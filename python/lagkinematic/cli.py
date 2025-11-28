from __future__ import annotations
import sys, pathlib, yaml, click
from datetime import datetime, timedelta
import numpy as np
import xarray as xr
# --- moduli locali ---
from lagkinematic.io.filespec import expand_files
from lagkinematic.io.regular_latlon import validate_regular_latlon
from lagkinematic.io.starts import read_starts
# ⬇️ usa SOLO il writer per-rank
from lagkinematic.writers.parquet_writer import write_manifest, write_initial_positions_rank
from lagkinematic.writers.steps_parquet_writer import StepsParquetWriter
from lagkinematic.io.timeaxis import compute_time_axis_regular_latlon
from lagkinematic.core.timegrid import build_time_grid, parse_duration_seconds
from lagkinematic.geometry import LonLatDomain
from lagkinematic.integration import ParticleState, ParticlePairState, EulerIntegrator
from lagkinematic.chunking import TimeChunkManager
from lagkinematic.utils import is_snapshot_time

from lagkinematic.sampling.regular_latlon import RegularLatLonSampler
from lagkinematic.sampling.mask_regular_latlon import RegularLatLonMaskSampler
from lagkinematic.sampling.bottom_regular_latlon import RegularLatLonBottomSampler

# --- MPI fallback ---
try:
    from mpi4py import MPI  # type: ignore
    COMM = MPI.COMM_WORLD
    RANK = COMM.Get_rank()
    SIZE = COMM.Get_size()
    def bcast(obj, root=0):
        return COMM.bcast(obj, root=root)
except Exception:
    COMM = None
    RANK = 0
    SIZE = 1
    def bcast(obj, root=0):
        return obj

def _partition_indices(n: int, size: int, rank: int) -> tuple[int, int]:
    q, r = divmod(n, size)
    start = rank * q + min(rank, r)
    stop = start + q + (1 if rank < r else 0)
    return start, stop

def _parse_duration_seconds(s: str) -> int:
    s = s.strip().lower()
    if s.endswith("d") and s[:-1].isdigit():
        return int(s[:-1]) * 86400
    raise ValueError(f"Formato duration non supportato: {s}")

@click.command()
@click.option(
    "--config", "-c",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help="File YAML di configurazione della simulazione.",
)
def main(config: str):
    if RANK == 0:
        click.echo("LAGKINEMATIC 3.0 — Python controller")
        with open(config, "r") as f:
            cfg = yaml.safe_load(f)
        outdir = pathlib.Path(cfg["run"]["output"]["dir"]).expanduser().absolute()
        outdir.mkdir(parents=True, exist_ok=True)
        cfg["run"]["output"]["dir"] = str(outdir)
    else:
        cfg = None

    cfg = bcast(cfg, root=0)

    # 1) Espansione file
    try:
        file_list = expand_files(cfg["run"], cfg["domain"]["files"])
    except Exception as e:
        if RANK == 0:
            click.echo(f"[errore] espansione files: {e}", err=True)
        sys.exit(1)
    cfg["domain"]["_expanded_files"] = file_list

    # 2) Validazione dominio
    try:
        if cfg["domain"]["type"] == "RegularLatLon":
            info = validate_regular_latlon(cfg)
        else:
            raise NotImplementedError(f"domain.type non supportato: {cfg['domain']['type']}")
    except Exception as e:
        if RANK == 0:
            click.echo(f"[errore] validazione dataset: {e}", err=True)
        raise SystemExit(2)
    
    # --- 3bis. Asse temporale dal dataset (+ loop temporali) ---
    tinfo = compute_time_axis_regular_latlon(cfg)
    t_all = tinfo["time"]                       # datetime64[ns]
    native_dt_s = tinfo["native_dt_seconds"]    # tipicamente 86400 per daily, ecc.
    n_cycles = tinfo["n_cycles"]

    # 3ter) Costruisci il sampler RegularLatLon per le velocità
    sampler_core = RegularLatLonSampler.from_config(cfg)

    # 3ter-bis) Costruzione della maschera se richiesta (NaN-based)
    mask_cfg = cfg.get("mask", {})
    if mask_cfg.get("enabled", False):
        mask_sampler = RegularLatLonMaskSampler.from_first_file(cfg)
        mask_threshold = float(mask_cfg.get("threshold", 0.5))
    else:
        mask_sampler = None
        mask_threshold = 0.5

    # Sampler batimetrico opzionale (se esiste un asse depth)
    bottom_sampler = None

    # Modalità di gestione spiaggiamento e fondale
    beaching_mode = mask_cfg.get("beaching_mode", "kill")  # "kill" o "bounce"
    bottom_mode = mask_cfg.get("bottom_mode", "kill")      # "kill" o "bounce"

    # 3) Lettura STARTS
    try:
        starts_data, starts_info = read_starts(cfg["starts"]["file"], expected_cols=8)
    except Exception as e:
        if RANK == 0:
            click.echo(f"[errore] lettura starts: {e}", err=True)
        raise SystemExit(3)

    # Partizione locale per rank
    n_pairs_total = int(starts_info["n_pairs"])
    lo, hi = _partition_indices(n_pairs_total, SIZE, RANK)
    starts_local = starts_data[lo:hi]
    local_count = starts_local.shape[0]
    ids_global = (np.arange(local_count, dtype=np.int64) + lo)  # 0-based globale per questo rank

    # 4) Calcoli temporali base
    dt_s = int(cfg["run"]["dt_seconds"])
    dur_s = _parse_duration_seconds(cfg["run"]["duration"])
    n_steps = (dur_s // dt_s) + 1
    snap_s = int(cfg["run"]["snapshots_seconds"])
    snaps = list(range(0, n_steps * dt_s + 1, snap_s))
    # 4bis) Griglia temporale globale di integrazione
    t_grid, t_offs = build_time_grid(cfg["run"]["start"], cfg["run"]["duration"], dt_s)

    # Indice di start per ciascuna coppia locale in base a t_delay
    # t_delay è colonna 6 degli starts (in secondi)
    tdelay_local = starts_local[:, 6].astype(np.float64)
    # idx_start = ceil(t_delay / dt) clampato a [0, N-1]
    idx_start = np.ceil(tdelay_local / dt_s).astype(np.int64)
    idx_start = np.clip(idx_start, 0, t_grid.size - 1)

    # warning se non multiplo di dt
    frac = (tdelay_local / dt_s) - np.floor(tdelay_local / dt_s)
    n_nonint = int(np.count_nonzero(frac > 1e-9))
    if n_nonint and RANK == 0:
        click.echo(f"[warning] {n_nonint} start delay non multipli di dt={dt_s}s → uso ceil (inizio al primo step utile)")

    # Salva un piccolo riassunto (solo rank 0)
    if RANK == 0:
        click.echo(f"  t_grid: N={t_grid.size}  dt={dt_s}s  first={str(t_grid[0])}  last={str(t_grid[-1])}")

    # 5) Manifest (nota 'partitioning')
    meta = {
        "n_steps": int(n_steps),
        "dt_seconds": dt_s,
        "snapshots_seconds": snap_s,
        "estimated_snapshots": len(snaps),
        "first_file_checked": info["file"],
        "dims_first_file": info["dims"],
        "dtypes_first_file": info["dtypes"],
        "n_pairs": n_pairs_total,
        "parquet_schema": ['id','time','X1','Y1','Z1','X2','Y2','Z2','age'],
        "partitioning": "per-rank-files",
        "time_axis": {
            "native_dt_seconds": native_dt_s,
            "n_cycles": n_cycles,
            "cycle_length_seconds": tinfo["cycle_length_seconds"],
            "n_points": int(t_all.size),
            "first_iso": str(t_all[0]) if t_all.size else None,
            "last_iso":  str(t_all[-1]) if t_all.size else None,
        },
        "time_loops": int(cfg["run"].get("time_loops", 0)),
        "integration_time_grid": {
            "n_points": int(t_grid.size),
            "dt_seconds": dt_s,
            "first_iso": str(t_grid[0]),
            "last_iso": str(t_grid[-1]),
        },
        "id_offset_rank": {"rank": RANK, "offset": int(lo), "count": int(local_count)},
        "output": {
            "mode": cfg["run"]["output"].get("mode", "steps"),
            "steps_layout": cfg["run"]["output"].get("steps_layout", "time"),
            "steps_chunk_size": int(cfg["run"]["output"].get("steps_chunk_size", 1000)),
        },
    }
    write_manifest(cfg["run"]["output"]["dir"], cfg, meta)

    # 6) Scrivi posizioni iniziali per-rank (ogni rank il suo file)
    written_path = write_initial_positions_rank(
        cfg["run"]["output"]["dir"],
        RANK,
        starts_local,
        cfg["run"]["start"],
    )

    # 7) Stampa riepilogo (rank 0)
    if RANK == 0:
        click.echo(f"[lagk] ranks={SIZE}")
        click.echo(f"  start={cfg['run']['start']} duration={cfg['run']['duration']} dt={cfg['run']['dt_seconds']}s")
        click.echo(f"  snapshots_every={cfg['run']['snapshots_seconds']}s back={cfg['run']['backtrajectory']}")
        click.echo(f"  domain={cfg['domain']['type']}")
        click.echo(f"  coords: lon={cfg['domain']['coords']['lon']} lat={cfg['domain']['coords']['lat']} depth={cfg['domain']['coords'].get('depth')}")
        uvw = cfg['domain']['velocity']
        click.echo(f"  vars: u={uvw['u']['var']} v={uvw['v']['var']} w={uvw.get('w',{}).get('var')}")
        mode = cfg['domain']['files'].get('mode', 'list')
        click.echo(f"  files.mode={mode}")
        click.echo(f"  #files={len(cfg['domain']['_expanded_files'])}")
        click.echo(f"  starts={cfg['starts']['file']} outdir={cfg['run']['output']['dir']}")
        click.echo(f"  primo_file={info['file']}")
        click.echo(f"  dims={info['dims']}")
        click.echo(f"  dtypes={info['dtypes']}")
        click.echo(f"  steps={n_steps} (dt={dt_s}s, duration={cfg['run']['duration']})")
        click.echo(f"  snapshots≈{len(snaps)} (ogni {snap_s}s)")
        click.echo(f"  wrote initial positions → {written_path} (per-rank)")
        click.echo(f"  pairs_total={n_pairs_total}  distribuzione=~{n_pairs_total//SIZE}±1 per rank")
        click.echo(f"  time_axis: native_dt≈{native_dt_s:.1f}s  points={t_all.size}  cycles={n_cycles}")
        if n_cycles > 1:
            click.echo(f"  time_loops: {cfg['run'].get('time_loops', 0)}  cycle_length≈{tinfo['cycle_length_seconds']}s")

        # Anteprima 2 righe
        prev = starts_data[:2]
        for i, r in enumerate(prev):
            lon1, lat1, dep1, lon2, lat2, dep2, tdelay, col8 = r
            click.echo(
                f"    start[{i}] p1=({lon1:.6f},{lat1:.6f},{dep1:.3f})  "
                f"p2=({lon2:.6f},{lat2:.3f},{dep2:.3f})  t_delay={tdelay:.1f}  col8={col8:.3f}"
            )

    # 8) Integrazione Euler minimale di test
    # ---------------------------------------------------------------
    # Costruisco le coppie di particelle per questo rank
    pairs: list[ParticlePairState] = []
    for i in range(local_count):
        lon1, lat1, dep1, lon2, lat2, dep2, tdelay, col8 = starts_local[i]

        pid_pair = int(ids_global[i])
        p1 = ParticleState(
            id=int(2 * pid_pair),
            lon=float(lon1),
            lat=float(lat1),
            depth=float(dep1),
        )
        p2 = ParticleState(
            id=int(2 * pid_pair + 1),
            lon=float(lon2),
            lat=float(lat2),
            depth=float(dep2),
        )
        pair = ParticlePairState(id=pid_pair, p1=p1, p2=p2, t=0.0)
        pairs.append(pair)

    # Adapter: da (t_sec) → np.datetime64 per usare RegularLatLonSampler.sample_uv
    class SamplerAdapter:
        def __init__(self, core_sampler, start_iso: str):
            self.core = core_sampler
            self.t0 = np.datetime64(start_iso)

        def sample(self, lon_deg: float, lat_deg: float, depth_m: float, t_sec: float):
            t_ns = self.t0 + np.timedelta64(int(round(t_sec)), "s")
            return self.core.sample_uv(lon_deg, lat_deg, depth_m, t_ns)

    sampler = SamplerAdapter(sampler_core, cfg["run"]["start"])
    
    # Dominio lon/lat dal primo file NetCDF del dominio
    first_file = cfg["domain"]["_expanded_files"][0]
    lon_name = cfg["domain"]["coords"]["lon"]
    lat_name = cfg["domain"]["coords"]["lat"]

    with xr.open_dataset(first_file, decode_times=False) as ds:
        lons = ds[lon_name].values
        lats = ds[lat_name].values

    domain = LonLatDomain(
        lon_min=float(np.nanmin(lons)),
        lon_max=float(np.nanmax(lons)),
        lat_min=float(np.nanmin(lats)),
        lat_max=float(np.nanmax(lats)),
    )

    # Profondità massima globale (se c'è un asse depth)
    depth_name = cfg["domain"]["coords"].get("depth", None)
    max_depth = None
    if depth_name is not None:
        bottom_sampler = RegularLatLonBottomSampler.from_first_file(cfg)
        max_depth = float(np.nanmax(bottom_sampler.bottom.values))



    # Modello subgrid (per ora assente)
    subgrid_model = None

    # Integratore Euler con sampler reale (adapter) + maschera opzionale
    integrator = EulerIntegrator(
        sampler=sampler,
        domain=domain,
        dt=dt_s,
        mask=mask_sampler,
        mask_threshold=mask_threshold,
        subgrid=subgrid_model,
        beaching_mode=beaching_mode,
        bottom_mode=bottom_mode,
        bottom=bottom_sampler,
        max_depth=max_depth,
    )

    # 🔴 Kill iniziale: uccidi le coppie che partono su terra PRIMA di scrivere qualsiasi output
    integrator.apply_initial_mask(pairs)

    # Writer reale per gli steps: salva i chunk su Parquet per-rank
    steps_writer = StepsParquetWriter(
        base_dir=cfg["run"]["output"]["dir"],
        rank=RANK,
        layout=cfg["run"]["output"].get("steps_layout", "time"),
    )

    # TimeChunkManager: usa steps_chunk_size come dimensione in secondi del chunk
    steps_chunk_size = int(cfg["run"]["output"]["steps_chunk_size"])
    chunk_manager = TimeChunkManager(
        steps_writer=steps_writer,
        chunk_size_seconds=steps_chunk_size,
        t0=0.0,
    )

    # Snapshot writer dummy: per ora stampa solo qualcosa
    def dummy_snapshot_writer(pairs_local, t_sec: float) -> None:
        # usa solo la prima coppia viva, se esiste
        live = [p for p in pairs_local if p.p1.alive and p.p2.alive]
        if not live:
            return
        p0 = live[0]
        click.echo(
            f"[rank {RANK}] SNAPSHOT t={t_sec}s "
            f"pair0.id={p0.id} "
            f"p1=({p0.p1.lon:.4f},{p0.p1.lat:.4f}) "
            f"p2=({p0.p2.lon:.4f},{p0.p2.lat:.4f})"
        )

    # Mini-loop di integrazione: per ora solo pochi step di test
    n_test_steps = min(3, n_steps)

    if RANK == 0:
        click.echo(f"[lagk] Eseguo integrazione di test: n_test_steps={n_test_steps} dt={dt_s}s")

    for step_idx in range(n_test_steps):
        t_sec = step_idx * dt_s

        # 1) OUTPUT allo stato corrente (prima dell'integrazione)
        if (t_sec % snap_s) == 0:
            for i, pair in enumerate(pairs):
                if step_idx < idx_start[i]:
                    continue

                # se la coppia è morta, non scriviamo
                if (not pair.p1.alive) or (not pair.p2.alive):
                    continue

                # calcolo dell'età
                t_delay = float(tdelay_local[i])
                pair.age = max(0.0, t_sec - t_delay)

                # tempo attuale
                pair.t = float(t_sec)

                # scrittura singola
                chunk_manager.add_step(pair)

            # snapshot (stampa) – passa tutte le coppie, ma la funzione filtra le vive
            dummy_snapshot_writer(pairs, t_sec)

        # 2) INTEGRAZIONE da t → t + dt
        for i, pair in enumerate(pairs):
            if step_idx < idx_start[i]:
                continue

            # se la coppia è già morta, non integriamo
            if (not pair.p1.alive) or (not pair.p2.alive):
                continue

            integrator.step_pair(pair, t_sec)

    # Flush finale
    chunk_manager.finalize()

    if RANK == 0:
        click.echo("[lagk] Integrazione Euler di test completata (RegularLatLonSampler + parquet writer)")
