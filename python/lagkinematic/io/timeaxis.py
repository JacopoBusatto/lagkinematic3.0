# python/lagkinematic/io/timeaxis.py
from __future__ import annotations
from pathlib import Path
from typing import Dict, Any
import numpy as np
import xarray as xr

def _as_ns(arr) -> np.ndarray:
    """Converte un array di datetime64 (qualsiasi unità) in datetime64[ns]."""
    a = np.array(arr)
    if np.issubdtype(a.dtype, np.datetime64):
        return a.astype('datetime64[ns]')
    # se fosse numerico con units strane, qui si potrebbe usare decode_cf=False e units...
    # per ora richiediamo decode_cf True nel reader.
    raise TypeError("L'asse temporale non è datetime64: serve decode_cf=True nei NetCDF.")

def _median_dt_seconds(t_ns: np.ndarray) -> float:
    if t_ns.size < 2:
        return np.nan
    diffs = np.diff(t_ns).astype('timedelta64[s]').astype(np.float64)
    return float(np.median(diffs))

def compute_time_axis_regular_latlon(cfg: Dict[str, Any]) -> Dict[str, Any]:
    """
    Legge e concatena l'asse temporale dai NetCDF (solo la coordinata tempo),
    e se run.time_loops > 0 replica il ciclo temporale aggiungendo un offset
    pari alla durata del ciclo (per avere un asse monotonico).

    Ritorna:
      {
        'time': np.ndarray(datetime64[ns]),
        'native_dt_seconds': float,
        'cycle_length_seconds': float | None,
        'n_cycles': int
      }
    """
    files: list[str] = cfg["domain"]["_expanded_files"]
    tname: str = cfg["domain"]["coords"]["time"]
    loops: int = int(cfg["run"].get("time_loops", 0))

    # 1) leggi time da ogni file (decode_cf=True per avere datetime64)
    times = []
    for fp in files:
        p = Path(fp)
        if not p.exists():
            raise FileNotFoundError(f"File non trovato durante la lettura dell'asse tempo: {p}")
        with xr.open_dataset(p, decode_times=True) as ds:
            if tname not in ds.coords:
                raise KeyError(f"Coordinata tempo '{tname}' non trovata in {p.name}. Disponibili: {list(ds.coords)}")
            t = _as_ns(ds[tname].values)
            times.append(t)

    if not times:
        raise ValueError("Nessun asse temporale letto dai file.")

    # concat singolo ciclo
    t_cycle = np.concatenate(times)
    # assicurati ordinamento e univocità soft (non error se duplicati)
    t_cycle = np.unique(t_cycle)

    native_dt = _median_dt_seconds(t_cycle)

    # --- controllo regolarità temporale ---
    if t_cycle.size > 2:
        diffs = np.diff(t_cycle).astype("timedelta64[s]").astype(float)
        med = np.nanmedian(diffs)
        mad = np.nanmedian(np.abs(diffs - med))  # robust MAD
        rel_spread = mad / med if med != 0 else np.inf
        if rel_spread > 1e-3:  # ~0.1% di variazione
            print(
                f"[warning] intervallo temporale irregolare: "
                f"Δt_med={med:.1f}s  MAD={mad:.1f}s ({rel_spread*100:.2f}% spread)"
            )
    # 2) se loops > 0, replica il ciclo con offset temporale
    if loops > 0 and t_cycle.size >= 2:
        # lunghezza ciclo = (last - first) + native_dt (così il 1° tempo del ciclo successivo non collide)
        cycle_len_s = (
            (t_cycle[-1] - t_cycle[0]).astype('timedelta64[s]').astype(np.int64)
            + int(round(native_dt if np.isfinite(native_dt) else 0))
        )
        chunks = [t_cycle]
        for k in range(1, loops + 1):
            offset = np.timedelta64(cycle_len_s * k, 's')
            chunks.append(t_cycle + offset)
        time_all = np.concatenate(chunks)
        n_cycles = 1 + loops
        out = {
            "time": time_all,
            "native_dt_seconds": native_dt,
            "cycle_length_seconds": float(cycle_len_s),
            "n_cycles": n_cycles,
        }
    else:
        out = {
            "time": t_cycle,
            "native_dt_seconds": native_dt,
            "cycle_length_seconds": None,
            "n_cycles": 1,
        }

    return out
