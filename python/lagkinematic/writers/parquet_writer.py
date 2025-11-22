# python/lagkinematic/writers/parquet_writer.py
from __future__ import annotations
from pathlib import Path
import json
import pyarrow as pa
import pyarrow.parquet as pq

def write_initial_positions_rank(
    out_dir: str,
    rank: int,
    starts_local: "np.ndarray",
    start_date_iso: str,
) -> str:
    """
    Scrive le posizioni iniziali del SOLO rank in:
      {out_dir}/trajectories/part-r{rank:04d}-init.parquet

    Schema colonne:
      ['id','time','X1','Y1','Z1','X2','Y2','Z2','age']
    """
    import numpy as np
    from datetime import datetime, timedelta

    traj_dir = Path(out_dir) / "trajectories"
    traj_dir.mkdir(parents=True, exist_ok=True)
    out_path = traj_dir / f"part-r{rank:04d}-init.parquet"

    # id locale al rank
    ids = np.arange(starts_local.shape[0], dtype=np.int64)

    # tempi reali: start + t_delay (colonna 6)
    base = datetime.fromisoformat(start_date_iso)
    times_py = [base + timedelta(seconds=float(row[6])) for row in starts_local]
    times = np.array(times_py, dtype="datetime64[ns]")  # datetime64[ns]

    # posizioni
    X1, Y1, Z1, X2, Y2, Z2 = [starts_local[:, i] for i in range(6)]
    age = np.zeros(starts_local.shape[0], dtype=np.int64)  # secondi dall'inizio

    table = pa.table({
        "id": ids,
        "time": times,
        "X1": X1, "Y1": Y1, "Z1": Z1,
        "X2": X2, "Y2": Y2, "Z2": Z2,
        "age": age,
    })
    pq.write_table(table, out_path)  # file unico per rank
    return str(out_path)

def write_manifest(out_dir: str, cfg: dict, meta: dict) -> None:
    out = Path(out_dir) / "manifest.json"
    payload = {
        "lagkinematic_version": "3.0.0",
        "run": cfg.get("run", {}),
        "domain": {
            "type": cfg.get("domain", {}).get("type"),
            "files_mode": cfg.get("domain", {}).get("files", {}).get("mode", "list"),
            "num_files": len(cfg.get("domain", {}).get("_expanded_files", [])),
        },
        "meta": meta,
    }
    out.write_text(json.dumps(payload, indent=2))

def create_empty_parquet(out_dir: str) -> str:
    """
    Crea 'trajectories.parquet' VUOTO con schema:
      id:int64, time:timestamp[ns], X1:float64, Y1:float64, Z1:float64,
      X2:float64, Y2:float64, Z2:float64, age:int64 (secondi)
    """
    out_path = str(Path(out_dir) / "trajectories.parquet")
    schema = pa.schema([
        pa.field("id", pa.int64()),
        pa.field("time", pa.timestamp("ns")),
        pa.field("X1", pa.float64()),
        pa.field("Y1", pa.float64()),
        pa.field("Z1", pa.float64()),
        pa.field("X2", pa.float64()),
        pa.field("Y2", pa.float64()),
        pa.field("Z2", pa.float64()),
        pa.field("age", pa.int64()),  # secondi dall'inizio della traiettoria
    ])
    empty = pa.table({name: pa.array([], type=field.type) for name, field in zip(schema.names, schema)})
    pq.write_table(empty, out_path)
    return out_path