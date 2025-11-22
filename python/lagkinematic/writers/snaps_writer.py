# python/lagkinematic/writers/snaps_writer.py
from __future__ import annotations
from pathlib import Path
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

SNAP_COLUMNS = ['id','age','X','Y','Z','start_time','time','X0','Y0']

def write_snapshot_parquet_rank(
    out_dir: str,
    rank: int,
    step_idx: int,
    ids: np.ndarray,             # [M] int64 (globali)
    time_ns: np.ndarray,         # [M] datetime64[ns]
    age_s: np.ndarray,           # [M] int64
    X: np.ndarray, Y: np.ndarray, Z: np.ndarray,       # [M] float64
    start_time_ns: np.ndarray,   # [M] datetime64[ns]
    X0: np.ndarray, Y0: np.ndarray,                   # [M] float64
) -> str:
    """
    Scrive lo snapshot di questo rank in:
      {out_dir}/trajectories/snaps/rank-{rank:05d}/step-{step_idx:06d}.parquet
    Colonne: id, age, X, Y, Z, start_time, time, X0, Y0
    """
    snap_dir = Path(out_dir) / "trajectories" / "snaps" / f"rank-{rank:05d}"
    snap_dir.mkdir(parents=True, exist_ok=True)
    out_path = snap_dir / f"step-{step_idx:06d}.parquet"

    table = pa.table({
        "id": pa.array(ids, type=pa.int64()),
        "age": pa.array(age_s, type=pa.int64()),
        "X": pa.array(X, type=pa.float64()),
        "Y": pa.array(Y, type=pa.float64()),
        "Z": pa.array(Z, type=pa.float64()),
        "start_time": pa.array(start_time_ns.astype("datetime64[ns]")),
        "time": pa.array(time_ns.astype("datetime64[ns]")),
        "X0": pa.array(X0, type=pa.float64()),
        "Y0": pa.array(Y0, type=pa.float64()),
    })
    pq.write_table(table, out_path)
    return str(out_path)