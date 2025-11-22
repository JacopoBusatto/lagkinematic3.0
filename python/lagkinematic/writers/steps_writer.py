# python/lagkinematic/writers/steps_writer.py
from __future__ import annotations
from pathlib import Path
from typing import Dict, Any, List
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

STEP_COLUMNS = ['id','time','X1','Y1','Z1','X2','Y2','Z2','age']

def _steps_dir(out_dir: str, rank: int) -> Path:
    p = Path(out_dir) / "trajectories" / "steps" / f"rank-{rank:05d}"
    p.mkdir(parents=True, exist_ok=True)
    return p

def write_steps_chunk_rank(
    out_dir: str,
    rank: int,
    first_step_idx: int,
    last_step_idx: int,
    rows: Dict[str, np.ndarray],  # chiavi = STEP_COLUMNS
) -> str:
    """
    Scrive un chunk di step per questo rank:
      {out_dir}/trajectories/steps/rank-00000/part-r00000-steps-XXXXXX-YYYYYY.parquet
    rows: dizionario di array numpy ALLINEATI sulle colonne STEP_COLUMNS.
    """
    steps_dir = _steps_dir(out_dir, rank)
    out_path = steps_dir / f"part-r{rank:05d}-steps-{first_step_idx:06d}-{last_step_idx:06d}.parquet"
    table = pa.table({k: pa.array(rows[k]) for k in STEP_COLUMNS})
    pq.write_table(table, out_path)
    return str(out_path)

class TimeChunkBuffer:
    """
    Accumula righe per step e flush-a a file ogni 'chunk_size' step.
    Uso:
      buf = TimeChunkBuffer(out_dir, rank, chunk_size)
      buf.add(step_idx, batch_rows_dict)  # batch_rows_dict con colonne STEP_COLUMNS
      buf.flush_final()  # a fine run
    """
    def __init__(self, out_dir: str, rank: int, chunk_size: int):
        self.out_dir = out_dir
        self.rank = rank
        self.chunk_size = int(chunk_size)
        self._current_first: int | None = None
        self._current_last: int | None = None
        self._cols: Dict[str, List] = {k: [] for k in STEP_COLUMNS}

    def add(self, step_idx: int, rows: Dict[str, np.ndarray]) -> list[str]:
        # allinea i batch (append lista Python → meno allocazioni frequenti)
        for k in STEP_COLUMNS:
            self._cols[k].append(rows[k])
        if self._current_first is None:
            self._current_first = step_idx
        self._current_last = step_idx

        # se ho raggiunto un confine di chunk (per indice di step), flush
        written: list[str] = []
        if ((self._current_last - self._current_first + 1) >= self.chunk_size):
            written.append(self._flush())
        return written

    def _flush(self) -> str:
        assert self._current_first is not None and self._current_last is not None
        # concatena ogni colonna
        cat = {k: np.concatenate(self._cols[k], axis=0) if len(self._cols[k]) > 0 else np.empty((0,)) for k in STEP_COLUMNS}
        path = write_steps_chunk_rank(self.out_dir, self.rank, self._current_first, self._current_last, cat)
        # reset
        self._current_first, self._current_last = None, None
        self._cols = {k: [] for k in STEP_COLUMNS}
        return path

    def flush_final(self) -> list[str]:
        if self._current_first is None:
            return []
        return [self._flush()]