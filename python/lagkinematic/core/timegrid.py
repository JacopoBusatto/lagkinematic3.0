# python/lagkinematic/core/timegrid.py
from __future__ import annotations
import numpy as np
from datetime import datetime, timedelta

def parse_duration_seconds(s: str) -> int:
    s = s.strip().lower()
    if s.endswith("d") and s[:-1].isdigit():
        return int(s[:-1]) * 86400
    raise ValueError(f"Formato duration non supportato: {s}")

def build_time_grid(start_iso: str, duration_str: str, dt_seconds: int):
    """
    Ritorna:
      t_grid    : np.ndarray(datetime64[ns]) [N]
      t_offsets : np.ndarray(int64 seconds)  [N] (0, dt, 2dt, ...)
    """
    start = np.datetime64(datetime.fromisoformat(start_iso), 'ns')
    dur_s = parse_duration_seconds(duration_str)
    # N inclusivo: t = 0, dt, ..., dur
    N = dur_s // dt_seconds + 1
    offs = np.arange(N, dtype=np.int64) * int(dt_seconds)
    t_grid = start + offs.astype(f'timedelta64[s]').astype('timedelta64[ns]')
    return t_grid, offs
