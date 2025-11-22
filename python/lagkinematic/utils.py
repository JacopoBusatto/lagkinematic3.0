# python/lagkinematic/utils.py
from __future__ import annotations


def is_snapshot_time(t: float, snapshot_seconds: float, *, t0: float = 0.0) -> bool:
    """
    Restituisce True se al tempo t va fatto uno snapshot, usando logica tipo Global::Snapshot:

        - snapshot a t == 0
        - poi ogni multiplo di snapshot_seconds (confronto in interi).

    Nota: assumiamo che t e snapshot_seconds siano in secondi.
    """
    if t <= t0:
        return True
    t_int = int(round(t - t0))
    s_int = int(round(snapshot_seconds))
    if s_int <= 0:
        return False
    return (t_int % s_int) == 0
