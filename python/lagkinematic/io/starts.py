# python/lagkinematic/io/starts.py
from __future__ import annotations
from pathlib import Path
import numpy as np
from typing import Tuple, Dict, Any

def read_starts(path: str, expected_cols: int = 8) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Legge file starts:
      - Prima riga opzionale: numero totale coppie (IGNORATO ai fini logici)
      - Dati riga: lon1 lat1 dep1  lon2 lat2 dep2  t_delay  col8

    Ritorna:
      data: ndarray shape (N, expected_cols)
      info: {'declared_pairs': int|None, 'n_pairs': int, 'has_col8': bool}
    """
    p = Path(path).expanduser()
    if not p.exists():
        raise FileNotFoundError(f"File starts non trovato: {p}")

    with p.open("r") as f:
        lines = f.readlines()

    def _is_comment_or_empty(s: str) -> bool:
        s = s.strip()
        return (not s) or s.startswith("#") or s.startswith("//")

    # rileva eventuale riga dichiarata (ignorata se presente)
    declared = None
    first_payload_idx = 0
    for i, ln in enumerate(lines):
        if _is_comment_or_empty(ln):
            continue
        toks = ln.strip().split()
        if len(toks) == 1:
            try:
                declared = int(toks[0])
                first_payload_idx = i + 1
                break
            except ValueError:
                first_payload_idx = i
                break
        else:
            first_payload_idx = i
            break

    # parse dati
    rows = []
    for ln in lines[first_payload_idx:]:
        if _is_comment_or_empty(ln):
            continue
        vals = ln.strip().split()
        if len(vals) < expected_cols:
            vals = vals + ["0.0"] * (expected_cols - len(vals))
        if len(vals) != expected_cols:
            raise ValueError(f"Riga con {len(vals)} colonne ma attese {expected_cols}: '{ln.strip()}'")
        rows.append([float(x) for x in vals])

    data = np.asarray(rows, dtype=np.float64)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    info = {
        "declared_pairs": declared,    # solo informativa
        "n_pairs": int(data.shape[0]),
        "has_col8": data.shape[1] >= 8,
    }
    return data, info