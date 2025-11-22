# python/lagkinematic/io/filespec.py
from __future__ import annotations
from pathlib import Path
from datetime import datetime, timedelta
from typing import List, Dict

def _date_range(start: str, end: str | None, freq: str, fallback_duration_days: int | None) -> List[datetime]:
    if end is None and fallback_duration_days is None:
        raise ValueError("Fornisci 'end' oppure run.duration (es. '180d').")
    ds = datetime.fromisoformat(start)
    de = datetime.fromisoformat(end) if end else ds + timedelta(days=fallback_duration_days or 0)

    if not freq.endswith("D"):
        raise NotImplementedError("Per ora supportiamo solo frequenze a giorni: es. '1D'")
    step_days = int(freq[:-1]) if freq[:-1] else 1
    step = timedelta(days=step_days)

    out, t = [], ds
    while t <= de:
        out.append(t)
        t += step
    return out

def expand_files(run_cfg: Dict, files_cfg: Dict) -> List[str]:
    mode = files_cfg.get("mode", "list")
    if mode == "list":
        return list(files_cfg.get("list", []))
    if mode == "listfile":
        p = Path(files_cfg["path"]).expanduser()
        with p.open() as f:
            return [ln.strip() for ln in f if ln.strip()]
    if mode == "daterange":
        base = Path(files_cfg["base_dir"]).expanduser()
        template = files_cfg["template"]  # es. "dataset_{date:%Y%m%d}.nc"
        # ✅ usa run.start se 'start' non è specificato nel blocco files
        start = files_cfg.get("start", run_cfg.get("start"))
        if start is None:
            raise ValueError("Devi fornire 'start' in domain.files oppure 'run.start' nel config.")
        end  = files_cfg.get("end")
        freq = files_cfg.get("freq", "1D")
        dur = run_cfg.get("duration", "")
        fallback_days = int(dur[:-1]) if isinstance(dur, str) and dur.endswith("d") and dur[:-1].isdigit() else None
        dates = _date_range(start, end, freq, fallback_days)
        return [str((base / template.format(date=dt)).absolute()) for dt in dates]

    raise ValueError(f"files.mode sconosciuto: {mode}")
