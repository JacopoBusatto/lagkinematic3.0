from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Sequence
import json


def write_manifest(
    output_dir: str,
    cfg: Dict[str, Any],
    meta: Dict[str, Any],
    steps_compacted: Sequence[str] | None = None,
) -> str:
    """
    Write a JSON manifest describing the run configuration and outputs.

    Returns the path to the manifest file.
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    payload = {
        "config_path": meta.get("config_path"),
        "run": {
            "start": cfg["run"].get("start"),
            "duration": cfg["run"].get("duration"),
            "dt_seconds": cfg["run"].get("dt_seconds"),
            "n_steps": meta.get("n_steps"),
        },
        "output": {
            "dir": output_dir,
            "parquet_compression": cfg["run"].get("output", {}).get(
                "parquet_compression", "snappy"
            ),
            "steps_compacted": list(steps_compacted or []),
        },
        "pairs": {
            "n_pairs_total": meta.get("n_pairs"),
        },
        "meta": meta,
    }

    manifest_path = out_dir / "manifest.json"
    manifest_path.write_text(json.dumps(payload, indent=2))
    return str(manifest_path)
