from __future__ import annotations

from pathlib import Path
import glob
from typing import List

import pyarrow as pa
import pyarrow.parquet as pq


def _discover_chunk_files(rank_dir: Path) -> list[str]:
    pattern = str(rank_dir / "part-r*-c*.parquet")
    files = sorted(glob.glob(pattern))
    return files


def _build_output_path(rank_dir: Path, output_file: str | None) -> Path:
    if output_file is not None:
        return Path(output_file)
    rank_name = rank_dir.name
    return rank_dir / f"steps_{rank_name}.parquet"


def compact_rank_steps(
    rank_dir: str, output_file: str | None = None, compression: str | None = "snappy"
) -> str:
    """
    Concatenate all step chunk Parquet files within a rank directory.

    Parameters
    ----------
    rank_dir:
        Directory containing ``part-r*-c*.parquet`` files for a rank.
    output_file:
        Optional output path. If omitted, writes ``steps_<rank>.parquet``
        inside ``rank_dir``.

    Returns
    -------
    str
        Path to the compacted Parquet file.

    Raises
    ------
    ValueError
        If no chunk files are discovered in ``rank_dir``.
    """
    rank_path = Path(rank_dir)
    chunk_files = _discover_chunk_files(rank_path)
    if not chunk_files:
        raise ValueError(f"No chunk Parquet files found in {rank_dir}")

    tables: List[pa.Table] = [pq.read_table(path) for path in chunk_files]
    combined = pa.concat_tables(tables, promote_options="default")

    if {"id", "time"}.issubset(set(combined.column_names)):
        combined = combined.sort_by([("id", "ascending"), ("time", "ascending")])

    output_path = _build_output_path(rank_path, output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(combined, output_path, compression=compression)

    return str(output_path)


def compact_all_ranks(steps_root: str, compression: str | None = "snappy") -> list[str]:
    """
    Compact step chunk files for every rank directory under ``steps_root``.

    Parameters
    ----------
    steps_root:
        Path containing per-rank subdirectories named like ``rank0000``.

    Returns
    -------
    list[str]
        Paths to the compacted Parquet files for each rank directory.
    """
    root_path = Path(steps_root)
    rank_dirs = sorted(
        path for path in root_path.iterdir() if path.is_dir() and path.name.startswith("rank")
    )

    compacted_files: list[str] = []
    for rank_dir in rank_dirs:
        compacted_files.append(compact_rank_steps(str(rank_dir), compression=compression))

    return compacted_files
