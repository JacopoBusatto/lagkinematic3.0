from __future__ import annotations

from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from lagkinematic.writers.steps_compactor import compact_all_ranks, compact_rank_steps


COLUMNS = ["id", "time", "X1", "Y1", "Z1", "X2", "Y2", "Z2", "age"]


def _write_chunk(path: Path, rows: list[dict]) -> None:
    df = pd.DataFrame(rows, columns=COLUMNS)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=False)


def test_compact_rank_steps(tmp_path):
    rank_dir = tmp_path / "rank0000"

    chunk1 = [
        {"id": 2, "time": 2, "X1": 0.0, "Y1": 0.1, "Z1": 0.2, "X2": 0.3, "Y2": 0.4, "Z2": 0.5, "age": 1},
        {"id": 1, "time": 3, "X1": 1.0, "Y1": 1.1, "Z1": 1.2, "X2": 1.3, "Y2": 1.4, "Z2": 1.5, "age": 2},
    ]
    chunk2 = [
        {"id": 1, "time": 1, "X1": 2.0, "Y1": 2.1, "Z1": 2.2, "X2": 2.3, "Y2": 2.4, "Z2": 2.5, "age": 3},
    ]

    _write_chunk(rank_dir / "part-r0000-c00000.parquet", chunk1)
    _write_chunk(rank_dir / "part-r0000-c00001.parquet", chunk2)

    output_path = Path(compact_rank_steps(rank_dir))

    assert output_path.exists()

    compacted = pq.read_table(output_path)
    assert compacted.num_rows == len(chunk1) + len(chunk2)

    ids = compacted.column("id").to_pylist()
    times = compacted.column("time").to_pylist()
    assert list(zip(ids, times)) == sorted(zip(ids, times))


def test_compact_all_ranks(tmp_path):
    rank0 = tmp_path / "rank0000"
    rank1 = tmp_path / "rank0001"

    _write_chunk(
        rank0 / "part-r0000-c00000.parquet",
        [{"id": 2, "time": 2, "X1": 0, "Y1": 0, "Z1": 0, "X2": 0, "Y2": 0, "Z2": 0, "age": 0}],
    )
    _write_chunk(
        rank1 / "part-r0001-c00000.parquet",
        [{"id": 1, "time": 1, "X1": 1, "Y1": 1, "Z1": 1, "X2": 1, "Y2": 1, "Z2": 1, "age": 1}],
    )

    outputs = compact_all_ranks(tmp_path)

    expected_paths = {
        str(rank0 / "steps_rank0000.parquet"),
        str(rank1 / "steps_rank0001.parquet"),
    }
    assert set(outputs) == expected_paths
    for out in outputs:
        assert Path(out).exists()
