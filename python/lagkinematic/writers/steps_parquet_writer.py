# python/lagkinematic/writers/steps_parquet_writer.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Any

import pandas as pd
import click

from lagkinematic.chunking import StepsWriter


@dataclass
class StepsParquetWriter(StepsWriter):
    """
    Writer per gli 'steps' su file Parquet, per-rank e per-chunk.

    Layout attuale (per run con outdir=...):

        outdir/
          steps/
            rank0000/
              part-r0000-c00000.parquet
              part-r0000-c00001.parquet
            rank0001/
              part-r0001-c00000.parquet
              ...

    Ogni file contiene tutte le righe accumulate in quel chunk temporale
    per il rank corrente.
    """

    base_dir: str
    rank: int
    layout: str = "time"  # per ora solo "time", ma lasciamo il campo per estensioni future

    def __post_init__(self) -> None:
        self.base_path = Path(self.base_dir).expanduser().absolute()
        self.steps_dir = self.base_path / "steps" / f"rank{self.rank:04d}"
        self.steps_dir.mkdir(parents=True, exist_ok=True)

    def _chunk_path(self, chunk_index: int) -> Path:
        """
        Costruisce il path del file parquet per un chunk dato.
        """
        fname = f"part-r{self.rank:04d}-c{chunk_index:05d}.parquet"
        return self.steps_dir / fname

    def write_chunk(self, chunk_index: int, rows: Iterable[dict[str, Any]]) -> None:
        """
        Scrive un chunk su parquet.
        """
        rows_list = list(rows)
        if not rows_list:
            return

        path = self._chunk_path(chunk_index)
        df = pd.DataFrame.from_records(rows_list)

        # Scrittura parquet (usando engine di default: pyarrow se disponibile)
        df.to_parquet(path, index=False)

        click.echo(
            f"[rank {self.rank}] StepsParquetWriter: wrote chunk={chunk_index} "
            f"rows={len(df)} → {path}",
            err=False,
        )
