# python/lagkinematic/chunking.py
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Protocol, Any, Iterable

from .integration import ParticlePairState


class StepsWriter(Protocol):
    """
    Interfaccia minimale per il writer degli 'steps' (parquet o altro).
    L'implementazione concreta la avrai da qualche parte nel tuo codice.
    """

    def write_chunk(self, chunk_index: int, rows: Iterable[dict[str, Any]]) -> None:
        ...


@dataclass
class TimeChunkBuffer:
    """
    Buffer in memoria per un singolo chunk temporale.
    Accumula righe (dizionari) e, quando flushato, le passa al writer.
    """
    chunk_index: int
    writer: StepsWriter
    rows: list[dict[str, Any]] = field(default_factory=list)

    def add_pair_step(self, pair: ParticlePairState) -> None:
        """
        Aggiunge una riga per la coppia allo step corrente.
        Usa pair.t come tempo (secondi dall'inizio del run).
        """
        t = pair.t
        row = {
            "id": pair.id,
            "time": t,
            "X1": pair.p1.lon,
            "Y1": pair.p1.lat,
            "Z1": pair.p1.depth,
            "X2": pair.p2.lon,
            "Y2": pair.p2.lat,
            "Z2": pair.p2.depth,
            # age per ora fissato a 0.0 (lo definiremo meglio in futuro)
            "age": 0.0,
        }
        self.rows.append(row)

    def flush(self) -> None:
        """
        Scrive il contenuto corrente tramite il writer e svuota il buffer.
        """
        if not self.rows:
            return
        self.writer.write_chunk(self.chunk_index, self.rows)
        self.rows.clear()


@dataclass
class TimeChunkManager:
    """
    Gestisce i chunk temporali:

    - decide quale sia il chunk corrente dato t
    - crea/distrugge il TimeChunkBuffer quando si passa a un nuovo chunk
    - fornisce un metodo comodo 'add_step' da usare nel loop di integrazione
    """
    steps_writer: StepsWriter
    chunk_size_seconds: float
    t0: float = 0.0  # tempo di riferimento (inizio run)
    _current_chunk_index: int | None = field(init=False, default=None)
    _buffer: TimeChunkBuffer | None = field(init=False, default=None)

    def _compute_chunk_index(self, t: float) -> int:
        """
        Determina l'indice di chunk dato il tempo t.
        """
        return int((t - self.t0) // self.chunk_size_seconds)

    def _ensure_buffer_for_time(self, t: float) -> None:
        """
        Verifica che esista un buffer per il chunk corrispondente a t.
        Se il chunk cambia, fa flush del buffer vecchio e ne crea uno nuovo.
        """
        new_index = self._compute_chunk_index(t)

        if self._current_chunk_index is None:
            # Primo chunk
            self._current_chunk_index = new_index
            self._buffer = TimeChunkBuffer(chunk_index=new_index, writer=self.steps_writer)
            return

        if new_index != self._current_chunk_index:
            # Cambio chunk: flush del vecchio buffer e creazione del nuovo
            assert self._buffer is not None
            self._buffer.flush()
            self._current_chunk_index = new_index
            self._buffer = TimeChunkBuffer(chunk_index=new_index, writer=self.steps_writer)

    def add_step(self, pair: ParticlePairState) -> None:
        """
        Aggiunge la riga per la coppia nel buffer del chunk corretto,
        usando pair.t come tempo globale.
        """
        t = pair.t
        self._ensure_buffer_for_time(t)
        assert self._buffer is not None
        self._buffer.add_pair_step(pair)

    def finalize(self) -> None:
        """
        Da chiamare a fine simulazione per flushare l'ultimo chunk.
        """
        if self._buffer is not None:
            self._buffer.flush()
