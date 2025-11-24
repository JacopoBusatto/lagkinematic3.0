# python/lagkinematic/integration.py
from __future__ import annotations
import numpy as np

from dataclasses import dataclass
from typing import Protocol, Optional

from .geometry import displace_lonlat, wrap_longitude, LonLatDomain


# ---------- Interfacce / protocolli ----------
def _is_longitude_periodic(domain: LonLatDomain) -> bool:
    """
    Ritorna True se il dominio è 'globale' in longitudine e quindi
    ha senso applicare condizioni periodiche.
    Usiamo una soglia larga ( >350° ) per coprire casi tipo 0.5–359.5.
    """
    lon_range = domain.lon_max - domain.lon_min
    return lon_range > 350.0

class VelocitySampler(Protocol):
    """
    Interfaccia minimale per il nostro RegularLatLonSampler.
    Deve restituire (u, v, w) in m/s dato (lon, lat, depth, time).
    """

    def sample(self, lon_deg: float, lat_deg: float, depth_m: float, t: float) -> tuple[float, float, float]:
        ...


class MaskSampler(Protocol):
    """
    Interfaccia per la maschera terra/mare.
    Deve restituire un valore scalare (es. 0/1 o fra 0 e 1) dato (lon, lat).
    """

    def sample_mask(self, lon_deg: float, lat_deg: float) -> float:
        ...


class SubgridModel(Protocol):
    """
    Interfaccia generica per il modello cinematico (2D/3D).
    Oggi: restituisce solo una velocità subgrid (u_sgs, v_sgs, w_sgs).
    Domani: qui dentro potrai implementare la logica
    'spegni scala se troppo vicino alla costa', usando la maschera.
    """

    def velocity(self, lon_deg: float, lat_deg: float, depth_m: float, t: float) -> tuple[float, float, float]:
        ...


# ---------- Stato delle particelle ----------

@dataclass
class ParticleState:
    """
    Stato di una singola particella in coordinate geografiche + profondità.
    """
    id: int
    lon: float  # gradi
    lat: float  # gradi
    depth: float  # metri (negativi in acqua, se segui il legacy)
    alive: bool = True

    def kill(self) -> None:
        self.alive = False


@dataclass
class ParticlePairState:
    """
    Stato di una coppia di particelle + tempo corrente.
    """
    id: int
    p1: ParticleState
    p2: ParticleState
    t: float = 0.0           # tempo in secondi dall'inizio del run
    age: float = 0.0         # età della coppia rispetto al rilascio (t - t_delay)


# ---------- Integratore Euler ----------
@dataclass
class EulerIntegrator:
    sampler: VelocitySampler
    domain: LonLatDomain
    dt: float
    mask: Optional[MaskSampler] = None
    mask_threshold: float = 0.5
    subgrid: Optional[SubgridModel] = None

    def apply_initial_mask(self, pairs: list[ParticlePairState]) -> None:
        """
        Uccide le particelle che iniziano sulla terra.
        Se almeno UNA particella in una coppia è su terra, la coppia viene
        considerata morta: uccidiamo ENTRAMBE.
        Va chiamato PRIMA del primo step, prima del primo writer.
        """
        if self.mask is None:
            return

        for pair in pairs:
            # calcola il valore di maschera per entrambe
            m1 = self.mask.sample_mask(pair.p1.lon, pair.p1.lat, pair.p1.depth)
            m2 = self.mask.sample_mask(pair.p2.lon, pair.p2.lat, pair.p2.depth)

            if (m1 < self.mask_threshold) or (m2 < self.mask_threshold):
                pair.p1.kill()
                pair.p2.kill()


    def _step_one_particle(self, p: ParticleState, t: float) -> None:
        if not p.alive:
            return

        # 1) Controllo maschera terra/mare
        if self.mask is not None:
            m_val = self.mask.sample_mask(p.lon, p.lat, p.depth)
            if m_val < self.mask_threshold:
                p.kill()
                return

        # 2) Velocità risolta (u,v,w) dal sampler
        u_res, v_res, w_res = self.sampler.sample(p.lon, p.lat, p.depth, t)

        # 3) Subgrid
        if self.subgrid is not None:
            u_sgs, v_sgs, w_sgs = self.subgrid.velocity(p.lon, p.lat, p.depth, t)
        else:
            u_sgs = v_sgs = w_sgs = 0.0

        # Velocità totale
        u_tot = u_res + u_sgs
        v_tot = v_res + v_sgs
        w_tot = w_res + w_sgs

        # 4) Spostamento in metri
        dx = u_tot * self.dt
        dy = v_tot * self.dt
        dz = w_tot * self.dt

        # 5) Aggiornamento lon/lat con formula legacy (via helper)
        lon_new, lat_new = displace_lonlat(
            lon_deg=p.lon,
            lat_deg=p.lat,
            dx_m=dx,
            dy_m=dy,
        )

        # 5bis) Trattamento bordo longitudinale:
        #       - se dominio "globale" → condizioni periodiche
        #       - se dominio regionale → kill se esce dal range lon/lat
        if _is_longitude_periodic(self.domain):
            # wrapping periodico della longitudine
            lon_new = wrap_longitude(lon_new, self.domain)
        else:
            # dominio non periodico: se esce, la particella muore
            if lon_new < self.domain.lon_min or lon_new > self.domain.lon_max:
                p.kill()
                return
        if lat_new < self.domain.lat_min or lat_new > self.domain.lat_max:
            p.kill()
            return
        
        # 6) Aggiornamento profondità (qui semplicemente Z + w*dt)
        depth_new = p.depth + dz

        # 7) Scrittura stato aggiornato
        p.lon = lon_new
        p.lat = lat_new
        p.depth = depth_new

    def step_pair(self, pair: ParticlePairState, t: float) -> None:
        """
        Esegue un passo di integrazione per la coppia al tempo globale t (secondi).
        Regola FSLE-style: se UNA delle due particelle muore (terra o bordo),
        la coppia viene considerata morta → entrambe killate e nessun update di t.
        """
        # se la coppia è già "morta" (una delle due non è viva), non facciamo nulla
        if (not pair.p1.alive) or (not pair.p2.alive):
            return

        # integriamo entrambe
        self._step_one_particle(pair.p1, t)
        self._step_one_particle(pair.p2, t)

        # dopo lo step, se una è morta → kill anche l'altra
        if (not pair.p1.alive) or (not pair.p2.alive):
            pair.p1.kill()
            pair.p2.kill()
            # non aggiorniamo pair.t: la coppia è finita qui
            return

        # se entrambe ancora vive → aggiorniamo il tempo della coppia
        pair.t = t + self.dt
