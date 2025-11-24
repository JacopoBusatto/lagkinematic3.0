# python/lagkinematic/integration.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, Optional

from .geometry import displace_lonlat, wrap_longitude, LonLatDomain


# ---------- Interfacce / protocolli ----------

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
    """
    Integra le traiettorie con schema di Eulero esplicito in lon/lat/depth.

    - sampler: RegularLatLonSampler (u, v, w) risolto
    - domain: dominio longitudinale per il wrapping
    - dt: passo temporale in secondi
    - mask: maschera terra/mare (opzionale)
    - mask_threshold: soglia minima per considerare "mare"
    - subgrid: modello cinematico (opzionale)
    """
    sampler: VelocitySampler
    domain: LonLatDomain
    dt: float
    mask: Optional[MaskSampler] = None
    mask_threshold: float = 0.5
    subgrid: Optional[SubgridModel] = None

    def _step_one_particle(self, p: ParticleState, t: float) -> None:
        """
        Aggiorna IN-PLACE la particella p di un passo dt.
        """
        if not p.alive:
            return

        # 1) Controllo maschera terra/mare (semplice 0-1 o continuo)
        if self.mask is not None:
            m_val = self.mask.sample_mask(p.lon, p.lat)
            if m_val < self.mask_threshold:
                # Particella "assorbita" dalla terra: la marchiamo come morta.
                p.kill()
                return

        # 2) Velocità risolta (u,v,w) dal sampler
        u_res, v_res, w_res = self.sampler.sample(p.lon, p.lat, p.depth, t)

        # 3) Velocità subgrid (cinematica), se presente
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
        # wrapping periodico della longitudine in base al dominio
        lon_new = wrap_longitude(lon_new, self.domain)

        # 6) Aggiornamento profondità (qui semplicemente Z + w*dt)
        depth_new = p.depth + dz

        # 7) Scrittura stato aggiornato
        p.lon = lon_new
        p.lat = lat_new
        p.depth = depth_new

    def step_pair(self, pair: ParticlePairState, t: float) -> None:
        """
        Esegue un passo di integrazione per la coppia al tempo globale t (secondi).
        Dopo il passo, pair.t viene aggiornato a t + dt.
        """
        self._step_one_particle(pair.p1, t)
        self._step_one_particle(pair.p2, t)
        pair.t = t + self.dt
