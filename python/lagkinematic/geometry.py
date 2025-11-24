# python/lagkinematic/geometry.py
from __future__ import annotations

import math
from dataclasses import dataclass

METERS_PER_DEGREE_LAT = 111_039.0  # ~metri per grado di latitudine (coerente con legacy C++)

@dataclass(frozen=True)
class LonLatDomain:
    """
    Dominio geografico in longitudine e latitudine.
    - lon_min, lon_max permettono wrapping periodico solo in lon
    - lat_min, lat_max sono sempre rigidi (nessun wrapping ai poli)
    """
    lon_min: float
    lon_max: float
    lat_min: float
    lat_max: float

    @property
    def width(self) -> float:
        return self.lon_max - self.lon_min



def displace_lonlat(
    lon_deg: float,
    lat_deg: float,
    dx_m: float,
    dy_m: float,
) -> tuple[float, float]:
    """
    Sposta un punto (lon, lat) di (dx, dy) in metri usando la convenzione legacy:

        dlat = dy / 111039
        dlon = dx / (111039 * cos(lat * pi/180))

    Restituisce (lon_new, lat_new) in gradi, SENZA wrapping periodico.
    """
    # gradi di latitudine
    dlat = dy_m / METERS_PER_DEGREE_LAT

    # metri per grado di longitudine, dipende dalla latitudine
    meters_per_degree_lon = METERS_PER_DEGREE_LAT * math.cos(math.radians(lat_deg))

    # per lat ≈ ±90° meters_per_degree_lon → 0; protezione minimale
    if abs(meters_per_degree_lon) < 1e-6:
        dlon = 0.0
    else:
        dlon = dx_m / meters_per_degree_lon

    return lon_deg + dlon, lat_deg + dlat


def wrap_longitude(lon: float, domain: LonLatDomain) -> float:
    """
    Mappa la longitudine nel range [lon_min, lon_max) usando condizioni periodiche.

    È indipendente dalla risoluzione della griglia, funziona per:
      - [-180, 180)
      - [0, 360)
      - qualunque altro intervallo continuo di ampiezza 360 (o multipli).
    """
    return ((lon - domain.lon_min) % domain.width) + domain.lon_min
