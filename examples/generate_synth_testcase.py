"""
Genera un dataset sintetico per testare LAGKINEMATIC 3.0.

Crea:
- 10 NetCDF giornalieri in examples/synth/dataset_synth_dayXX.nc
- una griglia RegularLatLon con:
    lon: 10–20E, step 0.5°
    lat: 35–40N, step 0.5°
    depth: 0,1,2,3,4 m (profondità positiva verso il basso)
- velocità (ugo, vgo) diagonali, decrescenti con la profondità
- una "barriera" di terra dove u=v=NaN (usata come mask)
- un file starts in examples/starts_synth.dat con:
    - 1 coppia che parte su terra (deve essere uccisa subito dal mask)
    - 4 coppie alla stessa posizione orizzontale ma profondità 0–4 m
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import xarray as xr


def make_grid():
    lon = np.arange(10.0, 20.0 + 1e-6, 0.5)   # 10 → 20 deg E
    lat = np.arange(35.0, 40.0 + 1e-6, 0.5)   # 35 → 40 deg N
    depth = np.arange(0.0, 5.0, 1.0)          # 0,1,2,3,4 m

    return lon, lat, depth


def make_velocity_field(lon, lat, depth, time_value):
    """
    Costruisce ugo, vgo con:
      - direzione diagonale (u>0, v>0)
      - modulo che decresce linearmente con la profondità
      - barriera "terra" con NaN in una fascia di coordinate
    """
    # dimensioni
    nlon = lon.size
    nlat = lat.size
    ndep = depth.size

    # velocità di riferimento alla superficie (m/s)
    u0 = 0.10  # ~ 0.1 m/s verso Est
    v0 = 0.05  # ~ 0.05 m/s verso Nord

    # fattore di decadimento con profondità: 1 a z=0, 0 a z=5
    depth_factor = 1.0 - depth / 5.0  # shape (ndep,)
    depth_factor = np.clip(depth_factor, 0.0, 1.0)

    # broadcasting a (depth, lat, lon)
    u3d = np.zeros((ndep, nlat, nlon), dtype=np.float64)
    v3d = np.zeros((ndep, nlat, nlon), dtype=np.float64)

    for k in range(ndep):
        f = depth_factor[k]
        u3d[k, :, :] = u0 * f
        v3d[k, :, :] = v0 * f

    # Barriera "terra": mettiamo NaN in una fascia
    # Esempio: lon >= 16 e 37 <= lat <= 38
    lon2d, lat2d = np.meshgrid(lon, lat)
    land_mask = (lon2d >= 16.0) & (lat2d >= 37.0) & (lat2d <= 38.0)

    # applica la maschera a tutti i livelli
    for k in range(ndep):
        u3d[k, land_mask] = np.nan
        v3d[k, land_mask] = np.nan

    # aggiungi dimensione time (1)
    u4d = u3d[None, ...]  # (time, depth, lat, lon)
    v4d = v3d[None, ...]

    ds = xr.Dataset(
        data_vars={
            "ugo": (("time", "depth", "lat", "lon"), u4d),
            "vgo": (("time", "depth", "lat", "lon"), v4d),
        },
        coords={
            "time": np.array([np.datetime64(time_value)]),
            "depth": depth,
            "lat": lat,
            "lon": lon,
        },
    )

    return ds


def generate_netcdf_series(outdir: Path, n_days: int = 10):
    outdir.mkdir(parents=True, exist_ok=True)
    lon, lat, depth = make_grid()

    base_date = np.datetime64("2016-06-01")
    paths = []

    for i in range(n_days):
        time_value = base_date + np.timedelta64(i, "D")
        ds = make_velocity_field(lon, lat, depth, time_value)

        fname = outdir / f"dataset_synth_day{i+1:02d}.nc"
        ds.to_netcdf(fname)
        paths.append(fname)

        print(f"[synth] scritto {fname}")

    return paths, (lon, lat, depth)


def generate_starts(file_path: Path, lon, lat, depth):
    """
    Crea uno starts_synth.dat con 5 coppie:

    Formato colonne:
        lon1, lat1, dep1, lon2, lat2, dep2, t_delay, col8
    """

    rows = []

    # 1) Coppia su terra (nella fascia di NaN): deve essere uccisa subito dal mask
    # scegliamo lon=16.5, lat=37.5, profondità 0m
    rows.append(
        [16.5, 37.5, 0.0, 16.6, 37.5, 0.0, 0.0, 0.0]
    )

    # 2–5) Coppie alla stessa posizione orizzontale ma profondità 0–4 m
    # posizione in mare: longe ~ 12, late ~ 36
    lon_ocean = 12.0
    lat_ocean = 36.0
    for dep in [0.0, 1.0, 2.0, 4.0]:
        rows.append(
            [lon_ocean, lat_ocean, dep, lon_ocean + 0.1, lat_ocean, dep, 0.0, 0.0]
        )

    data = np.array(rows, dtype=np.float64)
    file_path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(file_path, data, fmt="%.6f")

    print(f"[synth] scritto starts → {file_path}")
    print(f"[synth] n_coppie = {data.shape[0]}")


def main():
    root = Path(__file__).resolve().parent.parent  # cartella Lagkinematic3.0
    synth_dir = root / "examples" / "synth"
    starts_path = root / "examples" / "starts_synth.dat"

    print(f"[synth] root = {root}")
    print(f"[synth] outdir netcdf = {synth_dir}")

    nc_paths, (lon, lat, depth) = generate_netcdf_series(synth_dir, n_days=10)
    generate_starts(starts_path, lon, lat, depth)

    print("[synth] Fatto.")
    print("[synth] NetCDF generati:")
    for p in nc_paths:
        print(f"   - {p}")


if __name__ == "__main__":
    main()
