# python/lagkinematic/io/regular_latlon.py
from __future__ import annotations
from pathlib import Path
import xarray as xr

def _vars_and_coords(ds):
    return sorted(list(ds.variables) + list(ds.coords))

def validate_regular_latlon(cfg: dict) -> dict:
    files = cfg["domain"]["_expanded_files"]
    if not files:
        raise RuntimeError("Nessun file NetCDF trovato dopo l'espansione.")

    first = Path(files[0])
    if not first.exists():
        raise FileNotFoundError(f"Primo file non trovato: {first}")

    coords = cfg["domain"]["coords"]
    vel    = cfg["domain"]["velocity"]

    with xr.open_dataset(first, decode_cf=False) as ds:
        avail = _vars_and_coords(ds)

        name_lon = coords["lon"]
        name_lat = coords["lat"]
        name_time = coords["time"]
        name_depth = coords.get("depth")
        name_u = vel["u"]["var"]
        name_v = vel["v"]["var"]
        name_w = (vel.get("w") or {}).get("var")

        for cname, label in [(name_lon,"lon"), (name_lat,"lat"), (name_time,"time")]:
            if cname not in ds.variables and cname not in ds.coords:
                raise ValueError(
                    f"Coordinata richiesta '{label}'='{cname}' non trovata in {first.name}.\n"
                    f"Disponibili: {avail}"
                )

        if name_w or name_depth:
            if not name_depth:
                raise ValueError("Dominio 3D richiesto ma 'coords.depth' non è definito nel config.")
            if name_depth not in ds.variables and name_depth not in ds.coords:
                raise ValueError(
                    f"Coordinata 'depth'='{name_depth}' non trovata in {first.name}.\n"
                    f"Disponibili: {avail}"
                )

        for vname, label in [(name_u,"u"), (name_v,"v")]:
            if vname not in ds.variables:
                raise ValueError(
                    f"Variabile velocità '{label}'='{vname}' non trovata in {first.name}.\n"
                    f"Disponibili: {avail}"
                )
        if name_w and name_w not in ds.variables:
            raise ValueError(
                f"Variabile 'w'='{name_w}' non trovata in {first.name}.\n"
                f"Disponibili: {avail}"
            )

        dims = {k: int(v) for k, v in ds.sizes.items()}
        dtypes = {"u": str(ds[name_u].dtype), "v": str(ds[name_v].dtype)}
        if name_w:
            dtypes["w"] = str(ds[name_w].dtype)

    return {"file": str(first), "dims": dims, "dtypes": dtypes}
