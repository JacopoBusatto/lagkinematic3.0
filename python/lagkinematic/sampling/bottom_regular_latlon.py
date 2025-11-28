from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict

import numpy as np
import xarray as xr


@dataclass
class RegularLatLonBottomSampler:
    """
    Stima la profondità locale H(x,y) su una griglia regolare lon/lat.

    - Se è presente un asse di profondità, usa il primo file NetCDF e il campo U
      per individuare il livello valido più profondo (non-NaN) per ciascuna cella.
    - In assenza dell'asse di profondità, restituisce 0 ovunque.

    Convenzione: profondità positive verso il basso (0 = superficie, >0 fondale).
    """

    bottom: xr.DataArray
    lon_name: str
    lat_name: str
    _values: np.ndarray | None = None
    _lon: np.ndarray | None = None
    _lat: np.ndarray | None = None

    @classmethod
    def from_first_file(cls, cfg: Dict[str, Any]) -> "RegularLatLonBottomSampler":
        file0 = cfg["domain"]["_expanded_files"][0]

        coords = cfg["domain"]["coords"]
        lon_name = coords["lon"]
        lat_name = coords["lat"]
        depth_name = coords.get("depth")
        time_name = coords.get("time")

        uvar = cfg["domain"]["velocity"]["u"]["var"]

        with xr.open_dataset(file0, decode_times=False) as ds:
            da_u = ds[uvar]

            if time_name and (time_name in da_u.dims):
                da_u = da_u.isel({time_name: 0})

            if depth_name and (depth_name in da_u.dims):
                da_u = da_u.transpose(depth_name, lat_name, lon_name)
                depth_vals = ds[depth_name].values
                data = da_u.values
                valid = np.isfinite(data)
                depth_grid = np.broadcast_to(depth_vals.reshape(-1, 1, 1), data.shape)
                bottom_vals = np.where(valid, depth_grid, -np.inf).max(axis=0)
                bottom_vals[~np.isfinite(bottom_vals)] = 0.0
            else:
                # nessun asse di profondità → profondità nulla ovunque
                shape_lat = ds[lat_name].shape
                shape_lon = ds[lon_name].shape
                bottom_vals = np.zeros((shape_lat[0], shape_lon[0]), dtype=float)

            lat_vals = ds[lat_name].values
            lon_vals = ds[lon_name].values

        bottom_da = xr.DataArray(
            bottom_vals,
            coords={lat_name: lat_vals, lon_name: lon_vals},
            dims=(lat_name, lon_name),
        )

        sampler = cls(bottom=bottom_da, lon_name=lon_name, lat_name=lat_name)

        bottom_std = bottom_da.transpose(lat_name, lon_name)
        sampler._values = np.nan_to_num(bottom_std.values, nan=0.0)
        sampler._lon = bottom_std[lon_name].values
        sampler._lat = bottom_std[lat_name].values

        return sampler

    def sample_bottom(self, lon: float, lat: float) -> float:
        if self._lon is None or self._lat is None or self._values is None:
            return float(self.bottom)

        if (lon < self._lon[0]) or (lon > self._lon[-1]) or (lat < self._lat[0]) or (lat > self._lat[-1]):
            return 0.0

        def _find_bracketing_idx(vec: np.ndarray, x: float) -> tuple[int, float]:
            i = int(np.searchsorted(vec, x) - 1)
            if i < 0:
                return 0, 0.0
            if i >= len(vec) - 1:
                return len(vec) - 2, 1.0
            x0, x1 = float(vec[i]), float(vec[i + 1])
            a = 0.0 if x1 == x0 else (x - x0) / (x1 - x0)
            return i, float(np.clip(a, 0.0, 1.0))

        def _bilinear_xy(A: np.ndarray, lon_val: float, lat_val: float) -> float:
            j, ay = _find_bracketing_idx(self._lat, lat_val)
            i, ax = _find_bracketing_idx(self._lon, lon_val)
            f00 = float(A[j, i])
            f10 = float(A[j, i + 1])
            f01 = float(A[j + 1, i])
            f11 = float(A[j + 1, i + 1])
            return (
                (1 - ax) * (1 - ay) * f00
                + ax * (1 - ay) * f10
                + (1 - ax) * ay * f01
                + ax * ay * f11
            )

        return float(_bilinear_xy(self._values, lon, lat))
