import numpy as np
import xarray as xr

class RegularLatLonMaskSampler:
    """
    Maschera terra/mare basata sui NaN del campo U.
    - 1 = mare (U non è NaN)
    - 0 = terra (U è NaN)

    depth segue la convenzione positiva verso il basso (0 = superficie, >0 sotto).
    L'uso della maschera nell'integratore può essere controllato con
    beaching_mode: "off" (ignora la maschera), "kill" (uccidi) o "bounce"
    (riflessione).
    """

    def __init__(self, mask_da: xr.DataArray, lon_name: str, lat_name: str, depth_name: str | None):
        self.mask = mask_da
        self.lon_name = lon_name
        self.lat_name = lat_name
        self.depth_name = depth_name

        # cache coordinate e valori per interp manuale
        dims_order = [lat_name, lon_name] if depth_name is None else [depth_name, lat_name, lon_name]
        mask_std = mask_da.transpose(*dims_order)
        self._values = np.nan_to_num(mask_std.values, nan=0.0)
        self._lon = mask_std[lon_name].values
        self._lat = mask_std[lat_name].values
        self._depth = mask_std[depth_name].values if depth_name is not None else None

    @classmethod
    def from_first_file(cls, cfg):
        """
        Legge il PRIMO file NetCDF del dominio e costruisce la maschera.
        """
        file0 = cfg["domain"]["_expanded_files"][0]

        coords = cfg["domain"]["coords"]
        lon_name = coords["lon"]
        lat_name = coords["lat"]
        depth_name = coords.get("depth")
        time_name = coords.get("time")

        uvar = cfg["domain"]["velocity"]["u"]["var"]

        with xr.open_dataset(file0, decode_times=False) as ds:
            da_u = ds[uvar]

            # estrai time=0 se serve
            if time_name and (time_name in da_u.dims):
                da_u = da_u.isel({time_name: 0})

            # maschera 3D (depth, lat, lon)
            if depth_name and (depth_name in da_u.dims):
                mask = xr.where(da_u.notnull(), 1.0, 0.0)
            else:
                # maschera 2D (lat, lon)
                mask2D = xr.where(da_u.notnull(), 1.0, 0.0)
                mask = mask2D

        return cls(mask_da=mask, lon_name=lon_name, lat_name=lat_name, depth_name=depth_name)

    # -----------------------------------------------------------------------------

    def sample_mask(self, lon, lat, depth=0.0):
        """
        Restituisce la maschera interpolata in (lon,lat[,depth]).
        Ritorna 0.0 fuori dominio.
        """
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

        if self._depth is None:
            val = _bilinear_xy(self._values, lon, lat)
        elif self._values.shape[0] == 1:
            val = _bilinear_xy(self._values[0, ...], lon, lat)
        else:
            k, az = _find_bracketing_idx(self._depth, depth)
            v0 = _bilinear_xy(self._values[k, ...], lon, lat)
            v1 = _bilinear_xy(self._values[k + 1, ...], lon, lat)
            val = (1 - az) * v0 + az * v1

        return float(val)
