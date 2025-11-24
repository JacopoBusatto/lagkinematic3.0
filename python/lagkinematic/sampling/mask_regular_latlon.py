import numpy as np
import xarray as xr

class RegularLatLonMaskSampler:
    """
    Maschera terra/mare basata sui NaN del campo U.
    - 1 = mare (U non è NaN)
    - 0 = terra (U è NaN)
    """

    def __init__(self, mask_da: xr.DataArray, lon_name: str, lat_name: str, depth_name: str | None):
        self.mask = mask_da
        self.lon_name = lon_name
        self.lat_name = lat_name
        self.depth_name = depth_name

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
        coords = {
            self.lon_name: lon,
            self.lat_name: lat,
        }

        if self.depth_name is not None:
            coords[self.depth_name] = depth

        out = self.mask.interp(coords, method="linear").fillna(0.0)
        values = out.values

        if np.isscalar(values):
            return float(values)
        return values.astype(float)
