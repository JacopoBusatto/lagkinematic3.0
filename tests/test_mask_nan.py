import numpy as np
import yaml
import pathlib

from lagkinematic.sampling.mask_regular_latlon import RegularLatLonMaskSampler
from lagkinematic.io.filespec import expand_files

def test_mask_nan():

    # Path robusto al config
    CONFIG_PATH = pathlib.Path(__file__).parent.parent / "examples" / "config_list.yml"

    with open(CONFIG_PATH) as f:
        cfg = yaml.safe_load(f)

    # Espandi i file come fa il CLI
    cfg["domain"]["_expanded_files"] = expand_files(cfg["run"], cfg["domain"]["files"])

    # Costruisci la maschera
    mask = RegularLatLonMaskSampler.from_first_file(cfg)

    # Scegli due punti (esempi)
    lon_sea, lat_sea = 12.0, 36.5
    lon_land, lat_land = 12.0, 44.0  # un punto sicuramente fuori mare
    depth = 0.0

    m_sea = mask.sample_mask(lon_sea, lat_sea, depth)
    m_land = mask.sample_mask(lon_land, lat_land, depth)

    print("SEA mask:", m_sea)
    print("LAND mask:", m_land)

    assert m_sea > 0.5
    assert m_land < 0.5


if __name__ == "__main__":
    test_mask_nan()
