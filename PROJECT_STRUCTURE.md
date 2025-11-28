# Project Structure

```
lagkinematic3.0/
├── python/
│   └── lagkinematic/
│       ├── __init__.py
│       ├── __main__.py
│       ├── chunking.py
│       ├── cli.py
│       ├── geometry.py
│       ├── integration.py
│       ├── utils.py
│       ├── core/
│       │   └── timegrid.py
│       ├── io/
│       │   ├── filespec.py
│       │   ├── regular_latlon.py
│       │   ├── starts.py
│       │   └── timeaxis.py
│       ├── sampling/
│       │   ├── __init__.py
│       │   ├── mask_regular_latlon.py
│       │   └── regular_latlon.py
│       └── writers/
│           ├── parquet_writer.py
│           ├── snaps_writer.py
│           ├── steps_parquet_writer.py
│           └── steps_writer.py
├── tests/
│   ├── test_mask_integration.py
│   └── test_mask_nan.py
├── examples/
│   ├── config.yml
│   ├── config_daterange.yml
│   ├── config_list.yml
│   ├── config_listfile.yml
│   ├── dataset_day01.nc
│   ├── dataset_day02.nc
│   ├── dataset_day03.nc
│   ├── elenco_dummy.txt
│   └── starts_dummy.dat
├── OUT_DEV_LIST/
│   ├── manifest.json
│   ├── trajectories.parquet
│   └── trajectories/
│       └── part-r0000-init.parquet
├── legacy_cpp/ (reference C++ implementation; read-only)
│   ├── CAMPI/
│   │   ├── 4DMED
│   │   ├── 4DMED_ANCHOVIES
│   │   ├── CGLORIS_1-12
│   │   ├── GLORAN_cglo
│   │   ├── GLORAN_cglo_2010-2019
│   │   ├── Mercator
│   │   └── RADAR
│   ├── NRc++/ (Numerical Recipes headers and utilities)
│   └── src/ (original Lagrangian model sources)
├── README.md
├── pyproject.toml
├── to_run.txt
└── lagkinematic3.0.7z
```
