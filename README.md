# 🌊 Lagkinematic 3.0  
**A modern, modular, Python-based Lagrangian kinematic model**

Lagkinematic 3.0 is a full Python rewrite of the original *LAGKINEMATIC 1.0* C++ model.  
The goal is to provide a **clean, modular, extensible and HPC-friendly** framework for Lagrangian simulations, with a focus on:

- oceanic and atmospheric velocity fields  
- pair dispersion  
- synthetic tracers  
- time-chunked execution  
- offline coupling with reanalyses (e.g., CMEMS/Mercator)

This version does **not** include any C++ code:  
🟦 *everything is being rewritten from scratch in Python*.

---

## 🚀 Features (current & planned)

### ✔ Implemented
- Modular package structure (`lagkinematic/`)
- Configurable CLI (`lagkinematic-cli` entry point)
- Time-chunked trajectory integration
- Snapshotting system for long simulations
- Writers for:
  - Parquet trajectory chunks
  - Per-step outputs (debug / diagnostics)
  - Snapshot saving
- Sampling of gridded velocity fields  
  (regular lat/lon via `RegularLatLonSampler`)
- Minimal configuration system (`config_min.yml`)

### 🔧 In progress
- Vector field abstraction
- Kinematic subgrid model
- Geometry utilities
- IO layer cleanup (grids, NetCDF, etc.)

### 🧪 Planned
- Full replacement of all C++ logic
- Parallel execution (Dask / multiprocessing)
- Integration with CMEMS, HYCOM, ECCO
- Random walk / turbulent diffusion
- Particle classes (buoyant, anchovy eggs, plume sources)

---

## 🗂 Repository Structure

lagkinematic3.0/
│
├── examples/ # Example configs & datasets
├── python/
│ └── lagkinematic/
│ ├── core/ # Core trajectory integrator
│ ├── io/ # Field loaders, readers
│ ├── sampling/ # Spatial/temporal interpolators
│ ├── writers/ # Output writers (parquet, snapshots...)
│ ├── geometry.py # Spherical distances, utilities
│ ├── integration.py # Time-stepping logic
│ ├── utils.py # Helpers
│ └── main.py # CLI entry point
│
├── config_min.yml # Minimal working configuration
├── pyproject.toml # Build system & dependencies
├── .gitignore
└── README.md


---

## ⚙️ Installation

### 🔹 Create a virtual environment

```bash
python -m venv .venv
source .venv/bin/activate      # Linux/macOS
.venv\Scripts\activate         # Windows

🔹 Install the package in editable mode
pip install -e .

▶️ Running a Simulation

Assuming you have config_min.yml in your folder:

lagkinematic-run config_min.yml


This will:

read input parameters

load time-chunked velocity fields

integrate trajectories step-by-step

save outputs under OUT_DEV/ or OUT_DEV_LIST/

🧪 Development

Format code:

black .


Lint:

ruff check .


Run tests (optional):

pytest tests/

📦 Packaging & Distribution

Build package:

pip install build
python -m build

🤝 Contributing

The project is under active development.
Feel free to open issues or feature requests, or propose improvements.

📜 License

MIT License (or another license — update as needed).

👤 Author

Jacopo Busatto
CNR-ISMAR
Oceanographic modelling & Lagrangian dynamics


---