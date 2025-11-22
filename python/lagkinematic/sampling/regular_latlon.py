# python/lagkinematic/sampling/regular_latlon.py
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple, Dict, Any
import numpy as np
import xarray as xr

@dataclass
class FieldNames:
    lon: str
    lat: str
    depth: Optional[str]  # può essere None in 2D
    time: str
    u: str
    v: str
    w: Optional[str]      # opzionale

class RegularLatLonSampler:
    """
    Sampler per griglie regolari lat/lon:
    - Interpolazione temporale: lineare tra due snapshot adiacenti
    - Interpolazione spaziale: bilineare in XY (+ lineare in Z se presente)
    NOTE:
      * Minimizza I/O caricando solo 2 file alla volta (t_k, t_{k+1})
      * Assunzione: coordinate monotone crescenti
    """

    def __init__(self, files: list[str], names: FieldNames):
        self.files = [str(Path(f)) for f in files]
        self.nf = len(self.files)
        self.names = names

        # cache “doppia” di dataset: idx_k e idx_k+1
        self._ds_k: Optional[xr.Dataset] = None
        self._ds_kp1: Optional[xr.Dataset] = None
        self._ik: Optional[int] = None
        self._ikp1: Optional[int] = None

        # indici veloci (lon/lat/depth/time) per il file corrente
        self._lon = self._lat = self._depth = None
        self._time_k = self._time_kp1 = None

    @classmethod
    def from_config(cls, cfg: Dict[str, Any]) -> "RegularLatLonSampler":
        fcfg = cfg["domain"]
        files = fcfg["_expanded_files"]
        coords = fcfg["coords"]
        vel = fcfg["velocity"]
        names = FieldNames(
            lon=coords["lon"],
            lat=coords["lat"],
            depth=coords.get("depth"),
            time=coords["time"],
            u=vel["u"]["var"],
            v=vel["v"]["var"],
            w=(vel.get("w") or {}).get("var"),
        )
        return cls(files, names)

    def _open_ds(self, i: int) -> xr.Dataset:
        p = Path(self.files[i])
        if not p.exists():
            raise FileNotFoundError(f"File non trovato: {p}")
        # decode_times=True per avere datetime64; non carichiamo in memoria (lazy arrays)
        return xr.open_dataset(p, decode_times=True)

    def _ensure_time_bracket(self, t_ns: np.datetime64) -> Tuple[int, int, float]:
        """
        Trova due file consecutivi (k, k+1) tali che t_k <= t <= t_{k+1},
        e restituisce (k, k+1, alpha) con alpha in [0,1] per l'interpolazione temporale.
        Strategia: scorri file fino a trovare il bracketing usando il primo/ultimo time di ciascun file.
        """
        # quick pass: se abbiamo già un bracket valido, riusalo
        if self._ik is not None and self._ikp1 is not None and self._time_k is not None and self._time_kp1 is not None:
            if self._time_k[0] <= t_ns <= self._time_kp1[-1]:
                # alpha verrà calcolata internamente quando campioniamo (intra-file)
                return self._ik, self._ikp1, 0.0  # alpha placeholder

        # cerca bracketing k, k+1
        # apri sequenzialmente e usa min/max time del file
        last_ds = None
        for k in range(self.nf):
            ds = self._open_ds(k)
            t = ds[self.names.time].values.astype("datetime64[ns]")
            t0, t1 = t[0], t[-1]
            if t_ns < t0 and k == 0:
                # prima dell’inizio: clamp al primo file
                self._assign_pair(k, min(k + 1, self.nf - 1))
                ds.close() if hasattr(ds, "close") else None
                break
            if t0 <= t_ns <= t1:
                # dentro k — bracket (k,k) o (k,k+1) se esiste il successivo
                self._assign_pair(k, min(k + 1, self.nf - 1))
                ds.close() if hasattr(ds, "close") else None
                break
            last_ds = ds
        if self._ik is None:
            # oltre l'ultimo: clamp all’ultimo pair possibile
            self._assign_pair(self.nf - 2, self.nf - 1)

        # alpha non calcolata qui (dipende dall’indice temporale interno)
        return self._ik, self._ikp1, 0.0

    def _assign_pair(self, k: int, kp1: int) -> None:
        # chiudi vecchi
        for d in (self._ds_k, self._ds_kp1):
            try:
                d.close()
            except Exception:
                pass
        # apri nuovi
        self._ds_k = self._open_ds(k)
        self._ds_kp1 = self._open_ds(kp1)
        self._ik, self._ikp1 = k, kp1

        # cache coordinate
        self._lon = self._ds_k[self.names.lon].values
        self._lat = self._ds_k[self.names.lat].values
        self._depth = self._ds_k[self.names.depth].values if (self.names.depth and self.names.depth in self._ds_k) else None
        self._time_k = self._ds_k[self.names.time].values.astype("datetime64[ns]")
        self._time_kp1 = self._ds_kp1[self.names.time].values.astype("datetime64[ns]")

    # ----------------- INTERPOLAZIONI DI BASE -----------------

    @staticmethod
    def _find_bracketing_idx(vec: np.ndarray, x: float) -> Tuple[int, float]:
        """
        Ritorna (i, a) tali che x è tra vec[i] e vec[i+1] e a in [0,1] è la frazione.
        Clamp ai bordi se fuori dominio.
        """
        i = int(np.searchsorted(vec, x) - 1)
        if i < 0:
            return 0, 0.0
        if i >= len(vec) - 1:
            return len(vec) - 2, 1.0
        x0, x1 = float(vec[i]), float(vec[i + 1])
        a = 0.0 if x1 == x0 else (x - x0) / (x1 - x0)
        return i, float(np.clip(a, 0.0, 1.0))

    def _bilinear_xy(self, A: np.ndarray, lon: float, lat: float) -> float:
        """
        A: 2D array [lat, lon] o 3D con leading depth -> verrà indicizzato fuori.
        Bilineare semplice su griglia regolare.
        """
        j, ay = self._find_bracketing_idx(self._lat, lat)
        i, ax = self._find_bracketing_idx(self._lon, lon)
        # valori ai 4 vertici
        f00 = float(A[j    , i    ])
        f10 = float(A[j    , i + 1])
        f01 = float(A[j + 1, i    ])
        f11 = float(A[j + 1, i + 1])
        # bilineare
        return ( (1-ax)*(1-ay)*f00
               +    ax *(1-ay)*f10
               + (1-ax)*   ay *f01
               +    ax *   ay *f11 )

    def _linear_z(self, A3: np.ndarray, zvec: np.ndarray, z: float, lon: float, lat: float) -> float:
        """
        Interpolazione lineare in Z + bilineare XY.
        A3: 3D array [depth, lat, lon]
        """
        if zvec is None or zvec.size < 2:
            # fallback: niente z -> bilineare sul livello unico
            return self._bilinear_xy(A3[0, ...], lon, lat)
        k, az = self._find_bracketing_idx(zvec, z)
        v0 = self._bilinear_xy(A3[k    , ...], lon, lat)
        v1 = self._bilinear_xy(A3[k + 1, ...], lon, lat)
        return (1 - az) * v0 + az * v1

    # ----------------- API PUBBLICA -----------------

    def sample_uv(self, lon: float, lat: float, depth: float, t_ns: np.datetime64) -> Tuple[float, float, float]:
        """
        Restituisce (u,v,w) a (lon,lat,depth,t) con:
          - time-linear tra file k e k+1
          - bilinear in XY
          - linear in Z se depth è presente
        """
        k, kp1, _ = self._ensure_time_bracket(t_ns)
        assert self._ds_k is not None and self._ds_kp1 is not None

        # trova indici temporali interni ai due file e alpha temporale
        tk = self._time_k
        tkp1 = self._time_kp1
        # clamp t dentro il range di k..kp1
        t0 = tk[0]; t1 = tk[-1]; t2 = tkp1[-1]
        t_clamp = np.clip(t_ns, t0, t2)

        # se t è dentro k, usa bracketing intra-file (k)
        if t_clamp <= t1 and tk.size >= 2:
            idx, a = self._find_bracketing_idx(tk.astype("datetime64[ns]").astype("int64"), t_clamp.astype("int64"))
            # estrai slice temporali per u,v (2D o 3D)
            uA = self._ds_k[self.names.u].isel({self.names.time: idx}).values
            uB = self._ds_k[self.names.u].isel({self.names.time: idx+1}).values
            vA = self._ds_k[self.names.v].isel({self.names.time: idx}).values
            vB = self._ds_k[self.names.v].isel({self.names.time: idx+1}).values
            # W opzionale
            wA = wB = None
            if self.names.w and self.names.w in self._ds_k:
                wA = self._ds_k[self.names.w].isel({self.names.time: idx}).values
                wB = self._ds_k[self.names.w].isel({self.names.time: idx+1}).values
            # interp spazio
            if self._depth is not None and self.names.depth in self._ds_k[self.names.u].dims:
                u0 = self._linear_z(uA, self._depth, depth, lon, lat)
                u1 = self._linear_z(uB, self._depth, depth, lon, lat)
                v0 = self._linear_z(vA, self._depth, depth, lon, lat)
                v1 = self._linear_z(vB, self._depth, depth, lon, lat)
                if wA is not None:
                    w0 = self._linear_z(wA, self._depth, depth, lon, lat)
                    w1 = self._linear_z(wB, self._depth, depth, lon, lat)
                else:
                    w0 = w1 = 0.0
            else:
                u0 = self._bilinear_xy(uA, lon, lat); u1 = self._bilinear_xy(uB, lon, lat)
                v0 = self._bilinear_xy(vA, lon, lat); v1 = self._bilinear_xy(vB, lon, lat)
                w0 = w1 = 0.0 if wA is None else (
                    self._bilinear_xy(wA, lon, lat), self._bilinear_xy(wB, lon, lat)  # type: ignore
                )
            # time-linear
            u = (1 - a) * u0 + a * u1
            v = (1 - a) * v0 + a * v1
            w = (1 - a) * (w0 if isinstance(w0, float) else w0) + a * (w1 if isinstance(w1, float) else w1)
            return float(u), float(v), float(w)

        # altrimenti bracketing tra file k (ultimo tempo) e k+1 (primo tempo)
        uA = self._ds_k[self.names.u].isel({self.names.time: -1}).values
        uB = self._ds_kp1[self.names.u].isel({self.names.time: 0}).values
        vA = self._ds_k[self.names.v].isel({self.names.time: -1}).values
        vB = self._ds_kp1[self.names.v].isel({self.names.time: 0}).values
        wA = wB = None
        if self.names.w and (self.names.w in self._ds_k) and (self.names.w in self._ds_kp1):
            wA = self._ds_k[self.names.w].isel({self.names.time: -1}).values
            wB = self._ds_kp1[self.names.w].isel({self.names.time: 0}).values

        # alpha tra t1(k) e t0(k+1)
        t1 = self._time_k[-1].astype("int64")
        t2 = self._time_kp1[0].astype("int64")
        tt = t_ns.astype("int64")
        a = 0.0 if t2 == t1 else float(np.clip((tt - t1) / (t2 - t1), 0.0, 1.0))

        if self._depth is not None and self.names.depth in self._ds_k[self.names.u].dims:
            u0 = self._linear_z(uA, self._depth, depth, lon, lat)
            u1 = self._linear_z(uB, self._depth, depth, lon, lat)
            v0 = self._linear_z(vA, self._depth, depth, lon, lat)
            v1 = self._linear_z(vB, self._depth, depth, lon, lat)
            if wA is not None and wB is not None:
                w0 = self._linear_z(wA, self._depth, depth, lon, lat)
                w1 = self._linear_z(wB, self._depth, depth, lon, lat)
            else:
                w0 = w1 = 0.0
        else:
            u0 = self._bilinear_xy(uA, lon, lat); u1 = self._bilinear_xy(uB, lon, lat)
            v0 = self._bilinear_xy(vA, lon, lat); v1 = self._bilinear_xy(vB, lon, lat)
            w0 = w1 = 0.0 if wA is None else (
                self._bilinear_xy(wA, lon, lat), self._bilinear_xy(wB, lon, lat)  # type: ignore
            )

        u = (1 - a) * u0 + a * u1
        v = (1 - a) * v0 + a * v1
        w = (1 - a) * (w0 if isinstance(w0, float) else w0) + a * (w1 if isinstance(w1, float) else w1)
        return float(u), float(v), float(w)
