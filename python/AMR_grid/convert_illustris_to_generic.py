#!/usr/bin/env python3
"""
convert_illustris_to_generic.py
================================

Convert Illustris/IllustrisTNG snapshot or cutout data into the generic AMR
format understood by LaRT v2.00+.

Illustris/TNG uses an unstructured Voronoi mesh (Arepo code).  LaRT requires
an octree AMR grid.  This converter:

  1. Reads PartType0 (gas) cells from an HDF5 snapshot or API cutout.
  2. Converts Illustris comoving+h code units to physical (kpc, cm^-3, K, km/s).
  3. Builds an adaptive octree via ``AMRGrid`` with density/velocity gradient
     refinement and optional resolution matching to the local Voronoi cell size.
  4. Assigns leaf-cell properties via nearest-neighbor KD-tree lookup
     (physically exact for Voronoi tessellation).
  5. Optionally computes physics columns (xHI, n_e, emissivity, ndust).
  6. Writes the generic AMR output (HDF5, FITS, or text).

Dependencies
------------
- numpy, scipy (cKDTree), h5py
- requests (optional, for API download mode)
- astropy  (optional, for FITS output)

Usage
-----
Local file::

    python convert_illustris_to_generic.py cutout_99.hdf5 -o galaxy.h5 \\
        --output-unit kpc --level-max 8 --compute-physics

API download::

    python convert_illustris_to_generic.py \\
        --api-key YOUR_KEY --simulation TNG50-1 --snap 99 --subhalo-id 0 \\
        -o galaxy.h5 --compute-physics
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
MSUN_CGS   = 1.989e33                    # g
KPC_CM     = 3.0856775814913673e21       # cm
PC_CM      = 3.0856775814913673e18       # cm
MASS_H_CGS = 1.6726e-24                  # g  (hydrogen atom mass)
MASS_P_CGS = 1.6726e-24                  # g  (proton mass, ~m_H)
KB_CGS     = 1.3807e-16                  # erg K^-1
X_H        = 0.76                        # primordial hydrogen mass fraction
GAMMA      = 5.0 / 3.0                   # adiabatic index

# Illustris internal unit: 10^10 Msun/h
ILLUSTRIS_MASS_CGS = 1.0e10 * MSUN_CGS   # g / h


# ---------------------------------------------------------------------------
# Physics functions (imported from convert_ramses_to_generic if available,
#                     otherwise defined here as fallback)
# ---------------------------------------------------------------------------
try:
    from convert_ramses_to_generic import (
        cie_neutral_fraction,
        cie_neutral_fraction_table,
        electron_density,
        laursen09_ndust,
        caseB_lya_emissivity,
        cie_civ_emissivity,
        cie_ovi_emissivity,
    )
except ImportError:
    def cie_neutral_fraction(T):
        """CIE neutral fraction (simple formula)."""
        T = np.maximum(T, 10.0)
        T4 = T / 1.0e4
        k_ion = 5.84862e-9 * np.sqrt(T4) * np.exp(-15.78215 / T4)
        k_rec = 4.13e-13 * T4**(-0.7131 - 0.0115 * np.log(T4))
        return k_rec / (k_ion + k_rec)

    def cie_neutral_fraction_table(nH, T):
        return cie_neutral_fraction(T)

    def electron_density(nH, xHI):
        return nH * (1.0 - xHI)

    def laursen09_ndust(nH, xHI, Z, Z_ref=0.0134, f_ion=0.01):
        nHI  = nH * xHI
        nHII = nH * (1.0 - xHI)
        return (Z / np.maximum(Z_ref, 1e-30)) * (nHI + f_ion * nHII)

    def caseB_lya_emissivity(nH, T, xHI, ne):
        T = np.maximum(T, 10.0)
        lam = 315614.0 / T
        aB = 2.753e-14 * lam**1.5 / (1.0 + (lam / 2.74)**0.407)**2.242
        Ta = np.maximum(T, 100.0)
        PB = 0.686 - 0.106 * np.log10(Ta / 1e4) - 0.009 * (Ta / 1e4)**(-0.44)
        nHI  = nH * xHI
        nHII = nH * (1.0 - xHI)
        e_rec = PB * aB * ne * nHII
        qc = (6.58e-18 / T**0.185) * np.exp(-4.86e4 / T**0.895)
        e_col = nHI * ne * qc
        return e_rec + e_col

    # CIE metal-line emissivity (fallback)
    _K_BOLTZMANN_EV = 8.6173e-5
    def _cie_metal_emiss(nH, T, ne, Z, Z_ref, A_X, fit, Omega, DE_eV):
        T = np.maximum(np.asarray(T, dtype=float), 10.0)
        logT = np.log10(T)
        f_ion = np.clip(fit['f_peak'] * np.exp(
            -0.5 * ((logT - fit['logT_peak']) / fit['sigma'])**2), 0, 1)
        n_ion = nH * (Z / np.maximum(Z_ref, 1e-30)) * A_X * f_ion
        q = (8.629e-6 / (2.0 * np.sqrt(T))) * Omega * np.exp(
            -DE_eV / (_K_BOLTZMANN_EV * T))
        return n_ion * ne * q

    def cie_civ_emissivity(nH, T, ne, Z, Z_ref=0.0134):
        return _cie_metal_emiss(nH, T, ne, Z, Z_ref, 2.692e-4,
            dict(logT_peak=5.05, f_peak=0.29, sigma=0.20), 7.5, 8.00)

    def cie_ovi_emissivity(nH, T, ne, Z, Z_ref=0.0134):
        return _cie_metal_emiss(nH, T, ne, Z, Z_ref, 4.898e-4,
            dict(logT_peak=5.45, f_peak=0.20, sigma=0.18), 3.8, 11.98)


# ---------------------------------------------------------------------------
# Illustris data container
# ---------------------------------------------------------------------------
@dataclass
class IllustrisData:
    """Container for Illustris/TNG gas cell data in physical units."""
    # Positions [kpc, physical, centered]
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    # Physical quantities
    nH: np.ndarray             # cm^-3
    T: np.ndarray              # K
    vx: np.ndarray             # km/s
    vy: np.ndarray             # km/s
    vz: np.ndarray             # km/s
    # Optional Illustris-native quantities
    metallicity: Optional[np.ndarray] = None    # mass fraction
    xHI: Optional[np.ndarray] = None            # NeutralHydrogenAbundance
    electron_abundance: Optional[np.ndarray] = None
    sfr: Optional[np.ndarray] = None            # Msun/yr
    volume_kpc3: Optional[np.ndarray] = None    # kpc^3
    # Header
    hubble_param: float = 0.6774
    redshift: float = 0.0
    scale_factor: float = 1.0
    boxsize_ckpc_h: float = 0.0   # original comoving boxsize


# ---------------------------------------------------------------------------
# Data loader — local HDF5
# ---------------------------------------------------------------------------
def _read_header(f) -> dict:
    """Read Illustris/TNG HDF5 header attributes."""
    hdr = {}
    h = f["Header"]
    for key in h.attrs:
        val = h.attrs[key]
        hdr[key] = val
    return hdr


def load_illustris_hdf5(filepath: str,
                        center_kpc: Optional[np.ndarray] = None,
                        boxsize_kpc: Optional[float] = None) -> IllustrisData:
    """
    Load gas cell data from an Illustris/TNG HDF5 file (cutout or snapshot chunk).

    Parameters
    ----------
    filepath : str
        Path to HDF5 file.
    center_kpc : (3,) array or None
        If given, recenter coordinates so this point becomes the origin [physical kpc].
    boxsize_kpc : float or None
        If given, override the box size for coordinate wrapping [physical kpc].
        Only needed for full periodic snapshots when recentering.
    """
    import h5py

    with h5py.File(filepath, "r") as f:
        hdr = _read_header(f)

        h = float(hdr.get("HubbleParam", 0.6774))
        z = float(hdr.get("Redshift", 0.0))
        a = float(hdr.get("Time", 1.0 / (1.0 + z)))
        boxsize_ckpc_h = float(hdr.get("BoxSize", 0.0))

        pt0 = f["PartType0"]

        # --- Raw data ---
        coords   = pt0["Coordinates"][:]          # (N, 3) ckpc/h
        density  = pt0["Density"][:]               # (N,)   (10^10 Msun/h)/(ckpc/h)^3
        int_en   = pt0["InternalEnergy"][:]        # (N,)   (km/s)^2
        vel      = pt0["Velocities"][:]            # (N, 3) km*sqrt(a)/s

        # Optional fields
        metallicity = pt0["GFM_Metallicity"][:] if "GFM_Metallicity" in pt0 else None
        xHI_native  = pt0["NeutralHydrogenAbundance"][:] if "NeutralHydrogenAbundance" in pt0 else None
        xe          = pt0["ElectronAbundance"][:] if "ElectronAbundance" in pt0 else None
        sfr_arr     = pt0["StarFormationRate"][:] if "StarFormationRate" in pt0 else None
        vol_raw     = pt0["Volume"][:] if "Volume" in pt0 else None

        # If Volume not available, estimate from Masses/Density
        if vol_raw is None and "Masses" in pt0:
            masses = pt0["Masses"][:]   # 10^10 Msun/h
            vol_raw = masses / np.maximum(density, 1e-30)  # (ckpc/h)^3

    # --- Unit conversions ---
    # Positions: ckpc/h → kpc (physical)
    x_kpc = coords[:, 0] * a / h
    y_kpc = coords[:, 1] * a / h
    z_kpc = coords[:, 2] * a / h

    # Recenter
    if center_kpc is not None:
        cx, cy, cz = center_kpc
        x_kpc -= cx
        y_kpc -= cy
        z_kpc -= cz
        # Periodic wrapping (only meaningful for full snapshots)
        if boxsize_kpc is not None:
            half = boxsize_kpc / 2.0
            x_kpc = np.where(x_kpc >  half, x_kpc - boxsize_kpc, x_kpc)
            x_kpc = np.where(x_kpc < -half, x_kpc + boxsize_kpc, x_kpc)
            y_kpc = np.where(y_kpc >  half, y_kpc - boxsize_kpc, y_kpc)
            y_kpc = np.where(y_kpc < -half, y_kpc + boxsize_kpc, y_kpc)
            z_kpc = np.where(z_kpc >  half, z_kpc - boxsize_kpc, z_kpc)
            z_kpc = np.where(z_kpc < -half, z_kpc + boxsize_kpc, z_kpc)

    # Density: (10^10 Msun/h) / (ckpc/h)^3 → g/cm^3 → nH [cm^-3]
    rho_cgs = density * (ILLUSTRIS_MASS_CGS * h**2) / (KPC_CM * a)**3
    nH = rho_cgs * X_H / MASS_H_CGS

    # Temperature: from InternalEnergy + ElectronAbundance
    if xe is not None:
        mu = 4.0 / (1.0 + 3.0 * X_H + 4.0 * X_H * xe) * MASS_P_CGS
    else:
        # Assume fully neutral as fallback
        mu = 4.0 / (1.0 + 3.0 * X_H) * MASS_P_CGS
    # InternalEnergy is in (km/s)^2 = 1e10 erg/g
    T = (GAMMA - 1.0) * int_en * 1.0e10 * mu / KB_CGS

    # Velocities: km*sqrt(a)/s → km/s (physical peculiar)
    vx_kms = vel[:, 0] * np.sqrt(a)
    vy_kms = vel[:, 1] * np.sqrt(a)
    vz_kms = vel[:, 2] * np.sqrt(a)

    # Volume: (ckpc/h)^3 → kpc^3 (physical)
    vol_kpc3 = vol_raw * (a / h)**3 if vol_raw is not None else None

    return IllustrisData(
        x=x_kpc, y=y_kpc, z=z_kpc,
        nH=nH, T=T, vx=vx_kms, vy=vy_kms, vz=vz_kms,
        metallicity=metallicity,
        xHI=xHI_native,
        electron_abundance=xe,
        sfr=sfr_arr,
        volume_kpc3=vol_kpc3,
        hubble_param=h,
        redshift=z,
        scale_factor=a,
        boxsize_ckpc_h=boxsize_ckpc_h,
    )


# ---------------------------------------------------------------------------
# API download (optional)
# ---------------------------------------------------------------------------
def download_tng_cutout(simulation: str, snap_num: int, subhalo_id: int,
                        api_key: str, cache_dir: str = ".") -> Path:
    """
    Download a subhalo cutout HDF5 from the TNG public API.

    Parameters
    ----------
    simulation : str
        Simulation name, e.g. ``'TNG50-1'``, ``'TNG100-1'``.
    snap_num : int
        Snapshot number (e.g. 99 for z=0).
    subhalo_id : int
        Subhalo ID within the snapshot.
    api_key : str
        TNG API key (register at tng-project.org).
    cache_dir : str
        Directory to save the downloaded file.

    Returns
    -------
    Path to the downloaded HDF5 file.
    """
    try:
        import requests
    except ImportError:
        raise ImportError("API download requires the 'requests' package.  "
                          "Install with: pip install requests")

    base_url = "https://www.tng-project.org/api"
    headers = {"api-key": api_key}

    # Get subhalo metadata to find cutout URL
    sub_url = f"{base_url}/{simulation}/snapshots/{snap_num}/subhalos/{subhalo_id}"
    r = requests.get(sub_url, headers=headers)
    r.raise_for_status()
    meta = r.json()

    # Request cutout with gas fields
    cutout_url = meta["cutouts"]["subhalo"]
    params = {
        "gas": "Coordinates,Density,InternalEnergy,Velocities,"
               "GFM_Metallicity,NeutralHydrogenAbundance,"
               "ElectronAbundance,Masses,StarFormationRate"
    }
    r = requests.get(cutout_url, headers=headers, params=params)
    r.raise_for_status()

    outpath = Path(cache_dir) / f"cutout_{simulation}_snap{snap_num}_sub{subhalo_id}.hdf5"
    outpath.write_bytes(r.content)
    print(f"Downloaded cutout: {outpath}  ({outpath.stat().st_size / 1e6:.1f} MB)")
    return outpath


# ---------------------------------------------------------------------------
# VoronoiInterpolator — nearest-neighbor from Voronoi cells
# ---------------------------------------------------------------------------
class VoronoiInterpolator:
    """
    Resample Voronoi cell data onto query points using a KD-tree.

    Two methods are supported:

    ``method='nearest'`` (default)
        Each query point takes the value of the nearest Voronoi generator.
        Physically exact for a Voronoi tessellation (each point in space
        belongs to exactly one cell), but produces blocky / staircase fields.

    ``method='gaussian'``
        Volume-weighted, normalized Gaussian gather over the K nearest cells::

            A(r) = Sum_i V_i A_i exp(-d_i^2 / 2 h_i^2)
                 / Sum_i V_i      exp(-d_i^2 / 2 h_i^2)

        The smoothing length is adaptive by default, ``h_i = factor * r_eff,i``
        with ``r_eff,i = (3 V_i / 4 pi)^(1/3)`` the effective radius of cell i,
        so the kernel tracks the local Voronoi resolution.  A fixed ``h`` (in
        the coordinate unit) can be forced with ``kernel_size``.
    """

    def __init__(self, positions: np.ndarray, data: IllustrisData,
                 method: str = "nearest", kernel_factor: float = 1.0,
                 kernel_neighbors: int = 32, kernel_size=None):
        """
        Parameters
        ----------
        positions : (N, 3) array
            Voronoi generator positions [kpc].
        data : IllustrisData
            Container with all physical arrays (indexed 1:1 with positions).
        method : {'nearest', 'gaussian'}
            Value-assignment scheme (see class docstring).
        kernel_factor : float
            Adaptive smoothing length in units of the cell effective radius
            (``h_i = kernel_factor * r_eff,i``).  Used when ``kernel_size`` is
            None.
        kernel_neighbors : int
            Number of nearest cells gathered per query point (gaussian only).
        kernel_size : float or None
            If set (> 0), a fixed smoothing length for all points, overriding
            the adaptive ``kernel_factor * r_eff`` rule.
        """
        from scipy.spatial import cKDTree
        self.tree = cKDTree(positions)
        self.data = data
        self.positions = positions
        self.method = method
        self.kernel_factor = float(kernel_factor)
        self.kernel_neighbors = int(kernel_neighbors)
        self.kernel_size = kernel_size
        self._wcache = None   # (key, idx, wnorm) -- reuse weights across fields

        if method == "gaussian":
            n = len(data.x)
            self.kernel_neighbors = max(1, min(self.kernel_neighbors, n))
            # Effective cell radius and volume weights.
            if data.volume_kpc3 is not None:
                vol = np.maximum(np.asarray(data.volume_kpc3, dtype=float), 1e-30)
            else:
                # Fallback: estimate cell size from nearest-neighbor spacing.
                dnn, _ = self.tree.query(positions, k=2)
                r_nn = np.maximum(dnn[:, 1], 1e-30)
                vol = (4.0 / 3.0) * np.pi * r_nn**3
            self.cell_volume = vol
            self.r_eff = (3.0 * vol / (4.0 * np.pi))**(1.0 / 3.0)

    def query(self, x, y, z):
        """Return indices of nearest Voronoi cells for query points."""
        pts = np.column_stack([np.atleast_1d(x),
                               np.atleast_1d(y),
                               np.atleast_1d(z)])
        _, idx = self.tree.query(pts)
        return idx

    def _gather_weights(self, x, y, z):
        """Return (idx (P,K), wnorm (P,K)) normalized Gaussian kernel weights."""
        x = np.atleast_1d(x); y = np.atleast_1d(y); z = np.atleast_1d(z)
        key = (id(x), id(y), id(z), x.shape[0])
        if self._wcache is not None and self._wcache[0] == key:
            return self._wcache[1], self._wcache[2]
        pts = np.column_stack([x, y, z])
        K = self.kernel_neighbors
        d, idx = self.tree.query(pts, k=K)
        if K == 1:
            d = d[:, None]; idx = idx[:, None]
        if self.kernel_size is not None and self.kernel_size > 0:
            h = float(self.kernel_size)              # fixed, scalar
        else:
            h = self.kernel_factor * self.r_eff[idx]  # adaptive, per-neighbor
        w = self.cell_volume[idx] * np.exp(-0.5 * (d / h)**2)
        wsum = w.sum(axis=1)
        zero = (wsum == 0.0)
        safe = np.where(zero, 1.0, wsum)
        wnorm = w / safe[:, None]
        if zero.any():                                # fall back to nearest cell
            wnorm[zero] = 0.0
            wnorm[zero, 0] = 1.0
        self._wcache = (key, idx, wnorm)
        return idx, wnorm

    def sample(self, field, x, y, z):
        """Resample a scalar ``field`` (indexed 1:1 with cells) at query points."""
        if self.method == "gaussian":
            idx, wnorm = self._gather_weights(x, y, z)
            return (wnorm * field[idx]).sum(axis=1)
        return field[self.query(x, y, z)]

    def density_fn(self, x, y, z):
        """Return nH [cm^-3] at query points."""
        return self.sample(self.data.nH, x, y, z)

    def temperature_fn(self, x, y, z):
        """Return T [K] at query points."""
        return self.sample(self.data.T, x, y, z)

    def velocity_fn(self, x, y, z):
        """Return (vx, vy, vz) [km/s] at query points."""
        if self.method == "gaussian":
            idx, wnorm = self._gather_weights(x, y, z)
            return ((wnorm * self.data.vx[idx]).sum(axis=1),
                    (wnorm * self.data.vy[idx]).sum(axis=1),
                    (wnorm * self.data.vz[idx]).sum(axis=1))
        idx = self.query(x, y, z)
        return self.data.vx[idx], self.data.vy[idx], self.data.vz[idx]

    def properties_fn(self, x, y, z):
        """Return (nH, T, vx, vy, vz) at query points (for set_properties)."""
        if self.method == "gaussian":
            idx, wnorm = self._gather_weights(x, y, z)
            g = lambda f: (wnorm * f[idx]).sum(axis=1)
            return (g(self.data.nH), g(self.data.T),
                    g(self.data.vx), g(self.data.vy), g(self.data.vz))
        idx = self.query(x, y, z)
        return (self.data.nH[idx], self.data.T[idx],
                self.data.vx[idx], self.data.vy[idx], self.data.vz[idx])


# ---------------------------------------------------------------------------
# HDF5 output writer (same schema as convert_ramses_to_generic.py)
# ---------------------------------------------------------------------------
def _auto_bitpix(arr: np.ndarray) -> int:
    """Choose float32 or float64 based on dynamic range."""
    a = np.abs(arr)
    nz = a[a > 0]
    if nz.size == 0:
        return -32
    ratio = nz.max() / nz.min()
    return -32 if ratio < 1e6 else -64


def write_generic_hdf5(
    filename: Path,
    x: np.ndarray, y: np.ndarray, z: np.ndarray,
    level: np.ndarray,
    nH: np.ndarray, T: np.ndarray,
    vx: np.ndarray, vy: np.ndarray, vz: np.ndarray,
    boxlen: float,
    unit_l_cgs: float,
    xHI: Optional[np.ndarray] = None,
    n_e: Optional[np.ndarray] = None,
    emissivity: Optional[np.ndarray] = None,
    ndust: Optional[np.ndarray] = None,
    metallicity: Optional[np.ndarray] = None,
) -> None:
    """Write generic AMR HDF5 file matching LaRT read_generic_amr.f90 schema."""
    import h5py

    cols = [
        ("x", x), ("y", y), ("z", z),
        ("level", level),
        ("gasDen", nH), ("T", T),
        ("vx", vx), ("vy", vy), ("vz", vz),
    ]
    for name, arr in [("xHI", xHI), ("n_e", n_e), ("emissivity", emissivity),
                      ("ndust", ndust), ("metallicity", metallicity)]:
        if arr is not None:
            cols.append((name, arr))

    p = Path(filename)
    if p.exists():
        p.unlink()
    with h5py.File(p, "w", libver="latest", track_order=True) as f:
        g = f.create_group("AMRGRID", track_order=True)
        for name, arr in cols:
            if name == "level":
                darr = arr.astype(np.int32)
            else:
                bp = _auto_bitpix(arr)
                darr = arr.astype(np.float32 if bp == -32 else np.float64)
            kw = {}
            if darr.size > 4096:
                kw["chunks"] = (min(darr.size, 4096),)
                kw["compression"] = "gzip"
                kw["compression_opts"] = 4
            g.create_dataset(name, data=darr, **kw)

        g.attrs["EXTNAME"]  = "AMRGRID"
        g.attrs["CREATOR"]  = "convert_illustris_to_generic.py"
        g.attrs["DATE"]     = datetime.now().strftime("%Y-%m-%dT%H:%M")
        g.attrs["BOXLEN"]   = float(boxlen)
        g.attrs["UNITLCGS"] = float(unit_l_cgs)
        g.attrs["ORIGINX"]  = -0.5 * float(boxlen)
        g.attrs["ORIGINY"]  = -0.5 * float(boxlen)
        g.attrs["ORIGINZ"]  = -0.5 * float(boxlen)
        g.attrs["NLEAF"]    = np.int32(len(x))
        g.attrs["NAXIS2"]   = np.int32(len(x))


def write_generic_fits(
    filename: Path,
    x: np.ndarray, y: np.ndarray, z: np.ndarray,
    level: np.ndarray,
    nH: np.ndarray, T: np.ndarray,
    vx: np.ndarray, vy: np.ndarray, vz: np.ndarray,
    boxlen: float,
    unit_l_cgs: float,
    xHI: Optional[np.ndarray] = None,
    n_e: Optional[np.ndarray] = None,
    emissivity: Optional[np.ndarray] = None,
    ndust: Optional[np.ndarray] = None,
    metallicity: Optional[np.ndarray] = None,
) -> None:
    """Write generic AMR FITS binary table."""
    from astropy.io import fits as pyfits

    def _bp(arr):
        bp = _auto_bitpix(arr)
        return "E" if bp == -32 else "D"

    cols = [
        pyfits.Column(name="x",      format=_bp(x),  array=x),
        pyfits.Column(name="y",      format=_bp(y),  array=y),
        pyfits.Column(name="z",      format=_bp(z),  array=z),
        pyfits.Column(name="level",  format="J",     array=level.astype(np.int32)),
        pyfits.Column(name="gasDen", format=_bp(nH), array=nH),
        pyfits.Column(name="T",      format=_bp(T),  array=T),
        pyfits.Column(name="vx",     format=_bp(vx), array=vx),
        pyfits.Column(name="vy",     format=_bp(vy), array=vy),
        pyfits.Column(name="vz",     format=_bp(vz), array=vz),
    ]
    for name, arr in [("xHI", xHI), ("n_e", n_e), ("emissivity", emissivity),
                      ("ndust", ndust), ("metallicity", metallicity)]:
        if arr is not None:
            cols.append(pyfits.Column(name=name, format=_bp(arr), array=arr))

    hdu = pyfits.BinTableHDU.from_columns(cols, name="AMRGRID")
    hdu.header["CREATOR"] = "convert_illustris_to_generic.py"
    hdu.header["DATE"]    = datetime.now().strftime("%Y-%m-%dT%H:%M")
    hdu.header["BOXLEN"]  = float(boxlen)
    hdu.header["UNITLCGS"] = float(unit_l_cgs)
    hdu.header["ORIGINX"] = -0.5 * float(boxlen)
    hdu.header["ORIGINY"] = -0.5 * float(boxlen)
    hdu.header["ORIGINZ"] = -0.5 * float(boxlen)
    hdu.header["NLEAF"]   = np.int32(len(x))

    p = Path(filename)
    pyfits.HDUList([pyfits.PrimaryHDU(), hdu]).writeto(p, overwrite=True)


# ---------------------------------------------------------------------------
# Cartesian grid output writers (uniform 3D arrays)
# ---------------------------------------------------------------------------
def _cart_dtype(arr: np.ndarray):
    """Choose float32 or float64 based on dynamic range."""
    return np.float32 if _auto_bitpix(arr) == -32 else np.float64


def write_cartesian_hdf5(
    filename: Path,
    nH_3d: np.ndarray, T_3d: np.ndarray,
    vx_3d: np.ndarray, vy_3d: np.ndarray, vz_3d: np.ndarray,
    boxlen: float,
    unit_l_cgs: float,
    xHI_3d: Optional[np.ndarray] = None,
    n_e_3d: Optional[np.ndarray] = None,
    emissivity_3d: Optional[np.ndarray] = None,
    ndust_3d: Optional[np.ndarray] = None,
    metallicity_3d: Optional[np.ndarray] = None,
) -> None:
    """Write uniform Cartesian grid HDF5 matching LaRT iofile_mod convention.

    Layout (one HDF5 group per quantity, each containing a 'data' dataset):
      /gasDen/data [nz, ny, nx]   (Fortran column-major → stored as [nz,ny,nx] in C order)
      /T/data      [nz, ny, nx]
      /vx/data     [nz, ny, nx]
      /vy/data     [nz, ny, nx]
      /vz/data     [nz, ny, nx]
      ... optional groups ...

    Header attributes are on the first group (/gasDen/).
    """
    import h5py

    nx, ny, nz = nH_3d.shape  # Python: (nx, ny, nz) but stored (nz, ny, nx) for Fortran
    p = Path(filename)
    if p.exists():
        p.unlink()

    mandatory = [
        ("gasDen", nH_3d), ("T", T_3d),
        ("vx", vx_3d), ("vy", vy_3d), ("vz", vz_3d),
    ]
    optional = []
    for name, arr in [("xHI", xHI_3d), ("n_e", n_e_3d),
                      ("emissivity", emissivity_3d), ("ndust", ndust_3d),
                      ("metallicity", metallicity_3d)]:
        if arr is not None:
            optional.append((name, arr))

    with h5py.File(p, "w", libver="latest", track_order=True) as f:
        first_group = True
        for name, arr in mandatory + optional:
            darr = arr.astype(_cart_dtype(arr))
            # Transpose to Fortran order: Python (nx,ny,nz) → stored (nz,ny,nx)
            darr_f = np.asfortranarray(darr).T
            g = f.create_group(name)
            kw = {}
            if darr_f.size > 4096:
                chunk_z = min(darr_f.shape[0], 64)
                chunk_y = min(darr_f.shape[1], 64)
                chunk_x = min(darr_f.shape[2], 64)
                kw["chunks"] = (chunk_z, chunk_y, chunk_x)
                kw["compression"] = "gzip"
                kw["compression_opts"] = 4
            g.create_dataset("data", data=darr_f, **kw)
            # EXTNAME = group name (LaRT reads this)
            g.attrs["EXTNAME"] = np.bytes_(name)

            # Header attributes on the first group only
            if first_group:
                g.attrs["CREATOR"]  = np.bytes_("convert_illustris_to_generic.py")
                g.attrs["DATE"]     = np.bytes_(datetime.now().strftime("%Y-%m-%dT%H:%M"))
                g.attrs["NX"]       = np.int32(nx)
                g.attrs["NY"]       = np.int32(ny)
                g.attrs["NZ"]       = np.int32(nz)
                g.attrs["BOXLEN"]   = np.float64(boxlen)
                g.attrs["UNITLCGS"] = np.float64(unit_l_cgs)
                g.attrs["ORIGINX"]  = np.float64(-0.5 * boxlen)
                g.attrs["ORIGINY"]  = np.float64(-0.5 * boxlen)
                g.attrs["ORIGINZ"]  = np.float64(-0.5 * boxlen)
                first_group = False


def write_cartesian_fits(
    filename: Path,
    nH_3d: np.ndarray, T_3d: np.ndarray,
    vx_3d: np.ndarray, vy_3d: np.ndarray, vz_3d: np.ndarray,
    boxlen: float,
    unit_l_cgs: float,
    xHI_3d: Optional[np.ndarray] = None,
    n_e_3d: Optional[np.ndarray] = None,
    emissivity_3d: Optional[np.ndarray] = None,
    ndust_3d: Optional[np.ndarray] = None,
    metallicity_3d: Optional[np.ndarray] = None,
) -> None:
    """Write uniform Cartesian grid as multi-extension FITS (one ImageHDU per quantity)."""
    from astropy.io import fits as pyfits

    nx, ny, nz = nH_3d.shape

    primary = pyfits.PrimaryHDU()
    primary.header["CREATOR"] = "convert_illustris_to_generic.py"
    primary.header["DATE"]    = datetime.now().strftime("%Y-%m-%dT%H:%M")
    primary.header["NX"]      = nx
    primary.header["NY"]      = ny
    primary.header["NZ"]      = nz
    primary.header["BOXLEN"]  = float(boxlen)
    primary.header["UNITLCGS"] = float(unit_l_cgs)
    primary.header["ORIGINX"] = -0.5 * float(boxlen)
    primary.header["ORIGINY"] = -0.5 * float(boxlen)
    primary.header["ORIGINZ"] = -0.5 * float(boxlen)

    hdus = [primary]
    all_arrays = [
        ("gasDen", nH_3d), ("T", T_3d),
        ("vx", vx_3d), ("vy", vy_3d), ("vz", vz_3d),
    ]
    for name, arr in [("xHI", xHI_3d), ("n_e", n_e_3d),
                      ("emissivity", emissivity_3d), ("ndust", ndust_3d),
                      ("metallicity", metallicity_3d)]:
        if arr is not None:
            all_arrays.append((name, arr))

    for name, arr in all_arrays:
        darr = arr.astype(_cart_dtype(arr))
        hdu = pyfits.ImageHDU(data=darr, name=name)
        hdu.header["BOXLEN"]  = float(boxlen)
        hdu.header["UNITLCGS"] = float(unit_l_cgs)
        hdus.append(hdu)

    p = Path(filename)
    pyfits.HDUList(hdus).writeto(p, overwrite=True)


# ---------------------------------------------------------------------------
# Converter pipeline
# ---------------------------------------------------------------------------
def apply_sfr_treatment(data: IllustrisData, mode: str, cap_temp: float) -> IllustrisData:
    """Handle star-forming cells whose InternalEnergy is on an effective EOS."""
    if data.sfr is None or mode == "as-is":
        return data

    sfr_mask = data.sfr > 0

    if mode == "cap-temperature":
        data.T = np.where(sfr_mask, np.minimum(data.T, cap_temp), data.T)
        n_capped = np.count_nonzero(sfr_mask)
        if n_capped > 0:
            print(f"  SFR cap: {n_capped} cells capped at T <= {cap_temp:.0f} K")
        return data

    elif mode == "exclude":
        keep = ~sfr_mask
        n_excl = np.count_nonzero(sfr_mask)
        print(f"  SFR exclude: removed {n_excl} cells with SFR > 0")
        return IllustrisData(
            x=data.x[keep], y=data.y[keep], z=data.z[keep],
            nH=data.nH[keep], T=data.T[keep],
            vx=data.vx[keep], vy=data.vy[keep], vz=data.vz[keep],
            metallicity=data.metallicity[keep] if data.metallicity is not None else None,
            xHI=data.xHI[keep] if data.xHI is not None else None,
            electron_abundance=data.electron_abundance[keep] if data.electron_abundance is not None else None,
            sfr=data.sfr[keep],
            volume_kpc3=data.volume_kpc3[keep] if data.volume_kpc3 is not None else None,
            hubble_param=data.hubble_param,
            redshift=data.redshift,
            scale_factor=data.scale_factor,
            boxsize_ckpc_h=data.boxsize_ckpc_h,
        )

    return data


def convert(args):
    """Main conversion pipeline."""

    # ---- 1. Load data ----
    if args.input_file is not None:
        filepath = args.input_file
    elif args.api_key:
        filepath = download_tng_cutout(
            args.simulation, args.snap, args.subhalo_id,
            args.api_key, cache_dir=str(Path(args.output).parent),
        )
    else:
        print("ERROR: either provide an input file or --api-key with "
              "--simulation/--snap/--subhalo-id", file=sys.stderr)
        sys.exit(1)

    print(f"Loading {filepath} ...")
    center = np.array(args.center, dtype=float) if args.center else None
    data = load_illustris_hdf5(str(filepath), center_kpc=center)

    ncells = len(data.x)
    print(f"  {ncells:,} gas cells loaded")
    print(f"  z = {data.redshift:.4f},  a = {data.scale_factor:.4f},  h = {data.hubble_param:.4f}")

    # ---- 2. ISM treatment ----
    data = apply_sfr_treatment(data, args.sfr_treatment, args.sfr_cap_temp)
    ncells = len(data.x)

    # ---- 3. Determine box size ----
    if args.boxsize is not None:
        boxlen = args.boxsize
    else:
        # Auto: bounding box of data + 10% margin
        extent = max(data.x.max() - data.x.min(),
                     data.y.max() - data.y.min(),
                     data.z.max() - data.z.min())
        boxlen = extent * 1.1
    print(f"  Box size: {boxlen:.2f} kpc")

    # Recenter to box center if not already centered
    if center is None:
        cx = 0.5 * (data.x.max() + data.x.min())
        cy = 0.5 * (data.y.max() + data.y.min())
        cz = 0.5 * (data.z.max() + data.z.min())
        data.x -= cx
        data.y -= cy
        data.z -= cz
        print(f"  Auto-centered at ({cx:.2f}, {cy:.2f}, {cz:.2f}) kpc")

    # Clip to box (cells outside [-boxlen/2, boxlen/2])
    half = boxlen / 2.0
    inside = ((np.abs(data.x) < half) & (np.abs(data.y) < half) &
              (np.abs(data.z) < half))
    n_outside = ncells - np.count_nonzero(inside)
    if n_outside > 0:
        print(f"  Clipped {n_outside} cells outside box")
        data = IllustrisData(
            x=data.x[inside], y=data.y[inside], z=data.z[inside],
            nH=data.nH[inside], T=data.T[inside],
            vx=data.vx[inside], vy=data.vy[inside], vz=data.vz[inside],
            metallicity=data.metallicity[inside] if data.metallicity is not None else None,
            xHI=data.xHI[inside] if data.xHI is not None else None,
            electron_abundance=data.electron_abundance[inside] if data.electron_abundance is not None else None,
            sfr=data.sfr[inside] if data.sfr is not None else None,
            volume_kpc3=data.volume_kpc3[inside] if data.volume_kpc3 is not None else None,
            hubble_param=data.hubble_param,
            redshift=data.redshift,
            scale_factor=data.scale_factor,
            boxsize_ckpc_h=data.boxsize_ckpc_h,
        )
        ncells = len(data.x)

    # ---- 4. Print data summary ----
    print(f"\n  Data summary ({ncells:,} cells):")
    print(f"    nH [cm^-3]:  min={data.nH.min():.3e}  max={data.nH.max():.3e}  "
          f"median={np.median(data.nH):.3e}")
    print(f"    T  [K]:      min={data.T.min():.3e}  max={data.T.max():.3e}  "
          f"median={np.median(data.T):.3e}")
    vmag = np.sqrt(data.vx**2 + data.vy**2 + data.vz**2)
    print(f"    |v| [km/s]:  min={vmag.min():.1f}  max={vmag.max():.1f}  "
          f"median={np.median(vmag):.1f}")
    if data.metallicity is not None:
        print(f"    Z:           min={data.metallicity.min():.3e}  "
              f"max={data.metallicity.max():.3e}")
    if data.xHI is not None:
        print(f"    xHI (TNG):   min={data.xHI.min():.3e}  max={data.xHI.max():.3e}")

    # ---- 5. Build KD-tree interpolator ----
    print("\nBuilding KD-tree ...")
    positions = np.column_stack([data.x, data.y, data.z])
    interp = VoronoiInterpolator(
        positions, data,
        method=args.resample_method,
        kernel_factor=args.kernel_size_factor,
        kernel_neighbors=args.kernel_neighbors,
        kernel_size=args.kernel_size)
    print(f"  KD-tree built with {ncells:,} points")
    if args.resample_method == "gaussian":
        if args.kernel_size and args.kernel_size > 0:
            print(f"  Resampling: Gaussian kernel, fixed h = {args.kernel_size} "
                  f"{args.output_unit} (K={interp.kernel_neighbors})")
        else:
            print(f"  Resampling: Gaussian kernel, adaptive h = "
                  f"{args.kernel_size_factor} x r_eff (K={interp.kernel_neighbors})")
    else:
        print("  Resampling: nearest-neighbor (Voronoi-exact)")

    # ---- 6–9. Build grid, assign properties, compute physics, write output ----
    outpath = Path(args.output)
    unit_l_cgs = KPC_CM  # default: kpc
    if args.output_unit.lower() == "pc":
        unit_l_cgs = PC_CM
    elif args.output_unit.lower() == "cm":
        unit_l_cgs = 1.0

    # The diagnostic plot (if --plot) is produced inside the conversion
    # routine, where the octree / Cartesian grid arrays are available.
    if args.grid_type == "cartesian":
        _convert_cartesian(args, data, interp, boxlen, outpath, unit_l_cgs)
    else:
        _convert_amr(args, data, interp, boxlen, outpath, unit_l_cgs)


# ---------------------------------------------------------------------------
# Shared physics computation
# ---------------------------------------------------------------------------
def _compute_physics(args, data, interp, nH_arr, T_arr, cx_arr, cy_arr, cz_arr):
    """Compute optional physics columns (xHI, n_e, emissivity, ndust, metallicity).

    Returns (xHI_arr, ne_arr, emiss_arr, ndust_arr, Z_arr), any of which may be None.
    """
    xHI_arr = None
    ne_arr  = None
    emiss_arr = None
    ndust_arr = None
    Z_arr   = None

    if not (args.compute_physics or args.ionization != "none"):
        return xHI_arr, ne_arr, emiss_arr, ndust_arr, Z_arr

    # Native data-file fields (metallicity, xHI, ElectronAbundance) are
    # resampled with the SAME method as nH/T/v (nearest or Gaussian kernel).
    # Metallicity
    if data.metallicity is not None:
        Z_arr = interp.sample(data.metallicity, cx_arr, cy_arr, cz_arr)

    # Ionization
    if args.ionization == "from_illustris" and data.xHI is not None:
        xHI_arr = interp.sample(data.xHI, cx_arr, cy_arr, cz_arr)
        print("  xHI: from Illustris NeutralHydrogenAbundance")
    elif args.ionization == "cie":
        xHI_arr = cie_neutral_fraction(T_arr)
        print("  xHI: CIE neutral fraction (simple formula)")
    elif args.ionization == "full_neutral":
        xHI_arr = np.ones_like(T_arr)
        print("  xHI: full neutral (xHI = 1)")

    # Electron density
    if xHI_arr is not None:
        if args.ionization == "from_illustris" and data.electron_abundance is not None:
            ne_arr = nH_arr * interp.sample(data.electron_abundance,
                                            cx_arr, cy_arr, cz_arr)
            print("  n_e: from Illustris ElectronAbundance")
        else:
            ne_arr = electron_density(nH_arr, xHI_arr)
            print("  n_e: from xHI")

    # Emissivity and dust
    if args.compute_physics and xHI_arr is not None and ne_arr is not None:
        eline = getattr(args, 'emissivity_line', 'lya')
        if eline == 'civ':
            if Z_arr is None:
                print("  WARNING: --emissivity-line civ requires metallicity; skipping")
            else:
                emiss_arr = cie_civ_emissivity(nH_arr, T_arr, ne_arr, Z_arr)
                print(f"  emissivity: CIE C IV 1548+1550, "
                      f"max = {emiss_arr.max():.3e} cm^-3 s^-1")
        elif eline == 'ovi':
            if Z_arr is None:
                print("  WARNING: --emissivity-line ovi requires metallicity; skipping")
            else:
                emiss_arr = cie_ovi_emissivity(nH_arr, T_arr, ne_arr, Z_arr)
                print(f"  emissivity: CIE O VI 1032+1038, "
                      f"max = {emiss_arr.max():.3e} cm^-3 s^-1")
        else:
            emiss_arr = caseB_lya_emissivity(nH_arr, T_arr, xHI_arr, ne_arr)
            print(f"  emissivity: Case B Lya, "
                  f"max = {emiss_arr.max():.3e} cm^-3 s^-1")

        if Z_arr is not None:
            ndust_arr = laursen09_ndust(nH_arr, xHI_arr, Z_arr)
            print(f"  ndust: Laursen+09, max = {ndust_arr.max():.3e} cm^-3")

    return xHI_arr, ne_arr, emiss_arr, ndust_arr, Z_arr


# ---------------------------------------------------------------------------
# AMR grid conversion (existing path, refactored)
# ---------------------------------------------------------------------------
def _convert_amr(args, data, interp, boxlen, outpath, unit_l_cgs):
    """Build adaptive octree and write generic AMR format."""

    print(f"\nBuilding octree (level_min={args.level_min}, level_max={args.level_max}) ...")

    # Import AMRGrid
    amr_grid_dir = Path(__file__).parent
    if str(amr_grid_dir) not in sys.path:
        sys.path.insert(0, str(amr_grid_dir))
    from AMR_grid import AMRGrid

    grid = AMRGrid(boxlen)

    # Resolution-matching refinement criterion
    if args.match_resolution and data.volume_kpc3 is not None:
        r_eff = (3.0 * data.volume_kpc3 / (4.0 * np.pi))**(1.0 / 3.0)
        factor = args.resolution_factor

        def resolution_criterion(cell):
            idx = interp.query(cell.cx, cell.cy, cell.cz)
            return cell.h > factor * r_eff[idx]

        grid.refine_by_physics(
            dens_fn=interp.density_fn,
            vel_fn=interp.velocity_fn,
            dens_threshold=args.dens_threshold,
            vel_threshold=args.vel_threshold,
            level_max=args.level_max,
            level_min=args.level_min,
            nprobe=2,
        )
        grid.refine(resolution_criterion, args.level_max)
    else:
        grid.refine_by_physics(
            dens_fn=interp.density_fn,
            vel_fn=interp.velocity_fn,
            dens_threshold=args.dens_threshold,
            vel_threshold=args.vel_threshold,
            level_max=args.level_max,
            level_min=args.level_min,
            nprobe=2,
        )

    leaflist = grid.leaves()
    nleaf = len(leaflist)
    print(f"  Octree: {nleaf:,} leaf cells")
    print(f"  {grid.info()}")

    # Assign leaf-cell properties
    print("\nAssigning leaf-cell properties via nearest-neighbor ...")
    grid.set_properties(interp.properties_fn)

    # Extract leaf arrays
    cx_arr = np.array([lf.cx for lf in leaflist])
    cy_arr = np.array([lf.cy for lf in leaflist])
    cz_arr = np.array([lf.cz for lf in leaflist])
    level_arr = np.array([lf.level for lf in leaflist], dtype=np.int32)
    nH_arr = np.array([lf.dens for lf in leaflist])
    T_arr  = np.array([lf.T for lf in leaflist])
    vx_arr = np.array([lf.vx for lf in leaflist])
    vy_arr = np.array([lf.vy for lf in leaflist])
    vz_arr = np.array([lf.vz for lf in leaflist])

    # Physics
    xHI_arr, ne_arr, emiss_arr, ndust_arr, Z_arr = \
        _compute_physics(args, data, interp, nH_arr, T_arr, cx_arr, cy_arr, cz_arr)

    # Write output
    print(f"\nWriting {outpath} ({nleaf:,} leaf cells) ...")
    name_lower = outpath.name.lower()
    if name_lower.endswith(".h5") or name_lower.endswith(".hdf5"):
        write_generic_hdf5(outpath, cx_arr, cy_arr, cz_arr, level_arr,
                           nH_arr, T_arr, vx_arr, vy_arr, vz_arr,
                           boxlen, unit_l_cgs,
                           xHI=xHI_arr, n_e=ne_arr, emissivity=emiss_arr,
                           ndust=ndust_arr, metallicity=Z_arr)
    elif name_lower.endswith(".fits") or name_lower.endswith(".fits.gz"):
        write_generic_fits(outpath, cx_arr, cy_arr, cz_arr, level_arr,
                           nH_arr, T_arr, vx_arr, vy_arr, vz_arr,
                           boxlen, unit_l_cgs,
                           xHI=xHI_arr, n_e=ne_arr, emissivity=emiss_arr,
                           ndust=ndust_arr, metallicity=Z_arr)
    else:
        header = f"{nleaf}  {boxlen:.10e}\n"
        header += "# x  y  z  level  gasDen  T  vx  vy  vz"
        extra_cols = []
        for nm, ar in [("metallicity", Z_arr), ("xHI", xHI_arr), ("n_e", ne_arr),
                       ("emissivity", emiss_arr), ("ndust", ndust_arr)]:
            if ar is not None:
                header += f"  {nm}"
                extra_cols.append(ar)
        data_cols = np.column_stack(
            [cx_arr, cy_arr, cz_arr, level_arr, nH_arr, T_arr,
             vx_arr, vy_arr, vz_arr] + extra_cols
        )
        np.savetxt(outpath, data_cols, header=header, comments="",
                   fmt="%.8e", delimiter="  ")

    print(f"  Done: {outpath}")

    # Suggested LaRT input
    dist_str = f"par%distance_unit = '{args.output_unit}'"
    if args.output_unit.lower() == "cm":
        dist_str = "par%distance2cm = 1.0"
    print(f"""
{'─' * 60}
Suggested LaRT input snippet:
{'─' * 60}
 &parameters
  par%use_amr_grid    = .true.
  par%amr_type        = 'generic'
  par%amr_file        = '{outpath.name}'
  {dist_str}
  par%no_photons      = 1e+06
  par%spectral_type   = 'voigt'
  par%source_geometry = 'diffuse_emissivity'
  par%nxfreq  = 121
  par%nxim    = 100
  par%nyim    = 100
 /
{'─' * 60}""")

    if args.plot:
        make_diagnostic_plot(data, interp, boxlen, outpath,
                             grid_kind="amr",
                             leaf_cx=cx_arr, leaf_cy=cy_arr, leaf_cz=cz_arr,
                             leaf_level=level_arr, leaf_nH=nH_arr, leaf_T=T_arr,
                             scatter=args.plot_scatter)


# ---------------------------------------------------------------------------
# Cartesian grid conversion (new path)
# ---------------------------------------------------------------------------
def _convert_cartesian(args, data, interp, boxlen, outpath, unit_l_cgs):
    """Build uniform N^3 grid and write Cartesian grid format."""

    N = args.ngrid
    print(f"\nBuilding uniform Cartesian grid ({N} x {N} x {N}) ...")

    half = boxlen / 2.0
    dx = boxlen / N

    # Cell centers: centered convention [-half+dx/2, half-dx/2]
    edges = np.linspace(-half, half, N + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    # 3D meshgrid
    cx_1d, cy_1d, cz_1d = centers, centers, centers
    CX, CY, CZ = np.meshgrid(cx_1d, cy_1d, cz_1d, indexing='ij')
    # Flatten for KD-tree query
    cx_flat = CX.ravel()
    cy_flat = CY.ravel()
    cz_flat = CZ.ravel()

    ntotal = N * N * N
    print(f"  {ntotal:,} cells, dx = {dx:.4f} kpc")

    # Assign properties (nearest-neighbor or Gaussian kernel, per --resample-method)
    print("  Querying KD-tree ...")
    nH_flat, T_flat, vx_flat, vy_flat, vz_flat = \
        interp.properties_fn(cx_flat, cy_flat, cz_flat)

    # Reshape to 3D (nx, ny, nz)
    nH_3d = nH_flat.reshape((N, N, N))
    T_3d  = T_flat.reshape((N, N, N))
    vx_3d = vx_flat.reshape((N, N, N))
    vy_3d = vy_flat.reshape((N, N, N))
    vz_3d = vz_flat.reshape((N, N, N))

    print(f"  nH:  min={nH_3d.min():.3e}  max={nH_3d.max():.3e}")
    print(f"  T:   min={T_3d.min():.3e}  max={T_3d.max():.3e}")

    # Physics
    xHI_arr, ne_arr, emiss_arr, ndust_arr, Z_arr = \
        _compute_physics(args, data, interp, nH_flat, T_flat, cx_flat, cy_flat, cz_flat)

    # Reshape optional arrays to 3D
    xHI_3d   = xHI_arr.reshape((N, N, N)) if xHI_arr is not None else None
    ne_3d    = ne_arr.reshape((N, N, N)) if ne_arr is not None else None
    emiss_3d = emiss_arr.reshape((N, N, N)) if emiss_arr is not None else None
    ndust_3d = ndust_arr.reshape((N, N, N)) if ndust_arr is not None else None
    Z_3d     = Z_arr.reshape((N, N, N)) if Z_arr is not None else None

    # Write output
    print(f"\nWriting {outpath} ({ntotal:,} cells, {N}^3 Cartesian grid) ...")
    name_lower = outpath.name.lower()
    if name_lower.endswith(".h5") or name_lower.endswith(".hdf5"):
        write_cartesian_hdf5(outpath, nH_3d, T_3d, vx_3d, vy_3d, vz_3d,
                             boxlen, unit_l_cgs,
                             xHI_3d=xHI_3d, n_e_3d=ne_3d,
                             emissivity_3d=emiss_3d, ndust_3d=ndust_3d,
                             metallicity_3d=Z_3d)
    elif name_lower.endswith(".fits") or name_lower.endswith(".fits.gz"):
        write_cartesian_fits(outpath, nH_3d, T_3d, vx_3d, vy_3d, vz_3d,
                             boxlen, unit_l_cgs,
                             xHI_3d=xHI_3d, n_e_3d=ne_3d,
                             emissivity_3d=emiss_3d, ndust_3d=ndust_3d,
                             metallicity_3d=Z_3d)
    else:
        print("  ERROR: Cartesian grid only supports .h5 or .fits/.fits.gz output",
              file=sys.stderr)
        sys.exit(1)

    print(f"  Done: {outpath}")

    # Suggested LaRT input
    dist_str = f"par%distance_unit = '{args.output_unit}'"
    if args.output_unit.lower() == "cm":
        dist_str = "par%distance2cm = 1.0"
    src_str = "par%source_geometry = 'diffuse_emissivity'"
    if emiss_3d is None:
        src_str = "par%source_geometry = 'point'"
    print(f"""
{'─' * 60}
Suggested LaRT input snippet:
{'─' * 60}
 &parameters
  par%cart_file       = '{outpath.name}'
  {dist_str}
  {src_str}
  par%no_photons      = 1e+06
  par%spectral_type   = 'voigt'
  par%nxfreq  = 121
  par%nxim    = 100
  par%nyim    = 100
 /
{'─' * 60}""")

    if args.plot:
        make_diagnostic_plot(data, interp, boxlen, outpath,
                             grid_kind="cartesian",
                             cart_nH_3d=nH_3d, cart_T_3d=T_3d, cart_centers=centers,
                             leaf_cx=cx_flat, leaf_cy=cy_flat, leaf_cz=cz_flat,
                             leaf_nH=nH_flat, leaf_T=T_flat,
                             scatter=args.plot_scatter)


def _draw_octree_slice(ax, plt, cx, cy, cz, level, values, boxlen,
                       cmap, vmin, vmax, label):
    """Fill the z=0 plane with true octree leaf squares (area-filled slice).

    A leaf cell of side ``boxlen / 2^level`` is intersected by the z=0 plane
    when ``|cz| <= side/2``; each such cell is drawn as a filled square.
    """
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle
    size = boxlen / (2.0 ** level)
    halfc = size / 2.0
    sel = np.where(np.abs(cz) <= halfc)[0]
    rects = [Rectangle((cx[i] - halfc[i], cy[i] - halfc[i]), size[i], size[i])
             for i in sel]
    pc = PatchCollection(rects, cmap=cmap)
    pc.set_array(values[sel])
    pc.set_clim(vmin, vmax)
    ax.add_collection(pc)
    plt.colorbar(pc, ax=ax, label=label)


def make_diagnostic_plot(data, interp, boxlen, outpath, *,
                         grid_kind="amr",
                         leaf_cx=None, leaf_cy=None, leaf_cz=None,
                         leaf_level=None, leaf_nH=None, leaf_T=None,
                         cart_nH_3d=None, cart_T_3d=None, cart_centers=None,
                         scatter=False):
    """Side-by-side original Voronoi vs resampled grid, for nH and T.

    Default (``scatter=False``): true area-filled z=0 slices. The left column
    is the exact Voronoi slice (nearest-neighbor of the 3D mesh sampled on a
    pixel grid at z=0); the right column is the resampled grid -- octree leaf
    squares (``grid_kind='amr'``) or a Cartesian image (``grid_kind='cartesian'``).

    With ``scatter=True`` the old behavior is reproduced: a scatter of cell
    centers lying in a thin ``|z| < boxlen/50`` slab.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib not available; skipping plot")
        return

    half = boxlen / 2.0
    extent = [-half, half, -half, half]
    nH_label = r"$\log_{10}(n_H$ [cm$^{-3}$])"
    T_label  = r"$\log_{10}(T$ [K])"
    grid_name = "Cartesian" if grid_kind == "cartesian" else "Octree"

    def _cosmetics(ax, title):
        ax.set_xlim(-half, half)
        ax.set_ylim(-half, half)
        ax.set_title(title)
        ax.set_xlabel("x [kpc]")
        ax.set_ylabel("y [kpc]")
        ax.set_aspect("equal")

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    if scatter:
        dz = boxlen / 50.0
        mask_v = np.abs(data.z) < dz
        mask_o = np.abs(leaf_cz) < dz
        for (ax, mx, my, mc, cmap, vmin, vmax, label, title) in [
            (axes[0, 0], data.x, data.y, np.log10(np.maximum(data.nH, 1e-10)),
             "viridis", -6, 2, nH_label, "Voronoi: nH (z=0 slab)"),
            (axes[0, 1], leaf_cx, leaf_cy, np.log10(np.maximum(leaf_nH, 1e-10)),
             "viridis", -6, 2, nH_label, f"{grid_name}: nH (z=0 slab)"),
            (axes[1, 0], data.x, data.y, np.log10(np.maximum(data.T, 10.0)),
             "hot", 3, 7, T_label, "Voronoi: T (z=0 slab)"),
            (axes[1, 1], leaf_cx, leaf_cy, np.log10(np.maximum(leaf_T, 10.0)),
             "hot", 3, 7, T_label, f"{grid_name}: T (z=0 slab)"),
        ]:
            m = mask_v if mx is data.x else mask_o
            if m.any():
                sc = ax.scatter(mx[m], my[m], c=mc[m], s=1,
                                cmap=cmap, vmin=vmin, vmax=vmax)
                plt.colorbar(sc, ax=ax, label=label)
            _cosmetics(ax, title)
    else:
        # LEFT: exact Voronoi slice via nearest-neighbor on a pixel grid at z=0
        npix = 512
        px = np.linspace(-half, half, npix)
        PX, PY = np.meshgrid(px, px, indexing="ij")
        idx = interp.query(PX.ravel(), PY.ravel(), np.zeros(PX.size))
        vor_nH = np.log10(np.maximum(data.nH[idx].reshape(npix, npix), 1e-10))
        vor_T  = np.log10(np.maximum(data.T[idx].reshape(npix, npix), 10.0))

        im = axes[0, 0].imshow(vor_nH.T, origin="lower", extent=extent,
                               cmap="viridis", vmin=-6, vmax=2)
        plt.colorbar(im, ax=axes[0, 0], label=nH_label)
        _cosmetics(axes[0, 0], "Voronoi: nH (z=0 slice)")

        im = axes[1, 0].imshow(vor_T.T, origin="lower", extent=extent,
                               cmap="hot", vmin=3, vmax=7)
        plt.colorbar(im, ax=axes[1, 0], label=T_label)
        _cosmetics(axes[1, 0], "Voronoi: T (z=0 slice)")

        # RIGHT: resampled grid slice
        if grid_kind == "cartesian":
            k = int(np.argmin(np.abs(cart_centers)))
            grid_nH = np.log10(np.maximum(cart_nH_3d[:, :, k], 1e-10))
            grid_T  = np.log10(np.maximum(cart_T_3d[:, :, k], 10.0))
            im = axes[0, 1].imshow(grid_nH.T, origin="lower", extent=extent,
                                   cmap="viridis", vmin=-6, vmax=2)
            plt.colorbar(im, ax=axes[0, 1], label=nH_label)
            im = axes[1, 1].imshow(grid_T.T, origin="lower", extent=extent,
                                   cmap="hot", vmin=3, vmax=7)
            plt.colorbar(im, ax=axes[1, 1], label=T_label)
        else:
            _draw_octree_slice(axes[0, 1], plt, leaf_cx, leaf_cy, leaf_cz,
                               leaf_level, np.log10(np.maximum(leaf_nH, 1e-10)),
                               boxlen, "viridis", -6, 2, nH_label)
            _draw_octree_slice(axes[1, 1], plt, leaf_cx, leaf_cy, leaf_cz,
                               leaf_level, np.log10(np.maximum(leaf_T, 10.0)),
                               boxlen, "hot", 3, 7, T_label)
        _cosmetics(axes[0, 1], f"{grid_name}: nH (z=0 slice)")
        _cosmetics(axes[1, 1], f"{grid_name}: T (z=0 slice)")

    plt.tight_layout()
    figpath = outpath.with_suffix(".png")
    plt.savefig(figpath, dpi=150)
    print(f"  Diagnostic plot: {figpath}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Convert Illustris/TNG data to LaRT generic AMR format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:

  # Local cutout file:
  %(prog)s cutout.hdf5 -o galaxy.h5 --compute-physics

  # API download:
  %(prog)s --api-key YOUR_KEY --simulation TNG50-1 --snap 99 \\
           --subhalo-id 0 -o galaxy.h5 --compute-physics

  # Full options:
  %(prog)s cutout.hdf5 -o galaxy.h5 --output-unit kpc \\
           --level-max 8 --level-min 3 --dens-threshold 0.3 \\
           --match-resolution --compute-physics \\
           --ionization from_illustris --sfr-treatment cap-temperature \\
           --plot
""",
    )

    # Input
    parser.add_argument("input_file", nargs="?", default=None,
                        help="Illustris/TNG HDF5 file (cutout or snapshot chunk)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output file (.h5, .fits.gz, or .dat)")
    parser.add_argument("--output-unit", default="kpc",
                        help="Unit of output coordinates (default: kpc)")

    # Grid type
    parser.add_argument("--grid-type", default="amr", choices=["amr", "cartesian"],
                        help="Output grid type: amr (adaptive octree, default) or "
                             "cartesian (uniform N^3 grid)")
    parser.add_argument("--ngrid", type=int, default=128,
                        help="Grid cells per side for --grid-type cartesian (default: 128)")

    # Octree
    parser.add_argument("--level-max", type=int, default=8,
                        help="Maximum AMR level (default: 8)")
    parser.add_argument("--level-min", type=int, default=3,
                        help="Minimum uniform refinement level (default: 3)")
    parser.add_argument("--dens-threshold", type=float, default=0.3,
                        help="Density gradient threshold (default: 0.3)")
    parser.add_argument("--vel-threshold", type=float, default=0.3,
                        help="Velocity gradient threshold (default: 0.3)")
    parser.add_argument("--match-resolution", action="store_true",
                        help="Also refine to match local Voronoi cell size")
    parser.add_argument("--resolution-factor", type=float, default=1.0,
                        help="Factor for resolution matching (default: 1.0)")

    # Resampling (Voronoi -> grid value assignment)
    parser.add_argument("--resample-method", default="nearest",
                        choices=["nearest", "gaussian"],
                        help="How cell values are assigned to grid points: "
                             "nearest (Voronoi-exact, default) or gaussian "
                             "(volume-weighted kernel smoothing)")
    parser.add_argument("--kernel-size-factor", type=float, default=1.0,
                        help="Gaussian mode: adaptive smoothing length "
                             "h = factor x cell effective radius r_eff "
                             "(default: 1.0)")
    parser.add_argument("--kernel-neighbors", type=int, default=32,
                        help="Gaussian mode: number of nearest cells gathered "
                             "per grid point (default: 32)")
    parser.add_argument("--kernel-size", type=float, default=None,
                        help="Gaussian mode: fixed smoothing length in output "
                             "units, overriding the adaptive --kernel-size-factor")

    # Physics
    parser.add_argument("--compute-physics", action="store_true",
                        help="Compute xHI, n_e, emissivity, ndust")
    parser.add_argument("--ionization", default="from_illustris",
                        choices=["from_illustris", "cie", "full_neutral", "none"],
                        help="Ionization source (default: from_illustris)")
    parser.add_argument("--emissivity-line", default="lya",
                        choices=["lya", "civ", "ovi"],
                        help="Emission line for emissivity column: "
                             "lya (Case B Lyman-alpha), "
                             "civ (CIE C IV 1548+1550), "
                             "ovi (CIE O VI 1032+1038). Default: lya")

    # ISM
    parser.add_argument("--sfr-treatment", default="cap-temperature",
                        choices=["cap-temperature", "exclude", "as-is"],
                        help="Treatment of star-forming cells (default: cap-temperature)")
    parser.add_argument("--sfr-cap-temp", type=float, default=2.0e4,
                        help="Temperature cap for SFR>0 cells (default: 2e4 K)")

    # Geometry
    parser.add_argument("--center", nargs=3, type=float, default=None,
                        metavar=("X", "Y", "Z"),
                        help="Center coordinates in physical kpc (default: auto)")
    parser.add_argument("--boxsize", type=float, default=None,
                        help="Box side length in kpc (default: auto from data)")

    # API
    parser.add_argument("--api-key", default=None,
                        help="TNG API key for downloading cutouts")
    parser.add_argument("--simulation", default="TNG50-1",
                        help="Simulation name (default: TNG50-1)")
    parser.add_argument("--snap", type=int, default=99,
                        help="Snapshot number (default: 99 = z~0)")
    parser.add_argument("--subhalo-id", type=int, default=0,
                        help="Subhalo ID (default: 0)")

    # Misc
    parser.add_argument("--plot", action="store_true",
                        help="Generate diagnostic slice plot (area-filled slices "
                             "by default)")
    parser.add_argument("--plot-scatter", action="store_true",
                        help="With --plot, draw the old scatter-of-cell-centers "
                             "style (thin |z|<boxlen/50 slab) instead of "
                             "area-filled slices")

    args = parser.parse_args()
    convert(args)


if __name__ == "__main__":
    main()
