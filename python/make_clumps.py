#!/usr/bin/env python3
"""make_clumps.py -- Python port of LaRT's clump_mod.f90 / make_clumps.x

Generate a population of spherical clumps (positions, bulk velocities,
clump physical parameters) and write a FITS or HDF5 file with the
same schema as the Fortran writer in clump_mod.f90. The output file is
consumed by LaRT via ``par%clump_input_file = '<file>'``.

Algorithm
---------
Random Sequential Addition (RSA) with a uniform-grid hash for O(N) overlap
checks. Two sampling paths, both supporting an optional biconical mask
along the z-axis:

  - uniform (default): direct r,theta,phi sampling within the spherical
    shell [rmin, rmax] (or the bicone subset of it).
  - radial profile (any of clump_*_profile != 'constant'): inverse-CDF
    radial draw weighted by shape_number(r) * r^2; isotropic phi; cos(theta)
    direct in the cone window; clump radius / opacity / temperature
    derived from shape_radius / shape_density / file table.

The cone-fully-inside check uses
    cos(theta_min) = |cos(theta)| * sqrt(1 - (r_cl/r)^2) - sin(theta) * (r_cl/r),
matching clump_mod.f90:953.

Output schema
-------------
EXTNAME='Clumps' BinTable (FITS) or '/Clumps' group (HDF5):

  columns:   X Y Z                [code units]
             VX VY VZ             [km/s]
             R_CLUMP, RHOKAP, TEMP  (only when non-constant; const_tol = 1e-3)
  header:    N_CLUMPS, SPHERE_R, RMIN, CL_RAD, F_VOL, F_COV,
             TAU0, SIGMA_V, TEMP_CL, RHOKAP, CL_DFREQ, VTHERM, VOIGT_A,
             RMAX, IN_FCOV, IN_FVOL, IN_NCL, IN_NHI, IN_NH, IN_TEMP,
             DISTUNIT, DIST_CM

Usage
-----
  # uniform sphere, 5000 clumps, r_cl=0.02, tau0=1e3, Lya, T=1e4
  python make_clumps.py --rmax 1.0 --clump_radius 0.02 \\
         --clump_N_clumps 5000 --clump_tau0 1e3 --temperature 1e4 \\
         -o uniform_clumps.fits.gz

  # by covering factor instead of N
  python make_clumps.py --rmax 1.0 --clump_radius 1e-3 \\
         --clump_f_cov 5 --clump_tau0 1e3 -o fcov5_clumps.fits.gz

  # biconical outflow, 45 deg half-opening, Hubble flow at 200 km/s
  python make_clumps.py --rmax 1.0 --clump_radius 0.02 --clump_f_cov 2 \\
         --cone_opening 45 --clump_tau0 1e3 \\
         --velocity_type hubble --Vexp 200 \\
         -o bicone_clumps.fits.gz

  # powerlaw radial number density profile
  python make_clumps.py --rmax 1.0 --clump_radius 0.02 --clump_f_vol 0.1 \\
         --clump_number_profile powerlaw --clump_number_alpha 2.0 \\
         --clump_number_r0 0.3 --clump_tau0 1e3 -o profile_clumps.fits.gz
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from typing import Dict, Tuple

import numpy as np

# numba is used to JIT the RSA inner loop (overlap test + linked-list
# insert). Falls back to a pure-Python path if numba is unavailable.
try:
    from numba import njit
    _HAVE_NUMBA = True
except Exception:
    _HAVE_NUMBA = False
    def njit(*a, **k):
        def deco(fn): return fn
        return deco if (len(a) == 0 or not callable(a[0])) else a[0]

# ---------------------------------------------------------------------------
# Physical constants (mirror line_mod.f90 / define.f90)
# ---------------------------------------------------------------------------
SIGMA_0     = 0.026540083434       # pi * e^2 / (m_e c) [cm^2 Hz]
VTHERM1_AMU = 0.12895319011972164  # sqrt(2k_B/m_amu) [km/s @ 1 K, 1 amu]
UM2KM       = 1.0e-9               # 1 micron in km

PI     = np.pi
TWOPI  = 2.0 * np.pi
FOURPI = 4.0 * np.pi

# Line-specific atomic data (singlet f12 used to compute cross0; for doublets
# this is the strong component, matching f12(1) in setup_resonance_line()).
LINE_DATA: Dict[str, Dict[str, float]] = {
    'ly_alpha':  dict(wavelength0=0.1215668237310, damping=6.2649e8, f12=0.4164,  mass_amu=1.00797),
    'CIV_1548':  dict(wavelength0=0.1548187,       damping=2.647e8,  f12=0.190,   mass_amu=12.011),
    'NV_1239':   dict(wavelength0=0.1238821,       damping=3.390e8,  f12=0.156,   mass_amu=14.007),
    'OVI_1032':  dict(wavelength0=0.1031912,       damping=4.137e8,  f12=0.1325,  mass_amu=15.999),
    'NaI_D':     dict(wavelength0=0.5891583253,    damping=6.153e7,  f12=0.641,   mass_amu=22.990),
    'CaII_HK':   dict(wavelength0=0.3934777,       damping=1.446667e8, f12=0.682, mass_amu=40.078),
    'MgII_2796': dict(wavelength0=0.2796352,       damping=2.590e8,  f12=0.6080,  mass_amu=24.305),
    'SiIV_1394': dict(wavelength0=0.1393755,       damping=8.743e8,  f12=0.513,   mass_amu=28.085),
    'AlII_1671': dict(wavelength0=0.1670787,       damping=1.39e9,   f12=1.77,    mass_amu=26.982),
    'CII_1334':  dict(wavelength0=0.1334532,       damping=2.880e8,  f12=0.129,   mass_amu=12.011),
}

# Distance unit -> cm (mirrors define.f90 unit table).
DIST_UNIT_CM: Dict[str, float] = {
    '':       0.0,
    'cm':     1.0,
    'm':      1.0e2,
    'km':     1.0e5,
    'au':     1.495978707e13,
    'pc':     3.0856775814913673e18,
    'kpc':    3.0856775814913673e21,
    'Mpc':    3.0856775814913673e24,
    'Rsun':   6.957e10,
}


# ---------------------------------------------------------------------------
# Voigt function (Faddeeva-based) -- voigt(0, a) is needed for tau0 -> rhokap
# ---------------------------------------------------------------------------
def voigt(x: float, a: float) -> float:
    """Normalized Voigt profile (same as voigt_mod.f90:voigt).

    voigt(x, a) = Re(w(x + i a)) / sqrt(pi), where w is the Faddeeva
    function. voigt(0, 0) = 1/sqrt(pi).
    """
    from scipy.special import wofz
    z = complex(x, a)
    return float(np.real(wofz(z))) / np.sqrt(PI)


# ---------------------------------------------------------------------------
# Radial profile shapes (mirror clump_mod.f90:profile_shape)
# ---------------------------------------------------------------------------
def profile_shape(name: str, alpha: float, r0: float, r: np.ndarray) -> np.ndarray:
    name = (name or 'constant').strip().lower()
    if name in ('constant', ''):
        return np.ones_like(r)
    if name in ('powerlaw', 'power_law'):
        if r0 <= 0.0:
            return np.ones_like(r)
        r_floor = 0.05 * r0
        r_eff   = np.maximum(r, r_floor)
        r0_eff  = max(r0, r_floor)
        return (r_eff / r0_eff) ** (-alpha)
    if name == 'gaussian':
        if r0 <= 0.0:
            return np.ones_like(r)
        return np.exp(-(r / r0) ** 2)
    if name == 'exponential':
        if r0 <= 0.0:
            return np.ones_like(r)
        return np.exp(-r / r0)
    raise ValueError(f"Unknown profile shape: {name!r}")


# ---------------------------------------------------------------------------
# Line atomic helpers
# ---------------------------------------------------------------------------
def line_reference(line_id: str, temperature: float) -> Dict[str, float]:
    """Return reference Voigt / Doppler / cross-section for the given line
    at the given temperature, in the same units as cl_*_ref in clump_mod.f90.
    """
    if line_id not in LINE_DATA:
        raise ValueError(
            f"Unknown line_id {line_id!r}. Supported: {sorted(LINE_DATA)}"
        )
    L = LINE_DATA[line_id]
    cross0  = SIGMA_0 / np.sqrt(PI) * L['f12']              # cm^2 Hz
    vtherm1 = VTHERM1_AMU / np.sqrt(L['mass_amu'])          # km/s @ 1 K
    vtherm  = vtherm1 * np.sqrt(temperature)                # km/s
    Dfreq   = vtherm / (L['wavelength0'] * UM2KM)           # Hz
    voigt_a = (L['damping'] / FOURPI) / Dfreq               # dimensionless
    return dict(
        wavelength0=L['wavelength0'],
        damping=L['damping'],
        f12=L['f12'],
        mass_amu=L['mass_amu'],
        cross0=cross0,
        vtherm1=vtherm1,
        vtherm=vtherm,
        Dfreq=Dfreq,
        voigt_a=voigt_a,
    )


# ---------------------------------------------------------------------------
# Radial CDF table for non-uniform sampling
#   P(r) ∝ shape_number(r) * r^2  on [r_min_clump, sphere_R].
# ---------------------------------------------------------------------------
class RadialProfile:
    NPROF = 4001

    def __init__(self, args, sphere_R: float, r_min_clump: float, base_radius: float):
        self.r       = np.linspace(0.0, sphere_R, self.NPROF)
        self.r_min   = r_min_clump

        # shape factors on the radial grid
        self.shape_number  = profile_shape(args.clump_number_profile,
                                           args.clump_number_alpha,
                                           args.clump_number_r0 or sphere_R, self.r)
        self.shape_radius  = profile_shape(args.clump_radius_profile,
                                           args.clump_radius_alpha,
                                           args.clump_radius_r0 or sphere_R, self.r)
        self.shape_density = profile_shape(args.clump_density_profile,
                                           args.clump_density_alpha,
                                           args.clump_density_r0 or sphere_R, self.r)

        # inner-cavity mask on the radial number density (matches
        # clump_mod.f90:341 -- only shape_number is zeroed inside r_min)
        self.shape_number = np.where(self.r < r_min_clump, 0.0, self.shape_number)

        integrand = self.shape_number * self.r * self.r
        # cumulative trapezoidal
        cdf = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1])
                                               * (self.r[1:] - self.r[:-1]))])
        total = cdf[-1]
        if total > 0.0:
            cdf /= total
        else:
            # degenerate case: fall back to uniform-in-volume
            cdf = (self.r / sphere_R) ** 3
        self.cdf = cdf

        # cl_radius_max: largest clump radius that the profile may produce
        rcl_local = base_radius * self.shape_radius
        rcl_local = np.where(self.shape_number > 0, rcl_local, 0.0)
        self.r_cl_max = float(max(base_radius, rcl_local.max()))

    # ---- sampling helpers (vectorized over n) ----
    def sample_r(self, n: int, rng: np.random.Generator) -> np.ndarray:
        u = rng.random(n)
        return np.interp(u, self.cdf, self.r)

    def eval_radius(self, r: np.ndarray, base: float) -> np.ndarray:
        return base * np.interp(r, self.r, self.shape_radius)

    def eval_density(self, r: np.ndarray) -> np.ndarray:
        return np.interp(r, self.r, self.shape_density)

    # ---- analytic integrals used for N <-> f_vol / f_cov conversion ----
    def total_count(self, A: float) -> float:
        # ∫ A * shape_number(r) * 4πr^2 dr  (total number, unnormalized)
        integrand = A * self.shape_number * FOURPI * self.r * self.r
        return float(np.trapz(integrand, self.r))

    def f_vol(self, A: float, sphere_R: float) -> float:
        # ∫ A * shape_number * (4π/3) r_cl(r)^3 * 4πr^2 dr / V_shell
        rcl = self.r_cl_max  # use shape_radius
        # use the actual local r_cl(r)
        rcl_r = self.shape_radius   # base factor; multiply by base outside? No -- this stores shape only
        # we'll be more careful: f_vol is sum of clump volumes / shell volume.
        # numerical: N(r) dr = A * shape_number(r) * 4πr^2 dr; V_clump(r) = (4π/3) (base * shape_radius)^3
        # but build_radial_profile_tables only uses *shape*, not base. We carry base via the multiplication.
        raise NotImplementedError  # see helper below

    def f_cov_LOS(self, A: float, base_radius: float, sphere_R: float, r_min_clump: float) -> float:
        # LOS covering factor (radial sightline):
        #   f_cov = ∫_{rmin}^{R} A * shape_number * π * (base * shape_radius)^2 dr
        rcl_local = base_radius * self.shape_radius
        mask = self.r >= r_min_clump
        integrand = A * self.shape_number * PI * rcl_local ** 2
        integrand = np.where(mask, integrand, 0.0)
        return float(np.trapz(integrand, self.r))

    def vol_int(self, A: float, base_radius: float, sphere_R: float, r_min_clump: float) -> float:
        # ∫ N(r) V_clump(r) dr = ∫ A * shape_number(r) * 4πr^2 * (4π/3)(base*shape_radius)^3 dr
        rcl_local = base_radius * self.shape_radius
        mask = self.r >= r_min_clump
        integrand = A * self.shape_number * FOURPI * self.r * self.r \
                    * (FOURPI / 3.0) * rcl_local ** 3
        integrand = np.where(mask, integrand, 0.0)
        return float(np.trapz(integrand, self.r))


# ---------------------------------------------------------------------------
# RSA placement: numba-accelerated CSR linked-list hash grid
# ---------------------------------------------------------------------------
def _hash_grid_params(sphere_R: float, n_expected: int, r_cl_max: float
                      ) -> Tuple[int, float]:
    """Pick a (rg, rg_cell) such that cell size >= 2 * r_cl_max.

    Same heuristic as the legacy ``HashGrid`` class: target ~N^(1/3)
    cells per side, clamped to [32, 512], then enforce the minimum
    cell size so the 27-neighbor stencil captures every overlap.
    """
    rg = min(512, max(32, int(round(n_expected ** (1.0 / 3.0))) + 1))
    rg_cell = max(2.0 * sphere_R / rg, 2.0 * r_cl_max)
    rg = max(2, int((2.0 * sphere_R) / rg_cell) + 1)
    rg_cell = (2.0 * sphere_R) / rg
    return rg, rg_cell


@njit(cache=True, fastmath=True)
def _rsa_kernel(N_target: int,
                xb: np.ndarray, yb: np.ndarray, zb: np.ndarray,
                rcl_batch: np.ndarray, ok_geom: np.ndarray,
                cl_x: np.ndarray, cl_y: np.ndarray, cl_z: np.ndarray,
                cl_r: np.ndarray,
                head: np.ndarray, nxt: np.ndarray,
                rg: int, rg_cell: float, sphere_R: float,
                placed_start: int,
                profiles_active: bool, min_sep2_uniform: float,
                allow_overlap: bool) -> int:
    """Place candidates from one batch into the linked-list hash grid.

    Returns the new ``placed`` count. Mutates ``cl_x/y/z/r``, ``head``,
    ``nxt`` in place.
    """
    placed = placed_start
    n_batch = xb.shape[0]
    inv_cell = 1.0 / rg_cell
    rgm1 = rg - 1
    rg2  = rg * rg

    for k in range(n_batch):
        if placed >= N_target:
            return placed
        if not ok_geom[k]:
            continue
        x = xb[k]; y = yb[k]; z = zb[k]; rcl = rcl_batch[k]

        ig = int((x + sphere_R) * inv_cell)
        if ig < 0: ig = 0
        elif ig > rgm1: ig = rgm1
        jg = int((y + sphere_R) * inv_cell)
        if jg < 0: jg = 0
        elif jg > rgm1: jg = rgm1
        kg = int((z + sphere_R) * inv_cell)
        if kg < 0: kg = 0
        elif kg > rgm1: kg = rgm1

        overlap = False
        if not allow_overlap:
            kg_lo = kg - 1 if kg > 0 else 0
            kg_hi = kg + 2 if kg + 2 <= rg else rg
            jg_lo = jg - 1 if jg > 0 else 0
            jg_hi = jg + 2 if jg + 2 <= rg else rg
            ig_lo = ig - 1 if ig > 0 else 0
            ig_hi = ig + 2 if ig + 2 <= rg else rg
            for kk in range(kg_lo, kg_hi):
                for jj in range(jg_lo, jg_hi):
                    for ii in range(ig_lo, ig_hi):
                        cidx = kk * rg2 + jj * rg + ii
                        idx = head[cidx]
                        while idx >= 0:
                            dx = x - cl_x[idx]
                            dy = y - cl_y[idx]
                            dz = z - cl_z[idx]
                            d2 = dx * dx + dy * dy + dz * dz
                            if profiles_active:
                                sep = rcl + cl_r[idx]
                                if d2 < sep * sep:
                                    overlap = True
                                    break
                            else:
                                if d2 < min_sep2_uniform:
                                    overlap = True
                                    break
                            idx = nxt[idx]
                        if overlap: break
                    if overlap: break
                if overlap: break

        if overlap:
            continue

        cl_x[placed] = x
        cl_y[placed] = y
        cl_z[placed] = z
        cl_r[placed] = rcl
        cidx = kg * rg2 + jg * rg + ig
        nxt[placed]  = head[cidx]
        head[cidx]   = placed
        placed += 1

    return placed


# ---------------------------------------------------------------------------
# Position sampling (vectorized batches)
# ---------------------------------------------------------------------------
def sample_positions_uniform(n: int, r_min: float, r_max: float,
                             cos_cone: float, rng: np.random.Generator,
                             ) -> Tuple[np.ndarray, np.ndarray, np.ndarray,
                                        np.ndarray, np.ndarray, np.ndarray]:
    """Direct sampling: r in shell, cos(theta) in cone window, phi uniform.

    Returns (x, y, z, r, cos_theta, sin_theta).
    """
    # r^3 uniform between rmin^3 and rmax^3 -> uniform volume
    u = rng.random(n)
    r3 = r_min ** 3 + (r_max ** 3 - r_min ** 3) * u
    r = r3 ** (1.0 / 3.0)

    if cos_cone > 0.0:
        cos_t = cos_cone + (1.0 - cos_cone) * rng.random(n)
        sign = rng.random(n) < 0.5
        cos_t = np.where(sign, -cos_t, cos_t)
    else:
        cos_t = 2.0 * rng.random(n) - 1.0
    sin_t = np.sqrt(np.maximum(0.0, 1.0 - cos_t * cos_t))
    phi = TWOPI * rng.random(n)
    x = r * sin_t * np.cos(phi)
    y = r * sin_t * np.sin(phi)
    z = r * cos_t
    return x, y, z, r, cos_t, sin_t


def sample_positions_profile(n: int, prof: RadialProfile, cos_cone: float,
                             rng: np.random.Generator,
                             ) -> Tuple[np.ndarray, np.ndarray, np.ndarray,
                                        np.ndarray, np.ndarray, np.ndarray]:
    r = prof.sample_r(n, rng)
    if cos_cone > 0.0:
        cos_t = cos_cone + (1.0 - cos_cone) * rng.random(n)
        sign = rng.random(n) < 0.5
        cos_t = np.where(sign, -cos_t, cos_t)
    else:
        cos_t = 2.0 * rng.random(n) - 1.0
    sin_t = np.sqrt(np.maximum(0.0, 1.0 - cos_t * cos_t))
    phi = TWOPI * rng.random(n)
    x = r * sin_t * np.cos(phi)
    y = r * sin_t * np.sin(phi)
    z = r * cos_t
    return x, y, z, r, cos_t, sin_t


def cone_fully_inside_ok(r: np.ndarray, cos_t: np.ndarray, sin_t: np.ndarray,
                         rcl: np.ndarray, cos_cone: float) -> np.ndarray:
    """Vectorized cone fully-inside test.

    cos(theta_min) = |cos(theta)| * sqrt(1 - (r_cl/r)^2) - sin(theta) * (r_cl/r)
    A candidate is accepted if cos(theta_min) >= cos_cone.

    When r <= rcl, the clump straddles the origin; accept unconditionally
    (matches the Fortran "if (r_trial > rcl_trial)" guard).
    """
    if cos_cone <= 0.0:
        return np.ones_like(r, dtype=bool)
    ok = np.ones_like(r, dtype=bool)
    big = r > rcl
    ratio = np.where(big, rcl / np.maximum(r, 1e-300), 0.0)
    cos_min = np.abs(cos_t) * np.sqrt(np.maximum(0.0, 1.0 - ratio * ratio)) \
              - sin_t * ratio
    return np.where(big, cos_min >= cos_cone, ok)


# ---------------------------------------------------------------------------
# Velocity assignment (mirror assign_clump_velocities_from_type)
# ---------------------------------------------------------------------------
def assign_velocity_systematic(args, x: np.ndarray, y: np.ndarray, z: np.ndarray,
                               sphere_R: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    rr = np.sqrt(x * x + y * y + z * z)
    rcyl = np.sqrt(x * x + y * y)
    vx = np.zeros_like(x);  vy = np.zeros_like(x);  vz = np.zeros_like(x)
    vt = (args.velocity_type or '').strip().lower()

    if vt == 'hubble':
        vx = args.Vexp * x / sphere_R
        vy = args.Vexp * y / sphere_R
        vz = args.Vexp * z / sphere_R
    elif vt == 'constant_radial':
        m = rr > 0
        vx[m] = args.Vexp * x[m] / rr[m]
        vy[m] = args.Vexp * y[m] / rr[m]
        vz[m] = args.Vexp * z[m] / rr[m]
    elif vt == 'power_law':
        m = rr > 0
        scale = args.Vexp * (rr[m] / sphere_R) ** args.velocity_alpha
        vx[m] = scale * x[m] / rr[m]
        vy[m] = scale * y[m] / rr[m]
        vz[m] = scale * z[m] / rr[m]
    elif vt == 'linear_decelerate':
        m = rr > 0
        denom = sphere_R - max(args.rmin if args.rmin is not None else 0.0, 0.0)
        scale = args.Vexp * np.maximum(0.0, (sphere_R - rr[m]) / max(denom, 1e-300))
        vx[m] = scale * x[m] / rr[m]
        vy[m] = scale * y[m] / rr[m]
        vz[m] = scale * z[m] / rr[m]
    elif vt == 'parallel_velocity':
        vx[:] = args.Vx;  vy[:] = args.Vy;  vz[:] = args.Vz
    elif vt == 'ssh':
        # Song, Seon & Hwang (2020)
        m = rr > 0
        inside = (rr < args.rpeak) & m
        outside = (rr >= args.rpeak) & m
        scale_i = args.Vpeak / args.rpeak
        vx[inside] = scale_i * x[inside]
        vy[inside] = scale_i * y[inside]
        vz[inside] = scale_i * z[inside]
        scale_o = args.Vpeak + args.DeltaV * (rr[outside] - args.rpeak) \
                  / (sphere_R - args.rpeak)
        vx[outside] = scale_o * x[outside] / rr[outside]
        vy[outside] = scale_o * y[outside] / rr[outside]
        vz[outside] = scale_o * z[outside] / rr[outside]
    elif vt == 'rotating_solid_body':
        vx = -args.Vrot * y / sphere_R
        vy =  args.Vrot * x / sphere_R
    elif vt == 'rotating_galaxy_halo':
        m_in  = (rcyl < args.rinner) & (rcyl > 0)
        m_out = rcyl >= args.rinner
        vx[m_in]  = -args.Vrot * y[m_in]  / args.rinner
        vy[m_in]  =  args.Vrot * x[m_in]  / args.rinner
        vx[m_out] = -args.Vrot * y[m_out] / rcyl[m_out]
        vy[m_out] =  args.Vrot * x[m_out] / rcyl[m_out]
    elif vt in ('', 'none'):
        pass
    else:
        print(f'WARNING: unknown velocity_type {vt!r}; no systematic flow applied.')

    return vx, vy, vz


# ---------------------------------------------------------------------------
# RSA driver
# ---------------------------------------------------------------------------
def derive_N_clumps(args, sphere_R: float, r_min_clump: float,
                    base_radius: float, prof: RadialProfile | None):
    """Convert the user's normalization choice to N_clumps + A_norm.

    Returns (N_clumps, A_norm, f_vol_realized, f_cov_realized).
    For the uniform case A_norm is None.
    """
    if prof is None:
        if args.clump_N_clumps and args.clump_N_clumps > 0:
            N = int(args.clump_N_clumps)
        elif args.clump_f_vol and args.clump_f_vol > 0:
            N = int(round(args.clump_f_vol * (sphere_R**3 - r_min_clump**3)
                          / base_radius**3))
        elif args.clump_f_cov and args.clump_f_cov > 0:
            N = int(round((4.0 / 3.0) * args.clump_f_cov
                          * (sphere_R**2 + sphere_R * r_min_clump + r_min_clump**2)
                          / base_radius**2))
        else:
            sys.exit('ERROR: specify --clump_N_clumps, --clump_f_vol, or --clump_f_cov')
        f_vol = N * base_radius**3 / max(sphere_R**3 - r_min_clump**3, 1e-300)
        f_cov = 0.75 * N * base_radius**2 \
                / max(sphere_R**2 + sphere_R * r_min_clump + r_min_clump**2, 1e-300)
        return N, None, f_vol, f_cov

    # ---- profile case ----
    # total_count(A=1) integrates A * shape_number * 4πr^2 (analogue of f90 routine)
    tc1 = prof.total_count(1.0)
    if args.clump_N_clumps and args.clump_N_clumps > 0:
        N = int(args.clump_N_clumps)
        A = N / max(tc1, 1e-300)
    elif args.clump_f_vol and args.clump_f_vol > 0:
        v_int1 = prof.vol_int(1.0, base_radius, sphere_R, r_min_clump)
        V_shell = (FOURPI / 3.0) * (sphere_R**3 - r_min_clump**3)
        fvol_unit = v_int1 / max(V_shell, 1e-300)
        A = args.clump_f_vol / max(fvol_unit, 1e-300)
        N = int(round(prof.total_count(A)))
    elif args.clump_f_cov and args.clump_f_cov > 0:
        fcov_unit = prof.f_cov_LOS(1.0, base_radius, sphere_R, r_min_clump)
        A = args.clump_f_cov / max(fcov_unit, 1e-300)
        N = int(round(prof.total_count(A)))
    else:
        sys.exit('ERROR: specify --clump_N_clumps, --clump_f_vol, or --clump_f_cov')
    v_int = prof.vol_int(A, base_radius, sphere_R, r_min_clump)
    V_shell = (FOURPI / 3.0) * (sphere_R**3 - r_min_clump**3)
    f_vol = v_int / max(V_shell, 1e-300)
    f_cov = prof.f_cov_LOS(A, base_radius, sphere_R, r_min_clump)
    return N, A, f_vol, f_cov


def derive_rhokap(args, line_ref: Dict[str, float], base_radius: float,
                  sphere_R: float, r_min_clump: float, N: int, A: float | None,
                  prof: RadialProfile | None, distance2cm: float) -> float:
    """Back-solve the peak (reference) rhokap from the opacity input mode."""
    voigt0 = voigt(0.0, line_ref['voigt_a'])
    if args.clump_tau0 and args.clump_tau0 > 0:
        return args.clump_tau0 / (voigt0 * base_radius)
    if args.clump_NHI and args.clump_NHI > 0:
        return args.clump_NHI * line_ref['cross0'] / (line_ref['Dfreq'] * base_radius)
    if args.clump_nH and args.clump_nH > 0:
        if distance2cm <= 0.0:
            sys.exit('ERROR: --clump_nH requires --distance_unit')
        return args.clump_nH * line_ref['cross0'] * distance2cm / line_ref['Dfreq']
    if (args.taumax and args.taumax > 0) or (args.N_HImax and args.N_HImax > 0):
        if prof is None:
            GF = N * base_radius**3 / max(sphere_R**2 + sphere_R*r_min_clump + r_min_clump**2, 1e-300)
        else:
            # LOS geometric factor for the profile case:
            # GF = ∫_{rmin}^{R} A * shape_number(r) * (4/3) r_cl(r)^3 dr  (radial sightline)
            rcl_local = base_radius * prof.shape_radius
            mask = prof.r >= r_min_clump
            integrand = A * prof.shape_number * (4.0/3.0) * rcl_local**3
            integrand = np.where(mask, integrand, 0.0)
            GF = float(np.trapz(integrand, prof.r))
        if GF <= 0:
            sys.exit('ERROR: cannot back-solve opacity (geometric factor is zero)')
        if args.taumax and args.taumax > 0:
            return args.taumax / (GF * voigt0)
        return args.N_HImax * line_ref['cross0'] / (GF * line_ref['Dfreq'])
    sys.exit('ERROR: specify one of --clump_tau0, --clump_NHI, --clump_nH, --taumax, --N_HImax')


def generate_clumps(args, sphere_R: float, r_min_clump: float,
                    base_radius: float, N: int, A: float | None,
                    prof: RadialProfile | None, cos_cone: float,
                    rng: np.random.Generator) -> Dict[str, np.ndarray]:
    """Place N clumps via RSA. Returns dict of column arrays."""
    cl_x  = np.zeros(N, dtype=np.float64)
    cl_y  = np.zeros(N, dtype=np.float64)
    cl_z  = np.zeros(N, dtype=np.float64)
    cl_r  = np.full(N, base_radius if prof is None else prof.r_cl_max, dtype=np.float64)

    if prof is None:
        r_cl_max = base_radius
    else:
        r_cl_max = prof.r_cl_max

    if args.clump_fully_inside:
        if r_cl_max >= sphere_R:
            sys.exit('ERROR: --clump_fully_inside but max clump radius >= sphere_R')
        if r_min_clump + 2.0 * r_cl_max > sphere_R:
            sys.exit(f'ERROR: --clump_fully_inside but rmin + 2*r_cl_max > rmax '
                     f'(rmin={r_min_clump:g}, r_cl_max={r_cl_max:g}, rmax={sphere_R:g})')

    # shell limits for uniform path
    if args.clump_fully_inside:
        r_max_center = sphere_R - base_radius
        r_min_center = r_min_clump + base_radius
    else:
        r_max_center = sphere_R
        r_min_center = r_min_clump

    rg, rg_cell = _hash_grid_params(sphere_R, N, r_cl_max)
    head = np.full(rg * rg * rg, -1, dtype=np.int64)
    nxt  = np.full(N,             -1, dtype=np.int64)
    min_sep2_uniform = (2.0 * base_radius) ** 2
    profiles_active = prof is not None
    allow_overlap = bool(args.clump_allow_overlap)

    print(f' RSA: grid {rg}^3, placing {N} clumps...  '
          f'(numba={"on" if _HAVE_NUMBA else "off (fallback)"})')
    print(f' RSA: clump_fully_inside = {args.clump_fully_inside}')
    if r_min_clump > 0:
        print(f' RSA: shell rmin/rmax    = {r_min_clump:.5f} {sphere_R:.5f}')
    if cos_cone > 0:
        print(f' RSA: cone_opening       = {args.cone_opening:.3f} deg')

    placed = 0
    attempts = 0
    batch = max(1024, min(4 * N, 200000))
    t0 = time.time()
    next_report = 100000

    while placed < N:
        attempts += batch

        if prof is None:
            xb, yb, zb, rb, ct, st = sample_positions_uniform(
                batch, r_min_center, r_max_center, cos_cone, rng)
            rcl_batch = np.full(batch, base_radius, dtype=np.float64)
            ok_geom = np.ones(batch, dtype=np.bool_)
            if args.clump_fully_inside and cos_cone > 0:
                ok_geom = cone_fully_inside_ok(rb, ct, st, rcl_batch, cos_cone)
        else:
            xb, yb, zb, rb, ct, st = sample_positions_profile(batch, prof, cos_cone, rng)
            rcl_batch = prof.eval_radius(rb, base_radius).astype(np.float64, copy=False)
            ok_geom = np.ones(batch, dtype=np.bool_)
            if args.clump_fully_inside:
                ok_geom &= (rb + rcl_batch <= sphere_R)
                ok_geom &= (rb - rcl_batch >= r_min_clump)
                if cos_cone > 0:
                    ok_geom &= cone_fully_inside_ok(rb, ct, st, rcl_batch, cos_cone)

        new_placed = _rsa_kernel(
            N, xb, yb, zb, rcl_batch, ok_geom,
            cl_x, cl_y, cl_z, cl_r,
            head, nxt, rg, rg_cell, sphere_R,
            placed, profiles_active, min_sep2_uniform, allow_overlap)

        if new_placed >= next_report:
            dt = time.time() - t0
            print(f'   placed {new_placed:>10d} / {N}   '
                  f'(attempts={attempts}, {new_placed/max(dt,1e-9):.1e} clumps/s)')
            next_report = ((new_placed // 100000) + 1) * 100000

        if new_placed == placed:
            # No progress in this batch -- rare but possible if cone+overlap
            # acceptance is extremely low. Continue; outer while loop will
            # try another batch.
            if attempts > 200 * max(N, 1):
                sys.exit(f'ERROR: RSA stalled at {new_placed}/{N} after '
                         f'{attempts} attempts; loosen geometry or N.')
        placed = new_placed

    acc = N / max(attempts, 1) * 100.0
    if args.clump_allow_overlap:
        print(' Random placement done (overlap allowed).')
    else:
        print(f' RSA done, acceptance rate = {acc:.1f}%')

    return dict(X=cl_x, Y=cl_y, Z=cl_z, R_CLUMP=cl_r)


# ---------------------------------------------------------------------------
# Clump physics (rhokap, voigt_a, Dfreq, vtherm, temperature)
# ---------------------------------------------------------------------------
def assign_perclump_physics(cols: Dict[str, np.ndarray], prof: RadialProfile | None,
                            ref: Dict[str, float], rhokap_ref: float,
                            args, base_radius: float) -> None:
    N = len(cols['X'])
    r = np.sqrt(cols['X']**2 + cols['Y']**2 + cols['Z']**2)
    if prof is None:
        cols['R_CLUMP'] = np.full(N, base_radius)
        cols['RHOKAP']  = np.full(N, rhokap_ref)
        cols['TEMP']    = np.full(N, args.temperature if args.clump_temperature < 0
                                  else args.clump_temperature)
        cols['_VOIGT_A'] = np.full(N, ref['voigt_a'])
        cols['_DFREQ']   = np.full(N, ref['Dfreq'])
        cols['_VTHERM']  = np.full(N, ref['vtherm'])
    else:
        # clump radius from shape_radius * base
        rcl = prof.eval_radius(r, base_radius)
        cols['R_CLUMP'] = rcl
        # temperature: clump_mod uses ref T uniformly when no tab_temperature
        # (file-based clump T is not supported here).
        T0 = args.temperature if args.clump_temperature < 0 else args.clump_temperature
        T  = np.full(N, T0)
        cols['TEMP'] = T
        vtherm = (VTHERM1_AMU / np.sqrt(LINE_DATA[args.line_id]['mass_amu'])) * np.sqrt(T)
        cols['_VTHERM'] = vtherm
        Dfreq = vtherm / (ref['wavelength0'] * UM2KM)
        cols['_DFREQ'] = Dfreq
        cols['_VOIGT_A'] = (ref['damping'] / FOURPI) / Dfreq
        # local opacity: rhokap = base_rhokap * shape_density(r) * (Dfreq_ref / Dfreq)
        dens = prof.eval_density(r)
        cols['RHOKAP'] = rhokap_ref * dens * (ref['Dfreq'] / Dfreq)


# ---------------------------------------------------------------------------
# FITS / HDF5 writer (matches clump_mod.f90:write_clumps_info schema)
# ---------------------------------------------------------------------------
def write_clumps_file(path: str, cols: Dict[str, np.ndarray],
                      header: Dict[str, object]) -> None:
    # decide format from extension
    p_low = path.lower()
    is_h5 = p_low.endswith('.h5') or p_low.endswith('.hdf5')

    # decide which optional clump columns to keep (const_tol = 1e-3,
    # matching the Fortran writer).
    const_tol = 1.0e-3
    def varying(arr):
        amax, amin = float(arr.max()), float(arr.min())
        amean = float(arr.mean()) if arr.size > 0 else 0.0
        return (amax - amin) > const_tol * max(abs(amean), 1e-300)

    write_radius = varying(cols['R_CLUMP'])
    write_rhokap = varying(cols['RHOKAP'])
    write_temp   = varying(cols['TEMP'])

    # Fortran serializes columns as float32 (bitpix=-32).
    def to32(arr): return np.asarray(arr, dtype=np.float32)

    out_cols = dict(
        X  = to32(cols['X']),
        Y  = to32(cols['Y']),
        Z  = to32(cols['Z']),
        VX = to32(cols['VX']),
        VY = to32(cols['VY']),
        VZ = to32(cols['VZ']),
    )
    if write_radius: out_cols['R_CLUMP'] = to32(cols['R_CLUMP'])
    if write_rhokap: out_cols['RHOKAP']  = to32(cols['RHOKAP'])
    if write_temp:   out_cols['TEMP']    = to32(cols['TEMP'])

    if is_h5:
        _write_hdf5(path, out_cols, header)
    else:
        _write_fits(path, out_cols, header)

    extra = {'R_CLUMP': write_radius, 'RHOKAP': write_rhokap, 'TEMP': write_temp}
    print(' Clumps: clump columns saved (R_CLUMP/RHOKAP/TEMP) = '
          + ' '.join(str(extra[k]) for k in ('R_CLUMP', 'RHOKAP', 'TEMP')))
    print(f' Clumps saved to {path}')


def _write_fits(path: str, cols: Dict[str, np.ndarray],
                header: Dict[str, object]) -> None:
    from astropy.io import fits
    # build BinTable
    fits_cols = [fits.Column(name=k, array=arr, format='E')
                 for k, arr in cols.items()]
    hdu_table = fits.BinTableHDU.from_columns(fits_cols, name='Clumps')
    # populate header keywords
    for k, v in header.items():
        if v is None: continue
        try:
            if isinstance(v, str):
                hdu_table.header[k[:8].upper()] = v
            else:
                hdu_table.header[k[:8].upper()] = v
        except Exception:
            pass
    hdul = fits.HDUList([fits.PrimaryHDU(), hdu_table])
    if os.path.exists(path):
        os.remove(path)
    hdul.writeto(path)


def _write_hdf5(path: str, cols: Dict[str, np.ndarray],
                header: Dict[str, object]) -> None:
    import h5py
    if os.path.exists(path):
        os.remove(path)
    with h5py.File(path, 'w', libver='latest', track_order=True) as f:
        g = f.create_group('Clumps', track_order=True)
        for k, arr in cols.items():
            kwargs = {}
            if arr.size > 4096:
                kwargs.update(dict(chunks=(min(arr.size, 4096),),
                                   compression='gzip', compression_opts=4))
            g.create_dataset(k, data=arr, **kwargs)
        for k, v in header.items():
            if v is None: continue
            g.attrs[k.upper()] = v


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description='Python port of clump_mod.f90 / make_clumps.x',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # geometry
    p.add_argument('--rmax', type=float, required=True,
                   help='outer sphere radius [code units]')
    p.add_argument('--rmin', type=float, default=0.0,
                   help='inner placement radius [code units]')
    p.add_argument('--cone_opening', type=float, default=0.0,
                   help='biconical half-opening angle [deg]; 0 = full sphere')
    p.add_argument('--distance_unit', type=str, default='',
                   choices=list(DIST_UNIT_CM.keys()),
                   help='unit of rmax/rmin (only needed for --clump_nH)')

    # clump radius and counts
    p.add_argument('--clump_radius', type=float, required=True,
                   help='base clump radius [code units]')
    p.add_argument('--clump_N_clumps', type=float, default=0,
                   help='specify exact number of clumps (one of: N_clumps, f_vol, f_cov)')
    p.add_argument('--clump_f_vol', type=float, default=0,
                   help='target volume filling factor')
    p.add_argument('--clump_f_cov', type=float, default=0,
                   help='target covering factor (radial sightline)')

    # opacity (one of)
    p.add_argument('--clump_tau0', type=float, default=0,
                   help='clump line-center tau (center to surface)')
    p.add_argument('--clump_NHI',  type=float, default=0,
                   help='clump HI column density [cm^-2]')
    p.add_argument('--clump_nH',   type=float, default=0,
                   help='clump HI number density [cm^-3]; requires --distance_unit')
    p.add_argument('--taumax',     type=float, default=0,
                   help='system-level radial sightline tau (back-solves rhokap)')
    p.add_argument('--N_HImax',    type=float, default=0,
                   help='system-level HI column [cm^-2] (back-solves rhokap)')

    # thermal
    p.add_argument('--temperature', type=float, default=1.0e4,
                   help='global temperature [K] (used for clumps unless --clump_temperature set)')
    p.add_argument('--clump_temperature', type=float, default=-1.0,
                   help='clump temperature [K]; <0 means use --temperature')
    p.add_argument('--line_id', type=str, default='ly_alpha',
                   choices=list(LINE_DATA.keys()),
                   help='atomic line for cross-section / Doppler factors')

    # RSA controls
    p.add_argument('--clump_fully_inside', action=argparse.BooleanOptionalAction,
                   default=True,
                   help='reject clumps protruding past rmax, into rmin, or outside the cone')
    p.add_argument('--clump_allow_overlap', action='store_true', default=False,
                   help='skip overlap rejection (LaRT will use overlap-aware raytrace)')

    # systematic velocity field
    p.add_argument('--velocity_type', type=str, default='',
                   choices=['', 'none', 'hubble', 'constant_radial', 'power_law',
                            'linear_decelerate', 'parallel_velocity', 'ssh',
                            'rotating_solid_body', 'rotating_galaxy_halo'],
                   help='systematic velocity field type')
    p.add_argument('--Vexp', type=float, default=0.0, help='expansion speed [km/s]')
    p.add_argument('--Vrot', type=float, default=0.0, help='rotation speed [km/s]')
    p.add_argument('--velocity_alpha', type=float, default=1.0,
                   help='exponent for power_law velocity field')
    p.add_argument('--rinner', type=float, default=0.0,
                   help='inner radius for rotating_galaxy_halo [code units]')
    p.add_argument('--rpeak',  type=float, default=0.0,
                   help='SSH model peak radius [code units]')
    p.add_argument('--Vpeak',  type=float, default=0.0,
                   help='SSH model peak velocity [km/s]')
    p.add_argument('--DeltaV', type=float, default=0.0,
                   help='SSH model outer linear delta [km/s]')
    p.add_argument('--Vx', type=float, default=0.0)
    p.add_argument('--Vy', type=float, default=0.0)
    p.add_argument('--Vz', type=float, default=0.0)
    p.add_argument('--clump_sigma_v', type=float, default=0.0,
                   help='Gaussian sigma of random clump bulk velocity [km/s]')

    # radial profile knobs (for each axis: 'constant' | 'powerlaw' | 'gaussian' | 'exponential')
    p.add_argument('--clump_radius_profile',  type=str, default='constant')
    p.add_argument('--clump_density_profile', type=str, default='constant')
    p.add_argument('--clump_number_profile',  type=str, default='constant')
    p.add_argument('--clump_radius_alpha',  type=float, default=0.0)
    p.add_argument('--clump_density_alpha', type=float, default=0.0)
    p.add_argument('--clump_number_alpha',  type=float, default=0.0)
    p.add_argument('--clump_radius_r0',  type=float, default=0.0)
    p.add_argument('--clump_density_r0', type=float, default=0.0)
    p.add_argument('--clump_number_r0',  type=float, default=0.0)

    # output
    p.add_argument('-o', '--output', type=str, required=True,
                   help='output filename (.fits, .fits.gz, .h5, .hdf5)')
    p.add_argument('--seed', type=int, default=12345,
                   help='RNG seed for reproducibility')
    return p


def main(argv=None) -> int:
    args = build_argparser().parse_args(argv)
    rng  = np.random.default_rng(args.seed)

    sphere_R    = args.rmax
    r_min_clump = max(0.0, args.rmin)
    base_radius = args.clump_radius

    if base_radius <= 0:
        sys.exit('ERROR: --clump_radius must be > 0')
    if r_min_clump >= sphere_R:
        sys.exit(f'ERROR: --rmin must be < --rmax (got {r_min_clump} >= {sphere_R})')

    cos_cone = -1.0
    if 0.0 < args.cone_opening < 90.0:
        cos_cone = float(np.cos(np.deg2rad(args.cone_opening)))

    distance2cm = DIST_UNIT_CM.get(args.distance_unit, 0.0)

    # profile setup
    profiles_active = any(
        (p or 'constant').strip().lower() != 'constant'
        for p in (args.clump_radius_profile,
                  args.clump_density_profile,
                  args.clump_number_profile)
    )
    prof = RadialProfile(args, sphere_R, r_min_clump, base_radius) \
           if profiles_active else None

    # line atomic data
    temp_cl = args.temperature if args.clump_temperature < 0 else args.clump_temperature
    ref = line_reference(args.line_id, temp_cl)

    # population counts and reference rhokap
    N, A, fvol_real, fcov_real = derive_N_clumps(
        args, sphere_R, r_min_clump, base_radius, prof)
    if N <= 0: N = 1
    rhokap_ref = derive_rhokap(args, ref, base_radius, sphere_R, r_min_clump,
                               N, A, prof, distance2cm)

    print(f' Clumps: N_clumps  = {N}')
    print(f' Clumps: f_vol     = {fvol_real:.6f}')
    print(f' Clumps: f_cov     = {fcov_real:.5f}')
    print(f' Clumps: rmin/rmax = {r_min_clump:.5f} {sphere_R:.5f}')
    print(f' Clumps: cl_rhokap = {rhokap_ref:.4e}')
    print(f' Clumps: voigt_a   = {ref["voigt_a"]:.5f}')
    print(f' Clumps: cl_Dfreq  = {ref["Dfreq"]:.4e}')
    if profiles_active:
        print(f' Clumps: radius profile  = {args.clump_radius_profile}')
        print(f' Clumps: density profile = {args.clump_density_profile}')
        print(f' Clumps: number profile  = {args.clump_number_profile}')

    # RSA placement
    cols = generate_clumps(args, sphere_R, r_min_clump, base_radius,
                           N, A, prof, cos_cone, rng)

    # clump physical params
    assign_perclump_physics(cols, prof, ref, rhokap_ref, args, base_radius)

    # velocities -- systematic + clump random sigma_v.
    vx_sys, vy_sys, vz_sys = assign_velocity_systematic(
        args, cols['X'], cols['Y'], cols['Z'], sphere_R)
    if args.clump_sigma_v > 0:
        vx_sys = vx_sys + rng.normal(0.0, args.clump_sigma_v, size=N)
        vy_sys = vy_sys + rng.normal(0.0, args.clump_sigma_v, size=N)
        vz_sys = vz_sys + rng.normal(0.0, args.clump_sigma_v, size=N)
    cols['VX'] = vx_sys
    cols['VY'] = vy_sys
    cols['VZ'] = vz_sys

    # realized f_vol / f_cov for the header
    if prof is None:
        f_vol_actual = N * base_radius**3 / max(sphere_R**3 - r_min_clump**3, 1e-300)
        f_cov_actual = 0.75 * N * base_radius**2 \
                       / max(sphere_R**2 + sphere_R*r_min_clump + r_min_clump**2, 1e-300)
    else:
        f_vol_actual = float(np.sum(cols['R_CLUMP']**3)) \
                       / max(sphere_R**3 - r_min_clump**3, 1e-300)
        f_cov_actual = fcov_real  # quadrature value (same convention as Fortran)

    # header keywords (mirror write_clumps_info)
    cl_radius_max = base_radius if prof is None else float(prof.r_cl_max)
    header = dict(
        N_CLUMPS = int(N),
        SPHERE_R = float(sphere_R),
        RMIN     = float(r_min_clump),
        CL_RAD   = float(cl_radius_max),
        F_VOL    = float(f_vol_actual),
        F_COV    = float(f_cov_actual),
        TAU0     = float(args.clump_tau0),
        SIGMA_V  = float(args.clump_sigma_v),
        TEMP_CL  = float(temp_cl),
        RHOKAP   = float(rhokap_ref),
        CL_DFREQ = float(ref['Dfreq']),
        VTHERM   = float(ref['vtherm']),
        VOIGT_A  = float(ref['voigt_a']),
        RMAX     = float(args.rmax),
        IN_FCOV  = float(args.clump_f_cov),
        IN_FVOL  = float(args.clump_f_vol),
        IN_NCL   = float(args.clump_N_clumps),
        IN_NHI   = float(args.clump_NHI),
        IN_NH    = float(args.clump_nH),
        IN_TEMP  = float(args.clump_temperature),
        DISTUNIT = args.distance_unit,
        DIST_CM  = float(distance2cm),
        LINE_ID  = args.line_id,
        CONE_OP  = float(args.cone_opening),
    )

    write_clumps_file(args.output, cols, header)

    # diagnostics
    rcl = cols['R_CLUMP']
    print(f' Clumps: realized f_vol      = {f_vol_actual:.4e}')
    print(f' Clumps: realized f_cov      = {f_cov_actual:.4e}')
    print(f' Clumps: r_cl min/mean/max   = {rcl.min():.4e} {rcl.mean():.4e} {rcl.max():.4e}')
    tau_mean = float(np.mean(cols['RHOKAP'] * voigt(0.0, ref['voigt_a']) * rcl))
    print(f' Clumps: mean tau_per_clump  = {tau_mean:.4e}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
