#!/usr/bin/env python3
"""
make_amr_sphere_radial.py
============================
Build a spherical AMR grid whose refinement level increases toward the center.

Refinement is purely geometric (not gradient-driven): each concentric shell
is assigned a fixed AMR level, with the coarsest cells at the outer edge and
the finest cells near the center.  This is useful when the density or
velocity field has radial structure that benefits from a uniform radial
resolution hierarchy independent of local gradients.

Shell spacing modes (--spacing):
  linear  — equal-width shells; each shell has the same radial thickness
  log     — geometrically shrinking shells; finer shells near the center

Boundary refinement (--refine_boundary):
  Forces every cell intersecting the sphere surface to --boundary_level_max
  (defaults to --level_max).  Use --boundary_level_max to set the surface
  resolution independently of the interior level_max.

Velocity laws (--velocity):
  none           — zero velocity everywhere (default)
  hubble         — Hubble-like expansion: v_r = V_exp * r/rmax  (v ∝ r)
  constant_radial — constant radial outflow: v_r = V_exp  (uniform speed)
  power_law      — power-law radial: v_r = V_exp * (r/rmax)^velocity_alpha

  For all laws, (vx,vy,vz) are the Cartesian projections of the radial
  velocity.  Negative V_exp gives inflow.  Velocities are zero outside rmax.
  The velocity is embedded in the grid data file; no par%%velocity_type is
  needed in the LaRT input.

Usage
-----
    # Basic: linear shells, uniform density, boxlen=2
    python make_amr_sphere_radial.py

    # Log spacing with ratio=0.5
    python make_amr_sphere_radial.py --spacing log --ratio 0.5

    # Custom level range and density profile
    python make_amr_sphere_radial.py --level_min 2 --level_max 6 --density gaussian

    # Add sharp sphere-surface refinement (surface = level_max)
    python make_amr_sphere_radial.py --refine_boundary

    # Surface at level 3, interior up to level 5
    python make_amr_sphere_radial.py --refine_boundary --boundary_level_max 3

    # Hubble-flow outflow (V_exp=200 km/s at rmax)
    python make_amr_sphere_radial.py --velocity hubble --vexp 200

    # Constant radial outflow
    python make_amr_sphere_radial.py --velocity constant_radial --vexp 100

    # Power-law v_r = 200*(r/rmax)^2
    python make_amr_sphere_radial.py --velocity power_law --vexp 200 --velocity_alpha 2.0

    # Solid-body rotation about z-axis (Garavito-Camargo+2014)
    python make_amr_sphere_radial.py --velocity rotating_solid_body --vrot 150

    # Flat rotation curve above rinner (Kim et al.)
    python make_amr_sphere_radial.py --velocity rotating_galaxy_halo --vrot 200 --rinner 0.3

    # Compare radial-only vs radial+boundary, save slice plots
    python make_amr_sphere_radial.py --compare --refine_boundary --plot

    # Log spacing + boundary + log-scale density plot
    python make_amr_sphere_radial.py --spacing log --refine_boundary \\
        --boundary_level_max 3 --density gaussian --plot --plot_log

    python make_amr_sphere_radial.py --help
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from AMR_grid import AMRGrid, _OPTIONAL_COLUMNS


# ---------------------------------------------------------------------------
# Density profiles
# ---------------------------------------------------------------------------

def make_density_fn(profile, n0, rmax, sigma=None, r_scale=None,
                    density_alpha=0.0, rmin=0.0, cone_opening=0.0):
    # Precompute cone mask threshold
    if 0.0 < cone_opening < 90.0:
        cos_cone = np.cos(np.radians(cone_opening))
    else:
        cos_cone = None  # no cone mask

    def _apply_cone(dens, x, y, z, r):
        """Zero density outside the bicone (z-axis)."""
        if cos_cone is None:
            return dens
        cos_theta = np.where(r > 0, np.abs(np.asarray(z)) / r, 1.0)
        return np.where(cos_theta >= cos_cone, dens, 0.0)

    if profile == 'uniform':
        def fn(x, y, z):
            r = np.sqrt(np.asarray(x)**2 + np.asarray(y)**2 + np.asarray(z)**2)
            d = np.where((r <= rmax) & (r >= rmin), n0, 0.0)
            return _apply_cone(d, x, y, z, r)
    elif profile == 'gaussian':
        sig = sigma if sigma else rmax * 0.4
        def fn(x, y, z):
            r2 = np.asarray(x)**2 + np.asarray(y)**2 + np.asarray(z)**2
            r  = np.sqrt(r2)
            d = np.where((r <= rmax) & (r >= rmin),
                         n0 * np.exp(-r2 / (2.0 * sig**2)), 0.0)
            return _apply_cone(d, x, y, z, r)
    elif profile == 'exponential':
        rs = r_scale if r_scale else rmax * 0.3
        def fn(x, y, z):
            r = np.sqrt(np.asarray(x)**2 + np.asarray(y)**2 + np.asarray(z)**2)
            d = np.where((r <= rmax) & (r >= rmin), n0 * np.exp(-r / rs), 0.0)
            return _apply_cone(d, x, y, z, r)
    elif profile == 'power_law':
        # n(r) = n0 * (rmax/r)^density_alpha
        # Matches LaRT Cartesian convention in grid_mod_car.f90.
        def fn(x, y, z):
            r = np.sqrt(np.asarray(x)**2 + np.asarray(y)**2 + np.asarray(z)**2)
            safe_r = np.where(r > 0, r, 1.0)
            dens = n0 * (rmax / safe_r) ** density_alpha
            d = np.where((r <= rmax) & (r >= rmin) & (r > 0), dens, 0.0)
            return _apply_cone(d, x, y, z, r)
    else:
        raise ValueError(f'Unknown density profile: {profile!r}')
    return fn


# ---------------------------------------------------------------------------
# Velocity laws
# ---------------------------------------------------------------------------

def make_velocity_fn(law, rmax, v_exp=0.0, velocity_alpha=1.0, vrot=0.0, rinner=None,
                     rmin=0.0):
    """
    Return a vectorised velocity function  fn(x, y, z) -> (vx, vy, vz) [km/s].
    Velocity is zero outside the sphere of radius rmax.

    Parameters
    ----------
    law : str
        'none', 'hubble', 'constant_radial', 'power_law',
        'rotating_solid_body', or 'rotating_galaxy_halo'.
    rmax : float
        Sphere radius (code units).
    v_exp : float
        Characteristic speed [km/s] for radial laws.  Negative = inflow.
    velocity_alpha : float
        Exponent for 'power_law'  (v_r = v_exp * (r/rmax)^velocity_alpha).
    vrot : float
        Maximum rotation speed [km/s] for rotation laws.
    rinner : float or None
        Inner solid-body radius for 'rotating_galaxy_halo' (code units).
    """
    if law == 'none':
        return None

    # ── radial laws ────────────────────────────────────────────────────────
    if law in ('hubble', 'constant_radial', 'power_law', 'linear_decelerate'):
        def fn(x, y, z):
            x_ = np.asarray(x, float); y_ = np.asarray(y, float)
            z_ = np.asarray(z, float)
            r  = np.sqrt(x_**2 + y_**2 + z_**2)
            inside   = (r <= rmax) & (r >= rmin)
            safe_r   = np.where(r > 0.0, r, 1.0)

            if law == 'hubble':
                scale = np.where(inside, v_exp / rmax, 0.0)
                return scale * x_, scale * y_, scale * z_

            elif law == 'constant_radial':
                scale = np.where(inside & (r > 0.0), v_exp / safe_r, 0.0)
                return scale * x_, scale * y_, scale * z_

            elif law == 'linear_decelerate':
                # v_r = v_exp * (rmax - r) / (rmax - rmin)
                denom = rmax - rmin if rmax > rmin else 1.0
                vr = np.where(inside & (r > 0.0),
                              v_exp * (rmax - safe_r) / denom, 0.0)
                scale = np.where(r > 0.0, vr / safe_r, 0.0)
                return scale * x_, scale * y_, scale * z_

            else:  # power_law
                vr    = np.where(inside & (r > 0.0),
                                 v_exp * (safe_r / rmax) ** velocity_alpha, 0.0)
                scale = np.where(r > 0.0, vr / safe_r, 0.0)
                return scale * x_, scale * y_, scale * z_
        return fn

    # ── rotation laws ──────────────────────────────────────────────────────
    if law == 'rotating_solid_body':
        # Garavito-Camargo et al. (2014): solid-body rotation about z-axis.
        # vx = -Vrot * y/rmax,  vy = Vrot * x/rmax,  vz = 0
        scale = vrot / rmax
        def fn(x, y, z):
            x_ = np.asarray(x, float); y_ = np.asarray(y, float)
            z_ = np.asarray(z, float)
            r  = np.sqrt(x_**2 + y_**2 + z_**2)
            inside = r <= rmax
            return (np.where(inside, -scale * y_, 0.0),
                    np.where(inside,  scale * x_, 0.0),
                    np.zeros_like(x_))
        return fn

    if law == 'rotating_galaxy_halo':
        # Kim et al.: flat rotation curve above rinner (rotation about z-axis).
        # rxy < rinner : solid-body  vx=-Vrot*y/rinner, vy=Vrot*x/rinner
        # rxy >= rinner: flat curve  vx=-Vrot*y/rxy,    vy=Vrot*x/rxy
        # vz = 0;  velocity zero outside sphere
        if rinner is None:
            raise ValueError('--rinner is required for rotating_galaxy_halo')
        if rinner <= 0 or rinner >= rmax:
            raise ValueError(
                f'rinner={rinner} must satisfy 0 < rinner < rmax={rmax}')
        scale_in = vrot / rinner
        def fn(x, y, z):
            x_  = np.asarray(x, float); y_ = np.asarray(y, float)
            z_  = np.asarray(z, float)
            r   = np.sqrt(x_**2 + y_**2 + z_**2)
            rxy = np.sqrt(x_**2 + y_**2)
            safe_rxy = np.where(rxy > 0, rxy, 1.0)
            inner  = rxy < rinner
            inside = r   <= rmax
            vx_ = np.where(inner, -scale_in * y_, -vrot * y_ / safe_rxy)
            vy_ = np.where(inner,  scale_in * x_,  vrot * x_ / safe_rxy)
            return (np.where(inside, vx_, 0.0),
                    np.where(inside, vy_, 0.0),
                    np.zeros_like(x_))
        return fn

    raise ValueError(f'Unknown velocity law: {law!r}')


# ---------------------------------------------------------------------------
# Vectorized fast path (no per-cell Cell objects)
# ---------------------------------------------------------------------------
def _shell_radii_array(rmax, level_min, level_max, spacing, ratio, custom_radii):
    """Outer-boundary radius of each refinement shell (outermost first)."""
    n = level_max - level_min + 1
    if custom_radii is not None:
        return np.asarray(custom_radii, dtype=float)
    if spacing == 'linear':
        return np.array([rmax * (n - i) / n for i in range(n)], dtype=float)
    if spacing == 'log':
        return np.array([rmax * ratio ** i for i in range(n)], dtype=float)
    raise ValueError(f"spacing must be 'linear' or 'log', got {spacing!r}")


def _target_level(d, radii_inc, level_min):
    """Vectorized AMRGrid._target(): the finest shell level whose ball encloses
    distance ``d``.  ``radii_inc`` is the shell radii sorted ascending."""
    n_gt = radii_inc.size - np.searchsorted(radii_inc, d, side='right')
    return np.where(n_gt > 0, level_min + n_gt - 1, 0)


def radial_leaves_fast(boxlen, rmax, level_min, level_max,
                       spacing='log', ratio=0.5, custom_radii=None,
                       refine_boundary=False, boundary_level_max=None):
    """Generate the leaf ``(cx, cy, cz, level)`` arrays of a centered
    radial-shell sphere octree by vectorized breadth-first refinement.

    This reproduces ``AMRGrid.refine_sphere_radial`` **bit-for-bit** (verified
    against the object tree) without ever instantiating a ``Cell`` object, which
    makes very deep grids (level >= 8, millions of leaves) feasible in seconds
    and a few hundred MB instead of ~minutes and several GB.

    Centered convention: sphere center = box center = origin.  Refinement
    criterion, applied in a single pass and re-evaluated at every level:

    * radial: split while ``max(target(d_center), target(d_closest)) > level``
      where ``d_center`` is the center-to-origin distance and ``d_closest`` is
      the closest distance from the (axis-aligned) cell to the origin;
    * boundary (optional): cells that straddle the sphere surface — i.e. they
      intersect the sphere (``d_closest <= rmax``) but are not fully inside it
      (``d_farthest > rmax``) — are forced up to ``boundary_level_max``.
    """
    half0 = boxlen / 2.0
    radii = _shell_radii_array(rmax, level_min, level_max, spacing, ratio, custom_radii)
    radii_inc = np.sort(radii)
    blmax = (boundary_level_max if (refine_boundary and boundary_level_max is not None)
             else level_max)
    cx = np.zeros(1); cy = np.zeros(1); cz = np.zeros(1)
    level = 0
    chunks = []
    while cx.size:
        h = half0 / (2 ** level)
        d_cen = np.sqrt(cx * cx + cy * cy + cz * cz)
        ax = np.maximum(np.abs(cx) - h, 0.0)
        ay = np.maximum(np.abs(cy) - h, 0.0)
        az = np.maximum(np.abs(cz) - h, 0.0)
        d_cls = np.sqrt(ax * ax + ay * ay + az * az)
        t = np.maximum(_target_level(d_cen, radii_inc, level_min),
                       _target_level(d_cls, radii_inc, level_min))
        if refine_boundary:
            d_far = np.sqrt((np.abs(cx) + h) ** 2 +
                            (np.abs(cy) + h) ** 2 +
                            (np.abs(cz) + h) ** 2)
            straddle = (d_cls <= rmax) & (d_far > rmax)
            t = np.where(straddle, np.maximum(t, blmax), t)
        refine = level < t
        leaf = ~refine
        if leaf.any():
            chunks.append((cx[leaf], cy[leaf], cz[leaf],
                           np.full(int(leaf.sum()), level, dtype=np.int32)))
        if not refine.any():
            break
        pcx = cx[refine]; pcy = cy[refine]; pcz = cz[refine]
        q = h / 2.0
        s = np.array([-q, q])
        ox, oy, oz = (a.ravel() for a in np.meshgrid(s, s, s, indexing='ij'))
        cx = (pcx[:, None] + ox).ravel()
        cy = (pcy[:, None] + oy).ravel()
        cz = (pcz[:, None] + oz).ravel()
        level += 1
    return (np.concatenate([c[0] for c in chunks]),
            np.concatenate([c[1] for c in chunks]),
            np.concatenate([c[2] for c in chunks]),
            np.concatenate([c[3] for c in chunks]))


class FastAMRGrid(AMRGrid):
    """Array-backed :class:`AMRGrid`: leaf columns are held as numpy arrays
    instead of a ``Cell`` octree.  Produced by the vectorized generators for
    fast writing of very large grids.  It reuses the inherited ``write()``
    machinery (HDF5/FITS) and supports ``info()`` and ``len(grid.leaves())``;
    per-cell iteration and ``slice_plot()`` are intentionally unavailable —
    read the written file back with ``AMRGrid.read()`` for analysis."""

    def __init__(self, boxlen, cx, cy, cz, level, dens, T, vx, vy, vz,
                 origin=None, opt=None, metadata=None):
        self.boxlen = float(boxlen)
        h = self.boxlen * 0.5
        self.origin = (-h, -h, -h) if origin is None else tuple(float(v) for v in origin)
        self.metadata = dict(metadata) if metadata else {}
        n = len(cx)
        self._n = n
        dt = [('x', 'f8'), ('y', 'f8'), ('z', 'f8'), ('level', 'i4'),
              ('dens', 'f8'), ('T', 'f8'), ('vx', 'f8'), ('vy', 'f8'), ('vz', 'f8')]
        opt = opt or {}
        opt_present = [a for a in _OPTIONAL_COLUMNS if opt.get(a) is not None]
        for a in opt_present:
            dt.append((a, 'f8'))
        data = np.empty(n, dtype=dt)
        data['x'] = cx; data['y'] = cy; data['z'] = cz
        data['level'] = level
        data['dens'] = dens; data['T'] = T
        data['vx'] = vx; data['vy'] = vy; data['vz'] = vz
        for a in opt_present:
            data[a] = opt[a]
        self._table = data

    # -- overrides so the inherited write() works without a Cell tree --
    def _leaf_table_array(self):
        return self._table

    def leaves(self):
        # Only len() is used by callers (leaf-count reporting); per-cell
        # iteration is intentionally unsupported on the array-backed grid.
        return range(self._n)

    def info(self):
        d = self._table
        levels = d['level']
        hist = {int(l): int((levels == l).sum()) for l in np.unique(levels)}
        lines = [f'AMRGrid  boxlen={self.boxlen}  nleaf={self._n}',
                 f'  Level distribution: {dict(sorted(hist.items()))}']
        dens = d['dens']
        if dens.max() > 0:
            nz = dens[dens > 0]
            lines.append(f'  dens [cm^-3]: min={nz.min():.3e}  '
                         f'max={dens.max():.3e}  mean={dens.mean():.3e}')
        lines.append(f'  T  [K]:     min={d["T"].min():.3e}  max={d["T"].max():.3e}')
        return '\n'.join(lines)


# ---------------------------------------------------------------------------
# Grid builder
# ---------------------------------------------------------------------------

def build_radial_grid(boxlen, rmax, level_min, level_max, spacing, ratio,
                      density, n0, temperature, sigma=None, r_scale=None,
                      density_alpha=0.0, rmin=0.0,
                      custom_radii=None, refine_boundary=False,
                      boundary_level_max=None,
                      velocity='none', v_exp=0.0, velocity_alpha=1.0,
                      vrot=0.0, rinner=None,
                      cone_opening=0.0, engine='fast'):
    """
    Build and return an AMRGrid with radial shell refinement.

    Returns
    -------
    grid : AMRGrid
    radii : list of float
        Shell boundary radii used (outermost first).
    """
    dens_fn = make_density_fn(density, n0, rmax, sigma=sigma, r_scale=r_scale,
                              density_alpha=density_alpha, rmin=rmin,
                              cone_opening=cone_opening)
    vel_fn  = make_velocity_fn(velocity, rmax,
                               v_exp=v_exp, velocity_alpha=velocity_alpha,
                               vrot=vrot, rinner=rinner, rmin=rmin)

    has_cone = 0.0 < cone_opening < 90.0

    # Build cone region check for restricting radial refinement to the bicone.
    # The bicone interior is |z|/r >= cos(cone_opening), equivalently
    # z^2 * sin^2 >= cos^2 * (x^2 + y^2).  A cell intersects the bicone iff
    # max(|z|)^2 * sin^2 >= cos^2 * min(x^2+y^2) over the cell.
    cone_check = None
    if has_cone:
        cos_cone = np.cos(np.radians(cone_opening))
        sin2 = 1.0 - cos_cone**2
        cos2 = cos_cone**2
        def cone_check(c):
            """True if cell intersects the bicone (exact box-cone test)."""
            z_max = max(abs(c.zmin), abs(c.zmax))
            xm2 = 0.0 if c.xmin <= 0 <= c.xmax else min(c.xmin**2, c.xmax**2)
            ym2 = 0.0 if c.ymin <= 0 <= c.ymax else min(c.ymin**2, c.ymax**2)
            return z_max * z_max * sin2 >= cos2 * (xm2 + ym2)

    # -----------------------------------------------------------------
    # Fast path: vectorized leaf generation (no Cell-object octree).
    # Reproduces the object tree bit-for-bit for the plain radial sphere.
    # Cone geometry (region_check) is not vectorized -> fall back to tree.
    # -----------------------------------------------------------------
    if engine == 'fast' and not has_cone:
        cx, cy, cz, lev = radial_leaves_fast(
            boxlen, rmax, level_min, level_max,
            spacing=spacing, ratio=ratio, custom_radii=custom_radii,
            refine_boundary=refine_boundary, boundary_level_max=boundary_level_max)
        ones = np.ones_like(cx)
        dens = np.asarray(dens_fn(cx, cy, cz), dtype=float) * ones
        if vel_fn is not None:
            vx, vy, vz = vel_fn(cx, cy, cz)
            vx = np.asarray(vx, dtype=float) * ones
            vy = np.asarray(vy, dtype=float) * ones
            vz = np.asarray(vz, dtype=float) * ones
        else:
            vx = np.zeros_like(cx); vy = np.zeros_like(cx); vz = np.zeros_like(cx)
        zero = dens == 0.0          # zero velocity where density is zero
        vx = np.where(zero, 0.0, vx)
        vy = np.where(zero, 0.0, vy)
        vz = np.where(zero, 0.0, vz)
        T = np.full(cx.shape, float(temperature))
        grid = FastAMRGrid(boxlen, cx, cy, cz, lev, dens, T, vx, vy, vz)
        blmax_used = boundary_level_max if boundary_level_max is not None else level_max
        grid.metadata.update({
            'LMIN':   (level_min,       'Outermost shell AMR level'),
            'LMAX':   (level_max,       'Innermost shell AMR level'),
            'BLMAX':  (blmax_used,      'Boundary refinement max level'),
            'REFBND': (refine_boundary, 'Sphere boundary refinement applied'),
        })
        radii = list(_shell_radii_array(rmax, level_min, level_max,
                                        spacing, ratio, custom_radii))
        return grid, radii

    grid = AMRGrid(boxlen)
    grid.refine_sphere_radial(
        0.0, 0.0, 0.0, rmax,
        level_min          = level_min,
        level_max          = level_max,
        spacing            = spacing,
        ratio              = ratio,
        custom_radii       = custom_radii,
        refine_boundary    = refine_boundary and not has_cone,
        boundary_level_max = boundary_level_max,
        region_check       = cone_check,
    )

    # Bicone geometry: combined boundary refinement for the density region
    # (sphere ∩ bicone).  This replaces separate sphere-surface and
    # cone-surface passes, avoiding wasted refinement in density=0 regions
    # outside the bicone.
    if refine_boundary and has_cone:
        cos_cone = np.cos(np.radians(cone_opening))
        blmax = boundary_level_max if boundary_level_max is not None else level_max
        r2 = rmax * rmax

        def _region_inside(c):
            """All 8 corners inside sphere AND inside bicone."""
            corners = c.corners()
            d2 = corners[:,0]**2 + corners[:,1]**2 + corners[:,2]**2
            if not bool(np.all(d2 <= r2)):
                return False
            r = np.sqrt(d2)
            r = np.where(r > 0, r, 1.0)
            cos_th = np.abs(corners[:,2]) / r
            return bool(np.all(cos_th >= cos_cone))

        def _region_intersects(c):
            """Cell overlaps the density region (sphere ∩ bicone)."""
            # Quick reject: must overlap sphere (exact sphere-box test)
            dx = max(c.xmin, 0.0, -c.xmax)
            dy = max(c.ymin, 0.0, -c.ymax)
            dz = max(c.zmin, 0.0, -c.zmax)
            if dx*dx + dy*dy + dz*dz > r2:
                return False
            # At least one corner must be in the density region
            corners = c.corners()
            d2 = corners[:,0]**2 + corners[:,1]**2 + corners[:,2]**2
            in_sphere = d2 <= r2
            r = np.sqrt(d2)
            r_safe = np.where(r > 0, r, 1.0)
            cos_th = np.abs(corners[:,2]) / r_safe
            in_cone = cos_th >= cos_cone
            return bool(np.any(in_sphere & in_cone))

        grid.refine_geometry(_region_inside, _region_intersects, blmax,
                             criterion_inside=None)

    grid.set_density(dens_fn)
    grid.set_temperature(temperature)
    if vel_fn is not None:
        grid.set_velocity(vel_fn)

    # Zero velocity where density is zero
    for lf in grid.leaves():
        if lf.dens == 0.0:
            lf.vx = 0.0; lf.vy = 0.0; lf.vz = 0.0

    # Compute the shell radii actually used (for reporting)
    n_shells = level_max - level_min + 1
    if custom_radii is not None:
        radii = list(custom_radii)
    elif spacing == 'linear':
        radii = [rmax * (n_shells - i) / n_shells for i in range(n_shells)]
    else:
        radii = [rmax * (ratio ** i) for i in range(n_shells)]

    return grid, radii


def print_shell_table(radii, level_min, level_max):
    n = len(radii)
    print()
    print('  Shell structure (outermost → center):')
    print(f'  {"Outer radius":>14}  {"Inner radius":>14}  {"Level":>5}')
    print('  ' + '-' * 38)
    for i, (r_out, lv) in enumerate(zip(radii, range(level_min, level_max + 1))):
        r_in = radii[i + 1] if i + 1 < n else 0.0
        print(f'  {r_out:14.4f}  {r_in:14.4f}  {lv:5d}')
    print()


# ---------------------------------------------------------------------------
# Filename helpers
# ---------------------------------------------------------------------------

def _suffix_filename(path, suffix):
    """Insert ``suffix`` before the recognized AMR-file extension.

    ``foo.h5`` + ``_bar`` → ``foo_bar.h5``; same for ``.hdf5``, ``.fits.gz``,
    ``.fits``, ``.dat``, ``.txt``.  Unknown extensions: append after the path.
    """
    lower = path.lower()
    for ext in ('.fits.gz', '.hdf5', '.fits', '.h5', '.dat', '.txt'):
        if lower.endswith(ext):
            return path[: -len(ext)] + suffix + path[-len(ext):]
    return path + suffix


# ---------------------------------------------------------------------------
# LaRT input snippet
# ---------------------------------------------------------------------------

def print_lart_hint(outfile, boxlen, temperature, taumax=1e4, nphotons=1e6,
                    velocity='none', v_exp=0.0, velocity_alpha=1.0,
                    vrot=0.0, rinner=None):
    print()
    print('─' * 60)
    print('Suggested LaRT input snippet:')
    print('─' * 60)
    print(' &parameters')
    print(f"  par%use_amr_grid    = .true.")
    print(f"  par%amr_type        = 'generic'")
    print(f"  par%amr_file        = '{outfile}'")
    print(f"  par%distance_unit   = 1.0")
    print(f"  par%no_photons      = {nphotons:.0e}")
    print(f"  par%taumax          = {taumax:.1e}")
    print(f"  par%temperature     = {temperature:.4e}")
    print(f"  par%spectral_type   = 'voigt'")
    print(f"  par%source_geometry = 'point'")
    print(f"  par%xs_point        = 0.0")
    print(f"  par%ys_point        = 0.0")
    print(f"  par%zs_point        = 0.0")
    print(f"  par%nxfreq  = 121")
    print(f"  par%nxim    = 100")
    print(f"  par%nyim    = 100")
    if velocity != 'none':
        if velocity in ('hubble', 'constant_radial'):
            vel_info = f"V_exp={v_exp} km/s"
        elif velocity == 'power_law':
            vel_info = f"V_exp={v_exp} km/s, power={velocity_alpha}"
        elif velocity == 'rotating_solid_body':
            vel_info = f"Vrot={vrot} km/s"
        elif velocity == 'rotating_galaxy_halo':
            vel_info = f"Vrot={vrot} km/s, rinner={rinner}"
        else:
            vel_info = ''
        print(f"  ! velocity '{velocity}' ({vel_info}) "
              "embedded in grid file — no par%velocity_type needed")
    print(' /')
    print('─' * 60)


# ---------------------------------------------------------------------------
# Slice plot helpers
# ---------------------------------------------------------------------------

def _slice_velocity_data(grid, slice_val=0.0):
    """
    Collect leaf cells intersecting the z=slice_val plane and return
    arrays of cell centers, half-sizes, in-plane velocity (vx,vy),
    out-of-plane velocity (vz), and speed magnitude |v|.
    """
    leaves = []
    grid._collect_slice_leaves(grid.root, 'cz', slice_val, leaves)
    if not leaves:
        return None
    cx   = np.array([lf.cx  for lf in leaves])
    cy   = np.array([lf.cy  for lf in leaves])
    h    = np.array([lf.h   for lf in leaves])
    vx   = np.array([lf.vx  for lf in leaves])
    vy   = np.array([lf.vy  for lf in leaves])
    vz   = np.array([lf.vz  for lf in leaves])
    vmag = np.sqrt(vx**2 + vy**2 + vz**2)
    return dict(cx=cx, cy=cy, h=h, vx=vx, vy=vy, vz=vz, vmag=vmag)


def _add_quiver(ax, vdata, grid, color='white', alpha=0.55, arrow_frac=0.06):
    """
    Overlay quiver arrows on *ax* for cells with non-zero velocity.
    Arrow length is proportional to speed; max arrow = arrow_frac * boxlen.
    """
    mask = vdata['vmag'] > 0
    if not mask.any():
        return
    cx, cy = vdata['cx'][mask], vdata['cy'][mask]
    vx, vy = vdata['vx'][mask], vdata['vy'][mask]
    vmag   = vdata['vmag'][mask]
    vmax   = vmag.max()
    if vmax == 0:
        return
    scale = vmax / (arrow_frac * grid.boxlen)
    ax.quiver(cx, cy, vx, vy,
              scale=scale, scale_units='xy',
              pivot='mid',
              color=color, alpha=alpha,
              width=0.002, headwidth=4, headlength=4)


def _plot_vmag_panel(ax, vdata, grid, slice_val):
    """
    Draw a velocity-magnitude panel: PolyCollection colored by |v| [km/s]
    plus quiver arrows showing in-plane direction.
    """
    from matplotlib.collections import PolyCollection
    import matplotlib.pyplot as plt

    cx, cy = vdata['cx'], vdata['cy']
    h, vmag = vdata['h'], vdata['vmag']

    # Build axis-aligned rectangle polygons (n, 4, 2)
    verts = np.stack([
        np.column_stack([cx - h, cy - h]),
        np.column_stack([cx + h, cy - h]),
        np.column_stack([cx + h, cy + h]),
        np.column_stack([cx - h, cy + h]),
    ], axis=1)

    # Only color cells with v > 0; background shows vacuum / zero-v region
    ax.set_facecolor('0.15')
    keep = vmag > 0
    if keep.any():
        pc = PolyCollection(verts[keep], array=vmag[keep],
                            cmap='plasma', zorder=1)
        ax.add_collection(pc)
        cb = plt.colorbar(pc, ax=ax, pad=0.02)
        cb.set_label('|v|  [km/s]', fontsize=9)

    # AMR cell boundaries (all cells in slice)
    bc = PolyCollection(verts, facecolor='none',
                        edgecolor='white', linewidth=0.2, alpha=0.4, zorder=2)
    ax.add_collection(bc)

    # Quiver arrows
    _add_quiver(ax, vdata, grid, color='white', alpha=0.75)

    ax.set_xlim(grid.origin[0], grid.origin[0] + grid.boxlen)
    ax.set_ylim(grid.origin[1], grid.origin[1] + grid.boxlen)
    ax.set_aspect('equal')
    ax.set_xlabel('x  [code unit]', fontsize=9)
    ax.set_ylabel('y  [code unit]', fontsize=9)


# ---------------------------------------------------------------------------
# Slice plot
# ---------------------------------------------------------------------------

def make_slice_plot(grids, labels, outfile, log=False, show_velocity=False):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib not available — skipping plot')
        return

    ncols  = len(grids)
    nrows  = 2 if show_velocity else 1
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(6.2 * ncols, 5.5 * nrows),
                             squeeze=False)

    slice_val = 0.0

    for col, (grid, label) in enumerate(zip(grids, labels)):
        ax_dens = axes[0, col]

        # ── density panel ────────────────────────────────────────────────
        grid.slice_plot('z', slice_val, 'dens', ax=ax_dens,
                        log=log, cmap='viridis',
                        show_leaf_boundaries=True,
                        boundary_color='white',
                        boundary_lw=0.25,
                        boundary_alpha=0.5)
        ax_dens.set_title(label, fontsize=10)

        if show_velocity:
            vdata = _slice_velocity_data(grid, slice_val)
            if vdata is not None and vdata['vmag'].max() > 0:
                # quiver overlay on density panel
                _add_quiver(ax_dens, vdata, grid,
                            color='white', alpha=0.55)

                # ── velocity magnitude panel ──────────────────────────────
                ax_vel = axes[1, col]
                _plot_vmag_panel(ax_vel, vdata, grid, slice_val)
                ax_vel.set_title(
                    f'|v|  [km/s]  (arrows: in-plane  vx, vy)', fontsize=10)

    dens_label = ('log$_{10}$ nH  [cm$^{-3}$]' if log else 'nH  [cm$^{-3}$]')
    suptitle = f'Radial shell refinement — row 1: {dens_label}'
    if show_velocity:
        suptitle += '  (arrows: velocity)    row 2: |v| [km/s]'
    fig.suptitle(suptitle, fontsize=10)
    fig.tight_layout()

    png = outfile.replace('.dat', '').replace('.fits.gz', '').replace('.fits', '').replace('.hdf5', '').replace('.h5', '')
    png += '_radial_slice.png'
    plt.savefig(png, dpi=140, bbox_inches='tight')
    print(f'Slice plot saved to {png}')
    plt.show()


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    p.add_argument('--boxlen', type=float, default=2.0,
                   help='Box side length [code unit] (default: 2)')
    p.add_argument('--rmax', type=float, default=None,
                   help='Sphere radius [code unit] (default: boxlen/2)')

    p.add_argument('--density',
                   choices=['uniform', 'gaussian', 'exponential', 'power_law'],
                   default='uniform',
                   help='Density profile (default: uniform)')
    p.add_argument('--density_alpha', type=float, default=0.0,
                   help='[power_law] exponent: n(r) = n0 * (rmax/r)^alpha '
                        '(default: 0; 2 = isothermal)')
    p.add_argument('--rmin', type=float, default=0.0,
                   help='Inner shell radius: density and velocity = 0 for r < rmin '
                        '(default: 0)')
    p.add_argument('--n0', type=float, default=1.0,
                   help='Peak density [cm^-3] (default: 1.0)')
    p.add_argument('--sigma', type=float, default=None,
                   help='Gaussian scale radius (default: 0.4 * rmax)')
    p.add_argument('--r_scale', type=float, default=None,
                   help='Exponential scale radius (default: 0.3 * rmax)')
    p.add_argument('--temperature', type=float, default=1e4,
                   help='Gas temperature [K] (default: 1e4)')

    p.add_argument('--level_min', type=int, default=1,
                   help='Outermost shell refinement level (default: 1)')
    p.add_argument('--level_max', type=int, default=5,
                   help='Innermost (center) refinement level (default: 5)')
    p.add_argument('--spacing', choices=['linear', 'log'], default='linear',
                   help='Shell spacing mode (default: linear)')
    p.add_argument('--ratio', type=float, default=0.5,
                   help='[log] radius reduction factor per level (default: 0.5)')

    p.add_argument('--refine_boundary', action='store_true',
                   help='Also force sphere-surface cells to boundary_level_max '
                        '(sharp geometric boundary on top of radial shells)')
    p.add_argument('--boundary_level_max', type=int, default=None,
                   help='Max AMR level for sphere-surface cells when '
                        '--refine_boundary is set (default: same as --level_max)')
    p.add_argument('--compare', action='store_true',
                   help='Build without and with refine_boundary side-by-side '
                        '(overrides --spacing to show both)')

    vg = p.add_argument_group('velocity law')
    vg.add_argument('--velocity',
                    choices=['none', 'hubble', 'constant_radial', 'power_law',
                             'linear_decelerate',
                             'rotating_solid_body', 'rotating_galaxy_halo'],
                    default='none',
                    help='Velocity law (default: none)\n'
                         '  hubble               : v_r = V_exp * r/rmax\n'
                         '  constant_radial      : v_r = V_exp\n'
                         '  power_law            : v_r = V_exp * (r/rmax)^velocity_alpha\n'
                         '  linear_decelerate    : v_r = V_exp * (rmax-r)/(rmax-rmin)\n'
                         '  rotating_solid_body  : solid-body rotation about z (Garavito-Camargo+2014)\n'
                         '  rotating_galaxy_halo : flat rotation above rinner (Kim et al.)')
    vg.add_argument('--vexp', type=float, default=0.0,
                    help='[hubble/constant_radial/power_law] '
                         'characteristic speed [km/s] (default: 0). Negative = inflow.')
    vg.add_argument('--velocity_alpha', type=float, default=1.0,
                    help='[power_law] exponent (default: 1.0)')
    vg.add_argument('--cone_opening', type=float, default=0.0,
                    help='Bicone half-opening angle [deg]; 0 = full sphere (default: 0)')
    vg.add_argument('--vrot', type=float, default=0.0,
                    help='[rotating_*] maximum rotation speed [km/s] (default: 0)')
    vg.add_argument('--rinner', type=float, default=None,
                    help='[rotating_galaxy_halo] inner solid-body radius in code units')

    p.add_argument('-o', '--output', default='sphere_radial.h5',
                   help='Output AMR data file (default: sphere_radial.h5)')
    p.add_argument('--taumax', type=float, default=1e4,
                   help='taumax for the LaRT input snippet (default: 1e4)')
    p.add_argument('--plot', action='store_true',
                   help='Save a z=0 slice plot with AMR cell boundaries')
    p.add_argument('--plot_log', action='store_true',
                   help='Use log scale for density in slice plot')

    return p.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    rmax = args.rmax if args.rmax is not None else args.boxlen * 0.5

    common = dict(
        boxlen             = args.boxlen,
        rmax               = rmax,
        level_min          = args.level_min,
        level_max          = args.level_max,
        spacing            = args.spacing,
        ratio              = args.ratio,
        density            = args.density,
        n0                 = args.n0,
        temperature        = args.temperature,
        sigma              = args.sigma,
        r_scale            = args.r_scale,
        density_alpha      = args.density_alpha,
        rmin               = args.rmin,
        boundary_level_max = args.boundary_level_max,
        velocity           = args.velocity,
        v_exp              = args.vexp,
        velocity_alpha            = args.velocity_alpha,
        vrot               = args.vrot,
        rinner             = args.rinner,
        cone_opening       = args.cone_opening,
    )

    # The vectorized fast engine builds no Cell tree, so slice plots (which
    # need per-cell traversal) require the object-tree engine.
    _engine = 'tree' if args.plot else 'fast'

    if args.compare:
        # ── radial only vs. radial + boundary ────────────────────────────
        print('=' * 60)
        print('Mode: radial shells only')
        print('=' * 60)
        grid_base, radii_base = build_radial_grid(
            **common, refine_boundary=False, engine=_engine)
        print(grid_base.info())
        print_shell_table(radii_base, args.level_min, args.level_max)

        out_base = _suffix_filename(args.output, '_no_boundary')
        grid_base.write(out_base)

        print()
        print('=' * 60)
        print('Mode: radial shells + boundary refinement')
        print('=' * 60)
        grid_bnd, radii_bnd = build_radial_grid(
            **common, refine_boundary=True, engine=_engine)
        print(grid_bnd.info())
        print_shell_table(radii_bnd, args.level_min, args.level_max)

        out_bnd = _suffix_filename(args.output, '_with_boundary')
        grid_bnd.write(out_bnd)

        if args.plot:
            make_slice_plot(
                [grid_base, grid_bnd],
                [f'radial only — {grid_base.level_counts()}',
                 f'+ boundary  — {grid_bnd.level_counts()}'],
                args.output,
                log=args.plot_log,
                show_velocity=(args.velocity != 'none'),
            )
    else:
        # ── build a single grid ───────────────────────────────────────────
        print(f'Building AMR grid  boxlen={args.boxlen}  rmax={rmax}')
        print(f'  density  : {args.density}  n0={args.n0:.2e} cm^-3')
        print(f'  T        : {args.temperature:.2e} K')
        if args.velocity != 'none':
            if args.velocity in ('rotating_solid_body', 'rotating_galaxy_halo'):
                vel_desc = f'{args.velocity}  Vrot={args.vrot} km/s'
                if args.velocity == 'rotating_galaxy_halo':
                    vel_desc += f'  rinner={args.rinner}'
            else:
                vel_desc = f'{args.velocity}  V_exp={args.vexp} km/s'
                if args.velocity == 'power_law':
                    vel_desc += f'  power={args.velocity_alpha}'
            print(f'  velocity : {vel_desc}')
        blmax = args.boundary_level_max if args.boundary_level_max is not None else args.level_max
        print(f'  AMR      : level_min={args.level_min}  level_max={args.level_max}'
              f'  spacing={args.spacing}'
              + (f'  ratio={args.ratio}' if args.spacing == 'log' else '')
              + (f'  refine_boundary=True  boundary_level_max={blmax}'
                 if args.refine_boundary else ''))

        grid, radii = build_radial_grid(
            **common, refine_boundary=args.refine_boundary, engine=_engine)

        print(grid.info())
        print_shell_table(radii, args.level_min, args.level_max)

        grid.write(args.output)

        if args.plot:
            label = (f'{args.spacing} + boundary' if args.refine_boundary
                     else args.spacing)
            if args.velocity == 'none':
                vel_tag = ''
            elif args.velocity in ('rotating_solid_body', 'rotating_galaxy_halo'):
                vel_tag = f' + {args.velocity}(Vrot={args.vrot} km/s)'
            else:
                vel_tag = f' + {args.velocity}({args.vexp} km/s)'
            make_slice_plot([grid],
                            [f'{label} spacing — {args.density} density{vel_tag}'],
                            args.output,
                            log=args.plot_log,
                            show_velocity=(args.velocity != 'none'))

        print_lart_hint(args.output, args.boxlen, args.temperature,
                        taumax=args.taumax,
                        velocity=args.velocity, v_exp=args.vexp,
                        velocity_alpha=args.velocity_alpha,
                        vrot=args.vrot, rinner=args.rinner)


if __name__ == '__main__':
    main()
