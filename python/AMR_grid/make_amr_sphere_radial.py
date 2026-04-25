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
  power_law      — power-law radial: v_r = V_exp * (r/rmax)^v_power

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
    python make_amr_sphere_radial.py --velocity power_law --vexp 200 --v_power 2.0

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
from AMR_grid import AMRGrid


# ---------------------------------------------------------------------------
# Density profiles
# ---------------------------------------------------------------------------

def make_density_fn(profile, n0, rmax, sigma=None, r_scale=None):
    if profile == 'uniform':
        def fn(x, y, z):
            r = np.sqrt(np.asarray(x)**2 + np.asarray(y)**2 + np.asarray(z)**2)
            return np.where(r <= rmax, n0, 0.0)
    elif profile == 'gaussian':
        sig = sigma if sigma else rmax * 0.4
        def fn(x, y, z):
            r2 = np.asarray(x)**2 + np.asarray(y)**2 + np.asarray(z)**2
            return np.where(np.sqrt(r2) <= rmax,
                            n0 * np.exp(-r2 / (2.0 * sig**2)), 0.0)
    elif profile == 'exponential':
        rs = r_scale if r_scale else rmax * 0.3
        def fn(x, y, z):
            r = np.sqrt(np.asarray(x)**2 + np.asarray(y)**2 + np.asarray(z)**2)
            return np.where(r <= rmax, n0 * np.exp(-r / rs), 0.0)
    else:
        raise ValueError(f'Unknown density profile: {profile!r}')
    return fn


# ---------------------------------------------------------------------------
# Velocity laws
# ---------------------------------------------------------------------------

def make_velocity_fn(law, rmax, v_exp=0.0, v_power=1.0, vrot=0.0, rinner=None):
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
    v_power : float
        Exponent for 'power_law'  (v_r = v_exp * (r/rmax)^v_power).
    vrot : float
        Maximum rotation speed [km/s] for rotation laws.
    rinner : float or None
        Inner solid-body radius for 'rotating_galaxy_halo' (code units).
    """
    if law == 'none':
        return None

    # ── radial laws ────────────────────────────────────────────────────────
    if law in ('hubble', 'constant_radial', 'power_law'):
        def fn(x, y, z):
            x_ = np.asarray(x, float); y_ = np.asarray(y, float)
            z_ = np.asarray(z, float)
            r  = np.sqrt(x_**2 + y_**2 + z_**2)
            inside   = r <= rmax
            safe_r   = np.where(r > 0.0, r, 1.0)

            if law == 'hubble':
                # v_r = v_exp * r/rmax  →  v_i = v_exp/rmax * x_i
                scale = np.where(inside, v_exp / rmax, 0.0)
                return scale * x_, scale * y_, scale * z_

            elif law == 'constant_radial':
                # v_r = v_exp (unit radial vector)
                scale = np.where(inside & (r > 0.0), v_exp / safe_r, 0.0)
                return scale * x_, scale * y_, scale * z_

            else:  # power_law
                # v_r = v_exp * (r/rmax)^v_power
                vr    = np.where(inside & (r > 0.0),
                                 v_exp * (safe_r / rmax) ** v_power, 0.0)
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
# Grid builder
# ---------------------------------------------------------------------------

def build_radial_grid(boxlen, rmax, level_min, level_max, spacing, ratio,
                      density, n0, temperature, sigma=None, r_scale=None,
                      custom_radii=None, refine_boundary=False,
                      boundary_level_max=None,
                      velocity='none', v_exp=0.0, v_power=1.0,
                      vrot=0.0, rinner=None):
    """
    Build and return an AMRGrid with radial shell refinement.

    Returns
    -------
    grid : AMRGrid
    radii : list of float
        Shell boundary radii used (outermost first).
    """
    dens_fn = make_density_fn(density, n0, rmax, sigma=sigma, r_scale=r_scale)
    vel_fn  = make_velocity_fn(velocity, rmax,
                               v_exp=v_exp, v_power=v_power,
                               vrot=vrot, rinner=rinner)

    grid = AMRGrid(boxlen)
    grid.refine_sphere_radial(
        0.0, 0.0, 0.0, rmax,
        level_min          = level_min,
        level_max          = level_max,
        spacing            = spacing,
        ratio              = ratio,
        custom_radii       = custom_radii,
        refine_boundary    = refine_boundary,
        boundary_level_max = boundary_level_max,
    )

    grid.set_density(dens_fn)
    grid.set_temperature(temperature)
    if vel_fn is not None:
        grid.set_velocity(vel_fn)

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
    print('  Shell structure (outermost → centre):')
    print(f'  {"Outer radius":>14}  {"Inner radius":>14}  {"Level":>5}')
    print('  ' + '-' * 38)
    for i, (r_out, lv) in enumerate(zip(radii, range(level_min, level_max + 1))):
        r_in = radii[i + 1] if i + 1 < n else 0.0
        print(f'  {r_out:14.4f}  {r_in:14.4f}  {lv:5d}')
    print()


# ---------------------------------------------------------------------------
# LaRT input snippet
# ---------------------------------------------------------------------------

def print_lart_hint(outfile, boxlen, temperature, taumax=1e4, nphotons=1e6,
                    velocity='none', v_exp=0.0, v_power=1.0,
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
            vel_info = f"V_exp={v_exp} km/s, power={v_power}"
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
    Draw a velocity-magnitude panel: PolyCollection coloured by |v| [km/s]
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

    # Only colour cells with v > 0; background shows vacuum / zero-v region
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

    png = outfile.replace('.dat', '').replace('.fits.gz', '').replace('.fits', '')
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
                   choices=['uniform', 'gaussian', 'exponential'],
                   default='uniform',
                   help='Density profile (default: uniform)')
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
                   help='Innermost (centre) refinement level (default: 5)')
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
                             'rotating_solid_body', 'rotating_galaxy_halo'],
                    default='none',
                    help='Velocity law (default: none)\n'
                         '  hubble               : v_r = V_exp * r/rmax\n'
                         '  constant_radial      : v_r = V_exp\n'
                         '  power_law            : v_r = V_exp * (r/rmax)^v_power\n'
                         '  rotating_solid_body  : solid-body rotation about z (Garavito-Camargo+2014)\n'
                         '  rotating_galaxy_halo : flat rotation above rinner (Kim et al.)')
    vg.add_argument('--vexp', type=float, default=0.0,
                    help='[hubble/constant_radial/power_law] '
                         'characteristic speed [km/s] (default: 0). Negative = inflow.')
    vg.add_argument('--v_power', type=float, default=1.0,
                    help='[power_law] exponent (default: 1.0)')
    vg.add_argument('--vrot', type=float, default=0.0,
                    help='[rotating_*] maximum rotation speed [km/s] (default: 0)')
    vg.add_argument('--rinner', type=float, default=None,
                    help='[rotating_galaxy_halo] inner solid-body radius in code units')

    p.add_argument('-o', '--output', default='sphere_radial.dat',
                   help='Output AMR data file (default: sphere_radial.dat)')
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
        boundary_level_max = args.boundary_level_max,
        velocity           = args.velocity,
        v_exp              = args.vexp,
        v_power            = args.v_power,
        vrot               = args.vrot,
        rinner             = args.rinner,
    )

    if args.compare:
        # ── radial only vs. radial + boundary ────────────────────────────
        print('=' * 60)
        print('Mode: radial shells only')
        print('=' * 60)
        grid_base, radii_base = build_radial_grid(
            **common, refine_boundary=False)
        print(grid_base.info())
        print_shell_table(radii_base, args.level_min, args.level_max)

        out_base = args.output.replace('.dat', '_no_boundary.dat')
        grid_base.write(out_base)

        print()
        print('=' * 60)
        print('Mode: radial shells + boundary refinement')
        print('=' * 60)
        grid_bnd, radii_bnd = build_radial_grid(
            **common, refine_boundary=True)
        print(grid_bnd.info())
        print_shell_table(radii_bnd, args.level_min, args.level_max)

        out_bnd = args.output.replace('.dat', '_with_boundary.dat')
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
                    vel_desc += f'  power={args.v_power}'
            print(f'  velocity : {vel_desc}')
        blmax = args.boundary_level_max if args.boundary_level_max is not None else args.level_max
        print(f'  AMR      : level_min={args.level_min}  level_max={args.level_max}'
              f'  spacing={args.spacing}'
              + (f'  ratio={args.ratio}' if args.spacing == 'log' else '')
              + (f'  refine_boundary=True  boundary_level_max={blmax}'
                 if args.refine_boundary else ''))

        grid, radii = build_radial_grid(
            **common, refine_boundary=args.refine_boundary)

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
                        v_power=args.v_power,
                        vrot=args.vrot, rinner=args.rinner)


if __name__ == '__main__':
    main()
