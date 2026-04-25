#!/usr/bin/env python3
"""
make_amr_grid.py
===============
Generate a LaRT v2.00 generic AMR data file for a spherical medium
with radial density profiles and various velocity fields.

Supported density profiles (radial, centred at origin):
  uniform       dens = n0                          (constant inside sphere)
  gaussian      dens = n0 * exp(-r^2 / (2*sigma^2))
  exponential   dens = n0 * exp(-r / r_scale)

Supported velocity fields [km/s] (matching grid_mod_car.f90):
  none                zero velocity
  hubble              v ∝ r  (vx=Vexp*x/rmax, vy=..., vz=...; linear expansion)
  constant_radial     |v| = Vexp everywhere (radially outward)
  ssh                 Song, Seon & Hwang (2020) two-zone radial:
                        r < rpeak : v = Vpeak*r/rpeak  (linear core)
                        r >= rpeak: v = Vpeak + DeltaV*(r-rpeak)/(rmax-rpeak)
  rotating_solid_body solid-body rotation about z-axis (Garavito-Camargo+2014):
                        vx = -Vrot*y/rmax,  vy = Vrot*x/rmax,  vz = 0
  rotating_galaxy_halo flat rotation curve above rinner (Kim et al.):
                        rxy = sqrt(x^2+y^2)
                        rxy < rinner: vx = -Vrot*y/rinner, vy = Vrot*x/rinner
                        rxy >= rinner: vx = -Vrot*y/rxy,   vy = Vrot*x/rxy
                        vz = 0

Usage examples
--------------
# Gaussian density + Hubble expansion
python make_amr_grid.py --density gaussian --n0 1e-2 --sigma 20 \\
    --velocity hubble --Vexp 200 \\
    --boxlen 100 --rmax 50 --temperature 1e4 --level_max 6 -o sphere_hubble.dat

# Exponential density + SSH velocity
python make_amr_grid.py --density exponential --n0 1e-3 --r_scale 15 \\
    --velocity ssh --Vpeak 200 --rpeak 20 --DeltaV 100 \\
    --boxlen 100 --rmax 50 --level_max 5 -o sphere_ssh.dat

# Uniform density + rotating solid body
python make_amr_grid.py --density uniform --n0 1e-2 \\
    --velocity rotating_solid_body --Vrot 150 \\
    --boxlen 100 --rmax 50 --level_max 5 -o sphere_rot.dat

# Gaussian density + rotating galaxy halo
python make_amr_grid.py --density gaussian --n0 1e-2 --sigma 25 \\
    --velocity rotating_galaxy_halo --Vrot 200 --rinner 15 \\
    --boxlen 100 --rmax 50 --level_max 6 -o sphere_gal.dat
"""

import argparse
import os
import sys

import numpy as np

# Allow running from any working directory
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from AMR_grid import AMRGrid  # noqa: E402


# ---------------------------------------------------------------------------
# Density profile factory
# ---------------------------------------------------------------------------

def make_density_fn(profile, n0, rmax, sigma=None, r_scale=None):
    """
    Return a vectorisable density function  fn(x, y, z) -> nH [cm^-3].
    Density is zero outside the sphere of radius rmax.
    """
    if profile == 'uniform':
        def fn(x, y, z):
            r = np.sqrt(np.asarray(x, float)**2 +
                        np.asarray(y, float)**2 +
                        np.asarray(z, float)**2)
            return np.where(r <= rmax, n0, 0.0)

    elif profile == 'gaussian':
        if not sigma or sigma <= 0:
            raise ValueError('--sigma must be > 0 for gaussian density')
        def fn(x, y, z):
            r2 = (np.asarray(x, float)**2 +
                  np.asarray(y, float)**2 +
                  np.asarray(z, float)**2)
            return np.where(np.sqrt(r2) <= rmax,
                            n0 * np.exp(-r2 / (2.0 * sigma**2)),
                            0.0)

    elif profile == 'exponential':
        if not r_scale or r_scale <= 0:
            raise ValueError('--r_scale must be > 0 for exponential density')
        def fn(x, y, z):
            r = np.sqrt(np.asarray(x, float)**2 +
                        np.asarray(y, float)**2 +
                        np.asarray(z, float)**2)
            return np.where(r <= rmax,
                            n0 * np.exp(-r / r_scale),
                            0.0)
    else:
        raise ValueError(f'Unknown density profile: {profile!r}')

    return fn


# ---------------------------------------------------------------------------
# Velocity field factory
# ---------------------------------------------------------------------------

def make_velocity_fn(vtype, rmax, Vexp=0.0, Vpeak=0.0, rpeak=None,
                     DeltaV=0.0, Vrot=0.0, rinner=None):
    """
    Return a vectorisable velocity function  fn(x, y, z) -> (vx, vy, vz) [km/s].
    Matches the implementations in grid_mod_car.f90.
    """
    if vtype == 'none':
        def fn(x, y, z):
            z0 = np.zeros_like(np.asarray(x, float))
            return z0, z0.copy(), z0.copy()

    elif vtype == 'hubble':
        # v = Vexp * r / rmax  (linear Hubble-like expansion)
        # vx = Vexp/rmax * x,  vy = Vexp/rmax * y,  vz = Vexp/rmax * z
        scale = Vexp / rmax
        def fn(x, y, z):
            x_ = np.asarray(x, float)
            y_ = np.asarray(y, float)
            z_ = np.asarray(z, float)
            return scale * x_, scale * y_, scale * z_

    elif vtype == 'constant_radial':
        # |v| = Vexp radially outward; v=0 at origin
        eps = rmax * 1e-6
        def fn(x, y, z):
            x_ = np.asarray(x, float)
            y_ = np.asarray(y, float)
            z_ = np.asarray(z, float)
            rr = np.sqrt(x_**2 + y_**2 + z_**2)
            safe = np.where(rr > eps, rr, 1.0)
            mask = rr > eps
            return (np.where(mask, Vexp * x_ / safe, 0.0),
                    np.where(mask, Vexp * y_ / safe, 0.0),
                    np.where(mask, Vexp * z_ / safe, 0.0))

    elif vtype == 'ssh':
        # Song, Seon & Hwang (2020) two-zone radial velocity.
        # Inner zone (r < rpeak): solid Hubble-like, v ∝ r
        #   vx = Vpeak/rpeak * x, ...
        # Outer zone (r >= rpeak): linear ramp from Vpeak to Vpeak+DeltaV
        #   Vscale = Vpeak + DeltaV*(r-rpeak)/(rmax-rpeak)
        #   vx = Vscale * x/r, ...
        if rpeak is None:
            raise ValueError('--rpeak is required for ssh velocity')
        if rpeak <= 0 or rpeak >= rmax:
            raise ValueError(f'--rpeak must satisfy 0 < rpeak < rmax={rmax}')
        scale_inner = Vpeak / rpeak
        denom_outer  = max(rmax - rpeak, 1e-30)
        def fn(x, y, z):
            x_ = np.asarray(x, float)
            y_ = np.asarray(y, float)
            z_ = np.asarray(z, float)
            rr   = np.sqrt(x_**2 + y_**2 + z_**2)
            safe = np.where(rr > 0, rr, 1.0)
            Vscale_out = Vpeak + DeltaV * (rr - rpeak) / denom_outer
            inner = rr < rpeak
            vx = np.where(inner, scale_inner * x_, Vscale_out * x_ / safe)
            vy = np.where(inner, scale_inner * y_, Vscale_out * y_ / safe)
            vz = np.where(inner, scale_inner * z_, Vscale_out * z_ / safe)
            return vx, vy, vz

    elif vtype == 'rotating_solid_body':
        # Garavito-Camargo et al. (2014): solid-body rotation about z-axis.
        # vx = -Vrot * y/rmax,  vy = Vrot * x/rmax,  vz = 0
        scale = Vrot / rmax
        def fn(x, y, z):
            x_ = np.asarray(x, float)
            y_ = np.asarray(y, float)
            return -scale * y_, scale * x_, np.zeros_like(x_)

    elif vtype == 'rotating_galaxy_halo':
        # Kim et al.: flat rotation curve above rinner.
        # rxy = sqrt(x^2+y^2)
        # rxy < rinner: solid-body  vx=-Vrot*y/rinner, vy=Vrot*x/rinner
        # rxy >= rinner: flat curve vx=-Vrot*y/rxy,    vy=Vrot*x/rxy
        # vz = 0
        if rinner is None:
            raise ValueError('--rinner is required for rotating_galaxy_halo')
        if rinner <= 0:
            raise ValueError(f'--rinner must be > 0, got {rinner}')
        scale_in = Vrot / rinner
        def fn(x, y, z):
            x_ = np.asarray(x, float)
            y_ = np.asarray(y, float)
            rxy  = np.sqrt(x_**2 + y_**2)
            safe = np.where(rxy > 0, rxy, 1.0)
            inner = rxy < rinner
            vx = np.where(inner, -scale_in * y_, -Vrot * y_ / safe)
            vy = np.where(inner,  scale_in * x_,  Vrot * x_ / safe)
            return vx, vy, np.zeros_like(x_)

    else:
        raise ValueError(f'Unknown velocity type: {vtype!r}')

    return fn


# ---------------------------------------------------------------------------
# Grid builder
# ---------------------------------------------------------------------------

def build_grid(args):
    rmax = args.rmax if args.rmax is not None else args.boxlen * 0.5

    dens_fn = make_density_fn(args.density, args.n0, rmax,
                               sigma=args.sigma, r_scale=args.r_scale)
    vel_fn  = make_velocity_fn(
        args.velocity, rmax,
        Vexp   = args.Vexp,
        Vpeak  = args.Vpeak,
        rpeak  = args.rpeak,
        DeltaV = args.DeltaV,
        Vrot   = args.Vrot,
        rinner = args.rinner,
    )

    print(f'Building AMR grid  boxlen={args.boxlen}  rmax={rmax}')
    print(f'  density  : {args.density}  n0={args.n0:.2e} cm^-3'
          + (f'  sigma={args.sigma}'    if args.density == 'gaussian'    else '')
          + (f'  r_scale={args.r_scale}' if args.density == 'exponential' else ''))
    print(f'  velocity : {args.velocity}'
          + (_vel_summary(args)))
    print(f'  T        : {args.temperature:.2e} K')
    print(f'  AMR      : level_min={args.level_min}  level_max={args.level_max}')

    grid = AMRGrid(args.boxlen)

    # For uniform density the physics gradient only fires at the sphere
    # boundary (step from n0 → 0); force the boundary cells to level_max
    # so the surface is sharply resolved.
    force_boundary = (args.density == 'uniform') or args.refine_boundary

    grid.refine_sphere_by_physics(
        0.0, 0.0, 0.0, rmax,
        dens_fn        = dens_fn,
        vel_fn         = (vel_fn if args.velocity != 'none' else None),
        dens_threshold = args.dens_threshold,
        vel_threshold  = args.vel_threshold,
        level_max      = args.level_max,
        level_min      = args.level_min,
        refine_boundary= force_boundary,
    )

    grid.set_density(dens_fn)
    grid.set_temperature(args.temperature)
    grid.set_velocity(vel_fn)

    return grid, rmax


def _vel_summary(args):
    v = args.velocity
    if v == 'none':               return ''
    if v in ('hubble', 'constant_radial'):
        return f'  Vexp={args.Vexp} km/s'
    if v == 'ssh':
        return (f'  Vpeak={args.Vpeak} km/s'
                f'  rpeak={args.rpeak}'
                f'  DeltaV={args.DeltaV} km/s')
    if v == 'rotating_solid_body':
        return f'  Vrot={args.Vrot} km/s'
    if v == 'rotating_galaxy_halo':
        return f'  Vrot={args.Vrot} km/s  rinner={args.rinner}'
    return ''


# ---------------------------------------------------------------------------
# LaRT input file hint
# ---------------------------------------------------------------------------

def print_lart_hint(args, rmax, outfile):
    bu = args.distance_unit or ''
    print()
    print('─' * 60)
    print('Suggested LaRT input snippet:')
    print('─' * 60)
    print(' &parameters')
    print(f'  par%use_amr_grid    = .true.')
    print(f'  par%amr_type        = \'generic\'')
    print(f'  par%amr_file        = \'{outfile}\'')
    if bu:
        print(f'  par%distance_unit   = \'{bu}\'')
    print(f'  par%taumax          = 1.0e4    ! tau: box centre -> +z edge')
    print(f'  par%temperature     = {args.temperature:.4e}')
    print(f'  par%source_geometry = \'point\'')
    print(f'  par%xs_point        = 0.0      ! box centre')
    print(f'  par%ys_point        = 0.0')
    print(f'  par%zs_point        = 0.0')
    if args.velocity != 'none':
        print(f'  ! velocity set in AMR data file; par%velocity_type is ignored')
    print(' /')
    print('─' * 60)


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # ── geometry ─────────────────────────────────────────────────────────────
    g = p.add_argument_group('Grid geometry')
    g.add_argument('--boxlen', type=float, default=2.0,
                   help='Box side length in code units (default: 2)')
    g.add_argument('--rmax', type=float, default=None,
                   help='Sphere radius in code units (default: boxlen/2)')
    g.add_argument('--distance_unit', default='',
                   help="Unit of coordinates for LaRT: 'kpc', 'pc', 'au', "
                        "or '' for dimensionless (default: '')")

    # ── density ──────────────────────────────────────────────────────────────
    d = p.add_argument_group('Density profile')
    d.add_argument('--density',
                   choices=['uniform', 'gaussian', 'exponential'],
                   default='gaussian',
                   help='Profile type (default: gaussian)')
    d.add_argument('--n0', type=float, default=1.0,
                   help='Central/peak number density [cm^-3] (default: 1.0)')
    d.add_argument('--sigma', type=float, default=None,
                   help='Gaussian scale radius (code units; required for gaussian)')
    d.add_argument('--r_scale', type=float, default=None,
                   help='Exponential scale radius (code units; required for exponential)')

    # ── temperature ──────────────────────────────────────────────────────────
    p.add_argument('--temperature', type=float, default=1e4,
                   help='Gas temperature [K] (uniform throughout; default: 1e4)')

    # ── velocity ─────────────────────────────────────────────────────────────
    v = p.add_argument_group('Velocity field')
    v.add_argument('--velocity',
                   choices=['none', 'hubble', 'constant_radial', 'ssh',
                            'rotating_solid_body', 'rotating_galaxy_halo'],
                   default='none',
                   help='Velocity field type (default: none)')
    v.add_argument('--Vexp', type=float, default=0.0,
                   help='[hubble / constant_radial]  expansion speed [km/s]')
    v.add_argument('--Vpeak', type=float, default=0.0,
                   help='[ssh]  peak radial speed at rpeak [km/s]')
    v.add_argument('--rpeak', type=float, default=None,
                   help='[ssh]  transition radius (code units)')
    v.add_argument('--DeltaV', type=float, default=0.0,
                   help='[ssh]  additional speed from rpeak to rmax [km/s]; '
                        'v(rmax) = Vpeak + DeltaV')
    v.add_argument('--Vrot', type=float, default=0.0,
                   help='[rotating_*]  maximum rotation speed [km/s]')
    v.add_argument('--rinner', type=float, default=None,
                   help='[rotating_galaxy_halo]  inner solid-body radius (code units)')

    # ── AMR refinement ────────────────────────────────────────────────────────
    r = p.add_argument_group('AMR refinement')
    r.add_argument('--level_max', type=int, default=5,
                   help='Maximum AMR level (default: 5 → 32^3 finest cells)')
    r.add_argument('--level_min', type=int, default=2,
                   help='Minimum AMR level before gradient scan (default: 2)')
    r.add_argument('--dens_threshold', type=float, default=0.1,
                   help='Density gradient threshold [0,1) for refinement (default: 0.1)')
    r.add_argument('--vel_threshold', type=float, default=0.1,
                   help='Velocity gradient threshold [0,1) for refinement (default: 0.1)')
    r.add_argument('--refine_boundary', action='store_true',
                   help='Force sphere surface cells to level_max '
                        '(auto-enabled for uniform density)')

    # ── output ────────────────────────────────────────────────────────────────
    p.add_argument('-o', '--output', default='amr_grid.dat',
                   help='Output file name (default: amr_grid.dat)')
    p.add_argument('--plot', action='store_true',
                   help='Save a z=0 slice plot (density + velocity)')
    p.add_argument('--plot_log', action='store_true',
                   help='Use log scale for density in slice plot')

    return p.parse_args()


# ---------------------------------------------------------------------------
# Slice plot
# ---------------------------------------------------------------------------

def make_slice_plot(grid, args, outfile):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib not available — skipping plot')
        return

    vnames = {'hubble': ('vx', 'vy'),
              'constant_radial': ('vx', 'vy'),
              'ssh': ('vx', 'vy'),
              'rotating_solid_body': ('vx', 'vy'),
              'rotating_galaxy_halo': ('vx', 'vy'),
              'none': None}
    vcomps = vnames.get(args.velocity)

    ncols = 3 if vcomps else 1
    fig, axes = plt.subplots(1, ncols, figsize=(5 * ncols * 1.07, 5))
    if ncols == 1:
        axes = [axes]

    grid.slice_plot('z', 0.0, 'dens', ax=axes[0],
                    log=args.plot_log, cmap='viridis')
    axes[0].set_title('log$_{10}$ nH [cm$^{-3}$]' if args.plot_log
                       else 'nH [cm$^{-3}$]')

    if vcomps:
        for ax, comp in zip(axes[1:], vcomps):
            grid.slice_plot('z', 0.0, comp, ax=ax,
                            log=False, cmap='RdBu_r',
                            show_leaf_boundaries=True,
                            background_color='white')
            ax.set_title(f'{comp} [km/s]')

    plt.suptitle(
        f'{args.density} density + {args.velocity} velocity\n'
        f'boxlen={args.boxlen}  rmax={args.rmax}  '
        f'level_max={args.level_max}',
        fontsize=10
    )
    plt.tight_layout()
    png = outfile.replace('.dat', '').replace('.fits.gz', '').replace('.fits', '')
    png += '_slice.png'
    plt.savefig(png, dpi=120, bbox_inches='tight')
    plt.show()
    print(f'Slice plot saved to {png}')


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    # Defaults that depend on other arguments
    if args.rmax is None:
        args.rmax = args.boxlen * 0.5
    if args.density == 'gaussian' and args.sigma is None:
        args.sigma = args.rmax * 0.4   # sensible default: 40% of rmax
        print(f'Note: --sigma not specified; defaulting to {args.sigma:.4g} '
              f'(= 0.4 * rmax)')
    if args.density == 'exponential' and args.r_scale is None:
        args.r_scale = args.rmax * 0.3
        print(f'Note: --r_scale not specified; defaulting to {args.r_scale:.4g} '
              f'(= 0.3 * rmax)')

    grid, rmax = build_grid(args)
    print(grid.info())

    # velocity stats
    leaves = grid.leaves()
    vx_arr = np.array([lf.vx for lf in leaves])
    vy_arr = np.array([lf.vy for lf in leaves])
    vz_arr = np.array([lf.vz for lf in leaves])
    vmag   = np.sqrt(vx_arr**2 + vy_arr**2 + vz_arr**2)
    if vmag.max() > 0:
        print(f'  |v| [km/s]: min={vmag.min():.2f}  max={vmag.max():.2f}  '
              f'mean={vmag.mean():.2f}')

    grid.write(args.output)

    if args.plot:
        make_slice_plot(grid, args, args.output)

    print_lart_hint(args, rmax, args.output)


if __name__ == '__main__':
    main()
