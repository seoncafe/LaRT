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
# Grid builder
# ---------------------------------------------------------------------------

def build_radial_grid(boxlen, rmax, level_min, level_max, spacing, ratio,
                      density, n0, temperature, sigma=None, r_scale=None,
                      custom_radii=None, refine_boundary=False,
                      boundary_level_max=None):
    """
    Build and return an AMRGrid with radial shell refinement.

    Returns
    -------
    grid : AMRGrid
    radii : list of float
        Shell boundary radii used (outermost first).
    """
    dens_fn = make_density_fn(density, n0, rmax, sigma=sigma, r_scale=r_scale)

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

def print_lart_hint(outfile, boxlen, temperature, taumax=1e4, nphotons=1e6):
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
    print(' /')
    print('─' * 60)


# ---------------------------------------------------------------------------
# Slice plot
# ---------------------------------------------------------------------------

def make_slice_plot(grids, labels, outfile, log=False):
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
    except ImportError:
        print('matplotlib not available — skipping plot')
        return

    ncols = len(grids)
    fig, axes = plt.subplots(1, ncols, figsize=(6.2 * ncols, 5.5))
    if ncols == 1:
        axes = [axes]

    for ax, grid, label in zip(axes, grids, labels):
        # Density with AMR cell boundaries overlaid
        grid.slice_plot('z', 0.0, 'dens', ax=ax,
                        log=log, cmap='viridis',
                        show_leaf_boundaries=True,
                        boundary_color='white',
                        boundary_lw=0.25,
                        boundary_alpha=0.5)
        ax.set_title(label, fontsize=10)

    title = ('log$_{10}$ nH [cm$^{-3}$]' if log else 'nH [cm$^{-3}$]')
    fig.suptitle(f'Radial shell refinement — {title}', fontsize=11)

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
            )
    else:
        # ── build a single grid ───────────────────────────────────────────
        print(f'Building AMR grid  boxlen={args.boxlen}  rmax={rmax}')
        print(f'  density  : {args.density}  n0={args.n0:.2e} cm^-3')
        print(f'  T        : {args.temperature:.2e} K')
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
            make_slice_plot([grid],
                            [f'{label} spacing — {args.density} density'],
                            args.output,
                            log=args.plot_log)

        print_lart_hint(args.output, args.boxlen, args.temperature,
                        taumax=args.taumax)


if __name__ == '__main__':
    main()
