"""
example_amr_sphere.py
=====================
Create an AMR data file for a spherical HI cloud.

Two refinement levels are used:
  - Level 2 (coarse, half-width = boxlen/8)  : outer region
  - Level 3 (fine,   half-width = boxlen/16) : inner sphere

Density profile  : Gaussian   nH(r) = nH0 * exp(-r^2 / (2*r_s^2))
Temperature      : uniform    T = 1e4 K
Velocity options : static  |  Hubble-like expansion  |  solid-body rotation

Usage
-----
    python example_amr_sphere.py                     # static sphere
    python example_amr_sphere.py --velocity hubble   # expanding sphere
    python example_amr_sphere.py --velocity rotate   # rotating sphere
    python example_amr_sphere.py --help
"""

import argparse
import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))
from make_amr_grid import AMRGrid


def make_sphere(boxlen=100.0,
                r_refine=40.0,
                r_sigma=30.0,
                nH0=3.2e-4,
                T=1e4,
                level_coarse=2,
                level_fine=3,
                velocity='static',
                v_amp=20.0,
                outfile='sphere_amr.dat'):
    """
    Build and write an AMR sphere grid.

    Parameters
    ----------
    boxlen : float
        Box side length [kpc].  Box occupies [0, boxlen]^3.
    r_refine : float
        Refinement sphere radius [kpc] from box centre.
        Cells inside are resolved to ``level_fine``; outside to ``level_coarse``.
    r_sigma : float
        Gaussian density scale radius [kpc].
    nH0 : float
        Peak hydrogen number density at r=0 [cm^-3].
    T : float
        Gas temperature [K] (uniform).
    level_coarse : int
        AMR level for the outer region (root = 0).
    level_fine : int
        AMR level for the inner sphere.
    velocity : str
        'static'  : v = 0
        'hubble'  : v_r = v_amp * r / (boxlen/2)  [km/s]  (outflow)
        'rotate'  : solid-body rotation about z-axis with amplitude v_amp [km/s]
    v_amp : float
        Velocity amplitude [km/s].
    outfile : str
        Output filename.
    """
    cx0 = cy0 = cz0 = boxlen / 2.0   # box centre

    # ---- Build octree --------------------------------------------------
    grid = AMRGrid(boxlen)

    # Coarse level everywhere first
    grid.refine(lambda c: True, level_max=level_coarse)

    # Fine level inside the refinement sphere
    grid.refine_sphere(cx0, cy0, cz0, r_refine, level_max=level_fine)

    # ---- Density (Gaussian) -------------------------------------------
    def nH_fn(x, y, z):
        r2 = (x-cx0)**2 + (y-cy0)**2 + (z-cz0)**2
        return nH0 * np.exp(-r2 / (2.0 * r_sigma**2))

    grid.set_density(nH_fn)

    # ---- Temperature (uniform) ----------------------------------------
    grid.set_temperature(T)

    # ---- Velocity ------------------------------------------------------
    if velocity == 'static':
        pass   # default vx=vy=vz=0

    elif velocity == 'hubble':
        # Hubble-like: v_r = v_amp * r / r_box, directed radially outward
        r_box = boxlen / 2.0
        def vel_hubble(x, y, z):
            dx, dy, dz = x-cx0, y-cy0, z-cz0
            r = np.sqrt(dx**2 + dy**2 + dz**2)
            if r < 1e-10:
                return 0.0, 0.0, 0.0
            scale = v_amp * r / r_box / r
            return scale*dx, scale*dy, scale*dz
        grid.set_velocity(vel_hubble)

    elif velocity == 'rotate':
        # Solid-body rotation about z-axis: v = v_amp * (r_perp / r_box) * phi_hat
        r_box = boxlen / 2.0
        def vel_rotate(x, y, z):
            dx, dy = x-cx0, y-cy0
            r_perp = np.sqrt(dx**2 + dy**2)
            if r_perp < 1e-10:
                return 0.0, 0.0, 0.0
            scale = v_amp * r_perp / r_box / r_perp
            return -scale*dy, scale*dx, 0.0
        grid.set_velocity(vel_rotate)

    else:
        raise ValueError(f"Unknown velocity type: {velocity!r}. "
                         "Choose 'static', 'hubble', or 'rotate'.")

    # ---- Write and report ---------------------------------------------
    grid.write(outfile)
    print(grid.info())
    return grid


def make_laRT_input(outfile_dat, tau=1e4, nphotons=1e6, outfile_in=None):
    """Write a matching LaRT input file."""
    if outfile_in is None:
        base = os.path.splitext(outfile_dat)[0]
        outfile_in = base + '.in'
    boxlen = 100.0
    center = boxlen / 2.0
    content = (
        f'&parameters\n'
        f' par%use_amr_grid  = .true.\n'
        f' par%amr_type      = \'generic\'\n'
        f' par%amr_file      = \'{outfile_dat}\'\n'
        f' par%distance_unit = \'kpc\'\n'
        f' par%no_photons    = {nphotons:.0e}\n'
        f' par%taumax        = {tau:.1e}\n'
        f' par%DGR           = 0.0\n'
        f' par%spectral_type = \'voigt\'\n'
        f' par%source_geometry = \'point\'\n'
        f' par%xs_point      = {center:.1f}\n'
        f' par%ys_point      = {center:.1f}\n'
        f' par%zs_point      = {center:.1f}\n'
        f' par%nxfreq  = 121\n'
        f' par%nxim    = 64\n'
        f' par%nyim    = 64\n'
        f' par%nprint  = 100000\n'
        f' par%out_merge = .false.\n'
        f'/\n'
    )
    with open(outfile_in, 'w') as f:
        f.write(content)
    print(f'LaRT input file: {outfile_in}')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Create an AMR sphere data file for LaRT v2.00.')
    parser.add_argument('--boxlen',   type=float, default=100.0,
                        help='Box side length [kpc] (default: 100)')
    parser.add_argument('--r-refine', type=float, default=40.0,
                        help='Refinement sphere radius [kpc] (default: 40)')
    parser.add_argument('--r-sigma',  type=float, default=30.0,
                        help='Gaussian density scale radius [kpc] (default: 30)')
    parser.add_argument('--nH0',      type=float, default=3.2e-4,
                        help='Peak HI density [cm^-3] (default: 3.2e-4)')
    parser.add_argument('--temperature', type=float, default=1e4,
                        help='Gas temperature [K] (default: 1e4)')
    parser.add_argument('--level-coarse', type=int, default=2,
                        help='Coarse AMR level (default: 2)')
    parser.add_argument('--level-fine',   type=int, default=3,
                        help='Fine AMR level inside refinement sphere (default: 3)')
    parser.add_argument('--velocity', choices=['static', 'hubble', 'rotate'],
                        default='static',
                        help='Velocity field type (default: static)')
    parser.add_argument('--v-amp',    type=float, default=20.0,
                        help='Velocity amplitude [km/s] (default: 20)')
    parser.add_argument('--outfile',  default='sphere_amr.dat',
                        help='Output AMR data file (default: sphere_amr.dat)')
    parser.add_argument('--tau',      type=float, default=1e4,
                        help='taumax for the LaRT input file (default: 1e4)')
    parser.add_argument('--nphotons', type=float, default=1e6,
                        help='Number of photons (default: 1e6)')
    parser.add_argument('--write-input', action='store_true',
                        help='Also write a matching LaRT .in file')
    args = parser.parse_args()

    grid = make_sphere(
        boxlen      = args.boxlen,
        r_refine    = args.r_refine,
        r_sigma     = args.r_sigma,
        nH0         = args.nH0,
        T           = args.temperature,
        level_coarse= args.level_coarse,
        level_fine  = args.level_fine,
        velocity    = args.velocity,
        v_amp       = args.v_amp,
        outfile     = args.outfile,
    )

    if args.write_input:
        make_laRT_input(args.outfile, tau=args.tau, nphotons=args.nphotons)
