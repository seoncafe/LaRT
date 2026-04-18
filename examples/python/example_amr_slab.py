"""
example_amr_slab.py
===================
Create an AMR data file for a plane-parallel slab (or multi-slab) geometry.

The slab is oriented perpendicular to the z-axis and centred at z = boxlen/2.
Cells within the slab thickness are refined to a higher AMR level.

Density profiles available:
  uniform    : nH = const inside the slab, 0 outside
  gaussian   : nH(z) = nH0 * exp(-(z - z0)^2 / (2 * h_s^2))
  exponential: nH(z) = nH0 * exp(-|z - z0| / h_s)

Velocity options:
  static     : v = 0
  gradient   : vz = v_amp * (z - z0) / (boxlen/2)   (shear / outflow along z)
  turbulent  : random Gaussian velocity for each leaf  (seed can be set)

Usage
-----
    python example_amr_slab.py                         # basic uniform slab
    python example_amr_slab.py --profile gaussian      # Gaussian slab
    python example_amr_slab.py --velocity gradient --v-amp 50
    python example_amr_slab.py --help
"""

import argparse
import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))
from make_amr_grid import AMRGrid


def make_slab(boxlen=100.0,
              slab_half=20.0,
              profile='uniform',
              nH0=1e-2,
              h_sigma=10.0,
              T=1e4,
              level_coarse=2,
              level_fine=3,
              velocity='static',
              v_amp=30.0,
              turb_seed=42,
              outfile='slab_amr.dat'):
    """
    Build and write a plane-parallel slab AMR grid.

    Parameters
    ----------
    boxlen : float
        Box side length [kpc].  Box occupies [0, boxlen]^3.
    slab_half : float
        Half-thickness of the slab refinement region [kpc].
        Cells with |z - boxlen/2| < slab_half are refined to ``level_fine``.
    profile : str
        Density profile: 'uniform', 'gaussian', or 'exponential'.
    nH0 : float
        Peak hydrogen number density [cm^-3].
    h_sigma : float
        Scale height [kpc] for Gaussian or exponential profile.
    T : float or callable
        Gas temperature [K].  Can be a callable ``T(x, y, z) -> float``.
    level_coarse : int
        AMR level for the region outside the slab (root = 0).
    level_fine : int
        AMR level inside the slab.
    velocity : str
        'static'   : v = 0
        'gradient' : vz proportional to z displacement from slab centre
        'turbulent': each leaf gets a random Gaussian (vx, vy, vz)
    v_amp : float
        Velocity amplitude [km/s].
    turb_seed : int
        Random seed for turbulent velocity field.
    outfile : str
        Output filename.
    """
    z0 = boxlen / 2.0   # slab centre

    # ---- Build octree --------------------------------------------------
    grid = AMRGrid(boxlen)

    # Coarse level everywhere
    grid.refine(lambda c: True, level_max=level_coarse)

    # Fine level inside the slab
    grid.refine_slab('z', z0 - slab_half, z0 + slab_half, level_max=level_fine)

    # ---- Density -------------------------------------------------------
    if profile == 'uniform':
        def nH_fn(x, y, z):
            return nH0 if abs(z - z0) <= slab_half else 0.0

    elif profile == 'gaussian':
        def nH_fn(x, y, z):
            return nH0 * np.exp(-0.5 * ((z - z0) / h_sigma)**2)

    elif profile == 'exponential':
        def nH_fn(x, y, z):
            return nH0 * np.exp(-abs(z - z0) / h_sigma)

    else:
        raise ValueError(f"Unknown profile: {profile!r}. "
                         "Choose 'uniform', 'gaussian', or 'exponential'.")

    grid.set_density(nH_fn)

    # ---- Temperature ---------------------------------------------------
    grid.set_temperature(T)

    # ---- Velocity ------------------------------------------------------
    if velocity == 'static':
        pass

    elif velocity == 'gradient':
        # vz = v_amp * (z - z0) / (boxlen/2), vx = vy = 0
        half = boxlen / 2.0
        def vel_grad(x, y, z):
            vz = v_amp * (z - z0) / half
            return 0.0, 0.0, vz
        grid.set_velocity(vel_grad)

    elif velocity == 'turbulent':
        # Each leaf draws (vx,vy,vz) from N(0, v_amp/sqrt(3))
        rng    = np.random.default_rng(turb_seed)
        leaves = grid.leaves()
        sigma  = v_amp / np.sqrt(3.0)
        vels   = rng.normal(0.0, sigma, size=(len(leaves), 3))
        for lf, (vx, vy, vz) in zip(leaves, vels):
            lf.vx, lf.vy, lf.vz = vx, vy, vz

    else:
        raise ValueError(f"Unknown velocity type: {velocity!r}. "
                         "Choose 'static', 'gradient', or 'turbulent'.")

    # ---- Write and report ----------------------------------------------
    grid.write(outfile)
    print(grid.info())
    return grid


def make_laRT_input_slab(outfile_dat, tau=1e4, nphotons=1e6, outfile_in=None):
    """Write a matching LaRT input file for the slab (uniform source)."""
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
        description='Create an AMR slab data file for LaRT v2.00.')
    parser.add_argument('--boxlen',      type=float, default=100.0,
                        help='Box side length [kpc] (default: 100)')
    parser.add_argument('--slab-half',   type=float, default=20.0,
                        help='Slab half-thickness for refinement [kpc] (default: 20)')
    parser.add_argument('--profile',     choices=['uniform', 'gaussian', 'exponential'],
                        default='uniform',
                        help='Density profile shape (default: uniform)')
    parser.add_argument('--nH0',         type=float, default=1e-2,
                        help='Peak HI density [cm^-3] (default: 1e-2)')
    parser.add_argument('--h-sigma',     type=float, default=10.0,
                        help='Density scale height [kpc] (default: 10)')
    parser.add_argument('--temperature', type=float, default=1e4,
                        help='Gas temperature [K] (default: 1e4)')
    parser.add_argument('--level-coarse', type=int, default=2,
                        help='Coarse AMR level (default: 2)')
    parser.add_argument('--level-fine',   type=int, default=3,
                        help='Fine AMR level inside slab (default: 3)')
    parser.add_argument('--velocity',    choices=['static', 'gradient', 'turbulent'],
                        default='static',
                        help='Velocity field type (default: static)')
    parser.add_argument('--v-amp',       type=float, default=30.0,
                        help='Velocity amplitude [km/s] (default: 30)')
    parser.add_argument('--turb-seed',   type=int, default=42,
                        help='Random seed for turbulent velocity (default: 42)')
    parser.add_argument('--outfile',     default='slab_amr.dat',
                        help='Output AMR data file (default: slab_amr.dat)')
    parser.add_argument('--tau',         type=float, default=1e4,
                        help='taumax for the LaRT input file (default: 1e4)')
    parser.add_argument('--nphotons',    type=float, default=1e6,
                        help='Number of photons (default: 1e6)')
    parser.add_argument('--write-input', action='store_true',
                        help='Also write a matching LaRT .in file')
    args = parser.parse_args()

    grid = make_slab(
        boxlen       = args.boxlen,
        slab_half    = args.slab_half,
        profile      = args.profile,
        nH0          = args.nH0,
        h_sigma      = args.h_sigma,
        T            = args.temperature,
        level_coarse = args.level_coarse,
        level_fine   = args.level_fine,
        velocity     = args.velocity,
        v_amp        = args.v_amp,
        turb_seed    = args.turb_seed,
        outfile      = args.outfile,
    )

    if args.write_input:
        make_laRT_input_slab(args.outfile, tau=args.tau, nphotons=args.nphotons)
