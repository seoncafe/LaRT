"""
example_amr_sphere.py
=====================
Create an AMR data file for a spherical HI cloud.

Two refinement levels are used:
  - Level 2 (coarse, half-width = boxlen/8)  : background level
  - Level 3 (fine,   half-width = boxlen/16) : cells not fully contained in
    the refinement sphere, plus any fully contained cells that satisfy the
    density/velocity refinement criterion.

Density profile  : Gaussian   dens(r) = dens0 * exp(-r^2 / (2*r_s^2))
Temperature      : uniform    T = 1e4 K
Velocity options : static  |  Hubble-like expansion  |  solid-body rotation

Usage
-----
    python example_amr_sphere.py
    python example_amr_sphere.py --velocity hubble
    python example_amr_sphere.py --velocity rotate --dens-threshold 0.05 --nprobe 4
    python example_amr_sphere.py --help
"""

import argparse
import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))
from AMR_grid import AMRGrid


def make_sphere(boxlen=2.0,
                r_refine=1.0,
                r_sigma=1.0,
                dens0=1.0,
                T=1e4,
                level_coarse=2,
                level_fine=3,
                velocity='static',
                v_amp=1.0,
                dens_threshold=0.1,
                vel_threshold=0.1,
                nprobe=2,
                outfile='sphere_amr.dat'):
    """
    Build and write an AMR sphere grid.

    Parameters
    ----------
    boxlen : float
        Box side length [code unit].  Box occupies [-boxlen/2, boxlen/2]^3.
    r_refine : float
        Refinement sphere radius [code unit] from box centre.
        Cells not fully contained in the sphere are always refined to
        ``level_fine``. Cells fully contained in the sphere are refined
        toward ``level_fine`` only when the physics criterion is satisfied.
    r_sigma : float
        Gaussian density scale radius [code unit].
    dens0 : float
        Peak gas density (including hydrogen) at r=0 [cm^-3].
    T : float
        Gas temperature [K] (uniform).
    level_coarse : int
        AMR level for the background region (root = 0).
    level_fine : int
        Maximum AMR level associated with the sphere.
    velocity : str
        'static'  : v = 0
        'hubble'  : v_r = v_amp * r / (boxlen/2)  [km/s]  (outflow)
        'rotate'  : solid-body rotation about z-axis with amplitude v_amp [km/s]
    v_amp : float
        Velocity amplitude [km/s].
    dens_threshold : float
        Relative density-gradient threshold used for refinement of cells fully
        contained inside the sphere.
    vel_threshold : float
        Relative velocity-gradient threshold used for refinement of cells fully
        contained inside the sphere.
    nprobe : int
        Number of probe points per axis used to evaluate the refinement
        criteria inside each cell.
    outfile : str
        Output filename.
    """
    cx0 = cy0 = cz0 = 0.0

    def dens_fn(x, y, z):
        r2 = x**2 + y**2 + z**2
        return dens0 * np.exp(-r2 / (2.0 * r_sigma**2)) if r2 < boxlen**2 and r_sigma > 0.0 else 0.0

    vel_fn = None
    if velocity == 'static':
        pass
    elif velocity == 'hubble':
        r_box = boxlen / 2.0
        def vel_hubble(x, y, z):
            dx, dy, dz = x-cx0, y-cy0, z-cz0
            r = np.sqrt(dx**2 + dy**2 + dz**2)
            if r < 1e-10:
                return 0.0, 0.0, 0.0
            scale = v_amp * r / r_box / r
            return scale*dx, scale*dy, scale*dz
        vel_fn = vel_hubble
    elif velocity == 'rotate':
        r_box = boxlen / 2.0
        def vel_rotate(x, y, z):
            dx, dy = x-cx0, y-cy0
            r_perp = np.sqrt(dx**2 + dy**2)
            if r_perp < 1e-10:
                return 0.0, 0.0, 0.0
            scale = v_amp * r_perp / r_box / r_perp
            return -scale*dy, scale*dx, 0.0
        vel_fn = vel_rotate
    else:
        raise ValueError(f"Unknown velocity type: {velocity!r}. Choose 'static', 'hubble', or 'rotate'.")

    grid = AMRGrid(boxlen)
    grid.refine(lambda c: True, level_max=level_coarse)
    grid.refine_sphere_by_physics(
        cx0, cy0, cz0, r_refine,
        dens_fn=dens_fn,
        vel_fn=vel_fn,
        dens_threshold=dens_threshold,
        vel_threshold=vel_threshold,
        level_max=level_fine,
        nprobe=nprobe,
    )

    grid.set_density(dens_fn)
    grid.set_temperature(T)
    if vel_fn is not None:
        grid.set_velocity(vel_fn)

    grid.write(outfile)
    print(grid.info())
    return grid


def make_laRT_input(outfile_dat, tau=1e4, nphotons=1e6, outfile_in=None):
    """Write a matching LaRT input file."""
    if outfile_in is None:
        base = os.path.splitext(outfile_dat)[0]
        outfile_in = base + '.in'
    content = (
        f"&parameters\n"
        f" par%use_amr_grid  = .true.\n"
        f" par%amr_type      = 'generic'\n"
        f" par%amr_file      = '{outfile_dat}'\n"
        f" par%distance_unit = 1.0\n"
        f" par%no_photons    = {nphotons:.0e}\n"
        f" par%taumax        = {tau:.1e}\n"
        f" par%DGR           = 0.0\n"
        f" par%spectral_type = 'voigt'\n"
        f" par%source_geometry = 'point'\n"
        f" par%xs_point      = 0.0\n"
        f" par%ys_point      = 0.0\n"
        f" par%zs_point      = 0.0\n"
        f" par%nxfreq  = 121\n"
        f" par%nxim    = 100\n"
        f" par%nyim    = 100\n"
        f" par%nprint  = 100000\n"
        f" par%out_merge = .false.\n"
        f"/\n"
    )
    with open(outfile_in, 'w') as f:
        f.write(content)
    print(f'LaRT input file: {outfile_in}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Create an AMR sphere data file for LaRT v2.00.')
    parser.add_argument('--boxlen',   type=float, default=2.0,
                        help='Box side length (default: 2.0)')
    parser.add_argument('--r-refine', type=float, default=1.0,
                        help='Refinement sphere radius (default: 1.0)')
    parser.add_argument('--r-sigma',  type=float, default=1.0,
                        help='Gaussian density scale radius (default: 1.0)')
    parser.add_argument('--dens0',        type=float, default=1.0,
                        help='Peak gas density [cm^-3] (default: 1.0)')
    parser.add_argument('--temperature', type=float, default=1e4,
                        help='Gas temperature [K] (default: 1e4)')
    parser.add_argument('--level-coarse', type=int, default=2,
                        help='Coarse AMR level (default: 2)')
    parser.add_argument('--level-fine',   type=int, default=3,
                        help='Fine AMR level inside refinement sphere (default: 3)')
    parser.add_argument('--velocity', choices=['static', 'hubble', 'rotate'],
                        default='static',
                        help='Velocity field type (default: static)')
    parser.add_argument('--v-amp',    type=float, default=1.0,
                        help='Velocity amplitude [km/s] (default: 1.0)')
    parser.add_argument('--dens-threshold', type=float, default=0.1,
                        help='Density-gradient refinement threshold (default: 0.1)')
    parser.add_argument('--vel-threshold',  type=float, default=0.1,
                        help='Velocity-gradient refinement threshold (default: 0.1)')
    parser.add_argument('--nprobe', type=int, default=2,
                        help='Probe count per axis for refinement criteria (default: 2)')
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
        dens0       = args.dens0,
        T           = args.temperature,
        level_coarse= args.level_coarse,
        level_fine  = args.level_fine,
        velocity    = args.velocity,
        v_amp       = args.v_amp,
        dens_threshold = args.dens_threshold,
        vel_threshold  = args.vel_threshold,
        nprobe      = args.nprobe,
        outfile     = args.outfile,
    )

    if args.write_input:
        make_laRT_input(args.outfile, tau=args.tau, nphotons=args.nphotons)
