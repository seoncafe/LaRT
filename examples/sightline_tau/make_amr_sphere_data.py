#!/usr/bin/env python3
"""
Generate `amr_sphere.fits.gz` for the sight-line tau AMR example.

A uniform sphere of radius 1.0 (code units) inside a box of side 2.0,
with the boundary refined to level_max so that the sphere surface is
captured cleanly.  Density and temperature are constant inside the
sphere and zero outside.

Run once:
    python make_amr_sphere_data.py
The resulting `amr_sphere.fits.gz` is consumed by `amr_sightline.in`.
"""
import os
import sys
import numpy as np

# pull in the AMR_grid helper from the python tree
HERE     = os.path.dirname(os.path.abspath(__file__))
PYDIR    = os.path.normpath(os.path.join(HERE, '..', '..', 'python', 'AMR_grid'))
sys.path.insert(0, PYDIR)
from make_amr_grid import AMRGrid


def make_uniform_sphere(level_min=3, level_max=6):
    boxlen   = 2.0
    cx0 = cy0 = cz0 = 0.0
    r_sphere = 1.0
    dens0    = 1.0     # arbitrary; LaRT renormalises by par%taumax

    grid = AMRGrid(boxlen)

    def dens_fn(x, y, z):
        if x*x + y*y + z*z > r_sphere * r_sphere:
            return 0.0
        return dens0

    def vel_fn(x, y, z):
        return 0.0, 0.0, 0.0

    grid.refine_sphere_by_physics(
        cx0, cy0, cz0, r_sphere,
        dens_fn=dens_fn,
        vel_fn=vel_fn,
        dens_threshold=0.1,
        vel_threshold=0.1,
        level_min=level_min,
        level_max=level_max,
        nprobe=2,
        refine_boundary=True,
    )

    grid.set_density(dens_fn)
    grid.set_temperature(1.0e4)
    grid.set_velocity(vel_fn)
    return grid


if __name__ == '__main__':
    grid = make_uniform_sphere(level_min=3, level_max=6)
    print('level counts:', grid.level_counts())
    out = os.path.join(HERE, 'amr_sphere.fits.gz')
    grid.write(out)
    print('wrote', out)
