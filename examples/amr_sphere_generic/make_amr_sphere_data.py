import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append('../python')

from make_amr_grid import AMRGrid

def make_sphere_grid(r_scale = 1.0, level_min=1, level_max=3, refine_boundary=False):
    boxlen = 2.0
    cx0 = cy0 = cz0 = 0.0
    r_sphere = boxlen/2.0
    dens0 = 1.0

    grid = AMRGrid(boxlen)

    # Gaussian (r_scale > 0) or constant (r_scale <= 0) density inside the sphere
    def dens_fn(x, y, z):
        r2 = x**2 + y**2 + z**2
        if r2 > r_sphere**2:
            return 0.0
        return dens0 * np.exp(-r2 / (2.0 * r_scale**2)) if r_scale > 0.0 else dens0

    # static velocity
    def vel_fn(x, y, z):
        return 0.0, 0.0, 0.0

    # Geometry-aware refinement with forced boundary resolution.
    #
    # refine_boundary=True is required here because the density field is a
    # hard step function (1 inside, 0 outside).  At high refinement levels
    # (level_max >= 5) cells near the surface become small enough that all
    # 8 probe points of a boundary cell can fall on the same side of the
    # sphere, making the gradient criterion miss that cell.  Forcing every
    # cell that geometrically intersects the surface guarantees a sharp,
    # gap-free boundary in the slice plot.
    #
    # Use refine_boundary=False (the default) when the density field varies
    # smoothly and the gradient criterion is sufficient to capture the region
    # of interest without special treatment of the geometric boundary.
    grid.refine_sphere_by_physics(
        cx0, cy0, cz0, r_sphere,
        dens_fn=dens_fn,
        vel_fn=vel_fn,
        dens_threshold=0.1,
        vel_threshold=0.1,
        level_min=level_min,
        level_max=level_max,
        nprobe=2,
        refine_boundary=refine_boundary,
    )

    grid.set_density(dens_fn)
    grid.set_temperature(1.0e4)
    grid.set_velocity(vel_fn)

    return grid


r_scale = 0.0
# ---- level 3 -----------------------------------------------------------
#grid = make_sphere_grid(r_scale=r_scale, level_max=3)
#print("level=3:", grid.level_counts())
#outfile = 'uniform_amr_sphere.fits.gz'
#grid.write(outfile)

# ---- level 4 -----------------------------------------------------------
#grid = make_sphere_grid(r_scale=r_scale, level_max=4)
#print("level=4:", grid.level_counts())
#outfile = 'uniform_amr_sphere.fits.gz'
#grid.write(outfile)

# ---- level 7  -----------------------------------------------------------
grid = make_sphere_grid(r_scale=r_scale, level_min=3,level_max=7,refine_boundary=True)
print("level=7:", grid.level_counts())
outfile = 'uniform_amr_sphere.fits.gz'
grid.write(outfile)
