import numpy as np
import matplotlib.pyplot as plt

from make_amr_grid import AMRGrid


def make_uniform_sphere_grid(level_min=1, level_max=3, refine_boundary=False):
    boxlen = 100.0
    cx0 = cy0 = cz0 = boxlen / 2.0
    r_sphere = boxlen/2.0
    #r_sphere = 40.0

    grid = AMRGrid(boxlen)

    # uniform-density sphere
    def gasDen_fn(x, y, z):
        r2 = (x - cx0)**2 + (y - cy0)**2 + (z - cz0)**2
        return 1.0 if r2 <= r_sphere**2 else 0.0

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
        gasDen_fn=gasDen_fn,
        vel_fn=vel_fn,
        dens_threshold=0.1,
        vel_threshold=0.1,
        level_min=level_min,
        level_max=level_max,
        nprobe=2,
        refine_boundary=refine_boundary,
    )

    grid.set_density(gasDen_fn)
    grid.set_temperature(1.0e4)
    grid.set_velocity(vel_fn)

    return grid


# ---- level 3 -----------------------------------------------------------
grid3 = make_uniform_sphere_grid(level_max=3)
print("level=3:", grid3.level_counts())

fig, ax = plt.subplots(figsize=(6.5, 6))
grid3.slice_plot(axis="z", value=50.0, quantity="gasDen", ax=ax, log=False, show_leaf_boundaries=True, show_leaf_centers=True)
ax.set_title("Uniform sphere slice (level_max=3)")
#plt.tight_layout()
plt.savefig("uniform_sphere_slice_l3.pdf", dpi=600)
plt.close(fig)


# ---- level 4 -----------------------------------------------------------
grid4 = make_uniform_sphere_grid(level_max=4)
print("level=4:", grid4.level_counts())

fig, ax = plt.subplots(figsize=(6.5, 6))
grid4.slice_plot(axis="z", value=50.0, quantity="gasDen", ax=ax, log=False, show_leaf_boundaries=True, show_leaf_centers=True)
ax.set_title("Uniform sphere slice (level_max=4)")
#plt.tight_layout()
plt.savefig("uniform_sphere_slice_l4.pdf", dpi=600)
plt.close(fig)

# ---- level 7  -----------------------------------------------------------
grid = make_uniform_sphere_grid(level_max=6,level_min=3,refine_boundary=True)
outfile = 'uniform_100kpc.fits.gz'
grid.write(outfile)
