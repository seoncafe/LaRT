import numpy as np
import matplotlib.pyplot as plt

from make_amr_grid import AMRGrid


def make_uniform_sphere_grid(level_fine=3):
    boxlen = 100.0
    cx0 = cy0 = cz0 = boxlen / 2.0
    r_sphere = 40.0

    grid = AMRGrid(boxlen)

    # uniform-density sphere
    def nH_fn(x, y, z):
        r2 = (x - cx0)**2 + (y - cy0)**2 + (z - cz0)**2
        return 1.0 if r2 <= r_sphere**2 else 0.0

    # static velocity
    def vel_fn(x, y, z):
        return 0.0, 0.0, 0.0

    # geometry-aware refinement
    # - fully inside sphere  -> refine only if physics criterion is met
    # - not fully inside     -> refine unconditionally to level_fine
    grid.refine_sphere_by_physics(
        cx0, cy0, cz0, r_sphere,
        nH_fn=nH_fn,
        vel_fn=vel_fn,
        dens_threshold=0.1,
        vel_threshold=0.1,
        level_max=level_fine,
        nprobe=2,
    )

    grid.set_density(nH_fn)
    grid.set_temperature(1.0e4)
    grid.set_velocity(vel_fn)

    return grid


# ---- level 3 -----------------------------------------------------------
grid3 = make_uniform_sphere_grid(level_fine=3)
print("level=3:", grid3.level_counts())

fig, ax = plt.subplots(figsize=(6, 6))
grid3.slice_plot(axis="z", value=50.0, quantity="nH", ax=ax, log=False, show_leaf_boundaries=True, show_leaf_centers=True)
ax.set_title("Uniform sphere slice (level_max=3)")
plt.tight_layout()
plt.savefig("uniform_sphere_slice_l3.png", dpi=150)
plt.close(fig)


# ---- level 4 -----------------------------------------------------------
grid4 = make_uniform_sphere_grid(level_fine=4)
print("level=4:", grid4.level_counts())

fig, ax = plt.subplots(figsize=(6, 6))
grid4.slice_plot(axis="z", value=50.0, quantity="nH", ax=ax, log=False, show_leaf_boundaries=True, show_leaf_centers=True)
ax.set_title("Uniform sphere slice (level_max=4)")
plt.tight_layout()
plt.savefig("uniform_sphere_slice_l4.png", dpi=150)
plt.close(fig)

grid5 = make_uniform_sphere_grid(level_fine=5)
outfile = 'uniform_100kpc.dat'
grid5.write(outfile)
