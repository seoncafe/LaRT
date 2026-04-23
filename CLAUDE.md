# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Directory Is

`LaRT_v2.00/` is the AMR (Adaptive Mesh Refinement) extension of `LaRT_combined`. It adds an octree grid alongside the existing Cartesian grid, selectable at runtime via `par%use_amr_grid = .true.`. The broader project context (build system, input format, Cartesian architecture) is documented in the parent `CLAUDE.md` at `../CLAUDE.md`.

## Build & Run

```bash
# Build (auto-detects compiler: mpiifort > mpiifx > mpif90)
make              # produces LaRT.x
make clean        # removes .o and .mod files
make cleanall     # removes .o, .mod, and executables
make DEBUG=1      # full debug flags + runtime checks
make F90=mpif90   # force GNU compiler

# Auxiliary converter: RAMSES output → generic text AMR format
make ramses2fits  # produces convert_ramses_to_generic_fits.x

# Run
mpirun -np 8 LaRT.x input.in

# AMR workflow: generate grid file first, then run simulation
cd examples/python
python example_amr_sphere.py      # writes sphere_amr.dat (sphere geometry)
python example_amr_slab.py        # writes slab AMR .dat files
cd ../amr_sphere_generic
mpirun -np 4 ../../LaRT.x sphere_amr_tau4_point.in
```

The Makefile cleans `.o/.mod` before linking (`default: clean main`), so there are no incremental builds — every `make` recompiles everything.

## Architecture: Dual-Grid Dispatch

The key architectural pattern is **module-level procedure pointers** declared in `define.f90` (lines ~551–704). `setup_procedure()` in `setup.f90` assigns them at startup based on `par%use_amr_grid` and other flags. All callers (`run_simulation_mod.f90`, `main.f90`, peel-off routines) call these pointers, never the concrete functions directly.

```
define.f90          → declares procedure pointers (raytrace_to_tau, scatter_dust, ...)
setup.f90           → setup_procedure() assigns them based on par%use_amr_grid
main.f90            → calls grid_create / grid_create_amr, then procedure pointers
run_simulation_mod  → calls run_simulation (pointer), raytrace_to_tau (pointer), ...
```

### Cartesian path (`par%use_amr_grid = .false.`)
- Grid data: `grid_type` in `define.f90`
- Grid build: `grid_mod_car.f90` → `grid_create(grid)`
- Raytrace: `raytrace_car.f90`
- Scatter: `scattering_car.f90`
- Peel-off: `peelingoff_rect.f90` / `peelingoff_heal.f90`

### AMR path (`par%use_amr_grid = .true.`)
- Grid data: `amr_grid_type` in `octree_mod.f90` (global `amr_grid` instance)
- Grid build: `grid_mod_amr.f90` → `grid_create_amr(grid)` + `amr_sync_to_grid(grid)`
- Octree: `octree_mod.f90` — `amr_build_tree`, `amr_build_neighbors`, `amr_next_leaf`
- Raytrace: `raytrace_amr.f90` — O(1) cell crossings via precomputed face-neighbor table
- Scatter: `scattering_amr.f90` — per-cell Voigt parameters (`voigt_a(il)`, `Dfreq(il)`); selects Stokes vs. no-Stokes variant based on `par%use_stokes`
- Peel-off: `peelingoff_amr.f90`
- Reader: `read_ramses_amr.f90` — RAMSES binary and generic text formats

`amr_sync_to_grid(grid)` copies scalar grid geometry (box size, frequency axes) from `amr_grid` into the Cartesian `grid` variable so that shared output/observer routines work unchanged.

## AMR Data Structures

`amr_grid_type` (defined in `octree_mod.f90`) uses flat arrays indexed by cell index (1-based). Internal cells and leaf cells share the same index space:
- `ileaf(icell)` > 0 → leaf index; 0 → internal cell
- `icell_of_leaf(il)` → cell index for leaf `il`
- `neighbor(iface, icell)`: precomputed face neighbors; face convention: 1=+x, 2=-x, 3=+y, 4=-y, 5=+z, 6=-z; 0 = outside box
- Physical data (`rhokap`, `Dfreq`, `voigt_a`, `vfx/y/z`) indexed by leaf index `il`

All large arrays use MPI-3 shared memory (one copy per node) via `create_shared_mem`. Output arrays (`Jout`, `Jin`, `Jabs`) are per-rank regular allocations.

## AMR Input Parameters

```fortran
&parameters
 par%use_amr_grid  = .true.
 par%amr_type      = 'generic'    ! 'generic' or 'ramses' (code default: 'ramses' — always specify explicitly)
 par%amr_file      = 'grid.dat'  ! for 'ramses': path to snapshot directory; for 'generic': path to .dat file
 par%distance_unit = 'kpc'        ! unit of leaf-cell coordinates in data file
 par%taumax        = 1.0e4        ! tau from box center to +z edge along z-axis
 par%source_geometry = 'point'
 par%xs_point = 50.0              ! point source coords in code units (box center)
 par%ys_point = 50.0
 par%zs_point = 50.0
/
```

`taumax` convention for AMR: optical depth from box center to the +z boundary, traversed along the z-axis (matches Cartesian convention where `taumax` = tau from center to sphere surface at `rmax`).

## Generic AMR Data File Format

Plain text consumed by `read_ramses_amr.f90`:
```
# comment lines start with #
BOXLEN  <box_length_in_code_units>
NLEAF   <number_of_leaf_cells>
<x> <y> <z> <level> <nH[cm^-3]> <T[K]> <vx[km/s]> <vy[km/s]> <vz[km/s]>
...
```

Leaf-cell coordinates are in the same unit as `BOXLEN` (physical unit given by `par%distance_unit`). Origin is at the box corner; range is `[0, BOXLEN]`.

## Generating AMR Data Files (Python)

`examples/python/make_amr_grid.py` provides the `AMRGrid` class:

```python
from make_amr_grid import AMRGrid
import numpy as np

grid = AMRGrid(boxlen=100.0)   # kpc
grid.refine(lambda c: c.dist(50, 50, 50) < 40, level_max=3)
grid.set_density(lambda x, y, z: 1e-3 * np.exp(-((x-50)**2+(y-50)**2+(z-50)**2)/(2*20**2)))
grid.set_temperature(lambda x, y, z: 1e4)
grid.write('my_grid.dat')
```

`set_density/temperature/velocity` accept vectorised functions `fn(cx_array, cy_array, cz_array)` for performance on large grids. `make_amr_grid.py` is a library — import `AMRGrid`; `example_amr_sphere.py` and `example_amr_slab.py` are ready-made generators.

## Refactoring Plan (In Progress)

The CLAUDE.md in the parent directory (`../CLAUDE.md`) contains the full class-based OOP refactoring plan. The goal is to replace the global procedure pointer pattern with Fortran 2003 `class(grid_base_type), allocatable` dispatch. The `par%grid_type = 'car' | 'amr'` string parameter will replace `par%use_amr_grid`. Implementation order and design are detailed in the parent CLAUDE.md under "Class-Based Refactoring Plan".

## Key Files Quick Reference

| File | Role |
|------|------|
| `define.f90` | All global types + procedure pointer declarations |
| `setup.f90` | `read_input()` (namelist) + `setup_procedure()` (pointer assignments) |
| `main.f90` | MPI init, grid creation branch (AMR vs. Cartesian), output |
| `octree_mod.f90` | `amr_grid_type`, tree build, neighbor table, cell lookup |
| `grid_mod_amr.f90` | AMR grid setup: reads data, builds tree, normalizes opacity |
| `raytrace_amr.f90` | Photon propagation through AMR cells |
| `scattering_amr.f90` | Resonance/dust scattering with Stokes support |
| `peelingoff_amr.f90` | Peel-off observer calculation for AMR |
| `read_ramses_amr.f90` | Reader for RAMSES binary and generic text AMR formats |
| `run_simulation_mod.f90` | Main photon loop (master-slave or equal-split) |

Supporting resources:
- `docs/LaRT_AMR_description.pdf` — technical description of the AMR implementation
- `docs/RAMSES_data_structure.pdf` — RAMSES AMR data structure reference
- `examples/amr_sphere_generic/` — generic-format sphere validation cases (AMR vs. Cartesian)
- `examples/amr_ramses/` — RAMSES snapshot example

## Validation

AMR results are validated against equivalent Cartesian runs. Cartesian comparison uses `rmax = L_box/2` (no `distance_unit`). The ~3–5% amplitude difference in `Jout_max` at tau=1e2 and tau=1e4 is expected from octree volume-averaging at the sphere boundary, not a normalization error. Peak velocity and `N(HI)_pole` agree exactly.
