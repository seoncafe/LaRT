# LaRT v2.0

LaRT v2.0 is a radiative transfer code for resonance-line scattering. In addition to Lyα, it is designed to handle resonance and fluorescence scattering of a variety of metallic lines.

This version adds support for adaptive mesh refinement (AMR) grids. Simulation data from RAMSES or Illustris/IllustrisTNG are first converted to a generic AMR file (HDF5, gzipped FITS, or plain text) via standalone converters (`convert_ramses_to_generic.{x,py}` and `convert_illustris_to_generic.py`). **The AMR implementation is currently under active testing and validation.**

## Status

- AMR support: octree grid with leaf physical quantities; testing and validation are ongoing.
- **cell-by-cell core-skip**: `xcrit` is now recomputed at every scattering from the local cell's `voigt_a × rhokap × dl_face` ([Smith+15](https://ui.adsabs.harvard.edu/abs/2015MNRAS.449.4336S/abstract) Eq.35). This replaces the volume-averaged `xcrit` for inhomogeneous (e.g. galaxy-scale) boxes where the old prescription gave `xcrit = 0` and photons trapped in dense regions never escaped. The old behavior is retained behind `par%core_skip_global = .true.` for benchmarking.
- Clump overlap handling added: file-loaded or internally generated overlapping clumps are handled with an event-based multi-component raytrace; enabled via `par%clump_allow_overlap = .true.` for internally generated populations.
- **HDF5 I/O support (output and input)** is now available alongside FITS for all output streams (spectrum, peel-off, sight-line tau, CALCJ/CALCP arrays) and all gridded inputs (density, temperature, velocity, emissivity, clump files). The default is HDF5 (`par%file_format = 'hdf5'`); set `par%file_format = 'fits'` to fall back to `.fits.gz`. Build with `make HDF5=1` (default) or `make HDF5=0` to skip the HDF5 link dependency. See `python/lart_io.py` for a format-agnostic Python reader/converter.
- **CALCJ / CALCP / CALCPnew on the AMR grid** (2026-07-02): the internal mean-intensity spectrum `J(x)` and the scattering rate per atom `P_alpha` are now accumulated on the octree grid as well. Deposits are **position-binned** (segment midpoint / scattering position) with volume-weighted normalization from the true leaf/bin overlap, so radial profiles stay correct even in coarsely refined cores. Output sections: `Jx_AMR`/`Pa_AMR` (per leaf) or `Jx_1D`/`Pa_1D` etc. (binned profiles) with bin-axis keywords (`nr/rmax/dr`, `nz/zmin/dz`); `CALCPnew` sections carry the `_new` suffix. AMR radial/cylindrical binning requires `par%rmax > 0` in the input file. Validated against Cartesian runs to 1-2% (rms) at matched resolution. `python/read_lart.py` loads all CALC sections and provides `plot_J_profile()` / `plot_Pa_profile()`.
- **Bug fix** (2026-07-02): inputs without an explicit `par%geometry` were silently promoted to `'sphere'`, which inflated `nx=ny=1` `xy_periodic` slabs to full N^3 cubes with sphere density zeroing. The sphere dimension-forcing is now skipped for `xy_periodic` runs; re-run any Cartesian slab results produced since late April 2026.
- The generic AMR **plain-text format** now accepts `#` comment lines and `BOXLEN`/`NLEAF` keyword headers (the legacy `"<nleaf> <boxlen>"` first line still works).
- The usage guide is still being prepared.

## References

For scientific background and related methodology, see:

- [Seon & Kim (2020), ApJS, 250, 9](https://ui.adsabs.harvard.edu/abs/2020ApJS..250....9S/abstract)
- [Seon et al. (2022), ApJS, 259, 3](https://ui.adsabs.harvard.edu/abs/2022ApJS..259....3S/abstract)
- [Seon (2024), ApJ, 971, 184](https://ui.adsabs.harvard.edu/abs/2024ApJ...971..184S/abstract)
- [Yan et al. (2022), ApJ, 936, 177](https://ui.adsabs.harvard.edu/abs/2022ApJ...936..177Y/abstract) — Lyα radiative transfer in exoplanet atmospheres (WASP-52b)

## Documentation

For compilation instructions, input parameters, and usage examples, see [README_HOWTO.md](README_HOWTO.md).
Detailed usage instructions are in preparation.

---

Last updated: 2026-07-10 00:07
