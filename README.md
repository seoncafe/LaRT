# LaRT v2.0

LaRT v2.0 is a radiative transfer code for resonance-line scattering. In addition to Lyα, it is designed to handle resonance and fluorescence scattering of a variety of metallic lines.

This version adds support for adaptive mesh refinement (AMR) grids. RAMSES snapshots are not read directly by `LaRT.x`; instead, they are first converted to a generic AMR file (HDF5, gzipped FITS, or plain text) via a standalone converter (`convert_ramses_to_generic.{x,py}`). **The AMR implementation is currently under active testing and validation.**

## Status

- AMR support: octree grid with per-leaf physical quantities; testing and validation are ongoing.
- **per-cell core-skip**: `xcrit` is now recomputed at every scattering from the local cell's `voigt_a × rhokap × dl_face` (Smith+15 Eq.35). This replaces the volume-averaged `xcrit` for inhomogeneous (e.g. galaxy-scale) boxes where the old prescription gave `xcrit = 0` and photons trapped in dense regions never escaped. The old behavior is retained behind `par%core_skip_global = .true.` for benchmarking.
- Clump overlap handling added: file-loaded or internally generated overlapping clumps are handled with an event-based multi-component raytrace; enabled via `par%clump_allow_overlap = .true.` for internally generated populations.
- **HDF5 I/O support (output and input)** is now available alongside FITS for all output streams (spectrum, peel-off, sight-line tau, CALCJ/CALCP arrays) and all gridded inputs (density, temperature, velocity, emissivity, clump files). The default is HDF5 (`par%file_format = 'hdf5'`); set `par%file_format = 'fits'` to fall back to `.fits.gz`. Build with `make HDF5=1` (default) or `make HDF5=0` to skip the HDF5 link dependency. See `python/lart_io.py` for a format-agnostic Python reader/converter.
- The usage guide is still being prepared.

## References

For scientific background and related methodology, see:

- [Seon & Kim (2020), ApJS, 250, 9](https://ui.adsabs.harvard.edu/abs/2020ApJS..250....9S/abstract)
- [Seon et al. (2022), ApJS, 259, 3](https://ui.adsabs.harvard.edu/abs/2022ApJS..259....3S/abstract)
- [Seon (2024), ApJ, 971, 184](https://ui.adsabs.harvard.edu/abs/2024ApJ...971..184S/abstract)

## Documentation

Detailed usage instructions are in preparation.
