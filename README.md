# LaRT v2.0

LaRT v2.0 is a radiative transfer code for resonance-line scattering. In addition to Ly-alpha, it is designed to handle resonance and fluorescence scattering of a variety of metallic lines.

This version adds support for adaptive mesh refinement (AMR) grids, including both RAMSES-based AMR outputs and more general AMR grid data. **The AMR implementation is currently under active testing and validation, and its behavior should be checked carefully before being used as fully verified production software.**

## Status

- AMR support has been added; testing and validation are ongoing.
- Clump overlap handling added: file-loaded or internally generated overlapping clumps are handled with an event-based multi-component raytrace; enabled via `par%clump_allow_overlap = .true.` for internally generated populations.
- Both RAMSES AMR snapshots and a generic text/FITS/HDF5 AMR format are supported.
- **HDF5 I/O support (output and input)** is now available alongside FITS for all output streams (spectrum, peel-off, sight-line tau, CALCJ/CALCP arrays) and all gridded inputs (density, temperature, velocity, emissivity, clump files). The default is HDF5 (`par%file_format = 'hdf5'`); set `par%file_format = 'fits'` to fall back to `.fits.gz`. Build with `make HDF5=1` (default) or `make HDF5=0` to skip the HDF5 link dependency. See `HDF5_CHANGES.txt` and `python/lart_io.py` for the schema and a format-agnostic Python reader/converter.
- The usage guide is still being prepared.

## References

For scientific background and related methodology, see:

- [Seon & Kim (2020), ApJS, 250, 9](https://ui.adsabs.harvard.edu/abs/2020ApJS..250....9S/abstract)
- [Seon et al. (2022), ApJS, 259, 3](https://ui.adsabs.harvard.edu/abs/2022ApJS..259....3S/abstract)
- [Seon (2024), ApJ, 971, 184](https://ui.adsabs.harvard.edu/abs/2024ApJ...971..184S/abstract)

## Documentation

Detailed usage instructions are in preparation.
