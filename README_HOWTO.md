# LaRT User Guide

## How to Compile and Run

### Prerequisites

1. **Fortran/C compilers**: Fortran 2003 or later is required (e.g., Intel oneAPI Toolkit or GNU compilers).
   `gfortran` v4 does not support Fortran 2003. You may need to modify the Makefile for GNU `gfortran`.

2. **CFITSIO library**: <https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html>

3. **HDF5 library** (optional): For HDF5 output support, install HDF5 with Fortran bindings compiled
   with the same compiler as LaRT (e.g., ifx-built HDF5 for Intel).
   Build with `make HDF5=1` (the default). To skip HDF5, use `make HDF5=0`.

4. **MPI library**: Either MPICH or OpenMPI.
   - Intel oneAPI HPC toolkit for Linux contains an MPI library.
   - MPICH: <https://www.mpich.org>
   - OpenMPI: <https://www.open-mpi.org>

### Build Options

Before compiling, edit the Makefile to set the preprocessor options. For usual purposes, set all to 0.

| Flag | Description |
|------|-------------|
| `CALCPnew=1` | Calculate "scattering rate" using the faster method (Seon & Kim 2020) |
| `CALCJ=1` | Calculate "mean intensity" within the medium (Seon & Kim 2020) |
| `CALCP=1` | Calculate "scattering rate" using the slower method (Seon & Kim 2020) |
| `FINE_STRUCTURE=1` | Consider fine-structure levels of the n=2 state (much slower) |

`CALCPnew=1`, `CALCJ=1`, and/or `CALCP=1` will require a large amount of RAM.

### Compile and Run

```bash
make
cd examples/sphere
mpirun -np 8 ../../LaRT.x t1tau4.in
```

Use the number of CPU threads your system has in place of `8`.
See `run.sh` in each example directory.

### Reference

See `params_type` in `define.f90` for the default values of all input parameters.

---

## Input Parameters

> **Note:** All parameters require the `par%` prefix in the namelist input file.

### Grid Geometry

| Parameter | Default | Description |
|-----------|---------|-------------|
| `nx`, `ny`, `nz` | 1, 1, 11 | Number of cells in x, y, z |
| `nr` | -999 | Shorthand: sets `nx = ny = nz = nr` |
| `xmax`, `ymax`, `zmax` | 1.0 | Maximum positive extent in each direction |
| `rmax` | -999 | Sphere radius; density = 0 for r > rmax when rmax > 0 |
| `rmin` | -999 | Inner radius for shell geometry |
| `geometry` | `''` | `'sphere'`, `'rectangle'`/`'box'`, `'cylinder'`, `'plane_atmosphere'`, `'spherical_atmosphere'` |

**Symmetry options** (reduce RAM usage):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `xyz_symmetry` | `.false.` | Use only 1/8 of box. Set to `.false.` when peel-off is performed. |
| `xy_symmetry` | `.false.` | Use only 1/4 of box |
| `xy_periodic` | `.false.` | Infinitely periodic slab in xy |

Box size:
- `(-xmax, xmax)` if `xyz_symmetry = .false.`
- `(0, xmax)` if `xyz_symmetry = .true.` and `nx` is even
- `(-dx/2, xmax)` if `xyz_symmetry = .true.` and `nx` is odd, where `dx = xmax/(nx-0.5)`

The system center is always at (0, 0, 0).
If `par%input_field` is given, `(nx,ny,nz)` and `(xmax,ymax,zmax)` are read from the density file header.

### Grid Type

| Parameter | Default | Description |
|-----------|---------|-------------|
| `grid_type` | `'car'` | `'car'` (Cartesian) or `'amr'` (AMR octree; same as `use_amr_grid = .true.`) |
| `use_clump_medium` | `.false.` | Clumpy medium mode (see [Clumpy Medium](#clumpy-medium)) |

### Distance Unit

| Parameter | Default | Description |
|-----------|---------|-------------|
| `distance_unit` | `''` | `'kpc'`, `'pc'`, `'Mpc'`, `'au'`, or `''` (dimensionless code units) |
| `distance2cm` | 1.0 | Distance unit in cm (alternative to `distance_unit`) |

When `distance_unit = ''`, the model is dimensionless and the density is determined by
`par%taumax`, `par%tauhomo`, `par%N_HImax`, or `par%N_HIhomo`.

### Luminosity

LaRT gives outputs for a unit luminosity.
The output must be multiplied by the total luminosity of your system.

### Peel-Off Observers

If `nxim` and `nyim` are set to non-zero integers, the peel-off process is performed to
calculate spectral images. If `dxim` and `dyim` are not given, they are automatically
determined to cover the whole grid system.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `save_peeloff` | `.false.` | Enable peel-off observer images |
| `save_peeloff_2D` | `.false.` | Save 2D spatial map (integrated over frequency) |
| `save_peeloff_3D` | `.true.` | Save 3D data cube (x, y, frequency) |
| `nxim`, `nyim` | 0 | Image pixel counts |
| `dxim`, `dyim` | nan | Pixel size in degrees (`CD1_1`, `CD2_2` in FITS) |
| `distance` | nan | Distance from system center to observer (same units as grid) |
| `save_radial_profile` | `.false.` | Save radial intensity profile from peel-off |

**Observer direction** (choose one set):

```fortran
par%obsx = ..., par%obsy = ..., par%obsz = ..., par%distance = ...
! or
par%alpha = ..., par%beta = ..., par%gamma = ..., par%distance = ...
```

If not given, the observer is at (0, 0, 1). If `par%distance` is not given,
it defaults to 100 times the maximum grid size.

**Multiple observers:**

```fortran
par%alpha =  0.0 10.0 20.0
par%beta  = 10.0 30.0 40.0
```

The maximum number of observers is `MAX_OBSERVERS = 99` (defined in `define.f90`).
A large number of observers will require a huge amount of RAM.
The peel-off raytrace is performed for every observer at each scattering event,
so computation time also increases linearly with the number of observers.

### HEALPix Inside Observer (All-Sky Maps)

Setting `par%nside > 0` enables HEALPix all-sky observer mode.
`par%observer_located_inside` is automatically set to `.true.`,
switching the entire output chain to HEALPix all-sky maps.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `nside` | 0 | HEALPix resolution parameter (0 = disabled) |

Stokes polarization is not yet supported for the HEALPix path.

### Angular Spectrum (Jmu)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `save_Jmu` | `.false.` | Save escaped spectrum binned by mu = cos(theta_z) |
| `nmu` | 11 | Number of mu bins (odd default places mu=0 at bin center) |

`Jmu(nxfreq, nmu)` is written as a 2D image extension with full WCS.
The mu range is [0, +1] if `xyz_symmetry = .true.` (folded to +z hemisphere)
and [-1, +1] otherwise. Normalization: `mean(Jmu, axis=mu) = Jout`.

### Simulation Control

| Parameter | Default | Description |
|-----------|---------|-------------|
| `no_photons` | 1e5 | Number of photon packets (double precision, converted to int64) |
| `use_master_slave` | `.true.` | Master-worker algorithm for load balancing. If `.false.`, equal photon counts per process. |

### Output Control

| Parameter | Default | Description |
|-----------|---------|-------------|
| `save_Jin` | `.true.` | Save input spectrum |
| `save_Jabs` | `.true.` | Save spectrum absorbed by dust |
| `save_all` | `.false.` | Save 3D J(nu,x,y,z) and Pa(x,y,z). Equivalent to `geometry_JPa = 3`. |
| `geometry_JPa` | auto | 1 (spherical 1D), 2 (cylindrical 2D), or 3 (full 3D) |
| `save_backup` | `.false.` | Save backup copies before `out_merge` |
| `save_direc0` | `.false.` | Save direct light component |
| `save_sightline_tau` | `.false.` | Save sightline optical depth map |
| `save_dust_scattered` | `.false.` | Save dust-scattered component separately |
| `save_input_grid` | `.false.` | Save input grid data to the output file |
| `out_file` | `''` | Output file name (defaults to input file base name) |
| `out_merge` | `.false.` | Merge new results into pre-existing output. Input parameters must match. |
| `out_bitpix` | 0 | Output precision: 0 = auto, -32 = float32, -64 = float64 |
| `file_format` | `'hdf5'` | Output format: `'hdf5'` or `'fits'`. Build with `make HDF5=0` to disable HDF5. |

| Parameter | Default | Description |
|-----------|---------|-------------|
| `use_reduced_wgt` | `.false.` | If `.true.`, photon weight is reduced by albedo instead of destroying absorbed photons |
| `use_cie_condition` | `.false.` | Use CIE for neutral fraction. In AMR mode, use `ionization_model` instead. |

### Scattering Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `recoil` | `.false.` | Include the recoil effect |
| `core_skip` | `.false.` | Per-cell core-skipping: `xcrit = (voigt_a * rhokap * dl_face)^(1/3) / 5` (Smith+15 Eq.35). Required for inhomogeneous boxes. |
| `core_skip_global` | `.false.` | Fall back to old volume-averaged xcrit (benchmarking only) |

### Optical Depth

The following parameters scale the density field.
They can be used even when a realistic density is given.

| Parameter | Description |
|-----------|-------------|
| `taumax` | Optical depth at line center, from center to outer boundary along z-axis |
| `tauhomo` | Same, assuming constant (homogeneous) density |
| `N_HImax` | HI column density from center to outer boundary |
| `N_HIhomo` | Same, assuming constant density |
| `N_HI` | Alias for `N_HImax` (deprecated) |
| `N_gasmax` | Ion column density [cm^-2] (for metal lines with `ion_model='none'`) |

### Density

| Parameter | Default | Description |
|-----------|---------|-------------|
| `rmin` | -999 | Inner radius (density = 0 for r < rmin) |
| `rmax` | -999 | Outer radius (density = 0 for r > rmax) |
| `density_rscale` | -999 | Scale length for exponential density profile |
| `density_alpha` | 0.0 | Power-law index: n(r) = n0 * (rmax/r)^density_alpha. 0 = uniform, 2 = isothermal. |
| `cone_opening` | 0.0 | Biconical outflow half-opening angle [degrees] along z-axis. 0 = full sphere. When > 0, density is set to zero outside the cone (both hemispheres). Works with Cartesian, AMR, and clumpy-medium modes. |

### Velocity

| `velocity_type` | Formula |
|-----------------|---------|
| `'hubble'` | V(r) = Vexp * (r / r_max) |
| `'constant_radial'` | V(r) = Vexp * r / \|r\| |
| `'power_law'` | V(r) = Vexp * (r / rmax)^velocity_alpha |
| `'linear_decelerate'` | V(r) = Vexp * (rmax - r) / (rmax - rmin); V=Vexp at rmin, V=0 at rmax |
| `'ssh'` | Song, Seon, Hwang (2020); two-zone radial profile |
| `'rotating_solid_body'` | V = Vrot * r_cyl / rmax |
| `'rotating_galaxy_halo'` | Solid-body core (r < rinner) + flat rotation curve |
| `'parallel_velocity'` | V = (Vx, Vy, Vz) constant bulk velocity |

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `Vexp` | 0.0 | Maximum velocity [km/s]; >0 = outflow, <0 = infall |
| `velocity_alpha` | 1.0 | Power-law exponent for `'power_law'` velocity type |
| `Vpeak` | 0.0 | Peak velocity for `'ssh'` [km/s] |
| `rpeak` | 0.0 | Peak radius for `'ssh'` |
| `DeltaV` | 0.0 | Velocity difference for `'ssh'` [km/s] |
| `Vrot` | 0.0 | Rotation velocity [km/s] |
| `rinner` | 0.0 | Inner radius for `'rotating_galaxy_halo'` |
| `Vx`, `Vy`, `Vz` | 0.0 | Bulk velocity components for `'parallel_velocity'` [km/s] |

### Input Data Files

| Parameter | Description |
|-----------|-------------|
| `input_field` | Base name for density/temperature/velocity files (e.g., `'m1'` for `m1.dens.fits.gz`, `m1.temp.fits.gz`, `m1.velo.fits.gz`) |
| `dens_file` | Density file (H number density in cm^-3) |
| `temp_file` | Temperature file [K] |
| `velo_file` | Velocity file [km/s] |
| `emiss_file` | Emissivity file [photons/s/cm^3] |
| `star_file` | Text file with columns (x, y, z, luminosity); only relative luminosities are used |

Input files can be FITS (`.fits`, `.fits.gz`) or HDF5 (`.h5`, `.hdf5`).
The density will be rescaled if `par%taumax`, `par%N_HImax`, etc. are given.

### Source Geometry

| Parameter | Default | Description |
|-----------|---------|-------------|
| `source_geometry` | `'point'` | `'point'`, `'uniform'`, `'uniform_xy'`, `'gaussian'`, `'exponential'`, `'ssh'`/`'sersic'`, `'star_file'`, `'diffuse_emissivity'` |
| `source_zscale` | 0.0 | z-scale height for Gaussian or exponential source |
| `xs_point`, `ys_point`, `zs_point` | 0.0 | Point source location |

### Line Profile

| Parameter | Default | Description |
|-----------|---------|-------------|
| `comoving_source` | `.true.` | Photon source is comoving with the medium |
| `spectral_type` | `'voigt'` | `'voigt'`, `'voigt0'`, `'monochromatic'`, `'continuum'`, `'gaussian'`, `'continuum+gaussian'` |
| `gaussian_sigma_vel` | 12.84 | Gaussian sigma [km/s] for `spectral_type = 'gaussian'` |
| `gaussian_FWHM_vel` | -999 | Gaussian FWHM [km/s]; overrides sigma if > 0. Used by `'gaussian'` and `'continuum+gaussian'`. |
| `EW_line` | 0.0 | Equivalent width [Angstrom] for `'continuum+gaussian'` |
| `continuum_normalize` | `.true.` | Normalize output spectrum so continuum level = 1 |
| `line_prof_file` | `''` | Custom line profile from file |
| `line_prof_file_type` | 0 | 0 = (frequency, strength), 1 = (wavelength, strength) |

### Frequency / Wavelength / Velocity Range

| Parameter | Default | Description |
|-----------|---------|-------------|
| `xfreq0` | 0.0 | Initial photon frequency offset (x parameter) |
| `nxfreq` | 121 | Number of frequency bins for output spectrum |
| `xfreq_min`, `xfreq_max` | auto | Frequency range (auto-determined from mean temperature if not given) |
| `nvelocity` | 0 | Number of velocity bins (overrides xfreq if > 0) |
| `velocity_min`, `velocity_max` | nan | Velocity range [km/s] |
| `nwavelength` | 0 | Number of wavelength bins (overrides xfreq if > 0) |
| `wavelength_min`, `wavelength_max` | nan | Wavelength range |

### Dust

| Parameter | Default | Description |
|-----------|---------|-------------|
| `DGR` | 0.0 | Dust-to-gas ratio relative to Milky Way (0 = no dust) |
| `hgg` | 0.6761 | Henyey-Greenstein asymmetry parameter |
| `albedo` | 0.3253 | Dust single-scattering albedo |
| `cext_dust` | 1.6059e-21 | Dust extinction cross-section per H atom [cm^2] |

If `use_stokes = .true.`, then `hgg` and `albedo` are obtained from the scattering matrix file.

### Polarization

| Parameter | Default | Description |
|-----------|---------|-------------|
| `scatt_mat_file` | `''` | Mueller matrix file (e.g., Weingartner-Draine dust model) |
| `use_stokes` | `.false.` | Enable full Stokes polarization (I, Q, U, V) |

### Shearing Box (TIGRESS)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `Omega` | 0.0 | Angular velocity [(km/s)/kpc] |
| `q` | 1.0 | Shear parameter |

### Line ID and Per-Line Options

| `par%line_id` | Ion | Wavelength | Type |
|---------------|-----|------------|------|
| `'ly_alpha'` (default) | H I | 1216 | singlet / doublet |
| `'ly_alpha_HD'` | H I + D I | 1216 + 1215 | dual resonance |
| `'CII_1334'` | C II | 1334 (+1336 fluorescence) | resonance + fluorescence |
| `'CIV_1548'` | C IV | 1548, 1551 | doublet |
| `'NV_1239'` | N V | 1239, 1243 | doublet |
| `'OVI_1032'` | O VI | 1032, 1038 | doublet |
| `'NaI_D'` | Na I | 5892, 5898 | doublet |
| `'CaII_HK'` | Ca II | 3935, 3970 | doublet |
| `'MgII_2796'` | Mg II | 2796, 2804 | doublet |
| `'SiIV_1394'` | Si IV | 1394, 1403 | doublet |
| `'AlII_1671'` | Al II | 1671 | singlet |
| `'SiII_1260'` | Si II | 1260 | resonance + fluorescence |
| `'SiII_1193'` / `'SiII_1190'` | Si II | 1193, 1190 | 2 resonances + fluorescence |
| `'SiII_1304'` | Si II | 1304 | resonance + fluorescence |
| `'SiII_1527'` | Si II | 1527 | resonance + fluorescence |
| `'FeII_2250'` | Fe II | 2250 | resonance + fluorescence |
| `'FeII_2261'` | Fe II | 2261 | resonance + fluorescence |
| `'FeII_UV3'` / `'FeII_2344'` | Fe II | 2344 | resonance + 2 fluorescence |
| `'FeII_UV2'` / `'FeII_2383'` | Fe II | 2383 | 2 resonances + fluorescence |
| `'FeII_UV1'` / `'FeII_2600'` | Fe II | 2600 | 2 resonances + fluorescence |
| `'HeI_10833'` | He I | 10833 | triplet |

**Per-line options:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `fine_structure` | `.false.` | For `'ly_alpha'`: include H I 2P1/2-2P3/2 splitting (doublet) |
| `include_deuterium` | `.false.` | With `'ly_alpha'`, promotes to `'ly_alpha_HD'` |
| `D_to_H_ratio` | 1.5e-5 | Deuterium-to-hydrogen number ratio (cosmic primordial) |
| `HeI_coherent` | `.false.` | For `'HeI_10833'`: use frequency-dependent (E1, E2, E3) from the Real-Phi polynomial form instead of per-component incoherent values. See `examples/HeI_coherent_test/`. |

**Atomic data** (do not edit unless you know what you are doing):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `f12` | 0.4126 | Oscillator strength |
| `A21` | 6.265e8 | Einstein A coefficient |

---

## AMR Mode

Activated by `par%grid_type = 'amr'` or `par%use_amr_grid = .true.`.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `amr_type` | `'generic'` | Only `'generic'` is supported. Direct RAMSES reading was moved to a standalone converter. |
| `amr_file` | `''` | Path to generic AMR file (`.h5`, `.hdf5`, `.fits.gz`, or `.dat`) |
| `distance_unit` | `''` | Unit of position columns: `'kpc'`, `'pc'`, `'Mpc'`, `'au'`, or `''` (dimensionless code units; density determined by `taumax`, `tauhomo`, `N_HImax`, or `N_HIhomo`) |

### Generic AMR File Format

**Mandatory columns (9):**

| Column | Unit | Description |
|--------|------|-------------|
| x, y, z | distance_unit | Leaf-cell center positions |
| level | -- | AMR refinement level |
| nH | cm^-3 | Hydrogen number density |
| T | K | Temperature |
| vx, vy, vz | km/s | Velocity |

**Optional columns** (detected by name in HDF5/FITS):

| Column | Unit | Description |
|--------|------|-------------|
| metallicity | mass fraction | Z from simulation |
| xHI | dimensionless | Neutral hydrogen fraction |
| n_e | cm^-3 | Electron density |
| n_ion | cm^-3 | Scattering ion density (for metal lines) |
| emissivity | cm^-3 s^-1 | Lya emissivity rate |
| ndust | cm^-3 | Dust pseudo-number density |

When a column is present, LaRT uses it directly.
When absent, LaRT computes the quantity using the model selected below.

### AMR Physics Models

| Parameter | Default | Values |
|-----------|---------|--------|
| `ionization_model` | `'cie_formula'` | `'cie_formula'`, `'cie_table'` (Voronov + Verner), `'full_neutral'`, `'from_file'` |
| `dust_model` | `'global_dgr'` | `'global_dgr'`, `'laursen09'` (per-cell: ndust = Z/Z_ref * (nHI + f_ion * nHII)), `'from_file'` |
| `emissivity_model` | `'none'` | `'none'`, `'caseB'` (Case B recomb + collisional), `'from_file'` |
| `ion_model` | `'none'` | `'none'`, `'solar_cie'` (Asplund+09 * Gnat-Sternberg CIE), `'from_file'` |
| `metallicity_global` | -1.0 | Global Z for `'laursen09'`/`'solar_cie'` when no metallicity column (negative = unset) |
| `Z_ref` | 0.0134 | Reference solar metallicity (Asplund+09) |
| `f_ion_dust` | 0.01 | Dust survival fraction in ionized gas (Laursen+09) |

### Coordinate Convention

The file is expected in centered convention `[-boxlen/2, +boxlen/2]`
(i.e., `ORIGINX = ORIGINY = ORIGINZ = -boxlen/2`).
Corner-based files (`ORIGINX = 0`) are still readable;
LaRT will warn if the data does not cover the box center.

---

## Clumpy Medium

Activated by `par%use_clump_medium = .true.`.

### Basic Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `rmax` | 1.0 | Outer sphere radius (code units) |
| `clump_radius` | -1.0 | Individual clump radius |
| `clump_f_cov` | -1.0 | Covering factor (specify one of: `f_cov`, `f_vol`, `clump_N_clumps`) |
| `clump_f_vol` | -1.0 | Volume filling factor |
| `clump_N_clumps` | -1.0 | Number of clumps |
| `clump_tau0` | -1.0 | Line-center tau from clump center to surface |
| `clump_NHI` | -1.0 | HI column density through clump center |
| `clump_nH` | -1.0 | Hydrogen number density inside clumps |
| `clump_temperature` | -1.0 | Clump temperature [K] (defaults to `par%temperature`) |
| `clump_sigma_v` | 0.0 | Gaussian sigma of random bulk velocity [km/s] |
| `clump_fully_inside` | `.true.` | All clumps fully inside the sphere |
| `clump_allow_overlap` | `.false.` | Allow overlapping clumps (overlap-aware raytrace) |
| `save_clump_info` | `.false.` | Write clump positions/velocities to output |
| `clump_input_file` | `''` | Load clump positions from external file (`.h5`, `.fits.gz`, or `.dat`). Generate with `make make_clumps && ./make_clumps.x <in>` (Fortran) or `python python/make_clumps.py ... -o file.h5` (Python). |

### Spatially-Varying Clump Profiles

| Parameter | Default | Description |
|-----------|---------|-------------|
| `clump_radius_profile` | `'constant'` | Clump radius profile r_cl(r) |
| `clump_density_profile` | `'constant'` | n_H(r) inside clumps |
| `clump_number_profile` | `'constant'` | n_cl(r) spatial number density |
| `clump_radius_alpha`, `clump_radius_r0` | 0.0 | Profile shape parameters |
| `clump_density_alpha`, `clump_density_r0` | 0.0 | Profile shape parameters |
| `clump_number_alpha`, `clump_number_r0` | 0.0 | Profile shape parameters |
| `clump_radius_min` | -1.0 | Minimum clump radius |
| `clump_radius_max_in` | -1.0 | Maximum clump radius (input) |
| `clump_profile_file` | `''` | Tabulated profile file |

**Supported velocity types in clump mode:**
`'hubble'`, `'constant_radial'`, `'power_law'`, `'linear_decelerate'`,
`'parallel_velocity'`, `'ssh'`, `'rotating_solid_body'`, `'rotating_galaxy_halo'`

**Supported density/number profiles in clump mode:**
`'constant'`, `'power_law'` (alias `'powerlaw'`), `'gaussian'`, `'exponential'`, `'file'`

---

## Converters and Tools

### RAMSES Converters (`python/AMR_grid/`)

| Tool | Description |
|------|-------------|
| `convert_ramses_to_generic.x` (Fortran) | RAMSES snapshot to generic AMR file. `--compute-physics` adds xHI, n_e, Case B emissivity, Laursen+09 dust columns. |
| `convert_ramses_to_generic.py` (Python) | Same as above, Python version. |
| `extract_amr_subset.py` | Cut a cubic sub-region from a generic AMR file (centered convention). |
| `extract_amr_region.py` | Like `extract_amr_subset`, but keeps original BOXLEN. Streams HDF5 for large files. |
| `recenter_amr.py` | Shift coordinates so box is centered on origin (`ORIGINX = -BOXLEN/2`). |

**Example workflow:**

```bash
# 1) RAMSES -> generic
python convert_ramses_to_generic.py <repo> <snapnum> \
       -o sim.h5 --output-unit kpc --compute-physics

# 2) Extract a centered cube
python extract_amr_subset.py sim.h5 \
       --center 15 15 15 --size 30 \
       -o sim_subset.h5

# 3) Run LaRT
mpirun -np N LaRT.x input.in   # par%amr_file = 'sim_subset.h5'
```

### Python Analysis Tools (`python/`)

| Tool | Description |
|------|-------------|
| `read_lart.py` | High-level reader for LaRT outputs (`.fits.gz` or `.h5`). `read_lart('input.in')` returns a `LaRTOutput` with `summary()`, `plot_spectrum()`, `plot_peeling_map()`, `plot_peeling_spectrum()`, `plot_peel_jmu_compare()`, `plot_velocity_moment_map()`, `plot_clump_slice()`, `plot_jmu()`, etc. |
| `lart_io.py` | Format-agnostic reader/converter. CLI: `python lart_io.py info <file>` or `python lart_io.py convert <src> <dst>`. |
| `AMR_grid/AMR_grid.py` | In-memory `AMRGrid` builder + `slice_plot()`. Quantities: `nH`, `T`, `vx/vy/vz`, `v_mag`, `emissivity`, `xHI`, `n_e`, `n_ion`, `metallicity`, `ndust`, or a callable. |
