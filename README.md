# LaRT (Ly&alpha; Radiative Transfer)

## What is LaRT?
LaRT is a three-dimensional Monte Carlo Ly&alpha; radiative transfer code, which is developed to study the Ly&alpha; RT. LaRT is capable of predicting (1) the emergent Ly&alpha; spectrum and surface brightness profile, (2) the Wouthuysen-Field effect, and (3) the polarization of Ly&alpha; radiation. The code is capable of treating arbitrary geometries, density distributions and source distributions. LaRT uses the "peeling-off" (next event estimation) technique to obtain high signal-to-noise spectro=polarimetric images in a detector plane. LaRT can place an observer in an arbitrary location and make a detector plane to have an arbitrary orientation in the sky. The current version of LaRT uses a Cartesian grid to model the density distribution of hydrogen gas. LaRT has been benchmarked for a number of standard cases.

### Versions
LaRT is provided in two versions.
  - LaRT_v1.34 : pure MPI version
  - LaRT_v1.34_hybrid : MPI + openmp hybrid version (The hybrid versioin uses a less amount of RAM than the MPI version.)

### References
  - [Seon & Kim (2020 ApJS, 250, 9)](https://ui.adsabs.harvard.edu/abs/2020ApJS..250....9S/abstract)
  - [Seon, Song, & Chang (2021, ApJS, submitted)](https://ui.adsabs.harvard.edu/abs/2021arXiv210515062S%2F/abstract)

## How to compile and run:

1. The code is written in modern Fortran (Fortran 2003 or later).\
      You need to install fortran/C compilers (e.g., Intel oneAPI Toolkit or GNU compilers).\
      Note that gfortran v4 does not supprt fortran 2003.\
      You may need to modify Makefile if you want to use GNU gfortran.
2. You need to install CFITSIO library.\
(https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html)
3. You need to install an MPI library. (Either mpich or openmpi is okay.)
  - Intel oneAPI HPC toolkit for Linux contains an MPI library, called openmpi.
  - MPICH   -> www.mpich.org
  - OPENMPI -> www.open-mpi.org \
   In order to install openmpi on MacOSX, you may need to do "setenv TMPDIR /tmp" in tcsh shell or "export TMPDIR=/tmp" in bash.
4. How to compile and run:
  - unix> cd LaRT_v1.34 ; make
  - unix> cd examples/sphere
  - unix> mpirun -np 8 ../LaRT_calcP.x t1tau4.in \
      Use the number of threads that your system has in the place of "8". \
      Note that a Quad-core CPU has 8 threads. \
      Number of threads = 2 x number of cores. \
      Please refer to "run.sh" or "run_hybrid.sh" in each directory.
5. See "params_type" in define_v2.f90 for the default values of the input parameters.

### Setup a model & Examples:
  - Read README_HOWTO in each directory
  - See the example directory in each version. Examples are located under the directories, sphere, slab, etc.
  - Feel free to email me if you have any questions.
