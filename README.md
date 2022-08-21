 # LaRT (Ly&alpha; Radiative Transfer)

## What is LaRT?
LaRT is a three-dimensional Monte Carlo Ly&alpha; radiative transfer code, developed to study the Ly&alpha; radiative transfer. LaRT is capable of predicting (1) the Ly&alpha; spectrum and surface brightness profile, (2) the Wouthuysen-Field effect, and (3) the polarization of Ly&alpha; radiation. The code is capable of treating arbitrary geometries, density distributions, and source distributions. LaRT uses the "peeling-off" (next event estimation) technique to obtain high signal-to-noise spectro-polarimetric images in a detector plane. LaRT is superb, compared to the preexisting codes, in that it uses a smoothly and seamlessly varying phase function as frequency changes to deal with the polarization of Ly&alpha; radiation. LaRT can place an observer in an arbitrary location and make a detector plane have an arbitrary orientation in the sky. The current version of LaRT uses a Cartesian grid to model the density distribution of hydrogen gas. LaRT has been benchmarked for a number of standard cases.

### Versions
LaRT is provided in two versions.
  - LaRT_v1.34b : pure MPI version
  - LaRT_v1.34b_hybrid : MPI + openmp hybrid version (The hybrid version uses a smaller amount of RAM than the MPI version.)

### References
If you use LaRT, please acknowledge the following two papers.
  - [Seon & Kim (2020 ApJS, 250, 9)](https://ui.adsabs.harvard.edu/abs/2020ApJS..250....9S/abstract)
  - [Seon, Song, & Chang (2022, ApJS, 259, 3)](https://ui.adsabs.harvard.edu/abs/2022ApJS..259....3S/abstract)

## How to compile and run:

1. The code is written in modern Fortran (Fortran 2003 or later).\
      You need to install fortran/C compilers (e.g., [Intel oneAPI HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#hpc-kit) or GNU compilers).\
      Note that gfortran v4 does not supprt fortran 2003.\
      You may need to modify Makefile if you want to use GNU gfortran.

2. You need to install CFITSIO library.\
(https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html)

3. You need to install an MPI library. (Either mpich or openmpi is okay.)
  - Intel oneAPI HPC toolkit for Linux contains an MPI library, called intel mpi.
  - MPICH   -> www.mpich.org
  - OPENMPI -> www.open-mpi.org \
   In order to install openmpi on MacOSX, you may need to do "setenv TMPDIR /tmp" in tcsh shell or "export TMPDIR=/tmp" in bash.

4. Edit Makefile
  - unix> cd LaRT_v1.34b
  - Before you compile the code, edit Makefile to set the four preprocessor options (CALCPnew, CALCJ, CALCP, FINE_STRUCTURE).
   For the usual purpose, it would be good to set all four options to 0.

   |option | explanation |
   |---------------------------|---------------------------------------|
   |CALCPnew       = 1| to calculate the "scattering rate" using the faster second method (Seon & Kim 2020).|
   |CALCJ          = 1| to calculate the "radiation field strength (mean intensity)" within the medium (Seon & Kim 2020).|
   |CALCP          = 1| to calculate the "scattering rate" using the slower first method (Seon & Kim 2020).|
   |FINE_STRUCTURE = 1| to consider the fine structure levels of the n = 2 state.|

   "CALCPnew = 1" "CALCJ = 1" and/or "CALCP = 1" will reauire a large amount of RAM memories.\
   "FINE STRUCTURE = 1" will make the code much slower.

5. Compile and run:
  - unix> make
  - unix> cd examples/sphere
  - unix> mpirun -np 8 ../../LaRT.x t1tau4.in \
      Use the number of threads that your system has in the place of "8". \
      Note that a Quad-core CPU has 8 threads. \
      Number of threads = 2 x number of cores. \
      Please refer to "run.sh" or "run_hybrid.sh" in each directory.

6. See "params_type" in define_v2.f90 for the default values of the input parameters.

### Setup a model & Examples:
  - Read README_HOWTO in each directory
  - See the example directory in each version. Examples are located under the directories, sphere, slab, etc.
  - Feel free to email me if you have any questions.


## Author
Dr./Prof. Kwang-Il Seon \
[https://seoncafe.github.io](https://seoncafe.github.io) \
[KASI](http://www.kasi.re.kr)/[UST](http://ust.ac.kr)
