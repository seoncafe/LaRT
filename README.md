**LaRT (Ly&alpha; Radiative Transfer)**

**What is LaRT?**

LaRT is a three-dimensional Monte Carlo Ly&alpha; radiative transfer code, which is developed to study the Ly&alpha; RT, the Wouthuysen-Field effect, and the polarization of Ly&alpha; radiation.

**How to compile and run:**

1. The code is written in modern Fortran (Fortran 2003 or later).\
      You need to install fortran/C compilers (e.g., Inten oneAPI Toolkit or GNU compilers).\
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
5. examples are located under the directories, sphere, slab, etc.
6. See "params_type" in define_v2.f90 for the default values of the input parameters.

**Versions:**
LaRT is provided in two versions.
  - LaRT_v1.34 : pure MPI version
  - LaRT_v1.34_hybrid : MPI + openmp hybrid version (The hybrid versioin uses a less amount of RAM than the MPI version.)

**How to setup a model:**
  - Read README_HOWTO in each directory
  - Feel free to email me if you have any questions.

**Examples:**
  - See the example directory in each version.
