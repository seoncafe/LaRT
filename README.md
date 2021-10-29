LaRT (Lyman-alpha Radiative Transfer)

How to compile and run:

   1) The code is written in modern Fortran (Fortran 2003 or later is required).
      You need to install fortran/C compilers (for instance, Inten oneAPI Toolkit or GNU compilers).
      (gfortran v4 does not supprt fortran 2003.)
      (You need to modify Makefile if you want to use GNU gfortran.)

   2) You need to install CFITSIO library. (https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html)
      The version 3.410 does not work on macosx Sierra.

   3) You need to install an MPI library. (Either mpich or openmpi is okay.)
      a: Intel oneAPI HPC toolkit for Linux contains an MPI library, called openmpi.
      b: MPICH   -> www.mpich.org
         OPENMPI -> www.open-mpi.org
         In order to install openmpi on MacOSX, you may need to do "setenv TMPDIR /tmp" in tcsh shell or "export TMPDIR=/tmp" in bash.

   4) How to compile and run:
      unix> cd LaRT_v1.34 ; make
      unix> cd examples/sphere
      unix> mpirun -np 8 ../LaRT_peel_calcP.x t1tau4.in
      (Use the number of threads that your system has in the place of "8".)
      (A Quad-core CPU has 8 threads. Number of threads = 2 x number of cores.)
      Please refer to "run.sh" or "run_hybrid.sh" in each directory.

   5) examples are located under the directories, sphere, slab, etc.

   6) See "params_type" in define_v2.f90 for the default values of the input parameters.

Input parameters:

<<---------------------------------->>
<< note par% in front of parameters >>
<<---------------------------------->>

!-- Grid geometry.
  par%nx      = 101, number of cells in x-direction
  par%ny      = 101, number of cells in y-direction
  par%nz      = 101, number of cells in z-direction
  par%xmax    = 1.0, maximum positive x-value
  par%ymax    = 1.0, maximum positive y-value
  par%zmax    = 1.0, maximum positive z-value
  par%rmax    = -999, maximum radius of sphere in cartesian space.
                A spherical cloud will be constructed in cartesian space if par%rmax > 0
                (i.e., density = 0 for r = sqrt(x^2 + y^2 + z^2) > rmax if rmax > 0)

  The following options are helpful to use a less RAM memory by utilizing the symmetry of the system.
  par%xyz_symmetry = .false. (use only 1/8 part of a box if xyz_symmetry = .true.)
                     This will be set to be .false. when the peelingoff is performed.
  par%xy_periodic  = .false. (infinitely periodic slab in xy dimension)
  par%xy_symmetry  = .false. (use only 1/4 part of a box if xy_symmetry = .true.)

       box size = (-xmax, xmax) if xyz_symmetry = .false.
                = (0,     xmax) if xyz_symmetry = .true. and nx is even number.
                = (-dx/2, xmax) if xyz_symmetry = .true. and nx is odd  number, where dx = xmax/(nx-0.5).
                The system center is always at (x,y,z) = (0,0,0).
       If par%input_field is given, (nx,ny,nz) and (xmax,ymax,zmax) are obtained from header of density FITS file.

!-- Distance Unit
  par%distance_unit = '' (default) dimensionless unit for xmax, ymax, zmax of density.
                                   In this case, the density is determined by the optical depth.
                      'kpc' if distance unit for grid is killo-parsec.
                      'pc'  if distance unit for grid is parsec.
                      'au'  if distance unit for grid is AU.
  par%distance2cm   = 1.0 (default) distance unit in cm, or = 1 if the model is dimensionless.
                      Instead of using par%distance_unit, the distance unit in cm can be specified by par%distance2cm.

!-- Luminosity
  LaRT gives outputs for a unit luminosity.
  The output must be multiplied by a total luminosity of your system, that are generated inside the grid system
  or incident onto the grid system.

!-- Peelingoff & Observers
   If nxim and nyim are set to be non-zero integers, the peel-off process is performed to calculate spectral images that
   are supposed to observed by observer(s). If dxim and dyim are not given, they are automatically determined to cover
   the whole grid system. On the other hand, if dxim and dyim are given, the image size will be determined by nxim*dxim and
   nyim*dyim. The parameters dxim and dyim must be given in degree.

     par%nxim = 0
     par%nyim = 0
     par%dxim = nan64 (corresponds to 'CD1_1', as defined in the FITS standard.)
     par%dyim = nan64 (corresponds to 'CD2_2', as defined in the FITS standard.)
     par%distance = nan64 (distance between the system center and the observer)

   If par%distance is not given in the input file, the distance is assumed to be 100 x (maximum size of the grid system).
   The following parameters can be used to specify the location(s) of observer(s).
   If the following parameters are not given, the observer is assumed to be located at the direction of (obsx,obsy,obsz) = (0,0,1).
     par%obsx     = nan64
     par%obsy     = nan64
     par%obsz     = nan64
     par%distance = nan64
   or
     par%alpha    = nan64
     par%beta     = nan64
     par%distance = nan64

   The following parameter determines the orientation of the detector plane of the observer(s).
     par%gamma = nan64

   Multiple observers can be assumed by defining the above parameters as follows:
     par%alpha =  0.0 10.0 20.0
     par%beta  = 10.0 30.0 40.0

   The maximum number of observers are preset to be MAX_OBSERVERS = 99 in define_oo.f90.
   One may want to increase the maximum number of observers.
   But, a too large number of observers will requires a huge RAM.

!--
  par%use_master_slave  = .true.  (use master-worker algorithm for the load balancing between processes)
                                   if .false., then an equal number of photon packets is assigned to each processes.
  par%save_Jin          = .false. (save input spectrum if set to .true.)
  par%save_Jabs         = .false. (save the absorbed spectrum by dust grains if set to .true.)
  par%save_all          = .false. (save 3D J(nu,x,y,z) and Pa(x,y,z) arrays if set to .true.)
                                   By default, only 2D cylindrical or 1D spherical data is saved, depending on the geometry.
                                   par%geometry_JPa = 1, 2, or 3 can be set explicitly.
                                   par%save_all = .true. is equivalent to par%geometry_JPa = 3.
  par%use_reduced_wgt   = .false. (photons are destroyed when absorbed by dust grains)
                                  (One the other hand, if par%use_reduced_wgt = .true., the photon weight is reduced by an albedo.)
  par%use_cie_condition = .false. (Use Collisional Ionization Equilibrium Condition to calculate Neutral Hydrogen Density)

  par%out_file          = ''      (output FITS file name; The outputs will use the same name as the input file,
                                   if par%out_file is not given.)
  par%out_merge         = .false. (pre-existing output FITS file will be removed and new output file will be made if .false.)
                                  (On the other hand, if .true., new output results will be merged to pre-existing FITS file.)
                                  (Make sure that the input parameters are exactly the same as those of the previous run.)
                                   If you want to reduce noise, but still do not want to loose the previous results,
                                   then, set par%out_merge = .true.
  par%out_bitpix        = -32     (precision of the output FITS files)
                                  (single precision FITS files will be made if -32)
                                  (double precision FITS files will be made if -64)

  par%recoil     = .false. (the recoil effect is not considered)
  par%core_skip  = .false. (the core skipping algorithm is not used.)
                           (the core skipping algorithm will be updated later)
  par%no_photons = the number of photon packets to be used in calculation.
                   (this is a double precision floating number and is transformed to int64.)

!-- Optical Depth
  The following parameters are used to scale up or down the density field.
  Tese parameters can be used even for the case where a realistic density are given to scale it arbitrarily.

  par%taumax   = optical depth at the line center, measured from the center to outer boundary along the z-axis.
  par%tauhomo  = optical depth at the line center, measured from the center to outer boundary along the z-axis,
                 when the density is assumed to be constant.
  par%N_HImax  = column density of the HI gas measured from the center to outer boundary
  par%N_HIhomo = column density of the HI gas measured from the center to outer boundary,
                 when the density is assumed to be constant.
  par%N_HI     = par%N_HImax

!-- Density
  (1) If the following parameters are specified, then density is set to zero for r < rmin or r > rmax.
      These parameters can be used to simulate a spherical or shell geometry.
     par%rmin = -999 (default)
     par%rmax = -999 (default)
  (2) To simulate an exponential density profile, use the following parameter together with par%taumax, etc.
     par%density_rscale = -999 (default)

!-- Velocity
  There are three types of the velocity field implemented.
  par%velocity_type = 'hubble', 'constant_radial', or 'ssh'

  (1) par%velocity_type = 'hubble'
      V(r)     = Vexp * (r/r_max)
      par%Vexp = 0.0 (maximum velocity in km/s)
                   r     = radial vector from the cloud center
                   r_max = radius of the cloud
                   Vexp  > 0 : expanding medium (outflow)
                   Vexp  < 0 : contracting medium (infall)
  (2) par%velocity_type = 'contant_radial'
      V(r)     = Vexp * r / |r|
      par%Vexp = 0.0
  (3) par%velocity_type = 'ssh'  (Here, ssh stands for Song, Seon, Hwang 2020, ApJ)
      The radial velocity field is defined as follows:
      V(r) = Vpeak * r / rpeak                                             if |r| < rpeak
           = (Vpeak + DeltaV * (|r| - rpeak) / (rmax - rpeak)) * r / |r|   if |r| > rpeak
      par%Vpeak  = 0.0
      par%rpeak  = 0.0
      par%DeltaV = 0.0

!-- Density, Temperature, Velocity, and/or Emissivity input files.
  par%input_field = ... (base name for input density, temperature, and velocity fiels.)
                  = 'm1' if density file = 'm1.dens.fits.gz', temperature file = 'm1.temp.fits.gz',
                            and velocity file = 'm1.velo.fits.gz'
  par%dens_file   = a density FITS file nanme (for instance, 'm1_dens.fits.gz')
  par%temp_file   = a temperature FITS file name (e.g., 'm1_temp.fits.gz')
  par%velo_file   = a velocity FITS file name (e.g., 'm1_velo.fits.gz')
  par%emiss_file  = a emissivity file name (e.g., 'm1_emiss.fits.gz')
                  ---> density must be given in units of H number/cm^3.
                       (However, the density will be rescale if par%taumax, par%N_HImax, etc are given as an input.)
                  ---> temperature in K.
                  ---> velocity in units of km/s.
                  ---> emissivity in photons/s/cm^3.
  par%star_file   = a text file that describes information on individual stars.
                    The file must contain four columns (x, y, z, and luminosity).
                    The luminosity is not required to have an absolute physical unit.
                    Only their relative strengths are used.

!-- Source geometry and Source type
  par%source_geometry = 'point'   (point source at (x,y,z) = (0,0,0) or r = 0)
                      = 'uniform' (uniformly distributed source over the whole system)
                      = 'uniform_xy'
                      = 'gaussian'
                      = 'exponential'
                      = 'ssh' or 'sersic'
                      = 'star_file' -> Use this option to simulate multiple individual stars.
                      = 'diffuse_emissivity'
  par%source_zscale   = 0.0
                        (z-scale height for gaussian or exponential distribution)
  par%emiss_file      = ''
                        (The file name must be specified if par%source_geometry = 'diffuse_emissivity')

!-- Line profile
  par%comoving_sopurce = .true. (whether the photon source is comoving with the medium or not)
  par%spectral_type    = 'monochromatic' (delta function at the Lyman-alpha central frequency, meaning xfreq = 0)
                         'continuum'     (continnum spectrum over xfreq_min to xfreq_max)
                         'voigt'         (Voigt profile for the temperature at the photon emission location)
                         'voigt0'        (Voigt profile for the temperature of par%temperature0,
                                          independently of the temperaure of the photon emission location)
                         'gaussian'      (Gaussian line profile)
                                         The line width is specified by par%gassian_width_vel in units of km/s.

  The following parameters can be used to simulate an arbitrary line profile.
     par%line_prof_file = ''
     par%line_prof_file_type = 0 (if the file contains two columns, frequency and line strength)
                               1 (if the file contains wavelength and line strength).

!-- Frequency (wavelength or velocity) range
  par%xfreq0    = 0.0, initial photon frequency (defined by x parameter)
                  The initial photon frequency will be shifted by this quantity if it is set to a non-zero.
  par%nxfreq    = 121, Number of frequencies for the output spectral binning.
  par%xfreq_min = minimum frequency
  par%xfreq_max = maximum frequency
                  x = (nu-nu_0)/Doppler width (nu_0 = Lyman alpha frequency)
                  freqeuncy limits are automatically calculated when no input is given.

  If the following parameters are specified in the input file, the xfreq range is determined accordingly.
     par%nlambda      = 0
     par%lambda_min   = nan64
     par%lambda_max   = nan64
     par%nvelocity    = 0
     par%velocity_min = nan64
     par%velocity_max = nan64

  The frequency range is automatically determined by the mean temperature of the medium.
  But, this autumatically determined range can be inappropriate sometimes.
  If one wants to change the frequency or wavelength bins, use the above parameters.

!-- Dust
  par%DGR    = 0.0,     dust-to-gas ratio, dust abundance relative to Milky-Way dust-to-gas ratio.
  par%hgg    = 0.67617, asymmetry scattering phase factor of dust scattering (Heyney-Greenstein phase function)
  par%albedo = 0.32536, dust albedo (scattering cross-section/extinction cross-section)
    If par%use_stokes = .true., then hgg and albedo are not required to be explicitly specified.

!-- Polarization
  Set up the following parameters to simulation the Ly-alpha polarization.
  The scattering matrix files, calculated using the Weingartner-Draine dust model for the MilkyWay, are included in the "data" directory.
    par%scatt_mat_file = ''
    par%use_stokes     = .false.

!-- TIGRESS simulation: Use par%omega to simulate the shear effect for the TIGRESS simulation data.
  par%Omega       = 28.0 (Omega value for TIGRESS data in (km/s)/kpc)

!-- Atomic data for Lya,
  Do not edit the following parameters, unless you know what you are doing.
    par%f12 = 0.4126,     oscillator strength
    par%A21 = 6.265x10^8, Einstein A coeffiecient

