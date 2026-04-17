module define
  use, intrinsic :: iso_fortran_env, only: real32, real64, real128, int32, int64
public

!! kind parameter for single and double precision (6 and 15 significant decimal digits)
!! precision of 32-, 64-, and 128-bit reals
!  integer, parameter :: sp = selected_real_kind(6, 37)
!  integer, parameter :: dp = selected_real_kind(15, 307)
!  integer, parameter :: qp = selected_real_kind(33, 4931)
!! precision of 32-, 64- bit integers
!  integer, parameter :: i4b = selected_int_kind(r=9)
!  integer, parameter :: i8b = selected_int_kind(r=18)
  integer, parameter :: sp  = real32
  integer, parameter :: dp  = real64
  integer, parameter :: qp  = real128
  integer, parameter :: i4b = int32
  integer, parameter :: i8b = int64

! working precision
! Note that (1) iwp (integer working precision) applies only to seed of KISS random number generator.
!           (2) nphotons is always 64-bit integer.
!           (3) other integer variables are always 32-bit integers.
!  integer, parameter :: wp  = sp
!  integer, parameter :: iwp = i4b
  integer, parameter :: wp  = dp
  integer, parameter :: iwp = i8b

! tinest = the smallest positive number
! eps    = the least positive number
!          that added to 1 returns a number that is greater than 1
! hugest = the hugest positive number
  real(kind=wp), parameter :: tinest = tiny(0.0_wp)
  real(kind=wp), parameter :: eps    = epsilon(0.0_wp)
  real(kind=sp), parameter :: eps_sp = epsilon(0.0_sp)
  real(kind=dp), parameter :: eps_dp = epsilon(0.0_dp)
  real(kind=wp), parameter :: hugest = huge(1.0_wp)

! define NaN
  real(real64),  parameter :: nan64 = transfer(-2251799813685248_int64, 1._real64)
  real(real32),  parameter :: nan32 = transfer(-4194304_int32, 1._real32)

! numerical values
  real(kind=wp), parameter :: pi      = 3.141592653589793238462643383279502884197_wp
  real(kind=wp), parameter :: twopi   = 6.283185307179586476925286766559005768394_wp
  real(kind=wp), parameter :: fourpi  = 12.56637061435917295385057353311801153679_wp
  real(kind=wp), parameter :: halfpi  = pi/2.0_wp
  real(kind=wp), parameter :: sqrttwo = sqrt(2.0_wp)

! conversion factor (radian to degree, degree to radian)
  real(kind=wp), parameter :: rad2deg = 180.0_wp/pi
  real(kind=wp), parameter :: deg2rad = pi/180.0_wp

! distances
  real(kind=wp), parameter :: pc2cm  = 3.0856776e18_wp
  real(kind=wp), parameter :: kpc2cm = pc2cm * 1e3_wp
  real(kind=wp), parameter :: au2cm  = 1.4960e13_wp
  real(kind=wp), parameter :: ang2m  = 1.0e-10_wp
  real(kind=wp), parameter :: ang2km = 1.0e-13_wp
  real(kind=wp), parameter :: um2m   = 1.0e-6_wp
  real(kind=wp), parameter :: um2km  = 1.0e-9_wp

! speed of light in km/s
  real(kind=wp), parameter :: speedc  = 2.99792458e5_wp

! h_planck    = Planck's constant in m^2 kg/s
! massH       = Hydrogen mass in kg
! lambda_LyaH = Lya H line (2S1/2-2P1/2) wavelength in um.
! DnuHK_Hz    = Frequency difference between 2P1/2 - 2P3/2 of Hydrogen atom in Hz.
  real(kind=wp), parameter :: h_planck    = 6.62607004e-34_wp
  real(kind=wp), parameter :: massH       = 1.6737236e-27_wp
  real(kind=wp), parameter :: lambda_LyaH = 0.1215673123130_wp
  real(kind=wp), parameter :: g_recoil0   = (h_planck/massH)/(lambda_LyaH*um2m)**2
  !real(kind=wp), parameter :: DnuHK_Hz      = 1.1e10_wp
  real(kind=wp), parameter :: DnuHK_Hz      = 1.08e10_wp
  real(kind=wp), parameter :: DnuHK_Hz_half = DnuHK_Hz/2.0_wp

! maximum number of observers
  integer, parameter :: MAX_OBSERVERS = 99

! dust mass per hydrogen nuclei (g/H)
  real(kind=wp), parameter :: h2dust = 1.87e-26_wp

! kappa_v : dust extinction cross-section (cm^2/g) at V-band
  real(kind=wp), parameter :: kappa_v = 25992.6_wp

! photon type
  type photon_type
     integer(int64):: id
     real(kind=wp) :: nscatt_HI
     real(kind=wp) :: nscatt_dust
     real(kind=wp) :: x,y,z
     real(kind=wp) :: kx,ky,kz
     real(kind=wp) :: mx,my,mz
     real(kind=wp) :: nx,ny,nz
     integer       :: icell,jcell,kcell
     real(kind=wp) :: xfreq
     real(kind=wp) :: xfreq_ref
     real(kind=wp) :: wgt
     real(kind=wp) :: vfy_shear
     logical       :: inside
     ! Stokes parameters
     real(kind=wp) :: I,Q,U,V
  end type photon_type

! cylindrical grid type
! nr, np, nz          - number of cells
! rface, pface, zface - locations of cell faces
! opacity             - opacity per unit length
  type grid_type
     integer :: nx,ny,nz
     integer :: nxfreq
     integer :: i0 = 0, j0 = 0, k0 = 0
     real(kind=wp) :: rmin = -999.0
     real(kind=wp) :: rmax
     real(kind=wp) :: xmin,xmax,xrange
     real(kind=wp) :: ymin,ymax,yrange
     real(kind=wp) :: zmin,zmax,zrange
     real(kind=wp) :: dx,dy,dz
     real(kind=wp) :: xfreq_min, xfreq_max, dxfreq
     real(kind=wp) :: Dfreq_ref, DnuHK_ref, DnuHK_ref_half
     real(kind=wp) :: voigt_amean, Dfreq_mean
     real(kind=wp) :: xcrit, xcrit2
     real(kind=wp), pointer :: xfreq(:)       => null()
     real(kind=wp), pointer :: velocity(:)    => null()
     real(kind=wp), pointer :: dlambda(:)     => null()
     real(kind=wp), pointer :: xface(:)       => null()
     real(kind=wp), pointer :: yface(:)       => null()
     real(kind=wp), pointer :: zface(:)       => null()
     real(kind=wp), pointer :: vfx(:,:,:)     => null()
     real(kind=wp), pointer :: vfy(:,:,:)     => null()
     real(kind=wp), pointer :: vfz(:,:,:)     => null()
     real(kind=wp), pointer :: Dfreq(:,:,:)   => null()
     real(kind=wp), pointer :: voigt_a(:,:,:) => null()
     real(kind=wp), pointer :: rhokap(:,:,:)  => null()
     real(kind=wp), pointer :: Jin(:)         => null()
     real(kind=wp), pointer :: Jout(:)        => null()
     real(kind=wp), pointer :: rhokapD(:,:,:) => null()
     real(kind=wp), pointer :: Jabs(:)        => null()
     real(kind=wp), pointer :: Dfreq_ratio(:,:,:) => null()
     !--- for exoplanet atmosphere model.
     integer,       pointer :: mask(:,:,:)    => null()
     real(kind=wp), pointer :: Jabs2(:)       => null()
     !--- emissivity
     real(kind=wp), pointer :: Pem(:,:,:)     => null()
#if defined CALCJ || defined CALCP || defined CALCPnew
     ! cell indices for radial coordinate (spherical or cylindrical).
     integer, pointer      :: rind_sph(:,:,:) => null()
     integer, pointer      :: rind_cyl(:,:,:) => null()
#endif
#ifdef CALCJ
     real(kind=dp), pointer :: J(:,:,:,:)     => null()
     ! 1D profiles (z or radial profiles)
     real(kind=wp), pointer :: J1(:,:)        => null()
     ! 2D profiles (z and radial profiles for a cylindrically symmetric case)
     real(kind=wp), pointer :: J2(:,:,:)      => null()
#endif
#ifdef CALCP
     real(kind=dp), pointer :: Pa(:,:,:)      => null()
     ! 1D profiles (z or radial profiles)
     real(kind=wp), pointer :: P1(:)          => null()
     ! 2D profiles (z and radial profiles for a cylindrically symmetric case)
     real(kind=wp), pointer :: P2(:,:)        => null()
#endif
#ifdef CALCPnew
     real(kind=dp), pointer :: Pa_new(:,:,:)  => null()
     ! 1D profiles (z or radial profiles)
     real(kind=wp), pointer :: P1_new(:)      => null()
     ! 2D profiles (z and radial profiles for a cylindrically symmetric case)
     real(kind=wp), pointer :: P2_new(:,:)    => null()
#endif
     ! Pa and J should be always defined as double precision (2017-04-25, 05-24).
     ! Pa: single precision results in saturation.
     ! J : scattering distance is too short so that the disatnce is virtually zero with single precision.
  end type grid_type

! input parameters
  type params_type
     integer(kind=int64) :: nphotons = 1e5
     integer       :: nprint       = 1e7
     real(kind=wp) :: no_photons   = 1e5
     real(kind=wp) :: no_print     = 1e7
     integer       :: iseed        = 0
     real(kind=wp) :: luminosity   = 1.0
     real(kind=wp) :: temperature  = 1d4
     real(kind=wp) :: temperature0 = -999.0
     real(kind=wp) :: Dfreq0       = -999.0  ! Doppler frequency for the photon source
     real(kind=wp) :: voigt_a0     = -999.0  ! voigt_a parameter for the photon source
     !--- modified 2020.10.22
     !--- tau0 and N_HI are reagared to be taumax and N_HImax, respectively.
     real(kind=wp) :: taumax       = -999.0
     real(kind=wp) :: tauhomo      = -999.0
     real(kind=wp) :: tau0         = -999.0
     real(kind=wp) :: N_HImax      = -999.0
     real(kind=wp) :: N_HIhomo     = -999.0
     real(kind=wp) :: N_HI         = -999.0
     real(kind=wp) :: atau3
     real(kind=wp) :: f12          = 0.4126_wp
     real(kind=wp) :: A21          = 6.265e8_wp
     integer       :: atom_no      = 1
     real(kind=wp) :: cross0       = 0.02656_wp/sqrt(pi)*0.4126
     !--- wavelength for Lyman-alpha H (1S1/2 - 2P1/2) line in um.
     real(kind=wp) :: lambda0      = lambda_LyaH
     !--- parameter for Hubble-like expansion
     real(kind=wp) :: Vexp            = 0.0_wp
     !--- parameters for Song, Seon, & Hwang (2020) model
     real(kind=wp) :: Vpeak         = 0.0_wp
     real(kind=wp) :: rpeak         = 0.0_wp
     real(kind=wp) :: DeltaV        = 0.0_wp
     real(kind=wp) :: source_rscale = 0.0_wp
     real(kind=wp) :: sersic_m      = 1.0_wp
     real(kind=wp) :: Reff          = 0.0_wp
     !---
     logical       :: comoving_source  = .true.
     logical       :: recoil           = .false.
     logical       :: core_skip        = .false.
     logical       :: core_skip_global = .false.
     !--- grid geometry parameters (default geometry is a rectagle, 2020.08.28)
     logical       :: xyz_symmetry    = .false.
     logical       :: xy_periodic     = .false.
     logical       :: z_symmetry      = .false.
     logical       :: save_sightline_tau = .false.
     !character(len=128) :: geometry = 'sphere'
     character(len=128) :: geometry = ''
     character(len=128) :: velocity_type = ''
     !--- grid parameters
     integer       :: nx   = 1
     integer       :: ny   = 1
     integer       :: nz   = 11
     integer       :: nr   = -999
     real(kind=wp) :: xmax = 1.0_wp
     real(kind=wp) :: ymax = 1.0_wp
     real(kind=wp) :: zmax = 1.0_wp
     real(kind=wp) :: xmin = nan64
     real(kind=wp) :: ymin = nan64
     real(kind=wp) :: zmin = nan64
     real(kind=wp) :: rmin = -999.0
     real(kind=wp) :: rmax = -999.0
     real(kind=wp) :: source_rmax = -999.0_wp
     !--- density gamma
     real(kind=wp) :: density_rscale = -999.9_wp
     !--- location of a point source
     real(kind=wp) :: xs_point = 0.0_wp
     real(kind=wp) :: ys_point = 0.0_wp
     real(kind=wp) :: zs_point = 0.0_wp
     !--- frequency parameters
     real(kind=wp) :: xfreq0    = 0.0_wp
     real(kind=wp) :: xfreq_min = 0.0_wp
     real(kind=wp) :: xfreq_max = 0.0_wp
     integer       :: nxfreq    = 121
     !--- TIGRESS model
     !--- background y-velocity due to shearing box
     !--- vy0 = -q * Omega * x (Omega in km/s/kpc and x in kpc)
     !--- q   = 1 (flat rotation).
     real(kind=wp) :: q         = 1.0_wp
     real(kind=wp) :: Omega     = 0.0_wp
     !--- continuum parameters
     real(kind=wp) :: continuum_slope     = 0.0_wp
     real(kind=wp) :: continuum_wgt0      = 1.0_wp
     real(kind=wp) :: continuum_abs_depth = 0.0_wp
     real(kind=wp) :: continuum_abs_width = 1.0_wp
     real(kind=wp) :: source_zscale       = 0.0_wp
     real(kind=wp) :: distance2cm          = -999.9_wp
     real(kind=wp) :: gaussian_width_vel   = 12.843374_wp
     character(len=128) :: distance_unit   = 'kpc'
     character(len=128) :: source_geometry = 'point'
     !character(len=128) :: spectral_type   = 'monochromatic'
     character(len=128) :: spectral_type   = 'voigt'
     !--- input/output filename related
     character(len=128) :: base_name       = ''
     character(len=128) :: out_file        = ''
     logical            :: out_merge       = .false.
     integer            :: out_bitpix      = -32
     !--- master-slave algorithm
     integer       :: num_send_at_once     = 1
     logical       :: use_master_slave     = .true.
     logical       :: use_dynamic_schedule = .false.
     logical       :: use_reduced_wgt      = .false.
     logical       :: use_stokes           = .false.
     !--- output-related parameters
     logical       :: save_all         = .false.
     logical       :: save_Jin         = .true.
     logical       :: save_Jabs        = .true.
     logical       :: save_backup      = .false.
     logical       :: save_direc0      = .false.
     logical       :: save_all_photons = .false.
     logical       :: save_input_grid  = .false.
     !--- Collisional Ionization Equilibrium
     logical       :: use_cie_condition = .false.
     !--- Dust-related parameters
     real(kind=wp) :: hgg           = 0.6761
     real(kind=wp) :: albedo        = 0.3253
     real(kind=wp) :: cext_dust     = 1.6059e-21
     real(kind=wp) :: metalZ        = 0.0
     real(kind=wp) :: DGR           = 0.0
     character(len=128) :: scatt_mat_file = ''
     !--- external density file
     character(len=128) :: input_field = ''
     character(len=128) :: dens_file  = ''
     character(len=128) :: temp_file  = ''
     character(len=128) :: velo_file  = ''
     character(len=128) :: emiss_file = ''
     integer            :: reduce_factor = 1
     integer            :: centering     = 0
     !--- number of scatterings, excutiona time, these are not input parameters.
     real(kind=wp) :: nscatt_HI   = 0.0_wp
     real(kind=wp) :: nscatt_dust = 0.0_wp
     real(kind=wp) :: nscatt_tot  = 0.0_wp
     real(kind=wp) :: exetime     = 0.0_wp
#ifdef PEELINGOFF
     !--- parameters for PEELING-OFF
     integer :: nobs = 1
     integer :: nxim = 129
     integer :: nyim = 129
     logical :: peel_average = .false.
     real(kind=wp) :: distance                         = nan64
     real(kind=wp) :: inclination_angle(MAX_OBSERVERS) = nan64
     real(kind=wp) :: position_angle(MAX_OBSERVERS)    = nan64
     real(kind=wp) :: phase_angle(MAX_OBSERVERS)       = nan64
     real(kind=wp) :: alpha(MAX_OBSERVERS) = nan64
     real(kind=wp) :: beta(MAX_OBSERVERS)  = nan64
     real(kind=wp) :: gamma(MAX_OBSERVERS) = nan64
     real(kind=wp) :: obsx(MAX_OBSERVERS)  = nan64
     real(kind=wp) :: obsy(MAX_OBSERVERS)  = nan64
     real(kind=wp) :: obsz(MAX_OBSERVERS)  = nan64
     real(kind=wp) :: dxim = nan64
     real(kind=wp) :: dyim = nan64
#endif
  end type params_type

#ifdef PEELINGOFF
  !--- observer type
  type observer_type
     real(kind=wp) :: x,y,z
     real(kind=wp) :: inclination_angle = 0.0_wp
     real(kind=wp) :: position_angle    = 0.0_wp
     real(kind=wp) :: phase_angle       = 0.0_wp
     real(kind=wp) :: distance          = -999.9_wp
     real(kind=wp) :: rmatrix(3,3)
     integer       :: nxim, nyim
     integer       :: nxfreq
     real(kind=wp) :: dxim,dyim
     real(kind=wp) :: steradian_pix
     real(kind=wp), pointer :: tau_HI(:,:,:) => null()
     real(kind=wp), pointer :: N_HI(:,:)     => null()
     real(kind=wp), pointer :: tau_dust(:,:) => null()
     real(kind=wp), pointer :: scatt(:,:,:)  => null()
     real(kind=wp), pointer :: direc(:,:,:)  => null()
     real(kind=wp), pointer :: direc0(:,:,:) => null()
     real(kind=wp), pointer :: I(:,:,:)      => null()
     real(kind=wp), pointer :: Q(:,:,:)      => null()
     real(kind=wp), pointer :: U(:,:,:)      => null()
     real(kind=wp), pointer :: V(:,:,:)      => null()
     real(kind=wp), pointer :: radial_spec(:,:) => null()
     real(kind=wp), pointer :: radial_pol(:)    => null()
     real(kind=wp), pointer :: radial_r(:)      => null()
     real(kind=wp), pointer :: radial_I(:)      => null()
     real(kind=wp), pointer :: radial_Q(:)      => null()
     real(kind=wp), pointer :: radial_U(:)      => null()
     real(kind=wp), pointer :: radial_V(:)      => null()
  end type observer_type
#endif

  type all_photons_type
     real(kind=wp), pointer :: rp0(:)         => null()  !- radius of input   photon on a projected plane
     real(kind=wp), pointer :: rp(:)          => null()  !- radius of output  photon on a projected plane
     real(kind=wp), pointer :: xfreq1(:)      => null()  !- initial frequency
     real(kind=wp), pointer :: xfreq2(:)      => null()  !- final   frequency
     real(kind=wp), pointer :: nscatt_HI(:)   => null()  !- number of scatterings due to HI
     real(kind=wp), pointer :: nscatt_dust(:) => null()  !- number of scatterings due to dust
     real(kind=wp), pointer :: I(:)           => null()  !- Stokes I, weight just in case, when dust exists.
     real(kind=wp), pointer :: Q(:)           => null()  !- Stokes Q
     real(kind=wp), pointer :: U(:)           => null()  !- Stokes U
     real(kind=wp), pointer :: V(:)           => null()  !- Stokes V
  end type all_photons_type

! scattering related global variables
  type scattering_matrix_type
     integer :: nPDF
     real(kind=wp), pointer :: coss(:)      => null()
     real(kind=wp), pointer :: S11(:)       => null()
     real(kind=wp), pointer :: S12(:)       => null()
     real(kind=wp), pointer :: S33(:)       => null()
     real(kind=wp), pointer :: S34(:)       => null()
     real(kind=wp), pointer :: cosp(:)      => null()
     real(kind=wp), pointer :: phase_CDF(:) => null()
  end type scattering_matrix_type

  type params_mpi
     !--- p_rank    = process (thread) rank
     !--- h_rank    = host rank (rank defined for each process of a node)
     !--- nproc     = number of total processes
     !--- hostcomm  = communicator for host
     !--- SAME_HRANK_COMM  = communicator for processes with the same h_rank
     !--- SAME_HRANK_NPROC = number of processes with the same h_rank
     integer :: p_rank
     integer :: h_rank
     integer :: nproc
     integer :: hostcomm
     integer :: SAME_HRANK_COMM
     integer :: SAME_HRANK_NPROC
  end type params_mpi

  ! continuum parameters
  type params_cont
     real(kind=wp) :: slope
     real(kind=wp) :: slope_max
     real(kind=wp) :: x0, wgt0
     real(kind=wp) :: abs_depth, abs_width
  end type params_cont

  ! globar parameters
  type(params_type) :: par
  type(params_mpi)  :: mpar
  type(params_cont) :: cpar
#ifdef PEELINGOFF
  !type(observer_type) :: observer
  type(observer_type), allocatable :: observer(:)
#endif
  type(scattering_matrix_type) :: scatt_mat
  type(all_photons_type) :: allph

  !------------------------------------------------------------
  ! define procedure pointers
  procedure(raytrace_tau), pointer :: raytrace_to_tau => null()
  abstract interface
     subroutine raytrace_tau(photon,grid,tau_in)
     import
     type(photon_type), intent(inout) :: photon
     type(grid_type),   intent(inout) :: grid
     real(kind=wp),     intent(in)    :: tau_in
     end subroutine raytrace_tau
  end interface

  procedure(raytrace_edge), pointer :: raytrace_to_edge => null()
  abstract interface
     subroutine raytrace_edge(photon,grid,tau_out)
     import
     type(photon_type), intent(in)  :: photon
     type(grid_type),   intent(in)  :: grid
     real(kind=wp),     intent(out) :: tau_out
     end subroutine raytrace_edge
  end interface

  procedure(raytrace_edge_tau_HI), pointer :: raytrace_to_edge_tau_HI => null()
  abstract interface
     subroutine raytrace_edge_tau_HI(photon,grid,tau_HI)
     import
     type(photon_type), intent(in)  :: photon
     type(grid_type),   intent(in)  :: grid
     real(kind=wp),     intent(out) :: tau_HI
     end subroutine raytrace_edge_tau_HI
  end interface

  procedure(raytrace_edge_column), pointer :: raytrace_to_edge_column => null()
  abstract interface
     subroutine raytrace_edge_column(photon,grid,N_HI,tau_dust)
     import
     type(photon_type), intent(in)  :: photon
     type(grid_type),   intent(in)  :: grid
     real(kind=wp),     intent(out) :: N_HI
     real(kind=wp),     intent(out) :: tau_dust
     end subroutine raytrace_edge_column
  end interface

  procedure(scatter_routine), pointer :: scatter_dust      => null()
  procedure(scatter_routine), pointer :: scatter_resonance => null()
  abstract interface
     subroutine scatter_routine(photon,grid)
     import
     type(photon_type), intent(inout) :: photon
     type(grid_type),   intent(inout) :: grid
   end subroutine scatter_routine
  end interface

  procedure(run_sim), pointer :: run_simulation => null()
  abstract interface
     subroutine run_sim(grid)
     import
     type(grid_type),   intent(inout) :: grid
     end subroutine run_sim
  end interface

  procedure(write_out), pointer :: write_output => null()
  abstract interface
     subroutine write_out(filename,grid)
     import
     character(len=*),    intent(in)    :: filename
     type(grid_type),     intent(inout) :: grid
     end subroutine write_out
  end interface
  !------------------------------------------------------------
end module define
