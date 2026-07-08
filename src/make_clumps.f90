!=========================================================================
! make_clumps.f90
!
! Standalone Fortran program to generate a clump population (positions,
! bulk velocities, clump physical parameters) and write a FITS or
! HDF5 file with the same schema as LaRT's write_clumps_info().
!
! This is the standalone Fortran equivalent of python/make_clumps.py and
! is a CLI-driven replacement for the old MPI/namelist-based make_clumps
! driver.  Output is consumed by LaRT via par%clump_input_file = '<file>'.
!
! Build:
!   make make_clumps        (see Makefile target)
!
! Usage:
!   ./make_clumps.x [options]
!
! Options match python/make_clumps.py.  See --help for the full list.
!=========================================================================
program make_clumps
  use define,     only: wp, pi, twopi, fourpi, deg2rad
  use iofile_mod
  use voigt_mod,  only: voigt
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  implicit none

  !-----------------------------------------------------------------------
  ! Physical constants (mirror line_mod.f90 / define.f90)
  !-----------------------------------------------------------------------
  real(wp), parameter :: SIGMA_0     = 0.026540083434_wp
  real(wp), parameter :: VTHERM1_AMU = 0.12895319011972164_wp
  real(wp), parameter :: UM2KM       = 1.0e-9_wp
  real(wp), parameter :: CONST_TOL   = 1.0e-3_wp  ! constant-column threshold

  !-----------------------------------------------------------------------
  ! Atomic data table for each line (parallel arrays)
  !-----------------------------------------------------------------------
  integer, parameter :: N_LINES = 10
  character(len=12), parameter :: LINE_NAME(N_LINES) = [character(len=12) :: &
       'ly_alpha   ', 'CIV_1548   ', 'NV_1239    ', 'OVI_1032   ', &
       'NaI_D      ', 'CaII_HK    ', 'MgII_2796  ', 'SiIV_1394  ', &
       'AlII_1671  ', 'CII_1334   ']
  real(wp), parameter :: LINE_WL0(N_LINES) = [ &
       0.1215668237310_wp, 0.1548187_wp, 0.1238821_wp, 0.1031912_wp, &
       0.5891583253_wp,    0.3934777_wp, 0.2796352_wp, 0.1393755_wp, &
       0.1670787_wp,       0.1334532_wp ]
  real(wp), parameter :: LINE_DAMP(N_LINES) = [ &
       6.2649e8_wp, 2.647e8_wp, 3.390e8_wp, 4.137e8_wp, &
       6.153e7_wp, 1.446667e8_wp, 2.590e8_wp, 8.743e8_wp, &
       1.39e9_wp,  2.880e8_wp ]
  real(wp), parameter :: LINE_F12(N_LINES) = [ &
       0.4164_wp, 0.190_wp, 0.156_wp, 0.1325_wp, &
       0.641_wp,  0.682_wp, 0.6080_wp, 0.513_wp, &
       1.77_wp,   0.129_wp ]
  real(wp), parameter :: LINE_MASS(N_LINES) = [ &
       1.00797_wp, 12.011_wp, 14.007_wp, 15.999_wp, &
       22.990_wp,  40.078_wp, 24.305_wp, 28.085_wp, &
       26.982_wp,  12.011_wp ]

  !-----------------------------------------------------------------------
  ! Distance unit -> cm (subset of define.f90 unit table)
  !-----------------------------------------------------------------------
  integer, parameter :: N_UNITS = 8
  character(len=8), parameter :: UNIT_NAME(N_UNITS) = [character(len=8) :: &
       'cm   ', 'm    ', 'km   ', 'au   ', 'pc   ', 'kpc  ', 'Mpc  ', 'Rsun ']
  real(wp), parameter :: UNIT_CM(N_UNITS) = [ &
       1.0_wp, 1.0e2_wp, 1.0e5_wp, 1.495978707e13_wp, &
       3.0856775814913673e18_wp, 3.0856775814913673e21_wp, &
       3.0856775814913673e24_wp, 6.957e10_wp ]

  !-----------------------------------------------------------------------
  ! CLI parameters
  !-----------------------------------------------------------------------
  real(wp) :: rmax_in      = 0.0_wp
  real(wp) :: rmin_in      = 0.0_wp
  real(wp) :: cone_opening = 0.0_wp
  character(len=16) :: distance_unit = ''

  real(wp) :: clump_radius   = 0.0_wp
  real(wp) :: clump_N_clumps = 0.0_wp
  real(wp) :: clump_f_vol    = 0.0_wp
  real(wp) :: clump_f_cov    = 0.0_wp

  real(wp) :: clump_tau0 = 0.0_wp
  real(wp) :: clump_NHI  = 0.0_wp
  real(wp) :: clump_nH   = 0.0_wp
  real(wp) :: taumax     = 0.0_wp
  real(wp) :: N_HImax    = 0.0_wp

  real(wp) :: temperature       = 1.0e4_wp
  real(wp) :: clump_temperature = -1.0_wp
  character(len=16) :: line_id = 'ly_alpha'

  logical :: clump_fully_inside  = .true.
  logical :: clump_allow_overlap = .false.

  character(len=32) :: velocity_type = ''
  real(wp) :: Vexp = 0.0_wp
  real(wp) :: Vrot = 0.0_wp
  real(wp) :: velocity_alpha = 1.0_wp
  real(wp) :: rinner = 0.0_wp
  real(wp) :: rpeak  = 0.0_wp
  real(wp) :: Vpeak  = 0.0_wp
  real(wp) :: DeltaV = 0.0_wp
  real(wp) :: Vxin = 0.0_wp, Vyin = 0.0_wp, Vzin = 0.0_wp
  real(wp) :: clump_sigma_v = 0.0_wp

  character(len=32) :: clump_radius_profile  = 'constant'
  character(len=32) :: clump_density_profile = 'constant'
  character(len=32) :: clump_number_profile  = 'constant'
  real(wp) :: clump_radius_alpha  = 0.0_wp
  real(wp) :: clump_density_alpha = 0.0_wp
  real(wp) :: clump_number_alpha  = 0.0_wp
  real(wp) :: clump_radius_r0  = 0.0_wp
  real(wp) :: clump_density_r0 = 0.0_wp
  real(wp) :: clump_number_r0  = 0.0_wp

  character(len=512) :: output_file = ''
  integer(int64)     :: seed_in     = 12345_int64

  !-----------------------------------------------------------------------
  ! Derived runtime state
  !-----------------------------------------------------------------------
  integer  :: line_idx
  real(wp) :: wavelength0, damping, f12, mass_amu
  real(wp) :: cross0, vtherm1, vtherm_ref, Dfreq_ref, voigt_a_ref, voigt0_ref
  real(wp) :: temp_cl, distance2cm
  real(wp) :: sphere_R, r_min_clump, base_radius
  real(wp) :: cos_cone
  logical  :: profiles_active

  ! Radial profile tables (4001 points, allocated only if profiles_active)
  integer, parameter :: NPROF = 4001
  real(wp), allocatable :: prof_r(:), prof_shape_n(:), prof_shape_r(:), &
                           prof_shape_d(:), prof_cdf(:)
  real(wp) :: r_cl_max

  ! Population
  integer(int64) :: N_target, N_placed
  real(wp), allocatable :: cl_x(:), cl_y(:), cl_z(:), cl_r(:)
  real(wp), allocatable :: cl_vx(:), cl_vy(:), cl_vz(:)
  real(wp), allocatable :: cl_rhokap(:), cl_temp(:)
  real(wp) :: A_norm, rhokap_ref, fvol_input_est, fcov_input_est

  ! RNG state
  integer(int64) :: rng_state

  ! Driver
  integer :: status

  !-----------------------------------------------------------------------
  ! Parse CLI
  !-----------------------------------------------------------------------
  call parse_args()
  call validate_args()
  call rng_init(seed_in)

  sphere_R    = rmax_in
  r_min_clump = max(0.0_wp, rmin_in)
  base_radius = clump_radius

  cos_cone = -1.0_wp
  if (cone_opening > 0.0_wp .and. cone_opening < 90.0_wp) &
       cos_cone = cos(cone_opening * deg2rad)

  call resolve_line_id(line_id, line_idx)
  wavelength0 = LINE_WL0(line_idx)
  damping     = LINE_DAMP(line_idx)
  f12         = LINE_F12(line_idx)
  mass_amu    = LINE_MASS(line_idx)

  temp_cl = temperature
  if (clump_temperature >= 0.0_wp) temp_cl = clump_temperature

  distance2cm = unit_to_cm(distance_unit)

  call line_reference(temp_cl, cross0, vtherm1, vtherm_ref, Dfreq_ref, voigt_a_ref)
  voigt0_ref = voigt(0.0_wp, voigt_a_ref)

  !-----------------------------------------------------------------------
  ! Radial profile tables (only if any axis is non-constant)
  !-----------------------------------------------------------------------
  profiles_active = .not. ( is_constant(clump_radius_profile)  .and. &
                            is_constant(clump_density_profile) .and. &
                            is_constant(clump_number_profile) )
  if (profiles_active) call build_radial_profile()
  if (.not. profiles_active) r_cl_max = base_radius

  !-----------------------------------------------------------------------
  ! N_clumps + reference rhokap
  !-----------------------------------------------------------------------
  call derive_N_clumps_and_norm()
  call derive_rhokap_ref()

  write(*,'(a,i0)')     ' Clumps: N_clumps  = ', N_target
  write(*,'(a,f10.6)')  ' Clumps: f_vol     = ', fvol_input_est
  write(*,'(a,f10.5)')  ' Clumps: f_cov     = ', fcov_input_est
  write(*,'(a,2f10.5)') ' Clumps: rmin/rmax = ', r_min_clump, sphere_R
  write(*,'(a,es12.4)') ' Clumps: cl_rhokap = ', rhokap_ref
  write(*,'(a,f10.5)')  ' Clumps: voigt_a   = ', voigt_a_ref
  write(*,'(a,es12.4)') ' Clumps: cl_Dfreq  = ', Dfreq_ref
  if (profiles_active) then
     write(*,'(a,a)') ' Clumps: radius profile  = ', trim(clump_radius_profile)
     write(*,'(a,a)') ' Clumps: density profile = ', trim(clump_density_profile)
     write(*,'(a,a)') ' Clumps: number profile  = ', trim(clump_number_profile)
  end if

  !-----------------------------------------------------------------------
  ! RSA placement
  !-----------------------------------------------------------------------
  allocate(cl_x(N_target), cl_y(N_target), cl_z(N_target), cl_r(N_target))
  call rsa_place()

  !-----------------------------------------------------------------------
  ! Clump physics (R_CLUMP already set in cl_r; assign RHOKAP, TEMP)
  !-----------------------------------------------------------------------
  allocate(cl_rhokap(N_target), cl_temp(N_target))
  call assign_perclump_physics()

  !-----------------------------------------------------------------------
  ! Velocities
  !-----------------------------------------------------------------------
  allocate(cl_vx(N_target), cl_vy(N_target), cl_vz(N_target))
  call assign_velocities()

  !-----------------------------------------------------------------------
  ! Output
  !-----------------------------------------------------------------------
  call write_output()

  call print_summary()

contains

  !=======================================================================
  ! line atomic data
  !=======================================================================
  subroutine resolve_line_id(name, idx)
    character(len=*), intent(in)  :: name
    integer,          intent(out) :: idx
    integer :: i
    idx = -1
    do i = 1, N_LINES
       if (trim(LINE_NAME(i)) == trim(name)) then
          idx = i
          return
       end if
    end do
    write(*,'(a,a)') 'ERROR: unknown line_id ', trim(name)
    write(*,'(a)')   '       Supported lines:'
    do i = 1, N_LINES
       write(*,'(a,a)') '         ', trim(LINE_NAME(i))
    end do
    stop 1
  end subroutine

  subroutine line_reference(T, cross0_o, vtherm1_o, vtherm_o, Dfreq_o, voigt_a_o)
    real(wp), intent(in)  :: T
    real(wp), intent(out) :: cross0_o, vtherm1_o, vtherm_o, Dfreq_o, voigt_a_o
    cross0_o   = SIGMA_0 / sqrt(pi) * f12
    vtherm1_o  = VTHERM1_AMU / sqrt(mass_amu)
    vtherm_o   = vtherm1_o * sqrt(T)
    Dfreq_o    = vtherm_o / (wavelength0 * UM2KM)
    voigt_a_o  = (damping / fourpi) / Dfreq_o
  end subroutine

  function unit_to_cm(s) result(d2cm)
    character(len=*), intent(in) :: s
    real(wp) :: d2cm
    integer :: i
    d2cm = 0.0_wp
    if (len_trim(s) == 0) return
    do i = 1, N_UNITS
       if (trim(UNIT_NAME(i)) == trim(s)) then
          d2cm = UNIT_CM(i)
          return
       end if
    end do
    write(*,'(a,a)') 'ERROR: unknown distance_unit ', trim(s)
    stop 1
  end function

  logical function is_constant(s)
    character(len=*), intent(in) :: s
    is_constant = (len_trim(s) == 0 .or. trim(adjustl(lowercase(s))) == 'constant')
  end function

  function lowercase(s) result(t)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: t
    integer :: i, c
    t = s
    do i = 1, len_trim(t)
       c = iachar(t(i:i))
       if (c >= iachar('A') .and. c <= iachar('Z')) t(i:i) = achar(c + 32)
    end do
  end function

  !=======================================================================
  ! Profile shape (mirror clump_mod.f90:profile_shape, python:profile_shape)
  !=======================================================================
  function profile_shape(name, alpha, r0, r_arr) result(s)
    character(len=*), intent(in) :: name
    real(wp), intent(in) :: alpha, r0
    real(wp), intent(in) :: r_arr(:)
    real(wp) :: s(size(r_arr))
    character(len=len(name)) :: nm
    real(wp) :: r_floor, r0_eff
    integer  :: i

    nm = trim(adjustl(lowercase(name)))
    select case (trim(nm))
    case ('', 'constant')
       s = 1.0_wp
    case ('powerlaw', 'power_law')
       if (r0 <= 0.0_wp) then
          s = 1.0_wp
       else
          r_floor = 0.05_wp * r0
          r0_eff  = max(r0, r_floor)
          do i = 1, size(r_arr)
             s(i) = (max(r_arr(i), r_floor) / r0_eff) ** (-alpha)
          end do
       end if
    case ('gaussian')
       if (r0 <= 0.0_wp) then
          s = 1.0_wp
       else
          s = exp(-(r_arr / r0)**2)
       end if
    case ('exponential')
       if (r0 <= 0.0_wp) then
          s = 1.0_wp
       else
          s = exp(-r_arr / r0)
       end if
    case default
       write(*,'(a,a)') 'ERROR: unknown profile shape ', trim(name)
       stop 1
    end select
  end function

  subroutine build_radial_profile()
    integer  :: i
    real(wp) :: r0_n, r0_r, r0_d
    real(wp) :: integrand(NPROF), cdf(NPROF)
    real(wp) :: total
    real(wp), allocatable :: rcl_local(:)

    allocate(prof_r(NPROF), prof_shape_n(NPROF), prof_shape_r(NPROF), &
             prof_shape_d(NPROF), prof_cdf(NPROF))

    do i = 1, NPROF
       prof_r(i) = sphere_R * real(i-1, wp) / real(NPROF-1, wp)
    end do

    r0_n = clump_number_r0;  if (r0_n <= 0.0_wp) r0_n = sphere_R
    r0_r = clump_radius_r0;  if (r0_r <= 0.0_wp) r0_r = sphere_R
    r0_d = clump_density_r0; if (r0_d <= 0.0_wp) r0_d = sphere_R

    prof_shape_n = profile_shape(clump_number_profile,  clump_number_alpha,  r0_n, prof_r)
    prof_shape_r = profile_shape(clump_radius_profile,  clump_radius_alpha,  r0_r, prof_r)
    prof_shape_d = profile_shape(clump_density_profile, clump_density_alpha, r0_d, prof_r)

    ! Zero shape_number inside the inner cavity (matches clump_mod.f90)
    do i = 1, NPROF
       if (prof_r(i) < r_min_clump) prof_shape_n(i) = 0.0_wp
    end do

    integrand = prof_shape_n * prof_r * prof_r
    cdf(1) = 0.0_wp
    do i = 2, NPROF
       cdf(i) = cdf(i-1) + 0.5_wp * (integrand(i) + integrand(i-1)) &
                         * (prof_r(i) - prof_r(i-1))
    end do
    total = cdf(NPROF)
    if (total > 0.0_wp) then
       cdf = cdf / total
    else
       ! degenerate -> uniform-in-volume fallback
       do i = 1, NPROF
          cdf(i) = (prof_r(i) / sphere_R)**3
       end do
    end if
    prof_cdf = cdf

    ! r_cl_max: largest clump radius the profile may produce
    allocate(rcl_local(NPROF))
    rcl_local = base_radius * prof_shape_r
    do i = 1, NPROF
       if (prof_shape_n(i) <= 0.0_wp) rcl_local(i) = 0.0_wp
    end do
    r_cl_max = max(base_radius, maxval(rcl_local))
    deallocate(rcl_local)
  end subroutine

  ! Linear interpolation: y = interp(x, prof_r, prof_<table>)
  function interp_prof(x, table) result(y)
    real(wp), intent(in) :: x
    real(wp), intent(in) :: table(:)
    real(wp) :: y
    real(wp) :: t
    integer  :: i
    if (x <= prof_r(1)) then
       y = table(1); return
    end if
    if (x >= prof_r(NPROF)) then
       y = table(NPROF); return
    end if
    i = int(x / sphere_R * real(NPROF-1, wp)) + 1
    if (i < 1) i = 1
    if (i >= NPROF) i = NPROF - 1
    do while (prof_r(i+1) < x .and. i < NPROF-1); i = i + 1; end do
    do while (prof_r(i) > x   .and. i > 1);       i = i - 1; end do
    t = (x - prof_r(i)) / (prof_r(i+1) - prof_r(i))
    y = table(i) + t * (table(i+1) - table(i))
  end function

  ! Inverse-CDF sample: given u in [0,1], return r
  function sample_r_from_cdf(u) result(r)
    real(wp), intent(in) :: u
    real(wp) :: r, t
    integer  :: lo, hi, mid
    lo = 1; hi = NPROF
    do while (hi - lo > 1)
       mid = (lo + hi) / 2
       if (prof_cdf(mid) <= u) then
          lo = mid
       else
          hi = mid
       end if
    end do
    if (prof_cdf(hi) > prof_cdf(lo)) then
       t = (u - prof_cdf(lo)) / (prof_cdf(hi) - prof_cdf(lo))
    else
       t = 0.0_wp
    end if
    r = prof_r(lo) + t * (prof_r(hi) - prof_r(lo))
  end function

  !=======================================================================
  ! N_clumps + A normalization
  !=======================================================================
  subroutine derive_N_clumps_and_norm()
    real(wp) :: V_shell, tc1, fcov_unit, vint1, fvol_unit
    real(wp) :: rcl_local(NPROF), integrand(NPROF)
    integer  :: i

    A_norm = 0.0_wp

    if (.not. profiles_active) then
       if (clump_N_clumps > 0.0_wp) then
          N_target = nint(clump_N_clumps, int64)
       else if (clump_f_vol > 0.0_wp) then
          N_target = nint(clump_f_vol * (sphere_R**3 - r_min_clump**3) &
                          / base_radius**3, int64)
       else if (clump_f_cov > 0.0_wp) then
          N_target = nint((4.0_wp/3.0_wp) * clump_f_cov &
                  * (sphere_R**2 + sphere_R*r_min_clump + r_min_clump**2) &
                  / base_radius**2, int64)
       else
          write(*,'(a)') 'ERROR: specify --clump_N_clumps, --clump_f_vol, or --clump_f_cov'
          stop 1
       end if
       if (N_target <= 0_int64) N_target = 1_int64
       fvol_input_est = real(N_target,wp) * base_radius**3 &
                       / max(sphere_R**3 - r_min_clump**3, tiny(1.0_wp))
       fcov_input_est = 0.75_wp * real(N_target,wp) * base_radius**2 &
                       / max(sphere_R**2 + sphere_R*r_min_clump + r_min_clump**2, &
                             tiny(1.0_wp))
       return
    end if

    ! Profile case
    integrand = prof_shape_n * fourpi * prof_r * prof_r
    tc1 = trapz(integrand)

    if (clump_N_clumps > 0.0_wp) then
       N_target = nint(clump_N_clumps, int64)
       A_norm   = real(N_target, wp) / max(tc1, tiny(1.0_wp))
    else if (clump_f_vol > 0.0_wp) then
       ! v_int1 = ∫ A=1 * shape_n * 4πr^2 * (4π/3)(base*shape_r)^3 dr
       rcl_local = base_radius * prof_shape_r
       integrand = prof_shape_n * fourpi * prof_r * prof_r &
                 * (fourpi/3.0_wp) * rcl_local**3
       do i = 1, NPROF
          if (prof_r(i) < r_min_clump) integrand(i) = 0.0_wp
       end do
       vint1     = trapz(integrand)
       V_shell   = (fourpi/3.0_wp) * (sphere_R**3 - r_min_clump**3)
       fvol_unit = vint1 / max(V_shell, tiny(1.0_wp))
       A_norm    = clump_f_vol / max(fvol_unit, tiny(1.0_wp))
       N_target  = nint(A_norm * tc1, int64)
    else if (clump_f_cov > 0.0_wp) then
       rcl_local = base_radius * prof_shape_r
       integrand = prof_shape_n * pi * rcl_local**2
       do i = 1, NPROF
          if (prof_r(i) < r_min_clump) integrand(i) = 0.0_wp
       end do
       fcov_unit = trapz(integrand)
       A_norm    = clump_f_cov / max(fcov_unit, tiny(1.0_wp))
       N_target  = nint(A_norm * tc1, int64)
    else
       write(*,'(a)') 'ERROR: specify --clump_N_clumps, --clump_f_vol, or --clump_f_cov'
       stop 1
    end if
    if (N_target <= 0_int64) N_target = 1_int64

    ! Realized f_vol, f_cov estimates from the profile
    rcl_local = base_radius * prof_shape_r
    integrand = A_norm * prof_shape_n * fourpi * prof_r * prof_r &
              * (fourpi/3.0_wp) * rcl_local**3
    do i = 1, NPROF
       if (prof_r(i) < r_min_clump) integrand(i) = 0.0_wp
    end do
    V_shell        = (fourpi/3.0_wp) * (sphere_R**3 - r_min_clump**3)
    fvol_input_est = trapz(integrand) / max(V_shell, tiny(1.0_wp))

    integrand = A_norm * prof_shape_n * pi * rcl_local**2
    do i = 1, NPROF
       if (prof_r(i) < r_min_clump) integrand(i) = 0.0_wp
    end do
    fcov_input_est = trapz(integrand)
  end subroutine

  function trapz(y) result(s)
    real(wp), intent(in) :: y(:)
    real(wp) :: s
    integer  :: i
    s = 0.0_wp
    do i = 2, size(y)
       s = s + 0.5_wp * (y(i) + y(i-1)) * (prof_r(i) - prof_r(i-1))
    end do
  end function

  !=======================================================================
  ! Reference opacity from one of tau0 / NHI / nH / taumax / N_HImax
  !=======================================================================
  subroutine derive_rhokap_ref()
    real(wp) :: GF
    real(wp) :: rcl_local(NPROF), integrand(NPROF)
    integer  :: i

    if (clump_tau0 > 0.0_wp) then
       rhokap_ref = clump_tau0 / (voigt0_ref * base_radius)
       return
    end if
    if (clump_NHI > 0.0_wp) then
       rhokap_ref = clump_NHI * cross0 / (Dfreq_ref * base_radius)
       return
    end if
    if (clump_nH > 0.0_wp) then
       if (distance2cm <= 0.0_wp) then
          write(*,'(a)') 'ERROR: --clump_nH requires --distance_unit'
          stop 1
       end if
       rhokap_ref = clump_nH * cross0 * distance2cm / Dfreq_ref
       return
    end if
    if (taumax > 0.0_wp .or. N_HImax > 0.0_wp) then
       if (.not. profiles_active) then
          GF = real(N_target,wp) * base_radius**3 / &
               max(sphere_R**2 + sphere_R*r_min_clump + r_min_clump**2, tiny(1.0_wp))
       else
          rcl_local = base_radius * prof_shape_r
          integrand = A_norm * prof_shape_n * (4.0_wp/3.0_wp) * rcl_local**3
          do i = 1, NPROF
             if (prof_r(i) < r_min_clump) integrand(i) = 0.0_wp
          end do
          GF = trapz(integrand)
       end if
       if (GF <= 0.0_wp) then
          write(*,'(a)') 'ERROR: cannot back-solve opacity (geometric factor is zero)'
          stop 1
       end if
       if (taumax > 0.0_wp) then
          rhokap_ref = taumax / (GF * voigt0_ref)
       else
          rhokap_ref = N_HImax * cross0 / (GF * Dfreq_ref)
       end if
       return
    end if
    write(*,'(a)') 'ERROR: specify one of --clump_tau0, --clump_NHI, --clump_nH, --taumax, --N_HImax'
    stop 1
  end subroutine

  !=======================================================================
  ! xorshift64* RNG
  !=======================================================================
  subroutine rng_init(seed)
    integer(int64), intent(in) :: seed
    if (seed == 0_int64) then
       rng_state = 88172645463325252_int64
    else
       rng_state = seed
    end if
  end subroutine

  function rng_u64() result(x)
    integer(int64) :: x
    x = rng_state
    x = ieor(x, ishft(x,  13))
    x = ieor(x, ishft(x, -7))
    x = ieor(x, ishft(x,  17))
    rng_state = x
    x = x * 2685821657736338717_int64
  end function

  function rng_uniform() result(u)
    real(wp) :: u
    integer(int64) :: x
    ! 53 high bits of x -> uniform in [0,1)
    x = rng_u64()
    x = iand(ishft(x, -11), int(z'001FFFFFFFFFFFFF', int64))
    u = real(x, wp) / 9007199254740992.0_wp    ! 2^53
  end function

  function rng_normal() result(g)
    real(wp) :: g
    real(wp) :: u1, u2
    u1 = rng_uniform();  if (u1 < 1.0e-300_wp) u1 = 1.0e-300_wp
    u2 = rng_uniform()
    g = sqrt(-2.0_wp * log(u1)) * cos(twopi * u2)
  end function

  !=======================================================================
  ! Position sampling
  !=======================================================================
  subroutine sample_one(x, y, z, r, ct, st, rcl)
    real(wp), intent(out) :: x, y, z, r, ct, st, rcl
    real(wp) :: u, r3, phi, sgn
    real(wp) :: r_min_eff, r_max_eff

    if (profiles_active) then
       u = rng_uniform()
       r = sample_r_from_cdf(u)
       rcl = base_radius * interp_prof(r, prof_shape_r)
    else
       if (clump_fully_inside) then
          r_min_eff = r_min_clump + base_radius
          r_max_eff = sphere_R - base_radius
       else
          r_min_eff = r_min_clump
          r_max_eff = sphere_R
       end if
       u = rng_uniform()
       r3 = r_min_eff**3 + (r_max_eff**3 - r_min_eff**3) * u
       r  = r3 ** (1.0_wp/3.0_wp)
       rcl = base_radius
    end if

    if (cos_cone > 0.0_wp) then
       ct = cos_cone + (1.0_wp - cos_cone) * rng_uniform()
       sgn = rng_uniform()
       if (sgn < 0.5_wp) ct = -ct
    else
       ct = 2.0_wp * rng_uniform() - 1.0_wp
    end if
    st  = sqrt(max(0.0_wp, 1.0_wp - ct*ct))
    phi = twopi * rng_uniform()
    x   = r * st * cos(phi)
    y   = r * st * sin(phi)
    z   = r * ct
  end subroutine

  ! Cone "fully inside" test
  function cone_inside_ok(r, ct, st, rcl) result(ok)
    real(wp), intent(in) :: r, ct, st, rcl
    logical :: ok
    real(wp) :: ratio, cos_min
    if (cos_cone <= 0.0_wp) then
       ok = .true.; return
    end if
    if (r <= rcl) then
       ok = .true.; return
    end if
    ratio   = rcl / r
    cos_min = abs(ct) * sqrt(max(0.0_wp, 1.0_wp - ratio*ratio)) - st * ratio
    ok = (cos_min >= cos_cone)
  end function

  !=======================================================================
  ! RSA placement with hash-grid + linked lists
  !=======================================================================
  subroutine rsa_place()
    integer(int64), allocatable :: head(:)
    integer(int64), allocatable :: nxt(:)
    integer :: rg
    real(wp) :: rg_cell, inv_cell
    integer :: ig, jg, kg, ii, jj, kk
    integer :: ig_lo, ig_hi, jg_lo, jg_hi, kg_lo, kg_hi
    integer(int64) :: idx, cidx
    integer(int64) :: placed, attempts
    integer(int64) :: max_attempts
    real(wp) :: x, y, z, rr, ct, st, rcl
    real(wp) :: dx, dy, dz, d2, sep
    real(wp) :: min_sep2_uniform
    logical  :: overlap, ok_geom
    real(wp) :: t0, t1, acc
    integer(int64) :: next_report
    integer(int32) :: tval
    integer :: rgm1, rg2

    rg = nint( real(N_target, wp) ** (1.0_wp/3.0_wp) ) + 1
    if (rg < 32) rg = 32
    if (rg > 512) rg = 512
    rg_cell = max(2.0_wp * sphere_R / real(rg, wp), 2.0_wp * r_cl_max)
    rg = int(2.0_wp * sphere_R / rg_cell) + 1
    if (rg < 2) rg = 2
    rg_cell = 2.0_wp * sphere_R / real(rg, wp)
    inv_cell = 1.0_wp / rg_cell
    rgm1 = rg - 1
    rg2  = rg * rg

    allocate(head(int(rg, int64) * int(rg, int64) * int(rg, int64)))
    allocate(nxt(N_target))
    head = -1_int64
    nxt  = -1_int64

    min_sep2_uniform = (2.0_wp * base_radius) ** 2

    write(*,'(a,i0,a,i0,a)') ' RSA: grid ', rg, '^3, placing ', N_target, ' clumps...'
    write(*,'(a,l1)')        ' RSA: clump_fully_inside = ', clump_fully_inside
    if (r_min_clump > 0.0_wp) &
         write(*,'(a,2f10.5)') ' RSA: shell rmin/rmax    = ', r_min_clump, sphere_R
    if (cos_cone > 0.0_wp) &
         write(*,'(a,f10.3,a)') ' RSA: cone_opening       = ', cone_opening, ' deg'

    call cpu_time(t0)

    placed = 0_int64
    attempts = 0_int64
    max_attempts = 200_int64 * max(N_target, 1_int64)
    next_report = 100000_int64

    do while (placed < N_target)
       attempts = attempts + 1_int64
       call sample_one(x, y, z, rr, ct, st, rcl)

       ! Geometric acceptance (fully-inside / cone-inside)
       ok_geom = .true.
       if (clump_fully_inside) then
          if (profiles_active) then
             if (rr + rcl > sphere_R)     ok_geom = .false.
             if (rr - rcl < r_min_clump)  ok_geom = .false.
          end if
          if (cos_cone > 0.0_wp) then
             if (.not. cone_inside_ok(rr, ct, st, rcl)) ok_geom = .false.
          end if
       end if
       if (.not. ok_geom) then
          if (attempts > max_attempts) then
             write(*,'(a,i0,a,i0,a)') 'ERROR: RSA stalled at ', placed, '/', N_target, &
                  '; loosen geometry or N.'
             stop 1
          end if
          cycle
       end if

       ! Hash cell
       ig = int((x + sphere_R) * inv_cell)
       if (ig < 0) ig = 0;  if (ig > rgm1) ig = rgm1
       jg = int((y + sphere_R) * inv_cell)
       if (jg < 0) jg = 0;  if (jg > rgm1) jg = rgm1
       kg = int((z + sphere_R) * inv_cell)
       if (kg < 0) kg = 0;  if (kg > rgm1) kg = rgm1

       overlap = .false.
       if (.not. clump_allow_overlap) then
          ig_lo = max(ig-1, 0);    ig_hi = min(ig+1, rgm1)
          jg_lo = max(jg-1, 0);    jg_hi = min(jg+1, rgm1)
          kg_lo = max(kg-1, 0);    kg_hi = min(kg+1, rgm1)
          outer: do kk = kg_lo, kg_hi
             do jj = jg_lo, jg_hi
                do ii = ig_lo, ig_hi
                   cidx = int(kk, int64) * int(rg2, int64) &
                        + int(jj, int64) * int(rg,  int64) &
                        + int(ii, int64) + 1_int64
                   idx = head(cidx)
                   do while (idx >= 1_int64)
                      dx = x - cl_x(idx)
                      dy = y - cl_y(idx)
                      dz = z - cl_z(idx)
                      d2 = dx*dx + dy*dy + dz*dz
                      if (profiles_active) then
                         sep = rcl + cl_r(idx)
                         if (d2 < sep*sep) then
                            overlap = .true.; exit outer
                         end if
                      else
                         if (d2 < min_sep2_uniform) then
                            overlap = .true.; exit outer
                         end if
                      end if
                      idx = nxt(idx)
                   end do
                end do
             end do
          end do outer
       end if

       if (overlap) then
          if (attempts > max_attempts) then
             write(*,'(a,i0,a,i0,a)') 'ERROR: RSA stalled at ', placed, '/', N_target, &
                  '; loosen geometry or N.'
             stop 1
          end if
          cycle
       end if

       placed = placed + 1_int64
       cl_x(placed) = x
       cl_y(placed) = y
       cl_z(placed) = z
       cl_r(placed) = rcl
       cidx = int(kg, int64) * int(rg2, int64) &
            + int(jg, int64) * int(rg,  int64) &
            + int(ig, int64) + 1_int64
       nxt(placed)  = head(cidx)
       head(cidx)   = placed

       if (placed >= next_report) then
          call cpu_time(t1)
          write(*,'(a,i10,a,i0,a,i0,a,es10.2,a)') &
               '   placed ', placed, ' / ', N_target, &
               '   (attempts=', attempts, ', ', &
               real(placed, wp) / max(t1 - t0, 1.0e-9_wp), ' clumps/s)'
          next_report = ((placed / 100000_int64) + 1_int64) * 100000_int64
       end if
    end do

    call cpu_time(t1)
    acc = 100.0_wp * real(N_target, wp) / max(real(attempts, wp), 1.0_wp)
    if (clump_allow_overlap) then
       write(*,'(a)') ' Random placement done (overlap allowed).'
    else
       write(*,'(a,f6.2,a)') ' RSA done, acceptance rate = ', acc, '%'
    end if

    N_placed = placed
    deallocate(head, nxt)
  end subroutine

  !=======================================================================
  ! Clump RHOKAP, TEMP (R_CLUMP already set during placement)
  !=======================================================================
  subroutine assign_perclump_physics()
    integer(int64) :: i
    real(wp) :: r, T, vtherm_i, Dfreq_i, dens

    if (.not. profiles_active) then
       cl_rhokap = rhokap_ref
       cl_temp   = temp_cl
       return
    end if

    do i = 1_int64, N_target
       r        = sqrt(cl_x(i)**2 + cl_y(i)**2 + cl_z(i)**2)
       cl_temp(i) = temp_cl
       T        = temp_cl
       vtherm_i = vtherm1 * sqrt(T)
       Dfreq_i  = vtherm_i / (wavelength0 * UM2KM)
       dens     = interp_prof(r, prof_shape_d)
       cl_rhokap(i) = rhokap_ref * dens * (Dfreq_ref / Dfreq_i)
    end do
  end subroutine

  !=======================================================================
  ! Velocity assignment
  !=======================================================================
  subroutine assign_velocities()
    integer(int64) :: i
    real(wp) :: x, y, z, rr, rcyl, scale, denom, sx, sy, sz
    character(len=32) :: vt
    vt = trim(adjustl(lowercase(velocity_type)))

    do i = 1_int64, N_target
       x = cl_x(i);  y = cl_y(i);  z = cl_z(i)
       rr = sqrt(x*x + y*y + z*z)
       rcyl = sqrt(x*x + y*y)
       sx = 0.0_wp; sy = 0.0_wp; sz = 0.0_wp

       select case (trim(vt))
       case ('', 'none')
          ! no flow
       case ('hubble')
          sx = Vexp * x / sphere_R
          sy = Vexp * y / sphere_R
          sz = Vexp * z / sphere_R
       case ('constant_radial')
          if (rr > 0.0_wp) then
             scale = Vexp / rr
             sx = scale * x; sy = scale * y; sz = scale * z
          end if
       case ('power_law')
          if (rr > 0.0_wp) then
             scale = Vexp * (rr / sphere_R) ** velocity_alpha
             sx = scale * x / rr
             sy = scale * y / rr
             sz = scale * z / rr
          end if
       case ('linear_decelerate')
          if (rr > 0.0_wp) then
             denom = sphere_R - max(r_min_clump, 0.0_wp)
             if (denom <= 0.0_wp) denom = 1.0_wp
             scale = Vexp * max(0.0_wp, (sphere_R - rr) / denom)
             sx = scale * x / rr
             sy = scale * y / rr
             sz = scale * z / rr
          end if
       case ('parallel_velocity')
          sx = Vxin;  sy = Vyin;  sz = Vzin
       case ('ssh')
          if (rr > 0.0_wp .and. rpeak > 0.0_wp) then
             if (rr < rpeak) then
                scale = Vpeak / rpeak
                sx = scale * x;  sy = scale * y;  sz = scale * z
             else
                if (sphere_R > rpeak) then
                   scale = Vpeak + DeltaV * (rr - rpeak) / (sphere_R - rpeak)
                else
                   scale = Vpeak
                end if
                sx = scale * x / rr
                sy = scale * y / rr
                sz = scale * z / rr
             end if
          end if
       case ('rotating_solid_body')
          sx = -Vrot * y / sphere_R
          sy =  Vrot * x / sphere_R
       case ('rotating_galaxy_halo')
          if (rcyl > 0.0_wp) then
             if (rcyl < rinner .and. rinner > 0.0_wp) then
                sx = -Vrot * y / rinner
                sy =  Vrot * x / rinner
             else
                sx = -Vrot * y / rcyl
                sy =  Vrot * x / rcyl
             end if
          end if
       case default
          write(*,'(a,a,a)') 'WARNING: unknown velocity_type ', trim(vt), &
                             '; no systematic flow applied.'
       end select

       if (clump_sigma_v > 0.0_wp) then
          sx = sx + clump_sigma_v * rng_normal()
          sy = sy + clump_sigma_v * rng_normal()
          sz = sz + clump_sigma_v * rng_normal()
       end if

       cl_vx(i) = sx;  cl_vy(i) = sy;  cl_vz(i) = sz
    end do
  end subroutine

  !=======================================================================
  ! Output via iofile_mod (matches clump_mod.f90:write_clumps_info schema)
  !=======================================================================
  subroutine write_output()
    type(io_file_type) :: iofh
    integer :: bitpix
    logical :: write_radius, write_rhokap, write_temp
    real(wp) :: rcl_min, rcl_max, rcl_mean
    real(wp) :: kap_min, kap_max, kap_mean
    real(wp) :: T_min, T_max, T_mean
    real(wp) :: f_vol_actual, f_cov_actual
    real(wp) :: cl_radius_max_out
    integer(int64) :: i, ncl

    ncl    = N_target
    bitpix = -32

    rcl_min  = huge(1.0_wp); rcl_max = 0.0_wp; rcl_mean = 0.0_wp
    kap_min  = huge(1.0_wp); kap_max = 0.0_wp; kap_mean = 0.0_wp
    T_min    = huge(1.0_wp); T_max   = 0.0_wp; T_mean   = 0.0_wp
    do i = 1_int64, ncl
       rcl_min  = min(rcl_min, cl_r(i));  rcl_max = max(rcl_max, cl_r(i))
       rcl_mean = rcl_mean + cl_r(i)
       kap_min  = min(kap_min, cl_rhokap(i));  kap_max = max(kap_max, cl_rhokap(i))
       kap_mean = kap_mean + cl_rhokap(i)
       T_min    = min(T_min, cl_temp(i));  T_max   = max(T_max, cl_temp(i))
       T_mean   = T_mean   + cl_temp(i)
    end do
    if (ncl > 0_int64) then
       rcl_mean = rcl_mean / real(ncl, wp)
       kap_mean = kap_mean / real(ncl, wp)
       T_mean   = T_mean   / real(ncl, wp)
    end if

    write_radius = (rcl_max - rcl_min) > CONST_TOL * max(abs(rcl_mean), tiny(1.0_wp))
    write_rhokap = (kap_max - kap_min) > CONST_TOL * max(abs(kap_mean), tiny(1.0_wp))
    write_temp   = (T_max   - T_min  ) > CONST_TOL * max(abs(T_mean  ), tiny(1.0_wp))

    cl_radius_max_out = rcl_max

    if (profiles_active) then
       f_vol_actual = 0.0_wp
       do i = 1_int64, ncl
          f_vol_actual = f_vol_actual + cl_r(i)**3
       end do
       f_vol_actual = f_vol_actual &
                    / max(sphere_R**3 - r_min_clump**3, tiny(1.0_wp))
       f_cov_actual = fcov_input_est
    else
       f_vol_actual = real(ncl, wp) * cl_radius_max_out**3 &
                    / max(sphere_R**3 - r_min_clump**3, tiny(1.0_wp))
       f_cov_actual = 0.75_wp * real(ncl, wp) * cl_radius_max_out**2 &
                    / max(sphere_R**2 + sphere_R*r_min_clump + r_min_clump**2, &
                          tiny(1.0_wp))
    end if

    status = 0
    call io_open_new(iofh, trim(output_file), status)
    if (status /= 0) then
       write(*,'(2a)') 'ERROR: cannot create output file ', trim(output_file)
       stop 1
    end if

    call io_write_table_column(iofh, 'X',  cl_x(1:ncl), status, bitpix)
    call io_write_table_column(iofh, 'Y',  cl_y(1:ncl), status, bitpix)
    call io_write_table_column(iofh, 'Z',  cl_z(1:ncl), status, bitpix)
    call io_write_table_column(iofh, 'VX', cl_vx(1:ncl), status, bitpix)
    call io_write_table_column(iofh, 'VY', cl_vy(1:ncl), status, bitpix)
    call io_write_table_column(iofh, 'VZ', cl_vz(1:ncl), status, bitpix)

    if (write_radius) &
         call io_write_table_column(iofh, 'R_CLUMP', cl_r(1:ncl),     status, bitpix)
    if (write_rhokap) &
         call io_write_table_column(iofh, 'RHOKAP',  cl_rhokap(1:ncl),status, bitpix)
    if (write_temp) &
         call io_write_table_column(iofh, 'TEMP',    cl_temp(1:ncl),  status, bitpix)

    call io_put_keyword(iofh, 'N_CLUMPS', ncl,                   'number of clumps (realized)',         status)
    call io_put_keyword(iofh, 'SPHERE_R', sphere_R,              'outer sphere radius [code units]',     status)
    call io_put_keyword(iofh, 'RMIN',     r_min_clump,           'inner placement radius [code units]',  status)
    call io_put_keyword(iofh, 'CL_RAD',   cl_radius_max_out,     'clump radius (max) [code units]',      status)
    call io_put_keyword(iofh, 'F_VOL',    f_vol_actual,          'volume filling factor (realized)',     status)
    call io_put_keyword(iofh, 'F_COV',    f_cov_actual,          'covering factor (realized)',           status)
    call io_put_keyword(iofh, 'TAU0',     clump_tau0,            'line-center tau (center to surface)',  status)
    call io_put_keyword(iofh, 'SIGMA_V',  clump_sigma_v,         'bulk velocity sigma [km/s]',           status)
    call io_put_keyword(iofh, 'TEMP_CL',  temp_cl,               'clump temperature [K] (reference)',    status)
    call io_put_keyword(iofh, 'RHOKAP',   rhokap_ref,            'opacity/code-unit inside clumps (ref)',status)
    call io_put_keyword(iofh, 'CL_DFREQ', Dfreq_ref,             'Doppler frequency [Hz] (reference)',   status)
    call io_put_keyword(iofh, 'VTHERM',   vtherm_ref,            'thermal velocity [km/s] (reference)',  status)
    call io_put_keyword(iofh, 'VOIGT_A',  voigt_a_ref,           'Voigt damping parameter (reference)',  status)
    call io_put_keyword(iofh, 'RMAX',     rmax_in,               'outer sphere radius input',            status)
    call io_put_keyword(iofh, 'IN_FCOV',  clump_f_cov,           'covering factor (input)',              status)
    call io_put_keyword(iofh, 'IN_FVOL',  clump_f_vol,           'volume filling factor (input)',        status)
    call io_put_keyword(iofh, 'IN_NCL',   clump_N_clumps,        'N_clumps (input)',                     status)
    call io_put_keyword(iofh, 'IN_NHI',   clump_NHI,             'clump NHI input [cm^-2]',          status)
    call io_put_keyword(iofh, 'IN_NH',    clump_nH,              'clump nH density input [cm^-3]',       status)
    call io_put_keyword(iofh, 'IN_TEMP',  clump_temperature,     'clump temperature input [K]',          status)
    call io_put_keyword(iofh, 'DISTUNIT', trim(distance_unit),   'Distance Unit',                        status)
    call io_put_keyword(iofh, 'DIST_CM',  distance2cm,           'Distance Unit (cm)',                   status)
    call io_put_keyword(iofh, 'LINE_ID',  trim(line_id),         'line atomic data',                     status)
    call io_put_keyword(iofh, 'CONE_OP',  cone_opening,          'bicone half-opening angle [deg]',      status)

    call io_close(iofh, status)
    if (status /= 0) then
       write(*,'(a)') 'ERROR: closing output file failed'
       stop 1
    end if

    write(*,'(a,3l2)') ' Clumps: clump columns saved (R_CLUMP/RHOKAP/TEMP) = ', &
         write_radius, write_rhokap, write_temp
    write(*,'(2a)')    ' Clumps saved to ', trim(output_file)
  end subroutine

  !=======================================================================
  ! Diagnostic summary
  !=======================================================================
  subroutine print_summary()
    integer(int64) :: i
    real(wp) :: rcl_min, rcl_max, rcl_mean, tau_mean
    real(wp) :: f_vol_actual, f_cov_actual

    rcl_min  = huge(1.0_wp); rcl_max = 0.0_wp; rcl_mean = 0.0_wp; tau_mean = 0.0_wp
    do i = 1_int64, N_target
       rcl_min  = min(rcl_min,  cl_r(i))
       rcl_max  = max(rcl_max,  cl_r(i))
       rcl_mean = rcl_mean + cl_r(i)
       tau_mean = tau_mean + cl_rhokap(i) * voigt(0.0_wp, voigt_a_ref) * cl_r(i)
    end do
    if (N_target > 0_int64) then
       rcl_mean = rcl_mean / real(N_target, wp)
       tau_mean = tau_mean / real(N_target, wp)
    end if

    if (.not. profiles_active) then
       f_vol_actual = real(N_target,wp) * base_radius**3 &
                    / max(sphere_R**3 - r_min_clump**3, tiny(1.0_wp))
       f_cov_actual = 0.75_wp * real(N_target,wp) * base_radius**2 &
                    / max(sphere_R**2 + sphere_R*r_min_clump + r_min_clump**2, &
                          tiny(1.0_wp))
    else
       f_vol_actual = 0.0_wp
       do i = 1_int64, N_target
          f_vol_actual = f_vol_actual + cl_r(i)**3
       end do
       f_vol_actual = f_vol_actual &
                    / max(sphere_R**3 - r_min_clump**3, tiny(1.0_wp))
       f_cov_actual = fcov_input_est
    end if

    write(*,'(a,es12.4)') ' Clumps: realized f_vol      = ', f_vol_actual
    write(*,'(a,es12.4)') ' Clumps: realized f_cov      = ', f_cov_actual
    write(*,'(a,3es12.4)')' Clumps: r_cl min/mean/max   = ', rcl_min, rcl_mean, rcl_max
    write(*,'(a,es12.4)') ' Clumps: mean tau_per_clump  = ', tau_mean
  end subroutine

  !=======================================================================
  ! CLI parser
  !=======================================================================
  subroutine validate_args()
    if (rmax_in <= 0.0_wp) then
       write(*,'(a)') 'ERROR: --rmax must be > 0'
       stop 1
    end if
    if (clump_radius <= 0.0_wp) then
       write(*,'(a)') 'ERROR: --clump_radius must be > 0'
       stop 1
    end if
    if (rmin_in >= rmax_in) then
       write(*,'(a)') 'ERROR: --rmin must be < --rmax'
       stop 1
    end if
    if (len_trim(output_file) == 0) then
       write(*,'(a)') 'ERROR: -o/--output is required'
       stop 1
    end if
  end subroutine

  subroutine parse_args()
    integer :: nargs, i
    character(len=256) :: arg, val

    nargs = command_argument_count()
    if (nargs == 0) then
       call print_usage(); stop
    end if
    i = 1
    do while (i <= nargs)
       call get_command_argument(i, arg)
       select case (trim(arg))
       case ('--rmax');          call get_next(i,val); read(val,*) rmax_in
       case ('--rmin');          call get_next(i,val); read(val,*) rmin_in
       case ('--cone_opening');  call get_next(i,val); read(val,*) cone_opening
       case ('--distance_unit'); call get_next(i,val); distance_unit = trim(val)

       case ('--clump_radius');   call get_next(i,val); read(val,*) clump_radius
       case ('--clump_N_clumps'); call get_next(i,val); read(val,*) clump_N_clumps
       case ('--clump_f_vol');    call get_next(i,val); read(val,*) clump_f_vol
       case ('--clump_f_cov');    call get_next(i,val); read(val,*) clump_f_cov

       case ('--clump_tau0'); call get_next(i,val); read(val,*) clump_tau0
       case ('--clump_NHI');  call get_next(i,val); read(val,*) clump_NHI
       case ('--clump_nH');   call get_next(i,val); read(val,*) clump_nH
       case ('--taumax');     call get_next(i,val); read(val,*) taumax
       case ('--N_HImax');    call get_next(i,val); read(val,*) N_HImax

       case ('--temperature');       call get_next(i,val); read(val,*) temperature
       case ('--clump_temperature'); call get_next(i,val); read(val,*) clump_temperature
       case ('--line_id');           call get_next(i,val); line_id = trim(val)

       case ('--clump_fully_inside');    clump_fully_inside  = .true.
       case ('--no-clump_fully_inside'); clump_fully_inside  = .false.
       case ('--clump_allow_overlap');   clump_allow_overlap = .true.

       case ('--velocity_type');  call get_next(i,val); velocity_type = trim(val)
       case ('--Vexp');           call get_next(i,val); read(val,*) Vexp
       case ('--Vrot');           call get_next(i,val); read(val,*) Vrot
       case ('--velocity_alpha'); call get_next(i,val); read(val,*) velocity_alpha
       case ('--rinner');         call get_next(i,val); read(val,*) rinner
       case ('--rpeak');          call get_next(i,val); read(val,*) rpeak
       case ('--Vpeak');          call get_next(i,val); read(val,*) Vpeak
       case ('--DeltaV');         call get_next(i,val); read(val,*) DeltaV
       case ('--Vx');             call get_next(i,val); read(val,*) Vxin
       case ('--Vy');             call get_next(i,val); read(val,*) Vyin
       case ('--Vz');             call get_next(i,val); read(val,*) Vzin
       case ('--clump_sigma_v');  call get_next(i,val); read(val,*) clump_sigma_v

       case ('--clump_radius_profile');  call get_next(i,val); clump_radius_profile  = trim(val)
       case ('--clump_density_profile'); call get_next(i,val); clump_density_profile = trim(val)
       case ('--clump_number_profile');  call get_next(i,val); clump_number_profile  = trim(val)
       case ('--clump_radius_alpha');    call get_next(i,val); read(val,*) clump_radius_alpha
       case ('--clump_density_alpha');   call get_next(i,val); read(val,*) clump_density_alpha
       case ('--clump_number_alpha');    call get_next(i,val); read(val,*) clump_number_alpha
       case ('--clump_radius_r0');       call get_next(i,val); read(val,*) clump_radius_r0
       case ('--clump_density_r0');      call get_next(i,val); read(val,*) clump_density_r0
       case ('--clump_number_r0');       call get_next(i,val); read(val,*) clump_number_r0

       case ('-o', '--output'); call get_next(i,val); output_file = trim(val)
       case ('--seed');         call get_next(i,val); read(val,*) seed_in

       case ('-h', '--help');   call print_usage(); stop

       case default
          write(*,'(2a)') 'Unknown argument: ', trim(arg)
          call print_usage(); stop 1
       end select
       i = i + 1
    end do
  end subroutine

  subroutine get_next(i, val)
    integer, intent(inout) :: i
    character(len=*), intent(out) :: val
    i = i + 1
    if (i > command_argument_count()) then
       write(*,'(a)') 'ERROR: missing value for last argument'
       stop 1
    end if
    call get_command_argument(i, val)
  end subroutine

  subroutine print_usage()
    write(*,'(a)') 'Usage: make_clumps.x [options]'
    write(*,'(a)') ''
    write(*,'(a)') 'Geometry:'
    write(*,'(a)') '  --rmax VAL                  outer sphere radius (REQUIRED)'
    write(*,'(a)') '  --rmin VAL                  inner placement radius (default 0)'
    write(*,'(a)') '  --cone_opening VAL          bicone half-opening angle [deg]'
    write(*,'(a)') '  --distance_unit U           cm|m|km|au|pc|kpc|Mpc|Rsun (for --clump_nH)'
    write(*,'(a)') ''
    write(*,'(a)') 'Clump size & count (one of N_clumps/f_vol/f_cov required):'
    write(*,'(a)') '  --clump_radius VAL          base clump radius (REQUIRED)'
    write(*,'(a)') '  --clump_N_clumps VAL        exact number of clumps'
    write(*,'(a)') '  --clump_f_vol VAL           volume filling factor'
    write(*,'(a)') '  --clump_f_cov VAL           covering factor (radial sightline)'
    write(*,'(a)') ''
    write(*,'(a)') 'Opacity (one of these required):'
    write(*,'(a)') '  --clump_tau0 VAL            clump line-center tau'
    write(*,'(a)') '  --clump_NHI VAL             clump HI column [cm^-2]'
    write(*,'(a)') '  --clump_nH VAL              HI density [cm^-3] (needs --distance_unit)'
    write(*,'(a)') '  --taumax VAL                system-level radial tau (back-solve)'
    write(*,'(a)') '  --N_HImax VAL               system-level HI column (back-solve)'
    write(*,'(a)') ''
    write(*,'(a)') 'Thermal & line:'
    write(*,'(a)') '  --temperature VAL           [K] (default 1e4)'
    write(*,'(a)') '  --clump_temperature VAL     override clump T (<0 = use --temperature)'
    write(*,'(a)') '  --line_id NAME              ly_alpha|CIV_1548|NV_1239|OVI_1032|'
    write(*,'(a)') '                              NaI_D|CaII_HK|MgII_2796|SiIV_1394|'
    write(*,'(a)') '                              AlII_1671|CII_1334'
    write(*,'(a)') ''
    write(*,'(a)') 'RSA controls:'
    write(*,'(a)') '  --clump_fully_inside        require clumps to lie entirely inside shell+cone'
    write(*,'(a)') '  --no-clump_fully_inside     allow clumps to protrude'
    write(*,'(a)') '  --clump_allow_overlap       skip overlap rejection'
    write(*,'(a)') ''
    write(*,'(a)') 'Systematic velocity:'
    write(*,'(a)') '  --velocity_type T           none|hubble|constant_radial|power_law|'
    write(*,'(a)') '                              linear_decelerate|parallel_velocity|ssh|'
    write(*,'(a)') '                              rotating_solid_body|rotating_galaxy_halo'
    write(*,'(a)') '  --Vexp VAL                  expansion speed [km/s]'
    write(*,'(a)') '  --Vrot VAL                  rotation speed [km/s]'
    write(*,'(a)') '  --velocity_alpha VAL        power_law exponent'
    write(*,'(a)') '  --rinner VAL                inner radius for rotating_galaxy_halo'
    write(*,'(a)') '  --rpeak VAL --Vpeak VAL --DeltaV VAL    SSH model parameters'
    write(*,'(a)') '  --Vx VAL --Vy VAL --Vz VAL              parallel_velocity components'
    write(*,'(a)') '  --clump_sigma_v VAL         Gaussian random sigma [km/s]'
    write(*,'(a)') ''
    write(*,'(a)') 'Radial profile knobs (constant|powerlaw|gaussian|exponential):'
    write(*,'(a)') '  --clump_radius_profile/--clump_radius_alpha/--clump_radius_r0'
    write(*,'(a)') '  --clump_density_profile/--clump_density_alpha/--clump_density_r0'
    write(*,'(a)') '  --clump_number_profile/--clump_number_alpha/--clump_number_r0'
    write(*,'(a)') ''
    write(*,'(a)') 'Output:'
    write(*,'(a)') '  -o, --output FILE           output filename (.fits, .fits.gz, .h5) (REQUIRED)'
    write(*,'(a)') '  --seed N                    RNG seed (default 12345)'
    write(*,'(a)') '  -h, --help                  show this help'
  end subroutine

end program make_clumps
