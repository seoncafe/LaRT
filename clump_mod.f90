module clump_mod
!---------------------------------------------------------------------------
! Clumpy medium support for LaRT_v2.00.
!
! Geometry: N_clumps spherical clumps of (per-clump) radius cl_radius(:)
! placed uniformly at random (non-overlapping via RSA) inside a sphere of
! radius sphere_R. Each clump has its own opacity cl_rhokap(:), temperature
! cl_temperature(:), Doppler frequency cl_Dfreq(:), Voigt parameter
! cl_voigt_a(:), and independent Gaussian random bulk velocities.
! In Phase 1 every entry is uniform; later phases populate the arrays
! from radial profiles.
!
! Algorithms:
!   Placement : linked-list RSA, O(N_cl) expected time for f_vol < 35%.
!   Raytrace  : DDA through CSR acceleration grid (Amanatides & Woo 1987).
!   Memory    : MPI-3 shared memory (one copy per node).
!---------------------------------------------------------------------------
  use define
  use memory_mod
  use voigt_mod
  use random
  use mpi
  implicit none
  public

  !--- Physical properties of the clump population
  integer(int64), save :: N_clumps   = 0_int64
  real(kind=wp),  save :: sphere_R   = 0.0_wp   ! outer sphere radius [code units]

  !--- Per-clump physical properties (MPI shared memory, dimension N_clumps).
  !    In Phase 1 every entry is uniform (filled with the corresponding _ref
  !    scalar below); later phases (radial profiles) populate them as
  !    functions of clump-centre radius.
  real(kind=wp), pointer, save :: cl_radius(:)      => null()  ! clump radius [code units]
  real(kind=wp), pointer, save :: cl_radius2(:)     => null()  ! cl_radius(icl)**2
  real(kind=wp), pointer, save :: cl_rhokap(:)      => null()  ! opacity/code-unit inside clump
  real(kind=wp), pointer, save :: cl_voigt_a(:)     => null()  ! Voigt damping parameter
  real(kind=wp), pointer, save :: cl_Dfreq(:)       => null()  ! Doppler frequency [Hz]
  real(kind=wp), pointer, save :: cl_vtherm(:)      => null()  ! thermal velocity [km/s]
  real(kind=wp), pointer, save :: cl_temperature(:) => null()  ! clump temperature [K]

  !--- Reference / representative scalars used during initialisation, grid
  !    setup, and FITS-header reporting. In Phase 1 these equal every
  !    cl_*(icl) entry. cl_radius_max is also used to size the RSA / CSR
  !    acceleration grids so that the 27-neighbour overlap search remains
  !    complete when r_cl varies between clumps.
  real(kind=wp), save :: cl_radius_max      = 0.0_wp
  real(kind=wp), save :: cl_rhokap_ref      = 0.0_wp
  real(kind=wp), save :: cl_voigt_a_ref     = 0.0_wp
  real(kind=wp), save :: cl_Dfreq_ref       = 0.0_wp
  real(kind=wp), save :: cl_vtherm_ref      = 0.0_wp
  real(kind=wp), save :: cl_temperature_ref = 0.0_wp

  !--- Clump positions [code units] and bulk velocities – MPI shared memory.
  !
  !    cl_vx/y/z are stored DIMENSIONLESSLY as v / cl_vtherm(icl) (matching the
  !    Cartesian/AMR grid%vfx convention). Each setter routine
  !    (generate_clumps and assign_clump_velocities_from_type) divides by
  !    cl_vtherm(icl) inline before writing, so no separate post-pass rescaling
  !    is needed. write_clumps_fits() multiplies by cl_vtherm(icl) to keep the
  !    user-facing FITS output in km/s.
  real(kind=dp), pointer, save :: cl_x(:)  => null()
  real(kind=dp), pointer, save :: cl_y(:)  => null()
  real(kind=dp), pointer, save :: cl_z(:)  => null()
  real(kind=dp), pointer, save :: cl_vx(:) => null()
  real(kind=dp), pointer, save :: cl_vy(:) => null()
  real(kind=dp), pointer, save :: cl_vz(:) => null()

  !--- Radial-profile machinery (Phase 3+4+5).
  !    base_* values are the user-supplied peak values that the shape factors
  !    multiply. r_cl(r) = base_radius_in * shape_radius(r), etc.
  !    Profiles defined on r in [0, sphere_R].
  real(kind=wp), save :: base_radius_in   = 0.0_wp     ! par%clump_radius
  real(kind=wp), save :: base_rhokap_in   = 0.0_wp     ! density-equivalent peak rhokap
  real(kind=wp), save :: base_nH_in       = 0.0_wp     ! peak n_H [cm^-3] (when known)
  logical,       save :: profiles_active  = .false.    ! .true. if any profile != 'constant'
  logical,       save :: clumps_from_file = .false.    ! .true. if init_clumps loaded the population
                                                       !  from par%clump_input_file (via read_clumps_fits)

  !--- Tabulated radial CDF for inverse-CDF sampling of clump positions.
  integer, parameter :: NPROF = 4001                   ! 1-D radial table size
  real(kind=wp), save :: prof_dr     = 0.0_wp
  real(kind=wp), save :: prof_r(NPROF)
  real(kind=wp), save :: prof_shape_number(NPROF)      ! n_cl(r) shape (unnormalised)
  real(kind=wp), save :: prof_shape_radius(NPROF)      ! r_cl(r) shape
  real(kind=wp), save :: prof_shape_density(NPROF)     ! n_H(r) shape (rhokap shape)
  real(kind=wp), save :: prof_cdf_pos(NPROF)           ! CDF for sampling r ~ shape_number * r^2

  !--- Optional tabulated profile from `clump_profile_file`. When loaded,
  !    table_r/r_cl_in_table/n_H_in_table/T_in_table/n_cl_shape_in_table are
  !    interpolated linearly. Per-axis profile selection still uses the
  !    string flag; 'file' on any axis pulls from the corresponding column.
  integer, save :: NTAB = 0
  real(kind=wp), allocatable, save :: tab_r(:)
  real(kind=wp), allocatable, save :: tab_radius(:)
  real(kind=wp), allocatable, save :: tab_density(:)
  real(kind=wp), allocatable, save :: tab_temperature(:)
  real(kind=wp), allocatable, save :: tab_number(:)

  !--- CSR acceleration grid (MPI shared memory)
  !    Cell (i,j,k), 0-based → 1-based index = 1 + i + j*cgx + k*cgx*cgy
  !    cg_start(icell) .. cg_start(icell+1)-1 → entries in cg_list for that cell
  integer,       save :: cgx = 0, cgy = 0, cgz = 0
  real(kind=wp), save :: cg_xmin, cg_ymin, cg_zmin
  real(kind=wp), save :: cg_dx, cg_dy, cg_dz
  real(kind=wp), save :: cg_inv_dx, cg_inv_dy, cg_inv_dz
  integer(int32), pointer, save :: cg_start(:) => null()   ! size ncells+1
  integer(int32), pointer, save :: cg_list(:)  => null()   ! size total_registrations

contains

  !===========================================================================
  pure integer function cg_cell_idx(i, j, k)
  integer, intent(in) :: i, j, k
  cg_cell_idx = 1 + i + j*cgx + k*cgx*cgy
  end function cg_cell_idx
  !===========================================================================

  !===========================================================================
  ! Radial-profile evaluator. Returns the multiplicative shape factor
  ! associated with a profile name at radius r (in code units, 0 <= r <= R_box).
  ! The shape is normalised so that the user-supplied base value
  ! (par%clump_radius / par%clump_tau0 / etc.) is multiplied by this factor
  ! to obtain the local quantity. Built-in shapes:
  !   'constant'    : f = 1
  !   'powerlaw'    : f = (max(r, r_floor) / max(r0, r_floor)) ** (-alpha)
  !                   r_floor = 0.05 * r0 prevents divergence at r=0.
  !   'gaussian'    : f = exp(-(r/r0)**2)
  !   'exponential' : f = exp(-r/r0)
  !   'file'        : linear interpolation from the relevant tabulated column.
  !
  ! `axis_id` selects which file column when shape_name == 'file':
  !   1 = radius, 2 = density, 3 = number-density (shape only), 4 = temperature
  !===========================================================================
  pure real(kind=wp) function profile_shape(shape_name, alpha, r0, r, axis_id) &
       result(f)
  character(len=*), intent(in) :: shape_name
  real(kind=wp),    intent(in) :: alpha, r0, r
  integer,          intent(in) :: axis_id
  real(kind=wp) :: r_eff, r0_eff, r_floor

  select case (trim(shape_name))
  case ('constant', '')
     f = 1.0_wp
  case ('powerlaw')
     if (r0 > 0.0_wp) then
        r_floor = 0.05_wp * r0
        r_eff   = max(r,  r_floor)
        r0_eff  = max(r0, r_floor)
        f = (r_eff / r0_eff) ** (-alpha)
     else
        f = 1.0_wp
     end if
  case ('gaussian')
     if (r0 > 0.0_wp) then
        f = exp(-(r/r0)**2)
     else
        f = 1.0_wp
     end if
  case ('exponential')
     if (r0 > 0.0_wp) then
        f = exp(-r/r0)
     else
        f = 1.0_wp
     end if
  case ('file')
     f = profile_file_interp(r, axis_id)
  case default
     f = 1.0_wp
  end select
  end function profile_shape
  !===========================================================================

  !===========================================================================
  ! Linear interpolation of the tabulated `clump_profile_file`.
  ! Returns 1.0 if the table is not loaded or axis column is missing.
  !===========================================================================
  pure real(kind=wp) function profile_file_interp(r, axis_id) result(f)
  real(kind=wp), intent(in) :: r
  integer,       intent(in) :: axis_id
  integer :: lo, hi, mid
  real(kind=wp) :: t

  if (NTAB <= 1) then
     f = 1.0_wp;  return
  end if

  if (r <= tab_r(1)) then
     lo = 1;  hi = 2
  else if (r >= tab_r(NTAB)) then
     lo = NTAB-1;  hi = NTAB
  else
     lo = 1;  hi = NTAB
     do while (hi - lo > 1)
        mid = (lo + hi) / 2
        if (r >= tab_r(mid)) then
           lo = mid
        else
           hi = mid
        end if
     end do
  end if

  t = (r - tab_r(lo)) / max(tab_r(hi) - tab_r(lo), tiny(1.0_wp))
  t = max(0.0_wp, min(1.0_wp, t))

  select case (axis_id)
  case (1)
     if (allocated(tab_radius)) then
        f = (1.0_wp - t) * tab_radius(lo) + t * tab_radius(hi)
     else
        f = 1.0_wp
     end if
  case (2)
     if (allocated(tab_density)) then
        f = (1.0_wp - t) * tab_density(lo) + t * tab_density(hi)
     else
        f = 1.0_wp
     end if
  case (3)
     if (allocated(tab_number)) then
        f = (1.0_wp - t) * tab_number(lo) + t * tab_number(hi)
     else
        f = 1.0_wp
     end if
  case (4)
     if (allocated(tab_temperature)) then
        f = (1.0_wp - t) * tab_temperature(lo) + t * tab_temperature(hi)
     else
        f = 1.0_wp
     end if
  case default
     f = 1.0_wp
  end select
  end function profile_file_interp
  !===========================================================================

  !===========================================================================
  ! Convenience wrappers selecting the user-specified profile per axis.
  !===========================================================================
  pure real(kind=wp) function shape_radius(r) result(f)
  real(kind=wp), intent(in) :: r
  f = profile_shape(par%clump_radius_profile,  par%clump_radius_alpha, &
                    par%clump_radius_r0, r, 1)
  end function shape_radius

  pure real(kind=wp) function shape_density(r) result(f)
  real(kind=wp), intent(in) :: r
  f = profile_shape(par%clump_density_profile, par%clump_density_alpha, &
                    par%clump_density_r0, r, 2)
  end function shape_density

  pure real(kind=wp) function shape_number(r) result(f)
  real(kind=wp), intent(in) :: r
  f = profile_shape(par%clump_number_profile,  par%clump_number_alpha, &
                    par%clump_number_r0, r, 3)
  end function shape_number
  !===========================================================================

  !===========================================================================
  ! Build the radial CDF used for inverse-CDF sampling of clump positions.
  ! P(r) ∝ shape_number(r) * r^2  on [0, sphere_R].
  ! Trapezoidal integration over NPROF=4001 points; the tabulated CDF is
  ! used to convert a uniform U(0,1) draw into a radius via linear interp.
  !
  ! Sets module-level prof_r, prof_shape_*, prof_cdf_pos, prof_dr.
  ! Also computes cl_radius_max as the maximum of base_radius_in *
  ! shape_radius(r) over the radial table — needed for grid sizing before
  ! clumps are placed.
  !===========================================================================
  subroutine build_radial_profile_tables(R_box)
  real(kind=wp), intent(in) :: R_box
  integer        :: i
  real(kind=wp)  :: r, integrand, last_integrand, total, rcl_local

  prof_dr = R_box / real(NPROF - 1, wp)
  cl_radius_max = base_radius_in   ! starting estimate: base value

  prof_r(1)             = 0.0_wp
  prof_shape_number(1)  = shape_number(0.0_wp)
  prof_shape_radius(1)  = shape_radius(0.0_wp)
  prof_shape_density(1) = shape_density(0.0_wp)
  prof_cdf_pos(1)       = 0.0_wp

  rcl_local = base_radius_in * prof_shape_radius(1)
  if (rcl_local > cl_radius_max) cl_radius_max = rcl_local

  last_integrand = prof_shape_number(1) * 0.0_wp     ! r^2 = 0 at r=0

  do i = 2, NPROF
     r = real(i-1, wp) * prof_dr
     prof_r(i)             = r
     prof_shape_number(i)  = shape_number(r)
     prof_shape_radius(i)  = shape_radius(r)
     prof_shape_density(i) = shape_density(r)
     integrand             = prof_shape_number(i) * r * r
     prof_cdf_pos(i) = prof_cdf_pos(i-1) + 0.5_wp * (last_integrand + integrand) * prof_dr
     last_integrand = integrand

     rcl_local = base_radius_in * prof_shape_radius(i)
     if (rcl_local > cl_radius_max) cl_radius_max = rcl_local
  end do

  total = prof_cdf_pos(NPROF)
  if (total > 0.0_wp) then
     prof_cdf_pos(:) = prof_cdf_pos(:) / total
  else
     ! Degenerate (shape ≡ 0): fall back to uniform-in-volume sampling
     do i = 1, NPROF
        r = real(i-1, wp) * prof_dr
        prof_cdf_pos(i) = (r / R_box) ** 3
     end do
  end if
  end subroutine build_radial_profile_tables
  !===========================================================================

  !===========================================================================
  ! Sample a clump-centre radius from the tabulated inverse CDF.
  !===========================================================================
  real(kind=wp) function sample_clump_radius() result(r_out)
  real(kind=wp) :: u, t
  integer       :: lo, hi, mid

  u = rand_number()
  if (u <= prof_cdf_pos(1)) then
     r_out = prof_r(1);  return
  end if
  if (u >= prof_cdf_pos(NPROF)) then
     r_out = prof_r(NPROF);  return
  end if

  lo = 1;  hi = NPROF
  do while (hi - lo > 1)
     mid = (lo + hi) / 2
     if (u >= prof_cdf_pos(mid)) then
        lo = mid
     else
        hi = mid
     end if
  end do

  t = (u - prof_cdf_pos(lo)) / max(prof_cdf_pos(hi) - prof_cdf_pos(lo), tiny(1.0_wp))
  t = max(0.0_wp, min(1.0_wp, t))
  r_out = (1.0_wp - t) * prof_r(lo) + t * prof_r(hi)
  end function sample_clump_radius
  !===========================================================================

  !===========================================================================
  ! Trapezoidal integration over the radial table.
  !===========================================================================
  pure real(kind=wp) function integrate_table(integrand) result(total)
  real(kind=wp), intent(in) :: integrand(NPROF)
  integer :: i
  total = 0.0_wp
  do i = 2, NPROF
     total = total + 0.5_wp * (integrand(i-1) + integrand(i)) * prof_dr
  end do
  end function integrate_table
  !===========================================================================

  !===========================================================================
  ! Numerical f_cov via the radial line-of-sight integral (LaRT convention,
  ! see docs/covering_factor_definitions.tex). f_cov is the expected number
  ! of clump intersections along a radial sightline from the centre to the
  ! sphere surface, evaluated for the *current* normalisation A_norm:
  !
  !     f_cov = ∫₀ᴿ A_norm * shape_number(r) * π * r_cl(r)^2 dr
  !
  ! For uniform shapes this reduces to the closed form (3/4) N (r_cl/R)^2.
  !===========================================================================
  pure real(kind=wp) function f_cov_LOS_quad(A_norm) result(fcov)
  real(kind=wp), intent(in) :: A_norm
  real(kind=wp) :: rcl, integrand(NPROF)
  integer :: i
  do i = 1, NPROF
     rcl = base_radius_in * prof_shape_radius(i)
     integrand(i) = A_norm * prof_shape_number(i) * pi * rcl * rcl
  end do
  fcov = integrate_table(integrand)
  end function f_cov_LOS_quad
  !===========================================================================

  !===========================================================================
  ! Numerical volume filling factor:
  !     f_vol = (4π/R³) * ∫₀ᴿ A_norm * shape_number(r) * r_cl(r)³ * r² dr
  ! Derived from f_vol = (sum clump volumes)/V_box with the spatial number
  ! density A * shape(r) and per-clump volume (4/3)π r_cl(r)³ integrated
  ! over the shell 4π r² dr. For uniform shapes this reduces to
  ! N (r_cl/R)³.
  !===========================================================================
  pure real(kind=wp) function f_vol_quad(A_norm, R_box) result(fvol)
  real(kind=wp), intent(in) :: A_norm, R_box
  real(kind=wp) :: rcl, r, integrand(NPROF)
  integer :: i
  do i = 1, NPROF
     r   = prof_r(i)
     rcl = base_radius_in * prof_shape_radius(i)
     integrand(i) = A_norm * prof_shape_number(i) * rcl**3 * r * r
  end do
  fvol = fourpi * integrate_table(integrand) / (R_box**3)
  end function f_vol_quad
  !===========================================================================

  !===========================================================================
  ! Total clump count for normalisation A_norm:
  !     N = ∫₀ᴿ A_norm * shape_number(r) * 4π r^2 dr
  !===========================================================================
  pure real(kind=wp) function total_count_quad(A_norm) result(Ntot)
  real(kind=wp), intent(in) :: A_norm
  real(kind=wp) :: integrand(NPROF), r
  integer :: i
  do i = 1, NPROF
     r = prof_r(i)
     integrand(i) = A_norm * prof_shape_number(i) * fourpi * r * r
  end do
  Ntot = integrate_table(integrand)
  end function total_count_quad
  !===========================================================================

  !===========================================================================
  ! Load tabulated radial profile from `par%clump_profile_file`.
  ! Expected ASCII format (whitespace-separated, '#' = comment):
  !
  !   # r [code units]   r_cl_factor   n_H_factor   T [K]   n_cl_shape
  !   <data rows ...>
  !
  ! All five columns are required. Set to 1.0 (factors) or a constant (T) on
  ! axes where no variation is desired. The reader stores the columns so
  ! that profile_file_interp() can do a linear lookup. The radii must be
  ! monotonically increasing.
  !===========================================================================
  subroutine load_clump_profile_file(R_box)
  real(kind=wp), intent(in) :: R_box
  integer :: u, ios, nread, count_lines, ierr
  character(len=512) :: linebuf
  real(kind=wp) :: r_v, rc_v, nH_v, T_v, nc_v

  ! Load on every rank so that the table is available wherever
  ! profile_file_interp() is called.
  if (len_trim(par%clump_profile_file) == 0) return

  open(newunit=u, file=trim(par%clump_profile_file), status='old', &
       action='read', iostat=ios)
  if (ios /= 0) then
     write(*,*) 'ERROR: cannot open clump_profile_file: ', &
                trim(par%clump_profile_file)
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end if

  ! Pass 1: count valid data lines
  count_lines = 0
  do
     read(u, '(a)', iostat=ios) linebuf
     if (ios /= 0) exit
     linebuf = adjustl(linebuf)
     if (len_trim(linebuf) == 0) cycle
     if (linebuf(1:1) == '#') cycle
     count_lines = count_lines + 1
  end do
  rewind(u)

  if (count_lines < 2) then
     write(*,*) 'ERROR: clump_profile_file needs at least 2 data rows'
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end if

  NTAB = count_lines
  if (allocated(tab_r))           deallocate(tab_r)
  if (allocated(tab_radius))      deallocate(tab_radius)
  if (allocated(tab_density))     deallocate(tab_density)
  if (allocated(tab_temperature)) deallocate(tab_temperature)
  if (allocated(tab_number))      deallocate(tab_number)
  allocate(tab_r(NTAB), tab_radius(NTAB), tab_density(NTAB), &
           tab_temperature(NTAB), tab_number(NTAB))

  nread = 0
  do
     read(u, '(a)', iostat=ios) linebuf
     if (ios /= 0) exit
     linebuf = adjustl(linebuf)
     if (len_trim(linebuf) == 0) cycle
     if (linebuf(1:1) == '#') cycle
     read(linebuf, *, iostat=ios) r_v, rc_v, nH_v, T_v, nc_v
     if (ios /= 0) cycle
     nread = nread + 1
     if (nread > NTAB) exit
     tab_r(nread)           = r_v
     tab_radius(nread)      = rc_v
     tab_density(nread)     = nH_v
     tab_temperature(nread) = T_v
     tab_number(nread)      = nc_v
  end do
  close(u)

  if (nread /= NTAB) then
     write(*,*) 'WARNING: clump_profile_file: read ', nread, ' rows, expected ', NTAB
     NTAB = nread
  end if

  if (mpar%p_rank == 0) then
     write(*,'(a,i6,a)') ' Clumps: loaded ', NTAB, &
                         ' rows from clump_profile_file'
     write(*,'(a,2es12.4)') ' Clumps: profile r range = ', tab_r(1), tab_r(NTAB)
     if (tab_r(NTAB) < R_box) write(*,'(a)') &
        '   NOTE: profile radius range does not cover full sphere; ' // &
        'extrapolating at edges.'
  end if
  end subroutine load_clump_profile_file
  !===========================================================================

  !===========================================================================
  subroutine init_clumps(R_sphere)
  !---------------------------------------------------------------------------
  ! Derive N_clumps, compute physical parameters, allocate MPI shared memory,
  ! run RSA placement, and build the CSR acceleration grid.
  !---------------------------------------------------------------------------
  implicit none
  real(kind=wp), intent(in) :: R_sphere
  real(kind=wp) :: temp_cl, vtherm, A_norm, fvol_unit, fcov_unit, fvol_realised, fcov_realised
  integer       :: ierr

  sphere_R       = R_sphere
  base_radius_in = par%clump_radius
  if (base_radius_in <= 0.0_wp) then
     if (mpar%p_rank == 0) write(*,*) 'ERROR: clump_radius must be > 0'
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end if

  !--- If the user supplied a pre-built clump FITS file, load it and skip the
  !    internal profile / RSA generation entirely. Sets cl_*, sphere_R,
  !    cl_radius_max and reference scalars; builds the CSR acceleration grid.
  if (len_trim(par%clump_input_file) > 0) then
     call read_clumps_fits(trim(par%clump_input_file), R_sphere)
     return
  end if

  !--- Phase 2/3: detect any active radial profile.
  profiles_active = (trim(par%clump_radius_profile)  /= 'constant') .or. &
                    (trim(par%clump_density_profile) /= 'constant') .or. &
                    (trim(par%clump_number_profile)  /= 'constant')

  !--- If 'file' is requested on any axis, load the tabulated profile.
  if (profiles_active .and. ( &
       trim(par%clump_radius_profile)  == 'file' .or. &
       trim(par%clump_density_profile) == 'file' .or. &
       trim(par%clump_number_profile)  == 'file')) then
     call load_clump_profile_file(R_sphere)
  end if

  !--- Build radial tables. For uniform case we still set cl_radius_max
  !    from base_radius_in and skip the table to save a small amount of work.
  if (profiles_active) then
     call build_radial_profile_tables(R_sphere)
  else
     cl_radius_max = base_radius_in
  end if

  !--- clump temperature → Doppler frequency and Voigt parameter (peak / reference)
  temp_cl = par%clump_temperature
  if (temp_cl < 0.0_wp) temp_cl = par%temperature
  vtherm             = line%vtherm1 * sqrt(temp_cl)
  cl_Dfreq_ref       = vtherm / (line%wavelength0 * um2km)
  cl_vtherm_ref      = vtherm                          ! km/s
  cl_temperature_ref = temp_cl                         ! [K]
  cl_voigt_a_ref     = (line%damping / fourpi) / cl_Dfreq_ref

  !--- clump opacity (peak value) from tau0, NHI (column density), or nH.
  !    Whichever input is given fixes base_rhokap_in. When density_profile
  !    is non-constant, this acts as the peak and shape_density(r)
  !    multiplies it.
  base_nH_in = 0.0_wp
  if (par%clump_tau0 > 0.0_wp) then
     cl_rhokap_ref = par%clump_tau0 / (voigt(0.0_wp, cl_voigt_a_ref) * base_radius_in)
  else if (par%clump_NHI > 0.0_wp) then
     ! clump_NHI = per-clump column density [cm^-2] from clump centre to surface
     ! (peak value; shape_density(r) modulates per clump).
     cl_rhokap_ref = par%clump_NHI * line%cross0 / (cl_Dfreq_ref * base_radius_in)
  else if (par%clump_nH > 0.0_wp) then
     if (par%distance2cm <= 0.0_wp) then
        if (mpar%p_rank == 0) write(*,*) 'ERROR: clump_nH requires distance_unit'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
     cl_rhokap_ref = par%clump_nH * line%cross0 * par%distance2cm / cl_Dfreq_ref
     base_nH_in    = par%clump_nH
  else
     if (mpar%p_rank == 0) write(*,*) 'ERROR: specify clump_tau0, clump_NHI, or clump_nH'
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end if
  base_rhokap_in = cl_rhokap_ref

  !--- derive N_clumps
  !    Uniform case: closed-form (preserves Phase 1 behaviour byte-for-byte).
  !    Profile case: numerical quadrature with the LOS f_cov definition
  !                  (see docs/covering_factor_definitions.tex).
  if (.not. profiles_active) then
     if (par%clump_N_clumps > 0.0_wp) then
        N_clumps = int(par%clump_N_clumps, int64)
     else if (par%clump_f_vol > 0.0_wp) then
        N_clumps = nint(par%clump_f_vol * (R_sphere/cl_radius_max)**3, int64)
     else if (par%clump_f_cov > 0.0_wp) then
        N_clumps = nint((4.0_wp/3.0_wp)*par%clump_f_cov*(R_sphere/cl_radius_max)**2, int64)
     else
        if (mpar%p_rank == 0) &
           write(*,*) 'ERROR: specify clump_N_clumps, clump_f_vol, or clump_f_cov'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
     fvol_realised = real(N_clumps, wp) * (cl_radius_max/R_sphere)**3
     fcov_realised = 0.75_wp * real(N_clumps, wp) * (cl_radius_max/R_sphere)**2
  else
     !--- non-uniform: integrate the shape and back-solve for the density
     !    normalisation A_norm of n_cl(r) = A_norm * shape_number(r).
     if (par%clump_N_clumps > 0.0_wp) then
        N_clumps = int(par%clump_N_clumps, int64)
        A_norm   = real(N_clumps, wp) / max(total_count_quad(1.0_wp), tiny(1.0_wp))
     else if (par%clump_f_vol > 0.0_wp) then
        fvol_unit = f_vol_quad(1.0_wp, R_sphere)
        A_norm    = par%clump_f_vol / max(fvol_unit, tiny(1.0_wp))
        N_clumps  = nint(total_count_quad(A_norm), int64)
     else if (par%clump_f_cov > 0.0_wp) then
        fcov_unit = f_cov_LOS_quad(1.0_wp)
        A_norm    = par%clump_f_cov / max(fcov_unit, tiny(1.0_wp))
        N_clumps  = nint(total_count_quad(A_norm), int64)
     else
        if (mpar%p_rank == 0) &
           write(*,*) 'ERROR: specify clump_N_clumps, clump_f_vol, or clump_f_cov'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
     fvol_realised = f_vol_quad(A_norm, R_sphere)
     fcov_realised = f_cov_LOS_quad(A_norm)
  end if
  if (N_clumps <= 0_int64) N_clumps = 1_int64

  if (mpar%p_rank == 0) then
     write(*,'(a,i14)')   ' Clumps: N_clumps  = ', N_clumps
     write(*,'(a,f12.6)') ' Clumps: f_vol     = ', fvol_realised
     write(*,'(a,f12.5)') ' Clumps: f_cov     = ', fcov_realised
     write(*,'(a,es12.4)') ' Clumps: cl_rhokap = ', cl_rhokap_ref
     write(*,'(a,f12.5)')  ' Clumps: voigt_a   = ', cl_voigt_a_ref
     write(*,'(a,es12.4)') ' Clumps: cl_Dfreq  = ', cl_Dfreq_ref
     if (profiles_active) then
        write(*,'(a,a)') ' Clumps: radius profile  = ', trim(par%clump_radius_profile)
        write(*,'(a,a)') ' Clumps: density profile = ', trim(par%clump_density_profile)
        write(*,'(a,a)') ' Clumps: number profile  = ', trim(par%clump_number_profile)
        write(*,'(a,f12.6,a,f12.6)') ' Clumps: r_cl range      = ', &
              base_radius_in*minval(prof_shape_radius), ' to ', cl_radius_max
     end if
  end if

  !--- shared memory for clump data
  call create_shared_mem(cl_x,  [int(N_clumps)])
  call create_shared_mem(cl_y,  [int(N_clumps)])
  call create_shared_mem(cl_z,  [int(N_clumps)])
  call create_shared_mem(cl_vx, [int(N_clumps)])
  call create_shared_mem(cl_vy, [int(N_clumps)])
  call create_shared_mem(cl_vz, [int(N_clumps)])

  !--- shared memory for per-clump physical properties
  call create_shared_mem(cl_radius,      [int(N_clumps)])
  call create_shared_mem(cl_radius2,     [int(N_clumps)])
  call create_shared_mem(cl_rhokap,      [int(N_clumps)])
  call create_shared_mem(cl_voigt_a,     [int(N_clumps)])
  call create_shared_mem(cl_Dfreq,       [int(N_clumps)])
  call create_shared_mem(cl_vtherm,      [int(N_clumps)])
  call create_shared_mem(cl_temperature, [int(N_clumps)])

  !--- Phase 1: fill per-clump arrays uniformly with the reference values
  !    so that behaviour is identical to the previous scalar implementation.
  !    Phase 5 will replace this loop with profile-driven assignments after
  !    clump positions are placed.
  if (mpar%h_rank == 0) then
     cl_radius(:)      = cl_radius_max
     cl_radius2(:)     = cl_radius_max * cl_radius_max
     cl_rhokap(:)      = cl_rhokap_ref
     cl_voigt_a(:)     = cl_voigt_a_ref
     cl_Dfreq(:)       = cl_Dfreq_ref
     cl_vtherm(:)      = cl_vtherm_ref
     cl_temperature(:) = cl_temperature_ref
  end if
  call MPI_BARRIER(mpar%hostcomm, ierr)

  !--- p_rank=0 generates positions; broadcast to all h_rank=0 (one per node)
  !--- so every node has an identical clump layout in its shared memory.
  if (mpar%h_rank == 0) then
     if (mpar%p_rank == 0) then
        call generate_clumps()
        if (len_trim(par%velocity_type) > 0) call assign_clump_velocities_from_type()
     end if
     call MPI_BCAST(cl_x,  int(N_clumps), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_y,  int(N_clumps), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_z,  int(N_clumps), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_vx, int(N_clumps), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_vy, int(N_clumps), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_vz, int(N_clumps), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     ! cl_vx/y/z are already stored as v / cl_vtherm(icl) (normalised inline
     ! by generate_clumps and assign_clump_velocities_from_type).
  end if
  call MPI_BARRIER(mpar%hostcomm, ierr)

  !--- h_rank=0 builds CSR grid; barrier before use
  call build_clump_csr()
  call MPI_BARRIER(mpar%hostcomm, ierr)

  end subroutine init_clumps
  !===========================================================================

  !===========================================================================
  subroutine generate_clumps()
  !---------------------------------------------------------------------------
  ! RSA with linked-list grid acceleration. Called only on h_rank=0.
  ! Two sampling paths:
  !  - profiles_active = .false.: uniform-in-sphere box rejection (legacy).
  !  - profiles_active = .true. : inverse-CDF sampling on r, isotropic angles,
  !    per-clump radius from shape_radius(r), per-clump opacity / temperature
  !    / Voigt parameter from the active profiles.
  ! Per-pair RSA overlap test handles non-uniform r_cl correctly.
  !---------------------------------------------------------------------------
  implicit none
  integer        :: rg, ncells_rsa
  real(kind=wp)  :: rg_cell, min_sep2_uniform
  real(kind=wp)  :: xc, yc, zc, dx, dy, dz, d2, sep_pair
  real(kind=wp)  :: r_trial, rcl_trial, cos_theta, sin_theta, phi_az
  real(kind=wp)  :: temp_loc, vth_loc, Df_loc, va_loc, kap_loc, dens_factor
  integer        :: ig, jg, kg, ig2, jg2, kg2, icell_rsa, jnb
  integer(int64) :: icl, n_attempts
  logical        :: overlap
  integer, allocatable :: head(:), nxt(:)

  rg         = min(512, max(32, int(real(N_clumps,wp)**(1.0_wp/3.0_wp)) + 1))
  ncells_rsa = rg**3
  ! RSA grid spans the bounding sphere; cells must be at least 2*cl_radius_max
  ! across so the 27-neighbour search captures every possible overlap.
  rg_cell    = max(2.0_wp * sphere_R / real(rg, wp), 2.0_wp * cl_radius_max)
  rg         = max(2, int((2.0_wp * sphere_R) / rg_cell) + 1)
  rg_cell    = (2.0_wp * sphere_R) / real(rg, wp)
  ncells_rsa = rg**3
  min_sep2_uniform = (2.0_wp * cl_radius_max)**2

  allocate(head(ncells_rsa),    source=-1)
  allocate(nxt(int(N_clumps)), source=-1)

  if (mpar%p_rank == 0) &
     write(*,'(a,i5,a,i14,a)') ' RSA: grid ', rg, '^3, placing ', N_clumps, ' clumps...'

  n_attempts = 0_int64
  icl        = 0_int64
  rcl_trial  = base_radius_in   ! default for uniform path

  do while (icl < N_clumps)
     n_attempts = n_attempts + 1_int64

     if (profiles_active) then
        !--- inverse-CDF radial draw + isotropic angles
        r_trial    = sample_clump_radius()
        rcl_trial  = base_radius_in * shape_radius(r_trial)
        cos_theta  = 2.0_wp * rand_number() - 1.0_wp
        sin_theta  = sqrt(max(0.0_wp, 1.0_wp - cos_theta*cos_theta))
        phi_az     = twopi * rand_number()
        xc = r_trial * sin_theta * cos(phi_az)
        yc = r_trial * sin_theta * sin(phi_az)
        zc = r_trial * cos_theta
     else
        !--- uniform random point inside sphere (legacy box rejection)
        do
           xc = (2.0_wp * rand_number() - 1.0_wp) * sphere_R
           yc = (2.0_wp * rand_number() - 1.0_wp) * sphere_R
           zc = (2.0_wp * rand_number() - 1.0_wp) * sphere_R
           if (xc*xc + yc*yc + zc*zc <= sphere_R*sphere_R) exit
        end do
     end if

     !--- cell in RSA grid (0-based)
     ig = min(rg-1, max(0, int((xc + sphere_R) / rg_cell)))
     jg = min(rg-1, max(0, int((yc + sphere_R) / rg_cell)))
     kg = min(rg-1, max(0, int((zc + sphere_R) / rg_cell)))

     !--- check 27 neighbours for overlap
     overlap = .false.
     outer: do kg2 = max(0,kg-1), min(rg-1,kg+1)
        do jg2 = max(0,jg-1), min(rg-1,jg+1)
           do ig2 = max(0,ig-1), min(rg-1,ig+1)
              icell_rsa = 1 + ig2 + jg2*rg + kg2*rg*rg
              jnb = head(icell_rsa)
              do while (jnb > 0)
                 dx = xc - cl_x(jnb);  dy = yc - cl_y(jnb);  dz = zc - cl_z(jnb)
                 d2 = dx*dx + dy*dy + dz*dz
                 if (profiles_active) then
                    sep_pair = rcl_trial + cl_radius(int(jnb,int64))
                    if (d2 < sep_pair*sep_pair) then
                       overlap = .true.;  exit outer
                    end if
                 else
                    if (d2 < min_sep2_uniform) then
                       overlap = .true.;  exit outer
                    end if
                 end if
                 jnb = nxt(jnb)
              end do
           end do
        end do
     end do outer
     if (overlap) cycle

     !--- accept
     icl = icl + 1_int64
     cl_x(icl) = real(xc, dp);  cl_y(icl) = real(yc, dp);  cl_z(icl) = real(zc, dp)

     if (profiles_active) then
        !--- per-clump physical assignments from profiles (Phase 5)
        cl_radius(icl)     = rcl_trial
        cl_radius2(icl)    = rcl_trial * rcl_trial
        if (allocated(tab_temperature)) then
           temp_loc = profile_file_interp(r_trial, 4)
           if (temp_loc <= 0.0_wp) temp_loc = cl_temperature_ref
        else
           temp_loc = cl_temperature_ref
        end if
        cl_temperature(icl) = temp_loc
        vth_loc = line%vtherm1 * sqrt(temp_loc)
        cl_vtherm(icl)      = vth_loc
        Df_loc  = vth_loc / (line%wavelength0 * um2km)
        cl_Dfreq(icl)       = Df_loc
        va_loc  = (line%damping / fourpi) / Df_loc
        cl_voigt_a(icl)     = va_loc
        ! Local opacity scaling: rhokap ∝ n_H / Dfreq.
        ! base_rhokap_in encodes the peak rhokap at Dfreq_ref.
        ! shape_density(r) gives the n_H ratio relative to peak.
        dens_factor   = shape_density(r_trial)
        kap_loc       = base_rhokap_in * dens_factor * (cl_Dfreq_ref / Df_loc)
        cl_rhokap(icl) = kap_loc
     end if

     if (par%clump_sigma_v > 0.0_wp) then
        ! store v / cl_vtherm(icl) directly (dimensionless, matches AMR convention)
        cl_vx(icl) = real(par%clump_sigma_v / cl_vtherm(icl) * rand_gauss(), dp)
        cl_vy(icl) = real(par%clump_sigma_v / cl_vtherm(icl) * rand_gauss(), dp)
        cl_vz(icl) = real(par%clump_sigma_v / cl_vtherm(icl) * rand_gauss(), dp)
     end if

     !--- insert into RSA linked list
     icell_rsa       = 1 + ig + jg*rg + kg*rg*rg
     nxt(int(icl))   = head(icell_rsa)
     head(icell_rsa) = int(icl)

     if (mod(icl, 1000000_int64) == 0_int64 .and. mpar%p_rank == 0) &
        write(*,'(a,i14,a,i14)') '   placed ', icl, ' / ', N_clumps
  end do

  if (mpar%p_rank == 0) write(*,'(a,f6.1,a)') &
     ' RSA done, acceptance rate = ', &
     real(N_clumps,wp)/real(n_attempts,wp)*100.0_wp, '%'

  deallocate(head, nxt)
  end subroutine generate_clumps
  !===========================================================================

  !===========================================================================
  subroutine assign_clump_velocities_from_type()
  !---------------------------------------------------------------------------
  ! Add a systematic velocity component to each clump based on par%velocity_type
  ! and the clump's centre position.  Called only on h_rank=0 after
  ! generate_clumps(); adds to any existing sigma_v random component.
  !
  ! Supported types (same parameter names as grid_mod_car):
  !   'hubble'              : v = Vexp * (r / sphere_R)       [linear expansion]
  !   'constant_radial'     : v = Vexp * (r / |r|)           [uniform outflow]
  !   'parallel_velocity'   : v = (Vx, Vy, Vz)               [uniform bulk]
  !   'ssh'                 : Song, Seon & Hwang (2020) galaxy model
  !   'rotating_solid_body' : v = Vrot * (-y, x, 0) / sphere_R
  !   'rotating_galaxy_halo': flat rotation curve (Vrot, rinner)
  !---------------------------------------------------------------------------
  implicit none
  integer(int64) :: icl
  real(kind=wp)  :: xc, yc, zc, rr, rr_cyl, Vscale, vx, vy, vz

  do icl = 1_int64, N_clumps
     xc = real(cl_x(icl), wp)
     yc = real(cl_y(icl), wp)
     zc = real(cl_z(icl), wp)
     rr = sqrt(xc**2 + yc**2 + zc**2)
     vx = 0.0_wp;  vy = 0.0_wp;  vz = 0.0_wp

     select case (trim(par%velocity_type))

     case ('hubble')
        ! v_i = Vexp * r_i / sphere_R  [km/s]
        vx = par%Vexp * xc / sphere_R
        vy = par%Vexp * yc / sphere_R
        vz = par%Vexp * zc / sphere_R

     case ('constant_radial')
        ! v = Vexp * r_hat  [km/s]
        if (rr > 0.0_wp) then
           vx = par%Vexp * xc / rr
           vy = par%Vexp * yc / rr
           vz = par%Vexp * zc / rr
        end if

     case ('parallel_velocity')
        ! uniform bulk velocity  [km/s]
        vx = par%Vx
        vy = par%Vy
        vz = par%Vz

     case ('ssh')
        ! Song, Seon & Hwang (2020): linear inside rpeak, then Vpeak + DeltaV*(r-rpeak)/(R-rpeak)
        if (rr > 0.0_wp) then
           if (rr < par%rpeak) then
              Vscale = par%Vpeak / par%rpeak
              vx = Vscale * xc
              vy = Vscale * yc
              vz = Vscale * zc
           else
              Vscale = par%Vpeak + par%DeltaV * (rr - par%rpeak) / (sphere_R - par%rpeak)
              vx = Vscale * xc / rr
              vy = Vscale * yc / rr
              vz = Vscale * zc / rr
           end if
        end if

     case ('rotating_solid_body')
        ! solid-body rotation about z-axis: v = Vrot * (-y, x, 0) / sphere_R
        vx = -par%Vrot * yc / sphere_R
        vy =  par%Vrot * xc / sphere_R

     case ('rotating_galaxy_halo')
        ! flat rotation curve: Vrot inside rinner, then Vrot * rinner / rr
        rr_cyl = sqrt(xc**2 + yc**2)
        if (rr_cyl > 0.0_wp) then
           if (rr_cyl < par%rinner) then
              vx = -par%Vrot * yc / par%rinner
              vy =  par%Vrot * xc / par%rinner
           else
              vx = -par%Vrot * yc / rr_cyl
              vy =  par%Vrot * xc / rr_cyl
           end if
        end if

     end select

     ! store v / cl_vtherm(icl) directly (dimensionless, matches AMR convention)
     cl_vx(icl) = cl_vx(icl) + real(vx / cl_vtherm(icl), dp)
     cl_vy(icl) = cl_vy(icl) + real(vy / cl_vtherm(icl), dp)
     cl_vz(icl) = cl_vz(icl) + real(vz / cl_vtherm(icl), dp)
  end do

  if (mpar%p_rank == 0) &
     write(*,'(2a)') ' Clumps: velocity_type   = ', trim(par%velocity_type)

  end subroutine assign_clump_velocities_from_type
  !===========================================================================

  !===========================================================================
  subroutine build_clump_csr()
  !---------------------------------------------------------------------------
  ! Two-pass CSR construction. All ranks allocate shared memory;
  ! h_rank=0 fills it. Caller handles barrier after return.
  !---------------------------------------------------------------------------
  implicit none
  integer        :: ncells, total_regs
  integer        :: i, j, k, imin, imax, jmin, jmax, kmin, kmax, icell
  integer(int64) :: icl
  integer,       allocatable :: cnt(:)
  integer        :: ierr

  !--- CSR grid size: ~1 clump per cell
  cgx = min(512, max(32, int(real(N_clumps,wp)**(1.0_wp/3.0_wp)) + 1))
  cgy = cgx;  cgz = cgx
  ncells = cgx * cgy * cgz

  !--- bounding box: sphere + one clump-radius margin (use cl_radius_max so
  !    the largest clump still fits cleanly into the CSR cell scheme)
  cg_xmin = -(sphere_R + cl_radius_max);  cg_ymin = cg_xmin;  cg_zmin = cg_xmin
  cg_dx   = (2.0_wp*(sphere_R + cl_radius_max)) / real(cgx, wp)
  cg_dy   = cg_dx;  cg_dz = cg_dx
  cg_inv_dx = 1.0_wp / cg_dx
  cg_inv_dy = 1.0_wp / cg_dy
  cg_inv_dz = 1.0_wp / cg_dz

  !--- allocate cg_start on all ranks (h_rank=0 fills)
  call create_shared_mem(cg_start, [ncells + 1])

  !--- pass 1: count (h_rank=0 only)
  if (mpar%h_rank == 0) then
     allocate(cnt(ncells), source=0)
     do icl = 1_int64, N_clumps
        call clump_cell_range(icl, imin, imax, jmin, jmax, kmin, kmax)
        do k = kmin, kmax
           do j = jmin, jmax
              do i = imin, imax
                 icell = cg_cell_idx(i, j, k)
                 cnt(icell) = cnt(icell) + 1
              end do
           end do
        end do
     end do
     !--- prefix sum → cg_start (1-based: cg_start(1)=1, cg_start(ncells+1)=total+1)
     cg_start(1) = 1
     do icell = 1, ncells
        cg_start(icell+1) = cg_start(icell) + cnt(icell)
     end do
  end if
  call MPI_BARRIER(mpar%hostcomm, ierr)

  total_regs = cg_start(ncells+1) - 1
  call create_shared_mem(cg_list, [total_regs])

  !--- pass 2: fill (h_rank=0 only)
  if (mpar%h_rank == 0) then
     cnt(:) = 0
     do icl = 1_int64, N_clumps
        call clump_cell_range(icl, imin, imax, jmin, jmax, kmin, kmax)
        do k = kmin, kmax
           do j = jmin, jmax
              do i = imin, imax
                 icell = cg_cell_idx(i, j, k)
                 cg_list(cg_start(icell) + cnt(icell)) = int(icl, int32)
                 cnt(icell) = cnt(icell) + 1
              end do
           end do
        end do
     end do
     deallocate(cnt)
     if (mpar%p_rank == 0) write(*,'(a,i12,a,i5,a)') &
           ' CSR grid: ', total_regs, ' registrations in ', cgx, '^3 cells.'
  end if

  end subroutine build_clump_csr
  !===========================================================================

  !===========================================================================
  pure subroutine clump_cell_range(icl, imin, imax, jmin, jmax, kmin, kmax)
  !---------------------------------------------------------------------------
  ! CSR cells whose box overlaps clump icl's sphere.
  !---------------------------------------------------------------------------
  integer(int64), intent(in)  :: icl
  integer,        intent(out) :: imin, imax, jmin, jmax, kmin, kmax
  real(kind=wp) :: rcl
  rcl = cl_radius(icl)
  imin = max(0, int((cl_x(icl) - cg_xmin - rcl) * cg_inv_dx))
  imax = min(cgx-1, int((cl_x(icl) - cg_xmin + rcl) * cg_inv_dx))
  jmin = max(0, int((cl_y(icl) - cg_ymin - rcl) * cg_inv_dy))
  jmax = min(cgy-1, int((cl_y(icl) - cg_ymin + rcl) * cg_inv_dy))
  kmin = max(0, int((cl_z(icl) - cg_zmin - rcl) * cg_inv_dz))
  kmax = min(cgz-1, int((cl_z(icl) - cg_zmin + rcl) * cg_inv_dz))
  end subroutine clump_cell_range
  !===========================================================================

  !===========================================================================
  pure subroutine ray_sphere_isect(ox, oy, oz, kx, ky, kz, cx, cy, cz, &
                                    icl, t_entry, t_exit, hit)
  !---------------------------------------------------------------------------
  ! Ray–sphere intersection.  Ray: P(t) = origin + t*dir, t arbitrary.
  ! Sphere: centre (cx,cy,cz), squared-radius cl_radius2(icl).
  ! Returns hit=.true. and t_entry <= t_exit when t_exit > 0.
  !---------------------------------------------------------------------------
  real(kind=wp),  intent(in)  :: ox, oy, oz, kx, ky, kz, cx, cy, cz
  integer(int64), intent(in)  :: icl
  real(kind=wp),  intent(out) :: t_entry, t_exit
  logical,        intent(out) :: hit
  real(kind=wp) :: rx, ry, rz, b, disc
  rx = ox - cx;  ry = oy - cy;  rz = oz - cz
  b    = rx*kx + ry*ky + rz*kz
  disc = b*b - (rx*rx + ry*ry + rz*rz) + cl_radius2(icl)
  if (disc < 0.0_wp) then
     hit = .false.;  t_entry = 0.0_wp;  t_exit = 0.0_wp
  else
     disc    = sqrt(disc)
     t_entry = -b - disc
     t_exit  = -b + disc
     hit     = (t_exit > 0.0_wp)
  end if
  end subroutine ray_sphere_isect
  !===========================================================================

  !===========================================================================
  subroutine find_next_clump(xp, yp, zp, kx, ky, kz, skip_icl, t_max, &
                              t_entry, t_exit, icl_found, found)
  !---------------------------------------------------------------------------
  ! DDA through CSR grid to find the nearest clump hit along the ray
  !   P(t) = (xp,yp,zp) + t*(kx,ky,kz),  0 < t <= t_max.
  ! Clump skip_icl (> 0) is excluded (the one just exited).
  !
  ! Pattern: Amanatides & Woo (1987).
  !   tx,ty,tz = absolute path length to NEXT x/y/z face crossing.
  !   d        = path length at current cell entry.
  ! Stopping: once d > best_te, no earlier hit can be found.
  !---------------------------------------------------------------------------
  real(kind=wp),  intent(in)  :: xp, yp, zp, kx, ky, kz
  integer(int64), intent(in)  :: skip_icl
  real(kind=wp),  intent(in)  :: t_max
  real(kind=wp),  intent(out) :: t_entry, t_exit
  integer(int64), intent(out) :: icl_found
  logical,        intent(out) :: found

  integer  :: ci, cj, ck, si, sj, sk, icell, ip
  integer(int64) :: icl
  real(kind=wp) :: tx, ty, tz, delx, dely, delz, d
  real(kind=wp) :: te, tx2, best_te, best_tx2
  integer(int64):: best_icl
  logical  :: hit

  found    = .false.
  best_te  = hugest
  best_tx2 = 0.0_wp
  best_icl = 0_int64
  d        = 0.0_wp

  !--- starting cell (clamped to grid)
  ci = max(0, min(cgx-1, int((xp - cg_xmin) * cg_inv_dx)))
  cj = max(0, min(cgy-1, int((yp - cg_ymin) * cg_inv_dy)))
  ck = max(0, min(cgz-1, int((zp - cg_zmin) * cg_inv_dz)))

  !--- DDA setup: absolute path lengths to first face crossing
  if (kx > 0.0_wp) then
     si   =  1
     tx   = ((cg_xmin + real(ci+1,wp)*cg_dx) - xp) / kx
     delx =  cg_dx / kx
  else if (kx < 0.0_wp) then
     si   = -1
     tx   = ((cg_xmin + real(ci,wp)*cg_dx) - xp) / kx
     delx = -cg_dx / kx
  else
     si = 0;  tx = hugest;  delx = hugest
  end if

  if (ky > 0.0_wp) then
     sj   =  1
     ty   = ((cg_ymin + real(cj+1,wp)*cg_dy) - yp) / ky
     dely =  cg_dy / ky
  else if (ky < 0.0_wp) then
     sj   = -1
     ty   = ((cg_ymin + real(cj,wp)*cg_dy) - yp) / ky
     dely = -cg_dy / ky
  else
     sj = 0;  ty = hugest;  dely = hugest
  end if

  if (kz > 0.0_wp) then
     sk   =  1
     tz   = ((cg_zmin + real(ck+1,wp)*cg_dz) - zp) / kz
     delz =  cg_dz / kz
  else if (kz < 0.0_wp) then
     sk   = -1
     tz   = ((cg_zmin + real(ck,wp)*cg_dz) - zp) / kz
     delz = -cg_dz / kz
  else
     sk = 0;  tz = hugest;  delz = hugest
  end if

  do while(.true.)
     !--- stopping: current cell starts beyond the best hit or past t_max
     if (d > best_te .or. d > t_max) exit

     !--- check all clumps in current cell
     icell = cg_cell_idx(ci, cj, ck)
     do ip = cg_start(icell), cg_start(icell+1) - 1
        icl = int(cg_list(ip), int64)
        if (icl == skip_icl) cycle
        call ray_sphere_isect(xp, yp, zp, kx, ky, kz, &
             real(cl_x(icl),wp), real(cl_y(icl),wp), real(cl_z(icl),wp), &
             icl, te, tx2, hit)
        if (hit .and. tx2 > 0.0_wp .and. te < best_te .and. &
            (te > 0.0_wp .or. icl /= skip_icl)) then
           best_te  = te
           best_tx2 = tx2
           best_icl = icl
        end if
     end do

     !--- advance to next cell (Amanatides & Woo pattern)
     if (tx <= ty .and. tx <= tz) then
        d  = tx
        ci = ci + si;  if (ci < 0 .or. ci >= cgx) exit
        tx = tx + delx
     else if (ty <= tz) then
        d  = ty
        cj = cj + sj;  if (cj < 0 .or. cj >= cgy) exit
        ty = ty + dely
     else
        d  = tz
        ck = ck + sk;  if (ck < 0 .or. ck >= cgz) exit
        tz = tz + delz
     end if
  end do

  if (best_icl > 0_int64 .and. best_te <= t_max) then
     found     = .true.
     t_entry   = best_te
     t_exit    = min(best_tx2, t_max)
     icl_found = best_icl
  end if

  end subroutine find_next_clump
  !===========================================================================

  !===========================================================================
  pure real(kind=wp) function clump_exit_dist(xp, yp, zp, kx, ky, kz, icl)
  !---------------------------------------------------------------------------
  ! Distance from (xp,yp,zp) moving in direction (kx,ky,kz) to the exit point
  ! of clump icl. Assumes the photon is currently inside the clump.
  !---------------------------------------------------------------------------
  real(kind=wp),  intent(in) :: xp, yp, zp, kx, ky, kz
  integer(int64), intent(in) :: icl
  real(kind=wp) :: rx, ry, rz, b, disc
  rx = xp - real(cl_x(icl),wp)
  ry = yp - real(cl_y(icl),wp)
  rz = zp - real(cl_z(icl),wp)
  b    = rx*kx + ry*ky + rz*kz
  disc = b*b - (rx*rx + ry*ry + rz*rz) + cl_radius2(icl)
  if (disc < 0.0_wp) disc = 0.0_wp
  clump_exit_dist = max(0.0_wp, -b + sqrt(disc))
  end function clump_exit_dist
  !===========================================================================

  !===========================================================================
  pure real(kind=wp) function sphere_exit_dist(xp, yp, zp, kx, ky, kz)
  !---------------------------------------------------------------------------
  ! Distance to the exit of the outer sphere (radius sphere_R) from (xp,yp,zp).
  !---------------------------------------------------------------------------
  real(kind=wp), intent(in) :: xp, yp, zp, kx, ky, kz
  real(kind=wp) :: b, disc
  b    = xp*kx + yp*ky + zp*kz
  disc = b*b - (xp*xp + yp*yp + zp*zp) + sphere_R*sphere_R
  if (disc < 0.0_wp) disc = 0.0_wp
  sphere_exit_dist = max(0.0_wp, -b + sqrt(disc))
  end function sphere_exit_dist
  !===========================================================================

  !===========================================================================
  subroutine destroy_clumps()
  implicit none
  ! destroy_mem fails on non-root ranks because get_window compares against
  ! rank-0 base addresses (MPI_WIN_SHARED_QUERY rank=0), which differ from
  ! local pointers on other ranks.  Use destroy_shared_mem_all() so that
  ! MPI_WIN_FREE is called collectively, then nullify all pointers manually.
  call destroy_shared_mem_all()
  nullify(cl_x, cl_y, cl_z, cl_vx, cl_vy, cl_vz)
  nullify(cl_radius, cl_radius2, cl_rhokap, cl_voigt_a, cl_Dfreq, &
          cl_vtherm, cl_temperature)
  nullify(cg_start, cg_list)
  N_clumps = 0_int64
  end subroutine destroy_clumps
  !===========================================================================

  !===========================================================================
  subroutine write_clumps_fits(fname)
  !---------------------------------------------------------------------------
  ! Write clump positions, velocities, and physical parameters to a FITS
  ! binary table.  Called only on p_rank == 0 after init_clumps().
  !
  ! HDU 1 (primary): empty image, all clump scalars as header keywords.
  ! HDU 2 (BinTable): columns X, Y, Z [code units], VX, VY, VZ [km/s].
  !---------------------------------------------------------------------------
  use define
  use fitsio_mod
  implicit none
  character(len=*), intent(in) :: fname

  integer :: unit, status, bitpix
  real(kind=real64), allocatable :: tmp(:)
  integer(int64) :: ncl, i
  real(kind=wp)  :: f_vol_actual, f_cov_actual, total_HI_mass
  real(kind=wp)  :: rcl_min, rcl_max, rcl_mean, tau_mean
  real(kind=wp)  :: kap_min_w, kap_max_w, kap_mean_w
  real(kind=wp)  :: T_min_w, T_max_w, T_mean_w
  logical        :: write_radius, write_rhokap, write_temp
  real(kind=wp), parameter :: const_tol = 1.0e-3_wp  ! relative spread threshold

  status = 0
  ncl    = N_clumps
  bitpix = -32   ! save as float32 (sufficient precision for positions/velocities)

  !--- Compute realized f_vol and f_cov from actual placed clumps.
  !    Uniform case: closed-form (preserves Phase 1 numbers).
  !    Profile case: 4π/3 * Σ r_cl^3 / V_box and (LOS f_cov from quadrature).
  !    From-file case: per-clump sums (no profile tables available).
  if (clumps_from_file) then
     f_vol_actual = 0.0_wp
     f_cov_actual = 0.0_wp
     do i = 1_int64, ncl
        f_vol_actual = f_vol_actual + cl_radius(i)**3
        f_cov_actual = f_cov_actual + cl_radius(i)**2
     end do
     f_vol_actual = f_vol_actual / sphere_R**3
     f_cov_actual = 0.75_wp * f_cov_actual / sphere_R**2
  else if (.not. profiles_active) then
     f_vol_actual = real(ncl,wp) * (cl_radius_max / sphere_R)**3
     f_cov_actual = 0.75_wp * real(ncl,wp) * (cl_radius_max / sphere_R)**2
  else
     !--- realised f_vol from sum of individual clump volumes
     f_vol_actual = 0.0_wp
     do i = 1_int64, ncl
        f_vol_actual = f_vol_actual + cl_radius(i)**3
     end do
     f_vol_actual = f_vol_actual / sphere_R**3
     !--- realised f_cov from radial quadrature (LOS integral)
     f_cov_actual = f_cov_LOS_quad(real(ncl,wp) / max(total_count_quad(1.0_wp), tiny(1.0_wp)))
  end if

  !--- per-clump diagnostics (peak, mean, range)
  rcl_min  = huge(1.0_wp);  rcl_max = 0.0_wp;  rcl_mean = 0.0_wp;  tau_mean = 0.0_wp
  do i = 1_int64, ncl
     rcl_min  = min(rcl_min,  cl_radius(i))
     rcl_max  = max(rcl_max,  cl_radius(i))
     rcl_mean = rcl_mean + cl_radius(i)
     tau_mean = tau_mean + cl_rhokap(i) * voigt(0.0_wp, cl_voigt_a(i)) * cl_radius(i)
  end do
  if (ncl > 0_int64) then
     rcl_mean = rcl_mean / real(ncl, wp)
     tau_mean = tau_mean / real(ncl, wp)
  end if

  !--- estimate of total HI mass [in proton masses, when distance2cm > 0]:
  !    M_HI ∝ Σ n_H_clump * V_clump = Σ (rhokap_i * Dfreq_i / cross0) * V_i (in code units),
  !    converted to proton masses by * distance2cm**3 * m_proton.
  total_HI_mass = 0.0_wp
  do i = 1_int64, ncl
     total_HI_mass = total_HI_mass + (cl_rhokap(i) * cl_Dfreq(i) / line%cross0) &
                                  * (4.0_wp/3.0_wp) * pi * cl_radius(i)**3
  end do

  if (mpar%p_rank == 0) then
     write(*,'(a,es12.4)') ' Clumps: realised f_vol      = ', f_vol_actual
     write(*,'(a,es12.4)') ' Clumps: realised f_cov      = ', f_cov_actual
     write(*,'(a,3es12.4)')' Clumps: r_cl min/mean/max   = ', rcl_min, rcl_mean, rcl_max
     write(*,'(a,es12.4)') ' Clumps: mean tau_per_clump  = ', tau_mean
     write(*,'(a,es12.4)') ' Clumps: integrated HI col   = ', total_HI_mass
  end if

  call fits_open_new(unit, trim(fname), status)
  if (status /= 0) then
     write(*,*) 'WARNING: write_clumps_fits: cannot open ', trim(fname)
     return
  end if

  !--- BinTable HDU: X, Y, Z, VX, VY, VZ + (conditionally) R_CLUMP, RHOKAP, TEMP
  !
  !    The three per-clump physical columns are only written when their
  !    spread (max-min) exceeds const_tol * |mean|. Otherwise the value
  !    is fully captured by the corresponding header keyword (CL_RAD,
  !    RHOKAP, TEMP_CL) and the column would be a 4-byte-per-row constant
  !    waste. read_clumps_fits() falls back to the header keyword when the
  !    column is absent, so old 6-column files (predating Phase 6) and new
  !    constant-case files are read the same way.
  allocate(tmp(ncl))

  tmp = cl_x(1:ncl)
  call fits_write_table_column(unit, 'X',  tmp, status, bitpix)
  tmp = cl_y(1:ncl)
  call fits_write_table_column(unit, 'Y',  tmp, status, bitpix)
  tmp = cl_z(1:ncl)
  call fits_write_table_column(unit, 'Z',  tmp, status, bitpix)
  !--- cl_v* are stored as v/cl_vtherm(icl) internally (since 2026-04-27);
  !    multiply back by cl_vtherm(icl) to write velocities in km/s.
  tmp = cl_vx(1:ncl) * cl_vtherm(1:ncl)
  call fits_write_table_column(unit, 'VX', tmp, status, bitpix)
  tmp = cl_vy(1:ncl) * cl_vtherm(1:ncl)
  call fits_write_table_column(unit, 'VY', tmp, status, bitpix)
  tmp = cl_vz(1:ncl) * cl_vtherm(1:ncl)
  call fits_write_table_column(unit, 'VZ', tmp, status, bitpix)

  !--- per-clump physical columns (Phase 6 additions): write only if non-constant.
  kap_min_w  = minval(cl_rhokap(1:ncl))
  kap_max_w  = maxval(cl_rhokap(1:ncl))
  kap_mean_w = sum(cl_rhokap(1:ncl)) / max(real(ncl,wp), 1.0_wp)
  T_min_w    = minval(cl_temperature(1:ncl))
  T_max_w    = maxval(cl_temperature(1:ncl))
  T_mean_w   = sum(cl_temperature(1:ncl)) / max(real(ncl,wp), 1.0_wp)

  write_radius = (rcl_max  - rcl_min ) > const_tol * max(abs(rcl_mean ), tiny(1.0_wp))
  write_rhokap = (kap_max_w - kap_min_w) > const_tol * max(abs(kap_mean_w), tiny(1.0_wp))
  write_temp   = (T_max_w   - T_min_w  ) > const_tol * max(abs(T_mean_w  ), tiny(1.0_wp))

  if (write_radius) then
     tmp = cl_radius(1:ncl)
     call fits_write_table_column(unit, 'R_CLUMP', tmp, status, bitpix)
  end if
  if (write_rhokap) then
     tmp = cl_rhokap(1:ncl)
     call fits_write_table_column(unit, 'RHOKAP', tmp, status, bitpix)
  end if
  if (write_temp) then
     tmp = cl_temperature(1:ncl)
     call fits_write_table_column(unit, 'TEMP', tmp, status, bitpix)
  end if

  deallocate(tmp)

  if (mpar%p_rank == 0) then
     write(*,'(a,3l2)') ' Clumps: per-clump columns saved (R_CLUMP/RHOKAP/TEMP) = ', &
          write_radius, write_rhokap, write_temp
  end if

  !--- Write keywords into the BinTable HDU (current HDU after fits_write_table_column)
  call fits_put_keyword(unit, 'N_CLUMPS',  ncl,                'number of clumps (realized)',          status)
  call fits_put_keyword(unit, 'SPHERE_R',  sphere_R,           'outer sphere radius [code units]',      status)
  call fits_put_keyword(unit, 'CL_RAD',    cl_radius_max,      'clump radius (max) [code units]',       status)
  call fits_put_keyword(unit, 'F_VOL',     f_vol_actual,       'volume filling factor (realized)',      status)
  call fits_put_keyword(unit, 'F_COV',     f_cov_actual,       'covering factor (realized)',            status)
  call fits_put_keyword(unit, 'TAU0',      par%clump_tau0,     'line-center tau (center to surface)',   status)
  call fits_put_keyword(unit, 'SIGMA_V',   par%clump_sigma_v,  'bulk velocity sigma [km/s]',            status)
  call fits_put_keyword(unit, 'TEMP_CL',   cl_temperature_ref, 'clump temperature [K] (reference)',     status)
  call fits_put_keyword(unit, 'RHOKAP',    cl_rhokap_ref,      'opacity/code-unit inside clumps (ref)', status)
  call fits_put_keyword(unit, 'CL_DFREQ',  cl_Dfreq_ref,       'Doppler frequency [Hz] (reference)',    status)
  call fits_put_keyword(unit, 'VTHERM',    cl_vtherm_ref,      'thermal velocity [km/s] (reference)',   status)
  call fits_put_keyword(unit, 'VOIGT_A',   cl_voigt_a_ref,     'Voigt damping parameter (reference)',   status)
  call fits_put_keyword(unit, 'RMAX',      par%rmax,      'outer sphere radius input (par%rmax)',   status)
  call fits_put_keyword(unit, 'IN_FCOV',   par%clump_f_cov,   'covering factor (input)',            status)
  call fits_put_keyword(unit, 'IN_FVOL',   par%clump_f_vol,   'volume filling factor (input)',      status)
  call fits_put_keyword(unit, 'IN_NCL',    par%clump_N_clumps,'N_clumps (input)',                   status)
  call fits_put_keyword(unit, 'IN_NHI',    par%clump_NHI,     'per-clump NHI input [cm^-2] (clump center to surface)',status)
  call fits_put_keyword(unit, 'IN_NH',     par%clump_nH,      'clump nH density input [cm^-3]',    status)
  call fits_put_keyword(unit, 'IN_TEMP',   par%clump_temperature,'clump temperature input [K]',     status)

  call fits_close(unit, status)

  if (mpar%p_rank == 0) write(*,'(2a)') ' Clumps saved to ', trim(fname)
  end subroutine write_clumps_fits
  !===========================================================================

  !===========================================================================
  subroutine read_perclump_or_keyword(unit, colnames, keyname, scratch, ncl, dst)
  !---------------------------------------------------------------------------
  ! Helper used by read_clumps_fits.  Try each comma-separated entry in
  ! `colnames` in order; the first one that exists in the FITS table is
  ! used.  Allows alternative spellings to be accepted (e.g. RHOKAP,
  ! DENSITY, DENS for the per-clump opacity proxy).  If none of the
  ! columns are present, fall back to the value of header keyword
  ! `keyname` and fill `dst` uniformly with it.  This is how
  !   - a pre-Phase-6 file without the R_CLUMP/RHOKAP/TEMP columns, and
  !   - a post-Phase-8 column-skipping file where the column was dropped
  !     because it was effectively constant
  ! are both handled with a single read path.
  !---------------------------------------------------------------------------
  use fitsio_mod
  implicit none
  integer,          intent(in)    :: unit
  character(len=*), intent(in)    :: colnames, keyname
  real(real32),     intent(inout) :: scratch(:)
  integer(int64),   intent(in)    :: ncl
  real(kind=wp),    intent(out)   :: dst(:)

  integer            :: status, colnum, ierr, ks, ke, lstr
  real(kind=wp)      :: keyval
  character(len=80)  :: cmt
  character(len=64)  :: candidate
  character(len=64)  :: tried
  logical            :: found

  found = .false.
  tried = ''
  lstr  = len_trim(colnames)
  ks    = 1
  do while (ks <= lstr .and. .not. found)
     ke = index(colnames(ks:lstr), ',')
     if (ke == 0) then
        candidate = adjustl(colnames(ks:lstr))
        ks = lstr + 1
     else
        candidate = adjustl(colnames(ks:ks+ke-2))
        ks = ks + ke
     end if
     if (len_trim(candidate) == 0) cycle
     tried = trim(tried)//' '//trim(candidate)
     status = 0
     call fits_get_column_number(unit, trim(candidate), colnum, status)
     if (status == 0) then
        call fits_read_table_column(unit, colnum, scratch, status)
        dst(1:ncl) = real(scratch(1:ncl), wp)
        found = .true.
     end if
  end do

  if (.not. found) then
     status = 0
     call fits_get_keyword(unit, trim(keyname), keyval, status, cmt)
     if (status /= 0) then
        write(*,'(5a)') 'ERROR: clump_input_file is missing every column in {', &
             trim(adjustl(tried)), ' } and header keyword ', trim(keyname), &
             '; cannot fall back to a constant value'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
     dst(1:ncl) = keyval
  end if
  end subroutine read_perclump_or_keyword
  !===========================================================================

  !===========================================================================
  subroutine read_clumps_fits(fname, R_sphere)
  !---------------------------------------------------------------------------
  ! Read a clump population from a FITS file produced by write_clumps_fits()
  ! (standalone make_clumps.x driver, or a previous LaRT run with
  ! par%save_clump_info = .true.). The HDU layout follows that subroutine:
  !   HDU 2 (BinTable) columns:
  !     X, Y, Z          [code units]
  !     VX, VY, VZ       [km/s]
  !     R_CLUMP, RHOKAP, TEMP    (per-clump physical props)
  !
  ! Allocates the same MPI shared-memory arrays as the generate path,
  ! distributes via the standard hostcomm/SAME_HRANK_COMM pattern, and
  ! builds the CSR acceleration grid afterwards.  Sets clumps_from_file=.true.
  ! so write_clumps_fits and other diagnostic branches can detect this case.
  !---------------------------------------------------------------------------
  use define
  use line_mod
  use fitsio_mod
  use voigt_mod, only: voigt
  use mpi
  implicit none
  character(len=*), intent(in) :: fname
  real(kind=wp),    intent(in) :: R_sphere

  integer        :: unit, status, ierr, colnum
  integer(int64) :: ncl, i
  real(kind=real32), allocatable :: tmp(:)
  character(len=80) :: cmt

  status   = 0
  sphere_R = R_sphere
  base_radius_in  = par%clump_radius
  profiles_active = .false.
  clumps_from_file = .true.

  if (mpar%p_rank == 0) then
     write(*,'(2a)') ' Clumps: reading from ', trim(fname)
     call fits_open_old(unit, trim(fname), status)
     if (status /= 0) then
        write(*,'(2a)') 'ERROR: cannot open clump_input_file: ', trim(fname)
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
     !--- HDU 1 is the empty primary; the binary table is HDU 2.
     call fits_move_to_next_hdu(unit, status)
     if (status /= 0) then
        write(*,'(a)') 'ERROR: clump_input_file: cannot reach binary-table HDU 2'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
     call fits_get_keyword(unit, 'N_CLUMPS', ncl, status, cmt)
     if (status /= 0 .or. ncl <= 0_int64) then
        write(*,'(a)') 'ERROR: clump_input_file: N_CLUMPS keyword missing or <= 0'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
  end if
  call MPI_BCAST(ncl, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
  N_clumps = ncl

  !--- Allocate the same MPI shared-memory arrays the generate path uses.
  call create_shared_mem(cl_x,  [int(N_clumps)])
  call create_shared_mem(cl_y,  [int(N_clumps)])
  call create_shared_mem(cl_z,  [int(N_clumps)])
  call create_shared_mem(cl_vx, [int(N_clumps)])
  call create_shared_mem(cl_vy, [int(N_clumps)])
  call create_shared_mem(cl_vz, [int(N_clumps)])

  call create_shared_mem(cl_radius,      [int(N_clumps)])
  call create_shared_mem(cl_radius2,     [int(N_clumps)])
  call create_shared_mem(cl_rhokap,      [int(N_clumps)])
  call create_shared_mem(cl_voigt_a,     [int(N_clumps)])
  call create_shared_mem(cl_Dfreq,       [int(N_clumps)])
  call create_shared_mem(cl_vtherm,      [int(N_clumps)])
  call create_shared_mem(cl_temperature, [int(N_clumps)])

  !--- Read columns on p_rank=0 (which is also h_rank=0 of node 0); the
  !    other h_rank=0 ranks receive copies via MPI_BCAST below.
  if (mpar%p_rank == 0) then
     allocate(tmp(ncl))

     call fits_get_column_number(unit, 'X', colnum, status)
     call fits_read_table_column(unit, colnum, tmp, status)
     cl_x = real(tmp, dp)
     call fits_get_column_number(unit, 'Y', colnum, status)
     call fits_read_table_column(unit, colnum, tmp, status)
     cl_y = real(tmp, dp)
     call fits_get_column_number(unit, 'Z', colnum, status)
     call fits_read_table_column(unit, colnum, tmp, status)
     cl_z = real(tmp, dp)
     call fits_get_column_number(unit, 'VX', colnum, status)
     call fits_read_table_column(unit, colnum, tmp, status)
     cl_vx = real(tmp, dp)        ! km/s; converted to dimensionless below
     call fits_get_column_number(unit, 'VY', colnum, status)
     call fits_read_table_column(unit, colnum, tmp, status)
     cl_vy = real(tmp, dp)
     call fits_get_column_number(unit, 'VZ', colnum, status)
     call fits_read_table_column(unit, colnum, tmp, status)
     cl_vz = real(tmp, dp)

     !--- R_CLUMP / RHOKAP / TEMP are optional. write_clumps_fits omits
     !    each one when its spread is below const_tol (= 1e-3 of mean).
     !    Old (pre-Phase-6) FITS files do not carry these columns at all.
     !    In both cases we fall back to the corresponding header keyword.
     call read_perclump_or_keyword(unit, 'R_CLUMP',                'CL_RAD',  tmp, ncl, cl_radius)
     call read_perclump_or_keyword(unit, 'RHOKAP,DENSITY,DENS',    'RHOKAP',  tmp, ncl, cl_rhokap)
     call read_perclump_or_keyword(unit, 'TEMP',                   'TEMP_CL', tmp, ncl, cl_temperature)
     deallocate(tmp)
     call fits_close(unit, status)

     !--- Derive per-clump Doppler frequency, thermal velocity, Voigt damping
     !    from the per-clump temperature so that the line-data choice (par%line_id)
     !    used in the LaRT run can differ from the one that generated the file.
     do i = 1_int64, ncl
        cl_vtherm(i)  = line%vtherm1 * sqrt(cl_temperature(i))
        cl_Dfreq(i)   = cl_vtherm(i) / (line%wavelength0 * um2km)
        cl_voigt_a(i) = (line%damping / fourpi) / cl_Dfreq(i)
        cl_radius2(i) = cl_radius(i) * cl_radius(i)
        !--- VX/VY/VZ in the FITS file are stored in km/s (write_clumps_fits
        !    multiplies by cl_vtherm). LaRT internally uses v / cl_vtherm(icl).
        cl_vx(i) = cl_vx(i) / cl_vtherm(i)
        cl_vy(i) = cl_vy(i) / cl_vtherm(i)
        cl_vz(i) = cl_vz(i) / cl_vtherm(i)
     end do
  end if

  !--- Distribute the populated arrays to every node-local h_rank=0; other
  !    h_rank values share via MPI-3 shared memory after the host barrier.
  call MPI_BARRIER(mpar%hostcomm, ierr)
  if (mpar%h_rank == 0) then
     call MPI_BCAST(cl_x,           int(ncl), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_y,           int(ncl), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_z,           int(ncl), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_vx,          int(ncl), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_vy,          int(ncl), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_vz,          int(ncl), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_radius,      int(ncl), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_radius2,     int(ncl), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_rhokap,      int(ncl), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_voigt_a,     int(ncl), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_Dfreq,       int(ncl), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_vtherm,      int(ncl), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_temperature, int(ncl), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
  end if
  call MPI_BARRIER(mpar%hostcomm, ierr)

  !--- Reference scalars from the loaded data. cl_radius_max sets the RSA cell
  !    size; the *_ref scalars are reported in the diagnostic header and
  !    consumed by grid_create_clump for grid-array fill-in.
  cl_radius_max      = maxval(cl_radius(1:ncl))
  cl_temperature_ref = cl_temperature(1)
  cl_vtherm_ref      = cl_vtherm(1)
  cl_Dfreq_ref       = cl_Dfreq(1)
  cl_voigt_a_ref     = cl_voigt_a(1)
  cl_rhokap_ref      = cl_rhokap(1)
  base_rhokap_in     = cl_rhokap_ref

  if (mpar%p_rank == 0) then
     write(*,'(a,i14)')    ' Clumps: N_clumps  = ', N_clumps
     write(*,'(a,es12.4)') ' Clumps: cl_rhokap = ', cl_rhokap_ref
     write(*,'(a,f12.5)')  ' Clumps: voigt_a   = ', cl_voigt_a_ref
     write(*,'(a,es12.4)') ' Clumps: cl_Dfreq  = ', cl_Dfreq_ref
     write(*,'(a,2es12.4)')' Clumps: r_cl min/max  = ', minval(cl_radius(1:ncl)), cl_radius_max
  end if

  call build_clump_csr()
  call MPI_BARRIER(mpar%hostcomm, ierr)

  end subroutine read_clumps_fits
  !===========================================================================

end module clump_mod
