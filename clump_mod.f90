module clump_mod
!---------------------------------------------------------------------------
! Clumpy medium support for LaRT_v2.00.
!
! Geometry: N_clumps spherical clumps of (per-clump) radius cl_radius(:)
! placed uniformly at random (non-overlapping via RSA) inside a sphere of
! radius sphere_R. Each clump has its own opacity cl_rhokap(:), temperature
! cl_temperature(:), Doppler frequency cl_Dfreq(:), Voigt parameter
! cl_voigt_a(:), and independent Gaussian random bulk velocities.
! When all radial profiles are 'constant', every per-clump entry equals
! the corresponding reference scalar; otherwise the arrays are populated
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
  integer(int64), save :: N_clumps    = 0_int64
  real(kind=wp),  save :: sphere_R    = 0.0_wp   ! outer sphere radius [code units]
  real(kind=wp),  save :: r_min_clump = 0.0_wp   ! inner placement radius [code units]
                                                 ! (= max(0, par%rmin); 0 -> filled sphere)

  !--- Per-clump physical properties (MPI shared memory, dimension N_clumps).
  !    When all radial profiles are 'constant', every entry is uniform and
  !    equals the corresponding _ref scalar below; when a radial profile is
  !    active, they are populated as functions of clump-center radius.
  real(kind=wp), pointer, save :: cl_radius(:)      => null()  ! clump radius [code units]
  real(kind=wp), pointer, save :: cl_radius2(:)     => null()  ! cl_radius(icl)**2
  real(kind=wp), pointer, save :: cl_rhokap(:)      => null()  ! opacity/code-unit inside clump
  real(kind=wp), pointer, save :: cl_rhokapD(:)     => null()  ! dust opacity/code-unit (allocated only if DGR>0)
  real(kind=wp), pointer, save :: cl_voigt_a(:)     => null()  ! Voigt damping parameter
  real(kind=wp), pointer, save :: cl_Dfreq(:)       => null()  ! Doppler frequency [Hz]
  real(kind=wp), pointer, save :: cl_vtherm(:)      => null()  ! thermal velocity [km/s]
  real(kind=wp), pointer, save :: cl_temperature(:) => null()  ! clump temperature [K]

  !--- Reference / representative scalars used during initialisation, grid
  !    setup, and FITS-header reporting. In the uniform case these equal
  !    every cl_*(icl) entry. cl_radius_max is also used to size the RSA /
  !    CSR acceleration grids so that the 27-neighbor overlap search
  !    remains complete when r_cl varies between clumps.
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
  !    is needed. write_clumps_info() multiplies by cl_vtherm(icl) to keep the
  !    user-facing FITS output in km/s.
  real(kind=dp), pointer, save :: cl_x(:)  => null()
  real(kind=dp), pointer, save :: cl_y(:)  => null()
  real(kind=dp), pointer, save :: cl_z(:)  => null()
  real(kind=dp), pointer, save :: cl_vx(:) => null()
  real(kind=dp), pointer, save :: cl_vy(:) => null()
  real(kind=dp), pointer, save :: cl_vz(:) => null()

  !--- Radial-profile machinery.
  !    base_* values are the user-supplied peak values that the shape factors
  !    multiply. r_cl(r) = base_radius_in * shape_radius(r), etc.
  !    Profiles defined on r in [0, sphere_R].
  real(kind=wp), save :: base_radius_in   = 0.0_wp     ! par%clump_radius
  real(kind=wp), save :: base_rhokap_in   = 0.0_wp     ! density-equivalent peak rhokap
  real(kind=wp), save :: base_nH_in       = 0.0_wp     ! peak n_H [cm^-3] (when known)
  logical,       save :: profiles_active  = .false.    ! .true. if any profile != 'constant'
  logical,       save :: clumps_from_file = .false.    ! .true. if init_clumps loaded the population
  !--- When .true., kappa_clump returns the GAS line opacity only (no dust).
  !    Set transiently by the raytrace_to_*_tau_gas_clump wrappers so that
  !    gas-only sightline tau matches the Cartesian _tau_gas routines.
  !    Safe as a module scalar because v2.00 is MPI-only (one thread/rank).
  logical,       save :: clump_gas_only   = .false.
                                                       !  from par%clump_input_file (via read_clumps_info)
  logical,       save :: has_overlap      = .false.    ! .true. if file-loaded clumps contain overlapping pairs

  !--- Tabulated radial CDF for inverse-CDF sampling of clump positions.
  integer, parameter :: NPROF = 4001                   ! 1-D radial table size
  real(kind=wp), save :: prof_dr     = 0.0_wp
  real(kind=wp), save :: prof_r(NPROF)
  real(kind=wp), save :: prof_shape_number(NPROF)      ! n_cl(r) shape (unnormalized)
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
  ! Per-clump unit-conversion helpers.  photon%xfreq and the Jout / Jin / peel
  ! arrays use REF Doppler units (cl_Dfreq_ref) globally.  cl_voigt_a(icl) is
  ! defined with cl_Dfreq(icl), and cl_vx/y/z(icl) is stored as v/cl_vtherm(icl).
  ! When cl_Dfreq(icl) /= cl_Dfreq_ref (per-clump T), the Voigt argument and
  ! bulk-velocity Doppler shift need rescaling.  For uniform T both helpers
  ! reduce to the previous expressions (ratio = 1).
  !===========================================================================
  function voigt_clump(xfreq, icl) result(v)
  real(kind=wp),  intent(in) :: xfreq
  integer(int64), intent(in) :: icl
  real(kind=wp)              :: v, xloc, Dnu, a_ratio, f_ratio
  integer                    :: i
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt
  !--- Multiplet profile (exact mirror of calc_voigt3) in this clump's LOCAL
  !    Doppler units.  photon%xfreq is carried in REF units, so rescale once
  !    to local units (xloc); the member offsets delE_Hz(i)/cl_Dfreq(icl) are
  !    likewise in local units.  For single-line types (nup=1) the loop is
  !    skipped -> identical to the previous single-voigt behaviour.
  xloc = xfreq * (cl_Dfreq_ref/cl_Dfreq(icl))
  v    = voigt(xloc, cl_voigt_a(icl))
  do i = 2, line%nup
     Dnu     = line%delE_Hz(i)   / cl_Dfreq(icl)
     a_ratio = line%b(i)%damping / line%b(1)%damping
     f_ratio = line%f12(i)       / line%f12(1)
     v = v + voigt(xloc + Dnu, cl_voigt_a(icl)*a_ratio) * f_ratio
  end do
  end function voigt_clump
  !===========================================================================
  ! Total opacity per code-unit length for a clump leaf: gas line opacity
  ! (multiplet-aware via voigt_clump) plus the co-located dust continuum
  ! opacity when DGR > 0.  Mirrors the Cartesian raytrace
  !   rhokap = grid%rhokap*calc_voigt(...) + grid%rhokapD
  ! so that the dust/resonance split and total tau are consistent.
  !===========================================================================
  function kappa_clump(xfreq, icl) result(kap)
  real(kind=wp),  intent(in) :: xfreq
  integer(int64), intent(in) :: icl
  real(kind=wp)              :: kap
  kap = cl_rhokap(icl) * voigt_clump(xfreq, icl)
  if (par%DGR > 0.0_wp .and. .not. clump_gas_only) kap = kap + cl_rhokapD(icl)
  end function kappa_clump
  !===========================================================================
  pure function ulos_clump(icl, kx, ky, kz) result(u)
  integer(int64), intent(in) :: icl
  real(kind=wp),  intent(in) :: kx, ky, kz
  real(kind=wp)              :: u
  u = (real(cl_vx(icl),wp)*kx + real(cl_vy(icl),wp)*ky + real(cl_vz(icl),wp)*kz) &
      * (cl_Dfreq(icl)/cl_Dfreq_ref)
  end function ulos_clump

  !===========================================================================
  pure integer function cg_cell_idx(i, j, k)
  integer, intent(in) :: i, j, k
  cg_cell_idx = 1 + i + j*cgx + k*cgx*cgy
  end function cg_cell_idx
  !===========================================================================

  !===========================================================================
  ! Radial-profile evaluator. Returns the multiplicative shape factor
  ! associated with a profile name at radius r (in code units, 0 <= r <= R_box).
  ! The shape is normalized so that the user-supplied base value
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
  case ('powerlaw', 'power_law')
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
  real(kind=wp) :: r0_use
  r0_use = par%clump_radius_r0
  if (r0_use <= 0.0_wp) r0_use = sphere_R
  f = profile_shape(par%clump_radius_profile,  par%clump_radius_alpha, &
                    r0_use, r, 1)
  end function shape_radius

  pure real(kind=wp) function shape_density(r) result(f)
  real(kind=wp), intent(in) :: r
  real(kind=wp) :: r0_use
  r0_use = par%clump_density_r0
  if (r0_use <= 0.0_wp) r0_use = sphere_R
  f = profile_shape(par%clump_density_profile, par%clump_density_alpha, &
                    r0_use, r, 2)
  end function shape_density

  pure real(kind=wp) function shape_number(r) result(f)
  real(kind=wp), intent(in) :: r
  real(kind=wp) :: r0_use
  r0_use = par%clump_number_r0
  if (r0_use <= 0.0_wp) r0_use = sphere_R
  f = profile_shape(par%clump_number_profile,  par%clump_number_alpha, &
                    r0_use, r, 3)
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

  ! Suppress the radial number profile inside the inner cavity so that
  ! position sampling, total counts, and integrated f_vol/f_cov all see
  ! n_cl(r) = 0 for r < r_min_clump.  shape_radius/shape_density are kept
  ! intact for diagnostic output; they have no effect when shape_number = 0.
  if (prof_r(1) < r_min_clump) prof_shape_number(1) = 0.0_wp

  if (prof_shape_number(1) > 0.0_wp) then
     rcl_local = base_radius_in * prof_shape_radius(1)
     if (rcl_local > cl_radius_max) cl_radius_max = rcl_local
  end if

  last_integrand = prof_shape_number(1) * 0.0_wp     ! r^2 = 0 at r=0

  do i = 2, NPROF
     r = real(i-1, wp) * prof_dr
     prof_r(i)             = r
     prof_shape_number(i)  = shape_number(r)
     prof_shape_radius(i)  = shape_radius(r)
     prof_shape_density(i) = shape_density(r)
     if (r < r_min_clump) prof_shape_number(i) = 0.0_wp
     integrand             = prof_shape_number(i) * r * r
     prof_cdf_pos(i) = prof_cdf_pos(i-1) + 0.5_wp * (last_integrand + integrand) * prof_dr
     last_integrand = integrand

     if (prof_shape_number(i) > 0.0_wp) then
        rcl_local = base_radius_in * prof_shape_radius(i)
        if (rcl_local > cl_radius_max) cl_radius_max = rcl_local
     end if
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
  ! Sample a clump-center radius from the tabulated inverse CDF.
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
  ! of clump intersections along a radial sightline from the center to the
  ! sphere surface, evaluated for the *current* normalization A_norm:
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
  real(kind=wp) :: rcl, r, integrand(NPROF), V_shell
  integer :: i
  do i = 1, NPROF
     r   = prof_r(i)
     rcl = base_radius_in * prof_shape_radius(i)
     integrand(i) = A_norm * prof_shape_number(i) * rcl**3 * r * r
  end do
  V_shell = max(R_box**3 - r_min_clump**3, tiny(1.0_wp))
  fvol = fourpi * integrate_table(integrand) / V_shell
  end function f_vol_quad
  !===========================================================================

  !===========================================================================
  ! Total clump count for normalization A_norm:
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
  ! Geometric factor GF for the radial-sightline taumax / N_HImax integrals
  ! through the clumpy medium:
  !
  !   GF(A_norm) = (4*pi/3) * A_norm * base_radius_in^3
  !                * integral_rmin^R shape_number(r) * shape_radius(r)^3
  !                                * shape_density(r) dr
  !
  ! Given a uniform reference temperature (cl_Dfreq_ref, cl_voigt_a_ref) the
  ! line-center optical depth and HI column density, integrated radially from
  ! the box center to the outer surface, are
  !
  !   taumax  = GF(A_norm) * base_rhokap * voigt(0, cl_voigt_a_ref)
  !   N_HImax = GF(A_norm) * base_rhokap * cl_Dfreq_ref / line%cross0
  !
  ! These match the uniform closed-form taumax = (4/3) * f_cov * clump_tau0
  ! and N_HImax = (4/3) * f_cov * clump_NHI when shape_* = 1 on [rmin, R].
  ! The leading (4*pi/3) drops to f_vol_local(r) on combination with the
  ! shape_radius^3 factor; the remaining shape_density(r) carries the local
  ! n_HI scaling. Profile mode only -- the uniform path uses the closed
  ! form directly.
  !===========================================================================
  pure real(kind=wp) function LOS_geometric_quad(A_norm) result(GF)
  real(kind=wp), intent(in) :: A_norm
  real(kind=wp) :: integrand(NPROF)
  integer :: i
  do i = 1, NPROF
     integrand(i) = prof_shape_number(i) * prof_shape_radius(i)**3 * prof_shape_density(i)
  end do
  GF = (fourpi / 3.0_wp) * A_norm * base_radius_in**3 * integrate_table(integrand)
  end function LOS_geometric_quad
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
  real(kind=wp) :: temp_cl, vtherm, A_norm, fvol_unit, fcov_unit, fvol_realized, fcov_realized
  real(kind=wp) :: GF_los
  integer       :: ierr

  sphere_R       = R_sphere
  base_radius_in = par%clump_radius
  if (base_radius_in <= 0.0_wp) then
     if (mpar%p_rank == 0) write(*,*) 'ERROR: clump_radius must be > 0'
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end if

  !--- Inner placement radius. par%rmin = -999 (default) is treated as "no
  !    inner cavity". Negative values map to 0; rmin >= rmax is fatal.
  r_min_clump = max(0.0_wp, par%rmin)
  if (r_min_clump >= R_sphere) then
     if (mpar%p_rank == 0) write(*,'(a,2es12.4)') &
        'ERROR: par%rmin must be < par%rmax for clump placement; rmin, rmax = ', &
        r_min_clump, R_sphere
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end if

  !--- If the user supplied a pre-built clump FITS file, load it and skip the
  !    internal profile / RSA generation entirely. Sets cl_*, sphere_R,
  !    cl_radius_max and reference scalars; builds the CSR acceleration grid.
  if (len_trim(par%clump_input_file) > 0) then
     call read_clumps_info(trim(par%clump_input_file), R_sphere)
     call check_has_overlap()
     return
  end if

  !--- Detect any active radial profile.
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
  vtherm             = vtherm_total(temp_cl)
  cl_Dfreq_ref       = vtherm / (line%wavelength0 * um2km)
  cl_vtherm_ref      = vtherm                          ! km/s
  cl_temperature_ref = temp_cl                         ! [K]
  cl_voigt_a_ref     = (line%damping / fourpi) / cl_Dfreq_ref

  !--- derive N_clumps + f_vol/f_cov_realized first; opacity may need the
  !    realized population to back-solve from system-level taumax / N_HImax.
  !    None of the helpers below reference cl_rhokap_ref.
  !
  !    Uniform case: closed-form using shell volume V_shell = (4pi/3)(R^3 -
  !    rmin^3) and shell sightline factor (R^2 + R*rmin + rmin^2). With
  !    r_min_clump = 0 these reduce to the original full-sphere formulas.
  !    Profile case: numerical quadrature with the LOS f_cov definition
  !    (see docs/covering_factor_definitions.tex).
  A_norm = 0.0_wp
  if (.not. profiles_active) then
     if (par%clump_N_clumps > 0.0_wp) then
        N_clumps = int(par%clump_N_clumps, int64)
     else if (par%clump_f_vol > 0.0_wp) then
        N_clumps = nint(par%clump_f_vol * (R_sphere**3 - r_min_clump**3) / cl_radius_max**3, int64)
     else if (par%clump_f_cov > 0.0_wp) then
        N_clumps = nint((4.0_wp/3.0_wp)*par%clump_f_cov * &
                        (R_sphere**2 + R_sphere*r_min_clump + r_min_clump**2) / &
                        cl_radius_max**2, int64)
     else
        if (mpar%p_rank == 0) &
           write(*,*) 'ERROR: specify clump_N_clumps, clump_f_vol, or clump_f_cov'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
     fvol_realized = real(N_clumps, wp) * cl_radius_max**3 / &
                     max(R_sphere**3 - r_min_clump**3, tiny(1.0_wp))
     fcov_realized = 0.75_wp * real(N_clumps, wp) * cl_radius_max**2 / &
                     max(R_sphere**2 + R_sphere*r_min_clump + r_min_clump**2, tiny(1.0_wp))
  else
     !--- non-uniform: integrate the shape and back-solve for the density
     !    normalization A_norm of n_cl(r) = A_norm * shape_number(r).
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
     fvol_realized = f_vol_quad(A_norm, R_sphere)
     fcov_realized = f_cov_LOS_quad(A_norm)
  end if
  if (N_clumps <= 0_int64) N_clumps = 1_int64

  !--- clump opacity (peak value).
  !    Per-clump direct inputs (clump_tau0, clump_NHI, clump_nH) take
  !    priority over the system-level fallbacks (par%taumax, par%N_HImax),
  !    which back-solve for the peak opacity assuming the realized f_cov
  !    (uniform) or shape quadrature (profile) hits the requested radial
  !    sightline target.
  base_nH_in = 0.0_wp
  if (par%clump_tau0 > 0.0_wp) then
     cl_rhokap_ref = par%clump_tau0 / (voigt(0.0_wp, cl_voigt_a_ref) * base_radius_in)
  else if (par%clump_NHI > 0.0_wp) then
     ! clump_NHI = per-clump column density [cm^-2] from clump center to surface
     ! (peak value; shape_density(r) modulates per clump).
     cl_rhokap_ref = par%clump_NHI * line%cross0 / (cl_Dfreq_ref * base_radius_in)
  else if (par%clump_nH > 0.0_wp) then
     if (par%distance2cm <= 0.0_wp) then
        if (mpar%p_rank == 0) write(*,*) 'ERROR: clump_nH requires distance_unit'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
     cl_rhokap_ref = par%clump_nH * line%cross0 * par%distance2cm / cl_Dfreq_ref
     base_nH_in    = par%clump_nH
  else if (par%taumax > 0.0_wp .or. par%N_HImax > 0.0_wp) then
     !--- Back-solve from system-level radial sightline target.
     !    For uniform: GF = N * r_cl^3 / (R^2 + R*rmin + rmin^2)
     !    For profile: GF = LOS_geometric_quad(A_norm)
     !    Then taumax = GF * cl_rhokap_ref * voigt0   or
     !         N_HImax = GF * cl_rhokap_ref * cl_Dfreq_ref / cross0.
     if (profiles_active) then
        GF_los = LOS_geometric_quad(A_norm)
     else
        GF_los = real(N_clumps, wp) * cl_radius_max**3 / &
                 max(R_sphere**2 + R_sphere*r_min_clump + r_min_clump**2, tiny(1.0_wp))
     end if
     if (GF_los <= 0.0_wp) then
        if (mpar%p_rank == 0) write(*,*) &
           'ERROR: cannot back-solve clump opacity (geometric factor is zero)'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
     if (par%taumax > 0.0_wp) then
        cl_rhokap_ref = par%taumax / (GF_los * voigt(0.0_wp, cl_voigt_a_ref))
     else
        cl_rhokap_ref = par%N_HImax * line%cross0 / (GF_los * cl_Dfreq_ref)
     end if
  else
     if (mpar%p_rank == 0) write(*,*) &
        'ERROR: specify clump_tau0, clump_NHI, clump_nH, taumax, or N_HImax'
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end if
  base_rhokap_in = cl_rhokap_ref

  if (mpar%p_rank == 0) then
     write(*,'(a,i14)')    ' Clumps: N_clumps  = ', N_clumps
     write(*,'(a,f12.6)')  ' Clumps: f_vol     = ', fvol_realized
     write(*,'(a,f12.5)')  ' Clumps: f_cov     = ', fcov_realized
     write(*,'(a,2f12.5)') ' Clumps: rmin/rmax = ', r_min_clump, R_sphere
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
  if (par%DGR > 0.0_wp) call create_shared_mem(cl_rhokapD, [int(N_clumps)])
  call create_shared_mem(cl_voigt_a,     [int(N_clumps)])
  call create_shared_mem(cl_Dfreq,       [int(N_clumps)])
  call create_shared_mem(cl_vtherm,      [int(N_clumps)])
  call create_shared_mem(cl_temperature, [int(N_clumps)])

  !--- Fill per-clump arrays uniformly with the reference values. When
  !    radial profiles are active, the profile-driven assignments below
  !    overwrite these defaults once clump positions are placed.
  if (mpar%h_rank == 0) then
     cl_radius(:)      = cl_radius_max
     cl_radius2(:)     = cl_radius_max * cl_radius_max
     cl_rhokap(:)      = cl_rhokap_ref
     cl_voigt_a(:)     = cl_voigt_a_ref
     cl_Dfreq(:)       = cl_Dfreq_ref
     cl_vtherm(:)      = cl_vtherm_ref
     cl_temperature(:) = cl_temperature_ref
     !--- uniform-case dust opacity (overwritten per-clump below if profiles
     !    are active), matching the Cartesian rhokapD/rhokap ratio.
     if (par%DGR > 0.0_wp .and. associated(cl_rhokapD)) &
        cl_rhokapD(:) = cl_rhokap_ref * par%cext_dust * par%DGR * cl_Dfreq_ref / line%cross0
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
     ! cl_vx/y/z are already stored as v / cl_vtherm(icl) (normalized inline
     ! by generate_clumps and assign_clump_velocities_from_type).
  end if
  call MPI_BARRIER(mpar%hostcomm, ierr)

  !--- h_rank=0 builds CSR grid; barrier before use
  call build_clump_csr()
  call MPI_BARRIER(mpar%hostcomm, ierr)

  !--- For internally generated clumps with overlap allowed, run overlap
  !    detection so the overlap-aware raytrace path is engaged if needed.
  if (par%clump_allow_overlap) call check_has_overlap()

  end subroutine init_clumps
  !===========================================================================

  !===========================================================================
  subroutine generate_clumps()
  !---------------------------------------------------------------------------
  ! RSA with linked-list grid acceleration. Called only on h_rank=0.
  ! Two sampling paths:
  !  - profiles_active = .false.: uniform-in-sphere box rejection.
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
  real(kind=wp)  :: r_min_center, r_max_center, r_min_center2, r_max_center2
  real(kind=wp)  :: cos_cone_opening, r_trial2, cos_theta_cone, cos_theta_min
  integer        :: ig, jg, kg, ig2, jg2, kg2, icell_rsa, jnb
  integer        :: ierr
  integer(int64) :: icl, n_attempts
  logical        :: overlap
  integer, allocatable :: head(:), nxt(:)

  rg         = min(512, max(32, int(real(N_clumps,wp)**(1.0_wp/3.0_wp)) + 1))
  ncells_rsa = rg**3
  ! RSA grid spans the bounding sphere; cells must be at least 2*cl_radius_max
  ! across so the 27-neighbor search captures every possible overlap.
  rg_cell    = max(2.0_wp * sphere_R / real(rg, wp), 2.0_wp * cl_radius_max)
  rg         = max(2, int((2.0_wp * sphere_R) / rg_cell) + 1)
  rg_cell    = (2.0_wp * sphere_R) / real(rg, wp)
  ncells_rsa = rg**3
  min_sep2_uniform = (2.0_wp * cl_radius_max)**2

  !--- "fully-inside" mode: a clump is accepted only if it lies entirely in
  !    the radial shell [r_min_clump, sphere_R], i.e. both
  !    r_center + R_clump <= sphere_R   (no protrusion past the outer edge)
  !    r_center - R_clump >= r_min_clump (no protrusion into the inner cavity).
  !    This requires r_min_clump + 2*cl_radius_max <= sphere_R; abort otherwise.
  if (par%clump_fully_inside) then
     if (cl_radius_max >= sphere_R) then
        if (mpar%p_rank == 0) write(*,'(a)') &
           'ERROR: par%clump_fully_inside=.true. but max clump radius >= sphere_R; '// &
           'cannot fit any clump entirely inside the medium.'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
     if (r_min_clump + 2.0_wp*cl_radius_max > sphere_R) then
        if (mpar%p_rank == 0) write(*,'(a,3es12.4)') &
           'ERROR: par%clump_fully_inside=.true. but rmin + 2*cl_radius_max > rmax; '// &
           'no clump fits inside the shell. rmin, cl_radius_max, rmax = ', &
           r_min_clump, cl_radius_max, sphere_R
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
  end if

  allocate(head(ncells_rsa),    source=-1)
  allocate(nxt(int(N_clumps)), source=-1)

  if (mpar%p_rank == 0) then
     write(*,'(a,i5,a,i14,a)') ' RSA: grid ', rg, '^3, placing ', N_clumps, ' clumps...'
     if (par%clump_fully_inside) then
        write(*,'(a)') ' RSA: clump_fully_inside = .true.  (clumps must fit inside the shell)'
     else
        write(*,'(a)') ' RSA: clump_fully_inside = .false. (only centers inside the shell)'
     end if
     if (r_min_clump > 0.0_wp) &
        write(*,'(a,2f12.5)') ' RSA: shell rmin/rmax    = ', r_min_clump, sphere_R
  end if

  n_attempts = 0_int64
  icl        = 0_int64
  rcl_trial  = base_radius_in   ! default for uniform path

  !--- Precompute cone angle for biconical geometry
  if (par%cone_opening > 0.0_wp .and. par%cone_opening < 90.0_wp) then
     cos_cone_opening = cos(par%cone_opening * deg2rad)
  else
     cos_cone_opening = -1.0_wp   ! accept everything (full sphere)
  end if

  do while (icl < N_clumps)
     n_attempts = n_attempts + 1_int64

     if (profiles_active) then
        !--- inverse-CDF radial draw + isotropic angles. The CDF is built so
        !    that prof_shape_number(r) = 0 for r < r_min_clump, so the draw
        !    naturally avoids the inner cavity even when clump_fully_inside
        !    is .false.
        r_trial    = sample_clump_radius()
        rcl_trial  = base_radius_in * shape_radius(r_trial)
        !--- Direct sampling of angle within the cone (or full sphere)
        if (cos_cone_opening > 0.0_wp) then
           cos_theta = cos_cone_opening + (1.0_wp - cos_cone_opening) * rand_number()
           if (rand_number() < 0.5_wp) cos_theta = -cos_theta
        else
           cos_theta = 2.0_wp * rand_number() - 1.0_wp
        end if
        sin_theta  = sqrt(max(0.0_wp, 1.0_wp - cos_theta*cos_theta))
        phi_az     = twopi * rand_number()
        xc = r_trial * sin_theta * cos(phi_az)
        yc = r_trial * sin_theta * sin(phi_az)
        zc = r_trial * cos_theta
        !--- "fully-inside" rejection: retry if the clump would protrude
        !    past the outer sphere boundary, into the inner cavity,
        !    or outside the cone boundary.
        if (par%clump_fully_inside) then
           if (r_trial + rcl_trial > sphere_R)  cycle
           if (r_trial - rcl_trial < r_min_clump) cycle
           !--- Cone fully-inside: reject if the clump protrudes outside
           !    the cone. The angular radius of the clump seen from the
           !    origin is delta = asin(rcl/r). The minimum |cos(theta)|
           !    over the clump sphere is cos(theta_c + delta).
           if (cos_cone_opening > 0.0_wp .and. r_trial > rcl_trial) then
              cos_theta_min = abs(cos_theta) * sqrt(1.0_wp - (rcl_trial/r_trial)**2) &
                            - sin_theta * (rcl_trial / r_trial)
              if (cos_theta_min < cos_cone_opening) cycle
           end if
        end if
     else
        !--- uniform random point inside the placement shell (box-rejection
        !    sampling). When clump_fully_inside is set, the shell is inset
        !    by one base clump radius on both sides so the entire clump
        !    fits inside the medium.
        if (par%clump_fully_inside) then
           r_max_center = sphere_R    - base_radius_in
           r_min_center = r_min_clump + base_radius_in
        else
           r_max_center = sphere_R
           r_min_center = r_min_clump
        end if
        r_max_center2 = r_max_center * r_max_center
        r_min_center2 = r_min_center * r_min_center
        if (cos_cone_opening > 0.0_wp) then
           !--- Direct sampling inside the bicone
           do
              r_trial   = (r_min_center**3 + (r_max_center**3 - r_min_center**3) * rand_number())**(1.0_wp/3.0_wp)
              cos_theta = cos_cone_opening + (1.0_wp - cos_cone_opening) * rand_number()
              if (rand_number() < 0.5_wp) cos_theta = -cos_theta
              sin_theta = sqrt(max(0.0_wp, 1.0_wp - cos_theta*cos_theta))
              !--- fully-inside cone check
              if (par%clump_fully_inside .and. r_trial > base_radius_in) then
                 cos_theta_min = abs(cos_theta) * sqrt(1.0_wp - (base_radius_in/r_trial)**2) &
                               - sin_theta * (base_radius_in / r_trial)
                 if (cos_theta_min < cos_cone_opening) cycle
              end if
              exit
           end do
           phi_az    = twopi * rand_number()
           xc = r_trial * sin_theta * cos(phi_az)
           yc = r_trial * sin_theta * sin(phi_az)
           zc = r_trial * cos_theta
        else
           !--- Standard box-rejection for full sphere
           do
              xc = (2.0_wp * rand_number() - 1.0_wp) * r_max_center
              yc = (2.0_wp * rand_number() - 1.0_wp) * r_max_center
              zc = (2.0_wp * rand_number() - 1.0_wp) * r_max_center
              d2 = xc*xc + yc*yc + zc*zc
              if (d2 <= r_max_center2 .and. d2 >= r_min_center2) exit
           end do
        end if
     end if

     !--- cell in RSA grid (0-based)
     ig = min(rg-1, max(0, int((xc + sphere_R) / rg_cell)))
     jg = min(rg-1, max(0, int((yc + sphere_R) / rg_cell)))
     kg = min(rg-1, max(0, int((zc + sphere_R) / rg_cell)))

     !--- check 27 neighbors for overlap
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
     if (overlap .and. .not. par%clump_allow_overlap) cycle

     !--- accept
     icl = icl + 1_int64
     cl_x(icl) = real(xc, dp);  cl_y(icl) = real(yc, dp);  cl_z(icl) = real(zc, dp)

     if (profiles_active) then
        !--- per-clump physical assignments from radial profiles
        cl_radius(icl)     = rcl_trial
        cl_radius2(icl)    = rcl_trial * rcl_trial
        if (allocated(tab_temperature)) then
           temp_loc = profile_file_interp(r_trial, 4)
           if (temp_loc <= 0.0_wp) temp_loc = cl_temperature_ref
        else
           temp_loc = cl_temperature_ref
        end if
        cl_temperature(icl) = temp_loc
        vth_loc = vtherm_total(temp_loc)
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
        if (par%DGR > 0.0_wp .and. associated(cl_rhokapD)) &
           cl_rhokapD(icl) = kap_loc * par%cext_dust * par%DGR * Df_loc / line%cross0
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

  if (mpar%p_rank == 0) then
     if (par%clump_allow_overlap) then
        write(*,'(a)') ' Random placement done (overlap allowed).'
     else
        write(*,'(a,f6.1,a)') ' RSA done, acceptance rate = ', &
           real(N_clumps,wp)/real(n_attempts,wp)*100.0_wp, '%'
     end if
  end if

  deallocate(head, nxt)
  end subroutine generate_clumps
  !===========================================================================

  !===========================================================================
  subroutine assign_clump_velocities_from_type()
  !---------------------------------------------------------------------------
  ! Add a systematic velocity component to each clump based on par%velocity_type
  ! and the clump's center position.  Called only on h_rank=0 after
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

     case ('power_law')
        ! v(r) = Vexp * (r / sphere_R)^velocity_alpha  [km/s]
        if (rr > 0.0_wp) then
           Vscale = par%Vexp * (rr / sphere_R)**par%velocity_alpha
           vx = Vscale * xc / rr
           vy = Vscale * yc / rr
           vz = Vscale * zc / rr
        end if

     case ('linear_decelerate')
        ! v(r) = Vexp * (sphere_R - r) / (sphere_R - rmin)  [km/s]
        if (rr > 0.0_wp) then
           Vscale = par%Vexp * max(0.0_wp, (sphere_R - rr) / (sphere_R - max(par%rmin, 0.0_wp)))
           vx = Vscale * xc / rr
           vy = Vscale * yc / rr
           vz = Vscale * zc / rr
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
  ! Sphere: center (cx,cy,cz), squared-radius cl_radius2(icl).
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
  subroutine check_has_overlap()
  !---------------------------------------------------------------------------
  ! Scan all clump pairs using the CSR grid to detect any overlapping pair.
  ! Sets has_overlap = .true. on rank 0 if found, then broadcasts.
  ! O(N * k_neighbors) time; called once at init after read_clumps_info.
  !---------------------------------------------------------------------------
  implicit none
  integer        :: imin, imax, jmin, jmax, kmin, kmax
  integer        :: i, j, k, icell, ip
  integer(int64) :: icl_i, icl_j
  real(kind=dp)  :: dx, dy, dz, dist2, r_sum
  integer        :: ierr

  has_overlap = .false.
  if (mpar%h_rank == 0) then
     outer: do icl_i = 1_int64, N_clumps
        call clump_cell_range(icl_i, imin, imax, jmin, jmax, kmin, kmax)
        do k = kmin, kmax
           do j = jmin, jmax
              do i = imin, imax
                 icell = cg_cell_idx(i, j, k)
                 do ip = cg_start(icell), cg_start(icell+1) - 1
                    icl_j = int(cg_list(ip), int64)
                    if (icl_j <= icl_i) cycle
                    dx    = cl_x(icl_j) - cl_x(icl_i)
                    dy    = cl_y(icl_j) - cl_y(icl_i)
                    dz    = cl_z(icl_j) - cl_z(icl_i)
                    dist2 = dx*dx + dy*dy + dz*dz
                    r_sum = real(cl_radius(icl_i) + cl_radius(icl_j), dp)
                    if (dist2 < r_sum*r_sum) then
                       has_overlap = .true.
                       exit outer
                    end if
                 end do
              end do
           end do
        end do
     end do outer
  end if
  call MPI_BCAST(has_overlap, 1, MPI_LOGICAL, 0, mpar%SAME_HRANK_COMM, ierr)
  if (mpar%p_rank == 0) then
     if (has_overlap) then
        write(*,'(a)') ' Overlap check: overlapping clumps detected — using overlap-aware raytrace.'
     else
        write(*,'(a)') ' Overlap check: no overlaps found — using standard raytrace.'
     end if
  end if
  end subroutine check_has_overlap
  !===========================================================================

  !===========================================================================
  subroutine active_set_at_point(xp, yp, zp, active, n_active)
  !---------------------------------------------------------------------------
  ! Return all clump indices containing point (xp,yp,zp).
  ! active(1:n_active) are the 1-based clump indices.
  ! The caller must provide active(:) with size >= some upper bound.
  !---------------------------------------------------------------------------
  implicit none
  real(kind=wp),  intent(in)  :: xp, yp, zp
  integer(int64), intent(out) :: active(:)
  integer,        intent(out) :: n_active

  integer        :: ci, cj, ck, i, j, k, icell, ip
  integer(int64) :: icl
  real(kind=dp)  :: rx, ry, rz

  n_active = 0
  ci = max(0, min(cgx-1, int((xp - cg_xmin) * cg_inv_dx)))
  cj = max(0, min(cgy-1, int((yp - cg_ymin) * cg_inv_dy)))
  ck = max(0, min(cgz-1, int((zp - cg_zmin) * cg_inv_dz)))

  do k = max(0, ck-1), min(cgz-1, ck+1)
     do j = max(0, cj-1), min(cgy-1, cj+1)
        do i = max(0, ci-1), min(cgx-1, ci+1)
           icell = cg_cell_idx(i, j, k)
           do ip = cg_start(icell), cg_start(icell+1) - 1
              icl = int(cg_list(ip), int64)
              rx  = real(xp, dp) - cl_x(icl)
              ry  = real(yp, dp) - cl_y(icl)
              rz  = real(zp, dp) - cl_z(icl)
              if (rx*rx + ry*ry + rz*rz <= real(cl_radius2(icl), dp)) then
                 !--- avoid duplicates (a clump can register in multiple CSR cells)
                 if (n_active == 0 .or. all(active(1:n_active) /= icl)) then
                    n_active = n_active + 1
                    active(n_active) = icl
                 end if
              end if
           end do
        end do
     end do
  end do
  end subroutine active_set_at_point
  !===========================================================================

  !===========================================================================
  subroutine collect_ray_events_overlap(xp, yp, zp, kx, ky, kz, t_max, &
                                         ev_t, ev_icl, ev_type, n_ev)
  !---------------------------------------------------------------------------
  ! Full DDA scan: collect every clump entry/exit event along the ray
  !   P(t) = (xp,yp,zp) + t*(kx,ky,kz),   0 < t <= t_max.
  ! ev_type: +1 = ENTER, -1 = EXIT.
  ! Clumps that already contain the ray origin emit only an EXIT event.
  ! Events are returned sorted in ascending t order (insertion sort; n_ev small).
  ! The caller must provide arrays sized >= some upper bound (e.g. 2*N_clumps).
  !---------------------------------------------------------------------------
  implicit none
  real(kind=wp),  intent(in)  :: xp, yp, zp, kx, ky, kz, t_max
  real(kind=wp),  intent(out) :: ev_t(:)
  integer(int64), intent(out) :: ev_icl(:)
  integer,        intent(out) :: ev_type(:)
  integer,        intent(out) :: n_ev

  integer  :: ci, cj, ck, si, sj, sk, icell, ip, ie, ii
  integer(int64) :: icl
  real(kind=wp)  :: tx, ty, tz, delx, dely, delz, d
  real(kind=wp)  :: t_entry, t_exit
  real(kind=wp)  :: ev_t_tmp
  integer(int64) :: ev_icl_tmp
  integer        :: ev_type_tmp
  logical        :: hit

  n_ev = 0

  !--- DDA initialization (same pattern as find_next_clump)
  ci = max(0, min(cgx-1, int((xp - cg_xmin) * cg_inv_dx)))
  cj = max(0, min(cgy-1, int((yp - cg_ymin) * cg_inv_dy)))
  ck = max(0, min(cgz-1, int((zp - cg_zmin) * cg_inv_dz)))

  if (kx > 0.0_wp) then
     si   =  1;  delx = cg_dx / kx
     tx   = (cg_xmin + real(ci+1,wp)*cg_dx - xp) / kx
  else if (kx < 0.0_wp) then
     si   = -1;  delx = -cg_dx / kx
     tx   = (cg_xmin + real(ci,wp)*cg_dx - xp) / kx
  else
     si = 0;  delx = hugest;  tx = hugest
  end if
  if (ky > 0.0_wp) then
     sj   =  1;  dely = cg_dy / ky
     ty   = (cg_ymin + real(cj+1,wp)*cg_dy - yp) / ky
  else if (ky < 0.0_wp) then
     sj   = -1;  dely = -cg_dy / ky
     ty   = (cg_ymin + real(cj,wp)*cg_dy - yp) / ky
  else
     sj = 0;  dely = hugest;  ty = hugest
  end if
  if (kz > 0.0_wp) then
     sk   =  1;  delz = cg_dz / kz
     tz   = (cg_zmin + real(ck+1,wp)*cg_dz - zp) / kz
  else if (kz < 0.0_wp) then
     sk   = -1;  delz = -cg_dz / kz
     tz   = (cg_zmin + real(ck,wp)*cg_dz - zp) / kz
  else
     sk = 0;  delz = hugest;  tz = hugest
  end if
  d = 0.0_wp

  do
     if (d > t_max) exit
     if (ci < 0 .or. ci >= cgx .or. cj < 0 .or. cj >= cgy .or. ck < 0 .or. ck >= cgz) exit
     icell = cg_cell_idx(ci, cj, ck)
     do ip = cg_start(icell), cg_start(icell+1) - 1
        icl = int(cg_list(ip), int64)
        call ray_sphere_isect(xp, yp, zp, kx, ky, kz, &
             real(cl_x(icl),wp), real(cl_y(icl),wp), real(cl_z(icl),wp), &
             icl, t_entry, t_exit, hit)
        if (.not. hit) cycle
        if (t_exit <= 0.0_wp) cycle
        if (t_entry > t_max) cycle
        !--- avoid duplicate events (clump registered in multiple CSR cells)
        do ie = 1, n_ev
           if (ev_icl(ie) == icl) goto 10
        end do
        if (n_ev + 2 > size(ev_t)) goto 10
        if (t_entry > 0.0_wp) then
           !--- ENTER event
           n_ev = n_ev + 1
           ev_t(n_ev)    = t_entry
           ev_icl(n_ev)  = icl
           ev_type(n_ev) = +1
        end if
        !--- EXIT event (always)
        n_ev = n_ev + 1
        ev_t(n_ev)    = min(t_exit, t_max)
        ev_icl(n_ev)  = icl
        ev_type(n_ev) = -1
10      continue
     end do
     !--- advance DDA
     if (tx <= ty .and. tx <= tz) then
        d = tx;  tx = tx + delx;  ci = ci + si
     else if (ty <= tz) then
        d = ty;  ty = ty + dely;  cj = cj + sj
     else
        d = tz;  tz = tz + delz;  ck = ck + sk
     end if
  end do

  !--- insertion sort by ev_t (n_ev is typically small)
  do ie = 2, n_ev
     ev_t_tmp    = ev_t(ie)
     ev_icl_tmp  = ev_icl(ie)
     ev_type_tmp = ev_type(ie)
     ii = ie - 1
     do while (ii >= 1 .and. ev_t(ii) > ev_t_tmp)
        ev_t(ii+1)    = ev_t(ii)
        ev_icl(ii+1)  = ev_icl(ii)
        ev_type(ii+1) = ev_type(ii)
        ii = ii - 1
     end do
     ev_t(ii+1)    = ev_t_tmp
     ev_icl(ii+1)  = ev_icl_tmp
     ev_type(ii+1) = ev_type_tmp
  end do

  end subroutine collect_ray_events_overlap
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
  nullify(cl_radius, cl_radius2, cl_rhokap, cl_rhokapD, cl_voigt_a, cl_Dfreq, &
          cl_vtherm, cl_temperature)
  nullify(cg_start, cg_list)
  N_clumps = 0_int64
  end subroutine destroy_clumps
  !===========================================================================

  !===========================================================================
  subroutine write_clumps_info(fname)
  !---------------------------------------------------------------------------
  ! Write clump positions, velocities, and physical parameters to a FITS
  ! binary table.  Called only on p_rank == 0 after init_clumps().
  !
  ! HDU 1 (primary): empty image, all clump scalars as header keywords.
  ! HDU 2 (BinTable): columns X, Y, Z [code units], VX, VY, VZ [km/s].
  !---------------------------------------------------------------------------
  use define
  use iofile_mod
  implicit none
  character(len=*), intent(in) :: fname

  type(io_file_type) :: iofh
  integer :: status, bitpix
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
  !    Uniform case: closed-form expression.
  !    Profile case: 4π/3 * Σ r_cl^3 / V_shell and (LOS f_cov from quadrature).
  !    From-file case: per-clump sums (no profile tables available).
  !    All formulas use shell volume V_shell = (4π/3)(R^3 - rmin^3) and shell
  !    sightline factor (R^2 + R*rmin + rmin^2); rmin = 0 reduces them to the
  !    original full-sphere expressions.
  if (clumps_from_file) then
     f_vol_actual = 0.0_wp
     f_cov_actual = 0.0_wp
     do i = 1_int64, ncl
        f_vol_actual = f_vol_actual + cl_radius(i)**3
        f_cov_actual = f_cov_actual + cl_radius(i)**2
     end do
     f_vol_actual = f_vol_actual / max(sphere_R**3 - r_min_clump**3, tiny(1.0_wp))
     f_cov_actual = 0.75_wp * f_cov_actual / &
                    max(sphere_R**2 + sphere_R*r_min_clump + r_min_clump**2, tiny(1.0_wp))
  else if (.not. profiles_active) then
     f_vol_actual = real(ncl,wp) * cl_radius_max**3 / &
                    max(sphere_R**3 - r_min_clump**3, tiny(1.0_wp))
     f_cov_actual = 0.75_wp * real(ncl,wp) * cl_radius_max**2 / &
                    max(sphere_R**2 + sphere_R*r_min_clump + r_min_clump**2, tiny(1.0_wp))
  else
     !--- realized f_vol from sum of individual clump volumes
     f_vol_actual = 0.0_wp
     do i = 1_int64, ncl
        f_vol_actual = f_vol_actual + cl_radius(i)**3
     end do
     f_vol_actual = f_vol_actual / max(sphere_R**3 - r_min_clump**3, tiny(1.0_wp))
     !--- realized f_cov from radial quadrature (LOS integral)
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
     write(*,'(a,es12.4)') ' Clumps: realized f_vol      = ', f_vol_actual
     write(*,'(a,es12.4)') ' Clumps: realized f_cov      = ', f_cov_actual
     write(*,'(a,3es12.4)')' Clumps: r_cl min/mean/max   = ', rcl_min, rcl_mean, rcl_max
     write(*,'(a,es12.4)') ' Clumps: mean tau_per_clump  = ', tau_mean
     write(*,'(a,es12.4)') ' Clumps: integrated HI col   = ', total_HI_mass
  end if

  call io_open_new(iofh, trim(fname), status)
  if (status /= 0) then
     write(*,*) 'WARNING: write_clumps_info: cannot open ', trim(fname)
     return
  end if

  !--- BinTable HDU: X, Y, Z, VX, VY, VZ + (conditionally) R_CLUMP, RHOKAP, TEMP
  !
  !    The three per-clump physical columns are only written when their
  !    spread (max-min) exceeds const_tol * |mean|. Otherwise the value
  !    is fully captured by the corresponding header keyword (CL_RAD,
  !    RHOKAP, TEMP_CL) and the column would be a 4-byte-per-row constant
  !    waste. read_clumps_info() falls back to the header keyword when the
  !    column is absent, so legacy 6-column files and new constant-case
  !    files are read the same way.
  allocate(tmp(ncl))

  tmp = cl_x(1:ncl)
  call io_write_table_column(iofh, 'X',  tmp, status, bitpix)
  tmp = cl_y(1:ncl)
  call io_write_table_column(iofh, 'Y',  tmp, status, bitpix)
  tmp = cl_z(1:ncl)
  call io_write_table_column(iofh, 'Z',  tmp, status, bitpix)
  !--- cl_v* are stored as v/cl_vtherm(icl) internally (since 2026-04-27);
  !    multiply back by cl_vtherm(icl) to write velocities in km/s.
  tmp = cl_vx(1:ncl) * cl_vtherm(1:ncl)
  call io_write_table_column(iofh, 'VX', tmp, status, bitpix)
  tmp = cl_vy(1:ncl) * cl_vtherm(1:ncl)
  call io_write_table_column(iofh, 'VY', tmp, status, bitpix)
  tmp = cl_vz(1:ncl) * cl_vtherm(1:ncl)
  call io_write_table_column(iofh, 'VZ', tmp, status, bitpix)

  !--- per-clump physical columns: write only if non-constant.
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
     call io_write_table_column(iofh, 'R_CLUMP', tmp, status, bitpix)
  end if
  if (write_rhokap) then
     tmp = cl_rhokap(1:ncl)
     call io_write_table_column(iofh, 'RHOKAP', tmp, status, bitpix)
  end if
  if (write_temp) then
     tmp = cl_temperature(1:ncl)
     call io_write_table_column(iofh, 'TEMP', tmp, status, bitpix)
  end if

  deallocate(tmp)

  if (mpar%p_rank == 0) then
     write(*,'(a,3l2)') ' Clumps: per-clump columns saved (R_CLUMP/RHOKAP/TEMP) = ', &
          write_radius, write_rhokap, write_temp
  end if

  !--- Write keywords into the BinTable HDU (current HDU after io_write_table_column)
  call io_put_keyword(iofh, 'N_CLUMPS',  ncl,                'number of clumps (realized)',          status)
  call io_put_keyword(iofh, 'SPHERE_R',  sphere_R,           'outer sphere radius [code units]',      status)
  call io_put_keyword(iofh, 'RMIN',      r_min_clump,        'inner placement radius [code units]',   status)
  call io_put_keyword(iofh, 'CL_RAD',    cl_radius_max,      'clump radius (max) [code units]',       status)
  call io_put_keyword(iofh, 'F_VOL',     f_vol_actual,       'volume filling factor (realized)',      status)
  call io_put_keyword(iofh, 'F_COV',     f_cov_actual,       'covering factor (realized)',            status)
  call io_put_keyword(iofh, 'TAU0',      par%clump_tau0,     'line-center tau (center to surface)',   status)
  call io_put_keyword(iofh, 'SIGMA_V',   par%clump_sigma_v,  'bulk velocity sigma [km/s]',            status)
  call io_put_keyword(iofh, 'TEMP_CL',   cl_temperature_ref, 'clump temperature [K] (reference)',     status)
  call io_put_keyword(iofh, 'RHOKAP',    cl_rhokap_ref,      'opacity/code-unit inside clumps (ref)', status)
  call io_put_keyword(iofh, 'CL_DFREQ',  cl_Dfreq_ref,       'Doppler frequency [Hz] (reference)',    status)
  call io_put_keyword(iofh, 'VTHERM',    cl_vtherm_ref,      'thermal velocity [km/s] (reference)',   status)
  call io_put_keyword(iofh, 'VOIGT_A',   cl_voigt_a_ref,     'Voigt damping parameter (reference)',   status)
  call io_put_keyword(iofh, 'RMAX',      par%rmax,      'outer sphere radius input (par%rmax)',   status)
  call io_put_keyword(iofh, 'IN_FCOV',   par%clump_f_cov,   'covering factor (input)',            status)
  call io_put_keyword(iofh, 'IN_FVOL',   par%clump_f_vol,   'volume filling factor (input)',      status)
  call io_put_keyword(iofh, 'IN_NCL',    par%clump_N_clumps,'N_clumps (input)',                   status)
  call io_put_keyword(iofh, 'IN_NHI',    par%clump_NHI,     'per-clump NHI input [cm^-2] (clump center to surface)',status)
  call io_put_keyword(iofh, 'IN_NH',     par%clump_nH,      'clump nH density input [cm^-3]',    status)
  call io_put_keyword(iofh, 'IN_TEMP',   par%clump_temperature,'clump temperature input [K]',     status)
  call io_put_keyword(iofh, 'DISTUNIT',  par%distance_unit, 'Distance Unit',                     status)
  call io_put_keyword(iofh, 'DIST_CM',   par%distance2cm,   'Distance Unit (cm)',                status)

  call io_close(iofh, status)

  if (mpar%p_rank == 0) write(*,'(2a)') ' Clumps saved to ', trim(fname)
  end subroutine write_clumps_info
  !===========================================================================

  !===========================================================================
  subroutine read_perclump_or_keyword(iofh, colnames, keyname, scratch, ncl, dst, found_name)
  !---------------------------------------------------------------------------
  ! Helper used by read_clumps_info.  Try each comma-separated entry in
  ! `colnames` in order; the first one that exists in the FITS table is
  ! used.  Allows alternative spellings to be accepted (e.g. RHOKAP,
  ! DENSITY, DENS for the per-clump opacity proxy).  If none of the
  ! columns are present, fall back to the value of header keyword
  ! `keyname` and fill `dst` uniformly with it.  This is how
  !   - a legacy file without the R_CLUMP/RHOKAP/TEMP columns, and
  !   - a column-skipping file where the column was dropped because it
  !     was effectively constant
  ! are both handled with a single read path.
  !
  ! Optional output `found_name` returns the name of the column that was
  ! actually used (upper-case, trimmed), or empty if the header-keyword
  ! fallback was taken.  The caller can use it to disambiguate aliases --
  ! e.g. DENSITY/DENS need a units conversion that RHOKAP does not.
  !---------------------------------------------------------------------------
  use iofile_mod
  implicit none
  type(io_file_type), intent(in)             :: iofh
  character(len=*),   intent(in)             :: colnames, keyname
  real(real32),       intent(inout)          :: scratch(:)
  integer(int64),     intent(in)             :: ncl
  real(kind=wp),      intent(out)            :: dst(:)
  character(len=*),   intent(out),  optional :: found_name

  integer            :: status, colnum, ierr, ks, ke, lstr
  real(kind=wp)      :: keyval
  character(len=80)  :: cmt
  character(len=64)  :: candidate
  character(len=64)  :: tried
  logical            :: found

  found = .false.
  tried = ''
  if (present(found_name)) found_name = ''
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
     call io_get_column_number(iofh, trim(candidate), colnum, status)
     if (status == 0) then
        call io_read_table_column(iofh, colnum, scratch, status)
        dst(1:ncl) = real(scratch(1:ncl), wp)
        if (present(found_name)) found_name = trim(candidate)
        found = .true.
     end if
  end do

  if (.not. found) then
     status = 0
     call io_get_keyword(iofh, trim(keyname), keyval, status, cmt)
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
  subroutine read_clumps_info(fname, R_sphere)
  !---------------------------------------------------------------------------
  ! Read a clump population from a FITS file produced by write_clumps_info()
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
  ! so write_clumps_info and other diagnostic branches can detect this case.
  !---------------------------------------------------------------------------
  use define
  use line_mod
  use iofile_mod
  use voigt_mod, only: voigt
  use mpi
  implicit none
  character(len=*), intent(in) :: fname
  real(kind=wp),    intent(in) :: R_sphere

  type(io_file_type) :: iofh
  integer        :: status, ierr, colnum
  integer(int64) :: ncl, i
  real(kind=real32), allocatable :: tmp(:)
  character(len=80)  :: cmt
  character(len=64)  :: rhokap_colname
  character(len=128) :: distunit_file
  real(kind=wp)      :: distcm_file
  logical            :: rhokap_is_density, distance_unit_ok, user_set_distunit
  character(len=:), allocatable :: resolved_fname

  status   = 0
  sphere_R = R_sphere
  base_radius_in  = par%clump_radius
  profiles_active = .false.
  clumps_from_file = .true.
  r_min_clump = max(0.0_wp, par%rmin)

  ! If `fname` has no extension, try `.h5`, `.fits.gz`, `.fits`, `.dat`
  ! in that order via io_resolve_filename.
  resolved_fname = io_resolve_filename(fname)

  if (mpar%p_rank == 0) then
     write(*,'(2a)') ' Clumps: reading from ', resolved_fname
     call io_open_old(iofh, resolved_fname, status)
     if (status /= 0) then
        write(*,'(2a)') 'ERROR: cannot open clump_input_file: ', resolved_fname
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
     !--- HDU 1 is the empty primary; the binary table is HDU 2.
     call io_move_to_next_section(iofh, status)
     if (status /= 0) then
        write(*,'(a)') 'ERROR: clump_input_file: cannot reach binary-table HDU 2'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
     call io_get_keyword(iofh, 'N_CLUMPS', ncl, status, cmt)
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
  if (par%DGR > 0.0_wp) call create_shared_mem(cl_rhokapD, [int(N_clumps)])
  call create_shared_mem(cl_voigt_a,     [int(N_clumps)])
  call create_shared_mem(cl_Dfreq,       [int(N_clumps)])
  call create_shared_mem(cl_vtherm,      [int(N_clumps)])
  call create_shared_mem(cl_temperature, [int(N_clumps)])

  !--- Read columns on p_rank=0 (which is also h_rank=0 of node 0); the
  !    other h_rank=0 ranks receive copies via MPI_BCAST below.
  if (mpar%p_rank == 0) then
     allocate(tmp(ncl))

     call io_get_column_number(iofh, 'X', colnum, status)
     call io_read_table_column(iofh, colnum, tmp, status)
     cl_x = real(tmp, dp)
     call io_get_column_number(iofh, 'Y', colnum, status)
     call io_read_table_column(iofh, colnum, tmp, status)
     cl_y = real(tmp, dp)
     call io_get_column_number(iofh, 'Z', colnum, status)
     call io_read_table_column(iofh, colnum, tmp, status)
     cl_z = real(tmp, dp)
     call io_get_column_number(iofh, 'VX', colnum, status)
     call io_read_table_column(iofh, colnum, tmp, status)
     cl_vx = real(tmp, dp)        ! km/s; converted to dimensionless below
     call io_get_column_number(iofh, 'VY', colnum, status)
     call io_read_table_column(iofh, colnum, tmp, status)
     cl_vy = real(tmp, dp)
     call io_get_column_number(iofh, 'VZ', colnum, status)
     call io_read_table_column(iofh, colnum, tmp, status)
     cl_vz = real(tmp, dp)

     !--- R_CLUMP / RHOKAP / TEMP are optional. write_clumps_info omits
     !    each one when its spread is below const_tol (= 1e-3 of mean).
     !    Legacy FITS files do not carry these columns at all.
     !    In both cases we fall back to the corresponding header keyword.
     call read_perclump_or_keyword(iofh, 'R_CLUMP',                'CL_RAD',  tmp, ncl, cl_radius)
     call read_perclump_or_keyword(iofh, 'RHOKAP,DENSITY,DENS',    'RHOKAP',  tmp, ncl, cl_rhokap, &
                                   found_name=rhokap_colname)
     call read_perclump_or_keyword(iofh, 'TEMP',                   'TEMP_CL', tmp, ncl, cl_temperature)

     !--- Optional DISTUNIT / DIST_CM keywords (added by write_clumps_info
     !    from 2026-05-12 onward).  Cached locally; the adoption decision is
     !    made below, after closing the file.
     status        = 0
     distunit_file = ''
     call io_get_keyword(iofh, 'DISTUNIT', distunit_file, status, cmt)
     if (status /= 0) distunit_file = ''
     status        = 0
     distcm_file   = 0.0_wp
     call io_get_keyword(iofh, 'DIST_CM', distcm_file, status, cmt)
     if (status /= 0) distcm_file = 0.0_wp

     deallocate(tmp)
     call io_close(iofh, status)

     !--- Adopt DISTUNIT / DIST_CM from the file when the user did NOT
     !    supply par%distance_unit/par%distance2cm in the input.  The
     !    "user-supplied" criterion is
     !        par%distance_unit /= '' .AND. par%distance2cm > 1.0
     !    so the file values are picked up when EITHER condition fails
     !    (default state: distance_unit = '', distance2cm = 1.0).
     user_set_distunit = (len_trim(par%distance_unit) > 0 .and. par%distance2cm > 1.0_wp)
     if (.not. user_set_distunit) then
        if (len_trim(distunit_file) > 0) then
           par%distance_unit = trim(distunit_file)
           write(*,'(3a)') ' Clumps: adopting DISTUNIT = ''', &
                trim(par%distance_unit), ''' from clump_input_file'
        end if
        if (distcm_file > 0.0_wp) then
           par%distance2cm = distcm_file
           write(*,'(a,es12.5,a)') ' Clumps: adopting DIST_CM = ', &
                par%distance2cm, ' cm/code-unit from clump_input_file'
        end if
     end if

     !--- Derive per-clump Doppler frequency, thermal velocity, Voigt damping
     !    from the per-clump temperature so that the line-data choice (par%line_id)
     !    used in the LaRT run can differ from the one that generated the file.
     do i = 1_int64, ncl
        cl_vtherm(i)  = vtherm_total(cl_temperature(i))
        cl_Dfreq(i)   = cl_vtherm(i) / (line%wavelength0 * um2km)
        cl_voigt_a(i) = (line%damping / fourpi) / cl_Dfreq(i)
        cl_radius2(i) = cl_radius(i) * cl_radius(i)
        !--- VX/VY/VZ in the FITS file are stored in km/s (write_clumps_info
        !    multiplies by cl_vtherm). LaRT internally uses v / cl_vtherm(icl).
        cl_vx(i) = cl_vx(i) / cl_vtherm(i)
        cl_vy(i) = cl_vy(i) / cl_vtherm(i)
        cl_vz(i) = cl_vz(i) / cl_vtherm(i)
     end do

     !--- Column-name interpretation for the opacity field.  If the file
     !    used DENSITY or DENS rather than RHOKAP, the values are number
     !    densities [cm^-3] and need to be converted to line-center opacity
     !    per code unit via
     !        rhokap_i = n_HI_i * line%cross0 / cl_Dfreq(i) * par%distance2cm.
     !    The conversion only needs a meaningful distance2cm, so the check
     !    is OR (looser than the adoption check above) -- if EITHER
     !    distance_unit was supplied (so setup.f90 set distance2cm) OR a
     !    raw distance2cm was supplied / adopted from DIST_CM, the
     !    conversion has the data it needs.  Otherwise it is undefined and
     !    we emit a warning, treating the values as RHOKAP (the legacy
     !    interpretation).
     rhokap_is_density = (trim(rhokap_colname) == 'DENSITY' .or. &
                          trim(rhokap_colname) == 'DENS')
     distance_unit_ok  = (len_trim(par%distance_unit) > 0 .or. par%distance2cm > 1.0_wp)
     if (rhokap_is_density) then
        if (distance_unit_ok) then
           do i = 1_int64, ncl
              cl_rhokap(i) = cl_rhokap(i) * line%cross0 / cl_Dfreq(i) * par%distance2cm
           end do
           write(*,'(3a,es12.4,a)') ' Clumps: ', trim(rhokap_colname), &
                ' column converted to opacity using par%distance2cm = ', par%distance2cm, &
                ' cm/code-unit'
        else
           write(*,'(3a)')  ' Clumps: WARNING -- column ', trim(rhokap_colname), &
                ' found but par%distance_unit/distance2cm is not set;'
           write(*,'(a)')   '                    treating the values as RHOKAP'// &
                ' (opacity per code unit) without conversion.'
        end if
     end if
  end if

  !--- par%distance_unit / par%distance2cm may have been adopted from the
  !    clump file on p_rank=0 above; broadcast to keep all ranks consistent
  !    for downstream output (write_output_rect writes DISTUNIT/DIST_CM into
  !    the result FITS, and the peel-off/observer paths use distance2cm to
  !    convert code-unit areas to cm^2).
  call MPI_BCAST(par%distance_unit, len(par%distance_unit), MPI_CHARACTER,        0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(par%distance2cm,   1,                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

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

  !--- Optional rescaling: when par%taumax or par%N_HImax is supplied along
  !    with a pre-built clump file, multiply every cl_rhokap(i) by a single
  !    factor so the realized radial-sightline tau / NHI matches the target.
  !    The realized value uses the small-angle (V/d^2) approximation, capped
  !    at d^2 >= cl_radius^2 to avoid divergence when a clump straddles the
  !    origin (e.g. when par%rmin = 0 and a clump center lies near the box
  !    center).  Both cl_rhokap and cl_rhokap_ref are updated in lockstep
  !    so write_clumps_info and grid_create_clump see the rescaled values.
  call rescale_loaded_clumps_to_target()

  !--- Per-clump dust opacity from the finalized (rescaled) cl_rhokap, using
  !    the same Cartesian convention as the generate path:
  !      cl_rhokapD = cl_rhokap * cext_dust * DGR * cl_Dfreq / cross0.
  if (par%DGR > 0.0_wp .and. associated(cl_rhokapD)) then
     if (mpar%h_rank == 0) then
        do i = 1_int64, N_clumps
           cl_rhokapD(i) = cl_rhokap(i) * par%cext_dust * par%DGR * cl_Dfreq(i) / line%cross0
        end do
     end if
     call MPI_BARRIER(mpar%hostcomm, ierr)
  end if

  if (mpar%p_rank == 0) then
     write(*,'(a,i14)')    ' Clumps: N_clumps  = ', N_clumps
     write(*,'(a,es12.4)') ' Clumps: cl_rhokap = ', cl_rhokap_ref
     write(*,'(a,f12.5)')  ' Clumps: voigt_a   = ', cl_voigt_a_ref
     write(*,'(a,es12.4)') ' Clumps: cl_Dfreq  = ', cl_Dfreq_ref
     write(*,'(a,2es12.4)')' Clumps: r_cl min/max  = ', minval(cl_radius(1:ncl)), cl_radius_max
  end if

  call build_clump_csr()
  call MPI_BARRIER(mpar%hostcomm, ierr)

  end subroutine read_clumps_info
  !===========================================================================

  !===========================================================================
  subroutine compute_clump_scalars(tauhomo_out, taumax_out, N_gashomo_out, N_gasmax_out)
  !---------------------------------------------------------------------------
  ! Compute the four system-level scalars from the per-clump arrays
  ! cl_rhokap, cl_voigt_a, cl_Dfreq, cl_radius, cl_x/y/z (all assumed
  ! populated; works for both loaded and generated clumps, uniform or
  ! profile mode).
  !
  ! sphere_R and r_min_clump are module variables set in init_clumps;
  ! line%cross0 is the line-center cross section [cm^2 * Hz].
  !
  !   tauhomo  = Sum_i rhokap_i * voigt0_i * r_i^3 / (R^2 + R*r0 + r0^2)
  !              (uniform-smear of clump opacity over the shell, traversed
  !               radially)
  !
  !   taumax   = Sum_i rhokap_i * voigt0_i * r_i^3 / (3 * max(d_i^2, r_i^2))
  !              (small-angle expected line-center tau along a radial
  !               sightline from the origin; cap on d^2 avoids the singularity
  !               for a clump that straddles the origin)
  !
  !   N_gashomo / N_gasmax = same forms with cl_Dfreq(i) / line%cross0 in
  !              place of voigt0_i.
  !
  ! In the uniform-isotropic limit (all r_i = r_cl, all rhokap_i = rhokap_ref,
  ! and clumps distributed uniformly inside the shell) both _homo and _max
  ! reduce to the same closed-form value
  !
  !       N * r_cl^3 * rhokap_ref * voigt0 / (R^2 + R*r0 + r0^2)
  !     = (4/3) * f_cov * tau_per_clump
  !
  ! so this routine is backward compatible with the original closed-form
  ! formula used for procedurally generated uniform clumps.  Computed on
  ! p_rank == 0, then broadcast.
  !---------------------------------------------------------------------------
  use line_mod
  implicit none
  real(kind=wp), intent(out) :: tauhomo_out, taumax_out, N_gashomo_out, N_gasmax_out

  integer        :: ierr
  integer(int64) :: i
  real(kind=wp)  :: voigt0_i, di2, shell_R2_factor, w_homo, w_max
  real(kind=wp)  :: vals(4)

  shell_R2_factor = sphere_R**2 + sphere_R*r_min_clump + r_min_clump**2
  if (shell_R2_factor <= 0.0_wp) shell_R2_factor = sphere_R**2

  vals(:) = 0.0_wp
  if (mpar%p_rank == 0 .and. N_clumps > 0_int64) then
     do i = 1_int64, N_clumps
        voigt0_i = voigt(0.0_wp, cl_voigt_a(i))
        di2      = max(cl_x(i)**2 + cl_y(i)**2 + cl_z(i)**2, cl_radius(i)**2)
        w_homo   = cl_radius(i)**3 / shell_R2_factor
        w_max    = cl_radius(i)**3 / (3.0_wp * di2)
        vals(1)  = vals(1) + cl_rhokap(i) * voigt0_i    * w_homo
        vals(2)  = vals(2) + cl_rhokap(i) * voigt0_i    * w_max
        vals(3)  = vals(3) + cl_rhokap(i) * cl_Dfreq(i) * w_homo / line%cross0
        vals(4)  = vals(4) + cl_rhokap(i) * cl_Dfreq(i) * w_max  / line%cross0
     end do
  end if
  call MPI_BCAST(vals, 4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  tauhomo_out   = vals(1)
  taumax_out    = vals(2)
  N_gashomo_out = vals(3)
  N_gasmax_out  = vals(4)

  end subroutine compute_clump_scalars
  !===========================================================================

  !===========================================================================
  subroutine rescale_loaded_clumps_to_target()
  !---------------------------------------------------------------------------
  ! Helper for read_clumps_info.  When the user supplies a system-level
  ! target via par%taumax / par%tauhomo / par%N_gasmax / par%N_gashomo
  ! (priority order, highest first), rescales every loaded cl_rhokap(i) by
  ! a single multiplicative factor so the realized scalar (computed by
  ! compute_clump_scalars from the loaded distribution) matches the target.
  !
  ! Note: par%N_HImax / par%N_HIhomo are aliased to par%N_gasmax /
  ! par%N_gashomo in setup.f90 before this is called, so checking the
  ! N_gas* names alone is sufficient.
  !
  ! If none of the four targets are set, returns without modifying
  ! cl_rhokap -- the per-clump opacities loaded from the file are taken
  ! as-is and the four system-level scalars are derived from them in
  ! grid_mod_clump.
  !---------------------------------------------------------------------------
  use line_mod
  implicit none
  integer           :: ierr
  integer(int64)    :: i
  real(kind=wp)     :: realized, target_val, alpha
  real(kind=wp)     :: tauhomo_r, taumax_r, N_gashomo_r, N_gasmax_r
  character(len=12) :: which_target

  if (par%taumax    <= 0.0_wp .and. par%tauhomo   <= 0.0_wp .and. &
      par%N_gasmax  <= 0.0_wp .and. par%N_gashomo <= 0.0_wp) return

  call compute_clump_scalars(tauhomo_r, taumax_r, N_gashomo_r, N_gasmax_r)

  if (par%taumax > 0.0_wp) then
     realized     = taumax_r;    target_val = par%taumax;    which_target = 'taumax'
  else if (par%tauhomo > 0.0_wp) then
     realized     = tauhomo_r;   target_val = par%tauhomo;   which_target = 'tauhomo'
  else if (par%N_gasmax > 0.0_wp) then
     realized     = N_gasmax_r;  target_val = par%N_gasmax;  which_target = 'N_gasmax'
  else
     realized     = N_gashomo_r; target_val = par%N_gashomo; which_target = 'N_gashomo'
  end if

  if (realized <= 0.0_wp) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        ' Clumps: WARNING -- realized tau/NHI from loaded clumps is zero; '// &
        'cannot rescale. Leaving cl_rhokap unchanged.'
     return
  end if
  alpha = target_val / realized

  if (mpar%h_rank == 0) then
     do i = 1_int64, N_clumps
        cl_rhokap(i) = cl_rhokap(i) * alpha
     end do
  end if
  call MPI_BARRIER(mpar%hostcomm, ierr)

  cl_rhokap_ref  = cl_rhokap_ref * alpha
  base_rhokap_in = cl_rhokap_ref

  if (mpar%p_rank == 0) then
     write(*,'(3a,es12.4,a,es12.4)') &
        ' Clumps: rescaled to par%', trim(which_target), ' = ', target_val, &
        ', alpha = ', alpha
  end if

  end subroutine rescale_loaded_clumps_to_target
  !===========================================================================

end module clump_mod
