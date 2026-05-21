module raytrace_clump_mod
!---------------------------------------------------------------------------
! Raytrace routines for the clump medium.
!
! The medium consists of N spherical clumps (radii cl_radius(:)) placed inside
! a sphere of radius sphere_R. Outside the clumps the medium is vacuum.
! Each clump has its own opacity cl_rhokap(icl), Voigt parameter cl_voigt_a(icl),
! Doppler frequency cl_Dfreq(icl), and thermal velocity cl_vtherm(icl).
! With uniform clump properties every entry equals the corresponding _ref
! scalar; when a radial profile is active, init_clumps sets the entries
! as functions of clump-center radius.
!
! Frequency convention (same as Cartesian):
!   photon%xfreq is in units of cl_Dfreq_ref (= grid%Dfreq_ref).
!   cl_vx/y/z are stored as v / cl_vtherm(icl) (dimensionless), set in init_clumps.
!   At clump entry: xfreq -= v_los_clump   (v already in vtherm units)
!   At clump exit:  xfreq += v_los_clump
!
! photon%icell_clump: 0 = in vacuum; > 0 = index of current clump (1-based).
! photon%icell/jcell/kcell: maintained using floor() map from position.
!---------------------------------------------------------------------------
  use define,    only: wp
  use clump_mod
  use voigt_mod
  use line_mod
  implicit none
  private

  public :: raytrace_to_tau_clump
  public :: raytrace_to_edge_clump
  public :: raytrace_to_edge_tau_gas_clump
  public :: raytrace_to_edge_column_clump
  public :: raytrace_to_dist_clump
  public :: raytrace_to_dist_tau_gas_clump
  public :: raytrace_to_dist_column_clump
  public :: tau_huge_clump
  ! Capped variants: same algorithm as raytrace_to_edge_clump /
  ! raytrace_to_dist_clump but exit early once tau >= tau_max.  Use these
  ! from photon-transport callers (peeling-off) where exp(-tau) underflows
  ! above ~745; the uncapped originals stay around because their signature
  ! is locked by the raytrace_to_edge / raytrace_to_dist procedure-pointer
  ! interface in define.f90.
  public :: raytrace_to_edge_clump_capped
  public :: raytrace_to_dist_clump_capped
  ! Overlap-aware variants: used when has_overlap=.true. (file-loaded clumps)
  public :: raytrace_to_tau_clump_overlap
  public :: raytrace_to_edge_clump_overlap
  public :: raytrace_to_edge_clump_overlap_capped
  public :: raytrace_to_edge_tau_gas_clump_overlap
  public :: raytrace_to_edge_column_clump_overlap
  public :: raytrace_to_dist_clump_overlap
  public :: raytrace_to_dist_clump_overlap_capped
  public :: raytrace_to_dist_tau_gas_clump_overlap
  public :: raytrace_to_dist_column_clump_overlap

  ! Default early-exit threshold for photon-transport callers (peeling, ...):
  ! once tau exceeds this, exp(-tau) is zero in double precision so further
  ! accumulation only burns CPU.  Sight-line callers DO NOT pass tau_max --
  ! they need the full integrated tau, however large.
  real(kind=wp), parameter :: tau_huge_clump = 745.2_wp

contains

  !===========================================================================
  ! Internal helper: update photon%icell/jcell/kcell from position.
  !===========================================================================
  pure subroutine update_cell_idx(photon, grid)
  use define
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(in)    :: grid
  photon%icell = max(1, min(grid%nx, floor((photon%x - grid%xmin)/grid%dx) + 1))
  photon%jcell = max(1, min(grid%ny, floor((photon%y - grid%ymin)/grid%dy) + 1))
  photon%kcell = max(1, min(grid%nz, floor((photon%z - grid%zmin)/grid%dz) + 1))
  end subroutine update_cell_idx
  !===========================================================================

  !===========================================================================
  ! Internal helper: accumulate Jout contribution at current photon state.
  ! (Used by raytrace_to_edge_clump when photon exits sphere.)
  !===========================================================================

  !===========================================================================
  subroutine raytrace_to_tau_clump(photon, grid, tau_in)
  !---------------------------------------------------------------------------
  ! Advance photon through the clump medium until optical depth tau_in is
  ! accumulated or photon leaves the sphere.
  ! Updates: photon%x/y/z, photon%xfreq, photon%icell_clump, photon%inside.
  !---------------------------------------------------------------------------
  use define
  implicit none
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid
  real(kind=wp),     intent(in)    :: tau_in

  real(kind=wp)  :: kx, ky, kz, tau_rem, kap, ds, t_seg, t_sp
  real(kind=wp)  :: te, tx2, u_los
  integer(int64) :: icl, icl_found, last_icl
  logical        :: found
  integer        :: ix_jout

  kx = photon%kx;  ky = photon%ky;  kz = photon%kz
  tau_rem  = tau_in
  last_icl = 0_int64   ! index of most-recently exited clump (to avoid re-entry)

  do while (photon%inside)

     if (photon%icell_clump > 0) then
        !--- photon is inside clump icl
        icl   = int(photon%icell_clump, int64)
        t_seg = clump_exit_dist(photon%x, photon%y, photon%z, kx, ky, kz, icl)

        kap = cl_rhokap(icl) * voigt_clump(photon%xfreq, icl)

        if (tau_rem <= kap * t_seg) then
           !--- scatter inside this clump
           ds = tau_rem / max(kap, tiny(1.0_wp))
           photon%x = photon%x + ds * kx
           photon%y = photon%y + ds * ky
           photon%z = photon%z + ds * kz
           call update_cell_idx(photon, grid)
           return
        end if

        !--- traverse through this clump and exit
        tau_rem = tau_rem - kap * t_seg
        photon%x = photon%x + t_seg * kx
        photon%y = photon%y + t_seg * ky
        photon%z = photon%z + t_seg * kz

        !--- shift xfreq back to global frame at exit
        u_los = ulos_clump(icl, kx, ky, kz)
        photon%xfreq = photon%xfreq + u_los

        last_icl           = icl
        photon%icell_clump = 0

        !--- check sphere exit
        if (photon%x**2 + photon%y**2 + photon%z**2 >= sphere_R**2) then
           photon%inside = .false.
           call update_cell_idx(photon, grid)
           ix_jout = floor((photon%xfreq - grid%xfreq_min)/grid%dxfreq) + 1
           if (ix_jout >= 1 .and. ix_jout <= grid%nxfreq) then
              !$OMP ATOMIC UPDATE
              grid%Jout(ix_jout) = grid%Jout(ix_jout) + photon%wgt
              if (par%save_Jmu) call add_to_Jmu_clump(photon, grid, ix_jout)
           end if
           return
        end if

     else
        !--- photon in vacuum: find next clump or sphere exit
        t_sp = sphere_exit_dist(photon%x, photon%y, photon%z, kx, ky, kz)
        if (t_sp <= 0.0_wp) then
           photon%inside = .false.
           ix_jout = floor((photon%xfreq - grid%xfreq_min)/grid%dxfreq) + 1
           if (ix_jout >= 1 .and. ix_jout <= grid%nxfreq) then
              !$OMP ATOMIC UPDATE
              grid%Jout(ix_jout) = grid%Jout(ix_jout) + photon%wgt
              if (par%save_Jmu) call add_to_Jmu_clump(photon, grid, ix_jout)
           end if
           return
        end if

        call find_next_clump(photon%x, photon%y, photon%z, kx, ky, kz, &
             last_icl, t_sp, te, tx2, icl_found, found)

        if (found) then
           !--- advance to clump entry (vacuum, no tau)
           te = max(0.0_wp, te)
           photon%x = photon%x + te * kx
           photon%y = photon%y + te * ky
           photon%z = photon%z + te * kz

           !--- shift xfreq to clump frame at entry
           u_los = ulos_clump(icl_found, kx, ky, kz)
           photon%xfreq = photon%xfreq - u_los

           last_icl           = 0_int64   ! entering a new clump; no skip needed
           photon%icell_clump = int(icl_found)
           call update_cell_idx(photon, grid)
        else
           !--- no clump on this ray: skip to sphere boundary
           photon%x = photon%x + t_sp * kx
           photon%y = photon%y + t_sp * ky
           photon%z = photon%z + t_sp * kz
           photon%inside = .false.
           call update_cell_idx(photon, grid)
           ix_jout = floor((photon%xfreq - grid%xfreq_min)/grid%dxfreq) + 1
           if (ix_jout >= 1 .and. ix_jout <= grid%nxfreq) then
              !$OMP ATOMIC UPDATE
              grid%Jout(ix_jout) = grid%Jout(ix_jout) + photon%wgt
              if (par%save_Jmu) call add_to_Jmu_clump(photon, grid, ix_jout)
           end if
           return
        end if

     end if

  end do

  end subroutine raytrace_to_tau_clump
  !===========================================================================

  !===========================================================================
  subroutine raytrace_to_edge_clump(photon0, grid, tau)
  !---------------------------------------------------------------------------
  ! Accumulate total optical depth from photon0 to the sphere edge along k.
  ! photon0 is read-only. tau = integral of rhokap*H(x,a) ds through clumps.
  !---------------------------------------------------------------------------
  use define
  implicit none
  type(photon_type), intent(in)  :: photon0
  type(grid_type),   intent(in)  :: grid
  real(kind=wp),     intent(out) :: tau

  real(kind=wp)  :: xp, yp, zp, kx, ky, kz
  real(kind=wp)  :: xfreq_loc
  real(kind=wp)  :: t_sp, t_seg, te, tx2, u_los, kap
  integer(int64) :: icl, icl_cur, icl_found
  logical        :: found

  tau      = 0.0_wp
  xp       = photon0%x;  yp = photon0%y;  zp = photon0%z
  kx       = photon0%kx; ky = photon0%ky; kz = photon0%kz
  xfreq_loc = photon0%xfreq
  icl_cur   = int(photon0%icell_clump, int64)

  !--- first, if photon starts inside a clump, traverse it
  if (icl_cur > 0) then
     t_seg = clump_exit_dist(xp, yp, zp, kx, ky, kz, icl_cur)
     kap   = cl_rhokap(icl_cur) * voigt_clump(xfreq_loc, icl_cur)
     tau   = tau + kap * t_seg
     xp    = xp + t_seg * kx
     yp    = yp + t_seg * ky
     zp    = zp + t_seg * kz
     u_los = ulos_clump(icl_cur, kx, ky, kz)
     xfreq_loc = xfreq_loc + u_los
     ! keep icl_cur as skip index for first find_next_clump call
     if (xp**2 + yp**2 + zp**2 >= sphere_R**2) return
  end if

  !--- DDA through remaining clumps to sphere boundary
  do
     t_sp = sphere_exit_dist(xp, yp, zp, kx, ky, kz)
     if (t_sp <= 0.0_wp) exit

     call find_next_clump(xp, yp, zp, kx, ky, kz, icl_cur, t_sp, &
                           te, tx2, icl_found, found)
     if (.not. found) exit

     !--- advance to clump entry, accumulate tau through clump
     te = max(0.0_wp, te)
     xp = xp + te * kx;  yp = yp + te * ky;  zp = zp + te * kz

     !--- xfreq shift at entry
     u_los = ulos_clump(icl_found, kx, ky, kz)
     xfreq_loc = xfreq_loc - u_los

     !--- tau through this clump
     t_seg = clump_exit_dist(xp, yp, zp, kx, ky, kz, icl_found)
     kap   = cl_rhokap(icl_found) * voigt_clump(xfreq_loc, icl_found)
     tau   = tau + kap * t_seg
     xp    = xp + t_seg * kx;  yp = yp + t_seg * ky;  zp = zp + t_seg * kz
     xfreq_loc = xfreq_loc + u_los  ! restore global frame

     icl_cur = icl_found  ! skip this clump in next search (just exited)
     if (xp**2 + yp**2 + zp**2 >= sphere_R**2) exit
  end do

  end subroutine raytrace_to_edge_clump
  !===========================================================================

  !===========================================================================
  subroutine raytrace_to_edge_tau_gas_clump(photon0, grid, tau_gas)
  use define
  implicit none
  type(photon_type), intent(in)  :: photon0
  type(grid_type),   intent(in)  :: grid
  real(kind=wp),     intent(out) :: tau_gas
  ! For clump medium, all opacity is gas (no dust unless DGR>0 handled elsewhere)
  call raytrace_to_edge_clump(photon0, grid, tau_gas)
  end subroutine raytrace_to_edge_tau_gas_clump
  !===========================================================================

  !===========================================================================
  subroutine raytrace_to_edge_column_clump(photon0, grid, N_gas, tau_dust)
  use define
  implicit none
  type(photon_type), intent(in)  :: photon0
  type(grid_type),   intent(in)  :: grid
  real(kind=wp),     intent(out) :: N_gas, tau_dust

  real(kind=wp)  :: xp, yp, zp, kx, ky, kz
  real(kind=wp)  :: t_sp, t_seg, te, tx2, u_los
  integer(int64) :: icl, icl_cur, icl_found
  logical        :: found

  N_gas    = 0.0_wp
  tau_dust = 0.0_wp
  xp = photon0%x;  yp = photon0%y;  zp = photon0%z
  kx = photon0%kx; ky = photon0%ky; kz = photon0%kz
  icl_cur = int(photon0%icell_clump, int64)

  if (icl_cur > 0) then
     t_seg = clump_exit_dist(xp, yp, zp, kx, ky, kz, icl_cur)
     N_gas = N_gas + (cl_rhokap(icl_cur) / line%cross0) * cl_Dfreq(icl_cur) * t_seg
     xp = xp + t_seg*kx;  yp = yp + t_seg*ky;  zp = zp + t_seg*kz
     ! keep icl_cur as skip for first find_next_clump
     if (xp**2 + yp**2 + zp**2 >= sphere_R**2) return
  end if

  do
     t_sp = sphere_exit_dist(xp, yp, zp, kx, ky, kz)
     if (t_sp <= 0.0_wp) exit
     call find_next_clump(xp, yp, zp, kx, ky, kz, icl_cur, t_sp, &
                           te, tx2, icl_found, found)
     if (.not. found) exit
     te = max(0.0_wp, te)
     xp = xp + te*kx;  yp = yp + te*ky;  zp = zp + te*kz
     t_seg = clump_exit_dist(xp, yp, zp, kx, ky, kz, icl_found)
     N_gas = N_gas + (cl_rhokap(icl_found) / line%cross0) * cl_Dfreq(icl_found) * t_seg
     xp = xp + t_seg*kx;  yp = yp + t_seg*ky;  zp = zp + t_seg*kz
     icl_cur = icl_found
     if (xp**2 + yp**2 + zp**2 >= sphere_R**2) exit
  end do

  end subroutine raytrace_to_edge_column_clump
  !===========================================================================

  !===========================================================================
  subroutine raytrace_to_dist_clump(photon0, grid, dist_in, tau)
  !---------------------------------------------------------------------------
  ! Accumulate optical depth to distance dist_in along the photon ray.
  ! Used by peeling-off routines.
  !---------------------------------------------------------------------------
  use define
  implicit none
  type(photon_type), intent(in)  :: photon0
  type(grid_type),   intent(in)  :: grid
  real(kind=wp),     intent(in)  :: dist_in
  real(kind=wp),     intent(out) :: tau

  real(kind=wp)  :: xp, yp, zp, kx, ky, kz, xfreq_loc
  real(kind=wp)  :: t_rem, t_seg, te, tx2, u_los, kap
  integer(int64) :: icl_cur, icl_found
  logical        :: found

  tau = 0.0_wp
  xp  = photon0%x;  yp = photon0%y;  zp = photon0%z
  kx  = photon0%kx; ky = photon0%ky; kz = photon0%kz
  xfreq_loc = photon0%xfreq
  icl_cur   = int(photon0%icell_clump, int64)
  t_rem     = dist_in

  if (icl_cur > 0) then
     t_seg = min(clump_exit_dist(xp, yp, zp, kx, ky, kz, icl_cur), t_rem)
     kap   = cl_rhokap(icl_cur) * voigt_clump(xfreq_loc, icl_cur)
     tau   = tau + kap * t_seg
     t_rem = t_rem - t_seg
     if (t_rem <= 0.0_wp) return
     xp = xp + t_seg*kx;  yp = yp + t_seg*ky;  zp = zp + t_seg*kz
     u_los = ulos_clump(icl_cur, kx, ky, kz)
     xfreq_loc = xfreq_loc + u_los
     ! keep icl_cur as skip for first find_next_clump
  end if

  do while (t_rem > 0.0_wp)
     call find_next_clump(xp, yp, zp, kx, ky, kz, icl_cur, t_rem, &
                           te, tx2, icl_found, found)
     if (.not. found) exit

     te = max(0.0_wp, te)
     xp = xp + te*kx;  yp = yp + te*ky;  zp = zp + te*kz
     t_rem = t_rem - te

     u_los = ulos_clump(icl_found, kx, ky, kz)
     xfreq_loc = xfreq_loc - u_los

     t_seg = min(clump_exit_dist(xp, yp, zp, kx, ky, kz, icl_found), t_rem)
     kap   = cl_rhokap(icl_found) * voigt_clump(xfreq_loc, icl_found)
     tau   = tau + kap * t_seg
     t_rem = t_rem - t_seg
     if (t_rem <= 0.0_wp) return

     xp = xp + t_seg*kx;  yp = yp + t_seg*ky;  zp = zp + t_seg*kz
     xfreq_loc = xfreq_loc + u_los
     icl_cur = icl_found
  end do

  end subroutine raytrace_to_dist_clump
  !===========================================================================

  !===========================================================================
  subroutine raytrace_to_dist_tau_gas_clump(photon0, grid, dist_in, tau_gas)
  use define
  implicit none
  type(photon_type), intent(in)  :: photon0
  type(grid_type),   intent(in)  :: grid
  real(kind=wp),     intent(in)  :: dist_in
  real(kind=wp),     intent(out) :: tau_gas
  call raytrace_to_dist_clump(photon0, grid, dist_in, tau_gas)
  end subroutine raytrace_to_dist_tau_gas_clump
  !===========================================================================

  !===========================================================================
  subroutine raytrace_to_dist_column_clump(photon0, grid, dist_in, N_gas, tau_dust)
  use define
  implicit none
  type(photon_type), intent(in)  :: photon0
  type(grid_type),   intent(in)  :: grid
  real(kind=wp),     intent(in)  :: dist_in
  real(kind=wp),     intent(out) :: N_gas, tau_dust

  real(kind=wp)  :: xp, yp, zp, kx, ky, kz
  real(kind=wp)  :: t_rem, t_seg, te, tx2
  integer(int64) :: icl_cur, icl_found
  logical        :: found

  N_gas    = 0.0_wp
  tau_dust = 0.0_wp
  xp = photon0%x;  yp = photon0%y;  zp = photon0%z
  kx = photon0%kx; ky = photon0%ky; kz = photon0%kz
  icl_cur = int(photon0%icell_clump, int64)
  t_rem   = dist_in

  if (icl_cur > 0) then
     t_seg = min(clump_exit_dist(xp, yp, zp, kx, ky, kz, icl_cur), t_rem)
     N_gas = N_gas + (cl_rhokap(icl_cur) / line%cross0) * cl_Dfreq(icl_cur) * t_seg
     t_rem = t_rem - t_seg
     if (t_rem <= 0.0_wp) return
     xp = xp + t_seg*kx;  yp = yp + t_seg*ky;  zp = zp + t_seg*kz
     ! keep icl_cur as skip for first find_next_clump
  end if

  do while (t_rem > 0.0_wp)
     call find_next_clump(xp, yp, zp, kx, ky, kz, icl_cur, t_rem, &
                           te, tx2, icl_found, found)
     if (.not. found) exit
     te = max(0.0_wp, te)
     xp = xp + te*kx;  yp = yp + te*ky;  zp = zp + te*kz
     t_rem = t_rem - te
     t_seg = min(clump_exit_dist(xp, yp, zp, kx, ky, kz, icl_found), t_rem)
     N_gas = N_gas + (cl_rhokap(icl_found) / line%cross0) * cl_Dfreq(icl_found) * t_seg
     t_rem = t_rem - t_seg
     if (t_rem <= 0.0_wp) return
     xp = xp + t_seg*kx;  yp = yp + t_seg*ky;  zp = zp + t_seg*kz
     icl_cur = icl_found
  end do

  end subroutine raytrace_to_dist_column_clump
  !===========================================================================

!--- accumulate grid%Jmu(xfreq, mu) at escape (mu = cos(theta_z)).
  subroutine add_to_Jmu_clump(photon, grid, ix)
  use define
  type(photon_type), intent(in)    :: photon
  type(grid_type),   intent(inout) :: grid
  integer,           intent(in)    :: ix
  integer  :: imu
  real(wp) :: mu_val
  if (.not. associated(grid%Jmu)) return
  mu_val = photon%kz
  if (par%xyz_symmetry) mu_val = abs(mu_val)
  imu = floor((mu_val - par%mu_min)/par%dmu) + 1
  if (imu < 1)       imu = 1
  if (imu > par%nmu) imu = par%nmu
  !$OMP ATOMIC UPDATE
  grid%Jmu(ix, imu) = grid%Jmu(ix, imu) + photon%wgt
  end subroutine add_to_Jmu_clump
  !===========================================================================

  !===========================================================================
  ! Capped variant of raytrace_to_edge_clump for photon-transport callers
  ! (peeling-off).  Same algorithm; bails out as soon as the running tau
  ! exceeds tau_max, since exp(-tau) underflows in double precision above
  ! tau_max ~ 745 and any further accumulation contributes nothing.  The
  ! uncapped raytrace_to_edge_clump above keeps its 3-arg signature so that
  ! the raytrace_to_edge procedure pointer (define.f90 :: raytrace_edge
  ! interface) still matches.
  !===========================================================================
  subroutine raytrace_to_edge_clump_capped(photon0, grid, tau, tau_max)
  use define
  implicit none
  type(photon_type), intent(in)  :: photon0
  type(grid_type),   intent(in)  :: grid
  real(kind=wp),     intent(out) :: tau
  real(kind=wp),     intent(in)  :: tau_max

  real(kind=wp)  :: xp, yp, zp, kx, ky, kz
  real(kind=wp)  :: xfreq_loc
  real(kind=wp)  :: t_sp, t_seg, te, tx2, u_los, kap
  integer(int64) :: icl_cur, icl_found
  logical        :: found

  tau      = 0.0_wp
  xp       = photon0%x;  yp = photon0%y;  zp = photon0%z
  kx       = photon0%kx; ky = photon0%ky; kz = photon0%kz
  xfreq_loc = photon0%xfreq
  icl_cur   = int(photon0%icell_clump, int64)

  if (icl_cur > 0) then
     t_seg = clump_exit_dist(xp, yp, zp, kx, ky, kz, icl_cur)
     kap   = cl_rhokap(icl_cur) * voigt_clump(xfreq_loc, icl_cur)
     tau   = tau + kap * t_seg
     if (tau >= tau_max) return
     xp    = xp + t_seg * kx
     yp    = yp + t_seg * ky
     zp    = zp + t_seg * kz
     u_los = ulos_clump(icl_cur, kx, ky, kz)
     xfreq_loc = xfreq_loc + u_los
     if (xp**2 + yp**2 + zp**2 >= sphere_R**2) return
  end if

  do
     t_sp = sphere_exit_dist(xp, yp, zp, kx, ky, kz)
     if (t_sp <= 0.0_wp) exit
     call find_next_clump(xp, yp, zp, kx, ky, kz, icl_cur, t_sp, &
                           te, tx2, icl_found, found)
     if (.not. found) exit
     te = max(0.0_wp, te)
     xp = xp + te * kx;  yp = yp + te * ky;  zp = zp + te * kz
     u_los = ulos_clump(icl_found, kx, ky, kz)
     xfreq_loc = xfreq_loc - u_los
     t_seg = clump_exit_dist(xp, yp, zp, kx, ky, kz, icl_found)
     kap   = cl_rhokap(icl_found) * voigt_clump(xfreq_loc, icl_found)
     tau   = tau + kap * t_seg
     if (tau >= tau_max) return
     xp    = xp + t_seg * kx;  yp = yp + t_seg * ky;  zp = zp + t_seg * kz
     xfreq_loc = xfreq_loc + u_los
     icl_cur = icl_found
     if (xp**2 + yp**2 + zp**2 >= sphere_R**2) exit
  end do
  end subroutine raytrace_to_edge_clump_capped
  !===========================================================================

  !===========================================================================
  ! Capped variant of raytrace_to_dist_clump (same conventions as
  ! raytrace_to_edge_clump_capped above).
  !===========================================================================
  subroutine raytrace_to_dist_clump_capped(photon0, grid, dist_in, tau, tau_max)
  use define
  implicit none
  type(photon_type), intent(in)  :: photon0
  type(grid_type),   intent(in)  :: grid
  real(kind=wp),     intent(in)  :: dist_in
  real(kind=wp),     intent(out) :: tau
  real(kind=wp),     intent(in)  :: tau_max

  real(kind=wp)  :: xp, yp, zp, kx, ky, kz, xfreq_loc
  real(kind=wp)  :: t_rem, t_seg, te, tx2, u_los, kap
  integer(int64) :: icl_cur, icl_found
  logical        :: found

  tau = 0.0_wp
  xp  = photon0%x;  yp = photon0%y;  zp = photon0%z
  kx  = photon0%kx; ky = photon0%ky; kz = photon0%kz
  xfreq_loc = photon0%xfreq
  icl_cur   = int(photon0%icell_clump, int64)
  t_rem     = dist_in

  if (icl_cur > 0) then
     t_seg = min(clump_exit_dist(xp, yp, zp, kx, ky, kz, icl_cur), t_rem)
     kap   = cl_rhokap(icl_cur) * voigt_clump(xfreq_loc, icl_cur)
     tau   = tau + kap * t_seg
     if (tau >= tau_max) return
     t_rem = t_rem - t_seg
     if (t_rem <= 0.0_wp) return
     xp = xp + t_seg*kx;  yp = yp + t_seg*ky;  zp = zp + t_seg*kz
     u_los = ulos_clump(icl_cur, kx, ky, kz)
     xfreq_loc = xfreq_loc + u_los
  end if

  do while (t_rem > 0.0_wp)
     call find_next_clump(xp, yp, zp, kx, ky, kz, icl_cur, t_rem, &
                           te, tx2, icl_found, found)
     if (.not. found) exit
     te = max(0.0_wp, te)
     xp = xp + te*kx;  yp = yp + te*ky;  zp = zp + te*kz
     t_rem = t_rem - te
     u_los = ulos_clump(icl_found, kx, ky, kz)
     xfreq_loc = xfreq_loc - u_los
     t_seg = min(clump_exit_dist(xp, yp, zp, kx, ky, kz, icl_found), t_rem)
     kap   = cl_rhokap(icl_found) * voigt_clump(xfreq_loc, icl_found)
     tau   = tau + kap * t_seg
     if (tau >= tau_max) return
     t_rem = t_rem - t_seg
     if (t_rem <= 0.0_wp) return
     xp = xp + t_seg*kx;  yp = yp + t_seg*ky;  zp = zp + t_seg*kz
     xfreq_loc = xfreq_loc + u_los
     icl_cur = icl_found
  end do
  end subroutine raytrace_to_dist_clump_capped

  !=============================================================================
  !  OVERLAP-AWARE RAYTRACE ROUTINES
  !  Used when has_overlap=.true. (file-loaded clumps that may overlap).
  !
  !  Frequency convention: photon%xfreq stays in the GLOBAL rest frame throughout
  !  the spatial raytrace. Doppler shifts are applied only at scatter time by the
  !  scatter_resonance_clump_nostokes / _stokes wrappers in scattering_car.f90.
  !=============================================================================

  !===========================================================================
  ! Internal helper: compute kap_total = sum of rhokap*H(x_global, a) over all
  ! clumps in the active set, where each H is evaluated at the clump-frame
  ! frequency x_clump = x_global - u_los_i.
  !===========================================================================
  real(kind=wp) function sum_kap_active(active, n_active, xfreq_g, kx, ky, kz)
  use define, only: wp
  integer(int64), intent(in) :: active(:)
  integer,        intent(in) :: n_active
  real(kind=wp),  intent(in) :: xfreq_g, kx, ky, kz
  integer :: m
  integer(int64) :: icl
  real(kind=wp)  :: u_los, x_cl
  sum_kap_active = 0.0_wp
  do m = 1, n_active
     icl   = active(m)
     u_los = ulos_clump(icl, kx, ky, kz)
     x_cl  = xfreq_g - u_los
     sum_kap_active = sum_kap_active + cl_rhokap(icl) * voigt_clump(x_cl, icl)
  end do
  end function sum_kap_active
  !===========================================================================

  !===========================================================================
  ! Internal helper: opacity-weighted sampling of the owner clump from the
  ! active set.  Returns the 1-based clump index.
  ! Requires n_active >= 1 and kap_total > 0.
  !===========================================================================
  integer(int64) function sample_owner_clump(active, n_active, xfreq_g, kx, ky, kz)
  use define,  only: wp
  use random,  only: rand_number
  integer(int64), intent(in) :: active(:)
  integer,        intent(in) :: n_active
  real(kind=wp),  intent(in) :: xfreq_g, kx, ky, kz
  integer :: m
  integer(int64) :: icl
  real(kind=wp)  :: u_los, x_cl, kap, cumul, rnd
  cumul = 0.0_wp
  rnd   = rand_number()
  do m = 1, n_active
     icl   = active(m)
     u_los = ulos_clump(icl, kx, ky, kz)
     x_cl  = xfreq_g - u_los
     kap   = cl_rhokap(icl) * voigt_clump(x_cl, icl)
     cumul = cumul + kap
     sample_owner_clump = icl
     if (rnd * sum_kap_active(active, n_active, xfreq_g, kx, ky, kz) <= cumul) return
  end do
  end function sample_owner_clump
  !===========================================================================

  !===========================================================================
  subroutine raytrace_to_tau_clump_overlap(photon, grid, tau_in)
  !---------------------------------------------------------------------------
  ! Overlap-aware advance to optical depth tau_in.
  ! photon%xfreq is kept in the global frame; the active set is tracked across
  ! ENTER/EXIT events from collect_ray_events_overlap.
  !---------------------------------------------------------------------------
  use define
  implicit none
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid
  real(kind=wp),     intent(in)    :: tau_in

  integer, parameter :: MAX_EVT = 2048
  real(kind=wp)   :: ev_t(MAX_EVT)
  integer(int64)  :: ev_icl(MAX_EVT), active(MAX_EVT)
  integer         :: ev_type(MAX_EVT), n_ev, n_active

  real(kind=wp)  :: kx, ky, kz, tau_rem, kap_tot, ds, t_sp, t_cur, t_next, dt
  real(kind=wp)  :: xfreq_g
  integer        :: ie, ix_jout

  kx      = photon%kx;  ky = photon%ky;  kz = photon%kz
  xfreq_g = photon%xfreq    ! global frame — unchanged throughout raytrace
  tau_rem = tau_in

  do while (photon%inside)

     t_sp = sphere_exit_dist(photon%x, photon%y, photon%z, kx, ky, kz)
     if (t_sp <= 0.0_wp) then
        photon%inside = .false.
        call update_cell_idx(photon, grid)
        ix_jout = floor((xfreq_g - grid%xfreq_min)/grid%dxfreq) + 1
        if (ix_jout >= 1 .and. ix_jout <= grid%nxfreq) then
           !$OMP ATOMIC UPDATE
           grid%Jout(ix_jout) = grid%Jout(ix_jout) + photon%wgt
           if (par%save_Jmu) call add_to_Jmu_clump(photon, grid, ix_jout)
        end if
        return
     end if

     !--- collect all events along this ray segment to sphere boundary
     call collect_ray_events_overlap(photon%x, photon%y, photon%z, kx, ky, kz, &
          t_sp, ev_t, ev_icl, ev_type, n_ev)
     !--- initialize active set at current position
     call active_set_at_point(photon%x, photon%y, photon%z, active, n_active)

     t_cur = 0.0_wp

     do ie = 1, n_ev + 1
        if (ie <= n_ev) then
           t_next = ev_t(ie)
        else
           t_next = t_sp
        end if
        t_next = min(t_next, t_sp)
        dt = t_next - t_cur
        if (dt <= 0.0_wp) then
           if (ie <= n_ev) call apply_event(ev_icl(ie), ev_type(ie), active, n_active)
           cycle
        end if

        kap_tot = sum_kap_active(active, n_active, xfreq_g, kx, ky, kz)

        if (kap_tot > 0.0_wp .and. tau_rem <= kap_tot * dt) then
           !--- scatter inside this segment
           ds = tau_rem / kap_tot
           photon%x = photon%x + (t_cur + ds) * kx
           photon%y = photon%y + (t_cur + ds) * ky
           photon%z = photon%z + (t_cur + ds) * kz
           photon%xfreq       = xfreq_g
           photon%icell_clump = int(sample_owner_clump(active, n_active, xfreq_g, kx, ky, kz))
           call update_cell_idx(photon, grid)
           return
        end if

        tau_rem = tau_rem - kap_tot * dt
        t_cur   = t_next

        if (ie <= n_ev) call apply_event(ev_icl(ie), ev_type(ie), active, n_active)
     end do

     !--- photon has traversed to sphere exit
     photon%x = photon%x + t_sp * kx
     photon%y = photon%y + t_sp * ky
     photon%z = photon%z + t_sp * kz
     photon%xfreq       = xfreq_g
     photon%icell_clump = 0
     photon%inside = .false.
     call update_cell_idx(photon, grid)
     ix_jout = floor((xfreq_g - grid%xfreq_min)/grid%dxfreq) + 1
     if (ix_jout >= 1 .and. ix_jout <= grid%nxfreq) then
        !$OMP ATOMIC UPDATE
        grid%Jout(ix_jout) = grid%Jout(ix_jout) + photon%wgt
        if (par%save_Jmu) call add_to_Jmu_clump(photon, grid, ix_jout)
     end if
     return

  end do

  contains

    subroutine apply_event(icl_ev, etype, aset, na)
    integer(int64), intent(in)    :: icl_ev
    integer,        intent(in)    :: etype
    integer(int64), intent(inout) :: aset(:)
    integer,        intent(inout) :: na
    integer :: mm
    if (etype == +1) then
       if (na < size(aset)) then
          na = na + 1;  aset(na) = icl_ev
       end if
    else
       do mm = 1, na
          if (aset(mm) == icl_ev) then
             aset(mm) = aset(na);  na = na - 1;  return
          end if
       end do
    end if
    end subroutine apply_event

  end subroutine raytrace_to_tau_clump_overlap
  !===========================================================================

  !===========================================================================
  subroutine raytrace_to_edge_clump_overlap(photon0, grid, tau)
  !---------------------------------------------------------------------------
  ! Overlap-aware sight-line tau from photon0 position to sphere edge.
  ! photon0 read-only; xfreq kept in global frame.
  !---------------------------------------------------------------------------
  use define
  implicit none
  type(photon_type), intent(in)  :: photon0
  type(grid_type),   intent(in)  :: grid
  real(kind=wp),     intent(out) :: tau

  integer, parameter :: MAX_EVT = 2048
  real(kind=wp)   :: ev_t(MAX_EVT)
  integer(int64)  :: ev_icl(MAX_EVT), active(MAX_EVT)
  integer         :: ev_type(MAX_EVT), n_ev, n_active

  real(kind=wp)  :: kx, ky, kz, t_sp, t_cur, t_next, dt, kap_tot, xfreq_g
  integer        :: ie

  tau     = 0.0_wp
  kx      = photon0%kx;  ky = photon0%ky;  kz = photon0%kz
  xfreq_g = photon0%xfreq

  t_sp = sphere_exit_dist(photon0%x, photon0%y, photon0%z, kx, ky, kz)
  if (t_sp <= 0.0_wp) return

  call collect_ray_events_overlap(photon0%x, photon0%y, photon0%z, kx, ky, kz, &
       t_sp, ev_t, ev_icl, ev_type, n_ev)
  call active_set_at_point(photon0%x, photon0%y, photon0%z, active, n_active)

  t_cur = 0.0_wp
  do ie = 1, n_ev + 1
     t_next  = merge(ev_t(ie), t_sp, ie <= n_ev)
     t_next  = min(t_next, t_sp)
     dt      = t_next - t_cur
     if (dt > 0.0_wp) then
        kap_tot = sum_kap_active(active, n_active, xfreq_g, kx, ky, kz)
        tau     = tau + kap_tot * dt
     end if
     t_cur = t_next
     if (ie <= n_ev) then
        call apply_event_local(ev_icl(ie), ev_type(ie), active, n_active)
     end if
  end do

  contains

    subroutine apply_event_local(icl_ev, etype, aset, na)
    integer(int64), intent(in)    :: icl_ev
    integer,        intent(in)    :: etype
    integer(int64), intent(inout) :: aset(:)
    integer,        intent(inout) :: na
    integer :: mm
    if (etype == +1) then
       if (na < size(aset)) then; na = na + 1;  aset(na) = icl_ev; end if
    else
       do mm = 1, na
          if (aset(mm) == icl_ev) then; aset(mm) = aset(na);  na = na - 1;  return; end if
       end do
    end if
    end subroutine apply_event_local

  end subroutine raytrace_to_edge_clump_overlap
  !===========================================================================

  !===========================================================================
  subroutine raytrace_to_edge_clump_overlap_capped(photon0, grid, tau, tau_max)
  use define
  implicit none
  type(photon_type), intent(in)  :: photon0
  type(grid_type),   intent(in)  :: grid
  real(kind=wp),     intent(out) :: tau
  real(kind=wp),     intent(in)  :: tau_max

  integer, parameter :: MAX_EVT = 2048
  real(kind=wp)   :: ev_t(MAX_EVT)
  integer(int64)  :: ev_icl(MAX_EVT), active(MAX_EVT)
  integer         :: ev_type(MAX_EVT), n_ev, n_active

  real(kind=wp)  :: kx, ky, kz, t_sp, t_cur, t_next, dt, kap_tot, xfreq_g
  integer        :: ie

  tau     = 0.0_wp
  kx      = photon0%kx;  ky = photon0%ky;  kz = photon0%kz
  xfreq_g = photon0%xfreq

  t_sp = sphere_exit_dist(photon0%x, photon0%y, photon0%z, kx, ky, kz)
  if (t_sp <= 0.0_wp) return

  call collect_ray_events_overlap(photon0%x, photon0%y, photon0%z, kx, ky, kz, &
       t_sp, ev_t, ev_icl, ev_type, n_ev)
  call active_set_at_point(photon0%x, photon0%y, photon0%z, active, n_active)

  t_cur = 0.0_wp
  do ie = 1, n_ev + 1
     t_next  = merge(ev_t(ie), t_sp, ie <= n_ev)
     t_next  = min(t_next, t_sp)
     dt      = t_next - t_cur
     if (dt > 0.0_wp) then
        kap_tot = sum_kap_active(active, n_active, xfreq_g, kx, ky, kz)
        tau     = tau + kap_tot * dt
        if (tau >= tau_max) return
     end if
     t_cur = t_next
     if (ie <= n_ev) then
        call apply_event_local(ev_icl(ie), ev_type(ie), active, n_active)
     end if
  end do

  contains

    subroutine apply_event_local(icl_ev, etype, aset, na)
    integer(int64), intent(in)    :: icl_ev
    integer,        intent(in)    :: etype
    integer(int64), intent(inout) :: aset(:)
    integer,        intent(inout) :: na
    integer :: mm
    if (etype == +1) then
       if (na < size(aset)) then; na = na + 1;  aset(na) = icl_ev; end if
    else
       do mm = 1, na
          if (aset(mm) == icl_ev) then; aset(mm) = aset(na);  na = na - 1;  return; end if
       end do
    end if
    end subroutine apply_event_local

  end subroutine raytrace_to_edge_clump_overlap_capped
  !===========================================================================

  !===========================================================================
  subroutine raytrace_to_edge_tau_gas_clump_overlap(photon0, grid, tau_gas)
  use define
  implicit none
  type(photon_type), intent(in)  :: photon0
  type(grid_type),   intent(in)  :: grid
  real(kind=wp),     intent(out) :: tau_gas
  call raytrace_to_edge_clump_overlap(photon0, grid, tau_gas)
  end subroutine raytrace_to_edge_tau_gas_clump_overlap
  !===========================================================================

  !===========================================================================
  subroutine raytrace_to_edge_column_clump_overlap(photon0, grid, N_gas, tau_dust)
  !---------------------------------------------------------------------------
  ! Column density version: accumulates N_HI * ds summed over active clumps.
  !---------------------------------------------------------------------------
  use define
  implicit none
  type(photon_type), intent(in)  :: photon0
  type(grid_type),   intent(in)  :: grid
  real(kind=wp),     intent(out) :: N_gas, tau_dust

  integer, parameter :: MAX_EVT = 2048
  real(kind=wp)   :: ev_t(MAX_EVT)
  integer(int64)  :: ev_icl(MAX_EVT), active(MAX_EVT)
  integer         :: ev_type(MAX_EVT), n_ev, n_active

  real(kind=wp)  :: kx, ky, kz, t_sp, t_cur, t_next, dt
  integer        :: ie, m
  integer(int64) :: icl

  N_gas    = 0.0_wp
  tau_dust = 0.0_wp
  kx = photon0%kx;  ky = photon0%ky;  kz = photon0%kz

  t_sp = sphere_exit_dist(photon0%x, photon0%y, photon0%z, kx, ky, kz)
  if (t_sp <= 0.0_wp) return

  call collect_ray_events_overlap(photon0%x, photon0%y, photon0%z, kx, ky, kz, &
       t_sp, ev_t, ev_icl, ev_type, n_ev)
  call active_set_at_point(photon0%x, photon0%y, photon0%z, active, n_active)

  t_cur = 0.0_wp
  do ie = 1, n_ev + 1
     t_next = merge(ev_t(ie), t_sp, ie <= n_ev)
     t_next = min(t_next, t_sp)
     dt     = t_next - t_cur
     if (dt > 0.0_wp) then
        do m = 1, n_active
           icl   = active(m)
           N_gas = N_gas + (cl_rhokap(icl) / line%cross0) * cl_Dfreq(icl) * dt
        end do
     end if
     t_cur = t_next
     if (ie <= n_ev) then
        call apply_event_local(ev_icl(ie), ev_type(ie), active, n_active)
     end if
  end do

  contains

    subroutine apply_event_local(icl_ev, etype, aset, na)
    integer(int64), intent(in)    :: icl_ev
    integer,        intent(in)    :: etype
    integer(int64), intent(inout) :: aset(:)
    integer,        intent(inout) :: na
    integer :: mm
    if (etype == +1) then
       if (na < size(aset)) then; na = na + 1;  aset(na) = icl_ev; end if
    else
       do mm = 1, na
          if (aset(mm) == icl_ev) then; aset(mm) = aset(na);  na = na - 1;  return; end if
       end do
    end if
    end subroutine apply_event_local

  end subroutine raytrace_to_edge_column_clump_overlap
  !===========================================================================

  !===========================================================================
  subroutine raytrace_to_dist_clump_overlap(photon0, grid, dist_in, tau)
  !---------------------------------------------------------------------------
  ! Overlap-aware tau accumulation over a fixed path length dist_in.
  !---------------------------------------------------------------------------
  use define
  implicit none
  type(photon_type), intent(in)  :: photon0
  type(grid_type),   intent(in)  :: grid
  real(kind=wp),     intent(in)  :: dist_in
  real(kind=wp),     intent(out) :: tau

  integer, parameter :: MAX_EVT = 2048
  real(kind=wp)   :: ev_t(MAX_EVT)
  integer(int64)  :: ev_icl(MAX_EVT), active(MAX_EVT)
  integer         :: ev_type(MAX_EVT), n_ev, n_active

  real(kind=wp)  :: kx, ky, kz, t_cur, t_next, dt, kap_tot, xfreq_g
  integer        :: ie

  tau     = 0.0_wp
  kx      = photon0%kx;  ky = photon0%ky;  kz = photon0%kz
  xfreq_g = photon0%xfreq

  call collect_ray_events_overlap(photon0%x, photon0%y, photon0%z, kx, ky, kz, &
       dist_in, ev_t, ev_icl, ev_type, n_ev)
  call active_set_at_point(photon0%x, photon0%y, photon0%z, active, n_active)

  t_cur = 0.0_wp
  do ie = 1, n_ev + 1
     t_next  = merge(ev_t(ie), dist_in, ie <= n_ev)
     t_next  = min(t_next, dist_in)
     dt      = t_next - t_cur
     if (dt > 0.0_wp) then
        kap_tot = sum_kap_active(active, n_active, xfreq_g, kx, ky, kz)
        tau     = tau + kap_tot * dt
     end if
     t_cur = t_next
     if (ie <= n_ev) then
        call apply_event_local(ev_icl(ie), ev_type(ie), active, n_active)
     end if
  end do

  contains

    subroutine apply_event_local(icl_ev, etype, aset, na)
    integer(int64), intent(in)    :: icl_ev
    integer,        intent(in)    :: etype
    integer(int64), intent(inout) :: aset(:)
    integer,        intent(inout) :: na
    integer :: mm
    if (etype == +1) then
       if (na < size(aset)) then; na = na + 1;  aset(na) = icl_ev; end if
    else
       do mm = 1, na
          if (aset(mm) == icl_ev) then; aset(mm) = aset(na);  na = na - 1;  return; end if
       end do
    end if
    end subroutine apply_event_local

  end subroutine raytrace_to_dist_clump_overlap
  !===========================================================================

  !===========================================================================
  subroutine raytrace_to_dist_clump_overlap_capped(photon0, grid, dist_in, tau, tau_max)
  use define
  implicit none
  type(photon_type), intent(in)  :: photon0
  type(grid_type),   intent(in)  :: grid
  real(kind=wp),     intent(in)  :: dist_in
  real(kind=wp),     intent(out) :: tau
  real(kind=wp),     intent(in)  :: tau_max

  integer, parameter :: MAX_EVT = 2048
  real(kind=wp)   :: ev_t(MAX_EVT)
  integer(int64)  :: ev_icl(MAX_EVT), active(MAX_EVT)
  integer         :: ev_type(MAX_EVT), n_ev, n_active

  real(kind=wp)  :: kx, ky, kz, t_cur, t_next, dt, kap_tot, xfreq_g
  integer        :: ie

  tau     = 0.0_wp
  kx      = photon0%kx;  ky = photon0%ky;  kz = photon0%kz
  xfreq_g = photon0%xfreq

  call collect_ray_events_overlap(photon0%x, photon0%y, photon0%z, kx, ky, kz, &
       dist_in, ev_t, ev_icl, ev_type, n_ev)
  call active_set_at_point(photon0%x, photon0%y, photon0%z, active, n_active)

  t_cur = 0.0_wp
  do ie = 1, n_ev + 1
     t_next  = merge(ev_t(ie), dist_in, ie <= n_ev)
     t_next  = min(t_next, dist_in)
     dt      = t_next - t_cur
     if (dt > 0.0_wp) then
        kap_tot = sum_kap_active(active, n_active, xfreq_g, kx, ky, kz)
        tau     = tau + kap_tot * dt
        if (tau >= tau_max) return
     end if
     t_cur = t_next
     if (ie <= n_ev) then
        call apply_event_local(ev_icl(ie), ev_type(ie), active, n_active)
     end if
  end do

  contains

    subroutine apply_event_local(icl_ev, etype, aset, na)
    integer(int64), intent(in)    :: icl_ev
    integer,        intent(in)    :: etype
    integer(int64), intent(inout) :: aset(:)
    integer,        intent(inout) :: na
    integer :: mm
    if (etype == +1) then
       if (na < size(aset)) then; na = na + 1;  aset(na) = icl_ev; end if
    else
       do mm = 1, na
          if (aset(mm) == icl_ev) then; aset(mm) = aset(na);  na = na - 1;  return; end if
       end do
    end if
    end subroutine apply_event_local

  end subroutine raytrace_to_dist_clump_overlap_capped
  !===========================================================================

  !===========================================================================
  subroutine raytrace_to_dist_tau_gas_clump_overlap(photon0, grid, dist_in, tau_gas)
  use define
  implicit none
  type(photon_type), intent(in)  :: photon0
  type(grid_type),   intent(in)  :: grid
  real(kind=wp),     intent(in)  :: dist_in
  real(kind=wp),     intent(out) :: tau_gas
  call raytrace_to_dist_clump_overlap(photon0, grid, dist_in, tau_gas)
  end subroutine raytrace_to_dist_tau_gas_clump_overlap
  !===========================================================================

  !===========================================================================
  subroutine raytrace_to_dist_column_clump_overlap(photon0, grid, dist_in, N_gas, tau_dust)
  use define
  implicit none
  type(photon_type), intent(in)  :: photon0
  type(grid_type),   intent(in)  :: grid
  real(kind=wp),     intent(in)  :: dist_in
  real(kind=wp),     intent(out) :: N_gas, tau_dust

  integer, parameter :: MAX_EVT = 2048
  real(kind=wp)   :: ev_t(MAX_EVT)
  integer(int64)  :: ev_icl(MAX_EVT), active(MAX_EVT)
  integer         :: ev_type(MAX_EVT), n_ev, n_active

  real(kind=wp)  :: kx, ky, kz, t_cur, t_next, dt
  integer        :: ie, m
  integer(int64) :: icl

  N_gas    = 0.0_wp
  tau_dust = 0.0_wp
  kx = photon0%kx;  ky = photon0%ky;  kz = photon0%kz

  call collect_ray_events_overlap(photon0%x, photon0%y, photon0%z, kx, ky, kz, &
       dist_in, ev_t, ev_icl, ev_type, n_ev)
  call active_set_at_point(photon0%x, photon0%y, photon0%z, active, n_active)

  t_cur = 0.0_wp
  do ie = 1, n_ev + 1
     t_next = merge(ev_t(ie), dist_in, ie <= n_ev)
     t_next = min(t_next, dist_in)
     dt     = t_next - t_cur
     if (dt > 0.0_wp) then
        do m = 1, n_active
           icl   = active(m)
           N_gas = N_gas + (cl_rhokap(icl) / line%cross0) * cl_Dfreq(icl) * dt
        end do
     end if
     t_cur = t_next
     if (ie <= n_ev) then
        call apply_event_local(ev_icl(ie), ev_type(ie), active, n_active)
     end if
  end do

  contains

    subroutine apply_event_local(icl_ev, etype, aset, na)
    integer(int64), intent(in)    :: icl_ev
    integer,        intent(in)    :: etype
    integer(int64), intent(inout) :: aset(:)
    integer,        intent(inout) :: na
    integer :: mm
    if (etype == +1) then
       if (na < size(aset)) then; na = na + 1;  aset(na) = icl_ev; end if
    else
       do mm = 1, na
          if (aset(mm) == icl_ev) then; aset(mm) = aset(na);  na = na - 1;  return; end if
       end do
    end if
    end subroutine apply_event_local

  end subroutine raytrace_to_dist_column_clump_overlap
  !===========================================================================

end module raytrace_clump_mod
