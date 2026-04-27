module raytrace_clump_mod
!---------------------------------------------------------------------------
! Raytrace routines for the clump medium.
!
! The medium consists of N spherical clumps (radius cl_radius) placed inside
! a sphere of radius sphere_R. Outside the clumps the medium is vacuum.
! Each clump has uniform opacity cl_rhokap and Voigt parameter cl_voigt_a.
!
! Frequency convention (same as Cartesian):
!   photon%xfreq is in units of cl_Dfreq (= Dfreq_ref).
!   cl_vx/y/z are stored as v / cl_vtherm (dimensionless), set in init_clumps.
!   At clump entry: xfreq -= v_los_clump   (v already in vtherm units)
!   At clump exit:  xfreq += v_los_clump
!
! photon%icell_clump: 0 = in vacuum; > 0 = index of current clump (1-based).
! photon%icell/jcell/kcell: maintained using floor() map from position.
!---------------------------------------------------------------------------
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

        kap = cl_rhokap * voigt(photon%xfreq, cl_voigt_a)

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
        u_los = cl_vx(icl)*kx + cl_vy(icl)*ky + cl_vz(icl)*kz
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
           u_los = cl_vx(icl_found)*kx + cl_vy(icl_found)*ky + &
                   cl_vz(icl_found)*kz
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
     kap   = cl_rhokap * voigt(xfreq_loc, cl_voigt_a)
     tau   = tau + kap * t_seg
     xp    = xp + t_seg * kx
     yp    = yp + t_seg * ky
     zp    = zp + t_seg * kz
     u_los = cl_vx(icl_cur)*kx + cl_vy(icl_cur)*ky + cl_vz(icl_cur)*kz
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
     u_los = cl_vx(icl_found)*kx + cl_vy(icl_found)*ky + cl_vz(icl_found)*kz
     xfreq_loc = xfreq_loc - u_los

     !--- tau through this clump
     t_seg = clump_exit_dist(xp, yp, zp, kx, ky, kz, icl_found)
     kap   = cl_rhokap * voigt(xfreq_loc, cl_voigt_a)
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
     N_gas = N_gas + (cl_rhokap / line%cross0) * t_seg
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
     N_gas = N_gas + (cl_rhokap / line%cross0) * t_seg
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
     kap   = cl_rhokap * voigt(xfreq_loc, cl_voigt_a)
     tau   = tau + kap * t_seg
     t_rem = t_rem - t_seg
     if (t_rem <= 0.0_wp) return
     xp = xp + t_seg*kx;  yp = yp + t_seg*ky;  zp = zp + t_seg*kz
     u_los = cl_vx(icl_cur)*kx + cl_vy(icl_cur)*ky + cl_vz(icl_cur)*kz
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

     u_los = cl_vx(icl_found)*kx + cl_vy(icl_found)*ky + cl_vz(icl_found)*kz
     xfreq_loc = xfreq_loc - u_los

     t_seg = min(clump_exit_dist(xp, yp, zp, kx, ky, kz, icl_found), t_rem)
     kap   = cl_rhokap * voigt(xfreq_loc, cl_voigt_a)
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
     N_gas = N_gas + (cl_rhokap / line%cross0) * t_seg
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
     N_gas = N_gas + (cl_rhokap / line%cross0) * t_seg
     t_rem = t_rem - t_seg
     if (t_rem <= 0.0_wp) return
     xp = xp + t_seg*kx;  yp = yp + t_seg*ky;  zp = zp + t_seg*kz
     icl_cur = icl_found
  end do

  end subroutine raytrace_to_dist_column_clump
  !===========================================================================

end module raytrace_clump_mod
