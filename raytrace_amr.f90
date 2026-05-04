module raytrace_amr_mod
  !-----------------------------------------------------------------------
  ! AMR raytrace routines: octree-traversal equivalents of the Cartesian
  ! raytrace_car_v2h_refactored.f90 routines.
  !
  ! Coordinate convention (same as Cartesian):
  !   All positions and path lengths are in CODE UNITS (kpc, pc, etc.).
  !   rhokap is per code unit  →  tau = rhokap * ds   is dimensionless.
  !
  ! Cell-traversal convention:
  !   x, y, z  = running absolute position (updated at every face crossing).
  !   iface    = face index returned by amr_cell_exit:
  !               1=+x  2=-x  3=+y  4=-y  5=+z  6=-z
  !
  ! Cell lookup:
  !   amr_next_leaf(icell, iface, x, y, z) — O(1) via precomputed neighbor
  !   table; no position-based tree traversal from root at every step.
  !
  ! Frequency-shift convention (same as Cartesian):
  !   xfreq in the comoving frame of the current cell.
  !   At each cell crossing:
  !     xfreq_new = (xfreq_old + u_old) * Dfreq_old / Dfreq_new - u_new
  !   where u = dot(v_cell, k_hat)  [velocity in units of v_thermal].
  !
  ! Boundary handling at face exits when amr_next_leaf returns no neighbor:
  !
  !   par%xy_periodic       (e.g. plane_atmosphere):
  !     iface = 1 (+x) -> x = x - amr_grid%xrange,  re-find leaf
  !     iface = 2 (-x) -> x = x + amr_grid%xrange,  re-find leaf
  !     iface = 3 (+y) -> y = y - amr_grid%yrange,  re-find leaf
  !     iface = 4 (-y) -> y = y + amr_grid%yrange,  re-find leaf
  !     iface = 5,6   -> photon escapes
  !
  !   par%xy_symmetry       (mirror at x=xmin and y=ymin):
  !     iface = 2 (-x) -> kx -> -kx, stay in same leaf (specular reflection)
  !     iface = 4 (-y) -> ky -> -ky, stay in same leaf
  !     other faces   -> photon escapes
  !
  !   par%xyz_symmetry      (mirror at x=xmin, y=ymin, z=zmin):
  !     iface = 2,4,6 -> reflect kx/ky/kz, stay in same leaf
  !     iface = 1,3,5 -> photon escapes
  !
  ! Reflection assumption: no AMR leaf straddles the reflective plane (i.e.
  ! leaves are aligned so the symmetry plane coincides exactly with leaf
  ! faces at amr_grid%xmin / ymin / zmin).  Lab-frame frequency is
  ! preserved across reflection: u1 is recomputed with the new kx/ky/kz and
  ! photon%xfreq is shifted to xfreq_lab - u1_new (no Dfreq change).
  !
  ! Final photon position:
  !   raytrace_to_tau_amr uses the running x/y/z directly (advanced at face
  !   crossings and at the in-cell overshoot exit), so wraps and reflections
  !   do not need a separate offset accumulator.
  !-----------------------------------------------------------------------
  use octree_mod
  use voigt_mod
  implicit none
  private

  public :: raytrace_to_tau_amr
  public :: raytrace_to_edge_amr
  public :: raytrace_to_edge_tau_gas_amr
  public :: raytrace_to_edge_column_amr
  public :: raytrace_to_dist_amr
  public :: raytrace_to_dist_tau_gas_amr
  public :: raytrace_to_dist_column_amr

  real(wp), parameter :: tau_huge = 745.2_wp  ! exp(-tau_huge) ≈ 0 in double

contains

  !=========================================================================
  ! Walk photon through AMR grid to optical depth tau_in.
  ! Updates photon position, cell index, frequency, and inside flag.
  !=========================================================================
  subroutine raytrace_to_tau_amr(photon, grid, tau_in)
    use define
    implicit none
    type(photon_type), intent(inout) :: photon
    type(grid_type),   intent(inout) :: grid   ! unused for AMR physics; kept for interface
    real(wp),          intent(in)    :: tau_in

    integer  :: il, il_new, icell, ix
    real(wp) :: x, y, z, kx, ky, kz
    real(wp) :: tau, t_exit, d_step
    real(wp) :: rhokap, rhokapH
    real(wp) :: u1, u2, Df_old, Df_new, xfreq_ref
    integer  :: iface
    logical  :: reflected

    x  = photon%x;    y  = photon%y;    z  = photon%z
    kx = photon%kx;   ky = photon%ky;   kz = photon%kz
    il = photon%icell_amr

    if (il <= 0) then
      il = amr_find_leaf(x, y, z)
      if (il <= 0) then
        photon%inside = .false.
        return
      end if
    end if

    tau = 0.0_wp
    u1  = amr_grid%vfx(il)*kx + amr_grid%vfy(il)*ky + amr_grid%vfz(il)*kz

    do while (photon%inside)
      icell = amr_grid%icell_of_leaf(il)
      call amr_cell_exit(icell, x, y, z, kx, ky, kz, t_exit, iface)

      rhokapH = amr_grid%rhokap(il) * voigt(photon%xfreq, amr_grid%voigt_a(il))
      if (par%DGR > 0.0_wp .and. associated(amr_grid%rhokapD)) then
        rhokap = rhokapH + amr_grid%rhokapD(il)
      else
        rhokap = rhokapH
      end if

      if (tau + t_exit * rhokap >= tau_in) then
        ! Overshoot: stop inside this cell
        if (rhokap > 0.0_wp) then
          d_step = (tau_in - tau) / rhokap
        else
          d_step = t_exit
        end if
        x   = x + d_step * kx
        y   = y + d_step * ky
        z   = z + d_step * kz
        tau = tau_in
        exit
      end if

      ! Advance to face
      tau = tau + t_exit * rhokap
      x   = x + t_exit * kx
      y   = y + t_exit * ky
      z   = z + t_exit * kz

      ! Next leaf via precomputed neighbor table (O(1) + optional descent)
      il_new    = amr_next_leaf(icell, iface, x, y, z)
      reflected = .false.

      if (il_new <= 0) then
        if (par%xy_periodic) then
          select case (iface)
          case (1)
            x = x - amr_grid%xrange
          case (2)
            x = x + amr_grid%xrange
          case (3)
            y = y - amr_grid%yrange
          case (4)
            y = y + amr_grid%yrange
          end select
          if (iface <= 4) il_new = amr_find_leaf(x, y, z)
        else if (par%xyz_symmetry) then
          select case (iface)
          case (2);  kx = -kx; il_new = il; reflected = .true.
          case (4);  ky = -ky; il_new = il; reflected = .true.
          case (6);  kz = -kz; il_new = il; reflected = .true.
          end select
        else if (par%xy_symmetry) then
          select case (iface)
          case (2);  kx = -kx; il_new = il; reflected = .true.
          case (4);  ky = -ky; il_new = il; reflected = .true.
          end select
        end if
      end if

      if (il_new <= 0) then
        photon%inside = .false.
        exit
      end if

      ! Frequency shift across cell boundary (or reflection in same cell)
      Df_old = amr_grid%Dfreq(il)
      Df_new = amr_grid%Dfreq(il_new)
      u2     = amr_grid%vfx(il_new)*kx + amr_grid%vfy(il_new)*ky + amr_grid%vfz(il_new)*kz
      if (reflected) then
        ! Lab-frame frequency invariant: xfreq_lab = xfreq + u1 (old kx,ky,kz)
        ! After reflection, recompute u1 with new direction; same cell -> Df ratio = 1.
        photon%xfreq = (photon%xfreq + u1) - u2
      else
        photon%xfreq = (photon%xfreq + u1) * Df_old / Df_new - u2
      end if
      u1 = u2

      il = il_new
    end do

    ! Final state (running x/y/z and possibly-reflected kx/ky/kz)
    photon%x  = x;   photon%y  = y;   photon%z  = z
    photon%kx = kx;  photon%ky = ky;  photon%kz = kz
    photon%icell_amr = il

    if (.not. photon%inside) then
      ! Photon escaped: convert frequency to lab frame and record
      u1 = amr_grid%vfx(il)*kx + amr_grid%vfy(il)*ky + amr_grid%vfz(il)*kz
      photon%xfreq     = photon%xfreq + u1
      xfreq_ref        = photon%xfreq * (amr_grid%Dfreq(il) / amr_grid%Dfreq_ref)
      photon%xfreq_ref = xfreq_ref
      ix = floor((xfreq_ref - amr_grid%xfreq_min) / amr_grid%dxfreq) + 1
      if (ix >= 1 .and. ix <= amr_grid%nxfreq) then
        !$OMP ATOMIC UPDATE
        amr_grid%Jout(ix) = amr_grid%Jout(ix) + photon%wgt
        if (par%save_Jmu) call add_to_Jmu_amr(photon, ix)
      end if
    end if
  end subroutine raytrace_to_tau_amr

  !=========================================================================
  ! Integrate total (gas + dust) optical depth from photon to box edge.
  ! Photon state is NOT modified.
  !=========================================================================
  subroutine raytrace_to_edge_amr(photon0, grid, tau)
    use define
    implicit none
    type(photon_type), intent(in)  :: photon0
    type(grid_type),   intent(in)  :: grid
    real(wp),          intent(out) :: tau

    integer  :: il, il_new, icell
    real(wp) :: x, y, z, kx, ky, kz
    real(wp) :: t_exit
    real(wp) :: rhokap, u1, u2, Df_old, Df_new, xfreq_loc
    integer  :: iface
    logical  :: reflected

    x  = photon0%x;    y  = photon0%y;    z  = photon0%z
    kx = photon0%kx;   ky = photon0%ky;   kz = photon0%kz
    il = photon0%icell_amr

    tau = 0.0_wp
    if (il <= 0) then
      il = amr_find_leaf(x, y, z)
      if (il <= 0) return
    end if

    u1        = amr_grid%vfx(il)*kx + amr_grid%vfy(il)*ky + amr_grid%vfz(il)*kz
    xfreq_loc = photon0%xfreq

    do
      icell  = amr_grid%icell_of_leaf(il)
      call amr_cell_exit(icell, x, y, z, kx, ky, kz, t_exit, iface)
      rhokap = amr_grid%rhokap(il) * voigt(xfreq_loc, amr_grid%voigt_a(il))
      if (par%DGR > 0.0_wp .and. associated(amr_grid%rhokapD)) &
          rhokap = rhokap + amr_grid%rhokapD(il)
      tau = tau + t_exit * rhokap
      if (tau >= tau_huge) return

      x = x + t_exit * kx
      y = y + t_exit * ky
      z = z + t_exit * kz
      il_new    = amr_next_leaf(icell, iface, x, y, z)
      reflected = .false.

      if (il_new <= 0) then
        if (par%xy_periodic) then
          select case (iface)
          case (1); x = x - amr_grid%xrange
          case (2); x = x + amr_grid%xrange
          case (3); y = y - amr_grid%yrange
          case (4); y = y + amr_grid%yrange
          end select
          if (iface <= 4) il_new = amr_find_leaf(x, y, z)
        else if (par%xyz_symmetry) then
          select case (iface)
          case (2);  kx = -kx; il_new = il; reflected = .true.
          case (4);  ky = -ky; il_new = il; reflected = .true.
          case (6);  kz = -kz; il_new = il; reflected = .true.
          end select
        else if (par%xy_symmetry) then
          select case (iface)
          case (2);  kx = -kx; il_new = il; reflected = .true.
          case (4);  ky = -ky; il_new = il; reflected = .true.
          end select
        end if
      end if

      if (il_new <= 0) return

      Df_old    = amr_grid%Dfreq(il)
      Df_new    = amr_grid%Dfreq(il_new)
      u2        = amr_grid%vfx(il_new)*kx + amr_grid%vfy(il_new)*ky + amr_grid%vfz(il_new)*kz
      if (reflected) then
        xfreq_loc = (xfreq_loc + u1) - u2
      else
        xfreq_loc = (xfreq_loc + u1) * Df_old / Df_new - u2
      end if
      u1 = u2
      il = il_new
    end do
  end subroutine raytrace_to_edge_amr

  !=========================================================================
  ! Gas-only optical depth from photon to box edge.
  !=========================================================================
  subroutine raytrace_to_edge_tau_gas_amr(photon0, grid, tau_gas)
    use define
    implicit none
    type(photon_type), intent(in)  :: photon0
    type(grid_type),   intent(in)  :: grid
    real(wp),          intent(out) :: tau_gas

    integer  :: il, il_new, icell
    real(wp) :: x, y, z, kx, ky, kz
    real(wp) :: t_exit
    real(wp) :: rhokap, u1, u2, Df_old, Df_new, xfreq_loc
    integer  :: iface
    logical  :: reflected

    x  = photon0%x;    y  = photon0%y;    z  = photon0%z
    kx = photon0%kx;   ky = photon0%ky;   kz = photon0%kz
    il = photon0%icell_amr

    tau_gas = 0.0_wp
    if (il <= 0) then
      il = amr_find_leaf(x, y, z)
      if (il <= 0) return
    end if

    u1        = amr_grid%vfx(il)*kx + amr_grid%vfy(il)*ky + amr_grid%vfz(il)*kz
    xfreq_loc = photon0%xfreq

    do
      icell  = amr_grid%icell_of_leaf(il)
      call amr_cell_exit(icell, x, y, z, kx, ky, kz, t_exit, iface)
      rhokap  = amr_grid%rhokap(il) * voigt(xfreq_loc, amr_grid%voigt_a(il))
      tau_gas = tau_gas + t_exit * rhokap
      ! Note: no tau_huge early-out here — sightline tau output must be the
      ! full integrated value, not capped at exp(-tau)≈0.

      x = x + t_exit * kx
      y = y + t_exit * ky
      z = z + t_exit * kz
      il_new    = amr_next_leaf(icell, iface, x, y, z)
      reflected = .false.

      if (il_new <= 0) then
        if (par%xy_periodic) then
          select case (iface)
          case (1); x = x - amr_grid%xrange
          case (2); x = x + amr_grid%xrange
          case (3); y = y - amr_grid%yrange
          case (4); y = y + amr_grid%yrange
          end select
          if (iface <= 4) il_new = amr_find_leaf(x, y, z)
        else if (par%xyz_symmetry) then
          select case (iface)
          case (2);  kx = -kx; il_new = il; reflected = .true.
          case (4);  ky = -ky; il_new = il; reflected = .true.
          case (6);  kz = -kz; il_new = il; reflected = .true.
          end select
        else if (par%xy_symmetry) then
          select case (iface)
          case (2);  kx = -kx; il_new = il; reflected = .true.
          case (4);  ky = -ky; il_new = il; reflected = .true.
          end select
        end if
      end if

      if (il_new <= 0) return

      Df_old    = amr_grid%Dfreq(il)
      Df_new    = amr_grid%Dfreq(il_new)
      u2        = amr_grid%vfx(il_new)*kx + amr_grid%vfy(il_new)*ky + amr_grid%vfz(il_new)*kz
      if (reflected) then
        xfreq_loc = (xfreq_loc + u1) - u2
      else
        xfreq_loc = (xfreq_loc + u1) * Df_old / Df_new - u2
      end if
      u1 = u2
      il = il_new
    end do
  end subroutine raytrace_to_edge_tau_gas_amr

  !=========================================================================
  ! HI column density and dust optical depth from photon to box edge.
  !=========================================================================
  subroutine raytrace_to_edge_column_amr(photon0, grid, N_gas, tau_dust)
    use define
    implicit none
    type(photon_type), intent(in)  :: photon0
    type(grid_type),   intent(in)  :: grid
    real(wp),          intent(out) :: N_gas, tau_dust

    integer  :: il, il_new, icell
    real(wp) :: x, y, z, kx, ky, kz
    real(wp) :: t_exit
    real(wp) :: u1, u2, Df_old, Df_new, xfreq_loc
    integer  :: iface
    logical  :: reflected

    x  = photon0%x;    y  = photon0%y;    z  = photon0%z
    kx = photon0%kx;   ky = photon0%ky;   kz = photon0%kz
    il = photon0%icell_amr

    N_gas    = 0.0_wp
    tau_dust = 0.0_wp
    if (il <= 0) then
      il = amr_find_leaf(x, y, z)
      if (il <= 0) return
    end if

    u1        = amr_grid%vfx(il)*kx + amr_grid%vfy(il)*ky + amr_grid%vfz(il)*kz
    xfreq_loc = photon0%xfreq

    do
      icell  = amr_grid%icell_of_leaf(il)
      call amr_cell_exit(icell, x, y, z, kx, ky, kz, t_exit, iface)
      N_gas  = N_gas + t_exit * amr_grid%rhokap(il) * amr_grid%Dfreq(il) / line%cross0
      if (par%DGR > 0.0_wp .and. associated(amr_grid%rhokapD)) &
          tau_dust = tau_dust + t_exit * amr_grid%rhokapD(il)

      x = x + t_exit * kx
      y = y + t_exit * ky
      z = z + t_exit * kz
      il_new    = amr_next_leaf(icell, iface, x, y, z)
      reflected = .false.

      if (il_new <= 0) then
        if (par%xy_periodic) then
          select case (iface)
          case (1); x = x - amr_grid%xrange
          case (2); x = x + amr_grid%xrange
          case (3); y = y - amr_grid%yrange
          case (4); y = y + amr_grid%yrange
          end select
          if (iface <= 4) il_new = amr_find_leaf(x, y, z)
        else if (par%xyz_symmetry) then
          select case (iface)
          case (2);  kx = -kx; il_new = il; reflected = .true.
          case (4);  ky = -ky; il_new = il; reflected = .true.
          case (6);  kz = -kz; il_new = il; reflected = .true.
          end select
        else if (par%xy_symmetry) then
          select case (iface)
          case (2);  kx = -kx; il_new = il; reflected = .true.
          case (4);  ky = -ky; il_new = il; reflected = .true.
          end select
        end if
      end if

      if (il_new <= 0) return

      Df_old    = amr_grid%Dfreq(il)
      Df_new    = amr_grid%Dfreq(il_new)
      u2        = amr_grid%vfx(il_new)*kx + amr_grid%vfy(il_new)*ky + amr_grid%vfz(il_new)*kz
      if (reflected) then
        xfreq_loc = (xfreq_loc + u1) - u2
      else
        xfreq_loc = (xfreq_loc + u1) * Df_old / Df_new - u2
      end if
      u1 = u2
      il = il_new
    end do
  end subroutine raytrace_to_edge_column_amr

  !=========================================================================
  ! Total (gas + dust) optical depth from photon to distance dist_in.
  !=========================================================================
  subroutine raytrace_to_dist_amr(photon0, grid, dist_in, tau)
    use define
    implicit none
    type(photon_type), intent(in)  :: photon0
    type(grid_type),   intent(in)  :: grid
    real(wp),          intent(in)  :: dist_in
    real(wp),          intent(out) :: tau

    integer  :: il, il_new, icell
    real(wp) :: x, y, z, kx, ky, kz
    real(wp) :: d, t_exit
    real(wp) :: rhokap, u1, u2, Df_old, Df_new, xfreq_loc
    integer  :: iface
    logical  :: reflected

    x  = photon0%x;    y  = photon0%y;    z  = photon0%z
    kx = photon0%kx;   ky = photon0%ky;   kz = photon0%kz
    il = photon0%icell_amr

    tau = 0.0_wp
    if (il <= 0) then
      il = amr_find_leaf(x, y, z)
      if (il <= 0) return
    end if

    d         = 0.0_wp
    u1        = amr_grid%vfx(il)*kx + amr_grid%vfy(il)*ky + amr_grid%vfz(il)*kz
    xfreq_loc = photon0%xfreq

    do
      icell  = amr_grid%icell_of_leaf(il)
      call amr_cell_exit(icell, x, y, z, kx, ky, kz, t_exit, iface)
      rhokap = amr_grid%rhokap(il) * voigt(xfreq_loc, amr_grid%voigt_a(il))
      if (par%DGR > 0.0_wp .and. associated(amr_grid%rhokapD)) &
          rhokap = rhokap + amr_grid%rhokapD(il)

      if (d + t_exit >= dist_in) then
        tau = tau + (dist_in - d) * rhokap
        return
      end if
      tau = tau + t_exit * rhokap
      if (tau >= tau_huge) return

      d = d + t_exit
      x = x + t_exit * kx
      y = y + t_exit * ky
      z = z + t_exit * kz
      il_new    = amr_next_leaf(icell, iface, x, y, z)
      reflected = .false.

      if (il_new <= 0) then
        if (par%xy_periodic) then
          select case (iface)
          case (1); x = x - amr_grid%xrange
          case (2); x = x + amr_grid%xrange
          case (3); y = y - amr_grid%yrange
          case (4); y = y + amr_grid%yrange
          end select
          if (iface <= 4) il_new = amr_find_leaf(x, y, z)
        else if (par%xyz_symmetry) then
          select case (iface)
          case (2);  kx = -kx; il_new = il; reflected = .true.
          case (4);  ky = -ky; il_new = il; reflected = .true.
          case (6);  kz = -kz; il_new = il; reflected = .true.
          end select
        else if (par%xy_symmetry) then
          select case (iface)
          case (2);  kx = -kx; il_new = il; reflected = .true.
          case (4);  ky = -ky; il_new = il; reflected = .true.
          end select
        end if
      end if

      if (il_new <= 0) return

      Df_old    = amr_grid%Dfreq(il)
      Df_new    = amr_grid%Dfreq(il_new)
      u2        = amr_grid%vfx(il_new)*kx + amr_grid%vfy(il_new)*ky + amr_grid%vfz(il_new)*kz
      if (reflected) then
        xfreq_loc = (xfreq_loc + u1) - u2
      else
        xfreq_loc = (xfreq_loc + u1) * Df_old / Df_new - u2
      end if
      u1 = u2
      il = il_new
    end do
  end subroutine raytrace_to_dist_amr

  !=========================================================================
  ! Gas-only optical depth from photon to distance dist_in.
  !=========================================================================
  subroutine raytrace_to_dist_tau_gas_amr(photon0, grid, dist_in, tau_gas)
    use define
    implicit none
    type(photon_type), intent(in)  :: photon0
    type(grid_type),   intent(in)  :: grid
    real(wp),          intent(in)  :: dist_in
    real(wp),          intent(out) :: tau_gas

    integer  :: il, il_new, icell
    real(wp) :: x, y, z, kx, ky, kz
    real(wp) :: d, t_exit
    real(wp) :: rhokap, u1, u2, Df_old, Df_new, xfreq_loc
    integer  :: iface
    logical  :: reflected

    x  = photon0%x;    y  = photon0%y;    z  = photon0%z
    kx = photon0%kx;   ky = photon0%ky;   kz = photon0%kz
    il = photon0%icell_amr

    tau_gas = 0.0_wp
    if (il <= 0) then
      il = amr_find_leaf(x, y, z)
      if (il <= 0) return
    end if

    d         = 0.0_wp
    u1        = amr_grid%vfx(il)*kx + amr_grid%vfy(il)*ky + amr_grid%vfz(il)*kz
    xfreq_loc = photon0%xfreq

    do
      icell  = amr_grid%icell_of_leaf(il)
      call amr_cell_exit(icell, x, y, z, kx, ky, kz, t_exit, iface)
      rhokap = amr_grid%rhokap(il) * voigt(xfreq_loc, amr_grid%voigt_a(il))
      if (d + t_exit >= dist_in) then
        tau_gas = tau_gas + (dist_in - d) * rhokap
        return
      end if
      tau_gas = tau_gas + t_exit * rhokap

      d = d + t_exit
      x = x + t_exit * kx
      y = y + t_exit * ky
      z = z + t_exit * kz
      il_new    = amr_next_leaf(icell, iface, x, y, z)
      reflected = .false.

      if (il_new <= 0) then
        if (par%xy_periodic) then
          select case (iface)
          case (1); x = x - amr_grid%xrange
          case (2); x = x + amr_grid%xrange
          case (3); y = y - amr_grid%yrange
          case (4); y = y + amr_grid%yrange
          end select
          if (iface <= 4) il_new = amr_find_leaf(x, y, z)
        else if (par%xyz_symmetry) then
          select case (iface)
          case (2);  kx = -kx; il_new = il; reflected = .true.
          case (4);  ky = -ky; il_new = il; reflected = .true.
          case (6);  kz = -kz; il_new = il; reflected = .true.
          end select
        else if (par%xy_symmetry) then
          select case (iface)
          case (2);  kx = -kx; il_new = il; reflected = .true.
          case (4);  ky = -ky; il_new = il; reflected = .true.
          end select
        end if
      end if

      if (il_new <= 0) return

      Df_old    = amr_grid%Dfreq(il)
      Df_new    = amr_grid%Dfreq(il_new)
      u2        = amr_grid%vfx(il_new)*kx + amr_grid%vfy(il_new)*ky + amr_grid%vfz(il_new)*kz
      if (reflected) then
        xfreq_loc = (xfreq_loc + u1) - u2
      else
        xfreq_loc = (xfreq_loc + u1) * Df_old / Df_new - u2
      end if
      u1 = u2
      il = il_new
    end do
  end subroutine raytrace_to_dist_tau_gas_amr

  !=========================================================================
  ! HI column density and dust optical depth to distance dist_in.
  !=========================================================================
  subroutine raytrace_to_dist_column_amr(photon0, grid, dist_in, N_gas, tau_dust)
    use define
    implicit none
    type(photon_type), intent(in)  :: photon0
    type(grid_type),   intent(in)  :: grid
    real(wp),          intent(in)  :: dist_in
    real(wp),          intent(out) :: N_gas, tau_dust

    integer  :: il, il_new, icell
    real(wp) :: x, y, z, kx, ky, kz
    real(wp) :: d, t_exit, t_step
    real(wp) :: u1, u2, Df_old, Df_new, xfreq_loc
    integer  :: iface
    logical  :: reflected

    x  = photon0%x;    y  = photon0%y;    z  = photon0%z
    kx = photon0%kx;   ky = photon0%ky;   kz = photon0%kz
    il = photon0%icell_amr

    N_gas    = 0.0_wp
    tau_dust = 0.0_wp
    if (il <= 0) then
      il = amr_find_leaf(x, y, z)
      if (il <= 0) return
    end if

    d         = 0.0_wp
    u1        = amr_grid%vfx(il)*kx + amr_grid%vfy(il)*ky + amr_grid%vfz(il)*kz
    xfreq_loc = photon0%xfreq

    do
      icell  = amr_grid%icell_of_leaf(il)
      call amr_cell_exit(icell, x, y, z, kx, ky, kz, t_exit, iface)
      if (d + t_exit >= dist_in) then
        t_step = dist_in - d
      else
        t_step = t_exit
      end if
      N_gas    = N_gas + t_step * amr_grid%rhokap(il) * amr_grid%Dfreq(il) / line%cross0
      if (par%DGR > 0.0_wp .and. associated(amr_grid%rhokapD)) &
          tau_dust = tau_dust + t_step * amr_grid%rhokapD(il)
      if (d + t_exit >= dist_in) return

      d = d + t_exit
      x = x + t_exit * kx
      y = y + t_exit * ky
      z = z + t_exit * kz
      il_new    = amr_next_leaf(icell, iface, x, y, z)
      reflected = .false.

      if (il_new <= 0) then
        if (par%xy_periodic) then
          select case (iface)
          case (1); x = x - amr_grid%xrange
          case (2); x = x + amr_grid%xrange
          case (3); y = y - amr_grid%yrange
          case (4); y = y + amr_grid%yrange
          end select
          if (iface <= 4) il_new = amr_find_leaf(x, y, z)
        else if (par%xyz_symmetry) then
          select case (iface)
          case (2);  kx = -kx; il_new = il; reflected = .true.
          case (4);  ky = -ky; il_new = il; reflected = .true.
          case (6);  kz = -kz; il_new = il; reflected = .true.
          end select
        else if (par%xy_symmetry) then
          select case (iface)
          case (2);  kx = -kx; il_new = il; reflected = .true.
          case (4);  ky = -ky; il_new = il; reflected = .true.
          end select
        end if
      end if

      if (il_new <= 0) return

      Df_old    = amr_grid%Dfreq(il)
      Df_new    = amr_grid%Dfreq(il_new)
      u2        = amr_grid%vfx(il_new)*kx + amr_grid%vfy(il_new)*ky + amr_grid%vfz(il_new)*kz
      if (reflected) then
        xfreq_loc = (xfreq_loc + u1) - u2
      else
        xfreq_loc = (xfreq_loc + u1) * Df_old / Df_new - u2
      end if
      u1 = u2
      il = il_new
    end do
  end subroutine raytrace_to_dist_column_amr

!--- accumulate amr_grid%Jmu(xfreq, mu) at escape (mu = cos(theta_z)).
  subroutine add_to_Jmu_amr(photon, ix)
  use define
  use octree_mod, only: amr_grid
  type(photon_type), intent(in) :: photon
  integer,           intent(in) :: ix
  integer       :: imu
  real(kind=wp) :: mu_val
  if (.not. allocated(amr_grid%Jmu)) return
  mu_val = photon%kz
  if (par%xyz_symmetry) mu_val = abs(mu_val)
  imu = floor((mu_val - par%mu_min)/par%dmu) + 1
  if (imu < 1)       imu = 1
  if (imu > par%nmu) imu = par%nmu
  !$OMP ATOMIC UPDATE
  amr_grid%Jmu(ix, imu) = amr_grid%Jmu(ix, imu) + photon%wgt
  end subroutine add_to_Jmu_amr

end module raytrace_amr_mod
