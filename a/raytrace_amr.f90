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
  !   d        = cumulative distance from the ORIGINAL photon position
  !               (photon%x/y/z at entry); used only for the final
  !               photon%x = photon%x + d*kx  update.
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
    real(wp) :: tau, d, t_exit
    real(wp) :: rhokap, rhokapH, d_overshoot
    real(wp) :: u1, u2, Df_old, Df_new, xfreq_ref
    integer  :: iface

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
    d   = 0.0_wp
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

      tau = tau + t_exit * rhokap

      if (tau >= tau_in) then
        ! Overshoot: backtrack to exact tau_in
        if (rhokap > 0.0_wp) then
          d_overshoot = (tau - tau_in) / rhokap
          d           = d + t_exit - d_overshoot
        else
          d = d + t_exit
        end if
        exit
      end if

      ! Advance to face
      d = d + t_exit
      x = x + t_exit * kx
      y = y + t_exit * ky
      z = z + t_exit * kz

      ! Next leaf via precomputed neighbor table (O(1) + optional descent)
      il_new = amr_next_leaf(icell, iface, x, y, z)

      if (il_new <= 0) then
        photon%inside = .false.
        exit
      end if

      ! Frequency shift across cell boundary
      Df_old = amr_grid%Dfreq(il)
      Df_new = amr_grid%Dfreq(il_new)
      u2     = amr_grid%vfx(il_new)*kx + amr_grid%vfy(il_new)*ky + amr_grid%vfz(il_new)*kz
      photon%xfreq = (photon%xfreq + u1) * Df_old / Df_new - u2
      u1 = u2

      il = il_new
    end do

    ! Final position
    photon%x = photon%x + d * kx
    photon%y = photon%y + d * ky
    photon%z = photon%z + d * kz
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
      il_new = amr_next_leaf(icell, iface, x, y, z)
      if (il_new <= 0) return

      Df_old    = amr_grid%Dfreq(il)
      Df_new    = amr_grid%Dfreq(il_new)
      u2        = amr_grid%vfx(il_new)*kx + amr_grid%vfy(il_new)*ky + amr_grid%vfz(il_new)*kz
      xfreq_loc = (xfreq_loc + u1) * Df_old / Df_new - u2
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
      if (tau_gas >= tau_huge) return

      x = x + t_exit * kx
      y = y + t_exit * ky
      z = z + t_exit * kz
      il_new = amr_next_leaf(icell, iface, x, y, z)
      if (il_new <= 0) return

      Df_old    = amr_grid%Dfreq(il)
      Df_new    = amr_grid%Dfreq(il_new)
      u2        = amr_grid%vfx(il_new)*kx + amr_grid%vfy(il_new)*ky + amr_grid%vfz(il_new)*kz
      xfreq_loc = (xfreq_loc + u1) * Df_old / Df_new - u2
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
      il_new = amr_next_leaf(icell, iface, x, y, z)
      if (il_new <= 0) return

      Df_old    = amr_grid%Dfreq(il)
      Df_new    = amr_grid%Dfreq(il_new)
      u2        = amr_grid%vfx(il_new)*kx + amr_grid%vfy(il_new)*ky + amr_grid%vfz(il_new)*kz
      xfreq_loc = (xfreq_loc + u1) * Df_old / Df_new - u2
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
      il_new = amr_next_leaf(icell, iface, x, y, z)
      if (il_new <= 0) return

      Df_old    = amr_grid%Dfreq(il)
      Df_new    = amr_grid%Dfreq(il_new)
      u2        = amr_grid%vfx(il_new)*kx + amr_grid%vfy(il_new)*ky + amr_grid%vfz(il_new)*kz
      xfreq_loc = (xfreq_loc + u1) * Df_old / Df_new - u2
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
      il_new = amr_next_leaf(icell, iface, x, y, z)
      if (il_new <= 0) return

      Df_old    = amr_grid%Dfreq(il)
      Df_new    = amr_grid%Dfreq(il_new)
      u2        = amr_grid%vfx(il_new)*kx + amr_grid%vfy(il_new)*ky + amr_grid%vfz(il_new)*kz
      xfreq_loc = (xfreq_loc + u1) * Df_old / Df_new - u2
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
      il_new = amr_next_leaf(icell, iface, x, y, z)
      if (il_new <= 0) return

      Df_old    = amr_grid%Dfreq(il)
      Df_new    = amr_grid%Dfreq(il_new)
      u2        = amr_grid%vfx(il_new)*kx + amr_grid%vfy(il_new)*ky + amr_grid%vfz(il_new)*kz
      xfreq_loc = (xfreq_loc + u1) * Df_old / Df_new - u2
      u1 = u2
      il = il_new
    end do
  end subroutine raytrace_to_dist_column_amr

end module raytrace_amr_mod
