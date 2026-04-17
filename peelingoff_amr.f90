module peelingoff_amr_mod
  !-----------------------------------------------------------------------
  ! AMR peel-off routines: equivalents of peelingoff_rect.f90 routines
  ! that use amr_grid leaf-cell data instead of Cartesian grid arrays.
  !
  ! Key differences from Cartesian versions:
  !   grid%vfx(icell,jcell,kcell)   → amr_grid%vfx(photon%icell_amr)
  !   grid%Dfreq(icell,jcell,kcell) → amr_grid%Dfreq(photon%icell_amr)
  !   call raytrace_to_edge(...)     → procedure pointer → raytrace_to_edge_amr
  !
  ! grid%Dfreq_ref, grid%xfreq_min, grid%dxfreq, grid%nxfreq are
  ! filled from amr_grid by amr_sync_to_grid, so grid% lookups are fine.
  !-----------------------------------------------------------------------
  use define
  use mathlib
  use utility
  use line_mod
  use octree_mod
  implicit none
  private

  public :: peeling_direct_amr
  public :: peeling_dust_nostokes_amr
  public :: peeling_resonance_nostokes_amr

contains

  !=========================================================================
  subroutine peeling_direct_amr(photon, grid)
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  type(photon_type) :: pobs
  real(wp) :: r2, r, wgt, wgt0, tau
  real(wp) :: kx, ky, kz
  real(wp) :: xfreq_ref, u1, u2
  integer  :: ix, iy, ixf, il, i

  il = photon%icell_amr

  do i = 1, par%nobs
    pobs    = photon
    pobs%kx = observer(i)%x - photon%x
    pobs%ky = observer(i)%y - photon%y
    pobs%kz = observer(i)%z - photon%z
    r2      = pobs%kx**2 + pobs%ky**2 + pobs%kz**2
    r       = sqrt(r2)
    pobs%kx = pobs%kx / r
    pobs%ky = pobs%ky / r
    pobs%kz = pobs%kz / r

    kx = observer(i)%rmatrix(1,1)*pobs%kx + observer(i)%rmatrix(1,2)*pobs%ky + observer(i)%rmatrix(1,3)*pobs%kz
    ky = observer(i)%rmatrix(2,1)*pobs%kx + observer(i)%rmatrix(2,2)*pobs%ky + observer(i)%rmatrix(2,3)*pobs%kz
    kz = observer(i)%rmatrix(3,1)*pobs%kx + observer(i)%rmatrix(3,2)*pobs%ky + observer(i)%rmatrix(3,3)*pobs%kz

    ix = floor(atan2(-kx,kz)*rad2deg/observer(i)%dxim + observer(i)%nxim/2.0_wp) + 1
    iy = floor(atan2(-ky,kz)*rad2deg/observer(i)%dyim + observer(i)%nyim/2.0_wp) + 1

    if (.not. par%comoving_source) then
       u1        = amr_grid%vfx(il)*photon%kx + amr_grid%vfy(il)*photon%ky + amr_grid%vfz(il)*photon%kz
       xfreq_ref = photon%xfreq + u1
       u2        = amr_grid%vfx(il)*pobs%kx   + amr_grid%vfy(il)*pobs%ky   + amr_grid%vfz(il)*pobs%kz
       pobs%xfreq = xfreq_ref - u2
    else
       u1        = amr_grid%vfx(il)*pobs%kx + amr_grid%vfy(il)*pobs%ky + amr_grid%vfz(il)*pobs%kz
       xfreq_ref = photon%xfreq + u1
    end if

    xfreq_ref = xfreq_ref * (amr_grid%Dfreq(il) / grid%Dfreq_ref)
    ixf       = floor((xfreq_ref - grid%xfreq_min) / grid%dxfreq) + 1

    if (ix >= 1 .and. ix <= observer(i)%nxim .and. iy >= 1 .and. iy <= observer(i)%nyim) then
       call raytrace_to_edge(pobs, grid, tau)
       wgt0 = 1.0_wp / (fourpi*r2) * photon%wgt
       wgt  = exp(-tau) * wgt0

       if (par%save_peeloff_2D) then
          !$OMP ATOMIC UPDATE
          observer(i)%direc_2D(ix,iy) = observer(i)%direc_2D(ix,iy) + wgt
          if (par%save_direc0) then
             !$OMP ATOMIC UPDATE
             observer(i)%direc0_2D(ix,iy) = observer(i)%direc0_2D(ix,iy) + wgt0
          end if
       end if

       if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
          !$OMP ATOMIC UPDATE
          observer(i)%direc(ixf,ix,iy) = observer(i)%direc(ixf,ix,iy) + wgt
          if (par%save_direc0) then
             !$OMP ATOMIC UPDATE
             observer(i)%direc0(ixf,ix,iy) = observer(i)%direc0(ixf,ix,iy) + wgt0
          end if
       end if
    end if
  end do
  end subroutine peeling_direct_amr

  !=========================================================================
  subroutine peeling_dust_nostokes_amr(photon, grid)
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  type(photon_type) :: pobs
  real(wp) :: r2, r, kx, ky, kz
  real(wp) :: cosa, wgt, peel, tau
  real(wp) :: xfreq_ref, u1
  integer  :: ix, iy, ixf, il, i

  il = photon%icell_amr

  do i = 1, par%nobs
    pobs    = photon
    pobs%kx = observer(i)%x - photon%x
    pobs%ky = observer(i)%y - photon%y
    pobs%kz = observer(i)%z - photon%z
    r2      = pobs%kx**2 + pobs%ky**2 + pobs%kz**2
    r       = sqrt(r2)
    pobs%kx = pobs%kx / r
    pobs%ky = pobs%ky / r
    pobs%kz = pobs%kz / r

    kx = observer(i)%rmatrix(1,1)*pobs%kx + observer(i)%rmatrix(1,2)*pobs%ky + observer(i)%rmatrix(1,3)*pobs%kz
    ky = observer(i)%rmatrix(2,1)*pobs%kx + observer(i)%rmatrix(2,2)*pobs%ky + observer(i)%rmatrix(2,3)*pobs%kz
    kz = observer(i)%rmatrix(3,1)*pobs%kx + observer(i)%rmatrix(3,2)*pobs%ky + observer(i)%rmatrix(3,3)*pobs%kz

    ix = floor(atan2(-kx,kz)*rad2deg/observer(i)%dxim + observer(i)%nxim/2.0_wp) + 1
    iy = floor(atan2(-ky,kz)*rad2deg/observer(i)%dyim + observer(i)%nyim/2.0_wp) + 1

    if (ix >= 1 .and. ix <= observer(i)%nxim .and. iy >= 1 .and. iy <= observer(i)%nyim) then
       u1        = amr_grid%vfx(il)*pobs%kx + amr_grid%vfy(il)*pobs%ky + amr_grid%vfz(il)*pobs%kz
       xfreq_ref = photon%xfreq + u1
       xfreq_ref = xfreq_ref * (amr_grid%Dfreq(il) / grid%Dfreq_ref)
       ixf       = floor((xfreq_ref - grid%xfreq_min) / grid%dxfreq) + 1

       call raytrace_to_edge(pobs, grid, tau)
       cosa = photon%kx*pobs%kx + photon%ky*pobs%ky + photon%kz*pobs%kz
       peel = (1.0_wp - par%hgg**2) / ((1.0_wp + par%hgg**2) - 2.0_wp*par%hgg*cosa)**1.5_wp / fourpi
       wgt  = peel / r2 * exp(-tau) * photon%wgt

       if (par%save_peeloff_2D) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt_2D(ix,iy) = observer(i)%scatt_2D(ix,iy) + wgt
       end if

       if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt(ixf,ix,iy) = observer(i)%scatt(ixf,ix,iy) + wgt
          if (par%save_dust_scattered) then
             !$OMP ATOMIC UPDATE
             observer(i)%scatt_dust(ixf,ix,iy) = observer(i)%scatt_dust(ixf,ix,iy) + wgt
          end if
       end if
    end if
  end do
  end subroutine peeling_dust_nostokes_amr

  !=========================================================================
  subroutine peeling_resonance_nostokes_amr(photon, grid, xfreq_atom, vel_atom)
  use random
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  real(wp),          intent(in) :: xfreq_atom
  real(wp),          intent(in) :: vel_atom(3)
  type(photon_type) :: pobs
  real(wp) :: r2, r, kx, ky, kz
  real(wp) :: wgt, peel, tau
  real(wp) :: cost, cost2, sint, cosp, sinp, rho1, rho
  real(wp) :: xfreq_ref, xfreq, g_recoil, u1
  integer  :: ix, iy, ixf, il, i

  il = photon%icell_amr

  do i = 1, par%nobs
    pobs    = photon
    pobs%kx = observer(i)%x - photon%x
    pobs%ky = observer(i)%y - photon%y
    pobs%kz = observer(i)%z - photon%z
    r2      = pobs%kx**2 + pobs%ky**2 + pobs%kz**2
    r       = sqrt(r2)
    pobs%kx = pobs%kx / r
    pobs%ky = pobs%ky / r
    pobs%kz = pobs%kz / r

    kx = observer(i)%rmatrix(1,1)*pobs%kx + observer(i)%rmatrix(1,2)*pobs%ky + observer(i)%rmatrix(1,3)*pobs%kz
    ky = observer(i)%rmatrix(2,1)*pobs%kx + observer(i)%rmatrix(2,2)*pobs%ky + observer(i)%rmatrix(2,3)*pobs%kz
    kz = observer(i)%rmatrix(3,1)*pobs%kx + observer(i)%rmatrix(3,2)*pobs%ky + observer(i)%rmatrix(3,3)*pobs%kz

    ix = floor(atan2(-kx,kz)*rad2deg/observer(i)%dxim + observer(i)%nxim/2.0_wp) + 1
    iy = floor(atan2(-ky,kz)*rad2deg/observer(i)%dyim + observer(i)%nyim/2.0_wp) + 1

    if (ix >= 1 .and. ix <= observer(i)%nxim .and. iy >= 1 .and. iy <= observer(i)%nyim) then
       cost  = photon%kx*pobs%kx + photon%ky*pobs%ky + photon%kz*pobs%kz
       cost2 = cost**2
       sint  = sqrt(1.0_wp - cost2)
       rho1  = sqrt(1.0_wp - photon%kz**2) * sint

       if (rho1 == 0.0_wp) then
          cosp  = 1.0_wp
          sinp  = 0.0_wp
          xfreq = xfreq_atom + (vel_atom(1)*cosp + vel_atom(2)*sinp)*sint + vel_atom(3)*cost
       else
          rho   = 1.0_wp / rho1
          cosp  = rho * (cost*photon%kz - pobs%kz)
          sinp  = rho * (photon%kx*pobs%ky - pobs%kx*photon%ky)
          xfreq = xfreq_atom + (vel_atom(1)*cosp + vel_atom(2)*sinp)*sint + vel_atom(3)*cost
       end if

       if (par%recoil) then
          g_recoil = line%g_recoil0 / amr_grid%Dfreq(il)
          xfreq    = xfreq - g_recoil * (1.0_wp - cost)
       end if

       u1        = amr_grid%vfx(il)*pobs%kx + amr_grid%vfy(il)*pobs%ky + amr_grid%vfz(il)*pobs%kz
       xfreq_ref = xfreq + u1
       xfreq_ref = xfreq_ref * (amr_grid%Dfreq(il) / grid%Dfreq_ref)
       ixf       = floor((xfreq_ref - grid%xfreq_min) / grid%dxfreq) + 1

       pobs%xfreq = xfreq
       call raytrace_to_edge(pobs, grid, tau)
       peel = three_over_four * photon%E1 * (cost2 + 1.0_wp) + photon%E2
       wgt  = peel / (fourpi*r2) * exp(-tau) * photon%wgt

       if (par%save_peeloff_2D) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt_2D(ix,iy) = observer(i)%scatt_2D(ix,iy) + wgt
       end if

       if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt(ixf,ix,iy) = observer(i)%scatt(ixf,ix,iy) + wgt
       end if
    end if
  end do
  end subroutine peeling_resonance_nostokes_amr

end module peelingoff_amr_mod
