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
  public :: peeling_dust_stokes_amr
  public :: peeling_resonance_nostokes_amr
  public :: peeling_resonance_stokes_amr
  public :: peeling_direct_inside_amr
  public :: peeling_dust_nostokes_inside_amr
  public :: peeling_resonance_nostokes_inside_amr

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

  !=========================================================================
  subroutine peeling_dust_stokes_amr(photon, grid)
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  type(photon_type) :: pobs
  real(wp) :: r2, r
  real(wp) :: cost, sint, cosp, sinp, cos2p, sin2p, cosg, sing, cos2g, sin2g
  real(wp) :: S11, S12, S33, S34
  real(wp) :: Q0, U0
  real(wp) :: Iobs, Qobs, Uobs, Vobs, Idet, Qdet, Udet, Vdet
  real(wp) :: kx, ky, kz
  real(wp) :: tau, wgt
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

    u1        = amr_grid%vfx(il)*pobs%kx + amr_grid%vfy(il)*pobs%ky + amr_grid%vfz(il)*pobs%kz
    xfreq_ref = photon%xfreq + u1
    xfreq_ref = xfreq_ref * (amr_grid%Dfreq(il) / grid%Dfreq_ref)
    ixf       = floor((xfreq_ref - grid%xfreq_min) / grid%dxfreq) + 1

    if (ix >= 1 .and. ix <= observer(i)%nxim .and. iy >= 1 .and. iy <= observer(i)%nyim) then
       cost = photon%kx*pobs%kx + photon%ky*pobs%ky + photon%kz*pobs%kz
       sint = sqrt(1.0_wp - cost*cost)

       if (sint == 0.0_wp) then
          cosp  = 1.0_wp;  sinp  = 0.0_wp
          cos2p = 1.0_wp;  sin2p = 0.0_wp
       else
          cosp  = (pobs%kx*photon%mx + pobs%ky*photon%my + pobs%kz*photon%mz) / sint
          sinp  = (pobs%kx*photon%nx + pobs%ky*photon%ny + pobs%kz*photon%nz) / sint
          cos2p = 2.0_wp*cosp*cosp - 1.0_wp
          sin2p = 2.0_wp*cosp*sinp
       end if

       pobs%nx = -sinp*photon%mx + cosp*photon%nx
       pobs%ny = -sinp*photon%my + cosp*photon%ny
       pobs%nz = -sinp*photon%mz + cosp*photon%nz

       call interp_eq(scatt_mat%coss, scatt_mat%S11, cost, S11)
       call interp_eq(scatt_mat%coss, scatt_mat%S12, cost, S12)
       call interp_eq(scatt_mat%coss, scatt_mat%S33, cost, S33)
       call interp_eq(scatt_mat%coss, scatt_mat%S34, cost, S34)

       Q0 =  cos2p*photon%Q + sin2p*photon%U
       U0 = -sin2p*photon%Q + cos2p*photon%U

       Iobs = ( S11*photon%I + S12*Q0      ) / twopi
       Qobs = ( S12*photon%I + S11*Q0      ) / twopi
       Uobs = ( S33*U0       + S34*photon%V) / twopi
       Vobs = (-S34*U0       + S33*photon%V) / twopi

       cosg  = -(observer(i)%rmatrix(1,1)*pobs%nx + observer(i)%rmatrix(1,2)*pobs%ny + observer(i)%rmatrix(1,3)*pobs%nz)
       sing  =   observer(i)%rmatrix(2,1)*pobs%nx + observer(i)%rmatrix(2,2)*pobs%ny + observer(i)%rmatrix(2,3)*pobs%nz
       cos2g = 2.0_wp*cosg*cosg - 1.0_wp
       sin2g = 2.0_wp*cosg*sing

       Idet =  Iobs
       Qdet =  cos2g*Qobs + sin2g*Uobs
       Udet = -sin2g*Qobs + cos2g*Uobs
       Vdet =  Vobs

       call raytrace_to_edge(pobs, grid, tau)
       wgt = 1.0_wp / r2 * exp(-tau) * photon%wgt

       if (par%save_peeloff_2D) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt_2D(ix,iy) = observer(i)%scatt_2D(ix,iy) + wgt * Idet
          !$OMP ATOMIC UPDATE
          observer(i)%I_2D(ix,iy)     = observer(i)%I_2D(ix,iy)     + wgt * Idet
          !$OMP ATOMIC UPDATE
          observer(i)%Q_2D(ix,iy)     = observer(i)%Q_2D(ix,iy)     + wgt * Qdet
          !$OMP ATOMIC UPDATE
          observer(i)%U_2D(ix,iy)     = observer(i)%U_2D(ix,iy)     + wgt * Udet
          !$OMP ATOMIC UPDATE
          observer(i)%V_2D(ix,iy)     = observer(i)%V_2D(ix,iy)     + wgt * Vdet
       end if

       if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt(ixf,ix,iy) = observer(i)%scatt(ixf,ix,iy) + wgt * Idet
          if (par%save_dust_scattered) then
             !$OMP ATOMIC UPDATE
             observer(i)%scatt_dust(ixf,ix,iy) = observer(i)%scatt_dust(ixf,ix,iy) + wgt * Idet
          end if
          !$OMP ATOMIC UPDATE
          observer(i)%I(ixf,ix,iy) = observer(i)%I(ixf,ix,iy) + wgt * Idet
          !$OMP ATOMIC UPDATE
          observer(i)%Q(ixf,ix,iy) = observer(i)%Q(ixf,ix,iy) + wgt * Qdet
          !$OMP ATOMIC UPDATE
          observer(i)%U(ixf,ix,iy) = observer(i)%U(ixf,ix,iy) + wgt * Udet
          !$OMP ATOMIC UPDATE
          observer(i)%V(ixf,ix,iy) = observer(i)%V(ixf,ix,iy) + wgt * Vdet
       end if
    end if
  end do
  end subroutine peeling_dust_stokes_amr

  !=========================================================================
  subroutine peeling_resonance_stokes_amr(photon, grid, xfreq_atom, vel_atom)
  use random
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  real(wp),          intent(in) :: xfreq_atom
  real(wp),          intent(in) :: vel_atom(3)
  type(photon_type) :: pobs
  real(wp) :: r2, r
  real(wp) :: cost, cost2, sint, cosp, sinp, cos2p, sin2p, cosg, sing, cos2g, sin2g
  real(wp) :: S11, S12, S22, S33, S44
  real(wp) :: Q0, U0
  real(wp) :: Iobs, Qobs, Uobs, Vobs, Idet, Qdet, Udet, Vdet
  real(wp) :: kx, ky, kz
  real(wp) :: tau, wgt
  real(wp) :: xfreq, xfreq_ref, g_recoil, u1
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

       if (sint == 0.0_wp) then
          cosp  = 1.0_wp;  sinp  = 0.0_wp
          cos2p = 1.0_wp;  sin2p = 0.0_wp
       else
          cosp  = (pobs%kx*photon%mx + pobs%ky*photon%my + pobs%kz*photon%mz) / sint
          sinp  = (pobs%kx*photon%nx + pobs%ky*photon%ny + pobs%kz*photon%nz) / sint
          cos2p = 2.0_wp*cosp*cosp - 1.0_wp
          sin2p = 2.0_wp*cosp*sinp
       end if

       pobs%nx = -sinp*photon%mx + cosp*photon%nx
       pobs%ny = -sinp*photon%my + cosp*photon%ny
       pobs%nz = -sinp*photon%mz + cosp*photon%nz

       S22 = three_over_four * photon%E1 * (cost2 + 1.0_wp)
       S11 = S22 + photon%E2
       S12 = three_over_four * photon%E1 * (cost2 - 1.0_wp)
       S33 = three_over_two  * photon%E1 * cost
       S44 = three_over_two  * photon%E3 * cost

       xfreq = xfreq_atom + (vel_atom(1)*cosp + vel_atom(2)*sinp)*sint + vel_atom(3)*cost
       if (par%recoil) then
          g_recoil = line%g_recoil0 / amr_grid%Dfreq(il)
          xfreq    = xfreq - g_recoil * (1.0_wp - cost)
       end if

       u1        = amr_grid%vfx(il)*pobs%kx + amr_grid%vfy(il)*pobs%ky + amr_grid%vfz(il)*pobs%kz
       xfreq_ref = xfreq + u1
       xfreq_ref = xfreq_ref * (amr_grid%Dfreq(il) / grid%Dfreq_ref)
       ixf       = floor((xfreq_ref - grid%xfreq_min) / grid%dxfreq) + 1

       Q0 =  cos2p*photon%Q + sin2p*photon%U
       U0 = -sin2p*photon%Q + cos2p*photon%U

       Iobs = (S11 + S12*Q0) / fourpi
       Qobs = (S12 + S22*Q0) / fourpi
       Uobs = (S33*U0      ) / fourpi
       Vobs = (S44*photon%V) / fourpi

       cosg  = -(observer(i)%rmatrix(1,1)*pobs%nx + observer(i)%rmatrix(1,2)*pobs%ny + observer(i)%rmatrix(1,3)*pobs%nz)
       sing  =   observer(i)%rmatrix(2,1)*pobs%nx + observer(i)%rmatrix(2,2)*pobs%ny + observer(i)%rmatrix(2,3)*pobs%nz
       cos2g = 2.0_wp*cosg*cosg - 1.0_wp
       sin2g = 2.0_wp*cosg*sing

       Idet =  Iobs
       Qdet =  cos2g*Qobs + sin2g*Uobs
       Udet = -sin2g*Qobs + cos2g*Uobs
       Vdet =  Vobs

       pobs%xfreq = xfreq
       call raytrace_to_edge(pobs, grid, tau)
       wgt = 1.0_wp / r2 * exp(-tau) * photon%wgt

       if (par%save_peeloff_2D) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt_2D(ix,iy) = observer(i)%scatt_2D(ix,iy) + wgt * Idet
          !$OMP ATOMIC UPDATE
          observer(i)%I_2D(ix,iy)     = observer(i)%I_2D(ix,iy)     + wgt * Idet
          !$OMP ATOMIC UPDATE
          observer(i)%Q_2D(ix,iy)     = observer(i)%Q_2D(ix,iy)     + wgt * Qdet
          !$OMP ATOMIC UPDATE
          observer(i)%U_2D(ix,iy)     = observer(i)%U_2D(ix,iy)     + wgt * Udet
          !$OMP ATOMIC UPDATE
          observer(i)%V_2D(ix,iy)     = observer(i)%V_2D(ix,iy)     + wgt * Vdet
       end if

       if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt(ixf,ix,iy) = observer(i)%scatt(ixf,ix,iy) + wgt * Idet
          !$OMP ATOMIC UPDATE
          observer(i)%I(ixf,ix,iy) = observer(i)%I(ixf,ix,iy) + wgt * Idet
          !$OMP ATOMIC UPDATE
          observer(i)%Q(ixf,ix,iy) = observer(i)%Q(ixf,ix,iy) + wgt * Qdet
          !$OMP ATOMIC UPDATE
          observer(i)%U(ixf,ix,iy) = observer(i)%U(ixf,ix,iy) + wgt * Udet
          !$OMP ATOMIC UPDATE
          observer(i)%V(ixf,ix,iy) = observer(i)%V(ixf,ix,iy) + wgt * Vdet
       end if
    end if
  end do
  end subroutine peeling_resonance_stokes_amr

  !=========================================================================
  ! HEALPix inside-observer AMR variants.
  ! These mirror peelingoff_heal.f90's inside routines but use amr_grid leaf
  ! arrays via photon%icell_amr.  Stokes is intentionally not supported (the
  ! Cartesian _inside path does not support Stokes either).
  ! raytrace_to_dist is a procedure pointer already routed to raytrace_to_dist_amr
  ! when par%use_amr_grid is set, so those calls remain unchanged.
  !=========================================================================
  subroutine peeling_direct_inside_amr(photon,grid)
  use healpix
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  type(photon_type) :: pobs
  real(wp) :: r2, r, wgt, wgt0, tau
  real(wp) :: xfreq_ref, u1, u2
  integer  :: ipix, ixf, il, i

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

    call vec2pix(observer(i)%nside, -pobs%kx, -pobs%ky, -pobs%kz, ipix)

    if (.not. par%comoving_source) then
       u1 = amr_grid%vfx(il)*photon%kx + amr_grid%vfy(il)*photon%ky + amr_grid%vfz(il)*photon%kz
       xfreq_ref = photon%xfreq + u1
       u2 = amr_grid%vfx(il)*pobs%kx + amr_grid%vfy(il)*pobs%ky + amr_grid%vfz(il)*pobs%kz
       pobs%xfreq = xfreq_ref - u2
    else
       u1 = amr_grid%vfx(il)*pobs%kx + amr_grid%vfy(il)*pobs%ky + amr_grid%vfz(il)*pobs%kz
       xfreq_ref = photon%xfreq + u1
    endif

    xfreq_ref = xfreq_ref * (amr_grid%Dfreq(il) / grid%Dfreq_ref)
    ixf       = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq) + 1

    if (ipix >= 1 .and. ipix <= observer(i)%npix) then
       call raytrace_to_dist(pobs,grid,r,tau)
       wgt0 = 1.0_wp/(fourpi*r2) * photon%wgt
       wgt  = exp(-tau)*wgt0

       if (par%save_peeloff_2D) then
          !$OMP ATOMIC UPDATE
          observer(i)%direc_heal_2D(ipix) = observer(i)%direc_heal_2D(ipix) + wgt
          if (par%save_direc0) then
             !$OMP ATOMIC UPDATE
             observer(i)%direc0_heal_2D(ipix) = observer(i)%direc0_heal_2D(ipix) + wgt0
          endif
       endif

       if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
          !$OMP ATOMIC UPDATE
          observer(i)%direc_heal(ixf,ipix) = observer(i)%direc_heal(ixf,ipix) + wgt
          if (par%save_direc0) then
             !$OMP ATOMIC UPDATE
             observer(i)%direc0_heal(ixf,ipix) = observer(i)%direc0_heal(ixf,ipix) + wgt0
          endif
       endif
    endif
  end do
  end subroutine peeling_direct_inside_amr

  !=========================================================================
  subroutine peeling_dust_nostokes_inside_amr(photon,grid)
  use healpix
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  type(photon_type) :: pobs
  real(wp) :: r2, r, cosa, wgt, peel, tau
  real(wp) :: xfreq_ref, u1
  integer  :: ipix, ixf, il, i

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

    call vec2pix(observer(i)%nside, -pobs%kx, -pobs%ky, -pobs%kz, ipix)

    if (ipix >= 1 .and. ipix <= observer(i)%npix) then
       u1 = amr_grid%vfx(il)*pobs%kx + amr_grid%vfy(il)*pobs%ky + amr_grid%vfz(il)*pobs%kz
       xfreq_ref = photon%xfreq + u1

       xfreq_ref = xfreq_ref * (amr_grid%Dfreq(il) / grid%Dfreq_ref)
       ixf       = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq) + 1

       call raytrace_to_dist(pobs,grid,r,tau)
       cosa = photon%kx*pobs%kx + photon%ky*pobs%ky + photon%kz*pobs%kz
       peel = (1.0_wp - par%hgg**2)/((1.0_wp + par%hgg**2) - 2.0_wp*par%hgg*cosa)**1.5_wp / fourpi
       wgt  = peel/r2 * exp(-tau) * photon%wgt

       if (par%save_peeloff_2D) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt_heal_2D(ipix) = observer(i)%scatt_heal_2D(ipix) + wgt
       endif
       if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt_heal(ixf,ipix) = observer(i)%scatt_heal(ixf,ipix) + wgt
       endif
    endif
  end do
  end subroutine peeling_dust_nostokes_inside_amr

  !=========================================================================
  subroutine peeling_resonance_nostokes_inside_amr(photon,grid,xfreq_atom,vel_atom)
  use healpix
  use random
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  real(wp),          intent(in) :: xfreq_atom
  real(wp),          intent(in) :: vel_atom(3)
  type(photon_type) :: pobs
  real(wp) :: r2, r, wgt, peel, tau
  real(wp) :: cost, cost2, sint, cosp, sinp, rho1, rho
  real(wp) :: xfreq_ref, xfreq, g_recoil, u1
  integer  :: ipix, ixf, il, i

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

    call vec2pix(observer(i)%nside, -pobs%kx, -pobs%ky, -pobs%kz, ipix)

    if (ipix >= 1 .and. ipix <= observer(i)%npix) then
       cost  = photon%kx*pobs%kx + photon%ky*pobs%ky + photon%kz*pobs%kz
       cost2 = cost**2
       sint  = sqrt(1.0_wp - cost2)
       rho1  = sqrt(1.0_wp - photon%kz**2) * sint

       if (rho1 == 0.0_wp) then
          cosp  = 1.0_wp
          sinp  = 0.0_wp
          xfreq = xfreq_atom + (vel_atom(1)*cosp + vel_atom(2)*sinp)*sint + vel_atom(3)*cost
       else
          rho   = 1.0_wp/rho1
          cosp  = rho * (cost*photon%kz - pobs%kz)
          sinp  = rho * (photon%kx*pobs%ky - pobs%kx*photon%ky)
          xfreq = xfreq_atom + (vel_atom(1)*cosp + vel_atom(2)*sinp)*sint + vel_atom(3)*cost
       endif

       if (par%recoil) then
          g_recoil = line%g_recoil0 / amr_grid%Dfreq(il)
          xfreq    = xfreq - g_recoil * (1.0_wp - cost)
       endif

       u1 = amr_grid%vfx(il)*pobs%kx + amr_grid%vfy(il)*pobs%ky + amr_grid%vfz(il)*pobs%kz
       xfreq_ref = xfreq + u1

       xfreq_ref = xfreq_ref * (amr_grid%Dfreq(il) / grid%Dfreq_ref)
       ixf       = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq) + 1

       pobs%xfreq = xfreq
       call raytrace_to_dist(pobs,grid,r,tau)
       peel = three_over_four * photon%E1 * (cost2 + 1.0_wp) + photon%E2
       wgt  = peel/(fourpi*r2) * exp(-tau) * photon%wgt

       if (par%save_peeloff_2D) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt_heal_2D(ipix) = observer(i)%scatt_heal_2D(ipix) + wgt
       endif
       if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt_heal(ixf,ipix) = observer(i)%scatt_heal(ixf,ipix) + wgt
       endif
    endif
  end do
  end subroutine peeling_resonance_nostokes_inside_amr

end module peelingoff_amr_mod
