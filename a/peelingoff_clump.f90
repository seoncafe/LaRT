!-- Peeling-off routines for the clumpy medium model.
!--
!-- The Cartesian routines (peelingoff_rect) compute the lab-frame frequency for
!-- spectral binning using grid%vfx/vfy/vfz, which are identically 0 in clump
!-- mode (bulk velocity is stored per-clump in cl_vx/cl_vy/cl_vz).  These
!-- replacements apply the correct Doppler correction from the scattering
!-- clump's bulk velocity.
!--
!-- Convention:
!--   photon%icell_clump > 0  → photon is inside that clump; xfreq is in the
!--                              clump's rest frame.
!--   photon%icell_clump = 0  → photon is in vacuum; xfreq is in the lab frame.
!--
module peelingoff_clump_mod
  use define
  use mathlib
  use utility
  use line_mod
  use clump_mod
contains

  !---------------------------------------------------------------------------
  ! Lab-frame frequency correction at a scattering point inside clump icl.
  ! xfreq     : frequency in clump rest frame (x-units = Doppler widths)
  ! kx,ky,kz  : unit vector toward observer (peel-off direction)
  ! icl       : clump index (0 = vacuum, no correction)
  ! Returns frequency in the lab frame (still in x-units).
  !---------------------------------------------------------------------------
  real(wp) function peel_xfreq_lab(xfreq, kx, ky, kz, icl)
  implicit none
  real(wp), intent(in) :: xfreq, kx, ky, kz
  integer,  intent(in) :: icl
  real(wp) :: u_los
  if (icl > 0) then
     u_los = real((cl_vx(icl)*kx + cl_vy(icl)*ky + cl_vz(icl)*kz) / cl_vtherm, wp)
     peel_xfreq_lab = xfreq + u_los
  else
     peel_xfreq_lab = xfreq
  end if
  end function peel_xfreq_lab

  !--------------------------------------------------
  subroutine peeling_direct_clump(photon, grid)
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  ! local variables
  type(photon_type) :: pobs
  real(wp) :: r2, r, wgt, wgt0, tau
  real(wp) :: kx, ky, kz
  real(wp) :: xfreq_ref, u_los
  integer  :: ix, iy, ixf
  integer  :: i, icl

  do i = 1, par%nobs
    if ((trim(par%source_geometry) == 'plane_illumination') .and. &
        (observer(i)%x /= 0.0_wp .or. observer(i)%y /= 0.0_wp)) cycle

    pobs    = photon
    pobs%kx = observer(i)%x - photon%x
    pobs%ky = observer(i)%y - photon%y
    pobs%kz = observer(i)%z - photon%z
    r2      = pobs%kx**2 + pobs%ky**2 + pobs%kz**2
    r       = sqrt(r2)
    pobs%kx = pobs%kx / r
    pobs%ky = pobs%ky / r
    pobs%kz = pobs%kz / r

    !--- Transform the peeling-off vector to the observer's frame.
    kx = observer(i)%rmatrix(1,1)*pobs%kx + observer(i)%rmatrix(1,2)*pobs%ky + observer(i)%rmatrix(1,3)*pobs%kz
    ky = observer(i)%rmatrix(2,1)*pobs%kx + observer(i)%rmatrix(2,2)*pobs%ky + observer(i)%rmatrix(2,3)*pobs%kz
    kz = observer(i)%rmatrix(3,1)*pobs%kx + observer(i)%rmatrix(3,2)*pobs%ky + observer(i)%rmatrix(3,3)*pobs%kz

    ix = floor(atan2(-kx,kz)*rad2deg/observer(i)%dxim + observer(i)%nxim/2.0_wp) + 1
    iy = floor(atan2(-ky,kz)*rad2deg/observer(i)%dyim + observer(i)%nyim/2.0_wp) + 1

    !--- xfreq_ref = lab-frame frequency for spectral binning.
    !--- For direct light the source is always in vacuum (icell_clump=0),
    !--- so peel_xfreq_lab returns photon%xfreq unchanged.
    icl       = photon%icell_clump
    xfreq_ref = peel_xfreq_lab(photon%xfreq, pobs%kx, pobs%ky, pobs%kz, icl)
    xfreq_ref = xfreq_ref * (grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)
    ixf       = floor((xfreq_ref - grid%xfreq_min) / grid%dxfreq) + 1

    if (ix >= 1 .and. ix <= observer(i)%nxim .and. iy >= 1 .and. iy <= observer(i)%nyim) then
       call raytrace_to_edge(pobs, grid, tau)
       wgt0 = 1.0_wp / (fourpi*r2) * photon%wgt
       wgt  = exp(-tau) * wgt0

       if (par%save_peeloff_2D) then
          !$OMP ATOMIC UPDATE
          observer(i)%direc_2D(ix,iy) = observer(i)%direc_2D(ix,iy) + wgt
          if (par%use_stokes) then
             !$OMP ATOMIC UPDATE
             observer(i)%I_2D(ix,iy) = observer(i)%I_2D(ix,iy) + wgt
          endif
          if (par%save_direc0) then
             !$OMP ATOMIC UPDATE
             observer(i)%direc0_2D(ix,iy) = observer(i)%direc0_2D(ix,iy) + wgt0
          endif
       endif

       if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
          !$OMP ATOMIC UPDATE
          observer(i)%direc(ixf,ix,iy) = observer(i)%direc(ixf,ix,iy) + wgt
          if (par%use_stokes) then
             !$OMP ATOMIC UPDATE
             observer(i)%I(ixf,ix,iy) = observer(i)%I(ixf,ix,iy) + wgt
          endif
          if (par%save_direc0) then
             !$OMP ATOMIC UPDATE
             observer(i)%direc0(ixf,ix,iy) = observer(i)%direc0(ixf,ix,iy) + wgt0
          endif
       endif
    endif
  enddo
  end subroutine peeling_direct_clump

  !--------------------------------------------------
  subroutine peeling_dust_stokes_clump(photon, grid)
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
  real(wp) :: xfreq_ref
  real(wp) :: xx, yy, zz, rr, r_dot_k, det
  integer  :: ix, iy, ixf, icl
  integer  :: i

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

    if (trim(par%source_geometry) == 'stellar_illumination') then
       xx      = pobs%x
       yy      = pobs%y
       zz      = pobs%z + par%distance_star_to_planet
       rr      = sqrt(xx**2 + yy**2 + zz**2)
       r_dot_k = xx*pobs%kx + yy*pobs%ky + zz*pobs%kz
       det     = r_dot_k**2 - (rr**2 - par%stellar_radius**2)
       if (r_dot_k < 0.0_wp .and. det >= 0.0_wp) cycle
    endif

    kx = observer(i)%rmatrix(1,1)*pobs%kx + observer(i)%rmatrix(1,2)*pobs%ky + observer(i)%rmatrix(1,3)*pobs%kz
    ky = observer(i)%rmatrix(2,1)*pobs%kx + observer(i)%rmatrix(2,2)*pobs%ky + observer(i)%rmatrix(2,3)*pobs%kz
    kz = observer(i)%rmatrix(3,1)*pobs%kx + observer(i)%rmatrix(3,2)*pobs%ky + observer(i)%rmatrix(3,3)*pobs%kz

    ix = floor(atan2(-kx,kz)*rad2deg/observer(i)%dxim + observer(i)%nxim/2.0_wp) + 1
    iy = floor(atan2(-ky,kz)*rad2deg/observer(i)%dyim + observer(i)%nyim/2.0_wp) + 1

    !--- Lab-frame frequency: add clump bulk velocity along peel-off direction.
    icl       = photon%icell_clump
    xfreq_ref = peel_xfreq_lab(photon%xfreq, pobs%kx, pobs%ky, pobs%kz, icl)
    xfreq_ref = xfreq_ref * (grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)
    ixf       = floor((xfreq_ref - grid%xfreq_min) / grid%dxfreq) + 1

    if (ix >= 1 .and. ix <= observer(i)%nxim .and. iy >= 1 .and. iy <= observer(i)%nyim) then
       cost = photon%kx*pobs%kx + photon%ky*pobs%ky + photon%kz*pobs%kz
       sint = sqrt(1.0_wp - cost*cost)

       if (sint == 0.0_wp) then
          cosp  = 1.0_wp
          sinp  = 0.0_wp
          cos2p = 1.0_wp
          sin2p = 0.0_wp
       else
          cosp  = (pobs%kx*photon%mx + pobs%ky*photon%my + pobs%kz*photon%mz) / sint
          sinp  = (pobs%kx*photon%nx + pobs%ky*photon%ny + pobs%kz*photon%nz) / sint
          cos2p = 2.0_wp*cosp*cosp - 1.0_wp
          sin2p = 2.0_wp*cosp*sinp
       endif

       pobs%nx = -sinp*photon%mx + cosp*photon%nx
       pobs%ny = -sinp*photon%my + cosp*photon%ny
       pobs%nz = -sinp*photon%mz + cosp*photon%nz

       call interp_eq(scatt_mat%coss, scatt_mat%S11, cost, S11)
       call interp_eq(scatt_mat%coss, scatt_mat%S12, cost, S12)
       call interp_eq(scatt_mat%coss, scatt_mat%S33, cost, S33)
       call interp_eq(scatt_mat%coss, scatt_mat%S34, cost, S34)

       Q0 =  cos2p*photon%Q + sin2p*photon%U
       U0 = -sin2p*photon%Q + cos2p*photon%U

       Iobs = (S11*photon%I + S12*Q0      ) / twopi
       Qobs = (S12*photon%I + S11*Q0      ) / twopi
       Uobs = (S33*U0       + S34*photon%V) / twopi
       Vobs = (-S34*U0      + S33*photon%V) / twopi

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
       endif

       if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt(ixf,ix,iy) = observer(i)%scatt(ixf,ix,iy) + wgt * Idet
          if (par%save_dust_scattered) then
             !$OMP ATOMIC UPDATE
             observer(i)%scatt_dust(ixf,ix,iy) = observer(i)%scatt_dust(ixf,ix,iy) + wgt * Idet
          endif
          !$OMP ATOMIC UPDATE
          observer(i)%I(ixf,ix,iy) = observer(i)%I(ixf,ix,iy) + wgt * Idet
          !$OMP ATOMIC UPDATE
          observer(i)%Q(ixf,ix,iy) = observer(i)%Q(ixf,ix,iy) + wgt * Qdet
          !$OMP ATOMIC UPDATE
          observer(i)%U(ixf,ix,iy) = observer(i)%U(ixf,ix,iy) + wgt * Udet
          !$OMP ATOMIC UPDATE
          observer(i)%V(ixf,ix,iy) = observer(i)%V(ixf,ix,iy) + wgt * Vdet
       endif
    endif
  enddo
  end subroutine peeling_dust_stokes_clump

  !--------------------------------------------------
  !--- xfreq_atom = frequency in the scattering atom's rest frame (clump rest frame).
  !--- vel_atom   = velocity of the scattering atom in the frame comoving with the clump.
  subroutine peeling_resonance_stokes_clump(photon, grid, xfreq_atom, vel_atom)
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
  real(wp) :: xfreq, xfreq_ref, g_recoil, u_los
  real(wp) :: xx, yy, zz, rr, r_dot_k, det
  integer  :: ix, iy, ixf, icl
  integer  :: i

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

    if (trim(par%source_geometry) == 'stellar_illumination') then
       xx      = pobs%x
       yy      = pobs%y
       zz      = pobs%z + par%distance_star_to_planet
       rr      = sqrt(xx**2 + yy**2 + zz**2)
       r_dot_k = xx*pobs%kx + yy*pobs%ky + zz*pobs%kz
       det     = r_dot_k**2 - (rr**2 - par%stellar_radius**2)
       if (r_dot_k < 0.0_wp .and. det >= 0.0_wp) cycle
    endif

    kx = observer(i)%rmatrix(1,1)*pobs%kx + observer(i)%rmatrix(1,2)*pobs%ky + observer(i)%rmatrix(1,3)*pobs%kz
    ky = observer(i)%rmatrix(2,1)*pobs%kx + observer(i)%rmatrix(2,2)*pobs%ky + observer(i)%rmatrix(2,3)*pobs%kz
    kz = observer(i)%rmatrix(3,1)*pobs%kx + observer(i)%rmatrix(3,2)*pobs%ky + observer(i)%rmatrix(3,3)*pobs%kz

    ix = floor(atan2(-kx,kz)*rad2deg/observer(i)%dxim + observer(i)%nxim/2.0_wp) + 1
    iy = floor(atan2(-ky,kz)*rad2deg/observer(i)%dyim + observer(i)%nyim/2.0_wp) + 1

    if (ix >= 1 .and. ix <= observer(i)%nxim .and. iy >= 1 .and. iy <= observer(i)%nyim) then
       cost = photon%kx*pobs%kx + photon%ky*pobs%ky + photon%kz*pobs%kz
       sint = sqrt(1.0_wp - cost*cost)

       if (sint == 0.0_wp) then
          cosp  = 1.0_wp
          sinp  = 0.0_wp
          cos2p = 1.0_wp
          sin2p = 0.0_wp
       else
          cosp  = (pobs%kx*photon%mx + pobs%ky*photon%my + pobs%kz*photon%mz) / sint
          sinp  = (pobs%kx*photon%nx + pobs%ky*photon%ny + pobs%kz*photon%nz) / sint
          cos2p = 2.0_wp*cosp*cosp - 1.0_wp
          sin2p = 2.0_wp*cosp*sinp
       endif

       pobs%nx = -sinp*photon%mx + cosp*photon%nx
       pobs%ny = -sinp*photon%my + cosp*photon%ny
       pobs%nz = -sinp*photon%mz + cosp*photon%nz

       cost2 = cost**2
       S22   = three_over_four * photon%E1 * (cost2 + 1.0_wp)
       S11   = S22 + photon%E2
       S12   = three_over_four * photon%E1 * (cost2 - 1.0_wp)
       S33   = three_over_two  * photon%E1 * cost
       S44   = three_over_two  * photon%E3 * cost

       !--- Frequency in the clump rest frame toward the observer (peel-off direction).
       xfreq = xfreq_atom + (vel_atom(1)*cosp + vel_atom(2)*sinp)*sint + vel_atom(3)*cost
       if (par%recoil) then
          g_recoil = line%g_recoil0 / grid%Dfreq(photon%icell,photon%jcell,photon%kcell)
          xfreq    = xfreq - g_recoil * (1.0_wp - cost)
       endif

       !--- Convert to lab-frame frequency by adding the clump bulk velocity along the peel-off direction.
       icl       = photon%icell_clump
       xfreq_ref = peel_xfreq_lab(xfreq, pobs%kx, pobs%ky, pobs%kz, icl)
       xfreq_ref = xfreq_ref * (grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)
       ixf       = floor((xfreq_ref - grid%xfreq_min) / grid%dxfreq) + 1

       Q0 =  cos2p*photon%Q + sin2p*photon%U
       U0 = -sin2p*photon%Q + cos2p*photon%U
       Iobs = (S11 + S12*Q0) / fourpi
       Qobs = (S12 + S22*Q0) / fourpi
       Uobs = (S33*U0       ) / fourpi
       Vobs = (S44*photon%V ) / fourpi

       cosg  = -(observer(i)%rmatrix(1,1)*pobs%nx + observer(i)%rmatrix(1,2)*pobs%ny + observer(i)%rmatrix(1,3)*pobs%nz)
       sing  =   observer(i)%rmatrix(2,1)*pobs%nx + observer(i)%rmatrix(2,2)*pobs%ny + observer(i)%rmatrix(2,3)*pobs%nz
       cos2g = 2.0_wp*cosg*cosg - 1.0_wp
       sin2g = 2.0_wp*cosg*sing

       Idet =  Iobs
       Qdet =  cos2g*Qobs + sin2g*Uobs
       Udet = -sin2g*Qobs + cos2g*Uobs
       Vdet =  Vobs

       !--- Tau along peel-off: use clump-frame xfreq (raytrace_to_edge_clump handles clump exit shift).
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
       endif

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
       endif
    endif
  enddo
  end subroutine peeling_resonance_stokes_clump

  !--------------------------------------------------
  subroutine peeling_dust_nostokes_clump(photon, grid)
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  type(photon_type) :: pobs
  real(wp) :: r2, r, kx, ky, kz
  real(wp) :: cosa, wgt, peel, tau
  real(wp) :: xfreq_ref
  real(wp) :: xx, yy, zz, rr, r_dot_k, det
  integer  :: ix, iy, ixf, icl
  integer  :: i

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

    if (trim(par%source_geometry) == 'stellar_illumination') then
       xx      = pobs%x
       yy      = pobs%y
       zz      = pobs%z + par%distance_star_to_planet
       rr      = sqrt(xx**2 + yy**2 + zz**2)
       r_dot_k = xx*pobs%kx + yy*pobs%ky + zz*pobs%kz
       det     = r_dot_k**2 - (rr**2 - par%stellar_radius**2)
       if (r_dot_k < 0.0_wp .and. det >= 0.0_wp) cycle
    endif

    kx = observer(i)%rmatrix(1,1)*pobs%kx + observer(i)%rmatrix(1,2)*pobs%ky + observer(i)%rmatrix(1,3)*pobs%kz
    ky = observer(i)%rmatrix(2,1)*pobs%kx + observer(i)%rmatrix(2,2)*pobs%ky + observer(i)%rmatrix(2,3)*pobs%kz
    kz = observer(i)%rmatrix(3,1)*pobs%kx + observer(i)%rmatrix(3,2)*pobs%ky + observer(i)%rmatrix(3,3)*pobs%kz

    ix = floor(atan2(-kx,kz)*rad2deg/observer(i)%dxim + observer(i)%nxim/2.0_wp) + 1
    iy = floor(atan2(-ky,kz)*rad2deg/observer(i)%dyim + observer(i)%nyim/2.0_wp) + 1

    if (ix >= 1 .and. ix <= observer(i)%nxim .and. iy >= 1 .and. iy <= observer(i)%nyim) then
       !--- Lab-frame frequency: add clump bulk velocity along peel-off direction.
       icl       = photon%icell_clump
       xfreq_ref = peel_xfreq_lab(photon%xfreq, pobs%kx, pobs%ky, pobs%kz, icl)
       xfreq_ref = xfreq_ref * (grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)
       ixf       = floor((xfreq_ref - grid%xfreq_min) / grid%dxfreq) + 1

       call raytrace_to_edge(pobs, grid, tau)
       cosa = photon%kx*pobs%kx + photon%ky*pobs%ky + photon%kz*pobs%kz
       peel = (1.0_wp - par%hgg**2) / ((1.0_wp + par%hgg**2) - 2.0_wp*par%hgg*cosa)**1.5_wp / fourpi
       wgt  = peel / r2 * exp(-tau) * photon%wgt

       if (par%save_peeloff_2D) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt_2D(ix,iy) = observer(i)%scatt_2D(ix,iy) + wgt
       endif

       if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt(ixf,ix,iy) = observer(i)%scatt(ixf,ix,iy) + wgt
          if (par%save_dust_scattered) then
             !$OMP ATOMIC UPDATE
             observer(i)%scatt_dust(ixf,ix,iy) = observer(i)%scatt_dust(ixf,ix,iy) + wgt
          endif
       endif
    endif
  enddo
  end subroutine peeling_dust_nostokes_clump

  !--------------------------------------------------
  subroutine peeling_resonance_nostokes_clump(photon, grid, xfreq_atom, vel_atom)
  use random
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  real(wp),          intent(in) :: xfreq_atom
  real(wp),          intent(in) :: vel_atom(3)
  type(photon_type) :: pobs
  real(wp) :: r2, r
  real(wp) :: wgt, peel, tau
  real(wp) :: cost, cost2, sint, cosp, sinp, rho1, rho
  real(wp) :: kx, ky, kz
  real(wp) :: xfreq, xfreq_ref, g_recoil
  real(wp) :: xx, yy, zz, rr, r_dot_k, det
  integer  :: ix, iy, ixf, icl
  integer  :: i

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

    if (trim(par%source_geometry) == 'stellar_illumination') then
       xx      = pobs%x
       yy      = pobs%y
       zz      = pobs%z + par%distance_star_to_planet
       rr      = sqrt(xx**2 + yy**2 + zz**2)
       r_dot_k = xx*pobs%kx + yy*pobs%ky + zz*pobs%kz
       det     = r_dot_k**2 - (rr**2 - par%stellar_radius**2)
       if (r_dot_k < 0.0_wp .and. det >= 0.0_wp) cycle
    endif

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
          cosp  = rho * ( cost*photon%kz - pobs%kz)
          sinp  = rho * (photon%kx*pobs%ky - pobs%kx*photon%ky)
          xfreq = xfreq_atom + (vel_atom(1)*cosp + vel_atom(2)*sinp)*sint + vel_atom(3)*cost
       endif

       if (par%recoil) then
          g_recoil = line%g_recoil0 / grid%Dfreq(photon%icell,photon%jcell,photon%kcell)
          xfreq    = xfreq - g_recoil * (1.0_wp - cost)
       endif

       !--- Convert clump-frame xfreq to lab frame for spectral binning.
       icl       = photon%icell_clump
       xfreq_ref = peel_xfreq_lab(xfreq, pobs%kx, pobs%ky, pobs%kz, icl)
       xfreq_ref = xfreq_ref * (grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)
       ixf       = floor((xfreq_ref - grid%xfreq_min) / grid%dxfreq) + 1

       !--- Tau along peel-off: use clump-frame xfreq (raytrace_to_edge_clump handles clump exit shift).
       pobs%xfreq = xfreq
       call raytrace_to_edge(pobs, grid, tau)
       peel = three_over_four * photon%E1 * (cost2 + 1.0_wp) + photon%E2
       wgt  = peel / (fourpi*r2) * exp(-tau) * photon%wgt

       if (par%save_peeloff_2D) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt_2D(ix,iy) = observer(i)%scatt_2D(ix,iy) + wgt
       endif

       if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt(ixf,ix,iy) = observer(i)%scatt(ixf,ix,iy) + wgt
       endif
    endif
  enddo
  end subroutine peeling_resonance_nostokes_clump

end module peelingoff_clump_mod
