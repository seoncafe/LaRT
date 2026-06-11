!-- Modification History
!   2018-02-03, Frequency binning for peeling-off image should be done in observer's frame (lab frame).
!--
module peelingoff_heal
  use define
  use mathlib
  use utility
  use line_mod
  use healpix
contains
  !--------------------------------------------------
  subroutine peeling_direct_inside(photon,grid)
  implicit none
  type(photon_type),   intent(in)    :: photon
  type(grid_type),     intent(in)    :: grid
  ! local variables
  type(photon_type) :: pobs
  real(kind=wp) :: r2,r,wgt,wgt0,tau
  real(kind=wp) :: kx,ky,kz
  real(kind=wp) :: xfreq_ref, u1, u2
  integer :: ipix,ixf
  integer :: i

  do i=1,par%nobs
    pobs    = photon
    pobs%kx = (observer(i)%x-photon%x)
    pobs%ky = (observer(i)%y-photon%y)
    pobs%kz = (observer(i)%z-photon%z)
    r2      = pobs%kx*pobs%kx + pobs%ky*pobs%ky + pobs%kz*pobs%kz
    r       = sqrt(r2)
    pobs%kx = pobs%kx/r
    pobs%ky = pobs%ky/r
    pobs%kz = pobs%kz/r

    !--- location in healpix
    call vec2pix(observer(i)%nside, -pobs%kx, -pobs%ky, -pobs%kz, ipix)

    !--- xfreq_ref = lab (observer) frame frequency
     if (.not. par%comoving_source) then
        !--- transform the frequency back to the lab frame.
        u1 = grid%vfx(photon%icell,photon%jcell,photon%kcell)*photon%kx + &
             grid%vfy(photon%icell,photon%jcell,photon%kcell)*photon%ky + &
             grid%vfz(photon%icell,photon%jcell,photon%kcell)*photon%kz
        xfreq_ref = photon%xfreq + u1
        !--- frequency of the photon propagating toward the observer in the comoving frame
        u2 = grid%vfx(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kx + &
             grid%vfy(pobs%icell,pobs%jcell,pobs%kcell)*pobs%ky + &
             grid%vfz(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kz
        pobs%xfreq = xfreq_ref - u2
     else
        !--- for comoving source, the propagation vector of a peeled-off, direct ray is always toward the observer.
        u1 = grid%vfx(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kx + &
             grid%vfy(pobs%icell,pobs%jcell,pobs%kcell)*pobs%ky + &
             grid%vfz(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kz
        xfreq_ref = photon%xfreq + u1
     endif

     !--- frequency bin
     xfreq_ref = xfreq_ref * (grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)
     ixf       = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1

     if (ipix >= 1 .and. ipix <= observer(i)%npix) then
        call raytrace_to_dist(pobs,grid,r,tau)
        !-- bug-fixed (2022.04.26), direc was the same as direc0 when par%save_direc0 = .true.
        wgt0 = 1.0_wp/(fourpi*r2) * photon%wgt
        wgt  = exp(-tau)*wgt0

        !--- 2D image
        if (par%save_peeloff_2D) then
           !$OMP ATOMIC UPDATE
           observer(i)%direc_heal_2D(ipix) = observer(i)%direc_heal_2D(ipix) + wgt
           !if (par%use_stokes) then
           !   !$OMP ATOMIC UPDATE
           !   observer(i)%I_2D(ipix) = observer(i)%I_2D(ipix) + wgt
           !endif
           if (par%save_direc0) then
              !$OMP ATOMIC UPDATE
              observer(i)%direc0_heal_2D(ipix) = observer(i)%direc0_heal_2D(ipix) + wgt0
           endif
        endif

        !--- 3D spectral image
        if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
           !$OMP ATOMIC UPDATE
           observer(i)%direc_heal(ixf,ipix) = observer(i)%direc_heal(ixf,ipix) + wgt
           !if (par%use_stokes) then
           !   !$OMP ATOMIC UPDATE
           !   observer(i)%I_heal(ixf,ipix) = observer(i)%I_heal(ixf,ipix) + wgt
           !endif
           if (par%save_direc0) then
              !$OMP ATOMIC UPDATE
              observer(i)%direc0_heal(ixf,ipix) = observer(i)%direc0_heal(ixf,ipix) + wgt0
           endif
        endif
     endif
  enddo
  end subroutine peeling_direct_inside
  !--------------------------------------------------
  !--- We need to think about how to deal with Stokes parameters in Healpix coordinates.
  !--------------------------------------------------
  subroutine peeling_dust_nostokes_inside(photon,grid)
  implicit none
  type(photon_type),   intent(in)    :: photon
  type(grid_type),     intent(in)    :: grid
  ! local variables
  type (photon_type) :: pobs
  real(kind=wp) :: r2,r,kx,ky,kz
  real(kind=wp) :: cosa,wgt,peel,tau
  real(kind=wp) :: xfreq_ref, u1
  real(kind=wp) :: xx,yy,zz,rr,r_dot_k,det
  integer       :: ipix,ixf
  integer       :: i

  do i=1,par%nobs
    pobs    = photon
    pobs%kx = (observer(i)%x-photon%x)
    pobs%ky = (observer(i)%y-photon%y)
    pobs%kz = (observer(i)%z-photon%z)
    r2      = pobs%kx*pobs%kx + pobs%ky*pobs%ky + pobs%kz*pobs%kz
    r       = sqrt(r2)
    pobs%kx = pobs%kx/r
    pobs%ky = pobs%ky/r
    pobs%kz = pobs%kz/r

    !--- location in healpix
    call vec2pix(observer(i)%nside, -pobs%kx, -pobs%ky, -pobs%kz, ipix)

    if (ipix >= 1 .and. ipix <= observer(i)%npix) then
       !--- frequency for spectral binning.
       !--- xfreq_ref = lab (observer) frame frequency of peel-off photon (2020.08.28, bug-fixed)
       u1 = grid%vfx(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kx + &
            grid%vfy(pobs%icell,pobs%jcell,pobs%kcell)*pobs%ky + &
            grid%vfz(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kz
       xfreq_ref = photon%xfreq + u1

       !--- Note that frequency for the calculation of optical depth along the peel-off direction is
       !--- pobs%xfreq = photon%xfreq because dust grains are assumed to have no thermal motion.
       !--- comment added on 2020.08.23.

       !--- frequency bin
       xfreq_ref = xfreq_ref * (grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)
       ixf       = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1

       call raytrace_to_dist(pobs,grid,r,tau)
       cosa = photon%kx*pobs%kx+photon%ky*pobs%ky+photon%kz*pobs%kz
       peel = (1.0_wp - par%hgg**2)/((1.0_wp + par%hgg**2)-2.0_wp*par%hgg*cosa)**1.5_wp/fourpi
       !--- albedo was already multiplied before this routine is called.
       wgt  = peel/r2 * exp(-tau) * photon%wgt

       !--- 2D image
       if (par%save_peeloff_2D) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt_heal_2D(ipix) = observer(i)%scatt_heal_2D(ipix) + wgt
       endif

       !--- 3D spectral image
       if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt_heal(ixf,ipix) = observer(i)%scatt_heal(ixf,ipix) + wgt
       endif
    endif
  enddo
  end subroutine peeling_dust_nostokes_inside
  !--------------------------------------------------
  subroutine peeling_resonance_nostokes_inside(photon,grid,xfreq_atom,vel_atom)
  use random
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  real(kind=wp),     intent(in) :: xfreq_atom
  real(kind=wp),     intent(in) :: vel_atom(3)
  ! local variables
  type (photon_type) :: pobs
  real(kind=wp) :: r2,r
  real(kind=wp) :: wgt,peel,tau
  real(kind=wp) :: cost,cost2,sint,cosp,sinp,rho1,rho
  real(kind=wp) :: kx,ky,kz
  real(kind=wp) :: xfreq_ref, xfreq, g_recoil, u1
  real(kind=wp) :: xx,yy,zz,rr,r_dot_k,det
  integer       :: ipix,ixf
  integer       :: i

  do i=1,par%nobs
    pobs    = photon
    pobs%kx = (observer(i)%x-photon%x)
    pobs%ky = (observer(i)%y-photon%y)
    pobs%kz = (observer(i)%z-photon%z)
    r2      = pobs%kx*pobs%kx + pobs%ky*pobs%ky + pobs%kz*pobs%kz
    r       = sqrt(r2)
    pobs%kx = pobs%kx/r
    pobs%ky = pobs%ky/r
    pobs%kz = pobs%kz/r

    !--- location in healpix
    call vec2pix(observer(i)%nside, -pobs%kx, -pobs%ky, -pobs%kz, ipix)

    if (ipix >= 1 .and. ipix <= observer(i)%npix) then
       cost  = photon%kx * pobs%kx + photon%ky * pobs%ky + photon%kz * pobs%kz
       cost2 = cost**2
       sint  = sqrt(1.0_wp - cost2)
       rho1  = sqrt(1.0_wp - photon%kz**2) * sint

       !--- Calculate azimuthal scattering angle toward the observer.
       !--- bug-fixed (2021.05.03), the case where the photon direction and observer direction are coincident was missing.
       if (rho1 == 0.0_wp) then
          cosp  = 1.0_wp
          sinp  = 0.0_wp
          xfreq = xfreq_atom + (vel_atom(1)*cosp + vel_atom(2)*sinp)*sint + vel_atom(3)*cost
       else
          rho   = 1.0_wp/rho1
          cosp  = rho * ( cost*photon%kz -pobs%kz)
          sinp  = rho * (photon%kx*pobs%ky - pobs%kx*photon%ky)
          xfreq = xfreq_atom + (vel_atom(1)*cosp + vel_atom(2)*sinp)*sint + vel_atom(3)*cost
       endif

       if (par%recoil) then
          g_recoil = line%g_recoil0 /grid%Dfreq(photon%icell,photon%jcell,photon%kcell)
          xfreq    = xfreq - g_recoil * (1.0_wp - cost)
       endif

       !--- xfreq_ref = lab (observer) frame frequency of peel-off photon.
       u1 = grid%vfx(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kx + &
            grid%vfy(pobs%icell,pobs%jcell,pobs%kcell)*pobs%ky + &
            grid%vfz(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kz
       xfreq_ref = xfreq + u1

       !--- frequency bin
       xfreq_ref = xfreq_ref * (grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)
       ixf       = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1

       !----------------------------------------------------------------------------------
       pobs%xfreq = xfreq
       call raytrace_to_dist(pobs,grid,r,tau)
       peel = three_over_four * photon%E1 * (cost2 + 1.0_wp) + photon%E2
       wgt  = peel/(fourpi*r2) * exp(-tau) * photon%wgt

       !--- 2D image
       if (par%save_peeloff_2D) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt_heal_2D(ipix) = observer(i)%scatt_heal_2D(ipix) + wgt
       endif

       !--- 3D spectral image
       if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt_heal(ixf,ipix) = observer(i)%scatt_heal(ixf,ipix) + wgt
       endif
    endif
  enddo
  end subroutine peeling_resonance_nostokes_inside
  !--------------------------------------------------
end module peelingoff_heal
