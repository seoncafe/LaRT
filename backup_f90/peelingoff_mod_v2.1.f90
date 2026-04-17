!-- Modification History
!   2018-02-03, Frequency binning for peeling-off image should be done in observer's frame (lab frame).
!--
module peelingoff_mod
  use define
  use mathlib
  use utility
  use memory_mod
contains
#ifdef PEELINGOFF
  !--------------------------------------------------
  subroutine peeling_direct(photon,grid)
  implicit none
  type(photon_type),   intent(in)    :: photon
  type(grid_type),     intent(in)    :: grid
  ! local variables
  type(photon_type) :: pobs
  real(kind=wp) :: r2,r,wgt,tau
  real(kind=wp) :: kx,ky,kz,vdet(3)
  real(kind=wp) :: xfreq_ref, u1, u2
  integer :: ix,iy,ixf
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

    !--- Transform the peeling-off vector to the observer's frame.
    kx = observer(i)%rmatrix(1,1)*pobs%kx + observer(i)%rmatrix(1,2)*pobs%ky + observer(i)%rmatrix(1,3)*pobs%kz
    ky = observer(i)%rmatrix(2,1)*pobs%kx + observer(i)%rmatrix(2,2)*pobs%ky + observer(i)%rmatrix(2,3)*pobs%kz
    kz = observer(i)%rmatrix(3,1)*pobs%kx + observer(i)%rmatrix(3,2)*pobs%ky + observer(i)%rmatrix(3,3)*pobs%kz

    !--- Location in TAN image.
    ix = floor(atan2(-kx,kz)*rad2deg/observer(i)%dxim+observer(i)%nxim/2.0_wp) + 1
    iy = floor(atan2(-ky,kz)*rad2deg/observer(i)%dyim+observer(i)%nyim/2.0_wp) + 1

    !--- bug-fixed, 2020.08.19
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
     xfreq_ref = xfreq_ref * grid%Dfreq_ratio(photon%icell,photon%jcell,photon%kcell)
     ixf       = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1

     if (ix >= 1 .and. ix <= observer(i)%nxim .and. iy >= 1 .and. iy <= observer(i)%nyim &
         .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
        call raytrace_to_edge(pobs,grid,tau)
        wgt = exp(-tau)/(fourpi*r2) * photon%wgt
        observer(i)%direc(ixf,ix,iy) = observer(i)%direc(ixf,ix,iy) + wgt
        if (par%use_stokes) then
           observer(i)%I(ixf,ix,iy) = observer(i)%I(ixf,ix,iy) + wgt
        endif
        if (par%save_direc0) then
           wgt = 1.0_wp/(fourpi*r2) * photon%wgt
           observer(i)%direc0(ixf,ix,iy) = observer(i)%direc0(ixf,ix,iy) + wgt
        endif
     endif
  enddo
  end subroutine peeling_direct
  !--------------------------------------------------
  subroutine peeling_dust_stokes(photon,grid)
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  ! local variables
  type(photon_type) :: pobs
  real(kind=wp) :: r2,r
  real(kind=wp) :: cost,sint,cosp,sinp,cos2p,sin2p,cosg,sing,cos2g,sin2g
  real(kind=wp) :: S11,S12,S33,S34
  real(kind=wp) :: Q0,U0
  real(kind=wp) :: Iobs,Qobs,Uobs,Vobs,Idet,Qdet,Udet,Vdet
  real(kind=wp) :: kx,ky,kz,nx,ny
  real(kind=wp) :: tau,wgt
  real(kind=wp) :: xfreq_ref, u1
  integer :: ix,iy,ixf
  integer :: i

  do i=1,par%nobs
    !--- Calculate the propagation vector for the peeling-off.
    pobs    = photon
    pobs%kx = (observer(i)%x-photon%x)
    pobs%ky = (observer(i)%y-photon%y)
    pobs%kz = (observer(i)%z-photon%z)
    r2      = pobs%kx**2 + pobs%ky**2 + pobs%kz**2
    r       = sqrt(r2)
    pobs%kx = pobs%kx/r
    pobs%ky = pobs%ky/r
    pobs%kz = pobs%kz/r

    !--- Transform the peeling-off vector to the observer's frame.
    kx = observer(i)%rmatrix(1,1)*pobs%kx + observer(i)%rmatrix(1,2)*pobs%ky + observer(i)%rmatrix(1,3)*pobs%kz
    ky = observer(i)%rmatrix(2,1)*pobs%kx + observer(i)%rmatrix(2,2)*pobs%ky + observer(i)%rmatrix(2,3)*pobs%kz
    kz = observer(i)%rmatrix(3,1)*pobs%kx + observer(i)%rmatrix(3,2)*pobs%ky + observer(i)%rmatrix(3,3)*pobs%kz

    !--- Location in TAN image.
    ix = floor(atan2(-kx,kz)*rad2deg/observer(i)%dxim+observer(i)%nxim/2.0_wp) + 1
    iy = floor(atan2(-ky,kz)*rad2deg/observer(i)%dyim+observer(i)%nyim/2.0_wp) + 1

    !--- frequency for spectral binning in the detector plane.
    !--- xfreq_ref = lab (observer) frame frequency of peel-off photon (2020.08.28, bug-fixed)
    u1 = grid%vfx(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kx + &
         grid%vfy(pobs%icell,pobs%jcell,pobs%kcell)*pobs%ky + &
         grid%vfz(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kz
    xfreq_ref = photon%xfreq + u1

    !--- Note that frequency for the calculation of optical depth along the peel-off direction is
    !--- pobs%xfreq = photon%xfreq because dust grains are assumed to have no thermal motion.
    !--- comment added on 2020.08.23.

    !--- frequency bin
    xfreq_ref = xfreq_ref * grid%Dfreq_ratio(photon%icell,photon%jcell,photon%kcell)
    ixf       = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq) + 1

    if (ix >= 1 .and. ix <= observer(i)%nxim .and. iy >= 1 .and. iy <= observer(i)%nyim &
       .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
       !--- Calculate scattering angle toward the observer.
       cost = photon%kx * pobs%kx + photon%ky * pobs%ky + photon%kz * pobs%kz
       sint = sqrt(1.0_wp - cost*cost)

       !--- Calculate the reference normal vector for peeling-off.
       if (sint /= 0.0_wp) then
          pobs%nx = (photon%ky * pobs%kz - photon%kz * pobs%ky)/sint
          pobs%ny = (photon%kz * pobs%kx - photon%kx * pobs%kz)/sint
          pobs%nz = (photon%kx * pobs%ky - photon%ky * pobs%kx)/sint
       else
          pobs%nx = photon%nx
          pobs%ny = photon%ny
          pobs%nz = photon%nz
       endif

       !--- Calculate azimuthal scattering angle toward the observer.
       cosp  =   pobs%nx * photon%nx + pobs%ny * photon%ny + pobs%nz * photon%nz
       sinp  = -(pobs%nx * photon%mx + pobs%ny * photon%my + pobs%nz * photon%mz)
       cos2p = 2.0_wp*cosp*cosp - 1.0_wp
       sin2p = 2.0_wp*cosp*sinp

       !--- Calculate scattering matrix for the peeling-off vector for the scattering angles.
       call interp_eq(scatt_mat%coss,scatt_mat%S11,cost,S11)
       call interp_eq(scatt_mat%coss,scatt_mat%S12,cost,S12)
       call interp_eq(scatt_mat%coss,scatt_mat%S33,cost,S33)
       call interp_eq(scatt_mat%coss,scatt_mat%S34,cost,S34)

       !--- Calculate Stokes parameters for the peeling-off vector.
       !--- Stokes S11 is normalized by Integral(S11(cos\theta)) d(cos\theta).
       !--- Therefore, 1/2pi should be multiplied for 4pi factor.
       Q0 =  cos2p*photon%Q + sin2p*photon%U
       U0 = -sin2p*photon%Q + cos2p*photon%U

       Iobs = ( S11*photon%I + S12*Q0      )/twopi
       Qobs = ( S12*photon%I + S11*Q0      )/twopi
       Uobs = ( S33*U0       + S34*photon%V)/twopi
       Vobs = (-S34*U0       + S33*photon%V)/twopi

       !--- Rotate the Stokes parameters to the detector plane.
       !--- The Stokes parameters are defined by counter-clock rotation from the North to the South.
       !--- IAU 1974 recommendation.
       !nx = observer(i)%rmatrix(1,1)*pobs%nx + observer(i)%rmatrix(1,2)*pobs%ny + observer(i)%rmatrix(1,3)*pobs%nz
       !ny = observer(i)%rmatrix(2,1)*pobs%nx + observer(i)%rmatrix(2,2)*pobs%ny + observer(i)%rmatrix(2,3)*pobs%nz
       !nz = observer(i)%rmatrix(3,1)*pobs%nx + observer(i)%rmatrix(3,2)*pobs%ny + observer(i)%rmatrix(3,3)*pobs%nz
       !cosg  = -nx
       !sing  =  ny
       cosg  = -(observer(i)%rmatrix(1,1)*pobs%nx + observer(i)%rmatrix(1,2)*pobs%ny + observer(i)%rmatrix(1,3)*pobs%nz)
       sing  =   observer(i)%rmatrix(2,1)*pobs%nx + observer(i)%rmatrix(2,2)*pobs%ny + observer(i)%rmatrix(2,3)*pobs%nz
       cos2g = 2.0_wp*cosg*cosg - 1.0_wp
       sin2g = 2.0_wp*cosg*sing

       Idet = Iobs
       Qdet =  cos2g*Qobs + sin2g*Uobs
       Udet = -sin2g*Qobs + cos2g*Uobs
       Vdet = Vobs

       !--- Place a fraction to be peeled off to the output Stokes images.
       !--- note the angle-dependent, peeled fraction is took into account in Idet, Qdet, Udet, and Vdet.
       call raytrace_to_edge(pobs,grid,tau)
       wgt  = 1.0_wp/r2 * exp(-tau) * photon%wgt

       observer(i)%scatt(ixf,ix,iy) = observer(i)%scatt(ixf,ix,iy) + wgt * Idet
       observer(i)%I(ixf,ix,iy)     = observer(i)%I(ixf,ix,iy)     + wgt * Idet
       observer(i)%Q(ixf,ix,iy)     = observer(i)%Q(ixf,ix,iy)     + wgt * Qdet
       observer(i)%U(ixf,ix,iy)     = observer(i)%U(ixf,ix,iy)     + wgt * Udet
       observer(i)%V(ixf,ix,iy)     = observer(i)%V(ixf,ix,iy)     + wgt * Vdet
    endif
  enddo
  end subroutine peeling_dust_stokes
  !--------------------------------------------------
  !--- xfreq_atom = frequency in the scattering atom's rest frame.
  !--- vel_atom   = velocity of the scattering atom in the referenc frame carried with photon.
  subroutine peeling_resonance_stokes(photon,grid,xfreq_atom,vel_atom)
  use random
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  real(kind=wp),     intent(in) :: xfreq_atom
  real(kind=wp),     intent(in) :: vel_atom(3)
  ! local variables
  type(photon_type) :: pobs
  real(kind=wp) :: r2,r
  real(kind=wp) :: cost,cost2,sint,cosp,sinp,cos2p,sin2p,cosg,sing,cos2g,sin2g
  real(kind=wp) :: S11,S12,S22,S33,S44
  real(kind=wp) :: Q0,U0
  real(kind=wp) :: Iobs,Qobs,Uobs,Vobs,Idet,Qdet,Udet,Vdet
  real(kind=wp) :: kx,ky,kz,nx,ny
  real(kind=wp) :: tau,wgt
  real(kind=wp) :: xfreq, xfreq_ref, g_recoil, E1, u1
  integer :: ix,iy,ixf
  integer :: i
#ifdef FINE_STRUCTURE
  real(kind=wp) :: DnuHK,qK,qH,qH2
#endif
  real(kind=wp), parameter :: three_over_four = 3.0_wp/4.0_wp, &
                              three_over_two  = 3.0_wp/2.0_wp, &
                              one_over_three  = 1.0_wp/3.0_wp

  do i=1, par%nobs
    !--- Calculate the propagation vector for the peeling-off.
    pobs    = photon
    pobs%kx = (observer(i)%x-photon%x)
    pobs%ky = (observer(i)%y-photon%y)
    pobs%kz = (observer(i)%z-photon%z)
    r2      = pobs%kx**2 + pobs%ky**2 + pobs%kz**2
    r       = sqrt(r2)
    pobs%kx = pobs%kx/r
    pobs%ky = pobs%ky/r
    pobs%kz = pobs%kz/r

    !--- Transform the peeling-off vector to the observer's frame.
    kx = observer(i)%rmatrix(1,1)*pobs%kx + observer(i)%rmatrix(1,2)*pobs%ky + observer(i)%rmatrix(1,3)*pobs%kz
    ky = observer(i)%rmatrix(2,1)*pobs%kx + observer(i)%rmatrix(2,2)*pobs%ky + observer(i)%rmatrix(2,3)*pobs%kz
    kz = observer(i)%rmatrix(3,1)*pobs%kx + observer(i)%rmatrix(3,2)*pobs%ky + observer(i)%rmatrix(3,3)*pobs%kz

    !--- Location in TAN image.
    ix = floor(atan2(-kx,kz)*rad2deg/observer(i)%dxim+observer(i)%nxim/2.0_wp) + 1
    iy = floor(atan2(-ky,kz)*rad2deg/observer(i)%dyim+observer(i)%nyim/2.0_wp) + 1

    if (ix >= 1 .and. ix <= observer(i)%nxim .and. iy >= 1 .and. iy <= observer(i)%nyim) then
       !--- Calculate scattering angle toward the observer.
       cost = photon%kx * pobs%kx + photon%ky * pobs%ky + photon%kz * pobs%kz
       sint = sqrt(1.0_wp - cost*cost)

       !--- Calculate the reference normal vector for peeling-off.
       if (sint /= 0.0_wp) then
          pobs%nx = (photon%ky * pobs%kz - photon%kz * pobs%ky)/sint
          pobs%ny = (photon%kz * pobs%kx - photon%kx * pobs%kz)/sint
          pobs%nz = (photon%kx * pobs%ky - photon%ky * pobs%kx)/sint
       else
          pobs%nx = photon%nx
          pobs%ny = photon%ny
          pobs%nz = photon%nz
       endif

       !--- Calculate azimuthal scattering angle toward the observer.
       cosp  =   pobs%nx * photon%nx + pobs%ny * photon%ny + pobs%nz * photon%nz
       sinp  = -(pobs%nx * photon%mx + pobs%ny * photon%my + pobs%nz * photon%mz)
       cos2p = 2.0_wp*cosp*cosp - 1.0_wp
       sin2p = 2.0_wp*cosp*sinp

       !--- frequency in the fluid frame (comoving frame).
       xfreq = xfreq_atom + (vel_atom(1)*cosp + vel_atom(2)*sinp)*sint + vel_atom(3)*cost
       if (par%recoil) then
          g_recoil = g_recoil0 /(grid%Dfreq(photon%icell,photon%jcell,photon%kcell)*par%atom_no)
          xfreq    = xfreq - g_recoil * (1.0_wp - cost)
       endif

       !--- xfreq_ref = lab (observer) frame frequency of peel-off photon.
       u1 = grid%vfx(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kx + &
            grid%vfy(pobs%icell,pobs%jcell,pobs%kcell)*pobs%ky + &
            grid%vfz(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kz
       xfreq_ref = xfreq + u1

       !--- frequency bin
       xfreq_ref = xfreq_ref * grid%Dfreq_ratio(photon%icell,photon%jcell,photon%kcell)
       ixf       = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq) + 1
       if (ixf < 1 .or. ixf > grid%nxfreq) return
#ifdef FINE_STRUCTURE
       !--- Select a new scattering angle (cos(theta)).
       !--- E1 is a function of frequency at the atom's rest frame.
       !--- We are assuming nu_H = nu_K. The difference make no practical difference.
       !--- qH = xfreq_atom, qK = qH - DnuHK
       !=DnuHK = grid%DnuHK_ref/grid%Dfreq_ratio(photon%icell,photon%jcell,photon%kcell)
       !=qH2   = xfreq_atom**2
       !=qK    = xfreq_atom - DnuHK
       !=E1    = (2.0_wp*qK*xfreq_atom + qH2)/(qK**2 + 2.0_wp*qH2)
       DnuHK = grid%DnuHK_ref_half/grid%Dfreq_ratio(photon%icell,photon%jcell,photon%kcell)
       qH    = xfreq_atom + DnuHK
       qK    = xfreq_atom - DnuHK
       E1    = (2.0_wp*qK*qH + qH**2)/(qK**2 + 2.0_wp*qH**2)

       !--- Calculate Stokes parameters. (S11, S12, S22, S33, and S44)
       !--- Note S44 /= S33 and S34 = 0.
       !--- E2 = 1-E1 and E3 = (E1+2)/3. (3/2) E3 = (E1+2)/2.
       cost2 = cost**2
       S22   = three_over_four * E1 * (cost2 + 1.0_wp)
       S11   = S22 + (1.0_wp - E1)
       S12   = three_over_four * E1 * (cost2 - 1.0_wp)
       S33   = three_over_two  * E1 * cost
       S44   = (E1 + 2.0_wp)/2.0_wp * cost
#else
       E1    = 1.0_wp
       !--- Calculate Stokes parameters. (S11, S12, S33, and S44)
       !--- Note S44 /= S33 and S34 = 0.
       cost2 = cost**2
       S11   = three_over_four * (cost2 + 1.0_wp)
       S22   = S11
       S12   = three_over_four * (cost2 - 1.0_wp)
       S33   = three_over_two  * cost
       S44   = S33
#endif
       !--- Calculate Stokes parameters for the peeling-off vector.
       !--- Stokes S11 is normalized by Integral[S11(mu), {mu,-1,1}] = 2.
       !--- Therefore, 1/4pi should be multiplied.
       Q0 =  cos2p*photon%Q + sin2p*photon%U
       U0 = -sin2p*photon%Q + cos2p*photon%U
       !Iobs = (S11*photon%I + S12*Q0)/fourpi
       !Qobs = (S12*photon%I + S22*Q0)/fourpi
       Iobs = (S11 + S12*Q0)/fourpi
       Qobs = (S12 + S22*Q0)/fourpi
       Uobs = (S33*U0      )/fourpi
       Vobs = (S44*photon%V)/fourpi

       !--- Rotate the Stokes parameters to the detector plane.
       !--- The Stokes parameters are defined by counter-clock rotation from the North to the East.
       !--- (from Y to -X in detector plane)
       !--- IAU 1974 recommendation.
       !nx = observer(i)%rmatrix(1,1)*pobs%nx + observer(i)%rmatrix(1,2)*pobs%ny + observer(i)%rmatrix(1,3)*pobs%nz
       !ny = observer(i)%rmatrix(2,1)*pobs%nx + observer(i)%rmatrix(2,2)*pobs%ny + observer(i)%rmatrix(2,3)*pobs%nz
       !nz = observer(i)%rmatrix(3,1)*pobs%nx + observer(i)%rmatrix(3,2)*pobs%ny + observer(i)%rmatrix(3,3)*pobs%nz
       !cosg  = -nx
       !sing  =  ny
       cosg  = -(observer(i)%rmatrix(1,1)*pobs%nx + observer(i)%rmatrix(1,2)*pobs%ny + observer(i)%rmatrix(1,3)*pobs%nz)
       sing  =   observer(i)%rmatrix(2,1)*pobs%nx + observer(i)%rmatrix(2,2)*pobs%ny + observer(i)%rmatrix(2,3)*pobs%nz
       cos2g = 2.0_wp*cosg*cosg - 1.0_wp
       sin2g = 2.0_wp*cosg*sing

       Idet = Iobs
       Qdet =  cos2g*Qobs + sin2g*Uobs
       Udet = -sin2g*Qobs + cos2g*Uobs
       Vdet = Vobs

       !--- Place a fraction to be peeled off to the output Stokes images.
       !--- note the angle-dependent, peeled fraction is took into account in Idet, Qdet, Udet, and Vdet.
       pobs%xfreq = xfreq
       call raytrace_to_edge(pobs,grid,tau)
       wgt  = 1.0_wp/r2 * exp(-tau) * photon%wgt

       observer(i)%scatt(ixf,ix,iy) = observer(i)%scatt(ixf,ix,iy) + wgt * Idet
       observer(i)%I(ixf,ix,iy)     = observer(i)%I(ixf,ix,iy)     + wgt * Idet
       observer(i)%Q(ixf,ix,iy)     = observer(i)%Q(ixf,ix,iy)     + wgt * Qdet
       observer(i)%U(ixf,ix,iy)     = observer(i)%U(ixf,ix,iy)     + wgt * Udet
       observer(i)%V(ixf,ix,iy)     = observer(i)%V(ixf,ix,iy)     + wgt * Vdet
    endif
  enddo
  end subroutine peeling_resonance_stokes
  !--------------------------------------------------
  subroutine peeling_dust_nostokes(photon,grid)
  implicit none
  type(photon_type),   intent(in)    :: photon
  type(grid_type),     intent(in)    :: grid
  ! local variables
  real(kind=wp), parameter :: three_over_16pi = 3.0_wp/(16.0_wp*pi)
  type (photon_type) :: pobs
  real(kind=wp) :: r2,r,kx,ky,kz!,vdet(3)
  real(kind=wp) :: cosa,wgt,peel,tau
  real(kind=wp) :: xfreq_ref, u1
  integer       :: ix,iy,ixf
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

    !--- Transform the peeling-off vector to the observer's frame.
    kx = observer(i)%rmatrix(1,1)*pobs%kx + observer(i)%rmatrix(1,2)*pobs%ky + observer(i)%rmatrix(1,3)*pobs%kz
    ky = observer(i)%rmatrix(2,1)*pobs%kx + observer(i)%rmatrix(2,2)*pobs%ky + observer(i)%rmatrix(2,3)*pobs%kz
    kz = observer(i)%rmatrix(3,1)*pobs%kx + observer(i)%rmatrix(3,2)*pobs%ky + observer(i)%rmatrix(3,3)*pobs%kz

    !--- Location in TAN image.
    ix = floor(atan2(-kx,kz)*rad2deg/observer(i)%dxim+observer(i)%nxim/2.0_wp) + 1
    iy = floor(atan2(-ky,kz)*rad2deg/observer(i)%dyim+observer(i)%nyim/2.0_wp) + 1

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
    xfreq_ref = xfreq_ref * grid%Dfreq_ratio(photon%icell,photon%jcell,photon%kcell)
    ixf       = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1

    if (ix >= 1 .and. ix <= observer(i)%nxim .and. iy >= 1 .and. iy <= observer(i)%nyim &
        .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
       call raytrace_to_edge(pobs,grid,tau)
       cosa = photon%kx*pobs%kx+photon%ky*pobs%ky+photon%kz*pobs%kz
       peel = (1.0_wp - par%hgg**2)/((1.0_wp + par%hgg**2)-2.0_wp*par%hgg*cosa)**1.5_wp/fourpi
       !--- albedo was already multiplied before this routine is called.
       wgt  = peel/r2 * exp(-tau) * photon%wgt
       observer(i)%scatt(ixf,ix,iy) = observer(i)%scatt(ixf,ix,iy) + wgt
    endif
  enddo
  end subroutine peeling_dust_nostokes
  !--------------------------------------------------
  subroutine peeling_resonance_nostokes(photon,grid,xfreq_atom,vel_atom)
  use random
  implicit none
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  real(kind=wp),     intent(in) :: xfreq_atom
  real(kind=wp),     intent(in) :: vel_atom(3)
  ! local variables
  real(kind=wp), parameter :: three_over_16pi = 3.0_wp/(16.0_wp*pi)
  type (photon_type) :: pobs
  real(kind=wp) :: r2,r
  !real(kind=wp) :: vec(3)
  real(kind=wp) :: wgt,peel,tau
  real(kind=wp) :: cost,cost2,sint,cosp,sinp,rho1,rho
  real(kind=wp) :: kx,ky,kz
  real(kind=wp) :: xfreq_ref, xfreq, g_recoil, u1
  integer       :: ix,iy,ixf
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

    !----------------------------------------------------------------------------------
    !--- Transform the peeling-off vector to the observer's frame.
    kx = observer(i)%rmatrix(1,1)*pobs%kx + observer(i)%rmatrix(1,2)*pobs%ky + observer(i)%rmatrix(1,3)*pobs%kz
    ky = observer(i)%rmatrix(2,1)*pobs%kx + observer(i)%rmatrix(2,2)*pobs%ky + observer(i)%rmatrix(2,3)*pobs%kz
    kz = observer(i)%rmatrix(3,1)*pobs%kx + observer(i)%rmatrix(3,2)*pobs%ky + observer(i)%rmatrix(3,3)*pobs%kz

    !--- Location in TAN image.
    ix = floor(atan2(-kx,kz)*rad2deg/observer(i)%dxim+observer(i)%nxim/2.0_wp) + 1
    iy = floor(atan2(-ky,kz)*rad2deg/observer(i)%dyim+observer(i)%nyim/2.0_wp) + 1

    if (ix >= 1 .and. ix <= observer(i)%nxim .and. iy >= 1 .and. iy <= observer(i)%nyim) then
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
          g_recoil = g_recoil0 /(grid%Dfreq(photon%icell,photon%jcell,photon%kcell)*par%atom_no)
          xfreq    = xfreq - g_recoil * (1.0_wp - cost)
       endif

       !--- xfreq_ref = lab (observer) frame frequency of peel-off photon.
       u1 = grid%vfx(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kx + &
            grid%vfy(pobs%icell,pobs%jcell,pobs%kcell)*pobs%ky + &
            grid%vfz(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kz
       xfreq_ref = xfreq + u1

       !--- frequency bin
       xfreq_ref = xfreq_ref * grid%Dfreq_ratio(photon%icell,photon%jcell,photon%kcell)
       ixf       = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1

       !----------------------------------------------------------------------------------
       if (ixf >= 1 .and. ixf <= grid%nxfreq) then
          pobs%xfreq = xfreq
          call raytrace_to_edge(pobs,grid,tau)
          peel = three_over_16pi * (1.0_wp + cost2)
          wgt  = peel/r2 * exp(-tau) * photon%wgt
          observer(i)%scatt(ixf,ix,iy) = observer(i)%scatt(ixf,ix,iy) + wgt
       endif
    endif
  enddo
  end subroutine peeling_resonance_nostokes
  !--------------------------------------------------
  subroutine make_sightline_tau(grid)
  !--- calculate the optical depth along sight lines that are projected to a detector-plane pixel (2020/09/20).
  use mpi
  implicit none
  type(grid_type),  intent(in) :: grid
  !--- local variables
  type(photon_type) :: pobs
  real(kind=wp) :: delt(6),dist
  real(kind=wp) :: tau_HI,N_HI,tau_dust
  real(kind=wp) :: kx,ky,kz,kr,u1
  integer       :: ix,iy
  integer       :: loop,loop1,loop2,nsize,nsize2
  integer       :: ierr
  integer       :: i,jj,kk

  do i=1,par%nobs
    nsize = observer(i)%nxim * observer(i)%nyim
    call loop_divide(nsize,mpar%nproc,mpar%p_rank,loop1,loop2)
    do loop=loop1,loop2
       call array_2D_indices(observer(i)%nxim,observer(i)%nyim,loop,ix,iy)
       kx = tan((ix - (observer(i)%nxim+1.0_wp)/2.0_wp) * observer(i)%dxim/rad2deg)
       ky = tan((iy - (observer(i)%nyim+1.0_wp)/2.0_wp) * observer(i)%dyim/rad2deg)
       kz = -1.0_wp
       kr = sqrt(kx*kx + ky*ky + kz*kz)
       kx = kx/kr
       ky = ky/kr
       kz = kz/kr
       pobs%kx = observer(i)%rmatrix(1,1)*kx + observer(i)%rmatrix(2,1)*ky + observer(i)%rmatrix(3,1)*kz
       pobs%ky = observer(i)%rmatrix(1,2)*kx + observer(i)%rmatrix(2,2)*ky + observer(i)%rmatrix(3,2)*kz
       pobs%kz = observer(i)%rmatrix(1,3)*kx + observer(i)%rmatrix(2,3)*ky + observer(i)%rmatrix(3,3)*kz
       if (pobs%kx == 0.0_wp) then
          delt(1) = hugest
          delt(2) = hugest
       else
          delt(1) = (grid%xmax-observer(i)%x)/pobs%kx
          delt(2) = (grid%xmin-observer(i)%x)/pobs%kx
       endif
       if (pobs%ky == 0.0_wp) then
          delt(3) = hugest
          delt(4) = hugest
       else
          delt(3) = (grid%ymax-observer(i)%y)/pobs%ky
          delt(4) = (grid%ymin-observer(i)%y)/pobs%ky
       endif
       if (pobs%kz == 0.0_wp) then
          delt(5) = hugest
          delt(6) = hugest
       else
          delt(5) = (grid%zmax-observer(i)%z)/pobs%kz
          delt(6) = (grid%zmin-observer(i)%z)/pobs%kz
       endif

       !-- Find the closest boundary where the ray touches the grid system.
       !dist = hugest
       !do jj=1,6
       !   if (delt(jj) > 0.0_wp .and. delt(jj) < hugest) then
       !      pobs%x     = observer(i)%x + pobs%kx * delt(jj)
       !      pobs%y     = observer(i)%y + pobs%ky * delt(jj)
       !      pobs%z     = observer(i)%z + pobs%kz * delt(jj)
       !      pobs%icell = floor((pobs%x-grid%xmin)/grid%dx)+1
       !      pobs%jcell = floor((pobs%y-grid%ymin)/grid%dy)+1
       !      pobs%kcell = floor((pobs%z-grid%zmin)/grid%dz)+1
       !      if (pobs%icell >=0 .and. pobs%icell <= grid%nx+1 .and. &
       !          pobs%jcell >=0 .and. pobs%jcell <= grid%ny+1 .and. &
       !          pobs%kcell >=0 .and. pobs%kcell <= grid%nz+1) then
       !         if (delt(jj) < dist) dist = delt(jj)
       !      endif
       !   endif
       !enddo
       !-- Find the farthest boundary where the ray touches the grid system.
       !-- We measure optical depth from the distant universe toward Earth. (2020.10.20)
       dist = -999.9_wp
       do jj=1,6
          if (delt(jj) > 0.0_wp .and. delt(jj) < hugest) then
             pobs%x     = observer(i)%x + pobs%kx * delt(jj)
             pobs%y     = observer(i)%y + pobs%ky * delt(jj)
             pobs%z     = observer(i)%z + pobs%kz * delt(jj)
             pobs%icell = floor((pobs%x-grid%xmin)/grid%dx)+1
             pobs%jcell = floor((pobs%y-grid%ymin)/grid%dy)+1
             pobs%kcell = floor((pobs%z-grid%zmin)/grid%dz)+1
             if (pobs%icell >=0 .and. pobs%icell <= grid%nx+1 .and. &
                 pobs%jcell >=0 .and. pobs%jcell <= grid%ny+1 .and. &
                 pobs%kcell >=0 .and. pobs%kcell <= grid%nz+1) then
                if (delt(jj) > dist) dist = delt(jj)
             endif
          endif
       enddo

       !if (dist < hugest) then
       if (dist > 0.0_wp .and. dist < hugest) then
          pobs%x     = observer(i)%x + pobs%kx * dist
          pobs%y     = observer(i)%y + pobs%ky * dist
          pobs%z     = observer(i)%z + pobs%kz * dist
          pobs%icell = floor((pobs%x-grid%xmin)/grid%dx)+1
          pobs%jcell = floor((pobs%y-grid%ymin)/grid%dy)+1
          pobs%kcell = floor((pobs%z-grid%zmin)/grid%dz)+1

          !--- After finding the starting position of the ray, the direction vector should be revsersed. (2020.10.20)
          pobs%kx = -pobs%kx
          pobs%ky = -pobs%ky
          pobs%kz = -pobs%kz

          if (pobs%icell == 0) pobs%icell = 1
          if (pobs%jcell == 0) pobs%jcell = 1
          if (pobs%kcell == 0) pobs%kcell = 1
          if (pobs%icell == grid%nx+1) pobs%icell = grid%nx
          if (pobs%jcell == grid%ny+1) pobs%jcell = grid%ny
          if (pobs%kcell == grid%nz+1) pobs%kcell = grid%nz

          if (pobs%icell >=1 .and. pobs%icell <= grid%nx .and. &
              pobs%jcell >=1 .and. pobs%jcell <= grid%ny .and. &
              pobs%kcell >=1 .and. pobs%kcell <= grid%nz) then
             !--- Calculate the optical depth for the whole range of frequency (2020.10.19).
             do kk=1, observer(i)%nxfreq
                pobs%xfreq = grid%xfreq(kk)
                u1 = grid%vfx(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kx + &
                     grid%vfy(pobs%icell,pobs%jcell,pobs%kcell)*pobs%ky + &
                     grid%vfz(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kz
                pobs%xfreq = pobs%xfreq - u1

                !--- note that we are using raytrace_to_edge_tau_HI
                call raytrace_to_edge_tau_HI(pobs,grid,tau_HI)
                observer(i)%tau_HI(kk,ix,iy) = tau_HI
             enddo
             call raytrace_to_edge_column(pobs,grid,N_HI,tau_dust)
             observer(i)%N_HI(ix,iy)  = N_HI
             if (par%DGR > 0.0_wp) observer(i)%tau_dust(ix,iy) = tau_dust
          endif
       endif
    enddo

    nsize2 = observer(i)%nxim * observer(i)%nyim * observer(i)%nxfreq
    call reduce_mem(observer(i)%tau_HI)
    call reduce_mem(observer(i)%N_HI)
    if (par%DGR > 0.0_wp) call reduce_mem(observer(i)%tau_dust)
  enddo
  end subroutine make_sightline_tau
#endif
end module peelingoff_mod
