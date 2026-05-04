!-- Modification History
!   2018-02-03, Frequency binning for peeling-off image should be done in observer's frame (lab frame).
!--
module peelingoff_rect
  use define
  use mathlib
  use utility
  use line_mod
  use raytrace_clump_mod, only: raytrace_to_edge_clump_capped, tau_huge_clump
  !--- Clump-mode bulk velocity: when photon is inside a clump
  !    (par%use_clump_medium .and. photon%icell_clump > 0), grid%vfx/y/z = 0
  !    and the photon's xfreq is in the clump's rest frame. We must add the
  !    clump's bulk velocity to convert to the lab frame before binning.
  !    cl_v* are stored as v / cl_vtherm (dimensionless, see clump_mod.f90),
  !    so they can be used directly in u1 without further normalisation,
  !    matching the Cartesian grid%vfx convention. (Bug fix 2026-04-27:
  !    previously used grid%vfx=0, which caused peel-off spectra to be in
  !    the clump frame and disagree with the Jout spectrum, e.g., for
  !    rotating-halo or Hubble-flow clump runs.)
  use clump_mod, only: cl_vx, cl_vy, cl_vz
contains
  !--------------------------------------------------
  subroutine peeling_direct_outside(photon,grid)
  implicit none
  type(photon_type),   intent(in)    :: photon
  type(grid_type),     intent(in)    :: grid
  ! local variables
  type(photon_type) :: pobs
  real(kind=wp) :: r2,r,wgt,wgt0,tau
  real(kind=wp) :: kx,ky,kz
  real(kind=wp) :: xfreq_ref, u1, u2
  integer :: ix,iy,ixf
  integer :: i

  do i=1,par%nobs
    !--- For a plane-parallel illumination source, the direct and direc0 data are calcualted only when (obsx, obsy, obsz) = (0, 0, 1).
    !--- (2021.05.16) (2021.07.10)
    if ((trim(par%source_geometry) == 'plane_illumination') .and. &
        (observer(i)%x /= 0.0_wp .or. observer(i)%y /= 0.0_wp)) then
      cycle
    endif

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
     xfreq_ref = xfreq_ref * (grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)
     ixf       = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1

     if (ix >= 1 .and. ix <= observer(i)%nxim .and. iy >= 1 .and. iy <= observer(i)%nyim) then
        call peel_raytrace_to_edge(pobs, grid, tau)
        !-- bug-fixed (2022.04.26), direc was the same as direc0 when par%save_direc0 = .true.
        wgt0 = 1.0_wp/(fourpi*r2) * photon%wgt
        wgt  = exp(-tau)*wgt0
        !wgt  = exp(-tau)/(fourpi*r2) * photon%wgt

        !--- 2D image
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

        !--- 3D spectral image
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
  end subroutine peeling_direct_outside
  !--------------------------------------------------
  subroutine peeling_dust_stokes_outside(photon,grid)
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
  real(kind=wp) :: xx,yy,zz,rr,r_dot_k,det
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

    !---------------------------------------------------
    !-- Check the photon is blocked by the star or not.
    !-- Skip if the ray is blocked.
    if (trim(par%source_geometry) == 'stellar_illumination') then
       xx      = pobs%x
       yy      = pobs%y
       zz      = pobs%z + par%distance_star_to_planet
       rr      = sqrt(xx**2 + yy**2 + zz**2)
       r_dot_k = xx*pobs%kx + yy*pobs%ky + zz*pobs%kz
       det     = r_dot_k**2 - (rr**2 - par%stellar_radius**2)
       if (r_dot_k < 0.0_wp .and. det >= 0.0_wp) cycle
    endif
    !---------------------------------------------------

    !--- Transform the peeling-off vector to the observer's frame.
    kx = observer(i)%rmatrix(1,1)*pobs%kx + observer(i)%rmatrix(1,2)*pobs%ky + observer(i)%rmatrix(1,3)*pobs%kz
    ky = observer(i)%rmatrix(2,1)*pobs%kx + observer(i)%rmatrix(2,2)*pobs%ky + observer(i)%rmatrix(2,3)*pobs%kz
    kz = observer(i)%rmatrix(3,1)*pobs%kx + observer(i)%rmatrix(3,2)*pobs%ky + observer(i)%rmatrix(3,3)*pobs%kz

    !--- Location in TAN image.
    ix = floor(atan2(-kx,kz)*rad2deg/observer(i)%dxim+observer(i)%nxim/2.0_wp) + 1
    iy = floor(atan2(-ky,kz)*rad2deg/observer(i)%dyim+observer(i)%nyim/2.0_wp) + 1

    !--- frequency for spectral binning in the detector plane.
    !--- xfreq_ref = lab (observer) frame frequency of peel-off photon (2020.08.28, bug-fixed)
    if (par%use_clump_medium .and. photon%icell_clump > 0) then
       u1 = cl_vx(photon%icell_clump)*pobs%kx + &
            cl_vy(photon%icell_clump)*pobs%ky + &
            cl_vz(photon%icell_clump)*pobs%kz
    else
       u1 = grid%vfx(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kx + &
            grid%vfy(pobs%icell,pobs%jcell,pobs%kcell)*pobs%ky + &
            grid%vfz(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kz
    endif
    xfreq_ref = photon%xfreq + u1

    !--- Note that frequency for the calculation of optical depth along the peel-off direction is
    !--- pobs%xfreq = photon%xfreq because dust grains are assumed to have no thermal motion.
    !--- comment added on 2020.08.23.

    !--- frequency bin
    xfreq_ref = xfreq_ref * (grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)
    ixf       = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq) + 1

    if (ix >= 1 .and. ix <= observer(i)%nxim .and. iy >= 1 .and. iy <= observer(i)%nyim) then
       !--- Calculate polar scattering angle toward the observer.
       cost = photon%kx * pobs%kx + photon%ky * pobs%ky + photon%kz * pobs%kz
       sint = sqrt(1.0_wp - cost*cost)

       !--- Calculate azimuthal scattering angle toward the observer.
       if (sint == 0.0_wp) then
          cosp  = 1.0_wp
          sinp  = 0.0_wp
          cos2p = 1.0_wp
          sin2p = 0.0_wp
       else
          cosp  = (pobs%kx * photon%mx + pobs%ky * photon%my + pobs%kz * photon%mz)/sint
          sinp  = (pobs%kx * photon%nx + pobs%ky * photon%ny + pobs%kz * photon%nz)/sint
          cos2p = 2.0_wp*cosp*cosp - 1.0_wp
          sin2p = 2.0_wp*cosp*sinp
       endif

       !--- Calculate the reference normal vector for peeling-off.
       pobs%nx = -sinp * photon%mx + cosp * photon%nx
       pobs%ny = -sinp * photon%my + cosp * photon%ny
       pobs%nz = -sinp * photon%mz + cosp * photon%nz

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
       call peel_raytrace_to_edge(pobs, grid, tau)
       wgt  = 1.0_wp/r2 * exp(-tau) * photon%wgt

       !--- 2D image
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

       !--- 3D spectral image
       if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt(ixf,ix,iy) = observer(i)%scatt(ixf,ix,iy) + wgt * Idet
          if (par%save_dust_scattered) then
             !$OMP ATOMIC UPDATE
             observer(i)%scatt_dust(ixf,ix,iy) = observer(i)%scatt_dust(ixf,ix,iy) + wgt * Idet
          endif
          !$OMP ATOMIC UPDATE
          observer(i)%I(ixf,ix,iy)     = observer(i)%I(ixf,ix,iy)     + wgt * Idet
          !$OMP ATOMIC UPDATE
          observer(i)%Q(ixf,ix,iy)     = observer(i)%Q(ixf,ix,iy)     + wgt * Qdet
          !$OMP ATOMIC UPDATE
          observer(i)%U(ixf,ix,iy)     = observer(i)%U(ixf,ix,iy)     + wgt * Udet
          !$OMP ATOMIC UPDATE
          observer(i)%V(ixf,ix,iy)     = observer(i)%V(ixf,ix,iy)     + wgt * Vdet
       endif
    endif
  enddo
  end subroutine peeling_dust_stokes_outside
  !--------------------------------------------------
  !--- xfreq_atom = frequency in the scattering atom's rest frame.
  !--- vel_atom   = velocity of the scattering atom in the referenc frame carried with photon.
  subroutine peeling_resonance_stokes_outside(photon,grid,xfreq_atom,vel_atom)
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
  real(kind=wp) :: xx,yy,zz,rr,r_dot_k,det
  integer :: ix,iy,ixf
  integer :: i

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

    !---------------------------------------------------
    !-- Check the photon is blocked by the star or not.
    !-- Skip if the ray is blocked.
    if (trim(par%source_geometry) == 'stellar_illumination') then
       xx      = pobs%x
       yy      = pobs%y
       zz      = pobs%z + par%distance_star_to_planet
       rr      = sqrt(xx**2 + yy**2 + zz**2)
       r_dot_k = xx*pobs%kx + yy*pobs%ky + zz*pobs%kz
       det     = r_dot_k**2 - (rr**2 - par%stellar_radius**2)
       if (r_dot_k < 0.0_wp .and. det >= 0.0_wp) cycle
    endif
    !---------------------------------------------------

    !--- Transform the peeling-off vector to the observer's frame.
    kx = observer(i)%rmatrix(1,1)*pobs%kx + observer(i)%rmatrix(1,2)*pobs%ky + observer(i)%rmatrix(1,3)*pobs%kz
    ky = observer(i)%rmatrix(2,1)*pobs%kx + observer(i)%rmatrix(2,2)*pobs%ky + observer(i)%rmatrix(2,3)*pobs%kz
    kz = observer(i)%rmatrix(3,1)*pobs%kx + observer(i)%rmatrix(3,2)*pobs%ky + observer(i)%rmatrix(3,3)*pobs%kz

    !--- Location in TAN image.
    ix = floor(atan2(-kx,kz)*rad2deg/observer(i)%dxim+observer(i)%nxim/2.0_wp) + 1
    iy = floor(atan2(-ky,kz)*rad2deg/observer(i)%dyim+observer(i)%nyim/2.0_wp) + 1

    if (ix >= 1 .and. ix <= observer(i)%nxim .and. iy >= 1 .and. iy <= observer(i)%nyim) then
       !--- Calculate polar scattering angle toward the observer.
       cost = photon%kx * pobs%kx + photon%ky * pobs%ky + photon%kz * pobs%kz
       sint = sqrt(1.0_wp - cost*cost)

       !--- Calculate azimuthal scattering angle toward the observer.
       if (sint == 0.0_wp) then
          cosp  = 1.0_wp
          sinp  = 0.0_wp
          cos2p = 1.0_wp
          sin2p = 0.0_wp
       else
          cosp  = (pobs%kx * photon%mx + pobs%ky * photon%my + pobs%kz * photon%mz)/sint
          sinp  = (pobs%kx * photon%nx + pobs%ky * photon%ny + pobs%kz * photon%nz)/sint
          cos2p = 2.0_wp*cosp*cosp - 1.0_wp
          sin2p = 2.0_wp*cosp*sinp
       endif

       !--- Calculate the reference normal vector for peeling-off.
       pobs%nx = -sinp * photon%mx + cosp * photon%nx
       pobs%ny = -sinp * photon%my + cosp * photon%ny
       pobs%nz = -sinp * photon%mz + cosp * photon%nz

       !--- Calculate Stokes parameters. (S11, S12, S33, and S44)
       !--- Note S44 /= S33 and S34 = 0.
       cost2 = cost**2
       S22   = three_over_four * photon%E1 * (cost2 + 1.0_wp)
       S11   = S22 + photon%E2
       S12   = three_over_four * photon%E1 * (cost2 - 1.0_wp)
       S33   = three_over_two  * photon%E1 * cost
       S44   = three_over_two  * photon%E3 * cost

       !--- frequency in the fluid frame (comoving frame).
       xfreq = xfreq_atom + (vel_atom(1)*cosp + vel_atom(2)*sinp)*sint + vel_atom(3)*cost
       if (par%recoil) then
          if (line%line_type == 7 .and. line%selected_species_HD == 2) then
             g_recoil = line%g_recoil0_D / grid%Dfreq(photon%icell,photon%jcell,photon%kcell)
          else
             g_recoil = line%g_recoil0   / grid%Dfreq(photon%icell,photon%jcell,photon%kcell)
          endif
          xfreq    = xfreq - g_recoil * (1.0_wp - cost)
       endif

       !--- xfreq_ref = lab (observer) frame frequency of peel-off photon.
       if (par%use_clump_medium .and. photon%icell_clump > 0) then
          u1 = cl_vx(photon%icell_clump)*pobs%kx + &
               cl_vy(photon%icell_clump)*pobs%ky + &
               cl_vz(photon%icell_clump)*pobs%kz
       else
          u1 = grid%vfx(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kx + &
               grid%vfy(pobs%icell,pobs%jcell,pobs%kcell)*pobs%ky + &
               grid%vfz(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kz
       endif
       xfreq_ref = xfreq + u1

       !--- frequency bin
       xfreq_ref = xfreq_ref * (grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)
       ixf       = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq) + 1

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
       call peel_raytrace_to_edge(pobs, grid, tau)
       wgt  = 1.0_wp/r2 * exp(-tau) * photon%wgt

       !--- 2D image
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

       !--- 3D spectral image
       if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt(ixf,ix,iy) = observer(i)%scatt(ixf,ix,iy) + wgt * Idet
          !$OMP ATOMIC UPDATE
          observer(i)%I(ixf,ix,iy)     = observer(i)%I(ixf,ix,iy)     + wgt * Idet
          !$OMP ATOMIC UPDATE
          observer(i)%Q(ixf,ix,iy)     = observer(i)%Q(ixf,ix,iy)     + wgt * Qdet
          !$OMP ATOMIC UPDATE
          observer(i)%U(ixf,ix,iy)     = observer(i)%U(ixf,ix,iy)     + wgt * Udet
          !$OMP ATOMIC UPDATE
          observer(i)%V(ixf,ix,iy)     = observer(i)%V(ixf,ix,iy)     + wgt * Vdet
       endif
    endif
  enddo
  end subroutine peeling_resonance_stokes_outside
  !--------------------------------------------------
  subroutine peeling_dust_nostokes_outside(photon,grid)
  implicit none
  type(photon_type),   intent(in)    :: photon
  type(grid_type),     intent(in)    :: grid
  ! local variables
  type (photon_type) :: pobs
  real(kind=wp) :: r2,r,kx,ky,kz
  real(kind=wp) :: cosa,wgt,peel,tau
  real(kind=wp) :: xfreq_ref, u1
  real(kind=wp) :: xx,yy,zz,rr,r_dot_k,det
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

    !---------------------------------------------------
    !-- Check the photon is blocked by the star or not.
    !-- Skip if the ray is blocked.
    if (trim(par%source_geometry) == 'stellar_illumination') then
       xx      = pobs%x
       yy      = pobs%y
       zz      = pobs%z + par%distance_star_to_planet
       rr      = sqrt(xx**2 + yy**2 + zz**2)
       r_dot_k = xx*pobs%kx + yy*pobs%ky + zz*pobs%kz
       det     = r_dot_k**2 - (rr**2 - par%stellar_radius**2)
       if (r_dot_k < 0.0_wp .and. det >= 0.0_wp) cycle
    endif
    !---------------------------------------------------

    !--- Transform the peeling-off vector to the observer's frame.
    kx = observer(i)%rmatrix(1,1)*pobs%kx + observer(i)%rmatrix(1,2)*pobs%ky + observer(i)%rmatrix(1,3)*pobs%kz
    ky = observer(i)%rmatrix(2,1)*pobs%kx + observer(i)%rmatrix(2,2)*pobs%ky + observer(i)%rmatrix(2,3)*pobs%kz
    kz = observer(i)%rmatrix(3,1)*pobs%kx + observer(i)%rmatrix(3,2)*pobs%ky + observer(i)%rmatrix(3,3)*pobs%kz

    !--- Location in TAN image.
    ix = floor(atan2(-kx,kz)*rad2deg/observer(i)%dxim+observer(i)%nxim/2.0_wp) + 1
    iy = floor(atan2(-ky,kz)*rad2deg/observer(i)%dyim+observer(i)%nyim/2.0_wp) + 1

    if (ix >= 1 .and. ix <= observer(i)%nxim .and. iy >= 1 .and. iy <= observer(i)%nyim) then
       !--- frequency for spectral binning.
       !--- xfreq_ref = lab (observer) frame frequency of peel-off photon (2020.08.28, bug-fixed)
       if (par%use_clump_medium .and. photon%icell_clump > 0) then
          u1 = cl_vx(photon%icell_clump)*pobs%kx + &
               cl_vy(photon%icell_clump)*pobs%ky + &
               cl_vz(photon%icell_clump)*pobs%kz
       else
          u1 = grid%vfx(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kx + &
               grid%vfy(pobs%icell,pobs%jcell,pobs%kcell)*pobs%ky + &
               grid%vfz(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kz
       endif
       xfreq_ref = photon%xfreq + u1

       !--- Note that frequency for the calculation of optical depth along the peel-off direction is
       !--- pobs%xfreq = photon%xfreq because dust grains are assumed to have no thermal motion.
       !--- comment added on 2020.08.23.

       !--- frequency bin
       xfreq_ref = xfreq_ref * (grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)
       ixf       = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1

       call peel_raytrace_to_edge(pobs, grid, tau)
       cosa = photon%kx*pobs%kx+photon%ky*pobs%ky+photon%kz*pobs%kz
       peel = (1.0_wp - par%hgg**2)/((1.0_wp + par%hgg**2)-2.0_wp*par%hgg*cosa)**1.5_wp/fourpi
       !--- albedo was already multiplied before this routine is called.
       wgt  = peel/r2 * exp(-tau) * photon%wgt

       !--- 2D image
       if (par%save_peeloff_2D) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt_2D(ix,iy) = observer(i)%scatt_2D(ix,iy) + wgt
       endif

       !--- 3D spectral image
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
  end subroutine peeling_dust_nostokes_outside
  !--------------------------------------------------
  subroutine peeling_resonance_nostokes_outside(photon,grid,xfreq_atom,vel_atom)
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

    !---------------------------------------------------
    !-- Check the photon is blocked by the star or not.
    !-- Skip if the ray is blocked.
    if (trim(par%source_geometry) == 'stellar_illumination') then
       xx      = pobs%x
       yy      = pobs%y
       zz      = pobs%z + par%distance_star_to_planet
       rr      = sqrt(xx**2 + yy**2 + zz**2)
       r_dot_k = xx*pobs%kx + yy*pobs%ky + zz*pobs%kz
       det     = r_dot_k**2 - (rr**2 - par%stellar_radius**2)
       if (r_dot_k < 0.0_wp .and. det >= 0.0_wp) cycle
    endif
    !---------------------------------------------------

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
          if (line%line_type == 7 .and. line%selected_species_HD == 2) then
             g_recoil = line%g_recoil0_D / grid%Dfreq(photon%icell,photon%jcell,photon%kcell)
          else
             g_recoil = line%g_recoil0   / grid%Dfreq(photon%icell,photon%jcell,photon%kcell)
          endif
          xfreq    = xfreq - g_recoil * (1.0_wp - cost)
       endif

       !--- xfreq_ref = lab (observer) frame frequency of peel-off photon.
       if (par%use_clump_medium .and. photon%icell_clump > 0) then
          u1 = cl_vx(photon%icell_clump)*pobs%kx + &
               cl_vy(photon%icell_clump)*pobs%ky + &
               cl_vz(photon%icell_clump)*pobs%kz
       else
          u1 = grid%vfx(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kx + &
               grid%vfy(pobs%icell,pobs%jcell,pobs%kcell)*pobs%ky + &
               grid%vfz(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kz
       endif
       xfreq_ref = xfreq + u1

       !--- frequency bin
       xfreq_ref = xfreq_ref * (grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)
       ixf       = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1

       !----------------------------------------------------------------------------------
       pobs%xfreq = xfreq
       call peel_raytrace_to_edge(pobs, grid, tau)
       peel = three_over_four * photon%E1 * (cost2 + 1.0_wp) + photon%E2
       wgt  = peel/(fourpi*r2) * exp(-tau) * photon%wgt

       !--- 2D image
       if (par%save_peeloff_2D) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt_2D(ix,iy) = observer(i)%scatt_2D(ix,iy) + wgt
       endif

       !--- 3D spectral image
       if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
          !$OMP ATOMIC UPDATE
          observer(i)%scatt(ixf,ix,iy) = observer(i)%scatt(ixf,ix,iy) + wgt
       endif
    endif
  enddo
  end subroutine peeling_resonance_nostokes_outside
  !--------------------------------------------------
  subroutine peel_raytrace_to_edge(pobs, grid, tau)
  !---------------------------------------------------------------------------
  ! Mode-aware tau-to-edge raytrace used by every peel-off routine in this
  ! module.  In Cartesian / AMR mode the raytrace_to_edge procedure pointer
  ! (assigned by setup_procedure to one of the _car / _amr variants) handles
  ! the integral correctly.  In clump mode that pointer is set to
  ! raytrace_to_edge_clump, but the procedure-pointer interface has no
  ! tau_max argument -- so the integration runs uncapped and wastes CPU on
  ! photons whose contribution exp(-tau) has long underflowed.  Bypass the
  ! pointer in clump mode and call raytrace_to_edge_clump_capped directly.
  ! (The cap is correct here because peel-off only needs tau accurately up
  ! to the exp(-tau) underflow point; the cap is ABSENT in
  ! raytrace_to_edge_tau_gas_clump where the sight-line FITS output
  ! requires the full integrated tau.)
  !---------------------------------------------------------------------------
  type(photon_type), intent(in)  :: pobs
  type(grid_type),   intent(in)  :: grid
  real(kind=wp),     intent(out) :: tau
  if (par%use_clump_medium) then
     call raytrace_to_edge_clump_capped(pobs, grid, tau, tau_max=tau_huge_clump)
  else
     call raytrace_to_edge(pobs, grid, tau)
  end if
  end subroutine peel_raytrace_to_edge
  !--------------------------------------------------
end module peelingoff_rect
