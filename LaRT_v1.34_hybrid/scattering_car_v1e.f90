module scatter_mod
  use define
  use random
  use voigt_mod
  use peelingoff_mod
  use mathlib
  public
contains
  !=====================================
  subroutine scattering(photon,grid)
  implicit none
  ! Calculate new direction cosines and wavelength after scattering
  ! 2020/09/27, bug-fixed. Fine Structure Splitting was not taken into account.
  ! 2016/10/26, updated to use the acceleration scheme of Ahn et al. (2002) and Dijkstra et al. (2006)
  !             Laursen et al. (2009)
  ! Written by Kwang-il Seon, 2010/10/02
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid
  ! local variables
  real(kind=wp) :: p_dust
  integer :: i1,i2,i3
#ifdef FINE_STRUCTURE
  real(kind=wp) :: DnuHK, pH, pK
  real(kind=wp), parameter :: one_over_three  = 1.0_wp/3.0_wp, &
                              two_over_three  = 2.0_wp/3.0_wp
#endif
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt

  if (par%DGR > 0.0_wp) then
     i1 = photon%icell
     i2 = photon%jcell
     i3 = photon%kcell
#ifdef FINE_STRUCTURE
     !--- bug-fixed, 2020.09.27. The FS splitting effect had not been considered.
     DnuHK  = grid%DnuHK_ref_half/grid%Dfreq(i1,i2,i3)
     pH     = voigt(photon%xfreq+DnuHK, grid%voigt_a(i1,i2,i3)) * one_over_three
     pK     = voigt(photon%xfreq-DnuHK, grid%voigt_a(i1,i2,i3)) * two_over_three
     p_dust = grid%rhokapD(i1,i2,i3)/(grid%rhokap(i1,i2,i3) *(pH + pK) + grid%rhokapD(i1,i2,i3))
#else
     p_dust = grid%rhokapD(i1,i2,i3)/(grid%rhokap(i1,i2,i3)*voigt(photon%xfreq,grid%voigt_a(i1,i2,i3)) + grid%rhokapD(i1,i2,i3))
#endif
     if (rand_number() <= p_dust) then
        call scatter_dust(photon,grid)
     else
        call scatter_resonance(photon,grid)
     endif
  else
     call scatter_resonance(photon,grid)
  endif
  return
  end subroutine scattering
  !=====================================
  subroutine scatter_dust_stokes(photon,grid)
    use random
    implicit none
    ! Calculate new direction cosines
    ! Written by Kwang-il Seon, 2017/09/04
    type(photon_type), intent(inout) :: photon
    type(grid_type),   intent(inout) :: grid

    ! local variables
    real(kind=wp) :: xfreq_ref, uu1
    real(kind=wp) :: cost,sint,phi,phi1,cosp,sinp,cos2p,sin2p
    real(kind=wp) :: px,py,pz
    real(kind=wp) :: Prand, Pcomp, QoverI, UoverI, S12overS11
    real(kind=wp) :: S11,S12,S33,S34
    real(kind=wp) :: Q0,U0,I1,Q1,U1,V1
    integer :: ix

    photon%nscatt_dust = photon%nscatt_dust + 1

    !--- Add absorbed portion to Jabs.
    if (.not. par%use_reduced_wgt) then
       if (rand_number() > par%albedo) then
         ! absorbed by dust.
         !2018-02-04, Jabs is now calculated in lab frame.
         if (par%save_Jabs) then
            uu1 = grid%vfx(photon%icell,photon%jcell,photon%kcell)*photon%kx + &
                  grid%vfy(photon%icell,photon%jcell,photon%kcell)*photon%ky + &
                  grid%vfz(photon%icell,photon%jcell,photon%kcell)*photon%kz
            xfreq_ref = (photon%xfreq + uu1)*(grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)
            ix        = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1
            if (ix >= 1 .and. ix <= grid%nxfreq) then
               !$OMP ATOMIC UPDATE
               grid%Jabs(ix) = grid%Jabs(ix) + photon%wgt
            endif
         endif
         photon%inside = .false.
         if (par%save_all_photons) then
            !--- xfreq_ref is required to record the photon frequency when it is absorbed, measured in the lab frame.
            !--- (2020.09.27)
            if (.not. par%save_Jabs) then
               uu1 = grid%vfx(photon%icell,photon%jcell,photon%kcell)*photon%kx + &
                     grid%vfy(photon%icell,photon%jcell,photon%kcell)*photon%ky + &
                     grid%vfz(photon%icell,photon%jcell,photon%kcell)*photon%kz
               xfreq_ref = (photon%xfreq + uu1)*(grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)
            endif
            photon%xfreq_ref = xfreq_ref
            photon%wgt       = 0.0_wp
         endif
         return
       endif
    else
       !2018-02-04, Jabs is now calculated in lab frame.
       if (par%save_Jabs) then
          uu1 = grid%vfx(photon%icell,photon%jcell,photon%kcell)*photon%kx + &
                grid%vfy(photon%icell,photon%jcell,photon%kcell)*photon%ky + &
                grid%vfz(photon%icell,photon%jcell,photon%kcell)*photon%kz
          xfreq_ref = (photon%xfreq + uu1)*(grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)
          ix        = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1
          if (ix >= 1 .and. ix <= grid%nxfreq) then
             !$OMP ATOMIC UPDATE
             grid%Jabs(ix) = grid%Jabs(ix) + photon%wgt * (1.0_wp - par%albedo)
          endif
       endif
       photon%wgt = photon%wgt * par%albedo
    endif
    !--- Peeling-off
    if (par%save_peeloff) call peeling_dust_stokes(photon,grid)

    !--- Select a new cos(theta) from the numerical table using the alias method. (2021.08.30)
    cost  = rand_alias_linear(scatt_mat%phase_PDF, scatt_mat%alias, scatt_mat%coss, scatt_mat%S11)
    sint  = sqrt(1.0d0-cost**2)

    !--- Calculate Scattering Matrix.
    call interp_eq(scatt_mat%coss,scatt_mat%S11,cost,S11)
    call interp_eq(scatt_mat%coss,scatt_mat%S12,cost,S12)
    call interp_eq(scatt_mat%coss,scatt_mat%S33,cost,S33)
    call interp_eq(scatt_mat%coss,scatt_mat%S34,cost,S34)
    S12overS11 = S12/S11

    !--- Select a new phi using a rejection method.
    do while(.true.)
       phi    = twopi*rand_number()
       phi1   = 2.0_wp * phi
       ! Note: photon%Q = photon%Q/photon%I and photon%U = photon%U/photon%I
       Prand  = (1.0_wp + abs(S12overS11) * sqrt(photon%Q**2 + photon%U**2))*rand_number()
       Pcomp  =  1.0_wp + S12overS11 * (photon%Q*cos(phi1) + photon%U*sin(phi1))
       if (Prand <= Pcomp) exit
    enddo
    cosp = cos(phi)
    sinp = sin(phi)

    ! Which method is faster?
    !cos2p = cos(phi1)
    !sin2p = sin(phi1)
    cos2p = 2.0_wp*cosp*cosp - 1.0_wp
    sin2p = 2.0_wp*sinp*cosp

    !--- Calculate new Stokes parameters.
    Q0 =  cos2p*photon%Q + sin2p*photon%U
    U0 = -sin2p*photon%Q + cos2p*photon%U
    I1 =  S11 + S12*Q0
    Q1 =  S12 + S11*Q0
    U1 =  S33*U0 + S34*photon%V
    V1 = -S34*U0 + S33*photon%V

    photon%I  = 1.0_wp
    photon%Q  = Q1/I1
    photon%U  = U1/I1
    photon%V  = V1/I1

    !--- Calculate new reference vectors and propagation vector.
    !--- Do not change the calculation order.
    !--- Note (m), (n), (k) forms a right-handed triad.
    !--- Not sure the normalization of the reference vectors is really required.
    px = cosp * photon%mx + sinp * photon%nx
    py = cosp * photon%my + sinp * photon%ny
    pz = cosp * photon%mz + sinp * photon%nz

    photon%nx = cosp * photon%nx - sinp * photon%mx
    photon%ny = cosp * photon%ny - sinp * photon%my
    photon%nz = cosp * photon%nz - sinp * photon%mz

    photon%mx = cost * px - sint * photon%kx
    photon%my = cost * py - sint * photon%ky
    photon%mz = cost * pz - sint * photon%kz

    photon%kx = sint * px + cost * photon%kx
    photon%ky = sint * py + cost * photon%ky
    photon%kz = sint * pz + cost * photon%kz
  end subroutine scatter_dust_stokes
  !=====================================
  subroutine scatter_resonance_stokes(photon,grid)
    implicit none
    ! Written by Kwang-il Seon, 2017/09/04
    type(photon_type), intent(inout) :: photon
    type(grid_type),   intent(inout) :: grid
    ! local variables
    real(kind=wp) :: cost,cost2,sint,phi,phi1,phi2,cosp,sinp,cos2p,sin2p
    real(kind=wp) :: px,py,pz
    real(kind=wp) :: Prand, Pcomp, QoverI, UoverI, S12overS11
    real(kind=wp) :: S11,S12,S22,S33,S44
    real(kind=wp) :: Q0,U0,I1,Q1,U1,V1
    real(kind=wp) :: DnuHK,xfreq_atom,qK,qH,pH,pK
    real(kind=wp) :: ux,uy,uxy,uz,E1,g_recoil
    real(kind=wp), parameter :: three_over_four = 3.0_wp/4.0_wp, &
                                three_over_two  = 3.0_wp/2.0_wp, &
                                one_over_three  = 1.0_wp/3.0_wp, &
                                two_over_three  = 2.0_wp/3.0_wp, &
                                one_over_sqrt2  = 1.0_wp/sqrt(2.0_wp)
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt

    photon%nscatt_HI = photon%nscatt_HI + 1
#ifdef CALCP
    call add_to_Pa(photon,grid,photon%icell,photon%jcell,photon%kcell)
#endif

#ifdef FINE_STRUCTURE
    !--- Select an atom by which the photon is scattered.
    DnuHK = grid%DnuHK_ref_half/grid%Dfreq(photon%icell,photon%jcell,photon%kcell)
    pH    = voigt(photon%xfreq+DnuHK, grid%voigt_a(photon%icell,photon%jcell,photon%kcell)) * one_over_three
    pK    = voigt(photon%xfreq-DnuHK, grid%voigt_a(photon%icell,photon%jcell,photon%kcell)) * two_over_three
    pH    = pH/(pH + pK)
    if (rand_number() < pH) then
       !--- scattered by an atom at H state (2P1/2).
       uz = rand_resonance_vz(photon%xfreq + DnuHK, grid%voigt_a(photon%icell,photon%jcell,photon%kcell))
    else
       !--- scattered by an atom at K state (2P3/2).
       uz = rand_resonance_vz(photon%xfreq - DnuHK, grid%voigt_a(photon%icell,photon%jcell,photon%kcell))
    endif

    !--- xfreq_atom = frequency in the atom's rest frame.
    xfreq_atom = photon%xfreq - uz

    !--- Select a new scattering angle (cos(theta)).
    !--- E1 is a function of frequency in the atom's rest frame.
    !--- We are assuming nu_H = nu_K. The difference make no practical difference.
    !--- qH = nu - nu_H, qK = nu - nu_K, nu_K = nu_H + DnuHK_Hz
    !--- qH = xfreq_atom, qK = qH - DnuHK
    qH    = xfreq_atom + DnuHK
    qK    = xfreq_atom - DnuHK
    E1    = (2.0_wp*qK*qH + qH**2)/(qK**2 + 2.0_wp*qH**2)
    cost  = rand_resonance(E1)
    cost2 = cost**2
    sint  = sqrt(1.0_wp-cost2)

    !--- Calculate Scattering Matrix. (S11, S12, S33, and S44)
    !--- Note S44 /= S33 and S34 = 0.
    !--- E2 = 1-E1 and E3 = (E1+2)/3. (3/2) E3 = (E1+2)/2.
    S22 = three_over_four * E1 * (cost2 + 1.0_wp)
    S11 = S22 + (1.0_wp - E1)
    S12 = three_over_four * E1 * (cost2 - 1.0_wp)
    S33 = three_over_two  * E1 * cost
    S44 = (E1 + 2.0_wp)/2.0_wp * cost
    S12overS11 = S12/S11
#else
    !--- Select an atom by which the photon is scattered.
    uz = rand_resonance_vz(photon%xfreq, grid%voigt_a(photon%icell,photon%jcell,photon%kcell))

    !--- xfreq_atom = frequency in the atom's rest frame.
    xfreq_atom = photon%xfreq - uz

    !--- Select a new scattering angle (cos(theta)).
    cost  = rand_rayleigh()
    cost2 = cost**2
    sint  = sqrt(1.0_wp-cost2)

    !--- Calculate Scattering Matrix. (S11, S12, S33, and S44)
    S11 = three_over_four * (cost2 + 1.0_wp)
    S22 = S11
    S12 = three_over_four * (cost2 - 1.0_wp)
    S33 = three_over_two  * cost
    S44 = S33
    S12overS11 = S12/S11
#endif

    !--- Select a new scattering azimuthal angle (phi) using a rejection method.
    do while(.true.)
       phi    = twopi*rand_number()
       phi1   = 2.0_wp * phi
       ! Note: photon%Q = photon%Q/photon%I and photon%U = photon%U/photon%I
       Prand  = (1.0_wp + abs(S12overS11) * sqrt(photon%Q**2 + photon%U**2))*rand_number()
       Pcomp  =  1.0_wp + S12overS11 * (photon%Q*cos(phi1) + photon%U*sin(phi1))
       if (Prand <= Pcomp) exit
    enddo
    cosp = cos(phi)
    sinp = sin(phi)

    !--- Update photon's frequency after scattering.
    ! 2016-10-26, K.-I. Seon
    ! Note that the acceleration scheme uses a cut-off exponential distribution to obtain velocity amplitude, and
    ! divide it into x- and y- components, whereas the original scheme uses a gaussian distribution to obtain velocity components.
    ! This causes the difference in 1/sqrt(2) factor.
    ! new propagation direction: k' = (sint*cosp) e_x + (sint*sinp) e_y + cost e_z
    ! atom velociy             : u = ux e_x + uy e_y + uz e_z
    ! dot product              : u.k' = (ux*cosp + uy*sinp)sint + uz*cost
    if (par%core_skip .and. abs(photon%xfreq) < grid%xcrit) then
       phi2 = twopi * rand_number()
       uxy  = sqrt(grid%xcrit2 - log(rand_number()))
       ux   = uxy * cos(phi2)
       uy   = uxy * sin(phi2)
       photon%xfreq = xfreq_atom + uz*cost + (ux*cosp + uy*sinp)*sint
    else
       ux   = rand_gauss()*one_over_sqrt2
       uy   = rand_gauss()*one_over_sqrt2
       photon%xfreq = xfreq_atom + uz*cost + (ux*cosp + uy*sinp)*sint
    endif

    if (par%recoil) then
       g_recoil     = g_recoil0 /grid%Dfreq(photon%icell,photon%jcell,photon%kcell)
       photon%xfreq = photon%xfreq - g_recoil * (1.0_wp - cost)
    endif

    ! Peeling-off routine should be called before updating photon's triad vectors and Stokes parameters.
    ! But, after the selection of scattering atom.
    if (par%save_peeloff) call peeling_resonance_stokes(photon,grid,xfreq_atom,[ux,uy,uz])

    !--- Which method is faster?
    !cos2p = cos(phi1)
    !sin2p = sin(phi1)
    cos2p = 2.0_wp*cosp*cosp - 1.0_wp
    sin2p = 2.0_wp*sinp*cosp

    !--- Calculate new Stokes parameters.
    Q0 =  cos2p*photon%Q + sin2p*photon%U
    U0 = -sin2p*photon%Q + cos2p*photon%U
    I1 =  S11 + S12*Q0
    Q1 =  S12 + S22*Q0
    U1 =  S33*U0
    V1 =  S44*photon%V

    photon%I  = 1.0_wp
    photon%Q  = Q1/I1
    photon%U  = U1/I1
    photon%V  = V1/I1

    !--- Calculate new reference vectors and propagation vector.
    !--- Do not change the calculation order.
    !--- Note (m), (n), (k) forms a right-handed triad.
    px = cosp * photon%mx + sinp * photon%nx
    py = cosp * photon%my + sinp * photon%ny
    pz = cosp * photon%mz + sinp * photon%nz

    photon%nx = cosp * photon%nx - sinp * photon%mx
    photon%ny = cosp * photon%ny - sinp * photon%my
    photon%nz = cosp * photon%nz - sinp * photon%mz

    photon%mx = cost * px - sint * photon%kx
    photon%my = cost * py - sint * photon%ky
    photon%mz = cost * pz - sint * photon%kz

    photon%kx = sint * px + cost * photon%kx
    photon%ky = sint * py + cost * photon%ky
    photon%kz = sint * pz + cost * photon%kz
    return
  end subroutine scatter_resonance_stokes
  !=====================================
  subroutine scatter_dust_nostokes(photon,grid)
  implicit none
  ! Calculate new direction cosines and wavelength after scattering
  ! Written by Kwang-il Seon, 2010/10/02
  ! 2016/10/26, updated to use the acceleration scheme of Ahn et al. (2002) and Dijkstra et al. (2006)
  !             Laursen et al. (2009)
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid
  ! local variables
  real(kind=wp) :: xfreq_ref, uu1
  real(kind=wp) :: cost,sint,phi,cosp,sinp
  real(kind=wp) :: kx,ky,kz
  real(kind=wp) :: kx1,ky1,kz1,kr
  real(kind=wp) :: uz,ux,uy,uxy
  integer :: ix

  photon%nscatt_dust = photon%nscatt_dust + 1

  !--- Add absorbed portion to Jabs.
  if (.not. par%use_reduced_wgt) then
     if (rand_number() > par%albedo) then
       ! absorbed by dust.
       !2018-02-04, Jabs is now calculated in lab frame.
       if (par%save_Jabs) then
          uu1 = grid%vfx(photon%icell,photon%jcell,photon%kcell)*photon%kx + &
                grid%vfy(photon%icell,photon%jcell,photon%kcell)*photon%ky + &
                grid%vfz(photon%icell,photon%jcell,photon%kcell)*photon%kz
          xfreq_ref = (photon%xfreq + uu1)*(grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)
          ix        = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1
          if (ix >= 1 .and. ix <= grid%nxfreq) then
             !$OMP ATOMIC UPDATE
             grid%Jabs(ix) = grid%Jabs(ix) + photon%wgt
          endif
       endif
       photon%inside = .false.
       if (par%save_all_photons) then
          !--- xfreq_ref is required to record the photon frequency when it is absorbed, measured in the lab frame.
          !--- (2020.09.27)
          if (.not. par%save_Jabs) then
             uu1 = grid%vfx(photon%icell,photon%jcell,photon%kcell)*photon%kx + &
                   grid%vfy(photon%icell,photon%jcell,photon%kcell)*photon%ky + &
                   grid%vfz(photon%icell,photon%jcell,photon%kcell)*photon%kz
             xfreq_ref = (photon%xfreq + uu1)*(grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)
          endif
          photon%xfreq_ref = xfreq_ref
          photon%wgt       = 0.0_wp
       endif
       return
     endif
  else
     !2018-02-04, Jabs is now calculated in lab frame.
     if (par%save_Jabs) then
        uu1 = grid%vfx(photon%icell,photon%jcell,photon%kcell)*photon%kx + &
              grid%vfy(photon%icell,photon%jcell,photon%kcell)*photon%ky + &
              grid%vfz(photon%icell,photon%jcell,photon%kcell)*photon%kz
        xfreq_ref = (photon%xfreq + uu1)*(grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)
        ix        = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1
        if (ix >= 1 .and. ix <= grid%nxfreq) then
           !$OMP ATOMIC UPDATE
           grid%Jabs(ix) = grid%Jabs(ix) + photon%wgt * (1.0_wp - par%albedo)
        endif
     endif
     photon%wgt = photon%wgt * par%albedo
  endif
  if (par%save_peeloff) call peeling_dust_nostokes(photon,grid)

  !--- Select a new cos(theta) from the Henyey-Greenstein Function.
  cost = rand_henyey_greenstein(par%hgg)
  sint = sqrt(1.0d0-cost**2)

  ! New phi
  phi  = twopi * rand_number()
  cosp = cos(phi)
  sinp = sin(phi)

  ! New direction
  if (abs(photon%kz) >= 0.99999999999_wp) then
     !--- bug-fixed, 2021.04.22 (sometimes, we eject photons along the z-direction).
     photon%kx = sint*cosp
     photon%ky = sint*sinp
     photon%kz = cost
  else
     kx1 = photon%kx
     ky1 = photon%ky
     kz1 = photon%kz
     kr  = sqrt(kx1*kx1 + ky1*ky1)
     !--- 2018-02-03
     photon%kx = cost*kx1 + sint*(kz1*kx1*cosp - ky1*sinp)/kr
     photon%ky = cost*ky1 + sint*(kz1*ky1*cosp + kx1*sinp)/kr
     photon%kz = cost*kz1 - sint*cosp*kr
  endif
  return
  end subroutine scatter_dust_nostokes
  !=====================================
  subroutine scatter_resonance_nostokes(photon,grid)
  implicit none
  ! Calculate new direction cosines and wavelength after scattering
  ! Written by Kwang-il Seon, 2010/10/02
  ! 2016/10/26, updated to use the acceleration scheme of Ahn et al. (2002) and Dijkstra et al. (2006)
  !             Laursen et al. (2009)
  ! 2017-07-16, Photon package is removed and added to Jabs, instead of reducing weight, when the photon is absorbed by dust.
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid
  ! local variables
  real(kind=wp) :: cost,sint,phi,cosp,sinp
  real(kind=wp) :: kx,ky,kz
  real(kind=wp) :: kx1,ky1,kz1,kr
  real(kind=wp) :: phi2
  real(kind=wp) :: DnuHK,xfreq_atom,qK,qH,pH,pK
  real(kind=wp) :: ux,uy,uxy,uz,E1,g_recoil
  integer :: i1,i2,i3
  integer :: ix
  real(kind=wp), parameter :: one_over_three  = 1.0_wp/3.0_wp, &
                              two_over_three  = 2.0_wp/3.0_wp
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt

  i1 = photon%icell
  i2 = photon%jcell
  i3 = photon%kcell
  photon%nscatt_HI = photon%nscatt_HI + 1
#ifdef CALCP
  call add_to_Pa(photon,grid,photon%icell,photon%jcell,photon%kcell)
#endif

#ifdef FINE_STRUCTURE
  !--- Select an atom by which the photon is scattered and calculate new frequency.
  DnuHK = grid%DnuHK_ref_half/grid%Dfreq(i1,i2,i3)
  pH    = voigt(photon%xfreq+DnuHK, grid%voigt_a(i1,i2,i3)) * one_over_three
  pK    = voigt(photon%xfreq-DnuHK, grid%voigt_a(i1,i2,i3)) * two_over_three
  pH    = pH/(pH + pK)
  if (rand_number() < pH) then
     !--- scattered by an atom at H state (2P1/2).
     uz = rand_resonance_vz(photon%xfreq + DnuHK, grid%voigt_a(i1,i2,i3))
  else
     !--- scattered by an atom at K state (2P3/2).
     uz = rand_resonance_vz(photon%xfreq - DnuHK, grid%voigt_a(i1,i2,i3))
  endif

  !--- xfreq_atom = frequency in the atom's rest frame.
  xfreq_atom = photon%xfreq - uz

  !--- Select a new scattering angle (cos(theta)).
  !--- E1 is a function of frequency in the atom's rest frame.
  !--- We are assuming nu_H = nu_K. The difference make no practical difference.
  !--- qH = nu - nu_H, qK = nu - nu_K, nu_K = nu_H + DnuHK_Hz
  !--- qH = xfreq_atom, qK = qH - DnuHK
  qH    = xfreq_atom + DnuHK
  qK    = xfreq_atom - DnuHK
  E1    = (2.0_wp*qK*qH + qH**2)/(qK**2 + 2.0_wp*qH**2)
  cost  = rand_resonance(E1)
  sint  = sqrt(1.0_wp-cost**2)
#else
   !--- Select an atom by which the photon is scattered.
  uz = rand_resonance_vz(photon%xfreq, grid%voigt_a(i1,i2,i3))

  !--- xfreq_atom = frequency in the atom's rest frame.
  xfreq_atom = photon%xfreq - uz

  !--- Select a new scattering angle (cos(theta)).
  cost  = rand_rayleigh()
  sint  = sqrt(1.0_wp-cost**2)
#endif

  ! New phi
  phi  = twopi * rand_number()
  cosp = cos(phi)
  sinp = sin(phi)

  ! 2016-10-26, K.-I. Seon
  ! Note that the acceleration scheme uses a cut-off exponential distribution to obtain velocity amplitude, and
  ! divide it into x- and y- components, whereas the original scheme uses a gaussian distribution to obtain velocity components.
  ! This causes the difference in 1/sqrt(2) factor.
  if (par%core_skip .and. abs(photon%xfreq) < grid%xcrit) then
     phi2 = twopi * rand_number()
     uxy  = sqrt(grid%xcrit2 - log(rand_number()))
     ux   = uxy * cos(phi2)
     uy   = uxy * sin(phi2)
     photon%xfreq = xfreq_atom + uz*cost + (ux*cosp + uy*sinp)*sint
  else
     phi2 = twopi * rand_number()
     uxy  = sqrt(-log(rand_number()))
     ux   = uxy * cos(phi2)
     uy   = uxy * sin(phi2)
     photon%xfreq = xfreq_atom + uz*cost + (ux*cosp + uy*sinp)*sint
  endif

  if (par%recoil) then
     g_recoil     = g_recoil0 /grid%Dfreq(i1,i2,i3)
     photon%xfreq = photon%xfreq - g_recoil * (1.0_wp - cost)
  endif
  if (par%save_peeloff) call peeling_resonance_nostokes(photon,grid,xfreq_atom,[ux,uy,uz])

  ! New direction vector
  if (abs(photon%kz) >= 0.99999999999_wp) then
     !--- bug-fixed, 2021.04.22 (sometimes, we eject photons along the z-direction).
     photon%kx = sint*cosp
     photon%ky = sint*sinp
     photon%kz = cost
  else
     kx1 = photon%kx
     ky1 = photon%ky
     kz1 = photon%kz
     kr  = sqrt(kx1*kx1 + ky1*ky1)
     !--- 2018-02-03
     photon%kx = cost*kx1 + sint*(kz1*kx1*cosp - ky1*sinp)/kr
     photon%ky = cost*ky1 + sint*(kz1*ky1*cosp + kx1*sinp)/kr
     photon%kz = cost*kz1 - sint*cosp*kr
  endif
  return
  end subroutine scatter_resonance_nostokes

#ifdef CALCP
  subroutine add_to_Pa(photon,grid,icell,jcell,kcell)
  use define
  type(photon_type), intent(in)    :: photon
  type(grid_type),   intent(inout) :: grid
  integer,           intent(in)    :: icell,jcell,kcell
  real(kind=wp)                    :: rhokap
  if (icell > 0 .and. icell <= grid%nx .and. jcell > 0 .and. jcell <= grid%ny .and. kcell > 0 .and. kcell <= grid%nz) then
     if (grid%rhokap(icell,jcell,kcell) > 0.0_wp) then
        !--- Convert rhokap into density unit * distance2cm. (number/cm^3) * distance2cm.
        rhokap = grid%rhokap(icell,jcell,kcell) * grid%Dfreq(icell,jcell,kcell)/par%cross0
        select case (grid%geometry_JPa)
        case (3)
           !$OMP ATOMIC UPDATE
           grid%Pa(icell,jcell,kcell)               = grid%Pa(icell,jcell,kcell) + photon%wgt / rhokap
        case (2)
           !$OMP ATOMIC UPDATE
           grid%P2(grid%ind_cyl(icell,jcell),kcell) = grid%P2(grid%ind_cyl(icell,jcell),kcell) + photon%wgt / rhokap
        case (1)
           !$OMP ATOMIC UPDATE
           grid%P1(grid%ind_sph(icell,jcell,kcell)) = grid%P1(grid%ind_sph(icell,jcell,kcell)) + photon%wgt / rhokap
        case (-1)
           !$OMP ATOMIC UPDATE
           grid%P1(kcell)                           = grid%P1(kcell) + photon%wgt / rhokap
        end select
     endif
  endif
  end subroutine add_to_Pa
#endif

end module scatter_mod
