module scatter_mod
  use define
  use random
  use voigt_mod
  use line_mod
  !--- now, peelingoff subrotines are defined in define.f90 (2023.01.16).
  !use peelingoff_mod
  use mathlib
  use grid_mod, only: car_xcrit_local
  public
contains
  !=====================================
  subroutine scattering(photon,grid)
  use octree_mod, only: amr_grid, amr_line_profile
  use clump_mod,  only: cl_rhokap, cl_rhokapD, voigt_clump
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
  integer(int64) :: icl_s

  !--- band 2 (H-alpha, ly_beta line_type = 8): every interaction is a dust
  !--- event (dust-only opacity along the ray, see raytrace_to_tau). The
  !--- scatter_dust_Ha pointer selects the Cartesian or AMR band-2 dust routine.
  if (photon%iband /= 1) then
     call scatter_dust_Ha(photon,grid)
     return
  endif

  !--- AMR mode: use amr_grid arrays indexed by leaf index
  if (par%use_amr_grid) then
     i1 = photon%icell_amr
     if (par%DGR > 0.0_wp .and. associated(amr_grid%rhokapD)) then
        p_dust = amr_grid%rhokapD(i1) / &
            (amr_grid%rhokap(i1)*amr_line_profile(i1, photon%xfreq) + amr_grid%rhokapD(i1))
        if (rand_number() <= p_dust) then
           call scatter_dust(photon,grid)
        else
           call scatter_resonance(photon,grid)
        endif
     else
        call scatter_resonance(photon,grid)
     endif
     return
  endif

  !--- Clump mode: dust-vs-resonance split from per-clump opacities.
  !    Gas line opacity is multiplet-aware via voigt_clump; cl_rhokapD is
  !    the co-located dust continuum opacity (allocated only when DGR>0).
  if (par%use_clump_medium) then
     if (par%DGR > 0.0_wp .and. associated(cl_rhokapD)) then
        icl_s  = int(photon%icell_clump, int64)
        p_dust = cl_rhokapD(icl_s) / &
            (cl_rhokap(icl_s)*voigt_clump(photon%xfreq, icl_s) + cl_rhokapD(icl_s))
        if (rand_number() <= p_dust) then
           call scatter_dust(photon,grid)
        else
           call scatter_resonance(photon,grid)
        endif
     else
        call scatter_resonance(photon,grid)
     endif
     return
  endif

  if (par%DGR > 0.0_wp) then
     i1 = photon%icell
     i2 = photon%jcell
     i3 = photon%kcell
     p_dust = grid%rhokapD(i1,i2,i3)/(grid%rhokap(i1,i2,i3)*calc_voigt(grid,photon%xfreq,i1,i2,i3) + grid%rhokapD(i1,i2,i3))
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

    photon%nscatt_dust = photon%nscatt_dust + photon%wgt

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
    use clump_mod, only: cl_Dfreq, cl_Dfreq_ref
    implicit none
    ! Written by Kwang-il Seon, 2017/09/04
    type(photon_type), intent(inout) :: photon
    type(grid_type),   intent(inout) :: grid
    ! local variables
    !real(kind=wp) :: cost,cost2,sint,phi,phi1,phi2,cosp,sinp,cos2p,sin2p
    real(kind=wp) :: cost,sint,phi,phi1,phi2,cosp,sinp,cos2p,sin2p
    real(kind=wp) :: px,py,pz
    real(kind=wp) :: Prand, Pcomp, QoverI, UoverI, S12overS11
    real(kind=wp) :: S11,S12,S22,S33,S44
    real(kind=wp) :: Q0,U0,I1,Q1,U1,V1
    real(kind=wp) :: DnuHK,xfreq_atom,qK,qH,pH,pK
    real(kind=wp) :: ux,uy,uxy,uz,E1,g_recoil
    real(kind=wp) :: vth_ratio
    real(kind=wp) :: xcrit_cell, xcrit_cell2
    real(kind=wp), parameter :: three_over_four = 3.0_wp/4.0_wp, &
                                three_over_two  = 3.0_wp/2.0_wp, &
                                one_over_three  = 1.0_wp/3.0_wp, &
                                two_over_three  = 2.0_wp/3.0_wp, &
                                one_over_sqrt2  = 1.0_wp/sqrt(2.0_wp)
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt

    photon%nscatt_gas = photon%nscatt_gas + photon%wgt
#ifdef CALCP
    call add_to_Pa(photon,grid,photon%icell,photon%jcell,photon%kcell)
#endif

    call do_resonance(photon,grid,uz,xfreq_atom,cost,sint,S11,S22,S12,S33,S44)
    S12overS11  = S12/S11

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
    !-- Per-clump thermal-width ratio (see scatter_resonance_nostokes for rationale).
    if (par%use_clump_medium .and. photon%icell_clump > 0) then
       vth_ratio = cl_Dfreq(int(photon%icell_clump, int64)) / cl_Dfreq_ref
    else
       vth_ratio = 1.0_wp
    endif
    !-- RASCAS-style per-cell xcrit (Smith+15 Eq.35), see car_xcrit_local.
    if (par%core_skip) then
       call car_xcrit_local(grid, photon%icell, photon%jcell, photon%kcell, &
                            photon%x, photon%y, photon%z, xcrit_cell, xcrit_cell2)
    else
       xcrit_cell  = 0.0_wp
       xcrit_cell2 = 0.0_wp
    endif
    if (par%core_skip .and. abs(photon%xfreq) < xcrit_cell) then
       phi2 = twopi * rand_number()
       uxy  = sqrt(xcrit_cell2 - log(rand_number()))
       ux   = uxy * cos(phi2)
       uy   = uxy * sin(phi2)
       !-- For ly_alpha_HD D scatter: convert atom perpendicular velocity
       !-- from D Doppler units to caller's H-Doppler convention.
       !-- Same factor as for uz: divide by ratio_Dfreq_HD = Dfreq_H/Dfreq_D.
       if (line%line_type == 7 .and. line%selected_species_HD == 2) then
          ux = ux / line%ratio_Dfreq_HD
          uy = uy / line%ratio_Dfreq_HD
       endif
       ux = ux * vth_ratio
       uy = uy * vth_ratio
       photon%xfreq = xfreq_atom + uz*cost + (ux*cosp + uy*sinp)*sint
    else
       ux   = rand_gauss()*one_over_sqrt2
       uy   = rand_gauss()*one_over_sqrt2
       !-- For ly_alpha_HD D scatter: convert atom perpendicular velocity
       !-- from D Doppler units to caller's H-Doppler convention.
       !-- Same factor as for uz: divide by ratio_Dfreq_HD = Dfreq_H/Dfreq_D.
       if (line%line_type == 7 .and. line%selected_species_HD == 2) then
          ux = ux / line%ratio_Dfreq_HD
          uy = uy / line%ratio_Dfreq_HD
       endif
       ux = ux * vth_ratio
       uy = uy * vth_ratio
       photon%xfreq = xfreq_atom + uz*cost + (ux*cosp + uy*sinp)*sint
    endif

    if (par%recoil) then
       if (par%use_clump_medium .and. photon%icell_clump > 0) then
          if (line%line_type == 7 .and. line%selected_species_HD == 2) then
             g_recoil = line%g_recoil0_D / cl_Dfreq(int(photon%icell_clump, int64))
          else
             g_recoil = line%g_recoil0   / cl_Dfreq(int(photon%icell_clump, int64))
          endif
       else
          if (line%line_type == 7 .and. line%selected_species_HD == 2) then
             g_recoil = line%g_recoil0_D / grid%Dfreq(photon%icell,photon%jcell,photon%kcell)
          else
             g_recoil = line%g_recoil0   / grid%Dfreq(photon%icell,photon%jcell,photon%kcell)
          endif
       endif
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

  photon%nscatt_dust = photon%nscatt_dust + photon%wgt

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
       !--- ly_beta bookkeeping: band-1 dust-absorbed weight.
       if (line%line_type == 8) par%W_abs1 = par%W_abs1 + photon%wgt
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
     !--- ly_beta bookkeeping: band-1 dust-absorbed weight (reduced-weight scheme).
     if (line%line_type == 8) par%W_abs1 = par%W_abs1 + photon%wgt * (1.0_wp - par%albedo)
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
  subroutine scatter_dust_Ha_nostokes(photon,grid)
  implicit none
  !--- Band-2 (H-alpha, ly_beta line_type = 8) dust scattering event (Cartesian).
  !--- Mirrors scatter_dust_nostokes with the band-2 dust properties
  !--- (par%albedo_Ha, par%hgg_Ha) and band-2 tallies (Jabs_Ha, W_abs2,
  !--- peel_Ha). photon%xfreq is the LAB-frame frequency in REFERENCE Doppler
  !--- units (dust scattering is elastic) -> no fluid shift, no Dfreq scaling.
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid
  ! local variables
  real(kind=wp) :: cost,sint,phi,cosp,sinp
  real(kind=wp) :: kx1,ky1,kz1,kr
  integer :: ix

  photon%nscatt_dust = photon%nscatt_dust + photon%wgt

  !--- Add absorbed portion to Jabs_Ha.
  if (.not. par%use_reduced_wgt) then
     if (rand_number() > par%albedo_Ha) then
       ! absorbed by dust.
       if (par%save_Jabs .and. associated(grid%Jabs_Ha)) then
          ix = floor((photon%xfreq - grid%xfreq_min_Ha)/grid%dxfreq_Ha)+1
          if (ix >= 1 .and. ix <= grid%nxfreq_Ha) then
             grid%Jabs_Ha(ix) = grid%Jabs_Ha(ix) + photon%wgt
          endif
       endif
       photon%inside = .false.
       par%W_abs2    = par%W_abs2 + photon%wgt
       if (par%save_all_photons) then
          photon%xfreq_ref = photon%xfreq
          photon%wgt       = 0.0_wp
       endif
       return
     endif
  else
     if (par%save_Jabs .and. associated(grid%Jabs_Ha)) then
        ix = floor((photon%xfreq - grid%xfreq_min_Ha)/grid%dxfreq_Ha)+1
        if (ix >= 1 .and. ix <= grid%nxfreq_Ha) then
           grid%Jabs_Ha(ix) = grid%Jabs_Ha(ix) + photon%wgt * (1.0_wp - par%albedo_Ha)
        endif
     endif
     par%W_abs2 = par%W_abs2 + photon%wgt * (1.0_wp - par%albedo_Ha)
     photon%wgt = photon%wgt * par%albedo_Ha
  endif
  !--- Peel-off with the H-alpha HG phase function into peel_Ha
  !--- (inside-observer path is rejected for ly_beta in setup.f90).
  if (par%save_peeloff) call peeling_dust_Ha(photon,grid)

  !--- Select a new cos(theta) from the Henyey-Greenstein Function.
  cost = rand_henyey_greenstein(par%hgg_Ha)
  sint = sqrt(1.0d0-cost**2)

  ! New phi
  phi  = twopi * rand_number()
  cosp = cos(phi)
  sinp = sin(phi)

  ! New direction
  if (abs(photon%kz) >= 0.99999999999_wp) then
     photon%kx = sint*cosp
     photon%ky = sint*sinp
     photon%kz = cost
  else
     kx1 = photon%kx
     ky1 = photon%ky
     kz1 = photon%kz
     kr  = sqrt(kx1*kx1 + ky1*ky1)
     photon%kx = cost*kx1 + sint*(kz1*kx1*cosp - ky1*sinp)/kr
     photon%ky = cost*ky1 + sint*(kz1*ky1*cosp + kx1*sinp)/kr
     photon%kz = cost*kz1 - sint*cosp*kr
  endif
  return
  end subroutine scatter_dust_Ha_nostokes
  !=====================================
  subroutine scatter_resonance_nostokes(photon,grid)
  use clump_mod, only: cl_Dfreq, cl_Dfreq_ref
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
  real(kind=wp) :: vth_ratio
  real(kind=wp) :: xcrit_cell, xcrit_cell2
  real(kind=wp) :: u1_fluid
  logical :: converted
  integer :: i1,i2,i3
  integer :: ix
  real(kind=wp), parameter :: one_over_three  = 1.0_wp/3.0_wp, &
                              two_over_three  = 2.0_wp/3.0_wp
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt

  i1 = photon%icell
  i2 = photon%jcell
  i3 = photon%kcell
  photon%nscatt_gas = photon%nscatt_gas + photon%wgt
#ifdef CALCP
  call add_to_Pa(photon,grid,photon%icell,photon%jcell,photon%kcell)
#endif

  call do_resonance(photon,grid, uz,xfreq_atom,cost,sint)

  !--- ly_beta (line_type = 8): do_resonance8 sets photon%iband = 2 when the
  !--- 3p->2s conversion branch fires (transmutation to the H-alpha band).
  converted = (line%line_type == 8 .and. photon%iband == 2)
  if (converted) then
     par%W_conv = par%W_conv + photon%wgt
#ifdef CALCP
     call add_to_Pconv(photon,grid,i1,i2,i3)
#endif
  endif

  ! New phi
  phi  = twopi * rand_number()
  cosp = cos(phi)
  sinp = sin(phi)

  !-- Per-clump thermal-width ratio: scales perpendicular atom velocities and
  !-- recoil shift into REF Doppler units when clump T differs from ref T.
  !-- = 1 for uniform-T (cl_Dfreq(icl) == cl_Dfreq_ref), preserving prior behavior.
  if (par%use_clump_medium .and. photon%icell_clump > 0) then
     vth_ratio = cl_Dfreq(int(photon%icell_clump, int64)) / cl_Dfreq_ref
  else
     vth_ratio = 1.0_wp
  endif

  ! 2016-10-26, K.-I. Seon
  ! Note that the acceleration scheme uses a cut-off exponential distribution to obtain velocity amplitude, and
  ! divide it into x- and y- components, whereas the original scheme uses a gaussian distribution to obtain velocity components.
  ! This causes the difference in 1/sqrt(2) factor.
  !-- RASCAS-style per-cell xcrit (Smith+15 Eq.35), see car_xcrit_local.
  if (par%core_skip) then
     call car_xcrit_local(grid, i1, i2, i3, photon%x, photon%y, photon%z, &
                          xcrit_cell, xcrit_cell2)
  else
     xcrit_cell  = 0.0_wp
     xcrit_cell2 = 0.0_wp
  endif
  if (par%core_skip .and. abs(photon%xfreq) < xcrit_cell) then
     phi2 = twopi * rand_number()
     uxy  = sqrt(xcrit_cell2 - log(rand_number()))
     ux   = uxy * cos(phi2)
     uy   = uxy * sin(phi2)
     !-- For ly_alpha_HD D scatter: convert atom perpendicular velocity
     !-- from D Doppler units to caller's H-Doppler convention.
     !-- Same factor as for uz: divide by ratio_Dfreq_HD = Dfreq_H/Dfreq_D.
     if (line%line_type == 7 .and. line%selected_species_HD == 2) then
        ux = ux / line%ratio_Dfreq_HD
        uy = uy / line%ratio_Dfreq_HD
     endif
     ux = ux * vth_ratio
     uy = uy * vth_ratio
     photon%xfreq = xfreq_atom + uz*cost + (ux*cosp + uy*sinp)*sint
  else
     phi2 = twopi * rand_number()
     uxy  = sqrt(-log(rand_number()))
     ux   = uxy * cos(phi2)
     uy   = uxy * sin(phi2)
     !-- For ly_alpha_HD D scatter: convert atom perpendicular velocity
     !-- from D Doppler units to caller's H-Doppler convention.
     !-- Same factor as for uz: divide by ratio_Dfreq_HD = Dfreq_H/Dfreq_D.
     if (line%line_type == 7 .and. line%selected_species_HD == 2) then
        ux = ux / line%ratio_Dfreq_HD
        uy = uy / line%ratio_Dfreq_HD
     endif
     ux = ux * vth_ratio
     uy = uy * vth_ratio
     photon%xfreq = xfreq_atom + uz*cost + (ux*cosp + uy*sinp)*sint
  endif

  !--- Recoil is skipped for the conversion channel: the correction below is
  !--- defined for band-1 (Ly-beta) re-emission; par%recoil defaults to .false.
  !--- anyway. (Existing lines: converted is always .false. -> unchanged.)
  if (par%recoil .and. .not. converted) then
     if (par%use_clump_medium .and. photon%icell_clump > 0) then
        if (line%line_type == 7 .and. line%selected_species_HD == 2) then
           g_recoil = line%g_recoil0_D / cl_Dfreq(int(photon%icell_clump, int64))
        else
           g_recoil = line%g_recoil0   / cl_Dfreq(int(photon%icell_clump, int64))
        endif
     else
        if (line%line_type == 7 .and. line%selected_species_HD == 2) then
           g_recoil = line%g_recoil0_D / grid%Dfreq(i1,i2,i3)
        else
           g_recoil = line%g_recoil0   / grid%Dfreq(i1,i2,i3)
        endif
     endif
     photon%xfreq = photon%xfreq - g_recoil * (1.0_wp - cost)
  endif
  if (par%save_peeloff) then
     if (converted) then
        !--- direct fluorescent (H-alpha) peel of the newborn band-2 photon
        !--- (inside-observer path is rejected for ly_beta in setup.f90).
        call peeling_conversion_Ha(photon,grid,[ux,uy,uz])
     else
        call peeling_resonance_nostokes(photon,grid,xfreq_atom,[ux,uy,uz])
     endif
  endif

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

  !--- Band transmutation (ly_beta conversion): the H-alpha photon is emitted
  !--- at line center in the ATOM frame (plan section 8.4), so the emitted
  !--- frequency is the atom-velocity projection onto the NEW direction:
  !---   x_com = photon%xfreq - xfreq_atom = u_atom . k_new  (local Doppler units)
  !--- Transform to the GLOBAL (lab) frame in REFERENCE Doppler units, exactly
  !--- mirroring the band-1 escape transformation (add local fluid velocity
  !--- along k_new, then rescale by Dfreq_local/Dfreq_ref). v_th is identical
  !--- for both bands, so the numerical Doppler-unit scale carries over.
  !--- Band 2 NEVER updates xfreq again (dust scattering is elastic).
  if (converted) then
     u1_fluid = grid%vfx(i1,i2,i3)*photon%kx + grid%vfy(i1,i2,i3)*photon%ky + &
                grid%vfz(i1,i2,i3)*photon%kz
     photon%xfreq = (photon%xfreq - xfreq_atom + u1_fluid) * &
                    (grid%Dfreq(i1,i2,i3) / grid%Dfreq_ref)
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
        rhokap = grid%rhokap(icell,jcell,kcell) * grid%Dfreq(icell,jcell,kcell)/line%cross0
        select case (grid%geometry_JPa)
        case (3)
           !$OMP ATOMIC UPDATE
           grid%Pa(icell,jcell,kcell)               = grid%Pa(icell,jcell,kcell) + photon%wgt / rhokap
        case (2)
           !--- guard: corner cells can have cylindrical radius > rmax (ind_cyl > nr)
           if (grid%ind_cyl(icell,jcell) >= 1 .and. grid%ind_cyl(icell,jcell) <= grid%nr) then
           !$OMP ATOMIC UPDATE
           grid%P2(grid%ind_cyl(icell,jcell),kcell) = grid%P2(grid%ind_cyl(icell,jcell),kcell) + photon%wgt / rhokap
           endif
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

  subroutine add_to_Pconv(photon,grid,icell,jcell,kcell)
  !--- Ly-beta (line_type = 8) conversion-rate map: mirrors add_to_Pa (same
  !--- per-atom-rate convention, same binning, same guards), but accumulates
  !--- into the Pconv arrays at conversion events only.
  !--- Expectation: Pconv/Pa -> P_down(2) = 0.11834.
  use define
  type(photon_type), intent(in)    :: photon
  type(grid_type),   intent(inout) :: grid
  integer,           intent(in)    :: icell,jcell,kcell
  real(kind=wp)                    :: rhokap
  if (icell > 0 .and. icell <= grid%nx .and. jcell > 0 .and. jcell <= grid%ny .and. kcell > 0 .and. kcell <= grid%nz) then
     if (grid%rhokap(icell,jcell,kcell) > 0.0_wp) then
        !--- Convert rhokap into density unit * distance2cm. (number/cm^3) * distance2cm.
        rhokap = grid%rhokap(icell,jcell,kcell) * grid%Dfreq(icell,jcell,kcell)/line%cross0
        select case (grid%geometry_JPa)
        case (3)
           !$OMP ATOMIC UPDATE
           grid%Pc(icell,jcell,kcell)                = grid%Pc(icell,jcell,kcell) + photon%wgt / rhokap
        case (2)
           !--- guard: corner cells can have cylindrical radius > rmax (ind_cyl > nr)
           if (grid%ind_cyl(icell,jcell) >= 1 .and. grid%ind_cyl(icell,jcell) <= grid%nr) then
           !$OMP ATOMIC UPDATE
           grid%Pc2(grid%ind_cyl(icell,jcell),kcell) = grid%Pc2(grid%ind_cyl(icell,jcell),kcell) + photon%wgt / rhokap
           endif
        case (1)
           !$OMP ATOMIC UPDATE
           grid%Pc1(grid%ind_sph(icell,jcell,kcell)) = grid%Pc1(grid%ind_sph(icell,jcell,kcell)) + photon%wgt / rhokap
        case (-1)
           !$OMP ATOMIC UPDATE
           grid%Pc1(kcell)                           = grid%Pc1(kcell) + photon%wgt / rhokap
        end select
     endif
  endif
  end subroutine add_to_Pconv
#endif

  !=====================================
  subroutine scatter_resonance_clump_nostokes(photon, grid)
  !---------------------------------------------------------------------------
  ! Wrapper for the overlap-aware clump path (nostokes).
  ! photon%xfreq arrives in the GLOBAL frame.  Apply Doppler shift to the
  ! owner clump's frame, scatter, then shift back to global.
  !---------------------------------------------------------------------------
  use clump_mod, only: ulos_clump
  implicit none
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid
  integer(int64) :: icl
  real(kind=wp)  :: u_in, u_out
  icl  = int(photon%icell_clump, int64)
  u_in = ulos_clump(icl, photon%kx, photon%ky, photon%kz)
  photon%xfreq = photon%xfreq - u_in
  call scatter_resonance_nostokes(photon, grid)
  u_out = ulos_clump(icl, photon%kx, photon%ky, photon%kz)
  photon%xfreq = photon%xfreq + u_out
  end subroutine scatter_resonance_clump_nostokes
  !=====================================

  !=====================================
  subroutine scatter_resonance_clump_stokes(photon, grid)
  !---------------------------------------------------------------------------
  ! Wrapper for the overlap-aware clump path (stokes).
  ! Same Doppler logic as the nostokes variant.
  !---------------------------------------------------------------------------
  use clump_mod, only: ulos_clump
  implicit none
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(inout) :: grid
  integer(int64) :: icl
  real(kind=wp)  :: u_in, u_out
  icl  = int(photon%icell_clump, int64)
  u_in = ulos_clump(icl, photon%kx, photon%ky, photon%kz)
  photon%xfreq = photon%xfreq - u_in
  call scatter_resonance_stokes(photon, grid)
  u_out = ulos_clump(icl, photon%kx, photon%ky, photon%kz)
  photon%xfreq = photon%xfreq + u_out
  end subroutine scatter_resonance_clump_stokes
  !=====================================

end module scatter_mod
