module scattering_amr_mod
  !-----------------------------------------------------------------------
  ! AMR scattering routines.  These mirror the routines in
  ! scattering_car_v2a.f90 but read physical quantities from amr_grid
  ! instead of the Cartesian grid arrays.
  !
  ! photon%icell_amr = current leaf index in amr_grid.
  !
  ! In AMR mode, the procedure pointers are reassigned to these routines:
  !   scatter_dust      => scatter_dust_nostokes_amr / scatter_dust_stokes_amr
  !   scatter_resonance => scatter_resonance_nostokes_amr / scatter_resonance_stokes_amr
  !   do_resonance      => do_resonance1_amr  (or do_resonance2_amr etc.)
  !-----------------------------------------------------------------------
  use octree_mod
  use voigt_mod
  use random
  use mathlib
  implicit none
  private

  public :: scatter_dust_nostokes_amr
  public :: scatter_dust_stokes_amr
  public :: scatter_resonance_nostokes_amr
  public :: scatter_resonance_stokes_amr
  public :: do_resonance1_amr
  public :: do_resonance2_amr
  public :: do_resonance4_amr
  public :: do_resonance5_amr
  public :: do_resonance6_amr
  public :: scattering_amr_dispatch

contains

  !=========================================================================
  ! Top-level scattering dispatcher for AMR mode.
  ! Replaces the Cartesian scattering() routine when par%use_amr_grid=.true.
  !=========================================================================
  subroutine scattering_amr_dispatch(photon, grid)
    use define
    implicit none
    type(photon_type), intent(inout) :: photon
    type(grid_type),   intent(inout) :: grid
    integer  :: il
    real(wp) :: p_dust

    il = photon%icell_amr
    if (il <= 0) return

    if (par%DGR > 0.0_wp .and. associated(amr_grid%rhokapD)) then
      p_dust = amr_grid%rhokapD(il) / &
          (amr_grid%rhokap(il) * voigt(photon%xfreq, amr_grid%voigt_a(il)) + amr_grid%rhokapD(il))
      if (rand_number() <= p_dust) then
        call scatter_dust(photon, grid)
      else
        call scatter_resonance(photon, grid)
      end if
    else
      call scatter_resonance(photon, grid)
    end if
  end subroutine scattering_amr_dispatch

  !=========================================================================
  ! AMR dust scattering (no Stokes).
  !=========================================================================
  subroutine scatter_dust_nostokes_amr(photon, grid)
    use define
    implicit none
    type(photon_type), intent(inout) :: photon
    type(grid_type),   intent(inout) :: grid

    integer  :: il, ix
    real(wp) :: uu1, xfreq_ref
    real(wp) :: cost, sint, phi, cosp, sinp
    real(wp) :: kx1, ky1, kz1, kr

    il = photon%icell_amr
    photon%nscatt_dust = photon%nscatt_dust + photon%wgt

    ! Handle dust absorption (reduced-weight or stochastic)
    if (.not. par%use_reduced_wgt) then
      if (rand_number() > par%albedo) then
        ! Photon absorbed by dust
        if (par%save_Jabs .and. allocated(amr_grid%Jabs)) then
          uu1       = amr_grid%vfx(il)*photon%kx + amr_grid%vfy(il)*photon%ky + amr_grid%vfz(il)*photon%kz
          xfreq_ref = (photon%xfreq + uu1) * (amr_grid%Dfreq(il) / amr_grid%Dfreq_ref)
          ix        = floor((xfreq_ref - amr_grid%xfreq_min) / amr_grid%dxfreq) + 1
          if (ix >= 1 .and. ix <= amr_grid%nxfreq) then
            !$OMP ATOMIC UPDATE
            amr_grid%Jabs(ix) = amr_grid%Jabs(ix) + photon%wgt
          end if
        end if
        photon%inside = .false.
        return
      end if
    else
      if (par%save_Jabs .and. allocated(amr_grid%Jabs)) then
        uu1       = amr_grid%vfx(il)*photon%kx + amr_grid%vfy(il)*photon%ky + amr_grid%vfz(il)*photon%kz
        xfreq_ref = (photon%xfreq + uu1) * (amr_grid%Dfreq(il) / amr_grid%Dfreq_ref)
        ix        = floor((xfreq_ref - amr_grid%xfreq_min) / amr_grid%dxfreq) + 1
        if (ix >= 1 .and. ix <= amr_grid%nxfreq) then
          !$OMP ATOMIC UPDATE
          amr_grid%Jabs(ix) = amr_grid%Jabs(ix) + photon%wgt * (1.0_wp - par%albedo)
        end if
      end if
      photon%wgt = photon%wgt * par%albedo
    end if

    if (par%save_peeloff) call peeling_dust_nostokes(photon, grid)

    ! New scattering angle (Henyey–Greenstein or isotropic)
    if (par%hgg == 0.0_wp) then
      cost = 2.0_wp * rand_number() - 1.0_wp
    else
      cost = (1.0_wp + par%hgg**2 - &
              ((1.0_wp - par%hgg**2) / (1.0_wp - par%hgg + 2.0_wp*par%hgg*rand_number()))**2) &
              / (2.0_wp * par%hgg)
    end if
    sint = sqrt(max(0.0_wp, 1.0_wp - cost**2))
    phi  = twopi * rand_number()
    cosp = cos(phi);  sinp = sin(phi)

    ! Rotate direction
    kx1 = photon%kx;  ky1 = photon%ky;  kz1 = photon%kz
    if (abs(kz1) >= 0.99999999999_wp) then
      photon%kx = sint * cosp
      photon%ky = sint * sinp
      photon%kz = cost
    else
      kr         = sqrt(kx1**2 + ky1**2)
      photon%kx  = cost*kx1 + sint*(kz1*kx1*cosp - ky1*sinp)/kr
      photon%ky  = cost*ky1 + sint*(kz1*ky1*cosp + kx1*sinp)/kr
      photon%kz  = cost*kz1 - sint*cosp*kr
    end if
  end subroutine scatter_dust_nostokes_amr

  !=========================================================================
  ! AMR resonance scattering (no Stokes).
  !=========================================================================
  subroutine scatter_resonance_nostokes_amr(photon, grid)
    use define
    implicit none
    type(photon_type), intent(inout) :: photon
    type(grid_type),   intent(inout) :: grid

    integer  :: il
    real(wp) :: uz, xfreq_atom, cost, sint, phi, phi2, cosp, sinp
    real(wp) :: ux, uy, uxy
    real(wp) :: kx1, ky1, kz1, kr
    real(wp) :: g_recoil

    il = photon%icell_amr
    photon%nscatt_gas = photon%nscatt_gas + photon%wgt

    ! Sample scattering atom and scattering angle via procedure pointer do_resonance
    call do_resonance(photon, grid, uz, xfreq_atom, cost, sint)

    phi  = twopi * rand_number()
    cosp = cos(phi);  sinp = sin(phi)

    ! Atom thermal velocity component perpendicular to k (core-skip or full)
    if (par%core_skip .and. abs(photon%xfreq) < amr_grid%xcrit) then
      phi2 = twopi * rand_number()
      uxy  = sqrt(amr_grid%xcrit2 - log(rand_number()))
      ux   = uxy * cos(phi2);  uy = uxy * sin(phi2)
    else
      phi2 = twopi * rand_number()
      uxy  = sqrt(-log(rand_number()))
      ux   = uxy * cos(phi2);  uy = uxy * sin(phi2)
    end if
    photon%xfreq = xfreq_atom + uz*cost + (ux*cosp + uy*sinp)*sint

    if (par%recoil) then
      g_recoil     = line%g_recoil0 / amr_grid%Dfreq(il)
      photon%xfreq = photon%xfreq - g_recoil * (1.0_wp - cost)
    end if

    if (par%save_peeloff) call peeling_resonance_nostokes(photon, grid, xfreq_atom, [ux, uy, uz])

    ! New direction vector
    kx1 = photon%kx;  ky1 = photon%ky;  kz1 = photon%kz
    if (abs(kz1) >= 0.99999999999_wp) then
      photon%kx = sint * cosp
      photon%ky = sint * sinp
      photon%kz = cost
    else
      kr        = sqrt(kx1**2 + ky1**2)
      photon%kx = cost*kx1 + sint*(kz1*kx1*cosp - ky1*sinp)/kr
      photon%ky = cost*ky1 + sint*(kz1*ky1*cosp + kx1*sinp)/kr
      photon%kz = cost*kz1 - sint*cosp*kr
    end if
  end subroutine scatter_resonance_nostokes_amr

  !=========================================================================
  ! AMR dust scattering with Stokes polarization.
  ! Mirrors scatter_dust_stokes in scattering_car.f90.
  !=========================================================================
  subroutine scatter_dust_stokes_amr(photon, grid)
    use define
    implicit none
    type(photon_type), intent(inout) :: photon
    type(grid_type),   intent(inout) :: grid

    integer  :: il, ix
    real(wp) :: uu1, xfreq_ref
    real(wp) :: cost, sint, phi, phi1, cosp, sinp, cos2p, sin2p
    real(wp) :: px, py, pz
    real(wp) :: Prand, Pcomp, S12overS11
    real(wp) :: S11, S12, S33, S34
    real(wp) :: Q0, U0, I1, Q1, U1, V1

    il = photon%icell_amr
    photon%nscatt_dust = photon%nscatt_dust + photon%wgt

    if (.not. par%use_reduced_wgt) then
      if (rand_number() > par%albedo) then
        if (par%save_Jabs .and. allocated(amr_grid%Jabs)) then
          uu1       = amr_grid%vfx(il)*photon%kx + amr_grid%vfy(il)*photon%ky + amr_grid%vfz(il)*photon%kz
          xfreq_ref = (photon%xfreq + uu1) * (amr_grid%Dfreq(il) / amr_grid%Dfreq_ref)
          ix        = floor((xfreq_ref - amr_grid%xfreq_min) / amr_grid%dxfreq) + 1
          if (ix >= 1 .and. ix <= amr_grid%nxfreq) then
            !$OMP ATOMIC UPDATE
            amr_grid%Jabs(ix) = amr_grid%Jabs(ix) + photon%wgt
          end if
        end if
        photon%inside = .false.
        return
      end if
    else
      if (par%save_Jabs .and. allocated(amr_grid%Jabs)) then
        uu1       = amr_grid%vfx(il)*photon%kx + amr_grid%vfy(il)*photon%ky + amr_grid%vfz(il)*photon%kz
        xfreq_ref = (photon%xfreq + uu1) * (amr_grid%Dfreq(il) / amr_grid%Dfreq_ref)
        ix        = floor((xfreq_ref - amr_grid%xfreq_min) / amr_grid%dxfreq) + 1
        if (ix >= 1 .and. ix <= amr_grid%nxfreq) then
          !$OMP ATOMIC UPDATE
          amr_grid%Jabs(ix) = amr_grid%Jabs(ix) + photon%wgt * (1.0_wp - par%albedo)
        end if
      end if
      photon%wgt = photon%wgt * par%albedo
    end if

    if (par%save_peeloff) call peeling_dust_stokes(photon, grid)

    ! Sample new polar angle from scattering phase function via alias method
    cost = rand_alias_linear(scatt_mat%phase_PDF, scatt_mat%alias, scatt_mat%coss, scatt_mat%S11)
    sint = sqrt(1.0_wp - cost**2)

    call interp_eq(scatt_mat%coss, scatt_mat%S11, cost, S11)
    call interp_eq(scatt_mat%coss, scatt_mat%S12, cost, S12)
    call interp_eq(scatt_mat%coss, scatt_mat%S33, cost, S33)
    call interp_eq(scatt_mat%coss, scatt_mat%S34, cost, S34)
    S12overS11 = S12 / S11

    ! Sample azimuthal angle via rejection method (Stokes-weighted)
    do while (.true.)
      phi   = twopi * rand_number()
      phi1  = 2.0_wp * phi
      Prand = (1.0_wp + abs(S12overS11) * sqrt(photon%Q**2 + photon%U**2)) * rand_number()
      Pcomp =  1.0_wp + S12overS11 * (photon%Q*cos(phi1) + photon%U*sin(phi1))
      if (Prand <= Pcomp) exit
    end do
    cosp  = cos(phi);  sinp  = sin(phi)
    cos2p = 2.0_wp*cosp*cosp - 1.0_wp
    sin2p = 2.0_wp*cosp*sinp

    ! Update Stokes parameters
    Q0 =  cos2p*photon%Q + sin2p*photon%U
    U0 = -sin2p*photon%Q + cos2p*photon%U
    I1 =  S11 + S12*Q0
    Q1 =  S12 + S11*Q0
    U1 =  S33*U0 + S34*photon%V
    V1 = -S34*U0 + S33*photon%V
    photon%I = 1.0_wp
    photon%Q = Q1 / I1
    photon%U = U1 / I1
    photon%V = V1 / I1

    ! Update reference frame triad and propagation direction
    px = cosp*photon%mx + sinp*photon%nx
    py = cosp*photon%my + sinp*photon%ny
    pz = cosp*photon%mz + sinp*photon%nz

    photon%nx = cosp*photon%nx - sinp*photon%mx
    photon%ny = cosp*photon%ny - sinp*photon%my
    photon%nz = cosp*photon%nz - sinp*photon%mz

    photon%mx = cost*px - sint*photon%kx
    photon%my = cost*py - sint*photon%ky
    photon%mz = cost*pz - sint*photon%kz

    photon%kx = sint*px + cost*photon%kx
    photon%ky = sint*py + cost*photon%ky
    photon%kz = sint*pz + cost*photon%kz
  end subroutine scatter_dust_stokes_amr

  !=========================================================================
  ! AMR resonance scattering with Stokes polarization.
  ! Mirrors scatter_resonance_stokes in scattering_car.f90.
  !=========================================================================
  subroutine scatter_resonance_stokes_amr(photon, grid)
    use define
    implicit none
    type(photon_type), intent(inout) :: photon
    type(grid_type),   intent(inout) :: grid

    integer  :: il
    real(wp) :: uz, xfreq_atom, cost, sint, phi, phi1, phi2, cosp, sinp, cos2p, sin2p
    real(wp) :: px, py, pz
    real(wp) :: Prand, Pcomp, S12overS11
    real(wp) :: S11, S12, S22, S33, S44
    real(wp) :: Q0, U0, I1, Q1, U1, V1
    real(wp) :: ux, uy, uxy, g_recoil

    il = photon%icell_amr
    photon%nscatt_gas = photon%nscatt_gas + photon%wgt

    call do_resonance(photon, grid, uz, xfreq_atom, cost, sint, S11, S22, S12, S33, S44)
    S12overS11 = S12 / S11

    ! Sample azimuthal angle via rejection method (Stokes-weighted)
    do while (.true.)
      phi   = twopi * rand_number()
      phi1  = 2.0_wp * phi
      Prand = (1.0_wp + abs(S12overS11) * sqrt(photon%Q**2 + photon%U**2)) * rand_number()
      Pcomp =  1.0_wp + S12overS11 * (photon%Q*cos(phi1) + photon%U*sin(phi1))
      if (Prand <= Pcomp) exit
    end do
    cosp = cos(phi);  sinp = sin(phi)

    ! Update photon frequency (core-skip or full thermal sampling)
    if (par%core_skip .and. abs(photon%xfreq) < amr_grid%xcrit) then
      phi2 = twopi * rand_number()
      uxy  = sqrt(amr_grid%xcrit2 - log(rand_number()))
      ux   = uxy * cos(phi2);  uy = uxy * sin(phi2)
      photon%xfreq = xfreq_atom + uz*cost + (ux*cosp + uy*sinp)*sint
    else
      ux   = rand_gauss() * one_over_sqrt2
      uy   = rand_gauss() * one_over_sqrt2
      photon%xfreq = xfreq_atom + uz*cost + (ux*cosp + uy*sinp)*sint
    end if

    if (par%recoil) then
      g_recoil     = line%g_recoil0 / amr_grid%Dfreq(il)
      photon%xfreq = photon%xfreq - g_recoil * (1.0_wp - cost)
    end if

    ! Peeling-off before updating triad
    if (par%save_peeloff) call peeling_resonance_stokes(photon, grid, xfreq_atom, [ux, uy, uz])

    cos2p = 2.0_wp*cosp*cosp - 1.0_wp
    sin2p = 2.0_wp*cosp*sinp

    ! Update Stokes parameters
    Q0 =  cos2p*photon%Q + sin2p*photon%U
    U0 = -sin2p*photon%Q + cos2p*photon%U
    I1 =  S11 + S12*Q0
    Q1 =  S12 + S22*Q0
    U1 =  S33*U0
    V1 =  S44*photon%V
    photon%I = 1.0_wp
    photon%Q = Q1 / I1
    photon%U = U1 / I1
    photon%V = V1 / I1

    ! Update reference frame triad and propagation direction
    px = cosp*photon%mx + sinp*photon%nx
    py = cosp*photon%my + sinp*photon%ny
    pz = cosp*photon%mz + sinp*photon%nz

    photon%nx = cosp*photon%nx - sinp*photon%mx
    photon%ny = cosp*photon%ny - sinp*photon%my
    photon%nz = cosp*photon%nz - sinp*photon%mz

    photon%mx = cost*px - sint*photon%kx
    photon%my = cost*py - sint*photon%ky
    photon%mz = cost*pz - sint*photon%kz

    photon%kx = sint*px + cost*photon%kx
    photon%ky = sint*py + cost*photon%ky
    photon%kz = sint*pz + cost*photon%kz
  end subroutine scatter_resonance_stokes_amr

  !=========================================================================
  ! AMR do_resonance1: single-line resonance (Ly-alpha H I, FeII, etc.)
  ! Uses amr_grid%voigt_a(photon%icell_amr) instead of grid%voigt_a(i,j,k).
  !=========================================================================
  subroutine do_resonance1_amr(photon, grid, uz, xfreq_atom, cost, sint, &
                                S11, S22, S12, S33, S44)
    use define
    implicit none
    type(photon_type), intent(inout) :: photon
    type(grid_type),   intent(in)    :: grid
    real(wp),          intent(out)   :: uz, xfreq_atom, cost, sint
    real(wp), optional, intent(out)  :: S11, S22, S12, S33, S44

    integer  :: il
    real(wp) :: cost2

    il = photon%icell_amr
    uz         = rand_resonance_vz(photon%xfreq, amr_grid%voigt_a(il))
    xfreq_atom = photon%xfreq - uz

    photon%E1 = line%E1
    photon%E2 = line%E2
    photon%E3 = line%E3

    cost  = rand_resonance(photon%E1)
    cost2 = cost**2
    sint  = sqrt(1.0_wp - cost2)

    if (present(S44)) then
      S22 = three_over_four * photon%E1 * (cost2 + 1.0_wp)
      S11 = S22 + photon%E2
      S12 = three_over_four * photon%E1 * (cost2 - 1.0_wp)
      S33 = three_over_two  * photon%E1 * cost
      S44 = three_over_two  * photon%E3 * cost
    end if
  end subroutine do_resonance1_amr

  !=========================================================================
  ! AMR do_resonance2: two-line resonance (CII, Ca II H&K, etc.)
  !=========================================================================
  subroutine do_resonance2_amr(photon, grid, uz, xfreq_atom, cost, sint, &
                                S11, S22, S12, S33, S44)
    use define
    implicit none
    type(photon_type), intent(inout) :: photon
    type(grid_type),   intent(in)    :: grid
    real(wp),          intent(out)   :: uz, xfreq_atom, cost, sint
    real(wp), optional, intent(out)  :: S11, S22, S12, S33, S44

    integer  :: il
    real(wp) :: DnuHK, pH, pK, qH, qK, E1, cost2
    real(wp) :: va
    il   = photon%icell_amr
    va   = amr_grid%voigt_a(il)
    DnuHK = line%DnuHK_Hz / amr_grid%Dfreq(il)
    pH    = voigt(photon%xfreq + DnuHK, va) * one_over_three
    pK    = voigt(photon%xfreq,         va) * two_over_three
    pH    = pH / (pH + pK)

    if (rand_number() < pH) then
      ! scattered by an atom at H state (2P1/2)
      uz = rand_resonance_vz(photon%xfreq + DnuHK, va)
    else
      ! scattered by an atom at K state (2P3/2)
      uz = rand_resonance_vz(photon%xfreq, va)
    end if

    xfreq_atom = photon%xfreq - uz
    !--- E1 depends on frequency in the atom's rest frame
    !    qH = nu - nu_H, qK = nu - nu_K, nu_K = nu_H + DnuHK_Hz
    qH    = xfreq_atom + DnuHK
    qK    = xfreq_atom
    E1    = (2.0_wp*qK*qH + qH**2) / (qK**2 + 2.0_wp*qH**2)
    cost  = rand_resonance(E1)
    cost2 = cost**2
    sint  = sqrt(1.0_wp - cost2)

    if (present(S44)) then
      !--- Stokes scattering matrix: E2 = 1-E1, (3/2)E3 = (E1+2)/2
      S22 = three_over_four * E1 * (cost2 + 1.0_wp)
      S11 = S22 + (1.0_wp - E1)
      S12 = three_over_four * E1 * (cost2 - 1.0_wp)
      S33 = three_over_two  * E1 * cost
      S44 = (E1 + 2.0_wp)/2.0_wp * cost
    end if
    photon%E1 = E1
    photon%E2 = 1.0_wp - E1
    photon%E3 = (E1 + 2.0_wp)/3.0_wp
  end subroutine do_resonance2_amr

  !=========================================================================
  ! AMR do_resonance4: one upward transition, multiple downward transitions.
  ! (e.g. FeII 2374)
  !=========================================================================
  subroutine do_resonance4_amr(photon, grid, uz, xfreq_atom, cost, sint, &
                                S11, S22, S12, S33, S44)
    use define
    implicit none
    type(photon_type), intent(inout) :: photon
    type(grid_type),   intent(in)    :: grid
    real(wp),          intent(out)   :: uz, xfreq_atom, cost, sint
    real(wp), optional, intent(out)  :: S11, S22, S12, S33, S44

    integer  :: il, idown
    real(wp) :: va, cost2, del_xfreq

    il = photon%icell_amr
    va = amr_grid%voigt_a(il)

    uz         = rand_resonance_vz(photon%xfreq, va)
    xfreq_atom = photon%xfreq - uz

    idown = rand_alias_choise(line%b(1)%P_down, line%b(1)%A_down)
    if (idown /= 1) then
      del_xfreq  = line%b(1)%Elow_Hz(idown) / amr_grid%Dfreq(il)
      xfreq_atom = xfreq_atom - del_xfreq
    end if

    photon%E1 = line%b(1)%E1(idown)
    photon%E2 = line%b(1)%E2(idown)
    photon%E3 = line%b(1)%E3(idown)

    cost  = rand_resonance(photon%E1)
    cost2 = cost**2
    sint  = sqrt(1.0_wp - cost2)

    if (present(S44)) then
      S22 = three_over_four * photon%E1 * (cost2 + 1.0_wp)
      S11 = S22 + photon%E2
      S12 = three_over_four * photon%E1 * (cost2 - 1.0_wp)
      S33 = three_over_two  * photon%E1 * cost
      S44 = three_over_two  * photon%E3 * cost
    end if
  end subroutine do_resonance4_amr

  !=========================================================================
  ! AMR do_resonance5: two upward transitions, multiple downward transitions.
  ! (e.g. FeII multi-level)
  !=========================================================================
  subroutine do_resonance5_amr(photon, grid, uz, xfreq_atom, cost, sint, &
                                S11, S22, S12, S33, S44)
    use define
    implicit none
    type(photon_type), intent(inout) :: photon
    type(grid_type),   intent(in)    :: grid
    real(wp),          intent(out)   :: uz, xfreq_atom, cost, sint
    real(wp), optional, intent(out)  :: S11, S22, S12, S33, S44

    integer  :: il, iup, idown
    real(wp) :: va1, va2, Dx, p1, p2, cost2, del_xfreq

    il  = photon%icell_amr
    Dx  = line%delE_Hz(2) / amr_grid%Dfreq(il)
    va1 = amr_grid%voigt_a(il)
    va2 = amr_grid%voigt_a(il) * line%b(2)%damping / line%b(1)%damping
    p1  = voigt(photon%xfreq,      va1) * line%f12(1)
    p2  = voigt(photon%xfreq + Dx, va2) * line%f12(2)
    p1  = p1 / (p1 + p2)

    if (rand_number() < p1) then
      uz  = rand_resonance_vz(photon%xfreq,      va1)
      iup = 1
    else
      uz  = rand_resonance_vz(photon%xfreq + Dx, va2)
      iup = 2
    end if

    xfreq_atom = photon%xfreq - uz

    if (line%b(iup)%ndown > 1) then
      idown      = rand_alias_choise(line%b(iup)%P_down, line%b(iup)%A_down)
      del_xfreq  = line%b(iup)%Elow_Hz(idown) / amr_grid%Dfreq(il)
      xfreq_atom = xfreq_atom - del_xfreq
    else
      idown = 1
    end if

    photon%E1 = line%b(iup)%E1(idown)
    photon%E2 = line%b(iup)%E2(idown)
    photon%E3 = line%b(iup)%E3(idown)

    cost  = rand_resonance(photon%E1)
    cost2 = cost**2
    sint  = sqrt(1.0_wp - cost2)

    if (present(S44)) then
      S22 = three_over_four * photon%E1 * (cost2 + 1.0_wp)
      S11 = S22 + photon%E2
      S12 = three_over_four * photon%E1 * (cost2 - 1.0_wp)
      S33 = three_over_two  * photon%E1 * cost
      S44 = three_over_two  * photon%E3 * cost
    end if
  end subroutine do_resonance5_amr

  !=========================================================================
  ! AMR do_resonance6: three upward transitions, one downward each.
  ! (e.g. three-level system)
  !=========================================================================
  subroutine do_resonance6_amr(photon, grid, uz, xfreq_atom, cost, sint, &
                                S11, S22, S12, S33, S44)
    use define
    implicit none
    type(photon_type), intent(inout) :: photon
    type(grid_type),   intent(in)    :: grid
    real(wp),          intent(out)   :: uz, xfreq_atom, cost, sint
    real(wp), optional, intent(out)  :: S11, S22, S12, S33, S44

    integer  :: il, iup
    real(wp) :: va1, va2, va3, Dx2, Dx3
    real(wp) :: p1, p2, p3, ptot, pcum1, pcum2, xi, cost2

    il  = photon%icell_amr
    Dx2 = line%delE_Hz(2) / amr_grid%Dfreq(il)
    Dx3 = line%delE_Hz(3) / amr_grid%Dfreq(il)
    va1 = amr_grid%voigt_a(il)
    va2 = amr_grid%voigt_a(il) * line%b(2)%damping / line%b(1)%damping
    va3 = amr_grid%voigt_a(il) * line%b(3)%damping / line%b(1)%damping
    p1  = voigt(photon%xfreq,       va1) * line%f12(1)
    p2  = voigt(photon%xfreq + Dx2, va2) * line%f12(2)
    p3  = voigt(photon%xfreq + Dx3, va3) * line%f12(3)
    ptot  = p1 + p2 + p3
    pcum1 = p1 / ptot
    pcum2 = (p1 + p2) / ptot

    xi = rand_number()
    if (xi < pcum1) then
      uz  = rand_resonance_vz(photon%xfreq,       va1)
      iup = 1
    else if (xi < pcum2) then
      uz  = rand_resonance_vz(photon%xfreq + Dx2, va2)
      iup = 2
    else
      uz  = rand_resonance_vz(photon%xfreq + Dx3, va3)
      iup = 3
    end if

    xfreq_atom = photon%xfreq - uz

    photon%E1 = line%b(iup)%E1(1)
    photon%E2 = line%b(iup)%E2(1)
    photon%E3 = line%b(iup)%E3(1)

    cost  = rand_resonance(photon%E1)
    cost2 = cost**2
    sint  = sqrt(1.0_wp - cost2)

    if (present(S44)) then
      S22 = three_over_four * photon%E1 * (cost2 + 1.0_wp)
      S11 = S22 + photon%E2
      S12 = three_over_four * photon%E1 * (cost2 - 1.0_wp)
      S33 = three_over_two  * photon%E1 * cost
      S44 = three_over_two  * photon%E3 * cost
    end if
  end subroutine do_resonance6_amr

end module scattering_amr_mod
