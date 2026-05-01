module line_mod
  use define
  use voigt_mod
  use random
contains

!++++++++++++++++++++++++++++++++++++++
  function calc_voigt1(grid,xfreq,icell,jcell,kcell) result(voigt1)
  use define
  implicit none
  type(grid_type),   intent(in) :: grid
  real(kind=wp),     intent(in) :: xfreq
  integer,           intent(in) :: icell,jcell,kcell
  real(kind=wp) :: voigt1
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt
  voigt1 = voigt(xfreq, grid%voigt_a(icell,jcell,kcell))
  end function calc_voigt1
!---
  function calc_voigt2(grid,xfreq,icell,jcell,kcell) result(voigt2)
  use define
  implicit none
  type(grid_type),   intent(in) :: grid
  real(kind=wp),     intent(in) :: xfreq
  integer,           intent(in) :: icell,jcell,kcell
  real(kind=wp) :: DnuHK, voigt2
  real(kind=wp), parameter :: fraction1 = 1.0_wp/3.0_wp, &
                              fraction2 = 2.0_wp/3.0_wp
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt
  DnuHK  = line%DnuHK_Hz/grid%Dfreq(icell,jcell,kcell)
  voigt2 = voigt(xfreq+DnuHK, grid%voigt_a(icell,jcell,kcell))*fraction1 + &
           voigt(xfreq,       grid%voigt_a(icell,jcell,kcell))*fraction2
  end function calc_voigt2
!---
  function calc_voigt3(grid,xfreq,icell,jcell,kcell) result(voigt3)
  use define
  implicit none
  type(grid_type),   intent(in) :: grid
  real(kind=wp),     intent(in) :: xfreq
  integer,           intent(in) :: icell,jcell,kcell
  real(kind=wp) :: Dnu, voigt3, a_ratio, f_ratio
  integer       :: i
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt
  voigt3 = voigt(xfreq, grid%voigt_a(icell,jcell,kcell))
  do i=2, line%nup
     Dnu     = line%delE_Hz(i)   /grid%Dfreq(icell,jcell,kcell)
     a_ratio = line%b(i)%damping /line%b(1)%damping
     f_ratio = line%f12(i)       /line%f12(1)
     voigt3  = voigt3 + voigt(xfreq+Dnu, grid%voigt_a(icell,jcell,kcell)*a_ratio)*f_ratio
  enddo
  end function calc_voigt3
!++++++++++++++++++++++++++++++++++++++
  function calc_voigt_HD(grid,xfreq,icell,jcell,kcell) result(voigt_HD)
  !---
  !-- Combined H + D Lyman-α Voigt profile (line_type = 7).
  !-- Returns the dimensionless profile that, when multiplied by the
  !-- cell-stored rhokap (built from H atomic data), gives the total
  !-- (H + D) line opacity at the photon's frequency.
  !-- Photon xfreq is carried in H-frame Doppler units (H line center, H Dfreq).
  !-- The deuterium contribution is computed on-the-fly using the cross-species
  !-- constants precomputed at line setup; no per-cell D arrays are stored.
  !---
  use define
  implicit none
  type(grid_type),   intent(in) :: grid
  real(kind=wp),     intent(in) :: xfreq
  integer,           intent(in) :: icell,jcell,kcell
  real(kind=wp) :: voigt_HD
  real(kind=wp) :: voigt_a_H, dx_HD, xfreq_D, voigt_H, voigt_D
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt
  voigt_a_H = grid%voigt_a(icell,jcell,kcell)
  dx_HD     = line%delta_nu_HD_Hz / grid%Dfreq(icell,jcell,kcell)
  xfreq_D   = (xfreq - dx_HD) * line%ratio_Dfreq_HD
  voigt_H   = voigt(xfreq,   voigt_a_H)
  voigt_D   = voigt(xfreq_D, voigt_a_H * line%ratio_voigta_HD)
  voigt_HD  = voigt_H + line%nD_over_nH * line%ratio_Dfreq_HD * voigt_D
  end function calc_voigt_HD
!++++++++++++++++++++++++++++++++++++++
  subroutine do_resonance1(photon,grid,uz,xfreq_atom,cost,sint,S11,S22,S12,S33,S44)
  use define
  implicit none
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(in)    :: grid
  real(kind=wp),     intent(out)   :: uz,xfreq_atom, cost,sint
  real(kind=wp), optional, intent(out) :: S11,S22,S12,S33,S44
  real(kind=wp) :: cost2
  !--- Select an atom by which the photon is scattered.
  uz = rand_resonance_vz(photon%xfreq, grid%voigt_a(photon%icell,photon%jcell,photon%kcell))

  !--- xfreq_atom = frequency in the atom's rest frame.
  xfreq_atom = photon%xfreq - uz

  photon%E1 = line%E1
  photon%E2 = line%E2
  photon%E3 = line%E3

  !--- Select a new scattering angle (cos(theta)).
  cost  = rand_resonance(photon%E1)
  cost2 = cost**2
  sint  = sqrt(1.0_wp-cost2)

  !--- Calculate Scattering Matrix. (S11, S12, S33, and S44)
  if (present(S44)) then
     S22 = three_over_four * photon%E1 * (cost2 + 1.0_wp)
     S11 = S22 + photon%E2
     S12 = three_over_four * photon%E1 * (cost2 - 1.0_wp)
     S33 = three_over_two  * photon%E1 * cost
     S44 = three_over_two  * photon%E3 * cost
  endif
  end subroutine do_resonance1
!++++++++++++++++++++++++++++++++++++++
  subroutine do_resonance2(photon,grid,uz,xfreq_atom,cost,sint,S11,S22,S12,S33,S44)
  use define
  implicit none
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(in)    :: grid
  real(kind=wp),     intent(out)   :: uz,xfreq_atom, cost,sint
  real(kind=wp), optional, intent(out) :: S11,S22,S12,S33,S44
  real(kind=wp) :: DnuHK, pK, pH
  real(kind=wp) :: qH, qK, E1
  real(kind=wp) :: cost2

  !--- Select an atom by which the photon is scattered.
  DnuHK = line%DnuHK_Hz/grid%Dfreq(photon%icell,photon%jcell,photon%kcell)
  pH    = voigt(photon%xfreq+DnuHK, grid%voigt_a(photon%icell,photon%jcell,photon%kcell)) * one_over_three
  pK    = voigt(photon%xfreq,       grid%voigt_a(photon%icell,photon%jcell,photon%kcell)) * two_over_three
  pH    = pH/(pH + pK)
  if (rand_number() < pH) then
     !--- scattered by an atom at H state (2P1/2).
     uz = rand_resonance_vz(photon%xfreq + DnuHK, grid%voigt_a(photon%icell,photon%jcell,photon%kcell))
  else
     !--- scattered by an atom at K state (2P3/2).
     uz = rand_resonance_vz(photon%xfreq,         grid%voigt_a(photon%icell,photon%jcell,photon%kcell))
  endif

  !--- xfreq_atom = frequency in the atom's rest frame.
  xfreq_atom = photon%xfreq - uz

  !--- Select a new scattering angle (cos(theta)).
  !--- E1 is a function of frequency in the atom's rest frame.
  !--- qH = nu - nu_H, qK = nu - nu_K, nu_K = nu_H + DnuHK_Hz
  !--- qK = xfreq_atom, qH = qK + DnuHK
  qH    = xfreq_atom + DnuHK
  qK    = xfreq_atom
  E1    = (2.0_wp*qK*qH + qH**2)/(qK**2 + 2.0_wp*qH**2)
  cost  = rand_resonance(E1)
  cost2 = cost**2
  sint  = sqrt(1.0_wp-cost2)

  if (present(S44)) then
     !--- Calculate Stokes parameters. (S11, S12, S33, and S44)
     !--- Note S44 /= S33 and S34 = 0.
     !--- E2 = 1-E1 and E3 = (E1+2)/3. (3/2) E3 = (E1+2)/2.
     S22 = three_over_four * E1 * (cost2 + 1.0_wp)
     S11 = S22 + (1.0_wp - E1)
     S12 = three_over_four * E1 * (cost2 - 1.0_wp)
     S33 = three_over_two  * E1 * cost
     S44 = (E1 + 2.0_wp)/2.0_wp * cost
  endif
  photon%E1 = E1
  photon%E2 = 1.0_wp - E1
  photon%E3 = (E1 + 2.0_wp)/3.0_wp
  end subroutine do_resonance2
!++++++++++++++++++++++++++++++++++++++
  subroutine do_resonance3(photon,grid,uz,xfreq_atom,cost,sint,S11,S22,S12,S33,S44)
  use define
  implicit none
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(in)    :: grid
  real(kind=wp),     intent(out)   :: uz,xfreq_atom, cost,sint
  real(kind=wp), optional, intent(out) :: S11,S22,S12,S33,S44
  real(kind=wp) :: Dx, p1, p2, va1, va2
  real(kind=wp) :: cost2, E1, del_xfreq
  integer       :: iup

  !--- Select which upward transition will occur.
  Dx  = line%delE_Hz(2)/grid%Dfreq(photon%icell,photon%jcell,photon%kcell)
  va1 = grid%voigt_a(photon%icell,photon%jcell,photon%kcell)
  va2 = grid%voigt_a(photon%icell,photon%jcell,photon%kcell) * line%b(2)%damping/line%b(1)%damping
  p1  = voigt(photon%xfreq,    va1) * line%f12(1)
  p2  = voigt(photon%xfreq+Dx, va2) * line%f12(2)
  p1  = p1/(p2 + p1)

  !--- Select an atom by which the photon is scattered.
  if (rand_number() < p1) then
     !--- scattered by an atom at 1 state
     uz  = rand_resonance_vz(photon%xfreq,      va1)
     iup = 1
  else
     !--- scattered by an atom at 2 state
     uz  = rand_resonance_vz(photon%xfreq + Dx, va2)
     iup = 2
  endif

  !--- xfreq_atom = frequency in the atom's rest frame.
  xfreq_atom = photon%xfreq - uz

  !--- only one downward transition for each absorption.
  !--- This part have to be updated.
  photon%E1 = line%b(iup)%E1(1)
  photon%E2 = line%b(iup)%E2(1)
  photon%E3 = line%b(iup)%E3(1)

  !--- Select a new scattering angle (cos(theta)).
  cost  = rand_resonance(photon%E1)
  cost2 = cost**2
  sint  = sqrt(1.0_wp-cost2)

  if (present(S44)) then
     !--- Calculate Stokes parameters. (S11, S12, S33, and S44)
     !--- Note S44 /= S33 and S34 = 0.
     S22 = three_over_four * photon%E1 * (cost2 + 1.0_wp)
     S11 = S22 + photon%E2
     S12 = three_over_four * photon%E1 * (cost2 - 1.0_wp)
     S33 = three_over_two  * photon%E1 * cost
     S44 = three_over_two  * photon%E3 * cost
  endif
  end subroutine do_resonance3
!++++++++++++++++++++++++++++++++++++++
  subroutine do_resonance4(photon,grid,uz,xfreq_atom,cost,sint,S11,S22,S12,S33,S44)
  use define
  implicit none
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(in)    :: grid
  real(kind=wp),     intent(out)   :: uz,xfreq_atom, cost,sint
  real(kind=wp), optional, intent(out) :: S11,S22,S12,S33,S44
  real(kind=wp) :: cost2, del_xfreq
  integer :: idown
  !--- Select an atom by which the photon is scattered.
  uz = rand_resonance_vz(photon%xfreq, grid%voigt_a(photon%icell,photon%jcell,photon%kcell))

  !--- xfreq_atom = frequency in the atom's rest frame.
  xfreq_atom = photon%xfreq - uz

  idown      = rand_alias_choise(line%b(1)%P_down, line%b(1)%A_down)
  if (idown /= 1) then
     !--- bug-fixed (2023.07.25)
     !del_xfreq  = line%b(1)%Elow_Hz(2)/grid%Dfreq(photon%icell,photon%jcell,photon%kcell)
     del_xfreq  = line%b(1)%Elow_Hz(idown)/grid%Dfreq(photon%icell,photon%jcell,photon%kcell)
     xfreq_atom = xfreq_atom - del_xfreq
  endif

  photon%E1 = line%b(1)%E1(idown)
  photon%E2 = line%b(1)%E2(idown)
  photon%E3 = line%b(1)%E3(idown)

  !--- Select a new scattering angle (cos(theta)).
  cost  = rand_resonance(photon%E1)
  cost2 = cost**2
  sint  = sqrt(1.0_wp-cost2)

  !--- Calculate Scattering Matrix. (S11, S12, S33, and S44)
  if (present(S44)) then
     !--- Calculate Stokes parameters. (S11, S12, S33, and S44)
     !--- Note S44 /= S33 and S34 = 0.
     S22 = three_over_four * photon%E1* (cost2 + 1.0_wp)
     S11 = S22 + photon%E2
     S12 = three_over_four * photon%E1 * (cost2 - 1.0_wp)
     S33 = three_over_two  * photon%E1 * cost
     S44 = three_over_two  * photon%E3 * cost
  endif
  end subroutine do_resonance4
!++++++++++++++++++++++++++++++++++++++
  subroutine do_resonance5(photon,grid,uz,xfreq_atom,cost,sint,S11,S22,S12,S33,S44)
  use define
  implicit none
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(in)    :: grid
  real(kind=wp),     intent(out)   :: uz,xfreq_atom, cost,sint
  real(kind=wp), optional, intent(out) :: S11,S22,S12,S33,S44
  real(kind=wp) :: Dx, p1, p2, va1, va2
  real(kind=wp) :: cost2, E1, del_xfreq
  integer       :: iup, idown

  !--- Select which upward transition will occur.
  Dx  = line%delE_Hz(2)/grid%Dfreq(photon%icell,photon%jcell,photon%kcell)
  va1 = grid%voigt_a(photon%icell,photon%jcell,photon%kcell)
  va2 = grid%voigt_a(photon%icell,photon%jcell,photon%kcell) * line%b(2)%damping/line%b(1)%damping
  p1  = voigt(photon%xfreq,    va1) * line%f12(1)
  p2  = voigt(photon%xfreq+Dx, va2) * line%f12(2)
  p1  = p1/(p2 + p1)

  !--- Select an atom by which the photon is scattered.
  if (rand_number() < p1) then
     !--- scattered by an atom at 1 state
     uz  = rand_resonance_vz(photon%xfreq,      va1)
     iup = 1
  else
     !--- scattered by an atom at 2 state
     uz  = rand_resonance_vz(photon%xfreq + Dx, va2)
     iup = 2
  endif

  !--- xfreq_atom = frequency in the atom's rest frame.
  xfreq_atom = photon%xfreq - uz

  !--- Select which downward transition will occur.
  if (line%b(iup)%ndown > 1) then
     idown      = rand_alias_choise(line%b(iup)%P_down, line%b(iup)%A_down)
     del_xfreq  = line%b(iup)%Elow_Hz(idown)/grid%Dfreq(photon%icell,photon%jcell,photon%kcell)
     xfreq_atom = xfreq_atom - del_xfreq
  else
     idown = 1
  endif
  photon%E1 = line%b(iup)%E1(idown)
  photon%E2 = line%b(iup)%E2(idown)
  photon%E3 = line%b(iup)%E3(idown)

  !--- Select a new scattering angle (cos(theta)).
  cost  = rand_resonance(photon%E1)
  cost2 = cost**2
  sint  = sqrt(1.0_wp-cost2)

  if (present(S44)) then
     !--- Calculate Stokes parameters. (S11, S12, S33, and S44)
     !--- Note S44 /= S33 and S34 = 0.
     S22 = three_over_four * photon%E1 * (cost2 + 1.0_wp)
     S11 = S22 + photon%E2
     S12 = three_over_four * photon%E1 * (cost2 - 1.0_wp)
     S33 = three_over_two  * photon%E1 * cost
     S44 = three_over_two  * photon%E3 * cost
  endif
  end subroutine do_resonance5
!++++++++++++++++++++++++++++++++++++++
  subroutine do_resonance6(photon,grid,uz,xfreq_atom,cost,sint,S11,S22,S12,S33,S44)
  use define
  implicit none
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(in)    :: grid
  real(kind=wp),     intent(out)   :: uz,xfreq_atom, cost,sint
  real(kind=wp), optional, intent(out) :: S11,S22,S12,S33,S44
  real(kind=wp) :: Dx2, Dx3, p1, p2, p3, ptot, pcum1, pcum2, xi, va1, va2, va3
  real(kind=wp) :: cost2, E1, del_xfreq
  integer       :: iup, idown

  !--- Select which upward transition will occur.
  Dx2 = line%delE_Hz(2)/grid%Dfreq(photon%icell,photon%jcell,photon%kcell)
  Dx3 = line%delE_Hz(3)/grid%Dfreq(photon%icell,photon%jcell,photon%kcell)
  va1 = grid%voigt_a(photon%icell,photon%jcell,photon%kcell)
  va2 = grid%voigt_a(photon%icell,photon%jcell,photon%kcell) * line%b(2)%damping/line%b(1)%damping
  va3 = grid%voigt_a(photon%icell,photon%jcell,photon%kcell) * line%b(3)%damping/line%b(1)%damping
  p1  = voigt(photon%xfreq,     va1) * line%f12(1)
  p2  = voigt(photon%xfreq+Dx2, va2) * line%f12(2)
  p3  = voigt(photon%xfreq+Dx3, va3) * line%f12(3)
  ptot  = p1 + p2 + p3
  pcum1 = p1/ptot
  pcum2 = (p1 + p2)/ptot

  !--- Select an atom by which the photon is scattered.
  xi = rand_number()
  if (xi < pcum1) then
     !--- scattered by an atom at 1 state
     uz  = rand_resonance_vz(photon%xfreq,       va1)
     iup = 1
  else if (xi < pcum2) then
     !--- scattered by an atom at 2 state
     uz  = rand_resonance_vz(photon%xfreq + Dx2, va2)
     iup = 2
  else
     !--- scattered by an atom at 3 state
     uz  = rand_resonance_vz(photon%xfreq + Dx3, va3)
     iup = 3
  endif

  !--- xfreq_atom = frequency in the atom's rest frame.
  xfreq_atom = photon%xfreq - uz

  photon%E1 = line%b(iup)%E1(1)
  photon%E2 = line%b(iup)%E2(1)
  photon%E3 = line%b(iup)%E3(1)

  !--- Select a new scattering angle (cos(theta)).
  cost  = rand_resonance(photon%E1)
  cost2 = cost**2
  sint  = sqrt(1.0_wp-cost2)

  if (present(S44)) then
     !--- Calculate Stokes parameters. (S11, S12, S33, and S44)
     !--- Note S44 /= S33 and S34 = 0.
     S22 = three_over_four * photon%E1 * (cost2 + 1.0_wp)
     S11 = S22 + photon%E2
     S12 = three_over_four * photon%E1 * (cost2 - 1.0_wp)
     S33 = three_over_two  * photon%E1 * cost
     S44 = three_over_two  * photon%E3 * cost
  endif
  end subroutine do_resonance6
!++++++++++++++++++++++++++++++++++++++
  subroutine do_resonance_HD(photon,grid,uz,xfreq_atom,cost,sint,S11,S22,S12,S33,S44)
  !---
  !-- Combined H + D Lyman-α resonance scattering (line_type = 7).
  !-- Stochastically selects H or D based on the local partial opacity, then
  !-- samples the atom velocity in the chosen species' Doppler frame and
  !-- converts the result back to the H-frame Doppler units carried by
  !-- photon%xfreq.
  !-- v1 limitation: the caller's perpendicular-velocity sampling and recoil
  !-- step in scatter_resonance_*nostokes/stokes use H atomic constants for
  !-- both species. For D scatters this introduces a small error in the
  !-- thermal kick width and recoil shift; it is acceptable because D
  !-- scatters are rare at cosmic D/H abundance and the dominant D feature
  !-- (the absorption dip) depends on opacity, not on the re-emission shape.
  !---
  use define
  implicit none
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(in)    :: grid
  real(kind=wp),     intent(out)   :: uz, xfreq_atom, cost, sint
  real(kind=wp), optional, intent(out) :: S11, S22, S12, S33, S44
  real(kind=wp) :: voigt_a_H, dx_HD, xfreq_D
  real(kind=wp) :: pH, pD_term, pH_norm
  real(kind=wp) :: uz_D, xfreq_atom_D
  real(kind=wp) :: cost2
!DIR$ ATTRIBUTES INLINE :: voigt_mod_voigt
  voigt_a_H = grid%voigt_a(photon%icell,photon%jcell,photon%kcell)
  dx_HD     = line%delta_nu_HD_Hz / grid%Dfreq(photon%icell,photon%jcell,photon%kcell)
  xfreq_D   = (photon%xfreq - dx_HD) * line%ratio_Dfreq_HD

  !--- Partial opacity contributions: H term, D term (with abundance + Doppler-ratio prefactor).
  pH      = voigt(photon%xfreq, voigt_a_H)
  pD_term = line%nD_over_nH * line%ratio_Dfreq_HD * &
            voigt(xfreq_D, voigt_a_H * line%ratio_voigta_HD)
  pH_norm = pH / (pH + pD_term)

  if (rand_number() < pH_norm) then
     !--- Hydrogen scatter: identical to do_resonance1 path.
     uz         = rand_resonance_vz(photon%xfreq, voigt_a_H)
     xfreq_atom = photon%xfreq - uz
  else
     !--- Deuterium scatter: sample atom velocity in D's Doppler frame, then
     !-- convert uz and xfreq_atom back to H-frame Doppler units for downstream.
     uz_D         = rand_resonance_vz(xfreq_D, voigt_a_H * line%ratio_voigta_HD)
     xfreq_atom_D = xfreq_D - uz_D
     uz           = uz_D / line%ratio_Dfreq_HD
     xfreq_atom   = photon%xfreq - uz
  endif

  !--- Phase function: simple dipole (E1=1, E2=0, E3=1) for both species
  !--  in v1 (no fine-structure splitting).
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
  endif
  end subroutine do_resonance_HD
!++++++++++++++++++++++++++++++++++++++
  subroutine setup_resonance_line()
  use define
  use random
  implicit none
  real(kind=wp), parameter :: speedc_cm   = 2.99792458e10_wp
  real(kind=wp), parameter :: sigma_0     = 0.026540083434_wp
  real(kind=wp), parameter :: amu         = 1.67262192e-24_wp
  real(kind=wp), parameter :: vtherm1_amu = 0.12895319011972164_wp
  real(kind=wp) :: mass_amu
  integer :: i
  !-- line id is defined using the shorter wavelength of the doublet.
  !-- the central wavelength is the harmonic mean of the two wavelengths.
  !-- the data are obtained from physics.nist.gov.
  !-- cross0  = cross-section at line center in units of cm^2
  !-- sigma_0 = pi x e^2/(m_e c) = pi x c x R_e = 0.026540083434 (cm^2 Hz)
  !-- vtherm1_amu = sqrt(2*1K/1amu) = thermal velocity of a particle with 1 amu mass at 1 K, in units of km/s.
  !-- atomic mass: https://www.angelo.edu/faculty/kboudrea/periodic/structure_mass.htm

  !-- line types
  !   1 : singlet (one resonance)
  !   2 : doublet (two upward transitions, two resonances)
  !   3 : two resonance lines with different f12 and damping constants. (two upward transitions)
  !   4 : one upward transition + several downward transitions
  !       (one resonance + one or more fluorescences) (SiII 1527, FeII UV3 lines)
  !   5 : FeII UV1 or UV2 like lines (the most general case)
  !   6 : three upward transitions + one downward transition
  !   7 : H + D Lyman-α (two coexisting two-level resonance scatterers; no fine structure)

  !-- Convenience flag: par%include_deuterium = .true. with line_id = 'ly_alpha'
  !-- promotes the run to 'ly_alpha_HD'. Fine-structure splitting is not yet
  !-- supported for the H+D mode: warn and disable if both flags are set.
  if (par%include_deuterium .and. trim(par%line_id) == 'ly_alpha') then
     par%line_id = 'ly_alpha_HD'
  endif
  if (trim(par%line_id) == 'ly_alpha_HD' .and. par%fine_structure) then
     write(*,'(a)') 'WARNING: fine_structure is not supported for ly_alpha_HD; disabling.'
     par%fine_structure = .false.
  endif

  select case (trim(par%line_id))
  case ('CIV_1548')
     !-- vacuum wavelengths 1548.187 (short, K, S1/2-P3/2), 1550.772 (long, H, S1/2-P1/2)
     !-- observed wavelengths 1548.202, 1550.774
     !-- Aki = 2.65e8 (short), 2.64e8 (long), fik = 1.90e-1 (short), 9.52e-2 (long)
     line%ion_id    = 'C IV'
     line%line_type = 2
     line%wavelength0 = 0.1548187_wp
     !line%f12(1)    = 0.2852_wp
     line%f12(1:2)  = [0.190_wp, 0.0952_wp]
     line%damping   = 2.647e8_wp
     line%cross0    = sigma_0/sqrt(pi)*line%f12(1)
     mass_amu       = 12.011_wp
     line%vtherm1   = vtherm1_amu/sqrt(mass_amu)
     line%DnuHK_Hz  = speedc_cm * (64591.7_wp - 64484.0_wp)
     line%g_recoil0 = (h_planck/amu/line%mass_amu)/(line%wavelength0*um2m)**2
  case ('NV_1239')
     !-- vacuum wavelengths 1238.821 (short, K, S1/2-P3/2), 1242.804 (long, H, S1/2-P1/2)
     !-- observed wavelengths 1238.921, 1242.804
     !-- Aki = 3.40e8 (short), 3.37e8 (long), fik = 1.56e-1 (short), 7.80e-2 (long)
     line%ion_id    = 'N V'
     line%line_type = 2
     line%wavelength0 = 0.1238821_wp
     !line%f12(1)    = 0.2340_wp
     line%f12(1:2)  = [0.156_wp, 0.078_wp]
     line%damping   = 3.390e8_wp
     line%cross0    = sigma_0/sqrt(pi)*line%f12(1)
     mass_amu       = 14.0067_wp
     line%vtherm1   = vtherm1_amu/sqrt(mass_amu)
     line%DnuHK_Hz  = speedc_cm * (80721.9_wp - 80463.2_wp)
     line%g_recoil0 = (h_planck/amu/mass_amu)/(line%wavelength0*um2m)**2
  case ('OVI_1032')
     !-- vacuum wavelengths 1031.912 (short, K, S1/2-P3/2), 1037.613 (long, H, S1/2-P1/2)
     !-- observed wavelengths ??? (not shown in NIST)
     !-- Aki = 4.16e8 (short), 4.09e8 (long), fik = 1.33e-1 (short), 6.60e-2 (long)
     line%ion_id    = 'O VI'
     line%line_type = 2
     line%wavelength0 = 0.1031912_wp
     !line%f12(1)    = 0.1990_wp
     line%f12(1:2)  = [0.133_wp, 0.066_wp]
     line%damping   = 4.137e8_wp
     line%cross0    = sigma_0/sqrt(pi)*line%f12(1)
     mass_amu       = 15.9994_wp
     line%vtherm1   = vtherm1_amu/sqrt(mass_amu)
     line%DnuHK_Hz  = speedc_cm * (96907.5_wp - 96375.0_wp)
     line%g_recoil0 = (h_planck/amu/mass_amu)/(line%wavelength0*um2m)**2
  case ('NaI_D')
     !-- vacuum wavelengths 5891.583253 (short, K, S1/2-P3/2), 5897.558147 (long, H, S1/2-P1/2)
     !-- observed wavelengths 5891.583264 (short), 5897.558147 (long)
     !-- Aki = 6.16e7 (short), 6.14e7 (long), fik = 6.41e-1 (short), 3.20e-1 (long)
     line%ion_id    = 'Na I'
     line%line_type = 2
     line%wavelength0 = 0.5891583253_wp
     !line%f12(1)    = 0.9610_wp
     line%f12(1:2)  = [0.641_wp, 0.320_wp]
     line%damping   = 6.153e7_wp
     line%cross0    = sigma_0/sqrt(pi)*line%f12(1)
     mass_amu       = 22.98977_wp
     line%vtherm1   = vtherm1_amu/sqrt(mass_amu)
     line%DnuHK_Hz  = speedc_cm * (16973.36619_wp - 16956.17025_wp)
     line%g_recoil0 = (h_planck/amu/mass_amu)/(line%wavelength0*um2m)**2
  case ('CaII_HK')
     !-- vacuum wavelengths 3934.777 (short, K, S1/2-P3/2), 3969.591 (long, H, S1/2-P1/2)
     !-- observed wavelengths 3934.77 (short), 3969.59 (long)
     !-- Aki = 1.47e8 (short), 1.4e8 (long), fik = 6.82e-01 (short), 3.3e-01 (long)
     line%ion_id    = 'Ca II'
     line%line_type = 2
     line%wavelength0 = 0.3934777_wp
     !line%f12(1)    = 1.012_wp
     line%f12(1:2)  = [0.682_wp, 0.330_wp]
     line%damping   = 1.446667e8_wp
     line%cross0    = sigma_0/sqrt(pi)*line%f12(1)
     mass_amu       = 40.078_wp
     line%vtherm1   = vtherm1_amu/sqrt(mass_amu)
     line%DnuHK_Hz  = speedc_cm * (25414.40_wp - 25191.51_wp)
     line%g_recoil0 = (h_planck/amu/mass_amu)/(line%wavelength0*um2m)**2
  case ('MgII_2796')
     !-- vacuum wavelengths 2796.352 (short, K, S1/2-P3/2), 2803.531 (long, H, S1/2-P1/2)
     !-- observed wavelengths 2796.352 (short), 2803.530 (long)
     !-- Aki = 2.60e8 (short), 2.57e8 (long), fik = 6.08e-1 (short), 3.03e-1 (long)
     line%ion_id    = 'Mg II'
     line%line_type = 2
     line%wavelength0 = 0.2796352_wp
     !line%f12(1)    = 0.9110_wp
     line%f12(1:2)  = [0.608_wp, 0.303_wp]
     line%damping   = 2.590e8_wp
     line%cross0    = sigma_0/sqrt(pi)*line%f12(1)
     mass_amu       = 24.305_wp
     line%vtherm1   = vtherm1_amu/sqrt(mass_amu)
     line%DnuHK_Hz  = speedc_cm * (35760.88_wp - 35669.31_wp)
     line%g_recoil0 = (h_planck/amu/mass_amu)/(line%wavelength0*um2m)**2
  case ('SiIV_1394')
     !-- vacuum wavelengths 1393.755 (short, K, S1/2-P3/2), 1402.770 (long, H, S1/2-P1/2)
     !-- observed wavelengths 1393.76 (short), 1402.77 (long)
     !-- Aki = 8.80e8 (short), 8.63e8 (long), fik = 5.13e-1 (short), 2.55e-1 (long)
     line%ion_id    = 'Si IV'
     line%line_type = 2
     line%wavelength0 = 0.1393755_wp
     !line%f12(1)    = 0.768_wp
     line%f12(1:2)  = [0.513_wp, 0.255_wp]
     line%damping   = 8.743e8_wp
     line%cross0    = sigma_0/sqrt(pi)*line%f12(1)
     mass_amu       = 28.0855_wp
     line%vtherm1   = vtherm1_amu/sqrt(mass_amu)
     line%DnuHK_Hz  = speedc_cm * (71748.64_wp - 71287.54_wp)
     line%g_recoil0 = (h_planck/amu/mass_amu)/(line%wavelength0*um2m)**2
  case ('AlII_1671')
     !-- vacuum wavelengths 1670.7874 (singlet, S0-P1)
     !-- observed wavelengths 1670.7867
     !-- Aki = 1.41e9, fik = 1.77
     line%ion_id    = 'Al II'
     line%line_type = 1
     line%wavelength0 = 0.16707874_wp
     line%f12(1)    = 1.77_wp
     line%damping   = 1.41e9_wp
     line%cross0    = sigma_0/sqrt(pi)*line%f12(1)
     mass_amu       = 26.98154
     line%vtherm1   = vtherm1_amu/sqrt(mass_amu)
     line%DnuHK_Hz  = 0.0_wp
     line%g_recoil0 = (h_planck/amu/mass_amu)/(line%wavelength0*um2m)**2
     line%E1        = 1.0_wp
     line%E2        = 0.0_wp
     line%E3        = 1.0_wp
  case ('SiII_1527')
     !-- vacuum wavelengths 1527.707, 1533.431 (resonance + fluorescence)
     !-- observed wavelengths 1526.72, 1533.45
     !-- Aki = 3.81e8, fik = 0.133 / Aki = 7.52e8, fik = 0.133
     line%ion_id    = 'Si II'
     line%line_type = 4
     line%wavelength0 = 0.1526707_wp
     line%f12(1)    = 0.133_wp

     line%cross0    = sigma_0/sqrt(pi)*line%f12(1)
     mass_amu       = 28.0855_wp
     line%vtherm1   = vtherm1_amu/sqrt(mass_amu)
     line%g_recoil0 = (h_planck/amu/mass_amu)/(line%wavelength0*um2m)**2
     line%nup       = 1
     allocate(line%b(line%nup))

     line%b(1)%ndown = 2
     allocate(line%b(1)%A21(line%b(1)%ndown))
     allocate(line%b(1)%Elow_Hz(line%b(1)%ndown))
     allocate(line%b(1)%E1(line%b(1)%ndown))
     allocate(line%b(1)%E2(line%b(1)%ndown))
     allocate(line%b(1)%E3(line%b(1)%ndown))
     allocate(line%b(1)%P_down(line%b(1)%ndown))
     allocate(line%b(1)%A_down(line%b(1)%ndown))
     line%b(1)%A21(1:2)     = [3.81e8_wp,   7.52e8_wp]
     line%b(1)%Elow_Hz(1:2) = [0.0_wp,    287.24_wp] * speedc_cm
     line%b(1)%E1(1:2)      = [0.0_wp,         0.0_wp]
     line%b(1)%E2(1:2)      = [1.0_wp,         1.0_wp]
     line%b(1)%E3(1:2)      = [2.0_wp/3.0_wp, -1.0_wp/3.0_wp]
     line%b(1)%damping      = sum(line%b(1)%A21(:))
     line%b(1)%P_down(:)    = line%b(1)%A21(:)/line%b(1)%damping
     call random_alias_setup(line%b(1)%P_down, line%b(1)%A_down)
     line%damping   = line%b(1)%damping
  case ('SiII_1260')
     line%ion_id    = 'Si II'
     line%line_type = 4
     line%wavelength0 = 0.1260422_wp
     line%f12(1)    = 1.22_wp

     line%cross0    = sigma_0/sqrt(pi)*line%f12(1)
     mass_amu       = 28.0855_wp
     line%vtherm1   = vtherm1_amu/sqrt(mass_amu)
     line%g_recoil0 = (h_planck/amu/mass_amu)/(line%wavelength0*um2m)**2
     line%nup       = 1
     allocate(line%b(line%nup))

     line%b(1)%ndown = 2
     allocate(line%b(1)%A21(line%b(1)%ndown))
     allocate(line%b(1)%Elow_Hz(line%b(1)%ndown))
     allocate(line%b(1)%E1(line%b(1)%ndown))
     allocate(line%b(1)%E2(line%b(1)%ndown))
     allocate(line%b(1)%E3(line%b(1)%ndown))
     allocate(line%b(1)%P_down(line%b(1)%ndown))
     allocate(line%b(1)%A_down(line%b(1)%ndown))
     line%b(1)%A21(1:2)     = [2.57e9_wp,   4.73e8_wp]
     line%b(1)%Elow_Hz(1:2) = [0.0_wp,    287.24_wp] * speedc_cm
     line%b(1)%E1(1:2)      = [1.0_wp/2.0_wp, -2.0_wp/5.0_wp]
     line%b(1)%E2(1:2)      = [1.0_wp/2.0_wp,  7.0_wp/5.0_wp]
     line%b(1)%E3(1:2)      = [5.0_wp/6.0_wp,  1.0_wp/3.0_wp]
     line%b(1)%damping      = sum(line%b(1)%A21(:))
     line%b(1)%P_down(:)    = line%b(1)%A21(:)/line%b(1)%damping
     call random_alias_setup(line%b(1)%P_down, line%b(1)%A_down)
     line%damping   = line%b(1)%damping
  case ('SiII_1193', 'SiII_1190')
     line%ion_id    = 'Si II'
     line%line_type = 5
     line%wavelength0 = 0.1193290_wp
     ! the central frequency is assumed to be that with the highest f12.
     line%f12(1:2)  = [0.575_wp, 0.277_wp]
     line%cross0    = sigma_0/sqrt(pi)*line%f12(1)
     mass_amu       = 28.0855_wp
     line%vtherm1   = vtherm1_amu/sqrt(mass_amu)
     line%g_recoil0 = (h_planck/amu/mass_amu)/(line%wavelength0*um2m)**2

     !--- nbranch_up   = number of upward transitions
     !--- nbranch_down = number of downward transitions for each upward transition.
     line%nup = 2
     line%delE_Hz(1:2) = [83801.95_wp, 84004.26_wp] * speedc_cm
     line%delE_Hz(1:2) = line%delE_Hz(1) - line%delE_Hz(1:2)
     allocate(line%b(line%nup))
     line%b(1)%ndown = 2
     line%b(2)%ndown = 2
     do i=1, line%nup
        allocate(line%b(i)%A21(line%b(i)%ndown))
        allocate(line%b(i)%Elow_Hz(line%b(i)%ndown))
        allocate(line%b(i)%E1(line%b(i)%ndown))
        allocate(line%b(i)%E2(line%b(i)%ndown))
        allocate(line%b(i)%E3(line%b(i)%ndown))
     enddo

     line%b(1)%A21(1:2)     = [2.69e9_wp,   1.40e9_wp]
     line%b(1)%Elow_Hz(1:2) = [0.0_wp,    287.24_wp] * speedc_cm
     line%b(1)%E1(1:2)      = [0.0_wp,         0.0_wp]
     line%b(1)%E2(1:2)      = [1.0_wp,         1.0_wp]
     line%b(1)%E3(1:2)      = [2.0_wp/3.0_wp, -1.0_wp/3.0_wp]

     line%b(2)%A21(1:2)     = [6.53e8_wp,   3.45e9_wp]
     line%b(2)%Elow_Hz(1:2) = [0.0_wp,    287.24_wp] * speedc_cm
     line%b(2)%E1(1:2)      = [1.0_wp/2.0_wp, -2.0_wp/5.0_wp]
     line%b(2)%E2(1:2)      = [1.0_wp/2.0_wp,  7.0_wp/5.0_wp]
     line%b(2)%E3(1:2)      = [5.0_wp/6.0_wp, -1.0_wp/3.0_wp]

     do i=1, line%nup
        line%b(i)%damping = sum(line%b(i)%A21(:))
        if (line%b(i)%ndown > 1) then
           allocate(line%b(i)%P_down(line%b(i)%ndown))
           allocate(line%b(i)%A_down(line%b(i)%ndown))
           line%b(i)%P_down(:) = line%b(i)%A21(:) / line%b(i)%damping
           call random_alias_setup(line%b(i)%P_down, line%b(i)%A_down)
        endif
     enddo
     line%damping   = line%b(1)%damping
  case ('SiII_1304')
     line%ion_id    = 'Si II'
     line%line_type = 4
     line%wavelength0 = 0.1304370_wp
     ! the central frequency is assumed to be that with the highest f12.
     line%f12(1)    = 0.0928_wp
     line%cross0    = sigma_0/sqrt(pi)*line%f12(1)
     mass_amu       = 28.0855_wp
     line%vtherm1   = vtherm1_amu/sqrt(mass_amu)
     line%g_recoil0 = (h_planck/amu/mass_amu)/(line%wavelength0*um2m)**2

     !--- nbranch_up   = number of upward transitions
     !--- nbranch_down = number of downward transitions for each upward transition.
     line%nup = 1
     allocate(line%b(line%nup))
     line%b(1)%ndown = 2

     allocate(line%b(1)%A21(line%b(i)%ndown))
     allocate(line%b(1)%Elow_Hz(line%b(i)%ndown))
     allocate(line%b(1)%E1(line%b(i)%ndown))
     allocate(line%b(1)%E2(line%b(i)%ndown))
     allocate(line%b(1)%E3(line%b(i)%ndown))
     allocate(line%b(1)%P_down(line%b(1)%ndown))
     allocate(line%b(1)%A_down(line%b(1)%ndown))

     line%b(1)%A21(1:2)     = [3.64e8_wp,   6.23e8_wp]
     line%b(1)%Elow_Hz(1:2) = [0.0_wp,    287.24_wp] * speedc_cm
     line%b(1)%E1(1:2)      = [0.0_wp,         0.0_wp]
     line%b(1)%E2(1:2)      = [1.0_wp,         1.0_wp]
     line%b(1)%E3(1:2)      = [2.0_wp/3.0_wp, -1.0_wp/3.0_wp]
     line%b(1)%damping      = sum(line%b(1)%A21(:))
     line%b(1)%P_down(:)    = line%b(1)%A21(:)/line%b(1)%damping
     call random_alias_setup(line%b(1)%P_down, line%b(1)%A_down)
     line%damping   = line%b(1)%damping
  case ('FeII_2250')
     line%ion_id    = 'Fe II'
     line%line_type = 4
     line%wavelength0 = 0.224988_wp
     line%f12(1)    = 0.00182_wp

     line%cross0    = sigma_0/sqrt(pi)*line%f12(1)
     mass_amu       = 55.845_wp
     line%vtherm1   = vtherm1_amu/sqrt(mass_amu)
     line%g_recoil0 = (h_planck/amu/mass_amu)/(line%wavelength0*um2m)**2
     line%nup       = 1
     allocate(line%b(line%nup))

     line%b(1)%ndown = 2
     allocate(line%b(1)%A21(line%b(1)%ndown))
     allocate(line%b(1)%Elow_Hz(line%b(1)%ndown))
     allocate(line%b(1)%E1(line%b(1)%ndown))
     allocate(line%b(1)%E2(line%b(1)%ndown))
     allocate(line%b(1)%E3(line%b(1)%ndown))
     allocate(line%b(1)%P_down(line%b(1)%ndown))
     allocate(line%b(1)%A_down(line%b(1)%ndown))
     line%b(1)%A21(1:2)     = [3.00e6_wp,   4.00e5_wp]
     line%b(1)%Elow_Hz(1:2) = [0.0_wp,    384.7872_wp] * speedc_cm
     line%b(1)%E1(1:2)      = [  7.0_wp/150.0_wp, -2.0_wp/15.0_wp]
     line%b(1)%E2(1:2)      = [143.0_wp/150.0_wp, 17.0_wp/15.0_wp]
     line%b(1)%E3(1:2)      = [  7.0_wp/18.0_wp,  -1.0_wp/9.0_wp]
     line%b(1)%damping      = sum(line%b(1)%A21(:))
     line%b(1)%P_down(:)    = line%b(1)%A21(:)/line%b(1)%damping
     call random_alias_setup(line%b(1)%P_down, line%b(1)%A_down)
     line%damping   = line%b(1)%damping
  case ('FeII_2261')
     line%ion_id    = 'Fe II'
     line%line_type = 4
     line%wavelength0 = 0.226078_wp
     line%f12(1)    = 0.00244_wp

     line%cross0    = sigma_0/sqrt(pi)*line%f12(1)
     mass_amu       = 55.847_wp
     line%vtherm1   = vtherm1_amu/sqrt(mass_amu)
     line%g_recoil0 = (h_planck/amu/mass_amu)/(line%wavelength0*um2m)**2
     line%nup       = 1
     allocate(line%b(line%nup))

     line%b(1)%ndown = 2
     allocate(line%b(1)%A21(line%b(1)%ndown))
     allocate(line%b(1)%Elow_Hz(line%b(1)%ndown))
     allocate(line%b(1)%E1(line%b(1)%ndown))
     allocate(line%b(1)%E2(line%b(1)%ndown))
     allocate(line%b(1)%E3(line%b(1)%ndown))
     allocate(line%b(1)%P_down(line%b(1)%ndown))
     allocate(line%b(1)%A_down(line%b(1)%ndown))
     line%b(1)%A21(1:2)     = [3.18e6_wp,   4.49e6_wp]
     line%b(1)%Elow_Hz(1:2) = [0.0_wp,    384.7872_wp] * speedc_cm
     line%b(1)%E1(1:2)      = [ 64.0_wp/165.0_wp, -4.0_wp/15.0_wp]
     line%b(1)%E2(1:2)      = [101.0_wp/165.0_wp, 19.0_wp/15.0_wp]
     line%b(1)%E3(1:2)      = [  2.0_wp/99.0_wp,   1.0_wp/9.0_wp]
     line%b(1)%damping      = sum(line%b(1)%A21(:))
     line%b(1)%P_down(:)    = line%b(1)%A21(:)/line%b(1)%damping
     call random_alias_setup(line%b(1)%P_down, line%b(1)%A_down)
     line%damping   = line%b(1)%damping
  case ('FeII_UV3', 'FeII_2344')
     !-- vacuum wavelengths 2344.21274, 2365.55055, 2381.48756 (resonance + fluorescence)
     !-- Aki = 1.73e8, fik = 0.114 / Aki = 5.9e7, fik = 0.0495 / Aki = 3.1e7, fik = 0.0351
     line%ion_id    = 'Fe II'
     line%line_type = 4
     line%wavelength0 = 0.234421274_wp
     line%f12(1)    = 0.114_wp
     !line%damping   = 1.73e8_wp + 5.90e7_wp + 3.10e7_wp

     line%cross0    = sigma_0/sqrt(pi)*line%f12(1)
     mass_amu       = 55.847_wp
     line%vtherm1   = vtherm1_amu/sqrt(mass_amu)
     line%g_recoil0 = (h_planck/amu/mass_amu)/(line%wavelength0*um2m)**2

     line%nup       = 1
     allocate(line%b(line%nup))

     line%b(1)%ndown = 3
     allocate(line%b(1)%A21(line%b(1)%ndown))
     allocate(line%b(1)%Elow_Hz(line%b(1)%ndown))
     allocate(line%b(1)%E1(line%b(1)%ndown))
     allocate(line%b(1)%E2(line%b(1)%ndown))
     allocate(line%b(1)%E3(line%b(1)%ndown))
     allocate(line%b(1)%P_down(line%b(1)%ndown))
     allocate(line%b(1)%A_down(line%b(1)%ndown))

     line%b(1)%A21(1:3)     = [1.73e8_wp,   5.90e7_wp,   3.10e7_wp]
     line%b(1)%Elow_Hz(1:3) = [0.0_wp,    384.7872_wp, 667.6829_wp] * speedc_cm
     line%b(1)%E1(1:3)      = [7.0_wp/150.0_wp,  -2.0_wp/15.0_wp,  1.0_wp/10.0_wp]
     line%b(1)%E2(1:3)      = [143_wp/150.0_wp,  17.0_wp/15.0_wp,  9.0_wp/10.0_wp]
     line%b(1)%E3(1:3)      = [7.0_wp/18.0_wp,   -1.0_wp/9.0_wp,  -1.0_wp/2.0_wp]
     line%b(1)%damping      = sum(line%b(1)%A21(:))
     line%b(1)%P_down(:)    = line%b(1)%A21(:)/line%b(1)%damping
     call random_alias_setup(line%b(1)%P_down, line%b(1)%A_down)
     line%damping   = line%b(1)%damping
  case ('FeII_UV1', 'FeII_2600')
     line%ion_id    = 'Fe II'
     line%line_type = 5
     line%wavelength0 = 0.260017206_wp
     ! the central frequency is assumed to be that with the highest f12.
     line%f12(1:2)  = [0.239_wp, 0.0717_wp]
     line%cross0    = sigma_0/sqrt(pi)*line%f12(1)
     mass_amu       = 55.847_wp
     line%vtherm1   = vtherm1_amu/sqrt(mass_amu)
     line%g_recoil0 = (h_planck/amu/mass_amu)/(line%wavelength0*um2m)**2

     !--- nbranch_up   = number of upward transitions
     !--- nbranch_down = number of downward transitions for each upward transition.
     line%nup = 2
     line%delE_Hz(1:2) = [38458.9934_wp, 38660.0537_wp] * speedc_cm
     line%delE_Hz(1:2) = line%delE_Hz(1) - line%delE_Hz(1:2)
     allocate(line%b(line%nup))
     line%b(1)%ndown = 2
     line%b(2)%ndown = 3
     do i=1, line%nup
        allocate(line%b(i)%A21(line%b(i)%ndown))
        allocate(line%b(i)%Elow_Hz(line%b(i)%ndown))
        allocate(line%b(i)%E1(line%b(i)%ndown))
        allocate(line%b(i)%E2(line%b(i)%ndown))
        allocate(line%b(i)%E3(line%b(i)%ndown))
     enddo
#ifdef PROCHASKA
     line%b(1)%A21(1:2)     = [2.36e8_wp,   3.41e7_wp]
     line%b(1)%Elow_Hz(1:2) = [0.0_wp,    384.7872_wp] * speedc_cm
     line%b(1)%E1(1:2)      = [64.0_wp/165.0_wp,  -4.0_wp/15.0_wp]
     line%b(1)%E2(1:2)      = [101.0_wp/165.0_wp, 19.0_wp/15.0_wp]
     line%b(1)%E3(1:2)      = [2.0_wp/99.0_wp,     1.0_wp/9.0_wp]

     line%b(2)%A21(1:3)     = [8.61e7_wp,   1.23e8_wp,    6.21e7_wp]
     line%b(2)%Elow_Hz(1:3) = [0.0_wp,    384.7872_wp,  667.6829_wp] * speedc_cm
     line%b(2)%E1(1:3)      = [0.0_wp, 0.0_wp,  0.0_wp]
     line%b(2)%E2(1:3)      = [1.0_wp, 1.0_wp,  1.0_wp]
     line%b(2)%E3(1:3)      = [0.0_wp, 0.0_wp,  0.0_wp]
#else
     line%b(1)%A21(1:2)     = [2.35e8_wp,   3.52e8_wp]
     line%b(1)%Elow_Hz(1:2) = [0.0_wp,    384.7872_wp] * speedc_cm
     line%b(1)%E1(1:2)      = [64.0_wp/165.0_wp,  -4.0_wp/15.0_wp]
     line%b(1)%E2(1:2)      = [101.0_wp/165.0_wp, 19.0_wp/15.0_wp]
     line%b(1)%E3(1:2)      = [2.0_wp/99.0_wp,     1.0_wp/9.0_wp]

     line%b(2)%A21(1:3)     = [8.94e7_wp,   1.20e8_wp,    6.29e7_wp]
     line%b(2)%Elow_Hz(1:3) = [0.0_wp,    384.7872_wp,  667.6829_wp] * speedc_cm
     line%b(2)%E1(1:3)      = [7.0_wp/150.0_wp,   -2.0_wp/15.0_wp,  1.0_wp/10.0_wp]
     line%b(2)%E2(1:3)      = [143.0_wp/150.0_wp, 17.0_wp/15.0_wp,  9.0_wp/10.0_wp]
     line%b(2)%E3(1:3)      = [7.0_wp/18.0_wp,    -1.0_wp/9.0_wp,  -1.0_wp/2.0_wp]
#endif
     do i=1, line%nup
        line%b(i)%damping = sum(line%b(i)%A21(:))
        if (line%b(i)%ndown > 1) then
           allocate(line%b(i)%P_down(line%b(i)%ndown))
           allocate(line%b(i)%A_down(line%b(i)%ndown))
           line%b(i)%P_down(:) = line%b(i)%A21(:) / line%b(i)%damping
           call random_alias_setup(line%b(i)%P_down, line%b(i)%A_down)
        endif
     enddo
     line%damping   = line%b(1)%damping
  case ('FeII_UV2', 'FeII_2383')
     line%ion_id    = 'Fe II'
     line%line_type = 5
     line%wavelength0 = 0.238276386_wp
     ! the central frequency is assumed to be that with the highest f12.
     line%f12(1:2)  = [0.320_wp, 0.0359_wp]
     line%cross0    = sigma_0/sqrt(pi)*line%f12(1)
     mass_amu       = 55.847_wp
     line%vtherm1   = vtherm1_amu/sqrt(mass_amu)
     line%g_recoil0 = (h_planck/amu/mass_amu)/(line%wavelength0*um2m)**2

     !--- nbranch_up   = number of upward transitions
     !--- nbranch_down = number of downward transitions for each upward transition.
     line%nup = 2
     line%delE_Hz(1:2) = [41968.0698_wp, 42114.8380_wp] * speedc_cm
     line%delE_Hz(1:2) = line%delE_Hz(1) - line%delE_Hz(1:2)
     allocate(line%b(line%nup))
     line%b(1)%ndown = 1
     line%b(2)%ndown = 2
     do i=1, line%nup
        allocate(line%b(i)%A21(line%b(i)%ndown))
        allocate(line%b(i)%Elow_Hz(line%b(i)%ndown))
        allocate(line%b(i)%E1(line%b(i)%ndown))
        allocate(line%b(i)%E2(line%b(i)%ndown))
        allocate(line%b(i)%E3(line%b(i)%ndown))
     enddo
     line%b(1)%A21(1:1)     = [3.13e8_wp]
     line%b(1)%Elow_Hz(1:1) = [0.0_wp] * speedc_cm
     line%b(1)%E1(1:1)      = [91.0_wp/550.0_wp]
     line%b(1)%E2(1:1)      = [459.0_wp/550.0_wp]
     line%b(1)%E3(1:1)      = [12.0_wp/22.0_wp]

     line%b(2)%A21(1:2)     = [4.25e7_wp,   2.59e8_wp]
     line%b(2)%Elow_Hz(1:2) = [0.0_wp,    384.7872_wp] * speedc_cm
     line%b(2)%E1(1:2)      = [64.0_wp/165.0_wp,   -4.0_wp/15.0_wp]
     line%b(2)%E2(1:2)      = [101.0_wp/165.0_wp,  19.0_wp/15.0_wp]
     line%b(2)%E3(1:2)      = [2.0_wp/99.0_wp,      1.0_wp/9.0_wp]
     do i=1, line%nup
        line%b(i)%damping = sum(line%b(i)%A21(:))
        if (line%b(i)%ndown > 1) then
           allocate(line%b(i)%P_down(line%b(i)%ndown))
           allocate(line%b(i)%A_down(line%b(i)%ndown))
           line%b(i)%P_down(:) = line%b(i)%A21(:) / line%b(i)%damping
           call random_alias_setup(line%b(i)%P_down, line%b(i)%A_down)
        endif
     enddo
     line%damping   = line%b(1)%damping
  case ('HeI_10833')
     line%ion_id    = 'He I'
     line%line_type = 6
     line%wavelength0 = 1.0833306444
     ! the central frequency is assumed to be that with the highest f12.
     line%f12(1:3)  = [2.9958e-1, 1.797e-1, 5.9902e-2]
     line%cross0    = sigma_0/sqrt(pi)*line%f12(1)
     mass_amu       = 4.0026032545_wp
     line%vtherm1   = vtherm1_amu/sqrt(mass_amu)
     line%g_recoil0 = (h_planck/amu/mass_amu)/(line%wavelength0*um2m)**2

     !--- nbranch_up   = number of upward transitions
     !--- nbranch_down = number of downward transitions for each upward transition.
     line%nup = 3
     line%delE_Hz(1:3) = [169086.7664725_wp, 169086.8428979_wp, 169087.8308131_wp] * speedc_cm
     line%delE_Hz(1:3) = line%delE_Hz(1) - line%delE_Hz(1:3)
     allocate(line%b(line%nup))
     line%b(1)%ndown = 1
     line%b(2)%ndown = 1
     line%b(3)%ndown = 1
     do i=1, line%nup
        allocate(line%b(i)%A21(line%b(i)%ndown))
        allocate(line%b(i)%Elow_Hz(line%b(i)%ndown))
        allocate(line%b(i)%E1(line%b(i)%ndown))
        allocate(line%b(i)%E2(line%b(i)%ndown))
        allocate(line%b(i)%E3(line%b(i)%ndown))
     enddo
     line%b(1)%A21(1:1)     = [1.0216e7_wp]
     line%b(1)%Elow_Hz(1:1) = [0.0_wp] * speedc_cm
     line%b(1)%E1(1:1)      = [7.0_wp/20.0_wp]
     line%b(1)%E2(1:1)      = [13.0_wp/20.0_wp]
     line%b(1)%E3(1:1)      = [3.0_wp/4.0_wp]

     line%b(2)%A21(1:1)     = [1.0216e7_wp]
     line%b(2)%Elow_Hz(1:1) = [0.0_wp] * speedc_cm
     line%b(2)%E1(1:1)      = [1.0_wp/4.0_wp]
     line%b(2)%E2(1:1)      = [3.0_wp/4.0_wp]
     line%b(2)%E3(1:1)      = [1.0_wp/4.0_wp]

     line%b(3)%A21(1:1)     = [1.0216e7_wp]
     line%b(3)%Elow_Hz(1:1) = [0.0_wp] * speedc_cm
     line%b(3)%E1(1:1)      = [0.0_wp]
     line%b(3)%E2(1:1)      = [1.0_wp]
     line%b(3)%E3(1:1)      = [0.0_wp]

     do i=1, line%nup
        line%b(i)%damping = sum(line%b(i)%A21(:))
     enddo
     line%damping   = line%b(1)%damping
  case ('ly_alpha_HD')
     !-- Hydrogen + Deuterium Lyman-α (line_type = 7).
     !-- Two coexisting two-level resonance scatterers. Fine-structure
     !-- splitting is ignored for both species in this version. The deuterium
     !-- abundance defaults to par%D_to_H_ratio (cosmic 1.5e-5).
     line%ion_id    = 'H+D'
     line%line_type = 7
     !-- Hydrogen primary species (same as 'ly_alpha' non-fine-structure branch)
     line%wavelength0 = 0.1215668237310_wp
     line%damping     = 6.2649e8_wp
     line%f12(1:2)    = [0.27760_wp, 0.13881_wp]
     line%cross0      = sigma_0/sqrt(pi)*sum(line%f12)
     mass_amu         = 1.00797_wp
     line%mass_amu    = mass_amu
     line%vtherm1     = vtherm1_amu/sqrt(mass_amu)
     line%g_recoil0   = (h_planck/amu/mass_amu)/(line%wavelength0*um2m)**2
     line%DnuHK_Hz    = 0.0_wp
     line%E1          = 1.0_wp
     line%E2          = 0.0_wp
     line%E3          = 1.0_wp
     line%nup         = 1
     line%delE_Hz     = 0.0_wp
     !-- Deuterium secondary species (NIST vacuum λ; Γ ≈ Γ_H to <0.1%)
     line%wavelength0_D = 0.1215337431_wp
     line%mass_amu_D    = 2.01410177812_wp
     line%damping_D     = 6.2649e8_wp
     line%f12_D         = sum(line%f12)
     line%cross0_D      = sigma_0/sqrt(pi)*line%f12_D
     line%vtherm1_D     = vtherm1_amu/sqrt(line%mass_amu_D)
     line%g_recoil0_D   = (h_planck/amu/line%mass_amu_D)/(line%wavelength0_D*um2m)**2
     !-- Cross-species precomputed constants
     !-- delta_nu_HD_Hz = c [cm/s] * (1/λ_D − 1/λ_H) [1/cm], λ in cm.
     line%nD_over_nH      = par%D_to_H_ratio
     line%delta_nu_HD_Hz  = speedc_cm * (1.0_wp/(line%wavelength0_D*um2m*1.0e2_wp) - &
                                          1.0_wp/(line%wavelength0  *um2m*1.0e2_wp))
     line%ratio_Dfreq_HD  = (line%wavelength0_D/line%wavelength0) * sqrt(line%mass_amu_D/mass_amu)
     line%ratio_voigta_HD = (line%damping_D/line%damping) * line%ratio_Dfreq_HD
  case default
     !-- Lyman-alpha
     !-- line%wavelength0 is in units of micron.
     !-- vacuum wavelengths 1215.668237310 (short, K, S1/2-P3/2), 1215.673644608 (long, H, S1/2-P1/2)
     !-- observed wavelengths 1215.6699 (short), 1215.6699 (long)
     !-- Aki = 6.2648e8 (short), 6.2649e8 (long), fik = 2.7760e-1 (short), 1.3881e-1 (long)
     line%ion_id    = 'H  I'
     line%delE_Hz(1:2) = [82259.2850014_wp, 82258.9191133_wp] * speedc_cm
     line%delE_Hz(1:2) = line%delE_Hz(1) - line%delE_Hz(1:2)
     if (par%fine_structure) then
        line%line_type = 2
        !line%DnuHK_Hz  = speedc_cm * (82259.2850014_wp - 82258.9191133_wp)
        line%DnuHK_Hz  = line%delE_Hz(2)
     else
        line%line_type = 1
        line%DnuHK_Hz  = 0.0_wp
        line%E1        = 1.0_wp
        line%E2        = 0.0_wp
        line%E3        = 1.0_wp
     endif
     line%wavelength0 = 0.1215668237310_wp
     !line%f12(1)    = 0.41641_wp
     !line%cross0    = sigma_0/sqrt(pi)*line%f12(1)
     !line%damping   = 6.26483e8_wp
     line%damping   = 6.2649e8_wp
     line%f12(1:2)  = [0.27760_wp, 0.13881_wp]
     line%cross0    = sigma_0/sqrt(pi)*sum(line%f12)
     mass_amu       = 1.00797_wp
     line%vtherm1   = vtherm1_amu/sqrt(mass_amu)
     line%g_recoil0 = (h_planck/amu/mass_amu)/(line%wavelength0*um2m)**2
  endselect
  end subroutine setup_resonance_line
!++++++++++++++++++++++++++++++++++++++
end module line_mod
