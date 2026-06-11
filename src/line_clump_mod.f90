module line_clump_mod
!---------------------------------------------------------------------------
! Clump-mode resonance variants do_resonance{1,2,4,5,6,HD}_clump.
!
! Used when par%use_clump_medium = .true. and the photon is inside a clump
! (photon%icell_clump > 0).  They read the per-clump Voigt parameter
! cl_voigt_a(icl) and Doppler frequency cl_Dfreq(icl) instead of the (uniform)
! Cartesian grid arrays, and rescale photon%xfreq to the local clump
! Doppler-width units for sampling, then convert uz and xfreq_atom back to
! REF Doppler units (cl_Dfreq_ref) before returning, so the calling
! scatter_resonance_* code keeps using a single global xfreq convention.
!
! For uniform-T clumps (cl_Dfreq(icl) == cl_Dfreq_ref) these are numerically
! equivalent to the standard do_resonance{N} in line_mod.
!
! This module is separate from line_mod because clump_mod uses line_mod
! internally (read/write_clumps_info), and putting "use clump_mod" inside
! line_mod would create a module-level cycle.  line_clump_mod is compiled
! after both line_mod and clump_mod.
!---------------------------------------------------------------------------
  use define
  use voigt_mod
  use random
  use line_mod,  only: line, compute_HeI_E_coherent
  use clump_mod, only: cl_voigt_a, cl_Dfreq, cl_Dfreq_ref
  implicit none
  public
contains
!++++++++++++++++++++++++++++++++++++++
  subroutine do_resonance1_clump(photon,grid,uz,xfreq_atom,cost,sint,S11,S22,S12,S33,S44)
  implicit none
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(in)    :: grid
  real(kind=wp),     intent(out)   :: uz, xfreq_atom, cost, sint
  real(kind=wp), optional, intent(out) :: S11, S22, S12, S33, S44
  real(kind=wp) :: cost2, voigt_a_loc, xfreq_loc, uz_loc, scale, scale_inv
  integer(int64) :: icl
  icl         = int(photon%icell_clump, int64)
  scale       = cl_Dfreq_ref / cl_Dfreq(icl)
  scale_inv   = 1.0_wp / scale
  voigt_a_loc = cl_voigt_a(icl)
  xfreq_loc   = photon%xfreq * scale
  uz_loc      = rand_resonance_vz(xfreq_loc, voigt_a_loc)
  xfreq_atom  = (xfreq_loc - uz_loc) * scale_inv
  uz          = uz_loc * scale_inv
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
  end subroutine do_resonance1_clump
!++++++++++++++++++++++++++++++++++++++
  subroutine do_resonance2_clump(photon,grid,uz,xfreq_atom,cost,sint,S11,S22,S12,S33,S44)
  implicit none
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(in)    :: grid
  real(kind=wp),     intent(out)   :: uz, xfreq_atom, cost, sint
  real(kind=wp), optional, intent(out) :: S11, S22, S12, S33, S44
  real(kind=wp) :: DnuHK, pK, pH, qH, qK, E1, cost2
  real(kind=wp) :: voigt_a_loc, xfreq_loc, uz_loc, xfreq_atom_loc, scale, scale_inv
  integer(int64) :: icl
  icl         = int(photon%icell_clump, int64)
  scale       = cl_Dfreq_ref / cl_Dfreq(icl)
  scale_inv   = 1.0_wp / scale
  voigt_a_loc = cl_voigt_a(icl)
  xfreq_loc   = photon%xfreq * scale
  DnuHK       = line%DnuHK_Hz / cl_Dfreq(icl)
  pH          = voigt(xfreq_loc + DnuHK, voigt_a_loc) * one_over_three
  pK          = voigt(xfreq_loc,         voigt_a_loc) * two_over_three
  pH          = pH/(pH + pK)
  if (rand_number() < pH) then
     uz_loc = rand_resonance_vz(xfreq_loc + DnuHK, voigt_a_loc)
  else
     uz_loc = rand_resonance_vz(xfreq_loc,         voigt_a_loc)
  endif
  xfreq_atom_loc = xfreq_loc - uz_loc
  qH    = xfreq_atom_loc + DnuHK
  qK    = xfreq_atom_loc
  E1    = (2.0_wp*qK*qH + qH**2)/(qK**2 + 2.0_wp*qH**2)
  cost  = rand_resonance(E1)
  cost2 = cost**2
  sint  = sqrt(1.0_wp - cost2)
  if (present(S44)) then
     S22 = three_over_four * E1 * (cost2 + 1.0_wp)
     S11 = S22 + (1.0_wp - E1)
     S12 = three_over_four * E1 * (cost2 - 1.0_wp)
     S33 = three_over_two  * E1 * cost
     S44 = (E1 + 2.0_wp)/2.0_wp * cost
  endif
  photon%E1 = E1
  photon%E2 = 1.0_wp - E1
  photon%E3 = (E1 + 2.0_wp)/3.0_wp
  xfreq_atom = xfreq_atom_loc * scale_inv
  uz         = uz_loc * scale_inv
  end subroutine do_resonance2_clump
!++++++++++++++++++++++++++++++++++++++
  subroutine do_resonance4_clump(photon,grid,uz,xfreq_atom,cost,sint,S11,S22,S12,S33,S44)
  implicit none
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(in)    :: grid
  real(kind=wp),     intent(out)   :: uz, xfreq_atom, cost, sint
  real(kind=wp), optional, intent(out) :: S11, S22, S12, S33, S44
  real(kind=wp) :: cost2, del_xfreq
  real(kind=wp) :: voigt_a_loc, xfreq_loc, uz_loc, xfreq_atom_loc, scale, scale_inv
  integer :: idown
  integer(int64) :: icl
  icl         = int(photon%icell_clump, int64)
  scale       = cl_Dfreq_ref / cl_Dfreq(icl)
  scale_inv   = 1.0_wp / scale
  voigt_a_loc = cl_voigt_a(icl)
  xfreq_loc   = photon%xfreq * scale
  uz_loc      = rand_resonance_vz(xfreq_loc, voigt_a_loc)
  xfreq_atom_loc = xfreq_loc - uz_loc
  idown = rand_alias_choise(line%b(1)%P_down, line%b(1)%A_down)
  if (idown /= 1) then
     del_xfreq      = line%b(1)%Elow_Hz(idown) / cl_Dfreq(icl)
     xfreq_atom_loc = xfreq_atom_loc - del_xfreq
  endif
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
  endif
  xfreq_atom = xfreq_atom_loc * scale_inv
  uz         = uz_loc * scale_inv
  end subroutine do_resonance4_clump
!++++++++++++++++++++++++++++++++++++++
  subroutine do_resonance5_clump(photon,grid,uz,xfreq_atom,cost,sint,S11,S22,S12,S33,S44)
  implicit none
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(in)    :: grid
  real(kind=wp),     intent(out)   :: uz, xfreq_atom, cost, sint
  real(kind=wp), optional, intent(out) :: S11, S22, S12, S33, S44
  real(kind=wp) :: Dx, p1, p2, va1, va2, cost2, E1, del_xfreq
  real(kind=wp) :: voigt_a_loc, xfreq_loc, uz_loc, xfreq_atom_loc, scale, scale_inv
  integer :: iup, idown
  integer(int64) :: icl
  icl         = int(photon%icell_clump, int64)
  scale       = cl_Dfreq_ref / cl_Dfreq(icl)
  scale_inv   = 1.0_wp / scale
  voigt_a_loc = cl_voigt_a(icl)
  xfreq_loc   = photon%xfreq * scale
  Dx  = line%delE_Hz(2) / cl_Dfreq(icl)
  va1 = voigt_a_loc
  va2 = voigt_a_loc * line%b(2)%damping/line%b(1)%damping
  p1  = voigt(xfreq_loc,      va1) * line%f12(1)
  p2  = voigt(xfreq_loc + Dx, va2) * line%f12(2)
  p1  = p1/(p2 + p1)
  if (rand_number() < p1) then
     uz_loc = rand_resonance_vz(xfreq_loc,      va1)
     iup    = 1
  else
     uz_loc = rand_resonance_vz(xfreq_loc + Dx, va2)
     iup    = 2
  endif
  xfreq_atom_loc = xfreq_loc - uz_loc
  if (line%b(iup)%ndown > 1) then
     idown          = rand_alias_choise(line%b(iup)%P_down, line%b(iup)%A_down)
     del_xfreq      = line%b(iup)%Elow_Hz(idown) / cl_Dfreq(icl)
     xfreq_atom_loc = xfreq_atom_loc - del_xfreq
  else
     idown = 1
  endif
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
  endif
  xfreq_atom = xfreq_atom_loc * scale_inv
  uz         = uz_loc * scale_inv
  end subroutine do_resonance5_clump
!++++++++++++++++++++++++++++++++++++++
  subroutine do_resonance6_clump(photon,grid,uz,xfreq_atom,cost,sint,S11,S22,S12,S33,S44)
  implicit none
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(in)    :: grid
  real(kind=wp),     intent(out)   :: uz, xfreq_atom, cost, sint
  real(kind=wp), optional, intent(out) :: S11, S22, S12, S33, S44
  real(kind=wp) :: Dx2, Dx3, p1, p2, p3, ptot, pcum1, pcum2, xi
  real(kind=wp) :: va1, va2, va3, cost2, E1, del_xfreq
  real(kind=wp) :: voigt_a_loc, xfreq_loc, uz_loc, xfreq_atom_loc, scale, scale_inv
  integer :: iup, idown
  integer(int64) :: icl
  icl         = int(photon%icell_clump, int64)
  scale       = cl_Dfreq_ref / cl_Dfreq(icl)
  scale_inv   = 1.0_wp / scale
  voigt_a_loc = cl_voigt_a(icl)
  xfreq_loc   = photon%xfreq * scale
  Dx2 = line%delE_Hz(2) / cl_Dfreq(icl)
  Dx3 = line%delE_Hz(3) / cl_Dfreq(icl)
  va1 = voigt_a_loc
  va2 = voigt_a_loc * line%b(2)%damping/line%b(1)%damping
  va3 = voigt_a_loc * line%b(3)%damping/line%b(1)%damping
  p1  = voigt(xfreq_loc,       va1) * line%f12(1)
  p2  = voigt(xfreq_loc + Dx2, va2) * line%f12(2)
  p3  = voigt(xfreq_loc + Dx3, va3) * line%f12(3)
  ptot  = p1 + p2 + p3
  pcum1 = p1/ptot
  pcum2 = (p1 + p2)/ptot
  xi = rand_number()
  if (xi < pcum1) then
     uz_loc = rand_resonance_vz(xfreq_loc,       va1)
     iup    = 1
  else if (xi < pcum2) then
     uz_loc = rand_resonance_vz(xfreq_loc + Dx2, va2)
     iup    = 2
  else
     uz_loc = rand_resonance_vz(xfreq_loc + Dx3, va3)
     iup    = 3
  endif
  xfreq_atom_loc = xfreq_loc - uz_loc
  if (par%HeI_coherent) then
     call compute_HeI_E_coherent(xfreq_atom_loc, Dx2, Dx3, &
                                 photon%E1, photon%E2, photon%E3)
  else
     photon%E1 = line%b(iup)%E1(1)
     photon%E2 = line%b(iup)%E2(1)
     photon%E3 = line%b(iup)%E3(1)
  endif
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
  xfreq_atom = xfreq_atom_loc * scale_inv
  uz         = uz_loc * scale_inv
  end subroutine do_resonance6_clump
!++++++++++++++++++++++++++++++++++++++
  subroutine do_resonance_HD_clump(photon,grid,uz,xfreq_atom,cost,sint,S11,S22,S12,S33,S44)
  implicit none
  type(photon_type), intent(inout) :: photon
  type(grid_type),   intent(in)    :: grid
  real(kind=wp),     intent(out)   :: uz, xfreq_atom, cost, sint
  real(kind=wp), optional, intent(out) :: S11, S22, S12, S33, S44
  real(kind=wp) :: voigt_a_H, dx_HD, xfreq_D_loc
  real(kind=wp) :: pH, pD_term, pH_norm, uz_D, cost2
  real(kind=wp) :: xfreq_loc, uz_loc, xfreq_atom_loc, scale, scale_inv
  integer(int64) :: icl
  icl         = int(photon%icell_clump, int64)
  scale       = cl_Dfreq_ref / cl_Dfreq(icl)
  scale_inv   = 1.0_wp / scale
  voigt_a_H   = cl_voigt_a(icl)
  xfreq_loc   = photon%xfreq * scale
  dx_HD       = line%delta_nu_HD_Hz / cl_Dfreq(icl)
  xfreq_D_loc = (xfreq_loc - dx_HD) * line%ratio_Dfreq_HD
  pH      = voigt(xfreq_loc, voigt_a_H)
  pD_term = line%nD_over_nH * line%ratio_Dfreq_HD * &
            voigt(xfreq_D_loc, voigt_a_H * line%ratio_voigta_HD)
  pH_norm = pH / (pH + pD_term)
  if (rand_number() < pH_norm) then
     line%selected_species_HD = 1
     uz_loc         = rand_resonance_vz(xfreq_loc, voigt_a_H)
     xfreq_atom_loc = xfreq_loc - uz_loc
  else
     line%selected_species_HD = 2
     uz_D           = rand_resonance_vz(xfreq_D_loc, voigt_a_H * line%ratio_voigta_HD)
     uz_loc         = uz_D / line%ratio_Dfreq_HD
     xfreq_atom_loc = xfreq_loc - uz_loc
  endif
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
  xfreq_atom = xfreq_atom_loc * scale_inv
  uz         = uz_loc * scale_inv
  end subroutine do_resonance_HD_clump
!++++++++++++++++++++++++++++++++++++++
end module line_clump_mod
