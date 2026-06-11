module physics_amr_mod
  !-----------------------------------------------------------------------
  ! Physics utility functions for AMR grid initialization.
  !
  ! Shared by both LaRT.x (grid_mod_amr.f90) and the converter
  ! (convert_ramses_to_generic.f90), ensuring identical formulas.
  !
  ! Functions:
  !   cie_neutral_fraction_formula(T) -- single CIE formula (existing)
  !   laursen09_dust(nH, xHI, Z, Z_ref, f_ion) -- dust from metallicity
  !   caseB_emissivity(nH, T, xHI, ne) -- Lya emissivity
  !   solar_ion_density(nH, Z, T, ion_id) -- metal ion density
  !-----------------------------------------------------------------------
  use define, only: wp
  implicit none
  private

  public :: cie_neutral_fraction_formula
  public :: cie_neutral_fraction_table
  public :: laursen09_ndust
  public :: caseB_lya_emissivity
  public :: electron_density_from_xHI

  ! Physical constants
  real(wp), parameter :: massH_cgs     = 1.6726e-24_wp   ! hydrogen mass [g]
  real(wp), parameter :: boltzmann_cgs = 1.381e-16_wp    ! Boltzmann constant [erg/K]

contains

  !=========================================================================
  ! CIE neutral fraction from single-formula ionization/recombination rates.
  ! This is the existing formula used in grid_mod_amr.f90.
  !=========================================================================
  elemental function cie_neutral_fraction_formula(T) result(xHI)
    real(wp), intent(in) :: T
    real(wp) :: xHI
    real(wp) :: T4, k_ionize, k_recomb

    T4       = max(T, 10.0_wp) / 1.0e4_wp
    k_ionize = 5.84862e-9_wp * sqrt(T4) * exp(-15.78215_wp / T4)
    k_recomb = 4.13e-13_wp * T4**(-0.7131_wp - 0.0115_wp * log(T4))
    xHI      = k_recomb / (k_ionize + k_recomb)
  end function cie_neutral_fraction_formula

  !=========================================================================
  ! Electron density from neutral fraction (assuming primordial composition).
  ! ne ~ nH * (1 - xHI) for hydrogen-only; He contributes a small correction.
  !=========================================================================
  elemental function electron_density_from_xHI(nH, xHI) result(ne)
    real(wp), intent(in) :: nH, xHI
    real(wp) :: ne
    ne = nH * (1.0_wp - xHI)
  end function electron_density_from_xHI

  !=========================================================================
  ! Dust pseudo-number density from metallicity (Laursen+09).
  !
  ! ndust = (Z / Z_ref) * (nH * xHI + f_ion * nH * (1-xHI))
  !
  ! The returned ndust is multiplied by cext_dust * distance2cm in the
  ! caller to get rhokapD.
  !=========================================================================
  elemental function laursen09_ndust(nH, xHI, Z, Z_ref, f_ion) result(ndust)
    real(wp), intent(in) :: nH, xHI, Z, Z_ref, f_ion
    real(wp) :: ndust
    real(wp) :: nHI, nHII

    nHI   = nH * xHI
    nHII  = nH * (1.0_wp - xHI)
    ndust = (Z / max(Z_ref, 1.0e-30_wp)) * (nHI + f_ion * nHII)
  end function laursen09_ndust

  !=========================================================================
  ! Case B Lyman-alpha emissivity [cm^-3 s^-1].
  !
  ! emissivity = P_B * alpha_B * ne * nHII  (recombination)
  !            + q_coll * ne * nHI           (collisional excitation)
  !
  ! References:
  !   alpha_B : Hui & Gnedin (1997), MNRAS 292, 27
  !   P_B     : Cantalupo et al. (2008), ApJ 672, 48
  !   q_coll  : Harley/RASCAS formula
  !=========================================================================
  elemental function caseB_lya_emissivity(nH, T, xHI, ne) result(emiss)
    real(wp), intent(in) :: nH, T, xHI, ne
    real(wp) :: emiss
    real(wp) :: lambda_HG, alpha_B, P_B, Ta
    real(wp) :: nHI, nHII, emiss_recomb, emiss_coll, q_coll

    ! Case B recombination coefficient (Hui & Gnedin 1997)
    lambda_HG = 315614.0_wp / max(T, 10.0_wp)
    alpha_B   = 2.753e-14_wp * lambda_HG**1.5_wp &
              / (1.0_wp + (lambda_HG / 2.74_wp)**0.407_wp)**2.242_wp

    ! Case B Lya fraction (Cantalupo+ 2008)
    Ta  = max(T, 100.0_wp)
    P_B = 0.686_wp - 0.106_wp * log10(Ta / 1.0e4_wp) &
        - 0.009_wp * (Ta / 1.0e4_wp)**(-0.44_wp)

    nHI  = nH * xHI
    nHII = nH * (1.0_wp - xHI)
    emiss_recomb = P_B * alpha_B * ne * nHII

    ! Collisional excitation rate (Harley/RASCAS formula)
    q_coll = (6.58e-18_wp / max(T, 10.0_wp)**0.185_wp) &
           * exp(-4.86e4_wp / max(T, 10.0_wp)**0.895_wp)
    emiss_coll = nHI * ne * q_coll

    emiss = emiss_recomb + emiss_coll
  end function caseB_lya_emissivity

  !=========================================================================
  ! CIE neutral fraction from pre-computed table with log-linear
  ! interpolation. Uses Voronov (1997) collisional ionization rate and
  ! Verner & Ferland (1996) Case A recombination rate.
  !
  !   Gamma_HI(T)  = 5.85e-11 * T^{1/2} * exp(-157809.1/T)
  !                  / (1 + (T/1e5)^{1/2})                  [cm^3/s]
  !   alpha_A(T)   = 4.309e-13 * (T/1e4)^{-0.6166}
  !                  / (1 + 0.6703*(T/1e4)^{0.5300})        [cm^3/s]
  !   xHI          = alpha_A / (Gamma_HI + alpha_A)
  !
  ! Table: log10(T) from 3.0 to 8.0 in steps of 0.1 (51 points).
  ! Interpolation is linear in log10(T) and log10(xHI).
  !
  ! nH is unused (CIE is density-independent) but included in the
  ! interface for future density-dependent extensions.
  !=========================================================================
  function cie_neutral_fraction_table(nH, T) result(xHI)
    real(wp), intent(in) :: nH, T
    real(wp) :: xHI

    integer, parameter :: NTAB = 51
    real(wp), parameter :: LOGT_MIN = 3.0_wp
    real(wp), parameter :: LOGT_MAX = 8.0_wp
    real(wp), parameter :: DLOGT    = 0.1_wp

    real(wp), save :: log_xHI_tab(NTAB)
    logical, save  :: initialized = .false.

    real(wp) :: logT, frac, log_xHI_interp
    integer  :: idx

    ! --- Build the table on first call ---
    if (.not. initialized) then
      call build_cie_table(NTAB, LOGT_MIN, DLOGT, log_xHI_tab)
      initialized = .true.
    end if

    ! --- Interpolate ---
    if (T <= 10.0_wp**(LOGT_MIN)) then
      xHI = 1.0_wp
      return
    end if

    if (T >= 10.0_wp**(LOGT_MAX)) then
      xHI = 10.0_wp**(log_xHI_tab(NTAB))
      return
    end if

    logT = log10(max(T, 1.0_wp))
    frac = (logT - LOGT_MIN) / DLOGT
    idx  = int(frac) + 1
    idx  = max(1, min(idx, NTAB - 1))
    frac = frac - real(idx - 1, wp)

    log_xHI_interp = log_xHI_tab(idx) + frac * (log_xHI_tab(idx+1) - log_xHI_tab(idx))
    xHI = 10.0_wp**(log_xHI_interp)
    xHI = max(0.0_wp, min(1.0_wp, xHI))
  end function cie_neutral_fraction_table

  !--- Helper: build log10(xHI) table from Voronov + Verner rates ---
  subroutine build_cie_table(n, logt_min, dlogt, log_xHI)
    integer, intent(in)   :: n
    real(wp), intent(in)  :: logt_min, dlogt
    real(wp), intent(out) :: log_xHI(n)

    integer  :: i
    real(wp) :: logT_i, T_i, T4, sqrtT, Gamma_HI, alpha_A, xHI_i

    do i = 1, n
      logT_i = logt_min + real(i-1, wp) * dlogt
      T_i    = 10.0_wp**logT_i

      ! Voronov (1997) collisional ionization rate [cm^3/s]
      sqrtT    = sqrt(T_i)
      Gamma_HI = 5.85e-11_wp * sqrtT * exp(-157809.1_wp / T_i) &
               / (1.0_wp + sqrt(T_i / 1.0e5_wp))

      ! Verner & Ferland (1996) Case A recombination rate [cm^3/s]
      T4      = T_i / 1.0e4_wp
      alpha_A = 4.309e-13_wp * T4**(-0.6166_wp) &
              / (1.0_wp + 0.6703_wp * T4**(0.5300_wp))

      xHI_i = alpha_A / (Gamma_HI + alpha_A)
      ! Floor to avoid log10(0)
      xHI_i = max(xHI_i, 1.0e-30_wp)
      log_xHI(i) = log10(xHI_i)
    end do
  end subroutine build_cie_table

end module physics_amr_mod
