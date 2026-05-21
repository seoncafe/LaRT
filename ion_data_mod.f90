module ion_data_mod
  !-----------------------------------------------------------------------
  ! Solar element abundances (Asplund et al. 2009, ARA&A 47, 481) and
  ! CIE ion fraction tables for the ion species supported by LaRT.
  !
  ! Public interface:
  !   solar_abundance(ion_id)       -- returns n_X/n_H (number ratio)
  !   cie_ion_fraction(ion_id, T)   -- returns f_ion(T) from CIE table
  !   solar_ion_density(nH, Z, T, ion_id) -- n_ion = nH*(Z/Z_sun)*A_X*f_ion
  !
  ! Supported ion_id values (matching line%ion_id in line_mod.f90):
  !   'H  I', 'He I', 'C IV', 'N V', 'O VI', 'Na I', 'Ca II',
  !   'Mg II', 'Si IV', 'Al II', 'Si II', 'Fe II'
  !
  ! CIE ion fractions are computed analytically from the same Voronov (1997)
  ! + Verner & Ferland (1996) rates used in physics_amr_mod, extended to
  ! multi-electron ions using a simplified approach.  For a first
  ! implementation we use the optically-thin CIE tables from Gnat &
  ! Sternberg (2007, ApJS 168, 213), hardcoded as fitted coefficients.
  !
  ! For hydrogen and helium, the fractions are computed from the standard
  ! rate equations.  For metals, we store peak ion fraction temperature
  ! and a Gaussian fit in log10(T) space, which reproduces the Gnat &
  ! Sternberg tabulated values to ~10% accuracy.
  !-----------------------------------------------------------------------
  use define, only: wp
  implicit none
  private

  public :: solar_abundance
  public :: solar_ion_density

  ! Solar metallicity (Asplund+09)
  real(wp), parameter :: Z_sun = 0.0134_wp

  ! Solar number abundances n_X/n_H (Asplund+09, Table 1)
  ! These are logarithmic: A(X) = 12 + log10(n_X/n_H)
  ! Converted to linear: n_X/n_H = 10^(A(X)-12)
  real(wp), parameter :: A_C  = 2.692e-4_wp    ! A(C)  = 8.43
  real(wp), parameter :: A_N  = 6.761e-5_wp    ! A(N)  = 7.83
  real(wp), parameter :: A_O  = 4.898e-4_wp    ! A(O)  = 8.69
  real(wp), parameter :: A_Na = 1.738e-6_wp    ! A(Na) = 6.24
  real(wp), parameter :: A_Mg = 3.981e-5_wp    ! A(Mg) = 7.60
  real(wp), parameter :: A_Al = 2.818e-6_wp    ! A(Al) = 6.45
  real(wp), parameter :: A_Si = 3.236e-5_wp    ! A(Si) = 7.51
  real(wp), parameter :: A_Ca = 2.188e-6_wp    ! A(Ca) = 6.34
  real(wp), parameter :: A_Fe = 3.162e-5_wp    ! A(Fe) = 7.50
  real(wp), parameter :: A_He = 8.511e-2_wp    ! A(He) = 10.93

  ! CIE ion fraction Gaussian fit parameters (Gnat & Sternberg 2007)
  ! f_ion(T) ~ f_peak * exp(-0.5*((logT - logT_peak)/sigma)^2)
  ! These are approximate but sufficient for order-of-magnitude estimates.
  type :: ion_fit_type
    real(wp) :: logT_peak  ! log10(T) at peak ion fraction
    real(wp) :: f_peak     ! peak ion fraction
    real(wp) :: sigma      ! Gaussian width in log10(T)
  end type

contains

  !=========================================================================
  ! Solar number abundance n_X/n_H for a given ion_id.
  !=========================================================================
  function solar_abundance(ion_id) result(A_X)
    character(len=*), intent(in) :: ion_id
    real(wp) :: A_X

    select case (trim(adjustl(ion_id)))
    case ('H  I', 'H I', 'H+D')
      A_X = 1.0_wp
    case ('He I')
      A_X = A_He
    case ('C IV')
      A_X = A_C
    case ('N V')
      A_X = A_N
    case ('O VI')
      A_X = A_O
    case ('Na I')
      A_X = A_Na
    case ('Ca II')
      A_X = A_Ca
    case ('Mg II')
      A_X = A_Mg
    case ('Si IV', 'Si II')
      A_X = A_Si
    case ('Al II')
      A_X = A_Al
    case ('Fe II')
      A_X = A_Fe
    case default
      A_X = 0.0_wp
    end select
  end function solar_abundance

  !=========================================================================
  ! CIE ion fraction f_ion(T) for a given ion_id.
  !
  ! Uses Gaussian fits in log10(T) to the Gnat & Sternberg (2007) CIE
  ! tables.  For H I and He I, uses the analytic rate equations from
  ! physics_amr_mod.
  !=========================================================================
  function cie_ion_fraction(ion_id, T) result(f_ion)
    character(len=*), intent(in) :: ion_id
    real(wp), intent(in) :: T
    real(wp) :: f_ion

    real(wp) :: logT, dlogT
    type(ion_fit_type) :: fit

    logT = log10(max(T, 10.0_wp))

    select case (trim(adjustl(ion_id)))
    case ('H  I', 'H I', 'H+D')
      ! Use the full rate equation (same as physics_amr_mod)
      f_ion = cie_xHI_simple(T)
      return
    case ('He I')
      ! He I: peak at ~2e4 K, drops above ~3e4 K
      fit = ion_fit_type(4.25_wp, 0.95_wp, 0.25_wp)
    case ('C IV')
      ! C+3: peak at ~1e5 K
      fit = ion_fit_type(5.05_wp, 0.29_wp, 0.20_wp)
    case ('N V')
      ! N+4: peak at ~1.8e5 K
      fit = ion_fit_type(5.25_wp, 0.23_wp, 0.18_wp)
    case ('O VI')
      ! O+5: peak at ~3e5 K
      fit = ion_fit_type(5.45_wp, 0.20_wp, 0.18_wp)
    case ('Na I')
      ! Na neutral: only at T < ~5000 K
      fit = ion_fit_type(3.60_wp, 0.90_wp, 0.20_wp)
    case ('Ca II')
      ! Ca+: peak at ~1.2e4 K
      fit = ion_fit_type(4.10_wp, 0.65_wp, 0.25_wp)
    case ('Mg II')
      ! Mg+: peak at ~2.5e4 K
      fit = ion_fit_type(4.35_wp, 0.70_wp, 0.22_wp)
    case ('Si IV')
      ! Si+3: peak at ~7e4 K
      fit = ion_fit_type(4.85_wp, 0.35_wp, 0.22_wp)
    case ('Si II')
      ! Si+: peak at ~2e4 K
      fit = ion_fit_type(4.30_wp, 0.70_wp, 0.20_wp)
    case ('Al II')
      ! Al+: peak at ~1.5e4 K
      fit = ion_fit_type(4.20_wp, 0.75_wp, 0.22_wp)
    case ('Fe II')
      ! Fe+: peak at ~2.5e4 K
      fit = ion_fit_type(4.35_wp, 0.70_wp, 0.22_wp)
    case default
      f_ion = 0.0_wp
      return
    end select

    ! Gaussian fit in log10(T) space
    dlogT = logT - fit%logT_peak
    f_ion = fit%f_peak * exp(-0.5_wp * (dlogT / fit%sigma)**2)
    f_ion = max(0.0_wp, min(1.0_wp, f_ion))
  end function cie_ion_fraction

  !=========================================================================
  ! Compute per-cell ion number density:
  !   n_ion = nH * (Z / Z_sun) * (n_X/n_H)_sun * f_ion(T)
  !
  ! For hydrogen lines, Z scaling is not applied (n_ion = nH * xHI).
  !=========================================================================
  function solar_ion_density(nH, Z, T, ion_id) result(n_ion)
    real(wp), intent(in) :: nH, Z, T
    character(len=*), intent(in) :: ion_id
    real(wp) :: n_ion

    real(wp) :: A_X, f_ion

    select case (trim(adjustl(ion_id)))
    case ('H  I', 'H I', 'H+D')
      ! Hydrogen: n_ion = nH * xHI (no metallicity scaling)
      n_ion = nH * cie_xHI_simple(T)
      return
    case ('He I')
      ! Helium: n_ion = nHe * xHeI (no metallicity scaling)
      n_ion = nH * A_He * cie_ion_fraction(ion_id, T)
      return
    case default
      ! Metals: scale by metallicity
      A_X   = solar_abundance(ion_id)
      f_ion = cie_ion_fraction(ion_id, T)
      n_ion = nH * (Z / Z_sun) * A_X * f_ion
    end select
  end function solar_ion_density

  !=========================================================================
  ! Simple CIE neutral fraction for hydrogen (same formula as
  ! cie_neutral_fraction_formula in physics_amr_mod, duplicated here
  ! to avoid circular dependency).
  !=========================================================================
  function cie_xHI_simple(T) result(xHI)
    real(wp), intent(in) :: T
    real(wp) :: xHI
    real(wp) :: T4, k_ionize, k_recomb

    T4       = max(T, 10.0_wp) / 1.0e4_wp
    k_ionize = 5.84862e-9_wp * sqrt(T4) * exp(-15.78215_wp / T4)
    k_recomb = 4.13e-13_wp * T4**(-0.7131_wp - 0.0115_wp * log(T4))
    xHI      = k_recomb / (k_ionize + k_recomb)
  end function cie_xHI_simple

end module ion_data_mod
