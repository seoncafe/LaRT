module qmc_mod
!---------------------------------------------------------------------------
! Quasi-random (randomized quasi-Monte Carlo) photon launch:
! par%launch_sequence = 'sobol' draws the launch variables of a photon packet
! and its first-flight optical depth from an Owen-scrambled Sobol point, while
! every later transport draw stays on the Mersenne Twister.
! par%launch_sequence = 'random' (default) calls rand_number() in the same
! places and order as before, so runs are unchanged.
!
! The Sobol point of a packet is addressed by its global photon number
! (photon%id = 1..nphotons), so the launch set does not depend on the MPI task
! count or on the master-slave work assignment.  The Owen scramble is keyed by
! par%qmc_seed, identical on every rank; a different seed gives an independent
! randomized replicate (successive par%out_merge runs must use different seeds).
!
! Fixed dimension layout (a dimension keeps its meaning for every packet;
! a dimension a source does not use is left unused, which is harmless):
!   1, 2  frequency: Gaussian (Box-Muller), or the continuum position (dim 1)
!   3     frequency: Cauchy part of a Voigt source, or the line/continuum
!         choice of 'continuum+gaussian'
!   4     upward component of a multiplet (line_type 2, 3, 5, 6)
!   5, 6  isotropic emission direction (cos theta, phi)
!   7-9   source position ('uniform_sphere'/'sphere', 'uniform_cylinder'/
!         'cylinder', 'uniform', 'uniform_xy'); other geometries use the
!         Mersenne Twister for the position
!   10    first-flight optical depth (forced first scattering)
! Draws that use rejection or an alias table (line_prof_file, the downward
! branch of line_type 4/5, star_file, illumination and emissivity sources)
! stay on the Mersenne Twister.
!
! Sobol points and the Owen scramble are taken from MoCHII_v1.00
! (src/qmc_mod.f90, K. Seon): dimension 1 is van der Corput base 2,
! dimensions 2-12 use the Joe & Kuo (2008, new-joe-kuo-6) direction numbers
! (they reproduce scipy.stats.qmc.Sobol); the scramble is the hash-based Owen
! scramble of Burley (2020, JCGT 9, 1).
!---------------------------------------------------------------------------
  use, intrinsic :: iso_fortran_env, only : int64
  use define, only : wp, twopi, pi
  use random, only : rand_number, rand_gauss, rand_voigt
  implicit none
  private

  integer, parameter, public :: QMC_D_FREQ1 = 1
  integer, parameter, public :: QMC_D_FREQ2 = 2
  integer, parameter, public :: QMC_D_FREQ3 = 3
  integer, parameter, public :: QMC_D_COMP  = 4
  integer, parameter, public :: QMC_D_MU    = 5
  integer, parameter, public :: QMC_D_PHI   = 6
  integer, parameter, public :: QMC_D_POS1  = 7
  integer, parameter, public :: QMC_D_POS2  = 8
  integer, parameter, public :: QMC_D_POS3  = 9
  integer, parameter, public :: QMC_D_TAU1  = 10
  integer, parameter, public :: QMC_NDIM    = 10

  integer, parameter :: QMC_MAXDIM = 12
  integer, parameter :: NBITS = 32
  integer(int64), parameter :: MASK32 = int(z'FFFFFFFF', int64)
  real(wp),       parameter :: TWO_M32 = 1.0_wp / 4294967296.0_wp   ! 2^-32

  !--- .true. for par%launch_sequence = 'sobol' (set by qmc_setup).
  logical, public, save :: qmc_launch = .false.

  integer(int64), save :: Vdir(NBITS, QMC_MAXDIM)
  integer(int64), save :: dkey(QMC_MAXDIM)
  !--- scrambled Sobol point of the packet being launched (pure MPI: one
  !--- packet at a time on each rank).
  real(wp),       save :: u_launch(QMC_NDIM)

  public :: qmc_setup, qmc_launch_point, qmc_uniforms_raw
  public :: launch_uniform, launch_gauss, launch_voigt

contains

  !=========================================================================
  ! Direction numbers for dimensions 1..QMC_MAXDIM and the Owen scramble keys
  ! from the seed (the same on every rank).
  !=========================================================================
  subroutine qmc_setup(seed)
    implicit none
    integer, intent(in) :: seed
    !--- Joe & Kuo (2008) new-joe-kuo-6: sdeg = polynomial degree s, apoly =
    !--- interior-coefficient bit pattern a, minit = initial direction integers.
    integer, parameter :: sdeg(2:QMC_MAXDIM)  = [1, 2, 3, 3, 4, 4, 5, 5, 5, 5, 5]
    integer, parameter :: apoly(2:QMC_MAXDIM) = [0, 1, 1, 2, 1, 4, 2, 4, 7, 11, 13]
    integer :: minit(QMC_MAXDIM, 5)
    integer :: d, k, j, s, a
    integer(int64) :: v

    minit = 0
    minit(2, 1:1) = [1]
    minit(3, 1:2) = [1, 3]
    minit(4, 1:3) = [1, 3, 1]
    minit(5, 1:3) = [1, 1, 1]
    minit(6, 1:4) = [1, 1, 3, 3]
    minit(7, 1:4) = [1, 3, 5, 13]
    minit(8, 1:5) = [1, 1, 5, 5, 17]
    minit(9, 1:5) = [1, 1, 5, 5,  5]
    minit(10,1:5) = [1, 1, 7, 11, 19]
    minit(11,1:5) = [1, 1, 5, 1,  1]
    minit(12,1:5) = [1, 1, 1, 3, 11]

    !--- dimension 1: van der Corput base 2, V(k) = 2^(32-k).
    do k = 1, NBITS
       Vdir(k, 1) = ishft(1_int64, NBITS - k)
    end do

    !--- dimensions 2..12: V(k) = m_k << (32-k) for k <= s, and for k > s
    !---   V(k) = V(k-s) XOR (V(k-s) >> s) XOR [ sum_{j=1}^{s-1} a_j V(k-j) ],
    !--- with a_j = bit (s-1-j) of a.
    do d = 2, QMC_MAXDIM
       s = sdeg(d)
       a = apoly(d)
       do k = 1, s
          Vdir(k, d) = iand(ishft(int(minit(d, k), int64), NBITS - k), MASK32)
       end do
       do k = s + 1, NBITS
          v = ieor(Vdir(k-s, d), ishft(Vdir(k-s, d), -s))
          do j = 1, s - 1
             if (iand(ishft(int(a, int64), -(s-1-j)), 1_int64) == 1_int64) &
                v = ieor(v, Vdir(k-j, d))
          end do
          Vdir(k, d) = iand(v, MASK32)
       end do
    end do

    do d = 1, QMC_MAXDIM
       dkey(d) = hash32(iand(int(seed, int64) + int(d, int64) * int(z'9E3779B9', int64), MASK32))
    end do
    qmc_launch = .true.
  end subroutine qmc_setup

  !=========================================================================
  ! Store the scrambled Sobol point of global photon number id (1-based).
  ! Each coordinate is strictly inside (0,1).
  !=========================================================================
  subroutine qmc_launch_point(id)
    implicit none
    integer(int64), intent(in) :: id
    integer :: d
    integer(int64) :: code
    do d = 1, QMC_NDIM
       code = owen_scramble(sobol_code(id - 1_int64, d), dkey(d))
       u_launch(d) = (real(code, wp) + 0.5_wp) * TWO_M32
    end do
  end subroutine qmc_launch_point

  !=========================================================================
  ! Unscrambled uniforms of the 0-based index idx (check against scipy).
  !=========================================================================
  subroutine qmc_uniforms_raw(idx, u)
    implicit none
    integer(int64), intent(in)  :: idx
    real(wp),       intent(out) :: u(:)
    integer :: d
    do d = 1, size(u)
       u(d) = (real(sobol_code(idx, d), wp) + 0.5_wp) * TWO_M32
    end do
  end subroutine qmc_uniforms_raw

  !=========================================================================
  ! Launch draws: the Sobol coordinate of dimension d, or rand_number().
  !=========================================================================
  function launch_uniform(d) result(u)
    implicit none
    integer, intent(in) :: d
    real(wp) :: u
    if (qmc_launch) then
       u = u_launch(d)
    else
       u = rand_number()
    endif
  end function launch_uniform

  !--- standard normal deviate: Box-Muller on dimensions 1 and 2, or rand_gauss().
  function launch_gauss() result(g)
    implicit none
    real(wp) :: g
    if (qmc_launch) then
       g = sqrt(-2.0_wp*log(u_launch(QMC_D_FREQ1))) * cos(twopi*u_launch(QMC_D_FREQ2))
    else
       g = rand_gauss()
    endif
  end function launch_gauss

  !--- Voigt deviate a*Cauchy + Gauss/sqrt(2) (as rand_voigt): Cauchy by
  !--- inversion of dimension 3, Gaussian from dimensions 1 and 2.
  function launch_voigt(a) result(x)
    implicit none
    real(wp), intent(in) :: a
    real(wp) :: x
    real(wp), parameter :: one_over_sqrt2 = 1.0_wp/sqrt(2.0_wp)
    if (qmc_launch) then
       x = a*tan(pi*(u_launch(QMC_D_FREQ3) - 0.5_wp)) + launch_gauss()*one_over_sqrt2
    else
       x = rand_voigt(a)
    endif
  end function launch_voigt

  !=========================================================================
  ! Unscrambled Sobol integer: x = XOR over the set bits k of idx of V_d(k).
  !=========================================================================
  pure function sobol_code(idx, d) result(code)
    implicit none
    integer(int64), intent(in) :: idx
    integer,        intent(in) :: d
    integer(int64) :: code, ii
    integer :: k
    code = 0_int64
    ii   = idx
    k    = 1
    do while (ii /= 0_int64 .and. k <= NBITS)
       if (iand(ii, 1_int64) == 1_int64) code = ieor(code, Vdir(k, d))
       ii = ishft(ii, -1)
       k  = k + 1
    end do
  end function sobol_code

  !--- Owen (nested-uniform) scramble by a keyed hash (Burley 2020):
  !--- reverse bits, Laine-Karras permutation, reverse back.
  pure function owen_scramble(x, seed) result(y)
    implicit none
    integer(int64), intent(in) :: x, seed
    integer(int64) :: y
    y = reverse32(x)
    y = laine_karras(y, seed)
    y = reverse32(y)
  end function owen_scramble

  pure function laine_karras(x0, seed) result(x)
    implicit none
    integer(int64), intent(in) :: x0, seed
    integer(int64) :: x
    x = iand(x0 + seed, MASK32)
    x = iand(ieor(x, iand(x * int(z'6C50B47C', int64), MASK32)), MASK32)
    x = iand(ieor(x, iand(x * int(z'B82F1E52', int64), MASK32)), MASK32)
    x = iand(ieor(x, iand(x * int(z'C7AFE638', int64), MASK32)), MASK32)
    x = iand(ieor(x, iand(x * int(z'8D22F6E6', int64), MASK32)), MASK32)
  end function laine_karras

  pure function reverse32(n0) result(n)
    implicit none
    integer(int64), intent(in) :: n0
    integer(int64) :: n
    n = iand(n0, MASK32)
    n = iand(ior(ishft(n, 16), ishft(n, -16)), MASK32)
    n = iand(ior(ishft(iand(n, int(z'00FF00FF', int64)),  8),  &
                 ishft(iand(n, int(z'FF00FF00', int64)), -8)), MASK32)
    n = iand(ior(ishft(iand(n, int(z'0F0F0F0F', int64)),  4),  &
                 ishft(iand(n, int(z'F0F0F0F0', int64)), -4)), MASK32)
    n = iand(ior(ishft(iand(n, int(z'33333333', int64)),  2),  &
                 ishft(iand(n, int(z'CCCCCCCC', int64)), -2)), MASK32)
    n = iand(ior(ishft(iand(n, int(z'55555555', int64)),  1),  &
                 ishft(iand(n, int(z'AAAAAAAA', int64)), -1)), MASK32)
  end function reverse32

  !--- lowbias32 integer avalanche hash (Wellons), for the scramble keys.
  pure function hash32(x0) result(x)
    implicit none
    integer(int64), intent(in) :: x0
    integer(int64) :: x
    x = iand(x0, MASK32)
    x = iand(ieor(x, ishft(x, -16)), MASK32)
    x = iand(x * int(z'7FEB352D', int64), MASK32)
    x = iand(ieor(x, ishft(x, -15)), MASK32)
    x = iand(x * int(z'846CA68B', int64), MASK32)
    x = iand(ieor(x, ishft(x, -16)), MASK32)
  end function hash32

end module qmc_mod
