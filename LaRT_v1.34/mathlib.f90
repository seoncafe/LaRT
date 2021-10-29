module mathlib
use define, only: wp, sp, dp
use, intrinsic :: iso_fortran_env, only: real32, real64, real128, int32, int64
implicit none

private
public gauleg, dot_product, cross_product, &
       locate, interp, interp_eq, &
       calc_Integral, calc_CDF, &
       array_3D_indices, array_2D_indices

interface dot_product
   module procedure dot_product_r32,dot_product_r64
end interface dot_product
interface cross_product
   module procedure cross_product_r32,cross_product_r64
end interface cross_product

contains
  subroutine gauleg(x1,x2,x,w)
  ! Taken from Numerical Recipes
  real(kind=wp), intent(in) :: x1,x2
  real(kind=wp), dimension(:), intent(out) :: x,w

  real(kind=dp), parameter :: eps = 3.0e-14_dp
  real(kind=dp), parameter :: pi  = 3.141592653589793238462643383279502884197_dp
  integer,       parameter :: maxit = 10
  integer       :: its,j,m,n
  real(kind=dp) :: xl,xm
  real(kind=dp), dimension((size(x)+1)/2) :: p1,p2,p3,pp,z,z1
  logical,       dimension((size(x)+1)/2) :: unfinished

  n  = size(x)
  m  = (n+1)/2
  xm = 0.5_dp*(x2+x1)
  xl = 0.5_dp*(x2-x1)
  do j=1,m
     z(j) = cos(pi*(j-0.25_dp)/(n+0.5_dp))
  enddo
  unfinished = .true.

  do its=1,maxit
     where (unfinished)
        p1 = 1.0
        p2 = 0.0
     end where
     do j=1,n
        where (unfinished)
           p3 = p2
           p2 = p1
           p1 = ((2.0_dp*j-1.0_dp)*z*p2-(j-1.0_dp)*p3)/j
        end where
     end do
     where (unfinished)
        pp = n*(z*p1-p2)/(z*z-1.0_dp)
        z1 = z
        z  = z1-p1/pp
        unfinished = (abs(z-z1) > eps)
     end where
     if (.not. any(unfinished)) exit
  end do

  if (its == maxit+1) write(*,*) 'too many iterations in gauleg'
  x(1:m)        = xm-xl*z
  x(n:n-m+1:-1) = xm+xl*z
  w(1:m)        = 2.0_dp*xl/((1.0_dp-z**2)*pp**2)
  w(n:n-m+1:-1) = w(1:m)

  end subroutine gauleg

  subroutine interp_eq(x,y,xnew,ynew,ix)
  use define, only: wp
  implicit none
  real(kind=wp),     intent(in)  :: x(:),y(:),xnew
  real(kind=wp),     intent(out) :: ynew
  integer, optional, intent(out) :: ix
!---------------
! local variable
  integer :: n, i
  real(kind=wp) :: dx

  n  = size(x)
  dx = x(2) - x(1)
  i  = (xnew-x(1))/dx + 1
  if (i <= 0) then
     ynew = y(1)
  else if (i >= n) then
     ynew = y(n)
  else
     ynew = y(i) + (y(i+1)-y(i))*(xnew-x(i))/dx
  endif
  if (present(ix)) then
     ix = i
  endif

  return
  end subroutine interp_eq

  subroutine interp(x,y,xnew,ynew,ix)
  use define, only: wp
  implicit none
  real(kind=wp), intent(in)  :: x(:),y(:),xnew
  real(kind=wp), intent(out) :: ynew
  integer, optional, intent(out) :: ix
!---------------
! local variable
  integer :: n, i
  logical :: ascend
!
  n = size(x)
  ascend = (x(n) >= x(1))
  call locate(x,xnew,i)
  if (i == 0) then
     ynew = y(1)
  else if (i == n) then
     ynew = y(n)
  else
     if (ascend) then
        ynew = y(i)   + (y(i+1)-y(i))*(xnew-x(i))/(x(i+1)-x(i))
     else
        ynew = y(i+1) + (y(i)-y(i+1))*(xnew-x(i+1))/(x(i)-x(i+1))
     endif
  endif
  if (present(ix)) then
     ix = i
  endif

  return
  end subroutine interp
!=====================================================================
  SUBROUTINE locate(xx,x,j)
  use define, only: wp
  IMPLICIT NONE
  REAL(kind=wp), intent(in) :: xx(:), x
  INTEGER,      intent(out) :: j
! Taken from Numerical Recipes
! Given an array xx(1:n), and given a value x, returns a value j such that
! x is between xx(j) and xx(j+1). xx(1:n) must be monotonic, either increasing
! or decreasing. j = 0 or j =n is returned to indicate that x is out of range.

  INTEGER :: n,jl,jm,ju
  LOGICAL :: ascend

   n = size(xx)
   ascend = (xx(n) >= xx(1))
   jl = 0   ! Initialize lower
   ju = n+1 ! and upper limits.
   do
      if (ju-jl <= 1) exit  ! If we are done,
      jm = (ju+jl)/2        ! compute a midpoint,
      if (ascend .eqv. (x >= xx(jm))) then
         jl = jm            ! and replace either the lower limit
      else
         ju = jm            ! or the upper limit, as appropriate.
      endif
   enddo
   if (x == xx(1)) then ! Then set the output
     j = 1
   else if (x == xx(n)) then
     j = n-1
   else
     j = jl
   endif
   return
  END SUBROUTINE locate
!===========================================================
  !----------------------------
  function dot_product_r32(v1,v2) result(v1dotv2)
    real(kind=real32), intent(in) :: v1(3), v2(3)
    real(kind=real32) :: v1dotv2
    v1dotv2 = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
    return
  end function dot_product_r32
  !----------------------------
  function dot_product_r64(v1,v2) result(v1dotv2)
    real(kind=real64), intent(in) :: v1(3), v2(3)
    real(kind=real64) :: v1dotv2
    v1dotv2 = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
    return
  end function dot_product_r64
  !----------------------------
  function cross_product_r32(v1,v2) result(v1crossv2)
    real(kind=real32), intent(in) :: v1(3), v2(3)
    real(kind=real32) :: v1crossv2(3)
    v1crossv2(1) = v1(2)*v2(3) - v1(3)*v2(2)
    v1crossv2(2) = v1(3)*v2(1) - v1(1)*v2(3)
    v1crossv2(3) = v1(1)*v2(2) - v1(2)*v2(1)
    return
  end function cross_product_r32
  !----------------------------
  function cross_product_r64(v1,v2) result(v1crossv2)
    real(kind=real64), intent(in) :: v1(3), v2(3)
    real(kind=real64) :: v1crossv2(3)
    v1crossv2(1) = v1(2)*v2(3) - v1(3)*v2(2)
    v1crossv2(2) = v1(3)*v2(1) - v1(1)*v2(3)
    v1crossv2(3) = v1(1)*v2(2) - v1(2)*v2(1)
    return
  end function cross_product_r64
!+++++++++++++++++++++++++++++++++++++
  function calc_Integral(x,y)
!------------
! Numerical Integration of a tabular function y(x_i) = y_i for i = 1,...,n
!     2015-01-15, K.I.Seon
!------------
  real(kind=wp), intent(in) :: x(:),y(:)
  real(kind=wp) :: calc_Integral
! local variables
  integer :: n,i

  n = size(x)
  if (n < 2 .or. n > size(y)) then
     write(*,*) 'Error: size(y) >= size(x) >= 2 in function calc_Integral.'
     stop
  endif

  calc_Integral = 0.0_wp
  do i=2,n
    calc_Integral = calc_Integral + 0.5_wp * (y(i)+y(i-1))*abs(x(i)-x(i-1))
  enddo

  return
  end function calc_Integral
!+++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++
  subroutine calc_CDF(x,y,CDF)
  real(kind=wp), intent(in)  :: x(:),y(:)
  real(kind=wp), intent(out) :: CDF(:)
! local variables
  integer :: n,i

  n = size(x)
  if (n < 2 .or. n /= size(y) .or. n /= size(CDF)) then
     write(*,*) 'Error: size(x) = size(y) = size(CDF) >= 1 in function calc_CDF.'
     stop
  endif

  CDF(1) = 0.0_wp
  do i=2,n
    CDF(i) = CDF(i-1) + 0.5_wp*(y(i)+y(i-1))*abs(x(i)-x(i-1))
  enddo
  CDF(:) = CDF(:)/CDF(n)

  return
  end subroutine calc_CDF
!+++++++++++++++++++++++++++++++++++++

  subroutine newton_raphson(func,xroot,x1,x2)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Using the Newton-Raphson method, find the root of a function known to lie in the interval [x1, x2].
! The root xroot will be refined until its accuracy is known within +-xacc.
! func is a user-supplied subroutine that returns both the function value and
!              the first derivative of the function at the point x.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  external                     :: func
  real(kind=wp), intent(inout) :: xroot
  real(kind=wp),    intent(in) :: x1,x2

  integer,       parameter :: JMAX = 20
  real(kind=wp), parameter :: xacc = 1.0e-6_wp
  real(kind=wp) :: df,dx,f
  integer       :: j

  do j=1,JMAX
    call func(xroot,f,df)
    dx    = f/df
    xroot = xroot - dx
    if ((x1-xroot)*(xroot-x2) < 0.0_wp) then
       write(*,*) 'rtnewt - jumped out of brackets'
       stop
    endif
    if (abs(dx) < xacc) return
  enddo
  write(*,*) 'rtnew - exceeded maximum iterations, JMAX'

  end subroutine newton_raphson

  !+++ Logarithmic Module.
  !--- Returns log(a + b) when log(a) and log(b) are given.
  !----      Here, lna = log(a), lnb = log(b), and lnc = log(a+b).
  elemental function logadd(lna,lnb) result(lnc)
  real(kind=wp), intent(in) :: lna, lnb
  real(kind=wp) :: lnc
  if (lna >= lnb) then
     lnc = lna + log(1.0_wp + exp((lnb-lna)))
  else
     lnc = lnb + log(1.0_wp + exp((lna-lnb)))
  endif
  return
  end function logadd

  !--- Returns log(a - b) when log(a) and log(b) are given (a > b).
  elemental function logsub(lna,lnb) result(lnc)
  real(kind=wp), intent(in) :: lna, lnb
  real(kind=wp) :: lnc
  real(kind=wp), parameter :: hugest = huge(1.0_wp)
  if (lna > lnb) then
     lnc = lna + log(1.0_wp - exp(lnb - lna))
  else
     lnc = -hugest
  endif
  end function logsub
!+++++++++++++++++++++++++++++++++++++
  !--------------------
  subroutine array_3D_indices(nx,ny,nz,loc,i,j,k)
  !  This converts one-dimensional subscript (loc) of an array into corresponding 3D subscripts (i,j,k).
  implicit none
  integer, intent(in)  :: nx,ny,nz
  integer, intent(in)  :: loc
  integer, intent(out) :: i,j,k
  integer :: nxy
  nxy = nx*ny
  k   = (loc-1)/nxy+1
  j   = (loc-(k-1)*nxy-1)/nx+1
  i   = loc-(j-1)*nx-(k-1)*nxy
  return
  end subroutine array_3D_indices
  !--------------------
  subroutine array_2D_indices(nx,ny,loc,i,j)
  !  This converts one-dimensional subscript (loc) of an array into corresponding 2D subscripts (i,j).
  implicit none
  integer, intent(in)  :: nx,ny
  integer, intent(in)  :: loc
  integer, intent(out) :: i,j
  j   = (loc-1)/nx+1
  i   = loc-(j-1)*nx
  return
  end subroutine array_2D_indices
  !--------------------
end module mathlib
