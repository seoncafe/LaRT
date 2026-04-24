module random_sersic
  use define,  only : wp, dp, pi, halfpi
  use mathlib, only : interp, locate
  use random,  only : rand_number
contains
  !=============================
  ! Author: Kwang-Il Seon
  ! 2010-10-10
  !    included functions for special functions, such as gamma function.
  ! 2016-02-02
  !    bug fixed: sersic_dlnx should be declared as a save variable.
  ! 2016-01-26
  !    Optional rmax parameter is added. Interpolation at r ~ 0 is updated.
  ! 2015-12-11
  !     Correct way to make saved variables thread safe is to define them as THREADPRIVATE.
  !     http://sc.tamu.edu/IBM.Tutorial/docs/Compilers/xlf_8.1/html/lr225.HTM
  ! 2015-11-21
  !     Local variables with the SAVE attribute declared in procedures called from a parallel region are implicitly SHARED.
  !     http://www.openmp.org/mp-documents/finterp10.html#interp3
  !     Therefore, "SAVE"ed variables should have unique names,
  ! 2015-11-19
  !     Underestimation is caused by low-resolution radial bin in the final profile in rand_sersic.
  ! 2015-11-17
  !     Numerical method seems to underestimates at very small radia.
  !     In practical sense, underestimation is negligible.
  !
  ! Calculate deprojected 3D cumulative luminosity function from 2D Sersic profile.
  !     Sersic profile: I(R) = I_0 exp(-b(R/R_e)^(1/m)) where R is projected 2D radius.
  !     m       = Sersic Index (m >= 0.5)
  !     radius  = R/R_e (R_e = effective radius containing half of luminosity in projected 2D image)
  !     profile = normalized deprojected 3D luminosity profile
  !
  ! Probably, the algorithm will be described in Seon et al. (2021...).
  !
  subroutine sersic_cumulative_3D(m, radius, profile, rmax_in)
    implicit none
    real(wp), intent(in)    :: m
    real(wp), intent(inout) :: radius(:), profile(:)
    real(wp), optional, intent(in) :: rmax_in

    ! Local variables
    integer,  parameter :: nx   = 20001
    real(dp), parameter :: xmax = 10000
    real(dp), save, allocatable :: sersic_x(:), sersic_fx(:)
    logical,  save :: sersic_fx_saved = .false.
    real(dp), save :: sersic_dlnx
    !$OMP THREADPRIVATE(sersic_x,sersic_fx,sersic_fx_saved,sersic_dlnx)

    real(dp) :: m2
    real(dp) :: fr
    real(dp) :: sum1, sum2, norm
    real(dp) :: rmin,rmax,dlnr
    integer  :: nr,nr1
    integer  :: i,j

    real(dp), parameter :: bcoeff(4) = (/ 4d0/405d0, 46d0/25515d0, 131d0/1148175d0, -2194697d0/30690717750d0 /)
    real(dp) :: r, b

    !--- b as a function of sersic index (Ciotti & Bertin, 1999, A&A, 352, 447)
    b = 0.0d0
    !do i=4,1,-1
    do i=2,1,-1
       b = (b + bcoeff(i))/m
    enddo
    b  = b + 2d0*m - 1d0/3d0
    m2 = 2d0*m

    !--- Determine radius range, which includes most of luminosity (or mass).
    !--- Probably, nice formulae will be found later.
    !--- Radius ranges including 0.001 and 99.999% for m = 0.5-10 were numerically calculated using Mathematica,
    !---      and then the calculated rmin and rmax values were fitted using a functional form.
    nr   = size(radius)
    nr1  = nr - 1
    rmax = 2.5d0*((14.995674d0+4.0947738d0*m-0.052804581d0*m*m)/b)**m
    rmin = -0.27566682d0+0.21713972d0*m+0.037967891d0*m*m
    if (present(rmax_in)) then
       if (rmax_in > 0.0_wp) rmax = rmax_in
    endif
    if (rmin < 0.0d0) then
       rmin = rmax/1d4
    else
       rmin = 0.005d0*(rmin/b)**m
    endif
    dlnr = log(rmax/rmin)/(nr1-1)

    !--- Integration formula is obtained using Abel Transform
    !---      and changing the order of integration variables using Fubini's theorem.
    if (.not.(sersic_fx_saved)) then
       sersic_dlnx = log(xmax)/dble(nx-1)
       if (.not. allocated(sersic_x))   allocate(sersic_x(nx))
       if (.not. allocated(sersic_fx))  allocate(sersic_fx(nx))
       do i=1,nx
          sersic_x(i)   = exp((i-1)*sersic_dlnx)
          if (sersic_x(i) == 1.d0) then
             sersic_fx(i) = halfpi
          else
             sersic_fx(i) = -sqrt(1d0-1d0/sersic_x(i)**2) + sersic_x(i)*atan(1d0/sqrt(sersic_x(i)**2-1d0))
          endif
       enddo
       sersic_fx_saved = .true.
    endif

    radius(1)  = 0.0_dp
    profile(1) = 0.0_dp
    do j=1,nr1
       r    = rmin * exp((j-1)*dlnr)
       sum1 = gammp(m2+1d0, b*r**(1d0/m))

       sum2 = 0.0d0
       do i=1,nx
          fr = exp(-b*(r*sersic_x(i))**(1d0/m)) * (r*sersic_x(i))**(1d0/m)
          ! trapezoidal rule
          if (i == 1 .or. i == nx) then
             sum2  = sum2 + sersic_x(i)*fr*sersic_fx(i)*0.5d0
          else
             sum2  = sum2 + sersic_x(i)*fr*sersic_fx(i)
          endif
       enddo
       norm = (2d0/pi) * (b**(m2+1))/m * exp(-gammln(m2+1d0))
       sum2 = norm * r**2 * sum2*sersic_dlnx

       radius(j+1)  = r
       profile(j+1) = sum1 + sum2
    enddo

  end subroutine sersic_cumulative_3D
  !=============================
  function rand_sersic(m,rmax) result(ranv)
  !--- Random number generator for Sersic Profile
  !--- Optional rmax parameter should be given in unit of Reff.
     implicit none
     real(wp), intent(in) :: m
     real(wp), optional, intent(in) :: rmax
     real(wp) :: ranv

     integer  :: ip
     real(dp) :: pv
     real(dp), save, allocatable :: sersic_r(:), sersic_p(:)
     real(dp), save :: sersic_m_saved    = -huge(0.0_dp)
     integer,  save :: sersic_nr         = 100
     real(dp), save :: sersic_rmax_saved = -huge(0.0_wp)
     logical :: initialize
     !$OMP THREADPRIVATE(sersic_r,sersic_p,sersic_m_saved,sersic_rmax_saved,sersic_nr)

     initialize = .false.
     if (m /= sersic_m_saved) initialize = .true.
     if (present(rmax)) then
         if (sersic_rmax_saved /= rmax) initialize = .true.
     endif

     if (initialize) then
        if (.not. allocated(sersic_r))  allocate(sersic_r(sersic_nr))
        if (.not. allocated(sersic_p))  allocate(sersic_p(sersic_nr))
        call sersic_cumulative_3D(m,sersic_r,sersic_p,rmax_in=rmax)
        do ip=2,sersic_nr
           sersic_p(ip) = sersic_p(ip)/sersic_p(sersic_nr)
        enddo
        sersic_m_saved    = m
        sersic_rmax_saved = rmax
     endif

     !call random_number(pv)
     pv = rand_number()
     call locate(sersic_p,pv,ip)
     !--call locate(sersic_p,sersic_nr,pv,ip)
     !if (.not.(pv >= sersic_p(ip) .and. pv < sersic_p(ip+1))) print*,'something wrong in rand_sersic.'
     if (ip == 1) then
        ! linear interpolation
        !ranv = sersic_r(ip+1)/sersic_p(ip+1) * pv
        ! Note: sersic_p is proportional to r^(1/m + 2) as r -> 0.
        ranv = sersic_r(ip+1) * (pv/sersic_p(ip+1))**(1d0/(1d0/m+2.d0))
     else
        ! linear interpolation
        !ranv = (sersic_r(ip+1)-sersic_r(ip))/(sersic_p(ip+1)-sersic_p(ip)) * (pv - sersic_p(ip)) + sersic_r(ip)
        ! log-linear interpolation
        !ranv = log(sersic_r(ip+1)/sersic_r(ip))/(sersic_p(ip+1)-sersic_p(ip)) * (pv - sersic_p(ip)) + log(sersic_r(ip))
        ! log-log interpolation
        ranv = log(sersic_r(ip+1)/sersic_r(ip))/log(sersic_p(ip+1)/sersic_p(ip)) * log(pv/sersic_p(ip)) + log(sersic_r(ip))
        ranv = exp(ranv)
     endif

  end function rand_sersic
  !+++++++++++++++++++++++++++++++++++++
  !--------
  !  These routines are from Numerical Recipes.
  !    2015-11-16
  !-------- Gammln
  FUNCTION gammln(xx)
  IMPLICIT NONE
  REAL(WP), INTENT(IN) :: xx
  REAL(DP) :: gammln
  REAL(DP) :: tmp,x
  REAL(DP) :: ser,y
  REAL(DP) :: stp = 2.5066282746310005_dp
  REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
                                    -86.50532032941677_dp,&
                                     24.01409824083091_dp,&
                                     -1.231739572450155_dp,&
                                      0.1208650973866179e-2_dp,&
                                     -0.5395239384953e-5_dp/)
  INTEGER :: i
  x   = xx
  tmp = x+5.5_dp
  tmp = (x+0.5_dp)*log(tmp)-tmp
  ser = 1.000000000190015_dp
  y   = x
  do i=1, size(coef)
     y   = y+1.0_dp
     ser = ser+coef(i)/y
  end do
  gammln = tmp + log(stp*ser/x)
  END FUNCTION gammln
  !-------- Gcf
  FUNCTION gcf(a,x,gln)
  IMPLICIT NONE
  REAL(WP), INTENT(IN) :: a,x
  REAL(WP), OPTIONAL, INTENT(OUT) :: gln
  REAL(DP) :: gcf
  INTEGER,  PARAMETER :: ITMAX=100
  REAL(DP), PARAMETER :: EPS=epsilon(x), FPMIN=tiny(x)/EPS
  INTEGER :: i
  REAL(DP) :: an,b,c,d,del,h
  if (x == 0.0) then
     gcf = 1.0
     RETURN
  end if
  b = x+1.0_dp-a
  c = 1.0_dp/FPMIN
  d = 1.0_dp/b
  h = d
  do i=1,ITMAX
     an = -i*(i-a)
     b  = b+2.0_dp
     d  = an*d+b
     if (abs(d) < FPMIN) d = FPMIN
     c  = b+an/c
     if (abs(c) < FPMIN) c = FPMIN
     d   = 1.0_dp/d
     del = d*c
     h   = h*del
     if (abs(del-1.0_dp) <= EPS) exit
  end do
  if (i > ITMAX) then
     print*,'a too large, ITMAX too small in gcf'
     stop
  endif
  if (present(gln)) then
     gln = gammln(a)
     gcf = exp(-x+a*log(x)-gln)*h
  else
     gcf = exp(-x+a*log(x)-gammln(a))*h
  end if
  END FUNCTION gcf
  !-------- Gser
  FUNCTION gser(a,x,gln)
  IMPLICIT NONE
  REAL(WP), INTENT(IN) :: a,x
  REAL(WP), OPTIONAL, INTENT(OUT) :: gln
  REAL(DP) :: gser
  INTEGER,  PARAMETER :: ITMAX=100
  REAL(DP), PARAMETER :: EPS=epsilon(x)
  INTEGER  :: n
  REAL(DP) :: ap,del,summ
  !
  if (x == 0.0) then
     gser = 0.0
     RETURN
  end if
  ap   = a
  summ = 1.0_dp/a
  del  = summ
  do n=1, ITMAX
     ap   = ap+1.0_dp
     del  = del*x/ap
     summ = summ+del
     if (abs(del) < abs(summ)*EPS) exit
  end do
  if (n > ITMAX) then
     print*,'a too large, ITMAX too small in gser'
     stop
  endif
  if (present(gln)) then
     gln  = gammln(a)
     gser = summ*exp(-x+a*log(x)-gln)
  else
     gser = summ*exp(-x+a*log(x)-gammln(a))
  end if
  END FUNCTION gser
  !-------- GammaP
  FUNCTION gammp(a,x)
  IMPLICIT NONE
  REAL(WP), INTENT(IN) :: a,x
  REAL(DP) :: gammp
  if (x < a+1.0_dp) then
     gammp = gser(a,x)
  else
     gammp = 1.0_dp-gcf(a,x)
  end if
  END FUNCTION gammp
  !=============================
end module random_sersic
