module healpix
 use, intrinsic :: iso_fortran_env, only: real32, real64, real128, int32, int64
 integer, private, parameter :: wp = real64
contains
 function nside2npix(nside) result(npix)
 !=======================================================================
 ! given nside, returns npix such that npix = 12*nside^2
 !  nside should be a power of 2 smaller than 8192
 !  if not, -1 is returned
 ! EH, Feb-2000
 !=======================================================================
 integer, intent(in) :: nside
 integer :: npix
 integer :: listnside(14)
 !=======================================================================

 listnside = (/1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192/)
 npix      = 12*nside*nside

 if (minval(abs(listnside-nside)) > 0) then
    print*,"nside2npix : invalid nside ",nside
    print "(a,4(i2),3(i3),3(i4),4(i5))", " Nside should be among ",listnside
    npix = -1
 endif

 return
 end function nside2npix
!-----------------------------------------------------------------------
 subroutine vec2pix(nside,vx,vy,vz,ipix)
!=======================================================================
!     renders the pixel number ipix (RING scheme) for a pixel which contains
!     a point on a sphere at coordinate vector (=vx,vy,vz), given the map
!     resolution parameter nside
!
!     Translated from HEALPIX fortran 90 routine by Kwangil Seon
!=======================================================================
 implicit none
 integer, intent(in) :: nside
 real(kind=wp), intent(in) :: vx,vy,vz
 integer, intent(out) :: ipix

 integer :: nl2, nl4, ncap, npix, jp, jm, ipix1
 real(kind=wp) :: z, za, tt, tp, tmp, dnorm, phi
 integer :: ir, ip, kshift

 integer, parameter :: ns_max=8192

 real(kind=wp), parameter :: pi=3.141592653589793238462643383279502884197d0,&
                             twopi=6.283185307179586476925286766559005768394d0,&
                             fourpi=12.56637061435917295385057353311801153679d0,&
                             halfpi=1.570796326794896619231321691639751442099d0,&
                             twothird=0.6666666666666666666666666666666666666666d0
!-----------------------------------------------------------------------
 if (nside.lt.1 .or. nside.gt.ns_max) then
    write(*,*) 'nside out of range'
    stop
 endif

 dnorm = sqrt(vx**2+vy**2+vz**2)
 z = vz / dnorm
 phi = 0.0
 if (vx .ne. 0.0 .or. vy .ne. 0.0) then
    phi = atan2(vy,vx)                ! phi in ]-pi,pi]
 endif

 za = abs(z)
 if (phi .lt. 0.0)  phi = phi + twopi ! phi in [0,2pi[
 tt = phi / halfpi                    ! in [0,4)

 nl2  = 2*nside
 nl4  = 4*nside
 ncap = nl2*(nside-1)                 ! number of pixels in the north polar cap
 npix = 12*nside**2

 if ( za .le. twothird ) then           ! Equatorial region ------------------
    jp = int(nside*(0.5 + tt - z*0.75)) ! index of  ascending edge line
    jm = int(nside*(0.5 + tt + z*0.75)) ! index of descending edge line

    ir = nside + 1 + jp - jm            ! in {1,2n+1} (ring number counted from z=2/3)
    kshift = 0
    if (modulo(ir,2) .eq. 0) kshift = 1 ! kshift=1 if ir even, 0 otherwise
    ! note that MOD(A,B) = A - int(A/B) * B and MODULO(A,B) = A - floor(A/B) * B.

    ip = int( ( jp+jm - nside + kshift + 1 ) / 2 ) + 1 ! in {1,4n}
    if (ip .gt. nl4) ip = ip - nl4

    ipix1 = ncap + nl4*(ir-1) + ip
 else                                    ! North & South polar caps --------------
    tp = tt - int(tt)                    ! modulo(tt,1.0)
    tmp = nside * sqrt( 3.0*(1.0 - za) )

    jp = int( tp         * tmp )         ! increasing edge line index
    jm = int( (1.0 - tp) * tmp )         ! decreasing edge line index

    ir = jp + jm + 1                     ! ring number counted from the closest pole
    ip = int( tt * ir ) + 1              ! in {1,4*ir}
    if (ip .gt. 4*ir) ip = ip - 4*ir

    if (z .le. 0.0) then
       ipix1 = npix - 2*ir*(ir+1) + ip
    else
       ipix1 = 2*ir*(ir-1) + ip
    endif
 endif

! ipix = ipix1 - 1           ! in {0, npix-1}
 ipix = ipix1                ! in {1, npix}
! if (ipix.lt.1)    write(*,*) 'Something wrong ipix < 1',ipix
! if (ipix.gt.npix) write(*,*) 'Something wrong ipix > npix',ipix

 return
 end subroutine vec2pix
!-----------------------------------------------------------------------
 subroutine pix2vec(nside, ipix, vx,vy,vz)
!=======================================================================
!  renders vector (x,y,z) coordinates of the nominal pixel center
!  for the pixel number ipix (RING scheme)
!  given the map resolution parameter nside
!  also returns the (x,y,z) position of the 4 pixel vertices (=corners)
!  in the order N,W,S,E
!=======================================================================
 implicit none
 integer,       intent(in)  :: nside, ipix
 real(kind=wp), intent(out) :: vx,vy,vz

 integer       :: nl2, nl4, npix, ncap, iring, iphi, ip, ipix1
 real(kind=wp) :: fact1, fact2, fodd, hip, fihip, z, sth, phi
 real(kind=wp) :: phi_nv, phi_sv

 integer, parameter :: ns_max=8192
 real(kind=wp), parameter :: pi=3.141592653589793238462643383279502884197d0,&
                             twopi=6.283185307179586476925286766559005768394d0,&
                             fourpi=12.56637061435917295385057353311801153679d0,&
                             halfpi=1.570796326794896619231321691639751442099d0,&
                             twothird=0.6666666666666666666666666666666666666666d0

!-----------------------------------------------------------------------
 if (nside.lt.1 .or. nside.gt.ns_max) print*,"nside out of range"
 npix = 12*nside**2       ! total number of points
 if (ipix.lt.1 .or. ipix.gt.npix) print*, "ipix out of range"
  
! ipix1 = ipix + 1          ! in {0, npix-1}
 ipix1 = ipix              ! in {1, npix}
 nl2   = 2*nside
 nl4   = 4*nside
 ncap  = 2*nside*(nside-1) ! points in each polar cap, =0 for nside =1
 fact1 = 1.5*nside
 fact2 = 3.0*nside**2

 phi_nv = 0.0
 phi_sv = 0.0
 if (ipix1 .le. ncap) then ! North Polar cap -------------
    hip   = ipix1/2.0
    fihip = aint ( hip )
    iring = int( sqrt( max(0.0_wp, hip - sqrt(fihip)) ) ) + 1 ! counted from North pole
    iphi  = ipix1 - 2*iring*(iring - 1)

    z     =  1.0 - iring**2 / fact2
    phi   = (real(iphi) - 0.5) * PI/(2.0*iring)
 elseif (ipix1 .le. nl2*(5*nside+1)) then   ! Equatorial region ------
    ip    = ipix1 - ncap - 1
    iring = int( ip / nl4 ) + nside         ! counted from North pole
    iphi  = mod(ip,nl4) + 1

    fodd  = 0.5 * (1 + mod(iring+nside,2))  ! 1 if iring+nside is odd
                                            ! 1/2 otherwise
    z     = (nl2 - iring) / fact1
    phi   = (real(iphi) - fodd) * PI /(2.0*nside)
 else ! South Polar cap -----------------------------------
    ip    = npix - ipix1 + 1
    hip   = ip/2.0
    fihip = aint ( hip )
    iring = int( sqrt( max(0.0_wp, hip - sqrt(fihip)) ) ) + 1  ! counted from South pole
    iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1))

    z     = -1.0 + iring**2 / fact2
    phi   = (real(iphi) - 0.5) * PI/(2.0*iring)
 endif
  
 ! pixel center
 sth = sqrt((1.0-z)*(1.0+z))
 vx  = sth * cos(phi)
 vy  = sth * sin(phi)
 vz  = z

 return
 end subroutine pix2vec
end module healpix
