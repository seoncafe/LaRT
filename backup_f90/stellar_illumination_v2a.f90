module stellar_illumination_mod
  !--
  !-- 2021.09.14: use the composite method for the sampling of cos(vartheta), if par%sampling_method > 0.
  !-- 2021.09.01: Written by K.I. Seon.
  !--     random_stellar_illumination2: inspired by the method of Dongdong Yan (Yunnan Observatories, Chinese Academy of Sciences)
  !--
  public random_stellar_illumination
  public peeling_direct_stellar_illumination
  private
  interface random_stellar_illumination
     !-- 2020.09.02
     !-- '...1' version use the rejection method for (vartheta, varphi, theta, phi).
     !-- '...2' version choose vartheta and varphi using the alias method, and  use the rejection method for (theta, phi).
     !module procedure random_stellar_illumination1
     module procedure random_stellar_illumination2
  end interface random_stellar_illumination
  interface peeling_direct_stellar_illumination
     !-- '...1' version is for when a spherical medium is filled.
     !-- '...2' version is for a more general purpose.
     module procedure peeling_direct_stellar_illumination1
     !module procedure peeling_direct_stellar_illumination2
  end interface peeling_direct_stellar_illumination
contains
  !================================================
  subroutine random_stellar_illumination1(grid,photon)
  use define
  use random
  !TEST---use omp_lib
  implicit none
  type(grid_type),   intent(in)    :: grid
  type(photon_type), intent(inout) :: photon
  real(kind=wp) :: cos_ang
  real(kind=wp) :: cosvt, sinvt, vphi, cosvp, sinvp
  real(kind=wp) :: cost, sint, phi, cosp, sinp
  real(kind=wp) :: x0, y0, z0, x, y, z, kx0, ky0, kz0, rr, kr
  real(kind=wp) :: det, r_dot_k, dist
  logical,       save :: parameters_initialized = .false.
  real(kind=wp), save :: cosvt_max, cost_max, flux_fac1
  !$OMP THREADPRIVATE(parameters_initialized, cosvt_max, cost_max, flux_fac1)

  !TEST---integer, save :: my_threadid
  !TEST---!$OMP THREADPRIVATE(my_threadid)

  if (.not. parameters_initialized) then
     cosvt_max = (par%stellar_radius - par%rmax)/par%distance_star_to_planet
     cost_max  = sqrt(1.0_wp - (par%rmax/(par%distance_star_to_planet - par%stellar_radius))**2)
     flux_fac1 = (1.0_wp - cosvt_max) * (1.0_wp - cost_max)
     parameters_initialized = .true.
     !TEST---my_threadid = omp_get_thread_num()
  endif

  photon%wgt       = 1.0_wp
  photon%nrejected = 0.0_wp
  do while(.true.)
     !-- (vartheta, varphi) = angles where the photon is injected on the stellar surface (with respect to the star.)
     cosvt = (1.0_wp - cosvt_max) * rand_number() + cosvt_max
     sinvt = sqrt(1.0_wp - cosvt**2)
     vphi  = twopi * rand_number()
     cosvp = cos(vphi)
     sinvp = sin(vphi)

     !-- (x0, y0, z0) = unit vector for the starting position of the photon
     !-- (x,  y,  z)  = the starting position of the photon expressed in the grid system (planet coordinate system).
     x0  = sinvt*cosvp
     y0  = sinvt*sinvp
     z0  = cosvt
     x   = par%stellar_radius * x0
     y   = par%stellar_radius * y0
     z   = par%stellar_radius * z0 - par%distance_star_to_planet

     !-- (kx0, ky0, kz0) = unit vector connecting the planet center and the injection position on the stellar surface
     rr  = sqrt(x**2 + y**2 + z**2)
     kx0 = -x/rr
     ky0 = -y/rr
     kz0 = -z/rr

     !-- (theta, phi) = angles for the direction vector of the photon, expressed about the unit vector (kx0, ky0, kz0).
     cost     = (1.0_wp - cost_max) * rand_number() + cost_max
     sint     = sqrt(1.0_wp - cost**2)
     phi      = twopi * rand_number()
     cosp     = cos(phi)
     sinp     = sin(phi)

     !-- photon direction vector, expressed in the grid system (planet coordinate system).
     if (abs(kz0) >= 0.99999999999_wp) then
        photon%kx = sint*cosp
        photon%ky = sint*sinp
        photon%kz = cost
     else
        kr        = sqrt(kx0**2 + ky0**2)
        photon%kx = cost*kx0 + sint*(kz0*kx0*cosp - ky0*sinp)/kr
        photon%ky = cost*ky0 + sint*(kz0*ky0*cosp + kx0*sinp)/kr
        photon%kz = cost*kz0 - sint*cosp*kr
     endif

     !-- Find the location where the ray touches the sphere.
     r_dot_k = x*photon%kx + y*photon%ky + z*photon%kz
     det     = r_dot_k**2 - (rr**2 - par%rmax**2)

     !-- angle between "the unit vector of the starting-position" and "the propation vector."
     cos_ang = x0*photon%kx + y0*photon%ky + z0*photon%kz

     if (cos_ang >= 0.0_wp .and. det >= 0.0_wp) then
        dist         = -r_dot_k - sqrt(det)
        photon%x     = x + photon%kx * dist
        photon%y     = y + photon%ky * dist
        photon%z     = z + photon%kz * dist
        photon%icell = floor((photon%x-grid%xmin)/grid%dx)+1
        photon%jcell = floor((photon%y-grid%ymin)/grid%dy)+1
        photon%kcell = floor((photon%z-grid%zmin)/grid%dz)+1
        if (photon%icell == 0) photon%icell = 1
        if (photon%jcell == 0) photon%jcell = 1
        if (photon%kcell == 0) photon%kcell = 1
        if (photon%icell == grid%nx+1) photon%icell = grid%nx
        if (photon%jcell == grid%ny+1) photon%jcell = grid%ny
        if (photon%kcell == grid%nz+1) photon%kcell = grid%nz
        !--- a factor to calculate the flux (luminosity) incident onto the planet surface.
        x  = photon%x
        y  = photon%y
        z  = photon%z + par%distance_star_to_planet
        rr = sqrt(x**2 + y**2 + z**2)
        photon%flux_factor = flux_fac1 * (x*photon%kx + y*photon%ky + z*photon%kz) / rr
        !TEST---write((mpar%p_rank+1) * 1000 + my_threadid+1, '(4ES15.7)') cosvt, vphi, cost, phi
        exit
     else
        photon%nrejected = photon%nrejected + 1.0_wp
     endif
  enddo

  if (par%use_stokes) then
     !--- Set the polarization basis vectors (m) and (n) perpendicular to the propagation direction.
     photon%mx =  cost * cosp
     photon%my =  cost * sinp
     photon%mz = -sint
     photon%nx = -sinp
     photon%ny =  cosp
     photon%nz =  0.0_wp

     !--- Set the Stokes parameters (assume unpolarized light)
     photon%I = 1.0_wp
     photon%Q = 0.0_wp
     photon%U = 0.0_wp
     photon%V = 0.0_wp
  endif
  end subroutine random_stellar_illumination1
  !================================================
  subroutine random_stellar_illumination2(grid,photon)
  use define
  use random
  use mathlib
  !TEST---use omp_lib
  implicit none
  type(grid_type),   intent(in)    :: grid
  type(photon_type), intent(inout) :: photon
  real(wp) :: cos_ang
  real(wp) :: cosvt, sinvt, vphi, cosvp, sinvp
  real(wp) :: cost, sint, phi, cosp, sinp
  real(wp) :: x0, y0, z0, x, y, z, kx0, ky0, kz0, rr, kr
  real(wp) :: det, r_dot_k, dist

  integer  :: i, idx
  real(wp) :: tmp, mu1, cotc, cosp1, sinp1, phi1, mu_s, dcos_vt
  integer, parameter :: nvtheta = 1000
  logical,  save :: parameters_initialized = .false.
  real(wp), save :: r1, r2, tot_omega_S, cos_vt1, cos_vt2, cos_vt3
  real(wp), allocatable, save :: cos_vt(:), Prob(:), Palias(:), wgt(:)
  integer,  allocatable, save :: alias(:)
  !$OMP THREADPRIVATE(parameters_initialized, r1, r2, tot_omega_S, cos_vt1, cos_vt2, cos_vt3, cos_vt, Prob, Palias, wgt, alias)

  !TEST---integer, save :: my_threadid
  !TEST---!$OMP THREADPRIVATE(my_threadid)

  !--- Initialize Probability Distribution function for the polar angle on the stellar surface.
  if (.not. parameters_initialized) then
     r1      = par%stellar_radius / par%distance_star_to_planet
     r2      = par%rmax           / par%distance_star_to_planet
     cos_vt1 = r1 + r2 
     cos_vt2 = r1
     cos_vt3 = r1 - r2

     !TEST---my_threadid = omp_get_thread_num()

     if (.not.allocated(cos_vt)) allocate(cos_vt(nvtheta))
     if (.not.allocated(Prob))   allocate(Prob(nvtheta))

     dcos_vt = (1.0d0 - cos_vt3)/(nvtheta - 1.0d0)
     do i=1, nvtheta
        cos_vt(i)  = (i - 1.0d0)*dcos_vt + cos_vt3
        tmp        = r1**2 + 1.0d0 - 2.0d0*r1*cos_vt(i)
        mu1        = sqrt((tmp - r2**2)/tmp)
        if (cos_vt(i) >= cos_vt1) then
           Prob(i) = 2.0d0 * pi * (1.0d0 - mu1)
        else if (cos_vt(i) >= cos_vt2) then
           sinvt   = sqrt(1.0d0 - cos_vt(i)**2)
           cotc    = (cos_vt(i) - r1)/sinvt
           cosp1   = cotc * mu1/sqrt(1.0d0 - mu1**2)
           !-- to avoide a numerical error.
           if (cosp1 > 1.0d0) cosp1 = 1.0d0
           phi1    = acos(cosp1)
           sinp1   = sqrt(1.0d0 - cosp1**2)
           Prob(i) = 2.0d0*pi*(1.0d0-mu1) + 2.0d0*phi1*mu1 - 2.0d0*asin(sinp1/sqrt(1.0d0+cotc**2))
        else
           sinvt   = sqrt(1.0d0 - cos_vt(i)**2)
           cotc    = (cos_vt(i) - r1)/sinvt
           cosp1   = cotc * mu1/sqrt(1.0d0 - mu1**2)
           !-- to avoide a numerical error.
           if (cosp1 < -1.0d0) cosp1 = -1.0d0
           phi1    = acos(cosp1)
           sinp1   = sqrt(1.0d0 - cosp1**2)
           Prob(i) = -2.0d0*pi*mu1 + 2.0d0*phi1*mu1 + 2.0d0*asin(sinp1/sqrt(1.0d0+cotc**2))
        endif
        !TEST---if (mpar%p_rank == 0 .and. my_threadid == 1) then
        !TEST---   write(101,'(2ES15.7)') cos_vt(i), Prob(i)
        !TEST---endif
     enddo
     !-- probability distribution density function for cosvt
     tot_omega_S = sum(Prob) * dcos_vt
     Prob(:)     = Prob(:)/tot_omega_S

     !-- probability for the cosvt bins
     if (.not.allocated(Palias)) allocate(Palias(nvtheta-1))
     if (.not.allocated(alias))  allocate(alias(nvtheta-1))
     do i=1, nvtheta-1
        Palias(i) = (Prob(i)+Prob(i+1))/2.0d0
     enddo
     Palias(:) = Palias(:)/sum(Palias)

     !-- use the composite method (added on 2021.09.14).
     !-- Prob is the probability density function at a single value, but Palias is the probability integrated over a bin.
     if (par%sampling_method > 0) then
        if (.not.allocated(wgt)) allocate(wgt(nvtheta))
        Palias(:) = Palias(:) * (1.0_wp - par%f_composite) + par%f_composite / (nvtheta - 1.0_wp)
        wgt(:)    = Prob(:) / (Prob(:) * (1.0_wp - par%f_composite) + par%f_composite / (1.0_wp - cos_vt3))
        Prob(:)   = Prob(:) * (1.0_wp - par%f_composite) + par%f_composite / (1.0_wp - cos_vt3)
     endif

     call random_alias_setup(Palias, alias)
     parameters_initialized = .true.
  endif

  photon%nrejected = 0.0_wp

  !-- (vartheta, varphi) = angles where the photon is injected on the stellar surface (with respect to the star.)
  !-- assume that Prob is a linear function between cos_vt(i) and cos_vt(i+1).
  if (par%sampling_method > 0) then
     call random_alias_linear_wgt(Palias, alias, cos_vt, Prob, wgt, cosvt, photon%wgt)
  else
     !cosvt = rand_alias_linear(Palias, alias, cos_vt, Prob)
     call random_alias_linear(Palias, alias, cos_vt, Prob, cosvt)
     photon%wgt = 1.0_wp
  endif
  sinvt = sqrt(1.0_wp - cosvt**2)
  vphi  = twopi * rand_number()
  cosvp = cos(vphi)
  sinvp = sin(vphi)

  !-- (theta, phi) = angles for the direction vector of the photon, expressed about the unit vector (kx0, ky0, kz0).
  tmp = r1**2 + 1.0d0 - 2.0d0*r1*cosvt
  mu1 = sqrt((tmp - r2**2)/tmp)
  if (cosvt > cos_vt1) then
     phi  = twopi * rand_number()
     cosp = cos(phi)
     sinp = sin(phi)

     cost = (1.0_wp - mu1) * rand_number() + mu1
     sint = sqrt(1.0_wp - cost**2)
  else if (cosvt > cos_vt2) then
     cotc  = (cosvt - r1)/sinvt

     do while(.true.)
        phi  = twopi * rand_number()
        cosp = cos(phi)
        cost = (1.0d0 - mu1) *rand_number() + mu1
        sint = sqrt(1.0d0 - cost**2)
        if (cotc*cost/sint >= cosp) exit
     enddo
     sinp = sin(phi)
  else
     cotc  = (cosvt - r1)/sinvt
     cosp1 = cotc * mu1/sqrt(1.0d0 - mu1**2)
     !-- to avoide a numerical error.
     if (cosp1 < -1.0d0) cosp1 = -1.0d0
     phi1  = acos(cosp1)

     do while(.true.)
        phi   = 2.0_wp*(pi -  phi1) * rand_number() + phi1
        cosp  = cos(phi)
        cost  = (1.0d0 - mu1)*rand_number() + mu1
        sint  = sqrt(1.0d0 - cost**2)
        if (cotc*cost/sint >= cosp) exit
     enddo
     sinp  = sin(phi)
  endif

  !TEST---write((mpar%p_rank+1) * 1000 + my_threadid+1, '(5ES15.7)') cosvt, vphi, cost, phi, photon%wgt

  !-- (x0, y0, z0) = unit vector for the starting position of the photon
  !-- (x,  y,  z)  = the starting position of the photon expressed in the grid system (planet coordinate system).
  x0  = sinvt*cosvp
  y0  = sinvt*sinvp
  z0  = cosvt
  x   = par%stellar_radius * x0
  y   = par%stellar_radius * y0
  z   = par%stellar_radius * z0 - par%distance_star_to_planet

  !-- (kx0, ky0, kz0) = unit vector connecting the planet center and the injection position on the stellar surface
  rr  = sqrt(x**2 + y**2 + z**2)
  kx0 = -x/rr
  ky0 = -y/rr
  kz0 = -z/rr

  !-- photon direction vector, expressed in the grid system (planet coordinate system).
  if (abs(kz0) >= 0.99999999999_wp) then
     photon%kx = sint*cosp
     photon%ky = sint*sinp
     photon%kz = cost
  else
     kr        = sqrt(kx0**2 + ky0**2)
     photon%kx = cost*kx0 + sint*(kz0*kx0*cosp - ky0*sinp)/kr
     photon%ky = cost*ky0 + sint*(kz0*ky0*cosp + kx0*sinp)/kr
     photon%kz = cost*kz0 - sint*cosp*kr
  endif

  !-- Find the location where the ray touches the sphere.
  r_dot_k = x*photon%kx + y*photon%ky + z*photon%kz
  det     = r_dot_k**2 - (rr**2 - par%rmax**2)

  !-- angle between "the unit vector of the starting-position" and "the propation vector."
  cos_ang = x0*photon%kx + y0*photon%ky + z0*photon%kz

  if (cos_ang >= 0.0_wp .and. det >= 0.0_wp) then
     dist         = -r_dot_k - sqrt(det)
     photon%x     = x + photon%kx * dist
     photon%y     = y + photon%ky * dist
     photon%z     = z + photon%kz * dist
     photon%icell = floor((photon%x-grid%xmin)/grid%dx)+1
     photon%jcell = floor((photon%y-grid%ymin)/grid%dy)+1
     photon%kcell = floor((photon%z-grid%zmin)/grid%dz)+1
     if (photon%icell == 0) photon%icell = 1
     if (photon%jcell == 0) photon%jcell = 1
     if (photon%kcell == 0) photon%kcell = 1
     if (photon%icell == grid%nx+1) photon%icell = grid%nx
     if (photon%jcell == grid%ny+1) photon%jcell = grid%ny
     if (photon%kcell == grid%nz+1) photon%kcell = grid%nz
     !--- a factor to calculate the flux (luminosity) incident onto the planet surface.
     ! (x, y, z) = coordinates where the ray touches the planet, expressed in the star's frame.
     x  = photon%x
     y  = photon%y
     z  = photon%z + par%distance_star_to_planet
     rr = sqrt(x**2 + y**2 + z**2)
     photon%flux_factor = (x*photon%kx + y*photon%ky + z*photon%kz) / rr * (tot_omega_S/twopi) * photon%wgt
  else
     write(*,*) 'Something wrong in random_stellar_illumination...'
     stop
  endif

  if (par%use_stokes) then
     !--- Set the polarization basis vectors (m) and (n) perpendicular to the propagation direction.
     photon%mx =  cost * cosp
     photon%my =  cost * sinp
     photon%mz = -sint
     photon%nx = -sinp
     photon%ny =  cosp
     photon%nz =  0.0_wp

     !--- Set the Stokes parameters (assume unpolarized light)
     photon%I = 1.0_wp
     photon%Q = 0.0_wp
     photon%U = 0.0_wp
     photon%V = 0.0_wp
  endif
  end subroutine random_stellar_illumination2
  !================================================
  !--------------------------------------------------
  subroutine peeling_direct_stellar_illumination1(photon,grid)
  use define
  use random
  implicit none
  type(photon_type),   intent(in)    :: photon
  type(grid_type),     intent(in)    :: grid
  !-- local variables
  type(photon_type) :: pobs
  real(kind=wp) :: cosvt0,cosvt,sinvt,vphi,cosvp,sinvp
  real(kind=wp) :: r_dot_k,rr,det,dist
  real(kind=wp) :: r2,r,wgt,tau
  real(kind=wp) :: kx0,ky0,kz0,kr0,xx,yy,zz,kx,ky,kz
  real(kind=wp) :: xfreq_ref, u1, u2
  real(kind=wp) :: distance_star_to_obs2
  integer :: ix,iy,ixf
  integer :: i

  !-- take "frequency" from the input photon.
  !-- But, note that the input frequency is expressed in the local cell of the exoplanet.
  !-- Thus, transform the frequency back to the lab frame.
  pobs      = photon
  u1        = grid%vfx(photon%icell,photon%jcell,photon%kcell)*photon%kx + &
              grid%vfy(photon%icell,photon%jcell,photon%kcell)*photon%ky + &
              grid%vfz(photon%icell,photon%jcell,photon%kcell)*photon%kz
  xfreq_ref = (photon%xfreq + u1) * (grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)

  !-- reset photon%wgt.
  pobs%wgt = 1.0_wp

  !-- Location in spectral bins. (note that external radiation is emitted in a non-comoving frame).
  ixf       = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1

  !-- choose a random position on the stellar surface
  !-- assume that the distance between the star and the earch is fixed.
  cosvt0 = par%stellar_radius / par%distance
  !-- this is wrong. (2021.09.08)
  !cosvt  = (1.0_wp - cosvt0) * rand_number() + cosvt0
  cosvt  = (1.0_wp+cosvt0**2 - (1.0_wp-cosvt0)**2 * ((1.0_wp+cosvt0)/(1.0_wp-cosvt0))**rand_number()) / (2.0_wp*cosvt0)
  sinvt  = sqrt(1.0_wp - cosvt**2)
  vphi   = twopi*rand_number()
  cosvp  = cos(vphi)
  sinvp  = sin(vphi)

  do i=1,par%nobs
    !-- (kx0, ky0, kz0) = unit vector connecting the star and observer. a vector defining a lightcone toward the observer.
    kx0 = observer(i)%x
    ky0 = observer(i)%y
    kz0 = observer(i)%z + par%distance_star_to_planet
    kr0 = sqrt(kx0**2 + ky0**2 + kz0**2)
    kx0 = kx0/kr0
    ky0 = ky0/kr0
    kz0 = kz0/kr0

    !-- (xx, yy, zz) = a vector on the stellar surface, defined by (vartheta, varphi) about the star-observer sightline.
    if (abs(kz0) >= 0.99999999999_wp) then
       xx  = sinvt*cosvp
       yy  = sinvt*sinvp
       zz  = cosvt
    else
       kr0 = sqrt(kx0**2 + ky0**2)
       xx  = cosvt*kx0 + sinvt*(kz0*kx0*cosvp - ky0*sinvp)/kr0
       yy  = cosvt*ky0 + sinvt*(kz0*ky0*cosvp + kx0*sinvp)/kr0
       zz  = cosvt*kz0 - sinvt*cosvp*kr0
    endif
    !-- (xx, yy, zz) = photon location on the stellar surface, expressed in the grid system (planet coordinate system).
    xx = par%stellar_radius * xx
    yy = par%stellar_radius * yy
    zz = par%stellar_radius * zz - par%distance_star_to_planet

    pobs%kx = (observer(i)%x-xx)
    pobs%ky = (observer(i)%y-yy)
    pobs%kz = (observer(i)%z-zz)
    r2      = pobs%kx*pobs%kx + pobs%ky*pobs%ky + pobs%kz*pobs%kz
    r       = sqrt(r2)
    pobs%kx = pobs%kx/r
    pobs%ky = pobs%ky/r
    pobs%kz = pobs%kz/r

    !-- Transform the peeling-off vector to the observer's frame.
    kx = observer(i)%rmatrix(1,1)*pobs%kx + observer(i)%rmatrix(1,2)*pobs%ky + observer(i)%rmatrix(1,3)*pobs%kz
    ky = observer(i)%rmatrix(2,1)*pobs%kx + observer(i)%rmatrix(2,2)*pobs%ky + observer(i)%rmatrix(2,3)*pobs%kz
    kz = observer(i)%rmatrix(3,1)*pobs%kx + observer(i)%rmatrix(3,2)*pobs%ky + observer(i)%rmatrix(3,3)*pobs%kz

    !-- Location in TAN image.
    ix = floor(atan2(-kx,kz)*rad2deg/observer(i)%dxim+observer(i)%nxim/2.0_wp) + 1
    iy = floor(atan2(-ky,kz)*rad2deg/observer(i)%dyim+observer(i)%nyim/2.0_wp) + 1

    if (ix >= 1 .and. ix <= observer(i)%nxim .and. iy >= 1 .and. iy <= observer(i)%nyim) then
       distance_star_to_obs2 = observer(i)%x**2 + observer(i)%y**2 + (observer(i)%z + par%distance_star_to_planet)**2

       !-- Calculate direc0 (spectral) image.
       if (par%save_direc0) then
          !--- 2D image
          if (par%save_peeloff_2D) then
             !wgt = 1.0_wp/(fourpi*r2) * pobs%wgt
             wgt = 1.0_wp/(fourpi*distance_star_to_obs2) * pobs%wgt
             !$OMP ATOMIC UPDATE
             observer(i)%direc0_2D(ix,iy) = observer(i)%direc0_2D(ix,iy) + wgt
          endif

          !--- 3D spectral image
          if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
             !wgt = 1.0_wp/(fourpi*r2) * pobs%wgt
             wgt = 1.0_wp/(fourpi*distance_star_to_obs2) * pobs%wgt
             !$OMP ATOMIC UPDATE
             observer(i)%direc0(ixf,ix,iy) = observer(i)%direc0(ixf,ix,iy) + wgt
          endif
       endif

       !-- Calculate direc (spectral) image, which is attenuated by the intervening material, if the ray pass through the medium.
       !-- Find the location where the ray touches the sphere.
       r_dot_k = xx*pobs%kx + yy*pobs%ky + zz*pobs%kz
       rr      = sqrt(xx**2 + yy**2 + zz**2)
       det     = r_dot_k**2 - (rr**2 - par%rmax**2)

       !-- bug-fixed (2021.09.29).
       if (r_dot_k < 0.0_wp .and. det >= 0.0_wp) then
          dist       = -r_dot_k - sqrt(det)
          pobs%x     = xx + pobs%kx * dist
          pobs%y     = yy + pobs%ky * dist
          pobs%z     = zz + pobs%kz * dist
          pobs%icell = floor((pobs%x-grid%xmin)/grid%dx)+1
          pobs%jcell = floor((pobs%y-grid%ymin)/grid%dy)+1
          pobs%kcell = floor((pobs%z-grid%zmin)/grid%dz)+1
          if (pobs%icell == 0) pobs%icell = 1
          if (pobs%jcell == 0) pobs%jcell = 1
          if (pobs%kcell == 0) pobs%kcell = 1
          if (pobs%icell == grid%nx+1) pobs%icell = grid%nx
          if (pobs%jcell == grid%ny+1) pobs%jcell = grid%ny
          if (pobs%kcell == grid%nz+1) pobs%kcell = grid%nz

          !-- xfreq_ref is expressed in the lab frame. Here, we need to transform the frequency to the fluid rest frame.
          !-- note that pobs%icell /= photon%icell, etc.
          u2 = grid%vfx(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kx + &
               grid%vfy(pobs%icell,pobs%jcell,pobs%kcell)*pobs%ky + &
               grid%vfz(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kz
          pobs%xfreq = xfreq_ref * grid%Dfreq_ref/grid%Dfreq(pobs%icell,pobs%jcell,pobs%kcell) - u2

          call raytrace_to_edge(pobs,grid,tau)
          !wgt = exp(-tau)/(fourpi*r2) * pobs%wgt
          wgt = exp(-tau)/(fourpi*distance_star_to_obs2) * pobs%wgt
       else
          !wgt = 1.0_wp/(fourpi*r2) * pobs%wgt
          wgt = 1.0_wp/(fourpi*distance_star_to_obs2) * pobs%wgt
       endif

       !-- 2D image
       if (par%save_peeloff_2D) then
          !$OMP ATOMIC UPDATE
          observer(i)%direc_2D(ix,iy) = observer(i)%direc_2D(ix,iy) + wgt
          if (par%use_stokes) then
             !$OMP ATOMIC UPDATE
             observer(i)%I_2D(ix,iy) = observer(i)%I_2D(ix,iy) + wgt
          endif
       endif

       !-- 3D spectral image
       if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
          !$OMP ATOMIC UPDATE
          observer(i)%direc(ixf,ix,iy) = observer(i)%direc(ixf,ix,iy) + wgt
          if (par%use_stokes) then
             !$OMP ATOMIC UPDATE
             observer(i)%I(ixf,ix,iy) = observer(i)%I(ixf,ix,iy) + wgt
          endif
       endif
    endif
  enddo
  end subroutine peeling_direct_stellar_illumination1
  !--------------------------------------------------
  subroutine peeling_direct_stellar_illumination2(photon,grid)
  use define
  use random
  implicit none
  type(photon_type),   intent(in)    :: photon
  type(grid_type),     intent(in)    :: grid
  !-- local variables
  type(photon_type) :: pobs
  real(kind=wp) :: cosvt0,cosvt,sinvt,vphi,cosvp,sinvp
  real(kind=wp) :: x0,y0,z0,dist,delt(6)
  real(kind=wp) :: r2,r,wgt,tau
  real(kind=wp) :: kx0,ky0,kz0,kr0,xx,yy,zz,kx,ky,kz
  real(kind=wp) :: xfreq_ref, u1, u2
  real(kind=wp) :: distance_star_to_obs2
  integer :: ix,iy,ixf
  integer :: i, jj

  !-- take "frequency" and "weight" from the input photon.
  !-- But, note that the input frequency is expressed in the local cell of the exoplanet.
  !-- Thus, transform the frequency back to the lab frame.
  pobs      = photon
  u1        = grid%vfx(photon%icell,photon%jcell,photon%kcell)*photon%kx + &
              grid%vfy(photon%icell,photon%jcell,photon%kcell)*photon%ky + &
              grid%vfz(photon%icell,photon%jcell,photon%kcell)*photon%kz
  xfreq_ref = (photon%xfreq + u1) * (grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)

  !-- reset photon%wgt.
  pobs%wgt = 1.0_wp

  !-- Location in spectral bins. (note that external radiation is emitted in a non-comoving frame).
  ixf       = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1

  !-- choose a random position on the stellar surface
  !-- assume that the distance between the star and the earch is fixed.
  cosvt0 = par%stellar_radius / par%distance
  !-- this is wrong. (2021.09.08)
  !cosvt  = (1.0_wp - cosvt0) * rand_number() + cosvt0
  cosvt  = (1.0_wp+cosvt0**2 - (1.0_wp-cosvt0)**2 * ((1.0_wp+cosvt0)/(1.0_wp-cosvt0))**rand_number()) / (2.0_wp*cosvt0)
  sinvt  = sqrt(1.0_wp - cosvt**2)
  vphi   = 2.0_wp*pi*rand_number()
  cosvp  = cos(vphi)
  sinvp  = sin(vphi)

  do i=1,par%nobs
    !-- (kx0, ky0, kz0) = unit vector connecting the star and observer. a vector defining a lightcone toward the observer.
    kx0 = observer(i)%x
    ky0 = observer(i)%y
    kz0 = observer(i)%z + par%distance_star_to_planet
    kr0 = sqrt(kx0**2 + ky0**2 + kz0**2)
    kx0 = kx0/kr0
    ky0 = ky0/kr0
    kz0 = kz0/kr0

    !-- (xx, yy, zz) = a vector on the stellar surface, defined by (vartheta, varphi) about the star-observer sightline.
    if (abs(kz0) >= 0.99999999999_wp) then
       xx  = sinvt*cosvp
       yy  = sinvt*sinvp
       zz  = cosvt
    else
       kr0 = sqrt(kx0**2 + ky0**2)
       xx  = cosvt*kx0 + sinvt*(kz0*kx0*cosvp - ky0*sinvp)/kr0
       yy  = cosvt*ky0 + sinvt*(kz0*ky0*cosvp + kx0*sinvp)/kr0
       zz  = cosvt*kz0 - sinvt*cosvp*kr0
    endif
    !-- (xx, yy, zz) = photon location on the stellar surface, expressed in the grid system (planet coordinate system).
    xx = par%stellar_radius * xx
    yy = par%stellar_radius * yy
    zz = par%stellar_radius * zz - par%distance_star_to_planet

    pobs%kx = (observer(i)%x-xx)
    pobs%ky = (observer(i)%y-yy)
    pobs%kz = (observer(i)%z-zz)
    r2      = pobs%kx*pobs%kx + pobs%ky*pobs%ky + pobs%kz*pobs%kz
    r       = sqrt(r2)
    pobs%kx = pobs%kx/r
    pobs%ky = pobs%ky/r
    pobs%kz = pobs%kz/r

    !-- Transform the peeling-off vector to the observer's frame.
    kx = observer(i)%rmatrix(1,1)*pobs%kx + observer(i)%rmatrix(1,2)*pobs%ky + observer(i)%rmatrix(1,3)*pobs%kz
    ky = observer(i)%rmatrix(2,1)*pobs%kx + observer(i)%rmatrix(2,2)*pobs%ky + observer(i)%rmatrix(2,3)*pobs%kz
    kz = observer(i)%rmatrix(3,1)*pobs%kx + observer(i)%rmatrix(3,2)*pobs%ky + observer(i)%rmatrix(3,3)*pobs%kz

    !-- Location in TAN image.
    ix = floor(atan2(-kx,kz)*rad2deg/observer(i)%dxim+observer(i)%nxim/2.0_wp) + 1
    iy = floor(atan2(-ky,kz)*rad2deg/observer(i)%dyim+observer(i)%nyim/2.0_wp) + 1

    if (ix >= 1 .and. ix <= observer(i)%nxim .and. iy >= 1 .and. iy <= observer(i)%nyim) then
       distance_star_to_obs2 = observer(i)%x**2 + observer(i)%y**2 + (observer(i)%z + par%distance_star_to_planet)**2

       !-- Calculate direc0 (spectral) image.
       if (par%save_direc0) then
          !--- 2D image
          if (par%save_peeloff_2D) then
             !wgt = 1.0_wp/(fourpi*r2) * pobs%wgt
             wgt = 1.0_wp/(fourpi*distance_star_to_obs2) * pobs%wgt
             !$OMP ATOMIC UPDATE
             observer(i)%direc0_2D(ix,iy) = observer(i)%direc0_2D(ix,iy) + wgt
          endif

          !--- 3D spectral image
          if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
             !wgt = 1.0_wp/(fourpi*r2) * pobs%wgt
             wgt = 1.0_wp/(fourpi*distance_star_to_obs2) * pobs%wgt
             !$OMP ATOMIC UPDATE
             observer(i)%direc0(ixf,ix,iy) = observer(i)%direc0(ixf,ix,iy) + wgt
          endif
       endif

       !-- Find the closest boundary where the ray touches the grid system.
       if (pobs%kx == 0.0_wp) then
          delt(1) = hugest
          delt(2) = hugest
       else
          delt(1) = (grid%xmax-xx)/pobs%kx
          delt(2) = (grid%xmin-xx)/pobs%kx
       endif
       if (pobs%ky == 0.0_wp) then
          delt(3) = hugest
          delt(4) = hugest
       else
          delt(3) = (grid%ymax-yy)/pobs%ky
          delt(4) = (grid%ymin-yy)/pobs%ky
       endif
       if (pobs%kz == 0.0_wp) then
          delt(5) = hugest
          delt(6) = hugest
       else
          delt(5) = (grid%zmax-zz)/pobs%kz
          delt(6) = (grid%zmin-zz)/pobs%kz
       endif

       dist = hugest
       do jj=1,6
          if (delt(jj) > 0.0_wp .and. delt(jj) < hugest) then
             x0         = xx + pobs%kx * delt(jj)
             y0         = yy + pobs%ky * delt(jj)
             z0         = zz + pobs%kz * delt(jj)
             pobs%icell = floor((x0-grid%xmin)/grid%dx)+1
             pobs%jcell = floor((y0-grid%ymin)/grid%dy)+1
             pobs%kcell = floor((z0-grid%zmin)/grid%dz)+1
             if (pobs%icell >=0 .and. pobs%icell <= grid%nx+1 .and. &
                 pobs%jcell >=0 .and. pobs%jcell <= grid%ny+1 .and. &
                 pobs%kcell >=0 .and. pobs%kcell <= grid%nz+1) then
                if (delt(jj) < dist) dist = delt(jj)
             endif
          endif
       enddo
       !--

       !-- Calculate direc (spectral) image, which is attenuated by the intervening material, if the ray pass through the medium.
       if (dist < hugest) then
          pobs%x     = xx + pobs%kx * dist
          pobs%y     = yy + pobs%ky * dist
          pobs%z     = zz + pobs%kz * dist
          pobs%icell = floor((pobs%x-grid%xmin)/grid%dx)+1
          pobs%jcell = floor((pobs%y-grid%ymin)/grid%dy)+1
          pobs%kcell = floor((pobs%z-grid%zmin)/grid%dz)+1
          if (pobs%icell == 0) pobs%icell = 1
          if (pobs%jcell == 0) pobs%jcell = 1
          if (pobs%kcell == 0) pobs%kcell = 1
          if (pobs%icell == grid%nx+1) pobs%icell = grid%nx
          if (pobs%jcell == grid%ny+1) pobs%jcell = grid%ny
          if (pobs%kcell == grid%nz+1) pobs%kcell = grid%nz

          !-- xfreq_ref is expressed in the lab frame. Here, we need to transform the frequency to the fluid rest frame.
          !-- note that pobs%icell /= photon%icell, etc.
          u2 = grid%vfx(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kx + &
               grid%vfy(pobs%icell,pobs%jcell,pobs%kcell)*pobs%ky + &
               grid%vfz(pobs%icell,pobs%jcell,pobs%kcell)*pobs%kz
          pobs%xfreq = xfreq_ref * grid%Dfreq_ref/grid%Dfreq(pobs%icell,pobs%jcell,pobs%kcell) - u2

          call raytrace_to_edge(pobs,grid,tau)
          !wgt = exp(-tau)/(fourpi*r2) * pobs%wgt
          wgt = exp(-tau)/(fourpi*distance_star_to_obs2) * pobs%wgt
       else
          !wgt = 1.0_wp/(fourpi*r2) * pobs%wgt
          wgt = 1.0_wp/(fourpi*distance_star_to_obs2) * pobs%wgt
       endif

       !-- 2D image
       if (par%save_peeloff_2D) then
          !$OMP ATOMIC UPDATE
          observer(i)%direc_2D(ix,iy) = observer(i)%direc_2D(ix,iy) + wgt
          if (par%use_stokes) then
             !$OMP ATOMIC UPDATE
             observer(i)%I_2D(ix,iy) = observer(i)%I_2D(ix,iy) + wgt
          endif
       endif

       !-- 3D spectral image
       if (par%save_peeloff_3D .and. ixf >= 1 .and. ixf <= grid%nxfreq) then
          !$OMP ATOMIC UPDATE
          observer(i)%direc(ixf,ix,iy) = observer(i)%direc(ixf,ix,iy) + wgt
          if (par%use_stokes) then
             !$OMP ATOMIC UPDATE
             observer(i)%I(ixf,ix,iy) = observer(i)%I(ixf,ix,iy) + wgt
          endif
       endif
    endif
  enddo
  end subroutine peeling_direct_stellar_illumination2
  !--------------------------------------------------
end module stellar_illumination_mod
