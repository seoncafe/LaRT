module stellar_illumination_mod
  !--
  !-- 2021.09.30: Take into account the Limb Darkening effect. (par%stellar_limb_darkening = 0,1,2,3)
  !-- 2021.09.14: use the composite method for the sampling of cos(vartheta), if par%sampling_method > 0.
  !-- 2021.09.01: Written by K.I. Seon.
  !--     random_stellar_illumination2: inspired by the method of Dongdong Yan (Yunnan Observatories, Chinese Academy of Sciences)
  !--
  use define
  implicit none
  public random_stellar_illumination
  public peeling_direct_stellar_illumination
  private
  interface random_stellar_illumination
     !-- 2020.10.05
     !-- '...0' version : the most naive method. The results from other routines must be consistent with the results of this version.
     !-- '...1' version use the rejection method for (vartheta, varphi, theta, phi).
     !-- '...2' version choose vartheta and varphi using the alias method, and  use the rejection method for (theta, phi).
     !-- '...3' version : wgt = 1 even when the Limb Drarkening effect is taken into account.
     !module procedure random_stellar_illumination0
     !module procedure random_stellar_illumination1
     module procedure random_stellar_illumination2
     !module procedure random_stellar_illumination3
  end interface random_stellar_illumination
  interface peeling_direct_stellar_illumination
     !-- '...1' version is for when a spherical medium is filled.
     !-- '...2' version is for a more general purpose.
     module procedure peeling_direct_stellar_illumination1
     !module procedure peeling_direct_stellar_illumination2
  end interface peeling_direct_stellar_illumination

  !-- Coefficients for the Limb Darkening Function
  !-- Please update the coefficients to use a different limb darkening function.
  !-- Eddington approximation (Hubany & Mihalas, Theory of Stellar Atmospheres, p 572).
  !real(kind=wp), parameter :: limb_coeff(3) = [0.4d0, 0.6d0, 0.0d0]
  !-- The limb darkening function for the Sun at 550 nm. (Cox 2002, Allen's Astrophysical Quantities, p356).
  !real(kind=wp), parameter :: limb_coeff(3) = [0.3d0, 0.93d0, -0.23d0]
  !-- The limb darkening function for the Sun at 200 nm.
  real(kind=wp), parameter :: limb_coeff(3) = [0.55d0, 0.12d0, 0.33d0]
contains
  !================================================
  !-- Function to describe the limb darkening function, multiplied by cost (2021.09.27).
  !-- I(mu)/I(0) = c[0] + c[1]*mu + c[2]*mu^2
  !-- P(mu)      = I(mu)*mu = c[0]*mu + c[1]*mu^2 + c[2]*mu^3, ignoring a normalization constant.
  !-- mu         = cos(theta), theta = angle betwen the stellar radius vector and the line of sight.
  function  f_general_limb_darkening(cost) result(func)
  use define
  real(kind=wp), intent(in) :: cost
  real(kind=wp)       :: func
  integer             :: i
  logical,       save :: f_limb_init = .false.
  real(kind=wp), save :: f_sum
  !$OMP THREADPRIVATE(f_limb_init, f_sum)
  if (.not.f_limb_init) then
     f_sum = 0.0d0
     do i=1, size(limb_coeff)
        f_sum = f_sum + limb_coeff(i) / (i+1.0d0)
     enddo
     f_limb_init = .true.
  endif

  func = 0.0d0
  do i=1, size(limb_coeff)
     func = func + limb_coeff(i) * cost**i
  enddo
  func = func / f_sum
  end function f_general_limb_darkening
  !!================================================
  !!-- Weighting function to describe the limb darkening function, when the sampling of mu is performed by P(mu) = mu (2021.10.02).
  !!-- weight = P_general(mu)/(2.0*mu)
  !!-- mu     = cos(theta), theta = angle betwen the stellar radius vector and the line of sight.
  !function  w_general_limb_darkening(cost) result(func)
  !use define
  !real(kind=wp), intent(in) :: cost
  !real(kind=wp)       :: func
  !integer             :: i
  !logical,       save :: w_limb_init = .false.
  !real(kind=wp), save :: w_sum
  !!$OMP THREADPRIVATE(w_limb_init, w_sum)
  !if (.not.w_limb_init) then
  !   w_sum = 0.0d0
  !   do i=1, size(limb_coeff)
  !      w_sum = w_sum + limb_coeff(i) / (i+1.0d0)
  !   enddo
  !   w_limb_init = .true.
  !endif
  !func = 0.0d0
  !do i=1, size(limb_coeff)
  !   func = func + limb_coeff(i) * cost**(i-1)
  !enddo
  !func = func / w_sum / 2.0d0
  !end function w_general_limb_darkening
  !================================================
  !-- Random Number for the limb darkening function, which is described by a polynomial of mu. (2021.09.28)
  !-- I(mu)/I(0) = c[0] + c[1]*mu + c[2]*mu^2
  !-- P(mu)      = I(mu)*mu = c[0]*mu + c[1]*mu^2 + c[2]*mu^3, ignoring a normalization constant.
  !-- mu         = cos(theta), theta = angle betwen the stellar radius vector and the line of sight.
  function  rand_general_limb_darkening() result(cost)
  use define
  use random
  integer                     :: i, j
  integer,          parameter :: ncost = 201
  real(wp)                    :: cost
  logical,               save :: rand_limb_darkening_init = .false.
  real(wp), allocatable, save :: cost_limb(:), P_limb(:), Palias_limb(:)
  integer,  allocatable, save :: alias_limb(:)
  !$OMP THREADPRIVATE(rand_limb_darkening_init, cost_limb,P_limb,Palias_limb,alias_limb)
  if (.not.rand_limb_darkening_init) then
     if (.not.allocated(cost_limb))   allocate(cost_limb(ncost))
     if (.not.allocated(P_limb))      allocate(P_limb(ncost))
     if (.not.allocated(Palias_limb)) allocate(Palias_limb(ncost-1))
     if (.not.allocated(alias_limb))  allocate(alias_limb(ncost-1))
     do j=1, ncost
        cost_limb(j) = (j - 1.0_wp)/(ncost - 1.0_wp)
        P_limb(j)    = 0.0_wp
        do i=1, size(limb_coeff)
           P_limb(j) = P_limb(j) + limb_coeff(i) * cost_limb(j)**i
        enddo
     enddo
     P_limb = P_limb / (sum(P_limb) / (ncost - 1.0_wp))
     do j=1, ncost-1
        Palias_limb(j) = (P_limb(j) + P_limb(j+1))/2.0_wp
     enddo
     Palias_limb = Palias_limb / sum(Palias_limb)
     call random_alias_setup(Palias_limb, alias_limb)
     rand_limb_darkening_init = .true.
  endif

  cost = rand_alias_linear(Palias_limb, alias_limb, cost_limb, P_limb)
  end function rand_general_limb_darkening
  !================================================
  !-- Random Number of the Eddington limb darkening function for the Gray, LTE atmosphere (2021.09.30)
  !-- I(mu)/I(0) = 3/5 (mu + 2/3).
  !-- P(mu)      = I(mu)*mu, ignoring a normalization constant.
  !-- mu         = cos(theta), theta = angle betwen the stellar radius vector and the line of sight.
  function rand_eddington_limb_darkening() result(mu)
  use random
  real(wp) :: xi, Q, W, y, mu
  real(wp), parameter :: xi_c = 2.0_wp/27.0_wp
  xi = rand_number()
  Q  = 27.0_wp * xi - 1.0_wp
  if (xi <= xi_c) then
     y = cos(acos(Q)/3.0_wp)
  else
     W = (abs(Q) + sqrt(Q**2 - 1.0_wp))**(1.0_wp/3.0_wp)
     y = sign(1.0_wp,Q)/2.0_wp * (W + 1.0_wp/W)
  endif
  mu = (2.0_wp * y - 1.0_wp) / 3.0_wp
  end function rand_eddington_limb_darkening
  !================================================
  subroutine random_stellar_illumination0(grid,photon)
  use define
  use random
#ifdef TEST
  use omp_lib
#endif
  implicit none
  type(grid_type),   intent(in)    :: grid
  type(photon_type), intent(inout) :: photon
  real(kind=wp) :: cosvt, sinvt, vphi, cosvp, sinvp
  real(kind=wp) :: cost, sint, phi, cosp, sinp
  real(kind=wp) :: x, y, z, rr, kr, kx0, ky0, kz0
  real(kind=wp) :: det, r_dot_k, dist
  logical,       save :: parameters_initialized = .false.
  real(kind=wp), save :: flux_fac1
  !$OMP THREADPRIVATE(parameters_initialized, flux_fac1)

#ifdef TEST
  integer, save :: my_threadid
  !$OMP THREADPRIVATE(my_threadid)
#endif

  if (.not. parameters_initialized) then
     !-- In this method, flux_fac1 is independent of the limb darkening function.
     !flux_fac1              = 1.0_wp/2.0_wp
     flux_fac1              = 1.0_wp
     parameters_initialized = .true.
#ifdef TEST
     my_threadid = omp_get_thread_num()
     if (mpar%p_rank == 0 .and. my_threadid == 1) write(201,'(3ES15.7)') par%rmax,par%stellar_radius,par%distance_star_to_planet
#endif
  endif

  photon%wgt       = 1.0_wp
  photon%nrejected = 0.0_wp
  do while(.true.)
     !-- (vt, vphi) = angles for the central line connecting the stellar center and a point on the stellar surface.
     !-- use the following for flux_fac1 = 0.5_wp
     !cosvt = rand_number()
     !-- use the following for flux_fac1 = 1.0_wp
     cosvt = 2.0_wp * rand_number() - 1.0_wp
     sinvt = sqrt(1.0_wp - cosvt**2)
     vphi  = twopi * rand_number()
     cosvp = cos(vphi)
     sinvp = sin(vphi)

     !-- the Limb Darkening Function (2021.09.30).
     if (par%stellar_limb_darkening <= 0) then
        !-- constant "directional" flux
        cost = rand_number()
     else if (par%stellar_limb_darkening == 1) then
        !-- Lambert surface (constant intensity)
        cost = sqrt(rand_number())
     else if (par%stellar_limb_darkening == 2) then
        !-- Eddington Limb Darkening Function
        cost = rand_eddington_limb_darkening()
     else
        !-- a polynomial limb darkening function
        cost = rand_general_limb_darkening()
     endif

     sint = sqrt(1.0_wp - cost**2)
     phi  = twopi * rand_number()
     cosp = cos(phi)
     sinp = sin(phi)

     !-- rotation axis
     kx0 = sinvt * cosvp
     ky0 = sinvt * sinvp
     kz0 = cosvt

     !-- photon direction vector.
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

     !-- (x, y, z) = the starting position of the photon on the stellar surface expressed in the grid system (planet).
     x = par%stellar_radius * kx0
     y = par%stellar_radius * ky0
     z = par%stellar_radius * kz0 - par%distance_star_to_planet

     !-- Find the location where the ray meets the exosphere.
     r_dot_k = x*photon%kx + y*photon%ky + z*photon%kz
     rr      = sqrt(x**2 + y**2 + z**2)
     det     = r_dot_k**2 - (rr**2 - par%rmax**2)

     if (r_dot_k < 0.0_wp .and. det > 0.0_wp) then
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

        photon%flux_factor = flux_fac1 * photon%wgt
#ifdef TEST
        !-- Note that cosvt, vphi, cost, phi are different from those in illumination1 and 2 routines.
        write((mpar%p_rank+1)*1000 +my_threadid+1, '(11ES15.7)') cosvt, vphi, cost, phi, photon%wgt,&
                                                                 photon%x,photon%y,photon%z,photon%kx,photon%ky,photon%kz
#endif
        exit
     else
        photon%nrejected = photon%nrejected + 1.0_wp
     endif
  enddo

  if (par%use_stokes) then
     cost = photon%kz
     sint = sqrt(1.0_wp - cost**2)
     if (sint > 0.0_wp) then
        cosp = photon%kx / sint
        sinp = photon%ky / sint
     else
        cosp = 1.0_wp
        sinp = 0.0_wp
     endif

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
  end subroutine random_stellar_illumination0
  !================================================
  subroutine random_stellar_illumination1(grid,photon)
  use define
  use random
#ifdef TEST
  use omp_lib
#endif
  implicit none
  type(grid_type),   intent(in)    :: grid
  type(photon_type), intent(inout) :: photon
  integer       :: i
  real(kind=wp) :: cos_ang
  real(kind=wp) :: cosvt, sinvt, vphi, cosvp, sinvp
  real(kind=wp) :: cost, sint, phi, cosp, sinp
  real(kind=wp) :: x0, y0, z0, x, y, z, kx0, ky0, kz0, rr, kr
  real(kind=wp) :: det, r_dot_k, dist
  logical,       save :: parameters_initialized = .false.
  real(kind=wp), save :: cosvt_max, cost_max, flux_fac1
  !$OMP THREADPRIVATE(parameters_initialized, cosvt_max, cost_max, flux_fac1)
#ifdef USE_COSTHETA
  real(kind=wp) :: fsum1, fsum2
#endif

#ifdef TEST
  integer, save :: my_threadid
  !$OMP THREADPRIVATE(my_threadid)
#endif

  if (.not. parameters_initialized) then
     cosvt_max = (par%stellar_radius - par%rmax)/par%distance_star_to_planet
     cost_max  = sqrt(1.0_wp - (par%rmax/(par%distance_star_to_planet - par%stellar_radius))**2)

     flux_fac1 = (1.0_wp - cosvt_max) * (1.0_wp - cost_max) / 2.0_wp
#ifdef USE_COSTHETA
     !-- the default is not to use cos(theta) factor for flux but to use the number of photons (2021.10.05/10.07)
     !-- flux_fac1 depends on the limb darkening function if the flux is used to calculate luminosity.
     if (par%stellar_limb_darkening <= 0) then
        flux_fac1 = flux_fac1 * 2.0_wp
     else if (par%stellar_limb_darkening == 1) then
        flux_fac1 = flux_fac1 * 3.0_wp / 2.0_wp
     else if (par%stellar_limb_darkening == 2) then
        flux_fac1 = flux_fac1 * 24.0_wp / 17.0_wp
     else
        fsum1 = 0.0_wp
        fsum2 = 0.0_wp
        do i=1, size(limb_coeff)
           fsum1 = fsum1 + limb_coeff(i) / (i+1.0_wp)
           fsum2 = fsum2 + limb_coeff(i) / (i+2.0_wp)
        enddo
        flux_fac1 = flux_fac1 * (fsum1 / fsum2)
     endif
#endif
     parameters_initialized = .true.
#ifdef TEST
     my_threadid = omp_get_thread_num()
     if (mpar%p_rank == 0 .and. my_threadid == 1) write(201,'(3ES15.7)') par%rmax,par%stellar_radius,par%distance_star_to_planet
#endif
  endif

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

        if (par%stellar_limb_darkening <= 0) then
           !-- constant "directional" flux
           photon%wgt = 1.0_wp
        else if (par%stellar_limb_darkening == 1) then
           !-- Lambertian surface (isotropic intensity, flat disk when viewed by an observer)
           photon%wgt = 2.0d0 * cos_ang
        else if (par%stellar_limb_darkening == 2) then
           !-- Eddington Limb Darkening
           photon%wgt = cos_ang*(1.5d0*cos_ang + 1.0d0)
        else
           !-- Limb Darkening Function defined by a polynomial.
           photon%wgt = f_general_limb_darkening(cos_ang)
        endif

        photon%flux_factor = flux_fac1 * photon%wgt
#ifdef USE_COSTHETA
        !-- the default is not to use cos(theta) factor for flux but to use the number of photons (2021.10.05/10.07)
        photon%flux_factor = photon%flux_factor * cos_ang
#endif

#ifdef TEST
        write((mpar%p_rank+1)*1000 +my_threadid+1, '(11ES15.7)') cosvt, vphi, cost, phi, photon%wgt,&
                                                                 photon%x,photon%y,photon%z,photon%kx,photon%ky,photon%kz
#endif
        exit
     else
        photon%nrejected = photon%nrejected + 1.0_wp
     endif
  enddo

  if (par%use_stokes) then
     cost = photon%kz
     sint = sqrt(1.0_wp - cost**2)
     if (sint > 0.0_wp) then
        cosp = photon%kx / sint
        sinp = photon%ky / sint
     else
        cosp = 1.0_wp
        sinp = 0.0_wp
     endif

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
#ifdef TEST
  use omp_lib
#endif
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
  real(wp), save :: r1, r2, tot_omega_S, flux_fac1, cos_vt1, cos_vt2, cos_vt3
  real(wp), allocatable, save :: cos_vt(:), Prob(:), Palias(:), wgt(:)
  integer,  allocatable, save :: alias(:)
  !$OMP THREADPRIVATE(parameters_initialized, r1,r2, flux_fac1, cos_vt1,cos_vt2,cos_vt3, cos_vt, Prob, Palias, wgt, alias)
#ifdef USE_COSTHETA
  real(kind=wp) :: fsum1, fsum2
#endif

#ifdef TEST
  integer, save :: my_threadid
  !$OMP THREADPRIVATE(my_threadid)
#endif

  !--- Initialize Probability Distribution function for the polar angle on the stellar surface.
  if (.not. parameters_initialized) then
     r1      = par%stellar_radius / par%distance_star_to_planet
     r2      = par%rmax           / par%distance_star_to_planet
     cos_vt1 = r1 + r2 
     cos_vt2 = r1
     cos_vt3 = r1 - r2
#ifdef TEST
     my_threadid = omp_get_thread_num()
     if (mpar%p_rank == 0 .and. my_threadid == 1) write(201,'(3ES15.7)') par%rmax,par%stellar_radius,par%distance_star_to_planet
#endif

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
!#ifdef TEST
!        if (mpar%p_rank == 0 .and. my_threadid == 1) write(101,'(2ES15.7)') cos_vt(i), Prob(i)
!#endif
     enddo
     !-- probability distribution density function for cosvt
     tot_omega_S = sum(Prob) * dcos_vt
     Prob(:)     = Prob(:)/tot_omega_S

     flux_fac1 = tot_omega_S / fourpi
#ifdef USE_COSTHETA
     !-- the default is not to use cos(theta) factor for flux but to use the number of photons (2021.10.05/10.07)
     !-- flux_fac1 depends on the limb darkening function if the flux is used to calculate luminosity.
     if (par%stellar_limb_darkening <= 0) then
        flux_fac1 = flux_fac1 * 2.0_wp
     else if (par%stellar_limb_darkening == 1) then
        flux_fac1 = flux_fac1 * 3.0_wp / 2.0_wp
     else if (par%stellar_limb_darkening == 2) then
        flux_fac1 = flux_fac1 * 24.0_wp / 17.0_wp
     else
        fsum1 = 0.0_wp
        fsum2 = 0.0_wp
        do i=1, size(limb_coeff)
           fsum1 = fsum1 + limb_coeff(i) / (i+1.0_wp)
           fsum2 = fsum2 + limb_coeff(i) / (i+2.0_wp)
        enddo
        flux_fac1 = flux_fac1 * (fsum1 / fsum2)
     endif
#endif

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
     cosvt = rand_alias_linear(Palias, alias, cos_vt, Prob)
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

     !-- No weight change is required for the Isotropic "directional" flux case (par%stellar_limb_darkening <= 0)
     if (par%stellar_limb_darkening == 1) then
        !-- Lambertian surface (isotropic intensity, flat disk when viewed by an observer)
        photon%wgt = photon%wgt * 2.0d0 * cos_ang
     else if (par%stellar_limb_darkening == 2) then
        !-- Eddington Limb Darkening
        photon%wgt = photon%wgt * cos_ang*(1.5d0*cos_ang + 1.0d0)
     else if (par%stellar_limb_darkening > 2) then
        !-- Limb Darkening Function defined by a polynomial.
        photon%wgt = photon%wgt * f_general_limb_darkening(cos_ang)
     endif

     !-- the default is not to use cos(theta) factor for flux but to use the number of photons (2021.10.05/10.07)
     photon%flux_factor = flux_fac1 * photon%wgt
#ifdef USE_COSTHETA
     photon%flux_factor = photon%flux_factor * cos_ang
#endif

#ifdef TEST
     write((mpar%p_rank+1)*1000 +my_threadid+1, '(11ES15.7)') cosvt, vphi, cost, phi, photon%wgt,&
                                                              photon%x,photon%y,photon%z,photon%kx,photon%ky,photon%kz
#endif
  else
     write(*,*) 'Something wrong in random_stellar_illumination...'
     stop
  endif

  if (par%use_stokes) then
     cost = photon%kz
     sint = sqrt(1.0_wp - cost**2)
     if (sint > 0.0_wp) then
        cosp = photon%kx / sint
        sinp = photon%ky / sint
     else
        cosp = 1.0_wp
        sinp = 0.0_wp
     endif

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
  subroutine random_stellar_illumination3(grid,photon)
  use define
  use random
#ifdef TEST
  use omp_lib
#endif
  implicit none
  type(grid_type),   intent(in)    :: grid
  type(photon_type), intent(inout) :: photon
  real(kind=wp) :: sin_alpha
  real(kind=wp) :: cosvt_c, sinvt_c, vphi_c, cosvp_c, sinvp_c
  real(kind=wp) :: kx_c, ky_c, kz_c, x_c, y_c, z_c
  real(kind=wp) :: cosvt, sinvt, vphi, cosvp, sinvp
  real(kind=wp) :: cost, sint, phi, cosp, sinp
  real(kind=wp) :: x0, y0, z0, x, y, z, rr, kr
  real(kind=wp) :: det, r_dot_k, dist
  real(kind=wp) :: r_star, r_atmo, cos_alpha, tan_omega
  logical,       save :: parameters_initialized = .false.
  real(kind=wp), save :: distance_o, sin_omega, cosvt0, flux_fac1
  !$OMP THREADPRIVATE(parameters_initialized, distance_o, sin_omega, cosvt0, flux_fac1)

#ifdef TEST
  integer, save :: my_threadid
  !$OMP THREADPRIVATE(my_threadid)
#endif

  if (.not. parameters_initialized) then
     !-- initial version
     distance_o = par%distance_star_to_planet
     sin_omega  = par%stellar_radius/distance_o
     sin_alpha  = (par%stellar_radius + par%rmax)/par%distance_star_to_planet
     cosvt0     = sqrt(1.0_wp - sin_omega**2)*sqrt(1.0_wp - sin_alpha**2) + sin_omega * sin_alpha

     !-- optimal radius (?) of a virtual sphere to minimize the rejection rate (2021.10.10).
     !r_star     = par%stellar_radius / par%distance_star_to_planet
     !r_atmo     = par%rmax           / par%distance_star_to_planet
     !sin_alpha  = r_star + r_atmo
     !cos_alpha  = sqrt(1.0_wp - sin_alpha**2)
     !tan_omega  = (-cos_alpha + sqrt(cos_alpha**2 + 4.0*r_star*r_atmo)) / (2.0_wp*r_atmo)
     !sin_omega  = tan_omega/sqrt(1.0_wp+tan_omega**2)
     !distance_o = r_star / sin_omega * par%distance_star_to_planet
     !cosvt0     = sqrt(1.0_wp - sin_omega**2)*sqrt(1.0_wp - sin_alpha**2) + sin_omega * sin_alpha

     !-- flux_fac1 is independent of the limb darkening function if the number of photons are used to estimated
     !-- the fraction of luminosity incident on the exosphere.
     flux_fac1  = (1.0_wp - cosvt0) / 2.0_wp
     parameters_initialized = .true.
#ifdef TEST
     my_threadid = omp_get_thread_num()
     if (mpar%p_rank == 0 .and. my_threadid == 1) write(201,'(3ES15.7)') par%rmax,par%stellar_radius,par%distance_star_to_planet
#endif
  endif

  photon%wgt       = 1.0_wp
  photon%nrejected = 0.0_wp
  do while(.true.)
     !-- (vt_c, vphi_c) = angles for the central line connecting the stellar center and a point on the spherical surface a distanct_o
     cosvt_c = (1.0_wp - cosvt0) * rand_number() + cosvt0
     sinvt_c = sqrt(1.0_wp - cosvt_c**2)
     vphi_c  = twopi * rand_number()
     cosvp_c = cos(vphi_c)
     sinvp_c = sin(vphi_c)

     !-- unit vector for the central line.
     kx_c = sinvt_c*cosvp_c
     ky_c = sinvt_c*sinvp_c
     kz_c = cosvt_c

     !--
     x_c = distance_o * kx_c
     y_c = distance_o * ky_c
     z_c = distance_o * kz_c - par%distance_star_to_planet

     !-- the Limb Darkening Function (2021.09.30).
     if (par%stellar_limb_darkening <= 0) then
        !-- constant "directional" flux
        cost = rand_number()
     else if (par%stellar_limb_darkening == 1) then
        !-- Lambertian surface (isotropic intensity, flat disk when viewed by an observer)
        cost = sqrt(rand_number())
     else if (par%stellar_limb_darkening == 2) then
        !-- Eddington Limb Darkening Function
        cost = rand_eddington_limb_darkening()
     else
        !-- A polynomial Limb Darkening Function
        cost = rand_general_limb_darkening()
     endif

     !-- angles about the central line.
     cosvt = cost*sqrt(1.0_wp - sin_omega**2 + (sin_omega*cost)**2) + sin_omega*(1.0_wp - cost**2)
     sinvt = sqrt(1.0_wp - cosvt**2)
     vphi  = twopi*rand_number()
     cosvp = cos(vphi)
     sinvp = sin(vphi)

     !-- (x0, y0, z0) = unit vector for the position on the stellar surface.
     if (abs(kz_c) >= 0.99999999999_wp) then
        x0 = sinvt*cosvp
        y0 = sinvt*sinvp
        z0 = cosvt
     else
        kr = sqrt(kx_c**2 + ky_c**2)
        x0 = cosvt*kx_c + sinvt*(kz_c*kx_c*cosvp - ky_c*sinvp)/kr
        y0 = cosvt*ky_c + sinvt*(kz_c*ky_c*cosvp + kx_c*sinvp)/kr
        z0 = cosvt*kz_c - sinvt*cosvp*kr
     endif

     !-- (x,  y,  z) = the starting position of the photon on the stellar surface expressed in the grid system (planet).
     x = par%stellar_radius * x0
     y = par%stellar_radius * y0
     z = par%stellar_radius * z0 - par%distance_star_to_planet

     !--  photon direction vector
     photon%kx = x_c - x
     photon%ky = y_c - y
     photon%kz = z_c - z
     kr        = sqrt(photon%kx**2 + photon%ky**2 + photon%kz**2)
     photon%kx = photon%kx / kr
     photon%ky = photon%ky / kr
     photon%kz = photon%kz / kr

     !-- Find the location where the ray meets the exosphere.
     r_dot_k = x*photon%kx + y*photon%ky + z*photon%kz
     rr      = sqrt(x**2 + y**2 + z**2)
     det     = r_dot_k**2 - (rr**2 - par%rmax**2)

     if (det >= 0.0_wp) then
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

        photon%flux_factor = flux_fac1 * photon%wgt
#ifdef TEST
        !-- Note that cosvt, vphi, cost, phi are different from those in illumination1 and 2 routines.
        write((mpar%p_rank+1)*1000 +my_threadid+1, '(11ES15.7)') cosvt, vphi, cost, phi, photon%wgt,&
                                                                 photon%x,photon%y,photon%z,photon%kx,photon%ky,photon%kz
#endif
        exit
     else
        photon%nrejected = photon%nrejected + 1.0_wp
     endif
  enddo

  if (par%use_stokes) then
     cost = photon%kz
     sint = sqrt(1.0_wp - cost**2)
     if (sint > 0.0_wp) then
        cosp = photon%kx / sint
        sinp = photon%ky / sint
     else
        cosp = 1.0_wp
        sinp = 0.0_wp
     endif

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
  end subroutine random_stellar_illumination3
  !================================================
  subroutine peeling_direct_stellar_illumination1(photon,grid)
  use define
  use random
  implicit none
  type(photon_type),   intent(in)    :: photon
  type(grid_type),     intent(in)    :: grid
  !-- local variables
  type(photon_type) :: pobs
  real(kind=wp) :: cost
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
  ixf      = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1

  !-- choose a random position on the stellar surface
  !-- assume that the distance between the star and the earch is fixed.
  cosvt0 = par%stellar_radius / par%distance

  !-- the Limb Darkening Function (2021.09.30).
  if (par%stellar_limb_darkening <= 0) then
     !-- constant "directional" flux
     cost = rand_number()
  else if (par%stellar_limb_darkening == 1) then
     !-- Lambertian surface (isotropic intensity, flat disk when viewed by an observer)
     cost = sqrt(rand_number())
     !!-- test2
     !!cost = rand_number()
     !!pobs%wgt = 2.0d0 * cost
  else if (par%stellar_limb_darkening == 2) then
     !-- Eddington Limb Darkening Function
     cost = rand_eddington_limb_darkening()
     !!-- test1
     !!cost = sqrt(rand_number())
     !!pobs%wgt = 0.75d0*cost + 0.5d0
     !!-- test2
     !!cost = rand_number()
     !!pobs%wgt = cost*(1.5d0*cost + 1.0d0)
  else
     !-- A polynomial Limb Darkening Function
     cost = rand_general_limb_darkening()
     !!-- test1
     !!cost = sqrt(rand_number())
     !!pobs%wgt = w_general_limb_darkening(cost)
     !!-- test2
     !!cost = rand_number()
     !!pobs%wgt = f_general_limb_darkening(cost)
  endif

  cosvt = cost*sqrt(1.0_wp - cosvt0**2 + (cosvt0*cost)**2) + cosvt0*(1.0_wp - cost**2)
  sinvt = sqrt(1.0_wp - cosvt**2)
  vphi  = twopi*rand_number()
  cosvp = cos(vphi)
  sinvp = sin(vphi)
  !write(*,'(3ES15.7)') cosvt, cosvt0

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
  !================================================
  subroutine peeling_direct_stellar_illumination2(photon,grid)
  use define
  use random
  implicit none
  type(photon_type),   intent(in)    :: photon
  type(grid_type),     intent(in)    :: grid
  !-- local variables
  type(photon_type) :: pobs
  real(kind=wp) :: cost
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
  ixf      = floor((xfreq_ref - grid%xfreq_min)/grid%dxfreq)+1

  !-- choose a random position on the stellar surface
  !-- assume that the distance between the star and the earch is fixed.
  cosvt0 = par%stellar_radius / par%distance

  !-- the Limb Darkening Function (2021.09.30).
  if (par%stellar_limb_darkening <= 0) then
     !-- constant "directional" flux
     cost = rand_number()
  else if (par%stellar_limb_darkening == 1) then
     !-- Lambertian surface (isotropic intensity, flat disk when viewed by an observer)
     cost = sqrt(rand_number())
  else if (par%stellar_limb_darkening == 2) then
     !-- Eddington Limb Darkening Function
     cost = rand_eddington_limb_darkening()
  else
     !-- A polynomial Limb Darkening Function
     cost = rand_general_limb_darkening()
  endif

  cosvt = cost*sqrt(1.0_wp - cosvt0**2 + (cosvt0*cost)**2) + cosvt0*(1.0_wp - cost**2)
  sinvt = sqrt(1.0_wp - cosvt**2)
  vphi  = 2.0_wp*pi*rand_number()
  cosvp = cos(vphi)
  sinvp = sin(vphi)

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
  !================================================
end module stellar_illumination_mod
