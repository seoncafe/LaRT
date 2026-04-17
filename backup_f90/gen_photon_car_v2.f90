module photon_mod
contains
  subroutine gen_photon(grid,photon)

  use define
  use random
  use random_sersic
  use mathlib
#ifdef PEELINGOFF
  use peelingoff_mod
#endif
  implicit none
  type(grid_type),   intent(inout) :: grid
  type(photon_type), intent(inout) :: photon
  !--- local variables
  real(kind=wp) :: sint,cost,phi,phi2,sinp,cosp,rp
  real(kind=wp) :: xfreq_lab, u1
  real(kind=wp) :: Prand
  real(kind=wp) :: DnuHK
  integer       :: ix
  real(kind=wp), parameter :: one_over_three = 1.0_wp/3.0_wp

  !=== set up photon's position vector.
  !--- Central source
  if ((trim(par%source_geometry) == 'uniform' .or. trim(par%source_geometry) == 'sphere') .and. par%source_rmax > 0.0_wp) then
     rp   = (rand_number())**(1.0_wp/3.0_wp) * par%source_rmax
     cost = 2.0_wp*rand_number()-1.0_wp
     sint = sqrt(1.0_wp-cost*cost)
     phi  = twopi*rand_number()
     photon%x = rp*sint*cos(phi)
     photon%y = rp*sint*sin(phi)
     photon%z = rp*cost
  else if (trim(par%source_geometry) == 'uniform' .and. par%source_rmax <= 0.0_wp) then
     photon%x = grid%xrange*rand_number() + grid%xmin
     photon%y = grid%yrange*rand_number() + grid%ymin
     photon%z = grid%zrange*rand_number() + grid%zmin
  else if (trim(par%source_geometry) == 'uniform_xy' .and. par%source_rmax <= 0.0_wp) then
     photon%x = grid%xrange*rand_number()+grid%xmin
     photon%y = grid%yrange*rand_number()+grid%ymin
     photon%z = 0.0_wp
  else if (trim(par%source_geometry) == 'gaussian') then
     photon%x = grid%xrange*rand_number()+grid%xmin
     photon%y = grid%yrange*rand_number()+grid%ymin
     photon%z = par%source_zscale/sqrt(2.0_wp)*rand_gauss()
  else if (trim(par%source_geometry) == 'exponential') then
     photon%x = grid%xrange*rand_number()+grid%xmin
     photon%y = grid%yrange*rand_number()+grid%ymin
     photon%z = par%source_zscale*rand_zexp(par%zmax/par%source_zscale)
  else if (trim(par%source_geometry) == 'sersic' .or. trim(par%source_geometry) == 'ssh') then
     !--- galaxy model in Song, Seon, & Hwang (2020).
     rp   = rand_sersic(par%sersic_m, par%source_rmax/par%Reff) * par%Reff
     cost = 2.0_wp*rand_number()-1.0_wp
     sint = sqrt(1.0_wp-cost*cost)
     phi  = twopi*rand_number()
     photon%x = rp*sint*cos(phi)
     photon%y = rp*sint*sin(phi)
     photon%z = rp*cost
  else if (trim(par%source_geometry) == 'plane_illumination') then
     !--- plane-parallel illumination from z = zmax
     if (trim(par%geometry) == 'plane_atmosphere') then
        photon%x = 0.0_wp
        photon%y = 0.0_wp
        photon%z = par%zmax
     !--- plane-parallel illumination from z = zmin
     else if (trim(par%geometry) == 'spherical_atmosphere') then
        rp       = grid%rmax * sqrt(rand_number())
        if (par%xy_symmetry) then
           phi   = halfpi*rand_number()
        else
           phi   = twopi*rand_number()
        endif
        photon%x = rp*cos(phi)
        photon%y = rp*sin(phi)
        photon%z = grid%zmin
     endif
  else if (trim(par%source_geometry) == 'diffuse_emissivity') then
     !--- diffuse emission, using emissivity
     if (par%sampling_method == 1) then
        call random_emiss_composite(grid,photon)
     else
        call random_emiss_naive(grid,photon)
     endif
  else
     photon%x = par%xs_point
     photon%y = par%ys_point
     photon%z = par%zs_point
  endif

  if (par%xyz_symmetry) then
     if (photon%x < grid%xmin) photon%x = -photon%x
     if (photon%y < grid%ymin) photon%y = -photon%y
     if (photon%z < grid%zmin) photon%z = -photon%z
  endif

  !=== set up photon's propagation and direction vetor.
  if (trim(par%source_geometry) == 'plane_illumination') then
     !--- 1D plane-parallel illumination from z = zmax for exoplanet atmosphere
     if (trim(par%geometry) == 'plane_atmosphere') then
        cost = -1.0_wp
        sint =  0.0_wp
     else if (trim(par%geometry) == 'spherical_atmosphere') then
        cost = 1.0_wp
        sint = 0.0_wp
     endif
     cosp =  1.0_wp
     sinp =  0.0_wp
  else
     !--- isotropic emission
     cost = 2.0_wp*rand_number()-1.0_wp
     sint = sqrt(1.0_wp-cost*cost)
     phi2 = twopi*rand_number()
     cosp = cos(phi2)
     sinp = sin(phi2)
  endif

  !--- Set propagation direction vectors of the photon
  photon%kx = sint*cosp
  photon%ky = sint*sinp
  photon%kz = cost

  !--- Cell Index
  photon%icell = floor((photon%x-grid%xmin)/grid%dx)+1
  photon%jcell = floor((photon%y-grid%ymin)/grid%dy)+1
  photon%kcell = floor((photon%z-grid%zmin)/grid%dz)+1
  !--- This treatment is necessary for the external illuminating source (2021.04.25).
  if (photon%kx < 0.0_wp .and. photon%icell == grid%nx+1) photon%icell = grid%nx
  if (photon%ky < 0.0_wp .and. photon%jcell == grid%ny+1) photon%jcell = grid%ny
  if (photon%kz < 0.0_wp .and. photon%kcell == grid%nz+1) photon%kcell = grid%nz

  if (par%use_stokes) then
     !--- Set the reference normal vectors (m) and (n) perpendicular to the propagation direction.
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

  if (.not.(trim(par%source_geometry) == 'diffuse_emissivity' .and. par%sampling_method == 1)) then
     photon%wgt      = 1.0_wp
  endif

  photon%nscatt_HI   = 0.0_wp
  photon%nscatt_dust = 0.0_wp
  photon%inside      = .true.
  photon%xfreq       = par%xfreq0

#ifdef FINE_STRUCTURE
  !DnuHK = grid%DnuHK_ref_half/grid%Dfreq_ratio(photon%icell,photon%jcell,photon%kcell)
  DnuHK = grid%DnuHK_ref_half/grid%Dfreq(photon%icell,photon%jcell,photon%kcell)
  if (rand_number() > one_over_three) then
     photon%xfreq = photon%xfreq + DnuHK
  else
     photon%xfreq = photon%xfreq - DnuHK
  endif
#endif

  !--- This variable is defined to take into account the shear effect of TIGRESS simulation data.
  !--- 2017-07-14
  photon%vfy_shear = 0.0_wp

  !--- xfreq should be expressed in units of the Doppler frequency for the reference temperature (par%temperature).
  if (trim(par%spectral_type) == 'continuum') then
     photon%xfreq = rand_number() * (grid%xfreq_max - grid%xfreq_min) + grid%xfreq_min
     if (cpar%slope /= 0.0_wp) then
        photon%wgt = cpar%wgt0 + cpar%slope*(photon%xfreq - cpar%x0)
     else if (cpar%abs_depth > 0.0_wp) then
        photon%wgt = cpar%wgt0 - cpar%abs_depth*exp(-(photon%xfreq/cpar%abs_width)**2)
        if (photon%wgt < 0.0_wp) photon%wgt = 0.0_wp
     endif
     !--- xfreq should be transformed into units of local cell (2021-05-18).
     !photon%xfreq = photon%xfreq / grid%Dfreq_ratio(photon%icell,photon%jcell,photon%kcell)
     photon%xfreq = photon%xfreq / (grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)
  else if (trim(par%spectral_type) == 'voigt0') then
     !--- voigt0 is intended to make photons from a temperature gas that is independent of the cell temperature.
     !--- xfreq should be transformed into units of local cell (2019-08-18).
     photon%xfreq = photon%xfreq + &
                    rand_voigt(par%voigt_a0) * par%Dfreq0 / grid%Dfreq(photon%icell, photon%jcell, photon%kcell)
  else if (trim(par%spectral_type) == 'voigt') then
     photon%xfreq = photon%xfreq + rand_voigt(grid%voigt_a(photon%icell,photon%jcell,photon%kcell))
  else if (trim(par%spectral_type) == 'gaussian') then
     !photon%xfreq = photon%xfreq + 1.0_wp/sqrttwo * rand_gauss() * (par%gaussian_width_vel / 12.843374_wp)
     photon%xfreq = photon%xfreq + rand_gauss() * (par%gaussian_width_vel / 12.843374_wp)
  else if (trim(par%spectral_type) == 'line_prof_file') then
     Prand = rand_number()
     call interp_eq(line_prof%CDF, line_prof%xfreq, Prand, photon%xfreq)
     !--- xfreq should be transformed into units of local cell.
     !photon%xfreq = photon%xfreq / grid%Dfreq_ratio(photon%icell,photon%jcell,photon%kcell)
     photon%xfreq = photon%xfreq / (grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)
  endif

  !--- if the photon is generated in the lab frame, then the frequency should be transformed to the fluid rest frame.
  if (.not. par%comoving_source) then
     u1 = grid%vfx(photon%icell,photon%jcell,photon%kcell)*photon%kx + &
          grid%vfy(photon%icell,photon%jcell,photon%kcell)*photon%ky + &
          grid%vfz(photon%icell,photon%jcell,photon%kcell)*photon%kz
     photon%xfreq = photon%xfreq - u1
  endif

  !--- Jin is calculated in lab frame frequency (2020-10-17)
  if (par%save_Jin) then
     u1 = grid%vfx(photon%icell,photon%jcell,photon%kcell)*photon%kx + &
          grid%vfy(photon%icell,photon%jcell,photon%kcell)*photon%ky + &
          grid%vfz(photon%icell,photon%jcell,photon%kcell)*photon%kz
     !xfreq_lab = (photon%xfreq + u1) *  grid%Dfreq_ratio(photon%icell, photon%jcell, photon%kcell)
     xfreq_lab = (photon%xfreq + u1) * (grid%Dfreq(photon%icell, photon%jcell, photon%kcell) / grid%Dfreq_ref)
     ix        = floor((xfreq_lab - grid%xfreq_min)/grid%dxfreq)+1
     if (ix >= 1 .and. ix <= grid%nxfreq) then
        grid%Jin(ix) = grid%Jin(ix) + photon%wgt
     endif
  endif

#ifdef PEELINGOFF
  call peeling_direct(photon,grid)
#endif

  return
  end subroutine gen_photon
  !================================================
  !--- sampling routines according to emissivity
  subroutine random_emiss_naive(grid,photon)
  use define
  use random
  implicit none
  type(grid_type),   intent(in)    :: grid
  type(photon_type), intent(inout) :: photon
  do while(.true.)
     photon%x     = grid%xrange*rand_number() + grid%xmin
     photon%y     = grid%yrange*rand_number() + grid%ymin
     photon%z     = grid%zrange*rand_number() + grid%zmin
     photon%icell = floor((photon%x-grid%xmin)/grid%dx)+1
     photon%jcell = floor((photon%y-grid%ymin)/grid%dy)+1
     photon%kcell = floor((photon%z-grid%zmin)/grid%dz)+1
     if (rand_number() <= grid%Pem(photon%icell,photon%jcell,photon%kcell)) exit
  enddo
  photon%wgt = 1.0_wp
  end subroutine random_emiss_naive
  !--- Composit biasing method (added on 2021.05.18)
  subroutine random_emiss_composite(grid,photon)
  use define
  use random
  implicit none
  type(grid_type),   intent(in)    :: grid
  type(photon_type), intent(inout) :: photon

  logical, save       :: setup_Qem = .false.
  real(kind=wp), save :: Qem

  if (.not. setup_Qem) then
     Qem       = sum(grid%Pem) / count(grid%Pem > 0.0_wp)
     setup_Qem = .true.
  endif

  if (rand_number() < par%f_composite) then
     do while(.true.)
        photon%x     = grid%xrange*rand_number() + grid%xmin
        photon%y     = grid%yrange*rand_number() + grid%ymin
        photon%z     = grid%zrange*rand_number() + grid%zmin
        photon%icell = floor((photon%x-grid%xmin)/grid%dx)+1
        photon%jcell = floor((photon%y-grid%ymin)/grid%dy)+1
        photon%kcell = floor((photon%z-grid%zmin)/grid%dz)+1
        if (grid%Pem(photon%icell,photon%jcell,photon%kcell) > 0.0_wp) exit
     enddo
  else
     do while(.true.)
        photon%x     = grid%xrange*rand_number() + grid%xmin
        photon%y     = grid%yrange*rand_number() + grid%ymin
        photon%z     = grid%zrange*rand_number() + grid%zmin
        photon%icell = floor((photon%x-grid%xmin)/grid%dx)+1
        photon%jcell = floor((photon%y-grid%ymin)/grid%dy)+1
        photon%kcell = floor((photon%z-grid%zmin)/grid%dz)+1
        if (rand_number() <= grid%Pem(photon%icell,photon%jcell,photon%kcell)) exit
     enddo
  endif
  photon%wgt = 1.0_wp / ((1.0_wp-par%f_composite) + par%f_composite * Qem/grid%Pem(photon%icell,photon%jcell,photon%kcell))
  end subroutine random_emiss_composite
  !================================================
end module photon_mod
