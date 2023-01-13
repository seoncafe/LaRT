module photon_mod
contains
  subroutine gen_photon(grid,photon)

  use define
  use random
  use random_sersic
  use mathlib
  use stellar_illumination_mod
  !TEST---use omp_lib
  use peelingoff_mod
  implicit none
  type(grid_type),   intent(inout) :: grid
  type(photon_type), intent(inout) :: photon
  !--- local variables
  real(kind=wp) :: sint,cost,phi,rp
  real(kind=wp) :: xfreq_lab, u1
  real(kind=wp) :: DnuHK
  integer       :: ix, idx
  real(kind=wp), parameter :: one_over_three = 1.0_wp/3.0_wp

  !TEST---logical, save :: threadid_initialized = .false.
  !TEST---integer, save :: my_threadid
  !TEST---!$OMP THREADPRIVATE(threadid_initialized,my_threadid)
  !TEST---if (.not. threadid_initialized) then
  !TEST---   my_threadid          = omp_get_thread_num()
  !TEST---   threadid_initialized = .true.
  !TEST---endif

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
     call setup_isotropic_injection(grid,photon)
  else if (trim(par%source_geometry) == 'uniform_cylinder' .and. par%source_rmax > 0.0_wp) then
     rp   = sqrt(rand_number()) * par%source_rmax
     phi  = twopi*rand_number()
     photon%x = rp*cos(phi)
     photon%y = rp*sin(phi)
     photon%z = grid%zrange*rand_number() + grid%zmin
     call setup_isotropic_injection(grid,photon)
  else if (trim(par%source_geometry) == 'uniform' .and. par%source_rmax <= 0.0_wp) then
     photon%x = grid%xrange*rand_number() + grid%xmin
     photon%y = grid%yrange*rand_number() + grid%ymin
     photon%z = grid%zrange*rand_number() + grid%zmin
     call setup_isotropic_injection(grid,photon)
  else if (trim(par%source_geometry) == 'uniform_xy' .and. par%source_rmax <= 0.0_wp) then
     photon%x = grid%xrange*rand_number()+grid%xmin
     photon%y = grid%yrange*rand_number()+grid%ymin
     photon%z = 0.0_wp
     call setup_isotropic_injection(grid,photon)
  else if (trim(par%source_geometry) == 'gaussian') then
     photon%x = grid%xrange*rand_number()+grid%xmin
     photon%y = grid%yrange*rand_number()+grid%ymin
     photon%z = par%source_zscale/sqrt(2.0_wp)*rand_gauss()
     call setup_isotropic_injection(grid,photon)
  else if (trim(par%source_geometry) == 'exponential') then
     photon%x = grid%xrange*rand_number()+grid%xmin
     photon%y = grid%yrange*rand_number()+grid%ymin
     photon%z = par%source_zscale*rand_zexp(par%zmax/par%source_zscale)
     call setup_isotropic_injection(grid,photon)
  else if (trim(par%source_geometry) == 'exponential_sphere') then
     rp   = rand_r2exp(par%source_rmax/par%source_rscale) * par%source_rscale
     cost = 2.0_wp*rand_number()-1.0_wp
     sint = sqrt(1.0_wp-cost*cost)
     phi  = twopi*rand_number()
     photon%x = rp*sint*cos(phi)
     photon%y = rp*sint*sin(phi)
     photon%z = rp*cost
     call setup_isotropic_injection(grid,photon)
  else if (trim(par%source_geometry) == 'exponential_cylinder') then
     rp   = rand_r1exp(par%source_rmax/par%source_rscale) * par%source_rscale
     phi  = twopi*rand_number()
     photon%x = rp*cos(phi)
     photon%y = rp*sin(phi)
     photon%z = grid%zrange*rand_number() + grid%zmin
     call setup_isotropic_injection(grid,photon)
  else if (trim(par%source_geometry) == 'sersic' .or. trim(par%source_geometry) == 'ssh') then
     !--- galaxy model in Song, Seon, & Hwang (2020).
     rp   = rand_sersic(par%sersic_m, par%source_rmax/par%Reff) * par%Reff
     cost = 2.0_wp*rand_number()-1.0_wp
     sint = sqrt(1.0_wp-cost*cost)
     phi  = twopi*rand_number()
     photon%x = rp*sint*cos(phi)
     photon%y = rp*sint*sin(phi)
     photon%z = rp*cost
     call setup_isotropic_injection(grid,photon)
  else if (trim(par%source_geometry) == 'star_file') then
     idx      = rand_alias_choise(star%prob, star%alias)
     photon%x = star%x(idx)
     photon%y = star%y(idx)
     photon%z = star%z(idx)
     call setup_isotropic_injection(grid,photon)
     if (par%sampling_method > 0) photon%wgt = star%wgt(idx)
  else if (trim(par%source_geometry) == 'plane_illumination') then
     call random_plane_illumination(grid,photon)
  else if (trim(par%source_geometry) == 'stellar_illumination') then
     call random_stellar_illumination(grid,photon)
  else if (trim(par%source_geometry) == 'diffuse_emissivity') then
     if (trim(par%geometry) == 'plane_atmosphere') then
        photon%x = grid%xrange * rand_number() + grid%xmin
        photon%y = grid%yrange * rand_number() + grid%ymin
        if (par%sampling_method > 0) then
           call random_alias_linear_wgt(emiss_prof%prob_alias, emiss_prof%alias, &
                                        emiss_prof%axis, emiss_prof%prob, emiss_prof%wgt, photon%z, photon%wgt)
        else
           call random_alias_linear(emiss_prof%prob_alias, emiss_prof%alias, emiss_prof%axis, emiss_prof%prob, photon%z)
           photon%wgt = 1.0_wp
        endif
        call setup_isotropic_injection(grid,photon)
     else if (trim(par%geometry) == 'spherical_atmosphere') then
        if (par%sampling_method > 0) then
           call random_alias_linear_wgt(emiss_prof%prob_alias, emiss_prof%alias, &
                                        emiss_prof%axis, emiss_prof%prob, emiss_prof%wgt, rp, photon%wgt)
        else
           call random_alias_linear(emiss_prof%prob_alias, emiss_prof%alias, emiss_prof%axis, emiss_prof%prob, rp)
           photon%wgt = 1.0_wp
        endif
        cost       = 2.0_wp*rand_number()-1.0_wp
        sint       = sqrt(1.0_wp-cost*cost)
        phi        = twopi*rand_number()
        photon%x   = rp*sint*cos(phi)
        photon%y   = rp*sint*sin(phi)
        photon%z   = rp*cost
        call setup_isotropic_injection(grid,photon)
        !TEST---write((mpar%rank+1) * 1000 + my_threadid+1, '(2ES15.7)') rp, photon%wgt
     else
        !--- diffuse emission, using emissivity
        if (par%sampling_method == 0) then
           call random_emiss_alias(grid,photon)
        else if (par%sampling_method == 1) then
           call random_emiss_composite_alias(grid,photon)
        else if (par%sampling_method == 2) then
           call random_emiss_naive(grid,photon)
        else
           call random_emiss_composite(grid,photon)
        endif
        call setup_isotropic_injection(grid,photon)
     endif
  else
     photon%x = par%xs_point
     photon%y = par%ys_point
     photon%z = par%zs_point
     call setup_isotropic_injection(grid,photon)
  endif

  photon%nscatt_HI   = 0.0_wp
  photon%nscatt_dust = 0.0_wp
  photon%inside      = .true.
  photon%xfreq       = par%xfreq0

  !=== set up photon's frequency.
#ifdef FINE_STRUCTURE
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
     !--- xfreq should be transformed into units of local cell (2021-05-18).
     photon%xfreq = photon%xfreq / (grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)
  else if (trim(par%spectral_type) == 'voigt0') then
     !--- voigt0 is intended to make photons from a temperature that is independent of the cell temperature.
     !--- xfreq should be transformed into units of local cell (2019-08-18).
     photon%xfreq = photon%xfreq + &
                    rand_voigt(par%voigt_a0) * par%Dfreq0 / grid%Dfreq(photon%icell, photon%jcell, photon%kcell)
  else if (trim(par%spectral_type) == 'voigt') then
     photon%xfreq = photon%xfreq + rand_voigt(grid%voigt_a(photon%icell,photon%jcell,photon%kcell))
  else if (trim(par%spectral_type) == 'gaussian') then
     !--- bug-fixed (2022.05.28)
     photon%xfreq = photon%xfreq + rand_gauss() * (par%gaussian_width_vel / (0.12843374_wp * sqrt(par%temperature)))
     photon%xfreq = photon%xfreq / (grid%Dfreq(photon%icell,photon%jcell,photon%kcell) / grid%Dfreq_ref)
  else if (trim(par%spectral_type) == 'line_prof_file') then
     photon%xfreq = rand_alias_constant(line_prof%PDF, line_prof%alias, line_prof%xfreq)
     !--- xfreq should be transformed into units of local cell.
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
     xfreq_lab = (photon%xfreq + u1) * (grid%Dfreq(photon%icell, photon%jcell, photon%kcell) / grid%Dfreq_ref)
     ix        = floor((xfreq_lab - grid%xfreq_min)/grid%dxfreq)+1
     if (ix >= 1 .and. ix <= grid%nxfreq) then
        !$OMP ATOMIC UPDATE
        grid%Jin(ix) = grid%Jin(ix) + photon%wgt
     endif
  endif

  if (par%save_peeloff) then
     if (trim(par%source_geometry) == 'stellar_illumination') then
        call peeling_direct_stellar_illumination(photon,grid)
     else
        call peeling_direct(photon,grid)
     endif
  endif

  return
  end subroutine gen_photon

  !================================================
  subroutine setup_isotropic_injection(grid,photon)
  use define
  use random
  implicit none
  type(grid_type),   intent(in)    :: grid
  type(photon_type), intent(inout) :: photon
  real(kind=wp) :: cost, sint, phi, cosp, sinp

  !-- photon weight
  if (.not.(trim(par%source_geometry) == 'diffuse_emissivity')) then
     photon%wgt = 1.0_wp
  endif

  if (par%xyz_symmetry) then
     if (photon%x < grid%xmin) photon%x = -photon%x
     if (photon%y < grid%ymin) photon%y = -photon%y
     if (photon%z < grid%zmin) photon%z = -photon%z
  endif

  !--- isotropic emission
  cost = 2.0_wp*rand_number()-1.0_wp
  sint = sqrt(1.0_wp-cost*cost)
  phi  = twopi*rand_number()
  cosp = cos(phi)
  sinp = sin(phi)

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
  end subroutine setup_isotropic_injection

  !================================================
  !--- Alias sampling routine according to emissivity
  subroutine random_emiss_alias(grid,photon)
  use define
  use random
  use mathlib
  implicit none
  type(grid_type),   intent(in)    :: grid
  type(photon_type), intent(inout) :: photon
  integer :: idx

  !-- photon weight
  photon%wgt       = 1.0_wp
  photon%nrejected = 0.0_wp

  idx      = rand_alias_choise(grid%Pem1D, grid%alias)
  call array_3D_indices(grid%nx,grid%ny,grid%nz,idx,photon%icell,photon%jcell,photon%kcell)
  photon%x = (grid%xface(photon%icell+1) - grid%xface(photon%icell))*rand_number() + grid%xface(photon%icell)
  photon%y = (grid%yface(photon%jcell+1) - grid%yface(photon%jcell))*rand_number() + grid%yface(photon%jcell)
  photon%z = (grid%zface(photon%kcell+1) - grid%zface(photon%kcell))*rand_number() + grid%zface(photon%kcell)
  end subroutine random_emiss_alias

  !================================================
  !--- Alias sampling routine according to emissivity + composite method
  subroutine random_emiss_composite_alias(grid,photon)
  use define
  use random
  use mathlib
  implicit none
  type(grid_type),   intent(in)    :: grid
  type(photon_type), intent(inout) :: photon
  integer :: idx

  !-- photon weight
  photon%nrejected = 0.0_wp

  idx        = rand_alias_choise(grid%Pem1D, grid%alias)
  call array_3D_indices(grid%nx,grid%ny,grid%nz,idx,photon%icell,photon%jcell,photon%kcell)
  photon%wgt = grid%Pwgt(photon%icell,photon%jcell,photon%kcell)
  photon%x   = (grid%xface(photon%icell+1) - grid%xface(photon%icell))*rand_number() + grid%xface(photon%icell)
  photon%y   = (grid%yface(photon%jcell+1) - grid%yface(photon%jcell))*rand_number() + grid%yface(photon%jcell)
  photon%z   = (grid%zface(photon%kcell+1) - grid%zface(photon%kcell))*rand_number() + grid%zface(photon%kcell)
  end subroutine random_emiss_composite_alias

  !================================================
  !--- A simple rejection routine for sampling according to emissivity
  subroutine random_emiss_naive(grid,photon)
  use define
  use random
  implicit none
  type(grid_type),   intent(in)    :: grid
  type(photon_type), intent(inout) :: photon

  !-- photon weight
  photon%wgt       = 1.0_wp
  photon%nrejected = 0.0_wp

  do while(.true.)
     photon%x     = grid%xrange*rand_number() + grid%xmin
     photon%y     = grid%yrange*rand_number() + grid%ymin
     photon%z     = grid%zrange*rand_number() + grid%zmin
     photon%icell = floor((photon%x-grid%xmin)/grid%dx)+1
     photon%jcell = floor((photon%y-grid%ymin)/grid%dy)+1
     photon%kcell = floor((photon%z-grid%zmin)/grid%dz)+1
     if (rand_number() <= grid%Pem(photon%icell,photon%jcell,photon%kcell)) then
        exit
     else
        photon%nrejected = photon%nrejected + 1.0_wp
     endif
  enddo
  end subroutine random_emiss_naive

  !================================================
  !--- Composit biasing method (added on 2021.05.18)
  subroutine random_emiss_composite(grid,photon)
  use define
  use random
  use mathlib
  implicit none
  type(grid_type),   intent(in)    :: grid
  type(photon_type), intent(inout) :: photon

  logical, save       :: setup_Qem = .false.
  real(kind=wp), save :: Qem

  if (.not. setup_Qem) then
     Qem       = sum(grid%Pem) / count(grid%Pem > 0.0_wp)
     setup_Qem = .true.
  endif

  photon%nrejected = 0.0_wp
  if (rand_number() < par%f_composite) then
     do while(.true.)
        photon%x     = grid%xrange*rand_number() + grid%xmin
        photon%y     = grid%yrange*rand_number() + grid%ymin
        photon%z     = grid%zrange*rand_number() + grid%zmin
        photon%icell = floor((photon%x-grid%xmin)/grid%dx)+1
        photon%jcell = floor((photon%y-grid%ymin)/grid%dy)+1
        photon%kcell = floor((photon%z-grid%zmin)/grid%dz)+1
        if (grid%Pem(photon%icell,photon%jcell,photon%kcell) > 0.0_wp) exit
        !if (grid%Pem(photon%icell,photon%jcell,photon%kcell) > 0.0_wp) then
        !   exit
        !else
        !   photon%nrejected = photon%nrejected + 1.0_wp
        !endif
     enddo
  else
     do while(.true.)
        photon%x     = grid%xrange*rand_number() + grid%xmin
        photon%y     = grid%yrange*rand_number() + grid%ymin
        photon%z     = grid%zrange*rand_number() + grid%zmin
        photon%icell = floor((photon%x-grid%xmin)/grid%dx)+1
        photon%jcell = floor((photon%y-grid%ymin)/grid%dy)+1
        photon%kcell = floor((photon%z-grid%zmin)/grid%dz)+1
        if (rand_number() <= grid%Pem(photon%icell,photon%jcell,photon%kcell)) then
           exit
        else
           photon%nrejected = photon%nrejected + 1.0_wp
        endif
     enddo
  endif

  !-- photon weight
  photon%wgt = 1.0_wp / ((1.0_wp-par%f_composite) + par%f_composite * Qem/grid%Pem(photon%icell,photon%jcell,photon%kcell))
  end subroutine random_emiss_composite

  !================================================
  subroutine random_plane_illumination(grid,photon)
  use define
  use random
  implicit none
  type(grid_type),   intent(in)    :: grid
  type(photon_type), intent(inout) :: photon
  real(kind=wp) :: cost, sint, phi, cosp, sinp, rp

  !-- photon weight
  photon%wgt = 1.0_wp

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

  if (par%xyz_symmetry) then
     if (photon%x < grid%xmin) photon%x = -photon%x
     if (photon%y < grid%ymin) photon%y = -photon%y
     if (photon%z < grid%zmin) photon%z = -photon%z
  endif

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
  end subroutine random_plane_illumination
  !================================================
end module photon_mod
