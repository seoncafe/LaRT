module generate_photon_mod
contains
  subroutine generate_photon(grid,photon)

  use define
  use random
  use random_sersic
  use mathlib
  use stellar_illumination_mod
  use point_illumination_mod
  use octree_mod, only: amr_grid, amr_find_leaf
  !--- now, peelingoff subrotines are defined in define.f90 (2023.01.16).
  !use peelingoff_mod
  implicit none
  type(grid_type),   intent(inout) :: grid
  type(photon_type), intent(inout) :: photon
  !--- local variables
  real(kind=wp) :: sint,cost,phi,rp
  real(kind=wp) :: xfreq_lab, u1
  real(kind=wp) :: DnuHK
  real(kind=wp) :: del_xfreq, p1,p2,p3,ptot,pcum1,pcum2,xi
  integer       :: iup, idown
  integer       :: ix, idx, il
  real(kind=wp) :: Dfreq_local, voigt_a_local, vfx_local, vfy_local, vfz_local

  !=== set up photon's position vector.
  !--- grid%xmin/xrange/etc. are valid for both Cartesian and AMR (amr_sync_to_grid
  !--- copies amr_grid geometry into grid before photon generation begins).
  select case(trim(par%source_geometry))
  case ('uniform_sphere','sphere')
     rp   = (rand_number())**(1.0_wp/3.0_wp) * par%source_rmax
     cost = 2.0_wp*rand_number()-1.0_wp
     sint = sqrt(1.0_wp-cost*cost)
     phi  = twopi*rand_number()
     photon%x = rp*sint*cos(phi)
     photon%y = rp*sint*sin(phi)
     photon%z = rp*cost
     call setup_isotropic_injection(grid,photon)
  case ('uniform_cylinder','cylinder')
     rp   = sqrt(rand_number()) * par%source_rmax
     phi  = twopi*rand_number()
     photon%x = rp*cos(phi)
     photon%y = rp*sin(phi)
     photon%z = grid%zrange*rand_number() + grid%zmin
     call setup_isotropic_injection(grid,photon)
  case ('uniform')
     photon%x = grid%xrange*rand_number() + grid%xmin
     photon%y = grid%yrange*rand_number() + grid%ymin
     photon%z = grid%zrange*rand_number() + grid%zmin
     call setup_isotropic_injection(grid,photon)
  case ('uniform_xy')
     if (par%source_rmax > 0.0_wp) then
        rp   = sqrt(rand_number()) * par%source_rmax
        phi  = twopi*rand_number()
        photon%x = rp*cos(phi)
        photon%y = rp*sin(phi)
     else
        photon%x = grid%xrange*rand_number()+grid%xmin
        photon%y = grid%yrange*rand_number()+grid%ymin
     endif
     photon%z = 0.0_wp
     call setup_isotropic_injection(grid,photon)
  case ('gaussian')
     photon%x = grid%xrange*rand_number()+grid%xmin
     photon%y = grid%yrange*rand_number()+grid%ymin
     photon%z = par%source_zscale/sqrt(2.0_wp)*rand_gauss()
     call setup_isotropic_injection(grid,photon)
  case ('exponential')
     photon%x = grid%xrange*rand_number()+grid%xmin
     photon%y = grid%yrange*rand_number()+grid%ymin
     photon%z = par%source_zscale*rand_zexp(par%zmax/par%source_zscale)
     call setup_isotropic_injection(grid,photon)
  case ('exponential_sphere')
     rp   = rand_r2exp(par%source_rmax/par%source_rscale) * par%source_rscale
     cost = 2.0_wp*rand_number()-1.0_wp
     sint = sqrt(1.0_wp-cost*cost)
     phi  = twopi*rand_number()
     photon%x = rp*sint*cos(phi)
     photon%y = rp*sint*sin(phi)
     photon%z = rp*cost
     call setup_isotropic_injection(grid,photon)
  case ('exponential_cylinder')
     rp   = rand_r1exp(par%source_rmax/par%source_rscale) * par%source_rscale
     phi  = twopi*rand_number()
     photon%x = rp*cos(phi)
     photon%y = rp*sin(phi)
     if (par%source_zscale > 0.0) then
        photon%z = par%source_zscale*rand_zexp(par%zmax/par%source_zscale)
     else
        photon%z = grid%zrange*rand_number() + grid%zmin
     endif
     call setup_isotropic_injection(grid,photon)
  case ('sersic', 'ssh')
     !--- galaxy model in Song, Seon, & Hwang (2020).
     rp   = rand_sersic(par%sersic_m, par%source_rmax/par%Reff) * par%Reff
     cost = 2.0_wp*rand_number()-1.0_wp
     sint = sqrt(1.0_wp-cost*cost)
     phi  = twopi*rand_number()
     photon%x = rp*sint*cos(phi)
     photon%y = rp*sint*sin(phi)
     photon%z = rp*cost
     call setup_isotropic_injection(grid,photon)
  case ('star_file')
     idx      = rand_alias_choise(star%prob, star%alias)
     photon%x = star%x(idx)
     photon%y = star%y(idx)
     photon%z = star%z(idx)
     call setup_isotropic_injection(grid,photon)
     if (par%sampling_method > 0) photon%wgt = star%wgt(idx)
  case ('plane_illumination')
     call random_plane_illumination(grid,photon)
  case ('stellar_illumination')
     call random_stellar_illumination(grid,photon)
  case ('point_illumination')
     call random_point_illumination(grid,photon)
  case ('diffuse_emissivity')
     call setup_diffuse_emissivity(grid,photon)
  case default
     photon%x = par%xs_point
     photon%y = par%ys_point
     photon%z = par%zs_point
     call setup_isotropic_injection(grid,photon)
  endselect

  photon%nscatt_gas  = 0.0_wp
  photon%nscatt_dust = 0.0_wp
  photon%inside      = .true.
  photon%xfreq       = par%xfreq0
  photon%vfy_shear   = 0.0_wp

  !=== get local cell data (Dfreq, voigt_a, velocity)
  if (par%use_amr_grid) then
     il            = photon%icell_amr
     Dfreq_local   = amr_grid%Dfreq(il)
     voigt_a_local = amr_grid%voigt_a(il)
     vfx_local     = amr_grid%vfx(il)
     vfy_local     = amr_grid%vfy(il)
     vfz_local     = amr_grid%vfz(il)
  else
     Dfreq_local   = grid%Dfreq(photon%icell, photon%jcell, photon%kcell)
     voigt_a_local = grid%voigt_a(photon%icell, photon%jcell, photon%kcell)
     vfx_local     = grid%vfx(photon%icell, photon%jcell, photon%kcell)
     vfy_local     = grid%vfy(photon%icell, photon%jcell, photon%kcell)
     vfz_local     = grid%vfz(photon%icell, photon%jcell, photon%kcell)
  endif

  !=== set up photon's frequency.
  if (trim(par%spectral_type) /= 'continuum') then
     select case(line%line_type)
     case (2)
        !--- two upward transitions
        DnuHK = line%DnuHK_Hz / Dfreq_local
        if (rand_number() <= one_over_three) then
           photon%xfreq = photon%xfreq - DnuHK
        endif
     case (3)
        !--- two upward transitions + one downward transitions.
        !--- Select which upward transition will occur.
        p1 = line%f12(1)
        p2 = line%f12(2)
        p1 = p1/(p2 + p1)

        !--- Select an atom by which the photon is scattered.
        if (rand_number() > p1) then
           !--- scattered by an atom at 2 state
           iup = 2
           del_xfreq    = line%delE_Hz(2) / Dfreq_local
           photon%xfreq = photon%xfreq - del_xfreq
        endif
     case (4)
        !--- one upward transition + multiple downward transitions.
        idown = rand_alias_choise(line%b(1)%P_down, line%b(1)%A_down)
        if (idown /= 1) then
           del_xfreq    = line%b(1)%Elow_Hz(idown) / Dfreq_local
           photon%xfreq = photon%xfreq - del_xfreq
        endif
     case (5)
        !--- two upward transitions + multiple downward transitions.
        !--- Select which upward transition will occur.
        p1 = line%f12(1)
        p2 = line%f12(2)
        p1 = p1/(p2 + p1)

        !--- Select an atom by which the photon is scattered.
        if (rand_number() < p1) then
           !--- scattered by an atom at 1 state
           iup = 1
        else
           !--- scattered by an atom at 2 state
           iup = 2
           del_xfreq    = line%delE_Hz(2) / Dfreq_local
           photon%xfreq = photon%xfreq - del_xfreq
        endif

        !--- Select which downward transition will occur.
        if (line%b(iup)%ndown > 1) then
           idown        = rand_alias_choise(line%b(iup)%P_down, line%b(iup)%A_down)
           del_xfreq    = line%b(iup)%Elow_Hz(idown) / Dfreq_local
           photon%xfreq = photon%xfreq - del_xfreq
        else
           idown        = 1
        endif
     case (6)
        p1    = line%f12(1)
        p2    = line%f12(2)
        p3    = line%f12(3)
        ptot  = p1 + p2 + p3
        pcum1 = p1/ptot
        pcum2 = (p1 + p2)/ptot

        !--- Select an atom by which the photon is scattered.
        xi = rand_number()
        if (xi < pcum1) then
           !--- scattered by an atom at 1 state
           iup = 1
        else if (xi < pcum2) then
           !--- scattered by an atom at 2 state
           iup = 2
           del_xfreq    = line%delE_Hz(2) / Dfreq_local
           photon%xfreq = photon%xfreq - del_xfreq
        else
           !--- scattered by an atom at 3 state
           iup = 3
           del_xfreq    = line%delE_Hz(3) / Dfreq_local
           photon%xfreq = photon%xfreq - del_xfreq
        endif
     endselect
  endif

  !--- xfreq should be expressed in units of the Doppler frequency for the reference temperature (par%temperature).
  select case(trim(par%spectral_type))
  case ('continuum')
     photon%xfreq = rand_number() * (grid%xfreq_max - grid%xfreq_min) + grid%xfreq_min
     !--- xfreq should be transformed into units of local cell (2021-05-18).
     photon%xfreq = photon%xfreq / (Dfreq_local / grid%Dfreq_ref)
  case ('voigt0')
     !--- voigt0 is intended to make photons from a temperature that is independent of the cell temperature.
     !--- xfreq should be transformed into units of local cell (2019-08-18).
     photon%xfreq = photon%xfreq + rand_voigt(par%voigt_a0) * par%Dfreq0 / Dfreq_local
  case ('voigt')
     photon%xfreq = photon%xfreq + rand_voigt(voigt_a_local)
  case ('gaussian')
     !--- bug-fixed (2022.05.28)
     photon%xfreq = photon%xfreq + rand_gauss() * (par%gaussian_width_vel / (line%vtherm1 * sqrt(par%temperature)))
     photon%xfreq = photon%xfreq / (Dfreq_local / grid%Dfreq_ref)
  case ('line_prof_file')
     photon%xfreq = rand_alias_constant(line_prof%PDF, line_prof%alias, line_prof%xfreq)
     !--- xfreq should be transformed into units of local cell.
     photon%xfreq = photon%xfreq / (Dfreq_local / grid%Dfreq_ref)
  endselect

  !--- if the photon is generated in the lab frame, then the frequency should be transformed to the fluid rest frame.
  if (.not. par%comoving_source) then
     u1 = vfx_local*photon%kx + vfy_local*photon%ky + vfz_local*photon%kz
     photon%xfreq = photon%xfreq - u1
  endif

  !--- Jin is calculated in lab frame frequency (2020-10-17)
  if (par%save_Jin) then
     u1        = vfx_local*photon%kx + vfy_local*photon%ky + vfz_local*photon%kz
     xfreq_lab = (photon%xfreq + u1) * (Dfreq_local / grid%Dfreq_ref)
     ix        = floor((xfreq_lab - grid%xfreq_min)/grid%dxfreq)+1
     if (ix >= 1 .and. ix <= grid%nxfreq) then
        !$OMP ATOMIC UPDATE
        grid%Jin(ix) = grid%Jin(ix) + photon%wgt
     endif
  endif

  if (par%save_peeloff) then
     call peeling_direct(photon,grid)
  endif

  return
  end subroutine generate_photon

  !================================================
  subroutine setup_isotropic_injection(grid,photon)
  use define
  use random
  use octree_mod, only: amr_find_leaf
  implicit none
  type(grid_type),   intent(in)    :: grid
  type(photon_type), intent(inout) :: photon
  real(kind=wp) :: cost, sint, phi, cosp, sinp

  !-- photon weight
  if (trim(par%source_geometry) /= 'diffuse_emissivity') then
     photon%wgt = 1.0_wp
  endif

  if (.not. par%use_amr_grid .and. par%xyz_symmetry) then
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

  !--- Cell index
  if (par%use_amr_grid) then
     photon%icell_amr = amr_find_leaf(photon%x, photon%y, photon%z)
     photon%icell = 1
     photon%jcell = 1
     photon%kcell = 1
  else
     photon%icell = floor((photon%x-grid%xmin)/grid%dx)+1
     photon%jcell = floor((photon%y-grid%ymin)/grid%dy)+1
     photon%kcell = floor((photon%z-grid%zmin)/grid%dz)+1
     !--- This treatment is necessary for the external illuminating source (2021.04.25).
     if (photon%kx < 0.0_wp .and. photon%icell == grid%nx+1) photon%icell = grid%nx
     if (photon%ky < 0.0_wp .and. photon%jcell == grid%ny+1) photon%jcell = grid%ny
     if (photon%kz < 0.0_wp .and. photon%kcell == grid%nz+1) photon%kcell = grid%nz
     if (photon%kx > 0.0_wp .and. photon%icell < 1) photon%icell = 1
     if (photon%ky > 0.0_wp .and. photon%jcell < 1) photon%jcell = 1
     if (photon%kz > 0.0_wp .and. photon%kcell < 1) photon%kcell = 1
  endif

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
  subroutine setup_diffuse_emissivity(grid,photon)
  use define
  use random
  use octree_mod, only: amr_grid
  implicit none
  type(grid_type),   intent(in)    :: grid
  type(photon_type), intent(inout) :: photon
  real(kind=wp) :: rp, cost, sint, phi

  select case(trim(par%geometry))
  case ('plane_atmosphere')
     photon%x = grid%xrange * rand_number() + grid%xmin
     photon%y = grid%yrange * rand_number() + grid%ymin
     if (par%sampling_method > 0) then
        call random_alias_linear_wgt(emiss_prof%prob_alias, emiss_prof%alias, &
                                     emiss_prof%axis, emiss_prof%prob, emiss_prof%wgt, photon%z, photon%wgt)
     else
        call random_alias_linear(emiss_prof%prob_alias, emiss_prof%alias, &
                                 emiss_prof%axis, emiss_prof%prob, photon%z)
        photon%wgt = 1.0_wp
     endif
  case ('spherical_atmosphere', 'sphere')
     if (par%sampling_method > 0) then
        call random_alias_linear_wgt(emiss_prof%prob_alias, emiss_prof%alias, &
                                     emiss_prof%axis, emiss_prof%prob, emiss_prof%wgt, rp, photon%wgt)
     else
        call random_alias_linear(emiss_prof%prob_alias, emiss_prof%alias, &
                                 emiss_prof%axis, emiss_prof%prob, rp)
        photon%wgt = 1.0_wp
     endif
     cost     = 2.0_wp*rand_number()-1.0_wp
     sint     = sqrt(1.0_wp-cost*cost)
     phi      = twopi*rand_number()
     photon%x = rp*sint*cos(phi)
     photon%y = rp*sint*sin(phi)
     photon%z = rp*cost
  case default
     if (par%use_amr_grid) then
        if (par%sampling_method == 0) then
           call random_emiss_alias_amr(photon)
        else if (par%sampling_method == 1) then
           call random_emiss_composite_alias_amr(photon)
        else if (par%sampling_method == 2) then
           call random_emiss_naive_amr(photon)
        else
           call random_emiss_composite_amr(photon)
        endif
     else
        if (par%sampling_method == 0) then
           call random_emiss_alias(grid,photon)
        else if (par%sampling_method == 1) then
           call random_emiss_composite_alias(grid,photon)
        else if (par%sampling_method == 2) then
           call random_emiss_naive(grid,photon)
        else
           call random_emiss_composite(grid,photon)
        endif
     endif
  endselect

  call setup_isotropic_injection(grid,photon)
  end subroutine setup_diffuse_emissivity

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
  subroutine random_emiss_alias_amr(photon)
  use define
  use random
  use octree_mod, only: amr_grid
  implicit none
  type(photon_type), intent(inout) :: photon
  integer :: il, icell
  real(kind=wp) :: h

  photon%wgt       = 1.0_wp
  photon%nrejected = 0.0_wp

  il               = rand_alias_choise(amr_grid%Pem, amr_grid%alias)
  icell            = amr_grid%icell_of_leaf(il)
  h                = amr_grid%ch(icell)
  photon%icell_amr = il
  photon%x         = 2.0_wp*h*rand_number() + (amr_grid%cx(icell) - h)
  photon%y         = 2.0_wp*h*rand_number() + (amr_grid%cy(icell) - h)
  photon%z         = 2.0_wp*h*rand_number() + (amr_grid%cz(icell) - h)
  end subroutine random_emiss_alias_amr

  !================================================
  subroutine random_emiss_composite_alias_amr(photon)
  use define
  use random
  use octree_mod, only: amr_grid
  implicit none
  type(photon_type), intent(inout) :: photon
  integer :: il, icell
  real(kind=wp) :: h

  photon%nrejected = 0.0_wp

  il               = rand_alias_choise(amr_grid%Pem, amr_grid%alias)
  icell            = amr_grid%icell_of_leaf(il)
  h                = amr_grid%ch(icell)
  photon%icell_amr = il
  photon%wgt       = amr_grid%Pwgt(il)
  photon%x         = 2.0_wp*h*rand_number() + (amr_grid%cx(icell) - h)
  photon%y         = 2.0_wp*h*rand_number() + (amr_grid%cy(icell) - h)
  photon%z         = 2.0_wp*h*rand_number() + (amr_grid%cz(icell) - h)
  end subroutine random_emiss_composite_alias_amr

  !================================================
  subroutine random_emiss_naive_amr(photon)
  use define
  use random
  use octree_mod, only: amr_grid, amr_find_leaf
  implicit none
  type(photon_type), intent(inout) :: photon
  integer :: il

  photon%wgt       = 1.0_wp
  photon%nrejected = 0.0_wp

  do while(.true.)
     photon%x = amr_grid%xrange*rand_number() + amr_grid%xmin
     photon%y = amr_grid%yrange*rand_number() + amr_grid%ymin
     photon%z = amr_grid%zrange*rand_number() + amr_grid%zmin
     il = amr_find_leaf(photon%x, photon%y, photon%z)
     if (il <= 0) cycle
     photon%icell_amr = il
     if (rand_number() <= amr_grid%Pem(il)) then
        exit
     else
        photon%nrejected = photon%nrejected + 1.0_wp
     endif
  enddo
  end subroutine random_emiss_naive_amr

  !================================================
  subroutine random_emiss_composite_amr(photon)
  use define
  use random
  use octree_mod, only: amr_grid, amr_find_leaf
  implicit none
  type(photon_type), intent(inout) :: photon
  logical, save       :: setup_Qem = .false.
  real(kind=wp), save :: Qem
  integer :: il
  real(kind=wp) :: positive_volume, cell_volume

  if (.not. setup_Qem) then
     Qem = 0.0_wp
     positive_volume = 0.0_wp
     do il = 1, amr_grid%nleaf
        if (amr_grid%Pem(il) > 0.0_wp) then
           cell_volume = (2.0_wp * amr_grid%ch(amr_grid%icell_of_leaf(il)))**3
           Qem = Qem + amr_grid%Pem(il) * cell_volume
           positive_volume = positive_volume + cell_volume
        endif
     enddo
     if (positive_volume > 0.0_wp) Qem = Qem / positive_volume
     setup_Qem = .true.
  endif

  photon%nrejected = 0.0_wp
  if (rand_number() < par%f_composite) then
     do while(.true.)
        photon%x = amr_grid%xrange*rand_number() + amr_grid%xmin
        photon%y = amr_grid%yrange*rand_number() + amr_grid%ymin
        photon%z = amr_grid%zrange*rand_number() + amr_grid%zmin
        il = amr_find_leaf(photon%x, photon%y, photon%z)
        if (il > 0 .and. amr_grid%Pem(il) > 0.0_wp) exit
     enddo
  else
     do while(.true.)
        photon%x = amr_grid%xrange*rand_number() + amr_grid%xmin
        photon%y = amr_grid%yrange*rand_number() + amr_grid%ymin
        photon%z = amr_grid%zrange*rand_number() + amr_grid%zmin
        il = amr_find_leaf(photon%x, photon%y, photon%z)
        if (il <= 0) cycle
        if (rand_number() <= amr_grid%Pem(il)) then
           exit
        else
           photon%nrejected = photon%nrejected + 1.0_wp
        endif
     enddo
  endif

  photon%icell_amr = il
  photon%wgt = 1.0_wp / ((1.0_wp-par%f_composite) + par%f_composite * Qem/amr_grid%Pem(il))
  end subroutine random_emiss_composite_amr

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
end module generate_photon_mod
