module grid_mod
  use memory_mod
  use read_mod
  use read_atmosphere_mod
  use fitsio_mod
  use utility
  implicit none
contains
  !----------------------------------------
  subroutine grid_create(grid)
  use define
  use voigt_mod
  use mpi
  implicit none
  type(grid_type), intent(inout) :: grid

  !--- local variables
  integer :: i,j,k
  real(kind=wp) :: x,y,z,opac,taupole,taupole_dust, tauhomo,tauhomo_dust, N_HIpole,N_HIhomo
  real(kind=wp) :: dens, Vscale
  real(kind=wp) :: Pem_max
  real(kind=wp) :: opac_norm, opac_avg, opac_length
  real(kind=wp) :: radial_cell_width, density_rmin, density_rmax
  real(kind=wp) :: atau0, atau0_cell, atau1, xi, chi
  real(kind=wp) :: k_ionize, k_recomb, T4
  real(kind=wp), allocatable :: Temp(:,:,:), opacity(:,:,:), vtherm(:,:,:)
  real(kind=wp), allocatable :: xx(:), yy(:), zz(:), rr(:,:,:)
  integer :: nxcen,nycen,nzcen

  integer :: ierr
  integer :: unit,status

  !--- grid's cell faces.
  grid%nx   = par%nx
  grid%ny   = par%ny
  grid%nz   = par%nz
  grid%xmax = par%xmax
  grid%ymax = par%ymax
  grid%zmax = par%zmax
  grid%rmin = par%rmin
  grid%rmax = par%rmax

  ! grid%i0, grid%j0, grid%k0 = cell indices when the photon met a boundary and refelected into the grid system.
  ! center of grid system is always located at (x,y,z) = (0,0,0).
  if (par%xyz_symmetry) then
     if ((par%nx/2)*2 == par%nx) then
        grid%dx   = grid%xmax/grid%nx
        grid%xmin = 0.0_wp
        grid%i0   = 1
     else
        grid%dx   = grid%xmax/(grid%nx-0.5_wp)
        grid%xmin = -grid%dx/2.0_wp
        grid%i0   = 2
     endif
     if ((par%ny/2)*2 == par%ny) then
        grid%dy   = grid%ymax/grid%ny
        grid%ymin = 0.0_wp
        grid%j0   = 1
     else
        grid%dy   = grid%ymax/(grid%ny-0.5_wp)
        grid%ymin = -grid%dy/2.0_wp
        grid%j0   = 2
     endif
     if ((par%nz/2)*2 == par%nz) then
        grid%dz   = grid%zmax/grid%nz
        grid%zmin = 0.0_wp
        grid%k0   = 1
     else
        grid%dz   = grid%zmax/(grid%nz-0.5_wp)
        grid%zmin = -grid%dz/2.0_wp
        grid%k0   = 2
     endif
  else if (par%xy_symmetry) then
     if ((par%nx/2)*2 == par%nx) then
        grid%dx   = grid%xmax/grid%nx
        grid%xmin = 0.0_wp
        grid%i0   = 1
     else
        grid%dx   = grid%xmax/(grid%nx-0.5_wp)
        grid%xmin = -grid%dx/2.0_wp
        grid%i0   = 2
     endif
     if ((par%ny/2)*2 == par%ny) then
        grid%dy   = grid%ymax/grid%ny
        grid%ymin = 0.0_wp
        grid%j0   = 1
     else
        grid%dy   = grid%ymax/(grid%ny-0.5_wp)
        grid%ymin = -grid%dy/2.0_wp
        grid%j0   = 2
     endif
     grid%dz   = 2.0_wp*par%zmax/grid%nz
     grid%zmin = -grid%zmax
     grid%k0   = 0
  else if (par%z_symmetry) then
     grid%dx   = 2.0_wp*par%xmax/par%nx
     grid%dy   = 2.0_wp*par%ymax/par%ny
     grid%xmin = -grid%xmax
     grid%ymin = -grid%ymax
     grid%i0   = 0
     grid%j0   = 0
     if ((par%nz/2)*2 == par%nz) then
        grid%dz   = grid%zmax/par%nz
        grid%zmin = 0.0_wp
        grid%k0   = 1
     else
        grid%dz   = grid%zmax/(par%nz-0.5_wp)
        grid%zmin = -grid%dz/2.0_wp
        grid%k0   = 2
     endif
  else if (trim(par%geometry) == 'plane_atmosphere') then
     !--- 1D plane-parallel, exoplanet atmosphere (2021.04.21)
     grid%dx   = 2.0_wp*par%xmax/par%nx
     grid%dy   = 2.0_wp*par%ymax/par%ny
     grid%xmin = -grid%xmax
     grid%ymin = -grid%ymax
     grid%i0   = 0
     grid%j0   = 0
     !---
     if (is_finite(par%zmin)) then
        grid%zmin = par%zmin
     else
        grid%zmin = 0.0_wp
     endif
     grid%dz   = (grid%zmax-grid%zmin)/par%nz
     grid%k0   = 0
  else
     grid%dx   = 2.0_wp*par%xmax/par%nx
     grid%dy   = 2.0_wp*par%ymax/par%ny
     grid%dz   = 2.0_wp*par%zmax/par%nz
     grid%xmin = -grid%xmax
     grid%ymin = -grid%ymax
     grid%zmin = -grid%zmax
     grid%i0   = 0
     grid%j0   = 0
     grid%k0   = 0
  endif

  grid%xrange = grid%xmax - grid%xmin
  grid%yrange = grid%ymax - grid%ymin
  grid%zrange = grid%zmax - grid%zmin

  !--- allocate shared memories for xface, yface, zface
  call create_shared_mem(grid%xface, [grid%nx+1])
  call create_shared_mem(grid%yface, [grid%ny+1])
  call create_shared_mem(grid%zface, [grid%nz+1])
  if (mpar%h_rank == 0) then
     grid%xface(:) = [ ((i-1)*grid%dx + grid%xmin, i=1,grid%nx+1) ]
     grid%yface(:) = [ ((j-1)*grid%dy + grid%ymin, j=1,grid%ny+1) ]
     grid%zface(:) = [ ((k-1)*grid%dz + grid%zmin, k=1,grid%nz+1) ]
  endif

  !--- velocity field
  call create_shared_mem(grid%vfx, [grid%nx,grid%ny,grid%nz])
  call create_shared_mem(grid%vfy, [grid%nx,grid%ny,grid%nz])
  call create_shared_mem(grid%vfz, [grid%nx,grid%ny,grid%nz])

  !--- Dfreq_ratio = Ratio of local Doppler frequency to reference Doppler frequency.
  !--- Dfreq       = local Doppler frequency.
  !--- voigt_a     = A21/4pi/Dfreq (Width of Voigt profile)
  !--- rhokap      = neutral hydrogen density (rho) x kappa (absorption coefficient, per length)
  call create_shared_mem(grid%Dfreq_ratio, [grid%nx,grid%ny,grid%nz])
  call create_shared_mem(grid%Dfreq,       [grid%nx,grid%ny,grid%nz])
  call create_shared_mem(grid%voigt_a,     [grid%nx,grid%ny,grid%nz])
  call create_shared_mem(grid%rhokap,      [grid%nx,grid%ny,grid%nz])
  if (mpar%h_rank == 0) then
     grid%Dfreq_ratio(:,:,:) = 1.0_wp
  endif

  !--- rhokapD = dust density x dust kappa (dust extinction coefficient, per length)
  if (par%DGR > 0.0_wp) then
     call create_shared_mem(grid%rhokapD, [grid%nx,grid%ny,grid%nz])
  endif

  if (mpar%h_rank == 0) then
     if (.not. allocated(Temp))    allocate(Temp(grid%nx,grid%ny,grid%nz))
     if (.not. allocated(vtherm))  allocate(vtherm(grid%nx,grid%ny,grid%nz))
     if (.not. allocated(opacity)) allocate(opacity(grid%nx,grid%ny,grid%nz))
     Temp(:,:,:)    = 0.0_wp
     vtherm(:,:,:)  = 0.0_wp
     opacity(:,:,:) = 0.0_wp

     if (.not. allocated(xx)) allocate(xx(grid%nx))
     if (.not. allocated(yy)) allocate(yy(grid%ny))
     if (.not. allocated(zz)) allocate(zz(grid%nz))
     if (.not. allocated(rr)) allocate(rr(grid%nx,grid%ny,grid%nz))
     do i=1,grid%nx
        xx(i) = (grid%xface(i)+grid%xface(i+1))/2.0_wp
     enddo
     do j=1,grid%ny
        yy(j) = (grid%yface(j)+grid%yface(j+1))/2.0_wp
     enddo
     do k=1,grid%nz
        zz(k) = (grid%zface(k)+grid%zface(k+1))/2.0_wp
     enddo
     do k=1,grid%nz
     do j=1,grid%ny
     do i=1,grid%nx
        if (trim(par%geometry) == 'cylinder') then
           rr(i,j,k) = sqrt(xx(i)**2 + yy(j)**2)
        else
           rr(i,j,k) = sqrt(xx(i)**2 + yy(j)**2 + zz(k)**2)
        endif
     enddo
     enddo
     enddo
  endif

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  !--- reference temperature and Doppler frequencey.
  grid%Dfreq_ref      = 0.12843374_wp*sqrt(par%temperature)/sqrt(1.0_wp*par%atom_no)/(par%lambda0*um2km)
#ifdef FINE_STRUCTURE
  grid%DnuHK_ref      = DnuHK_Hz/grid%Dfreq_ref
  grid%DnuHK_ref_half = DnuHK_Hz_half/grid%Dfreq_ref
#endif

  !--- To mask locations where Lya photons will be completely destroyed.
  if (trim(par%geometry) == 'spherical_atmosphere') then
     call create_shared_mem(grid%mask, [grid%nx,grid%ny,grid%nz])
  endif
 
  !---
  !--- (1) setup temperature, Doppler-frequency, Voigt-a parameter.
  !---
  if (mpar%h_rank == 0) then
     if (len_trim(par%temp_file) > 0) then
        if (trim(par%geometry) == 'plane_atmosphere') then
           call read_plane_atmosphere_data(trim(par%temp_file),Temp,grid)
        else if (trim(par%geometry) == 'spherical_atmosphere') then
           call read_spherical_atmosphere_data(trim(par%temp_file),Temp,grid)
        else
           call read_3D(trim(par%temp_file),Temp,reduce_factor=par%reduce_factor)
        endif
     else
        Temp(:,:,:) = par%temperature
     endif
     !--- Temp should have positive values everywhere in the cells. (Do not remove this.)
     if (any(Temp <= 0.0_wp)) then
        where (Temp <= 0.0_wp)
           Temp = par%temperature
        endwhere
     endif

     do k=1,grid%nz
     do j=1,grid%ny
     do i=1,grid%nx
        !--- Dfreq should have non-zero values even in zero-density cells.
        vtherm(i,j,k)           = 0.12843374_wp*sqrt(Temp(i,j,k))/sqrt(1.0_wp*par%atom_no)
        grid%Dfreq(i,j,k)       = vtherm(i,j,k)/(par%lambda0*um2km)
        grid%Dfreq_ratio(i,j,k) = grid%Dfreq(i,j,k)/grid%Dfreq_ref
        grid%voigt_a(i,j,k)     = (par%A21/fourpi)/grid%Dfreq(i,j,k)
     enddo
     enddo
     enddo
  endif

  !---
  !--- (2) setup density
  !---
  if (mpar%h_rank == 0) then
     if (len_trim(par%dens_file) > 0) then
        if (trim(par%geometry) == 'plane_atmosphere') then
           call read_plane_atmosphere_data(trim(par%dens_file),grid%rhokap,grid)
        else if (trim(par%geometry) == 'spherical_atmosphere') then
           call read_spherical_atmosphere_data(trim(par%dens_file),grid%rhokap,grid)
           do k=1,grid%nz
           do j=1,grid%ny
           do i=1,grid%nx
              if (rr(i,j,k) <= grid%rmin) grid%mask(i,j,k) = -1
           enddo
           enddo
           enddo
        else
           call read_3D(trim(par%dens_file),grid%rhokap,reduce_factor=par%reduce_factor,centering=par%centering)
        endif
        !grid%rhokap(:,:,:) = grid%rhokap(:,:,:) * par%distance2cm
        !if (par%DGR > 0.0_wp) grid%rhokapD(:,:,:) = grid%rhokap(:,:,:)*par%cext_dust * par%DGR
        !--- in polaris machine, the above statement causes segmentaiton error. (2020-10-02).
        do k=1,grid%nz
           grid%rhokap(:,:,k) = grid%rhokap(:,:,k) * par%distance2cm
           if (par%DGR > 0.0_wp) grid%rhokapD(:,:,k) = grid%rhokap(:,:,k)*par%cext_dust * par%DGR
        enddo

        !--- Shear Effect for TIGRESS data.
        !--- The shear effect should be considered in raytrace module. (2017-07-14). bug-fixed (2017-07-15).
        if (par%Omega /= 0.0_wp) then
           if (par%distance_unit /= 'kpc') par%Omega = par%Omega * (par%distance2cm/kpc2cm)
           par%Omega = par%q * par%Omega * grid%xrange
        endif
     else
        par%distance_unit  = ''
        par%distance2cm    = 1.0_wp
        grid%rhokap(:,:,:) = 1.0_wp
        if (par%DGR > 0.0_wp) grid%rhokapD(:,:,:) = par%cext_dust * par%DGR
     endif

     if (par%rmax > 0.0_wp) then
        radial_cell_width = minval([grid%dx,grid%dy,grid%dz])
        if (par%rmin > 0.0_wp) then
           density_rmin   = par%rmin + radial_cell_width/2.1_wp
        else
           density_rmin   = par%rmin
        endif
        density_rmax      = par%rmax - radial_cell_width/2.1_wp
        do k=1,grid%nz
        do j=1,grid%ny
        do i=1,grid%nx
           !if (rr(i,j,k) < par%rmin .or. rr(i,j,k) > par%rmax) then
           !--- 2020.08.31. we are using an approximate sphere. a cell has a radial width of dz.
           if (rr(i,j,k) < density_rmin .or. rr(i,j,k) > density_rmax) then
              grid%rhokap(i,j,k) = 0.0_wp
              if (par%DGR > 0.0_wp) grid%rhokapD(i,j,k) = 0.0_wp
           endif
           if (par%density_rscale > 0.0_wp) then
              grid%rhokap(i,j,k) = grid%rhokap(i,j,k)*exp(-rr(i,j,k)/par%density_rscale)
           endif
        enddo
        enddo
        enddo
     endif

     !--- Calculate Neutral Fraction of Hydrogen in Collisional Ionization Equilibtrium (2017-07-14)
     !--- Note: Dust amount is assumed to be proportional to total (ionized + neutral) Hydrogen.
     if (par%use_cie_condition) then
        do k=1,grid%nz
        do j=1,grid%ny
        do i=1,grid%nx
           ! Collisional Ionization rate and case A Recombination rate from Draine (p. 134 and p. 139)
           T4       = Temp(i,j,k)/1.0d4
           k_ionize = 5.84862e-9_wp * sqrt(t4)  * exp(-15.78215d0/t4)
           k_recomb = 4.13e-13_wp   * (T4)**(-0.7131_wp-0.0115_wp*log(T4))
           grid%rhokap(i,j,k) = grid%rhokap(i,j,k)*(k_recomb/(k_ionize + k_recomb))
           !if (par%DGR > 0.0_wp) grid%rhokapD(i,j,k) = grid%rhokapD(i,j,k)*(k_recomb/(k_ionize + k_recomb))
        enddo
        enddo
        enddo
     endif

     !--- opacity = optical depth  per unit distance at x = 0.
     do k=1,grid%nz
     do j=1,grid%ny
     do i=1,grid%nx
        grid%rhokap(i,j,k) = grid%rhokap(i,j,k)/grid%Dfreq(i,j,k)*par%cross0
        opacity(i,j,k)     = grid%rhokap(i,j,k)*voigt(0.0_wp,grid%voigt_a(i,j,k))
     enddo
     enddo
     enddo

     if (par%rmax > 0.0_wp .and. par%rmin > 0.0_wp) then
        opac_length = par%rmax - par%rmin
     else if (par%rmax > 0.0_wp) then
        opac_length = par%rmax
     else if (grid%zmax == -grid%zmin) then
        opac_length = grid%zrange/2.0_wp
     else
        opac_length = grid%zrange
     endif

     !--- center indices.
     if (par%xyz_symmetry) then
        nxcen = 1
        nycen = 1
        nzcen = 1
     else
        nxcen = (grid%nx+1)/2
        nycen = (grid%ny+1)/2
        nzcen = (grid%nz+1)/2
     endif

     !--- We need to set tauhomoe = -999.0, taumax = -999.0 and N_HIhomo = -999.0 in define.f90 (2020.08.27)
     if (par%taumax > 0.0_wp) then
        !--- the case when taumax is given.
        !--- taumax is the radial or z-direction optical depth.
        if (par%xyz_symmetry) then
           if ((grid%nz/2)*2 == grid%nz) then
              opac_norm = par%taumax/(sum(opacity(nxcen,nycen,:))*grid%dz)
           else
              opac_norm = par%taumax/((sum(opacity(nxcen,nycen,:))-opacity(nxcen,nycen,1)/2.0_wp)*grid%dz)
           endif
        else if (grid%zmax == -grid%zmin) then
           opac_norm = (2.0_wp*par%taumax)/(sum(opacity(nxcen,nycen,:))*grid%dz)
        else
           opac_norm = par%taumax/(sum(opacity(nxcen,nycen,:))*grid%dz)
        endif
        opacity(:,:,:)     = opacity(:,:,:)     * opac_norm
        grid%rhokap(:,:,:) = grid%rhokap(:,:,:) * opac_norm
        if (par%DGR > 0.0_wp) grid%rhokapD(:,:,:) = grid%rhokapD(:,:,:) * opac_norm
     else if (par%tauhomo > 0.0_wp) then
        !--- the case when tauhomo is given.
        !--- tauhomo is the optical depth obtained for an homogenous medium that contains the same mass.
        opac_avg           = sum(opacity)/count(opacity > 0.0_wp) * opac_length
        opac_norm          = par%tauhomo/opac_avg
        opacity(:,:,:)     = opacity(:,:,:)     * opac_norm
        grid%rhokap(:,:,:) = grid%rhokap(:,:,:) * opac_norm
        if (par%DGR > 0.0_wp) grid%rhokapD(:,:,:) = grid%rhokapD(:,:,:) * opac_norm
     else if (par%N_HImax > 0.0_wp) then
        !--- the case when N(HI)_max is given.
        !--- N_HImax is the radial or z-direction column density.
        if (par%xyz_symmetry) then
           if ((grid%nz/2)*2 == grid%nz) then
              opac_norm = par%N_HImax/(sum(grid%rhokap(nxcen,nycen,:)*grid%Dfreq(nxcen,nycen,:))*grid%dz/par%cross0)
           else
              opac_norm = par%N_HImax/((sum(grid%rhokap(nxcen,nycen,:)*grid%Dfreq(nxcen,nycen,:)) &
                                              -grid%rhokap(nxcen,nycen,1)*grid%Dfreq(nxcen,nycen,1)/2.0_wp)*grid%dz/par%cross0)
           endif
        else if (grid%zmax == -grid%zmin) then
           opac_norm = (2.0_wp*par%N_HImax)/(sum(grid%rhokap(nxcen,nycen,:)*grid%Dfreq(nxcen,nycen,:))*grid%dz/par%cross0)
        else
           opac_norm = par%N_HImax/(sum(grid%rhokap(nxcen,nycen,:)*grid%Dfreq(nxcen,nycen,:))*grid%dz/par%cross0)
        endif
        opacity(:,:,:)     = opacity(:,:,:)     * opac_norm
        grid%rhokap(:,:,:) = grid%rhokap(:,:,:) * opac_norm
        if (par%DGR > 0.0_wp) grid%rhokapD(:,:,:) = grid%rhokapD(:,:,:) * opac_norm
     else if (par%N_HIhomo > 0.0_wp) then
        !--- the case when N(HI)_homo is given.
        !--- N_HIhomo is the column density obtained for an homogenous medium that contains the same mass.
        dens               = sum(grid%rhokap*grid%Dfreq)/count(grid%rhokap > 0.0_wp)/par%cross0
        opac_norm          = par%N_HIhomo/(dens * opac_length)
        opacity(:,:,:)     = opacity(:,:,:)     * opac_norm
        grid%rhokap(:,:,:) = grid%rhokap(:,:,:) * opac_norm
        if (par%DGR > 0.0_wp) grid%rhokapD(:,:,:) = grid%rhokapD(:,:,:) * opac_norm
     endif

     !--- calculate optical depth
     tauhomo = sum(opacity)/count(opacity > 0.0_wp) * opac_length
     if (par%xyz_symmetry) then
        taupole = sum(opacity(nxcen,nycen,:))*grid%dz
        if ((grid%nz/2)*2 /= grid%nz) taupole = taupole - opacity(nxcen,nycen,1)*grid%dz/2.0_wp
     else if (grid%zmax == -grid%zmin) then
        taupole = sum(opacity(nxcen,nycen,:))*grid%dz/2.0_wp
     else
        taupole = sum(opacity(nxcen,nycen,:))*grid%dz
     endif

     !--- calculate N(HI)
     N_HIhomo  = sum(grid%rhokap*grid%Dfreq)/count(grid%rhokap > 0.0_wp)/par%cross0 * opac_length
     if (par%xyz_symmetry) then
        N_HIpole = sum(grid%rhokap(nxcen,nycen,:)*grid%Dfreq(nxcen,nycen,:))*grid%dz/par%cross0
        if ((grid%nz/2)*2 /= grid%nz) N_HIpole = N_HIpole &
                                                  - grid%rhokap(nxcen,nycen,1)*grid%Dfreq(nxcen,nycen,1)*grid%dz/2.0_wp/par%cross0
     else if (grid%zmax == -grid%zmin) then
        N_HIpole = sum(grid%rhokap(nxcen,nycen,:)*grid%Dfreq(nxcen,nycen,:))*grid%dz/2.0_wp/par%cross0
     else
        N_HIpole = sum(grid%rhokap(nxcen,nycen,:)*grid%Dfreq(nxcen,nycen,:))*grid%dz/par%cross0
     endif

     if (par%DGR  > 0.0_wp) then
        tauhomo_dust = sum(grid%rhokapD)/count(grid%rhokapD > 0.0_wp) * opac_length
        if (par%xyz_symmetry) then
           taupole_dust = sum(grid%rhokapD(nxcen,nycen,:))*grid%dz
           if ((grid%nz/2)*2 /= grid%nz) taupole_dust = taupole_dust - grid%rhokapD(nxcen,nycen,1)*grid%dz/2.0_wp
        else if (grid%zmax == -grid%zmin) then
           taupole_dust = sum(grid%rhokapD(nxcen,nycen,:))*grid%dz/2.0_wp
        else
           taupole_dust = sum(grid%rhokapD(nxcen,nycen,:))*grid%dz
        endif
     else
        tauhomo_dust = 0.0_wp
        taupole_dust = 0.0_wp
     endif

     if (par%taumax   <= 0.0_wp) par%taumax   = taupole
     if (par%tauhomo  <= 0.0_wp) par%tauhomo  = tauhomo
     if (par%N_HImax  <= 0.0_wp) par%N_HImax  = N_HIpole
     if (par%N_HIhomo <= 0.0_wp) par%N_HIhomo = N_HIhomo
  endif
  call MPI_BCAST(par%N_HImax,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(par%N_HIhomo, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(par%taumax,   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(par%tauhomo,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  !---
  !--- (3) setup velocity field
  !---
  if (mpar%h_rank == 0) then
     if (len_trim(par%velo_file) > 0) then
        if (trim(par%geometry) == 'plane_atmosphere') then
           call read_plane_atmosphere_data(trim(par%velo_file),grid%vfz,grid)
           grid%vfx(:,:,:) = 0.0_wp
           grid%vfy(:,:,:) = 0.0_wp
        else if (trim(par%geometry) == 'spherical_atmosphere') then
           call read_spherical_atmosphere_velocity(trim(par%velo_file),grid%vfx,grid%vfy,grid%vfz,grid)
        else
           call read_velocity(trim(par%velo_file),grid%vfx,grid%vfy,grid%vfz,reduce_factor=par%reduce_factor)
        endif
        grid%vfx(:,:,:) = grid%vfx(:,:,:)/vtherm(:,:,:)
        grid%vfy(:,:,:) = grid%vfy(:,:,:)/vtherm(:,:,:)
        grid%vfz(:,:,:) = grid%vfz(:,:,:)/vtherm(:,:,:)
     else if (trim(par%velocity_type) == 'hubble') then
        if (par%rmax <= 0.0_wp) then
           par%rpeak = maxval([grid%xmax,grid%ymax,grid%zmax])
        else
           par%rpeak = par%rmax
        endif
        do k=1,grid%nz
        do j=1,grid%ny
        do i=1,grid%nx
           if (grid%rhokap(i,j,k) > 0.0_wp) then
              grid%vfx(i,j,k) = (par%Vexp/vtherm(i,j,k)) * xx(i) / par%rpeak
              grid%vfy(i,j,k) = (par%Vexp/vtherm(i,j,k)) * yy(j) / par%rpeak
              grid%vfz(i,j,k) = (par%Vexp/vtherm(i,j,k)) * zz(k) / par%rpeak
           endif
        enddo
        enddo
        enddo
     else if (trim(par%velocity_type) == 'ssh') then
        !--- galaxy model in Song, Seon, & Hwang (2020).
        do k=1,grid%nz
        do j=1,grid%ny
        do i=1,grid%nx
           if (grid%rhokap(i,j,k) > 0.0_wp) then
              if (rr(i,j,k) < par%rpeak) then
                 Vscale          = par%Vpeak / par%rpeak
                 grid%vfx(i,j,k) = Vscale/vtherm(i,j,k) * xx(i)
                 grid%vfy(i,j,k) = Vscale/vtherm(i,j,k) * yy(j)
                 grid%vfz(i,j,k) = Vscale/vtherm(i,j,k) * zz(k)
              else
                 Vscale          = par%Vpeak + par%DeltaV * (rr(i,j,k) - par%rpeak) / (par%rmax - par%rpeak)
                 grid%vfx(i,j,k) = Vscale/vtherm(i,j,k) * xx(i)/rr(i,j,k)
                 grid%vfy(i,j,k) = Vscale/vtherm(i,j,k) * yy(j)/rr(i,j,k)
                 grid%vfz(i,j,k) = Vscale/vtherm(i,j,k) * zz(k)/rr(i,j,k)
              endif
           endif
        enddo
        enddo
        enddo
     else if (trim(par%velocity_type) == 'constant_radial') then
        do k=1,grid%nz
        do j=1,grid%ny
        do i=1,grid%nx
           !--- 2020.08.30.
           !--- Here, it is safe to apply "rho > 0" condition for thin shell model with a constant velocity (Dijkstra & Loeb 2008)
           !--- For an odd number of cells, "rr" at (x,y,z)=(0,0,0) can be non-zero due to numerical instability.
           !--- Then, the constant-velocity model will give a constant velocity at (0,0,0) if rho > 0 condition is not applied.
           !---       the velocity direction will be randomly drawn due to numerical instability.
           !--- In this case, if comoving_source = .true. (default) is assumed, this will yield a weird result
           !---       (asymmetric peeled-off image).
           !--- To avoid this issue, we need to assign velocity only in non-zero-density cells, or set comoving_source = .false.
           !---       it would be even much safer to assume v = 0 for rr < dz.
           !if (rr(i,j,k) > 0.0_wp) then
           !if (rr(i,j,k) > 0.0_wp .and. grid%rhokap(i,j,k) > 0.0_wp) then
           !---
           if (rr(i,j,k) > grid%dz/10.0_wp .and. grid%rhokap(i,j,k) > 0.0_wp) then
              grid%vfx(i,j,k) = par%Vexp/vtherm(i,j,k) * xx(i) / rr(i,j,k)
              grid%vfy(i,j,k) = par%Vexp/vtherm(i,j,k) * yy(j) / rr(i,j,k)
              grid%vfz(i,j,k) = par%Vexp/vtherm(i,j,k) * zz(k) / rr(i,j,k)
           endif
        enddo
        enddo
        enddo
     endif
  endif

  !---
  !--- (4) Emissivity
  !---
  if (len_trim(par%emiss_file) > 0) then
     call create_shared_mem(grid%Pem, [grid%nx,grid%ny,grid%nz])
     if (mpar%h_rank == 0) then
        !--- grid%Pem : emissivity and probability distribution function for rejection method.
        if (trim(par%geometry) == 'plane_atmosphere') then
           call read_plane_atmosphere_data(trim(par%emiss_file),grid%Pem,grid)
        else if (trim(par%geometry) == 'spherical_atmosphere') then
           call read_spherical_atmosphere_data(trim(par%emiss_file),grid%Pem,grid)
        else
           call read_3D(trim(par%emiss_file),grid%Pem,reduce_factor=par%reduce_factor,centering=par%centering)
        endif
        Pem_max  = maxval(grid%Pem)
        grid%Pem = grid%Pem/Pem_max
     endif
  endif

  !--- save input grid, if you want to check the inputs were correctly given.
  if (par%save_input_grid .and. mpar%p_rank == 0) then
     !-- density * cross-section / unit length at x = 0.
     call fits_open_new(unit,trim(par%base_name)//'_opac.fits.gz',status)
     call fits_append_image(unit,opacity,status, bitpix=-32)
     call fits_close(unit,status)

     !-- density or density per unit length
     opacity(:,:,:) = grid%rhokap(:,:,:)*grid%Dfreq(:,:,:)/par%cross0 / par%distance2cm
     call fits_open_new(unit,trim(par%base_name)//'_dens.fits.gz',status)
     call fits_append_image(unit,opacity,status,bitpix=-32)
     call fits_close(unit,status)

     !-- temperature in K.
     call fits_open_new(unit,trim(par%base_name)//'_temp.fits.gz',status)
     call fits_append_image(unit,Temp,status,bitpix=-32)
     call fits_close(unit,status)

     !-- velocity
     !-- let's use "opacity" array temperarily.
     call fits_open_new(unit,trim(par%base_name)//'_vfx.fits.gz',status)
     opacity = grid%vfx*vtherm
     call fits_append_image(unit,opacity,status,bitpix=-32)
     call fits_close(unit,status)

     call fits_open_new(unit,trim(par%base_name)//'_vfy.fits.gz',status)
     opacity = grid%vfy*vtherm
     call fits_append_image(unit,opacity,status,bitpix=-32)
     call fits_close(unit,status)

     call fits_open_new(unit,trim(par%base_name)//'_vfz.fits.gz',status)
     opacity = grid%vfz*vtherm
     call fits_append_image(unit,opacity,status,bitpix=-32)
     call fits_close(unit,status)

     !-- emissivity
     if (associated(grid%Pem)) then
        !-- let's use "opacity" array temperarily.
        opacity = grid%Pem*Pem_max
        call fits_open_new(unit,trim(par%base_name)//'_emiss.fits.gz',status)
        call fits_append_image(unit,opacity,status,bitpix=-32)
        call fits_close(unit,status)
     endif
  endif

  !--- The following arrays are no longer needed. The arrays are deallocated.
  if (mpar%h_rank == 0) then
     if (allocated(Temp))    deallocate(Temp)
     if (allocated(vtherm))  deallocate(vtherm)
     if (allocated(opacity)) deallocate(opacity)
     if (allocated(xx))      deallocate(xx)
     if (allocated(yy))      deallocate(yy)
     if (allocated(zz))      deallocate(zz)
     if (allocated(rr))      deallocate(rr)
  endif

  if (par%save_all_photons) then
     !-- Note that allph arrays are all declared to be shared memory (2020-11-02).
     call create_shared_mem(allph%rp,          [int(par%nphotons)])
     call create_shared_mem(allph%xfreq1,      [int(par%nphotons)])
     call create_shared_mem(allph%xfreq2,      [int(par%nphotons)])
     call create_shared_mem(allph%nscatt_HI,   [int(par%nphotons)])
     call create_shared_mem(allph%nscatt_dust, [int(par%nphotons)])
     if (trim(par%source_geometry) /= 'point') then
        call create_shared_mem(allph%rp0, [int(par%nphotons)])
     endif
     if (par%use_stokes) then
        call create_shared_mem(allph%I,   [int(par%nphotons)])
        call create_shared_mem(allph%Q,   [int(par%nphotons)])
        call create_shared_mem(allph%U,   [int(par%nphotons)])
        call create_shared_mem(allph%V,   [int(par%nphotons)])
     endif
  endif
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  ! mean values
  ! updated 2017-06-20.
  !grid%voigt_amean = sum(grid%voigt_a, grid%voigt_a > 0.0_wp)/count(grid%voigt_a > 0.0_wp)
  !grid%Dfreq_mean  = sum(grid%Dfreq,   grid%Dfreq   > 0.0_wp)/count(grid%Dfreq   > 0.0_wp)
  grid%voigt_amean = (par%A21/fourpi)/grid%Dfreq_ref
  grid%Dfreq_mean  = grid%Dfreq_ref
  atau0            = grid%voigt_amean * par%tauhomo
  par%atau3        = (atau0)**(1.0_wp/3.0_wp)

  ! set up frequency grid.
  grid%nxfreq = par%nxfreq
  if (par%xfreq_max == 0.0_wp .and. par%xfreq_min == 0.0_wp) then
     if (par%Vexp == 0.0_wp) then
        par%xfreq_max  = floor(4.5_wp * par%atau3)+1
        par%xfreq_min  = -par%xfreq_max
     else if (par%Vexp > 0.0_wp) then    ! expanding
        par%xfreq_max  = floor(4.5_wp * par%atau3)+1
        par%xfreq_min  = -(floor(10.5_wp * par%atau3)+1)
     else if (par%Vexp < 0.0_wp) then    ! infalling
        par%xfreq_max  = floor(10.5_wp * par%atau3)+1
        par%xfreq_min  = -(floor(4.5_wp * par%atau3)+1)
     endif
     if (trim(par%spectral_type) == 'continuum') then
        if (par%Vexp > 0.0_wp) then         ! expanding
           par%xfreq_max  = floor(10.5_wp * par%atau3)+1
           par%xfreq_min  = -(floor(4.5_wp * par%atau3)+1)
        else if (par%Vexp < 0.0_wp) then    ! infalling
           par%xfreq_max  = floor(4.5_wp * par%atau3)+1
           par%xfreq_min  = -(floor(10.5_wp * par%atau3)+1)
        endif
     endif
     !if (mpar%p_rank == 0) write(*,*) 'xfreq = ', par%xfreq_min, par%xfreq_max
  endif
  if ((par%nxfreq/2)*2 == par%nxfreq .or. par%nxfreq == 1) then
     grid%dxfreq    = (par%xfreq_max - par%xfreq_min)/par%nxfreq
     grid%xfreq_max = par%xfreq_max
     grid%xfreq_min = par%xfreq_min
  else
     grid%dxfreq    = (par%xfreq_max - par%xfreq_min)/(par%nxfreq-1)
     grid%xfreq_min = par%xfreq_min - grid%dxfreq/2.0_wp
     grid%xfreq_max = par%xfreq_max + grid%dxfreq/2.0_wp
     par%xfreq_min  = grid%xfreq_min
     par%xfreq_max  = grid%xfreq_max
  endif

  !--- allocate frequency, velocity, and wavelength arrays for spectral outputs.
  call create_shared_mem(grid%xfreq,    [grid%nxfreq])
  call create_shared_mem(grid%velocity, [grid%nxfreq])
  call create_shared_mem(grid%dlambda,  [grid%nxfreq])
  if (mpar%h_rank == 0) then
     grid%xfreq(:) = [ ((i-0.5_wp)*grid%dxfreq + grid%xfreq_min, i=1, grid%nxfreq) ]
     !--- numerically, to make zero real zero.
     if ( (grid%xfreq(grid%nxfreq) == abs(grid%xfreq(1))) .and. &
          ((grid%nxfreq/2)*2 /= grid%nxfreq) ) grid%xfreq((grid%nxfreq+1)/2) = 0.0_wp
     grid%velocity = -0.12843374_wp*sqrt(par%temperature) * grid%xfreq
     grid%dlambda  =  grid%velocity * (par%lambda0 * 1e4  / speedc)
  endif

  !--- J(freq,x,y,z)  = mean intensity
  !--- Pa(x,y,z)      = number of scattering (P_alpha) per atom.
  !--- Jin(xfreq)     = input spectrum (intensity unit)
  !--- Jout(xfreq)    = escaped spectrum (intensity unit)
  call create_mem(grid%Jout, [grid%nxfreq])
  if (par%save_Jin) then
     call create_mem(grid%Jin,  [grid%nxfreq])
  endif
  !--- Jabs(xfreq)    = spectrum absorbed by dust (intensity unit)
  if (par%DGR > 0.0_wp .and. par%save_Jabs) then
     call create_mem(grid%Jabs, [grid%nxfreq])
  endif
  !--- Jabs2(xfreq)   = spectrum absorbed by molecular zone of the exoplanet atmosphere (intensity unit)
  !    2021.04.22
  if (trim(par%geometry) == 'plane_atmosphere' .or. trim(par%geometry) == 'spherical_atmosphere') then
     call create_mem(grid%Jabs2, [grid%nxfreq])
  endif

  ! continuum parameters
  if (trim(par%spectral_type) == 'continuum') then
     cpar%x0        = (grid%xfreq_min + grid%xfreq_max)/2.0_wp
     cpar%slope_max = 2.0_wp*cpar%wgt0/(grid%xfreq_max - grid%xfreq_min)
     if (abs(cpar%slope) > cpar%slope_max) then
        if (cpar%slope > 0.0_wp) cpar%slope =  cpar%slope_max
        if (cpar%slope < 0.0_wp) cpar%slope = -cpar%slope_max
        print*,'slope = ',cpar%slope
     endif
  endif

  ! parameters for core skipping.
  if (par%core_skip_global) then
     ! let's use global cell value for atau0.
     if (atau0 <= 1.0_wp) then
        grid%xcrit  = 0.0_wp
        grid%xcrit2 = 0.0_wp
     else
        if (atau0 <= 60.0_wp) then
           xi  = 0.6_wp
           chi = 1.2_wp
        else
           xi  = 1.4_wp
           chi = 0.6_wp
        endif
        grid%xcrit  = 0.02_wp * exp(xi*(log(atau0))**chi)
        grid%xcrit2 = grid%xcrit ** 2
     endif
  else
     ! let's use local cell value for atau0. (2017-06-12)
     atau0_cell = atau0/(grid%xmax/grid%dx)
     if (atau0_cell <= 1.0_wp) then
        grid%xcrit  = 0.0_wp
        grid%xcrit2 = 0.0_wp
     else
        if (atau0_cell <= 60.0_wp) then
           xi  = 0.6_wp
           chi = 1.2_wp
        else
           xi  = 1.4_wp
           chi = 0.6_wp
        endif
        grid%xcrit  = 0.02_wp * exp(xi*(log(atau0_cell))**chi)
        grid%xcrit2 = grid%xcrit ** 2
     endif
  endif

#if defined (CALCJ) || defined (CALCP) || defined (CALCPnew)
  !--- Jx and Pa arrays should be created after grid%nxfreq, grid%nx,ny,nz parameters are all defined.
  call create_JPa_mem(grid)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif

  !--- print out important parameters.
  if (mpar%p_rank == 0) then
     write(6,'(a,es16.3)') 'voigt_a          : ',grid%voigt_amean
     write(6,'(a,es16.3)') 'temperature (K)  : ',par%temperature
     write(6,'(a,es16.3)') 'N(HI)_pole       : ',N_HIpole
     write(6,'(a,es16.3)') 'N(HI)_homo       : ',N_HIhomo
     write(6,'(a,es16.3)') 'tau_pole  (HI)   : ',taupole
     write(6,'(a,es16.3)') 'tau_homo  (HI)   : ',tauhomo
     write(6,'(a,es16.3)') 'tau_pole  (dust) : ',taupole_dust
     write(6,'(a,es16.3)') 'tau_homo  (dust) : ',tauhomo_dust
  endif
  end subroutine grid_create

!++++++++++++++++++++++++++++++++++++++++++
#if defined (CALCJ) || defined (CALCP) || defined (CALCPnew)
  subroutine create_JPa_mem(grid)
  use define
  type(grid_type), intent(inout) :: grid
  ! local variables
  real(kind=wp), allocatable :: xx(:),yy(:),zz(:),rr_sph(:,:,:),rr_cyl(:,:)
  integer       :: i,j,k,nadd
  real(kind=wp) :: dr,roff

  grid%geometry_JPa = par%geometry_JPa
  if (grid%geometry_JPa == 3 .and. par%xyz_symmetry) then
     call create_shared_mem(grid%ncount3D, [grid%nx,grid%ny,grid%nz])
  else if (grid%geometry_JPa == 2) then
     grid%nr   = maxval([grid%nx,   grid%ny])
     if (.not. par%xy_symmetry) then
        if ((grid%nr/2)*2 == grid%nr) then
           grid%nr = grid%nr/2
        else
           grid%nr = (grid%nr-1)/2 + 1
        endif
     endif
     grid%rmax = minval([grid%xmax, grid%ymax])
     call create_shared_mem(grid%ind_cyl,    [grid%nx,grid%ny])
     call create_shared_mem(grid%ncount_cyl, [grid%nr,grid%nz])
  else if (grid%geometry_JPa == 1) then
     grid%nr   = maxval([grid%nx, grid%ny, grid%nz])
     if (.not. par%xyz_symmetry) then
        if ((grid%nr/2)*2 == grid%nr) then
           grid%nr = grid%nr/2
        else
           grid%nr = (grid%nr-1)/2 + 1
        endif
     endif
     grid%rmax = minval([grid%xmax, grid%ymax, grid%zmax])
     call create_shared_mem(grid%ind_sph,    [grid%nx,grid%ny,grid%nz])
     call create_shared_mem(grid%ncount_sph, [grid%nr])
  else if (grid%geometry_JPa == -1) then
     call create_shared_mem(grid%ncount_plane, [grid%nz])
  endif

  !-- minor-bug fixed (2021.05.26)
  if (grid%geometry_JPa == 1 .or. grid%geometry_JPa == 2) then
     if ((grid%nr/2)*2 == grid%nr) then
        dr   = grid%rmax/grid%nr
        roff = 0.0_wp
     else
        dr   = grid%rmax/(grid%nr-0.5_wp)
        roff = -dr/2.0_wp
     endif
  endif

  if (mpar%h_rank == 0) then
     if (grid%geometry_JPa == 1 .or. grid%geometry_JPa == 2) then
        if (.not. allocated(xx)) allocate(xx(grid%nx))
        if (.not. allocated(yy)) allocate(yy(grid%ny))
        do i=1,grid%nx
           xx(i) = (grid%xface(i)+grid%xface(i+1))/2.0_wp
        enddo
        do j=1,grid%ny
           yy(j) = (grid%yface(j)+grid%yface(j+1))/2.0_wp
        enddo
        if (grid%geometry_JPa == 1) then
           if (.not. allocated(zz)) allocate(zz(grid%nz))
           do k=1,grid%nz
              zz(k) = (grid%zface(k)+grid%zface(k+1))/2.0_wp
           enddo
           if (.not. allocated(rr_sph)) allocate(rr_sph(grid%nx,grid%ny,grid%nz))
           do k=1,grid%nz
           do j=1,grid%ny
           do i=1,grid%nx
              rr_sph(i,j,k) = sqrt(xx(i)**2 + yy(j)**2 + zz(k)**2)
           enddo
           enddo
           enddo
        else if (grid%geometry_JPa == 2) then
           if (.not. allocated(rr_cyl)) allocate(rr_cyl(grid%nx,grid%ny))
           do j=1,grid%ny
           do i=1,grid%nx
              rr_cyl(i,j) = sqrt(xx(i)**2 + yy(j)**2)
           enddo
           enddo
        endif
     endif
     if (par%geometry_JPa == 1) then
        do k=1,grid%nz
        do j=1,grid%ny
        do i=1,grid%nx
           grid%ind_sph(i,j,k) = floor((rr_sph(i,j,k)-roff)/dr) + 1
           !--- Internal radiation field and scattering rate are defined only in cells with positive density. (2021.05.10)
           !--- Counting only the cells with density > 0 gives the result that is consistent with old code. (at r = rmax cell).
           !if (grid%ind_sph(i,j,k) >= 1 .and. grid%ind_sph(i,j,k) <= grid%nr) then
           if (grid%rhokap(i,j,k) > 0.0_wp .and. grid%ind_sph(i,j,k) >= 1 .and. grid%ind_sph(i,j,k) <= grid%nr) then
              if (par%xyz_symmetry) then
                 nadd = 8
                 if (i == 1 .and. (grid%nx/2)*2 /= grid%nx) nadd = nadd/2
                 if (j == 1 .and. (grid%ny/2)*2 /= grid%ny) nadd = nadd/2
                 if (k == 1 .and. (grid%nz/2)*2 /= grid%nz) nadd = nadd/2
              else
                 nadd = 1
              endif
              grid%ncount_sph(grid%ind_sph(i,j,k)) = grid%ncount_sph(grid%ind_sph(i,j,k)) + nadd
           endif
        enddo
        enddo
        enddo
     else if (par%geometry_JPa == 2) then
        do k=1,grid%nz
        do j=1,grid%ny
        do i=1,grid%nx
           grid%ind_cyl(i,j) = floor((rr_cyl(i,j)-roff)/dr) + 1
           !--- Internal radiation field and scattering rate are defined only in cells with positive density. (2021.05.10)
           !--- Counting only the cells with density > 0 gives the result that is consistent with old code. (at r = rmax cell).
           !if (grid%ind_cyl(i,j) >= 1 .and. grid%ind_cyl(i,j) <= grid%nr) then
           if (grid%rhokap(i,j,k) > 0.0_wp .and. grid%ind_cyl(i,j) >= 1 .and. grid%ind_cyl(i,j) <= grid%nr) then
              if (par%xy_symmetry) then
                 nadd = 4
                 if (i == 1 .and. (grid%nx/2)*2 /= grid%nx) nadd = nadd/2
                 if (j == 1 .and. (grid%ny/2)*2 /= grid%ny) nadd = nadd/2
              else
                 nadd = 1
              endif
              grid%ncount_cyl(grid%ind_cyl(i,j),k) = grid%ncount_cyl(grid%ind_cyl(i,j),k) + nadd
           endif
        enddo
        enddo
        enddo
     else if (par%geometry_JPa == -1) then
        do k=1,grid%nz
        do j=1,grid%ny
        do i=1,grid%nx
           !--- Internal radiation field and scattering rate are defined only in cells with positive density. (2021.05.10)
           if (grid%rhokap(i,j,k) > 0.0_wp) then
              grid%ncount_plane(k) = grid%ncount_plane(k) + 1
           endif
        enddo
        enddo
        enddo
     else if (par%geometry_JPa == 3 .and. par%xyz_symmetry) then
        do k=1,grid%nz
        do j=1,grid%ny
        do i=1,grid%nx
           nadd = 8
           if (i == 1 .and. (grid%nx/2)*2 /= grid%nx) nadd = nadd/2
           if (j == 1 .and. (grid%ny/2)*2 /= grid%ny) nadd = nadd/2
           if (k == 1 .and. (grid%nz/2)*2 /= grid%nz) nadd = nadd/2
           grid%ncount3D(i,j,k) = grid%ncount3D(i,j,k) + nadd
        enddo
        enddo
        enddo
     endif
  endif
  if (mpar%h_rank == 0) then
     if (allocated(xx))      deallocate(xx)
     if (allocated(yy))      deallocate(yy)
     if (allocated(zz))      deallocate(zz)
     if (allocated(rr_sph))  deallocate(rr_sph)
     if (allocated(rr_cyl))  deallocate(rr_cyl)
  endif
#ifdef CALCP
  select case (grid%geometry_JPa)
     case (3)
        call create_mem(grid%Pa, [grid%nx,grid%ny,grid%nz])
     case (2)
        call create_mem(grid%P2, [grid%nr,grid%nz])
     case (-1)
        call create_mem(grid%P1, [grid%nz])
     case default
        call create_mem(grid%P1, [grid%nr])
  end select
#endif
#ifdef CALCPnew
  select case (grid%geometry_JPa)
     case (3)
        call create_mem(grid%Pa_new, [grid%nx,grid%ny,grid%nz])
     case (2)
        call create_mem(grid%P2_new, [grid%nr,grid%nz])
     case (-1)
        call create_mem(grid%P1_new, [grid%nz])
     case default
        call create_mem(grid%P1_new, [grid%nr])
  end select
#endif
#ifdef CALCJ
  select case (grid%geometry_JPa)
     case (3)
        call create_mem(grid%J,  [grid%nxfreq, grid%nx,grid%ny,grid%nz])
     case (2)
        call create_mem(grid%J2, [grid%nxfreq, grid%nr,grid%nz])
     case (-1)
        call create_mem(grid%J1, [grid%nxfreq, grid%nz])
     case default
        call create_mem(grid%J1, [grid%nxfreq, grid%nr])
  end select
#endif
  end subroutine create_JPa_mem
#endif
!++++++++++++++++++++++++++++++++++++++++++

  !-----------------
  subroutine grid_destroy(grid)
  use define
  use mpi
  implicit none
  type(grid_type), intent(inout) :: grid
  integer :: ierr

  call destroy_shared_mem_all()

  if (associated(grid%Jout))  deallocate(grid%Jout)
  if (associated(grid%Jin))   deallocate(grid%Jin)
  if (associated(grid%Jabs))  deallocate(grid%Jabs)
  if (associated(grid%Jabs2)) deallocate(grid%Jabs2)

#ifdef CALCJ
  if (associated(grid%J))    deallocate(grid%J)
  if (associated(grid%J1))   deallocate(grid%J1)
  if (associated(grid%J2))   deallocate(grid%J2)
#endif

#ifdef CALCP
  if (associated(grid%Pa))   deallocate(grid%Pa)
  if (associated(grid%P1))   deallocate(grid%P1)
  if (associated(grid%P2))   deallocate(grid%P2)
#endif

#ifdef CALCPnew
  if (associated(grid%Pa_new)) deallocate(grid%Pa_new)
  if (associated(grid%P1_new)) deallocate(grid%P1_new)
  if (associated(grid%P2_new)) deallocate(grid%P2_new)
#endif
  end subroutine grid_destroy
end module grid_mod
