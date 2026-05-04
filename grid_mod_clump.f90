module grid_mod_clump
!---------------------------------------------------------------------------
! Grid setup for the clump medium mode (par%use_clump_medium = .true.).
!
! The medium consists of N spherical clumps of radius par%clump_radius
! placed uniformly at random (non-overlapping) inside a sphere of radius
! par%rmax. Outside the clumps the medium is vacuum.
!
! Approach: call grid_create (Cartesian) with the sphere box geometry, then
! override Dfreq/voigt_a/rhokap/vf arrays with uniform clump values, and
! call init_clumps to build the clump population.
!
! After grid_create_clump returns:
!   grid%rhokap      = 0 everywhere (opacity only inside clumps via clump_mod)
!   grid%voigt_a     = cl_voigt_a_ref everywhere (for do_resonance compatibility)
!   grid%Dfreq       = cl_Dfreq_ref everywhere (for frequency calculations)
!   grid%vfx/vfy/vfz = 0 everywhere
!   grid%Dfreq_ref   = cl_Dfreq_ref
!   par%rmax         = outer sphere radius (from input)
!   par%xmax = par%ymax = par%zmax = par%rmax (box enclosing sphere)
!---------------------------------------------------------------------------
  use grid_mod
  use clump_mod
  use voigt_mod, only: voigt
  implicit none

  public :: grid_create_clump, grid_destroy_clump

contains

  !===========================================================================
  subroutine grid_create_clump(grid)
  !---------------------------------------------------------------------------
  ! Set up the Cartesian shell grid and initialize the clump population.
  !---------------------------------------------------------------------------
  use define
  use line_mod
  use mpi
  implicit none
  type(grid_type), intent(inout) :: grid

  integer       :: ierr
  real(kind=wp) :: f_vol_realized, f_cov_realized, voigt0, tau_per_clump_lc

  if (par%rmax <= 0.0_wp) then
     if (mpar%p_rank == 0) write(*,*) 'ERROR: par%rmax must be > 0 for clump medium'
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end if

  !--- Force the Cartesian box to enclose the sphere [-rmax, +rmax]^3
  par%xmax = par%rmax
  par%ymax = par%rmax
  par%zmax = par%rmax

  !--- Ensure reasonable grid resolution (default 11 if not set by user)
  if (par%nx < 1) par%nx = 11
  if (par%ny < 1) par%ny = 11
  if (par%nz < 1) par%nz = 11

  !--- Force standard (no symmetry) geometry so grid%xmin = -par%rmax
  par%xyz_symmetry = .false.
  par%xy_symmetry  = .false.
  par%xy_periodic  = .false.
  par%z_symmetry   = .false.

  !--- Initialize clumps FIRST: computes cl_Dfreq_ref, cl_voigt_a_ref, cl_rhokap_ref,
  !    determines N_clumps, places clumps via RSA, builds CSR acceleration grid.
  !    Done before grid_create so that par%tauhomo can be set correctly for
  !    the auto-detection of the frequency range inside car_setup_freq_grid.
  call init_clumps(par%rmax)

  !--- Set par%tauhomo / par%taumax / par%N_gashomo / par%N_gasmax for the
  !    clumpy medium, only if the user did not specify them in the input.
  !
  !    tauhomo (homogeneous-equivalent line-center optical depth): the
  !    optical depth of a uniform medium of the same total HI mass spread
  !    over the sphere volume:
  !
  !        tauhomo = f_vol * cl_rhokap_ref * voigt(0, cl_voigt_a_ref) * R_max
  !                = (4/3) * f_cov * tau_per_clump
  !
  !    where f_vol = N_clumps * (R_cl/R_max)^3 (volume filling fraction),
  !          f_cov = (3/4) * N_clumps * (R_cl/R_max)^2 (covering factor;
  !                  expected number of clump crossings per radial path),
  !          tau_per_clump = cl_rhokap_ref * voigt(0, cl_voigt_a_ref) * R_cl
  !                        (line-center tau from clump center to surface).
  !
  !    For a clumpy medium the EXPECTED radial taumax along a random LOS
  !    equals tauhomo, although a specific realization can deviate widely
  !    (no clumps -> 0; many clumps stacked -> >> tauhomo).
  !    N_gashomo and N_gasmax are set similarly via N_HI_per_clump =
  !        cl_rhokap_ref / line%cross0 * cl_Dfreq_ref * R_cl  (column from clump
  !        center to surface; reduces to par%clump_NHI when that input is used).
  voigt0           = voigt(0.0_wp, cl_voigt_a_ref)
  f_vol_realized   = real(N_clumps, wp) * (par%clump_radius / par%rmax)**3
  f_cov_realized   = 0.75_wp * real(N_clumps, wp) * (par%clump_radius / par%rmax)**2
  tau_per_clump_lc = cl_rhokap_ref * voigt0 * par%clump_radius

  if (par%tauhomo  <= 0.0_wp) par%tauhomo  = (4.0_wp/3.0_wp) * f_cov_realized * tau_per_clump_lc
  if (par%taumax   <= 0.0_wp) par%taumax   = (4.0_wp/3.0_wp) * f_cov_realized * tau_per_clump_lc
  if (par%N_gashomo <= 0.0_wp) &
       par%N_gashomo = (4.0_wp/3.0_wp) * f_cov_realized * (cl_rhokap_ref / line%cross0) * cl_Dfreq_ref * par%clump_radius
  if (par%N_gasmax  <= 0.0_wp) &
       par%N_gasmax  = (4.0_wp/3.0_wp) * f_cov_realized * (cl_rhokap_ref / line%cross0) * cl_Dfreq_ref * par%clump_radius

  if (mpar%p_rank == 0) then
     write(*,'(a,es14.5)') ' Clump derived: tauhomo  = ', par%tauhomo
     write(*,'(a,es14.5)') ' Clump derived: taumax   = ', par%taumax
     write(*,'(a,es14.5)') ' Clump derived: N_gashomo= ', par%N_gashomo
  end if

  !--- Create the base Cartesian grid (allocates shared-memory arrays,
  !    sets up Jout/Jin/xfreq axes, face arrays, xcrit, etc.)
  !    grid%rhokap will be normalised by par%tauhomo here (waste but harmless),
  !    then overridden to 0 below.
  call grid_create(grid)

  !--- Override grid arrays with uniform clump values (h_rank=0 only).
  if (mpar%h_rank == 0) then
     grid%rhokap(:,:,:)   = 0.0_wp    ! opacity is only inside clumps
     grid%Dfreq(:,:,:)    = cl_Dfreq_ref  ! uniform thermal Doppler frequency
     grid%voigt_a(:,:,:)  = cl_voigt_a_ref  ! uniform Voigt damping parameter
     grid%vfx(:,:,:)      = 0.0_wp    ! no background bulk velocity
     grid%vfy(:,:,:)      = 0.0_wp
     grid%vfz(:,:,:)      = 0.0_wp
  end if
  call MPI_BARRIER(mpar%hostcomm, ierr)

  !--- Update scalar reference values to clump Doppler frequency
  grid%Dfreq_ref   = cl_Dfreq_ref
  grid%voigt_amean = cl_voigt_a_ref
  grid%Dfreq_mean  = cl_Dfreq_ref

  if (mpar%p_rank == 0) then
     write(*,'(a,f12.5)')  ' Clump grid: rmax      = ', par%rmax
     write(*,'(a,3i5)')    ' Clump grid: nx,ny,nz  = ', grid%nx, grid%ny, grid%nz
     write(*,'(a,es12.4)') ' Clump grid: cl_Dfreq_ref  = ', cl_Dfreq_ref
     write(*,'(a,f12.5,a)')' Clump grid: cl_vtherm_ref = ', cl_vtherm_ref, ' km/s'
     write(*,'(a,f12.5)')  ' Clump grid: voigt_a   = ', cl_voigt_a_ref
  end if

  !--- Optionally save clump positions/velocities to FITS (p_rank=0 only)
  if (par%save_clump_info .and. mpar%p_rank == 0) then
     call write_clumps_fits(trim(par%base_name)//'_clumps.fits.gz')
  end if

  end subroutine grid_create_clump
  !===========================================================================

  !===========================================================================
  subroutine grid_destroy_clump(grid)
  use define
  implicit none
  type(grid_type), intent(inout) :: grid
  call grid_destroy(grid)
  call destroy_clumps()
  end subroutine grid_destroy_clump
  !===========================================================================

end module grid_mod_clump
