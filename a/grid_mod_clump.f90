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
!   grid%voigt_a     = cl_voigt_a everywhere (for do_resonance compatibility)
!   grid%Dfreq       = cl_Dfreq everywhere (for frequency calculations)
!   grid%vfx/vfy/vfz = 0 everywhere
!   grid%Dfreq_ref   = cl_Dfreq
!   par%rmax         = outer sphere radius (from input)
!   par%xmax = par%ymax = par%zmax = par%rmax (box enclosing sphere)
!---------------------------------------------------------------------------
  use grid_mod
  use clump_mod
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

  integer :: ierr

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

  !--- Create the base Cartesian grid (allocates shared-memory arrays,
  !    sets up Jout/Jin/xfreq axes, face arrays, xcrit, etc.)
  !    rhokap will be overridden to 0 below.
  call grid_create(grid)

  !--- Initialize clumps: computes cl_Dfreq, cl_voigt_a, cl_rhokap,
  !    places N clumps via RSA, builds CSR acceleration grid.
  call init_clumps(par%rmax)

  !--- Override grid arrays with uniform clump values (h_rank=0 only).
  if (mpar%h_rank == 0) then
     grid%rhokap(:,:,:)   = 0.0_wp    ! opacity is only inside clumps
     grid%Dfreq(:,:,:)    = cl_Dfreq  ! uniform thermal Doppler frequency
     grid%voigt_a(:,:,:)  = cl_voigt_a  ! uniform Voigt damping parameter
     grid%vfx(:,:,:)      = 0.0_wp    ! no background bulk velocity
     grid%vfy(:,:,:)      = 0.0_wp
     grid%vfz(:,:,:)      = 0.0_wp
  end if
  call MPI_BARRIER(mpar%hostcomm, ierr)

  !--- Update scalar reference values to clump Doppler frequency
  grid%Dfreq_ref   = cl_Dfreq
  grid%voigt_amean = cl_voigt_a
  grid%Dfreq_mean  = cl_Dfreq

  if (mpar%p_rank == 0) then
     write(*,'(a,f12.5)')  ' Clump grid: rmax      = ', par%rmax
     write(*,'(a,3i5)')    ' Clump grid: nx,ny,nz  = ', grid%nx, grid%ny, grid%nz
     write(*,'(a,es12.4)') ' Clump grid: cl_Dfreq  = ', cl_Dfreq
     write(*,'(a,f12.5,a)')' Clump grid: cl_vtherm = ', cl_vtherm, ' km/s'
     write(*,'(a,f12.5)')  ' Clump grid: voigt_a   = ', cl_voigt_a
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
