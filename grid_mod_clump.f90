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
  use iofile_mod, only: io_file_extension
  implicit none

  public :: grid_create_clump, grid_destroy_clump

contains

  !===========================================================================
  subroutine grid_create_clump(grid)
  !---------------------------------------------------------------------------
  ! Set up the Cartesian shell grid and initialize the clump population.
  !---------------------------------------------------------------------------
  use define
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

  !--- Initialize clumps FIRST: computes cl_Dfreq_ref, cl_voigt_a_ref, cl_rhokap_ref,
  !    determines N_clumps, places clumps via RSA, builds CSR acceleration grid.
  !    Done before grid_create so that par%tauhomo can be set correctly for
  !    the auto-detection of the frequency range inside car_setup_freq_grid.
  call init_clumps(par%rmax)

  !--- Compute the four system-level scalars (tauhomo, taumax, N_gashomo,
  !    N_gasmax) from the per-clump arrays populated by init_clumps.
  !
  !    For loaded clumps (par%clump_input_file given): rescale_loaded_clumps_to_target
  !    inside read_clumps_info has already multiplied cl_rhokap(:) by alpha
  !    so the realized scalar matches whichever target (taumax / tauhomo /
  !    N_gasmax / N_gashomo) the user supplied; the other three scalars are
  !    derived consistently from the same distribution.
  !
  !    For generated uniform clumps (cl_radius(:) = const, cl_rhokap(:) = const),
  !    the per-clump sums collapse to the original closed-form expressions
  !    (4/3) * f_cov_shell * tau_per_clump_lc and the equivalent
  !    column-density form, so backward compatibility is preserved.  See
  !    compute_clump_scalars in clump_mod.f90 for the exact formulas.
  call compute_clump_scalars(par%tauhomo, par%taumax, par%N_gashomo, par%N_gasmax)

  if (mpar%p_rank == 0) then
     write(*,'(a,es14.5)') ' Clump derived: tauhomo   = ', par%tauhomo
     write(*,'(a,es14.5)') ' Clump derived: taumax    = ', par%taumax
     write(*,'(a,es14.5)') ' Clump derived: N_gashomo = ', par%N_gashomo
     write(*,'(a,es14.5)') ' Clump derived: N_gasmax  = ', par%N_gasmax
  end if

  !--- Create the base Cartesian grid (allocates shared-memory arrays,
  !    sets up Jout/Jin/xfreq axes, face arrays, xcrit, etc.)
  !    grid%rhokap will be normalized by par%tauhomo here (waste but harmless),
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

  !--- Optionally save clump positions/velocities (p_rank=0 only).  The
  !    extension follows par%file_format (HDF5 by default since 2026-05),
  !    matching the make_clumps.x driver.
  if (par%save_clump_info .and. mpar%p_rank == 0) then
     call write_clumps_info(trim(par%base_name)//'_clumps'// &
                            trim(io_file_extension(par%file_format)))
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
