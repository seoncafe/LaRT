program make_clumps
!---------------------------------------------------------------------------
! Standalone driver that builds a clump population using the same Fortran
! routines LaRT uses internally, then writes the result to a FITS file
! readable by LaRT (set par%clump_input_file = '<file>' in the LaRT input
! to consume it).
!
! Usage:
!     mpirun -np N make_clumps.x clumps_input.in
!
! The input file uses the same namelist syntax as LaRT.x.  The driver only
! reads the parameters that init_clumps / generate_clumps need, namely:
!     par%rmax                      outer sphere radius [code units]
!     par%clump_radius              base clump radius
!     par%clump_tau0 / clump_NHI / clump_nH    opacity choice
!     par%temperature / clump_temperature      gas temperature
!     par%clump_N_clumps / clump_f_cov / clump_f_vol     normalisation choice
!     par%clump_sigma_v             optional velocity sigma
!     par%clump_*_profile / *_alpha / *_r0      optional radial profiles
!     par%clump_profile_file        optional tabulated profile
!     par%iseed                     RNG seed (use a fixed value for
!                                   reproducible round-trip tests)
!     par%out_file (optional)       used only to derive the output base name;
!                                   the FITS is always <base>_clumps.fits.gz
!     par%line_id                   line atomic data (defaults to ly_alpha)
!
! Output: <base_name>_clumps.fits.gz with the same schema as the file
! produced by LaRT when par%save_clump_info = .true.  The same schema is
! consumed by read_clumps_fits() inside LaRT.
!---------------------------------------------------------------------------
  use mpi
  use define
  use setup_mod, only: read_input
  use clump_mod, only: init_clumps, write_clumps_fits
  use random,    only: init_random_seed
  use utility,   only: time_stamp, get_date_time
  implicit none

  integer       :: ierr
  real(kind=wp) :: dtime

  call MPI_INIT(ierr)
  call time_stamp(dtime)

  !--- Read the namelist (also populates mpar via the standard MPI setup
  !    block in setup.f90, and initialises line atomic data).
  call read_input()

  !--- Independent RNG init: setup_procedure() does this for LaRT, but the
  !    standalone driver does not need the rest of setup_procedure (calc_voigt
  !    / do_resonance pointers are only used by the photon loop).
  call init_random_seed(par%iseed)

  if (par%rmax <= 0.0_wp) then
     if (mpar%p_rank == 0) write(*,'(a)') &
          'ERROR: par%rmax must be > 0 for the clump generator'
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end if

  if (len_trim(par%clump_input_file) > 0 .and. mpar%p_rank == 0) then
     write(*,'(a)') &
          'WARNING: par%clump_input_file is set; make_clumps.x will REGENERATE'
     write(*,'(a)') &
          '         clumps from scratch and OVERWRITE the named file.'
  end if

  !--- Build the clump population.  init_clumps performs the same work it
  !    does inside grid_create_clump: derives N_clumps, allocates shared
  !    memory, places clumps via RSA, and builds the CSR acceleration grid.
  call init_clumps(par%rmax)

  !--- Save to FITS exactly as LaRT does when par%save_clump_info = .true.
  if (mpar%p_rank == 0) then
     call write_clumps_fits(trim(par%base_name)//'_clumps.fits.gz')
     call time_stamp(dtime)
     write(*,'(a,f8.3,a)') 'Total Execution Time         : ', dtime/60.0_wp, ' mins'
     write(*,'(2a)')       ' >>> STOP  @ ', get_date_time()
  end if

  call MPI_FINALIZE(ierr)
end program make_clumps
