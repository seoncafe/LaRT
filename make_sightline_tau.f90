program make_sightline_tau_main
!---------------------------------------------------------------------------
! Standalone driver that produces, from the peel-off observer vantage,
! the radial sight-line maps:
!     tau_gas[nxfreq, nxim, nyim]   gas optical depth as f(velocity)
!     N_gas  [nxim,   nyim]         gas (HI) column density
!     tau_dust[nxim,  nyim]         dust optical depth (only when DGR > 0)
!
! Same algorithm and FITS layout as a normal LaRT run with
! par%save_sightline_tau = .true., but the photon-transport loop is skipped.
! Supports all three grid modes (Cartesian, AMR, clump) via the same
! procedure-pointer dispatch used inside LaRT.
!
! Usage:
!     mpirun -np N make_sightline_tau.x sightline.in
!
! `sightline.in` uses the same namelist syntax as LaRT.x.  The driver
! forces par%save_sightline_tau = .true. and par%save_peeloff = .true.
! (so that the observer object is created), and otherwise relies on the
! same parameter handling as LaRT.x.
!---------------------------------------------------------------------------
  use mpi
  use define
  use setup_mod
  use grid_mod
  use grid_mod_amr
  use grid_mod_clump
  use utility, only: time_stamp, get_date_time
  implicit none

  type(grid_type) :: grid
  real(kind=wp)   :: dtime
  integer         :: ierr

  call MPI_INIT(ierr)
  call time_stamp(dtime)

  call read_input()
  call setup_procedure()

  !--- Force the sight-line + observer paths on; the rest of LaRT's photon
  !    transport flags stay at the user's values (and are simply not used).
  par%save_peeloff       = .true.
  par%save_sightline_tau = .true.

  !--- Build the grid for the user-selected mode.  The procedure pointers
  !    set in setup_procedure handle the per-mode raytrace dispatch later.
  if (par%use_clump_medium) then
     call grid_create_clump(grid)
  else if (par%use_amr_grid) then
     call grid_create_amr(grid)
     call amr_sync_to_grid(grid)
  else
     call grid_create(grid)
  end if

  call observer_create()
  call make_sightline_tau(grid)
  call observer_destroy()

  if (par%use_clump_medium) then
     call grid_destroy_clump(grid)
  else if (par%use_amr_grid) then
     call grid_destroy_amr()
  else
     call grid_destroy(grid)
  end if

  if (mpar%p_rank == 0) then
     call time_stamp(dtime)
     write(*,'(a,f8.3,a)') 'Total Execution Time         : ', dtime/60.0_wp, ' mins'
     write(*,'(2a)')       ' >>> STOP  @ ', get_date_time()
  end if

  call MPI_FINALIZE(ierr)
end program make_sightline_tau_main
