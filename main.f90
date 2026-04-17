  program main

  use mpi
  use define
  use setup_mod
  use grid_mod
  use grid_mod_amr
  use utility

  !--- Parameter
  implicit none
  type(grid_type) :: grid
  real(kind=wp)   :: dtime

  !--- for MPI
  integer :: ierr

  call MPI_INIT(ierr)
  call time_stamp(dtime)

  !--- Initial set up
  call read_input()
  call setup_procedure()
  if (par%use_amr_grid) then
     call grid_create_amr(grid)
     call amr_sync_to_grid(grid)
  else
     call grid_create(grid)
  end if
  if (par%save_peeloff) then
     call observer_create()
     if (par%save_sightline_tau) call make_sightline_tau(grid)
  endif

  !--- Run Main Calculation
  call time_stamp(dtime)
  if (mpar%p_rank == 0) write(6,'(a,f8.3,a)') '---> Now starting simulation...  @ ', dtime/60.0_wp, ' mins'
  call run_simulation(grid)

  !--- Final Output
  if (mpar%p_rank == 0) write(6,'(a)') '---> Now Gathering Results...'
  call output_reduce(grid)

  call time_stamp(dtime)
  if (mpar%p_rank == 0) then
     call output_normalize(grid)
     par%exetime = dtime/60.0_wp
     call write_output(trim(par%out_file),grid)

     write(6,'(a,es12.4)')  'Average Number of scattering : ', par%nscatt_tot
     write(6,'(a,2es12.4)') 'Nscatt_dust, Nscatt_gas      : ', par%nscatt_dust, par%nscatt_gas
     write(6,'(a,f8.3,a)')  'Total Excution Time          : ', par%exetime,' mins'
     write(6,'(2a)')        ' >>> STOP  @ ', get_date_time()
  endif

  if (par%use_amr_grid) then
     call grid_destroy_amr()
  else
     call grid_destroy(grid)
  end if
  if (par%save_peeloff) call observer_destroy()
  call MPI_FINALIZE(ierr)
  stop
  end program main
