  program main

  use mpi
  use define
  use setup_mod
  use observer_mod
  use grid_mod
  use output_sum
  use utility
  use peelingoff_mod
  use sightline_tau_mod

  !--- Parameter
  implicit none
  type(grid_type)   :: grid
  real(kind=wp)     :: dtime

  !--- for MPI
  integer :: ierr

  call MPI_INIT(ierr)
  call time_stamp(dtime)

  !--- Initial set up
  call read_input()
  call setup_procedure()
  call grid_create(grid)
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
     write(6,'(a,2es12.4)') 'Nscatt_dust, Nscatt_HI       : ', par%nscatt_dust, par%nscatt_HI
     write(6,'(a,f8.3,a)')  'Total Excution Time          : ', par%exetime,' mins'
     write(6,'(2a)')        ' >>> STOP                    : ', get_date_time()
  endif

  call grid_destroy(grid)
  if (par%save_peeloff) call observer_destroy()
  call MPI_FINALIZE(ierr)
  stop
  end program main
