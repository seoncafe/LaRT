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
  type(grid_type) :: grid
  real(kind=wp)   :: dtime

  !--- for MPI
  integer :: ierr

  !--- Initialize MPI + Openmp
  !--- call MPI_INIT_THREAD(MPI_THREAD_FUNNELED, MPI_provided,ierr)
  !--- Please use MPI_THREAD_FUNNELED if MPI_provided < MPI_THREAD_MULTIPLE (comment added on 2020-11-08).
  !--- Then, par%use_master_slave will be set to false in setup_MPI_openMP_parameters.
  call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, mpar%MPI_provided,ierr)
  call time_stamp(dtime)

  !--- Initial set up
  call read_input()
  call setup_procedure()

  !===== Start OPENMP
  !--- It is safe to do global settings outside of the parallel construct. (comment added, 2020.10.19).
  !--- Warning: Don't even try to parallelize global settings. You will be messed up with the settings.
  call grid_create(grid)
  if (par%save_peeloff) then
     call observer_create()
     if (par%save_sightline_tau) call make_sightline_tau(grid)
  endif

  !--- Run Main Calculation
  call time_stamp(dtime)
  if (mpar%p_rank == 0) write(6,'(a,f8.3,a)') '---> Now starting simulation...  @ ', dtime/60.0_wp, ' mins'
  !$OMP PARALLEL default(shared)
  call run_simulation(grid)
  !$OMP END PARALLEL

  !--- Final Output
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
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
     write(6,'(2a)')        ' >>> STOP  @ ', get_date_time()
  endif

  call grid_destroy(grid)
  if (par%save_peeloff) call observer_destroy()
  call MPI_FINALIZE(ierr)
  stop
  end program main
