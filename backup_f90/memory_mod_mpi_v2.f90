module memory_mod
  use, intrinsic :: iso_c_binding,   only: C_PTR, C_F_POINTER, C_LOC, C_ASSOCIATED
  use, intrinsic :: iso_fortran_env, only: character_storage_size, real32, real64, real128, int8, int16, int32, int64
  use mpi
  !-- 2021-06-29. Now, the routines work for both "real32" and "real64" types.
  !--             The routines for "int8", "int16", "int32", and "inte64" type array are added.
  !-- 2021-05-03. bug-fixed
  !-- 2020-11-07. Now, "window" can be obtained using "get_window." for shared memory.
  !--             individual shared memory can be destroy using "destroy_mem".
  !--             all shared memory can be destroyed using "destroy_shared_mem_all".
  !-- Rewrote the code (2020-11-02). Now, "window" is internally stored.
  !
  !-- stroage_size returns the storage size in bits and is a standard in Fortran 2008.
  !-- character_storage_size is defined in module iso_fortran_env. (character_storage_size = 8 in most cases, but not always.)
  !-- sizeof returns the size in bytes, but is not a standard. (2017-06-10)
  implicit none
  private
  integer(kind=MPI_ADDRESS_KIND), parameter :: real32_size = int(storage_size(0.0_real32)/character_storage_size, MPI_ADDRESS_KIND)
  integer(kind=MPI_ADDRESS_KIND), parameter :: real64_size = int(storage_size(0.0_real64)/character_storage_size, MPI_ADDRESS_KIND)
  integer(kind=MPI_ADDRESS_KIND), parameter :: int8_size   = int(storage_size(0_int8 )/character_storage_size, MPI_ADDRESS_KIND)
  integer(kind=MPI_ADDRESS_KIND), parameter :: int16_size  = int(storage_size(0_int16)/character_storage_size, MPI_ADDRESS_KIND)
  integer(kind=MPI_ADDRESS_KIND), parameter :: int32_size  = int(storage_size(0_int32)/character_storage_size, MPI_ADDRESS_KIND)
  integer(kind=MPI_ADDRESS_KIND), parameter :: int64_size  = int(storage_size(0_int64)/character_storage_size, MPI_ADDRESS_KIND)
  integer, parameter :: MAX_NUM_WINDOWS = 1000
  integer :: windows(MAX_NUM_WINDOWS) = MPI_WIN_NULL
  integer :: num_windows = 0
  logical :: mpi_parameters_defined = .false.
  integer :: nprocs, p_rank, h_rank, hostcomm, same_hrank_comm, num_hosts

  public create_shared_mem, create_mem, reduce_mem, &
         destroy_shared_mem_window, destroy_shared_mem_all, get_window, destroy_mem

  interface create_shared_mem
     module procedure create_shared_mem_1D_real64, create_shared_mem_2D_real64, create_shared_mem_3D_real64, &
                      create_shared_mem_4D_real64, create_shared_mem_5D_real64, &
                      create_shared_mem_1D_real32, create_shared_mem_2D_real32, create_shared_mem_3D_real32, &
                      create_shared_mem_4D_real32, create_shared_mem_5D_real32, &
                      create_shared_mem_1D_int64, create_shared_mem_2D_int64, create_shared_mem_3D_int64, &
                      create_shared_mem_1D_int32, create_shared_mem_2D_int32, create_shared_mem_3D_int32, &
                      create_shared_mem_1D_int16, create_shared_mem_2D_int16, create_shared_mem_3D_int16, &
                      create_shared_mem_1D_int8,  create_shared_mem_2D_int8,  create_shared_mem_3D_int8
  end interface create_shared_mem
  interface create_mem
     module procedure create_mem_1D_real64, create_mem_2D_real64, create_mem_3D_real64, create_mem_4D_real64, &
                      create_mem_1D_real32, create_mem_2D_real32, create_mem_3D_real32, create_mem_4D_real32, &
                      create_mem_5D_real64, create_mem_5D_real32, &
                      create_mem_1D_int64, create_mem_2D_int64, create_mem_3D_int64, &
                      create_mem_1D_int32, create_mem_2D_int32, create_mem_3D_int32, &
                      create_mem_1D_int16, create_mem_2D_int16, create_mem_3D_int16, &
                      create_mem_1D_int8,  create_mem_2D_int8,  create_mem_3D_int8
  end interface create_mem
  interface reduce_mem
     module procedure reduce_mem_1D_real64, reduce_mem_2D_real64, reduce_mem_3D_real64, reduce_mem_4D_real64, &
                      reduce_mem_1D_real32, reduce_mem_2D_real32, reduce_mem_3D_real32, reduce_mem_4D_real32, &
                      reduce_mem_5D_real64, reduce_mem_5D_real32, &
                      reduce_mem_1D_int64, reduce_mem_2D_int64, reduce_mem_3D_int64, &
                      reduce_mem_1D_int32, reduce_mem_2D_int32, reduce_mem_3D_int32, &
                      reduce_mem_1D_int16, reduce_mem_2D_int16, reduce_mem_3D_int16, &
                      reduce_mem_1D_int8,  reduce_mem_2D_int8,  reduce_mem_3D_int8
  end interface reduce_mem
  interface get_window
     module procedure get_window_1D_real64, get_window_2D_real64, get_window_3D_real64, get_window_4D_real64, &
                      get_window_1D_real32, get_window_2D_real32, get_window_3D_real32, get_window_4D_real32, &
                      get_window_5D_real64, get_window_5D_real32, &
                      get_window_1D_int64, get_window_2D_int64, get_window_3D_int64, &
                      get_window_1D_int32, get_window_2D_int32, get_window_3D_int32, &
                      get_window_1D_int16, get_window_2D_int16, get_window_3D_int16, &
                      get_window_1D_int8,  get_window_2D_int8,  get_window_3D_int8
  end interface get_window
  interface destroy_mem
     module procedure destroy_mem_1D_real64, destroy_mem_2D_real64, destroy_mem_3D_real64, destroy_mem_4D_real64, &
                      destroy_mem_1D_real32, destroy_mem_2D_real32, destroy_mem_3D_real32, destroy_mem_4D_real32, &
                      destroy_mem_5D_real64, destroy_mem_5D_real32, &
                      destroy_mem_1D_int64, destroy_mem_2D_int64, destroy_mem_3D_int64, &
                      destroy_mem_1D_int32, destroy_mem_2D_int32, destroy_mem_3D_int32, &
                      destroy_mem_1D_int16, destroy_mem_2D_int16, destroy_mem_3D_int16, &
                      destroy_mem_1D_int8,  destroy_mem_2D_int8,  destroy_mem_3D_int8
  end interface destroy_mem

contains
  !---------------------------
  subroutine init_mpi_parameters()
  integer :: ierr
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, p_rank, ierr)
  call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, hostcomm, ierr)
  call MPI_COMM_RANK(hostcomm, h_rank, ierr)
  call MPI_COMM_SPLIT(MPI_COMM_WORLD, h_rank, p_rank, same_hrank_comm, ierr)
  call MPI_COMM_SIZE(same_hrank_comm, num_hosts, ierr)
  mpi_parameters_defined = .true.
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  end subroutine init_mpi_parameters
  !---------------------------
  subroutine create_shared_mem_1D_real64(array,arrshape)
  implicit none
  real(kind=real64), pointer, intent(inout) :: array(:)
  integer,                    intent(in)    :: arrshape(1)
  integer :: ierr
 
  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit
  
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(arrshape(1), MPI_ADDRESS_KIND)*real64_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif
 
  disp_unit   = real64_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0.0_real64
  endif
  call MPI_BARRIER(hostcomm, ierr)
  end subroutine create_shared_mem_1D_real64
  !---------------------------
  subroutine create_shared_mem_2D_real64(array,arrshape)
  implicit none
  real(kind=real64), pointer, intent(inout) :: array(:,:)
  integer,                    intent(in)    :: arrshape(2)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(int(arrshape(1),int64)*arrshape(2), MPI_ADDRESS_KIND)*real64_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = real64_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0.0_real64
  endif
  call MPI_BARRIER(hostcomm, ierr)
  end subroutine create_shared_mem_2D_real64
  !---------------------------
  subroutine create_shared_mem_3D_real64(array,arrshape)
  implicit none
  real(kind=real64), pointer, intent(inout) :: array(:,:,:)
  integer,                    intent(in)    :: arrshape(3)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(int(arrshape(1),int64)*arrshape(2)*arrshape(3), MPI_ADDRESS_KIND)*real64_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = real64_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0.0_real64
  endif
  call MPI_BARRIER(hostcomm, ierr)
  end subroutine create_shared_mem_3D_real64
  !---------------------------
  subroutine create_shared_mem_4D_real64(array,arrshape)
  implicit none
  real(kind=real64), pointer, intent(inout) :: array(:,:,:,:)
  integer,                    intent(in)    :: arrshape(4)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(int(arrshape(1),int64)*arrshape(2)*arrshape(3)*arrshape(4), MPI_ADDRESS_KIND)*real64_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = real64_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0.0_real64
  endif
  call MPI_BARRIER(hostcomm, ierr)
  end subroutine create_shared_mem_4D_real64
  !---------------------------
  subroutine create_shared_mem_5D_real64(array,arrshape)
  implicit none
  real(kind=real64), pointer, intent(inout) :: array(:,:,:,:,:)
  integer,                    intent(in)    :: arrshape(5)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(int(arrshape(1),int64)*arrshape(2)*arrshape(3)*arrshape(4)*arrshape(5), MPI_ADDRESS_KIND)*real64_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = real64_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0.0_real64
  endif
  call MPI_BARRIER(hostcomm, ierr)
  end subroutine create_shared_mem_5D_real64
  !---------------------------
  subroutine create_mem_1D_real64(array,arrshape,shared_memory)
  implicit none
  real(kind=real64), pointer, intent(inout) :: array(:)
  integer,                    intent(in)    :: arrshape(1)
  logical, optional,          intent(in)    :: shared_memory
  logical :: use_shared_memory
 
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_1D_real64(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1)))
     array = 0.0_real64
  endif
  end subroutine create_mem_1D_real64
  !---------------------------
  subroutine create_mem_2D_real64(array,arrshape,shared_memory)
  implicit none
  real(kind=real64), pointer, intent(inout) :: array(:,:)
  integer,                    intent(in)    :: arrshape(2)
  logical, optional,          intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_2D_real64(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2)))
     array = 0.0_real64
  endif
  end subroutine create_mem_2D_real64
  !---------------------------
  subroutine create_mem_3D_real64(array,arrshape,shared_memory)
  implicit none
  real(kind=real64), pointer, intent(inout) :: array(:,:,:)
  integer,                    intent(in)    :: arrshape(3)
  logical, optional,          intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_3D_real64(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2),arrshape(3)))
     array = 0.0_real64
  endif
  end subroutine create_mem_3D_real64
  !---------------------------
  subroutine create_mem_4D_real64(array,arrshape,shared_memory)
  implicit none
  real(kind=real64), pointer, intent(inout) :: array(:,:,:,:)
  integer,                    intent(in)    :: arrshape(4)
  logical, optional,          intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_4D_real64(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2),arrshape(3),arrshape(4)))
     array = 0.0_real64
  endif
  end subroutine create_mem_4D_real64
  !---------------------------
  subroutine create_mem_5D_real64(array,arrshape,shared_memory)
  implicit none
  real(kind=real64), pointer, intent(inout) :: array(:,:,:,:,:)
  integer,                    intent(in)    :: arrshape(5)
  logical, optional,          intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_5D_real64(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2),arrshape(3),arrshape(4),arrshape(5)))
     array = 0.0_real64
  endif
  end subroutine create_mem_5D_real64
  !---------------------------
  subroutine reduce_mem_1D_real64(array, shared_memory)
  implicit none
  real(kind=real64), pointer,  intent(inout) :: array(:)
  logical,           optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, nsize
  logical :: do_reduce, use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory

  do_reduce = .false.
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:),nsize,MPI_REAL8,MPI_SUM,0,communicator,ierr)
     else
        call MPI_REDUCE(array(:),           0,nsize,MPI_REAL8,MPI_SUM,0,communicator,ierr)
     endif
  endif
  end subroutine reduce_mem_1D_real64
  !---------------------------
  subroutine reduce_mem_2D_real64(array, shared_memory)
  implicit none
  real(kind=real64), pointer,  intent(inout) :: array(:,:)
  logical,           optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, nsize
  logical :: do_reduce, use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory

  do_reduce = .false.
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:,:),nsize,MPI_REAL8,MPI_SUM,0,communicator,ierr)
     else
        call MPI_REDUCE(array(:,:),           0,nsize,MPI_REAL8,MPI_SUM,0,communicator,ierr)
     endif
  endif
  end subroutine reduce_mem_2D_real64
  !---------------------------
  subroutine reduce_mem_3D_real64(array, shared_memory)
  implicit none
  real(kind=real64), pointer,  intent(inout) :: array(:,:,:)
  logical,           optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, k, nsize
  logical :: do_reduce, use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory

  do_reduce = .false.
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array(:,:,1))
     do k=1, size(array(1,1,:))
        if (p_rank == 0) then
           call MPI_REDUCE(MPI_IN_PLACE,array(:,:,k),nsize,MPI_REAL8,MPI_SUM,0,communicator,ierr)
        else
           call MPI_REDUCE(array(:,:,k),           0,nsize,MPI_REAL8,MPI_SUM,0,communicator,ierr)
        endif
     enddo
  endif
  end subroutine reduce_mem_3D_real64
  !---------------------------
  subroutine reduce_mem_4D_real64(array, shared_memory)
  implicit none
  real(kind=real64), pointer,  intent(inout) :: array(:,:,:,:)
  logical,           optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, j, k, nsize
  logical :: do_reduce, use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory

  do_reduce = .false.
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array(:,:,1,1))
     do k=1, size(array(1,1,1,:))
     do j=1, size(array(1,1,:,1))
        if (p_rank == 0) then
           call MPI_REDUCE(MPI_IN_PLACE,array(:,:,j,k),nsize,MPI_REAL8,MPI_SUM,0,communicator,ierr)
        else
           call MPI_REDUCE(array(:,:,j,k),           0,nsize,MPI_REAL8,MPI_SUM,0,communicator,ierr)
        endif
     enddo
     enddo
  endif
  end subroutine reduce_mem_4D_real64
  !---------------------------
  subroutine reduce_mem_5D_real64(array, shared_memory)
  implicit none
  real(kind=real64), pointer,  intent(inout) :: array(:,:,:,:,:)
  logical,           optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, i, j, k, nsize
  logical :: do_reduce, use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory

  do_reduce = .false.
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array(:,:,1,1,1))
     do k=1, size(array(1,1,1,1,:))
     do j=1, size(array(1,1,1,:,1))
     do i=1, size(array(1,1,:,1,1))
        if (p_rank == 0) then
           call MPI_REDUCE(MPI_IN_PLACE,array(:,:,i,j,k),nsize,MPI_REAL8,MPI_SUM,0,communicator,ierr)
        else
           call MPI_REDUCE(array(:,:,i,j,k),           0,nsize,MPI_REAL8,MPI_SUM,0,communicator,ierr)
        endif
     enddo
     enddo
     enddo
  endif
  end subroutine reduce_mem_5D_real64
  !---------------------------
  subroutine get_window_1D_real64(arr, win)
  implicit none
  real(kind=real64), pointer, intent(in)  :: arr(:)
  integer,                    intent(out) :: win
  integer(MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr_arr, ptr
  integer     :: disp_unit, ierr, i

  ptr_arr   = C_LOC(arr)
  win       = MPI_WIN_NULL
  disp_unit = real64_size
  do i=1, num_windows
     if (windows(i) /= MPI_WIN_NULL) then
        !-- note that destroyed window causes a segmentation error in calling MPI_WIN_SHARED_QUERY.
        call MPI_WIN_SHARED_QUERY(windows(i), 0, wsize, disp_unit, ptr, ierr)
        if (C_ASSOCIATED(ptr, ptr_arr)) then
           win  = windows(i)
           exit
        endif
     endif
  enddo
  end subroutine get_window_1D_real64
  !---------------------------
  subroutine get_window_2D_real64(arr, win)
  implicit none
  real(kind=real64), pointer, intent(in)  :: arr(:,:)
  integer,                    intent(out) :: win
  integer(MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr_arr, ptr
  integer     :: disp_unit, ierr, i

  ptr_arr   = C_LOC(arr)
  win       = MPI_WIN_NULL
  disp_unit = real64_size
  do i=1, num_windows
     if (windows(i) /= MPI_WIN_NULL) then
        call MPI_WIN_SHARED_QUERY(windows(i), 0, wsize, disp_unit, ptr, ierr)
        if (C_ASSOCIATED(ptr, ptr_arr)) then
           win = windows(i)
           exit
        endif
     endif
  enddo
  end subroutine get_window_2D_real64
  !---------------------------
  subroutine get_window_3D_real64(arr, win)
  implicit none
  real(kind=real64), pointer, intent(in)  :: arr(:,:,:)
  integer,                    intent(out) :: win
  integer(MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr_arr, ptr
  integer     :: disp_unit, ierr, i

  ptr_arr   = C_LOC(arr)
  win       = MPI_WIN_NULL
  disp_unit = real64_size
  do i=1, num_windows
     if (windows(i) /= MPI_WIN_NULL) then
        call MPI_WIN_SHARED_QUERY(windows(i), 0, wsize, disp_unit, ptr, ierr)
        if (C_ASSOCIATED(ptr, ptr_arr)) then
           win = windows(i)
           exit
        endif
     endif
  enddo
  end subroutine get_window_3D_real64
  !---------------------------
  subroutine get_window_4D_real64(arr, win)
  implicit none
  real(kind=real64), pointer, intent(in)  :: arr(:,:,:,:)
  integer,                    intent(out) :: win
  integer(MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr_arr, ptr
  integer     :: disp_unit, ierr, i

  ptr_arr   = C_LOC(arr)
  win       = MPI_WIN_NULL
  disp_unit = real64_size
  do i=1, num_windows
     if (windows(i) /= MPI_WIN_NULL) then
        call MPI_WIN_SHARED_QUERY(windows(i), 0, wsize, disp_unit, ptr, ierr)
        if (C_ASSOCIATED(ptr, ptr_arr)) then
           win = windows(i)
           exit
        endif
     endif
  enddo
  end subroutine get_window_4D_real64
  !---------------------------
  subroutine get_window_5D_real64(arr, win)
  implicit none
  real(kind=real64), pointer, intent(in)  :: arr(:,:,:,:,:)
  integer,                    intent(out) :: win
  integer(MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr_arr, ptr
  integer     :: disp_unit, ierr, i

  ptr_arr   = C_LOC(arr)
  win       = MPI_WIN_NULL
  disp_unit = real64_size
  do i=1, num_windows
     if (windows(i) /= MPI_WIN_NULL) then
        call MPI_WIN_SHARED_QUERY(windows(i), 0, wsize, disp_unit, ptr, ierr)
        if (C_ASSOCIATED(ptr, ptr_arr)) then
           win = windows(i)
           exit
        endif
     endif
  enddo
  end subroutine get_window_5D_real64
  !---------------------------
  subroutine destroy_mem_1D_real64(arr)
  implicit none
  real(kind=real64), pointer, intent(inout) :: arr(:)
  integer :: win, ierr, i
  call get_window_1D_real64(arr, win)
  if (win /= MPI_WIN_NULL) then
     call destroy_shared_mem_window(win)
     nullify(arr)
  else
     if (associated(arr)) deallocate(arr)
  endif
  end subroutine destroy_mem_1D_real64
  !---------------------------
  subroutine destroy_mem_2D_real64(arr)
  implicit none
  real(kind=real64), pointer, intent(inout) :: arr(:,:)
  integer :: win
  call get_window_2D_real64(arr, win)
  if (win /= MPI_WIN_NULL) then
     call destroy_shared_mem_window(win)
     nullify(arr)
  else
     if (associated(arr)) deallocate(arr)
  endif
  end subroutine destroy_mem_2D_real64
  !---------------------------
  subroutine destroy_mem_3D_real64(arr)
  implicit none
  real(kind=real64), pointer, intent(inout) :: arr(:,:,:)
  integer :: win
  call get_window_3D_real64(arr, win)
  if (win /= MPI_WIN_NULL) then
     call destroy_shared_mem_window(win)
     nullify(arr)
  else
     if (associated(arr)) deallocate(arr)
  endif
  end subroutine destroy_mem_3D_real64
  !---------------------------
  subroutine destroy_mem_4D_real64(arr)
  implicit none
  real(kind=real64), pointer, intent(inout) :: arr(:,:,:,:)
  integer :: win
  call get_window_4D_real64(arr, win)
  if (win /= MPI_WIN_NULL) then
     call destroy_shared_mem_window(win)
     nullify(arr)
  else
     if (associated(arr)) deallocate(arr)
  endif
  end subroutine destroy_mem_4D_real64
  !---------------------------
  subroutine destroy_mem_5D_real64(arr)
  implicit none
  real(kind=real64), pointer, intent(inout) :: arr(:,:,:,:,:)
  integer :: win
  call get_window_5D_real64(arr, win)
  if (win /= MPI_WIN_NULL) then
     call destroy_shared_mem_window(win)
     nullify(arr)
  else
     if (associated(arr)) deallocate(arr)
  endif
  end subroutine destroy_mem_5D_real64
  !---------------------------
  subroutine create_shared_mem_1D_real32(array,arrshape)
  implicit none
  real(kind=real32), pointer, intent(inout) :: array(:)
  integer,                    intent(in)    :: arrshape(1)
  integer :: ierr
 
  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit
  
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(arrshape(1), MPI_ADDRESS_KIND)*real32_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif
 
  disp_unit   = real32_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0.0_real32
  endif
  call MPI_BARRIER(hostcomm, ierr)
  end subroutine create_shared_mem_1D_real32
  !---------------------------
  subroutine create_shared_mem_2D_real32(array,arrshape)
  implicit none
  real(kind=real32), pointer, intent(inout) :: array(:,:)
  integer,                    intent(in)    :: arrshape(2)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(int(arrshape(1),int64)*arrshape(2), MPI_ADDRESS_KIND)*real32_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = real32_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0.0_real32
  endif
  call MPI_BARRIER(hostcomm, ierr)
  end subroutine create_shared_mem_2D_real32
  !---------------------------
  subroutine create_shared_mem_3D_real32(array,arrshape)
  implicit none
  real(kind=real32), pointer, intent(inout) :: array(:,:,:)
  integer,                    intent(in)    :: arrshape(3)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(int(arrshape(1),int64)*arrshape(2)*arrshape(3), MPI_ADDRESS_KIND)*real32_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = real32_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0.0_real32
  endif
  call MPI_BARRIER(hostcomm, ierr)
  end subroutine create_shared_mem_3D_real32
  !---------------------------
  subroutine create_shared_mem_4D_real32(array,arrshape)
  implicit none
  real(kind=real32), pointer, intent(inout) :: array(:,:,:,:)
  integer,                    intent(in)    :: arrshape(4)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(int(arrshape(1),int64)*arrshape(2)*arrshape(3)*arrshape(4), MPI_ADDRESS_KIND)*real32_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = real32_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0.0_real32
  endif
  call MPI_BARRIER(hostcomm, ierr)
  end subroutine create_shared_mem_4D_real32
  !---------------------------
  subroutine create_shared_mem_5D_real32(array,arrshape)
  implicit none
  real(kind=real32), pointer, intent(inout) :: array(:,:,:,:,:)
  integer,                    intent(in)    :: arrshape(5)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(int(arrshape(1),int64)*arrshape(2)*arrshape(3)*arrshape(4)*arrshape(5), MPI_ADDRESS_KIND)*real32_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = real32_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0.0_real32
  endif
  call MPI_BARRIER(hostcomm, ierr)
  end subroutine create_shared_mem_5D_real32
  !---------------------------
  subroutine create_mem_1D_real32(array,arrshape,shared_memory)
  implicit none
  real(kind=real32), pointer, intent(inout) :: array(:)
  integer,                    intent(in)    :: arrshape(1)
  logical, optional,          intent(in)    :: shared_memory
  logical :: use_shared_memory
 
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_1D_real32(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1)))
     array = 0.0_real32
  endif
  end subroutine create_mem_1D_real32
  !---------------------------
  subroutine create_mem_2D_real32(array,arrshape,shared_memory)
  implicit none
  real(kind=real32), pointer, intent(inout) :: array(:,:)
  integer,                    intent(in)    :: arrshape(2)
  logical, optional,          intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_2D_real32(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2)))
     array = 0.0_real32
  endif
  end subroutine create_mem_2D_real32
  !---------------------------
  subroutine create_mem_3D_real32(array,arrshape,shared_memory)
  implicit none
  real(kind=real32), pointer, intent(inout) :: array(:,:,:)
  integer,                    intent(in)    :: arrshape(3)
  logical, optional,          intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_3D_real32(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2),arrshape(3)))
     array = 0.0_real32
  endif
  end subroutine create_mem_3D_real32
  !---------------------------
  subroutine create_mem_4D_real32(array,arrshape,shared_memory)
  implicit none
  real(kind=real32), pointer, intent(inout) :: array(:,:,:,:)
  integer,                    intent(in)    :: arrshape(4)
  logical, optional,          intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_4D_real32(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2),arrshape(3),arrshape(4)))
     array = 0.0_real32
  endif
  end subroutine create_mem_4D_real32
  !---------------------------
  subroutine create_mem_5D_real32(array,arrshape,shared_memory)
  implicit none
  real(kind=real32), pointer, intent(inout) :: array(:,:,:,:,:)
  integer,                    intent(in)    :: arrshape(5)
  logical, optional,          intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_5D_real32(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2),arrshape(3),arrshape(4),arrshape(5)))
     array = 0.0_real32
  endif
  end subroutine create_mem_5D_real32
  !---------------------------
  subroutine reduce_mem_1D_real32(array, shared_memory)
  implicit none
  real(kind=real32), pointer,  intent(inout) :: array(:)
  logical,           optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, nsize
  logical :: do_reduce, use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory

  do_reduce = .false.
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:),nsize,MPI_REAL4,MPI_SUM,0,communicator,ierr)
     else
        call MPI_REDUCE(array(:),           0,nsize,MPI_REAL4,MPI_SUM,0,communicator,ierr)
     endif
  endif
  end subroutine reduce_mem_1D_real32
  !---------------------------
  subroutine reduce_mem_2D_real32(array, shared_memory)
  implicit none
  real(kind=real32), pointer,  intent(inout) :: array(:,:)
  logical,           optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, nsize
  logical :: do_reduce, use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory

  do_reduce = .false.
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:,:),nsize,MPI_REAL4,MPI_SUM,0,communicator,ierr)
     else
        call MPI_REDUCE(array(:,:),           0,nsize,MPI_REAL4,MPI_SUM,0,communicator,ierr)
     endif
  endif
  end subroutine reduce_mem_2D_real32
  !---------------------------
  subroutine reduce_mem_3D_real32(array, shared_memory)
  implicit none
  real(kind=real32), pointer,  intent(inout) :: array(:,:,:)
  logical,           optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, k, nsize
  logical :: do_reduce, use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory

  do_reduce = .false.
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array(:,:,1))
     do k=1, size(array(1,1,:))
        if (p_rank == 0) then
           call MPI_REDUCE(MPI_IN_PLACE,array(:,:,k),nsize,MPI_REAL4,MPI_SUM,0,communicator,ierr)
        else
           call MPI_REDUCE(array(:,:,k),           0,nsize,MPI_REAL4,MPI_SUM,0,communicator,ierr)
        endif
     enddo
  endif
  end subroutine reduce_mem_3D_real32
  !---------------------------
  subroutine reduce_mem_4D_real32(array, shared_memory)
  implicit none
  real(kind=real32), pointer,  intent(inout) :: array(:,:,:,:)
  logical,           optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, j, k, nsize
  logical :: do_reduce, use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory

  do_reduce = .false.
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array(:,:,1,1))
     do k=1, size(array(1,1,1,:))
     do j=1, size(array(1,1,:,1))
        if (p_rank == 0) then
           call MPI_REDUCE(MPI_IN_PLACE,array(:,:,j,k),nsize,MPI_REAL4,MPI_SUM,0,communicator,ierr)
        else
           call MPI_REDUCE(array(:,:,j,k),           0,nsize,MPI_REAL4,MPI_SUM,0,communicator,ierr)
        endif
     enddo
     enddo
  endif
  end subroutine reduce_mem_4D_real32
  !---------------------------
  subroutine reduce_mem_5D_real32(array, shared_memory)
  implicit none
  real(kind=real32), pointer,  intent(inout) :: array(:,:,:,:,:)
  logical,           optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, i, j, k, nsize
  logical :: do_reduce, use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory

  do_reduce = .false.
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array(:,:,1,1,1))
     do k=1, size(array(1,1,1,1,:))
     do j=1, size(array(1,1,1,:,1))
     do i=1, size(array(1,1,:,1,1))
        if (p_rank == 0) then
           call MPI_REDUCE(MPI_IN_PLACE,array(:,:,i,j,k),nsize,MPI_REAL4,MPI_SUM,0,communicator,ierr)
        else
           call MPI_REDUCE(array(:,:,i,j,k),           0,nsize,MPI_REAL4,MPI_SUM,0,communicator,ierr)
        endif
     enddo
     enddo
     enddo
  endif
  end subroutine reduce_mem_5D_real32
  !---------------------------
  subroutine get_window_1D_real32(arr, win)
  implicit none
  real(kind=real32), pointer, intent(in)  :: arr(:)
  integer,                    intent(out) :: win
  integer(MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr_arr, ptr
  integer     :: disp_unit, ierr, i

  ptr_arr   = C_LOC(arr)
  win       = MPI_WIN_NULL
  disp_unit = real32_size
  do i=1, num_windows
     if (windows(i) /= MPI_WIN_NULL) then
        !-- note that destroyed window causes a segmentation error in calling MPI_WIN_SHARED_QUERY.
        call MPI_WIN_SHARED_QUERY(windows(i), 0, wsize, disp_unit, ptr, ierr)
        if (C_ASSOCIATED(ptr, ptr_arr)) then
           win  = windows(i)
           exit
        endif
     endif
  enddo
  end subroutine get_window_1D_real32
  !---------------------------
  subroutine get_window_2D_real32(arr, win)
  implicit none
  real(kind=real32), pointer, intent(in)  :: arr(:,:)
  integer,                    intent(out) :: win
  integer(MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr_arr, ptr
  integer     :: disp_unit, ierr, i

  ptr_arr   = C_LOC(arr)
  win       = MPI_WIN_NULL
  disp_unit = real32_size
  do i=1, num_windows
     if (windows(i) /= MPI_WIN_NULL) then
        call MPI_WIN_SHARED_QUERY(windows(i), 0, wsize, disp_unit, ptr, ierr)
        if (C_ASSOCIATED(ptr, ptr_arr)) then
           win = windows(i)
           exit
        endif
     endif
  enddo
  end subroutine get_window_2D_real32
  !---------------------------
  subroutine get_window_3D_real32(arr, win)
  implicit none
  real(kind=real32), pointer, intent(in)  :: arr(:,:,:)
  integer,                    intent(out) :: win
  integer(MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr_arr, ptr
  integer     :: disp_unit, ierr, i

  ptr_arr   = C_LOC(arr)
  win       = MPI_WIN_NULL
  disp_unit = real32_size
  do i=1, num_windows
     if (windows(i) /= MPI_WIN_NULL) then
        call MPI_WIN_SHARED_QUERY(windows(i), 0, wsize, disp_unit, ptr, ierr)
        if (C_ASSOCIATED(ptr, ptr_arr)) then
           win = windows(i)
           exit
        endif
     endif
  enddo
  end subroutine get_window_3D_real32
  !---------------------------
  subroutine get_window_4D_real32(arr, win)
  implicit none
  real(kind=real32), pointer, intent(in)  :: arr(:,:,:,:)
  integer,                    intent(out) :: win
  integer(MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr_arr, ptr
  integer     :: disp_unit, ierr, i

  ptr_arr   = C_LOC(arr)
  win       = MPI_WIN_NULL
  disp_unit = real32_size
  do i=1, num_windows
     if (windows(i) /= MPI_WIN_NULL) then
        call MPI_WIN_SHARED_QUERY(windows(i), 0, wsize, disp_unit, ptr, ierr)
        if (C_ASSOCIATED(ptr, ptr_arr)) then
           win = windows(i)
           exit
        endif
     endif
  enddo
  end subroutine get_window_4D_real32
  !---------------------------
  subroutine get_window_5D_real32(arr, win)
  implicit none
  real(kind=real32), pointer, intent(in)  :: arr(:,:,:,:,:)
  integer,                    intent(out) :: win
  integer(MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr_arr, ptr
  integer     :: disp_unit, ierr, i

  ptr_arr   = C_LOC(arr)
  win       = MPI_WIN_NULL
  disp_unit = real32_size
  do i=1, num_windows
     if (windows(i) /= MPI_WIN_NULL) then
        call MPI_WIN_SHARED_QUERY(windows(i), 0, wsize, disp_unit, ptr, ierr)
        if (C_ASSOCIATED(ptr, ptr_arr)) then
           win = windows(i)
           exit
        endif
     endif
  enddo
  end subroutine get_window_5D_real32
  !---------------------------
  subroutine destroy_mem_1D_real32(arr)
  implicit none
  real(kind=real32), pointer, intent(inout) :: arr(:)
  integer :: win, ierr, i
  call get_window_1D_real32(arr, win)
  if (win /= MPI_WIN_NULL) then
     call destroy_shared_mem_window(win)
     nullify(arr)
  else
     if (associated(arr)) deallocate(arr)
  endif
  end subroutine destroy_mem_1D_real32
  !---------------------------
  subroutine destroy_mem_2D_real32(arr)
  implicit none
  real(kind=real32), pointer, intent(inout) :: arr(:,:)
  integer :: win
  call get_window_2D_real32(arr, win)
  if (win /= MPI_WIN_NULL) then
     call destroy_shared_mem_window(win)
     nullify(arr)
  else
     if (associated(arr)) deallocate(arr)
  endif
  end subroutine destroy_mem_2D_real32
  !---------------------------
  subroutine destroy_mem_3D_real32(arr)
  implicit none
  real(kind=real32), pointer, intent(inout) :: arr(:,:,:)
  integer :: win
  call get_window_3D_real32(arr, win)
  if (win /= MPI_WIN_NULL) then
     call destroy_shared_mem_window(win)
     nullify(arr)
  else
     if (associated(arr)) deallocate(arr)
  endif
  end subroutine destroy_mem_3D_real32
  !---------------------------
  subroutine destroy_mem_4D_real32(arr)
  implicit none
  real(kind=real32), pointer, intent(inout) :: arr(:,:,:,:)
  integer :: win
  call get_window_4D_real32(arr, win)
  if (win /= MPI_WIN_NULL) then
     call destroy_shared_mem_window(win)
     nullify(arr)
  else
     if (associated(arr)) deallocate(arr)
  endif
  end subroutine destroy_mem_4D_real32
  !---------------------------
  subroutine destroy_mem_5D_real32(arr)
  implicit none
  real(kind=real32), pointer, intent(inout) :: arr(:,:,:,:,:)
  integer :: win
  call get_window_5D_real32(arr, win)
  if (win /= MPI_WIN_NULL) then
     call destroy_shared_mem_window(win)
     nullify(arr)
  else
     if (associated(arr)) deallocate(arr)
  endif
  end subroutine destroy_mem_5D_real32
  !---------------------------
  subroutine create_shared_mem_1D_int64(array,arrshape)
  implicit none
  integer(int64), pointer, intent(inout) :: array(:)
  integer,                 intent(in)    :: arrshape(1)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(arrshape(1), MPI_ADDRESS_KIND)*int64_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = int64_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0_int64
  endif
  call MPI_BARRIER(hostcomm, ierr)
  end subroutine create_shared_mem_1D_int64
  !---------------------------
  subroutine create_shared_mem_2D_int64(array,arrshape)
  implicit none
  integer(int64), pointer, intent(inout) :: array(:,:)
  integer,                 intent(in)    :: arrshape(2)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(int(arrshape(1),int64)*arrshape(2), MPI_ADDRESS_KIND)*int64_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = int64_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0_int64
  endif
  call MPI_BARRIER(hostcomm, ierr)
  end subroutine create_shared_mem_2D_int64
  !---------------------------
  subroutine create_shared_mem_3D_int64(array,arrshape)
  implicit none
  integer(int64), pointer, intent(inout) :: array(:,:,:)
  integer,                 intent(in)    :: arrshape(3)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(int(arrshape(1),int64)*arrshape(2)*arrshape(3), MPI_ADDRESS_KIND)*int64_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = int64_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0_int64
  endif
  call MPI_BARRIER(hostcomm, ierr)
  end subroutine create_shared_mem_3D_int64
  !---------------------------
  subroutine create_mem_1D_int64(array,arrshape,shared_memory)
  implicit none
  integer(int64), pointer,  intent(inout) :: array(:)
  integer,                  intent(in)    :: arrshape(1)
  logical,        optional, intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_1D_int64(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1)))
     array = 0_int64
  endif
  end subroutine create_mem_1D_int64
  !---------------------------
  subroutine create_mem_2D_int64(array,arrshape,shared_memory)
  implicit none
  integer(int64), pointer,  intent(inout) :: array(:,:)
  integer,                  intent(in)    :: arrshape(2)
  logical,        optional, intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_2D_int64(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2)))
     array = 0_int64
  endif
  end subroutine create_mem_2D_int64
  !---------------------------
  subroutine create_mem_3D_int64(array,arrshape,shared_memory)
  implicit none
  integer(int64), pointer,  intent(inout) :: array(:,:,:)
  integer,                  intent(in)    :: arrshape(3)
  logical,        optional, intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_3D_int64(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2),arrshape(3)))
     array = 0_int64
  endif
  end subroutine create_mem_3D_int64
  !---------------------------
  subroutine reduce_mem_1D_int64(array, shared_memory)
  implicit none
  integer(int64), pointer,  intent(inout) :: array(:)
  logical,        optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, nsize
  logical :: do_reduce, use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory

  do_reduce = .false.
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:),nsize,MPI_INTEGER8,MPI_SUM,0,communicator,ierr)
     else
        call MPI_REDUCE(array(:),           0,nsize,MPI_INTEGER8,MPI_SUM,0,communicator,ierr)
     endif
  endif
  end subroutine reduce_mem_1D_int64
  !---------------------------
  subroutine reduce_mem_2D_int64(array, shared_memory)
  implicit none
  integer(int64), pointer,  intent(inout) :: array(:,:)
  logical,        optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, nsize
  logical :: do_reduce, use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory

  do_reduce = .false.
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:,:),nsize,MPI_INTEGER8,MPI_SUM,0,communicator,ierr)
     else
        call MPI_REDUCE(array(:,:),           0,nsize,MPI_INTEGER8,MPI_SUM,0,communicator,ierr)
     endif
  endif
  end subroutine reduce_mem_2D_int64
  !---------------------------
  subroutine reduce_mem_3D_int64(array, shared_memory)
  implicit none
  integer(int64), pointer,  intent(inout) :: array(:,:,:)
  logical,        optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, k, nsize
  logical :: do_reduce, use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory

  do_reduce = .false.
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array(:,:,1))
     do k=1, size(array(1,1,:))
        if (p_rank == 0) then
           call MPI_REDUCE(MPI_IN_PLACE,array(:,:,k),nsize,MPI_INTEGER8,MPI_SUM,0,communicator,ierr)
        else
           call MPI_REDUCE(array(:,:,k),           0,nsize,MPI_INTEGER8,MPI_SUM,0,communicator,ierr)
        endif
     enddo
  endif
  end subroutine reduce_mem_3D_int64
  !---------------------------
  subroutine get_window_1D_int64(arr, win)
  implicit none
  integer(int64), pointer, intent(in)  :: arr(:)
  integer,                 intent(out) :: win
  integer(MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr_arr, ptr
  integer     :: disp_unit, ierr, i

  ptr_arr   = C_LOC(arr)
  win       = MPI_WIN_NULL
  disp_unit = int64_size
  do i=1, num_windows
     if (windows(i) /= MPI_WIN_NULL) then
        call MPI_WIN_SHARED_QUERY(windows(i), 0, wsize, disp_unit, ptr, ierr)
        if (C_ASSOCIATED(ptr, ptr_arr)) then
           win  = windows(i)
           exit
        endif
     endif
  enddo
  end subroutine get_window_1D_int64
  !---------------------------
  subroutine get_window_2D_int64(arr, win)
  implicit none
  integer(int64), pointer, intent(in)  :: arr(:,:)
  integer,                 intent(out) :: win
  integer(MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr_arr, ptr
  integer     :: disp_unit, ierr, i

  ptr_arr   = C_LOC(arr)
  win       = MPI_WIN_NULL
  disp_unit = int64_size
  do i=1, num_windows
     if (windows(i) /= MPI_WIN_NULL) then
        call MPI_WIN_SHARED_QUERY(windows(i), 0, wsize, disp_unit, ptr, ierr)
        if (C_ASSOCIATED(ptr, ptr_arr)) then
           win  = windows(i)
           exit
        endif
     endif
  enddo
  end subroutine get_window_2D_int64
  !---------------------------
  subroutine get_window_3D_int64(arr, win)
  implicit none
  integer(int64), pointer, intent(in)  :: arr(:,:,:)
  integer,                 intent(out) :: win
  integer(MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr_arr, ptr
  integer     :: disp_unit, ierr, i

  ptr_arr   = C_LOC(arr)
  win       = MPI_WIN_NULL
  disp_unit = int64_size
  do i=1, num_windows
     if (windows(i) /= MPI_WIN_NULL) then
        call MPI_WIN_SHARED_QUERY(windows(i), 0, wsize, disp_unit, ptr, ierr)
        if (C_ASSOCIATED(ptr, ptr_arr)) then
           win  = windows(i)
           exit
        endif
     endif
  enddo
  end subroutine get_window_3D_int64
  !---------------------------
  subroutine destroy_mem_1D_int64(arr)
  implicit none
  integer(int64), pointer, intent(inout) :: arr(:)
  integer :: win, ierr
  call get_window_1D_int64(arr, win)
  if (win /= MPI_WIN_NULL) then
     call destroy_shared_mem_window(win)
     nullify(arr)
  else
     if (associated(arr)) deallocate(arr)
  endif
  end subroutine destroy_mem_1D_int64
  !---------------------------
  subroutine destroy_mem_2D_int64(arr)
  implicit none
  integer(int64), pointer, intent(inout) :: arr(:,:)
  integer :: win, ierr
  call get_window_2D_int64(arr, win)
  if (win /= MPI_WIN_NULL) then
     call destroy_shared_mem_window(win)
     nullify(arr)
  else
     if (associated(arr)) deallocate(arr)
  endif
  end subroutine destroy_mem_2D_int64
  !---------------------------
  subroutine destroy_mem_3D_int64(arr)
  implicit none
  integer(int64), pointer, intent(inout) :: arr(:,:,:)
  integer :: win, ierr
  call get_window_3D_int64(arr, win)
  if (win /= MPI_WIN_NULL) then
     call destroy_shared_mem_window(win)
     nullify(arr)
  else
     if (associated(arr)) deallocate(arr)
  endif
  end subroutine destroy_mem_3D_int64
  !---------------------------
  subroutine create_shared_mem_1D_int32(array,arrshape)
  implicit none
  integer(int32), pointer, intent(inout) :: array(:)
  integer,                 intent(in)    :: arrshape(1)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(arrshape(1), MPI_ADDRESS_KIND)*int32_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = int32_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0_int32
  endif
  call MPI_BARRIER(hostcomm, ierr)
  end subroutine create_shared_mem_1D_int32
  !---------------------------
  subroutine create_shared_mem_2D_int32(array,arrshape)
  implicit none
  integer(int32), pointer, intent(inout) :: array(:,:)
  integer,                 intent(in)    :: arrshape(2)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(int(arrshape(1),int64)*arrshape(2), MPI_ADDRESS_KIND)*int32_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = int32_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0_int32
  endif
  call MPI_BARRIER(hostcomm, ierr)
  end subroutine create_shared_mem_2D_int32
  !---------------------------
  subroutine create_shared_mem_3D_int32(array,arrshape)
  implicit none
  integer(int32), pointer, intent(inout) :: array(:,:,:)
  integer,                 intent(in)    :: arrshape(3)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(int(arrshape(1),int64)*arrshape(2)*arrshape(3), MPI_ADDRESS_KIND)*int32_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = int32_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0_int32
  endif
  call MPI_BARRIER(hostcomm, ierr)
  end subroutine create_shared_mem_3D_int32
  !---------------------------
  subroutine create_mem_1D_int32(array,arrshape,shared_memory)
  implicit none
  integer(int32), pointer,  intent(inout) :: array(:)
  integer,                  intent(in)    :: arrshape(1)
  logical,        optional, intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_1D_int32(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1)))
     array = 0_int32
  endif
  end subroutine create_mem_1D_int32
  !---------------------------
  subroutine create_mem_2D_int32(array,arrshape,shared_memory)
  implicit none
  integer(int32), pointer,  intent(inout) :: array(:,:)
  integer,                  intent(in)    :: arrshape(2)
  logical,        optional, intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_2D_int32(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2)))
     array = 0_int32
  endif
  end subroutine create_mem_2D_int32
  !---------------------------
  subroutine create_mem_3D_int32(array,arrshape,shared_memory)
  implicit none
  integer(int32), pointer,  intent(inout) :: array(:,:,:)
  integer,                  intent(in)    :: arrshape(3)
  logical,        optional, intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_3D_int32(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2),arrshape(3)))
     array = 0_int32
  endif
  end subroutine create_mem_3D_int32
  !---------------------------
  subroutine reduce_mem_1D_int32(array, shared_memory)
  implicit none
  integer(int32), pointer,  intent(inout) :: array(:)
  logical,        optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, nsize
  logical :: do_reduce, use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory

  do_reduce = .false.
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:),nsize,MPI_INTEGER4,MPI_SUM,0,communicator,ierr)
     else
        call MPI_REDUCE(array(:),           0,nsize,MPI_INTEGER4,MPI_SUM,0,communicator,ierr)
     endif
  endif
  end subroutine reduce_mem_1D_int32
  !---------------------------
  subroutine reduce_mem_2D_int32(array, shared_memory)
  implicit none
  integer(int32), pointer,  intent(inout) :: array(:,:)
  logical,        optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, nsize
  logical :: do_reduce, use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory

  do_reduce = .false.
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:,:),nsize,MPI_INTEGER4,MPI_SUM,0,communicator,ierr)
     else
        call MPI_REDUCE(array(:,:),           0,nsize,MPI_INTEGER4,MPI_SUM,0,communicator,ierr)
     endif
  endif
  end subroutine reduce_mem_2D_int32
  !---------------------------
  subroutine reduce_mem_3D_int32(array, shared_memory)
  implicit none
  integer(int32), pointer,  intent(inout) :: array(:,:,:)
  logical,        optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, k, nsize
  logical :: do_reduce, use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory

  do_reduce = .false.
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array(:,:,1))
     do k=1, size(array(1,1,:))
        if (p_rank == 0) then
           call MPI_REDUCE(MPI_IN_PLACE,array(:,:,k),nsize,MPI_INTEGER4,MPI_SUM,0,communicator,ierr)
        else
           call MPI_REDUCE(array(:,:,k),           0,nsize,MPI_INTEGER4,MPI_SUM,0,communicator,ierr)
        endif
     enddo
  endif
  end subroutine reduce_mem_3D_int32
  !---------------------------
  subroutine get_window_1D_int32(arr, win)
  implicit none
  integer(int32), pointer, intent(in)  :: arr(:)
  integer,                 intent(out) :: win
  integer(MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr_arr, ptr
  integer     :: disp_unit, ierr, i

  ptr_arr   = C_LOC(arr)
  win       = MPI_WIN_NULL
  disp_unit = int32_size
  do i=1, num_windows
     if (windows(i) /= MPI_WIN_NULL) then
        call MPI_WIN_SHARED_QUERY(windows(i), 0, wsize, disp_unit, ptr, ierr)
        if (C_ASSOCIATED(ptr, ptr_arr)) then
           win  = windows(i)
           exit
        endif
     endif
  enddo
  end subroutine get_window_1D_int32
  !---------------------------
  subroutine get_window_2D_int32(arr, win)
  implicit none
  integer(int32), pointer, intent(in)  :: arr(:,:)
  integer,                 intent(out) :: win
  integer(MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr_arr, ptr
  integer     :: disp_unit, ierr, i

  ptr_arr   = C_LOC(arr)
  win       = MPI_WIN_NULL
  disp_unit = int32_size
  do i=1, num_windows
     if (windows(i) /= MPI_WIN_NULL) then
        call MPI_WIN_SHARED_QUERY(windows(i), 0, wsize, disp_unit, ptr, ierr)
        if (C_ASSOCIATED(ptr, ptr_arr)) then
           win  = windows(i)
           exit
        endif
     endif
  enddo
  end subroutine get_window_2D_int32
  !---------------------------
  subroutine get_window_3D_int32(arr, win)
  implicit none
  integer(int32), pointer, intent(in)  :: arr(:,:,:)
  integer,                 intent(out) :: win
  integer(MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr_arr, ptr
  integer     :: disp_unit, ierr, i

  ptr_arr   = C_LOC(arr)
  win       = MPI_WIN_NULL
  disp_unit = int32_size
  do i=1, num_windows
     if (windows(i) /= MPI_WIN_NULL) then
        call MPI_WIN_SHARED_QUERY(windows(i), 0, wsize, disp_unit, ptr, ierr)
        if (C_ASSOCIATED(ptr, ptr_arr)) then
           win  = windows(i)
           exit
        endif
     endif
  enddo
  end subroutine get_window_3D_int32
  !---------------------------
  subroutine destroy_mem_1D_int32(arr)
  implicit none
  integer(int32), pointer, intent(inout) :: arr(:)
  integer :: win, ierr
  call get_window_1D_int32(arr, win)
  if (win /= MPI_WIN_NULL) then
     call destroy_shared_mem_window(win)
     nullify(arr)
  else
     if (associated(arr)) deallocate(arr)
  endif
  end subroutine destroy_mem_1D_int32
  !---------------------------
  subroutine destroy_mem_2D_int32(arr)
  implicit none
  integer(int32), pointer, intent(inout) :: arr(:,:)
  integer :: win, ierr
  call get_window_2D_int32(arr, win)
  if (win /= MPI_WIN_NULL) then
     call destroy_shared_mem_window(win)
     nullify(arr)
  else
     if (associated(arr)) deallocate(arr)
  endif
  end subroutine destroy_mem_2D_int32
  !---------------------------
  subroutine destroy_mem_3D_int32(arr)
  implicit none
  integer(int32), pointer, intent(inout) :: arr(:,:,:)
  integer :: win, ierr
  call get_window_3D_int32(arr, win)
  if (win /= MPI_WIN_NULL) then
     call destroy_shared_mem_window(win)
     nullify(arr)
  else
     if (associated(arr)) deallocate(arr)
  endif
  end subroutine destroy_mem_3D_int32
  !---------------------------
  subroutine create_shared_mem_1D_int16(array,arrshape)
  implicit none
  integer(int16), pointer, intent(inout) :: array(:)
  integer,                 intent(in)    :: arrshape(1)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(arrshape(1), MPI_ADDRESS_KIND)*int16_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = int16_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0_int16
  endif
  call MPI_BARRIER(hostcomm, ierr)
  end subroutine create_shared_mem_1D_int16
  !---------------------------
  subroutine create_shared_mem_2D_int16(array,arrshape)
  implicit none
  integer(int16), pointer, intent(inout) :: array(:,:)
  integer,                 intent(in)    :: arrshape(2)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(int(arrshape(1),int64)*arrshape(2), MPI_ADDRESS_KIND)*int16_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = int16_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0_int16
  endif
  call MPI_BARRIER(hostcomm, ierr)
  end subroutine create_shared_mem_2D_int16
  !---------------------------
  subroutine create_shared_mem_3D_int16(array,arrshape)
  implicit none
  integer(int16), pointer, intent(inout) :: array(:,:,:)
  integer,                 intent(in)    :: arrshape(3)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(int(arrshape(1),int64)*arrshape(2)*arrshape(3), MPI_ADDRESS_KIND)*int16_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = int16_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0_int16
  endif
  call MPI_BARRIER(hostcomm, ierr)
  end subroutine create_shared_mem_3D_int16
  !---------------------------
  subroutine create_mem_1D_int16(array,arrshape,shared_memory)
  implicit none
  integer(int16), pointer,  intent(inout) :: array(:)
  integer,                  intent(in)    :: arrshape(1)
  logical,        optional, intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_1D_int16(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1)))
     array = 0_int16
  endif
  end subroutine create_mem_1D_int16
  !---------------------------
  subroutine create_mem_2D_int16(array,arrshape,shared_memory)
  implicit none
  integer(int16), pointer,  intent(inout) :: array(:,:)
  integer,                  intent(in)    :: arrshape(2)
  logical,        optional, intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_2D_int16(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2)))
     array = 0_int16
  endif
  end subroutine create_mem_2D_int16
  !---------------------------
  subroutine create_mem_3D_int16(array,arrshape,shared_memory)
  implicit none
  integer(int16), pointer,  intent(inout) :: array(:,:,:)
  integer,                  intent(in)    :: arrshape(3)
  logical,        optional, intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_3D_int16(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2),arrshape(3)))
     array = 0_int16
  endif
  end subroutine create_mem_3D_int16
  !---------------------------
  subroutine reduce_mem_1D_int16(array, shared_memory)
  implicit none
  integer(int16), pointer,  intent(inout) :: array(:)
  logical,        optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, nsize
  logical :: do_reduce, use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory

  do_reduce = .false.
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:),nsize,MPI_INTEGER2,MPI_SUM,0,communicator,ierr)
     else
        call MPI_REDUCE(array(:),           0,nsize,MPI_INTEGER2,MPI_SUM,0,communicator,ierr)
     endif
  endif
  end subroutine reduce_mem_1D_int16
  !---------------------------
  subroutine reduce_mem_2D_int16(array, shared_memory)
  implicit none
  integer(int16), pointer,  intent(inout) :: array(:,:)
  logical,        optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, nsize
  logical :: do_reduce, use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory

  do_reduce = .false.
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:,:),nsize,MPI_INTEGER2,MPI_SUM,0,communicator,ierr)
     else
        call MPI_REDUCE(array(:,:),           0,nsize,MPI_INTEGER2,MPI_SUM,0,communicator,ierr)
     endif
  endif
  end subroutine reduce_mem_2D_int16
  !---------------------------
  subroutine reduce_mem_3D_int16(array, shared_memory)
  implicit none
  integer(int16), pointer,  intent(inout) :: array(:,:,:)
  logical,        optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, k, nsize
  logical :: do_reduce, use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory

  do_reduce = .false.
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array(:,:,1))
     do k=1, size(array(1,1,:))
        if (p_rank == 0) then
           call MPI_REDUCE(MPI_IN_PLACE,array(:,:,k),nsize,MPI_INTEGER2,MPI_SUM,0,communicator,ierr)
        else
           call MPI_REDUCE(array(:,:,k),           0,nsize,MPI_INTEGER2,MPI_SUM,0,communicator,ierr)
        endif
     enddo
  endif
  end subroutine reduce_mem_3D_int16
  !---------------------------
  subroutine get_window_1D_int16(arr, win)
  implicit none
  integer(int16), pointer, intent(in)  :: arr(:)
  integer,                 intent(out) :: win
  integer(MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr_arr, ptr
  integer     :: disp_unit, ierr, i

  ptr_arr   = C_LOC(arr)
  win       = MPI_WIN_NULL
  disp_unit = int16_size
  do i=1, num_windows
     if (windows(i) /= MPI_WIN_NULL) then
        call MPI_WIN_SHARED_QUERY(windows(i), 0, wsize, disp_unit, ptr, ierr)
        if (C_ASSOCIATED(ptr, ptr_arr)) then
           win  = windows(i)
           exit
        endif
     endif
  enddo
  end subroutine get_window_1D_int16
  !---------------------------
  subroutine get_window_2D_int16(arr, win)
  implicit none
  integer(int16), pointer, intent(in)  :: arr(:,:)
  integer,                 intent(out) :: win
  integer(MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr_arr, ptr
  integer     :: disp_unit, ierr, i

  ptr_arr   = C_LOC(arr)
  win       = MPI_WIN_NULL
  disp_unit = int16_size
  do i=1, num_windows
     if (windows(i) /= MPI_WIN_NULL) then
        call MPI_WIN_SHARED_QUERY(windows(i), 0, wsize, disp_unit, ptr, ierr)
        if (C_ASSOCIATED(ptr, ptr_arr)) then
           win  = windows(i)
           exit
        endif
     endif
  enddo
  end subroutine get_window_2D_int16
  !---------------------------
  subroutine get_window_3D_int16(arr, win)
  implicit none
  integer(int16), pointer, intent(in)  :: arr(:,:,:)
  integer,                 intent(out) :: win
  integer(MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr_arr, ptr
  integer     :: disp_unit, ierr, i

  ptr_arr   = C_LOC(arr)
  win       = MPI_WIN_NULL
  disp_unit = int16_size
  do i=1, num_windows
     if (windows(i) /= MPI_WIN_NULL) then
        call MPI_WIN_SHARED_QUERY(windows(i), 0, wsize, disp_unit, ptr, ierr)
        if (C_ASSOCIATED(ptr, ptr_arr)) then
           win  = windows(i)
           exit
        endif
     endif
  enddo
  end subroutine get_window_3D_int16
  !---------------------------
  subroutine destroy_mem_1D_int16(arr)
  implicit none
  integer(int16), pointer, intent(inout) :: arr(:)
  integer :: win, ierr
  call get_window_1D_int16(arr, win)
  if (win /= MPI_WIN_NULL) then
     call destroy_shared_mem_window(win)
     nullify(arr)
  else
     if (associated(arr)) deallocate(arr)
  endif
  end subroutine destroy_mem_1D_int16
  !---------------------------
  subroutine destroy_mem_2D_int16(arr)
  implicit none
  integer(int16), pointer, intent(inout) :: arr(:,:)
  integer :: win, ierr
  call get_window_2D_int16(arr, win)
  if (win /= MPI_WIN_NULL) then
     call destroy_shared_mem_window(win)
     nullify(arr)
  else
     if (associated(arr)) deallocate(arr)
  endif
  end subroutine destroy_mem_2D_int16
  !---------------------------
  subroutine destroy_mem_3D_int16(arr)
  implicit none
  integer(int16), pointer, intent(inout) :: arr(:,:,:)
  integer :: win, ierr
  call get_window_3D_int16(arr, win)
  if (win /= MPI_WIN_NULL) then
     call destroy_shared_mem_window(win)
     nullify(arr)
  else
     if (associated(arr)) deallocate(arr)
  endif
  end subroutine destroy_mem_3D_int16
  !---------------------------
  subroutine create_shared_mem_1D_int8(array,arrshape)
  implicit none
  integer(int8), pointer, intent(inout) :: array(:)
  integer,                intent(in)    :: arrshape(1)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(arrshape(1), MPI_ADDRESS_KIND)*int8_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = int8_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0_int8
  endif
  call MPI_BARRIER(hostcomm, ierr)
  end subroutine create_shared_mem_1D_int8
  !---------------------------
  subroutine create_shared_mem_2D_int8(array,arrshape)
  implicit none
  integer(int8), pointer, intent(inout) :: array(:,:)
  integer,                intent(in)    :: arrshape(2)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(int(arrshape(1),int64)*arrshape(2), MPI_ADDRESS_KIND)*int8_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = int8_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0_int8
  endif
  call MPI_BARRIER(hostcomm, ierr)
  end subroutine create_shared_mem_2D_int8
  !---------------------------
  subroutine create_shared_mem_3D_int8(array,arrshape)
  implicit none
  integer(int8), pointer, intent(inout) :: array(:,:,:)
  integer,                intent(in)    :: arrshape(3)
  integer :: ierr

  integer(kind=MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr
  integer     :: disp_unit

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (h_rank == 0) then
     wsize = int(int(arrshape(1),int64)*arrshape(2)*arrshape(3), MPI_ADDRESS_KIND)*int8_size
  else
     wsize = 0_MPI_ADDRESS_KIND
  endif

  disp_unit   = int8_size
  num_windows = num_windows + 1
  if (num_windows > MAX_NUM_WINDOWS) then
     stop 'Please increase MAX_NUM_WINDOWS'
  endif
  call MPI_WIN_ALLOCATE_SHARED(wsize, disp_unit, MPI_INFO_NULL, hostcomm, ptr, windows(num_windows), ierr)
  if (h_rank /= 0) then
     call MPI_WIN_SHARED_QUERY(windows(num_windows), 0, wsize, disp_unit, ptr, ierr)
  endif
  call C_F_POINTER(ptr, array, arrshape)
  call MPI_WIN_FENCE(0, windows(num_windows), ierr)

  if (h_rank == 0) then
     array = 0_int8
  endif
  call MPI_BARRIER(hostcomm, ierr)
  end subroutine create_shared_mem_3D_int8
  !---------------------------
  subroutine create_mem_1D_int8(array,arrshape,shared_memory)
  implicit none
  integer(int8), pointer,  intent(inout) :: array(:)
  integer,                 intent(in)    :: arrshape(1)
  logical,       optional, intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_1D_int8(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1)))
     array = 0_int8
  endif
  end subroutine create_mem_1D_int8
  !---------------------------
  subroutine create_mem_2D_int8(array,arrshape,shared_memory)
  implicit none
  integer(int8), pointer,  intent(inout) :: array(:,:)
  integer,                 intent(in)    :: arrshape(2)
  logical,       optional, intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_2D_int8(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2)))
     array = 0_int8
  endif
  end subroutine create_mem_2D_int8
  !---------------------------
  subroutine create_mem_3D_int8(array,arrshape,shared_memory)
  implicit none
  integer(int8), pointer,  intent(inout) :: array(:,:,:)
  integer,                 intent(in)    :: arrshape(3)
  logical,       optional, intent(in)    :: shared_memory
  logical :: use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory
  if (use_shared_memory) then
     call create_shared_mem_3D_int8(array,arrshape)
  else
     if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2),arrshape(3)))
     array = 0_int8
  endif
  end subroutine create_mem_3D_int8
  !---------------------------
  subroutine reduce_mem_1D_int8(array, shared_memory)
  implicit none
  integer(int8), pointer,  intent(inout) :: array(:)
  logical,       optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, nsize
  logical :: do_reduce, use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory

  do_reduce = .false.
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:),nsize,MPI_INTEGER1,MPI_SUM,0,communicator,ierr)
     else
        call MPI_REDUCE(array(:),           0,nsize,MPI_INTEGER1,MPI_SUM,0,communicator,ierr)
     endif
  endif
  end subroutine reduce_mem_1D_int8
  !---------------------------
  subroutine reduce_mem_2D_int8(array, shared_memory)
  implicit none
  integer(int8), pointer,  intent(inout) :: array(:,:)
  logical,       optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, nsize
  logical :: do_reduce, use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory

  do_reduce = .false.
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:,:),nsize,MPI_INTEGER1,MPI_SUM,0,communicator,ierr)
     else
        call MPI_REDUCE(array(:,:),           0,nsize,MPI_INTEGER1,MPI_SUM,0,communicator,ierr)
     endif
  endif
  end subroutine reduce_mem_2D_int8
  !---------------------------
  subroutine reduce_mem_3D_int8(array, shared_memory)
  implicit none
  integer(int8), pointer,  intent(inout) :: array(:,:,:)
  logical,       optional, intent(in)    :: shared_memory
  integer :: communicator, ierr, k, nsize
  logical :: do_reduce, use_shared_memory

  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  use_shared_memory = .false.
  if (present(shared_memory)) use_shared_memory = shared_memory

  do_reduce = .false.
  if (use_shared_memory) then
     communicator = same_hrank_comm
     if (h_rank == 0 .and. num_hosts > 1) do_reduce = .true.
  else
     communicator = MPI_COMM_WORLD
     do_reduce = .true.
  endif
  if (do_reduce) then
     nsize = size(array(:,:,1))
     do k=1, size(array(1,1,:))
        if (p_rank == 0) then
           call MPI_REDUCE(MPI_IN_PLACE,array(:,:,k),nsize,MPI_INTEGER1,MPI_SUM,0,communicator,ierr)
        else
           call MPI_REDUCE(array(:,:,k),           0,nsize,MPI_INTEGER1,MPI_SUM,0,communicator,ierr)
        endif
     enddo
  endif
  end subroutine reduce_mem_3D_int8
  !---------------------------
  subroutine get_window_1D_int8(arr, win)
  implicit none
  integer(int8), pointer, intent(in)  :: arr(:)
  integer,                intent(out) :: win
  integer(MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr_arr, ptr
  integer     :: disp_unit, ierr, i

  ptr_arr   = C_LOC(arr)
  win       = MPI_WIN_NULL
  disp_unit = int8_size
  do i=1, num_windows
     if (windows(i) /= MPI_WIN_NULL) then
        call MPI_WIN_SHARED_QUERY(windows(i), 0, wsize, disp_unit, ptr, ierr)
        if (C_ASSOCIATED(ptr, ptr_arr)) then
           win  = windows(i)
           exit
        endif
     endif
  enddo
  end subroutine get_window_1D_int8
  !---------------------------
  subroutine get_window_2D_int8(arr, win)
  implicit none
  integer(int8), pointer, intent(in)  :: arr(:,:)
  integer,                intent(out) :: win
  integer(MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr_arr, ptr
  integer     :: disp_unit, ierr, i

  ptr_arr   = C_LOC(arr)
  win       = MPI_WIN_NULL
  disp_unit = int8_size
  do i=1, num_windows
     if (windows(i) /= MPI_WIN_NULL) then
        call MPI_WIN_SHARED_QUERY(windows(i), 0, wsize, disp_unit, ptr, ierr)
        if (C_ASSOCIATED(ptr, ptr_arr)) then
           win  = windows(i)
           exit
        endif
     endif
  enddo
  end subroutine get_window_2D_int8
  !---------------------------
  subroutine get_window_3D_int8(arr, win)
  implicit none
  integer(int8), pointer, intent(in)  :: arr(:,:,:)
  integer,                intent(out) :: win
  integer(MPI_ADDRESS_KIND) :: wsize
  type(C_PTR) :: ptr_arr, ptr
  integer     :: disp_unit, ierr, i

  ptr_arr   = C_LOC(arr)
  win       = MPI_WIN_NULL
  disp_unit = int8_size
  do i=1, num_windows
     if (windows(i) /= MPI_WIN_NULL) then
        call MPI_WIN_SHARED_QUERY(windows(i), 0, wsize, disp_unit, ptr, ierr)
        if (C_ASSOCIATED(ptr, ptr_arr)) then
           win  = windows(i)
           exit
        endif
     endif
  enddo
  end subroutine get_window_3D_int8
  !---------------------------
  subroutine destroy_mem_1D_int8(arr)
  implicit none
  integer(int8), pointer, intent(inout) :: arr(:)
  integer :: win, ierr
  call get_window_1D_int8(arr, win)
  if (win /= MPI_WIN_NULL) then
     call destroy_shared_mem_window(win)
     nullify(arr)
  else
     if (associated(arr)) deallocate(arr)
  endif
  end subroutine destroy_mem_1D_int8
  !---------------------------
  subroutine destroy_mem_2D_int8(arr)
  implicit none
  integer(int8), pointer, intent(inout) :: arr(:,:)
  integer :: win, ierr
  call get_window_2D_int8(arr, win)
  if (win /= MPI_WIN_NULL) then
     call destroy_shared_mem_window(win)
     nullify(arr)
  else
     if (associated(arr)) deallocate(arr)
  endif
  end subroutine destroy_mem_2D_int8
  !---------------------------
  subroutine destroy_mem_3D_int8(arr)
  implicit none
  integer(int8), pointer, intent(inout) :: arr(:,:,:)
  integer :: win, ierr
  call get_window_3D_int8(arr, win)
  if (win /= MPI_WIN_NULL) then
     call destroy_shared_mem_window(win)
     nullify(arr)
  else
     if (associated(arr)) deallocate(arr)
  endif
  end subroutine destroy_mem_3D_int8
  !---------------------------
  subroutine destroy_shared_mem_all()
  implicit none
  integer, parameter :: assert = 0
  integer :: ierr, i
  do i=1, num_windows
     if (windows(i) /= MPI_WIN_NULL) then
        call MPI_WIN_FENCE(assert,windows(i),ierr)
        call MPI_WIN_FREE(windows(i),ierr)
     endif
  enddo
  !-- there is no way to nullify fortran array pointer within this module.
  !-- num_windows should be reset.
  num_windows = 0
  end subroutine destroy_shared_mem_all
  !---------------------------
  subroutine destroy_shared_mem_window(win)
  implicit none
  integer, intent(in) :: win
  integer, parameter  :: assert = 0
  integer :: ierr, i
  do i=1, num_windows
     if (win == windows(i) .and.  windows(i) /= MPI_WIN_NULL) then
        call MPI_WIN_FENCE(assert,windows(i),ierr)
        call MPI_WIN_FREE(windows(i),ierr)
        exit
     endif
  enddo
  end subroutine destroy_shared_mem_window
  !---------------------------
end module memory_mod
