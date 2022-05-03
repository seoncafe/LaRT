module memory_mod
  use, intrinsic :: iso_fortran_env, only: real32, real64, real128, int8, int16, int32, int64
#ifdef MPI
  use mpi
#endif
  !--
  !-- 2021-06-29. Now, the routines work for both "real32" and "real64" types.
  !--             The routines for "int8", "int16", "int32", and "inte64" type array are added.
  !--
  implicit none
  private
  logical :: mpi_parameters_defined = .false.
  integer :: p_rank, num_nodes

  public create_mem, reduce_mem
  interface create_mem
     module procedure create_mem_1D_real64, create_mem_2D_real64, create_mem_3D_real64, create_mem_4D_real64, create_mem_5D_real64, &
                      create_mem_1D_real32, create_mem_2D_real32, create_mem_3D_real32, create_mem_4D_real32, create_mem_5D_real32, &
                      create_mem_1D_int64, create_mem_2D_int64, create_mem_3D_int64, &
                      create_mem_1D_int32, create_mem_2D_int32, create_mem_3D_int32, &
                      create_mem_1D_int16, create_mem_2D_int16, create_mem_3D_int16, &
                      create_mem_1D_int8,  create_mem_2D_int8,  create_mem_3D_int8
  end interface create_mem
  interface reduce_mem
     module procedure reduce_mem_1D_real64, reduce_mem_2D_real64, reduce_mem_3D_real64, reduce_mem_4D_real64, reduce_mem_5D_real64, &
                      reduce_mem_1D_real32, reduce_mem_2D_real32, reduce_mem_3D_real32, reduce_mem_4D_real32, reduce_mem_5D_real32, &
                      reduce_mem_1D_int64, reduce_mem_2D_int64, reduce_mem_3D_int64, &
                      reduce_mem_1D_int32, reduce_mem_2D_int32, reduce_mem_3D_int32, &
                      reduce_mem_1D_int16, reduce_mem_2D_int16, reduce_mem_3D_int16, &
                      reduce_mem_1D_int8,  reduce_mem_2D_int8,  reduce_mem_3D_int8
  end interface reduce_mem

contains
  !---------------------------
  subroutine init_mpi_parameters()
  integer :: ierr
#ifdef MPI
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_nodes, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, p_rank,    ierr)
  mpi_parameters_defined = .true.
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#else
  num_nodes = 1
  p_rank    = 0
  mpi_parameters_defined = .true.
#endif
  end subroutine init_mpi_parameters
  !---------------------------
  subroutine create_mem_1D_real64(array,arrshape)
  implicit none
  real(kind=real64), pointer, intent(inout) :: array(:)
  integer,                    intent(in)    :: arrshape(1)
 
  if (.not. associated(array)) allocate(array(arrshape(1)))
  array = 0.0_real64
  end subroutine create_mem_1D_real64
  !---------------------------
  subroutine create_mem_2D_real64(array,arrshape)
  implicit none
  real(kind=real64), pointer, intent(inout) :: array(:,:)
  integer,                    intent(in)    :: arrshape(2)

  if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2)))
  array = 0.0_real64
  end subroutine create_mem_2D_real64
  !---------------------------
  subroutine create_mem_3D_real64(array,arrshape)
  implicit none
  real(kind=real64), pointer, intent(inout) :: array(:,:,:)
  integer,                    intent(in)    :: arrshape(3)

  if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2),arrshape(3)))
  array = 0.0_real64
  end subroutine create_mem_3D_real64
  !---------------------------
  subroutine create_mem_4D_real64(array,arrshape)
  implicit none
  real(kind=real64), pointer, intent(inout) :: array(:,:,:,:)
  integer,                    intent(in)    :: arrshape(4)

  if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2),arrshape(3),arrshape(4)))
  array = 0.0_real64
  end subroutine create_mem_4D_real64
  !---------------------------
  subroutine create_mem_5D_real64(array,arrshape)
  implicit none
  real(kind=real64), pointer, intent(inout) :: array(:,:,:,:,:)
  integer,                    intent(in)    :: arrshape(5)

  if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2),arrshape(3),arrshape(4),arrshape(5)))
  array = 0.0_real64
  end subroutine create_mem_5D_real64
  !---------------------------
  subroutine reduce_mem_1D_real64(array)
  implicit none
  real(kind=real64), pointer, intent(inout) :: array(:)
  integer :: ierr, nsize

#ifdef MPI
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (num_nodes > 1) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:),nsize,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     else
        call MPI_REDUCE(array(:),           0,nsize,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     endif
  endif
#endif
  end subroutine reduce_mem_1D_real64
  !---------------------------
  subroutine reduce_mem_2D_real64(array)
  implicit none
  real(kind=real64), pointer, intent(inout) :: array(:,:)
  integer :: ierr, nsize

#ifdef MPI
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (num_nodes > 1) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:,:),nsize,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     else
        call MPI_REDUCE(array(:,:),           0,nsize,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     endif
  endif
#endif
  end subroutine reduce_mem_2D_real64
  !---------------------------
  subroutine reduce_mem_3D_real64(array)
  implicit none
  real(kind=real64), pointer, intent(inout) :: array(:,:,:)
  integer :: ierr, k, nsize

#ifdef MPI
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (num_nodes > 1) then
     nsize = size(array(:,:,1))
     do k=1, size(array(1,1,:))
        if (p_rank == 0) then
           call MPI_REDUCE(MPI_IN_PLACE,array(:,:,k),nsize,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        else
           call MPI_REDUCE(array(:,:,k),           0,nsize,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        endif
     enddo
  endif
#endif
  end subroutine reduce_mem_3D_real64
  !---------------------------
  subroutine reduce_mem_4D_real64(array)
  implicit none
  real(kind=real64), pointer, intent(inout) :: array(:,:,:,:)
  integer :: ierr, j, k, nsize

#ifdef MPI
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (num_nodes > 1) then
     nsize = size(array(:,:,1,1))
     do k=1, size(array(1,1,1,:))
     do j=1, size(array(1,1,:,1))
        if (p_rank == 0) then
           call MPI_REDUCE(MPI_IN_PLACE,array(:,:,j,k),nsize,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        else
           call MPI_REDUCE(array(:,:,j,k),           0,nsize,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        endif
     enddo
     enddo
  endif
#endif
  end subroutine reduce_mem_4D_real64
  !---------------------------
  subroutine reduce_mem_5D_real64(array)
  implicit none
  real(kind=real64), pointer, intent(inout) :: array(:,:,:,:,:)
  integer :: ierr, i, j, k, nsize

#ifdef MPI
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (num_nodes > 1) then
     nsize = size(array(:,:,1,1,1))
     do k=1, size(array(1,1,1,1,:))
     do j=1, size(array(1,1,1,:,1))
     do i=1, size(array(1,1,:,1,1))
        if (p_rank == 0) then
           call MPI_REDUCE(MPI_IN_PLACE,array(:,:,i,j,k),nsize,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        else
           call MPI_REDUCE(array(:,:,i,j,k),           0,nsize,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        endif
     enddo
     enddo
     enddo
  endif
#endif
  end subroutine reduce_mem_5D_real64
  !---------------------------
  subroutine create_mem_1D_real32(array,arrshape)
  implicit none
  real(kind=real32), pointer, intent(inout) :: array(:)
  integer,                    intent(in)    :: arrshape(1)
 
  if (.not. associated(array)) allocate(array(arrshape(1)))
  array = 0.0_real32
  end subroutine create_mem_1D_real32
  !---------------------------
  subroutine create_mem_2D_real32(array,arrshape)
  implicit none
  real(kind=real32), pointer, intent(inout) :: array(:,:)
  integer,                    intent(in)    :: arrshape(2)

  if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2)))
  array = 0.0_real32
  end subroutine create_mem_2D_real32
  !---------------------------
  subroutine create_mem_3D_real32(array,arrshape)
  implicit none
  real(kind=real32), pointer, intent(inout) :: array(:,:,:)
  integer,                    intent(in)    :: arrshape(3)

  if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2),arrshape(3)))
  array = 0.0_real32
  end subroutine create_mem_3D_real32
  !---------------------------
  subroutine create_mem_4D_real32(array,arrshape)
  implicit none
  real(kind=real32), pointer, intent(inout) :: array(:,:,:,:)
  integer,                    intent(in)    :: arrshape(4)

  if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2),arrshape(3),arrshape(4)))
  array = 0.0_real32
  end subroutine create_mem_4D_real32
  !---------------------------
  subroutine create_mem_5D_real32(array,arrshape)
  implicit none
  real(kind=real32), pointer, intent(inout) :: array(:,:,:,:,:)
  integer,                    intent(in)    :: arrshape(5)

  if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2),arrshape(3),arrshape(4),arrshape(5)))
  array = 0.0_real32
  end subroutine create_mem_5D_real32
  !---------------------------
  subroutine reduce_mem_1D_real32(array)
  implicit none
  real(kind=real32), pointer, intent(inout) :: array(:)
  integer :: ierr, nsize

#ifdef MPI
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (num_nodes > 1) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:),nsize,MPI_REAL4,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     else
        call MPI_REDUCE(array(:),           0,nsize,MPI_REAL4,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     endif
  endif
#endif
  end subroutine reduce_mem_1D_real32
  !---------------------------
  subroutine reduce_mem_2D_real32(array)
  implicit none
  real(kind=real32), pointer, intent(inout) :: array(:,:)
  integer :: ierr, nsize

#ifdef MPI
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (num_nodes > 1) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:,:),nsize,MPI_REAL4,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     else
        call MPI_REDUCE(array(:,:),           0,nsize,MPI_REAL4,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     endif
  endif
#endif
  end subroutine reduce_mem_2D_real32
  !---------------------------
  subroutine reduce_mem_3D_real32(array)
  implicit none
  real(kind=real32), pointer, intent(inout) :: array(:,:,:)
  integer :: ierr, k, nsize

#ifdef MPI
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (num_nodes > 1) then
     nsize = size(array(:,:,1))
     do k=1, size(array(1,1,:))
        if (p_rank == 0) then
           call MPI_REDUCE(MPI_IN_PLACE,array(:,:,k),nsize,MPI_REAL4,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        else
           call MPI_REDUCE(array(:,:,k),           0,nsize,MPI_REAL4,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        endif
     enddo
  endif
#endif
  end subroutine reduce_mem_3D_real32
  !---------------------------
  subroutine reduce_mem_4D_real32(array)
  implicit none
  real(kind=real32), pointer, intent(inout) :: array(:,:,:,:)
  integer :: ierr, j, k, nsize

#ifdef MPI
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (num_nodes > 1) then
     nsize = size(array(:,:,1,1))
     do k=1, size(array(1,1,1,:))
     do j=1, size(array(1,1,:,1))
        if (p_rank == 0) then
           call MPI_REDUCE(MPI_IN_PLACE,array(:,:,j,k),nsize,MPI_REAL4,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        else
           call MPI_REDUCE(array(:,:,j,k),           0,nsize,MPI_REAL4,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        endif
     enddo
     enddo
  endif
#endif
  end subroutine reduce_mem_4D_real32
  !---------------------------
  subroutine reduce_mem_5D_real32(array)
  implicit none
  real(kind=real32), pointer, intent(inout) :: array(:,:,:,:,:)
  integer :: ierr, i, j, k, nsize

#ifdef MPI
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (num_nodes > 1) then
     nsize = size(array(:,:,1,1,1))
     do k=1, size(array(1,1,1,1,:))
     do j=1, size(array(1,1,1,:,1))
     do i=1, size(array(1,1,:,1,1))
        if (p_rank == 0) then
           call MPI_REDUCE(MPI_IN_PLACE,array(:,:,i,j,k),nsize,MPI_REAL4,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        else
           call MPI_REDUCE(array(:,:,i,j,k),           0,nsize,MPI_REAL4,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        endif
     enddo
     enddo
     enddo
  endif
#endif
  end subroutine reduce_mem_5D_real32
  !---------------------------
  subroutine create_mem_1D_int64(array,arrshape)
  implicit none
  integer(int64), pointer, intent(inout) :: array(:)
  integer,                 intent(in)    :: arrshape(1)

  if (.not. associated(array)) allocate(array(arrshape(1)))
  array = 0_int64
  end subroutine create_mem_1D_int64
  !---------------------------
  subroutine create_mem_2D_int64(array,arrshape)
  implicit none
  integer(int64), pointer, intent(inout) :: array(:,:)
  integer,                 intent(in)    :: arrshape(2)

  if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2)))
  array = 0_int64
  end subroutine create_mem_2D_int64
  !---------------------------
  subroutine create_mem_3D_int64(array,arrshape)
  implicit none
  integer(int64), pointer, intent(inout) :: array(:,:,:)
  integer,                 intent(in)    :: arrshape(3)

  if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2),arrshape(3)))
  array = 0_int64
  end subroutine create_mem_3D_int64
  !---------------------------
  subroutine reduce_mem_1D_int64(array)
  implicit none
  integer(int64), pointer, intent(inout) :: array(:)
  integer :: ierr, nsize

#ifdef MPI
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (num_nodes > 1) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:),nsize,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     else
        call MPI_REDUCE(array(:),           0,nsize,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     endif
  endif
#endif
  end subroutine reduce_mem_1D_int64
  !---------------------------
  subroutine reduce_mem_2D_int64(array)
  implicit none
  integer(int64), pointer, intent(inout) :: array(:,:)
  integer :: ierr, nsize

#ifdef MPI
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (num_nodes > 1) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:,:),nsize,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     else
        call MPI_REDUCE(array(:,:),           0,nsize,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     endif
  endif
#endif
  end subroutine reduce_mem_2D_int64
  !---------------------------
  subroutine reduce_mem_3D_int64(array)
  implicit none
  integer(int64), pointer, intent(inout) :: array(:,:,:)
  integer :: ierr, k, nsize

#ifdef MPI
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (num_nodes > 1) then
     nsize = size(array(:,:,1))
     do k=1, size(array(1,1,:))
        if (p_rank == 0) then
           call MPI_REDUCE(MPI_IN_PLACE,array(:,:,k),nsize,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        else
           call MPI_REDUCE(array(:,:,k),           0,nsize,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        endif
     enddo
  endif
#endif
  end subroutine reduce_mem_3D_int64
  !---------------------------
  subroutine create_mem_1D_int32(array,arrshape)
  implicit none
  integer(int32), pointer, intent(inout) :: array(:)
  integer,                 intent(in)    :: arrshape(1)

  if (.not. associated(array)) allocate(array(arrshape(1)))
  array = 0_int32
  end subroutine create_mem_1D_int32
  !---------------------------
  subroutine create_mem_2D_int32(array,arrshape)
  implicit none
  integer(int32), pointer, intent(inout) :: array(:,:)
  integer,                 intent(in)    :: arrshape(2)

  if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2)))
  array = 0_int32
  end subroutine create_mem_2D_int32
  !---------------------------
  subroutine create_mem_3D_int32(array,arrshape)
  implicit none
  integer(int32), pointer, intent(inout) :: array(:,:,:)
  integer,                 intent(in)    :: arrshape(3)

  if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2),arrshape(3)))
  array = 0_int32
  end subroutine create_mem_3D_int32
  !---------------------------
  subroutine reduce_mem_1D_int32(array)
  implicit none
  integer(int32), pointer, intent(inout) :: array(:)
  integer :: ierr, nsize

#ifdef MPI
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (num_nodes > 1) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:),nsize,MPI_INTEGER4,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     else
        call MPI_REDUCE(array(:),           0,nsize,MPI_INTEGER4,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     endif
  endif
#endif
  end subroutine reduce_mem_1D_int32
  !---------------------------
  subroutine reduce_mem_2D_int32(array)
  implicit none
  integer(int32), pointer, intent(inout) :: array(:,:)
  integer :: ierr, nsize

#ifdef MPI
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (num_nodes > 1) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:,:),nsize,MPI_INTEGER4,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     else
        call MPI_REDUCE(array(:,:),           0,nsize,MPI_INTEGER4,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     endif
  endif
#endif
  end subroutine reduce_mem_2D_int32
  !---------------------------
  subroutine reduce_mem_3D_int32(array)
  implicit none
  integer(int32), pointer, intent(inout) :: array(:,:,:)
  integer :: ierr, k, nsize

#ifdef MPI
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (num_nodes > 1) then
     nsize = size(array(:,:,1))
     do k=1, size(array(1,1,:))
        if (p_rank == 0) then
           call MPI_REDUCE(MPI_IN_PLACE,array(:,:,k),nsize,MPI_INTEGER4,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        else
           call MPI_REDUCE(array(:,:,k),           0,nsize,MPI_INTEGER4,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        endif
     enddo
  endif
#endif
  end subroutine reduce_mem_3D_int32
  !---------------------------
  subroutine create_mem_1D_int16(array,arrshape)
  implicit none
  integer(int16), pointer, intent(inout) :: array(:)
  integer,                 intent(in)    :: arrshape(1)

  if (.not. associated(array)) allocate(array(arrshape(1)))
  array = 0_int16
  end subroutine create_mem_1D_int16
  !---------------------------
  subroutine create_mem_2D_int16(array,arrshape)
  implicit none
  integer(int16), pointer, intent(inout) :: array(:,:)
  integer,                 intent(in)    :: arrshape(2)

  if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2)))
  array = 0_int16
  end subroutine create_mem_2D_int16
  !---------------------------
  subroutine create_mem_3D_int16(array,arrshape)
  implicit none
  integer(int16), pointer, intent(inout) :: array(:,:,:)
  integer,                 intent(in)    :: arrshape(3)

  if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2),arrshape(3)))
  array = 0_int16
  end subroutine create_mem_3D_int16
  !---------------------------
  subroutine reduce_mem_1D_int16(array)
  implicit none
  integer(int16), pointer, intent(inout) :: array(:)
  integer :: ierr, nsize

#ifdef MPI
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (num_nodes > 1) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:),nsize,MPI_INTEGER2,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     else
        call MPI_REDUCE(array(:),           0,nsize,MPI_INTEGER2,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     endif
  endif
#endif
  end subroutine reduce_mem_1D_int16
  !---------------------------
  subroutine reduce_mem_2D_int16(array)
  implicit none
  integer(int16), pointer, intent(inout) :: array(:,:)
  integer :: ierr, nsize

#ifdef MPI
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (num_nodes > 1) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:,:),nsize,MPI_INTEGER2,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     else
        call MPI_REDUCE(array(:,:),           0,nsize,MPI_INTEGER2,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     endif
  endif
#endif
  end subroutine reduce_mem_2D_int16
  !---------------------------
  subroutine reduce_mem_3D_int16(array)
  implicit none
  integer(int16), pointer, intent(inout) :: array(:,:,:)
  integer :: ierr, k, nsize

#ifdef MPI
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (num_nodes > 1) then
     nsize = size(array(:,:,1))
     do k=1, size(array(1,1,:))
        if (p_rank == 0) then
           call MPI_REDUCE(MPI_IN_PLACE,array(:,:,k),nsize,MPI_INTEGER2,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        else
           call MPI_REDUCE(array(:,:,k),           0,nsize,MPI_INTEGER2,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        endif
     enddo
  endif
#endif
  end subroutine reduce_mem_3D_int16
  !---------------------------
  subroutine create_mem_1D_int8(array,arrshape)
  implicit none
  integer(int8), pointer, intent(inout) :: array(:)
  integer,                intent(in)    :: arrshape(1)

  if (.not. associated(array)) allocate(array(arrshape(1)))
  array = 0_int8
  end subroutine create_mem_1D_int8
  !---------------------------
  subroutine create_mem_2D_int8(array,arrshape)
  implicit none
  integer(int8), pointer, intent(inout) :: array(:,:)
  integer,                intent(in)    :: arrshape(2)

  if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2)))
  array = 0_int8
  end subroutine create_mem_2D_int8
  !---------------------------
  subroutine create_mem_3D_int8(array,arrshape)
  implicit none
  integer(int8), pointer, intent(inout) :: array(:,:,:)
  integer,                intent(in)    :: arrshape(3)

  if (.not. associated(array)) allocate(array(arrshape(1),arrshape(2),arrshape(3)))
  array = 0_int8
  end subroutine create_mem_3D_int8
  !---------------------------
  subroutine reduce_mem_1D_int8(array)
  implicit none
  integer(int8), pointer, intent(inout) :: array(:)
  integer :: ierr, nsize

#ifdef MPI
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (num_nodes > 1) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:),nsize,MPI_INTEGER1,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     else
        call MPI_REDUCE(array(:),           0,nsize,MPI_INTEGER1,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     endif
  endif
#endif
  end subroutine reduce_mem_1D_int8
  !---------------------------
  subroutine reduce_mem_2D_int8(array)
  implicit none
  integer(int8), pointer, intent(inout) :: array(:,:)
  integer :: ierr, nsize

#ifdef MPI
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (num_nodes > 1) then
     nsize = size(array)
     if (p_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,array(:,:),nsize,MPI_INTEGER1,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     else
        call MPI_REDUCE(array(:,:),           0,nsize,MPI_INTEGER1,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     endif
  endif
#endif
  end subroutine reduce_mem_2D_int8
  !---------------------------
  subroutine reduce_mem_3D_int8(array)
  implicit none
  integer(int8), pointer, intent(inout) :: array(:,:,:)
  integer :: ierr, k, nsize

#ifdef MPI
  if (.not.mpi_parameters_defined) call init_mpi_parameters()
  if (num_nodes > 1) then
     nsize = size(array(:,:,1))
     do k=1, size(array(1,1,:))
        if (p_rank == 0) then
           call MPI_REDUCE(MPI_IN_PLACE,array(:,:,k),nsize,MPI_INTEGER1,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        else
           call MPI_REDUCE(array(:,:,k),           0,nsize,MPI_INTEGER1,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        endif
     enddo
  endif
#endif
  end subroutine reduce_mem_3D_int8
  !---------------------------
end module memory_mod
