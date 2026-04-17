module read_atmosphere_mod
  use define
  use mathlib
  use utility
contains
  !=================================================================
  !--- for par%geometry == 'plane_atmosphere'
  subroutine read_plane_atmosphere_data(fname,array,grid)
  implicit none
  character(len=*), intent(in)    :: fname
  real(kind=wp),    intent(inout) :: array(:,:,:)
  type(grid_type),  intent(in)    :: grid

  !--- local variables
  real(kind=wp), allocatable :: xarr(:), yarr(:)
  real(kind=wp) :: x1,y1,xnew,ynew,x2
  integer :: unit, info
  integer :: n, i, loop1, loop2

  open(newunit=unit,file=trim(fname),status='old')
  if (.not.allocated(xarr)) allocate(xarr(0))
  if (.not.allocated(yarr)) allocate(yarr(0))
  n = 0
  do while(.true.)
     read(unit,*,iostat=info) x1,y1
     if (info /= 0) exit
     n    = n + 1
     xarr = [xarr, x1]
     yarr = [yarr, y1]
  enddo
  close(unit)

  x1 = minval(xarr)
  x2 = maxval(xarr)
  call loop_divide(grid%nz,mpar%h_nproc,mpar%h_rank,loop1,loop2)
  do i=loop1, loop2
     xnew = (grid%zface(i)+grid%zface(i+1))/2.0
     if (xnew >= x1 .and. xnew <= x2) then
        call interp(xarr,yarr,xnew,ynew)
     endif
     array(:,:,i) = ynew
  enddo
  if (allocated(xarr)) deallocate(xarr)
  if (allocated(yarr)) deallocate(yarr)
  end subroutine read_plane_atmosphere_data
  !=================================================================
  !--- for par%geometry == 'spherical_atmosphere'
  subroutine read_spherical_atmosphere_data(fname,array,grid)
  implicit none
  character(len=*), intent(in)    :: fname
  real(kind=wp),    intent(inout) :: array(:,:,:)
  type(grid_type),  intent(in)    :: grid

  !--- local variables
  real(kind=wp), allocatable :: xarr(:), yarr(:)
  real(kind=wp) :: x1,y1,x,y,z,r,ynew,x2
  integer :: unit, info
  integer :: n, i, j, k, loop, loop1, loop2, nsize

  open(newunit=unit,file=trim(fname),status='old')
  if (.not.allocated(xarr)) allocate(xarr(0))
  if (.not.allocated(yarr)) allocate(yarr(0))
  n = 0
  do while(.true.)
     read(unit,*,iostat=info) x1,y1
     if (info /= 0) exit
     n    = n + 1
     xarr = [xarr, x1]
     yarr = [yarr, y1]
  enddo
  close(unit)

  x1 = minval(xarr)
  x2 = maxval(xarr)
  nsize = grid%nx * grid%ny * grid%nz
  call loop_divide(nsize,mpar%h_nproc,mpar%h_rank,loop1,loop2)
  do loop=loop1, loop2
     call array_3D_indices(grid%nx,grid%ny,grid%nz,loop,i,j,k)
     x = (grid%xface(i)+grid%xface(i+1))/2.0
     y = (grid%yface(j)+grid%yface(j+1))/2.0
     z = (grid%zface(k)+grid%zface(k+1))/2.0
     r = sqrt(x**2 + y**2 + z**2)
     if (r >= x1 .and. r <= x2) then
        call interp(xarr,yarr,r,ynew)
        array(i,j,k) = ynew
     endif
  enddo
  if (allocated(xarr)) deallocate(xarr)
  if (allocated(yarr)) deallocate(yarr)
  end subroutine read_spherical_atmosphere_data
  !=================================================================
  !--- for par%geometry == 'spherical_atmosphere'
  subroutine read_spherical_atmosphere_velocity(fname,vx,vy,vz,grid)
  implicit none
  character(len=*), intent(in)    :: fname
  real(kind=wp),    intent(inout) :: vx(:,:,:), vy(:,:,:), vz(:,:,:)
  type(grid_type),  intent(in)    :: grid

  !--- local variables
  real(kind=wp), allocatable :: xarr(:), yarr(:)
  real(kind=wp) :: x1,y1,x,y,z,r,ynew,x2
  integer :: unit, info
  integer :: n, i, j, k, loop, loop1, loop2, nsize

  open(newunit=unit,file=trim(fname),status='old')
  if (.not.allocated(xarr)) allocate(xarr(0))
  if (.not.allocated(yarr)) allocate(yarr(0))
  n = 0
  do while(.true.)
     read(unit,*,iostat=info) x1,y1
     if (info /= 0) exit
     n    = n + 1
     xarr = [xarr, x1]
     yarr = [yarr, y1]
  enddo
  close(unit)

  x1 = minval(xarr)
  x2 = maxval(xarr)
  nsize = grid%nx * grid%ny * grid%nz
  call loop_divide(nsize,mpar%h_nproc,mpar%h_rank,loop1,loop2)
  do loop=loop1, loop2
     call array_3D_indices(grid%nx,grid%ny,grid%nz,loop,i,j,k)
     x = (grid%xface(i)+grid%xface(i+1))/2.0
     y = (grid%yface(j)+grid%yface(j+1))/2.0
     z = (grid%zface(k)+grid%zface(k+1))/2.0
     r = sqrt(x**2 + y**2 + z**2)
     call interp(xarr,yarr,r,ynew)
     if (r > 0.0_wp .and. r >= x1 .and. r <= x2) then
        vx(i,j,k) = x/r * ynew
        vy(i,j,k) = y/r * ynew
        vz(i,j,k) = z/r * ynew
     endif
  enddo
  if (allocated(xarr)) deallocate(xarr)
  if (allocated(yarr)) deallocate(yarr)
  end subroutine read_spherical_atmosphere_velocity
  !=================================================================
end module read_atmosphere_mod
