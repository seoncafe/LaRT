module read_grid_data
  !
  ! Reads gridded input data (density/temperature/velocity fields, plus
  ! generic 1-D / 3-D / 4-D arrays) through the iofile_mod facade so that
  ! both FITS (.fits / .fits.gz) and HDF5 (.h5 / .hdf5) inputs work
  ! transparently.  Format is detected from the filename extension by
  ! io_open_old.
  !
  ! HDF5 input convention (matches what LaRT writes itself):
  !   array data is at /<section>/data  (e.g. /Spectrum/data, /section_001/data),
  !   or at /<section>/<colname>        (table-style; not used by these readers).
  ! Shape queries that target FITS-style 'NAXIS', 'NAXIS1', ... keywords are
  ! transparently answered from the dataset shape on the HDF5 backend.
  !
  ! Renamed from `read_fits_data` (2026-05) once the module learned HDF5.
  !
  use define
  use iofile_mod
contains
  !=================================================================
  subroutine read_velocity(fname,vfx,vfy,vfz,reduce_factor)
  implicit none
  character(len=*),    intent(in) :: fname
  real(kind=wp),    intent(inout) :: vfx(:,:,:), vfy(:,:,:), vfz(:,:,:)
  integer, optional,   intent(in) :: reduce_factor
  real(kind=wp), allocatable :: velo(:,:,:,:)
  integer :: nx,ny,nz

  nx = size(vfx,1)
  ny = size(vfx,2)
  nz = size(vfx,3)

  if (.not. allocated(velo)) allocate(velo(3,nx,ny,nz))
  if (present(reduce_factor)) then
     call read_4D(fname,velo,reduce_factor=reduce_factor)
  else
     call read_4D(fname,velo)
  endif
  vfx(:,:,:) = velo(1,:,:,:)
  vfy(:,:,:) = velo(2,:,:,:)
  vfz(:,:,:) = velo(3,:,:,:)
  if (allocated(velo)) deallocate(velo)
  end subroutine read_velocity

  !=================================================================
  subroutine read_3D(fname,array,reduce_factor,centering)
  implicit none
  character(len=*),    intent(in) :: fname
  real(kind=wp),    intent(inout) :: array(:,:,:)
  integer, optional,   intent(in) :: reduce_factor, centering

  !--- local variables
  type(io_file_type) :: iofh
  integer :: status = 0
  integer :: n1,n2,n3
  integer :: loc(3),n1cen,n2cen,n3cen
  integer :: i,j,k,i1,i2,j1,j2,k1,k2
  !--- note the following should be always a single precision variable. (not any more, 2020.09.02).
  real(kind=wp), allocatable :: arr(:,:,:)

  call io_open_old(iofh,trim(fname),status)
  if (status == 0) then
     ! HDF5: io_open_old auto-positions at the first image-style section.
     ! FITS: cur HDU is 1 (Primary); the NAXIS keywords live there if the
     !       data was written as a single image.
     call io_get_keyword(iofh,'NAXIS1',n1,status)
     call io_get_keyword(iofh,'NAXIS2',n2,status)
     call io_get_keyword(iofh,'NAXIS3',n3,status)
     if (.not. allocated(arr)) allocate(arr(n1,n2,n3))
     call io_read_image(iofh,arr,status)
     call io_close(iofh,status)
  else
     write(*,*) 'Error in reading the data file :', trim(fname),' status = ',status
  endif

  if (present(reduce_factor)) then
     if (reduce_factor > 1) then
        n1 = n1/reduce_factor
        n2 = n2/reduce_factor
        n3 = n3/reduce_factor
        do k=1,n3
        do j=1,n2
        do i=1,n1
           i2 = reduce_factor*i
           i1 = i2 - (reduce_factor-1)
           j2 = reduce_factor*j
           j1 = j2 - (reduce_factor-1)
           k2 = reduce_factor*k
           k1 = k2 - (reduce_factor-1)
           array(i,j,k) = sum(arr(i1:i2,j1:j2,k1:k2))/reduce_factor**3
        enddo
        enddo
        enddo
     else
        array(:,:,:) = arr(:,:,:)
     endif
  else
     array(:,:,:) = arr(:,:,:)
  endif
  if (allocated(arr)) deallocate(arr)

  if (present(centering)) then
     if (centering == 1 .or. centering == 2) then
        if (.not.allocated(arr)) allocate(arr(n1,n2,n3))
        n1cen = n1/2
        n2cen = n2/2
        n3cen = n3/2
        if (centering == 1) then
           loc = minloc(array)
        else
           loc = maxloc(array)
        endif
        if (loc(1) < n1cen) then
           arr(n1cen-loc(1)+1:n1,:,:)   = array(1:loc(1)+n1cen,:,:)
           arr(1:n1cen-loc(1),:,:)      = array(loc(1)+n1cen+1:n1,:,:)
           array(:,:,:) = arr(:,:,:)
        else if (loc(1) > n1cen) then
           arr(1:3*n1cen-loc(1),:,:)    = array(loc(1)-n1cen+1:n1,:,:)
           arr(3*n1cen-loc(1)+1:n1,:,:) = array(1:loc(1)-n1cen,:,:)
           array(:,:,:) = arr(:,:,:)
        endif
        if (loc(2) < n2cen) then
           arr(:,n2cen-loc(2)+1:n2,:)   = array(:,1:loc(2)+n2cen,:)
           arr(:,1:n2cen-loc(2),:)      = array(:,loc(2)+n2cen+1:n2,:)
           array(:,:,:) = arr(:,:,:)
        else if (loc(2) > n2cen) then
           arr(:,1:3*n2cen-loc(2),:)    = array(:,loc(2)-n2cen+1:n2,:)
           arr(:,3*n2cen-loc(2)+1:n2,:) = array(:,1:loc(2)-n2cen,:)
           array(:,:,:) = arr(:,:,:)
        endif
        if (loc(3) < n3cen) then
           arr(:,:,n3cen-loc(3)+1:n3)   = array(:,:,1:loc(3)+n3cen)
           arr(:,:,1:n3cen-loc(3))      = array(:,:,loc(3)+n2cen+1:n3)
           array(:,:,:) = arr(:,:,:)
        else if (loc(3) > n3cen) then
           arr(:,:,1:3*n3cen-loc(3))    = array(:,:,loc(3)-n3cen+1:n3)
           arr(:,:,3*n3cen-loc(3)+1:n3) = array(:,:,1:loc(3)-n3cen)
           array(:,:,:) = arr(:,:,:)
        endif
        if (allocated(arr)) deallocate(arr)
        if (centering == 1) then
           loc = minloc(array)
        else
           loc = maxloc(array)
        endif
        if (loc(1) /= n1cen .or. loc(2) /= n2cen .or. loc(3) /= n3cen) then
           write(6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
           write(6,*) 'Something wrong. Please check the centering algorithm.'
           write(6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        endif
     endif
  endif

  end subroutine read_3D
!--------------------------------------------------------------
  subroutine read_4D(fname,array,reduce_factor)
  implicit none
  character(len=*), intent(in)    :: fname
  real(kind=wp),    intent(inout) :: array(:,:,:,:)
  integer, optional,   intent(in) :: reduce_factor

  !--- local variables
  type(io_file_type) :: iofh
  integer :: status = 0
  integer :: n1,n2,n3,n4
  integer :: i,j,k,i1,i2,j1,j2,k1,k2
  !--- note the following should be always a single precision variable. (2020.09.02, not any more).
  real(kind=wp), allocatable :: arr(:,:,:,:)

  call io_open_old(iofh,trim(fname),status)
  if (status == 0) then
     call io_get_keyword(iofh,'NAXIS1',n1,status)
     call io_get_keyword(iofh,'NAXIS2',n2,status)
     call io_get_keyword(iofh,'NAXIS3',n3,status)
     call io_get_keyword(iofh,'NAXIS4',n4,status)
     if (.not. allocated(arr)) allocate(arr(n1,n2,n3,n4))
     call io_read_image(iofh,arr,status)
  else
     write(*,*) 'Error in reading the data file :', trim(fname),' status = ',status
  endif
  call io_close(iofh,status)

  if (present(reduce_factor)) then
     if (reduce_factor > 1) then
        n2 = n2/reduce_factor
        n3 = n3/reduce_factor
        n4 = n4/reduce_factor
        do k=1,n4
        do j=1,n3
        do i=1,n2
           i2 = reduce_factor*i
           i1 = i2 - (reduce_factor-1)
           j2 = reduce_factor*j
           j1 = j2 - (reduce_factor-1)
           k2 = reduce_factor*k
           k1 = k2 - (reduce_factor-1)
           array(:,i,j,k) = sum(arr(:,i1:i2,j1:j2,k1:k2))/reduce_factor**3
        enddo
        enddo
        enddo
     else
        array(:,:,:,:) = arr(:,:,:,:)
     endif
  else
     array(:,:,:,:) = arr(:,:,:,:)
  endif
  if (allocated(arr)) deallocate(arr)

  end subroutine read_4D
  !=================================================================
  ! get 1D dimeisions from density fits/hdf5 file.
  ! updated to deal with xmax,ymax,zmax (2023.01.20).
  subroutine get_dimension_1D(fname,nz,zmax,status)
  implicit none
  character(len=*), intent(in)    :: fname
  integer,          intent(out)   :: nz
  real(kind=wp),    intent(inout) :: zmax
  integer,          intent(inout) :: status

  !--- local variables
  type(io_file_type) :: iofh
  real(kind=wp) :: zmax1

  status = 0
  call io_open_old(iofh,trim(fname),status)
  if (status == 0) then
     call io_get_keyword(iofh,'NAXIS1',nz,status)
     call io_get_keyword(iofh,'zmax',zmax1,status)
     if (status == 0) then
        zmax = zmax1
     else
        status = 0
     endif
  else
     write(*,*) 'Error in reading the data file :', trim(fname),' status = ',status
  endif
  call io_close(iofh,status)
  end subroutine get_dimension_1D
  !=================================================================
  ! get 3D dimeisions from density fits/hdf5 file.
  ! updated to deal with xmax,ymax,zmax (2023.01.20).
  subroutine get_dimension(fname,nx,ny,nz,xmax,ymax,zmax,status,reduce_factor)
  implicit none
  character(len=*),  intent(in)    :: fname
  integer,           intent(out)   :: nx,ny,nz
  real(kind=wp),     intent(inout) :: xmax,ymax,zmax
  integer,           intent(inout) :: status
  integer, optional, intent(inout) :: reduce_factor

  !--- local variables
  type(io_file_type) :: iofh
  real(kind=wp) :: xmax1, ymax1, zmax1

  status = 0
  call io_open_old(iofh,trim(fname),status)
  if (status == 0) then
     call io_get_keyword(iofh,'NAXIS1',nx,status)
     call io_get_keyword(iofh,'NAXIS2',ny,status)
     call io_get_keyword(iofh,'NAXIS3',nz,status)
     call io_get_keyword(iofh,'xmax',xmax1,status)
     if (status == 0) then
        xmax = xmax1
     else
        status = 0
     endif
     call io_get_keyword(iofh,'ymax',ymax1,status)
     if (status == 0)then
        ymax = ymax1
     else
        status = 0
     endif
     call io_get_keyword(iofh,'zmax',zmax1,status)
     if (status == 0) then
        zmax = zmax1
     else
        status = 0
     endif
     if (present(reduce_factor)) then
        if (reduce_factor > 1) then
           if ((nx/reduce_factor) > 0 .and. (ny/reduce_factor) > 0 .and. (nz/reduce_factor) > 0) then
              nx = nx/reduce_factor
              ny = ny/reduce_factor
              nz = nz/reduce_factor
           else
              reduce_factor = 1
              write(*,*) 'Oops, Using the original size...'
           endif
        endif
     endif
  else
     write(*,*) 'Error in reading the data file :', trim(fname),' status = ',status
  endif
  call io_close(iofh,status)
  end subroutine get_dimension
  !=================================================================
end module read_grid_data
