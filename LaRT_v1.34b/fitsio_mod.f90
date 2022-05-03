module fitsio_mod
!
!-- 2020-11-10, added fits_append_table_column_int and fits_read_table_column_int
!-- 2020-10.18, added fits_get_column_number
!-- 2018-01-15, Null values are now NaN.
!-- 2017-07-13, Added Bin Table routines.
!-- Written by Kwang-il Seon, 2017-06-28
!      Include only some FITS routines that are required in MoCafe and LaRT.
!
  use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int64
  use, intrinsic :: ieee_arithmetic
  implicit none
  private
  public fits_open_new, fits_open_old, fits_close, fits_move_to_next_hdu, fits_move_to_nth_hdu, &
         fits_write_image, fits_append_image, fits_read_image, fits_put_keyword, fits_get_keyword, &
         fits_write_table_column, fits_append_table_column, fits_read_table_column, &
         fits_hdu_copy, fits_get_column_number

  ! note that the default output bitpix is the same as that of input array.
  interface fits_write_image
     module procedure fits_append_1D_real32, fits_append_2D_real32, fits_append_3D_real32, fits_append_4D_real32, &
                      fits_append_1D_real64, fits_append_2D_real64, fits_append_3D_real64, fits_append_4D_real64
  end interface fits_write_image
  interface fits_append_image
     module procedure fits_append_1D_real32, fits_append_2D_real32, fits_append_3D_real32, fits_append_4D_real32, &
                      fits_append_1D_real64, fits_append_2D_real64, fits_append_3D_real64, fits_append_4D_real64, &
                      fits_append_1D_int32,  fits_append_2D_int32,  fits_append_3D_int32,  fits_append_4D_int32
  end interface fits_append_image
  interface fits_read_image
     module procedure fits_read_1D_real32, fits_read_2D_real32, fits_read_3D_real32, fits_read_4D_real32, &
                      fits_read_1D_real64, fits_read_2D_real64, fits_read_3D_real64, fits_read_4D_real64
  end interface fits_read_image

  interface fits_write_table_column
     module procedure fits_append_table_column_real32, fits_append_table_column_real64, &
                      fits_append_table_column_int32,  fits_append_table_column_int64
  end interface fits_write_table_column
  interface fits_append_table_column
     module procedure fits_append_table_column_real32, fits_append_table_column_real64, &
                      fits_append_table_column_int32,  fits_append_table_column_int64
  end interface fits_append_table_column
  interface fits_read_table_column
     module procedure fits_read_table_column_real32, fits_read_table_column_real64, &
                      fits_read_table_column_int32,  fits_read_table_column_int64
  end interface fits_read_table_column

  interface fits_put_keyword
     module procedure fits_put_key_logical, fits_put_key_real32, fits_put_key_real64, &
                      fits_put_key_string,  fits_put_key_int32,  fits_put_key_int64
  end interface fits_put_keyword
  interface fits_get_keyword
     module procedure fits_get_key_logical, fits_get_key_real32, fits_get_key_real64, &
                      fits_get_key_string,  fits_get_key_int32,  fits_get_key_int64
  end interface fits_get_keyword
contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine fits_open_new(unit,fname,status)
  integer,            intent(out) :: unit
  character(len=*),    intent(in) :: fname
  integer,          intent(inout) :: status
  integer :: blocksize=1
  call unlink(trim(fname))
  status = 0
  call ftgiou(unit,status)
  call ftinit(unit,trim(fname),blocksize,status)
  end subroutine fits_open_new
!--------------------
  subroutine fits_open_old(unit,fname,status,overwrite)
  integer,            intent(out) :: unit
  character(len=*),    intent(in) :: fname
  integer,          intent(inout) :: status
  logical, optional,   intent(in) :: overwrite
  integer :: blocksize=1, readwrite=0
  if (present(overwrite)) then
     if (overwrite) readwrite = 1
  endif
  status = 0
  call ftgiou(unit,status)
  call ftopen(unit,trim(fname),readwrite,blocksize,status)
  end subroutine fits_open_old
!--------------------
  subroutine fits_close(unit,status)
  integer, intent(in)    :: unit
  integer, intent(inout) :: status
  call ftclos(unit,status)
  call ftfiou(unit,status)
  end subroutine fits_close
!--------------------
  subroutine fits_hdu_copy(unit1,unit2,status)
  integer, intent(in)    :: unit1,unit2
  integer, intent(inout) :: status
  integer :: morekeys = 0
  call ftcopy(unit1,unit2,morekeys,status)
  end subroutine fits_hdu_copy
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine fits_move_to_next_hdu(unit,status,nstep)
  integer,           intent(in)    :: unit
  integer,           intent(inout) :: status
  integer, optional, intent(in)    :: nstep
  integer :: move_step, hdutype
  move_step = 1
  if (present(nstep)) then
     if (nstep >= 0) move_step = nstep
  endif
  call ftmrhd(unit,move_step,hdutype,status)
  end subroutine fits_move_to_next_hdu
!--------------------
  subroutine fits_move_to_nth_hdu(unit,nth,status)
  integer,    intent(in) :: unit
  integer,    intent(in) :: nth
  integer, intent(inout) :: status
  integer :: hdutype
  call ftmahd(unit,nth,hdutype,status)
  end subroutine fits_move_to_nth_hdu
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine fits_put_key_logical(unit,keyname,keyvalue,comment,status)
  integer,             intent(in) :: unit
  character(len=*),    intent(in) :: keyname
  logical,             intent(in) :: keyvalue
  character(len=*),    intent(in) :: comment
  integer,          intent(inout) :: status
  call ftpkyl(unit,trim(keyname),keyvalue,trim(comment),status)
  end subroutine fits_put_key_logical
!--------------------
  subroutine fits_put_key_real32(unit,keyname,keyvalue,comment,status)
  integer,             intent(in) :: unit
  character(len=*),    intent(in) :: keyname
  real(real32),        intent(in) :: keyvalue
  character(len=*),    intent(in) :: comment
  integer,          intent(inout) :: status
  integer :: decimal=-8
  call ftpkye(unit,trim(keyname),keyvalue,decimal,trim(comment),status)
  end subroutine fits_put_key_real32
!--------------------
  subroutine fits_put_key_real64(unit,keyname,keyvalue,comment,status)
  integer,             intent(in) :: unit
  character(len=*),    intent(in) :: keyname
  real(real64),        intent(in) :: keyvalue
  character(len=*),    intent(in) :: comment
  integer,          intent(inout) :: status
  integer :: decimal=-8
  call ftpkyd(unit,trim(keyname),keyvalue,decimal,trim(comment),status)
  end subroutine fits_put_key_real64
!--------------------
  subroutine fits_put_key_int32(unit,keyname,keyvalue,comment,status)
  integer,             intent(in) :: unit
  character(len=*),    intent(in) :: keyname
  integer(int32),      intent(in) :: keyvalue
  character(len=*),    intent(in) :: comment
  integer,          intent(inout) :: status
  call ftpkyj(unit,trim(keyname),keyvalue,trim(comment),status)
  end subroutine fits_put_key_int32
!--------------------
  subroutine fits_put_key_int64(unit,keyname,keyvalue,comment,status)
  integer,             intent(in) :: unit
  character(len=*),    intent(in) :: keyname
  integer(int64),      intent(in) :: keyvalue
  character(len=*),    intent(in) :: comment
  integer,          intent(inout) :: status
  call ftpkyk(unit,trim(keyname),keyvalue,trim(comment),status)
  end subroutine fits_put_key_int64
!--------------------
  subroutine fits_put_key_string(unit,keyname,keyvalue,comment,status)
  integer,             intent(in) :: unit
  character(len=*),    intent(in) :: keyname
  character(len=*),    intent(in) :: keyvalue
  character(len=*),    intent(in) :: comment
  integer,          intent(inout) :: status
  call ftpkys(unit,trim(keyname),trim(keyvalue),trim(comment),status)
  end subroutine fits_put_key_string
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine fits_get_key_logical(unit,keyname,keyvalue,status,comment)
  integer,                     intent(in) :: unit
  character(len=*),            intent(in) :: keyname
  logical,                    intent(out) :: keyvalue
  integer,                  intent(inout) :: status
  character(len=*), optional, intent(out) :: comment
  character(len=72) :: comment0
  call ftgkyl(unit,trim(keyname),keyvalue,comment0,status)
  if (present(comment)) comment = trim(comment0)
  end subroutine fits_get_key_logical
!--------------------
  subroutine fits_get_key_real32(unit,keyname,keyvalue,status,comment)
  integer,                     intent(in) :: unit
  character(len=*),            intent(in) :: keyname
  real(real32),               intent(out) :: keyvalue
  integer,                  intent(inout) :: status
  character(len=*), optional, intent(out) :: comment
  character(len=72) :: comment0
  call ftgkye(unit,trim(keyname),keyvalue,comment0,status)
  if (present(comment)) comment = trim(comment0)
  end subroutine fits_get_key_real32
!--------------------
  subroutine fits_get_key_real64(unit,keyname,keyvalue,status,comment)
  integer,                     intent(in) :: unit
  character(len=*),            intent(in) :: keyname
  real(real64),               intent(out) :: keyvalue
  integer,                  intent(inout) :: status
  character(len=*), optional, intent(out) :: comment
  character(len=72) :: comment0
  call ftgkyd(unit,trim(keyname),keyvalue,comment0,status)
  if (present(comment)) comment = trim(comment0)
  end subroutine fits_get_key_real64
!--------------------
  subroutine fits_get_key_int32(unit,keyname,keyvalue,status,comment)
  integer,                     intent(in) :: unit
  character(len=*),            intent(in) :: keyname
  integer(int32),             intent(out) :: keyvalue
  integer,                  intent(inout) :: status
  character(len=*), optional, intent(out) :: comment
  character(len=72) :: comment0
  call ftgkyj(unit,trim(keyname),keyvalue,comment0,status)
  if (present(comment)) comment = trim(comment0)
  end subroutine fits_get_key_int32
!--------------------
  subroutine fits_get_key_int64(unit,keyname,keyvalue,status,comment)
  integer,                     intent(in) :: unit
  character(len=*),            intent(in) :: keyname
  integer(int64),             intent(out) :: keyvalue
  integer,                  intent(inout) :: status
  character(len=*), optional, intent(out) :: comment
  character(len=72) :: comment0
  call ftgkyk(unit,trim(keyname),keyvalue,comment0,status)
  if (present(comment)) comment = trim(comment0)
  end subroutine fits_get_key_int64
!--------------------
  subroutine fits_get_key_string(unit,keyname,keyvalue,status,comment)
  integer,                     intent(in) :: unit
  character(len=*),            intent(in) :: keyname
  character(len=*),           intent(out) :: keyvalue
  integer,                  intent(inout) :: status
  character(len=*), optional, intent(out) :: comment
  character(len=72) :: comment0
  call ftgkys(unit,trim(keyname),keyvalue,comment0,status)
  if (present(comment)) comment = trim(comment0)
  end subroutine fits_get_key_string
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine fits_read_1D_real32(unit,array,status)
  integer,      intent(in)  :: unit
  real(real32), intent(out) :: array(:)
  integer,    intent(inout) :: status

  !--- local variables
  real(real32) :: nullval32
  real(real64) :: nullval64
  logical  :: ltemp
  integer  :: bitpix, group=1, fpixel, n1
  character(len=72) :: comment
  real(real64), allocatable :: arr(:)

  nullval32 = ieee_value(0.0_real32,ieee_quiet_nan)
  nullval64 = ieee_value(0.0_real64,ieee_quiet_nan)
  fpixel = 1
  n1     = size(array)

  call ftgkyj(unit,'bitpix',bitpix,comment,status)
  if (bitpix /= -32) then
     if (.not. allocated(arr)) allocate(arr(n1))
     call ftgpvd(unit,group,fpixel,n1,nullval64,arr,ltemp,status)
     array = arr
     if (allocated(arr)) deallocate(arr)
  else
     call ftgpve(unit,group,fpixel,n1,nullval32,array,ltemp,status)
  endif
  end subroutine fits_read_1D_real32
!--------------------
  subroutine fits_read_2D_real32(unit,array,status)
  integer,      intent(in)  :: unit
  real(real32), intent(out) :: array(:,:)
  integer,    intent(inout) :: status

  !--- local variables
  real(real32) :: nullval32
  real(real64) :: nullval64
  logical  :: ltemp
  integer  :: bitpix, group=1, fpixel, n1, n2, j
  character(len=72) :: comment
  real(real64), allocatable :: arr(:,:)

  nullval32 = ieee_value(0.0_real32,ieee_quiet_nan)
  nullval64 = ieee_value(0.0_real64,ieee_quiet_nan)
  fpixel = 1
  n1     = size(array,1)
  n2     = size(array,2)

  call ftgkyj(unit,'bitpix',bitpix,comment,status)
  if (bitpix /= -32) then
     if (.not. allocated(arr)) allocate(arr(n1,n2))
     do j=1, n2
        call ftgpvd(unit,group,fpixel,n1,nullval64,arr(:,j),ltemp,status)
        fpixel = fpixel + n1
     enddo
     array = arr
     if (allocated(arr)) deallocate(arr)
  else
     do j=1, n2
        call ftgpve(unit,group,fpixel,n1,nullval32,array(:,j),ltemp,status)
        fpixel = fpixel + n1
     enddo
  endif
  end subroutine fits_read_2D_real32
!--------------------
  subroutine fits_read_3D_real32(unit,array,status)
  integer,      intent(in)  :: unit
  real(real32), intent(out) :: array(:,:,:)
  integer,    intent(inout) :: status

  !--- local variables
  real(real32) :: nullval32
  real(real64) :: nullval64
  logical  :: ltemp
  integer  :: bitpix, group=1, fpixel, n1, n2, n3, j, k
  character(len=72) :: comment
  real(real64), allocatable :: arr(:,:,:)

  nullval32 = ieee_value(0.0_real32,ieee_quiet_nan)
  nullval64 = ieee_value(0.0_real64,ieee_quiet_nan)
  fpixel = 1
  n1     = size(array,1)
  n2     = size(array,2)
  n3     = size(array,3)

  call ftgkyj(unit,'bitpix',bitpix,comment,status)
  if (bitpix /= -32) then
     if (.not. allocated(arr)) allocate(arr(n1,n2,n3))
     do k=1, n3
     do j=1, n2
        call ftgpvd(unit,group,fpixel,n1,nullval64,arr(:,j,k),ltemp,status)
        fpixel = fpixel + n1
     enddo
     enddo
     array = arr
     if (allocated(arr)) deallocate(arr)
  else
     do k=1, n3
     do j=1, n2
        call ftgpve(unit,group,fpixel,n1,nullval32,array(:,j,k),ltemp,status)
        fpixel = fpixel + n1
     enddo
     enddo
  endif
  end subroutine fits_read_3D_real32
!--------------------
  subroutine fits_read_4D_real32(unit,array,status)
  integer,      intent(in)  :: unit
  real(real32), intent(out) :: array(:,:,:,:)
  integer,    intent(inout) :: status

  !--- local variables
  real(real32) :: nullval32
  real(real64) :: nullval64
  logical  :: ltemp
  integer  :: bitpix, group=1, fpixel, n1, n2, n3, n4, j, k, l
  character(len=72) :: comment
  real(real64), allocatable :: arr(:,:,:,:)

  nullval32 = ieee_value(0.0_real32,ieee_quiet_nan)
  nullval64 = ieee_value(0.0_real64,ieee_quiet_nan)
  fpixel = 1
  n1     = size(array,1)
  n2     = size(array,2)
  n3     = size(array,3)
  n4     = size(array,4)

  call ftgkyj(unit,'bitpix',bitpix,comment,status)
  if (bitpix /= -32) then
     if (.not. allocated(arr)) allocate(arr(n1,n2,n3,n4))
     do l=1, n4
     do k=1, n3
     do j=1, n2
        call ftgpvd(unit,group,fpixel,n1,nullval64,arr(:,j,k,l),ltemp,status)
        fpixel = fpixel + n1
     enddo
     enddo
     enddo
     array = arr
     if (allocated(arr)) deallocate(arr)
  else
     do l=1, n4
     do k=1, n3
     do j=1, n2
        call ftgpve(unit,group,fpixel,n1,nullval32,array(:,j,k,l),ltemp,status)
        fpixel = fpixel + n1
     enddo
     enddo
     enddo
  endif
  end subroutine fits_read_4D_real32
!--------------------
  subroutine fits_read_1D_real64(unit,array,status)
  integer,      intent(in)  :: unit
  real(real64), intent(out) :: array(:)
  integer,    intent(inout) :: status

  !--- local variables
  real(real32) :: nullval32
  real(real64) :: nullval64
  logical  :: ltemp
  integer  :: bitpix, group=1, fpixel, n1
  character(len=72) :: comment
  real(real32), allocatable :: arr(:)

  nullval32 = ieee_value(0.0_real32,ieee_quiet_nan)
  nullval64 = ieee_value(0.0_real64,ieee_quiet_nan)
  fpixel = 1
  n1     = size(array)

  call ftgkyj(unit,'bitpix',bitpix,comment,status)
  if (bitpix /= -64) then
     if (.not. allocated(arr)) allocate(arr(n1))
     call ftgpve(unit,group,fpixel,n1,nullval32,arr,ltemp,status)
     array = arr
     if (allocated(arr)) deallocate(arr)
  else
     call ftgpvd(unit,group,fpixel,n1,nullval64,array,ltemp,status)
  endif
  end subroutine fits_read_1D_real64
!--------------------
  subroutine fits_read_2D_real64(unit,array,status)
  integer,      intent(in)  :: unit
  real(real64), intent(out) :: array(:,:)
  integer,    intent(inout) :: status

  !--- local variables
  real(real32) :: nullval32
  real(real64) :: nullval64
  logical  :: ltemp
  integer  :: bitpix, group=1, fpixel, n1, n2, j
  character(len=72) :: comment
  real(real32), allocatable :: arr(:,:)

  nullval32 = ieee_value(0.0_real32,ieee_quiet_nan)
  nullval64 = ieee_value(0.0_real64,ieee_quiet_nan)
  fpixel = 1
  n1     = size(array,1)
  n2     = size(array,2)

  call ftgkyj(unit,'bitpix',bitpix,comment,status)
  if (bitpix /= -64) then
     if (.not. allocated(arr)) allocate(arr(n1,n2))
     do j=1, n2
        call ftgpve(unit,group,fpixel,n1,nullval32,arr(:,j),ltemp,status)
        fpixel = fpixel + n1
     enddo
     array = arr
     if (allocated(arr)) deallocate(arr)
  else
     do j=1, n2
        call ftgpvd(unit,group,fpixel,n1,nullval64,array(:,j),ltemp,status)
        fpixel = fpixel + n1
     enddo
  endif
  end subroutine fits_read_2D_real64
!--------------------
  subroutine fits_read_3D_real64(unit,array,status)
  integer,      intent(in)  :: unit
  real(real64), intent(out) :: array(:,:,:)
  integer,    intent(inout) :: status

  !--- local variables
  real(real32) :: nullval32
  real(real64) :: nullval64
  logical  :: ltemp
  integer  :: bitpix, group=1, fpixel, n1, n2, n3, j, k
  character(len=72) :: comment
  real(real32), allocatable :: arr(:,:,:)

  nullval32 = ieee_value(0.0_real32,ieee_quiet_nan)
  nullval64 = ieee_value(0.0_real64,ieee_quiet_nan)
  fpixel = 1
  n1     = size(array,1)
  n2     = size(array,2)
  n3     = size(array,3)

  call ftgkyj(unit,'bitpix',bitpix,comment,status)
  if (bitpix /= -64) then
     if (.not. allocated(arr)) allocate(arr(n1,n2,n3))
     do k=1, n3
     do j=1, n2
        call ftgpve(unit,group,fpixel,n1,nullval32,arr(:,j,k),ltemp,status)
        fpixel = fpixel + n1
     enddo
     enddo
     array = arr
     if (allocated(arr)) deallocate(arr)
  else
     do k=1, n3
     do j=1, n2
        call ftgpvd(unit,group,fpixel,n1,nullval64,array(:,j,k),ltemp,status)
        fpixel = fpixel + n1
     enddo
     enddo
  endif
  end subroutine fits_read_3D_real64
!--------------------
  subroutine fits_read_4D_real64(unit,array,status)
  integer,      intent(in)  :: unit
  real(real64), intent(out) :: array(:,:,:,:)
  integer,    intent(inout) :: status

  !--- local variables
  real(real32) :: nullval32
  real(real64) :: nullval64
  logical  :: ltemp
  integer  :: bitpix, group=1, fpixel, n1, n2, n3, n4, j, k, l
  character(len=72) :: comment
  real(real32), allocatable :: arr(:,:,:,:)

  nullval32 = ieee_value(0.0_real32,ieee_quiet_nan)
  nullval64 = ieee_value(0.0_real64,ieee_quiet_nan)
  fpixel = 1
  n1     = size(array,1)
  n2     = size(array,2)
  n3     = size(array,3)
  n4     = size(array,4)

  call ftgkyj(unit,'bitpix',bitpix,comment,status)
  if (bitpix /= -64) then
     if (.not. allocated(arr)) allocate(arr(n1,n2,n3,n4))
     do l=1, n4
     do k=1, n3
     do j=1, n2
        call ftgpve(unit,group,fpixel,n1,nullval32,arr(:,j,k,l),ltemp,status)
        fpixel = fpixel + n1
     enddo
     enddo
     enddo
     array = arr
     if (allocated(arr)) deallocate(arr)
  else
     do l=1, n4
     do k=1, n3
     do j=1, n2
        call ftgpvd(unit,group,fpixel,n1,nullval64,array(:,j,k,l),ltemp,status)
        fpixel = fpixel + n1
     enddo
     enddo
     enddo
  endif
  end subroutine fits_read_4D_real64
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine fits_append_1D_real32(unit,array,status,bitpix)
  integer,      intent(in) :: unit
  real(real32), intent(in) :: array(:)
  integer,           intent(inout) :: status
  integer, optional, intent(in)    :: bitpix

  !--- local variables
  integer :: bitpix_out
  integer :: naxis
  integer(int64) :: naxes(1)
  integer :: group=1, felem=1, nelem
  real(real64), allocatable :: arr(:)

  if (.not. present(bitpix)) then
     bitpix_out = -32
  else
     bitpix_out = bitpix
  endif

  naxis    = 1
  naxes(1) = size(array,1)
  nelem    = size(array)
  call ftiimgll(unit,bitpix_out,naxis,naxes(1:naxis),status)
  if (bitpix_out /= -32) then
     if (.not. allocated(arr)) allocate(arr(naxes(1)))
     arr = array
     call ftpprd(unit,group,felem,nelem,arr,status)
     if (allocated(arr)) deallocate(arr)
  else
     call ftppre(unit,group,felem,nelem,array,status)
  endif
  end subroutine fits_append_1D_real32
!--------------------
  subroutine fits_append_2D_real32(unit,array,status,bitpix)
  integer,      intent(in) :: unit
  real(real32), intent(in) :: array(:,:)
  integer,           intent(inout) :: status
  integer, optional, intent(in)    :: bitpix

  !--- local variables
  integer :: bitpix_out
  integer :: naxis
  integer(int64) :: naxes(2)
  integer :: group=1, felem=1, nelem
  real(real64), allocatable :: arr(:,:)

  if (.not. present(bitpix)) then
     bitpix_out = -32
  else
     bitpix_out = bitpix
  endif

  naxis    = 2
  naxes(1) = size(array,1)
  naxes(2) = size(array,2)
  nelem    = size(array)
  call ftiimgll(unit,bitpix_out,naxis,naxes(1:naxis),status)
  if (bitpix_out /= -32) then
     if (.not. allocated(arr)) allocate(arr(naxes(1),naxes(2)))
     arr = array
     call ftpprd(unit,group,felem,nelem,arr,status)
     if (allocated(arr)) deallocate(arr)
  else
     call ftppre(unit,group,felem,nelem,array,status)
  endif
  end subroutine fits_append_2D_real32
!--------------------
  subroutine fits_append_3D_real32(unit,array,status,bitpix)
  integer,      intent(in) :: unit
  real(real32), intent(in) :: array(:,:,:)
  integer,           intent(inout) :: status
  integer, optional, intent(in)    :: bitpix

  !--- local variables
  integer :: bitpix_out
  integer :: naxis
  integer(int64) :: naxes(3)
  integer :: group=1, felem=1, nelem
  real(real64), allocatable :: arr(:,:,:)

  if (.not. present(bitpix)) then
     bitpix_out = -32
  else
     bitpix_out = bitpix
  endif

  naxis    = 3
  naxes(1) = size(array,1)
  naxes(2) = size(array,2)
  naxes(3) = size(array,3)
  if (naxes(1) == 1 .and. naxes(2) == 1) then
     naxis      = 1
     naxes(1)   = naxes(3)
     naxes(2:3) = 1
  endif
  nelem    = size(array)
  call ftiimgll(unit,bitpix_out,naxis,naxes(1:naxis),status)
  if (bitpix_out /= -32) then
     if (.not. allocated(arr)) allocate(arr(naxes(1),naxes(2),naxes(3)))
     arr = array
     call ftpprd(unit,group,felem,nelem,arr,status)
     if (allocated(arr)) deallocate(arr)
  else
     call ftppre(unit,group,felem,nelem,array,status)
  endif
  end subroutine fits_append_3D_real32
!--------------------
  subroutine fits_append_4D_real32(unit,array,status,bitpix)
  integer,      intent(in) :: unit
  real(real32), intent(in) :: array(:,:,:,:)
  integer,           intent(inout) :: status
  integer, optional, intent(in)    :: bitpix

  !--- local variables
  integer :: bitpix_out
  integer :: naxis
  integer(int64) :: naxes(4)
  integer :: group=1, felem=1, nelem
  real(real64), allocatable :: arr(:,:,:,:)

  if (.not. present(bitpix)) then
     bitpix_out = -32
  else
     bitpix_out = bitpix
  endif

  naxis    = 4
  naxes(1) = size(array,1)
  naxes(2) = size(array,2)
  naxes(3) = size(array,3)
  naxes(4) = size(array,4)
  if (naxes(2) == 1 .and. naxes(3) == 1) then
     naxis      = 2
     naxes(2)   = naxes(4)
     naxes(3:4) = 1
  endif
  nelem    = size(array)
  call ftiimgll(unit,bitpix_out,naxis,naxes(1:naxis),status)
  if (bitpix_out /= -32) then
     if (.not. allocated(arr)) allocate(arr(naxes(1),naxes(2),naxes(3),naxes(4)))
     arr = array
     call ftpprd(unit,group,felem,nelem,arr,status)
     if (allocated(arr)) deallocate(arr)
  else
     call ftppre(unit,group,felem,nelem,array,status)
  endif
  end subroutine fits_append_4D_real32
!--------------------
  subroutine fits_append_1D_real64(unit,array,status,bitpix)
  integer,      intent(in) :: unit
  real(real64), intent(in) :: array(:)
  integer,           intent(inout) :: status
  integer, optional, intent(in)    :: bitpix

  !--- local variables
  integer :: bitpix_out
  integer :: naxis
  integer(int64) :: naxes(1)
  integer :: group=1, felem=1, nelem
  real(real32), allocatable :: arr(:)

  if (.not. present(bitpix)) then
     bitpix_out = -64
  else
     bitpix_out = bitpix
  endif
  naxis    = 1
  naxes(1) = size(array,1)
  nelem    = size(array)
  call ftiimgll(unit,bitpix_out,naxis,naxes(1:naxis),status)
  if (bitpix_out /= -32) then
     call ftpprd(unit,group,felem,nelem,array,status)
  else
     if (.not. allocated(arr)) allocate(arr(naxes(1)))
     arr = array
     call ftppre(unit,group,felem,nelem,arr,status)
     if (allocated(arr)) deallocate(arr)
  endif
  end subroutine fits_append_1D_real64
!--------------------
  subroutine fits_append_2D_real64(unit,array,status,bitpix)
  integer,      intent(in) :: unit
  real(real64), intent(in) :: array(:,:)
  integer,           intent(inout) :: status
  integer, optional, intent(in)    :: bitpix

  !--- local variables
  integer :: bitpix_out
  integer :: naxis
  integer(int64) :: naxes(2)
  integer :: group=1, felem=1, nelem
  real(real32), allocatable :: arr(:,:)

  if (.not. present(bitpix)) then
     bitpix_out = -64
  else
     bitpix_out = bitpix
  endif

  naxis    = 2
  naxes(1) = size(array,1)
  naxes(2) = size(array,2)
  nelem    = size(array)
  call ftiimgll(unit,bitpix_out,naxis,naxes(1:naxis),status)
  if (bitpix_out /= -32) then
     call ftpprd(unit,group,felem,nelem,array,status)
  else
     if (.not. allocated(arr)) allocate(arr(naxes(1),naxes(2)))
     arr = array
     call ftppre(unit,group,felem,nelem,arr,status)
     if (allocated(arr)) deallocate(arr)
  endif
  end subroutine fits_append_2D_real64
!--------------------
  subroutine fits_append_3D_real64(unit,array,status,bitpix)
  integer,      intent(in) :: unit
  real(real64), intent(in) :: array(:,:,:)
  integer,           intent(inout) :: status
  integer, optional, intent(in)    :: bitpix

  !--- local variables
  integer :: bitpix_out
  integer :: naxis
  integer(int64) :: naxes(3)
  integer :: group=1, felem=1, nelem
  real(real32), allocatable :: arr(:,:,:)

  if (.not. present(bitpix)) then
     bitpix_out = -64
  else
     bitpix_out = bitpix
  endif

  naxis    = 3
  naxes(1) = size(array,1)
  naxes(2) = size(array,2)
  naxes(3) = size(array,3)
  if (naxes(1) == 1 .and. naxes(2) == 1) then
     naxis      = 1
     naxes(1)   = naxes(3)
     naxes(2:3) = 1
  endif
  nelem    = size(array)
  call ftiimgll(unit,bitpix_out,naxis,naxes(1:naxis),status)
  if (bitpix_out /= -32) then
     call ftpprd(unit,group,felem,nelem,array,status)
  else
     if (.not. allocated(arr)) allocate(arr(naxes(1),naxes(2),naxes(3)))
     arr = array
     call ftppre(unit,group,felem,nelem,arr,status)
     if (allocated(arr)) deallocate(arr)
  endif
  end subroutine fits_append_3D_real64
!--------------------
  subroutine fits_append_4D_real64(unit,array,status,bitpix)
  integer,      intent(in) :: unit
  real(real64), intent(in) :: array(:,:,:,:)
  integer,           intent(inout) :: status
  integer, optional, intent(in)    :: bitpix

  !--- local variables
  integer :: bitpix_out
  integer :: naxis
  integer(int64) :: naxes(4)
  integer :: group=1, felem=1, nelem
  real(real32), allocatable :: arr(:,:,:,:)

  if (.not. present(bitpix)) then
     bitpix_out = -64
  else
     bitpix_out = bitpix
  endif

  naxis    = 4
  naxes(1) = size(array,1)
  naxes(2) = size(array,2)
  naxes(3) = size(array,3)
  naxes(4) = size(array,4)
  if (naxes(2) == 1 .and. naxes(3) == 1) then
     naxis      = 2
     naxes(2)   = naxes(4)
     naxes(3:4) = 1
  endif
  nelem    = size(array)
  call ftiimgll(unit,bitpix_out,naxis,naxes(1:naxis),status)
  if (bitpix_out /= -32) then
     call ftpprd(unit,group,felem,nelem,array,status)
  else
     if (.not. allocated(arr)) allocate(arr(naxes(1),naxes(2),naxes(3),naxes(4)))
     arr = array
     call ftppre(unit,group,felem,nelem,arr,status)
     if (allocated(arr)) deallocate(arr)
  endif
  end subroutine fits_append_4D_real64
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine fits_append_1D_int32(unit,array,status)
  integer, intent(in) :: unit
  integer, intent(in) :: array(:)
  integer, intent(inout) :: status

  !--- local variables
  integer :: bitpix_out
  integer :: naxis
  integer(int64) :: naxes(1)
  integer :: group=1, felem=1, nelem

  bitpix_out = 32

  naxis    = 1
  naxes(1) = size(array,1)
  nelem    = size(array)
  call ftiimgll(unit,bitpix_out,naxis,naxes(1:naxis),status)
  call ftpprj(unit,group,felem,nelem,array,status)
  end subroutine fits_append_1D_int32
!--------------------
  subroutine fits_append_2D_int32(unit,array,status)
  integer, intent(in) :: unit
  integer, intent(in) :: array(:,:)
  integer, intent(inout) :: status

  !--- local variables
  integer :: bitpix_out
  integer :: naxis
  integer(int64) :: naxes(2)
  integer :: group=1, felem=1, nelem

  bitpix_out = 32

  naxis    = 2
  naxes(1) = size(array,1)
  naxes(2) = size(array,2)
  nelem    = size(array)
  call ftiimgll(unit,bitpix_out,naxis,naxes(1:naxis),status)
  call ftpprj(unit,group,felem,nelem,array,status)
  end subroutine fits_append_2D_int32
!--------------------
  subroutine fits_append_3D_int32(unit,array,status)
  integer, intent(in) :: unit
  integer, intent(in) :: array(:,:,:)
  integer, intent(inout) :: status

  !--- local variables
  integer :: bitpix_out
  integer :: naxis
  integer(int64) :: naxes(3)
  integer :: group=1, felem=1, nelem

  bitpix_out = 32

  naxis    = 3
  naxes(1) = size(array,1)
  naxes(2) = size(array,2)
  naxes(3) = size(array,3)
  nelem    = size(array)
  call ftiimgll(unit,bitpix_out,naxis,naxes(1:naxis),status)
  call ftpprj(unit,group,felem,nelem,array,status)
  end subroutine fits_append_3D_int32
!--------------------
  subroutine fits_append_4D_int32(unit,array,status)
  integer, intent(in) :: unit
  integer, intent(in) :: array(:,:,:,:)
  integer, intent(inout) :: status

  !--- local variables
  integer :: bitpix_out
  integer :: naxis
  integer(int64) :: naxes(4)
  integer :: group=1, felem=1, nelem

  bitpix_out = 32

  naxis    = 4
  naxes(1) = size(array,1)
  naxes(2) = size(array,2)
  naxes(3) = size(array,3)
  naxes(4) = size(array,4)
  nelem    = size(array)
  call ftiimgll(unit,bitpix_out,naxis,naxes(1:naxis),status)
  call ftpprj(unit,group,felem,nelem,array,status)
  end subroutine fits_append_4D_int32
!--------------------
  subroutine fits_get_column_number(unit,colname,colnum,status)
  integer,           intent(in) :: unit
  character(len=*),  intent(in) :: colname
  integer,          intent(out) :: colnum, status
  logical :: casesen
  casesen = .false.
  call ftgcno(unit, casesen, trim(colname), colnum, status)
  end subroutine fits_get_column_number
!--------------------
  subroutine fits_append_table_column_int32(unit,ttype,array,status)
  integer,           intent(in) :: unit
  character(len=*),  intent(in) :: ttype
  integer(int32),    intent(in) :: array(:)
  integer,        intent(inout) :: status

  integer :: bitpix_out, colnum, nrows, frow, felem
  character(len=16) :: tform, tunit = ' ', extname = ' '
  integer :: tfields = 1, varidat = 0
  integer :: hdutype

  tform  = '1J'
  nrows  = size(array)

  call ftghdt(unit,hdutype,status)
  if (hdutype /= 2) then
     ! create a new BinTable HDU and move to it.
     call ftibin(unit,nrows,tfields,ttype,tform,tunit,extname,varidat,status)
     colnum = 1
  else
     call ftgncl(unit,colnum,status)
     colnum = colnum + 1
     call fticol(unit,colnum,ttype,tform,status)
  endif

  frow  = 1
  felem = 1
  call ftpclj(unit,colnum,frow,felem,nrows,array,status)
  end subroutine fits_append_table_column_int32
!--------------------
  subroutine fits_append_table_column_int64(unit,ttype,array,status)
  integer,           intent(in) :: unit
  character(len=*),  intent(in) :: ttype
  integer(int64),    intent(in) :: array(:)
  integer,        intent(inout) :: status

  integer :: bitpix_out, colnum, nrows, frow, felem
  character(len=16) :: tform, tunit = ' ', extname = ' '
  integer :: tfields = 1, varidat = 0
  integer :: hdutype

  tform  = '1K'
  nrows  = size(array)

  call ftghdt(unit,hdutype,status)
  if (hdutype /= 2) then
     ! create a new BinTable HDU and move to it.
     call ftibin(unit,nrows,tfields,ttype,tform,tunit,extname,varidat,status)
     colnum = 1
  else
     call ftgncl(unit,colnum,status)
     colnum = colnum + 1
     call fticol(unit,colnum,ttype,tform,status)
  endif

  frow  = 1
  felem = 1
  call ftpclk(unit,colnum,frow,felem,nrows,array,status)
  end subroutine fits_append_table_column_int64
!--------------------
  subroutine fits_append_table_column_real64(unit,ttype,array,status,bitpix)
  integer,           intent(in) :: unit
  character(len=*),  intent(in) :: ttype
  real(real64),      intent(in) :: array(:)
  integer,        intent(inout) :: status
  integer, optional, intent(in) :: bitpix

  integer :: bitpix_out, colnum, nrows, frow, felem
  character(len=16) :: tform, tunit = ' ', extname = ' '
  integer :: tfields = 1, varidat = 0
  real(real32), allocatable :: arr(:)
  integer :: hdutype

  if (.not. present(bitpix)) then
     bitpix_out = -64
  else
     bitpix_out = bitpix
  endif
  if (bitpix_out == -64) then
     tform      = '1D'
  else
     tform      = '1E'
  endif
  nrows  = size(array)

  call ftghdt(unit,hdutype,status)
  if (hdutype /= 2) then
     ! create a new BinTable HDU and move to it.
     call ftibin(unit,nrows,tfields,ttype,tform,tunit,extname,varidat,status)
     colnum = 1
  else
     call ftgncl(unit,colnum,status)
     colnum = colnum + 1
     call fticol(unit,colnum,ttype,tform,status)
  endif

  frow  = 1
  felem = 1
  if (bitpix_out == -64) then
     call ftpcld(unit,colnum,frow,felem,nrows,array,status)
  else
     if (.not. allocated(arr)) allocate(arr(nrows))
     arr = array
     call ftpcle(unit,colnum,frow,felem,nrows,arr,status)
     if (allocated(arr)) deallocate(arr)
  endif
  end subroutine fits_append_table_column_real64
!--------------------
  subroutine fits_append_table_column_real32(unit,ttype,array,status,bitpix)
  integer,           intent(in) :: unit
  character(len=*),  intent(in) :: ttype
  real(real32),      intent(in) :: array(:)
  integer,        intent(inout) :: status
  integer, optional, intent(in) :: bitpix

  integer :: bitpix_out, colnum, nrows, frow, felem
  character(len=16) :: tform, tunit = ' ', extname = ' '
  integer :: tfields = 1, varidat = 0
  real(real64), allocatable :: arr(:)
  integer :: hdutype

  if (.not. present(bitpix)) then
     bitpix_out = -32
  else
     bitpix_out = bitpix
  endif
  if (bitpix_out == -32) then
     tform      = '1E'
  else
     tform      = '1D'
  endif
  nrows  = size(array)

  call ftghdt(unit,hdutype,status)
  if (hdutype /= 2) then
     ! create a new BinTable HDU and move to it.
     call ftibin(unit,nrows,tfields,ttype,tform,tunit,extname,varidat,status)
     colnum = 1
  else
     call ftgncl(unit,colnum,status)
     colnum = colnum + 1
     call fticol(unit,colnum,ttype,tform,status)
  endif

  frow  = 1
  felem = 1
  if (bitpix_out == -32) then
     call ftpcle(unit,colnum,frow,felem,nrows,array,status)
  else
     if (.not. allocated(arr)) allocate(arr(nrows))
     arr = array
     call ftpcld(unit,colnum,frow,felem,nrows,arr,status)
     if (allocated(arr)) deallocate(arr)
  endif
  end subroutine fits_append_table_column_real32
!--------------------
  subroutine fits_read_table_column_real32(unit,colnum,array,status)
  integer,           intent(in) :: unit
  integer,           intent(in) :: colnum
  real(real32),     intent(out) :: array(:)
  integer,        intent(inout) :: status

  real(real64), allocatable :: arr(:)
  character(len=16) :: tform
  integer :: numkeys=1, nfound, frow=1, felem=1, nrows
  logical :: anynull
  real(real32) :: null32
  real(real64) :: null64

  null32 = ieee_value(0.0_real32,ieee_quiet_nan)
  null64 = ieee_value(0.0_real64,ieee_quiet_nan)
  call ftgkns(unit,'TFORM',colnum,numkeys,tform,nfound,status)
  nrows = size(array)
  if (trim(tform) == '1E') then
     call ftgcve(unit,colnum,frow,felem,nrows,null32,array,anynull,status)
  else
     if (.not. allocated(arr)) allocate(arr(nrows))
     call ftgcvd(unit,colnum,frow,felem,nrows,null64,arr,anynull,status)
     array = arr
     if (allocated(arr)) deallocate(arr)
  endif
  end subroutine fits_read_table_column_real32
!--------------------
  subroutine fits_read_table_column_real64(unit,colnum,array,status)
  integer,           intent(in) :: unit
  integer,           intent(in) :: colnum
  real(real64),     intent(out) :: array(:)
  integer,        intent(inout) :: status

  real(real32), allocatable :: arr(:)
  character(len=16) :: tform
  integer :: numkeys=1, nfound, frow=1, felem=1, nrows
  logical :: anynull
  real(real32) :: null32
  real(real64) :: null64

  null32 = ieee_value(0.0_real32,ieee_quiet_nan)
  null64 = ieee_value(0.0_real64,ieee_quiet_nan)
  call ftgkns(unit,'TFORM',colnum,numkeys,tform,nfound,status)
  nrows = size(array)
  if (trim(tform) == '1D') then
     call ftgcvd(unit,colnum,frow,felem,nrows,null32,array,anynull,status)
  else
     if (.not. allocated(arr)) allocate(arr(nrows))
     call ftgcve(unit,colnum,frow,felem,nrows,null64,arr,anynull,status)
     array = arr
     if (allocated(arr)) deallocate(arr)
  endif
  end subroutine fits_read_table_column_real64
!--------------------
  subroutine fits_read_table_column_int32(unit,colnum,array,status)
  integer,           intent(in) :: unit
  integer,           intent(in) :: colnum
  integer(int32),   intent(out) :: array(:)
  integer,        intent(inout) :: status

  real(real64), allocatable :: arr(:)
  character(len=16) :: tform
  integer :: numkeys=1, nfound, frow=1, felem=1, nrows
  logical :: anynull
  integer :: nullval = 0

  call ftgkns(unit,'TFORM',colnum,numkeys,tform,nfound,status)
  nrows = size(array)
  if (trim(tform) == '1J') then
     call ftgcvj(unit,colnum,frow,felem,nrows,nullval,array,anynull,status)
  else
     if (.not. allocated(arr)) allocate(arr(nrows))
     call ftgcvk(unit,colnum,frow,felem,nrows,nullval,arr,anynull,status)
     array = arr
     if (allocated(arr)) deallocate(arr)
  endif
  end subroutine fits_read_table_column_int32
!--------------------
  subroutine fits_read_table_column_int64(unit,colnum,array,status)
  integer,           intent(in) :: unit
  integer,           intent(in) :: colnum
  integer(int64),   intent(out) :: array(:)
  integer,        intent(inout) :: status

  real(real32), allocatable :: arr(:)
  character(len=16) :: tform
  integer :: numkeys=1, nfound, frow=1, felem=1, nrows
  logical :: anynull
  integer :: nullval = 0

  call ftgkns(unit,'TFORM',colnum,numkeys,tform,nfound,status)
  nrows = size(array)
  if (trim(tform) == '1K') then
     call ftgcvk(unit,colnum,frow,felem,nrows,nullval,array,anynull,status)
  else
     if (.not. allocated(arr)) allocate(arr(nrows))
     call ftgcvj(unit,colnum,frow,felem,nrows,nullval,arr,anynull,status)
     array = arr
     if (allocated(arr)) deallocate(arr)
  endif
  end subroutine fits_read_table_column_int64
!--------------------
end module fitsio_mod
