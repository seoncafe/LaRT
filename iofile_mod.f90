module iofile_mod
!
! I/O facade for LaRT. Callers never reference FITS or HDF5 libraries
! directly. Backend selection is per io_file_type instance, set by
! io_open_new/old based on filename extension (`.h5`/`.hdf5` → HDF5,
! else FITS). HDF5 backend is compiled in when -DHDF5 is set; otherwise
! HDF5 paths return IO_ERR_NO_HDF5 via the hdf5io_mod stubs.
!
! See docs/io_layout.md for the full specification.
!
  use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int64
  use fitsio_mod
  use hdf5io_mod
  implicit none
  private

  integer, parameter, public :: IO_FMT_FITS = 1
  integer, parameter, public :: IO_FMT_HDF5 = 2

  integer, parameter, public :: IO_ERR_NO_HDF5     = -8001
  integer, parameter, public :: IO_ERR_UNKNOWN_FMT = -8002

  type, public :: io_file_type
     integer            :: format    = IO_FMT_FITS
     integer            :: unit      = -1
     integer            :: cur_hdu   = 1
     character(len=256) :: filename  = ''
     logical            :: writable  = .false.
     type(hdf5_state_type), allocatable :: hdf5
  end type io_file_type

  public :: io_open_new, io_open_old, io_close, io_move_to_next_section
  public :: io_put_keyword, io_get_keyword, io_get_column_number
  public :: io_read_image, io_append_image, io_write_image
  public :: io_read_table_column, io_append_table_column, io_write_table_column
  public :: io_detect_format, io_extension, io_file_extension

  interface io_put_keyword
     module procedure io_put_key_logical, io_put_key_real32, io_put_key_real64, &
                      io_put_key_string,  io_put_key_int32,  io_put_key_int64
  end interface
  interface io_get_keyword
     module procedure io_get_key_logical, io_get_key_real32, io_get_key_real64, &
                      io_get_key_string,  io_get_key_int32,  io_get_key_int64
  end interface

  interface io_read_image
     module procedure io_read_1D_real32, io_read_2D_real32, io_read_3D_real32, io_read_4D_real32, &
                      io_read_1D_real64, io_read_2D_real64, io_read_3D_real64, io_read_4D_real64
  end interface
  interface io_append_image
     module procedure io_append_1D_real32, io_append_2D_real32, io_append_3D_real32, io_append_4D_real32, &
                      io_append_1D_real64, io_append_2D_real64, io_append_3D_real64, io_append_4D_real64, &
                      io_append_1D_int32,  io_append_2D_int32,  io_append_3D_int32,  io_append_4D_int32
  end interface
  interface io_write_image
     module procedure io_append_1D_real32, io_append_2D_real32, io_append_3D_real32, io_append_4D_real32, &
                      io_append_1D_real64, io_append_2D_real64, io_append_3D_real64, io_append_4D_real64, &
                      io_append_1D_int32,  io_append_2D_int32,  io_append_3D_int32,  io_append_4D_int32
  end interface

  interface io_read_table_column
     module procedure io_read_table_column_real32, io_read_table_column_real64, &
                      io_read_table_column_int32,  io_read_table_column_int64
  end interface
  interface io_append_table_column
     module procedure io_append_table_column_real32, io_append_table_column_real64, &
                      io_append_table_column_int32,  io_append_table_column_int64
  end interface
  interface io_write_table_column
     module procedure io_append_table_column_real32, io_append_table_column_real64, &
                      io_append_table_column_int32,  io_append_table_column_int64
  end interface

contains

!=============================================================================
! Format detection helpers
!=============================================================================
  function io_detect_format(fname) result(fmt)
    character(len=*), intent(in) :: fname
    integer :: fmt, ln
    ln = len_trim(fname)
    fmt = IO_FMT_FITS
    if (ln >= 3) then
       if (fname(ln-2:ln) == '.h5')    fmt = IO_FMT_HDF5
    endif
    if (ln >= 5) then
       if (fname(ln-4:ln) == '.hdf5')  fmt = IO_FMT_HDF5
    endif
  end function io_detect_format

  function io_extension(format) result(ext)
    integer, intent(in), optional :: format
    character(len=8) :: ext
    integer :: f
    f = IO_FMT_FITS
    if (present(format)) f = format
    if (f == IO_FMT_HDF5) then
       ext = '.h5'
    else
       ext = '.fits.gz'
    endif
  end function io_extension

  !---------------------------------------------------------------------------
  ! io_file_extension(par%file_format) → '.h5' or '.fits.gz'.
  ! Convenience wrapper used by write_output_*, sightline_tau_* etc. so they
  ! don't have to repeat the same case-insensitive comparison.
  !---------------------------------------------------------------------------
  function io_file_extension(format_str) result(ext)
    character(len=*), intent(in) :: format_str
    character(len=8) :: ext
    character(len=len(format_str)) :: s
    integer :: i
    s = adjustl(format_str)
    do i = 1, len_trim(s)
       if (s(i:i) >= 'A' .and. s(i:i) <= 'Z') s(i:i) = char(ichar(s(i:i)) + 32)
    enddo
    if (trim(s) == 'hdf5' .or. trim(s) == 'h5') then
       ext = '.h5'
    else
       ext = '.fits.gz'
    endif
  end function io_file_extension

!=============================================================================
! Open / close
!=============================================================================
  subroutine io_open_new(file, fname, status)
    type(io_file_type), intent(out)   :: file
    character(len=*),   intent(in)    :: fname
    integer,            intent(inout) :: status
    file%filename = fname
    file%writable = .true.
    file%cur_hdu  = 1
    file%format   = io_detect_format(fname)
    select case (file%format)
    case (IO_FMT_FITS)
       call fits_open_new(file%unit, fname, status)
    case (IO_FMT_HDF5)
       if (.not. allocated(file%hdf5)) allocate(file%hdf5)
       call hdf5_open_new(file%hdf5, fname, status)
    case default
       status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_open_new

  subroutine io_open_old(file, fname, status, overwrite)
    type(io_file_type), intent(out)   :: file
    character(len=*),   intent(in)    :: fname
    integer,            intent(inout) :: status
    logical, optional,  intent(in)    :: overwrite
    logical :: rw
    file%filename = fname
    file%cur_hdu  = 1
    file%format   = io_detect_format(fname)
    rw = .false.
    if (present(overwrite)) rw = overwrite
    file%writable = rw
    select case (file%format)
    case (IO_FMT_FITS)
       if (present(overwrite)) then
          call fits_open_old(file%unit, fname, status, overwrite)
       else
          call fits_open_old(file%unit, fname, status)
       endif
    case (IO_FMT_HDF5)
       if (.not. allocated(file%hdf5)) allocate(file%hdf5)
       call hdf5_open_old(file%hdf5, fname, status, rw)
    case default
       status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_open_old

  subroutine io_close(file, status)
    type(io_file_type), intent(inout) :: file
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS)
       call fits_close(file%unit, status)
    case (IO_FMT_HDF5)
       if (allocated(file%hdf5)) then
          call hdf5_close(file%hdf5, status)
          deallocate(file%hdf5)
       endif
    case default
       status = IO_ERR_UNKNOWN_FMT
    end select
    file%unit = -1
  end subroutine io_close

  subroutine io_move_to_next_section(file, status, nstep)
    type(io_file_type), intent(inout) :: file
    integer,            intent(inout) :: status
    integer, optional,  intent(in)    :: nstep
    integer :: step
    step = 1
    if (present(nstep)) then
       if (nstep >= 0) step = nstep
    endif
    select case (file%format)
    case (IO_FMT_FITS)
       call fits_move_to_next_hdu(file%unit, status, step)
       file%cur_hdu = file%cur_hdu + step
    case (IO_FMT_HDF5)
       call hdf5_move_to_next_section(file%hdf5, status, step)
    case default
       status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_move_to_next_section

  subroutine io_get_column_number(file, colname, colnum, status)
    type(io_file_type), intent(in)    :: file
    character(len=*),   intent(in)    :: colname
    integer,            intent(out)   :: colnum
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS)
       call fits_get_column_number(file%unit, colname, colnum, status)
    case (IO_FMT_HDF5)
       call hdf5_get_column_number(file%hdf5, colname, colnum, status)
    case default
       colnum = 0; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_get_column_number

!=============================================================================
! Keyword put
!=============================================================================
  subroutine io_put_key_logical(file, keyname, keyvalue, comment, status)
    type(io_file_type), intent(inout) :: file
    character(len=*),   intent(in)    :: keyname
    logical,            intent(in)    :: keyvalue
    character(len=*),   intent(in)    :: comment
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS); call fits_put_keyword(file%unit, keyname, keyvalue, comment, status)
    case (IO_FMT_HDF5); call hdf5_put_key_logical(file%hdf5, keyname, keyvalue, comment, status)
    case default;       status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_put_key_logical

  subroutine io_put_key_real32(file, keyname, keyvalue, comment, status)
    type(io_file_type), intent(inout) :: file
    character(len=*),   intent(in)    :: keyname
    real(real32),       intent(in)    :: keyvalue
    character(len=*),   intent(in)    :: comment
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS); call fits_put_keyword(file%unit, keyname, keyvalue, comment, status)
    case (IO_FMT_HDF5); call hdf5_put_key_real32(file%hdf5, keyname, keyvalue, comment, status)
    case default;       status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_put_key_real32

  subroutine io_put_key_real64(file, keyname, keyvalue, comment, status)
    type(io_file_type), intent(inout) :: file
    character(len=*),   intent(in)    :: keyname
    real(real64),       intent(in)    :: keyvalue
    character(len=*),   intent(in)    :: comment
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS); call fits_put_keyword(file%unit, keyname, keyvalue, comment, status)
    case (IO_FMT_HDF5); call hdf5_put_key_real64(file%hdf5, keyname, keyvalue, comment, status)
    case default;       status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_put_key_real64

  subroutine io_put_key_int32(file, keyname, keyvalue, comment, status)
    type(io_file_type), intent(inout) :: file
    character(len=*),   intent(in)    :: keyname
    integer(int32),     intent(in)    :: keyvalue
    character(len=*),   intent(in)    :: comment
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS); call fits_put_keyword(file%unit, keyname, keyvalue, comment, status)
    case (IO_FMT_HDF5); call hdf5_put_key_int32(file%hdf5, keyname, keyvalue, comment, status)
    case default;       status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_put_key_int32

  subroutine io_put_key_int64(file, keyname, keyvalue, comment, status)
    type(io_file_type), intent(inout) :: file
    character(len=*),   intent(in)    :: keyname
    integer(int64),     intent(in)    :: keyvalue
    character(len=*),   intent(in)    :: comment
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS); call fits_put_keyword(file%unit, keyname, keyvalue, comment, status)
    case (IO_FMT_HDF5); call hdf5_put_key_int64(file%hdf5, keyname, keyvalue, comment, status)
    case default;       status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_put_key_int64

  subroutine io_put_key_string(file, keyname, keyvalue, comment, status)
    type(io_file_type), intent(inout) :: file
    character(len=*),   intent(in)    :: keyname
    character(len=*),   intent(in)    :: keyvalue
    character(len=*),   intent(in)    :: comment
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS); call fits_put_keyword(file%unit, keyname, keyvalue, comment, status)
    case (IO_FMT_HDF5); call hdf5_put_key_string(file%hdf5, keyname, keyvalue, comment, status)
    case default;       status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_put_key_string

!=============================================================================
! Keyword get
!=============================================================================
  subroutine io_get_key_logical(file, keyname, keyvalue, status, comment)
    type(io_file_type),         intent(in)    :: file
    character(len=*),           intent(in)    :: keyname
    logical,                    intent(out)   :: keyvalue
    integer,                    intent(inout) :: status
    character(len=*), optional, intent(out)   :: comment
    select case (file%format)
    case (IO_FMT_FITS)
       if (present(comment)) then
          call fits_get_keyword(file%unit, keyname, keyvalue, status, comment)
       else
          call fits_get_keyword(file%unit, keyname, keyvalue, status)
       endif
    case (IO_FMT_HDF5)
       if (present(comment)) then
          call hdf5_get_key_logical(file%hdf5, keyname, keyvalue, status, comment)
       else
          call hdf5_get_key_logical(file%hdf5, keyname, keyvalue, status)
       endif
    case default; keyvalue = .false.; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_get_key_logical

  subroutine io_get_key_real32(file, keyname, keyvalue, status, comment)
    type(io_file_type),         intent(in)    :: file
    character(len=*),           intent(in)    :: keyname
    real(real32),               intent(out)   :: keyvalue
    integer,                    intent(inout) :: status
    character(len=*), optional, intent(out)   :: comment
    select case (file%format)
    case (IO_FMT_FITS)
       if (present(comment)) then
          call fits_get_keyword(file%unit, keyname, keyvalue, status, comment)
       else
          call fits_get_keyword(file%unit, keyname, keyvalue, status)
       endif
    case (IO_FMT_HDF5)
       if (present(comment)) then
          call hdf5_get_key_real32(file%hdf5, keyname, keyvalue, status, comment)
       else
          call hdf5_get_key_real32(file%hdf5, keyname, keyvalue, status)
       endif
    case default; keyvalue = 0.0_real32; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_get_key_real32

  subroutine io_get_key_real64(file, keyname, keyvalue, status, comment)
    type(io_file_type),         intent(in)    :: file
    character(len=*),           intent(in)    :: keyname
    real(real64),               intent(out)   :: keyvalue
    integer,                    intent(inout) :: status
    character(len=*), optional, intent(out)   :: comment
    select case (file%format)
    case (IO_FMT_FITS)
       if (present(comment)) then
          call fits_get_keyword(file%unit, keyname, keyvalue, status, comment)
       else
          call fits_get_keyword(file%unit, keyname, keyvalue, status)
       endif
    case (IO_FMT_HDF5)
       if (present(comment)) then
          call hdf5_get_key_real64(file%hdf5, keyname, keyvalue, status, comment)
       else
          call hdf5_get_key_real64(file%hdf5, keyname, keyvalue, status)
       endif
    case default; keyvalue = 0.0_real64; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_get_key_real64

  subroutine io_get_key_int32(file, keyname, keyvalue, status, comment)
    type(io_file_type),         intent(in)    :: file
    character(len=*),           intent(in)    :: keyname
    integer(int32),             intent(out)   :: keyvalue
    integer,                    intent(inout) :: status
    character(len=*), optional, intent(out)   :: comment
    select case (file%format)
    case (IO_FMT_FITS)
       if (present(comment)) then
          call fits_get_keyword(file%unit, keyname, keyvalue, status, comment)
       else
          call fits_get_keyword(file%unit, keyname, keyvalue, status)
       endif
    case (IO_FMT_HDF5)
       if (present(comment)) then
          call hdf5_get_key_int32(file%hdf5, keyname, keyvalue, status, comment)
       else
          call hdf5_get_key_int32(file%hdf5, keyname, keyvalue, status)
       endif
    case default; keyvalue = 0_int32; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_get_key_int32

  subroutine io_get_key_int64(file, keyname, keyvalue, status, comment)
    type(io_file_type),         intent(in)    :: file
    character(len=*),           intent(in)    :: keyname
    integer(int64),             intent(out)   :: keyvalue
    integer,                    intent(inout) :: status
    character(len=*), optional, intent(out)   :: comment
    select case (file%format)
    case (IO_FMT_FITS)
       if (present(comment)) then
          call fits_get_keyword(file%unit, keyname, keyvalue, status, comment)
       else
          call fits_get_keyword(file%unit, keyname, keyvalue, status)
       endif
    case (IO_FMT_HDF5)
       if (present(comment)) then
          call hdf5_get_key_int64(file%hdf5, keyname, keyvalue, status, comment)
       else
          call hdf5_get_key_int64(file%hdf5, keyname, keyvalue, status)
       endif
    case default; keyvalue = 0_int64; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_get_key_int64

  subroutine io_get_key_string(file, keyname, keyvalue, status, comment)
    type(io_file_type),         intent(in)    :: file
    character(len=*),           intent(in)    :: keyname
    character(len=*),           intent(out)   :: keyvalue
    integer,                    intent(inout) :: status
    character(len=*), optional, intent(out)   :: comment
    select case (file%format)
    case (IO_FMT_FITS)
       if (present(comment)) then
          call fits_get_keyword(file%unit, keyname, keyvalue, status, comment)
       else
          call fits_get_keyword(file%unit, keyname, keyvalue, status)
       endif
    case (IO_FMT_HDF5)
       if (present(comment)) then
          call hdf5_get_key_string(file%hdf5, keyname, keyvalue, status, comment)
       else
          call hdf5_get_key_string(file%hdf5, keyname, keyvalue, status)
       endif
    case default; keyvalue = ''; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_get_key_string

!=============================================================================
! Image read (real32 / real64; 1D-4D)
!=============================================================================
  subroutine io_read_1D_real32(file, array, status)
    type(io_file_type), intent(in)    :: file
    real(real32),       intent(out)   :: array(:)
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS); call fits_read_image(file%unit, array, status)
    case (IO_FMT_HDF5); call hdf5_read_1D_real32(file%hdf5, array, status)
    case default;       array = 0.0_real32; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_read_1D_real32

  subroutine io_read_2D_real32(file, array, status)
    type(io_file_type), intent(in)    :: file
    real(real32),       intent(out)   :: array(:,:)
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS); call fits_read_image(file%unit, array, status)
    case (IO_FMT_HDF5); call hdf5_read_2D_real32(file%hdf5, array, status)
    case default;       array = 0.0_real32; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_read_2D_real32

  subroutine io_read_3D_real32(file, array, status)
    type(io_file_type), intent(in)    :: file
    real(real32),       intent(out)   :: array(:,:,:)
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS); call fits_read_image(file%unit, array, status)
    case (IO_FMT_HDF5); call hdf5_read_3D_real32(file%hdf5, array, status)
    case default;       array = 0.0_real32; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_read_3D_real32

  subroutine io_read_4D_real32(file, array, status)
    type(io_file_type), intent(in)    :: file
    real(real32),       intent(out)   :: array(:,:,:,:)
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS); call fits_read_image(file%unit, array, status)
    case (IO_FMT_HDF5); call hdf5_read_4D_real32(file%hdf5, array, status)
    case default;       array = 0.0_real32; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_read_4D_real32

  subroutine io_read_1D_real64(file, array, status)
    type(io_file_type), intent(in)    :: file
    real(real64),       intent(out)   :: array(:)
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS); call fits_read_image(file%unit, array, status)
    case (IO_FMT_HDF5); call hdf5_read_1D_real64(file%hdf5, array, status)
    case default;       array = 0.0_real64; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_read_1D_real64

  subroutine io_read_2D_real64(file, array, status)
    type(io_file_type), intent(in)    :: file
    real(real64),       intent(out)   :: array(:,:)
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS); call fits_read_image(file%unit, array, status)
    case (IO_FMT_HDF5); call hdf5_read_2D_real64(file%hdf5, array, status)
    case default;       array = 0.0_real64; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_read_2D_real64

  subroutine io_read_3D_real64(file, array, status)
    type(io_file_type), intent(in)    :: file
    real(real64),       intent(out)   :: array(:,:,:)
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS); call fits_read_image(file%unit, array, status)
    case (IO_FMT_HDF5); call hdf5_read_3D_real64(file%hdf5, array, status)
    case default;       array = 0.0_real64; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_read_3D_real64

  subroutine io_read_4D_real64(file, array, status)
    type(io_file_type), intent(in)    :: file
    real(real64),       intent(out)   :: array(:,:,:,:)
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS); call fits_read_image(file%unit, array, status)
    case (IO_FMT_HDF5); call hdf5_read_4D_real64(file%hdf5, array, status)
    case default;       array = 0.0_real64; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_read_4D_real64

!=============================================================================
! Image append
!=============================================================================
  subroutine io_append_1D_real32(file, array, status, bitpix)
    type(io_file_type), intent(inout) :: file
    real(real32),       intent(in)    :: array(:)
    integer,            intent(inout) :: status
    integer, optional,  intent(in)    :: bitpix
    select case (file%format)
    case (IO_FMT_FITS)
       if (present(bitpix)) then; call fits_append_image(file%unit, array, status, bitpix)
       else;                       call fits_append_image(file%unit, array, status); endif
    case (IO_FMT_HDF5)
       if (present(bitpix)) then; call hdf5_append_1D_real32(file%hdf5, array, status, bitpix)
       else;                       call hdf5_append_1D_real32(file%hdf5, array, status); endif
    case default; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_append_1D_real32

  subroutine io_append_2D_real32(file, array, status, bitpix)
    type(io_file_type), intent(inout) :: file
    real(real32),       intent(in)    :: array(:,:)
    integer,            intent(inout) :: status
    integer, optional,  intent(in)    :: bitpix
    select case (file%format)
    case (IO_FMT_FITS)
       if (present(bitpix)) then; call fits_append_image(file%unit, array, status, bitpix)
       else;                       call fits_append_image(file%unit, array, status); endif
    case (IO_FMT_HDF5)
       if (present(bitpix)) then; call hdf5_append_2D_real32(file%hdf5, array, status, bitpix)
       else;                       call hdf5_append_2D_real32(file%hdf5, array, status); endif
    case default; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_append_2D_real32

  subroutine io_append_3D_real32(file, array, status, bitpix)
    type(io_file_type), intent(inout) :: file
    real(real32),       intent(in)    :: array(:,:,:)
    integer,            intent(inout) :: status
    integer, optional,  intent(in)    :: bitpix
    select case (file%format)
    case (IO_FMT_FITS)
       if (present(bitpix)) then; call fits_append_image(file%unit, array, status, bitpix)
       else;                       call fits_append_image(file%unit, array, status); endif
    case (IO_FMT_HDF5)
       if (present(bitpix)) then; call hdf5_append_3D_real32(file%hdf5, array, status, bitpix)
       else;                       call hdf5_append_3D_real32(file%hdf5, array, status); endif
    case default; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_append_3D_real32

  subroutine io_append_4D_real32(file, array, status, bitpix)
    type(io_file_type), intent(inout) :: file
    real(real32),       intent(in)    :: array(:,:,:,:)
    integer,            intent(inout) :: status
    integer, optional,  intent(in)    :: bitpix
    select case (file%format)
    case (IO_FMT_FITS)
       if (present(bitpix)) then; call fits_append_image(file%unit, array, status, bitpix)
       else;                       call fits_append_image(file%unit, array, status); endif
    case (IO_FMT_HDF5)
       if (present(bitpix)) then; call hdf5_append_4D_real32(file%hdf5, array, status, bitpix)
       else;                       call hdf5_append_4D_real32(file%hdf5, array, status); endif
    case default; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_append_4D_real32

  subroutine io_append_1D_real64(file, array, status, bitpix)
    type(io_file_type), intent(inout) :: file
    real(real64),       intent(in)    :: array(:)
    integer,            intent(inout) :: status
    integer, optional,  intent(in)    :: bitpix
    select case (file%format)
    case (IO_FMT_FITS)
       if (present(bitpix)) then; call fits_append_image(file%unit, array, status, bitpix)
       else;                       call fits_append_image(file%unit, array, status); endif
    case (IO_FMT_HDF5)
       if (present(bitpix)) then; call hdf5_append_1D_real64(file%hdf5, array, status, bitpix)
       else;                       call hdf5_append_1D_real64(file%hdf5, array, status); endif
    case default; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_append_1D_real64

  subroutine io_append_2D_real64(file, array, status, bitpix)
    type(io_file_type), intent(inout) :: file
    real(real64),       intent(in)    :: array(:,:)
    integer,            intent(inout) :: status
    integer, optional,  intent(in)    :: bitpix
    select case (file%format)
    case (IO_FMT_FITS)
       if (present(bitpix)) then; call fits_append_image(file%unit, array, status, bitpix)
       else;                       call fits_append_image(file%unit, array, status); endif
    case (IO_FMT_HDF5)
       if (present(bitpix)) then; call hdf5_append_2D_real64(file%hdf5, array, status, bitpix)
       else;                       call hdf5_append_2D_real64(file%hdf5, array, status); endif
    case default; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_append_2D_real64

  subroutine io_append_3D_real64(file, array, status, bitpix)
    type(io_file_type), intent(inout) :: file
    real(real64),       intent(in)    :: array(:,:,:)
    integer,            intent(inout) :: status
    integer, optional,  intent(in)    :: bitpix
    select case (file%format)
    case (IO_FMT_FITS)
       if (present(bitpix)) then; call fits_append_image(file%unit, array, status, bitpix)
       else;                       call fits_append_image(file%unit, array, status); endif
    case (IO_FMT_HDF5)
       if (present(bitpix)) then; call hdf5_append_3D_real64(file%hdf5, array, status, bitpix)
       else;                       call hdf5_append_3D_real64(file%hdf5, array, status); endif
    case default; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_append_3D_real64

  subroutine io_append_4D_real64(file, array, status, bitpix)
    type(io_file_type), intent(inout) :: file
    real(real64),       intent(in)    :: array(:,:,:,:)
    integer,            intent(inout) :: status
    integer, optional,  intent(in)    :: bitpix
    select case (file%format)
    case (IO_FMT_FITS)
       if (present(bitpix)) then; call fits_append_image(file%unit, array, status, bitpix)
       else;                       call fits_append_image(file%unit, array, status); endif
    case (IO_FMT_HDF5)
       if (present(bitpix)) then; call hdf5_append_4D_real64(file%hdf5, array, status, bitpix)
       else;                       call hdf5_append_4D_real64(file%hdf5, array, status); endif
    case default; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_append_4D_real64

  subroutine io_append_1D_int32(file, array, status)
    type(io_file_type), intent(inout) :: file
    integer,            intent(in)    :: array(:)
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS); call fits_append_image(file%unit, array, status)
    case (IO_FMT_HDF5); call hdf5_append_1D_int32(file%hdf5, array, status)
    case default;       status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_append_1D_int32

  subroutine io_append_2D_int32(file, array, status)
    type(io_file_type), intent(inout) :: file
    integer,            intent(in)    :: array(:,:)
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS); call fits_append_image(file%unit, array, status)
    case (IO_FMT_HDF5); call hdf5_append_2D_int32(file%hdf5, array, status)
    case default;       status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_append_2D_int32

  subroutine io_append_3D_int32(file, array, status)
    type(io_file_type), intent(inout) :: file
    integer,            intent(in)    :: array(:,:,:)
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS); call fits_append_image(file%unit, array, status)
    case (IO_FMT_HDF5); call hdf5_append_3D_int32(file%hdf5, array, status)
    case default;       status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_append_3D_int32

  subroutine io_append_4D_int32(file, array, status)
    type(io_file_type), intent(inout) :: file
    integer,            intent(in)    :: array(:,:,:,:)
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS); call fits_append_image(file%unit, array, status)
    case (IO_FMT_HDF5); call hdf5_append_4D_int32(file%hdf5, array, status)
    case default;       status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_append_4D_int32

!=============================================================================
! Bin-table column read
!=============================================================================
  subroutine io_read_table_column_real32(file, colnum, array, status)
    type(io_file_type), intent(in)    :: file
    integer,            intent(in)    :: colnum
    real(real32),       intent(out)   :: array(:)
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS); call fits_read_table_column(file%unit, colnum, array, status)
    case (IO_FMT_HDF5); call hdf5_read_table_column_real32(file%hdf5, colnum, array, status)
    case default;       array = 0.0_real32; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_read_table_column_real32

  subroutine io_read_table_column_real64(file, colnum, array, status)
    type(io_file_type), intent(in)    :: file
    integer,            intent(in)    :: colnum
    real(real64),       intent(out)   :: array(:)
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS); call fits_read_table_column(file%unit, colnum, array, status)
    case (IO_FMT_HDF5); call hdf5_read_table_column_real64(file%hdf5, colnum, array, status)
    case default;       array = 0.0_real64; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_read_table_column_real64

  subroutine io_read_table_column_int32(file, colnum, array, status)
    type(io_file_type), intent(in)    :: file
    integer,            intent(in)    :: colnum
    integer(int32),     intent(out)   :: array(:)
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS); call fits_read_table_column(file%unit, colnum, array, status)
    case (IO_FMT_HDF5); call hdf5_read_table_column_int32(file%hdf5, colnum, array, status)
    case default;       array = 0_int32; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_read_table_column_int32

  subroutine io_read_table_column_int64(file, colnum, array, status)
    type(io_file_type), intent(in)    :: file
    integer,            intent(in)    :: colnum
    integer(int64),     intent(out)   :: array(:)
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS); call fits_read_table_column(file%unit, colnum, array, status)
    case (IO_FMT_HDF5); call hdf5_read_table_column_int64(file%hdf5, colnum, array, status)
    case default;       array = 0_int64; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_read_table_column_int64

!=============================================================================
! Bin-table column append
!=============================================================================
  subroutine io_append_table_column_real32(file, ttype, array, status, bitpix)
    type(io_file_type), intent(inout) :: file
    character(len=*),   intent(in)    :: ttype
    real(real32),       intent(in)    :: array(:)
    integer,            intent(inout) :: status
    integer, optional,  intent(in)    :: bitpix
    select case (file%format)
    case (IO_FMT_FITS)
       if (present(bitpix)) then; call fits_append_table_column(file%unit, ttype, array, status, bitpix)
       else;                       call fits_append_table_column(file%unit, ttype, array, status); endif
    case (IO_FMT_HDF5)
       if (present(bitpix)) then; call hdf5_append_table_column_real32(file%hdf5, ttype, array, status, bitpix)
       else;                       call hdf5_append_table_column_real32(file%hdf5, ttype, array, status); endif
    case default; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_append_table_column_real32

  subroutine io_append_table_column_real64(file, ttype, array, status, bitpix)
    type(io_file_type), intent(inout) :: file
    character(len=*),   intent(in)    :: ttype
    real(real64),       intent(in)    :: array(:)
    integer,            intent(inout) :: status
    integer, optional,  intent(in)    :: bitpix
    select case (file%format)
    case (IO_FMT_FITS)
       if (present(bitpix)) then; call fits_append_table_column(file%unit, ttype, array, status, bitpix)
       else;                       call fits_append_table_column(file%unit, ttype, array, status); endif
    case (IO_FMT_HDF5)
       if (present(bitpix)) then; call hdf5_append_table_column_real64(file%hdf5, ttype, array, status, bitpix)
       else;                       call hdf5_append_table_column_real64(file%hdf5, ttype, array, status); endif
    case default; status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_append_table_column_real64

  subroutine io_append_table_column_int32(file, ttype, array, status)
    type(io_file_type), intent(inout) :: file
    character(len=*),   intent(in)    :: ttype
    integer(int32),     intent(in)    :: array(:)
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS); call fits_append_table_column(file%unit, ttype, array, status)
    case (IO_FMT_HDF5); call hdf5_append_table_column_int32(file%hdf5, ttype, array, status)
    case default;       status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_append_table_column_int32

  subroutine io_append_table_column_int64(file, ttype, array, status)
    type(io_file_type), intent(inout) :: file
    character(len=*),   intent(in)    :: ttype
    integer(int64),     intent(in)    :: array(:)
    integer,            intent(inout) :: status
    select case (file%format)
    case (IO_FMT_FITS); call fits_append_table_column(file%unit, ttype, array, status)
    case (IO_FMT_HDF5); call hdf5_append_table_column_int64(file%hdf5, ttype, array, status)
    case default;       status = IO_ERR_UNKNOWN_FMT
    end select
  end subroutine io_append_table_column_int64

end module iofile_mod
