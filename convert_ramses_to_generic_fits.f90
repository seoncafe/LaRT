program convert_ramses_to_generic_fits
  use define, only: wp, kpc2cm, pc2cm, au2cm
  use read_ramses_amr_mod, only: ramses_read_leaf_cells
  use fitsio_mod
  use, intrinsic :: iso_fortran_env, only: int32
  implicit none

  character(len=512) :: repository, output_file
  character(len=32)  :: output_unit
  character(len=64)  :: arg
  real(wp), allocatable :: xleaf(:), yleaf(:), zleaf(:)
  real(wp), allocatable :: nH_cgs(:), T_cgs(:), vx_kms(:), vy_kms(:), vz_kms(:)
  integer,  allocatable :: leaf_level(:)
  integer(int32), allocatable :: leaf_level_i4(:)
  integer :: snapnum, nleaf
  real(wp) :: boxlen_cm, unit2cm, unit_l_cgs
  integer :: status, unit

  call parse_args(repository, snapnum, output_file, output_unit)
  call read_ramses_unit_l(trim(repository), snapnum, unit_l_cgs)

  call ramses_read_leaf_cells(trim(repository), snapnum, &
      xleaf, yleaf, zleaf, leaf_level, &
      nH_cgs, T_cgs, vx_kms, vy_kms, vz_kms, &
      nleaf, boxlen_cm)

  call get_unit_scale(trim(output_unit), unit2cm)

  xleaf = xleaf / unit2cm
  yleaf = yleaf / unit2cm
  zleaf = zleaf / unit2cm
  boxlen_cm = boxlen_cm / unit2cm

  allocate(leaf_level_i4(nleaf))
  leaf_level_i4 = int(leaf_level, int32)

  status = 0
  call fits_open_new(unit, trim(output_file), status)
  if (status /= 0) stop 'convert_ramses_to_generic_fits: cannot create output FITS file'

  call fits_write_table_column(unit, 'x',      xleaf,         status)
  call fits_write_table_column(unit, 'y',      yleaf,         status)
  call fits_write_table_column(unit, 'z',      zleaf,         status)
  call fits_write_table_column(unit, 'level',  leaf_level_i4, status)
  call fits_write_table_column(unit, 'gasDen', nH_cgs,        status)
  call fits_write_table_column(unit, 'T',      T_cgs,         status)
  call fits_write_table_column(unit, 'vx',     vx_kms,        status)
  call fits_write_table_column(unit, 'vy',     vy_kms,        status)
  call fits_write_table_column(unit, 'vz',     vz_kms,        status)

  call fits_put_keyword(unit, 'BOXLEN',  boxlen_cm,               'Simulation box length', status)
  call fits_put_keyword(unit, 'ORIGINX', 0.0_wp,                  'Box origin x',          status)
  call fits_put_keyword(unit, 'ORIGINY', 0.0_wp,                  'Box origin y',          status)
  call fits_put_keyword(unit, 'ORIGINZ', 0.0_wp,                  'Box origin z',          status)
  call fits_put_keyword(unit, 'NLEAF',   int(nleaf, int32),       'Number of AMR leaf cells', status)
  call fits_put_keyword(unit, 'UNITLCGS', unit_l_cgs,             'RAMSES unit_l in cm',   status)
  call fits_put_keyword(unit, 'UNITPOS', trim(output_unit),       'Position unit in table', status)

  call fits_close(unit, status)
  if (status /= 0) stop 'convert_ramses_to_generic_fits: error while closing FITS file'

  write(*,'(a)') 'Converted RAMSES snapshot to generic FITS'
  write(*,'(a)') '  repository : '//trim(repository)
  write(*,'(a,i0)') '  snapnum    : ', snapnum
  write(*,'(a)') '  output     : '//trim(output_file)
  write(*,'(a,i0)') '  nleaf      : ', nleaf
  write(*,'(a,es12.4,1x,a)') '  boxlen     : ', boxlen_cm, trim(output_unit)

contains

  subroutine parse_args(repository, snapnum, output_file, output_unit)
    character(len=*), intent(out) :: repository, output_file, output_unit
    integer,          intent(out) :: snapnum
    integer :: nargs

    nargs = command_argument_count()
    if (nargs < 3 .or. nargs > 4) then
      call print_usage()
      stop 1
    end if

    call get_command_argument(1, repository)
    call get_command_argument(2, arg)
    read(arg, *, err=100) snapnum
    call get_command_argument(3, output_file)
    output_unit = 'kpc'
    if (nargs >= 4) call get_command_argument(4, output_unit)
    output_unit = trim(str_lower(output_unit))
    if (output_unit /= 'cm' .and. output_unit /= 'kpc' .and. output_unit /= 'pc' .and. output_unit /= 'au') then
      write(*,'(a)') 'Unsupported output unit: '//trim(output_unit)
      call print_usage()
      stop 1
    end if
    return
100 continue
    write(*,'(a)') 'Invalid snapshot number: '//trim(arg)
    call print_usage()
    stop 1
  end subroutine parse_args

  subroutine print_usage()
    write(*,'(a)') 'Usage: convert_ramses_to_generic_fits.x <repository> <snapnum> <output.fits|output.fits.gz> [kpc|pc|au|cm]'
  end subroutine print_usage

  subroutine get_unit_scale(unit_name, unit2cm)
    character(len=*), intent(in)  :: unit_name
    real(wp),         intent(out) :: unit2cm

    select case (trim(unit_name))
    case ('cm')
      unit2cm = 1.0_wp
    case ('kpc')
      unit2cm = kpc2cm
    case ('pc')
      unit2cm = pc2cm
    case ('au')
      unit2cm = au2cm
    case default
      stop 'convert_ramses_to_generic_fits: unsupported unit'
    end select
  end subroutine get_unit_scale

  subroutine read_ramses_unit_l(repository, snapnum, unit_l_cgs)
    character(len=*), intent(in)  :: repository
    integer,          intent(in)  :: snapnum
    real(wp),         intent(out) :: unit_l_cgs

    character(len=512) :: filename, line, key
    integer :: ios, info_unit, eqpos

    write(filename, '(a,"/output_",i5.5,"/info_",i5.5,".txt")') &
        trim(repository), snapnum, snapnum

    unit_l_cgs = 1.0_wp
    info_unit = 77
    open(info_unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) return

    do
      read(info_unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      eqpos = index(line, '=')
      if (eqpos <= 0) cycle
      key = adjustl(trim(line(:eqpos-1)))
      if (trim(key) /= 'unit_l') cycle
      read(line(eqpos+1:), *, iostat=ios) unit_l_cgs
      exit
    end do

    close(info_unit)
  end subroutine read_ramses_unit_l

  pure function str_lower(str) result(out)
    character(len=*), intent(in) :: str
    character(len=len(str)) :: out
    integer :: i, code

    out = str
    do i = 1, len(str)
      code = iachar(out(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) out(i:i) = achar(code + 32)
    end do
  end function str_lower

end program convert_ramses_to_generic_fits
