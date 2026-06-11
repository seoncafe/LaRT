program convert_ramses_to_generic
  use define, only: wp, kpc2cm, pc2cm, au2cm
  use read_ramses_amr_mod, only: ramses_read_leaf_cells, &
      ramses_metal_var, ramses_xHII_var, ramses_xHeII_var, ramses_xHeIII_var
  use physics_amr_mod
  use iofile_mod
  use, intrinsic :: iso_fortran_env, only: int32
  implicit none

  character(len=512) :: repository, output_file
  character(len=32)  :: output_unit
  character(len=64)  :: arg
  real(wp), allocatable :: xleaf(:), yleaf(:), zleaf(:)
  real(wp), allocatable :: nH_cgs(:), T_cgs(:), vx_kms(:), vy_kms(:), vz_kms(:)
  integer,  allocatable :: leaf_level(:)
  integer(int32), allocatable :: leaf_level_i4(:)
  ! Extended physics arrays
  real(wp), allocatable :: xHI_arr(:), ne_arr(:), ndust_arr(:), emiss_arr(:)
  ! RAMSES hydro extra fields
  real(wp), allocatable :: Z_ramses(:), xHII_ramses(:), xHeII_ramses(:), xHeIII_ramses(:)
  integer :: snapnum, nleaf
  real(wp) :: boxlen_cm, boxlen_eff, unit2cm, unit_l_cgs
  integer :: status, il
  type(io_file_type) :: iofh
  logical :: compute_physics
  real(wp) :: Z_global
  integer  :: nx_base
  integer  :: base_ix, base_iy, base_iz
  integer  :: ix_min, iy_min, iz_min, ix_max, iy_max, iz_max
  integer  :: sub_ix, sub_iy, sub_iz, extent_cells, m_sub, level_bump
  real(wp) :: base_cell_cm

  call parse_args(repository, snapnum, output_file, output_unit, compute_physics, Z_global)
  call read_ramses_unit_l(trim(repository), snapnum, unit_l_cgs)

  call ramses_read_leaf_cells(trim(repository), snapnum, &
      xleaf, yleaf, zleaf, leaf_level, &
      nH_cgs, T_cgs, vx_kms, vy_kms, vz_kms, &
      nleaf, boxlen_cm, &
      metallicity=Z_ramses, xHII=xHII_ramses, &
      xHeII=xHeII_ramses, xHeIII=xHeIII_ramses, &
      nx_base=nx_base)

  call get_unit_scale(trim(output_unit), unit2cm)

  xleaf = xleaf / unit2cm
  yleaf = yleaf / unit2cm
  zleaf = zleaf / unit2cm
  boxlen_cm = boxlen_cm / unit2cm

  ! Compute extended physics if requested
  if (compute_physics) then
    allocate(xHI_arr(nleaf), ne_arr(nleaf), ndust_arr(nleaf), emiss_arr(nleaf))
    do il = 1, nleaf
      xHI_arr(il)   = cie_neutral_fraction_formula(T_cgs(il))
      ne_arr(il)    = electron_density_from_xHI(nH_cgs(il), xHI_arr(il))
      emiss_arr(il) = caseB_lya_emissivity(nH_cgs(il), T_cgs(il), xHI_arr(il), ne_arr(il))
      if (Z_global >= 0.0_wp) then
        ndust_arr(il) = laursen09_ndust(nH_cgs(il), xHI_arr(il), Z_global, 0.0134_wp, 0.01_wp)
      else
        ndust_arr(il) = 0.0_wp
      end if
    end do
  end if

  allocate(leaf_level_i4(nleaf))
  leaf_level_i4 = int(leaf_level, int32)

  ! Octree-alignment fix for RAMSES nx_base > 1.  See the Python converter
  ! for the full rationale.  Briefly: native cell size at level L is
  ! boxlen/(nx_base * 2^L), which cannot match a LaRT octree of step
  ! boxlen/2^L' unless nx_base is a power of 2.  We re-anchor to the
  ! smallest power-of-2 cubic sub-block of the nx_base^3 base grid that
  ! encloses every populated leaf and bump levels by log2(m_sub) so the
  ! native cell size remains exact on the shrunk box.
  if (nx_base > 1) then
    base_cell_cm = boxlen_cm / real(nx_base, wp)
    ix_min = nx_base; iy_min = nx_base; iz_min = nx_base
    ix_max = -1;     iy_max = -1;     iz_max = -1
    do il = 1, nleaf
      base_ix = int(floor(xleaf(il) / base_cell_cm))
      base_iy = int(floor(yleaf(il) / base_cell_cm))
      base_iz = int(floor(zleaf(il) / base_cell_cm))
      if (base_ix < 0) base_ix = 0
      if (base_iy < 0) base_iy = 0
      if (base_iz < 0) base_iz = 0
      if (base_ix > nx_base - 1) base_ix = nx_base - 1
      if (base_iy > nx_base - 1) base_iy = nx_base - 1
      if (base_iz > nx_base - 1) base_iz = nx_base - 1
      if (base_ix < ix_min) ix_min = base_ix
      if (base_iy < iy_min) iy_min = base_iy
      if (base_iz < iz_min) iz_min = base_iz
      if (base_ix > ix_max) ix_max = base_ix
      if (base_iy > iy_max) iy_max = base_iy
      if (base_iz > iz_max) iz_max = base_iz
    end do
    extent_cells = max(ix_max - ix_min + 1, iy_max - iy_min + 1, iz_max - iz_min + 1)
    m_sub = 1
    do while (m_sub < extent_cells)
      m_sub = m_sub * 2
    end do
    if (m_sub > nx_base) then
      write(*,'(a,i0,a,i0,a)') 'ERROR: populated extent ', extent_cells, &
        ' base cells exceeds nx_base=', nx_base, &
        '; cannot fit on a power-of-2 sub-block.'
      stop 1
    end if
    sub_ix = min(ix_min, nx_base - m_sub)
    sub_iy = min(iy_min, nx_base - m_sub)
    sub_iz = min(iz_min, nx_base - m_sub)
    sub_ix = min(sub_ix, ix_max - m_sub + 1)
    sub_iy = min(sub_iy, iy_max - m_sub + 1)
    sub_iz = min(sub_iz, iz_max - m_sub + 1)
    sub_ix = max(0, sub_ix); sub_iy = max(0, sub_iy); sub_iz = max(0, sub_iz)
    do il = 1, nleaf
      base_ix = int(floor(xleaf(il) / base_cell_cm))
      base_iy = int(floor(yleaf(il) / base_cell_cm))
      base_iz = int(floor(zleaf(il) / base_cell_cm))
      if (base_ix < sub_ix .or. base_ix >= sub_ix + m_sub .or. &
          base_iy < sub_iy .or. base_iy >= sub_iy + m_sub .or. &
          base_iz < sub_iz .or. base_iz >= sub_iz + m_sub) then
        write(*,'(a)') 'ERROR: cell outside chosen sub-block (internal error).'
        stop 1
      end if
    end do
    xleaf(:) = xleaf(:) - real(sub_ix, wp) * base_cell_cm
    yleaf(:) = yleaf(:) - real(sub_iy, wp) * base_cell_cm
    zleaf(:) = zleaf(:) - real(sub_iz, wp) * base_cell_cm
    boxlen_eff = real(m_sub, wp) * base_cell_cm
    level_bump = 0
    do while (ishft(1, level_bump) < m_sub)
      level_bump = level_bump + 1
    end do
    if (level_bump > 0) then
      leaf_level_i4(:) = leaf_level_i4(:) + int(level_bump, int32)
    end if
    write(*,'(a,i0,a,i0,a,es12.4,a,es12.4,1x,a)') &
        '  RAMSES nx_base=', nx_base, &
        ' > 1: shrinking BOXLEN to a ', m_sub, '-base-cell-wide block (', &
        boxlen_cm, ' -> ', boxlen_eff, trim(output_unit)
    write(*,'(a,6i4)') '  populated bbox (base cells): ', &
        ix_min, ix_max, iy_min, iy_max, iz_min, iz_max
    write(*,'(a,3i4,a,i0)') '  sub-block origin: ', sub_ix, sub_iy, sub_iz, &
        '; level shift: +', level_bump
  else
    boxlen_eff = boxlen_cm
  end if

  ! Shift cell positions from RAMSES corner-based [0, boxlen_eff] to LaRT's
  ! centered convention [-boxlen_eff/2, +boxlen_eff/2].
  xleaf(:) = xleaf(:) - 0.5_wp * boxlen_eff
  yleaf(:) = yleaf(:) - 0.5_wp * boxlen_eff
  zleaf(:) = zleaf(:) - 0.5_wp * boxlen_eff

  status = 0
  call io_open_new(iofh, trim(output_file), status)
  if (status /= 0) stop 'convert_ramses_to_generic: cannot create output file'

  ! Mandatory 9 columns
  call io_write_table_column(iofh, 'x',      xleaf,         status)
  call io_write_table_column(iofh, 'y',      yleaf,         status)
  call io_write_table_column(iofh, 'z',      zleaf,         status)
  call io_write_table_column(iofh, 'level',  leaf_level_i4, status)
  call io_write_table_column(iofh, 'gasDen', nH_cgs,        status)
  call io_write_table_column(iofh, 'T',      T_cgs,         status)
  call io_write_table_column(iofh, 'vx',     vx_kms,        status)
  call io_write_table_column(iofh, 'vy',     vy_kms,        status)
  call io_write_table_column(iofh, 'vz',     vz_kms,        status)

  ! Extended columns (when physics computed)
  if (compute_physics) then
    call io_write_table_column(iofh, 'xHI',        xHI_arr,   status)
    call io_write_table_column(iofh, 'n_e',        ne_arr,    status)
    call io_write_table_column(iofh, 'emissivity', emiss_arr, status)
    if (Z_global >= 0.0_wp) then
      call io_write_table_column(iofh, 'ndust', ndust_arr, status)
    end if
  end if

  ! RAMSES hydro extra columns (when detected in descriptor)
  if (allocated(Z_ramses)) then
    call io_write_table_column(iofh, 'metallicity', Z_ramses, status)
  end if
  if (allocated(xHII_ramses)) then
    call io_write_table_column(iofh, 'xHII', xHII_ramses, status)
  end if
  if (allocated(xHeII_ramses)) then
    call io_write_table_column(iofh, 'xHeII', xHeII_ramses, status)
  end if
  if (allocated(xHeIII_ramses)) then
    call io_write_table_column(iofh, 'xHeIII', xHeIII_ramses, status)
  end if

  call io_put_keyword(iofh, 'BOXLEN',  boxlen_eff,                 'Simulation box length',     status)
  call io_put_keyword(iofh, 'ORIGINX', -0.5_wp * boxlen_eff,       'Box origin x (centered)',   status)
  call io_put_keyword(iofh, 'ORIGINY', -0.5_wp * boxlen_eff,       'Box origin y (centered)',   status)
  call io_put_keyword(iofh, 'ORIGINZ', -0.5_wp * boxlen_eff,       'Box origin z (centered)',   status)
  call io_put_keyword(iofh, 'NLEAF',   int(nleaf, int32),       'Number of AMR leaf cells', status)
  call io_put_keyword(iofh, 'UNITLCGS', unit_l_cgs,             'RAMSES unit_l in cm',   status)
  call io_put_keyword(iofh, 'UNITPOS', trim(output_unit),       'Position unit in table', status)

  call io_close(iofh, status)
  if (status /= 0) stop 'convert_ramses_to_generic: error while closing output file'

  write(*,'(a)') 'Converted RAMSES snapshot to generic format'
  write(*,'(a)') '  repository : '//trim(repository)
  write(*,'(a,i0)') '  snapnum    : ', snapnum
  write(*,'(a)') '  output     : '//trim(output_file)
  write(*,'(a,i0)') '  nleaf      : ', nleaf
  write(*,'(a,es12.4,1x,a)') '  boxlen     : ', boxlen_eff, trim(output_unit)
  if (compute_physics) then
    write(*,'(a)') '  physics    : xHI, n_e, emissivity computed (CIE formula)'
    if (Z_global >= 0.0_wp) write(*,'(a,es10.3)') '  metallicity: ', Z_global
  end if
  if (allocated(Z_ramses))      write(*,'(a)') '  metallicity: extracted from RAMSES'
  if (allocated(xHII_ramses))   write(*,'(a)') '  xHII       : extracted from RAMSES'
  if (allocated(xHeII_ramses))  write(*,'(a)') '  xHeII      : extracted from RAMSES'
  if (allocated(xHeIII_ramses)) write(*,'(a)') '  xHeIII     : extracted from RAMSES'

contains

  subroutine parse_args(repository, snapnum, output_file, output_unit, compute_physics, Z_global)
    character(len=*), intent(out) :: repository, output_file, output_unit
    integer,          intent(out) :: snapnum
    logical,          intent(out) :: compute_physics
    real(wp),         intent(out) :: Z_global
    integer :: nargs, iarg
    character(len=64) :: cur_arg

    compute_physics = .false.
    Z_global = -1.0_wp
    output_unit = 'kpc'

    nargs = command_argument_count()
    if (nargs < 3) then
      call print_usage()
      stop 1
    end if

    ! Positional arguments
    call get_command_argument(1, repository)
    call get_command_argument(2, arg)
    read(arg, *, err=100) snapnum
    call get_command_argument(3, output_file)

    ! Optional arguments
    iarg = 4
    do while (iarg <= nargs)
      call get_command_argument(iarg, cur_arg)
      select case (trim(cur_arg))
      case ('--compute-physics')
        compute_physics = .true.
      case ('--metallicity')
        iarg = iarg + 1
        if (iarg > nargs) then
          write(*,'(a)') 'ERROR: --metallicity requires a value'
          call print_usage()
          stop 1
        end if
        call get_command_argument(iarg, cur_arg)
        read(cur_arg, *, err=110) Z_global
      case ('cm', 'kpc', 'pc', 'au')
        output_unit = trim(str_lower(cur_arg))
      case default
        output_unit = trim(str_lower(cur_arg))
        if (output_unit /= 'cm' .and. output_unit /= 'kpc' .and. &
            output_unit /= 'pc' .and. output_unit /= 'au') then
          write(*,'(a)') 'Unknown argument: '//trim(cur_arg)
          call print_usage()
          stop 1
        end if
      end select
      iarg = iarg + 1
    end do
    return
100 continue
    write(*,'(a)') 'Invalid snapshot number: '//trim(arg)
    call print_usage()
    stop 1
110 continue
    write(*,'(a)') 'Invalid metallicity value: '//trim(cur_arg)
    call print_usage()
    stop 1
  end subroutine parse_args

  subroutine print_usage()
    write(*,'(a)') 'Usage: convert_ramses_to_generic.x <repository> <snapnum> <output> [unit] [options]'
    write(*,'(a)') '  unit: kpc (default), pc, au, cm'
    write(*,'(a)') '  --compute-physics      Compute xHI, n_e, emissivity from CIE'
    write(*,'(a)') '  --metallicity <Z>      Global metallicity (mass fraction) for dust'
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
      stop 'convert_ramses_to_generic: unsupported unit'
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

end program convert_ramses_to_generic
