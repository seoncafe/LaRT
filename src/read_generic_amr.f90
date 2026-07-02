module read_generic_amr_mod
  !-----------------------------------------------------------------------
  ! Reads generic AMR input files (text, FITS, or HDF5) and returns flat
  ! arrays of leaf-cell data.
  !
  ! Public interface:
  !   generic_amr_read  -- reads generic AMR file, returns leaf data
  !
  ! This module handles ONLY the generic AMR format.  RAMSES binary files
  ! are handled by read_ramses_amr_mod in the standalone converter
  ! (convert_ramses_to_generic.x), NOT in the main LaRT executable.
  !
  ! Supported generic formats:
  !   - Plain text (.dat, .txt)
  !   - FITS binary table (.fits, .fits.gz)
  !   - HDF5 (.h5, .hdf5)
  !
  ! Mandatory columns (9): x, y, z, level, nH, T, vx, vy, vz
  ! Optional columns (detected by name when present):
  !   metallicity, xHI, n_e, n_ion, emissivity, ndust
  !-----------------------------------------------------------------------
  use define
  use iofile_mod
  use, intrinsic :: iso_fortran_env, only: int32
  implicit none
  private

  public :: generic_amr_read

contains

  !=========================================================================
  ! Read generic AMR input in text, FITS, or HDF5 format.
  !
  ! Mandatory outputs (always allocated):
  !   xleaf, yleaf, zleaf  -- leaf center coordinates [code units]
  !   leaf_level           -- AMR level
  !   nH_cgs               -- total hydrogen number density [cm^-3]
  !   T_cgs                -- gas temperature [K]
  !   vx_cgs, vy_cgs, vz_cgs -- gas velocity [km/s]
  !   nleaf                -- number of leaf cells
  !   boxlen_phys          -- box length [code units]
  !
  ! Optional outputs (allocated only if corresponding column found in file):
  !   metallicity          -- Z (mass fraction)
  !   xHI                  -- neutral hydrogen fraction
  !   n_e                  -- electron density [cm^-3]
  !   n_ion                -- scattering ion density [cm^-3]
  !   emissivity           -- Lya emissivity [cm^-3 s^-1]
  !   ndust                -- dust pseudo-number density [cm^-3]
  !=========================================================================
  subroutine generic_amr_read(filename, &
      xleaf, yleaf, zleaf, leaf_level, &
      nH_cgs, T_cgs, vx_cgs, vy_cgs, vz_cgs, &
      nleaf, boxlen_phys, &
      metallicity, xHI, n_e, n_ion, emissivity, ndust, &
      origin_x, origin_y, origin_z)
    character(len=*), intent(in)       :: filename
    real(wp), allocatable, intent(out) :: xleaf(:), yleaf(:), zleaf(:)
    integer,  allocatable, intent(out) :: leaf_level(:)
    real(wp), allocatable, intent(out) :: nH_cgs(:), T_cgs(:)
    real(wp), allocatable, intent(out) :: vx_cgs(:), vy_cgs(:), vz_cgs(:)
    integer,               intent(out) :: nleaf
    real(wp),              intent(out) :: boxlen_phys
    ! Optional extended columns
    real(wp), allocatable, intent(out), optional :: metallicity(:)
    real(wp), allocatable, intent(out), optional :: xHI(:)
    real(wp), allocatable, intent(out), optional :: n_e(:)
    real(wp), allocatable, intent(out), optional :: n_ion(:)
    real(wp), allocatable, intent(out), optional :: emissivity(:)
    real(wp), allocatable, intent(out), optional :: ndust(:)
    ! Box origin (lower corner) in code units
    real(wp),              intent(out), optional :: origin_x, origin_y, origin_z

    ! Reads either:
    !   1) a simple text format
    !        header: nleaf  boxlen
    !        rows  : x, y, z, level, nH [cm^-3], T [K], vx [km/s], vy [km/s], vz [km/s]
    !        optionally more columns if a column-name header is present
    !   2) a FITS/HDF5 binary table with named columns
    !
    ! Position unit is set by par%distance_unit / par%distance2cm (same as
    ! Cartesian mode):
    !   distance_unit = 'kpc'  : x,y,z and boxlen are in kpc  -> distance2cm = kpc2cm
    !   distance_unit = 'pc'   : in pc                        -> distance2cm = pc2cm
    !   distance_unit = 'au'   : in au                        -> distance2cm = au2cm
    !   distance_unit = ''     : already in cm                -> distance2cm = 1
    !   distance2cm = <value>  : 1 data unit = value cm  (overrides distance_unit)
    ! par%distance2cm is always set by setup_v2c before this routine is called.
    real(wp) :: x, y, z, nH, T, vx, vy, vz
    integer  :: lv, n, unit, ios
    real(wp) :: bl
    character(len=512) :: linebuf
    character(len=:), allocatable :: resolved

    ! If `filename` carries no extension, try `.h5`, `.fits.gz`, `.fits`,
    ! `.dat` in that order and use the first one that exists. This mirrors
    ! the Python AMRGrid.read() fallback so users can drop the extension
    ! in par%amr_file.
    resolved = io_resolve_filename(filename)

    if (is_binary_amr_file(resolved)) then
      call generic_amr_read_binary(resolved, &
          xleaf, yleaf, zleaf, leaf_level, &
          nH_cgs, T_cgs, vx_cgs, vy_cgs, vz_cgs, &
          nleaf, boxlen_phys, &
          metallicity, xHI, n_e, n_ion, emissivity, ndust, &
          origin_x, origin_y, origin_z)
      return
    end if

    ! --- Text format ---
    unit = 51
    open(unit, file=resolved, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(6,'(a,a)') 'generic_amr_read: cannot open file: ', resolved
      stop 'generic_amr_read: cannot open file'
    end if

    !--- header: lines starting with '#' (and blank lines) are skipped.
    !--- Accepts either the legacy single header line "<nleaf> <boxlen>"
    !--- or keyword header lines "BOXLEN <value>" / "NLEAF <n>" in any order.
    n  = 0
    bl = -1.0_wp
    do
      read(unit, '(a)', iostat=ios) linebuf
      if (ios /= 0) exit
      linebuf = adjustl(linebuf)
      if (len_trim(linebuf) == 0 .or. linebuf(1:1) == '#') cycle
      if (linebuf(1:6) == 'BOXLEN' .or. linebuf(1:6) == 'boxlen') then
        read(linebuf(7:), *, iostat=ios) bl
        if (ios /= 0) stop 'generic_amr_read: cannot parse BOXLEN header'
      else if (linebuf(1:5) == 'NLEAF' .or. linebuf(1:5) == 'nleaf') then
        read(linebuf(6:), *, iostat=ios) n
        if (ios /= 0) stop 'generic_amr_read: cannot parse NLEAF header'
      else
        !--- legacy header line: "<nleaf> <boxlen>"
        read(linebuf, *, iostat=ios) n, bl
        if (ios /= 0) stop 'generic_amr_read: cannot parse text header'
      end if
      if (n > 0 .and. bl > 0.0_wp) exit
    end do
    if (n <= 0 .or. bl <= 0.0_wp) &
      stop 'generic_amr_read: text header must provide NLEAF and BOXLEN'
    nleaf       = n
    boxlen_phys = bl          ! in code units (kpc, pc, au, or cm) as written in the file

    allocate(xleaf(nleaf), yleaf(nleaf), zleaf(nleaf), leaf_level(nleaf))
    allocate(nH_cgs(nleaf), T_cgs(nleaf), vx_cgs(nleaf), vy_cgs(nleaf), vz_cgs(nleaf))

    !--- data rows (comment / blank lines allowed between rows)
    n = 0
    do while (n < nleaf)
      read(unit, '(a)', iostat=ios) linebuf
      if (ios /= 0) exit
      linebuf = adjustl(linebuf)
      if (len_trim(linebuf) == 0 .or. linebuf(1:1) == '#') cycle
      read(linebuf, *, iostat=ios) x, y, z, lv, nH, T, vx, vy, vz
      if (ios /= 0) exit
      n             = n + 1
      xleaf(n)      = x   ! code units; caller converts to cm via distance2cm
      yleaf(n)      = y
      zleaf(n)      = z
      leaf_level(n) = lv
      nH_cgs(n)     = nH
      T_cgs(n)      = T
      vx_cgs(n)     = vx
      vy_cgs(n)     = vy
      vz_cgs(n)     = vz
    end do
    if (n /= nleaf) then
      write(6,'(a,i0,a,i0)') 'generic_amr_read: expected ', nleaf, ' rows, got ', n
      stop 'generic_amr_read: truncated text file'
    end if
    close(unit)

    ! Text format does not carry origin information; default to a centered
    ! box (origin = -boxlen/2) so legacy text files keep their old behavior.
    if (present(origin_x)) origin_x = -0.5_wp * boxlen_phys
    if (present(origin_y)) origin_y = -0.5_wp * boxlen_phys
    if (present(origin_z)) origin_z = -0.5_wp * boxlen_phys

    ! Text format does not yet support optional columns.
    ! (Future: parse column-name header line and read extra columns.)
  end subroutine generic_amr_read

  !=========================================================================
  ! Read a FITS/HDF5 binary table generic AMR file.
  ! Mandatory columns are read by index (1-9) for backward compatibility.
  ! Optional columns are detected by name.
  !=========================================================================
  subroutine generic_amr_read_binary(filename, &
      xleaf, yleaf, zleaf, leaf_level, &
      nH_cgs, T_cgs, vx_cgs, vy_cgs, vz_cgs, &
      nleaf, boxlen_phys, &
      metallicity, xHI, n_e, n_ion, emissivity, ndust, &
      origin_x_out, origin_y_out, origin_z_out)
    character(len=*), intent(in)       :: filename
    real(wp), allocatable, intent(out) :: xleaf(:), yleaf(:), zleaf(:)
    integer,  allocatable, intent(out) :: leaf_level(:)
    real(wp), allocatable, intent(out) :: nH_cgs(:), T_cgs(:)
    real(wp), allocatable, intent(out) :: vx_cgs(:), vy_cgs(:), vz_cgs(:)
    integer,               intent(out) :: nleaf
    real(wp),              intent(out) :: boxlen_phys
    real(wp), allocatable, intent(out), optional :: metallicity(:)
    real(wp), allocatable, intent(out), optional :: xHI(:)
    real(wp), allocatable, intent(out), optional :: n_e(:)
    real(wp), allocatable, intent(out), optional :: n_ion(:)
    real(wp), allocatable, intent(out), optional :: emissivity(:)
    real(wp), allocatable, intent(out), optional :: ndust(:)
    real(wp),              intent(out), optional :: origin_x_out, origin_y_out, origin_z_out

    type(io_file_type) :: iofh
    integer :: status, colnum
    integer(int32) :: nleaf_i4
    integer(int32), allocatable :: leaf_level_i4(:)
    real(wp) :: origin_x, origin_y, origin_z
    integer :: status_ox, status_oy, status_oz

    status = 0
    call io_open_old(iofh, trim(filename), status)
    if (status /= 0) stop 'generic_amr_read_binary: cannot open file'

    call io_move_to_next_section(iofh, status)
    if (status /= 0) stop 'generic_amr_read_binary: cannot move to binary table HDU'

    call io_get_keyword(iofh, 'NAXIS2', nleaf_i4, status)
    if (status /= 0) stop 'generic_amr_read_binary: cannot read NAXIS2'

    boxlen_phys = 0.0_wp
    call io_get_keyword(iofh, 'BOXLEN',  boxlen_phys, status)
    if (status /= 0) stop 'generic_amr_read_binary: cannot read BOXLEN'

    ! Default origin = centered box (-boxlen/2) when the header is missing.
    origin_x  = -0.5_wp * boxlen_phys
    origin_y  = -0.5_wp * boxlen_phys
    origin_z  = -0.5_wp * boxlen_phys
    status_ox = 0; status_oy = 0; status_oz = 0
    call io_get_keyword(iofh, 'ORIGINX', origin_x, status_ox)
    call io_get_keyword(iofh, 'ORIGINY', origin_y, status_oy)
    call io_get_keyword(iofh, 'ORIGINZ', origin_z, status_oz)
    if (status_ox /= 0) origin_x = -0.5_wp * boxlen_phys
    if (status_oy /= 0) origin_y = -0.5_wp * boxlen_phys
    if (status_oz /= 0) origin_z = -0.5_wp * boxlen_phys

    nleaf = int(nleaf_i4)
    allocate(xleaf(nleaf), yleaf(nleaf), zleaf(nleaf), leaf_level(nleaf))
    allocate(nH_cgs(nleaf), T_cgs(nleaf), vx_cgs(nleaf), vy_cgs(nleaf), vz_cgs(nleaf))
    allocate(leaf_level_i4(nleaf))

    ! --- Read mandatory columns by index (backward compat) ---
    call io_read_table_column(iofh, 1, xleaf,         status)
    call io_read_table_column(iofh, 2, yleaf,         status)
    call io_read_table_column(iofh, 3, zleaf,         status)
    call io_read_table_column(iofh, 4, leaf_level_i4, status)
    call io_read_table_column(iofh, 5, nH_cgs,        status)
    call io_read_table_column(iofh, 6, T_cgs,         status)
    call io_read_table_column(iofh, 7, vx_cgs,        status)
    call io_read_table_column(iofh, 8, vy_cgs,        status)
    call io_read_table_column(iofh, 9, vz_cgs,        status)
    if (status /= 0) stop 'generic_amr_read_binary: cannot read mandatory table columns'

    if (present(origin_x_out)) origin_x_out = origin_x
    if (present(origin_y_out)) origin_y_out = origin_y
    if (present(origin_z_out)) origin_z_out = origin_z

    leaf_level(:) = int(leaf_level_i4(:))
    deallocate(leaf_level_i4)

    ! --- Read optional columns by name ---
    call try_read_optional_column(iofh, nleaf, 'metallicity', metallicity)
    call try_read_optional_column(iofh, nleaf, 'xHI',         xHI)
    call try_read_optional_column(iofh, nleaf, 'n_e',         n_e)
    call try_read_optional_column(iofh, nleaf, 'n_ion',       n_ion)
    call try_read_optional_column(iofh, nleaf, 'emissivity',  emissivity)
    call try_read_optional_column(iofh, nleaf, 'ndust',       ndust)

    call io_close(iofh, status)
  end subroutine generic_amr_read_binary

  !=========================================================================
  ! Try to read an optional column by name.  If the column exists in the
  ! file and the caller passed the optional argument, allocate and fill it.
  ! Otherwise do nothing.
  !=========================================================================
  subroutine try_read_optional_column(iofh, nleaf, colname, arr)
    type(io_file_type), intent(inout) :: iofh
    integer,            intent(in)    :: nleaf
    character(len=*),   intent(in)    :: colname
    real(wp), allocatable, intent(out), optional :: arr(:)

    integer :: colnum, status

    if (.not. present(arr)) return

    ! Try to find the column by name
    status = 0
    call io_get_column_number(iofh, colname, colnum, status)
    if (status /= 0 .or. colnum <= 0) return

    ! Column exists: allocate and read
    allocate(arr(nleaf))
    status = 0
    call io_read_table_column(iofh, colnum, arr, status)
    if (status /= 0) then
      deallocate(arr)
      return
    end if
  end subroutine try_read_optional_column

  !=========================================================================
  ! Detect binary AMR file formats from file extension.
  !=========================================================================
  logical function is_binary_amr_file(filename)
    character(len=*), intent(in) :: filename
    character(len=:), allocatable :: name
    integer :: n

    name = trim(filename)
    n = len(name)
    is_binary_amr_file = .false.
    if (n >= 5) then
      if (name(n-4:n) == '.fits' .or. name(n-4:n) == '.FITS') then
        is_binary_amr_file = .true.
        return
      end if
      if (name(n-4:n) == '.hdf5' .or. name(n-4:n) == '.HDF5') then
        is_binary_amr_file = .true.
        return
      end if
    end if
    if (n >= 8) then
      if (name(n-7:n) == '.fits.gz' .or. name(n-7:n) == '.FITS.GZ') then
        is_binary_amr_file = .true.
        return
      end if
    end if
    if (n >= 3) then
      if (name(n-2:n) == '.h5' .or. name(n-2:n) == '.H5') then
        is_binary_amr_file = .true.
        return
      end if
    end if
  end function is_binary_amr_file

end module read_generic_amr_mod
