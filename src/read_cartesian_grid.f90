module read_cartesian_grid_mod
  !-----------------------------------------------------------------------
  ! Reads an all-in-one Cartesian grid file (HDF5 or FITS) produced by
  ! convert_illustris_to_generic.py --grid-type cartesian.
  !
  ! File layout (HDF5 -- iofile_mod "image section" convention):
  !   /gasDen/data  [nz, ny, nx]   (3D image, mandatory)
  !   /T/data       [nz, ny, nx]   (3D image, mandatory)
  !   /vx/data      [nz, ny, nx]   (3D image, mandatory)
  !   /vy/data      [nz, ny, nx]   (3D image, mandatory)
  !   /vz/data      [nz, ny, nx]   (3D image, mandatory)
  !   + optional sections: xHI, n_e, emissivity, ndust, metallicity
  !   Attributes on first section: NX, NY, NZ, BOXLEN, UNITLCGS, ORIGINX/Y/Z
  !
  ! File layout (FITS -- multi-extension image):
  !   Primary HDU: header with NX, NY, NZ, BOXLEN, UNITLCGS, ORIGINX/Y/Z
  !   Extension 1: gasDen image [nx, ny, nz]
  !   Extension 2: T image
  !   Extension 3: vx image
  !   Extension 4: vy image
  !   Extension 5: vz image
  !   + optional extensions with EXTNAME = xHI, n_e, ...
  !
  ! Public interface:
  !   cartesian_grid_read  -- reads the file, returns 3D arrays
  !-----------------------------------------------------------------------
  use define
  use iofile_mod
  implicit none
  private

  public :: cartesian_grid_read
  public :: cartesian_grid_read_header

contains

  !=========================================================================
  ! Read only the header (dimensions + box size) without loading data.
  ! Safe to call from every MPI rank (light I/O).
  !=========================================================================
  subroutine cartesian_grid_read_header(filename, nx, ny, nz, boxlen, unit_l_cgs)
    character(len=*), intent(in)  :: filename
    integer,          intent(out) :: nx, ny, nz
    real(wp),         intent(out) :: boxlen, unit_l_cgs

    type(io_file_type) :: iofh
    integer :: status
    logical :: found

    status = 0
    call io_open_old(iofh, trim(filename), status)
    if (status /= 0) then
       write(6,'(a,a)') 'cartesian_grid_read_header: cannot open file: ', trim(filename)
       stop 'cartesian_grid_read_header: cannot open file'
    endif

    ! Scan sections to find the one carrying BOXLEN (order-independent).
    ! HDF5 auto-positions at section 1 on open; FITS starts at Primary.
    ! Try current position first, then scan forward.
    boxlen     = 0.0_wp
    unit_l_cgs = 0.0_wp
    found      = .false.

    ! Try current position (HDF5 may be auto-positioned at section 1)
    status = 0
    call io_get_keyword(iofh, 'BOXLEN', boxlen, status)
    if (status /= 0) then
       ! Not found — scan forward through sections (FITS path, or HDF5 non-primary)
       do
          status = 0
          call io_move_to_next_section(iofh, status)
          if (status /= 0) exit
          status = 0
          call io_get_keyword(iofh, 'BOXLEN', boxlen, status)
          if (status == 0) exit
       end do
    endif
    if (boxlen > 0.0_wp) then
       found = .true.
       status = 0; call io_get_keyword(iofh, 'UNITLCGS', unit_l_cgs, status)
       if (status /= 0) unit_l_cgs = 0.0_wp
       status = 0; call io_get_keyword(iofh, 'NAXIS1', nx, status)
       if (status /= 0) stop 'cartesian_grid_read_header: cannot read NAXIS1'
       status = 0; call io_get_keyword(iofh, 'NAXIS2', ny, status)
       if (status /= 0) stop 'cartesian_grid_read_header: cannot read NAXIS2'
       status = 0; call io_get_keyword(iofh, 'NAXIS3', nz, status)
       if (status /= 0) stop 'cartesian_grid_read_header: cannot read NAXIS3'
    endif

    call io_close(iofh, status)

    if (.not. found) stop 'cartesian_grid_read_header: BOXLEN not found in any section'
  end subroutine cartesian_grid_read_header

  subroutine cartesian_grid_read(filename, &
      nx, ny, nz, boxlen, unit_l_cgs, &
      gasDen, T_arr, vx_arr, vy_arr, vz_arr, &
      xHI, n_e, emissivity, ndust, metallicity)
    character(len=*), intent(in)        :: filename
    integer,          intent(out)       :: nx, ny, nz
    real(wp),         intent(out)       :: boxlen, unit_l_cgs
    real(wp), allocatable, intent(out)  :: gasDen(:,:,:)
    real(wp), allocatable, intent(out)  :: T_arr(:,:,:)
    real(wp), allocatable, intent(out)  :: vx_arr(:,:,:)
    real(wp), allocatable, intent(out)  :: vy_arr(:,:,:)
    real(wp), allocatable, intent(out)  :: vz_arr(:,:,:)
    real(wp), allocatable, intent(out), optional :: xHI(:,:,:)
    real(wp), allocatable, intent(out), optional :: n_e(:,:,:)
    real(wp), allocatable, intent(out), optional :: emissivity(:,:,:)
    real(wp), allocatable, intent(out), optional :: ndust(:,:,:)
    real(wp), allocatable, intent(out), optional :: metallicity(:,:,:)

    type(io_file_type) :: iofh
    integer :: status, n1, n2, n3
    character(len=64) :: extname
    logical :: got_header

    status = 0
    call io_open_old(iofh, trim(filename), status)
    if (status /= 0) then
       write(6,'(a,a)') 'cartesian_grid_read: cannot open file: ', trim(filename)
       stop 'cartesian_grid_read: cannot open file'
    endif

    ! Scan all sections by EXTNAME (order-independent).
    ! HDF5 auto-positions at section 1 on open; FITS starts at Primary.
    got_header = .false.
    boxlen     = 0.0_wp
    unit_l_cgs = 0.0_wp
    nx = 0; ny = 0; nz = 0

    ! Try to read from the current position first (HDF5 auto-positioned),
    ! then loop via io_move_to_next_section for remaining sections.
    ! Use a flag to track whether we need the initial move.
    status = 0
    call io_get_keyword(iofh, 'EXTNAME', extname, status)
    if (status /= 0) then
       ! Current position has no EXTNAME (FITS Primary or not positioned).
       ! Move to first real section.
       status = 0
       call io_move_to_next_section(iofh, status)
       if (status /= 0) stop 'cartesian_grid_read: no sections in file'
    endif

    ! Now we are at the first valid section.  Scan all sections.
    do
       extname = ''
       status = 0
       call io_get_keyword(iofh, 'EXTNAME', extname, status)

       ! Read header attributes from whatever section carries BOXLEN
       ! (normally the gasDen section, but accept any).
       if (.not. got_header) then
          status = 0
          call io_get_keyword(iofh, 'BOXLEN', boxlen, status)
          if (status == 0) then
             got_header = .true.
             status = 0; call io_get_keyword(iofh, 'UNITLCGS', unit_l_cgs, status)
             if (status /= 0) unit_l_cgs = 0.0_wp
             status = 0; call io_get_keyword(iofh, 'NAXIS1', n1, status)
             if (status /= 0) stop 'cartesian_grid_read: cannot read NAXIS1'
             status = 0; call io_get_keyword(iofh, 'NAXIS2', n2, status)
             if (status /= 0) stop 'cartesian_grid_read: cannot read NAXIS2'
             status = 0; call io_get_keyword(iofh, 'NAXIS3', n3, status)
             if (status /= 0) stop 'cartesian_grid_read: cannot read NAXIS3'
             nx = n1; ny = n2; nz = n3
          endif
       endif

       ! Dispatch by EXTNAME (order-independent)
       if (nx > 0) then
          select case (trim(extname))
          case ('gasDen')
             allocate(gasDen(nx, ny, nz))
             call io_read_image(iofh, gasDen, status)
          case ('T')
             allocate(T_arr(nx, ny, nz))
             call io_read_image(iofh, T_arr, status)
          case ('vx')
             allocate(vx_arr(nx, ny, nz))
             call io_read_image(iofh, vx_arr, status)
          case ('vy')
             allocate(vy_arr(nx, ny, nz))
             call io_read_image(iofh, vy_arr, status)
          case ('vz')
             allocate(vz_arr(nx, ny, nz))
             call io_read_image(iofh, vz_arr, status)
          case ('xHI')
             if (present(xHI)) then
                allocate(xHI(nx, ny, nz))
                call io_read_image(iofh, xHI, status)
             endif
          case ('n_e')
             if (present(n_e)) then
                allocate(n_e(nx, ny, nz))
                call io_read_image(iofh, n_e, status)
             endif
          case ('emissivity')
             if (present(emissivity)) then
                allocate(emissivity(nx, ny, nz))
                call io_read_image(iofh, emissivity, status)
             endif
          case ('ndust')
             if (present(ndust)) then
                allocate(ndust(nx, ny, nz))
                call io_read_image(iofh, ndust, status)
             endif
          case ('metallicity')
             if (present(metallicity)) then
                allocate(metallicity(nx, ny, nz))
                call io_read_image(iofh, metallicity, status)
             endif
          end select
       endif

       ! Move to next section
       status = 0
       call io_move_to_next_section(iofh, status)
       if (status /= 0) exit   ! no more sections
    end do

    call io_close(iofh, status)

    ! Verify mandatory arrays were found
    if (.not. got_header)      stop 'cartesian_grid_read: BOXLEN not found in any section'
    if (.not. allocated(gasDen)) stop 'cartesian_grid_read: gasDen section not found'
    if (.not. allocated(T_arr))  stop 'cartesian_grid_read: T section not found'
    if (.not. allocated(vx_arr)) stop 'cartesian_grid_read: vx section not found'
    if (.not. allocated(vy_arr)) stop 'cartesian_grid_read: vy section not found'
    if (.not. allocated(vz_arr)) stop 'cartesian_grid_read: vz section not found'
  end subroutine cartesian_grid_read

end module read_cartesian_grid_mod
