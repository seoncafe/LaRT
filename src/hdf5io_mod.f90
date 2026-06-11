module hdf5io_mod
!
! Serial HDF5 backend for the iofile_mod facade. Mirrors fitsio_mod's verb
! set: hdf5_open_new/old/close, hdf5_move_to_next_section,
! hdf5_put_keyword/get_keyword, hdf5_get_column_number,
! hdf5_read_image/append_image, hdf5_read_table_column/append_table_column.
!
! Schema (matches docs/io_layout.md):
!   /section_001/data           -- N-D image dataset (image-style HDU)
!   /section_002/<colname>      -- 1-D table columns (BinTable-style HDU)
!   /section_NNN/@<key>         -- attributes (FITS keywords)
!   /                           -- root attributes (DATE, etc.)
!
! EXTNAME keyword written to a section renames its group to /<EXTNAME>/.
!
! Current limitations:
!   - out_merge (read-old + sum + write-new) is not supported when both
!     ends are HDF5; io_open_old works for plain reads only.
!   - hsize_t / hid_t are assumed to be 8-byte integers (matches HDF5 1.14).
!
#ifdef HDF5
  use hdf5
#endif
  use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int64
  implicit none
  private

  integer, parameter, public :: HDF5_MAX_SECTIONS = 32
  integer, parameter, public :: HDF5_MAX_COLS     = 32
  integer, parameter, public :: HDF5_NAME_LEN     = 64

  type, public :: hdf5_state_type
     integer(int64) :: file_id    = -1_int64
     integer(int64) :: cur_group  = -1_int64
     integer        :: nsections  = 0
     integer        :: cur_idx    = 0    ! 1-based active section
     character(len=HDF5_NAME_LEN) :: section_names(HDF5_MAX_SECTIONS) = ''
     integer        :: ncols(HDF5_MAX_SECTIONS) = 0
     character(len=HDF5_NAME_LEN) :: col_names(HDF5_MAX_COLS, HDF5_MAX_SECTIONS) = ''
     logical        :: writable   = .false.
  end type hdf5_state_type

  public :: hdf5_init, hdf5_finalize
  public :: hdf5_open_new, hdf5_open_old, hdf5_close, hdf5_move_to_next_section
  public :: hdf5_get_column_number
  public :: hdf5_put_key_logical, hdf5_put_key_real32, hdf5_put_key_real64
  public :: hdf5_put_key_int32,   hdf5_put_key_int64,  hdf5_put_key_string
  public :: hdf5_get_key_logical, hdf5_get_key_real32, hdf5_get_key_real64
  public :: hdf5_get_key_int32,   hdf5_get_key_int64,  hdf5_get_key_string
  public :: hdf5_read_1D_real32,  hdf5_read_2D_real32,  hdf5_read_3D_real32,  hdf5_read_4D_real32
  public :: hdf5_read_1D_real64,  hdf5_read_2D_real64,  hdf5_read_3D_real64,  hdf5_read_4D_real64
  public :: hdf5_append_1D_real32,hdf5_append_2D_real32,hdf5_append_3D_real32,hdf5_append_4D_real32
  public :: hdf5_append_1D_real64,hdf5_append_2D_real64,hdf5_append_3D_real64,hdf5_append_4D_real64
  public :: hdf5_append_1D_int32, hdf5_append_2D_int32, hdf5_append_3D_int32, hdf5_append_4D_int32
  public :: hdf5_read_table_column_real32, hdf5_read_table_column_real64
  public :: hdf5_read_table_column_int32,  hdf5_read_table_column_int64
  public :: hdf5_append_table_column_real32, hdf5_append_table_column_real64
  public :: hdf5_append_table_column_int32,  hdf5_append_table_column_int64

  integer, parameter, public :: HDF5_ERR_OK         = 0
  integer, parameter, public :: HDF5_ERR_NOT_BUILT  = -9001
  integer, parameter, public :: HDF5_ERR_TOO_MANY   = -9002
  integer, parameter, public :: HDF5_ERR_NOT_FOUND  = -9003
  integer, parameter, public :: HDF5_ERR_UNSUPPORTED= -9004
  integer, parameter, public :: HDF5_ERR_BACKEND    = -9005

#ifdef HDF5
  logical :: hdf5_initialized = .false.
#endif

contains

#ifdef HDF5
!=============================================================================
!  Lifecycle
!=============================================================================
  subroutine hdf5_init(status)
    integer, intent(inout) :: status
    integer :: ierr
    integer(hid_t), parameter :: H5E_DEFAULT_ID = 0_hid_t
    if (.not. hdf5_initialized) then
       call h5open_f(ierr)
       if (ierr /= 0) then
          status = HDF5_ERR_BACKEND
          return
       endif
       !-- Silence HDF5's automatic stack-trace dumps on errors.  We check
       !   the Fortran ierr from every call ourselves and either recover
       !   (e.g. fall back from CRT_ORDER iteration to NAME iteration when
       !   the file wasn't written with creation-order tracking) or return
       !   HDF5_ERR_BACKEND to the caller.  Without this, every miss spews
       !   a 6-line H5E diagnostic to stderr.
       call h5eset_auto_f(0, ierr)
       hdf5_initialized = .true.
    endif
    status = 0
  end subroutine hdf5_init

  subroutine hdf5_finalize(status)
    integer, intent(inout) :: status
    integer :: ierr
    if (hdf5_initialized) then
       call h5close_f(ierr)
       hdf5_initialized = .false.
       if (ierr /= 0) status = HDF5_ERR_BACKEND
    endif
  end subroutine hdf5_finalize

!=============================================================================
!  Internal helpers
!=============================================================================
  function next_section_name(state) result(nm)
    type(hdf5_state_type), intent(in) :: state
    character(len=HDF5_NAME_LEN) :: nm
    write(nm, '(a,i3.3)') 'section_', state%nsections + 1
  end function next_section_name

  subroutine create_group_tracked(parent_id, path, gid, ierr)
    !-------------------------------------------------------------------------
    ! Create a group whose link-table records the creation order of its
    ! children.  Without this, the alphabetical default index leaks out into
    ! Python (`h5py` iteration), and table-column readers see Jin/Jout/Xfreq
    ! in alphabetical order instead of the order LaRT wrote them.
    !-------------------------------------------------------------------------
    integer(hid_t),    intent(in)  :: parent_id
    character(len=*),  intent(in)  :: path
    integer(hid_t),    intent(out) :: gid
    integer,           intent(out) :: ierr
    integer(hid_t) :: gcpl
    integer :: ierr_aux
    call h5pcreate_f(H5P_GROUP_CREATE_F, gcpl, ierr)
    if (ierr == 0) then
       call h5pset_link_creation_order_f(gcpl, &
            ior(H5P_CRT_ORDER_TRACKED_F, H5P_CRT_ORDER_INDEXED_F), ierr_aux)
    endif
    if (ierr == 0) then
       call h5gcreate_f(parent_id, trim(path), gid, ierr, gcpl_id=gcpl)
    else
       call h5gcreate_f(parent_id, trim(path), gid, ierr)
    endif
    call h5pclose_f(gcpl, ierr_aux)
  end subroutine create_group_tracked

  subroutine begin_new_section(state, status, group_id_out, name_out)
    !-------------------------------------------------------------------------
    ! Create a new group at the file root with an auto-generated name,
    ! register it in the section list, and set it as the current group.
    !-------------------------------------------------------------------------
    type(hdf5_state_type), intent(inout) :: state
    integer, intent(inout) :: status
    integer(hid_t), intent(out), optional :: group_id_out
    character(len=*), intent(out), optional :: name_out
    integer :: ierr
    integer(hid_t) :: gid
    character(len=HDF5_NAME_LEN) :: nm

    if (state%nsections >= HDF5_MAX_SECTIONS) then
       status = HDF5_ERR_TOO_MANY
       return
    endif

    nm = next_section_name(state)
    call create_group_tracked(int(state%file_id, hid_t), '/'//trim(nm), gid, ierr)
    if (ierr /= 0) then
       status = HDF5_ERR_BACKEND
       return
    endif

    ! Close prior current group if it was a section group (not the file root)
    if (state%cur_group > 0 .and. state%cur_group /= state%file_id) then
       call h5gclose_f(int(state%cur_group, hid_t), ierr)
    endif

    state%nsections = state%nsections + 1
    state%section_names(state%nsections) = trim(nm)
    state%ncols(state%nsections) = 0
    state%cur_idx   = state%nsections
    state%cur_group = int(gid, int64)

    if (present(group_id_out)) group_id_out = gid
    if (present(name_out))     name_out     = trim(nm)
    status = 0
  end subroutine begin_new_section

  subroutine ensure_table_section(state, status)
    !-------------------------------------------------------------------------
    ! For table-column appends: if the current section already has at least
    ! one column, reuse it; otherwise create a new section.
    !-------------------------------------------------------------------------
    type(hdf5_state_type), intent(inout) :: state
    integer, intent(inout) :: status
    integer(hid_t) :: gid
    integer :: ierr
    if (state%cur_idx >= 1 .and. state%cur_idx <= state%nsections) then
       if (state%ncols(state%cur_idx) > 0) then
          status = 0
          return
       endif
       ! Current section exists but has no columns yet — reuse it.
       status = 0
       return
    endif
    call begin_new_section(state, status, gid)
  end subroutine ensure_table_section

  function tdims_to_hsize(s) result(h)
    integer, intent(in) :: s
    integer(hsize_t) :: h
    h = int(s, hsize_t)
  end function tdims_to_hsize

  !---------------------------------------------------------------------------
  ! pick_bitpix64: resolve a bitpix=0 (auto) request given the array's
  ! dynamic range AND absolute magnitude. Rules:
  !   max_abs <= 0                  → -32 (all-zero array, nothing to lose)
  !   max_abs < tiny_f32            → -64 (even max would underflow float32)
  !   max_abs > huge_f32            → -64 (max would overflow float32)
  !   min_nz  < tiny_f32            → -64 (smallest non-zero underflows)
  !   max_abs / min_nz >= 1e6       → -64 (dynamic range exceeds float32)
  !   otherwise                     → -32 (float32 storage is enough)
  ! Here tiny_f32 = 1.175e-38 (smallest normal float32) and huge_f32 = 3.40e+38.
  !---------------------------------------------------------------------------
  function pick_bitpix64(max_abs, min_nz) result(bp)
    real(real64), intent(in) :: max_abs, min_nz
    integer :: bp
    real(real64), parameter :: tiny_f32 = 1.175494e-38_real64
    real(real64), parameter :: huge_f32 = 3.402823e+38_real64
    if (max_abs <= 0.0_real64) then
       bp = -32
    else if (max_abs < tiny_f32 .or. max_abs > huge_f32) then
       bp = -64
    else if (min_nz <= 0.0_real64) then
       bp = -32
    else if (min_nz < tiny_f32) then
       bp = -64
    else if (max_abs / min_nz >= 1.0e6_real64) then
       bp = -64
    else
       bp = -32
    endif
  end function pick_bitpix64

  function pick_bitpix32(max_abs, min_nz) result(bp)
    real(real32), intent(in) :: max_abs, min_nz
    integer :: bp
    real(real32), parameter :: tiny_f32 = 1.175494e-38_real32
    real(real32), parameter :: huge_f32 = 3.402823e+38_real32
    if (max_abs <= 0.0_real32) then
       bp = -32
    else if (max_abs < tiny_f32 .or. max_abs > huge_f32) then
       bp = -64
    else if (min_nz <= 0.0_real32) then
       bp = -32
    else if (min_nz < tiny_f32) then
       bp = -64
    else if (max_abs / min_nz >= 1.0e6_real32) then
       bp = -64
    else
       bp = -32
    endif
  end function pick_bitpix32

  !---------------------------------------------------------------------------
  ! Build a dataset creation property list with chunking + gzip(level 4).
  ! Falls back to contiguous (no DCPL) when total elements <= 4096 — small
  ! datasets gain nothing from chunked storage and pay metadata overhead.
  ! Caller must close the returned dcpl with h5pclose_f when > 0.
  !---------------------------------------------------------------------------
  subroutine build_chunked_dcpl(ndim, dims, dcpl_id, use_chunked, ierr)
    integer,           intent(in)  :: ndim
    integer(hsize_t),  intent(in)  :: dims(:)
    integer(hid_t),    intent(out) :: dcpl_id
    logical,           intent(out) :: use_chunked
    integer,           intent(out) :: ierr
    integer(hsize_t) :: chunk_dims(ndim), total
    integer :: i
    integer, parameter :: CHUNK_MAX = 64
    integer, parameter :: SMALL_LIMIT = 4096

    ierr = 0
    dcpl_id = -1_hid_t
    use_chunked = .false.

    total = 1_hsize_t
    do i = 1, ndim
       total = total * dims(i)
    enddo
    if (total <= int(SMALL_LIMIT, hsize_t)) return     ! contiguous

    ! Pick chunk dims: min(dims(i), CHUNK_MAX) per axis.
    ! For 1-D large arrays, use a generous chunk (4096 elts).
    if (ndim == 1) then
       chunk_dims(1) = min(dims(1), int(SMALL_LIMIT, hsize_t))
    else
       do i = 1, ndim
          chunk_dims(i) = min(dims(i), int(CHUNK_MAX, hsize_t))
          if (chunk_dims(i) < 1_hsize_t) chunk_dims(i) = 1_hsize_t
       enddo
    endif

    call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, ierr)
    if (ierr /= 0) return
    call h5pset_chunk_f(dcpl_id, ndim, chunk_dims, ierr)
    if (ierr /= 0) return
    call h5pset_deflate_f(dcpl_id, 4, ierr)   ! gzip level 4
    if (ierr /= 0) return
    use_chunked = .true.
  end subroutine build_chunked_dcpl

  subroutine write_dataset_real32_rank(group_id, name, array, dims, ndim, bitpix_out, status)
    !-------------------------------------------------------------------------
    ! Create an N-D dataset and write `array` to it.
    ! Datasets larger than 4096 elements are chunked + gzip(4); smaller ones
    ! are stored contiguous.  bitpix_out picks the on-disk type.
    !-------------------------------------------------------------------------
    integer(hid_t),    intent(in)    :: group_id
    character(len=*),  intent(in)    :: name
    real(real32),      intent(in)    :: array(*)
    integer(hsize_t),  intent(in)    :: dims(:)
    integer,           intent(in)    :: ndim
    integer,           intent(in)    :: bitpix_out
    integer,           intent(inout) :: status
    integer(hid_t) :: space_id, dset_id, dtype_id, dcpl_id
    integer :: ierr, ierr_dcpl
    logical :: chunked
    if (bitpix_out == -64) then
       dtype_id = H5T_NATIVE_DOUBLE
    else
       dtype_id = H5T_NATIVE_REAL
    endif
    call h5screate_simple_f(ndim, dims(1:ndim), space_id, ierr)
    call build_chunked_dcpl(ndim, dims, dcpl_id, chunked, ierr_dcpl)
    if (chunked) then
       call h5dcreate_f(group_id, trim(name), dtype_id, space_id, dset_id, ierr, dcpl_id)
    else
       call h5dcreate_f(group_id, trim(name), dtype_id, space_id, dset_id, ierr)
    endif
    call h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, dims(1:ndim), ierr)
    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(space_id, ierr)
    if (chunked) call h5pclose_f(dcpl_id, ierr_dcpl)
    if (ierr /= 0) status = HDF5_ERR_BACKEND
  end subroutine write_dataset_real32_rank

  subroutine write_dataset_real64_rank(group_id, name, array, dims, ndim, bitpix_out, status)
    integer(hid_t),    intent(in)    :: group_id
    character(len=*),  intent(in)    :: name
    real(real64),      intent(in)    :: array(*)
    integer(hsize_t),  intent(in)    :: dims(:)
    integer,           intent(in)    :: ndim
    integer,           intent(in)    :: bitpix_out
    integer,           intent(inout) :: status
    integer(hid_t) :: space_id, dset_id, dtype_id, dcpl_id
    integer :: ierr, ierr_dcpl
    logical :: chunked
    if (bitpix_out == -32) then
       dtype_id = H5T_NATIVE_REAL
    else
       dtype_id = H5T_NATIVE_DOUBLE
    endif
    call h5screate_simple_f(ndim, dims(1:ndim), space_id, ierr)
    call build_chunked_dcpl(ndim, dims, dcpl_id, chunked, ierr_dcpl)
    if (chunked) then
       call h5dcreate_f(group_id, trim(name), dtype_id, space_id, dset_id, ierr, dcpl_id)
    else
       call h5dcreate_f(group_id, trim(name), dtype_id, space_id, dset_id, ierr)
    endif
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, dims(1:ndim), ierr)
    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(space_id, ierr)
    if (chunked) call h5pclose_f(dcpl_id, ierr_dcpl)
    if (ierr /= 0) status = HDF5_ERR_BACKEND
  end subroutine write_dataset_real64_rank

  subroutine write_dataset_int32_rank(group_id, name, array, dims, ndim, status)
    integer(hid_t),    intent(in)    :: group_id
    character(len=*),  intent(in)    :: name
    integer(int32),    intent(in)    :: array(*)
    integer(hsize_t),  intent(in)    :: dims(:)
    integer,           intent(in)    :: ndim
    integer,           intent(inout) :: status
    integer(hid_t) :: space_id, dset_id, dcpl_id
    integer :: ierr, ierr_dcpl
    logical :: chunked
    call h5screate_simple_f(ndim, dims(1:ndim), space_id, ierr)
    call build_chunked_dcpl(ndim, dims, dcpl_id, chunked, ierr_dcpl)
    if (chunked) then
       call h5dcreate_f(group_id, trim(name), H5T_NATIVE_INTEGER, space_id, dset_id, ierr, dcpl_id)
    else
       call h5dcreate_f(group_id, trim(name), H5T_NATIVE_INTEGER, space_id, dset_id, ierr)
    endif
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, array, dims(1:ndim), ierr)
    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(space_id, ierr)
    if (chunked) call h5pclose_f(dcpl_id, ierr_dcpl)
    if (ierr /= 0) status = HDF5_ERR_BACKEND
  end subroutine write_dataset_int32_rank

!=============================================================================
!  Open / close
!=============================================================================
  subroutine hdf5_open_new(state, fname, status)
    type(hdf5_state_type), intent(out) :: state
    character(len=*),      intent(in)  :: fname
    integer,             intent(inout) :: status
    integer(hid_t) :: file_id, fcpl_id
    integer :: ierr, ios, tmp_unit

    call hdf5_init(status)
    if (status /= 0) return

    ! Delete pre-existing file (matches fits_open_new semantics)
    open(newunit=tmp_unit, iostat=ios, file=trim(fname), status='old')
    if (ios == 0) close(tmp_unit, status='delete')

    !-- Track link creation order on the root group so readers can iterate
    !   sections in the order they were written rather than alphabetically.
    call h5pcreate_f(H5P_FILE_CREATE_F, fcpl_id, ierr)
    if (ierr == 0) then
       call h5pset_link_creation_order_f(fcpl_id, &
            ior(H5P_CRT_ORDER_TRACKED_F, H5P_CRT_ORDER_INDEXED_F), ierr)
    endif
    if (ierr == 0) then
       call h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, ierr, &
                        creation_prp=fcpl_id)
    else
       call h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, ierr)
    endif
    call h5pclose_f(fcpl_id, ierr)
    if (ierr /= 0) then
       status = HDF5_ERR_BACKEND
       return
    endif
    state%file_id  = int(file_id, int64)
    state%cur_group= state%file_id     ! attributes go on root until first section
    state%cur_idx  = 0
    state%nsections= 0
    state%writable = .true.
    status = 0
  end subroutine hdf5_open_new

  subroutine hdf5_open_old(state, fname, status, writable)
    type(hdf5_state_type), intent(out) :: state
    character(len=*),      intent(in)  :: fname
    integer,             intent(inout) :: status
    logical, optional,     intent(in)  :: writable
    integer(hid_t) :: file_id
    integer :: ierr
    integer :: access_flag

    call hdf5_init(status)
    if (status /= 0) return

    access_flag = H5F_ACC_RDONLY_F
    if (present(writable)) then
       if (writable) access_flag = H5F_ACC_RDWR_F
    endif

    call h5fopen_f(trim(fname), access_flag, file_id, ierr)
    if (ierr /= 0) then
       status = HDF5_ERR_BACKEND
       return
    endif
    state%file_id = int(file_id, int64)

    ! Build the section list by iterating root group children in creation order.
    call collect_sections(state, status)
    if (status /= 0) return

    state%writable  = (access_flag == H5F_ACC_RDWR_F)

    !--- Position the cursor to mirror FITS open_old semantics.
    !
    ! In FITS, after open, cur HDU = 1 (Primary). For image-primary writes
    ! (peel-off, sight-line), CFITSIO folds the first image into Primary, so
    ! HDU 1 already holds image data + keywords. For table-primary writes
    ! (spectrum), Primary stays empty and HDU 2 holds the BinTable+keywords;
    ! the caller does an explicit fits_move_to_next_hdu first.
    !
    ! To reproduce this in HDF5, auto-position cur_group at section 1 IFF
    ! section 1 looks image-primary (single dataset named 'data'); otherwise
    ! leave cur_group at root so the spectrum-style caller's first
    ! io_move_to_next_section call still lands on /Spectrum.
    if (state%nsections >= 1 .and. state%ncols(1) == 1 .and. &
        trim(state%col_names(1,1)) == 'data') then
       call h5gopen_f(int(state%file_id, hid_t), &
                      '/'//trim(state%section_names(1)), state%cur_group, ierr)
       state%cur_idx = 1
    else
       state%cur_idx   = 0
       state%cur_group = state%file_id
    endif
    status = 0
  end subroutine hdf5_open_old

  subroutine collect_sections(state, status)
    !-------------------------------------------------------------------------
    ! Walk the root group's children in creation order so the section list
    ! reflects the order they were written (matches FITS HDU sequence).
    !-------------------------------------------------------------------------
    type(hdf5_state_type), intent(inout) :: state
    integer, intent(inout) :: status
    integer(hsize_t) :: idx, nlinks_h
    integer :: ierr
    integer(size_t) :: name_size
    integer :: storage_type, nlinks, max_corder
    character(len=HDF5_NAME_LEN) :: nm
    integer(hid_t) :: gid, oid
    integer :: type_

    call h5gget_info_f(int(state%file_id, hid_t), storage_type, nlinks, max_corder, ierr)
    if (ierr /= 0) then; status = HDF5_ERR_BACKEND; return; endif

    state%nsections = 0
    do idx = 0_hsize_t, int(nlinks - 1, hsize_t)
       name_size = HDF5_NAME_LEN
       nm = ''
       call h5lget_name_by_idx_f(int(state%file_id, hid_t), '/', &
                                 H5_INDEX_CRT_ORDER_F, H5_ITER_INC_F, &
                                 idx, nm, ierr, size=name_size)
       if (ierr /= 0) then
          ! Fall back: alphabetical index if creation order is unavailable.
          call h5lget_name_by_idx_f(int(state%file_id, hid_t), '/', &
                                    H5_INDEX_NAME_F, H5_ITER_INC_F, &
                                    idx, nm, ierr, size=name_size)
          if (ierr /= 0) cycle
       endif

       ! Skip non-group children at root (datasets, soft links, etc.)
       call h5oopen_f(int(state%file_id, hid_t), trim(nm), oid, ierr)
       if (ierr /= 0) cycle
       call h5iget_type_f(oid, type_, ierr)
       call h5oclose_f(oid, ierr)
       if (type_ /= H5I_GROUP_F) cycle

       if (state%nsections >= HDF5_MAX_SECTIONS) exit
       state%nsections = state%nsections + 1
       state%section_names(state%nsections) = trim(nm)

       call h5gopen_f(int(state%file_id, hid_t), trim(nm), gid, ierr)
       if (ierr == 0) then
          call list_datasets_in_group(gid, state, state%nsections)
          call h5gclose_f(gid, ierr)
       endif
    enddo
    status = 0
  end subroutine collect_sections

  subroutine list_datasets_in_group(gid, state, sec_idx)
    integer(hid_t),        intent(in)    :: gid
    type(hdf5_state_type), intent(inout) :: state
    integer,               intent(in)    :: sec_idx
    integer(hsize_t) :: k
    integer(size_t) :: name_size
    integer :: storage_type, nlinks, max_corder
    integer :: ierr, type_
    character(len=HDF5_NAME_LEN) :: dname
    integer(hid_t) :: oid

    call h5gget_info_f(gid, storage_type, nlinks, max_corder, ierr)
    if (ierr /= 0) return

    state%ncols(sec_idx) = 0
    do k = 0_hsize_t, int(nlinks - 1, hsize_t)
       name_size = HDF5_NAME_LEN
       dname = ''
       call h5lget_name_by_idx_f(gid, '.', H5_INDEX_CRT_ORDER_F, H5_ITER_INC_F, &
                                 k, dname, ierr, size=name_size)
       if (ierr /= 0) then
          call h5lget_name_by_idx_f(gid, '.', H5_INDEX_NAME_F, H5_ITER_INC_F, &
                                    k, dname, ierr, size=name_size)
          if (ierr /= 0) cycle
       endif
       call h5oopen_f(gid, trim(dname), oid, ierr)
       if (ierr /= 0) cycle
       call h5iget_type_f(oid, type_, ierr)
       call h5oclose_f(oid, ierr)
       if (type_ /= H5I_DATASET_F) cycle
       if (state%ncols(sec_idx) >= HDF5_MAX_COLS) exit
       state%ncols(sec_idx) = state%ncols(sec_idx) + 1
       state%col_names(state%ncols(sec_idx), sec_idx) = trim(dname)
    enddo
  end subroutine list_datasets_in_group

  subroutine hdf5_close(state, status)
    type(hdf5_state_type), intent(inout) :: state
    integer,               intent(inout) :: status
    integer :: ierr
    if (state%cur_group > 0 .and. state%cur_group /= state%file_id) then
       call h5gclose_f(int(state%cur_group, hid_t), ierr)
    endif
    if (state%file_id > 0) then
       call h5fclose_f(int(state%file_id, hid_t), ierr)
       if (ierr /= 0) status = HDF5_ERR_BACKEND
    endif
    state%file_id   = -1_int64
    state%cur_group = -1_int64
    state%nsections = 0
    state%cur_idx   = 0
    state%writable  = .false.
  end subroutine hdf5_close

  subroutine hdf5_move_to_next_section(state, status, nstep)
    type(hdf5_state_type), intent(inout) :: state
    integer,               intent(inout) :: status
    integer, optional,     intent(in)    :: nstep
    integer :: step, target_idx, ierr
    integer(hid_t) :: gid
    step = 1
    if (present(nstep)) then
       if (nstep >= 0) step = nstep
    endif
    target_idx = state%cur_idx + step
    if (target_idx < 1) target_idx = 1
    if (target_idx > state%nsections) then
       status = HDF5_ERR_NOT_FOUND
       return
    endif
    if (state%cur_group > 0 .and. state%cur_group /= state%file_id) then
       call h5gclose_f(int(state%cur_group, hid_t), ierr)
    endif
    call h5gopen_f(int(state%file_id, hid_t), &
                   '/'//trim(state%section_names(target_idx)), gid, ierr)
    if (ierr /= 0) then; status = HDF5_ERR_BACKEND; return; endif
    state%cur_group = int(gid, int64)
    state%cur_idx   = target_idx
    status = 0
  end subroutine hdf5_move_to_next_section

!=============================================================================
!  Column lookup
!=============================================================================
  subroutine hdf5_get_column_number(state, colname, colnum, status)
    type(hdf5_state_type), intent(in)    :: state
    character(len=*),      intent(in)    :: colname
    integer,               intent(out)   :: colnum
    integer,               intent(inout) :: status
    integer :: i
    colnum = 0
    if (state%cur_idx < 1 .or. state%cur_idx > state%nsections) then
       status = HDF5_ERR_NOT_FOUND; return
    endif
    do i = 1, state%ncols(state%cur_idx)
       if (trim(state%col_names(i, state%cur_idx)) == trim(colname)) then
          colnum = i
          status = 0
          return
       endif
    enddo
    status = HDF5_ERR_NOT_FOUND
  end subroutine hdf5_get_column_number

!=============================================================================
!  Attribute write (= FITS put_keyword)
!=============================================================================
  subroutine put_attr_scalar(loc_id, name, dtype, buf, status)
    integer(hid_t),   intent(in)    :: loc_id
    character(len=*), intent(in)    :: name
    integer(hid_t),   intent(in)    :: dtype
    integer(hid_t),   intent(in)    :: buf  ! dummy address, callers pass through
    integer,          intent(inout) :: status
    ! unused — real attribute writes are below.
    status = 0
  end subroutine put_attr_scalar

  subroutine attr_open_or_create(loc_id, name, dtype, attr_id, ierr)
    integer(hid_t),    intent(in)  :: loc_id
    character(len=*),  intent(in)  :: name
    integer(hid_t),    intent(in)  :: dtype
    integer(hid_t),    intent(out) :: attr_id
    integer,           intent(out) :: ierr
    integer(hid_t)   :: space_id
    integer(hsize_t) :: dims(1)
    logical :: exists
    dims(1) = 1
    call h5aexists_f(loc_id, trim(name), exists, ierr)
    if (exists) then
       call h5adelete_f(loc_id, trim(name), ierr)
    endif
    call h5screate_simple_f(1, dims, space_id, ierr)
    call h5acreate_f(loc_id, trim(name), dtype, space_id, attr_id, ierr)
    call h5sclose_f(space_id, ierr)
  end subroutine attr_open_or_create

  subroutine hdf5_put_key_logical(state, name, value, comment, status)
    type(hdf5_state_type), intent(in)    :: state
    character(len=*),      intent(in)    :: name
    logical,               intent(in)    :: value
    character(len=*),      intent(in)    :: comment
    integer,               intent(inout) :: status
    integer(hid_t) :: aid
    integer(hsize_t) :: dims(1)
    integer :: ierr, val_i
    val_i = 0
    if (value) val_i = 1
    dims(1) = 1
    call attr_open_or_create(int(state%cur_group, hid_t), trim(name), H5T_NATIVE_INTEGER, aid, ierr)
    call h5awrite_f(aid, H5T_NATIVE_INTEGER, val_i, dims, ierr)
    call h5aclose_f(aid, ierr)
    if (ierr /= 0) status = HDF5_ERR_BACKEND
  end subroutine hdf5_put_key_logical

  subroutine hdf5_put_key_real32(state, name, value, comment, status)
    type(hdf5_state_type), intent(in)    :: state
    character(len=*),      intent(in)    :: name
    real(real32),          intent(in)    :: value
    character(len=*),      intent(in)    :: comment
    integer,               intent(inout) :: status
    integer(hid_t) :: aid
    integer(hsize_t) :: dims(1)
    integer :: ierr
    dims(1) = 1
    call attr_open_or_create(int(state%cur_group, hid_t), trim(name), H5T_NATIVE_REAL, aid, ierr)
    call h5awrite_f(aid, H5T_NATIVE_REAL, value, dims, ierr)
    call h5aclose_f(aid, ierr)
    if (ierr /= 0) status = HDF5_ERR_BACKEND
  end subroutine hdf5_put_key_real32

  subroutine hdf5_put_key_real64(state, name, value, comment, status)
    type(hdf5_state_type), intent(in)    :: state
    character(len=*),      intent(in)    :: name
    real(real64),          intent(in)    :: value
    character(len=*),      intent(in)    :: comment
    integer,               intent(inout) :: status
    integer(hid_t) :: aid
    integer(hsize_t) :: dims(1)
    integer :: ierr
    dims(1) = 1
    call attr_open_or_create(int(state%cur_group, hid_t), trim(name), H5T_NATIVE_DOUBLE, aid, ierr)
    call h5awrite_f(aid, H5T_NATIVE_DOUBLE, value, dims, ierr)
    call h5aclose_f(aid, ierr)
    if (ierr /= 0) status = HDF5_ERR_BACKEND
  end subroutine hdf5_put_key_real64

  subroutine hdf5_put_key_int32(state, name, value, comment, status)
    type(hdf5_state_type), intent(in)    :: state
    character(len=*),      intent(in)    :: name
    integer(int32),        intent(in)    :: value
    character(len=*),      intent(in)    :: comment
    integer,               intent(inout) :: status
    integer(hid_t) :: aid
    integer(hsize_t) :: dims(1)
    integer :: ierr
    dims(1) = 1
    call attr_open_or_create(int(state%cur_group, hid_t), trim(name), H5T_NATIVE_INTEGER, aid, ierr)
    call h5awrite_f(aid, H5T_NATIVE_INTEGER, value, dims, ierr)
    call h5aclose_f(aid, ierr)
    if (ierr /= 0) status = HDF5_ERR_BACKEND
  end subroutine hdf5_put_key_int32

  subroutine hdf5_put_key_int64(state, name, value, comment, status)
    type(hdf5_state_type), intent(in)    :: state
    character(len=*),      intent(in)    :: name
    integer(int64),        intent(in)    :: value
    character(len=*),      intent(in)    :: comment
    integer,               intent(inout) :: status
    integer(hid_t) :: aid
    integer(hsize_t) :: dims(1)
    integer :: ierr
    dims(1) = 1
    call attr_open_or_create(int(state%cur_group, hid_t), trim(name), H5T_STD_I64LE, aid, ierr)
    call h5awrite_f(aid, H5T_STD_I64LE, value, dims, ierr)
    call h5aclose_f(aid, ierr)
    if (ierr /= 0) status = HDF5_ERR_BACKEND
  end subroutine hdf5_put_key_int64

  subroutine hdf5_put_key_string(state, name, value, comment, status)
    type(hdf5_state_type), intent(inout) :: state
    character(len=*),      intent(in)    :: name
    character(len=*),      intent(in)    :: value
    character(len=*),      intent(in)    :: comment
    integer,               intent(inout) :: status
    integer(hid_t) :: aid, tid
    integer(hsize_t) :: dims(1)
    integer :: ierr
    integer(size_t) :: slen

    ! Special-case EXTNAME: rename the current section group.
    if (trim(name) == 'EXTNAME' .or. trim(name) == 'extname') then
       call rename_current_section(state, trim(value), status)
       if (status /= 0) return
    endif

    slen = max(len_trim(value), 1)
    call h5tcopy_f(H5T_NATIVE_CHARACTER, tid, ierr)
    call h5tset_size_f(tid, slen, ierr)
    dims(1) = 1
    call attr_open_or_create(int(state%cur_group, hid_t), trim(name), tid, aid, ierr)
    call h5awrite_f(aid, tid, trim(value), dims, ierr)
    call h5aclose_f(aid, ierr)
    call h5tclose_f(tid, ierr)
    if (ierr /= 0) status = HDF5_ERR_BACKEND
  end subroutine hdf5_put_key_string

  subroutine rename_current_section(state, newname, status)
    type(hdf5_state_type), intent(inout) :: state
    character(len=*),      intent(in)    :: newname
    integer,               intent(inout) :: status
    integer :: ierr
    character(len=HDF5_NAME_LEN) :: oldname
    if (state%cur_idx < 1 .or. state%cur_idx > state%nsections) then
       status = HDF5_ERR_NOT_FOUND; return
    endif
    oldname = state%section_names(state%cur_idx)
    if (trim(oldname) == trim(newname)) return
    ! Close current group, rename via H5Lmove, reopen.
    if (state%cur_group > 0) call h5gclose_f(int(state%cur_group, hid_t), ierr)
    call h5lmove_f(int(state%file_id, hid_t), '/'//trim(oldname), &
                   int(state%file_id, hid_t), '/'//trim(newname), ierr)
    if (ierr /= 0) then
       status = HDF5_ERR_BACKEND
       ! reopen the original to keep state consistent
       call h5gopen_f(int(state%file_id, hid_t), '/'//trim(oldname), state%cur_group, ierr)
       return
    endif
    state%section_names(state%cur_idx) = trim(newname)
    call h5gopen_f(int(state%file_id, hid_t), '/'//trim(newname), state%cur_group, ierr)
  end subroutine rename_current_section

!=============================================================================
!  Attribute read (= FITS get_keyword)
!=============================================================================
  subroutine attr_open(loc_id, name, attr_id, ierr)
    integer(hid_t),   intent(in)  :: loc_id
    character(len=*), intent(in)  :: name
    integer(hid_t),   intent(out) :: attr_id
    integer,          intent(out) :: ierr
    call h5aopen_f(loc_id, trim(name), attr_id, ierr)
  end subroutine attr_open

  !---------------------------------------------------------------------------
  ! try_naxis_fallback: when a FITS keyword query asks for NAXIS / NAXIS1 /
  ! NAXIS2 / NAXIS3 / NAXIS4 and there is no such attribute on the current
  ! group, infer the answer from the 'data' dataset's shape.  This lets the
  ! existing read_fits_data code work unchanged when fed an HDF5 input file.
  ! Returns .true. if the name matched and dim_out was populated.
  !---------------------------------------------------------------------------
  function try_naxis_fallback(state, name, dim_out) result(handled)
    type(hdf5_state_type), intent(in)  :: state
    character(len=*),      intent(in)  :: name
    integer(int64),        intent(out) :: dim_out
    logical :: handled
    integer(hid_t) :: dset_id, space_id
    integer(hsize_t) :: dims(7), maxdims(7)
    integer :: ndim, ierr, axis
    character(len=16) :: ukey
    integer :: i

    handled = .false.
    dim_out = 0_int64

    ! Uppercase the requested name for comparison.
    ukey = ''
    do i = 1, min(len_trim(name), 16)
       if (name(i:i) >= 'a' .and. name(i:i) <= 'z') then
          ukey(i:i) = char(ichar(name(i:i)) - 32)
       else
          ukey(i:i) = name(i:i)
       endif
    enddo

    if (trim(ukey) == 'NAXIS') then
       axis = 0   ! return rank
    else if (trim(ukey) == 'NAXIS1') then
       axis = 1
    else if (trim(ukey) == 'NAXIS2') then
       axis = 2
    else if (trim(ukey) == 'NAXIS3') then
       axis = 3
    else if (trim(ukey) == 'NAXIS4') then
       axis = 4
    else
       return   ! not a NAXIS query
    endif

    call h5dopen_f(int(state%cur_group, hid_t), 'data', dset_id, ierr)
    if (ierr /= 0) return   ! no fallback available; treat as NOT_FOUND
    call h5dget_space_f(dset_id, space_id, ierr)
    if (ierr /= 0) then
       call h5dclose_f(dset_id, ierr)
       return
    endif
    call h5sget_simple_extent_dims_f(space_id, dims, maxdims, ndim)
    call h5sclose_f(space_id, ierr)
    call h5dclose_f(dset_id, ierr)
    if (ndim < 0) return

    if (axis == 0) then
       dim_out = int(ndim, int64)
       handled = .true.
    else if (axis >= 1 .and. axis <= ndim) then
       dim_out = int(dims(axis), int64)
       handled = .true.
    endif
  end function try_naxis_fallback

  subroutine hdf5_get_key_logical(state, name, value, status, comment)
    type(hdf5_state_type),      intent(in)    :: state
    character(len=*),           intent(in)    :: name
    logical,                    intent(out)   :: value
    integer,                    intent(inout) :: status
    character(len=*), optional, intent(out)   :: comment
    integer(hid_t) :: aid
    integer(hsize_t) :: dims(1)
    integer :: ierr, val_i
    dims(1) = 1
    call attr_open(int(state%cur_group, hid_t), trim(name), aid, ierr)
    if (ierr /= 0) then; value = .false.; status = HDF5_ERR_NOT_FOUND; return; endif
    call h5aread_f(aid, H5T_NATIVE_INTEGER, val_i, dims, ierr)
    call h5aclose_f(aid, ierr)
    value = (val_i /= 0)
    if (present(comment)) comment = ''
  end subroutine hdf5_get_key_logical

  subroutine hdf5_get_key_real32(state, name, value, status, comment)
    type(hdf5_state_type),      intent(in)    :: state
    character(len=*),           intent(in)    :: name
    real(real32),               intent(out)   :: value
    integer,                    intent(inout) :: status
    character(len=*), optional, intent(out)   :: comment
    integer(hid_t) :: aid
    integer(hsize_t) :: dims(1)
    integer :: ierr
    dims(1) = 1
    call attr_open(int(state%cur_group, hid_t), trim(name), aid, ierr)
    if (ierr /= 0) then; value = 0.0_real32; status = HDF5_ERR_NOT_FOUND; return; endif
    call h5aread_f(aid, H5T_NATIVE_REAL, value, dims, ierr)
    call h5aclose_f(aid, ierr)
    if (present(comment)) comment = ''
  end subroutine hdf5_get_key_real32

  subroutine hdf5_get_key_real64(state, name, value, status, comment)
    type(hdf5_state_type),      intent(in)    :: state
    character(len=*),           intent(in)    :: name
    real(real64),               intent(out)   :: value
    integer,                    intent(inout) :: status
    character(len=*), optional, intent(out)   :: comment
    integer(hid_t) :: aid
    integer(hsize_t) :: dims(1)
    integer :: ierr
    dims(1) = 1
    call attr_open(int(state%cur_group, hid_t), trim(name), aid, ierr)
    if (ierr /= 0) then; value = 0.0_real64; status = HDF5_ERR_NOT_FOUND; return; endif
    call h5aread_f(aid, H5T_NATIVE_DOUBLE, value, dims, ierr)
    call h5aclose_f(aid, ierr)
    if (present(comment)) comment = ''
  end subroutine hdf5_get_key_real64

  subroutine hdf5_get_key_int32(state, name, value, status, comment)
    type(hdf5_state_type),      intent(in)    :: state
    character(len=*),           intent(in)    :: name
    integer(int32),             intent(out)   :: value
    integer,                    intent(inout) :: status
    character(len=*), optional, intent(out)   :: comment
    integer(hid_t) :: aid
    integer(hsize_t) :: dims(1)
    integer :: ierr
    integer(int64) :: fallback_val
    dims(1) = 1
    call attr_open(int(state%cur_group, hid_t), trim(name), aid, ierr)
    if (ierr /= 0) then
       ! Fall back to dataset-shape lookup for NAXIS-style queries.
       if (try_naxis_fallback(state, name, fallback_val)) then
          value  = int(fallback_val, int32)
          status = 0
          if (present(comment)) comment = ''
          return
       endif
       value = 0_int32
       status = HDF5_ERR_NOT_FOUND
       return
    endif
    call h5aread_f(aid, H5T_NATIVE_INTEGER, value, dims, ierr)
    call h5aclose_f(aid, ierr)
    if (present(comment)) comment = ''
  end subroutine hdf5_get_key_int32

  subroutine hdf5_get_key_int64(state, name, value, status, comment)
    type(hdf5_state_type),      intent(in)    :: state
    character(len=*),           intent(in)    :: name
    integer(int64),             intent(out)   :: value
    integer,                    intent(inout) :: status
    character(len=*), optional, intent(out)   :: comment
    integer(hid_t) :: aid
    integer(hsize_t) :: dims(1)
    integer :: ierr
    integer(int64) :: fallback_val
    dims(1) = 1
    call attr_open(int(state%cur_group, hid_t), trim(name), aid, ierr)
    if (ierr /= 0) then
       if (try_naxis_fallback(state, name, fallback_val)) then
          value  = fallback_val
          status = 0
          if (present(comment)) comment = ''
          return
       endif
       value = 0_int64
       status = HDF5_ERR_NOT_FOUND
       return
    endif
    call h5aread_f(aid, H5T_STD_I64LE, value, dims, ierr)
    call h5aclose_f(aid, ierr)
    if (present(comment)) comment = ''
  end subroutine hdf5_get_key_int64

  subroutine hdf5_get_key_string(state, name, value, status, comment)
    type(hdf5_state_type),      intent(in)    :: state
    character(len=*),           intent(in)    :: name
    character(len=*),           intent(out)   :: value
    integer,                    intent(inout) :: status
    character(len=*), optional, intent(out)   :: comment
    integer(hid_t)   :: aid, tid, atype
    integer(hsize_t) :: dims(1)
    integer          :: ierr, ii
    integer(size_t)  :: slen
    logical          :: is_vlen
    dims(1) = 1
    value = ''
    call attr_open(int(state%cur_group, hid_t), trim(name), aid, ierr)
    if (ierr /= 0) then; status = HDF5_ERR_NOT_FOUND; return; endif

    !--- Detect variable-length string attributes (e.g. h5py's default
    !    string_dtype).  Reading vlen into a Fortran character buffer via
    !    h5aread_f returns garbage on most Fortran/HDF5 bindings, so we
    !    treat vlen as unreadable here and let the caller treat the value
    !    as missing.  External pipelines should write fixed-length strings
    !    (h5py: dtype=h5py.string_dtype(encoding='ascii', length=N) ).
    call h5aget_type_f(aid, atype, ierr)
    call h5tis_variable_str_f(atype, is_vlen, ierr)
    call h5tclose_f(atype, ierr)
    if (is_vlen) then
       write(*,'(3a)') ' HDF5: WARNING -- string attribute ''', trim(name), &
            ''' is variable-length; LaRT requires fixed-length strings ' // &
            '(use h5py.string_dtype(length=N) when writing).  Skipping.'
       call h5aclose_f(aid, ierr)
       value = ''
       status = HDF5_ERR_NOT_FOUND
       if (present(comment)) comment = ''
       return
    end if

    !--- Fixed-length read: build a Fortran-buffer-sized fixed type.
    !    HDF5 NULL-pads strings shorter than the buffer; we strip those
    !    padding NULs below so callers can trim() normally.
    call h5tcopy_f(H5T_NATIVE_CHARACTER, tid, ierr)
    slen = int(len(value), size_t)
    if (slen < 1) slen = 1
    call h5tset_size_f(tid, slen, ierr)
    call h5aread_f(aid, tid, value, dims, ierr)
    call h5tclose_f(tid, ierr)
    call h5aclose_f(aid, ierr)
    if (ierr /= 0) status = HDF5_ERR_BACKEND

    do ii = 1, len(value)
       if (value(ii:ii) == char(0)) value(ii:ii) = ' '
    end do

    if (present(comment)) comment = ''
  end subroutine hdf5_get_key_string

!=============================================================================
!  Image append (always goes into /section_NNN/data)
!=============================================================================
#define HDF5_APPEND_IMAGE_BODY(rank_, kind_macro_, h5kind_) \
    integer :: ierr; \
    integer(hsize_t) :: dims(rank_); \
    integer(hid_t) :: gid; \
    call begin_new_section(state, status, gid); \
    if (status /= 0) return;

  subroutine hdf5_append_1D_real32(state, array, status, bitpix)
    type(hdf5_state_type), intent(inout) :: state
    real(real32),          intent(in)    :: array(:)
    integer,               intent(inout) :: status
    integer, optional,     intent(in)    :: bitpix
    integer(hsize_t) :: dims(1)
    integer(hid_t) :: gid
    integer :: bp
    real(real32) :: max_abs, min_nz
    bp = -32; if (present(bitpix)) bp = bitpix
    if (bp == 0) then
       max_abs = maxval(abs(array))
       if (max_abs > 0.0_real32) then
          min_nz = minval(abs(array), mask=abs(array) > 0.0_real32)
       else
          min_nz = 0.0_real32
       endif
       bp = pick_bitpix32(max_abs, min_nz)
    endif
    call begin_new_section(state, status, gid); if (status /= 0) return
    dims(1) = size(array)
    call write_dataset_real32_rank(gid, 'data', array, dims, 1, bp, status)
  end subroutine hdf5_append_1D_real32

  subroutine hdf5_append_2D_real32(state, array, status, bitpix)
    type(hdf5_state_type), intent(inout) :: state
    real(real32),          intent(in)    :: array(:,:)
    integer,               intent(inout) :: status
    integer, optional,     intent(in)    :: bitpix
    integer(hsize_t) :: dims(2)
    integer(hid_t) :: gid
    integer :: bp
    real(real32) :: max_abs, min_nz
    bp = -32; if (present(bitpix)) bp = bitpix
    if (bp == 0) then
       max_abs = maxval(abs(array))
       if (max_abs > 0.0_real32) then
          min_nz = minval(abs(array), mask=abs(array) > 0.0_real32)
       else
          min_nz = 0.0_real32
       endif
       bp = pick_bitpix32(max_abs, min_nz)
    endif
    call begin_new_section(state, status, gid); if (status /= 0) return
    dims(1) = size(array,1); dims(2) = size(array,2)
    call write_dataset_real32_rank(gid, 'data', array, dims, 2, bp, status)
  end subroutine hdf5_append_2D_real32

  subroutine hdf5_append_3D_real32(state, array, status, bitpix)
    type(hdf5_state_type), intent(inout) :: state
    real(real32),          intent(in)    :: array(:,:,:)
    integer,               intent(inout) :: status
    integer, optional,     intent(in)    :: bitpix
    integer(hsize_t) :: dims(3)
    integer(hid_t) :: gid
    integer :: bp
    real(real32) :: max_abs, min_nz
    bp = -32; if (present(bitpix)) bp = bitpix
    if (bp == 0) then
       max_abs = maxval(abs(array))
       if (max_abs > 0.0_real32) then
          min_nz = minval(abs(array), mask=abs(array) > 0.0_real32)
       else
          min_nz = 0.0_real32
       endif
       bp = pick_bitpix32(max_abs, min_nz)
    endif
    call begin_new_section(state, status, gid); if (status /= 0) return
    dims(1) = size(array,1); dims(2) = size(array,2); dims(3) = size(array,3)
    call write_dataset_real32_rank(gid, 'data', array, dims, 3, bp, status)
  end subroutine hdf5_append_3D_real32

  subroutine hdf5_append_4D_real32(state, array, status, bitpix)
    type(hdf5_state_type), intent(inout) :: state
    real(real32),          intent(in)    :: array(:,:,:,:)
    integer,               intent(inout) :: status
    integer, optional,     intent(in)    :: bitpix
    integer(hsize_t) :: dims(4)
    integer(hid_t) :: gid
    integer :: bp
    real(real32) :: max_abs, min_nz
    bp = -32; if (present(bitpix)) bp = bitpix
    if (bp == 0) then
       max_abs = maxval(abs(array))
       if (max_abs > 0.0_real32) then
          min_nz = minval(abs(array), mask=abs(array) > 0.0_real32)
       else
          min_nz = 0.0_real32
       endif
       bp = pick_bitpix32(max_abs, min_nz)
    endif
    call begin_new_section(state, status, gid); if (status /= 0) return
    dims(1) = size(array,1); dims(2) = size(array,2)
    dims(3) = size(array,3); dims(4) = size(array,4)
    call write_dataset_real32_rank(gid, 'data', array, dims, 4, bp, status)
  end subroutine hdf5_append_4D_real32

  subroutine hdf5_append_1D_real64(state, array, status, bitpix)
    type(hdf5_state_type), intent(inout) :: state
    real(real64),          intent(in)    :: array(:)
    integer,               intent(inout) :: status
    integer, optional,     intent(in)    :: bitpix
    integer(hsize_t) :: dims(1)
    integer(hid_t) :: gid
    integer :: bp
    real(real64) :: max_abs, min_nz
    bp = -64; if (present(bitpix)) bp = bitpix
    if (bp == 0) then
       max_abs = maxval(abs(array))
       if (max_abs > 0.0_real64) then
          min_nz = minval(abs(array), mask=abs(array) > 0.0_real64)
       else
          min_nz = 0.0_real64
       endif
       bp = pick_bitpix64(max_abs, min_nz)
    endif
    call begin_new_section(state, status, gid); if (status /= 0) return
    dims(1) = size(array)
    call write_dataset_real64_rank(gid, 'data', array, dims, 1, bp, status)
  end subroutine hdf5_append_1D_real64

  subroutine hdf5_append_2D_real64(state, array, status, bitpix)
    type(hdf5_state_type), intent(inout) :: state
    real(real64),          intent(in)    :: array(:,:)
    integer,               intent(inout) :: status
    integer, optional,     intent(in)    :: bitpix
    integer(hsize_t) :: dims(2)
    integer(hid_t) :: gid
    integer :: bp
    real(real64) :: max_abs, min_nz
    bp = -64; if (present(bitpix)) bp = bitpix
    if (bp == 0) then
       max_abs = maxval(abs(array))
       if (max_abs > 0.0_real64) then
          min_nz = minval(abs(array), mask=abs(array) > 0.0_real64)
       else
          min_nz = 0.0_real64
       endif
       bp = pick_bitpix64(max_abs, min_nz)
    endif
    call begin_new_section(state, status, gid); if (status /= 0) return
    dims(1) = size(array,1); dims(2) = size(array,2)
    call write_dataset_real64_rank(gid, 'data', array, dims, 2, bp, status)
  end subroutine hdf5_append_2D_real64

  subroutine hdf5_append_3D_real64(state, array, status, bitpix)
    type(hdf5_state_type), intent(inout) :: state
    real(real64),          intent(in)    :: array(:,:,:)
    integer,               intent(inout) :: status
    integer, optional,     intent(in)    :: bitpix
    integer(hsize_t) :: dims(3)
    integer(hid_t) :: gid
    integer :: bp
    real(real64) :: max_abs, min_nz
    bp = -64; if (present(bitpix)) bp = bitpix
    if (bp == 0) then
       max_abs = maxval(abs(array))
       if (max_abs > 0.0_real64) then
          min_nz = minval(abs(array), mask=abs(array) > 0.0_real64)
       else
          min_nz = 0.0_real64
       endif
       bp = pick_bitpix64(max_abs, min_nz)
    endif
    call begin_new_section(state, status, gid); if (status /= 0) return
    dims(1) = size(array,1); dims(2) = size(array,2); dims(3) = size(array,3)
    call write_dataset_real64_rank(gid, 'data', array, dims, 3, bp, status)
  end subroutine hdf5_append_3D_real64

  subroutine hdf5_append_4D_real64(state, array, status, bitpix)
    type(hdf5_state_type), intent(inout) :: state
    real(real64),          intent(in)    :: array(:,:,:,:)
    integer,               intent(inout) :: status
    integer, optional,     intent(in)    :: bitpix
    integer(hsize_t) :: dims(4)
    integer(hid_t) :: gid
    integer :: bp
    real(real64) :: max_abs, min_nz
    bp = -64; if (present(bitpix)) bp = bitpix
    if (bp == 0) then
       max_abs = maxval(abs(array))
       if (max_abs > 0.0_real64) then
          min_nz = minval(abs(array), mask=abs(array) > 0.0_real64)
       else
          min_nz = 0.0_real64
       endif
       bp = pick_bitpix64(max_abs, min_nz)
    endif
    call begin_new_section(state, status, gid); if (status /= 0) return
    dims(1) = size(array,1); dims(2) = size(array,2)
    dims(3) = size(array,3); dims(4) = size(array,4)
    call write_dataset_real64_rank(gid, 'data', array, dims, 4, bp, status)
  end subroutine hdf5_append_4D_real64

  subroutine hdf5_append_1D_int32(state, array, status)
    type(hdf5_state_type), intent(inout) :: state
    integer(int32),        intent(in)    :: array(:)
    integer,               intent(inout) :: status
    integer(hsize_t) :: dims(1)
    integer(hid_t) :: gid
    call begin_new_section(state, status, gid); if (status /= 0) return
    dims(1) = size(array)
    call write_dataset_int32_rank(gid, 'data', array, dims, 1, status)
  end subroutine hdf5_append_1D_int32

  subroutine hdf5_append_2D_int32(state, array, status)
    type(hdf5_state_type), intent(inout) :: state
    integer(int32),        intent(in)    :: array(:,:)
    integer,               intent(inout) :: status
    integer(hsize_t) :: dims(2)
    integer(hid_t) :: gid
    call begin_new_section(state, status, gid); if (status /= 0) return
    dims(1) = size(array,1); dims(2) = size(array,2)
    call write_dataset_int32_rank(gid, 'data', array, dims, 2, status)
  end subroutine hdf5_append_2D_int32

  subroutine hdf5_append_3D_int32(state, array, status)
    type(hdf5_state_type), intent(inout) :: state
    integer(int32),        intent(in)    :: array(:,:,:)
    integer,               intent(inout) :: status
    integer(hsize_t) :: dims(3)
    integer(hid_t) :: gid
    call begin_new_section(state, status, gid); if (status /= 0) return
    dims(1) = size(array,1); dims(2) = size(array,2); dims(3) = size(array,3)
    call write_dataset_int32_rank(gid, 'data', array, dims, 3, status)
  end subroutine hdf5_append_3D_int32

  subroutine hdf5_append_4D_int32(state, array, status)
    type(hdf5_state_type), intent(inout) :: state
    integer(int32),        intent(in)    :: array(:,:,:,:)
    integer,               intent(inout) :: status
    integer(hsize_t) :: dims(4)
    integer(hid_t) :: gid
    call begin_new_section(state, status, gid); if (status /= 0) return
    dims(1) = size(array,1); dims(2) = size(array,2)
    dims(3) = size(array,3); dims(4) = size(array,4)
    call write_dataset_int32_rank(gid, 'data', array, dims, 4, status)
  end subroutine hdf5_append_4D_int32

!=============================================================================
!  Image read (reads /<cursection>/data)
!=============================================================================
  subroutine read_dataset_real32_rank(state, array, dims, ndim, status)
    type(hdf5_state_type), intent(in)    :: state
    real(real32),          intent(out)   :: array(*)
    integer(hsize_t),      intent(in)    :: dims(:)
    integer,               intent(in)    :: ndim
    integer,               intent(inout) :: status
    integer(hid_t) :: dset_id
    integer :: ierr
    call h5dopen_f(int(state%cur_group, hid_t), 'data', dset_id, ierr)
    if (ierr /= 0) then; status = HDF5_ERR_NOT_FOUND; return; endif
    call h5dread_f(dset_id, H5T_NATIVE_REAL, array, dims(1:ndim), ierr)
    call h5dclose_f(dset_id, ierr)
    if (ierr /= 0) status = HDF5_ERR_BACKEND
  end subroutine read_dataset_real32_rank

  subroutine read_dataset_real64_rank(state, array, dims, ndim, status)
    type(hdf5_state_type), intent(in)    :: state
    real(real64),          intent(out)   :: array(*)
    integer(hsize_t),      intent(in)    :: dims(:)
    integer,               intent(in)    :: ndim
    integer,               intent(inout) :: status
    integer(hid_t) :: dset_id
    integer :: ierr
    call h5dopen_f(int(state%cur_group, hid_t), 'data', dset_id, ierr)
    if (ierr /= 0) then; status = HDF5_ERR_NOT_FOUND; return; endif
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dims(1:ndim), ierr)
    call h5dclose_f(dset_id, ierr)
    if (ierr /= 0) status = HDF5_ERR_BACKEND
  end subroutine read_dataset_real64_rank

  subroutine hdf5_read_1D_real32(state, array, status)
    type(hdf5_state_type), intent(in)    :: state
    real(real32),          intent(out)   :: array(:)
    integer,               intent(inout) :: status
    integer(hsize_t) :: dims(1)
    dims(1) = size(array)
    call read_dataset_real32_rank(state, array, dims, 1, status)
  end subroutine hdf5_read_1D_real32

  subroutine hdf5_read_2D_real32(state, array, status)
    type(hdf5_state_type), intent(in)    :: state
    real(real32),          intent(out)   :: array(:,:)
    integer,               intent(inout) :: status
    integer(hsize_t) :: dims(2)
    dims(1) = size(array,1); dims(2) = size(array,2)
    call read_dataset_real32_rank(state, array, dims, 2, status)
  end subroutine hdf5_read_2D_real32

  subroutine hdf5_read_3D_real32(state, array, status)
    type(hdf5_state_type), intent(in)    :: state
    real(real32),          intent(out)   :: array(:,:,:)
    integer,               intent(inout) :: status
    integer(hsize_t) :: dims(3)
    dims(1) = size(array,1); dims(2) = size(array,2); dims(3) = size(array,3)
    call read_dataset_real32_rank(state, array, dims, 3, status)
  end subroutine hdf5_read_3D_real32

  subroutine hdf5_read_4D_real32(state, array, status)
    type(hdf5_state_type), intent(in)    :: state
    real(real32),          intent(out)   :: array(:,:,:,:)
    integer,               intent(inout) :: status
    integer(hsize_t) :: dims(4)
    dims(1) = size(array,1); dims(2) = size(array,2)
    dims(3) = size(array,3); dims(4) = size(array,4)
    call read_dataset_real32_rank(state, array, dims, 4, status)
  end subroutine hdf5_read_4D_real32

  subroutine hdf5_read_1D_real64(state, array, status)
    type(hdf5_state_type), intent(in)    :: state
    real(real64),          intent(out)   :: array(:)
    integer,               intent(inout) :: status
    integer(hsize_t) :: dims(1)
    dims(1) = size(array)
    call read_dataset_real64_rank(state, array, dims, 1, status)
  end subroutine hdf5_read_1D_real64

  subroutine hdf5_read_2D_real64(state, array, status)
    type(hdf5_state_type), intent(in)    :: state
    real(real64),          intent(out)   :: array(:,:)
    integer,               intent(inout) :: status
    integer(hsize_t) :: dims(2)
    dims(1) = size(array,1); dims(2) = size(array,2)
    call read_dataset_real64_rank(state, array, dims, 2, status)
  end subroutine hdf5_read_2D_real64

  subroutine hdf5_read_3D_real64(state, array, status)
    type(hdf5_state_type), intent(in)    :: state
    real(real64),          intent(out)   :: array(:,:,:)
    integer,               intent(inout) :: status
    integer(hsize_t) :: dims(3)
    dims(1) = size(array,1); dims(2) = size(array,2); dims(3) = size(array,3)
    call read_dataset_real64_rank(state, array, dims, 3, status)
  end subroutine hdf5_read_3D_real64

  subroutine hdf5_read_4D_real64(state, array, status)
    type(hdf5_state_type), intent(in)    :: state
    real(real64),          intent(out)   :: array(:,:,:,:)
    integer,               intent(inout) :: status
    integer(hsize_t) :: dims(4)
    dims(1) = size(array,1); dims(2) = size(array,2)
    dims(3) = size(array,3); dims(4) = size(array,4)
    call read_dataset_real64_rank(state, array, dims, 4, status)
  end subroutine hdf5_read_4D_real64

!=============================================================================
!  Table column append (multiple 1-D datasets under one section group)
!=============================================================================
  subroutine register_column(state, colname)
    type(hdf5_state_type), intent(inout) :: state
    character(len=*),      intent(in)    :: colname
    integer :: n
    if (state%cur_idx < 1) return
    n = state%ncols(state%cur_idx) + 1
    if (n > HDF5_MAX_COLS) return
    state%col_names(n, state%cur_idx) = trim(colname)
    state%ncols(state%cur_idx) = n
  end subroutine register_column

  subroutine hdf5_append_table_column_real32(state, colname, array, status, bitpix)
    type(hdf5_state_type), intent(inout) :: state
    character(len=*),      intent(in)    :: colname
    real(real32),          intent(in)    :: array(:)
    integer,               intent(inout) :: status
    integer, optional,     intent(in)    :: bitpix
    integer(hsize_t) :: dims(1)
    integer :: bp
    real(real32) :: max_abs, min_nz
    call ensure_table_section(state, status); if (status /= 0) return
    bp = -32; if (present(bitpix)) bp = bitpix
    if (bp == 0) then
       max_abs = maxval(abs(array))
       if (max_abs > 0.0_real32) then
          min_nz = minval(abs(array), mask=abs(array) > 0.0_real32)
       else
          min_nz = 0.0_real32
       endif
       bp = pick_bitpix32(max_abs, min_nz)
    endif
    dims(1) = size(array)
    call write_dataset_real32_rank(int(state%cur_group, hid_t), trim(colname), &
                                   array, dims, 1, bp, status)
    if (status == 0) call register_column(state, colname)
  end subroutine hdf5_append_table_column_real32

  subroutine hdf5_append_table_column_real64(state, colname, array, status, bitpix)
    type(hdf5_state_type), intent(inout) :: state
    character(len=*),      intent(in)    :: colname
    real(real64),          intent(in)    :: array(:)
    integer,               intent(inout) :: status
    integer, optional,     intent(in)    :: bitpix
    integer(hsize_t) :: dims(1)
    integer :: bp
    real(real64) :: max_abs, min_nz
    call ensure_table_section(state, status); if (status /= 0) return
    bp = -64; if (present(bitpix)) bp = bitpix
    if (bp == 0) then
       max_abs = maxval(abs(array))
       if (max_abs > 0.0_real64) then
          min_nz = minval(abs(array), mask=abs(array) > 0.0_real64)
       else
          min_nz = 0.0_real64
       endif
       bp = pick_bitpix64(max_abs, min_nz)
    endif
    dims(1) = size(array)
    call write_dataset_real64_rank(int(state%cur_group, hid_t), trim(colname), &
                                   array, dims, 1, bp, status)
    if (status == 0) call register_column(state, colname)
  end subroutine hdf5_append_table_column_real64

  subroutine hdf5_append_table_column_int32(state, colname, array, status)
    type(hdf5_state_type), intent(inout) :: state
    character(len=*),      intent(in)    :: colname
    integer(int32),        intent(in)    :: array(:)
    integer,               intent(inout) :: status
    integer(hsize_t) :: dims(1)
    call ensure_table_section(state, status); if (status /= 0) return
    dims(1) = size(array)
    call write_dataset_int32_rank(int(state%cur_group, hid_t), trim(colname), &
                                  array, dims, 1, status)
    if (status == 0) call register_column(state, colname)
  end subroutine hdf5_append_table_column_int32

  subroutine hdf5_append_table_column_int64(state, colname, array, status)
    type(hdf5_state_type), intent(inout) :: state
    character(len=*),      intent(in)    :: colname
    integer(int64),        intent(in)    :: array(:)
    integer,               intent(inout) :: status
    integer(hsize_t) :: dims(1)
    integer(hid_t) :: space_id, dset_id
    integer :: ierr
    call ensure_table_section(state, status); if (status /= 0) return
    dims(1) = size(array)
    call h5screate_simple_f(1, dims, space_id, ierr)
    call h5dcreate_f(int(state%cur_group, hid_t), trim(colname), &
                     H5T_STD_I64LE, space_id, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_STD_I64LE, array, dims, ierr)
    call h5dclose_f(dset_id, ierr)
    call h5sclose_f(space_id, ierr)
    if (ierr /= 0) then; status = HDF5_ERR_BACKEND; return; endif
    call register_column(state, colname)
  end subroutine hdf5_append_table_column_int64

!=============================================================================
!  Table column read (by 1-based column index in current section)
!=============================================================================
  subroutine hdf5_read_table_column_real32(state, colnum, array, status)
    type(hdf5_state_type), intent(in)    :: state
    integer,               intent(in)    :: colnum
    real(real32),          intent(out)   :: array(:)
    integer,               intent(inout) :: status
    integer(hsize_t) :: dims(1)
    integer(hid_t) :: dset_id
    integer :: ierr
    character(len=HDF5_NAME_LEN) :: nm
    if (colnum < 1 .or. colnum > state%ncols(state%cur_idx)) then
       status = HDF5_ERR_NOT_FOUND; return
    endif
    nm = state%col_names(colnum, state%cur_idx)
    dims(1) = size(array)
    call h5dopen_f(int(state%cur_group, hid_t), trim(nm), dset_id, ierr)
    if (ierr /= 0) then; status = HDF5_ERR_NOT_FOUND; return; endif
    call h5dread_f(dset_id, H5T_NATIVE_REAL, array, dims, ierr)
    call h5dclose_f(dset_id, ierr)
    if (ierr /= 0) status = HDF5_ERR_BACKEND
  end subroutine hdf5_read_table_column_real32

  subroutine hdf5_read_table_column_real64(state, colnum, array, status)
    type(hdf5_state_type), intent(in)    :: state
    integer,               intent(in)    :: colnum
    real(real64),          intent(out)   :: array(:)
    integer,               intent(inout) :: status
    integer(hsize_t) :: dims(1)
    integer(hid_t) :: dset_id
    integer :: ierr
    character(len=HDF5_NAME_LEN) :: nm
    if (colnum < 1 .or. colnum > state%ncols(state%cur_idx)) then
       status = HDF5_ERR_NOT_FOUND; return
    endif
    nm = state%col_names(colnum, state%cur_idx)
    dims(1) = size(array)
    call h5dopen_f(int(state%cur_group, hid_t), trim(nm), dset_id, ierr)
    if (ierr /= 0) then; status = HDF5_ERR_NOT_FOUND; return; endif
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, ierr)
    call h5dclose_f(dset_id, ierr)
    if (ierr /= 0) status = HDF5_ERR_BACKEND
  end subroutine hdf5_read_table_column_real64

  subroutine hdf5_read_table_column_int32(state, colnum, array, status)
    type(hdf5_state_type), intent(in)    :: state
    integer,               intent(in)    :: colnum
    integer(int32),        intent(out)   :: array(:)
    integer,               intent(inout) :: status
    integer(hsize_t) :: dims(1)
    integer(hid_t) :: dset_id
    integer :: ierr
    character(len=HDF5_NAME_LEN) :: nm
    if (colnum < 1 .or. colnum > state%ncols(state%cur_idx)) then
       status = HDF5_ERR_NOT_FOUND; return
    endif
    nm = state%col_names(colnum, state%cur_idx)
    dims(1) = size(array)
    call h5dopen_f(int(state%cur_group, hid_t), trim(nm), dset_id, ierr)
    if (ierr /= 0) then; status = HDF5_ERR_NOT_FOUND; return; endif
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, array, dims, ierr)
    call h5dclose_f(dset_id, ierr)
    if (ierr /= 0) status = HDF5_ERR_BACKEND
  end subroutine hdf5_read_table_column_int32

  subroutine hdf5_read_table_column_int64(state, colnum, array, status)
    type(hdf5_state_type), intent(in)    :: state
    integer,               intent(in)    :: colnum
    integer(int64),        intent(out)   :: array(:)
    integer,               intent(inout) :: status
    integer(hsize_t) :: dims(1)
    integer(hid_t) :: dset_id
    integer :: ierr
    character(len=HDF5_NAME_LEN) :: nm
    if (colnum < 1 .or. colnum > state%ncols(state%cur_idx)) then
       status = HDF5_ERR_NOT_FOUND; return
    endif
    nm = state%col_names(colnum, state%cur_idx)
    dims(1) = size(array)
    call h5dopen_f(int(state%cur_group, hid_t), trim(nm), dset_id, ierr)
    if (ierr /= 0) then; status = HDF5_ERR_NOT_FOUND; return; endif
    call h5dread_f(dset_id, H5T_STD_I64LE, array, dims, ierr)
    call h5dclose_f(dset_id, ierr)
    if (ierr /= 0) status = HDF5_ERR_BACKEND
  end subroutine hdf5_read_table_column_int64

#else
!=============================================================================
!  Stubs — compiled when HDF5 support is OFF.  All calls return NOT_BUILT.
!=============================================================================
  subroutine hdf5_init(status);     integer, intent(inout) :: status; status = HDF5_ERR_NOT_BUILT; end subroutine
  subroutine hdf5_finalize(status); integer, intent(inout) :: status; status = HDF5_ERR_NOT_BUILT; end subroutine

  subroutine hdf5_open_new(state, fname, status)
    type(hdf5_state_type), intent(out) :: state
    character(len=*),      intent(in)  :: fname
    integer,             intent(inout) :: status
    status = HDF5_ERR_NOT_BUILT
  end subroutine
  subroutine hdf5_open_old(state, fname, status, writable)
    type(hdf5_state_type), intent(out) :: state
    character(len=*),      intent(in)  :: fname
    integer,             intent(inout) :: status
    logical, optional,     intent(in)  :: writable
    status = HDF5_ERR_NOT_BUILT
  end subroutine
  subroutine hdf5_close(state, status)
    type(hdf5_state_type), intent(inout) :: state
    integer,               intent(inout) :: status
    status = HDF5_ERR_NOT_BUILT
  end subroutine
  subroutine hdf5_move_to_next_section(state, status, nstep)
    type(hdf5_state_type), intent(inout) :: state
    integer,               intent(inout) :: status
    integer, optional,     intent(in)    :: nstep
    status = HDF5_ERR_NOT_BUILT
  end subroutine
  subroutine hdf5_get_column_number(state, colname, colnum, status)
    type(hdf5_state_type), intent(in)    :: state
    character(len=*),      intent(in)    :: colname
    integer,               intent(out)   :: colnum
    integer,               intent(inout) :: status
    colnum = 0; status = HDF5_ERR_NOT_BUILT
  end subroutine

#define STUB_PUT_KEY(NAME, TYPE) \
  subroutine NAME(state, name, value, comment, status); \
    type(hdf5_state_type), intent(in)    :: state; \
    character(len=*),      intent(in)    :: name; \
    TYPE,                  intent(in)    :: value; \
    character(len=*),      intent(in)    :: comment; \
    integer,               intent(inout) :: status; \
    status = HDF5_ERR_NOT_BUILT; \
  end subroutine

  STUB_PUT_KEY(hdf5_put_key_logical, logical)
  STUB_PUT_KEY(hdf5_put_key_real32,  real(real32))
  STUB_PUT_KEY(hdf5_put_key_real64,  real(real64))
  STUB_PUT_KEY(hdf5_put_key_int32,   integer(int32))
  STUB_PUT_KEY(hdf5_put_key_int64,   integer(int64))

  subroutine hdf5_put_key_string(state, name, value, comment, status)
    type(hdf5_state_type), intent(inout) :: state
    character(len=*),      intent(in)    :: name, value, comment
    integer,               intent(inout) :: status
    status = HDF5_ERR_NOT_BUILT
  end subroutine

#define STUB_GET_KEY(NAME, TYPE, INIT) \
  subroutine NAME(state, name, value, status, comment); \
    type(hdf5_state_type),      intent(in)    :: state; \
    character(len=*),           intent(in)    :: name; \
    TYPE,                       intent(out)   :: value; \
    integer,                    intent(inout) :: status; \
    character(len=*), optional, intent(out)   :: comment; \
    value = INIT; status = HDF5_ERR_NOT_BUILT; \
    if (present(comment)) comment = ''; \
  end subroutine

  STUB_GET_KEY(hdf5_get_key_logical, logical,        .false.)
  STUB_GET_KEY(hdf5_get_key_real32,  real(real32),   0.0_real32)
  STUB_GET_KEY(hdf5_get_key_real64,  real(real64),   0.0_real64)
  STUB_GET_KEY(hdf5_get_key_int32,   integer(int32), 0_int32)
  STUB_GET_KEY(hdf5_get_key_int64,   integer(int64), 0_int64)

  subroutine hdf5_get_key_string(state, name, value, status, comment)
    type(hdf5_state_type),      intent(in)    :: state
    character(len=*),           intent(in)    :: name
    character(len=*),           intent(out)   :: value
    integer,                    intent(inout) :: status
    character(len=*), optional, intent(out)   :: comment
    value = ''; status = HDF5_ERR_NOT_BUILT
    if (present(comment)) comment = ''
  end subroutine

#define STUB_RW_IMAGE(NAME, RANKDECL, TYPE) \
  subroutine NAME(state, array, status, bitpix); \
    type(hdf5_state_type), intent(inout) :: state; \
    TYPE,                  intent(in)    :: array RANKDECL; \
    integer,               intent(inout) :: status; \
    integer, optional,     intent(in)    :: bitpix; \
    status = HDF5_ERR_NOT_BUILT; \
  end subroutine

  STUB_RW_IMAGE(hdf5_append_1D_real32, (:),       real(real32))
  STUB_RW_IMAGE(hdf5_append_2D_real32, (:,:),     real(real32))
  STUB_RW_IMAGE(hdf5_append_3D_real32, (:,:,:),   real(real32))
  STUB_RW_IMAGE(hdf5_append_4D_real32, (:,:,:,:), real(real32))
  STUB_RW_IMAGE(hdf5_append_1D_real64, (:),       real(real64))
  STUB_RW_IMAGE(hdf5_append_2D_real64, (:,:),     real(real64))
  STUB_RW_IMAGE(hdf5_append_3D_real64, (:,:,:),   real(real64))
  STUB_RW_IMAGE(hdf5_append_4D_real64, (:,:,:,:), real(real64))

#define STUB_IMG_INT(NAME, RANKDECL) \
  subroutine NAME(state, array, status); \
    type(hdf5_state_type), intent(inout) :: state; \
    integer(int32),        intent(in)    :: array RANKDECL; \
    integer,               intent(inout) :: status; \
    status = HDF5_ERR_NOT_BUILT; \
  end subroutine

  STUB_IMG_INT(hdf5_append_1D_int32, (:))
  STUB_IMG_INT(hdf5_append_2D_int32, (:,:))
  STUB_IMG_INT(hdf5_append_3D_int32, (:,:,:))
  STUB_IMG_INT(hdf5_append_4D_int32, (:,:,:,:))

#define STUB_READ_IMAGE(NAME, RANKDECL, TYPE, INIT) \
  subroutine NAME(state, array, status); \
    type(hdf5_state_type), intent(in)    :: state; \
    TYPE,                  intent(out)   :: array RANKDECL; \
    integer,               intent(inout) :: status; \
    array = INIT; status = HDF5_ERR_NOT_BUILT; \
  end subroutine

  STUB_READ_IMAGE(hdf5_read_1D_real32, (:),       real(real32), 0.0_real32)
  STUB_READ_IMAGE(hdf5_read_2D_real32, (:,:),     real(real32), 0.0_real32)
  STUB_READ_IMAGE(hdf5_read_3D_real32, (:,:,:),   real(real32), 0.0_real32)
  STUB_READ_IMAGE(hdf5_read_4D_real32, (:,:,:,:), real(real32), 0.0_real32)
  STUB_READ_IMAGE(hdf5_read_1D_real64, (:),       real(real64), 0.0_real64)
  STUB_READ_IMAGE(hdf5_read_2D_real64, (:,:),     real(real64), 0.0_real64)
  STUB_READ_IMAGE(hdf5_read_3D_real64, (:,:,:),   real(real64), 0.0_real64)
  STUB_READ_IMAGE(hdf5_read_4D_real64, (:,:,:,:), real(real64), 0.0_real64)

#define STUB_APPEND_TC(NAME, TYPE) \
  subroutine NAME(state, colname, array, status, bitpix); \
    type(hdf5_state_type), intent(inout) :: state; \
    character(len=*),      intent(in)    :: colname; \
    TYPE,                  intent(in)    :: array(:); \
    integer,               intent(inout) :: status; \
    integer, optional,     intent(in)    :: bitpix; \
    status = HDF5_ERR_NOT_BUILT; \
  end subroutine

  STUB_APPEND_TC(hdf5_append_table_column_real32, real(real32))
  STUB_APPEND_TC(hdf5_append_table_column_real64, real(real64))

  subroutine hdf5_append_table_column_int32(state, colname, array, status)
    type(hdf5_state_type), intent(inout) :: state
    character(len=*),      intent(in)    :: colname
    integer(int32),        intent(in)    :: array(:)
    integer,               intent(inout) :: status
    status = HDF5_ERR_NOT_BUILT
  end subroutine

  subroutine hdf5_append_table_column_int64(state, colname, array, status)
    type(hdf5_state_type), intent(inout) :: state
    character(len=*),      intent(in)    :: colname
    integer(int64),        intent(in)    :: array(:)
    integer,               intent(inout) :: status
    status = HDF5_ERR_NOT_BUILT
  end subroutine

#define STUB_READ_TC(NAME, TYPE, INIT) \
  subroutine NAME(state, colnum, array, status); \
    type(hdf5_state_type), intent(in)    :: state; \
    integer,               intent(in)    :: colnum; \
    TYPE,                  intent(out)   :: array(:); \
    integer,               intent(inout) :: status; \
    array = INIT; status = HDF5_ERR_NOT_BUILT; \
  end subroutine

  STUB_READ_TC(hdf5_read_table_column_real32, real(real32),   0.0_real32)
  STUB_READ_TC(hdf5_read_table_column_real64, real(real64),   0.0_real64)
  STUB_READ_TC(hdf5_read_table_column_int32,  integer(int32), 0_int32)
  STUB_READ_TC(hdf5_read_table_column_int64,  integer(int64), 0_int64)
#endif

end module hdf5io_mod
