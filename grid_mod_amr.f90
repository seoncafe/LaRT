module grid_mod_amr
  !-----------------------------------------------------------------------
  ! AMR grid creation: reads data (RAMSES or generic format), builds the
  ! octree, and normalizes physical quantities to LaRT's internal units.
  !
  ! Called from setup_v2c.f90 when par%use_amr_grid = .true.
  !
  ! After grid_create_amr returns:
  !   amr_grid%rhokap(il)   = HI opacity per unit length at line center
  !   amr_grid%rhokapD(il)  = dust opacity per unit length  (if DGR > 0)
  !   amr_grid%Dfreq(il)    = local Doppler frequency
  !   amr_grid%voigt_a(il)  = Voigt damping parameter
  !   amr_grid%vfx/vfy/vfz  = velocity / v_thermal (dimensionless)
  !   amr_grid%xfreq(*)     = frequency array
  !   amr_grid%Jout(*), Jin(*), Jabs(*) allocated
  !-----------------------------------------------------------------------
  use octree_mod
  use read_generic_amr_mod
  use physics_amr_mod
  use ion_data_mod
  use read_text_data
  use random
  use voigt_mod
  use utility
  implicit none

  public :: grid_create_amr, grid_destroy_amr, amr_sync_to_grid

  private :: assign_amr_velocities_from_type

contains

  !=========================================================================
  subroutine grid_create_amr(grid)
    use define
    use mpi
    use line_mod
    implicit none
    type(grid_type), intent(inout) :: grid   ! Cartesian grid (kept for observer/output compat.)

    ! leaf-cell arrays
    real(wp), allocatable :: xleaf(:), yleaf(:), zleaf(:)
    real(wp), allocatable :: nH_cgs(:), T_cgs(:)
    real(wp), allocatable :: vx_kms(:), vy_kms(:), vz_kms(:)
    integer,  allocatable :: leaf_level(:)
    integer  :: nleaf
    real(wp) :: boxlen_code   ! box length in code units (kpc, pc, cm, etc.)
    real(wp) :: box_origin_x, box_origin_y, box_origin_z   ! lower-corner origin

    ! Optional extended arrays (allocated only if present in generic file)
    real(wp), allocatable :: Z_arr(:)       ! metallicity
    real(wp), allocatable :: xHI_arr(:)     ! neutral fraction from file
    real(wp), allocatable :: ne_arr(:)      ! electron density from file
    real(wp), allocatable :: nion_arr(:)    ! ion density from file
    real(wp), allocatable :: emiss_arr(:)   ! emissivity from file
    real(wp), allocatable :: ndust_arr(:)   ! dust density from file
    logical :: have_Z, have_xHI, have_ne, have_nion, have_emiss, have_ndust

    ! Working arrays
    real(wp), allocatable :: nHI_frac(:)   ! nHI/nH neutral fraction
    integer  :: il, icell
    real(wp) :: vtherm
    real(wp) :: nH, opacity_sum, nopac, opac_norm, opac_length, cos_cone, rr
    real(wp) :: xscale, atau0, atau0_cell, xi, chi
    real(wp) :: taupole, N_HIpole, tauhomo, N_HIhomo, NHI_pole_raw
    real(wp) :: atau1
    integer  :: ix, ierr

    opac_norm = 1.0_wp
    !--- Step 1: Read leaf-cell data ----------------------------------
    if (mpar%p_rank == 0) then
      ! Generic AMR reader returns code units (e.g., kpc).
      ! par%distance2cm (set by distance_unit in setup_v2c) converts to cm.
      ! Optional columns are allocated only if present in the file.
      call generic_amr_read(trim(par%amr_file), &
          xleaf, yleaf, zleaf, leaf_level, &
          nH_cgs, T_cgs, vx_kms, vy_kms, vz_kms, &
          nleaf, boxlen_code, &
          metallicity=Z_arr, xHI=xHI_arr, n_e=ne_arr, &
          n_ion=nion_arr, emissivity=emiss_arr, ndust=ndust_arr, &
          origin_x=box_origin_x, origin_y=box_origin_y, origin_z=box_origin_z)
      write(6,'(a,i0)')    'AMR nleaf    : ', nleaf
      write(6,'(a,es12.4,a)') 'AMR boxlen   : ', boxlen_code, ' (code units)'
      write(6,'(a,es12.4,a)') 'AMR boxlen   : ', boxlen_code*par%distance2cm, ' cm'
      if (allocated(Z_arr))     write(6,'(a)') 'AMR optional : metallicity column found'
      if (allocated(xHI_arr))   write(6,'(a)') 'AMR optional : xHI column found'
      if (allocated(ne_arr))    write(6,'(a)') 'AMR optional : n_e column found'
      if (allocated(nion_arr))  write(6,'(a)') 'AMR optional : n_ion column found'
      if (allocated(emiss_arr)) write(6,'(a)') 'AMR optional : emissivity column found'
      if (allocated(ndust_arr)) write(6,'(a)') 'AMR optional : ndust column found'
    end if

    ! Broadcast nleaf, boxlen_code, and origin to all ranks
    call MPI_BCAST(nleaf,        1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(boxlen_code,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(box_origin_x, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(box_origin_y, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(box_origin_z, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    if (mpar%p_rank /= 0) then
      allocate(xleaf(nleaf), yleaf(nleaf), zleaf(nleaf))
      allocate(leaf_level(nleaf))
      allocate(nH_cgs(nleaf), T_cgs(nleaf))
      allocate(vx_kms(nleaf), vy_kms(nleaf), vz_kms(nleaf))
    end if
    call MPI_BCAST(xleaf(1),      nleaf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(yleaf(1),      nleaf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(zleaf(1),      nleaf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(leaf_level(1), nleaf, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nH_cgs(1),     nleaf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(T_cgs(1),      nleaf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(vx_kms(1),     nleaf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(vy_kms(1),     nleaf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(vz_kms(1),     nleaf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Broadcast optional column presence flags and data
    have_Z     = allocated(Z_arr)
    have_xHI   = allocated(xHI_arr)
    have_ne    = allocated(ne_arr)
    have_nion  = allocated(nion_arr)
    have_emiss = allocated(emiss_arr)
    have_ndust = allocated(ndust_arr)
    call MPI_BCAST(have_Z,     1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(have_xHI,   1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(have_ne,    1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(have_nion,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(have_emiss, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(have_ndust, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

    if (mpar%p_rank /= 0) then
      if (have_Z)     allocate(Z_arr(nleaf))
      if (have_xHI)   allocate(xHI_arr(nleaf))
      if (have_ne)    allocate(ne_arr(nleaf))
      if (have_nion)  allocate(nion_arr(nleaf))
      if (have_emiss) allocate(emiss_arr(nleaf))
      if (have_ndust) allocate(ndust_arr(nleaf))
    end if
    if (have_Z)     call MPI_BCAST(Z_arr(1),     nleaf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (have_xHI)   call MPI_BCAST(xHI_arr(1),   nleaf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (have_ne)    call MPI_BCAST(ne_arr(1),     nleaf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (have_nion)  call MPI_BCAST(nion_arr(1),   nleaf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (have_emiss) call MPI_BCAST(emiss_arr(1),  nleaf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (have_ndust) call MPI_BCAST(ndust_arr(1),  nleaf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    !--- Step 2: Build octree in code units ---------------------------
    ! Box extent is [origin, origin+boxlen] in the same frame as the input
    ! cell positions.  For data with default centered origin (-boxlen/2),
    ! box runs [-boxlen/2, +boxlen/2].  For data with corner-based origin
    ! (ORIGINX=0), box runs [0, boxlen].
    if (mpar%p_rank == 0) then
      write(6,'(a,3es12.4)') 'AMR box origin: ', box_origin_x, box_origin_y, box_origin_z
      write(6,'(a,3es12.4)') 'AMR box extent: ', box_origin_x + boxlen_code, &
                                                  box_origin_y + boxlen_code, &
                                                  box_origin_z + boxlen_code
      write(6,'(a,3es12.4)') 'AMR data x rng: ', minval(xleaf), maxval(xleaf)
      write(6,'(a,3es12.4)') 'AMR data y rng: ', minval(yleaf), maxval(yleaf)
      write(6,'(a,3es12.4)') 'AMR data z rng: ', minval(zleaf), maxval(zleaf)
      block
        real(wp) :: cx_box, cy_box, cz_box
        cx_box = box_origin_x + 0.5_wp * boxlen_code
        cy_box = box_origin_y + 0.5_wp * boxlen_code
        cz_box = box_origin_z + 0.5_wp * boxlen_code
        if (minval(xleaf) > cx_box .or. maxval(xleaf) < cx_box .or. &
            minval(yleaf) > cy_box .or. maxval(yleaf) < cy_box .or. &
            minval(zleaf) > cz_box .or. maxval(zleaf) < cz_box) then
          write(6,'(a)') '------------------------------------------------------------'
          write(6,'(a)') 'WARNING: AMR data does not cover the box center.'
          write(6,'(a,3es12.4)') '  box center:  ', cx_box, cy_box, cz_box
          write(6,'(a)') '  pole traversal (from box center +z) will return N(HI)_pole=0,'
          write(6,'(a)') '  and tau normalization will fall back to tauhomo.'
          write(6,'(a)') '------------------------------------------------------------'
        end if
      end block
    end if
    call amr_build_tree(xleaf, yleaf, zleaf, leaf_level, nleaf, &
        box_origin_x, box_origin_x + boxlen_code, &
        box_origin_y, box_origin_y + boxlen_code, &
        box_origin_z, box_origin_z + boxlen_code)

    ! Precompute face-neighbor table for O(1) cell-boundary crossings.
    call amr_build_neighbors

    !--- Set par box dimensions from AMR grid (overrides any user input).
    !    sphere  : rmax = xmax = ymax = zmax = boxlen/2
    !    cylinder: rmax = xmax = ymax = boxlen/2; zmax = boxlen/2 (cubic box)
    !    rectangle: rmax = -1 (undefined); xmax/ymax/zmax = boxlen/2
    select case(trim(par%geometry))
    case ('sphere')
       par%rmax = amr_grid%L_box * 0.5_wp
       par%xmax = par%rmax;  par%ymax = par%rmax;  par%zmax = par%rmax
    case ('cylinder')
       par%rmax = amr_grid%L_box * 0.5_wp
       par%xmax = par%rmax;  par%ymax = par%rmax
       par%zmax = amr_grid%L_box * 0.5_wp
    case default   ! 'rectangle'
       par%rmax = -1.0_wp
       par%xmax = amr_grid%xmax
       par%ymax = amr_grid%ymax
       par%zmax = amr_grid%zmax
    end select

    !--- Step 3: Set reference Doppler frequency ----------------------
    amr_grid%Dfreq_ref  = line%vtherm1 * sqrt(par%temperature) / (line%wavelength0 * um2km)
    amr_grid%voigt_amean = (line%damping / fourpi) / amr_grid%Dfreq_ref
    amr_grid%Dfreq_mean  = amr_grid%Dfreq_ref

    !--- Step 4: Allocate physical arrays and fill --------------------
    allocate(nHI_frac(nleaf))
    ! Allocate rhokapD whenever dust is in play: global DGR, ndust column
    ! in the file, or laursen09 dust_model with a metallicity source.
    call amr_alloc_phys( &
        par%DGR > 0.0_wp .or. &
        have_ndust .or. &
        (trim(par%dust_model) == 'laursen09' .and. &
         (have_Z .or. par%metallicity_global >= 0.0_wp)))

    ! Shared arrays: only rank 0 on each node fills; others wait at barrier.
    if (mpar%h_rank == 0) then
    do il = 1, nleaf
      T_cgs(il) = max(T_cgs(il), 10.0_wp)
      vtherm    = line%vtherm1 * sqrt(T_cgs(il))
      amr_grid%Dfreq(il)   = vtherm / (line%wavelength0 * um2km)
      amr_grid%voigt_a(il) = (line%damping / fourpi) / amr_grid%Dfreq(il)

      ! --- Neutral fraction: file column > ionization_model > use_cie_condition ---
      if (have_xHI) then
        nHI_frac(il) = xHI_arr(il)
      else
        select case (trim(par%ionization_model))
        case ('cie_formula')
          if (par%use_cie_condition) then
            nHI_frac(il) = cie_neutral_fraction_formula(T_cgs(il))
          else
            nHI_frac(il) = 1.0_wp
          end if
        case ('cie_table')
          nHI_frac(il) = cie_neutral_fraction_table(nH_cgs(il), T_cgs(il))
        case ('full_neutral')
          nHI_frac(il) = 1.0_wp
        case ('from_file')
          write(6,'(a)') 'ERROR: ionization_model=from_file but xHI column not found in file.'
          stop 'ionization_model=from_file requires xHI column'
        case default
          if (par%use_cie_condition) then
            nHI_frac(il) = cie_neutral_fraction_formula(T_cgs(il))
          else
            nHI_frac(il) = 1.0_wp
          end if
        end select
      end if

      ! --- Line opacity: rhokap = n_scatterer * cross0 / Dfreq * distance2cm ---
      if (have_nion) then
        nH = nion_arr(il)
      else
        select case (trim(par%ion_model))
        case ('solar_cie')
          if (have_Z) then
            nH = solar_ion_density(nH_cgs(il), Z_arr(il), T_cgs(il), trim(line%ion_id))
          else if (par%metallicity_global >= 0.0_wp) then
            nH = solar_ion_density(nH_cgs(il), par%metallicity_global, T_cgs(il), trim(line%ion_id))
          else
            nH = nH_cgs(il) * nHI_frac(il)
          end if
        case ('from_file')
          write(6,'(a)') 'ERROR: ion_model=from_file but n_ion column not found.'
          stop 'ion_model=from_file requires n_ion column'
        case default  ! 'none'
          nH = nH_cgs(il) * nHI_frac(il)
        end select
      end if
      amr_grid%rhokap(il) = nH * line%cross0 / amr_grid%Dfreq(il) * par%distance2cm

      ! --- Dust opacity: file ndust > dust_model > global DGR ---
      if (have_ndust) then
        amr_grid%rhokapD(il) = ndust_arr(il) * par%cext_dust * par%distance2cm
      else
        select case (trim(par%dust_model))
        case ('laursen09')
          if (have_Z) then
            amr_grid%rhokapD(il) = laursen09_ndust(nH_cgs(il), nHI_frac(il), &
                Z_arr(il), par%Z_ref, par%f_ion_dust) * par%cext_dust * par%distance2cm
          else if (par%metallicity_global >= 0.0_wp) then
            amr_grid%rhokapD(il) = laursen09_ndust(nH_cgs(il), nHI_frac(il), &
                par%metallicity_global, par%Z_ref, par%f_ion_dust) * par%cext_dust * par%distance2cm
          end if
        case ('from_file')
          write(6,'(a)') 'ERROR: dust_model=from_file but ndust column not found in file.'
          stop 'dust_model=from_file requires ndust column'
        case default  ! 'global_dgr'
          if (par%DGR > 0.0_wp) then
            amr_grid%rhokapD(il) = nH_cgs(il) * par%cext_dust * par%DGR * par%distance2cm
          end if
        end select
      end if

      ! Velocity in units of v_thermal (same convention as Cartesian code)
      amr_grid%vfx(il) = vx_kms(il) / vtherm
      amr_grid%vfy(il) = vy_kms(il) / vtherm
      amr_grid%vfz(il) = vz_kms(il) / vtherm
    end do

    !--- Placeholder leaves (filled by amr_build_tree to cover gaps in the
    !    octree) get safe defaults: zero opacity (already zero-initialized
    !    by create_shared_mem) and a non-zero Dfreq/voigt_a to avoid NaNs
    !    in any divide-by-Dfreq path.
    if (amr_grid%nleaf > nleaf) then
      do il = nleaf + 1, amr_grid%nleaf
        amr_grid%Dfreq(il)   = amr_grid%Dfreq_ref
        amr_grid%voigt_a(il) = amr_grid%voigt_amean
      end do
    end if

    !--- Optional override: apply analytic velocity_type model on top of the
    !    generic AMR data velocities.  This lets the user reuse the same
    !    octree structure (density, temperature, refinement) while testing
    !    different velocity scalings without regenerating the AMR file.
    !    Overrides file velocities with analytic model (Hubble, SSH, etc.).
    if (len_trim(par%velocity_type) > 0) then
       call assign_amr_velocities_from_type(T_cgs)
    end if
    end if
    call MPI_BARRIER(mpar%hostcomm, ierr)

    ! Keep nHI_frac and emiss_arr alive if needed by amr_setup_emissivity.
    if (trim(par%source_geometry) /= 'diffuse_emissivity') then
      deallocate(nHI_frac)
      if (allocated(emiss_arr)) deallocate(emiss_arr)
    end if
    if (allocated(Z_arr))     deallocate(Z_arr)
    if (allocated(xHI_arr))   deallocate(xHI_arr)
    if (allocated(ne_arr))    deallocate(ne_arr)
    if (allocated(nion_arr))  deallocate(nion_arr)
    if (allocated(ndust_arr)) deallocate(ndust_arr)

    !--- Step 4b: Biconical mask — zero opacity outside cone (z-axis) --
    if (par%cone_opening > 0.0_wp .and. par%cone_opening < 90.0_wp) then
      cos_cone = cos(par%cone_opening * deg2rad)
      do il = 1, nleaf
        icell = amr_grid%icell_of_leaf(il)
        rr = sqrt(amr_grid%cx(icell)**2 + amr_grid%cy(icell)**2 + amr_grid%cz(icell)**2)
        if (rr > 0.0_wp) then
          if (abs(amr_grid%cz(icell)) / rr < cos_cone) then
            amr_grid%rhokap(il) = 0.0_wp
            if (par%DGR > 0.0_wp .and. associated(amr_grid%rhokapD)) &
                amr_grid%rhokapD(il) = 0.0_wp
          end if
        end if
      end do
    end if

    !--- Step 5: Normalize to input optical depth / column density -----
    ! tau_pole: traverse from box center in +z direction to compute actual
    ! optical depth from center to z+ boundary (matches Cartesian taumax convention).
    opacity_sum = 0.0_wp
    do il = 1, nleaf
      opacity_sum = opacity_sum + amr_grid%rhokap(il) * voigt(0.0_wp, amr_grid%voigt_a(il))
    end do
    nopac       = real(count(amr_grid%rhokap > 0.0_wp), wp)
    ! Half-box = center-to-boundary distance; matches Cartesian convention where
    ! opac_length = par%rmax (center to sphere surface) for a sphere model.
    opac_length = boxlen_code / 2.0_wp
    if (nopac > 0.0_wp) then
      tauhomo = (opacity_sum / nopac) * opac_length
    else
      tauhomo = 0.0_wp
    end if

    ! Pole traversal from box center (+z direction): compute tau and N(HI) simultaneously.
    ! This matches the Cartesian convention (tau from center to sphere surface).
    ! When the octree only partially covers the box (e.g. when the user passes a
    ! subset of a larger simulation), the traversal may step into a region where
    ! no leaf cell exists.  We skip across that gap one virtual sub-cell at a
    ! time (zero opacity contribution) so the traversal can resume on the far
    ! side and reach the box boundary, giving a correct N(HI)_pole.
    NHI_pole_raw = 0.0_wp
    block
      integer  :: il_cur, icell_cur, iface_cur, icell_gap
      real(wp) :: xc, yc, zc, t_exit_cur, tau_raw_half
      integer  :: niter
      integer, parameter :: NITER_MAX = 10000000
      xc = amr_grid%cx(1)
      yc = amr_grid%cy(1)
      zc = amr_grid%cz(1)
      tau_raw_half = 0.0_wp
      niter        = 0
      do
        niter = niter + 1
        if (niter > NITER_MAX) exit
        if (zc >= amr_grid%zmax) exit

        il_cur = amr_find_leaf(xc, yc, zc)
        if (il_cur > 0) then
          icell_cur = amr_grid%icell_of_leaf(il_cur)
          call amr_cell_exit(icell_cur, xc, yc, zc, 0.0_wp, 0.0_wp, 1.0_wp, t_exit_cur, iface_cur)
          tau_raw_half = tau_raw_half + amr_grid%rhokap(il_cur) &
                         * voigt(0.0_wp, amr_grid%voigt_a(il_cur)) * t_exit_cur
          NHI_pole_raw = NHI_pole_raw + amr_grid%rhokap(il_cur) &
                         * amr_grid%Dfreq(il_cur) / line%cross0 * t_exit_cur
          zc = zc + t_exit_cur
        else
          icell_gap = amr_find_enclosing_cell(xc, yc, zc)
          if (icell_gap <= 0) exit
          call amr_gap_exit(icell_gap, xc, yc, zc, 0.0_wp, 0.0_wp, 1.0_wp, t_exit_cur, iface_cur)
          zc = zc + t_exit_cur
        end if
      end do
      if (tau_raw_half > 0.0_wp) then
        taupole = tau_raw_half
      else
        taupole = tauhomo
      end if
    end block

    N_HIpole = NHI_pole_raw
    N_HIhomo = N_HIpole

    ! Shared arrays: only rank 0 on each node normalizes; others wait at barrier.
    if (mpar%h_rank == 0) then
    if (par%taumax > 0.0_wp) then
      if (taupole > 0.0_wp) then
        opac_norm           = par%taumax / taupole
        amr_grid%rhokap(:) = amr_grid%rhokap(:) * opac_norm
        if (par%DGR > 0.0_wp .and. associated(amr_grid%rhokapD)) &
            amr_grid%rhokapD(:) = amr_grid%rhokapD(:) * opac_norm
      end if
    else if (par%N_HImax > 0.0_wp) then
      if (N_HIpole > 0.0_wp) then
        opac_norm           = par%N_HImax / N_HIpole
        amr_grid%rhokap(:) = amr_grid%rhokap(:) * opac_norm
        if (par%DGR > 0.0_wp .and. associated(amr_grid%rhokapD)) &
            amr_grid%rhokapD(:) = amr_grid%rhokapD(:) * opac_norm
      end if
    else if (par%N_gasmax > 0.0_wp) then
      ! N_gasmax: ion/gas column density from center to +z boundary.
      ! N_HIpole was computed as sum(rhokap * Dfreq / cross0 * t_exit) along +z pole,
      ! which gives the column density in the same units as N_gasmax.
      if (N_HIpole > 0.0_wp) then
        opac_norm           = par%N_gasmax / N_HIpole
        amr_grid%rhokap(:) = amr_grid%rhokap(:) * opac_norm
        if (par%DGR > 0.0_wp .and. associated(amr_grid%rhokapD)) &
            amr_grid%rhokapD(:) = amr_grid%rhokapD(:) * opac_norm
      end if
    end if
    end if
    call MPI_BARRIER(mpar%hostcomm, ierr)

    ! Recompute tauhomo and taupole from normalized rhokap.
    opacity_sum = 0.0_wp
    do il = 1, nleaf
      opacity_sum = opacity_sum + amr_grid%rhokap(il) * voigt(0.0_wp, amr_grid%voigt_a(il))
    end do
    if (nopac > 0.0_wp) then
      tauhomo = (opacity_sum / nopac) * opac_length
    end if
    ! taupole = taumax after normalization (pole traversal gave exact normalization)
    if (par%taumax > 0.0_wp) then
      taupole = par%taumax
    else
      taupole = tauhomo
    end if

    N_HIpole = NHI_pole_raw * opac_norm
    N_HIhomo = N_HIpole

    ! Update par values (used by xcrit and frequency grid)
    if (par%taumax  <= 0.0_wp) par%taumax  = taupole
    if (par%tauhomo <= 0.0_wp) par%tauhomo = tauhomo
    if (par%N_HImax  <= 0.0_wp) par%N_HImax  = N_HIpole
    if (par%N_HIhomo <= 0.0_wp) par%N_HIhomo = N_HIhomo
    if (par%N_gasmax <= 0.0_wp) par%N_gasmax = N_HIpole
    par%tauhomo = tauhomo
    par%taumax  = taupole

    !--- Step 6: Core-skip parameters ---------------------------------
    atau0      = amr_grid%voigt_amean * par%tauhomo
    ! ch(1) = boxlen_code/2 in code units; number of root cells along one axis = 1.
    atau0_cell = atau0 / (boxlen_code / (2.0_wp * amr_grid%ch(1)))
    if (atau0_cell <= 1.0_wp) then
      amr_grid%xcrit  = 0.0_wp
      amr_grid%xcrit2 = 0.0_wp
    else
      if (atau0_cell <= 60.0_wp) then; xi = 0.6_wp; chi = 1.2_wp
      else;                            xi = 1.4_wp; chi = 0.6_wp
      end if
      amr_grid%xcrit  = 0.02_wp * exp(xi * (log(atau0_cell))**chi)
      amr_grid%xcrit2 = amr_grid%xcrit ** 2
    end if

    !--- Step 7: Diffuse emissivity setup ------------------------------
    grid%xmin   = amr_grid%xmin
    grid%xmax   = amr_grid%xmax
    grid%ymin   = amr_grid%ymin
    grid%ymax   = amr_grid%ymax
    grid%zmin   = amr_grid%zmin
    grid%zmax   = amr_grid%zmax
    grid%xrange = amr_grid%xrange
    grid%yrange = amr_grid%yrange
    grid%zrange = amr_grid%zrange
    ! par%rmax/xmax/ymax/zmax already set by geometry block in Step 2.

    call amr_setup_emissivity(grid, nH_cgs, T_cgs, nHI_frac, emiss_arr)
    if (allocated(emiss_arr))  deallocate(emiss_arr)
    if (allocated(nHI_frac))   deallocate(nHI_frac)

    !--- Step 8: Set up frequency grid --------------------------------
    call amr_setup_freq_grid(atau0, boxlen_code * par%distance2cm)

    !--- Step 9: Allocate output arrays -------------------------------
    call amr_alloc_output(par%save_Jin, par%save_Jabs, par%DGR > 0.0_wp, par%save_Jmu, par%nmu)

    !--- Print summary ------------------------------------------------
    if (mpar%p_rank == 0) then
      write(6,'(a,es12.4)') 'AMR voigt_a      : ', amr_grid%voigt_amean
      write(6,'(a,es12.4)') 'AMR temperature  : ', par%temperature
      write(6,'(a,es12.4)') 'AMR N(HI)_pole   : ', N_HIpole
      write(6,'(a,es12.4)') 'AMR tau_pole      : ', taupole
    end if

    deallocate(xleaf, yleaf, zleaf, leaf_level, nH_cgs, T_cgs, vx_kms, vy_kms, vz_kms)
  end subroutine grid_create_amr

  !=========================================================================
  ! Set up AMR frequency grid (mirrors grid_mod_car logic).
  !=========================================================================
  subroutine amr_setup_freq_grid(atau0, boxlen_cm)
    use define
    use mpi
    use memory_mod
    use line_mod
    implicit none
    real(wp), intent(in) :: atau0, boxlen_cm

    real(wp) :: vtherm, xscale, atau1, atau0_arg
    integer  :: i, ierr

    vtherm = line%vtherm1 * sqrt(par%temperature)

    ! Translate wavelength/velocity range inputs to x-frequency range (mirrors grid_mod_car)
    if (is_finite(par%wavelength_min) .and. is_finite(par%wavelength_max)) then
      if (par%nwavelength == 0 .and. par%nxfreq > 0) par%nwavelength = par%nxfreq
      if (par%nwavelength > 0) par%nxfreq = par%nwavelength
      par%xfreq_min = -(par%wavelength_max - line%wavelength0*1.0e4_wp) / &
                       (line%wavelength0*1.0e4_wp) * (speedc/vtherm)
      par%xfreq_max = -(par%wavelength_min - line%wavelength0*1.0e4_wp) / &
                       (line%wavelength0*1.0e4_wp) * (speedc/vtherm)
    else if (is_finite(par%velocity_min) .and. is_finite(par%velocity_max)) then
      if (par%nvelocity == 0 .and. par%nxfreq > 0) par%nvelocity = par%nxfreq
      if (par%nvelocity > 0) par%nxfreq = par%nvelocity
      par%xfreq_min = -par%velocity_max / vtherm
      par%xfreq_max = -par%velocity_min / vtherm
    end if

    atau1 = (amr_grid%voigt_amean * par%tauhomo)**(1.0_wp/3.0_wp)

    if (.not. (is_finite(par%xfreq_max) .and. is_finite(par%xfreq_min))) then
      if (par%taumax <= 5.0e1_wp) then;  xscale = 25.0_wp
      else if (par%taumax <= 5.0e2_wp) then; xscale = 14.0_wp
      else if (par%taumax <= 5.0e3_wp) then; xscale = 10.0_wp
      else;                               xscale =  5.0_wp
      end if
      if (par%Vexp == 0.0_wp) then
        par%xfreq_max = floor(xscale * atau1) + 1.0_wp
        par%xfreq_min = -(floor(xscale * atau1 + line%DnuHK_Hz/amr_grid%Dfreq_ref) + 1.0_wp)
      else if (par%Vexp > 0.0_wp) then
        par%xfreq_max = floor(xscale * atau1) + 1.0_wp
        par%xfreq_min = -(floor(xscale * atau1 + abs(par%Vexp)/vtherm &
                            + line%DnuHK_Hz/amr_grid%Dfreq_ref) + 1.0_wp)
      else
        par%xfreq_max = floor(xscale * atau1 + abs(par%Vexp)/vtherm) + 1.0_wp
        par%xfreq_min = -(floor(xscale * atau1 + line%DnuHK_Hz/amr_grid%Dfreq_ref) + 1.0_wp)
      end if
      if (trim(par%spectral_type) == 'continuum') then
        xscale        = 4.0_wp * xscale
        par%xfreq_max = floor(xscale * atau1 + abs(par%Vexp)/vtherm) + 1.0_wp
        par%xfreq_min = -(floor(xscale * atau1 + abs(par%Vexp)/vtherm &
                            + line%DnuHK_Hz/amr_grid%Dfreq_ref) + 1.0_wp)
      end if
    end if

    amr_grid%nxfreq   = par%nxfreq
    amr_grid%xfreq_min = par%xfreq_min
    amr_grid%xfreq_max = par%xfreq_max
    amr_grid%dxfreq    = (par%xfreq_max - par%xfreq_min) / par%nxfreq

    call create_shared_mem(amr_grid%xfreq,     [par%nxfreq])
    call create_shared_mem(amr_grid%velocity,  [par%nxfreq])
    call create_shared_mem(amr_grid%wavelength,[par%nxfreq])
    ! Only h_rank==0 writes to shared arrays; all ranks compute the scalar dwave.
    if (mpar%h_rank == 0) then
      do i = 1, par%nxfreq
        amr_grid%xfreq(i)      = (i - 0.5_wp) * amr_grid%dxfreq + amr_grid%xfreq_min
        amr_grid%velocity(i)   = -vtherm * amr_grid%xfreq(i)
        amr_grid%wavelength(i) = (amr_grid%velocity(i)/speedc + 1.0_wp) * (line%wavelength0 * 1.0e4_wp)
      end do
    end if
    amr_grid%dwave = vtherm / speedc * (line%wavelength0 * 1.0e4_wp) * amr_grid%dxfreq

    !--- Compute line photon fraction for 'continuum+gaussian' spectral type.
    !    (Same logic as grid_mod_car.f90; needed by output_normalize_amr.)
    if (trim(par%spectral_type) == 'continuum+gaussian' .and. par%EW_line > 0.0_wp) then
       block
          real(kind=wp) :: EW_vel, dv_range
          real(kind=wp), parameter :: speedc_kms = 2.99792458e5_wp
          EW_vel   = par%EW_line / (line%wavelength0 * 1.0e4_wp) * speedc_kms
          dv_range = (amr_grid%xfreq_max - amr_grid%xfreq_min) * line%vtherm1 * sqrt(par%temperature)
          par%f_line = EW_vel / (EW_vel + dv_range)
          if (mpar%p_rank == 0) then
             write(*,'(a,f8.4)') ' continuum+gaussian: f_line = ', par%f_line
             write(*,'(a,f8.2,a)') '   EW_line = ', par%EW_line, ' A'
             write(*,'(a,f8.2,a)') '   dv_range = ', dv_range, ' km/s'
          endif
       end block
    endif

    call MPI_BARRIER(mpar%hostcomm, ierr)
  end subroutine amr_setup_freq_grid



  !=========================================================================
  ! Prepare diffuse-emissivity sampling data for AMR runs.
  !
  ! Only runs if par%source_geometry == 'diffuse_emissivity'. Otherwise
  ! returns immediately so that par%emiss_file (which may double as a
  ! mode keyword like 'density1') has no side effects when not requested.
  !=========================================================================
  subroutine amr_setup_emissivity(grid, nH_cgs, T_cgs, nHI_frac, emiss_arr)
    use define
    use physics_amr_mod, only: electron_density_from_xHI, caseB_lya_emissivity
    use mpi
    implicit none
    type(grid_type), intent(in) :: grid
    real(wp),              intent(in) :: nH_cgs(:), T_cgs(:), nHI_frac(:)
    real(wp), allocatable, intent(in) :: emiss_arr(:)

    integer  :: il, ierr
    real(wp) :: cell_volume, total_positive_volume, norm
    real(wp) :: f_comp, f_comp1
    real(wp) :: ne_il
    logical  :: have_emiss
    logical  :: did_setup

    ! Nothing to do if the source is not diffuse_emissivity.
    if (trim(par%source_geometry) /= 'diffuse_emissivity') return

    have_emiss = allocated(emiss_arr)
    did_setup  = .false.

    ! --- Emissivity model dispatch ------------------------------------
    ! Priority order:
    !   1) explicit par%emissivity_model parameter
    !   2) par%emiss_file (filename or keyword)
    !   3) auto-detect emissivity column in the grid file
    select case (trim(par%emissivity_model))
    case ('caseB')
      ! Case B recombination + collisional excitation from nH, T, xHI.
      call create_shared_mem(amr_grid%Pem, [amr_grid%nleaf])
      if (mpar%h_rank == 0) then
        do il = 1, amr_grid%nleaf
          ne_il = electron_density_from_xHI(nH_cgs(il), nHI_frac(il))
          amr_grid%Pem(il) = caseB_lya_emissivity(nH_cgs(il), T_cgs(il), &
                                                   nHI_frac(il), ne_il)
        end do
      end if
      did_setup = .true.
    case ('from_file')
      if (.not. have_emiss) then
        write(6,'(a)') 'ERROR: emissivity_model=from_file but the generic AMR file'
        write(6,'(a)') '       does not contain an emissivity column.'
        stop 'amr_setup_emissivity: emissivity_model=from_file requires emissivity column'
      end if
      call create_shared_mem(amr_grid%Pem, [amr_grid%nleaf])
      if (mpar%h_rank == 0) amr_grid%Pem(:) = emiss_arr(:)
      did_setup = .true.
    case ('', 'none')
      ! Fall through to emiss_file-based logic, or auto-detect emissivity column.
      if (len_trim(par%emiss_file) <= 0) then
        if (have_emiss) then
          ! Auto-use emissivity column from the grid file.
          if (mpar%p_rank == 0) &
            write(6,'(a)') 'AMR emissivity: using emissivity column from grid file'
          call create_shared_mem(amr_grid%Pem, [amr_grid%nleaf])
          if (mpar%h_rank == 0) amr_grid%Pem(:) = emiss_arr(:)
          did_setup = .true.
        else
          write(6,'(a)') 'ERROR: source_geometry=diffuse_emissivity but no emissivity source.'
          write(6,'(a)') '       Set par%emissivity_model to caseB/from_file, or set'
          write(6,'(a)') '       par%emiss_file, or include an emissivity column in'
          write(6,'(a)') '       the AMR grid file.'
          stop 'amr_setup_emissivity: missing emissivity configuration'
        end if
      end if

      if (.not. did_setup) then
        select case(trim(get_extension(par%emiss_file)))
        case ('txt', 'dat')
          if (trim(par%geometry) == 'plane_atmosphere') then
            call setup_plane_emissivity(trim(par%emiss_file), emiss_prof, grid, &
                                        par%sampling_method, par%f_composite)
          else
            call setup_spherical_emissivity(trim(par%emiss_file), emiss_prof, grid, &
                                            par%sampling_method, par%f_composite)
          end if
          return   ! 1D profile path: emiss_prof is populated, skip cell alias setup.
        case ('density1')
          call create_shared_mem(amr_grid%Pem, [amr_grid%nleaf])
          if (mpar%h_rank == 0) amr_grid%Pem(:) = amr_grid%rhokap(:)
          did_setup = .true.
        case ('density2')
          call create_shared_mem(amr_grid%Pem, [amr_grid%nleaf])
          if (mpar%h_rank == 0) amr_grid%Pem(:) = amr_grid%rhokap(:)**2
          did_setup = .true.
        case default
          write(6,'(a)') 'AMR diffuse emissivity supports: caseB, from_file,'
          write(6,'(a)') 'txt/dat profile, density1, density2, or grid-file column.'
          stop 'amr_setup_emissivity: unsupported emissivity input for AMR mode'
        end select
      end if
    case default
      write(6,'(a,a)') 'ERROR: unknown par%emissivity_model = ', trim(par%emissivity_model)
      stop 'amr_setup_emissivity: unsupported emissivity_model'
    end select

    call MPI_BARRIER(mpar%hostcomm, ierr)
    if (.not. did_setup) return
    select case(par%sampling_method)
    case (0)
      call create_shared_mem(amr_grid%alias, [amr_grid%nleaf])
      if (mpar%h_rank == 0) then
        do il = 1, amr_grid%nleaf
          cell_volume = (2.0_wp * amr_grid%ch(amr_grid%icell_of_leaf(il)))**3
          amr_grid%Pem(il) = amr_grid%Pem(il) * cell_volume
        end do
        norm = sum(amr_grid%Pem)
        if (norm > 0.0_wp) then
          amr_grid%Pem(:) = amr_grid%Pem(:) / norm
          call random_alias_setup(amr_grid%Pem, amr_grid%alias)
        else
          write(6,'(a)') 'ERROR: total emissivity is zero across all AMR leaf cells.'
          stop 'amr_setup_emissivity: emissivity sum is zero (sampling_method=0)'
        end if
      end if
    case (1)
      call create_shared_mem(amr_grid%Pwgt,  [amr_grid%nleaf])
      call create_shared_mem(amr_grid%alias, [amr_grid%nleaf])
      if (mpar%h_rank == 0) then
        do il = 1, amr_grid%nleaf
          cell_volume = (2.0_wp * amr_grid%ch(amr_grid%icell_of_leaf(il)))**3
          amr_grid%Pem(il) = amr_grid%Pem(il) * cell_volume
        end do

        norm = sum(amr_grid%Pem)
        if (norm > 0.0_wp) then
          amr_grid%Pem(:) = amr_grid%Pem(:) / norm
          total_positive_volume = 0.0_wp
          do il = 1, amr_grid%nleaf
            if (amr_grid%Pem(il) > 0.0_wp) then
              total_positive_volume = total_positive_volume + &
                  (2.0_wp * amr_grid%ch(amr_grid%icell_of_leaf(il)))**3
            end if
          end do

          f_comp  = par%f_composite
          f_comp1 = 1.0_wp - f_comp
          do il = 1, amr_grid%nleaf
            if (amr_grid%Pem(il) > 0.0_wp) then
              cell_volume = (2.0_wp * amr_grid%ch(amr_grid%icell_of_leaf(il)))**3
              amr_grid%Pwgt(il) = amr_grid%Pem(il) / &
                  (f_comp1 * amr_grid%Pem(il) + f_comp * cell_volume / total_positive_volume)
              amr_grid%Pem(il) = f_comp1 * amr_grid%Pem(il) + &
                  f_comp * cell_volume / total_positive_volume
            else
              amr_grid%Pwgt(il) = 0.0_wp
            end if
          end do
          call random_alias_setup(amr_grid%Pem, amr_grid%alias)
        else
          write(6,'(a)') 'ERROR: total emissivity is zero across all AMR leaf cells.'
          stop 'amr_setup_emissivity: emissivity sum is zero (sampling_method=1)'
        end if
      end if
    case default
      if (mpar%h_rank == 0) then
        norm = maxval(amr_grid%Pem)
        if (norm > 0.0_wp) amr_grid%Pem(:) = amr_grid%Pem(:) / norm
      end if
    end select
    call MPI_BARRIER(mpar%hostcomm, ierr)
  end subroutine amr_setup_emissivity

  !=========================================================================
  ! Point grid spectral pointer arrays at amr_grid data and copy scalar
  ! metadata.  Must be called after grid_create_amr and before any output
  ! routine that takes a grid_type argument.
  !=========================================================================
  subroutine amr_sync_to_grid(grid)
    use define
    implicit none
    type(grid_type), intent(inout) :: grid
    integer :: n

    n = amr_grid%nxfreq
    grid%nxfreq    = n
    grid%xfreq_min = amr_grid%xfreq_min
    grid%xfreq_max = amr_grid%xfreq_max
    grid%dxfreq    = amr_grid%dxfreq
    grid%dwave     = amr_grid%dwave
    grid%Dfreq_ref = amr_grid%Dfreq_ref
    grid%voigt_amean = amr_grid%voigt_amean
    ! nx/ny/nz set to 1 so nx==ny==nz check passes but rmax=0 forces box area
    grid%nx = 1
    grid%ny = 1
    grid%nz = 1

    ! Allocate grid pointer arrays (filled from amr_grid after MPI reduction)
    allocate(grid%xfreq(n),     source=amr_grid%xfreq)
    allocate(grid%velocity(n),  source=amr_grid%velocity)
    allocate(grid%wavelength(n),source=amr_grid%wavelength)
    allocate(grid%Jout(n));     grid%Jout = 0.0_wp
    if (allocated(amr_grid%Jin)) then
      allocate(grid%Jin(n));    grid%Jin  = 0.0_wp
    end if
    if (allocated(amr_grid%Jabs)) then
      allocate(grid%Jabs(n));   grid%Jabs = 0.0_wp
    end if
    if (allocated(amr_grid%Jmu)) then
      allocate(grid%Jmu(n, par%nmu)); grid%Jmu = 0.0_wp
    end if
  end subroutine amr_sync_to_grid

  !=========================================================================
  ! Deallocate all AMR grid arrays.
  !=========================================================================
  subroutine grid_destroy_amr()
    use define
    use memory_mod
    implicit none
    ! Free every shared-memory window in one shot.  Mirrors the convention
    ! used by grid_destroy (Cartesian) and observer_destroy_outside: any of
    ! those callers may run before us and zero out num_windows, in which case
    ! the per-array destroy_mem() route would fall through to a Fortran
    ! deallocate() of MPI-window memory and trip ifort severe (173).  Calling
    ! destroy_shared_mem_all() is idempotent (no-op when num_windows == 0).
    call destroy_shared_mem_all()
    nullify(amr_grid%neighbor, amr_grid%parent, amr_grid%children, &
            amr_grid%level, amr_grid%cx, amr_grid%cy, amr_grid%cz, amr_grid%ch, &
            amr_grid%ileaf, amr_grid%icell_of_leaf, &
            amr_grid%rhokap, amr_grid%rhokapD, amr_grid%Dfreq, amr_grid%voigt_a, &
            amr_grid%vfx, amr_grid%vfy, amr_grid%vfz, &
            amr_grid%Pem, amr_grid%Pwgt, amr_grid%alias, &
            amr_grid%xfreq, amr_grid%velocity, amr_grid%wavelength)
    ! Per-rank output arrays (allocatable): freed normally
    if (allocated(amr_grid%Jout))  deallocate(amr_grid%Jout)
    if (allocated(amr_grid%Jin))   deallocate(amr_grid%Jin)
    if (allocated(amr_grid%Jabs))  deallocate(amr_grid%Jabs)
    if (allocated(amr_grid%Jmu))   deallocate(amr_grid%Jmu)
  end subroutine grid_destroy_amr

  !=========================================================================
  ! assign_amr_velocities_from_type: override per-leaf vfx/vfy/vfz using the same
  ! analytic velocity-field models supported by Cartesian mode.  Called from
  ! grid_create_amr only when amr_type='generic' and par%velocity_type is set.
  !
  ! Supported par%velocity_type values (identical to Cartesian semantics):
  !   'hubble'              : v = Vexp * r / rmax              (linear expansion)
  !   'constant_radial'     : v = Vexp * r_hat                  (uniform outflow)
  !   'parallel_velocity'   : v = (Vx, Vy, Vz)                  (uniform bulk)
  !   'ssh'                 : Song, Seon & Hwang (2020) galaxy model
  !   'rotating_solid_body' : v = Vrot * (-y, x, 0) / rmax
  !   'rotating_galaxy_halo': flat rotation curve (Vrot, rinner)
  !
  ! Only h_rank==0 calls this (vfx/vfy/vfz are MPI shared-memory windows).
  ! Velocities are stored as v / v_thermal (dimensionless), per AMR convention.
  !=========================================================================
  subroutine assign_amr_velocities_from_type(T_cgs)
    use define
    use line_mod
    implicit none
    real(wp), intent(in) :: T_cgs(:)

    integer  :: il, icell
    real(wp) :: xc, yc, zc, rr, rr_cyl, vtherm, Vscale, rmax_eff

    rmax_eff = par%rmax
    if (rmax_eff <= 0.0_wp) rmax_eff = amr_grid%L_box * 0.5_wp

    if (mpar%p_rank == 0) then
       write(6,'(2a)') 'AMR generic: overriding input velocities with velocity_type = ', &
                       trim(par%velocity_type)
    end if

    select case (trim(par%velocity_type))

    case ('hubble')
       do il = 1, amr_grid%nleaf
          icell  = amr_grid%icell_of_leaf(il)
          vtherm = line%vtherm1 * sqrt(T_cgs(il))
          amr_grid%vfx(il) = (par%Vexp / vtherm) * amr_grid%cx(icell) / rmax_eff
          amr_grid%vfy(il) = (par%Vexp / vtherm) * amr_grid%cy(icell) / rmax_eff
          amr_grid%vfz(il) = (par%Vexp / vtherm) * amr_grid%cz(icell) / rmax_eff
       end do

    case ('constant_radial')
       do il = 1, amr_grid%nleaf
          icell  = amr_grid%icell_of_leaf(il)
          xc = amr_grid%cx(icell);  yc = amr_grid%cy(icell);  zc = amr_grid%cz(icell)
          rr = sqrt(xc*xc + yc*yc + zc*zc)
          if (rr > amr_grid%ch(icell) * 0.1_wp) then
             vtherm = line%vtherm1 * sqrt(T_cgs(il))
             amr_grid%vfx(il) = (par%Vexp / vtherm) * xc / rr
             amr_grid%vfy(il) = (par%Vexp / vtherm) * yc / rr
             amr_grid%vfz(il) = (par%Vexp / vtherm) * zc / rr
          else
             amr_grid%vfx(il) = 0.0_wp
             amr_grid%vfy(il) = 0.0_wp
             amr_grid%vfz(il) = 0.0_wp
          end if
       end do

    case ('parallel_velocity')
       do il = 1, amr_grid%nleaf
          vtherm = line%vtherm1 * sqrt(T_cgs(il))
          amr_grid%vfx(il) = par%Vx / vtherm
          amr_grid%vfy(il) = par%Vy / vtherm
          amr_grid%vfz(il) = par%Vz / vtherm
       end do

    case ('ssh')
       do il = 1, amr_grid%nleaf
          icell  = amr_grid%icell_of_leaf(il)
          xc = amr_grid%cx(icell);  yc = amr_grid%cy(icell);  zc = amr_grid%cz(icell)
          rr = sqrt(xc*xc + yc*yc + zc*zc)
          vtherm = line%vtherm1 * sqrt(T_cgs(il))
          if (rr < par%rpeak) then
             Vscale = par%Vpeak / par%rpeak
             amr_grid%vfx(il) = (Vscale / vtherm) * xc
             amr_grid%vfy(il) = (Vscale / vtherm) * yc
             amr_grid%vfz(il) = (Vscale / vtherm) * zc
          else if (rr > 0.0_wp) then
             Vscale = par%Vpeak + par%DeltaV * (rr - par%rpeak) / (rmax_eff - par%rpeak)
             amr_grid%vfx(il) = (Vscale / vtherm) * xc / rr
             amr_grid%vfy(il) = (Vscale / vtherm) * yc / rr
             amr_grid%vfz(il) = (Vscale / vtherm) * zc / rr
          else
             amr_grid%vfx(il) = 0.0_wp
             amr_grid%vfy(il) = 0.0_wp
             amr_grid%vfz(il) = 0.0_wp
          end if
       end do

    case ('rotating_solid_body')
       do il = 1, amr_grid%nleaf
          icell  = amr_grid%icell_of_leaf(il)
          vtherm = line%vtherm1 * sqrt(T_cgs(il))
          amr_grid%vfx(il) = -(par%Vrot / vtherm) * amr_grid%cy(icell) / rmax_eff
          amr_grid%vfy(il) =  (par%Vrot / vtherm) * amr_grid%cx(icell) / rmax_eff
          amr_grid%vfz(il) = 0.0_wp
       end do

    case ('rotating_galaxy_halo')
       do il = 1, amr_grid%nleaf
          icell  = amr_grid%icell_of_leaf(il)
          xc = amr_grid%cx(icell);  yc = amr_grid%cy(icell)
          rr_cyl = sqrt(xc*xc + yc*yc)
          vtherm = line%vtherm1 * sqrt(T_cgs(il))
          if (rr_cyl < par%rinner) then
             amr_grid%vfx(il) = -(par%Vrot / vtherm) * yc / par%rinner
             amr_grid%vfy(il) =  (par%Vrot / vtherm) * xc / par%rinner
          else if (rr_cyl > 0.0_wp) then
             amr_grid%vfx(il) = -(par%Vrot / vtherm) * yc / rr_cyl
             amr_grid%vfy(il) =  (par%Vrot / vtherm) * xc / rr_cyl
          else
             amr_grid%vfx(il) = 0.0_wp
             amr_grid%vfy(il) = 0.0_wp
          end if
          amr_grid%vfz(il) = 0.0_wp
       end do

    case default
       if (mpar%p_rank == 0) then
          write(6,'(3a)') 'AMR generic: WARNING velocity_type "', &
                          trim(par%velocity_type), &
                          '" not recognized - keeping input data velocities'
       end if
    end select

  end subroutine assign_amr_velocities_from_type

end module grid_mod_amr
