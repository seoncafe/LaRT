module grid_mod_amr
  !-----------------------------------------------------------------------
  ! AMR grid creation: reads data (RAMSES or generic format), builds the
  ! octree, and normalises physical quantities to LaRT's internal units.
  !
  ! Called from setup_v2c.f90 when par%use_amr_grid = .true.
  !
  ! After grid_create_amr returns:
  !   amr_grid%rhokap(il)   = HI opacity per unit length at line centre
  !   amr_grid%rhokapD(il)  = dust opacity per unit length  (if DGR > 0)
  !   amr_grid%Dfreq(il)    = local Doppler frequency
  !   amr_grid%voigt_a(il)  = Voigt damping parameter
  !   amr_grid%vfx/vfy/vfz  = velocity / v_thermal (dimensionless)
  !   amr_grid%xfreq(*)     = frequency array
  !   amr_grid%Jout(*), Jin(*), Jabs(*) allocated
  !-----------------------------------------------------------------------
  use octree_mod
  use read_ramses_amr_mod
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

    ! Working arrays
    real(wp), allocatable :: nHI_frac(:)   ! nHI/nH neutral fraction
    integer  :: il
    real(wp) :: vtherm, T4, k_ionize, k_recomb
    real(wp) :: nH, opacity_sum, nopac, opac_norm, opac_length
    real(wp) :: xscale, atau0, atau0_cell, xi, chi
    real(wp) :: taupole, N_HIpole, tauhomo, N_HIhomo, NHI_pole_raw
    real(wp) :: atau1
    integer  :: ix, ierr

    opac_norm = 1.0_wp
    !--- Step 1: Read leaf-cell data ----------------------------------
    if (mpar%p_rank == 0) then
      if (trim(par%amr_type) == 'ramses') then
        ! RAMSES reader returns positions in cm (unit_l from info.txt).
        ! For RAMSES, par%distance2cm is not meaningful; set it = 1 if unset.
        call ramses_read_leaf_cells(trim(par%amr_file), par%amr_snapnum, &
            xleaf, yleaf, zleaf, leaf_level, &
            nH_cgs, T_cgs, vx_kms, vy_kms, vz_kms, &
            nleaf, boxlen_code)
        if (par%distance2cm <= 0.0_wp) par%distance2cm = 1.0_wp
      else
        ! Generic reader returns code units (e.g., kpc).
        ! par%distance2cm (set by distance_unit in setup_v2c) converts to cm.
        call generic_amr_read(trim(par%amr_file), &
            xleaf, yleaf, zleaf, leaf_level, &
            nH_cgs, T_cgs, vx_kms, vy_kms, vz_kms, &
            nleaf, boxlen_code)
      end if
      write(6,'(a,i0)')    'AMR nleaf    : ', nleaf
      write(6,'(a,es12.4,a)') 'AMR boxlen   : ', boxlen_code, ' (code units)'
      write(6,'(a,es12.4,a)') 'AMR boxlen   : ', boxlen_code*par%distance2cm, ' cm'
    end if

    ! Broadcast nleaf and boxlen_code to all ranks
    call MPI_BCAST(nleaf,       1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(boxlen_code, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

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

    !--- Step 2: Build octree in code units ---------------------------
    ! Tree positions stay in code units (kpc, pc, etc.), exactly like the
    ! Cartesian grid.  Raytrace path lengths are in code units; rhokap is
    ! scaled to per-code-unit via *= distance2cm (same as Cartesian).
    ! (RAMSES reader returns cm, so for RAMSES distance2cm = 1 → no change.)
    ! Coordinates are in [-boxlen/2, +boxlen/2] (box centred at origin).
    call amr_build_tree(xleaf, yleaf, zleaf, leaf_level, nleaf, &
        -0.5_wp*boxlen_code, 0.5_wp*boxlen_code, &
        -0.5_wp*boxlen_code, 0.5_wp*boxlen_code, &
        -0.5_wp*boxlen_code, 0.5_wp*boxlen_code)

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
    call amr_alloc_phys(par%DGR > 0.0_wp)

    ! Shared arrays: only rank 0 on each node fills; others wait at barrier.
    if (mpar%h_rank == 0) then
    do il = 1, nleaf
      T_cgs(il) = max(T_cgs(il), 10.0_wp)
      vtherm    = line%vtherm1 * sqrt(T_cgs(il))
      amr_grid%Dfreq(il)   = vtherm / (line%wavelength0 * um2km)
      amr_grid%voigt_a(il) = (line%damping / fourpi) / amr_grid%Dfreq(il)

      ! Neutral fraction
      if (par%use_cie_condition) then
        T4       = T_cgs(il) / 1.0e4_wp
        k_ionize = 5.84862e-9_wp * sqrt(T4) * exp(-15.78215_wp / T4)
        k_recomb = 4.13e-13_wp * T4**(-0.7131_wp - 0.0115_wp*log(T4))
        nHI_frac(il) = k_recomb / (k_ionize + k_recomb)
      else
        nHI_frac(il) = 1.0_wp  ! assume fully neutral; user can modify
      end if

      ! rhokap = nHI * cross0 / Dfreq  [cm^-1]
      ! *= distance2cm converts to per-code-unit (same as Cartesian mode).
      nH = nH_cgs(il) * nHI_frac(il)
      amr_grid%rhokap(il) = nH * line%cross0 / amr_grid%Dfreq(il) * par%distance2cm

      ! Dust opacity (also per code unit)
      if (par%DGR > 0.0_wp) then
        amr_grid%rhokapD(il) = amr_grid%rhokap(il) * amr_grid%Dfreq(il) / line%cross0 &
                                * par%cext_dust * par%DGR
      end if

      ! Velocity in units of v_thermal (same convention as Cartesian code)
      amr_grid%vfx(il) = vx_kms(il) / vtherm
      amr_grid%vfy(il) = vy_kms(il) / vtherm
      amr_grid%vfz(il) = vz_kms(il) / vtherm
    end do

    !--- Optional override: apply analytic velocity_type model on top of the
    !    generic AMR data velocities.  This lets the user reuse the same
    !    octree structure (density, temperature, refinement) while testing
    !    different velocity scalings without regenerating the AMR file.
    !    Only enabled for amr_type='generic' (RAMSES velocities are physical).
    if (trim(par%amr_type) == 'generic' .and. len_trim(par%velocity_type) > 0) then
       call assign_amr_velocities_from_type(T_cgs)
    end if
    end if
    call MPI_BARRIER(mpar%hostcomm, ierr)

    deallocate(nHI_frac)

    !--- Step 5: Normalise to input optical depth / column density -----
    ! tau_pole: traverse from box center in +z direction to compute actual
    ! optical depth from center to z+ boundary (matches Cartesian taumax convention).
    opacity_sum = 0.0_wp
    do il = 1, nleaf
      opacity_sum = opacity_sum + amr_grid%rhokap(il) * voigt(0.0_wp, amr_grid%voigt_a(il))
    end do
    nopac       = real(count(amr_grid%rhokap > 0.0_wp), wp)
    ! Half-box = centre-to-boundary distance; matches Cartesian convention where
    ! opac_length = par%rmax (centre to sphere surface) for a sphere model.
    opac_length = boxlen_code / 2.0_wp
    if (nopac > 0.0_wp) then
      tauhomo = (opacity_sum / nopac) * opac_length
    else
      tauhomo = 0.0_wp
    end if

    ! Pole traversal from box center (+z direction): compute tau and N(HI) simultaneously.
    ! This matches the Cartesian convention (tau from center to sphere surface).
    NHI_pole_raw = 0.0_wp
    block
      integer  :: il_cur, icell_cur, il_next, iface_cur
      real(wp) :: xc, yc, zc, t_exit_cur, tau_raw_half
      xc = amr_grid%cx(1)
      yc = amr_grid%cy(1)
      zc = amr_grid%cz(1)
      il_cur = amr_find_leaf(xc, yc, zc)
      tau_raw_half = 0.0_wp
      if (il_cur > 0) then
        do
          icell_cur = amr_grid%icell_of_leaf(il_cur)
          call amr_cell_exit(icell_cur, xc, yc, zc, 0.0_wp, 0.0_wp, 1.0_wp, t_exit_cur, iface_cur)
          tau_raw_half  = tau_raw_half + amr_grid%rhokap(il_cur) &
                          * voigt(0.0_wp, amr_grid%voigt_a(il_cur)) * t_exit_cur
          NHI_pole_raw  = NHI_pole_raw + amr_grid%rhokap(il_cur) &
                          * amr_grid%Dfreq(il_cur) / line%cross0 * t_exit_cur
          zc = zc + t_exit_cur
          il_next = amr_next_leaf(icell_cur, iface_cur, xc, yc, zc)
          if (il_next <= 0) exit
          il_cur = il_next
        end do
      end if
      if (tau_raw_half > 0.0_wp) then
        taupole = tau_raw_half
      else
        taupole = tauhomo
      end if
    end block

    N_HIpole = NHI_pole_raw
    N_HIhomo = N_HIpole

    ! Shared arrays: only rank 0 on each node normalises; others wait at barrier.
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

    call amr_setup_emissivity(grid)

    !--- Step 8: Set up frequency grid --------------------------------
    call amr_setup_freq_grid(atau0, boxlen_code * par%distance2cm)

    !--- Step 9: Allocate output arrays -------------------------------
    call amr_alloc_output(par%save_Jin, par%save_Jabs, par%DGR > 0.0_wp)

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
    call MPI_BARRIER(mpar%hostcomm, ierr)
  end subroutine amr_setup_freq_grid



  !=========================================================================
  ! Prepare diffuse-emissivity sampling data for AMR runs.
  !=========================================================================
  subroutine amr_setup_emissivity(grid)
    use define
    use mpi
    implicit none
    type(grid_type), intent(in) :: grid

    integer  :: il, ierr
    real(wp) :: cell_volume, total_positive_volume, norm
    real(wp) :: f_comp, f_comp1

    if (len_trim(par%emiss_file) <= 0) return

    select case(trim(get_extension(par%emiss_file)))
    case ('txt', 'dat')
      if (trim(par%geometry) == 'plane_atmosphere') then
        call setup_plane_emissivity(trim(par%emiss_file), emiss_prof, grid, &
                                    par%sampling_method, par%f_composite)
      else
        call setup_spherical_emissivity(trim(par%emiss_file), emiss_prof, grid, &
                                        par%sampling_method, par%f_composite)
      end if
      return
    case ('density1')
      call create_shared_mem(amr_grid%Pem, [amr_grid%nleaf])
      if (mpar%h_rank == 0) amr_grid%Pem(:) = amr_grid%rhokap(:)
    case ('density2')
      call create_shared_mem(amr_grid%Pem, [amr_grid%nleaf])
      if (mpar%h_rank == 0) amr_grid%Pem(:) = amr_grid%rhokap(:)**2
    case default
      write(6,'(a)') 'AMR diffuse emissivity currently supports txt/dat, density1, and density2.'
      stop 'amr_setup_emissivity: unsupported emissivity input for AMR mode'
    end select
    call MPI_BARRIER(mpar%hostcomm, ierr)

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
  end subroutine amr_sync_to_grid

  !=========================================================================
  ! Deallocate all AMR grid arrays.
  !=========================================================================
  subroutine grid_destroy_amr()
    use define
    use memory_mod
    implicit none
    ! Shared-memory (pointer) arrays: freed via destroy_mem (MPI window release)
    if (associated(amr_grid%neighbor))      call destroy_mem(amr_grid%neighbor)
    if (associated(amr_grid%parent))        call destroy_mem(amr_grid%parent)
    if (associated(amr_grid%children))      call destroy_mem(amr_grid%children)
    if (associated(amr_grid%level))         call destroy_mem(amr_grid%level)
    if (associated(amr_grid%cx))            call destroy_mem(amr_grid%cx)
    if (associated(amr_grid%cy))            call destroy_mem(amr_grid%cy)
    if (associated(amr_grid%cz))            call destroy_mem(amr_grid%cz)
    if (associated(amr_grid%ch))            call destroy_mem(amr_grid%ch)
    if (associated(amr_grid%ileaf))         call destroy_mem(amr_grid%ileaf)
    if (associated(amr_grid%icell_of_leaf)) call destroy_mem(amr_grid%icell_of_leaf)
    if (associated(amr_grid%rhokap))        call destroy_mem(amr_grid%rhokap)
    if (associated(amr_grid%rhokapD))       call destroy_mem(amr_grid%rhokapD)
    if (associated(amr_grid%Dfreq))         call destroy_mem(amr_grid%Dfreq)
    if (associated(amr_grid%voigt_a))       call destroy_mem(amr_grid%voigt_a)
    if (associated(amr_grid%vfx))           call destroy_mem(amr_grid%vfx)
    if (associated(amr_grid%vfy))           call destroy_mem(amr_grid%vfy)
    if (associated(amr_grid%vfz))           call destroy_mem(amr_grid%vfz)
    if (associated(amr_grid%Pem))           call destroy_mem(amr_grid%Pem)
    if (associated(amr_grid%Pwgt))          call destroy_mem(amr_grid%Pwgt)
    if (associated(amr_grid%alias))         call destroy_mem(amr_grid%alias)
    if (associated(amr_grid%xfreq))         call destroy_mem(amr_grid%xfreq)
    if (associated(amr_grid%velocity))      call destroy_mem(amr_grid%velocity)
    if (associated(amr_grid%wavelength))    call destroy_mem(amr_grid%wavelength)
    ! Per-rank output arrays (allocatable): freed normally
    if (allocated(amr_grid%Jout))  deallocate(amr_grid%Jout)
    if (allocated(amr_grid%Jin))   deallocate(amr_grid%Jin)
    if (allocated(amr_grid%Jabs))  deallocate(amr_grid%Jabs)
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
                          '" not recognised - keeping input data velocities'
       end if
    end select

  end subroutine assign_amr_velocities_from_type

end module grid_mod_amr
