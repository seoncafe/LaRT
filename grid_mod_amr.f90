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
  use voigt_mod
  use utility
  implicit none

  public :: grid_create_amr, grid_destroy_amr, amr_sync_to_grid

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
    real(wp) :: taupole, N_HIpole, tauhomo, N_HIhomo
    real(wp) :: atau1
    integer  :: ix, ierr

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

    !--- Step 2: Build octree in cm -----------------------------------
    ! Convert code-unit positions to cm, then build the tree.
    ! Raytrace uses cm path lengths; rhokap will be scaled to cm^-1.
    ! (RAMSES: positions already in cm, distance2cm = 1; generic: positions in
    !  code units, distance2cm = kpc2cm etc.)
    xleaf(:) = xleaf(:) * par%distance2cm
    yleaf(:) = yleaf(:) * par%distance2cm
    zleaf(:) = zleaf(:) * par%distance2cm
    call amr_build_tree(xleaf, yleaf, zleaf, leaf_level, nleaf, &
        0.0_wp, boxlen_code*par%distance2cm, &
        0.0_wp, boxlen_code*par%distance2cm, &
        0.0_wp, boxlen_code*par%distance2cm)

    !--- Step 3: Set reference Doppler frequency ----------------------
    amr_grid%Dfreq_ref  = line%vtherm1 * sqrt(par%temperature) / (line%wavelength0 * um2km)
    amr_grid%voigt_amean = (line%damping / fourpi) / amr_grid%Dfreq_ref
    amr_grid%Dfreq_mean  = amr_grid%Dfreq_ref

    !--- Step 4: Allocate physical arrays and fill --------------------
    allocate(nHI_frac(nleaf))
    call amr_alloc_phys(par%DGR > 0.0_wp)

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

    deallocate(nHI_frac)

    !--- Step 5: Normalise to input optical depth / column density -----
    ! Approximate tau_pole: mean(rhokap_percode * voigt) * boxlen_code.
    ! rhokap is per-code-unit (= rhokap_cm * distance2cm), opac_length is
    ! boxlen in code units — same convention as Cartesian mode.
    opacity_sum = sum(amr_grid%rhokap * voigt_array(amr_grid%voigt_a, 0.0_wp))
    nopac       = real(count(amr_grid%rhokap > 0.0_wp), wp)
    if (nopac > 0.0_wp) then
      opac_length = boxlen_code    ! box length in code units
      tauhomo     = (opacity_sum / nopac) * opac_length
    else
      tauhomo = 0.0_wp
    end if
    taupole  = tauhomo  ! approximation for AMR

    ! N_HIpole [cm^-2]: rhokap_percode / distance2cm × Dfreq / cross0 × boxlen_code × distance2cm
    !                 = rhokap_percode × Dfreq / cross0 × boxlen_code
    N_HIpole = sum(amr_grid%rhokap * amr_grid%Dfreq / line%cross0) &
               / amr_grid%nleaf * boxlen_code
    N_HIhomo = N_HIpole

    if (par%taumax > 0.0_wp) then
      if (taupole > 0.0_wp) then
        opac_norm           = par%taumax / taupole
        amr_grid%rhokap(:) = amr_grid%rhokap(:) * opac_norm
        if (par%DGR > 0.0_wp .and. allocated(amr_grid%rhokapD)) &
            amr_grid%rhokapD(:) = amr_grid%rhokapD(:) * opac_norm
      end if
    else if (par%N_HImax > 0.0_wp) then
      if (N_HIpole > 0.0_wp) then
        opac_norm           = par%N_HImax / N_HIpole
        amr_grid%rhokap(:) = amr_grid%rhokap(:) * opac_norm
        if (par%DGR > 0.0_wp .and. allocated(amr_grid%rhokapD)) &
            amr_grid%rhokapD(:) = amr_grid%rhokapD(:) * opac_norm
      end if
    end if

    ! Recompute taupole and tauhomo from normalized rhokap (same as Cartesian).
    opacity_sum = sum(amr_grid%rhokap * voigt_array(amr_grid%voigt_a, 0.0_wp))
    if (nopac > 0.0_wp) then
      tauhomo = (opacity_sum / nopac) * opac_length
    end if
    taupole  = tauhomo

    N_HIpole = sum(amr_grid%rhokap * amr_grid%Dfreq / line%cross0) &
               / amr_grid%nleaf * boxlen_code
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
    ! amr_grid%ch(1) is the root cell half-width (= boxlen_cm/2 = boxlen_code*distance2cm/2).
    ! boxlen_code / (2*ch(1)/distance2cm) = number of root cells along one axis = 1 for cubic.
    atau0_cell = atau0 / (boxlen_code / (2.0_wp * amr_grid%ch(1) / par%distance2cm))
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

    !--- Step 7: Set up frequency grid --------------------------------
    call amr_setup_freq_grid(atau0, boxlen_code * par%distance2cm)

    !--- Step 8: Allocate output arrays -------------------------------
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
    use line_mod
    implicit none
    real(wp), intent(in) :: atau0, boxlen_cm

    real(wp) :: vtherm, xscale, atau1, atau0_arg
    integer  :: i

    atau1 = (amr_grid%voigt_amean * par%tauhomo)**(1.0_wp/3.0_wp)

    if (.not. (is_finite(par%xfreq_max) .and. is_finite(par%xfreq_min))) then
      if (par%taumax <= 5.0e1_wp) then;  xscale = 25.0_wp
      else if (par%taumax <= 5.0e2_wp) then; xscale = 14.0_wp
      else if (par%taumax <= 5.0e3_wp) then; xscale = 10.0_wp
      else;                               xscale =  5.0_wp
      end if
      vtherm = line%vtherm1 * sqrt(par%temperature)
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
    end if

    amr_grid%nxfreq   = par%nxfreq
    amr_grid%xfreq_min = par%xfreq_min
    amr_grid%xfreq_max = par%xfreq_max
    amr_grid%dxfreq    = (par%xfreq_max - par%xfreq_min) / par%nxfreq

    allocate(amr_grid%xfreq(par%nxfreq))
    allocate(amr_grid%velocity(par%nxfreq))
    allocate(amr_grid%wavelength(par%nxfreq))

    vtherm = line%vtherm1 * sqrt(par%temperature)
    do i = 1, par%nxfreq
      amr_grid%xfreq(i)     = (i - 0.5_wp) * amr_grid%dxfreq + amr_grid%xfreq_min
      amr_grid%velocity(i)  = -vtherm * amr_grid%xfreq(i)
      amr_grid%wavelength(i) = (amr_grid%velocity(i)/speedc + 1.0_wp) * (line%wavelength0 * 1.0e4_wp)
    end do
    amr_grid%dwave = vtherm / speedc * (line%wavelength0 * 1.0e4_wp) * amr_grid%dxfreq
  end subroutine amr_setup_freq_grid

  !=========================================================================
  ! Compute Voigt function for an array of voigt_a values at x = 0.
  !=========================================================================
  function voigt_array(va, x) result(v)
    use define
    implicit none
    real(wp), intent(in) :: va(:), x
    real(wp) :: v(size(va))
    integer  :: i
    do i = 1, size(va)
      v(i) = voigt(x, va(i))
    end do
  end function voigt_array

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
    implicit none
    if (allocated(amr_grid%parent))       deallocate(amr_grid%parent)
    if (allocated(amr_grid%children))     deallocate(amr_grid%children)
    if (allocated(amr_grid%level))        deallocate(amr_grid%level)
    if (allocated(amr_grid%cx))           deallocate(amr_grid%cx)
    if (allocated(amr_grid%cy))           deallocate(amr_grid%cy)
    if (allocated(amr_grid%cz))           deallocate(amr_grid%cz)
    if (allocated(amr_grid%ch))           deallocate(amr_grid%ch)
    if (allocated(amr_grid%ileaf))        deallocate(amr_grid%ileaf)
    if (allocated(amr_grid%icell_of_leaf)) deallocate(amr_grid%icell_of_leaf)
    if (allocated(amr_grid%rhokap))       deallocate(amr_grid%rhokap)
    if (allocated(amr_grid%rhokapD))      deallocate(amr_grid%rhokapD)
    if (allocated(amr_grid%Dfreq))        deallocate(amr_grid%Dfreq)
    if (allocated(amr_grid%voigt_a))      deallocate(amr_grid%voigt_a)
    if (allocated(amr_grid%vfx))          deallocate(amr_grid%vfx)
    if (allocated(amr_grid%vfy))          deallocate(amr_grid%vfy)
    if (allocated(amr_grid%vfz))          deallocate(amr_grid%vfz)
    if (allocated(amr_grid%xfreq))        deallocate(amr_grid%xfreq)
    if (allocated(amr_grid%velocity))     deallocate(amr_grid%velocity)
    if (allocated(amr_grid%wavelength))   deallocate(amr_grid%wavelength)
    if (allocated(amr_grid%Jout))         deallocate(amr_grid%Jout)
    if (allocated(amr_grid%Jin))          deallocate(amr_grid%Jin)
    if (allocated(amr_grid%Jabs))         deallocate(amr_grid%Jabs)
  end subroutine grid_destroy_amr

end module grid_mod_amr
