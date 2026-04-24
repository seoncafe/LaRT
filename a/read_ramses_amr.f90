module read_ramses_amr_mod
  !-----------------------------------------------------------------------
  ! Reads RAMSES AMR binary output and returns flat arrays of leaf cells.
  !
  ! Public interface:
  !   ramses_read_leaf_cells  -- reads AMR + hydro files, returns leaf data
  !   generic_amr_read        -- placeholder for other AMR formats
  !
  ! RAMSES convention assumed here:
  !   - Non-cosmological run (cosmo = .false.) OR cosmological (boxlen in Mpc/h)
  !   - Unit scales read from output_{snapnum}/info_{snapnum}.txt
  !   - If present, hydro_file_descriptor.txt is used to detect non-standard
  !     layouts such as RAMSES RMHD variants with velocity_x/y/z and
  !     thermal_pressure fields.
  !   - nHI computed from nH + T using either CIE or supplied neutral fraction
  !-----------------------------------------------------------------------
  use define
  use fitsio_mod
  use, intrinsic :: iso_fortran_env, only: int32
  implicit none
  private

  public :: ramses_read_leaf_cells
  public :: generic_amr_read

  ! Precision flag for RAMSES hydro output (4 = single, 8 = double)
  integer :: ramses_hydro_prec = 8
  integer :: ramses_density_var = 1
  integer :: ramses_velocity_var(3) = [2, 3, 4]
  integer :: ramses_thermo_var = 5
  character(len=16) :: ramses_velocity_layout = 'momentum'
  character(len=16) :: ramses_thermo_mode    = 'energy'
  logical :: ramses_descriptor_found = .false.

  real(kind=8), parameter :: massH_cgs = 1.6726d-24
  real(kind=8), parameter :: boltzmann_cgs = 1.381d-16
  real(kind=8), parameter :: mu_neutral = 1.22d0

contains

  !=========================================================================
  ! Main reader: scans all CPU output files and collects leaf cell data.
  !
  ! Outputs (allocatable arrays, caller must deallocate):
  !   xleaf, yleaf, zleaf  -- leaf centre coordinates [physical units set by unit_l]
  !   leaf_level           -- AMR level (0 = root)
  !   nH_cgs               -- total hydrogen number density [cm^-3]
  !   T_cgs                -- gas temperature [K]
  !   vx_cgs, vy_cgs, vz_cgs -- gas velocity [km/s]
  !   nleaf                -- number of leaf cells found
  !   boxlen_cm            -- box physical size [cm]
  !   nHI_frac             -- neutral fraction nHI/nH (allocated, caller fills or use CIE)
  !=========================================================================
  subroutine ramses_read_leaf_cells(repository, snapnum, &
      xleaf, yleaf, zleaf, leaf_level, &
      nH_cgs, T_cgs, vx_cgs, vy_cgs, vz_cgs, &
      nleaf, boxlen_cm)

    character(len=*), intent(in)               :: repository
    integer,          intent(in)               :: snapnum
    real(wp), allocatable, intent(out)         :: xleaf(:), yleaf(:), zleaf(:)
    integer,  allocatable, intent(out)         :: leaf_level(:)
    real(wp), allocatable, intent(out)         :: nH_cgs(:), T_cgs(:)
    real(wp), allocatable, intent(out)         :: vx_cgs(:), vy_cgs(:), vz_cgs(:)
    integer,               intent(out)         :: nleaf
    real(wp),              intent(out)         :: boxlen_cm

    ! RAMSES unit scales
    real(kind=8) :: unit_l, unit_d, unit_t, unit_v, boxlen_code
    real(kind=8) :: gamma_eos

    integer :: ncpu, nleaf_total
    integer :: nleaf_est
    integer :: il

    ! Temp storage (allocated to nleaf_total)
    real(wp), allocatable :: xl(:), yl(:), zl(:)
    real(wp), allocatable :: nH(:), Tgas(:)
    real(wp), allocatable :: vel_x(:), vel_y(:), vel_z(:)
    integer,  allocatable :: lvl(:)

    ! Read unit scales and number of CPUs from info file
    call ramses_read_info(repository, snapnum, ncpu, &
        unit_l, unit_d, unit_t, boxlen_code, gamma_eos)
    call ramses_detect_layout(repository, snapnum)

    unit_v    = unit_l / unit_t  ! cm/s
    boxlen_cm = boxlen_code * unit_l

    ! First pass: count leaf cells
    nleaf_total = ramses_count_leaves(repository, snapnum, ncpu)

    ! Allocate output arrays
    allocate(xl(nleaf_total), yl(nleaf_total), zl(nleaf_total))
    allocate(lvl(nleaf_total))
    allocate(nH(nleaf_total), Tgas(nleaf_total))
    allocate(vel_x(nleaf_total), vel_y(nleaf_total), vel_z(nleaf_total))

    ! Second pass: read actual data
    call ramses_read_all_cpus(repository, snapnum, ncpu, &
        unit_l, unit_d, unit_t, unit_v, gamma_eos, boxlen_code, &
        xl, yl, zl, lvl, nH, Tgas, vel_x, vel_y, vel_z, &
        nleaf_total)

    ! Convert positions from code units (fraction of boxlen) to physical [cm]
    xl = xl * boxlen_cm
    yl = yl * boxlen_cm
    zl = zl * boxlen_cm

    nleaf       = nleaf_total
    xleaf       = xl
    yleaf       = yl
    zleaf       = zl
    leaf_level  = lvl
    nH_cgs      = nH
    T_cgs       = Tgas
    vx_cgs      = vel_x / 1.0e5_wp  ! cm/s → km/s
    vy_cgs      = vel_y / 1.0e5_wp
    vz_cgs      = vel_z / 1.0e5_wp

    deallocate(xl, yl, zl, lvl, nH, Tgas, vel_x, vel_y, vel_z)
  end subroutine ramses_read_leaf_cells

  !=========================================================================
  ! Read the RAMSES info file to get unit scales and ncpu.
  !=========================================================================
  subroutine ramses_read_info(repository, snapnum, ncpu, &
      unit_l, unit_d, unit_t, boxlen, gamma_eos)
    character(len=*), intent(in)  :: repository
    integer,          intent(in)  :: snapnum
    integer,          intent(out) :: ncpu
    real(kind=8),     intent(out) :: unit_l, unit_d, unit_t, boxlen, gamma_eos

    character(len=512) :: filename, line
    integer :: ios, unit
    character(len=20) :: key
    real(kind=8) :: val

    write(filename, '(a,"/output_",i5.5,"/info_",i5.5,".txt")') &
        trim(repository), snapnum, snapnum

    ncpu      = 1
    unit_l    = 1.0d0
    unit_d    = 1.0d0
    unit_t    = 1.0d0
    boxlen    = 1.0d0
    gamma_eos = 5.0d0/3.0d0

    unit = 50
    open(unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(6,'(a)') 'read_ramses_amr: WARNING: cannot open info file: '//trim(filename)
      return
    end if

    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      ! Parse lines of the form "key = value"
      if (index(line, '=') > 0) then
        read(line, *, iostat=ios) key, val
        if (ios /= 0) then
          ! Try string format with '=' separator
          key = adjustl(line(1:index(line,'=')-1))
        end if
        key = adjustl(key)
        if     (trim(key) == 'ncpu')   then; read(line(index(line,'=')+1:), *, iostat=ios) ncpu
        else if (trim(key) == 'boxlen') then; read(line(index(line,'=')+1:), *, iostat=ios) boxlen
        else if (trim(key) == 'unit_l') then; read(line(index(line,'=')+1:), *, iostat=ios) unit_l
        else if (trim(key) == 'unit_d') then; read(line(index(line,'=')+1:), *, iostat=ios) unit_d
        else if (trim(key) == 'unit_t') then; read(line(index(line,'=')+1:), *, iostat=ios) unit_t
        else if (trim(key) == 'gamma')  then; read(line(index(line,'=')+1:), *, iostat=ios) gamma_eos
        end if
      end if
    end do
    close(unit)
  end subroutine ramses_read_info

  !=========================================================================
  ! Count total leaf cells by scanning AMR files.
  !=========================================================================
  integer function ramses_count_leaves(repository, snapnum, ncpu) result(nleaf_total)
    character(len=*), intent(in) :: repository
    integer,          intent(in) :: snapnum, ncpu

    character(len=512) :: filename
    integer :: icpu, iu
    integer :: ncpu_f, ndim, nx, ny, nz, nlevelmax, ngridmax, nboundary, ngrid_current
    integer :: ilevel, j, ngrida
    integer, allocatable :: ngridlevel(:,:), ngridfile(:,:), ngridbound(:,:)

    nleaf_total = 0
    iu = 99

    do icpu = 1, ncpu
      write(filename,'(a,"/output_",i5.5,"/amr_",i5.5,".out",i5.5)') &
          trim(repository), snapnum, snapnum, icpu
      open(iu, file=trim(filename), form='unformatted', status='old', action='read')
      read(iu) ncpu_f
      read(iu) ndim
      read(iu) nx, ny, nz
      read(iu) nlevelmax
      read(iu) ngridmax
      read(iu) nboundary
      read(iu) ngrid_current
      read(iu)  ! boxlen
      call skip_records(iu, 13)

      allocate(ngridlevel(ncpu_f, nlevelmax))
      allocate(ngridfile(ncpu_f + nboundary, nlevelmax))
      if (nboundary > 0) allocate(ngridbound(nboundary, nlevelmax))
      read(iu) ngridlevel
      ngridfile(1:ncpu_f, :) = ngridlevel
      read(iu)
      if (nboundary > 0) then
        read(iu); read(iu)
        read(iu) ngridbound
        ngridfile(ncpu_f+1:ncpu_f+nboundary, :) = ngridbound
        deallocate(ngridbound)
      end if
      call skip_records(iu, 6)

      ! Count leaves for this cpu
      do ilevel = 1, nlevelmax
        ngrida = ngridfile(icpu, ilevel)
        do j = 1, nboundary + ncpu_f
          if (ngridfile(j, ilevel) > 0) then
            call skip_records(iu, 3)          ! grid/next/prev indices
            call skip_records(iu, ndim)       ! grid positions
            call skip_records(iu, 1)          ! father
            call skip_records(iu, 2*ndim)     ! nbor
            if (j == icpu) then
              ! Read son array to count leaves
              call count_leaves_in_oct(iu, 2**ndim, ngrida, nleaf_total)
            else
              call skip_records(iu, 2**ndim)  ! son
            end if
            call skip_records(iu, 2 * 2**ndim)  ! cpu_map + ref_map
          end if
        end do
      end do

      deallocate(ngridlevel, ngridfile)
      close(iu)
    end do
  end function ramses_count_leaves

  !=========================================================================
  ! Read son array for one grid level/CPU and count leaf cells.
  !=========================================================================
  subroutine count_leaves_in_oct(iu, twotondim, ngrida, nleaf_total)
    integer, intent(in)    :: iu, twotondim, ngrida
    integer, intent(inout) :: nleaf_total
    integer, allocatable   :: son(:,:)
    integer :: ind
    allocate(son(ngrida, twotondim))
    do ind = 1, twotondim
      read(iu) son(:, ind)
    end do
    ! A cell is a leaf if son = 0
    nleaf_total = nleaf_total + count(son == 0)
    deallocate(son)
  end subroutine count_leaves_in_oct

  !=========================================================================
  ! Read all CPU files and fill the output arrays.
  !=========================================================================
  subroutine ramses_read_all_cpus(repository, snapnum, ncpu, &
      unit_l, unit_d, unit_t, unit_v, gamma_eos, boxlen_code, &
      xleaf, yleaf, zleaf, leaf_level, &
      nH_arr, T_arr, vx_arr, vy_arr, vz_arr, &
      nleaf_total)

    character(len=*), intent(in) :: repository
    integer,          intent(in) :: snapnum, ncpu
    real(kind=8),     intent(in) :: unit_l, unit_d, unit_t, unit_v, gamma_eos, boxlen_code
    real(wp), intent(out) :: xleaf(nleaf_total), yleaf(nleaf_total), zleaf(nleaf_total)
    integer,  intent(out) :: leaf_level(nleaf_total)
    real(wp), intent(out) :: nH_arr(nleaf_total), T_arr(nleaf_total)
    real(wp), intent(out) :: vx_arr(nleaf_total), vy_arr(nleaf_total), vz_arr(nleaf_total)
    integer,  intent(in)  :: nleaf_total

    character(len=512) :: amr_file, hydro_file
    integer :: iu_amr = 97, iu_hyd = 98
    integer :: icpu, ncpu_f, ndim, nx, ny, nz, nlevelmax, ngridmax, nboundary
    integer :: ngrid_current, nvar
    integer :: ilevel, j, ngrida, ind, ivar, il_out
    real(kind=8) :: boxlen_f, gamma_f

    integer,  allocatable :: ngridlevel(:,:), ngridfile(:,:), ngridbound(:,:)
    real(kind=8), allocatable :: xg(:,:), son_r(:,:)
    integer,      allocatable :: son(:,:)
    real(kind=8), allocatable :: var(:,:,:)
    real(kind=4), allocatable :: var_sp(:)
    real(kind=8), allocatable :: xc(:,:)  ! (twotondim, ndim) cell offsets
    integer :: twotondim
    real(kind=8) :: dx
    real(kind=8) :: xbound(3)
    real(kind=8) :: rho_cgs, dens_code, px, py, pz, eint, pressure_cgs
    real(kind=8) :: kinetic_code
    integer :: dens_ivar, vel_ivar(3), thermo_ivar
    character(len=16) :: velocity_layout, thermo_mode
    il_out = 0

    do icpu = 1, ncpu
      write(amr_file,  '(a,"/output_",i5.5,"/amr_",i5.5,".out",i5.5)') &
          trim(repository), snapnum, snapnum, icpu
      write(hydro_file,'(a,"/output_",i5.5,"/hydro_",i5.5,".out",i5.5)') &
          trim(repository), snapnum, snapnum, icpu

      ! --- Open and read AMR header ---
      open(iu_amr, file=trim(amr_file),   form='unformatted', status='old', action='read')
      open(iu_hyd, file=trim(hydro_file), form='unformatted', status='old', action='read')

      read(iu_amr) ncpu_f
      read(iu_amr) ndim
      read(iu_amr) nx, ny, nz
      read(iu_amr) nlevelmax
      read(iu_amr) ngridmax
      read(iu_amr) nboundary
      read(iu_amr) ngrid_current
      read(iu_amr) boxlen_f
      call skip_records(iu_amr, 13)

      twotondim = 2**ndim
      xbound = [dble(nx/2), dble(ny/2), dble(nz/2)]

      allocate(ngridlevel(ncpu_f, nlevelmax))
      allocate(ngridfile(ncpu_f + max(nboundary,1), nlevelmax))
      if (nboundary > 0) allocate(ngridbound(nboundary, nlevelmax))

      read(iu_amr) ngridlevel
      ngridfile(1:ncpu_f, :) = ngridlevel
      read(iu_amr)
      if (nboundary > 0) then
        read(iu_amr); read(iu_amr)
        read(iu_amr) ngridbound
        ngridfile(ncpu_f+1:ncpu_f+nboundary, :) = ngridbound
        deallocate(ngridbound)
      end if
      call skip_records(iu_amr, 6)

      ! --- Hydro header ---
      read(iu_hyd)          ! ncpu
      read(iu_hyd) nvar
      read(iu_hyd)          ! ndim
      read(iu_hyd)          ! nlevelmax
      read(iu_hyd)          ! nboundary
      read(iu_hyd) gamma_f

      dens_ivar = ramses_density_var
      vel_ivar  = ramses_velocity_var
      thermo_ivar = ramses_thermo_var
      velocity_layout = ramses_velocity_layout
      thermo_mode     = ramses_thermo_mode

      if (dens_ivar < 1 .or. dens_ivar > nvar) dens_ivar = 1
      if (any(vel_ivar < 1) .or. any(vel_ivar > nvar)) vel_ivar = [2, 3, 4]
      if (thermo_ivar < 1 .or. thermo_ivar > nvar) then
        thermo_ivar = min(5, nvar)
        if (thermo_mode /= 'pressure' .and. thermo_mode /= 'energy') thermo_mode = 'constant'
      end if

      allocate(xc(twotondim, ndim))

      ! --- Loop over levels ---
      do ilevel = 1, nlevelmax
        dx = 0.5d0**ilevel
        ! Cell offsets within the oct (RAMSES convention)
        do ind = 1, twotondim
          xc(ind, 3) = dble((ind-1)/4)
          xc(ind, 2) = dble(mod((ind-1)/2, 2))
          xc(ind, 1) = dble(mod(ind-1, 2))
          xc(ind, :) = (xc(ind, :) - 0.5d0) * dx
        end do

        if (allocated(xg))    deallocate(xg)
        if (allocated(son))   deallocate(son)
        if (allocated(var))   deallocate(var)
        if (allocated(var_sp)) deallocate(var_sp)

        ngrida = ngridfile(icpu, ilevel)
        if (ngrida > 0) then
          allocate(xg(ngrida, ndim))
          allocate(son(ngrida, twotondim))
          allocate(var(ngrida, twotondim, nvar))
          if (ramses_hydro_prec == 4) allocate(var_sp(ngrida))
        end if

        do j = 1, nboundary + ncpu_f
          if (ngridfile(j, ilevel) > 0) then
            call skip_records(iu_amr, 3)         ! grid/next/prev
            do ivar = 1, ndim
              if (j == icpu) then
                read(iu_amr) xg(:, ivar)
              else
                call skip_records(iu_amr, 1)
              end if
            end do
            call skip_records(iu_amr, 1)         ! father
            call skip_records(iu_amr, 2*ndim)    ! nbor
            do ind = 1, twotondim
              if (j == icpu) then
                read(iu_amr) son(:, ind)
              else
                call skip_records(iu_amr, 1)
              end if
            end do
            call skip_records(iu_amr, 2*twotondim) ! cpu_map + ref_map
          end if

          ! Hydro: skip level/domain headers
          call skip_records(iu_hyd, 2)
          if (ngridfile(j, ilevel) > 0) then
            do ind = 1, twotondim
              do ivar = 1, nvar
                if (j == icpu) then
                  if (ramses_hydro_prec == 4) then
                    read(iu_hyd) var_sp
                    var(:, ind, ivar) = dble(var_sp)
                  else
                    read(iu_hyd) var(:, ind, ivar)
                  end if
                else
                  call skip_records(iu_hyd, 1)
                end if
              end do
            end do
          end if
        end do

        ! Collect leaf cells for this cpu at this level
        if (ngrida > 0) then
          do ind = 1, twotondim
            do j = 1, ngrida
              if (son(j, ind) == 0) then
                ! This is a leaf cell
                il_out = il_out + 1
                if (il_out > nleaf_total) cycle  ! safety: shouldn't happen

                ! Position in code units (fraction of boxlen)
                xleaf(il_out) = (xg(j, 1) + xc(ind, 1) - xbound(1)) / dble(nx) + 0.5d0
                yleaf(il_out) = (xg(j, 2) + xc(ind, 2) - xbound(2)) / dble(ny) + 0.5d0
                zleaf(il_out) = (xg(j, 3) + xc(ind, 3) - xbound(3)) / dble(nz) + 0.5d0
                leaf_level(il_out) = ilevel

                ! Physical quantities
                ! var(:,:,1) = mass density [code units]
                ! var(:,:,2:4) = momentum density or velocity [code units]
                ! var(:,:,5) = total energy or thermal pressure [code units]
                dens_code = var(j, ind, dens_ivar)
                rho_cgs   = dens_code * unit_d
                nH_arr(il_out) = real(rho_cgs / massH_cgs, wp)

                kinetic_code = 0.0d0
                if (trim(velocity_layout) == 'velocity') then
                  px = var(j, ind, vel_ivar(1)) * unit_v
                  py = var(j, ind, vel_ivar(2)) * unit_v
                  pz = var(j, ind, vel_ivar(3)) * unit_v
                  kinetic_code = 0.5d0 * dens_code * ( &
                      var(j,ind,vel_ivar(1))**2 + &
                      var(j,ind,vel_ivar(2))**2 + &
                      var(j,ind,vel_ivar(3))**2 )
                else
                  ! Standard RAMSES: stored as momentum density rho*v
                  if (dens_code > 0.0d0) then
                    px = var(j, ind, vel_ivar(1)) / dens_code * unit_v
                    py = var(j, ind, vel_ivar(2)) / dens_code * unit_v
                    pz = var(j, ind, vel_ivar(3)) / dens_code * unit_v
                    kinetic_code = 0.5d0 * dens_code * ( &
                        (var(j,ind,vel_ivar(1))/max(dens_code,1d-40))**2 + &
                        (var(j,ind,vel_ivar(2))/max(dens_code,1d-40))**2 + &
                        (var(j,ind,vel_ivar(3))/max(dens_code,1d-40))**2 )
                  else
                    px = 0.0d0;  py = 0.0d0;  pz = 0.0d0
                  end if
                end if
                vx_arr(il_out) = real(px, wp)
                vy_arr(il_out) = real(py, wp)
                vz_arr(il_out) = real(pz, wp)

                if (trim(thermo_mode) == 'pressure') then
                  pressure_cgs = var(j, ind, thermo_ivar) * unit_d * unit_v**2
                  if (rho_cgs > 0.0d0) then
                    T_arr(il_out) = real(pressure_cgs * mu_neutral * massH_cgs / &
                        max(rho_cgs, 1.0d-40) / boltzmann_cgs, wp)
                  else
                    T_arr(il_out) = 10.0_wp
                  end if
                else if (trim(thermo_mode) == 'energy') then
                  eint = (var(j, ind, thermo_ivar) - kinetic_code) / max(dens_code, 1d-40)
                  eint = max(eint, 0.0d0)
                  T_arr(il_out) = real((gamma_f - 1.0d0) * eint * unit_v**2 * &
                      mu_neutral * massH_cgs / boltzmann_cgs, wp)
                else
                  T_arr(il_out) = 1.0e4_wp
                end if
                T_arr(il_out) = max(T_arr(il_out), 10.0_wp)
              end if
            end do
          end do
        end if
      end do  ! ilevel

      deallocate(ngridlevel, ngridfile, xc)
      if (allocated(xg))    deallocate(xg)
      if (allocated(son))   deallocate(son)
      if (allocated(var))   deallocate(var)
      if (allocated(var_sp)) deallocate(var_sp)
      close(iu_amr)
      close(iu_hyd)
    end do  ! icpu
  end subroutine ramses_read_all_cpus

  !=========================================================================
  ! Placeholder interface for arbitrary AMR input formats.
  ! Users supply their own routine that fills these arrays.
  !
  ! The caller must fill:
  !   xleaf, yleaf, zleaf  [cm or kpc]
  !   leaf_level
  !   nH_cgs [cm^-3], T_cgs [K], vx/vy/vz_cgs [km/s]
  !
  ! Then call grid_create_amr to build the octree from these arrays.
  !=========================================================================
  subroutine generic_amr_read(filename, &
      xleaf, yleaf, zleaf, leaf_level, &
      nH_cgs, T_cgs, vx_cgs, vy_cgs, vz_cgs, &
      nleaf, boxlen_phys)
    character(len=*), intent(in)       :: filename
    real(wp), allocatable, intent(out) :: xleaf(:), yleaf(:), zleaf(:)
    integer,  allocatable, intent(out) :: leaf_level(:)
    real(wp), allocatable, intent(out) :: nH_cgs(:), T_cgs(:)
    real(wp), allocatable, intent(out) :: vx_cgs(:), vy_cgs(:), vz_cgs(:)
    integer,               intent(out) :: nleaf
    real(wp),              intent(out) :: boxlen_phys

    ! Reads either:
    !   1) a simple text format
    !        header: nleaf  boxlen
    !        rows  : x, y, z, level, nH [cm^-3], T [K], vx [km/s], vy [km/s], vz [km/s]
    !   2) an AMRGrid FITS binary table when filename ends with .fits or .fits.gz
    !
    ! For FITS input, this routine reads the AMRGrid binary table by column index
    ! (1:x, 2:y, 3:z, 4:level, 5:gas density, 6:T, 7:vx, 8:vy, 9:vz), rather
    ! than depending on column names.
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

    if (filename_is_fits(filename)) then
      call generic_amr_read_fits(filename, &
          xleaf, yleaf, zleaf, leaf_level, &
          nH_cgs, T_cgs, vx_cgs, vy_cgs, vz_cgs, &
          nleaf, boxlen_phys)
      return
    end if

    unit = 51
    open(unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'generic_amr_read: cannot open file'

    read(unit, *) n, bl
    nleaf       = n
    boxlen_phys = bl          ! in code units (kpc, pc, au, or cm) as written in the file

    allocate(xleaf(n), yleaf(n), zleaf(n), leaf_level(n))
    allocate(nH_cgs(n), T_cgs(n), vx_cgs(n), vy_cgs(n), vz_cgs(n))

    do n = 1, nleaf
      read(unit, *, iostat=ios) x, y, z, lv, nH, T, vx, vy, vz
      if (ios /= 0) exit
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
    close(unit)
  end subroutine generic_amr_read

  subroutine generic_amr_read_fits(filename, &
      xleaf, yleaf, zleaf, leaf_level, &
      nH_cgs, T_cgs, vx_cgs, vy_cgs, vz_cgs, &
      nleaf, boxlen_phys)
    character(len=*), intent(in)       :: filename
    real(wp), allocatable, intent(out) :: xleaf(:), yleaf(:), zleaf(:)
    integer,  allocatable, intent(out) :: leaf_level(:)
    real(wp), allocatable, intent(out) :: nH_cgs(:), T_cgs(:)
    real(wp), allocatable, intent(out) :: vx_cgs(:), vy_cgs(:), vz_cgs(:)
    integer,               intent(out) :: nleaf
    real(wp),              intent(out) :: boxlen_phys

    integer :: unit, status
    integer(int32) :: nleaf_i4
    integer(int32), allocatable :: leaf_level_i4(:)
    real(wp) :: origin_x, origin_y, origin_z

    status = 0
    call fits_open_old(unit, trim(filename), status)
    if (status /= 0) stop 'generic_amr_read_fits: cannot open FITS file'

    call fits_move_to_next_hdu(unit, status)
    if (status /= 0) stop 'generic_amr_read_fits: cannot move to binary table HDU'

    call fits_get_keyword(unit, 'NAXIS2', nleaf_i4, status)
    if (status /= 0) stop 'generic_amr_read_fits: cannot read NAXIS2'

    boxlen_phys = 0.0_wp
    origin_x    = 0.0_wp
    origin_y    = 0.0_wp
    origin_z    = 0.0_wp
    call fits_get_keyword(unit, 'BOXLEN',  boxlen_phys, status)
    if (status /= 0) stop 'generic_amr_read_fits: cannot read BOXLEN'

    status = 0
    call fits_get_keyword(unit, 'ORIGINX', origin_x, status)
    if (status /= 0) status = 0
    call fits_get_keyword(unit, 'ORIGINY', origin_y, status)
    if (status /= 0) status = 0
    call fits_get_keyword(unit, 'ORIGINZ', origin_z, status)
    if (status /= 0) status = 0

    nleaf = int(nleaf_i4)
    allocate(xleaf(nleaf), yleaf(nleaf), zleaf(nleaf), leaf_level(nleaf))
    allocate(nH_cgs(nleaf), T_cgs(nleaf), vx_cgs(nleaf), vy_cgs(nleaf), vz_cgs(nleaf))
    allocate(leaf_level_i4(nleaf))

    call fits_read_table_column(unit, 1, xleaf,         status)
    call fits_read_table_column(unit, 2, yleaf,         status)
    call fits_read_table_column(unit, 3, zleaf,         status)
    call fits_read_table_column(unit, 4, leaf_level_i4, status)
    call fits_read_table_column(unit, 5, nH_cgs,        status)
    call fits_read_table_column(unit, 6, T_cgs,         status)
    call fits_read_table_column(unit, 7, vx_cgs,        status)
    call fits_read_table_column(unit, 8, vy_cgs,        status)
    call fits_read_table_column(unit, 9, vz_cgs,        status)
    if (status /= 0) stop 'generic_amr_read_fits: cannot read table columns'

    call fits_close(unit, status)

    xleaf(:) = xleaf(:) - origin_x
    yleaf(:) = yleaf(:) - origin_y
    zleaf(:) = zleaf(:) - origin_z
    leaf_level(:) = int(leaf_level_i4(:))

    deallocate(leaf_level_i4)
  end subroutine generic_amr_read_fits

  subroutine ramses_detect_layout(repository, snapnum)
    character(len=*), intent(in) :: repository
    integer,          intent(in) :: snapnum

    character(len=512) :: filename, line
    character(len=128) :: name, names(256)
    integer :: ios, unit, ivar, ipos1, ipos2

    ramses_descriptor_found = .false.
    ramses_density_var = 1
    ramses_velocity_var = [2, 3, 4]
    ramses_thermo_var = 5
    ramses_velocity_layout = 'momentum'
    ramses_thermo_mode = 'energy'
    names = ''

    write(filename, '(a,"/output_",i5.5,"/hydro_file_descriptor.txt")') &
        trim(repository), snapnum

    unit = 52
    open(unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) return

    ramses_descriptor_found = .true.
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      ipos1 = index(line, '#')
      ipos2 = index(line, ':')
      if (ipos1 <= 0 .or. ipos2 <= ipos1) cycle
      read(line(ipos1+1:ipos2-1), *, iostat=ios) ivar
      if (ios /= 0) cycle
      if (ivar < 1 .or. ivar > size(names)) cycle
      name = str_lower(adjustl(line(ipos2+1:)))
      names(ivar) = trim(name)
    end do
    close(unit)

    do ivar = 1, size(names)
      if (trim(names(ivar)) == 'density') then
        ramses_density_var = ivar
        exit
      end if
    end do

    if (index(names(2), 'velocity_x') > 0 .and. index(names(3), 'velocity_y') > 0 .and. &
        index(names(4), 'velocity_z') > 0) then
      ramses_velocity_layout = 'velocity'
      ramses_velocity_var = [2, 3, 4]
    else if (index(names(2), 'momentum_x') > 0 .and. index(names(3), 'momentum_y') > 0 .and. &
             index(names(4), 'momentum_z') > 0) then
      ramses_velocity_layout = 'momentum'
      ramses_velocity_var = [2, 3, 4]
    end if

    do ivar = 1, size(names)
      if (trim(names(ivar)) == 'thermal_pressure') then
        ramses_thermo_mode = 'pressure'
        ramses_thermo_var = ivar
        return
      end if
    end do

    do ivar = 1, size(names)
      select case (trim(names(ivar)))
      case ('pressure', 'gas_pressure', 'thermal_energy_density')
        ramses_thermo_mode = 'pressure'
        ramses_thermo_var = ivar
        return
      case ('total_energy', 'energy')
        ramses_thermo_mode = 'energy'
        ramses_thermo_var = ivar
        return
      end select
    end do
  end subroutine ramses_detect_layout

  pure function str_lower(str) result(out)
    character(len=*), intent(in) :: str
    character(len=len(str)) :: out
    integer :: i, code

    out = str
    do i = 1, len(str)
      code = iachar(out(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        out(i:i) = achar(code + 32)
      end if
    end do
  end function str_lower

  logical function filename_is_fits(filename)
    character(len=*), intent(in) :: filename
    character(len=:), allocatable :: name
    integer :: n

    name = trim(filename)
    n = len(name)
    filename_is_fits = .false.
    if (n >= 5) then
      if (name(n-4:n) == '.fits' .or. name(n-4:n) == '.FITS') then
        filename_is_fits = .true.
        return
      end if
    end if
    if (n >= 8) then
      if (name(n-7:n) == '.fits.gz' .or. name(n-7:n) == '.FITS.GZ') then
        filename_is_fits = .true.
        return
      end if
    end if
  end function filename_is_fits

  !=========================================================================
  ! Skip n Fortran unformatted records.
  !=========================================================================
  subroutine skip_records(unit, n)
    integer, intent(in) :: unit, n
    integer :: i
    do i = 1, n
      read(unit)
    end do
  end subroutine skip_records

end module read_ramses_amr_mod
