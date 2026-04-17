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
  !   - Standard conservative hydro variables (density, rho*v, energy/pressure, metallicity)
  !   - Unit scales read from output_{snapnum}/info_{snapnum}.txt
  !   - nHI computed from nH + T using either CIE or supplied neutral fraction
  !-----------------------------------------------------------------------
  use define
  implicit none
  private

  public :: ramses_read_leaf_cells
  public :: generic_amr_read

  ! Precision flag for RAMSES hydro output (4 = single, 8 = double)
  integer :: ramses_hydro_prec = 8

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
      call skip_records(iu, 5)

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
    real(kind=8) :: dx, massH_cgs
    real(kind=8) :: xbound(3)
    real(kind=8) :: rho_cgs, dens_code, px, py, pz, eint

    massH_cgs = 1.6726e-24_8
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
      call skip_records(iu_amr, 5)

      ! --- Hydro header ---
      read(iu_hyd)          ! ncpu
      read(iu_hyd) nvar
      read(iu_hyd)          ! ndim
      read(iu_hyd)          ! nlevelmax
      read(iu_hyd)          ! nboundary
      read(iu_hyd) gamma_f

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
                dens_code = var(j, ind, 1)
                rho_cgs   = dens_code * unit_d
                nH_arr(il_out) = real(rho_cgs / massH_cgs, wp)

                ! Velocity: check if RAMSES stored momentum or velocity
                ! Standard RAMSES: stored as momentum density ρ*v
                if (dens_code > 0.0d0) then
                  px = var(j, ind, 2) / dens_code * unit_v
                  py = var(j, ind, 3) / dens_code * unit_v
                  pz = var(j, ind, 4) / dens_code * unit_v
                else
                  px = 0.0d0;  py = 0.0d0;  pz = 0.0d0
                end if
                vx_arr(il_out) = real(px, wp)
                vy_arr(il_out) = real(py, wp)
                vz_arr(il_out) = real(pz, wp)

                ! Temperature from internal energy or thermal pressure
                ! Standard RAMSES: var5 = total energy density
                ! T = (gamma-1) * eint / (rho/mu/mH) * mH/kB  — use unit conversion
                ! A simpler approximation: T = (gamma-1)*eint_per_mass * mu*mH/kB
                ! We use unit_T2 = (unit_v)^2 * mu * mH / kB  where mu ~ 1.22 (fully neutral)
                if (nvar >= 5) then
                  eint = (var(j, ind, 5) - 0.5d0 * dens_code * (                &
                          (var(j,ind,2)/max(dens_code,1d-40))**2 +              &
                          (var(j,ind,3)/max(dens_code,1d-40))**2 +              &
                          (var(j,ind,4)/max(dens_code,1d-40))**2)) / dens_code
                  eint = max(eint, 0.0d0)
                  ! T in K: (gamma-1)*eint[erg/g] * mu*mH/kB  (mu*mH/kB in K/(erg/g))
                  ! unit_v^2 has units [cm^2/s^2 = erg/g]; kB = 1.381e-16 erg/K
                  T_arr(il_out) = real((gamma_f - 1.0d0) * eint * unit_v**2 * &
                      1.22d0 * 1.6726d-24 / 1.381d-16, wp)
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

    ! Reads a simple text format.
    ! Format (one leaf per line):
    !   x, y, z, level, nH [cm^-3], T [K], vx [km/s], vy [km/s], vz [km/s]
    ! Header line: nleaf  boxlen
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
