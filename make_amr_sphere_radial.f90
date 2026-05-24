!=========================================================================
! make_amr_sphere_radial.f90
!
! Standalone Fortran program to generate a generic AMR grid file for a
! spherically symmetric medium with radial shell refinement.
!
! This is the Fortran equivalent of python/AMR_grid/make_amr_sphere_radial.py.
! It avoids the Python h5py dependency and writes HDF5 or FITS files
! directly via the LaRT iofile_mod facade.
!
! Build:
!   make make_amr_sphere     (see Makefile target)
!
! Usage:
!   ./make_amr_sphere_radial.x [options]
!
! Options (all optional, defaults shown):
!   --boxlen 2.0            Box side length (code units)
!   --rmax   1.0            Sphere radius  (code units)
!   --rmin   0.0            Inner cavity radius
!   --level_min 5           Outermost shell refinement level
!   --level_max 10          Innermost shell (center) refinement level
!   --spacing log           Shell spacing: 'linear' or 'log'
!   --ratio   0.5           Geometric ratio for log spacing
!   --density uniform       Density profile: uniform|gaussian|exponential|power_law
!   --n0      1.0           Reference number density [cm^-3]
!   --sigma   0.0           Gaussian sigma (default 0.4*rmax)
!   --r_scale 0.0           Exponential scale radius (default 0.3*rmax)
!   --density_alpha 0.0     Power-law exponent: n(r) = n0*(rmax/r)^alpha
!   --temperature 1e4       Gas temperature [K]
!   --velocity none         Velocity law: none|hubble|constant_radial|
!                           power_law|linear_decelerate
!   --v_exp   0.0           Expansion speed [km/s]
!   --v_power 1.0           Power-law velocity exponent
!   --refine_boundary       Refine cells at sphere surface
!   --boundary_level_max -1 Boundary level (default = level_max)
!   -o output.h5            Output file (.h5, .fits.gz, or .dat)
!=========================================================================
program make_amr_sphere_radial
  use define, only: wp, deg2rad
  use iofile_mod
  use, intrinsic :: iso_fortran_env, only: int32
  implicit none

  ! Parameters
  real(wp) :: boxlen       = 2.0_wp
  real(wp) :: rmax         = 1.0_wp
  real(wp) :: rmin_val     = 0.0_wp
  integer  :: level_min    = 5
  integer  :: level_max    = 10
  character(len=32) :: spacing = 'log'
  real(wp) :: ratio        = 0.5_wp
  character(len=32) :: density_type = 'uniform'
  real(wp) :: n0           = 1.0_wp
  real(wp) :: sigma        = 0.0_wp    ! 0 -> auto = 0.4*rmax
  real(wp) :: r_scale      = 0.0_wp    ! 0 -> auto = 0.3*rmax
  real(wp) :: density_alpha = 0.0_wp
  real(wp) :: temperature  = 1.0e4_wp
  character(len=32) :: velocity_law = 'none'
  real(wp) :: v_exp        = 0.0_wp
  real(wp) :: v_power      = 1.0_wp
  real(wp) :: cone_opening  = 0.0_wp    ! bicone half-opening angle [deg]; 0 = full sphere
  logical  :: refine_boundary = .false.
  integer  :: boundary_level_max_in = -1   ! -1 -> use level_max
  character(len=512) :: output_file = 'amr_sphere.h5'

  ! Shell radii
  integer  :: n_shells, boundary_lmax
  real(wp), allocatable :: shell_radii(:)

  ! Leaf cell storage (dynamically grown)
  integer  :: nleaf, nleaf_alloc
  real(wp), allocatable :: xleaf(:), yleaf(:), zleaf(:)
  real(wp), allocatable :: dens(:), Tgas(:), vx(:), vy(:), vz(:)
  integer(int32), allocatable :: leaf_level(:)

  ! I/O
  type(io_file_type) :: iofh
  integer :: status

  ! Parse command-line arguments
  call parse_args()

  ! Set defaults
  if (sigma <= 0.0_wp)   sigma   = 0.4_wp * rmax
  if (r_scale <= 0.0_wp) r_scale = 0.3_wp * rmax
  boundary_lmax = boundary_level_max_in
  if (boundary_lmax < 0) boundary_lmax = level_max

  ! Build shell radii
  n_shells = level_max - level_min + 1
  allocate(shell_radii(n_shells))
  call build_shell_radii(rmax, n_shells, spacing, ratio, shell_radii)

  ! Generate leaf cells by recursive subdivision
  nleaf = 0
  nleaf_alloc = 100000
  allocate(xleaf(nleaf_alloc), yleaf(nleaf_alloc), zleaf(nleaf_alloc))
  allocate(dens(nleaf_alloc), Tgas(nleaf_alloc))
  allocate(vx(nleaf_alloc), vy(nleaf_alloc), vz(nleaf_alloc))
  allocate(leaf_level(nleaf_alloc))

  call subdivide_box()

  ! Print summary
  call print_summary()

  ! Write output
  status = 0
  call io_open_new(iofh, trim(output_file), status)
  if (status /= 0) stop 'make_amr_sphere_radial: cannot create output file'

  call io_write_table_column(iofh, 'x',      xleaf(1:nleaf),      status)
  call io_write_table_column(iofh, 'y',      yleaf(1:nleaf),      status)
  call io_write_table_column(iofh, 'z',      zleaf(1:nleaf),      status)
  call io_write_table_column(iofh, 'level',  leaf_level(1:nleaf), status)
  call io_write_table_column(iofh, 'gasDen', dens(1:nleaf),       status)
  call io_write_table_column(iofh, 'T',      Tgas(1:nleaf),       status)
  call io_write_table_column(iofh, 'vx',     vx(1:nleaf),         status)
  call io_write_table_column(iofh, 'vy',     vy(1:nleaf),         status)
  call io_write_table_column(iofh, 'vz',     vz(1:nleaf),         status)

  call io_put_keyword(iofh, 'BOXLEN',  boxlen,              'Simulation box length', status)
  call io_put_keyword(iofh, 'ORIGINX', -0.5_wp * boxlen,    'Box origin x',         status)
  call io_put_keyword(iofh, 'ORIGINY', -0.5_wp * boxlen,    'Box origin y',         status)
  call io_put_keyword(iofh, 'ORIGINZ', -0.5_wp * boxlen,    'Box origin z',         status)
  call io_put_keyword(iofh, 'NLEAF',   int(nleaf, int32),   'Number of leaf cells', status)

  call io_close(iofh, status)
  if (status /= 0) stop 'make_amr_sphere_radial: error closing output file'

  write(*,'(a,i0,a,a)') 'Written ', nleaf, ' leaf cells to ', trim(output_file)

contains

  !-----------------------------------------------------------------------
  ! Build shell boundary radii (outermost first)
  !-----------------------------------------------------------------------
  subroutine build_shell_radii(rmax, n, sp, rat, radii)
    real(wp), intent(in)  :: rmax, rat
    integer,  intent(in)  :: n
    character(len=*), intent(in) :: sp
    real(wp), intent(out) :: radii(n)
    integer :: i

    if (trim(sp) == 'linear') then
      do i = 1, n
        radii(i) = rmax * real(n - i + 1, wp) / real(n, wp)
      end do
    else  ! log
      do i = 1, n
        radii(i) = rmax * rat ** (i - 1)
      end do
    end if
  end subroutine

  !-----------------------------------------------------------------------
  ! Determine target refinement level from distance to center
  !-----------------------------------------------------------------------
  integer function target_level(d)
    real(wp), intent(in) :: d
    integer :: i
    target_level = 0
    do i = 1, n_shells
      if (d < shell_radii(i)) target_level = level_min + i - 1
    end do
  end function

  !-----------------------------------------------------------------------
  ! Recursive cell subdivision
  !-----------------------------------------------------------------------
  subroutine subdivide_box()
    real(wp) :: half
    half = 0.5_wp * boxlen
    call subdivide_cell(0.0_wp, 0.0_wp, 0.0_wp, half, 0)
  end subroutine

  recursive subroutine subdivide_cell(cx, cy, cz, ch, lev)
    real(wp), intent(in) :: cx, cy, cz, ch  ! cell center and half-width
    integer,  intent(in) :: lev
    real(wp) :: d_center, d_close, dx, dy, dz
    integer  :: t_cen, t_cls, target_lev
    real(wp) :: ch2, ox, oy, oz
    integer  :: ix, iy, iz

    ! Distance from cell center to origin (sphere center)
    d_center = sqrt(cx*cx + cy*cy + cz*cz)
    t_cen = target_level(d_center)

    ! Distance from closest point in cell to origin
    dx = max(abs(cx) - ch, 0.0_wp)
    dy = max(abs(cy) - ch, 0.0_wp)
    dz = max(abs(cz) - ch, 0.0_wp)
    d_close = sqrt(dx*dx + dy*dy + dz*dz)
    t_cls = target_level(d_close)

    target_lev = max(t_cen, t_cls)

    ! Check boundary refinement (sphere surface + cone boundary)
    if (refine_boundary .and. lev < boundary_lmax) then
      if (cell_intersects_sphere(cx, cy, cz, ch, rmax)) then
        target_lev = max(target_lev, boundary_lmax)
      end if
      if (cone_opening > 0.0_wp .and. cone_opening < 90.0_wp) then
        if (cell_intersects_cone(cx, cy, cz, ch, cos(cone_opening * deg2rad))) then
          target_lev = max(target_lev, boundary_lmax)
        end if
      end if
    end if

    ! If this cell needs further refinement, subdivide
    if (target_lev > lev .and. lev < level_max) then
      ch2 = ch * 0.5_wp
      do iz = 0, 1
        do iy = 0, 1
          do ix = 0, 1
            ox = cx + ch2 * (2*ix - 1)
            oy = cy + ch2 * (2*iy - 1)
            oz = cz + ch2 * (2*iz - 1)
            call subdivide_cell(ox, oy, oz, ch2, lev + 1)
          end do
        end do
      end do
      return
    end if

    ! This cell is a leaf — add it
    call add_leaf(cx, cy, cz, ch, lev)
  end subroutine

  !-----------------------------------------------------------------------
  ! Check if a cell intersects the sphere surface at radius r
  !-----------------------------------------------------------------------
  logical function cell_intersects_sphere(cx, cy, cz, ch, r)
    real(wp), intent(in) :: cx, cy, cz, ch, r
    real(wp) :: d_close, d_far, dx, dy, dz
    ! Closest point in cell to origin
    dx = max(abs(cx) - ch, 0.0_wp)
    dy = max(abs(cy) - ch, 0.0_wp)
    dz = max(abs(cz) - ch, 0.0_wp)
    d_close = sqrt(dx*dx + dy*dy + dz*dz)
    ! Farthest point in cell from origin
    dx = abs(cx) + ch
    dy = abs(cy) + ch
    dz = abs(cz) + ch
    d_far = sqrt(dx*dx + dy*dy + dz*dz)
    ! Cell intersects the sphere surface if d_close <= r <= d_far
    cell_intersects_sphere = (d_close <= r .and. d_far >= r)
  end function

  !-----------------------------------------------------------------------
  ! Check if a cell intersects the cone boundary surface.
  ! The cone boundary is |cos(theta)| = cos(cone_opening), i.e. the set
  ! of points where abs(z)/r = cos_cone.  A cell intersects this surface
  ! if the 8 corners span both inside and outside the cone.
  !-----------------------------------------------------------------------
  logical function cell_intersects_cone(cx, cy, cz, ch, cos_c)
    real(wp), intent(in) :: cx, cy, cz, ch, cos_c
    real(wp) :: x0, y0, z0, r2, cos_th
    integer  :: ix, iy, iz, n_in, n_out

    n_in  = 0
    n_out = 0
    do iz = 0, 1
    do iy = 0, 1
    do ix = 0, 1
      x0 = cx + ch * (2*ix - 1)
      y0 = cy + ch * (2*iy - 1)
      z0 = cz + ch * (2*iz - 1)
      r2 = x0*x0 + y0*y0 + z0*z0
      if (r2 > 0.0_wp) then
        cos_th = abs(z0) / sqrt(r2)
        if (cos_th >= cos_c) then
          n_in = n_in + 1
        else
          n_out = n_out + 1
        end if
      else
        n_in = n_in + 1  ! origin is inside cone
      end if
    end do
    end do
    end do
    ! Cell intersects cone boundary if some corners are inside and some outside
    cell_intersects_cone = (n_in > 0 .and. n_out > 0)
  end function

  !-----------------------------------------------------------------------
  ! Add a leaf cell with computed physical quantities
  !-----------------------------------------------------------------------
  subroutine add_leaf(cx, cy, cz, ch, lev)
    real(wp), intent(in) :: cx, cy, cz, ch
    integer,  intent(in) :: lev
    real(wp) :: r, d, vr, scale

    ! Grow arrays if needed
    if (nleaf >= nleaf_alloc) call grow_arrays()

    nleaf = nleaf + 1
    xleaf(nleaf) = cx
    yleaf(nleaf) = cy
    zleaf(nleaf) = cz
    leaf_level(nleaf) = int(lev, int32)
    Tgas(nleaf) = temperature

    r = sqrt(cx*cx + cy*cy + cz*cz)

    ! Density (with cone mask)
    if (r > rmax .or. r < rmin_val) then
      dens(nleaf) = 0.0_wp
    else if (cone_opening > 0.0_wp .and. cone_opening < 90.0_wp .and. r > 0.0_wp) then
      if (abs(cz) / r < cos(cone_opening * deg2rad)) then
        dens(nleaf) = 0.0_wp
      else
        dens(nleaf) = compute_density(r)
      end if
    else
      dens(nleaf) = compute_density(r)
    end if

    ! Velocity (zero where density is zero)
    if (dens(nleaf) > 0.0_wp) then
      call compute_velocity(cx, cy, cz, r, vx(nleaf), vy(nleaf), vz(nleaf))
    else
      vx(nleaf) = 0.0_wp; vy(nleaf) = 0.0_wp; vz(nleaf) = 0.0_wp
    end if
  end subroutine

  !-----------------------------------------------------------------------
  ! Compute density at radius r
  !-----------------------------------------------------------------------
  real(wp) function compute_density(r)
    real(wp), intent(in) :: r
    select case (trim(density_type))
    case ('uniform')
      compute_density = n0
    case ('gaussian')
      compute_density = n0 * exp(-r*r / (2.0_wp * sigma*sigma))
    case ('exponential')
      compute_density = n0 * exp(-r / r_scale)
    case ('power_law')
      if (r > 0.0_wp) then
        compute_density = n0 * (rmax / r) ** density_alpha
      else
        compute_density = n0
      end if
    case default
      compute_density = n0
    end select
  end function

  !-----------------------------------------------------------------------
  ! Compute velocity at position (x,y,z) with radius r
  !-----------------------------------------------------------------------
  subroutine compute_velocity(x, y, z, r, vxo, vyo, vzo)
    real(wp), intent(in)  :: x, y, z, r
    real(wp), intent(out) :: vxo, vyo, vzo
    real(wp) :: vr, scale, denom

    vxo = 0.0_wp; vyo = 0.0_wp; vzo = 0.0_wp
    if (r <= 0.0_wp .or. r > rmax .or. r < rmin_val) return

    select case (trim(velocity_law))
    case ('none')
      return
    case ('hubble')
      scale = v_exp / rmax
      vxo = scale * x; vyo = scale * y; vzo = scale * z
    case ('constant_radial')
      scale = v_exp / r
      vxo = scale * x; vyo = scale * y; vzo = scale * z
    case ('power_law')
      vr = v_exp * (r / rmax) ** v_power
      scale = vr / r
      vxo = scale * x; vyo = scale * y; vzo = scale * z
    case ('linear_decelerate')
      denom = rmax - rmin_val
      if (denom <= 0.0_wp) denom = 1.0_wp
      vr = v_exp * (rmax - r) / denom
      scale = vr / r
      vxo = scale * x; vyo = scale * y; vzo = scale * z
    end select
  end subroutine

  !-----------------------------------------------------------------------
  ! Grow dynamic arrays (double capacity)
  !-----------------------------------------------------------------------
  subroutine grow_arrays()
    integer :: new_size
    new_size = nleaf_alloc * 2
    call grow_real(xleaf, nleaf_alloc, new_size)
    call grow_real(yleaf, nleaf_alloc, new_size)
    call grow_real(zleaf, nleaf_alloc, new_size)
    call grow_real(dens,  nleaf_alloc, new_size)
    call grow_real(Tgas,  nleaf_alloc, new_size)
    call grow_real(vx,    nleaf_alloc, new_size)
    call grow_real(vy,    nleaf_alloc, new_size)
    call grow_real(vz,    nleaf_alloc, new_size)
    call grow_int4(leaf_level, nleaf_alloc, new_size)
    nleaf_alloc = new_size
  end subroutine

  subroutine grow_real(arr, old_size, new_size)
    real(wp), allocatable, intent(inout) :: arr(:)
    integer, intent(in) :: old_size, new_size
    real(wp), allocatable :: tmp(:)
    allocate(tmp(new_size))
    tmp(1:old_size) = arr(1:old_size)
    call move_alloc(tmp, arr)
  end subroutine

  subroutine grow_int4(arr, old_size, new_size)
    integer(int32), allocatable, intent(inout) :: arr(:)
    integer, intent(in) :: old_size, new_size
    integer(int32), allocatable :: tmp(:)
    allocate(tmp(new_size))
    tmp(1:old_size) = arr(1:old_size)
    call move_alloc(tmp, arr)
  end subroutine

  !-----------------------------------------------------------------------
  ! Print summary
  !-----------------------------------------------------------------------
  subroutine print_summary()
    integer :: i
    real(wp) :: r_in

    write(*,'(a)')       '============================================='
    write(*,'(a)')       '  make_amr_sphere_radial (Fortran)'
    write(*,'(a)')       '============================================='
    write(*,'(a,f10.4)') '  boxlen          = ', boxlen
    write(*,'(a,f10.4)') '  rmax            = ', rmax
    write(*,'(a,f10.4)') '  rmin            = ', rmin_val
    write(*,'(a,i6)')    '  level_min       = ', level_min
    write(*,'(a,i6)')    '  level_max       = ', level_max
    write(*,'(a,a)')     '  spacing         = ', trim(spacing)
    if (trim(spacing) == 'log') write(*,'(a,f10.4)') '  ratio           = ', ratio
    write(*,'(a,a)')     '  density         = ', trim(density_type)
    write(*,'(a,es12.4)')'  n0              = ', n0
    if (trim(density_type) == 'power_law') &
      write(*,'(a,f10.4)') '  density_alpha   = ', density_alpha
    write(*,'(a,es12.4)')'  temperature     = ', temperature
    write(*,'(a,a)')     '  velocity        = ', trim(velocity_law)
    if (trim(velocity_law) /= 'none') &
      write(*,'(a,f10.2)') '  v_exp           = ', v_exp
    if (refine_boundary) then
      write(*,'(a,i6)')  '  boundary_lmax   = ', boundary_lmax
    end if
    if (cone_opening > 0.0_wp) then
      write(*,'(a,f10.2)') '  cone_opening    = ', cone_opening
    end if
    write(*,'(a,i0)')    '  nleaf           = ', nleaf
    write(*,'(a,a)')     '  output          = ', trim(output_file)

    write(*,'(a)')       ''
    write(*,'(a)')       '  Shell structure (outermost -> center):'
    write(*,'(a)')       '    Outer radius   Inner radius   Level'
    write(*,'(a)')       '  ----------------------------------------'
    do i = 1, n_shells
      if (i < n_shells) then
        r_in = shell_radii(i+1)
      else
        r_in = 0.0_wp
      end if
      write(*,'(4x,f12.4,2x,f12.4,2x,i5)') shell_radii(i), r_in, level_min + i - 1
    end do
    write(*,'(a)')       ''
  end subroutine

  !-----------------------------------------------------------------------
  ! Command-line argument parser
  !-----------------------------------------------------------------------
  subroutine parse_args()
    integer :: nargs, i
    character(len=256) :: arg, val

    nargs = command_argument_count()
    i = 1
    do while (i <= nargs)
      call get_command_argument(i, arg)
      select case (trim(arg))
      case ('--boxlen')
        call get_next_arg(i, val); read(val, *) boxlen
      case ('--rmax')
        call get_next_arg(i, val); read(val, *) rmax
      case ('--rmin')
        call get_next_arg(i, val); read(val, *) rmin_val
      case ('--level_min')
        call get_next_arg(i, val); read(val, *) level_min
      case ('--level_max')
        call get_next_arg(i, val); read(val, *) level_max
      case ('--spacing')
        call get_next_arg(i, val); spacing = trim(val)
      case ('--ratio')
        call get_next_arg(i, val); read(val, *) ratio
      case ('--density')
        call get_next_arg(i, val); density_type = trim(val)
      case ('--n0')
        call get_next_arg(i, val); read(val, *) n0
      case ('--sigma')
        call get_next_arg(i, val); read(val, *) sigma
      case ('--r_scale')
        call get_next_arg(i, val); read(val, *) r_scale
      case ('--density_alpha')
        call get_next_arg(i, val); read(val, *) density_alpha
      case ('--temperature')
        call get_next_arg(i, val); read(val, *) temperature
      case ('--velocity')
        call get_next_arg(i, val); velocity_law = trim(val)
      case ('--v_exp')
        call get_next_arg(i, val); read(val, *) v_exp
      case ('--v_power')
        call get_next_arg(i, val); read(val, *) v_power
      case ('--refine_boundary')
        refine_boundary = .true.
      case ('--boundary_level_max')
        call get_next_arg(i, val); read(val, *) boundary_level_max_in
      case ('--cone_opening')
        call get_next_arg(i, val); read(val, *) cone_opening
      case ('-o', '--output')
        call get_next_arg(i, val); output_file = trim(val)
      case ('-h', '--help')
        call print_usage(); stop
      case default
        write(*,'(a,a)') 'Unknown argument: ', trim(arg)
        call print_usage(); stop
      end select
      i = i + 1
    end do
  end subroutine

  subroutine get_next_arg(i, val)
    integer, intent(inout) :: i
    character(len=*), intent(out) :: val
    i = i + 1
    if (i > command_argument_count()) then
      write(*,'(a)') 'Error: missing value for argument'
      stop
    end if
    call get_command_argument(i, val)
  end subroutine

  subroutine print_usage()
    write(*,'(a)') 'Usage: make_amr_sphere_radial.x [options]'
    write(*,'(a)') ''
    write(*,'(a)') 'Options:'
    write(*,'(a)') '  --boxlen VAL          Box side length (default: 2.0)'
    write(*,'(a)') '  --rmax VAL            Sphere radius (default: 1.0)'
    write(*,'(a)') '  --rmin VAL            Inner cavity radius (default: 0.0)'
    write(*,'(a)') '  --level_min N         Outermost shell level (default: 5)'
    write(*,'(a)') '  --level_max N         Innermost shell level (default: 10)'
    write(*,'(a)') '  --spacing TYPE        linear or log (default: log)'
    write(*,'(a)') '  --ratio VAL           Log spacing ratio (default: 0.5)'
    write(*,'(a)') '  --density TYPE        uniform|gaussian|exponential|power_law'
    write(*,'(a)') '  --n0 VAL              Reference density (default: 1.0)'
    write(*,'(a)') '  --sigma VAL           Gaussian sigma (default: 0.4*rmax)'
    write(*,'(a)') '  --r_scale VAL         Exponential scale (default: 0.3*rmax)'
    write(*,'(a)') '  --density_alpha VAL   Power-law exponent (default: 0.0)'
    write(*,'(a)') '  --temperature VAL     Gas temperature [K] (default: 1e4)'
    write(*,'(a)') '  --velocity TYPE       none|hubble|constant_radial|'
    write(*,'(a)') '                        power_law|linear_decelerate'
    write(*,'(a)') '  --v_exp VAL           Expansion speed [km/s] (default: 0)'
    write(*,'(a)') '  --v_power VAL         Velocity power-law exponent (default: 1)'
    write(*,'(a)') '  --refine_boundary     Refine cells at sphere surface'
    write(*,'(a)') '  --boundary_level_max N  Surface refinement level'
    write(*,'(a)') '  --cone_opening VAL    Bicone half-opening angle [deg] (0=sphere)'
    write(*,'(a)') '  -o FILE               Output file (default: amr_sphere.h5)'
    write(*,'(a)') '  -h, --help            Show this help'
  end subroutine

end program make_amr_sphere_radial
