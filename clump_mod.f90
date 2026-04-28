module clump_mod
!---------------------------------------------------------------------------
! Clumpy medium support for LaRT_v2.00.
!
! Geometry: N_clumps spherical clumps of radius cl_radius placed uniformly
! at random (non-overlapping via RSA) inside a sphere of radius sphere_R.
! All clumps have uniform opacity (cl_rhokap), temperature (cl_Dfreq,
! cl_voigt_a), and independent Gaussian random bulk velocities.
!
! Algorithms:
!   Placement : linked-list RSA, O(N_cl) expected time for f_vol < 35%.
!   Raytrace  : DDA through CSR acceleration grid (Amanatides & Woo 1987).
!   Memory    : MPI-3 shared memory (one copy per node).
!---------------------------------------------------------------------------
  use define
  use memory_mod
  use voigt_mod
  use random
  use mpi
  implicit none
  public

  !--- Physical properties of the clump population
  integer(int64), save :: N_clumps   = 0_int64
  real(kind=wp),  save :: sphere_R   = 0.0_wp   ! outer sphere radius [code units]
  real(kind=wp),  save :: cl_radius  = 0.0_wp   ! clump radius [code units]
  real(kind=wp),  save :: cl_radius2 = 0.0_wp   ! cl_radius**2
  real(kind=wp),  save :: cl_rhokap  = 0.0_wp   ! opacity/code-unit inside clumps
  real(kind=wp),  save :: cl_voigt_a = 0.0_wp   ! Voigt damping parameter
  real(kind=wp),  save :: cl_Dfreq   = 0.0_wp   ! Doppler frequency [Hz]
  real(kind=wp),  save :: cl_vtherm      = 0.0_wp   ! thermal velocity [km/s] = cl_Dfreq * lambda0 * um2km
  real(kind=wp),  save :: cl_temperature = 0.0_wp   ! actual clump temperature [K]

  !--- Clump positions [code units] and bulk velocities – MPI shared memory.
  !
  !    cl_vx/y/z are stored DIMENSIONLESSLY as v / cl_vtherm (matching the
  !    Cartesian/AMR grid%vfx convention). Each setter routine
  !    (generate_clumps and assign_clump_velocities_from_type) divides by
  !    cl_vtherm inline before writing, so no separate post-pass rescaling
  !    is needed. write_clumps_fits() multiplies by cl_vtherm to keep the
  !    user-facing FITS output in km/s.
  real(kind=dp), pointer, save :: cl_x(:)  => null()
  real(kind=dp), pointer, save :: cl_y(:)  => null()
  real(kind=dp), pointer, save :: cl_z(:)  => null()
  real(kind=dp), pointer, save :: cl_vx(:) => null()
  real(kind=dp), pointer, save :: cl_vy(:) => null()
  real(kind=dp), pointer, save :: cl_vz(:) => null()

  !--- CSR acceleration grid (MPI shared memory)
  !    Cell (i,j,k), 0-based → 1-based index = 1 + i + j*cgx + k*cgx*cgy
  !    cg_start(icell) .. cg_start(icell+1)-1 → entries in cg_list for that cell
  integer,       save :: cgx = 0, cgy = 0, cgz = 0
  real(kind=wp), save :: cg_xmin, cg_ymin, cg_zmin
  real(kind=wp), save :: cg_dx, cg_dy, cg_dz
  real(kind=wp), save :: cg_inv_dx, cg_inv_dy, cg_inv_dz
  integer(int32), pointer, save :: cg_start(:) => null()   ! size ncells+1
  integer(int32), pointer, save :: cg_list(:)  => null()   ! size total_registrations

contains

  !===========================================================================
  pure integer function cg_cell_idx(i, j, k)
  integer, intent(in) :: i, j, k
  cg_cell_idx = 1 + i + j*cgx + k*cgx*cgy
  end function cg_cell_idx
  !===========================================================================

  !===========================================================================
  subroutine init_clumps(R_sphere)
  !---------------------------------------------------------------------------
  ! Derive N_clumps, compute physical parameters, allocate MPI shared memory,
  ! run RSA placement, and build the CSR acceleration grid.
  !---------------------------------------------------------------------------
  implicit none
  real(kind=wp), intent(in) :: R_sphere
  real(kind=wp) :: temp_cl, vtherm
  integer       :: ierr

  sphere_R   = R_sphere
  cl_radius  = par%clump_radius
  if (cl_radius <= 0.0_wp) then
     if (mpar%p_rank == 0) write(*,*) 'ERROR: clump_radius must be > 0'
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end if
  cl_radius2 = cl_radius * cl_radius

  !--- derive N_clumps
  if (par%clump_N_clumps > 0.0_wp) then
     N_clumps = int(par%clump_N_clumps, int64)
  else if (par%clump_f_vol > 0.0_wp) then
     N_clumps = nint(par%clump_f_vol * (R_sphere/cl_radius)**3, int64)
  else if (par%clump_f_cov > 0.0_wp) then
     N_clumps = nint((4.0_wp/3.0_wp)*par%clump_f_cov*(R_sphere/cl_radius)**2, int64)
  else
     if (mpar%p_rank == 0) &
        write(*,*) 'ERROR: specify clump_N_clumps, clump_f_vol, or clump_f_cov'
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end if
  if (N_clumps <= 0_int64) N_clumps = 1_int64

  !--- clump temperature → Doppler frequency and Voigt parameter
  temp_cl = par%clump_temperature
  if (temp_cl < 0.0_wp) temp_cl = par%temperature
  vtherm     = line%vtherm1 * sqrt(temp_cl)
  cl_Dfreq   = vtherm / (line%wavelength0 * um2km)
  cl_vtherm      = vtherm                          ! km/s; = cl_Dfreq * wavelength0 * um2km
  cl_temperature = temp_cl                         ! actual temperature used [K]
  cl_voigt_a = (line%damping / fourpi) / cl_Dfreq

  !--- clump opacity from tau0, NHI (column density), or nH
  if (par%clump_tau0 > 0.0_wp) then
     cl_rhokap = par%clump_tau0 / (voigt(0.0_wp, cl_voigt_a) * cl_radius)
  else if (par%clump_NHI > 0.0_wp) then
     ! clump_NHI = per-clump column density [cm^-2] from individual clump center to its surface
     ! (distinct from total system column density N_HImax along sightlines through the clumpy sphere)
     ! line%cross0 [cm^2*Hz]; divide by cl_Dfreq [Hz] to get actual cross-section [cm^2]
     ! => cl_rhokap = clump_NHI * (cross0/cl_Dfreq) / cl_radius  [code_unit^-1]
     cl_rhokap = par%clump_NHI * line%cross0 / (cl_Dfreq * cl_radius)
  else if (par%clump_nH > 0.0_wp) then
     if (par%distance2cm <= 0.0_wp) then
        if (mpar%p_rank == 0) write(*,*) 'ERROR: clump_nH requires distance_unit'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
     ! line%cross0 [cm^2*Hz]; divide by cl_Dfreq [Hz] to get actual cross-section [cm^2]
     cl_rhokap = par%clump_nH * line%cross0 * par%distance2cm / cl_Dfreq
  else
     if (mpar%p_rank == 0) write(*,*) 'ERROR: specify clump_tau0, clump_NHI, or clump_nH'
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end if

  if (mpar%p_rank == 0) then
     write(*,'(a,i14)')   ' Clumps: N_clumps  = ', N_clumps
     write(*,'(a,f12.6)') ' Clumps: f_vol     = ', &
           real(N_clumps,wp) * (cl_radius/R_sphere)**3
     write(*,'(a,f12.5)') ' Clumps: f_cov     = ', &
           0.75_wp * real(N_clumps,wp) * (cl_radius/R_sphere)**2
     write(*,'(a,es12.4)') ' Clumps: cl_rhokap = ', cl_rhokap
     write(*,'(a,f12.5)')  ' Clumps: voigt_a   = ', cl_voigt_a
     write(*,'(a,es12.4)') ' Clumps: cl_Dfreq  = ', cl_Dfreq
  end if

  !--- shared memory for clump data
  call create_shared_mem(cl_x,  [int(N_clumps)])
  call create_shared_mem(cl_y,  [int(N_clumps)])
  call create_shared_mem(cl_z,  [int(N_clumps)])
  call create_shared_mem(cl_vx, [int(N_clumps)])
  call create_shared_mem(cl_vy, [int(N_clumps)])
  call create_shared_mem(cl_vz, [int(N_clumps)])

  !--- p_rank=0 generates positions; broadcast to all h_rank=0 (one per node)
  !--- so every node has an identical clump layout in its shared memory.
  if (mpar%h_rank == 0) then
     if (mpar%p_rank == 0) then
        call generate_clumps()
        if (len_trim(par%velocity_type) > 0) call assign_clump_velocities_from_type()
     end if
     call MPI_BCAST(cl_x,  int(N_clumps), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_y,  int(N_clumps), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_z,  int(N_clumps), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_vx, int(N_clumps), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_vy, int(N_clumps), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     call MPI_BCAST(cl_vz, int(N_clumps), MPI_DOUBLE_PRECISION, 0, mpar%SAME_HRANK_COMM, ierr)
     ! cl_vx/y/z are already stored as v / cl_vtherm (normalised inline by
     ! generate_clumps and assign_clump_velocities_from_type).
  end if
  call MPI_BARRIER(mpar%hostcomm, ierr)

  !--- h_rank=0 builds CSR grid; barrier before use
  call build_clump_csr()
  call MPI_BARRIER(mpar%hostcomm, ierr)

  end subroutine init_clumps
  !===========================================================================

  !===========================================================================
  subroutine generate_clumps()
  !---------------------------------------------------------------------------
  ! RSA with linked-list grid acceleration. Called only on h_rank=0.
  !---------------------------------------------------------------------------
  implicit none
  integer        :: rg, ncells_rsa
  real(kind=wp)  :: rg_cell, min_sep2
  real(kind=wp)  :: xc, yc, zc, dx, dy, dz, d2
  integer        :: ig, jg, kg, ig2, jg2, kg2, icell_rsa, jnb
  integer(int64) :: icl, n_attempts
  logical        :: overlap
  integer, allocatable :: head(:), nxt(:)

  rg         = min(512, max(32, int(real(N_clumps,wp)**(1.0_wp/3.0_wp)) + 1))
  ncells_rsa = rg**3
  rg_cell    = 2.0_wp * sphere_R / real(rg, wp)
  min_sep2   = (2.0_wp * cl_radius)**2

  allocate(head(ncells_rsa),    source=-1)
  allocate(nxt(int(N_clumps)), source=-1)

  if (mpar%p_rank == 0) &
     write(*,'(a,i5,a,i14,a)') ' RSA: grid ', rg, '^3, placing ', N_clumps, ' clumps...'

  n_attempts = 0_int64
  icl        = 0_int64

  do while (icl < N_clumps)
     n_attempts = n_attempts + 1_int64

     !--- uniform random point inside sphere (box rejection)
     do
        xc = (2.0_wp * rand_number() - 1.0_wp) * sphere_R
        yc = (2.0_wp * rand_number() - 1.0_wp) * sphere_R
        zc = (2.0_wp * rand_number() - 1.0_wp) * sphere_R
        if (xc*xc + yc*yc + zc*zc <= sphere_R*sphere_R) exit
     end do

     !--- cell in RSA grid (0-based)
     ig = min(rg-1, max(0, int((xc + sphere_R) / rg_cell)))
     jg = min(rg-1, max(0, int((yc + sphere_R) / rg_cell)))
     kg = min(rg-1, max(0, int((zc + sphere_R) / rg_cell)))

     !--- check 27 neighbours for overlap
     overlap = .false.
     outer: do kg2 = max(0,kg-1), min(rg-1,kg+1)
        do jg2 = max(0,jg-1), min(rg-1,jg+1)
           do ig2 = max(0,ig-1), min(rg-1,ig+1)
              icell_rsa = 1 + ig2 + jg2*rg + kg2*rg*rg
              jnb = head(icell_rsa)
              do while (jnb > 0)
                 dx = xc - cl_x(jnb);  dy = yc - cl_y(jnb);  dz = zc - cl_z(jnb)
                 d2 = dx*dx + dy*dy + dz*dz
                 if (d2 < min_sep2) then
                    overlap = .true.;  exit outer
                 end if
                 jnb = nxt(jnb)
              end do
           end do
        end do
     end do outer
     if (overlap) cycle

     !--- accept
     icl = icl + 1_int64
     cl_x(icl) = real(xc, dp);  cl_y(icl) = real(yc, dp);  cl_z(icl) = real(zc, dp)
     if (par%clump_sigma_v > 0.0_wp) then
        ! store v / cl_vtherm directly (dimensionless, matches AMR convention)
        cl_vx(icl) = real(par%clump_sigma_v / cl_vtherm * rand_gauss(), dp)
        cl_vy(icl) = real(par%clump_sigma_v / cl_vtherm * rand_gauss(), dp)
        cl_vz(icl) = real(par%clump_sigma_v / cl_vtherm * rand_gauss(), dp)
     end if

     !--- insert into RSA linked list
     icell_rsa       = 1 + ig + jg*rg + kg*rg*rg
     nxt(int(icl))   = head(icell_rsa)
     head(icell_rsa) = int(icl)

     if (mod(icl, 1000000_int64) == 0_int64 .and. mpar%p_rank == 0) &
        write(*,'(a,i14,a,i14)') '   placed ', icl, ' / ', N_clumps
  end do

  if (mpar%p_rank == 0) write(*,'(a,f6.1,a)') &
     ' RSA done, acceptance rate = ', &
     real(N_clumps,wp)/real(n_attempts,wp)*100.0_wp, '%'

  deallocate(head, nxt)
  end subroutine generate_clumps
  !===========================================================================

  !===========================================================================
  subroutine assign_clump_velocities_from_type()
  !---------------------------------------------------------------------------
  ! Add a systematic velocity component to each clump based on par%velocity_type
  ! and the clump's centre position.  Called only on h_rank=0 after
  ! generate_clumps(); adds to any existing sigma_v random component.
  !
  ! Supported types (same parameter names as grid_mod_car):
  !   'hubble'              : v = Vexp * (r / sphere_R)       [linear expansion]
  !   'constant_radial'     : v = Vexp * (r / |r|)           [uniform outflow]
  !   'parallel_velocity'   : v = (Vx, Vy, Vz)               [uniform bulk]
  !   'ssh'                 : Song, Seon & Hwang (2020) galaxy model
  !   'rotating_solid_body' : v = Vrot * (-y, x, 0) / sphere_R
  !   'rotating_galaxy_halo': flat rotation curve (Vrot, rinner)
  !---------------------------------------------------------------------------
  implicit none
  integer(int64) :: icl
  real(kind=wp)  :: xc, yc, zc, rr, rr_cyl, Vscale, vx, vy, vz

  do icl = 1_int64, N_clumps
     xc = real(cl_x(icl), wp)
     yc = real(cl_y(icl), wp)
     zc = real(cl_z(icl), wp)
     rr = sqrt(xc**2 + yc**2 + zc**2)
     vx = 0.0_wp;  vy = 0.0_wp;  vz = 0.0_wp

     select case (trim(par%velocity_type))

     case ('hubble')
        ! v_i = Vexp * r_i / sphere_R  [km/s]
        vx = par%Vexp * xc / sphere_R
        vy = par%Vexp * yc / sphere_R
        vz = par%Vexp * zc / sphere_R

     case ('constant_radial')
        ! v = Vexp * r_hat  [km/s]
        if (rr > 0.0_wp) then
           vx = par%Vexp * xc / rr
           vy = par%Vexp * yc / rr
           vz = par%Vexp * zc / rr
        end if

     case ('parallel_velocity')
        ! uniform bulk velocity  [km/s]
        vx = par%Vx
        vy = par%Vy
        vz = par%Vz

     case ('ssh')
        ! Song, Seon & Hwang (2020): linear inside rpeak, then Vpeak + DeltaV*(r-rpeak)/(R-rpeak)
        if (rr > 0.0_wp) then
           if (rr < par%rpeak) then
              Vscale = par%Vpeak / par%rpeak
              vx = Vscale * xc
              vy = Vscale * yc
              vz = Vscale * zc
           else
              Vscale = par%Vpeak + par%DeltaV * (rr - par%rpeak) / (sphere_R - par%rpeak)
              vx = Vscale * xc / rr
              vy = Vscale * yc / rr
              vz = Vscale * zc / rr
           end if
        end if

     case ('rotating_solid_body')
        ! solid-body rotation about z-axis: v = Vrot * (-y, x, 0) / sphere_R
        vx = -par%Vrot * yc / sphere_R
        vy =  par%Vrot * xc / sphere_R

     case ('rotating_galaxy_halo')
        ! flat rotation curve: Vrot inside rinner, then Vrot * rinner / rr
        rr_cyl = sqrt(xc**2 + yc**2)
        if (rr_cyl > 0.0_wp) then
           if (rr_cyl < par%rinner) then
              vx = -par%Vrot * yc / par%rinner
              vy =  par%Vrot * xc / par%rinner
           else
              vx = -par%Vrot * yc / rr_cyl
              vy =  par%Vrot * xc / rr_cyl
           end if
        end if

     end select

     ! store v / cl_vtherm directly (dimensionless, matches AMR convention)
     cl_vx(icl) = cl_vx(icl) + real(vx / cl_vtherm, dp)
     cl_vy(icl) = cl_vy(icl) + real(vy / cl_vtherm, dp)
     cl_vz(icl) = cl_vz(icl) + real(vz / cl_vtherm, dp)
  end do

  if (mpar%p_rank == 0) &
     write(*,'(2a)') ' Clumps: velocity_type   = ', trim(par%velocity_type)

  end subroutine assign_clump_velocities_from_type
  !===========================================================================

  !===========================================================================
  subroutine build_clump_csr()
  !---------------------------------------------------------------------------
  ! Two-pass CSR construction. All ranks allocate shared memory;
  ! h_rank=0 fills it. Caller handles barrier after return.
  !---------------------------------------------------------------------------
  implicit none
  integer        :: ncells, total_regs
  integer        :: i, j, k, imin, imax, jmin, jmax, kmin, kmax, icell
  integer(int64) :: icl
  integer,       allocatable :: cnt(:)
  integer        :: ierr

  !--- CSR grid size: ~1 clump per cell
  cgx = min(512, max(32, int(real(N_clumps,wp)**(1.0_wp/3.0_wp)) + 1))
  cgy = cgx;  cgz = cgx
  ncells = cgx * cgy * cgz

  !--- bounding box: sphere + one clump-radius margin
  cg_xmin = -(sphere_R + cl_radius);  cg_ymin = cg_xmin;  cg_zmin = cg_xmin
  cg_dx   = (2.0_wp*(sphere_R + cl_radius)) / real(cgx, wp)
  cg_dy   = cg_dx;  cg_dz = cg_dx
  cg_inv_dx = 1.0_wp / cg_dx
  cg_inv_dy = 1.0_wp / cg_dy
  cg_inv_dz = 1.0_wp / cg_dz

  !--- allocate cg_start on all ranks (h_rank=0 fills)
  call create_shared_mem(cg_start, [ncells + 1])

  !--- pass 1: count (h_rank=0 only)
  if (mpar%h_rank == 0) then
     allocate(cnt(ncells), source=0)
     do icl = 1_int64, N_clumps
        call clump_cell_range(icl, imin, imax, jmin, jmax, kmin, kmax)
        do k = kmin, kmax
           do j = jmin, jmax
              do i = imin, imax
                 icell = cg_cell_idx(i, j, k)
                 cnt(icell) = cnt(icell) + 1
              end do
           end do
        end do
     end do
     !--- prefix sum → cg_start (1-based: cg_start(1)=1, cg_start(ncells+1)=total+1)
     cg_start(1) = 1
     do icell = 1, ncells
        cg_start(icell+1) = cg_start(icell) + cnt(icell)
     end do
  end if
  call MPI_BARRIER(mpar%hostcomm, ierr)

  total_regs = cg_start(ncells+1) - 1
  call create_shared_mem(cg_list, [total_regs])

  !--- pass 2: fill (h_rank=0 only)
  if (mpar%h_rank == 0) then
     cnt(:) = 0
     do icl = 1_int64, N_clumps
        call clump_cell_range(icl, imin, imax, jmin, jmax, kmin, kmax)
        do k = kmin, kmax
           do j = jmin, jmax
              do i = imin, imax
                 icell = cg_cell_idx(i, j, k)
                 cg_list(cg_start(icell) + cnt(icell)) = int(icl, int32)
                 cnt(icell) = cnt(icell) + 1
              end do
           end do
        end do
     end do
     deallocate(cnt)
     if (mpar%p_rank == 0) write(*,'(a,i12,a,i5,a)') &
           ' CSR grid: ', total_regs, ' registrations in ', cgx, '^3 cells.'
  end if

  end subroutine build_clump_csr
  !===========================================================================

  !===========================================================================
  pure subroutine clump_cell_range(icl, imin, imax, jmin, jmax, kmin, kmax)
  !---------------------------------------------------------------------------
  ! CSR cells whose box overlaps clump icl's sphere.
  !---------------------------------------------------------------------------
  integer(int64), intent(in)  :: icl
  integer,        intent(out) :: imin, imax, jmin, jmax, kmin, kmax
  imin = max(0, int((cl_x(icl) - cg_xmin - cl_radius) * cg_inv_dx))
  imax = min(cgx-1, int((cl_x(icl) - cg_xmin + cl_radius) * cg_inv_dx))
  jmin = max(0, int((cl_y(icl) - cg_ymin - cl_radius) * cg_inv_dy))
  jmax = min(cgy-1, int((cl_y(icl) - cg_ymin + cl_radius) * cg_inv_dy))
  kmin = max(0, int((cl_z(icl) - cg_zmin - cl_radius) * cg_inv_dz))
  kmax = min(cgz-1, int((cl_z(icl) - cg_zmin + cl_radius) * cg_inv_dz))
  end subroutine clump_cell_range
  !===========================================================================

  !===========================================================================
  pure subroutine ray_sphere_isect(ox, oy, oz, kx, ky, kz, cx, cy, cz, &
                                    t_entry, t_exit, hit)
  !---------------------------------------------------------------------------
  ! Ray–sphere intersection.  Ray: P(t) = origin + t*dir, t arbitrary.
  ! Sphere: centre (cx,cy,cz), squared-radius cl_radius2.
  ! Returns hit=.true. and t_entry <= t_exit when t_exit > 0.
  !---------------------------------------------------------------------------
  real(kind=wp), intent(in)  :: ox, oy, oz, kx, ky, kz, cx, cy, cz
  real(kind=wp), intent(out) :: t_entry, t_exit
  logical,       intent(out) :: hit
  real(kind=wp) :: rx, ry, rz, b, disc
  rx = ox - cx;  ry = oy - cy;  rz = oz - cz
  b    = rx*kx + ry*ky + rz*kz
  disc = b*b - (rx*rx + ry*ry + rz*rz) + cl_radius2
  if (disc < 0.0_wp) then
     hit = .false.;  t_entry = 0.0_wp;  t_exit = 0.0_wp
  else
     disc    = sqrt(disc)
     t_entry = -b - disc
     t_exit  = -b + disc
     hit     = (t_exit > 0.0_wp)
  end if
  end subroutine ray_sphere_isect
  !===========================================================================

  !===========================================================================
  subroutine find_next_clump(xp, yp, zp, kx, ky, kz, skip_icl, t_max, &
                              t_entry, t_exit, icl_found, found)
  !---------------------------------------------------------------------------
  ! DDA through CSR grid to find the nearest clump hit along the ray
  !   P(t) = (xp,yp,zp) + t*(kx,ky,kz),  0 < t <= t_max.
  ! Clump skip_icl (> 0) is excluded (the one just exited).
  !
  ! Pattern: Amanatides & Woo (1987).
  !   tx,ty,tz = absolute path length to NEXT x/y/z face crossing.
  !   d        = path length at current cell entry.
  ! Stopping: once d > best_te, no earlier hit can be found.
  !---------------------------------------------------------------------------
  real(kind=wp),  intent(in)  :: xp, yp, zp, kx, ky, kz
  integer(int64), intent(in)  :: skip_icl
  real(kind=wp),  intent(in)  :: t_max
  real(kind=wp),  intent(out) :: t_entry, t_exit
  integer(int64), intent(out) :: icl_found
  logical,        intent(out) :: found

  integer  :: ci, cj, ck, si, sj, sk, icell, ip
  integer(int64) :: icl
  real(kind=wp) :: tx, ty, tz, delx, dely, delz, d
  real(kind=wp) :: te, tx2, best_te, best_tx2
  integer(int64):: best_icl
  logical  :: hit

  found    = .false.
  best_te  = hugest
  best_tx2 = 0.0_wp
  best_icl = 0_int64
  d        = 0.0_wp

  !--- starting cell (clamped to grid)
  ci = max(0, min(cgx-1, int((xp - cg_xmin) * cg_inv_dx)))
  cj = max(0, min(cgy-1, int((yp - cg_ymin) * cg_inv_dy)))
  ck = max(0, min(cgz-1, int((zp - cg_zmin) * cg_inv_dz)))

  !--- DDA setup: absolute path lengths to first face crossing
  if (kx > 0.0_wp) then
     si   =  1
     tx   = ((cg_xmin + real(ci+1,wp)*cg_dx) - xp) / kx
     delx =  cg_dx / kx
  else if (kx < 0.0_wp) then
     si   = -1
     tx   = ((cg_xmin + real(ci,wp)*cg_dx) - xp) / kx
     delx = -cg_dx / kx
  else
     si = 0;  tx = hugest;  delx = hugest
  end if

  if (ky > 0.0_wp) then
     sj   =  1
     ty   = ((cg_ymin + real(cj+1,wp)*cg_dy) - yp) / ky
     dely =  cg_dy / ky
  else if (ky < 0.0_wp) then
     sj   = -1
     ty   = ((cg_ymin + real(cj,wp)*cg_dy) - yp) / ky
     dely = -cg_dy / ky
  else
     sj = 0;  ty = hugest;  dely = hugest
  end if

  if (kz > 0.0_wp) then
     sk   =  1
     tz   = ((cg_zmin + real(ck+1,wp)*cg_dz) - zp) / kz
     delz =  cg_dz / kz
  else if (kz < 0.0_wp) then
     sk   = -1
     tz   = ((cg_zmin + real(ck,wp)*cg_dz) - zp) / kz
     delz = -cg_dz / kz
  else
     sk = 0;  tz = hugest;  delz = hugest
  end if

  do while(.true.)
     !--- stopping: current cell starts beyond the best hit or past t_max
     if (d > best_te .or. d > t_max) exit

     !--- check all clumps in current cell
     icell = cg_cell_idx(ci, cj, ck)
     do ip = cg_start(icell), cg_start(icell+1) - 1
        icl = int(cg_list(ip), int64)
        if (icl == skip_icl) cycle
        call ray_sphere_isect(xp, yp, zp, kx, ky, kz, &
             real(cl_x(icl),wp), real(cl_y(icl),wp), real(cl_z(icl),wp), &
             te, tx2, hit)
        if (hit .and. tx2 > 0.0_wp .and. te < best_te .and. &
            (te > 0.0_wp .or. icl /= skip_icl)) then
           best_te  = te
           best_tx2 = tx2
           best_icl = icl
        end if
     end do

     !--- advance to next cell (Amanatides & Woo pattern)
     if (tx <= ty .and. tx <= tz) then
        d  = tx
        ci = ci + si;  if (ci < 0 .or. ci >= cgx) exit
        tx = tx + delx
     else if (ty <= tz) then
        d  = ty
        cj = cj + sj;  if (cj < 0 .or. cj >= cgy) exit
        ty = ty + dely
     else
        d  = tz
        ck = ck + sk;  if (ck < 0 .or. ck >= cgz) exit
        tz = tz + delz
     end if
  end do

  if (best_icl > 0_int64 .and. best_te <= t_max) then
     found     = .true.
     t_entry   = best_te
     t_exit    = min(best_tx2, t_max)
     icl_found = best_icl
  end if

  end subroutine find_next_clump
  !===========================================================================

  !===========================================================================
  pure real(kind=wp) function clump_exit_dist(xp, yp, zp, kx, ky, kz, icl)
  !---------------------------------------------------------------------------
  ! Distance from (xp,yp,zp) moving in direction (kx,ky,kz) to the exit point
  ! of clump icl. Assumes the photon is currently inside the clump.
  !---------------------------------------------------------------------------
  real(kind=wp),  intent(in) :: xp, yp, zp, kx, ky, kz
  integer(int64), intent(in) :: icl
  real(kind=wp) :: rx, ry, rz, b, disc
  rx = xp - real(cl_x(icl),wp)
  ry = yp - real(cl_y(icl),wp)
  rz = zp - real(cl_z(icl),wp)
  b    = rx*kx + ry*ky + rz*kz
  disc = b*b - (rx*rx + ry*ry + rz*rz) + cl_radius2
  if (disc < 0.0_wp) disc = 0.0_wp
  clump_exit_dist = max(0.0_wp, -b + sqrt(disc))
  end function clump_exit_dist
  !===========================================================================

  !===========================================================================
  pure real(kind=wp) function sphere_exit_dist(xp, yp, zp, kx, ky, kz)
  !---------------------------------------------------------------------------
  ! Distance to the exit of the outer sphere (radius sphere_R) from (xp,yp,zp).
  !---------------------------------------------------------------------------
  real(kind=wp), intent(in) :: xp, yp, zp, kx, ky, kz
  real(kind=wp) :: b, disc
  b    = xp*kx + yp*ky + zp*kz
  disc = b*b - (xp*xp + yp*yp + zp*zp) + sphere_R*sphere_R
  if (disc < 0.0_wp) disc = 0.0_wp
  sphere_exit_dist = max(0.0_wp, -b + sqrt(disc))
  end function sphere_exit_dist
  !===========================================================================

  !===========================================================================
  subroutine destroy_clumps()
  implicit none
  ! destroy_mem fails on non-root ranks because get_window compares against
  ! rank-0 base addresses (MPI_WIN_SHARED_QUERY rank=0), which differ from
  ! local pointers on other ranks.  Use destroy_shared_mem_all() so that
  ! MPI_WIN_FREE is called collectively, then nullify all pointers manually.
  call destroy_shared_mem_all()
  nullify(cl_x, cl_y, cl_z, cl_vx, cl_vy, cl_vz)
  nullify(cg_start, cg_list)
  N_clumps = 0_int64
  end subroutine destroy_clumps
  !===========================================================================

  !===========================================================================
  subroutine write_clumps_fits(fname)
  !---------------------------------------------------------------------------
  ! Write clump positions, velocities, and physical parameters to a FITS
  ! binary table.  Called only on p_rank == 0 after init_clumps().
  !
  ! HDU 1 (primary): empty image, all clump scalars as header keywords.
  ! HDU 2 (BinTable): columns X, Y, Z [code units], VX, VY, VZ [km/s].
  !---------------------------------------------------------------------------
  use define
  use fitsio_mod
  implicit none
  character(len=*), intent(in) :: fname

  integer :: unit, status, bitpix
  real(kind=real64), allocatable :: tmp(:)
  integer(int64) :: ncl
  real(kind=wp)  :: f_vol_actual, f_cov_actual

  status = 0
  ncl    = N_clumps
  bitpix = -32   ! save as float32 (sufficient precision for positions/velocities)

  !--- Compute realized f_vol and f_cov from actual placed clumps
  f_vol_actual  = real(ncl,wp) * (cl_radius / sphere_R)**3
  f_cov_actual  = 0.75_wp * real(ncl,wp) * (cl_radius / sphere_R)**2

  call fits_open_new(unit, trim(fname), status)
  if (status /= 0) then
     write(*,*) 'WARNING: write_clumps_fits: cannot open ', trim(fname)
     return
  end if

  !--- BinTable HDU: X, Y, Z, VX, VY, VZ (first column creates the BinTable HDU)
  allocate(tmp(ncl))

  tmp = cl_x(1:ncl)
  call fits_write_table_column(unit, 'X',  tmp, status, bitpix)
  tmp = cl_y(1:ncl)
  call fits_write_table_column(unit, 'Y',  tmp, status, bitpix)
  tmp = cl_z(1:ncl)
  call fits_write_table_column(unit, 'Z',  tmp, status, bitpix)
  !--- cl_v* are stored as v/vtherm internally (since 2026-04-27);
  !    multiply back by cl_vtherm to write velocities in km/s in the FITS file.
  tmp = cl_vx(1:ncl) * cl_vtherm
  call fits_write_table_column(unit, 'VX', tmp, status, bitpix)
  tmp = cl_vy(1:ncl) * cl_vtherm
  call fits_write_table_column(unit, 'VY', tmp, status, bitpix)
  tmp = cl_vz(1:ncl) * cl_vtherm
  call fits_write_table_column(unit, 'VZ', tmp, status, bitpix)

  deallocate(tmp)

  !--- Write keywords into the BinTable HDU (current HDU after fits_write_table_column)
  call fits_put_keyword(unit, 'N_CLUMPS',  ncl,           'number of clumps (realized)',          status)
  call fits_put_keyword(unit, 'SPHERE_R',  sphere_R,      'outer sphere radius [code units]',      status)
  call fits_put_keyword(unit, 'CL_RAD',    cl_radius,     'clump radius [code units]',             status)
  call fits_put_keyword(unit, 'F_VOL',     f_vol_actual,  'volume filling factor (realized)',       status)
  call fits_put_keyword(unit, 'F_COV',     f_cov_actual,  'covering factor (realized)',             status)
  call fits_put_keyword(unit, 'TAU0',      par%clump_tau0,    'line-center tau (center to surface)',status)
  call fits_put_keyword(unit, 'SIGMA_V',   par%clump_sigma_v, 'bulk velocity sigma [km/s]',        status)
  call fits_put_keyword(unit, 'TEMP_CL',   cl_temperature,'clump temperature [K] (actual)',         status)
  call fits_put_keyword(unit, 'RHOKAP',    cl_rhokap,     'opacity/code-unit inside clumps',       status)
  call fits_put_keyword(unit, 'CL_DFREQ',  cl_Dfreq,      'Doppler frequency [Hz]',                status)
  call fits_put_keyword(unit, 'VTHERM',    cl_vtherm,     'thermal velocity [km/s]',               status)
  call fits_put_keyword(unit, 'VOIGT_A',   cl_voigt_a,    'Voigt damping parameter',               status)
  call fits_put_keyword(unit, 'RMAX',      par%rmax,      'outer sphere radius input (par%rmax)',   status)
  call fits_put_keyword(unit, 'IN_FCOV',   par%clump_f_cov,   'covering factor (input)',            status)
  call fits_put_keyword(unit, 'IN_FVOL',   par%clump_f_vol,   'volume filling factor (input)',      status)
  call fits_put_keyword(unit, 'IN_NCL',    par%clump_N_clumps,'N_clumps (input)',                   status)
  call fits_put_keyword(unit, 'IN_NHI',    par%clump_NHI,     'per-clump NHI input [cm^-2] (clump center to surface)',status)
  call fits_put_keyword(unit, 'IN_NH',     par%clump_nH,      'clump nH density input [cm^-3]',    status)
  call fits_put_keyword(unit, 'IN_TEMP',   par%clump_temperature,'clump temperature input [K]',     status)

  call fits_close(unit, status)

  if (mpar%p_rank == 0) write(*,'(2a)') ' Clumps saved to ', trim(fname)
  end subroutine write_clumps_fits
  !===========================================================================

end module clump_mod
