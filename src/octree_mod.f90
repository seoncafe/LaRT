module octree_mod
  use define
  use memory_mod
  use voigt_mod
  implicit none
  public

  !-----------------------------------------------------------------------
  ! AMR grid type: octree data structure for Adaptive Mesh Refinement.
  !
  ! Tree storage is a flat array of cells (internal + leaf), 1-indexed.
  ! Root = cell 1.  Octant ordering: child index = 1 + ix + 2*iy + 4*iz
  ! where ix=1 if x >= cell_center_x, else 0 (similarly for iy, iz).
  !
  ! Tree structure arrays and physical data arrays are allocated as
  ! MPI-3 shared memory (one copy per node), via create_shared_mem.
  ! Output arrays (Jout, Jin, Jabs) are per-rank (regular allocate).
  !-----------------------------------------------------------------------
  type amr_grid_type

    ! ------ domain geometry ------
    real(wp) :: xmin = 0.0_wp, xmax = 1.0_wp, xrange = 1.0_wp
    real(wp) :: ymin = 0.0_wp, ymax = 1.0_wp, yrange = 1.0_wp
    real(wp) :: zmin = 0.0_wp, zmax = 1.0_wp, zrange = 1.0_wp
    real(wp) :: L_box = 1.0_wp    ! box length (= xrange; kept for internal octree use)

    ! ------ octree structure ------
    integer  :: ncells   = 0      ! total cells (internal + leaf)
    integer  :: nleaf    = 0      ! number of leaf cells
    integer  :: levelmax = 0      ! maximum AMR level present

    ! Arrays of size ncells (tree topology) -- shared memory (pointer)
    integer,  pointer :: parent(:)      => null()  ! parent index; 0 for root
    integer,  pointer :: children(:,:)  => null()  ! (8, ncells): child cell indices; 0 = no child
    integer,  pointer :: level(:)       => null()  ! cell AMR level (0 = root)
    real(wp), pointer :: cx(:)          => null()  ! cell center x
    real(wp), pointer :: cy(:)          => null()  ! cell center y
    real(wp), pointer :: cz(:)          => null()  ! cell center z
    real(wp), pointer :: ch(:)          => null()  ! cell half-width

    ! Leaf-index mapping -- shared memory (pointer)
    integer,  pointer :: ileaf(:)           => null()  ! ncells: >0 is leaf index; 0 = internal
    integer,  pointer :: icell_of_leaf(:)   => null()  ! nleaf -> cell index

    ! ------ face-neighbor lookup table (size 6 x ncells) -- shared memory ------
    ! neighbor(iface, icell) = cell index of the face neighbor in direction iface:
    !   iface 1=+x  2=-x  3=+y  4=-y  5=+z  6=-z
    ! 0 => outside the box.
    integer, pointer :: neighbor(:,:)  => null()  ! (6, ncells)

    ! ------ physical data at leaf cells (size nleaf) -- shared memory ------
    real(wp), pointer :: rhokap(:)   => null()  ! HI opacity per unit length at line center
    real(wp), pointer :: rhokapD(:)  => null()  ! dust opacity per unit length
    real(wp), pointer :: Dfreq(:)    => null()  ! local Doppler frequency
    real(wp), pointer :: voigt_a(:)  => null()  ! Voigt damping parameter
    real(wp), pointer :: vfx(:)      => null()  ! x-velocity / v_thermal(local)
    real(wp), pointer :: vfy(:)      => null()  ! y-velocity / v_thermal(local)
    real(wp), pointer :: vfz(:)      => null()  ! z-velocity / v_thermal(local)
    real(wp), pointer :: Pem(:)      => null()  ! diffuse emissivity helper array
    real(wp), pointer :: Pwgt(:)     => null()  ! composite-bias emission weights
    integer,  pointer :: alias(:)    => null()  ! alias table for Pem

    ! ------ frequency / spectral grid -- shared memory ------
    integer  :: nxfreq    = 0
    real(wp) :: xfreq_min = 0.0_wp
    real(wp) :: xfreq_max = 0.0_wp
    real(wp) :: dxfreq    = 0.0_wp
    real(wp) :: dwave     = 0.0_wp
    real(wp) :: Dfreq_ref  = 0.0_wp  ! reference Doppler frequency
    real(wp) :: voigt_amean = 0.0_wp
    real(wp) :: Dfreq_mean  = 0.0_wp
    real(wp) :: xcrit  = 0.0_wp      ! core-skip threshold
    real(wp) :: xcrit2 = 0.0_wp

    real(wp), pointer :: xfreq(:)     => null()
    real(wp), pointer :: velocity(:)  => null()
    real(wp), pointer :: wavelength(:)=> null()

    ! ------ spectral output arrays (per frequency bin) -- per-rank (NOT shared) ------
    real(wp), allocatable :: Jout(:)  ! escaped spectrum
    real(wp), allocatable :: Jin(:)   ! input (injected) spectrum
    real(wp), allocatable :: Jabs(:)  ! absorbed by dust spectrum
    real(wp), allocatable :: Jmu(:,:) ! escaped spectrum binned by mu = cos(theta_z)

    ! ------ Ly-beta (line_type = 8): band-2 (H-alpha) spectral grid + spectra ------
    ! Band-2 axis is in H-alpha Doppler units (v_th identical to band 1, so the
    ! velocity range maps 1:1).  Mirrors the Cartesian grid_type Ha members.
    integer  :: nxfreq_Ha    = 0
    real(wp) :: xfreq_min_Ha = 0.0_wp, xfreq_max_Ha = 0.0_wp
    real(wp) :: dxfreq_Ha    = 0.0_wp, dwave_Ha     = 0.0_wp
    real(wp), pointer :: Jout_Ha(:) => null()  ! escaped H-alpha spectrum (band 2)
    real(wp), pointer :: Jabs_Ha(:) => null()  ! dust-absorbed H-alpha spectrum (band 2)

    ! ------ CALCJ / CALCP / CALCPnew : mean-intensity & scattering-rate maps ------
    ! Per-leaf storage (the natural AMR counterpart of the Cartesian (i,j,k) arrays).
    ! geometry_JPa: 3 = per-leaf 3D, 2 = cylinder (r,z), 1 = sphere (radial),
    !              -1 = plane-parallel (z bins).  Set from par%geometry_JPa.
    ! Binning is POSITION-based: each path segment (J, Pnew) or scattering event
    ! (Pa) deposits into the bin containing its position, independent of the
    ! leaf size, so bins narrower than the local leaf are resolved correctly.
    ! Normalization is volume-weighted: vol_* holds the gas (rhokap > 0) volume
    ! overlapping each bin, computed by sub-sampling the leaves at setup.
    ! Octant multiplicity for xyz_symmetry / xy_symmetry is folded into vol_*
    ! via the factor nadd (8/4), mirroring grid_mod_car.
#if defined (CALCJ) || defined (CALCP) || defined (CALCPnew)
    integer  :: geometry_JPa = huge(1)
    integer  :: nr_JPa = 0        ! number of radial bins (geom 1, 2)
    integer  :: nz_JPa = 0        ! number of z bins      (geom -1, 2)
    real(wp) :: rmax_JPa = 0.0_wp ! outer radius for radial binning
    real(wp) :: dr_JPa = 0.0_wp   ! radial bin width
    real(wp) :: dz_JPa = 0.0_wp   ! z bin width
    real(wp), pointer :: vol_leaf(:)  => null()  ! (nleaf) leaf volume * nadd (geom 3)
    real(wp), pointer :: vol_sph(:)   => null()  ! (nr)    gas volume per radial bin
    real(wp), pointer :: vol_cyl(:,:) => null()  ! (nr,nz) gas volume per (r,z) bin
    real(wp), pointer :: vol_plane(:) => null()  ! (nz)    gas volume per z bin
#endif
#ifdef CALCJ
    real(dp), pointer :: J(:,:)    => null()  ! (nxfreq, nleaf)  geom 3
    real(wp), pointer :: J1(:,:)   => null()  ! (nxfreq, nr) geom 1 / (nxfreq, nz) geom -1
    real(wp), pointer :: J2(:,:,:) => null()  ! (nxfreq, nr, nz) geom 2
#endif
#ifdef CALCP
    real(dp), pointer :: Pa(:)   => null()    ! (nleaf) geom 3
    real(wp), pointer :: P1(:)   => null()    ! (nr) geom 1 / (nz) geom -1
    real(wp), pointer :: P2(:,:) => null()    ! (nr,nz) geom 2
    !--- Ly-beta (line_type = 8) conversion-rate maps: same shapes/binning/
    !--- normalization as Pa/P1/P2, accumulated at 3p->2s conversion events only.
    real(dp), pointer :: Pc(:)   => null()    ! (nleaf) geom 3
    real(wp), pointer :: Pc1(:)  => null()    ! (nr) geom 1 / (nz) geom -1
    real(wp), pointer :: Pc2(:,:) => null()   ! (nr,nz) geom 2
#endif
#ifdef CALCPnew
    real(dp), pointer :: Pa_new(:)   => null()
    real(wp), pointer :: P1_new(:)   => null()
    real(wp), pointer :: P2_new(:,:) => null()
#endif

  end type amr_grid_type

  ! Global AMR grid instance (accessible from any module that uses octree_mod)
  type(amr_grid_type) :: amr_grid

contains

  !=========================================================================
  ! Find the leaf index for point (x, y, z).
  ! Returns 0 if the point is outside the domain or not covered by a leaf.
  !=========================================================================
  integer function amr_find_leaf(x, y, z) result(ileaf)
    real(wp), intent(in) :: x, y, z
    integer :: icell, ioct

    ileaf = 0
    if (x < amr_grid%xmin .or. x > amr_grid%xmax .or. &
        y < amr_grid%ymin .or. y > amr_grid%ymax .or. &
        z < amr_grid%zmin .or. z > amr_grid%zmax) return

    icell = 1
    do
      if (amr_grid%ileaf(icell) > 0) then
        ileaf = amr_grid%ileaf(icell)
        return
      end if
      ioct = 1
      if (x >= amr_grid%cx(icell)) ioct = ioct + 1
      if (y >= amr_grid%cy(icell)) ioct = ioct + 2
      if (z >= amr_grid%cz(icell)) ioct = ioct + 4
      icell = amr_grid%children(ioct, icell)
      if (icell == 0) return
    end do
  end function amr_find_leaf

  !=========================================================================
  ! Like amr_find_leaf but returns the cell index rather than the leaf index.
  !=========================================================================
  integer function amr_find_cell(x, y, z) result(icell_out)
    real(wp), intent(in) :: x, y, z
    integer :: icell, ioct

    icell_out = 0
    if (x < amr_grid%xmin .or. x > amr_grid%xmax .or. &
        y < amr_grid%ymin .or. y > amr_grid%ymax .or. &
        z < amr_grid%zmin .or. z > amr_grid%zmax) return

    icell = 1
    do
      if (amr_grid%ileaf(icell) > 0) then
        icell_out = icell
        return
      end if
      ioct = 1
      if (x >= amr_grid%cx(icell)) ioct = ioct + 1
      if (y >= amr_grid%cy(icell)) ioct = ioct + 2
      if (z >= amr_grid%cz(icell)) ioct = ioct + 4
      icell = amr_grid%children(ioct, icell)
      if (icell == 0) then
        icell_out = 0
        return
      end if
    end do
  end function amr_find_cell

  !=========================================================================
  ! Returns the deepest cell (leaf or internal) whose volume contains
  ! (x,y,z).  Unlike amr_find_cell, when the descent reaches an internal
  ! cell that has no child in the relevant octant (a "gap"), this returns
  ! that internal cell rather than 0.  Useful for pole traversal and any
  ! code that needs to step across gaps with the geometry of the gap-cell.
  !=========================================================================
  integer function amr_find_enclosing_cell(x, y, z) result(icell_out)
    real(wp), intent(in) :: x, y, z
    integer :: icell, ioct, child

    icell_out = 0
    if (x < amr_grid%xmin .or. x > amr_grid%xmax .or. &
        y < amr_grid%ymin .or. y > amr_grid%ymax .or. &
        z < amr_grid%zmin .or. z > amr_grid%zmax) return

    icell = 1
    do
      if (amr_grid%ileaf(icell) > 0) then
        icell_out = icell
        return
      end if
      ioct = 1
      if (x >= amr_grid%cx(icell)) ioct = ioct + 1
      if (y >= amr_grid%cy(icell)) ioct = ioct + 2
      if (z >= amr_grid%cz(icell)) ioct = ioct + 4
      child = amr_grid%children(ioct, icell)
      if (child == 0) then
        icell_out = icell   ! deepest internal cell that contains (x,y,z)
        return
      end if
      icell = child
    end do
  end function amr_find_enclosing_cell

  !=========================================================================
  ! RASCAS-style per-cell core-skip threshold (Smith+15 Eq.35).
  !   atau_cell = voigt_a(il) * rhokap(il) * dl_face
  !   xcrit     = (atau_cell)^(1/3) / 5   if atau_cell > 1, else 0
  ! dl_face is the minimum distance from (x,y,z) to any face of leaf cell il,
  ! i.e. the radius of the largest sphere centered at the photon position that
  ! is still contained in the cell.  When par%core_skip_global = .true. the
  ! routine returns the volume-averaged amr_grid%xcrit / xcrit2 instead.
  ! Reads from module-level amr_grid.
  !=========================================================================
  subroutine amr_xcrit_local(il, x, y, z, xcrit_out, xcrit2_out)
    integer,  intent(in)  :: il
    real(wp), intent(in)  :: x, y, z
    real(wp), intent(out) :: xcrit_out, xcrit2_out
    integer  :: icell
    real(wp) :: h, dl_face, atau_cell
    real(wp), parameter :: third = 1.0_wp/3.0_wp

    if (par%core_skip_global) then
      xcrit_out  = amr_grid%xcrit
      xcrit2_out = amr_grid%xcrit2
      return
    end if

    xcrit_out  = 0.0_wp
    xcrit2_out = 0.0_wp
    if (il <= 0) return

    icell = amr_grid%icell_of_leaf(il)
    h     = amr_grid%ch(icell)
    dl_face = min( h - abs(x - amr_grid%cx(icell)), &
                   h - abs(y - amr_grid%cy(icell)), &
                   h - abs(z - amr_grid%cz(icell)) )
    if (dl_face <= 0.0_wp) return

    atau_cell = amr_grid%voigt_a(il) * amr_grid%rhokap(il) * dl_face
    if (atau_cell > 1.0_wp) then
      xcrit_out  = atau_cell**third / 5.0_wp
      xcrit2_out = xcrit_out * xcrit_out
    end if
  end subroutine amr_xcrit_local

#if defined (CALCJ) || defined (CALCP) || defined (CALCPnew)
  !-----------------------------------------------------------------------
  ! Position -> bin index for the CALC* profile shapes (geometry_JPa 1/2/-1).
  ! Returns ibin1 (radial or z bin) and ibin2 (z bin, geom 2 only);
  ! ok = .false. when the position falls outside the binned range.
  !-----------------------------------------------------------------------
  subroutine amr_JPa_bin(amr, x0, y0, z0, ibin1, ibin2, ok)
    type(amr_grid_type), intent(in)  :: amr
    real(wp),            intent(in)  :: x0, y0, z0
    integer,             intent(out) :: ibin1, ibin2
    logical,             intent(out) :: ok
    real(wp) :: rr
    ibin1 = 0;  ibin2 = 0;  ok = .false.
    select case (amr%geometry_JPa)
    case (2)
       rr    = sqrt(x0*x0 + y0*y0)
       ibin1 = floor(rr/amr%dr_JPa) + 1
       ibin2 = floor((z0 - amr%zmin)/amr%dz_JPa) + 1
       if (ibin2 < 1)          ibin2 = 1
       if (ibin2 > amr%nz_JPa) ibin2 = amr%nz_JPa
       ok = (ibin1 >= 1 .and. ibin1 <= amr%nr_JPa)
    case (1)
       rr    = sqrt(x0*x0 + y0*y0 + z0*z0)
       ibin1 = floor(rr/amr%dr_JPa) + 1
       ok = (ibin1 >= 1 .and. ibin1 <= amr%nr_JPa)
    case (-1)
       ibin1 = floor((z0 - amr%zmin)/amr%dz_JPa) + 1
       if (ibin1 < 1)          ibin1 = 1
       if (ibin1 > amr%nz_JPa) ibin1 = amr%nz_JPa
       ok = .true.
    end select
  end subroutine amr_JPa_bin
#endif

  !=========================================================================
  ! Line profile for an AMR leaf, dispatched by line%line_type to mirror the
  ! Cartesian calc_voigt1/2/3/HD EXACTLY (same select case as setup.f90):
  !   type 2      -> calc_voigt2  (DnuHK doublet, 1/3 + 2/3)
  !   type 5, 6   -> calc_voigt3  (general multiplet, sum over line%nup)
  !   type 7      -> calc_voigt_HD (combined H + D Lyman-alpha)
  !   default 1,4 -> calc_voigt1  (single Voigt)
  ! Uses the per-leaf amr_grid arrays; xfreq is carried in this leaf's Doppler
  ! units by the AMR raytrace.
  !=========================================================================
  function amr_line_profile(il, xfreq) result(v)
    integer,  intent(in) :: il
    real(wp), intent(in) :: xfreq
    real(wp) :: v, Dnu, a_ratio, f_ratio, xD
    integer  :: i
    select case (line%line_type)
    case (2)
      Dnu = line%DnuHK_Hz / amr_grid%Dfreq(il)
      v = voigt(xfreq + Dnu, amr_grid%voigt_a(il)) * (1.0_wp/3.0_wp) &
        + voigt(xfreq,       amr_grid%voigt_a(il)) * (2.0_wp/3.0_wp)
    case (5, 6)
      v = voigt(xfreq, amr_grid%voigt_a(il))
      do i = 2, line%nup
        Dnu     = line%delE_Hz(i)   / amr_grid%Dfreq(il)
        a_ratio = line%b(i)%damping / line%b(1)%damping
        f_ratio = line%f12(i)       / line%f12(1)
        v = v + voigt(xfreq + Dnu, amr_grid%voigt_a(il)*a_ratio) * f_ratio
      end do
    case (7)
      Dnu = line%delta_nu_HD_Hz / amr_grid%Dfreq(il)
      xD  = (xfreq - Dnu) * line%ratio_Dfreq_HD
      v = voigt(xfreq, amr_grid%voigt_a(il)) &
        + line%nD_over_nH * line%ratio_Dfreq_HD &
          * voigt(xD, amr_grid%voigt_a(il) * line%ratio_voigta_HD)
    case default
      v = voigt(xfreq, amr_grid%voigt_a(il))
    end select
  end function amr_line_profile

  !=========================================================================
  ! Step through a virtual sub-cell of an internal cell: half the size of
  ! `icell`, centered in the sub-octant containing (x, y, z).  Used by
  ! pole traversal and any code that needs to advance across gaps one
  ! sub-octant at a time so it can find leaves that live in other
  ! sub-octants of the same internal cell.
  !=========================================================================
  subroutine amr_gap_exit(icell, x, y, z, kx, ky, kz, t_exit, iface)
    integer,  intent(in)  :: icell
    real(wp), intent(in)  :: x, y, z, kx, ky, kz
    real(wp), intent(out) :: t_exit
    integer,  intent(out) :: iface

    real(wp) :: cx, cy, cz, h_sub
    real(wp) :: cx_sub, cy_sub, cz_sub
    real(wp) :: t(6)
    integer  :: ix, iy, iz

    cx = amr_grid%cx(icell);  cy = amr_grid%cy(icell);  cz = amr_grid%cz(icell)
    h_sub = amr_grid%ch(icell) * 0.5_wp

    ix = 0;  if (x >= cx) ix = 1
    iy = 0;  if (y >= cy) iy = 1
    iz = 0;  if (z >= cz) iz = 1

    cx_sub = cx + real(2*ix - 1, wp) * h_sub
    cy_sub = cy + real(2*iy - 1, wp) * h_sub
    cz_sub = cz + real(2*iz - 1, wp) * h_sub

    if (kx > 0.0_wp) then
      t(1) = (cx_sub + h_sub - x) / kx;  t(2) = hugest
    else if (kx < 0.0_wp) then
      t(1) = hugest;  t(2) = (cx_sub - h_sub - x) / kx
    else
      t(1) = hugest;  t(2) = hugest
    end if
    if (ky > 0.0_wp) then
      t(3) = (cy_sub + h_sub - y) / ky;  t(4) = hugest
    else if (ky < 0.0_wp) then
      t(3) = hugest;  t(4) = (cy_sub - h_sub - y) / ky
    else
      t(3) = hugest;  t(4) = hugest
    end if
    if (kz > 0.0_wp) then
      t(5) = (cz_sub + h_sub - z) / kz;  t(6) = hugest
    else if (kz < 0.0_wp) then
      t(5) = hugest;  t(6) = (cz_sub - h_sub - z) / kz
    else
      t(5) = hugest;  t(6) = hugest
    end if

    iface  = minloc(t, dim=1)
    t_exit = t(iface)
  end subroutine amr_gap_exit

  !=========================================================================
  ! Compute the distance from (x,y,z) traveling (kx,ky,kz) to the nearest
  ! face of cell icell.  iface: 1=+x, 2=-x, 3=+y, 4=-y, 5=+z, 6=-z.
  !=========================================================================
  subroutine amr_cell_exit(icell, x, y, z, kx, ky, kz, t_exit, iface)
    integer,  intent(in)  :: icell
    real(wp), intent(in)  :: x, y, z, kx, ky, kz
    real(wp), intent(out) :: t_exit
    integer,  intent(out) :: iface

    real(wp) :: h, cx, cy, cz, t(6)

    cx = amr_grid%cx(icell)
    cy = amr_grid%cy(icell)
    cz = amr_grid%cz(icell)
    h  = amr_grid%ch(icell)

    if (kx > 0.0_wp) then
      t(1) = (cx + h - x) / kx;  t(2) = hugest
    else if (kx < 0.0_wp) then
      t(1) = hugest;              t(2) = (cx - h - x) / kx
    else
      t(1) = hugest;              t(2) = hugest
    end if

    if (ky > 0.0_wp) then
      t(3) = (cy + h - y) / ky;  t(4) = hugest
    else if (ky < 0.0_wp) then
      t(3) = hugest;              t(4) = (cy - h - y) / ky
    else
      t(3) = hugest;              t(4) = hugest
    end if

    if (kz > 0.0_wp) then
      t(5) = (cz + h - z) / kz;  t(6) = hugest
    else if (kz < 0.0_wp) then
      t(5) = hugest;              t(6) = (cz - h - z) / kz
    else
      t(5) = hugest;              t(6) = hugest
    end if

    iface  = minloc(t, dim=1)
    t_exit = t(iface)
  end subroutine amr_cell_exit

  !=========================================================================
  ! Build an octree from flat arrays of leaf-cell positions and AMR levels.
  !
  ! Tree is built on every rank using local allocatable temporaries (since
  ! all ranks receive identical broadcast data), then transferred to
  ! MPI-3 shared memory (one physical copy per node).
  !=========================================================================
  subroutine amr_build_tree(xleaf, yleaf, zleaf, leaf_level, nleaf, &
                              box_xmin, box_xmax, box_ymin, box_ymax, box_zmin, box_zmax)
    use mpi
    integer,  intent(in) :: nleaf
    real(wp), intent(in) :: xleaf(nleaf), yleaf(nleaf), zleaf(nleaf)
    integer,  intent(in) :: leaf_level(nleaf)
    real(wp), intent(in) :: box_xmin, box_xmax, box_ymin, box_ymax, box_zmin, box_zmax

    ! Local temporaries for dynamic build (move_alloc requires allocatable)
    integer,  allocatable :: t_parent(:), t_children(:,:), t_level(:), t_ileaf(:)
    real(wp), allocatable :: t_cx(:), t_cy(:), t_cz(:), t_ch(:)
    integer,  allocatable :: t_icell_of_leaf(:)

    integer  :: ncells_max, ncells
    integer  :: il, icell, lev, ioct, nc, ix, iy, iz
    integer  :: ierr

    ncells_max = max(nleaf * 2, 16)
    allocate(t_parent(ncells_max),    source=0)
    allocate(t_children(8,ncells_max),source=0)
    allocate(t_level(ncells_max),     source=0)
    allocate(t_ileaf(ncells_max),     source=0)
    allocate(t_cx(ncells_max),   source=0.0_wp)
    allocate(t_cy(ncells_max),   source=0.0_wp)
    allocate(t_cz(ncells_max),   source=0.0_wp)
    allocate(t_ch(ncells_max),   source=0.0_wp)
    allocate(t_icell_of_leaf(nleaf), source=0)

    amr_grid%xmin  = box_xmin;  amr_grid%xmax  = box_xmax;  amr_grid%xrange = box_xmax - box_xmin
    amr_grid%ymin  = box_ymin;  amr_grid%ymax  = box_ymax;  amr_grid%yrange = box_ymax - box_ymin
    amr_grid%zmin  = box_zmin;  amr_grid%zmax  = box_zmax;  amr_grid%zrange = box_zmax - box_zmin
    amr_grid%L_box = box_xmax - box_xmin

    ! Initialise root cell (index 1)
    ncells          = 1
    t_parent(1)     = 0
    t_children(:,1) = 0
    t_level(1)      = 0
    t_ileaf(1)      = 0
    t_cx(1) = (box_xmin + box_xmax) * 0.5_wp
    t_cy(1) = (box_ymin + box_ymax) * 0.5_wp
    t_cz(1) = (box_zmin + box_zmax) * 0.5_wp
    t_ch(1) = amr_grid%L_box * 0.5_wp

    ! Insert each leaf, creating internal nodes as needed
    do il = 1, nleaf
      icell = 1
      do lev = 0, leaf_level(il) - 1
        ioct = 1
        if (xleaf(il) >= t_cx(icell)) ioct = ioct + 1
        if (yleaf(il) >= t_cy(icell)) ioct = ioct + 2
        if (zleaf(il) >= t_cz(icell)) ioct = ioct + 4
        if (t_children(ioct, icell) == 0) then
          ncells = ncells + 1
          if (ncells > ncells_max) then
            call grow_local_arrays(ncells_max*2, ncells_max, &
                t_parent, t_children, t_level, t_ileaf, t_cx, t_cy, t_cz, t_ch)
            ncells_max = size(t_parent)
          end if
          nc = ncells
          t_parent(nc)      = icell
          t_children(:, nc) = 0
          t_level(nc)       = lev + 1
          t_ileaf(nc)       = 0
          t_ch(nc)          = t_ch(icell) * 0.5_wp
          ix = mod(ioct - 1, 2)
          iy = mod((ioct - 1) / 2, 2)
          iz = (ioct - 1) / 4
          t_cx(nc) = t_cx(icell) + real(2*ix - 1, wp) * t_ch(nc)
          t_cy(nc) = t_cy(icell) + real(2*iy - 1, wp) * t_ch(nc)
          t_cz(nc) = t_cz(icell) + real(2*iz - 1, wp) * t_ch(nc)
          t_children(ioct, icell) = nc
        end if
        icell = t_children(ioct, icell)
      end do
      t_ileaf(icell)      = il
      t_icell_of_leaf(il) = icell
    end do

    amr_grid%ncells   = ncells
    amr_grid%nleaf    = nleaf
    amr_grid%levelmax = maxval(leaf_level)

    ! Allocate shared memory and copy from local temps (h_rank==0 writes)
    call create_shared_mem(amr_grid%parent,        [ncells])
    call create_shared_mem(amr_grid%children,      [8, ncells])
    call create_shared_mem(amr_grid%level,         [ncells])
    call create_shared_mem(amr_grid%ileaf,         [ncells])
    call create_shared_mem(amr_grid%cx,            [ncells])
    call create_shared_mem(amr_grid%cy,            [ncells])
    call create_shared_mem(amr_grid%cz,            [ncells])
    call create_shared_mem(amr_grid%ch,            [ncells])
    call create_shared_mem(amr_grid%icell_of_leaf, [nleaf])

    if (mpar%h_rank == 0) then
      amr_grid%parent(1:ncells)        = t_parent(1:ncells)
      amr_grid%children(:, 1:ncells)   = t_children(:, 1:ncells)
      amr_grid%level(1:ncells)         = t_level(1:ncells)
      amr_grid%ileaf(1:ncells)         = t_ileaf(1:ncells)
      amr_grid%cx(1:ncells)            = t_cx(1:ncells)
      amr_grid%cy(1:ncells)            = t_cy(1:ncells)
      amr_grid%cz(1:ncells)            = t_cz(1:ncells)
      amr_grid%ch(1:ncells)            = t_ch(1:ncells)
      amr_grid%icell_of_leaf(1:nleaf)  = t_icell_of_leaf(1:nleaf)
    end if
    call MPI_BARRIER(mpar%hostcomm, ierr)

    deallocate(t_parent, t_children, t_level, t_ileaf, &
               t_cx, t_cy, t_cz, t_ch, t_icell_of_leaf)
  end subroutine amr_build_tree

  !=========================================================================
  ! Allocate physical property arrays as shared memory.
  ! h_rank==0 must fill the arrays and call MPI_BARRIER(mpar%hostcomm)
  ! after this subroutine returns.
  !=========================================================================
  subroutine amr_alloc_phys(use_dust)
    logical, intent(in) :: use_dust
    integer :: n
    n = amr_grid%nleaf
    call create_shared_mem(amr_grid%rhokap,  [n])
    call create_shared_mem(amr_grid%Dfreq,   [n])
    call create_shared_mem(amr_grid%voigt_a, [n])
    call create_shared_mem(amr_grid%vfx,     [n])
    call create_shared_mem(amr_grid%vfy,     [n])
    call create_shared_mem(amr_grid%vfz,     [n])
    if (use_dust) call create_shared_mem(amr_grid%rhokapD, [n])
  end subroutine amr_alloc_phys

  !=========================================================================
  ! Allocate spectral-output arrays (per-rank, NOT shared).
  !=========================================================================
  subroutine amr_alloc_output(save_jin, save_jabs, use_dust, save_jmu, nmu)
    logical, intent(in) :: save_jin, save_jabs, use_dust
    logical, intent(in), optional :: save_jmu
    integer, intent(in), optional :: nmu
    integer :: n
    n = amr_grid%nxfreq
    allocate(amr_grid%Jout(n), source=0.0_wp)
    if (save_jin)                 allocate(amr_grid%Jin(n),  source=0.0_wp)
    if (save_jabs .and. use_dust) allocate(amr_grid%Jabs(n), source=0.0_wp)
    if (present(save_jmu)) then
       if (save_jmu .and. present(nmu)) then
          allocate(amr_grid%Jmu(n, nmu), source=0.0_wp)
       endif
    endif
    !--- Ly-beta (line_type = 8): band-2 (H-alpha) escaped/absorbed spectra.
    !--- nxfreq_Ha is finalized in amr_setup_freq_grid (called before this).
    if (line%line_type == 8) then
       allocate(amr_grid%Jout_Ha(amr_grid%nxfreq_Ha), source=0.0_wp)
       if (save_jabs .and. use_dust) &
          allocate(amr_grid%Jabs_Ha(amr_grid%nxfreq_Ha), source=0.0_wp)
    endif
  end subroutine amr_alloc_output

  !=========================================================================
  ! Precompute the face-neighbor table neighbor(6, ncells) as shared memory.
  ! Must be called after amr_build_tree.
  !=========================================================================
  subroutine amr_build_neighbors
    use mpi
    integer  :: icell, ierr, iface
    real(wp) :: cx, cy, cz, h, hp
    integer  :: ineigh

    call create_shared_mem(amr_grid%neighbor, [6, amr_grid%ncells])

    if (mpar%h_rank == 0) then
      amr_grid%neighbor = 0
      do icell = 1, amr_grid%ncells
        cx  = amr_grid%cx(icell)
        cy  = amr_grid%cy(icell)
        cz  = amr_grid%cz(icell)
        h   = amr_grid%ch(icell)
        ! Query the would-be same-level neighbor's CENTER (1 cell width past icell's face),
        ! NOT just past the face itself.  Querying at cx + h lands EXACTLY on the parent's
        ! +x face when icell is on the +x side of its parent, and amr_find_cell_at_level's
        ! `>=` descent then routes the query back into icell (self-loop in the neighbor
        ! table).  Querying at cx + 2*h is unambiguously inside the +x neighbor.
        hp  = 2.0_wp * h

        if (cx + hp <= amr_grid%xmax) &
          amr_grid%neighbor(1, icell) = amr_find_cell_at_level(cx+hp, cy,   cz,   amr_grid%level(icell))
        if (cx - hp >= amr_grid%xmin) &
          amr_grid%neighbor(2, icell) = amr_find_cell_at_level(cx-hp, cy,   cz,   amr_grid%level(icell))
        if (cy + hp <= amr_grid%ymax) &
          amr_grid%neighbor(3, icell) = amr_find_cell_at_level(cx,    cy+hp, cz,  amr_grid%level(icell))
        if (cy - hp >= amr_grid%ymin) &
          amr_grid%neighbor(4, icell) = amr_find_cell_at_level(cx,    cy-hp, cz,  amr_grid%level(icell))
        if (cz + hp <= amr_grid%zmax) &
          amr_grid%neighbor(5, icell) = amr_find_cell_at_level(cx,    cy,   cz+hp, amr_grid%level(icell))
        if (cz - hp >= amr_grid%zmin) &
          amr_grid%neighbor(6, icell) = amr_find_cell_at_level(cx,    cy,   cz-hp, amr_grid%level(icell))

        ! If amr_find_cell_at_level returned an ANCESTOR of icell, the descent
        ! ran out of refinement (sparse-data octree gap).  Such an ancestor would
        ! route a descent in amr_next_leaf back into icell itself (self-loop
        ! equivalent).  Mark these as "no neighbor" so the raytrace treats them
        ! as escapes.
        do iface = 1, 6
          ineigh = amr_grid%neighbor(iface, icell)
          if (ineigh > 0 .and. ineigh /= icell .and. is_ancestor(ineigh, icell)) then
            amr_grid%neighbor(iface, icell) = 0
          end if
        end do
      end do

      block
        integer :: nself, iface_dbg
        nself = 0
        do icell = 1, amr_grid%ncells
          do iface_dbg = 1, 6
            if (amr_grid%neighbor(iface_dbg, icell) == icell) nself = nself + 1
          end do
        end do
        if (nself /= 0) then
          write(6,'(a,i0)') 'amr_build_neighbors WARNING: self-loops detected: ', nself
          flush(6)
        end if
      end block
    end if
    call MPI_BARRIER(mpar%hostcomm, ierr)
  end subroutine amr_build_neighbors

  ! Return .true. if 'anc' is an ancestor of 'desc' in the octree.
  pure logical function is_ancestor(anc, desc)
    integer, intent(in) :: anc, desc
    integer :: c
    is_ancestor = .false.
    c = desc
    do while (c > 0)
      c = amr_grid%parent(c)
      if (c == anc) then
        is_ancestor = .true.
        return
      end if
    end do
  end function is_ancestor

  !=========================================================================
  ! Return the leaf index a photon enters after crossing face iface of
  ! cell icell.  (x, y, z) is the photon position AT the face (code units),
  ! used ONLY for the two transverse sub-octant bits; the iface-normal
  ! sub-octant bit is determined topologically from iface itself so that
  ! a photon sitting exactly on the face cannot be re-routed by floating-
  ! point round-off in the normal coordinate.
  !
  ! Face convention: 1=+x, 2=-x, 3=+y, 4=-y, 5=+z, 6=-z (source cell).
  ! The photon enters the destination cell through the opposite face, so:
  !   iface=1 -> dest -x side -> x sub-octant bit = 0
  !   iface=2 -> dest +x side -> x sub-octant bit = 1
  !   iface=3 -> dest -y side -> y sub-octant bit = 0
  !   iface=4 -> dest +y side -> y sub-octant bit = 1
  !   iface=5 -> dest -z side -> z sub-octant bit = 0
  !   iface=6 -> dest +z side -> z sub-octant bit = 1
  !
  ! Returns 0 if the photon has left the box.
  !=========================================================================
  integer function amr_next_leaf(icell, iface, x, y, z) result(il_new)
    integer,  intent(in) :: icell, iface
    real(wp), intent(in) :: x, y, z
    integer :: ineigh, child, ioct

    il_new = 0
    ineigh = amr_grid%neighbor(iface, icell)
    if (ineigh == 0) return   ! outside box

    do while (amr_grid%ileaf(ineigh) == 0)
      ioct = 1
      select case (iface)
      case (1)                                                    ! dest -x side
        if (y >= amr_grid%cy(ineigh)) ioct = ioct + 2
        if (z >= amr_grid%cz(ineigh)) ioct = ioct + 4
      case (2)                                                    ! dest +x side
        ioct = ioct + 1
        if (y >= amr_grid%cy(ineigh)) ioct = ioct + 2
        if (z >= amr_grid%cz(ineigh)) ioct = ioct + 4
      case (3)                                                    ! dest -y side
        if (x >= amr_grid%cx(ineigh)) ioct = ioct + 1
        if (z >= amr_grid%cz(ineigh)) ioct = ioct + 4
      case (4)                                                    ! dest +y side
        if (x >= amr_grid%cx(ineigh)) ioct = ioct + 1
        ioct = ioct + 2
        if (z >= amr_grid%cz(ineigh)) ioct = ioct + 4
      case (5)                                                    ! dest -z side
        if (x >= amr_grid%cx(ineigh)) ioct = ioct + 1
        if (y >= amr_grid%cy(ineigh)) ioct = ioct + 2
      case (6)                                                    ! dest +z side
        if (x >= amr_grid%cx(ineigh)) ioct = ioct + 1
        if (y >= amr_grid%cy(ineigh)) ioct = ioct + 2
        ioct = ioct + 4
      end select
      child = amr_grid%children(ioct, ineigh)
      if (child == 0) exit
      ineigh = child
    end do

    il_new = amr_grid%ileaf(ineigh)
  end function amr_next_leaf

  !=========================================================================
  ! Same as amr_next_leaf but tells the caller whether a returned il_new==0
  ! is "outside the box" or "inside a gap" (internal cell with no child at
  ! the relevant octant).  When in a gap, icell_gap and h_gap give the
  ! internal cell index and its half-width.
  !
  ! Output codes:
  !   il_new  > 0 :  valid leaf found
  !   il_new == 0 .and. icell_gap > 0 :  in a gap; h_gap = gap-cell half-width
  !   il_new == 0 .and. icell_gap == 0:  truly outside the box
  !=========================================================================
  subroutine amr_next_leaf_or_gap(icell, iface, x, y, z, il_new, icell_gap, h_gap)
    integer,  intent(in)  :: icell, iface
    real(wp), intent(in)  :: x, y, z
    integer,  intent(out) :: il_new, icell_gap
    real(wp), intent(out) :: h_gap
    integer :: ineigh, child, ioct

    il_new    = 0
    icell_gap = 0
    h_gap     = 0.0_wp
    ineigh = amr_grid%neighbor(iface, icell)
    if (ineigh == 0) return   ! truly outside box

    do while (amr_grid%ileaf(ineigh) == 0)
      ioct = 1
      select case (iface)
      case (1)
        if (y >= amr_grid%cy(ineigh)) ioct = ioct + 2
        if (z >= amr_grid%cz(ineigh)) ioct = ioct + 4
      case (2)
        ioct = ioct + 1
        if (y >= amr_grid%cy(ineigh)) ioct = ioct + 2
        if (z >= amr_grid%cz(ineigh)) ioct = ioct + 4
      case (3)
        if (x >= amr_grid%cx(ineigh)) ioct = ioct + 1
        if (z >= amr_grid%cz(ineigh)) ioct = ioct + 4
      case (4)
        if (x >= amr_grid%cx(ineigh)) ioct = ioct + 1
        ioct = ioct + 2
        if (z >= amr_grid%cz(ineigh)) ioct = ioct + 4
      case (5)
        if (x >= amr_grid%cx(ineigh)) ioct = ioct + 1
        if (y >= amr_grid%cy(ineigh)) ioct = ioct + 2
      case (6)
        if (x >= amr_grid%cx(ineigh)) ioct = ioct + 1
        if (y >= amr_grid%cy(ineigh)) ioct = ioct + 2
        ioct = ioct + 4
      end select
      child = amr_grid%children(ioct, ineigh)
      if (child == 0) then
        icell_gap = ineigh
        h_gap     = amr_grid%ch(ineigh)
        return
      end if
      ineigh = child
    end do

    il_new = amr_grid%ileaf(ineigh)
  end subroutine amr_next_leaf_or_gap

  ! ---- private helpers ----

  !=========================================================================
  ! Descend the octree to find the cell at exactly `target_level` that
  ! contains (x, y, z).
  !=========================================================================
  integer function amr_find_cell_at_level(x, y, z, target_level) result(icell)
    real(wp), intent(in) :: x, y, z
    integer,  intent(in) :: target_level
    integer :: ioct, child

    icell = 0
    if (x < amr_grid%xmin .or. x > amr_grid%xmax .or. &
        y < amr_grid%ymin .or. y > amr_grid%ymax .or. &
        z < amr_grid%zmin .or. z > amr_grid%zmax) return

    icell = 1
    do
      if (amr_grid%level(icell) >= target_level) return
      if (amr_grid%ileaf(icell) > 0)             return
      ioct  = 1
      if (x >= amr_grid%cx(icell)) ioct = ioct + 1
      if (y >= amr_grid%cy(icell)) ioct = ioct + 2
      if (z >= amr_grid%cz(icell)) ioct = ioct + 4
      child = amr_grid%children(ioct, icell)
      if (child == 0) return
      icell = child
    end do
  end function amr_find_cell_at_level

  !=========================================================================
  ! Grow local allocatable build arrays to new_max elements.
  !=========================================================================
  subroutine grow_local_arrays(new_max, old_max, &
      t_par, t_chi, t_lev, t_ile, t_cx, t_cy, t_cz, t_ch)
    integer,  intent(in)    :: new_max, old_max
    integer,  allocatable, intent(inout) :: t_par(:), t_chi(:,:)
    integer,  allocatable, intent(inout) :: t_lev(:), t_ile(:)
    real(wp), allocatable, intent(inout) :: t_cx(:), t_cy(:), t_cz(:), t_ch(:)

    integer,  allocatable :: tmp_i1(:), tmp_i2(:,:)
    real(wp), allocatable :: tmp_r(:)

    allocate(tmp_i1(new_max))
    tmp_i1(1:old_max) = t_par;  tmp_i1(old_max+1:) = 0
    call move_alloc(tmp_i1, t_par)

    allocate(tmp_i2(8, new_max))
    tmp_i2(:, 1:old_max) = t_chi;  tmp_i2(:, old_max+1:) = 0
    call move_alloc(tmp_i2, t_chi)

    allocate(tmp_i1(new_max))
    tmp_i1(1:old_max) = t_lev;  tmp_i1(old_max+1:) = 0
    call move_alloc(tmp_i1, t_lev)

    allocate(tmp_i1(new_max))
    tmp_i1(1:old_max) = t_ile;  tmp_i1(old_max+1:) = 0
    call move_alloc(tmp_i1, t_ile)

    allocate(tmp_r(new_max))
    tmp_r(1:old_max) = t_cx;  tmp_r(old_max+1:) = 0.0_wp
    call move_alloc(tmp_r, t_cx)

    allocate(tmp_r(new_max))
    tmp_r(1:old_max) = t_cy;  tmp_r(old_max+1:) = 0.0_wp
    call move_alloc(tmp_r, t_cy)

    allocate(tmp_r(new_max))
    tmp_r(1:old_max) = t_cz;  tmp_r(old_max+1:) = 0.0_wp
    call move_alloc(tmp_r, t_cz)

    allocate(tmp_r(new_max))
    tmp_r(1:old_max) = t_ch;  tmp_r(old_max+1:) = 0.0_wp
    call move_alloc(tmp_r, t_ch)
  end subroutine grow_local_arrays

end module octree_mod
