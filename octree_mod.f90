module octree_mod
  use define
  implicit none
  public

  !-----------------------------------------------------------------------
  ! AMR grid type: octree data structure for Adaptive Mesh Refinement.
  !
  ! Tree storage is a flat array of cells (internal + leaf), 1-indexed.
  ! Root = cell 1.  Octant ordering: child index = 1 + ix + 2*iy + 4*iz
  ! where ix=1 if x >= cell_center_x, else 0 (similarly for iy, iz).
  !-----------------------------------------------------------------------
  type amr_grid_type

    ! ------ domain geometry ------
    real(wp) :: xmin = 0.0_wp, xmax = 1.0_wp
    real(wp) :: ymin = 0.0_wp, ymax = 1.0_wp
    real(wp) :: zmin = 0.0_wp, zmax = 1.0_wp
    real(wp) :: L_box = 1.0_wp    ! box length (same units as xmin..xmax)

    ! ------ octree structure ------
    integer  :: ncells   = 0      ! total cells (internal + leaf)
    integer  :: nleaf    = 0      ! number of leaf cells
    integer  :: levelmax = 0      ! maximum AMR level present

    ! Arrays of size ncells (tree topology)
    integer,  allocatable :: parent(:)      ! parent index; 0 for root
    integer,  allocatable :: children(:,:)  ! (8, ncells): child cell indices; 0 = no child
    integer,  allocatable :: level(:)       ! cell AMR level (0 = root)
    real(wp), allocatable :: cx(:), cy(:), cz(:) ! cell center coordinates
    real(wp), allocatable :: ch(:)          ! cell half-width

    ! Leaf-index mapping
    integer,  allocatable :: ileaf(:)           ! ncells: >0 is leaf index; 0 = internal
    integer,  allocatable :: icell_of_leaf(:)   ! nleaf → cell index

    ! ------ face-neighbor lookup table (size 6 × ncells) ------
    ! neighbor(iface, icell) = cell index of the face neighbor in direction iface:
    !   iface 1=+x  2=-x  3=+y  4=-y  5=+z  6=-z
    ! 0  → outside the box.
    ! Leaf cell  (ileaf > 0) → direct O(1) result.
    ! Internal cell (ileaf = 0, neighbor is at finer level) → descent at runtime
    !   using the actual photon position; see amr_next_leaf.
    integer, allocatable :: neighbor(:,:)  ! (6, ncells)

    ! ------ physical data at leaf cells (size nleaf) ------
    real(wp), allocatable :: rhokap(:)   ! HI opacity per unit length at line centre
    real(wp), allocatable :: rhokapD(:)  ! dust opacity per unit length
    real(wp), allocatable :: Dfreq(:)    ! local Doppler frequency
    real(wp), allocatable :: voigt_a(:)  ! Voigt damping parameter a = A21/(4pi*Dfreq)
    real(wp), allocatable :: vfx(:)      ! x-velocity / v_thermal(local)
    real(wp), allocatable :: vfy(:)      ! y-velocity / v_thermal(local)
    real(wp), allocatable :: vfz(:)      ! z-velocity / v_thermal(local)

    ! ------ frequency / spectral grid ------
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

    real(wp), allocatable :: xfreq(:), velocity(:), wavelength(:)

    ! ------ spectral output arrays (per frequency bin) ------
    real(wp), allocatable :: Jout(:)  ! escaped spectrum
    real(wp), allocatable :: Jin(:)   ! input (injected) spectrum
    real(wp), allocatable :: Jabs(:)  ! absorbed by dust spectrum

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
  ! Input:
  !   xleaf, yleaf, zleaf  -- leaf cell centre coordinates (same units as box)
  !   leaf_level           -- AMR refinement level of each leaf (0 = root)
  !   nleaf                -- number of leaf cells
  !   box_*                -- domain extent (physical units, arbitrary but consistent)
  !
  ! Post-condition: amr_grid%{parent,children,level,cx,cy,cz,ch,ileaf,icell_of_leaf}
  !                 are allocated and filled.  Physical data arrays are NOT yet
  !                 allocated; call amr_alloc_phys after this.
  !=========================================================================
  subroutine amr_build_tree(xleaf, yleaf, zleaf, leaf_level, nleaf, &
                              box_xmin, box_xmax, box_ymin, box_ymax, box_zmin, box_zmax)
    integer,  intent(in) :: nleaf
    real(wp), intent(in) :: xleaf(nleaf), yleaf(nleaf), zleaf(nleaf)
    integer,  intent(in) :: leaf_level(nleaf)
    real(wp), intent(in) :: box_xmin, box_xmax
    real(wp), intent(in) :: box_ymin, box_ymax
    real(wp), intent(in) :: box_zmin, box_zmax

    integer  :: ncells_max, ncells
    integer  :: il, icell, lev, ioct, nc, ix, iy, iz

    ncells_max = max(nleaf * 2, 16)
    call alloc_tree_arrays(ncells_max, nleaf)

    amr_grid%xmin  = box_xmin;  amr_grid%xmax  = box_xmax
    amr_grid%ymin  = box_ymin;  amr_grid%ymax  = box_ymax
    amr_grid%zmin  = box_zmin;  amr_grid%zmax  = box_zmax
    amr_grid%L_box = box_xmax - box_xmin

    ! Initialise root cell (index 1)
    ncells                   = 1
    amr_grid%parent(1)       = 0
    amr_grid%children(:, 1)  = 0
    amr_grid%level(1)        = 0
    amr_grid%ileaf(1)        = 0
    amr_grid%cx(1)           = (box_xmin + box_xmax) * 0.5_wp
    amr_grid%cy(1)           = (box_ymin + box_ymax) * 0.5_wp
    amr_grid%cz(1)           = (box_zmin + box_zmax) * 0.5_wp
    amr_grid%ch(1)           = amr_grid%L_box * 0.5_wp

    ! Insert each leaf, creating internal nodes as needed
    do il = 1, nleaf
      icell = 1
      do lev = 0, leaf_level(il) - 1
        ioct = 1
        if (xleaf(il) >= amr_grid%cx(icell)) ioct = ioct + 1
        if (yleaf(il) >= amr_grid%cy(icell)) ioct = ioct + 2
        if (zleaf(il) >= amr_grid%cz(icell)) ioct = ioct + 4
        if (amr_grid%children(ioct, icell) == 0) then
          ncells = ncells + 1
          if (ncells > ncells_max) then
            ncells_max = ncells_max * 2
            call grow_tree_arrays(ncells_max)
          end if
          nc = ncells
          amr_grid%parent(nc)      = icell
          amr_grid%children(:, nc) = 0
          amr_grid%level(nc)       = lev + 1
          amr_grid%ileaf(nc)       = 0
          amr_grid%ch(nc)          = amr_grid%ch(icell) * 0.5_wp
          ix = mod(ioct - 1, 2)
          iy = mod((ioct - 1) / 2, 2)
          iz = (ioct - 1) / 4
          amr_grid%cx(nc) = amr_grid%cx(icell) + real(2*ix - 1, wp) * amr_grid%ch(nc)
          amr_grid%cy(nc) = amr_grid%cy(icell) + real(2*iy - 1, wp) * amr_grid%ch(nc)
          amr_grid%cz(nc) = amr_grid%cz(icell) + real(2*iz - 1, wp) * amr_grid%ch(nc)
          amr_grid%children(ioct, icell) = nc
        end if
        icell = amr_grid%children(ioct, icell)
      end do
      ! icell now sits at the target level and is the leaf cell
      amr_grid%ileaf(icell)       = il
      amr_grid%icell_of_leaf(il)  = icell
    end do

    amr_grid%ncells   = ncells
    amr_grid%nleaf    = nleaf
    amr_grid%levelmax = maxval(leaf_level)
  end subroutine amr_build_tree

  !=========================================================================
  ! Allocate physical property arrays (rhokap, Dfreq, voigt_a, vf*, rhokapD).
  !=========================================================================
  subroutine amr_alloc_phys(use_dust)
    logical, intent(in) :: use_dust
    integer :: n
    n = amr_grid%nleaf
    allocate(amr_grid%rhokap(n),   source=0.0_wp)
    allocate(amr_grid%Dfreq(n),    source=amr_grid%Dfreq_ref)
    allocate(amr_grid%voigt_a(n),  source=amr_grid%voigt_amean)
    allocate(amr_grid%vfx(n),      source=0.0_wp)
    allocate(amr_grid%vfy(n),      source=0.0_wp)
    allocate(amr_grid%vfz(n),      source=0.0_wp)
    if (use_dust) allocate(amr_grid%rhokapD(n), source=0.0_wp)
  end subroutine amr_alloc_phys

  !=========================================================================
  ! Allocate spectral-output arrays (Jout, Jin, Jabs).
  !=========================================================================
  subroutine amr_alloc_output(save_jin, save_jabs, use_dust)
    logical, intent(in) :: save_jin, save_jabs, use_dust
    integer :: n
    n = amr_grid%nxfreq
    allocate(amr_grid%Jout(n), source=0.0_wp)
    if (save_jin)                  allocate(amr_grid%Jin(n),  source=0.0_wp)
    if (save_jabs .and. use_dust)  allocate(amr_grid%Jabs(n), source=0.0_wp)
  end subroutine amr_alloc_output

  !=========================================================================
  ! Precompute the face-neighbor table neighbor(6, ncells).
  ! Must be called after amr_build_tree (so cx,cy,cz,ch,level,ileaf are set).
  ! Stored neighbors are at the SAME level or coarser. If the neighbor is
  ! an internal cell (ileaf==0) the photon crossed into a finer region;
  ! amr_next_leaf descends to the correct leaf using the crossing position.
  !=========================================================================
  subroutine amr_build_neighbors
    integer  :: icell, lev
    real(wp) :: cx, cy, cz, h, hp

    if (allocated(amr_grid%neighbor)) deallocate(amr_grid%neighbor)
    allocate(amr_grid%neighbor(6, amr_grid%ncells), source=0)

    do icell = 1, amr_grid%ncells
      cx  = amr_grid%cx(icell)
      cy  = amr_grid%cy(icell)
      cz  = amr_grid%cz(icell)
      h   = amr_grid%ch(icell)
      lev = amr_grid%level(icell)
      ! Probe just past each face: 1e-8 * h is > fp-epsilon at scale h
      ! yet << h, so it reliably lands inside the neighboring cell.
      hp  = h * (1.0_wp + 1.0e-8_wp)

      if (cx + hp <= amr_grid%xmax) &
        amr_grid%neighbor(1, icell) = amr_find_cell_at_level(cx+hp, cy,   cz,   lev)
      if (cx - hp >= amr_grid%xmin) &
        amr_grid%neighbor(2, icell) = amr_find_cell_at_level(cx-hp, cy,   cz,   lev)
      if (cy + hp <= amr_grid%ymax) &
        amr_grid%neighbor(3, icell) = amr_find_cell_at_level(cx,    cy+hp, cz,  lev)
      if (cy - hp >= amr_grid%ymin) &
        amr_grid%neighbor(4, icell) = amr_find_cell_at_level(cx,    cy-hp, cz,  lev)
      if (cz + hp <= amr_grid%zmax) &
        amr_grid%neighbor(5, icell) = amr_find_cell_at_level(cx,    cy,   cz+hp, lev)
      if (cz - hp >= amr_grid%zmin) &
        amr_grid%neighbor(6, icell) = amr_find_cell_at_level(cx,    cy,   cz-hp, lev)
    end do
  end subroutine amr_build_neighbors

  !=========================================================================
  ! Return the leaf index a photon enters after crossing face iface of
  ! cell icell.  (x, y, z) is the photon position AT the face (code units).
  ! Returns 0 if the photon has left the box.
  !
  ! Algorithm:
  !   1. Look up neighbor(iface, icell) from the precomputed table (O(1)).
  !   2. If the result is an internal cell (finer refinement), descend using
  !      the actual crossing position until a leaf is found (O(Δlevel)).
  !=========================================================================
  integer function amr_next_leaf(icell, iface, x, y, z) result(il_new)
    integer,  intent(in) :: icell, iface
    real(wp), intent(in) :: x, y, z
    integer :: ineigh, child, ioct

    il_new = 0
    ineigh = amr_grid%neighbor(iface, icell)
    if (ineigh == 0) return   ! outside box

    ! Descend into finer cells if the neighbor is an internal cell.
    do while (amr_grid%ileaf(ineigh) == 0)
      ioct  = 1
      if (x >= amr_grid%cx(ineigh)) ioct = ioct + 1
      if (y >= amr_grid%cy(ineigh)) ioct = ioct + 2
      if (z >= amr_grid%cz(ineigh)) ioct = ioct + 4
      child = amr_grid%children(ioct, ineigh)
      if (child == 0) exit   ! grid gap (should not happen in a valid grid)
      ineigh = child
    end do

    il_new = amr_grid%ileaf(ineigh)
  end function amr_next_leaf

  ! ---- private helpers ----

  !=========================================================================
  ! Descend the octree to find the cell at exactly `target_level` that
  ! contains (x, y, z).  Returns the coarsest leaf if the grid is not
  ! refined to target_level at that position.  Returns 0 if outside box.
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
      if (amr_grid%level(icell) >= target_level) return  ! reached target level
      if (amr_grid%ileaf(icell) > 0)             return  ! leaf at coarser level
      ioct  = 1
      if (x >= amr_grid%cx(icell)) ioct = ioct + 1
      if (y >= amr_grid%cy(icell)) ioct = ioct + 2
      if (z >= amr_grid%cz(icell)) ioct = ioct + 4
      child = amr_grid%children(ioct, icell)
      if (child == 0) return   ! grid gap
      icell = child
    end do
  end function amr_find_cell_at_level

  subroutine alloc_tree_arrays(ncells_max, nleaf)
    integer, intent(in) :: ncells_max, nleaf
    allocate(amr_grid%parent(ncells_max))
    allocate(amr_grid%children(8, ncells_max))
    allocate(amr_grid%level(ncells_max))
    allocate(amr_grid%cx(ncells_max), amr_grid%cy(ncells_max), amr_grid%cz(ncells_max))
    allocate(amr_grid%ch(ncells_max))
    allocate(amr_grid%ileaf(ncells_max))
    amr_grid%parent   = 0
    amr_grid%children = 0
    amr_grid%level    = 0
    amr_grid%ileaf    = 0
    amr_grid%cx       = 0.0_wp
    amr_grid%cy       = 0.0_wp
    amr_grid%cz       = 0.0_wp
    amr_grid%ch       = 0.0_wp
    allocate(amr_grid%icell_of_leaf(nleaf))
    amr_grid%icell_of_leaf = 0
  end subroutine alloc_tree_arrays

  subroutine grow_tree_arrays(new_max)
    integer, intent(in) :: new_max
    integer  :: n_old
    integer,  allocatable :: tmp_i1(:), tmp_i2(:,:)
    real(wp), allocatable :: tmp_r(:)

    n_old = size(amr_grid%parent)

    allocate(tmp_i1(new_max))
    tmp_i1(1:n_old)     = amr_grid%parent;    tmp_i1(n_old+1:) = 0
    call move_alloc(tmp_i1, amr_grid%parent)

    allocate(tmp_i2(8, new_max))
    tmp_i2(:, 1:n_old)  = amr_grid%children;  tmp_i2(:, n_old+1:) = 0
    call move_alloc(tmp_i2, amr_grid%children)

    allocate(tmp_i1(new_max))
    tmp_i1(1:n_old)     = amr_grid%level;     tmp_i1(n_old+1:) = 0
    call move_alloc(tmp_i1, amr_grid%level)

    allocate(tmp_i1(new_max))
    tmp_i1(1:n_old)     = amr_grid%ileaf;     tmp_i1(n_old+1:) = 0
    call move_alloc(tmp_i1, amr_grid%ileaf)

    allocate(tmp_r(new_max))
    tmp_r(1:n_old) = amr_grid%cx;  tmp_r(n_old+1:) = 0.0_wp
    call move_alloc(tmp_r, amr_grid%cx)

    allocate(tmp_r(new_max))
    tmp_r(1:n_old) = amr_grid%cy;  tmp_r(n_old+1:) = 0.0_wp
    call move_alloc(tmp_r, amr_grid%cy)

    allocate(tmp_r(new_max))
    tmp_r(1:n_old) = amr_grid%cz;  tmp_r(n_old+1:) = 0.0_wp
    call move_alloc(tmp_r, amr_grid%cz)

    allocate(tmp_r(new_max))
    tmp_r(1:n_old) = amr_grid%ch;  tmp_r(n_old+1:) = 0.0_wp
    call move_alloc(tmp_r, amr_grid%ch)
  end subroutine grow_tree_arrays

end module octree_mod
