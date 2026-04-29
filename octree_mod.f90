module octree_mod
  use define
  use memory_mod
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
    real(wp), pointer :: rhokap(:)   => null()  ! HI opacity per unit length at line centre
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
  end subroutine amr_alloc_output

  !=========================================================================
  ! Precompute the face-neighbor table neighbor(6, ncells) as shared memory.
  ! Must be called after amr_build_tree.
  !=========================================================================
  subroutine amr_build_neighbors
    use mpi
    integer  :: icell, ierr
    real(wp) :: cx, cy, cz, h, hp

    call create_shared_mem(amr_grid%neighbor, [6, amr_grid%ncells])

    if (mpar%h_rank == 0) then
      amr_grid%neighbor = 0
      do icell = 1, amr_grid%ncells
        cx  = amr_grid%cx(icell)
        cy  = amr_grid%cy(icell)
        cz  = amr_grid%cz(icell)
        h   = amr_grid%ch(icell)
        hp  = h * (1.0_wp + 1.0e-8_wp)

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
      end do
    end if
    call MPI_BARRIER(mpar%hostcomm, ierr)
  end subroutine amr_build_neighbors

  !=========================================================================
  ! Return the leaf index a photon enters after crossing face iface of
  ! cell icell.  (x, y, z) is the photon position AT the face (code units).
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
      ioct  = 1
      if (x >= amr_grid%cx(ineigh)) ioct = ioct + 1
      if (y >= amr_grid%cy(ineigh)) ioct = ioct + 2
      if (z >= amr_grid%cz(ineigh)) ioct = ioct + 4
      child = amr_grid%children(ioct, ineigh)
      if (child == 0) exit
      ineigh = child
    end do

    il_new = amr_grid%ileaf(ineigh)
  end function amr_next_leaf

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
