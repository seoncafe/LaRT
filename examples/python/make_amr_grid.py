"""
make_amr_grid.py
================
Library for building generic-format AMR data files for LaRT v2.00.

Output file format (amr_type = 'generic')
------------------------------------------
    nleaf  boxlen
    x  y  z  level  gasDen[cm^-3]  T[K]  vx[km/s]  vy[km/s]  vz[km/s]
    ...

  - x, y, z   : leaf-cell centre coordinates in the same unit as boxlen
                 (origin at box corner, range [0, boxlen])
  - level      : AMR refinement level (root = 0; minimum leaf level = 1)
  - gasDen     : total hydrogen number density [cm^-3]
                 (neutral fraction applied in LaRT via par%use_cie_condition)
  - T          : gas temperature [K]
  - vx/vy/vz   : bulk gas velocity [km/s]

Corresponding LaRT input file
------------------------------
    par%use_amr_grid  = .true.
    par%amr_type      = 'generic'
    par%amr_file      = 'my_grid.dat'
    par%distance_unit = 'kpc'    ! unit of x,y,z and boxlen in the data file

Quick start
-----------
    from make_amr_grid import AMRGrid
    import numpy as np

    boxlen = 100.0                          # kpc
    grid   = AMRGrid(boxlen)

    # Refine cells inside a sphere of radius 40 kpc centred on the box
    grid.refine(lambda c: c.dist(50, 50, 50) < 40, level_max=3)

    # Assign physical properties
    def gasDen_fn(x, y, z):
        r = np.sqrt((x-50)**2 + (y-50)**2 + (z-50)**2)
        return 1e-3 * np.exp(-r**2 / (2*20.0**2))

    grid.set_density(gasDen_fn)
    grid.set_temperature(lambda x, y, z: 1e4)
    # velocity defaults to zero

    grid.write('my_grid.dat')
    print(grid.info())

Performance optimizations
--------------------------
The following optimizations have been applied to handle large grids
(level_max >= 6, millions of leaf cells) efficiently.

1. AMRGrid.read / _build_from_rows  (O(N·depth) tree insertion)
   Previous implementation called grid.refine() once per row in the data
   file.  refine() traverses the entire current tree, giving O(N·tree_size)
   ≈ O(N²) total cost — prohibitively slow for level_max = 7 with ~10⁶
   leaves.  Replaced with _insert_cell(), which navigates from the root
   directly to the target depth using the octree child-index formula
   (iz*4 + iy*2 + ix), splitting nodes on the way.  Cost per cell is
   O(depth) = O(level_max), so total cost is O(N·level_max).
   Physical data are stored on the leaf at insertion time, eliminating
   the separate O(N log N) sort-and-assign pass.

2. leaves()  (iterative stack traversal, no Python recursion)
   The original recursive _collect_leaves() incurred one Python function-
   call per tree node (~0.1–0.3 µs overhead each).  For millions of
   nodes that amounts to seconds of pure call overhead.  Replaced with an
   explicit stack loop (stack.extend / stack.pop), which avoids all
   function-call cost and is typically 2–3× faster on large trees.

3. set_density / set_temperature / set_velocity / set_properties
   (vectorised physics-function evaluation)
   The original code called fn(scalar, scalar, scalar) once per leaf
   cell.  For numpy-based functions (e.g. Gaussian profiles built with
   np.exp / np.sqrt) this creates per-call numpy overhead N times.
   The new implementation collects all leaf coordinates into numpy arrays
   (cx, cy, cz) and calls fn(cx_array, cy_array, cz_array) once.  numpy
   then evaluates the expression in a single C-level vectorised pass,
   giving 10–100× speedup for typical physics functions.  A try/except
   fallback to the scalar loop handles non-vectorisable callables.

4. _make_dens_criterion / _make_vel_criterion
   (vectorised probe-point evaluation)
   Each candidate leaf during refinement is probed at 8 (nprobe=2) or
   nprobe³ sub-cell centres.  Previously these were evaluated with a
   Python zip-loop (N_probe individual fn calls per leaf).  Now the
   function is probed once at closure creation to detect array support;
   if supported, all probe points are evaluated with a single fn(xx, yy,
   zz) array call, reducing per-leaf fn calls from N_probe to 1.

5. info()  (single-pass leaf traversal)
   The original info() called self.leaves() three times (once inside
   level_counts(), twice for gasDen/T statistics).  Replaced with a
   single pass that collects level counts and physical-data lists
   simultaneously.

6. _leaf_table_array  (column-wise structured-array filling)
   Filling a numpy structured array row-by-row with data[i] = tuple is
   slower than assigning whole columns at once (data['x'] = list).
   Switched to column-wise assignment, which lets numpy convert each
   Python list to a contiguous typed array in one C-level operation.

7. write (text format)  (batched output via np.savetxt)
   The original code called f.write(f'...') individually for every leaf,
   incurring per-line Python format overhead.  Replaced with np.column_stack
   to build a 2-D array of all data, then np.savetxt for bulk formatted
   output.

8. slice_plot  (tree-pruned leaf collection + numpy polygon array)
   The original code collected all leaves with self.leaves() and then
   filtered by the slice plane, visiting O(N_total) nodes regardless of
   how many intersect the slice.  Replaced with _collect_slice_leaves(),
   a recursive traversal that prunes any subtree whose bounding box does
   not intersect the slice plane, reducing traversal cost to
   O(N_slice·depth).  Polygon vertices are then built as a single
   (N, 4, 2) numpy array instead of a Python list-of-lists, eliminating
   per-cell Python object creation before passing to PolyCollection.

Geometry boundary refinement  (refine_boundary parameter)
-----------------------------------------------------------
The ``refine_sphere_adaptive``, ``refine_slab_adaptive``,
``refine_box_adaptive`` methods and their ``_by_physics`` variants all
accept a ``refine_boundary`` keyword argument (default ``False``).

When ``refine_boundary=False`` (default):
    Only the physics criterion (density and/or velocity gradient) drives
    refinement.  No special treatment is applied to the geometric boundary
    of the sphere/slab/box.  This is sufficient — and faster — when the
    density field itself captures the boundary (e.g. density drops sharply
    at the sphere edge so the gradient criterion naturally produces fine
    cells there).

    Internally this is equivalent to calling ``grid.refine(criterion,
    level_max)`` directly; the geometry shape plays no role.

When ``refine_boundary=True``:
    In addition to the physics criterion, every leaf cell that intersects
    the geometric surface is forced all the way to ``level_max``
    regardless of the local gradient.  Use this when the boundary is not
    captured by the density/velocity field alone (e.g. a sharp density
    step that lies exactly between two coarse cells, or a smooth field
    whose gradient falls below the threshold near the boundary).

    This was the only available behaviour before the ``refine_boundary``
    parameter was introduced.

Example usage::

    # Default: physics gradient only (recommended starting point)
    grid.refine_sphere_by_physics(50, 50, 50, radius=40, gasDen_fn=dens,
                                  level_max=7)

    # Forced boundary resolution on top of the physics criterion
    grid.refine_sphere_by_physics(50, 50, 50, radius=40, gasDen_fn=dens,
                                  level_max=7, refine_boundary=True)

Trade-off:
    ``refine_boundary=True`` guarantees a sharp geometric surface but
    typically generates ~20–40 % more leaf cells near the boundary,
    increasing both memory usage and LaRT run time.  Prefer
    ``refine_boundary=False`` unless you observe staircase artefacts at
    the geometry surface in your output.

Minimum refinement level  (level_min parameter)
-------------------------------------------------
All physics-based refinement methods (``refine_by_density``,
``refine_by_velocity``, ``refine_by_physics``, ``refine_sphere_by_physics``,
``refine_slab_by_physics``, ``refine_box_by_physics``) accept a
``level_min`` keyword argument (default ``2``).

Before applying the gradient criterion, the entire grid is uniformly split
down to ``level_min`` by calling ``refine_uniform(level_min)``.  This
guarantees a coarse baseline resolution everywhere in the box and prevents
the criterion from being evaluated only on the very coarse root cell (where
all probe points may happen to land on the same side of a feature, missing
it entirely).

    level_min=2  →  8² = 64 cells before the gradient scan starts
    level_min=3  →  8³ = 512 cells before the gradient scan starts
    level_min=0  →  no pre-refinement (original behaviour)

Set ``level_min=0`` to restore the old behaviour and let the gradient
criterion drive all splitting from the root.

The ``refine_uniform(level)`` method is also exposed as a public API and
can be called directly to establish any desired baseline resolution before
applying a custom ``refine()`` criterion::

    grid = AMRGrid(100.0)
    grid.refine_uniform(2)                        # 64 equal cells
    grid.refine(my_criterion, level_max=7)        # physics refinement on top
"""

import numpy as np

try:
    from astropy.io import fits
except ImportError:
    fits = None


# ---------------------------------------------------------------------------
# Cell
# ---------------------------------------------------------------------------

class Cell:
    """
    Single node in the octree.

    Parameters
    ----------
    cx, cy, cz : float
        Cell centre coordinates (same unit as boxlen).
    h : float
        Cell half-width (cell spans [c-h, c+h] along each axis).
    level : int
        AMR refinement level (root = 0).

    Attributes
    ----------
    children : list of Cell or None
        Eight children if refined; None if this is a leaf.
    gasDen, T, vx, vy, vz : float
        Physical data (only meaningful for leaves).
    """

    __slots__ = ('cx', 'cy', 'cz', 'h', 'level', 'children',
                 'gasDen', 'T', 'vx', 'vy', 'vz')

    def __init__(self, cx, cy, cz, h, level):
        self.cx    = cx
        self.cy    = cy
        self.cz    = cz
        self.h     = h
        self.level = level
        self.children = None
        # Physical defaults
        self.gasDen = 0.0
        self.T  = 1e4
        self.vx = 0.0
        self.vy = 0.0
        self.vz = 0.0

    # ------------------------------------------------------------------ #
    @property
    def is_leaf(self):
        return self.children is None

    def dist(self, x0, y0, z0):
        """Distance from cell centre to point (x0, y0, z0)."""
        return np.sqrt((self.cx - x0)**2 +
                       (self.cy - y0)**2 +
                       (self.cz - z0)**2)

    @property
    def xmin(self): return self.cx - self.h
    @property
    def xmax(self): return self.cx + self.h
    @property
    def ymin(self): return self.cy - self.h
    @property
    def ymax(self): return self.cy + self.h
    @property
    def zmin(self): return self.cz - self.h
    @property
    def zmax(self): return self.cz + self.h

    def corners(self):
        """Return the 8 cell-corner coordinates as an (8, 3) array."""
        offs = (-self.h, self.h)
        return np.array([
            (self.cx + dx, self.cy + dy, self.cz + dz)
            for dz in offs
            for dy in offs
            for dx in offs
        ], dtype=float)

    # ------------------------------------------------------------------ #
    def split(self):
        """Replace this leaf with 8 child leaves at half the size."""
        h2 = self.h * 0.5
        lv = self.level + 1
        self.children = [
            Cell(self.cx + (2*ix - 1)*h2,
                 self.cy + (2*iy - 1)*h2,
                 self.cz + (2*iz - 1)*h2,
                 h2, lv)
            for iz in (0, 1)
            for iy in (0, 1)
            for ix in (0, 1)
        ]

    def __repr__(self):
        tag = 'leaf' if self.is_leaf else f'{len(self.children)} children'
        return (f'Cell(lv={self.level}, '
                f'c=({self.cx:.2f},{self.cy:.2f},{self.cz:.2f}), '
                f'h={self.h:.2f}, {tag})')


# ---------------------------------------------------------------------------
# AMRGrid
# ---------------------------------------------------------------------------

class AMRGrid:
    """
    Octree AMR grid for LaRT v2.00.

    Parameters
    ----------
    boxlen : float
        Side length of the cubic domain (same unit as leaf coordinates).
    origin : tuple of float, optional
        Box corner coordinates (default (0, 0, 0)).  The root cell centre
        is at origin + (boxlen/2, boxlen/2, boxlen/2).

    Examples
    --------
    >>> grid = AMRGrid(100.0)
    >>> grid.refine(lambda c: c.dist(50, 50, 50) < 40, level_max=3)
    >>> grid.set_density(lambda x, y, z: 1e-3)
    >>> grid.write('uniform_100kpc.dat')
    """

    def __init__(self, boxlen, origin=(0.0, 0.0, 0.0)):
        self.boxlen = float(boxlen)
        self.origin = tuple(float(v) for v in origin)
        ox, oy, oz  = self.origin
        h           = self.boxlen * 0.5
        self.root   = Cell(ox + h, oy + h, oz + h, h, level=0)

    # ------------------------------------------------------------------ #
    # Refinement
    # ------------------------------------------------------------------ #

    def refine(self, criterion, level_max):
        """
        Recursively split leaf cells that satisfy ``criterion`` up to
        ``level_max``.

        Parameters
        ----------
        criterion : callable
            ``criterion(cell) -> bool``.  Receives a :class:`Cell` instance;
            should return True if the cell should be refined.
        level_max : int
            Cells at this level or deeper are never split.

        Notes
        -----
        Calling ``refine`` multiple times with different criteria is safe and
        additive (new criteria are applied to the current leaf set).
        """
        self._refine_recursive(self.root, criterion, level_max)

    def _refine_recursive(self, cell, criterion, level_max):
        if cell.is_leaf:
            if cell.level < level_max and criterion(cell):
                cell.split()
                for child in cell.children:
                    self._refine_recursive(child, criterion, level_max)
        else:
            for child in cell.children:
                self._refine_recursive(child, criterion, level_max)

    def refine_sphere(self, cx, cy, cz, radius, level_max):
        """
        Refine all cells whose centres lie within ``radius`` of
        (``cx``, ``cy``, ``cz``).
        """
        self.refine(lambda c: c.dist(cx, cy, cz) < radius, level_max)

    def refine_slab(self, axis, lo, hi, level_max):
        """
        Refine all cells whose centre along ``axis`` ('x', 'y', or 'z')
        falls in the interval [``lo``, ``hi``].
        """
        attr = {'x': 'cx', 'y': 'cy', 'z': 'cz'}[axis.lower()]
        self.refine(lambda c: lo <= getattr(c, attr) <= hi, level_max)

    def refine_box(self, xlo, xhi, ylo, yhi, zlo, zhi, level_max):
        """Refine all cells whose centre lies inside the given rectangular box."""
        self.refine(
            lambda c: (xlo <= c.cx <= xhi and
                       ylo <= c.cy <= yhi and
                       zlo <= c.cz <= zhi),
            level_max
        )

    def refine_uniform(self, level):
        """
        Unconditionally split all leaf cells below ``level``.

        Every leaf whose level is less than ``level`` is split, regardless
        of any physics criterion.  After this call the grid has at least
        ``8**level`` cells, each at level ≥ ``level``.

        This is called automatically by the physics-based refinement methods
        via their ``level_min`` parameter, but can also be invoked directly
        to establish a coarse baseline before applying a custom criterion.

        Parameters
        ----------
        level : int
            Target minimum level.  ``level=0`` is a no-op.
        """
        if level > 0:
            self.refine(lambda c: True, level)

    def refine_geometry(self, inside_fn, intersects_fn, level_max, criterion_inside=None):
        """
        Refine a geometry-aware region up to ``level_max``.

        Cells that intersect the geometry boundary are always refined until
        ``level_max``.  Cells that are fully contained inside the geometry
        are refined only when ``criterion_inside(cell)`` is True.  Cells
        completely outside the geometry are left untouched.

        Parameters
        ----------
        inside_fn : callable
            ``inside_fn(cell) -> bool``.  Must return True when the whole
            cell is contained in the geometry.
        intersects_fn : callable
            ``intersects_fn(cell) -> bool``.  Must return True when the cell
            overlaps the geometry at all (including fully-contained cells).
        level_max : int
            Maximum refinement level.
        criterion_inside : callable or None
            Optional refinement criterion for cells fully inside the geometry.
            If None, fully-contained cells are not further refined.
        """
        self._refine_geometry_recursive(
            self.root, inside_fn, intersects_fn, criterion_inside, level_max
        )

    def _refine_geometry_recursive(self, cell, inside_fn, intersects_fn,
                                   criterion_inside, level_max):
        if cell.is_leaf:
            if cell.level >= level_max:
                return

            fully_inside = inside_fn(cell)
            if fully_inside:
                if criterion_inside is None or not criterion_inside(cell):
                    return
            elif not intersects_fn(cell):
                return

            cell.split()
            for child in cell.children:
                self._refine_geometry_recursive(
                    child, inside_fn, intersects_fn, criterion_inside, level_max
                )
        else:
            for child in cell.children:
                self._refine_geometry_recursive(
                    child, inside_fn, intersects_fn, criterion_inside, level_max
                )

    def refine_sphere_adaptive(self, cx, cy, cz, radius, level_max,
                               criterion_inside=None, refine_boundary=False):
        """
        Refinement for a sphere, with optional geometry-forced boundary cells.

        Parameters
        ----------
        cx, cy, cz : float
            Centre of the sphere.
        radius : float
            Sphere radius (same unit as boxlen).
        level_max : int
            Maximum refinement level.
        criterion_inside : callable or None
            Additional refinement criterion applied to cells fully inside the
            sphere.  Ignored when ``refine_boundary=False``.
        refine_boundary : bool, optional (default False)
            If False (default), apply ``criterion_inside`` to all leaf cells
            with no special treatment of the geometric boundary — equivalent
            to calling ``refine(criterion_inside, level_max)``.
            If True, force every cell that intersects the sphere surface all
            the way to ``level_max`` regardless of the physics criterion, so
            that the geometric boundary is sharply resolved.
        """
        if not refine_boundary:
            if criterion_inside is not None:
                self.refine(criterion_inside, level_max)
            return

        r2 = radius**2

        def inside(c):
            corners = c.corners()
            dist2 = ((corners[:, 0] - cx)**2 +
                     (corners[:, 1] - cy)**2 +
                     (corners[:, 2] - cz)**2)
            return np.all(dist2 <= r2)

        def intersects(c):
            dx = max(c.xmin - cx, 0.0, cx - c.xmax)
            dy = max(c.ymin - cy, 0.0, cy - c.ymax)
            dz = max(c.zmin - cz, 0.0, cz - c.zmax)
            return dx*dx + dy*dy + dz*dz <= r2

        self.refine_geometry(inside, intersects, level_max, criterion_inside)

    def refine_slab_adaptive(self, axis, lo, hi, level_max,
                             criterion_inside=None, refine_boundary=False):
        """
        Refinement for a slab, with optional geometry-forced boundary cells.

        Parameters
        ----------
        axis : {'x', 'y', 'z'}
            Normal axis of the slab.
        lo, hi : float
            Slab extent along ``axis``.
        level_max : int
            Maximum refinement level.
        criterion_inside : callable or None
            Additional refinement criterion for cells fully inside the slab.
        refine_boundary : bool, optional (default False)
            If False (default), apply ``criterion_inside`` globally with no
            forced boundary refinement.
            If True, force every cell that overlaps the slab boundary to
            ``level_max``.
        """
        if not refine_boundary:
            if criterion_inside is not None:
                self.refine(criterion_inside, level_max)
            return

        axis = axis.lower()
        min_attr = {'x': 'xmin', 'y': 'ymin', 'z': 'zmin'}[axis]
        max_attr = {'x': 'xmax', 'y': 'ymax', 'z': 'zmax'}[axis]

        def inside(c):
            return lo <= getattr(c, min_attr) and getattr(c, max_attr) <= hi

        def intersects(c):
            return getattr(c, max_attr) >= lo and getattr(c, min_attr) <= hi

        self.refine_geometry(inside, intersects, level_max, criterion_inside)

    def refine_box_adaptive(self, xlo, xhi, ylo, yhi, zlo, zhi, level_max,
                            criterion_inside=None, refine_boundary=False):
        """
        Refinement for a box, with optional geometry-forced boundary cells.

        Parameters
        ----------
        xlo, xhi, ylo, yhi, zlo, zhi : float
            Box extents.
        level_max : int
            Maximum refinement level.
        criterion_inside : callable or None
            Additional refinement criterion for cells fully inside the box.
        refine_boundary : bool, optional (default False)
            If False (default), apply ``criterion_inside`` globally with no
            forced boundary refinement.
            If True, force every cell that overlaps the box boundary to
            ``level_max``.
        """
        if not refine_boundary:
            if criterion_inside is not None:
                self.refine(criterion_inside, level_max)
            return

        def inside(c):
            return (xlo <= c.xmin and c.xmax <= xhi and
                    ylo <= c.ymin and c.ymax <= yhi and
                    zlo <= c.zmin and c.zmax <= zhi)

        def intersects(c):
            return (c.xmax >= xlo and c.xmin <= xhi and
                    c.ymax >= ylo and c.ymin <= yhi and
                    c.zmax >= zlo and c.zmin <= zhi)

        self.refine_geometry(inside, intersects, level_max, criterion_inside)

    def refine_by_density(self, gasDen_fn, threshold=0.1, level_max=6,
                          level_min=2, nprobe=2, floor=1e-30):
        """
        Refine cells where the density gradient exceeds ``threshold``.

        Follows the RASCAS (module_mesh.f90) criterion:

            delta = (gasDen_max - gasDen_min) / (gasDen_max + gasDen_min + floor)

        The density is sampled at a uniform ``nprobe × nprobe × nprobe``
        sub-grid inside the cell plus the cell centre.  ``nprobe=2``
        samples the 8 cell corners (fast); ``nprobe=4`` (64 points) gives
        a more thorough probe.  With ``nprobe > 2`` the function also tries
        ``nprobe=2`` first and returns early if the gradient already exceeds
        the threshold (same early-exit logic as the Fortran reference).

        Parameters
        ----------
        gasDen_fn : callable
            ``gasDen_fn(x, y, z) -> gasDen [cm^-3]`` in the same coordinate units
            as the grid (e.g. kpc).
        threshold : float
            Gradient threshold ∈ [0, 1).  0.1 means 10 % relative variation.
        level_max : int
            Maximum refinement level.
        level_min : int, optional (default 2)
            Minimum refinement level.  All cells are unconditionally split to
            this level before the gradient criterion is applied.  This ensures
            a coarse baseline resolution everywhere and prevents the criterion
            from being evaluated only on the (very coarse) root cell.
            Set to 0 to disable uniform pre-refinement.
        nprobe : int
            Sub-grid sampling density per axis (≥ 2).  Default 2 (corners).
        floor : float
            Density floor to avoid division by zero.
        """
        self.refine_uniform(level_min)
        self.refine(
            self._make_dens_criterion(gasDen_fn, threshold, floor, nprobe),
            level_max
        )

    def refine_by_velocity(self, vel_fn, gasDen_fn=None, threshold=0.1,
                           level_max=6, level_min=2, nprobe=2, floor=1e-30):
        """
        Refine cells where the velocity magnitude gradient exceeds ``threshold``.

        Follows the RASCAS criterion applied to |v|:

            delta = (|v|_max - |v|_min) / (|v|_max + |v|_min + floor)

        Sampling is done at the same sub-grid as :meth:`refine_by_density`.
        If ``gasDen_fn`` is supplied, sub-cells with ``gasDen == 0`` are skipped
        (avoids refining voids where the velocity field is undefined).

        Parameters
        ----------
        vel_fn : callable
            ``vel_fn(x, y, z) -> (vx, vy, vz) [km/s]``.
        gasDen_fn : callable or None
            Optional density function.  Sub-cells with ``gasDen_fn(x,y,z) <= 0``
            do not contribute to the velocity gradient.
        threshold : float
            Gradient threshold ∈ [0, 1).  Default 0.1.
        level_max : int
            Maximum refinement level.
        level_min : int, optional (default 2)
            Minimum refinement level.  All cells are unconditionally split to
            this level before the gradient criterion is applied.
            Set to 0 to disable uniform pre-refinement.
        nprobe : int
            Sub-grid sampling density per axis.  Default 2 (corners).
        floor : float
            Velocity floor to avoid division by zero.
        """
        self.refine_uniform(level_min)
        self.refine(
            self._make_vel_criterion(vel_fn, gasDen_fn, threshold, floor, nprobe),
            level_max
        )

    def refine_by_physics(self, gasDen_fn, vel_fn=None,
                          dens_threshold=0.1, vel_threshold=0.1,
                          level_max=6, level_min=2, nprobe=2, floor=1e-30):
        """
        Refine cells satisfying either the density or velocity gradient criterion.

        This is a convenience wrapper that calls both
        :meth:`refine_by_density` and :meth:`refine_by_velocity` with a
        single pass over the leaf set (cells that satisfy *either* criterion
        are refined).

        Parameters
        ----------
        gasDen_fn : callable
            ``gasDen_fn(x, y, z) -> gasDen [cm^-3]``.
        vel_fn : callable or None
            ``vel_fn(x, y, z) -> (vx, vy, vz) [km/s]``.  If None, only the
            density criterion is applied.
        dens_threshold : float or None
            Density gradient threshold.  Pass ``None`` to skip.
        vel_threshold : float or None
            Velocity gradient threshold.  Pass ``None`` to skip (or if
            ``vel_fn`` is None).
        level_max : int
            Maximum refinement level.
        level_min : int, optional (default 2)
            Minimum refinement level.  All cells are unconditionally split to
            this level before the gradient criterion is applied.
            Set to 0 to disable uniform pre-refinement.
        nprobe : int
            Sub-grid sampling density per axis.  Default 2.
        floor : float
            Floor for division.
        """
        self.refine_uniform(level_min)
        criteria = []
        if dens_threshold is not None:
            criteria.append(self._make_dens_criterion(gasDen_fn, dens_threshold, floor, nprobe))
        if vel_fn is not None and vel_threshold is not None:
            criteria.append(self._make_vel_criterion(vel_fn, gasDen_fn, vel_threshold, floor, nprobe))
        if not criteria:
            return

        def combined(c):
            return any(crit(c) for crit in criteria)

        self.refine(combined, level_max)

    def _make_physics_criterion(self, gasDen_fn, vel_fn=None,
                                dens_threshold=0.1, vel_threshold=0.1,
                                nprobe=2, floor=1e-30):
        """Return the combined density/velocity refinement criterion."""
        criteria = []
        if dens_threshold is not None:
            criteria.append(self._make_dens_criterion(gasDen_fn, dens_threshold, floor, nprobe))
        if vel_fn is not None and vel_threshold is not None:
            criteria.append(self._make_vel_criterion(vel_fn, gasDen_fn, vel_threshold, floor, nprobe))
        if not criteria:
            return None
        return (lambda c: any(crit(c) for crit in criteria))

    def refine_sphere_by_physics(self, cx, cy, cz, radius, gasDen_fn, vel_fn=None,
                                 dens_threshold=0.1, vel_threshold=0.1,
                                 level_max=6, level_min=2, nprobe=2, floor=1e-30,
                                 refine_boundary=False):
        """
        Refine inside (or around) a sphere using physics-based gradient criteria.

        Parameters
        ----------
        cx, cy, cz : float
            Centre of the sphere.
        radius : float
            Sphere radius (same unit as boxlen).
        gasDen_fn : callable
            ``gasDen_fn(x, y, z) -> gasDen [cm^-3]``.
        vel_fn : callable or None
            ``vel_fn(x, y, z) -> (vx, vy, vz) [km/s]``.
        dens_threshold : float or None
            Density gradient threshold.  Pass ``None`` to skip.
        vel_threshold : float or None
            Velocity gradient threshold.  Pass ``None`` to skip.
        level_max : int
            Maximum refinement level.
        level_min : int, optional (default 2)
            Minimum refinement level.  All cells are unconditionally split to
            this level before the physics/geometry criterion is applied.
            Set to 0 to disable uniform pre-refinement.
        nprobe : int
            Sub-grid sampling density per axis.  Default 2.
        floor : float
            Floor for gradient denominator.
        refine_boundary : bool, optional (default False)
            If False (default), refine using only the density/velocity
            gradient criterion — no special treatment of the sphere surface.
            This is sufficient when the density field itself captures the
            boundary (e.g. density drops sharply at the sphere edge).
            If True, additionally force every cell that intersects the sphere
            surface to ``level_max``, ensuring sharp geometric resolution of
            the boundary regardless of the local gradient.
        """
        self.refine_uniform(level_min)
        criterion = self._make_physics_criterion(
            gasDen_fn, vel_fn, dens_threshold, vel_threshold, nprobe, floor
        )
        self.refine_sphere_adaptive(cx, cy, cz, radius, level_max, criterion,
                                    refine_boundary=refine_boundary)

    def refine_slab_by_physics(self, axis, lo, hi, gasDen_fn, vel_fn=None,
                               dens_threshold=0.1, vel_threshold=0.1,
                               level_max=6, level_min=2, nprobe=2, floor=1e-30,
                               refine_boundary=False):
        """
        Refine inside (or around) a slab using physics-based gradient criteria.

        Parameters
        ----------
        axis : {'x', 'y', 'z'}
            Normal axis of the slab.
        lo, hi : float
            Slab extent along ``axis``.
        gasDen_fn : callable
            ``gasDen_fn(x, y, z) -> gasDen [cm^-3]``.
        vel_fn : callable or None
            ``vel_fn(x, y, z) -> (vx, vy, vz) [km/s]``.
        dens_threshold, vel_threshold : float or None
            Gradient thresholds.
        level_max : int
            Maximum refinement level.
        level_min : int, optional (default 2)
            Minimum refinement level.  All cells are unconditionally split to
            this level before the physics/geometry criterion is applied.
            Set to 0 to disable uniform pre-refinement.
        nprobe : int
            Sub-grid sampling density per axis.  Default 2.
        floor : float
            Floor for gradient denominator.
        refine_boundary : bool, optional (default False)
            If False (default), use only the physics criterion.
            If True, also force cells that overlap the slab boundary to
            ``level_max``.
        """
        self.refine_uniform(level_min)
        criterion = self._make_physics_criterion(
            gasDen_fn, vel_fn, dens_threshold, vel_threshold, nprobe, floor
        )
        self.refine_slab_adaptive(axis, lo, hi, level_max, criterion,
                                  refine_boundary=refine_boundary)

    def refine_box_by_physics(self, xlo, xhi, ylo, yhi, zlo, zhi, gasDen_fn,
                              vel_fn=None, dens_threshold=0.1,
                              vel_threshold=0.1, level_max=6, level_min=2,
                              nprobe=2, floor=1e-30, refine_boundary=False):
        """
        Refine inside (or around) a box using physics-based gradient criteria.

        Parameters
        ----------
        xlo, xhi, ylo, yhi, zlo, zhi : float
            Box extents.
        gasDen_fn : callable
            ``gasDen_fn(x, y, z) -> gasDen [cm^-3]``.
        vel_fn : callable or None
            ``vel_fn(x, y, z) -> (vx, vy, vz) [km/s]``.
        dens_threshold, vel_threshold : float or None
            Gradient thresholds.
        level_max : int
            Maximum refinement level.
        level_min : int, optional (default 2)
            Minimum refinement level.  All cells are unconditionally split to
            this level before the physics/geometry criterion is applied.
            Set to 0 to disable uniform pre-refinement.
        nprobe : int
            Sub-grid sampling density per axis.  Default 2.
        floor : float
            Floor for gradient denominator.
        refine_boundary : bool, optional (default False)
            If False (default), use only the physics criterion.
            If True, also force cells that overlap the box boundary to
            ``level_max``.
        """
        self.refine_uniform(level_min)
        criterion = self._make_physics_criterion(
            gasDen_fn, vel_fn, dens_threshold, vel_threshold, nprobe, floor
        )
        self.refine_box_adaptive(xlo, xhi, ylo, yhi, zlo, zhi, level_max,
                                 criterion, refine_boundary=refine_boundary)

    # ------------------------------------------------------------------ #
    # Private gradient criterion builders
    # ------------------------------------------------------------------ #

    @staticmethod
    def _subcell_centers(c, np_):
        """
        Return arrays (xx, yy, zz) of ``np_^3`` sub-cell centres inside
        cell ``c`` at sub-grid density ``np_`` per axis.
        """
        dh = 2.0 * c.h / np_
        ox = c.cx - c.h + (np.arange(np_) + 0.5) * dh
        oy = c.cy - c.h + (np.arange(np_) + 0.5) * dh
        oz = c.cz - c.h + (np.arange(np_) + 0.5) * dh
        gx, gy, gz = np.meshgrid(ox, oy, oz, indexing='ij')
        return gx.ravel(), gy.ravel(), gz.ravel()

    @staticmethod
    def _make_dens_criterion(gasDen_fn, threshold, floor, nprobe):
        """
        Return a cell criterion function for density gradient refinement.

        Probes at np=2 first (8 sub-cell centres); if no refinement found and
        nprobe > 2, also probes at np=nprobe (nprobe^3 sub-cells).
        Accumulates running (rhomin, rhomax) across all probe levels,
        matching the RASCAS multi-scale approach.
        Numpy-vectorisable gasDen_fn is evaluated with a single array call
        per probe level (fast path); scalar functions use a loop fallback.
        """
        _probe = np.array([0.0, 1.0])
        try:
            _r = gasDen_fn(_probe, _probe, _probe)
            _vec = np.ndim(_r) >= 1
        except Exception:
            _vec = False

        def criterion(c):
            rhomin = rhomax = max(float(gasDen_fn(c.cx, c.cy, c.cz)), 0.0)

            for np_ in ([2] if nprobe <= 2 else [2, nprobe]):
                xx, yy, zz = AMRGrid._subcell_centers(c, np_)
                if _vec:
                    vals = np.maximum(gasDen_fn(xx, yy, zz), 0.0)
                    v_max = float(vals.max())
                    v_min = float(vals.min())
                    if v_max > rhomax: rhomax = v_max
                    if v_min < rhomin: rhomin = v_min
                else:
                    for x, y, z in zip(xx, yy, zz):
                        rho2 = max(float(gasDen_fn(x, y, z)), 0.0)
                        if rho2 > rhomax: rhomax = rho2
                        if rho2 < rhomin: rhomin = rho2
                delta = (rhomax - rhomin) / (rhomax + rhomin + floor)
                if delta >= threshold:
                    return True
            return False
        return criterion

    @staticmethod
    def _make_vel_criterion(vel_fn, gasDen_fn, threshold, floor, nprobe):
        """
        Return a cell criterion function for velocity gradient refinement.

        Same multi-scale probing as ``_make_dens_criterion``, but applied
        to |v|.  Sub-cells with gasDen == 0 are skipped when gasDen_fn is given.
        Numpy-vectorisable vel_fn is evaluated with a single array call per
        probe level (fast path).
        """
        _probe = np.array([0.0, 1.0])
        try:
            _rv = vel_fn(_probe, _probe, _probe)
            _vec_vel = np.ndim(_rv[0]) >= 1
        except Exception:
            _vec_vel = False

        if gasDen_fn is not None:
            try:
                _rd = gasDen_fn(_probe, _probe, _probe)
                _vec_den = np.ndim(_rd) >= 1
            except Exception:
                _vec_den = False
        else:
            _vec_den = False

        def speed_arr(xx, yy, zz):
            vx, vy, vz = vel_fn(xx, yy, zz)
            return np.sqrt(np.asarray(vx)**2 + np.asarray(vy)**2 + np.asarray(vz)**2)

        def criterion(c):
            if gasDen_fn is not None and float(gasDen_fn(c.cx, c.cy, c.cz)) <= 0.0:
                return False
            vx0, vy0, vz0 = vel_fn(c.cx, c.cy, c.cz)
            vmin = vmax = float(np.sqrt(float(vx0)**2 + float(vy0)**2 + float(vz0)**2))

            for np_ in ([2] if nprobe <= 2 else [2, nprobe]):
                xx, yy, zz = AMRGrid._subcell_centers(c, np_)
                if _vec_vel:
                    if gasDen_fn is not None and _vec_den:
                        mask = np.asarray(gasDen_fn(xx, yy, zz)) > 0.0
                        if not np.any(mask):
                            continue
                        spd = speed_arr(xx[mask], yy[mask], zz[mask])
                    else:
                        spd = speed_arr(xx, yy, zz)
                    v_max = float(spd.max())
                    v_min = float(spd.min())
                    if v_max > vmax: vmax = v_max
                    if v_min < vmin: vmin = v_min
                else:
                    for x, y, z in zip(xx, yy, zz):
                        if gasDen_fn is not None and float(gasDen_fn(x, y, z)) <= 0.0:
                            continue
                        vx2, vy2, vz2 = vel_fn(x, y, z)
                        v2 = float(np.sqrt(float(vx2)**2 + float(vy2)**2 + float(vz2)**2))
                        if v2 > vmax: vmax = v2
                        if v2 < vmin: vmin = v2
                delta = (vmax - vmin) / (vmax + vmin + floor)
                if delta >= threshold:
                    return True
            return False
        return criterion

    # ------------------------------------------------------------------ #
    # Physical property assignment
    # ------------------------------------------------------------------ #

    @staticmethod
    def _coords(leaflist):
        """Return (cx, cy, cz) numpy arrays for a list of leaves."""
        return (np.array([lf.cx for lf in leaflist]),
                np.array([lf.cy for lf in leaflist]),
                np.array([lf.cz for lf in leaflist]))

    @staticmethod
    def _broadcast(val, n):
        """Ensure val is a 1-D array of length n (broadcast scalar if needed)."""
        v = np.asarray(val, dtype=float)
        return np.full(n, float(v)) if v.ndim == 0 else v

    def set_density(self, fn):
        """
        Assign hydrogen number density to every leaf.

        Parameters
        ----------
        fn : callable or float
            ``fn(cx, cy, cz) -> gasDen [cm^-3]``.
            If a scalar is given, all leaves get that value.
            Numpy-vectorisable callables are evaluated in one call (fast path).
        """
        if callable(fn):
            leaflist = self.leaves()
            try:
                cx, cy, cz = self._coords(leaflist)
                vals = self._broadcast(fn(cx, cy, cz), len(leaflist))
            except Exception:
                for lf in leaflist:
                    lf.gasDen = float(fn(lf.cx, lf.cy, lf.cz))
                return
            for lf, v in zip(leaflist, vals):
                lf.gasDen = float(v)
        else:
            val = float(fn)
            for lf in self.leaves():
                lf.gasDen = val

    def set_temperature(self, fn):
        """
        Assign temperature to every leaf.

        Parameters
        ----------
        fn : callable or float
            ``fn(cx, cy, cz) -> T [K]``.
            Numpy-vectorisable callables are evaluated in one call (fast path).
        """
        if callable(fn):
            leaflist = self.leaves()
            try:
                cx, cy, cz = self._coords(leaflist)
                vals = self._broadcast(fn(cx, cy, cz), len(leaflist))
            except Exception:
                for lf in leaflist:
                    lf.T = float(fn(lf.cx, lf.cy, lf.cz))
                return
            for lf, v in zip(leaflist, vals):
                lf.T = float(v)
        else:
            val = float(fn)
            for lf in self.leaves():
                lf.T = val

    def set_velocity(self, fn):
        """
        Assign bulk velocity to every leaf.

        Parameters
        ----------
        fn : callable or tuple
            ``fn(cx, cy, cz) -> (vx, vy, vz) [km/s]``.
            A 3-tuple ``(vx, vy, vz)`` sets uniform velocity.
            Numpy-vectorisable callables are evaluated in one call (fast path).
        """
        if callable(fn):
            leaflist = self.leaves()
            try:
                cx, cy, cz = self._coords(leaflist)
                res = fn(cx, cy, cz)
                n = len(leaflist)
                vx_a = self._broadcast(res[0], n)
                vy_a = self._broadcast(res[1], n)
                vz_a = self._broadcast(res[2], n)
            except Exception:
                for lf in leaflist:
                    lf.vx, lf.vy, lf.vz = fn(lf.cx, lf.cy, lf.cz)
                return
            for lf, vx, vy, vz in zip(leaflist, vx_a, vy_a, vz_a):
                lf.vx = float(vx); lf.vy = float(vy); lf.vz = float(vz)
        else:
            vx, vy, vz = (float(v) for v in fn)
            for lf in self.leaves():
                lf.vx, lf.vy, lf.vz = vx, vy, vz

    def set_properties(self, fn):
        """
        Assign all physical properties at once.

        Parameters
        ----------
        fn : callable
            ``fn(cx, cy, cz) -> (gasDen, T, vx, vy, vz)``.
            Numpy-vectorisable callables are evaluated in one call (fast path).
        """
        leaflist = self.leaves()
        try:
            cx, cy, cz = self._coords(leaflist)
            res = fn(cx, cy, cz)
            n = len(leaflist)
            gd_a = self._broadcast(res[0], n)
            T_a  = self._broadcast(res[1], n)
            vx_a = self._broadcast(res[2], n)
            vy_a = self._broadcast(res[3], n)
            vz_a = self._broadcast(res[4], n)
        except Exception:
            for lf in leaflist:
                lf.gasDen, lf.T, lf.vx, lf.vy, lf.vz = fn(lf.cx, lf.cy, lf.cz)
            return
        for lf, gd, T, vx, vy, vz in zip(leaflist, gd_a, T_a, vx_a, vy_a, vz_a):
            lf.gasDen = float(gd); lf.T = float(T)
            lf.vx = float(vx); lf.vy = float(vy); lf.vz = float(vz)

    # ------------------------------------------------------------------ #
    # Tree traversal
    # ------------------------------------------------------------------ #

    def leaves(self):
        """Return a list of all leaf :class:`Cell` objects."""
        result = []
        stack = [self.root]
        while stack:
            cell = stack.pop()
            if cell.children is None:
                result.append(cell)
            else:
                stack.extend(cell.children)
        return result

    def level_counts(self):
        """Return a dict {level: count} of leaf cells per AMR level."""
        counts = {}
        for lf in self.leaves():
            counts[lf.level] = counts.get(lf.level, 0) + 1
        return dict(sorted(counts.items()))

    def info(self):
        """Return a summary string of the grid."""
        leaflist = self.leaves()
        counts = {}
        gasDen_list = []
        T_list = []
        for lf in leaflist:
            counts[lf.level] = counts.get(lf.level, 0) + 1
            gasDen_list.append(lf.gasDen)
            T_list.append(lf.T)
        lv = dict(sorted(counts.items()))

        lines = [
            f'AMRGrid  boxlen={self.boxlen}  nleaf={len(leaflist)}',
            f'  Level distribution: {lv}',
        ]
        gasDen_vals = np.array(gasDen_list)
        T_vals      = np.array(T_list)
        if gasDen_vals.max() > 0:
            lines.append(
                f'  gasDen [cm^-3]: min={gasDen_vals.min():.3e}  '
                f'max={gasDen_vals.max():.3e}  mean={gasDen_vals.mean():.3e}'
            )
        lines.append(
            f'  T  [K]:     min={T_vals.min():.3e}  '
            f'max={T_vals.max():.3e}'
        )
        return '\n'.join(lines)

    # ------------------------------------------------------------------ #
    # I/O
    # ------------------------------------------------------------------ #

    @staticmethod
    def _is_fits_filename(filename):
        name = str(filename).lower()
        return name.endswith('.fits') or name.endswith('.fits.gz')

    @staticmethod
    def _require_fits():
        if fits is None:
            raise ImportError(
                "FITS I/O requires astropy. Install astropy to read/write "
                ".fits or .fits.gz AMR files."
            )

    def _leaf_table_array(self):
        leaflist = self.leaves()
        data = np.empty(len(leaflist), dtype=[
            ('x', 'f8'), ('y', 'f8'), ('z', 'f8'), ('level', 'i4'),
            ('gasDen', 'f8'), ('T', 'f8'), ('vx', 'f8'), ('vy', 'f8'), ('vz', 'f8'),
        ])
        # Column-wise assignment is much faster than per-row structured-array writes
        data['x']      = [lf.cx     for lf in leaflist]
        data['y']      = [lf.cy     for lf in leaflist]
        data['z']      = [lf.cz     for lf in leaflist]
        data['level']  = [lf.level  for lf in leaflist]
        data['gasDen'] = [lf.gasDen for lf in leaflist]
        data['T']      = [lf.T      for lf in leaflist]
        data['vx']     = [lf.vx     for lf in leaflist]
        data['vy']     = [lf.vy     for lf in leaflist]
        data['vz']     = [lf.vz     for lf in leaflist]
        return data

    def _write_fits(self, filename):
        self._require_fits()
        data = self._leaf_table_array()
        primary = fits.PrimaryHDU()
        table = fits.BinTableHDU(data=data, name='AMRGRID')
        table.header['BOXLEN'] = (self.boxlen, 'Simulation box length')
        table.header['ORIGINX'] = (self.origin[0], 'Box origin x')
        table.header['ORIGINY'] = (self.origin[1], 'Box origin y')
        table.header['ORIGINZ'] = (self.origin[2], 'Box origin z')
        table.header['NLEAF'] = (len(data), 'Number of AMR leaf cells')
        fits.HDUList([primary, table]).writeto(filename, overwrite=True)
        print(f'Written {len(data)} leaf cells to {filename}')

    @classmethod
    def _read_fits(cls, filename):
        cls._require_fits()
        with fits.open(filename) as hdul:
            table_hdu = None
            for hdu in hdul:
                if isinstance(hdu, fits.BinTableHDU):
                    table_hdu = hdu
                    break
            if table_hdu is None:
                raise ValueError(f'No FITS binary table found in {filename}')

            header = table_hdu.header
            boxlen = float(header['BOXLEN'])
            origin = (
                float(header.get('ORIGINX', 0.0)),
                float(header.get('ORIGINY', 0.0)),
                float(header.get('ORIGINZ', 0.0)),
            )
            table = table_hdu.data
            gasDen_key = 'gasDen' if 'gasDen' in table.names else 'nH'
            data = np.column_stack([
                np.asarray(table['x'], dtype=float),
                np.asarray(table['y'], dtype=float),
                np.asarray(table['z'], dtype=float),
                np.asarray(table['level'], dtype=int),
                np.asarray(table[gasDen_key], dtype=float),
                np.asarray(table['T'], dtype=float),
                np.asarray(table['vx'], dtype=float),
                np.asarray(table['vy'], dtype=float),
                np.asarray(table['vz'], dtype=float),
            ])

        return cls._build_from_rows(boxlen, data, origin=origin)

    def _insert_cell(self, cx, cy, cz, target_level,
                     gasDen=0.0, T=1e4, vx=0.0, vy=0.0, vz=0.0):
        """Navigate from root to the leaf at (cx,cy,cz)/target_level, splitting as needed.

        O(target_level) per call — avoids the O(tree_size) full traversal of refine().
        Child index convention matches Cell.split(): iz*4 + iy*2 + ix.
        """
        cell = self.root
        for _ in range(target_level):
            if cell.is_leaf:
                cell.split()
            ix = 1 if cx >= cell.cx else 0
            iy = 1 if cy >= cell.cy else 0
            iz = 1 if cz >= cell.cz else 0
            cell = cell.children[iz * 4 + iy * 2 + ix]
        cell.gasDen = gasDen
        cell.T  = T
        cell.vx = vx
        cell.vy = vy
        cell.vz = vz

    @classmethod
    def _build_from_rows(cls, boxlen, data, origin=(0.0, 0.0, 0.0)):
        grid = cls(boxlen, origin=origin)
        for row in data:
            grid._insert_cell(
                float(row[0]), float(row[1]), float(row[2]), int(row[3]),
                float(row[4]), float(row[5]), float(row[6]), float(row[7]), float(row[8]),
            )
        return grid

    def write(self, filename):
        """
        Write the AMR grid to a LaRT generic file.

        Text output is used by default. If ``filename`` ends with ``.fits``
        or ``.fits.gz``, the grid is written as a FITS binary table.

        Parameters
        ----------
        filename : str
            Output file path.
        """
        if self._is_fits_filename(filename):
            self._write_fits(filename)
            return

        leaflist = self.leaves()
        nleaf    = len(leaflist)
        mat = np.column_stack([
            [lf.cx     for lf in leaflist],
            [lf.cy     for lf in leaflist],
            [lf.cz     for lf in leaflist],
            [lf.level  for lf in leaflist],
            [lf.gasDen for lf in leaflist],
            [lf.T      for lf in leaflist],
            [lf.vx     for lf in leaflist],
            [lf.vy     for lf in leaflist],
            [lf.vz     for lf in leaflist],
        ])
        with open(filename, 'w') as f:
            f.write(f'{nleaf}  {self.boxlen:.6f}\n')
            np.savetxt(f, mat,
                       fmt='%.6f  %.6f  %.6f  %d  %.6e  %.4e  %.4f  %.4f  %.4f')
        print(f'Written {nleaf} leaf cells to {filename}')

    @classmethod
    def read(cls, filename):
        """
        Read a LaRT generic AMR file back into an :class:`AMRGrid`.

        Returns a grid whose octree is rebuilt by inserting each leaf cell.
        Physical data (gasDen, T, velocity) is stored on the leaf cells.

        Parameters
        ----------
        filename : str
            Path to the AMR data file.
        """
        if cls._is_fits_filename(filename):
            return cls._read_fits(filename)

        data = np.loadtxt(filename, skiprows=1)
        with open(filename) as f:
            first = f.readline().split()
        nleaf, boxlen = int(first[0]), float(first[1])
        if data.ndim == 1:
            data = data[np.newaxis, :]

        if len(data) != nleaf:
            raise ValueError(
                f'Header nleaf={nleaf} does not match row count {len(data)} in {filename}'
            )

        return cls._build_from_rows(boxlen, data)

    # ------------------------------------------------------------------ #
    # Visualisation helpers
    # ------------------------------------------------------------------ #

    def _collect_slice_leaves(self, cell, na, value, result):
        """Traverse the tree, pruning subtrees entirely outside the slice plane."""
        if abs(getattr(cell, na) - value) > cell.h:
            return
        if cell.is_leaf:
            result.append(cell)
        else:
            for child in cell.children:
                self._collect_slice_leaves(child, na, value, result)

    def slice_plot(self, axis='z', value=None, quantity='gasDen',
                   ax=None, log=False, cmap='viridis',
                   background_color='white',
                   show_leaf_boundaries=False,
                   boundary_color='k', boundary_lw=0.3,
                   boundary_alpha=0.7,
                   show_leaf_centers=False,
                   center_color='k', center_marker='o',
                   center_size=8, center_alpha=0.9, **kwargs):
        """
        Quick 2-D slice plot of a physical quantity through the grid.

        Uses tree pruning (skips subtrees outside the slice plane entirely)
        and numpy-vectorised polygon construction for speed at high refinement
        levels (e.g. ``level_max >= 6``).

        Parameters
        ----------
        axis : {'x', 'y', 'z'}
            Normal axis of the slice plane.
        value : float, optional
            Position of the slice along ``axis``.  Defaults to box centre.
        quantity : {'gasDen', 'T', 'vx', 'vy', 'vz'}
            Quantity to plot.
        ax : matplotlib.axes.Axes, optional
        log : bool
            Plot log10 of the quantity if True.
        cmap : str
        background_color : str or None, optional (default 'white')
            Colour of the axes background — the region outside (or between)
            the cell polygons.  Any matplotlib colour string is accepted
            (e.g. ``'white'``, ``'black'``, ``'lightgray'``, ``'#eeeeee'``).
            Pass ``None`` to leave the axes background unchanged (inherits
            from the current matplotlib style).
        show_leaf_boundaries : bool
            If True, overlay boundaries of the leaf cells that intersect
            the slice plane.
        boundary_color : str
            Edge colour for the leaf-cell boundary overlay.
        boundary_lw : float
            Line width for the leaf-cell boundary overlay.
        boundary_alpha : float
            Transparency for the leaf-cell boundary overlay.
        show_leaf_centers : bool
            If True, overlay markers at the centres of leaf cells that
            intersect the slice plane.
        center_color : str
            Marker colour for the leaf-cell centre overlay.
        center_marker : str
            Matplotlib marker style for the leaf-cell centre overlay.
        center_size : float
            Marker size for the leaf-cell centre overlay.
        center_alpha : float
            Transparency for the leaf-cell centre overlay.
        **kwargs
            Passed to ``matplotlib.collections.PolyCollection`` for the filled
            cells.

        Returns
        -------
        im : matplotlib.collections.PolyCollection
        """
        import matplotlib.pyplot as plt
        from matplotlib.collections import PolyCollection
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes

        if ax is None:
            _, ax = plt.subplots()

        if background_color is not None:
            ax.set_facecolor(background_color)

        axis = axis.lower()
        if value is None:
            value = self.origin[{'x': 0, 'y': 1, 'z': 2}[axis]] + self.boxlen * 0.5

        ax_map = {'x': ('y', 'z', 'cy', 'cz', 'cx'),
                  'y': ('x', 'z', 'cx', 'cz', 'cy'),
                  'z': ('x', 'y', 'cx', 'cy', 'cz')}
        ha, va, ca, da, na = ax_map[axis]
        quantity_key = 'gasDen' if quantity == 'nH' else quantity

        # Collect only leaves that intersect the slice plane (tree pruning).
        # This visits O(N_slice * depth) nodes instead of all N_total leaves.
        slice_leaves = []
        self._collect_slice_leaves(self.root, na, value, slice_leaves)

        n = len(slice_leaves)
        if n == 0:
            return None

        # Build numpy arrays in a single pass over the (small) slice leaf list.
        h_arr  = np.empty(n, dtype=float)
        cx_arr = np.empty(n, dtype=float)
        cy_arr = np.empty(n, dtype=float)
        v_arr  = np.empty(n, dtype=float)
        for i, lf in enumerate(slice_leaves):
            h_arr[i]  = lf.h
            cx_arr[i] = getattr(lf, ca)
            cy_arr[i] = getattr(lf, da)
            v_arr[i]  = getattr(lf, quantity_key)

        xmin = cx_arr - h_arr;  xmax = cx_arr + h_arr
        ymin = cy_arr - h_arr;  ymax = cy_arr + h_arr

        # Build (N, 4, 2) polygon array without any Python-level list nesting.
        polys = np.empty((n, 4, 2), dtype=float)
        polys[:, 0, 0] = xmin;  polys[:, 0, 1] = ymin
        polys[:, 1, 0] = xmax;  polys[:, 1, 1] = ymin
        polys[:, 2, 0] = xmax;  polys[:, 2, 1] = ymax
        polys[:, 3, 0] = xmin;  polys[:, 3, 1] = ymax

        values = np.log10(np.maximum(v_arr, 1e-100)) if log else v_arr

        pc = PolyCollection(
            polys,
            cmap=cmap,
            edgecolors='none',
            closed=True,
            **kwargs,
        )
        pc.set_array(values)
        ax.add_collection(pc)

        if show_leaf_boundaries:
            pc_boundary = PolyCollection(
                polys,
                facecolors='none',
                edgecolors=boundary_color,
                linewidths=boundary_lw,
                alpha=boundary_alpha,
                closed=True,
            )
            ax.add_collection(pc_boundary)

        if show_leaf_centers:
            ax.scatter(
                cx_arr, cy_arr,
                c=center_color, marker=center_marker, s=center_size,
                alpha=center_alpha
            )

        ax.set_xlim(self.origin[0], self.origin[0] + self.boxlen)
        ax.set_ylim(self.origin[1], self.origin[1] + self.boxlen)
        ax.set_aspect('equal')
        ax.set_xlabel(ha)
        ax.set_ylabel(va)
        label = f'log10({quantity_key})' if log else quantity_key

        cax = inset_axes(
            ax,
            width="4%",
            height="100%",
            loc='lower left',
            bbox_to_anchor=(1.02, 0.0, 1.0, 1.0),
            bbox_transform=ax.transAxes,
            borderpad=0,
        )
        ax.figure.colorbar(pc, cax=cax, label=label)
        return pc
