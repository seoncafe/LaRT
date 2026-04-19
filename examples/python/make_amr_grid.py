"""
make_amr_grid.py
================
Library for building generic-format AMR data files for LaRT v2.00.

Output file format (amr_type = 'generic')
------------------------------------------
    nleaf  boxlen
    x  y  z  level  nH[cm^-3]  T[K]  vx[km/s]  vy[km/s]  vz[km/s]
    ...

  - x, y, z   : leaf-cell centre coordinates in the same unit as boxlen
                 (origin at box corner, range [0, boxlen])
  - level      : AMR refinement level (root = 0; minimum leaf level = 1)
  - nH         : total hydrogen number density [cm^-3]
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
    def nH_fn(x, y, z):
        r = np.sqrt((x-50)**2 + (y-50)**2 + (z-50)**2)
        return 1e-3 * np.exp(-r**2 / (2*20.0**2))

    grid.set_density(nH_fn)
    grid.set_temperature(lambda x, y, z: 1e4)
    # velocity defaults to zero

    grid.write('my_grid.dat')
    print(grid.info())
"""

import numpy as np


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
    nH, T, vx, vy, vz : float
        Physical data (only meaningful for leaves).
    """

    __slots__ = ('cx', 'cy', 'cz', 'h', 'level', 'children',
                 'nH', 'T', 'vx', 'vy', 'vz')

    def __init__(self, cx, cy, cz, h, level):
        self.cx    = cx
        self.cy    = cy
        self.cz    = cz
        self.h     = h
        self.level = level
        self.children = None
        # Physical defaults
        self.nH = 0.0
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
                               criterion_inside=None):
        """
        Geometry-aware refinement for a sphere.

        Boundary cells (cells that intersect the sphere but are not fully
        contained in it) are always refined until ``level_max``.  Cells
        fully inside the sphere are refined only when ``criterion_inside``
        evaluates to True.
        """
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
                             criterion_inside=None):
        """
        Geometry-aware refinement for a slab.

        The slab spans ``[lo, hi]`` along ``axis`` ('x', 'y', or 'z').
        Boundary cells are always refined until ``level_max``; cells fully
        inside the slab are refined only when ``criterion_inside`` is True.
        """
        axis = axis.lower()
        min_attr = {'x': 'xmin', 'y': 'ymin', 'z': 'zmin'}[axis]
        max_attr = {'x': 'xmax', 'y': 'ymax', 'z': 'zmax'}[axis]

        def inside(c):
            return lo <= getattr(c, min_attr) and getattr(c, max_attr) <= hi

        def intersects(c):
            return getattr(c, max_attr) >= lo and getattr(c, min_attr) <= hi

        self.refine_geometry(inside, intersects, level_max, criterion_inside)

    def refine_box_adaptive(self, xlo, xhi, ylo, yhi, zlo, zhi, level_max,
                            criterion_inside=None):
        """
        Geometry-aware refinement for an axis-aligned rectangular box.

        Boundary cells are always refined until ``level_max``; cells fully
        inside the box are refined only when ``criterion_inside`` is True.
        """
        def inside(c):
            return (xlo <= c.xmin and c.xmax <= xhi and
                    ylo <= c.ymin and c.ymax <= yhi and
                    zlo <= c.zmin and c.zmax <= zhi)

        def intersects(c):
            return (c.xmax >= xlo and c.xmin <= xhi and
                    c.ymax >= ylo and c.ymin <= yhi and
                    c.zmax >= zlo and c.zmin <= zhi)

        self.refine_geometry(inside, intersects, level_max, criterion_inside)

    def refine_by_density(self, nH_fn, threshold=0.1, level_max=6,
                          nprobe=2, floor=1e-30):
        """
        Refine cells where the density gradient exceeds ``threshold``.

        Follows the RASCAS (module_mesh.f90) criterion:

            delta = (nH_max - nH_min) / (nH_max + nH_min + floor)

        The density is sampled at a uniform ``nprobe × nprobe × nprobe``
        sub-grid inside the cell plus the cell centre.  ``nprobe=2``
        samples the 8 cell corners (fast); ``nprobe=4`` (64 points) gives
        a more thorough probe.  With ``nprobe > 2`` the function also tries
        ``nprobe=2`` first and returns early if the gradient already exceeds
        the threshold (same early-exit logic as the Fortran reference).

        Parameters
        ----------
        nH_fn : callable
            ``nH_fn(x, y, z) -> nH [cm^-3]`` in the same coordinate units
            as the grid (e.g. kpc).
        threshold : float
            Gradient threshold ∈ [0, 1).  0.1 means 10 % relative variation.
        level_max : int
            Maximum refinement level.
        nprobe : int
            Sub-grid sampling density per axis (≥ 2).  Default 2 (corners).
        floor : float
            Density floor to avoid division by zero.
        """
        self.refine(
            self._make_dens_criterion(nH_fn, threshold, floor, nprobe),
            level_max
        )

    def refine_by_velocity(self, vel_fn, nH_fn=None, threshold=0.1,
                           level_max=6, nprobe=2, floor=1e-30):
        """
        Refine cells where the velocity magnitude gradient exceeds ``threshold``.

        Follows the RASCAS criterion applied to |v|:

            delta = (|v|_max - |v|_min) / (|v|_max + |v|_min + floor)

        Sampling is done at the same sub-grid as :meth:`refine_by_density`.
        If ``nH_fn`` is supplied, sub-cells with ``nH == 0`` are skipped
        (avoids refining voids where the velocity field is undefined).

        Parameters
        ----------
        vel_fn : callable
            ``vel_fn(x, y, z) -> (vx, vy, vz) [km/s]``.
        nH_fn : callable or None
            Optional density function.  Sub-cells with ``nH_fn(x,y,z) <= 0``
            do not contribute to the velocity gradient.
        threshold : float
            Gradient threshold ∈ [0, 1).  Default 0.1.
        level_max : int
            Maximum refinement level.
        nprobe : int
            Sub-grid sampling density per axis.  Default 2 (corners).
        floor : float
            Velocity floor to avoid division by zero.
        """
        self.refine(
            self._make_vel_criterion(vel_fn, nH_fn, threshold, floor, nprobe),
            level_max
        )

    def refine_by_physics(self, nH_fn, vel_fn=None,
                          dens_threshold=0.1, vel_threshold=0.1,
                          level_max=6, nprobe=2, floor=1e-30):
        """
        Refine cells satisfying either the density or velocity gradient criterion.

        This is a convenience wrapper that calls both
        :meth:`refine_by_density` and :meth:`refine_by_velocity` with a
        single pass over the leaf set (cells that satisfy *either* criterion
        are refined).

        Parameters
        ----------
        nH_fn : callable
            ``nH_fn(x, y, z) -> nH [cm^-3]``.
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
        nprobe : int
            Sub-grid sampling density per axis.  Default 2.
        floor : float
            Floor for division.
        """
        criteria = []
        if dens_threshold is not None:
            criteria.append(self._make_dens_criterion(nH_fn, dens_threshold, floor, nprobe))
        if vel_fn is not None and vel_threshold is not None:
            criteria.append(self._make_vel_criterion(vel_fn, nH_fn, vel_threshold, floor, nprobe))
        if not criteria:
            return

        def combined(c):
            return any(crit(c) for crit in criteria)

        self.refine(combined, level_max)

    def _make_physics_criterion(self, nH_fn, vel_fn=None,
                                dens_threshold=0.1, vel_threshold=0.1,
                                nprobe=2, floor=1e-30):
        """Return the combined density/velocity refinement criterion."""
        criteria = []
        if dens_threshold is not None:
            criteria.append(self._make_dens_criterion(nH_fn, dens_threshold, floor, nprobe))
        if vel_fn is not None and vel_threshold is not None:
            criteria.append(self._make_vel_criterion(vel_fn, nH_fn, vel_threshold, floor, nprobe))
        if not criteria:
            return None
        return (lambda c: any(crit(c) for crit in criteria))

    def refine_sphere_by_physics(self, cx, cy, cz, radius, nH_fn, vel_fn=None,
                                 dens_threshold=0.1, vel_threshold=0.1,
                                 level_max=6, nprobe=2, floor=1e-30):
        """
        Sphere refinement with boundary cells forced to ``level_max``.

        A leaf cell fully inside the sphere is refined only when the density
        and/or velocity gradient criterion is satisfied.  Any cell that
        intersects the sphere but is not fully contained in it is refined all
        the way to ``level_max`` so that the geometric boundary is resolved.
        """
        criterion = self._make_physics_criterion(
            nH_fn, vel_fn, dens_threshold, vel_threshold, nprobe, floor
        )
        self.refine_sphere_adaptive(cx, cy, cz, radius, level_max, criterion)

    def refine_slab_by_physics(self, axis, lo, hi, nH_fn, vel_fn=None,
                               dens_threshold=0.1, vel_threshold=0.1,
                               level_max=6, nprobe=2, floor=1e-30):
        """
        Slab refinement with boundary cells forced to ``level_max``.

        A leaf cell fully inside the slab is refined only when the density
        and/or velocity gradient criterion is satisfied.  Any cell that
        overlaps the slab boundary is refined all the way to ``level_max``.
        """
        criterion = self._make_physics_criterion(
            nH_fn, vel_fn, dens_threshold, vel_threshold, nprobe, floor
        )
        self.refine_slab_adaptive(axis, lo, hi, level_max, criterion)

    def refine_box_by_physics(self, xlo, xhi, ylo, yhi, zlo, zhi, nH_fn,
                              vel_fn=None, dens_threshold=0.1,
                              vel_threshold=0.1, level_max=6, nprobe=2,
                              floor=1e-30):
        """
        Box refinement with boundary cells forced to ``level_max``.

        A leaf cell fully inside the box is refined only when the density
        and/or velocity gradient criterion is satisfied.  Any cell that
        overlaps the box boundary is refined all the way to ``level_max``.
        """
        criterion = self._make_physics_criterion(
            nH_fn, vel_fn, dens_threshold, vel_threshold, nprobe, floor
        )
        self.refine_box_adaptive(xlo, xhi, ylo, yhi, zlo, zhi, level_max, criterion)

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
    def _make_dens_criterion(nH_fn, threshold, floor, nprobe):
        """
        Return a cell criterion function for density gradient refinement.

        Probes at np=2 first (8 corners); if no refinement found and
        nprobe > 2, also probes at np=nprobe (nprobe^3 sub-cells).
        Accumulates running (rhomin, rhomax) across all probe levels,
        matching the RASCAS multi-scale approach.
        """
        def criterion(c):
            rhomin = rhomax = max(nH_fn(c.cx, c.cy, c.cz), 0.0)

            for np_ in ([2] if nprobe <= 2 else [2, nprobe]):
                xx, yy, zz = AMRGrid._subcell_centers(c, np_)
                for x, y, z in zip(xx, yy, zz):
                    rho2 = max(nH_fn(x, y, z), 0.0)
                    if rho2 > rhomax: rhomax = rho2
                    if rho2 < rhomin: rhomin = rho2
                delta = (rhomax - rhomin) / (rhomax + rhomin + floor)
                if delta >= threshold:
                    return True
            return False
        return criterion

    @staticmethod
    def _make_vel_criterion(vel_fn, nH_fn, threshold, floor, nprobe):
        """
        Return a cell criterion function for velocity gradient refinement.

        Same multi-scale probing as ``_make_dens_criterion``, but applied
        to |v|.  Sub-cells with nH == 0 are skipped when nH_fn is given.
        """
        def speed(x, y, z):
            vx, vy, vz = vel_fn(x, y, z)
            return np.sqrt(vx**2 + vy**2 + vz**2)

        def has_gas(x, y, z):
            return (nH_fn is None) or (nH_fn(x, y, z) > 0.0)

        def criterion(c):
            if not has_gas(c.cx, c.cy, c.cz):
                return False
            vmin = vmax = speed(c.cx, c.cy, c.cz)

            for np_ in ([2] if nprobe <= 2 else [2, nprobe]):
                xx, yy, zz = AMRGrid._subcell_centers(c, np_)
                for x, y, z in zip(xx, yy, zz):
                    if not has_gas(x, y, z):
                        continue
                    v2 = speed(x, y, z)
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

    def set_density(self, fn):
        """
        Assign hydrogen number density to every leaf.

        Parameters
        ----------
        fn : callable or float
            ``fn(cx, cy, cz) -> nH [cm^-3]``.
            If a scalar is given, all leaves get that value.
        """
        if callable(fn):
            for lf in self.leaves():
                lf.nH = float(fn(lf.cx, lf.cy, lf.cz))
        else:
            val = float(fn)
            for lf in self.leaves():
                lf.nH = val

    def set_temperature(self, fn):
        """
        Assign temperature to every leaf.

        Parameters
        ----------
        fn : callable or float
            ``fn(cx, cy, cz) -> T [K]``.
        """
        if callable(fn):
            for lf in self.leaves():
                lf.T = float(fn(lf.cx, lf.cy, lf.cz))
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
        """
        if callable(fn):
            for lf in self.leaves():
                lf.vx, lf.vy, lf.vz = fn(lf.cx, lf.cy, lf.cz)
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
            ``fn(cx, cy, cz) -> (nH, T, vx, vy, vz)``.
        """
        for lf in self.leaves():
            lf.nH, lf.T, lf.vx, lf.vy, lf.vz = fn(lf.cx, lf.cy, lf.cz)

    # ------------------------------------------------------------------ #
    # Tree traversal
    # ------------------------------------------------------------------ #

    def leaves(self):
        """Return a list of all leaf :class:`Cell` objects."""
        result = []
        self._collect_leaves(self.root, result)
        return result

    def _collect_leaves(self, cell, result):
        if cell.is_leaf:
            result.append(cell)
        else:
            for child in cell.children:
                self._collect_leaves(child, result)

    def level_counts(self):
        """Return a dict {level: count} of leaf cells per AMR level."""
        counts = {}
        for lf in self.leaves():
            counts[lf.level] = counts.get(lf.level, 0) + 1
        return dict(sorted(counts.items()))

    def info(self):
        """Return a summary string of the grid."""
        lv = self.level_counts()
        lines = [
            f'AMRGrid  boxlen={self.boxlen}  nleaf={sum(lv.values())}',
            f'  Level distribution: {lv}',
        ]
        nH_vals = np.array([lf.nH for lf in self.leaves()])
        T_vals  = np.array([lf.T  for lf in self.leaves()])
        if nH_vals.max() > 0:
            lines.append(
                f'  nH [cm^-3]: min={nH_vals.min():.3e}  '
                f'max={nH_vals.max():.3e}  mean={nH_vals.mean():.3e}'
            )
        lines.append(
            f'  T  [K]:     min={T_vals.min():.3e}  '
            f'max={T_vals.max():.3e}'
        )
        return '\n'.join(lines)

    # ------------------------------------------------------------------ #
    # I/O
    # ------------------------------------------------------------------ #

    def write(self, filename):
        """
        Write the AMR grid to a LaRT generic text file.

        Parameters
        ----------
        filename : str
            Output file path.
        """
        leaflist = self.leaves()
        nleaf    = len(leaflist)
        with open(filename, 'w') as f:
            f.write(f'{nleaf}  {self.boxlen:.6f}\n')
            for lf in leaflist:
                f.write(
                    f'{lf.cx:.6f}  {lf.cy:.6f}  {lf.cz:.6f}  '
                    f'{lf.level}  '
                    f'{lf.nH:.6e}  {lf.T:.4e}  '
                    f'{lf.vx:.4f}  {lf.vy:.4f}  {lf.vz:.4f}\n'
                )
        print(f'Written {nleaf} leaf cells to {filename}')

    @classmethod
    def read(cls, filename):
        """
        Read a LaRT generic AMR file back into an :class:`AMRGrid`.

        Returns a grid whose octree is rebuilt by inserting each leaf cell.
        Physical data (nH, T, velocity) is stored on the leaf cells.

        Parameters
        ----------
        filename : str
            Path to the AMR data file.
        """
        data = np.loadtxt(filename, skiprows=1)
        with open(filename) as f:
            first = f.readline().split()
        nleaf, boxlen = int(first[0]), float(first[1])

        grid = cls(boxlen)

        # Insert leaves by descending the tree to the correct level
        for row in data:
            cx, cy, cz = row[0], row[1], row[2]
            level       = int(row[3])
            nH, T       = row[4], row[5]
            vx, vy, vz  = row[6], row[7], row[8]
            # Refine until the cell at (cx,cy,cz) exists at `level`
            grid.refine(
                lambda c, _cx=cx, _cy=cy, _cz=cz, _lv=level:
                    c.level < _lv and
                    c.xmin <= _cx <= c.xmax and
                    c.ymin <= _cy <= c.ymax and
                    c.zmin <= _cz <= c.zmax,
                level_max=level
            )
        # Now assign physical data
        leaflist = grid.leaves()
        for lf, row in zip(
            sorted(leaflist, key=lambda c: (round(c.cx,6), round(c.cy,6), round(c.cz,6))),
            sorted(data,     key=lambda r: (round(r[0],6), round(r[1],6), round(r[2],6)))
        ):
            lf.nH = row[4]; lf.T = row[5]
            lf.vx = row[6]; lf.vy = row[7]; lf.vz = row[8]
        return grid

    # ------------------------------------------------------------------ #
    # Visualisation helpers
    # ------------------------------------------------------------------ #

    def slice_plot(self, axis='z', value=None, quantity='nH',
                   ax=None, log=True, cmap='viridis',
                   show_leaf_boundaries=False,
                   boundary_color='k', boundary_lw=0.3,
                   boundary_alpha=0.7,
                   show_leaf_centers=False,
                   center_color='k', center_marker='o',
                   center_size=8, center_alpha=0.9, **kwargs):
        """
        Quick 2-D slice plot of a physical quantity through the grid.

        Parameters
        ----------
        axis : {'x', 'y', 'z'}
            Normal axis of the slice plane.
        value : float, optional
            Position of the slice along ``axis``.  Defaults to box centre.
        quantity : {'nH', 'T', 'vx', 'vy', 'vz'}
            Quantity to plot.
        ax : matplotlib.axes.Axes, optional
        log : bool
            Plot log10 of the quantity if True.
        cmap : str
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
            Passed to ``matplotlib.patches.Rectangle`` for the filled cells.

        Returns
        -------
        im : matplotlib.collections.PatchCollection
        """
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib.collections import PatchCollection
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes

        if ax is None:
            _, ax = plt.subplots()

        if value is None:
            value = self.origin[{'x': 0, 'y': 1, 'z': 2}[axis]] + self.boxlen * 0.5

        ax_map = {'x': ('y', 'z', 'cy', 'cz', 'cx'),
                  'y': ('x', 'z', 'cx', 'cz', 'cy'),
                  'z': ('x', 'y', 'cx', 'cy', 'cz')}
        ha, va, ca, da, na = ax_map[axis.lower()]

        patches, values = [], []
        boundary_patches = []
        center_x, center_y = [], []
        for lf in self.leaves():
            c_normal = getattr(lf, na)
            if abs(c_normal - value) <= lf.h:
                cx_ = getattr(lf, ca)
                cy_ = getattr(lf, da)
                rect = mpatches.Rectangle(
                    (cx_ - lf.h, cy_ - lf.h), 2*lf.h, 2*lf.h, **kwargs)
                patches.append(rect)
                if show_leaf_boundaries:
                    boundary_patches.append(
                        mpatches.Rectangle(
                            (cx_ - lf.h, cy_ - lf.h), 2*lf.h, 2*lf.h
                        )
                    )
                if show_leaf_centers:
                    center_x.append(cx_)
                    center_y.append(cy_)
                v = getattr(lf, quantity)
                values.append(np.log10(max(v, 1e-100)) if log else v)

        pc = PatchCollection(patches, cmap=cmap, linewidths=0)
        pc.set_array(np.array(values))
        ax.add_collection(pc)

        if show_leaf_boundaries and boundary_patches:
            pc_boundary = PatchCollection(
                boundary_patches,
                facecolor='none',
                edgecolor=boundary_color,
                linewidths=boundary_lw,
                alpha=boundary_alpha,
            )
            ax.add_collection(pc_boundary)

        if show_leaf_centers and center_x:
            ax.scatter(
                center_x, center_y,
                c=center_color, marker=center_marker, s=center_size,
                alpha=center_alpha
            )

        ax.set_xlim(self.origin[0], self.origin[0] + self.boxlen)
        ax.set_ylim(self.origin[1], self.origin[1] + self.boxlen)
        ax.set_aspect('equal')
        ax.set_xlabel(ha)
        ax.set_ylabel(va)
        label = f'log10({quantity})' if log else quantity

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
