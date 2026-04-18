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
                   ax=None, log=True, cmap='viridis', **kwargs):
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
        **kwargs
            Passed to ``matplotlib.patches.Rectangle``.

        Returns
        -------
        im : matplotlib.collections.PatchCollection
        """
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib.collections import PatchCollection

        if ax is None:
            _, ax = plt.subplots()

        if value is None:
            value = self.origin[{'x': 0, 'y': 1, 'z': 2}[axis]] + self.boxlen * 0.5

        ax_map = {'x': ('y', 'z', 'cy', 'cz', 'cx'),
                  'y': ('x', 'z', 'cx', 'cz', 'cy'),
                  'z': ('x', 'y', 'cx', 'cy', 'cz')}
        ha, va, ca, da, na = ax_map[axis.lower()]

        patches, values = [], []
        for lf in self.leaves():
            c_normal = getattr(lf, na)
            if abs(c_normal - value) <= lf.h:
                cx_ = getattr(lf, ca)
                cy_ = getattr(lf, da)
                rect = mpatches.Rectangle(
                    (cx_ - lf.h, cy_ - lf.h), 2*lf.h, 2*lf.h, **kwargs)
                patches.append(rect)
                v = getattr(lf, quantity)
                values.append(np.log10(max(v, 1e-100)) if log else v)

        pc = PatchCollection(patches, cmap=cmap, linewidths=0)
        pc.set_array(np.array(values))
        ax.add_collection(pc)
        ax.set_xlim(self.origin[0], self.origin[0] + self.boxlen)
        ax.set_ylim(self.origin[1], self.origin[1] + self.boxlen)
        ax.set_aspect('equal')
        ax.set_xlabel(ha)
        ax.set_ylabel(va)
        label = f'log10({quantity})' if log else quantity
        plt.colorbar(pc, ax=ax, label=label)
        return pc
