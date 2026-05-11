#!/usr/bin/env python3
"""Read a LaRT output FITS file given the input (*.in) filename.

Returns spectral arrays (wavelength, velocity, xfreq, Jout, Jin, Jabs, Jabs2)
and the Jmu output (mu array + Jmu) when present.

Peel-off observer outputs (`*_obs.fits.gz`, `*_obs_NNN.fits.gz`) are also
loaded when present.  Each peel-off file is one observer (one (alpha,beta,
gamma) direction); the routine collects them all into a list of
PeelObservation objects, and provides a plotting method to compare each
observer's spatially-averaged spectrum with the Jmu(mu = cos(beta)) curve.

Usage (module):
    from read_lart import read_lart
    out = read_lart('t1tau2.in')
    print(out.Jout.shape, out.mu, out.Jmu.shape)
    if out.peelings:
        out.plot_peel_jmu_compare()

Usage (CLI):
    python read_lart.py t1tau2.in
"""
from __future__ import annotations

import glob
import os
import re
import sys
from dataclasses import dataclass, field
from typing import Any, List, Optional, Tuple

import numpy as np
from lart_io import load_lart, find_lart_file, glob_lart


# ---------------------------------------------------------------------------
# Containers
# ---------------------------------------------------------------------------

@dataclass
class PeelObservation:
    """One peel-off observer (one (alpha, beta, gamma) direction).

    The 3D spectral cube is stored in numpy convention ``(nyim, nxim, nxfreq)``
    (Fortran-side shape was ``(nxfreq, nxim, nyim)``).  Both the scattered
    and the direct-light components are kept separately; ``cube`` is the sum.
    """

    file_name: str
    alpha:     float
    beta:      float
    gamma:     float
    distance:  float            # observer distance in code units
    dist_cm:   float            # distance unit in cm (1.0 if dimensionless)
    nphotons:  float
    sr_pix:    float            # solid angle of one pixel (steradian)
    nxim:      int
    nyim:      int
    cube:      np.ndarray       # scattered + direct, shape (nyim, nxim, nxfreq)
    scatt:     np.ndarray
    direc:     np.ndarray
    direc0:    Optional[np.ndarray] = None
    header:    dict = field(default_factory=dict)

    @property
    def mu(self) -> float:
        """Direction cosine of the observer (cos(beta) along the +z axis)."""
        return float(np.cos(np.deg2rad(self.beta)))

    def average_spectrum(self) -> np.ndarray:
        """1D spectrum: mean specific intensity over the 2D image pixels."""
        return self.cube.mean(axis=(0, 1))

    def velocity_moment_map(self,
                            velocity: np.ndarray,
                            order: int = 1,
                            component: str = 'all',
                            vel_range: Optional[Tuple[float, float]] = None
                            ) -> np.ndarray:
        r"""Velocity moment map, one value per spatial image pixel.

        For each pixel the spectrum :math:`I(v)` is integrated along the
        velocity axis.  Higher moments (order >= 1) are always
        intensity-normalized -- the raw, unnormalized integrals carry
        units that depend on the spectral binning and have no
        intuitive interpretation, so the helper does not return them.

        ======  =================================================================
        order   result for each pixel
        ======  =================================================================
        0       :math:`\int I\,dv`  (integrated intensity / mom-0)
        1       :math:`\langle v\rangle = \int I\,v\,dv \,/\, \int I\,dv`
                  intensity-weighted mean velocity [km/s]
        2       :math:`\sigma_v = \sqrt{\int I (v-\langle v\rangle)^2 dv
                  \,/\, \int I\,dv}` velocity dispersion [km/s]
        ======  =================================================================

        Parameters
        ----------
        velocity : (nxfreq,) array, km/s
            Velocity grid -- typically pass ``LaRTOutput.velocity``.
        order : {0, 1, 2}
            Moment order (default 1).
        component : {'all', 'scatt', 'direct'}
            Which spectrum to integrate.  ``'all'`` (default) is
            scattered + direct, matching ``self.cube``.
        vel_range : (lo, hi), optional
            Restrict integration to ``lo <= v <= hi`` (km/s).  Either
            bound may be None to leave it unbounded.

        Returns
        -------
        (nyim, nxim) ndarray
            Pixels with zero integrated intensity are returned as NaN
            for order >= 1.
        """
        if component == 'all':
            cube = self.cube
        elif component == 'scatt':
            cube = self.scatt
        elif component == 'direct':
            cube = self.direc
        else:
            raise ValueError(f"component must be 'all', 'scatt', or "
                             f"'direct', got {component!r}")
        if order not in (0, 1, 2):
            raise ValueError(f"order must be 0, 1, or 2, got {order}")

        v = np.asarray(velocity, dtype=float)
        if v.ndim != 1 or v.size != cube.shape[-1]:
            raise ValueError(
                f"velocity shape {v.shape} does not match cube nxfreq "
                f"axis {cube.shape[-1]}")

        if vel_range is not None:
            lo, hi = vel_range
            if lo is None: lo = -np.inf
            if hi is None: hi = +np.inf
            mask = (v >= lo) & (v <= hi)
            v = v[mask]
            cube = cube[..., mask]
            if v.size == 0:
                raise ValueError("vel_range excluded every velocity bin")

        dv = float(np.abs(v[1] - v[0])) if v.size >= 2 else 1.0
        m0 = cube.sum(axis=-1) * dv
        if order == 0:
            return m0

        m1 = (cube * v).sum(axis=-1) * dv
        with np.errstate(invalid='ignore', divide='ignore'):
            vmean = np.where(m0 > 0, m1 / m0, np.nan)
        if order == 1:
            return vmean

        # order == 2
        with np.errstate(invalid='ignore', divide='ignore'):
            m2c = (cube * (v - np.nan_to_num(vmean)[..., None])**2
                   ).sum(axis=-1) * dv
            return np.where(m0 > 0,
                            np.sqrt(np.maximum(m2c / m0, 0.0)),
                            np.nan)


# ---------------------------------------------------------------------------
# Clump-input dataclass + analysis
# ---------------------------------------------------------------------------

@dataclass
class ClumpsOutput:
    """A LaRT clump-input file loaded into per-clump arrays + attributes.

    Built by :func:`read_clumps` from either format (HDF5 or FITS).  Carries
    the analysis / plotting methods that operate purely on the clump
    population, so it can be used standalone (no LaRT output required) as
    well as embedded inside a :class:`LaRTOutput` via the ``clumps`` field.
    """

    clumps_file: str
    input_file:  str  = ''
    params:      dict = field(default_factory=dict)
    # core per-clump arrays (None if the column is absent from the file)
    x:           Optional[np.ndarray] = None
    y:           Optional[np.ndarray] = None
    z:           Optional[np.ndarray] = None
    vx:          Optional[np.ndarray] = None
    vy:          Optional[np.ndarray] = None
    vz:          Optional[np.ndarray] = None
    radius:      Optional[np.ndarray] = None        # R_CLUMP column
    rhokap:      Optional[np.ndarray] = None
    temp:        Optional[np.ndarray] = None
    # group / HDU attributes -- case-insensitive lookup via ``attr()``
    attrs:       dict = field(default_factory=dict)

    # ------------------------------------------------------------------
    # Convenience accessors
    # ------------------------------------------------------------------
    @property
    def pos(self) -> np.ndarray:
        """(N, 3) array of clump centres in code units."""
        return np.column_stack([self.x, self.y, self.z])

    @property
    def vel(self) -> Optional[np.ndarray]:
        """(N, 3) bulk velocities in km/s, or None if the VX/VY/VZ columns
        are absent."""
        if self.vx is None or self.vy is None or self.vz is None:
            return None
        return np.column_stack([self.vx, self.vy, self.vz])

    def attr(self, name: str, default: Any = None) -> Any:
        """Case-insensitive attribute lookup."""
        if name in self.attrs:
            return self.attrs[name]
        lname = name.lower()
        for k, v in self.attrs.items():
            if k.lower() == lname:
                return v
        return default

    @property
    def n_clumps(self) -> int:
        v = self.attr('N_CLUMPS')
        return int(v) if v is not None else (len(self.x) if self.x is not None else 0)

    @property
    def sphere_r(self) -> float:
        return float(self.attr('SPHERE_R', self.attr('RMAX', 1.0)))

    @property
    def r_min(self) -> float:
        """Inner placement radius (code units).  Reads the ``RMIN``
        attribute; defaults to 0."""
        return float(self.attr('RMIN', 0.0))

    def _radii_array(self) -> Optional[np.ndarray]:
        """Per-clump radii used by the f_vol / f_cov formulas: prefer the
        ``R_CLUMP`` column, fall back to the scalar ``CL_RAD`` attribute
        broadcast to every clump.  Returns ``None`` if neither is available.
        """
        if self.radius is not None:
            return np.asarray(self.radius, dtype=float)
        r0 = self.attr('CL_RAD')
        if r0 is not None and self.x is not None:
            return np.full(len(self.x), float(r0))
        return None

    def compute_f_vol(self) -> Optional[float]:
        """Recompute the realised volume filling factor from the loaded
        clump radii.  Mirrors LaRT's ``write_clumps_fits`` from-file branch:

        .. math::  f_{\\rm vol} = \\frac{\\sum r_i^3}{R^3 - r_{\\rm min}^3}

        where :math:`R` is ``sphere_r`` and :math:`r_{\\rm min}` is ``r_min``
        (both read from the file attributes).  Returns ``None`` if the
        clump radii cannot be determined.
        """
        radii = self._radii_array()
        if radii is None:
            return None
        R, rmin = self.sphere_r, self.r_min
        denom = max(R**3 - rmin**3, np.finfo(float).tiny)
        return float(np.sum(radii**3) / denom)

    def compute_f_cov(self) -> Optional[float]:
        r"""Recompute the realised line-of-sight covering factor.  Mirrors
        LaRT's ``write_clumps_fits`` from-file branch:

        .. math::  f_{\rm cov} = \frac{3}{4}
            \frac{\sum r_i^2}{R^2 + R\,r_{\rm min} + r_{\rm min}^2}.

        Reduces to :math:`\tfrac34 \sum r_i^2 / R^2` when ``r_min = 0``.
        Returns ``None`` if the clump radii cannot be determined.
        """
        radii = self._radii_array()
        if radii is None:
            return None
        R, rmin = self.sphere_r, self.r_min
        denom = max(R*R + R*rmin + rmin*rmin, np.finfo(float).tiny)
        return float(0.75 * np.sum(radii**2) / denom)

    @property
    def f_vol(self) -> Optional[float]:
        """Volume filling factor.  Returns the ``F_VOL`` attribute when
        present, otherwise falls back to :meth:`compute_f_vol`."""
        v = self.attr('F_VOL')
        return float(v) if v is not None else self.compute_f_vol()

    @property
    def f_cov(self) -> Optional[float]:
        """LOS covering factor.  Returns the ``F_COV`` attribute when
        present, otherwise falls back to :meth:`compute_f_cov`."""
        v = self.attr('F_COV')
        return float(v) if v is not None else self.compute_f_cov()

    def summary(self) -> str:
        lines = [f"Clumps file: {self.clumps_file}",
                 f"Input file : {self.input_file or '(none)'}",
                 f"N_clumps   : {self.n_clumps}",
                 f"sphere_R   : {self.sphere_r:.4g}"]
        if self.f_vol is not None:
            lines.append(f"f_vol      : {self.f_vol:.4g}")
        if self.f_cov is not None:
            lines.append(f"f_cov      : {self.f_cov:.4g}")
        if self.radius is not None:
            lines.append(f"R_clump    : min/max = {self.radius.min():.3e} / "
                         f"{self.radius.max():.3e}")
        if self.rhokap is not None:
            lines.append(f"RHOKAP     : min/max = {self.rhokap.min():.3e} / "
                         f"{self.rhokap.max():.3e}")
        if self.temp is not None:
            lines.append(f"TEMP       : min/max = {self.temp.min():.3e} / "
                         f"{self.temp.max():.3e}")
        return '\n'.join(lines)

    # ------------------------------------------------------------------
    # Plotting: clump cross-section slice
    # ------------------------------------------------------------------
    def plot_clump_slice(self,
                         axis: str = 'z',
                         value: float = 0.0,
                         colorby: Optional[str] = None,
                         fill: bool = True,
                         cmap: str = 'viridis',
                         vmin: Optional[float] = None,
                         vmax: Optional[float] = None,
                         linewidth: float = 0.5,
                         alpha: Optional[float] = None,
                         show_sphere: bool = True,
                         ax=None,
                         figsize: Tuple[float, float] = (7.0, 7.0),
                         title: Optional[str] = None,
                         show: bool = False,
                         savefig: Optional[str] = None):
        r"""Plot the cross-section of clumps that intersect a coordinate
        slice plane.

        Operates purely on the in-memory clump arrays loaded by
        :func:`read_clumps`; nothing about the on-disk format matters by
        the time this method is called.

        Parameters
        ----------
        axis : {'x', 'y', 'z'}, default 'z'
        value : float, default 0
        colorby : {'rhokap', 'radius', 'rcross', 'temp', 'none'}, optional
            Default ``'rhokap'`` if the RHOKAP column is loaded, else
            ``'radius'``.
        See :meth:`LaRTOutput.plot_clump_slice` for the rest of the
        parameters; the API is identical apart from the absent
        ``clumps_file`` override (already fixed by ``read_clumps``).
        """
        if axis not in ('x', 'y', 'z'):
            raise ValueError(f"axis must be 'x', 'y', or 'z', got {axis!r}")
        if colorby is not None and colorby not in (
                'none', 'rhokap', 'radius', 'rcross', 'temp'):
            raise ValueError(f"colorby must be one of 'none', 'rhokap', "
                             f"'radius', 'rcross', 'temp'; got {colorby!r}")
        if self.x is None or self.radius is None:
            raise RuntimeError(
                f"{self.clumps_file!r} did not contain the X/Y/Z and "
                f"R_CLUMP columns required for slice plotting.")

        try:
            from plot_clump_slice import slice_clumps, AXIS_TRIPLE
        except ImportError:
            here = os.path.dirname(os.path.abspath(__file__))
            if here not in sys.path:
                sys.path.insert(0, here)
            from plot_clump_slice import slice_clumps, AXIS_TRIPLE

        import matplotlib.pyplot as plt
        from matplotlib.patches import Circle
        from matplotlib.collections import PatchCollection

        a, b, rcross, idx = slice_clumps(self.pos, self.radius, axis, value)

        if colorby is None:
            colorby = 'rhokap' if self.rhokap is not None else 'radius'
        if colorby == 'none':
            color_vals, color_label = None, None
        elif colorby == 'rhokap':
            if self.rhokap is None:
                raise ValueError(
                    f"colorby='rhokap' requested but {os.path.basename(self.clumps_file)} "
                    f"has no RHOKAP column.")
            color_vals = np.asarray(self.rhokap)[idx]
            color_label = 'RHOKAP'
        elif colorby == 'radius':
            color_vals = np.asarray(self.radius)[idx]
            color_label = r'$R_{\rm clump}$'
        elif colorby == 'rcross':
            color_vals = rcross
            color_label = r'$r_{\rm cross}$'
        else:  # 'temp'
            if self.temp is None:
                raise ValueError(
                    f"colorby='temp' requested but {os.path.basename(self.clumps_file)} "
                    f"has no TEMP column.")
            color_vals = np.asarray(self.temp)[idx]
            color_label = 'Temperature  [K]'

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure

        R_sphere = self.sphere_r
        if show_sphere and abs(value) < R_sphere:
            r_outer = np.sqrt(R_sphere**2 - value**2)
            ax.add_patch(Circle((0.0, 0.0), r_outer, fill=False,
                                edgecolor='black', linestyle='--',
                                linewidth=1.0,
                                label=f'sphere @ {axis}={value:g}'))

        patches = [Circle((ai, bi), ri) for ai, bi, ri in zip(a, b, rcross)]
        fill_alpha = (alpha if alpha is not None
                      else (0.6 if fill else 1.0))

        if color_vals is None:
            col = PatchCollection(
                patches,
                facecolors=('tab:blue' if fill else 'none'),
                edgecolors='tab:blue',
                linewidths=linewidth,
                alpha=fill_alpha,
            )
            ax.add_collection(col)
        else:
            col = PatchCollection(
                patches, cmap=cmap,
                edgecolors=('face' if fill else None),
                linewidths=linewidth,
                alpha=fill_alpha,
            )
            col.set_array(np.asarray(color_vals))
            if vmin is not None or vmax is not None:
                col.set_clim(vmin=vmin, vmax=vmax)
            if not fill:
                col.set_facecolor('none')
            ax.add_collection(col)
            cb = fig.colorbar(col, ax=ax, fraction=0.046, pad=0.04)
            cb.set_label(color_label)

        _, _, _, lab_a, lab_b = AXIS_TRIPLE[axis]
        ax.set_xlabel(lab_a)
        ax.set_ylabel(lab_b)
        ax.set_aspect('equal')
        lim = R_sphere * 1.05
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)

        n_in_plane = len(rcross)
        if title is None:
            title = (f'{os.path.basename(self.clumps_file)}\n'
                     f'slice {axis} = {value:g}  |  '
                     f'{n_in_plane} / {len(self.pos)} clumps cross plane')
            if self.f_cov is not None and self.f_vol is not None:
                title += (f'\nf_cov={self.f_cov:.3f}  '
                          f'f_vol={self.f_vol:.4f}')
        ax.set_title(title, fontsize=10)
        if show_sphere and abs(value) < R_sphere:
            ax.legend(loc='lower right', fontsize=8)

        fig.tight_layout()
        if savefig:
            fig.savefig(savefig, dpi=150)
        if show:
            plt.show()
        return fig, ax


@dataclass
class LaRTOutput:
    """Holds arrays + headers from a LaRT output file (FITS or HDF5).

    In the "clumps-only" pre-run state (LaRT has not been executed yet but
    the namelist's ``par%clump_input_file`` exists on disk) the spectrum
    arrays are ``None`` and only ``input_file``, ``params``, plus the
    ``_clumps_file_path``-driven plotting hooks are usable.
    """

    fits_file:   str
    input_file:  str
    # --- Spectrum table (None in the clumps-only pre-run state) ---
    xfreq:       Optional[np.ndarray] = None
    velocity:    Optional[np.ndarray] = None
    wavelength:  Optional[np.ndarray] = None
    Jout:        Optional[np.ndarray] = None
    Jin:         Optional[np.ndarray] = None
    Jabs:        Optional[np.ndarray] = None
    Jabs2:       Optional[np.ndarray] = None
    # --- Jmu (optional) ---
    mu:          Optional[np.ndarray] = None       # bin centers
    mu_edges:    Optional[np.ndarray] = None       # bin edges (length nmu+1)
    Jmu:         Optional[np.ndarray] = None       # shape (nmu, nxfreq)
    nmu:         Optional[int]        = None
    mu_min:      Optional[float]      = None
    dmu:         Optional[float]      = None
    # --- Peel-off observers (optional) ---
    peelings:    List[PeelObservation] = field(default_factory=list)
    # --- Metadata ---
    spectrum_header: dict = field(default_factory=dict)
    jmu_header:      dict = field(default_factory=dict)
    params:          dict = field(default_factory=dict)   # parsed *.in namelist
    # --- Clump-input population (optional, lazy-loaded by plot_clump_slice) ---
    clumps:          Optional['ClumpsOutput'] = None

    @property
    def is_clumps_only(self) -> bool:
        """True if no LaRT output has been read (pre-run state with a
        clump input file)."""
        return self.Jout is None and self.xfreq is None

    def summary(self) -> str:
        nxfreq_str = str(len(self.xfreq)) if self.xfreq is not None else '(absent)'
        lines = [f"FITS file: {self.fits_file}",
                 f"Input    : {self.input_file}",
                 f"nxfreq   : {nxfreq_str}"]
        if self.is_clumps_only:
            cfile = self.params.get('clump_input_file') if self.params else None
            lines.insert(0,
                "Mode     : clump-input only (LaRT output not yet generated)")
            if cfile:
                lines.append(f"Clumps   : {cfile}")
            return '\n'.join(lines)
        for name in ('Jout', 'Jin', 'Jabs', 'Jabs2'):
            arr = getattr(self, name)
            lines.append(f"{name:8s} : {'present' if arr is not None else '(absent)'}")
        if self.Jmu is not None:
            lines.append(f"Jmu      : present  (nmu={self.nmu}, "
                         f"mu_min={self.mu_min:+.3f}, dmu={self.dmu:.4f})")
        else:
            lines.append("Jmu      : (absent)")
        if self.peelings:
            lines.append(f"peelings : {len(self.peelings)} observer(s)")
            for i, p in enumerate(self.peelings, 1):
                lines.append(f"   #{i:02d}: alpha={p.alpha:+6.1f} "
                             f"beta={p.beta:+6.1f} gamma={p.gamma:+6.1f}  "
                             f"mu={p.mu:+.4f}  ({os.path.basename(p.file_name)})")
        else:
            lines.append("peelings : (none)")
        return '\n'.join(lines)

    # ------------------------------------------------------------------
    # Plotting: Jmu only
    # ------------------------------------------------------------------
    def plot_jmu(self,
                 ax=None,
                 x: str = 'velocity',
                 kind: str = 'lines',
                 overplot_jout: bool = True,
                 cmap: str = 'viridis',
                 title: Optional[str] = None,
                 xlim: Optional[tuple] = None,
                 ylim: Optional[tuple] = None,
                 xmin: Optional[float] = None,
                 xmax: Optional[float] = None,
                 ymin: Optional[float] = None,
                 ymax: Optional[float] = None,
                 z_symmetry: bool = False,
                 show: bool = False,
                 savefig: Optional[str] = None):
        """Plot Jmu as a function of frequency/velocity/wavelength.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Existing axes to draw on; created if None.
        x : {'velocity', 'xfreq', 'wavelength'}
            Quantity for the x-axis (default 'velocity', km/s).
        kind : {'lines', 'image'}
            'lines'  -> one curve per mu bin, colored by mu;
            'image'  -> 2-D heatmap (x vs mu).
        overplot_jout : bool
            For kind='lines', also draw Jout in black for reference.
        cmap : str
            Colormap name.
        xlim, ylim : (lo, hi) tuple, optional
            Axis limits.  Either bound may be None to leave it autoscaled.
        xmin, xmax, ymin, ymax : float, optional
            Per-bound shorthand; merged with xlim/ylim (the per-bound
            value wins if both are given).
        z_symmetry : bool
            If True, average Jmu(x; +mu) and Jmu(x; -mu) before plotting,
            collapsing the mu axis onto |mu| in [0, +1].  Useful for
            simulations where the medium is symmetric about the z=0
            plane but the run was launched with par%xyz_symmetry=.false.;
            folding doubles the photon count per mu bin.  No effect when
            the data are already on [0, +1] (par%xyz_symmetry was on).
        title, show, savefig : convenience options.

        Returns
        -------
        ax : matplotlib.axes.Axes
        """
        if self.Jmu is None or self.mu is None:
            raise RuntimeError(
                "Jmu is not present in this output -- run with "
                "par%save_Jmu = .true. to enable.")
        # lazy import so the module is usable without matplotlib
        import matplotlib.pyplot as plt
        from matplotlib.cm import ScalarMappable
        from matplotlib.colors import Normalize

        xlim_eff = _resolve_lim(xlim, xmin, xmax)
        ylim_eff = _resolve_lim(ylim, ymin, ymax)

        xvals, xlabel, yvar, yunit = _x_axis_pick(self, x)
        factor = _spectral_jacobian(self, xvals)
        Jmu_y  = self.Jmu * factor
        Jout_y = (self.Jout * factor) if self.Jout is not None else None
        mu_arr = self.mu
        mu_edges_arr = self.mu_edges

        # --- Optional fold of +/-mu pairs: assume z=0 plane symmetry ---
        if z_symmetry and self.mu_min is not None and self.mu_min < -1e-6:
            # Group bins by |mu| (rounded to 4 decimals so near-symmetric
            # +/- pairs merge cleanly).
            abs_mu_key = np.round(np.abs(mu_arr), 4)
            keys, inv  = np.unique(abs_mu_key, return_inverse=True)
            new_jmu    = np.zeros((len(keys), Jmu_y.shape[1]))
            counts     = np.zeros(len(keys))
            for i, k in enumerate(inv):
                new_jmu[k] += Jmu_y[i]
                counts[k]  += 1
            new_jmu /= counts[:, None]
            Jmu_y  = new_jmu
            mu_arr = keys
            # rebuild mu_edges (uniform spacing from folded centres)
            if len(mu_arr) >= 2:
                dmu_new = float(mu_arr[1] - mu_arr[0])
            else:
                dmu_new = 1.0
            mu_edges_arr = np.concatenate(
                [[max(0.0, mu_arr[0] - 0.5*dmu_new)],
                 mu_arr + 0.5*dmu_new])
            mu_label_text = r'$|\mu| = |\cos\theta_z|$'
        else:
            mu_label_text = r'$\mu = \cos\theta_z$'

        ylabel_lines = rf'$J({yvar};\mu){yunit}$'
        ylabel_img   = ylabel_lines

        if ax is None:
            _, ax = plt.subplots(figsize=(7.0, 4.5))

        if kind == 'lines':
            norm = Normalize(vmin=float(mu_arr.min()),
                             vmax=float(mu_arr.max()))
            sm   = ScalarMappable(norm=norm, cmap=cmap)
            for i, mu in enumerate(mu_arr):
                ax.plot(xvals, Jmu_y[i, :], color=sm.to_rgba(mu),
                        lw=1.0, label=f'$\\mu={mu:+.2f}$')
            if overplot_jout and Jout_y is not None:
                ax.plot(xvals, Jout_y, color='black', lw=1.5, ls='--',
                        label=r'$J_{\rm out}$ (all $\mu$)')
            sm.set_array([])
            cbar = ax.figure.colorbar(sm, ax=ax, label=mu_label_text)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel_lines)
            if overplot_jout and Jout_y is not None:
                ax.legend(loc='best', fontsize=11, ncol=1,
                          handles=[ax.lines[-1]], frameon=False)

        elif kind == 'image':
            # build edges along x-axis (assume uniform spacing)
            dx       = float(xvals[1] - xvals[0])
            x_edges  = np.concatenate([xvals - 0.5*dx,
                                       [xvals[-1] + 0.5*dx]])
            mesh     = ax.pcolormesh(x_edges, mu_edges_arr, Jmu_y,
                                     cmap=cmap, shading='flat')
            ax.figure.colorbar(mesh, ax=ax, label=ylabel_img)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(mu_label_text)
        else:
            raise ValueError(f"kind must be 'lines' or 'image', got {kind!r}")

        if xlim_eff is not None:
            ax.set_xlim(*xlim_eff)
        if ylim_eff is not None:
            ax.set_ylim(*ylim_eff)
        if title is None:
            title = os.path.basename(self.fits_file)
        ax.set_title(title)
        ax.figure.tight_layout()
        if savefig:
            ax.figure.savefig(savefig, dpi=150)
        if show:
            plt.show()
        return ax

    # ------------------------------------------------------------------
    # Geometry helpers (used by Jpeel conversion)
    # ------------------------------------------------------------------
    def _source_surface_area(self) -> Tuple[float, str]:
        """Return ``(A_surface_in_code_units, label)`` matching the Fortran
        ``area`` factor used in Jout/Jmu normalization, but **without** the
        ``distance2cm**2`` part (which cancels with the peel-off scaling).

        Resolution order:
          1. AMR run (``par%use_amr_grid = .true.``)  -> A = 4π (rmax = 1).
          2. ``xy_periodic`` slab                      -> A = 2 * (2xmax)(2ymax).
          3. ``par%geometry = 'sphere'``               -> A = 4π·rmax².
          4. ``par%geometry = 'rectangle'/'box'``      -> A = 8(xy + yz + zx).
          5. Fallback (no input file): assume sphere with rmax = max(xmax,ymax,zmax).
        """
        h = self.spectrum_header
        xmax = float(h.get('XMAX', 1.0))
        ymax = float(h.get('YMAX', xmax))
        zmax = float(h.get('ZMAX', xmax))
        xy_per = bool(h.get('XY_PER', False))

        params = self.params or {}
        geometry = (params.get('geometry') or '').strip().lower()
        rmax = float(params.get('rmax') or 0.0)
        use_amr = bool(params.get('use_amr_grid', False))

        if use_amr:
            return 4.0 * np.pi, 'AMR (rmax=1)'
        if xy_per:
            return 2.0 * (2.0 * xmax) * (2.0 * ymax), 'slab'
        if geometry == 'sphere':
            if rmax <= 0.0:
                rmax = max(xmax, ymax, zmax)
            return 4.0 * np.pi * rmax**2, f'sphere (rmax={rmax:g})'
        if geometry in ('rectangle', 'box'):
            return 8.0 * (xmax*ymax + ymax*zmax + zmax*xmax), 'box'
        # No geometry from input file: prefer sphere when xmax==ymax==zmax,
        # otherwise box.
        if abs(xmax - ymax) < 1e-12 and abs(xmax - zmax) < 1e-12:
            rmax_eff = rmax if rmax > 0.0 else xmax
            return 4.0 * np.pi * rmax_eff**2, f'sphere (rmax={rmax_eff:g}, inferred)'
        return 8.0 * (xmax*ymax + ymax*zmax + zmax*xmax), 'box (inferred)'

    # ------------------------------------------------------------------
    # Plotting: peel-off vs Jmu comparison
    # ------------------------------------------------------------------
    def plot_peel_jmu_compare(self,
                              x: str = 'velocity',
                              component: str = 'all',
                              ncols: Optional[int] = None,
                              figsize: Optional[Tuple[float, float]] = None,
                              xlim: Optional[tuple] = None,
                              ylim: Optional[tuple] = None,
                              xmin: Optional[float] = None,
                              xmax: Optional[float] = None,
                              ymin: Optional[float] = None,
                              ymax: Optional[float] = None,
                              yscale: str = 'linear',
                              show: bool = False,
                              savefig: Optional[str] = None):
        r"""Compare each peel-off observer's spatially-summed spectrum with
        Jmu(mu = cos(beta)), in identical Jout/Jmu units (no arbitrary
        rescaling).

        For each peel-off observer in ``self.peelings`` this routine plots
        two curves on a separate axes panel:

        1. ``Jpeel(nu) = sum_pixels( I_peel ) * scale_for_peel`` where::

               scale_for_peel = (4 pi D^2) * dOmega_pix / (A_surface * 2 pi)

           This matches ``examples/sphere/plot_Jpeel_conversion.ipynb``.
           The bin_unit (dxfreq vs dwave) cancels because both Jmu and the
           peel-off cube are already normalized per-bin_unit by the
           Fortran code (see ``output_sum_rect.f90``).
        2. Jmu spectrum at the corresponding ``mu = cos(beta)``, linearly
           interpolated between the closest mu bin centres.

        The conversion is geometric and exact -- there is no fitted scale.

        Parameters
        ----------
        x : {'velocity', 'xfreq', 'wavelength'}
            Quantity for the x-axis (default 'velocity', km/s).
        component : {'all', 'scatt', 'direct'}
            Which part of the peel-off cube to sum.  ``'all'`` (default) is
            ``Scattered + Direct`` and matches the photon population that
            Jmu accumulates.  ``'scatt'`` uses only the scattered light.
            ``'direct'`` plots only the direct (unscattered) light.
        yscale : {'linear', 'log'}
            Y-axis scale.
        ncols : int, optional
            Number of subplot columns.  Default = ceil(sqrt(nobs)).
        figsize : (w, h), optional
            Figure size; default scales with the number of panels.
        xlim, ylim : tuples, optional
            Axis limits applied to every panel.  Either bound may be None.
        xmin, xmax, ymin, ymax : float, optional
            Per-bound overrides (per-bound value wins if both forms given).
        show, savefig : convenience options.

        Returns
        -------
        fig, axes : matplotlib.figure.Figure, np.ndarray of Axes
        """
        if not self.peelings:
            raise RuntimeError("No peel-off observers were found for this run.")
        if self.Jmu is None or self.mu is None:
            raise RuntimeError(
                "Jmu is not present in this output -- run with "
                "par%save_Jmu = .true. to enable Jmu output.")
        if component not in ('all', 'scatt', 'direct'):
            raise ValueError(f"component must be 'all', 'scatt', or "
                             f"'direct', got {component!r}")
        if yscale not in ('linear', 'log'):
            raise ValueError(f"yscale must be 'linear' or 'log', got "
                             f"{yscale!r}")

        import matplotlib.pyplot as plt

        xvals, xlabel, yvar, yunit = _x_axis_pick(self, x)
        factor = _spectral_jacobian(self, xvals)

        # Source surface area (code units squared) for Jpeel conversion.
        A_surface, geom_label = self._source_surface_area()

        nobs = len(self.peelings)
        if ncols is None:
            ncols = int(np.ceil(np.sqrt(nobs)))
        nrows = int(np.ceil(nobs / ncols))
        if figsize is None:
            figsize = (5.0 * ncols, 3.6 * nrows)
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
        axes_flat = axes.ravel()

        xlim_eff = _resolve_lim(xlim, xmin, xmax)
        ylim_eff = _resolve_lim(ylim, ymin, ymax)

        for i, peel in enumerate(self.peelings):
            ax = axes_flat[i]
            mu_target = peel.mu
            jmu_at_mu = _interp_jmu_at_mu(self.mu, self.Jmu, mu_target)
            jmu_y     = jmu_at_mu * factor
            # select cube component to sum
            if component == 'all':
                cube_use = peel.cube
                peel_label = r'$J_{\rm peel}$'
            elif component == 'scatt':
                cube_use = peel.scatt
                peel_label = r'$J_{\rm peel}^{\rm scatt}$'
            else:  # 'direct'
                cube_use = peel.direc
                peel_label = r'$J_{\rm peel}^{\rm direct}$'

            # Geometric conversion factor: no arbitrary scaling.
            #   scale_for_peel = 4*pi*D^2 * dOmega_pix / (A_surface * 2*pi)
            scale_for_peel = ((4.0 * np.pi * peel.distance**2) * peel.sr_pix
                              / (A_surface * 2.0 * np.pi))
            peel_y = cube_use.sum(axis=(0, 1)) * scale_for_peel * factor

            # Plot
            ax.plot(xvals, jmu_y,  color='black', lw=1.6,
                    label=rf'$J(\mu={mu_target:+.3f})$')
            ax.plot(xvals, peel_y, color='C3', lw=1.2, ls='--',
                    label=peel_label)

            ax.set_xlabel(xlabel)
            ax.set_ylabel(rf'$J({yvar};\mu){yunit}$')
            if xlim_eff is not None: ax.set_xlim(*xlim_eff)
            if ylim_eff is not None: ax.set_ylim(*ylim_eff)
            if yscale == 'log': ax.set_yscale('log')
            ax.legend(loc='best', fontsize=9, frameon=False)
            ax.set_title(f'$\\beta = {peel.beta:.1f}^{{\\circ}}$, '
                         f'$\\mu = {mu_target:+.3f}$',
                         fontsize=11)

        # blank out unused panels
        for j in range(nobs, nrows*ncols):
            axes_flat[j].axis('off')

        fig.suptitle(f'{os.path.basename(self.fits_file)}  [{geom_label}]',
                     fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        if savefig:
            fig.savefig(savefig, dpi=150)
        if show:
            plt.show()
        return fig, axes

    # ------------------------------------------------------------------
    # Plotting: velocity moment map(s) of peel-off image(s)
    # ------------------------------------------------------------------
    def plot_velocity_moment_map(self,
                                 observer: Optional[int] = None,
                                 order: int = 1,
                                 component: str = 'all',
                                 vel_range: Optional[Tuple[float, float]] = None,
                                 cmap: Optional[str] = None,
                                 vmin: Optional[float] = None,
                                 vmax: Optional[float] = None,
                                 symmetric: Optional[bool] = None,
                                 share_color: bool = True,
                                 transpose: bool = False,
                                 ncols: Optional[int] = None,
                                 figsize: Optional[Tuple[float, float]] = None,
                                 title: Optional[str] = None,
                                 show: bool = False,
                                 savefig: Optional[str] = None):
        r"""Plot velocity moment map(s) of the peel-off observer image(s).

        For every spatial pixel of each peel-off cube ``I(x, y, v)`` this
        computes the velocity moment via
        :py:meth:`PeelObservation.velocity_moment_map`.  Order 1 (default)
        produces the intensity-weighted mean velocity field
        :math:`\langle v\rangle = \int I\,v\,dv / \int I\,dv` [km/s];
        order 2 produces the velocity dispersion :math:`\sigma_v` [km/s];
        order 0 produces the integrated intensity :math:`\int I\,dv`.

        Parameters
        ----------
        observer : int, optional
            Index into ``self.peelings``.  If None (default), every
            peel-off observer is drawn as a subplot panel.
        order, component, vel_range :
            Forwarded to ``PeelObservation.velocity_moment_map``.
        cmap : str, optional
            Default 'RdBu_r' for order=1, 'inferno' for order=0/2.
        vmin, vmax : float, optional
            Color-scale bounds.  When both are None and
            ``symmetric=True`` (auto for order=1), a symmetric range
            about zero is chosen.
        symmetric : bool, optional
            Force symmetric color limits about 0.  Auto-set to True
            for order=1, otherwise False.
        share_color : bool
            Multi-observer mode: share one color scale across panels
            (default True).
        transpose : bool
            If True, swap image x and y axes (transpose the map and
            swap the extent + labels).  Default False.
        ncols, figsize, title, show, savefig : convenience options.

        Returns
        -------
        fig, axes : matplotlib.figure.Figure, Axes (or array of Axes)
        """
        if not self.peelings:
            raise RuntimeError("No peel-off observers were found for this run.")
        if self.velocity is None:
            raise RuntimeError(
                "velocity grid not available -- the FITS Spectrum table "
                "has no 'velocity' column.")

        import matplotlib.pyplot as plt
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        def _attach_cbar(im, ax, label):
            """Append a colorbar matched to the image axes height."""
            divider = make_axes_locatable(ax)
            cax     = divider.append_axes('right', size='4%', pad=0.08)
            return ax.figure.colorbar(im, cax=cax, label=label)

        if cmap is None:
            cmap = 'RdBu_r' if order == 1 else 'inferno'
        if symmetric is None:
            symmetric = (order == 1)

        if order == 0:
            cbar_label = r'$\int I\,dv$'
        elif order == 1:
            cbar_label = r'$\langle v\rangle$ [km s$^{-1}$]'
        else:
            cbar_label = r'$\sigma_v$ [km s$^{-1}$]'

        xlabel_str = 'image y' if transpose else 'image x'
        ylabel_str = 'image x' if transpose else 'image y'

        def _imshow_one(ax, mmap, peel, vmin_, vmax_):
            xlo, xhi, ylo, yhi = _peel_extent_code(peel)
            if transpose:
                arr = mmap.T
                ext = (ylo, yhi, xlo, xhi)
            else:
                arr = mmap
                ext = (xlo, xhi, ylo, yhi)
            return ax.imshow(arr, origin='lower', extent=ext,
                             cmap=cmap, vmin=vmin_, vmax=vmax_,
                             interpolation='nearest', aspect='equal')

        def _auto_limits(maps):
            if not symmetric:
                lo = min(float(np.nanmin(m)) for m in maps)
                hi = max(float(np.nanmax(m)) for m in maps)
                return lo, hi
            vabs = max(float(np.nanmax(np.abs(m))) for m in maps)
            if not np.isfinite(vabs) or vabs == 0:
                vabs = 1.0
            return -vabs, +vabs

        # --- single observer path --------------------------------------
        if observer is not None:
            peel = self.peelings[observer]
            mmap = peel.velocity_moment_map(self.velocity, order=order,
                                            component=component,
                                            vel_range=vel_range)
            lo_auto, hi_auto = _auto_limits([mmap])
            if vmin is None: vmin = lo_auto
            if vmax is None: vmax = hi_auto
            if figsize is None:
                figsize = (5.5, 4.6)
            fig, ax = plt.subplots(figsize=figsize)
            im = _imshow_one(ax, mmap, peel, vmin, vmax)
            _attach_cbar(im, ax, cbar_label)
            ax.set_xlabel(xlabel_str); ax.set_ylabel(ylabel_str)
            if title is None:
                title = (f"{os.path.basename(peel.file_name)}\n"
                         f"$\\alpha$={peel.alpha:+.1f}, "
                         f"$\\beta$={peel.beta:+.1f}, "
                         f"$\\gamma$={peel.gamma:+.1f}")
            ax.set_title(title, fontsize=11)
            fig.tight_layout()
            if savefig:
                fig.savefig(savefig, dpi=150)
            if show:
                plt.show()
            return fig, ax

        # --- multi-observer grid ---------------------------------------
        nobs = len(self.peelings)
        if ncols is None:
            ncols = int(np.ceil(np.sqrt(nobs)))
        nrows = int(np.ceil(nobs / ncols))
        if figsize is None:
            figsize = (4.6 * ncols, 3.9 * nrows)

        maps = [p.velocity_moment_map(self.velocity, order=order,
                                      component=component,
                                      vel_range=vel_range)
                for p in self.peelings]

        if share_color:
            lo_auto, hi_auto = _auto_limits(maps)
            vmin_eff = lo_auto if vmin is None else vmin
            vmax_eff = hi_auto if vmax is None else vmax
        else:
            vmin_eff = vmax_eff = None  # per-panel auto

        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
        axes_flat = axes.ravel()

        for i, (peel, mmap) in enumerate(zip(self.peelings, maps)):
            ax = axes_flat[i]
            if share_color:
                im = _imshow_one(ax, mmap, peel, vmin_eff, vmax_eff)
            else:
                lo, hi = _auto_limits([mmap])
                im = _imshow_one(ax, mmap, peel, lo, hi)
            _attach_cbar(im, ax, cbar_label)
            ax.set_title(f'$\\alpha$={peel.alpha:+.1f}, '
                         f'$\\beta$={peel.beta:+.1f}, '
                         f'$\\gamma$={peel.gamma:+.1f}', fontsize=10)
            ax.set_xlabel(xlabel_str); ax.set_ylabel(ylabel_str)

        for j in range(nobs, nrows * ncols):
            axes_flat[j].axis('off')

        if title is None:
            title = os.path.basename(self.fits_file)
        fig.suptitle(title, fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        if savefig:
            fig.savefig(savefig, dpi=150)
        if show:
            plt.show()
        return fig, axes

    # ------------------------------------------------------------------
    # Plotting: peel-off integrated-intensity map(s)
    # ------------------------------------------------------------------
    def plot_peeling_map(self,
                         observer: Optional[int] = None,
                         component: str = 'all',
                         vel_range: Optional[Tuple[float, float]] = None,
                         cmap: str = 'inferno',
                         vmin: Optional[float] = None,
                         vmax: Optional[float] = None,
                         log: bool = False,
                         transpose: bool = False,
                         share_color: bool = True,
                         ncols: Optional[int] = None,
                         figsize: Optional[Tuple[float, float]] = None,
                         title: Optional[str] = None,
                         show: bool = False,
                         savefig: Optional[str] = None):
        r"""Plot the spatial map of frequency-integrated intensity from each
        peel-off observer.

        For every spatial pixel of each peel-off cube ``I(x, y, v)`` this
        integrates over the velocity (frequency) axis and draws the
        resulting 2D map :math:`\int I\,dv`.

        Parameters
        ----------
        observer : int, optional
            Index into ``self.peelings``.  If None (default), every
            peel-off observer is drawn as a subplot panel.
        component : {'all', 'scatt', 'direct'}
            Which part of the peel-off cube to use.
        vel_range : (lo, hi), optional
            Restrict integration to ``lo <= v <= hi`` (km/s).
        cmap : str
            Matplotlib colormap (default 'inferno').
        vmin, vmax : float, optional
            Color-scale bounds.  If both None, set from the data.
        log : bool
            If True, use logarithmic color scaling (zero/negative pixels
            are masked).
        transpose : bool
            If True, swap image x and y axes.
        share_color : bool
            Multi-observer mode: share one color scale across panels.
        ncols, figsize, title, show, savefig : convenience options.

        Returns
        -------
        fig, ax(es) : matplotlib.figure.Figure, Axes (or array of Axes)
        """
        if not self.peelings:
            raise RuntimeError("No peel-off observers were found for this run.")
        if self.velocity is None:
            raise RuntimeError(
                "velocity grid not available -- the FITS Spectrum table "
                "has no 'velocity' column.")

        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm, Normalize
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        cbar_label = r'$\int I\,dv$'

        def _attach_cbar(im, ax, label):
            divider = make_axes_locatable(ax)
            cax     = divider.append_axes('right', size='4%', pad=0.08)
            return ax.figure.colorbar(im, cax=cax, label=label)

        xlabel_str = 'image y' if transpose else 'image x'
        ylabel_str = 'image x' if transpose else 'image y'

        def _imshow_one(ax, mmap, peel, norm_):
            xlo, xhi, ylo, yhi = _peel_extent_code(peel)
            if transpose:
                arr = mmap.T
                ext = (ylo, yhi, xlo, xhi)
            else:
                arr = mmap
                ext = (xlo, xhi, ylo, yhi)
            return ax.imshow(arr, origin='lower', extent=ext,
                             cmap=cmap, norm=norm_,
                             interpolation='nearest', aspect='equal')

        def _make_norm(maps):
            if log:
                # ignore zero/negative for log scaling
                vals = np.concatenate([m[m > 0].ravel() for m in maps])
                if vals.size == 0:
                    raise RuntimeError("Cannot use log scale: no positive "
                                       "pixel values found.")
                lo = float(vals.min()) if vmin is None else vmin
                hi = float(vals.max()) if vmax is None else vmax
                return LogNorm(vmin=lo, vmax=hi)
            lo = (min(float(np.nanmin(m)) for m in maps)
                  if vmin is None else vmin)
            hi = (max(float(np.nanmax(m)) for m in maps)
                  if vmax is None else vmax)
            return Normalize(vmin=lo, vmax=hi)

        def _moment0(peel):
            return peel.velocity_moment_map(self.velocity, order=0,
                                            component=component,
                                            vel_range=vel_range)

        # --- single observer ------------------------------------------
        if observer is not None:
            peel = self.peelings[observer]
            mmap = _moment0(peel)
            norm_ = _make_norm([mmap])
            if figsize is None:
                figsize = (5.5, 4.6)
            fig, ax = plt.subplots(figsize=figsize)
            im = _imshow_one(ax, mmap, peel, norm_)
            _attach_cbar(im, ax, cbar_label)
            ax.set_xlabel(xlabel_str); ax.set_ylabel(ylabel_str)
            if title is None:
                title = (f"{os.path.basename(peel.file_name)}\n"
                         f"$\\alpha$={peel.alpha:+.1f}, "
                         f"$\\beta$={peel.beta:+.1f}, "
                         f"$\\gamma$={peel.gamma:+.1f}")
            ax.set_title(title, fontsize=11)
            fig.tight_layout()
            if savefig:
                fig.savefig(savefig, dpi=150)
            if show:
                plt.show()
            return fig, ax

        # --- multi-observer grid --------------------------------------
        nobs = len(self.peelings)
        if ncols is None:
            ncols = int(np.ceil(np.sqrt(nobs)))
        nrows = int(np.ceil(nobs / ncols))
        if figsize is None:
            figsize = (4.6 * ncols, 3.9 * nrows)

        maps = [_moment0(p) for p in self.peelings]

        if share_color:
            shared_norm = _make_norm(maps)
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
        axes_flat = axes.ravel()

        for i, (peel, mmap) in enumerate(zip(self.peelings, maps)):
            ax = axes_flat[i]
            norm_ = shared_norm if share_color else _make_norm([mmap])
            im = _imshow_one(ax, mmap, peel, norm_)
            _attach_cbar(im, ax, cbar_label)
            ax.set_title(f'$\\alpha$={peel.alpha:+.1f}, '
                         f'$\\beta$={peel.beta:+.1f}, '
                         f'$\\gamma$={peel.gamma:+.1f}', fontsize=10)
            ax.set_xlabel(xlabel_str); ax.set_ylabel(ylabel_str)

        for j in range(nobs, nrows * ncols):
            axes_flat[j].axis('off')

        if title is None:
            title = os.path.basename(self.fits_file)
        fig.suptitle(title, fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        if savefig:
            fig.savefig(savefig, dpi=150)
        if show:
            plt.show()
        return fig, axes

    # ------------------------------------------------------------------
    # Plotting: peel-off spatially-integrated 1D spectrum
    # ------------------------------------------------------------------
    def plot_peeling_spectrum(self,
                              observer: Optional[int] = None,
                              x: str = 'velocity',
                              component: str = 'all',
                              mode: str = 'sum',
                              yscale: str = 'linear',
                              ax=None,
                              ncols: Optional[int] = None,
                              figsize: Optional[Tuple[float, float]] = None,
                              xlim: Optional[tuple] = None,
                              ylim: Optional[tuple] = None,
                              xmin: Optional[float] = None,
                              xmax: Optional[float] = None,
                              ymin: Optional[float] = None,
                              ymax: Optional[float] = None,
                              title: Optional[str] = None,
                              show: bool = False,
                              savefig: Optional[str] = None):
        """Plot the spatially-integrated (or averaged) 1D spectrum of each
        peel-off observer.

        Parameters
        ----------
        observer : int, optional
            Index into ``self.peelings``.  If None (default), every
            observer is drawn -- on a single set of axes (``ax`` given
            or only one observer) or as a grid of subplots otherwise.
        x : {'velocity', 'xfreq', 'wavelength'}
            X-axis quantity.
        component : {'all', 'scatt', 'direct'}
            Which cube to use.
        mode : {'sum', 'mean'}
            Spatial reduction over the (y, x) image axes.  ``'sum'``
            (default) integrates over pixels; ``'mean'`` returns the
            spatial average specific intensity.
        yscale : {'linear', 'log'}
        ax : matplotlib.axes.Axes, optional
            If provided, plot all observers as separate curves on this
            single axes (overrides the multi-panel grid).
        xlim, ylim, xmin/xmax/ymin/ymax : axis limit overrides.
        ncols, figsize, title, show, savefig : convenience options.

        Returns
        -------
        fig, ax(es) : matplotlib.figure.Figure, Axes (or array of Axes)
        """
        if not self.peelings:
            raise RuntimeError("No peel-off observers were found for this run.")
        if mode not in ('sum', 'mean'):
            raise ValueError(f"mode must be 'sum' or 'mean', got {mode!r}")
        if component not in ('all', 'scatt', 'direct'):
            raise ValueError(f"component must be 'all', 'scatt', or "
                             f"'direct', got {component!r}")
        if yscale not in ('linear', 'log'):
            raise ValueError(f"yscale must be 'linear' or 'log', got {yscale!r}")

        import matplotlib.pyplot as plt

        xvals, xlabel, _, _ = _x_axis_pick(self, x)
        xlim_eff = _resolve_lim(xlim, xmin, xmax)
        ylim_eff = _resolve_lim(ylim, ymin, ymax)

        ylabel = (r'$\sum_{\rm pix} I\,(\mathrm{pixel}\cdot \mathrm{flux})$'
                  if mode == 'sum'
                  else r'$\langle I\rangle_{\rm pix}$')

        def _spec(peel):
            if component == 'all':
                cube = peel.cube
            elif component == 'scatt':
                cube = peel.scatt
            else:
                cube = peel.direc
            if mode == 'sum':
                return cube.sum(axis=(0, 1))
            return cube.mean(axis=(0, 1))

        def _label(peel):
            return (rf'$\alpha={peel.alpha:+.0f}^\circ,\ '
                    rf'\beta={peel.beta:+.0f}^\circ,\ '
                    rf'\gamma={peel.gamma:+.0f}^\circ$')

        # --- single observer ------------------------------------------
        if observer is not None:
            peel = self.peelings[observer]
            if ax is None:
                if figsize is None:
                    figsize = (7.0, 4.5)
                fig, ax = plt.subplots(figsize=figsize)
            else:
                fig = ax.figure
            ax.plot(xvals, _spec(peel), lw=1.4)
            ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
            if yscale == 'log': ax.set_yscale('log')
            if xlim_eff is not None: ax.set_xlim(*xlim_eff)
            if ylim_eff is not None: ax.set_ylim(*ylim_eff)
            if title is None:
                title = (f"{os.path.basename(peel.file_name)}\n"
                         f"{_label(peel)}")
            ax.set_title(title, fontsize=11)
            fig.tight_layout()
            if savefig:
                fig.savefig(savefig, dpi=150)
            if show:
                plt.show()
            return fig, ax

        # --- all observers, single axes (overlay) --------------------
        if ax is not None:
            fig = ax.figure
            for peel in self.peelings:
                ax.plot(xvals, _spec(peel), lw=1.2, label=_label(peel))
            ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
            if yscale == 'log': ax.set_yscale('log')
            if xlim_eff is not None: ax.set_xlim(*xlim_eff)
            if ylim_eff is not None: ax.set_ylim(*ylim_eff)
            ax.legend(fontsize=9, frameon=False)
            if title is None:
                title = os.path.basename(self.fits_file)
            ax.set_title(title, fontsize=11)
            fig.tight_layout()
            if savefig:
                fig.savefig(savefig, dpi=150)
            if show:
                plt.show()
            return fig, ax

        # --- multi-observer grid --------------------------------------
        nobs = len(self.peelings)
        if ncols is None:
            ncols = int(np.ceil(np.sqrt(nobs)))
        nrows = int(np.ceil(nobs / ncols))
        if figsize is None:
            figsize = (4.8 * ncols, 3.4 * nrows)
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
        axes_flat = axes.ravel()

        for i, peel in enumerate(self.peelings):
            ax_i = axes_flat[i]
            ax_i.plot(xvals, _spec(peel), lw=1.4, color='C0')
            ax_i.set_xlabel(xlabel); ax_i.set_ylabel(ylabel)
            if yscale == 'log': ax_i.set_yscale('log')
            if xlim_eff is not None: ax_i.set_xlim(*xlim_eff)
            if ylim_eff is not None: ax_i.set_ylim(*ylim_eff)
            ax_i.set_title(_label(peel), fontsize=10)

        for j in range(nobs, nrows * ncols):
            axes_flat[j].axis('off')

        if title is None:
            title = os.path.basename(self.fits_file)
        fig.suptitle(title, fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        if savefig:
            fig.savefig(savefig, dpi=150)
        if show:
            plt.show()
        return fig, axes

    # convenience alias
    plot_peeling_spec = plot_peeling_spectrum

    # ------------------------------------------------------------------
    # Plotting: peel-off azimuthally-averaged radial profile
    # ------------------------------------------------------------------
    def plot_peeling_radial_profile(self,
                                    observer: Optional[int] = None,
                                    component: str = 'all',
                                    vel_range: Optional[Tuple[float, float]] = None,
                                    nbins: Optional[int] = None,
                                    rmax: Optional[float] = None,
                                    yscale: str = 'linear',
                                    xscale: str = 'linear',
                                    ax=None,
                                    ncols: Optional[int] = None,
                                    figsize: Optional[Tuple[float, float]] = None,
                                    xlim: Optional[tuple] = None,
                                    ylim: Optional[tuple] = None,
                                    xmin: Optional[float] = None,
                                    xmax: Optional[float] = None,
                                    ymin: Optional[float] = None,
                                    ymax: Optional[float] = None,
                                    title: Optional[str] = None,
                                    show: bool = False,
                                    savefig: Optional[str] = None):
        r"""Azimuthally averaged radial profile of the frequency-integrated
        peel-off map :math:`\int I\,dv`.

        For each peel-off observer the cube is integrated over the velocity
        axis to produce a 2D map; pixel distances from the image centre
        are then binned and the mean intensity within each annulus is
        computed.

        Parameters
        ----------
        observer : int, optional
            Index into ``self.peelings``.  If None (default), every
            observer is drawn -- on a single set of axes (``ax`` given)
            or as a grid of subplots otherwise.
        component : {'all', 'scatt', 'direct'}
        vel_range : (lo, hi), optional
            Restrict integration to ``lo <= v <= hi`` (km/s).
        nbins : int, optional
            Number of radial bins.  Default = ``nxim // 2``.
        rmax : float, optional
            Maximum radius (image-axis units).  Default = half the
            shorter image side.
        yscale, xscale : {'linear', 'log'}
        ax : matplotlib.axes.Axes, optional
            If given, plot all observers as separate curves on this single
            axes (overrides the multi-panel grid).
        xlim, ylim, xmin/xmax/ymin/ymax : axis limits.
        ncols, figsize, title, show, savefig : convenience options.

        Returns
        -------
        fig, ax(es) : matplotlib.figure.Figure, Axes (or array of Axes)
        """
        if not self.peelings:
            raise RuntimeError("No peel-off observers were found for this run.")
        if self.velocity is None:
            raise RuntimeError(
                "velocity grid not available -- the FITS Spectrum table "
                "has no 'velocity' column.")
        if yscale not in ('linear', 'log'):
            raise ValueError(f"yscale must be 'linear' or 'log', got {yscale!r}")
        if xscale not in ('linear', 'log'):
            raise ValueError(f"xscale must be 'linear' or 'log', got {xscale!r}")

        import matplotlib.pyplot as plt

        xlim_eff = _resolve_lim(xlim, xmin, xmax)
        ylim_eff = _resolve_lim(ylim, ymin, ymax)

        def _profile(peel: PeelObservation):
            mmap = peel.velocity_moment_map(self.velocity, order=0,
                                            component=component,
                                            vel_range=vel_range)
            xs, ys, _, _ = _peel_pixel_coords_code(peel)
            X, Y = np.meshgrid(xs, ys)
            R = np.sqrt(X*X + Y*Y)
            rmax_eff = (min(float(np.abs(xs).max()), float(np.abs(ys).max()))
                        if rmax is None else float(rmax))
            nb = peel.nxim // 2 if nbins is None else int(nbins)
            if nb < 1:
                raise ValueError("nbins must be >= 1")
            edges = np.linspace(0.0, rmax_eff, nb + 1)
            centers = 0.5 * (edges[:-1] + edges[1:])
            counts, _ = np.histogram(R.ravel(), bins=edges)
            weighted, _ = np.histogram(R.ravel(), bins=edges,
                                       weights=mmap.ravel())
            with np.errstate(invalid='ignore', divide='ignore'):
                prof = np.where(counts > 0, weighted / counts, np.nan)
            return centers, prof

        def _label(peel: PeelObservation):
            return (rf'$\alpha={peel.alpha:+.0f}^\circ,\ '
                    rf'\beta={peel.beta:+.0f}^\circ,\ '
                    rf'\gamma={peel.gamma:+.0f}^\circ$')

        xlabel = 'image radius'
        ylabel = r'$\langle\int I\,dv\rangle$ (annular avg)'

        def _finalize(ax_, with_legend=False, title_=None):
            ax_.set_xlabel(xlabel); ax_.set_ylabel(ylabel)
            if xscale == 'log': ax_.set_xscale('log')
            if yscale == 'log': ax_.set_yscale('log')
            if xlim_eff is not None: ax_.set_xlim(*xlim_eff)
            if ylim_eff is not None: ax_.set_ylim(*ylim_eff)
            if with_legend: ax_.legend(fontsize=9, frameon=False)
            if title_ is not None: ax_.set_title(title_, fontsize=11)

        # --- single observer ------------------------------------------
        if observer is not None:
            peel = self.peelings[observer]
            r, prof = _profile(peel)
            if ax is None:
                if figsize is None: figsize = (7.0, 4.5)
                fig, ax = plt.subplots(figsize=figsize)
            else:
                fig = ax.figure
            ax.plot(r, prof, lw=1.4)
            t = (f"{os.path.basename(peel.file_name)}\n{_label(peel)}"
                 if title is None else title)
            _finalize(ax, with_legend=False, title_=t)
            fig.tight_layout()
            if savefig: fig.savefig(savefig, dpi=150)
            if show: plt.show()
            return fig, ax

        # --- all observers, single axes (overlay) --------------------
        if ax is not None:
            fig = ax.figure
            for peel in self.peelings:
                r, prof = _profile(peel)
                ax.plot(r, prof, lw=1.2, label=_label(peel))
            t = os.path.basename(self.fits_file) if title is None else title
            _finalize(ax, with_legend=True, title_=t)
            fig.tight_layout()
            if savefig: fig.savefig(savefig, dpi=150)
            if show: plt.show()
            return fig, ax

        # --- multi-observer grid --------------------------------------
        nobs = len(self.peelings)
        if ncols is None:
            ncols = int(np.ceil(np.sqrt(nobs)))
        nrows = int(np.ceil(nobs / ncols))
        if figsize is None:
            figsize = (4.8 * ncols, 3.4 * nrows)
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
        axes_flat = axes.ravel()

        for i, peel in enumerate(self.peelings):
            r, prof = _profile(peel)
            ax_i = axes_flat[i]
            ax_i.plot(r, prof, lw=1.4, color='C0')
            _finalize(ax_i, with_legend=False, title_=_label(peel))

        for j in range(nobs, nrows * ncols):
            axes_flat[j].axis('off')

        if title is None:
            title = os.path.basename(self.fits_file)
        fig.suptitle(title, fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        if savefig: fig.savefig(savefig, dpi=150)
        if show: plt.show()
        return fig, axes

    # convenience alias
    plot_peeling_rprof = plot_peeling_radial_profile

    # ------------------------------------------------------------------
    # Plotting: clump cross-section slice
    # ------------------------------------------------------------------
    def _clumps_file_path(self) -> str:
        """Return the path of the corresponding clump file (HDF5 or FITS).

        Resolution order:
          1. ``par%clump_input_file`` from the namelist (if set) — resolved
             relative to the input-file directory when not absolute.
          2. ``<output-stem>_clumps.{h5,hdf5,fits.gz,fits}`` siblings of the
             main LaRT output file.
        Returns the first match that exists on disk; if nothing exists, the
        HDF5 default ``<output-stem>_clumps.h5`` is returned so a downstream
        FileNotFoundError points at the expected location.
        """
        cfile = self.params.get('clump_input_file') if self.params else None
        if isinstance(cfile, str) and cfile.strip():
            cfile = cfile.strip()
            if not os.path.isabs(cfile):
                indir = (os.path.dirname(os.path.abspath(self.input_file))
                         if self.input_file else
                         os.path.dirname(os.path.abspath(self.fits_file or '.')))
                cfile = os.path.join(indir or '.', cfile)
            if os.path.exists(cfile):
                return cfile
        # fall back to siblings of the output file
        base = self.fits_file
        for ext in ('.fits.gz', '.hdf5', '.fits', '.h5'):
            if base.lower().endswith(ext):
                base = base[:-len(ext)]
                break
        for ext in ('_clumps.h5', '_clumps.hdf5', '_clumps.fits.gz', '_clumps.fits'):
            cand = base + ext
            if os.path.exists(cand):
                return cand
        return base + '_clumps.h5'

    def plot_clump_slice(self, *args, clumps_file: Optional[str] = None, **kwargs):
        r"""Plot the cross-section of clumps that intersect a coordinate
        slice plane.

        Thin wrapper around :meth:`ClumpsOutput.plot_clump_slice`:
        loads the matching clump file via :func:`read_clumps` on first
        call (cached in ``self.clumps`` for reuse), then forwards all
        arguments.  Pass ``clumps_file=PATH`` to override the auto-located
        path.

        See :meth:`ClumpsOutput.plot_clump_slice` for the full parameter
        list.
        """
        if clumps_file is not None:
            co = read_clumps(clumps_file)
            if self.clumps is None:
                self.clumps = co
        else:
            if self.clumps is None:
                cfile = self._clumps_file_path()
                if not os.path.exists(cfile):
                    raise FileNotFoundError(
                        f"Clumps file not found: {cfile}\n"
                        f"  (derived from {self.fits_file}; pass clumps_file=... "
                        f"to override)")
                self.clumps = read_clumps(cfile)
            co = self.clumps
        return co.plot_clump_slice(*args, **kwargs)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _peel_pixel_coords_code(p: 'PeelObservation'):
    """Pixel-centre coordinates of a peel-off image in **code units**.

    The peel-off cube's FITS WCS stores ``CD2_2`` / ``CD3_3`` in
    degrees/pixel and the observer distance in code units (``DISTANCE``
    header / ``peel.distance``).  The transverse pixel size projected
    back to the source plane is therefore ``cd_deg * (pi/180) * D`` in
    the small-angle limit (always valid for peel-off geometry where the
    observer is far from the source).

    Returns
    -------
    xs : (nxim,) ndarray of pixel-centre x in code units
    ys : (nyim,) ndarray of pixel-centre y in code units
    dx, dy : float
        Pixel pitch along x / y in code units.
    """
    h    = p.header
    cd1d = abs(float(h.get('CD2_2', 1.0)))
    cd2d = abs(float(h.get('CD3_3', cd1d)))
    cx   = float(h.get('CRPIX2', 0.5*(p.nxim + 1)))
    cy   = float(h.get('CRPIX3', 0.5*(p.nyim + 1)))
    dx   = cd1d * (np.pi/180.0) * p.distance
    dy   = cd2d * (np.pi/180.0) * p.distance
    xs   = (np.arange(p.nxim) + 1.0 - cx) * dx
    ys   = (np.arange(p.nyim) + 1.0 - cy) * dy
    return xs, ys, dx, dy


def _peel_extent_code(p: 'PeelObservation'):
    """``(xlo, xhi, ylo, yhi)`` extent of a peel-off image in code units."""
    xs, ys, dx, dy = _peel_pixel_coords_code(p)
    return (xs[0] - 0.5*dx, xs[-1] + 0.5*dx,
            ys[0] - 0.5*dy, ys[-1] + 0.5*dy)


def _resolve_lim(lim, lo, hi):
    """Combine (lo, hi) tuple with explicit per-bound overrides."""
    if lim is None:
        lo_eff, hi_eff = None, None
    else:
        lo_eff, hi_eff = lim
    if lo is not None: lo_eff = lo
    if hi is not None: hi_eff = hi
    if lo_eff is None and hi_eff is None:
        return None
    return (lo_eff, hi_eff)


def _x_axis_pick(out: 'LaRTOutput', x: str):
    """Return (array, xlabel, yvar, yunit)."""
    x_axes = {
        'velocity':  (out.velocity,   r'velocity [km s$^{-1}$]',
                      'v', r'\,(\mathrm{km/s})^{-1}'),
        'xfreq':     (out.xfreq,      r'$x = (\nu - \nu_0)/\Delta\nu_D$',
                      'x', ''),
        'wavelength':(out.wavelength, r'wavelength [$\mathrm{\AA}$]',
                      r'\lambda', r'\,\mathrm{\AA}^{-1}'),
    }
    if x not in x_axes:
        raise ValueError(f"x must be one of {list(x_axes)}")
    xvals, xlabel, yvar, yunit = x_axes[x]
    if xvals is None:
        raise RuntimeError(f"'{x}' column not in the FITS file.")
    return xvals, xlabel, yvar, yunit


def _spectral_jacobian(out: 'LaRTOutput', xvals: np.ndarray) -> float:
    """Multiplier that converts the stored spectrum (per-d(orig)) to
    per-d(requested), so that the area under the curve is preserved."""
    i_unit = int(out.spectrum_header.get('I_unit', 0) or 0)
    orig_x = out.wavelength if i_unit == 1 else out.xfreq
    orig_dx = float(np.abs(orig_x[1] - orig_x[0]))
    req_dx  = float(np.abs(xvals[1]  - xvals[0]))
    return orig_dx / req_dx


def _interp_jmu_at_mu(mu_arr: np.ndarray,
                      Jmu:    np.ndarray,
                      mu_target: float) -> np.ndarray:
    """Linearly interpolate Jmu(mu, x) along mu, returning a 1-D spectrum."""
    mu_target = float(np.clip(mu_target, mu_arr.min(), mu_arr.max()))
    idx_above = int(np.searchsorted(mu_arr, mu_target))
    if idx_above >= len(mu_arr):
        return Jmu[-1].copy()
    if idx_above == 0:
        return Jmu[0].copy()
    a, b = idx_above - 1, idx_above
    w = (mu_target - mu_arr[a]) / (mu_arr[b] - mu_arr[a])
    return (1.0 - w) * Jmu[a] + w * Jmu[b]


# ---------------------------------------------------------------------------
# Input-file parsing
# ---------------------------------------------------------------------------

# match: par%key = value
#   - key is alphanumeric/underscore
#   - value can be quoted string, number, or .true./.false.
_PAR_RE = re.compile(
    r"par\s*%\s*(\w+)\s*=\s*(.+?)(?:\s*!|\s*$|,)",
    re.IGNORECASE,
)


def _strip_comment(line: str) -> str:
    """Return the part of `line` before any namelist comment ('!')."""
    in_str = False
    quote  = ''
    for i, ch in enumerate(line):
        if in_str:
            if ch == quote:
                in_str = False
        else:
            if ch in ("'", '"'):
                in_str = True
                quote  = ch
            elif ch == '!':
                return line[:i]
    return line


def _parse_value(raw: str):
    """Best-effort namelist value parser (string / float / int / bool)."""
    s = raw.strip().rstrip(',').strip()
    # quoted string
    if (s.startswith("'") and s.endswith("'")) or \
       (s.startswith('"') and s.endswith('"')):
        return s[1:-1]
    # logical
    low = s.lower().strip('.')
    if low in ('true', 't'):  return True
    if low in ('false', 'f'): return False
    # Fortran double-precision exponent (1.0d3 -> 1.0e3)
    s_num = s.replace('d', 'e').replace('D', 'e')
    # number
    try:
        if '.' in s_num or 'e' in s_num.lower():
            return float(s_num)
        return int(s_num)
    except ValueError:
        return s


def parse_input_file(infile: str) -> dict:
    """Parse a LaRT *.in namelist into a flat dict of par% values."""
    params = {}
    with open(infile, 'r') as f:
        for raw_line in f:
            line = _strip_comment(raw_line).strip()
            if not line or line.startswith('&') or line.startswith('/'):
                continue
            # there can be multiple par% entries on one line
            for m in re.finditer(r"par\s*%\s*(\w+)\s*=", line, re.IGNORECASE):
                key   = m.group(1).lower()
                start = m.end()
                # find end: next "par%" or end-of-line
                nxt = re.search(r",?\s*par\s*%", line[start:], re.IGNORECASE)
                end = start + nxt.start() if nxt else len(line)
                params[key] = _parse_value(line[start:end])
    return params


def resolve_input_file(name: str) -> str:
    """Accept 't1tau2.in', 't1tau2', or 't1tau2.fits.gz' and return the
    corresponding *.in path.  Adds the '.in' suffix if missing.  If the
    user passes a FITS path (or a stem that maps to one), returns that
    same path so read_lart can fall through to read it directly."""
    if os.path.isfile(name) and name.endswith('.in'):
        return name
    # already a FITS file?
    if name.endswith('.fits') or name.endswith('.fits.gz'):
        return name
    # try adding .in
    if not name.endswith('.in'):
        cand = name + '.in'
        if os.path.isfile(cand):
            return cand
        # also check if the FITS exists directly (stem only)
        for ext in ('.fits.gz', '.fits'):
            if os.path.isfile(name + ext):
                return name + ext
    return name  # let downstream raise FileNotFoundError with the original


def fits_path_for(infile: str) -> str:
    """Determine the output file path for a given LaRT input file.

    Accepts either format the simulation may have produced (HDF5 or FITS).
    Honours par%out_file if set in the namelist; otherwise falls back to
    ``<basename>.{h5,hdf5,fits.gz,fits}`` in the same directory as the
    input file, with HDF5 preferred (the LaRT v2 default).  The name is
    kept for backwards compatibility — it returns whatever exists,
    regardless of the actual on-disk format.
    """
    params = parse_input_file(infile)
    out_file = params.get('out_file', '').strip() if isinstance(
        params.get('out_file', ''), str) else ''
    indir = os.path.dirname(os.path.abspath(infile)) or '.'
    if out_file:
        if not os.path.isabs(out_file):
            out_file = os.path.join(indir, out_file)
        # If the namelist's value points at an existing file, return it
        # as-is.  Otherwise try sibling files with the LaRT-recognised
        # extensions: the simulation may have run with a different
        # par%file_format than the namelist string suggests.
        if os.path.exists(out_file):
            return out_file
        stem = out_file
        for ext in ('.fits.gz', '.hdf5', '.fits', '.h5'):
            if stem.lower().endswith(ext):
                stem = stem[: -len(ext)]
                break
        resolved = find_lart_file(stem)
        return resolved if resolved is not None else out_file
    # default: <basename>.{h5,fits.gz,...} in the same dir as the .in file
    base = os.path.splitext(os.path.basename(infile))[0]
    stem = os.path.join(indir, base)
    resolved = find_lart_file(stem)
    if resolved is not None:
        return resolved
    # nothing exists yet; return the HDF5 default so the caller's
    # FileNotFoundError points at the expected location.
    return stem + '.h5'


# ---------------------------------------------------------------------------
# Backend-agnostic reading via lart_io (handles both FITS and HDF5)
# ---------------------------------------------------------------------------

def _peel_file_list(main_file: str) -> List[str]:
    """Return sorted list of peel-off observer files corresponding to a
    main output file.

    Accepts both legacy (single observer): ``<base>_obs.{ext}``, and
    multi-observer: ``<base>_obs_001.{ext}``, ``<base>_obs_002.{ext}``
    (3-digit zero-padded suffix written by Fortran ``write(...,'(a,i3.3)')``).
    Both FITS and HDF5 extensions are searched.
    """
    # strip the recognised output extensions
    base = main_file
    for ext in ('.fits.gz', '.hdf5', '.fits', '.h5'):
        if base.lower().endswith(ext):
            base = base[: -len(ext)]
            break

    # multi-observer form: <base>_obs_NNN.<ext>
    candidates = sorted(glob_lart(base + '_obs_???'))
    # single-observer form: <base>_obs.<ext> (only if no NNN matches)
    if not candidates:
        single = find_lart_file(base, suffix='_obs')
        if single is not None:
            candidates.append(single)
    return candidates


def _read_peel_observation(fname: str) -> Optional[PeelObservation]:
    """Read one peel-off file (FITS or HDF5).  Returns None if the file
    is missing or does not contain a ``Scattered`` section."""
    if not os.path.exists(fname):
        return None
    lf = load_lart(fname)
    scattered = lf.section('Scattered')
    if scattered is None or scattered.data is None:
        return None
    scatt = np.asarray(scattered.data)               # (nyim, nxim, nxfreq)
    direct = lf.section('Direct')
    direc = (np.asarray(direct.data) if direct is not None and direct.data is not None
             else np.zeros_like(scatt))
    direct0 = lf.section('Direct0')
    direc0 = (np.asarray(direct0.data)
              if direct0 is not None and direct0.data is not None else None)
    nyim, nxim, _ = scatt.shape
    cdx = abs(float(scattered.attr('CD2_2', 0.0)))
    cdy = abs(float(scattered.attr('CD3_3', cdx)))
    sr_pix = cdx * cdy * (np.pi/180.0)**2
    # Build a header-ish dict for backwards-compat callers that use
    # ``po.header['KEY']``.  Provide case-variant aliases for the few
    # mixed-case keywords that LaRT writes (alpha/beta/gamma/nphotons).
    hdr = dict(scattered.attrs)
    for key in ('alpha', 'beta', 'gamma', 'nphotons'):
        if key in hdr and key.upper() not in hdr:
            hdr[key.upper()] = hdr[key]
    return PeelObservation(
        file_name = fname,
        alpha     = float(scattered.attr('alpha', 0.0)),
        beta      = float(scattered.attr('beta', 0.0)),
        gamma     = float(scattered.attr('gamma', 0.0)),
        distance  = float(scattered.attr('DISTANCE', 1.0)),
        dist_cm   = float(scattered.attr('DIST_CM', 1.0)),
        nphotons  = float(scattered.attr('nphotons', 0.0)),
        sr_pix    = sr_pix,
        nxim      = nxim,
        nyim      = nyim,
        cube      = scatt + direc,
        scatt     = scatt,
        direc     = direc,
        direc0    = direc0,
        header    = hdr,
    )


def _resolve_clump_path(name: str) -> Tuple[str, str, dict]:
    """Resolve a user-supplied name to ``(clumps_file, input_file, params)``.

    Accepts:
      * an existing clump file with an HDF5/FITS extension --- returned
        as-is (input_file/params blank).
      * an input file (``run.in`` or stem) --- ``par%clump_input_file`` is
        read from the namelist and resolved relative to the input dir.
      * a stem with a matching ``<stem>_clumps.{h5,hdf5,fits.gz,fits}`` ---
        returned with empty input metadata.

    Raises ``FileNotFoundError`` if no clump file can be located.
    """
    # 1. Direct clump file given?
    if name.lower().endswith(('.h5', '.hdf5', '.fits', '.fits.gz')) \
            and os.path.exists(name):
        return name, '', {}

    # 2. Try to resolve as an input file
    infile = resolve_input_file(name)
    params: dict = {}
    if os.path.isfile(infile) and infile.endswith('.in'):
        try:
            params = parse_input_file(infile)
        except OSError:
            params = {}
        cval = params.get('clump_input_file', '')
        if isinstance(cval, str) and cval.strip():
            cpath = cval.strip()
            if not os.path.isabs(cpath):
                cpath = os.path.join(
                    os.path.dirname(os.path.abspath(infile)) or '.', cpath)
            if os.path.exists(cpath):
                return cpath, infile, params

    # 3. Try '<stem>_clumps.*' siblings
    stem = name
    for ext in ('.in', '.fits.gz', '.hdf5', '.fits', '.h5'):
        if stem.lower().endswith(ext):
            stem = stem[: -len(ext)]
            break
    resolved = find_lart_file(stem, suffix='_clumps')
    if resolved is not None:
        return resolved, (infile if infile.endswith('.in') and os.path.isfile(infile) else ''), params

    raise FileNotFoundError(
        f"No clump file could be located from {name!r} "
        f"(checked direct path, par%clump_input_file in {infile!r}, "
        f"and <stem>_clumps.{{h5,hdf5,fits.gz,fits}}).")


def read_clumps(name: str) -> ClumpsOutput:
    """Load a LaRT clump-input file (HDF5 or FITS) into a :class:`ClumpsOutput`.

    ``name`` may be:
      * a clump file path: ``'foo_clumps.h5'`` or ``'foo_clumps.fits.gz'``
      * an input file: ``'run.in'`` -- reads ``par%clump_input_file`` and
        loads it (relative paths resolved against the input-file directory)
      * a stem: ``'run'`` -- tries ``'run.in'`` first, then a sibling
        ``run_clumps.{h5,hdf5,fits.gz,fits}``.

    Returns a :class:`ClumpsOutput` with the per-clump arrays and group
    attributes populated.
    """
    cpath, infile, params = _resolve_clump_path(name)
    lf = load_lart(cpath)
    sec = lf.section('Clumps') or (lf.sections[0] if lf.sections else None)
    if sec is None:
        raise RuntimeError(f"No clump section found in {cpath!r}; "
                           f"available sections: {[s.name for s in lf.sections]}")

    def _col(col_name: str) -> Optional[np.ndarray]:
        v = sec.col(col_name)
        return np.asarray(v) if v is not None else None

    attrs = dict(sec.attrs)
    for k in list(attrs):
        if k.upper() not in attrs:
            attrs[k.upper()] = attrs[k]

    return ClumpsOutput(
        clumps_file = cpath,
        input_file  = infile,
        params      = params,
        x           = _col('X'),
        y           = _col('Y'),
        z           = _col('Z'),
        vx          = _col('VX'),
        vy          = _col('VY'),
        vz          = _col('VZ'),
        radius      = _col('R_CLUMP'),
        rhokap      = _col('RHOKAP'),
        temp        = _col('TEMP'),
        attrs       = attrs,
    )


# Backward-compat singular alias.
read_clump = read_clumps


def read_lart(name: str) -> LaRTOutput:
    """Read a LaRT output file (FITS or HDF5).

    `name` may be:
      - a *.in input file       ('t1tau2.in')
      - the stem of one         ('t1tau2')        — '.in' is added automatically
      - a LaRT output file path ('t1tau2.fits.gz' or 't1tau2.h5')

    Peel-off files matching ``<base>_obs*.<ext>`` (with `.<ext>` in
    `.h5`, `.hdf5`, `.fits.gz`, `.fits`) are auto-loaded into
    ``LaRTOutput.peelings`` when present.
    """
    infile = resolve_input_file(name)
    lower  = infile.lower()
    if lower.endswith(('.fits', '.fits.gz', '.h5', '.hdf5')):
        # caller passed the output file directly
        fits_file = infile
        infile    = ''                              # no namelist available
        params    = {}
    else:
        fits_file = fits_path_for(infile)
        try:
            params = parse_input_file(infile)
        except OSError:
            params = {}
    if not os.path.exists(fits_file):
        # Clumps-only fallback: LaRT has not been run yet, but the namelist
        # points at an existing par%clump_input_file.  Delegate to
        # read_clumps() so we get the full ClumpsOutput (per-clump arrays
        # + attrs + plot_clump_slice method) without duplicating logic.
        try:
            clumps = read_clumps(infile) if infile else None
        except (FileNotFoundError, OSError):
            clumps = None
        if clumps is not None:
            print(f"read_lart: output {fits_file} not yet present; "
                  f"loading clump input only ({clumps.clumps_file}).")
            return LaRTOutput(
                fits_file = fits_file,
                input_file = infile,
                params = params,
                clumps = clumps,
            )
        raise FileNotFoundError(
            f"Output file not found: {fits_file}\n"
            f"  (input was {infile})"
        )

    lf = load_lart(fits_file)
    spec = lf.section('Spectrum')
    if spec is None:
        raise RuntimeError(
            f"No 'Spectrum' section found in {fits_file}.  "
            f"Available sections: {[s.name for s in lf.sections]}"
        )

    xfreq      = spec.col('Xfreq')
    velocity   = spec.col('velocity')
    wavelength = spec.col('wavelength')
    Jout       = spec.col('Jout')
    Jin        = spec.col('Jin')
    Jabs       = spec.col('Jabs')
    Jabs2      = spec.col('Jabs2')
    # Backwards-compat header dict: include both original-case attrs and
    # uppercase aliases so legacy callers using either style still work.
    spec_hdr = dict(spec.attrs)
    for k in list(spec_hdr):
        if k.upper() not in spec_hdr:
            spec_hdr[k.upper()] = spec_hdr[k]

    # Jmu (optional)
    jmu_sec = lf.section('Jmu')
    mu_arr      = None
    mu_edges    = None
    Jmu_arr     = None
    nmu         = None
    mu_min      = None
    dmu         = None
    jmu_hdr_d   = {}

    if jmu_sec is not None and jmu_sec.data is not None:
        jmu_hdr_d = dict(jmu_sec.attrs)
        for k in list(jmu_hdr_d):
            if k.upper() not in jmu_hdr_d:
                jmu_hdr_d[k.upper()] = jmu_hdr_d[k]
        data = np.asarray(jmu_sec.data)
        # Fortran-order (nxfreq, nmu) is read by astropy as numpy shape
        # (nmu, nxfreq).  HDF5 follows the same convention thanks to the
        # writer-side dimension ordering in iofile_mod.
        nmu = jmu_sec.attr('nmu', None)
        nmu = int(nmu) if nmu is not None else None
        nxfreq = len(xfreq) if xfreq is not None else None
        if nmu is None:
            if nxfreq is not None and data.shape[0] == nxfreq:
                nmu = data.shape[1]
            elif nxfreq is not None and data.shape[1] == nxfreq:
                nmu = data.shape[0]
            else:
                nmu = min(data.shape)
        if data.shape == (nmu, data.size // nmu):
            Jmu_arr = data
        elif data.shape == (data.size // nmu, nmu):
            Jmu_arr = data.T
        else:
            raise RuntimeError(
                f"Unexpected Jmu shape {data.shape}; nmu={nmu}, "
                f"nxfreq={nxfreq}."
            )
        mu_min  = float(jmu_sec.attr('mu_min', -1.0))
        dmu     = float(jmu_sec.attr('dmu', 2.0/nmu))
        mu_arr  = mu_min + (np.arange(nmu) + 0.5) * dmu
        mu_edges = mu_min + np.arange(nmu + 1) * dmu

    # Peel-off observers (optional)
    peelings: List[PeelObservation] = []
    for fname in _peel_file_list(fits_file):
        po = _read_peel_observation(fname)
        if po is not None:
            peelings.append(po)
    # sort by beta for stable ordering when plotted
    peelings.sort(key=lambda p: (p.beta, p.alpha, p.gamma))

    return LaRTOutput(
        fits_file = fits_file,
        input_file = infile,
        xfreq = xfreq,
        velocity = velocity,
        wavelength = wavelength,
        Jout = Jout,
        Jin  = Jin,
        Jabs = Jabs,
        Jabs2 = Jabs2,
        mu = mu_arr,
        mu_edges = mu_edges,
        Jmu = Jmu_arr,
        nmu = nmu,
        mu_min = mu_min,
        dmu = dmu,
        peelings = peelings,
        spectrum_header = spec_hdr,
        jmu_header = jmu_hdr_d,
        params = params,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _main(argv):
    if len(argv) < 2:
        print(__doc__)
        sys.exit(1)
    out = read_lart(argv[1])
    print(out.summary())
    if out.Jmu is not None and out.Jout is not None:
        # show per-bin Jmu_max vs Jout_max as a sanity check.
        Jout_max = float(np.nanmax(out.Jout))
        print(f"\nJout_max = {Jout_max:.4e}")
        print("Jmu_max per mu bin / Jout_max:")
        for i, mu in enumerate(out.mu):
            jmax = float(np.nanmax(out.Jmu[i, :]))
            ratio = jmax / Jout_max if Jout_max > 0 else float('nan')
            print(f"  mu={mu:+.3f}  Jmu_max={jmax:.4e}  ratio={ratio:.3f}")
    if out.peelings:
        print("\nPeel-off observers:")
        for p in out.peelings:
            avg = p.average_spectrum()
            print(f"  beta={p.beta:+6.1f}  mu={p.mu:+.4f}  "
                  f"<I_peel>_max = {float(np.nanmax(avg)):.4e}  "
                  f"({os.path.basename(p.file_name)})")


if __name__ == '__main__':
    _main(sys.argv)
