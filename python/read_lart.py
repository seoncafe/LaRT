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
from typing import List, Optional, Tuple

import numpy as np
from astropy.io import fits


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


@dataclass
class LaRTOutput:
    """Holds arrays + headers from a LaRT FITS output."""

    fits_file:   str
    input_file:  str
    # --- Spectrum table (always present) ---
    xfreq:       np.ndarray
    velocity:    np.ndarray
    wavelength:  np.ndarray
    Jout:        np.ndarray
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

    def summary(self) -> str:
        nxfreq_str = str(len(self.xfreq)) if self.xfreq is not None else '(absent)'
        lines = [f"FITS file: {self.fits_file}",
                 f"Input    : {self.input_file}",
                 f"nxfreq   : {nxfreq_str}"]
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


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

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
    """Determine the output FITS file path for a given input file.

    Honours par%out_file if set; otherwise falls back to <basename>.fits.gz
    in the same directory as the input file.
    """
    params = parse_input_file(infile)
    out_file = params.get('out_file', '').strip() if isinstance(
        params.get('out_file', ''), str) else ''
    indir = os.path.dirname(os.path.abspath(infile)) or '.'
    if out_file:
        # may be absolute or relative
        if not os.path.isabs(out_file):
            out_file = os.path.join(indir, out_file)
        return out_file
    # default: base name of *.in -> .fits.gz (in same dir)
    base = os.path.splitext(os.path.basename(infile))[0]
    return os.path.join(indir, base + '.fits.gz')


# ---------------------------------------------------------------------------
# FITS reading
# ---------------------------------------------------------------------------

def _find_hdu(hdul, extname: str):
    """Return the first HDU with matching EXTNAME (case-insensitive), else None."""
    for h in hdul:
        ext = h.header.get('EXTNAME', '')
        if ext.strip().lower() == extname.lower():
            return h
    return None


def _peel_file_list(fits_file: str) -> List[str]:
    """Return sorted list of peel-off observer files corresponding to a
    main output FITS file.

    Accepts both legacy (single observer): ``<base>_obs.fits.gz``, and
    multi-observer: ``<base>_obs_001.fits.gz``, ``<base>_obs_002.fits.gz``
    (3-digit zero-padded suffix written by Fortran `write(...,'(a,i3.3)')`).
    """
    # strip .fits / .fits.gz
    base = fits_file
    for ext in ('.fits.gz', '.fits'):
        if base.endswith(ext):
            base = base[:-len(ext)]
            break

    candidates = []
    # multi-observer form: <base>_obs_NNN.fits.gz
    candidates.extend(sorted(glob.glob(base + '_obs_???.fits.gz')))
    candidates.extend(sorted(glob.glob(base + '_obs_???.fits')))
    # single-observer form: <base>_obs.fits.gz (only if no NNN matches)
    if not candidates:
        for cand in (base + '_obs.fits.gz', base + '_obs.fits'):
            if os.path.exists(cand):
                candidates.append(cand)
    return candidates


def _read_peel_observation(fname: str) -> Optional[PeelObservation]:
    """Read one peel-off FITS file. Returns None if the file is missing or
    does not contain a 'Scattered' image extension."""
    if not os.path.exists(fname):
        return None
    with fits.open(fname) as hdul:
        scatt_hdu = _find_hdu(hdul, 'Scattered')
        if scatt_hdu is None:
            return None
        h = scatt_hdu.header
        scatt = np.asarray(scatt_hdu.data)            # (nyim, nxim, nxfreq)
        direc_hdu = _find_hdu(hdul, 'Direct')
        direc = (np.asarray(direc_hdu.data)
                 if direc_hdu is not None else np.zeros_like(scatt))
        direc0_hdu = _find_hdu(hdul, 'Direct0')
        direc0 = (np.asarray(direc0_hdu.data)
                  if direc0_hdu is not None else None)
        nyim, nxim, _ = scatt.shape
        cdx = abs(float(h.get('CD2_2', 0.0)))
        cdy = abs(float(h.get('CD3_3', cdx)))
        sr_pix = cdx * cdy * (np.pi/180.0)**2
        return PeelObservation(
            file_name = fname,
            alpha     = float(h.get('alpha', 0.0)),
            beta      = float(h.get('beta', 0.0)),
            gamma     = float(h.get('gamma', 0.0)),
            distance  = float(h.get('DISTANCE', 1.0)),
            dist_cm   = float(h.get('DIST_CM', 1.0)),
            nphotons  = float(h.get('nphotons', 0.0)),
            sr_pix    = sr_pix,
            nxim      = nxim,
            nyim      = nyim,
            cube      = scatt + direc,
            scatt     = scatt,
            direc     = direc,
            direc0    = direc0,
            header    = dict(h),
        )


def read_lart(name: str) -> LaRTOutput:
    """Read a LaRT output FITS file.

    `name` may be:
      - a *.in input file ('t1tau2.in')
      - the stem of one     ('t1tau2')      — '.in' is added automatically
      - the FITS file path  ('t1tau2.fits.gz')

    Peel-off files matching `<base>_obs*.fits.gz` are auto-loaded into
    ``LaRTOutput.peelings`` when present.
    """
    infile = resolve_input_file(name)
    if infile.endswith('.fits') or infile.endswith('.fits.gz'):
        # caller passed the FITS file directly
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
        raise FileNotFoundError(
            f"Output FITS file not found: {fits_file}\n"
            f"  (input was {infile})"
        )

    with fits.open(fits_file) as hdul:
        spec = _find_hdu(hdul, 'Spectrum')
        if spec is None:
            raise RuntimeError(
                f"No 'Spectrum' HDU found in {fits_file}.  "
                f"Available HDUs: "
                f"{[h.header.get('EXTNAME','') for h in hdul]}"
            )

        cols  = spec.data.columns.names
        get   = lambda name: (np.asarray(spec.data[name])
                              if name in cols else None)
        xfreq      = get('Xfreq')
        velocity   = get('velocity')
        wavelength = get('wavelength')
        Jout       = get('Jout')
        Jin        = get('Jin')
        Jabs       = get('Jabs')
        Jabs2      = get('Jabs2')
        spec_hdr   = dict(spec.header)

        # Jmu (optional)
        jmu_hdu  = _find_hdu(hdul, 'Jmu')
        mu_arr      = None
        mu_edges    = None
        Jmu_arr     = None
        nmu         = None
        mu_min      = None
        dmu         = None
        jmu_hdr_d   = {}

        if jmu_hdu is not None:
            jmu_hdr_d = dict(jmu_hdu.header)
            data      = np.asarray(jmu_hdu.data)
            # Fortran-order (nxfreq, nmu) is read by astropy as numpy
            # shape (nmu, nxfreq).
            nmu = jmu_hdu.header.get('nmu', None)
            nmu = int(nmu) if nmu is not None else None
            nxfreq = len(xfreq) if xfreq is not None else None
            if nmu is None:
                # try to infer from shape + nxfreq
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
            mu_min  = float(jmu_hdu.header.get('mu_min', -1.0))
            dmu     = float(jmu_hdu.header.get('dmu', 2.0/nmu))
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
