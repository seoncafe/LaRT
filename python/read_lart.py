#!/usr/bin/env python3
"""Read a LaRT output FITS file given the input (*.in) filename.

Returns spectral arrays (wavelength, velocity, xfreq, Jout, Jin, Jabs, Jabs2)
and the Jmu output (mu array + Jmu) when present.

Peel-off observer outputs (`*_obs.fits.gz`, `*_obs2D.fits.gz`) are not
loaded here.

Usage (module):
    from read_lart import read_lart
    out = read_lart('t1tau2.in')
    print(out.Jout.shape, out.mu, out.Jmu.shape)

Usage (CLI):
    python read_lart.py t1tau2.in
"""
from __future__ import annotations

import os
import re
import sys
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from astropy.io import fits


# ---------------------------------------------------------------------------
# Container
# ---------------------------------------------------------------------------

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
    # --- Metadata ---
    spectrum_header: dict = field(default_factory=dict)
    jmu_header:      dict = field(default_factory=dict)

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
        return '\n'.join(lines)

    # ------------------------------------------------------------------
    # Plotting
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

        def _resolve(lim, lo, hi):
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
        xlim_eff = _resolve(xlim, xmin, xmax)
        ylim_eff = _resolve(ylim, ymin, ymax)

        # axis name -> (array, x-axis label, y-label variable, y-unit
        #                inside math mode)
        x_axes = {
            'velocity':  (self.velocity,   r'velocity [km s$^{-1}$]',
                          'v', r'\,(\mathrm{km/s})^{-1}'),
            'xfreq':     (self.xfreq,      r'$x = (\nu - \nu_0)/\Delta\nu_D$',
                          'x', ''),
            'wavelength':(self.wavelength, r'wavelength [$\mathrm{\AA}$]',
                          r'\lambda', r'\,\mathrm{\AA}^{-1}'),
        }
        if x not in x_axes:
            raise ValueError(f"x must be one of {list(x_axes)}")
        xvals, xlabel, yvar, yunit = x_axes[x]
        if xvals is None:
            raise RuntimeError(f"'{x}' column not in the FITS file.")

        # --- Convert spectrum normalisation from "per d(orig)" to
        #     "per d(requested)" so that area under the curve is preserved.
        i_unit  = int(self.spectrum_header.get('I_unit', 0) or 0)
        if i_unit == 1:
            orig_x = self.wavelength      # stored per dwavelength
        else:
            orig_x = self.xfreq           # stored per dxfreq (default)
        orig_dx = float(np.abs(orig_x[1] - orig_x[0]))
        req_dx  = float(np.abs(xvals[1]  - xvals[0]))
        factor  = orig_dx / req_dx        # multiply Jmu / Jout by this
        Jmu_y   = self.Jmu * factor
        Jout_y  = (self.Jout * factor) if self.Jout is not None else None

        ylabel_lines = rf'$J({yvar};\mu){yunit}$'
        ylabel_img   = ylabel_lines

        if ax is None:
            _, ax = plt.subplots(figsize=(7.0, 4.5))

        if kind == 'lines':
            norm = Normalize(vmin=float(self.mu.min()),
                             vmax=float(self.mu.max()))
            sm   = ScalarMappable(norm=norm, cmap=cmap)
            for i, mu in enumerate(self.mu):
                ax.plot(xvals, Jmu_y[i, :], color=sm.to_rgba(mu),
                        lw=1.0, label=f'$\\mu={mu:+.2f}$')
            if overplot_jout and Jout_y is not None:
                ax.plot(xvals, Jout_y, color='black', lw=1.5, ls='--',
                        label=r'$J_{\rm out}$ (all $\mu$)')
            sm.set_array([])
            cbar = ax.figure.colorbar(sm, ax=ax,
                                      label=r'$\mu = \cos\theta_z$')
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
            mu_edges = self.mu_edges
            mesh     = ax.pcolormesh(x_edges, mu_edges, Jmu_y,
                                     cmap=cmap, shading='flat')
            ax.figure.colorbar(mesh, ax=ax, label=ylabel_img)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(r'$\mu = \cos\theta_z$')
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


def read_lart(name: str) -> LaRTOutput:
    """Read a LaRT output FITS file.

    `name` may be:
      - a *.in input file ('t1tau2.in')
      - the stem of one     ('t1tau2')      — '.in' is added automatically
      - the FITS file path  ('t1tau2.fits.gz')
    """
    infile = resolve_input_file(name)
    if infile.endswith('.fits') or infile.endswith('.fits.gz'):
        # caller passed the FITS file directly
        fits_file = infile
        infile    = ''                              # no namelist available
    else:
        fits_file = fits_path_for(infile)
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
        spectrum_header = spec_hdr,
        jmu_header = jmu_hdr_d,
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


if __name__ == '__main__':
    _main(sys.argv)
