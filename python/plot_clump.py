#!/usr/bin/env python3
"""
plot_clump.py  --  Quick analysis plots for the LaRT_v2.00 clumpy-sphere examples.

Usage
-----
  python plot_clump.py                  # plots all available output files
  python plot_clump.py clump_tau3_fcov5 # plots a single run by base name

Requires: astropy, matplotlib, numpy
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from lart_io import load_lart, find_lart_file

# -----------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------

def load_spectrum(fname):
    """Return (velocity, Jout, Jin) arrays from a LaRT spectrum file
    (FITS .fits.gz or HDF5 .h5)."""
    lf = load_lart(fname)
    spec = lf.section('Spectrum')
    vel  = spec.col('velocity')
    Jout = spec.col('Jout')
    Jin  = spec.col('Jin')
    return vel, Jout, Jin


def load_clumps(fname):
    """Return (x, y, z, vx, vy, vz) arrays from a *_clumps.* file."""
    lf = load_lart(fname)
    # The clumps file has a single table section (named 'Clumps', or just
    # 'section_001' if the writer didn't set EXTNAME — fall back to the
    # first section in either case).
    sec = lf.section('Clumps') or (lf.sections[0] if lf.sections else None)
    if sec is None:
        raise ValueError(f'No section found in {fname!r}')
    return sec.col('X'), sec.col('Y'), sec.col('Z'), \
           sec.col('VX'), sec.col('VY'), sec.col('VZ')


def load_peeloff(fname):
    """Return peel-off image cube (nyim, nxim, nxfreq) from a peel-off file."""
    lf = load_lart(fname)
    sec = lf.section('Scattered') or (lf.sections[0] if lf.sections else None)
    if sec is None or sec.data is None:
        raise ValueError(f'No scattered cube found in {fname!r}')
    return sec.data   # shape: (nyim, nxim, nxfreq) or (nyim, nxim)


def load_clump_header(fname):
    """Return attrs dict (case-insensitive lookup) of the clumps file."""
    lf = load_lart(fname)
    sec = lf.section('Clumps') or (lf.sections[0] if lf.sections else None)
    return sec.attrs if sec is not None else {}


# -----------------------------------------------------------------------
# Plots for each run
# -----------------------------------------------------------------------

def plot_run(base, save_pdf=False):
    spec_file  = find_lart_file(base)
    clump_file = find_lart_file(base, suffix='_clumps')
    obs_file   = find_lart_file(base, suffix='_obs')

    if spec_file is None:
        print(f'  [skip] no spectrum file found for {base!r}')
        return

    fig, axes = plt.subplots(1, 3 if clump_file is not None else 2,
                             figsize=(14, 4))
    fig.suptitle(os.path.basename(base), fontsize=12)

    # --- Spectrum ---
    ax = axes[0]
    vel, Jout, Jin = load_spectrum(spec_file)
    ax.plot(vel, Jout, lw=1.5, label='Jout (escaped)')
    if Jin is not None:
        ax.plot(vel, Jin,  lw=1.0, ls='--', label='Jin (total interior)')
    ax.set_xlabel('Velocity  [km/s]')
    ax.set_ylabel('J  [arbitrary]')
    ax.set_title('Spectrum')
    ax.legend(fontsize=8)

    # --- Peel-off image (collapsed along frequency) ---
    ax = axes[1]
    if obs_file is not None:
        cube = load_peeloff(obs_file)
        if cube.ndim == 3:
            #im = cube.sum(axis=0)   # collapse frequency axis
            im = cube.sum(axis=2)   # collapse frequency axis
        else:
            im = cube
        masked = np.ma.masked_where(im == 0, im)
        cmap = cm.viridis.copy()
        cmap.set_bad('white')
        img = ax.imshow(masked, origin='lower', cmap=cmap)
        plt.colorbar(img, ax=ax, fraction=0.046)
        ax.set_title('Peel-off image (sum over frequency)')
    else:
        ax.text(0.5, 0.5, 'No peel-off file', ha='center', va='center',
                transform=ax.transAxes)
        ax.set_title('Peel-off image')
    ax.set_xlabel('x pixel');  ax.set_ylabel('y pixel')

    # --- Clump positions (z-slice |z| < 0.1) ---
    if clump_file is not None and len(axes) > 2:
        ax = axes[2]
        x, y, z, vx, vy, vz = load_clumps(clump_file)
        # Case-insensitive header access (works for both FITS and HDF5).
        lf_cl = load_lart(clump_file)
        cl_sec = lf_cl.section('Clumps') or (lf_cl.sections[0] if lf_cl.sections else None)
        R      = float(cl_sec.attr('SPHERE_R', 1.0))
        r_cl   = float(cl_sec.attr('CL_RAD',   0.01))
        n_cl   = int  (cl_sec.attr('N_CLUMPS', len(x)))
        f_vol  = float(cl_sec.attr('F_VOL',    0.0))
        n_cov  = float(cl_sec.attr('N_COV',    0.0))
        mask = np.abs(z) < 0.1
        ax.scatter(x[mask], y[mask], s=max(1, int(300 * r_cl)), alpha=0.5)
        ax.set_xlim(-R * 1.05, R * 1.05)
        ax.set_ylim(-R * 1.05, R * 1.05)
        ax.set_aspect('equal')
        ax.set_xlabel('x  [code units]');  ax.set_ylabel('y  [code units]')
        ax.set_title(f'Clump positions (|z|<0.1)\n'
                     f'N={n_cl}, f_vol={f_vol:.3f}, N_cov={n_cov:.2f}')

    plt.tight_layout()
    if save_pdf:
       out_pdf = base + '_plot.pdf'
       plt.savefig(out_pdf, dpi=150)
       print(f'  saved  {out_pdf}')
       plt.close()


# -----------------------------------------------------------------------
# Comparison plot: spectra from multiple runs
# -----------------------------------------------------------------------

def plot_spectra_comparison(bases, title='Clumpy sphere spectra', save_pdf=False):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.set_title(title)
    for base in bases:
        fname = find_lart_file(base)
        if fname is None:
            continue
        vel, Jout, _ = load_spectrum(fname)
        ax.plot(vel, Jout, lw=1.5, label=os.path.basename(base))
    ax.set_xlabel('Velocity  [km/s]')
    ax.set_ylabel('J  [arbitrary]')
    ax.legend(fontsize=8)
    plt.tight_layout()
    if save_pdf:
       out_pdf = 'clump_spectra_comparison.pdf'
       plt.savefig(out_pdf, dpi=150)
       print(f'  saved  {out_pdf}')
       plt.close()


# -----------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------

DEFAULT_BASES = [
    'clump_tau3_fcov5',
    'clump_tau4_fcov5',
    'clump_tau3_fcov20',
    'clump_tau3_fcov5_hubble',
    'clump_tau3_fcov5_sigv',
    'clump_tau3_fcov5_peel',
]

if __name__ == '__main__':
    bases = sys.argv[1:] if len(sys.argv) > 1 else DEFAULT_BASES

    print('Individual run plots:')
    for base in bases:
        print(f'  processing {base} ...')
        plot_run(base)

    print('Comparison spectrum plot:')
    plot_spectra_comparison(bases)
