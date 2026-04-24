#!/usr/bin/env python3
"""
plot_clump.py  --  Quick analysis plots for the LaRT_v2.00 clumpy-sphere examples.

Usage
-----
  python plot_clump.py                  # plots all available output files
  python plot_clump.py clump_tau3_ncov5 # plots a single run by base name

Requires: astropy, matplotlib, numpy
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy.io import fits

# -----------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------

def load_spectrum(fname):
    """Return (velocity, Jout, Jin) arrays from a LaRT FITS spectrum file."""
    with fits.open(fname) as hdul:
        tbl = hdul[1].data
        vel  = tbl['velocity']
        Jout = tbl['Jout']
        Jin  = tbl.field('Jin') if 'Jin' in tbl.names else None
    return vel, Jout, Jin


def load_clumps(fname):
    """Return (x, y, z, vx, vy, vz) arrays from a *_clumps.fits.gz file."""
    with fits.open(fname) as hdul:
        tbl = hdul[1].data
        x  = tbl['X'];   y  = tbl['Y'];   z  = tbl['Z']
        vx = tbl['VX'];  vy = tbl['VY'];  vz = tbl['VZ']
    return x, y, z, vx, vy, vz


def load_peeloff(fname):
    """Return peel-off image cube (nyim, nxim, nxfreq) from _obs.fits.gz."""
    with fits.open(fname) as hdul:
        return hdul[0].data   # shape: (nyim, nxim, nxfreq) or (nyim, nxim)


# -----------------------------------------------------------------------
# Per-run plots
# -----------------------------------------------------------------------

def plot_run(base, save_pdf=False):
    spec_file   = base + '.fits.gz'
    clump_file  = base + '_clumps.fits.gz'
    obs_file    = base + '_obs.fits.gz'

    if not os.path.exists(spec_file):
        print(f'  [skip] {spec_file} not found')
        return

    fig, axes = plt.subplots(1, 3 if os.path.exists(clump_file) else 2,
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
    if os.path.exists(obs_file):
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
    if os.path.exists(clump_file) and len(axes) > 2:
        ax = axes[2]
        x, y, z, vx, vy, vz = load_clumps(clump_file)
        with fits.open(clump_file) as hdul:
            hdr    = hdul[1].header
            R      = hdr.get('SPHERE_R', 1.0)
            r_cl   = hdr.get('CL_RAD',   0.01)
            n_cl   = hdr.get('N_CLUMPS', len(x))
            f_vol  = hdr.get('F_VOL',    0.0)
            n_cov  = hdr.get('N_COV',    0.0)
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
        fname = base + '.fits.gz'
        if not os.path.exists(fname):
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
    'clump_tau3_ncov5',
    'clump_tau4_ncov5',
    'clump_tau3_ncov20',
    'clump_tau3_ncov5_hubble',
    'clump_tau3_ncov5_sigv',
    'clump_tau3_ncov5_peel',
]

if __name__ == '__main__':
    bases = sys.argv[1:] if len(sys.argv) > 1 else DEFAULT_BASES

    print('Individual run plots:')
    for base in bases:
        print(f'  processing {base} ...')
        plot_run(base)

    print('Comparison spectrum plot:')
    plot_spectra_comparison(bases)
