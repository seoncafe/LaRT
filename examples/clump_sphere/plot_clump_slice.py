#!/usr/bin/env python3
"""
plot_clump_slice.py  --  Plot the boundaries of clumps cut by a coordinate
slice (x = const, y = const, or z = const) using the *_clumps.fits.gz file
written by LaRT.

For each clump intersected by the slice plane, the displayed boundary is the
exact 2D cross-section of the spherical clump:

    boundary radius  r_cross = sqrt(R_CLUMP^2 - (slice - center)^2)

Clumps whose centre lies farther from the slice than their R_CLUMP do not
appear (they don't cross the plane).

Usage
-----
  # default: z = 0 slice, no fill, white background
  python plot_clump_slice.py clump_powerlaw_density_clumps.fits.gz

  # x = 0.2, colour outline by full clump radius
  python plot_clump_slice.py file_clumps.fits.gz --axis x --value 0.2 \
      --colorby radius

  # y = 0, fill the cross-section disks (alpha=0.4) coloured by temperature
  python plot_clump_slice.py file_clumps.fits.gz --axis y --value 0 \
      --colorby temp --fill

  # save to file instead of showing
  python plot_clump_slice.py file_clumps.fits.gz -o slice.png

Color-by options:
  none     -- single colour outline (default tab:blue)
  radius   -- colour by full 3D R_CLUMP (uniform per clump)
  rcross   -- colour by 2D intersection radius (depends on slice)
  temp     -- per-clump temperature
  rhokap   -- per-clump rhokap

The full sphere boundary (par%rmax / SPHERE_R from header) is drawn as a
dashed circle for reference.

Requires: numpy, matplotlib, astropy
"""

import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from astropy.io import fits


AXIS_TRIPLE = {
    'x': (1, 2, 0, 'y', 'z'),   # plotted axes index, slice axis index, labels
    'y': (0, 2, 1, 'x', 'z'),
    'z': (0, 1, 2, 'x', 'y'),
}


def load_clumps(fname):
    """Return dict of clump arrays + header metadata."""
    with fits.open(fname) as hdul:
        tbl = hdul[1].data
        hdr = hdul[1].header
        pos = np.column_stack([tbl['X'], tbl['Y'], tbl['Z']])
        out = {
            'pos':    pos,
            'radius': tbl['R_CLUMP'] if 'R_CLUMP' in tbl.names else
                      np.full(len(pos), hdr.get('CL_RAD', 0.01)),
            'rhokap': tbl['RHOKAP']  if 'RHOKAP'  in tbl.names else None,
            'temp':   tbl['TEMP']    if 'TEMP'    in tbl.names else None,
            'sphere_R': hdr.get('SPHERE_R', hdr.get('RMAX', 1.0)),
            'n_clumps': hdr.get('N_CLUMPS', len(pos)),
            'f_vol':   hdr.get('F_VOL',  None),
            'f_cov':   hdr.get('F_COV',  None),
        }
    return out


def slice_clumps(pos, radius, axis, value):
    """Return (a, b, rcross, idx) for clumps that cross the slice plane.

    a, b -- coordinates in the 2D plot plane
    rcross -- 2D circle radius at the cut
    idx -- indices into the original arrays (for colour-by-property)
    """
    ia, ib, ic, _, _ = AXIS_TRIPLE[axis]
    delta = value - pos[:, ic]
    mask  = np.abs(delta) < radius
    rcross = np.sqrt(radius[mask] ** 2 - delta[mask] ** 2)
    return pos[mask, ia], pos[mask, ib], rcross, np.where(mask)[0]


def make_color_array(clumps, idx, key, rcross):
    """Return a 1D array of values to colour clumps by, plus a label."""
    if key == 'none':
        return None, None
    if key == 'radius':
        return clumps['radius'][idx], 'R_clump  [code]'
    if key == 'rcross':
        return rcross, 'r_cross  [code]'
    if key == 'temp':
        if clumps['temp'] is None:
            sys.exit("Error: TEMP column not found in clumps file")
        return clumps['temp'][idx], 'Temperature  [K]'
    if key == 'rhokap':
        if clumps['rhokap'] is None:
            sys.exit("Error: RHOKAP column not found in clumps file")
        return clumps['rhokap'][idx], 'RHOKAP  [code]'
    sys.exit(f"Unknown --colorby option: {key}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('clumps_file',
                    help='*_clumps.fits.gz produced by LaRT')
    ap.add_argument('--axis', choices=('x', 'y', 'z'), default='z',
                    help="slice axis (default: z)")
    ap.add_argument('--value', type=float, default=0.0,
                    help="slice coordinate (default: 0.0, i.e. through origin)")
    ap.add_argument('--colorby',
                    choices=('none', 'radius', 'rcross', 'temp', 'rhokap'),
                    default='none',
                    help='colour outlines by per-clump quantity (default: none)')
    ap.add_argument('--fill', action='store_true',
                    help='fill cross-section disks (alpha=0.4) instead of outline only')
    ap.add_argument('--cmap', default='viridis',
                    help='matplotlib colormap when --colorby is set')
    ap.add_argument('--linewidth', type=float, default=0.7,
                    help='outline line width')
    ap.add_argument('--alpha', type=float, default=None,
                    help='alpha for filled disks (default 0.4 with --fill, 1.0 outline)')
    ap.add_argument('--no-sphere', action='store_true',
                    help="do not draw the outer sphere boundary")
    ap.add_argument('-o', '--output',
                    help='save figure to this path instead of showing')
    ap.add_argument('--dpi', type=int, default=150)
    ap.add_argument('--figsize', type=float, nargs=2, default=(7.0, 7.0))
    args = ap.parse_args()

    if not os.path.exists(args.clumps_file):
        sys.exit(f"File not found: {args.clumps_file}")

    clumps = load_clumps(args.clumps_file)
    a, b, rcross, idx = slice_clumps(clumps['pos'], clumps['radius'],
                                     args.axis, args.value)

    fig, ax = plt.subplots(figsize=tuple(args.figsize))

    # outer sphere boundary at this slice
    R_sphere = clumps['sphere_R']
    if not args.no_sphere and abs(args.value) < R_sphere:
        r_outer = np.sqrt(R_sphere ** 2 - args.value ** 2)
        ax.add_patch(Circle((0.0, 0.0), r_outer, fill=False,
                            edgecolor='black', linestyle='--', linewidth=1.0,
                            label=f'sphere @ {args.axis}={args.value:g}'))

    color_vals, color_label = make_color_array(clumps, idx,
                                               args.colorby, rcross)

    patches = [Circle((ai, bi), ri) for ai, bi, ri in zip(a, b, rcross)]

    fill_alpha = (args.alpha if args.alpha is not None
                  else (0.4 if args.fill else 1.0))

    if color_vals is None:
        col = PatchCollection(
            patches,
            facecolors=('tab:blue' if args.fill else 'none'),
            edgecolors='tab:blue',
            linewidths=args.linewidth,
            alpha=(fill_alpha if args.fill else 1.0),
        )
        ax.add_collection(col)
    else:
        col = PatchCollection(
            patches,
            cmap=args.cmap,
            edgecolors=('face' if args.fill else None),
            linewidths=args.linewidth,
            alpha=(fill_alpha if args.fill else 1.0),
        )
        col.set_array(np.asarray(color_vals))
        if not args.fill:
            col.set_facecolor('none')
        ax.add_collection(col)
        cb = fig.colorbar(col, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label(color_label)

    # axis labels
    _, _, _, lab_a, lab_b = AXIS_TRIPLE[args.axis]
    ax.set_xlabel(f'{lab_a}  [code]')
    ax.set_ylabel(f'{lab_b}  [code]')
    ax.set_aspect('equal')

    lim = R_sphere * 1.05
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    n_in_plane = len(rcross)
    title = (f'{os.path.basename(args.clumps_file)}\n'
             f'slice {args.axis} = {args.value:g}  |  '
             f'{n_in_plane} / {len(clumps["pos"])} clumps cross plane')
    if clumps['f_cov'] is not None and clumps['f_vol'] is not None:
        title += (f'\nf_cov={clumps["f_cov"]:.3f}  '
                  f'f_vol={clumps["f_vol"]:.4f}')
    ax.set_title(title, fontsize=10)

    if not args.no_sphere and abs(args.value) < R_sphere:
        ax.legend(loc='lower right', fontsize=8)

    plt.tight_layout()

    if args.output:
        plt.savefig(args.output, dpi=args.dpi)
        print(f'  saved {args.output}')
    else:
        plt.show()


if __name__ == '__main__':
    main()
