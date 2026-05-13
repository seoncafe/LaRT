#!/usr/bin/env python
"""Fit Mueller matrix elements S11, S12, S33, S34 with simple analytic models.

Python port of fit_all.pro / fit_s11.pro / fit_s12.pro and the
func_s11.pro, func_s12.pro, func_s33.pro, func_s34.pro model functions.
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# ---------------------------------------------------------------------
# Model functions (faithful ports of func_s11/12/33/34.pro)
# ---------------------------------------------------------------------

def hg(cost, g):
    """Henyey-Greenstein phase function."""
    return 0.5 * (1.0 - g**2) / (1.0 + g**2 - 2.0 * g * cost)**1.5


def func_s11(cost, scale, g1, g2):
    """Two-component HG mixture (func_s11.pro, np=3 branch)."""
    s = abs(scale)
    return s * hg(cost, g1) + (1.0 - s) * hg(cost, g2)


def func_s12(cost, scale):
    """func_s12.pro: scale * (1 - cos^2 theta)."""
    return scale * (1.0 - cost**2)


def func_s33(cost, p0):
    """func_s33.pro: 2 cos / (1 + p0 * sqrt(|cos|))."""
    return 2.0 * cost / (1.0 + p0 * np.sqrt(np.abs(cost)))


def func_s34(cost, pcir, s, se, th0, ex):
    """func_s34.pro: 5-parameter polarization model."""
    theta = np.degrees(np.arccos(cost))
    with np.errstate(invalid='ignore'):
        cc = np.cos(np.radians(theta + s * 3.13 * theta
                               * np.exp(-se * (theta - th0)**ex / 180.0)))
    return pcir * (1.0 - cc**2) / (1.0 + cc**2)


# ---------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------

def read_mueller(fname):
    with open(fname) as f:
        f.readline()
        meta = f.readline().split()
    wavl, _, albedo, hgg = (float(x) for x in meta[:4])
    cost, s11, s12, s33, s34 = np.loadtxt(fname, skiprows=3, unpack=True)
    return wavl, albedo, hgg, cost, s11, s12, s33, s34


# ---------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('-f', '--file', default='mueller_Lyalpha.dat',
                   help='Mueller matrix data file')
    p.add_argument('-o', '--out', default=None,
                   help='Output figure file. If omitted, show on screen.')
    args = p.parse_args()

    wavl, albedo, hgg, cost, s11, s12, s33, s34 = read_mueller(args.file)
    print(f'wavl = {wavl}, albedo = {albedo}, hgg (g) = {hgg}')

    # --- S11
    p11, _ = curve_fit(func_s11, cost, s11, p0=[0.7, 0.8, 0.4])
    print('S11 params (scale, g1, g2):', p11)
    yfit_s11 = func_s11(cost, *p11)

    # --- S12/S11
    p12, _ = curve_fit(func_s12, cost, s12 / s11, p0=[-0.4])
    print('S12 params (scale):', p12)
    yfit_s12 = func_s12(cost, *p12)

    # --- S33/S11
    p33, _ = curve_fit(func_s33, cost, s33 / s11, p0=[1.0])
    print('S33 params (p0):', p33)
    yfit_s33 = func_s33(cost, *p33)

    # --- S34/S11
    p34, _ = curve_fit(func_s34, cost, s34 / s11,
                       p0=[0.39, 1.0, 3.0, 0.0, 1.0])
    print('S34 params (pcir, s, se, th0, ex):', p34)
    yfit_s34 = func_s34(cost, *p34)

    # --- Plot 2x4: top row = ratio (Sij/S11), bottom row = absolute Sij
    fig, ax = plt.subplots(2, 4, figsize=(12, 6))

    # row 0: ratios
    ax[0, 0].plot(cost, s11, 'k', label=r'$S_{11}$')
    ax[0, 0].plot(cost, yfit_s11, 'r', label='fit')
    ax[0, 0].set_yscale('log')
    ax[0, 0].set_ylabel(r'$S_{11}$')
    ax[0, 0].legend()

    ax[0, 1].plot(cost, s12 / s11, 'k')
    ax[0, 1].plot(cost, yfit_s12, 'r')
    ax[0, 1].set_ylabel(r'$S_{12}/S_{11}$')

    ax[0, 2].plot(cost, s33 / s11, 'k')
    ax[0, 2].plot(cost, yfit_s33, 'r')
    ax[0, 2].set_ylabel(r'$S_{33}/S_{11}$')

    ax[0, 3].plot(cost, s34 / s11, 'k')
    ax[0, 3].plot(cost, yfit_s34, 'r')
    ax[0, 3].set_ylabel(r'$S_{34}/S_{11}$')

    # row 1: absolute matrix elements with fit * yfit_s11
    ax[1, 0].plot(cost, s11, 'k')
    ax[1, 0].plot(cost, yfit_s11, 'r')
    ax[1, 0].set_ylabel(r'$S_{11}$')

    ax[1, 1].plot(cost, s12, 'k')
    ax[1, 1].plot(cost, yfit_s12 * yfit_s11, 'r')
    ax[1, 1].set_ylabel(r'$S_{12}$')

    ax[1, 2].plot(cost, s33, 'k')
    ax[1, 2].plot(cost, yfit_s33 * yfit_s11, 'r')
    ax[1, 2].set_ylabel(r'$S_{33}$')

    ax[1, 3].plot(cost, s34, 'k')
    ax[1, 3].plot(cost, yfit_s34 * yfit_s11, 'r')
    ax[1, 3].set_ylabel(r'$S_{34}$')

    for a in ax.ravel():
        a.set_xlim(-1, 1)
        a.set_xlabel(r'$\cos\theta$')

    fig.tight_layout()
    if args.out:
        fig.savefig(args.out)
    else:
        plt.show()


if __name__ == '__main__':
    main()
