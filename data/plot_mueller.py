#!/usr/bin/env python
"""Plot Mueller matrix elements S11, S12, S33, S34 vs cos(theta).

Python port of plot_mueller.pro.
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt


def read_mueller(fname):
    """Read a Mueller matrix data file.

    Header layout (3 header lines):
      line 1: column description
      line 2: wavl, Cext, albedo, <cos> (hgg), nang
      line 3: column description for the table
      line 4+: cos(theta), S11, S12, S33, S34
    """
    with open(fname) as f:
        f.readline()
        meta = f.readline().split()
    wavl   = float(meta[0])
    cext   = float(meta[1])
    albedo = float(meta[2])
    hgg    = float(meta[3])
    nang   = int(meta[4])

    cost, s11, s12, s33, s34 = np.loadtxt(fname, skiprows=3, unpack=True)
    return dict(wavl=wavl, cext=cext, albedo=albedo, hgg=hgg, nang=nang,
                cost=cost, s11=s11, s12=s12, s33=s33, s34=s34)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('-f', '--file', default='mueller_Lyalpha.dat',
                   help='Mueller matrix data file')
    p.add_argument('-o', '--out', default=None,
                   help='Output figure (e.g. mueller_Lya.pdf). If omitted, show on screen.')
    args = p.parse_args()

    d = read_mueller(args.file)
    cost = d['cost']

    fig, ax = plt.subplots(2, 2, figsize=(8, 6))
    ax = ax.ravel()

    ax[0].plot(cost, d['s11'])
    ax[0].set_ylabel(r'$S_{11}$')

    ax[1].plot(cost, d['s12'])
    ax[1].set_ylim(-0.5, 0)
    ax[1].set_ylabel(r'$S_{12}$')

    ax[2].plot(cost, d['s33'])
    ax[2].set_ylabel(r'$S_{33}$')

    ax[3].plot(cost, d['s34'])
    ax[3].set_ylabel(r'$S_{34}$')

    for a in ax:
        a.set_xlim(-1, 1)
        a.set_xlabel(r'$\cos\theta$')

    fig.suptitle(r'$\lambda = ' + f'{d["wavl"]:.4f}' + r'\,\mu$m, '
                 + r'$g = ' + f'{d["hgg"]:.4f}' + r'$, '
                 + r'albedo $= ' + f'{d["albedo"]:.4f}' + r'$')
    fig.tight_layout()

    if args.out:
        fig.savefig(args.out)
    else:
        plt.show()


if __name__ == '__main__':
    main()
