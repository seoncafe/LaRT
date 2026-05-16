#!/usr/bin/env python3
"""Generate input .in files for the HeI 10833 coherent-vs-incoherent test.

Grid:
    - source geometry: central point ('point') vs filled sphere ('uniform_sphere')
    - taumax: 1, 10, 100, 1000
    - HeI_coherent: .false. vs .true.

Total: 4 tau x 2 sources x 2 coherence = 16 cases.

Output file name convention:
    <source>_tau<tau>_<inc|coh>.in
where source = 'pt' (point) | 'un' (uniform).
"""
from __future__ import annotations
import os

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

TEMPLATE = """\
&parameters
 par%line_id      = 'HeI_10833'
 par%HeI_coherent = {coh}
 par%no_photons   = {nphot:.1e}
 par%temperature  = 1.0e4
 par%taumax       = {tau:.4e}
 par%DGR             = 0.0
 par%comoving_source = .false.
 par%save_all        = .false.
 par%recoil          = .false.
 par%use_stokes      = .true.
 par%geometry        = 'sphere'
 par%source_geometry = '{src_geom}'
 par%source_rmax     = 1.0
 par%spectral_type   = 'voigt'
 par%nx               = 101
 par%ny               = 101
 par%nz               = 101
 par%rmax             = 1.0
 par%xmax             = 1.0
 par%ymax             = 1.0
 par%zmax             = 1.0
 par%nvelocity        = 201
 par%velocity_min     = {vmin:.1f}
 par%velocity_max     = {vmax:.1f}
 par%save_peeloff     = .true.
 par%save_peeloff_3D  = .true.
 par%nxim             = 101
 par%nyim             = 101
 par%distance         = 100.0
 par%nobs             = 1
 par%alpha(1)         = 0.0
 par%beta(1)          = 0.0
 par%nprint           = 1.0e7
/
"""

# tau-dependent velocity window (km/s) - wider for higher tau
VWINDOW = {
    0.1:  (-50.0,  30.0),
    1:    (-50.0,  30.0),
    10:   (-80.0,  40.0),
    100:  (-120.0, 60.0),
    1000: (-200.0, 100.0),
}

# tau-dependent photon count - higher tau needs more photons for clean wings
NPHOT = {
    0.1:  5.0e6,
    1:    5.0e6,
    10:   5.0e6,
    100:  5.0e6,
    1000: 5.0e6,
}

# source geometry mapping: tag -> LaRT keyword
SOURCES = {
    'pt': 'point',
    'un': 'uniform_sphere',
}


def main() -> None:
    written = []
    for src_tag, src_geom in SOURCES.items():
        for tau, (vmin, vmax) in VWINDOW.items():
            for coh_tag, coh_val in [('inc', '.false.'), ('coh', '.true.')]:
                fname = f'{src_tag}_tau{tau}_{coh_tag}.in'
                fpath = os.path.join(THIS_DIR, fname)
                content = TEMPLATE.format(
                    coh=coh_val,
                    nphot=NPHOT[tau],
                    tau=float(tau),
                    src_geom=src_geom,
                    vmin=vmin,
                    vmax=vmax,
                )
                with open(fpath, 'w') as f:
                    f.write(content)
                written.append(fname)
    print(f'Wrote {len(written)} input files:')
    for f in written:
        print(' ', f)


if __name__ == '__main__':
    main()
