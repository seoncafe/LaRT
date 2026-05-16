# He I 10833 Coherent vs. Incoherent Scattering Test

This directory contains a test suite that compares two treatments of the
He I 10833 Å triplet resonance scattering in LaRT:

| Toggle                     | $E_1, E_2, E_3$ source |
|----------------------------|------------------------|
| `par%HeI_coherent = .false.` (default, legacy) | Per-component incoherent values $(E_1, E_2, E_3) = (0, 1, 0)$, $(\tfrac14, \tfrac34, \tfrac14)$, $(\tfrac{7}{20}, \tfrac{13}{20}, \tfrac{3}{4})$ for $P_0, P_1, P_2$ (Seon `scatter_matrix_v1.1` §9.4) — assigned channel-by-channel based on which upper state was sampled. |
| `par%HeI_coherent = .true.`  | Frequency-dependent $E_1(\nu), E_2(\nu), E_3(\nu)$ from the Real-$\Phi$ polynomial form of Seon `scatter_matrix_HeI` memo, eqs (29) & (30) bottom line. The upper-state selection probability is unchanged because the total cross-section is the same incoherent $1{:}3{:}5$ Lorentzian sum (interference cross-terms cancel in $W_{11}+2W_{14}$). |

## Test grid

* Geometry: `sphere`, $T = 10^4$ K
* Source: `point` (central) and `uniform_sphere` (filled)
* Optical depth (pole-to-edge): $\tau = 0.1, 1, 10, 100, 1000$
* Both `HeI_coherent = .false.` and `.true.`

Total: **20 runs**.

## Reproducing the test

```bash
# 1) Generate all 20 input files
python3 generate_inputs.py

# 2) Run them (defaults to 90% of available cores).
#    Each run writes <basename>.h5 (spectrum), <basename>_obs.h5 (peel-off cube),
#    <basename>_stokes.h5 (I,Q,U,V), and a log under logs/.
bash run_all.sh

# Override MPI ranks if desired:
NP=32 bash run_all.sh

# 3) Open the notebook for comparison plots
jupyter notebook compare.ipynb
```

The notebook produces:

1. Angle-integrated escape spectra ($J_{\rm out}$ vs. velocity) — coherent vs. incoherent overlaid, for each $(\text{source}, \tau)$.
2. Fractional difference $(J_{\rm coh} - J_{\rm inc}) / \bar J$.
3. Frequency-integrated 2-D peel-off surface brightness images.
4. Azimuthally averaged radial surface-brightness profiles $\langle I(b)\rangle$ vs. impact parameter $b$.
5. Linear polarization fraction profiles $p(b) = \sqrt{Q^2 + U^2} / I$.
6. The $E_1(\nu), E_2(\nu), E_3(\nu)$ curves themselves (Real-$\Phi$ polynomial form), with vertical lines at $P_0, P_1, P_2$ line centers.

## Expected differences

Because the upper-state selection probabilities are unchanged, the photon-by-photon
re-emission frequency redistribution is statistically identical between the two
treatments. The differences come entirely from the **phase function**:

* `E_1` controls the linear (Rayleigh) anisotropy. Changes in `E_1(\nu)` redistribute
  photons in $\cos\theta$ at each scattering, which propagates into the surface-brightness
  profile and the integrated escape spectrum at high $\tau$.
* `E_3` controls the circular polarizability — visible only in the Stokes-V signature.

For low $\tau$ ($\tau \lesssim 10$) the two treatments should give nearly identical
escape spectra (few scatterings) but already-different linear polarization profiles.
At $\tau = 1000$ the cumulative phase-function difference may produce measurable
shifts in the line profile and a noticeably different $p(b)$.

## File layout

```
HeI_coherent_test/
├── README.md                 # this file
├── generate_inputs.py        # generates the 20 .in files
├── run_all.sh                # batch runner (mpirun all 20)
├── compare.ipynb             # comparison plots
├── pt_tau*_inc.in            # central point source, incoherent
├── pt_tau*_coh.in            # central point source, coherent
├── un_tau*_inc.in            # uniform sphere source, incoherent
├── un_tau*_coh.in            # uniform sphere source, coherent
├── *.h5                      # outputs (spectrum)
├── *_obs.h5                  # outputs (peel-off cube)
├── *_stokes.h5               # outputs (Stokes I,Q,U,V cube)
└── logs/                     # stdout/stderr of each LaRT run
```

## Notes

* The same identical input files (and `run_all.sh` / notebook) work in
  `LaRT_v2.00/examples/HeI_coherent_test/`.
* Output is HDF5 (default in v2.10; equivalent FITS output via
  `par%file_format = 'fits'`).
* Stokes output (`*_stokes.h5`) requires `par%use_stokes = .true.` (set in the
  generated input files).
