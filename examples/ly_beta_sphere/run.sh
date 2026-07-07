#!/bin/bash
# Ly-beta fluorescence smoke tests (Phase 1a: Cartesian, outside observer).
# Static uniform sphere, T = 1e4 K, taumax = 1e4, central point source,
# monochromatic injection at line center.
#   t4tau1e4.in       : dust-free   (band2_esc == W_conv exactly)
#   t4tau1e4_dust.in  : DGR = 100   (visible dust absorption in both bands)
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NP=${NP:-8}
mpirun -np $NP "$BASE"/../../LaRT.x "$BASE"/t4tau1e4.in
mpirun -np $NP "$BASE"/../../LaRT.x "$BASE"/t4tau1e4_dust.in
