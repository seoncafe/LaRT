#!/bin/bash
# Generate the AMR data file (once) and run the standalone
# make_sightline_tau.x driver for all three grid modes.
set -e

NP=${NP:-4}
LART=../../make_sightline_tau.x

# 1) AMR data file (skip if already present)
if [ ! -f amr_sphere.fits.gz ]; then
    echo "=== generating amr_sphere.fits.gz ==="
    python3 make_amr_sphere_data.py
fi

# 2) Cartesian
echo "=== Cartesian sightline tau ==="
mpirun -np ${NP} ${LART} car_sightline.in

# 3) AMR
echo "=== AMR sightline tau ==="
mpirun -np ${NP} ${LART} amr_sightline.in

# 4) Clump
echo "=== Clump sightline tau ==="
mpirun -np ${NP} ${LART} clump_sightline.in

# 5) Cartesian (inside observer, HEALPix all-sky)
echo "=== Cartesian sightline tau (inside observer / HEALPix) ==="
mpirun -np ${NP} ${LART} car_inside.in

# 6) AMR (inside observer, HEALPix all-sky)
echo "=== AMR sightline tau (inside observer / HEALPix) ==="
mpirun -np ${NP} ${LART} amr_inside.in
