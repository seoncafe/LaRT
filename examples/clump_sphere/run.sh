#!/bin/bash
# Run the clumpy-sphere examples for LaRT_v2.00 (use_clump_medium = .true.).
# Adjust NP (MPI rank count) and LART path as needed.

LART=../../LaRT.x
NP=4

echo "=== Clumpy sphere: tau0=1e3, N_cov=5, static ==="
mpirun -np ${NP} ${LART} clump_tau3_ncov5.in

echo "=== Clumpy sphere: tau0=1e4, N_cov=5, static ==="
mpirun -np ${NP} ${LART} clump_tau4_ncov5.in

echo "=== Clumpy sphere: tau0=1e3, N_cov=20, static (dense covering) ==="
mpirun -np ${NP} ${LART} clump_tau3_ncov20.in

echo "=== Clumpy sphere: tau0=1e3, N_cov=5, Hubble expansion Vexp=200 km/s ==="
mpirun -np ${NP} ${LART} clump_tau3_ncov5_hubble.in

echo "=== Clumpy sphere: tau0=1e3, N_cov=5, clump sigma_v=50 km/s ==="
mpirun -np ${NP} ${LART} clump_tau3_ncov5_sigv.in

echo "=== Clumpy sphere: tau0=1e3, N_cov=5, peel-off observer ==="
mpirun -np ${NP} ${LART} clump_tau3_ncov5_peel.in
