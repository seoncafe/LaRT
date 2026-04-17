#!/bin/bash
# Run the generic-format AMR sphere example.
# Adjust the path to LaRT.x and the number of MPI ranks as needed.

LART=../../LaRT.x
NP=4

echo "=== AMR sphere (generic format), tau=1e4 ==="
mpirun -np ${NP} ${LART} sphere_amr_tau4.in

echo "=== AMR sphere (generic format), tau=1e4 + peel-off ==="
mpirun -np ${NP} ${LART} sphere_amr_tau4_peel.in
