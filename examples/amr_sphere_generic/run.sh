#!/bin/bash
# Run the generic-format AMR sphere example.
# Adjust the path to LaRT.x and the number of MPI ranks as needed.
iexec < /dev/null 2>&1
trap "" HUP

LART=../../LaRT.x
NP=70

mpirun -np ${NP} ${LART} sphere_car_tau2_point.in
mpirun -np ${NP} ${LART} sphere_car_tau4_point.in

#mpirun -np ${NP} ${LART} sphere_amr_tau2_point.in
#mpirun -np ${NP} ${LART} sphere_amr_tau4_point.in
