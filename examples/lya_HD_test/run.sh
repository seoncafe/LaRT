#!/bin/bash
# Run all H+D Lyα reference cases. Uses 90% of available CPU threads.
# After completion, open inspect_HD.ipynb to inspect the spectra.

NP=72

#mpirun -np "$NP" ../../LaRT.x "sphere_HD_ref.in"
#mpirun -np "$NP" ../../LaRT.x "sphere_HD_test.in"
#mpirun -np "$NP" ../../LaRT.x "sphere_HD_test0.in"
#mpirun -np "$NP" ../../LaRT.x "sphere_HD_test2.in"
mpirun -np "$NP" ../../LaRT.x "sphere_HD_test3.in"
