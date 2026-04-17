#!/bin/bash
exec < /dev/null 2>&1
trap "" HUP

EXEC=../../LaRT.x

NTHREADS=72

mpirun -np $NTHREADS $EXEC t1e5tau1e0.in
mpirun -np $NTHREADS $EXEC t1e5tau1e0_V050.in
mpirun -np $NTHREADS $EXEC t1e5tau1e0_V100.in

mpirun -np $NTHREADS $EXEC t1e5tau2e0.in
mpirun -np $NTHREADS $EXEC t1e5tau2e0_V050.in
mpirun -np $NTHREADS $EXEC t1e5tau2e0_V100.in

mpirun -np $NTHREADS $EXEC t1e5tau5e0.in
mpirun -np $NTHREADS $EXEC t1e5tau5e0_V050.in
mpirun -np $NTHREADS $EXEC t1e5tau5e0_V100.in

mpirun -np $NTHREADS $EXEC t1e5tau1e1.in
mpirun -np $NTHREADS $EXEC t1e5tau1e1_V050.in
mpirun -np $NTHREADS $EXEC t1e5tau1e1_V100.in

mpirun -np $NTHREADS $EXEC t1e5tau2e1.in
mpirun -np $NTHREADS $EXEC t1e5tau2e1_V050.in
mpirun -np $NTHREADS $EXEC t1e5tau2e1_V100.in

#mpirun -machinefile $NTHREADS $EXEC t1e5tau1e0.in
#mpirun -machinefile $NTHREADS $EXEC t1e5tau1e0_V050.in
#mpirun -machinefile $NTHREADS $EXEC t1e5tau1e0_V100.in
#
#mpirun -machinefile $NTHREADS $EXEC t1e5tau2e0.in
#mpirun -machinefile $NTHREADS $EXEC t1e5tau2e0_V050.in
#mpirun -machinefile $NTHREADS $EXEC t1e5tau2e0_V100.in
#
#mpirun -machinefile $NTHREADS $EXEC t1e5tau5e0.in
#mpirun -machinefile $NTHREADS $EXEC t1e5tau5e0_V050.in
#mpirun -machinefile $NTHREADS $EXEC t1e5tau5e0_V100.in
#
#mpirun -machinefile $NTHREADS $EXEC t1e5tau1e1.in
#mpirun -machinefile $NTHREADS $EXEC t1e5tau1e1_V050.in
#mpirun -machinefile $NTHREADS $EXEC t1e5tau1e1_V100.in
#
#mpirun -machinefile $NTHREADS $EXEC t1e5tau2e1.in
#mpirun -machinefile $NTHREADS $EXEC t1e5tau2e1_V050.in
#mpirun -machinefile $NTHREADS $EXEC t1e5tau2e1_V100.in
