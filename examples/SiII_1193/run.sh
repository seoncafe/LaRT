#!/bin/bash
exec < /dev/null 2>&1
trap "" HUP

EXEC=../LaRT.x
HOST=all_hosts

mpirun -machinefile $HOST $EXEC tau1e+0_V000.in
mpirun -machinefile $HOST $EXEC tau1e+0_V050.in
mpirun -machinefile $HOST $EXEC tau1e+0_V100.in
mpirun -machinefile $HOST $EXEC tau1e+0_V200.in
mpirun -machinefile $HOST $EXEC tau1e+0_V400.in
mpirun -machinefile $HOST $EXEC tau2e+0_V000.in
mpirun -machinefile $HOST $EXEC tau2e+0_V050.in
mpirun -machinefile $HOST $EXEC tau2e+0_V100.in
mpirun -machinefile $HOST $EXEC tau2e+0_V200.in
mpirun -machinefile $HOST $EXEC tau2e+0_V400.in
mpirun -machinefile $HOST $EXEC tau5e+0_V000.in
mpirun -machinefile $HOST $EXEC tau5e+0_V050.in
mpirun -machinefile $HOST $EXEC tau5e+0_V100.in
mpirun -machinefile $HOST $EXEC tau5e+0_V200.in
mpirun -machinefile $HOST $EXEC tau5e+0_V400.in
mpirun -machinefile $HOST $EXEC tau1e+1_V000.in
mpirun -machinefile $HOST $EXEC tau1e+1_V050.in
mpirun -machinefile $HOST $EXEC tau1e+1_V100.in
mpirun -machinefile $HOST $EXEC tau1e+1_V200.in
mpirun -machinefile $HOST $EXEC tau1e+1_V400.in
mpirun -machinefile $HOST $EXEC tau2e+1_V000.in
mpirun -machinefile $HOST $EXEC tau2e+1_V050.in
mpirun -machinefile $HOST $EXEC tau2e+1_V100.in
mpirun -machinefile $HOST $EXEC tau2e+1_V200.in
mpirun -machinefile $HOST $EXEC tau2e+1_V400.in
mpirun -machinefile $HOST $EXEC tau5e+1_V000.in
mpirun -machinefile $HOST $EXEC tau5e+1_V050.in
mpirun -machinefile $HOST $EXEC tau5e+1_V100.in
mpirun -machinefile $HOST $EXEC tau5e+1_V200.in
mpirun -machinefile $HOST $EXEC tau5e+1_V400.in
mpirun -machinefile $HOST $EXEC tau1e+2_V000.in
mpirun -machinefile $HOST $EXEC tau1e+2_V050.in
mpirun -machinefile $HOST $EXEC tau1e+2_V100.in
mpirun -machinefile $HOST $EXEC tau1e+2_V200.in
mpirun -machinefile $HOST $EXEC tau1e+2_V400.in
mpirun -machinefile $HOST $EXEC tau2e+2_V000.in
mpirun -machinefile $HOST $EXEC tau2e+2_V050.in
mpirun -machinefile $HOST $EXEC tau2e+2_V100.in
mpirun -machinefile $HOST $EXEC tau2e+2_V200.in
mpirun -machinefile $HOST $EXEC tau2e+2_V400.in
mpirun -machinefile $HOST $EXEC tau5e+2_V000.in
mpirun -machinefile $HOST $EXEC tau5e+2_V050.in
mpirun -machinefile $HOST $EXEC tau5e+2_V100.in
mpirun -machinefile $HOST $EXEC tau5e+2_V200.in
mpirun -machinefile $HOST $EXEC tau5e+2_V400.in
mpirun -machinefile $HOST $EXEC tau1e+3_V000.in
mpirun -machinefile $HOST $EXEC tau1e+3_V050.in
mpirun -machinefile $HOST $EXEC tau1e+3_V100.in
mpirun -machinefile $HOST $EXEC tau1e+3_V200.in
mpirun -machinefile $HOST $EXEC tau1e+3_V400.in
