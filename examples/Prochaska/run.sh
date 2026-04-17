#!/bin/bash
exec < /dev/null 2>&1
trap "" HUP

EXEC=../LaRT.x
EXEC_P=../LaRT_PROCHASKA.x
HOST=all_hosts
#HOST=lart12

mpirun -machinefile $HOST $EXEC_P FeII_UV1_Pc.in
mpirun -machinefile $HOST $EXEC_P FeII_UV1_Pb.in
mpirun -machinefile $HOST $EXEC_P FeII_UV1_Pa.in

mpirun -machinefile $HOST $EXEC FeII_UV1_c.in
mpirun -machinefile $HOST $EXEC FeII_UV1_b.in
mpirun -machinefile $HOST $EXEC FeII_UV1_a.in

mpirun -machinefile $HOST $EXEC FeII_UV2_c.in
mpirun -machinefile $HOST $EXEC FeII_UV2_b.in
mpirun -machinefile $HOST $EXEC FeII_UV2_a.in

mpirun -machinefile $HOST $EXEC FeII_UV3_c.in
mpirun -machinefile $HOST $EXEC FeII_UV3_b.in
mpirun -machinefile $HOST $EXEC FeII_UV3_a.in

mpirun -machinefile $HOST $EXEC MgII_b.in
mpirun -machinefile $HOST $EXEC MgII_a.in
