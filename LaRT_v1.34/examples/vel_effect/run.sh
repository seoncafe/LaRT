#!/bin/bash
exec < /dev/null 2>&1
trap "" HUP

#EXEC=../../LaRT_carJ.x
EXEC=../../LaRT.x
HOST=all_hosts

mpirun -machinefile $HOST $EXEC t4NHI2_20_V2000.in
mpirun -machinefile $HOST $EXEC t4NHI2_20_V0200.in
mpirun -machinefile $HOST $EXEC t4NHI2_20_V0020.in
mpirun -machinefile $HOST $EXEC t4NHI2_20_V0000.in
