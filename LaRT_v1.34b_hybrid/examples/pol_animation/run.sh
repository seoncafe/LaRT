#!/bin/bash
exec < /dev/null 2>&1
trap "" HUP

#EXEC=../../LaRT_FS.x
EXEC=../../LaRT.x
HOST=all_hosts

mpirun -machinefile $HOST $EXEC t1tau3.in
mpirun -machinefile $HOST $EXEC t1tau3_cub111.in
mpirun -machinefile $HOST $EXEC t1tau3_cub111b.in
mpirun -machinefile $HOST $EXEC t1tau3_cub111c.in
mpirun -machinefile $HOST $EXEC t1tau3_cub111_Vexp20.in
