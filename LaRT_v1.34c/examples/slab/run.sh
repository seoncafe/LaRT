#!/bin/bash
exec < /dev/null 2>&1
trap "" HUP

HOST=all_hosts

#EXEC=../../LaRT_FS.x
#EXEC=../../LaRT_calcP.x
EXEC=../../LaRT.x

mpirun -machinefile $HOST $EXEC t1tau4.in
mpirun -machinefile $HOST $EXEC t1tau5.in
mpirun -machinefile $HOST $EXEC t1tau6.in
mpirun -machinefile $HOST $EXEC t1tau7.in

mpirun -machinefile $HOST $EXEC t4tau4.in
mpirun -machinefile $HOST $EXEC t4tau5.in
mpirun -machinefile $HOST $EXEC t4tau6.in
mpirun -machinefile $HOST $EXEC t4tau7.in
