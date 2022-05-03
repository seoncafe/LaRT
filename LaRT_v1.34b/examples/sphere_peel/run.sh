#!/bin/bash
exec < /dev/null 2>&1
trap "" HUP

#EXEC=../../LaRT_calcP.x
EXEC=../../LaRT.x

HOST=all_hosts

mpirun -machinefile $HOST $EXEC t1tau3.in
mpirun -machinefile $HOST $EXEC t1tau4.in
mpirun -machinefile $HOST $EXEC t1tau5.in
mpirun -machinefile $HOST $EXEC t1tau6.in
mpirun -machinefile $HOST $EXEC t1tau7.in

mpirun -machinefile $HOST $EXEC t4tau3.in
mpirun -machinefile $HOST $EXEC t4tau4.in
mpirun -machinefile $HOST $EXEC t4tau5.in
mpirun -machinefile $HOST $EXEC t4tau6.in
mpirun -machinefile $HOST $EXEC t4tau7.in

#mpirun -machinefile $HOST $EXEC t4tau6dust.in
#mpirun -machinefile $HOST $EXEC t1tau7dust.in
#mpirun -machinefile $HOST $EXEC t4tau7dust.in
