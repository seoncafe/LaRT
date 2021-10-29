#!/bin/bash
exec < /dev/null 2>&1
trap "" HUP

#EXEC=../../LaRT_FS.x
EXEC=../../LaRT.x
HOST=all_hosts

#----- These are the best parameters for DL2008 model.
mpirun -machinefile $HOST $EXEC DL19e.in
mpirun -machinefile $HOST $EXEC DL20e.in

#mpirun -machinefile $HOST $EXEC DL19e_dust.in
#mpirun -machinefile $HOST $EXEC DL20e_dust.in
