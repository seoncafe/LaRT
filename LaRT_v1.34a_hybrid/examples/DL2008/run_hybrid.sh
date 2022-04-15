#!/bin/bash
exec < /dev/null 2>&1
trap "" HUP

#EXEC=../../LaRT_FS.x
EXEC=../../LaRT.x
HOSTS=mocafe,lart1,lart2,lart3

#----- These are the best parameters for DL2008 model.
mpirun -hosts $HOSTS -ppn 1 $EXEC DL19e.in
mpirun -hosts $HOSTS -ppn 1 $EXEC DL20e.in

#mpirun -hosts $HOSTS -ppn 1 $EXEC DL19e_dust.in
#mpirun -hosts $HOSTS -ppn 1 $EXEC DL20e_dust.in
