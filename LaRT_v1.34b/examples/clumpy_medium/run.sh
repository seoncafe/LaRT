#!/bin/bash
exec < /dev/null 2>&1
trap "" HUP

#EXEC=../../LaRT_calcJP.x
#EXEC=../../LaRT_calcP.x
EXEC=../../LaRT.x

HOST=all_hosts

mpirun -machinefile $HOST $EXEC clumpy.in
