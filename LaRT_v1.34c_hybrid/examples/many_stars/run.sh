#!/bin/bash
exec < /dev/null 2>&1
trap "" HUP

EXEC=../../LaRT.x
HOST=all_hosts

#mpirun -machinefile $HOST $EXEC stars1.in
mpirun -machinefile $HOST $EXEC stars2.in
