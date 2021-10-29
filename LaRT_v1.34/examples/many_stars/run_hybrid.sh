#!/bin/bash
exec < /dev/null 2>&1
trap "" HUP

EXEC=../../LaRT.x
HOSTS=mocafe,lart1,lart2,lart3

mpirun -hosts $HOSTS -ppn 1 $EXEC stars1.in
mpirun -hosts $HOSTS -ppn 1 $EXEC stars2.in
