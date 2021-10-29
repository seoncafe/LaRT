#!/bin/bash
exec < /dev/null 2>&1
trap "" HUP

#EXEC=../../LaRT_peel.x
#EXEC=../../LaRT_peel_calcJP.x
#EXEC=../../LaRT_calcP.x
EXEC=../../LaRT_peel_calcP.x
HOSTS=mocafe,lart1,lart2,lart3

mpirun -hosts $HOSTS -ppn 1 $EXEC clumpy.in
