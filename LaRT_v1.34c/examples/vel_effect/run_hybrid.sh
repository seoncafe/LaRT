#!/bin/bash
exec < /dev/null 2>&1
trap "" HUP

#EXEC=../../LaRT_carJ.x
EXEC=../../LaRT.x
HOSTS=mocafe,lart1,lart2,lart3

mpirun -hosts $HOSTS -ppn 1 $EXEC t4NHI2_20_V2000.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4NHI2_20_V0200.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4NHI2_20_V0020.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4NHI2_20_V0000.in
