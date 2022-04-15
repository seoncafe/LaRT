#!/bin/bash
exec < /dev/null 2>&1
trap "" HUP

#EXEC=../../LaRT_FS.x
#EXEC=../../LaRT_calcP.x
#EXEC=../../LaRT_calcJP.x
EXEC=../../LaRT.x

HOSTS=mocafe,lart1,lart2,lart3

mpirun -hosts $HOSTS -ppn 1 $HOST $EXEC t1tau3.in
mpirun -hosts $HOSTS -ppn 1 $HOST $EXEC t1tau4.in
mpirun -hosts $HOSTS -ppn 1 $HOST $EXEC t1tau5.in
mpirun -hosts $HOSTS -ppn 1 $HOST $EXEC t1tau6.in
mpirun -hosts $HOSTS -ppn 1 $HOST $EXEC t1tau7.in

mpirun -hosts $HOSTS -ppn 1 $EXEC t4tau3.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4tau4.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4tau5.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4tau6.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4tau7.in

#mpirun -hosts $HOSTS -ppn 1 $EXEC t4tau6dust.in
#mpirun -hosts $HOSTS -ppn 1 $EXEC t1tau7dust.in
#mpirun -hosts $HOSTS -ppn 1 $EXEC t4tau7dust.in
