#!/bin/bash
exec < /dev/null 2>&1
trap "" HUP

#EXEC=../../LaRT_FS.x
EXEC=../../LaRT.x

HOSTS=mocafe,lart1,lart2,lart3

mpirun -hosts $HOSTS -ppn 1 $EXEC t1tau4.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t1tau5.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t1tau6.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t1tau7.in

mpirun -hosts $HOSTS -ppn 1 $EXEC t4tau4.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4tau5.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4tau6.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4tau7.in

mpirun -hosts $HOSTS -ppn 1 $EXEC t2tau3.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t2tau4.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t2tau5.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t2tau6.in

mpirun -hosts $HOSTS -ppn 1 $EXEC t3tau3.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t3tau4.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t3tau5.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t3tau6.in
