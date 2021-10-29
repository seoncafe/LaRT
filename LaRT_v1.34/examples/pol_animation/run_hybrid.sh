#!/bin/bash
exec < /dev/null 2>&1
trap "" HUP

#EXEC=../../LaRT_FS.x
EXEC=../../LaRT.x
HOSTS=mocafe,lart1,lart2,lart3

mpirun -hosts $HOSTS -ppn 1 $EXEC t1tau3.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t1tau3_cub111.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t1tau3_cub111b.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t1tau3_cub111c.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t1tau3_cub111_Vexp20.in
