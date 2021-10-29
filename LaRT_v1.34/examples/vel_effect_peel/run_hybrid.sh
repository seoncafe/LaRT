#!/bin/bash
exec < /dev/null 2>&1
trap "" HUP

#EXEC=../../LaRT_FS.x
EXEC=../../LaRT.x

HOSTS=mocafe,lart1,lart2,lart3

mpirun -hosts $HOSTS -ppn 1 $EXEC t4NHI2_18_V0020.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4NHI2_18_V0200.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4NHI2_18_V2000.in

mpirun -hosts $HOSTS -ppn 1 $EXEC t4NHI2_19_V0020.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4NHI2_19_V0100.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4NHI2_19_V0200.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4NHI2_19_V0500.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4NHI2_19_V1000.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4NHI2_19_V1500.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4NHI2_19_V2000.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4NHI2_19_V3000.in


mpirun -hosts $HOSTS -ppn 1 $EXEC t4NHI2_20_V0020.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4NHI2_20_V0100.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4NHI2_20_V0200.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4NHI2_20_V0500.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4NHI2_20_V1000.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4NHI2_20_V1500.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4NHI2_20_V2000.in
mpirun -hosts $HOSTS -ppn 1 $EXEC t4NHI2_20_V3000.in

#mpirun -hosts $HOSTS -ppn 1 $EXEC t4NHI2_20_V0020n.in
