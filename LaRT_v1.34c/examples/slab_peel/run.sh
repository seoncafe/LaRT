#!/bin/bash

#EXEC=../../LaRT_FS.x
EXEC=../../LaRT.x
HOST=all_hosts

mpirun -machinefile $HOST $EXEC t1tau4.in
mpirun -machinefile $HOST $EXEC t1tau5.in
mpirun -machinefile $HOST $EXEC t1tau6.in
mpirun -machinefile $HOST $EXEC t1tau7.in

mpirun -machinefile $HOST $EXEC t4tau4.in
mpirun -machinefile $HOST $EXEC t4tau5.in
mpirun -machinefile $HOST $EXEC t4tau6.in
mpirun -machinefile $HOST $EXEC t4tau7.in

mpirun -machinefile $HOST $EXEC t2tau4.in
mpirun -machinefile $HOST $EXEC t2tau5.in
mpirun -machinefile $HOST $EXEC t2tau6.in
mpirun -machinefile $HOST $EXEC t2tau7.in

mpirun -machinefile $HOST $EXEC t3tau4.in
mpirun -machinefile $HOST $EXEC t3tau5.in
mpirun -machinefile $HOST $EXEC t3tau6.in
mpirun -machinefile $HOST $EXEC t3tau7.in

#mpirun -machinefile $HOST $EXEC t1tau4r.in
