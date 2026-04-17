#!/bin/bash
exec < /dev/null 2>&1
trap "" HUP

EXEC=../../LaRT.x

HOSTS=mocafe
#HOSTS=mocafe,lart2
#HOSTS=mocafe,lart1,lart2,lart3

#====== Do not touch starting from here =========
host_file=/tmp/host_file_$RANDOM

array=$(echo $HOSTS | tr "," "\n")
for host in $array
do
   if [[ $host = "mocafe" ]]
   then
      num=88
   else
      num=72
   fi
   echo $host:$num >> $host_file
done
#====== Do not touch up to here =========

#----- These are the best parameters for DL2008 model.
mpirun -machinefile $host_file $EXEC DL19e.in
mpirun -machinefile $host_file $EXEC DL20e.in

#mpirun -machinefile $host_file $EXEC DL19e_dust.in
#mpirun -machinefile $host_file $EXEC DL20e_dust.in
