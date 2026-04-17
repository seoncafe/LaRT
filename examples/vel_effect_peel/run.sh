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

echo ""
echo "Running $EXEC on $HOSTS with $NTHREADS threads."
echo "   with the machinefile $host_file"
#====== Do not touch up to here =========

mpirun -machinefile $host_file $EXEC t4NHI2_18_V0020.in
mpirun -machinefile $host_file $EXEC t4NHI2_18_V0200.in
mpirun -machinefile $host_file $EXEC t4NHI2_18_V2000.in

mpirun -machinefile $host_file $EXEC t4NHI2_19_V0020.in
mpirun -machinefile $host_file $EXEC t4NHI2_19_V0100.in
mpirun -machinefile $host_file $EXEC t4NHI2_19_V0200.in
mpirun -machinefile $host_file $EXEC t4NHI2_19_V0500.in
mpirun -machinefile $host_file $EXEC t4NHI2_19_V1000.in
mpirun -machinefile $host_file $EXEC t4NHI2_19_V1500.in
mpirun -machinefile $host_file $EXEC t4NHI2_19_V2000.in
mpirun -machinefile $host_file $EXEC t4NHI2_19_V3000.in

mpirun -machinefile $host_file $EXEC t4NHI2_20_V0020.in
mpirun -machinefile $host_file $EXEC t4NHI2_20_V0100.in
mpirun -machinefile $host_file $EXEC t4NHI2_20_V0200.in
mpirun -machinefile $host_file $EXEC t4NHI2_20_V0500.in
mpirun -machinefile $host_file $EXEC t4NHI2_20_V1000.in
mpirun -machinefile $host_file $EXEC t4NHI2_20_V1500.in
mpirun -machinefile $host_file $EXEC t4NHI2_20_V2000.in
mpirun -machinefile $host_file $EXEC t4NHI2_20_V3000.in

#mpirun -machinefile $host_file $EXEC t4NHI2_20_V0020n.in
