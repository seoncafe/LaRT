#!/bin/bash
exec < /dev/null 2>&1
trap "" HUP

EXEC=../../LaRT.x
HOSTS=lart4,lart3,lart2,lart1,mocafe

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
echo "Running $EXEC on $HOSTS"
echo "   with the machinefile $host_file"
#====== Do not touch up to here =========

mpirun -machinefile $host_file $EXEC t1tau3.in
mpirun -machinefile $host_file $EXEC t1tau2.in
mpirun -machinefile $host_file $EXEC t1tau1.in
mpirun -machinefile $host_file $EXEC t1tau0.in
mpirun -machinefile $host_file $EXEC t1tau01.in
mpirun -machinefile $host_file $EXEC t1tau001.in
mpirun -machinefile $host_file $EXEC t1tau0001.in

mpirun -machinefile $host_file $EXEC t4tau3.in
mpirun -machinefile $host_file $EXEC t4tau2.in
mpirun -machinefile $host_file $EXEC t4tau1.in
mpirun -machinefile $host_file $EXEC t4tau0.in
mpirun -machinefile $host_file $EXEC t4tau01.in
mpirun -machinefile $host_file $EXEC t4tau001.in
mpirun -machinefile $host_file $EXEC t4tau0001.in
