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

mpirun -machinefile $host_file $EXEC t1tau3.in
mpirun -machinefile $host_file $EXEC t1tau3_cub111.in
mpirun -machinefile $host_file $EXEC t1tau3_cub111b.in
mpirun -machinefile $host_file $EXEC t1tau3_cub111c.in
mpirun -machinefile $host_file $EXEC t1tau3_cub111_Vexp20.in
