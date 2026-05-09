#!/bin/bash
exec < /dev/null 2>&1
trap "" HUP

EXEC=../../LaRT.x

#HOSTS=mocafe
#HOSTS=mocafe,lart2
#HOSTS=mocafe,lart1,lart2,lart3,lart4
#HOSTS=lart4,lart3,lart2,lart1
HOSTS=lart4
#HOSTS=lart4,lart3,lart2
#HOSTS=lart4,lart3,lart1

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
echo "Running $EXEC on $HOSTS."
echo "   with the machinefile $host_file"
#====== Do not touch up to here =================
#================================================

#mpirun -machinefile ${host_file} ${EXEC} clump_NHI18_fcov1.in
mpirun -machinefile ${host_file} ${EXEC} clump_NHI18_fcov2.in
#mpirun -machinefile ${host_file} ${EXEC} clump_NHI18_fcov5.in
#mpirun -machinefile ${host_file} ${EXEC} clump_NHI18_fcov20.in
mpirun -machinefile ${host_file} ${EXEC} clump_NHI18_fcov50.in
mpirun -machinefile ${host_file} ${EXEC} clump_NHI18_fcov100.in

#mpirun -machinefile ${host_file} ${EXEC} clump_NHI18_fcov1_gauss1.in
#mpirun -machinefile ${host_file} ${EXEC} clump_NHI18_fcov2_gauss1.in
#mpirun -machinefile ${host_file} ${EXEC} clump_NHI18_fcov5_gauss1.in
#mpirun -machinefile ${host_file} ${EXEC} clump_NHI18_fcov20_gauss1.in
#mpirun -machinefile ${host_file} ${EXEC} clump_NHI18_fcov50_gauss1.in
#mpirun -machinefile ${host_file} ${EXEC} clump_NHI18_fcov100_gauss1.in

mpirun -machinefile ${host_file} ${EXEC} clump_NHI20_fcov1.in
mpirun -machinefile ${host_file} ${EXEC} clump_NHI20_fcov2.in
mpirun -machinefile ${host_file} ${EXEC} clump_NHI20_fcov5.in
mpirun -machinefile ${host_file} ${EXEC} clump_NHI20_fcov20.in
mpirun -machinefile ${host_file} ${EXEC} clump_NHI20_fcov50.in
mpirun -machinefile ${host_file} ${EXEC} clump_NHI20_fcov100.in
#mpirun -machinefile ${host_file} ${EXEC} clump_NHI20_fcov1_gauss1.in
#mpirun -machinefile ${host_file} ${EXEC} clump_NHI20_fcov2_gauss1.in
#mpirun -machinefile ${host_file} ${EXEC} clump_NHI20_fcov5_gauss1.in
#mpirun -machinefile ${host_file} ${EXEC} clump_NHI20_fcov20_gauss1.in
#mpirun -machinefile ${host_file} ${EXEC} clump_NHI20_fcov50_gauss1.in
#mpirun -machinefile ${host_file} ${EXEC} clump_NHI20_fcov100_gauss1.in
