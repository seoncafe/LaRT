#!/bin/bash
# Run the clumpy-sphere examples for LaRT_v2.00 (use_clump_medium = .true.).
# Adjust NP (MPI rank count) and LART path as needed.

#LART=../../LaRT.x
#NP=4
#
#echo "=== Clumpy sphere: tau0=1e3, N_cov=5, static ==="
#mpirun -np ${NP} ${LART} clump_tau3_fcov5.in
#
#echo "=== Clumpy sphere: tau0=1e4, N_cov=5, static ==="
#mpirun -np ${NP} ${LART} clump_tau4_fcov5.in
#
#echo "=== Clumpy sphere: tau0=1e3, N_cov=20, static (dense covering) ==="
#mpirun -np ${NP} ${LART} clump_tau3_fcov20.in
#
#echo "=== Clumpy sphere: tau0=1e3, N_cov=5, Hubble expansion Vexp=200 km/s ==="
#mpirun -np ${NP} ${LART} clump_tau3_fcov5_hubble.in
#
#echo "=== Clumpy sphere: tau0=1e3, N_cov=5, clump sigma_v=50 km/s ==="
#mpirun -np ${NP} ${LART} clump_tau3_fcov5_sigv.in
#
#echo "=== Clumpy sphere: tau0=1e3, N_cov=5, peel-off observer ==="
#mpirun -np ${NP} ${LART} clump_tau3_fcov5_peel.in
#
#================================================
#================================================
exec < /dev/null 2>&1
trap "" HUP

EXEC=../../LaRT.x

#HOSTS=mocafe
#HOSTS=mocafe,lart2
#HOSTS=mocafe,lart1,lart2,lart3,lart4
#HOSTS=lart4,lart3,lart2,lart1
#HOSTS=lart4
HOSTS=lart4,lart3,lart2

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

#echo "=== Clumpy sphere: tau0=1e3, N_cov=5, static ==="
#mpirun -machinefile ${host_file} ${EXEC} clump_tau3_fcov5.in
#
#echo "=== Clumpy sphere: tau0=1e4, N_cov=5, static ==="
#mpirun -machinefile ${host_file} ${EXEC} clump_tau4_fcov5.in
#
#echo "=== Clumpy sphere: tau0=1e3, N_cov=20, static (dense covering) ==="
#mpirun -machinefile ${host_file} ${EXEC} clump_tau3_fcov20.in
#
#echo "=== Clumpy sphere: tau0=1e3, N_cov=5, Hubble expansion Vexp=200 km/s ==="
#mpirun -machinefile ${host_file} ${EXEC} clump_tau3_fcov5_hubble.in
#
#echo "=== Clumpy sphere: tau0=1e3, N_cov=5, clump sigma_v=50 km/s ==="
#mpirun -machinefile ${host_file} ${EXEC} clump_tau3_fcov5_sigv.in

#echo "=== Clumpy sphere: tau0=1e3, N_cov=5, peel-off observer ==="
#mpirun -machinefile ${host_file} ${EXEC} clump_tau3_fcov5_peel.in

mpirun -machinefile ${host_file} ${EXEC} clump_NHI18_fcov1.in
mpirun -machinefile ${host_file} ${EXEC} clump_NHI18_fcov5.in
mpirun -machinefile ${host_file} ${EXEC} clump_NHI18_fcov20.in
