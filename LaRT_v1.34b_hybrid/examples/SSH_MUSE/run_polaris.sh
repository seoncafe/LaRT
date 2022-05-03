#!/bin/bash

#SBATCH --job-name="MUSE"          # Job name
#SBATCH --output=out_polaris.txt   # Name of stdout output file (%j expands to jobId)
#SBATCH --open-mode=append         # open the output and error files using append output mode (default is truncate mode)

##SBATCH --time=4-00:00:00     # time limit, 1 hour, for this job.
##SBATCH --partition=sm        # select sm partition (or queue), which is composed of supermicro servers
##SBATCH --nodes=7             # requested number of nodes
##SBATCH --ntasks-per-node=40  # requested number of mpi tasks per node

#SBATCH --time=7-00:00:00     # time limit, 1 hour, for this job.
#SBATCH --partition=ibm       # select sm partition (or queue), which is composed of supermicro servers
#SBATCH --nodes=20            # requested number of nodes
#SBATCH --ntasks-per-node=16  # requested number of mpi tasks per node

module purge                  # remove all the modules
#module load intel/2018.5.274   # load intel compiler

#export I_MPI_FABRICS=shm:dapl # use shared memory and infiniband as interconnects of the MPI job
export I_MPI_FABRICS=shm:ofi # use shared memory and infiniband as interconnects of the MPI job
unset I_MPI_SHM_LMT

cd $SLURM_SUBMIT_DIR          # change your working directory where you launch your job.

EXEC=../LaRT_peel_FS.x

mpirun $EXEC MUSE1185.in
mpirun $EXEC MUSE1185.in

mpirun $EXEC MUSE0082.in
mpirun $EXEC MUSE0082.in

mpirun $EXEC MUSE6905.in
mpirun $EXEC MUSE6905.in

#mpirun $EXEC MUSE1343.in
#mpirun $EXEC MUSE0053.in
#mpirun $EXEC MUSE0171.in
#mpirun $EXEC MUSE0547.in
#mpirun $EXEC MUSE0364.in
