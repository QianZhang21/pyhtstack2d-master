#!/bin/bash -x
#SBATCH -J HT
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -p short
#SBATCH -o %j.log
#SBATCH -e %j.err

echo Time is `date`
echo Directory is $PWD
echo This job runs on the following nodes:
echo $SLURM_JOB_NODELIST
echo This job has allocated $SLURM_JOB_CPUS_PER_NODE cpu cores.

cd $SLURM_SUBMIT_DIR
module load vasp/5.4.4-intel2019

MPI=mpirun
VASP=vasp_std_2D

for file in `cat entropylist.txt`
do
    cd $file
    $MPI $VASP > runlog
    cd ../../..
done
