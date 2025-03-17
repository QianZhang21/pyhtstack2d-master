#!/bin/bash -x
#SBATCH -J HT
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p short
#SBATCH -o %j.log
#SBATCH -e %j.err

echo Time is `date`
echo Directory is $PWD
echo This job runs on the following nodes:
echo $SLURM_JOB_NODELIST
echo This job has allocated $SLURM_JOB_CPUS_PER_NODE cpu cores.

cd $SLURM_SUBMIT_DIR

for file in `cat mass_mid.txt`
do
    cd $file
    cd band
    sumo-bandstats
    if [ `grep -c "ERROR" sumo-bandstats.log` -eq 0 ];then
        echo "=============" $file >> ../../../../../mass_sumo.txt
        nstar=$(grep -n -m 1 "Hole effective masses" sumo-bandstats.log | awk -F: '{print $1}')
        sed -n "${nstar},\$p" sumo-bandstats.log >> ../../../../../mass_sumo.txt
    fi
    cd ..
    cd ../../..
done
