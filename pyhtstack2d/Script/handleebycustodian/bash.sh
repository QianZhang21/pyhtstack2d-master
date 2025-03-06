#!/bin/bash -x
#SBATCH -J HT3
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
module load vasp/5.4.4-intel2019 # Load the vasp module
runvasp="path/script/custodian_vasp.py" # Specify the path of custodian_vasp.py

# ====================================================================================================
# For 1-level directory
for file in `ls $1`
do
    if [ ! -d $file ];then
        continue
    fi
    cd $file
    ln -nfs $runvasp custodian_vasp.py
    python custodian_vasp.py
    cd ..
done
# ====================================================================================================


## ====================================================================================================
## For 2-level directory
#for file in `ls $1`
#do
#    if [ ! -d $file ];then
#        continue
#    fi
#    cd $file
#    for file_0 in `ls $1`
#    do
#        if [ ! -d $file_0 ]; then
#            continue
#        fi
#        cd $file_0
#        ln -nfs $runvasp custodian_vasp.py
#        python custodian_vasp.py
#        cd ..
#    done
#    cd ..
#done
## ====================================================================================================


# ====================================================================================================
## For 3-level directory
#for file in `ls $1`
#do
#    if [ ! -d $file ];then
#        continue
#    fi
#    cd $file
#    for file_0 in `ls $1`
#    do
#        if [ ! -d $file_0 ]; then
#            continue
#        fi
#        cd $file_0
#        for file_1 in `ls $1`
#        do
#            if [ ! -d $file_1 ]; then
#                continue
#            fi
#            cd $file_1
#            ln -nfs $runvasp custodian_vasp.py
#            python custodian_vasp.py
#            cd ..
#        done
#        cd ..
#    done
#    cd ..
#done
## ====================================================================================================