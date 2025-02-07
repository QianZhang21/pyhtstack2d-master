#!/bin/bash -x
#SBATCH -J band
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -p long
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

jobconv='reached'
jobfinish='Voluntary'

for file in `cat nomaglist.txt`
do
    cd $file
    if [ ! -d scf ]; then
        mkdir scf
    fi
    cd scf
    ln -nfs ../CONTCAR POSCAR
    ln -nfs ../POTCAR POTCAR
    ln -nfs ../../../../INCAR-scf INCAR
    ln -nfs ../../../../KPOINTS-scf KPOINTS
    $MPI $VASP > runlog
    if [ ! -f OUTCAR ] || [ `grep -c $jobconv OUTCAR` -eq 0 ] || [ `grep -c $jobfinish OUTCAR` -eq 0 ]; then
        echo $file",scf-failed" >> ../../../../failed.txt
        cd ../../../..
        continue
    fi
    cd ..

    if [ ! -d band ]; then
        mkdir band
    fi
    cd band
    ln -nfs ../CONTCAR POSCAR
    ln -nfs ../POTCAR POTCAR
    ln -nfs ../../../../INCAR-band INCAR
    ln -nfs ../../../../KPOINTS-band KPOINTS
    ln -nfs ../scf/CHG CHG
    ln -nfs ../scf/CHGCAR CHGCAR
    if [ ! -f OUTCAR ] || [ `grep -c $jobconv OUTCAR` -eq 0 ] || [ `grep -c $jobfinish OUTCAR` -eq 0 ]; then
        echo $file",band-failed" >> ../../../../failed.txt
    fi
    cd ../../../..
done
