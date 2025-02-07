#!/bin/bash -x
#SBATCH -J HT2
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

echo "#################" >> failed.txt

jobconv='reached'
jobfinish='Voluntary'


for file in `cat maglist.txt`
do
    cd $file
    
    if [ ! -d scf ]; then
        mkdir scf
    fi
    cd scf
    ln -nfs ../POSCAR POSCAR
    ln -nfs ../POTCAR POTCAR

    if [ -f ../POSCAR221 ]; then
        ln -nfs ../AFM1/KPOINTS 
    else
        ln -nfs ../../../../KPOINTS-scf KPOINTS
    fi

    if [ -f ../MAGMOM ]; then
        cp ../../../../INCAR-scf INCAR
        cat ../MAGMOM >> INCAR
    else
        ln -nfs ../../../../INCAR-scf INCAR
    fi

    $MPI $VASP > runlog
    if [ ! -f OUTCAR ] || [ `grep -c $jobconv OUTCAR` -eq 0 ] || [ `grep -c $jobfinish OUTCAR` -eq 0 ]; then
        cd ../../../../
	echo $file","scf-failed >> failed.txt
	continue
    fi
    cd ..

    if [ ! -d band ]; then
        mkdir band
    fi
    cd band
    ln -nfs ../POSCAR POSCAR
    ln -nfs ../POTCAR POTCAR

    if [ -f ../MAGMOM ]; then
        cp ../../../../INCAR-band INCAR
        cat ../MAGMOM >> INCAR
    else
        ln -nfs ../../../../INCAR-band INCAR
    fi

    ln -nfs ../../../../KPOINTS-band KPOINTS
    ln -nfs ../scf/CHG CHG
    ln -nfs ../scf/CHGCAR CHGCAR
    $MPI $VASP > runlog
    if [ ! -f OUTCAR ] || [ `grep -c $jobconv OUTCAR` -eq 0 ] || [ `grep -c $jobfinish OUTCAR` -eq 0 ]; then
        echo $file","band-failed >> ../../../../failed.txt
    fi
    cd ../../../..
done
