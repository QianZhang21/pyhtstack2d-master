#!/bin/bash -x
#SBATCH -J HT
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



for file in `cat mass_mid.txt`
do
    if [ ! -d $file ];then
        continue
    fi
    cd $file
	  mkdir mass
    cd mass
    ln -nfs ../band/POTCAR
    ln -nfs ../band/POSCAR
    ln -nfs ../../../../INCAR-mass INCAR
    mv ../VPKIT.in ./
	  echo -e "913\n2\n0.02" | vaspkit >/dev/null 2>&1
    $MPI $VASP > runlog
    sed -i '1s/.*/2/' VPKIT.in
    echo -e "913" | vaspkit >> mass.txt
    cd ..
    if [ -f cbm_VPKIT.in ]; then
        mkdir mass_cbm
        cd mass_cbm
        ln -nfs ../band/POTCAR
        ln -nfs ../band/POSCAR
        ln -nfs ../mass/INCAR
        mv ../cbm_VPKIT.in VPKIT.in
        echo -e "913\n2\n0.02" | vaspkit >/dev/null 2>&1
        $MPI $VASP > runlog
        sed -i '1s/.*/2/' VPKIT.in
        echo -e "913" | vaspkit >> mass.txt
        cd ..
    fi
    cd ../../..
done


