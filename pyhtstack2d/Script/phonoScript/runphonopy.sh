#!/bin/bash -x
#PBS -N jobname
#PBS -l nodes=1:ppn=28
#PBS -j oe
#PBS -q normal

# Use the Phonopy-DFPT
# Reference: “First-principles Phonon Calculations with Phonopy and Phono3py”, Atsushi Togo, J. Phys. Soc. Jpn., 92, 012001-1-21 (2023)
# 			 “Implementation strategies in phonopy and phono3py”, Atsushi Togo, Laurent Chaput, Terumasa Tadano, and Isao Tanaka, J. Phys. Condens. Matter 35, 353001-1-22 (2023)

echo "Working directory: $PBS_O_WORKDIR"
NCPUS=`wc -l $PBS_NODEFILE | awk '{print $1}'`
echo "Number of CPUs: $NCPUS"

cd $PBS_O_WORKDIR
module load vasp/5.4.4-intel2019

MPI=mpirun
VASP=vasp_std_2D

echo "#################" >> $HOME/log/failed.txt

jobconv='reached'
jobfinish='Voluntary'

spdim="3 3 1"
INCARpath="$HOME/Input/INCAR-phono"
failedpath="$HOME/log/failed.txt"
bandconfpath="$HOME/Input/band.conf"

start_dir=$(pwd)

for file in `cat phono.txt`
do
    if [ ! -d "$file" ]; then
        continue
    fi
    cd "$file"

    mkdir -p phono
    cd phono
    cp ../POSCAR POSCAR
    ln -nfs ../POTCAR POTCAR
    cp $INCARpath INCAR
    phonopy -d --dim=$spdim
    mv POSCAR POSCAR-unitcell
    mv SPOSCAR POSCAR
	  echo -e "102\n2\n0.04" | vaspkit >/dev/null 2>&1
    $MPI $VASP > runlog
    if [ -f OUTCAR ] && grep -q "$jobconv" OUTCAR && grep -q "$jobfinish" OUTCAR; then
        phonopy --fc vasprun.xml > phono.log
        if [ -f FORCE_CONSTANTS ]; then
            phonopy -c POSCAR-unitcell $bandconfpath -p -s
        else
            echo "FORCE_CONSTANTS not generated for $file" >> $failedpath
        fi
    else
        echo "VASP job failed for $file" >> $failedpath
    fi
    cd $start_dir
done
