#!/bin/bash -x
#PBS -N jobname
#PBS -l nodes=1:ppn=28
#PBS -j oe
#PBS -q normal

# Use the Phonopy
# Reference: “First-principles Phonon Calculations with Phonopy and Phono3py”, Atsushi Togo, J. Phys. Soc. Jpn., 92, 012001-1-21 (2023)
# 			 “Implementation strategies in phonopy and phono3py”, Atsushi Togo, Laurent Chaput, Terumasa Tadano, and Isao Tanaka, J. Phys. Condens. Matter 35, 353001-1-22 (2023)

echo $PBS_O_WORKDIR
NCPUS=`wc -l $PBS_NODEFILE | awk '{print $1}'`
echo $NCPUS

cd $PBS_O_WORKDIR
module load vasp/5.4.4-intel2019

MPI=mpirun
VASP=vasp_std_2D

jobconv='reached'
jobfinish='Voluntary'


spdim="3 3 1"
INCARpath="$HOME/Input/INCAR-phono"
failedpath="$HOME/log/failed.txt"
bandconfpath="$HOME/Input/band.conf"
echo "#################" >> $failedpath
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

    num_poscars=$(ls -d POSCAR-0* 2>/dev/null | wc -l)
    count=1

    cp POSCAR-001 POSCAR
    echo -e "102\n2\n0.04" | vaspkit >/dev/null 2>&1

    for poscar_file in POSCAR-0*; do
        cp "$poscar_file" POSCAR
        $MPI $VASP > runlog

        if [ -f runlog ] && grep -q $jobconv OUTCAR && grep -q $jobfinish OUTCAR; then
            mv vasprun.xml vasprun.xml-$count
        fi
        count=$((count+1))
    done

    file_list=$(eval echo vasprun.xml-{1..$num_poscars})
    missing_files=$(ls $file_list 2>/dev/null | wc -l)

    if [ "$missing_files" -eq "$num_poscars" ]; then
        phonopy -f vasprun.xml-{1..$num_poscars}
        phonopy -c POSCAR-unitcell $bandconfpath -p -s
    else
        echo "Some vasprun.xml files are missing."
        echo "$file" >> $failedpath
    fi

    cd $start_dir
done

