#!/bin/bash -x


for file in `cat mass_mid.txt`
do
    if [ ! -d $file ];then
        continue
    fi
    cd $file
    cd mass
    echo "=============" $file >> ../../../../../mass.txt
    #echo "**************VBM***********" >> ../../../../../mass.txt
    nstar=$(grep -n -m 1 "Band Index" mass.txt | awk -F: '{print $1}')
    nfinal=$((nstar + 4))
    sed -n "${nstar},${nfinal}p" mass.txt >> ../../../../../mass.txt
    cd ..
    if [ -d mass_cbm ]; then
        cd mass_cbm
        #echo "**************CBM***********" >> ../../../../../mass.txt
        nstar=$(grep -n -m 1 "Band Index" mass.txt | awk -F: '{print $1}')
        nfinal=$((nstar + 4))
        sed -n "${nstar},${nfinal}p" mass.txt >> ../../../../../mass.txt
        cd ..
    fi
    cd ../../..
done

