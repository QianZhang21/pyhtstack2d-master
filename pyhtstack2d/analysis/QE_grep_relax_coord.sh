#!/bin/bash


# QE relax Resulting Atomic Position Extraction Script

rm -f relax_coord.txt

nat=`grep 'nat' C6H6.0_relax.in | awk -F "=" '{print $2}' | sed 's/[^0-9]//g'`
finalstr=`cat C6H6.0_relax.out | grep -n 'Begin final coordinates' | awk -F ":" '{print $1}'`
addn=3
let nbegin=$finalstr+$addn
let nfinal=$nbegin+$nat

for ((i=$nbegin; i<$nfinal; i++))
do
element=`awk 'NR=='$i  C6H6.0_relax.out | awk '{print $1}'`
coord1=`awk 'NR=='$i  C6H6.0_relax.out | awk '{print $2}'`
coord2=`awk 'NR=='$i  C6H6.0_relax.out | awk '{print $3}'`
coord3=`awk 'NR=='$i  C6H6.0_relax.out | awk '{print $4}'`
echo $element'        '$coord1'   '$coord2'   '$coord3 >> relax_coord.txt
done


