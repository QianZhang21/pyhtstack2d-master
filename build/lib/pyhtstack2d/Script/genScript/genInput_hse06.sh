#!/bin/bash

for rows in `cat material.txt`
do
cd $rows
cd ./scf
vaspkit -task 302 >/dev/null 2>&1
echo -e "251\n2\n0.03\n0.03" | vaspkit >/dev/null 2>&1
cp ~/work/material_design/high_throughput/INCAR_file/PBE_INCAR INCAR
cd ..
if [ ! -d hse ];then
mkdir hse
fi
cd ./hse
cp ~/work/material_design/high_throughput/INCAR_file/HSE_INCAR INCAR
cp ../scf/KPOINTS ./
cp ../scf/PO* ./
cp ../scf/KPATH.in ./
cd ..
cd ../../..
done











