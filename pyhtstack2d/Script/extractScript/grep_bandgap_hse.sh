#!/bin/bash

# Batch data extraction via shell script
# Band structure calculation for TMD bilayers

# Remove the old output file if it exists, and add the header
rm -f hse_gap_inf.txt
echo "material cella gaptype gapvalue E VBMLx VBMLy cordprefix element elementnum hm1 hm2 h1 h2 h3 h4 Ef Evacuum" >> hse_gap_inf.txt

# Loop through each row in the select_list_gt3.txt file
while read -r rows; do
    cd "$rows" || { echo "Failed to change directory to $rows"; continue; }

    # Extract values from the path and CONTCAR file
    material=$(echo "$rows" | awk -F '/' '{print $1}')
    cordprefix=$(echo "$rows" | awk -F '/' '{print $2}')
    cordf=$(echo "$rows" | awk -F '/' '{print $3}')
    cella=$(awk 'NR==3 {print $2}' CONTCAR)
    element=$(sed -n '6p' CONTCAR | sed 's/[ ][ ]*/-/g')
    elementnum=$(sed -n '7p' CONTCAR | sed 's/[ ][ ]*/-/g')
    hm1=$(awk 'NR==9 {print $4}' CONTCAR)
    hm2=$(awk 'NR==10 {print $4}' CONTCAR)
    h1=$(awk 'NR==11 {print $4}' CONTCAR)
    h2=$(awk 'NR==12 {print $4}' CONTCAR)
    h3=$(awk 'NR==13 {print $4}' CONTCAR)
    h4=$(awk 'NR==14 {print $4}' CONTCAR)

    # Move to the 'hse' directory and check for BAND_GAP file
    cd ./hse || { echo "Failed to change directory to ./hse in $rows"; continue; }
    if [ ! -f BAND_GAP ]; then
        cd ../hse2 || { echo "Failed to change directory to ../hse2 in $rows"; continue; }
    fi

    # Run VASPkit tasks
    vaspkit -task 252 >/dev/null 2>&1
    vaspkit -task 927 > Vacuum.txt

    # Extract data from BAND_GAP and OUTCAR
    gaptype=$(awk '/Character/ {print $3}' BAND_GAP)
    gapvalue=$(awk '/Gap/ {print $4}' BAND_GAP)
    VBMLx=$(awk '/Location of VBM/ {print $4}' BAND_GAP)
    VBMLy=$(awk '/Location of VBM/ {print $5}' BAND_GAP)
    Ef=$(awk '/E-fermi/ {print $3}' OUTCAR)
    Evacuum=$(awk '/Vacuum Level/ {print $4}' Vacuum.txt)
    E=$(awk 'END {print $5}' OSZICAR)

    # Move back to the top directory and append the results
    cd ../../../.. || { echo "Failed to change directory back to root"; continue; }
    echo "$material $cella $gaptype $gapvalue $E $VBMLx $VBMLy $cordprefix $element $elementnum $hm1 $hm2 $h1 $h2 $h3 $h4 $Ef $Evacuum" >> hse_gap_inf.txt

done < material.txt


