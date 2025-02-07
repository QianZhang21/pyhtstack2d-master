#!/bin/bash -x

# Create the maglist.txt file if all the materials in the work folder are magnetic systems.

rm -f maglist.txt
path=$(pwd -P) 
for file in `ls $path`; do  
if [ -d "$file" ]; then
materialpath="$file"
echo $materialpath >> maglist.txt
fi
done




### For nested directories, V1
#multilevel=2
#path=$(pwd -P)
#find "$path" -mindepth $multilevel -type d | while read -r dir; do
#    echo "$dir" >> maglist.txt
#done




### For nested directories, V2
#path=$(pwd -P)
#for file in `ls $path`; do
#if [ -d "$file" ]; then
#for file2 in `ls $file`; do
#materialpath="$file/$file2"
#if [ -d $materialpath ]; then
#echo $materialpath >> maglist.txt
#fi
#done
#fi
#done

