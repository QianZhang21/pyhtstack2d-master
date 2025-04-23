#!/bin/sh

for file in `ls $1`
do
    if [ -d $file/AA_7 ];then
        echo $file
        mkdir $file"_1"
        mv $file/AA_{7..12} $file"_1"
    fi
done
