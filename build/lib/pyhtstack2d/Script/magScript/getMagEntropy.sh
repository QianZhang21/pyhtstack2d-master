#!/bin/bash -x

rm mag.txt entropy.txt

for file in `ls $1`
do
    if [ ! -d $file ];then
        continue
    fi
    cd $file
    for file_0 in `ls $1`
    do
        if [ ! -d $file_0 ]; then
            continue
        fi
        cd $file_0
        for file_1 in `ls $1`
        do
            if [ ! -d $file_1 ]; then
                continue
            fi
            cd $file_1
            cd opt
            if [ ! -f runlog ] || [ `grep -c "reached required accuracy" runlog` -eq 0 ];then
                echo $file/$file_0/$file_1
            fi
            mag=`tail -1 OSZICAR | awk '{print $10}'`
            echo $file/$file_0/$file_1","$mag >> ../../../../mag.txt
            entropy=`grep "entropy T" OUTCAR | tail -1 | awk '{print $5}'`
            echo $file/$file_0/$file_1","$entropy  >> ../../../../entropy.txt
            cd ../
            cd ..
        done
        cd ..
    done
    cd ..
done
