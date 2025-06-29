#!/bin/sh

rm mid.txt
for file in `ls $1`
do
    if [ ! -d $file ];then
        continue
    fi
    if [ `grep -c "reached required accuracy" $file"/AA/cord1/runlog"` -ne 0 ]; then
        mag=`tail -1 $file"/AA/cord1/OSZICAR" | awk '{print $10}'` 
        E0=`tail -1 $file"/AA/cord1/OSZICAR" | awk '{print $5}'`
        echo $file" "$mag" "$E0 >> mid.txt
    else
        echo $file
    fi
done
