#!/bin/bash
file=cadata.txt
count=`wc -l $file | cut -f1 -d" "`;
for i in `seq 1 $count`
do
    echo $i
    num=`sed -n "${i}p" $file`;
    ./radialdistribution.py $num >move$num.txt
done
