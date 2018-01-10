#!/bin/bash
k=0.0625
sub='tt'
for i in `seq 7`
do
    echo $i' s='$k
    rm -f data.* global.dat
    OMP_STACKSIZE=100000 ./nbody.out -i 2 -t 100.0 -R 0.1 -r 0.3 -b 0.3 -o 0.0625 -D $k -e 0.0 input/p3t.input 1>result/nbody.log.p3t.$sub$i 2>&1 
    #echo 'input/p3t.input' >datalist
    ls data.*|sort -n -k1.6 > datalist
    #./getdata -i 1 input/p3t.input 1>getdata.log 2>&1
    ./getdata -i 3 -f datalist 1>getdata.log 2>&1
    mv global.dat result/nbody.dat.p3t.$sub$i
    k=`echo 'scale=10; '$k'/2'|bc`
done
