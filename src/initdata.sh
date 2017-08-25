#!/bin/bash

#Initial data file generator, read data file (m,r[3],v[3])

read -p "Data filename: " fname
[ -z $fname ] && echo 'No file name found!' && exit
read -p "Ignore line number (0): " igline
[ -z $igline ] && igline=0

n=`wc -l $fname|cut -d' ' -f1`
n=`expr $n - $igline`
awk -v n=$n -v ig=$igline 'BEGIN{print 0,n,n,0,0,0} {if(NR>ig) print $LINE,0,0,NR-ig,0,0,0,0,0,0}' $fname >$fname.input
