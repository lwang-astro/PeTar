#!/bin/bash

#Initial data file generator, read data file (m,r[3],v[3])

read -p "Data filename: " fname
[ -z $fname ] && echo 'No file name found!' && exit
read -p "output filename (fname.input): " fout
[ -z $fout ] && fout=$fname.input
read -p "Ignore line number (0): " igline
[ -z $igline ] && igline=0
read -p "With stellar evolution y/n (y): " seflag
[ -z $seflag ] && seflag=y

n=`wc -l $fname|cut -d' ' -f1`
n=`expr $n - $igline`
if [[ $seflag -eq 'y' ]]; then
    read -p "Stellar radius (0)" radius
    [ -z $radius ] && radius=0
    awk -v n=$n -v ig=$igline 'BEGIN{print 0,n,0} {if(NR>ig) print $LINE,'$radius',0,0,0,0,NR-ig,0,0,0,0,0,0,0,0,0}' $fname >$fout
else
    awk -v n=$n -v ig=$igline 'BEGIN{print 0,n,0} {if(NR>ig) print $LINE,0,NR-ig,0,0,0,0,0,0,0,0,0}' $fname >$fout
fi
