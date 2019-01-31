#!/bin/bash
source ~/.uftools.sh

lst=`powseq -n 13 -s 11 0.5`

r=0.01

for kd in kdk kdkdk2 kdkdk4
do
    #rm -f $kd.$r.err
    for dt in $lst
    do
	echo $kd.$dt.$r
	fac=`echo $dt'/0.00048828125*4'|bc -l`
	./nbody.out.$kd -n 1000 -r $r -t 0.25 -o $dt -s $dt --dt-max-factor $fac -T 0.0 &>$kd.$dt.$r.log
	egrep Enow-Einit $kd.$dt.$r.log |awk -v dt=$dt 'BEGIN{t=0.0;h=0.0} {terr=$8>0?$8:-$8; herr=$11>0?$11:-$11; if (t<terr) t=terr; if (h<herr) h=herr;} END{print dt,t,h}' >>$kd.$r.err
    done
done
