#!/bin/bash
source ~/.uftools.sh

lst=`powseq -n 11 -s 5 0.5`
npath=`pwd`


T=0.3
r=0.01

#n=5000
#rdir=n$n.r$r.T$T.grad

input=$npath'/input/Plum.N5k.input'
rdir=n5000.r$r.T$T.imf.grad

[ -d $rdir ] || mkdir $rdir
cd $rdir

for kd in kdk kdkdk2
do
    rm -f $kd.err
    for dt in $lst
    do
	echo $kd.$dt
	fac=`echo $dt'/0.00048828125*8'|bc -l`
	#$npath/nbody.out.$kd -n $n -r $r -t 0.125 -o $dt -s $dt --dt-max-factor $fac -T $T &>$kd.$dt.log
	$npath/nbody.out.$kd -r $r -t 0.125 -o $dt -s $dt --dt-max-factor $fac -T $T --r-bin 1e-10 $input &>$kd.$dt.log
	egrep Enow-Einit $kd.$dt.log |awk -v dt=$dt 'BEGIN{t=0.0;h=0.0} {terr=$8>0?$8:-$8; herr=$11>0?$11:-$11; if (t<terr) t=terr; if (h<herr) h=herr;} END{print dt,t,h}' >>$kd.err
    done
done

cd $npath
