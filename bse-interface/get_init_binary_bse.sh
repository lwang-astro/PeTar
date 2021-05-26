#!/bin/bash
unset fout
unset nbin
unset fname

until [[ `echo x$1` == 'x' ]]
do
    case $1 in
	-h) shift;
	    echo 'Read PeTar initial data and generate the input file for binary stellar evolution (petar.[mo]bse)';
	    echo 'The keplerorbit tool from SDAR/sample/Kepler needs to be installed first'
	    echo 'Usage: petar.get.init.binary.bse [options] [input data filename]';
	    echo 'Input data file should have the astronomical unit: Msun, pc, pc/Myr';
	    echo 'Options:';
	    echo '  -f: output file (petar input data) name (default: intput file name + ".bin0")';
	    echo '      Each line: m1[Msun] m2[Msun] 0 0 period[Myr] ecc 0'
	    echo '  -h: help';
	    echo '  -b: number of binaries (1)';
	    exit;;
	-f) shift; fout=$1; shift;;
	-b) shift; nbin=$1; shift;;
	*) fname=$1;shift;;
    esac
done

if [ ! -e $fname ] | [ -z $fname ] ; then
    echo 'Error, file name not provided' 
    exit
fi

[ -z $fout ] && fout=$fname.bin0
[ -z $nbin ] && nbin=0

echo 'generate binary file: '$fout'; number of binaries: '$nbin


awk '{if(NR>1) print $1,$2,$3,$4,$5,$6,$7}' $fname >$fname.tmp_

echo $nbin >$fout
keplerorbit -i -n $nbin -u 4 $fname.tmp_ |awk '{print $16,$17,0,0,$14,$9,0}' >>$fout
rm -f $fname.tmp_
