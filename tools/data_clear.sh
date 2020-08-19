#!/bin/bash

suffixes=(esc group sse bse status prof.rank)
tindices=(1 3 21 13 1 2)
# sn single 12 binary 14
unset tcrit

until [[ `echo x$1` == 'x' ]]
do
    case $1 in
	-h) shift;
	    echo 'Remove data after a given time criterion if users want to restart simulations at middle and overwrite the output files (file suffixes: '$suffixes')';
	    echo 'Usage: petar.gether [options] [data filename prefix]';
	    echo '       data filename prefix is defined by "petar -f", defaulted case is "data".'
	    echo 'Options:';
	    echo '  -t: time criterion (default: none)'
	    echo '  -n: MPI processes number (default: auto detect)';
	    echo '  -i: replace backup files (default: false)';
	    exit;;
	-n) shift; nmpi=$1; shift;;
	-t) shift; tcrit=$1; shift;;
	-i) rmi=1; shift;;
	*) fname=$1;shift;;
    esac
done

if [ ! -e $fname ] | [ -z $fname ] ; then
    echo 'Error, file name not provided' 
    exit
fi

if [ -z $tcrit ]; then
    echo 'Time criterion not provided'
    exit
fi

echo 'data filename prefix: '$fname
echo 'time criterion: '$tcrit

# esc group
for ((i=0;i<6;i=i+1))
do
    s=${suffixes[$i]}
    tindex=${tindices[$i]}
    file=$fname.$s
    echo $file'; time column: '$tindex
    if [ ! -z $nmpi ]; then
	nend=`expr $nmpi - 1`
	lst=`seq 0 $nend|awk -v f=$file '{printf("$s.$d\n", f,$1)}`
    else
	lst=`ls |egrep $file'.[0-9]+$'`
    fi
    for f in $lst
    do
	echo 'process '$f
	if [ -e $f.bk ]; then
	    if [ ! -z $rmi ]; then
		rm -i $f.bk
		mv $f $f.bk
	    fi
	else
	    mv $f $f.bk
	fi
	if [ $s == 'sse' ]; then
	    awk -v t=$tcrit -v ti=$tindex '{if (($1!="SN_kick" && $ti<=t) || ($1=="SN_kick" && $12<=t)) print $LINE}' $f.bk >$f
	elif [ $s == 'bse' ]; then
	    awk -v t=$tcrit -v ti=$tindex '{if ($1=="SN_kick") {if ($14<=t) print $LINE} else if ($1=="Dynamic_merge:") {if($35<=t) print $LINE} else if ($ti<=t) print $LINE;}' $f.bk >$f
	else
	    awk -v t=$tcrit -v ti=$tindex '{if ($ti<=t) print $LINE}' $f.bk >$f
	fi
    done
done

#



