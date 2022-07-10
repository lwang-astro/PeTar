#!/bin/bash

# sn single 12 binary 14
tindex_sse_type_change=21
tindex_sse_sn_kick=12
tindex_bse_type_change=23
tindex_bse_dyn_merge=38
tindex_bse_sn_kick=14

ncol_sse_type_change=22
ncol_sse_sn_kick=13
ncol_bse_type_change=46
ncol_bse_dyn_merge=49
ncol_bse_sn_kick=15

unset tcrit

until [[ `echo x$1` == 'x' ]]
do
    case $1 in
	-h) shift;
	    echo 'Remove data after a given time criterion if users want to restart simulations at middle and overwrite the output files (file suffixes: '${suffixes[@]}')';
	    echo 'Usage: petar.data.clear [options] [data filename prefix]';
	    echo '       data filename prefix is defined by "petar -f", defaulted case is "data".'
	    echo 'Options:';
	    echo '  -t: time criterion (default: none)'
	    echo '  -n: MPI processes number (default: auto detect)';
	    echo '  -b: use previous backup files instead of replacing (default: replacing)';
	    echo '  --less-dyn-merge: less output mode for [prefix].*bse*.dynmical_merge (three columns less, for the PeTar version before Sep 10, 2020)'
	    echo '  --less-type-change: less output mode for [prefix].*bse*.type_change (20 columns less, for the PeTar version before Jun 2, 2022)'
	    exit;;
	-n) shift; nmpi=$1; shift;;
	-t) shift; tcrit=$1; shift;;
	-b) rmi=1; shift;;
	--less-dyn-merge) tindex_bse_dyn_merge=35; ncol_bse_dyn_merge=46;shift;; 
	--less-type-change) tindex_bse_type_change=13;ncol_bse_type_change=26;shift;;
	*) fname=$1;shift;;
    esac
done

suffixes=(esc group sse mosse sseEmp bse mobse bseEmp status prof.rank)
tindices=(1 3 0 0 0 0 0 0 1 2)
nsuffixes=${#suffixes[@]}

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

# check consistence for the number of columns
ncol=`egrep -m 1 'SN_kick' $fname.*sse*.[0-9]|wc -w`
if [[ $ncol != $ncol_sse_sn_kick ]] && [[ $ncol != 0 ]]; then
    echo 'Error! column number not matches for SSE SN kick, should be '$ncol_sse_sn_kick', the file has '$ncol'.'
    exit
fi
ncol=`egrep -v -m 1 'SN_kick' $fname.*sse*.[0-9]|wc -w`
if [[ $ncol -ne $ncol_sse_type_change ]] && [[ $ncol -ne 0 ]]; then
    echo 'Error! column number not matches for SSE Type Change, should be '$ncol_sse_type_change', the file has '$ncol
    exit
fi

ncol=`egrep -m 1 'Dynamic_merge' $fname.*bse*.[0-9]|wc -w`
if [[ $ncol -ne $ncol_bse_dyn_merge ]] && [[ $ncol -ne 0 ]]; then
    echo 'Error! column number not matches for BSE dynamical merger, should be '$ncol_bse_dyn_merge', the file has '$ncol
    exit
fi

ncol=`egrep -m 1 'SN_kick' $fname.*bse*.[0-9]|wc -w`
if [[ $ncol -ne $ncol_bse_sn_kick ]] && [[ $ncol -ne 0 ]]; then
    echo 'Error! column number not matches for BSE SN kick, should be '$ncol_bse_sn_kick', the file has '$ncol
    exit
fi

ncol=`egrep -v -m 1 '(SN_kick|Dynamic_merge)' $fname.*bse*.[0-9]|wc -w`
if [[ $ncol -ne $ncol_bse_type_change ]] && [[ $ncol -ne 0 ]]; then
    echo 'Error! column number not matches for BSE Type Change, should be '$ncol_bse_type_change', the file has '$ncol
    exit
fi


# clear each file
for ((i=0;i<$nsuffixes;i=i+1))
do
    s=${suffixes[$i]}
    tindex=${tindices[$i]}
    file=$fname.$s
    echo $file'; time column: '$tindex
    if [ $s == 'status' ]; then
	lst=$file
    elif [ ! -z $nmpi ]; then
	nend=`expr $nmpi - 1`
	lst=`seq 0 $nend|awk -v f=$file '{printf("$s.$d\n", f,$1)}`
    else
	lst=`ls |egrep $file'.[0-9]+$'`
    fi
    for f in $lst
    do
	echo 'process '$f
	if [ -e $f.bk ]; then
	    if [ -z $rmi ]; then
		rm -i $f.bk
		mv $f $f.bk
		echo 'backup '$f' to '$f.bk
	    fi
	else
	    mv $f $f.bk
	    echo 'backup '$f' to '$f.bk
	fi
	if [[ $s == *'sse'* ]]; then
	    awk -v t=$tcrit -v tsn=$tindex_sse_sn_kick -v ttch=$tindex_sse_type_change '{if (($1!="SN_kick" && $ttch<=t) || ($1=="SN_kick" && $tsn<=t)) print $LINE}' $f.bk >$f
	elif [[ $s == *'bse'* ]]; then
	    awk -v t=$tcrit -v tsn=$tindex_bse_sn_kick -v ttch=$tindex_bse_type_change -v tdyn=$tindex_bse_dyn_merge '{if ($1=="SN_kick") {if ($tsn<=t) print $LINE} else if ($1=="Dynamic_merge:") {if($tdyn<=t) print $LINE} else if ($ttch<=t) print $LINE;}' $f.bk >$f
	else
	    awk -v t=$tcrit -v ti=$tindex '{if ($ti<=t) print $LINE}' $f.bk >$f
	fi
    done
done

#



