#!/bin/bash

suffixes='esc sse bse mosse mobse sseEmp bseEmp'
unset rmi
unset onlylist
unset groupflag

until [[ `echo x$1` == 'x' ]]
do
    case $1 in
	-h) shift;
	    echo 'A tool for organizing petar output files (for version before 1708e).';
	    echo 'Functionality:';
	    echo '    1) Combine separated output data from multiple MPI processes with filename suffixes: '$suffixes' group';
	    echo '    2) Split legacy SSE/BSE output files into different files with suffixes "type_change", "sn_kick", "gw_kick" and "dynamic_merge".';
	    echo '    3) If the option "-g" is used, combine group files from mutliple MPI processes:';
	    echo '         [prefix].group.[rank].n[N] -> [output].group.n[N]';
	    echo '       (works for both ASCII and BINARY group files)';
	    echo 'Usage: petar.data.gether [options] [data filename prefix]';
	    echo '       The data filename prefix is defined by "petar -f"; the default case is "data".';
            echo 'Options (default arguments shown in parentheses at the end):';
	    echo '  -f [S] Output filename prefix (default: [data filename prefix])';
	    echo '  -n [I] Number of MPI processes (default: auto-detect)';
	    echo '  -i     Ask before removing existing combined files (default: no ask)';
	    echo '  -l     Only generate a list of snapshot data files';
	    echo '  -g     Combine group files (can be slow if file sizes are large and there are many MPI processes)';
	    exit;;
	-f) shift; fout=$1; shift;;
	-n) shift; nmpi=$1; shift;;
	-i) rmi=1; shift;;
	-l) onlylist=1; shift;;
	-g) groupflag=1; shift;;
	*) fname=$1;shift;;
    esac
done

if [ ! -e $fname ] | [ -z $fname ] ; then
    echo 'Error, file name not provided'
    exit
fi
[ -z $fout ] && fout=$fname

echo 'data filename prefix: '$fout

if [ ! -e $fout.snap.lst ]; then
	flen=`expr ${#fname} + 2`
	ls|egrep '^'$fname'.[0-9]+$' |sort -n -k 1.${flen} >$fout.snap.lst
fi

[ ! -z $onlylist ] && exit

for s in $suffixes
do
    file=$fname.$s
    echo $file
    if [ -e $file.0 ]; then
	if [ -e $fout.$s ]; then
	    if [ -z $rmi ]; then
		rm -f $fout.$s
	    else
		rm -i $fout.$s
	    fi
	fi

	echo 'gether '$file'.* to '$fout.$s
	if [ ! -z $nmpi ]; then
	    nend=`expr $nmpi - 1`
	    lst=`seq 0 $nend`
	    for i in $lst
	    do
		cat $file.$i >>$fout.$s
	    done
	else
	    flist=`ls |egrep $file'.[0-9]+$'`
	    echo $flist
	    cat $flist >$fout.$s
	fi
    fi
done

if [ ! -z $groupflag ]; then
	# New-format group files: [prefix].group.[rank].n[N]
	nlist=`ls | egrep '^'$fname'.group.[0-9]+.n[0-9]+$' | sed 's/.*\.n/n/' | sort -u`
	for ns in $nlist
	do
		out_group=$fout.group.$ns
		if [ -e $out_group ]; then
			if [ -z $rmi ]; then
				rm -f $out_group
			else
				rm -i $out_group
			fi
		fi

		echo 'gether '$fname'.group.*.'$ns' to '$out_group
		if [ ! -z $nmpi ]; then
			nend=`expr $nmpi - 1`
			lst=`seq 0 $nend`
			for i in $lst
			do
				in_group=$fname.group.$i.$ns
				[ -e $in_group ] && cat $in_group >>$out_group
			done
		else
			flist=`ls | egrep '^'$fname'.group.[0-9]+\.'$ns'$'`
			[ ! -z "$flist" ] && cat $flist >$out_group
		fi
	done

fi

sse_opt='.sse .mosse .sseEmp'
for s in $sse_opt
do
    if [ -e $fout$s ]; then
	echo 'get '$s' type_change, sn_kick'
	egrep '^Type_change ' $fout$s |sed 's/Type_change//g' >$fout$s.type_change
	egrep '^SN_kick ' $fout$s |sed 's/SN_kick//g' >$fout$s.sn_kick
    fi
done

bse_opt='.bse .mobse .bseEmp'
for s in $bse_opt
do
    if [ -e $fout$s ]; then
	echo 'get '$s' type_change, sn_kick, gw_kick, dynamic_merge, GW_tide_merge, hyperbolic_tde, binary_tde, tide'
	egrep '^Dynamic_merge' $fout$s |sed 's/Dynamic_merge://g' >$fout$s.dynamic_merge.tmp
	awk '{if (NF==45) {for (i=1; i<=5; i++) printf("%s ", $i); printf("0 0 0 "); for (i=6; i<=NF; i++) printf("%s ", $i); printf("\n");} else print $LINE}' $fout$s.dynamic_merge.tmp > $fout$s.dynamic_merge
	rm -f $fout$s.dynamic_merge.tmp
	egrep '^Binary_merge' $fout$s |sed 's/Binary_merge://g' >$fout$s.binary_merge.tmp
	awk '{if (NF==45) {for (i=1; i<=5; i++) printf("%s ", $i); printf("0 0 0 "); for (i=6; i<=NF; i++) printf("%s ", $i); printf("\n");} else print $LINE}' $fout$s.binary_merge.tmp > $fout$s.binary_merge
	rm -f $fout$s.binary_merge.tmp
	egrep '^Hyperbolic_TDE' $fout$s |sed 's/Hyperbolic_TDE://g' >$fout$s.hyperbolic_tde.tmp
	awk '{if (NF==45) {for (i=1; i<=5; i++) printf("%s ", $i); printf("0 0 0 "); for (i=6; i<=NF; i++) printf("%s ", $i); printf("\n");} else print $LINE}' $fout$s.hyperbolic_tde.tmp > $fout$s.hyperbolic_tde
	rm -f $fout$s.hyperbolic_tde.tmp
	egrep '^Binary_TDE' $fout$s |sed 's/Binary_TDE://g' >$fout$s.binary_tde.tmp
	awk '{if (NF==45) {for (i=1; i<=5; i++) printf("%s ", $i); printf("0 0 0 "); for (i=6; i<=NF; i++) printf("%s ", $i); printf("\n");} else print $LINE}' $fout$s.binary_tde.tmp > $fout$s.binary_tde
	rm -f $fout$s.binary_tde.tmp
	egrep '^SN_kick' $fout$s |sed 's/SN_kick//g' >$fout$s.sn_kick
	egrep '^Hyperbolic_TDE' $fout$s |sed 's/Hyperbolic_TDE://g' >$fout$s.hyperbolic_tde
	egrep '^Tide' $fout$s |sed 's/Tide//g' >$fout$s.tide
	egrep '^GW_kick' $fout$s |sed 's/GW_kick//g' >$fout$s.gw_kick
	egrep -v '^(Dynamic_merge|GW_tide_merge|Hyperbolic_TDE|Binary_TDE|SN_kick|Tide|GW_kick)' $fout$s |awk '{for (i=2;i<=NF;i++) printf("%s ", $i); printf("\n")}' >$fout$s.type_change
    fi
done
