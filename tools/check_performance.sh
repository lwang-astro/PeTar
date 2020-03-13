#!/bin/bash

until [[ `echo x$1` == 'x' ]]
do
    case $1 in
	-h) shift;
	    echo 'Check the tree soft step for best performance, ';
	    echo 'Usage: petar.init [options] [data file name]';
	    echo 'Options:';
	    echo '  -p: petar commander name (default: petar)';
	    echo '  -s: base tree step size (default: auto)';
	    echo '  -b: binary number (default: 0)';
	    echo '  -m: number of MPI processors (default: 1)';
	    echo '  -o: number of OpenMP processors (default: auto)';
	    exit;;
	-p) shift; pbin=$1; shift;;
	-s) shift; dt_base=$1; shift;;
	-b) shift; bnum=$1; shift;;
	-m) shift; nmpi=$1; shift;;
	-o) shift; nomp=$1; shift;;
	*) fname=$1;shift;;
    esac
done

if [ ! -e $fname ] | [ -z $fname ] ; then
    echo 'Error, file name not provided' 
    exit
fi

[ -z $pbin ] && pbin=petar
[ -z $nomp ] || prefix='env OMP_NUM_THREADS='$nomp
[ -z $nmpi ] || prefix=$prefix' mpiexec -n '$nmpi

if [ -z $dt_base ]; then
    $prefix $pbin -w 0 -t 0.0 $fname &>.check.perf.test.log
    dt_base=`egrep dt_soft .check.perf.test.log |awk '{OFMT="%.14g"; print $3/8}'`
    rm -f .check.perf.test.log
else
    dt_base=`echo $dt_base|awk '{OFMT="%.14g"; print $1/2}'`
fi

tperf_pre=1e10
dt=$dt_base
dt_min=$dt
check_flag=true
while [[ $check_flag == true ]]
do
    dt=`echo $dt|awk '{OFMT="%.14g"; print $1*2.0}'`
    tend=`echo $dt |awk '{print $1*6.01}' `
    ${prefix} $pbin -w 0 -t $tend -s $dt -o $dt $fname &>.check.perf.$dt.log
    tperf_list=`egrep 'Wallclock' -A 3 .check.perf.$dt.log |sed -n '/^\ *[0-9]/ p'|awk '{if (NR>1 && NR%2==0) print $1}'`
    #tperf=`echo $tperf_list|awk -v dt=$dt 'BEGIN{t=1e10;ns=1.0/dt} {if (NR>1 && NR%2==0) t=(t<$1)?t:$1} END{print t*ns}'`
    tperf=`echo $tperf_list|awk -v dt=$dt '{t=1e10;ns=1.0/dt; for(i=1;i<=NF;i++) t=(t<$i)?t:$i; print t*ns}'`
    de=`grep 'dE(SD)' .check.perf.$dt.log|awk '{print $11}'`
    echo 'check tree step: '$dt', wallclock time for one time unit: '$tperf
    echo '  wallclock time first 6 steps: '$tperf_list
    echo '  cumulative error: '$de
    check_flag=`echo $tperf |awk -v pre=$tperf_pre '{if ($1<pre) print "true"; else print "false"}'`
    [[ $check_flag == true ]] && dt_min=$dt
    tperf_pre=$tperf
    rm -f .check.perf.$dt.log
done

#index_min=`echo $tperf_collect|awk '{tmin=1e10; imin=0; for (i=1;i<=NF;i++) {if (tmin>$i) {tmin=$i; imin=i}}; print imin}'`
#ds_min=`echo $dt_list|awk -v imin=$index_min '{print $imin}'`
echo 'Best performance choice: tree step: '$dt_min
echo ${prefix}' '$pbin' -s '$dt_min' '$fname
