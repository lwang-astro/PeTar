#!/bin/bash

unset run
unset pbin
unset dt_base
unset bnum
unset nmpi
unset nomp
unset fname
unset prefix_flag
unset opts
unset unit
unset G

until [[ `echo x$1` == 'x' ]]
do
    case $1 in
	-h) shift;
	    echo 'Check the tree soft step for best performance, ';
	    echo 'Usage: petar.init [options] [data file name]';
	    echo 'Options:';
	    echo '  -r: user defined prefix before program name, please use " " to enclose the options.'
            echo '      If used, -m and -o are suppressed (default: "")'
	    echo '  -a: user defined options used for petar commander, please use " " to enclose the options.'
	    echo '      For example, "-b [binary number], -u [unit set] -G [gravitaitonal constant]" (default: "")'
	    echo '  -p: petar commander name (default: petar)';
	    echo '  -s: base tree step size (default: auto)';
#	    echo '  -b: binary number (default: 0)';
	    echo '  -m: number of MPI processors (default: 1)';
	    echo '  -o: number of OpenMP processors (default: auto)';
#	    echo '  -u: petar unit set (default: 0)';
#	    echo '  -G: gravitational constant (default: 1)'
	    exit;;
	-r) shift; run=$1; prefix_flag=yes; shift;;
	-a) shift; opts=$1; shift;;
	-p) shift; pbin=$1; shift;;
	-s) shift; dt_base=$1; shift;;
#	-b) shift; bnum=$1; shift;;
	-m) shift; nmpi=$1; shift;;
	-o) shift; nomp=$1; shift;;
#	-u) shift; unit=$1; shift;;
#	-G) shift; G=$1; shift;;
	*) fname=$1;shift;;
    esac
done

if [ ! -e $fname ] | [ -z $fname ] ; then
    echo 'Error, file name not provided' 
    exit
fi

#[ -z $opts ] && opts=''
[ -z $pbin ] && pbin=petar
[ -z $nomp ] || prefix='env OMP_NUM_THREADS='$nomp
[ -z $nmpi ] || prefix=$prefix' mpiexec -n '$nmpi
#[ -z $bnum ] || opts=$opts' -b '$bnum
[ -z $prefix_flag ] || prefix=$run
#[ -z $unit ] || opts=$opts' -u '$unit
#[ -z $G ] || opts=$opts' -G '$G

echo 'commander: '$prefix' '$pbin' '$opts

if [ -z $dt_base ]; then
    $prefix $pbin -w 0 $opts -t 0.0 $fname &>.check.perf.test.log
    dt_base=`egrep dt_soft .check.perf.test.log |awk '{OFMT="%.14g"; print $3/4}'`
    echo 'Auto determine dt_base: '$dt_base
    rm -f .check.perf.test.log
else
    dt_base=`echo $dt_base|awk '{OFMT="%.14g"; print $1/2}'`
fi

tperf_pre=1e10
tzero=`head -1 $fname|awk '{print $3}'`
dt=$dt_base
dt_min=$dt
check_flag=true
tcum=10000 # sec
while [[ $check_flag == true ]]
do
    dt=`echo $dt|awk '{OFMT="%.14g"; print $1*2.0}'`
    tend=`echo $dt |awk -v tzero=$tzero '{print tzero+$1*6.01}' `
    ${prefix} timeout $tcum $pbin -w 0 -t $tend $opts  -s $dt -o $dt $fname &>.check.perf.$dt.log
    egrep 'Wallclock' -A 3 .check.perf.$dt.log |sed -n '/^\ *[0-9]/ p'|awk '{if (NR>1 && NR%2==0) print $1,$2+$3+$4+$8,$7+$9,$6+$10+$11,$12+$13}' >.check.perf.$dt.tperf
    if [[ `wc -l .check.perf.$dt.tperf|awk '{print $1}'` -ge 6 ]]; then 
	dt_reg=`egrep dt_soft .check.perf.$dt.log|awk '{OFMT="%.14g"; print $3}'`
	tperf=(`awk -v dt=$dt_reg 'BEGIN{t=1e10; ns=1.0/dt; th=0; ts=0; tc=0; td=0;} {if (t>$1) {t=$1; th=$2; ts=$3; tc=$4; td=$5;}} END{print t*ns,th*ns,ts*ns,tc*ns,td*ns,ns}' .check.perf.$dt.tperf`)
	de=`grep 'dE(SD)' .check.perf.$dt.log|awk '{print $11}'`
	echo 'check tree step: '$dt_reg', wallclock time for one time unit: '${tperf[0]}'  hard: '${tperf[1]}'  soft: '${tperf[2]}'   clustering: '${tperf[3]}'   domain: '${tperf[4]}'   step number: '${tperf[5]}
	tcum=`awk 'BEGIN{t=0} {t=t+$1} END{if (t<1) {t=1}; print t*3}' .check.perf.$dt.tperf`
	echo '  wallclock time first 6 steps: '`awk '{print $1}' .check.perf.$dt.tperf`' sum: '$tcum
	echo '  cumulative error: '$de
	check_flag=`echo $tperf |awk -v pre=$tperf_pre '{if ($1<pre) print "true"; else print "false"}'`
	[[ $check_flag == true ]] && dt_min=$dt
	tperf_pre=$tperf
    else
	check_flag=false
    fi
    #rm -f .check.perf.$dt.log
    #rm -f .check.perf.$dt.tperf
done

#index_min=`echo $tperf_collect|awk '{tmin=1e10; imin=0; for (i=1;i<=NF;i++) {if (tmin>$i) {tmin=$i; imin=i}}; print imin}'`
#ds_min=`echo $dt_list|awk -v imin=$index_min '{print $imin}'`
echo 'Best performance choice: tree step: '$dt_min
echo ${prefix}' '$pbin' '$opts' -s '$dt_min' '$fname
