#!/bin/bash

unset run
unset pbin
unset dt_base
unset nmpi
unset nomp
unset fname
unset prefix_flag
unset opts
unset sfmt

until [[ `echo x$1` == 'x' ]]
do
    case $1 in
        -h) shift;
            echo 'A tool to find the suitable tree time step for a simulation.';
            echo '   By running simulations with the given snapshot file, find the tree time step that results in the best performance for the initial time steps.';
            echo 'Usage: petar.find.dt [options] [petar snapshot filename]';
            echo 'Options (default arguments shown in parentheses at the end):';
            echo '  -r [S] User-defined prefix before the petar commander (default: not used).'
            echo '         If mpiexec is not used for the MPI launcher or OMP_NUM_THREADS cannot set the thread number, users can manually define the prefix to set up MPI and thread numbers.';
            echo '         For example, if the syntax to use two MPI processors is "srun -N 2 petar ...", users can specify: -r "srun -N 2".';
	    echo '         Please make sure to enclose the prefix string in "".';
            echo '         If this option is used, -m and -o are suppressed.';
            echo '  -a [S] User-defined options used for the petar commander (default: not used).'
            echo '         For example, to set the binary number to 10 and use astronomical units, specify -a "-b 10 -u 1".';
	    echo '         Please make sure to enclose the options in "".';
            echo '         Note that the "-o", "-w", "-t", "-i", and "-s" options from the petar commander are defined here and cannot be used inside this -a block.';
            echo '  -p [S] petar commander name (default: petar)';
            echo '  -s [F] base tree step size (default: auto)';
            echo '  -m [I] number of MPI processors used for "mpiexec -n " (default: MPI is not used)';
            echo '  -o [I] number of OpenMP processors (default: auto)';
            echo '  -i [I] format of snapshot: 0 for BINARY, 1 for ASCII (default: 1)';
            echo 'PS: 1) Since the test is only based on the first 6 steps, the suggested tree step may not be the best for a long-term simulation.';
            echo '       Sometimes, the best tree step may result in a large changeover in radii, causing the wallclock time for the hard part to increase after some steps.';
            echo '       In such cases, the next smallest choice of tree step (0.5*best one) may be better.';
            echo '    2) For a new simulation, changeover radii (-r), minimum neighbor search radius (--r-search-min), and binary criterion (-r-bin) are automatically determined.';
            echo '       However, for a restarting simulation where input parameter files (input.par[.bse/.galpy]) exist, these parameters may need to be updated.';
            echo '       In this case, using -a "-p input.par -r 0 --r-bin 0 --r-search-min 0" can enable updating based on the new tree step size.';
            echo '    3) Files named check.perf.[time step].log and check.perf.test.log will be generated.';
            echo '       If the tool does not work correctly, checking these files can help identify the problems.';
	    exit;;
	-r) shift; run=$1; prefix_flag=yes; shift;;
	-a) shift; opts=$1; shift;;
	-p) shift; pbin=$1; shift;;
	-s) shift; dt_base=$1; shift;;
	-m) shift; nmpi=$1; shift;;
	-o) shift; nomp=$1; shift;;
	-i) shift; sfmt=$1; shift;;
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
[ -z $sfmt ] && sfmt=1
[ -z $prefix_flag ] || prefix=$run

echo 'commander: '$prefix' '$pbin' '$opts

if ! command -v $pbin &> /dev/null
then
    echo "Error, "$pbin" could not be found!"
    exit
fi

if [ -z $dt_base ]; then
    $prefix $pbin -w 0 $opts -i $sfmt -t 0.0 $fname &>check.perf.test.log
    dt_base=`egrep '^ dt_soft' check.perf.test.log |awk '{OFMT="%.14g"; print $3/4}'`
    echo 'Auto determine dt_base: '$dt_base
    #rm -f check.perf.test.log
else
    dt_base=`echo $dt_base|awk '{OFMT="%.14g"; print $1/2}'`
fi

tperf_pre=1e10
if [[ $sfmt == 1 ]]; then
    tzero=`head -1 $fname|awk '{print $3}'`
else
    petar.format.transfer -b -H -f $fname
    tzero=`head -1 $fname.H|awk '{print $3}'`
fi
echo 'time start: '$tzero
dt=$dt_base
dt_min=$dt
check_flag=true
success_flag=false
tcum=10000 # sec
while [[ $check_flag == true ]]
do
    dt=`echo $dt|awk '{OFMT="%.14g"; print $1*2.0}'`
    tend=`echo $dt |awk -v tzero=$tzero '{OFMT="%.14g"; print tzero+$1*6.01}' `
    echo 'test ending time '$tend' dt '$dt
    ${prefix} timeout $tcum $pbin $opts -w 0 -i $sfmt -t $tend  -s $dt -o $dt $fname &>check.perf.$dt.log
    egrep 'Wallclock' -A 3 check.perf.$dt.log |sed -n '/^\ *[0-9]/ p'|awk '{if (NR>1 && NR%2==0) print $1,$2+$3+$4+$8,$7+$9,$6+$10+$11,$12+$13}' >check.perf.$dt.tperf
    if [[ `wc -l check.perf.$dt.tperf|awk '{print $1}'` -ge 6 ]]; then 
	dt_reg=`egrep '^\ +dt_soft' check.perf.$dt.log|awk '{OFMT="%.14g"; print $3}'`
	tperf=(`awk -v dt=$dt_reg 'BEGIN{t=1e10; ns=1.0/dt; th=0; ts=0; tc=0; td=0;} {if (t>$1) {t=$1; th=$2; ts=$3; tc=$4; td=$5;}} END{print t*ns,th*ns,ts*ns,tc*ns,td*ns,ns}' check.perf.$dt.tperf`)
	de=`grep 'Slowdown:' check.perf.$dt.log|awk '{print $2}'`
	echo 'check tree step: '$dt_reg', wallclock time for one time unit: '${tperf[0]}'  hard: '${tperf[1]}'  soft: '${tperf[2]}'   clustering: '${tperf[3]}'   domain: '${tperf[4]}'   step number: '${tperf[5]}
	tcum=`awk 'BEGIN{t=0} {t=t+$1} END{if (t<1.0) {t=1.0}; print t*3}' check.perf.$dt.tperf`
	echo '  wallclock time first 6 steps: '`awk '{print $1}' check.perf.$dt.tperf`' next timeout check: '$tcum
	echo '  Slowdown relative energy error: '$de
	check_flag=`echo $tperf |awk -v pre=$tperf_pre '{if ($1<pre) print "true"; else print "false"}'`
	[[ $check_flag == true ]] && dt_min=$dt
	tperf_pre=$tperf
	success_flag=true
    else
	check_flag=false
    fi
    #rm -f check.perf.$dt.log
    #rm -f check.perf.$dt.tperf
done

#index_min=`echo $tperf_collect|awk '{tmin=1e10; imin=0; for (i=1;i<=NF;i++) {if (tmin>$i) {tmin=$i; imin=i}}; print imin}'`
#ds_min=`echo $dt_list|awk -v imin=$index_min '{print $imin}'`
if [[ $success_flag == true ]];then
    echo 'Best performance choice: tree step: '$dt_min
    echo ${prefix}' '$pbin' '$opts' -s '$dt_min' '$fname
else
    echo 'Fail to check the performance, make sure the options are correct'
    echo 'See check.perf.test.log and check.perf.'$dt'.log for details from the output of petar'
fi

