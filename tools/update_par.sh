#!/bin/bash
# update old version of data to new version

unset mode
app_mode=1
se_mode=0
orb_mode=0
ehard_mode=1
agroup_mode=1
 
until [[ `echo x$1` == 'x' ]]
do
    case $1 in
 	-h) shift;
 	    echo 'Update input parameter from version before Oct 18, 2020 to the new version';
 	    echo 'Usage: petar.update.par [options] [input parameter filename]';
 	    echo 'Options:';
 	    echo '   -p: initial parameters of petar';
 	    echo '   -a: no append option (in versions before Aug 8, 2020 on Github; for -p)';
 	    echo '   -s: stellar evolution is used (for -p)';
 	    echo '   -o: orbit-samping is used (for -p)';
 	    echo '   -h: no hard-check-energy (for -p)';
 	    echo '   -g: no adjust group print (for -p)';
 	    echo 'PS: -o, -h, -g are used for non-default case (when users switch on/off special features in configure)';
 	    echo '   -b: initial parameters of bse';
 	    echo '   -t: initial parameters of galpy';
 	    exit;;
 	-p) mode='p'; shift;;
 	-b) mode='b'; shift;;
 	-t) mode='t'; shift;;
 	-a) app_mode=0; shift;;
 	-s) se_mode=1; shift;;
 	-o) orb_mode=1; shift;;
 	-h) ehard_mode=0; shift;;
 	-g) agroup_mode=0; shift;;
 	*) fname=$1; shift;;
    esac
done
 
if [ ! -e $fname ] | [ -z $fname ] ; then
    echo 'Error, file name not provided' 
    exit
fi
 
if  [ -z $mode ]; then
    echo 'Please provide which type of the input data is (-p, -b, -g)'
    exit
fi
 
# input par
if [ $mode == 'p' ]; then
    value=(`cat $fname`)
    ic=0
    echo 'F r-ratio '${value[$ic]} >$fname.new
    ic=`expr $ic + 1`
    echo 'F T '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'I number-step-tt '${value[$ic]%.*} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F t '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F hermite-eta '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F G '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F s '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F o '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F search-vel-factor '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F search-peri-factor '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F dt-max-factor '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    #echo 'F dt-error-pert '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    #echo 'F dt-error-iso '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F energy-err-ar '${value[$ic]} >>$fname.new
    if [ $ehard_mode == 1 ];then
	ic=`expr $ic + 1`
	echo 'F energy-err-hard '${value[$ic]} >>$fname.new
    fi
    ic=`expr $ic + 1`
    echo 'F soft-eps '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F r '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F r-bin '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F r-search-min '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F r-escape '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F slowdown-factor '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'I b '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'I n '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'I id-offset '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'I number-leaf-limit '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'I number-group-limit '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'I number-interrupt-limit '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'I number-sample-average '${value[$ic]} >>$fname.new
    if [[ $orb_mode == 1 ]];then
	ic=`expr $ic + 1`
	echo 'I number-split '${value[$ic]} >>$fname.new
    fi	
    ic=`expr $ic + 1`
    echo 'I u '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'I dt-min-hermite '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    #echo 'F dt-min-ar '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'I step-limit-ar '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'I i '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'I w '${value[$ic]} >>$fname.new
    if [ $se_mode == 1 ];then
	ic=`expr $ic + 1`
	echo 'I stellar-evolution '${value[$ic]} >>$fname.new
    fi
    ic=`expr $ic + 1`
    echo 'I detect-interrupt '${value[$ic]} >>$fname.new
    if [ $agroup_mode == 1 ]; then
	ic=`expr $ic + 1`
	echo 'I write-group-info '${value[$ic]} >>$fname.new
    fi
    if [ $app_mode == 1 ];then
	ic=`expr $ic + 1`
	echo 'I a '${value[$ic]} >>$fname.new
    fi
    ic=`expr $ic + 1`
    echo 'S f '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'S p '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'S snap-filename '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    [ ${#value[@]} != $ic ] && echo 'Error: number of values not match, read '${#value[@]}', should be '$ic
    echo 'generate new input parameter to '$fname.new
fi

# input par for bse
if [ $mode == 'b' ]; then
    value=(`cat $fname`)
    ic=0
    echo 'F bse-neta '${value[$ic]} >$fname.new
    ic=`expr $ic + 1`
    echo 'F bse-bwind '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F bse-hewind '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F bse-alpha '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F bse-lambda '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F bse-beta '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F bse-xi '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F bse-bhwacc '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F bse-epsnov '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F bse-eddfac '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F bse-gamma '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F bse-sigma '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F bse-pts1 '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F bse-pts2 '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F bse-pts3 '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F bse-tscale '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F bse-rscale '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F bse-mscale '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F bse-vscale '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F bse-metallicity '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'I bse-ceflag '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'I bse-tflag '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'I bse-wdflag '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'I bse-bhflag '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'I bse-nsflag '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'I bse-psflag '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'I bse-kmech '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'I bse-ecflag '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'I bse-idum '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    [ ${#value[@]} != $ic ] && echo 'Error: number of values not match, read '${#value[@]}', should be '$ic
    echo 'generate new input parameter to '$fname.new
fi

# input par for galpy
if [ $mode == 't' ]; then
    value=(`cat $fname`)
    ic=0
    echo 'F galpy-rscale '${value[$ic]} >$fname.new
    ic=`expr $ic + 1`
    echo 'F galpy-tscale '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F galpy-vscale '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F galpy-fscale '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    echo 'F galpy-pscale '${value[$ic]} >>$fname.new
    ic=`expr $ic + 1`
    if [[ ${value[$ic]} == 'None' ]]; then
	echo 'S galpy-type-arg __NONE__' >>$fname.new
    else
	echo 'S galpy-type-arg '${value[$ic]} >>$fname.new
    fi
    ic=`expr $ic + 1`
    if [[ ${value[$ic]} == 'None' ]]; then
	echo 'S galpy-set __NONE__' >>$fname.new
    else
	echo 'S galpy-set '${value[$ic]} >>$fname.new
    fi
    ic=`expr $ic + 1`
    if [[ ${value[$ic]} == 'None' ]]; then
	echo 'S galpy-conf-file __NONE__' >>$fname.new
    else
	echo 'S galpy-conf-file '${value[$ic]} >>$fname.new
    fi
    ic=`expr $ic + 1`
    [ ${#value[@]} != $ic ] && echo 'Error: number of values not match, read '${#value[@]}', should be '$ic
    echo 'generate new input parameter to '$fname.new
fi
