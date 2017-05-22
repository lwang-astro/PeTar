if [ -f .check ]; then
    data=(`cat .check`)
    dsn=${data[0]}
    tendn=${data[1]}
    nloopn=${data[2]}
else
    dsn=1e-8
    tendn=0.00390625
    nloopn=1001
fi

read -p 'ds ('$dsn'): ' ds
read -p 'tend ('$tendn'): ' tend
read -p 'nloop ('$nloopn'): ' nloop
read -p 'itermax ('$iter'): ' iter

[[ -z $ds ]] && ds=$dsn
[[ -z $tend ]] && tend=$tendn
[[ -z $nloop ]] && nloop=$nloopn

[[ -z $iter ]] || iter='-I '$iter


echo $ds' '$tend' '$nloop' '$iter >.check

data=ARC_dump.dat
pdata=ARC_dump.par

./ARC_debug.out -s $ds -t $tend -n 1 $iter 2>log 
egrep 'Loop' log |head -1001 >tloop 
./ARC_debug.out -s `calc \* $ds 0.001` $iter -n 1000 2>log 1>tint
