## plot semi-ecc-binary for PeTar snapshots in parallel
##-----------------------------##

read -p "Parameter group name (test)": logname
[[ -z $logname ]] && logname="test"
read -p "time suffix list (tlst):" dlist
[[ -z $dlist ]] && dlist=tlst
[ -e $dlist ] || (echo $dlist" not exist" && exit 1) || exit 

#read -p "output file with scale factor (output): " ofile
#[[ -z $ofile ]] && ofile=output
#scale $ofile
nb_tscale=1
nb_vscale=1
nb_mscale=1
nb_rscale=1

ntime=`date +%F-%H-%M-%S`
echo dlist $dlist>> semi_ecc_binary.$logname.log.$ntime
read -p "time scale to NB ($nb_tscale):" tscale
[[ -z $tscale ]] && tscale=$nb_tscale
echo tscale $tscale >> semi_ecc_binary.$logname.log.$ntime

read -p "distance scale to NB ($nb_rscale):" rscale
[[ -z $rscale ]] && rscale=$nb_rscale
echo rscale $rscale >> semi_ecc_binary.$logname.log.$ntime

read -p "velocit scale to NB ($nb_vscale):" vscale
[[ -z $vscale ]] && vscale=$nb_vscale
echo vscale $vscale >> semi_ecc_binary.$logname.log.$ntime

read -p "mass scale to NB ($nb_mscale):" mscale
[[ -z $mscale ]] && mscale=$nb_mscale
echo mscale $mscale >> semi_ecc_binary.$logname.log.$ntime

read -p "output figure name (semi_ecc_binary):" dout
[[ -z $dout ]] && dout=semi_ecc_binary
echo dout $dout >> semi_ecc_binary.$logname.log.$ntime

read -p "semi min (1e-7):" amin
[[ -z $amin ]] && amin=1e-7
echo amin $amin >> semi_ecc_binary.$logname.log.$ntime

read -p "semi max (5):" amax
[[ -z $amax ]] && amax=5
echo amax $amax >> semi_ecc_binary.$logname.log.$ntime

read -p "r min (1e-2):" rmin
[[ -z $rmin ]] && rmin=1e-2
echo rmin $rmin >> semi_ecc_binary.$logname.log.$ntime

read -p "r max (100):" rmax
[[ -z $rmax ]] && rmax=100
echo rmax $rmax >> semi_ecc_binary.$logname.log.$ntime

read -p "ekin min (1e-7):" ekmin
[[ -z $ekmin ]] && ekmin=1e-7
echo ekmin $ekmin >> semi_ecc_binary.$logname.log.$ntime

read -p "ekin max (1000):" ekmax
[[ -z $ekmax ]] && ekmax=1000
echo ekmax $ekmax >> semi_ecc_binary.$logname.log.$ntime

read -p "Number of parallel threads (4):" nloop
[[ -z $nloop ]] && nloop=4

ic=`wc -l $dlist|cut -d' ' -f1`
lnum=`expr $ic / $nloop`
[ `expr $lnum \* $nloop` -lt $ic ] && let lnum++
rm -f $dlist"_split."*
split -a 2 -d -l $lnum $dlist $dlist"_split."

lst=`ls $dlist"_split."*`
for kk in $lst
do
    rm -f semi-ecc-binary.${kk}
    python semi-ecc-binary-movie.py $kk $dout $tscale $rscale $vscale $mscale $amin $amax $rmin $rmax $ekmin $ekmax &> semi-ecc-binary.${kk##*.}.$logname &
done
