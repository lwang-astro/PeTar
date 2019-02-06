source ~/.uftools.sh

flist=hard_dump.list

nlst=`wc -l $flist|coln c1`
for ((i=1;i<=$nlst;i=i+1)) 
do
    fline=`sed -n $i' p' $flist`
    fname=`echo $fline|cut -d':' -f1`
    finfo=`echo $fline|cut -d':' -f2`
    opt=`echo $fline|cut -d':' -f3`
    echo $fline
    flog=hard.debug.$fname.log
    fdata=hard.debug.$fname.data
    ./hard_debug.out $opt hard_dump.$fname &>$flog
    mv hard_debug_ptcl $fdata
    tail -10 $flog |egrep 'Tot Energy:'
done
