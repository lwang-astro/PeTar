source ~/.uftools.sh

flist=hard_dump.list

nlst=`wc -l $flist|coln c1`
for ((i=1;i<=$nlst;i=i+1)) 
do
    fline=`sed -n $i' p' $flist`
    fname=`echo $fline|cut -d':' -f1`
    if [ -e hard_dump.$fname ]; then
	finfo=`echo $fline|cut -d':' -f2`
	opt=`echo $fline|cut -d':' -f3`
	echo $fline
	flog=hard.debug.$fname.log
	fdata=hard.debug.$fname.data
	./hard_debug.out $opt hard_dump.$fname 1>$fdata 2>$flog
	tail -10 $flog |egrep 'Hard Energy:'
    fi
done
