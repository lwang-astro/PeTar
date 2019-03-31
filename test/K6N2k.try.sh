source ~/.uftools.sh
#lst=`powseq -s 6 -n 10 0.5`
#for i in $lst
#do
#    ./nbody.out -s $i -t 10.0 input/K6N2k.input >K6N2k.$i.log
#    egrep '^Time=' K6N2k.$i.log >K6N2k.$i.de
#done

lst='1 2 3 4 5'
for i in $lst
do
    ./nbody.out -t 10.0 --search-factor $i input/K6N2k.input >K6N2k.sf.$i.log
    egrep '^Time=' K6N2k.sf.$i.log >K6N2k.sf.$i.de
done
