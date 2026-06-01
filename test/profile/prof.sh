make -f ../Makefile nbody.out
for k in 6 5 4 3 2 1
#for k in 6
do
    for i in 1 10 100
#    for i in 100
    do
	echo 'Now OMP='$k' N='$i'00'
	rm -f data.*
#	echo OMP_STACKSIZE=100000 CPUPROFILE_FREQUENCY=500 CPUPROFILE=./prof.out.${i}k.omp${k} OMP_NUM_THREADS=$k ./nbody.out -i 2 -B ${i}00 -t 0.0625 -R 0.006 input/p${i}kb10.input 1>nbody.log.${i}k.omp${k} 2>&1
	OMP_STACKSIZE=100000 CPUPROFILE_FREQUENCY=500 CPUPROFILE=./prof.out.${i}k.omp${k} OMP_NUM_THREADS=$k ./nbody.out -i 2 -B ${i}00  -t 0.0625 -R 0.006 input/p${i}kb10.input 1>nbody.log.${i}k.omp${k} 2>&1
	ls data.*|sort -n -k1.6 > datalist
	./getdata -i 3 -f datalist 1>getdata.log 2>&1
	mv global.dat result/nbody.dat.p${i}kb10.omp${k}
    done
done

#pprof nbody.out ./prof.out 
