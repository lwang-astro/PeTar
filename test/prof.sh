make -f ../Makefile nbody.out 
CPUPROFILE_FREQUENCY=1000 CPUPROFILE=./prof.out OMP_NUM_THREADS=6 ./nbody.out -i 1 -B 100 -t 10.0 input/plum1000b100.input 1>nbody.log 2>&1
pprof nbody.out ./prof.out 
