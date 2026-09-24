#!/bin/bash
# end-to-end PeTar demo comparison: base / autovec / NEON, t=10 Myr, 24 MPI x 1 OMP
ROOT=$HOME/tsv110-bench/e2e
SRC=$HOME/tests/PeTar-project/demo/cluster
BIN_BASE=$HOME/tests/PeTar-project/install/bin/petar.mpi.omp
BIN_AUTO=$ROOT/PeTar-auto/build/petar.mpi.omp
BIN_NEON=$ROOT/PeTar-neon/build/petar.mpi.omp

run_one(){
  name=$1; bin=$2
  mkdir -p $ROOT/run_$name
  cp -f $SRC/input $ROOT/run_$name/input
  for f in $SRC/input.par*; do cp -f $f $ROOT/run_$name/ 2>/dev/null; done
  cd $ROOT/run_$name
  rm -f data.* run.log
  start=$(date +%s)
  OMP_NUM_THREADS=1 OMP_STACKSIZE=128M mpirun -np 24 --bind-to none $bin -u 1 -b 500 -t 10 -o 2 input > run.log 2>&1
  rc=$?
  end=$(date +%s)
  echo "RESULT $name rc=$rc wall_s=$((end-start))"
  tail -4 run.log
  ls data.10 2>/dev/null && wc -l data.10 || true
}

run_one base $BIN_BASE
run_one auto $BIN_AUTO
run_one neon $BIN_NEON
echo ALL_DONE
