#!/bin/bash
# repeated E2E timing: 3 reps each, t=10 Myr
ROOT=$HOME/tsv110-bench/e2e
SRC=$HOME/tests/PeTar-project/demo/cluster
BIN_BASE=$HOME/tests/PeTar-project/install/bin/petar.mpi.omp
BIN_AUTO=$ROOT/PeTar-auto/build/petar.mpi.omp
BIN_NEON=$ROOT/PeTar-neon/build/petar.mpi.omp

run_one(){
  name=$1; bin=$2; rep=$3
  d=$ROOT/tm_$name
  mkdir -p $d
  cp -f $SRC/input $d/input
  for f in $SRC/input.par*; do cp -f $f $d/ 2>/dev/null; done
  cd $d
  rm -f data.* run.log
  start=$(date +%s.%N)
  OMP_NUM_THREADS=1 OMP_STACKSIZE=128M mpirun -np 24 --bind-to none $bin -u 1 -b 500 -t 10 -o 2 input > run.log 2>&1
  rc=$?
  end=$(date +%s.%N)
  w=$(awk "BEGIN{printf \"%.1f\", $end-$start}")
  ok=$(grep -c "successfully finished" run.log)
  echo "RESULT $name rep$rep rc=$rc wall_s=$w finished=$ok"
}

for r in 1 2 3; do
  run_one base $BIN_BASE $r
  run_one auto $BIN_AUTO $r
  run_one neon $BIN_NEON $r
done
echo ALL_DONE
