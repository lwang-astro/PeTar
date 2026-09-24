#!/bin/bash
# N-scaling: base vs NEON on existing single-star inputs, t=2 Myr
ROOT=$HOME/tsv110-bench/e2e
BIN_BASE=$HOME/tests/PeTar-project/install/bin/petar.mpi.omp
BIN_NEON=$ROOT/PeTar-neon/build/petar.mpi.omp
SRC=$HOME/tests/PeTar-project/gravothermal1/scan

for N in 2000 5000 10000; do
  for name in base neon; do
    d=$ROOT/scal_${N}_${name}
    mkdir -p $d
    cp -f $SRC/N$N/input $d/input
    cd $d; rm -f data.* run.log
    bin=$BIN_BASE; [ "$name" = "neon" ] && bin=$BIN_NEON
    start=$(date +%s.%N)
    OMP_NUM_THREADS=1 OMP_STACKSIZE=128M mpirun -np 24 --bind-to none $bin -u 1 -b 0 -t 2 -o 1 -w 3 input > run.log 2>&1
    rc=$?
    end=$(date +%s.%N)
    w=$(awk "BEGIN{printf \"%.1f\", $end-$start}")
    ok=$(grep -c "successfully finished" run.log)
    steps=$(grep 'Tree step number' run.log | awk '{s+=$4} END{print s+0}')
    per=$(grep -A2 'Wallclock time per step' run.log | tail -1 | awk '{print $1}')
    tf=$(grep -A2 'Wallclock time per step' run.log | tail -1 | awk '{print $7}')
    tnb=$(grep -A2 'Wallclock time per step' run.log | tail -1 | awk '{print $6}')
    cf=$(grep -A2 'FDPS tree soft force time profile' run.log | tail -1 | awk '{print $9}')
    cfnb=$(grep -A2 'Tree neighbor time profile' run.log | tail -1 | awk '{print $9}')
    echo "SCAL N=$N $name wall=${w}s rc=$rc ok=$ok steps=$steps total_step=${per}s tree_force=${tf}s tree_nb=${tnb}s fdps_calcforce=${cf}s nb_calcforce=${cfnb}s"
  done
done
echo SCAL_DONE
