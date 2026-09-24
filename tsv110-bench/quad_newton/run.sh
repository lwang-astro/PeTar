#!/bin/bash
# run the quad Newton study and store raw CSV in data/
set -e
cd "$(dirname "$0")"
mkdir -p data
NACC=${NACC:-2000000}

echo "== rsinv accuracy/throughput =="
./rsinv cubic  $NACC > data/rsinv_cubic.csv
./rsinv newton $NACC > data/rsinv_newton.csv

echo "== SP quad force/potential errors =="
: > data/errors.csv
for cfg in "1 0" "10 0" "1 1000000"; do
  set -- $cfg; scale=$1; off=$2
  for seed in 1 2 3; do
    ./qn_cubic  err 1000 2000 $scale $off $seed >> data/errors.csv
    ./qn_newton err 1000 2000 $scale $off $seed >> data/errors.csv
  done
done
cat data/errors.csv

echo "== per-particle error dump (CDF input), 3 seeds merged =="
: > data/errdump_cubic.csv; : > data/errdump_newton.csv
for seed in 1 2 3; do
  ./qn_cubic  errdump 1000 2000 1 0 $seed >  data/errdump_cubic_$seed.csv
  ./qn_newton errdump 1000 2000 1 0 $seed >  data/errdump_newton_$seed.csv
done

echo "== quad kernel timing (3 scans per variant) =="
: > data/timing.csv
for scan in 1 2 3; do
  for ni in 4 16 64 256 1024; do
    for ns in 8 32 128 512 2048; do
      ./qn_cubic  time $ni $ns >> data/timing.csv
      ./qn_newton time $ni $ns >> data/timing.csv
    done
  done
done

echo "== mass>0 filtering check (EP-EP) =="
: > data/masszero.csv
./qn_cubic  masszero 1000 2000 7 >> data/masszero.csv
./qn_newton masszero 1000 2000 7 >> data/masszero.csv
cat data/masszero.csv
echo RUN_DONE
