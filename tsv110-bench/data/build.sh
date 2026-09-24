#!/bin/bash
# build the kernel benchmark variants (AArch64 + NEON required)
set -e
# PETAR must point to the PeTar source tree; FDPS/SDAR default to its siblings.
PETAR=${PETAR:?please set PETAR to the PeTar source directory}
FDPS=${FDPS:-$PETAR/../FDPS}
SDAR=${SDAR:-$PETAR/../SDAR}
INC="-I$PETAR/src -I$FDPS/src -I$SDAR/src"
DEF="-D USE_QUAD"
COMMON="-O3 -Wall -std=c++17 $DEF $INC -fopenmp"
echo "== bench_g (generic -O3) =="
g++ $COMMON -o bench_g bench_kernels.cxx
echo "== bench_a (-mcpu=tsv110 autovec) =="
g++ $COMMON -mcpu=tsv110 -o bench_a bench_kernels.cxx
echo "== bench_n (-mcpu=tsv110 + NEON) =="
g++ $COMMON -mcpu=tsv110 -D USE_NEON_KERNEL -o bench_n bench_kernels.cxx
echo built
