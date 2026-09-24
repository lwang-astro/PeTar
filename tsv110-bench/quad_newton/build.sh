#!/bin/bash
# build the quad Newton study (NEON/AArch64 required)
set -e
PETAR=${PETAR:-$HOME/tests/PeTar-project/PeTar}
FDPS=${FDPS:-$HOME/tests/PeTar-project/FDPS}
SDAR=${SDAR:-$HOME/tests/PeTar-project/SDAR}
INC="-I$PETAR/src -I$FDPS/src -I$SDAR/src"
COMMON="-O3 -Wall -std=c++17 -mcpu=tsv110 -D USE_QUAD -D USE_NEON_KERNEL $INC -fopenmp"

gcc -O3 -mcpu=tsv110 -o rsinv rsinv_test.c -lm

g++ $COMMON -DNEON_QUAD_NEWTON=0 -o qn_cubic  quad_newton.cxx
g++ $COMMON -DNEON_QUAD_NEWTON=1 -o qn_newton quad_newton.cxx
echo built
