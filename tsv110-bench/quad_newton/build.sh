#!/bin/bash
# build the quad Newton study (AArch64 + NEON required)
set -e
# PETAR must point to the PeTar source tree; FDPS/SDAR default to its siblings.
PETAR=${PETAR:?please set PETAR to the PeTar source directory}
FDPS=${FDPS:-$PETAR/../FDPS}
SDAR=${SDAR:-$PETAR/../SDAR}
INC="-I$PETAR/src -I$FDPS/src -I$SDAR/src"
COMMON="-O3 -Wall -std=c++17 -mcpu=tsv110 -D USE_QUAD -D USE_NEON_KERNEL $INC -fopenmp"

gcc -O3 -mcpu=tsv110 -o rsinv rsinv_test.c -lm

g++ $COMMON -DNEON_QUAD_NEWTON=0 -o qn_cubic  quad_newton.cxx
g++ $COMMON -DNEON_QUAD_NEWTON=1 -o qn_newton quad_newton.cxx
echo built
