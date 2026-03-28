#!/usr/bin/env bash
set -euo pipefail

# Reusable PeTar launch template based on particle scale and primordial-binary richness.
# Usage:
#   N_SCALE=1e3 BIN_RICHNESS=none PETAR_BIN=petar \
#   PETAR_OPTS='-u 1 -t 10 -o 1' INPUT_SNAPSHOT=input OUTPUT_LOG=output \
#   bash .github/skills/petar-nbody-simulation/assets/launch-template-by-scale.sh
#
# Required env vars:
#   N_SCALE         : 1e3 | 1e4 | larger
#   BIN_RICHNESS    : none | few | many
#   PETAR_BIN       : solver binary (e.g. petar, petar.mpi.omp.avx2.bse)
#   PETAR_OPTS      : run options except snapshot filename
#   INPUT_SNAPSHOT  : input snapshot path
#   OUTPUT_LOG      : output log filename
# Optional env vars:
#   USE_MPI         : 0 | 1 (default 0)
#   OMP_THREADS     : force OMP threads (override heuristic)
#   MPI_RANKS       : force MPI ranks (override heuristic)
#   OUTPUT_PREFIX   : PeTar output prefix (-f). Auto-selected when empty.
#   RESTART_MODE    : 0 | 1 (default 0). If 1, avoid silently reusing default `data`.

N_SCALE=${N_SCALE:-1e3}
BIN_RICHNESS=${BIN_RICHNESS:-none}
PETAR_BIN=${PETAR_BIN:-petar}
PETAR_OPTS=${PETAR_OPTS:-"-u 1 -t 10 -o 1"}
INPUT_SNAPSHOT=${INPUT_SNAPSHOT:-input}
OUTPUT_LOG=${OUTPUT_LOG:-output}
USE_MPI=${USE_MPI:-0}
OUTPUT_PREFIX=${OUTPUT_PREFIX:-}
RESTART_MODE=${RESTART_MODE:-0}

heuristic_threads=1
heuristic_ranks=1

case "${N_SCALE}:${BIN_RICHNESS}" in
  1e3:none|1e3:few)
    heuristic_threads=1
    heuristic_ranks=1
    ;;
  1e3:many)
    heuristic_threads=4
    heuristic_ranks=1
    ;;
  1e4:none|1e4:few|1e4:many)
    heuristic_threads=4
    heuristic_ranks=1
    ;;
  larger:none|larger:few|larger:many)
    heuristic_threads=8
    heuristic_ranks=2
    ;;
  *)
    echo "Unsupported combination: N_SCALE='${N_SCALE}', BIN_RICHNESS='${BIN_RICHNESS}'" >&2
    exit 2
    ;;
esac

OMP_NUM_THREADS=${OMP_THREADS:-${heuristic_threads}}
MPI_RANKS=${MPI_RANKS:-${heuristic_ranks}}

if [[ -z "${OUTPUT_PREFIX}" ]]; then
  if [[ "${RESTART_MODE}" == "0" ]] && ! compgen -G "data*" > /dev/null; then
    OUTPUT_PREFIX="data"
  else
    OUTPUT_PREFIX="data.$(date +%Y%m%d_%H%M%S)"
  fi
fi

if [[ " ${PETAR_OPTS} " != *" -f "* ]]; then
  PETAR_OPTS="${PETAR_OPTS} -f ${OUTPUT_PREFIX}"
fi

echo "[PeTar launch template]"
echo "  N_SCALE=${N_SCALE}"
echo "  BIN_RICHNESS=${BIN_RICHNESS}"
echo "  PETAR_BIN=${PETAR_BIN}"
echo "  OMP_NUM_THREADS=${OMP_NUM_THREADS}"
echo "  MPI_RANKS=${MPI_RANKS}"
echo "  OUTPUT_PREFIX=${OUTPUT_PREFIX}"

echo "[run]"
if [[ "${USE_MPI}" == "1" ]]; then
  OMP_STACKSIZE=128M OMP_NUM_THREADS=${OMP_NUM_THREADS} \
    mpiexec -n ${MPI_RANKS} "${PETAR_BIN}" ${PETAR_OPTS} "${INPUT_SNAPSHOT}" &> "${OUTPUT_LOG}"
else
  OMP_STACKSIZE=128M OMP_NUM_THREADS=${OMP_NUM_THREADS} \
    "${PETAR_BIN}" ${PETAR_OPTS} "${INPUT_SNAPSHOT}" &> "${OUTPUT_LOG}"
fi

echo "Done: ${OUTPUT_LOG}"
