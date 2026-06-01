#!/usr/bin/env bash
set -u

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT_DIR"

SCENARIO="test/validation/scenarios/t3_binary_hard_switch_longterm.json"
REPORT="test/out/report.t3.binary.json"
OUTDIR="test/out/validation_t3_binary"
FINAL_HTML="test/out/t3_binary_hard_switch_summary.html"

BIN_NON64="${1:-/home/lwang/bin/petar.mpi.omp.avx2.kdkdk4}"

echo "[T3] non64b binary: $BIN_NON64"

echo "[T3] Run binary hard-switching scenario..."
set +e
python3 test/validation/run_validation.py \
  --scenario "$SCENARIO" \
  --out-dir "$OUTDIR" \
  --report "$REPORT" \
  --var petar_bin_switch="$BIN_NON64"
RC=$?
set -e

echo "[T3] Build HTML summary..."
python3 test/validation/t3_binary_hard_switch_report.py \
  --report "$REPORT" \
  --output "$FINAL_HTML"

echo "[T3] Done."
echo "  - report json:  $REPORT (exit=$RC)"
echo "  - summary html: $FINAL_HTML"

if [[ $RC -eq 0 ]]; then
  exit 0
fi

exit 1
