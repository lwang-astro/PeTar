#!/usr/bin/env bash
set -u

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT_DIR"

SCENARIO="test/validation/scenarios/t2_binary_conservation_longterm.json"
REPORT="test/validation/out/report.t2.binary.json"
OUTDIR="test/validation/out_t2_binary"
FINAL_HTML="test/validation/out/t2_binary_conservation_summary.html"

BIN_NON64="${1:-/home/lwang/bin/petar.mpi.omp.avx2.kdkdk4}"

echo "[T2] non64b binary: $BIN_NON64"

echo "[T2] Run long-term binary conservation scenario..."
set +e
python3 test/validation/run_validation.py \
  --scenario "$SCENARIO" \
  --out-dir "$OUTDIR" \
  --report "$REPORT" \
  --var petar_bin_switch="$BIN_NON64"
RC=$?
set -e

echo "[T2] Build HTML summary..."
python3 test/validation/t2_binary_conservation_report.py \
  --report "$REPORT" \
  --output "$FINAL_HTML"

echo "[T2] Done."
echo "  - report json:  $REPORT (exit=$RC)"
echo "  - summary html: $FINAL_HTML"

if [[ $RC -eq 0 ]]; then
  exit 0
fi

exit 1
