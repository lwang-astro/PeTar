#!/usr/bin/env bash
set -u

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT_DIR"

SCENARIO="test/validation/scenarios/t1_high_ecc_changeover.json"
REPORT_64="test/validation/out/report.t1.64b.overlay.json"
REPORT_NON64="test/validation/out/report.t1.non64b.overlay.json"
OUTDIR_64="test/validation/out64b"
OUTDIR_NON64="test/validation/outnon64b"
FINAL_HTML="test/validation/out/t1_kdkdk4_summary.html"

BIN_64="${1:-/home/lwang/bin/petar.mpi.omp.avx2.64b.kdkdk4}"
BIN_NON64="${2:-/home/lwang/bin/petar.mpi.omp.avx2.kdkdk4}"

echo "[T1] 64b binary:    $BIN_64"
echo "[T1] non64b binary: $BIN_NON64"

echo "[T1] Run 64b scenario..."
set +e
python3 test/validation/run_validation.py \
  --scenario "$SCENARIO" \
  --out-dir "$OUTDIR_64" \
  --report "$REPORT_64" \
  --var petar_bin_order4="$BIN_64"
RC_64=$?
set -e

echo "[T1] Run non64b scenario..."
set +e
python3 test/validation/run_validation.py \
  --scenario "$SCENARIO" \
  --out-dir "$OUTDIR_NON64" \
  --report "$REPORT_NON64" \
  --var petar_bin_order4="$BIN_NON64"
RC_NON64=$?
set -e

echo "[T1] Build merged summary HTML..."
python3 test/validation/t1_kdkdk4_report.py \
  --report "$REPORT_64" \
  --primary-tag 64b \
  --compare-report "$REPORT_NON64" \
  --compare-tag non64b \
  --output "$FINAL_HTML"

echo "[T1] Done."
echo "  - 64b report:    $REPORT_64 (exit=$RC_64)"
echo "  - non64b report: $REPORT_NON64 (exit=$RC_NON64)"
echo "  - summary html:  $FINAL_HTML"

if [[ $RC_64 -eq 0 && $RC_NON64 -eq 0 ]]; then
  exit 0
fi

exit 1
