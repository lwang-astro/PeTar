#!/usr/bin/env bash
set -u

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT_DIR"

SCENARIO="test/validation/scenarios/t4_tree_hard_from_triple.json"
REPORT="test/out/report.t4.triple.json"
OUTDIR="test/out/validation_t4_triple"
FINAL_HTML="test/out/t4_triple_mode_summary.html"

BIN_NON64="${1:-petar}"

echo "[T4] non64b binary: $BIN_NON64"

echo "[T4] Run five-mode triple scenario..."
set +e
python3 test/validation/run_validation.py \
  --scenario "$SCENARIO" \
  --out-dir "$OUTDIR" \
  --report "$REPORT" \
  --var petar_bin_switch="$BIN_NON64"
RC=$?
set -e

echo "[T4] Build HTML summary..."
python3 test/validation/t4_triple_mode_report.py \
  --report "$REPORT" \
  --output "$FINAL_HTML"

echo "[T4] Done."
echo "  - report json:  $REPORT (exit=$RC)"
echo "  - summary html: $FINAL_HTML"

if [[ $RC -eq 0 ]]; then
  exit 0
fi

exit 1
