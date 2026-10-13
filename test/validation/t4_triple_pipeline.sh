#!/usr/bin/env bash
set -u

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT_DIR"

SCENARIO="test/validation/scenarios/t4_tree_hard_from_triple.json"
REPORT="test/out/report.t4.triple.json"
OUTDIR="test/out/validation_t4_triple"
FINAL_HTML="test/out/t4_triple_mode_summary.html"

# Rule: select the required family via petar.select; if not found or the
# installed binary version does not match the repository VERSION, reinstall.
ver_repo="$(head -c 16 VERSION | tr -d '[:space:]')"
sel_bin="$(petar.select --require bse --optional mpi,omp,avx2 >/dev/null 2>&1 && readlink -f "$(command -v petar)")"
need_install=0
if [[ -z "$sel_bin" || ! -x "$sel_bin" ]]; then
    need_install=1
else
    ver_bin="$("$sel_bin" -h 2>&1 | grep -m1 '^Version:' | awk '{print $2}' | cut -d_ -f1)"
    [[ "$ver_bin" != "$ver_repo" ]] && need_install=1
fi
if [[ "$need_install" -eq 1 ]]; then
    echo "[T4] bse-family binary missing or stale (repo VERSION=$ver_repo); reinstalling (--with-interrupt=bse)..."
    ./configure --with-interrupt=bse
    make -j"$(nproc)" install
    petar.select --require bse --optional mpi,omp,avx2
    sel_bin="$(readlink -f "$(command -v petar)")"
fi

BIN_NON64="${1:-$sel_bin}"

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
