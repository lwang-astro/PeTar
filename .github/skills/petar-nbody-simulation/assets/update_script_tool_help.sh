#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
OUT_DIR="$ROOT_DIR/assets/tool-help"
mkdir -p "$OUT_DIR"

TOOLS=(
  petar.init
  petar.find.dt
  petar.update.par
  petar.data.clear
  petar.data.process
  petar.movie
  petar.data.gether
  petar.get.object.snap
  petar.format.transfer.post
  petar.galev.process
  petar.external.pot.movie
  petar.galpy.pot.movie
  petar.get.init.binary
)

rm -f "$OUT_DIR"/*.help.txt

for tool in "${TOOLS[@]}"; do
  out="$OUT_DIR/${tool}.help.txt"
  if "$tool" -h > "$out" 2>&1; then
    echo "$tool OK"
  else
    echo "$tool NONZERO"
  fi
done
