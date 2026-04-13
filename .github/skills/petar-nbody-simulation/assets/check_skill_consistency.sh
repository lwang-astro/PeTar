#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../../../" && pwd)"
cd "$REPO_ROOT"

FAIL=0

check() {
  local name="$1"
  local pattern="$2"
  shift 2
  local files=("$@")

  echo "[CHECK] $name"
  if grep -nE "$pattern" "${files[@]}"; then
    echo "[PASS] $name"
  else
    echo "[FAIL] $name"
    FAIL=$((FAIL+1))
  fi
  echo
}

check "SKILL sources include latest samples" \
  'star_cluster_plummer_N1k_binaries\.sh|star_cluster_plummer_N1k_GalpyMWPot\.sh|star_cluster_plummer_N1k_AgamaMWPotHunter24\.sh' \
  .github/skills/petar-nbody-simulation/SKILL.md

check "SKILL contains external tools and map rule" \
  'petar\.galpy\.help|petar\.external\.galpy|petar\.external\.agama|petar\.external\.<galpy\|agama>' \
  .github/skills/petar-nbody-simulation/SKILL.md

check "Agama post-processing aligned in SKILL and default-postprocessing" \
  'petar\.data\.process -t agama --r-escape tidal -G 0\.00449830997959438|petar\.data\.process -i bse -t agama --r-escape tidal -G 0\.00449830997959438' \
  .github/skills/petar-nbody-simulation/SKILL.md \
  .github/skills/petar-nbody-simulation/assets/default-postprocessing.md

check "script-tools contains external map workflow" \
  'petar\.external\.galpy|petar\.external\.agama|pot_conf|petar\.external\.pot\.movie|petar\.movie --ext-pot' \
  .github/skills/petar-nbody-simulation/assets/script-tools.md

check "prompt-starters contains external visualization template" \
  '外势势场图与可视化链路|petar\.external\.agama|petar\.external\.pot\.movie|petar\.movie --ext-pot' \
  .github/skills/petar-nbody-simulation/assets/prompt-starters.md

if [[ "$FAIL" -ne 0 ]]; then
  echo "Consistency check finished with $FAIL failure(s)."
  exit 1
fi

echo "All consistency checks passed."
