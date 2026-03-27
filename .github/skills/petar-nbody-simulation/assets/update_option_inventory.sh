#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
ASSET_DIR="$ROOT_DIR/assets"
BIN_GLOB="/home/lwang/bin/petar.mpi.omp.avx512*"

mkdir -p "$ASSET_DIR"

# Remove stale inventory files from deleted binaries before regenerating.
rm -f "$ASSET_DIR"/petar.mpi.omp.avx512*.options.txt

extract_opts() {
    local bin="$1"
    local out="$2"
        if ! "$bin" -h 2>/dev/null \
      | awk '{for(i=1;i<=NF;i++){if($i ~ /^-[A-Za-z0-9]$/ || $i ~ /^--[A-Za-z0-9][A-Za-z0-9-]*$/) print $i}}' \
      | sed 's/[,:;]$//' \
            | sort -u > "$out"; then
                echo "Warning: failed to parse help from $(basename "$bin")" >&2
                : > "$out"
        fi
}

# 1) Extract per-binary options
for b in $BIN_GLOB; do
    [[ -x "$b" ]] || continue
    bn="$(basename "$b")"
    extract_opts "$b" "$ASSET_DIR/${bn}.options.txt"
done

# 2) Build matrix for solver binaries and list helper tools separately.
mapfile -t solver_files < <(ls "$ASSET_DIR"/petar.mpi.omp.avx512*.options.txt 2>/dev/null \
  | grep -v 'format.transfer' \
  | grep -v 'hard.debug' \
  | sort)

mapfile -t helper_files < <(ls "$ASSET_DIR"/petar.mpi.omp.avx512*.options.txt 2>/dev/null \
  | grep -E 'format.transfer|hard.debug' \
  | sort)

matrix="$ASSET_DIR/option-matrix.md"
{
    echo "# PeTar AVX512 Option Matrix (Auto-generated)"
    echo
    echo "This file is generated from <binary> -h option extraction."
    echo
    echo "## Solver Binaries Included"
    echo
    for f in "${solver_files[@]}"; do
        bn="$(basename "$f" .options.txt)"
        n="$(wc -l < "$f")"
        echo "- $bn ($n options)"
    done

    echo
    echo "## Helper Tool Binaries Excluded From Simulation Commands"
    echo
    for f in "${helper_files[@]}"; do
        bn="$(basename "$f" .options.txt)"
        n="$(wc -l < "$f")"
        echo "- $bn ($n options)"
    done

    echo
    echo "## Common Core Options"
    echo

    if [[ ${#solver_files[@]} -gt 0 ]]; then
        tmp_common="$(mktemp)"
        cp "${solver_files[0]}" "$tmp_common"
        for f in "${solver_files[@]:1}"; do
            comm -12 "$tmp_common" "$f" > "$tmp_common.next"
            mv "$tmp_common.next" "$tmp_common"
        done
        sed 's/^/- /' "$tmp_common"
        echo
        echo "## Feature-specific Options By Binary"
        echo
        for f in "${solver_files[@]}"; do
            bn="$(basename "$f" .options.txt)"
            echo "### $bn"
            echo
            extras="$(comm -23 "$f" "$tmp_common")"
            if [[ -n "$extras" ]]; then
                echo "$extras" | sed 's/^/- /'
            else
                echo "- (no extra options beyond core)"
            fi
            echo
        done
        rm -f "$tmp_common"
    fi
} > "$matrix"

echo "Updated option inventory in: $ASSET_DIR"
