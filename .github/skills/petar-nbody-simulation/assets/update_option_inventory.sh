#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------------------
# Machine-specific binary option inventory.
#
# This script scans PATH for installed PeTar binaries and runs <binary> -h
# on each to produce a snapshot of the currently installed binary families
# and their options. This is a MACHINE-SPECIFIC snapshot.
#
# For the complete theoretical option space across all configure variants,
# use generate_option_reference.py instead, which parses source headers
# without requiring any binary installation.
#
# See also: option-reference.md  (source-generated, complete)
#           option-matrix.md     (this script's output, machine snapshot)
# ---------------------------------------------------------------------------

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
ASSET_DIR="$ROOT_DIR/assets"

mkdir -p "$ASSET_DIR"

# Remove stale inventory files before regenerating.
rm -f "$ASSET_DIR"/petar*.options.txt

declare -a SCRIPT_TOOL_NAMES=(
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

discover_petar_bins() {
    local dir
    local bin
    declare -A seen=()
    IFS=':' read -r -a path_dirs <<< "${PATH:-}"

    for dir in "${path_dirs[@]}"; do
        [[ -d "$dir" ]] || continue
        shopt -s nullglob
        for bin in "$dir"/petar*; do
            [[ -x "$bin" && -f "$bin" ]] || continue
            real_bin="$(readlink -f "$bin" 2>/dev/null || echo "$bin")"
            [[ -z "${seen[$real_bin]+x}" ]] || continue
            seen[$real_bin]=1
            echo "$real_bin"
        done
        shopt -u nullglob
    done | sort -u
}

is_script_tool() {
    local bn="$1"
    local tool
    for tool in "${SCRIPT_TOOL_NAMES[@]}"; do
        [[ "$bn" == "$tool" ]] && return 0
    done
    return 1
}

is_helper_tool() {
    local bn="$1"
    [[ "$bn" == *".hard.debug"* || "$bn" == *".format.transfer"* ]]
}

is_solver_by_options() {
    local opt_file="$1"
    grep -qx -- '-u' "$opt_file" && grep -qx -- '-t' "$opt_file" && grep -qx -- '-o' "$opt_file"
}

solver_perf_key() {
    local bn="$1"
    local opt_count="$2"
    local has_gpu=0
    local has_avx512=0
    local has_avx2=0
    local has_omp=0
    local has_mpi=0

    [[ "$bn" == *".gpu"* ]] && has_gpu=1
    [[ "$bn" == *".avx512"* ]] && has_avx512=1
    [[ "$bn" == *".avx2"* ]] && has_avx2=1
    [[ "$bn" == *".omp"* ]] && has_omp=1
    [[ "$bn" == *".mpi"* ]] && has_mpi=1

    printf "%d%d%d%d%d|%04d" "$has_gpu" "$has_avx512" "$has_avx2" "$has_omp" "$has_mpi" "$opt_count"
}

solver_perf_tags() {
    local bn="$1"
    local -a tags=()

    [[ "$bn" == *".gpu"* ]] && tags+=("gpu")
    [[ "$bn" == *".avx512"* ]] && tags+=("avx512")
    [[ "$bn" == *".avx2"* ]] && tags+=("avx2")
    [[ "$bn" == *".omp"* ]] && tags+=("omp")
    [[ "$bn" == *".mpi"* ]] && tags+=("mpi")

    if [[ ${#tags[@]} -eq 0 ]]; then
        echo "none"
    else
        (IFS=','; echo "${tags[*]}")
    fi
}

solver_has_interrupt() {
    local bn="$1"
    [[ "$bn" =~ \.(base|bse|mobse|bseEmp|dsm)(\.|$) ]]
}

solver_has_stellar_evolution() {
    local bn="$1"
    [[ "$bn" =~ \.(bse|mobse|bseEmp)(\.|$) ]]
}

solver_has_external_galpy() {
    local bn="$1"
    [[ "$bn" == *".galpy"* ]]
}

solver_has_external_agama() {
    local bn="$1"
    [[ "$bn" == *".agama"* ]]
}

solver_has_external_hard() {
    local bn="$1"
    [[ "$bn" == *".gasdrag"* ]]
}

solver_has_pn() {
    local bn="$1"
    [[ "$bn" == *".pn"* ]]
}

solver_has_64b() {
    local bn="$1"
    [[ "$bn" =~ (^|\.)64b(\.|$) ]]
}

solver_has_mpfrc() {
    local bn="$1"
    [[ "$bn" =~ (^|\.)mp(\.|$) ]]
}

solver_is_debug_g() {
    local bn="$1"
    [[ "$bn" =~ (^|\.)g(\.|$) ]]
}

solver_is_debug_d() {
    local bn="$1"
    [[ "$bn" =~ (^|\.)d(\.|$) ]]
}

emit_filtered_candidates() {
    local title="$1"
    local matcher_fn="$2"
    local max_count="$3"
    local shown=0
    local entry
    local key
    local bn
    local n

    echo "### $title"
    echo

    for entry in "${ranked_solver_entries[@]}"; do
        IFS=$'\t' read -r key bn n <<< "$entry"
        if "$matcher_fn" "$bn"; then
            shown=$((shown + 1))
            echo "- $bn"
            echo "  - perf-tags: $(solver_perf_tags "$bn")"
            echo "  - options: $n"
            echo "  - path: ${solver_path_map[$bn]}"
            [[ $shown -lt $max_count ]] || break
        fi
    done

    if [[ $shown -eq 0 ]]; then
        echo "- (no matching solver detected)"
    fi

    echo
}

token_group_from_configure_feature() {
    local token="$1"
    case "$token" in
        mpi|omp|gpu)
            echo "parallel-runtime"
            ;;
        avx|avx2|avx512|fugaku|64b)
            echo "architecture"
            ;;
        kdk|kdkdk|kdkdk4)
            echo "step-mode"
            ;;
        base|bse|mobse|bseEmp|dsm)
            echo "interrupt-mode"
            ;;
        galpy|agama)
            echo "external-potential"
            ;;
        gasdrag)
            echo "external-hard"
            ;;
        pnall|pnhermite|pnsdar|pnpnall|pnpnhermite)
            echo "post-newtonian"
            ;;
        tt2nd|tt3rd)
            echo "tidal-tensor"
            ;;
        os|pm)
            echo "orbit-method"
            ;;
        mp)
            echo "mpfrc"
            ;;
        g|d)
            echo "debug"
            ;;
        *)
            echo ""
            ;;
    esac
}

extract_opts() {
    local bin="$1"
    local out="$2"
    local raw_help
    local rc=0

    raw_help="$(mktemp)"
    if command -v timeout >/dev/null 2>&1; then
        timeout 8s "$bin" -h > "$raw_help" 2>&1 || rc=$?
    else
        "$bin" -h > "$raw_help" 2>&1 || rc=$?
    fi

        tr '()' '  ' < "$raw_help" \
            | awk '{for(i=1;i<=NF;i++){tok=$i; gsub(/[,:;]+$/, "", tok); if(tok ~ /^-[A-Za-z0-9]$/ || tok ~ /^--[A-Za-z0-9][A-Za-z0-9-]*$/) print tok}}' \
      | sort -u > "$out" || true

    if [[ ! -s "$out" ]]; then
        echo "Warning: no option tokens extracted from $(basename "$bin") (exit=$rc)" >&2
    fi

    rm -f "$raw_help"
}

# 1) Discover commands from PATH and extract options.
mapfile -t discovered_bins < <(discover_petar_bins)

declare -a solver_files=()
declare -a helper_files=()
declare -a script_files=()
declare -a other_files=()
declare -A solver_path_map=()
declare -A helper_path_map=()
declare -A script_path_map=()
declare -A other_path_map=()

for b in "${discovered_bins[@]}"; do
    bn="$(basename "$b")"
    opt_file="$ASSET_DIR/${bn}.options.txt"

    if is_script_tool "$bn"; then
        extract_opts "$b" "$opt_file"
        script_files+=("$opt_file")
        script_path_map["$bn"]="$b"
    elif is_helper_tool "$bn"; then
        # Some hard-debug helpers abort on -h; reuse the corresponding solver help instead.
        if [[ "$bn" == *".hard.debug" ]]; then
            base_bn="${bn%.hard.debug}"
            base_bin="$(command -v "$base_bn" 2>/dev/null || true)"
            if [[ -n "$base_bin" ]]; then
                extract_opts "$base_bin" "$opt_file"
            else
                : > "$opt_file"
                echo "Warning: cannot infer options for $bn (missing $base_bn in PATH)" >&2
            fi
        else
            extract_opts "$b" "$opt_file"
        fi
        helper_files+=("$opt_file")
        helper_path_map["$bn"]="$b"
    else
        extract_opts "$b" "$opt_file"
        if is_solver_by_options "$opt_file"; then
            solver_files+=("$opt_file")
            solver_path_map["$bn"]="$b"
        else
            other_files+=("$opt_file")
            other_path_map["$bn"]="$b"
        fi
    fi
done

IFS=$'\n' solver_files=($(printf '%s\n' "${solver_files[@]}" | sort -u))
IFS=$'\n' helper_files=($(printf '%s\n' "${helper_files[@]}" | sort -u))
IFS=$'\n' script_files=($(printf '%s\n' "${script_files[@]}" | sort -u))
IFS=$'\n' other_files=($(printf '%s\n' "${other_files[@]}" | sort -u))

mapfile -t ranked_solver_entries < <(
    for f in "${solver_files[@]}"; do
        bn="$(basename "$f" .options.txt)"
        n="$(wc -l < "$f")"
        key="$(solver_perf_key "$bn" "$n")"
        printf "%s\t%s\t%s\n" "$key" "$bn" "$n"
    done | sort -r
)

# 2) Summarize configure.ac-related suffix token families seen in solver names.
declare -A token_group_map=()
declare -A token_count_map=()

for f in "${solver_files[@]}"; do
    bn="$(basename "$f" .options.txt)"
    IFS='.' read -r -a parts <<< "$bn"
    for token in "${parts[@]:1}"; do
        group="$(token_group_from_configure_feature "$token")"
        [[ -n "$group" ]] || continue
        token_group_map["$token"]="$group"
        token_count_map["$token"]=$(( ${token_count_map["$token"]:-0} + 1 ))
    done
done

matrix="$ASSET_DIR/option-matrix.md"
{
    echo "# PeTar Option Matrix (Auto-generated)"
    echo
    echo "This file is generated from commands discovered in PATH via <binary> -h option extraction."
    echo
    echo "Discovery rule:"
    echo "- scan PATH for executable names starting with 'petar'"
    echo "- classify helper tools by suffix (*.hard.debug, *.format.transfer)"
    echo "- classify solver binaries by required core options (-u, -t, -o)"
    echo
    echo "## Solver Binaries Included"
    echo
    for f in "${solver_files[@]}"; do
        bn="$(basename "$f" .options.txt)"
        n="$(wc -l < "$f")"
        echo "- $bn ($n options)"
        echo "  - path: ${solver_path_map[$bn]}"
    done

    echo
    echo "## Requirement-Driven Solver Filters"
    echo
    echo "No fixed priority across physics scenarios. Select filters by user requirements first, then choose by performance."
    echo
    echo "Configure mapping:"
    echo "- --with-interrupt -> suffix tokens: base | bse | mobse | bseEmp | dsm"
    echo "- --with-external -> suffix tokens: galpy | agama"
    echo "- --with-external-hard -> suffix token: gasdrag"
    echo "- --with-pn -> suffix tokens: pn*"
    echo
    emit_filtered_candidates "Interrupt Module (--with-interrupt)" solver_has_interrupt 8
    emit_filtered_candidates "Stellar Evolution (bse/mobse/bseEmp)" solver_has_stellar_evolution 8
    emit_filtered_candidates "Long-timescale External: Galpy (--with-external=galpy)" solver_has_external_galpy 8
    emit_filtered_candidates "Long-timescale External: Agama (--with-external=agama)" solver_has_external_agama 8
    emit_filtered_candidates "Short-timescale External Hard: Gas Drag (--with-external-hard=gasdrag)" solver_has_external_hard 8
    emit_filtered_candidates "Post-Newtonian Relativity (--with-pn)" solver_has_pn 8

    echo "### Special-purpose Features (Non-default Recommendations)"
    echo
    echo "These are for special user requests and should not be selected by default."
    echo
    emit_filtered_candidates "High-precision Tree Force (--enable-64b)" solver_has_64b 6
    emit_filtered_candidates "High-precision Position Representation (--enable-mpfrc)" solver_has_mpfrc 6
    emit_filtered_candidates "Debug Build: g (--with-debug=g)" solver_is_debug_g 6
    emit_filtered_candidates "Debug Build: assert (--with-debug=assert)" solver_is_debug_d 6

    echo
    echo "## Recommended Solver Candidates (Performance-first)"
    echo
    echo "Priority order: gpu > avx512 > avx2 > omp > mpi"
    echo "Apply this ordering after choosing a requirement-driven filter."
    if [[ ${#ranked_solver_entries[@]} -gt 0 ]]; then
        rank=0
        for entry in "${ranked_solver_entries[@]}"; do
            IFS=$'\t' read -r key bn n <<< "$entry"
            rank=$((rank + 1))
            [[ $rank -le 12 ]] || break
            echo "- #$rank $bn"
            echo "  - perf-tags: $(solver_perf_tags "$bn")"
            echo "  - options: $n"
            echo "  - path: ${solver_path_map[$bn]}"
        done
    else
        echo "- (no solver binaries detected in PATH)"
    fi

    echo
    echo "## Helper Tool Binaries Excluded From Simulation Commands"
    echo
    for f in "${helper_files[@]}"; do
        bn="$(basename "$f" .options.txt)"
        n="$(wc -l < "$f")"
        echo "- $bn ($n options)"
        echo "  - path: ${helper_path_map[$bn]}"
    done

    echo
    echo "## Script Tools (Non-solver Workflow Commands)"
    echo
    for f in "${script_files[@]}"; do
        bn="$(basename "$f" .options.txt)"
        n="$(wc -l < "$f")"
        echo "- $bn ($n options)"
        echo "  - path: ${script_path_map[$bn]}"
    done

    echo
    echo "## Other Petar Commands In PATH"
    echo
    for f in "${other_files[@]}"; do
        bn="$(basename "$f" .options.txt)"
        n="$(wc -l < "$f")"
        echo "- $bn ($n options)"
        echo "  - path: ${other_path_map[$bn]}"
    done

    echo
    echo "## Detected Configure-Feature Suffix Tokens"
    echo
    if [[ ${#token_group_map[@]} -gt 0 ]]; then
        for token in $(printf "%s\n" "${!token_group_map[@]}" | sort); do
            echo "- .$token"
            echo "  - group: ${token_group_map[$token]}"
            echo "  - binaries: ${token_count_map[$token]}"
        done
    else
        echo "- (none detected)"
    fi

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
    else
        echo "- (no solver binaries detected in PATH)"
    fi
} > "$matrix"

echo "Updated option inventory in: $ASSET_DIR"
