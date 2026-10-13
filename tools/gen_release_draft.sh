#!/bin/bash
# Generate a categorized release.note draft from git history.
#
# Usage: bash tools/gen_release_draft.sh [range]
#   range: git commit range; default master..HEAD. Use the commit of the
#   last release.note update as the base (git log -1 --format=%h -- release.note).
#
# Output (stdout): commit subjects grouped by topic. This is a mechanical
# draft only -- refine into release.note items: merge related commits into
# one thematic item, extend an existing item when a commit continues its
# topic, and skip internal-only changes (lessons, agent infra, doc hygiene).

RANGE="${1:-master..HEAD}"

classify() {
    local s="$1"
    case "$s" in
        hard:*|Hermite*|AR:*) echo "Hard integration / AR" ;;
        bse:*|SEVN*|SSE*) echo "Stellar evolution" ;;
        tools:*|python*|update\ ipynb*) echo "Tools / analysis" ;;
        validation:*|functional:*|test*) echo "Validation / tests" ;;
        build:*|configure*|Makefile*) echo "Build system" ;;
        docs:*|README*|skill*|lesson*|agent*) echo "Documentation / infra" ;;
        *changeover*|*group*|*tidal*|*slowdown*|*orbit*) echo "Hard integration / AR" ;;
        *) echo "Other" ;;
    esac
}

ORDER=("Hard integration / AR" "Stellar evolution" "Tools / analysis" \
       "Validation / tests" "Build system" "Documentation / infra" "Other")
declare -A CAT_LINES
while read -r line; do
    hash="${line%% *}"
    subject="${line#* }"
    [[ "$subject" == Merge* ]] && continue
    CAT_LINES["$(classify "$subject")"]+="$hash $subject"$'\n'
done < <(git log --format='%h %s' "$RANGE")

for cat in "${ORDER[@]}"; do
    [[ -n "${CAT_LINES[$cat]}" ]] || continue
    echo "===== $cat ====="
    printf '%s' "${CAT_LINES[$cat]}"
    echo
done
