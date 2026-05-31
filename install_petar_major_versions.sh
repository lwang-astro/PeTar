#!/usr/bin/env bash
set -euo pipefail

# Rebuild/install controlled PeTar binary families so option inventory is deterministic.
# This script intentionally covers the feature matrix used by validation + skill docs.

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT_DIR"

JOBS="${JOBS:-$(nproc)}"
PURGE_OLD_BINARIES="${PURGE_OLD_BINARIES:-0}"
BIN_DIR="${BIN_DIR:-$HOME/bin}"
INSTALL_PREFIX="${INSTALL_PREFIX:-$HOME}"
STRICT_BUILD="${STRICT_BUILD:-0}"
BUILD_TAGS="${BUILD_TAGS:-all}"

echo "[INFO] root: $ROOT_DIR"
echo "[INFO] jobs: $JOBS"
echo "[INFO] bin dir: $BIN_DIR"
echo "[INFO] install prefix: $INSTALL_PREFIX"
echo "[INFO] strict build: $STRICT_BUILD"
echo "[INFO] build tags: $BUILD_TAGS"

HAS_MPFR=0
if ldconfig -p 2>/dev/null | grep -q "libmpfr.so"; then
	HAS_MPFR=1
fi

HAS_NVCC=0
if command -v nvcc >/dev/null 2>&1; then
	HAS_NVCC=1
fi

echo "[INFO] has mpfr: $HAS_MPFR"
echo "[INFO] has nvcc: $HAS_NVCC"

FAILED_TAGS=()
SKIPPED_TAGS=()
REQUESTED_TAGS=()

current_source_version() {
	local petar_version=""
	local sdar_version=""

	if [[ -f "$ROOT_DIR/VERSION" ]]; then
		petar_version="$(tr -d '[:space:]' < "$ROOT_DIR/VERSION")"
	fi
	if [[ -f "$ROOT_DIR/../SDAR/VERSION" ]]; then
		sdar_version="$(tr -d '[:space:]' < "$ROOT_DIR/../SDAR/VERSION")"
	fi

	if [[ -z "$petar_version" || -z "$sdar_version" ]]; then
		return 1
	fi

	printf '%s_%s\n' "$petar_version" "$sdar_version"
}

solver_matches_tag() {
	local solver_name="$1"
	local tag="$2"
	local token
	local core
	local parts_csv
	local parts_str
	local core_tokens=" merger dsm bse galpy agama "

	if [[ -z "$tag" || "$tag" == "std" ]]; then
		# For std, require a solver without any core physics tags.
		parts_csv="${solver_name#petar}"
		parts_csv="${parts_csv#.}"
		parts_str=" ${parts_csv//./ } "
		for core in merger dsm bse galpy agama; do
			if [[ "$parts_str" == *" $core "* ]]; then
				return 1
			fi
		done
		return 0
	fi

	# Parse solver features from dot-separated executable suffixes.
	# Example: petar.mpi.omp.avx2.bse.agama -> " mpi omp avx2 bse agama "
	parts_csv="${solver_name#petar}"
	parts_csv="${parts_csv#.}"
	parts_str=" ${parts_csv//./ } "

	IFS='-' read -r -a tokens <<< "$tag"
	for token in "${tokens[@]}"; do
		[[ -z "$token" ]] && continue
		if [[ "$parts_str" != *" $token "* ]]; then
			return 1
		fi
	done

	# Prevent mixed-physics false matches, e.g. treat bse-agama as NOT matching bse.
	for core in merger dsm bse galpy agama; do
		if [[ "$parts_str" == *" $core "* && "$tag" != *"$core"* ]]; then
			return 1
		fi
	done

	return 0
}

extract_solver_version() {
	local solver_path="$1"
	"$solver_path" -h 2>&1 | sed -n 's/^Version:[[:space:]]*//p' | head -n 1
}

has_matching_installed_version() {
	local tag="$1"
	local expected_version="$2"
	local solver_path=""
	local solver_name=""
	local installed_version=""

	shopt -s nullglob
	for solver_path in "$BIN_DIR"/petar.*; do
		solver_name="$(basename "$solver_path")"

		if [[ "$solver_name" == *.hard.debug || "$solver_name" == *.format.transfer ]]; then
			continue
		fi
		if [[ ! -x "$solver_path" ]]; then
			continue
		fi
		if [[ ! -f "$BIN_DIR/$solver_name.hard.debug" || ! -f "$BIN_DIR/$solver_name.format.transfer" ]]; then
			continue
		fi
		if ! solver_matches_tag "$solver_name" "$tag"; then
			continue
		fi

		installed_version="$(extract_solver_version "$solver_path")"
		if [[ "$installed_version" == "$expected_version" ]]; then
			shopt -u nullglob
			return 0
		fi
		done
	shopt -u nullglob
	return 1
}

if [[ "$PURGE_OLD_BINARIES" == "1" ]]; then
	if [[ -d "$BIN_DIR" ]]; then
		echo "[INFO] removing previous petar* executables in $BIN_DIR"
		find "$BIN_DIR" -maxdepth 1 -type f -name 'petar*' -print -delete || true
		find "$BIN_DIR" -maxdepth 1 -type l -name 'petar*' -print -delete || true
	fi
fi

build_install() {
	local tag="$1"
	shift
	local source_version=""

	echo "[INFO] ===== build: $tag ====="
	echo "[INFO] configure args: $*"

	if source_version="$(current_source_version)"; then
		echo "[INFO] source version: $source_version"
		if has_matching_installed_version "$tag" "$source_version"; then
			echo "[INFO] skip rebuild for $tag: installed version matches source version $source_version"
			return 0
		fi
		echo "[INFO] rebuild $tag: no installed solver with matching version $source_version"
	else
		echo "[WARN] source version detection failed; continue with rebuild for $tag"
	fi

	if ! ./configure --prefix="$INSTALL_PREFIX" "$@"; then
		echo "[WARN] configure failed for $tag"
		FAILED_TAGS+=("$tag:configure")
		if [[ "$STRICT_BUILD" == "1" ]]; then
			return 1
		fi
		return 0
	fi
	make clean
	if ! make -j"$JOBS"; then
		echo "[WARN] make failed for $tag"
		FAILED_TAGS+=("$tag:make")
		if [[ "$STRICT_BUILD" == "1" ]]; then
			return 1
		fi
		return 0
	fi
	if ! make install; then
		echo "[WARN] make install failed for $tag"
		FAILED_TAGS+=("$tag:install")
		if [[ "$STRICT_BUILD" == "1" ]]; then
			return 1
		fi
		return 0
	fi
}

if [[ "$BUILD_TAGS" != "all" ]]; then
	IFS=',' read -r -a REQUESTED_TAGS <<< "$BUILD_TAGS"
fi

normalize_tag() {
	local tag="$1"
	echo "$tag" | tr '[:upper:]' '[:lower:]' | xargs
}

want_tag() {
	local tag
	tag="$(normalize_tag "$1")"
	if [[ "$BUILD_TAGS" == "all" ]]; then
		return 0
	fi
	for req in "${REQUESTED_TAGS[@]}"; do
		if [[ "$(normalize_tag "$req")" == "$tag" ]]; then
			return 0
		fi
	done
	return 1
}

# Merger + major physics families
if want_tag std; then build_install std; fi
if want_tag merger || want_tag base; then build_install merger --with-interrupt=merger; fi
if want_tag dsm; then build_install dsm --with-interrupt=dsm; fi
if want_tag bse; then build_install bse --with-interrupt=bse; fi
if want_tag galpy; then build_install galpy --with-external=galpy; fi
if want_tag bse-galpy; then build_install bse-galpy --with-interrupt=bse --with-external=galpy; fi
if want_tag agama; then build_install agama --with-external=agama; fi
if want_tag bse-agama; then build_install bse-agama --with-interrupt=bse --with-external=agama; fi
if want_tag mpfrc-galpy; then
	if [[ "$HAS_MPFR" == "1" ]]; then
		build_install mpfrc-galpy --with-external=galpy --enable-mpfrc
	else
		echo "[WARN] skip mpfrc-galpy: libmpfr not found"
		SKIPPED_TAGS+=("mpfrc-galpy:no-mpfr")
	fi
fi

# External-hard families
if want_tag gasdrag; then build_install gasdrag --with-external-hard=gasdrag; fi
if want_tag bse-gasdrag; then build_install bse-gasdrag --with-interrupt=bse --with-external-hard=gasdrag; fi
if want_tag bse-galpy-gasdrag; then build_install bse-galpy-gasdrag --with-interrupt=bse --with-external=galpy --with-external-hard=gasdrag; fi
if want_tag dsm-gasdrag; then build_install dsm-gasdrag --with-interrupt=dsm --with-external-hard=gasdrag; fi
if want_tag dsm-galpy-gasdrag; then build_install dsm-galpy-gasdrag --with-interrupt=dsm --with-external=galpy --with-external-hard=gasdrag; fi

# CUDA families
if want_tag cuda || want_tag cuda-bse || want_tag cuda-bse-galpy || want_tag cuda-bse-galpy-gasdrag || want_tag cuda-dsm-galpy-gasdrag; then
	if [[ "$HAS_NVCC" == "1" ]]; then
		if want_tag cuda; then build_install cuda --enable-cuda; fi
		if want_tag cuda-bse; then build_install cuda-bse --enable-cuda --with-interrupt=bse; fi
		if want_tag cuda-bse-galpy; then build_install cuda-bse-galpy --enable-cuda --with-interrupt=bse --with-external=galpy; fi
		if want_tag cuda-bse-galpy-gasdrag; then build_install cuda-bse-galpy-gasdrag --enable-cuda --with-interrupt=bse --with-external=galpy --with-external-hard=gasdrag; fi
		if want_tag cuda-dsm-galpy-gasdrag; then build_install cuda-dsm-galpy-gasdrag --enable-cuda --with-interrupt=dsm --with-external=galpy --with-external-hard=gasdrag; fi
	else
		echo "[WARN] skip requested cuda variants: nvcc not found"
		SKIPPED_TAGS+=("cuda-family:no-nvcc")
	fi
fi

# Precision / PN / validation-specific step modes
if want_tag 64b; then build_install 64b --enable-64b; fi
if want_tag pnall; then build_install pnall --with-pn=all; fi
if want_tag bse-pnhermite; then build_install bse-pnhermite --with-interrupt=bse --with-pn=hermite; fi

# Explicit T1 comparison binaries: order2=kdk, order4=kdkdk4 (debug g)
if want_tag kdk-g; then build_install kdk-g --with-step-mode=kdk --with-debug=g --with-mpi=yes --disable-omp; fi
if want_tag kdkdk4-g; then build_install kdkdk4-g --with-step-mode=kdkdk4 --with-debug=g --with-mpi=yes --disable-omp; fi

if [[ "${#SKIPPED_TAGS[@]}" -gt 0 ]]; then
	echo "[INFO] skipped variants: ${SKIPPED_TAGS[*]}"
fi

if [[ "${#FAILED_TAGS[@]}" -gt 0 ]]; then
	echo "[WARN] failed variants: ${FAILED_TAGS[*]}"
	if [[ "$STRICT_BUILD" == "1" ]]; then
		exit 1
	fi
fi

echo "[INFO] build/install pass completed"
