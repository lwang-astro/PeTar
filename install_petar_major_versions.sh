#!/usr/bin/env bash
set -euo pipefail

# Rebuild/install controlled PeTar binary families so option inventory is deterministic.
# This script intentionally covers the feature matrix used by validation + skill docs.

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT_DIR"

JOBS="${JOBS:-$(nproc)}"
PURGE_OLD_BINARIES="${PURGE_OLD_BINARIES:-1}"
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

	echo "[INFO] ===== build: $tag ====="
	echo "[INFO] configure args: $*"

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

# Base + major physics families
if want_tag base; then build_install base; fi
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
