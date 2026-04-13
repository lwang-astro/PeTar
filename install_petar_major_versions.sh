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

echo "[INFO] root: $ROOT_DIR"
echo "[INFO] jobs: $JOBS"
echo "[INFO] bin dir: $BIN_DIR"
echo "[INFO] install prefix: $INSTALL_PREFIX"
echo "[INFO] strict build: $STRICT_BUILD"

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

# Base + major physics families
build_install base
build_install bse --with-interrupt=bse
build_install galpy --with-external=galpy
build_install bse-galpy --with-interrupt=bse --with-external=galpy
build_install agama --with-external=agama
build_install bse-agama --with-interrupt=bse --with-external=agama
if [[ "$HAS_MPFR" == "1" ]]; then
	build_install mpfrc-galpy --with-external=galpy --enable-mpfrc
else
	echo "[WARN] skip mpfrc-galpy: libmpfr not found"
	SKIPPED_TAGS+=("mpfrc-galpy:no-mpfr")
fi

# External-hard families
build_install gasdrag --with-external-hard=gasdrag
build_install bse-gasdrag --with-interrupt=bse --with-external-hard=gasdrag
build_install bse-galpy-gasdrag --with-interrupt=bse --with-external=galpy --with-external-hard=gasdrag
build_install dsm-gasdrag --with-interrupt=dsm --with-external-hard=gasdrag
build_install dsm-galpy-gasdrag --with-interrupt=dsm --with-external=galpy --with-external-hard=gasdrag

# CUDA families
if [[ "$HAS_NVCC" == "1" ]]; then
	build_install cuda --enable-cuda
	build_install cuda-bse --enable-cuda --with-interrupt=bse
	build_install cuda-bse-galpy --enable-cuda --with-interrupt=bse --with-external=galpy
	build_install cuda-bse-galpy-gasdrag --enable-cuda --with-interrupt=bse --with-external=galpy --with-external-hard=gasdrag
	build_install cuda-dsm-galpy-gasdrag --enable-cuda --with-interrupt=dsm --with-external=galpy --with-external-hard=gasdrag
else
	echo "[WARN] skip all cuda variants: nvcc not found"
	SKIPPED_TAGS+=("cuda-family:no-nvcc")
fi

# Precision / PN / validation-specific step modes
build_install 64b --enable-64b
build_install pnall --with-pn=all
build_install bse-pnhermite --with-interrupt=bse --with-pn=hermite

# Explicit T1 comparison binaries: order2=kdk, order4=kdkdk4 (debug g)
build_install kdk-g --with-step-mode=kdk --with-debug=g --with-mpi=yes --disable-omp
build_install kdkdk4-g --with-step-mode=kdkdk4 --with-debug=g --with-mpi=yes --disable-omp

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
