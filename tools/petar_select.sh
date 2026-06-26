#!/usr/bin/env bash

set -euo pipefail

usage() {
	cat <<'EOF'
Usage:
  petar.select <suffix|binary-name> [bin-dir]
  petar.select --list [bin-dir]
	petar.select --require feat1,feat2 [--optional feat3,feat4] [bin-dir]
	petar.select --require=feat1,feat2 [--optional=feat3,feat4] [bin-dir]
	petar.select --optional feat1,feat2 [bin-dir]

Description:
  Rebuild symlinks in the install bin directory:
	petar
	petar.hard.debug
	petar.format.transfer
	petar.dump2test

  The target can be either:
	1) a full binary name, e.g. petar.mpi.omp.avx2.agama
	2) a suffix, e.g. .mpi.omp.avx2.agama or mpi.omp.avx2.agama

Examples:
  petar.select .mpi.omp.avx2.agama
  petar.select petar.mpi.omp.avx2.agama
  petar.select --list
	petar.select --optional mpi,omp,avx512,avx2
	petar.select --require bse,agama --optional mpi,omp,avx2

Feature tokens and related configure options:
	merger|base|bse|mobse|bseEmp|dsm -> --with-interrupt=<token>
	galpy|agama           -> --with-external=<token>
	gasdrag               -> --with-external-hard=gasdrag
	mpi                   -> --with-mpi=yes
	omp                   -> OpenMP enabled by default (avoid --disable-omp)
	avx|avx2|avx512       -> --with-simd=<token>
	gpu                   -> --enable-cuda [--with-cuda-prefix=<CUDA_PREFIX>]
	64b                   -> --enable-64b
	mp|mpfrc              -> --enable-mpfrc
	pn* (e.g. pnsdar)     -> --with-pn=<mode>
	g|d                   -> --with-debug=g | --with-debug=assert
EOF
}

resolve_bindir() {
	local script_path
	script_path="$(readlink -f "$0")"
	dirname "$script_path"
}

normalize_target() {
	local input="$1"
	if [[ "$input" == petar* ]]; then
		printf '%s\n' "$input"
	elif [[ "$input" == .* ]]; then
		printf 'petar%s\n' "$input"
	else
		printf 'petar.%s\n' "$input"
	fi
}

split_csv_to_array() {
	local csv="$1"
	local -n out_arr=$2
	out_arr=()
	if [[ -z "$csv" ]]; then
		return 0
	fi
	IFS=',' read -r -a out_arr <<< "$csv"
}

feature_to_configure_hint() {
	local token="$1"
	case "$token" in
		base|merger|bse|mobse|bseEmp|dsm)
			canonical_interrupt_token="$(canonical_feature_token "$token")"
			printf '%s\n' "  - $token -> ./configure --with-interrupt=$canonical_interrupt_token"
			return 0
			;;
		galpy|agama)
			printf '%s\n' "  - $token -> ./configure --with-external=$token"
			return 0
			;;
		gasdrag)
			printf '%s\n' "  - $token -> ./configure --with-external-hard=gasdrag"
			return 0
			;;
		mpi)
			printf '%s\n' "  - mpi -> ./configure --with-mpi=yes"
			return 0
			;;
		omp)
			printf '%s\n' "  - omp -> keep OpenMP enabled (default); do not pass --disable-omp"
			return 0
			;;
		avx|avx2|avx512)
			printf '%s\n' "  - $token -> ./configure --with-simd=$token"
			return 0
			;;
		gpu)
			printf '%s\n' "  - gpu -> ./configure --enable-cuda [--with-cuda-prefix=<CUDA_PREFIX>]"
			return 0
			;;
		64b)
			printf '%s\n' "  - 64b -> ./configure --enable-64b"
			return 0
			;;
		mp|mpfrc)
			printf '%s\n' "  - $token -> ./configure --enable-mpfrc"
			return 0
			;;
		pn*)
			local pn_mode="${token#pn}"
			if [[ -n "$pn_mode" ]]; then
				printf '%s\n' "  - $token -> ./configure --with-pn=$pn_mode"
			else
				printf '%s\n' "  - $token -> ./configure --with-pn=<hermite|sdar|all>"
			fi
			return 0
			;;
		g)
			printf '%s\n' "  - g -> ./configure --with-debug=g"
			return 0
			;;
		d)
			printf '%s\n' "  - d -> ./configure --with-debug=assert"
			return 0
			;;
		*)
			return 1
			;;
	esac
}

is_supported_feature_token() {
	local token="$1"
	case "$token" in
		base|merger|bse|mobse|bseEmp|dsm|galpy|agama|gasdrag|mpi|omp|avx|avx2|avx512|gpu|64b|mp|mpfrc|g|d)
			return 0
			;;
		pn*)
			return 0
			;;
		*)
			return 1
			;;
	esac
}

canonical_feature_token() {
	local token="$1"
	case "$token" in
		base)
			printf '%s\n' "merger"
			;;
		mpfrc)
			printf '%s\n' "mp"
			;;
		*)
			printf '%s\n' "$token"
			;;
	esac
}

is_interrupt_feature_token() {
	local token="$1"
	case "$token" in
		base|merger|bse|mobse|bseEmp|dsm)
			return 0
			;;
		*)
			return 1
			;;
	esac
}

is_require_only_feature_token() {
	local token="$1"
	case "$token" in
		base|merger|bse|mobse|bseEmp|dsm|galpy|agama|gasdrag|64b|gpu|mp|mpfrc)
			return 0
			;;
		pn*)
			return 0
			;;
		*)
			return 1
			;;
	esac
}

array_contains_token() {
	local needle="$1"
	shift
	local item
	for item in "$@"; do
		if [[ "$item" == "$needle" ]]; then
			return 0
		fi
	done
	return 1
}

is_debug_feature_token() {
	local token="$1"
	case "$token" in
		g|d)
			return 0
			;;
		*)
			return 1
			;;
	esac
}

is_performance_feature_token() {
	local token="$1"
	case "$token" in
		gpu|avx|avx2|avx512|omp|mpi)
			return 0
			;;
		*)
			return 1
			;;
	esac
}

print_supported_features() {
	echo "Supported features: merger,base,bse,mobse,bseEmp,dsm,galpy,agama,gasdrag,mpi,omp,avx,avx2,avx512,gpu,64b,mp,mpfrc,pn*,g,d" >&2
}

validate_feature_csv() {
	local csv="$1"
	local argname="$2"
	local tokens=()
	local token
	local cleaned

	split_csv_to_array "$csv" tokens
	for token in "${tokens[@]}"; do
		cleaned="${token//[[:space:]]/}"
		if [[ -z "$cleaned" ]]; then
			continue
		fi
		if ! is_supported_feature_token "$cleaned"; then
			echo "Error: unsupported feature '$cleaned' in $argname" >&2
			print_supported_features
			return 1
		fi
	done
}

filter_supported_feature_csv() {
	local csv="$1"
	local argname="$2"
	local tokens=()
	local token
	local cleaned
	local canonical
	local out_csv=""
	declare -A seen

	split_csv_to_array "$csv" tokens
	for token in "${tokens[@]}"; do
		cleaned="${token//[[:space:]]/}"
		if [[ -z "$cleaned" ]]; then
			continue
		fi
		if ! is_supported_feature_token "$cleaned"; then
			echo "Warning: unsupported feature '$cleaned' in $argname; ignored" >&2
			continue
		fi
		canonical="$(canonical_feature_token "$cleaned")"
		if [[ -n "${seen[$canonical]+x}" ]]; then
			continue
		fi
		seen[$canonical]=1
		if [[ -z "$out_csv" ]]; then
			out_csv="$canonical"
		else
			out_csv="$out_csv,$canonical"
		fi
	done

	printf '%s\n' "$out_csv"
}

print_configure_hints_for_csv() {
	local csv="$1"
	local tokens=()
	local token
	local cleaned
	local line
	local printed=0
	declare -A seen

	split_csv_to_array "$csv" tokens
	for token in "${tokens[@]}"; do
		cleaned="${token//[[:space:]]/}"
		if [[ -z "$cleaned" ]]; then
			continue
		fi
		if [[ -n "${seen[$cleaned]+x}" ]]; then
			continue
		fi
		seen[$cleaned]=1
		if line="$(feature_to_configure_hint "$cleaned")"; then
			if [[ "$printed" -eq 0 ]]; then
				echo "Hint: required features correspond to configure options such as:" >&2
				printed=1
			fi
			echo "$line" >&2
		fi
	done

	if [[ "$printed" -eq 1 ]]; then
		echo "Hint: rebuild and install, then run petar.select again." >&2
	fi
}

has_token() {
	local name="$1"
	local token="$2"
	local part
	IFS='.' read -r -a parts <<< "$name"
	for part in "${parts[@]}"; do
		if [[ "$part" == "$token" ]]; then
			return 0
		fi
	done
	return 1
}

count_optional_hits() {
	local name="$1"
	shift
	local token
	local score=0
	for token in "$@"; do
		if [[ -n "$token" ]] && has_token "$name" "$token"; then
			score=$((score+1))
		fi
	done
	printf '%d\n' "$score"
}

has_interrupt_token_in_name() {
	local name="$1"
	local token
	for token in base bse mobse bseEmp dsm; do
		if has_token "$name" "$token"; then
			return 0
		fi
	done
	return 1
}

is_dsm_requested() {
	local requested_physical_tokens=("$@")
	if array_contains_token "dsm" "${requested_physical_tokens[@]}"; then
		return 0
	fi
	return 1
}

has_unrequested_physical_token_in_name() {
	local name="$1"
	shift
	local requested_physical_tokens=("$@")
	local part
	IFS='.' read -r -a parts <<< "$name"
	for part in "${parts[@]}"; do
		if ! is_require_only_feature_token "$part"; then
			continue
		fi
		if ! array_contains_token "$part" "${requested_physical_tokens[@]}"; then
			return 0
		fi
	done
	return 1
}

has_debug_token_in_name() {
	local name="$1"
	if has_token "$name" g || has_token "$name" d; then
		return 0
	fi
	return 1
}

performance_score() {
	local name="$1"
	local score=0

	if has_token "$name" gpu; then
		score=$((score+10000))
	fi

	if has_token "$name" avx512; then
		score=$((score+1000))
	elif has_token "$name" avx2; then
		score=$((score+800))
	elif has_token "$name" avx; then
		score=$((score+500))
	fi

	if has_token "$name" omp; then
		score=$((score+100))
	fi
	if has_token "$name" mpi; then
		score=$((score+80))
	fi

	printf '%d\n' "$score"
}

count_optional_performance_hits() {
	local name="$1"
	shift
	local token
	local score=0
	for token in "$@"; do
		if [[ -z "$token" ]]; then
			continue
		fi
		if ! is_performance_feature_token "$token"; then
			continue
		fi
		if has_token "$name" "$token"; then
			score=$((score+20))
		fi
	done
	printf '%d\n' "$score"
}

select_by_features() {
	local bindir="$1"
	local require_csv="$2"
	local optional_csv="$3"

	local require_tokens=()
	local optional_tokens=()
	split_csv_to_array "$require_csv" require_tokens
	split_csv_to_array "$optional_csv" optional_tokens

	local candidates=()
	local name
	while IFS= read -r name; do
		candidates+=("$name")
	done < <(list_versions "$bindir")

	local best_name=""
	local best_score=-1
	local token
	local miss
	local score
	local require_physical_tokens=()
	local require_has_debug=0
	local dsm_requested=0

	for token in "${require_tokens[@]}"; do
		if [[ -z "$token" ]]; then
			continue
		fi
		if is_require_only_feature_token "$token"; then
			require_physical_tokens+=("$token")
		fi
		if is_debug_feature_token "$token"; then
			require_has_debug=1
		fi
	done

	if is_dsm_requested "${require_physical_tokens[@]}"; then
		dsm_requested=1
	fi

	for name in "${candidates[@]}"; do
		miss=0
		for token in "${require_tokens[@]}"; do
			if [[ -z "$token" ]]; then
				continue
			fi
			if ! has_token "$name" "$token"; then
				miss=1
				break
			fi
		done
		if [[ "$miss" -eq 1 ]]; then
			continue
		fi

		if has_unrequested_physical_token_in_name "$name" "${require_physical_tokens[@]}"; then
			continue
		fi

		score="$(performance_score "$name")"
		score=$((score + $(count_optional_performance_hits "$name" "${optional_tokens[@]}")))
		if [[ "$require_has_debug" -eq 0 ]] && ! has_debug_token_in_name "$name"; then
			score=$((score+100000))
		fi
		if [[ "$dsm_requested" -eq 0 ]] && has_token "$name" dsm; then
			score=$((score-50000))
		fi
		if [[ "$score" -gt "$best_score" ]]; then
			best_score="$score"
			best_name="$name"
		fi
	done

	if [[ -z "$best_name" ]]; then
		echo "Error: no installed petar version matches required features: ${require_csv}" >&2
		echo "Hint: require-only feature families are included only when explicitly required." >&2
		echo "Hint: this applies to interrupt (base/bse/mobse/bseEmp/dsm), external (galpy/agama), external-hard (gasdrag), pn*, and mpfrc/mp features." >&2
		print_configure_hints_for_csv "$require_csv"
		echo "Available versions in $bindir:" >&2
		list_versions "$bindir" >&2 || true
		return 1
	fi

	printf '%s\n' "$best_name"
}

update_links() {
	local bindir="$1"
	local target="$2"

	local req
	for req in "$target" "$target.hard.debug" "$target.format.transfer" "$target.dump2test"; do
		if [[ ! -f "$bindir/$req" ]]; then
			echo "Error: required file not found: $bindir/$req" >&2
			if [[ "$target" == petar.* ]]; then
				print_configure_hints_for_csv "${target#petar.}"
			fi
			exit 1
		fi
	done

	ln -sfn "$target" "$bindir/petar"
	ln -sfn "$target.hard.debug" "$bindir/petar.hard.debug"
	ln -sfn "$target.format.transfer" "$bindir/petar.format.transfer"
	ln -sfn "$target.dump2test" "$bindir/petar.dump2test"

	echo "Updated links in $bindir"
	echo "  petar -> $target"
	echo "  petar.hard.debug -> $target.hard.debug"
	echo "  petar.format.transfer -> $target.format.transfer"
	echo "  petar.dump2test -> $target.dump2test"
}

list_versions() {
	local bindir="$1"
	local has_any=0

	shopt -s nullglob
	for f in "$bindir"/petar.*; do
		local name
		name="$(basename "$f")"

		if [[ "$name" == *.hard.debug || "$name" == *.format.transfer || "$name" == *.dump2test ]]; then
			continue
		fi

		if [[ -f "$bindir/$name.hard.debug" && -f "$bindir/$name.format.transfer" ]]; then
			has_any=1
			printf '%s\n' "$name"
		fi
	done
	shopt -u nullglob

	if [[ "$has_any" -eq 0 ]]; then
		echo "No selectable petar version found in $bindir" >&2
		return 1
	fi
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
	usage
	exit 0
fi

mode="select"
require_csv=""
optional_csv=""
bindir=""
target_arg=""

while [[ $# -gt 0 ]]; do
	case "$1" in
		--list)
			mode="list"
			shift
			;;
		--require)
			mode="auto"
			if [[ $# -lt 2 ]]; then
				echo "Error: --require needs an argument" >&2
				exit 1
			fi
			require_csv="$2"
			shift 2
			;;
		--require=*)
			mode="auto"
			require_csv="${1#--require=}"
			shift
			;;
		--optional)
			mode="auto"
			if [[ $# -lt 2 ]]; then
				echo "Error: --optional needs an argument" >&2
				exit 1
			fi
			optional_csv="$2"
			shift 2
			;;
		--optional=*)
			mode="auto"
			optional_csv="${1#--optional=}"
			shift
			;;
		-h|--help)
			usage
			exit 0
			;;
		--)
			shift
			break
			;;
		-*)
			echo "Error: unknown option: $1" >&2
			usage >&2
			exit 1
			;;
		*)
			if [[ "$mode" == "select" ]]; then
				if [[ -z "$target_arg" ]]; then
					target_arg="$1"
				elif [[ -z "$bindir" ]]; then
					bindir="$1"
				else
					echo "Error: too many positional arguments" >&2
					exit 1
				fi
			else
				if [[ -z "$bindir" ]]; then
					bindir="$1"
				else
					echo "Error: too many positional arguments" >&2
					exit 1
				fi
			fi
			shift
			;;
	esac
done

if [[ -z "$bindir" ]]; then
	bindir="$(resolve_bindir)"
fi

if [[ ! -d "$bindir" ]]; then
	echo "Error: bin directory does not exist: $bindir" >&2
	if [[ "$bindir" == *,* ]]; then
		echo "Hint: this looks like a feature list. Use --optional <feat1,feat2,...> (or --optional=<...>)" >&2
	fi
	exit 1
fi

if [[ "$mode" == "list" ]]; then
	list_versions "$bindir"
	exit 0
fi

if [[ "$mode" == "auto" ]]; then
	if [[ -z "$require_csv" && -z "$optional_csv" ]]; then
		echo "Error: at least one of --require or --optional must be provided in auto mode" >&2
		exit 1
	fi
	if [[ -n "$require_csv" ]]; then
		validate_feature_csv "$require_csv" "--require"
		require_csv="$(filter_supported_feature_csv "$require_csv" "--require")"
	fi
	optional_csv="$(filter_supported_feature_csv "$optional_csv" "--optional")"
	target="$(select_by_features "$bindir" "$require_csv" "$optional_csv")"
	update_links "$bindir" "$target"
	exit 0
fi

if [[ -z "$target_arg" ]]; then
	usage
	echo
	echo "Available versions in $bindir:" >&2
	list_versions "$bindir" || true
	exit 1
fi

target="$(normalize_target "$target_arg")"
update_links "$bindir" "$target"
