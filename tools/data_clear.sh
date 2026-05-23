#!/bin/bash

set -o pipefail

# fallback time-column indices (1-based, include event-label token)
tindex_sse_type_change=25
tindex_sse_sn_kick=14
tindex_bse_type_change=27
tindex_bse_dyn_merge=44
tindex_bse_sn_kick=16

ncol_sse_type_change=26
ncol_sse_sn_kick=15
ncol_bse_type_change=54
ncol_bse_dyn_merge=57
ncol_bse_sn_kick=17
less_dyn_merge=0
less_type_change=0

unset tcrit
interrupt_mode=none
external_mode=none
status_n_particle=0
unset use_mpfrc
petar_binary=petar
auto_mode=1
ignore_binary_warnings=0
interrupt_mode_set=0
external_mode_set=0
clear_tmp=0

detect_modes_from_petar() {
	local pbin=$1
	local ppath pname htxt

	ppath=$(command -v "$pbin" 2>/dev/null)
	if [ -z "$ppath" ]; then
		echo "Warning! cannot find petar binary '$pbin' in PATH, skip auto-detection." >&2
		return 1
	fi

	if [ -L "$ppath" ]; then
		ppath=$(readlink -f "$ppath")
	fi
	pname=$(basename "$ppath")

	# Priority 1: infer from executable suffix tokens (most reliable for compiled features)
	if [ $interrupt_mode_set -eq 0 ]; then
		if [[ $pname == *'.bseEmp'* ]]; then
			interrupt_mode=bseEmp
		elif [[ $pname == *'.mobse'* ]]; then
			interrupt_mode=mobse
		elif [[ $pname == *'.bse'* ]]; then
			interrupt_mode=bse
		elif [[ $pname == *'.dsm'* ]]; then
			interrupt_mode=dsm
		elif [[ $pname == *'.base'* ]]; then
			interrupt_mode=base
		fi
	fi

	if [ $external_mode_set -eq 0 ]; then
		if [[ $pname == *'.galpy'* ]]; then
			external_mode=galpy
		elif [[ $pname == *'.agama'* ]]; then
			external_mode=agama
		fi
	fi

	# Priority 2: fallback to help-option probing when suffixes are not informative.
	if [ $interrupt_mode_set -eq 0 ] || [ $external_mode_set -eq 0 ]; then
		htxt=$($pbin -h 2>&1)
		if [ $interrupt_mode_set -eq 0 ] && [ "$interrupt_mode" = "none" ]; then
			if echo "$htxt" | grep -Eq -- '--dsm-'; then
				interrupt_mode=dsm
			elif echo "$htxt" | grep -Eq -- '--bse-|--stellar-evolution'; then
				interrupt_mode=bse
			fi
		fi
		if [ $external_mode_set -eq 0 ] && [ "$external_mode" = "none" ]; then
			if echo "$htxt" | grep -Eq -- '--galpy-'; then
				external_mode=galpy
			elif echo "$htxt" | grep -Eq -- '--agama-'; then
				external_mode=agama
			fi
		fi
	fi

	echo "Auto-detected parser modes from $pname: interrupt=$interrupt_mode, external=$external_mode" >&2
	return 0
}

has_binary_status_or_group() {
	local prefix=$1
	local nmpi_local=$2
	local f

	if [ -e "$prefix.status" ] && is_binary_file "$prefix.status"; then
		return 0
	fi

	if [ ! -z "$nmpi_local" ]; then
		local nend ir
		nend=`expr $nmpi_local - 1`
		for ir in `seq 0 $nend`
		do
			for f in `ls 2>/dev/null | egrep '^'$prefix'.group.'$ir'(\.n[0-9]+)?$'`
			do
				if is_binary_file "$f"; then
					return 0
				fi
			done
		done
	else
		for f in `ls 2>/dev/null | egrep '^'$prefix'.group.[0-9]+(\.n[0-9]+)?$'`
		do
			if is_binary_file "$f"; then
				return 0
			fi
		done
	fi

	return 1
}

is_binary_file() {
	local f=$1
	if [ ! -s "$f" ]; then
		return 1
	fi
	if python3 - "$f" <<'PY'
import sys

path = sys.argv[1]
with open(path, 'rb') as handle:
    while True:
        chunk = handle.read(4096)
        if not chunk:
            sys.exit(1)
        if b'\x00' in chunk:
            sys.exit(0)
PY
	then
		return 0
	fi
	if ! LC_ALL=C head -c 1024 "$f" | LC_ALL=C grep -q '^[[:print:][:space:]]*$'; then
		return 0
	fi
	return 1
}

filter_binary_file_by_time() {
	local kind=$1
	local fin=$2
	local fout=$3
	local tcut=$4
	local n_group_member=${5:-0}
	local tmpout="${fout}.tmp.$$"

	if ! python3 - "$kind" "$fin" "$tmpout" "$tcut" "$n_group_member" "$interrupt_mode" "$external_mode" "$status_n_particle" "${use_mpfrc:-0}" "$ignore_binary_warnings" <<'PY'
import sys
import warnings

kind = sys.argv[1]
fin = sys.argv[2]
fout = sys.argv[3]
tcrit = float(sys.argv[4])
n_group = int(sys.argv[5])
interrupt_mode = sys.argv[6]
external_mode = sys.argv[7]
status_n_particle = int(sys.argv[8])
use_mpfrc = bool(int(sys.argv[9]))
ignore_binary_warnings = bool(int(sys.argv[10]))

import petar

common = {
	"interrupt_mode": interrupt_mode,
	"external_mode": external_mode,
	"use_mpfrc": use_mpfrc,
}

if kind == "status":
	if status_n_particle > 0:
		common["N_particle"] = status_n_particle
	data = petar.Status(**common)
elif kind == "group":
	common["N"] = n_group
	data = petar.GroupInfo(**common)
else:
	raise ValueError(f"Unsupported kind: {kind}")

try:
	with warnings.catch_warnings(record=True) as caught:
		warnings.simplefilter("always")
		data.fromfile(fin)
	if caught:
		messages = [str(item.message) for item in caught]
		if not ignore_binary_warnings:
			raise RuntimeError(f"Abort clearing {kind}: warning raised while reading {fin}: {' | '.join(messages)}")
		for message in messages:
			print(f"Warning! ignore binary warning for {kind} file {fin}: {message}", file=sys.stderr)
	selected = data[data.time <= tcrit]
	selected.tofile(fout)
except Warning as warning:
	if not ignore_binary_warnings:
		raise RuntimeError(f"Abort clearing {kind}: warning raised while reading {fin}: {warning}")
	print(f"Warning! ignore binary warning for {kind} file {fin}: {warning}", file=sys.stderr)
except RuntimeError:
	raise
PY
	then
		rm -f "$tmpout"
		return 1
	fi
	mv "$tmpout" "$fout"
}

initialize_event_column_counts() {
	local pyout
	if ! pyout=$(python3 - "$less_type_change" "$less_dyn_merge" <<'PY'
import sys

less_type_change = bool(int(sys.argv[1]))
less_dyn_merge = bool(int(sys.argv[2]))

import petar

def ncols_of_type(tp):
	if isinstance(tp, tuple):
		return int(tp[1])
	if hasattr(tp, 'ncols'):
		return int(tp.ncols)
	try:
		obj = tp()
		if hasattr(obj, 'ncols'):
			return int(obj.ncols)
	except Exception:
		pass
	return 1

def time_index_with_label(obj, path):
	parts = path.split('.')
	cur = obj
	idx = 1
	for i, part in enumerate(parts):
		found = False
		for key, tp in cur.keys:
			if key == part:
				if i == len(parts) - 1:
					return idx + 1
				cur = getattr(cur, key)
				found = True
				break
			idx += ncols_of_type(tp)
		if not found and i != len(parts) - 1:
			raise KeyError(path)
	raise KeyError(path)

sse_type_change = petar.SSETypeChange()
sse_sn_kick = petar.SSESNKick()
bse_type_change = petar.BSETypeChange(base_output=less_type_change)
bse_dyn_merge = petar.BSEDynamicMerge(less_output=less_dyn_merge)
bse_sn_kick = petar.BSEKick()

values = {
	'ncol_sse_type_change': sse_type_change.ncols + 1,
	'ncol_sse_sn_kick': sse_sn_kick.ncols + 1,
	'ncol_bse_type_change': bse_type_change.ncols + 1,
	'ncol_bse_dyn_merge': bse_dyn_merge.ncols + 1,
	'ncol_bse_sn_kick': bse_sn_kick.ncols + 1,
	'tindex_sse_type_change': time_index_with_label(sse_type_change, 'final.time'),
	'tindex_sse_sn_kick': time_index_with_label(sse_sn_kick, 'star.time'),
	'tindex_bse_type_change': time_index_with_label(bse_type_change, 'final.time'),
	'tindex_bse_dyn_merge': time_index_with_label(bse_dyn_merge, 'final.p1.time'),
	'tindex_bse_sn_kick': time_index_with_label(bse_sn_kick, 'star.time'),
}

for key, value in values.items():
    print(f"{key}={value}")
PY
	); then
		echo 'Warning! failed to derive SSE/BSE column counts from petar analysis classes; fallback to built-in defaults.' >&2
		return 1
	fi
	eval "$pyout"
	echo 'ASCII event layout from petar analysis: sse(type_change ncol='${ncol_sse_type_change}', t='${tindex_sse_type_change}'; sn_kick ncol='${ncol_sse_sn_kick}', t='${tindex_sse_sn_kick}') bse(type_change ncol='${ncol_bse_type_change}', t='${tindex_bse_type_change}'; dyn_merge ncol='${ncol_bse_dyn_merge}', t='${tindex_bse_dyn_merge}'; sn_kick ncol='${ncol_bse_sn_kick}', t='${tindex_bse_sn_kick}')' >&2
	return 0
}

until [[ `echo x$1` == 'x' ]]
do
    case $1 in
	-h) shift;
	    echo 'A tool for clearing data after a specified time.';
	    echo '   Remove data after a specified time criterion for output files with name suffixes: '${suffixes[@]};
	    echo '   When users wish to restart a simulation from an outputted snapshot file, other output files may contain events that occurred after the time of this snapshot file.';
	    echo '   This tool can help clear up these events before restarting.';
	    echo '   Subsequently, when restarting the simulation, the same events will not be recorded twice.';
	    echo 'Usage: petar.data.clear [options] [data filename prefix]';
	    echo '       The data filename prefix is defined by "petar -f"; the default is "data".';
	    echo 'Options (default arguments shown in parentheses at the end):';
	    echo '  -t [F]: time criterion for clearing up data, must be provided (default: none)';
	    echo '  -n [I]: number of MPI processes (default: auto)';
	    echo '  -b    : use previous backup files instead of replacing (default: replacing)';
	    echo '  --interrupt-mode [S]: mode for binary status/group parsing: none, base, bse, bseEmp, mobse, dsm (default: auto-detect from active petar binary)';
	    echo '  --external-mode [S]: mode for binary status/group parsing: none, galpy, agama (default: auto-detect from active petar binary)';
	    echo '  --petar-binary [S]: petar command name/path used for auto-detection (default: petar)';
	    echo '  --no-auto-mode: disable auto-detection and keep user-provided/default mode values';
	    echo '  --status-n-particle [I]: particle count for binary status when petar uses -w 2 (default: 0, i.e. no embedded particles)';
	    echo '  --use-mpfrc: enable mpfrc layout for binary status/group parsing (default: off)';
	    echo '  --ignore-binary-warnings: continue clearing binary status/group even if Python readers emit layout warnings (default: abort on warnings)';
	    echo '  --clear-tmp: remove residual transactional tmp files after data clearing (prefix.*.tmp and object_*.tmp)';
	    echo '  --less-dyn-merge:   reduced output mode for [prefix].*bse*.dynmical_merge (three columns less)';
	    echo '                      Only applicable for PeTar versions before Sep 10, 2020)';
	    echo '  --less-type-change: reduced output mode for [prefix].*bse*.type_change (20 columns less)';
	    echo '                      Only applicable for PeTar versions before Jun 2, 2022)';
	    exit;;
	-n) shift; nmpi=$1; shift;;
	-t) shift; tcrit=$1; shift;;
	-b) rmi=1; shift;;
	--interrupt-mode) shift; interrupt_mode=$1; interrupt_mode_set=1; shift;;
	--external-mode) shift; external_mode=$1; external_mode_set=1; shift;;
	--petar-binary) shift; petar_binary=$1; shift;;
	--no-auto-mode) auto_mode=0; shift;;
	--status-n-particle) shift; status_n_particle=$1; shift;;
	--use-mpfrc) use_mpfrc=1; shift;;
	--ignore-binary-warnings) ignore_binary_warnings=1; shift;;
	--clear-tmp) clear_tmp=1; shift;;
	--less-dyn-merge) less_dyn_merge=1; shift;; 
	--less-type-change) less_type_change=1; shift;;
	*) fname=$1;shift;;
    esac
done

suffixes=(esc group sse mosse sseEmp bse mobse bseEmp interrupt status prof.rank)
tindices=(1 3 0 0 0 0 0 0 1 1 2)
nsuffixes=${#suffixes[@]}

if [ ! -e $fname ] | [ -z $fname ] ; then
    echo 'Error, file name not provided' 
	exit 1
fi

if [ -z $tcrit ]; then
    echo 'Time criterion not provided'
	exit 1
fi

if [ $auto_mode -eq 1 ]; then
	if has_binary_status_or_group "$fname" "$nmpi"; then
		detect_modes_from_petar "$petar_binary"
	else
		echo "Skip auto-detection: no binary status/group files detected for prefix '$fname'." >&2
	fi
fi

echo 'data filename prefix: '$fname
echo 'time criterion: '$tcrit
echo 'binary parse mode: interrupt='${interrupt_mode}', external='${external_mode}', status_n_particle='${status_n_particle}', use_mpfrc='${use_mpfrc:-0}', petar_binary='${petar_binary}', auto_mode='${auto_mode}', ignore_binary_warnings='${ignore_binary_warnings}

initialize_event_column_counts

# check consistence for the number of columns
# SSE/BSE event files in current PeTar are ASCII outputs.
ncol=`egrep -m 1 'SN_kick' $fname.*sse*.0|wc -w`
if [[ $ncol != $ncol_sse_sn_kick ]] && [[ $ncol != 0 ]]; then
    echo 'Error! column number not matches for SSE SN kick, should be '$ncol_sse_sn_kick', the file has '$ncol'.'
	exit 1
fi
ncol=`egrep -v -m 1 'SN_kick' $fname.*sse*.0|wc -w`
if [[ $ncol -ne $ncol_sse_type_change ]] && [[ $ncol -ne 0 ]]; then
    echo 'Error! column number not matches for SSE Type Change, should be '$ncol_sse_type_change', the file has '$ncol
	exit 1
fi

ncol=`egrep -m 1 'Dynamic_merge' $fname.*bse*.0|wc -w`
if [[ $ncol -ne $ncol_bse_dyn_merge ]] && [[ $ncol -ne 0 ]]; then
    echo 'Error! column number not matches for BSE dynamical merger, should be '$ncol_bse_dyn_merge', the file has '$ncol
	exit 1
fi

ncol=`egrep -m 1 'SN_kick' $fname.*bse*.0|wc -w`
if [[ $ncol -ne $ncol_bse_sn_kick ]] && [[ $ncol -ne 0 ]]; then
    echo 'Error! column number not matches for BSE SN kick, should be '$ncol_bse_sn_kick', the file has '$ncol
	exit 1
fi

ncol=`egrep -v -m 1 '(SN_kick|Dynamic_merge)' $fname.*bse*.0|wc -w`
if [[ $ncol -ne $ncol_bse_type_change ]] && [[ $ncol -ne 0 ]]; then
    echo 'Error! column number not matches for BSE Type Change, should be '$ncol_bse_type_change', the file has '$ncol
	exit 1
fi


# clear each file
for ((i=0;i<$nsuffixes;i=i+1))
do
    s=${suffixes[$i]}
    tindex=${tindices[$i]}
    file=$fname.$s
    echo $file'; time column: '$tindex
    if [ $s == 'status' ]; then
	# status output is a single file (no rank splitting)
	if [ -e "$file" ]; then
	    lst="$file"
	else
	    lst=''
	fi
    elif [ $s == 'interrupt' ]; then
	# interrupt output may appear as a single file or rank-suffixed files
	lst=''
	if [ -e "$file" ]; then
	    lst="$file"
	fi
	if [ ! -z $nmpi ]; then
	    nend=`expr $nmpi - 1`
	    for ir in `seq 0 $nend`
	    do
		if [ -e "$file.$ir" ]; then
		    lst="$lst $file.$ir"
		fi
	    done
	else
	    for itf in `ls 2>/dev/null | egrep '^'$file'.[0-9]+$'`
	    do
		lst="$lst $itf"
	    done
	fi
    elif [ $s == 'group' ]; then
	if [ ! -z $nmpi ]; then
	    nend=`expr $nmpi - 1`
	    lst=''
	    for ir in `seq 0 $nend`
	    do
		if [ -e $file.$ir ]; then
		    lst="$lst $file.$ir"
		fi
		for gf in `ls 2>/dev/null | egrep '^'$file'.'$ir'.n[0-9]+$'`
		do
		    lst="$lst $gf"
		done
	    done
	else
	    lst=`ls | egrep '^'$file'.[0-9]+(\.n[0-9]+)?$'`
	fi
    elif [ ! -z $nmpi ]; then
	nend=`expr $nmpi - 1`
	lst=`seq 0 $nend|awk -v f=$file '{printf("$s.$d\n", f,$1)}`
    else
	lst=`ls |egrep $file'.[0-9]+$'`
    fi
    for f in $lst
    do
	echo 'process '$f
	if [ -e $f.bk ]; then
	    if [ -z $rmi ]; then
		echo 'Start to remove previous backup file, if you do not want this, please quit petar.data.clear (ctrl c) immediately!'
		rm -i $f.bk
		mv $f $f.bk
		echo 'backup '$f' to '$f.bk
	    fi
	else
	    mv $f $f.bk
	    echo 'backup '$f' to '$f.bk
	fi
	if [[ $s == *'sse'* ]]; then
	    if is_binary_file $f.bk; then
		echo 'Warning! unexpected binary SSE event file '$f' detected; current PeTar SSE event outputs are treated as ASCII. Keep original backup content.'
		cp $f.bk $f
		continue
	    fi
	    awk -v t=$tcrit -v tsn=$tindex_sse_sn_kick -v ttch=$tindex_sse_type_change '{if (($1!="SN_kick" && $ttch<=t) || ($1=="SN_kick" && $tsn<=t)) print $LINE}' $f.bk >$f
	elif [[ $s == *'bse'* ]]; then
	    if is_binary_file $f.bk; then
		echo 'Warning! unexpected binary BSE event file '$f' detected; current PeTar BSE event outputs are treated as ASCII. Keep original backup content.'
		cp $f.bk $f
		continue
	    fi
	    awk -v t=$tcrit -v tsn=$tindex_bse_sn_kick -v ttch=$tindex_bse_type_change -v tdyn=$tindex_bse_dyn_merge '{if ($1=="SN_kick") {if ($tsn<=t) print $LINE} else if ($1=="Dynamic_merge:") {if($tdyn<=t) print $LINE} else if ($ttch<=t) print $LINE;}' $f.bk >$f
	elif [ $s == 'interrupt' ]; then
	    if is_binary_file $f.bk; then
		echo 'Warning! unexpected binary interrupt file '$f' detected; interrupt outputs are treated as ASCII. Keep original backup content.'
		cp $f.bk $f
		continue
	    fi
	    awk -v t=$tcrit -v ti=$tindex '{if ($ti<=t) print $LINE}' $f.bk >$f
	elif is_binary_file $f.bk; then
	    if [ $s == 'status' ]; then
		echo 'binary '$s': use Python parser to clear by time criterion'
		if ! filter_binary_file_by_time status $f.bk $f $tcrit; then
		    echo 'Error! failed to clear binary status file '$f
		    exit 1
		fi
	    elif [ $s == 'group' ]; then
		n_group=0
		if [[ $f =~ \.n([0-9]+)$ ]]; then
		    n_group=${BASH_REMATCH[1]}
		fi
		if [ $n_group -le 0 ]; then
		    echo 'Warning! cannot determine N from group filename '$f', skip this file.'
		    cp $f.bk $f
		    continue
		fi
		echo 'binary '$s': use Python parser to clear by time criterion (N='$n_group')'
		if ! filter_binary_file_by_time group $f.bk $f $tcrit $n_group; then
		    echo 'Error! failed to clear binary group file '$f
		    exit 1
		fi
	    else
		echo 'Warning! binary file '$f' is not supported in this script branch, keep original backup content.'
		cp $f.bk $f
	    fi
	else
	    awk -v t=$tcrit -v ti=$tindex '{if ($ti<=t) print $LINE}' $f.bk >$f
	fi
    done
done

#

if [ $clear_tmp -eq 1 ]; then
	echo 'remove transactional tmp files'
	tmp_list=''
	for tf in `ls 2>/dev/null | egrep '^'$fname'.*\.tmp$'`
	do
		tmp_list="$tmp_list $tf"
	done
	for tf in `ls 2>/dev/null | egrep '^object_.*\.tmp$'`
	do
		tmp_list="$tmp_list $tf"
	done

	n_tmp=0
	for tf in $tmp_list
	do
		if [ -e "$tf" ]; then
			rm -f "$tf"
			n_tmp=`expr $n_tmp + 1`
		fi
	done
	echo 'removed tmp files: '$n_tmp
fi



