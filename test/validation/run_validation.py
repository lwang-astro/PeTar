#!/usr/bin/env python3
import argparse
import json
import math
import os
import re
import shlex
import struct
import subprocess
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

from metrics import (
    check_convergence_ratio,
    check_loglog_slope,
    check_max_threshold,
    check_metric_ratio,
    check_min_threshold,
    summarize_results,
)


ENERGY_HEADER_RE = re.compile(r"^Energy:\s+(.*)$")
ENERGY_PHYSIC_RE = re.compile(r"^Physic:\s+(.*)$")
TIME_RE = re.compile(r"^Time:\s*([0-9eE+\-.]+)")
COUNT_RE = re.compile(
    r"^Time:\s*([0-9eE+\-.]+)\s+N_real\(loc\):\s*(\d+)\s+N_real\(glb\):\s*(\d+)\s+N_all\(loc\):\s*(\d+)\s+N_all\(glb\):\s*(\d+)"
)
G_MSUN_PC_MYR = 0.00449830997959438


@dataclass
class RunResult:
    run_id: str
    command: str
    output_path: Path
    metrics: Dict[str, Any]


def load_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def load_criteria(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    return load_json(path)


def resolve_threshold(
    check: Dict[str, Any],
    scenario_name: str,
    criteria: Dict[str, Any],
    metric: str,
    key: str,
    fallback: float,
) -> float:
    section = criteria.get(scenario_name, {})
    metric_map = section.get(metric, {})
    if isinstance(metric_map, dict) and key in metric_map:
        return float(metric_map[key])

    ref = check.get("criteria_ref")
    if isinstance(ref, list) and len(ref) == 3:
        sname, mname, kname = ref
        if sname in criteria and isinstance(criteria[sname], dict):
            cmap = criteria[sname].get(mname, {})
            if isinstance(cmap, dict) and kname in cmap:
                return float(cmap[kname])

    return fallback


def resolve_convergence(check: Dict[str, Any], scenario_name: str, criteria: Dict[str, Any]) -> Tuple[float, float]:
    expected = float(check["expected_ratio"])
    tolerance = float(check["relative_tolerance"])

    ckey = check.get("criteria_key")
    if ckey:
        section = criteria.get(scenario_name, {})
        conv = section.get("convergence", {})
        if isinstance(conv, dict) and ckey in conv and isinstance(conv[ckey], dict):
            citem = conv[ckey]
            expected = float(citem.get("expected", expected))
            tolerance = float(citem.get("relative_tolerance", tolerance))

    return expected, tolerance


def resolve_slope(check: Dict[str, Any], scenario_name: str, criteria: Dict[str, Any]) -> Tuple[float, float]:
    expected = float(check["expected_slope"])
    tolerance = float(check["slope_tolerance"])
    ckey = check.get("criteria_key")
    if ckey:
        section = criteria.get(scenario_name, {})
        slopes = section.get("loglog_slope", {})
        if isinstance(slopes, dict) and ckey in slopes and isinstance(slopes[ckey], dict):
            sitem = slopes[ckey]
            expected = float(sitem.get("expected", expected))
            tolerance = float(sitem.get("tolerance", tolerance))
    return expected, tolerance


def substitute_vars(command: str, variables: Dict[str, str]) -> str:
    for key, value in variables.items():
        command = command.replace("{" + key + "}", value)
    return command


def parse_key_value(items: List[str]) -> Dict[str, str]:
    parsed = {}
    for item in items:
        if "=" not in item:
            raise ValueError(f"Invalid --var entry '{item}'. Use key=value format.")
        key, value = item.split("=", 1)
        parsed[key.strip()] = value.strip()
    return parsed


def parse_command_options(command: str) -> Dict[str, float]:
    tokens = shlex.split(command)
    parsed: Dict[str, float] = {}
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if tok in {"-s", "-r", "-t", "-o"} and i + 1 < len(tokens):
            try:
                parsed[{"-s": "dt_soft", "-r": "rout", "-t": "tend", "-o": "dt_out"}[tok]] = float(tokens[i + 1])
            except ValueError:
                pass
            i += 2
            continue
        if tok in {"--r-group", "--r-search-group", "--tt-switch", "--tt-nstep", "--r-ratio"} and i + 1 < len(tokens):
            key_map = {
                "--r-group": "r_group",
                "--r-search-group": "r_search_group",
                "--tt-switch": "tt_switch",
                "--tt-nstep": "tt_nstep",
                "--r-ratio": "r_ratio",
            }
            try:
                parsed[key_map[tok]] = float(tokens[i + 1])
            except ValueError:
                pass
            i += 2
            continue
        i += 1
    if "rout" in parsed:
        r_ratio = parsed.get("r_ratio", 0.1)
        parsed["rin"] = r_ratio * parsed["rout"]
    return parsed


def parse_time_counts(log_text: str) -> List[Dict[str, float]]:
    rows: List[Dict[str, float]] = []
    for line in log_text.splitlines():
        match = COUNT_RE.match(line.strip())
        if not match:
            continue
        rows.append(
            {
                "time": float(match.group(1)),
                "n_real_loc": float(match.group(2)),
                "n_real_glb": float(match.group(3)),
                "n_all_loc": float(match.group(4)),
                "n_all_glb": float(match.group(5)),
            }
        )
    return rows


def parse_energy_rows(log_text: str) -> List[Dict[str, float]]:
    rows = []
    lines = log_text.splitlines()
    header = None
    for idx, line in enumerate(lines):
        head_match = ENERGY_HEADER_RE.match(line.strip())
        if head_match:
            header = head_match.group(1).split()
            if idx + 1 >= len(lines):
                continue
            phys_line = lines[idx + 1].strip()
            phys_match = ENERGY_PHYSIC_RE.match(phys_line)
            if not phys_match:
                continue
            parts = phys_match.group(1).split()
            if len(parts) < len(header):
                continue
            values = [float(parts[i]) for i in range(len(header))]
            rows.append(dict(zip(header, values)))
    return rows


def parse_time_error_series(log_text: str, metric: str = "Error/Total") -> List[Tuple[float, float]]:
    lines = log_text.splitlines()
    series: List[Tuple[float, float]] = []
    current_time = None

    for idx, line in enumerate(lines):
        tmatch = TIME_RE.match(line.strip())
        if tmatch:
            try:
                current_time = float(tmatch.group(1))
            except ValueError:
                current_time = None
            continue

        hmatch = ENERGY_HEADER_RE.match(line.strip())
        if hmatch and idx + 1 < len(lines) and current_time is not None:
            pmatch = ENERGY_PHYSIC_RE.match(lines[idx + 1].strip())
            if not pmatch:
                continue
            header = hmatch.group(1).split()
            values = pmatch.group(1).split()
            if len(values) < len(header):
                continue
            rec = {header[j]: float(values[j]) for j in range(len(header))}
            series.append((current_time, abs(rec.get(metric, 0.0))))

    return series


def period_maxima(
    series: List[Tuple[float, float]],
    period: float,
    n_period: int,
) -> List[float]:
    if period <= 0.0 or n_period <= 0:
        return []

    period_max: Dict[int, float] = {}
    t_upper = period * n_period + 1e-12

    for t, value in series:
        if t < 0.0 or t > t_upper:
            continue
        idx = int(t / period)
        if idx >= n_period:
            idx = n_period - 1
        if idx < 0:
            continue
        period_max[idx] = max(period_max.get(idx, 0.0), value)

    return [period_max[k] for k in sorted(period_max.keys())]


def resolve_status_path_from_command(command: str, workdir: Path) -> Path | None:
    tokens = shlex.split(command)
    prefix = None
    run_dir = workdir

    for i, tok in enumerate(tokens):
        if tok == "-f" and i + 1 < len(tokens):
            prefix = tokens[i + 1]
            break

    if prefix is None:
        return None

    if len(tokens) >= 3 and tokens[0] == "cd" and tokens[2] == "&&":
        p = Path(tokens[1])
        run_dir = p if p.is_absolute() else (workdir / p)

    return run_dir.resolve() / f"{prefix}.status"


def extract_orbital_drift_from_status(status_path: Path, repo_root: Path, n_particle: int, interrupt_mode: str = "none") -> Dict[str, float]:
    if n_particle < 2:
        return {}

    tools_path = (repo_root / "tools").resolve()
    if str(tools_path) not in sys.path:
        sys.path.insert(0, str(tools_path))

    from analysis.status import Status  # pylint: disable=import-outside-toplevel  # pyright: ignore[reportMissingImports]

    # interrupt_mode must match the producing solver's feature set, otherwise
    # the binary status dtype misaligns (bse columns shift every record)
    st = Status(N_particle=n_particle, interrupt_mode=interrupt_mode)
    with warnings.catch_warnings(record=True) as warn_list:
        warnings.simplefilter("always", category=UserWarning)
        st.fromfile(str(status_path))
    if warn_list:
        msg = "; ".join(str(item.message) for item in warn_list)
        raise RuntimeError(
            f"Status parse warning for {status_path} with N_particle={n_particle}: {msg}. "
            "This usually indicates an N_particle mismatch."
        )

    p0 = st.particles.p0
    p1 = st.particles.p1
    m1 = float(p0.mass[0])
    m2 = float(p1.mass[0])
    mu = G_MSUN_PC_MYR * (m1 + m2)

    r = p1.pos - p0.pos
    v = p1.vel - p0.vel
    rnorm = np.linalg.norm(r, axis=1)
    v2 = np.sum(v * v, axis=1)
    eps = 0.5 * v2 - mu / rnorm
    a = -mu / (2.0 * eps)
    h = np.cross(r, v)
    h2 = np.sum(h * h, axis=1)
    e = np.sqrt(np.maximum(0.0, 1.0 + 2.0 * eps * h2 / (mu * mu)))

    a0 = float(a[0])
    e0 = float(e[0])
    da_rel = np.abs((a - a0) / a0)
    de_abs = np.abs(e - e0)

    return {
        "max_rel_drift_semi": float(np.max(da_rel)),
        "max_abs_drift_ecc": float(np.max(de_abs)),
        "status_sample_count": int(st.size),
    }


def extract_momentum_from_status(status_path: Path, repo_root: Path, n_particle: int, interrupt_mode: str = "none") -> Dict[str, float]:
    """Total linear momentum drift from per-particle status rows (-w 2).

    The tidal-tensor mean-leak bug (fixed 2026-10-06) manifests as a secular,
    one-directional drift of sum(m*v); the fix caps it at the baseline
    integration-drift level. Returns drift normalized to velocity units
    (pc/Myr) so thresholds are scale-free.
    """
    if n_particle < 2:
        return {}

    tools_path = (repo_root / "tools").resolve()
    if str(tools_path) not in sys.path:
        sys.path.insert(0, str(tools_path))

    from analysis.status import Status  # pylint: disable=import-outside-toplevel  # pyright: ignore[reportMissingImports]

    st = Status(N_particle=n_particle, interrupt_mode=interrupt_mode)
    with warnings.catch_warnings(record=True) as warn_list:
        warnings.simplefilter("always", category=UserWarning)
        st.fromfile(str(status_path))
    if warn_list:
        msg = "; ".join(str(item.message) for item in warn_list)
        raise RuntimeError(
            f"Status parse warning for {status_path} with N_particle={n_particle}: {msg}. "
            "This usually indicates an N_particle mismatch."
        )

    px = np.zeros(st.size)
    py = np.zeros(st.size)
    pz = np.zeros(st.size)
    m_tot = np.zeros(st.size)
    for k in range(n_particle):
        pk = getattr(st.particles, f"p{k}")
        px += pk.mass * pk.vel[:, 0]
        py += pk.mass * pk.vel[:, 1]
        pz += pk.mass * pk.vel[:, 2]
        m_tot += pk.mass

    p_norm = np.sqrt(px * px + py * py + pz * pz)
    drift = np.max(np.abs(p_norm - p_norm[0])) / m_tot[0]

    return {
        "max_momentum_drift_pcm": float(drift),
    }


def extract_regex_counts(log_text: str, patterns: Dict[str, str]) -> Dict[str, int]:
    counts = {}
    for key, pattern in patterns.items():
        counts[key] = len(re.findall(pattern, log_text, flags=re.MULTILINE))
    return counts


def extract_regex_floats(log_text: str, patterns: Dict[str, str]) -> Dict[str, float]:
    """First capture group of the LAST match per pattern (final report line)."""
    values: Dict[str, float] = {}
    for key, pattern in patterns.items():
        matches = re.findall(pattern, log_text, flags=re.MULTILINE)
        if not matches:
            continue
        last = matches[-1]
        text = last[0] if isinstance(last, tuple) else last
        try:
            values[key] = float(text)
        except ValueError:
            continue
    return values


def extract_metrics(output_path: Path, extract_cfg: Dict[str, Any], command: str, workdir: Path, repo_root: Path) -> Dict[str, Any]:
    text = output_path.read_text(encoding="utf-8", errors="replace")
    energy_rows = parse_energy_rows(text)
    count_rows = parse_time_counts(text)
    max_err_total = 0.0
    max_err_pp = 0.0
    if energy_rows:
        max_err_total = max(abs(r.get("Error/Total", 0.0)) for r in energy_rows)
        max_err_pp = max(abs(r.get("Error_PP", 0.0)) for r in energy_rows)

    regex_cfg = extract_cfg.get("regex_counts", {})
    regex_counts = extract_regex_counts(text, regex_cfg) if regex_cfg else {}

    floats_cfg = extract_cfg.get("regex_floats", {})
    regex_floats = extract_regex_floats(text, floats_cfg) if floats_cfg else {}

    period = float(extract_cfg.get("period", 0.0) or 0.0)
    n_period = int(extract_cfg.get("n_period", 0) or 0)
    period_metric = str(extract_cfg.get("period_metric", "Error/Total"))
    series = parse_time_error_series(text, metric=period_metric)
    pmax = period_maxima(series, period=period, n_period=n_period) if (period > 0.0 and n_period > 0) else []

    max_period_error = max(pmax) if pmax else 0.0
    sum_period_error = sum(pmax) if pmax else 0.0
    mean_period_error = (sum_period_error / len(pmax)) if pmax else 0.0
    # stationarity: max of the second half of the period maxima over the first
    # half (≈1 stationary; >>1 secular drift, <<1 anomalous decay)
    if len(pmax) >= 2:
        half = len(pmax) // 2
        early_half = max(pmax[:half]) if half > 0 else pmax[0]
        late_half = max(pmax[half:]) if pmax[half:] else early_half
        late_early_period_max_ratio = late_half / early_half if early_half > 0.0 else float("inf")
    else:
        late_early_period_max_ratio = float("nan")

    orbital_metrics: Dict[str, float] = {}
    status_n_particle = int(extract_cfg.get("status_n_particle", 0) or 0)
    if status_n_particle >= 2:
        status_path = resolve_status_path_from_command(command, workdir)
        if status_path is not None and status_path.exists():
            orbital_metrics = extract_orbital_drift_from_status(status_path, repo_root=repo_root, n_particle=status_n_particle, interrupt_mode=str(extract_cfg.get("interrupt_mode", "none")))
            if extract_cfg.get("momentum", False):
                orbital_metrics.update(
                    extract_momentum_from_status(
                        status_path,
                        repo_root=repo_root,
                        n_particle=status_n_particle,
                        interrupt_mode=str(extract_cfg.get("interrupt_mode", "none")),
                    )
                )

    # AR group census from the final statistics block (cumulative counters):
    # header " AR_tsyn_sum Stab_Ngroup Iso_Ngroup Form_Ngroup End_Ngroup ..."
    # NOTE: these profile counters reset per hard-integration segment, so they
    # only reflect the last segment; the authoritative cumulative census is
    # parsed from the group event dumps (group_dump_* below).
    ar_form_ngroup = -1
    ar_end_ngroup = -1
    log_lines = text.splitlines()
    for i, line in enumerate(log_lines):
        if line.startswith(" AR_tsyn_sum Stab_Ngroup"):
            parts = log_lines[i + 1].split()
            if len(parts) >= 5:
                ar_form_ngroup = int(float(parts[3]))
                ar_end_ngroup = int(float(parts[4]))

    # Cumulative group form/break events from the ADJUST_GROUP_PRINT binary
    # dumps <prefix>.group.n<m> next to the status file. Record header:
    # [int32 type (0=form, 1=break)][int32 n_member][float64 time]...
    # Records have fixed size per file (n_member fixed), so candidates are
    # filtered onto the dominant stride grid to reject payload bit patterns
    # that mimic the header.
    group_form_events = -1
    group_break_events = -1
    status_path = resolve_status_path_from_command(command, workdir)
    if status_path is not None:
        group_files = sorted(status_path.parent.glob(status_path.stem + ".group.n*"))
        n_form = 0
        n_break = 0
        for gfile in group_files:
            try:
                n_member_expected = int(gfile.name.rsplit(".n", 1)[1])
                raw = gfile.read_bytes()
            except (OSError, ValueError):
                continue
            a = np.frombuffer(raw[: len(raw) // 4 * 4], dtype=np.int32)
            cand = np.where(((a[:-1] == 0) | (a[:-1] == 1)))[0]
            cand = cand[(a[cand + 1] == n_member_expected) & (cand % 2 == 0)]
            keep = []
            for i in cand:
                off = i * 4 + 8
                if off + 8 > len(raw):
                    continue
                tval = struct.unpack_from("<d", raw, off)[0]
                if 0.0 <= tval <= 1.0e3:
                    keep.append(int(i))
            if len(keep) >= 3:
                diffs = np.diff(keep)
                stride = int(np.median(diffs))
                if stride > 0:
                    base = keep[0]
                    keep = [i for i in keep if (i - base) % stride == 0]
            for i in keep:
                if a[i] == 0:
                    n_form += 1
                else:
                    n_break += 1
        # no group dump files at all = no group events (petar deletes empty
        # dumps on exit); only an unresolvable status path keeps the -1 sentinel
        if group_files or status_path is not None:
            group_form_events = n_form
            group_break_events = n_break

    max_n_real_glb = max((row["n_real_glb"] for row in count_rows), default=0.0)
    max_n_all_glb = max((row["n_all_glb"] for row in count_rows), default=0.0)
    max_n_real_loc = max((row["n_real_loc"] for row in count_rows), default=0.0)
    max_n_all_loc = max((row["n_all_loc"] for row in count_rows), default=0.0)
    max_artificial_particles_glb = max((row["n_all_glb"] - row["n_real_glb"] for row in count_rows), default=0.0)
    max_artificial_particles_loc = max((row["n_all_loc"] - row["n_real_loc"] for row in count_rows), default=0.0)
    final_time = count_rows[-1]["time"] if count_rows else 0.0

    return {
        "max_abs_error_over_total": max_err_total,
        "max_abs_error_pp": max_err_pp,
        "energy_sample_count": len(energy_rows),
        "time_sample_count": len(count_rows),
        "max_n_real_glb": max_n_real_glb,
        "max_n_all_glb": max_n_all_glb,
        "max_n_real_loc": max_n_real_loc,
        "max_n_all_loc": max_n_all_loc,
        "max_artificial_particles_glb": max_artificial_particles_glb,
        "max_artificial_particles_loc": max_artificial_particles_loc,
        "ar_form_ngroup": float(ar_form_ngroup),
        "ar_end_ngroup": float(ar_end_ngroup),
        "group_form_events": float(group_form_events),
        "group_break_events": float(group_break_events),
        "final_time_log": final_time,
        "regex_counts": regex_counts,
        "regex_floats": regex_floats,
        **regex_floats,
        "period_bin_count": len(pmax),
        "max_abs_error_period_max": max_period_error,
        "sum_abs_error_period_max": sum_period_error,
        "mean_abs_error_period_max": mean_period_error,
        "late_early_period_max_ratio": late_early_period_max_ratio,
        **orbital_metrics,
    }


def check_monotonic_trend(
    values: List[float],
    direction: str,
    relative_tolerance: float,
) -> Dict[str, Any]:
    if len(values) < 2:
        return {
            "passed": False,
            "message": "need at least two points for monotonic trend check",
        }

    if direction not in {"nondecreasing", "nonincreasing"}:
        return {
            "passed": False,
            "message": f"unknown monotonic direction: {direction}",
        }

    for i in range(1, len(values)):
        prev_v = values[i - 1]
        curr_v = values[i]
        if direction == "nondecreasing":
            lower = prev_v * (1.0 - relative_tolerance)
            if curr_v < lower:
                return {
                    "passed": False,
                    "message": f"index {i}: {curr_v:.6e} < allowed lower {lower:.6e} (prev={prev_v:.6e})",
                }
        else:
            upper = prev_v * (1.0 + relative_tolerance)
            if curr_v > upper:
                return {
                    "passed": False,
                    "message": f"index {i}: {curr_v:.6e} > allowed upper {upper:.6e} (prev={prev_v:.6e})",
                }

    return {
        "passed": True,
        "message": f"{direction} trend satisfied with relative_tolerance={relative_tolerance:.3f}",
    }


def run_shell(command: str, cwd: Path, output_path: Path) -> int:
    with output_path.open("w", encoding="utf-8") as out:
        proc = subprocess.run(command, cwd=str(cwd), shell=True, stdout=out, stderr=subprocess.STDOUT)
    return proc.returncode


def _binary_provenance(variables: Dict[str, str]) -> Dict[str, Any]:
    """md5 + Build banner of every petar binary referenced by the run variables.

    Provenance is recorded per report so that error levels can always be traced
    back to the exact binaries that produced them (build-flag sensitivity made
    this necessary; see NOTE 2026-10-12/13 9b)."""
    import hashlib
    import subprocess as _sp

    seen: Dict[str, Dict[str, Any]] = {}
    for value in variables.values():
        cand = Path(value)
        try:
            if not cand.is_file() or not (cand.stat().st_mode & 0o111):
                continue
        except OSError:
            continue
        key = str(cand.resolve())
        if key in seen:
            continue
        entry: Dict[str, Any] = {"path": key}
        try:
            entry["md5"] = hashlib.md5(cand.read_bytes()).hexdigest()
        except OSError:
            pass
        try:
            r = _sp.run([str(cand), "-h"], capture_output=True, text=True, timeout=60)
            for line in (r.stdout + r.stderr).splitlines():
                if line.startswith("Build:"):
                    entry["build"] = line.strip()
                    break
        except Exception:
            pass
        seen[key] = entry
    return {"binaries": list(seen.values())}


def execute_scenario(
    scenario: Dict[str, Any],
    base_dir: Path,
    out_dir: Path,
    variables: Dict[str, str],
    criteria: Dict[str, Any],
    dry_run: bool,
) -> Tuple[List[RunResult], List[Dict[str, Any]]]:
    scenario_name = scenario["name"]
    workdir = Path(substitute_vars(scenario.get("workdir", "."), variables))
    if not workdir.is_absolute():
        workdir = (base_dir / workdir).resolve()

    scenario_out = out_dir / f"validation_{scenario_name}"
    scenario_out.mkdir(parents=True, exist_ok=True)

    setup_cmds = scenario.get("setup_commands", [])
    for idx, cmd in enumerate(setup_cmds):
        resolved = substitute_vars(cmd, variables)
        if dry_run:
            print(f"[dry-run] setup[{idx}] {resolved}")
            continue
        code = run_shell(resolved, workdir, scenario_out / f"setup.{idx}.log")
        if code != 0:
            raise RuntimeError(f"Setup command failed in {scenario_name}: {resolved}")

    run_results: List[RunResult] = []
    for run in scenario["runs"]:
        run_id = run["id"]
        cmd = substitute_vars(run["command"], variables)
        output_name = run.get("output", f"{run_id}.log")
        output_path = scenario_out / output_name

        if dry_run:
            print(f"[dry-run] run {scenario_name}/{run_id}: {cmd}")
            metrics = parse_command_options(cmd)
        else:
            code = run_shell(cmd, workdir, output_path)
            if code != 0:
                raise RuntimeError(f"Run command failed ({scenario_name}/{run_id}): {cmd}")
            metrics = extract_metrics(
                output_path,
                scenario.get("extract", {}),
                command=cmd,
                workdir=workdir,
                repo_root=base_dir,
            )
            metrics.update(parse_command_options(cmd))
            if "dt_soft" in metrics:
                dt_soft = float(metrics["dt_soft"])
                metrics["dt_soft_times_max_abs_error_period_max"] = dt_soft * float(metrics.get("max_abs_error_period_max", 0.0))
                metrics["dt_soft_times_sum_abs_error_period_max"] = dt_soft * float(metrics.get("sum_abs_error_period_max", 0.0))

        run_results.append(RunResult(run_id=run_id, command=cmd, output_path=output_path, metrics=metrics))

    check_results = []
    if not dry_run:
        per_run = {item.run_id: item.metrics for item in run_results}
        for check in scenario.get("checks", []):
            ctype = check["type"]
            skip = False
            req = check.get("require_var_contains")
            if req and req[1] not in str(variables.get(req[0], "")):
                skip = True
            reqn = check.get("require_var_not_contains")
            if reqn and reqn[1] in str(variables.get(reqn[0], "")):
                skip = True
            if skip:
                check_results.append(
                    {
                        "passed": True,
                        "skipped": True,
                        "scenario": scenario_name,
                        "check": check["name"],
                        "message": "skipped: build condition not met (require_var_[not_]contains)",
                    }
                )
                continue
            if ctype == "max_threshold":
                run_id = check["run_id"]
                metric = check["metric"]
                threshold = resolve_threshold(
                    check=check,
                    scenario_name=scenario_name,
                    criteria=criteria,
                    metric=metric,
                    key=run_id,
                    fallback=float(check["threshold"]),
                )
                value = float(per_run[run_id][metric])
                result = check_max_threshold(value, threshold)
                result.update({"scenario": scenario_name, "check": check["name"]})
                check_results.append(result)
            elif ctype == "min_threshold":
                run_id = check["run_id"]
                metric = check["metric"]
                threshold = float(check["threshold"])
                value = float(per_run[run_id][metric])
                result = check_min_threshold(value, threshold)
                result.update({"scenario": scenario_name, "check": check["name"]})
                check_results.append(result)
            elif ctype == "convergence_ratio":
                skip_vars = check.get("skip_if_vars_equal", [])
                if len(skip_vars) == 2:
                    var_a = variables.get(skip_vars[0], "")
                    var_b = variables.get(skip_vars[1], "")
                    if var_a and var_b and var_a == var_b:
                        check_results.append(
                            {
                                "passed": True,
                                "skipped": True,
                                "scenario": scenario_name,
                                "check": check["name"],
                                "message": f"skipped: {skip_vars[0]} and {skip_vars[1]} are identical ({var_a})",
                            }
                        )
                        continue

                coarse_id = check["coarse_run_id"]
                fine_id = check["fine_run_id"]
                metric = check["metric"]
                expected_ratio, tolerance = resolve_convergence(check, scenario_name, criteria)
                coarse = float(per_run[coarse_id][metric])
                fine = float(per_run[fine_id][metric])

                floor = check.get("min_fine_error")
                if floor is not None and abs(fine) <= float(floor):
                    check_results.append(
                        {
                            "passed": True,
                            "skipped": True,
                            "scenario": scenario_name,
                            "check": check["name"],
                            "message": f"skipped: fine error {fine:.6e} <= floor {float(floor):.6e}",
                        }
                    )
                    continue

                result = check_convergence_ratio(coarse, fine, expected_ratio, tolerance)
                result.update({"scenario": scenario_name, "check": check["name"]})
                check_results.append(result)
            elif ctype == "regex_count_max":
                run_id = check["run_id"]
                regex_key = check["regex_key"]
                threshold = int(
                    resolve_threshold(
                        check=check,
                        scenario_name=scenario_name,
                        criteria=criteria,
                        metric="regex_count_max",
                        key=regex_key,
                        fallback=float(check["threshold"]),
                    )
                )
                value = int(per_run[run_id]["regex_counts"].get(regex_key, 0))
                result = check_max_threshold(float(value), float(threshold))
                result.update({"scenario": scenario_name, "check": check["name"], "value": value, "threshold": threshold})
                check_results.append(result)
            elif ctype == "metric_ratio":
                numerator_id = check["numerator_run_id"]
                denominator_id = check["denominator_run_id"]
                metric = check["metric"]
                missing = [rid for rid in (numerator_id, denominator_id) if rid not in per_run]
                if missing:
                    check_results.append(
                        {
                            "passed": False,
                            "scenario": scenario_name,
                            "check": check["name"],
                            "message": f"missing metrics for runs: {', '.join(missing)}",
                        }
                    )
                    continue
                numerator = float(per_run[numerator_id][metric])
                denominator = float(per_run[denominator_id][metric])
                result = check_metric_ratio(
                    numerator, denominator, float(check["min_ratio"]), float(check["max_ratio"])
                )
                result.update({"scenario": scenario_name, "check": check["name"]})
                check_results.append(result)
            elif ctype == "loglog_slope":
                run_ids = check["run_ids"]
                metric = check["metric"]
                x_metric = check.get("x_metric", "dt_soft")
                x_values = []
                y_values = []
                missing = []
                for rid in run_ids:
                    run_rec = next((r for r in run_results if r.run_id == rid), None)
                    if run_rec is None:
                        missing.append(rid)
                        continue
                    x_val = run_rec.metrics.get(x_metric)
                    y_val = run_rec.metrics.get(metric)
                    if x_val is None or y_val is None:
                        missing.append(rid)
                        continue
                    x_values.append(float(x_val))
                    y_values.append(float(y_val))
                if missing:
                    check_results.append(
                        {
                            "passed": False,
                            "scenario": scenario_name,
                            "check": check["name"],
                            "message": f"missing metrics for runs: {', '.join(missing)}",
                        }
                    )
                    continue
                expected, tolerance = resolve_slope(check, scenario_name, criteria)
                result = check_loglog_slope(
                    x_values, y_values, expected, tolerance, bound=str(check.get("bound", "both"))
                )
                result.update({"scenario": scenario_name, "check": check["name"]})
                check_results.append(result)
            elif ctype == "monotonic_trend":
                run_ids = check["run_ids"]
                metric = check["metric"]
                direction = str(check.get("direction", "nondecreasing"))
                relative_tolerance = float(check.get("relative_tolerance", 0.0))

                values = []
                missing = []
                for rid in run_ids:
                    run_rec = next((r for r in run_results if r.run_id == rid), None)
                    if run_rec is None:
                        missing.append(rid)
                        continue
                    mv = run_rec.metrics.get(metric)
                    if mv is None:
                        missing.append(rid)
                        continue
                    values.append(float(mv))

                if missing:
                    check_results.append(
                        {
                            "passed": False,
                            "scenario": scenario_name,
                            "check": check["name"],
                            "message": f"missing metrics for runs: {', '.join(missing)}",
                        }
                    )
                    continue

                result = check_monotonic_trend(values, direction=direction, relative_tolerance=relative_tolerance)
                result.update({"scenario": scenario_name, "check": check["name"]})
                check_results.append(result)
            else:
                raise ValueError(f"Unknown check type: {ctype}")

    return run_results, check_results


def _detect_petar_binary() -> str:
    """Try to auto-detect the active petar binary, preferring KDKDK4 variants."""
    import shutil
    # 1) if petar.select is available, try to get a KDKDK4 binary
    sel = shutil.which("petar.select")
    if sel is not None:
        import subprocess as _sp
        try:
            result = _sp.run(
                [sel, "--require", "kdkdk4", "--optional", "mpi,omp,avx2"],
                capture_output=True, text=True, timeout=10,
            )
            if result.returncode == 0:
                line = result.stdout.strip().split("\n")[0]
                if line:
                    path = shutil.which(line.split()[0] if " " not in line else line)
                    if path:
                        return path
        except Exception:
            pass
    # 2) check the current petar symlink
    petar_path = shutil.which("petar")
    if petar_path is not None:
        resolved = Path(petar_path).resolve()
        parent = resolved.parent
        stem = resolved.name
        # try kdkdk4 suffix variants
        for suffix in ["kdkdk4", "64b.kdkdk4"]:
            candidate = parent / f"{stem}.{suffix}"
            if candidate.exists():
                return str(candidate)
        # try sibling kdkdk4 binaries with similar prefix
        for suffix in ["kdkdk4", "64b.kdkdk4"]:
            for candidate in sorted(parent.glob(f"petar.*.{suffix}")):
                return str(candidate)
        return petar_path
    # 3) fallback: check bin/ relative to repo root
    candidates = [
        Path(__file__).resolve().parents[2] / "bin" / "petar",
    ]
    for c in candidates:
        if c.exists():
            return str(c)
    return "petar"


def _validate_variables(variables: Dict[str, str], scenario: Dict[str, Any]) -> None:
    """Check that all template variables used in scenario commands are defined.

    Only matches {var} patterns where var contains at least one underscore,
    to avoid false positives like awk's {print}.
    """
    import re as _re
    # require at least one underscore — real template vars are like {petar_bin_switch}, {repo_root}
    pattern = _re.compile(r"\{([a-zA-Z_][a-zA-Z0-9_]*_[a-zA-Z0-9_.]+)}")
    for run in scenario.get("runs", []):
        cmd = run.get("command", "")
        for m in pattern.findall(cmd):
            if m not in variables:
                raise ValueError(f"Undefined template variable '{{{m}}}' in run '{run.get('id')}'. "
                                 f"Use --var {m}=<value> to define it.")
    for cmd in scenario.get("setup_commands", []):
        for m in pattern.findall(cmd):
            if m not in variables:
                raise ValueError(f"Undefined template variable '{{{m}}}' in setup command. "
                                 f"Use --var {m}=<value> to define it.")


def _select_or_build_petar(require_64b: bool | None = None) -> str:
    """Select a petar binary via petar.select; reinstall when the selection
    fails or the installed binary version does not match the repository
    VERSION (see report_common.select_or_install_petar). Never falls back to
    hardcoded binary names or the default `petar` symlink.
    """
    from report_common import select_or_install_petar

    repo_root = Path(__file__).resolve().parents[2]
    if require_64b:
        return select_or_install_petar(require=["64b"], optional=["avx2", "omp"], repo_root=repo_root)
    return select_or_install_petar(optional=["avx2", "omp"], repo_root=repo_root)


def _pipeline(scenario_file: str, report_json: str, out_dir: str, petar_var: str, petar_bin: str) -> int:
    """Run a single scenario via run_validation.py subprocess."""
    cmd = [
        sys.executable, __file__,
        "--scenario", scenario_file,
        "--out-dir", out_dir,
        "--report", report_json,
        "--var", f"{petar_var}={petar_bin}",
    ]
    print(f"[pipeline] {' '.join(cmd)}")
    return subprocess.run(cmd, check=False).returncode


def main() -> int:
    parser = argparse.ArgumentParser(description="PeTar validation runner")
    parser.add_argument("--scenario", action="append", help="Scenario JSON path (can be repeated)")
    parser.add_argument("--scenario-dir", default="test/validation/scenarios", help="Default scenario directory")
    parser.add_argument("--out-dir", default="test/out", help="Output directory")
    parser.add_argument("--var", action="append", default=[], help="Template variables key=value")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running")
    parser.add_argument("--report", default="test/out/report.json", help="Report JSON path")
    parser.add_argument("--criteria", default="test/validation/criteria.json", help="Criteria JSON path")
    parser.add_argument("--continue-on-error", action="store_true",
                        help="Continue with remaining scenarios/runs even if one fails")
    parser.add_argument("--pipeline", choices=["t1", "t2", "t3"],
                        help="Run full pipeline: select/build binary, run scenario, generate HTML report")
    args = parser.parse_args()

    # --pipeline mode: single-test end-to-end with HTML report
    if args.pipeline == "t1":
        bin_64 = str(Path(_select_or_build_petar(True)).resolve())
        bin_non64 = str(Path(_select_or_build_petar(False)).resolve())
        print(f"[pipeline] 64b: {bin_64}")
        sfile = "test/validation/scenarios/t1_high_ecc_changeover.json"
        _pipeline(sfile, "test/out/report.t1.64b.overlay.json", "test/out/validation_t1_binary_changeover_64b", "petar_bin_order4", bin_64)
        _pipeline(sfile, "test/out/report.t1.non64b.overlay.json", "test/out/validation_t1_binary_changeover_non64b", "petar_bin_order4", bin_non64)
        subprocess.run([
            sys.executable, "test/validation/t1_binary_changeover_kdkdk4_report.py",
            "--report", "test/out/report.t1.64b.overlay.json",
            "--primary-tag", "64b",
            "--compare-report", "test/out/report.t1.non64b.overlay.json",
            "--compare-tag", "non64b",
            "--output", "test/out/t1_binary_changeover_kdkdk4_summary.html",
        ], check=False)
        return 0
    if args.pipeline == "t2":
        bin_pet = _select_or_build_petar()
        print(f"[pipeline] t2 binary: {bin_pet}")
        sfile = "test/validation/scenarios/t2_binary_conservation_longterm.json"
        _pipeline(sfile, "test/out/report.t2.binary.json", "test/out/validation_t2_binary", "petar_bin_switch", bin_pet)
        subprocess.run([
            sys.executable, "test/validation/t2_binary_longterm_autochangeover_report.py",
            "--report", "test/out/report.t2.binary.json",
            "--output", "test/out/t2_binary_longterm_autochangeover_summary.html",
        ], check=False)
        return 0
    if args.pipeline == "t3":
        bin_pet = _select_or_build_petar()
        sfile = "test/validation/scenarios/t3_binary_hard_switch_longterm.json"
        _pipeline(sfile, "test/out/report.t3.binary.json", "test/out/validation_t3_binary", "petar_bin_switch", bin_pet)
        subprocess.run([
            sys.executable, "test/validation/t3_binary_sdar_grouplink_report.py",
            "--report", "test/out/report.t3.binary.json",
            "--output", "test/out/t3_binary_sdar_grouplink_summary.html",
        ], check=False)
        return 0

    repo_root = Path(__file__).resolve().parents[2]
    scenario_dir = (repo_root / args.scenario_dir).resolve()
    out_dir = (repo_root / args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    criteria_path = Path(args.criteria)
    if not criteria_path.is_absolute():
        criteria_path = (repo_root / args.criteria).resolve()
    criteria = load_criteria(criteria_path)

    variables = parse_key_value(args.var)
    variables.setdefault("repo_root", str(repo_root))
    detected = _detect_petar_binary()
    variables.setdefault("petar_bin_order2", detected)
    variables.setdefault("petar_bin_order4", detected)
    variables.setdefault("petar_bin_switch", detected)
    variables.setdefault("petar_bin_blogh_a", detected)
    variables.setdefault("petar_bin_blogh_b", detected)

    scenario_paths = []
    if args.scenario:
        for item in args.scenario:
            p = Path(item)
            if not p.is_absolute():
                p = (repo_root / item).resolve()
            scenario_paths.append(p)
    else:
        scenario_paths = sorted(scenario_dir.glob("*.json"))

    if not scenario_paths:
        print("No scenario files found.", file=sys.stderr)
        return 2

    all_checks: List[Dict[str, Any]] = []
    all_runs: List[Dict[str, Any]] = []

    for spath in scenario_paths:
        scenario = load_json(spath)
        if not args.dry_run:
            _validate_variables(variables, scenario)
        try:
            run_results, check_results = execute_scenario(
                scenario=scenario,
                base_dir=repo_root,
                out_dir=out_dir,
                variables=variables,
                criteria=criteria,
                dry_run=args.dry_run,
            )
        except (RuntimeError, ValueError) as exc:
            print(f"[ERROR] Scenario '{scenario.get('name', '?')}' failed: {exc}", file=sys.stderr)
            if not args.continue_on_error:
                raise
            run_results = []
            check_results = [{
                "passed": False,
                "scenario": scenario.get("name", "?"),
                "check": "scenario_execution",
                "message": str(exc),
            }]

        for run_item in run_results:
            all_runs.append(
                {
                    "scenario": scenario["name"],
                    "run_id": run_item.run_id,
                    "command": run_item.command,
                    "output": str(run_item.output_path),
                    "metrics": run_item.metrics,
                }
            )
        all_checks.extend(check_results)

    summary = summarize_results(all_checks) if not args.dry_run else {"passed": True, "n_total": 0, "n_failed": 0}
    report = {
        "dry_run": args.dry_run,
        "petar_provenance": _binary_provenance(variables),
        "summary": summary,
        "runs": all_runs,
        "checks": all_checks,
    }

    report_path = Path(args.report)
    if not report_path.is_absolute():
        report_path = (repo_root / report_path).resolve()
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    if args.dry_run:
        print(f"Dry run finished. Report: {report_path}")
        return 0

    for item in all_checks:
        status = "PASS" if item["passed"] else "FAIL"
        print(f"[{status}] {item['scenario']}::{item['check']} - {item['message']}")

    print(
        f"Summary: passed={summary['passed']} total={summary['n_total']} failed={summary['n_failed']} report={report_path}"
    )
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())