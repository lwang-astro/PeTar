#!/usr/bin/env python3
import argparse
import json
import math
import os
import re
import shlex
import subprocess
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

from metrics import check_convergence_ratio, check_loglog_slope, check_max_threshold, check_min_threshold, summarize_results


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


def extract_orbital_drift_from_status(status_path: Path, repo_root: Path, n_particle: int) -> Dict[str, float]:
    if n_particle < 2:
        return {}

    tools_path = (repo_root / "tools").resolve()
    if str(tools_path) not in sys.path:
        sys.path.insert(0, str(tools_path))

    from analysis.status import Status  # pylint: disable=import-outside-toplevel  # pyright: ignore[reportMissingImports]

    st = Status(N_particle=n_particle)
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


def extract_regex_counts(log_text: str, patterns: Dict[str, str]) -> Dict[str, int]:
    counts = {}
    for key, pattern in patterns.items():
        counts[key] = len(re.findall(pattern, log_text, flags=re.MULTILINE))
    return counts


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

    period = float(extract_cfg.get("period", 0.0) or 0.0)
    n_period = int(extract_cfg.get("n_period", 0) or 0)
    period_metric = str(extract_cfg.get("period_metric", "Error/Total"))
    series = parse_time_error_series(text, metric=period_metric)
    pmax = period_maxima(series, period=period, n_period=n_period) if (period > 0.0 and n_period > 0) else []

    max_period_error = max(pmax) if pmax else 0.0
    sum_period_error = sum(pmax) if pmax else 0.0
    mean_period_error = (sum_period_error / len(pmax)) if pmax else 0.0

    orbital_metrics: Dict[str, float] = {}
    status_n_particle = int(extract_cfg.get("status_n_particle", 0) or 0)
    if status_n_particle >= 2:
        status_path = resolve_status_path_from_command(command, workdir)
        if status_path is not None and status_path.exists():
            orbital_metrics = extract_orbital_drift_from_status(status_path, repo_root=repo_root, n_particle=status_n_particle)

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
        "final_time_log": final_time,
        "regex_counts": regex_counts,
        "period_bin_count": len(pmax),
        "max_abs_error_period_max": max_period_error,
        "sum_abs_error_period_max": sum_period_error,
        "mean_abs_error_period_max": mean_period_error,
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
                result = check_loglog_slope(x_values, y_values, expected, tolerance)
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


def main() -> int:
    parser = argparse.ArgumentParser(description="PeTar validation runner (minimal framework)")
    parser.add_argument("--scenario", action="append", help="Scenario JSON path (can be repeated)")
    parser.add_argument("--scenario-dir", default="test/validation/scenarios", help="Default scenario directory")
    parser.add_argument("--out-dir", default="test/out", help="Output directory")
    parser.add_argument("--var", action="append", default=[], help="Template variables key=value")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running")
    parser.add_argument("--report", default="test/out/report.json", help="Report JSON path")
    parser.add_argument("--criteria", default="test/validation/criteria.json", help="Criteria JSON path")
    args = parser.parse_args()

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
    variables.setdefault("petar_bin_order2", "petar")
    variables.setdefault("petar_bin_order4", "petar")
    variables.setdefault("petar_bin_switch", "petar")
    variables.setdefault("petar_bin_blogh_a", variables.get("petar_bin_switch", "petar"))
    variables.setdefault("petar_bin_blogh_b", variables.get("petar_bin_switch", "petar"))

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
        run_results, check_results = execute_scenario(
            scenario=scenario,
            base_dir=repo_root,
            out_dir=out_dir,
            variables=variables,
            criteria=criteria,
            dry_run=args.dry_run,
        )

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