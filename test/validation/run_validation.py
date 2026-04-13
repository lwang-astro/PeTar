#!/usr/bin/env python3
import argparse
import json
import os
import re
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple

from metrics import check_convergence_ratio, check_max_threshold, summarize_results


ENERGY_HEADER_RE = re.compile(r"^Energy:\s+(.*)$")
ENERGY_PHYSIC_RE = re.compile(r"^Physic:\s+(.*)$")


@dataclass
class RunResult:
    run_id: str
    command: str
    output_path: Path
    metrics: Dict[str, Any]


def load_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


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


def extract_regex_counts(log_text: str, patterns: Dict[str, str]) -> Dict[str, int]:
    counts = {}
    for key, pattern in patterns.items():
        counts[key] = len(re.findall(pattern, log_text, flags=re.MULTILINE))
    return counts


def extract_metrics(output_path: Path, extract_cfg: Dict[str, Any]) -> Dict[str, Any]:
    text = output_path.read_text(encoding="utf-8", errors="replace")
    energy_rows = parse_energy_rows(text)
    max_err_total = 0.0
    max_err_pp = 0.0
    if energy_rows:
        max_err_total = max(abs(r.get("Error/Total", 0.0)) for r in energy_rows)
        max_err_pp = max(abs(r.get("Error_PP", 0.0)) for r in energy_rows)

    regex_cfg = extract_cfg.get("regex_counts", {})
    regex_counts = extract_regex_counts(text, regex_cfg) if regex_cfg else {}

    return {
        "max_abs_error_over_total": max_err_total,
        "max_abs_error_pp": max_err_pp,
        "energy_sample_count": len(energy_rows),
        "regex_counts": regex_counts,
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
    dry_run: bool,
) -> Tuple[List[RunResult], List[Dict[str, Any]]]:
    scenario_name = scenario["name"]
    workdir = Path(substitute_vars(scenario.get("workdir", "."), variables))
    if not workdir.is_absolute():
        workdir = (base_dir / workdir).resolve()

    scenario_out = out_dir / scenario_name
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
            metrics = {}
        else:
            code = run_shell(cmd, workdir, output_path)
            if code != 0:
                raise RuntimeError(f"Run command failed ({scenario_name}/{run_id}): {cmd}")
            metrics = extract_metrics(output_path, scenario.get("extract", {}))

        run_results.append(RunResult(run_id=run_id, command=cmd, output_path=output_path, metrics=metrics))

    check_results = []
    if not dry_run:
        per_run = {item.run_id: item.metrics for item in run_results}
        for check in scenario.get("checks", []):
            ctype = check["type"]
            if ctype == "max_threshold":
                run_id = check["run_id"]
                metric = check["metric"]
                threshold = float(check["threshold"])
                value = float(per_run[run_id][metric])
                result = check_max_threshold(value, threshold)
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
                expected_ratio = float(check["expected_ratio"])
                tolerance = float(check["relative_tolerance"])
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
                threshold = int(check["threshold"])
                value = int(per_run[run_id]["regex_counts"].get(regex_key, 0))
                result = check_max_threshold(float(value), float(threshold))
                result.update({"scenario": scenario_name, "check": check["name"], "value": value, "threshold": threshold})
                check_results.append(result)
            else:
                raise ValueError(f"Unknown check type: {ctype}")

    return run_results, check_results


def main() -> int:
    parser = argparse.ArgumentParser(description="PeTar validation runner (minimal framework)")
    parser.add_argument("--scenario", action="append", help="Scenario JSON path (can be repeated)")
    parser.add_argument("--scenario-dir", default="test/validation/scenarios", help="Default scenario directory")
    parser.add_argument("--out-dir", default="test/validation/out", help="Output directory")
    parser.add_argument("--var", action="append", default=[], help="Template variables key=value")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running")
    parser.add_argument("--report", default="test/validation/out/report.json", help="Report JSON path")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[2]
    scenario_dir = (repo_root / args.scenario_dir).resolve()
    out_dir = (repo_root / args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    variables = parse_key_value(args.var)
    variables.setdefault("repo_root", str(repo_root))
    variables.setdefault("petar_bin_order2", "petar")
    variables.setdefault("petar_bin_order4", "petar")
    variables.setdefault("petar_bin_switch", "petar")

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