#!/usr/bin/env python3
import argparse
import json
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List


def run_cmd(cmd: List[str], cwd: Path, env: Dict[str, str], log_path: Path) -> int:
        print(f"[DEBUG] run_cmd: {cmd} log_path: {log_path}")
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("a", encoding="utf-8") as log:
            log.write("\n$ " + " ".join(shlex.quote(x) for x in cmd) + "\n")
            log.flush()
            proc = subprocess.run(cmd, cwd=str(cwd), env=env, stdout=log, stderr=subprocess.STDOUT, check=False)
            log.write(f"[exit] {proc.returncode}\n")
        # After each step, scan for critical warnings
        with log_path.open("r", encoding="utf-8") as log:
            log_text = log.read()
            # Add more patterns as needed
            critical_patterns = [
                "UserWarning: Binary file size",
                "File may be truncated",
                "is not aligned with dtype itemsize",
                "RuntimeWarning",
                "Traceback (most recent call last):"
            ]
            for pat in critical_patterns:
                if pat in log_text:
                    # Use 100 as a special code for warning-triggered fail
                    return 100
        return proc.returncode


def load_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def apply_vars(text: str) -> str:
    out = text
    for key, value in os.environ.items():
        out = out.replace("{" + key + "}", value)
    return out


def resolve_tokens(tokens: List[str]) -> List[str]:
    return [apply_vars(tok) for tok in tokens]


def unresolved_vars(tokens: List[str]) -> List[str]:
    unresolved = []
    for tok in tokens:
        if "{" in tok and "}" in tok:
            unresolved.append(tok)
    return unresolved


def ensure_commands(commands: List[str]) -> List[str]:
    missing = []
    for cmd in commands:
        if subprocess.run(["bash", "-lc", f"command -v {shlex.quote(cmd)} >/dev/null 2>&1"], check=False).returncode != 0:
            missing.append(cmd)
    return missing


def read_first_snapshot_list(path: Path) -> str:
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                return line
    return ""


def run_case(repo_root: Path, case: Dict[str, Any], matrix: Dict[str, Any], out_root: Path) -> Dict[str, Any]:
    print(f"[DEBUG] run_case input: {case}")
    name = case["name"]
    case_dir = out_root / f"functional_{name}"
    case_dir.mkdir(parents=True, exist_ok=True)
    log_path = case_dir / "run.log"
    result = {
        "name": name,
        "status": "pass",
        "steps": [],
        "optional": bool(case.get("optional", False)),
        "skip_reason": "",
        "dir": str(case_dir),
        "log": str(log_path),
    }

    all_tokens = []
    for key in ["init_args", "run_args", "process_args"]:
        all_tokens.extend(case.get(key, []))

    skip_keys = case.get("skip_if_unresolved", [])
    if skip_keys:
        for key in skip_keys:
            if not os.environ.get(key):
                result["status"] = "skip"
                result["skip_reason"] = f"environment variable {key} is not set"
                return result

    unresolved = unresolved_vars(resolve_tokens(all_tokens))
    if unresolved:
        result["status"] = "skip" if result["optional"] else "fail"
        result["skip_reason"] = "unresolved tokens: " + ", ".join(unresolved)
        return result

    env = os.environ.copy()
    env.update(matrix.get("global", {}).get("run_env", {}))

    # 1) Generate deterministic IC in Msun, pc, pc/Myr.
    raw_ic = case_dir / "ic.raw"
    ic_case = case.get("ic_case", matrix.get("global", {}).get("ic_case", "t2"))
    rc = run_cmd(
        [
            "python3",
            "test/validation/make_ic.py",
            "--case",
            ic_case,
            "--output",
            str(raw_ic),
        ],
        cwd=repo_root,
        env=env,
        log_path=log_path,
    )
    result["steps"].append({"name": "make_ic", "code": rc})
    if rc != 0:
        if rc == 100:
            result["status"] = "fail"
            result["skip_reason"] = "Critical warning detected in log (see run.log)"
        else:
            result["status"] = "fail"
        return result

    # 2) Convert IC to PeTar format.
    init_args = resolve_tokens(case.get("init_args", []))
    rc = run_cmd(
        ["petar.init", *init_args, "-f", "input", str(raw_ic)],
        cwd=case_dir,
        env=env,
        log_path=log_path,
    )
    result["steps"].append({"name": "petar.init", "code": rc})
    if rc != 0:
        result["status"] = "fail"
        return result

    # 3) Switch to required installed binary family.
    require = case.get("require", "")
    optional = matrix.get("global", {}).get("optional_select", "mpi,omp,avx512,avx2")
    select_cmd = ["petar.select"]
    if require:
        select_cmd.extend(["--require", require])
    if optional:
        select_cmd.extend(["--optional", optional])
    rc = run_cmd(
        select_cmd,
        cwd=case_dir,
        env=env,
        log_path=log_path,
    )
    result["steps"].append({"name": "petar.select", "code": rc})
    if rc != 0:
        result["status"] = "skip" if result["optional"] else "fail"
        if result["status"] == "skip":
            result["skip_reason"] = "required binary family is unavailable"
        return result

    # 4) Fast run with tiny end time/output cadence.
    t_end = str(case.get("t_end", matrix.get("global", {}).get("t_end", 0.2)))
    dt_out = str(case.get("dt_out", matrix.get("global", {}).get("dt_out", 0.1)))
    run_args = resolve_tokens(case.get("run_args", []))
    run_cmdline = ["petar", "-u", "1", "-t", t_end, "-o", dt_out, *run_args, "input"]
    rc = run_cmd(run_cmdline, cwd=case_dir, env=env, log_path=log_path)
    result["steps"].append({"name": "petar", "code": rc})
    if rc != 0:
        result["status"] = "fail"
        return result

    # 5) Post-processing pipeline smoke checks.
    rc = run_cmd(["petar.data.gether", "data"], cwd=case_dir, env=env, log_path=log_path)
    result["steps"].append({"name": "petar.data.gether", "code": rc})
    if rc != 0:
        result["status"] = "fail"
        return result

    snap_list = case_dir / "data.snap.lst"
    if not snap_list.exists():
        result["status"] = "fail"
        result["skip_reason"] = "data.snap.lst not found"
        return result

    first_snap = read_first_snapshot_list(snap_list)
    if not first_snap:
        result["status"] = "fail"
        result["skip_reason"] = "data.snap.lst is empty"
        return result

    process_args = resolve_tokens(case.get("process_args", []))
    rc = run_cmd(["petar.data.process", *process_args, "data.snap.lst"], cwd=case_dir, env=env, log_path=log_path)
    result["steps"].append({"name": "petar.data.process", "code": rc})
    if rc != 0:
        result["status"] = "fail"
        return result


    # 6) petar.get.object.snap with correct -i/-t for dtype
    object_snap_args = case.get("object_snap_args", [])
    snap_cmd = ["petar.get.object.snap", *object_snap_args, "-p", "object", "-f", "origin", "-m", "id", "1", "data.snap.lst"]
    print(f"[DEBUG] petar.get.object.snap cmd: {snap_cmd}")
    rc = run_cmd(snap_cmd, cwd=case_dir, env=env, log_path=log_path)
    result["steps"].append({"name": "petar.get.object.snap", "code": rc})
    if rc != 0:
        result["status"] = "fail"
        return result

    rc = run_cmd(["petar.format.transfer.post", "-d", "single", "-s", "binary", "-o", "npy", "data.snap.lst"], cwd=case_dir, env=env, log_path=log_path)
    result["steps"].append({"name": "petar.format.transfer.post", "code": rc})
    if rc != 0:
        result["status"] = "fail"
        return result

    expected_files = [
        case_dir / first_snap,
        case_dir / "data.lagr",
        case_dir / "object.1",
    ]
    missing = [str(p) for p in expected_files if not p.exists()]
    if missing:
        result["status"] = "fail"
        result["skip_reason"] = "missing expected outputs: " + ", ".join(missing)

    return result


def run_build(repo_root: Path, matrix: Dict[str, Any], out_root: Path) -> Dict[str, Any]:
    build = matrix.get("build", {})
    script = build.get("script")
    if not script:
        return {"status": "skip", "reason": "build script is not configured"}

    script_path = repo_root / script
    if not script_path.exists():
        return {"status": "fail", "reason": f"build script not found: {script_path}"}

    env = os.environ.copy()
    env.update(build.get("strict_env", {}))

    log_path = out_root / "build.log"
    rc = run_cmd(["bash", str(script_path)], cwd=repo_root, env=env, log_path=log_path)
    return {
        "status": "pass" if rc == 0 else "fail",
        "code": rc,
        "script": str(script_path),
        "log": str(log_path),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="PeTar functional smoke test: build matrix + quick run + analysis pipeline")
    parser.add_argument("--matrix", default="test/functional/functional_matrix.json")
    parser.add_argument("--out-dir", default="test/functional/out")
    parser.add_argument("--phase", choices=["all", "build", "run"], default="all")
    parser.add_argument("--cases", default="all", help="comma-separated case names or 'all'")
    args = parser.parse_args()

    # 其余逻辑保持原缩进

    # 修正 matrix 路径为脚本同目录下的相对路径，确保加载的就是 test/functional/functional_matrix.json
    # 直接用参数路径 resolve，避免路径重复

    matrix = load_json(Path(args.matrix).resolve())
    out_root = Path(args.out_dir).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    repo_root = Path.cwd()


    required_cmds = [
        "python3",
        "petar.init",
        "petar.select",
        "petar.data.gether",
        "petar.data.process",
        "petar.get.object.snap",
        "petar.format.transfer.post",
    ]
    missing_cmds = ensure_commands(required_cmds)
    if missing_cmds:
        print("[ERROR] missing commands:", ", ".join(missing_cmds), file=sys.stderr)
        return 2


    selected = matrix.get("cases", [])
    if args.cases != "all":
        wanted = {x.strip() for x in args.cases.split(",") if x.strip()}
        selected = [c for c in selected if c.get("name") in wanted]


    if not selected and args.phase in {"all", "run"}:
        print("[ERROR] no cases selected", file=sys.stderr)
        return 2

    report: Dict[str, Any] = {
        "started_at": int(time.time()),
        "phase": args.phase,
        "build": None,
        "cases": [],
    }

    if args.phase in {"all", "build"}:
        report["build"] = run_build(repo_root, matrix, out_root)
        if args.phase == "build":
            report_path = out_root / "report.functional.json"
            report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
            print(json.dumps(report, indent=2))
            return 0 if report["build"]["status"] == "pass" else 1

    if args.phase in {"all", "run"}:
        for case in selected:
            report["cases"].append(run_case(repo_root, case, matrix, out_root))

    report["finished_at"] = int(time.time())
    report_path = out_root / "report.functional.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))

    failed = [c for c in report.get("cases", []) if c.get("status") == "fail"]
    build_failed = report.get("build", {}).get("status") == "fail" if report.get("build") else False
    return 1 if failed or build_failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
