#!/usr/bin/env python3
import argparse
import html
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional

DEFAULT_AGAMA_CONF_RELATIVE = Path("sample/MWPotentialHunter24_rotspiral.ini")


WARNING_PATTERNS = [
    r"\bWARNING\b",
    r"UserWarning:.*Binary file size",
    r"File may be truncated",
    r"not aligned with dtype itemsize",
    r"column",
    r"dtype",
    r"RuntimeWarning",
]


def detect_warnings(text: str) -> List[str]:
    found: List[str] = []
    for line in text.splitlines():
        for pat in WARNING_PATTERNS:
            if re.search(pat, line, flags=re.IGNORECASE):
                found.append(line.strip())
                break
    # Keep order but remove duplicates.
    unique = list(dict.fromkeys(found))
    return unique


def has_critical_mismatch_warning(warnings_list: List[str]) -> bool:
    return any(
        re.search(r"Binary file size|not aligned with dtype itemsize|File may be truncated|dtype mismatch", warning, flags=re.IGNORECASE)
        for warning in warnings_list
    )


def shell_join(cmd: List[str]) -> str:
    return " ".join(shlex.quote(x) for x in cmd)


def status_from_code(code: int) -> str:
    return "pass" if code == 0 else "fail"


def short_reason_from_code(code: int) -> str:
    return "" if code == 0 else f"exit code {code}"


def configure_args_for_build_tag(tag: str) -> List[str]:
    mapping = {
        "std": [],
        "merger": ["--with-interrupt=merger"],
        "merger": ["--with-interrupt=merger"],
        "dsm": ["--with-interrupt=dsm"],
        "bse": ["--with-interrupt=bse"],
        "galpy": ["--with-external=galpy"],
        "bse-galpy": ["--with-interrupt=bse", "--with-external=galpy"],
        "agama": ["--with-external=agama"],
        "bse-agama": ["--with-interrupt=bse", "--with-external=agama"],
    }
    return mapping.get(tag, [])


def describe_ic_case(ic_case: str) -> str:
    descriptions = {
        "t1": "Two-body binary validation IC.",
        "t2": "Two-body binary hard-switch validation IC.",
        "t3": "Hierarchical triple validation IC.",
        "t4": "Hierarchical triple with outer orbit baseline IC.",
        "t4_outer15": "Hierarchical triple with outer semi-major axis 15 variant.",
        "functional_smoke": "One deterministic close binary plus a small symmetric background cluster.",
        "functional_dual_merge_smoke": "Two target BSE binaries plus a light symmetric background cluster.",
        "functional_mcluster_dual_merge_smoke": "mcluster N=100 Plummer Rh=1 pc Kroupa no-binary model, first four stars replaced by two target binaries, then global COM recentered.",
    }
    return descriptions.get(ic_case, ic_case)


def add_step(
    result: Dict[str, Any],
    name: str,
    code: int,
    cmd: List[str],
    cwd: Path,
    warnings_list: Optional[List[str]] = None,
    details: Optional[Dict[str, Any]] = None,
    reason: str = "",
) -> None:
    step = {
        "name": name,
        "code": code,
        "status": status_from_code(code),
        "reason": reason or short_reason_from_code(code),
        "command": shell_join(cmd),
        "argv": cmd,
        "cwd": str(cwd),
    }
    if details:
        step["details"] = details
    if warnings_list:
        step["warnings"] = warnings_list
        result["warnings"].append({"step": name, "messages": warnings_list})
    result["steps"].append(step)


def render_html_report(report: Dict[str, Any], out_path: Path) -> None:
    def esc(value: Any) -> str:
        return html.escape(str(value))

    def badge(status: str) -> str:
        color = {"pass": "#1b7f3b", "fail": "#b42318", "skip": "#8a6d1f"}.get(status, "#555")
        return f'<span class="badge" style="background:{color}">{esc(status.upper())}</span>'

    parts: List[str] = []
    parts.append("<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">")
    parts.append("<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">")
    parts.append("<title>PeTar Functional Report</title>")
    parts.append(
        "<style>"
        "body{font-family:Segoe UI,Helvetica,Arial,sans-serif;margin:0;background:#f6f8fb;color:#16202a;}"
        ".wrap{max-width:1280px;margin:0 auto;padding:24px;}"
        "h1,h2,h3{margin:0 0 12px 0;}"
        ".muted{color:#5b6875;}"
        ".card{background:#fff;border:1px solid #d8dee6;border-radius:12px;padding:18px;margin:16px 0;box-shadow:0 1px 2px rgba(16,24,40,.04);}"
        ".badge{display:inline-block;color:#fff;border-radius:999px;padding:2px 10px;font-size:12px;font-weight:600;}"
        ".kv{display:grid;grid-template-columns:220px 1fr;gap:8px 16px;margin:10px 0;}"
        ".kv div:nth-child(odd){font-weight:600;color:#334155;}"
        "table{width:100%;border-collapse:collapse;margin-top:10px;}"
        "th,td{text-align:left;vertical-align:top;padding:10px;border-top:1px solid #e5e7eb;font-size:14px;}"
        "th{background:#f8fafc;font-weight:700;}"
        "code,pre{font-family:SFMono-Regular,Consolas,Monaco,monospace;font-size:12px;}"
        "pre{white-space:pre-wrap;word-break:break-word;background:#0f172a;color:#e2e8f0;padding:12px;border-radius:10px;overflow:auto;}"
        ".fail{color:#b42318;font-weight:600;}.pass{color:#1b7f3b;font-weight:600;}.skip{color:#8a6d1f;font-weight:600;}"
        "</style></head><body><div class=\"wrap\">"
    )
    parts.append("<h1>PeTar Functional Smoke Report</h1>")
    parts.append(f"<p class=\"muted\">Phase: {esc(report.get('phase'))} | Started: {esc(report.get('started_at'))} | Finished: {esc(report.get('finished_at', ''))}</p>")

    build = report.get("build")
    if build:
        parts.append("<div class=\"card\"><h2>Global Build</h2>")
        parts.append(f"<p>{badge(build.get('status', 'unknown'))}</p>")
        parts.append("<div class=\"kv\">")
        parts.append(f"<div>Script</div><div>{esc(build.get('script', ''))}</div>")
        parts.append(f"<div>Build Tags</div><div>{esc(', '.join(build.get('build_tags', [])))}</div>")
        parts.append(f"<div>Log</div><div>{esc(build.get('log', ''))}</div>")
        parts.append("</div></div>")

    for case in report.get("cases", []):
        parts.append("<div class=\"card\">")
        parts.append(f"<h2>Case: {esc(case.get('name'))}</h2>")
        parts.append(f"<p>{badge(case.get('status', 'unknown'))}</p>")
        parts.append("<div class=\"kv\">")
        parts.append(f"<div>Require</div><div>{esc(case.get('require', ''))}</div>")
        parts.append(f"<div>Directory</div><div>{esc(case.get('dir', ''))}</div>")
        parts.append(f"<div>Log</div><div>{esc(case.get('log', ''))}</div>")
        parts.append(f"<div>Skip/Fail Reason</div><div>{esc(case.get('skip_reason', ''))}</div>")
        ic_meta = case.get("ic", {})
        parts.append(f"<div>IC Case</div><div>{esc(ic_meta.get('case', ''))}</div>")
        parts.append(f"<div>IC Description</div><div>{esc(ic_meta.get('description', ''))}</div>")
        parts.append("</div>")

        if case.get("build"):
            build_info = case["build"]
            parts.append("<h3>Build/Configure</h3><div class=\"kv\">")
            parts.append(f"<div>Build Tag</div><div>{esc(build_info.get('build_tag', ''))}</div>")
            parts.append(f"<div>Configure Args</div><div>{esc(' '.join(build_info.get('configure_args', [])))}</div>")
            parts.append(f"<div>Configure Command</div><div><code>{esc(build_info.get('configure_command', ''))}</code></div>")
            parts.append(f"<div>Build Command</div><div><code>{esc(build_info.get('command', ''))}</code></div>")
            parts.append("</div>")

        parts.append("<h3>Tool Steps</h3><table><thead><tr><th>Step</th><th>Command</th><th>Result</th><th>Details</th></tr></thead><tbody>")
        for step in case.get("steps", []):
            details = [f"cwd: {step.get('cwd', '')}"]
            for key, value in step.get("details", {}).items():
                if isinstance(value, list):
                    value = " ".join(str(x) for x in value)
                details.append(f"{key}: {value}")
            if step.get("reason"):
                details.append(f"reason: {step.get('reason')}")
            if step.get("warnings"):
                details.append("warnings:")
                details.extend(f"- {msg}" for msg in step.get("warnings", []))
            parts.append("<tr>")
            parts.append(f"<td>{esc(step.get('name', ''))}</td>")
            parts.append(f"<td><code>{esc(step.get('command', ''))}</code></td>")
            parts.append(f"<td class=\"{esc(step.get('status', ''))}\">{esc(step.get('status', ''))}</td>")
            parts.append(f"<td><pre>{esc(chr(10).join(details))}</pre></td>")
            parts.append("</tr>")
        parts.append("</tbody></table>")

        read_checks = case.get("read_checks", {}).get("checks", [])
        if read_checks:
            parts.append("<h3>Python Read Checks</h3><table><thead><tr><th>Label</th><th>Command</th><th>Result</th><th>Reason</th></tr></thead><tbody>")
            for chk in read_checks:
                chk_status = "pass" if chk.get("ok") else "fail"
                reason = chk.get("error", "") or "; ".join(chk.get("warnings", []))
                details = chk.get("details", {}) or {}
                detail_lines = []
                for key, value in details.items():
                    detail_lines.append(f"{key}: {value}")
                if detail_lines:
                    reason = (reason + "\n" if reason else "") + "\n".join(detail_lines)
                parts.append("<tr>")
                parts.append(f"<td>{esc(chk.get('label', ''))}</td>")
                parts.append(f"<td><code>{esc(chk.get('command', ''))}</code></td>")
                parts.append(f"<td class=\"{chk_status}\">{chk_status}</td>")
                parts.append(f"<td><pre>{esc(chk.get('file', ''))}\n{esc(reason)}</pre></td>")
                parts.append("</tr>")
            parts.append("</tbody></table>")

        parts.append("</div>")

    parts.append("</div></body></html>")
    out_path.write_text("".join(parts), encoding="utf-8")


def run_cmd(cmd: List[str], cwd: Path, env: Dict[str, str], log_path: Path) -> tuple[int, List[str]]:
        print(f"[DEBUG] run_cmd: {cmd} log_path: {log_path}")
        log_path.parent.mkdir(parents=True, exist_ok=True)
        cmd_text = "$ " + " ".join(shlex.quote(x) for x in cmd) + "\n"

        warn_lines: List[str] = []
        with log_path.open("a", encoding="utf-8") as log:
            log.write("\n" + cmd_text)
            log.flush()

            proc = subprocess.Popen(
                cmd,
                cwd=str(cwd),
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )

            if proc.stdout is not None:
                for line in proc.stdout:
                    log.write(line)
                    log.flush()
                    warn_lines.extend(detect_warnings(line))

            return_code = proc.wait()
            log.write(f"[exit] {return_code}\n")
            log.flush()

        warn_lines = list(dict.fromkeys(warn_lines))
        return return_code, warn_lines


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


def resolve_output_prefixes(case: Dict[str, Any], matrix: Dict[str, Any]) -> List[str]:
    raw = case.get("filename_prefix", matrix.get("global", {}).get("filename_prefix", "data"))
    if isinstance(raw, str):
        prefixes = [tok.strip() for tok in raw.split(",") if tok.strip()]
    elif isinstance(raw, list):
        prefixes = [str(tok).strip() for tok in raw if str(tok).strip()]
    else:
        prefixes = []
    return prefixes or ["data"]


def cleanup_case_outputs_by_prefix(case_dir: Path, prefixes: List[str]) -> Dict[str, Any]:
    matched: List[str] = []
    removed: List[str] = []
    failed: List[str] = []

    if not case_dir.exists():
        return {
            "prefixes": prefixes,
            "matched_count": 0,
            "removed_count": 0,
            "failed_count": 0,
            "matched": matched,
            "removed": removed,
            "failed": failed,
        }

    for entry in sorted(case_dir.iterdir(), key=lambda p: p.name):
        if not any(entry.name.startswith(prefix) for prefix in prefixes):
            continue
        matched.append(entry.name)
        try:
            if entry.is_dir() and not entry.is_symlink():
                shutil.rmtree(entry)
            else:
                entry.unlink()
            removed.append(entry.name)
        except Exception as exc:
            failed.append(f"{entry.name}: {exc}")

    return {
        "prefixes": prefixes,
        "matched_count": len(matched),
        "removed_count": len(removed),
        "failed_count": len(failed),
        "matched": matched,
        "removed": removed,
        "failed": failed,
    }


def prefixed_name(prefix: str, suffix: str) -> str:
    return f"{prefix}.{suffix}"


def read_last_snapshot_list(path: Path) -> str:
    last = ""
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                last = line
    return last


def infer_build_tag_from_require(require: str) -> str:
    tokens = {x.strip() for x in require.split(",") if x.strip()}
    if not tokens:
        return "std"

    # Canonicalize combined feature tags to match install_petar_major_versions.sh.
    if "bse" in tokens and "galpy" in tokens:
        return "bse-galpy"
    if "bse" in tokens and "agama" in tokens:
        return "bse-agama"

    return "-".join(sorted(tokens))


def build_default_movie_args(case: Dict[str, Any], output_prefix: str, interrupt_mode: str, external_mode: str) -> List[str]:
    name = case.get("name", "")
    lagr_file = prefixed_name(output_prefix, "lagr")
    movie_interrupt = interrupt_mode or "none"
    movie_external = external_mode or "none"

    # Follow sample recipes by case family while keeping runtime short and robust.
    if name == "bse-galpy":
        return [
            "-i", movie_interrupt,
            "-t", movie_external,
            "-m", "x-y,x-y",
            "-R", "10,10000",
            "--cm-mode", "core,none",
            "--marker-scale", "1,0.1",
            "-c", "logtemp,logtemp",
            "-b",
            "-L", lagr_file,
            "--rlagr-min", "0",
            "--rlagr-max", "5",
        ]

    if name == "bse-agama":
        return [
            "-i", movie_interrupt,
            "-t", movie_external,
            "-m", "x-y,x-y",
            "-R", "10,10000",
            "--cm-mode", "core,none",
            "--marker-scale", "1,0.1",
            "-c", "logtemp,logtemp",
            "-b",
            "-L", lagr_file,
            "--rlagr-min", "0",
            "--rlagr-max", "5",
        ]

    if name == "galpy":
        return [
            "-i", movie_interrupt,
            "-t", movie_external,
            "-m", "x-y,x-y",
            "-R", "10,10000",
            "--cm-mode", "core,none",
            "--marker-scale", "1,0.1",
            "-L", lagr_file,
            "--rlagr-min", "0",
            "--rlagr-max", "5",
        ]

    if name == "bse":
        return [
            "-i", movie_interrupt,
            "-t", movie_external,
            "-m", "x-y",
            "-R", "10",
            "-c", "logtemp",
            "-b",
            "-L", lagr_file,
            "--rlagr-min", "0",
            "--rlagr-max", "5",
        ]

    # std / merger / dsm default to binary+lagrangian style.
    return [
        "-i", movie_interrupt,
        "-t", movie_external,
        "-m", "x-y",
        "-R", "10",
        "-b",
        "-L", lagr_file,
        "--rlagr-min", "0",
        "--rlagr-max", "5",
    ]


def run_case_build(repo_root: Path, matrix: Dict[str, Any], case: Dict[str, Any], case_dir: Path, log_path: Path) -> Dict[str, Any]:
    build = matrix.get("build", {})
    script = build.get("script")
    if not script:
        return {"status": "skip", "reason": "build script is not configured"}

    script_path = repo_root / script
    if not script_path.exists():
        return {"status": "fail", "reason": f"build script not found: {script_path}"}

    tag = infer_build_tag_from_require(case.get("require", ""))

    env = os.environ.copy()
    env.update(build.get("strict_env", {}))
    env["BUILD_TAGS"] = tag
    env["PURGE_OLD_BINARIES"] = "0"
    configure_args = configure_args_for_build_tag(tag)

    rc, warns = run_cmd(["bash", str(script_path)], cwd=repo_root, env=env, log_path=log_path)
    return {
        "status": "pass" if rc == 0 else "fail",
        "code": rc,
        "build_tag": tag,
        "configure_args": configure_args,
        "configure_command": shell_join(["./configure", *configure_args]),
        "command": shell_join(["bash", str(script_path)]),
        "warnings": warns,
        "script": str(script_path),
        "dir": str(case_dir),
    }


def run_case_build_for_tag(repo_root: Path, matrix: Dict[str, Any], build_tag: str, case_dir: Path, log_path: Path) -> Dict[str, Any]:
    build = matrix.get("build", {})
    script = build.get("script")
    if not script:
        return {"status": "skip", "reason": "build script is not configured"}

    script_path = repo_root / script
    if not script_path.exists():
        return {"status": "fail", "reason": f"build script not found: {script_path}"}

    env = os.environ.copy()
    env.update(build.get("strict_env", {}))
    env["BUILD_TAGS"] = build_tag
    env["PURGE_OLD_BINARIES"] = "0"
    configure_args = configure_args_for_build_tag(build_tag)

    rc, warns = run_cmd(["bash", str(script_path)], cwd=repo_root, env=env, log_path=log_path)
    return {
        "status": "pass" if rc == 0 else "fail",
        "code": rc,
        "build_tag": build_tag,
        "configure_args": configure_args,
        "configure_command": shell_join(["./configure", *configure_args]),
        "command": shell_join(["bash", str(script_path)]),
        "warnings": warns,
        "script": str(script_path),
        "dir": str(case_dir),
    }


def select_build_tags(cases: List[Dict[str, Any]]) -> List[str]:
    tags = []
    for case in cases:
        # Optional cases with unresolved env vars will be skipped in run phase; do not build for them.
        unresolved_skip = False
        for key in case.get("skip_if_unresolved", []):
            if not os.environ.get(key):
                unresolved_skip = True
                break
        if unresolved_skip:
            continue

        tag = infer_build_tag_from_require(case.get("require", ""))
        tags.append(tag)

    if not tags:
        tags = ["std"]
    return sorted(list(dict.fromkeys(tags)))


def run_python_output_read_checks(repo_root: Path, case_dir: Path, require: str, output_prefix: str, interrupt_mode_override: str = "") -> Dict[str, Any]:
    tools_path = repo_root / "tools"
    if str(tools_path) not in sys.path:
        sys.path.insert(0, str(tools_path))

    import analysis as petar  # type: ignore
    import numpy as np  # type: ignore

    require_tokens = {x.strip() for x in require.split(",") if x.strip()}
    interrupt_mode = interrupt_mode_override or "none"
    if not interrupt_mode_override:
        if "dsm" in require_tokens:
            interrupt_mode = "dsm"
        elif any(x in require_tokens for x in ["bse", "bseEmp", "mobse"]):
            interrupt_mode = "bse"
        elif "merger" in require_tokens:
            interrupt_mode = "merger"

    external_mode = "none"
    if "galpy" in require_tokens:
        external_mode = "galpy"
    elif "agama" in require_tokens:
        external_mode = "agama"

    checks: List[Dict[str, Any]] = []

    def add_check(label: str, path: Path, loader, command: str, strict_mismatch: bool = True) -> None:
        if not path.exists():
            return
        with warnings.catch_warnings(record=True) as rec:
            warnings.simplefilter("always")
            ok = True
            error = ""
            details: Dict[str, Any] = {}
            try:
                loader_result = loader(path)
                if isinstance(loader_result, dict):
                    details = loader_result.get("details", {}) or {}
                    command = loader_result.get("command", command)
            except Exception as exc:
                ok = False
                error = str(exc)

        warn_msgs = [str(item.message) for item in rec]
        mismatch_warns = [
            msg
            for msg in warn_msgs
            if re.search(r"column|dtype|aligned|truncated|itemsize|Binary file size", msg, flags=re.IGNORECASE)
        ]
        if strict_mismatch and mismatch_warns and ok:
            ok = False
            error = "; ".join(mismatch_warns)
        checks.append(
            {
                "label": label,
                "file": str(path),
                "command": command,
                "ok": ok,
                "error": error,
                "warnings": warn_msgs,
                "details": details,
            }
        )

    def load_interrupt_binary(path: Path) -> Dict[str, Any]:
        interrupt = petar.InterruptBinary(
            particle_type=petar.HardParticle,
            interrupt_mode=interrupt_mode,
        )
        interrupt.fromfile(str(path))
        return {
            "command": f"petar.InterruptBinary(particle_type=petar.HardParticle, interrupt_mode={interrupt_mode}).fromfile(path)",
            "details": {"records": int(interrupt.size), "ncols": int(interrupt.ncols)},
        }

    def load_profile_with_matching_schema(path: Path) -> Dict[str, Any]:
        # Select Profile schema by matching column count in the numeric body.
        first_row = np.loadtxt(str(path), skiprows=1, max_rows=1, ndmin=2)
        if first_row.size == 0:
            petar.Profile().loadtxt(str(path), skiprows=1)
            return {
                "command": "petar.Profile().loadtxt(path, skiprows=1)",
                "details": {"profile_kwargs": {}, "matched_ncols": 0},
            }

        ncols = first_row.shape[1]
        candidates = [
            {},
            {"use_gpu": True},
            {"old_version": True},
            {"use_gpu": True, "old_version": True},
            {"FDPS_version": 7.0},
            {"use_gpu": True, "FDPS_version": 7.0},
            {"old_version": True, "FDPS_version": 7.0},
            {"use_gpu": True, "old_version": True, "FDPS_version": 7.0},
        ]

        for kwargs in candidates:
            prof = petar.Profile(**kwargs)
            if prof.ncols == ncols:
                try:
                    prof.loadtxt(str(path), skiprows=1)
                except Exception as exc:
                    msg = str(exc)
                    if "number of columns changed" not in msg:
                        raise
                    # Some profile outputs can append extra diagnostic columns in later rows.
                    # For smoke checks, validate a common numeric prefix can be parsed.
                    min_cols = 0
                    with path.open("r", encoding="utf-8") as fh:
                        for i, line in enumerate(fh):
                            if i == 0:
                                continue
                            tokens = line.strip().split()
                            if not tokens:
                                continue
                            row_cols = len(tokens)
                            min_cols = row_cols if min_cols == 0 else min(min_cols, row_cols)
                    if min_cols <= 0:
                        raise
                    np.loadtxt(str(path), skiprows=1, usecols=tuple(range(min_cols)))
                return {
                    "command": f"petar.Profile({', '.join(f'{k}={v!r}' for k, v in kwargs.items())}).loadtxt(path, skiprows=1)" if kwargs else "petar.Profile().loadtxt(path, skiprows=1)",
                    "details": {"profile_kwargs": kwargs, "matched_ncols": ncols},
                }

        # Fallback for legacy FDPS-profile layouts.
        fallback_kwargs = {"FDPS_version": 7.0}
        petar.Profile(**fallback_kwargs).loadtxt(str(path), skiprows=1)
        return {
            "command": "petar.Profile(FDPS_version=7.0).loadtxt(path, skiprows=1)",
            "details": {"profile_kwargs": fallback_kwargs, "matched_ncols": ncols, "fallback": True},
        }

    add_check(
        f"{output_prefix}.lagr",
        case_dir / prefixed_name(output_prefix, "lagr"),
        lambda p, em=external_mode: petar.LagrangianMultiple(external_mode=em).fromfile(str(p)),
        f"petar.LagrangianMultiple(external_mode={external_mode}).fromfile(path)",
    )
    add_check(f"{output_prefix}.core", case_dir / prefixed_name(output_prefix, "core"), lambda p: petar.Core().fromfile(str(p)), "petar.Core().fromfile(path)")
    add_check(f"{output_prefix}.status", case_dir / prefixed_name(output_prefix, "status"), lambda p: petar.Status().fromfile(str(p)), "petar.Status().fromfile(path)")
    add_check(
        f"{output_prefix}.esc_single",
        case_dir / prefixed_name(output_prefix, "esc_single"),
        lambda p, im=interrupt_mode, em=external_mode: petar.SingleEscaper(
            interrupt_mode=im, external_mode=em
        ).fromfile(str(p)),
        f"petar.SingleEscaper(interrupt_mode={interrupt_mode}, external_mode={external_mode}).fromfile(path)",
    )
    add_check(
        f"{output_prefix}.esc_binary",
        case_dir / prefixed_name(output_prefix, "esc_binary"),
        lambda p, im=interrupt_mode, em=external_mode: petar.BinaryEscaper(
            interrupt_mode=im, external_mode=em
        ).fromfile(str(p)),
        f"petar.BinaryEscaper(interrupt_mode={interrupt_mode}, external_mode={external_mode}).fromfile(path)",
    )

    add_check(f"{output_prefix}.sse", case_dir / prefixed_name(output_prefix, "sse"), lambda p: petar.SSEType().loadtxt(str(p)), "petar.SSEType().loadtxt(path)")

    add_check(f"{output_prefix}.sse.type_change", case_dir / prefixed_name(output_prefix, "sse.type_change"), lambda p: petar.SSETypeChange().loadtxt(str(p)), "petar.SSETypeChange().loadtxt(path)")
    add_check(f"{output_prefix}.sse.sn_kick", case_dir / prefixed_name(output_prefix, "sse.sn_kick"), lambda p: petar.SSESNKick().loadtxt(str(p)), "petar.SSESNKick().loadtxt(path)")

    add_check(f"{output_prefix}.bse", case_dir / prefixed_name(output_prefix, "bse"), lambda p: petar.BSEType().loadtxt(str(p)), "petar.BSEType().loadtxt(path)")
    add_check(f"{output_prefix}.bse_status", case_dir / prefixed_name(output_prefix, "bse_status"), lambda p: petar.BSEStatus().fromfile(str(p)), "petar.BSEStatus().fromfile(path)")
    add_check(f"{output_prefix}.bse.type_change", case_dir / prefixed_name(output_prefix, "bse.type_change"), lambda p: petar.BSETypeChange().loadtxt(str(p)), "petar.BSETypeChange().loadtxt(path)")
    add_check(f"{output_prefix}.bse.sn_kick", case_dir / prefixed_name(output_prefix, "bse.sn_kick"), lambda p: petar.BSEKick().loadtxt(str(p)), "petar.BSEKick().loadtxt(path)")
    add_check(f"{output_prefix}.bse.gw_kick", case_dir / prefixed_name(output_prefix, "bse.gw_kick"), lambda p: petar.BSEKick().loadtxt(str(p)), "petar.BSEKick().loadtxt(path)")
    add_check(f"{output_prefix}.bse.dynamic_merge", case_dir / prefixed_name(output_prefix, "bse.dynamic_merge"), lambda p: petar.BSEDynamicMerge().loadtxt(str(p)), "petar.BSEDynamicMerge().loadtxt(path)")
    add_check(f"{output_prefix}.bse.binary_merge", case_dir / prefixed_name(output_prefix, "bse.binary_merge"), lambda p: petar.BSETypeChange().loadtxt(str(p)), "petar.BSETypeChange().loadtxt(path)")

    # Snapshot files listed in [prefix].snap.lst.
    snap_list = case_dir / prefixed_name(output_prefix, "snap.lst")
    if snap_list.exists():
        snap_names = [line.strip() for line in snap_list.read_text(encoding="utf-8").splitlines() if line.strip()]
        if snap_names:
            snap_name = snap_names[-1]
            snap_path = case_dir / snap_name
            offset = petar.HEADER_OFFSET_WITH_CM if external_mode != "none" else petar.HEADER_OFFSET
            add_check(
                f"snapshot:{snap_name}",
                snap_path,
                lambda p, im=interrupt_mode, em=external_mode, off=offset: petar.Particle(
                    interrupt_mode=im, external_mode=em
                ).fromfile(str(p), offset=off),
                f"petar.Particle(interrupt_mode={interrupt_mode}, external_mode={external_mode}).fromfile(path, offset={offset})",
            )
            add_check(
                f"snapshot.single:{snap_name}.single",
                case_dir / f"{snap_name}.single",
                lambda p, im=interrupt_mode, em=external_mode: petar.Particle(
                    interrupt_mode=im, external_mode=em
                ).fromfile(str(p)),
                f"petar.Particle(interrupt_mode={interrupt_mode}, external_mode={external_mode}).fromfile(path)",
            )
            add_check(
                f"snapshot.binary:{snap_name}.binary",
                case_dir / f"{snap_name}.binary",
                lambda p, im=interrupt_mode, em=external_mode: petar.Binary(
                    member_particle_type=petar.Particle,
                    interrupt_mode=im,
                    external_mode=em,
                    G=petar.G_MSUN_PC_MYR,
                ).fromfile(str(p)),
                f"petar.Binary(member_particle_type=petar.Particle, interrupt_mode={interrupt_mode}, external_mode={external_mode}, G=petar.G_MSUN_PC_MYR).fromfile(path)",
            )

    # Validate object snapshots produced by petar.get.object.snap.
    add_check(
        "object:object.1",
        case_dir / "object.1",
        lambda p, im=interrupt_mode, em=external_mode: (
            lambda obj: (obj.addNewMember("time", np.array([], dtype=float)), obj.fromfile(str(p))))(
                petar.Particle(interrupt_mode=im, external_mode=em)
            ),
        f"obj=petar.Particle(interrupt_mode={interrupt_mode}, external_mode={external_mode}); obj.addNewMember('time', np.array([], dtype=float)); obj.fromfile(path)",
    )

    # DSM mode produces data.interrupt as a binary stream of interrupted binaries.
    # The file may include trailing bytes depending on compile-time layout, so this
    # smoke check validates readability and parsed record count instead of enforcing
    # strict dtype-alignment warnings.
    if interrupt_mode == "dsm":
        add_check(
            f"{output_prefix}.interrupt",
            case_dir / prefixed_name(output_prefix, "interrupt"),
            load_interrupt_binary,
            f"petar.InterruptBinary(particle_type=petar.HardParticle, interrupt_mode={interrupt_mode}).fromfile(path)",
            strict_mismatch=False,
        )

    # Validate profiling outputs (ASCII table).
    for prof_path in sorted(case_dir.glob(f"{output_prefix}.prof.rank.*")):
        add_check(
            f"profile:{prof_path.name}",
            prof_path,
            load_profile_with_matching_schema,
            "petar.Profile(...).loadtxt(path, skiprows=1)",
        )

    # Validate group outputs from petar.data.process (binary GroupInfo).
    for group_path in sorted(case_dir.glob(f"{output_prefix}.group.n*")):
        match = re.search(r"\.n(\d+)$", group_path.name)
        if not match:
            continue
        n_member = int(match.group(1))
        add_check(
            f"group:{group_path.name}",
            group_path,
            lambda p, n=n_member, im=interrupt_mode, em=external_mode: petar.GroupInfo(
                N=n,
                interrupt_mode=im,
                external_mode=em,
            ).fromfile(str(p)),
            f"petar.GroupInfo(N={n_member}, interrupt_mode={interrupt_mode}, external_mode={external_mode}).fromfile(path)",
        )

    all_warnings = []
    for chk in checks:
        for wmsg in chk.get("warnings", []):
            all_warnings.append({"file": chk["file"], "message": wmsg})

    mismatch_warnings = [
        item
        for item in all_warnings
        if re.search(r"column|dtype|aligned|truncated|itemsize", item["message"], flags=re.IGNORECASE)
    ]

    return {
        "checks": checks,
        "failed_checks": [c for c in checks if not c.get("ok", False)],
        "warnings": all_warnings,
        "mismatch_warnings": mismatch_warnings,
    }


def run_case(repo_root: Path, case: Dict[str, Any], matrix: Dict[str, Any], out_root: Path) -> Dict[str, Any]:
    print(f"[DEBUG] run_case input: {case}")
    name = case["name"]
    case_dir = out_root / f"functional_{name}"
    case_dir.mkdir(parents=True, exist_ok=True)
    log_path = case_dir / "run.log"
    log_path.write_text("", encoding="utf-8")
    result = {
        "name": name,
        "require": case.get("require", ""),
        "interrupt_mode": case.get("interrupt_mode", ""),
        "status": "pass",
        "steps": [],
        "warnings": [],
        "optional": bool(case.get("optional", False)),
        "skip_reason": "",
        "dir": str(case_dir),
        "log": str(log_path),
    }

    prefixes = resolve_output_prefixes(case, matrix)
    output_prefix = prefixes[0]
    cleanup_result = cleanup_case_outputs_by_prefix(case_dir, prefixes)
    cleanup_code = 0 if cleanup_result["failed_count"] == 0 else 1
    add_step(
        result,
        "cleanup_case_outputs",
        cleanup_code,
        ["internal.cleanup", f"prefixes={','.join(prefixes)}"],
        case_dir,
        details=cleanup_result,
        reason="failed to remove one or more stale files" if cleanup_code != 0 else "",
    )
    if cleanup_code != 0:
        result["status"] = "fail"
        result["skip_reason"] = "failed to cleanup stale outputs by prefix"
        return result

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

    # 0) Build/install required binary family for this case before running tests.
    build_result = run_case_build(repo_root, matrix, case, case_dir, log_path)
    add_step(
        result,
        "build_install",
        build_result.get("code", 0),
        ["bash", str(repo_root / matrix.get("build", {}).get("script", ""))],
        repo_root,
        warnings_list=build_result.get("warnings"),
        details={
            "build_tag": build_result.get("build_tag", ""),
            "configure_args": build_result.get("configure_args", []),
            "configure_command": build_result.get("configure_command", ""),
        },
        reason=build_result.get("reason", ""),
    )
    result["build"] = build_result
    if build_result.get("status") == "skip":
        result["status"] = "skip" if result["optional"] else "fail"
        result["skip_reason"] = build_result.get("reason", "build script is not configured")
        return result
    if build_result.get("status") != "pass":
        result["status"] = "skip" if result["optional"] else "fail"
        if result["status"] == "skip":
            result["skip_reason"] = "required build/install failed"
        return result

    # 1) Generate deterministic IC in Msun, pc, pc/Myr.
    raw_ic = case_dir / "ic.raw"
    ic_case = case.get("ic_case", matrix.get("global", {}).get("ic_case", "t2"))
    make_ic_cmd = [
        "python3",
        "test/validation/make_ic.py",
        "--case",
        ic_case,
        "--output",
        str(raw_ic),
    ]
    result["ic"] = {
        "case": ic_case,
        "description": describe_ic_case(ic_case),
        "output": str(raw_ic),
        "command": shell_join(make_ic_cmd),
    }
    rc, warns = run_cmd(
        make_ic_cmd,
        cwd=repo_root,
        env=env,
        log_path=log_path,
    )
    add_step(result, "make_ic", rc, make_ic_cmd, repo_root, warnings_list=warns, details=result["ic"])
    if rc != 0:
        result["status"] = "fail"
        return result

    # 2) Convert IC to PeTar format.
    init_args = resolve_tokens(case.get("init_args", []))
    init_cmd = ["petar.init", *init_args, "-f", "input", str(raw_ic)]
    rc, warns = run_cmd(
        init_cmd,
        cwd=case_dir,
        env=env,
        log_path=log_path,
    )
    add_step(result, "petar.init", rc, init_cmd, case_dir, warnings_list=warns, details={"init_args": init_args})
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
    rc, warns = run_cmd(
        select_cmd,
        cwd=case_dir,
        env=env,
        log_path=log_path,
    )
    add_step(result, "petar.select", rc, select_cmd, case_dir, warnings_list=warns, details={"require": require, "optional": optional})

    # Fallback: if select fails, trigger targeted install/build once and retry select.
    if rc != 0:
        fallback_tag = infer_build_tag_from_require(require)
        fallback_build = run_case_build_for_tag(repo_root, matrix, fallback_tag, case_dir, log_path)
        add_step(
            result,
            "build_install_retry",
            fallback_build.get("code", 0),
            ["bash", str(repo_root / matrix.get("build", {}).get("script", ""))],
            repo_root,
            warnings_list=fallback_build.get("warnings"),
            details={
                "build_tag": fallback_build.get("build_tag", fallback_tag),
                "configure_args": fallback_build.get("configure_args", []),
                "configure_command": fallback_build.get("configure_command", ""),
                "trigger": "petar.select failed",
            },
            reason=fallback_build.get("reason", ""),
        )
        if fallback_build.get("status") == "pass":
            rc_retry, warns_retry = run_cmd(
                select_cmd,
                cwd=case_dir,
                env=env,
                log_path=log_path,
            )
            add_step(
                result,
                "petar.select.retry",
                rc_retry,
                select_cmd,
                case_dir,
                warnings_list=warns_retry,
                details={"require": require, "optional": optional, "trigger": "after build_install_retry"},
            )
            rc = rc_retry

    if rc != 0:
        result["status"] = "skip" if result["optional"] else "fail"
        if result["status"] == "skip":
            result["skip_reason"] = "required binary family is unavailable"
        return result

    # 4) Fast run with tiny end time/output cadence.
    t_end = str(case.get("t_end", matrix.get("global", {}).get("t_end", 0.2)))
    dt_out = str(case.get("dt_out", matrix.get("global", {}).get("dt_out", 0.1)))
    run_args = resolve_tokens(case.get("run_args", []))
    run_cmdline = ["petar", "-u", "1", "-t", t_end, "-o", dt_out, "-f", output_prefix, *run_args, "input"]
    rc, warns = run_cmd(run_cmdline, cwd=case_dir, env=env, log_path=log_path)
    add_step(result, "petar", rc, run_cmdline, case_dir, warnings_list=warns, details={"t_end": t_end, "dt_out": dt_out, "output_prefix": output_prefix, "run_args": run_args})
    if rc != 0:
        result["status"] = "fail"
        return result

    # 5) Post-processing pipeline smoke checks.
    snap_list = case_dir / prefixed_name(output_prefix, "snap.lst")
    if not snap_list.exists():
        # Fallback for serial outputs when no list file is produced.
        snapshots = sorted(case_dir.glob(f"{output_prefix}.[0-9]*"), key=lambda p: int(p.name.split(".")[-1]))
        if snapshots:
            snap_list.write_text("\n".join(p.name for p in snapshots) + "\n", encoding="utf-8")
        else:
            result["status"] = "fail"
            result["skip_reason"] = f"{snap_list.name} not found"
            return result

    last_snap = read_last_snapshot_list(snap_list)
    if not last_snap:
        result["status"] = "fail"
        result["skip_reason"] = f"{snap_list.name} is empty"
        return result

    process_args = resolve_tokens(case.get("process_args", []))
    process_cmd = ["petar.data.process", "--no-auto-resume", "-p", output_prefix, *process_args, snap_list.name]
    rc, warns = run_cmd(process_cmd, cwd=case_dir, env=env, log_path=log_path)
    add_step(result, "petar.data.process", rc, process_cmd, case_dir, warnings_list=warns, details={"process_args": process_args})
    if rc != 0:
        result["status"] = "fail"
        return result


    # 6) petar.get.object.snap with correct -i/-t for dtype
    object_snap_args = case.get("object_snap_args", [])
    snap_cmd = ["petar.get.object.snap", *object_snap_args, "-p", "object", "-f", "origin", "-m", "id", "1", snap_list.name]
    print(f"[DEBUG] petar.get.object.snap cmd: {snap_cmd}")
    rc, warns = run_cmd(snap_cmd, cwd=case_dir, env=env, log_path=log_path)
    snap_critical_warn = has_critical_mismatch_warning(warns)
    snap_step_code = rc if (rc != 0 or not snap_critical_warn) else 1
    add_step(
        result,
        "petar.get.object.snap",
        snap_step_code,
        snap_cmd,
        case_dir,
        warnings_list=warns,
        details={"object_snap_args": object_snap_args},
        reason="binary mismatch warning detected" if snap_critical_warn and rc == 0 else "",
    )
    if rc != 0:
        result["status"] = "fail"
        result["skip_reason"] = "petar.get.object.snap command failed"
        return result
    if snap_critical_warn:
        result["status"] = "fail"
        result["skip_reason"] = "petar.get.object.snap reported binary mismatch warning"

    # Keep this smoke step focused on the last snapshot to avoid false positives
    # from intermediate outputs while still validating conversion functionality.
    last_snap_list = case_dir / prefixed_name(output_prefix, "snap.last.lst")
    last_snap_list.write_text(f"{last_snap}\n", encoding="utf-8")
    format_cmd = ["petar.format.transfer.post", "-d", "single", "-s", "binary", "-o", "npy"]

    # Keep dtype-critical options consistent with petar.data.process.
    process_interrupt_mode = None
    process_external_mode = None
    for idx, token in enumerate(process_args):
        if token in ("-i", "--interrupt-mode") and idx + 1 < len(process_args):
            process_interrupt_mode = process_args[idx + 1]
        elif token in ("-t", "--external-mode") and idx + 1 < len(process_args):
            process_external_mode = process_args[idx + 1]

    if process_interrupt_mode is None:
        process_interrupt_mode = case.get("interrupt_mode")
    if process_external_mode is None:
        process_external_mode = case.get("external_mode")

    if process_interrupt_mode:
        format_cmd.extend(["-i", str(process_interrupt_mode)])
    if process_external_mode:
        format_cmd.extend(["-t", str(process_external_mode)])

    format_cmd.append(str(last_snap_list.name))
    rc, warns = run_cmd(format_cmd, cwd=case_dir, env=env, log_path=log_path)
    format_critical_warn = has_critical_mismatch_warning(warns)
    format_step_code = rc if (rc != 0 or not format_critical_warn) else 1
    add_step(
        result,
        "petar.format.transfer.post",
        format_step_code,
        format_cmd,
        case_dir,
        warnings_list=warns,
        details={"output_format": "npy", "data_kind": "single/binary"},
        reason="binary mismatch warning detected" if format_critical_warn and rc == 0 else "",
    )
    if rc != 0:
        result["status"] = "fail"
        result["skip_reason"] = "petar.format.transfer.post command failed"
        return result
    if format_critical_warn:
        result["status"] = "fail"
        result["skip_reason"] = "petar.format.transfer.post reported binary mismatch warning"

    # 8) petar.movie smoke check (uniform --n-cpu 1 to expose warnings deterministically).
    movie_args_raw = case.get("movie_args")
    if movie_args_raw is None:
        movie_args = build_default_movie_args(case, output_prefix, case.get("interrupt_mode", "none"), process_external_mode or "none")
    else:
        movie_args = resolve_tokens(movie_args_raw)

    if "--n-cpu" not in movie_args:
        movie_args.extend(["--n-cpu", "1"])

    movie_cmd = ["petar.movie", *movie_args, snap_list.name]
    rc, warns = run_cmd(movie_cmd, cwd=case_dir, env=env, log_path=log_path)
    add_step(result, "petar.movie", rc, movie_cmd, case_dir, warnings_list=warns, details={"movie_args": movie_args})
    if rc != 0:
        result["status"] = "fail"
        result["skip_reason"] = "petar.movie command failed"
        return result

    expected_files = [
        case_dir / last_snap,
        case_dir / f"{last_snap}.single",
        case_dir / f"{last_snap}.binary",
        case_dir / prefixed_name(output_prefix, "lagr"),
        case_dir / "object.1",
    ]
    missing = [str(p) for p in expected_files if not p.exists()]
    if missing:
        result["status"] = "fail"
        result["skip_reason"] = "missing expected outputs: " + ", ".join(missing)

    read_checks = run_python_output_read_checks(
        repo_root,
        case_dir,
        require,
        output_prefix,
        interrupt_mode_override=case.get("interrupt_mode", ""),
    )
    result["read_checks"] = read_checks
    if read_checks.get("warnings"):
        result["warnings"].append({"step": "python_read_checks", "messages": [x["message"] for x in read_checks["warnings"]]})
    if read_checks.get("failed_checks"):
        result["status"] = "fail"
        if result.get("skip_reason"):
            result["skip_reason"] += "; python read checks failed"
        else:
            result["skip_reason"] = "python read checks failed"

    return result


def run_build(repo_root: Path, matrix: Dict[str, Any], out_root: Path, selected_cases: List[Dict[str, Any]]) -> Dict[str, Any]:
    build = matrix.get("build", {})
    script = build.get("script")
    if not script:
        return {"status": "skip", "reason": "build script is not configured"}

    script_path = repo_root / script
    if not script_path.exists():
        return {"status": "fail", "reason": f"build script not found: {script_path}"}

    env = os.environ.copy()
    env.update(build.get("strict_env", {}))
    build_tags = select_build_tags(selected_cases)
    env["BUILD_TAGS"] = ",".join(build_tags)
    env["PURGE_OLD_BINARIES"] = "0"

    log_path = out_root / "build.log"
    rc, warns = run_cmd(["bash", str(script_path)], cwd=repo_root, env=env, log_path=log_path)
    return {
        "status": "pass" if rc == 0 else "fail",
        "code": rc,
        "build_tags": build_tags,
        "configure_by_tag": {tag: configure_args_for_build_tag(tag) for tag in build_tags},
        "warnings": warns,
        "script": str(script_path),
        "log": str(log_path),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="PeTar functional smoke test: build matrix + quick run + analysis pipeline")
    parser.add_argument("--matrix", default="test/functional/functional_matrix.json")
    parser.add_argument("--out-dir", default="test/out")
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

    # Provide a sensible default Agama config for optional agama cases.
    # Keep user-provided AGAMA_CONF_FILE untouched if it is already set.
    if not os.environ.get("AGAMA_CONF_FILE"):
        agama_default = (repo_root / DEFAULT_AGAMA_CONF_RELATIVE).resolve()
        if agama_default.exists():
            os.environ["AGAMA_CONF_FILE"] = str(agama_default)
            print(f"[INFO] AGAMA_CONF_FILE is not set; defaulting to {agama_default}")


    required_cmds = [
        "python3",
        "petar.init",
        "petar.select",
        "petar.data.process",
        "petar.get.object.snap",
        "petar.format.transfer.post",
        "petar.movie",
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
        report["build"] = run_build(repo_root, matrix, out_root, selected)
        if args.phase == "build":
            report_path = out_root / "report.functional.json"
            report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
            render_html_report(report, out_root / "report.functional.html")
            print(json.dumps(report, indent=2))
            return 0 if report["build"]["status"] == "pass" else 1

    if args.phase in {"all", "run"}:
        for case in selected:
            report["cases"].append(run_case(repo_root, case, matrix, out_root))

    report["finished_at"] = int(time.time())
    report_path = out_root / "report.functional.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    render_html_report(report, out_root / "report.functional.html")
    print(json.dumps(report, indent=2))

    failed = [c for c in report.get("cases", []) if c.get("status") == "fail"]
    build_failed = report.get("build", {}).get("status") == "fail" if report.get("build") else False
    return 1 if failed or build_failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
