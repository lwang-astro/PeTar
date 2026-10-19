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
import threading
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
    # mismatch-shaped only — bare "column"/"dtype" matched informational tool
    # echoes (e.g. petar.init "Add stellar evolution columns: no")
    r"number of columns",
    r"column.{0,40}mismatch|mismatch.{0,40}column",
    r"dtype.{0,40}mismatch|mismatch.{0,40}dtype",
    r"RuntimeWarning",
]


def detect_warnings(text: str) -> List[str]:
    found: List[str] = []
    for line in text.splitlines():
        # petar.data.clear prints one informational banner per suffix
        # ("<prefix>.<suffix>; time column: N", a configured fallback index,
        # mostly placeholders) — not a warning.
        if re.search(r"; time column: \d+$", line):
            continue
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
        "dsm": ["--with-interrupt=dsm"],
        "bse": ["--with-interrupt=bse"],
        "galpy": ["--with-external=galpy"],
        "bse-galpy": ["--with-interrupt=bse", "--with-external=galpy"],
        "agama": ["--with-external=agama"],
        "bse-agama": ["--with-interrupt=bse", "--with-external=agama"],
        "sevn": ["--with-interrupt=sevn"],
        "gasdrag": ["--with-external-hard=gasdrag"],
        "bse-gasdrag": ["--with-interrupt=bse", "--with-external-hard=gasdrag"],
        "bse-galpy-gasdrag": ["--with-interrupt=bse", "--with-external=galpy", "--with-external-hard=gasdrag"],
        "dsm-gasdrag": ["--with-interrupt=dsm", "--with-external-hard=gasdrag"],
        "dsm-galpy-gasdrag": ["--with-interrupt=dsm", "--with-external=galpy", "--with-external-hard=gasdrag"],
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
        "functional_mcluster_smoke": "mcluster N=100 Plummer Rh=1 pc Kroupa model with 10 primordial hard binaries (members written first as adjacent pairs); for modes without a stellar-evolution merge channel (dual-merge contact binaries would stall them).",
        "functional_mcluster_pair_smoke": "Binary-free mcluster N=100 cluster plus one inbound hyperbolic 10+10 Msun collision pair written as the FIRST member pair (petar -b 1); used by merger and dsm — the interrupt paths hang at t=0 when the pair is not pre-marked as a member pair.",
        "functional_binary_cluster": "mcluster N=200 Plummer Rh=0.5 pc Kroupa 95% primordial-binary model (seed 11), members written first as adjacent pairs.",
        "functional_sevn_cluster": "mcluster N=60 Plummer Rh=1 pc Kroupa single-star model (seed 7), masses restricted to the SEVN gold-table grid 0.7-80 Msun.",
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
        "tr.step-details td{background:#fbfcfe;border-top:none;padding-top:2px;padding-bottom:10px;}"
        "tr.step-details pre{margin:0;background:#f1f5f9;color:#334155;}"
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

    preflight = report.get("preflight")
    if preflight:
        parts.append("<div class=\"card\"><h2>Preflight (binary families via petar.select)</h2>")
        parts.append("<table><thead><tr><th>Tag</th><th>Require</th><th>Status</th><th>Selected</th><th>Version</th><th>Expected</th></tr></thead><tbody>")
        for tag, entry in preflight.items():
            parts.append("<tr>")
            parts.append(f"<td>{esc(tag)}</td><td>{esc(entry.get('require', ''))}</td>")
            parts.append(f"<td class=\"{esc(entry.get('status', ''))}\">{esc(entry.get('status', ''))}</td>")
            parts.append(f"<td><code>{esc(entry.get('selected', ''))}</code></td>")
            parts.append(f"<td>{esc(entry.get('version', ''))}</td><td>{esc(entry.get('expected_version', ''))}</td>")
            parts.append("</tr>")
        parts.append("</tbody></table></div>")

    for case in report.get("cases", []):
        parts.append("<div class=\"card\">")
        case_name_html = esc(case.get('name'))
        per_case_html = out_path.parent / f"report.functional.{case.get('name')}.html"
        if per_case_html.exists() and per_case_html.resolve() != out_path.resolve():
            case_name_html = f'<a href="{esc(per_case_html.name)}">{case_name_html}</a>'
        parts.append(f"<h2>Case: {case_name_html}</h2>")
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

        parts.append("<h3>Tool Steps</h3><table><thead><tr><th>Step</th><th>Command</th><th>Result</th></tr></thead><tbody>")
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
            parts.append("</tr>")
            if details:
                parts.append(f"<tr class=\"step-details\"><td colspan=\"3\"><pre>{esc(chr(10).join(details))}</pre></td></tr>")
        parts.append("</tbody></table>")

        read_checks = case.get("read_checks", {}).get("checks", [])
        if read_checks:
            parts.append("<h3>Python Read Checks</h3><table><thead><tr><th>Label</th><th>Command</th><th>Result</th></tr></thead><tbody>")
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
                parts.append("</tr>")
                parts.append(f"<tr class=\"step-details\"><td colspan=\"3\"><pre>{esc(chk.get('file', ''))}\n{esc(reason)}</pre></td></tr>")
            parts.append("</tbody></table>")

        parts.append("</div>")

    parts.append("</div></body></html>")
    out_path.write_text("".join(parts), encoding="utf-8")


def run_cmd(cmd: List[str], cwd: Path, env: Dict[str, str], log_path: Path, timeout_s: Optional[float] = None) -> tuple[int, List[str]]:
        """Run a command, streaming output to the log; kill it after timeout_s seconds.

        A killed step returns -9 (SIGKILL); the log carries a '[timeout]' marker
        line — zero-output silent hangs otherwise take ~20 min to notice by hand.
        """
        log_path.parent.mkdir(parents=True, exist_ok=True)
        cmd_text = "$ " + " ".join(shlex.quote(x) for x in cmd) + "\n"

        warn_lines: List[str] = []
        timed_out = False
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

            timer = None
            if timeout_s is not None and timeout_s > 0:
                def _kill_on_timeout() -> None:
                    nonlocal timed_out
                    timed_out = True
                    proc.kill()

                timer = threading.Timer(timeout_s, _kill_on_timeout)
                timer.daemon = True
                timer.start()

            try:
                if proc.stdout is not None:
                    for line in proc.stdout:
                        log.write(line)
                        log.flush()
                        warn_lines.extend(detect_warnings(line))

                return_code = proc.wait()
            finally:
                if timer is not None:
                    timer.cancel()

            if timed_out:
                log.write(f"[timeout] killed after {timeout_s}s (exit code -9)\n")
                return_code = -9
            else:
                log.write(f"[exit] {return_code}\n")
            log.flush()

        warn_lines = list(dict.fromkeys(warn_lines))
        return return_code, warn_lines


def _step_timeout(matrix: Dict[str, Any]) -> float:
    """Per-step wall-clock limit (matrix global 'step_timeout_s', default 1800)."""
    try:
        return float(matrix.get("global", {}).get("step_timeout_s", 1800))
    except (TypeError, ValueError):
        return 1800.0


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


def resolve_sevn_tables() -> Optional[str]:
    """Resolve the SEVN table set the same way as sample/star_cluster_plummer_N1k_binaries_sevn.sh.

    Order: SEVN_TABLES env > gold table set next to the --tables default baked
    into the selected petar binary > SEVN_DIR source checkout.
    Returns None when no table set is found (case should skip).
    """
    env_tables = os.environ.get("SEVN_TABLES")
    if env_tables and Path(env_tables).is_dir():
        return env_tables

    try:
        proc = subprocess.run(["petar", "-h"], capture_output=True, text=True, check=False)
    except OSError:
        return None
    default_tables = None
    for line in (proc.stdout + proc.stderr).splitlines():
        match = re.match(r"\s*--tables\s+(\S+)", line)
        if match:
            default_tables = match.group(1)
            break
    if default_tables:
        candidate = Path(default_tables).parent / "SEVNtracks_MIST_AGBrobust_gold"
        if candidate.is_dir():
            return str(candidate)

    sevn_dir = os.environ.get("SEVN_DIR")
    if sevn_dir:
        candidate = Path(sevn_dir) / "tables" / "SEVNtracks_MIST_AGBrobust_gold"
        if candidate.is_dir():
            return str(candidate)
    return None


def expected_source_stamp(repo_root: Path) -> str:
    """Version stamp pair expected from a current build: <PeTar VERSION>_<SDAR VERSION>."""
    petar_version = ""
    sdar_version = ""
    vf = repo_root / "VERSION"
    if vf.exists():
        petar_version = vf.read_text(encoding="utf-8").strip()
    sf = (repo_root / ".." / "SDAR" / "VERSION").resolve()
    if sf.exists():
        sdar_version = sf.read_text(encoding="utf-8").strip()
    return f"{petar_version}_{sdar_version}"


def petar_stamp(solver: str = "petar") -> str:
    """Version stamp printed by '<solver> -h' (empty when unreadable)."""
    try:
        proc = subprocess.run([solver, "-h"], capture_output=True, text=True, check=False)
    except OSError:
        return ""
    for line in (proc.stdout + proc.stderr).splitlines():
        match = re.match(r"\s*Version:\s*(\S+)", line)
        if match:
            return match.group(1)
    return ""


def active_petar_realpath() -> str:
    proc = subprocess.run(["bash", "-lc", "command -v petar"], capture_output=True, text=True, check=False)
    path = proc.stdout.strip().splitlines()[0] if proc.stdout.strip() else ""
    if not path:
        return ""
    return str(Path(path).resolve())


def estimate_group_records(repo_root: Path, group_path: Path) -> Optional[int]:
    """Estimate the record count of a binary group file via the petar package (None when unavailable)."""
    match = re.search(r"\.n(\d+)$", group_path.name)
    if not match:
        return None
    n_member = int(match.group(1))
    try:
        tools_path = repo_root / "tools"
        if str(tools_path) not in sys.path:
            sys.path.insert(0, str(tools_path))
        import numpy as _np  # noqa: WPS433
        import analysis as _petar  # noqa: WPS433

        dt = _np.dtype(_petar.GroupInfo(N=n_member).collectDtype())
    except Exception:
        return None
    if dt.itemsize <= 0:
        return None
    return group_path.stat().st_size // dt.itemsize


def duplicate_group_records(repo_root: Path, group_path: Path) -> Optional[int]:
    """Count duplicated group records by (time, id1, id2) — different groups may
    legitimately share a time, so time alone is not a duplicate signal."""
    import collections

    match = re.search(r"\.n(\d+)$", group_path.name)
    if not match:
        return None
    n_member = int(match.group(1))
    try:
        tools_path = repo_root / "tools"
        if str(tools_path) not in sys.path:
            sys.path.insert(0, str(tools_path))
        import numpy as _np  # noqa: WPS433
        import analysis as _petar  # noqa: WPS433

        gi = _petar.GroupInfo(N=n_member)
        gi.fromfile(str(group_path))
        keys = list(
            zip(
                map(float, _np.asarray(gi.time)),
                map(int, _np.asarray(gi.bin0.p1.id)),
                map(int, _np.asarray(gi.bin0.p2.id)),
            )
        )
    except Exception:
        return None
    counts = collections.Counter(keys)
    return sum(c - 1 for c in counts.values() if c > 1)


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
    to = _step_timeout(matrix)

    rc, warns = run_cmd(["bash", str(script_path)], cwd=repo_root, env=env, log_path=log_path, timeout_s=to)
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
    to = _step_timeout(matrix)
    to = _step_timeout(matrix)

    rc, warns = run_cmd(["bash", str(script_path)], cwd=repo_root, env=env, log_path=log_path, timeout_s=to)
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


def prune_stale_case_reports(out_root: Path, case_names: List[str]) -> List[str]:
    """Remove per-case report HTMLs for cases no longer in the matrix."""
    keep = {"report.functional.html"} | {f"report.functional.{n}.html" for n in case_names}
    removed = []
    for path in out_root.glob("report.functional.*.html"):
        if path.name not in keep:
            path.unlink()
            removed.append(path.name)
    return removed


def preflight_families(repo_root: Path, matrix: Dict[str, Any], selected: List[Dict[str, Any]], out_root: Path) -> Dict[str, Dict[str, Any]]:
    """Resolve every binary family before any case runs.

    Mirrors the validation-suite rule: petar.select per family; when the
    selection fails or the selected binary's version stamp does not match the
    repository (PeTar VERSION + SDAR VERSION pair), reinstall the family via
    the build script, then re-select. Prevents the stale-binary failure class
    (an old installed family carrying abort paths already removed from
    current source). The per-case petar.select step remains the authority for
    availability; this preflight guarantees freshness up front.
    """
    log_path = out_root / "preflight.log"
    log_path.write_text("", encoding="utf-8")
    optional = matrix.get("global", {}).get("optional_select", "mpi,omp,avx512,avx2")
    expected = expected_source_stamp(repo_root)
    to = _step_timeout(matrix)

    def select(require: str) -> int:
        cmd = ["petar.select"]
        if require:
            cmd.extend(["--require", require])
        if optional:
            cmd.extend(["--optional", optional])
        rc, _ = run_cmd(cmd, cwd=repo_root, env=os.environ.copy(), log_path=log_path, timeout_s=to)
        return rc

    tags: List[str] = []
    require_by_tag: Dict[str, str] = {}
    for case in selected:
        if any(not os.environ.get(key) for key in case.get("skip_if_unresolved", [])):
            continue
        tag = infer_build_tag_from_require(case.get("require", ""))
        if tag not in require_by_tag:
            require_by_tag[tag] = case.get("require", "")
            tags.append(tag)

    results: Dict[str, Dict[str, Any]] = {}
    for tag in tags:
        require = require_by_tag[tag]
        entry: Dict[str, Any] = {"tag": tag, "require": require, "expected_version": expected}
        rc = select(require)
        if rc == 0:
            entry["selected"] = active_petar_realpath()
            entry["version"] = petar_stamp()
            if entry["version"] == expected:
                entry["status"] = "current"
                results[tag] = entry
                print(f"[preflight] {tag}: current ({entry['version']})")
                continue
            entry["status"] = "stale"
        else:
            entry["status"] = "missing"
            entry["selected"] = ""
            entry["version"] = ""

        print(f"[preflight] {tag}: {entry['status']} -> reinstalling")
        rebuild = run_case_build_for_tag(repo_root, matrix, tag, out_root, log_path)
        entry["rebuild_status"] = rebuild.get("status", "")
        rc2 = select(require)
        if rc2 == 0:
            entry["selected"] = active_petar_realpath()
            entry["version"] = petar_stamp()
            entry["status"] = "rebuilt" if entry["version"] == expected else "stale-after-rebuild"
        else:
            entry["selected"] = ""
            entry["version"] = ""
            entry["status"] = "unavailable"
        print(f"[preflight] {tag}: {entry['status']} ({entry['version']})")
        results[tag] = entry
    return results


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
        elif "sevn" in require_tokens:
            interrupt_mode = "sevn"
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
            loaded = None
            try:
                loader_result = loader(path)
                if isinstance(loader_result, dict):
                    details = loader_result.get("details", {}) or {}
                    command = loader_result.get("command", command)
                    loaded = loader_result.get("_object")
                else:
                    loaded = loader_result
            except Exception as exc:
                ok = False
                error = str(exc)

        # Record count of the loaded data (DictNpArrayMix .size = number of rows).
        if loaded is not None and "records" not in details:
            size = getattr(loaded, "size", None)
            if size is not None:
                try:
                    details["records"] = int(size)
                except (TypeError, ValueError):
                    pass

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
        # data.interrupt is an ASCII single-line-per-record stream
        # (printColumnAscii), NOT binary — fromfile() misparses the text bytes
        # and the old strict_mismatch=False tolerance turned that into a bogus
        # pass (2140 B file, 744 B dtype, 652 B "trailing").
        interrupt = petar.InterruptBinary(
            particle_type=petar.HardParticle,
            interrupt_mode=interrupt_mode,
        )
        interrupt.loadtxt(str(path))
        return {
            "command": f"petar.InterruptBinary(particle_type=petar.HardParticle, interrupt_mode={interrupt_mode}).loadtxt(path)",
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
                    "details": {"profile_kwargs": kwargs, "matched_ncols": ncols, "records": int(prof.size)},
                }

        # Fallback for legacy FDPS-profile layouts.
        fallback_kwargs = {"FDPS_version": 7.0}
        prof_fb = petar.Profile(**fallback_kwargs)
        prof_fb.loadtxt(str(path), skiprows=1)
        return {
            "command": "petar.Profile(FDPS_version=7.0).loadtxt(path, skiprows=1)",
            "details": {"profile_kwargs": fallback_kwargs, "matched_ncols": ncols, "fallback": True, "records": int(prof_fb.size)},
        }

    # Stellar-event record loaders switch dtype (SEVN columns) via interrupt_mode.
    def se_aware(cls, im: str):
        def load(path: Path):
            obj = cls(interrupt_mode=im)
            obj.loadtxt(str(path))
            return obj

        return load

    # Load and return the object so add_check can report its record count.
    def obj_loader(ctor, method: str = "fromfile", **kwargs):
        def load(path: Path):
            obj = ctor(**kwargs)
            getattr(obj, method)(str(path))
            return obj

        return load

    add_check(
        f"{output_prefix}.lagr",
        case_dir / prefixed_name(output_prefix, "lagr"),
        obj_loader(petar.LagrangianMultiple, external_mode=external_mode),
        f"petar.LagrangianMultiple(external_mode={external_mode}).fromfile(path)",
    )
    add_check(f"{output_prefix}.core", case_dir / prefixed_name(output_prefix, "core"), obj_loader(petar.Core), "petar.Core().fromfile(path)")
    add_check(f"{output_prefix}.status", case_dir / prefixed_name(output_prefix, "status"), obj_loader(petar.Status), "petar.Status().fromfile(path)")
    add_check(
        f"{output_prefix}.esc_single",
        case_dir / prefixed_name(output_prefix, "esc_single"),
        obj_loader(petar.SingleEscaper, interrupt_mode=interrupt_mode, external_mode=external_mode),
        f"petar.SingleEscaper(interrupt_mode={interrupt_mode}, external_mode={external_mode}).fromfile(path)",
    )
    add_check(
        f"{output_prefix}.esc_binary",
        case_dir / prefixed_name(output_prefix, "esc_binary"),
        obj_loader(petar.BinaryEscaper, interrupt_mode=interrupt_mode, external_mode=external_mode),
        f"petar.BinaryEscaper(interrupt_mode={interrupt_mode}, external_mode={external_mode}).fromfile(path)",
    )

    add_check(f"{output_prefix}.sse", case_dir / prefixed_name(output_prefix, "sse"), obj_loader(petar.SSEType, "loadtxt"), "petar.SSEType().loadtxt(path)")

    add_check(
        f"{output_prefix}.sse.type_change",
        case_dir / prefixed_name(output_prefix, "sse.type_change"),
        se_aware(petar.SSETypeChange, interrupt_mode),
        f"petar.SSETypeChange(interrupt_mode={interrupt_mode}).loadtxt(path)",
    )
    add_check(
        f"{output_prefix}.sse.sn_kick",
        case_dir / prefixed_name(output_prefix, "sse.sn_kick"),
        se_aware(petar.SSESNKick, interrupt_mode),
        f"petar.SSESNKick(interrupt_mode={interrupt_mode}).loadtxt(path)",
    )

    add_check(f"{output_prefix}.bse", case_dir / prefixed_name(output_prefix, "bse"), obj_loader(petar.BSEType, "loadtxt"), "petar.BSEType().loadtxt(path)")
    add_check(f"{output_prefix}.bse_status", case_dir / prefixed_name(output_prefix, "bse_status"), obj_loader(petar.BSEStatus), "petar.BSEStatus().fromfile(path)")
    add_check(
        f"{output_prefix}.bse.type_change",
        case_dir / prefixed_name(output_prefix, "bse.type_change"),
        se_aware(petar.BSETypeChange, interrupt_mode),
        f"petar.BSETypeChange(interrupt_mode={interrupt_mode}).loadtxt(path)",
    )
    add_check(
        f"{output_prefix}.bse.sn_kick",
        case_dir / prefixed_name(output_prefix, "bse.sn_kick"),
        se_aware(petar.BSEKick, interrupt_mode),
        f"petar.BSEKick(interrupt_mode={interrupt_mode}).loadtxt(path)",
    )
    add_check(
        f"{output_prefix}.bse.gw_kick",
        case_dir / prefixed_name(output_prefix, "bse.gw_kick"),
        se_aware(petar.BSEKick, interrupt_mode),
        f"petar.BSEKick(interrupt_mode={interrupt_mode}).loadtxt(path)",
    )
    add_check(
        f"{output_prefix}.bse.dynamic_merge",
        case_dir / prefixed_name(output_prefix, "bse.dynamic_merge"),
        se_aware(petar.BSEDynamicMerge, interrupt_mode),
        f"petar.BSEDynamicMerge(interrupt_mode={interrupt_mode}).loadtxt(path)",
    )
    add_check(
        f"{output_prefix}.bse.binary_merge",
        case_dir / prefixed_name(output_prefix, "bse.binary_merge"),
        se_aware(petar.BSETypeChange, interrupt_mode),
        f"petar.BSETypeChange(interrupt_mode={interrupt_mode}).loadtxt(path)",
    )

    # Snapshot files listed in [prefix].snap.lst.
    snap_list = case_dir / prefixed_name(output_prefix, "snap.lst")
    if snap_list.exists():
        snap_names = [line.strip() for line in snap_list.read_text(encoding="utf-8").splitlines() if line.strip()]
        if snap_names:
            snap_name = snap_names[-1]
            snap_path = case_dir / snap_name
            offset = petar.HEADER_OFFSET_WITH_CM if external_mode != "none" else petar.HEADER_OFFSET

            def snap_particle(path: Path, off=offset, im=interrupt_mode, em=external_mode):
                obj = petar.Particle(interrupt_mode=im, external_mode=em)
                obj.fromfile(str(path), offset=off)
                return obj

            add_check(
                f"snapshot:{snap_name}",
                snap_path,
                snap_particle,
                f"petar.Particle(interrupt_mode={interrupt_mode}, external_mode={external_mode}).fromfile(path, offset={offset})",
            )
            add_check(
                f"snapshot.single:{snap_name}.single",
                case_dir / f"{snap_name}.single",
                obj_loader(petar.Particle, interrupt_mode=interrupt_mode, external_mode=external_mode),
                f"petar.Particle(interrupt_mode={interrupt_mode}, external_mode={external_mode}).fromfile(path)",
            )
            add_check(
                f"snapshot.binary:{snap_name}.binary",
                case_dir / f"{snap_name}.binary",
                obj_loader(
                    petar.Binary,
                    member_particle_type=petar.Particle,
                    interrupt_mode=interrupt_mode,
                    external_mode=external_mode,
                    G=petar.G_MSUN_PC_MYR,
                ),
                f"petar.Binary(member_particle_type=petar.Particle, interrupt_mode={interrupt_mode}, external_mode={external_mode}, G=petar.G_MSUN_PC_MYR).fromfile(path)",
            )
            # Early-interval binary readback: the last snapshot can be binary-free
            # by design (injected compact binaries merge via BSE; std/merger ICs
            # carry none), which would make the binary reader check trivial.
            if len(snap_names) >= 2:
                early_name = snap_names[1]
                add_check(
                    f"snapshot.binary.early:{early_name}.binary",
                    case_dir / f"{early_name}.binary",
                    obj_loader(
                        petar.Binary,
                        member_particle_type=petar.Particle,
                        interrupt_mode=interrupt_mode,
                        external_mode=external_mode,
                        G=petar.G_MSUN_PC_MYR,
                    ),
                    f"petar.Binary(member_particle_type=petar.Particle, interrupt_mode={interrupt_mode}, external_mode={external_mode}, G=petar.G_MSUN_PC_MYR).fromfile(path)",
                )

    # Validate object snapshots produced by petar.get.object.snap
    # (naming: [prefix].[object_typename].[origin|single|binary]).
    object_file = case_dir / "object.1.origin"

    def object_particle(path: Path, im=interrupt_mode, em=external_mode):
        obj = petar.Particle(interrupt_mode=im, external_mode=em)
        obj.addNewMember("time", np.array([], dtype=float))
        obj.fromfile(str(path))
        return obj

    add_check(
        "object:object.1.origin",
        object_file,
        object_particle,
        f"obj=petar.Particle(interrupt_mode={interrupt_mode}, external_mode={external_mode}); obj.addNewMember('time', np.array([], dtype=float)); obj.fromfile(path)",
    )

    # Interrupt modes write data.interrupt as an ASCII stream of interrupted
    # binaries (merger/dsm); strict column matching — a mismatch fails the case.
    if interrupt_mode in ("dsm", "merger"):
        add_check(
            f"{output_prefix}.interrupt",
            case_dir / prefixed_name(output_prefix, "interrupt"),
            load_interrupt_binary,
            f"petar.InterruptBinary(particle_type=petar.HardParticle, interrupt_mode={interrupt_mode}).loadtxt(path)",
        )

    # Validate profiling outputs (ASCII table).
    for prof_path in sorted(case_dir.glob(f"{output_prefix}.prof.rank.*")):
        add_check(
            f"profile:{prof_path.name}",
            prof_path,
            load_profile_with_matching_schema,
            "petar.Profile(...).loadtxt(path, skiprows=1)",
        )

    # Validate group outputs (binary GroupInfo) written by the petar run itself
    # when --write-group-info 2 is used: [prefix].group[.MPI rank].n[N_member]
    # (the rank segment is omitted for single-rank runs).
    for group_path in sorted(case_dir.glob(f"{output_prefix}.group*.n*")):
        match = re.search(r"\.group(?:\.\d+)?\.n(\d+)$", group_path.name)
        if not match:
            continue
        n_member = int(match.group(1))
        add_check(
            f"group:{group_path.name}",
            group_path,
            obj_loader(petar.GroupInfo, N=n_member, interrupt_mode=interrupt_mode, external_mode=external_mode),
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


def run_tool_smokes(result: Dict[str, Any], case: Dict[str, Any], case_dir: Path, output_prefix: str, log_path: Path, env: Dict[str, str], t_end: str, run_args: List[str], timeout_s: Optional[float] = None) -> None:
    """Extra per-case tool coverage (matrix key 'tool_smokes').

    Each listed tool runs against this case's outputs so every tool a case's
    workflow involves gets exercised. Tools that mutate inputs (update.par,
    data.clear, data.gether, format.transfer, find.dt, standalone bse) run in
    a scratch copy; readers run in place. A failed smoke fails the case.
    """
    to = timeout_s if timeout_s is not None and timeout_s > 0 else _step_timeout({})
    names = case.get("tool_smokes", [])
    if not names:
        return

    scratch = case_dir / "tool_smoke_scratch"
    if scratch.exists():
        shutil.rmtree(scratch)
    scratch.mkdir(parents=True)

    def copy_to_scratch(pattern: str) -> int:
        n = 0
        for src in case_dir.glob(pattern):
            if src.is_file():
                shutil.copy2(src, scratch / src.name)
                n += 1
        return n

    for tool in names:
        if tool == "read.par":
            for par in [f"{output_prefix}.par", f"{output_prefix}.par.hard"]:
                if (case_dir / par).exists():
                    cmd = ["petar.read.par", par]
                    rc, warns = run_cmd(cmd, cwd=case_dir, env=env, log_path=log_path, timeout_s=to)
                    add_step(result, f"petar.read.par:{par}", rc, cmd, case_dir, warnings_list=warns)
        elif tool == "update.par":
            copy_to_scratch(f"{output_prefix}.par")
            cmd = ["petar.update.par", "-p", f"{output_prefix}.par"]
            rc, warns = run_cmd(cmd, cwd=scratch, env=env, log_path=log_path, timeout_s=to)
            add_step(result, "petar.update.par", rc, cmd, scratch, warnings_list=warns)
        elif tool == "data.clear":
            copied = copy_to_scratch(f"{output_prefix}.group*.n*") + copy_to_scratch(f"{output_prefix}.status")
            clear_t = str(float(t_end) / 2.0)
            cmd = ["petar.data.clear", "-t", clear_t, output_prefix]
            rc, warns = run_cmd(cmd, cwd=scratch, env=env, log_path=log_path, timeout_s=to)
            add_step(result, "petar.data.clear", rc, cmd, scratch, warnings_list=warns, details={"files_copied": copied, "t": clear_t})
        elif tool == "data.gether":
            copy_to_scratch(f"{output_prefix}.status")
            copy_to_scratch(f"{output_prefix}.snap.lst")
            copy_to_scratch(f"{output_prefix}.group*.n*")
            cmd = ["petar.data.gether", output_prefix]
            rc, warns = run_cmd(cmd, cwd=scratch, env=env, log_path=log_path, timeout_s=to)
            add_step(result, "petar.data.gether", rc, cmd, scratch, warnings_list=warns)
        elif tool == "format.transfer":
            ft_binary = Path(active_petar_realpath())
            ft_tool = ft_binary.parent / (ft_binary.name + ".format.transfer")
            if not ft_tool.exists():
                add_step(result, "petar.format.transfer", 1, ["petar.format.transfer"], case_dir, reason=f"family tool not found: {ft_tool.name}")
            else:
                copy_to_scratch(f"{output_prefix}.snap.lst")
                for line in (case_dir / f"{output_prefix}.snap.lst").read_text(encoding="utf-8").splitlines():
                    if line.strip():
                        copy_to_scratch(line.strip())
                cmd = [str(ft_tool), "-b", f"{output_prefix}.snap.lst"]
                rc, warns = run_cmd(cmd, cwd=scratch, env=env, log_path=log_path, timeout_s=to)
                add_step(result, "petar.format.transfer", rc, cmd, scratch, warnings_list=warns)
        elif tool == "find.dt":
            copy_to_scratch("input")
            core_args = []
            skip_next = False
            for tok in run_args:
                if skip_next:
                    skip_next = False
                    continue
                if tok == "--write-group-info":
                    skip_next = True
                    continue
                core_args.append(tok)
            a_args = "-u 1" + ("".join(f" {t}" for t in core_args))
            cmd = ["petar.find.dt", "-a", a_args, "-i", "1", "input"]
            rc, warns = run_cmd(cmd, cwd=scratch, env=env, log_path=log_path, timeout_s=to)
            add_step(result, "petar.find.dt", rc, cmd, scratch, warnings_list=warns, details={"a_args": a_args})
        elif tool == "bse":
            cmd = ["petar.bse", "-n", "5", "--mmin", "0.8", "--mmax", "20", "-t", "100"]
            rc, warns = run_cmd(cmd, cwd=scratch, env=env, log_path=log_path, timeout_s=to)
            add_step(result, "petar.bse", rc, cmd, scratch, warnings_list=warns)
        else:
            add_step(result, f"tool_smoke:{tool}", 1, ["internal.tool_smoke", tool], case_dir, reason="unknown tool smoke name")

    failed = [s for s in result["steps"] if s["name"].startswith(("petar.read.par", "petar.update.par", "petar.data.clear", "petar.data.gether", "petar.format.transfer", "petar.find.dt", "petar.bse")) and s["status"] != "pass"]
    if failed:
        result["status"] = "fail"
        result["skip_reason"] = "tool smoke failed: " + ", ".join(s["name"] for s in failed)


def run_case(repo_root: Path, case: Dict[str, Any], matrix: Dict[str, Any], out_root: Path) -> Dict[str, Any]:
    name = case["name"]
    case_dir = out_root / f"functional_{name}"
    case_dir.mkdir(parents=True, exist_ok=True)
    log_path = case_dir / "run.log"
    log_path.write_text("", encoding="utf-8")
    to = _step_timeout(matrix)
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

    unresolved = [tok for tok in unresolved_vars(resolve_tokens(all_tokens)) if "{SEVN_TABLES}" not in tok]
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
    select_details = {"require": require, "optional": optional}
    if rc == 0:
        select_details["selected"] = active_petar_realpath()
        select_details["petar_version"] = petar_stamp()
    add_step(result, "petar.select", rc, select_cmd, case_dir, warnings_list=warns, details=select_details)

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
    run_args_raw = case.get("run_args", [])

    if any("{SEVN_TABLES}" in tok for tok in run_args_raw):
        sevn_tables = resolve_sevn_tables()
        if not sevn_tables:
            result["status"] = "skip"
            result["skip_reason"] = "SEVN table set not found (checked SEVN_TABLES, SEVN_DIR, baked --tables default)"
            return result
        os.environ["SEVN_TABLES"] = sevn_tables
        env["SEVN_TABLES"] = sevn_tables

    run_args = resolve_tokens(run_args_raw)
    snap_lst_path = case_dir / prefixed_name(output_prefix, "snap.lst")
    restart_cfg = case.get("restart") or matrix.get("global", {}).get("restart")

    def snap_entries() -> List[str]:
        if not snap_lst_path.exists():
            return []
        return [x.strip() for x in snap_lst_path.read_text(encoding="utf-8").splitlines() if x.strip()]

    # ASCII event-stream time columns, verified against the writers (2026-10-19):
    # sse/sevn/mosse rows = id + StarParameter blocks -> tphys of the first
    #   block = col 12 (bse-interface StarParameter: kw,m0,mt,r,mc,rc,ospin*3,
    #   epoch,tphys,lum);
    # bse/mobse/bseEmp/sevnB rows = type + BinaryEvent blocks -> record[0]
    #   (=min(p1.tphys,p2.tphys)) at col 2;
    # interrupt rows = col 1; prof.rank rows = output counter at col 2
    #   (matches petar.data.clear semantics; exact when dt_out=1);
    # bse.dynamic_merge has a label-led layout -> not parsed (lines kind).
    # Keep in sync with tools/data_clear.sh tindex_* assignments.
    def ascii_rule(name: str) -> Optional[str]:
        if ".sevnB" in name:
            return "col2"
        if ".prof.rank" in name:
            return "col2"
        if ".interrupt" in name:
            return "col1"
        if "dynamic_merge" in name:
            return None
        if ".bse" in name or ".mobse" in name or ".bseEmp" in name:
            return "col2"
        if ".sse" in name or ".sevn" in name or ".mosse" in name:
            return "col12"
        return None

    def ascii_line_time(line: str, rule: str) -> Optional[float]:
        cols = line.split()
        idx = {"col1": 1, "col2": 2, "col12": 12}[rule]
        try:
            return float(cols[idx - 1])
        except (IndexError, ValueError):
            return None

    def stream_stats(t_cut: float) -> Dict[str, Dict[str, Any]]:
        """Per-stream record stats at a time cut.

        Reader-parseable streams (status, group, esc solver stream, lagr/core)
        return exact kept/beyond counts by record time; ASCII event streams
        (sse/bse/mosse/sevn/interrupt/prof) return line totals (their record
        time column is label-dependent — kept by petar.data.clear itself).
        """
        stats: Dict[str, Dict[str, Any]] = {}
        try:
            tools_path = repo_root / "tools"
            if str(tools_path) not in sys.path:
                sys.path.insert(0, str(tools_path))
            import numpy as _np  # noqa: WPS433
            import analysis as _petar  # noqa: WPS433

        except Exception:
            _petar = None
            _np = None

        def add_reader(p: Path, ctor, **kw) -> None:
            if _petar is None:
                stats[p.name] = {"kind": "unreadable"}
                return
            try:
                obj = ctor(**kw)
                obj.fromfile(str(p))
                times = [float(x) for x in _np.asarray(obj.time)]
                stats[p.name] = {
                    "kind": "reader",
                    "kept": sum(1 for x in times if x <= t_cut + 1e-9),
                    "beyond": sum(1 for x in times if x > t_cut + 1e-9),
                }
            except Exception:
                stats[p.name] = {"kind": "unreadable"}

        mode_i = case.get("interrupt_mode", "")
        ext_e = "galpy" if "galpy" in {x.strip() for x in require.split(",") if x.strip()} else ("agama" if "agama" in {x.strip() for x in require.split(",") if x.strip()} else "none")

        for pattern in [f"{output_prefix}.sse*", f"{output_prefix}.bse*", f"{output_prefix}.mosse*", f"{output_prefix}.sevn*", f"{output_prefix}.interrupt*", f"{output_prefix}.prof.rank*"]:
            for p in sorted(case_dir.glob(pattern)):
                if p.is_file():
                    rule = ascii_rule(p.name)
                    if rule is None:
                        stats[p.name] = {"kind": "lines", "total": sum(1 for _ in p.open("r", encoding="utf-8", errors="ignore"))}
                        continue
                    kept = beyond = unknown = 0
                    with p.open("r", encoding="utf-8", errors="ignore") as fh:
                        for line in fh:
                            if not line.strip():
                                continue
                            tv = ascii_line_time(line, rule)
                            if tv is None:
                                unknown += 1
                            elif tv <= t_cut + 1e-9:
                                kept += 1
                            else:
                                beyond += 1
                    stats[p.name] = {"kind": "ascii", "kept": kept, "beyond": beyond, "unknown": unknown}
        if (case_dir / f"{output_prefix}.status").exists():
            add_reader(case_dir / f"{output_prefix}.status", _petar.Status)
        for p in sorted(case_dir.glob(f"{output_prefix}.esc*")):
            if re.fullmatch(rf"{re.escape(output_prefix)}\.esc(\.\d+)?", p.name) and p.is_file():
                add_reader(p, _petar.SingleEscaper, interrupt_mode=mode_i or "none", external_mode=ext_e)
        for p in sorted(case_dir.glob(f"{output_prefix}.group*.n*")):
            m = re.search(r"\.n(\d+)$", p.name)
            if m and p.is_file():
                add_reader(p, _petar.GroupInfo, N=int(m.group(1)), interrupt_mode=mode_i or "none", external_mode=ext_e)
        for name, ctor in [(f"{output_prefix}.lagr", None), (f"{output_prefix}.core", None)]:
            p = case_dir / name
            if not p.exists():
                continue
            if name.endswith(".lagr"):
                add_reader(p, _petar.LagrangianMultiple, external_mode=ext_e)
            else:
                add_reader(p, _petar.Core)
        return stats

    def status_times() -> Optional[List[float]]:
        try:
            tools_path = repo_root / "tools"
            if str(tools_path) not in sys.path:
                sys.path.insert(0, str(tools_path))
            import numpy as _np  # noqa: WPS433
            import analysis as _petar  # noqa: WPS433

            st = _petar.Status()
            st.fromfile(str(case_dir / prefixed_name(output_prefix, "status")))
            return [float(x) for x in _np.asarray(st.time)]
        except Exception:
            return None

    def duplicate_lines(path: Path) -> int:
        seen: Dict[str, int] = {}
        with path.open("r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                key = line.strip()
                if key:
                    seen[key] = seen.get(key, 0) + 1
        return sum(c - 1 for c in seen.values())

    if restart_cfg:
        # Universal restart flow: fresh run to t_end_first, trim event streams
        # with petar.data.clear (tests over-deletion), then the canonical
        # restart from the cleared-to snapshot (tests duplicate records).
        t_first = str(restart_cfg.get("t_end_first", 4))
        clear_t = float(restart_cfg.get("clear_t", 3))
        restart_snap = f"{output_prefix}.{int(clear_t / float(dt_out))}"

        cmd_a = ["petar", "-u", "1", "-t", t_first, "-o", dt_out, "-f", output_prefix, *run_args, "input"]
        rc, warns = run_cmd(cmd_a, cwd=case_dir, env=env, log_path=log_path, timeout_s=to)
        add_step(result, "petar.initial", rc, cmd_a, case_dir, warnings_list=warns, details={"t_end": t_first, "dt_out": dt_out, "output_prefix": output_prefix, "run_args": run_args})
        if rc != 0:
            result["status"] = "fail"
            return result

        pre_stats = stream_stats(clear_t)

        clear_cmd = ["petar.data.clear", "-t", str(clear_t), output_prefix]
        ext_tokens = {x.strip() for x in require.split(",") if x.strip()}
        if "galpy" in ext_tokens:
            clear_cmd += ["--external-mode", "galpy"]
        elif "agama" in ext_tokens:
            clear_cmd += ["--external-mode", "agama"]
        mode_i = case.get("interrupt_mode", "")
        if mode_i and mode_i != "none":
            clear_cmd += ["--interrupt-mode", mode_i]
        rc, warns = run_cmd(clear_cmd, cwd=case_dir, env=env, log_path=log_path, timeout_s=to)
        layout = ""
        with log_path.open("r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("ASCII event layout"):
                    layout = line.strip()
        add_step(result, "petar.data.clear", rc, clear_cmd, case_dir, warnings_list=warns, details={"clear_t": clear_t, "restart_snapshot": restart_snap, "event_layout": layout})
        if rc != 0:
            result["status"] = "fail"
            return result

        # Over-deletion / residue check. Reader-parseable streams: kept records
        # (time<=t) must be preserved exactly and beyond-t records removed.
        # ASCII streams: no growth and no fabricated lines versus the .bk
        # backup the tool leaves behind. snap.lst must hold exactly the
        # entries up to the trim index and the stale snapshot files must be
        # gone.
        post_stats = stream_stats(clear_t)
        clear_issues = []
        for name, before in pre_stats.items():
            after = post_stats.get(name)
            if after is None:
                clear_issues.append(f"{name}: file disappeared after clear")
                continue
            if before.get("kind") in ("reader", "ascii") and after.get("kind") == before.get("kind"):
                if after["beyond"] != 0:
                    clear_issues.append(f"{name}: {after['beyond']} beyond-t records survived")
                if after["kept"] != before["kept"]:
                    clear_issues.append(f"{name}: kept {after['kept']} records, expected {before['kept']} (over-deletion)")
            elif before.get("kind") == "lines":
                if after.get("total", 0) > before.get("total", 0):
                    clear_issues.append(f"{name}: grew {before.get('total')}->{after.get('total')}")
                bk = case_dir / (name + ".bk")
                if bk.exists():
                    import collections

                    pre_c = collections.Counter(l.strip() for l in bk.open("r", encoding="utf-8", errors="ignore") if l.strip())
                    post_c = collections.Counter(l.strip() for l in (case_dir / name).open("r", encoding="utf-8", errors="ignore") if l.strip())
                    fabricated = sum(v for k, v in post_c.items() if pre_c.get(k, 0) < v)
                    if fabricated:
                        clear_issues.append(f"{name}: {fabricated} fabricated line(s) after clear")
        expect_idx = int(clear_t / float(dt_out))
        expect_entries = [f"{output_prefix}.{i}" for i in range(expect_idx + 1)]
        if snap_entries() != expect_entries:
            clear_issues.append(f"{snap_lst_path.name}: entries != expected 0..{expect_idx} after clear")
        stale = [n for n in [f"{output_prefix}.{i}" for i in range(expect_idx + 1, expect_idx + 64)] if (case_dir / n).exists()]
        if stale:
            clear_issues.append(f"stale snapshot files remain: {stale[:3]}")
        add_step(
            result,
            "restart_clear_check",
            0 if not clear_issues else 1,
            ["internal.clear_check", f"t={clear_t}"],
            case_dir,
            details={"streams_before": pre_stats, "streams_after": post_stats, "issues": clear_issues[:10]},
            reason="" if not clear_issues else "; ".join(clear_issues[:3]),
        )
        if clear_issues:
            result["status"] = "fail"
            return result

        # Canonical minimal restart: -p carries every stored option (companions
        # included); only the changed end time is overridden. Re-passing -b or
        # run args would re-apply IC-time settings against the post-restart
        # particle layout (e.g. -b re-marking member pairs that no longer
        # match) — observed to trigger pathological hard-integration grinding.
        cmd_b = ["petar", "-p", f"{output_prefix}.par", "-i", "0", "-t", t_end, restart_snap]
        rc, warns = run_cmd(cmd_b, cwd=case_dir, env=env, log_path=log_path, timeout_s=to)
        add_step(result, "petar.restart", rc, cmd_b, case_dir, warnings_list=warns, details={"t_end": t_end, "dt_out": dt_out, "restart_snapshot": restart_snap})
        if rc != 0:
            result["status"] = "fail"
            return result
        entries_b = snap_entries()
        restart_details = {
            "t_end": t_end,
            "dt_out": dt_out,
            "output_prefix": output_prefix,
            "run_args": run_args,
            "restart_snapshot": restart_snap,
            "snapshots_after_restart": len(entries_b),
        }
        last_idx = int(entries_b[-1].rsplit(".", 1)[-1]) if entries_b else -1
        need_idx = int(float(t_end) / float(dt_out)) - 1
        ok_progress = last_idx >= need_idx
        add_step(
            result,
            "restart_progress_check",
            0 if ok_progress else 1,
            ["internal.restart_check", snap_lst_path.name],
            case_dir,
            details={**restart_details, "last_snapshot_index": last_idx, "required_min_index": need_idx},
            reason="" if ok_progress else f"last snapshot index {last_idx} below required {need_idx}",
        )
        if not ok_progress:
            result["status"] = "fail"
            return result

        # Duplicate-record check after the restart: no repeated snapshot list
        # entries, strictly increasing unique status times, no duplicated
        # event-stream lines, no duplicated group records (time, id1, id2).
        dup_issues = []
        names = snap_entries()
        if len(names) != len(set(names)):
            dup_issues.append(f"{snap_lst_path.name}: repeated entries")
        st_times = status_times() or []
        if len(st_times) != len(set(st_times)) or st_times != sorted(st_times):
            dup_issues.append("status: duplicated or non-monotonic times")
        for pattern in [f"{output_prefix}.sse*", f"{output_prefix}.bse*", f"{output_prefix}.mosse*", f"{output_prefix}.sevn*", f"{output_prefix}.esc_*", f"{output_prefix}.esc.[0-9]*", f"{output_prefix}.interrupt*", f"{output_prefix}.prof.rank*"]:
            for p in sorted(case_dir.glob(pattern)):
                if p.is_file():
                    n_dup = duplicate_lines(p)
                    if n_dup:
                        dup_issues.append(f"{p.name}: {n_dup} duplicated record line(s)")
        for p in sorted(case_dir.glob(f"{output_prefix}.group*.n*")):
            n_dup = duplicate_group_records(repo_root, p)
            if n_dup:
                dup_issues.append(f"{p.name}: {n_dup} duplicated group record(s)")
        add_step(
            result,
            "restart_duplicates_check",
            0 if not dup_issues else 1,
            ["internal.duplicates_check"],
            case_dir,
            details={"issues": dup_issues[:10]},
            reason="" if not dup_issues else "; ".join(dup_issues[:3]),
        )
        if dup_issues:
            result["status"] = "fail"
            return result
    else:
        run_cmdline = ["petar", "-u", "1", "-t", t_end, "-o", dt_out, "-f", output_prefix, *run_args, "input"]
        rc, warns = run_cmd(run_cmdline, cwd=case_dir, env=env, log_path=log_path, timeout_s=to)
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
    rc, warns = run_cmd(process_cmd, cwd=case_dir, env=env, log_path=log_path, timeout_s=to)
    add_step(result, "petar.data.process", rc, process_cmd, case_dir, warnings_list=warns, details={"process_args": process_args})
    if rc != 0:
        result["status"] = "fail"
        return result

    # 5b) petar.group.dedup smoke on the first binary group file (writes
    # <input>.dedup; input is never modified). The tool's layout validation
    # needs >=50 records, so tiny files are recorded as skipped, not failed.
    group_files = sorted(case_dir.glob(f"{output_prefix}.group*.n*"))
    if group_files:
        target = group_files[0]
        records_est = estimate_group_records(repo_root, target)
        dedup_cmd = ["petar.group.dedup", target.name]
        if records_est is not None and records_est < 50:
            add_step(
                result,
                "petar.group.dedup",
                0,
                dedup_cmd,
                case_dir,
                details={"skipped": True, "group_file": target.name, "records_estimate": records_est, "reason": "fewer than 50 records; dedup layout validation requires >=50"},
                reason="group file too small for dedup layout validation",
            )
        else:
            rc, warns = run_cmd(dedup_cmd, cwd=case_dir, env=env, log_path=log_path, timeout_s=to)
            dedup_out = Path(str(target) + ".dedup")
            dedup_ok = rc == 0 and dedup_out.exists() and dedup_out.stat().st_size > 0
            add_step(
                result,
                "petar.group.dedup",
                0 if dedup_ok else 1,
                dedup_cmd,
                case_dir,
                warnings_list=warns,
                details={"group_file": target.name, "dedup_output": dedup_out.name, "output_exists": dedup_out.exists(), "records_estimate": records_est},
                reason="" if dedup_ok else "dedup failed or produced no output",
            )
            if not dedup_ok:
                result["status"] = "fail"
                return result


    # 6) petar.get.object.snap with correct -i/-t for dtype
    object_snap_args = case.get("object_snap_args", [])
    # -C: the tool's core-file default is hardcoded to 'data.core'; cases with a
    # different output prefix (e.g. bse-galpy 'out') must point at <prefix>.core
    snap_cmd = ["petar.get.object.snap", *object_snap_args, "-p", "object", "-C", f"{output_prefix}.core", "-f", "origin", "-m", "id", "1", snap_list.name]
    rc, warns = run_cmd(snap_cmd, cwd=case_dir, env=env, log_path=log_path, timeout_s=to)
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
    rc, warns = run_cmd(format_cmd, cwd=case_dir, env=env, log_path=log_path, timeout_s=to)
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

    # movie's core-file default is hardcoded to 'data.core'; point it at the
    # case's actual output prefix (needed by --cm-mode core)
    if "--core-file" not in movie_args:
        movie_args.extend(["--core-file", f"{output_prefix}.core"])

    movie_cmd = ["petar.movie", *movie_args, snap_list.name]
    rc, warns = run_cmd(movie_cmd, cwd=case_dir, env=env, log_path=log_path, timeout_s=to)
    add_step(result, "petar.movie", rc, movie_cmd, case_dir, warnings_list=warns, details={"movie_args": movie_args})
    if rc != 0:
        result["status"] = "fail"
        result["skip_reason"] = "petar.movie command failed"
        return result

    # 9) Extra per-case tool coverage (read.par / update.par / data.clear /
    # data.gether / format.transfer / find.dt / standalone bse).
    run_tool_smokes(result, case, case_dir, output_prefix, log_path, env, t_end, run_args, timeout_s=to)
    if result["status"] == "fail":
        return result

    expected_files = [
        case_dir / last_snap,
        case_dir / f"{last_snap}.single",
        case_dir / f"{last_snap}.binary",
        case_dir / prefixed_name(output_prefix, "lagr"),
        case_dir / "object.1.origin",
    ]
    for extra in case.get("expect_files", []):
        expected_files.append(case_dir / extra.replace("{prefix}", output_prefix))
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
    to = _step_timeout(matrix)

    log_path = out_root / "build.log"
    rc, warns = run_cmd(["bash", str(script_path)], cwd=repo_root, env=env, log_path=log_path, timeout_s=to)
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
        # Prune by the full matrix case list, not the selection: running a case
        # subset must not delete the other cases' reports.
        stale = prune_stale_case_reports(out_root, [c["name"] for c in matrix.get("cases", [])])
        if stale:
            print(f"[cleanup] removed {len(stale)} stale per-case report(s): {', '.join(stale)}")
        report["preflight"] = preflight_families(repo_root, matrix, selected, out_root)
        for case in selected:
            case_result = run_case(repo_root, case, matrix, out_root)
            report["cases"].append(case_result)
            # Per-case pages show only this case's family preflight entry.
            case_tag = infer_build_tag_from_require(case.get("require", ""))
            case_preflight = {t: e for t, e in report["preflight"].items() if t == case_tag}
            render_html_report(
                {
                    "phase": args.phase,
                    "started_at": report["started_at"],
                    "preflight": case_preflight,
                    "build": None,
                    "cases": [case_result],
                },
                out_root / f"report.functional.{case['name']}.html",
            )

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
