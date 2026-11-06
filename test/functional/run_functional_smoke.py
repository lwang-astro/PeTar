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
        "table.mini{width:100%;margin:0;border-collapse:collapse;}"
        "table.mini th{background:transparent;font-weight:600;white-space:nowrap;width:180px;text-align:left;vertical-align:top;padding:2px 10px 2px 0;font-size:13px;color:#334155;}"
        "table.mini td{background:transparent;border-top:none;padding:2px 0;font-size:13px;}"
        "tr.block-title td{background:#eef2f7;border-top:1px solid #d8dee6;font-weight:700;padding:6px 10px;font-size:13px;color:#16202a;}"
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

    selftest = report.get("schema_selftest")
    if selftest:
        parts.append("<div class=\"card\"><h2>Schema Auto-Detection Selftest (synthetic)</h2>")
        parts.append(f"<p>{badge(selftest.get('status', 'unknown'))} "
                     "<span class=\"muted\">same-tier ambiguity raises; correct rare-flag sidecar accepted; "
                     "wrong sidecar rejected; explicit kwargs win over auto-fill "
                     "(see doc/output-schema-tagging-plan.md)</span></p>")
        parts.append("<table><thead><tr><th>Check</th><th>Result</th><th>Detail</th></tr></thead><tbody>")
        for chk in selftest.get("checks", []):
            chk_status = "pass" if chk.get("ok") else "fail"
            parts.append("<tr>")
            parts.append(f"<td>{esc(chk.get('name', ''))}</td>")
            parts.append(f"<td class=\"{chk_status}\">{chk_status}</td>")
            parts.append(f"<td><code>{esc(chk.get('detail', ''))}</code></td>")
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
                # 'variants'/'matrix' get their own tables below
                # (lagr.variants / format.matrix steps)
                if key in ("variants", "matrix"):
                    continue
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
        schema_checks = [chk for chk in read_checks if str(chk.get("label", "")).startswith("schema:")]
        other_checks = [chk for chk in read_checks if not str(chk.get("label", "")).startswith("schema:")]

        def _shorten_value(value: Any) -> str:
            """Compact a detail value: numpy class reprs -> short names,
            dict-style layout reprs -> k=v lists."""
            text = str(value)
            text = re.sub(r"<class 'numpy\.([^']+)'>", r"\1", text)
            text = re.sub(r"<class '([^']+)'>", r"\1", text)
            return text

        def _parse_layout(text: str) -> Dict[str, str]:
            layout: Dict[str, str] = {}
            for key, raw, sval in re.findall(r"'(\w+)':\s*('([^']*)'|<class '[^']*'>|True|False|[0-9.]+)", text):
                if sval:
                    layout[key] = sval
                elif "class" in raw:
                    layout[key] = raw.split("'")[1].replace("numpy.", "")
                else:
                    layout[key] = raw
            return layout

        def _fmt_schema_warning(msg: str) -> str:
            """One readable line for the tier-resolution schema notice:
            evidence, the resolved default layout, and each rarer
            alternative as ONLY the flags differing from the default."""
            m = re.match(
                r"Schema detection of (\S+) \(([^)]*)\) resolved to the default-build layout (\{.*?\}); "
                r"the file could also match rarer builds (\[.*?\])\.(.*)$",
                msg,
            )
            if not m:
                return _shorten_value(msg)
            path, evidence, resolved_s, alts_s, tail = m.groups()
            resolved = _parse_layout(resolved_s)
            alt_summaries = []
            for alt_s in re.findall(r"\{.*?\}", alts_s):
                alt = _parse_layout(alt_s)
                diff = ["%s=%s" % (k, alt[k]) for k in sorted(alt) if alt.get(k) != resolved.get(k)]
                alt_summaries.append("+".join(diff) if diff else "identical")
            head = "Schema: %s (%s) → default %s" % (
                Path(path).name,
                evidence,
                ", ".join("%s=%s" % kv for kv in sorted(resolved.items())),
            )
            if alt_summaries:
                head += "; rarer alternative: " + " | ".join(alt_summaries)
            tail = tail.strip()
            return head + ("" if not tail else " — " + tail)

        def render_check_rows(chk_list: List[Dict[str, Any]]) -> None:
            for chk in chk_list:
                chk_status = "pass" if chk.get("ok") else "fail"
                details = chk.get("details", {}) or {}
                parts.append("<tr>")
                parts.append(f"<td>{esc(chk.get('label', ''))}</td>")
                parts.append(f"<td><code>{esc(chk.get('command', ''))}</code></td>")
                parts.append(f"<td class=\"{chk_status}\">{chk_status}</td>")
                parts.append("</tr>")
                # structured details: mini key/value table instead of a
                # flat pre-formatted dump; warnings deduplicated and
                # schema notices compacted to their decisive content
                parts.append('<tr class="step-details"><td colspan="3"><table class="mini">')
                if chk.get("file"):
                    note = chk.get("file_note")
                    note_html = (" &mdash; " + esc(str(note))) if note else ""
                    parts.append(f'<tr><th>file</th><td><span title="{esc(chk["file"])}">{esc(Path(chk["file"]).name)}</span>{note_html}</td></tr>')
                if chk.get("error"):
                    parts.append(f'<tr><th>error</th><td class="fail">{esc(_shorten_value(chk["error"]))}</td></tr>')
                for key, value in details.items():
                    text = _shorten_value(value)
                    cmd_like = key in ("explicit", "auto", "command", "autofill", "snapshot") or "petar." in text[:24]
                    cell = f'<code>{esc(text)}</code>' if cmd_like else esc(text)
                    cls = ' class="pass"' if key == "comparison" and text.startswith("identical") else ""
                    parts.append(f"<tr><th>{esc(key)}</th><td{cls}>{cell}</td></tr>")
                warn_msgs = list(dict.fromkeys(chk.get("warnings", []) or []))
                for wmsg in warn_msgs:
                    parts.append(f'<tr><th>warning</th><td>{esc(_fmt_schema_warning(wmsg))}</td></tr>')
                parts.append("</table></td></tr>")

        if schema_checks:
            parts.append("<h3>Schema Auto-Detection Checks (data.meta sidecar + structure sniff)</h3><table><thead><tr><th>Label</th><th>Command</th><th>Result</th></tr></thead><tbody>")
            render_check_rows(schema_checks)
            parts.append("</tbody></table>")

        if other_checks:
            parts.append("<h3>Python Read Checks</h3><table><thead><tr><th>Label</th><th>Command</th><th>Result</th></tr></thead><tbody>")
            render_check_rows(other_checks)
            parts.append("</tbody></table>")

        # Variant-carrying steps go last: they build on the default
        # data.process run and the read checks above
        for vstep in case.get("steps", []):
            if vstep.get("name") == "format.matrix":
                fmt_matrix = (vstep.get("details") or {}).get("matrix")
                if fmt_matrix:
                    parts.append("<h3>Format Coverage (per file: binary / ASCII / npy)</h3>")
                    parts.append("<table><thead><tr><th>File</th><th>Binary</th><th>ASCII</th><th>npy</th></tr></thead><tbody>")
                    row_labels = [
                        ("snapshot", "snapshot <code>data.&lt;N&gt;</code>"),
                        ("lagr", "Lagrangian <code>data.lagr</code>"),
                        ("core", "core <code>data.core</code>"),
                    ]
                    fmt_order = ("binary", "ascii", "npy")
                    for key, label in row_labels:
                        parts.append(f"<tr><td>{label}</td>")
                        for fmt in fmt_order:
                            c = (fmt_matrix.get(key) or {}).get(fmt)
                            if c is None:
                                parts.append('<td><span class="muted">—</span></td>')
                                continue
                            status = "pass" if c.get("ok") else "fail"
                            cmd = esc(c.get("command", ""))
                            parts.append(
                                f'<td class="{status}"><b>{status}</b> · {esc(c.get("detail", ""))}'
                                + (f"<br><code>{cmd}</code>" if cmd else "")
                                + "</td>"
                            )
                        parts.append("</tr>")
                    parts.append("</tbody></table>")
                continue
            heading = {
                "lagr.variants": "Lagrangian Variant Checks (petar.data.process with non-default "
                                 "--mass-fraction / --add-star-type / --add-mass-range / --calc-energy / -o ascii)",
            }.get(vstep.get("name", ""))
            if heading is None:
                continue
            lagr_variant_checks = (vstep.get("details") or {}).get("variants", [])
            if not lagr_variant_checks:
                continue
            parts.append(f"<h3>{esc(heading)}</h3>")
            # the check name is "<variant>: <check>": the variant is the test
            # subject and becomes a block TITLE row (full width); the table
            # columns then carry only the short per-check remainder
            parts.append("<table><thead><tr><th>Check</th><th>Command</th><th>Result</th><th>Detail</th></tr></thead><tbody>")
            blocks: Dict[str, List[Dict[str, Any]]] = {}
            block_order: List[str] = []
            for chk in lagr_variant_checks:
                name = str(chk.get("name", ""))
                vname, _, sub = name.partition(":")
                if vname not in blocks:
                    blocks[vname] = []
                    block_order.append(vname)
                entry = dict(chk)
                entry["_sub"] = sub.strip()
                blocks[vname].append(entry)
            for vname in block_order:
                parts.append(f'<tr class="block-title"><td colspan="4">{esc(vname)}</td></tr>')
                for chk in blocks[vname]:
                    chk_status = "pass" if chk.get("ok") else "fail"
                    parts.append("<tr>")
                    parts.append(f"<td>{esc(chk.get('_sub', ''))}</td>")
                    parts.append(f"<td><code>{esc(chk.get('command', ''))}</code></td>")
                    parts.append(f"<td class=\"{chk_status}\">{chk_status}</td>")
                    parts.append(f"<td><code>{esc(chk.get('detail', ''))}</code></td>")
                    parts.append("</tr>")
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

    # Format-contract assertion helper (shared with the SDAR binary reader):
    # fromfile() uses the same sniff to auto-dispatch ASCII files to loadtxt().
    detect_ascii = getattr(petar, "detectAsciiFile", None)
    if detect_ascii is None:
        raise RuntimeError(
            "petar.detectAsciiFile is unavailable: reinstall the SDAR Python tools "
            "(make -C SDAR/tools) so the functional format assertions can run"
        )

    def add_check(label: str, path: Path, loader, command: str, strict_mismatch: bool = True, expect_format: str = "binary") -> None:
        if not path.exists():
            return
        with warnings.catch_warnings(record=True) as rec:
            warnings.simplefilter("always")
            ok = True
            error = ""
            details: Dict[str, Any] = {}
            loaded = None
            try:
                # Format contract assertion: the physical format (text vs binary)
                # must match the writer's contract, so a writer format flip or a
                # rename fails loudly here instead of surfacing as reader
                # confusion downstream. Empty files carry no format signal.
                sniffed_text = None
                if path.stat().st_size > 0:
                    sniffed_text = detect_ascii(str(path))
                    details["format"] = "ascii" if sniffed_text else "binary"
                    if expect_format == "binary" and sniffed_text:
                        raise ValueError(f"expected BINARY output, but the file sniffs as printable ASCII text (writer format flip?): {path.name}")
                    if expect_format == "ascii" and not sniffed_text:
                        raise ValueError(f"expected ASCII text output, but the file sniffs as binary content (writer format flip?): {path.name}")
                loader_result = loader(path)
                if isinstance(loader_result, dict):
                    details.update(loader_result.get("details", {}) or {})
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

        # Make the read concrete for the report reader: the 'path' token in
        # the API-call strings IS the checked file itself (see the file row)
        fname = "'%s'" % Path(str(path)).name
        command = re.sub(r"\bpath\b", fname, command)
        for mkey in ("explicit", "auto", "snapshot"):
            mval = details.get(mkey)
            if isinstance(mval, str):
                details[mkey] = re.sub(r"\bpath\b", fname, mval)

        warn_msgs = [str(item.message) for item in rec]
        mismatch_warns = [
            msg
            for msg in warn_msgs
            if re.search(r"column|dtype|aligned|truncated|itemsize|Binary file size", msg, flags=re.IGNORECASE)
            # the schema tier-resolution notice ("Schema detection of ... resolved
            # to the default-build layout") is informational, not a misalignment
            and not msg.startswith("Schema detection of")
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
        # Uniform concise-constructor read (same pattern as every other
        # stream): petar.Profile(path) detects the column layout from the
        # numeric body; legacy files carrying a title line are handled by
        # the generic path-constructor title skip (no prof-specific logic).
        prof = petar.Profile(str(path))
        return {
            "command": "petar.Profile(path)",
            "details": {"records": int(prof.size), "ncols": int(prof.ncols)},
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

    def compare_members(a, b):
        """Member-wise comparison of two DictNpArrayMix reads (walks sub-members)."""
        if a.size != b.size:
            return 'size mismatch: explicit %d vs auto %d' % (a.size, b.size)
        diffs = []

        def walk(x, y, prefix=""):
            for key, _par in x.keys:
                xa = getattr(x, key)
                ya = getattr(y, key)
                if isinstance(xa, np.ndarray):
                    if not np.array_equal(xa, ya):
                        diffs.append(prefix + key)
                elif hasattr(xa, "keys"):
                    walk(xa, ya, prefix + key + ".")

        walk(a, b)
        if not diffs:
            return 'identical (all %d columns x %d records)' % (a.ncols, a.size)
        return 'DIFFER: ' + ','.join(diffs[:6])

    def paired_loader(desc, load_explicit, make_auto, auto_desc, make_snapshot=None, snapshot_desc=None):
        """Read the same file twice (three times for snapshot-family files)
        and report the side-by-side comparison: the legacy explicit-kwargs
        read vs the concise auto-detection read (`petar.<Reader>(path)`),
        plus - when ``make_snapshot`` is given - the dedicated
        `petar.Snapshot(path)` path (particle + header in one object). The
        comparisons (detected kwargs source + member-wise equality) are
        carried in the check details, rendered in the HTML report."""

        def load(path: Path):
            obj_a = load_explicit(str(path))
            details: Dict[str, Any] = {}
            # when the dedicated Snapshot path is provided it is the PRIMARY
            # (default) reading way: it owns the check's command field and
            # the first details rows; the legacy explicit and auto reads
            # follow as comparison methods
            if make_snapshot is not None:
                try:
                    snap_obj = make_snapshot(str(path), {})
                    header_ok = (int(snap_obj.header.n) == int(snap_obj.particle.size))
                    details["snapshot"] = snapshot_desc
                    details["snapshot_result"] = "ok, %d records, time=%.17g" % (snap_obj.particle.size, snap_obj.header.time)
                    details["snapshot_comparison"] = "%s%s" % (
                        compare_members(obj_a, snap_obj.particle),
                        "; header N == records" if header_ok else "; HEADER N MISMATCH: %d vs %d" % (snap_obj.header.n, snap_obj.particle.size))
                    if not header_ok:
                        raise ValueError("snapshot header N (%d) != particle records (%d)" % (snap_obj.header.n, snap_obj.particle.size))
                except Exception as exc:  # noqa: BLE001
                    details["snapshot"] = snapshot_desc
                    details["snapshot_result"] = "failed: %s" % exc
                    details["snapshot_comparison"] = "SNAPSHOT READ FAILED"
                    raise
            details["explicit"] = desc
            details["explicit_result"] = "ok, %d records, %d columns" % (obj_a.size, obj_a.ncols)
            try:
                layout_keys = ("interrupt_mode", "external_mode", "spin_3d", "use_mpfrc", "collect_sp_acc", "float_type")
                # kwargs come from the READER CLASS's own resolve machinery
                # (not a fixed model class), so event/group/escaper readers
                # detect with their own LAYOUT_PARAMS declarations
                _d, kw, _p = type(obj_a).resolvePathInput(str(path), {})
                provenance = 'reader machinery (%s)' % type(obj_a).__name__
                obj_b = make_auto(str(path), kw)
                details["auto"] = auto_desc
                details["kwargs_source"] = (
                    "%s - getDataKwargs(path) locates <prefix>.meta from the data path "
                    "(here the sidecar of the producing run); when absent it falls back "
                    "to structure sniffing (header width + column/record layout)" % provenance)
                details["detected_kwargs"] = ", ".join("%s=%s" % (k, kw[k]) for k in layout_keys if k in kw)
                details["auto_result"] = "ok, %d records" % obj_b.size
                details["comparison"] = compare_members(obj_a, obj_b)
            except Exception as exc:  # noqa: BLE001 - report, do not mask the explicit read
                details["auto"] = auto_desc
                details["auto_result"] = "failed: %s" % exc
                details["comparison"] = "AUTO READ FAILED"
                raise
            details["records"] = int(obj_a.size)
            # the command field shows the PRIMARY reading way (the dedicated
            # Snapshot path when provided, else the legacy explicit read);
            # every method is also reported as its own details row - a
            # single "|" concatenated line wraps unreadably in the HTML table
            return {
                "command": snapshot_desc if make_snapshot is not None else desc,
                "details": details,
                "_object": obj_a,
            }

        return load

    add_check(
        f"{output_prefix}.lagr",
        case_dir / prefixed_name(output_prefix, "lagr"),
        paired_loader(
            f"petar.LagrangianMultiple(external_mode={external_mode}).fromfile(path)",
            obj_loader(petar.LagrangianMultiple, external_mode=external_mode),
            lambda p, kw: petar.LagrangianMultiple(p),
            "petar.LagrangianMultiple(path)",
        ),
        f"petar.LagrangianMultiple(external_mode={external_mode}).fromfile(path)",
    )
    add_check(
        f"{output_prefix}.core",
        case_dir / prefixed_name(output_prefix, "core"),
        paired_loader(
            "petar.Core().fromfile(path)",
            obj_loader(petar.Core),
            lambda p, kw: petar.Core(p),
            "petar.Core(path)",
        ),
        "petar.Core().fromfile(path)",
    )
    add_check(f"{output_prefix}.status", case_dir / prefixed_name(output_prefix, "status"), obj_loader(petar.Status), "petar.Status().fromfile(path)")
    add_check(
        f"{output_prefix}.esc_single",
        case_dir / prefixed_name(output_prefix, "esc_single"),
        paired_loader(
            f"petar.SingleEscaper(interrupt_mode={interrupt_mode}, external_mode={external_mode}).fromfile(path)",
            obj_loader(petar.SingleEscaper, interrupt_mode=interrupt_mode, external_mode=external_mode),
            lambda p, kw: petar.SingleEscaper(p),
            "petar.SingleEscaper(path)",
        ),
        f"petar.SingleEscaper(interrupt_mode={interrupt_mode}, external_mode={external_mode}).fromfile(path)",
    )
    add_check(
        f"{output_prefix}.esc_binary",
        case_dir / prefixed_name(output_prefix, "esc_binary"),
        paired_loader(
            f"petar.BinaryEscaper(interrupt_mode={interrupt_mode}, external_mode={external_mode}).fromfile(path)",
            obj_loader(petar.BinaryEscaper, interrupt_mode=interrupt_mode, external_mode=external_mode),
            lambda p, kw: petar.BinaryEscaper(p),
            "petar.BinaryEscaper(path)",
        ),
        f"petar.BinaryEscaper(interrupt_mode={interrupt_mode}, external_mode={external_mode}).fromfile(path)",
    )

    add_check(
        f"{output_prefix}.sse",
        case_dir / prefixed_name(output_prefix, "sse"),
        paired_loader(
            "petar.SSEType().loadtxt(path)",
            obj_loader(petar.SSEType, "loadtxt"),
            lambda p, kw: petar.SSEType(p),
            "petar.SSEType(path)",
        ),
        "petar.SSEType().loadtxt(path)",
        expect_format="ascii",
    )

    add_check(
        f"{output_prefix}.sse.type_change",
        case_dir / prefixed_name(output_prefix, "sse.type_change"),
        paired_loader(
            f"petar.SSETypeChange(interrupt_mode={interrupt_mode}).loadtxt(path)",
            se_aware(petar.SSETypeChange, interrupt_mode),
            lambda p, kw: petar.SSETypeChange(p),
            "petar.SSETypeChange(path)",
        ),
        f"petar.SSETypeChange(interrupt_mode={interrupt_mode}).loadtxt(path)",
        expect_format="ascii",
    )
    add_check(
        f"{output_prefix}.sse.sn_kick",
        case_dir / prefixed_name(output_prefix, "sse.sn_kick"),
        paired_loader(
            f"petar.SSESNKick(interrupt_mode={interrupt_mode}).loadtxt(path)",
            se_aware(petar.SSESNKick, interrupt_mode),
            lambda p, kw: petar.SSESNKick(p),
            "petar.SSESNKick(path)",
        ),
        f"petar.SSESNKick(interrupt_mode={interrupt_mode}).loadtxt(path)",
        expect_format="ascii",
    )

    add_check(
        f"{output_prefix}.bse",
        case_dir / prefixed_name(output_prefix, "bse"),
        paired_loader(
            "petar.BSEType().loadtxt(path)",
            obj_loader(petar.BSEType, "loadtxt"),
            lambda p, kw: petar.BSEType(p),
            "petar.BSEType(path)",
        ),
        "petar.BSEType().loadtxt(path)",
        expect_format="ascii",
    )
    add_check(
        f"{output_prefix}.bse_status",
        case_dir / prefixed_name(output_prefix, "bse_status"),
        paired_loader(
            "petar.BSEStatus().fromfile(path)",
            obj_loader(petar.BSEStatus),
            lambda p, kw: petar.BSEStatus(p),
            "petar.BSEStatus(path)",
        ),
        "petar.BSEStatus().fromfile(path)",
    )
    add_check(
        f"{output_prefix}.bse.type_change",
        case_dir / prefixed_name(output_prefix, "bse.type_change"),
        paired_loader(
            f"petar.BSETypeChange(interrupt_mode={interrupt_mode}).loadtxt(path)",
            se_aware(petar.BSETypeChange, interrupt_mode),
            lambda p, kw: petar.BSETypeChange(p),
            "petar.BSETypeChange(path)",
        ),
        f"petar.BSETypeChange(interrupt_mode={interrupt_mode}).loadtxt(path)",
        expect_format="ascii",
    )
    add_check(
        f"{output_prefix}.bse.sn_kick",
        case_dir / prefixed_name(output_prefix, "bse.sn_kick"),
        paired_loader(
            f"petar.BSEKick(interrupt_mode={interrupt_mode}).loadtxt(path)",
            se_aware(petar.BSEKick, interrupt_mode),
            lambda p, kw: petar.BSEKick(p),
            "petar.BSEKick(path)",
        ),
        f"petar.BSEKick(interrupt_mode={interrupt_mode}).loadtxt(path)",
        expect_format="ascii",
    )
    add_check(
        f"{output_prefix}.bse.gw_kick",
        case_dir / prefixed_name(output_prefix, "bse.gw_kick"),
        paired_loader(
            f"petar.BSEKick(interrupt_mode={interrupt_mode}).loadtxt(path)",
            se_aware(petar.BSEKick, interrupt_mode),
            lambda p, kw: petar.BSEKick(p),
            "petar.BSEKick(path)",
        ),
        f"petar.BSEKick(interrupt_mode={interrupt_mode}).loadtxt(path)",
        expect_format="ascii",
    )
    add_check(
        f"{output_prefix}.bse.dynamic_merge",
        case_dir / prefixed_name(output_prefix, "bse.dynamic_merge"),
        paired_loader(
            f"petar.BSEDynamicMerge(interrupt_mode={interrupt_mode}).loadtxt(path)",
            se_aware(petar.BSEDynamicMerge, interrupt_mode),
            lambda p, kw: petar.BSEDynamicMerge(p),
            "petar.BSEDynamicMerge(path)",
        ),
        f"petar.BSEDynamicMerge(interrupt_mode={interrupt_mode}).loadtxt(path)",
        expect_format="ascii",
    )
    add_check(
        f"{output_prefix}.bse.binary_merge",
        case_dir / prefixed_name(output_prefix, "bse.binary_merge"),
        paired_loader(
            f"petar.BSETypeChange(interrupt_mode={interrupt_mode}).loadtxt(path)",
            se_aware(petar.BSETypeChange, interrupt_mode),
            lambda p, kw: petar.BSETypeChange(p),
            "petar.BSETypeChange(path)",
        ),
        f"petar.BSETypeChange(interrupt_mode={interrupt_mode}).loadtxt(path)",
        expect_format="ascii",
    )

    # Snapshot files listed in [prefix].snap.lst.
    snap_list = case_dir / prefixed_name(output_prefix, "snap.lst")
    if snap_list.exists():
        snap_names = [line.strip() for line in snap_list.read_text(encoding="utf-8").splitlines() if line.strip()]
        if snap_names:


            def read_particle_explicit(p, im=interrupt_mode, em=external_mode):
                obj = petar.Particle(interrupt_mode=im, external_mode=em)
                obj.fromfile(p)
                return obj

            def read_particle_auto(p, kw):
                return petar.Particle(p)

            def read_binary_explicit(p, im=interrupt_mode, em=external_mode):
                obj = petar.Binary(member_particle_type=petar.Particle, interrupt_mode=im, external_mode=em, G=petar.G_MSUN_PC_MYR)
                obj.fromfile(p)
                return obj

            def read_binary_auto(p, kw):
                return petar.Binary(p, G=petar.G_MSUN_PC_MYR)

            snap_name = snap_names[-1]
            snap_path = case_dir / snap_name
            offset = petar.HEADER_OFFSET_WITH_CM if external_mode != "none" else petar.HEADER_OFFSET
            offset_name = "HEADER_OFFSET_WITH_CM" if external_mode != "none" else "HEADER_OFFSET"

            def snap_particle(path: Path, off=offset, im=interrupt_mode, em=external_mode):
                obj = petar.Particle(interrupt_mode=im, external_mode=em)
                obj.fromfile(str(path), offset=off)
                return obj

            def snap_particle_auto(p, kw):
                return petar.Particle(p)

            def snap_snapshot(p, kw):
                return petar.Snapshot(p)

            add_check(
                f"snapshot:{snap_name}",
                snap_path,
                paired_loader(
                    f"petar.Particle(interrupt_mode={interrupt_mode}, external_mode={external_mode}).fromfile(path, offset=petar.{offset_name})",
                    lambda p: snap_particle(p),
                    snap_particle_auto,
                    "petar.Particle(path)",
                    make_snapshot=snap_snapshot,
                    snapshot_desc="petar.Snapshot(path) — default (particle + header)",
                ),
                "petar.Snapshot(path)",
            )
            add_check(
                f"snapshot.single:{snap_name}.single",
                case_dir / f"{snap_name}.single",
                paired_loader(
                    f"petar.Particle(interrupt_mode={interrupt_mode}, external_mode={external_mode}).fromfile(path)",
                    lambda p, im=interrupt_mode, em=external_mode: read_particle_explicit(p, im, em),
                    read_particle_auto,
                    "petar.Particle(path)",
                ),
                f"petar.Particle(interrupt_mode={interrupt_mode}, external_mode={external_mode}).fromfile(path)",
            )
            add_check(
                f"snapshot.binary:{snap_name}.binary",
                case_dir / f"{snap_name}.binary",
                paired_loader(
                    f"petar.Binary(member_particle_type=petar.Particle, interrupt_mode={interrupt_mode}, external_mode={external_mode}, G=petar.G_MSUN_PC_MYR).fromfile(path)",
                    lambda p, im=interrupt_mode, em=external_mode: read_binary_explicit(p, im, em),
                    read_binary_auto,
                    "petar.Binary(path, G=petar.G_MSUN_PC_MYR)",
                ),
                f"petar.Binary(member_particle_type=petar.Particle, interrupt_mode={interrupt_mode}, external_mode={external_mode}, G=petar.G_MSUN_PC_MYR).fromfile(path)",
            )
            # ASCII twin of the last snapshot: exercises the auto reader's
            # format sniffing end to end (byte sniff -> header-line parse ->
            # loadtxt(skiprows=1)) against the same explicit-kwargs read in
            # text form. Written to a temp dir removed right after the check;
            # getDataKwargs then has no sidecar, so this also proves the
            # structure-sniff fallback path (ASCII header item count + row
            # column count).
            if snap_path.exists():
                import tempfile as _tempfile
                with _tempfile.TemporaryDirectory(dir=str(case_dir)) as _atd:
                    ascii_path = Path(_atd) / f"{snap_name}.ascii"
                    try:
                        _hdr = petar.PeTarDataHeader(str(snap_path), snapshot_format="binary", external_mode=external_mode)
                        _prt = petar.Particle(interrupt_mode=interrupt_mode, external_mode=external_mode)
                        _prt.fromfile(str(snap_path), offset=offset)
                        with ascii_path.open("w", encoding="utf-8") as _fh:
                            _hdr.savetxt(_fh)
                            _prt.savetxt(_fh)
                    except Exception as _exc:  # noqa: BLE001 - report as failed check, keep going
                        add_check(
                            f"snapshot.ascii:{snap_name}.ascii",
                            ascii_path,
                            (lambda _p, _e=_exc: (_ for _ in ()).throw(_e)),
                            f"ASCII twin generation failed: {_e}",
                            expect_format="ascii",
                        )
                    else:

                        def ascii_explicit(p, im=interrupt_mode, em=external_mode):
                            obj = petar.Particle(interrupt_mode=im, external_mode=em)
                            obj.loadtxt(p, skiprows=1)
                            return obj

                        add_check(
                            f"snapshot.ascii:{snap_name}.ascii",
                            ascii_path,
                            paired_loader(
                                f"petar.Particle(interrupt_mode={interrupt_mode}, external_mode={external_mode}).loadtxt(path, skiprows=1)",
                                ascii_explicit,
                                lambda p, kw: petar.Particle(p),
                                "petar.Particle(path)",
                            ),
                            f"petar.Particle(interrupt_mode={interrupt_mode}, external_mode={external_mode}).loadtxt(path, skiprows=1)",
                            expect_format="ascii",
                        )

            # Early-interval binary readback: the last snapshot can be binary-free
            # by design (injected compact binaries merge via BSE; std/merger ICs
            # carry none), which would make the binary reader check trivial.
            if len(snap_names) >= 2:
                early_name = snap_names[1]
                add_check(
                    f"snapshot.binary.early:{early_name}.binary",
                    case_dir / f"{early_name}.binary",
                    paired_loader(
                        f"petar.Binary(member_particle_type=petar.Particle, interrupt_mode={interrupt_mode}, external_mode={external_mode}, G=petar.G_MSUN_PC_MYR).fromfile(path)",
                        lambda p, im=interrupt_mode, em=external_mode: read_binary_explicit(p, im, em),
                        read_binary_auto,
                        "petar.Binary(path, G=petar.G_MSUN_PC_MYR)",
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
            expect_format="ascii",
        )

    # Validate profiling outputs (ASCII table).
    for prof_path in sorted(case_dir.glob(f"{output_prefix}.prof.rank.*")):
        if prof_path.name.endswith(".bk"):
            # petar.data.clear backup of the trimmed stream - not a live output
            continue
        add_check(
            f"profile:{prof_path.name}",
            prof_path,
            load_profile_with_matching_schema,
            "petar.Profile(path)",
            expect_format="ascii",
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
            paired_loader(
                f"petar.GroupInfo(N={n_member}, interrupt_mode={interrupt_mode}, external_mode={external_mode}).fromfile(path)",
                obj_loader(petar.GroupInfo, N=n_member, interrupt_mode=interrupt_mode, external_mode=external_mode),
                lambda p, kw, n=n_member: petar.GroupInfo(p, N=n),
                f"petar.GroupInfo(path, N={n_member})",
            ),
            f"petar.GroupInfo(N={n_member}, interrupt_mode={interrupt_mode}, external_mode={external_mode}).fromfile(path)",
        )

    # Schema sidecar + auto-detection checks (doc/output-schema-tagging-plan.md):
    # 1) the solver-written <prefix>.meta must exist and agree with the case's
    #    configure modes;
    # 2) structure sniffing (getDataKwargs) must recover the same schema;
    # 3) the CLI auto-fill path (autoFillKwargs from empty options) must fill
    #    the same layout kwargs;
    # 4) Snapshot (auto-detected) must equal the explicit-kwargs read.
    schema_api = {
        "getDataKwargs": getattr(petar, "getDataKwargs", None),
        "readMetaData": getattr(petar, "readMetaData", None),
        "autoFillKwargs": getattr(petar, "autoFillKwargs", None),
        "Snapshot": getattr(petar, "Snapshot", None),
    }
    bse_family = ("bse", "mobse", "bseEmp")

    def interrupt_family_matches(value: str, expected: str) -> bool:
        if value == expected:
            return True
        return expected == "bse" and value in bse_family

    def last_snapshot_path() -> Optional[Path]:
        candidates = [
            p
            for p in case_dir.glob(f"{output_prefix}.[0-9]*")
            if re.fullmatch(rf"{re.escape(output_prefix)}\.\d+", p.name)
        ]
        if not candidates:
            return None
        return max(candidates, key=lambda p: int(p.name.rsplit(".", 1)[1]))

    snap_for_schema = last_snapshot_path()
    if None in schema_api.values():
        # Old installed petar package: report as a failed check (environment
        # gap) rather than silently skipping, but only when snapshots exist.
        if snap_for_schema is not None:
            checks.append(
                {
                    "label": "schema:api",
                    "file": str(snap_for_schema),
                    "command": "petar.getDataKwargs/readMetaData/autoFillKwargs/Snapshot",
                    "ok": False,
                    "error": "installed petar package predates the schema-detection API; reinstall with make install",
                    "warnings": [],
                    "details": {},
                }
            )
    elif snap_for_schema is not None:
        snap_fmt = "ascii" if detect_ascii(str(snap_for_schema)) else "binary"

        def load_schema_meta(path: Path):
            meta = schema_api["readMetaData"](str(path))
            if meta is None:
                raise ValueError(f"{output_prefix}.meta not found next to {path.name}: the solver build predates the sidecar writer")
            details = {"meta": {k: meta[k] for k in sorted(meta)}}
            meta_im = meta.get("layout.interrupt_mode", meta.get("interrupt_mode", "none"))
            meta_ext = meta.get("layout.external_mode", meta.get("external_mode", "none"))
            if not interrupt_family_matches(meta_im, interrupt_mode):
                raise ValueError(f"data.meta interrupt_mode={meta_im} does not match the case mode {interrupt_mode}")
            if (meta_ext != "none") != (external_mode != "none"):
                raise ValueError(f"data.meta external_mode={meta_ext} does not match the case mode {external_mode}")
            return {
                "command": "petar.readMetaData(path)",
                "details": details,
                "_object": None,
            }

        def load_schema_detect(path: Path):
            kw = schema_api["getDataKwargs"](str(path))
            det_im = kw.get("interrupt_mode", "none")
            det_ext = kw.get("external_mode", "none")
            if not interrupt_family_matches(det_im, interrupt_mode):
                raise ValueError(f"detected interrupt_mode={det_im} does not match the case mode {interrupt_mode}")
            if (det_ext != "none") != (external_mode != "none"):
                raise ValueError(f"detected external_mode={det_ext} does not match the case mode {external_mode}")
            # round-trip: auto-detected read vs explicit-kwargs read
            auto_snap = schema_api["Snapshot"](str(path))
            auto_p = auto_snap.particle
            explicit = petar.Particle(interrupt_mode=interrupt_mode, external_mode=external_mode)
            if snap_fmt == "ascii":
                explicit.loadtxt(str(path), skiprows=1)
            else:
                explicit.fromfile(
                    str(path),
                    offset=petar.HEADER_OFFSET_WITH_CM if external_mode != "none" else petar.HEADER_OFFSET,
                )
            if auto_p.size != explicit.size:
                raise ValueError(f"detected-kwargs read size {auto_p.size} != explicit-kwargs read size {explicit.size}")
            for key, _ in explicit.keys:
                a, b = getattr(auto_p, key), getattr(explicit, key)
                if isinstance(a, np.ndarray) and not np.array_equal(a, b):
                    raise ValueError(f"detected-kwargs read differs from explicit-kwargs read in member '{key}'")
            return {
                "command": "petar.getDataKwargs(path) + petar.Snapshot(path) vs explicit kwargs",
                "details": {"detected": {k: str(v) for k, v in kw.items()}, "records": int(auto_p.size)},
                "_object": auto_p,
            }

        def load_schema_autofill(path: Path):
            kw: Dict[str, Any] = {}
            schema_api["autoFillKwargs"](str(path), kw, verbose=False)
            if not kw:
                raise ValueError("autoFillKwargs filled nothing from an empty option set")
            fill_im = kw.get("interrupt_mode", "none")
            fill_ext = kw.get("external_mode", "none")
            if not interrupt_family_matches(fill_im, interrupt_mode):
                raise ValueError(f"auto-filled interrupt_mode={fill_im} does not match the case mode {interrupt_mode}")
            if (fill_ext != "none") != (external_mode != "none"):
                raise ValueError(f"auto-filled external_mode={fill_ext} does not match the case mode {external_mode}")
            return {
                "command": "petar.autoFillKwargs(path, {})",
                "details": {"filled": {k: str(v) for k, v in kw.items()}},
                "_object": None,
            }

        add_check("schema:meta", snap_for_schema, load_schema_meta, "petar.readMetaData(path)", expect_format=snap_fmt)
        add_check("schema:detect", snap_for_schema, load_schema_detect, "petar.getDataKwargs(path)", expect_format=snap_fmt)
        add_check("schema:autofill", snap_for_schema, load_schema_autofill, "petar.autoFillKwargs(path, {})", expect_format=snap_fmt)

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


LAGR_DEFAULT_MASS_FRACTION = [0.1, 0.3, 0.5, 0.7, 0.9]


def _lagr_variant_expected_labels(args: List[str]) -> Dict[str, Any]:
    """Per-run Lagrangian labels a petar.data.process variant passes on the CLI."""
    expected: Dict[str, Any] = {}
    for idx, tok in enumerate(args):
        if tok in ("-m", "--mass-fraction") and idx + 1 < len(args):
            expected["mass_fraction"] = [float(x) for x in args[idx + 1].split(",")]
        elif tok == "--add-star-type" and idx + 1 < len(args):
            expected["add_star_type"] = [x for x in args[idx + 1].split(",") if x]
        elif tok == "--add-mass-range" and idx + 1 < len(args):
            expected["add_mass_range"] = [x for x in args[idx + 1].split(",") if x]
    return expected


def run_lagr_variant_checks(
    repo_root: Path,
    case: Dict[str, Any],
    case_dir: Path,
    env: Dict[str, str],
    log_path: Path,
    process_args: List[str],
    snap_list_name: str,
    timeout_s: Optional[float] = None,
    variants: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    """Lagrangian read/write matrix beyond the default petar.data.process run.

    Variants come from the matrix global 'lagr_variants' (universal defaults,
    run for EVERY case) merged with the case key of the same name (case entries
    override global ones by name, otherwise append). Each entry re-runs
    petar.data.process on the case's snapshot list under a scratch prefix
    (lagrvar<N>) with adjusted
    --mass-fraction / --add-star-type / --add-mass-range (plus the case's own
    mode flags inherited from process_args), then checks:

    1. the tool run itself — or, for an 'expect_fail' entry, the loud failure
       (e.g. the auto escape-radius guard when 0.5 is not among the mass
       fractions and no --r-escape is given),
    2. the zero-warning auto read petar.LagrangianMultiple(path) recovering
       the exact per-run labels (mass_fraction, add_star_type,
       add_mass_range, external_mode) from the <prefix>.lagr.meta sidecar,
       and a shape-strict explicit-kwargs read agreeing on record count
       and the time column,
    3. the legacy read of the bare table copied without its sidecar: a shape
       mismatch must raise (with the no-sidecar guidance warning), while a
       same-shape/different-value table must read successfully but carry the
       no-sidecar mislabel warning (labels fall back to the defaults).

    Returns a step dict; the caller appends it and fails the case when its
    status is not 'pass'.
    """
    import tempfile  # local, mirrors run_schema_selftest

    if variants is None:
        variants = case.get("lagr_variants", [])
    variant_names = [v.get("name", "") for v in variants]
    if len(variant_names) != len(set(variant_names)):
        dupes = sorted({n for n in variant_names if variant_names.count(n) > 1})
        return {
            "name": "lagr.variants",
            "code": 1,
            "status": "fail",
            "reason": f"duplicate lagr_variants names (override-by-name would silently drop work): {', '.join(dupes)}",
            "command": "internal.lagr_variants_validation",
            "cwd": str(case_dir),
            "details": {"variants": []},
        }
    step: Dict[str, Any] = {
        "name": "lagr.variants",
        "code": 0,
        "status": "pass",
        "reason": "",
        "command": f"petar.data.process -p lagrvar<N> (x{len(variants)}; scratch variants under lagrvar_scratch/) ... {snap_list_name}",
        "cwd": str(case_dir),
        "details": {"variants": []},
    }
    if not variants:
        step["details"]["skipped"] = True
        return step

    to = timeout_s if timeout_s is not None and timeout_s > 0 else _step_timeout({})
    tools_path = repo_root / "tools"
    if str(tools_path) not in sys.path:
        sys.path.insert(0, str(tools_path))
    import analysis as petar  # type: ignore
    import numpy as np  # type: ignore

    require_tokens = {x.strip() for x in case.get("require", "").split(",") if x.strip()}
    external_mode = "galpy" if "galpy" in require_tokens else ("agama" if "agama" in require_tokens else "none")

    snap_list = case_dir / snap_list_name
    n_snapshots = sum(1 for line in snap_list.read_text(encoding="utf-8").splitlines() if line.strip())

    # scratch products of earlier runs must not leak into this one
    for stale in case_dir.glob("lagrvar*"):
        if stale.is_file():
            stale.unlink()
    stale_scratch = case_dir / "lagrvar_scratch"
    if stale_scratch.exists():
        shutil.rmtree(stale_scratch)

    def check(name: str, ok: bool, detail: str = "", command: str = "") -> None:
        entry = {"name": name, "ok": bool(ok), "detail": detail}
        if command:
            entry["command"] = command
        step["details"]["variants"].append(entry)
        if not ok:
            step["code"] = 1
            step["status"] = "fail"
            step["reason"] = step["reason"] or f"lagr variant check failed: {name}"

    for vi, variant in enumerate(variants):
        vname = variant.get("name", f"variant{vi}")
        vprefix = f"lagrvar{vi}"
        vargs = resolve_tokens(variant.get("args", []))
        expect_fail_marker = variant.get("expect_fail", "")
        # scan the full effective CLI (case mode flags first, variant args last,
        # matching argparse precedence) so expected labels cannot drift from the
        # actual run if a case ever puts -m/--add-* into process_args
        eff_args = list(process_args) + vargs
        expected = _lagr_variant_expected_labels(eff_args)
        wants_ascii = False
        wants_calc_energy = False
        for idx, tok in enumerate(eff_args):
            if tok in ("-o", "--output-format") and idx + 1 < len(eff_args):
                wants_ascii = eff_args[idx + 1] == "ascii"
            elif tok in ("-e", "--calc-energy"):
                wants_calc_energy = True
        # -o/--output-format also rewrites the DERIVED snapshot files
        # (data.<N>.single/.binary) in place: such variants must run on a
        # scratch copy so the case's canonical outputs stay untouched
        # (matrix flag 'scratch').
        run_dir = case_dir
        if variant.get("scratch", False):
            run_dir = case_dir / "lagrvar_scratch"
            run_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(snap_list, run_dir / snap_list_name)
            for line in snap_list.read_text(encoding="utf-8").splitlines():
                if line.strip():
                    src = case_dir / line.strip()
                    if src.is_file():
                        shutil.copy2(src, run_dir / src.name)
            # solver sidecar is <output-prefix>.meta (e.g. data.snap.lst ->
            # data.meta), NOT <snap-list stem>.meta - derive from the prefix
            # so the scratch run keeps the sidecar-first detection path
            solver_meta = case_dir / (snap_list_name.split(".")[0] + ".meta")
            if solver_meta.is_file():
                shutil.copy2(solver_meta, run_dir / solver_meta.name)
        cmd = ["petar.data.process", "--no-auto-resume", "-p", vprefix, *process_args, *vargs, snap_list_name]
        log_pos = log_path.stat().st_size if log_path.exists() else 0
        rc, warns = run_cmd(cmd, cwd=run_dir, env=env, log_path=log_path, timeout_s=to)
        if warns:
            step.setdefault("warnings", []).extend(warns)
            step.setdefault("_warn_steps", []).append({"step": f"lagr.variants:{vname}", "messages": warns})

        if expect_fail_marker:
            log_segment = ""
            if log_path.exists():
                with log_path.open("r", encoding="utf-8", errors="ignore") as fh:
                    fh.seek(log_pos)
                    log_segment = fh.read()
            found = expect_fail_marker in log_segment
            check(f"{vname}: expected failure with guidance", rc != 0 and found, f"rc={rc}, marker found={found}", command="$ " + " ".join(cmd))
            continue
        if rc != 0:
            check(f"{vname}: petar.data.process exit", False, f"exit code {rc}", command="$ " + " ".join(cmd))
            continue
        check(f"{vname}: petar.data.process exit", True, "rc=0", command="$ " + " ".join(cmd))

        lagr_path = run_dir / f"{vprefix}.lagr"
        meta_path = run_dir / f"{vprefix}.lagr.meta"
        if not lagr_path.exists():
            check(f"{vname}: output table", False, f"{lagr_path.name} missing")
            continue
        check(f"{vname}: output table", True, lagr_path.name)
        check(f"{vname}: sidecar written", meta_path.exists(), meta_path.name if meta_path.exists() else f"{meta_path.name} missing")
        sniffed = petar.detectAsciiFile(str(lagr_path)) if lagr_path.stat().st_size > 0 else None
        check(f"{vname}: output format honored", sniffed == wants_ascii, f"wants ascii={wants_ascii}, sniffed={sniffed}")

        # 2) auto read: new data must resolve cleanly from the sidecar
        disp = str(lagr_path.relative_to(case_dir)) if lagr_path.is_relative_to(case_dir) else str(lagr_path)
        meta_disp = str(meta_path.relative_to(case_dir)) if meta_path.is_relative_to(case_dir) else str(meta_path)
        auto_cmd = f"petar.LagrangianMultiple('{disp}')"
        try:
            with warnings.catch_warnings(record=True) as rec:
                warnings.simplefilter("always")
                auto = petar.LagrangianMultiple(str(lagr_path))
            guidance = [str(w.message) for w in rec if "sidecar" in str(w.message) or "pass them explicitly" in str(w.message)]
            check(f"{vname}: auto read zero-warning", not guidance, "; ".join(guidance) if guidance else f"{int(auto.size)} records", command=auto_cmd)
            labels = dict(getattr(auto, "initargs", {}) or {})
            if "mass_fraction" in expected:
                got = labels.get("mass_fraction")
                got_vals = np.asarray(got, dtype=float) if got is not None else None
                check(
                    f"{vname}: mass_fraction label recovered",
                    got_vals is not None and np.allclose(got_vals, expected["mass_fraction"]),
                    f"got={got}",
                    command=f"auto.initargs['mass_fraction']",
                )
            for key in ("add_star_type", "add_mass_range"):
                if key in expected:
                    got = list(labels.get(key) or [])
                    check(f"{vname}: {key} label recovered", got == expected[key], f"got={got}", command=f"auto.initargs['{key}']")
            check(f"{vname}: external_mode label recovered", str(labels.get("external_mode", external_mode)) == external_mode, f"got={labels.get('external_mode')}", command="auto.initargs['external_mode']")
            if wants_calc_energy:
                check(f"{vname}: calc_energy label recovered", bool(labels.get("calc_energy")), f"got={labels.get('calc_energy')}", command="auto.initargs['calc_energy']")
            check(f"{vname}: record count == snapshots", int(auto.size) == n_snapshots, f"records={int(auto.size)} snapshots={n_snapshots}", command=auto_cmd)

            # explicit read with the known labels: shape-strict (fromfile
            # raises on mismatch) - compare record count and the time column
            explicit = petar.LagrangianMultiple(
                mass_fraction=labels.get("mass_fraction"),
                add_star_type=labels.get("add_star_type"),
                add_mass_range=labels.get("add_mass_range"),
                calc_energy=labels.get("calc_energy", wants_calc_energy),
                external_mode=external_mode,
            )
            explicit.fromfile(str(lagr_path))
            kw_txt = ", ".join(
                f"{k}={labels.get(k)!r}"
                for k in ("mass_fraction", "add_star_type", "add_mass_range", "calc_energy", "external_mode")
                if labels.get(k) is not None
            )
            explicit_cmd = f"petar.LagrangianMultiple({kw_txt}).fromfile('{disp}')"
            agrees = int(explicit.size) == int(auto.size) and np.allclose(np.asarray(auto.time, dtype=float), np.asarray(explicit.time, dtype=float))
            check(f"{vname}: explicit read agrees", agrees, f"records={int(explicit.size)}", command=explicit_cmd)

            # meta= explicit injection: the documented escape hatch for files
            # copied away from their run directory must behave like the auto walk
            if meta_path.exists():
                try:
                    with warnings.catch_warnings(record=True) as rec:
                        warnings.simplefilter("always")
                        injected = petar.LagrangianMultiple(str(lagr_path), meta=str(meta_path))
                    guidance_inj = [str(w.message) for w in rec if "sidecar" in str(w.message) or "pass them explicitly" in str(w.message)]
                    same = int(injected.size) == int(auto.size) and np.allclose(np.asarray(injected.time, dtype=float), np.asarray(auto.time, dtype=float))
                    inj_cmd = f"petar.LagrangianMultiple('{disp}', meta='{meta_disp}')"
                    check(f"{vname}: meta= injection read", not guidance_inj and same, "; ".join(guidance_inj) if guidance_inj else f"{int(injected.size)} records", command=inj_cmd)
                except Exception as exc:  # noqa: BLE001 - report, do not mask
                    check(f"{vname}: meta= injection read", False, str(exc))
        except Exception as exc:  # noqa: BLE001 - report, do not mask
            check(f"{vname}: auto read", False, str(exc))
            continue

        # 3) legacy probe: the bare table without its sidecar
        with tempfile.TemporaryDirectory(dir=str(lagr_path.parent)) as td:
            legacy = Path(td) / f"{vprefix}.legacy.lagr"
            shutil.copy2(lagr_path, legacy)
            exp_fracs = expected.get("mass_fraction", LAGR_DEFAULT_MASS_FRACTION)
            shape_differs = (
                len(exp_fracs) != len(LAGR_DEFAULT_MASS_FRACTION)
                or bool(expected.get("add_star_type"))
                or bool(expected.get("add_mass_range"))
            )

            captured: List[str] = []
            legacy_obj = None
            raised = ""
            try:
                with warnings.catch_warnings(record=True) as rec:
                    warnings.simplefilter("always")
                    legacy_obj = petar.LagrangianMultiple(str(legacy))
                captured = [str(w.message) for w in rec]
            except Exception as exc:  # noqa: BLE001 - any failure must fail the case, not the suite
                raised = f"{type(exc).__name__}: {exc}"
                captured = [str(w.message) for w in rec]
            guided = [m for m in captured if "no sidecar found" in m]
            legacy_cmd = f"petar.LagrangianMultiple('{legacy.name}')  # bare copy, no sidecar"
            if shape_differs:
                check(f"{vname}: legacy read raises (shape differs)", bool(raised), raised or "read succeeded without the sidecar", command=legacy_cmd)
                check(f"{vname}: legacy raise carries guidance", bool(guided), "; ".join(guided) or "no no-sidecar warning", command=legacy_cmd)
            else:
                check(f"{vname}: legacy read succeeds (same shape)", not raised, raised, command=legacy_cmd)
                check(f"{vname}: legacy read warns mislabel", bool(guided), "; ".join(guided) or "no no-sidecar warning", command=legacy_cmd)
                if legacy_obj is not None:
                    legacy_labels = dict(getattr(legacy_obj, "initargs", {}) or {})
                    raw = legacy_labels.get("mass_fraction")
                    mis = np.asarray(raw, dtype=float) if raw is not None else np.array([])
                    check(
                        f"{vname}: legacy labels fall back to defaults",
                        mis.size == len(LAGR_DEFAULT_MASS_FRACTION) and np.allclose(mis, LAGR_DEFAULT_MASS_FRACTION),
                        f"got={mis}",
                        command="legacy.initargs['mass_fraction']",
                    )
                    # meta=False is the documented explicit opt-out: structure /
                    # explicit kwargs only, no guidance warning
                    try:
                        with warnings.catch_warnings(record=True) as rec2:
                            warnings.simplefilter("always")
                            optout = petar.LagrangianMultiple(str(legacy), meta=False)
                        guided2 = [str(w.message) for w in rec2 if "no sidecar found" in str(w.message)]
                        optout_cmd = f"petar.LagrangianMultiple('{legacy.name}', meta=False)"
                        check(f"{vname}: meta=False opt-out is silent", not guided2 and int(optout.size) == n_snapshots, "; ".join(guided2) or f"{int(optout.size)} records", command=optout_cmd)
                    except Exception as exc:  # noqa: BLE001 - report, do not mask
                        check(f"{vname}: meta=False opt-out is silent", False, str(exc))

    return step


def run_format_matrix_checks(
    repo_root: Path,
    case: Dict[str, Any],
    case_dir: Path,
    env: Dict[str, str],
    log_path: Path,
    process_args: List[str],
    snap_list_name: str,
    timeout_s: Optional[float] = None,
) -> Dict[str, Any]:
    """Output-format coverage matrix: solver snapshots in BINARY and ASCII,
    petar.data.process products in npy / ascii / binary.

    In one scratch copy of the case snapshots: (1) a data.process run with
    -s binary -o npy (lagr/core .npy + per-snapshot .single.npy products,
    all read back through the concise constructors), (2) the snapshots are
    converted to ASCII with the family petar.<family>.format.transfer tool
    (solver ASCII snapshot format) and read back, (3) a data.process run
    with -s ascii -o ascii, its ASCII lagr/core read back, and the lagr
    time series cross-checked against the npy run. Entries carry commands
    like the lagr-variant step and render in the same HTML table family.
    """
    import tempfile

    step: Dict[str, Any] = {
        "name": "format.matrix",
        "code": 0,
        "status": "pass",
        "reason": "",
        "command": f"petar.data.process -s binary -o npy / -s ascii -o ascii (scratch) ... {snap_list_name}",
        "cwd": str(case_dir),
        "details": {"variants": []},
    }
    to = timeout_s if timeout_s is not None and timeout_s > 0 else _step_timeout({})
    tools_path = repo_root / "tools"
    if str(tools_path) not in sys.path:
        sys.path.insert(0, str(tools_path))
    import analysis as petar  # type: ignore
    import numpy as np  # type: ignore

    def check(name: str, ok: bool, detail: str = "", command: str = "") -> None:
        entry = {"name": name, "ok": bool(ok), "detail": detail}
        if command:
            entry["command"] = command
        step["details"]["variants"].append(entry)
        if not ok:
            step["code"] = 1
            step["status"] = "fail"
            step["reason"] = step["reason"] or f"format matrix check failed: {name}"

    # per-file x per-format grid for the HTML report (binary / ascii / npy)
    matrix: Dict[str, Dict[str, Dict[str, Any]]] = {"snapshot": {}, "lagr": {}, "core": {}}

    def cell(file_key: str, fmt: str, ok: bool, detail: str, command: str = "") -> None:
        matrix[file_key][fmt] = {"ok": bool(ok), "detail": detail, "command": command}

    # scratch: snap list + listed snapshots + solver sidecar
    scratch = case_dir / "fmt_scratch"
    if scratch.exists():
        shutil.rmtree(scratch)
    scratch.mkdir(parents=True)
    snap_list = case_dir / snap_list_name
    shutil.copy2(snap_list, scratch / snap_list_name)
    snap_files = [line.strip() for line in snap_list.read_text(encoding="utf-8").splitlines() if line.strip()]
    for name in snap_files:
        src = case_dir / name
        if src.is_file():
            shutil.copy2(src, scratch / name)
    solver_meta = case_dir / (snap_list_name.split(".")[0] + ".meta")
    if solver_meta.is_file():
        shutil.copy2(solver_meta, scratch / solver_meta.name)
    n_snapshots = len(snap_files)
    first_snap = snap_files[0]
    last_snap = snap_files[-1]

    # ---- 0) binary baseline: solver binary snapshots + the main run's
    # binary lagr/core products, read through the concise constructors
    try:
        bin_obj = petar.Particle(str(case_dir / first_snap))
        cell("snapshot", "binary", True, f"{int(bin_obj.size)} particles", command=f"petar.Particle('{first_snap}')")
    except Exception as exc:  # noqa: BLE001
        bin_obj = None
        cell("snapshot", "binary", False, str(exc)[:90], command=f"petar.Particle('{first_snap}')")
    for path, ctor, label in ((case_dir / "data.lagr", petar.LagrangianMultiple, "data.lagr"), (case_dir / "data.core", petar.Core, "data.core")):
        if not path.exists():
            continue
        try:
            obj = ctor(str(path))
            cell("lagr" if "lagr" in label else "core", "binary", int(obj.size) == n_snapshots, f"{int(obj.size)} records", command=f"{ctor.__name__}('{label}')")
        except Exception as exc:  # noqa: BLE001
            cell("lagr" if "lagr" in label else "core", "binary", False, str(exc)[:90], command=f"{ctor.__name__}('{label}')")

    # ---- 1) data.process: -s binary -o npy (binary snapshots, npy products)
    cmd_npy = ["petar.data.process", "--no-auto-resume", "-p", "fmtnpy", "-s", "binary", "-o", "npy", *process_args, snap_list_name]
    rc, warns = run_cmd(cmd_npy, cwd=scratch, env=env, log_path=log_path, timeout_s=to)
    if warns:
        step.setdefault("warnings", []).extend(warns)
        step.setdefault("_warn_steps", []).append({"step": "format.matrix:npy", "messages": warns})
    check("npy: petar.data.process exit", rc == 0, f"rc={rc}", command="$ " + " ".join(cmd_npy))
    lagr_npy_time = None
    if rc == 0:
        lagr_npy = scratch / "fmtnpy.lagr.npy"
        core_npy = scratch / "fmtnpy.core.npy"
        single_npy = scratch / f"{last_snap}.single.npy"
        for path, ctor, label in ((lagr_npy, petar.LagrangianMultiple, "fmtnpy.lagr.npy"), (core_npy, petar.Core, "fmtnpy.core.npy")):
            if not path.exists():
                check(f"npy: {label} exists", False, f"{label} missing")
                continue
            try:
                obj = ctor(str(path))
                check(f"npy: read {label}", True, f"{int(obj.size)} records", command=f"{ctor.__name__}('{label}')")
                cell("lagr" if "lagr" in label else "core", "npy", int(obj.size) == n_snapshots, f"{int(obj.size)} records", command=f"{ctor.__name__}('{label}')")
                if "lagr" in label:
                    lagr_npy_time = np.asarray(obj.time, dtype=float)
            except Exception as exc:  # noqa: BLE001
                check(f"npy: read {label}", False, str(exc)[:90], command=f"{ctor.__name__}('{label}')")
                cell("lagr" if "lagr" in label else "core", "npy", False, str(exc)[:90], command=f"{ctor.__name__}('{label}')")
        if single_npy.exists():
            try:
                p = petar.Particle(str(single_npy))
                check("npy: read snapshot .single.npy", True, f"{int(p.size)} particles", command=f"petar.Particle('{last_snap}.single.npy')")
                cell("snapshot", "npy", True, f"{int(p.size)} particles", command=f"petar.Particle('{last_snap}.single.npy')")
            except Exception as exc:  # noqa: BLE001
                check("npy: read snapshot .single.npy", False, str(exc)[:90], command=f"petar.Particle('{last_snap}.single.npy')")
                cell("snapshot", "npy", False, str(exc)[:90], command=f"petar.Particle('{last_snap}.single.npy')")
        else:
            check("npy: read snapshot .single.npy", False, f"{single_npy.name} missing")
            cell("snapshot", "npy", False, f"{single_npy.name} missing")

    # ---- 2) convert snapshots to ASCII (solver ASCII format) + read back
    ft_tool = Path(active_petar_realpath()).parent / (Path(active_petar_realpath()).name + ".format.transfer")
    if not ft_tool.exists():
        check("ascii: family format.transfer tool", False, f"not found: {ft_tool.name}")
    else:
        cmd_ft = [str(ft_tool), "-b", "-r", snap_list_name]
        rc, warns = run_cmd(cmd_ft, cwd=scratch, env=env, log_path=log_path, timeout_s=to)
        if warns:
            step.setdefault("warnings", []).extend(warns)
        check("ascii: format.transfer exit", rc == 0, f"rc={rc}", command="$ " + " ".join(cmd_ft))
        if rc == 0:
            ascii_first = scratch / first_snap
            sniffed = petar.detectAsciiFile(str(ascii_first)) if ascii_first.exists() else None
            check("ascii: snapshot format honored", sniffed is True, f"sniffed={sniffed}")
            try:
                ascii_obj = petar.Particle(str(ascii_first))
                ok = sniffed is True and bin_obj is not None and int(ascii_obj.size) == int(bin_obj.size)
                check("ascii: snapshot read (Particle(path))", ok, f"ascii={int(ascii_obj.size)} binary={int(bin_obj.size) if bin_obj is not None else '?'} particles", command=f"petar.Particle('{first_snap}')")
                cell("snapshot", "ascii", ok, f"{int(ascii_obj.size)} particles", command=f"petar.Particle('{first_snap}')")
            except Exception as exc:  # noqa: BLE001
                check("ascii: snapshot read (Particle(path))", False, str(exc)[:90], command=f"petar.Particle('{first_snap}')")
                cell("snapshot", "ascii", False, str(exc)[:90], command=f"petar.Particle('{first_snap}')")

        # ---- 3) data.process: -s ascii -o ascii + read back + cross-format
        cmd_ascii = ["petar.data.process", "--no-auto-resume", "-p", "fmtascii", "-s", "ascii", "-o", "ascii", *process_args, snap_list_name]
        rc, warns = run_cmd(cmd_ascii, cwd=scratch, env=env, log_path=log_path, timeout_s=to)
        if warns:
            step.setdefault("warnings", []).extend(warns)
            step.setdefault("_warn_steps", []).append({"step": "format.matrix:ascii", "messages": warns})
        check("ascii: petar.data.process exit", rc == 0, f"rc={rc}", command="$ " + " ".join(cmd_ascii))
        if rc == 0:
            lagr_ascii = scratch / "fmtascii.lagr"
            core_ascii = scratch / "fmtascii.core"
            lagr_ascii_time = None
            for path, ctor, label in ((lagr_ascii, petar.LagrangianMultiple, "fmtascii.lagr"), (core_ascii, petar.Core, "fmtascii.core")):
                if not path.exists():
                    check(f"ascii: {label} exists", False, f"{label} missing")
                    continue
                sniffed = petar.detectAsciiFile(str(path))
                try:
                    obj = ctor(str(path))
                    ok = sniffed is True and int(obj.size) == n_snapshots
                    check(f"ascii: read {label}", ok, f"ascii={sniffed}, {int(obj.size)} records", command=f"{ctor.__name__}('{label}')")
                    cell("lagr" if "lagr" in label else "core", "ascii", ok, f"{int(obj.size)} records", command=f"{ctor.__name__}('{label}')")
                    if "lagr" in label:
                        lagr_ascii_time = np.asarray(obj.time, dtype=float)
                except Exception as exc:  # noqa: BLE001
                    check(f"ascii: read {label}", False, str(exc)[:90], command=f"{ctor.__name__}('{label}')")
                    cell("lagr" if "lagr" in label else "core", "ascii", False, str(exc)[:90], command=f"{ctor.__name__}('{label}')")
            if lagr_npy_time is not None and lagr_ascii_time is not None:
                same = lagr_npy_time.size == lagr_ascii_time.size and np.allclose(lagr_npy_time, lagr_ascii_time)
                check("cross-format: lagr time series agrees (npy vs ascii)", same, f"{lagr_npy_time.size} vs {lagr_ascii_time.size} records")

    step["details"]["matrix"] = matrix
    return step


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


def append_packed_container_read_check(result: Dict[str, Any], sel_issues: List[str], sel_details: Dict[str, Any], idx_path: Path) -> None:
    """Mirror the packed-mode.select step into the read-checks table as a
    'packed:container-read' row, so the container read + selection forms
    appear next to the packed-extract rows (readback of the extracted
    standalone files) instead of only in the steps table."""
    result.setdefault("read_checks", {"checks": [], "failed_checks": [], "warnings": [], "mismatch_warnings": []})
    idx_name = idx_path.name
    idx_stem = idx_name[:-len(".pack.idx")] if idx_name.endswith(".pack.idx") else idx_name
    entry = {
        "label": "packed:container-read",
        "file": str(idx_path),
        "file_note": "index of the packed container; snapshots are seeked out of the %s.pack.<id> VOLUMES via this index (no copy)" % idx_stem,
        "command": "pack = petar.PackedSnapshots('%s'); pack.readSnapshot(time=/nfile=/row) / pack[i]; petar.pack extract %s --time/--nfile/--last/--all" % (idx_stem, idx_stem),
        "ok": not sel_issues,
        "error": "; ".join(sel_issues) if sel_issues else "",
        "warnings": [],
        "details": sel_details,
    }
    result["read_checks"]["checks"].append(entry)
    if sel_issues:
        result["read_checks"]["failed_checks"].append(entry)


def packed_selection_checks(pack) -> tuple:
    """Exercise and document every snapshot-selection form of a
    ``petar.PackedSnapshots`` container (discover times / by time / by
    snapshot id / by row / range mask / CLI materialization); returns
    (issues, details) for a ``packed-mode.select`` step. Every form is
    actually executed and its exact call form recorded for the report."""
    import numpy as np

    issues: List[str] = []
    details: Dict[str, Any] = {}
    idx = pack.index
    try:
        times = np.asarray(idx['time'], dtype=float)
        shown = [round(float(t), 6) for t in times[:8]]
        details['discover times'] = 'pack.index["time"] -> %s%s' % (shown, ' ...' if len(times) > 8 else '')
        t_mid = float(times[1]) if len(times) > 1 else float(times[0])
        s_time = pack.readSnapshot(time=t_mid)
        details['by time'] = 'pack.readSnapshot(time=%g) -> first snapshot at/after t: nfile %d, time %g' % (
            t_mid, int(s_time.header.file_id), float(s_time.header.time))
        expect_mid = int(idx['nfile'][1]) if len(times) > 1 else int(idx['nfile'][0])
        if int(s_time.header.file_id) != expect_mid:
            issues.append('time selection landed on nfile %d, expected %d' % (s_time.header.file_id, expect_mid))
        s_nfile = pack.readSnapshot(nfile=int(idx['nfile'][-1]))
        details['by snapshot id'] = 'pack.readSnapshot(nfile=%d) -> nfile %d, time %g' % (
            int(idx['nfile'][-1]), int(s_nfile.header.file_id), float(s_nfile.header.time))
        s_last = pack[-1]
        details['by row'] = 'pack[-1] / pack.readSnapshot(-1) -> nfile %d, time %g' % (
            int(s_last.header.file_id), float(s_last.header.time))
        if int(s_nfile.header.file_id) != int(s_last.header.file_id) or abs(float(s_nfile.header.time) - float(s_last.header.time)) > 1e-10:
            issues.append('nfile selection differs from last-row selection')
        details['range query'] = 'rows = np.where((pack.index["time"] >= a) & (pack.index["time"] <= b)); [pack.readSnapshot(row=int(r)) for r in rows]'
        details['materialize (CLI)'] = 'petar.pack extract <prefix> --time <t> | --nfile <n> | --last | --all --snap-lst'
    except Exception as exc:  # noqa: BLE001
        issues.append('selection-method error: %r' % (exc,))
        details['error'] = repr(exc)
    return issues, details


def run_packed_mode_coverage(result: Dict[str, Any], case: Dict[str, Any], case_dir: Path, repo_root: Path, env: Dict[str, str], log_path: Path, dt_out: str, run_args: List[str], require: str, timeout_s: Optional[float] = None) -> None:
    """Packed-output coverage for every solver mode: a short scratch run with
    ``-w 4`` under the ``packcov`` prefix (the main case run keeps the default
    ``-w 1`` single-file output and its readback untouched), index validation,
    ``petar.pack`` extraction, and the standard paired-read readback on the
    extracted snapshots - so every interrupt/external column layout also goes
    through the packed path. Appears in the HTML report as ``packed-mode.*``
    steps (container seek read + selection forms) and
    ``packed-extract:``-prefixed read-check rows (readback of the extracted
    standalone files). The dedicated std-packed
    case covers the deep restart/trim flow and skips this stage."""
    scratch_prefix = "packcov"
    # scratch artifacts from an earlier run must not be re-read
    for stale in case_dir.glob(f"{scratch_prefix}*"):
        if stale.is_file():
            stale.unlink()

    dt_out_f = float(dt_out)
    t_packed = repr(2 * dt_out_f)
    expected_records = 3  # t = 0, dt_out, 2*dt_out

    run_cmdline = ["petar", "-u", "1", "-t", t_packed, "-o", dt_out, "-f", scratch_prefix, *run_args, "-w", "4", "input"]
    rc, warns = run_cmd(run_cmdline, cwd=case_dir, env=env, log_path=log_path, timeout_s=timeout_s)
    add_step(result, "packed-mode.run (-w 4)", rc, run_cmdline, case_dir, warnings_list=warns,
             details={"t_end": t_packed, "dt_out": dt_out, "output_prefix": scratch_prefix, "expected_records": expected_records})
    if rc != 0:
        result["status"] = "fail"
        return

    list_cmd = ["petar.pack", "list", scratch_prefix]
    rc, warns = run_cmd(list_cmd, cwd=case_dir, env=env, log_path=log_path, timeout_s=timeout_s)
    add_step(result, "packed-mode.index (petar.pack list)", rc, list_cmd, case_dir, warnings_list=warns)
    if rc != 0:
        result["status"] = "fail"
        return

    extract_cmd = ["petar.pack", "extract", scratch_prefix, "--all", "--snap-lst"]
    rc, warns = run_cmd(extract_cmd, cwd=case_dir, env=env, log_path=log_path, timeout_s=timeout_s)
    add_step(result, "packed-mode.extract", rc, extract_cmd, case_dir, warnings_list=warns)
    if rc != 0:
        result["status"] = "fail"
        return

    # In-process validation: index shape, per-volume offset contiguity,
    # record-size consistency, and the direct-seek PackedSnapshots read
    # against the extracted standalone file.
    tools_path = repo_root / "tools"
    if str(tools_path) not in sys.path:
        sys.path.insert(0, str(tools_path))
    import analysis as petar  # type: ignore
    import numpy as np  # type: ignore

    issues: List[str] = []
    details: Dict[str, Any] = {}
    try:
        pack = petar.PackedSnapshots(str(case_dir / scratch_prefix))
        idx = pack.index
        details["records"] = int(len(idx))
        if len(idx) != expected_records:
            issues.append(f"record count {len(idx)} != expected {expected_records}")
        if len(idx) and (int(idx["nfile"][0]) != 0 or int(idx["nfile"][-1]) != len(idx) - 1):
            issues.append("nfile sequence does not start at 0 / is not contiguous")
        times = np.asarray(idx["time"], dtype=float)
        if len(times) > 1 and not np.all(np.diff(times) > 0):
            issues.append("times are not strictly ascending")
        itemsize = np.dtype(petar.Particle(**pack.kwargs).collectDtype()).itemsize
        cursor: Dict[int, int] = {}
        for rec in idx:
            vol = int(rec["volume"])
            off = int(rec["offset"])
            if vol in cursor and off != cursor[vol]:
                issues.append(f"volume {vol}: offset {off} != contiguous {cursor[vol]}")
            cursor[vol] = off + int(rec["nbytes"])
            if int(rec["nbytes"]) != pack.header_size + int(rec["n_body"]) * itemsize:
                issues.append(f"nfile {int(rec['nfile'])}: nbytes {int(rec['nbytes'])} != header {pack.header_size} + N*{itemsize}")
        snap_pack = pack.readSnapshot(row=-1)
        p_pack, h_pack = snap_pack.particle, snap_pack.header
        last_name = f"{scratch_prefix}.{int(idx['nfile'][-1])}"
        snap_ref = petar.Snapshot(str(case_dir / last_name))
        p_ref, h_ref = snap_ref.particle, snap_ref.header
        same = int(p_pack.size) == int(p_ref.size)
        for entry in p_pack.keys:
            m_pack, m_ref = getattr(p_pack, entry[0]), getattr(p_ref, entry[0])
            if hasattr(m_pack, "keys") or hasattr(m_ref, "keys"):
                continue  # nested containers: the byte-identical record already covers them
            same = same and np.array_equal(np.asarray(m_pack), np.asarray(m_ref))
        same = same and (h_pack.file_id == h_ref.file_id) and (h_pack.n == h_ref.n) \
            and (abs(h_pack.time - h_ref.time) <= 1e-10 * max(1.0, abs(h_ref.time)))
        if not same:
            issues.append("direct pack read differs from the extracted snapshot")
        details["last_snapshot"] = last_name
        details["volumes"] = sorted(cursor.keys())
    except Exception as exc:  # noqa: BLE001 - report as a failed step
        issues.append(f"python packed validation error: {exc!r}")
        details["error"] = repr(exc)
    # Container read row in the READ-CHECKS table (NOT a step: index
    # validation, the direct-seek read and the selection forms are Python
    # data-analysis reads, not executed shell commands - they belong next
    # to the packed-extract rows, not in the command steps). One row carries
    # the index validation, the direct-read comparison and every selection
    # form with its exact call form.
    sel_issues, sel_details = packed_selection_checks(pack)
    container_details = dict(sel_details)
    if 'error' in details:
        container_details['index validation'] = 'failed: %s' % details['error']
    else:
        container_details['index validation'] = 'ok: %d records, volumes %s, offsets contiguous, record sizes = header + N*itemsize' % (
            details.get('records', 0), details.get('volumes', []))
        container_details['direct read'] = 'pack.readSnapshot(-1) == extracted %s (member-wise identical)' % details.get('last_snapshot', '?')
    append_packed_container_read_check(result, issues + sel_issues, container_details, case_dir / (scratch_prefix + ".pack.idx"))
    if issues or sel_issues:
        result["status"] = "fail"
        return

    # Standard paired-read readback on the extracted snapshots (the same
    # machinery as the main chain; labels prefixed for the HTML report).
    # Packing changes SNAPSHOT storage only - status/escaper/event/profile
    # streams are written identically to the unpacked mode and are already
    # covered by the case's main readback - so only snapshot-family and
    # schema rows are kept in the packed report. (The dedicated std-packed
    # case runs entirely in packed mode and keeps its full readback.)
    try:
        scratch_checks = run_python_output_read_checks(
            repo_root,
            case_dir,
            require,
            scratch_prefix,
            interrupt_mode_override=case.get("interrupt_mode", ""),
        )
    except RuntimeError as exc:
        # environment-level failure: record it as a failed read-check row
        # (Python analysis, not a shell command), then fail the case
        result.setdefault("read_checks", {"checks": [], "failed_checks": [], "warnings": [], "mismatch_warnings": []})
        err_entry = {
            "label": "packed-extract:readback",
            "file": str(case_dir),
            "command": "paired reads of extracted packcov.* standalone snapshots",
            "ok": False,
            "error": str(exc),
            "warnings": [],
            "details": {},
        }
        result["read_checks"]["checks"].append(err_entry)
        result["read_checks"]["failed_checks"].append(err_entry)
        result["status"] = "fail"
        return
    packed_rows = [
        entry for entry in scratch_checks.get("checks", [])
        if str(entry.get("label", "")).startswith(("schema:", "snapshot"))
    ]
    for entry in packed_rows:
        # 'packed-extract:' marks readback of the EXTRACTED standalone files
        # (packcov.N) - the container read itself is the
        # packed:container-read row above (seek read + selection forms).
        # Annotate each extracted snapshot with its container origin
        # (record nfile/time/volume) so the file cell explains what it IS.
        entry["label"] = "packed-extract:" + str(entry.get("label", ""))
        base = Path(str(entry.get("file", ""))).name
        m = re.match(r"^%s\.(\d+)$" % re.escape(scratch_prefix), base)
        if m:
            n = int(m.group(1))
            row = np.where(idx["nfile"] == n)[0]
            if len(row):
                rec = idx[int(row[0])]
                entry["file_note"] = "standalone snapshot extracted from the packed container (petar.pack extract --all) = record nfile %d (t=%g) of %s.pack.%04d" % (
                    n, float(rec["time"]), scratch_prefix, int(rec["volume"]))
    merged = result.setdefault("read_checks", {"checks": [], "failed_checks": [], "warnings": [], "mismatch_warnings": []})
    merged.setdefault("checks", []).extend(packed_rows)
    merged.setdefault("failed_checks", []).extend([e for e in packed_rows if not e.get("ok")])
    merged.setdefault("warnings", []).extend(scratch_checks.get("warnings", []))
    n_failed = sum(1 for e in packed_rows if not e.get("ok"))
    if n_failed:
        result["status"] = "fail"


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
    # Boundary conventions at t == cut (2026-10-28): group streams and the
    # sse/bse-family rules below differ — group records are state prints
    # from BOTH sides of the boundary (phase-A substep end + restart
    # creation record at t=restart), so petar.data.clear trims them
    # STRICTLY; sse/bse events at exactly t==cut belong to the hard step
    # containing them and are baked into the post-event snapshot state, so
    # the restart does NOT re-emit them and the trim stays INCLUSIVE (a
    # strict trim would silently lose such events). If a boundary duplicate
    # ever fires the restart duplicate check in these streams, flip this
    # family to strict — the check is the arbiter.
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

        def add_reader(p: Path, ctor, strict_boundary: bool = False, **kw) -> None:
            if _petar is None:
                stats[p.name] = {"kind": "unreadable"}
                return
            try:
                obj = ctor(**kw)
                obj.fromfile(str(p))
                times = [float(x) for x in _np.asarray(obj.time)]
                # group streams: strict < (the restart re-emits the boundary
                # record for snapshot groups, so petar.data.clear trims it;
                # see the group branch there)
                if strict_boundary:
                    stats[p.name] = {
                        "kind": "reader",
                        "kept": sum(1 for x in times if x < t_cut - 1e-9),
                        "beyond": sum(1 for x in times if x >= t_cut - 1e-9),
                    }
                else:
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
                add_reader(p, _petar.GroupInfo, strict_boundary=True, N=int(m.group(1)), interrupt_mode=mode_i or "none", external_mode=ext_e)
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

        # Packed-output restart hook 1/3: materialize the phase-1 packed
        # snapshots as standalone files so the shared restart machinery
        # (snap.lst, data.clear, snapshot picking) works unchanged.
        if case.get("packed_output"):
            pack_extract_cmd = ["petar.pack", "extract", output_prefix, "--all", "--snap-lst"]
            rc, warns = run_cmd(pack_extract_cmd, cwd=case_dir, env=env, log_path=log_path, timeout_s=to)
            add_step(result, "petar.pack.extract", rc, pack_extract_cmd, case_dir, warnings_list=warns)
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

        # Packed-output restart hook 2/3: data.clear trims the extracted
        # copies; the container itself needs the matching petar.pack trim so
        # the restart appends without stale trailing records.
        if case.get("packed_output"):
            pack_trim_cmd = ["petar.pack", "trim", output_prefix, "-t", str(clear_t)]
            rc, warns = run_cmd(pack_trim_cmd, cwd=case_dir, env=env, log_path=log_path, timeout_s=to)
            add_step(result, "petar.pack.trim", rc, pack_trim_cmd, case_dir, warnings_list=warns, details={"clear_t": clear_t})
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

        # Packed-output restart hook 3/3: the restarted run appends new packed
        # records only (snap.lst is not touched in packed mode); re-extract so
        # the progress/duplicate checks see the full snapshot set.
        if case.get("packed_output"):
            pack_extract_cmd = ["petar.pack", "extract", output_prefix, "--all", "--snap-lst"]
            rc, warns = run_cmd(pack_extract_cmd, cwd=case_dir, env=env, log_path=log_path, timeout_s=to)
            add_step(result, "petar.pack.extract:restart", rc, pack_extract_cmd, case_dir, warnings_list=warns)
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

    # Packed-output case hook (matrix key 'packed_output', solver -w 4):
    # materialize the packed snapshots as standalone files (byte-identical to
    # the -w 1 layout) so the shared downstream pipeline (data.process /
    # movie / readback checks) exercises them, then verify the container
    # itself: the Python direct-seek read against the extracted copy.
    if case.get("packed_output"):
        tools_path = repo_root / "tools"
        if str(tools_path) not in sys.path:
            sys.path.insert(0, str(tools_path))
        import analysis as petar  # type: ignore
        import numpy as np  # type: ignore
        pack_extract_cmd = ["petar.pack", "extract", output_prefix, "--all", "--snap-lst"]
        rc, warns = run_cmd(pack_extract_cmd, cwd=case_dir, env=env, log_path=log_path, timeout_s=to)
        add_step(result, "petar.pack.extract", rc, pack_extract_cmd, case_dir, warnings_list=warns)
        if rc != 0:
            result["status"] = "fail"
            return result
        rb_cmd = ["petar.pack", "rebuild", output_prefix]
        rc, warns = run_cmd(rb_cmd, cwd=case_dir, env=env, log_path=log_path, timeout_s=to)
        add_step(result, "petar.pack.rebuild", rc, rb_cmd, case_dir, warnings_list=warns)
        if rc != 0:
            result["status"] = "fail"
            return result
        try:
            pack_obj = petar.PackedSnapshots(str(case_dir / output_prefix))
            snap_pack = pack_obj.readSnapshot(row=-1)
            last_name = f"{output_prefix}.{int(pack_obj.index[-1]['nfile'])}"
            snap_ref = petar.Snapshot(str(case_dir / last_name))
            p_pack, h_pack = snap_pack.particle, snap_pack.header
            p_ref, h_ref = snap_ref.particle, snap_ref.header
            same = int(p_pack.size) == int(p_ref.size)
            for entry in p_pack.keys:
                m_pack, m_ref = getattr(p_pack, entry[0]), getattr(p_ref, entry[0])
                if hasattr(m_pack, "keys") or hasattr(m_ref, "keys"):
                    continue  # nested containers: the byte-identical record already covers them
                same = same and np.array_equal(np.asarray(m_pack), np.asarray(m_ref))
            same = same and (h_pack.file_id == h_ref.file_id) and (h_pack.n == h_ref.n) \
                and (abs(h_pack.time - h_ref.time) <= 1e-10 * max(1.0, abs(h_ref.time)))
            n_idx = len(pack_obj.index)
        except Exception as exc:  # noqa: BLE001 - report as a failed container-read row
            result["_pending_packed_container_read"] = ([repr(exc)], {"error": repr(exc)}, case_dir / (output_prefix + ".pack.idx"))
            result["status"] = "fail"
            return result

        # container read + selection forms go to the READ-CHECKS table (not
        # steps: Python analysis reads, not shell commands); stashed because
        # this hook runs before the main read_checks assignment below
        sel_issues, sel_details = packed_selection_checks(pack_obj)
        sel_details["direct read"] = "pack.readSnapshot(-1) == extracted %s (member-wise %s)" % (last_name, "identical" if same else "DIFFERS")
        sel_details["index validation"] = "ok: %d records (index consistency verified by petar.pack rebuild)" % n_idx
        result["_pending_packed_container_read"] = (
            sel_issues + ([] if same else ["direct pack read differs from the extracted snapshot"]),
            sel_details,
            case_dir / (output_prefix + ".pack.idx"),
        )
        if not same or sel_issues:
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

    # 5a) Lagrangian variants: petar.data.process with non-default
    # --mass-fraction / --add-star-type / --add-mass-range under scratch
    # prefixes, plus the sidecar/legacy read matrix (see
    # run_lagr_variant_checks). Universal defaults come from the matrix
    # global 'lagr_variants'; case entries override by name or append.
    global_lagr_variants = matrix.get("global", {}).get("lagr_variants", [])
    effective_lagr_variants = list(global_lagr_variants)
    for case_variant in case.get("lagr_variants", []):
        for gi, gvi in enumerate(effective_lagr_variants):
            if gvi.get("name") == case_variant.get("name"):
                effective_lagr_variants[gi] = case_variant
                break
        else:
            effective_lagr_variants.append(case_variant)
    lagr_variant_step = run_lagr_variant_checks(repo_root, case, case_dir, env, log_path, process_args, snap_list.name, timeout_s=to, variants=effective_lagr_variants)
    if not lagr_variant_step["details"].get("skipped"):
        result["steps"].append(lagr_variant_step)
        for entry in lagr_variant_step.pop("_warn_steps", []):
            result["warnings"].append(entry)
        if lagr_variant_step["status"] != "pass":
            result["status"] = "fail"
            result["skip_reason"] = lagr_variant_step["reason"]
            return result

    # 5b) Output-format coverage matrix: solver snapshots in binary and
    # ASCII, petar.data.process products in npy / ascii / binary, read back
    # through the concise constructors (see run_format_matrix_checks).
    format_matrix_step = run_format_matrix_checks(repo_root, case, case_dir, env, log_path, process_args, snap_list.name, timeout_s=to)
    result["steps"].append(format_matrix_step)
    for entry in format_matrix_step.pop("_warn_steps", []):
        result["warnings"].append(entry)
    if format_matrix_step["status"] != "pass":
        result["status"] = "fail"
        result["skip_reason"] = format_matrix_step["reason"]
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

    try:
        read_checks = run_python_output_read_checks(
            repo_root,
            case_dir,
            require,
            output_prefix,
            interrupt_mode_override=case.get("interrupt_mode", ""),
        )
    except RuntimeError as exc:
        # Environment-level failure (e.g. a stale installed sdar package missing
        # petar.detectAsciiFile): fail THIS case loudly but keep the run going so
        # the report still lands, instead of aborting the whole matrix.
        read_checks = {
            "checks": [],
            "failed_checks": [{"label": "python_read_checks", "ok": False, "error": str(exc)}],
            "warnings": [],
            "mismatch_warnings": [],
        }
    result["read_checks"] = read_checks
    if read_checks.get("warnings"):
        result["warnings"].append({"step": "python_read_checks", "messages": [x["message"] for x in read_checks["warnings"]]})
    if read_checks.get("failed_checks"):
        result["status"] = "fail"
        if result.get("skip_reason"):
            result["skip_reason"] += "; python read checks failed"
        else:
            result["skip_reason"] = "python read checks failed"

    # the std-packed hook stashed its container-read row before this
    # assignment - append it to the final read-checks table now
    pending_container_read = result.pop("_pending_packed_container_read", None)
    if pending_container_read is not None:
        append_packed_container_read_check(result, *pending_container_read)

    # 10) Packed-output coverage for every solver mode (matrix-global): a
    # short -w 4 scratch run under the 'packcov' prefix + extraction + the
    # standard readback, reported as 'packed-mode.*' steps (container read +
    # selection forms) and 'packed-extract:' read-check rows (readback of
    # the extracted standalone files). The main -w 1 chain above is untouched.
    # The dedicated std-packed case covers the deep restart/trim flow and
    # skips this stage.
    if not case.get("packed_output"):
        run_packed_mode_coverage(result, case, case_dir, repo_root, env, log_path, dt_out, run_args, require, timeout_s=to)

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


def run_schema_selftest(repo_root: Path, out_root: Path) -> Dict[str, Any]:
    """One-time synthetic checks for the schema auto-detection layer.

    Covers the cases the per-case checks cannot (no rare-flag build exists in
    the matrix): same-tier ambiguity must raise, a correct rare-flag sidecar
    must be accepted even when the tier-0 structure sniff resolves elsewhere,
    and a wrong-family sidecar must be rejected.
    """
    tools_path = repo_root / "tools"
    if str(tools_path) not in sys.path:
        sys.path.insert(0, str(tools_path))
    import analysis as petar  # type: ignore
    import numpy as np  # type: ignore
    import tempfile

    result: Dict[str, Any] = {"status": "pass", "checks": []}

    def record(name: str, ok: bool, detail: str = "") -> None:
        result["checks"].append({"name": name, "ok": ok, "detail": detail})
        if not ok:
            result["status"] = "fail"

    def write_ascii_table(path: Path, ncol: int, rows: int = 5, ext: bool = False) -> None:
        with path.open("w", encoding="utf-8") as fh:
            fh.write(("0 %d 1.0 " % rows + " ".join(["1.0"] * 6) + "\n") if ext else ("0 %d 1.0\n" % rows))
            for _ in range(rows):
                fh.write(" ".join(["1.0"] * ncol) + "\n")

    with tempfile.TemporaryDirectory(dir=str(out_root)) as td:
        tdp = Path(td)
        # 1) same-tier collision (23 columns: none+sp_acc vs none+mpfrc) must raise
        f_amb = tdp / "data.0"
        write_ascii_table(f_amb, 23)
        try:
            petar.getDataKwargs(str(f_amb))
            record("same-tier-ambiguity-raises", False, "AmbiguousSchemaError not raised")
        except petar.AmbiguousSchemaError:
            record("same-tier-ambiguity-raises", True)
        # 2) correct tier-2 sidecar accepted although the tier-0 sniff resolves elsewhere
        #    (36 columns: dsm+mpfrc+sp_acc collides with bse/spin_3d)
        write_ascii_table(f_amb, 36)
        petar.writeMetaData(str(tdp / "data.meta"), interrupt_mode="dsm", external_mode="none",
                            use_mpfrc=True, collect_sp_acc=True)
        try:
            kw = petar.getDataKwargs(str(f_amb))
            ok = kw.get("interrupt_mode") == "dsm" and kw.get("use_mpfrc") and kw.get("collect_sp_acc")
            record("rare-flag-sidecar-accepted", ok, str(kw))
        except Exception as exc:  # noqa: BLE001
            record("rare-flag-sidecar-accepted", False, str(exc))
        # 3) wrong-family sidecar rejected
        petar.writeMetaData(str(tdp / "data.meta"), interrupt_mode="merger", external_mode="none")
        try:
            petar.getDataKwargs(str(f_amb))
            record("wrong-sidecar-rejected", False, "ValueError not raised")
        except ValueError:
            record("wrong-sidecar-rejected", True)
        # restore a CORRECT sidecar after check 3 deliberately corrupted
        # it: later checks validate clean paths (the sidecar itself is now
        # cross-checked, not just the merged kwargs)
        write_ascii_table(f_amb, 36)  # bse + spin_3d, no external
        petar.writeMetaData(str(tdp / "data.meta"), interrupt_mode="bse", external_mode="none", spin_3d=True)
        # 4) explicit-flag precedence in the CLI auto-fill helper
        kw = petar.autoFillKwargs(str(f_amb), {"interrupt_mode": "bse"}, verbose=False)
        record("autofill-explicit-wins", kw.get("interrupt_mode") == "bse", str(kw.get("interrupt_mode")))
        # 5) an explicit kwarg contradicting BOTH sidecar and structure is
        #    attributed to the explicit flag (e.g. --spin-1d on 3D-spin
        #    data), not misreported as a sidecar-from-another-run error
        write_ascii_table(f_amb, 36)  # bse + spin_3d, no external
        petar.writeMetaData(str(tdp / "data.meta"), interrupt_mode="bse", external_mode="none", spin_3d=True)
        try:
            petar.autoFillKwargs(str(f_amb), {"spin_3d": False}, verbose=False)
            record("explicit-wrong-kwarg-diagnosed", False, "no error raised")
        except ValueError as exc:  # noqa: BLE001
            msg = str(exc)
            record("explicit-wrong-kwarg-diagnosed", msg.startswith("explicitly given layout parameter") and "spin_3d=False" in msg, msg[:80])
        # 6) the same file without conflicting flags fills cleanly
        kw = petar.autoFillKwargs(str(f_amb), {}, verbose=False)
        record("explicit-clean-fill", kw.get("interrupt_mode") == "bse" and kw.get("spin_3d") is True, str(kw))

    return result


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
        report["schema_selftest"] = run_schema_selftest(repo_root, out_root)
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
                    "schema_selftest": report.get("schema_selftest"),
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
    selftest_failed = report.get("schema_selftest", {}).get("status") == "fail"
    return 1 if failed or build_failed or selftest_failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
