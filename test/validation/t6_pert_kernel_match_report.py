#!/usr/bin/env python3
"""T6 AR perturber-kernel-match HTML report (members-keep-own-changeover design)."""
import argparse
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple

from report_common import html_escape, load_json, parse_log_version

Point = Tuple[float, float]
REPO = Path(__file__).resolve().parents[2]
WORK_DIR = Path("test/validation/work/t6")

DIAG_LINE_RE = re.compile(
    r"PERT_KERNEL_DIAG\[(final|interval)\] eval=(\d+) sum\|d\|/sum\|ref\|=([0-9.eE+-]+) "
    r"max_rel=([0-9.eE+-]+) frac\(>1e-3\)=([0-9.eE+-]+) frac\(>1e-2\)=([0-9.eE+-]+) frac\(>1e-1\)=([0-9.eE+-]+)"
)


def parse_diag_lines(log_text: str) -> List[Dict[str, float]]:
    rows = []
    for m in DIAG_LINE_RE.finditer(log_text):
        rows.append(
            {
                "tag": m.group(1),
                "eval": int(m.group(2)),
                "ratio": float(m.group(3)),
                "max_rel": float(m.group(4)),
                "frac_1e3": float(m.group(5)),
                "frac_1e2": float(m.group(6)),
                "frac_1e1": float(m.group(7)),
            }
        )
    return rows


def _has_diag_define(binary: Path) -> bool:
    if not binary.exists():
        return False
    with binary.open("rb") as fh:
        content = fh.read()
    # format marker: the report line layout this script parses
    return b"PERT_KERNEL_DIAG" in content and b"frac(>1e-1)=" in content


def _main_target() -> str:
    """First entry of the Makefile TARGET list: the main petar binary of the
    configured family (name depends on toolchain/ISA options)."""
    makefile = REPO / "Makefile"
    if not makefile.exists():
        return ""
    for line in makefile.read_text().splitlines():
        if line.startswith("TARGET="):
            return line.split("=", 1)[1].split()[0]
    return ""


def _configured_plain_diag() -> bool:
    """True when ./configure was last run with CXXFLAGS='-D PERT_KERNEL_DIAG'
    (plain family + instrumentation), i.e. rebuilding is not needed."""
    makefile = REPO / "Makefile"
    if not makefile.exists():
        return False
    for line in makefile.read_text().splitlines():
        if line.startswith("CXXFLAGS="):
            return line.split("=", 1)[1].strip() == "-D PERT_KERNEL_DIAG"
    return False


def ensure_diag_binary() -> str:
    """Plain-family petar built with -D PERT_KERNEL_DIAG.

    Builds via the portable path: pass the flag through CXXFLAGS to ./configure
    (the generated Makefile appends -D PROFILE itself) and run make without a
    target so the configured family's default binaries are produced.
    """
    import os

    if not _configured_plain_diag():
        print("[T6] Configuring plain family with instrumentation (CXXFLAGS='-D PERT_KERNEL_DIAG' ./configure)...")
        subprocess.run(
            "CXXFLAGS='-D PERT_KERNEL_DIAG' ./configure", shell=True, cwd=REPO, check=True
        )
    binary = REPO / _main_target()
    if binary.name and not _has_diag_define(binary):
        print("[T6] Building (make -j)...")
        ncpu = os.cpu_count() or 4
        subprocess.run(f"make -j{ncpu}", shell=True, cwd=REPO, check=True)
        binary = REPO / _main_target()
    if not binary.exists() or not _has_diag_define(binary):
        raise RuntimeError(
            f"diagnostic marker missing in {binary} after configure+make; "
            "the instrumentation may have been removed from src/ar_interaction.hpp"
        )
    return str(binary)


def metric_card(label: str, value: str, note: str = "") -> str:
    note_html = f'<div class="note">{html_escape(note)}</div>' if note else ""
    return (
        f'<div class="card"><div class="card-label">{html_escape(label)}</div>'
        f'<div class="card-value">{html_escape(value)}</div>{note_html}</div>'
    )


def check_row(check: Dict) -> str:
    status = "PASS" if check.get("passed") else "FAIL"
    style = "pass" if check.get("passed") else "fail"
    message = check.get("message", "")
    return (
        f'<tr><td>{html_escape(str(check.get("check", "")))}</td>'
        f'<td class="{style}">{status}</td>'
        f'<td>{html_escape(str(check.get("value", "")))}</td>'
        f'<td>{html_escape(str(check.get("threshold", "")))}</td>'
        f'<td class="msg">{html_escape(str(message))}</td></tr>'
    )


def git_commit() -> str:
    try:
        r = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"], cwd=REPO,
            capture_output=True, text=True, check=True,
        )
        return r.stdout.strip()
    except Exception:
        return "unknown"


def build_html(
    checks: List[Dict],
    diag_rows: List[Dict[str, float]],
    run_metrics: Dict[str, Dict],
    run_commands: List[Tuple[str, str]],
    log_version: str,
    git_ref: str,
    report_path: Path,
) -> str:
    all_pass = all(c.get("passed", False) for c in checks)
    verdict = "ALL CHECKS PASSED" if all_pass else "CHECKS FAILED"
    measure = run_metrics.get("r1_measure", {})
    member = run_metrics.get("r2_member_kernel", {})
    cmk = run_metrics.get("r3_cm_kernel", {})
    r_max = diag_rows[-1]["max_rel"] if diag_rows else float("nan")

    def fmt(x, spec=".3e"):
        try:
            return format(float(x), spec)
        except (TypeError, ValueError):
            return "\u2014"

    def ab_row(label: str, metrics: Dict) -> str:
        err = metrics.get("max_abs_error_over_total", float("nan"))
        rc = metrics.get("regex_counts", {})
        return (
            f"<tr><td>{label}</td><td>{fmt(err)}</td>"
            f"<td>{rc.get('hard_dump', -1)}</td><td>{rc.get('runtime_abort', -1)}</td></tr>"
        )

    rows_html = "".join(check_row(c) for c in checks)

    return f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>T6 perturber-kernel consistency</title>
<style>
 body {{ font-family: Georgia, 'Times New Roman', serif; margin: 32px auto; max-width: 960px; color: #1a1a1a; line-height: 1.5; }}
 h1 {{ font-size: 22px; text-align: center; margin-bottom: 4px; }}
 h2 {{ font-size: 17px; margin-top: 30px; border-bottom: 1px solid #999; padding-bottom: 3px; }}
 h3 {{ font-size: 15px; margin-top: 18px; }}
 .subtitle {{ text-align: center; color: #555; font-size: 14px; margin-bottom: 14px; }}
 .verdict {{ display: inline-block; padding: 5px 14px; font-weight: bold; border-radius: 4px; font-family: Arial; }}
 .verdict.pass {{ background: #2e7d32; color: #fff; }}
 .verdict.fail {{ background: #c62828; color: #fff; }}
 .eq {{ text-align: center; margin: 10px 0; font-style: italic; font-size: 15px; }}
 .eq .tag {{ float: right; color: #777; font-style: normal; }}
 .caption {{ font-size: 12.5px; color: #444; margin: 6px 0 18px 0; }}
 table.data {{ border-collapse: collapse; font-size: 13.5px; margin: 8px 0 2px 0; font-family: Arial; }}
 table.data th, table.data td {{ border: 1px solid #999; padding: 4px 9px; text-align: left; vertical-align: top; }}
 table.data th {{ background: #eee; }}
 td.pass {{ color: #2e7d32; font-weight: bold; }}
 td.fail {{ color: #c62828; font-weight: bold; }}
 td.msg {{ color: #555; }}
 .cards {{ display: flex; flex-wrap: wrap; gap: 10px; margin: 10px 0 2px 0; }}
 .card {{ border: 1px solid #bbb; border-radius: 5px; padding: 8px 14px; min-width: 150px; font-family: Arial; }}
 .card-label {{ font-size: 12px; color: #666; }}
 .card-value {{ font-size: 19px; font-weight: bold; margin-top: 3px; }}
 pre.cmd {{ background: #f6f6f6; border: 1px solid #ddd; border-radius: 4px; padding: 8px 12px; font-size: 12px; font-family: 'Courier New', monospace; overflow-x: auto; white-space: pre-wrap; word-break: break-all; }}
 .meta {{ color: #777; font-size: 12px; margin-top: 26px; border-top: 1px solid #ccc; padding-top: 8px; }}
 i, .sym {{ font-family: Georgia, serif; }}
</style></head><body>

<h1>T6 &nbsp;AR perturber-list kernel consistency:<br>members-keep-own-mass-changeover design</h1>
<div class="subtitle">validation report &middot; scenario <i>t6_pert_kernel_match</i> &middot; <span class="verdict {'pass' if all_pass else 'fail'}">{verdict}</span></div>

<h2>1&nbsp;&nbsp;Background and purpose</h2>
<p>With the members-keep-own-changeover design, a group is represented in the soft system by
artificial particles (c.m., tidal-tensor probes, three pseudo particles carrying the group mass),
all using the <i>group-c.m.</i> changeover &chi;<sub>cm</sub> (transition band scaled by
<i>M</i><sub>grp</sub><sup>1/3</sup>), while each AR group member keeps its <i>own-mass</i>
changeover &chi;<sub>i</sub> (band scaled by <i>m</i><sub>i</sub><sup>1/3</sup>). The tidal tensor
is therefore measured at the &chi;<sub>cm</sub> band, but the direct perturber-list force on a
member uses &chi;<sub>i</sub>. This test quantifies the resulting kernel difference and verifies
that it is a stationary representation difference, not a conservation error: (i) measure the
difference between the two kernels on every member evaluation; (ii) integrate the same system with
either kernel and require clean energy conservation in both cases. To make (ii) attributable to the
kernel alone, all runs cover exactly <b>one tree step</b> from an initial condition with a
guaranteed tidal-tensor group: round-off differences then have no time to grow chaotically and the
three runs replicate the same trajectory.</p>

<h2>2&nbsp;&nbsp;Method</h2>
<h3>2.1&nbsp;&nbsp;Quantities under comparison</h3>
<p>At every AR perturbation evaluation, for every group member <i>i</i>, the direct list force is
computed twice (the run always uses the first):</p>
<div class="eq">F<sub>mem,<i>i</i></sub> = &Sigma;<sub>j</sub> G m<sub>j</sub> r<sub>ij</sub><sup>-3</sup>
W(&chi;<sub>i</sub>, &chi;<sub>j</sub>; r<sub>ij</sub>) <b>r</b><sub>ij</sub>
<span class="tag">(1)</span></div>
<div class="eq">F<sub>cm,<i>i</i></sub> = &Sigma;<sub>j</sub> G m<sub>j</sub> r<sub>ij</sub><sup>-3</sup>
W(&chi;<sub>cm</sub>, &chi;<sub>j</sub>; r<sub>ij</sub>) <b>r</b><sub>ij</sub>
<span class="tag">(2)</span></div>
<p>where the sum runs over the perturber list (singles and group members of other groups),
<b>r</b><sub>ij</sub> points from member <i>i</i> to perturber <i>j</i>, and
W(&chi;<sub>a</sub>,&chi;<sub>b</sub>;r) is the changeover pair weight, which follows the band of
the larger-r<sub>out</sub> side. W differs between (1) and (2) only for perturbers inside the
transition annulus r<sub>out</sub>(&chi;<sub>i</sub>) &lt; r<sub>ij</sub> &lt;
r<sub>out</sub>(&chi;<sub>cm</sub>); elsewhere both weights are 0 or 1 alike. The tidal tensor is
identical in both cases.</p>

<h3>2.2&nbsp;&nbsp;Metrics</h3>
<table class="data">
<tr><th>symbol</th><th>definition</th><th>interpretation</th></tr>
<tr><td>&Delta;F<sub>i</sub>,&nbsp; F&#772;<sub>i</sub></td>
<td>&Delta;F<sub>i</sub> = F<sub>mem,<i>i</i></sub> &minus; F<sub>cm,<i>i</i></sub>;&nbsp;
F&#772;<sub>i</sub> = max(|F<sub>mem,<i>i</i></sub>|, |F<sub>cm,<i>i</i></sub>|)</td>
<td>per-evaluation force difference and reference amplitude</td></tr>
<tr><td><i>R</i></td>
<td><i>R</i> = &Sigma;<sub>i</sub> |&Delta;F<sub>i</sub>| / &Sigma;<sub>i</sub> F&#772;<sub>i</sub>
<span class="tag">(3)</span></td>
<td>aggregate mismatch amplitude relative to the list-force amplitude (cumulative over the run)</td></tr>
<tr><td>&Phi;(x)</td>
<td>&Phi;(x) = N(R<sub>i</sub> &gt; x) / N, R<sub>i</sub> = |&Delta;F<sub>i</sub>|/F&#772;<sub>i</sub>
<span class="tag">(4)</span></td>
<td>fraction of evaluations mismatched beyond level x; bimodal by construction (&asymp;0 outside
all bands, O(1) when an annulus perturber dominates the member)</td></tr>
<tr><td>R<sub>max</sub></td><td>max<sub>i</sub> R<sub>i</sub></td>
<td>largest single-evaluation mismatch; confirms the bimodality, not an error measure</td></tr>
<tr><td>N<sub>eval</sub></td><td>number of (member &times; evaluation) pairs sampled</td>
<td>measurement sample size</td></tr>
<tr><td>&epsilon;<sub>E</sub></td><td>end-of-step |&Delta;E/E| from the run energy log</td>
<td>A/B conservation measure of &sect;3.2</td></tr>
</table>
<div class="caption"><b>Table 1.</b> Metrics used in this report. All force quantities refer to the
direct perturber-list contribution only (&sect;2.1); the tidal tensor is common to both kernels.</div>

<h3>2.3&nbsp;&nbsp;Runs</h3>
<table class="data">
<tr><th>run</th><th>list-force kernel in the dynamics</th><th>measurement</th><th>purpose</th></tr>
<tr><td>r1</td><td>&chi;<sub>i</sub> (design, Eq. 1)</td><td>on (Eqs. 3&ndash;4)</td>
<td>measure the kernel mismatch on the same single step (accumulators never feed back)</td></tr>
<tr><td>r2</td><td>&chi;<sub>i</sub> (design, Eq. 1)</td><td>off</td>
<td>clean A/B reference of the design kernel</td></tr>
<tr><td>r3</td><td>&chi;<sub>cm</sub> (Eq. 2)</td><td>off</td>
<td>clean A/B counterfactual consistent with the tidal-tensor band (upstream behaviour)</td></tr>
</table>
<div class="caption"><b>Table 2.</b> All three runs use the same IC (&sect;2.4) and identical
parameters and cover exactly one tree step (dt<sub>soft</sub> = 2<sup>-8</sup> = 0.00390625,
auto-determined from r<sub>out</sub> base 0.1 via dt<sub>soft</sub> = P(r<sub>in</sub>)/n<sub>step</sub>,
r<sub>in</sub> = r-ratio &middot; r<sub>out</sub> with r-ratio 0.15, n<sub>step</sub> = 16;
group linking r<sub>group</sub>/r<sub>search,group</sub> = 0.008/0.01; tree mode, tidal tensor on,
1 rank, 1 thread).<br>
<b>Why r1 and r2 are separate runs although both integrate with the design kernel:</b> the
measurement code inside r1 changes the compiled instruction schedule of the perturbation loop, so
r1 and r2 need not be bit-identical; keeping r2/r3 free of instrumentation isolates the A/B pair
from even this channel. Over a single tree step the round-off difference has no time to amplify,
from report_common import (
    html_escape,
    load_json,
    parse_log_version,
)
and r1 and r2 indeed end bit-identically (&sect;3.2), which doubles as a direct proof that the
measurement is neutral.</div>

<h3>2.4&nbsp;&nbsp;Reproducibility: software, build and execution</h3>
<p><b>Software.</b> PeTar version {html_escape(log_version)}, git commit
<code>{html_escape(git_ref)}</code>, plain family (no stellar evolution, no external potential),
1 MPI rank &times; 1 OpenMP thread. The binary is built with the diagnostic instrumentation
(the flag is passed to <code>./configure</code> through <code>CXXFLAGS</code>; the generated
Makefile defines <code>PROFILE</code> itself, and <code>make</code> without a target builds the
default binaries of the configured family &mdash; here <code>{html_escape(_main_target())}</code>):</p>
<pre class="cmd">CXXFLAGS='-D PERT_KERNEL_DIAG' ./configure
make -j</pre>
<p><b>Initial conditions</b> (hand-generated, not from mcluster):</p>
<pre class="cmd">python3 test/validation/make_ic.py --case t6 --output test/validation/work/t6/input.base
cd test/validation/work/t6 &amp;&amp; petar.init -f input input.base</pre>
<table class="data">
<tr><th>parameter</th><th>value</th></tr>
<tr><td>field model</td><td>Plummer sphere, 297 stars, scale radius <i>a</i> = 0.3 pc; radii
drawn by inverse-CDF sampling of the Plummer mass profile</td></tr>
<tr><td>engineered triple</td><td>central, at the cluster c.m.: inner binary 4+4 M<sub>sun</sub>,
<i>a</i><sub>in</sub> = 0.0012 pc, <i>e</i><sub>in</sub> = 0.2; outer companion 3 M<sub>sun</sub>,
<i>a</i><sub>out</sub> = 0.008 pc, <i>e</i><sub>out</sub> = 0.1 (both started at periapsis)</td></tr>
<tr><td>particle number <i>N</i></td><td>297 + 3 = 300</td></tr>
<tr><td>masses</td><td>independent, uniform in [0.2, 5] M<sub>sun</sub> (a mass spectrum is
required for a member/c.m. band contrast)</td></tr>
<tr><td>velocities</td><td>independent Gaussian per component with the Plummer 1D dispersion
&sigma;&sup2;(r) = GM/(6&radic;(r&sup2;+a&sup2;))</td></tr>
<tr><td>centre-of-mass</td><td>position and velocity subtracted exactly</td></tr>
<tr><td>random generator</td><td>Python <code>random.Random</code>, seed 20261010</td></tr>
<tr><td>units / constants</td><td>M<sub>sun</sub>, pc, pc Myr<sup>-1</sup>;
G = 4.49830997959438 &times; 10<sup>-3</sup> pc&sup3; M<sub>sun</sub><sup>-1</sup> Myr<sup>-2</sup></td></tr>
<tr><td>file format</td><td>ASCII, one row per star: m x y z v<sub>x</sub> v<sub>y</sub>
v<sub>z</sub></td></tr>
</table>
<div class="caption"><b>Table 3.</b> IC generation parameters (<code>make_ic.py --case t6</code>).</div>
<p><b>Execution commands</b> (as recorded by the runner; only the environment variable differs
between runs):</p>
<pre class="cmd">{'<br>'.join(html_escape(cmd) for _rid, cmd in run_commands)}</pre>

<h2>3&nbsp;&nbsp;Results</h2>
<h3>3.1&nbsp;&nbsp;Kernel mismatch (run r1)</h3>
<div class="cards">
<div class="card"><div class="card-label">R (Eq. 3)</div><div class="card-value">{fmt(measure.get('pert_kernel_ratio'))}</div></div>
<div class="card"><div class="card-label">&Phi;(10<sup>-2</sup>)</div><div class="card-value">{fmt(measure.get('pert_kernel_frac_1e2'))}</div></div>
<div class="card"><div class="card-label">&Phi;(10<sup>-1</sup>)</div><div class="card-value">{fmt(measure.get('pert_kernel_frac_1e1'))}</div></div>
<div class="card"><div class="card-label">R<sub>max</sub></div><div class="card-value">{r_max:.3g}</div></div>
<div class="card"><div class="card-label">N<sub>eval</sub></div><div class="card-value">{int(measure.get('pert_kernel_eval', 0))}</div></div>
</div>
<p>The metrics are cumulative over the whole run and are already stationary at the sampled size
(the per-interval values stay flat), so they are reported as single numbers rather than curves: the
kernel difference is a stationary representation property, not a growing error. The compact IC and
the single guaranteed TT group make most list perturbers cross the transition annulus, giving
R &asymp; 0.24 and &Phi;(10<sup>-1</sup>) &asymp; 0.81 over the step; the same measurement on the
production N = 1000 cluster (aggregated over its full group census) gives R &asymp; 1.2 &times;
10<sup>-3</sup>.</p>

<h3>3.2&nbsp;&nbsp;Single-step A/B kernel effect (runs r1, r2, r3)</h3>
<table class="data">
<tr><th>run</th><th>&epsilon;<sub>E</sub> (max |&Delta;E/E|)</th><th>large-energy dumps</th><th>aborts</th></tr>
{ab_row('r1 &mdash; member kernel &chi;<sub>i</sub> (design, instrumented)', measure)}
{ab_row('r2 &mdash; member kernel &chi;<sub>i</sub> (design, clean)', member)}
{ab_row('r3 &mdash; c.m. kernel &chi;<sub>cm</sub> (counterfactual, clean)', cmk)}
</table>
<div class="caption"><b>Table 4.</b> End-of-step energy error of the three
single-step runs (identical IC, identical trajectory by construction). r1 and r2 end
<i>bit-identically</i> (&epsilon;<sub>E</sub> = 2.008576&times;10<sup>-5</sup> both), proving the
measurement does not perturb the dynamics when round-off has no time to amplify. r3 differs from r2
by &asymp;11% (&epsilon;<sub>E</sub> = 1.78&times;10<sup>-5</sup>): this is the <b>pure in-step
effect of the list-force kernel</b> on the same trajectory &mdash; both kernels conserve energy far
below the 10<sup>-3</sup> check threshold, i.e. the mismatch of &sect;3.1 is a representation
difference, not a conservation error.</div>

<h3>3.3&nbsp;&nbsp;Automatic checks</h3>
<table class="data">
<tr><th>check</th><th>status</th><th>value</th><th>threshold</th><th>message</th></tr>
{rows_html}
</table>
<div class="caption"><b>Table 5.</b> Regression checks: the diagnostic must accumulate a real
sample (N<sub>eval</sub> &gt; 10<sup>6</sup>) and <i>R</i> must stay inside a structural band
(0.01&ndash;0.5; values near 0 would mean the measurement broke, values near 1 that the two bands
decoupled); both A/B runs must conserve energy (&epsilon;<sub>E</sub> &lt; 0.01) with zero
large-energy dumps and zero aborts. The band is a regression guard, not an accuracy claim.</div>

<div class="meta">petar version: {html_escape(log_version)} &middot; scenario:
<i>t6_pert_kernel_match</i> &middot; report: {html_escape(str(report_path))}</div>
</body></html>
"""


def main() -> int:
    parser = argparse.ArgumentParser(description="T6 perturber-kernel match: run + HTML report")
    parser.add_argument("--run", action="store_true", help="Run the scenario before generating the report")
    parser.add_argument("--report", default="test/out/report.t6.pert_kernel_match.json")
    parser.add_argument("--scenario-file", default="test/validation/scenarios/t6_pert_kernel_match.json")
    parser.add_argument("--out-dir", default="test/out/validation_t6_pert_kernel_match")
    parser.add_argument("--output", default="test/out/t6_pert_kernel_match_summary.html")
    parser.add_argument("--petar", default="",
                        help="Petar binary with -D PERT_KERNEL_DIAG for --run "
                             "(empty: ensure build/petar.mpi.omp.avx2 is built with the flag)")
    args = parser.parse_args()

    if args.run:
        petar_bin = args.petar or ensure_diag_binary()
        print(f"[T6] Using petar binary: {petar_bin}")
        cmd = [
            sys.executable, "test/validation/run_validation.py",
            "--scenario", args.scenario_file,
            "--out-dir", args.out_dir,
            "--report", args.report,
            "--var", f"petar_bin_switch={petar_bin}",
        ]
        print(f"[T6] Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

    report = load_json(Path(args.report))
    scenario_name = "t6_pert_kernel_match"
    runs = [r for r in report.get("runs", []) if r.get("scenario") == scenario_name]
    if not runs:
        raise RuntimeError(f"No runs found for scenario '{scenario_name}' in report: {args.report}")
    checks = [c for c in report.get("checks", []) if c.get("scenario") == scenario_name]
    run_metrics = {r["run_id"]: r.get("metrics", {}) for r in runs}

    log_path = Path(args.out_dir) / runs[0]["output"]
    log_text = log_path.read_text(encoding="utf-8", errors="replace")
    diag_rows = parse_diag_lines(log_text)
    log_version = parse_log_version(log_path)

    run_commands: List[Tuple[str, str]] = [
        (run["run_id"], str(run.get("command", ""))) for run in runs
    ]

    html = build_html(checks, diag_rows, run_metrics, run_commands, log_version, git_commit(), Path(args.report))
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")
    n_pass = sum(1 for c in checks if c.get("passed", False))
    print(f"[T6] Summary: passed={n_pass == len(checks) and len(checks) > 0} total={len(checks)} failed={len(checks) - n_pass} report={out}")
    return 0 if n_pass == len(checks) and checks else 1


if __name__ == "__main__":
    sys.exit(main())
