#!/usr/bin/env python3
"""T8 LogH vs BTLogH gauge comparison (extracted event classes): HTML report
(paper style, see REPORT_DESIGN.md; structured after the T7 report)."""
import argparse
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

from report_common import (
    PAPER_CSS,
    check_row,
    git_commit,
    html_escape,
    load_json,
    parse_log_version,
    select_or_install_petar,
)

REPO = Path(__file__).resolve().parents[2]

RUN_LABELS = {
    "r1_cascade_logh": ("r1: cascade, LogH gauge", "#d81b60"),
    "r2_cascade_btlogh": ("r2: cascade, BTLogH gauge (--ar-g-func 1)", "#7b1fa2"),
    "r3_chaotic_logh": ("r3: chaotic 6-body, LogH gauge", "#0b74de"),
    "r4_chaotic_btlogh": ("r4: chaotic 6-body, BTLogH gauge", "#2e7d32"),
    "r5_control_logh": ("r5: control (detached pairs), LogH gauge", "#e65100"),
    "r6_control_btlogh": ("r6: control, BTLogH gauge", "#5d4037"),
}


def main() -> int:
    parser = argparse.ArgumentParser(description="T8 gauge comparison: run + HTML report")
    parser.add_argument("--run", action="store_true", help="Run the scenario before generating the report")
    parser.add_argument("--report", default="test/out/report.t8.gauge_logh_btlogh_extracted.json")
    parser.add_argument("--scenario-file", default="test/validation/scenarios/t8_gauge_logh_btlogh_extracted.json")
    parser.add_argument("--out-dir", default="test/out/validation_t8_gauge_logh_btlogh_extracted")
    parser.add_argument("--output", default="test/out/t8_gauge_logh_btlogh_extracted_summary.html")
    parser.add_argument("--petar", default="",
                        help="Petar binary for --run (empty: petar.select --optional mpi,omp,avx2). "
                             "The ICs are plain-format, so a plain-family build is required.")
    args = parser.parse_args()

    if args.run:
        petar_bin = args.petar or select_or_install_petar(optional=["mpi", "omp", "avx2"], repo_root=REPO)
        print(f"[T8] Using petar binary: {petar_bin}")
        cmd = [
            sys.executable, "test/validation/run_validation.py",
            "--scenario", args.scenario_file,
            "--out-dir", args.out_dir,
            "--report", args.report,
            "--var", f"petar_bin_switch={petar_bin}",
        ]
        print(f"[T8] Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

    report = load_json(Path(args.report))
    scenario_name = "t8_gauge_logh_btlogh_extracted"
    runs = [r for r in report.get("runs", []) if r.get("scenario") == scenario_name]
    if not runs:
        raise RuntimeError(f"No runs found for scenario '{scenario_name}' in report: {args.report}")
    checks = [c for c in report.get("checks", []) if c.get("scenario") == scenario_name]

    table_rows: List[str] = []
    run_records: List[Dict[str, Any]] = []
    version = "unknown"
    metrics_by_run: Dict[str, Dict[str, Any]] = {}
    for run in runs:
        run_id = run["run_id"]
        log_path = Path(args.out_dir) / str(run.get("output", f"{run_id}.log"))
        if log_path.exists():
            version = parse_log_version(log_path)
        metrics = run.get("metrics", {})
        metrics_by_run[run_id] = metrics
        counts = metrics.get("regex_counts", {})
        floats = metrics.get("regex_floats", {})
        label = RUN_LABELS.get(run_id, (run_id, "#555"))[0]
        ar_steps = floats.get("ar_step_sum", float("nan"))
        ar_cell = f"{ar_steps:.3e}" if ar_steps == ar_steps else "&mdash;"
        table_rows.append(
            f"<tr><td>{html_escape(run_id)}</td><td>{html_escape(label)}</td>"
            f"<td>{ar_cell}</td>"
            f"<td>{counts.get('hard_energy_sig', -1)}</td>"
            f"<td>{counts.get('step_warning', -1)}</td><td>{counts.get('ds_floor_abort', -1)}</td>"
            f"<td>{counts.get('runtime_abort', -1)}</td></tr>"
        )
        run_records.append({"run_id": run_id, "command": run["command"]})

    def _steps(run_id: str) -> float:
        return metrics_by_run.get(run_id, {}).get("regex_floats", {}).get("ar_step_sum", float("nan"))

    def _fmt(x: float) -> str:
        return f"{x:.3e}" if x == x else "&mdash;"

    ratio_rows = []
    for name, a, b in [
        ("cascade (r2/r1)", "r2_cascade_btlogh", "r1_cascade_logh"),
        ("chaotic (r4/r3)", "r4_chaotic_btlogh", "r3_chaotic_logh"),
        ("control (r6/r5)", "r6_control_btlogh", "r5_control_logh"),
    ]:
        sa, sb = _steps(a), _steps(b)
        ratio = sa / sb if (sa == sa and sb == sb and sb > 0) else float("nan")
        ratio_cell = f"{ratio:.3f}" if ratio == ratio else "&mdash;"
        ratio_rows.append(f"<tr><td>{html_escape(name)}</td><td>{_fmt(sb)}</td><td>{_fmt(sa)}</td><td>{ratio_cell}</td></tr>")
    ratio_html = "".join(ratio_rows)

    all_pass = all(c.get("passed", False) for c in checks)
    verdict = "ALL CHECKS PASSED" if all_pass else "CHECKS FAILED"

    command_rows_html = "".join(
        f"<tr><td>{html_escape(r['run_id'])}</td><td><code>{html_escape(r['command'])}</code></td></tr>"
        for r in run_records
    )
    checks_html = "".join(check_row(c) for c in checks)

    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>T8 LogH vs BTLogH gauge comparison (extracted event classes)</title>
<style>{PAPER_CSS}</style></head><body>

<h1>T8 &nbsp;LogH vs BTLogH gauge comparison<br>on extracted multiple-system event classes</h1>
<div class="subtitle">validation report &middot; scenario <i>t8_gauge_logh_btlogh_extracted</i> &middot; <span class="verdict {'pass' if all_pass else 'fail'}">{verdict}</span></div>

<h2>1&nbsp;&nbsp;Background</h2>
<p>PeTar's AR (algorithmic-regularization) integrator supports two
time-transform gauges: the default <b>LogH</b> sum gauge,
&Gamma;<sub>LogH</sub> = ln[(E<sub>kin</sub> &minus; E<sub>0</sub>)/(−U)],
and the <b>BTLogH</b> product gauge (<code>--ar-g-func 1</code>), in which
each binary-tree level contributes a factor U<sub>node</sub>/&kappa; to a
product &Pi;g (see <code>SDAR/docs/hierarchical_blogh_impl_notes.md</code>).
The two gauges integrate the same physical trajectories but price resolution
differently: for <i>hierarchical</i> multiples the product form decouples the
levels, while for non-hierarchical chaotic complexes it can over-resolve. A
2026-10-22 dump-replay sweep on the N=1000 100%-binary benchmark
(<code>~/sdar_run/bench-n1k-arfix/</code>) measured per-event step ratios of
38&ndash;934&times; (BTLogH cheaper) on hierarchical binary&ndash;binary
cascades, parity on ordinary blocks, and a 4&times; step <i>cost</i> with a
cleaner energy monitor (0 vs 16 hard-energy events) on one hard chaotic
cluster. This suite pins those three signatures on self-contained
deterministic ICs so the relative gauge behavior is guarded across future
code changes. No stellar evolution and no external potential are involved.</p>

<h2>2&nbsp;&nbsp;Purpose</h2>
<ol>
<li><b>Hierarchical-advantage guard (r1, r2):</b> on the cascade IC the LogH
run must pay a real event (N<sub>AR</sub> &ge; 3&times;10<sup>3</sup> floor)
inside its ceiling (10<sup>5</sup>), and the BTLogH run must stay
&le; 6&times;10<sup>3</sup> steps with the gauge ratio
N<sub>AR,BTLogH</sub>/N<sub>AR,LogH</sub> &isin; [0.02, 0.5] &mdash; i.e. the
product gauge retains its multi-level advantage (dump-replay per-block
values 1/934&ndash;1/38; standalone-window value ~0.11).</li>
<li><b>Chaotic-accuracy guard (r3, r4):</b> on the 6-body exchange complex
both gauges complete within 10<sup>5</sup> steps, the LogH event floor
(3&times;10<sup>3</sup>) confirms the resonance triggers, and the BTLogH arm
must show <b>zero</b> hard-energy-significant events (measured: LogH 2,
BTLogH 0; dump-replay parent cluster: 16 vs 0) &mdash; the documented
accuracy-not-waste signature of the extra resolution the product gauge spends
on non-hierarchical complexes.</li>
<li><b>Parity / false-positive guard (r5, r6):</b> on two detached binaries
(no encounter within the window) both gauges stay &le; 5&times;10<sup>3</sup>
steps with zero hard-energy events &mdash; the cascade/chaotic assertions do
not fire on a quiescent configuration, and neither gauge pays a penalty
there.</li>
<li><b>Stability guards (all runs):</b> zero runtime aborts (assert/segfault)
and zero adaptive-ds time-resolution aborts in every run &mdash; both gauges
must remain inside the post-2026-10-21 fixed controller regime.</li>
<li><b>Scope note:</b> this suite guards <i>per-event</i> gauge signatures.
The full-cluster N=1000 BTLogH run aborts (ds denormal collapse at
t&asymp;18.5 Myr, estimator-channel suspicion) &mdash; tracked as SDAR plan
open item 5 (<code>etot-sd-ref-tracking-corruption-plan.md</code>), out of
T8's scope.</li>
</ol>

<h2>3&nbsp;&nbsp;Method</h2>
<h3>3.1&nbsp;&nbsp;Initial conditions</h3>
<p>All three ICs are deterministic full-double-precision particle tables
(hard-coded in <code>make_ic.py</code>, no RNG), extracted from
hard-integrator dump states of the 2026-10-22 benchmark runs, shifted to the
center-of-mass frame (position and velocity), and reduced to the physically
important members:</p>
<table class="data">
<tr><th>case</th><th>members</th><th>key orbital parameters (at extraction)</th></tr>
<tr><td>cascade (r1/r2)<br>11 bodies, 133.3 M<sub>sun</sub></td>
<td>the 4-body AR group of the n30 cluster core at mutual pericenter
&mdash; a 105 M<sub>sun</sub> tight binary + a 0.3 M<sub>sun</sub> binary
&mdash; plus the 7 bystanders carrying 95% of the cluster's tidal field on
the group (minimal set by G&middot;m/r<sup>3</sup>)</td>
<td>massive binary: m = 21.63+83.74, a = 6.48&times;10<sup>-7</sup> pc,
e = 0.016, P = 4.76&times;10<sup>-9</sup> Myr; small binary: m = 0.20+0.10,
a = 3.24&times;10<sup>-5</sup> pc, e = 0.65; root (CM&ndash;CM):
q &asymp; 6.3&times;10<sup>-3</sup> pc, e &asymp; 0.94 &mdash; the dump
state is at the mutual pericenter passage. Without the 7 bystanders the
external field is too weak to disengage the slowdown and both gauges
integrate the 4-body extract in ~300 steps &mdash; the disengagement, and
hence the cascade cost, is externally driven.</td></tr>
<tr><td>chaotic (r3/r4)<br>6 bodies, 107.3 M<sub>sun</sub></td>
<td>the step-dominant 6-body AR group of a dense 31-star core (the same
massive binary in a 4+2 low-mass exchange complex)</td>
<td>massive binary: same 105 M<sub>sun</sub> pair; four low-mass members
(0.11&ndash;1.33 M<sub>sun</sub>) on crossing orbits within
&sim;10<sup>-2</sup> pc of the binary</td></tr>
<tr><td>control (r5/r6)<br>4 bodies, 105.6 M<sub>sun</sub></td>
<td>two detached binaries (CM separation 0.093 pc, no encounter within the
window)</td>
<td>massive binary + one low-mass binary at large separation; both gauges
resolve the two independent binaries with slowdown engaged &mdash; the
quiescent reference point</td></tr>
</table>
<div class="caption"><b>Table 1.</b> Initial conditions. Source dumps:
<code>current.dump_large_step_h4_n30_g10_t5.929199</code> (cascade),
<code>baseline.dump_large_step_h4_n31_g10_t0.359375</code> (chaotic),
<code>baseline.dump_large_step_h4_n32_g13_t0.336914</code> (control), all in
<code>~/sdar_run/bench-n1k-arfix/</code>.</div>

<h3>3.2&nbsp;&nbsp;Metrics</h3>
<table class="data">
<tr><th>symbol</th><th>source</th><th>role in the checks</th></tr>
<tr><td>N<sub>AR</sub></td><td>cumulative AR substeps over the window
(<code>AR_step_sum</code>, final statistics block of the run log)</td>
<td>quantitative budget per run (floors guard that the event actually
triggers; ceilings guard against pathological cost), and the cross-gauge
ratio for the cascade case</td></tr>
<tr><td>N<sub>HES</sub></td><td>count of &ldquo;Hard energy significant&rdquo;
monitor events (|dE<sub>SD</sub>/E<sub>SD</sub>| above the hard-monitor
threshold)</td><td>must be 0 for BTLogH on the chaotic IC and for both gauges
on the control IC; informational (allowed &gt; 0) for the LogH chaotic and
cascade arms &mdash; transient reductions during violent encounters are
legitimate</td></tr>
<tr><td>N<sub>warn</sub></td><td>count of AR step-cap warnings</td>
<td>informational; must stay 0 (steps are far below the cap)</td></tr>
<tr><td>N<sub>floor</sub></td><td>count of adaptive-ds time-resolution
aborts</td><td>must be 0 everywhere</td></tr>
<tr><td>N<sub>abort</sub></td><td>runtime aborts (segfault/assertion)</td>
<td>must be 0 everywhere</td></tr>
</table>
<div class="caption"><b>Table 2.</b> Metrics and their roles.</div>

<h3>3.3&nbsp;&nbsp;Reproducibility</h3>
<p><b>Software.</b> PeTar version {html_escape(version)}, git commit
<code>{html_escape(git_commit(REPO))}</code>, 1 rank &times; 1 thread,
plain-family binary (no stellar-evolution or external-potential columns),
selected via <code>petar.select --optional mpi,omp,avx2</code>. Solver radii
are scaled to the IC (<code>--r 0.05 --r-ratio 0.15 --r-search-group
0.005</code>): the hard manager asserts
r<sub>search,group</sub> &le; r<sub>in</sub>, and the compact extracted
systems require compact radii. Units: G = 4.4983&times;10<sup>-3</sup>
pc<sup>3</sup> M<sub>sun</sub><sup>-1</sup> Myr<sup>-2</sup> (<code>-u 1</code>);
soft step s = 9.77&times;10<sup>-4</sup> Myr; windows: cascade/chaotic
t = 0.01 Myr, control t = 0.05 Myr.</p>
<p><b>Initial conditions</b> (deterministic, hardcoded numbers, full double
precision, no RNG):</p>
<pre class="cmd">python3 test/validation/make_ic.py --case t8_cascade --output .../input.cascade
python3 test/validation/make_ic.py --case t8_chaotic --output .../input.chaotic
python3 test/validation/make_ic.py --case t8_control --output .../input.control
cd test/validation/work/t8 &amp;&amp; petar.init -f input_cascade input.cascade \
  &amp;&amp; petar.init -f input_chaotic input.chaotic &amp;&amp; petar.init -f input_control input.control</pre>
<p><b>Execution commands</b> (as recorded by the runner):</p>
<table class="data">
<tr><th>run_id</th><th>command</th></tr>
{command_rows_html}
</table>

<h2>4&nbsp;&nbsp;Results</h2>
<h3>4.1&nbsp;&nbsp;Run summary</h3>
<table class="data">
<tr><th>run</th><th>configuration</th><th>N<sub>AR</sub></th><th>N<sub>HES</sub></th><th>N<sub>warn</sub></th><th>N<sub>floor</sub></th><th>N<sub>abort</sub></th></tr>
{''.join(table_rows)}
</table>
<div class="caption"><b>Table 3.</b> Per-run quantities. r1/r2 are the
hierarchical cascade pair (primary gauge-ratio guard); r3/r4 the chaotic
exchange complex (accuracy signature); r5/r6 the quiescent control. The
dump-replay reference values on the source events (full cluster environment)
are 2.2&times;10<sup>6</sup>&ndash;1.1&times;10<sup>7</sup> (LogH) vs
1.2&times;10<sup>4</sup>&ndash;3.5&times;10<sup>4</sup> (BTLogH) per hard
block for the cascade, and 6.6&times;10<sup>4</sup> vs 2.6&times;10<sup>5</sup>
steps (16 vs 0 HES events) for the chaotic parent cluster &mdash; documented
for context; the suite runs exclusively against the standalone ICs.</div>

<h3>4.2&nbsp;&nbsp;Gauge step ratios</h3>
<table class="data">
<tr><th>case</th><th>N<sub>AR</sub> (LogH)</th><th>N<sub>AR</sub> (BTLogH)</th><th>ratio BTLogH/LogH</th></tr>
{ratio_html}
</table>
<div class="caption"><b>Table 4.</b> Cross-gauge step ratios. Cascade:
BTLogH multi-level advantage (guarded band [0.02, 0.5]). Chaotic: near-parity
to mild advantage &mdash; the signature here is the monitor cleanliness, not
the step count. Control: parity within run-to-run variation.</div>

<h3>4.3&nbsp;&nbsp;Automatic checks</h3>
<table class="data">
<tr><th>check</th><th>status</th><th>value</th><th>threshold</th><th>message</th></tr>
{checks_html}
</table>
<div class="caption"><b>Table 5.</b> Automatic checks, calibrated against the
current binary (all values measured on these ICs with ~3&times; margin on
floors/ceilings; the cascade ratio band is widened to [0.02, 0.5] to absorb
gauge-cost evolution while still excluding parity).</div>

<div class="footer">version {html_escape(version)} &middot; commit {git_commit(REPO)} &middot; generated by t8_gauge_logh_btlogh_extracted_report.py</div>
</body></html>
"""
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")
    n_pass = sum(1 for c in checks if c.get("passed", False))
    print(f"[T8] Summary: passed={n_pass == len(checks) and len(checks) > 0} total={len(checks)} failed={len(checks) - n_pass} report={out}")
    return 0 if n_pass == len(checks) and checks else 1


if __name__ == "__main__":
    raise SystemExit(main())
