#!/usr/bin/env python3
"""T7 changeover rate-type ds-controller HTML report (paper style, see REPORT_DESIGN.md)."""
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
    "r1_pert_band_mid": ("r1: r_p = 0.15 (intruder band, mid)", "#d81b60"),
    "r2_pert_inside_rin": ("r2: r_p = 0.05 (inside r_in, k = 1)", "#0b74de"),
    "r3_pert_band_far": ("r3: r_p = 0.30 (intruder band, far side)", "#2e7d32"),
    "r4_pert_band_mid_btlogh": ("r4: r_p = 0.15, BTLogH gauge (--ar-g-func 1)", "#7b1fa2"),
    "r5_pert_outer": ("r5: r_p = 0.60 (outside r_out)", "#e65100"),
    "r6_pert_abort_band_near": ("r6: r_p = 0.35 (former abort band, near)", "#00838f"),
    "r7_pert_abort_band_far": ("r7: r_p = 0.50 (former abort band, far)", "#5d4037"),
}


def main() -> int:
    parser = argparse.ArgumentParser(description="T7 changeover rate-type drift: run + HTML report")
    parser.add_argument("--run", action="store_true", help="Run the scenario before generating the report")
    parser.add_argument("--report", default="test/out/report.t7.triple_perturber_ar_ds_adjustment.json")
    parser.add_argument("--scenario-file", default="test/validation/scenarios/t7_triple_perturber_ar_ds_adjustment.json")
    parser.add_argument("--out-dir", default="test/out/validation_t7_triple_perturber_ar_ds_adjustment")
    parser.add_argument("--output", default="test/out/t7_triple_perturber_ar_ds_adjustment_summary.html")
    parser.add_argument("--petar", default="",
                        help="Petar binary for --run (empty: petar.select --optional mpi,omp,avx2). "
                             "The IC is plain-format, so a plain-family build is required.")
    args = parser.parse_args()

    if args.run:
        petar_bin = args.petar or select_or_install_petar(optional=["mpi", "omp", "avx2"], repo_root=REPO)
        print(f"[T7] Using petar binary: {petar_bin}")
        cmd = [
            sys.executable, "test/validation/run_validation.py",
            "--scenario", args.scenario_file,
            "--out-dir", args.out_dir,
            "--report", args.report,
            "--var", f"petar_bin_switch={petar_bin}",
        ]
        print(f"[T7] Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

    report = load_json(Path(args.report))
    scenario_name = "t7_triple_perturber_ar_ds_adjustment"
    runs = [r for r in report.get("runs", []) if r.get("scenario") == scenario_name]
    if not runs:
        raise RuntimeError(f"No runs found for scenario '{scenario_name}' in report: {args.report}")
    checks = [c for c in report.get("checks", []) if c.get("scenario") == scenario_name]

    table_rows: List[str] = []
    run_records: List[Dict[str, Any]] = []
    version = "unknown"
    for run in runs:
        run_id = run["run_id"]
        log_path = Path(args.out_dir) / str(run.get("output", f"{run_id}.log"))
        if log_path.exists():
            version = parse_log_version(log_path)
        metrics = run.get("metrics", {})
        counts = metrics.get("regex_counts", {})
        floats = metrics.get("regex_floats", {})
        label = RUN_LABELS.get(run_id, (run_id, "#555"))[0]
        ar_steps = floats.get("ar_step_sum", float("nan"))
        ar_cell = f"{ar_steps:.3e}" if ar_steps == ar_steps else "&mdash;"
        energy = float(metrics.get("max_abs_error_over_total", float("nan")))
        energy_cell = f"{energy:.2e}" if energy == energy else "&mdash;"
        table_rows.append(
            f"<tr><td>{html_escape(run_id)}</td><td>{html_escape(label)}</td>"
            f"<td>{ar_cell}</td><td>{energy_cell}</td>"
            f"<td>{counts.get('step_warning', -1)}</td><td>{counts.get('ds_floor_abort', -1)}</td>"
            f"<td>{counts.get('runtime_abort', -1)}</td><td>{counts.get('large_energy', -1)}</td></tr>"
        )
        run_records.append({"run_id": run_id, "command": run["command"]})

    all_pass = all(c.get("passed", False) for c in checks)
    verdict = "ALL CHECKS PASSED" if all_pass else "CHECKS FAILED"

    command_rows_html = "".join(
        f"<tr><td>{html_escape(r['run_id'])}</td><td><code>{html_escape(r['command'])}</code></td></tr>"
        for r in run_records
    )
    checks_html = "".join(check_row(c) for c in checks)

    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>T7 AR adaptive-step adjustment: perturbed triple</title>
<style>{PAPER_CSS}</style></head><body>

<h1>T7 &nbsp;AR adaptive-step (ds) adjustment for a<br>hierarchical triple with an external perturber</h1>
<div class="subtitle">validation report &middot; scenario <i>t7_triple_perturber_ar_ds_adjustment</i> &middot; <span class="verdict {'pass' if all_pass else 'fail'}">{verdict}</span></div>

<h2>1&nbsp;&nbsp;Background</h2>
<p>The AR integrator's adaptive step size <code>ds</code> is controlled by a
per-step error audit of the transformed Hamiltonian
H = ln[(E<sub>kin,sd</sub> &minus; E<sub>ref,sd</sub>) / (&minus;E<sub>pot,sd</sub>)],
built only from the group's <i>internal</i> (slowdown-weighted) energies. When
an AR group is perturbed from outside, energy physically flows through the group
boundary, and the per-step H change acquires a component that scales with the
time step (a <i>rate</i>), not with ds<sup>order</sup> (truncation). If the
external perturber force is additionally modulated by the members' changeover
weight k(r) &mdash; which varies continuously as the group moves &mdash; this
rate component becomes continuous and large: with a tight hierarchical trio and
a perturber inside the changeover transition band it reaches ~3.6&times;10<sup>-5</sup>
per step at the natural ds. The root cause (SDAR 2026-10-21, commit
<code>6a8a357</code>) is a <i>bookkeeping</i> inconsistency, not a controller
misclassification: the SD reference E<sub>ref,sd</sub> absorbed the perturbation
work at <i>full</i> weight while the slowed system takes the &kappa;-weighted
part, so the residual R = S &minus; E<sub>ref,sd</sub> drifted at
&minus;(1&minus;1/&kappa;)&middot;W per unit work &mdash; and the controller
collapsed ds by ~5 orders of magnitude (10<sup>6</sup>-substep caps; the effect
was originally found as a pathological slowdown of hard integration in cluster
runs) in response to a genuinely drifting audit. The fix credits the work at
the slowdown system's own velocity-gradient (v<sub>sd</sub>) rate, removes the
n&ge;3 <code>fix_step_option=later</code> auto-grant that froze ds through
resolution swings, and adds physically-gated sd-gauge repairs; an interim
rate-drift EMA compensation (2026-10-19) was removed once the root cause
landed. This suite pins the trigger in a minimal pure-Newtonian 4-body system
&mdash; no stellar evolution, no external potential, standard
plain-format input.</p>

<h2>2&nbsp;&nbsp;Purpose</h2>
<ol>
<li><b>Regression guard (r1, r3):</b> with a perturber placed inside the
intruder's changeover transition band, the cumulative AR substep count stays
within budget (N<sub>AR</sub> &le; 10<sup>6</sup>; current binary ~1.3&times;10<sup>4</sup>)
and no step-cap warnings, ds-floor aborts or runtime aborts occur. A
reintroduced full-weight work credit trips these checks long before the
warning cap.</li>
<li><b>False-positive guard (r2):</b> with the perturber deep inside r<sub>in</sub>
(k = 1, no modulation) the same assertions hold &mdash; the checks do not fire on
an untriggered configuration.</li>
<li><b>Energy guard (r1):</b> the physical energy error stays within a gross
bound, guarding against the step-size budget being met by sacrificing accuracy.</li>
<li><b>Documentation of the zone structure (r1, r2, r3, r5, r6, r7):</b> the
six placements sample the changeover band (k &isin; (0,1)), its inside,
beyond r<sub>out</sub>, and the former ds-floor-abort band. r5 (beyond
r<sub>out</sub>) was formerly excluded because it aborted from
<code>fix_step_option=later</code> freezing ds &mdash; a gap removed by the
same 2026-10-21 root fix (n&ge;3 auto-grant removed). r6/r7 (0.35/0.50)
pin the oscillatory ds-floor-abort zone (pre-fix: bit-identical NaN
aborts old vs interim fix, from the frozen-ds one-step residual)
likewise fixed by the same commit; all three now carry abort-to-complete
guards (no step-cap warning, no ds-floor abort). The T7 trigger runs
with <code>Fix: none</code>, i.e. the adaptive path under guard is
active (verified from the pre-fix warning context).</li>
<li><b>Gauge comparison (r4):</b> the r1 configuration re-run with the BTLogH
g-function (<code>--ar-g-func 1</code>) exercises the same &kappa;-weighted
work-credit path under the alternative time-transform gauge; it must equally
stay within budget (~1.1&times;10<sup>5</sup> steps measured &mdash; consistently
~9&times; the LogH count on this IC, an expected gauge-cost difference, not a
failure mode) with zero warnings/aborts.</li>
</ol>

<h2>3&nbsp;&nbsp;Method</h2>
<h3>3.1&nbsp;&nbsp;System and zone placement</h3>
<table class="data">
<tr><th>component</th><th>parameters</th></tr>
<tr><td>inner binary</td><td>m = 23.66 + 19.04 M<sub>sun</sub>;
a<sub>in</sub> = 8.34&times;10<sup>-7</sup> pc, e<sub>in</sub> &asymp; 5&times;10<sup>-4</sup>
(near-circular), P<sub>in</sub> = 1.09&times;10<sup>-8</sup> Myr</td></tr>
<tr><td>intruder (trio member)</td><td>m = 0.084 M<sub>sun</sub>; barely unbound
orbit around the binary (e &asymp; 1.001), lingering at r &asymp; 0.036 pc
(&asymp; 4&times;10<sup>4</sup> a<sub>in</sub>) with v &asymp; v<sub>esc</sub>
&mdash; the slowly evolving third body whose orbital motion continuously
modulates the group state</td></tr>
<tr><td>perturber (external)</td><td>m = 1.0 M<sub>sun</sub>, static, on +x at
r<sub>p</sub> &isin; {{0.05, 0.15, 0.30, 0.35, 0.50, 0.60}} pc; its plain-Newtonian force on the
group members is weighted by the changeover factor k(r)</td></tr>
<tr><td>changeover radii (from <code>--r 0.6 --r-ratio 0.15</code>, mass-scaled)</td>
<td>intruder band [r<sub>in</sub>, r<sub>out</sub>] = [0.082, 0.545] pc;
binary-member band &asymp; [0.5, 3.6] pc &rarr; r1/r3 (0.15/0.30) sit inside the
intruder band (k &isin; (0,1), modulated); r2 (0.05) sits inside r<sub>in</sub>
(k = 1, unmodulated); r5 (0.60) sits beyond r<sub>out</sub> (k = 0,
unperturbed AR bookkeeping)</td></tr>
<tr><td>group</td><td>the trio forms one AR group (r-search-group = 0.05 pc);
the perturber stays external and enters the group's perturber list</td></tr>
<tr><td>units / window</td><td>G = 4.4983&times;10<sup>-3</sup> pc&sup3;
M<sub>sun</sub><sup>-1</sup> Myr<sup>-2</sup>; tree step 0.0156 Myr, integration
window t = 0.0078 Myr (the trigger develops within the first tree step)</td></tr>
</table>
<div class="caption"><b>Table 1.</b> Initial conditions. The trio is the exact
c.m.-subtracted state of a hierarchical triple extracted from a cluster
simulation snapshot (full double precision &mdash; the knife-edge intruder orbit
is precision-sensitive); the perturber is generated analytically.</div>

<h3>3.2&nbsp;&nbsp;Metrics</h3>
<table class="data">
<tr><th>symbol</th><th>source</th><th>role in the checks</th></tr>
<tr><td>N<sub>AR</sub></td><td>cumulative AR substeps over the window (final statistics block of the run log)</td>
<td>quantitative budget, r1: &le; 10<sup>6</sup> (current ~1.25&times;10<sup>4</sup> LogH, ~1.15&times;10<sup>5</sup> BTLogH r4; pre-fix behaviour on this IC: 2.27&times;10<sup>7</sup> with a step-cap warning &mdash; the threshold sits 80&times; (r1) and 8.7&times; (r4) above the current values and 23&times; below the pre-fix one, without needing a pre-fix binary)</td></tr>
<tr><td>|dE/E|<sub>max</sub></td><td>max over outputs of the physical energy error (Error/Total) from the run log</td>
<td>accuracy guard, r1: &le; 0.5 (measured ~0.09 on this violent-encounter IC; the bound is a gross-breakage guard, not a precision claim)</td></tr>
<tr><td>N<sub>warn</sub></td><td>count of 10<sup>6</sup>-substep step-cap warnings</td>
<td>must be 0 everywhere</td></tr>
<tr><td>N<sub>floor</sub></td><td>count of adaptive-ds time-resolution aborts</td>
<td>must be 0 everywhere</td></tr>
<tr><td>N<sub>abort</sub></td><td>runtime aborts (segfault/assertion)</td>
<td>must be 0 everywhere</td></tr>
<tr><td>N<sub>LEE</sub></td><td>Large_energy_error events (error-reduction actions)</td>
<td>informational only; genuine transient spikes may legitimately occur during the close encounter</td></tr>
</table>
<div class="caption"><b>Table 2.</b> Metrics and their roles.</div>

<h3>3.3&nbsp;&nbsp;Reproducibility</h3>
<p><b>Software.</b> PeTar version {html_escape(version)}, git commit
<code>{html_escape(git_commit(REPO))}</code>, 1 rank &times; 1 thread,
plain-family binary (no stellar-evolution or external-potential columns),
selected via <code>petar.select --optional mpi,omp,avx2</code>.</p>
<p><b>Initial conditions</b> (deterministic, hardcoded numbers, full double
precision, no RNG):</p>
<pre class="cmd">python3 test/validation/make_ic.py --case t7         --output .../input.band    # r_p = 0.15
python3 test/validation/make_ic.py --case t7_band2   --output .../input.band2   # r_p = 0.30
python3 test/validation/make_ic.py --case t7_control_in --output .../input.inside # r_p = 0.05
python3 test/validation/make_ic.py --case t7_outer      --output .../input.outer   # r_p = 0.60
python3 test/validation/make_ic.py --case t7_abortzone  --output .../input.abortzone  # r_p = 0.35
python3 test/validation/make_ic.py --case t7_abortzone2 --output .../input.abortzone2 # r_p = 0.50
cd test/validation/work/t7 &amp;&amp; petar.init -f input_band input.band   # (+ the two others)</pre>
<p><b>Execution commands</b> (as recorded by the runner):</p>
<table class="data">
<tr><th>run_id</th><th>command</th></tr>
{command_rows_html}
</table>

<h2>4&nbsp;&nbsp;Results</h2>
<h3>4.1&nbsp;&nbsp;Run summary</h3>
<table class="data">
<tr><th>run</th><th>placement</th><th>N<sub>AR</sub></th><th>|dE/E|<sub>max</sub></th><th>N<sub>warn</sub></th><th>N<sub>floor</sub></th><th>N<sub>abort</sub></th><th>N<sub>LEE</sub></th></tr>
{''.join(table_rows)}
</table>
<div class="caption"><b>Table 3.</b> Per-run quantities. The r1 row is the
primary regression guard (budget + zero warnings/aborts + energy bound); r3
repeats it at a second in-band point; r2 demonstrates the checks stay silent
without the changeover modulation; r4 repeats r1 under the BTLogH gauge; r5
extends coverage beyond r<sub>out</sub>; r6/r7 pin the former ds-floor-abort
band (pre-fix NaN-abort zone) as abort-to-complete guards. N<sub>LEE</sub> &gt; 0 in the in-band runs
reflects genuine transient error-reduction actions during the close encounter
&mdash; allowed by design.</div>

<h3>4.2&nbsp;&nbsp;Automatic checks</h3>
<table class="data">
<tr><th>check</th><th>status</th><th>value</th><th>threshold</th><th>message</th></tr>
{checks_html}
</table>
<div class="caption"><b>Table 4.</b> Automatic checks. Calibrated against the
current binary (all values measured on this IC); the pre-fix reference values
(N<sub>AR</sub> &asymp; 2.27&times;10<sup>7</sup>, 1 warning) are documented for
context only &mdash; the suite runs exclusively against the current build.</div>

<div class="meta">petar version: {html_escape(version)} &middot; scenario:
<i>t7_triple_perturber_ar_ds_adjustment</i> &middot; report: {html_escape(args.report)}</div>
</body></html>
"""

    Path(args.output).write_text(html, encoding="utf-8")
    n_pass = sum(1 for c in checks if c.get("passed", False))
    print(f"[T7] Summary: passed={n_pass == len(checks) and len(checks) > 0} total={len(checks)} failed={len(checks) - n_pass} report={args.output}")
    return 0 if n_pass == len(checks) and checks else 1


if __name__ == "__main__":
    raise SystemExit(main())
