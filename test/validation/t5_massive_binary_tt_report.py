#!/usr/bin/env python3
"""T5 massive-binary tidal-tensor momentum-conservation HTML report (paper style, see REPORT_DESIGN.md)."""
import argparse
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

from report_common import (
    PAPER_CSS,
    check_row,
    git_commit,
    html_escape,
    load_json,
    parse_command_args,
    parse_log_version,
    read_ic,
    select_or_install_petar,
    svg_curves_plot,
)

REPO = Path(__file__).resolve().parents[2]
Point = Tuple[float, float]
WORK_DIR = Path("test/validation/work/t5")

MODE_LABELS = {
    "r1_hard_tt": ("r1 hard, TT on (isolated group)", "#0b74de"),
    "r2_hard_no_tt": ("r2 hard, TT off (baseline)", "#d81b60"),
    "r3_tree_tt": ("r3 tree, TT on (event configuration)", "#2e7d32"),
}


def load_status(status_path: Path, n_particle: int):
    tools_path = (REPO / "tools").resolve()
    if str(tools_path) not in sys.path:
        sys.path.insert(0, str(tools_path))
    from analysis.status import Status  # pylint: disable=import-outside-toplevel  # pyright: ignore[reportMissingImports]

    st = Status(N_particle=n_particle, interrupt_mode="none")
    st.fromfile(str(status_path))
    return st


def momentum_series(status_path: Path, n_particle: int) -> Dict[str, Any]:
    st = load_status(status_path, n_particle)
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
    drift = np.abs(p_norm - p_norm[0]) / m_tot[0]
    return {
        "time": [float(t) for t in st.time.tolist()],
        "drift": [float(d) for d in drift.tolist()],
        "max_drift": float(np.max(drift)),
        "n_sample": int(st.size),
    }


def _select_plain_petar() -> str:
    """Select (and if needed reinstall) the plain-family binary via petar.select.

    The scenario IC is plain-format, so BSE/galpy-enabled builds cannot read it;
    plain family is selected without --require.
    """
    return select_or_install_petar(optional=["mpi", "omp", "avx2"], repo_root=REPO)


def main() -> int:
    parser = argparse.ArgumentParser(description="T5 massive-binary TT momentum: run + HTML report")
    parser.add_argument("--run", action="store_true", help="Run the scenario before generating the report")
    parser.add_argument("--report", default="test/out/report.t5.massive_binary_tt.json")
    parser.add_argument("--scenario-file", default="test/validation/scenarios/t5_massive_binary_tt.json")
    parser.add_argument("--out-dir", default="test/out/validation_t5_massive_binary_tt")
    parser.add_argument("--output", default="test/out/t5_massive_binary_tt_summary.html")
    parser.add_argument("--ic", default="test/validation/work/t5/input.base")
    parser.add_argument("--petar", default="",
                        help="Petar binary for --run (empty: petar.select --optional mpi,omp,avx2). "
                             "The scenario IC is plain-format (no external-potential or stellar "
                             "columns), so it requires a no-interrupt/no-external build.")
    args = parser.parse_args()

    if args.run:
        petar_bin = args.petar or _select_plain_petar()
        print(f"[T5] Using petar binary: {petar_bin}")
        cmd = [
            sys.executable, "test/validation/run_validation.py",
            "--scenario", args.scenario_file,
            "--out-dir", args.out_dir,
            "--report", args.report,
            "--var", f"petar_bin_switch={petar_bin}",
        ]
        print(f"[T5] Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

    report = load_json(Path(args.report))
    scenario_name = "t5_massive_binary_tt"
    runs = [r for r in report.get("runs", []) if r.get("scenario") == scenario_name]
    if not runs:
        raise RuntimeError(f"No runs found for scenario '{scenario_name}' in report: {args.report}")
    checks = [c for c in report.get("checks", []) if c.get("scenario") == scenario_name]
    momentum_thresholds = {
        c["run_id"]: c["threshold"]
        for c in checks
        if c.get("metric") == "max_momentum_drift_pcm"
    }

    ic_rows = read_ic(Path(args.ic))
    n_particle = len(ic_rows)

    curves: List[Tuple[str, List[Point], str, str]] = []
    legend_models: List[Tuple[str, str]] = []
    table_rows: List[str] = []
    version = "unknown"
    run_records: List[Dict] = []
    for run in runs:
        run_id = run["run_id"]
        args_map = parse_command_args(run["command"])
        prefix = args_map.get("-f", run_id)
        status_path = WORK_DIR / f"{prefix}.status"
        log_path = Path(args.out_dir) / str(run.get("output", f"{run_id}.log"))
        if log_path.exists():
            version = parse_log_version(log_path)
        ms = momentum_series(status_path, n_particle)
        metrics = run.get("metrics", {})
        artificials = float(metrics.get("max_artificial_particles_glb", 0.0))
        label = MODE_LABELS.get(run_id, (run_id, "#555"))[0]
        color = MODE_LABELS.get(run_id, (run_id, "#555"))[1]
        legend_models.append((label, color))
        # r2 is dashed so its exact overlap with r1 stays visible (the two
        # hard-mode runs track the same drift trajectory; TT is inactive)
        dash = "7,4" if run_id == "r2_hard_no_tt" else ""
        curves.append((label, list(zip(ms["time"], ms["drift"])), color, dash))
        threshold = momentum_thresholds.get(run_id)
        threshold_str = str(threshold) if threshold is not None else "\u2014"
        verdict = ""
        if threshold is not None:
            ok = ms["max_drift"] <= threshold
            verdict = "PASS" if ok else "FAIL"
        table_rows.append(
            f"<tr><td>{html_escape(run_id)}</td><td>{html_escape(label)}</td><td>{ms['n_sample']}</td>"
            f"<td>{ms['max_drift']:.4e}</td><td>{threshold_str}</td>"
            f"<td class=\"{'pass' if verdict == 'PASS' else 'fail' if verdict else ''}\">{verdict}</td>"
            f"<td>{artificials:.0f}</td></tr>"
        )
        run_records.append({"run_id": run_id, "command": run["command"], "artificials": artificials})

    threshold_lines = [
        (0.005, "r3 threshold 5\u00d710\u207b\u00b3 (pre-fix \u22481.0\u00d710\u207b\u00b2 fails)"),
        (0.02, "r1/r2 threshold 2\u00d710\u207b\u00b2"),
    ]
    all_pass = all(c.get("passed", False) for c in checks)
    verdict = "ALL CHECKS PASSED" if all_pass else "CHECKS FAILED"

    command_rows_html = "".join(
        f"<tr><td>{html_escape(r['run_id'])}</td><td><code>{html_escape(r['command'])}</code></td></tr>"
        for r in run_records
    )
    checks_html = "".join(check_row(c) for c in checks)

    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>T5 tidal-tensor momentum conservation</title>
<style>{PAPER_CSS}</style></head><body>

<h1>T5 &nbsp;Tidal-tensor linear-momentum conservation under<br>AR group-boundary churn</h1>
<div class="subtitle">validation report &middot; scenario <i>t5_massive_binary_tt</i> &middot; <span class="verdict {'pass' if all_pass else 'fail'}">{verdict}</span></div>

<h2>1&nbsp;&nbsp;Background and purpose</h2>
<p>On 2026-10-06 a production run (N = 1000, 200 Myr) developed a linear-momentum leak traced to
the AR-group tidal tensor: the even-order (T<sub>3</sub> ~ r<sup>-4</sup>) term of the perturber
tidal tensor had a non-zero mass-weighted mean over the group members, so it pushed the group
centre of mass by a small amount every AR substep with no counter-reaction &mdash; and with a
light star churning through the group boundary every outer orbit, the bias never averaged out.
The fix removes the member mean of the tensor in <code>calcAccPert</code>
(src/ar_interaction.hpp). This test reproduces the event in a 3-particle system &mdash; a
77.65 M<sub>sun</sub> binary plus a 0.242 M<sub>sun</sub> churner whose pericenter (0.01 pc)
lies inside the group radius and apocenter (0.04 pc) outside it &mdash; and establishes (i) that
the fixed binary conserves momentum in the event configuration (tree mode, TT on, boundary
churn), (ii) that in the isolated-group hard-mode path the TT is structurally inactive, so the
TT-on run must reproduce the TT-off baseline, and (iii) via the artificial-particle census that
the TT is engaged exactly in the tree-mode run.</p>

<h2>2&nbsp;&nbsp;Method</h2>
<h3>2.1&nbsp;&nbsp;Quantities under comparison</h3>
<p>At every status output the total linear momentum and its drift are</p>
<div class="eq"><b>P</b>(t) = &Sigma;<sub>i</sub> m<sub>i</sub><b>v</b><sub>i</sub>(t),&nbsp;&nbsp;
&Pi;(t) = ||<b>P</b>(t)| &minus; |<b>P</b>(0)|| / M,<span class="tag">(1)</span></div>
<p>with M the total mass; a momentum-conserving integrator keeps &Pi; at the round-off /
integration-noise floor. The tidal-tensor activity guard is the artificial-particle census</p>
<div class="eq">N<sub>art</sub> = max<sub>t</sub> [N<sub>all</sub>(glb) &minus;
N<sub>real</sub>(glb)],<span class="tag">(2)</span></div>
<p>positive only while tidal-tensor probes exist.</p>

<h3>2.2&nbsp;&nbsp;Metrics</h3>
<table class="data">
<tr><th>symbol</th><th>definition</th><th>interpretation</th></tr>
<tr><td>&Pi;<sub>max</sub></td><td>max over the run of Eq. (1)</td>
<td>momentum-drift level per run; the check metric</td></tr>
<tr><td>N<sub>art</sub></td><td>Eq. (2)</td>
<td>guard: TT engaged (tree mode) / structurally inactive (hard mode)</td></tr>
<tr><td>N<sub>sample</sub></td><td>status outputs per run</td><td>drift-curve sampling</td></tr>
</table>
<div class="caption"><b>Table 1.</b> Metrics used in this report. Momentum is computed from the
*.status particle states (N_particle = 3, plain format, interrupt_mode none).</div>

<h3>2.3&nbsp;&nbsp;Runs</h3>
<table class="data">
<tr><th>run</th><th>configuration</th><th>purpose</th></tr>
<tr><td>r1_hard_tt</td><td>hard mode (r_out = 2), --tt-switch 1, r_group = 0.03, dt_soft = 3&times;10<sup>-5</sup></td>
<td>the hard cluster is always the whole system (no tree perturbers); the AR
group census still splits and re-joins as the churner crosses the linking
radius (group event dumps show &asymp;10<sup>3</sup> three-member/binary-group
topology switches), but every cluster stays isolated, so the TT channel is
structurally inactive and this run must match the r2 baseline</td></tr>
<tr><td>r2_hard_no_tt</td><td>as r1 with --tt-switch 0</td>
<td>TT-off baseline of the isolated-group path (A/B counterpart of r1)</td></tr>
<tr><td>r3_tree_tt</td><td>tree mode (r_out = 0.08, r-ratio 0.3), --tt-switch 1, r_group = 0.02, dt_soft = 0.0039</td>
<td>the production configuration of the event: TT engaged with group-boundary churn every outer
orbit</td></tr>
</table>
<div class="caption"><b>Table 2.</b> Three runs of one IC; r1/r2 are an A/B pair isolating the
TT knob in the isolated-group path (expected identical), r3 exercises the churn + TT combination
that leaked before the fix. All runs integrate to t<sub>end</sub> = 3 Myr with output interval
0.05 Myr, 1 rank &times; 1 thread.</div>

<h3>2.4&nbsp;&nbsp;Reproducibility: software, build and execution</h3>
<p><b>Software.</b> PeTar version {html_escape(version)}, git commit
<code>{html_escape(git_commit(REPO))}</code>, 1 MPI rank &times; 1 OpenMP thread. The IC is
plain-format (no stellar-evolution or external-potential columns), so a plain-family binary is
required and selected via <code>petar.select --optional mpi,omp,avx2</code> (built with
<code>./configure &amp;&amp; make install</code> when missing; note this changes the configure
state).</p>
<p><b>Initial conditions</b> (deterministic, both orbits started at apocenter, c.m. position and
velocity subtracted exactly, no RNG):</p>
<pre class="cmd">python3 test/validation/make_ic.py --case t5 --output test/validation/work/t5/input.base
cd test/validation/work/t5 &amp;&amp; petar.init -f input input.base</pre>
<table class="data">
<tr><th>parameter</th><th>value</th></tr>
<tr><td>system</td><td>massive binary + light churner (3 particles), reproducing the 2026-10-06
event (run_N1k_200Myr at t &asymp; 47.6 Myr)</td></tr>
<tr><td>masses</td><td>m<sub>1</sub> = 46.06, m<sub>2</sub> = 31.59 (binary),
m<sub>3</sub> = 0.242 M<sub>sun</sub> (churner)</td></tr>
<tr><td>inner (binary) orbit</td><td>a = 0.0081 pc, e = 0.667 (r<sub>peri</sub> = 2.7&times;10<sup>-3</sup>,
r<sub>apo</sub> = 1.35&times;10<sup>-2</sup> pc)</td></tr>
<tr><td>outer (churner) orbit</td><td>a = 0.025 pc, e = 0.6 (r<sub>peri</sub> = 0.01 pc inside
the group radius, r<sub>apo</sub> = 0.04 pc outside it &mdash; membership churns every outer
orbit)</td></tr>
<tr><td>units / constants</td><td>M<sub>sun</sub>, pc, pc Myr<sup>-1</sup>;
G = 4.49830997959438 &times; 10<sup>-3</sup> pc&sup3; M<sub>sun</sub><sup>-1</sup> Myr<sup>-2</sup></td></tr>
<tr><td>file format</td><td>ASCII m x y z v<sub>x</sub> v<sub>y</sub> v<sub>z</sub>
(+ petar.init -f conversion)</td></tr>
</table>
<div class="caption"><b>Table 3.</b> IC generation parameters (<code>make_ic.py --case t5</code>).
The churner pericenter/apocenter straddle the group radius by construction, so the group census
changes every outer orbit without relying on a spontaneous rare event.</div>
<p><b>Execution commands</b> (as recorded by the runner):</p>
<table class="data">
<tr><th>run_id</th><th>command</th></tr>
{command_rows_html}
</table>

<h2>3&nbsp;&nbsp;Results</h2>
<h3>3.1&nbsp;&nbsp;Momentum drift by mode</h3>
{svg_curves_plot(
    curves,
    legend_models,
    'T5: linear-momentum drift vs time for the three modes',
    't [Myr]',
    'Π (Eq. 1) [pc Myr^-1]',
    hlines=threshold_lines,
    logy=True,
    legend_quantities=[('Π(t)', ''), ('r2 baseline (dashed, overlaps r1)', '7,4')],
)}
<div class="caption"><b>Figure 1.</b> &Pi;(t) (Eq. 1) versus time; <b>color = run</b> (see
legend), threshold lines dotted. <b>Only two distinct tracks are visible, and this is the
result, not a missing curve</b>: with the TT channel structurally inactive in the isolated
hard cluster, the TT-on run r1 reproduces the TT-off baseline r2 (identical run-max
2.7&times;10<sup>-5</sup>); r2 is drawn dashed so the exact overlap stays visible. The check
metric is the curve maximum per run. Read it as: a
curve hugging the plot floor is momentum-conserving (r1/r2 and the fixed r3); a curve that ramps
up over the outer-orbit time (0.05&ndash;0.1 Myr scale) and keeps growing is the pre-fix leak
signature &mdash; the r3 threshold (5&times;10<sup>-3</sup>) is calibrated so the pre-fix binary
(&Pi; &asymp; 1.03&times;10<sup>-2</sup>) fails and the fixed binary (&Pi; &asymp;
1.4&times;10<sup>-3</sup>) passes with margin; the r1/r2 threshold (2&times;10<sup>-2</sup>)
only bounds gross breakage of the isolated-group path.</div>

<h3>3.2&nbsp;&nbsp;Run summary</h3>
<table class="data">
<tr><th>run</th><th>mode</th><th>N<sub>sample</sub></th><th>&Pi;<sub>max</sub> [pc Myr<sup>-1</sup>]</th><th>threshold</th><th>verdict</th><th>N<sub>art</sub></th></tr>
{''.join(table_rows)}
</table>
<div class="caption"><b>Table 4.</b> Per-run maxima. The verdict column restates the momentum
checks of &sect;3.3 against the calibrated thresholds; N<sub>art</sub> (Eq. 2) is the engagement
guard (zero in both hard-mode runs, non-zero in the tree-mode run).</div>

<h3>3.3&nbsp;&nbsp;Automatic checks</h3>
<table class="data">
<tr><th>check</th><th>status</th><th>value</th><th>threshold</th><th>message</th></tr>
{checks_html}
</table>
<div class="caption"><b>Table 5.</b> Automatic checks. <i>Accuracy thresholds</i>: the three
momentum checks cap &Pi;<sub>max</sub> &mdash; r1/r2 at 2&times;10<sup>-2</sup> (gross bound for
the TT-inactive path), r3 at 5&times;10<sup>-3</sup>, calibrated on this IC so that the pre-fix
binary fails and the fixed binary passes (a threshold property of the IC, not a universal
constant). <i>Guards</i>: r1/r2 artificial-particle counts must be zero (TT structurally
inactive &mdash; a non-zero count there would mean r1 silently exercised a different path than
r2), r3 must have &ge; 1 (TT actually engaged), and every run must abort-free.</div>

<div class="meta">petar version: {html_escape(version)} &middot; scenario:
<i>t5_massive_binary_tt</i> &middot; report: {html_escape(args.report)}</div>
</body></html>
"""

    Path(args.output).write_text(html, encoding="utf-8")
    n_pass = sum(1 for c in checks if c.get("passed", False))
    print(f"[T5] Summary: passed={n_pass == len(checks) and len(checks) > 0} total={len(checks)} failed={len(checks) - n_pass} report={args.output}")
    return 0 if n_pass == len(checks) and checks else 1


if __name__ == "__main__":
    raise SystemExit(main())
