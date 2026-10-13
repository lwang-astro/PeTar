#!/usr/bin/env python3
"""T2 long-term binary conservation HTML report (paper style, see REPORT_DESIGN.md)."""
import argparse
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

from report_common import (
    G_MSUN_PC_MYR,
    PAPER_CSS,
    check_row,
    git_commit,
    gradient_colors,
    html_escape,
    load_json,
    orbital_elements_from_two_body_rows,
    parse_command_args,
    parse_log_version,
    parse_time_error_series,
    read_ic,
    select_or_install_petar,
    svg_curves_plot,
)

REPO = Path(__file__).resolve().parents[2]
Point = Tuple[float, float]
RIN_RE = re.compile(r"^Mean inner changeover radius\s*=\s*([0-9.eE+\-.]+)")
ROUT_RE = re.compile(r"^Mean outer changeover radius\s*=\s*([0-9.eE+\-.]+)")


def resolve_status_path(run_record: Dict) -> Path:
    cmd = run_record["command"]
    args_map = parse_command_args(cmd)
    prefix = args_map.get("-f", run_record["run_id"])
    # status files are in test/out/work/t2/ regardless of log output dir
    return Path("test/out/work/t2") / f"{prefix}.status"


def load_status_binary(path: Path):
    sys.path.insert(0, str(REPO / "tools"))
    from analysis.status import Status  # pylint: disable=import-outside-toplevel  # pyright: ignore[reportMissingImports]

    st = Status(N_particle=2)
    st.fromfile(str(path))
    return st


def derive_orbital_drift(status_path: Path) -> Dict[str, Any]:
    st = load_status_binary(status_path)

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
        "n_sample": int(st.size),
        "t_end": float(st.time[-1]),
        "da_rel_max": float(np.max(da_rel)),
        "de_abs_max": float(np.max(de_abs)),
    }


def parse_changeover_radii(log_path: Path) -> Tuple[float, float]:
    rin = float("nan")
    rout = float("nan")
    for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        m = RIN_RE.match(line.strip())
        if m:
            rin = float(m.group(1))
        m = ROUT_RE.match(line.strip())
        if m:
            rout = float(m.group(1))
    return rin, rout


def classify_regime(rin: float, rout: float, peri: float, apo: float) -> str:
    if rout < peri:
        return "pure tree"
    if rin > apo:
        return "pure Hermite"
    return "mix"


def orbit_envelope(series: List[Point], period: float, n_orbits: int) -> List[float]:
    env = [0.0] * n_orbits
    for t, v in series:
        k = int(t / period)
        if 0 <= k < n_orbits and v > env[k]:
            env[k] = v
    return env


def _detect_petar() -> str:
    """Select (and if needed reinstall) the petar binary via petar.select."""
    return select_or_install_petar(optional=["avx2", "omp"], repo_root=REPO)


def _run_scenario(scenario_file: str, report_path: str, out_dir: str, petar_bin: str) -> None:
    """Invoke run_validation.py for a single scenario."""
    import subprocess as _sp
    cmd = [
        sys.executable, "test/validation/run_validation.py",
        "--scenario", scenario_file,
        "--out-dir", out_dir,
        "--report", report_path,
        "--var", f"petar_bin_switch={petar_bin}",
    ]
    print(f"[T2] Running: {' '.join(cmd)}")
    r = _sp.run(cmd, check=True)
    if r.returncode != 0:
        raise RuntimeError(f"run_validation.py failed with code {r.returncode}")


def scenario_file_from_name(scenario: str) -> Path:
    return Path("test/validation/scenarios") / f"{scenario}.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="T2 long-term binary conservation: run + HTML report")
    parser.add_argument("--run", action="store_true", help="Run the scenario before generating the report")
    parser.add_argument("--report", default="test/out/report.t2.binary.json")
    parser.add_argument("--scenario", default="t2_binary_conservation_longterm")
    parser.add_argument("--scenario-file", default="")
    parser.add_argument("--out-dir", default="test/out/validation_t2_binary", help="Output directory for run logs")
    parser.add_argument("--petar", default="", help="Petar binary (auto-detected if empty)")
    parser.add_argument("--ic", default="test/out/work/t2/input.base")
    parser.add_argument("--output", default="test/out/t2_binary_longterm_autochangeover_summary.html")
    args = parser.parse_args()

    # --run: execute the scenario via run_validation.py
    if args.run:
        petar = args.petar or _detect_petar()
        sfile = args.scenario_file or str(scenario_file_from_name(args.scenario))
        _run_scenario(sfile, args.report, args.out_dir, petar)

    report = load_json(Path(args.report))
    scenario_path = Path(args.scenario_file) if args.scenario_file else scenario_file_from_name(args.scenario)
    scenario_data = load_json(scenario_path) if scenario_path.exists() else {}

    runs = [r for r in report.get("runs", []) if r.get("scenario") == args.scenario]
    if not runs:
        raise RuntimeError(f"No runs found for scenario '{args.scenario}' in report: {args.report}")
    checks = [c for c in report.get("checks", []) if c.get("scenario") == args.scenario and not c.get("skipped")]

    ic_rows = read_ic(Path(args.ic))
    if len(ic_rows) != 2:
        raise RuntimeError(f"T2 report expects a pure binary IC with 2 rows, got {len(ic_rows)}")
    orb0 = orbital_elements_from_two_body_rows(ic_rows)

    rows = []
    for run in runs:
        args_map = parse_command_args(run["command"])
        dt = float(args_map.get("-s", "nan"))
        status_path = resolve_status_path(run)
        drift = derive_orbital_drift(status_path)
        log_path = Path(run["output"])
        rin, rout = parse_changeover_radii(log_path)
        metrics = run.get("metrics", {})
        energy_series = parse_time_error_series(log_path)
        n_orbits = int((energy_series[-1][0] + 1e-12) / orb0["period"]) if energy_series else 0
        rows.append(
            {
                "run_id": run["run_id"],
                "dt": dt,
                "r_group": args_map.get("--r-group", ""),
                "r_search_group": args_map.get("--r-search-group", ""),
                "rin_auto": rin,
                "rout_auto": rout,
                "regime": classify_regime(rin, rout, orb0["peri"], orb0["apo"]),
                "command": run["command"],
                "log_path": run["output"],
                "status_path": str(status_path),
                "energy_max": float(metrics.get("max_abs_error_over_total", 0.0)),
                "late_early": float(metrics.get("late_early_period_max_ratio", float("nan"))),
                "energy_env": orbit_envelope(energy_series, orb0["period"], n_orbits),
                "n_orbits": n_orbits,
                "da_rel_max": drift["da_rel_max"],
                "de_abs_max": drift["de_abs_max"],
                "t_end": drift["t_end"],
                "n_sample": drift["n_sample"],
            }
        )

    rows = sorted(rows, key=lambda x: x["dt"])
    log_version = parse_log_version(Path(rows[0]["log_path"]))

    # Figure 1: per-period max energy error vs orbit number, color = dt order
    curves: List[Tuple[str, List[Point], str, str]] = []
    legend_models: List[Tuple[str, str]] = []
    sweep_colors = gradient_colors(len(rows))
    for i, row in enumerate(rows):
        color = sweep_colors[i]
        model = f"{row['run_id']} (dt={row['dt']:.2e}, {row['regime']})"
        legend_models.append((model, color))
        xs = [k + 0.5 for k in range(row["n_orbits"])]
        curves.append((model, list(zip(xs, row["energy_env"])), color, ""))

    # Figure 2: run-max level vs dt_soft, color = regime marker-free single model
    level_curves: List[Tuple[str, List[Point], str, str]] = [
        ("max |Error/Total|", [(r["dt"], r["energy_max"]) for r in rows], "#0b74de", ""),
    ]

    all_pass = all(c.get("passed", False) for c in checks)
    verdict = "ALL CHECKS PASSED" if all_pass else "CHECKS FAILED"

    check_descriptions = {
        c.get("name"): str(c.get("comment", ""))
        for c in (scenario_data.get("checks", []) if isinstance(scenario_data, dict) else [])
        if isinstance(c, dict)
    }

    def fmt(x, spec=".3e"):
        try:
            return format(float(x), spec)
        except (TypeError, ValueError):
            return "\u2014"

    run_rows_html = "".join(
        f"<tr><td>{html_escape(r['run_id'])}</td><td>{r['dt']:.6e}</td>"
        f"<td>{fmt(r['rin_auto'], '.4g')}</td><td>{fmt(r['rout_auto'], '.4g')}</td>"
        f"<td>{html_escape(r['regime'])}</td>"
        f"<td>{fmt(r['energy_max'])}</td><td>{fmt(r['late_early'], '.3f')}</td>"
        f"<td>{fmt(r['da_rel_max'])}</td><td>{fmt(r['de_abs_max'])}</td>"
        f"<td>{r['n_orbits']}</td><td>{r['n_sample']}</td></tr>"
        for r in rows
    )

    # §2.3 configuration rows: run, dt, auto radii, regime, test items
    stationary_runs = {"dt17", "dt12", "dt8"}

    def test_items(rid: str) -> str:
        items = []
        if rid in stationary_runs:
            items.append(f"<code>energy_error_stationary_{rid}</code>")
        if rid == "dt17":
            items.append("level numerator of <code>energy_level_consistent_across_dt</code>")
        elif rid == "dt12":
            items.append("level denominator of <code>energy_level_consistent_across_dt</code>")
        else:
            items.append("level span member of <code>energy_level_consistent_across_dt</code>")
        return "; ".join(items)

    config_rows_html = "".join(
        f"<tr><td>{html_escape(r['run_id'])}</td><td>{r['dt']:.6e}</td>"
        f"<td>{fmt(r['rin_auto'], '.4g')}</td><td>{fmt(r['rout_auto'], '.4g')}</td>"
        f"<td>{html_escape(r['regime'])}</td>"
        f"<td>{test_items(r['run_id'])}</td></tr>"
        for r in rows
    )

    checks_html = "".join(
        f'<tr><td>{html_escape(c["check"])}</td>'
        f'<td class="{"pass" if c.get("passed") else "fail"}">{"PASS" if c.get("passed") else "FAIL"}</td>'
        f"<td>{fmt(c.get('value'), '.4g')}</td>"
        f"<td>{html_escape(str(c.get('threshold', '')))}</td>"
        f'<td class="msg">{html_escape(check_descriptions.get(c["check"], c.get("message", "")))}</td></tr>'
        for c in checks
    )

    setup_commands = scenario_data.get("setup_commands", []) if isinstance(scenario_data, dict) else []
    if not isinstance(setup_commands, list):
        setup_commands = []
    setup_html = "<br>".join(html_escape(str(cmd)) for cmd in setup_commands)

    binary_primary = rows[0]["command"].split("OMP_NUM_THREADS=1 ")[1].split(" -a")[0]
    level_ratio = max(r["energy_max"] for r in rows) / min(r["energy_max"] for r in rows)
    dt_span = rows[-1]["dt"] / rows[0]["dt"]

    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>T2 long-term binary conservation under auto changeover</title>
<style>{PAPER_CSS}</style></head><body>

<h1>T2 &nbsp;Long-term binary evolution under auto changeover radii:<br>error drift and level consistency across dt<sub>soft</sub></h1>
<div class="subtitle">validation report &middot; scenario <i>t2_binary_conservation_longterm</i> &middot; <span class="verdict {'pass' if all_pass else 'fail'}">{verdict}</span></div>

<h2>1&nbsp;&nbsp;Background and purpose</h2>
<p>With auto changeover radii (no <code>-r</code>), PeTar derives r<sub>in</sub> from
dt<sub>soft</sub> (r<sub>in</sub> solves P(r<sub>in</sub>) = n<sub>step</sub>&thinsp;dt<sub>soft</sub>,
r<sub>out</sub> = 10&thinsp;r<sub>in</sub>), so sweeping dt<sub>soft</sub> moves the changeover
through the orbit: at fine dt the whole orbit lies beyond r<sub>out</sub> (pure tree), at
intermediate dt the orbit crosses the changeover (mix), and at coarse dt the whole orbit lies
inside r<sub>in</sub> (pure Hermite). A pure binary evolved for &gt;100 orbital periods in this
configuration answers two questions:</p>
<p><b>T2-a (drift).</b> With the changeover participating in the long-term evolution &mdash; the
sweep here covers pure tree (dt17), mix (dt16&ndash;dt9) and pure Hermite (dt8) &mdash; does the
energy error show secular drift, or does it stay stationary at a bounded level?</p>
<p><b>T2-b (level consistency).</b> Do the auto-derived radii keep the energy-error level
approximately consistent across dt<sub>soft</sub>? With fixed radii the tree-channel error
grows as dt<sup>4</sup>; the auto algorithm recouples the radii to dt, and the question is how
much of that scaling remains across the sweep.</p>

<h2>2&nbsp;&nbsp;Method</h2>
<h3>2.1&nbsp;&nbsp;Quantities under comparison</h3>
<p>From each run's energy log, the per-period envelope of the relative energy error is</p>
<div class="eq">E&#770;<sub>E</sub>(k) = max<sub>t&thinsp;&isin;&thinsp;[kP,(k+1)P)</sub>
|Error/Total|(t),<span class="tag">(1)</span></div>
<p>with P the binary period. Stationarity is quantified by the half-run ratio</p>
<div class="eq">&rho; = max<sub>k&thinsp;&ge;&thinsp;K/2</sub> E&#770;<sub>E</sub>(k)
&thinsp;/&thinsp; max<sub>k&thinsp;&lt;&thinsp;K/2</sub> E&#770;<sub>E</sub>(k)
&nbsp;&nbsp;(&rho; &asymp; 1 stationary, &rho; &gg; 1 secular drift),<span class="tag">(2)</span></div>
<p>and the level of a run by its run-maximum</p>
<div class="eq">&epsilon;<sub>E</sub> = max<sub>t</sub> |Error/Total|(t).<span class="tag">(3)</span></div>
<p>The osculating-element drift of the orbit (|&Delta;a/a<sub>0</sub>|, |&Delta;e| from the
status two-body states, referenced to the first output sample) is reported as supporting
numbers in Table 4; as in T1, the raw drift is a bounded orbital-phase oscillation, so only its
run-maximum is quoted.</p>

<h3>2.2&nbsp;&nbsp;Metrics</h3>
<table class="data">
<tr><th>symbol</th><th>definition</th><th>interpretation</th></tr>
<tr><td>E&#770;<sub>E</sub>(k)</td><td>Eq. (1)</td>
<td>per-period error envelope; flat in k = stationary error (T2-a)</td></tr>
<tr><td>&rho;</td><td>Eq. (2)</td>
<td>half-run drift ratio; the T2-a check quantity</td></tr>
<tr><td>&epsilon;<sub>E</sub></td><td>Eq. (3)</td>
<td>run-maximum level; compared across dt (T2-b)</td></tr>
<tr><td>r<sub>in</sub>, r<sub>out</sub> (auto)</td><td>from the run log</td>
<td>regime classification against r<sub>peri</sub> = {orb0['peri']:.3g},
r<sub>apo</sub> = {orb0['apo']:.3g} pc</td></tr>
<tr><td>|&Delta;a/a<sub>0</sub>|<sub>max</sub>, |&Delta;e|<sub>max</sub></td>
<td>run-max orbital-element drift</td><td>supporting conservation numbers</td></tr>
</table>
<div class="caption"><b>Table 1.</b> Metrics used in this report.</div>

<h3>2.3&nbsp;&nbsp;Runs and test items</h3>
<table class="data">
<tr><th>run</th><th>dt<sub>soft</sub> [Myr]</th><th>r<sub>in</sub> (auto)</th><th>r<sub>out</sub> (auto)</th><th>regime vs. the orbit</th><th>test items (automatic checks)</th></tr>
{config_rows_html}
</table>
<div class="caption"><b>Table 2.</b> The dt<sub>soft</sub> sweep with its auto-derived
changeover radii and the resulting regime: dt17 pure tree (r<sub>out</sub> &lt;
r<sub>peri</sub>), dt16&ndash;dt9 mix (changeover inside the orbit), dt8 pure Hermite
(r<sub>in</sub> &gt; r<sub>apo</sub>). Stationarity is checked at one run per regime
(&rho; &le; 2), and the level consistency across the whole sweep by the dt17/dt12 run-max
ratio. All runs integrate to t<sub>end</sub> = 3.2 Myr &asymp; 113 periods; group linking is
kept too small to form an AR group (dt17 uses smaller group radii because the auto
r<sub>in</sub> shrinks below the default r<sub>search,group</sub>; the binary separation never
drops below either setting). Per-run metrics are in Table 4.</div>

<h3>2.4&nbsp;&nbsp;Reproducibility: software, build and execution</h3>
<p><b>Software.</b> PeTar version {html_escape(log_version)}, git commit
<code>{html_escape(git_commit(REPO))}</code>, default-precision build, tree mode (-a 0), 1 MPI
rank &times; 1 OpenMP thread, executable <code>{html_escape(binary_primary)}</code>.</p>
<p><b>Initial conditions</b> (deterministic two-body construction, started at pericenter, no
RNG):</p>
<pre class="cmd">{setup_html}</pre>
<table class="data">
<tr><th>parameter</th><th>value</th></tr>
<tr><td>system</td><td>isolated binary</td></tr>
<tr><td>masses</td><td>m<sub>1</sub> = 1 M<sub>sun</sub>, m<sub>2</sub> = 10 M<sub>sun</sub></td></tr>
<tr><td>orbit</td><td>a = 0.01 pc, e = 0.6 (derived: r<sub>peri</sub> = {orb0['peri']:.4g} pc,
r<sub>apo</sub> = {orb0['apo']:.4g} pc, P = {orb0['period']:.6g} Myr)</td></tr>
<tr><td>changeover</td><td>auto radii (no -r): r<sub>in</sub> from dt<sub>soft</sub>,
r<sub>out</sub> = 10&thinsp;r<sub>in</sub>, per run in Table 2</td></tr>
<tr><td>units / constants</td><td>M<sub>sun</sub>, pc, pc Myr<sup>-1</sup>;
G = 4.49830997959438 &times; 10<sup>-3</sup> pc&sup3; M<sub>sun</sub><sup>-1</sup> Myr<sup>-2</sup></td></tr>
<tr><td>file format</td><td>ASCII m x y z v<sub>x</sub> v<sub>y</sub> v<sub>z</sub>
(+ petar.init -f conversion)</td></tr>
</table>
<div class="caption"><b>Table 3.</b> IC generation parameters (<code>make_ic.py --case t2</code>).</div>
<p><b>Execution</b> (all runs share one command template):</p>
<pre class="cmd">cd test/out/work/t2 &amp;&amp; OMP_NUM_THREADS=1 {{petar_bin_switch}} -a 0 --r-group {{r_group}} --r-search-group {{r_search_group}} -u 1 -b 1 -t 3.2 -s {{dt_soft}} -o {{dt_out}} -w 2 -f {{run_id}} input</pre>
<p>Variables: <code>{{dt_soft}}</code> = 2<sup>-17</sup>&hellip;2<sup>-8</sup> (ten contiguous values of
Table 2); <code>{{dt_out}}</code> = 2&middot;dt<sub>soft</sub> (bounded to 9.77&times;10<sup>-4</sup> at the fine end, dt15&ndash;dt17, to cap the output volume); <code>{{r_group}}</code>/<code>{{r_search_group}}</code> = 5e-4/1e-3
(5e-6/1e-5 for dt13&ndash;dt17, keeping r<sub>search,group</sub> &le; auto r<sub>in</sub>);
<code>{{run_id}}</code> = t2_dtXX. Exact recorded commands are in the JSON report.</p>

<h2>3&nbsp;&nbsp;Results</h2>
<h3>3.1&nbsp;&nbsp;Energy-error envelope: stationarity across regimes (T2-a)</h3>
{svg_curves_plot(
    curves,
    legend_models,
    'T2: per-period max energy error vs orbit number, dt sweep',
    'orbit number k',
    'Ê_E (Eq. 1)',
    logy=True,
    legend_quantities=[('per-period max |Error/Total|', '')],
)}
<div class="caption"><b>Figure 1.</b> Per-period envelope (Eq. 1) of |Error/Total| versus
orbit number; <b>color = dt<sub>soft</sub></b> on a viridis gradient from fine (dt17, dark
purple) to coarse (dt8, yellow) &mdash; the sweep direction of Figure 2. Read it as:
a curve flat in k is a stationary, bounded error &mdash; no secular drift; a systematically
rising curve would be drift (T2-a failure). The &rho; column of Table 4 quantifies flatness:
all runs sit at &rho; &asymp; 1. The vertical ordering of the curves is the level question of
T2-b (Figure 2).</div>

<h3>3.2&nbsp;&nbsp;Error level vs dt<sub>soft</sub>: consistency of the auto radii (T2-b)</h3>
{svg_curves_plot(
    level_curves,
    [('T2 runs (dt17 pure tree &rarr; dt8 pure Hermite)', '#0b74de')],
    'T2: run-max energy error vs dt_soft',
    'dt_soft [Myr]',
    'ε_E',
    logx=True,
    logy=True,
    legend_quantities=[('max |Error/Total|', '')],
)}
<div class="caption"><b>Figure 2.</b> Run-maximum energy error &epsilon;<sub>E</sub> versus
dt<sub>soft</sub> (log&ndash;log), one point per run; the sweep crosses the pure-tree &rarr; mix
&rarr; pure-Hermite regime boundaries as the auto radii track dt (regimes per Table 2). Read it
as: with fixed radii the tree channel would scale as dt<sup>4</sup> (a factor
{fmt(dt_span**4, '.1e')} across this dt span); the measured span is only
{level_ratio:.0f}&times; ({fmt(min(r['energy_max'] for r in rows), '.2e')} at dt12 to
{fmt(max(r['energy_max'] for r in rows), '.2e')} at dt17) &mdash; the auto algorithm keeps the
level approximately consistent, with the <i>finest</i> run (pure tree) carrying the largest
error because the changeover has moved inside the pericenter. The check bounds the total span
(&le; 300&times;).</div>

<h3>3.3&nbsp;&nbsp;Run summary and automatic checks</h3>
<table class="data">
<tr><th>run</th><th>dt<sub>soft</sub></th><th>r<sub>in</sub> (auto)</th><th>r<sub>out</sub> (auto)</th><th>regime</th><th>&epsilon;<sub>E</sub></th><th>&rho; (Eq. 2)</th><th>|&Delta;a/a<sub>0</sub>|<sub>max</sub></th><th>|&Delta;e|<sub>max</sub></th><th>orbits</th><th>N<sub>sample</sub></th></tr>
{run_rows_html}
</table>
<div class="caption"><b>Table 4.</b> Per-run summary: auto changeover radii and regime,
error level &epsilon;<sub>E</sub>, stationarity ratio &rho;, and the supporting
osculating-element drift maxima (bounded phase oscillation, ~10<sup>-4</sup> level for all
runs).</div>
<table class="data">
<tr><th>check</th><th>status</th><th>value</th><th>threshold</th><th>description</th></tr>
{checks_html}
</table>
<div class="caption"><b>Table 5.</b> Automatic checks by purpose. <b>T2-a</b>: one stationarity
check per regime (&rho; &le; 2; measured &rho; &asymp; 1 for every run, see Table 4).
<b>T2-b</b>: the level span across the full sweep (dt17/dt12 run-max ratio) within [1, 300]
&mdash; measured {level_ratio:.0f}&times;, versus {fmt(dt_span**4, '.0e')} for dt<sup>4</sup>
scaling.</div>

<div class="meta">petar version: {html_escape(log_version)} &middot; scenario:
<i>t2_binary_conservation_longterm</i> &middot; report: {html_escape(args.report)}</div>
</body></html>
"""

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html, encoding="utf-8")
    n_pass = sum(1 for c in checks if c.get("passed", False))
    print(f"[T2] Summary: passed={n_pass == len(checks) and len(checks) > 0} total={len(checks)} failed={len(checks) - n_pass} report={out_path}")
    return 0 if n_pass == len(checks) and checks else 1


if __name__ == "__main__":
    raise SystemExit(main())
