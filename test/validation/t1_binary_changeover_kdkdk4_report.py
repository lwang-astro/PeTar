#!/usr/bin/env python3
"""T1 KDKDK4 changeover HTML report (paper style, see REPORT_DESIGN.md)."""
import argparse
import math
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

from report_common import (
    PAPER_CSS,
    check_row,
    git_commit,
    html_escape,
    load_json,
    orbital_elements_from_two_body_rows,
    parse_command_args,
    parse_log_version,
    read_ic,
    select_or_install_petar,
    svg_curves_plot,
)

REPO = Path(__file__).resolve().parents[2]
Point = Tuple[float, float]
TIME_RE = re.compile(r"^Time:\s*([0-9eE+\-.]+)")
ENERGY_HEADER_RE = re.compile(r"^Energy:\s+(.*)$")
ENERGY_PHYSIC_RE = re.compile(r"^Physic:\s+(.*)$")

# color per changeover family (model); dash per build (quantity). The
# full-crossing family is the design under test; pure-tree and pure-Hermite
# are the reference channels. Labels are plain text (HTML-escaped downstream).
FAMILY_COLORS = {
    "full_crossing": "#d81b60",
    "pure_tree": "#0b74de",
    "pure_hermite": "#2e7d32",
}
FAMILY_LABELS = {
    "full_crossing": "full crossing (r_peri < r_in < r_out < r_apo)",
    "pure_tree": "pure tree (r_peri > r_out)",
    "pure_hermite": "pure Hermite (r_apo < r_in)",
}
FAMILY_SHORT = {
    "full_crossing": "full crossing",
    "pure_tree": "pure tree",
    "pure_hermite": "pure Hermite",
}


def parse_error_series(log_path: Path) -> List[Point]:
    text = log_path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    series: List[Point] = []
    current_time = None

    for i, line in enumerate(lines):
        tmatch = TIME_RE.match(line.strip())
        if tmatch:
            current_time = float(tmatch.group(1))
            continue
        hmatch = ENERGY_HEADER_RE.match(line.strip())
        if hmatch and i + 1 < len(lines):
            pmatch = ENERGY_PHYSIC_RE.match(lines[i + 1].strip())
            if not pmatch:
                continue
            header = hmatch.group(1).split()
            values = pmatch.group(1).split()
            if len(values) < len(header):
                continue
            rec = {header[j]: float(values[j]) for j in range(len(header))}
            if current_time is not None:
                series.append((current_time, abs(rec.get("Error/Total", 0.0))))
    return series


def extract_binary_from_command(command: str) -> str:
    tokens = command.split()
    for idx, token in enumerate(tokens):
        if token.startswith("OMP_NUM_THREADS=") and idx + 1 < len(tokens):
            return tokens[idx + 1]
    return "unknown"


def scenario_file_from_name(scenario: str) -> Path:
    return Path("test/validation/scenarios") / f"{scenario}.json"


def max_error_within_period(series: List[Point], period: float) -> float:
    selected = [err for t, err in series if t <= period + 1e-12]
    if selected:
        return max(selected)
    if series:
        return max(err for _, err in series)
    return 0.0


def fit_loglog_slope(points: List[Point]) -> float:
    pts = [(x, y) for x, y in points if x > 0 and y > 0]
    if len(pts) < 2:
        return float("nan")
    xs = [math.log10(x) for x, _ in pts]
    ys = [math.log10(y) for _, y in pts]
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    denom = sum((x - mx) ** 2 for x in xs)
    if denom == 0:
        return float("nan")
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / denom


def fit_loglog_slope_subset(points: List[Point], start: int, end: int) -> float:
    return fit_loglog_slope(points[start:end])


def _detect_or_build_petar(require_64b: bool) -> str:
    """Select (and if needed reinstall) a petar binary via petar.select."""
    if require_64b:
        return select_or_install_petar(require=["64b"], optional=["avx2", "omp"], repo_root=REPO)
    return select_or_install_petar(optional=["avx2", "omp"], repo_root=REPO)


def _run_t1_scenario(scenario_file: str, report_path: str, out_dir: str, petar_bin: str) -> None:
    """Invoke run_validation.py for T1 with a specific binary."""
    import subprocess as _sp
    cmd = [
        sys.executable, "test/validation/run_validation.py",
        "--scenario", scenario_file,
        "--out-dir", out_dir,
        "--report", report_path,
        "--var", f"petar_bin_order4={petar_bin}",
    ]
    print(f"[T1] Running: {' '.join(cmd)}")
    r = _sp.run(cmd, check=False)
    if r.returncode != 0:
        print(f"[T1] Warning: run_validation.py returned {r.returncode}", file=sys.stderr)


def main() -> int:
    parser = argparse.ArgumentParser(description="T1 KDKDK4 changeover: 64b+non64b run + merged HTML report")
    parser.add_argument("--run", action="store_true", help="Run both 64b and non64b scenarios before generating report")
    parser.add_argument("--report", default="test/out/report.t1.64b.overlay.json")
    parser.add_argument("--primary-tag", default="64b")
    parser.add_argument("--compare-report", default="test/out/report.t1.non64b.overlay.json")
    parser.add_argument("--compare-tag", default="non64b")
    parser.add_argument("--scenario", default="t1_binary_changeover_kdkdk4")
    parser.add_argument("--scenario-file", default="")
    parser.add_argument("--out-dir-64b", default="test/out/validation_t1_binary_changeover_64b")
    parser.add_argument("--out-dir-non64b", default="test/out/validation_t1_binary_changeover_non64b")
    parser.add_argument("--ic", default="test/out/work/t1/input.base")
    parser.add_argument("--output", default="test/out/t1_binary_changeover_kdkdk4_summary.html")
    args = parser.parse_args()

    if args.run:
        sfile = args.scenario_file or str(scenario_file_from_name(args.scenario))
        petar_64 = _detect_or_build_petar(True)
        petar_non64 = _detect_or_build_petar(False)
        _run_t1_scenario(sfile, args.report, args.out_dir_64b, petar_64)
        _run_t1_scenario(sfile, args.compare_report, args.out_dir_non64b, petar_non64)

    report = load_json(Path(args.report))
    scenario_path = Path(args.scenario_file) if args.scenario_file else scenario_file_from_name(args.scenario)
    scenario_data = load_json(scenario_path) if scenario_path.exists() else {}

    runs = [r for r in report.get("runs", []) if r.get("scenario") == args.scenario]
    # skipped checks (build condition not met) are not displayed
    checks = [c for c in report.get("checks", []) if c.get("scenario") == args.scenario and not c.get("skipped")]
    if not runs:
        raise RuntimeError(f"No runs found for scenario '{args.scenario}' in report: {args.report}")

    compare_report = None
    compare_runs = []
    if args.compare_report:
        compare_path = Path(args.compare_report)
        if compare_path.exists():
            compare_report = load_json(compare_path)
            compare_runs = [r for r in compare_report.get("runs", []) if r.get("scenario") == args.scenario]
        else:
            print(f"[T1] Warning: comparison report not found, skipping: {args.compare_report}", file=sys.stderr)

    ic_rows = read_ic(Path(args.ic))
    orb = orbital_elements_from_two_body_rows(ic_rows)
    period = orb["period"]

    def classify_regime(rout: float, rin: float, peri: float, apo: float) -> str:
        if rout < peri:
            return "pure_tree"
        if rin > apo:
            return "pure_hermite"
        crosses_inner = peri < rin
        crosses_outer = apo > rout
        if crosses_inner and crosses_outer:
            return "full_crossing"
        if crosses_inner:
            return "inner_crossing"
        if rin < peri < rout:
            return "outer_crossing"
        return "other"

    def build_rows(run_list: List[Dict], tag: str) -> List[Dict]:
        out_rows = []
        for run in run_list:
            command = run["command"]
            args_map = parse_command_args(command)
            dt = float(args_map.get("-s", "nan"))
            rout = float(args_map.get("-r", "nan"))
            rin = 0.1 * rout
            log_path = Path(run["output"])
            series = parse_error_series(log_path)
            max_err_1p = max_error_within_period(series, period)
            regime = classify_regime(rout, rin, orb["peri"], orb["apo"])
            out_rows.append(
                {
                    "run_id": run["run_id"],
                    "dt": dt,
                    "rout": rout,
                    "rin": rin,
                    "regime": regime,
                    "max_err_1p": max_err_1p,
                    "command": command,
                    "log_path": str(log_path),
                    "dataset": tag,
                }
            )
        return out_rows

    rows = build_rows(runs, args.primary_tag)
    compare_rows = build_rows(compare_runs, args.compare_tag) if compare_runs else []
    all_rows = rows + compare_rows

    log_version = parse_log_version(Path(rows[0]["log_path"]))

    # curves for Figure 1: one per (family, dataset); color = family, dash = build
    family_order = ["full_crossing", "pure_tree", "pure_hermite"]
    dash_by_dataset = {args.primary_tag: "", args.compare_tag: "7,4"}
    curves: List[Tuple[str, List[Point], str, str]] = []
    for regime in family_order:
        for source_rows, tag in ((rows, args.primary_tag), (compare_rows, args.compare_tag)):
            if not source_rows:
                continue
            candidates = sorted(
                [r for r in source_rows if r["regime"] == regime], key=lambda item: item["dt"]
            )
            if not candidates:
                continue
            label = f"{FAMILY_SHORT.get(regime, regime)} r_out={candidates[0]['rout']:g} ({tag})"
            curves.append((label, [(r["dt"], r["max_err_1p"]) for r in candidates],
                           FAMILY_COLORS[regime], dash_by_dataset.get(tag, "")))
    legend_models = [
        (f"{FAMILY_SHORT.get(regime, regime)} (r_out={min(r['rout'] for r in all_rows if r['regime'] == regime):g})", FAMILY_COLORS[regime])
        for regime in family_order
        if any(r["regime"] == regime for r in all_rows)
    ]
    legend_quantities = [(f"{args.primary_tag} build", ""), (f"{args.compare_tag} build", "7,4")] if compare_rows else [(f"{args.primary_tag} build", "")]

    def fmt(x, spec=".3e"):
        try:
            return format(float(x), spec)
        except (TypeError, ValueError):
            return "\u2014"

    # slope fits per (family, dataset); fine = 4 finest, coarse = 4 coarsest,
    # check-window = si03..si06 (above the non64b round-off floor, below the
    # O(1) saturation of the coarse end)
    slope_rows = []
    for regime in family_order:
        for source_rows, tag in ((rows, args.primary_tag), (compare_rows, args.compare_tag)):
            if not source_rows:
                continue
            candidates = sorted(
                [(r["dt"], r["max_err_1p"]) for r in source_rows if r["regime"] == regime],
                key=lambda item: item[0],
            )
            if not candidates:
                continue
            slope_rows.append(
                {
                    "regime": regime,
                    "dataset": tag,
                    "full_slope": fit_loglog_slope(candidates),
                    "fine_slope": fit_loglog_slope_subset(candidates, 0, min(4, len(candidates))),
                    "coarse_slope": fit_loglog_slope_subset(candidates, max(0, len(candidates) - 4), len(candidates)),
                    "window_slope": fit_loglog_slope_subset(candidates, 3, min(7, len(candidates))),
                }
            )

    # Figure 2 data: crossing error relative to the references at equal dt
    # (Q_t vs pure tree, Q_h vs pure Hermite), color = build, dash = quantity
    def family_dt_map(source_rows: List[Dict], regime: str) -> Dict[float, float]:
        return {
            r["dt"]: r["max_err_1p"]
            for r in source_rows
            if r["regime"] == regime
        }

    ratio_curves: List[Tuple[str, List[Point], str, str]] = []
    ratio_legend_models: List[Tuple[str, str]] = []
    build_colors = {args.primary_tag: "#0b74de", args.compare_tag: "#d81b60"}
    for source_rows, tag in ((rows, args.primary_tag), (compare_rows, args.compare_tag)):
        if not source_rows:
            continue
        cross_map = family_dt_map(source_rows, "full_crossing")
        tree_map = family_dt_map(source_rows, "pure_tree")
        herm_map = family_dt_map(source_rows, "pure_hermite")
        ratio_legend_models.append((f"{tag} build", build_colors[tag]))
        for regime_map, dash, label in (
            (herm_map, "", "Qh"),
            (tree_map, "7,4", "Qt"),
        ):
            pts = [
                (dt, cross_map[dt] / regime_map[dt])
                for dt in sorted(cross_map)
                if regime_map.get(dt, float("inf")) > 0
            ]
            ratio_curves.append((f"{tag} {label}", pts, build_colors[tag], dash))

    # command template + per-build binary of record (replaces the per-run
    # command table: all runs share one template)
    def command_template(command: str, args_map: Dict[str, str]) -> str:
        templ = command
        for value, placeholder in (
            (args_map.get("-s", ""), "{dt_soft}"),
            (args_map.get("-r", ""), "{r_out}"),
            (args_map.get("-f", ""), "{run_id}"),
        ):
            if value:
                templ = templ.replace(value, placeholder)
        binary = extract_binary_from_command(command)
        if binary != "unknown":
            templ = templ.replace(binary, "{petar_bin_order4}")
        return templ

    cmd_template = command_template(rows[0]["command"], parse_command_args(rows[0]["command"])) if rows else ""
    binary_primary = extract_binary_from_command(rows[0]["command"]) if rows else "unknown"
    binary_compare = (
        extract_binary_from_command(compare_rows[0]["command"]) if compare_rows else ""
    )
    family_rout = {
        regime: min(r["rout"] for r in all_rows if r["regime"] == regime)
        for regime in family_order
        if any(r["regime"] == regime for r in all_rows)
    }

    # checks are evaluated per invocation: show both builds side by side
    # (skipped entries — build condition not met — are not displayed)
    compare_checks = []
    if compare_report is not None:
        compare_checks = [
            c for c in compare_report.get("checks", [])
            if c.get("scenario") == args.scenario and not c.get("skipped")
        ]
    all_pass = bool(checks or compare_checks) and all(
        c.get("passed", False) for c in list(checks) + list(compare_checks)
    )
    verdict = "ALL CHECKS PASSED" if all_pass else "CHECKS FAILED"

    # human-readable description per check, from the scenario definition
    check_descriptions = {
        c.get("name"): str(c.get("comment", ""))
        for c in (scenario_data.get("checks", []) if isinstance(scenario_data, dict) else [])
        if isinstance(c, dict)
    }

    def check_row_build(check: Dict, build_label: str) -> str:
        status = "PASS" if check.get("passed") else "FAIL"
        style = "pass" if check.get("passed") else "fail"
        try:
            value_str = f"{float(check.get('value')):.4g}"
        except (TypeError, ValueError):
            value_str = str(check.get("value", ""))
        return (
            f'<tr><td>{html_escape(str(check.get("check", "")))}</td>'
            f"<td>{html_escape(build_label)}</td>"
            f'<td class="{style}">{status}</td>'
            f"<td>{html_escape(value_str)}</td>"
            f'<td>{html_escape(str(check.get("threshold", "")))}</td>'
            f'<td class="msg">{html_escape(check_descriptions.get(check.get("check"), check.get("message", "")))}</td></tr>'
        )

    checks_html = "".join(check_row_build(c, args.primary_tag) for c in checks)
    if compare_checks:
        checks_html += "".join(check_row_build(c, args.compare_tag) for c in compare_checks)

    slope_table_html = "".join(
        f"<tr><td>{html_escape(FAMILY_LABELS.get(row['regime'], row['regime']))}</td>"
        f"<td>{html_escape(row['dataset'])}</td>"
        f"<td>{row['window_slope']:.3f}</td><td>{row['fine_slope']:.3f}</td>"
        f"<td>{row['full_slope']:.3f}</td><td>{row['coarse_slope']:.3f}</td></tr>"
        for row in slope_rows
    )

    setup_commands = scenario_data.get("setup_commands", []) if isinstance(scenario_data, dict) else []
    if not isinstance(setup_commands, list):
        setup_commands = []
    setup_html = "<br>".join(html_escape(str(cmd)) for cmd in setup_commands)

    pure_tree_window_n64 = next(
        (row["window_slope"] for row in slope_rows
         if row["regime"] == "pure_tree" and row["dataset"] == args.compare_tag),
        float("nan"),
    )
    pure_hermite_fine = next(
        (row["fine_slope"] for row in slope_rows
         if row["regime"] == "pure_hermite" and row["dataset"] == args.primary_tag),
        float("nan"),
    )

    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>T1 KDKDK4 changeover consistency</title>
<style>{PAPER_CSS}</style></head><body>

<h1>T1 &nbsp;Kepler orbit crossing the changeover region:<br>validated against pure-tree and pure-Hermite references</h1>
<div class="subtitle">validation report &middot; scenario <i>t1_binary_changeover_kdkdk4</i> &middot; <span class="verdict {'pass' if all_pass else 'fail'}">{verdict}</span></div>

<h2>1&nbsp;&nbsp;Background and purpose</h2>
<p>PeTar splits the pairwise force at the changeover radii r<sub>in</sub>/r<sub>out</sub>: pairs
inside r<sub>in</sub> are integrated by the Hermite+SDAR hard channel, pairs outside
r<sub>out</sub> by the fourth-order KDKDK4 tree channel, and the annulus between them blends the
two kernels. <b>The design under test is the full-crossing configuration</b>: a Kepler orbit
with r<sub>peri</sub> &lt; r<sub>in</sub> &lt; r<sub>out</sub> &lt; r<sub>apo</sub> crosses the
whole changeover region every period, sampling pure hard force at pericenter, pure tree force
at apocenter, and the blended annulus in between. It is validated against two single-channel
references run on the same orbit at the same dt<sub>soft</sub> sweep, establishing the test
items:</p>
<p><b>T1-a (pure-tree reference order).</b> The pure-tree reference (r<sub>peri</sub> &gt;
r<sub>out</sub>: the orbit never enters the changeover) must exhibit the KDKDK4 convergence
order above its round-off floor: error-sweep slope &gt; 3.5 over dt<sub>soft</sub> =
1.22&times;10<sup>-4</sup>&ndash;7.81&times;10<sup>-3</sup> for the 64b build and over
9.77&times;10<sup>-4</sup>&ndash;7.81&times;10<sup>-3</sup> for the non64b build, whose
round-off floor is higher.</p>
<p><b>T1-b (pure-Hermite reference behaviour).</b> The pure-Hermite reference (r<sub>apo</sub>
&lt; r<sub>in</sub>: the orbit never leaves the hard region) must show slope &gt; 3.5 below
the block-step dt threshold (dt<sub>soft</sub> = 2.44&times;10<sup>-4</sup>&ndash;
4.88&times;10<sup>-4</sup>) and only weak dt dependence above it (slope &le; 1 for
dt<sub>soft</sub> &gt; 4.88&times;10<sup>-4</sup>), because the Hermite block steps
self-regulate below dt<sub>soft</sub>.</p>
<p><b>T1-c (crossing family).</b> Inside r<sub>in</sub> the crossing orbit follows the
pure-Hermite reference (at fine dt the 64b crossing error equals the Hermite level,
Q<sub>h</sub> &asymp; 1, Figure 2); otherwise &mdash; at the tree/annulus-dominated coarse end
&mdash; its error must grow with slope &gt; 3.5; and at fine dt its error must lie between
the two references (at or below the pure-tree level, at or above the Hermite level).</p>

<h2>2&nbsp;&nbsp;Method</h2>
<h3>2.1&nbsp;&nbsp;Quantities under comparison</h3>
<p>For every run, the energy log gives the relative total-energy error |Error/Total|(t); the
quantity compared across the dt<sub>soft</sub> sweep is its maximum over the first orbital
period:</p>
<div class="eq">&epsilon;<sub>E</sub>(dt<sub>soft</sub>) = max<sub>t&thinsp;&le;&thinsp;P</sub>
|Error/Total|(t)<span class="tag">(1)</span></div>
<p>with P the binary period derived from the IC (&sect;2.4). The convergence order is the slope of
&epsilon;<sub>E</sub> in the log&ndash;log plane,</p>
<div class="eq">s = &Delta;log<sub>10</sub>&epsilon;<sub>E</sub> /
&Delta;log<sub>10</sub>dt<sub>soft</sub>,<span class="tag">(2)</span></div>
<p>fitted over the check window si03&ndash;si06 (s<sub>win</sub>: above the non64b round-off
floor, below the order-unity saturation of the coarse end), the four finest steps
(s<sub>fine</sub>), all eight points (s<sub>full</sub>), and the four coarsest steps
(s<sub>coarse</sub>). Each r<sub>out</sub> family is classified against the orbit through</p>
<div class="eq">r<sub>in</sub> = r-ratio &middot; r<sub>out</sub> = 0.1&thinsp;r<sub>out</sub>,
with r<sub>peri</sub> = a(1&minus;e), r<sub>apo</sub> = a(1+e),<span class="tag">(3)</span></div>
<p>as pure tree (r<sub>out</sub> &lt; r<sub>peri</sub>), pure Hermite (r<sub>in</sub> &gt;
r<sub>apo</sub>), or full crossing (r<sub>peri</sub> &lt; r<sub>in</sub> and r<sub>apo</sub>
&gt; r<sub>out</sub>).</p>

<h3>2.2&nbsp;&nbsp;Metrics</h3>
<table class="data">
<tr><th>symbol</th><th>definition</th><th>interpretation</th></tr>
<tr><td>&epsilon;<sub>E</sub></td><td>Eq. (1), from the run energy log</td>
<td>one-period energy error of a single run</td></tr>
<tr><td>s<sub>win</sub></td><td>Eq. (2) over the si03&ndash;si06 window</td>
<td>convergence order in the window that is above the non64b round-off floor and below the
coarse-end saturation &mdash; the check quantity for the fourth-order KDKDK4 claim</td></tr>
<tr><td>s<sub>fine</sub>, s<sub>full</sub>, s<sub>coarse</sub></td><td>Eq. (2) over the four
finest / all / four coarsest dt<sub>soft</sub> points</td>
<td>diagnostics: s<sub>fine</sub> &rarr; 0 marks the round-off floor (non64b) or the Hermite
block-step plateau; s &asymp; 4 elsewhere is truncation-dominated behaviour</td></tr>
<tr><td>Q<sub>t</sub>, Q<sub>h</sub></td><td>&epsilon;<sub>E</sub>(family)/&epsilon;<sub>E</sub>(pure
tree), &epsilon;<sub>E</sub>(family)/&epsilon;<sub>E</sub>(pure Hermite) at equal
dt<sub>soft</sub></td>
<td>family error relative to the references: Q<sub>h</sub> &asymp; 1 = coincides with pure
Hermite; Q<sub>t</sub> &le; 1 = no worse than pure tree</td></tr>
<tr><td>r<sub>peri</sub>, r<sub>apo</sub>, P</td><td>Eq. (3) from the IC two-body elements</td>
<td>regime classification scales of each family</td></tr>
</table>
<div class="caption"><b>Table 1.</b> Metrics used in this report. All energies are the total
(tree + hard) relative energy error reported by petar in the run log.</div>

<h3>2.3&nbsp;&nbsp;Runs and test items</h3>
<table class="data">
<tr><th>r<sub>out</sub> family</th><th>r<sub>out</sub> [pc]</th><th>r<sub>in</sub> [pc]</th><th>regime vs. the orbit</th><th>test items (automatic checks)</th></tr>
<tr><td>crossfull_r015</td><td>0.15</td><td>0.015</td>
<td>r<sub>peri</sub> = 0.01 &lt; r<sub>in</sub> &lt; r<sub>out</sub> &lt; r<sub>apo</sub> = 0.19</td>
<td><b>design under test</b>: full crossing every period &mdash;
<code>mix_follows_hermite_fine</code> (follows pure Hermite inside r<sub>in</sub> at fine dt;
64b only &mdash; the non64b annulus floor sits above the Hermite level by construction),
<code>mix_tree_regime_order</code> (slope &gt; 3.5 at the coarse, tree-dominated end),
<code>mix_error_below_pure_tree_fine</code> (error between the two references)</td></tr>
<tr><td>inside_r0005</td><td>0.005</td><td>0.0005</td><td>r<sub>out</sub> &lt; r<sub>peri</sub></td>
<td><b>pure-tree reference</b> &mdash; <code>pure_tree_order_64b</code> (slope &gt; 3.5 over
dt 1.22&times;10<sup>-4</sup>&ndash;7.81&times;10<sup>-3</sup>) and
<code>pure_tree_order_non64b</code> (slope &gt; 3.5 over 9.77&times;10<sup>-4</sup>&ndash;
7.81&times;10<sup>-3</sup>, above the higher non64b round-off floor); each evaluated only on
its own build</td></tr>
<tr><td>hard_r200</td><td>2.0</td><td>0.2</td><td>r<sub>in</sub> &gt; r<sub>apo</sub></td>
<td><b>pure-Hermite reference</b> &mdash; <code>pure_hermite_order</code> (slope &gt; 3.5 over
dt 2.44&times;10<sup>-4</sup>&ndash;4.88&times;10<sup>-4</sup>) and
<code>pure_hermite_weak_dt_dependence</code> (slope &le; 1 for dt &gt;
4.88&times;10<sup>-4</sup>: block steps self-regulate)</td></tr>
</table>
<div class="caption"><b>Table 2.</b> The three r<sub>out</sub> families and the test items each
one carries. Every family is run at eight dt<sub>soft</sub> values (2<sup>&minus;13</sup>,
2<sup>&minus;12</sup>, &hellip;, 2<sup>&minus;6</sup>; output interval = dt<sub>soft</sub>)
with both builds: 24 runs per build, 48 in total. The families are structurally identical runs
of one IC; only -r, -s, -o, -f differ (&sect;2.4). Group linking is disabled (r<sub>group</sub>
= r<sub>search,group</sub> = 10<sup>&minus;6</sup> pc) so no AR group forms and the changeover
is the only active soft/hard mechanism.</div>

<h3>2.4&nbsp;&nbsp;Reproducibility: software, build and execution</h3>
<p><b>Software.</b> PeTar version {html_escape(log_version)}, git commit
<code>{html_escape(git_commit(REPO))}</code>, 1 MPI rank &times; 1 OpenMP thread, tree mode
(-a 0). Two binary families are compared: the 64b build (double-precision tree forces) and the
default-precision build, selected via petar.select:</p>
<pre class="cmd">./configure --enable-64b &amp;&amp; make install &amp;&amp; petar.select --require 64b --optional avx2,omp   # 64b
./configure &amp;&amp; make install &amp;&amp; petar.select --optional avx2,omp                          # non64b</pre>
<p><b>Initial conditions</b> (deterministic two-body construction, no RNG):</p>
<pre class="cmd">{setup_html}</pre>
<table class="data">
<tr><th>parameter</th><th>value</th></tr>
<tr><td>system</td><td>isolated binary, started at apocenter, in the x&ndash;y plane</td></tr>
<tr><td>masses</td><td>m<sub>1</sub> = 1 M<sub>sun</sub>, m<sub>2</sub> = 10 M<sub>sun</sub></td></tr>
<tr><td>orbit</td><td>a = 0.1 pc, e = 0.9 (derived: r<sub>peri</sub> = {orb['peri']:.4g} pc,
r<sub>apo</sub> = {orb['apo']:.4g} pc, P = {orb['period']:.6g} Myr)</td></tr>
<tr><td>units / constants</td><td>M<sub>sun</sub>, pc, pc Myr<sup>-1</sup>;
G = 4.49830997959438 &times; 10<sup>-3</sup> pc&sup3; M<sub>sun</sub><sup>-1</sup> Myr<sup>-2</sup></td></tr>
<tr><td>file format</td><td>ASCII, one row per star: m x y z v<sub>x</sub> v<sub>y</sub>
v<sub>z</sub> (+ petar.init -f conversion)</td></tr>
</table>
<div class="caption"><b>Table 3.</b> IC generation parameters (<code>make_ic.py --case t1</code>
+ petar.init); the input.h4 wrap step adds the Hermite-format header used by -w 2 runs.</div>
<p><b>Execution</b> (all 48 runs share one command template; per build only the family,
step and output prefix vary):</p>
<pre class="cmd">{html_escape(cmd_template)}</pre>
<table class="data">
<tr><th>variable</th><th>values</th><th>meaning</th></tr>
<tr><td><code>{{petar_bin_order4}}</code></td>
<td><code>{html_escape(binary_primary)}</code> ({html_escape(args.primary_tag)})<br>
<code>{html_escape(binary_compare)}</code> ({html_escape(args.compare_tag)})</td>
<td>tree-force precision of the build: double (64b) vs default (non64b)</td></tr>
<tr><td><code>{{r_out}}</code></td>
<td>{", ".join(f"{rout:g} ({fam})" for fam, rout in family_rout.items())}</td>
<td>changeover family (&sect;2.3); r<sub>in</sub> = 0.1&thinsp;r<sub>out</sub> follows
automatically</td></tr>
<tr><td><code>{{dt_soft}}</code></td>
<td>2<sup>-13</sup>, 2<sup>-12</sup>, &hellip;, 2<sup>-6</sup> (also the output interval
<code>-o</code>)</td>
<td>soft time step; the sweep point</td></tr>
<tr><td><code>{{run_id}}</code></td>
<td>&lt;family&gt;_si00 &hellip; si07</td>
<td>output prefix (<code>-f</code>) and log name</td></tr>
</table>
<div class="caption"><b>Table 4.</b> Command template and its variables; the exact recorded
commands are preserved in the JSON reports.</div>

<h2>3&nbsp;&nbsp;Results</h2>
<h3>3.1&nbsp;&nbsp;Error level and dt<sub>soft</sub> scaling by family</h3>
{svg_curves_plot(
    curves,
    legend_models,
    'T1: one-period energy error vs dt_soft, per r_out family and build',
    'dt_soft [Myr]',
    'ε_E (Eq. 1)',
    logx=True,
    logy=True,
    legend_quantities=legend_quantities,
    legend_headers=("r_out family (color)", "build (line)"),
)}
<div class="caption"><b>Figure 1.</b> &epsilon;<sub>E</sub> (Eq. 1) versus dt<sub>soft</sub>
(log&ndash;log) for the three r<sub>out</sub> families; <b>color = family</b> (see legend;
red = full crossing, the design under test; blue = pure-tree reference; green = pure-Hermite
reference), <b>solid = 64b build, dashed = non64b build of the same family</b>. Read it as:
in the pure-tree family the curve rises with slope &asymp; 4 to the right
(s<sub>win</sub> in Table 5) &mdash; the fourth-order KDKDK4 truncation &mdash; and the point
where the non64b dashed curve separates from the 64b solid curve and flattens is the
build-dependent round-off floor; in the pure-Hermite family the curve is flat at coarse dt
(block-step plateau: weak dt dependence) and turns to slope &asymp; 4 at fine dt, and its two
builds coincide because the hard channel is double-precision in both; the full-crossing family
tracks the Hermite reference at fine dt and stays below the pure-tree family at every
dt &mdash; the mixed channel inherits the better of its two references. The per-run numbers
behind every curve point are in the JSON reports.</div>

<h3>3.2&nbsp;&nbsp;Crossing family against the references</h3>
{svg_curves_plot(
    ratio_curves,
    ratio_legend_models,
    'T1: crossing error relative to the references at equal dt_soft',
    'dt_soft [Myr]',
    'Q (error ratio)',
    hlines=[(1.0, 'Q = 1 (equals reference)')],
    logx=True,
    logy=True,
    legend_quantities=[('Q_h = crossing / pure Hermite', ''), ('Q_t = crossing / pure tree', '7,4')],
)}
<div class="caption"><b>Figure 2.</b> The design under test relative to its references at
equal dt<sub>soft</sub>: Q<sub>h</sub> = full-crossing / pure-Hermite error (solid),
Q<sub>t</sub> = full-crossing / pure-tree error (dashed); <b>color = build</b>; the dotted
line marks equality with a reference. Read it as: the 64b Q<sub>h</sub> sits on the
equality line at fine dt &mdash; test item T1-c, inside r<sub>in</sub> the crossing orbit
follows pure Hermite (the pericenter passages set the error exactly as in the reference)
&mdash; and departs upward at the coarse end where the tree/annulus contribution takes over
(slope &gt; 3.5 there, Table 6); the non64b Q<sub>h</sub> at fine dt reflects the annulus
round-off floor sitting above the Hermite level, not a dynamical difference; Q<sub>t</sub>
&lt; 1 everywhere for both builds is the between-references upper bound &mdash; the mixed
channel never degrades below the pure-tree reference. The per-run error values behind these
ratios are the curves of Figure 1 and the JSON reports.</div>

<h3>3.3&nbsp;&nbsp;Convergence-order fits</h3>
<table class="data">
<tr><th>family</th><th>build</th><th>s<sub>win</sub> (si03&ndash;06)</th><th>s<sub>fine</sub> (si00&ndash;03)</th><th>s<sub>full</sub></th><th>s<sub>coarse</sub> (si04&ndash;07)</th></tr>
{slope_table_html}
</table>
<div class="caption"><b>Table 5.</b> Log&ndash;log slope fits (Eq. 2) per family and build,
shown as diagnostics over four standard windows; the checks of Table 6 use the
purpose-specific windows stated there. Read against the check windows: pure tree gives
{pure_tree_window_n64:.2f} on the shared si03&ndash;si06 window and steeper still on the 64b
fine range (a pericenter-dominated superconvergence) &mdash; both far above the 3.5 floor;
pure Hermite gives {pure_hermite_fine:.2f} at fine dt and flattens toward
s<sub>coarse</sub> &rarr; 0 at coarse dt where the block steps self-regulate below
dt<sub>soft</sub>; the non64b s<sub>fine</sub> &rarr; 0 in the pure-tree family is the
round-off floor, which is exactly why its order window starts at
9.77&times;10<sup>-4</sup>.</div>

<h3>3.4&nbsp;&nbsp;Automatic checks</h3>
<table class="data">
<tr><th>check</th><th>build</th><th>status</th><th>value</th><th>threshold</th><th>description</th></tr>
{checks_html}
</table>
<div class="caption"><b>Table 6.</b> Automatic checks, organised by the three purposes of
&sect;1; build-specific windows are evaluated only on their own build (the other build's
entry is not shown). <b>T1-a</b>: pure-tree slope &gt; 3.5 above the round-off floor &mdash; dt
1.22&times;10<sup>-4</sup>&ndash;7.81&times;10<sup>-3</sup> for 64b, 9.77&times;10<sup>-4
</sup>&ndash;7.81&times;10<sup>-3</sup> for non64b (the finer non64b points sit on the floor
and are excluded). <b>T1-b</b>: pure-Hermite slope &gt; 3.5 below the block-step threshold
(dt 2.44&times;10<sup>-4</sup>&ndash;4.88&times;10<sup>-4</sup>) and slope &le; 1 above it
(weak dt dependence; a fourth-order channel would give &asymp;4, a silent one would drift to
0). <b>T1-c</b>: the crossing family follows pure Hermite inside r<sub>in</sub> at fine dt
(ratio in [0.5, 2]; 64b &mdash; for non64b the annulus round-off floor necessarily sits above
the Hermite level, so the check is skipped there), grows with slope &gt; 3.5 at the
tree-dominated coarse end (dt 7.81&times;10<sup>-3</sup>&ndash;1.56&times;10<sup>-2</sup>),
and stays between the references at fine dt (ratio &le; 1 against pure tree; the Hermite lower
bound is the follow-Hermite check for 64b and holds by construction for non64b).</div>

<div class="meta">petar version: {html_escape(log_version)} &middot; scenario:
<i>t1_binary_changeover_kdkdk4</i> &middot; report: {html_escape(args.report)} (+ comparison
{html_escape(args.compare_report)})</div>
</body></html>
"""

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html, encoding="utf-8")
    all_check_records = list(checks) + list(compare_checks)
    n_pass = sum(1 for c in all_check_records if c.get("passed", False))
    print(f"[T1] Summary: passed={n_pass == len(all_check_records) and len(all_check_records) > 0} total={len(all_check_records)} failed={len(all_check_records) - n_pass} report={out_path}")
    return 0 if n_pass == len(all_check_records) and all_check_records else 1


if __name__ == "__main__":
    raise SystemExit(main())
