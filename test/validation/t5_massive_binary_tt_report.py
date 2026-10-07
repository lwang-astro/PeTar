#!/usr/bin/env python3
"""T5 massive-binary tidal-tensor momentum-conservation HTML report."""
import argparse
import json
import math
import re
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

Point = Tuple[float, float]
VERSION_RE = re.compile(r"^Version:\s*(.*)$")
WORK_DIR = Path("test/validation/work/t5")


def load_json(path: Path) -> Dict:
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def html_escape(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def parse_command_args(command: str) -> Dict[str, str]:
    tokens = shlex.split(command)
    out: Dict[str, str] = {}
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if tok in {"-s", "-t", "-o", "-f", "-b", "-u", "-w", "-a", "-r"} and i + 1 < len(tokens):
            out[tok] = tokens[i + 1]
            i += 2
            continue
        if tok.startswith("--") and i + 1 < len(tokens) and not tokens[i + 1].startswith("-"):
            out[tok] = tokens[i + 1]
            i += 2
            continue
        i += 1
    return out


def format_tick(value: float) -> str:
    return f"{value:.3g}"


def svg_multi_plot(
    series_map: Dict[str, List[Point]],
    title: str,
    x_label: str,
    y_label: str,
    hlines: List[Tuple[float, str]] = [],
    max_x_ticks: int = 8,
) -> str:
    width, height = 940, 500
    left, right, top, bottom = 88, 24, 44, 70
    plot_w = width - left - right
    plot_h = height - top - bottom

    all_points = [p for points in series_map.values() for p in points]
    if not all_points:
        return "<p>No data available.</p>"

    floor = 1e-9
    xs = [p[0] for p in all_points]
    ys = [max(p[1], floor) for p in all_points]

    xmin, xmax = min(xs), max(xs)
    if xmax == xmin:
        xmax += 1.0
    ymin, ymax = math.log10(min(ys)), math.log10(max(ys))
    if ymax == ymin:
        ymax += 1.0

    def sx(v: float) -> float:
        return left + (v - xmin) / (xmax - xmin) * plot_w

    def sy(v: float) -> float:
        return top + (ymax - math.log10(max(v, floor))) / (ymax - ymin) * plot_h

    palette = ["#0b74de", "#d81b60", "#2e7d32", "#8e24aa", "#f4511e", "#00897b"]

    parts: List[str] = []
    legends: List[str] = []
    for idx, (label, points) in enumerate(series_map.items()):
        color = palette[idx % len(palette)]
        poly = " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y in points)
        parts.append(f'<polyline points="{poly}" fill="none" stroke="{color}" stroke-width="2"/>')
        ly = top + 18 + idx * 18
        legends.append(f'<line x1="{width-right-320}" y1="{ly}" x2="{width-right-300}" y2="{ly}" stroke="{color}" stroke-width="2"/>')
        legends.append(f'<text x="{width-right-295}" y="{ly+4}" font-size="11" font-family="Arial">{html_escape(label)}</text>')

    for idx, (value, label) in enumerate(hlines):
        py = sy(value)
        parts.append(f'<line x1="{left}" y1="{py:.2f}" x2="{width-right}" y2="{py:.2f}" stroke="#555" stroke-width="1.2" stroke-dasharray="6,4"/>')
        ly = height - bottom - 10 - idx * 16
        parts.append(f'<text x="{left+6}" y="{ly}" font-size="11" font-family="Arial">{html_escape(label)}</text>')

    x_unique = sorted(set(xs))
    if len(x_unique) <= max_x_ticks:
        x_tick_values = x_unique
    else:
        step = (x_unique[-1] - x_unique[0]) / (max_x_ticks - 1)
        x_tick_values = [x_unique[0] + i * step for i in range(max_x_ticks)]

    p0 = math.floor(ymin)
    p1 = math.ceil(ymax)
    y_tick_values = [10.0 ** p for p in range(int(p0), int(p1) + 1)]

    x_marks: List[str] = []
    for value in x_tick_values:
        px = sx(value)
        x_marks.append(f'<line x1="{px:.2f}" y1="{top}" x2="{px:.2f}" y2="{height-bottom}" stroke="#e0e0e0" stroke-dasharray="3,3"/>')
        x_marks.append(f'<line x1="{px:.2f}" y1="{height-bottom}" x2="{px:.2f}" y2="{height-bottom+6}" stroke="#333"/>')
        x_marks.append(f'<text x="{px:.2f}" y="{height-bottom+22}" text-anchor="middle" font-size="11" font-family="Arial">{format_tick(value)}</text>')

    y_marks: List[str] = []
    for value in y_tick_values:
        py = sy(value)
        y_marks.append(f'<line x1="{left}" y1="{py:.2f}" x2="{width-right}" y2="{py:.2f}" stroke="#e0e0e0" stroke-dasharray="3,3"/>')
        y_marks.append(f'<line x1="{left-6}" y1="{py:.2f}" x2="{left}" y2="{py:.2f}" stroke="#333"/>')
        y_marks.append(f'<text x="{left-10}" y="{py+4:.2f}" text-anchor="end" font-size="11" font-family="Arial">{format_tick(value)}</text>')

    return f"""
<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg">
  <rect x="0" y="0" width="{width}" height="{height}" fill="white"/>
  <text x="{width/2:.1f}" y="24" text-anchor="middle" font-size="16" font-family="Arial">{html_escape(title)}</text>
  {' '.join(x_marks)}
  {' '.join(y_marks)}
  <line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="#333"/>
  <line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="#333"/>
  {' '.join(parts)}
  {' '.join(legends)}
  <text x="{width/2:.1f}" y="{height-18}" text-anchor="middle" font-size="13" font-family="Arial">{html_escape(x_label)}</text>
  <text x="22" y="{height/2:.1f}" transform="rotate(-90 22,{height/2:.1f})" text-anchor="middle" font-size="13" font-family="Arial">{html_escape(y_label)}</text>
</svg>
"""


def read_ic(ic_path: Path) -> List[List[float]]:
    rows: List[List[float]] = []
    with ic_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.strip().split()
            if parts:
                rows.append([float(x) for x in parts])
    return rows


def load_status(status_path: Path, n_particle: int):
    tools_path = (Path(__file__).resolve().parents[2] / "tools").resolve()
    if str(tools_path) not in sys.path:
        sys.path.insert(0, str(tools_path))
    from analysis.status import Status  # pylint: disable=import-outside-toplevel  # pyright: ignore[reportMissingImports]

    st = Status(N_particle=n_particle, interrupt_mode="none")
    st.fromfile(str(status_path))
    return st


def momentum_series(status_path: Path, n_particle: int) -> Dict[str, object]:
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


def parse_log_version(log_path: Path) -> str:
    text = log_path.read_text(encoding="utf-8", errors="replace")
    for line in text.splitlines():
        m = VERSION_RE.match(line.strip())
        if m:
            return m.group(1).strip()
    return "unknown"


def _select_plain_petar() -> str:
    """Select the plain (no interrupt/external) family via petar.select.

    The scenario IC is plain-format, so BSE/galpy-enabled builds cannot read it;
    petar.select without --require picks the plain family. If no plain binary is
    installed, build one (./configure && make install), mirroring the T1-T3
    _select_or_build_petar behavior. Note this changes the configure state.
    """
    import os
    import shutil
    sel = shutil.which("petar.select")
    if not sel:
        raise RuntimeError("petar.select not found in PATH")
    r = subprocess.run([sel, "--optional", "mpi,omp,avx2"], capture_output=True, text=True)
    if r.returncode != 0:
        print("[T5] Plain-family binary not found; building (./configure && make install)...")
        subprocess.run("./configure", shell=True, check=True)
        subprocess.run(f"make -j{os.cpu_count() or 4} install", shell=True, check=True)
        r = subprocess.run([sel, "--optional", "mpi,omp,avx2"], capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(
                "petar.select still cannot select a plain-family binary after rebuild:\n"
                f"{r.stdout}{r.stderr}"
            )
    petar = shutil.which("petar")
    if not petar:
        raise RuntimeError("petar not found in PATH after petar.select")
    return petar


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
    check_by_name = {c["check"]: c for c in checks}
    momentum_thresholds = {
        c["run_id"]: c["threshold"]
        for c in checks
        if c.get("metric") == "max_momentum_drift_pcm"
    }

    ic_rows = read_ic(Path(args.ic))
    n_particle = len(ic_rows)

    mode_labels = {
        "r1_hard_tt": "hard, TT on (isolated group path)",
        "r2_hard_no_tt": "hard, TT off (baseline)",
        "r3_tree_tt": "tree, TT on (event configuration)",
    }

    series: Dict[str, List[Point]] = {}
    table_rows: List[str] = []
    command_rows: List[str] = []
    version = "unknown"
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
        label = mode_labels.get(run_id, run_id)
        series[label] = list(zip(ms["time"], ms["drift"]))
        threshold = momentum_thresholds.get(run_id)
        verdict = ""
        if threshold is not None:
            ok = ms["max_drift"] <= threshold
            verdict = (
                '<span style="color:#2e7d32">PASS</span>' if ok
                else '<span style="color:#d81b60">FAIL</span>'
            )
        table_rows.append(
            "<tr>"
            f"<td>{html_escape(run_id)}</td><td>{html_escape(label)}</td><td>{ms['n_sample']}</td>"
            f"<td>{ms['max_drift']:.4e}</td><td>{threshold if threshold is not None else '—'}</td><td>{verdict}</td>"
            f"<td>{artificials:.0f}</td>"
            "</tr>"
        )
        command_rows.append(
            f"<tr><td>{html_escape(run_id)}</td><td><code>{html_escape(run['command'])}</code></td></tr>"
        )

    ic_table = "\n".join(
        "<tr>"
        f"<td>{i + 1}</td><td>{row[0]:.6e}</td><td>{row[1]:.6e}</td><td>{row[2]:.6e}</td><td>{row[3]:.6e}</td><td>{row[4]:.6e}</td><td>{row[5]:.6e}</td><td>{row[6]:.6e}</td>"
        "</tr>"
        for i, row in enumerate(ic_rows)
    )

    check_descriptions = {
        "r1_hard_tt_momentum_conserved": "max ||P|-|P0||/M in pure hard mode with TT on. With 3 particles the churner orbit overlaps the group-linking radius, so the system stays one isolated AR group; drift must stay at integration-noise level.",
        "r2_hard_no_tt_momentum_conserved": "Same as r1 with TT off: the baseline drift level without tidal tensor.",
        "r3_tree_tt_momentum_conserved": "max ||P|-|P0||/M in tree mode with TT engaged and group-boundary churn (the production configuration of the 2026-10-06 event). Threshold 0.005 calibrated so the pre-fix binary (~1.0e-2) fails and the fixed binary (~1.4e-3) passes with margin.",
        "r1_hard_tt_matches_no_tt_baseline": "Artificial-particle count in r1 must be zero: TT is structurally inactive in the isolated-group path, so r1 must reproduce the r2 baseline.",
        "r3_tt_engaged": "Artificial-particle count in r3 must be >= 1, proving the tidal tensor was actually exercised.",
        "r2_no_tt_particles": "Artificial-particle count in r2 must be zero (TT off).",
        "r1_no_runtime_abort": "No segmentation faults, assertion failures, or aborts in r1.",
        "r2_no_runtime_abort": "No segmentation faults, assertion failures, or aborts in r2.",
        "r3_no_runtime_abort": "No segmentation faults, assertion failures, or aborts in r3.",
    }
    checks_html = "\n".join(
        "<li><b>{name}</b>: {status} — {msg}<br><i>{desc}</i></li>".format(
            name=html_escape(c["check"]),
            status='<span style="color:#2e7d32">PASS</span>' if c.get("passed") else '<span style="color:#d81b60">FAIL</span>',
            msg=html_escape(c.get("message", "")),
            desc=html_escape(check_descriptions.get(c["check"], "")),
        )
        for c in checks
    )

    summary = report.get("summary", {})
    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>T5 Massive-Binary Tidal-Tensor Momentum Report</title>
<style>
body{{font-family:Arial,Helvetica,sans-serif;margin:20px;line-height:1.45;}}
h1,h2,h3{{color:#1f2d3d;}}
table{{border-collapse:collapse;width:100%;margin:8px 0 16px 0;}}
th,td{{border:1px solid #cfd8dc;padding:6px 8px;font-size:13px;vertical-align:top;}}
th{{background:#f5f7fa;}}
.note{{background:#eef7ff;border-left:4px solid #0b74de;padding:8px 10px;margin:10px 0;}}
code{{white-space:pre-wrap;word-break:break-all;}}
</style></head><body>
<h1>T5: Massive Binary + Light Churner, Tidal-Tensor Momentum Conservation</h1>
<div class="note">
Regression test for the 2026-10-06 tidal-tensor momentum leak: the even-order (T3 ~ r^-4)
term of the perturber tidal tensor had a non-zero mass-weighted mean over AR group members,
pushing the group center of mass every AR substep with no counter-reaction. The fix removes
the tensor member mean in <code>calcAccPert</code> (src/ar_interaction.hpp). The IC reproduces
the event: a 77.65 Msun binary (a = 0.0081 pc, e = 0.667) plus a 0.242 Msun churner whose
pericenter (0.01 pc) lies inside the group radius and apocenter (0.04 pc) outside it, so AR
group membership churns every outer orbit. Threshold calibration: pre-fix binary drifts to
~1.0e-2 (fails the r3 threshold), fixed binary stays at ~1.4e-3.
</div>

<h2>Figure 1: Momentum Drift ||P|-|P0||/M vs Time</h2>
{svg_multi_plot(series, 'T5: momentum drift for TT-on/off and hard/tree modes', 't [Myr]', '||P|-|P0||/M [pc/Myr]', hlines=[(0.005, 'r3 threshold 5e-3 (pre-fix ~1.0e-2 fails)'), (0.02, 'r1/r2 threshold 2e-2')], max_x_ticks=8)}

<h2>Initial Conditions (make_ic.py --case t5)</h2>
<table><tr><th>#</th><th>m [Msun]</th><th>x</th><th>y</th><th>z</th><th>vx</th><th>vy</th><th>vz</th></tr>
{ic_table}
</table>

<h2>Run Summary</h2>
<table><tr><th>run_id</th><th>mode</th><th>status samples</th><th>max drift [pc/Myr]</th><th>threshold</th><th>verdict</th><th>max artificial particles</th></tr>
{"".join(table_rows)}
</table>

<h2>Automatic Checks</h2>
<ul>{checks_html}</ul>

<h2>PeTar Version</h2>
<ul><li>PeTar version: <code>{html_escape(version)}</code></li></ul>

<h2>Commands</h2>
<table><tr><th>run_id</th><th>command</th></tr>
{"".join(command_rows)}
</table>

<hr>
<p>Scenario summary: passed={summary.get("passed")}, total={summary.get("n_total", summary.get("total"))}, failed={summary.get("n_failed", summary.get("failed"))}. Report JSON: <code>{html_escape(args.report)}</code></p>
</body></html>
"""

    Path(args.output).write_text(html, encoding="utf-8")
    print(f"[T5] HTML report written: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
