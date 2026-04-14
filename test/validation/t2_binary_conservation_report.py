#!/usr/bin/env python3
import argparse
import json
import math
import re
import shlex
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

G_MSUN_PC_MYR = 0.00449830997959438
PLOT_LOG_FLOOR = 1e-15
VERSION_RE = re.compile(r"^Version:\s*(.*)$")

Point = Tuple[float, float]


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


def read_ic(ic_path: Path) -> List[List[float]]:
    rows: List[List[float]] = []
    with ic_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.strip().split()
            if not parts:
                continue
            rows.append([float(x) for x in parts])
    return rows


def orbital_elements_from_two_body_rows(rows: List[List[float]]) -> Dict[str, float]:
    p1, p2 = rows[0], rows[1]
    m1, m2 = p1[0], p2[0]
    r = np.array([p2[1] - p1[1], p2[2] - p1[2], p2[3] - p1[3]], dtype=float)
    v = np.array([p2[4] - p1[4], p2[5] - p1[5], p2[6] - p1[6]], dtype=float)

    rnorm = float(np.linalg.norm(r))
    v2 = float(np.dot(v, v))
    mu = G_MSUN_PC_MYR * (m1 + m2)
    eps = 0.5 * v2 - mu / rnorm
    a = -mu / (2.0 * eps)

    h = np.cross(r, v)
    h2 = float(np.dot(h, h))
    e = math.sqrt(max(0.0, 1.0 + 2.0 * eps * h2 / (mu * mu)))
    peri = a * (1.0 - e)
    apo = a * (1.0 + e)
    period = 2.0 * math.pi * math.sqrt(a * a * a / mu) if a > 0 else float("nan")

    return {
        "m1": m1,
        "m2": m2,
        "a": a,
        "e": e,
        "peri": peri,
        "apo": apo,
        "period": period,
    }


def parse_log_metadata(log_path: Path) -> Dict[str, object]:
    text = log_path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    version = "unknown"
    features: List[str] = []
    capture_features = False
    for line in lines:
        stripped = line.strip()
        vmatch = VERSION_RE.match(stripped)
        if vmatch:
            version = vmatch.group(1).strip()
        if stripped == "=====================================" and not capture_features:
            capture_features = True
            continue
        if capture_features:
            if stripped.startswith("Data file output mode"):
                break
            if (
                stripped.startswith("Use ")
                or stripped.startswith("Print ")
                or stripped.startswith("Check ")
                or stripped.startswith("Count ")
                or stripped.startswith("Calculate ")
            ):
                features.append(stripped)
    return {"version": version, "features": features}


def format_tick(value: float) -> str:
    return f"{value:.3g}"


def svg_multi_plot(
    series_map: Dict[str, List[Point]],
    title: str,
    x_label: str,
    y_label: str,
    logx: bool = True,
    logy: bool = True,
    logy_floor_from_nonzero_x: bool = False,
    max_x_ticks: int = 8,
) -> str:
    width, height = 900, 500
    left, right, top, bottom = 88, 24, 44, 70
    plot_w = width - left - right
    plot_h = height - top - bottom

    all_points = [p for points in series_map.values() for p in points]
    if not all_points:
        return "<p>No data available.</p>"

    xs = [p[0] for p in all_points]
    dynamic_floor = PLOT_LOG_FLOOR
    if logy and logy_floor_from_nonzero_x:
        candidates = [p[1] for p in all_points if p[0] > 0.0 and p[1] > 0.0]
        if candidates:
            dynamic_floor = max(min(candidates), PLOT_LOG_FLOOR)

    ys = [max(p[1], dynamic_floor) for p in all_points]

    def tx(v: float) -> float:
        return math.log10(v) if logx else v

    def ty(v: float) -> float:
        return math.log10(v) if logy else v

    xvals = [tx(x) for x in xs]
    yvals = [ty(y) for y in ys]
    xmin, xmax = min(xvals), max(xvals)
    ymin, ymax = min(yvals), max(yvals)
    if xmax == xmin:
        xmax += 1.0
    if ymax == ymin:
        ymax += 1.0

    def sx(v: float) -> float:
        return left + (tx(v) - xmin) / (xmax - xmin) * plot_w

    def sy(v: float) -> float:
        y = max(v, dynamic_floor) if logy else v
        return top + (ymax - ty(y)) / (ymax - ymin) * plot_h

    palette = ["#0b74de", "#d81b60", "#2e7d32", "#8e24aa", "#f4511e", "#00897b"]

    polylines: List[str] = []
    legends: List[str] = []

    for idx, (label, points) in enumerate(series_map.items()):
        color = palette[idx % len(palette)]
        poly = " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y in points)
        polylines.append(f'<polyline points="{poly}" fill="none" stroke="{color}" stroke-width="2"/>')
        ly = top + 18 + idx * 18
        legends.append(
            f'<line x1="{width-right-310}" y1="{ly}" x2="{width-right-290}" y2="{ly}" stroke="{color}" stroke-width="2"/>'
        )
        legends.append(
            f'<text x="{width-right-285}" y="{ly+4}" font-size="11" font-family="Arial">{html_escape(label)}</text>'
        )

    x_unique = sorted(set(xs))
    if len(x_unique) <= max_x_ticks:
        x_tick_values = x_unique
    elif logx:
        px0 = math.floor(min(math.log10(x) for x in x_unique if x > 0.0))
        px1 = math.ceil(max(math.log10(x) for x in x_unique if x > 0.0))
        cand = [10.0 ** p for p in range(int(px0), int(px1) + 1)]
        x_tick_values = [v for v in cand if x_unique[0] <= v <= x_unique[-1]]
        if len(x_tick_values) > max_x_ticks:
            idx_list = [round(i * (len(x_tick_values) - 1) / (max_x_ticks - 1)) for i in range(max_x_ticks)]
            x_tick_values = [x_tick_values[i] for i in idx_list]
    else:
        x0 = min(x_unique)
        x1 = max(x_unique)
        step = (x1 - x0) / (max_x_ticks - 1)
        x_tick_values = [x0 + i * step for i in range(max_x_ticks)]
    if logy:
        p0 = math.floor(min(math.log10(y) for y in ys))
        p1 = math.ceil(max(math.log10(y) for y in ys))
        y_tick_values = [10.0 ** p for p in range(int(p0), int(p1) + 1)]
    else:
        y_tick_values = [min(ys), (min(ys) + max(ys)) * 0.5, max(ys)]

    x_ticks: List[str] = []
    x_grid: List[str] = []
    for value in x_tick_values:
        px = sx(value)
        x_grid.append(f'<line x1="{px:.2f}" y1="{top}" x2="{px:.2f}" y2="{height-bottom}" stroke="#e0e0e0" stroke-dasharray="3,3"/>')
        x_ticks.append(f'<line x1="{px:.2f}" y1="{height-bottom}" x2="{px:.2f}" y2="{height-bottom+6}" stroke="#333"/>')
        x_ticks.append(
            f'<text x="{px:.2f}" y="{height-bottom+22}" text-anchor="middle" font-size="11" font-family="Arial">{format_tick(value)}</text>'
        )

    y_ticks: List[str] = []
    y_grid: List[str] = []
    for value in y_tick_values:
        py = sy(value)
        y_grid.append(f'<line x1="{left}" y1="{py:.2f}" x2="{width-right}" y2="{py:.2f}" stroke="#e0e0e0" stroke-dasharray="3,3"/>')
        y_ticks.append(f'<line x1="{left-6}" y1="{py:.2f}" x2="{left}" y2="{py:.2f}" stroke="#333"/>')
        y_ticks.append(
            f'<text x="{left-10}" y="{py+4:.2f}" text-anchor="end" font-size="11" font-family="Arial">{format_tick(value)}</text>'
        )

    return f"""
<svg width=\"{width}\" height=\"{height}\" viewBox=\"0 0 {width} {height}\" xmlns=\"http://www.w3.org/2000/svg\">
  <rect x=\"0\" y=\"0\" width=\"{width}\" height=\"{height}\" fill=\"white\"/>
  <text x=\"{width/2:.1f}\" y=\"24\" text-anchor=\"middle\" font-size=\"16\" font-family=\"Arial\">{html_escape(title)}</text>
  {' '.join(x_grid)}
  {' '.join(y_grid)}
  <line x1=\"{left}\" y1=\"{height-bottom}\" x2=\"{width-right}\" y2=\"{height-bottom}\" stroke=\"#333\"/>
  <line x1=\"{left}\" y1=\"{top}\" x2=\"{left}\" y2=\"{height-bottom}\" stroke=\"#333\"/>
  {' '.join(polylines)}
  {' '.join(legends)}
  {' '.join(x_ticks)}
  {' '.join(y_ticks)}
  <text x=\"{width/2:.1f}\" y=\"{height-18}\" text-anchor=\"middle\" font-size=\"13\" font-family=\"Arial\">{html_escape(x_label)}</text>
  <text x=\"22\" y=\"{height/2:.1f}\" transform=\"rotate(-90 22,{height/2:.1f})\" text-anchor=\"middle\" font-size=\"13\" font-family=\"Arial\">{html_escape(y_label)}</text>
</svg>
"""


def resolve_status_path(run_record: Dict) -> Path:
    cmd = run_record["command"]
    args_map = parse_command_args(cmd)
    prefix = args_map.get("-f", run_record["run_id"])
    log_path = Path(run_record["output"])
    return log_path.parent / ".." / ".." / "work" / "t2" / f"{prefix}.status"


def extract_binary_from_command(command: str) -> str:
    tokens = shlex.split(command)
    for idx, token in enumerate(tokens):
        if token.startswith("OMP_NUM_THREADS=") and idx + 1 < len(tokens):
            return tokens[idx + 1]
    return "unknown"


def load_status_binary(path: Path):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
    from analysis.status import Status  # pylint: disable=import-outside-toplevel  # pyright: ignore[reportMissingImports]

    st = Status(N_particle=2)
    st.fromfile(str(path))
    return st


def derive_orbital_drift(status_path: Path) -> Dict[str, object]:
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
        "time": [float(x) for x in st.time.tolist()],
        "a0": a0,
        "e0": e0,
        "a_end": float(a[-1]),
        "e_end": float(e[-1]),
        "da_rel": [float(x) for x in da_rel.tolist()],
        "de_abs": [float(x) for x in de_abs.tolist()],
        "da_rel_max": float(np.max(da_rel)),
        "de_abs_max": float(np.max(de_abs)),
    }


def scenario_file_from_name(scenario: str) -> Path:
    return Path("test/validation/scenarios") / f"{scenario}.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate T2 long-term binary conservation HTML report")
    parser.add_argument("--report", default="test/validation/out/report.t2.binary.json")
    parser.add_argument("--scenario", default="t2_binary_conservation_longterm")
    parser.add_argument("--scenario-file", default="")
    parser.add_argument("--ic", default="test/validation/work/t2/input.base")
    parser.add_argument("--output", default="test/validation/out/t2_binary_conservation_report.html")
    args = parser.parse_args()

    report = load_json(Path(args.report))
    scenario_path = Path(args.scenario_file) if args.scenario_file else scenario_file_from_name(args.scenario)
    scenario_data = load_json(scenario_path) if scenario_path.exists() else {}

    runs = [r for r in report.get("runs", []) if r.get("scenario") == args.scenario]
    if not runs:
        raise RuntimeError(f"No runs found for scenario '{args.scenario}' in report: {args.report}")
    checks = [c for c in report.get("checks", []) if c.get("scenario") == args.scenario]

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
        rows.append(
            {
                "run_id": run["run_id"],
                "dt": dt,
                "command": run["command"],
                "log_path": run["output"],
                "status_path": str(status_path),
                "energy_max": float(run["metrics"].get("max_abs_error_over_total", 0.0)),
                "da_rel_max": drift["da_rel_max"],
                "de_abs_max": drift["de_abs_max"],
                "t_end": drift["t_end"],
                "n_sample": drift["n_sample"],
                "a_end": drift["a_end"],
                "e_end": drift["e_end"],
                "time": drift["time"],
                "da_rel": drift["da_rel"],
                "de_abs": drift["de_abs"],
                "args_map": args_map,
            }
        )

    rows = sorted(rows, key=lambda x: x["dt"])

    log_path0 = Path(rows[0]["log_path"])
    metadata = parse_log_metadata(log_path0)

    semi_series = {"max |Δa/a0|": [(r["dt"], r["da_rel_max"]) for r in rows]}
    ecc_series = {"max |Δe|": [(r["dt"], r["de_abs_max"]) for r in rows]}

    semi_time_series: Dict[str, List[Point]] = {}
    ecc_time_series: Dict[str, List[Point]] = {}
    for row in rows:
        label = f"dt_soft={row['dt']:.6e}"
        time_period = [float(t) / orb0["period"] for t in row["time"]]
        semi_time_series[label] = list(zip(time_period, row["da_rel"]))
        ecc_time_series[label] = list(zip(time_period, row["de_abs"]))

    best_semi = min(rows, key=lambda item: item["da_rel_max"])
    worst_semi = max(rows, key=lambda item: item["da_rel_max"])
    best_ecc = min(rows, key=lambda item: item["de_abs_max"])
    worst_ecc = max(rows, key=lambda item: item["de_abs_max"])

    n_periods = [r["t_end"] / orb0["period"] for r in rows]
    min_periods = min(n_periods)

    run_rows_html = "\n".join(
        "<tr>"
        f"<td>{html_escape(r['run_id'])}</td><td>{r['dt']:.6e}</td><td>{r['t_end']:.6e}</td><td>{(r['t_end']/orb0['period']):.2f}</td><td>{r['da_rel_max']:.6e}</td><td>{r['de_abs_max']:.6e}</td><td>{r['n_sample']}</td>"
        "</tr>"
        for r in rows
    )

    command_rows_html = "\n".join(
        "<tr>"
        f"<td>{html_escape(r['run_id'])}</td><td>{r['dt']:.6e}</td><td><code>{html_escape(r['command'])}</code></td><td><code>{html_escape(r['status_path'])}</code></td><td><code>{html_escape(r['log_path'])}</code></td>"
        "</tr>"
        for r in rows
    )

    checks_html = "\n".join(
        f"<li><b>{html_escape(c['check'])}</b>: {'PASS' if c.get('passed') else 'FAIL'} - {html_escape(c.get('message', ''))}</li>"
        for c in checks
    )

    setup_commands = scenario_data.get("setup_commands", []) if isinstance(scenario_data, dict) else []
    if not isinstance(setup_commands, list):
        setup_commands = []
    setup_rows_html = "\n".join(
        f"<tr><td>{idx + 1}</td><td><code>{html_escape(str(cmd))}</code></td></tr>"
        for idx, cmd in enumerate(setup_commands)
    )

    representative_options = rows[0]["args_map"] if rows else {}
    option_rows_html = "\n".join(
        f"<tr><td>{html_escape(key)}</td><td>{html_escape(value)}</td></tr>"
        for key, value in sorted(representative_options.items())
    )

    feature_list = metadata.get("features", [])
    if not isinstance(feature_list, list):
        feature_list = []
    feature_html = "".join(f"<li>{html_escape(str(item))}</li>" for item in feature_list)

    ic_table = "\n".join(
        "<tr>"
        f"<td>{i + 1}</td><td>{row[0]:.6e}</td><td>{row[1]:.6e}</td><td>{row[2]:.6e}</td><td>{row[3]:.6e}</td><td>{row[4]:.6e}</td><td>{row[5]:.6e}</td><td>{row[6]:.6e}</td>"
        "</tr>"
        for i, row in enumerate(ic_rows)
    )

    html = f"""<!doctype html>
<html lang=\"en\"><head><meta charset=\"utf-8\"><title>T2 Binary Long-term Conservation Report</title>
<style>
body{{font-family:Arial,Helvetica,sans-serif;margin:20px;line-height:1.45;}}
h1,h2,h3{{color:#1f2d3d;}}
table{{border-collapse:collapse;width:100%;margin:8px 0 16px 0;}}
th,td{{border:1px solid #cfd8dc;padding:6px 8px;font-size:13px;vertical-align:top;}}
th{{background:#f5f7fa;}}
.note{{background:#eef7ff;border-left:4px solid #0b74de;padding:8px 10px;margin:10px 0;}}
code{{white-space:pre-wrap;word-break:break-all;}}
</style></head><body>
<h1>T2: Long-term Binary Conservation (KDKDK4)</h1>
<div class=\"note\">
Goal: validate long-term orbital conservation for a deterministic pure binary over at least 100 periods, sweeping <code>dt_soft</code> while using PeTar auto-changeover (no <code>-r</code>).<br/>
This scenario targets non64b KDKDK4 by default and keeps <code>r_group</code> small to avoid SDAR-dominated behavior.<br/>
Under auto changeover radii, no explicit 4th-order or monotonic error scaling with <code>dt_soft</code> is assumed; the primary objective is long-term <code>a</code>/<code>e</code> conservation.
</div>

<h2>Figure 1: |Δa/a0| Evolution vs Time (all dt_soft)</h2>
{svg_multi_plot(semi_time_series, 'T2: |Δa/a0|(t) for all dt_soft', 't / P_orb', '|Δa/a0|', logx=False, logy=True, logy_floor_from_nonzero_x=True)}

<h2>Figure 2: |Δe| Evolution vs Time (all dt_soft)</h2>
{svg_multi_plot(ecc_time_series, 'T2: |Δe|(t) for all dt_soft', 't / P_orb', '|Δe|', logx=False, logy=True, logy_floor_from_nonzero_x=True)}

<h2>Figure 3: Semi-major-axis Drift vs Step Size</h2>
{svg_multi_plot(semi_series, 'T2: dt_soft vs max |Δa/a0| over full run', 'dt_soft', 'max |Δa/a0|', logx=True, logy=True)}

<h2>Figure 4: Eccentricity Drift vs Step Size</h2>
{svg_multi_plot(ecc_series, 'T2: dt_soft vs max |Δe| over full run', 'dt_soft', 'max |Δe|', logx=True, logy=True)}

<h2>Initial Conditions</h2>
<table><tr><th>#</th><th>m</th><th>x</th><th>y</th><th>z</th><th>vx</th><th>vy</th><th>vz</th></tr>
{ic_table}
</table>
<ul>
<li>a0 = {orb0['a']:.6e} pc, e0 = {orb0['e']:.6f}, peri = {orb0['peri']:.6e} pc, apo = {orb0['apo']:.6e} pc</li>
<li>One binary period = {orb0['period']:.6e} Myr</li>
<li>Minimum covered periods across dt sweep = {min_periods:.2f}</li>
</ul>

<h2>Run Summary</h2>
<table><tr><th>run_id</th><th>dt_soft</th><th>t_end</th><th>periods covered</th><th>max |Δa/a0|</th><th>max |Δe|</th><th>samples</th></tr>
{run_rows_html}
</table>

<h2>Conservation Extremes</h2>
<table><tr><th>metric</th><th>best run</th><th>best value</th><th>worst run</th><th>worst value</th></tr>
<tr><td>max |Δa/a0|</td><td>{html_escape(best_semi['run_id'])} (dt={best_semi['dt']:.6e})</td><td>{best_semi['da_rel_max']:.6e}</td><td>{html_escape(worst_semi['run_id'])} (dt={worst_semi['dt']:.6e})</td><td>{worst_semi['da_rel_max']:.6e}</td></tr>
<tr><td>max |Δe|</td><td>{html_escape(best_ecc['run_id'])} (dt={best_ecc['dt']:.6e})</td><td>{best_ecc['de_abs_max']:.6e}</td><td>{html_escape(worst_ecc['run_id'])} (dt={worst_ecc['dt']:.6e})</td><td>{worst_ecc['de_abs_max']:.6e}</td></tr>
</table>

<h2>PeTar Version and Executable Info</h2>
<ul>
<li>Executable: <code>{html_escape(extract_binary_from_command(rows[0]['command']))}</code></li>
<li>PeTar version: <code>{html_escape(str(metadata.get('version', 'unknown')))}</code></li>
</ul>

<h3>Build/Runtime Features Parsed from Logs</h3>
<ul>{feature_html}</ul>

<h2>PeTar Commands and Options Used</h2>
<table><tr><th>option</th><th>value</th></tr>
{option_rows_html}
</table>

<h2>Initialization Commands</h2>
<table><tr><th>#</th><th>exact setup command</th></tr>
{setup_rows_html}
</table>

<h2>All Executed Run Commands</h2>
<table><tr><th>run_id</th><th>dt_soft</th><th>exact command</th><th>status path</th><th>log path</th></tr>
{command_rows_html}
</table>

<h2>Analysis Definition</h2>
<ul>
<li><code>max |Δa/a0|</code> and <code>max |Δe|</code> are derived from <code>*.status</code> with <code>N_particle=2</code>, using relative two-body state vectors at each output step.</li>
<li>Two-body element formulas: <code>eps = v^2/2 - mu/r</code>, <code>a = -mu/(2*eps)</code>, <code>e = sqrt(1 + 2*eps*h^2/mu^2)</code>, where <code>mu = G*(m1+m2)</code>.</li>
<li>Drift references are the first output sample (<code>a0</code>, <code>e0</code>).</li>
<li>Interpretation is based on absolute long-term drift magnitudes only; no enforced <code>dt_soft</code> scaling law is assumed in this auto-changeover setup.</li>
</ul>

<h2>Automated Checks</h2>
<ul>{checks_html}</ul>

</body></html>
"""

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html, encoding="utf-8")
    print(f"HTML report written to: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
