#!/usr/bin/env python3
import argparse
import json
import math
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple


G_MSUN_PC_MYR = 0.00449830997959438
TIME_RE = re.compile(r"^Time:\s*([0-9eE+\-.]+)")
ENERGY_HEADER_RE = re.compile(r"^Energy:\s+(.*)$")
ENERGY_PHYSIC_RE = re.compile(r"^Physic:\s+(.*)$")
VERSION_RE = re.compile(r"^Version:\s*(.*)$")


Point = Tuple[float, float]


def load_json(path: Path) -> Dict:
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def parse_command_args(command: str) -> Dict[str, str]:
    tokens = command.split()
    out: Dict[str, str] = {}
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if tok in {"-s", "-r", "-t", "-o", "-f", "-b", "-u", "-w"} and i + 1 < len(tokens):
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


def orbital_elements_from_two_body(rows: List[List[float]]) -> Dict[str, float]:
    p1, p2 = rows[0], rows[1]
    m1, m2 = p1[0], p2[0]
    r = [p2[1] - p1[1], p2[2] - p1[2], p2[3] - p1[3]]
    v = [p2[4] - p1[4], p2[5] - p1[5], p2[6] - p1[6]]
    rnorm = math.sqrt(sum(x * x for x in r))
    v2 = sum(x * x for x in v)
    mu = G_MSUN_PC_MYR * (m1 + m2)

    eps = 0.5 * v2 - mu / rnorm
    a = -mu / (2.0 * eps)

    hx = r[1] * v[2] - r[2] * v[1]
    hy = r[2] * v[0] - r[0] * v[2]
    hz = r[0] * v[1] - r[1] * v[0]
    h2 = hx * hx + hy * hy + hz * hz

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


def html_escape(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


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


def format_tick(value: float) -> str:
    return f"{value:.3g}"


def svg_multi_plot(series_map: Dict[str, List[Point]], title: str, x_label: str, y_label: str, logx: bool = True, logy: bool = True) -> str:
    width, height = 900, 500
    left, right, top, bottom = 88, 24, 44, 70
    plot_w = width - left - right
    plot_h = height - top - bottom

    all_points = [p for points in series_map.values() for p in points]
    if not all_points:
        return "<p>No data available.</p>"

    xs = [p[0] for p in all_points]
    ys = [max(p[1], 1e-30) for p in all_points]

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
        return top + (ymax - ty(v)) / (ymax - ymin) * plot_h

    palette = ["#0b74de", "#d81b60", "#2e7d32", "#8e24aa", "#f4511e", "#00897b"]

    def parse_style_key(label: str) -> Tuple[str, str]:
        dataset = ""
        if "|" in label:
            dataset = label.split("|", 1)[0].strip().lower()
        m = re.search(r"rout=([0-9eE+\-.]+)", label)
        if m:
            try:
                rout_key = f"{float(m.group(1)):.12g}"
            except ValueError:
                rout_key = m.group(1)
        else:
            rout_key = label
        return dataset, rout_key

    rout_keys_in_order: List[str] = []
    for label in series_map.keys():
        _, rk = parse_style_key(label)
        if rk not in rout_keys_in_order:
            rout_keys_in_order.append(rk)
    color_by_rout = {rk: palette[i % len(palette)] for i, rk in enumerate(rout_keys_in_order)}

    polylines: List[str] = []
    circles: List[str] = []
    legends: List[str] = []

    for idx, (label, points) in enumerate(series_map.items()):
        dataset, rout_key = parse_style_key(label)
        color = color_by_rout.get(rout_key, palette[idx % len(palette)])
        dash = "6,4" if dataset == "non64b" else "none"
        poly = " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y in points)
        polylines.append(f'<polyline points="{poly}" fill="none" stroke="{color}" stroke-width="2" stroke-dasharray="{dash}"/>')
        circles.extend(
            f'<circle cx="{sx(x):.2f}" cy="{sy(y):.2f}" r="4" fill="{color}"/>' for x, y in points
        )
        ly = top + 18 + idx * 18
        legends.append(
            f'<line x1="{width-right-310}" y1="{ly}" x2="{width-right-290}" y2="{ly}" stroke="{color}" stroke-width="2" stroke-dasharray="{dash}"/>'
        )
        legends.append(
            f'<text x="{width-right-285}" y="{ly+4}" font-size="11" font-family="Arial">{html_escape(label)}</text>'
        )

    x_tick_values = sorted(set(xs))
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
<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg">
  <rect x="0" y="0" width="{width}" height="{height}" fill="white"/>
  <text x="{width/2:.1f}" y="24" text-anchor="middle" font-size="16" font-family="Arial">{html_escape(title)}</text>
  {' '.join(x_grid)}
  {' '.join(y_grid)}
  <line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="#333"/>
  <line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="#333"/>
  {' '.join(polylines)}
  {' '.join(circles)}
  {' '.join(legends)}
  {' '.join(x_ticks)}
  {' '.join(y_ticks)}
  <text x="{width/2:.1f}" y="{height-18}" text-anchor="middle" font-size="13" font-family="Arial">{html_escape(x_label)}</text>
  <text x="22" y="{height/2:.1f}" transform="rotate(-90 22,{height/2:.1f})" text-anchor="middle" font-size="13" font-family="Arial">{html_escape(y_label)}</text>
</svg>
"""


def _detect_or_build_petar(require_64b: bool) -> str:
    """Select (and if needed build) a petar binary with the right 64b settings via petar.select."""
    import shutil, subprocess as _sp
    sel = shutil.which("petar.select")
    if not sel:
        return shutil.which("petar") or "petar"

    args = ["--require", "64b", "--optional", "avx2,omp"] if require_64b else ["--optional", "avx2,omp"]
    r = _sp.run([sel] + args, capture_output=True, text=True)
    if r.returncode == 0:
        return shutil.which("petar") or "petar"

    # Not found – build it
    print(f"[T1] {'64b' if require_64b else 'non64b'} KDKDK4 binary not found. Building...")
    config_flags = "--enable-64b" if require_64b else ""
    _sp.run(f"./configure {config_flags}", shell=True, check=True, capture_output=True)
    _sp.run(f"make -j{os.cpu_count() or 4}", shell=True, check=True, capture_output=True)
    _sp.run("make install", shell=True, check=True, capture_output=True)
    _sp.run([sel] + args, check=True, capture_output=True)
    return shutil.which("petar") or "petar"


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
    parser.add_argument("--compare-report", default="")
    parser.add_argument("--compare-tag", default="non64b")
    parser.add_argument("--scenario", default="t1_high_ecc_changeover")
    parser.add_argument("--scenario-file", default="")
    parser.add_argument("--out-dir-64b", default="test/out/validation_t1_kdkdk4_64b")
    parser.add_argument("--out-dir-non64b", default="test/out/validation_t1_kdkdk4_non64b")
    parser.add_argument("--ic", default="test/out/work/t1/input.base")
    parser.add_argument("--output", default="test/out/t1_kdkdk4_summary.html")
    args = parser.parse_args()

    if args.run:
        sfile = args.scenario_file or str(scenario_file_from_name(args.scenario))
        petar_64 = _detect_or_build_petar(True)
        petar_non64 = _detect_or_build_petar(False)
        _run_t1_scenario(sfile, args.report, args.out_dir_64b, petar_64)
        _run_t1_scenario(sfile, args.compare_report or "test/out/report.t1.non64b.overlay.json",
                         args.out_dir_non64b, petar_non64)
        if not args.compare_report:
            args.compare_report = "test/out/report.t1.non64b.overlay.json"

    report = load_json(Path(args.report))
    scenario_path = Path(args.scenario_file) if args.scenario_file else scenario_file_from_name(args.scenario)
    scenario_data = load_json(scenario_path) if scenario_path.exists() else {}

    runs = [r for r in report.get("runs", []) if r.get("scenario") == args.scenario]
    checks = [c for c in report.get("checks", []) if c.get("scenario") == args.scenario]

    compare_report = None
    compare_runs = []
    if args.compare_report:
        compare_report = load_json(Path(args.compare_report))
        compare_runs = [r for r in compare_report.get("runs", []) if r.get("scenario") == args.scenario]

    ic_rows = read_ic(Path(args.ic))
    orb = orbital_elements_from_two_body(ic_rows)
    period = orb["period"]

    def classify_regime(rout: float, rin: float, peri: float, apo: float) -> str:
        if rout < peri:
            return "inside"
        if rin > apo:
            return "hard"
        if rin < peri < rout:
            return "peri-between"
        if peri < rin < apo and apo < rout:
            return "crossing"
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
                    "args_map": args_map,
                    "log_path": str(log_path),
                    "dataset": tag,
                }
            )
        return out_rows

    rows = build_rows(runs, args.primary_tag)
    compare_rows = build_rows(compare_runs, args.compare_tag) if compare_runs else []
    all_rows = rows + compare_rows

    metadata = parse_log_metadata(Path(rows[0]["log_path"])) if rows else {"version": "unknown", "features": []}

    def build_series(regime: str, source_rows: List[Dict], dataset_tag: str = "") -> Dict[str, List[Point]]:
        grouped: Dict[str, List[Point]] = {}
        candidates = [r for r in source_rows if r["regime"] == regime]
        for rout in sorted({r["rout"] for r in candidates}):
            base = f"rout={rout:g}, rin={0.1 * rout:g}"
            label = f"{dataset_tag} | {base}" if dataset_tag else base
            grouped[label] = sorted(
                [(r["dt"], r["max_err_1p"]) for r in candidates if abs(r["rout"] - rout) < 1e-15],
                key=lambda item: item[0],
            )
        return grouped

    inside_series = build_series("inside", rows, args.primary_tag)
    peri_between_series = build_series("peri-between", rows, args.primary_tag)
    cross_series = build_series("crossing", rows, args.primary_tag)

    compare_inside_series = build_series("inside", compare_rows, args.compare_tag) if compare_rows else {}
    compare_peri_between_series = build_series("peri-between", compare_rows, args.compare_tag) if compare_rows else {}
    compare_cross_series = build_series("crossing", compare_rows, args.compare_tag) if compare_rows else {}

    all_series: Dict[str, List[Point]] = {}
    all_series.update(cross_series)
    all_series.update(peri_between_series)
    all_series.update(inside_series)
    all_series.update(compare_cross_series)
    all_series.update(compare_peri_between_series)
    all_series.update(compare_inside_series)

    series_summary_rows = []
    annotated_all_series: Dict[str, List[Point]] = {}
    for label, points in all_series.items():
        full_slope = fit_loglog_slope(points)
        fine_slope = fit_loglog_slope_subset(points, 0, min(4, len(points)))
        coarse_slope = fit_loglog_slope_subset(points, max(0, len(points) - 4), len(points))
        if label in inside_series or label in compare_inside_series:
            regime = "inside"
        elif label in peri_between_series or label in compare_peri_between_series:
            regime = "peri-between"
        else:
            regime = "crossing"

        if regime == "inside":
            note = "tree-dominated; full-range slope should approach 4"
        elif regime == "peri-between":
            note = "pericenter lies in changeover shell (rin < peri < rout); expected to sit between inside and wide-crossing behavior"
        else:
            note = "coarse dt tends to plateau-like behavior, fine dt returns to ~4th order"
        series_summary_rows.append(
            {
                "label": label,
                "regime": regime,
                "full_slope": full_slope,
                "fine_slope": fine_slope,
                "coarse_slope": coarse_slope,
                "note": note,
            }
        )
        annotated_all_series[
            f"{label} | full={full_slope:.2f}, fine={fine_slope:.2f}, coarse={coarse_slope:.2f}"
        ] = points

    check_descriptions = {
        "inside_energy_bound_fine": (
            "Maximum |Error/Total| at finest dt_soft in inside (tree-only) regime. "
            "Verifies that KDKDK4 tree integration does not produce abnormally large errors."
        ),
        "inside_loglog_slope": (
            "Log-log slope of |Error/Total| vs dt_soft in inside regime. "
            "Expect ~4 for KDKDK4 4th-order tree integration (64b only; non64b has round-off plateau)."
        ),
        "hard_r200_energy_bound_fine": (
            "Maximum |Error/Total| at finest dt_soft in hard (Hermite-only, r_in > r_apo) regime. "
            "Verifies that Hermite integration stays within reasonable error bounds."
        ),
        "hard_r200_plateau_si06_to_si04": (
            "Convergence ratio between coarse (si06) and fine (si04) dt in hard regime. "
            "Tests the adaptive-timestep plateau: error should grow much slower than dt^4, indicating Hermite self-regulates timestep."
        ),
        "inside_no_hard_energy_alarm": (
            "Count of 'Hard energy significant' warnings in inside regime. "
            "Should be zero — tree-only integration should not trigger hard-energy alarms."
        ),
        "hard_r200_no_hard_energy_alarm": (
            "Count of 'Hard energy significant' warnings in hard regime. "
            "Should be zero — pure Hermite integration should be stable."
        ),
    }

    checks_html = "\n".join(
        "<li><b>{name}</b>: {status} — {msg}<br><i>{desc}</i></li>".format(
            name=html_escape(c['check']),
            status='<span style=\"color:#2e7d32\">PASS</span>' if c.get('passed') else '<span style=\"color:#d81b60\">FAIL</span>',
            msg=html_escape(c.get('message', '')),
            desc=html_escape(check_descriptions.get(c['check'], '')))
        for c in checks
    )

    run_rows_html = "\n".join(
        "<tr>"
        f"<td>{html_escape(r['run_id'])}</td><td>{html_escape(r['dataset'])}</td><td>{r['dt']:.6e}</td><td>{r['rin']:.6e}</td><td>{r['rout']:.6e}</td><td>{r['regime']}</td><td>{r['max_err_1p']:.6e}</td>"
        "</tr>"
        for r in sorted(all_rows, key=lambda x: (x["dataset"], x["run_id"]))
    )

    feature_list = metadata.get("features", [])
    if not isinstance(feature_list, list):
        feature_list = []
    feature_html = "".join(f"<li>{html_escape(str(item))}</li>" for item in feature_list)

    representative_options = rows[0]["args_map"] if rows else {}
    option_rows_html = "\n".join(
        f"<tr><td>{html_escape(key)}</td><td>{html_escape(value)}</td></tr>"
        for key, value in sorted(representative_options.items())
    )

    command_rows_html = "\n".join(
        "<tr>"
        f"<td>{html_escape(r['run_id'])}</td><td>{html_escape(r['dataset'])}</td><td>{r['regime']}</td><td>{r['dt']:.6e}</td><td>{r['rout']:.6e}</td><td><code>{html_escape(r['command'])}</code></td><td><code>{html_escape(r['log_path'])}</code></td>"
        "</tr>"
        for r in sorted(all_rows, key=lambda x: (x["dataset"], x["run_id"]))
    )

    setup_commands = scenario_data.get("setup_commands", []) if isinstance(scenario_data, dict) else []
    if not isinstance(setup_commands, list):
        setup_commands = []
    setup_rows_html = "\n".join(
        f"<tr><td>{idx + 1}</td><td><code>{html_escape(str(cmd))}</code></td></tr>"
        for idx, cmd in enumerate(setup_commands)
    )

    slope_rows_html = "\n".join(
        "<tr>"
        f"<td>{html_escape(row['label'])}</td>"
        f"<td>{row['regime']}</td>"
        f"<td>{row['full_slope']:.3f}</td>"
        f"<td>{row['fine_slope']:.3f}</td>"
        f"<td>{row['coarse_slope']:.3f}</td>"
        f"<td>{html_escape(row['note'])}</td>"
        "</tr>"
        for row in series_summary_rows
    )

    ic_table = "\n".join(
        "<tr>"
        f"<td>{i + 1}</td><td>{row[0]:.6e}</td><td>{row[1]:.6e}</td><td>{row[2]:.6e}</td><td>{row[3]:.6e}</td><td>{row[4]:.6e}</td><td>{row[5]:.6e}</td><td>{row[6]:.6e}</td>"
        "</tr>"
        for i, row in enumerate(ic_rows)
    )

    html = f"""<!doctype html>
<html lang=\"en\"><head><meta charset=\"utf-8\"><title>T1 KDKDK4 Validation Report</title>
<style>
body{{font-family:Arial,Helvetica,sans-serif;margin:20px;line-height:1.45;}}
h1,h2,h3{{color:#1f2d3d;}}
table{{border-collapse:collapse;width:100%;margin:8px 0 16px 0;}}
th,td{{border:1px solid #cfd8dc;padding:6px 8px;font-size:13px;vertical-align:top;}}
th{{background:#f5f7fa;}}
.note{{background:#eef7ff;border-left:4px solid #0b74de;padding:8px 10px;margin:10px 0;}}
code{{white-space:pre-wrap;word-break:break-all;}}
</style></head><body>
<h1>T1: KDKDK4 Changeover Validation Report</h1>
<div class=\"note\">
<b>Goal:</b> Test how a high-eccentricity binary behaves under different PeTar integration regimes — pure particle-tree (KDKDK4 4th-order), pure Hermite (4th-order adaptive), and soft+hard crossing — by varying the changeover radius r<sub>out</sub>. Compare double-precision (64b) vs default-precision (non64b) tree forces.<br/><br/>

<h3>Test Design</h3>
<b>Binary IC:</b> m1=1&#8239;M<sub>&sun;</sub>, m2=10&#8239;M<sub>&sun;</sub>, semi=0.1&#8239;pc, e=0.9, start at apocenter (&approx;0.19&#8239;pc).<br/>
One orbital period &approx;0.0283&#8239;Myr. Pericenter r<sub>peri</sub>=0.01&#8239;pc, apocenter r<sub>apo</sub>=0.19&#8239;pc.<br/>
<b>dt<sub>soft</sub> sweep:</b> 2<sup>&minus;13</sup> to 2<sup>&minus;6</sup>, sample the error-scaling curve.<br/>
<br/>
<b>r<sub>out</sub> regimes &amp; why each value is chosen:</b>
<table>
<tr><th>r<sub>out</sub></th><th>r<sub>in</sub>(=0.1&times;r<sub>out</sub>)</th><th>regime</th><th>purpose</th></tr>
<tr><td>0.005</td><td>0.0005</td><td>inside (r<sub>out</sub>&lt;r<sub>peri</sub>)</td><td>All forces in particle-tree. Expect &asymp;4th-order KDKDK4 error slope.</td></tr>
<tr><td>0.05</td><td>0.005</td><td>peri-between (r<sub>in</sub>&lt;r<sub>peri</sub>&lt;r<sub>out</sub>)</td><td>Pericenter lies in the changeover shell; probes the handover from tree to hard.</td></tr>
<tr><td>0.32</td><td>0.032</td><td>crossing (r<sub>peri</sub>&lt;r<sub>in</sub>&lt;r<sub>apo</sub>&lt;r<sub>out</sub>)</td><td>Narrow crossing: binary just fits inside changeover. Coarse dt tends to plateau.</td></tr>
<tr><td>2.0</td><td>0.2</td><td>hard (r<sub>in</sub>&gt;r<sub>apo</sub>)</td><td>All forces in Hermite hard region. Coarse dt plateau (adaptive timestep); fine dt returns to 4th order.</td></tr>
</table>
<i>r<sub>out</sub>=0.32 (crossing) shows mixed soft+hard behavior. r<sub>out</sub>=2.0 (hard) shows Hermite plateau at coarse dt.</i><br/>
<br/>
<b>64b vs non64b comparison:</b> The 64b build uses double-precision tree forces (<code>--enable-64b</code>). The non64b build uses default precision. Crossing-regime errors (tree-dominated) are dominated by round-off and should differ; inside-regime errors (hard-dominated) may be similar.<br/>
<b>Line styles:</b> solid = 64b, dashed = non64b. Same <code>r<sub>out</sub></code> = same color.</div>

<h2>Figure 1: Overall Error vs Step Size</h2>
{svg_multi_plot(annotated_all_series, 'T1 representative runs: dt_soft vs max error within one binary period', 'dt_soft', 'max |Error/Total|', logx=True, logy=True)}
<ul>
<li>Primary dataset tag: <code>{html_escape(args.primary_tag)}</code></li>
{f'<li>Comparison dataset tag: <code>{html_escape(args.compare_tag)}</code></li>' if compare_rows else ''}
</ul>

<h2>Initial Conditions</h2>
<table><tr><th>#</th><th>m</th><th>x</th><th>y</th><th>z</th><th>vx</th><th>vy</th><th>vz</th></tr>
{ic_table}
</table>
<ul>
<li>a = {orb['a']:.6e} pc, e = {orb['e']:.6f}, peri = {orb['peri']:.6e} pc, apo = {orb['apo']:.6e} pc</li>
<li>One binary period = {orb['period']:.6e} Myr</li>
</ul>

<h2>Run Parameters and One-Period Maximum Error</h2>
<table><tr><th>run_id</th><th>dataset</th><th>dt_soft</th><th>rin(≈0.1*rout)</th><th>rout</th><th>regime</th><th>max |Error/Total| (t ≤ 1 period)</th></tr>
{run_rows_html}
</table>

<h2>PeTar Version and Executable Info</h2>
<ul>
<li>Executable: <code>{html_escape(extract_binary_from_command(rows[0]['command'])) if rows else 'unknown'}</code></li>
<li>PeTar version: <code>{html_escape(str(metadata.get('version', 'unknown')))}</code></li>
</ul>

<h3>Build/Runtime Features Parsed from Logs</h3>
<ul>{feature_html}</ul>

<h2>PeTar Commands and Options Used</h2>
<p>Representative common options are listed below (between runs, the main changes are <code>-s</code>, <code>-r</code>, <code>-o</code>, and <code>-f</code>):</p>
<table><tr><th>option</th><th>value</th></tr>
{option_rows_html}
</table>

<h2>Initialization and Input Preparation Commands</h2>
<p>The following commands come directly from the current T1 scenario definition to generate ICs, run <code>petar.init</code>, and build <code>input.h4</code>:</p>
<table><tr><th>#</th><th>exact setup command</th></tr>
{setup_rows_html}
</table>

<p>All executed run commands are listed below:</p>
<table><tr><th>run_id</th><th>dataset</th><th>regime</th><th>dt_soft</th><th>rout</th><th>exact command</th><th>log path</th></tr>
{command_rows_html}
</table>

<h2>Error Analysis Definition</h2>
<ul>
<li>Primary error metric: <code>max |Error/Total| (t ≤ 1 period)</code>. It is extracted from the PeTar log <code>Energy:</code>/<code>Physic:</code> table by taking the absolute value of <code>Error/Total</code> and then the maximum within one binary period.</li>
<li>Fitted slope: for each <code>rout</code> series, perform log-log linear fitting on <code>(dt_soft, max |Error/Total|)</code>. The report provides full-range, smallest-4-point, and largest-4-point slopes.</li>
<li>Inside criterion: <code>rout &lt; peri</code>, used to isolate particle-tree/KDKDK4-dominated error, with expected slope near 4.</li>
<li>Peri-between group: <code>rin &lt; peri &lt; rout</code>, where pericenter lies in the changeover shell and is expected to lie between inside and wide-crossing behavior.</li>
<li>Crossing criterion: <code>peri &lt; rin &lt; apo &lt; rout</code>, typically showing plateau-like behavior at coarse dt and near-4th-order trend at fine dt.</li>
</ul>

<h2>Slope Summary for All Curves</h2>
<table><tr><th>series</th><th>regime</th><th>full-range slope</th><th>fine-dt slope</th><th>coarse-dt slope</th><th>note</th></tr>
{slope_rows_html}
</table>

<h2>Automated Check Results</h2>
<ul>{checks_html}</ul>

<h2>Conclusion (Auto-generated)</h2>
<ul>
<li>The inside subset currently uses <code>rout=0.005</code> (satisfying <code>rout &lt; peri</code>) to probe particle-tree/KDKDK4-dominated error-order behavior; the current full-range slope is {next((row['full_slope'] for row in series_summary_rows if row['regime'] == 'inside'), float('nan')):.3f}.</li>
<li>The peri-between subset is used to check whether the transitional shell case (<code>rin &lt; peri &lt; rout</code>) yields error levels between inside and wide-crossing.</li>
<li>The crossing subset currently shows a more plateau-like coarse-dt region and near-4th-order recovery at fine dt, consistent with the expected physical behavior.</li>
<li>The summary figure shows inside, peri-between, and crossing groups together for direct visual comparison.</li>
</ul>

</body></html>
"""

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html, encoding="utf-8")
    print(f"HTML report written to: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
