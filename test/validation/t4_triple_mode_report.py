#!/usr/bin/env python3
import argparse
import json
import math
import re
import shlex
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

G_MSUN_PC_MYR = 0.00449830997959438
VERSION_RE = re.compile(r"^Version:\s*(.*)$")
TIME_RE = re.compile(r"^Time:\s*([0-9eE+\-.]+)")
ENERGY_HEADER_RE = re.compile(r"^Energy:\s+(.*)$")
ENERGY_PHYSIC_RE = re.compile(r"^Physic:\s+(.*)$")
INPUT_UNIT_RE = re.compile(r"^Input data unit:\s*(\d+)")
UNIT_SET_RE = re.compile(r"^----- Unit set\s*(\d+)\s*:\s*(.*)-----$")
TREE_DT_RE = re.compile(r"^Tree time step\s*=\s*([0-9eE+\-.]+)$")
OUTPUT_DT_RE = re.compile(r"^Output time step\s*=\s*([0-9eE+\-.]+)$")

Point = Tuple[float, float]


def format_tick(value: float) -> str:
    if value == 0.0:
        return "0"
    av = abs(value)
    if av >= 1e3 or av < 1e-3:
        return f"{value:.2e}"
    return f"{value:.3g}"


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


def parse_time_metric_series(log_text: str, metric: str) -> List[Point]:
    lines = log_text.splitlines()
    series: List[Point] = []
    current_time = None
    for idx, line in enumerate(lines):
        tmatch = TIME_RE.match(line.strip())
        if tmatch:
            try:
                current_time = float(tmatch.group(1))
            except ValueError:
                current_time = None
            continue

        hmatch = ENERGY_HEADER_RE.match(line.strip())
        if hmatch and idx + 1 < len(lines) and current_time is not None:
            pmatch = ENERGY_PHYSIC_RE.match(lines[idx + 1].strip())
            if not pmatch:
                continue
            header = hmatch.group(1).split()
            values = pmatch.group(1).split()
            if len(values) < len(header):
                continue
            rec = {header[j]: float(values[j]) for j in range(len(header))}
            series.append((current_time, abs(rec.get(metric, 0.0))))
    return series


def parse_unit_info(log_text: str) -> Dict[str, str]:
    input_unit = "unknown"
    unit_set = "unknown"
    for line in log_text.splitlines():
        sm = line.strip()
        im = INPUT_UNIT_RE.match(sm)
        if im:
            input_unit = im.group(1)
        um = UNIT_SET_RE.match(sm)
        if um:
            unit_set = f"{um.group(1)}: {um.group(2).strip()}"
    return {"input_unit": input_unit, "unit_set": unit_set}


def parse_runtime_steps(log_text: str) -> Dict[str, float]:
    tree_dt = math.nan
    output_dt = math.nan
    for line in log_text.splitlines():
        sm = line.strip()
        tm = TREE_DT_RE.match(sm)
        if tm:
            tree_dt = float(tm.group(1))
        om = OUTPUT_DT_RE.match(sm)
        if om:
            output_dt = float(om.group(1))
    return {"tree_dt": tree_dt, "output_dt": output_dt}


def resolve_status_path(run_record: Dict) -> Path:
    cmd = run_record["command"]
    args_map = parse_command_args(cmd)
    prefix = args_map.get("-f", run_record["run_id"])
    cd_match = re.search(r"\bcd\s+([^&;]+?)\s*&&", cmd)
    if cd_match:
        work_dir = Path(cd_match.group(1).strip())
        if not work_dir.is_absolute():
            work_dir = (Path.cwd() / work_dir).resolve()
        return work_dir / f"{prefix}.status"
    log_path = Path(run_record["output"])
    return log_path.parent / f"{prefix}.status"


def load_status_particles(path: Path, n_particle: int):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
    from analysis.status import Status  # pylint: disable=import-outside-toplevel  # pyright: ignore[reportMissingImports]

    st = Status(N_particle=n_particle)
    with warnings.catch_warnings(record=True) as warn_list:
        warnings.simplefilter("always", category=UserWarning)
        st.fromfile(str(path))
    if warn_list:
        msg = "; ".join(str(item.message) for item in warn_list)
        raise RuntimeError(
            f"Status parse warning for {path} with N_particle={n_particle}: {msg}. "
            "This usually indicates an N_particle mismatch."
        )
    return st


def derive_inner_binary_drift(status_path: Path, n_particle: int) -> Dict[str, object]:
    st = load_status_particles(status_path, n_particle)

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
        "da_rel": [float(x) for x in da_rel.tolist()],
        "de_abs": [float(x) for x in de_abs.tolist()],
        "da_rel_max": float(np.max(da_rel)),
        "de_abs_max": float(np.max(de_abs)),
    }


def derive_periods_from_status(status_path: Path, n_particle: int) -> Dict[str, float]:
    st = load_status_particles(status_path, n_particle)
    if n_particle < 3:
        raise RuntimeError(f"Need at least 3 particles for triple period analysis: {status_path}")

    p0 = st.particles.p0
    p1 = st.particles.p1
    p2 = st.particles.p2

    m0 = float(p0.mass[0])
    m1 = float(p1.mass[0])
    m2 = float(p2.mass[0])

    r01 = np.array([p1.pos[0, 0] - p0.pos[0, 0], p1.pos[0, 1] - p0.pos[0, 1], p1.pos[0, 2] - p0.pos[0, 2]], dtype=float)
    v01 = np.array([p1.vel[0, 0] - p0.vel[0, 0], p1.vel[0, 1] - p0.vel[0, 1], p1.vel[0, 2] - p0.vel[0, 2]], dtype=float)

    mu_in = G_MSUN_PC_MYR * (m0 + m1)
    eps_in = 0.5 * float(np.dot(v01, v01)) - mu_in / float(np.linalg.norm(r01))
    a_in = -mu_in / (2.0 * eps_in)
    p_in = 2.0 * math.pi * math.sqrt(a_in * a_in * a_in / mu_in)

    r_cm = (m0 * np.array([p0.pos[0, 0], p0.pos[0, 1], p0.pos[0, 2]]) + m1 * np.array([p1.pos[0, 0], p1.pos[0, 1], p1.pos[0, 2]])) / (m0 + m1)
    v_cm = (m0 * np.array([p0.vel[0, 0], p0.vel[0, 1], p0.vel[0, 2]]) + m1 * np.array([p1.vel[0, 0], p1.vel[0, 1], p1.vel[0, 2]])) / (m0 + m1)
    r_out = np.array([p2.pos[0, 0], p2.pos[0, 1], p2.pos[0, 2]], dtype=float) - r_cm
    v_out = np.array([p2.vel[0, 0], p2.vel[0, 1], p2.vel[0, 2]], dtype=float) - v_cm

    mu_out = G_MSUN_PC_MYR * (m0 + m1 + m2)
    eps_out = 0.5 * float(np.dot(v_out, v_out)) - mu_out / float(np.linalg.norm(r_out))
    a_out = -mu_out / (2.0 * eps_out)
    p_out = 2.0 * math.pi * math.sqrt(a_out * a_out * a_out / mu_out)

    return {
        "period_inner": float(p_in),
        "period_outer": float(p_out),
    }


def derive_orbital_series_from_status(status_path: Path, n_particle: int) -> Dict[str, List[float]]:
    st = load_status_particles(status_path, n_particle)
    if n_particle < 3:
        raise RuntimeError(f"Need at least 3 particles for triple orbital-series analysis: {status_path}")

    p0 = st.particles.p0
    p1 = st.particles.p1
    p2 = st.particles.p2

    m0 = float(p0.mass[0])
    m1 = float(p1.mass[0])
    m2 = float(p2.mass[0])

    r0 = np.asarray(p0.pos, dtype=float)
    r1 = np.asarray(p1.pos, dtype=float)
    r2 = np.asarray(p2.pos, dtype=float)
    v0 = np.asarray(p0.vel, dtype=float)
    v1 = np.asarray(p1.vel, dtype=float)
    v2 = np.asarray(p2.vel, dtype=float)

    r_in = r1 - r0
    v_in = v1 - v0

    m01 = m0 + m1
    r_cm = (m0 * r0 + m1 * r1) / m01
    v_cm = (m0 * v0 + m1 * v1) / m01
    r_out = r2 - r_cm
    v_out = v2 - v_cm

    def calc_orbit_series(r: np.ndarray, v: np.ndarray, mu: float):
        rnorm = np.linalg.norm(r, axis=1)
        v2 = np.sum(v * v, axis=1)
        eps = 0.5 * v2 - mu / rnorm
        a = -mu / (2.0 * eps)

        h = np.cross(r, v)
        hnorm = np.linalg.norm(h, axis=1)
        e_vec = np.cross(v, h) / mu - r / rnorm[:, None]
        e = np.linalg.norm(e_vec, axis=1)

        cos_i = np.clip(h[:, 2] / np.maximum(hnorm, 1e-30), -1.0, 1.0)
        inc = np.degrees(np.arccos(cos_i))
        return a, e, inc, h

    mu_in = G_MSUN_PC_MYR * m01
    mu_out = G_MSUN_PC_MYR * (m01 + m2)
    a_in, e_in, i_in, h_in = calc_orbit_series(r_in, v_in, mu_in)
    a_out, e_out, i_out, h_out = calc_orbit_series(r_out, v_out, mu_out)

    h_in_norm = np.linalg.norm(h_in, axis=1)
    h_out_norm = np.linalg.norm(h_out, axis=1)
    cos_mut = np.sum(h_in * h_out, axis=1) / np.maximum(h_in_norm * h_out_norm, 1e-30)
    cos_mut = np.clip(cos_mut, -1.0, 1.0)
    i_mutual = np.degrees(np.arccos(cos_mut))

    return {
        "time": [float(x) for x in st.time.tolist()],
        "a_in": [float(x) for x in a_in.tolist()],
        "e_in": [float(x) for x in e_in.tolist()],
        "i_in": [float(x) for x in i_in.tolist()],
        "a_out": [float(x) for x in a_out.tolist()],
        "e_out": [float(x) for x in e_out.tolist()],
        "i_out": [float(x) for x in i_out.tolist()],
        "i_mutual": [float(x) for x in i_mutual.tolist()],
    }


def svg_multi_plot(
    series_map: Dict[str, List[Point]],
    title: str,
    x_label: str,
    y_label: str,
    logx: bool = False,
    logy: bool = True,
) -> str:
    width, height = 900, 500
    left, right, top, bottom = 88, 24, 44, 70
    plot_w = width - left - right
    plot_h = height - top - bottom

    all_points = [p for points in series_map.values() for p in points]
    if not all_points:
        return "<p>No data available.</p>"

    xs = [p[0] for p in all_points]
    ys = [max(p[1], 1e-15) if logy else p[1] for p in all_points]

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
        y = max(v, 1e-15) if logy else v
        return top + (ymax - ty(y)) / (ymax - ymin) * plot_h

    palette = ["#0b74de", "#d81b60", "#2e7d32", "#8e24aa", "#f4511e", "#00897b"]
    polylines: List[str] = []
    legends: List[str] = []

    x0 = min(xs)
    x1 = max(xs)
    if x1 == x0:
        x1 = x0 + 1.0
    x_tick_values = [x0 + (x1 - x0) * i / 5.0 for i in range(6)]

    y0 = min(ys)
    y1 = max(ys)
    y_tick_values: List[float]
    if logy:
        p0 = int(math.floor(math.log10(max(y0, 1e-15))))
        p1 = int(math.ceil(math.log10(max(y1, 1e-15))))
        y_tick_values = [10.0 ** p for p in range(p0, p1 + 1)]
    else:
        if y1 == y0:
            y1 = y0 + 1.0
        y_tick_values = [y0 + (y1 - y0) * i / 5.0 for i in range(6)]

    x_grid: List[str] = []
    x_ticks: List[str] = []
    for value in x_tick_values:
        px = sx(value)
        x_grid.append(f'<line x1="{px:.2f}" y1="{top}" x2="{px:.2f}" y2="{height-bottom}" stroke="#e0e0e0" stroke-dasharray="3,3"/>')
        x_ticks.append(f'<line x1="{px:.2f}" y1="{height-bottom}" x2="{px:.2f}" y2="{height-bottom+6}" stroke="#333"/>')
        x_ticks.append(f'<text x="{px:.2f}" y="{height-bottom+22}" text-anchor="middle" font-size="11" font-family="Arial">{html_escape(format_tick(value))}</text>')

    y_grid: List[str] = []
    y_ticks: List[str] = []
    for value in y_tick_values:
        py = sy(value)
        y_grid.append(f'<line x1="{left}" y1="{py:.2f}" x2="{width-right}" y2="{py:.2f}" stroke="#e0e0e0" stroke-dasharray="3,3"/>')
        y_ticks.append(f'<line x1="{left-6}" y1="{py:.2f}" x2="{left}" y2="{py:.2f}" stroke="#333"/>')
        y_ticks.append(f'<text x="{left-10}" y="{py+4:.2f}" text-anchor="end" font-size="11" font-family="Arial">{html_escape(format_tick(value))}</text>')

    items = list(series_map.items())
    # Draw m1 last and thicker so it is always visible even if curves overlap.
    items = [it for it in items if "m1_pure_sdar_ref" not in it[0]] + [it for it in items if "m1_pure_sdar_ref" in it[0]]

    for idx, (label, points) in enumerate(items):
        color = palette[idx % len(palette)]
        stroke_w = "3.5" if "m1_pure_sdar_ref" in label else "2"
        poly = " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y in points)
        polylines.append(f'<polyline points="{poly}" fill="none" stroke="{color}" stroke-width="{stroke_w}"/>')
        ly = top + 18 + idx * 18
        legends.append(
            f'<line x1="{width-right-310}" y1="{ly}" x2="{width-right-290}" y2="{ly}" stroke="{color}" stroke-width="{stroke_w}"/>'
        )
        legends.append(
            f'<text x="{width-right-285}" y="{ly+4}" font-size="11" font-family="Arial">{html_escape(label)}</text>'
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


def extract_binary_from_command(command: str) -> str:
    tokens = shlex.split(command)
    for idx, token in enumerate(tokens):
        if token.startswith("OMP_NUM_THREADS=") and idx + 1 < len(tokens):
            return tokens[idx + 1]
    return "unknown"


def scenario_file_from_name(scenario: str) -> Path:
    return Path("test/validation/scenarios") / f"{scenario}.json"


def normalize_mode_id(run_id: str) -> str:
    return run_id[:-4] if run_id.endswith("_o15") else run_id


def collect_plot_series(rows: List[Dict]) -> Dict[str, object]:
    energy_total_series: Dict[str, List[Point]] = {}
    energy_pp_series: Dict[str, List[Point]] = {}
    a_in_series: Dict[str, List[Point]] = {}
    e_in_series: Dict[str, List[Point]] = {}
    i_in_series: Dict[str, List[Point]] = {}
    a_out_series: Dict[str, List[Point]] = {}
    e_out_series: Dict[str, List[Point]] = {}
    i_out_series: Dict[str, List[Point]] = {}
    i_mutual_series: Dict[str, List[Point]] = {}
    unit_rows = []

    for row in rows:
        label = f"{row['run_id']} | tt={row['tt_switch']}"
        log_text = Path(row["log_path"]).read_text(encoding="utf-8", errors="replace")
        energy_total_series[label] = parse_time_metric_series(log_text, "Error/Total")
        energy_pp_series[label] = parse_time_metric_series(log_text, "Error_PP")
        orbit = derive_orbital_series_from_status(Path(row["status_path"]), int(row["max_n_real_glb"]))
        tvals = orbit["time"]
        row["a_in0"] = float(orbit["a_in"][0]) if orbit["a_in"] else math.nan
        row["e_in0"] = float(orbit["e_in"][0]) if orbit["e_in"] else math.nan
        row["i_in0"] = float(orbit["i_in"][0]) if orbit["i_in"] else math.nan
        row["a_out0"] = float(orbit["a_out"][0]) if orbit["a_out"] else math.nan
        row["e_out0"] = float(orbit["e_out"][0]) if orbit["e_out"] else math.nan
        row["i_out0"] = float(orbit["i_out"][0]) if orbit["i_out"] else math.nan
        row["i_mutual0"] = float(orbit["i_mutual"][0]) if orbit["i_mutual"] else math.nan
        a_in_series[label] = list(zip(tvals, orbit["a_in"]))
        e_in_series[label] = list(zip(tvals, orbit["e_in"]))
        i_in_series[label] = list(zip(tvals, orbit["i_in"]))
        a_out_series[label] = list(zip(tvals, orbit["a_out"]))
        e_out_series[label] = list(zip(tvals, orbit["e_out"]))
        i_out_series[label] = list(zip(tvals, orbit["i_out"]))
        i_mutual_series[label] = list(zip(tvals, orbit["i_mutual"]))
        units = parse_unit_info(log_text)
        unit_rows.append((row["run_id"], units["input_unit"], units["unit_set"]))

    return {
        "energy_total_series": energy_total_series,
        "energy_pp_series": energy_pp_series,
        "a_in_series": a_in_series,
        "e_in_series": e_in_series,
        "i_in_series": i_in_series,
        "a_out_series": a_out_series,
        "e_out_series": e_out_series,
        "i_out_series": i_out_series,
        "i_mutual_series": i_mutual_series,
        "unit_rows": unit_rows,
    }


def build_orbit_figure_section(title_prefix: str, figure_offset: int, series: Dict[str, Dict[str, List[Point]]]) -> str:
    return f"""
<h2>{html_escape(title_prefix)}: Energy Error</h2>
<h3>Figure {figure_offset}: Total Relative Energy Error vs Time</h3>
{svg_multi_plot(series['energy_total_series'], f'{title_prefix}: |Error/Total|(t)', 'time [Myr]', '|Error/Total|', logx=False, logy=True)}

<h3>Figure {figure_offset + 1}: Short-range Energy Error (PP) vs Time</h3>
{svg_multi_plot(series['energy_pp_series'], f'{title_prefix}: |Error_PP|(t)', 'time [Myr]', '|Error_PP|', logx=False, logy=True)}

<h2>{html_escape(title_prefix)}: Orbital Evolution of Inner Orbit</h2>
<h3>Figure {figure_offset + 2}: Inner Semi-major Axis a_in(t)</h3>
{svg_multi_plot(series['a_in_series'], f'{title_prefix}: inner orbit a(t)', 'time [Myr]', 'a_in [pc]', logx=False, logy=False)}

<h3>Figure {figure_offset + 3}: Inner Eccentricity e_in(t)</h3>
{svg_multi_plot(series['e_in_series'], f'{title_prefix}: inner orbit e(t)', 'time [Myr]', 'e_in', logx=False, logy=False)}

<h3>Figure {figure_offset + 4}: Inner Inclination i_in(t)</h3>
{svg_multi_plot(series['i_in_series'], f'{title_prefix}: inner orbit inclination i(t)', 'time [Myr]', 'i_in [deg]', logx=False, logy=False)}

<h2>{html_escape(title_prefix)}: Orbital Evolution of Outer Orbit</h2>
<h3>Figure {figure_offset + 5}: Outer Semi-major Axis a_out(t)</h3>
{svg_multi_plot(series['a_out_series'], f'{title_prefix}: outer orbit a(t)', 'time [Myr]', 'a_out [pc]', logx=False, logy=False)}

<h3>Figure {figure_offset + 6}: Outer Eccentricity e_out(t)</h3>
{svg_multi_plot(series['e_out_series'], f'{title_prefix}: outer orbit e(t)', 'time [Myr]', 'e_out', logx=False, logy=False)}

<h3>Figure {figure_offset + 7}: Outer Inclination i_out(t)</h3>
{svg_multi_plot(series['i_out_series'], f'{title_prefix}: outer orbit inclination i(t)', 'time [Myr]', 'i_out [deg]', logx=False, logy=False)}

<h3>Figure {figure_offset + 8}: Mutual Inclination i_mutual(t)</h3>
{svg_multi_plot(series['i_mutual_series'], f'{title_prefix}: mutual inclination i_mutual(t)', 'time [Myr]', 'i_mutual [deg]', logx=False, logy=False)}
"""


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate T4 triple mode comparison HTML report")
    parser.add_argument("--report", default="test/validation/out/report.t4.triple.json")
    parser.add_argument("--report-control", default="", help="Optional control-group report (e.g. outer a=1.5)")
    parser.add_argument("--scenario", default="t4_tree_hard_from_triple")
    parser.add_argument("--scenario-file", default="")
    parser.add_argument("--output", default="test/validation/out/t4_triple_mode_summary.html")
    args = parser.parse_args()

    report = load_json(Path(args.report))
    scenario_path = Path(args.scenario_file) if args.scenario_file else scenario_file_from_name(args.scenario)
    scenario_data = load_json(scenario_path) if scenario_path.exists() else {}

    runs = [r for r in report.get("runs", []) if r.get("scenario") == args.scenario]
    if not runs:
        raise RuntimeError(f"No runs found for scenario '{args.scenario}' in report: {args.report}")
    checks = [c for c in report.get("checks", []) if c.get("scenario") == args.scenario]

    rows = []
    period_inner = math.nan
    period_outer = math.nan
    for run in runs:
        args_map = parse_command_args(run["command"])
        status_path = resolve_status_path(run)
        metrics = run.get("metrics", {})
        log_text = Path(run["output"]).read_text(encoding="utf-8", errors="replace")
        runtime_steps = parse_runtime_steps(log_text)
        dt_soft = float(metrics.get("dt_soft", math.nan))
        if not math.isfinite(dt_soft):
            dt_soft = runtime_steps["tree_dt"]
        tree_dt = runtime_steps["tree_dt"]
        dt_le_tree = (dt_soft <= tree_dt + 1e-15) if (math.isfinite(dt_soft) and math.isfinite(tree_dt)) else False
        n_particle = int(float(metrics.get("max_n_real_glb", 0.0))) if "max_n_real_glb" in metrics else 0
        if n_particle < 2:
            n_particle = 3
        drift = derive_inner_binary_drift(status_path, n_particle)
        if not math.isfinite(period_outer):
            periods = derive_periods_from_status(status_path, n_particle)
            period_inner = periods["period_inner"]
            period_outer = periods["period_outer"]
        t_end = float(drift["time"][-1]) if drift["time"] else 0.0
        rows.append(
            {
                "run_id": run["run_id"],
                "mode": run["run_id"],
                "tt_switch": int(float(metrics.get("tt_switch", 0.0))),
                "tt_nstep": int(float(metrics.get("tt_nstep", 0.0))) if "tt_nstep" in metrics else 0,
                "r_group": float(metrics.get("r_group", math.nan)),
                "r_search_group": float(metrics.get("r_search_group", math.nan)),
                "rout": float(metrics.get("rout", math.nan)),
                "dt_soft": dt_soft,
                "tree_dt": tree_dt,
                "dt_le_tree": dt_le_tree,
                "max_abs_error_over_total": float(metrics.get("max_abs_error_over_total", 0.0)),
                "max_abs_error_pp": float(metrics.get("max_abs_error_pp", 0.0)),
                "max_artificial_particles_glb": int(float(metrics.get("max_artificial_particles_glb", 0.0))),
                "max_n_real_glb": int(float(metrics.get("max_n_real_glb", 0.0))),
                "max_n_all_glb": int(float(metrics.get("max_n_all_glb", 0.0))),
                "da_rel_max": drift["da_rel_max"],
                "de_abs_max": drift["de_abs_max"],
                "t_end": t_end,
                "inner_period_covered": (t_end / period_inner) if (period_inner > 0.0) else math.nan,
                "outer_period_covered": (t_end / period_outer) if (period_outer > 0.0) else math.nan,
                "time": drift["time"],
                "da_rel": drift["da_rel"],
                "de_abs": drift["de_abs"],
                "command": run["command"],
                "log_path": run["output"],
                "status_path": str(status_path),
                "args_map": args_map,
            }
        )

    order = [
        "m1_pure_sdar_ref",
        "m2_hard_no_tt",
        "m3_hard_tt",
        "m4_tree_no_tt",
        "m5_tree_tt",
    ]
    rows = sorted(rows, key=lambda x: order.index(x["run_id"]) if x["run_id"] in order else 999)
    log_path0 = Path(rows[0]["log_path"])
    metadata = parse_log_metadata(log_path0)

    series = collect_plot_series(rows)
    energy_total_series = series["energy_total_series"]
    energy_pp_series = series["energy_pp_series"]
    a_in_series = series["a_in_series"]
    e_in_series = series["e_in_series"]
    i_in_series = series["i_in_series"]
    a_out_series = series["a_out_series"]
    e_out_series = series["e_out_series"]
    i_out_series = series["i_out_series"]
    i_mutual_series = series["i_mutual_series"]
    unit_rows = series["unit_rows"]

    run_rows_html = "\n".join(
        "<tr>"
        f"<td>{html_escape(r['run_id'])}</td><td>{r['tt_switch']}</td><td>{r['tt_nstep']}</td><td>{r['r_group']:.6e}</td><td>{r['r_search_group']:.6e}</td><td>{r['rout']:.6e}</td><td>{r['dt_soft']:.6e}</td><td>{r['tree_dt']:.6e}</td><td>{'YES' if r['dt_le_tree'] else 'NO'}</td><td>{r['a_in0']:.6e}</td><td>{r['e_in0']:.6e}</td><td>{r['i_in0']:.3f}</td><td>{r['a_out0']:.6e}</td><td>{r['e_out0']:.6e}</td><td>{r['i_out0']:.3f}</td><td>{r['i_mutual0']:.3f}</td><td>{r['t_end']:.6e}</td><td>{r['inner_period_covered']:.2f}</td><td>{r['outer_period_covered']:.2f}</td><td>{r['max_abs_error_over_total']:.6e}</td><td>{r['max_abs_error_pp']:.6e}</td><td>{r['da_rel_max']:.6e}</td><td>{r['de_abs_max']:.6e}</td><td>{r['max_n_real_glb']}</td><td>{r['max_n_all_glb']}</td><td>{r['max_artificial_particles_glb']}</td>"
        "</tr>"
        for r in rows
    )

    base_figures_html = build_orbit_figure_section("Base Group", 1, series)
    control_compare_html = ""
    control_figures_html = ""
    if args.report_control:
        control_report = load_json(Path(args.report_control))
        control_runs = control_report.get("runs", [])
        if isinstance(control_runs, list) and control_runs:
            control_rows = []
            for rr in control_runs:
                rmet = rr.get("metrics", {})
                status_path = resolve_status_path(rr)
                log_text = Path(rr["output"]).read_text(encoding="utf-8", errors="replace")
                runtime_steps = parse_runtime_steps(log_text)
                dt_soft = float(rmet.get("dt_soft", math.nan))
                if not math.isfinite(dt_soft):
                    dt_soft = runtime_steps["tree_dt"]
                tree_dt = runtime_steps["tree_dt"]
                dt_le_tree = (dt_soft <= tree_dt + 1e-15) if (math.isfinite(dt_soft) and math.isfinite(tree_dt)) else False
                control_rows.append(
                    {
                        "run_id": rr.get("run_id", ""),
                        "mode_id": normalize_mode_id(str(rr.get("run_id", ""))),
                        "tt_switch": int(float(rmet.get("tt_switch", 0.0))),
                        "tt_nstep": int(float(rmet.get("tt_nstep", 0.0))) if "tt_nstep" in rmet else 0,
                        "r_group": float(rmet.get("r_group", math.nan)),
                        "r_search_group": float(rmet.get("r_search_group", math.nan)),
                        "rout": float(rmet.get("rout", math.nan)),
                        "dt_soft": dt_soft,
                        "tree_dt": tree_dt,
                        "dt_le_tree": dt_le_tree,
                        "err_tot": float(rmet.get("max_abs_error_over_total", math.nan)),
                        "err_pp": float(rmet.get("max_abs_error_pp", math.nan)),
                        "da": float(rmet.get("max_abs_a_frac_drift", rmet.get("max_rel_drift_semi", rmet.get("da_rel_max", math.nan)))),
                        "de": float(rmet.get("max_abs_e_abs_drift", rmet.get("max_abs_drift_ecc", rmet.get("de_abs_max", math.nan)))),
                        "max_art": float(rmet.get("max_artificial_particles_glb", math.nan)),
                        "max_abs_error_over_total": float(rmet.get("max_abs_error_over_total", math.nan)),
                        "max_abs_error_pp": float(rmet.get("max_abs_error_pp", math.nan)),
                        "max_artificial_particles_glb": int(float(rmet.get("max_artificial_particles_glb", 0.0))),
                        "max_n_real_glb": int(float(rmet.get("max_n_real_glb", 0.0))),
                        "max_n_all_glb": int(float(rmet.get("max_n_all_glb", 0.0))),
                        "command": rr["command"],
                        "log_path": rr["output"],
                        "status_path": str(status_path),
                    }
                )

            control_rows = sorted(control_rows, key=lambda x: order.index(x["mode_id"]) if x["mode_id"] in order else 999)
            control_series = collect_plot_series(control_rows)
            control_figures_html = build_orbit_figure_section("Control Group (Outer a=1.5)", 10, control_series)

            base_by_mode = {
                normalize_mode_id(str(r["run_id"])): r
                for r in rows
            }
            ctrl_by_mode = {r["mode_id"]: r for r in control_rows}

            lines = []
            for mode in ["m1_pure_sdar_ref", "m2_hard_no_tt", "m3_hard_tt", "m4_tree_no_tt", "m5_tree_tt"]:
                if mode not in base_by_mode or mode not in ctrl_by_mode:
                    continue
                b = base_by_mode[mode]
                c = ctrl_by_mode[mode]
                lines.append(
                    "<tr>"
                    f"<td>{html_escape(mode)}</td>"
                    f"<td>{b['max_abs_error_over_total']:.6e}</td>"
                    f"<td>{c['err_tot']:.6e}</td>"
                    f"<td>{b['max_abs_error_pp']:.6e}</td>"
                    f"<td>{c['err_pp']:.6e}</td>"
                    f"<td>{b['da_rel_max']:.6e}</td>"
                    f"<td>{c['da']:.6e}</td>"
                    f"<td>{b['de_abs_max']:.6e}</td>"
                    f"<td>{c['de']:.6e}</td>"
                    f"<td>{b['max_artificial_particles_glb']}</td>"
                    f"<td>{c['max_art']:.0f}</td>"
                    "</tr>"
                )
            if lines:
                control_compare_html = (
                    "<h2>Control Group Comparison (Outer a=1.5)</h2>"
                    "<table><tr><th>mode</th><th>base |Error/Total|</th><th>control |Error/Total|</th><th>base |Error_PP|</th><th>control |Error_PP|</th><th>base |Δa/a0|</th><th>control |Δa/a0|</th><th>base |Δe|</th><th>control |Δe|</th><th>base max artificial</th><th>control max artificial</th></tr>"
                    + "\n".join(lines)
                    + "</table>"
                )

    command_rows_html = "\n".join(
        "<tr>"
        f"<td>{html_escape(r['run_id'])}</td><td><code>{html_escape(r['command'])}</code></td><td><code>{html_escape(r['status_path'])}</code></td><td><code>{html_escape(r['log_path'])}</code></td>"
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

    feature_list = metadata.get("features", [])
    if not isinstance(feature_list, list):
        feature_list = []
    feature_html = "".join(f"<li>{html_escape(str(item))}</li>" for item in feature_list)

    html = f"""<!doctype html>
<html lang=\"en\"><head><meta charset=\"utf-8\"><title>T4 Triple Mode Comparison Report</title>
<style>
body{{font-family:Arial,Helvetica,sans-serif;margin:20px;line-height:1.45;}}
h1,h2,h3{{color:#1f2d3d;}}
table{{border-collapse:collapse;width:100%;margin:8px 0 16px 0;}}
th,td{{border:1px solid #cfd8dc;padding:6px 8px;font-size:13px;vertical-align:top;}}
th{{background:#f5f7fa;}}
.note{{background:#eef7ff;border-left:4px solid #0b74de;padding:8px 10px;margin:10px 0;}}
code{{white-space:pre-wrap;word-break:break-all;}}
</style></head><body>
<h1>T4: Triple Algorithm-Mode Comparison (SDAR/Hermite/Tree with Tidal Tensor)</h1>
<div class=\"note\">
Five-mode comparison from notebook-style hierarchical triple scales: <code>m1_pure_sdar_ref</code>, <code>m2_hard_no_tt</code>, <code>m3_hard_tt</code>, <code>m4_tree_no_tt</code>, <code>m5_tree_tt</code>.
Tidal tensor activation is validated by checking <code>N_all(glb)-N_real(glb) &gt; 0</code> when <code>--tt-switch=1</code>.
</div>

{base_figures_html}

{control_figures_html}

<h2>Unit Consistency Check</h2>
<table><tr><th>run_id</th><th>Input data unit</th><th>Unit set line</th></tr>
{''.join(f'<tr><td>{html_escape(rid)}</td><td>{html_escape(iu)}</td><td>{html_escape(us)}</td></tr>' for rid, iu, us in unit_rows)}
</table>

<h2>Period Coverage</h2>
<ul>
<li>Inner period (from initial status): <code>{period_inner:.6e}</code> Myr</li>
<li>Outer period (from initial status): <code>{period_outer:.6e}</code> Myr</li>
<li>Configured end time: <code>{rows[0]['t_end']:.6e}</code> Myr, equivalent to <code>{rows[0]['outer_period_covered']:.2f}</code> outer periods</li>
</ul>

<h2>Run Summary</h2>
<table><tr><th>run_id</th><th>tt_switch</th><th>tt_nstep</th><th>r_group</th><th>r_search_group</th><th>r_out</th><th>dt_soft</th><th>tree_dt(auto)</th><th>dt_soft&lt;=tree_dt</th><th>a_in(0)</th><th>e_in(0)</th><th>i_in(0) [deg]</th><th>a_out(0)</th><th>e_out(0)</th><th>i_out(0) [deg]</th><th>i_mutual(0) [deg]</th><th>t_end</th><th>inner periods</th><th>outer periods</th><th>max |Error/Total|</th><th>max |Error_PP|</th><th>max |Δa/a0| (dynamic)</th><th>max |Δe| (dynamic)</th><th>max N_real</th><th>max N_all</th><th>max artificial</th></tr>
{run_rows_html}
</table>

{control_compare_html}

<h2>PeTar Version and Executable Info</h2>
<ul>
<li>Executable: <code>{html_escape(extract_binary_from_command(rows[0]['command']))}</code></li>
<li>PeTar version: <code>{html_escape(str(metadata.get('version', 'unknown')))}</code></li>
</ul>

<h3>Build/Runtime Features Parsed from Logs</h3>
<ul>{feature_html}</ul>

<h2>Initialization Commands</h2>
<table><tr><th>#</th><th>exact setup command</th></tr>
{setup_rows_html}
</table>

<h2>All Executed Run Commands</h2>
<table><tr><th>run_id</th><th>exact command</th><th>status path</th><th>log path</th></tr>
{command_rows_html}
</table>

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
