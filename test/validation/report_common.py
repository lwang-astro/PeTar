#!/usr/bin/env python3
"""Shared helpers for PeTar validation HTML reports.

Every report script imports its generic helpers from this module instead of
defining local copies; see test/validation/REPORT_DESIGN.md for the report
design rules (paper-style sections, figure encoding, caption conventions).
"""
import json
import math
import os
import re
import shlex
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np  # used by orbital_elements_from_two_body_rows

Point = Tuple[float, float]
VERSION_RE = re.compile(r"^Version:\s*(.*)$")
G_MSUN_PC_MYR = 0.00449830997959438  # pc^3 Msun^-1 Myr^-2 (same as make_ic.py)

# Paper-style stylesheet shared by all validation reports (see REPORT_DESIGN.md).
PAPER_CSS = """
 body { font-family: Georgia, 'Times New Roman', serif; margin: 32px auto; max-width: 960px; color: #1a1a1a; line-height: 1.5; }
 h1 { font-size: 22px; text-align: center; margin-bottom: 4px; }
 h2 { font-size: 17px; margin-top: 30px; border-bottom: 1px solid #999; padding-bottom: 3px; }
 h3 { font-size: 15px; margin-top: 18px; }
 .subtitle { text-align: center; color: #555; font-size: 14px; margin-bottom: 14px; }
 .verdict { display: inline-block; padding: 5px 14px; font-weight: bold; border-radius: 4px; font-family: Arial; }
 .verdict.pass { background: #2e7d32; color: #fff; }
 .verdict.fail { background: #c62828; color: #fff; }
 .eq { text-align: center; margin: 10px 0; font-style: italic; font-size: 15px; }
 .eq .tag { float: right; color: #777; font-style: normal; }
 .caption { font-size: 12.5px; color: #444; margin: 6px 0 18px 0; }
 table.data { border-collapse: collapse; font-size: 13.5px; margin: 8px 0 2px 0; font-family: Arial; }
 table.data th, table.data td { border: 1px solid #999; padding: 4px 9px; text-align: left; vertical-align: top; }
 table.data th { background: #eee; }
 td.pass { color: #2e7d32; font-weight: bold; }
 td.fail { color: #c62828; font-weight: bold; }
 td.msg { color: #555; }
 .cards { display: flex; flex-wrap: wrap; gap: 10px; margin: 10px 0 2px 0; }
 .card { border: 1px solid #bbb; border-radius: 5px; padding: 8px 14px; min-width: 150px; font-family: Arial; }
 .card-label { font-size: 12px; color: #666; }
 .card-value { font-size: 19px; font-weight: bold; margin-top: 3px; }
 .note { font-size: 13.5px; color: #333; margin: 8px 0 2px 0; }
 pre.cmd { background: #f6f6f6; border: 1px solid #ddd; border-radius: 4px; padding: 8px 12px; font-size: 12px; font-family: 'Courier New', monospace; overflow-x: auto; white-space: pre-wrap; word-break: break-all; }
 .meta { color: #777; font-size: 12px; margin-top: 26px; border-top: 1px solid #ccc; padding-top: 8px; }
 i, .sym { font-family: Georgia, serif; }
"""


def git_commit(repo_root: Path) -> str:
    import subprocess

    try:
        r = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=repo_root,
            capture_output=True,
            text=True,
            check=True,
        )
        return r.stdout.strip()
    except Exception:
        return "unknown"


# configure flags for the feature tokens used by the validation suite
_FEATURE_CONFIGURE_FLAGS = {
    "64b": "--enable-64b",
    "bse": "--with-interrupt=bse",
    "merger": "--with-interrupt=merger",
    "mobse": "--with-interrupt=mobse",
    "bseEmp": "--with-interrupt=bseEmp",
    "sevn": "--with-interrupt=sevn",
    "dsm": "--with-interrupt=dsm",
    "galpy": "--with-external=galpy",
    "agama": "--with-external=agama",
}


def _petar_select(require: List[str], optional: List[str]) -> "Optional[Path]":
    """Run petar.select and return the resolved selected binary (symlink
    target), or None if the selection failed."""
    import shutil
    import subprocess

    sel = shutil.which("petar.select")
    if sel is None:
        return None
    cmd = [sel]
    if require:
        cmd += ["--require", ",".join(require)]
    if optional:
        cmd += ["--optional", ",".join(optional)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        return None
    bin_dir = Path(sel).resolve().parent
    link = bin_dir / "petar"
    if not link.exists():
        return None
    target = bin_dir / os.readlink(link)
    return target if target.exists() else None


def _petar_binary_version(binary: Path) -> str:
    """Installed binary version tag (the prefix before '_'), e.g. '1804e'."""
    import subprocess

    try:
        r = subprocess.run(
            [str(binary), "-h"], capture_output=True, text=True, timeout=60
        )
    except Exception:
        return ""
    m = re.search(r"^Version:\s*(\S+)", r.stdout + r.stderr, re.MULTILINE)
    return m.group(1).split("_", 1)[0] if m else ""


def select_or_install_petar(
    require: Optional[List[str]] = None,
    optional: Optional[List[str]] = None,
    repo_root: Optional[Path] = None,
) -> str:
    """Select a petar binary by feature via petar.select, following the
    validation-suite rule: if the selection fails or the installed binary's
    version does not match the repository VERSION, rebuild and reinstall the
    configured family, then re-select. Returns the resolved selected binary
    path (never the bare `petar` symlink name).

    Tests must call this instead of hardcoding binary suffixes or using the
    default `petar`.
    """
    import shutil
    import subprocess

    require = require or []
    optional = optional or []
    repo_root = repo_root or Path(__file__).resolve().parents[2]

    def select() -> Optional[Path]:
        return _petar_select(require, optional)

    def repo_version() -> str:
        vf = repo_root / "VERSION"
        return vf.read_text().strip().split("_", 1)[0] if vf.exists() else ""

    def install() -> None:
        flags = []
        for feat in require:
            flag = _FEATURE_CONFIGURE_FLAGS.get(feat)
            if flag and flag not in flags:
                flags.append(flag)
        configure = "./configure" + ("" if not flags else " " + " ".join(flags))
        ncpu = os.cpu_count() or 4
        print(f"[petar.select] selection failed or version stale; reinstalling: {configure} && make -j{ncpu} install")
        subprocess.run(configure, shell=True, cwd=repo_root, check=True)
        subprocess.run(f"make -j{ncpu}", shell=True, cwd=repo_root, check=True)
        subprocess.run("make install", shell=True, cwd=repo_root, check=True)

    binary = select()
    want = repo_version()
    if binary is not None and want and _petar_binary_version(binary) == want:
        return str(binary)
    install()
    binary = select()
    if binary is None or (want and _petar_binary_version(binary) != want):
        raise RuntimeError(
            f"petar.select still cannot provide an up-to-date binary "
            f"(require={require}, optional={optional}, repo VERSION={want})"
        )
    return str(binary)


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


def parse_log_version(log_path: Path) -> str:
    text = log_path.read_text(encoding="utf-8", errors="replace")
    for line in text.splitlines():
        m = VERSION_RE.match(line.strip())
        if m:
            return m.group(1).strip()
    return "unknown"


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
            if parts:
                rows.append([float(x) for x in parts])
    return rows


def format_tick(value: float) -> str:
    if value == 0.0:
        return "0"
    av = abs(value)
    if av >= 1e3 or av < 1e-3:
        return f"{value:.2e}"
    return f"{value:.3g}"


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


TIME_RE = re.compile(r"^Time:\s*([0-9eE+\-.]+)")
ENERGY_HEADER_RE = re.compile(r"^Energy:\s+(.*)$")
ENERGY_PHYSIC_RE = re.compile(r"^Physic:\s+(.*)$")


def parse_time_error_series(log_path: Path, metric: str = "Error/Total") -> List[Point]:
    """(t, |metric|(t)) samples from a run log's Time/Energy/Physic table."""
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
        if hmatch and i + 1 < len(lines) and current_time is not None:
            pmatch = ENERGY_PHYSIC_RE.match(lines[i + 1].strip())
            if not pmatch:
                continue
            header = hmatch.group(1).split()
            values = pmatch.group(1).split()
            if len(values) < len(header):
                continue
            rec = {header[j]: float(values[j]) for j in range(len(header))}
            series.append((current_time, abs(rec.get(metric, 0.0))))
    return series



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


def _log_floor(ys: List[float]) -> float:
    """Clip floor for log-y plots: a decade below the smallest positive value,
    never above the historical 1e-9 default (data below 1e-9 stays visible)."""
    pos = [v for v in ys if v > 0.0]
    if not pos:
        return 1e-9
    return min(1e-9, 10.0 ** (math.floor(math.log10(min(pos))) - 1))


def _svg_frame(
    width: int,
    height: int,
    left: int,
    right: int,
    top: int,
    bottom: int,
    xs: List[float],
    ys: List[float],
    logx: bool = False,
    logy: bool = True,
):
    plot_w = width - left - right
    plot_h = height - top - bottom
    floor = _log_floor(ys) if logy else 1e-9
    tx = (lambda v: math.log10(max(v, floor))) if logx else (lambda v: v)
    ty = (lambda v: math.log10(max(v, floor))) if logy else (lambda v: v)
    xtv = [tx(v) for v in xs]
    xmin, xmax = min(xtv), max(xtv)
    if xmax == xmin:
        xmax = xmin + 1.0
    ytv = [ty(v) for v in ys]
    ymin, ymax = min(ytv), max(ytv)
    if ymax == ymin:
        ymax = ymin + 1.0

    def sx(v: float) -> float:
        return left + (tx(v) - xmin) / (xmax - xmin) * plot_w

    def sy(v: float) -> float:
        return top + (ymax - ty(v)) / (ymax - ymin) * plot_h

    return sx, sy, xmin, xmax, plot_w, plot_h


# 10-stop approximation of matplotlib's viridis ramp (perceptual, colorblind-safe).
_VIRIDIS_STOPS: List[Tuple[int, int, int]] = [
    (68, 1, 84), (72, 40, 120), (62, 74, 137), (49, 104, 142), (38, 130, 142),
    (31, 158, 137), (53, 183, 121), (109, 205, 89), (180, 222, 44), (253, 231, 37),
]


def gradient_colors(n: int) -> List[str]:
    """Return n hex colors sampled evenly along the viridis ramp.

    Used when a plot encodes an ordered sweep (e.g. dt from fine to coarse) in
    color: adjacent palette entries repeat for >6 runs, a gradient keeps every
    curve distinguishable while making the ordering visible at a glance.
    """
    if n <= 0:
        return []
    stops = _VIRIDIS_STOPS
    if n == 1:
        return ["#%02x%02x%02x" % stops[len(stops) // 2]]
    colors = []
    for i in range(n):
        pos = i * (len(stops) - 1) / (n - 1)
        k = min(int(pos), len(stops) - 2)
        f = pos - k
        c = tuple(round(stops[k][j] + (stops[k + 1][j] - stops[k][j]) * f) for j in range(3))
        colors.append("#%02x%02x%02x" % c)
    return colors


def svg_curves_plot(
    curves: List[Tuple[str, List[Point], str, str]],
    legend_models: List[Tuple[str, str]],
    title: str,
    x_label: str,
    y_label: str,
    hlines: List[Tuple[float, str]] = [],
    bands: List[Tuple[str, List[Tuple[float, float]], str]] = [],
    max_x_ticks: int = 8,
    logx: bool = False,
    logy: bool = True,
    legend_quantities: Optional[List[Tuple[str, str]]] = None,
    legend_headers: Tuple[str, str] = ("model (color)", "quantity (line)"),
) -> str:
    """Paper-style plot: color encodes the model (run) and dash pattern the
    quantity (current style): curves are (label, points, color, dash) with
    dash == "" for solid; legend_models lists (model label, color) and
    legend_quantities lists (quantity label, dash) for the second legend column
    (default: instantaneous/cumulative); legend_headers overrides the two
    column titles when the encoding roles differ. The legend is a dedicated
    strip below the plot whose height adapts to the entry count, so it can
    never overlap the curves. bands: (label, [(t0, t1)...], color) draw
    translucent full-height rectangles marking active intervals. logx/logy
    select log axes (defaults: linear x, log y).

    See REPORT_DESIGN.md for the figure-encoding rules.
    """
    width, height = 940, 500
    left, right, top, bottom = 88, 24, 44, 70

    all_points = [p for _label, points, _c, _d in curves for p in points]
    if not all_points:
        return "<p>No data available.</p>"
    xs = [p[0] for p in all_points]
    ys = [p[1] for p in all_points]
    sx, sy, _, _, _, plot_h = _svg_frame(
        width, height, left, right, top, bottom, xs, ys, logx=logx, logy=logy
    )
    rxmin, rxmax = min(xs), max(xs)
    rymin, rymax = min(ys), max(ys)

    parts: List[str] = []
    for _label, band_list, band_color in bands:
        for t0, t1 in band_list:
            if t1 < rxmin or t0 > rxmax:
                continue
            x0, x1 = sx(max(t0, rxmin)), sx(min(t1, rxmax))
            parts.append(
                f'<rect x="{x0:.2f}" y="{top}" width="{x1-x0:.2f}" height="{plot_h:.2f}" '
                f'fill="{band_color}" fill-opacity="0.13" stroke="none"/>'
            )
    for _label, points, color, dash in curves:
        poly = " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y in points)
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        parts.append(f'<polyline points="{poly}" fill="none" stroke="{color}" stroke-width="2"{dash_attr}/>')

    for idx, (value, label) in enumerate(hlines):
        py = sy(value)
        parts.append(f'<line x1="{left}" y1="{py:.2f}" x2="{width-right}" y2="{py:.2f}" stroke="#333" stroke-width="1.2" stroke-dasharray="2,4"/>')
        ly = height - bottom - 10 - idx * 16
        parts.append(f'<text x="{left+6}" y="{ly}" font-size="11" font-family="Arial">{html_escape(label)}</text>')

    # legend strip below the plot: color = model, dash = quantity
    quantities = legend_quantities if legend_quantities is not None else [("instantaneous", ""), ("cumulative", "7,4")]
    n_tt_legend = 1 if bands else 0
    n_legend_rows = max(len(legend_models), len(quantities)) + n_tt_legend
    svg_h = height + 20 + (17 + 17 * max(n_legend_rows, 2) + 10)
    leg_y = height + 20
    col1_x = width - right - 500
    col2_x = width - right - 240
    legends: List[str] = [f'<text x="{col1_x}" y="{leg_y}" font-size="11" font-family="Arial" font-weight="bold">{html_escape(legend_headers[0])}</text>']
    for i, (label, color) in enumerate(legend_models):
        ly = leg_y + 17 + i * 17
        legends.append(f'<line x1="{col1_x}" y1="{ly}" x2="{col1_x+26}" y2="{ly}" stroke="{color}" stroke-width="3"/>')
        legends.append(f'<text x="{col1_x+31}" y="{ly+4}" font-size="11" font-family="Arial">{html_escape(label)}</text>')
    legends.append(f'<text x="{col2_x}" y="{leg_y}" font-size="11" font-family="Arial" font-weight="bold">{html_escape(legend_headers[1])}</text>')
    for i, (label, dash) in enumerate(quantities):
        ly = leg_y + 17 + i * 17
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        legends.append(f'<line x1="{col2_x}" y1="{ly}" x2="{col2_x+26}" y2="{ly}" stroke="#555" stroke-width="2"{dash_attr}/>')
        legends.append(f'<text x="{col2_x+31}" y="{ly+4}" font-size="11" font-family="Arial">{html_escape(label)}</text>')
    if bands:
        ly = leg_y + 17 + len(quantities) * 17 + 2
        legends.append(f'<rect x="{col2_x}" y="{ly-8}" width="26" height="11" fill="#555" fill-opacity="0.18"/>')
        legends.append(f'<text x="{col2_x+31}" y="{ly+2}" font-size="11" font-family="Arial">active intervals (run color)</text>')

    x_unique = sorted(set(xs))
    if len(x_unique) <= max_x_ticks:
        x_tick_values = x_unique
    elif logx:
        px0 = math.floor(math.log10(max(min(v for v in x_unique if v > 0.0), 1e-9)))
        px1 = math.ceil(math.log10(max(max(x_unique), 1e-9)))
        cand = [10.0 ** p for p in range(int(px0), int(px1) + 1)]
        x_tick_values = [v for v in cand if x_unique[0] <= v <= x_unique[-1]]
        if len(x_tick_values) > max_x_ticks:
            idx_list = [round(i * (len(x_tick_values) - 1) / (max_x_ticks - 1)) for i in range(max_x_ticks)]
            x_tick_values = [x_tick_values[i] for i in idx_list]
    else:
        step = (x_unique[-1] - x_unique[0]) / (max_x_ticks - 1)
        x_tick_values = [x_unique[0] + i * step for i in range(max_x_ticks)]
    if logy:
        floor = _log_floor(ys)
        p0 = math.floor(math.log10(max(min(ys), floor)))
        p1 = math.ceil(math.log10(max(max(ys), floor)))
        y_tick_values = [10.0 ** p for p in range(int(p0), int(p1) + 1)]
    else:
        y_tick_values = [rymin + (rymax - rymin) * i / 5.0 for i in range(6)]

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

    return (
        f'<svg width="{width}" height="{svg_h}" viewBox="0 0 {width} {svg_h}" xmlns="http://www.w3.org/2000/svg">\n'
        f'  <rect x="0" y="0" width="{width}" height="{svg_h}" fill="white"/>\n'
        f'  <text x="{width/2:.1f}" y="24" text-anchor="middle" font-size="16" font-family="Arial">{html_escape(title)}</text>\n'
        + " ".join(x_marks)
        + " ".join(y_marks)
        + f'  <line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="#333"/>\n'
        + f'  <line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="#333"/>\n'
        + " ".join(parts)
        + " ".join(legends)
        + f'  <text x="{width/2:.1f}" y="{height-18}" text-anchor="middle" font-size="13" font-family="Arial">{html_escape(x_label)}</text>\n'
        + f'  <text x="22" y="{height/2:.1f}" transform="rotate(-90 22,{height/2:.1f})" text-anchor="middle" font-size="13" font-family="Arial">{html_escape(y_label)}</text>\n'
        + "</svg>\n"
    )
