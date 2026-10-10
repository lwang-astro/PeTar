#!/usr/bin/env python3
"""Shared helpers for PeTar validation HTML reports.

Every report script imports its generic helpers from this module instead of
defining local copies; see test/validation/REPORT_DESIGN.md for the report
design rules (paper-style sections, figure encoding, caption conventions).
"""
import json
import math
import re
import shlex
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np  # used by orbital_elements_from_two_body_rows

Point = Tuple[float, float]
VERSION_RE = re.compile(r"^Version:\s*(.*)$")
G_MSUN_PC_MYR = 0.00449830997959438  # pc^3 Msun^-1 Myr^-2 (same as make_ic.py)


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



def metric_card(label: str, value: str, note: str = "") -> str:
    note_html = f'<div class="note">{html_escape(note)}</div>' if note else ""
    return (
        f'<div class="card"><div class="card-label">{html_escape(label)}</div>'
        f'<div class="card-value">{html_escape(value)}</div>{note_html}</div>'
    )


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


PALETTE = ["#0b74de", "#d81b60", "#2e7d32", "#8e24aa", "#f4511e", "#00897b"]


def _svg_frame(
    width: int,
    height: int,
    left: int,
    right: int,
    top: int,
    bottom: int,
    xs: List[float],
    ys: List[float],
):
    plot_w = width - left - right
    plot_h = height - top - bottom
    xmin, xmax = min(xs), max(xs)
    if xmax == xmin:
        xmax += 1.0
    floor = 1e-9
    ymin = math.log10(max(min(ys), floor))
    ymax = math.log10(max(max(ys), floor))
    if ymax == ymin:
        ymax += 1.0

    def sx(v: float) -> float:
        return left + (v - xmin) / (xmax - xmin) * plot_w

    def sy(v: float) -> float:
        return top + (ymax - math.log10(max(v, floor))) / (ymax - ymin) * plot_h

    return sx, sy, xmin, xmax, plot_w, plot_h


def svg_series_plot(
    series_map: Dict[str, List[Point]],
    title: str,
    x_label: str,
    y_label: str,
    hlines: List[Tuple[float, str]] = [],
    max_x_ticks: int = 8,
) -> str:
    """Log-y multi-series plot, auto palette, per-series legend on the plot
    (legacy style used by T1-T5)."""
    width, height = 940, 500
    left, right, top, bottom = 88, 24, 44, 70

    all_points = [p for points in series_map.values() for p in points]
    if not all_points:
        return "<p>No data available.</p>"
    xs = [p[0] for p in all_points]
    ys = [max(p[1], 1e-9) for p in all_points]
    sx, sy, _, _, _, _ = _svg_frame(width, height, left, right, top, bottom, xs, ys)

    parts: List[str] = []
    legends: List[str] = []
    for idx, (label, points) in enumerate(series_map.items()):
        color = PALETTE[idx % len(PALETTE)]
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
    p0 = math.floor(math.log10(max(min(ys), 1e-9)))
    p1 = math.ceil(math.log10(max(max(ys), 1e-9)))
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


def svg_curves_plot(
    curves: List[Tuple[str, List[Point], str, str]],
    legend_models: List[Tuple[str, str]],
    title: str,
    x_label: str,
    y_label: str,
    hlines: List[Tuple[float, str]] = [],
    bands: List[Tuple[str, List[Tuple[float, float]], str]] = [],
    max_x_ticks: int = 8,
) -> str:
    """Log-y plot where color encodes the model (run) and dash pattern the
    quantity (current style): curves are (label, points, color, dash) with
    dash == "" for solid; legend_models lists (model label, color). The legend
    is a dedicated strip below the plot whose height adapts to the entry
    count, so it can never overlap the curves. bands: (label, [(t0, t1)...],
    color) draw translucent full-height rectangles marking active intervals.

    See REPORT_DESIGN.md for when to use this versus svg_series_plot.
    """
    width, height = 940, 500
    left, right, top, bottom = 88, 24, 44, 70

    all_points = [p for _label, points, _c, _d in curves for p in points]
    if not all_points:
        return "<p>No data available.</p>"
    xs = [p[0] for p in all_points]
    ys = [max(p[1], 1e-9) for p in all_points]
    sx, sy, xmin, xmax, _, plot_h = _svg_frame(width, height, left, right, top, bottom, xs, ys)

    parts: List[str] = []
    for _label, band_list, band_color in bands:
        for t0, t1 in band_list:
            if t1 < xmin or t0 > xmax:
                continue
            x0, x1 = sx(max(t0, xmin)), sx(min(t1, xmax))
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
    n_tt_legend = 1 if bands else 0
    n_legend_rows = max(len(legend_models) + n_tt_legend, 2 + n_tt_legend)
    svg_h = height + 20 + (17 + 17 * n_legend_rows + 10)
    leg_y = height + 20
    col1_x = width - right - 500
    col2_x = width - right - 240
    legends: List[str] = [f'<text x="{col1_x}" y="{leg_y}" font-size="11" font-family="Arial" font-weight="bold">model (color)</text>']
    for i, (label, color) in enumerate(legend_models):
        ly = leg_y + 17 + i * 17
        legends.append(f'<line x1="{col1_x}" y1="{ly}" x2="{col1_x+26}" y2="{ly}" stroke="{color}" stroke-width="3"/>')
        legends.append(f'<text x="{col1_x+31}" y="{ly+4}" font-size="11" font-family="Arial">{html_escape(label)}</text>')
    legends.append(f'<text x="{col2_x}" y="{leg_y}" font-size="11" font-family="Arial" font-weight="bold">quantity (line)</text>')
    for i, (label, dash) in enumerate([("instantaneous", ""), ("cumulative", "7,4")]):
        ly = leg_y + 17 + i * 17
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        legends.append(f'<line x1="{col2_x}" y1="{ly}" x2="{col2_x+26}" y2="{ly}" stroke="#555" stroke-width="2"{dash_attr}/>')
        legends.append(f'<text x="{col2_x+31}" y="{ly+4}" font-size="11" font-family="Arial">{html_escape(label)}</text>')
    if bands:
        ly = leg_y + 17 + 2 * 17 + 2
        legends.append(f'<rect x="{col2_x}" y="{ly-8}" width="26" height="11" fill="#555" fill-opacity="0.18"/>')
        legends.append(f'<text x="{col2_x+31}" y="{ly+2}" font-size="11" font-family="Arial">active intervals (run color)</text>')

    x_unique = sorted(set(xs))
    if len(x_unique) <= max_x_ticks:
        x_tick_values = x_unique
    else:
        step = (x_unique[-1] - x_unique[0]) / (max_x_ticks - 1)
        x_tick_values = [x_unique[0] + i * step for i in range(max_x_ticks)]
    p0 = math.floor(math.log10(max(min(ys), 1e-9)))
    p1 = math.ceil(math.log10(max(max(ys), 1e-9)))
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
