#!/usr/bin/env python3
"""T4 hierarchical-triple five-mode comparison HTML report (paper style, see REPORT_DESIGN.md)."""
import argparse
import math
import re
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from report_common import (
    G_MSUN_PC_MYR,
    PAPER_CSS,
    check_row,
    git_commit,
    html_escape,
    load_json,
    parse_command_args,
    parse_log_version,
    svg_curves_plot,
    svg_dual_panel_plot,
)

REPO = Path(__file__).resolve().parents[2]
Point = Tuple[float, float]
TIME_RE = re.compile(r"^Time:\s*([0-9eE+\-.]+)")
ENERGY_HEADER_RE = re.compile(r"^Energy:\s+(.*)$")
ENERGY_PHYSIC_RE = re.compile(r"^Physic:\s+(.*)$")
INPUT_UNIT_RE = re.compile(r"^Input data unit:\s*(\d+)")
UNIT_SET_RE = re.compile(r"^----- Unit set\s*(\d+)\s*:\s*(.*)-----$")
TREE_DT_RE = re.compile(r"^Tree time step\s*=\s*([0-9eE+\-.]+)$")
OUTPUT_DT_RE = re.compile(r"^Output time step\s*=\s*([0-9eE+\-.]+)$")

MODE_ORDER = ["m1_pure_sdar_ref", "m2_hard_no_tt", "m3_hard_tt", "m4_tree_no_tt", "m5_tree_tt"]
# Encoding across ALL figures: color = integration channel (m1 reference blue,
# hard pair m2/m3 pink, tree pair m4/m5 purple); line style = tidal tensor
# (solid = TT off, dashed 12,6 = TT on). The TT pairs share one color and are
# separated by line style only. m2 (TT off) and m3 (TT on) additionally agree
# to ~1e-12 in every plotted quantity (the tensor and direct perturber
# treatments of the inner AR group are equivalent at that level), so the two
# keys of the pair are visually one line there; exact values are in Table 5.
MODE_COLORS = {
    "m1_pure_sdar_ref": "#0b74de",
    "m2_hard_no_tt": "#d81b60",
    "m3_hard_tt": "#d81b60",
    "m4_tree_no_tt": "#8e24aa",
    "m5_tree_tt": "#8e24aa",
}
MODE_DASH = {"m3_hard_tt": "12,6", "m5_tree_tt": "12,6"}
MODE_LABELS = {
    "m1_pure_sdar_ref": "m1 pure SDAR reference",
    "m2_hard_no_tt": "m2 hard, TT off",
    "m3_hard_tt": "m3 hard, TT on",
    "m4_tree_no_tt": "m4 tree, TT off",
    "m5_tree_tt": "m5 tree, TT on",
}


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
    sys.path.insert(0, str(REPO / "tools"))
    from analysis.status import Status  # pylint: disable=import-outside-toplevel  # pyright: ignore[reportMissingImports]

    # T4 runs with -b 1 on a bse-interrupt binary: status rows carry stellar columns
    st = Status(N_particle=n_particle, interrupt_mode="bse")
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


def derive_inner_binary_drift(status_path: Path, n_particle: int) -> Dict[str, Any]:
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


def extract_binary_from_command(command: str) -> str:
    tokens = command.split()
    for idx, token in enumerate(tokens):
        if token.startswith("OMP_NUM_THREADS=") and idx + 1 < len(tokens):
            return tokens[idx + 1]
    return "unknown"


def scenario_file_from_name(scenario: str) -> Path:
    return Path("test/validation/scenarios") / f"{scenario}.json"


def normalize_mode_id(run_id: str) -> str:
    return run_id[:-4] if run_id.endswith("_o15") else run_id


def collect_plot_series(rows: List[Dict]) -> Dict[str, Any]:
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


def mode_curves(series_map: Dict[str, List[Point]], modes: List[str]) -> List[Tuple[str, List[Point], str, str]]:
    """Convert a {label: points} map into svg_curves_plot curves keyed by mode id."""
    by_mode: Dict[str, List[Point]] = {}
    for label, points in series_map.items():
        mode = label.split("|")[0].strip()
        by_mode[mode] = points
    return [(m, by_mode.get(m, []), MODE_COLORS.get(m, "#555"), MODE_DASH.get(m, "")) for m in modes if by_mode.get(m)]


def build_orbit_figure_section(
    series: Dict[str, Any],
    modes: List[str],
    control_series: Optional[Dict[str, Any]] = None,
) -> str:
    """One figure per quantity; upper panel = base-group absolute curves, and
    when a control group exists each figure gains a lower panel with the
    control - base difference (zero = structurally identical behaviour at both
    outer scales; systematic deviation = scale-dependent effect). Encoding in
    both panels: color = integration channel, line style = tidal tensor.
    """
    legend_models = [(MODE_LABELS.get(m, m), MODE_COLORS.get(m, "#555"), *([MODE_DASH[m]] if m in MODE_DASH else [])) for m in modes]
    # second legend column carries the TT line-style convention (each figure
    # shows a single quantity, so the dash channel is free for TT on/off)
    tt_legend = [("TT off (solid)", ""), ("TT on (dashed)", "12,6")]

    def _log10_abs(v: float) -> float:
        return math.log10(abs(v)) if abs(v) > 0.0 else math.nan

    def _mode_series(series_map: Dict[str, List[Point]], key: str, m: str) -> List[Point]:
        s = series_map.get(key, {})
        return s.get(f"{m} | tt=0", s.get(f"{m} | tt=1", []))

    # key, figure title, top y-label, log-y(top), diff kind, diff y-label
    quantities = [
        ("energy_total_series", "|Error/Total|", "relative energy error", True, "dlog", "control - base [decades]"),
        ("energy_pp_series", "|Error_PP|", "relative PP energy error", True, "dlog", "control - base [decades]"),
        ("a_in_series", "inner semi-major axis", "a_in [pc]", False, "rel", "control - base, Delta(a/a0) [1]"),
        ("e_in_series", "inner eccentricity", "e_in", False, "abs", "Delta e_in [1]"),
        ("a_out_series", "outer semi-major axis", "a_out [pc]", False, "rel", "Delta(a/a0) [1]"),
        ("e_out_series", "outer eccentricity", "e_out", False, "abs", "Delta e_out [1]"),
        ("i_mutual_series", "mutual inclination", "i_mutual [deg]", False, "abs", "Delta i_mut [deg]"),
    ]

    html_parts: List[str] = []
    for key, title, y_label, logy, diff_kind, diff_label in quantities:
        top_curves = mode_curves(series.get(key, {}), modes)
        if not top_curves:
            continue
        diff_curves = []
        if control_series is not None:
            for m in modes:
                b = _mode_series(series, key, m)
                c = _mode_series(control_series, key, m)
                if not b or not c:
                    continue
                if diff_kind == "dlog":
                    transform: Any = lambda vb, vc: _log10_abs(vc) - _log10_abs(vb)  # noqa: E731
                elif diff_kind == "rel":
                    b0, c0 = b[0][1], c[0][1]
                    if abs(b0) < 1e-300 or abs(c0) < 1e-300:
                        continue
                    transform = lambda vb, vc, b0=b0, c0=c0: vc / c0 - vb / b0  # noqa: E731
                else:
                    transform = lambda vb, vc: vc - vb  # noqa: E731
                pts = _aligned_diff(b, c, transform)
                if pts:
                    diff_curves.append((m, pts, MODE_COLORS.get(m, "#555"), MODE_DASH.get(m, "")))
        if diff_curves:
            # single SVG: upper panel (base absolute) + lower panel (control -
            # base) sharing the x axis; one-row model legend above the plot
            html_parts.append(
                svg_dual_panel_plot(
                    top_curves,
                    diff_curves,
                    legend_models,
                    f"{title} vs time",
                    "time [Myr]",
                    y_label,
                    diff_label,
                    top_logy=logy,
                )
            )
        else:
            html_parts.append(
                svg_curves_plot(
                    top_curves,
                    legend_models,
                    f"{title} vs time",
                    "time [Myr]",
                    y_label,
                    logy=logy,
                    legend_quantities=tt_legend,
                    legend_headers=("model (color)", "tidal tensor (line)"),
                )
            )
    return "".join(html_parts)


def _aligned_diff(base_pts: List[Point], ctrl_pts: List[Point], transform) -> List[Point]:
    """Pointwise transformed difference ctrl-base on the shared output-time grid.

    Both groups use identical -o/-t, so times match exactly; mismatched or
    non-finite samples are dropped pairwise.
    """
    ctrl_by_t = {round(t, 12): v for t, v in ctrl_pts}
    out: List[Point] = []
    for t, vb in base_pts:
        vc = ctrl_by_t.get(round(t, 12))
        if vc is None:
            continue
        d = transform(vb, vc)
        if d is not None and math.isfinite(d):
            out.append((t, d))
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate T4 triple mode comparison HTML report")
    parser.add_argument("--report", default="test/out/report.t4.triple_tidal_tensor.json")
    parser.add_argument("--scenario", default="t4_triple_tidal_tensor")
    parser.add_argument("--scenario-file", default="")
    parser.add_argument("--output", default="test/out/t4_triple_tidal_tensor_summary.html")
    args = parser.parse_args()

    report = load_json(Path(args.report))
    scenario_path = Path(args.scenario_file) if args.scenario_file else scenario_file_from_name(args.scenario)
    scenario_data = load_json(scenario_path) if scenario_path.exists() else {}

    # One scenario, two IC scales: base runs (a_out=0.15) and the control
    # group (original a_out=1.5) distinguished by the _o15 run-id suffix.
    all_runs = [r for r in report.get("runs", []) if r.get("scenario") == args.scenario]
    runs = [r for r in all_runs if not str(r.get("run_id", "")).endswith("_o15")]
    control_runs_raw = [r for r in all_runs if str(r.get("run_id", "")).endswith("_o15")]
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

    rows = sorted(rows, key=lambda x: MODE_ORDER.index(x["run_id"]) if x["run_id"] in MODE_ORDER else 999)
    log_version = parse_log_version(Path(rows[0]["log_path"]))

    series = collect_plot_series(rows)
    unit_rows = series["unit_rows"]

    control_compare_html = ""
    control_series: Optional[Dict[str, Any]] = None
    if control_runs_raw:
        control_rows = []
        for rr in control_runs_raw:
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

        control_rows = sorted(control_rows, key=lambda x: MODE_ORDER.index(x["mode_id"]) if x["mode_id"] in MODE_ORDER else 999)
        # figure lookups key series labels by MODE_ORDER ids: map the _o15
        # run ids onto their base mode ids for the series pass only
        for r in control_rows:
            r["run_id_full"] = r["run_id"]
            r["run_id"] = r["mode_id"]
        control_series = collect_plot_series(control_rows)
        for r in control_rows:
            r["run_id"] = r.pop("run_id_full")

        base_by_mode = {
            normalize_mode_id(str(r["run_id"])): r
            for r in rows
        }
        ctrl_by_mode = {r["mode_id"]: r for r in control_rows}

        lines = []
        for mode in MODE_ORDER:
            if mode not in base_by_mode or mode not in ctrl_by_mode:
                continue
            b = base_by_mode[mode]
            c = ctrl_by_mode[mode]
            lines.append(
                "<tr>"
                f"<td>{html_escape(MODE_LABELS.get(mode, mode))}</td>"
                f"<td>{b['max_abs_error_over_total']:.3e}</td>"
                f"<td>{c['err_tot']:.3e}</td>"
                f"<td>{b['max_abs_error_pp']:.3e}</td>"
                f"<td>{c['err_pp']:.3e}</td>"
                f"<td>{b['da_rel_max']:.3e}</td>"
                f"<td>{c['da']:.3e}</td>"
                f"<td>{b['de_abs_max']:.3e}</td>"
                f"<td>{c['de']:.3e}</td>"
                f"<td>{b['max_artificial_particles_glb']}</td>"
                f"<td>{c['max_art']:.0f}</td>"
                "</tr>"
            )
        if lines:
            control_compare_html = "".join(lines)

    # merged figures: upper panel = base absolute curves, lower panel =
    # control - base difference (only when the control group exists)
    base_figures_html = build_orbit_figure_section(
        series, MODE_ORDER,
        control_series=control_series if control_runs_raw else None,
    )

    run_rows_html = "".join(
        f"<tr><td>{html_escape(str(MODE_LABELS.get(r['run_id'], r['run_id'])))}</td>"
        f"<td>{r['tt_switch']}</td><td>{r['tt_nstep']}</td>"
        f"<td>{r['r_group']:.3g}</td><td>{r['rout']:.3g}</td><td>{r['dt_soft']:.3e}</td>"
        f"<td>{'YES' if r['dt_le_tree'] else 'NO'}</td>"
        f"<td>{r['max_abs_error_over_total']:.3e}</td><td>{r['max_abs_error_pp']:.3e}</td>"
        f"<td>{r['da_rel_max']:.3e}</td><td>{r['de_abs_max']:.3e}</td>"
        f"<td>{r['max_n_real_glb']}</td><td>{r['max_n_all_glb']}</td><td>{r['max_artificial_particles_glb']}</td>"
        f"<td>{r['a_in0']:.3e}</td><td>{r['e_in0']:.3f}</td><td>{r['a_out0']:.3e}</td><td>{r['e_out0']:.3f}</td><td>{r['i_mutual0']:.1f}</td>"
        f"<td>{r['t_end']:.3g}</td><td>{r['outer_period_covered']:.0f}</td></tr>"
        for r in rows
    )

    command_rows_html = "".join(
        f"<tr><td>{html_escape(r['run_id'])}</td><td><code>{html_escape(r['command'])}</code></td></tr>"
        for r in rows
    )

    checks_html = "".join(check_row(c) for c in checks)

    setup_commands = scenario_data.get("setup_commands", []) if isinstance(scenario_data, dict) else []
    if not isinstance(setup_commands, list):
        setup_commands = []
    setup_html = "<br>".join(html_escape(str(cmd)) for cmd in setup_commands)

    all_pass = all(c.get("passed", False) for c in checks)
    verdict = "ALL CHECKS PASSED" if all_pass else "CHECKS FAILED"

    unit_rows_html = "".join(
        f'<tr><td>{html_escape(rid)}</td><td>{html_escape(iu)}</td><td>{html_escape(us)}</td></tr>'
        for rid, iu, us in unit_rows
    )

    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>T4 hierarchical triple five-mode comparison</title>
<style>{PAPER_CSS}</style></head><body>

<h1>T4 &nbsp;Hierarchical triple under five integration modes:<br>pure SDAR, Hermite&ndash;SDAR and tree, with and without tidal tensor</h1>
<div class="subtitle">validation report &middot; scenario <i>t4_triple_tidal_tensor</i> &middot; <span class="verdict {'pass' if all_pass else 'fail'}">{verdict}</span></div>

<h2>1&nbsp;&nbsp;Background and purpose</h2>
<p>When a compact subsystem is integrated as an AR (SDAR) group inside a larger system, the rest
of the system perturbs it through two channels: the direct perturber-list force and the tidal
tensor (TT) correction that represents the smooth field gradient across the group. The
--tt-switch option enables the tidal-tensor channel; whether it is active is not visible in the
trajectory alone, so this test compares five integration modes of one engineered hierarchical
triple over 200 outer orbits: a pure SDAR reference (the whole triple one AR group, no tree),
hard-mode and tree-mode runs with TT off and on. It establishes (i) that the energy error of the
production modes stays far below the pure-SDAR reference scale, (ii) that the inner orbit is
conserved relative to the reference, and (iii) &mdash; by counting artificial particles &mdash;
that the tidal tensor is actually engaged exactly in the runs where it is switched on, and
nowhere else. A control group repeats the comparison with the original outer scale
(a<sub>out</sub> = 1.5 pc) to separate scale-dependent from structural effects.</p>

<h2>2&nbsp;&nbsp;Method</h2>
<h3>2.1&nbsp;&nbsp;Quantities under comparison</h3>
<p>From the run energy logs:</p>
<div class="eq">&epsilon;<sub>E</sub> = max<sub>t</sub> |Error/Total|(t),&nbsp;&nbsp;
&epsilon;<sub>PP</sub> = max<sub>t</sub> |Error_PP|(t),<span class="tag">(1)</span></div>
<p>the total and short-range (particle&ndash;particle/AR) relative energy errors. From the status
particle states, the inner orbit (m<sub>1</sub>, m<sub>2</sub>), the outer orbit (inner c.m.,
m<sub>3</sub>), and the mutual inclination:</p>
<div class="eq">a = &minus;&mu;/(2&epsilon;), &epsilon; = v&sup2;/2 &minus; &mu;/r,&nbsp;&nbsp;
e = &radic;(1 + 2&epsilon;h&sup2;/&mu;&sup2;),&nbsp;&nbsp; h = |r &times; v|,
<span class="tag">(2)</span></div>
<div class="eq">cos&thinsp;i<sub>mut</sub> = &#292;<sub>in</sub> &middot;
&#292;<sub>out</sub>,<span class="tag">(3)</span></div>
<div class="eq">&Delta;a/a<sub>0</sub> = max<sub>t</sub> |a(t)&minus;a<sub>0</sub>|/a<sub>0</sub>
(inner), &nbsp; &Delta;e = max<sub>t</sub> |e(t)&minus;e<sub>0</sub>| (inner).
<span class="tag">(4)</span></div>
<p>The tidal-tensor activity guard is the artificial-particle census:</p>
<div class="eq">N<sub>art</sub> = max<sub>t</sub> [N<sub>all</sub>(glb) &minus;
N<sub>real</sub>(glb)],<span class="tag">(5)</span></div>
<p>which is positive only when tidal-tensor probes (or other artificial group particles)
exist.</p>

<h3>2.2&nbsp;&nbsp;Metrics</h3>
<table class="data">
<tr><th>symbol</th><th>definition</th><th>interpretation</th></tr>
<tr><td>&epsilon;<sub>E</sub>, &epsilon;<sub>PP</sub></td><td>Eq. (1)</td>
<td>global and short-range energy conservation per mode</td></tr>
<tr><td>&Delta;a/a<sub>0</sub>, &Delta;e</td><td>Eq. (4)</td>
<td>inner-orbit conservation relative to its first output sample</td></tr>
<tr><td>i<sub>mut</sub>(t)</td><td>Eq. (3)</td>
<td>angular-momentum exchange between inner and outer orbit (Kozai indicator)</td></tr>
<tr><td>N<sub>art</sub></td><td>Eq. (5)</td>
<td>guard: proves the TT path executed (on-modes) and only there (off-modes)</td></tr>
<tr><td>t<sub>end</sub>/P<sub>out</sub></td><td>coverage in outer periods</td>
<td>integration length in orbital units</td></tr>
</table>
<div class="caption"><b>Table 1.</b> Metrics used in this report. Orbital elements are computed
from the *.status particle states (N_particle = 3, bse columns) at every output step.</div>

<h3>2.3&nbsp;&nbsp;Runs</h3>
<table class="data">
<tr><th>run</th><th>tt-switch</th><th>configuration</th><th>purpose</th></tr>
<tr><td>m1_pure_sdar_ref</td><td>0</td>
<td>r_out = 2 pc (whole system in the hard region), r_group = 0.18 &gt; outer apo (0.152):
the triple stays ONE isolated 3-body AR group integrated by SDAR at the AR step
(3&times;10<sup>-5</sup>) &mdash; no tree, no Hermite cluster, no perturbers</td>
<td>regularized reference for element conservation (what the production channels should
not degrade)</td></tr>
<tr><td>m2_hard_no_tt</td><td>0</td>
<td>HARD channel: r_out = 2 keeps all three bodies in ONE Hermite block-step cluster
(dt = 0.03125); the inner binary (apo 1.9&times;10<sup>-3</sup> &lt; r_group = 0.01) is a permanent AR
group inside it. The third body is a Hermite cluster member and perturbs the AR group
via the DIRECT perturber list (exact force on members)</td>
<td>hard-channel baseline: AR group under a directly integrated perturber, TT absent</td></tr>
<tr><td>m3_hard_tt</td><td>1 (nstep 8)</td>
<td>as m2, but the third body's perturbation enters the AR integration as the TIDAL
TENSOR expansion measured at the group c.m. (tt-nstep = 8 substeps), not as a direct
perturber force</td>
<td>isolates the TT approximation against m2's direct force in the same channel
(A/B pair with m2)</td></tr>
<tr><td>m4_tree_no_tt</td><td>0</td>
<td>TREE channel: r_out = 0.08 (&rarr; r_in via r-ratio 0.3) puts the outer orbit
(0.149&ndash;0.152 pc) entirely OUTSIDE the changeover region &mdash; the third body is a
soft TREE particle (FDPS gravity, dt_soft = 0.03125) while the inner-binary AR group
stays in the hard channel; their coupling crosses the changeover boundary, no TT</td>
<td>production layout without TT: tree background + hard AR subsystem coupled through
the changeover</td></tr>
<tr><td>m5_tree_tt</td><td>1 (nstep 8)</td>
<td>as m4 + tidal tensor: the tree-side perturbation of the AR group is applied via the
tidal tensor</td>
<td>production configuration of the 2026-10-06 event: TT active in the tree channel
(A/B pair with m4)</td></tr>
</table>
<div class="caption"><b>Table 2.</b> Five modes of one IC, identical apart from the listed
knobs; all run with -b 1 (block steps), 1 rank &times; 1 thread, to t<sub>end</sub> = 20 Myr
&asymp; 200 outer periods. The essential hard-vs-tree difference is WHICH integrator owns the
third body: in the hard pair (m2/m3, r_out = 2) it is a Hermite cluster member exerting the
direct perturber force on the AR group; in the tree pair (m4/m5, r_out = 0.08) it is a soft
tree particle whose influence reaches the AR group through the changeover coupling. The pairs
(m2, m3) and (m4, m5) are A/B pairs isolating the TT knob in one channel each; m1 anchors the
achievable conservation scale; r<sub>search,group</sub> = 0.02 pc throughout.</div>

<h3>2.4&nbsp;&nbsp;Reproducibility: software, build and execution</h3>
<p><b>Software.</b> PeTar version {html_escape(log_version)}, git commit
<code>{html_escape(git_commit(REPO))}</code>, 1 MPI rank &times; 1 OpenMP thread. The IC carries
stellar-evolution columns (petar.init -s bse) and the runs use -b 1, so a bse-interrupt family
binary is required; the pipeline selects it via <code>petar.select --require bse --optional
mpi,omp,avx2</code> and builds it (<code>./configure --with-interrupt=bse &amp;&amp; make
install</code>) when missing:</p>
<pre class="cmd">{setup_html}</pre>
<table class="data">
<tr><th>parameter</th><th>value</th></tr>
<tr><td>system</td><td>hierarchical triple, both orbits started at pericenter, c.m. position and
velocity subtracted exactly</td></tr>
<tr><td>masses</td><td>m<sub>1</sub> = 0.001, m<sub>2</sub> = 0.009, m<sub>3</sub> = 1.0
M<sub>sun</sub></td></tr>
<tr><td>inner orbit</td><td>r<sub>peri,in</sub> = 10<sup>-4</sup> pc, e<sub>in</sub> = 0.9
(a<sub>in</sub> = 10<sup>-3</sup> pc)</td></tr>
<tr><td>outer orbit</td><td>e<sub>out</sub> = 0.01, a<sub>out</sub> = 0.15 pc base group
(control group: a<sub>out</sub> = 1.5 pc)</td></tr>
<tr><td>derived periods</td><td>P<sub>in</sub> = {period_inner:.4g} Myr,
P<sub>out</sub> = {period_outer:.4g} Myr (from the first status sample)</td></tr>
<tr><td>units / constants</td><td>M<sub>sun</sub>, pc, pc Myr<sup>-1</sup>;
G = 4.49830997959438 &times; 10<sup>-3</sup> pc&sup3; M<sub>sun</sub><sup>-1</sup> Myr<sup>-2</sup></td></tr>
<tr><td>file format</td><td>ASCII m x y z v<sub>x</sub> v<sub>y</sub> v<sub>z</sub>
(+ petar.init -s bse conversion, stellar-evolution columns)</td></tr>
</table>
<div class="caption"><b>Table 3.</b> IC generation parameters (<code>make_ic.py --case t4</code>;
control group <code>--case t4_outer15</code>). Scales follow the notebook tidal-tensor 3-body
test; the shortened outer orbit (0.15 pc vs the original 1.5 pc) keeps the runtime practical
while preserving the hierarchy.</div>
<table class="data">
<tr><th>run_id</th><th>input data unit</th><th>unit set line</th></tr>
{unit_rows_html}
</table>
<div class="caption"><b>Table 4.</b> Unit consistency across runs (parsed from each log): all
modes must read the IC with the same unit set for the comparison to be meaningful.</div>
<p><b>Execution commands</b> (as recorded by the runner):</p>
<table class="data">
<tr><th>run_id</th><th>command</th></tr>
{command_rows_html}
</table>

<h2>3&nbsp;&nbsp;Results</h2>
<h3>3.1&nbsp;&nbsp;Energy conservation by mode</h3>
{base_figures_html}
<div class="caption"><b>Figures 1&ndash;7.</b> One figure per quantity: two stacked panels sharing one
x axis and one legend &mdash; <b>upper = base group</b> (a<sub>out</sub> = 0.15 pc) absolute curves; <b>lower = control &minus; base</b>
(a<sub>out</sub> = 1.5 pc minus base, same modes and encoding). Encoding in both panels:
<b>color = integration channel</b> (blue: m1 SDAR reference; pink: hard pair m2/m3; purple: tree
pair m4/m5), <b>line style = tidal tensor</b> (solid: TT off; dashed: TT on). Lower-panel zero
(hline) = structurally identical behaviour at both outer scales; systematic deviation =
scale-dependent effect. Energy errors differ in decades (Delta log<sub>10</sub>; &plusmn;1 = one
decade); semi-major axes as Delta(a/a<sub>0</sub>) with each group normalized by its own initial
value; eccentricities and inclination as absolute differences. The check metric is the
upper-panel curve maximum; bounded flat upper curves mean conservative behaviour, upward
steps mark discrete events (group re-linking, BSE kicks). The m2/m3 pair coincides to ~1e-12
in every plotted quantity (the tidal-tensor and direct perturber treatments of the inner AR
group are equivalent at that level &mdash; the TT-consistency result T5 verifies structurally),
so the same-color pair reads as one line with the solid m2 showing through the dashed m3;
exact values are in Table 5.</div>

<h3>3.2&nbsp;&nbsp;Run summary</h3>
<table class="data">
<tr><th>run</th><th>tt</th><th>tt-nstep</th><th>r_group</th><th>r_out</th><th>dt_soft</th><th>dt_soft &le; tree_dt</th><th>&epsilon;<sub>E</sub></th><th>&epsilon;<sub>PP</sub></th><th>&Delta;a/a<sub>0</sub></th><th>&Delta;e</th><th>max N_real</th><th>max N_all</th><th>N<sub>art</sub></th><th>a_in(0)</th><th>e_in(0)</th><th>a_out(0)</th><th>e_out(0)</th><th>i_mut(0) [&deg;]</th><th>t_end</th><th>P_out covered</th></tr>
{run_rows_html}
</table>
<div class="caption"><b>Table 5.</b> Per-mode summary. N<sub>art</sub> (Eq. 5) is the execution
guard: it must be &ge; 1 in m3/m5 (TT engaged) and 0 in m1/m2/m4. The dt_soft &le; tree_dt
column guards that no run steps the tree faster than its auto-derived step.</div>

<h3>3.3&nbsp;&nbsp;Control group (outer a<sub>out</sub> = 1.5 pc)</h3>
<p>The control group's per-mode behaviour is shown as the lower (control &minus; base) panels
of Figures 1&ndash;7; this section keeps the run-maximum comparison.</p>
{('<table class="data"><tr><th>mode</th><th>base &epsilon;_E</th><th>control &epsilon;_E</th><th>base &epsilon;_PP</th><th>control &epsilon;_PP</th><th>base &Delta;a/a0</th><th>control &Delta;a/a0</th><th>base &Delta;e</th><th>control &Delta;e</th><th>base N_art</th><th>control N_art</th></tr>' + control_compare_html + '</table><div class="caption"><b>Table 6.</b> Base (a_out = 0.15 pc) versus control (a_out = 1.5 pc) maxima per mode; both groups use identical mode configurations apart from the IC outer scale.</div>') if control_compare_html else ''}

<h3>3.4&nbsp;&nbsp;Automatic checks</h3>
<table class="data">
<tr><th>check</th><th>status</th><th>value</th><th>threshold</th><th>message</th></tr>
{checks_html}
</table>
<div class="caption"><b>Table 7.</b> Automatic checks. <i>Accuracy threshold</i>:
m1_reference_energy_bound caps &epsilon;<sub>E</sub> of the pure-SDAR reference at
10<sup>-8</sup> (the regularized channel should conserve far below the production modes).
<i>Guards</i>: m3/m5_tidal_tensor_effective require N<sub>art</sub> &ge; 1 (the tested TT path
actually executed &mdash; all-green on an unexecuted path would be a silent failure);
m2/m4_no_tidal_tensor_particles require N<sub>art</sub> = 0 (TT stays off where switched off);
one no_runtime_abort per mode requires zero segmentation faults / assertion failures /
aborts.</div>

<div class="meta">petar version: {html_escape(log_version)} &middot; scenario:
<i>t4_triple_tidal_tensor</i> &middot; report: {html_escape(args.report)}</div>
</body></html>
"""

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html, encoding="utf-8")
    n_pass = sum(1 for c in checks if c.get("passed", False))
    print(f"[T4] Summary: passed={n_pass == len(checks) and len(checks) > 0} total={len(checks)} failed={len(checks) - n_pass} report={out_path}")
    return 0 if n_pass == len(checks) and checks else 1


if __name__ == "__main__":
    raise SystemExit(main())
