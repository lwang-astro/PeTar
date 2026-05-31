#!/usr/bin/env python3
import argparse
import math
import os
import re
import subprocess
import tempfile
from pathlib import Path
from typing import List, Tuple


def load_astro_units_constants() -> Tuple[float, float]:
    """Load G_ASTRO and KMS_TO_PCMYR from PeTar src/astro_units.hpp.

    Fallback defaults are used by caller if header is missing or malformed.
    """
    header = Path(__file__).resolve().parents[2] / "src" / "astro_units.hpp"
    text = header.read_text(encoding="utf-8")

    pattern = re.compile(r"^\s*#define\s+(\w+)\s+([+\-]?(?:\d+\.?\d*|\d*\.\d+)(?:[eE][+\-]?\d+)?)")
    constants = {}
    for line in text.splitlines():
        m = pattern.match(line)
        if m is not None:
            constants[m.group(1)] = float(m.group(2))

    return constants["G_ASTRO"], constants["KMS_TO_PCMYR"]


try:
    G_MSUN_PC_MYR, KM_S_TO_PC_MYR = load_astro_units_constants()
except Exception:
    # Conservative fallback in case repository layout changes.
    G_MSUN_PC_MYR = 0.00449830997959438
    KM_S_TO_PC_MYR = 1.022712165045695

PC_MYR_TO_KM_S = 1.0 / KM_S_TO_PC_MYR
MCLUSTER_SEED = 42
MCLUSTER_OUTPUT_VELOCITY_UNIT = "km/s"
FUNCTIONAL_RECENTER_AFTER_INJECTION = os.environ.get("FUNCTIONAL_RECENTER_AFTER_INJECTION", "0").strip().lower() in {"1", "true", "yes", "on"}

# Fixed binary parameters for functional_mcluster_dual_merge_smoke.
# Format per entry: (m1, m2, period_myr, ecc)
FUNCTIONAL_FIXED_BINARIES: List[Tuple[float, float, float, float]] = [
    (43.95648895255284, 59.062621526994754, 2.1688326622265456e-08, 0.7090993394052765),
    (120.80378755706538, 124.43730278955712, 9.475226105745696e-09, 0.02203162007734072),
    (19.069022579079945, 23.77114358019302, 2.0531407662661103e-08, 0.5909878636500059),
]


def env_float(name: str, default: float) -> float:
    value = os.environ.get(name)
    return float(value) if value is not None else default


def resolve_functional_binary_params() -> List[Tuple[float, float, float, float]]:
    """Resolve binary parameters for functional IC generation.

    The default is fully script-fixed values to avoid external file dependency.
    Optional env overrides are provided for controlled tuning experiments.
    """
    b1 = FUNCTIONAL_FIXED_BINARIES[0]
    b2 = FUNCTIONAL_FIXED_BINARIES[1]

    resolved = [
        (
            env_float("FUNCTIONAL_B1_M1", b1[0]),
            env_float("FUNCTIONAL_B1_M2", b1[1]),
            env_float("FUNCTIONAL_B1_PERIOD_MYR", b1[2]),
            env_float("FUNCTIONAL_B1_ECC", b1[3]),
        ),
        (
            env_float("FUNCTIONAL_B2_M1", b2[0]),
            env_float("FUNCTIONAL_B2_M2", b2[1]),
            env_float("FUNCTIONAL_B2_PERIOD_MYR", b2[2]),
            env_float("FUNCTIONAL_B2_ECC", b2[3]),
        ),
        (
            env_float("FUNCTIONAL_B3_M1", FUNCTIONAL_FIXED_BINARIES[2][0]),
            env_float("FUNCTIONAL_B3_M2", FUNCTIONAL_FIXED_BINARIES[2][1]),
            env_float("FUNCTIONAL_B3_PERIOD_MYR", FUNCTIONAL_FIXED_BINARIES[2][2]),
            env_float("FUNCTIONAL_B3_ECC", FUNCTIONAL_FIXED_BINARIES[2][3]),
        ),
    ]
    return resolved


def two_body_peri_state(m1: float, m2: float, peri: float, ecc: float, g_const: float) -> Tuple[List[float], List[float], List[float], List[float]]:
    mass_sum = m1 + m2
    rel_r = [peri, 0.0, 0.0]
    rel_vy = math.sqrt(g_const * mass_sum * (1.0 + ecc) / peri)
    rel_v = [0.0, rel_vy, 0.0]

    factor1 = -m2 / mass_sum
    factor2 = m1 / mass_sum

    r1 = [factor1 * x for x in rel_r]
    r2 = [factor2 * x for x in rel_r]
    v1 = [factor1 * x for x in rel_v]
    v2 = [factor2 * x for x in rel_v]
    return r1, v1, r2, v2


def two_body_apo_state(m1: float, m2: float, semi: float, ecc: float, g_const: float) -> Tuple[List[float], List[float], List[float], List[float]]:
    """Construct a two-body state at apo-center (true anomaly = pi)."""
    mass_sum = m1 + m2
    apo = semi * (1.0 + ecc)
    rel_r = [-apo, 0.0, 0.0]
    rel_vy = math.sqrt(g_const * mass_sum * (1.0 - ecc) / apo)
    rel_v = [0.0, -rel_vy, 0.0]

    factor1 = -m2 / mass_sum
    factor2 = m1 / mass_sum

    r1 = [factor1 * x for x in rel_r]
    r2 = [factor2 * x for x in rel_r]
    v1 = [factor1 * x for x in rel_v]
    v2 = [factor2 * x for x in rel_v]
    return r1, v1, r2, v2


def two_body_hyperbolic_inbound_state(
    m1: float,
    m2: float,
    semi: float,
    ecc: float,
    true_anomaly: float,
    g_const: float,
) -> Tuple[List[float], List[float], List[float], List[float]]:
    """Construct a hyperbolic two-body state with inbound relative motion.

    semi must be negative and ecc > 1. true_anomaly should be in
    (-arccos(-1/ecc), 0) for inbound motion before peri-center.
    """
    if semi >= 0.0:
        raise ValueError(f"Hyperbolic orbit requires negative semi-major axis, got {semi}")
    if ecc <= 1.0:
        raise ValueError(f"Hyperbolic orbit requires ecc>1, got {ecc}")

    f_inf = math.acos(-1.0 / ecc)
    if not (-f_inf < true_anomaly < 0.0):
        raise ValueError(
            f"Inbound true anomaly must be in (-acos(-1/e), 0), got {true_anomaly} for ecc={ecc}"
        )

    mu = g_const * (m1 + m2)
    p = semi * (1.0 - ecc * ecc)
    r = p / (1.0 + ecc * math.cos(true_anomaly))
    h = math.sqrt(mu * p)

    vr = mu / h * ecc * math.sin(true_anomaly)
    vt = mu / h * (1.0 + ecc * math.cos(true_anomaly))

    cosf = math.cos(true_anomaly)
    sinf = math.sin(true_anomaly)
    rel_r = [r * cosf, r * sinf, 0.0]
    rel_v = [vr * cosf - vt * sinf, vr * sinf + vt * cosf, 0.0]

    mass_sum = m1 + m2
    factor1 = -m2 / mass_sum
    factor2 = m1 / mass_sum

    r1 = [factor1 * x for x in rel_r]
    r2 = [factor2 * x for x in rel_r]
    v1 = [factor1 * x for x in rel_v]
    v2 = [factor2 * x for x in rel_v]
    return r1, v1, r2, v2


def orbital_elements_from_rel_state(
    rel_r: List[float],
    rel_v: List[float],
    m1: float,
    m2: float,
    g_const: float,
) -> Tuple[float, float]:
    """Return (period, ecc) from relative two-body phase-space state."""
    mu = g_const * (m1 + m2)
    rx, ry, rz = rel_r
    vx, vy, vz = rel_v
    r2 = rx * rx + ry * ry + rz * rz
    v2 = vx * vx + vy * vy + vz * vz
    r = math.sqrt(r2)

    hx = ry * vz - rz * vy
    hy = rz * vx - rx * vz
    hz = rx * vy - ry * vx

    ex = (vy * hz - vz * hy) / mu - rx / r
    ey = (vz * hx - vx * hz) / mu - ry / r
    ez = (vx * hy - vy * hx) / mu - rz / r
    ecc = math.sqrt(ex * ex + ey * ey + ez * ez)

    energy = 0.5 * v2 - mu / r
    semi = -mu / (2.0 * energy)
    period = 2.0 * math.pi * math.sqrt(semi * semi * semi / mu)
    return period, ecc


def validate_functional_binary_roundtrip(
    rows: List[Tuple[float, List[float], List[float]]],
    selected: List[Tuple[float, float, float, float]],
    velocity_unit: str,
) -> None:
    """Round-trip check from generated rows back to Kepler period/ecc.

    This guards against conversion drift when generating ic.raw.
    """
    if len(rows) < 6:
        raise ValueError("Need at least 6 rows to validate the first three binaries")

    def _to_pc_myr(vel: List[float]) -> List[float]:
        if velocity_unit == "km/s":
            return [x * KM_S_TO_PC_MYR for x in vel]
        return vel

    pair_indices = [(0, 1), (2, 3), (4, 5)]
    for pair_i, (ia, ib) in enumerate(pair_indices):
        m1, p_ref, e_ref = selected[pair_i][0], selected[pair_i][2], selected[pair_i][3]
        m2 = selected[pair_i][1]
        _, pa, va = rows[ia]
        _, pb, vb = rows[ib]

        rel_r = [pb[k] - pa[k] for k in range(3)]
        rel_v = [_to_pc_myr(vb)[k] - _to_pc_myr(va)[k] for k in range(3)]
        p_calc, e_calc = orbital_elements_from_rel_state(rel_r, rel_v, m1, m2, G_MSUN_PC_MYR)

        # Keep tolerance tight enough to catch real conversion mistakes, while
        # allowing tiny floating-point rounding differences.
        p_tol = max(1.0e-12, abs(p_ref) * 2.0e-4)
        e_tol = 2.0e-4
        if abs(p_calc - p_ref) > p_tol or abs(e_calc - e_ref) > e_tol:
            raise ValueError(
                "Binary round-trip check failed for pair "
                f"{pair_i + 1}: period(ref={p_ref:.16e}, calc={p_calc:.16e}), "
                f"ecc(ref={e_ref:.16e}, calc={e_calc:.16e}), "
                f"velocity_unit={velocity_unit}"
            )


def write_rows(output: Path, rows: List[Tuple[float, List[float], List[float]]]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as fh:
        for mass, pos, vel in rows:
            fh.write(
                "{:.14e} {:.14e} {:.14e} {:.14e} {:.14e} {:.14e} {:.14e}\n".format(
                    mass,
                    pos[0],
                    pos[1],
                    pos[2],
                    vel[0],
                    vel[1],
                    vel[2],
                )
            )


def build_t1(output: Path) -> None:
    m1 = 1.0
    m2 = 10.0
    semi = 0.1
    ecc = 0.9
    apo = semi * (1.0 + ecc)
    mass_sum = m1 + m2
    rel_r = [apo, 0.0, 0.0]
    rel_vy = math.sqrt(G_MSUN_PC_MYR * mass_sum * (1.0 - ecc) / (semi * (1.0 + ecc)))
    rel_v = [0.0, rel_vy, 0.0]

    factor1 = -m2 / mass_sum
    factor2 = m1 / mass_sum
    r1 = [factor1 * x for x in rel_r]
    r2 = [factor2 * x for x in rel_r]
    v1 = [factor1 * x for x in rel_v]
    v2 = [factor2 * x for x in rel_v]
    rows = [(m1, r1, v1), (m2, r2, v2)]
    write_rows(output, rows)


def build_t2(output: Path) -> None:
    m1 = 1.0
    m2 = 10.0
    semi = 0.01
    ecc = 0.6

    peri = semi * (1.0 - ecc)
    r1, v1, r2, v2 = two_body_peri_state(m1, m2, peri, ecc, G_MSUN_PC_MYR)

    rows = [(m1, r1, v1), (m2, r2, v2)]
    write_rows(output, rows)


def build_t3(output: Path) -> None:
    m1 = 1.0
    m2 = 0.5
    m3 = 0.05

    rin_peri = 0.002
    ein = 0.7
    r1_rel, v1_rel, r2_rel, v2_rel = two_body_peri_state(m1, m2, rin_peri, ein, G_MSUN_PC_MYR)

    m12 = m1 + m2
    rout_peri = 0.05
    eout = 0.6
    r12, v12, r3, v3 = two_body_peri_state(m12, m3, rout_peri, eout, G_MSUN_PC_MYR)

    p1 = [r12[i] + r1_rel[i] for i in range(3)]
    p2 = [r12[i] + r2_rel[i] for i in range(3)]
    v1 = [v12[i] + v1_rel[i] for i in range(3)]
    v2 = [v12[i] + v2_rel[i] for i in range(3)]

    rows = [(m1, p1, v1), (m2, p2, v2), (m3, r3, v3)]
    write_rows(output, rows)


def build_t4_with_outer_a(output: Path, outer_a: float) -> None:
    # Notebook-inspired hierarchical triple scales from "Tidal tensor test: 3-body".
    # Inner: (0.001, 0.009) with a=1e-3, e=0.9; outer: (0.01, 1.0) with e=0.01.
    # outer_a is configurable to support both the shortened and original-scale controls.
    m1 = 0.001
    m2 = 0.009
    m3 = 1.0

    rin_peri = 1.0e-4
    ein = 0.9
    r1_rel, v1_rel, r2_rel, v2_rel = two_body_peri_state(m1, m2, rin_peri, ein, G_MSUN_PC_MYR)

    m12 = m1 + m2
    rout_peri = outer_a * (1.0 - 0.01)
    eout = 0.01
    r12, v12, r3, v3 = two_body_peri_state(m12, m3, rout_peri, eout, G_MSUN_PC_MYR)

    p1 = [r12[i] + r1_rel[i] for i in range(3)]
    p2 = [r12[i] + r2_rel[i] for i in range(3)]
    v1 = [v12[i] + v1_rel[i] for i in range(3)]
    v2 = [v12[i] + v2_rel[i] for i in range(3)]

    rows: List[Tuple[float, List[float], List[float]]] = [(m1, p1, v1), (m2, p2, v2), (m3, r3, v3)]

    mass_tot = sum(item[0] for item in rows)
    com_pos = [sum(item[0] * item[1][k] for item in rows) / mass_tot for k in range(3)]
    com_vel = [sum(item[0] * item[2][k] for item in rows) / mass_tot for k in range(3)]

    recentered: List[Tuple[float, List[float], List[float]]] = []
    for mass, pos, vel in rows:
        pos_new = [pos[k] - com_pos[k] for k in range(3)]
        vel_new = [vel[k] - com_vel[k] for k in range(3)]
        recentered.append((mass, pos_new, vel_new))

    write_rows(output, recentered)


def build_t4(output: Path) -> None:
    # Shortened outer orbit for practical runtime.
    build_t4_with_outer_a(output, outer_a=0.15)


def build_t4_outer15(output: Path) -> None:
    # Original notebook-scale outer orbit.
    build_t4_with_outer_a(output, outer_a=1.5)


def build_functional_smoke(output: Path) -> None:
    rows: List[Tuple[float, List[float], List[float]]] = []

    # One deterministic close binary first, so binary-capable runs can use -b 1.
    m1 = 1.0
    m2 = 0.8
    semi = 0.01
    ecc = 0.2
    peri = semi * (1.0 - ecc)
    r1, v1, r2, v2 = two_body_peri_state(m1, m2, peri, ecc, G_MSUN_PC_MYR)
    rows.append((m1, r1, v1))
    rows.append((m2, r2, v2))

    # Add a small deterministic background cluster with zero net rotation.
    background = [
        (0.6, [0.20, 0.00, 0.00], [0.00, 0.10, 0.00]),
        (0.7, [-0.20, 0.00, 0.00], [0.00, -0.10, 0.00]),
        (0.5, [0.00, 0.22, 0.00], [-0.10, 0.00, 0.00]),
        (0.9, [0.00, -0.22, 0.00], [0.10, 0.00, 0.00]),
        (0.4, [0.00, 0.00, 0.18], [0.00, 0.05, -0.02]),
        (0.4, [0.00, 0.00, -0.18], [0.00, -0.05, 0.02]),
        (0.3, [0.15, 0.15, 0.00], [-0.06, 0.04, 0.00]),
        (0.3, [-0.15, -0.15, 0.00], [0.06, -0.04, 0.00]),
        (0.35, [0.15, -0.15, 0.00], [0.05, 0.03, 0.00]),
        (0.35, [-0.15, 0.15, 0.00], [-0.05, -0.03, 0.00]),
        (0.25, [0.10, 0.00, 0.12], [0.00, 0.04, 0.03]),
        (0.25, [-0.10, 0.00, -0.12], [0.00, -0.04, -0.03]),
        (0.20, [0.00, 0.10, -0.12], [-0.03, 0.00, 0.02]),
        (0.20, [0.00, -0.10, 0.12], [0.03, 0.00, -0.02]),
    ]
    rows.extend(background)

    mass_tot = sum(item[0] for item in rows)
    com_pos = [sum(item[0] * item[1][k] for item in rows) / mass_tot for k in range(3)]
    com_vel = [sum(item[0] * item[2][k] for item in rows) / mass_tot for k in range(3)]

    recentered: List[Tuple[float, List[float], List[float]]] = []
    for mass, pos, vel in rows:
        pos_new = [pos[k] - com_pos[k] for k in range(3)]
        vel_new = [vel[k] - com_vel[k] for k in range(3)]
        recentered.append((mass, pos_new, vel_new))

    write_rows(output, recentered)


def build_functional_bse_period_binary(output: Path, m1: float, m2: float, period_myr: float, ecc: float) -> None:
    """Build a two-body IC from BSE test-table parameters.

    The b0.dat binary table stores period in Myr and eccentricity.
    We reconstruct semi-major axis via Kepler's law in Msun/pc/Myr units,
    then place the binary at peri-center using the same deterministic setup
    as other two-body test builders.
    """
    semi = (G_MSUN_PC_MYR * (m1 + m2) * (period_myr / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0)
    peri = semi * (1.0 - ecc)
    r1, v1, r2, v2 = two_body_peri_state(m1, m2, peri, ecc, G_MSUN_PC_MYR)
    rows = [(m1, r1, v1), (m2, r2, v2)]
    write_rows(output, rows)


def build_functional_dual_merge_smoke(output: Path) -> None:
    """Two selected BSE binaries in one shared tiny model for all functional cases.

    Binary A (non-BH-BH merger candidate from data.bse.binary_merge):
      ids (17,18), m1=85.294690901635846, m2=23.528452122895846,
      period=4.594336e-09 Myr, ecc=0.499774915487773

    Binary B (GW-kick candidate from data.bse.gw_kick):
      ids (61,62), m1=35.947378474297203, m2=48.425671915280390,
      period=8.3652155e-08 Myr, ecc=0.427779218694198
    """
    rows: List[Tuple[float, List[float], List[float]]] = []

    # Build two binaries from period/ecc and place them apart in x to avoid overlap.
    # Binary A center offset
    b1_m1 = 85.294690901635846
    b1_m2 = 23.528452122895846
    b1_period = 4.594336e-09
    b1_ecc = 0.499774915487773
    b1_r1, b1_v1, b1_r2, b1_v2 = two_body_peri_state(
        b1_m1,
        b1_m2,
        (G_MSUN_PC_MYR * (b1_m1 + b1_m2) * (b1_period / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0) * (1.0 - b1_ecc),
        b1_ecc,
        G_MSUN_PC_MYR,
    )
    c1 = [-0.03, 0.0, 0.0]
    rows.append((b1_m1, [b1_r1[0] + c1[0], b1_r1[1] + c1[1], b1_r1[2] + c1[2]], b1_v1))
    rows.append((b1_m2, [b1_r2[0] + c1[0], b1_r2[1] + c1[1], b1_r2[2] + c1[2]], b1_v2))

    # Binary B center offset
    b2_m1 = 35.947378474297203
    b2_m2 = 48.425671915280390
    b2_period = 8.3652155e-08
    b2_ecc = 0.427779218694198
    b2_r1, b2_v1, b2_r2, b2_v2 = two_body_peri_state(
        b2_m1,
        b2_m2,
        (G_MSUN_PC_MYR * (b2_m1 + b2_m2) * (b2_period / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0) * (1.0 - b2_ecc),
        b2_ecc,
        G_MSUN_PC_MYR,
    )
    c2 = [0.03, 0.0, 0.0]
    rows.append((b2_m1, [b2_r1[0] + c2[0], b2_r1[1] + c2[1], b2_r1[2] + c2[2]], b2_v1))
    rows.append((b2_m2, [b2_r2[0] + c2[0], b2_r2[1] + c2[1], b2_r2[2] + c2[2]], b2_v2))

    # Add a light symmetric background so base/galpy smoke remains robust.
    background = [
        (0.5, [0.20, 0.00, 0.00], [0.00, 0.06, 0.00]),
        (0.5, [-0.20, 0.00, 0.00], [0.00, -0.06, 0.00]),
        (0.4, [0.00, 0.20, 0.00], [-0.06, 0.00, 0.00]),
        (0.4, [0.00, -0.20, 0.00], [0.06, 0.00, 0.00]),
    ]
    rows.extend(background)

    # Recenter to zero COM position/velocity.
    mass_tot = sum(item[0] for item in rows)
    com_pos = [sum(item[0] * item[1][k] for item in rows) / mass_tot for k in range(3)]
    com_vel = [sum(item[0] * item[2][k] for item in rows) / mass_tot for k in range(3)]

    recentered: List[Tuple[float, List[float], List[float]]] = []
    for mass, pos, vel in rows:
        pos_new = [pos[k] - com_pos[k] for k in range(3)]
        vel_new = [vel[k] - com_vel[k] for k in range(3)]
        recentered.append((mass, pos_new, vel_new))

    write_rows(output, recentered)


def build_functional_mcluster_dual_merge_smoke(output: Path) -> None:
    """Generate a physical background with mcluster and inject target binaries.

    Steps:
     1) Use mcluster to generate N=100, Plummer, Rh=1 pc, Kroupa IMF, no binaries (fixed seed).
    2) Replace the first six stars by three binaries (6 stars total).
        Binary masses/period/ecc are script-fixed (with optional env overrides).
     3) Replace stars #7/#8 with a hyperbolic merger pair (m=10,10; a=-1e-8 pc; e=1.001)
         in inbound configuration so the relative distance is shrinking.
     4) Force the ninth star mass to 27 Msun for single-star SN-kick coverage.
     5) Align each injected pair center of mass with one of the outermost stars.
     6) Optionally recenter to zero net COM position/velocity (disabled by default).

    Output convention for this builder:
     - position: pc
    - velocity: km/s (configurable only via MCLUSTER_OUTPUT_VELOCITY_UNIT)
     so that petar.init should use -v=KM_S_TO_PC_MYR under -u 1.
    """
    if MCLUSTER_OUTPUT_VELOCITY_UNIT not in {"km/s", "pc/Myr"}:
        raise ValueError(
            f"Unsupported MCLUSTER_OUTPUT_VELOCITY_UNIT={MCLUSTER_OUTPUT_VELOCITY_UNIT}; "
            "expected 'km/s' or 'pc/Myr'"
        )

    with tempfile.TemporaryDirectory(prefix="petar_make_ic_") as tmpdir:
        tmp_prefix_name = "mcluster_seed"
        cmd = [
            "mcluster",
            "-N",
            "100",
            "-P",
            "0",
            "-R",
            "1",
            "-f",
            "1",
            "-b",
            "0",
            "-u",
            "1",
            "-C",
            "3",
            "-s",
            str(MCLUSTER_SEED),
            "-o",
            tmp_prefix_name,
        ]
        subprocess.run(
            cmd,
            cwd=tmpdir,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        table_path = Path(tmpdir) / f"{tmp_prefix_name}.txt"
        rows: List[Tuple[float, List[float], List[float]]] = []
        with table_path.open("r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 7:
                    continue
                mass = float(parts[0])
                pos = [float(parts[1]), float(parts[2]), float(parts[3])]
                vel_raw = [float(parts[4]), float(parts[5]), float(parts[6])]
                if MCLUSTER_OUTPUT_VELOCITY_UNIT == "km/s":
                    vel = vel_raw
                else:
                    vel = [x * KM_S_TO_PC_MYR for x in vel_raw]
                rows.append((mass, pos, vel))

    if len(rows) < 9:
        raise RuntimeError(
            "mcluster returned fewer than 9 stars; cannot inject three binaries, "
            "one hyperbolic merger pair, and one single-star mass override"
        )

    # Replace the outermost stars first so the injected binaries start in the
    # weak-field outskirts rather than near the cluster center.
    rows.sort(key=lambda item: math.sqrt(sum(coord * coord for coord in item[1])), reverse=True)

    selected = resolve_functional_binary_params()
    if len(selected) < 3:
        raise ValueError(f"Need at least three binaries in FUNCTIONAL_FIXED_BINARIES, got {len(selected)}")
    b1_m1, b1_m2, b1_period, b1_ecc = selected[0]
    b2_m1, b2_m2, b2_period, b2_ecc = selected[1]
    b3_m1, b3_m2, b3_period, b3_ecc = selected[2]

    def binary_rel_state(m1: float, m2: float, period_myr: float, ecc: float) -> Tuple[List[float], List[float], List[float], List[float]]:
        semi = (G_MSUN_PC_MYR * (m1 + m2) * (period_myr / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0)
        return two_body_apo_state(m1, m2, semi, ecc, G_MSUN_PC_MYR)

    b1_r1, b1_v1, b1_r2, b1_v2 = binary_rel_state(b1_m1, b1_m2, b1_period, b1_ecc)
    b2_r1, b2_v1, b2_r2, b2_v2 = binary_rel_state(b2_m1, b2_m2, b2_period, b2_ecc)
    b3_r1, b3_v1, b3_r2, b3_v2 = binary_rel_state(b3_m1, b3_m2, b3_period, b3_ecc)

    h_m1 = 10.0
    h_m2 = 10.0
    h_semi = -1.0e-8
    h_ecc = 1.001
    # Start far on inbound branch (near asymptote) to keep initial separation
    # larger while maintaining dr/dt < 0 and the same (a, e).
    h_true_anomaly = -3.0955
    h_r1, h_v1, h_r2, h_v2 = two_body_hyperbolic_inbound_state(
        h_m1,
        h_m2,
        h_semi,
        h_ecc,
        h_true_anomaly,
        G_MSUN_PC_MYR,
    )

    if MCLUSTER_OUTPUT_VELOCITY_UNIT == "km/s":
        # Convert binary velocities to km/s to match mcluster output unit.
        b1_v1 = [v * PC_MYR_TO_KM_S for v in b1_v1]
        b1_v2 = [v * PC_MYR_TO_KM_S for v in b1_v2]
        b2_v1 = [v * PC_MYR_TO_KM_S for v in b2_v1]
        b2_v2 = [v * PC_MYR_TO_KM_S for v in b2_v2]
        b3_v1 = [v * PC_MYR_TO_KM_S for v in b3_v1]
        b3_v2 = [v * PC_MYR_TO_KM_S for v in b3_v2]
        h_v1 = [v * PC_MYR_TO_KM_S for v in h_v1]
        h_v2 = [v * PC_MYR_TO_KM_S for v in h_v2]

    # Anchor binary COMs to three outermost stars directly.
    # Using pair COM can place a binary near the center when the anchors
    # are on opposite sides of the cluster.
    pcom1, vcom1 = list(rows[0][1]), list(rows[0][2])
    pcom2, vcom2 = list(rows[1][1]), list(rows[1][2])
    pcom3, vcom3 = list(rows[4][1]), list(rows[4][2])
    pcom4, vcom4 = list(rows[7][1]), list(rows[7][2])

    rows[0] = (b1_m1, [pcom1[k] + b1_r1[k] for k in range(3)], [vcom1[k] + b1_v1[k] for k in range(3)])
    rows[1] = (b1_m2, [pcom1[k] + b1_r2[k] for k in range(3)], [vcom1[k] + b1_v2[k] for k in range(3)])
    rows[2] = (b2_m1, [pcom2[k] + b2_r1[k] for k in range(3)], [vcom2[k] + b2_v1[k] for k in range(3)])
    rows[3] = (b2_m2, [pcom2[k] + b2_r2[k] for k in range(3)], [vcom2[k] + b2_v2[k] for k in range(3)])
    rows[4] = (b3_m1, [pcom3[k] + b3_r1[k] for k in range(3)], [vcom3[k] + b3_v1[k] for k in range(3)])
    rows[5] = (b3_m2, [pcom3[k] + b3_r2[k] for k in range(3)], [vcom3[k] + b3_v2[k] for k in range(3)])
    rows[6] = (h_m1, [pcom4[k] + h_r1[k] for k in range(3)], [vcom4[k] + h_v1[k] for k in range(3)])
    rows[7] = (h_m2, [pcom4[k] + h_r2[k] for k in range(3)], [vcom4[k] + h_v2[k] for k in range(3)])
    rows[8] = (27.0, list(rows[8][1]), list(rows[8][2]))

    rows_out = rows
    if FUNCTIONAL_RECENTER_AFTER_INJECTION:
        mass_tot = sum(item[0] for item in rows)
        com_pos = [sum(item[0] * item[1][k] for item in rows) / mass_tot for k in range(3)]
        com_vel = [sum(item[0] * item[2][k] for item in rows) / mass_tot for k in range(3)]

        recentered: List[Tuple[float, List[float], List[float]]] = []
        for mass, pos, vel in rows:
            recentered.append(
                (
                    mass,
                    [pos[k] - com_pos[k] for k in range(3)],
                    [vel[k] - com_vel[k] for k in range(3)],
                )
            )
        rows_out = recentered

    validate_functional_binary_roundtrip(rows_out, selected, MCLUSTER_OUTPUT_VELOCITY_UNIT)
    write_rows(output, rows_out)


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate deterministic initial conditions for PeTar validation scenarios")
    parser.add_argument(
        "--case",
        choices=[
            "t1",
            "t2",
            "t3",
            "t4",
            "t4_outer15",
            "functional_smoke",
            "functional_dual_merge_smoke",
            "functional_mcluster_dual_merge_smoke",
        ],
        required=True,
    )
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    out = Path(args.output)
    if args.case == "t1":
        build_t1(out)
    elif args.case == "t2":
        build_t2(out)
    elif args.case == "t3":
        build_t3(out)
    elif args.case == "t4":
        build_t4(out)
    elif args.case == "t4_outer15":
        build_t4_outer15(out)
    elif args.case == "functional_smoke":
        build_functional_smoke(out)
    elif args.case == "functional_dual_merge_smoke":
        build_functional_dual_merge_smoke(out)
    elif args.case == "functional_mcluster_dual_merge_smoke":
        build_functional_mcluster_dual_merge_smoke(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())