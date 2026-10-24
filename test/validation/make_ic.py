#!/usr/bin/env python3
import argparse
import json
import math
import os
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Tuple


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
# length conversions used by the T9 builder
AU_PC = 4.8481e-6
RSUN_PC = AU_PC / 215.032
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


def build_t5(output: Path) -> None:
    # Massive-binary + light-churner configuration reproducing the 2026-10-06
    # tidal-tensor momentum-leak event (run_N1k_200Myr at t~47.6):
    # a 77.65 Msun binary (46.06 + 31.59) with a=0.0081 pc, e=0.667, plus a
    # 0.242 Msun star on an eccentric orbit whose pericenter (0.01 pc) sits
    # inside the group radius and apocenter (0.04 pc) outside it, so the AR
    # group membership churns every outer orbit and close passages dominate
    # the tidal-tensor fit (the T3 ~ r^-4 regime).
    m1 = 46.06
    m2 = 31.59
    m3 = 0.242

    ein = 0.667
    ain = 0.0081
    # start the inner binary at apocenter
    r12, v12, r3_in, v3_in = two_body_apo_state(m1, m2, ain, ein, G_MSUN_PC_MYR)

    m12 = m1 + m2
    aout = 0.025
    eout = 0.6
    # start the outer orbit at apocenter
    r12o, v12o, r3o, v3o = two_body_apo_state(m12, m3, aout, eout, G_MSUN_PC_MYR)

    p1 = [r12o[i] + r12[i] for i in range(3)]
    p2 = [r12o[i] + r3_in[i] for i in range(3)]
    v1 = [v12o[i] + v12[i] for i in range(3)]
    v2 = [v12o[i] + v3_in[i] for i in range(3)]

    rows: List[Tuple[float, List[float], List[float]]] = [(m1, p1, v1), (m2, p2, v2), (m3, r3o, v3o)]

    mass_tot = sum(item[0] for item in rows)
    com_pos = [sum(item[0] * item[1][k] for item in rows) / mass_tot for k in range(3)]
    com_vel = [sum(item[0] * item[2][k] for item in rows) / mass_tot for k in range(3)]

    recentered: List[Tuple[float, List[float], List[float]]] = []
    for mass, pos, vel in rows:
        pos_new = [pos[k] - com_pos[k] for k in range(3)]
        vel_new = [vel[k] - com_vel[k] for k in range(3)]
        recentered.append((mass, pos_new, vel_new))

    write_rows(output, recentered)


def build_t7(output: Path, rp: float) -> None:
    # Knife-edge hierarchical trio + one static perturber, reproducing the
    # 2026-10-19/20 rate-type ds-controller collapse (record:
    # SDAR/docs/etot_sd_ref_corruption_fix_notes.md). The trio state is the exact
    # c.m.-subtracted configuration extracted from the functional-bse grind
    # dump at t=8.171875 (hard-dump records 0, 1, 68): an ultra-tight
    # 42.7 Msun inner binary (a = 8.34e-7, e ~ 5e-4) plus a 0.084 Msun
    # intruder barely unbound (e = 1.0012) lingering at ~4e4 a_in with
    # v ~ v_esc. A single 1.0 Msun perturber placed at distance rp on +x
    # (v = 0) activates the changeover-modulated perturbation channel when
    # rp falls inside the members' changeover transition band; with the
    # solver radii used by the scenario (mean outer ~2.17, intruder band
    # [0.082, 0.545]) rp = 0.15 sits mid-band of the intruder.
    trio = [
        (23.660928497887493,
         [-5.3337175601186892e-05, 4.6028170444190053e-05, -1.1970474815936605e-05],
         [31.05023827013645, -211.76010591288758, -0.0021811559681917392]),
        (19.044431988547881,
         [-5.2512178304042756e-05, 4.6149125184635764e-05, -1.1970474986910951e-05],
         [-38.587592282131624, 263.10022998850241, -0.0028319265168438856]),
        (0.08396199811790174,
         [0.026941613553335575, -0.023438617110016713, 0.0060885097654825859],
         [2.3976235834667845, -1.8307267677018166, 1.2570044748639835]),
    ]
    rows = list(trio) + [(1.0, [rp, 0.0, 0.0], [0.0, 0.0, 0.0])]
    # The knife-edge intruder state (barely unbound, lingering at ~4e4 a_in)
    # is precision-sensitive: a 15-digit rounding shifts the trajectory enough
    # to move the perturber across the trigger-zone boundary. Write full
    # double precision instead of the shared write_rows (.14e).
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as fh:
        for mass, pos, vel in rows:
            fh.write("{:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e}\n".format(
                mass, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]))


def _write_t8_rows(output: Path, rows) -> None:
    # Full double precision (mid-encounter states are precision-sensitive,
    # same rationale as build_t7).
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as fh:
        for mass, pos, vel in rows:
            fh.write("{:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e}\n".format(
                mass, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]))


def build_t8_cascade(output: Path) -> None:
    # T8 hierarchical-cascade case: the 4-body group (members 2,3,6,7) from
    # the 2026-10-22 N=1000 benchmark dump
    # current.dump_large_step_h4_n30_g10_t5.929199 (bench-n1k-arfix/) -- a
    # 105 Msun tight binary (a=6.5e-7 pc, e=0.016, P=4.8e-9 Myr) colliding
    # with a 0.3 Msun binary at mutual pericenter (root q~0.0063 pc,
    # e~0.94) -- PLUS the 7 bystanders carrying 95% of the cluster's tidal
    # field on the group (minimal-set selection by G*m/r^3; the 4-body-only
    # extract integrates cheaply in BOTH gauges because the slowdown
    # disengagement in the original event was driven by this external
    # field). Dump-replay gauge A/B: LogH 2.2M-11.4M AR steps per hard
    # block vs BTLogH 12k-35k (38-934x cheaper) at comparable energy
    # monitors -- the flagship hierarchical-multiple win for the product
    # (BTLogH) gauge.
    _write_t8_rows(output, [
        (2.16290647351530012e+01, [2.51456703316543795e-03, -2.06976591139084176e-02, -1.24765142852054228e-02], [5.34463233418060426e+01, 2.00370466881110872e+02, 6.36228792189857245e+02]),
        (8.37422646967909969e+01, [2.51505813667540767e-03, -2.06980872336284261e-02, -1.24764239685754286e-02], [-1.39708701042149563e+01, -5.15929694683241280e+01, -1.64618485403352736e+02]),
        (1.98992673330040004e-01, [-3.91263914796458057e-03, -1.86856872708184296e-02, -1.38420936072554263e-02], [-3.45386068613545705e+00, -2.21391618967252501e+00, -1.31479461483637401e+01]),
        (9.91972549158660066e-02, [-3.92396056996457965e-03, -1.86883438878684416e-02, -1.38191147460754171e-02], [-7.30576172128457269e-01, -5.78096101307962495e+00, -6.61668487769453950e+00]),
        (1.70233580713849997e+00, [-2.29160220959445660e-02, -7.88579022022842246e-03, 3.73382620227458162e-03], [-2.07882505443205723e+00, 5.81302177774244866e-01, 2.76371456685926020e+00]),
        (2.38615570113539999e+00, [3.63964259318254046e-02, -5.98916917792184367e-02, -1.55540272379354194e-02], [2.75358403890394277e+00, -8.21854161516582593e+00, 1.07623175087182599e+01]),
        (2.35168022272280008e-01, [-2.29970252705445910e-02, -8.20403977885841851e-03, 4.05274885988457512e-03], [-2.85936956868995695e+00, 1.15889856014367476e+00, 3.04059771164276027e+00]),
        (8.13137780849759961e-01, [3.64036620797854038e-02, -5.98930949320784500e-02, -1.55608487146554203e-02], [1.21219591270884286e+00, 1.69395577018388757e+01, -2.30428065234047388e+01]),
        (1.67223535414359997e-01, [2.88122460499053989e-02, -3.61299796138442320e-04, -1.45111445604654277e-02], [-2.62583949128975691e+00, -2.70095300223775148e-01, -1.33976819446123985e+00]),
        (1.15714268543820005e+01, [-1.52601078412545887e-02, 1.07215364418017572e-01, 6.10829613892335740e-02], [2.49651036736204279e+00, 1.72284990067220861e+02, 1.54326432434967245e+02]),
        (1.07534750297469994e+01, [-1.52604837545845806e-02, 1.07215514475368576e-01, 6.10827990160735801e-02], [-1.58794043906265703e+00, -1.86103864618319136e+02, -1.64612150153052738e+02]),
    ])


def build_t8_chaotic(output: Path) -> None:
    # T8 chaotic case: CM-recentered extract of the step-dominant 6-body
    # group (members 12,13,10,11,0,1) from the pre-fix pathological cluster
    # baseline.dump_large_step_h4_n31_g10_t0.359375: the same 105 Msun
    # binary embedded in a dense 4+2 exchange complex with low-mass
    # stars. Dump-replay gauge A/B on the parent cluster: LogH 65,791
    # steps with 16 hard-energy-significant events vs BTLogH 264,956
    # steps with 0 events (BTLogH pays ~4x steps for a clean monitor) --
    # the non-hierarchical counterpoint to the cascade case.
    _write_t8_rows(output, [
        (1.55182606135500001e-01, [3.62853972415677628e-03, 1.86153494924545804e-03, -1.27761703803046592e-03], [3.83911750441267730e+00, -1.21808834763359783e+01, 5.91010231527986463e+00]),
        (1.07671330855760000e-01, [3.79843945169677788e-03, 1.91898265402545842e-03, -1.32666923963446759e-03], [4.21031695253597782e+00, -1.37506708736609777e+01, 6.99897222938166497e+00]),
        (1.32761026988319997e+00, [-6.35014888314322512e-03, 3.55839192234547186e-03, 6.31244862344253199e-03], [8.66770225115347692e+00, 3.35966317437472117e+00, -2.55865346753363543e+00]),
        (3.52232555331649988e-01, [-6.38193367854322791e-03, 3.53395075165546801e-03, 6.31580711190753152e-03], [3.61388144984747717e+00, 1.05851285368287211e+01, -6.53098219820803561e+00]),
        (8.37422646967909969e+01, [9.20182112167378108e-05, -6.12672555245219463e-05, -9.74345360444661690e-05], [4.06546146316437795e+00, 6.22445584967110221e+01, 1.68150033174157358e+02]),
        (2.16290647351530012e+01, [9.24942905368153490e-05, -6.16655741745297803e-05, -9.73046052654691351e-05], [-1.63798253051127247e+01, -2.41217930288848976e+02, -6.50848206974502659e+02]),
    ])


def build_t8_control(output: Path) -> None:
    # T8 control case: two detached binaries (CM separation 0.093 pc, no
    # encounter within the window) from
    # baseline.dump_large_step_h4_n32_g13_t0.336914 members (2,3)+(6,7).
    # Both gauges must stay cheap and near-equal (parity/false-positive
    # guard for the cascade and chaotic assertions).
    _write_t8_rows(output, [
        (8.37422646967909969e+01, [2.15118834008598547e-04, -7.70364143333202023e-05, 2.89664571128216142e-05], [2.65766536269051379e+01, 4.26192580230517422e+01, 1.71832257060643343e+02]),
        (2.16290647351530012e+01, [2.15585544198648726e-04, -7.74666156233005943e-05, 2.89930402388206832e-05], [-1.02931874148616856e+02, -1.65006859015856264e+02, -6.65291905658366659e+02]),
        (1.55182606135500001e-01, [-8.62791762918013738e-02, 3.08621417984666913e-02, -1.15757898937789774e-02], [1.88065702779903843e+00, -4.98562568410358731e-01, 2.54293672240903879e-01]),
        (1.07671330855760000e-01, [-8.62664603868213975e-02, 3.09967994524367185e-02, -1.16692884659742781e-02], [4.07872020828773874e+00, -1.09434054840058748e-01, -1.86873478169156115e-01]),
    ])


def build_t4_outer15(output: Path) -> None:
    # Original notebook-scale outer orbit.
    build_t4_with_outer_a(output, outer_a=1.5)


# ----- T9: micro-TDE trigger criterion (BH/NS + star pair) -----
# Template masses follow Rastello et al. (2026) event 1 (BH 20.059 Msun,
# MS star 2.3934 Msun, Z = 0.02). All runs use --bse-metallicity 0.02, so the
# BSE ZAMS radius of the star is deterministic; T9_R_STAR_ZAMS is the measured
# value and only fixes the beta grid (the report re-derives r_t from the
# logged pre-merger star radius, so a BSE radius drift shifts design and
# logged beta together and is caught by the beta_ratio_dev check).
T9_M_BH = 20.059
T9_M_NS = 1.4
T9_M_STAR = 2.3934
T9_R_STAR_ZAMS = 1.7543175952
T9_SEMI_AU = 100.0          # bound grid semi-major axis
T9_V_INF_KMS = 5.0          # hyperbolic grid v_inf
T9_R0_AU = 100.0            # hyperbolic grid start separation
T9_MSMS_RPERI_RSUN = 10.0   # MS-MS control pericentre (no CO -> no TDE, no contact)
T9_BETA_GRID_BOUND = [0.3, 0.6, 0.9, 1.1, 1.5, 3.0, 5.0]
T9_BETA_GRID_HYPERBOLIC = [0.5, 0.9, 1.1, 2.93, 5.0]
T9_BETA_NS = 2.0
# perturbed-candidate scenario (Rastello comm.): wide high-e binary that is a
# TDE candidate by osculating r_p < r_t, perturbed by a stellar flyby during
# the approach so the pericentre is lifted above r_t before the plunge; the
# run must NOT trigger (old versions fired prematurely). The control (no
# perturber) must fire exactly at the analytic pericentre epoch.
T9_PC_BETA = 1.5
T9_PC_SEMI_AU = 380.0
T9_PC_PERT_M = 1.0
T9_PC_PERT_PERI_AU = 250.0
T9_PC_PERT_A_ABS_AU = 1500.0
T9_PC_PERT_TPERI_MYR = 4.0e-4  # flyby pericentre ~400 yr after start


def tde_tidal_radius_rsun(m_co: float, m_star: float, r_star_rsun: float) -> float:
    """Rees (1988) tidal radius r_t = R_* (M_CO/M_*)^(1/3) [Rsun]."""
    return r_star_rsun * (m_co / m_star) ** (1.0 / 3.0)


def write_rows_typed(output: Path, rows: List[Tuple[float, List[float], List[float], int]]) -> None:
    """Write mass/pos/vel (pc, pc/Myr) + BSE stellar type column (8th)."""
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as fh:
        for mass, pos, vel, kw in rows:
            fh.write(
                "{:.14e} {:.14e} {:.14e} {:.14e} {:.14e} {:.14e} {:.14e} {}\n".format(
                    mass, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2], kw
                )
            )


def _build_t9_perturbed_candidate(output: Path, kind: str) -> Dict[str, float]:
    """Perturbed-candidate scenario: BH-star wide binary (a=380 AU) with
    osculating r_p < r_t (beta=1.5), first pericentre at P/2 ~ 781.7 yr.

    pc_control   : binary + bystander only -> TDE must fire at the pericentre
    pc_perturbed : + 1 Msun flyby (pericentre 250 AU from the pair CM at
                   ~400 yr, tilted plane) lifting r_p above r_t before the
                   plunge -> NO TDE may fire (premature-trigger regression)
    """
    m_co, m_star = T9_M_BH, T9_M_STAR
    r_t = tde_tidal_radius_rsun(m_co, m_star, T9_R_STAR_ZAMS)
    r_peri_rsun = r_t / T9_PC_BETA
    semi_pc = T9_PC_SEMI_AU * AU_PC
    ecc = 1.0 - (r_peri_rsun * RSUN_PC) / semi_pc
    # first pericentre epoch from the apocentre start (exact two-body)
    t_peri_myr = math.pi * math.sqrt(semi_pc ** 3 / (G_MSUN_PC_MYR * (m_co + m_star)))

    r1, v1, r2, v2 = two_body_apo_state(m_co, m_star, semi_pc, ecc, G_MSUN_PC_MYR)
    rows = [(m_co, r1, v1, 14), (m_star, r2, v2, 1)]

    if kind == "pc_perturbed":
        # perturber flyby around the pair CM (at origin): solve the inbound
        # true anomaly so its pericentre epoch is T9_PC_PERT_TPERI_MYR
        mt = m_co + m_star
        q_pc = T9_PC_PERT_PERI_AU * AU_PC
        a_abs = T9_PC_PERT_A_ABS_AU * AU_PC
        ecc_p = 1.0 + q_pc / a_abs
        mu = G_MSUN_PC_MYR * (mt + T9_PC_PERT_M)

        def _tp(nu: float) -> float:
            p = a_abs * (ecc_p * ecc_p - 1.0)
            r = p / (1.0 + ecc_p * math.cos(nu))
            # cosh H = (1 - r/semi)/e with semi = -a_abs (positive magnitude here)
            h_anom = math.acosh(max(1.0, (1.0 + r / a_abs) / ecc_p))
            m_anom = ecc_p * math.sinh(h_anom) - h_anom
            return m_anom / math.sqrt(mu / a_abs ** 3)

        lo, hi = -math.acos(-1.0 / ecc_p) + 1.0e-6, -1.0e-4
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if _tp(mid) > T9_PC_PERT_TPERI_MYR:
                lo = mid
            else:
                hi = mid
        nu_p = 0.5 * (lo + hi)
        _r1, _v1, pr, pv = two_body_hyperbolic_inbound_state(mt, T9_PC_PERT_M, -a_abs, ecc_p, nu_p, G_MSUN_PC_MYR)
        # tilt the flyby plane 90 deg (binary lives in x-y; flyby in x-z)
        pr, pv = [pr[0], -pr[2], pr[1]], [pv[0], -pv[2], pv[1]]
        rows.append((T9_PC_PERT_M, pr, pv, 1))

    rows.append((0.5, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1))
    _assert_no_coincident_positions([(m, p, v) for m, p, v, _ in rows], f"t9_tde_trigger({kind})")
    write_rows_typed(output, rows)

    design = {
        "kind": kind,
        "beta_design": float(T9_PC_BETA),
        "m_co": m_co,
        "m_star": m_star,
        "r_star_ref_rsun": T9_R_STAR_ZAMS,
        "r_tidal_design_rsun": r_t,
        "r_peri_design_rsun": r_peri_rsun,
        "r_init_au": semi_pc / AU_PC * (1.0 + ecc),
        "t_pericentre_myr": t_peri_myr,
        "expect_tde": kind == "pc_control",
        "metallicity": 0.02,
    }
    design_path = output.with_suffix("").with_name(output.name + ".design.json")
    design_path.write_text(json.dumps(design, indent=1), encoding="utf-8")
    return design


# events from the author reference sample replayed from their last pre-event
# snapshot (see test/validation/data/micro_TDE/README.md). Only events whose
# snapshot state quantitatively reproduces the author trigger are encoded:
# event 9 (E_2+2, r_p to 0.2%) and event 10 (E_2+2, r_p to 0.6%, a
# phase-affected). Events 1/3/4/5/7 do not reproduce from the stripped
# few-body state (missing cluster field; event 1 additionally has an internal
# velocity/position inconsistency in the reference CSV), events 2/8 exceed the
# suite runtime budget, event 6 produces instant-merger artifacts.
T9_REPLAY_EVENTS = [9, 10]
T9_REPLAY_DIR = Path(__file__).resolve().parent / "data" / "micro_TDE"


def _t9_replay_reference(event_id: int) -> Tuple[Dict[str, str], List[Dict[str, str]], float, float]:
    """Author reference row + last pre-event snapshot rows for an event."""
    import csv

    ev_row = None
    with (T9_REPLAY_DIR / "micro_TDE_events.csv").open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if int(row["event"]) == event_id:
                ev_row = row
                break
    if ev_row is None:
        raise ValueError(f"micro_TDE event {event_id} not found")

    snaps: Dict[int, List[Dict[str, str]]] = {}
    with (T9_REPLAY_DIR / "micro_TDE_snapshots.csv").open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if int(row["event"]) == event_id:
                snaps.setdefault(int(row["snap"]), []).append(row)
    last = max(snaps)
    rows = snaps[last]
    t_snap = float(rows[0]["t_snap_Myr"])
    dt_before = float(rows[0]["dt_before_TDE_Myr"])
    return ev_row, rows, t_snap, dt_before


def build_t9_tde_event_replay(output: Path, event_id: int) -> Dict[str, Any]:
    """Replay IC for one author micro-TDE event from its last pre-event snapshot.

    Writes the 9-column base IC (mass/pos/vel [Msun, pc, pc/Myr] + BSE type +
    stellar age t_snap) plus a design.json carrying the author trigger
    reference (a, e, r_peri, beta). The converted `input` still needs the
    clock alignment patch (see fix_t9_replay_input).
    """
    ev_row, rows, t_snap, dt_before = _t9_replay_reference(event_id)
    ic_rows = []
    for r in rows:
        ic_rows.append((
            float(r["mass"]),
            [float(r[c]) for c in "xyz"],
            [float(r["v" + c]) for c in "xyz"],
            int(r["type"]),
        ))
    _assert_no_coincident_positions(ic_rows, f"t9_tde_trigger(replay{event_id})")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as fh:
        for mass, pos, vel, kw in ic_rows:
            fh.write(
                "{:.14e} {:.14e} {:.14e} {:.14e} {:.14e} {:.14e} {:.14e} {} {}\n".format(
                    mass, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2], kw, t_snap
                )
            )

    design = {
        "kind": f"replay{event_id}",
        "beta_design": float(ev_row["beta"]),
        "m_co": float(ev_row["M_BH_Msun"]),
        "m_star": float(ev_row["M_star_Msun"]),
        "r_star_ref_rsun": float(ev_row["R_star_Rsun"]),
        "r_tidal_design_rsun": float(ev_row["r_tidal_Rsun"]),
        "r_peri_design_rsun": float(ev_row["r_peri_Rsun"]),
        "r_init_au": float("nan"),  # set after conversion (needs pair state)
        "t_pericentre_myr": float("nan"),
        "expect_tde": True,
        "author_a_rsun": float(ev_row["a_Rsun"]),
        "author_ecc": float(ev_row["ecc"]),
        "author_t_snap_myr": t_snap,
        "author_dt_before_myr": dt_before,
        "metallicity": 0.02,
    }
    design_path = output.with_suffix("").with_name(output.name + ".design.json")
    design_path.write_text(json.dumps(design, indent=1), encoding="utf-8")
    return design


def fix_t9_replay_input(input_path: Path, t_snap: float) -> None:
    """Align a converted replay IC with the snapshot epoch.

    petar.init --time sets the stellar ages (tphys = t_snap), but the header
    time and the per-particle t_record / t_interrupt (columns 11/12) stay 0;
    the hard integrator then asserts star.tphys <= time_origin + window and
    the first SSE catch-up overshoots. Patch: header time = t_snap and
    t_record = t_interrupt = t_snap for every row.
    """
    path = Path(input_path)
    lines = path.read_text(encoding="utf-8").splitlines()
    n = len(lines) - 1
    lines[0] = f"0 {n} {t_snap:.10f}"
    for i in range(1, len(lines)):
        c = lines[i].split()
        if len(c) >= 12:
            c[10] = f"{t_snap:.10f}"
            c[11] = f"{t_snap:.10f}"
            lines[i] = " ".join(c)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_t9_tde_trigger(output: Path, kind: str, beta: float) -> Dict[str, float]:
    """Generate one T9 micro-TDE trigger-criterion IC (3 particles).

    Rows: (CO, star, bystander 0.5 Msun at 1 pc). Bound kinds mark the first
    two rows as the primordial binary (petar -b 1); the hyperbolic kind uses
    two singles (group detection links them). The 8th column is the BSE type
    (petar.init --type '$8'): 14 BH, 13 NS, 1 MS.

    Also writes <output>.design.json with the design parameters for the
    validation extract/report.
    """
    if kind not in {"bound", "hyperbolic", "ns", "msms", "pc_control", "pc_perturbed"}:
        raise ValueError(f"Unknown T9 kind: {kind}")

    if kind in ("pc_control", "pc_perturbed"):
        return _build_t9_perturbed_candidate(output, kind)

    m_co = T9_M_NS if kind == "ns" else T9_M_BH
    kw_co = 13 if kind == "ns" else 14
    m_star, kw_star = T9_M_STAR, 1
    r_t = tde_tidal_radius_rsun(m_co, m_star, T9_R_STAR_ZAMS)

    if kind == "msms":
        # no CO -> getTDERadius returns 0; use a pericentre far from contact
        r_peri_rsun = T9_MSMS_RPERI_RSUN
        beta_design = float("nan")
        m_co, kw_co = T9_M_STAR, 1  # second star in place of the CO
    else:
        r_peri_rsun = r_t / beta
        beta_design = float(beta)

    semi_pc = T9_SEMI_AU * AU_PC
    if kind == "hyperbolic":
        # unbound orbit from (v_inf, r0): a = -GM/v_inf^2, e = 1 + r_p/|a|,
        # inbound true anomaly solved for r(nu) = r0 (the near-parabolic
        # pericentre-scale orbit makes r0 = O(100 AU) already "far")
        mt = m_co + m_star
        v_inf = T9_V_INF_KMS * KM_S_TO_PC_MYR
        a_abs = G_MSUN_PC_MYR * mt / (v_inf * v_inf)
        ecc = 1.0 + (r_peri_rsun * RSUN_PC) / a_abs
        p = a_abs * (ecc * ecc - 1.0)
        r0_pc = T9_R0_AU * AU_PC
        cosnu = max(-1.0, min(1.0, (p / r0_pc - 1.0) / ecc))
        nu = -math.acos(cosnu)
        r1, v1, r2, v2 = two_body_hyperbolic_inbound_state(
            m_co, m_star, -a_abs, ecc, nu, G_MSUN_PC_MYR
        )
    else:
        ecc = 1.0 - (r_peri_rsun * RSUN_PC) / semi_pc
        r1, v1, r2, v2 = two_body_apo_state(m_co, m_star, semi_pc, ecc, G_MSUN_PC_MYR)

    rows = [
        (m_co, r1, v1, kw_co),
        (m_star, r2, v2, kw_star),
        (0.5, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1),
    ]
    _assert_no_coincident_positions([(m, p, v) for m, p, v, _ in rows], f"t9_tde_trigger({kind},{beta})")
    write_rows_typed(output, rows)

    design = {
        "kind": kind,
        "beta_design": beta_design,
        "m_co": m_co,
        "m_star": m_star,
        "r_star_ref_rsun": T9_R_STAR_ZAMS,
        "r_tidal_design_rsun": float("nan") if kind == "msms" else r_t,
        "r_peri_design_rsun": r_peri_rsun,
        "r_init_au": (T9_R0_AU if kind == "hyperbolic" else semi_pc / AU_PC * (1.0 + ecc)),
        "t_pericentre_myr": math.pi * math.sqrt(semi_pc ** 3 / (G_MSUN_PC_MYR * (m_co + m_star))) if kind in ("bound", "msms", "ns") else float("nan"),
        "expect_tde": kind != "msms" and beta_design > 1.0,
        "metallicity": 0.02,
    }
    design_path = output.with_suffix("").with_name(output.name + ".design.json")
    design_path.write_text(json.dumps(design, indent=1), encoding="utf-8")
    return design


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


def _assert_no_coincident_positions(rows: List[Tuple[float, List[float], List[float]]], name: str) -> None:
    """Reject ICs with two particles at identical positions (trips SDAR's r2>0 assert)."""
    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            if sum((a - b) ** 2 for a, b in zip(rows[i][1], rows[j][1])) == 0.0:
                raise RuntimeError(f"{name}: rows {i} and {j} have identical positions")


def _generate_mcluster_plain_rows(n_star: int, seed: int, binary_fraction: float = 0.0) -> Tuple[List[Tuple[float, List[float], List[float]]], int]:
    """mcluster Plummer Rh=1 pc Kroupa rows in Msun / pc / km/s.

    With binary_fraction > 0 mcluster generates primordial binaries (ordinary
    hard binaries, NOT contact pairs — no stellar-evolution merge channel may
    rely on them) and writes their members first as adjacent (2k, 2k+1) pairs.
    Returns (rows, n_bin); the parsed binary count is verified against the
    expected value so the petar '-b' argument stays consistent.
    """
    n_bin_expected = int(binary_fraction * n_star // 2)
    with tempfile.TemporaryDirectory(prefix="petar_make_ic_") as tmpdir:
        tmp_prefix_name = "plainclu"
        cmd = [
            "mcluster",
            "-N", str(n_star),
            "-P", "0",
            "-R", "1",
            "-f", "1",
            "-b", str(binary_fraction),
            "-u", "1",
            "-C", "3",
            "-s", str(seed),
            "-o", tmp_prefix_name,
        ]
        proc = subprocess.run(
            cmd,
            cwd=tmpdir,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if binary_fraction > 0:
            match = re.search(r"Creating (\d+) primordial binary systems", proc.stdout)
            n_bin_found = int(match.group(1)) if match else -1
            if n_bin_found != n_bin_expected:
                raise RuntimeError(
                    f"mcluster primordial binary count mismatch: expected {n_bin_expected}, got {n_bin_found}"
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
                rows.append((float(parts[0]), [float(parts[1]), float(parts[2]), float(parts[3])], [float(parts[4]), float(parts[5]), float(parts[6])]))
    if len(rows) != n_star:
        raise RuntimeError(f"mcluster returned {len(rows)} stars, expected {n_star}")
    return rows, n_bin_expected


# 10 primordial binaries for N=100: exercises group/binary artifacts in the
# no-stellar-evolution modes without contact-binary pathologies.
FUNCTIONAL_SMOKE_BINARY_FRACTION = 0.2


def build_functional_mcluster_smoke(output: Path) -> None:
    """Plummer cluster with primordial binaries for modes without a stellar-evolution merge channel.

    10 ordinary hard binaries (mcluster Kroupa eigenperiod distribution; members
    written first as adjacent pairs) so `<snap>.binary` and group files carry
    real records in std/merger/dsm/gasdrag/external-only families.
    functional_mcluster_dual_merge_smoke instead injects ultra-compact massive
    binaries that BSE merges at the first interrupt; without stellar evolution
    those persist as contact binaries and stall the AR integrator in
    million-step loops.

    Output convention (same as functional_mcluster_dual_merge_smoke):
    position pc, velocity km/s -> petar.init -v KM_S_TO_PC_MYR under -u 1.
    """
    rows, _n_bin = _generate_mcluster_plain_rows(100, MCLUSTER_SEED, binary_fraction=FUNCTIONAL_SMOKE_BINARY_FRACTION)
    _assert_no_coincident_positions(rows, "functional_mcluster_smoke")
    write_rows(output, rows)


def build_functional_mcluster_pair_smoke(output: Path) -> None:
    """Binary-free cluster plus one inbound hyperbolic collision pair as the FIRST
    member pair (merger mode only; petar '-b 1').

    The simple-merger interrupt path hangs at t=0 when the collision pair is NOT
    pre-marked as a binary member pair (dynamic grouping + merger interrupt:
    silent infinite loop, 2026-10-17; dumps preserved in the functional merger
    case dir). With the pair as marked members (file rows 1-2, '-b 1') the
    interrupt fires through the member-pair path and completes — the
    configuration all earlier passing runs used implicitly via '-b 4'.
    """
    rows, _n_bin = _generate_mcluster_plain_rows(100, MCLUSTER_SEED)

    h_m1 = 10.0
    h_m2 = 10.0
    h_semi = -1.0e-8
    h_ecc = 1.001
    h_r1, h_v1, h_r2, h_v2 = two_body_hyperbolic_inbound_state(h_m1, h_m2, h_semi, h_ecc, -3.0955, G_MSUN_PC_MYR)
    if MCLUSTER_OUTPUT_VELOCITY_UNIT == "km/s":
        h_v1 = [v * PC_MYR_TO_KM_S for v in h_v1]
        h_v2 = [v * PC_MYR_TO_KM_S for v in h_v2]

    order = sorted(range(len(rows)), key=lambda i: -math.sqrt(sum(c * c for c in rows[i][1])))
    slot_a, slot_b = order[0], order[1]
    pcom, vcom = list(rows[slot_b][1]), list(rows[slot_b][2])
    pair_rows = [
        (h_m1, [pcom[k] + h_r1[k] for k in range(3)], [vcom[k] + h_v1[k] for k in range(3)]),
        (h_m2, [pcom[k] + h_r2[k] for k in range(3)], [vcom[k] + h_v2[k] for k in range(3)]),
    ]
    rest = [rows[i] for i in range(len(rows)) if i not in (slot_a, slot_b)]
    final_rows = pair_rows + rest

    _assert_no_coincident_positions(final_rows, "functional_mcluster_pair_smoke")
    write_rows(output, final_rows)


def build_t6(output: Path) -> None:
    """Plummer cluster with an engineered central massive triple for the T6
    perturber-kernel-match check.

    The compact triple (4+4 Msun inner binary, a=0.0012 pc, e=0.2; 3 Msun outer
    companion, a=0.008 pc, e=0.1) guarantees an AR group with artificial
    particles (tidal tensor engaged) from the first tree step, independent of
    chaotic realization differences; the run must pass --r-group/--r-search-
    group large enough to keep the outer orbit linked (e.g. 0.008/0.01 with
    --r-ratio 0.15, since r_search_group must stay <= r_in). The
    surrounding Plummer sphere (297 stars, scale a=0.3 pc, masses uniform in
    [0.2, 5] Msun) supplies list perturbers that cross the member-r_out..
    c.m.-r_out transition annulus of the triple group. Fixed seed. Units:
    Msun, pc, pc/Myr (G = G_MSUN_PC_MYR, same as t1-t3).
    """
    import random as _random

    rng = _random.Random(20261010)
    n = 297
    a = 0.3
    m_tot = 0.0
    rows = []
    for _ in range(n):
        mass = 0.2 + 4.8 * rng.random()
        m_tot += mass
        # Plummer radius by inverse-CDF sampling
        x1 = rng.random()
        while x1 < 1e-12:
            x1 = rng.random()
        r = a / math.sqrt(x1 ** (-2.0 / 3.0) - 1.0)
        r = min(r, 20.0 * a)
        ct = 2.0 * rng.random() - 1.0
        st = math.sqrt(1.0 - ct * ct)
        phi = 2.0 * math.pi * rng.random()
        pos = [r * st * math.cos(phi), r * st * math.sin(phi), r * ct]
        sigma = math.sqrt(G_MSUN_PC_MYR * m_tot / (6.0 * math.sqrt(r * r + a * a)))
        vel = [rng.gauss(0.0, sigma) for _ in range(3)]
        rows.append((mass, pos, vel))

    # engineered central massive triple: inner binary plus outer companion
    m1, m2, m3 = 4.0, 4.0, 3.0
    ra, va, rb, vb = two_body_peri_state(m1, m2, 0.0012 * (1.0 - 0.2), 0.2, G_MSUN_PC_MYR)
    rc, vc, rd, vd = two_body_peri_state(m1 + m2, m3, 0.008 * (1.0 - 0.1), 0.1, G_MSUN_PC_MYR)
    rows.append((m1, [rc[i] + ra[i] for i in range(3)], [vc[i] + va[i] for i in range(3)]))
    rows.append((m2, [rc[i] + rb[i] for i in range(3)], [vc[i] + vb[i] for i in range(3)]))
    rows.append((m3, rd, vd))

    cm_r = [sum(row[1][k] for row in rows) / len(rows) for k in range(3)]
    cm_v = [sum(row[2][k] for row in rows) / len(rows) for k in range(3)]
    rows = [
        (mass, [pos[k] - cm_r[k] for k in range(3)], [vel[k] - cm_v[k] for k in range(3)])
        for mass, pos, vel in rows
    ]
    write_rows(output, rows)


def build_functional_binary_cluster(output: Path) -> None:
    """Binary-rich Plummer cluster for the functional changeover-invariant case.

    mcluster N=200, Plummer Rh=0.5 pc, Kroupa IMF, 95% primordial binaries,
    fixed seed 11. mcluster writes binary members first as adjacent pairs,
    then remaining singles; the builder verifies the pair adjacency and the
    binary count so petar's '-b <n_bin>' (member IDs are 1..2*n_bin) holds.

    Output convention (same as functional_mcluster_dual_merge_smoke):
    position pc, velocity km/s -> petar.init -v KM_S_TO_PC_MYR under -u 1.
    """
    n_star = 200
    binary_fraction = 0.95
    n_bin_expected = int(binary_fraction * n_star // 2)

    with tempfile.TemporaryDirectory(prefix="petar_make_ic_") as tmpdir:
        tmp_prefix_name = "binclu"
        cmd = [
            "mcluster",
            "-N", str(n_star),
            "-P", "0",
            "-R", "0.5",
            "-f", "1",
            "-b", str(binary_fraction),
            "-u", "1",
            "-C", "3",
            "-s", "11",
            "-o", tmp_prefix_name,
        ]
        proc = subprocess.run(
            cmd,
            cwd=tmpdir,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        match = re.search(r"Creating (\d+) primordial binary systems", proc.stdout)
        if not match or int(match.group(1)) != n_bin_expected:
            found = match.group(1) if match else "none"
            raise RuntimeError(
                f"mcluster primordial binary count mismatch: expected {n_bin_expected}, got {found}"
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
                vel = [float(parts[4]), float(parts[5]), float(parts[6])]
                rows.append((mass, pos, vel))

    if len(rows) != n_star:
        raise RuntimeError(f"mcluster returned {len(rows)} stars, expected {n_star}")

    # mcluster writes members first as adjacent (2k, 2k+1) pairs; verify so the
    # petar '-b n_bin' ID convention (members are IDs 1..2*n_bin) is valid.
    pair_sep_max_sq = 0.05 ** 2
    for k in range(n_bin_expected):
        a, b = rows[2 * k], rows[2 * k + 1]
        d2 = sum((a[1][i] - b[1][i]) ** 2 for i in range(3))
        if d2 > pair_sep_max_sq:
            raise RuntimeError(
                f"binary pair {k} separation {math.sqrt(d2):.3e} pc exceeds 0.05 pc; "
                "member-first adjacency assumption broken"
            )

    write_rows(output, rows)


def build_functional_sevn_cluster(output: Path) -> None:
    """Single-star Plummer cluster for the functional SEVN case.

    Constraints from the SEVN stellar-evolution interface (see
    sample/star_cluster_plummer_N1k_binaries_sevn.sh):
     - stellar masses must stay inside the SEVNtracks_MIST_AGBrobust_gold grid
       (0.7-80 Msun); masses outside abort stellar initialisation;
     - no primordial binaries: binary-rich SEVN ICs can trigger a
       dump_large_step_h4 abort at t=0.

    mcluster N=60, Plummer Rh=1 pc, Kroupa IMF 0.7-80 Msun, no binaries,
    fixed seed 7.

    Output convention (same as functional_mcluster_dual_merge_smoke):
    position pc, velocity km/s -> petar.init -v KM_S_TO_PC_MYR under -u 1.
    """
    n_star = 60

    with tempfile.TemporaryDirectory(prefix="petar_make_ic_") as tmpdir:
        tmp_prefix_name = "sevnclu"
        cmd = [
            "mcluster",
            "-N", str(n_star),
            "-P", "0",
            "-R", "1",
            "-f", "1",
            "-b", "0",
            "-m", "0.7",
            "-m", "80",
            "-u", "1",
            "-C", "3",
            "-s", "7",
            "-o", tmp_prefix_name,
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
                vel = [float(parts[4]), float(parts[5]), float(parts[6])]
                rows.append((mass, pos, vel))

    if len(rows) != n_star:
        raise RuntimeError(f"mcluster returned {len(rows)} stars, expected {n_star}")
    bad = [m for m, _, _ in rows if m < 0.7 or m > 80.0]
    if bad:
        raise RuntimeError(f"{len(bad)} star masses outside the SEVN gold-table grid 0.7-80 Msun (e.g. {bad[0]:.3f})")

    write_rows(output, rows)


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
            "t5",
            "t6",
            "t7",
            "t7_band2",
            "t7_control_in",
            "t7_outer",
            "t7_abortzone",
            "t7_abortzone2",
            "t8_cascade",
            "t8_chaotic",
            "t8_control",
            "t9_tde_trigger",
            "t9_replay",
            "t9_replay_fix",
            "functional_smoke",
            "functional_dual_merge_smoke",
            "functional_mcluster_dual_merge_smoke",
            "functional_mcluster_smoke",
            "functional_mcluster_pair_smoke",
            "functional_binary_cluster",
            "functional_sevn_cluster",
        ],
        required=True,
    )
    parser.add_argument("--output", default="", help="output IC path (not needed for t9_replay_fix)")
    parser.add_argument(
        "--kind",
        choices=["bound", "hyperbolic", "ns", "msms", "pc_control", "pc_perturbed"],
        default="bound",
        help="T9 tde_trigger IC variant",
    )
    parser.add_argument("--beta", type=float, default=2.93, help="T9 tde_trigger design beta")
    parser.add_argument("--event", type=int, default=9, help="T9 replay event id (micro_TDE reference)")
    parser.add_argument("--input-file", default="", help="T9 replay_fix: converted petar input file to patch")
    parser.add_argument("--time", type=float, default=0.0, help="T9 replay_fix: snapshot epoch [Myr]")
    args = parser.parse_args()

    if args.case != "t9_replay_fix" and not args.output:
        parser.error(f"--output is required for case {args.case}")

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
    elif args.case == "t5":
        build_t5(out)
    elif args.case == "t6":
        build_t6(out)
    elif args.case == "t7":
        build_t7(out, rp=0.15)
    elif args.case == "t7_band2":
        build_t7(out, rp=0.30)
    elif args.case == "t7_control_in":
        build_t7(out, rp=0.05)
    elif args.case == "t7_outer":
        build_t7(out, rp=0.60)
    elif args.case == "t7_abortzone":
        build_t7(out, rp=0.35)
    elif args.case == "t7_abortzone2":
        build_t7(out, rp=0.50)
    elif args.case == "t8_cascade":
        build_t8_cascade(out)
    elif args.case == "t8_chaotic":
        build_t8_chaotic(out)
    elif args.case == "t8_control":
        build_t8_control(out)
    elif args.case == "t9_tde_trigger":
        build_t9_tde_trigger(out, kind=args.kind, beta=args.beta)
    elif args.case == "t9_replay":
        build_t9_tde_event_replay(out, event_id=args.event)
    elif args.case == "t9_replay_fix":
        if not args.input_file:
            raise ValueError("t9_replay_fix requires --input-file and --time")
        fix_t9_replay_input(Path(args.input_file), t_snap=args.time)
    elif args.case == "functional_smoke":
        build_functional_smoke(out)
    elif args.case == "functional_dual_merge_smoke":
        build_functional_dual_merge_smoke(out)
    elif args.case == "functional_mcluster_dual_merge_smoke":
        build_functional_mcluster_dual_merge_smoke(out)
    elif args.case == "functional_mcluster_smoke":
        build_functional_mcluster_smoke(out)
    elif args.case == "functional_mcluster_pair_smoke":
        build_functional_mcluster_pair_smoke(out)
    elif args.case == "functional_binary_cluster":
        build_functional_binary_cluster(out)
    elif args.case == "functional_sevn_cluster":
        build_functional_sevn_cluster(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())