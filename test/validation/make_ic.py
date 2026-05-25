#!/usr/bin/env python3
import argparse
import math
import subprocess
import tempfile
from pathlib import Path
from typing import List, Tuple


G_MSUN_PC_MYR = 0.00449830997959438
KM_S_TO_PC_MYR = 1.022712165045695


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
    """Generate a physical background with mcluster and inject two target binaries.

    Steps:
    1) Use mcluster to generate N=100, Plummer, Rh=1 pc, Kroupa IMF, no binaries.
    2) Replace the first four stars by two binaries (4 stars total).
    3) Align each binary center of mass with one of the first two 2-star pair COMs.
    4) Recenter to zero net COM position/velocity in Msun/pc/pc-Myr units.
    """
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
            "42",
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
                vel = [
                    float(parts[4]) * KM_S_TO_PC_MYR,
                    float(parts[5]) * KM_S_TO_PC_MYR,
                    float(parts[6]) * KM_S_TO_PC_MYR,
                ]
                rows.append((mass, pos, vel))

    if len(rows) < 4:
        raise RuntimeError("mcluster returned fewer than 4 stars; cannot inject two binaries")

    # Two selected binaries from sample outcomes.
    b1_m1 = 85.294690901635846
    b1_m2 = 23.528452122895846
    b1_period = 4.594336e-09
    b1_ecc = 0.499774915487773

    b2_m1 = 35.947378474297203
    b2_m2 = 48.425671915280390
    b2_period = 8.3652155e-08
    b2_ecc = 0.427779218694198

    def binary_rel_state(m1: float, m2: float, period_myr: float, ecc: float) -> Tuple[List[float], List[float], List[float], List[float]]:
        semi = (G_MSUN_PC_MYR * (m1 + m2) * (period_myr / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0)
        peri = semi * (1.0 - ecc)
        return two_body_peri_state(m1, m2, peri, ecc, G_MSUN_PC_MYR)

    b1_r1, b1_v1, b1_r2, b1_v2 = binary_rel_state(b1_m1, b1_m2, b1_period, b1_ecc)
    b2_r1, b2_v1, b2_r2, b2_v2 = binary_rel_state(b2_m1, b2_m2, b2_period, b2_ecc)

    # Anchor binary COMs to COMs of the first two star pairs: (0,1) and (2,3).
    def pair_com(row_a: Tuple[float, List[float], List[float]], row_b: Tuple[float, List[float], List[float]]) -> Tuple[List[float], List[float]]:
        ma, pa, va = row_a
        mb, pb, vb = row_b
        mt = ma + mb
        pcom = [(ma * pa[k] + mb * pb[k]) / mt for k in range(3)]
        vcom = [(ma * va[k] + mb * vb[k]) / mt for k in range(3)]
        return pcom, vcom

    pcom1, vcom1 = pair_com(rows[0], rows[1])
    pcom2, vcom2 = pair_com(rows[2], rows[3])

    rows[0] = (b1_m1, [pcom1[k] + b1_r1[k] for k in range(3)], [vcom1[k] + b1_v1[k] for k in range(3)])
    rows[1] = (b1_m2, [pcom1[k] + b1_r2[k] for k in range(3)], [vcom1[k] + b1_v2[k] for k in range(3)])
    rows[2] = (b2_m1, [pcom2[k] + b2_r1[k] for k in range(3)], [vcom2[k] + b2_v1[k] for k in range(3)])
    rows[3] = (b2_m2, [pcom2[k] + b2_r2[k] for k in range(3)], [vcom2[k] + b2_v2[k] for k in range(3)])

    # Global recenter.
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

    write_rows(output, recentered)


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