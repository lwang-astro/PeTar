#!/usr/bin/env python3
import argparse
import math
import random
from pathlib import Path
from typing import List, Tuple


G_MSUN_PC_MYR = 0.00449830997959438


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


def build_t4(output: Path) -> None:
    m1 = 1.0
    m2 = 0.5
    m3 = 0.05

    rin_peri = 0.002
    ein = 0.6
    r1_rel, v1_rel, r2_rel, v2_rel = two_body_peri_state(m1, m2, rin_peri, ein, G_MSUN_PC_MYR)

    m12 = m1 + m2
    rout_peri = 0.04
    eout = 0.5
    r12, v12, r3, v3 = two_body_peri_state(m12, m3, rout_peri, eout, G_MSUN_PC_MYR)

    p1 = [r12[i] + r1_rel[i] for i in range(3)]
    p2 = [r12[i] + r2_rel[i] for i in range(3)]
    v1 = [v12[i] + v1_rel[i] for i in range(3)]
    v2 = [v12[i] + v2_rel[i] for i in range(3)]

    rows: List[Tuple[float, List[float], List[float]]] = [(m1, p1, v1), (m2, p2, v2), (m3, r3, v3)]

    rng = random.Random(20260413)
    n_bg = 24
    m_bg = 0.02
    rmax = 0.25
    vmax = 0.25
    for _ in range(n_bg):
        radius = rmax * (rng.random() ** (1.0 / 3.0))
        cost = 2.0 * rng.random() - 1.0
        sint = math.sqrt(max(0.0, 1.0 - cost * cost))
        phi = 2.0 * math.pi * rng.random()
        pos = [
            radius * sint * math.cos(phi),
            radius * sint * math.sin(phi),
            radius * cost,
        ]

        v_radius = vmax * rng.random()
        v_cost = 2.0 * rng.random() - 1.0
        v_sint = math.sqrt(max(0.0, 1.0 - v_cost * v_cost))
        v_phi = 2.0 * math.pi * rng.random()
        vel = [
            v_radius * v_sint * math.cos(v_phi),
            v_radius * v_sint * math.sin(v_phi),
            v_radius * v_cost,
        ]
        rows.append((m_bg, pos, vel))

    mass_tot = sum(item[0] for item in rows)
    com_pos = [sum(item[0] * item[1][k] for item in rows) / mass_tot for k in range(3)]
    com_vel = [sum(item[0] * item[2][k] for item in rows) / mass_tot for k in range(3)]

    recentered: List[Tuple[float, List[float], List[float]]] = []
    for mass, pos, vel in rows:
        pos_new = [pos[k] - com_pos[k] for k in range(3)]
        vel_new = [vel[k] - com_vel[k] for k in range(3)]
        recentered.append((mass, pos_new, vel_new))

    write_rows(output, recentered)


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate deterministic initial conditions for PeTar validation scenarios")
    parser.add_argument("--case", choices=["t1", "t2", "t3", "t4"], required=True)
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())