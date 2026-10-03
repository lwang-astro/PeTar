//
// Created by Daniel Marín Pina
//
#pragma once

#include <cmath>
#include <stdexcept>
#include <algorithm>

/**
 * Evolve a compact binary under GW radiation reaction (Peters 1964).
 *
 * Integrates da/dt and de/dt until either dt is elapsed or the
 * sticky-sphere merger condition a*(1-e) <= r1+r2 is met.
 *
 * @param[in,out] a: Semimajor axis [Rsun]
 * @param[in,out] e: Eccentricity
 * @param m1, m2:  Component masses [Msun]
 * @param r1, r2:  Stellar radii   [Rsun]
 * @param dt: Maximum integration time [Myr]
 * @param rtol: Relative tolerance per step (default 1e-10)
 * @return Elapsed time [Myr]; equals dt if no merger, or the estimated merger time if merger occurred.
 */
inline double evolve_bco(double& a, double& e, double m1, double m2, double r1, double r2, double dt,
                         double rtol = 1.0e-10) {
    // Physical constants
    constexpr double G = 3.947841760435743e13; // 4*pi**2*1e12
    constexpr double c = 6.3198e10;
    constexpr double c5 = c * c * c * c * c;
    constexpr double G3 = G * G * G;

    // Input validation
    if (m1 <= 0.0 || m2 <= 0.0) throw std::invalid_argument("masses must be positive");
    if (r1 < 0.0 || r2 < 0.0) throw std::invalid_argument("radii must be non-negative");
    if (a <= 0.0) throw std::invalid_argument("semimajor axis must be positive");
    if (e < 0.0 || e >= 1.0) throw std::invalid_argument("eccentricity must be in [0, 1)");
    if (dt < 0.0) throw std::invalid_argument("dt must be non-negative");

    const double r_sum = r1 + r2;

    // Already merged at t = 0
    if (a * (1.0 - e) <= r_sum) return 0.0;

    const double beta = G3 * m1 * m2 * (m1 + m2) / c5;

    //  Peters RHS as a local lambda
    auto rhs = [&](double av, double ev, double& dadt, double& dedt) {
        const double e2 = ev * ev;
        const double oMe2 = 1.0 - e2;
        dadt = -(64.0 / 5.0) * beta
            * (1.0 + (73.0 / 24.0) * e2 + (37.0 / 96.0) * e2 * e2)
            / (av * av * av * std::pow(oMe2, 3.5));
        dedt = -(304.0 / 15.0) * beta * ev
            * (1.0 + (121.0 / 304.0) * e2)
            / (av * av * av * av * std::pow(oMe2, 2.5));
    };

    // Single RK4 step
    auto rk4 = [&](double av, double ev, double h,
                   double& a_out, double& e_out) {
        double da1, de1, da2, de2, da3, de3, da4, de4;
        auto clamp_e = [](double x) { return std::max(0.0, std::min(x, 1.0 - 1e-12)); };
        rhs(av, clamp_e(ev), da1, de1);
        rhs(av + 0.5 * h * da1, clamp_e(ev + 0.5 * h * de1), da2, de2);
        rhs(av + 0.5 * h * da2, clamp_e(ev + 0.5 * h * de2), da3, de3);
        rhs(av + h * da3, clamp_e(ev + h * de3), da4, de4);
        a_out = av + (h / 6.0) * (da1 + 2 * da2 + 2 * da3 + da4);
        e_out = ev + (h / 6.0) * (de1 + 2 * de2 + 2 * de3 + de4);
    };

    // Auto initial step: small fraction of the Peters coalescence time
    const double T_merge_approx = (12.0 / 85.0) * a * a * a * a / beta;
    double h = std::min({T_merge_approx * 1.0e-4, dt * 1.0e-3, dt});
    h = std::max(h, 1.0); // floor at 1 s to avoid stalling

    // Adaptive RK4 integration (step-doubling error control)
    constexpr double safety = 0.9;
    constexpr double h_min = 1.0e-3; // [s]
    constexpr double grow_max = 5.0;
    constexpr double shrink_min = 0.1;
    constexpr long max_steps = 10'000'000;

    double t = 0.0;
    long nsteps = 0;

    while (t < dt) {
        if (nsteps++ > max_steps)
            throw std::runtime_error("evolve_binary: max_steps exceeded");

        h = std::min(h, dt - t); // don't overshoot end time
        if (h < h_min) h = h_min;

        // One full step and two half-steps (step-doubling)
        double a1, e1;
        rk4(a, e, h, a1, e1); // one step of h
        double am, em;
        rk4(a, e, h * 0.5, am, em); // first half-step
        double a2, e2;
        rk4(am, em, h * 0.5, a2, e2); // second half-step

        // Richardson error estimate
        const double err_a = std::abs(a2 - a1) / 15.0;
        const double err_e = std::abs(e2 - e1) / 15.0;

        // Mixed relative/absolute error norm
        const double scale_a = rtol * std::max(std::abs(a), std::abs(a2));
        const double scale_e = std::max(rtol * std::max(std::abs(e), std::abs(e2)), 1.0e-13);
        const double err = std::max(err_a / scale_a, err_e / scale_e);

        if (err > 1.0 && h > h_min) {
            // Reject step — shrink h and retry
            h *= std::max(shrink_min, safety * std::pow(err, -0.25));
            --nsteps;
            continue;
        }

        // Richardson-extrapolated best estimate (5th-order)
        a2 += (a2 - a1) / 15.0;
        e2 += (e2 - e1) / 15.0;
        e2 = std::max(0.0, std::min(e2, 1.0 - 1.0e-12));

        // Check for merger after this accepted step
        if (a2 * (1.0 - e2) <= r_sum) {
            // Bisect within [0, h] to locate the merger time precisely
            double t_lo = 0.0, t_hi = h;
            double a_hi = a2, e_hi = e2;
            for (int i = 0; i < 60; ++i) {
                double t_mid = 0.5 * (t_lo + t_hi);
                double a_mid, e_mid;
                rk4(a, e, t_mid, a_mid, e_mid);
                if (a_mid * (1.0 - e_mid) <= r_sum) {
                    t_hi = t_mid;
                    a_hi = a_mid;
                    e_hi = e_mid;
                } else {
                    t_lo = t_mid;
                }
                if (t_hi - t_lo < 1.0) break; // converged to ~ 1s
            }
            a = a_hi;
            e = e_hi;
            return t + t_hi; // merger time
        }

        // Accept step
        a = a2;
        e = e2;
        t += h;

        // Grow h for next step
        h *= (err < 1.0e-30)
                 ? grow_max
                 : std::min(grow_max, safety * std::pow(err, -0.25));
    }

    return t; // full dt elapsed, no merger
}
