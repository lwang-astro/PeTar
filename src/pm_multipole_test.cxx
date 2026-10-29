//! Unit test: pseudo-particle multipole (PM) placement invariance under binary
//! member-order swap.
//!
//! Background (2026-10-29): particleToOrbit stores the orbit under the
//! (member1, member2) convention of the binding chain; under a member-order
//! swap rot_self shifts by pi and the PM triangle point-reflects. The
//! quadrupole (design target) is invariant, but the spurious odd multipoles
//! of the 3-point discretization flip sign, giving every external particle a
//! discontinuous force jump (~(a/d)^3 relative) at each group reformation.
//! The canonical-orientation fix (see pseudoparticle_multipole.hpp) defines
//! the in-plane +x axis toward the orbit-averaged position of the heavier
//! member (masses decide; member IDs break exact ties), making the placement
//! label-swap invariant by construction.
//!
//! This test drives the REAL production path: two particle states ->
//! COMM::BinaryTree::particleToOrbit -> PseudoParticleMultipoleManager::
//! createSampleParticles, for random bound orbits with random orientations,
//! and asserts the two orderings give identical PM positions.
//!
//! Build (see Makefile target build/petar.pm.test):
//!   mpic++ -I./src -I../FDPS/src -I../SDAR/src -std=c++17 ... pm_multipole_test.cxx

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>

#include "particle_simulator.hpp"
#include "hard_ptcl.hpp"
#include "pseudoparticle_multipole.hpp"

namespace {

//! minimal particle for the test (mass/pos/vel/id as used by particleToOrbit)
struct TestPtcl {
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::S64 id;
};

using TestBinTree = COMM::BinaryTree<TestPtcl, COMM::Binary>;

TestPtcl makeParticle(const PS::F64 mass, const PS::F64vec& pos, const PS::F64vec& vel, const PS::S64 id) {
    TestPtcl p;
    p.mass = mass;
    p.pos = pos;
    p.vel = vel;
    p.id = id;
    return p;
}

//! generate a bound binary with random orientation: relative orbit with given
//! semi/ecc at a random true anomaly, then random 3D rotation of the pair
void makeBinary(std::mt19937_64& rng, TestPtcl& p1, TestPtcl& p2) {
    std::uniform_real_distribution<PS::F64> u01(0.0, 1.0);
    const PS::F64 G = 4.498502151469553e-3; // pc^3/(Msun Myr^2)
    // component masses: asymmetric in general (the common case; the
    // canonical rule uses masses), occasionally equal (ID tie-break branch)
    PS::F64 m1 = 0.1 + u01(rng)*9.9;
    PS::F64 m2 = (u01(rng) < 0.2) ? m1 : 0.1 + u01(rng)*9.9;
    const PS::F64 m12 = m1 + m2;
    const PS::F64 GM = G*m12;
    const PS::F64 semi = 0.002 + u01(rng)*0.006;
    const PS::F64 ecc = 0.05 + u01(rng)*0.85;
    // place at a random true anomaly
    const PS::F64 f = u01(rng)*2.0*M_PI;
    const PS::F64 cf = std::cos(f), sf = std::sin(f);
    // relative state in the perifocal frame (true-anomaly form):
    // r = p/(1+e cos f), v = sqrt(GM/p) (-sin f, e+cos f)
    const PS::F64 p = semi*(1.0-ecc*ecc);
    const PS::F64 den = 1.0 + ecc*cf;
    PS::F64vec rel_r(p*cf/den, p*sf/den, 0.0);
    const PS::F64 h = std::sqrt(GM*p);
    PS::F64vec rel_v(-GM/h*sf, GM/h*(ecc + cf), 0.0);
    // random 3D rotation via two Euler angles
    const PS::F64 a1 = u01(rng)*2.0*M_PI, a2 = u01(rng)*2.0*M_PI, c2 = std::cos(a2), s2 = std::sin(a2);
    auto rot = [&](PS::F64vec& v) {
        PS::F64vec w(v.x, v.y*c2 - v.z*s2, v.y*s2 + v.z*c2);
        const PS::F64 c1 = std::cos(a1), s1 = std::sin(a1);
        v = PS::F64vec(w.x*c1 - w.y*s1, w.x*s1 + w.y*c1, w.z);
    };
    rot(rel_r);
    rot(rel_v);
    // c.m. at origin with zero velocity
    p1 = makeParticle(m1, -m2/m12*rel_r, -m2/m12*rel_v, 101);
    p2 = makeParticle(m2,  m1/m12*rel_r,  m1/m12*rel_v, 202);
}

PS::F64 pmPositions(const TestPtcl& pa, const TestPtcl& pb, PS::F64vec out[3]) {
    // production-faithful tree construction: generateBinaryTree builds the
    // binding chain, sets member pointers, aggregate fields (mass/pos/vel)
    // and the orbit elements (via particleToOrbit) exactly as in a group
    // formation; createSampleParticles then runs on a fully-built node.
    TestPtcl ptcl[2] = {pa, pb};
    int ptcl_list[2] = {0, 1};
    TestBinTree bins[1];
    COMM::BinaryTree<TestPtcl, COMM::Binary>::generateBinaryTree(bins, ptcl_list, 2, ptcl, 4.498502151469553e-3);
    PseudoParticleMultipoleManager mgr;
    TestPtcl arti[3];
    mgr.createSampleParticles(arti, bins[0]);
    PS::F64 scale = 0.0;
    for (int k=0; k<3; k++) {
        out[k] = arti[k].pos;
        // NaN/Inf guard: a vacuous NaN==NaN comparison would silently pass
        if (!std::isfinite(out[k].x) || !std::isfinite(out[k].y) || !std::isfinite(out[k].z)) {
            std::printf("FAIL: non-finite PM position at particle %d\n", k);
            std::exit(2);
        }
        scale = std::max(scale, std::abs(out[k].x) + std::abs(out[k].y) + std::abs(out[k].z));
    }
    return scale;
}

} // namespace

int main() {
    std::mt19937_64 rng(20261029);
    const int n_trial = 5000;
    PS::F64 max_rel_diff = 0.0;
    int worst = -1;
    for (int i=0; i<n_trial; i++) {
        TestPtcl p1, p2;
        makeBinary(rng, p1, p2);
        PS::F64vec pmA[3], pmB[3];
        const PS::F64 scaleA = pmPositions(p1, p2, pmA);
        const PS::F64 scaleB = pmPositions(p2, p1, pmB);
        const PS::F64 scale = std::max(scaleA, scaleB);
        for (int k=0; k<3; k++) {
            const PS::F64vec d = pmA[k] - pmB[k];
            const PS::F64 rel = std::sqrt(d*d)/std::max(scale, 1e-300);
            if (rel > max_rel_diff) {
                max_rel_diff = rel;
                worst = i;
            }
        }
    }
    std::printf("PM member-order swap invariance: n_trial=%d  max relative |dpos| = %.3e (worst trial %d)\n",
                n_trial, max_rel_diff, worst);
    // bound: canonical orientation must give roundoff-level agreement.
    // (pre-fix code fails this at O(1): the triangle point-reflects, giving
    //  relative differences of order the point radius itself)
    const PS::F64 tol = 1e-12;
    if (max_rel_diff > tol) {
        std::printf("FAIL: max relative difference %.3e exceeds tolerance %.3e\n", max_rel_diff, tol);
        return 1;
    }
    std::printf("PASS (tolerance %.1e)\n", tol);
    return 0;
}
