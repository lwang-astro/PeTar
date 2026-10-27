# Tidal-tensor (TT) method — design principles and code mapping

Status: reference note — documents the design rationale, algorithm, and
paper-vs-code mapping of the artificial-particle tidal-tensor (TT) method
(Wang et al. 2020, MNRAS 497, 536; `doc/PeTar_code.pdf` §6.4, pp. 543–547) as
implemented in PeTar `src/`. Compiled 2026-10-27 as background for the open
investigation `plan-tt-tree-channel-inner-binary-degradation.md` (T4
inner-binary degradation); records mechanism and design invariants only —
failure analysis and hypotheses live in that plan.

## 1. The problem TT solves (Hamiltonian-split context)

PeTar splits `H = H_S + H_L` (paper §5–6): short-range interactions are
integrated exactly (Hermite + SDAR/AR, individual steps) inside "hard"
clusters; long-range interactions are computed by the particle tree and
applied as kicks on a step `dt_L`. The changeover function `W(r)` provides
the smooth transition between the two regimes.

The long-range force is therefore **piecewise constant in time**: it is
evaluated once per tree kick and applied as an impulse. For a binary with
period `T_bin ≪ dt_L`, that impulse samples a **random orbital phase**, so
the differential (tidal) part of the kick on the two components is wrong at
the O(tide) level and accumulates secularly. The counterforce on distant
perturbers has the mirror problem: they would see the binary at one random
phase instead of its orbit average.

The slow-down (SD) method lengthens the effective binary period by the
factor `κ` and so helps make `dt_L ≲ κ·T_bin`, but only weakly perturbed
binaries can reach a large `κ`; a tight binary with neighbours keeps
`T_bin,eff ≪ dt_L` (frequent in clusters).

The artificial-particle algorithm (paper §6.4) resolves both sides:

- **tide on the binary** — reconstruct the smooth *spatial* structure of the
  long-range field around the binary c.m. (the tidal tensor), and evaluate
  it at each member's instantaneous position at **every AR step** during the
  hard segment, instead of applying one phase-random kick;
- **counterforce on perturbers** — represent the binary by orbit-sampling or
  pseudo-particle multipole particles so distant particles feel the
  orbit-averaged field (§5 below).

Cost trade-off chosen by design: more (zero/small-mass) particles in the
tree — cheap to parallelize — instead of forcing small `dt_L`.

## 2. Tensor construction (measure points + fit)

3D Taylor expansion of the long-range acceleration about the binary c.m.
(paper eq. 25):

```
A_i(r) = A_i(rcm) + A_ij·r'_j + ½ r'_j·A_ijk·r'_k + …     r' = r − rcm
A_ij = ∂A_i/∂r_j          A_ijk = ∂²A_i/∂r_j∂r_k
```

For gravity `A_ij`, `A_ijk` are symmetric: 2nd-order method stores 6 (+3)
elements, 3rd-order 10 (+3+6). Measurement recipe (paper §6.4.1):

- zero-mass **measure particles** are placed near the binary c.m. and are
  ordinary members of the soft (tree) system — they receive the long-range
  force at each tree kick;
- `A(rcm)` is measured directly from the c.m. acceleration;
- `A'(r_k) = A(r_k) − A(rcm)` at each point gives 3 linear equations on the
  tensor elements; the linear map `M·T = A'` (paper eq. 26–28) is inverted
  once and for all — the fit is a fixed linear combination of the measured
  accelerations;
- geometry: **2nd order** = 4 points at the corners of a regular tetrahedron
  (condition number 2); **3rd order** = 8 points at the corners of a cube
  (the 6-point octahedron is rank-deficient; cube condition number 12.7).
  Paper Table 1 gives the exact point coordinates and fit formulae.

Code mapping:

- `TidalTensor` (`src/tidal_tensor.hpp`): stores `T1[3]` (kept zero — the
  c.m. part is subtracted), general `T2[9]` (symmetrized on fit), symmetric
  `T3[10]` when `TIDAL_TENSOR_3RD` is defined (it is, in all standard build
  families: `Makefile.in` MT_FLAGS; see `src/build_info.hpp`).
- `TidalTensor::fit()` hard-codes the `M⁻¹` rows (numerically identical to
  paper Table 1); `createTidalTensorMeasureParticles()` places 3rd-order
  points at cube corners with half-edge `dc = 0.16·r_tidal_tensor`
  (2nd-order tetrahedron scale `dc = 0.5·size`); the fit rescale is `1/dc`.
- `r_tidal_tensor` is set to `r_group` (`HardManager::initial`,
  `src/hard.hpp`): the tensor scale and measure-point extent track the group
  radius, not the binary semi-major axis.
- **Input invariant of `fit()`**: it assumes the measure accs are already
  c.m.-subtracted (T1≡0 by construction; the Axxx-style T3 fit rows have nonzero
  coefficient sums, so a residual common offset leaks into T3 with gain ~1/dc²,
  while T2 rows sum to zero and stay immune). The subtraction
  (`ArtificialParticleManager::correctArtficialParticleForce` →
  `TidalTensor::subtractCMForce`) belongs to the force-evaluation phase only:
  under `KDKDK_4TH` the mid-step gradient pass (`treeSoftGradient`, tree reuse)
  re-enters the same correction path and must not subtract again — the accs are
  refreshed only once per tree step, and the 8-point stencil cannot absorb a
  constant offset (free-T1 fit is rank-deficient: constant ↔ even-T3 degeneracy),
  so any gauge residue is amplified by 1/dc². The guard (`!_acorr_flag` around
  the subtraction call sites in `src/hard.hpp`) enforces exactly-once semantics
  (plan F13/F15, fixed 2026-10-27).
- Measure particles co-move with the c.m. velocity and carry zero mass, so
  they only ever sample the long-range force.

## 3. Application inside AR

- `TidalTensor::eval()` / `evalPot()`: evaluates the Taylor series at a
  c.m.-frame position; `evalPot` returns the linear-gauge potential
  `−r·A(r)`.
- In the AR interaction (`src/ar_interaction.hpp`, perturbation-force
  evaluation used for both `n_pert>0` and isolated branches): at **every AR
  perturbation evaluation**, `soft_pert->eval()` is applied to each member
  with `|r| < r_out` (changeover gate) — **with its member mean kept intact**:
  the mean of a correct tensor is the physical second-order c.m. coupling of
  the extended binary and drives the outer-orbit secular evolution (measured
  2026-10-27: the outer-orbit drift Δe_out ≈ 3.3e-5 over 20 Myr is physical —
  three independent integrators agree to 3 digits with e_out(t) correlation
  +1.000; subtracting the mean suppresses it ~100x with phase correlation
  −0.63). History: a member-mean subtraction was introduced at 1795e to cure
  a "momentum leak" (binary c.m. drift) that was in fact the mean of the
  F13-corrupted tensor (~8000x the physical mean); after the F13 producer
  fix (acorr-phase guard) the subtraction was reverted 2026-10-27 — momentum
  cost of the revert is negligible (T4: 2.6e-7 vs 2.3e-7 kept; T5 gate 5e-3
  passes). See the (closed) investigation plan and lessons 2026-10-27 for
  the full evidence chain.
- Slowdown coupling: `ARPerturber::calcSoftPertMin()` (`src/ar_perturber.hpp`)
  evaluates the tensor's differential acceleration at the apocentre
  configuration (`E = π`) and stores `soft_pert_min`, which is added to
  `pert_out` in the slowdown-`κ` criterion (`evalPertOuterForce` tail). The
  tensor therefore also *measures* the perturbation strength that limits the
  slow-down factor.
- Activation semantics: `--tt-nstep` sets `n_step_per_orbit`;
  `ar_manager.slowdown_timescale_max = dt_max_hermite · n_step_per_orbit`
  caps the isolated-binary slowdown timescale. A binary whose SD factor
  would need `κ·T_bin` above that cap is instead handled with the TT
  representation. See lessons 2026-10-27 (tt-nstep is a cap/threshold, not a
  tensor update frequency).

## 4. Freeze/measurement cadence — the central design invariant

Paper, §6.4.1: *"In the P3T method, the numerical long-range forces are
constant within one dt_L, so should be the tensor field. Thus, we only need
to measure the tensors once per dt_L."*

The tensor is **frozen in time over exactly one long-range-force interval**
— the same interval for which the tree kick it reconstructs is frozen. The
method is the smooth-in-space reconstruction of a piecewise-constant-in-time
force; freeze interval (tensor) ≡ freeze interval (tree force) is the design
invariant that makes the reconstruction consistent with the Hamiltonian
split.

Code cadence (`integrateOneCluster`, `src/hard.hpp`): per `hard()` call —

1. `HardIntegrator::initial()` — fits every group tensor **once**, from the
   measure particles' just-computed tree accelerations; binds
   `group.perturber.soft_pert` and `tidal_tensor_index[]` (member → tensor
   mapping); computes `soft_pert_min`.
2. `integrateToTime(dt)` — integrates the whole hard segment (interval
   `dt` = the tree/block step for that cluster); **no refit inside the
   loop**; when `adjustGroups` reforms a group, the `soft_pert` pointer is
   only re-bound if the new membership is consistent with the recorded
   tensor.
3. `clear()` — tensors discarded; measure particles are re-created around
   the current c.m. for the next segment by `updateArtificialParticles()`
   (`src/artificial_particles.hpp`; called from the group-update paths in
   `hard.hpp`) — the position refresh happens by re-measurement next kick,
   not by shifting the stored tensor.

`TidalTensor::shiftCM()` (Taylor re-centering of a frozen tensor) exists,
but the per-step c.m. update inside `integrateToTime` is disabled with the
comment "Not consistent, should not shift" — the stored `pos` remains the
fit-time c.m. for the whole segment.

So: **one fit per hard segment** = per tree-kick interval, exactly the
paper's "once per dt_L" — provided the hard segment and the tree-force
freeze interval are the same interval (see §7.1).

## 5. Counterforce side (perturbers receiving the binary field)

- **Orbit-sampling** (paper §6.4.2.1): sample particles along both
  component orbits at equal eccentric-anomaly interval `ΔE`, masses
  weighted by the mean-anomaly interval (eq. 29) so the particles represent
  the time-averaged mass distribution; `src/orbit_sampling.hpp`.
- **Pseudo-particle multipole** (Kawai & Makino 2001): 3 equal-mass
  pseudo-particles reproduce the orbit-averaged traceless quadrupole
  `⟨A'⟩` of the binary (paper eq. 30–38; positions eq. 38);
  `src/pseudoparticle_multipole.hpp`.
- Paper result (Fig. 6): the counterforce choice does **not** affect the
  orbital evolution (within resolution), only the quality of linear
  momentum conservation — OS-8 ≈ PM (3 pseudo-particles) ≪ OS-4 ≪ c.m.-only.

## 6. Paper validation conditions and retained limitations

What the published claim covers (paper §6.4.3, Table 2, Fig. 5):

- hierarchical triple, perturber **outside** the binary's changeover region;
  inner binary integrated by SDAR, outer orbit by leapfrog with the same
  `dt_L`;
- `dt_L = 3.90625e-3` vs `T_bin,in = 1.97e-3`, i.e. `dt_L ≈ 2·T_bin,in`
  (no-TT samples the phase about every 2 inner periods);
- TT-3 (+ pseudo-particle multipole) reproduces the SDAR reference's
  **secular** evolution of the inner *and* outer orbital elements; TT-2
  reproduces the inner orbit but fails the outer orbit's secular
  oscillation; no-TT fails both.

Retained limitations stated in the paper (§6.4.1 end):

- the tensor **cannot be constructed mid-drift** — measure points receive
  forces only at kicks, so a binary newly formed inside a drift interval
  has no tensor until the next kick;
- the Taylor field is **invalid once members separate beyond the expansion
  radius** (disrupted binary);
- "a short interval error within one dt_L is not very serious" — the design
  goal is the correct **long-term cumulative (secular)** effect, not
  per-interval exactness.

Not covered by the published validation: block-stepped tree channels
(perturber and binary c.m. on different step hierarchies), freeze windows
spanning many inner periods under a bound slowdown factor, and eccentricity
conservation metrics at the precision T4 now demands.

## 7. Design invariants → code anchors (verification checklist)

Neutral mapping for investigators; the open questions are in the plan (§3
Q1–Q5).

1. **Freeze-interval equality** (tensor frozen exactly as long as the tree
   force it reconstructs). Anchor: one `fit()` per `HardIntegrator::initial()`
   per hard segment. Verify: hard-segment duration vs the *perturber's* own
   kick/block-step interval — under block stepping the two need not
   coincide (the tree force on the measure points updates on the measure
   points' own kick schedule; plan Q1/Q2).
2. **Smooth-in-space, frozen-in-time field.** Anchor: `eval()` at
   instantaneous member positions each AR evaluation; coefficients constant
   through the segment. Verify: measured tensor vs analytic
   quadrupole+octupole of the perturber state at the fit instant
   (sign/normalization; plan step 3, `PERT_KERNEL_DIAG`-style dump).
3. **Differential-only application.** Anchor: member-mean subtraction in
   `ar_interaction.hpp`; `T1 = 0` by fit. Verify: gauge consistency between
   `acc_pert` and the linear-gauge `pot_pert`, and against the c.m. kick the
   tree delivers to the group c.m.
4. **Validity radius.** Anchor: measure points at `0.16·r_group`; Taylor
   expansion assumes perturber distance ≫ group extent. Verify: `a_out`
   vs `r_group` for the scenario in question.
5. **Perturber outside changeover.** Anchor: `eval` gated by
   `|r| < r_out`; tree-side perturber reach bounded by `r_search_group`.
   Verify: whether the scenario's third body is inside/outside each radius
   (plan F2/Q4 — without TT, tree particles beyond `r_search_group` never
   perturb the AR binary at all).
6. **Activation threshold semantics.** Anchor:
   `slowdown_timescale_max = dt_max_hermite · tt-nstep`; TT activates only
   when the isolated-binary SD factor would exceed the cap (plan F3, F5
   saturation; lessons 2026-10-27).

## 8. Pointers

- Investigation plan: `plan-tt-tree-channel-inner-binary-degradation.md`
- Paper: `PeTar_code.pdf` §6.4 (pp. 543–547) — §6.4.1 tensor construction,
  Table 1 measure-point coordinates + fit formulae, §6.4.2 counterforce,
  Table 2/Fig. 5 validation, Fig. 6 counterforce momentum; §6.2 SDAR
  slowdown; §5 P3T split.
- Code: `src/tidal_tensor.hpp` (TidalTensor), `src/artificial_particles.hpp`
  (ArtificialParticleManager), `src/ar_perturber.hpp` (ARPerturber),
  `src/ar_interaction.hpp` (application + slowdown link), `src/hard.hpp`
  (HardIntegrator::initial / integrateToTime, integrateOneCluster),
  `src/orbit_sampling.hpp`, `src/pseudoparticle_multipole.hpp`.
- Unit-test hook: `src/tidal_tensor_test.cxx`.
