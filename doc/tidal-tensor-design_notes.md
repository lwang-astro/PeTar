# Tidal-tensor (TT) method — design principles and code mapping

Status: reference note — documents the design rationale, algorithm, and
paper-vs-code mapping of the artificial-particle tidal-tensor (TT) method
(Wang et al. 2020, MNRAS 497, 536; `doc/PeTar_code.pdf` §6.4, pp. 543–547) as
implemented in PeTar `src/`. Compiled 2026-10-27 as background for the TT
inner-binary degradation investigation (now closed; the temporary plan was
deleted after distillation — its record lives in git history as
`doc/plan-tt-tree-channel-inner-binary-degradation.md`); records mechanism
and design invariants only — failure analysis and hypotheses live in git
history and lessons 2026-10-27.

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
- **Considered and rejected (2026-10-27): a debug-gated leak-detector assert in
  `fit()`** — `|T3|·dc² ≤ C·|T2|·dc` (C~10). Principle: T2 fit rows sum to zero
  (immune to a constant offset in the stencil accs) while T3 rows do not (a
  common offset leaks in with gain ~1/dc²), so the ratio |T3|dc/|T2| ~ dc/R for
  physical fields but explodes under any constant-offset pollution — it would
  have caught F13 at the first fit. Rejected because of false-positive risk in
  production-relevant regimes: when the physical |T2| is near zero (no nearby
  perturber, near-uniform external field) the ratio is ill-conditioned and the
  assert cannot distinguish pollution from a legitimately negligible tensor.
  If a constant-offset corruption class recurs, re-evaluate with the ratio
  guarded by a |T2| floor. The input invariant above stays the authoritative
  statement.
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
- **Activation structure (corrected twice 2026-10-28 — polarity verified against
  code + fine-output census)**: TT has TWO entry paths:
  1. **General path (the default)**: any stable binary subtree reaching the
     tail loop of `findGroupsAndCreateArtificialParticlesOneCluster`
     (`src/hard.hpp`, `createArtificialParticles` call) gets TT with NO
     tide-strength/κ/tt-nstep/cost condition. Gates to reach it: tt-switch
     on; members in one hard cluster; stability (`findStableTree`,
     `t_crit = cap` — binaries with `period > cap` are not stable-grouped
     at all); and not caught by the isolated opt-out below. In multi-particle
     clusters (≥3 particles) the opt-out CANNOT fire — every stable binary
     sharing its cluster with any companion gets TT unconditionally.
  2. **Isolated branch = an OPT-OUT for pair-only clusters** (`_n_ptcl==2`):
     a trial slowdown (timescale = 2·cap) gives κ_est; if
     `period·κ_est > dt_tree·tt-nstep` the cap-limited slowdown already
     resolves the slowed orbit ("it is fine without tidal tensor") → plain
     slowdown group, NO artificial particles; otherwise fall through to the
     general path → TT. NOTE the polarity: the printed decision boolean
     `P·κ > cap` TRUE = plain (no TT), FALSE = TT. `pert_out` uses the c.m.'s
     OWN r_search (r_out-scaled + velocity cushion; T4-m5: ≈0.11 pc vs the
     actual perturber at 0.149) — order-unity accurate (T4: hyp 243 vs real
     ~302 per binary mass).
  Empirically (audits 2026-10-28) the opt-out NEVER fired (all recorded
  booleans 0 → TT): the mcluster ultra-hard binary, the T4-m5 binary, and a
  sparse-cluster population (26 live TT groups, n_all 200→536 measured with
  fine outputs) all take the general path — TT-on is effectively the default
  for stable binaries, and the compute-cost question lives on this path
  (empirical extension in §6: cost 0–30% at CPU N≲1e3 scales).
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
- **Label-swap symmetry property (verified 2026-10-29, fixed same day)**: under
  a swap of the binary member order, `particleToOrbit` gives incline/rot_horizon
  invariant (am = r×v with both factors negated) but `rot_self →
  rot_self + π` exactly — the PM triangle is point-reflected in the orbital
  plane. The design target (quadrupole, all even moments) is invariant;
  the SPURIOUS odd multipoles of the 3-point discretization (the difference
  between 3 points and the continuous ellipse) flip sign. Net effect on a
  test particle: a discontinuous force jump of relative size ~(a_bin/d)³
  (equal-mass triangle ⇒ zero dipole; leading spurious odd moment is the
  octupole; measured 3.4e-8 at d/a=150, 1.4e-5 at d/a=20, 9.2e-4 at d/a=5 —
  a constant ×(a/d)³ fit to 3%). The
  member order within a pair comes from the binding-chain reordering in
  `generateBinaryTree` (not a stable convention), so the mirror can fire at
  group reformations. **Fix (2026-10-29)**: canonical orientation in
  `PseudoParticleMultipoleManager::createSampleParticles` — the triangle is
  point-reflected (equivalent to rot_self += π, quadrupole unchanged) unless
  member 1 is the heavier, with the member IDs breaking the exact
  equal-mass tie; the rule is a total canonical ordering on leaf members, so
  both member orders map to the identical particle set. Residual corner: at
  group FORMATION the stable-tree nodes carry no node ids (only the
  AR-info/Hermite refresh trees get them via `checkAndSetBinaryPairID`), so
  an exactly-equal-AGGREGATE-mass hierarchical pair falls through the tie
  break (id 0 = 0) and may point-reflect once at the first refresh —
  measure-zero outside exactly-equal-mass ICs. Unit test
  `src/pm_multipole_test.cxx` (`make build/petar.pm.test`): 5000 random
  bound binaries built through the production path
  (`generateBinaryTree` → `createSampleParticles`), asserting
  PM(p1,p2) ≡ PM(p2,p1) — max relative position difference 2.8e-16
  post-fix (roundoff) vs 1.74 pre-fix (the point-reflection itself).
- **Counterforce method ladder (available escalation, 2026-10-29 audit)**:
  the counterforce representation of a binary for external particles is
  selectable at build/run time, in increasing fidelity and cost:
  (1) c.m.-only — switch SDAR off (`--r-group 0`); (2) PM, 3
  pseudo-particles — the default (`typedef PseudoParticleMultipoleManager
  OrbitManager`, `src/artificial_particles.hpp`); (3) orbit sampling, N
  particles — compile-time `-D ORBIT_SAMPLING` + `--os-nsplit`
  (`src/orbit_sampling.hpp`, `OrbitalSamplingManager`); paper Fig. 6:
  OS-8 ≈ PM ≪ OS-4 ≪ c.m.-only for momentum conservation, none affect the
  orbital evolution. Replacing PM by the **instantaneous member components**
  ("改回双星本体") was considered and rejected: the external force would
  carry the phase-random quadrupole fluctuation, relative amplitude
  ~(a_bin/d)² per kick, aliased at the unsynchronized tree cadence;
  PM's residual (static spurious odd moments, see above) is one to two
  multipole orders smaller and frozen with the orbit frame. The correct
  escalation when PM resolution is insufficient is OS-8 or enlarging the
  changeover `r_out` (cluster capture), not exposing the members.
- **Group-radius chain (PM validity region is environment-adaptive)**: the
  annulus where PM represents the binary (perturber inside `r_search_group`,
  outside the changeover `r_out`) is not a fixed multiple of `r_in`.
  Chain (`HardManager::initial`, `src/hard.hpp`): the 1D velocity dispersion
  σ₁ = √(⟨v²⟩/3) is measured from the ICs (`src/petar.hpp`, passed into
  `initial`); the hard-soft boundary a_hs = G·⟨m⟩/σ₁² then gives
  `r_search_group`(auto) = min(r_in, safety·a_hs) (`--r-search-group-safety`,
  default 1.0) and `r_group`(auto) = 0.8·r_search_group, with assertions
  enforcing 0 ≤ r_group ≤ r_search_group ≤ r_in. Consequence: in cold
  (bound, clustered) systems the auto chain shrinks the PM annulus to the
  physical hard-soft boundary, keeping d/a_bin large where PM is used;
  the r_in-proportional fallback (no σ₁) is the legacy behavior. Restart
  caveat: `--r-search-group-safety` or changeover updates re-determine both
  radii unless they were explicitly given.
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

**Empirical extension (2026-10-28/29 campaign; raw evidence in
`test/validation/work/{ttgrid,ttaudit,ttcost}/`, distilled from the closed
activation-criterion plan)**:

- Accuracy surface (48-run grid, T4 triple geometry, relative tide
  2.1e-8 … 5.7e-3, dt_soft/P_in ∈ {0.26, 1.05, 4.21}): TT matches the
  hard-channel reference across the entire range (inner-Δa ratio ≤ 1.14),
  including dt/P = 4.2 (4 inner orbits per tree kick) and the strongest
  tide (perturber peri 0.0495 pc). The no-TT error is a tide-INDEPENDENT
  kick-aliasing floor that dominates weakward (33–957x the physical signal
  at tide ≤ 2.1e-4, dt/P ≥ 1); at dt/P = 0.26 aliasing is negligible and
  no-TT ≈ reference. Consequence: "weak tide → TT useless" is refuted on
  accuracy — weak tide is where no-TT is worst relative to the signal;
  TT's own floor is 100–1000x lower. TT over-activation exists only as
  compute waste in the (fine step × weak tide) corner.
- Strong-tide boundary: clean to relative tide 5.7e-3; tide ~0.21
  (m1000/a0.15 probe) is a chaotic disruption/exchange regime where
  grouping churns and drift metrics are meaningless — the practical upper
  boundary. Retained caveat from the inclined-orbit probe (tide 2.1e-2,
  incl 60°): secular-angle agreement has a phase-drift gap (ψ 3.9% total,
  inclination pointwise 21.4° as a phase-shifted oscillation) — consistent
  with the paper's own "Except ψin" caveat; element-amplitude agreement
  remains good.
- Coverage hole (perturber strong but outside r_search with small κ_est)
  could NOT be produced: raising the mean mass 33–333x still activated TT
  (census N_all=15) with TT/ref ≤ 1.0 in Δe — the criterion is
  conservative on this grid.
- Cost (plain 1860e binary, 3 reps, OMP=4): N=200/Rh=1pc (dense,
  hard-channel-dominated) TT +~30% wall; N=1000/Rh=2pc statistically
  indistinguishable. Attribution: Tree_Force +156% with TT — the analytic
  particle multiplier 1+6.5·f_bin from live TT groups (n_all 1000→2548;
  12 artificial particles per TT binary in the tree); PP_single (hard)
  −9% (slowdown pays back part of the tree overhead). Untested
  extrapolation: tree-dominated regimes (large-N with GPU offload) where
  the 2.5x Tree_Force particle term's wall share grows — judged
  low-priority (the activation-criterion research closed 2026-10-29 with
  no code change; revisit only if a GPU production profile shows tree
  dominance).

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

- Closed investigation (temporary plan, deleted after distillation): see git
  history `doc/plan-tt-tree-channel-inner-binary-degradation.md` and lessons
  2026-10-27 entries
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
