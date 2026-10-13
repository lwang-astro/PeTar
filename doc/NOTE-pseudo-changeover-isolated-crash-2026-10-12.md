# Handoff note: isolated-cluster crash from pseudo-particle c.m.-changeover switch — ATTRIBUTION CORRECTED (not a 3848bc8 regression)

Date: 2026-10-12 (original) · corrected 2026-10-12 (late) · Found by: validation-suite T3 rework (`t3_binary_sdar_grouplink`) ·
Status: **RESOLVED 2026-10-12 (late): root cause found and fixed (SDAR `block_time_step.h` + `hermite_integrator.h`, uncommitted); NOT a 3848bc8 regression; reproducer + eta variant run clean to t=3.2 at dE/E~1.3e-12**

## 0. CORRECTION (read first)

The original note attributed the crash to commit 3848bc8 ("pseudo particles take the group c.m. changeover").
**That attribution is wrong.** Controlled A/B tests on 2026-10-12 evening:

| experiment (all with identical input/params, forced full rebuilds) | result |
|---|---|
| HEAD + revert of the 3848bc8 pseudo-branch (12 lines) | **still crashes**, identical signature t=0.168945 |
| HEAD + SDAR rolled back d474412 (6020bfc excluded; `touch src/*.hpp; make` — the Makefile does NOT track `../SDAR/src` headers, so earlier rollback tests without forced rebuild were linking stale objects) | **still crashes**, identical signature |
| checkout 980e6ad (the claimed working commit) + SDAR d474412 | **still crashes**, identical signature |
| `--tt-switch 0` (tidal tensor off) | **still crashes**, identical signature and counters |
| `-O1` build / `-O3 -g` build / instrumented builds | **runs clean** — the trigger is marginal (FP/codegen sensitive) |

**Decisive observation:** `/tmp/t3_diag/A.log` — the run this note originally cited as
"same config on the pre-regression binary — completed cleanly" — **itself crashed with the
identical error** (A.log line 2347, `A_grouped_plain.hard_dump_t0.168945_*`, same abort stack).
There was never a clean baseline among the saved artifacts; only the never-grouped variant
(B.log, Hermite pair) completed. The "AR error 4e-12 over 113 orbits" clean run is not among
the saved logs. Bisecting across differently-built installed binaries conflated codegen/layout
differences with source changes.

## 1. Problem (unchanged)

A plain-family PeTar run whose whole system is one isolated hard cluster containing a single
2-member AR group aborts with `time step size too small (2.22e-16) < dt_min (4.44e-16)`
inside the hard integration.

Reproducer (unchanged):

```bash
python3 test/validation/make_ic.py --case t2 --output /tmp/repro/input.base
cd /tmp/repro && petar.init -f input input.base
OMP_NUM_THREADS=1 <plain-family petar> -a 0 --r-group 0.018 --r-search-group 0.02 \
  -u 1 -t 3.2 -s 9.765625e-4 -r 0.2 -o 1.953125e-3 -w 2 -f repro input
```

## 2. What the crash regime actually is

Per-step log counters (`Stab/Iso/Form/End_Ngroup`, hard log "Number per step") show, for
EVERY tree step of the whole run:

- `Stab_Ngroup = 0`: the find level (`findGroupsAndCreateArtificialParticlesOMP`) **never**
  links the pair into a group — consistent with the mass-weighted criterion lesson of
  2026-10-12 (effective threshold for the 1:10 binary ≈0.013-0.018 regime is marginal);
- `Form_Ngroup = 1, End_Ngroup = 0`: each tree-step segment the two singles enter `h4_int`,
  the pair is linked **inside the integrator** (`adjustGroups`→`addGroups`), integrated with
  AR (~15-23 steps), and dissolves again at the next find level.

So the whole run lives in the "never-find-grouped / always-in-integrator re-grouped" regime.
At t≈0.168945 (≈6 orbital periods, near pericenter) the first Hermite correction of a segment
aborts: `calcBlockDt4th` ← `correctAndCalcDt4thList` ← `integrateSingleOneStepAct`
(gdb backtrace, frame #2-#4). The dump written at abort is the segment-start backup
(2 singles, n_group=0 — clean pos/vel).

Characterization experiments:

- `--hermite-eta-init 0.0001` (10x smaller initial steps): crash **moves** to t=0.227051 —
  the divergence is marginal and timing shifts with integrator parameters;
  this is genuine numerics, not memory UB.
- Debug replayer (`petar.mpi.omp.avx2.hard.debug`, -O0 + full asserts) replays the crash
  segment **cleanly** — again consistent with a marginal trigger.

## 2b. ROOT CAUSE (resolved)

The abort is NOT a force/numerical divergence at all. gdb evidence (exact-entry breakpoints
on `BlockTimeStep4th::calcBlockDt4th`; the aborting call):

- all acceleration inputs are EXACTLY ZERO (`acc0=acc1=acc2=acc3=0`);
- `dt_limit == dt_min` (4.44e-16), `dt_ref = 4.207e-16 < dt_min` → one halving below the
  grid floor → abort printing `dt = 2.22e-16 = dt_min/2`.

With zero force, `calcDt4th` takes the `HERMITE_MONOPOLE_DT_LIMIT` branch (PeTar Makefile
enables it): `dt_ref = eta * |pos| / |vel|` — the flying-time bound added for the 2026-09-13
IMBH slingshot. The aborting particle is the **group c.m. of the in-integrator formed pair**:
the isolated cluster is integrated in its own c.m. frame, so the group c.m. sits at
roundoff-distance from the origin with roundoff velocity — `|pos|/|vel|` is a pure
cancellation artifact that can land anywhere, here just below the grid floor. This is why the
crash flipped with every codegen change (-O1/-O3/-g/instrumentation) and moved with
`--hermite-eta-init`: the artifact depends on ulp-level rounding of the c.m. state.

A second, compounding defect in SDAR `initialIntegration`: the wide-binary guard
`while(dt_limit > t_peri) dt_limit *= 0.5;` has no floor. At formations with a
near-parabolic/hyperbolic instantaneous Kepler estimate the boundary mean anomaly is ~0,
`t_peri ~ 1e-16`, and the loop crushes `info.dt_limit` from dt_max to the grid floor
(observed: dt_limit = dt_min at the abort). Alone it is survivable (zero-force early return
usually saves the step); combined with the monopole artifact it guarantees sub-dt_min steps.

**Fixes (SDAR, uncommitted, at 6020bfc):**

1. `block_time_step.h` `calcBlockDt2nd/4th`: floor `dt_ref` at `dt_min` when the jerk is
   exactly zero (the monopole branch) — a roundoff `|pos|/|vel|` is not a dynamical time;
   force-driven small `dt_ref` still aborts (policy preserved).
2. `hermite_integrator.h` `initialIntegration`: floor the t_peri halving target at
   `step.getDtMin()`.

**Validation:** reproducer runs to t=3.2 (113 orbits) clean, final |dE/E| ~ 1.3e-12
(the AR-grouped accuracy the note originally expected); `--hermite-eta-init 1e-4` variant
(previously crashed at t=0.227) also clean to t=3.2.

## 3. Real state-pollution bug found and fixed (hygiene, uncommitted)

`ParticleH4` fields `acc0/acc1/pot/time/dt` are **never initialized in the soft system**
(they are IC garbage; observed NaN). They flow:

- soft → hard: `setPtclFor[Isolated]MultiCluster[OMP]` `DataCopy` brings NaN into
  `ptcl_hard_` every segment (this is the "garbage doubles after the first few fields" seen
  in dumps: the 9 columns after `r_out` are exactly acc0[3]/acc1[3]/pot/time/dt);
- hard → soft: `writeBackPtclForOneCluster[OMP]` `DataCopy` copies the AR member state
  (untouched Hermite fields) back into `system_soft`.

The singles path re-initializes these fields (`initialIntegration`), which is why the garbage
is normally invisible — but it is undefined behavior feeding every non-resetting consumer
(a plausible source of the codegen sensitivity). `r_search` of ex-members also carries the
stale group-era value (0.252) back into the soft system.

**Fix applied in `src/hard.hpp` (+14 lines, the two soft→hard entry copies
`setPtclFor[Isolated]MultiCluster[OMP]`, uncommitted, HEAD=449718a):** zero
`acc0/acc1/pot/time/dt` on the entry copy (`ptcl_hard_` is `PtclH4`, which has the fields).
Exit-side zeroing in `writeBackPtclForOneCluster` was tried first but does not compile
generically (`Tsys` can be `FPSoft`, which has no `acc0`) and was dropped. **The crash
itself is NOT fixed by this** — verified with a compiling build, crash time unchanged
(t=0.168945) — it is a correctness/hygiene fix, kept for review. An earlier "fix test"
in this session was invalid (build errors left the old binary in place); the result above
is from the verified build.

## 4. Remaining hunt for the crash driver (suggested next steps)

The abort is a Hermite single-particle correction with astronomically large acc/jerk
derivatives developing within the first few H4 steps of a segment, in the
singles→in-integrator-group→(break?) cycle near pericenter. Suggested entry points:

1. Instrument `correctAndCalcDt4thOne` inputs (acc0..acc3, dt_ref) at the abort site
   (SDAR `block_time_step.h`) in a build that still crashes — note any instrumentation
   shifts codegen and can mask the crash; iterate flags until a crashing instrumented
   build is found.
2. Check the in-integrator formation near pericenter: `checkNewGroup`/`addGroups` initial
   dt and the first AR/H4 synchronization after formation (`integrateToTimeList`,
   `calcDt2ndList` with zero-neighbour group c.m. force + `acc0_offset_sq`).
3. Compare against the structurally identical find-level-grouped regime (force with a
   larger `--r-ratio` so r_in grows past the effective mass-weighted threshold;
   `r_search_group <= r_in` is asserted) — if the sym_int/TT path is clean for the same
   orbit, the bug is specific to the in-integrator formation/dissolution cycle.
4. The B variant (never grouped, pure Hermite pair) completes 113 orbits cleanly at 4.7e-7 —
   use it as the reference trajectory to localize where the crashing run's singles state
   starts to deviate.

## 5. What is NOT the cause (updated)

- 3848bc8 pseudo-branch change — excluded by direct revert A/B (§0).
- SDAR 6020bfc kappa guard — excluded with a forced-rebuild rollback (§0).
- Tidal tensor mode — excluded (`--tt-switch 0` identical crash).
- The unified group criterion (SDAR 5a62e4e), BSE, mass-weighted linking — unchanged from
  original §5 (and the mass-weighted criterion lesson stands).
- The NaN acc/time/dt flow — real bug, fixed (§3), but not the crash driver.

## 6. Blocked work (unchanged, plus caveat)

`test/validation/t3_binary_sdar_grouplink` (T3) remains blocked on this crash. NOTE: the T3
grouped regime (r_group 0.014/0.018) IS the crashing regime (§2); the earlier census lesson
(Form=1 via in-integrator formation, not find level) applies. Physics findings (a)/(b) of the
original note stand.

## 7. Session state

- PeTar at `experiment` HEAD 449718a; only `src/hard.hpp` modified (+26 hygiene lines, §3);
  pre-existing unrelated modifications in validation docs/metrics untouched.
- SDAR restored to 6020bfc, clean (the temporary DT_DIAG patch in `block_time_step.h` was
  reverted).
- Build directory rebuilt at HEAD + fix (crash reproducible, same t=0.168945).
- `/tmp/r1/` holds this session's runs (repro/fix/ei/ttoff logs+par); `/tmp/t3_diag/` as before.
