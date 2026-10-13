# Handoff note: isolated-cluster crash from pseudo-particle c.m.-changeover switch — ATTRIBUTION CORRECTED (not a 3848bc8 regression)

Date: 2026-10-12 (original) · corrected 2026-10-12 (late) · Found by: validation-suite T3 rework (`t3_binary_sdar_grouplink`) ·
Status: **RESOLVED 2026-10-12 (late): root cause found and fixed (SDAR `block_time_step.h` + `hermite_integrator.h`, uncommitted); NOT a 3848bc8 regression; reproducer + eta variant run clean to t=3.2 at dE/E~1.3e-12** ·
**RESOLVED (2026-10-13 late): §9b — the T2 "explosion" was a build-state artifact of one specific binary, NOT an SDAR regression; clean paired A/B (incl. a 6020bfc shim endpoint) and a full T2 rerun reproduce the desktop baseline exactly**

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

## 8. Follow-up (2026-10-13, laptop LWLaptop): T1 coarsest-dt full-crossing endpoint crashes

After the 5765150+SDAR-63e8308 fix (T3 now passes 10/10 with bit-identical values on both
machines), a remaining failure appears on the laptop (Intel, gcc 13.3, Ubuntu 24.04):
T1 `crossfull_r015_si07` (coarsest dt_soft = 1.5625e-2, full crossing r_out=0.15, both 64b and
non64b, avx2 and avx512) aborts with `time step size too small (2^-49 or 2^-48) < dt_min
(2^-47)` near the first pericenter (t~0.40). The demanded dt sits exactly on the block-step
grid floor. si00-si06 of the same family pass; the desktop (AMD, older gcc) completed this
config with a large-but-stable error (3.4e-5).

Verified: deterministic on the laptop; reproduces with clean builds of SDAR 459e/460e/461e
(provenance checked via the configure-time version stamp; see build-hygiene notes below).
Conclusion: this endpoint sits at a stability cliff of the Hermite block-step criterion where
compiler/CPU FP details decide pass/fail; it is not an SDAR-version regression in 459-461e.

Build-hygiene traps discovered (documented in lessons): (1) the Makefile does not track
../SDAR header changes - SDAR checkouts require `touch src/hard.hpp` or a full clean; (2)
the version stamp `_SDARver` is baked at configure time into src/get_version.hpp - only a
fresh ./configure re-stamps it, so `make clean` alone keeps a stale stamp and misattributes
provenance; (3) always show build output - several silent partial builds invalidated early
bisect attempts.

## 9. NEW (2026-10-13, laptop): T2 coarse-dt mix-regime error explosion — ~4 orders above desktop baseline

After the 2c3dd3e gated fix resolved the §8 endpoint crashes (T1 64b now 6/6 PASS, mix
follows Hermite exactly, ratio 1.0 at si00), the T2 rerun shows the coarse half of the dt
sweep exploded while the fine half is bit-consistent with the desktop:

| run | regime (auto radii) | desktop (SDAR 6020bfc) | laptop (2c3dd3e binary) |
|---|---|---|---|
| dt17-dt12 | tree / fine mix | 6e-8 … 6e-6 | **matches** (6.01e-8 … 6.10e-6) |
| dt11 | mix, r_in=4.245e-3 (peri=4e-3 crosses it) | 2.8e-7 | **7.98e-3** |
| dt10 | mix | 4.7e-7 | **1.58e-2** |
| dt9 | mix | 6.3e-7 | **2.41e-3** |
| dt8 | pure Hermite (m1 side still in annulus) | 7.5e-7 | **6.51e-3** |

The "error suddenly small late in the run" curve first reported for this figure is **dt11**:
its per-period envelope OSCILLATES orbit-to-orbit (drops 2700x at orbit 60, jumps back 2652x
at orbit 63, again at 109/110) — each changeover-annulus crossing injects a ~1e-2-level error
that swings back next orbit; the late-time "small" values are envelope sampling on individual
orbits, not a real error decrease (quarters decline smoothly 8.0e-3 → 2.9e-4 for the same
reason). No AR group ever forms in these runs (Form=0/End=0); the pair is integrated purely
by Hermite/tree with the annulus crossing between them.

Same-machine single-variable A/B on dt11 (PeTar ec1dc7a fixed, only SDAR source varied,
fresh configure+build each time):

| SDAR | dt11 max abs error |
|---|---|
| 63e8308 | aborts (known crash era, §8) |
| b51ffa6 | 7.98e-3 |
| 2c3dd3e (gated fix) | 7.98e-3, **bit-identical to b51ffa6** |

Conclusion for handoff:
- The gated fix is NOT the cause — the gate only fires on the degenerate (zero acc/jerk)
  criterion branch, which this run never takes; b51ffa6 and 2c3dd3e are bit-identical here.
- The regression entered in the 63e8308..b51ffa6 chain and SURVIVES the gate. Prime suspect:
  63e8308's revert of the in-estimator monopole crossing bound — the annulus-crossing
  regime is exactly what that estimator path exercised. Needs SDAR-side review.
- Caveat: the 1e-7 desktop baseline ran SDAR 6020bfc against older PeTar (API-incompatible
  with current HEAD), so it is cross-machine evidence, not a same-machine A/B endpoint.
- T2's current checks pass despite the explosion (only rho-stationarity for dt17/dt12/dt8
  plus a dt17/dt12 level ratio are checked); per-regime level bounds should be added once
  the physics verdict is in.

Reproducer (dt11, ~40 s single-thread on the laptop):

```bash
python3 test/validation/make_ic.py --case t2 --output input.base
petar.init -f input input.base
OMP_NUM_THREADS=1 petar.mpi.omp.avx512 -a 0 --r-group 5e-4 --r-search-group 1e-3 \
  -u 1 -b 1 -t 3.2 -s 4.8828125e-4 -o 9.765625e-4 -w 2 -f dt11 input
# desktop-era level: ~2.8e-7; current level: ~7.98e-3 (max |Error/Total|)
```

## 9b. RESOLUTION (2026-10-13 late): §9 was a build-state artifact, not an SDAR regression

Same-machine paired A/B (PeTar ec1dc7a tree fixed, only the SDAR checkout varied,
forced rebuilds `touch src/hard.hpp; make -j`, binary taken from `build/`):

| SDAR | dt11 max \|err/total\| |
|---|---|
| 6020bfc + no-op `setRDtLimitOne` shim (the missing endpoint; note current PeTar builds without `HERMITE_MONOPOLE_DT_LIMIT`, so the old estimator branch was inactive) | **2.748237e-07** |
| 2c3dd3e (experiment) | **2.748237e-07 — identical** |

Provenance corroboration: a fresh forced rebuild of SDAR 2c3dd3e under the current configure
state is md5-identical (`1cd98362...`) to the binary installed at 16:56, and independently
reproduces 2.748237e-07; the 6020bfc-shim and 2c3dd3e builds give the same value, so the
level is robust to the SDAR source difference under one build configuration.

Full T2 rerun with the current build: dt11/dt10/dt9/dt8 = 2.75e-7 / 4.69e-7 / 6.25e-7 / 7.54e-7 —
**matching the desktop baseline column to the digit**; the fine half unchanged. The "explosion"
is not reproducible from the current tree.

What produced the 7.98e-3: the binary installed at the time of the §9 runs (call it B1, since
overwritten at 16:56 by a reinstall). B1 was built in the bisect-heavy window where -O1 /
-O3-g / instrumented configure states were being cycled (§0 documents these); autoconf keeps
previous arguments in `config.status`, so a plain `./configure` re-run can silently inherit
non-default flags — B1's exact configure line is unrecoverable, but it is the only constant
shared by every 7.98e-3 observation (the §9 "b51ffa6 bit-identical to 2c3dd3e" A/B ran the
same B1-style build for both, which is why source variation showed nothing).

Why build flags can flip 4 orders here: §9's own envelope observation — each annulus crossing
injects a ~1e-2-level error that largely CANCELS the next orbit (2700x drops/swings) — means
the observable error is a small residual of large cancelling injections. Such marginal
cancellation is exactly where FMA/contraction codegen differences decide whether the residual
stays 1e-7 or sits at the 1e-3..1e-2 injection level. Consistent with §8's cliff sensitivity
and with §0's -O1/-O3 behavior flips.

Hygiene rules going forward: (1) A/B endpoints must record BOTH binaries' md5 + configure-time
version stamp before trusting a level difference; (2) after bisect/instrumented sessions,
`./configure` does not reset flags — explicitly re-configure with the intended options before
producing numbers; (3) level comparisons in cancelling-error regimes are build-sensitive —
prefer same-binary A/B or report both binaries' provenance.
