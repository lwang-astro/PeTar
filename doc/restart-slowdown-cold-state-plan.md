# Plan: AR ds-controller collapse on rate-type H error — the "restart slow" grind (cause found)

Status: plan — open items only. Found 2026-10-19 by the functional `bse` case
(6.3 min vs seconds). **Root cause identified 2026-10-19 by `petar.hard.debug`
dump replay + gdb sampling** (see "Root cause"); the original restart-cold-state
hypothesis is refuted. Remaining items: fix decision + harness stabilization.

## Symptom

Restarting from an intermediate snapshot makes the hard integration of
pre-existing tight binaries pathological, while the identical evolution
integrated continuously is cheap:

| | fresh run 0→10 | restart run (from data.3 → 10) |
|---|---|---|
| wall time | **1.4 s** | **5 m 38 s** (restart segment only) |
| `Warning: step count is signficiant large` (10⁶-step caps) | 0 | **20** |
| 121+124 Msun contact-binary merge trajectory | t≈3 → t≈6 | **identical** (t≈3 → t≈6) |
| final state / all functional assertions | pass | pass |

Physics, BSE evolution, and merge times are the same along both paths — only
the cost differs. The restarted process grinds tight-binary encounters with
repeated step-count-max caps.

## Root cause (confirmed 2026-10-19)

Two independent findings; the first is the pathology, the second explains why
it was first misattributed to restart.

### 1. The grind: SDAR AR ds-controller collapses under rate-type H error

Measured by replaying `data.dump_large_step_ar_n3_t8.171875_*` with
`petar.mpi.omp.avx512.bse.hard.debug` (bit-identical reproduction: 5.54e6 AR
substeps in one cluster step) and gdb sampling inside
`SDAR/src/AR/symplectic_integrator.h` `integrateToTime`:

- **Configuration**: ultra-tight BSE binary (m≈23.66+19.04, semi≈8.3e-7,
  ecc≈5e-4, P≈1.1e-8) + slow near-parabolic third body (e≈1.0012, a≈−1.66,
  v_rel≈v_esc) lingering at r≈0.036 ≈ 4×10⁴·a_in for ~1 time unit. The trio is
  a 3-member AR group inside an n=80 Hermite cluster.
- **Slowdown machinery works**: κ_in ≈ 1.3e6, pinned at its timescale cap
  (κ = ts_in/P_in), tracking the encounter timescale. NOT a cold slowdown.
- **But κ_in drifts** at relative rate ~50–100 per unit time (the intruder's
  orbital phase drives `ts_in`), and the sd-energy bookkeeping follows with a
  ~2 %-per-update residual (gdb two-sample diffs: `ekin_sd_`/`epot_sd_`/
  `etot_sd_ref_` all drift together, `etot_sd_ref_` lagging).
- **Per-step error is rate-type, not truncation**: |H−H_bk| = tracking residual
  ≈ rate × dt_step. Measured scaling is linear in ds (err 3.3e-6 @ ds 5.4e-6 →
  err 6.5e-9 @ ds 1.3e-8, constant slope ≈ 0.6); shrinking ds cannot reduce it
  below the rate-implied floor.
- **Controller misreads it as truncation error** and collapses ds from the
  natural 8.3e-4 (≈P_eff/34) to the equilibrium ds ≈ 1.3e-8 — 5 orders below —
  where per-step error ~6.5e-9 finally sits under `ar-max-error` (1e-8).
- **Dead band pins it there**: growth gate requires error < 0.5×max (5e-9,
  `symplectic_integrator.h:3470`) but reduce fires at 1×max (`:3332`); the
  equilibrium error lands in (5e-9, 1e-8) — all observed warnings' per-step
  errors fall in this band. Net: ~7×10⁵–10⁶ substeps per Hermite call → the
  10⁶-step caps and minutes of wall time.
- Physical energy exchange is ~10¹³× too small to explain the H drift
  (pert_out/pert_in ≈ 6e-17) — the drift is bookkeeping, consistent with a
  κ-time-dependence residual (cf. SDAR lessons 2026-09-28 `correctSlowDownEnergy`
  fix, which removed one such leak; a residual remains under time-dependent κ).

### 2. Why it looked restart-specific: `--rand-seed` defaults to wall clock

`rand-seed` default 0 = "current cpu time" (`parallel-random/rand_io.hpp`), so
every run (and every functional-harness execution) rolls a new BSE stochastic
realization. Whether a given realization wanders into a pathological encounter
is chaotic. Same-seed experiments (2026-10-19, /tmp/bse_fresh):

| run | seed | result |
|---|---|---|
| fresh 0→10 | functional's exact seed 1792391053081014136 | **1.41 s, 0 warnings**; `f.3` bit-identical to functional `data.3` |
| restart pipeline (0→4, clear, resume 3→10) | same seed | **1.16 s, 0 warnings** |
| fresh 0→10, different (wall-clock) seed | — | grinds at t≈8.9 (10 warnings, minutes) |
| functional restart segment | wall-clock seed | grinds at t≈8.14–8.19 and 8.87–8.92 (≈6 min) |

Additional perturbation sources confirmed: `--write-group-info 2` changes the
saved random-stream state (`a.1.randseeds` ≠ `b.1.randseeds`) while particles
stay bit-identical until the next BSE kick event; even `-t 4` vs `-t 10` (same
seed) produces 864 differing bytes by t=3. Any perturbation re-rolls the dice.
Conclusion: restart reconstruction is *a* perturbation source, not *the* cause;
the original A-vs-C comparison was seed-confounded (each run wall-clock-seeded).

## Ruled out (2026-10-19)

- Redundant restart CLI options (`-b`, `-o`, `-f`, `--write-group-info`
  re-passed): the minimal canonical `petar -p c.par -i 0 -t 10 c.3` is
  equally slow. (The functional runner now uses the minimal form regardless.)
- `petar.data.clear` interference: physics identical with/without it.
- BSE / stellar-evolution state reconstruction: merge trajectory identical.
- **Cold AR/slowdown state at restart** (original hypothesis): during the grind
  κ is re-derived and at its cap; the grinding groups formed *after* the
  restart (trio group records only t=8.1357–8.1934; first n=3 group in the
  segment at t=5.26), so no group state crossed the restart.
- **Group geometry reconstruction** (r_search/r_group): member `t_record` fields
  synchronized; r_search values mass-scaled, not anomalous.
- **Root-orbit Kepler conversion bug**: recorded root elements (a=−1.66,
  e=1.0012) are self-consistent with the measured member state (v_rel≈v_esc)
  once G-units are handled correctly.

## Reproduction of the *mechanism* (deterministic)

The dump replay is the deterministic reproducer (no seed involved):

```bash
mkdir scratch && cd scratch
cp <repo>/test/out/functional_bse/data.par* .
cp <repo>/test/out/functional_bse/data.dump_large_step_ar_n3_* .
setarch $(uname -m) -R petar.mpi.omp.avx512.bse.hard.debug -p data.par \
    data.dump_large_step_ar_n3_t8.171875_M0_O0_c3_s1792391106
# → "Large AR step cluster found: total step: 5540018"; warning binary states
#   show SD at Tscale cap with SD_org ~1e12 (slowdown warm, ds collapsed)
```

The old simulation-level recipe (fresh vs restart from `petar.init`) is a
**coin flip** unless `--rand-seed` is pinned; keep it only for demonstrating
trigger variability, not as a regression test.

## Open items (suggested order)

0. **Two anomaly classes from the clean-4-body scan (analysed 2026-10-20,
   follow-up code work open)**:
   (a) *"fixed worse than old"* (perturber at 0.2/0.25/1.5): chaos-mediated
   trajectory lottery, not a fix defect — the fix's changed early ds
   decisions butterfly-diverge this knife-edge scattering system; in the
   fixed realization the inner binary is DISRUPTED (warning context: inner
   semi −1.97e-7, e = 2.15, scrambled tree, nonsense ds_init ≈ 35) and the
   scrambled tree defeats BOTH controllers (rp=0.4: bit-identical aborts;
   rp=0.2: old also grinds). Conservation on completed runs is comparable
   order (3 % vs 9 %). No code action; T7's window precedes the scramble.
   (b) *"oscillatory" ds-floor aborts* (perturber at 0.35–0.5, 5.0,
   bit-identical old/fixed — pre-existing, fix-independent; precise
   mechanism nailed 2026-10-20 evening by abort-frame member dump +
   conditional watchpoint): the encounter pumps the inner binary to
   e = 0.9955 (pericenter 1.9×10⁻⁹); during the pericenter plunge the
   rapid κ evolution (~5×10⁴) breaks the `de_sd` tracking in
   `correctSlowDownEnergy` (`etot_sd_ref_ += de_sd`,
   `symplectic_integrator.h:1698`): the reference is driven from −946
   across zero to +988, violating the audit invariant
   `etot_sd_ref ≈ ekin_sd + epot_sd` (true value ≈ −1.09e4) → the SD-channel
   H numerator `ekin_sd − etot_sd_ref` = −21.3 < 0 → `H_sd = ln(negative)`
   = NaN → the NaN branch halves ds unconditionally, but the NaN is
   structural → ds-floor abort. NOTE: `epot` stays negative in both
   channels and both log denominators remain positive — an earlier
   "hyperbolicity flips epot sign" reading was wrong and is corrected;
   the physical-channel H stays finite throughout. Suspected root: κ
   decomposition inconsistency between `ekin_sd_` (calcEKin at root κ) and
   `epot_sd_` (calcAccPotAndGTKickInv at node κ) in the force-sync branch —
   same family as the 2026-09-28 snapshot-ordering fix, on the κ side.
   Fix direction: on persistent H-NaN / grossly violated reference
   invariant, fall back to the physical energy-error channel or re-set the
   sd reference instead of halving; independently verify de_sd consistency
   under rapid κ change. **Spun off 2026-10-20 into a dedicated research
   thread, now archived as the reference note
   `SDAR/docs/etot_sd_ref_corruption_fix_notes.md` (investigation record +
   implemented fix, root-caused 2026-10-21; open follow-ups in
   `SDAR/docs/etot-sd-ref-tracking-corruption-plan.md`)
   (measured facts, suspected κ-decomposition root, git-archaeology leads —
   incl. the 881abf3 syncTreeSlowDownAndDs refactor as prime suspect —
   deterministic reproducer with dump/pars/gdb scripts in
   `~/sdar_run/clean4/`, and open items incl. a commit bisect of the
   etot_sd_ref_ update history).**
1. ~~**SDAR controller fix**~~ **Executed 2026-10-19** — rate-type drift
   separation implemented in `SDAR/src/AR/symplectic_integrator.h`
   (`integrateToTime`): per-step error becomes
   `max(0, |H−H_bk| − 1.5·|R̂|·dt)` with `R̂` a signed drift-rate EMA
   (`h_drift_rate_ema_`, running-mean warm-up then α=0.1, active after 3
   samples); reduce and growth gates both use the effective error; warning
   output adds `|Int_err_eff/E|` and `|H_drift_rate|`. Verification:
   deterministic dump replay 5,540,018 → **16,393** substeps (0 warnings,
   physical dE same order, slightly better); second grinding realization
   (34.87 Msun binary dump) → **183** substeps per cluster step; clean-run
   same-seed rerun diverges only via adaptive-path bifurcation (expected) with
   cleaner conservation books; validation T1–T4 pass (T5–T6 see below).
   **Still open within this item**: the κ-time-dependence bookkeeping residual
   itself (~2 % per update) is now tolerated by the controller but not
   eliminated — extends the 2026-09-28 `correctSlowDownEnergy` work.
   **Stability campaign (2026-10-19, `~/sdar_run/triple-stability-fixcheck`)**:
   35 isolated-triple runs (20 on `ar.logh.ttl.sd.cm`, 14+1 on `hermite`)
   covering hierarchy 1e2–1e6, e_in 0–0.9, e_out 0→1.5 (circular→near-parabolic
   →hyperbolic), inclinations 0/90/180°, three mass configurations, plus the
   **exact PeTar grind-trio state** extracted from the dump and replayed
   standalone. Result: old and fixed controllers produce **bit-identical
   outputs (step counts and final energies) in all 35 runs** — the fix is
   provably inert outside the rate-type regime.
   **ROOT CAUSE FULLY CLOSED (2026-10-20)** — both remaining links resolved:
   (1) *Cap identity*: no cap switch — inner-binary stab=0 in BOTH contexts, so
   the stability cap is inactive and ts = √t_min² (encounter timescale) is the
   binding cap everywhere; κ is updated per-step by the SAME code path in both
   contexts (hardware-watchpoint backtrace: `calcSlowDownInnerBinary` ←
   `calcBinaryTreeSlowDown` ← in-loop `syncTreeSlowDownAndDs`, plus the
   call-entry sync; standalone ts drifts per-step even faster, +3.6 %/step
   during the flyby transient).
   (2) *Amplifier*: the **perturber list**. The H error estimator is built on
   internal energies only (no perturbation-potential term), so the physical
   energy flow W between the external perturbers (n_pert=18 in-cluster) and the
   group appears per-step as |ΔH| ≈ W·dt/|E| — a rate-type error no ds can
   reduce. Without perturbers W=0: the isolated trio and the 80-particle
   emulation with the DEFAULT neighbor radius (which left the group's neighbor
   list empty — that is why the emulation originally looked clean) show
   pure-truncation errors ~2e-11. Re-running the emulation with populated
   neighbor lists (`--r-neighbor-over-group 10/30`) reproduces the regime
   standalone: per-step errors jump to 3.5e-8 and `Large_energy_error` events
   appear. All magnitudes are consistent: in-cluster W≈0.2–14/unit gives
   err=3.6e-5 at formation (dt 3.7e-4) and 6.5e-9 at equilibrium (dt 5.6e-9).
   The κ drift is the same physical flow expressed in the slowdown gauge
   (κ tracks the perturber-driven ts_in; identity dE_sd/dt=−(dκ/κ/dt)·E_sd),
   and the earlier "external work too small by 10×" argument was a misreading
   — the controller reacts to the per-step increment, not the mean drift.
   Consequence: the EMA fix subtracts exactly the missing open-system energy
   flow — principled, not symptomatic. A future alternative is including
   epert in the conserved-quantity check, but that changes error semantics
   globally.

   **Self-made-IC verification (2026-10-20, HardDump crafting)**: by cloning
   the dump binary layout (PtclH4 stride 240B; mass@0/pos@8/vel@32/
   bin_state@56/r_in@224/r_out@232) and replacing particle records, the
   trigger decomposes into **seed × activator**:
   | configuration | old | fixed | verdict |
   |---|---|---|---|
   | real trio alone (no field) | 57 | 57 | trio alone clean — field required |
   | **real trio + 16-star synthetic field (self-made)** | **16,062,968 (3 warns, 398 LEE)** | **1,245 (0 warns)** | **self-made IC triggers; fix cures 12,900×** |
   | synthetic trio (near-parabolic) + real field | 764,891 | 764,891 | legitimate cost; fix correctly inert |
   | synthetic trio + synthetic field | 1,415,259 | 1,415,259 | legitimate pericenter dive; fix inert |
   Seed = the real trio's knife-edge state (intruder barely unbound,
   lingering at ~4×10⁴·a_in with v≈v_esc — produced by the chaotic cluster
   capture; nominal-geometry synthetic trios do not land in this corner,
   which is why 35 synthetic configurations never triggered). Activator =
   any adequate external field (real or synthetic). Tidal tensor ruled out
   twice (`--soft-perturbation 0` bit-identical).

   **AMPLIFIER = the changeover weighting of perturber forces (2026-10-20,
   neutralization experiment)**: same seed+field, members' r_in/r_out
   enlarged to 1.0/1.2 so every field star sits deep inside the k=1 zone
   (forces at FULL strength, only the r-dependent modulation removed) —
   the grind VANISHES: 16,062,968 → **650** substeps (24,712×), 0 warnings,
   with the perturber forces now STRONGER than in the grinding run. So the
   rate-type error floor is generated by `calcAcc0WTwo`-weighted
   perturbation forces (PeTar `ARInteraction`; the trio members' r_out
   ~0.65–1.1 places the cluster stars inside the transition zone, and each
   substep's CM motion modulates k(r), injecting a smooth energy-flow
   signature into the internal-only H audit). This also explains why the
   standalone sample can only partially reproduce: `sample/Hermite`'s
   interaction has NO changeover — with the exact seed+field it caps at
   err≈2e-7 spikes (`AR_KDK_PERT` alone gives 6×; TTL is inert in the
   sample), ~170× below the PeTar-level continuous 3.6e-5 and without a
   sustained collapse. Complete chain: knife-edge trio (seed) × external
   field inside the members' changeover transition zone (activator ×
   amplifier) → modulated open-system flow → per-step H error ≈ rate·dt →
   old controller collapses ds; the EMA fix removes exactly that rate
   (12,900× cure on the self-made trigger, inert on legitimate cost).

   **Minimal-system localization (2026-10-20, single-perturber scan)**:
   real trio + exactly ONE perturber (m=1, static, on +x), everything else
   bystandered — cost vs perturber distance is a clean bell over the
   members' changeover transition band, baseline outside both ends:
   | r_p | zone | old substeps | fixed |
   |---|---|---|---|
   | 0.03 / 0.05 | inside r_in (k=1) | 774 / 546 | — |
   | 0.15 | transition edge | 12,681 | 727 (17×) |
   | 0.3 | transition | 566,265 | 8,331 (68×) |
   | 0.5 | transition center | **2,065,841** | **12,375 (167×)** |
   | 0.9 | transition far edge | 447,789 | — |
   | 1.5 | outside r_out (k=0) | 53 | — (=trio-alone baseline) |
   | 0.3, m=4 | transition | 579,158 | mass-independent |
   | 0.3, neutral k(r) | weight removed | 43 | baseline |
   Precise trigger condition: knife-edge trio seed × ≥1 perturber inside
   the members' changeover transition band (the k(r) modulation channel);
   perturber mass nearly irrelevant (modulation, not strength); the fix
   cures across the band while preserving the mild legitimate cost gradient.
2. **PeTar trigger-side mitigation** (optional, cheap): consider breaking a
   3+-member group early when the inner binary's slowdown is already at its
   timescale cap and the extra member only drives κ drift without approaching
   (binary + perturber treatment is the cheaper decomposition).
3. ~~**Functional harness stabilization**~~ partly done — the deterministic
   validation case exists: **T7** (`test/validation/scenarios/t7_changeover_
   rate_drift.json`, `--pipeline t7`): minimal 4-body plain-format system
   (knife-edge trio + one static perturber at 0.15, inside the intruder's
   changeover band; full-precision seed; half-tree-step window), 5 checks
   (no step-cap warnings / ds-floor aborts / runtime aborts, plus an
   inside-r_in false-positive guard), verified discriminative (FAILS with
   the pre-fix binary, 5/5 PASS with the fix, ~1 s suite cost). Remaining:
   pin `--rand-seed` for the functional `bse` case (its runtime is still a
   random variable, 1.2 s … 6+ min).
4. ~~**Verify**~~ dump replay verification done (see item 1); full validation
   suite regression recorded with the fix (T1–T6; see `test/out/fixcheck_*`).

## Cross-references

- SDAR lessons-learned 2026-10-19 entry (AR ds-controller collapse, gdb
  measurement template) — mechanism owner.
- PeTar lessons-learned 2026-10-19 entry (wall-clock seed confound; same-seed
  comparison protocol).
- Deterministic reproducer: `test/out/functional_bse/data.dump_large_step_ar_n3_*`
  (regenerable from any grinding realization; see recipe above).
- Related src areas: `SDAR/src/AR/symplectic_integrator.h` (`integrateToTime`
  error gates :3332/:3470, `syncTreeSlowDownAndDs`, `correctSlowDownEnergy`),
  `PeTar/src/hard.hpp` (defaults `ar-max-error` 1e-8, `ar-max-nstep` 1e6),
  `PeTar/parallel-random/rand_io.hpp` (wall-clock default seed).
