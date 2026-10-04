# hard_large_energy dumps fired by SE mass loss inside few-body systems (handoff from Pal5-IMF production runs)

Source: BSCC-M9 Pal 5 IMF grid (`petar.mpi.omp.avx512.bse.galpy`, SDAR `649643e`+`7d12afb`, 2026-09-27 build), 7 runs, 1–12 Gyr in progress. Diagnosed 2026-10-02 with `*.hard.debug` replay + gdb. Repro materials on this machine: `/data/lwang/Pal5_IMF/<model>/debug/` (dump + current `data.par*`).

## Symptom

`data.hard_large_energy_h4_*` warning dumps accumulate at up to ~40/hour in dense clusters (~11.8k files across 7 runs in 3 days; worst: N730k_Rh3.8 7199, N210k 2787, N550k_Rh4.8 1063). All replayed events share one pattern: the trigger fires on `dE` vs the step-start reference while the true in-step integrator error is 4–5 orders smaller, and — the decisive detail — **`dE` is already at its final value at the FIRST substep print and stays constant across the whole step** while `Ekin`/`Epot` evolve smoothly: the offset is baked in at the step boundary, not grown by integration. Slowdown factor is 1 in all three (`Ekin_SD==Ekin`, `Epot_SD==Epot`), so `dE_SD==dE`.

## The three cases (replayed, values match run-side dump records)

**Case 1 — `N550k_Rh4.8_a-2.0` n4 g0, t=369.003296** (`..._h4_n4_g0_ke2942b7b0015a98c_t369.003296_M3_O13_c0_s1790918131.tmp`)
- Final substep: `Ekin 53.0942, Epot -21.9547, Etot_SD_ref 31.5407, dE = -0.401199 (dE/Etot_SD = -0.01288)`
- `dE_change = 3.03e-06`, `dE_change_single = -2.13e-10`; `H4_step_sum 744, H4_force_sum 2243, sdar_substep_sum 124, sdar_tsyn_step_sum 41`
- Substep sequence: `Ekin 54.13 → 53.09` over ~10 printed substeps while `dE` pinned at `-0.401196` the entire time (first print already has it). System has `Etot = 31.14 > 0` (unbound 4-body); `Etot_SD_ref = 31.5407` sits ABOVE the whole trajectory — opposite sign to cases 2/3, so this event LOWERED the energy vs the reference (interpretation left open: not every SE channel raises E).

**Case 2 — `N550k_Rh4.8_a-2.0` n6 g0, t=369.487061** (`..._h4_n6_g0_k36338a6ce0bb4424_t369.487061_M5_O3_c0_s1790918567.tmp`; the "dE_SD=16.4" event)
- Final substep: `Ekin 3974.39, Epot -4455.14, Etot_SD_ref -497.156, dE = +16.4051 (dE/Etot_SD = -0.03412)`
- `dE_change = 1.10e-04`, `dE_change_single = -2.89e-10`; `H4_step_sum 1073, H4_force_sum 5357, sdar_substep_sum 99, sdar_tsyn_step_sum 31`
- Same constant-offset pattern (`16.4053 → 16.4051` across substeps, substeps at dt ~9.5e-7). Tight, deeply bound 6-body (`|E| ~ 481`); positive dE = system became LESS bound than the reference — classic wind/AGB mass loss weakening the potential. t=369 Myr: the ~1.5–2 M☉ population is in AGB/WD-formation era, so envelope loss inside a compact few-body is the natural source.

**Case 3 — `N210k_Rh5.85_a-2.3` n4 g0, t=1199.444336** (`..._h4_n4_g0_k5a212fdc596a433a_t1199.444336_M4_O21_c0_s1790918754`, post-rename non-.tmp)
- Final substep: `Ekin 29.4307, Epot -20.4674, dE = +0.138423 (dE/Etot_SD = 0.01544)`
- `dE_change = 3.13e-06`, `dE_change_single = -6.46e-11`; `H4_step_sum 711, H4_force_sum 2149, sdar_substep_sum 191, sdar_tsyn_step_sum 62`
- Same pattern; t=1199 Myr again the WD-formation mass-loss era.

gdb stack (case 2): `main (hard_debug.cxx:528) → HardIntegrator::integrateToTime` — no anomalous path.

## Root-cause pointer

A mass-change energy correction is **designed but not fully wired** in `src/hard.hpp`:

- dm tracking is ACTIVE in `writeBackPtclForOneCluster(OMP)`: records `sys[adr].dm = mass_new - mass_bk`, fills `_mass_modify_list` (~lines 3170 / 3215).
- The consumer is DISABLED: `correctSoftPotMassChangeOneParticle(...)` calls are commented out (~lines 2697 and 4266, `//// correct soft potential energy ... due to mass change`).

Expected fix directions (for the dev session to decide): apply the recorded dm to the energy bookkeeping — either update `Etot_SD_ref` / add the correction term to the large_energy check, and/or restore the soft-potential correction calls. Alternative cheap mitigation: trigger on `dE_change` instead of the reference-relative `dE` (would remove ~99% of the dump volume without physics change).

Replay command (family-matched, ASan build):

```bash
cd <scratch with data.par*> && setarch $(uname -m) -R petar.mpi.omp.avx512.bse.galpy.hard.debug <dumpfile>
```

Note: the replay tool prints its built-in `G = 0.0044985`; the runs used `-G 0.00449830997959438`. Replayed energies still match the run-side records.

## Follow-up (2026-10-02, dev session): root cause corrected — group-event transition injection, NOT SE mass loss

Event-level replay of all three cases **falsifies the SE-mass-loss hypothesis**:

- gdb on case 2: 12 `evolveStar` calls in the step, ALL with dm = 0 (BHs) or ~1e-14 (low-mass stars); `.sse.*` event files empty. No in-step mass loss.
- Instead each step contains **group form/break events**: a massive pair (40.5+40.5 Msun BHs in cases 1/2) passes peri of an unbound / marginally bound encounter (case 1: semi -2.0e-3 e 1.017; case 2: semi -0.835 e 1.00001; case 3: semi +1.06 e 0.9998). At peri d < r_crit momentarily -> a 2-member group forms (to protect the close passage) -> breaks as soon as d crosses r_crit. Cases 1 and 3 show form -> break -> **re-form** -> break chatter around the boundary.
- Each transition rewrites the integrator state (CM<->singles, slowdown/ds/vcm_record re-init) and injects O(1e-3) of the interacting-pair energy: case 2 +16.4 (~6e-4 of the encounter Ekin), case 1 1.12, case 3 0.138. This is the transition-injection layer measured earlier with in-memory palindromes (SDAR lessons 2026-09-30), now confirmed as the dominant large_energy source in dense clusters (~40 dumps/hour; the per-interval suppression keeps it bounded but every distinct close flyby fires once).
- `dE_mod`/`dE_change` stay ~0 because transition injections were never booked into the change columns — which is why the step-start-reference `dE` carries the full offset and re-fires the warning.

Fix applied (PeTar `src/hard.hpp`, integrateToTime energy readout): when a cluster step had group form/break/merge events, **re-baseline the energy references** after the readout (`calcEnergySlowDown(true)`), move the injected offset into the `dE_change` bookkeeping columns, and let the warning check see only the residual integration error. Verified: all three replays now report `dE ~ 1e-10` with zero warning triggers while the injection stays visible in `dE_change` (case 2: 291.7 accumulated over its event steps); genuine error detection is preserved for non-event steps and for post-event drift.

**Decision 2026-10-03: the re-baseline fix is withdrawn.** `de_change_cum` was designed to record *physical* energy changes (stellar evolution mass loss, binary interrupts); folding algorithmic transition injection into it conflates physics with numerics and masks the error. The warning pressure is the correct signal that the transition layer injects real error. Root-cause direction adopted instead: make the SDAR group-transition layer per-event invertible with pure-state derivations (form/break as inverse canonical changes of variables, all auxiliary quantities derived from the current state by shared code, switches at common sync boundaries) — detailed executable plan: `SDAR/docs/transition_unification_plan.md`. Until that lands, the production dump flood is capped only by the per-interval suppression; the existing dumps in the seven runs are noise from this known source and can be cleaned. The dm-tracking items above (disabled `correctSoftPotMassChangeOneParticle` consumers) remain a separate, real but currently-inactive gap: SSE wind masses change inside the hard domain and are only corrected at the soft level (`correctSoftPotMassChange`), so the hard-side `de_modify_single` path records only tiny amounts today.

Remaining known issue (not addressed here): the injection itself is real energy error — bounded but nonzero per event. Reducing it requires the transition-layer unification (same state-rewrite function for form and break), tracked as the next implementation goal in the SDAR palindrome notes.

## Closure (2026-10-04): diagnosis complete, root cause refined, remaining actions

The SDAR-side investigation of the transition layer is complete (full evidence
chain and experiment record: `SDAR/docs/transition_unification_plan.md`, commits
`a2dff7a` + `d474412`; regression bench: `SDAR/sample/test/group_transition_bench.sh`).
Final state of every hypothesis raised above:

| Hypothesis | Verdict | Evidence |
|---|---|---|
| SE mass loss | falsified (2026-10-02) | dm=0 in all evolveStar calls; `.sse.*` empty |
| Transition bookkeeping asymmetry (#4/#5) | already fixed by the 2026-09 consolidation; no net injection (form/break one-shot bookings cancel pairwise) | per-event decomposition trace; measured `dE_rewrite` <= 1e-16 everywhere |
| Re-entry initialization split (#7) | already unified (same force/dt init path, block-aligned) | code-path audit + state-purity measurements |
| AR adaptive ds sequence (#6) | ruled out | fixed-ds (`-e 1e10`) run bit-identical to adaptive |
| Remaining per-event offset (0.138..16.4) | **intrinsic local truncation of switching discretizations at the switch position** — dominated by WHERE the switch happens | matched-control palindromes |

Key discovery: the **switch position is the dominant lever** (palindrome
round-trip energy error: 400x better switching early in isolated encounters,
~3x with an intermediate optimum under a strong perturber, flat for bound
orbits). The unified criterion's kappa branch already implements adaptive
placement and is calibrated at the measured optimum (kappa gate 1e-2): all
sample optima form exactly at the kappa crossing. **Production is
radius-capped**: all three Pal5 cases form exactly at the radius boundary with
kappa_org 7–13 orders above the gate, i.e. the kappa branch requests a much
earlier switch and `r_group/r_in = 0.00375` clamps it.

Also falsified: "forming groups under extreme perturbation is harmful" —
grouping wins in all six tested environments (6x..1e8x); the worth-forming
gate direction is closed.

### Remaining actions (in priority order)

1. **Production switch-radius test segment** (only open item from this
   diagnosis): short restart segments of a dense model (e.g. N210k) with
   `--r-group`/`--r-search-group` scaled x3 and x10; measure hard_large_energy
   dump rate, dE distribution, group counts, AR step profile and wall-clock
   cost. Adopt the multiple with the best dump-reduction-per-CPU. Note the
   `*.hard.debug` replay cannot substitute: criterion radii are embedded in
   the dump's per-particle changeover state and ignore par-file edits
   (measured: identical r_crit across par factors 0.1–10).
2. **dm-tracking gap** (independent, inactive): the disabled
   `correctSoftPotMassChangeOneParticle` consumers (see Root-cause pointer
   above) remain unwired; hard-side wind mass corrections only land via the
   soft level today.
3. **Warning-criterion semantics** (user decision, no physics change):
   trigger on `dE_change` instead of reference-relative `dE` to silence the
   transition-seam dumps entirely (~99% of volume); orthogonal to action 1.
