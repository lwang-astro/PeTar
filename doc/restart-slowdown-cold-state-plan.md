# Plan: AR ds-controller collapse — the "restart slow" grind (open items)

Status: plan — open items only. Opened 2026-10-19 (functional `bse` case:
6.3 min vs 1.4 s); root cause identified 2026-10-19; fix thread closed
2026-10-22 by the SDAR root fix (commit `6a8a357`, which superseded and
removed the interim 2026-10-19 EMA compensation). The investigation,
mechanism, and validation record is archived in the SDAR reference note
`SDAR/docs/etot_sd_ref_corruption_fix_notes.md`; this file keeps only the
PeTar-side open items and pointers (the full investigation text is preserved
in this file's git history).

## Open items

1. **Optional trigger-side mitigation** (cheap, undecided — pending
   discussion): break a 3+-member group early when the inner binary's
   slowdown is already at its timescale cap and the extra member only
   drives κ drift without approaching — the binary + perturber
   decomposition is cheaper.

## Where the completed content went

- Mechanism, implemented fix, deterministic g4e reproducer, production
  N=1000 cluster A/B validation: SDAR reference note
  `SDAR/docs/etot_sd_ref_corruption_fix_notes.md` (registered in SDAR
  AGENTS.md → Documentation).
- SDAR-side follow-ups (LogH flat-path κ-weighted kick credit, LogH runtime
  validation of the negative-dt repair routing, hyperbolic-root ds_init
  estimator): SDAR plan `docs/etot-sd-ref-tracking-corruption-plan.md`.
- Prevention rules & diagnostic methodology: SDAR lessons-learned
  2026-10-19 (#5/#6), 2026-10-20 (#8/#9), 2026-10-21; PeTar lessons
  2026-10-19 (wall-clock seed confound + same-seed protocol;
  cold-state-hypothesis correction).
- Regression coverage: validation **T7** scenario
  `t7_triple_perturber_ar_ds_adjustment` (`--pipeline t7`; 22 checks —
  LogH/BTLogH gauges, six zone placements: 0.05/0.15/0.30/0.35/0.50/0.60,
  i.e. band inside/mid/far, the former ds-floor-abort zone, and beyond
  r_out). Measured 2026-10-22 on the root-fix binary: 2.27×10⁷ substeps
  pre-fix → 1.25×10⁴ (LogH) / 1.15×10⁵ (BTLogH) at the trigger, all
  placements complete with zero warnings/aborts; full T1–T7 regression
  clean.
- Functional determinism (2026-10-22, executed): the functional `bse` case
  pins `--rand-seed 1792391053081014136` in
  `test/functional/functional_matrix.json` — verified non-grinding on the
  root-fix binary (case pass, 9.2 s wall, zero step-cap warnings); the
  wall-clock-seed runtime lottery (1.2 s … 6+ min) is closed.
- Pre-fix grind-trio dump replay (needs the `bse.hard.debug` family):
  `test/out/functional_bse/data.dump_large_step_ar_n3_*` +
  `setarch $(uname -m) -R petar.mpi.omp.avx512.bse.hard.debug -p data.par
  data.dump_large_step_ar_n3_t8.171875_*` — regenerable from any grinding
  realization; superseded for regression purposes by T7.
- Related src areas: `SDAR/src/AR/symplectic_integrator.h`
  (`integrateToTime`, `syncTreeSlowDownAndDs`, `correctSlowDownEnergy`),
  `SDAR/src/AR/information.h` (fix_step_option grants),
  `PeTar/src/hard.hpp` (ar-max-error 1e-8, ar-max-nstep 1e6),
  `PeTar/parallel-random/rand_io.hpp` (wall-clock default seed).


