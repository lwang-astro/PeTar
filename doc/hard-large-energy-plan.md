# Plan: hard_large_energy remaining actions (diagnosis closed 2026-10-04)

The diagnosis is closed and distilled into
`.github/skills/petar-nbody-simulation/assets/lessons-learned.md`
(Simulation Execution, 2026-10-04 entry): root cause is group form/break
transition injection at the criterion boundary (NOT SE mass loss — falsified
by dm=0 event replays); the switch position is the dominant lever; production
is radius-capped (`r_group/r_in` clamps the kappa branch's much-earlier switch
requests). SDAR-side investigation closed 2026-10-03 with zero code changes
(the transition layer is already pure, grouped beats ungrouped 32x; residual
is intrinsic) — findings in SDAR lessons-learned (2026-10-03 entries; full
record: git `d474412`); regression bench `SDAR/sample/test/group_transition_bench.sh`.

Source runs: BSCC-M9 Pal 5 IMF grid (`petar.mpi.omp.avx512.bse.galpy`, SDAR
`649643e`+`7d12afb`, 2026-09-27 build), 7 runs, 1–12 Gyr. The existing
hard_large_energy dumps in these runs are noise from this known source and
can be cleaned.

## Remaining actions (in priority order)

1. **Production switch-radius test segment** (only open item from this
   diagnosis): short restart segments of a dense model (e.g. N210k) with
   `--r-group`/`--r-search-group` scaled x3 and x10; measure hard_large_energy
   dump rate, dE distribution, group counts, AR step profile and wall-clock
   cost. Adopt the multiple with the best dump-reduction-per-CPU. Note the
   `*.hard.debug` replay cannot substitute: criterion radii are embedded in
   the dump's per-particle changeover state and ignore par-file edits
   (measured: identical r_crit across par factors 0.1–10).
2. **dm-tracking gap** (independent, inactive): dm tracking is active in
   `writeBackPtclForOneCluster(OMP)` (`src/hard.hpp`: records
   `sys[adr].dm = mass_new - mass_bk`, fills `_mass_modify_list`, ~lines
   3170/3215) but the consumers are disabled
   (`correctSoftPotMassChangeOneParticle` calls commented out, ~lines 2697
   and 4266); SSE wind mass changes inside the hard domain are only corrected
   at the soft level (`correctSoftPotMassChange`) today, so the hard-side
   `de_modify_single` path records only tiny amounts.
3. **Warning-criterion semantics** (user decision, no physics change):
   trigger on `dE_change` instead of reference-relative `dE` to silence the
   transition-seam dumps entirely (~99% of volume); orthogonal to action 1.
