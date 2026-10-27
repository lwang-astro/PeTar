# Plan: TT (tidal-tensor) tree-channel inner-binary degradation — CLOSED

Status: **CLOSED 2026-10-27** — root cause found (F13), fix implemented and
validated (full T1–T6 campaign green with the fixed production binaries);
this file now holds only the remaining open follow-ups plus pointers to where
the completed content was distilled. The investigation record (facts F1–F16,
discriminator program, fix evidence) lives in git history of this file prior
to 2026-10-27 if ever needed; the durable knowledge was distilled to:

- `doc/tidal-tensor-design_notes.md` — mechanism: fit input invariant
  (c.m.-subtracted accs, exactly once per force refresh), the KDKDK_4TH
  acorr-phase guard, and the member-mean-kept-intact application (1795e
  revert) with its physical justification.
- `.github/skills/petar-nbody-simulation/assets/lessons-learned.md`
  (2026-10-27 entries) — the cross-layer invariant lesson, the
  cap-co-variation sweep rule, and the amplitude-vs-correlation methodology
  lesson.

## Closed with this work (one-line pointers)

- Tree-channel TT inner-binary degradation (25x/140x) → fixed by the acorr
  guard (`src/hard.hpp`) + 1795e revert (`src/ar_interaction.hpp`);
  validated: T4 base m5 Δe 8.214e-3 vs truth 8.220e-3, ctrl m3/m5 at truth
  level, outer-orbit tracking corr ≥ 0.987, T1–T6 all green (2026-10-27).

## Open items (follow-up, separate commit)

1. **T4 scenario hardening (was Q4/Q5)**:
   - add a Delta_e bound for TT modes vs m1 (the current suite passes on
     census/aborts only — it did not catch the 25x regression);
   - add an outer-orbit tracking check (e_out(t) correlation vs m1 ≥ 0.9 —
     the discriminator that separated physical signal from suppression);
   - census guard: N_all counts artificial particles only; extend the census
     to cluster/direct perturbers so the m4 baseline's tide delivery is
     verified (F2 correction: m4 is NOT tide-free — its inner e(t) tracks m1
     at corr 0.9996 through the cluster path);
   - decide whether `r_search_group` semantics (tree-side perturber reach for
     AR groups) need a dedicated direct-perturber tree mode, or the scenario
     should enlarge r_search_group to make m4 a true direct-tide baseline.
2. **Optional hardening**: debug-gated leak-detector assert in the TT fit
   (`|T3|·dc² ≪ C·|T2|·dc`, C~10) — would have caught F13 instantly; not
   implemented.
