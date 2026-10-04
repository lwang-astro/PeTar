# SEVN Integration: Status, Defects, and Development Plan

Status: 2026-10-04 — SEVN runs work end to end for no-primordial-binary
clusters (N=1k, 12 Myr: 38 type-change + 9 sn_kick events recorded);
binary-rich scenarios are blocked by a pre-existing, family-independent IC
hard failure (`doc/ic-hard-failures-plan.md`). Resolved items and false
alarms are recorded in
`.github/skills/petar-nbody-simulation/assets/lessons-learned.md`, not here —
this document lists only open defects and questions.

## 1. Confirmed Defects (Require Development)

### D1 Random seed not wired into PeTar's parallel random system — High priority
- Symptom: the SEVN branch of `BSEManager::initial()` wrote the seed into a local
  `struct {int idum;} value3_` (SEVN declaration block of bse_interface.h), but
  **the SEVN library never reads that struct** (it uses its own `utilities`
  random numbers); meanwhile `construct_default_sevn_params()` drops `rseed`
  via `unneeded_options`, so users cannot pass a seed to SEVN at all.
  (The dead `value3_`/`rand3_` plumbing and the unused `-idum` option were
  removed on 2026-10-04 — see lessons-learned.)
- Consequences: (1) random processes such as SN kicks in SEVN mode are
  **not reproducible** and the seed is uncontrolled run-to-run; (2) `rand_manager`
  (data.par randseeds, per-thread seeds) still writes its files in SEVN mode,
  creating the illusion of seed management with no consumer; (3) under MPI the
  SEVN-internal random sequence may repeat across ranks (unverified).
- Plan: expose `rseed` in `default_params_sevn` and fill it in
  `BSEManager::initial()` from the per-rank/per-thread seeds of `rand_manager`.

### D2 GW merger recoil (gw_kick) not wired into the SEVN path — High priority
- Symptom: the experiment-branch GW recoil chain (`GWKick::calcKickVel/calcFinalMass`,
  `getCompactChiRandom`, `compactOspinToChi`; shared region bse_interface.h
  L1855-1890) is only called from the **non-SEVN branch** of `evolveBinary()`
  (the merger_event block, L2284-2335); the SEVN branch goes through
  `evolve_bco` (Peters 1964) and lets `merge_SEVN` handle the merger internally,
  **producing no** PeTar-side GW kick velocities, remnant mass ratios, or
  event_flag=6 events.
- Consequences: in SEVN mode, BBH/BNS mergers lose the experiment branch's
  self-consistent recoil physics and the `fout_bse_gw_kick` event stream;
  post-processing `BSEKick`/GW statistics are missing.
- Plan: after the `merge_SEVN` call site in the SEVN branch, reuse the shared
  gw_kick chain (SEVN's RemnantType/spin fields can be mapped); or confirm on
  the SEVN library side that its built-in kick handling applies and align the
  event record format.

### D3 `spin_3d` semantics inconsistent on the SEVN path — Medium priority
- C++: `StarParameter::ospin[3]` (experiment 3D spin) only has `ospin[0]`
  written under SEVN (all scalar uses in evolve_sevn.cpp; `ospin[1..2]` stay 0).
- Python: `SEVNStarParameter` nonetheless defaults to `spin_3d=True` (3 columns);
  `spin_3d=False` is the 1-column form — opposite in direction to the bse family
  (for bse data 3D is the new format; for SEVN data 1D is the only format).
- Consequences: snapshot column-width interpretation can be confused; the
  dimensionless BH chi (3D) simply does not exist in SEVN data.
- Plan: document + set class defaults so that "SEVN is always scalar spin", or
  ignore `spin_3d` for SEVN and assert.

### D4 Binary event file streams unverified under SEVN — Medium priority
- SSE streams are verified working (12 Myr no-binary run: 38 `type_change` +
  9 `sn_kick` events on `data.sevn.*`); all 9 `.sevnB` streams were opened but
  stayed empty — no primordial binaries, and dynamically formed binaries
  produced no logged BSE-binary events within 12 Myr.
- Plan: verify `.sevnB` streams with a binary-rich run (blocked until the IC
  hard failure in `doc/ic-hard-failures-plan.md` is resolved); then mark each
  BSE event type's SEVN behaviour (native / missing / equivalent). Note the
  gw_kick/binary_merge/tide write points sit on non-SEVN paths (D2) — expected
  empty under SEVN regardless.

### D5 TDE criterion field assumptions unverified under SEVN — Medium priority
- `getTDERadius()` (bse_interface.h L2421, shared region) relies on `kw`
  (SSE type), `mt/r/mco`; SEVN provides kw via the PhaseBSE mapping, but the
  semantics of `mco` (SEVN's MCO_SEVN vs BSE's mc) are **unchecked**; the actual
  triggering of the ar_interaction TDE interrupt chain (L1590-1612) under a SEVN
  build is untested.
- Plan: construct a CO+star close-encounter binary case and verify that TDE
  events are produced and recorded in SEVN mode.

### D6 Remaining BSE→SEVN capability comparison (initial survey, pending runtime confirmation)
| BSE-family feature | SEVN status |
|---|---|
| SN kick (velocity + direction) | Native (in-library); SSE `sn_kick` events verified (9 in the 12 Myr run); direction randomness affected by the D1 seed issue |
| GW merger recoil + remnant mass | **Missing** (D2) |
| Hyperbolic/binary TDE | Criterion shared but unverified (D5) |
| Mass transfer / CE / common envelope | Native in SEVN (evolv2_SEVN/merge_SEVN), semantics differ from BSE, events via the BEvent columns |
| Dynamical (hyperbolic) merger records | Unverified (depends on the event-flag chain working in the SEVN branch) |
| 3D spin/chi | **Missing** (D3, scalar only) |
| Metallicity / table coverage | Limited: MIST **gold** set only (verified 0.7–80 M☉ at Z=0.00142857, discrete Z ≤0.00452; usage details in SKILL.md Gate 4 and `sample/star_cluster_plummer_N1k_binaries_sevn.sh`) |
| Random-number reproducibility | **Missing** (D1) |

## 2. Tooling Leftovers (Low priority)
- `petar.select --require sevn` is supported, but the actual selection flow for
  the sevn binary family is untested (regression-run once after a sevn build is
  installed).
- The sevn variant of `petar.bse.get.init.binary` (`petar.sevn.get.init.binary`)
  has not been adjusted for SEVN parameters (still BSE-biased assumptions);
  review before use.

## 3. Open Questions

1. The installed SEVN library is recent (built with cmake 4.4.3); divergence
   from the version the upstream PR author developed against is unevaluated —
   if `evolve_sevn.cpp`'s use of the SEVN API has version-sensitive spots,
   check them against the SEVN changelog (low priority until a SEVN upgrade).
2. Binary-rich verification (D4) is blocked by the pre-existing
   primordial-binary IC hard failure — not a SEVN defect; deterministic
   reproduction ticket: `doc/ic-hard-failures-plan.md`.

## 4. Recommended Implementation Order
1. D1 seed wiring (reproducibility is a prerequisite for scientific output);
2. D2 GW recoil wiring (the core physics gap);
3. D4+D5 event-stream and TDE verification (D4 needs the binary scenario
   unblocked, `doc/ic-hard-failures-plan.md`);
4. D3 spin semantics, §2 tooling leftovers.
