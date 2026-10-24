# Validation Framework (T1–T8)

**Design rules**: all validation reports follow [REPORT_DESIGN.md](REPORT_DESIGN.md) (paper-style
sections, symbol definitions, figure encoding, threshold semantics, test-design principles, shared
`report_common.py` helpers). T6 is the reference implementation.

Scriptable validation scenarios for PeTar algorithm changes.

## Quick Start

**Default (run everything)**: `python3 test/validation/run_validation.py` with no arguments runs the full
campaign T1–T6 end-to-end:

1. **Preflight** — every binary kind used by the suites (plain-64b, plain, bse, PERT_KERNEL_DIAG) is
   selected via `petar.select`; a missing family **or an installed version different from the repository
   `VERSION`** triggers an automatic `./configure && make install` rebuild (T6's instrumented family is
   built last without install, since it reconfigures `CXXFLAGS`). Builds run sequentially before any
   simulation starts.
2. **Parallel suites** — suites run concurrently under a core budget (`--max-cores`, default
   `nproc`); each suite occupies the `cores` slots declared in its scenario JSON for its whole
   duration, so a 1-core suite starts whenever a core is free.
3. **HTML reports** — each suite finishes by regenerating its summary HTML; exit code 0 only when
   every suite's checks pass.

```bash
python3 test/validation/run_validation.py                # all suites, auto-install, parallel, HTML
python3 test/validation/run_validation.py --max-cores 8  # cap the parallel budget
python3 test/validation/run_validation.py --pipeline t4  # one suite end-to-end
python3 test/validation/run_validation.py --dry-run      # print all commands, build nothing
```

**Per-test core marking**: scenario JSONs carry `cores` (slots the suite occupies; its runs are also
executed up to `cores`-way in parallel when `parallel_runs` is `true` — all current suites have
independent single-threaded runs) and `parallel_runs`. Within a suite, jobs run sequentially (T1's
64b/non64b overlays share a workdir); across suites everything is parallel.

**Single suite**: `--pipeline t1` … `--pipeline t6` (same preflight/jobs/HTML steps for that suite).
`t5/t6_*_report.py --run` likewise still work standalone.

**Raw scenario mode** (JSON report, no per-suite HTML): pass one or more `--scenario <file>`
explicitly; `--jobs N` bounds run parallelism inside a scenario. The runner auto-detects `petar` via
`command -v petar` for this mode.

**Binary selection rule (all tests)**: tests never hardcode binary suffixes or use the default
`petar` symlink directly. They select by feature via `petar.select --require/--optional` and rebuild
automatically when stale (`report_common.select_or_install_petar`); the pipeline preflight applies
the same rule to every suite. The only exception is T6, which builds its instrumented binary from
source by construction.

**Exit codes**: report scripts (`t1`–`t6` `*_report.py`) exit 0 when all recorded checks pass
and 1 otherwise (also when re-rendering stored JSON reports — a failed check fails the script).
Chain regenerations with `;` instead of `&&` when a stored FAIL is expected.

---

## Scenarios

### T1: Kepler Orbit Crossing the Changeover Region
- Notebook-faithful binary IC: `m1=1`, `m2=10`, `semi=0.1 pc`, `e=0.9`, apocenter start; three families at the same `dt_soft` sweep (`2^[-13,-6]`), both 64b and non64b builds (report merges both)
- **7 checks organised by three purposes** (one-sided slope bounds; build-specific windows skip on the other build via `require_var_[not_]contains`):
  - **pure tree** (`r_out=0.005`, `r_peri > r_out`): KDKDK4 slope > 3.5 above the round-off floor — dt `1.22e-4–7.81e-3` for 64b, `9.77e-4–7.81e-3` for non64b
  - **pure Hermite** (`r_out=2.0`, `r_apo < r_in`): slope > 3.5 below the block-step threshold (dt `2.44e-4–4.88e-4`); slope ≤ 1 (weak dt dependence) above it
  - **mix / full crossing** (`r_out=0.15`, `r_peri < r_in < r_out < r_apo`): follows pure Hermite inside `r_in` at fine dt (ratio in [0.5, 2], 64b — non64b annulus floor sits above the Hermite level so it is skipped); slope > 3.5 at the tree-dominated coarse end; error between the two references at fine dt (≤ pure tree)
- The crossing-vs-references table (Q ratios per build) documents the follow-Hermite coincidence; commands are documented as one template + variable table

```bash
python3 test/validation/run_validation.py --pipeline t1
```

Output: `test/out/t1_binary_changeover_kdkdk4_summary.html`

### T2: Long-term Binary under Auto Changeover Radii
- Pure binary, >100 periods; no `-r`: the auto changeover tracks `dt_soft`, so the contiguous sweep (dt17→dt8, ten values, 2⁻¹⁷…2⁻⁸) naturally covers **pure-tree (dt17: `r_out < r_peri`) / mix (dt16–dt9) / pure-Hermite (dt8: `r_in > r_apo`)** regimes; per-run auto radii and regime are tabulated in the report
- **Two purposes, 4 checks**: (a) no secular energy-error drift — per-period error envelope stationary over ~113 orbits, half-run ratio ρ ≤ 2 checked once per regime; (b) the auto-derived changeover keeps the error level approximately consistent across `dt_soft` (span ≤ 300× over the 512× dt sweep; fixed radii would scale as dt⁴)
- Osculating-element drift (|Δa/a0|, |Δe| run-maxima) reported as supporting numbers in the run table

```bash
python3 test/validation/run_validation.py --pipeline t2
```

Output: `test/out/t2_binary_longterm_autochangeover_summary.html`

### T3: Binary AR Group Linking in the Hermite Channel (Census-Verified)
- Plain-family binary (no stellar evolution — the earlier bse IC's 1.45e-7 "AR error" was the BSE energy floor), force pinned in the hard channel (`r_in = 0.02 > r_apo`); census from whole-run group event dumps (`<prefix>.group.n<m>`)
- Three verified regimes: **never grouped** (`r_group=0.002`, pure Hermite pair, 4.7e-7 block-step truncation) / **stable group** (`r_group=0.018`, pure AR, never breaks, ~6e-12 — essentially exact) / **orbit-frequency Hermite↔AR switching** (`r_group=0.013`, criterion inside the separation range: ~1 break/re-form per orbit, measured 110/113; level between channels at ~1e4× AR but ~7× better than never grouping; stationary)
- 10 checks: 3 census + 1 anti-explosion guard (>2000 breaks flags the pre-fix 411k step-frequency chatter bug, fixed in PeTar 5765150 + SDAR 63e8308) + 3 level ratios + 3 stationarity

```bash
python3 test/validation/run_validation.py --pipeline t3
```

Output: `test/out/t3_binary_sdar_grouplink_summary.html`

### T4: Hierarchical Triple, Five-Mode Comparison
- Hierarchical triple IC (tidal-tensor 3-body test), evolved for 200 outer periods
- Five modes compared: pure SDAR reference, hard/tree with and without tidal tensor, and variants
- Automatic checks (10): reference quality — `m1` final **absolute** cumulative energy error `Error_cum` ≤ 1e-7 (the `Error/Total` ratio is not usable on this IC: |E_tot| ~1.7e-4 turns the BSE wind-accounting residual ~2.4e-10 absolute into ~1.4e-6 in ratio, the shared floor of every hard-channel mode); tidal-tensor census — artificial particles ≥ 1 with `--tt-switch 1`, = 0 with `--tt-switch 0`; runtime-abort regex counts per mode
- Control group `*_o15` runs are part of the same scenario: the five modes repeated with the original outer semi-major axis `a_out = 1.5`, distinguishing structural from scale-dependent effects (base-vs-control comparison table in the report)
- Requires a bse-interrupt binary (IC has `petar.init -s bse` columns, runs use `-b 1`); the pipeline selects it via `petar.select --require bse --optional mpi,omp,avx2` and builds it (`./configure --with-interrupt=bse && make install`) if missing

```bash
python3 test/validation/run_validation.py --pipeline t4
```

Output: `test/out/t4_triple_mode_summary.html`

### T5: Massive Binary + Light Churner, Tidal-Tensor Momentum Conservation
- IC from the 2026-10-06 momentum-leak event: 77.65 Msun binary (`a=0.0081` pc, `e=0.667`) plus a 0.242 Msun star whose pericenter (0.01 pc) lies inside the group radius and apocenter (0.04 pc) outside it, so AR group membership churns every outer orbit (the T3 ~ r^-4 close-passage regime)
- Three runs: pure hard mode with `--tt-switch 1/0` (the hard cluster is always the whole system — no tree perturbers — so however the AR group census splits as the churner crosses the linking radius (~10³ topology switches, verified from group event dumps), the TT channel stays structurally inactive and the TT-on run must match the TT-off baseline) plus tree mode with `--tt-switch 1` (production configuration of the event: TT engaged with group-boundary churn)
- Automatic checks: linear-momentum drift `max ||P|-|P0||/M` from `.status` readback (tree-mode threshold 0.005 chosen so the pre-fix binary fails at ~0.0103 while the fixed binary passes at ~0.0014), artificial-particle counts (TT engaged in tree mode, none in hard mode), runtime-abort regex counts
- Requires a plain-family binary (no interrupt/external): the IC is plain-format, so galpy- or BSE-enabled builds cannot read it; `t5_massive_binary_tt_report.py --run` selects it automatically via `petar.select --optional mpi,omp,avx2` (this changes the active binary selection)

```bash
python3 test/validation/run_validation.py --pipeline t5   # standalone legacy: t5_massive_binary_tt_report.py --run
```

Output: `test/out/t5_massive_binary_tt_summary.html` (momentum-drift figure, run summary, checks; JSON report at `test/out/report.t5.massive_binary_tt.json`)

### T7: AR Adaptive-Step Adjustment for a Perturbed Hierarchical Triple
- Minimal 4-body plain-format reproducer: an ultra-tight hierarchical trio (42.7 Msun binary, a = 8.34e-7 pc; barely-unbound 0.084 Msun intruder lingering at ~4e4 a_in with v ~ v_esc; exact c.m.-subtracted state at full double precision) plus one static 1.0 Msun perturber. With the scenario radii the intruder's changeover band is [0.082, 0.545]: the k(r)-weighted perturber force injects a continuous per-step energy flow into the internal-only H audit, and the SD reference's full-weight work credit made its residual R = S - etot_sd_ref_ drift at -(1-1/kappa)*W, collapsing ds ~5 orders (10^6-substep caps; pre-fix reference on this IC: AR steps 2.27e7 + warning). Root-fixed in SDAR 2026-10-21 (commit 6a8a357): kappa-weighted work credit at the slowdown system's own velocity gradient, n>=3 `fix_step_option=later` auto-grant removed, physically-gated sd-gauge repairs; the interim 2026-10-19 rate-type EMA separation was removed (current binary: ~1.25e4 steps)
- Seven runs: r1 perturber at 0.15 and r3 at 0.30 (intruder band — trigger guards: cumulative AR-substep budget `ar_step_sum ≤ 10⁶` calibrated 80× above the current LogH value and 23× below the pre-fix one, zero step-cap warnings / ds-floor aborts / runtime aborts, physical energy-error guard), r2 at 0.05 (deep inside r_in, k = 1 — false-positive guard), r4 = r1 with the BTLogH g-function (`--ar-g-func 1` — gauge comparison; ~1.15e5 steps, ~9× the LogH count on this IC, same budget + zero-warning guards), r5 at 0.60 (beyond r_out — formerly excluded for aborting via the `later` freeze, a gap removed by the same root fix), and r6/r7 at 0.35/0.50 across the former ds-floor-abort band (pre-fix NaN-abort zone; abort-to-complete guards). r5–r7 guard no step-cap warning / ds-floor abort / runtime abort. The T7 trigger itself runs with `Fix: none`, i.e. the adaptive path under guard is active
- 22 checks; half-tree-step window (`-t 0.0078`): the trigger develops within the first tree step, suite cost ~2 s; verified discriminative (scenario FAILS with the pre-fix binary)

```bash
python3 test/validation/run_validation.py --pipeline t7
```

Output: `test/out/t7_triple_perturber_ar_ds_adjustment_summary.html`

### T8: LogH vs BTLogH Gauge Comparison (Extracted Event Classes)
- Three event classes extracted (CM-recentered, full double precision) from the
  2026-10-22 N=1000 binary-rich cluster benchmark dump sweep (`bench-n1k-arfix/`),
  keeping the physically important members: r1/r2 **cascade** — the 105 Msun
  tight binary (a=6.5e-7 pc, e=0.016) + 0.3 Msun binary at mutual pericenter,
  PLUS the 7 bystanders carrying 95% of the cluster's tidal field on the group
  (minimal set by G*m/r^3; without them the slowdown stays engaged and both
  gauges integrate cheaply — the disengagement was externally driven); dump
  replay measured LogH 2.2M-11.4M AR steps/block vs BTLogH 12k-35k (38-934x).
  r3/r4 **chaotic** — the step-dominant 6-body group of a pre-fix pathological
  31-star cluster; dump replay: LogH 65.8k steps with 16 hard-energy events vs
  BTLogH 265k with 0. r5/r6 **control** — two detached binaries (parity guard).
- Six runs (each case in LogH and `--ar-g-func 1`), plain family, solver radii
  scaled to the IC (`--r 0.05 --r-search-group 0.005`; the hard manager asserts
  r_search_group <= r_in)
- 22 checks: zero runtime/ds-floor aborts everywhere; calibrated budgets
  (cascade LogH floor 3e3 guards the event actually triggers, ceiling 1e5;
  cascade BTLogH ceiling 6e3; cascade gauge ratio in [0.02, 0.5] — measured 0.11,
  dump-replay 38-934x — guarding the BTLogH hierarchical advantage; chaotic
  LogH floor 3e3 + ceilings 1e5; chaotic BTLogH zero hard-energy events; control
  ceilings 5e3 + zero events both gauges); suite cost ~3 s
- The N=1000 full-cluster BTLogH abort observed at suite creation (ds-denormal
  collapse, t=18.5 Myr) was root-caused and fixed 2026-10-22 (SDAR `96f0d39`,
  PeTar `d14ce63`; record in SDAR `docs/hierarchical_blogh_impl_notes.md`,
  2026-10-22 entry — plan item 5, closed); the full-cluster run now completes
  t=30 at energy parity with LogH. T8 remains a per-event gauge-signature guard

```bash
python3 test/validation/run_validation.py --pipeline t8
```

Output: `test/out/t8_gauge_logh_btlogh_extracted_summary.html`

### T9: Micro-TDE Trigger Criterion (Rastello et al. 2026)
- First validation of the micro-TDE module (commit a9978e1 + the delayed-interrupt precedence fix): BH 20.059 Msun (kw=14) + MS star 2.3934 Msun (Z=0.02, BSE ZAMS R=1.7543 Rsun, r_t=3.563 Rsun) with a 0.5 Msun bystander at 1 pc; bound grid (a=100 AU marked binary, apocentre start, `-b 1`) and hyperbolic grid (v_inf=5 km/s, 100 AU inbound start, group detection) sweep the pericentre through β = 0.3–5; NS (1.4 Msun, kw=13) covers the second compact-object branch; MS-MS control has no compact object
- 16 runs, 83 checks, suite cost ~10 s: detection boundary at β=1 (exactly one event for β≥1.1, zero for β≤0.9 and for MS-MS — measured β_logged matches β_design to ~7e-6); orbit-sign routing (semi<0 → `data .bse.hyperbolic_tde`, semi>0 → `.binary_tde`); disrupted star becomes kw=15 massless remnant; BH accretes `ftzacc`=0.1 of the stellar mass exactly (mix.f)
- **Perturbed-candidate pair (Rastello comm. 2026-10, premature-trigger regression)**: a wide high-e binary (a=380 AU, e=0.99997, r_p=2.38 Rsun, β=1.5, first pericentre at 781.7 yr) that is a TDE candidate by osculating elements from t=0 — the control fires at trigger+t_to_peri = the analytic pericentre epoch (timing deviation ~0.04 yr, threshold 1 yr); its twin with a 1 Msun flyby (pericentre 250 AU from the pair CM at ~400 yr) lifts r_p to 22–400 Rsun during the approach and must fire **nothing** (old versions activated the TDE before the pair reached the trigger region)
- **Author-event replays (micro_TDE reference)**: events 9 and 10 replayed from their last pre-event snapshot (aged stars via `petar.init --time`, clock aligned to the snapshot epoch through the `t9_replay_fix` patch: header time + t_record/t_interrupt) — logged trigger orbits reproduce the author values (event 9: a to 0.01%, r_p to 0.2%; event 10: r_p to 0.6%, a to 28% — phase-chaotic, looser band). Reference data tracked under `test/validation/data/micro_TDE/`. Survey of the full sample: events 1/3/4/5/7 do not re-trigger from the stripped few-body state (missing cluster field; event 1 additionally has an internal velocity/position inconsistency in the reference CSV — |dv| from the velocity columns is ~3–8× the finite-difference of the positions), events 2/8 exceed the suite runtime budget, event 6 gives instant-merger artifacts (a=37 Rsun, t_peri~0)
- Documented behavior: the mix.f accretion mod covers BH+MS only (`K1=14, K2<=1`) — an NS TDE disrupts the star with **zero** accretion (encoded as a regression check; prescription pending confirmation with Rastello et al.)
- Requires a bse-interrupt binary (`petar.select --require bse --optional mpi,omp,avx2`); event parsing lives in `run_validation.py` (`extract_tde_event_metrics`, wired by the scenario's `"tde_event": true` extract key); the scenario JSON is generated by the report script (`t9_tde_trigger_criterion_report.py --gen-scenario`; `--run` regenerates it first), keeping the grid as code in one place; 18 runs, 97 checks, suite cost ~13 s

```bash
python3 test/validation/run_validation.py --pipeline t9
```

Output: `test/out/t9_tde_trigger_criterion_summary.html`

### T6: AR Perturber-List Kernel Match (Members Keep Own-Mass Changeover)
- Design under test: group artificial particles (c.m./tidal-tensor/pseudo) carry the group c.m. changeover while AR group members keep their own-mass changeover; the tidal tensor is measured at the c.m. band but the direct perturber-list force on members uses the member band, so the two kernels differ in the member-r_out..c.m.-r_out transition annulus
- IC: engineered central massive triple (4+4+3 Msun, a_in=0.0012/a_out=0.008 pc) inside a compact Plummer field (297 stars, a=0.3 pc, masses 0.2-5 Msun, seed 20261010) — the triple guarantees a tidal-tensor-active AR group from the first tree step; group linking needs `--r-ratio 0.15 --r-group 0.008 --r-search-group 0.01` (r_search_group <= r_in constraint)
- All runs cover exactly one tree step (`-t 0.00390625` = dt_soft auto-derived from `-r 0.1`): round-off has no time to amplify chaotically, so r1/r2/r3 replicate the same trajectory — r1 measures the kernel mismatch (`PERT_KERNEL_DIFF=1`), r2/r3 are the member-kernel/c.m.-kernel A/B whose end-of-step energy-error difference is the pure in-step kernel effect (r1 and r2 end bit-identically, proving measurement neutrality)
- Checks (14): single step and TT-active group (>=12 artificial particles) in every run, diagnostic sample size (>=5000 evaluations), mismatch-ratio regression band (0.01-0.5; measured ~0.24 on this IC vs ~1.2e-3 aggregated at production N=1000), end-of-step |dE/E| < 1e-3 for r2/r3, zero large-energy dumps, zero aborts
- Requires a plain-family binary built with `-D PERT_KERNEL_DIAG` (portable: `CXXFLAGS='-D PERT_KERNEL_DIAG' ./configure && make -j` — the Makefile adds `-D PROFILE` itself and `make` without a target builds the configured family); `t6_pert_kernel_match_report.py --run` configures and builds it automatically when missing (this changes the configure state)

```bash
python3 test/validation/run_validation.py --pipeline t6   # standalone legacy: t6_pert_kernel_match_report.py --run
```

Output: `test/out/t6_pert_kernel_match_summary.html` (paper-style sections: background, method with metric definitions, single-step A/B results, checks; JSON report at `test/out/report.t6.pert_kernel_match.json`)

---

## Files

| File | Purpose |
|------|---------|
| `run_validation.py` | Scenario runner, metric extractor, pass/fail checker, `--pipeline t1/t2/t3` driver |
| `report_common.py` | Shared report helpers (parsing, escaping, plotting, cards); see REPORT_DESIGN.md §5 |
| `REPORT_DESIGN.md` | Design rules for validation scenarios and HTML reports (T6 = reference) |
| `metrics.py` | Common check functions |
| `make_ic.py` | Deterministic IC generator (`--case t1…t6`) |
| `criteria.json` | Centralized threshold reference |
| `t1_binary_changeover_kdkdk4_report.py` | T1 report generation |
| `t2_binary_longterm_autochangeover_report.py` | T2 report generation |
| `t3_binary_sdar_grouplink_report.py` | T3 report generation |
| `t4_triple_mode_report.py` | T4 report generation |
| `t5_massive_binary_tt_report.py` | T5 report generation (`--run` executes the scenario first) |
| `t6_pert_kernel_match_report.py` | T6 report generation (`--run` executes the scenario first; reuses T5 SVG helpers) |
| `scenarios/*.json` | Scenario definitions (commands, runs, checks) |

## Outputs

- Run logs, metric JSON, and HTML reports → `test/out/` (`--out-dir` / `--report` defaults; pass `--out-dir test/validation/out` to keep them under the validation directory)
- Temporary work directories → `test/validation/work/` (generated at runtime, not tracked)

## Dry Run

Verify scenario expansion without running simulations:

```bash
python3 test/validation/run_validation.py --dry-run
```

`--dry-run` applies to scenario mode; the `--pipeline` driver ignores it and runs for real.

## Advanced

Override auto-detected binaries (scenario mode only):

```bash
python3 test/validation/run_validation.py \
  --var petar_bin_switch=/custom/path/to/petar
```

> `--pipeline t1…t6` does **not** accept `--var` (and ignores `--dry-run`): it selects or builds the binaries it needs itself, then runs the scenarios and the report generator. To use custom binaries, run the scenario directly with `--scenario` + `--var`.

Continue on scenario failure:

```bash
python3 test/validation/run_validation.py \
  --scenario test/validation/scenarios/t1_binary_changeover_kdkdk4.json \
  --continue-on-error
```