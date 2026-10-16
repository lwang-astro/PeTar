# Validation Framework (T1–T6)

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