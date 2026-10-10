# Validation Framework (T1–T6)

**Design rules**: new validation reports follow [REPORT_DESIGN.md](REPORT_DESIGN.md) (paper-style
sections, symbol definitions, figure encoding, threshold semantics, test-design principles, shared
`report_common.py` helpers). T6 is the reference implementation; T1–T5 are legacy style pending
migration.

Scriptable validation scenarios for PeTar algorithm changes.

## Quick Start

**Prerequisites**: Run `make install` or `petar.select <suffix>` so `petar` is in PATH.

Run any full pipeline (binary selection → scenario → HTML report, all auto):

```bash
python3 test/validation/run_validation.py --pipeline t1   # 64b vs non64b KDKDK4
python3 test/validation/run_validation.py --pipeline t2   # long-term binary conservation
python3 test/validation/run_validation.py --pipeline t3   # Hermite-SDAR switching
```

T4 has its own wrapper (it is not a `--pipeline` choice), taking an optional solver binary argument:

```bash
bash test/validation/t4_triple_pipeline.sh [path/to/petar]
```

Run all scenarios only (JSON report, no HTML):

```bash
python3 test/validation/run_validation.py
```

The runner auto-detects `petar` via `command -v petar`.

---

## Scenarios

### T1: KDKDK4 Changeover (64b vs non64b)
- Notebook-faithful binary IC: `m1=1`, `m2=10`, `semi=0.1 pc`, `e=0.9`, apocenter start
- Sweeps `dt_soft = 2^[-13,-6]` across inside, peri-between, and crossing `rout` regimes
- HTML report with overlay comparison of 64b vs non64b KDKDK4

```bash
python3 test/validation/run_validation.py --pipeline t1
```

Output: `test/out/t1_kdkdk4_summary.html`

### T2: Long-term Binary Conservation
- Pure binary, ≥100 orbital periods, `dt_soft` sweep, auto changeover radii (no `-r`)
- Tracks `max |Δa/a0|`, `max |Δe|`, and period-binned error proxies
- Automatic checks: monotonic error growth from fine to coarse `dt_soft`

```bash
python3 test/validation/run_validation.py --pipeline t2
```

Output: `test/out/t2_binary_conservation_summary.html`

### T3: Binary Hermite–SDAR Switching
- Fixed `r_out` and `dt_soft`, force stays in hard region (`r_in > r_apo`)
- Three `r_group` regimes: outside apocenter ↔ between peri/apo ↔ inside pericenter
- Automatic checks: ordering of semi/ecc drift and cumulative error across regimes

```bash
python3 test/validation/run_validation.py --pipeline t3
```

Output: `test/out/t3_binary_hard_switch_summary.html`

### T4: Hierarchical Triple, Five-Mode Comparison
- Hierarchical triple IC (tidal-tensor 3-body test), evolved for 200 outer periods
- Five modes compared: pure SDAR reference, hard/tree with and without tidal tensor, and variants
- Automatic checks: energy-error ordering across modes, inner/outer orbital-element conservation, runtime-abort regex counts
- Control scenario `t4_tree_hard_from_triple_outer15.json` repeats the comparison with the original outer semi-major axis `a_out = 1.5`
- Requires a bse-interrupt binary (IC has `petar.init -s bse` columns, runs use `-b 1`); the pipeline selects it via `petar.select --require bse --optional mpi,omp,avx2` and builds it (`./configure --with-interrupt=bse && make install`) if missing

```bash
bash test/validation/t4_triple_pipeline.sh [path/to/petar]
```

Output: `test/out/t4_triple_mode_summary.html`

### T5: Massive Binary + Light Churner, Tidal-Tensor Momentum Conservation
- IC from the 2026-10-06 momentum-leak event: 77.65 Msun binary (`a=0.0081` pc, `e=0.667`) plus a 0.242 Msun star whose pericenter (0.01 pc) lies inside the group radius and apocenter (0.04 pc) outside it, so AR group membership churns every outer orbit (the T3 ~ r^-4 close-passage regime)
- Three runs: pure hard mode with `--tt-switch 1/0` (with 3 particles the churner orbit necessarily overlaps the group-linking radius, so the system stays one isolated AR group — TT structurally inactive, and the TT-on run must match the TT-off baseline) plus tree mode with `--tt-switch 1` (production configuration of the event: TT engaged with group-boundary churn)
- Automatic checks: linear-momentum drift `max ||P|-|P0||/M` from `.status` readback (tree-mode threshold 0.005 chosen so the pre-fix binary fails at ~0.0103 while the fixed binary passes at ~0.0014), artificial-particle counts (TT engaged in tree mode, none in hard mode), runtime-abort regex counts
- Requires a plain-family binary (no interrupt/external): the IC is plain-format, so galpy- or BSE-enabled builds cannot read it; `t5_massive_binary_tt_report.py --run` selects it automatically via `petar.select --optional mpi,omp,avx2` (this changes the active binary selection)

```bash
python3 test/validation/t5_massive_binary_tt_report.py --run
```

Output: `test/out/t5_massive_binary_tt_summary.html` (momentum-drift figure, run summary, checks; JSON report at `test/out/report.t5.massive_binary_tt.json`)

### T6: AR Perturber-List Kernel Match (Members Keep Own-Mass Changeover)
- Design under test: group artificial particles (c.m./tidal-tensor/pseudo) carry the group c.m. changeover while AR group members keep their own-mass changeover; the tidal tensor is measured at the c.m. band but the direct perturber-list force on members uses the member band, so the two kernels differ in the member-r_out..c.m.-r_out transition annulus
- IC: engineered central massive triple (4+4+3 Msun, a_in=0.0012/a_out=0.008 pc) inside a compact Plummer field (297 stars, a=0.3 pc, masses 0.2-5 Msun, seed 20261010) — the triple guarantees a tidal-tensor-active AR group from the first tree step; group linking needs `--r-ratio 0.15 --r-group 0.008 --r-search-group 0.01` (r_search_group <= r_in constraint)
- All runs cover exactly one tree step (`-t 0.00390625` = dt_soft auto-derived from `-r 0.1`): round-off has no time to amplify chaotically, so r1/r2/r3 replicate the same trajectory — r1 measures the kernel mismatch (`PERT_KERNEL_DIFF=1`), r2/r3 are the member-kernel/c.m.-kernel A/B whose end-of-step energy-error difference is the pure in-step kernel effect (r1 and r2 end bit-identically, proving measurement neutrality)
- Checks (14): single step and TT-active group (>=12 artificial particles) in every run, diagnostic sample size (>=5000 evaluations), mismatch-ratio regression band (0.01-0.5; measured ~0.24 on this IC vs ~1.2e-3 aggregated at production N=1000), end-of-step |dE/E| < 1e-3 for r2/r3, zero large-energy dumps, zero aborts
- Requires a plain-family binary built with `-D PERT_KERNEL_DIAG` (portable: `CXXFLAGS='-D PERT_KERNEL_DIAG' ./configure && make -j` — the Makefile adds `-D PROFILE` itself and `make` without a target builds the configured family); `t6_pert_kernel_match_report.py --run` configures and builds it automatically when missing (this changes the configure state)

```bash
python3 test/validation/t6_pert_kernel_match_report.py --run
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
| `t1_kdkdk4_report.py` | T1 report generation |
| `t2_binary_conservation_report.py` | T2 report generation |
| `t3_binary_hard_switch_report.py` | T3 report generation |
| `t4_triple_mode_report.py` / `t4_triple_pipeline.sh` | T4 report generation & pipeline |
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

> `--pipeline t1/t2/t3` does **not** accept `--var` (and ignores `--dry-run`): it selects or builds the binaries it needs itself, then runs the scenarios and the report generator. To use custom binaries, run the scenario directly with `--scenario` + `--var`.

Continue on scenario failure:

```bash
python3 test/validation/run_validation.py \
  --scenario test/validation/scenarios/t1_high_ecc_changeover.json \
  --continue-on-error
```

Override the T4 pipeline binary by passing its path as an argument:

```bash
bash test/validation/t4_triple_pipeline.sh /path/to/petar
```