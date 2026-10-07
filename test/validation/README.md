# Validation Framework (T1–T5)

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

---

## Files

| File | Purpose |
|------|---------|
| `run_validation.py` | Scenario runner, metric extractor, pass/fail checker, `--pipeline t1/t2/t3` driver |
| `metrics.py` | Common check functions |
| `make_ic.py` | Deterministic IC generator (`--case t1…t5`) |
| `criteria.json` | Centralized threshold reference |
| `t1_kdkdk4_report.py` | T1 report generation |
| `t2_binary_conservation_report.py` | T2 report generation |
| `t3_binary_hard_switch_report.py` | T3 report generation |
| `t4_triple_mode_report.py` / `t4_triple_pipeline.sh` | T4 report generation & pipeline |
| `t5_massive_binary_tt_report.py` | T5 report generation (`--run` executes the scenario first) |
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