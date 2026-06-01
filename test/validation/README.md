# Validation Framework (T1-T4 + blogh placeholder)

This directory provides a minimal, scriptable validation framework for major PeTar algorithm changes.

Current scope:

- `T1` notebook-faithful KDKDK4 high-eccentricity binary test:
  - exact IC from `PeTar_KDKDK4_method.ipynb`: `m1=1`, `m2=10`, `semi=0.1 pc`, `e=0.9`, start at apocenter
  - corrected inside regime for particle-tree behavior: `rout < peri`
  - exact tree-step grid from notebook: `dt_soft = 2^[-13,-6]`
  - exact representative `rout` values from notebook sweep:
    - `rout=0.005` (`rout < peri=0.01`): expect particle-tree dominated 4th-order trend
    - `rout=0.32`, `1.28`: crossing regime representatives, expect small-dt plateau behavior
- `T2` long-term pure-binary conservation with KDKDK4:
  - deterministic 2-body IC
  - integration duration >= 100 periods
  - dt_soft sweep with PeTar auto changeover radii (no manual `-r`)
  - focus on long-term `max |Δa/a0|` and `max |Δe|`
  - do not assume strict error scaling with `dt_soft` under auto changeover
  - automated checks include 100-period growth-law proxies from logs:
    - `sum_abs_error_period_max` (sum of per-period maxima of `|Error/Total|`)
    - `dt_soft_times_max_abs_error_period_max` (`dt_soft *` per-period max error)
    - both are checked for nondecreasing trend from fine to coarse `dt_soft`
  - non64b binary as default execution target
- `T3` long-term binary hard-switching test:
  - fixed `r_out` and fixed `dt_soft`
  - choose `r_out` such that `r_in=0.1*r_out > r_apo`, so force calculation stays in hard region
  - compare three `r_group/r-search-group` regimes for Hermite-SDAR switching:
    - `r_group > r_apo` (outside apocenter)
    - `r_group < r_peri` (inside pericenter)
    - `r_peri < r_group < r_apo` (between)
  - report long-term `max |Δa/a0|`, `max |Δe|`, and 100-period error accumulation proxies
  - automated checks enforce ordering across the three regimes (outside apo -> between -> inside peri):
    - `max_rel_drift_semi` nondecreasing
    - `max_abs_drift_ecc` nondecreasing
    - `sum_abs_error_period_max` nondecreasing
- `T4` notebook-scale hierarchical 3-body mode comparison:
  - notebook-inspired triple scales from `PeTar_function_test.ipynb` tidal tensor 3-body section
  - five modes in one scenario:
    - `m1_pure_sdar_ref`: pure SDAR reference
    - `m2_hard_no_tt`: inner SDAR + outer Hermite, no tidal tensor
    - `m3_hard_tt`: inner SDAR + outer Hermite, with tidal tensor
    - `m4_tree_no_tt`: inner SDAR + outer particle-tree, no tidal tensor
    - `m5_tree_tt`: inner SDAR + outer particle-tree, with tidal tensor
  - explicit tidal-tensor usage check via runtime counts: require `N_all(glb) - N_real(glb) > 0` when `--tt-switch=1`
- `T3-blogh` A/B placeholder (same IC, two configurable binaries for future blogh integration)

## Files

- `run_validation.py`: scenario runner, metric extractor, and pass/fail checker.
- `metrics.py`: common check functions.
- `t1_kdkdk4_report.py`: generate one-file HTML report with plots, IC table, run parameters and interpretation.
- `scenarios/t1_high_ecc_changeover.json`: T1 scenario definition.
- `scenarios/t2_binary_conservation_longterm.json`: T2 long-term binary conservation scenario definition.
- `scenarios/t3_binary_hard_switch_longterm.json`: T3 scenario definition.
- `scenarios/t4_tree_hard_from_triple.json`: T4 scenario definition.
- `t4_triple_mode_report.py`: T4 one-file HTML report for five-mode comparison.
- `scenarios/t3_blogh_ab_placeholder.json`: blogh A/B placeholder scenario.
- `criteria.json`: centralized threshold reference (now read by `run_validation.py`).

## 1) Dry-run first

From repository root:

```bash
python3 test/validation/run_validation.py --dry-run \
  --var petar_bin_order2=petar \
  --var petar_bin_order4=petar \
  --var petar_bin_switch=petar
```

This verifies scenario expansion and command generation without starting simulations.

## 2) Select binaries for order comparison

Scenario defaults now include concrete commands and deterministic initial condition generation.

To compare 2nd/4th order behavior in T1, set binaries explicitly based on your local installation, for example:

- `petar_bin_order2`: a KDK-style binary (example from notebook: `petar.mpi.omp.avx2.64b`)
- `petar_bin_order4`: a KDKDK4-style binary (example from notebook: `petar.mpi.omp.avx2.64b.kdkdk4`)
- `petar_bin_switch`: binary for T2 long-term binary conservation checks

Example run:

```bash
python3 test/validation/run_validation.py \
  --var petar_bin_order2=petar.mpi.omp.avx2.64b \
  --var petar_bin_order4=petar.mpi.omp.avx2.64b.kdkdk4 \
  --var petar_bin_switch=petar.mpi.omp.avx2.64b.kdkdk4
```

Run only one scenario:

```bash
python3 test/validation/run_validation.py \
  --scenario test/validation/scenarios/t3_binary_hard_switch_longterm.json
```

Run blogh placeholder A/B with explicit binaries:

```bash
python3 test/validation/run_validation.py \
  --scenario test/validation/scenarios/t3_blogh_ab_placeholder.json \
  --var petar_bin_blogh_a=petar.mpi.omp.avx2.64b.kdkdk4 \
  --var petar_bin_blogh_b=petar.mpi.omp.avx2.64b
```

## 3) Outputs

- Run logs: `test/out/validation_<scenario_name>/*.log`
- Final report: `test/out/report.json`

The report includes:

- extracted metrics (`max_abs_error_over_total`, `max_abs_error_pp`, `regex_counts`)
- check records
- pass/fail summary

Notes:

- `T1` uses the notebook-faithful deterministic 2-body setup generated by `make_ic.py --case t1`.
- `T2` uses a deterministic pure-binary setup generated by `make_ic.py --case t2`.
- `T3` uses a deterministic pure-binary setup generated by `make_ic.py --case t2`.
- `T4` uses a deterministic triple + background setup generated by `make_ic.py --case t4`.
- `T4` now uses a deterministic notebook-scale hierarchical triple generated by `make_ic.py --case t4` (3 real particles; no background particles).
- `T3-blogh` reuses deterministic `t3` IC and provides A/B placeholders for future blogh-vs-reference comparison.
- Both setups are in `Msun, pc, pc/Myr` and converted by `petar.init` into PeTar input snapshots.
- T1 convergence-ratio checks are auto-skipped when `petar_bin_order2` and `petar_bin_order4` are identical or when the fine-step error is below a numeric floor.
- Thresholds are resolved from `criteria.json` first; per-scenario inline values are fallback defaults.

Generate T1 HTML report (after running T1 scenario):

```bash
python3 test/validation/t1_kdkdk4_report.py \
  --report test/out/report.t1.json \
  --output test/out/t1_kdkdk4_report.html
```

One-command T1 pipeline (64b + non64b + merged summary HTML):

```bash
test/validation/t1_kdkdk4_pipeline.sh \
  /home/lwang/bin/petar.mpi.omp.avx2.64b.kdkdk4 \
  /home/lwang/bin/petar.mpi.omp.avx2.kdkdk4
```

This generates:

- `test/out/report.t1.64b.overlay.json`
- `test/out/report.t1.non64b.overlay.json`
- `test/out/t1_kdkdk4_summary.html`

Generate T2 HTML report (after running T2 scenario):

```bash
python3 test/validation/t2_binary_conservation_report.py \
  --report test/out/report.t2.binary.json \
  --output test/out/t2_binary_conservation_summary.html
```

One-command T2 pipeline (non64b long-term binary run + HTML):

```bash
test/validation/t2_binary_pipeline.sh \
  /home/lwang/bin/petar.mpi.omp.avx2.kdkdk4
```

This generates:

- `test/out/report.t2.binary.json`
- `test/out/t2_binary_conservation_summary.html`

Generate T3 HTML report (after running T3 scenario):

```bash
python3 test/validation/t3_binary_hard_switch_report.py \
  --report test/out/report.t3.binary.json \
  --output test/out/t3_binary_hard_switch_summary.html
```

One-command T3 pipeline (non64b long-term three-regime run + HTML):

```bash
test/validation/t3_binary_pipeline.sh \
  /home/lwang/bin/petar.mpi.omp.avx2.kdkdk4
```

This generates:

- `test/out/report.t3.binary.json`
- `test/out/t3_binary_hard_switch_summary.html`

Generate T4 HTML report (after running T4 scenario):

```bash
python3 test/validation/t4_triple_mode_report.py \
  --report test/out/report.t4.triple.json \
  --output test/out/t4_triple_mode_summary.html
```

One-command T4 pipeline (five-mode triple comparison + tidal tensor activation checks + HTML):

```bash
test/validation/t4_triple_pipeline.sh \
  /home/lwang/bin/petar.mpi.omp.avx2.kdkdk4
```

This generates:

- `test/out/report.t4.triple.json`
- `test/out/t4_triple_mode_summary.html`

## 4) Next extension targets

- Replace blogh placeholder binary mapping with real blogh-enabled runtime mode
- Add statistical regression mode for multi-seed runs
- Add explicit A/B delta metrics (relative error ratio between methods)
