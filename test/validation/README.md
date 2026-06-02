# Validation Framework (T1–T3)

Scriptable validation scenarios for PeTar algorithm changes.

## Quick Start

**Prerequisites**: Run `make install` or `petar.select <suffix>` so `petar` is in PATH.

Run any full pipeline (binary selection → scenario → HTML report, all auto):

```bash
python3 test/validation/run_validation.py --pipeline t1   # 64b vs non64b KDKDK4
python3 test/validation/run_validation.py --pipeline t2   # long-term binary conservation
python3 test/validation/run_validation.py --pipeline t3   # Hermite-SDAR switching
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
bash test/validation/t1_kdkdk4_pipeline.sh
```

Output: `test/validation/out/t1_kdkdk4_summary.html`

### T2: Long-term Binary Conservation
- Pure binary, ≥100 orbital periods, `dt_soft` sweep, auto changeover radii (no `-r`)
- Tracks `max |Δa/a0|`, `max |Δe|`, and period-binned error proxies
- Automatic checks: monotonic error growth from fine to coarse `dt_soft`

```bash
bash test/validation/t2_binary_pipeline.sh
```

Output: `test/validation/out/t2_binary_conservation_summary.html`

### T3: Binary Hermite–SDAR Switching
- Fixed `r_out` and `dt_soft`, force stays in hard region (`r_in > r_apo`)
- Three `r_group` regimes: outside apocenter ↔ between peri/apo ↔ inside pericenter
- Automatic checks: ordering of semi/ecc drift and cumulative error across regimes

```bash
bash test/validation/t3_binary_pipeline.sh
```

Output: `test/validation/out/t3_binary_hard_switch_summary.html`

---

## Files

| File | Purpose |
|------|---------|
| `run_validation.py` | Scenario runner, metric extractor, pass/fail checker |
| `metrics.py` | Common check functions |
| `make_ic.py` | Deterministic IC generator |
| `criteria.json` | Centralized threshold reference |
| `t1_kdkdk4_report.py` / `t1_kdkdk4_pipeline.sh` | T1 report generation & pipeline |
| `t2_binary_conservation_report.py` / `t2_binary_pipeline.sh` | T2 report generation & pipeline |
| `t3_binary_hard_switch_report.py` / `t3_binary_pipeline.sh` | T3 report generation & pipeline |

| `scenarios/*.json` | Scenario definitions (commands, runs, checks) |

## Outputs

- Run logs, metric JSON, and HTML reports → `test/validation/out/`
- Temporary work directories → `test/validation/work/` (generated at runtime, not tracked)

## Dry Run

Verify scenario expansion without running simulations:

```bash
python3 test/validation/run_validation.py --dry-run
```

## Advanced

Override auto-detected binaries:

```bash
python3 test/validation/run_validation.py \
  --var petar_bin_switch=/custom/path/to/petar
```

Continue on scenario failure:

```bash
python3 test/validation/run_validation.py \
  --scenario test/validation/scenarios/t1_high_ecc_changeover.json \
  --continue-on-error
```

Override pipeline binaries by passing paths as arguments:

```bash
bash test/validation/t1_kdkdk4_pipeline.sh /path/to/64b/petar /path/to/non64b/petar
bash test/validation/t2_binary_pipeline.sh /path/to/petar
bash test/validation/t3_binary_pipeline.sh /path/to/petar
