---
name: petar-nbody-simulation
description: "Use when: setting up or running PeTar N-body simulations, including petar.select binary-family switching, petar.init conversion, isolated clusters, BSE/SSE, Galpy or Agama external potential, MPI/OpenMP/GPU launch, restart/resume, petar.data post-processing, timestep tuning with petar.find.dt, snapshot extraction with petar.get.object.snap, format conversion with petar.format.transfer.post, galev processing, and movie generation."
---

# PeTar N-body Simulation Skill (Compact)

## Purpose

Provide strict, command-level guidance for PeTar workflows in this repository.
This compact version prioritizes execution safety, option correctness, and reproducibility.

## Non-Negotiable Rules

### Gate 1 — Selection and inputs

- **Infer required physics/runtime features first, then switch binary family with `petar.select`.** Never run whatever `petar` currently points to, and never switch binaries by hand-editing symlinks. See "Binary Selection and Capability Validation".
- **Ask only for missing required inputs; never guess a physics-defining parameter.** Do not proceed until every scenario-required input is resolved.
- Only two defaults may be applied without asking: output prefix `data`, unit mode `-u 1`. Everything physics-defining must be confirmed — see "Star-Cluster Generation Inputs" for the definitive checklist.
- Mandatory scenario confirmations: stellar-evolution module (`merger` | `bse` | `bseEmp` | `mobse` | `dsm`); external-potential mode (`galpy` | `agama`); metallicity for BSE-family runs; and, for Galpy/Agama, COM position **and** velocity in physical units as numerical values — never accept vague descriptions such as "Sun position".
- Exclude debug and GPU families unless the user explicitly asks for them.
- Enforce unit consistency across IC generation, `petar.init`, runtime options, and post-processing.

### Gate 2 — Working directory

- **Always `cd` into the designated simulation directory before any IC generation, solver, or post-processing command, and create no files outside it.** Test runs, debugging, and exploratory commands go there too — never the repository root or another project directory. Create the directory first if it does not exist.
- **If the user does not specify a working directory, ask — do not decide the path yourself.** A run produces ICs, logs, snapshots, parameter dumps, and movies. If the user has no preference, suggest a convention such as `~/petar_sim/<descriptive-name>` and confirm before creating.

### Gate 3 — Confirmation and traceability

- **Before any solver invocation (`petar`, `petar.find.dt`, `petar.bse`, …), present a structured parameter summary and wait for explicit confirmation.** This applies to "quick tests" and "just checking" too — a solver call costs at least seconds and may write files.
  Summary contents, in this order:
  1. Working directory and input snapshot
  2. Full command line (binary, flags, options, redirects)
  3. Estimated runtime hint (e.g., "N=1000 × 100 Myr → a few minutes")
  4. Expected output files (e.g., `output`, `data.*`, `commands.log`)
  5. A clear prompt: *"Proceed? (y/n)"* — stop and wait for the response.
- **Record every command to `commands.log`** in the working directory before executing each major step (IC generation, `petar.init`, `petar.find.dt`, solver, post-processing, movie), with a timestamp, so the exact sequence is reproducible:
  ```bash
  echo "# $(date): petar -u 1 --galpy-set MWPotential2014 -t 100.0 -o 1.0 -s 0.0009765625 input" >> commands.log
  ```
  Append when resuming or adding steps; create a new file at the start of a fresh simulation.
- **Redirect stdout + stderr of every major command to a file** — lost output cannot be reviewed later. Consistent naming: `mcluster >mc.log`, `petar &>output`, `petar.find.dt &>finddt.log`, `petar.data.process &>process.log`, `petar.movie &>movie.log`, `petar.external.* &>ext.log`, `petar.data.gether &>gether.log`. (`petar.init` output is minimal — terminal is fine.)

### Gate 4 — Build and binary integrity

- **Check binary version against source before running a binary**, and ask whether to rebuild (`./configure` + `make install`) on mismatch:
  ```bash
  binary_version=$(strings <binary> | grep -A1 "^Version: $" | tail -1)
  source_version=$(echo "$(cat VERSION)_$(cat ../SDAR/VERSION)")
  ```
- **When rebuilding, pass only the flags the current scenario needs** — never reuse the previous `config.status` command verbatim; it may carry unneeded features (e.g. `gasdrag`, `pnhermite`) into an unnecessarily large binary. Derive flags from the features selected for this scenario:
  - `--with-interrupt=<bse|mobse|bseEmp|dsm>` for stellar evolution / disk star merger
  - `--with-external=<galpy|agama>` for external potential
  - `--with-external-hard=<gasdrag>` — **only if gas drag is explicitly requested**
  - `--with-pn=<pnall|pnhermite|…>` — **only if post-Newtonian correction is explicitly requested**
  - `--enable-mpfrc` — **only if MPFRC precision is explicitly requested**
  - `--enable-64b`, `--enable-avx2`, `--enable-avx512`, `--enable-omp`, `--enable-mpi` — architecture options, minimal feature risk
- `*.hard.debug`, `*.format.transfer`, `petar.hard.test`, and `*.dump2test` are **not** production solvers. Source-level debugging requires a rebuild with `--with-debug=g`.
- **Validate every custom option against the selected binary's `-h` output before emitting a runnable command.**

### Gate 5 — Timestep

- **Run `petar.find.dt` before every production run, then halve the recommendation and pass it as `-s <value>`.** Overhead is negligible (seconds to a minute); speedup is 2–5×. Exception: explicit quick/debug runs (`-t 0`, visibly sub-minute runtime) — skip it and let the solver auto-estimate. See "Timestep Tuning" below.

### Workflow integrity

- Keep fresh IC generation separate from restart/resume paths. **Never use `petar.init` in a restart workflow.**
- Modern PeTar (tmp-file mechanism) generates `data.snap.lst` at runtime and commits stream outputs automatically — **no gathering step is needed**. `petar.data.gether` is legacy, required only for old MPI runs lacking `snap.lst`; when it is required, gather before downstream processing.
- Prefer installed PeTar tools over abstract advice whenever the request maps to one (analysis, conversion, movie generation, restart cleanup, external-potential maps).
- For post-processing tools, mode flags must match the producing solver family and snapshot format.
- For `petar.data.process` in smoke or repeated-run scenarios, use `--no-auto-resume` to prevent unintended auto-restart against stale outputs.
- For functional smoke automation, required test inputs must be git-tracked; output artifacts are not inputs.
- **After any `astropy` upgrade** (or when adding a constant to `astro_units.hpp`), regenerate the physical constants header with `python tools/generate_astro_units.py` so `G_ASTRO`, `PCMYR_TO_KMS`, `SPEED_OF_LIGHT`, etc. stay aligned with the latest IAU/CODATA definitions. See `README.md` (Compilation and Installation → Physical Constants).
- When modifying source headers that define CLI options (`src/petar.hpp`, `src/hard.hpp`, `bse-interface/*.h`, `galpy-interface/*.h`, `agama-interface/*.h`, `src/disk_star_merger.hpp`, `src/gas_drag.hpp`, `src/external_hard.hpp`, `parallel-random/rand_io.hpp`), re-run `assets/generate_option_reference.py` so `option-reference.md` stays synchronized.

### Environment requirements

- Before any PeTar launch, require `export OMP_STACKSIZE=128M` (and `ulimit -s unlimited` on Linux). Missing this causes segfault on large simulations.
- For MPI+OpenMP launches, recommend `--bind-to none` to prevent processes being restricted to a single core.
- For MPI launches in UCX environments, require `export UCX_VFS_ENABLE=n` to avoid a known UCX segfault.
- When `./configure` runs on a cluster login node, warn that auto-detected SIMD capabilities may differ from compute nodes; suggest verifying with `--enable-avx2` / `--enable-avx512` flags.

### bseEmp-specific

- `bseEmp` requires manually linking `ffbonn` or `ffgeneva` metal-poor track directories before first use. Without this, the binary crashes at initialization with a file-not-found error. Ask the user to set up these links before composing a bseEmp run command.

### DSM-specific

- DSM init uses `petar.init -s dsm --type <type> --radius <radius>` where `--radius` is the disk outer edge in simulation units (pc for `-u 1`).
- DSM runtime requires `--detect-interrupt 1 --dsm-new-star-mode 0` for stable smoke behavior.
- DSM post-processing and snap reads use `-i dsm` or `interrupt_mode="dsm"`.
- The file `data.interrupt` records interrupted DSMs; its binary layout may include trailing padding bytes — validate by record count, not strict dtype alignment.

### Coordinate origin

- Warn when the density center is far from `(0,0,0)` or includes very distant particles — this can trigger `n_jp<=pg.NJMAX` assertion failures due to tree cell resolution limits.

## Required Input Checklist

Collect the minimum required fields before composing run commands.

### Common

1. Working directory.
2. Scenario: isolated | BSE/SSE | DSM | Galpy | Agama | restart | standalone-bse.
3. End time `-t`.
4. Output interval `-o`.
5. **Launch mode** — confirm each sub-field explicitly (do not default without asking):
   - Mode: serial | OpenMP | MPI+OpenMP | GPU.
   - MPI process count (if MPI mode).
   - OpenMP thread count (if OpenMP/MPI+OpenMP mode) — scale with N, see "Parallel Launch Heuristics" in `assets/script-tools.md`; for N ≲ 10³ recommend 1 thread and say why.
   - Custom launcher prefix (e.g. `srun -N 2`), if any — otherwise use `mpiexec -n` for MPI.
6. Initial data source: existing PeTar snapshot | raw table | generator.
7. Extra user options to pass through (validated only).

### Scenario-Specific

- **BSE/mobse/bseEmp**: primordial binary count `-b`; metallicity (`--bse-metallicity` or family-equivalent).
- **DSM**: only the DSM controls the user asks to override (`--dsm-*`).
- **Galpy**: COM position and velocity in physical units; one of `--galpy-set` / `--galpy-conf-file` / `--galpy-type-arg`; optional `--galpy-rscale`, `--galpy-vscale`. PeTar supports Galpy ≤ 1.10.2 only.
- **Agama**: COM position and velocity in physical units; `--agama-conf-file`; optional `--agama-rscale`, `--agama-vscale`.
- **Restart**: restart snapshot; parameter file (`-p`, default `data.par`); append/overwrite behaviour (`-a 1` / `-a 0`). Companion files and `-i` format: see "Restart/resume" below.
- **standalone-bse**: use `petar.bse` / `petar.mobse` / `petar.bseEmp` directly, not the N-body solver — no IC generation and no `petar.init`. See `README.md`.

Per-scenario ask-minimum lists, including the "do not ask again if already provided" rule: `assets/minimal-question-sets.md`.

## IC and Unit Rules

### Generator Availability

When generator-based star-cluster IC creation is required:

1. Prefer `mcluster_gpu` if available.
2. Else use `mcluster` if available.
3. If neither exists, stop and ask user whether to proceed with installation guidance.

### mcluster to PeTar Unit Handling

- **Always pass `-C 5`** — it emits the headerless `test.dat.10` that `petar.init` expects; the default `test.txt` header is ingested as a zero-mass ghost particle (N+1). Verify N in the converted snapshot header afterwards.
- `mcluster` output unit depends on generator `-u`.
- Common case: `mcluster -u 1` gives mass `Msun`, position `pc`, velocity `km/s`.
- If target runtime is `petar -u 1`, convert velocity with:
  - `petar.init -v 1.022712165045695`
- Prefer normalizing units in `petar.init` unless user explicitly requests a native-unit workflow.
- For any raw IC source, confirm or infer source mass/length/velocity units before `petar.init`.
- Recommended path: normalize the IC into the target PeTar unit system during `petar.init` so that, for `petar -u 1`, the final IC is in `Msun`, `pc`, `pc/Myr`.
- Advanced native-unit paths require a matching `-G` and consistent conversion of all unit-sensitive runtime and post-processing options; do not recommend this path unless the user explicitly asks to preserve native units.

### petar.init Argument Order

`-f` specifies the **PeTar output snapshot**; the positional argument is the **raw input table**. Writing `petar.init -f <raw_table> <pe_tar_input>` overwrites the raw IC with PeTar-formatted data. Always use `petar.init [...] -f <pe_tar_input> <raw_table>` and verify the `Transfer "<raw_table>" to PeTar input data file "<pe_tar_input>"` message.

### Star-Cluster Generation Inputs (Complete Checklist)

If IC must be generated for a star cluster (e.g., via mcluster), every parameter below must be explicitly confirmed with the user. Do not assume defaults for any physics-defining parameter. If the user does not volunteer a value, ask; do not proceed until all fields are resolved.

| # | Parameter | mcluster flag | Why it matters | Default if not asked |
|---|-----------|---------------|----------------|----------------------|
| 1 | Size scale | `-N` or `-M` | Total particle count or total mass | N/A — must be specified |
| 2 | Density profile model + parameters | `-P` | Plummer, King (W0), fractal, etc. | N/A — must be specified |
| 3 | Half-mass radius | `-R` | Controls density and dynamical timescale; directly affects evolution speed and escape rate | **Not safe to default** — ask user |
| 4 | Virial ratio (Q) | `-Q` | Q=0.5 = virial equilibrium; Q=0 = cold collapse; Q>0.5 = expanding. Drives initial dynamical phase | **Not safe to default** — ask user |
| 5 | IMF model + parameters | `-f` | Kroupa, Salpeter, top-heavy, etc. Determines stellar mass distribution | **Not safe to default** — ask user |
| 6 | Stellar mass range (min, max) | `-m` | Lower (typically 0.08 Msun) and upper mass limits | **Not safe to default** — confirm with user |
| 7 | Binary fraction + distribution model | `-b` | Fraction of binaries, period/semi-major axis distribution, mass-ratio distribution | Must be confirmed even for zero binaries |
| 8 | Mass segregation | `-S` | S=0 = no primordial segregation; S>0 = segregated | mcluster default (0) — confirm with user |
| 9 | Random seed | `-s` | Reproducibility. seed=0 = automatic (non-reproducible) | mcluster default (0) — confirm with user |
| 10 | External potential COM position | `-c` (via petar.init) | Cluster center-of-mass initial position (x,y,z) in simulation units | **Not safe to default** — ask user |
| 11 | External potential COM velocity | `-c` (via petar.init) | Cluster center-of-mass initial velocity (vx,vy,vz) in simulation units | **Not safe to default** — ask user |

**Enforcement rule**: This table is the single source of truth for star-cluster IC parameters — `assets/minimal-question-sets.md` points here rather than restating it. Before composing any mcluster command, iterate all applicable rows above; for each row without a user-provided value, ask. Do not proceed to command generation until every applicable field is resolved.

## Performance Optimisation: Tree Time Step Tuning

### Why tune the tree time step?

PeTar auto-estimates a tree time step (`-s`) from the initial particle distribution.
However, this auto estimate is conservative and often suboptimal for production runs.
Tuning `-s` with `petar.find.dt` can give **2–5× speedup** while maintaining accuracy.

### Tool: petar.find.dt

`petar.find.dt` benchmarks the solver across several tree time step candidates and reports the fastest one; it creates `check.perf.<timestep>.log` files.

**Key constraint**: `-a` carries only scenario options — never `-o`, `-w`, `-t`, `-i`, or `-s`, which the tool sets internally. Full flag reference and command template: `assets/script-tools.md`.

**Unit-mode constraint**: `-a` must repeat the production unit mode (`-u <mode>`, plus `-G`/scale factors if used); omitting `-u 1` makes the benchmark read the snapshot with the wrong gravitational constant and abort with a misleading `neighbor address list full`.

### Critical Caveat: Halve the Recommended Step

`petar.find.dt` only benchmarks the **first 6 solver steps**. In many simulations, the optimal tree time step for the initial configuration leads to increasingly large changeover radii as the cluster evolves (e.g., half-mass radius grows from mass loss, or binaries harden), causing the direct-integration (hard) part to slow down significantly after some time.

**Rule**: Always use **`s_rec / 2`** (the next smaller regularized step, i.e. half the recommended value) for the production run. This provides a safety margin against late-time performance degradation. The cost is only ~2× more tree steps while avoiding the risk of the hard solver becoming a bottleneck.

### Workflow Integration

Insert the `petar.find.dt` step **between IC preparation and the production run**:

```
IC generation (mcluster) → petar.init → petar.select → petar.find.dt → petar (full production run)
```

Step by step:

1. Generate IC and convert with `petar.init` (as usual).
2. Run `petar.select` to pick the solver binary matching your scenario.
3. Run `petar.find.dt` directly on the input snapshot — it handles all internal test runs automatically:
   ```
   petar.find.dt -a "-u 1 <scenario-opts>" -i 1 input
   ```
   Use `-i 1` for ASCII-format input (default from `petar.init`), `-i 0` for binary.  
   Pass the unit mode plus all scenario-specific options (e.g. `--galpy-set`, `--bse-metallicity`, `-b`) inside `-a "..."`.
4. **Halve** the recommended `s` value — use the next regularized step (×0.5).
5. Launch the production run with `-s <halved-value>` plus the other options.

## Binary Selection and Capability Validation

### Binary Selection

1. Determine required tokens from scenario:
   - interrupt: `merger`, `bse`, `mobse`, `bseEmp`, `dsm`
   - external: `galpy`, `agama`
   - external-hard: `gasdrag`
   - PN: `pn*`
   - precision: `mp` from `--enable-mpfrc`
2. Select with `petar.select --require ... [--optional ...]`.
3. Keep physics/structure-sensitive tokens in `--require`.
4. Keep performance tokens in `--optional` (for example `avx512,avx2,omp,mpi`).
5. If no candidate matches, stop and provide exact reconfigure/build hints, then explicitly ask the user whether to proceed with `./configure` + `make install` and wait for confirmation before executing any build command.
6. **Before comparing results across builds** (e.g. KDK vs KDKDK4, 32-bit vs 64-bit), confirm both binaries share the same version string (Gate 4). Otherwise the differences reflect code evolution, not the variable under test.

Binary suffix → scenario inference (default mapping, helper-binary classification): `assets/binary-scenario-map.md`.

Physics- or structure-sensitive families should not be treated as ranking-only preferences.
Use `--require` for interruption modules, external-potential families, `gasdrag`, `pn*`, and `mp`/MPFRC style precision requirements.

### Option Validation

Before composing runnable commands:

1. `command -v <binary>`
2. `<binary> -h` and extract supported options.
3. Reject unsupported options and explain why.
4. If the binary is a helper (`.hard.debug`, `.format.transfer`, `.hard.test`, `.dump2test`), switch to the matching solver binary — see Gate 4.

## Changeover Radius and Tree Time Step

### Basic Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-r <r_out>` | 0.0 (auto) | Outer changeover radius |
| `--r-ratio <ratio>` | 0.1 | Inner/outer ratio: `r_in = r_ratio * r_out` |
| `-s <dt_soft>` | 0.0 (auto) | Tree (PT) time step, regularized to 0.5^n |

### `-s` ↔ `-r` Coupling: Two Switching Criteria

Controlled by `--dt-soft-sigma-factor`:

- **Freefall-based** (default, `--dt-soft-sigma-factor=0`): `dt_soft = P(r_in) / nstep`, where `P(r_in)` is the binary period at semi-major axis `r_in`. `nstep` is set by `--dt-soft-kepler-nstep` (default 16, recommended 64 for high accuracy). Mass-dependent via per-particle mass-weighted `r_in`. See Wang et al. 2026, ApJ, 998, 233.
- **σ-based** (`--dt-soft-sigma-factor > 0`, recommended 0.2): `dt_soft = alpha * r_in / (sqrt(3) * sigma_3D)`. Classic criterion from Iwasawa et al. 2015. When `-r=0`, `dt_soft` is first auto-estimated as `2.6e-4 * G*M / sigma_3D^3`, then `r_out` and `r_in` are derived.

**Criterion selection** (Wang et al. 2026):
- Low-σ / loose clusters, wide binaries, subvirial/fractal ICs → **freefall-based** (σ hard to measure, or σ-based gives too-small r_in).
- High-σ / dense clusters (N ≥ 10⁴) → **σ-based** (larger r_in for same dt_soft).
- Both similar at optimal dt_soft for virial equilibrium N≈10³. Transition: `dt_soft / t_ch ∝ 1/N`.

**Environment considerations:**
- **Star cluster** (default): σ well-defined; both criteria applicable; auto chain works.
- **Isolated binary / few-body scattering**: σ undefined → auto defaults will fail. Must set `-r`, `-s` explicitly. Use freefall-based criterion directly: choose `-s` and let the code derive `r_out`, or choose `-r` and set `--dt-soft-kepler-nstep` for accuracy. SDAR handles close encounters by default (`--r-group` auto from `r_in`).
- **Stellar disk / DSM**: rotation-dominated, σ inapplicable. Use freefall-based criterion. Set `-r` explicitly from disk scale height or the encounter region of interest. `--r-group` should be smaller than the disk scale height to avoid false SDAR triggers from shear.
- **Unknown / mixed**: fall back to **freefall-based criterion** (no σ dependency, works for any particle configuration). Set `-r` and `-s` explicitly when possible; avoid relying on auto defaults unless the system is clearly a virialized cluster.

### Related Parameters

- **`--dt-soft-kepler-nstep`** (default 16.0): Steps per orbit for freefall criterion. ns=16 matches pure Hermite final error; ns=64 matches peak error (recommended for precision). Active only with freefall criterion.
- **`--dt-soft-sigma-factor`** (default 0.0): 0 = use freefall criterion; > 0 = use σ-based criterion with this alpha value.
- **`--r-search-min`** (default 0.0, auto): Hard-cluster neighbor search radius. Auto: `max(search_vel_factor * sigma_1D * dt_soft + r_out, 1.2 * r_out)`. Override when auto-detected radius misses a specific binary population.
- **`--r-search-group`** (default -1.0, auto): Group candidate search radius. Auto: `1.0 * r_in`. Set to 0 to disable SDAR. Controls which particles are considered for SDAR group membership.
- **`--r-group`** (default -1.0, auto): Multi-body detection radius and tidal tensor box size. Auto: `0.8 * r_search_group`. Must satisfy `r_group < r_in` for SDAR to activate on binaries.
  - **Accuracy trade-off**: SDAR (LogH) is designed and most accurate for **2-body** (binary) systems. Multi-body SDAR groups (hierarchical triples, 3+ bodies) have significantly degraded accuracy — the LogH method loses its Kepler-solver property when a third body is present (Wang 2025, ApJ, 978, 65). The hybrid BlogH method improves accuracy for weakly perturbed triples but is not yet integrated into PeTar's runtime SDAR.
  - **r_group too large**: risks capturing multiple stars into one SDAR group → multi-body SDAR (lower accuracy, especially for secular evolution like Kozai–Lidov). In disk/DSM environments, excessive r_group may trigger false SDAR groups from shear — keep r_group smaller than the disk scale height.
  - **r_group too small**: misses real binaries → they remain in the P3T Hermite/leapfrog integrator (also less accurate and slower for tight binaries).
  - **Default auto value** (0.8 * r_search_group ≈ 0.8 * r_in) is a reasonable balance for typical star clusters. Override only when a specific binary population is being missed (increase r_group) or too many multi-body groups degrade accuracy (decrease r_group).

**Do not confuse** `--r-search-min` (hard-binary neighbor search) with the changeover boundary (`-r` / `--r-ratio`). `--r-group` and `--r-search-group` control SDAR group detection, not force switching.

**SDAR deep-dive routing**: When a task requires replicating a close-encounter / few-body subsystem in standalone SDAR, debugging SDAR group detection (`--r-group` / `--r-search-group`), or modifying `SDAR/src/`, read the SDAR repository's [sdar-fewbody-integration skill](../../../../SDAR/.github/skills/sdar-fewbody-integration/SKILL.md) and [SDAR AGENTS.md](../../../../SDAR/AGENTS.md) first — they carry strict standalone workflow rules (working directory, input format, command logging, `import sdar` post-processing, and version consistency).

### Auto-Detection Chain (All Defaults)

```
-s = 0, -r = 0, --dt-soft-sigma-factor = 0:
  dt_soft = 2.6e-4 * G*M / sigma^3          (auto-estimate)
  rout derived from dt_soft via freefall criterion
  r_in = r_ratio * rout
  r_search_min = max(vel_factor*sigma*dt_soft + rout, 1.2*rout)
  r_search_group = r_in
  r_group = 0.8 * r_search_group
```

To fix r_in and vary changeover width:
```bash
# r_in = 0.01 pc, r_ratio = 0.5 → r_out = 0.02 pc (narrow)
-r 0.02 --r-ratio 0.5

# r_in = 0.01 pc, r_ratio = 0.1 → r_out = 0.10 pc (wide)
-r 0.10 --r-ratio 0.1
```

Set `-s` and `-r` explicitly to take full manual control. Any parameter left at default continues to auto-derive.

## Workflow Patterns

Input-source branching in more detail (raw table / existing snapshot / restart / repository sample): `assets/input-source-workflows.md`.

### Raw table to new run

1. Prepare/generate raw IC.
2. Convert via `petar.init`.
3. Select solver via `petar.select`.
4. Run solver.
5. Run post-processing tools with mode-consistent flags.

### Existing PeTar snapshot

1. Skip `petar.init`.
2. Select solver.
3. Run solver and optional downstream tools.

### Restart/resume

1. Use `-p <par_file>` before overrides.
2. Final positional argument is restart snapshot.
3. Copy **all** `<prefix>.par.*` companions alongside `<prefix>.par`; a missing one aborts with `Cannot open file <prefix>.par.<feature>`.
4. Match the **read** half of `-i` to the snapshot being read: a binary snapshot needs `-i 0` (or `-i 3`) — the default `-i 2` reads ASCII and aborts with `cannot read header`. Details and worked example: `assets/input-source-workflows.md`.
5. Use `-a 0` only when overwrite requested.
6. **Duplicate-event prevention is now handled automatically** via the transactional `*.tmp` file mechanism: each output interval is first written to temporary files and only committed on successful completion. This eliminates the need for manual cleanup before restarts in modern PeTar.
7. For **legacy (old-version) runs only**: if restarting from an intermediate snapshot with older output files that may contain already-committed duplicate records, use `petar.data.clear -t <time> [prefix]` to trim event files back to the snapshot time.

## Controlled Experiment Checklist

Before running a parameter scan or cross-build comparison:

### General checks

- [ ] Compared builds share the **same PeTar version** (check `strings ... | grep Version`)
- [ ] Parameters under study are explicitly set (not relying on defaults); verify derived values match intent
- [ ] A control run with known-good parameters (e.g., defaults from sample scripts) behaves as expected before scanning
- [ ] Unit system is consistent across IC generation, `petar.init`, and runtime options
- [ ] Binary family selected via `petar.select` matches the required physics features

### Scenario-specific checks

Refer to `Required Input Checklist` above for scenario-specific mandatory parameters, and to `Environment considerations` under `Changeover Radius and Tree Time Step` for changeover/timestep strategy per system type.

## Post-Processing Tool Policy

Use installed tools directly when user intent matches:

- `petar.data.process`
- `petar.movie`
- `petar.format.transfer.post`
- `petar.get.object.snap`
- `petar.galev.process`
- `petar.data.gether` — legacy; for old MPI runs without automatic `snap.lst`
- `petar.external.pot.movie`
- `petar.galpy.pot.movie`

### Other tools

- `petar.update.par` — migrate legacy `.par` files from old PeTar versions
- `petar.find.dt` — search for suitable tree time step (see `Timestep Tuning` section)
- `petar.galpy.help` — query Galpy potential families and argument/config help
- `petar.external.galpy` / `petar.external.agama` — generate external potential map snapshots from run parameter files for visualization workflows

Prefer the corresponding installed tool command pattern over generic workflow prose for conversion, extraction, movie generation, restart cleanup, and external-potential maps. Per-scenario `petar.data.process` and `petar.movie` defaults: `assets/default-postprocessing.md`.

## petar.movie Usage

`petar.movie` generates movies from `data.snap.lst` (auto-generated at runtime). Mode-consistency flags must match the producing solver family, or snapshot read errors will occur (see `Snapshot Read-Mismatch Policy`).

### Required mode flags

| Flag | Value | When |
|------|-------|------|
| `-i` | `none`, `merger`, `bse`, `mobse` | Match solver interrupt mode |
| `-t` | `none`, `galpy`, `agama` | Match solver external mode |
| `-s` | `ascii`, `binary`, `npy` | `npy` for `petar.data.process` output; `binary` for raw solver output |
| `--snapshot-type` | `origin`, `post` | `post` for processed snapshots; `origin` for raw solver output |
| `-G` | float | Gravitational constant; 0.00449830997959438 for Msun/pc/Myr; 1.0 for Henon |

### Quick examples

```bash
# Particle distribution
petar.movie -m x-y -R 10 data.snap.lst

# HR diagram
petar.movie -H data.snap.lst

# Combined panels for BSE simulation
petar.movie -m x-y -R 10 -H -b -i bse -s npy data.snap.lst
```

Per-scenario argument templates (isolated, BSE, Galpy/Agama, DSM, and combined panels): `assets/default-postprocessing.md`.

### Snapshot format matching

- Raw solver output (`data.*`): use `-s binary --snapshot-type origin`
- Post-processed output (`data.*.single`, `data.*.binary`): use `-s npy --snapshot-type post`

Mismatched flags cause snapshot read errors (see `Snapshot Read-Mismatch Policy`). For full flag reference and comparison-mode usage, see `README.md`.

## Python Data Analysis Tools

### Hard Rule: Read the Authority First

Before writing any Python analysis code that reads PeTar output files, **you MUST read `assets/data-readback-patterns.md`** with `read_file`. Do not guess the API from memory — the class constructors, keyword arguments, header offsets, and `fromfile`/`loadtxt` method signatures vary by file type and solver configuration. Guessing produces read errors and incorrect results.

The primary reference is `sample/data_analysis.ipynb`; the verified readback patterns are in `assets/data-readback-patterns.md`.

### Unit Conversion Rule

**Never use approximate hardcoded constants** (e.g., `semi * 206265`). **Always prefer `astropy.units`** for all unit conversions — it uses IAU 2015/2019 definitions and is authoritative. Use `petar` module constants (e.g., `petar.G_MSUN_PC_MYR`) only when exact consistency with PeTar's C++ solver is required. See `assets/data-readback-patterns.md` (Unit Conversion section) for complete guidance and examples.

### File-Type → Reader-Class Quick Reference

One-line hints only — exact constructor signatures, header offsets, and worked examples are in `assets/data-readback-patterns.md` and must be read from there rather than guessed.

| File pattern | Reader class | Key hint |
|---|---|---|
| `data.<N>` (raw snapshot) | `petar.Particle` | Binary: `offset=` depends on build (see asset); ASCII: `skiprows=1`; header via `petar.PeTarDataHeader` |
| `data.lagr` | `petar.LagrangianMultiple` | `external_mode` controls COM-offset columns |
| `data.core` / `data.status` | `petar.Core` / `petar.Status` | Plain `.fromfile()` |
| `data.esc_single` / `data.esc_binary` | `petar.SingleEscaper` / `petar.BinaryEscaper` | Needs `interrupt_mode`, `external_mode` |
| `data.sse` / `data.bse` | `petar.SSEType` / `petar.BSEType` | ASCII — `.loadtxt()` |
| `data.sse.*` / `data.bse.*` event tables | `petar.SSETypeChange`, `petar.SSESNKick`, `petar.BSETypeChange`, `petar.BSEKick`, `petar.BSEDynamicMerge` | `.loadtxt()`; `gw_kick` → `BSEKick`, `binary_merge` → `BSETypeChange` |
| `data.bse_status` | `petar.BSEStatus` | Binary — `.fromfile()` |
| `data.<N>.single[.npy]` | `petar.Particle` | No header/offset; `.npy` → `.load()` |
| `data.<N>.binary[.npy]` | `petar.Binary` | No header/offset; needs `G=`, `member_particle_type`; `.npy` → `.load()` |
| `data.<N>.triple` / `.quadruple` / `data.group.n<N>` | `petar.GroupInfo(N=…)` | `N=` must match the file |
| `object.<N>` | `petar.Particle` | Call `addNewMember("time", …)` before `fromfile()` |
| `data.prof.rank.*` | `petar.Profile` | Column count varies by GPU/FDPS build |
| `data.*_h4_<N>_<rank>.log` | `petar.HermiteData` | PeTar reader only (not `sdar`); split by column count first; rerun with `OMP_NUM_THREADS=1` |

### Universal Keyword Arguments

These kwargs control column layout and **must match the producing solver**. Pass them consistently to `Particle`, `Binary`, escaper, and Lagrangian readers:

| Argument | Values | Purpose |
|----------|--------|---------|
| `interrupt_mode` | `none`, `merger`, `bse`, `mobse`, `bseEmp`, `dsm` | Must match solver `--with-interrupt` |
| `external_mode` | `none`, `galpy`, `agama` | Must match solver `--with-external` |

Mismatched kwargs cause column misalignment and read errors (see `Snapshot Read-Mismatch Policy`).

### Analysis Modules (import petar.<module>)

| Module | Purpose |
|--------|---------|
| `petar.profile` | Radial density/velocity profiles |
| `petar.lagrangian` | Lagrangian radii |
| `petar.escaper` | Escaper identification |
| `petar.bse` | BSE stellar evolution analysis |
| `petar.bse.pulsar` | Pulsar analysis |
| `petar.dsm` | Disk star merger analysis |
| `petar.external` | External potential analysis |
| `petar.tide` | Tidal analysis |
| `petar.galev` | Galev stellar population synthesis |
| `petar.agama` | Agama MW potential utilities |
| `petar.parallel_data_process` | Parallel multi-snapshot processing |

## Snapshot Read-Mismatch Policy

Treat the following as read failures (not harmless warnings):

- binary size/dtype alignment mismatch
- ASCII shape/column mismatch
- text decode errors caused by format mismatch

**Typical failure signature**: forgetting the header offset (binary) or `skiprows=1` (ASCII) reads the header as the first particle, yielding negative masses and nonsensical IDs. Post-processed snapshots (`data.<N>.single`, `data.<N>.binary`, …) have **no header** and must not use either offset. Correct per-file-type examples: `assets/data-readback-patterns.md` (Pattern 6 for raw snapshots).

Recovery sequence:

1. Stop current downstream command.
2. Re-check producing solver family and runtime output mode.
3. Re-check reader flags (`-i`, `-t`, format flags, snapshot type flags) and Python kwargs (`interrupt_mode`, `external_mode`, header offset).
4. Retry with corrected mode.
5. **Do not write a custom binary/ASCII parser as a workaround.** PeTar output formats are intentionally covered by the installed Python readers; a column mismatch almost always means the reader mode/flags are wrong. If the mode/flags are verified and the mismatch persists, treat it as a potential bug in the reader or the producing solver, stop the analysis, and report the issue to the user with the exact file pattern, reader class, kwargs, and error message.

## Preferred Sources in This Repository

Use these as primary references:

- `README.md`
- `test/functional/README.md`
- `test/validation/README.md`
- `sample/star_cluster_plummer_N1k.sh`
- `sample/star_cluster_plummer_N1k_binaries.sh`
- `sample/star_cluster_plummer_N1k_binaries_bse.sh`
- `sample/star_cluster_plummer_N1k_GalpyMWPot.sh`
- `sample/star_cluster_plummer_N1k_binaries_bse_GalpyMWPot.sh`
- `sample/star_cluster_plummer_N1k_AgamaMWPotHunter24.sh`
- `sample/data_analysis.ipynb`
- `.github/skills/petar-nbody-simulation/assets/option-matrix.md`
- `.github/skills/petar-nbody-simulation/assets/option-reference.md`
- `.github/skills/petar-nbody-simulation/assets/script-tools.md`

Technical background (key algorithms):

- PeTar code description: P3T hybrid method, Hamiltonian splitting, performance benchmarks (Wang et al. 2020, MNRAS, 497, 536) — [https://doi.org/10.1093/mnras/staa1915](https://doi.org/10.1093/mnras/staa1915)
- SDAR integrator: slow-down + time-transformed symplectic method for few-body systems (Wang et al. 2020, MNRAS, 493, 3398) — [https://doi.org/10.1093/mnras/staa480](https://doi.org/10.1093/mnras/staa480)
- BlogH hybrid method and LogH accuracy limits for hierarchical triples (Wang 2025, ApJ, 978, 65) — [https://doi.org/10.3847/1538-4357/ad98f3](https://doi.org/10.3847/1538-4357/ad98f3)
- Freefall-based P3T switching criterion vs σ-based criterion (Wang et al. 2026, ApJ, 998, 233) — [https://doi.org/10.3847/1538-4357/ae367c](https://doi.org/10.3847/1538-4357/ae367c)

## Reference Documents (Must Read)

The following asset files contain critical information not inlined in this document.
When a task falls into the corresponding category, **read the file explicitly** with `read_file` before proceeding — do not guess or rely on memory.

| When to read | File | What it contains |
|-------------|------|------------------|
| **Any simulation scenario** (after scenario is identified) | `assets/minimal-question-sets.md` | Per-scenario required-ask lists, including the "do not ask if already known" rule |
| **Binary selection** (choosing or explaining the solver family) | `assets/binary-scenario-map.md` | Binary-suffix → scenario inference, solver vs helper classification, nearest-match guidance |
| **Tool usage** (before composing commands) | `assets/script-tools.md` | Command templates and usage patterns for all PeTar tools (`petar.init`, `petar.find.dt`, `petar.data.process`, `petar.movie`, etc.) |
| **Post-processing / movie defaults** (full-workflow requests) | `assets/default-postprocessing.md` | Per-scenario `petar.data.process` and `petar.movie` argument templates |
| **Input-source branching** (raw table vs snapshot vs restart) | `assets/input-source-workflows.md` | Detailed workflow per input source, including the canonical restart command form |
| **Python data analysis** (before writing any analysis code) | `assets/data-readback-patterns.md` | **MUST READ before writing any analysis code.** Verified readback patterns for all file types, constructor kwargs, header offsets (including MPFRC variants), `fromfile` vs `loadtxt`/`load` choice, mode-matching rules, and blocking-warning classification |
| **DSM workflow** (when scenario = DSM) | `assets/dsm-workflow.md` | DSM-specific IC preparation, runtime parameters, and post-processing pipeline |
| **Skill maintenance / new-machine recovery** (never for simulation tasks) | `assets/prompt-starters.md`, `HANDOFF.md` | Recovery and regression prompts, and the current skill-maintenance state |

## Scope Notes

- **This skill is the authority for running PeTar simulations and analysing PeTar output, and it is self-sufficient** — completing a simulation task never requires reading `AGENTS.md` or an agent definition. Agent files cover PeTar *development* and reference this skill instead of restating it.
- Execution-correctness content (gates, mode consistency, required inputs, error-prone constraints) is kept in this file. Reference-style content (full flag lists, class APIs, per-scenario templates) lives in `assets/` and is reached through "Reference Documents (Must Read)" above.
- **SDAR internals are owned by the SDAR repository** (`SDAR/src/`, standalone AR/Hermite builds, the `sdar` Python readers) — see the [sdar-fewbody-integration skill](../../../../SDAR/.github/skills/sdar-fewbody-integration/SKILL.md) and [SDAR AGENTS.md](../../../../SDAR/AGENTS.md). PeTar-specific difference: PeTar exposes SDAR group detection as runtime options (`--r-group`, `--r-search-group`, auto-derived from `r_in`), whereas standalone SDAR takes `-r-group` directly on the command line.
- Functional smoke defaults are documented in [README.md](../../../README.md) and [test/functional/README.md](../../../test/functional/README.md).
- If a required detail is not in this file, query source docs rather than guessing.
