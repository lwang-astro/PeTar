---
name: petar-nbody-simulation
description: "Use when: setting up or running PeTar N-body simulations, including petar.select binary-family switching, petar.init conversion, isolated clusters, BSE/SSE, Galpy or Agama external potential, MPI/OpenMP/GPU launch, restart/resume, petar.data post-processing, timestep tuning with petar.find.dt, restart cleanup with petar.data.clear, snapshot extraction with petar.get.object.snap, format conversion with petar.format.transfer.post, galev processing, and movie generation."
---

# PeTar N-body Simulation Skill (Compact)

## Purpose

Provide strict, command-level guidance for PeTar workflows in this repository.
This compact version prioritizes execution safety, option correctness, and reproducibility.

## Non-Negotiable Rules

- Infer required physics/runtime features first, then switch binary family with `petar.select`.
- Do not run whatever `petar` currently points to unless selected for this scenario.
- Ask only for missing required inputs; do not guess mandatory physics parameters.
- Do not proceed with a run until all scenario-required inputs are present.
- Before execution, provide parameter summary + exact commands; run only after user confirmation.
- Default output prefix is `data` when user does not provide one.
- Default unit mode is `-u 1` when user does not provide one.
- Enforce unit consistency across IC generation, `petar.init`, runtime options, and post-processing.
- For stellar evolution requests, require module confirmation: `merger`, `bse`, `bseEmp`, `mobse` (`moBSE`), or `dsm`.
- For external potential requests, require mode confirmation: `galpy` or `agama`.
- For star-cluster IC generation, do not proceed until the IC parameter set is complete.
- For BSE-family runs, metallicity is mandatory.
- For galactic/external potential runs, potential model and cluster COM phase-space are mandatory.
- Exclude debug families unless user explicitly asks for debug.
- Exclude GPU families unless user explicitly asks for GPU.
- Validate all custom options against selected binary `-h` before emitting runnable commands.
- Never use `*.hard.debug` or `*.format.transfer` binaries as production solvers.
- For source-level debugging, require rebuild with `--with-debug=g`.
- Prefer `petar.select` over manual symlink edits.
- For current PeTar runs (tmp file mechanism), `data.snap.lst` is generated at runtime and stream outputs are automatically committed — no gathering step is needed. `petar.data.gether` is a legacy tool only required for old MPI runs that lack `snap.lst`.
- Do not use `petar.init` in restart/resume workflows.
- Do not answer analysis, conversion, movie generation, or restart-cleanup requests only with abstract advice when an installed PeTar tool exists for that task.
- For post-processing tools, ensure mode flags match producing solver family and snapshot format.
- For functional smoke automation, required test inputs must be git-tracked; output artifacts are not inputs.
- When modifying source headers that define CLI options (`src/petar.hpp`, `src/hard.hpp`, `bse-interface/*.h`, `galpy-interface/*.h`, `agama-interface/*.h`, `src/disk_star_merger.hpp`, `src/gas_drag.hpp`, `src/external_hard.hpp`, `parallel-random/rand_io.hpp`), remind the user to re-run `.github/skills/petar-nbody-simulation/assets/generate_option_reference.py` so that `option-reference.md` stays synchronized.

### Environment requirements

- Before any PeTar launch, require `export OMP_STACKSIZE=128M` (and `ulimit -s unlimited` on Linux). Missing this causes segfault on large simulations.
- For MPI+OpenMP launches, recommend `--bind-to none` to prevent processes being restricted to a single core.
- For MPI launches in UCX environments, require `export UCX_VFS_ENABLE=n` to avoid a known UCX segfault.
- When `./configure` runs on a cluster login node, warn that auto-detected SIMD capabilities may differ from compute nodes; suggest verifying with `--enable-avx2` / `--enable-avx512` flags.

### bseEmp-specific

- `bseEmp` requires manually linking `ffbonn` or `ffgeneva` metal-poor track directories before first use. Without this, the binary crashes at initialization with a file-not-found error. Ask the user to set up these links before composing a bseEmp run command.

### Coordinate origin

- Warn when the density center is far from `(0,0,0)` or includes very distant particles — this can trigger `n_jp<=pg.NJMAX` assertion failures due to tree cell resolution limits.

## Required Input Checklist

Collect the minimum required fields before composing run commands.

### Common

1. Working directory.
2. Scenario: isolated | BSE/SSE | DSM | Galpy | Agama | restart | standalone-bse.
3. End time `-t`.
4. Output interval `-o`.
5. Launch mode: serial | OpenMP | MPI+OpenMP | GPU.
6. Initial data source: existing PeTar snapshot | raw table | generator.
7. Extra user options to pass through (validated only).

### Scenario-Specific

- BSE/mobse/bseEmp:
  - primordial binary count `-b`
  - metallicity (`--bse-metallicity` or family-equivalent)
- DSM:
  - only requested DSM controls (`--dsm-*`)
- Galpy:
  - one of `--galpy-set` / `--galpy-conf-file` / `--galpy-type-arg`
  - optional `--galpy-rscale`, `--galpy-vscale`
  - note: PeTar only supports Galpy ≤ 1.10.2
- Agama:
  - `--agama-conf-file`
  - optional `--agama-rscale`, `--agama-vscale`
- Restart:
  - restart snapshot
  - parameter file (`-p`, default `data.par` if appropriate)
  - append/overwrite behavior (`-a 1`/`-a 0`)
- standalone-bse:
  - use `petar.bse` (or `petar.mobse`, `petar.bseEmp`) binary directly, not the N-body solver
  - no IC generation or `petar.init` needed
  - see `README.md` for usage

## IC and Unit Rules

### Generator Availability

When generator-based star-cluster IC creation is required:

1. Prefer `mcluster_gpu` if available.
2. Else use `mcluster` if available.
3. If neither exists, stop and ask user whether to proceed with installation guidance.

### mcluster to PeTar Unit Handling

- `mcluster` output unit depends on generator `-u`.
- Common case: `mcluster -u 1` gives mass `Msun`, position `pc`, velocity `km/s`.
- If target runtime is `petar -u 1`, convert velocity with:
  - `petar.init -v 1.022712165045695`
- Prefer normalizing units in `petar.init` unless user explicitly requests a native-unit workflow.
- For any raw IC source, confirm or infer source mass/length/velocity units before `petar.init`.
- Recommended path: normalize the IC into the target PeTar unit system during `petar.init` so that, for `petar -u 1`, the final IC is in `Msun`, `pc`, `pc/Myr`.
- Advanced native-unit paths require a matching `-G` and consistent conversion of all unit-sensitive runtime and post-processing options; do not recommend this path unless the user explicitly asks to preserve native units.

### Star-Cluster Generation Inputs

If IC must be generated for a star cluster, require:

1. Size scale (total mass or total star count).
2. Half-mass radius.
3. Density profile model (+ profile parameters if needed).
4. IMF model.
5. Binary fraction + distribution model.

Do not continue to command generation until this parameter set is complete.

## Timestep Tuning with petar.find.dt

`petar.find.dt` searches for a suitable tree time step for a given snapshot and launch configuration:

1. Run the solver for a short test duration with data output every step.
2. Run `petar.find.dt` on the output snapshots to get a recommended `-s` value.
3. Restart the simulation with the recommended time step.

This is especially useful for simulations where the automatic time step estimate is suboptimal.

## Binary Selection and Capability Validation

### Binary Selection

1. Determine required tokens from scenario:
   - interrupt: `base`, `bse`, `mobse`, `bseEmp`, `dsm`
   - external: `galpy`, `agama`
   - external-hard: `gasdrag`
   - PN: `pn*`
   - precision: `mp` from `--enable-mpfrc`
2. Select with `petar.select --require ... [--optional ...]`.
3. Keep physics/structure-sensitive tokens in `--require`.
4. Keep performance tokens in `--optional` (for example `avx512,avx2,omp,mpi`).
5. If no candidate matches, stop and provide exact reconfigure/build hints.

Physics- or structure-sensitive families should not be treated as ranking-only preferences.
Use `--require` for interruption modules, external-potential families, `gasdrag`, `pn*`, and `mp`/MPFRC style precision requirements.

### Option Validation

Before composing runnable commands:

1. `command -v <binary>`
2. `<binary> -h` and extract supported options.
3. Reject unsupported options and explain why.
4. If binary name ends with `.hard.debug` or `.format.transfer`, switch to matching solver binary.

## Workflow Patterns

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
3. Use `-a 0` only when overwrite requested.
4. If restarting from an intermediate snapshot that may duplicate events, use `petar.data.clear` first.

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

Do not give only abstract advice when a matching tool exists.
If the request is specifically about conversion, extraction, movie generation, restart cleanup, or external-potential maps, prefer the corresponding installed tool command pattern over generic workflow prose.

## petar.movie Usage

`petar.movie` generates movies from `data.snap.lst` (auto-generated at runtime). Mode-consistency flags must match the producing solver family, or snapshot read errors will occur (see `Snapshot Read-Mismatch Policy`).

### Required mode flags

| Flag | Value | When |
|------|-------|------|
| `-i` | `none`, `merger`, `base`, `bse`, `mobse` | Match solver interrupt mode |
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

### Snapshot format matching

- Raw solver output (`data.*`): use `-s binary --snapshot-type origin`
- Post-processed output (`data.*.single`, `data.*.binary`): use `-s npy --snapshot-type post`

Mismatched flags cause snapshot read errors (see `Snapshot Read-Mismatch Policy`). For full flag reference and comparison-mode usage, see `README.md`.

## Python Data Analysis Tools

PeTar installs a Python analysis library (`import petar`) at `tools/analysis/`. The primary reference is `sample/data_analysis.ipynb`.

### Quick start

```python
import petar

header = petar.PeTarDataHeader("data.1")        # read snapshot header
particle = petar.Particle()                      # read particles
particle.fromfile("data.1", offset=petar.HEADER_OFFSET)
```

### Mode-specific reading (required for correctness)

The `Particle` columns depend on the solver's compile-time configuration. Pass matching keyword arguments:

```python
# BSE simulation
particle = petar.Particle(interrupt_mode="bse")
particle.fromfile("data.1", offset=petar.HEADER_OFFSET)

# Galpy external potential
particle = petar.Particle(external_mode="galpy")
particle.fromfile("data.1", offset=petar.HEADER_OFFSET)
```

**Key keyword arguments for `Particle()`:**

| Argument | Values | Purpose |
|----------|--------|---------|
| `interrupt_mode` | `none`, `merger`, `base`, `bse`, `mobse`, `dsm` | Must match solver `--with-interrupt` |
| `external_mode` | `none`, `galpy`, `agama` | Must match solver `--with-external` |

Mismatched keyword arguments cause column misalignment and read errors (see `Snapshot Read-Mismatch Policy`).

### Analysis modules

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

For detailed class API (header offsets, keyword arguments, member lists), see `sample/data_analysis.ipynb`.

## Snapshot Read-Mismatch Policy

Treat the following as read failures (not harmless warnings):

- binary size/dtype alignment mismatch
- ASCII shape/column mismatch
- text decode errors caused by format mismatch

Recovery sequence:

1. Stop current downstream command.
2. Re-check producing solver family and runtime output mode.
3. Re-check reader flags (`-i`, `-t`, format flags, snapshot type flags).
4. Retry with corrected mode.
5. If mismatch persists, stop and report outputs as untrustworthy.

## Minimal Interaction Rule

When invoked:

1. Infer scenario and binary requirements from user request.
2. Ask only missing required inputs.
3. Present concise parameter summary + exact command block.
4. Ask confirmation.
5. Execute after confirmation.

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

## Scope Notes

- This skill prioritizes execution correctness: it keeps mode-consistency rules, required inputs, and error-prone workflow constraints in full detail.
- Reference-style content (full flag lists, class API docs, comparison-mode usage) is kept in repository docs (`README.md`, `sample/data_analysis.ipynb`) to avoid duplication.
- Functional smoke defaults are documented in [README.md](README.md) and [test/functional/README.md](test/functional/README.md).
- If a required detail is not in this file, query source docs rather than guessing.
