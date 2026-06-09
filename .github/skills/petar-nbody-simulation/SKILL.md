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

## Required Input Checklist

Collect the minimum required fields before composing run commands.

### Common

1. Working directory.
2. Scenario: isolated | BSE/SSE | DSM | Galpy | Agama | restart.
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
- Agama:
  - `--agama-conf-file`
  - optional `--agama-rscale`, `--agama-vscale`
- Restart:
  - restart snapshot
  - parameter file (`-p`, default `data.par` if appropriate)
  - append/overwrite behavior (`-a 1`/`-a 0`)

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

Do not give only abstract advice when a matching tool exists.
If the request is specifically about conversion, extraction, movie generation, restart cleanup, or external-potential maps, prefer the corresponding installed tool command pattern over generic workflow prose.

## petar.movie Usage

`petar.movie` generates movies from simulation snapshots. It reads a snapshot list file that is automatically generated at runtime as `data.snap.lst` (with the default prefix `data`). No separate gathering step is needed for current PeTar runs.

### Input file

`data.snap.lst` is created during simulation runtime and can be used directly by `petar.movie` and `petar.data.process`. For legacy runs without this file, use `ls | egrep '^data.[0-9]+$' | sort -n -k 1.6 > snap.lst` or `petar.data.gether` as a fallback.

### Common Scenarios

#### 1. Particle position distribution in x-y plane

```bash
petar.movie -m x-y -R 10 -o particle_movie.mp4 data.snap.lst
```

#### 2. HR diagram alone

```bash
petar.movie -H -o hr_movie.mp4 data.snap.lst
```

#### 3. Combined: particle + HR + semi-ecc diagram for BSE simulations

```bash
petar.movie -m x-y -R 10 -H -b -i bse -s npy -o combined.mp4 data.snap.lst
```

#### 4. Lagrangian radii evolution

Requires `data.lagr` from `petar.data.process --calc-lagrangian`:

```bash
petar.movie -L data.lagr -o lagr_movie.mp4 data.snap.lst
```

### Mode-consistency flags

These must match the producing solver family:

| Flag | Value | When |
|------|-------|------|
| `-i` | `none`, `merger`, `base`, `bse`, `mobse` | Match solver interrupt mode |
| `-t` | `none`, `galpy`, `agama` | Match solver external mode |
| `-s` | `ascii`, `binary`, `npy` | `npy` when snapshots are from `petar.data.process`; `binary` when directly from solver output |
| `--snapshot-type` | `origin`, `post` | `post` for `petar.data.process` output (single/binary files); `origin` for raw solver output |
| `-G` | float | Gravitational constant; default 0.00449830997959438 (Msun, pc, Myr); use 1.0 for Henon units |

### Comparison mode

Compare multiple simulations side-by-side by preparing a model list file and using `-l`:

```bash
petar.movie -m x-y -l models.lst -o comparison.mp4 snapshots.lst
```

### Overlay external potential map

Combine with `petar.external.pot.movie` output or overlay directly:

```bash
petar.movie -m x-y -R 20 --ext-pot -o with_pot.mp4 data.snap.lst
```

### Snapshot format matching for movie

- Raw solver output (`data.*`): use `-s binary --snapshot-type origin`
- Post-processed output (`data.*.single`, `data.*.binary`): use `-s npy --snapshot-type post`

Mismatched flags cause snapshot read errors; see `Snapshot Read-Mismatch Policy`.

## Python Data Analysis Tools

PeTar provides a Python analysis library at `tools/analysis/`, importable as `import petar` after installation. This is the primary interface for custom data analysis beyond what `petar.data.process` and `petar.movie` provide.

### Import Pattern

```python
import petar
import numpy as np
import matplotlib.pyplot as plt
```

Reference: `sample/data_analysis.ipynb`

### Core Classes

#### PeTarDataHeader — read snapshot header

```python
header = petar.PeTarDataHeader("data.1")
print(header.time, header.n, header.file_id)
```

For snapshots with external potential offsets (Galpy/Agama), pass `external_mode`:

```python
header = petar.PeTarDataHeader("data.1", external_mode="galpy")
# header.pos_offset, header.vel_offset are now populated
```

#### Particle — read snapshot particle data

The `Particle` class inherits from `HardParticle` → `BaseParticle` → `SimpleParticle`. The exact columns depend on compile-time configuration, controlled by keyword arguments:

```python
# Pure gravity
particle = petar.Particle()
particle.fromfile("data.1", offset=petar.HEADER_OFFSET)

# BSE simulation
particle = petar.Particle(interrupt_mode="bse")
particle.fromfile("data.1", offset=petar.HEADER_OFFSET)

# Galpy external potential
particle = petar.Particle(external_mode="galpy")
particle.fromfile("data.1", offset=petar.HEADER_OFFSET)
```

**Key keyword arguments for `Particle.__init__`:**

| Argument | Values | Purpose |
|----------|--------|---------|
| `interrupt_mode` | `none`, `merger`, `base`, `bse`, `mobse`, `dsm` | Controls stellar evolution columns |
| `external_mode` | `none`, `galpy`, `agama` | Controls external potential column |
| `use_mpfrc` | bool | High-precision position parts |
| `collect_sp_acc` | bool | Superparticle acceleration column |

**Common particle members** (vary by mode):
- `mass`, `pos`, `vel` — inherited from SimpleParticle
- `binary_state` — interruption state flag
- `radius`, `dm`, `star` — when interrupt_mode is set
- `r_search`, `id` — from HardParticle
- `acc_soft`, `pot`, `pot_soft` — from Particle
- `pot_ext` — when external_mode is not `none`

#### Binary — read binary/triple systems

After `petar.data.process` generates `data.*.binary` files:

```python
binary = petar.Binary()
binary.fromfile("data.1.binary")

# Access members: binary.semi, binary.ecc, binary.mass, binary.pos, binary.vel
# binary.p1, binary.p2 are the two component particles
```

#### PeTarData — load full processed snapshot

```python
data = petar.PeTarData("data.1", interrupt_mode="bse")
# data.single  — single star Particle data
# data.binary  — Binary data
# data.header  — PeTarDataHeader
```

### Snapshot offset constants

Use the correct header offset for the particle data:

| Constant | Value | When |
|----------|-------|------|
| `petar.HEADER_OFFSET` | 24 | Default binary format |
| `petar.HEADER_OFFSET_F128` | 32 | Float128 build |
| `petar.HEADER_OFFSET_WITH_CM` | 72 | With center-of-mass in header |
| `petar.HEADER_OFFSET_WITH_CM_F128` | 128 | Float128 + CM in header |

### Analysis Modules

After loading data, use sub-modules for specific analyses:

| Module | Purpose |
|--------|---------|
| `petar.profile` | Radial density/velocity profiles |
| `petar.lagrangian` | Lagrangian radii computation |
| `petar.escaper` | Escaper identification |
| `petar.bse` | BSE stellar evolution analysis |
| `petar.bse.pulsar` | Pulsar-specific analysis |
| `petar.dsm` | Disk star merger analysis |
| `petar.external` | External potential analysis |
| `petar.tide` | Tidal interaction analysis |
| `petar.galev` | Galev stellar population synthesis |
| `petar.agama` | Agama MW potential utilities |
| `petar.parallel_data_process` | Parallel processing of multiple snapshots |

### Reading ASCII snapshots

If snapshots were converted to ASCII via `petar.format.transfer`:

```python
particle = petar.Particle()
particle.fromfile("data.1.A", snapshot_format="ascii")
```

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

- This skill is intentionally compact.
- This compact version still keeps the mandatory execution guardrails from earlier long-form versions in `Non-Negotiable Rules` and the most error-prone unit/workflow sections.
- Use repository docs and sample scripts for long-form explanations and example-heavy guidance.
- Functional smoke defaults are documented in [README.md](README.md) and [test/functional/README.md](test/functional/README.md).
- If a required detail is not in this file, query source docs rather than guessing.
