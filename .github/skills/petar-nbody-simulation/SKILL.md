---
name: petar-nbody-simulation
description: "Use when: setting up or running PeTar N-body simulations, including petar.select binary-family switching, petar.init conversion, isolated clusters, BSE/SSE, Galpy or Agama external potential, MPI/OpenMP/GPU launch, restart/resume, petar.data post-processing, timestep tuning with petar.find.dt, output gathering with petar.data.gether, restart cleanup with petar.data.clear, snapshot extraction with petar.get.object.snap, format conversion with petar.format.transfer.post, galev processing, and movie generation."
---

# PeTar N-body Simulation Skill (Compact)

## Purpose

Provide strict, command-level guidance for PeTar workflows in this repository.
This compact version prioritizes execution safety, option correctness, and reproducibility.

## Non-Negotiable Rules

- Infer required physics/runtime features first, then switch binary family with `petar.select`.
- Do not run whatever `petar` currently points to unless selected for this scenario.
- Ask only for missing required inputs; do not guess mandatory physics parameters.
- Before execution, provide parameter summary + exact commands; run only after user confirmation.
- Default output prefix is `data` when user does not provide one.
- Default unit mode is `-u 1` when user does not provide one.
- Enforce unit consistency across IC generation, `petar.init`, runtime options, and post-processing.
- For stellar evolution requests, require module confirmation: `merger`, `bse`, `bseEmp`, `mobse` (`moBSE`), or `dsm`.
- For external potential requests, require mode confirmation: `galpy` or `agama`.
- For BSE-family runs, metallicity is mandatory.
- For galactic/external potential runs, potential model and cluster COM phase-space are mandatory.
- Exclude debug families unless user explicitly asks for debug.
- Exclude GPU families unless user explicitly asks for GPU.
- Validate all custom options against selected binary `-h` before emitting runnable commands.
- Never use `*.hard.debug` or `*.format.transfer` binaries as production solvers.
- For source-level debugging, require rebuild with `--with-debug=g`.
- Prefer `petar.select` over manual symlink edits.
- If MPI output must feed downstream post-processing/movie/extraction, run `petar.data.gether` first.
- Do not use `petar.init` in restart/resume workflows.
- For post-processing tools, ensure mode flags match producing solver family and snapshot format.
- For functional smoke automation, required test inputs must be git-tracked; output artifacts are not inputs.

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

### Star-Cluster Generation Inputs

If IC must be generated for a star cluster, require:

1. Size scale (total mass or total star count).
2. Half-mass radius.
3. Density profile model (+ profile parameters if needed).
4. IMF model.
5. Binary fraction + distribution model.

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
- `petar.data.gether`
- `petar.external.pot.movie`

Do not give only abstract advice when a matching tool exists.

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
- `.github/skills/petar-nbody-simulation/assets/option-matrix.md`
- `.github/skills/petar-nbody-simulation/assets/script-tools.md`

## Scope Notes

- This skill is intentionally compact.
- Use repository docs and sample scripts for long-form explanations and example-heavy guidance.
- Functional smoke defaults are documented in [README.md](README.md) and [test/functional/README.md](test/functional/README.md).
- If a required detail is not in this file, query source docs rather than guessing.
