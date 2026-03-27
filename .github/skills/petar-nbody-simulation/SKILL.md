---
name: petar-nbody-simulation
description: "Use when: setting up or running PeTar N-body simulations, including petar.init conversion, isolated clusters, BSE/SSE, Galpy or Agama external potential, MPI/OpenMP/GPU launch, restart/resume, petar.data post-processing, timestep tuning with petar.find.dt, output gathering with petar.data.gether, restart cleanup with petar.data.clear, snapshot extraction with petar.get.object.snap, format conversion with petar.format.transfer.post, galev processing, and movie generation."
---

# PeTar N-body Simulation Skill

## Purpose

Provide reliable, command-level guidance for running PeTar simulations with patterns already used in this repository.
Support multiple installed executables and only suggest options confirmed by the selected binary help output.

## Sources In This Repository

Prefer these examples and docs as the source of truth:

- `sample/star_cluster_plummer_N1k.sh`
- `sample/star_cluster_plummer_N1k_binaries_bse.sh`
- `sample/star_cluster_plummer_N1k_binaries_bse_GalpyMWPot.sh`
- `README.md` sections for OpenMP, MPI, GPU, restart, and options.
- `assets/option-matrix.md` (generated from installed binary help)
- `assets/script-tools.md` (installed script-tool inventory from Makefile.in)

## Installed Script Tools (Current Host)

The repository installs the following workflow tools via `install_script_tool` in `Makefile.in`.

- `petar.init`
- `petar.find.dt`
- `petar.update.par`
- `petar.data.clear`
- `petar.data.process`
- `petar.movie`
- `petar.data.gether`
- `petar.get.object.snap`
- `petar.format.transfer.post`
- `petar.galev.process`
- `petar.external.pot.movie`

These tools are part of the skill surface and should be suggested when the user's request matches their purpose.

## Tool Selection By Task

Use these tools proactively when the user intent matches the task.

- `petar.init`:
  convert raw particle tables into PeTar input snapshots.
- `petar.find.dt`:
  find a suitable tree time step for a given snapshot and launch configuration.
- `petar.update.par`:
  update legacy parameter files from older PeTar versions.
- `petar.data.gether`:
  merge MPI outputs and generate snapshot path lists.
- `petar.data.process`:
  post-process snapshots into single/binary/multiple data, core data, Lagrangian radii, escapers.
- `petar.data.clear`:
  remove post-snapshot events before restarting from an intermediate snapshot.
- `petar.get.object.snap`:
  extract time series of selected objects, IDs, binary IDs, types, or mass ranges from snapshots.
- `petar.format.transfer.post`:
  convert post-processed snapshot formats such as ascii, binary, npy.
- `petar.galev.process`:
  prepare or convert snapshots for Galev-based photometric workflows.
- `petar.movie`:
  generate movies of snapshots, HR diagrams, binary evolution, and Lagrangian evolution.
- `petar.external.pot.movie`:
  generate movies for external potential map evolution from `petar.external` outputs.

If the user asks for one of these tasks, do not answer only with general advice; provide the corresponding tool command pattern.

## Installed Binary Families (Current Host)

This skill is designed to work with installed binaries that start with:

- `petar.mpi.omp.avx512`

Current solver executables on this host are:

- `petar.mpi.omp.avx512`
- `petar.mpi.omp.avx512.agama`
- `petar.mpi.omp.avx512.bse`
- `petar.mpi.omp.avx512.bse.agama`
- `petar.mpi.omp.avx512.bse.galpy`
- `petar.mpi.omp.avx512.galpy`

Installed helper tools with the same prefix family include:

- `*.format.transfer`
- `*.hard.debug`

These helper tools are not simulation executables. Do not use them as the main binary for N-body production runs.

## Binary Capability Discovery (Required)

Before composing a command with custom options, always detect capabilities from help output.

1. Resolve binary path:

```bash
command -v <petar_binary>
```

2. Read help and extract options:

```bash
<petar_binary> -h | awk '{for(i=1;i<=NF;i++){if($i ~ /^-[A-Za-z0-9]$/ || $i ~ /^--[A-Za-z0-9][A-Za-z0-9-]*$/) print $i}}' | sed 's/[,:;]$//' | sort -u
```

3. Validate every requested user option against the extracted list.

4. If an option is missing:
- explain which binary does not support it;
- suggest feature-matched binaries (for example, options with `--bse-` require a `.bse` build; `--galpy-*` requires `.galpy`; `--agama-*` requires `.agama`).

4a. If the requested binary name ends with `.hard.debug` or `.format.transfer`, stop and explain that it is a helper tool rather than a simulation executable. Then switch validation to the matching solver binary without that suffix.

5. Only emit final runnable command after validation passes.

## Binary-Driven Scenario Inference

If the user provides a binary name, infer the default physics stack before asking follow-up questions.

Inference rules:

- `petar.mpi.omp.avx512`
  Default scenario: isolated cluster.
- `*.bse`
  Default scenario: stellar evolution enabled.
- `*.galpy`
  Default scenario: external Galactic potential via Galpy.
- `*.agama`
  Default scenario: external potential via Agama.
- `*.bse.galpy`
  Default scenario: stellar evolution plus Galpy.
- `*.bse.agama`
  Default scenario: stellar evolution plus Agama.

When inferring from binary name:

1. Do not ask whether BSE, Galpy, or Agama is needed if the suffix already determines it.
2. Only ask for the missing physical inputs required by that inferred scenario.
3. If the user explicitly asks for a different physics stack than the binary suffix implies, explain the mismatch and suggest the correct binary.

## Binary-To-Option Family Mapping

Use the selected binary suffix to decide which extra option families are reasonable to propose by default:

- base solver:
  `-u`, `-t`, `-o`, `-p`, `-a`, `-b`, `-r`, `-s`, `--r-search-min`, `--r-group`, `--soft-eps`, `--dt-soft-*`, Hermite and AR options.
- `.bse`:
  all base options plus `--bse-*`, `--stellar-evolution`, `--detect-interrupt`, `--rand-seed`, `--rand-seedfile`.
- `.galpy`:
  all base options plus `--galpy-set`, `--galpy-conf-file`, `--galpy-type-arg`, `--galpy-rscale`, `--galpy-vscale`.
- `.agama`:
  all base options plus `--agama-conf-file`, `--agama-rscale`, `--agama-vscale`.

Do not proactively suggest options from families that the binary cannot support.

## Required Inputs From User

Ask for the minimum set before generating commands:

1. Simulation type: isolated | BSE binaries | galactic tidal field (Galpy) | external potential (Agama) | restart.
2. Unit choice: Henon (`-u 0` default) or astrophysical (`-u 1`).
3. End time (`-t`) and output interval (`-o`).
4. Parallel mode: serial | OpenMP | MPI+OpenMP | GPU.
5. Whether initial data already exists in PeTar snapshot format.
6. Target binary name (for example `petar.mpi.omp.avx512.bse.galpy`).
7. Extra custom options to pass through verbatim after validation.

If the binary name is already given, infer item 1 when possible instead of asking again.

## Minimal Required Question Sets

After binary inference, ask only for the missing inputs required by that scenario.

### Base Isolated Cluster

Required if not already provided:

1. Initial condition source:
  existing PeTar snapshot | raw particle table | generator command such as `mcluster`.
2. Unit mode:
  `-u 0` or `-u 1`.
3. End time `-t`.
4. Output interval `-o`.
5. Parallel launch mode:
  serial | OpenMP | MPI+OpenMP.

Do not ask about BSE, Galpy, or Agama.

### BSE / SSE Scenario

Required if not already provided:

1. Initial condition source.
2. Unit mode.
3. End time `-t`.
4. Output interval `-o`.
5. Number of primordial binaries for `-b`.
6. Metallicity for `--bse-metallicity`.
7. Parallel launch mode.

Do not ask whether stellar evolution should be enabled if the binary already ends with `.bse`.

### Galpy Scenario

Required if not already provided:

1. Initial condition source.
2. Unit mode.
3. End time `-t`.
4. Output interval `-o`.
5. Cluster center phase-space coordinates for `petar.init -c x,y,z,vx,vy,vz` if starting from raw input.
6. Galactic potential choice:
  `--galpy-set`, `--galpy-conf-file`, or `--galpy-type-arg`.
7. If using custom Galpy scaling, `--galpy-rscale` and `--galpy-vscale`.
8. Parallel launch mode.

If the binary is also `.bse`, additionally ask for:

9. Number of primordial binaries.
10. BSE metallicity.

### Agama Scenario

Required if not already provided:

1. Initial condition source.
2. Unit mode.
3. End time `-t`.
4. Output interval `-o`.
5. Cluster center phase-space coordinates for `petar.init -c x,y,z,vx,vy,vz` if starting from raw input.
6. Agama configuration file for `--agama-conf-file`.
7. If using custom scaling, `--agama-rscale` and `--agama-vscale`.
8. Parallel launch mode.

If the binary is also `.bse`, additionally ask for:

9. Number of primordial binaries.
10. BSE metallicity.

### Restart Scenario

Required if not already provided:

1. Restart snapshot filename.
2. Parameter file path, usually `input.par`.
3. Any overridden options after `-p`.
4. Whether outputs should append or overwrite (`-a 1` vs `-a 0`).
5. Parallel launch mode.

## Missing-Input Collection Rule

When invoked:

1. Infer scenario from binary name if possible.
2. Build the minimal required input set for that scenario.
3. Compare it against information already present in the user request.
4. Ask only for missing fields.
5. If enough information is already present, skip questions and emit the command directly.

## Input-Source Driven Workflow Selection

Decide the command flow from the user's input source before generating commands.

### Case A: Raw Particle Table

Definition:

- User has a plain particle file with mass, position, velocity columns.
- Or user provides a generator command such as `mcluster` that produces such a file.

Workflow:

1. Generate or identify the raw file.
2. Use `petar.init` to convert it into a PeTar snapshot.
3. Run the selected solver binary on the generated snapshot.
4. Run scenario-appropriate post-processing.

Default conversion rules:

- If user works in astrophysical units from examples in this repository, prefer `petar.init -v kms2pcmyr -f input ...`.
- If scenario includes BSE, add `-s bse` to `petar.init`.
- If scenario includes Galpy or Agama and the initial condition is in a galactic frame, add `-c ... -t` to `petar.init`.

### Case B: Existing PeTar Snapshot

Definition:

- User already has a PeTar input snapshot such as `input`, `data.0`, or another valid snapshot file.

Workflow:

1. Skip `petar.init`.
2. Run the selected solver binary directly on that snapshot.
3. Run scenario-appropriate post-processing if the user asks for it or if a full workflow is requested.

### Case C: Restart / Resume From Previous Run

Definition:

- User provides a previous output snapshot plus parameter file, or explicitly asks to resume.

Workflow:

1. Use `-p <parameter_file>` before overridden options.
2. Pass any new options after `-p`.
3. Use the restart snapshot as the final positional argument.
4. If requested, set `-a 0` to overwrite outputs instead of appending.

Canonical pattern:

```bash
<launcher> <petar_binary> -p input.par [override options] <restart_snapshot>
```

### Case D: Existing Generator Script In Repository

Definition:

- User wants a workflow matching repository examples such as `mcluster` + `petar.init` + solver + `petar.data.process`.

Workflow:

1. Prefer repository sample script structure.
2. Replace only the scenario-specific placeholders.
3. Preserve repository flag conventions unless the user explicitly overrides them.

## Tool-Driven Workflow Extensions

In addition to solver execution, extend the workflow with the script tools below when the request implies them.

### Timestep Tuning

If the user asks for a good tree step, performance tuning, or "best dt":

```bash
petar.find.dt [options] <snapshot>
```

Use the selected solver name with `-p` or custom launcher prefix with `-r` when needed.

### Output Gathering

If the run used MPI or the user needs a snapshot list file:

```bash
petar.data.gether [options] <data_prefix>
```

Use `-l` when only a snapshot list is needed.

### Restart Cleanup

If the user restarts from a non-final snapshot and wants to avoid duplicated events:

```bash
petar.data.clear -t <restart_time> <data_prefix>
```

### Object Extraction

If the user wants histories of BHs, NSs, IDs, binaries, or selected mass ranges:

```bash
petar.get.object.snap [options] <mode_argument> [mode_argument_b] <snapshot_list>
```

### Format Conversion Of Post-processed Snapshots

If the user wants to convert `petar.data.process` outputs between ascii, binary, and npy:

```bash
petar.format.transfer.post [options] <snapshot_list>
```

### Galev Preparation

If the user wants Galev photometric inputs or conversion:

```bash
petar.galev.process [options] <snapshot_list>
```

### Movie Generation

If the user wants visualization products:

- simulation movie:

```bash
petar.movie [options] <snapshot_list>
```

- external potential movie:

```bash
petar.external.pot.movie [options] <petar.external_parameter_file>
```

## Workflow

1. Validate build/runtime prerequisites.
2. Resolve and verify target binary exists.
3. Infer scenario and option families from binary suffix.
4. Infer input source type: raw input | existing snapshot | restart.
5. Build the scenario-specific minimal required input list.
6. Ask only for missing inputs.
7. Extract supported option list from `<binary> -h`.
8. Validate user custom options.
9. Select one scenario template below.
10. Fill placeholders and output runnable command blocks.
11. Add safety notes for OpenMP/MPI/GPU.
12. Add scenario-appropriate post-processing commands.
13. Add tool-driven follow-up commands when the request includes tuning, gathering, cleanup, extraction, conversion, Galev, or movies.

## Scenario-Specific Default Post-Processing

When the user asks for a full workflow, include post-processing by default unless they explicitly ask for run command only.

### Isolated Cluster

Default post-processing:

```bash
petar.data.gether data
petar.data.process -G 0.00449830997959438 data.snap.lst
```

Optional follow-ups when relevant:

- timestep tuning: `petar.find.dt`
- movie generation: `petar.movie`
- object extraction: `petar.get.object.snap`

### BSE / SSE

Default post-processing:

```bash
petar.data.gether data
petar.data.process -i bse data.snap.lst
```

Optional follow-ups when relevant:

- Galev preparation: `petar.galev.process`
- movie generation with HR diagram: `petar.movie -H`
- object extraction by stellar type: `petar.get.object.snap -m type`

### Galpy

Default post-processing:

```bash
petar.data.gether data
petar.data.process -t galpy --r-escape tidal data.snap.lst
```

Optional follow-ups when relevant:

- external potential movie: `petar.external.pot.movie`
- simulation movie with external potential overlay: `petar.movie --ext-pot`

### BSE + Galpy

Default post-processing:

```bash
petar.data.gether data
petar.data.process -i bse -t galpy --r-escape tidal data.snap.lst
```

Optional follow-ups when relevant:

- Galev preparation: `petar.galev.process`
- movie generation: `petar.movie -H --ext-pot`

### Agama

Default post-processing:

```bash
petar.data.gether data
petar.data.process data.snap.lst
```

### BSE + Agama

Default post-processing:

```bash
petar.data.gether data
petar.data.process -i bse data.snap.lst
```

Optional follow-ups when relevant:

- Galev preparation: `petar.galev.process`
- movie generation: `petar.movie -H`

If the user asks only for the simulation launch command, omit post-processing commands.

## Scenario Templates

### 1) Isolated Cluster

If raw initial condition exists as 7-column mass/position/velocity text:

```bash
# Convert input velocity unit km/s -> pc/Myr and write PeTar snapshot "input"
petar.init -v kms2pcmyr -f input <raw_input_file>

# Run simulation
OMP_NUM_THREADS=<threads> OMP_STACKSIZE=128M petar -u 1 -t <t_end_myr> -o <dt_out_myr> input > output 2>&1

# Gather and process outputs
petar.data.gether data
petar.data.process -G 0.00449830997959438 data.snap.lst
```

If a non-default executable is requested, replace `petar` by the selected binary name.

If the user already has a PeTar snapshot, skip the `petar.init` line and keep the solver plus post-processing only.

### 2) Stellar Evolution With BSE/SSE

Require configure-time support:

```bash
./configure --with-interrupt=bse [other options]
make && make install
```

Run template:

```bash
petar.init -s bse -v kms2pcmyr -f input <raw_input_file>
OMP_STACKSIZE=128M petar -u 1 -b <n_primordial_binaries> --bse-metallicity <Z> -t <t_end_myr> -o <dt_out_myr> input > output 2>&1
petar.data.gether data
petar.data.process -i bse data.snap.lst
```

Compatibility rule:

- any option starting with `--bse-`, `--stellar-evolution`, `--detect-interrupt` requires a `.bse` executable.

If the input is already a PeTar snapshot prepared for stellar evolution, skip `petar.init`.

### 3) External Galactic Potential (Galpy)

Require configure-time support:

```bash
./configure --with-interrupt=bse --with-external=galpy [other options]
make && make install
```

Run template:

```bash
# -t in petar.init enables external potential info in initial snapshot
petar.init -c <x,y,z,vx,vy,vz_in_galactic_frame> -t -s bse -v kms2pcmyr -f input <raw_input_file>
OMP_STACKSIZE=128M petar -u 1 -b <n_primordial_binaries> --galpy-set MWPotential2014 --bse-metallicity <Z> -t <t_end_myr> -o <dt_out_myr> input > output 2>&1
petar.data.gether data
petar.data.process -i bse -t galpy --r-escape tidal data.snap.lst
```

Compatibility rule:

- any option starting with `--galpy-` requires a `.galpy` executable.

If the input is already a PeTar snapshot with external-potential context prepared, skip `petar.init` and keep the solver plus post-processing.

### 3b) External Potential (Agama)

Run template:

```bash
petar.init -c <x,y,z,vx,vy,vz_in_galactic_frame> -t -v kms2pcmyr -f input <raw_input_file>
OMP_STACKSIZE=128M <petar_binary_with_agama> -u 1 -t <t_end_myr> -o <dt_out_myr> --agama-conf-file <agama_config_file> input > output 2>&1
petar.data.gether data
petar.data.process data.snap.lst
```

Compatibility rule:

- any option starting with `--agama-` requires an `.agama` executable.

If the input is already a PeTar snapshot with external-potential context prepared, skip `petar.init`.

### 4) Restart / Resume

```bash
# Keep old parameters from input.par, override selected options after -p
petar -p input.par -t <new_end_time> <snapshot_file>
```

To overwrite instead of append on restart, include `-a 0`.

For restart workflows, do not emit `petar.init`.

If restarting from a non-final snapshot, recommend `petar.data.clear` before the restart command when duplicate events are a concern.

## Parallel And Performance Safety Notes

Always include these reminders when applicable:

- OpenMP: set `OMP_STACKSIZE=128M`; ensure `ulimit -s unlimited`.
- MPI+OpenMP launch pattern:
  `OMP_STACKSIZE=128M OMP_NUM_THREADS=<threads> mpiexec -n <n_mpi> petar [options] <snapshot>`
- If MPI binds one core per rank, try `mpiexec --bind-to none`.
- GPU: one MPI rank typically drives one GPU job; tune ranks/threads to avoid CUDA OOM.
- `*.hard.debug` binaries are diagnostic tools.
- `*.format.transfer` binaries are format-conversion helpers.
- Neither `*.hard.debug` nor `*.format.transfer` should be used as the main simulation executable.

## Response Format When Skill Is Invoked

Return output in this fixed structure:

1. Chosen scenario and assumptions.
2. Selected binary, inferred feature stack, and inferred input-source workflow.
3. Exact command block(s), ready to copy.
4. Why each critical option is used (`-u`, `-b`, `-t`, `-o`, `-p`, `--galpy-set`, `--agama-conf-file`, `--bse-metallicity`, etc.).
5. Tool-driven follow-up commands when relevant (`petar.find.dt`, `petar.data.gether`, `petar.data.clear`, `petar.movie`, etc.).
6. Validation checklist (expected output files and quick sanity checks).
7. Optional next-step optimization knobs (`-s`, `-r`, `--r-search-min`, `--r-bin`).

## Guardrails

- Do not invent unsupported PeTar options.
- If user intent conflicts with current build (for example, asks for Galpy but binary lacks it), provide reconfigure command first.
- Prefer repository sample scripts and README wording over generic N-body advice.
- Never pass through user-provided extra options without checking `<selected_binary> -h` first.
- If the user selects a helper tool binary, explicitly switch to the nearest matching solver binary before generating the run command.
- If the binary suffix already determines the physics stack, do not ask redundant questions about whether to enable that same feature.
- Do not ask the full generic input list when the scenario-specific minimal set is smaller.
- If enough inputs are already present in the user request, generate the command immediately instead of asking confirmation questions.
- Do not emit `petar.init` for existing snapshots or restart workflows.
- Do not omit post-processing when the user asks for a complete end-to-end workflow.
- When the user asks for analysis, conversion, movie generation, or restart cleanup, use the installed script tools instead of describing the task abstractly.
