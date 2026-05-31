---
name: petar-nbody-simulation
description: "Use when: setting up or running PeTar N-body simulations, including petar.select binary-family switching, petar.init conversion, isolated clusters, BSE/SSE, Galpy or Agama external potential, MPI/OpenMP/GPU launch, restart/resume, petar.data post-processing, timestep tuning with petar.find.dt, output gathering with petar.data.gether, restart cleanup with petar.data.clear, snapshot extraction with petar.get.object.snap, format conversion with petar.format.transfer.post, galev processing, and movie generation."
---

# PeTar N-body Simulation Skill

## Purpose

Provide reliable, command-level guidance for running PeTar simulations with patterns already used in this repository.
Support multiple installed executables and only suggest options confirmed by the selected binary help output.
Treat the guardrails in this file as mandatory, not advisory.

## Hard Constraints

- Before any simulation run, first infer the required physics and runtime features from the user description, then select the matching PeTar binary family with `petar.select`.
- Do not directly execute whatever `petar` currently resolves to unless it was chosen through the selection step for the requested simulation.
- Before execution, ensure all scenario-required inputs are present; if any required input is missing, ask for it explicitly instead of guessing defaults.
- After required inputs are complete, provide a parameter summary and the exact command(s), then ask user confirmation before executing any run command.
- If output prefix/model name is not provided, default to PeTar prefix `data`.
- If unit mode is not provided, default to astrophysical unit mode `-u 1` (Msun, pc, pc/Myr).
- Unit-system self-consistency is mandatory: PeTar runtime does not independently re-scale mass/length/time/velocity; only `G` changes with the selected unit mode, so IC units must be mutually consistent with the chosen `G`.
- If the user requests stellar evolution, confirm the exact interruption module (`bse`, `bseEmp`, `mobse`/`moBSE`, or `dsm`) before selection, and use that module as a required token in `petar.select`.
- If the user requests external potential, confirm `galpy` or `agama` before selection, and use the confirmed mode as a required token in `petar.select`.
- For star-cluster simulations that require IC generation, do not proceed until the star-cluster IC parameter set is complete.
- For star-cluster IC generation, verify generator availability first: prefer `mcluster_gpu` when available, otherwise use `mcluster`.
- If neither `mcluster_gpu` nor `mcluster` is available, stop and ask the user to install one before proceeding.
- For IC produced by `mcluster`, do not assume a fixed unit system independent of generator options: `mcluster` output depends on its own `-u` setting. In the common `mcluster -u 1` case, the output is astrophysical-style data with mass in `Msun`, position in `pc`, and velocity in `km/s`; if this IC is later used with `petar -u 1` (target unit: `Msun`, `pc`, `pc/Myr`), `petar.init -v 1.022712165045695` must be applied to convert `km/s -> pc/Myr`.
- For any raw IC source (including `mcluster`), explicitly confirm or infer source mass/length/velocity units before `petar.init`. Recommended path unless the user explicitly requests otherwise: normalize the IC into the target PeTar unit system during `petar.init` using conversion options such as `-r` and `-v` so that, for `petar -u 1`, the final IC is in `Msun`, `pc`, `pc/Myr` and uses the standard astrophysical `G`. Advanced path: keep the native length scale (for example `kpc`) and convert velocity consistently (for example `km/s -> kpc/Myr`), but then you must also provide the matching `G` via `petar -G` and ensure every unit-sensitive physics module and post-processing option is converted consistently as well (for example BSE, Galpy, Agama scaling options). Because this path is error-prone, do not recommend it unless the user explicitly asks to preserve the native units.
- For BSE-related stellar-evolution runs (`bse`, `bseEmp`, `mobse`), metallicity is mandatory and must be confirmed before execution.
- For external-potential runs, the potential model and the simulation-object center phase-space coordinates in that potential are mandatory; do not pick a potential model on behalf of the user.
- If the user does not explicitly request debug mode, selected solver candidates must exclude debug families (for example suffix token `g` or assert-debug builds).
- If the user does not explicitly request GPU, selected solver candidates must exclude `gpu` families.
- Do not emit a runnable solver command until the selected binary has been capability-checked with `-h`.
- Do not pass through custom user options unless they are validated against the selected binary help output.
- Do not use helper binaries such as `*.hard.debug` or `*.format.transfer` as the main simulation executable.
- If source-level debugging is requested (for example, gdb backtrace with file/line symbols), require a rebuild with `--with-debug=g`; do not treat `petar.hard.debug` as a substitute for debug-symbol builds.
- Do not recommend manual symlink edits when `petar.select` can perform the family switch.
- Do not skip `petar.data.gether` after MPI runs when downstream tools need merged outputs or a snapshot list.
- Do not emit `petar.init` for restart or resume workflows.
- Do not treat analysis, conversion, movie generation, or restart cleanup as abstract guidance when an installed PeTar tool exists for the task; use the tool command pattern.
- For post-processing tools (`petar.data.process`, `petar.movie`, `petar.format.transfer.post`, `petar.galev.process`, `petar.get.object.snap`), required mode parameters must match the currently selected solver family and data format settings.

## Sources In This Repository

Prefer these examples and docs as the source of truth:

- `sample/star_cluster_plummer_N1k.sh`
- `sample/star_cluster_plummer_N1k_binaries.sh`
- `sample/star_cluster_plummer_N1k_binaries_bse.sh`
- `sample/star_cluster_plummer_N1k_GalpyMWPot.sh`
- `sample/star_cluster_plummer_N1k_binaries_bse_GalpyMWPot.sh`
- `sample/star_cluster_plummer_N1k_AgamaMWPotHunter24.sh`
- `sample/data_analysis.ipynb` (Python post-processing and plotting patterns)
- `README.md` sections for OpenMP, MPI, GPU, restart, and options.
- `assets/option-matrix.md` (generated from installed binary help)
- `assets/script-tools.md` (installed script-tool inventory from Makefile.in)
- `HANDOFF.md` (cross-machine continuation notes)
- `assets/prompt-starters.md` (chat prompt templates for resuming work)

## Installed Script Tools (Current Host)

The repository installs the following workflow tools via `install_script_tool` in `Makefile.in`.

- `petar.init`
- `petar.select`
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

Additional workflow utilities frequently available on host installations include:

- `petar.galpy.pot.movie`
- `petar.get.init.binary`
- `petar.galpy.help`
- `petar.external.galpy`
- `petar.external.agama`

These tools are part of the skill surface and should be suggested when the user's request matches their purpose.

## Tool Selection By Task

Use these tools proactively when the user intent matches the task.

- `petar.init`:
  convert raw particle tables into PeTar input snapshots.
- `petar.select`:
  switch installed symlinks `petar`, `petar.hard.debug`, and `petar.format.transfer` to a selected installed binary family.
  Support direct target mode (`petar.select <suffix|binary-name>`), feature auto-select mode (`petar.select --require ... [--optional ...]`), and listing (`petar.select --list`).
  In feature mode, all required tokens must match; optional tokens are used for ranking.
  Module disambiguation is mandatory before feature mode selection:
  - stellar evolution request: ask or infer one of `bse`, `bseEmp`, `mobse` (accept user wording `moBSE` as `mobse`), `dsm`;
  - external potential request: ask or infer one of `galpy`, `agama`.
  Then pass the confirmed module tokens in `--require`.
  For every simulation request, determine the needed feature set first, then use `petar.select` to choose the matching family before generating any run command.
  For every binary-family request, prefer `petar.select` over manual symlink edits, and always stop with configure hints if no feature-matched family exists.
  Treat physics/structure-sensitive families as require-only: interrupt (`base`, `bse`, `mobse`, `bseEmp`, `dsm`), external (`galpy`, `agama`), external-hard (`gasdrag`), `pn*`, and `mpfrc` (suffix token `mp`).
  In other words, candidates containing these tokens are excluded unless explicitly requested in `--require`.
  If no match exists, surface configure hints mapped from required features (for example, `bse -> --with-interrupt=bse`, `galpy -> --with-external=galpy`).
  Unknown feature handling: `--require` must fail fast; `--optional` should warn and ignore unsupported tokens.
- `petar.hard.debug`:
  helper tool for replaying and diagnosing hard-integrator dump files (for example, `[output_prefix].hard_dump.*`), not a full N-body production solver.
- `petar.find.dt`:
  find a suitable tree time step for a given snapshot and launch configuration.
  This is useful for star-cluster scale systems, but for very small-N setups
  (for example isolated binaries or a few-body test) do not trust the automatic
  estimate blindly: the tree time step can be badly misjudged when the system
  contains only a few particles. In that regime, choose the tree time step
  manually from the orbital timescale and make sure the binary orbit is still
  resolved with several steps per orbit.
- `petar.update.par`:
  update legacy parameter files from older PeTar versions.
- `petar.data.gether`:
  merge MPI outputs and generate snapshot path lists.
  If the workflow used MPI and any downstream post-processing, movie, or extraction step needs consolidated outputs, `petar.data.gether` is mandatory before those steps.
  By default, do not assume group files should be gathered.
  If the user explicitly wants merged group outputs, use `petar.data.gether -g <prefix>`.
  This matters because MPI ranks may write separate `data.group.*.nX` files and merged group outputs can be large.
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
- `petar.galpy.pot.movie`:
  generate movies for Galpy potential map evolution.
- `petar.get.init.binary`:
  generate primordial binary pairing tables for BSE initialization workflows.
- `petar.galpy.help`:
  inspect Galpy potential models and generate option/configuration guidance.
- `petar.external.galpy` / `petar.external.agama`:
  generate external potential map snapshots from run parameters (for example before `petar.external.pot.movie` or `petar.movie --ext-pot`).

If the user asks for one of these tasks, do not answer only with general advice; provide the corresponding tool command pattern.

For binary-family switching requests, prefer `petar.select` over manually editing symlinks.
When the request contains physics requirements (for example, `bse`, `galpy`, `agama`), use `--require` first and keep performance-related tokens (for example, `mpi`, `omp`, `avx512`, `avx2`) in `--optional`.

## Installed Binary Families (Current Host)

Do not assume a fixed install directory or a fixed SIMD family.

Binary discovery rule:

1. Scan commands in `PATH` whose executable name starts with `petar`.
2. Treat binaries ending with `*.hard.debug` or `*.format.transfer` as helper tools.
3. Treat binaries exposing core runtime options `-u`, `-t`, and `-o` in `-h` as solver binaries.
4. Use `.github/skills/petar-nbody-simulation/assets/option-matrix.md` as the machine-local inventory snapshot.

Requirement-driven selection rule:

1. Do not assign a fixed priority across physics scenarios; scenario choice depends on user intent.
2. Map user intent to configure features first:
  - `--with-interrupt` controls interruption module family (`base`, `bse`, `mobse`, `bseEmp`, `dsm`).
  - `--with-external` controls long-timescale external potential in tree steps (`galpy`, `agama`).
  - `--with-external-hard` controls short-timescale external forces in hard integrators (`gasdrag`).
  - `--with-pn` controls post-Newtonian relativistic corrections (`pn*`).
  - `--enable-mpfrc` changes position representation precision and affects snapshot structure compatibility (`mp` suffix token).
3. Filter binaries by required feature suffixes first using `assets/option-matrix.md` section `Requirement-Driven Solver Filters`.
4. Apply mandatory pre-filters before ranking:
  - if user did not explicitly request debug mode, exclude debug families (`g`, assert-debug variants);
  - if user did not explicitly request GPU, exclude `gpu` families.
5. After filtering, rank candidates by performance suffix priority:
  - if GPU is explicitly requested: `gpu` > `avx512` > `avx2` > `omp` > `mpi`;
  - otherwise: `avx512` > `avx2` > `omp` > `mpi`.
6. If no candidate remains, do not emit a run command. Ask whether to proceed with `./configure` + `make install`, provide the exact suggested commands, and execute only after explicit user confirmation.
7. Still validate requested options against `<selected_binary> -h` before emitting commands.

Special-purpose modules (on-demand, not part of default recommendation ranking):

- `--enable-64b`:
  enable 64-bit floating-point particle-tree force calculation (default tree-force path is lower precision).
- `--enable-mpfrc`:
  represent particle positions with split high-precision components to reduce round-off errors in large dynamic-range setups (for example, star clusters embedded in galactic environments with extreme scale separation).
- `--enable-gperf`:
  enable gperftools profiling for runtime performance investigation.
- `--with-debug=g`:
  debug-focused build (`-O0`, `-g`) for debugger workflows such as `gdb`.
- `--with-debug=assert`:
  assertion-focused debug mode with optimization retained to avoid severe runtime slowdown.

Step-mode selection rule:

- `--with-step-mode` controls tree-step integration mode and must be validated from current `./configure -h` before suggesting reconfigure commands.
- Treat supported modes as configure-driven build features, not runtime flags.
- For the current T1 performance workflow in this repository, use explicit `kdk` vs `kdkdk4` binaries as the primary forced comparison.
- If `kdkdk4` is missing on the current machine, say so explicitly and only then suggest temporary fallback to `kdkdk`, together with a warning that the comparison target changed.

When these features are requested by the user, prioritize feature-matched binaries or provide reconfigure guidance; otherwise keep them outside default solver recommendations.

Suffix tokens should be interpreted as configure-driven feature combinations from `configure.ac`, including categories such as:

- parallel/runtime (`mpi`, `omp`, `gpu`)
- architecture (`avx`, `avx2`, `avx512`, `64b`)
- interruption mode (`base`, `bse`, `mobse`, `bseEmp`, `dsm`)
- external potential (`galpy`, `agama`)
- external hard-force (`gasdrag`)
- post-Newtonian (`pn*`)
- high-precision tree force (`64b`)
- high-precision position (`mp` from `--enable-mpfrc`)
- debug suffix (`g`, `d`)

Helper tools include (but are not limited to):

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

4a. If the requested binary name ends with `.hard.debug` or `.format.transfer`, stop and explain that it is a helper tool rather than a simulation executable. Then switch validation to the matching solver binary without that suffix for production-run command generation.

4b. If source-level debugging is requested and the current binary is not built with `--with-debug=g`, stop and request reconfigure + rebuild with `--with-debug=g` before running debugger workflows.

5. Only emit final runnable command after validation passes.

## Binary-Driven Scenario Inference

If the user provides a binary name, infer the default physics stack before asking follow-up questions.

Inference rules:

- `petar` (or any solver without explicit physics suffix)
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
- `*.gasdrag`
  Default scenario extension: short-timescale external hard force is enabled.
- `*.pn*`
  Default scenario extension: post-Newtonian relativistic corrections are enabled.

When inferring from binary name:

1. Do not ask whether BSE, DSM, Galpy, or Agama is needed if the suffix already determines it.
2. Only ask for the missing physical inputs required by that inferred scenario.
3. If the user explicitly asks for a different physics stack than the binary suffix implies, explain the mismatch and suggest the correct binary.

## Binary-To-Option Family Mapping

Use the selected binary suffix to decide which extra option families are reasonable to propose by default:

- base solver:
  `-u`, `-t`, `-o`, `-p`, `-a`, `-b`, `-r`, `-s`, `--r-search-min`, `--r-group`, `--soft-eps`, `--dt-soft-*`, Hermite and AR options.
- `.bse`:
  all base options plus `--bse-*`, `--stellar-evolution`, `--detect-interrupt`, `--rand-seed`, `--rand-seedfile`.
- `.dsm`:
  all base options plus `--dsm-*` option family and DSM-specific interrupt controls.
- `.galpy`:
  all base options plus `--galpy-set`, `--galpy-conf-file`, `--galpy-type-arg`, `--galpy-rscale`, `--galpy-vscale`.
- `.agama`:
  all base options plus `--agama-conf-file`, `--agama-rscale`, `--agama-vscale`.
- `.gasdrag`:
  hard-part external force is compiled in (short-timescale), typically from `--with-external-hard=gasdrag`.
- `.pn*`:
  post-Newtonian force terms are compiled in, typically from `--with-pn=<mode>`.

Do not proactively suggest options from families that the binary cannot support.

## Required Inputs From User

Ask for the minimum set before generating commands:

1. Simulation execution path (working directory where command will run).
2. Simulation type: isolated | BSE binaries | DSM mode | galactic tidal field (Galpy) | external potential (Agama) | restart.
3. Unit choice: Henon (`-u 0`) or astrophysical (`-u 1`, default when omitted).
4. End time (`-t`) and output interval (`-o`).
5. Parallel mode: serial | OpenMP | MPI+OpenMP | GPU.
6. Whether initial data already exists in PeTar snapshot format.
7. Target binary name (for example `petar.mpi.omp.avx2.bse.galpy`).
8. Extra custom options to pass through verbatim after validation.

If the binary name is already given, infer item 1 when possible instead of asking again.

Default fill-ins when omitted:

- output prefix/model name: `data`.
- unit mode: `-u 1` (Msun, pc, pc/Myr).

### Star-Cluster IC Parameters (Required When IC Must Be Generated)

If the simulation object is a star cluster and IC is not already an existing PeTar snapshot, collect all of the following before emitting the run command:

1. One size scale for the stellar system:
  total mass or total star count.
2. Half-mass radius.
3. Density profile model:
  for example Plummer or King; if a model requires extra parameters (for example King W0), collect them explicitly.
4. IMF model.
5. Binary fraction and binary distribution model.

Do not assume defaults for these items when the user requests a generated star-cluster setup.

### Star-Cluster Generator Availability Rule

When the workflow needs generator-based star-cluster IC creation:

1. Check generator availability in this order:
  `command -v mcluster_gpu` then `command -v mcluster`.
2. If `mcluster_gpu` exists, use it by default for faster generation.
3. If `mcluster_gpu` is unavailable but `mcluster` exists, use `mcluster`.
4. If both are unavailable, do not continue with generation commands; inform the user that generator installation is required and ask whether to proceed with installation guidance.

## Minimal Required Question Sets

After binary inference, ask only for the missing inputs required by that scenario.

### Base Isolated Cluster

Required if not already provided:

1. Initial condition source:
  existing PeTar snapshot | raw particle table | generator command such as `mcluster`.
1a. Execution path (working directory).
2. Unit mode:
  `-u 0` or `-u 1` (if omitted, use `-u 1`).
3. End time `-t`.
4. Output interval `-o`.
5. Parallel launch mode:
  serial | OpenMP | MPI+OpenMP.

Do not ask about BSE, Galpy, or Agama.

If the source is raw/generator and object is a star cluster, also require `Star-Cluster IC Parameters` before emitting run commands.

### BSE / SSE Scenario

Required if not already provided:

1. Initial condition source.
1a. Execution path.
2. Unit mode.
  If omitted, use `-u 1`.
3. End time `-t`.
4. Output interval `-o`.
5. Number of primordial binaries for `-b`.
6. Metallicity for `--bse-metallicity`.
7. Parallel launch mode.

Do not ask whether stellar evolution should be enabled if the binary already ends with `.bse`.

Primordial-binary counting rule (avoid a common misunderstanding):

- `N` in generator commands and the PeTar input header is the total number of stars (single stars + binary members), not the number of systems.
- A binary fraction like `mcluster -b 0.95` means about 95% of stars are arranged into pairs, so the primordial binary count is approximately `n_bin ~= 0.95 * N / 2`.
- Therefore, a header value such as `500` still means 500 stars in total; with 95% binary fraction this implies about 237-238 binaries (about 474-476 stars in binaries), not 500 binaries.
- Keep `petar -b <n_bin>` consistent with the actual generated pairing count (read from generator log if available), while keeping header `N` as the total star count.

### DSM Scenario

Required if not already provided:

1. Initial condition source.
1a. Execution path.
2. Unit mode.
  If omitted, use `-u 1`.
3. End time `-t`.
4. Output interval `-o`.
5. DSM key controls to override (only those requested), for example `--dsm-seed-mass`, `--dsm-he-disk`, `--dsm-lambda0`, `--dsm-dt-factor`.
6. Parallel launch mode.

If the binary is also `.galpy` or `.agama`, additionally ask for the corresponding external-potential inputs.

If raw input is converted with `petar.init`, use `-s dsm` for DSM-compatible stellar columns.

### Galpy Scenario

Required if not already provided:

1. Initial condition source.
1a. Execution path.
2. Unit mode.
  If omitted, use `-u 1`.
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

Do not select a Galpy potential model implicitly; require user-provided model choice.

### Agama Scenario

Required if not already provided:

1. Initial condition source.
1a. Execution path.
2. Unit mode.
  If omitted, use `-u 1`.
3. End time `-t`.
4. Output interval `-o`.
5. Cluster center phase-space coordinates for `petar.init -c x,y,z,vx,vy,vz` if starting from raw input.
6. Agama potential model/configuration for `--agama-conf-file`.
7. If using custom scaling, `--agama-rscale` and `--agama-vscale`.
8. Parallel launch mode.

If the binary is also `.bse`, additionally ask for:

9. Number of primordial binaries.
10. BSE metallicity.

Do not select an Agama potential model implicitly; require user-provided model choice.

### Restart Scenario

Required if not already provided:

1. Restart snapshot filename.
1a. Execution path.
2. Parameter file path, usually `<output_prefix>.par` (default `data.par`).
3. Any overridden options after `-p`.
4. Whether outputs should append or overwrite (`-a 1` vs `-a 0`).
5. Parallel launch mode.

## Missing-Input Collection Rule

When invoked:

1. Infer scenario from binary name if possible.
2. Build the minimal required input set for that scenario.
3. Compare it against information already present in the user request.
4. Ask only for missing fields.
5. If enough information is already present, output a parameter summary and runnable command block first, then ask for execution confirmation before running.

## Pre-Execution Confirmation Rule

Before any actual run execution:

1. Present a compact parameter summary (scenario, binary family, execution path, key physics toggles, unit mode, output prefix/model name, end time/output interval, and critical module inputs).
2. Present the exact command block(s) to run.
3. Ask for explicit user confirmation.
4. Execute only after confirmation.

## Output Prefix Selection Rule (`-f`)

Choose output prefix with the default-first rule.

1. If user explicitly provides `-f <prefix>`, use it unchanged.
2. If workflow is a restart/resume from previous outputs, do not force default `data`; keep the restart context and choose prefix explicitly if needed.
3. If user does not specify a prefix and workflow is not restart, default to prefix `data`.
4. If user does not specify a prefix and existing `data*` outputs are present, keep default `data` and explicitly warn about overwrite/append behavior; if user asks to avoid collisions, propose an alternative prefix (for example `data.<tag>`) and confirm.

Downstream commands must use the same chosen prefix consistently, for example:

```bash
petar ... -f <prefix> ...
petar.data.gether <prefix>
petar.data.process <prefix>.snap.lst
```

## Snapshot Read-Mismatch Detection And Recovery

Apply this rule to Python-based snapshot readers, especially:

- `petar.movie`
- `petar.data.process`
- `petar.format.transfer.post`
- `petar.galev.process`
- `petar.get.object.snap`

Treat the following messages as evidence that snapshot reading is misconfigured rather than merely noisy warnings:

1. Binary misalignment warnings, for example:
   `Binary file size is not aligned with dtype itemsize`
2. ASCII column mismatch warnings, for example:
   `The reading data shape or the number of columns mismatches the number of columns`
3. Text decoding errors such as `utf-8` decode failures.

Apply the same rule to `data.lagr` and other `petar.data.process` outputs when the reader reports binary size, dtype, or alignment mismatches: those are read failures, not harmless warnings.

Interpretation:

- binary misalignment:
  snapshot structure or reader mode is wrong.
- ASCII shape/column mismatch:
  snapshot format or particle schema is wrong.
- `utf-8` decode failure:
  reader is trying to read binary data as ascii/text, or wrong format was selected.

Recovery procedure:

1. Stop the current processing command first if it is still running.
2. Re-check the selected solver family and the producing command parameters (`interrupt mode`, `external mode`, `snapshot format`) before retrying.
3. Re-check the tool help output before retrying.
4. Retry by correcting reader parameters rather than ignoring the warning.

Most common fixes to try:

- snapshot format selection:
  `-s ascii | binary | npy`
- snapshot origin/type selection when available:
  `--snapshot-type origin | post | generate_binary`
- interrupt mode selection:
  `-i none | base | bse | mobse | bseEmp | dsm`
- external mode selection:
  `-t none | galpy | agama`
- use the exact spelling shown by current `-h`; do not assume aliases like `no` are valid when help says `none`.

Decision rule:

- If the warning disappears after parameter correction, continue with the corrected command only.
- If warnings persist after reasonable parameter corrections, stop and warn the user explicitly that snapshot reading is likely misaligned and downstream products may be invalid.
- If warnings persist, abort downstream analysis/conversion/movie generation in this turn and report the unresolved mismatch to the user; do not continue with partial outputs.
- Do not present movies, post-processing results, extracted objects, or converted files as trustworthy when these warnings remain unresolved.

## Input-Source Driven Workflow Selection

Decide the command flow from the user's input source before generating commands.

### Case A: Raw Particle Table

Definition:

- User has a plain particle file with mass, position, velocity columns.
- Or user provides a generator command such as `mcluster` that produces such a file.

Workflow:

1. Generate or identify the raw file.
  If generation is required, apply `Star-Cluster Generator Availability Rule` first and prefer `mcluster_gpu` when available.
2. Use `petar.init` to convert it into a PeTar snapshot.
3. Run the selected solver binary on the generated snapshot.
4. Run scenario-appropriate post-processing.

Default conversion rules:

- If user works in astrophysical units from examples in this repository, prefer `petar.init -v kms2pcmyr -f input ...`.
- If scenario includes BSE, add `-s bse` to `petar.init`.
- If scenario includes Galpy or Agama and the initial condition is in a galactic frame, add `-c ... -t` to `petar.init`.
- If the source was generated with a high binary fraction (for example `mcluster -b 0.95`), do not reinterpret header `N` as binary count; `N` remains total stars and `petar -b` should use the actual primordial binary pair count.

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
5. If the restart starts from a non-final snapshot and duplicate events are possible, recommend `petar.data.clear` before the restart command.

Canonical pattern:

```bash
<launcher> <petar_binary> -p data.par [override options] <restart_snapshot>
```

### Case D: Existing Generator Script In Repository

Definition:

- User wants a workflow matching repository examples such as `mcluster` + `petar.init` + solver + `petar.data.process`.

Workflow:

1. Prefer repository sample script structure.
2. Replace only the scenario-specific placeholders.
3. If generator command is part of the script, prefer `mcluster_gpu` over `mcluster` when available.
4. Preserve repository flag conventions unless the user explicitly overrides them.

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

Group-data rule:

- default behavior: do not add `-g` automatically.
- if user explicitly asks to gather group files, use:

```bash
petar.data.gether -g <data_prefix>
```

- explain that this merges per-rank group outputs such as `data.group.*.n2`, `data.group.*.n3`, etc., and may create large files.

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

- For BSE-enabled scenarios (`-i bse`), default to `-c logtemp` for particle color unless the user explicitly requests another color mode.

- simulation movie:

```bash
petar.movie [options] <snapshot_list>
```

- external potential movie:

```bash
petar.external.pot.movie [options] <petar.external_parameter_file>
```

### External Potential Map Generation

If the user asks to generate potential maps from runtime parameters before visualization:

```bash
petar.external.<galpy|agama> -p <output_prefix>.par -m <pot_conf>
```

Typical use:

1. create a map configuration file (`pot_conf`) describing time/grid ranges;
2. generate map snapshots with `petar.external.agama` or `petar.external.galpy`;
3. visualize with `petar.external.pot.movie` or overlay in `petar.movie --ext-pot`.

## Python Data Analysis Templates (from `sample/data_analysis.ipynb`)

When the user asks for analysis of PeTar outputs, provide runnable Python code snippets (not only shell commands), reusing these patterns.

### Minimal Setup

```python
import petar
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
```

### Read Original Snapshots

- pure gravity:

```python
filename = 'data.10'
header = petar.PeTarDataHeader(path + filename)
particle = petar.Particle()
particle.fromfile(path + filename, offset=petar.HEADER_OFFSET)
```

- SSE/BSE compiled:

```python
particle = petar.Particle(interrupt_mode='bse')
particle.fromfile(path + 'data.10', offset=petar.HEADER_OFFSET)
```

- Galpy external potential:

```python
header = petar.PeTarDataHeader(path + 'data.10', external_mode='galpy')
particle = petar.Particle(interrupt_mode='bse', external_mode='galpy')
particle.fromfile(path + 'data.10', offset=petar.HEADER_OFFSET_WITH_CM)
```

### Read Post-processed Outputs (`petar.data.process`)

```python
single = petar.Particle()
single.fromfile(path + 'data.10.single')

binary = petar.Binary(member_particle_type=petar.Particle,
            interrupt_mode='bse', external_mode='galpy',
            G=petar.G_MSUN_PC_MYR)
binary.fromfile(path + 'data.10.binary')

lagr = petar.LagrangianMultiple()
lagr.fromfile(path + 'data.lagr')
```

### Read Group Files (`petar.data.gether -g`)

```python
g2 = petar.GroupInfo(N=2); g2.fromfile(path + 'data.group.n2')
g3 = petar.GroupInfo(N=3); g3.fromfile(path + 'data.group.n3')
g4 = petar.GroupInfo(N=4); g4.fromfile(path + 'data.group.n4')
```

### Read Stellar-Evolution Event Files

```python
sse_type = petar.SSETypeChange();    sse_type.loadtxt(path + 'data.sse.type_change')
sse_kick = petar.SSESNKick();        sse_kick.loadtxt(path + 'data.sse.sn_kick')
bse_type = petar.BSETypeChange();    bse_type.loadtxt(path + 'data.bse.type_change')
bse_kick = petar.BSEKick();          bse_kick.loadtxt(path + 'data.bse.sn_kick')
bse_gw_kick = petar.BSEKick();       bse_gw_kick.loadtxt(path + 'data.bse.gw_kick')
bse_dyn = petar.BSEDynamicMerge();   bse_dyn.loadtxt(path + 'data.bse.dynamic_merge')

merger = petar.BSEMerge()
merger.combine(bse_type, bse_dyn)
```

### Read `petar.get.object.snap` Outputs (time series)

```python
bms = petar.Binary(member_particle_type=petar.Particle,
           interrupt_mode='bse', G=petar.G_MSUN_PC_MYR)
bms.addNewMember('time', np.array([], dtype=float))
bms.fromfile(path + 'object.MS.MS.binary')

p1 = petar.Particle(interrupt_mode='bse')
p1.addNewMember('time', np.array([], dtype=float))
p1.fromfile(path + 'object.1')
```

### Read `data.status` for Few-body `-w 2` Runs

```python
status = petar.Status(external_mode='galpy', N_particle=1)
status.fromfile(path + 'data.status')
pos = status.particles.p0.pos
```

### Plot Templates

- HR diagram (BSE):

```python
lum = particle.star.lum
temp = 5778 * (particle.star.lum / (particle.star.rad * particle.star.rad))**0.25
fig, ax = plt.subplots(1, 1)
norm = colors.BoundaryNorm(boundaries=np.arange(16), ncolors=256)
pt = ax.scatter(temp, lum, c=particle.star.type, cmap='rainbow', norm=norm)
plt.colorbar(pt, ax=ax, label='SSE stellar type')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(30000, 1000); ax.set_ylim(1e-5, 1e6)
ax.set_xlabel('Temperature'); ax.set_ylabel('Luminosity')
```

- Lagrangian radii evolution:

```python
fig, ax = plt.subplots(1, 1)
for i in range(5):
  ax.plot(lagr.time, lagr.all.r[:, i], '-', label=lagr.initargs['mass_fraction'][i])
ax.plot(lagr.time, lagr.all.r[:, -1], '--', label='Rc')
ax.set_xlabel('Time'); ax.set_ylabel('R'); ax.set_yscale('log'); ax.legend()
```

- binary semi-ecc scatter:

```python
x = binary.semi * 206265.0  # AU
y = binary.ecc
q = np.minimum(binary.p1.mass, binary.p2.mass) / np.maximum(binary.p1.mass, binary.p2.mass)
fig, ax = plt.subplots(1, 1)
pt = ax.scatter(x, y, s=binary.mass, c=q)
ax.set_xscale('log')
ax.set_xlabel('Semi-major axes [AU]'); ax.set_ylabel('eccentricity')
plt.colorbar(pt, ax=ax, label='mass ratio')
```

### Python Analysis Guardrails

- Always match `interrupt_mode`, `external_mode`, `snapshot_format`, and `G` to the simulation/data-processing setup.
- For any Python read path, ensure these settings are consistent with the currently selected `petar` binary family (for example `.bse`, `.galpy`, `.agama`, `.dsm`) and the actual producer command used to write the files.
- For original snapshots (`data.*`), use header offset (`petar.HEADER_OFFSET` or `petar.HEADER_OFFSET_WITH_CM` for external-mode snapshots).
- For `petar.data.process` outputs (`*.single`, `*.binary`, `data.lagr`, `data.core`), read directly without header offset.
- For `data.lagr`, use `petar.LagrangianMultiple()`; do not use `petar.Lagrangian()` for functional post-processing outputs.
- When reading outputs from `petar.get.object.snap`, add a `time` member before `fromfile`.
- If column-mismatch or decode warnings appear during Python reading and cannot be resolved by parameter correction, stop the workflow and report the issue instead of continuing.
- If user asks for “give me Python code for this analysis”, return directly runnable snippet(s) using these templates and the user’s path/filename choices.

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
petar.data.process -t agama --r-escape tidal -G 0.00449830997959438 data.snap.lst
```

### BSE + Agama

Default post-processing:

```bash
petar.data.gether data
petar.data.process -i bse -t agama --r-escape tidal -G 0.00449830997959438 data.snap.lst
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
petar.data.process -t agama --r-escape tidal -G 0.00449830997959438 data.snap.lst
```

Compatibility rule:

- any option starting with `--agama-` requires an `.agama` executable.

If the input is already a PeTar snapshot with external-potential context prepared, skip `petar.init`.

### 4) Restart / Resume

```bash
# Keep old parameters from <output_prefix>.par, override selected options after -p
petar -p data.par -t <new_end_time> <snapshot_file>
```

To overwrite instead of append on restart, include `-a 0`.

For restart workflows, do not emit `petar.init`.

If restarting from a non-final snapshot, recommend `petar.data.clear` before the restart command when duplicate events are a concern.

If restart or resume is requested, do not propose `petar.init` as a fallback.

## Parallel And Performance Safety Notes

Always include these reminders when applicable:

- OpenMP: set `OMP_STACKSIZE=128M`; ensure `ulimit -s unlimited`.
- MPI+OpenMP launch pattern:
  `OMP_STACKSIZE=128M OMP_NUM_THREADS=<threads> mpiexec -n <n_mpi> petar [options] <snapshot>`
- If MPI binds one core per rank, try `mpiexec --bind-to none`.
- GPU: one MPI rank typically drives one GPU job; tune ranks/threads to avoid CUDA OOM.
- `*.hard.debug` binaries are helper tools for hard-dump diagnostics (`[output_prefix].hard_dump.*` replay), not full simulation drivers.
- `*.format.transfer` binaries are format-conversion helpers.
- Neither `*.hard.debug` nor `*.format.transfer` should be used as the main simulation executable.
- For source-level debugging, reconfigure with `--with-debug=g` and rebuild before launching debugger workflows.

### Parallel Sizing Heuristic (N and Primordial Binaries)

This heuristic must stay consistent with the comments in `sample/star_cluster_plummer_N1k*.sh`.

Use this default sizing rule before emitting launch commands:

1. Infer workload scale from initial particle number `N` and primordial-binary richness.
2. Prefer under-subscription over over-subscription for small-`N` runs.

Default recommendations:

- `N ~ 10^3` and no primordial binaries (or very few):
  use single-process single-thread by default (`mpiexec -n 1`, `OMP_NUM_THREADS=1`).
- `N ~ 10^4`:
  multi-core OpenMP and/or MPI can be beneficial; start from moderate parallelism, then tune.
- `N ~ 10^3` with many primordial binaries:
  OpenMP multi-threading can be beneficial; consider `OMP_NUM_THREADS > 1`.

When user does not provide this information, ask only these missing items:

- approximate particle count scale (`10^3`, `10^4`, larger), and
- whether primordial binaries are none/few/many.

Command emission rule:

- Do not leave OpenMP thread count implicit.
- Explicitly set `OMP_NUM_THREADS` (and MPI rank count when used) in every launch command.
- For small isolated runs, default to:
  `OMP_STACKSIZE=128M OMP_NUM_THREADS=1 petar [options] <snapshot>`
  or
  `OMP_STACKSIZE=128M OMP_NUM_THREADS=1 mpiexec -n 1 petar [options] <snapshot>`.

## Response Format When Skill Is Invoked

Return output in this fixed structure:

1. Chosen scenario and assumptions.
2. Selected binary, inferred feature stack, and inferred input-source workflow.
3. Exact command block(s), ready to copy.
4. Why each critical option is used (`-u`, `-b`, `-t`, `-o`, `-p`, `--galpy-set`, `--agama-conf-file`, `--bse-metallicity`, etc.).
5. Tool-driven follow-up commands when relevant (`petar.find.dt`, `petar.data.gether`, `petar.data.clear`, `petar.movie`, etc.).
6. Validation checklist (expected output files and quick sanity checks).
7. Optional next-step optimization knobs (`-s`, `-r`, `--r-search-min`, `--r-bin`).

For interactive chat replies, use this three-stage template before full command details:

1. Requirement detection (physics first):
  - explicitly summarize user-required feature switches:
    `interrupt?` | `external long-timescale? (galpy/agama)` | `external-hard short-timescale? (gasdrag)` | `pn?`.
  - if any is missing and required, ask only those missing items.
2. Default candidate selection (performance second):
  - pick candidates from `Requirement-Driven Solver Filters` first.
  - then rank within matched candidates by `gpu > avx512 > avx2 > omp > mpi`.
3. Special-purpose optional modules (only if user asks):
  - mention `--enable-64b`, `--enable-mpfrc`, `--enable-gperf`, `--with-debug=g|assert` as optional toggles,
    and explain why they are not default choices.

After these three stages, provide the runnable command block(s) and validation checklist.

### Minimal Opening Checklist (Quick Ask)

When user intent is not fully specified, ask this compact checklist first (only missing items):

1. Physics requirements (toggle set):
  - `interrupt`: off | base | bse | mobse | bseEmp | dsm
  - `external (long-timescale tree step)`: off | galpy | agama
  - `external-hard (short-timescale hard integrator)`: off | gasdrag
  - `pn`: off | pnhermite | pnsdar | pnall (or other pn mode)
2. Runtime/performance target: GPU needed? If not, accept best CPU candidate.
3. Core run inputs: unit mode (`-u`), end time (`-t`), output interval (`-o`), initial source (raw | snapshot | restart).
4. Parallel sizing inputs: approximate particle count scale (`10^3` | `10^4` | larger), primordial binaries (none | few | many).
5. Optional special modules (only if user asks): `64b`, `mpfrc`, `gperf`, debug mode (`g` or `assert`).

If all checklist items are already inferable from user input and selected binary name, do not ask again; emit commands directly.

## Guardrails

- Do not invent unsupported PeTar options.
- If user intent conflicts with current build (for example, asks for Galpy but binary lacks it), provide reconfigure command first.
- If `petar.select` has no match, explicitly ask user whether to run configure/install, show the exact commands, and wait for confirmation before executing them.
- Prefer repository sample scripts and README wording over generic N-body advice.
- Never pass through user-provided extra options without checking `<selected_binary> -h` first.
- If the user selects a helper tool binary, explicitly switch to the nearest matching solver binary before generating the run command.
- If the binary suffix already determines the physics stack, do not ask redundant questions about whether to enable that same feature.
- Do not ask the full generic input list when the scenario-specific minimal set is smaller.
- If enough inputs are already present in the user request, generate the command immediately instead of asking confirmation questions.
- Do not emit `petar.init` for existing snapshots or restart workflows.
- Do not omit post-processing when the user asks for a complete end-to-end workflow.
- When the user asks for analysis, conversion, movie generation, or restart cleanup, use the installed script tools instead of describing the task abstractly.
