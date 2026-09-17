# PeTar Minimal Question Sets

This note documents the minimum missing information to collect before generating commands.

## Rule

Ask only for parameters that are both **missing** and **required** for the scenario.
If the request already contains enough information for a runnable command, do not ask follow-up questions.

The mandatory field list lives in `SKILL.md` ("Required Input Checklist" → Common) and must be
complete before anything below applies. Ask these four for **every** scenario, in addition to the
scenario-specific entries:

- unit mode (`-u 0` or `-u 1`)
- end time (`-t`)
- output interval (`-o`)
- parallel launch mode (serial | OpenMP | MPI+OpenMP | GPU, with thread/process counts)

## Scenario extras

| Scenario | Ask in addition |
|---|---|
| **Isolated cluster** | when IC generation is required (no existing snapshot): every row of the "Star-Cluster Generation Inputs" table in `SKILL.md` |
| **BSE / SSE** | primordial binary count (`-b`); metallicity (`--bse-metallicity`); and the IC table above when generating ICs |
| **DSM** | only the DSM controls the user wants to override (e.g. `--dsm-seed-mass`, `--dsm-he-disk`, `--dsm-lambda0`, `--dsm-dt-factor`); matching external-potential inputs if Galpy/Agama is also enabled |
| **Galpy** | cluster COM position and velocity for `petar.init -c`; one of `--galpy-set` / `--galpy-conf-file` / `--galpy-type-arg`; optional `--galpy-rscale`, `--galpy-vscale`; plus `-b` and metallicity if BSE is also enabled |
| **Agama** | cluster COM position and velocity for `petar.init -c`; `--agama-conf-file`; optional `--agama-rscale`, `--agama-vscale`; plus `-b` and metallicity if BSE is also enabled |
| **Restart** | restart snapshot filename; parameter file after `-p`; options being overridden; append/overwrite mode (`-a`); companion-file and `-i` requirements are owned by `input-source-workflows.md` → "Restart / resume" |
| **standalone-bse** | no IC generation and no `petar.init` — confirm the target binary (`petar.bse` / `petar.mobse` / `petar.bseEmp`) and its metallicity |
