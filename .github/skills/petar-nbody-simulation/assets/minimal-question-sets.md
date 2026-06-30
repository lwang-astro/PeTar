# PeTar Minimal Question Sets

This note documents the minimum missing information that the skill should ask for before generating commands.

## Isolated cluster

When IC generation is required (no existing snapshot), ask for ALL fields in the "Star-Cluster Generation Inputs" checklist (SKILL.md) before composing commands.

Runtime-only asks (when snapshot already exists):
- unit mode (`-u 0` or `-u 1`)
- end time (`-t`)
- output interval (`-o`)
- parallel launch mode

## BSE / SSE

When IC generation is required (no existing snapshot), ask for ALL fields in the "Star-Cluster Generation Inputs" checklist (SKILL.md) before composing commands. Pay special attention to:
- primordial binary count (`-b`)
- parallel launch mode

Runtime-only asks (when snapshot already exists):
- unit mode
- end time (`-t`)
- output interval (`-o`)
- primordial binary count (`-b`)
- metallicity (`--bse-metallicity`)
- parallel launch mode

## DSM

Ask only if missing:
- initial condition source
- unit mode
- end time (`-t`)
- output interval (`-o`)
- DSM key controls to override (for example `--dsm-seed-mass`, `--dsm-he-disk`, `--dsm-lambda0`, `--dsm-dt-factor`)
- parallel launch mode

If `.galpy` or `.agama` is also enabled, also ask for matching external-potential inputs.

## Galpy

Ask only if missing:
- initial condition source
- unit mode
- end time (`-t`)
- output interval (`-o`)
- cluster center coordinates and velocity for `petar.init -c`
- one of: `--galpy-set`, `--galpy-conf-file`, `--galpy-type-arg`
- optional scaling: `--galpy-rscale`, `--galpy-vscale`
- parallel launch mode

If `.bse` is also enabled, also ask for:
- primordial binary count
- metallicity

## Agama

Ask only if missing:
- initial condition source
- unit mode
- end time (`-t`)
- output interval (`-o`)
- cluster center coordinates and velocity for `petar.init -c`
- `--agama-conf-file`
- optional scaling: `--agama-rscale`, `--agama-vscale`
- parallel launch mode

If `.bse` is also enabled, also ask for:
- primordial binary count
- metallicity

## Restart

Ask only if missing:
- restart snapshot filename
- parameter file path after `-p`
- overridden options
- append or overwrite mode (`-a`)
- parallel launch mode

## Rule

If the request already contains enough information for a runnable command, do not ask follow-up questions.
