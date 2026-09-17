# PeTar Input Source Workflows

This note documents how the skill should switch workflow based on the user's input source.

## Raw particle table

Definition:
- plain particle data with mass, position, velocity columns
- or a generator command that produces such a file

Workflow:
1. generate or identify raw input
2. run `petar.init`
3. run selected solver binary
4. run scenario-specific post-processing

Common additions:
- use `-v kms2pcmyr` for astrophysical example workflows in this repository
- add `-s bse` when stellar evolution is enabled
- add `-c ... -t` when external potential context is required at initialization

## Existing PeTar snapshot

Definition:
- user already has an `input`, `data.*`, or other valid PeTar snapshot

Workflow:
1. skip `petar.init`
2. run selected solver binary directly
3. add post-processing if a full workflow is requested

## Restart / resume

Definition:
- user provides a restart snapshot and parameter file
- or explicitly asks to continue a previous run

Workflow:
1. use `-p <parameter_file>`
2. place overrides after `-p`
3. copy **all** `<prefix>.par.*` companion files (not just `<prefix>.par`) into the restart directory
4. match `-i` to the snapshot being **read**: default `-i 2` reads ASCII and writes binary, so a binary restart snapshot needs `-i 0` (binary read + write) or `-i 3` (binary read, ASCII write)
5. use restart snapshot as final positional argument
6. optionally use `-a 0` to overwrite outputs

Canonical form:
`<launcher> <petar_binary> -p <output_prefix>.par -i 0 [overrides] <restart_snapshot>` (default prefix: `data`)

Worked example:
```bash
# data.par, data.par.hard (and any other data.par.*) plus the binary snapshot data.60 present
petar -p data.par -i 0 -t 16.0 -o 0.25 -s 0.001953125 data.60
```

Failure modes when the above is skipped:

| Omission | Error message | Cause |
|----------|---------------|-------|
| Missing `<prefix>.par.*` companion | `Cannot open file <prefix>.par.<feature>` at startup | Feature parameters live in the companions, not in `<prefix>.par` |
| Missing `-i 0` on a binary snapshot | `cannot read header` (core dump under MPI) | Default `-i 2` **reads** ASCII; binary snapshots are the default write mode, so the read half must be switched to binary (`-i 0` or `-i 3`) |

Known-good reference: restarting an N=500 run from `data.60` (t=15 Myr) with the full companion set and `-i 0` reproduced the original run's energies **bitwise** (ΔE = 0.000e+00 over the restarted interval).

## Repository sample workflow

Definition:
- user wants something matching repository sample scripts

Workflow:
1. follow repository sample ordering
2. replace placeholders only where needed
3. keep repository defaults unless the user overrides them
