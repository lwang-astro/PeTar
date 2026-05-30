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
3. use restart snapshot as final positional argument
4. optionally use `-a 0` to overwrite outputs

Canonical form:
`<launcher> <petar_binary> -p <output_prefix>.par [overrides] <restart_snapshot>` (default: `data.par`)

## Repository sample workflow

Definition:
- user wants something matching repository sample scripts

Workflow:
1. follow repository sample ordering
2. replace placeholders only where needed
3. keep repository defaults unless the user overrides them
