# SKILL Content Index (Compact Version Coverage)

This file tracks where detailed content from earlier long-form skill versions is documented after the skill was compacted.

## Purpose

- Keep `SKILL.md` compact and execution-focused.
- Preserve traceability of previously documented details.
- Point agents and maintainers to the authoritative location for each topic.

## Coverage Map

### 0. Mandatory execution guardrails and hard constraints

- Primary authority: `.github/skills/petar-nbody-simulation/SKILL.md`
- Supporting examples and edge-case detail: repository sample scripts, `README.md`, and skill assets
- In-skill coverage now includes:
  - selection-before-run and no-implicit-`petar` execution
  - required-input completion before command generation
  - confirmation-before-execution
  - unit-consistency guardrails
  - IC-generation completion requirements
  - restart/resume and post-processing blocking rules

Status: restored-in-skill

### 1. Binary selection, capability checks, and option compatibility

- Primary rule (compact): `.github/skills/petar-nbody-simulation/SKILL.md`
- Detailed capability inventory: `.github/skills/petar-nbody-simulation/assets/option-matrix.md`
- Installed tool surfaces and help-derived descriptions: `.github/skills/petar-nbody-simulation/assets/script-tools.md`
- In-skill coverage now includes:
  - physics-sensitive tokens stay in `--require`
  - performance tokens stay in `--optional`
  - selected binary must be checked with `-h`
  - helper binaries are not valid production solvers

Status: covered, with core hard constraints kept in `SKILL.md`

### 2. Minimal question sets and input collection details

- Compact ask-minimum principle: `.github/skills/petar-nbody-simulation/SKILL.md`
- Detailed question templates: `.github/skills/petar-nbody-simulation/assets/minimal-question-sets.md`
- In-skill coverage now includes:
  - scenario-required inputs must be complete before run commands
  - BSE metallicity requirement
  - external-potential model and COM phase-space requirement
  - star-cluster IC parameter completeness requirement

Status: covered, with required-input guardrails kept in `SKILL.md`

### 3. Input-source workflow branching (raw/snapshot/restart)

- Compact workflow patterns: `.github/skills/petar-nbody-simulation/SKILL.md`
- Detailed branching reference: `.github/skills/petar-nbody-simulation/assets/input-source-workflows.md`
- In-skill coverage now includes:
  - raw IC unit confirmation before `petar.init`
  - normalized-unit recommended path
  - native-unit path treated as advanced and opt-in
  - restart workflows must not emit `petar.init`

Status: covered, with restart and `petar.init` guardrails kept in `SKILL.md`

### 4. Default post-processing flows by scenario

- Compact post-processing policy: `.github/skills/petar-nbody-simulation/SKILL.md`
- Detailed defaults: `.github/skills/petar-nbody-simulation/assets/default-postprocessing.md`
- In-skill coverage now includes:
  - use installed PeTar tools instead of only abstract advice
  - mode matching between producer family and reader/post-processing tool
  - MPI gather requirement before downstream tools when applicable

Status: covered, with mandatory tool-use and mode-matching rules kept in `SKILL.md`

### 4a. Installed tool command patterns and extended tool behaviors

- Compact command-availability rule: `.github/skills/petar-nbody-simulation/SKILL.md`
- Detailed tool behaviors and inventory: `.github/skills/petar-nbody-simulation/assets/script-tools.md`
- Sample workflow references: `sample/` and `README.md`

Status: split between `SKILL.md` guardrails and asset-level detail

### 4b. Post-processing tool command templates

- Compact tool-usage rule (use installed tools, mode matching): `.github/skills/petar-nbody-simulation/SKILL.md`
- Detailed command templates for `petar.get.object.snap` and `petar.format.transfer.post`: `.github/skills/petar-nbody-simulation/assets/script-tools.md`
- Per-scenario `petar.movie` default arguments: `.github/skills/petar-nbody-simulation/assets/default-postprocessing.md`

Status: covered, with tool-specific templates and movie defaults in asset files

### 4c. Python data readback patterns

- Compact reading guidance and keyword-argument table: `.github/skills/petar-nbody-simulation/SKILL.md` (Python Data Analysis Tools section)
- Detailed readback patterns (10 patterns: lagr, core, status, escaper, SSE/BSE, snapshot offsets, object snap, profile matching, GroupInfo, DSM interrupt): `.github/skills/petar-nbody-simulation/assets/data-readback-patterns.md`
- Primary tutorial reference: `sample/data_analysis.ipynb`

Status: covered, with SKILL.md providing the compact keyword-argument reference and asset providing full readback strategies

### 4d. DSM workflow reference

- Compact DSM rules (init, runtime, post-processing): `.github/skills/petar-nbody-simulation/SKILL.md` (DSM-specific section)
- Detailed DSM workflow (IC prep, radius calculation, parameter reasoning, artifact checklist): `.github/skills/petar-nbody-simulation/assets/dsm-workflow.md`

Status: covered, with hard rules in SKILL.md and detailed reference in asset

### 5. Functional smoke defaults and test-specific policies

- User-facing summary and current defaults: `README.md` (Automated Test Layers -> Functional smoke tests)
- Functional runner details and matrix scope: `test/functional/README.md`

Status: covered

### 6. Validation scenarios and thresholds

- User-facing validation entry: `README.md` (Automated Test Layers -> Validation scenarios)
- Detailed validation design and criteria: `test/validation/README.md`

Status: covered

### 7. Prompt templates and session recovery shortcuts

- Recovery and maintenance notes: `.github/skills/petar-nbody-simulation/HANDOFF.md`
- Prompt templates: `.github/skills/petar-nbody-simulation/assets/prompt-starters.md`

Status: covered

## What Stays Out Of SKILL.md

The compact skill should not re-expand into a full host-inventory or long procedural manual. These details remain outside the main skill unless they become mandatory safety rules:

- exhaustive installed-tool descriptions
- full binary inventory snapshots
- long prompt starter catalogs
- long-form workflow examples already covered by `sample/` and `README.md`
- host-specific continuation notes

## Notes For Maintainers

- When removing detail from `SKILL.md`, update this index in the same change.
- Do not move mandatory run-safety or execution-blocking constraints out of `SKILL.md` unless an equivalent hard-rule location replaces them.
- Keep this file as a stable map; do not place transient run results here.
- If a topic cannot be mapped to another file, either:
  - restore a concise version into `SKILL.md`, or
  - add a new dedicated doc and reference it here.
