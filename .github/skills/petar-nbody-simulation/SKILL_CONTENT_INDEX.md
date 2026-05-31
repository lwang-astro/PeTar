# SKILL Content Index (Compact Version Coverage)

This file tracks where detailed content from earlier long-form skill versions is documented after the skill was compacted.

## Purpose

- Keep `SKILL.md` compact and execution-focused.
- Preserve traceability of previously documented details.
- Point agents and maintainers to the authoritative location for each topic.

## Coverage Map

### 1. Binary selection, capability checks, and option compatibility

- Primary rule (compact): `.github/skills/petar-nbody-simulation/SKILL.md`
- Detailed capability inventory: `.github/skills/petar-nbody-simulation/assets/option-matrix.md`
- Installed tool surfaces and help-derived descriptions: `.github/skills/petar-nbody-simulation/assets/script-tools.md`

Status: covered

### 2. Minimal question sets and input collection details

- Compact ask-minimum principle: `.github/skills/petar-nbody-simulation/SKILL.md`
- Detailed question templates: `.github/skills/petar-nbody-simulation/assets/minimal-question-sets.md`

Status: covered

### 3. Input-source workflow branching (raw/snapshot/restart)

- Compact workflow patterns: `.github/skills/petar-nbody-simulation/SKILL.md`
- Detailed branching reference: `.github/skills/petar-nbody-simulation/assets/input-source-workflows.md`

Status: covered

### 4. Default post-processing flows by scenario

- Compact post-processing policy: `.github/skills/petar-nbody-simulation/SKILL.md`
- Detailed defaults: `.github/skills/petar-nbody-simulation/assets/default-postprocessing.md`

Status: covered

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

## Notes For Maintainers

- When removing detail from `SKILL.md`, update this index in the same change.
- Keep this file as a stable map; do not place transient run results here.
- If a topic cannot be mapped to another file, either:
  - restore a concise version into `SKILL.md`, or
  - add a new dedicated doc and reference it here.
