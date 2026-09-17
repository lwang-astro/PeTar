# PeTar Agent Guide

This file is the root entry point for agents working in this repository. Keep it short, project-wide, and stable.

## Start Here

1. [README.md](README.md) for project structure, build/test entry points, and user-facing workflow notes.
2. [.github/skills/petar-nbody-simulation/SKILL.md](.github/skills/petar-nbody-simulation/SKILL.md) — **authority for running simulations and analysing output**. Self-sufficient: no agent file needs to be read to complete a simulation task.
3. [test/functional/README.md](test/functional/README.md) for functional smoke design and current defaults.
4. [test/validation/README.md](test/validation/README.md) for validation scenarios and thresholds.

Skill-maintenance notes and session-recovery context live in [HANDOFF.md](.github/skills/petar-nbody-simulation/HANDOFF.md). Read it only when maintaining the skill — not for ordinary simulation or development tasks.

**Before changing any customization file** (`AGENTS.md`, `.github/agents/*.md`, `.github/skills/**`), read [.github/skills/README.md](.github/skills/README.md) — it owns the design goals, layering rules, update rules, and health-check procedure for this workspace.

## Repository Map

- `src/`: core PeTar implementation.
- `tools/`: installed command-line helpers and post-processing scripts.
- `sample/`: runnable examples and reference configurations.
- `test/functional/`: fast workflow smoke coverage.
- `test/validation/`: numerical and physics validation coverage.
- `doc/`: longer-form documentation.

## Scope Split

- **Running simulations / analysing output** → SKILL.md is authoritative and complete on its own.
- **Developing PeTar** (configure, Makefile/Makefile.in, `src/`, `tools/`, `test/`) → this file plus the agent suite below.
- **SDAR internals** (integrator, `SDAR/src/`, standalone few-body runs) → [SDAR/AGENTS.md](../SDAR/AGENTS.md) is authoritative. PeTar only routes to it; see the skill's "SDAR deep-dive routing" note.

## Development Rules

- Prefer repository-relative tracked files as inputs; never rely on untracked local files for committed test logic.
- Treat `test/out/` as generated output only.
- Keep restart/resume workflows separate from fresh IC generation.
- Simulation execution rules (binary selection via `petar.select`, option validation, unit consistency) live in SKILL.md — do not restate or second-guess them here.
- Do not add undocumented configure or runtime options without checking repository docs and `-h` output.
- Do not commit or push without user confirmation.
- Bump `VERSION` in every commit: `(cd tools && bash get_version.sh)`, then append the `e` experiment suffix and stage `VERSION`.

## Custom Agents

Start with **PeTar Developer** for every request. It performs planning and focused edits directly, and routes to **PeTar Simulation Engineer** (all execution) and **PeTar Reviewer** (pre-close review and numerical verdicts). Read-only discovery uses the built-in **Explore**.

Delegation criteria, model allocation, and repository constraints are defined in [`.github/agents/petar-developer.agent.md`](.github/agents/petar-developer.agent.md) — do not restate them here.

## Lessons-Learned Capture

The capture process is canonical in [SDAR/AGENTS.md](../SDAR/AGENTS.md). PeTar routing rule:

- Lessons about PeTar surfaces (configure/build, `petar.*` tools, clustering physics, `test/`) → `.github/skills/petar-nbody-simulation/assets/lessons-learned.md`
- Lessons rooted in SDAR code (`SDAR/src/`, `SDAR/sample/*`, `SDAR/tools/`) → SDAR's lessons file, with a one-line pointer here if PeTar context matters.

## Test Organization

- Use `test/functional` to check the toolchain, binary selection, readback paths, and quick post-processing coverage.
- Use `test/validation` for physics-oriented scenarios and regression thresholds.
- Keep functional smoke defaults documented in [README.md](README.md) and [test/functional/README.md](test/functional/README.md) rather than duplicating them here.
