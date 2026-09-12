# PeTar Agent Guide

This file is the root entry point for agents working in this repository. Keep it short, project-wide, and stable.

## Start Here

1. [README.md](README.md) for project structure, build/test entry points, and user-facing workflow notes.
2. [.github/skills/petar-nbody-simulation/SKILL.md](.github/skills/petar-nbody-simulation/SKILL.md) for the executable PeTar workflow rules used by the VS Code skill system.
3. [.github/skills/petar-nbody-simulation/HANDOFF.md](.github/skills/petar-nbody-simulation/HANDOFF.md) for continuation notes and longer-lived skill context.
4. [test/functional/README.md](test/functional/README.md) for functional smoke design and current defaults.
5. [test/validation/README.md](test/validation/README.md) for validation scenarios and thresholds.

## Repository Map

- `src/`: core PeTar implementation.
- `tools/`: installed command-line helpers and post-processing scripts.
- `sample/`: runnable examples and reference configurations.
- `test/functional/`: fast workflow smoke coverage.
- `test/validation/`: numerical and physics validation coverage.
- `doc/`: longer-form documentation.

## General Working Rules

- Prefer repository-relative tracked files as inputs.
- Do not depend on untracked local files for committed test logic.
- Treat `test/out/` as generated output only.
- Use `petar.select` before generating any run command that depends on a specific binary family.
- Validate custom options against the selected binary help output before using them.
- Keep restart/resume workflows separate from fresh IC generation.
- When a task is clearly a simulation or post-processing task, consult the skill and the relevant README before acting.

## Test Organization

- Use `test/functional` to check the toolchain, binary selection, readback paths, and quick post-processing coverage.
- Use `test/validation` for physics-oriented scenarios and regression thresholds.
- Keep functional smoke defaults documented in [README.md](README.md) and [test/functional/README.md](test/functional/README.md) rather than duplicating them here.

## Custom Agents

The suite was consolidated on 2026-09-12 (see `.github/skills/petar-nbody-simulation/assets/lessons-learned.md`, "Agent Workflow & Delegation"): delegation is for context isolation and long-output isolation, not for packaging already-derived context.

- [PeTar Developer](.github/agents/petar-developer.agent.md): Conductor. Orchestrates work, writes plans and makes focused code/doc edits directly, and delegates execution and review below.
- [PeTar Simulation Engineer](.github/agents/petar-simulation-engineer.agent.md): All execution - configure/build/install maintenance, binary-family selection, simulation runs, functional smoke and validation harness execution, post-processing, restart debugging (`assist`/`debug` modes). Absorbs the former Build and Test Maintainer.
- [PeTar Reviewer](.github/agents/petar-reviewer.agent.md): Pre-close review of non-trivial changes plus numerical validation analysis (validation-layer choice, T1-T3, threshold verdicts). Absorbs the former Validation Analyst.

Read-only discovery uses the built-in **Explore** agent. Planning, focused implementation, and documentation sync are performed directly by PeTar Developer.

## Delegation Pattern

- Start with [PeTar Developer](.github/agents/petar-developer.agent.md) for all requests; it works directly when the context is already established: planning, focused code/script edits, documentation sync, and lessons-learned entries.
- Use **Explore** (built-in) when owning code or workflow surfaces are unknown and read-heavy discovery is needed.
- Use [PeTar Simulation Engineer](.github/agents/petar-simulation-engineer.agent.md) for all simulation execution tasks — end-user workflow assistance, build/install checks, binary selection, command composition, functional smoke, post-processing, restart debugging, and runtime troubleshooting. It operates in `assist` mode (conversational, confirmation-driven) or `debug` mode (technical, execution-heavy) as appropriate.
- Use [PeTar Reviewer](.github/agents/petar-reviewer.agent.md) before closing non-trivial changes, and whenever a verdict on numerical behavior or validation thresholds is required.
- Delegate only for new-context acquisition or long-output isolation; never to repackage already-derived context.

## If You Need More Detail

- Scenario and command templates: [README.md](README.md)
- Workflow rules and option constraints: [.github/skills/petar-nbody-simulation/SKILL.md](.github/skills/petar-nbody-simulation/SKILL.md)
- Cross-session context: [.github/skills/petar-nbody-simulation/HANDOFF.md](.github/skills/petar-nbody-simulation/HANDOFF.md)
- SDAR close-encounter / few-body subsystem work: consult the SDAR repository's [AGENTS.md](../SDAR/AGENTS.md) and [sdar-fewbody-integration skill](../SDAR/.github/skills/sdar-fewbody-integration/SKILL.md) before replicating a subsystem in standalone SDAR, debugging SDAR group detection, or modifying `SDAR/src/`.
