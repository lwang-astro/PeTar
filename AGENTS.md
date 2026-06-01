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

- [PeTar Developer](.github/agents/petar-developer.agent.md): Full-stack PeTar development agent. Use for C++ core development, build system changes, physics modules, Python analysis tools, test creation, interface integration, and documentation maintenance.

## If You Need More Detail

- Scenario and command templates: [README.md](README.md)
- Workflow rules and option constraints: [.github/skills/petar-nbody-simulation/SKILL.md](.github/skills/petar-nbody-simulation/SKILL.md)
- Cross-session context: [.github/skills/petar-nbody-simulation/HANDOFF.md](.github/skills/petar-nbody-simulation/HANDOFF.md)
