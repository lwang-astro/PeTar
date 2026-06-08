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

- [PeTar Developer](.github/agents/petar-developer.agent.md): Conductor agent. Use for multi-step PeTar work or broad PeTar requests that should be routed across research, implementation, simulation, validation, and review specialists.
- [PeTar Planner](.github/agents/petar-planner.agent.md): Planning-first agent for writing phased implementation plans that can be handed to PeTar Developer.
- [PeTar Researcher](.github/agents/petar-researcher.agent.md): Read-heavy investigator for locating owning files, symbols, workflow controls, and required documentation.
- [PeTar Implementer](.github/agents/petar-implementer.agent.md): Focused editor for C++, Python, build, test, and documentation changes.
- [PeTar Build and Test Maintainer](.github/agents/petar-build-test-maintainer.agent.md): Configure/build, binary-family selection, smoke harness, and validation entry-point specialist.
- [PeTar Simulation Engineer](.github/agents/petar-simulation-engineer.agent.md): Unified simulation execution specialist. Handles end-user simulation setup (assist mode) and technical workflow/build/post-processing debugging (debug mode). Replaces Operator and Specialist.
- [PeTar Validation Analyst](.github/agents/petar-validation-analyst.agent.md): Numerical regression specialist for T1-T3 and scenario-based validation.
- [PeTar Reviewer](.github/agents/petar-reviewer.agent.md): Review specialist for changed-slice correctness, workflow regressions, and documentation/test gaps.
- [PeTar Documentation Maintainer](.github/agents/petar-documentation-maintainer.agent.md): README, SKILL, sample-script, and agent-guide synchronization specialist.

## Delegation Pattern

- Start with [PeTar Developer](.github/agents/petar-developer.agent.md) for broad or ambiguous requests.
- Use [PeTar Planner](.github/agents/petar-planner.agent.md) when the user wants a plan first or when a large task should be broken into explicit phases before implementation.
- Use [PeTar Researcher](.github/agents/petar-researcher.agent.md) when you first need code ownership or workflow mapping.
- Use [PeTar Implementer](.github/agents/petar-implementer.agent.md) for focused edits after the target surface is known.
- Use [PeTar Build and Test Maintainer](.github/agents/petar-build-test-maintainer.agent.md) for `configure`/`make` issues and for choosing or repairing the right automated check.
- Use [PeTar Simulation Engineer](.github/agents/petar-simulation-engineer.agent.md) for all simulation execution tasks — end-user workflow assistance, build/install checks, binary selection, command composition, functional smoke, post-processing, restart debugging, and runtime troubleshooting. It operates in `assist` mode (conversational, confirmation-driven) or `debug` mode (technical, execution-heavy) as appropriate.
- Use [PeTar Validation Analyst](.github/agents/petar-validation-analyst.agent.md) only when the task depends on numerical behavior or validation thresholds.
- Use [PeTar Documentation Maintainer](.github/agents/petar-documentation-maintainer.agent.md) when user-facing workflow guidance must change with the implementation.
- Use [PeTar Reviewer](.github/agents/petar-reviewer.agent.md) before closing non-trivial changes.

## If You Need More Detail

- Scenario and command templates: [README.md](README.md)
- Workflow rules and option constraints: [.github/skills/petar-nbody-simulation/SKILL.md](.github/skills/petar-nbody-simulation/SKILL.md)
- Cross-session context: [.github/skills/petar-nbody-simulation/HANDOFF.md](.github/skills/petar-nbody-simulation/HANDOFF.md)
