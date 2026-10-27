# PeTar Agent Guide

This file is the root entry point for agents working in this repository. Keep it short, project-wide, and stable.

## Start Here

1. [README.md](README.md) for project structure, build/test entry points, and user-facing workflow notes.
2. [.github/skills/petar-nbody-simulation/SKILL.md](.github/skills/petar-nbody-simulation/SKILL.md) — **authority for running simulations and analysing output**. Self-sufficient: no agent file needs to be read to complete a simulation task.
3. [test/functional/README.md](test/functional/README.md) for functional smoke design and current defaults.
4. [test/validation/README.md](test/validation/README.md) for validation scenarios and thresholds.

**Before changing any customization file** (`AGENTS.md`, `.github/agents/*.md`, `.github/skills/**`), read [.github/skills/README.md](.github/skills/README.md) — it owns the design goals, layering rules, update rules, and health-check procedure for this workspace.

## Repository Map

- `src/`: core PeTar implementation.
- `tools/`: installed command-line helpers and post-processing scripts.
- `sample/`: runnable examples and reference configurations.
- `test/functional/`: fast workflow smoke coverage.
- `test/validation/`: numerical and physics validation coverage.
- `doc/`: longer-form documentation.

## Documentation

- [doc/output-format-contract_notes.md](doc/output-format-contract_notes.md) — reference note: output-format contract (bulk BINARY vs event-stream ASCII), its reader/test enforcement (2026-10-19), and the shelved P2-lite/P3 decisions with revival preconditions.

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
- Never commit or push autonomously: every commit needs explicit user confirmation for that specific change first — reporting finished work is not authorization.
- Bump `VERSION` in every commit: `(cd tools && bash get_version.sh)` and stage
  `VERSION` — the script preserves the `e` experiment suffix automatically.
- User-facing changes update `release.note` in the same commit:
  - one entry = one change = one line <= 150 chars, starting with a standard verb: **Add / Fix / Change / Improve / Remove**;
  - prefix the introducing version (`66. (1790e) Fix ...`); major entries get a trailing `[major]` tag (never a leading `*` — it breaks `grep -E '^[0-9]+\. '`); carry SDAR requirements over from the commit message (`needs SDAR >= 63e8308`);
  - what changed and why it matters only — root causes, fix mechanics, and validation numbers go to the commit message / lessons; minor fixes, internal work (lessons, agent infra, doc hygiene), and resolved issues (one-clause strikethrough) stay out;
  - draft: `bash tools/gen_release_draft.sh <last-release.note-commit>..HEAD`; insert and verify with `python3 tools/release_note_tool.py insert|sort|check` (hardened sort key for bare `(V)` prefixes; `check` enforces numbering, sort order, verbs, version prefixes, and the length cap for entries >= 1848e).

## Custom Agents

Start with **PeTar Developer** for every request. It performs planning and focused edits directly, and routes to **PeTar Simulation Engineer** (all execution) and **PeTar Reviewer** (pre-close review and numerical verdicts). Read-only discovery uses the built-in **Explore**.

Delegation criteria, model allocation, and repository constraints are defined in [`.github/agents/petar-developer.agent.md`](.github/agents/petar-developer.agent.md) — do not restate them here.

## Lessons-Learned Capture

The capture process is canonical in [SDAR/AGENTS.md](../SDAR/AGENTS.md). PeTar routing rule:

- Lessons about PeTar surfaces (configure/build, `petar.*` tools, clustering physics, `test/`) → `.github/skills/petar-nbody-simulation/assets/lessons-learned.md`
- Lessons rooted in SDAR code (`SDAR/src/`, `SDAR/sample/*`, `SDAR/tools/`) → SDAR's lessons file, with a one-line pointer here if PeTar context matters.

## Doc & Plan Organization

- Every md in `doc/` belongs to exactly one class — unclassified notes are a defect:
  - **plan** (`plan` in the filename): follow-up work; open items only (see below).
  - **reference note** (e.g. `*_notes.md`): permanent documentation of current-code mechanism (implementation notes), theory, or removed legacy behavior; carries a one-line status header stating its nature, and is registered in the repo's documentation listing (SDAR: `AGENTS.md` → Documentation). Retained limitations/caveats of implemented schemes live here, not in plans.
  - **temporary** (debug/handoff notes): no permanent home — distill then delete (below).
  Reference material (papers, method PDFs) needs no classification.
- When a plan item completes, distill its experience into the owning note in the same change; the plan keeps only open items plus pointers to where the completed content went. Closed-item records never live in the plan (no "Closed items" section) — at most a one-line pointer per closed item.
- Routing rule — **lessons-learned is the default home**. Lessons take anything with the Mistake / Root cause / Prevention shape: "what went wrong, why, how to avoid" (operational mistakes, methodology gotchas, diagnostic fingerprints); hard rules get promoted into SKILL assets. A reference note takes **only** documentation of how implemented code works or worked — mechanism, algebra, invariants, measured behavior, design rationale, retained limitations — answering "how/why does the code behave this way", never "what went wrong". An investigation yielding both splits: prevention rules → lessons entry, mechanism description → reference note, cross-linked. Boundary unclear → lessons.
- One-off debug/handoff notes are temporary: distill the necessary experience into lessons first, then delete — tracked files keep git history, but untracked content is lost on `rm`, so lessons capture comes first.
- Moving/renaming/deleting any doc updates its cross-references in the same change (`grep -rn "<name>"` across the workspace, including SDAR docs).
- Keep `doc/` free of temp/backup artifacts (`*~`, `Zone.Identifier`, scratch dumps).

## Test Organization

- Use `test/functional` to check the toolchain, binary selection, readback paths, and quick post-processing coverage.
- Use `test/validation` for physics-oriented scenarios and regression thresholds.
- Keep functional smoke defaults documented in [README.md](README.md) and [test/functional/README.md](test/functional/README.md) rather than duplicating them here.
