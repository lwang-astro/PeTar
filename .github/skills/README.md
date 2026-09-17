# Skills Directory Guide

Repository-specific skill definitions, assets, and maintenance notes.

> **Customization-maintenance authority.** This file owns the design goals, layering rules,
> and update checklist for every agent-customization file in this workspace — `AGENTS.md`,
> `.github/agents/*.md`, `.github/skills/**` — in **both** the PeTar and SDAR repositories.
> Referenced from both `AGENTS.md` files. It is **not** part of the simulation read path.

## Design Goals

Every change to a customization file is judged against these four goals:

1. **Context efficiency** — minimize the always-loaded rule text so the read path stays cheap.
2. **Layered disclosure** — top layer holds only must-read rules; scenario-specific detail moves to subset files, so a long context never dilutes attention.
3. **Clear ownership** — `SKILL.md` fully covers *using* PeTar for simulation and data analysis, complete without any agent file or `AGENTS.md`. `AGENTS.md` and the sub-agents cover PeTar *code development* and must not restate what `SKILL.md` covers. SDAR owns SDAR internals; PeTar references and states only the difference.
4. **Single authoritative home** — one home per fact; every other occurrence is a path pointer.

## Simplification Questions (ask before adding or keeping any rule)

1. **Would a mainstream model do this anyway?** If the rule merely describes sensible generic behaviour ("use the right tool", "create the directory if missing", "be concise"), delete it or keep one occurrence. **Exception — never delete counter-intuitive domain rules** even when they look like common sense: halve the `petar.find.dt` result, log to `commands.log`, the pre-solver confirmation gate, `petar.init` argument order, `--r-ratio` and `r_in` moving in the same direction, `mcluster -C 5`, `petar.find.dt -a "-u 1"`.
2. **Is it duplicated or contradictory?** Merge into the owner and replace the other occurrences with pointers. Treat "file A says X, file B says Y" as a defect to fix in the same change — do not leave both.
3. **Is it PeTar+SDAR shared?** Push the shared rule down to SDAR (its `AGENTS.md` is canonical for shared processes such as lessons capture); PeTar references it and states only the PeTar-specific difference.

## Layering Rules

| Fact type | Authoritative home |
|---|---|
| Execution-correctness rules for running simulations / analysing output | `petar-nbody-simulation/SKILL.md` (subsequently loaded every simulation task — keep it to rules and pointers) |
| PeTar development rules (build, configure, `src/`, `tools/`, `test/`) | `AGENTS.md` + `.github/agents/*.md` |
| Reference detail (flag tables, class APIs, per-scenario templates, measured benchmarks) | `.github/skills/…/assets/*.md` (on demand only) |
| Where each former detail went | `SKILL_CONTENT_INDEX.md` (`Layering Contract` + coverage map) |
| Incident records and prevention rules | `assets/lessons-learned.md` (never in the read path) |
| SDAR internals and shared processes | `../SDAR/AGENTS.md` |

Rules about maintaining the customization files themselves belong **here** — never in `SKILL.md` or an agent definition.

## Update Rules

- **Always-loaded files hold rules and pointers only.** `AGENTS.md` must not copy agent-file content, SKILL rules, or HANDOFF material. Maintenance/recovery material (`HANDOFF.md`, `prompt-starters.md`) must never enter the always-loaded read path.
- **Never restate a fact in a second place.** After adding a rule, grep for it; convert other occurrences into `<path> → "<section>"` pointers.
- **State the budget in the change**: report the line/byte delta of `SKILL.md` (and of `AGENTS.md` and `.github/agents/*.md` when touched). Adding a reference-style table or benchmark to `SKILL.md` is a defect — move it to an asset and leave one line behind.
- **Scenario detail goes to a subset file**, not to the top layer.
- **Agent definitions** carry role, tool/scope constraints, and delegation rules only — they never restate `SKILL.md`. They load on every invocation of the agent, so they meet the same context-efficiency bar: rules plus at most one clause of counter-intuitive rationale — persuasion, history, and narrative justification belong in `assets/lessons-learned.md`.
- **Lessons entry titles must be self-contained one-line rules** (`### <date>: <one-line rule>`), so `grep -n '^### ' lessons-learned.md` is a complete index. Do **not** split `lessons-learned.md` into per-lesson files while it stays under the split trigger: the read cost is identical (~index + one entry either way) while splitting adds index drift on every capture. **Split trigger**: file > ~800 lines, one category > ~200 lines, or ambiguous titles.

## Health-Check Procedure

Run this as a standalone review (not while writing new content), then land the findings separately:

1. **Always-loaded budget** — measure the read path (a typical execution path, not just one file). Record it; compare against the previous check.
2. **Redundancy scan** — count occurrences of each hard rule across `AGENTS.md`, agents, `SKILL.md`, and assets. Any rule appearing more than once outside its owner is a defect.
3. **Contradiction scan** — read candidate pairs side by side (e.g. "never assume defaults" vs "two permitted defaults"); fix in the same change.
4. **Generic-behaviour scan** — flag rules that any mainstream model would follow anyway (see Simplification Questions). Keep the counter-intuitive exceptions.
5. **Dead-pointer scan** — every referenced file must exist, and every asset must be reachable from a live pointer.
6. **Ownership check** — confirm `SKILL.md` is self-sufficient for simulation tasks and that no agent file restates it.
7. **Land the fixes**, then update `SKILL_CONTENT_INDEX.md` and append a `Configuration Hygiene` lessons entry.

## Primary Entry Points

- `petar-nbody-simulation/SKILL.md`: compact, execution-focused skill rules. **Authority for running simulations and analysing output**; self-sufficient for simulation tasks.
- `petar-nbody-simulation/SKILL_CONTENT_INDEX.md`: coverage map and the layering contract (where each fact's authoritative home is).
- `petar-nbody-simulation/HANDOFF.md`: cross-session skill-maintenance state.

## Assets (Detailed References)

Read on demand from `SKILL.md`'s "Reference Documents (Must Read)" table:

- `petar-nbody-simulation/assets/minimal-question-sets.md` — per-scenario required-ask lists
- `petar-nbody-simulation/assets/binary-scenario-map.md` — binary-suffix → scenario inference
- `petar-nbody-simulation/assets/script-tools.md` — tool command templates and flag references
- `petar-nbody-simulation/assets/default-postprocessing.md` — per-scenario post-processing and movie defaults
- `petar-nbody-simulation/assets/input-source-workflows.md` — raw / snapshot / restart branching
- `petar-nbody-simulation/assets/data-readback-patterns.md` — Python readback patterns (read before writing analysis code)
- `petar-nbody-simulation/assets/dsm-workflow.md` — DSM scenario detail

Generated or maintenance material (not part of the simulation read path):

- `petar-nbody-simulation/assets/option-matrix.md` (machine-specific binary snapshot)
- `petar-nbody-simulation/assets/option-reference.md` (source-generated option catalog)
- `petar-nbody-simulation/assets/option-list.txt`
- `petar-nbody-simulation/assets/prompt-starters.md` (recovery and regression prompts)
- `petar-nbody-simulation/assets/lessons-learned.md`

## Maintenance Scripts

- `petar-nbody-simulation/assets/update_option_inventory.sh`
- `petar-nbody-simulation/assets/update_script_tool_help.sh`
- `petar-nbody-simulation/assets/check_skill_consistency.sh`
- `petar-nbody-simulation/assets/generate_option_reference.py`

## Verifying a Change

```bash
bash .github/skills/petar-nbody-simulation/assets/check_skill_consistency.sh
```

When compacting or reorganizing `SKILL.md`, update `SKILL_CONTENT_INDEX.md` in the same change so removed detail stays traceable and the layering contract stays accurate.
