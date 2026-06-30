# Lessons Learned — PeTar Agent Workflow

This file records mistakes, gotchas, and pitfalls discovered during PeTar agent sessions.
Each entry documents: what went wrong, why, and how to prevent recurrence.

**Lifecycle**: New entries are added by agents automatically after encountering issues.
Periodically reviewed → verified entries are elevated to `SKILL.md` as hard rules.

---

## Categories

- [Build & Configure](#build--configure)
- [Binary Selection](#binary-selection)
- [Simulation Execution](#simulation-execution)
- [Post-Processing](#post-processing)
- [Restart / Resume](#restart--resume)
- [Documentation](#documentation)

---

## Build & Configure

*(No entries yet)*

---

## Binary Selection

*(No entries yet)*

---

## Simulation Execution

### 2026-06-30: Silent defaults for physics-defining IC parameters
**Mistake**: When user requested a Plummer N=500 cluster without specifying half-mass radius, virial ratio, IMF, mass range, or tidal filling, the agent silently used mcluster defaults (Rh=0.8 pc, Q=0.5, Kroupa IMF, 0.08–150 Msun) without asking — violating reproducibility and potentially producing unintended physics.
**Root cause**: The SKILL.md "Star-Cluster Generation Inputs" checklist was too brief (5 items) and did not enumerate all parameters that must be confirmed. The agent interpreted "do not proceed until the IC parameter set is complete" as satisfied by the 5 listed items, ignoring unlisted parameters.
**Prevention rule**: The definitive checklist in SKILL.md must enumerate every mcluster parameter that affects the initial physical state. Agents must iterate through all applicable rows and ask for any missing value. A "not safe to default" classification must be applied to Rh, Q, IMF, mass range, and COM phase-space.

---

## Post-Processing

*(No entries yet)*

---

## Restart / Resume

*(No entries yet)*

---

## Documentation

*(No entries yet)*
