# Controlled Experiment Checklist

Read before running a parameter scan or a cross-build comparison. Pointed to from
`SKILL.md` → "Reference Documents (Must Read)".

## General checks

- [ ] Parameters under study are explicitly set (not relying on defaults); verify derived values
      match intent
- [ ] A control run with known-good parameters (e.g., defaults from sample scripts) behaves as
      expected before scanning

Version consistency (Gate 4), unit consistency (Gate 1), and binary-family matching are enforced by
the `SKILL.md` gates — re-verify them there rather than duplicating here.

## Scenario-specific checks

Refer to `SKILL.md` → "Required Input Checklist" for scenario-specific mandatory parameters, and to
`assets/changeover-tuning.md` → "Environment considerations" for changeover/timestep strategy per
system type.
