# Output schema tagging and auto-detection — plan

Status: plan (open items only) — S1–S5 executed 2026-10-30; mechanism
documentation distilled to
`.github/skills/petar-nbody-simulation/assets/data-readback-patterns.md` →
"Schema auto-detection" (precedence, uniqueness table, tier policy,
sidecar cross-validation rules, caveats). Implementation pointers:
`tools/analysis/schema.py` (detection + sidecar + `getDataKwargs`/
`readSnapshot`/`autoFillKwargs`), `src/io.hpp` `writeOutputMeta()` called
from the parameter-initialization block in `src/petar.hpp` (refreshed on
restart), CLI auto-fill in
`tools/{data_process,get_object_snapshot,format_transfer,movie}.py`,
per-case smoke checks `schema:meta|detect|autofill` plus the one-time
`run_schema_selftest()` in `test/functional/run_functional_smoke.py`.
Decisions taken at execution: O3 → tiered policy (rarest-flag count;
lowest tier wins with a warning, same-tier collision raises); O4 →
resolved (all producing tools write inside the run directory); O5 →
resolved without configure changes (BSE variants via `BSEEMP`/`MOBSE`
macros, version via `GetVersion()`).

## Open items

- **O1** MPFR/float128 detection: the binary path carries the 32/128-byte
  header variants and couples `use_mpfrc` with float128 records; a
  hypothetical float128 build *without* `pos_high` columns is not in the
  candidate space and raises no-match (explicit kwargs or the sidecar
  cover it). Verify once such a build is actually exercised. Related
  residual: binary float64/float128 header aliasing at contrived N
  (documented as a caveat in the readback asset).
- **O2** Extend sniffing to `data.esc_*`, `data.group.n*`, `data.lagr`,
  `object.*` (weaker signals: headerless, N unknown) or keep them
  explicit-kwargs + meta-only (current state).
- **O3 (residual)** No functional case exercises a real rare-flag build
  (`COLLECT_SP_ACC`/MPFR compiled binaries); the synthetic
  `run_schema_selftest()` covers the collision logic, but an end-to-end
  rare-flag case would be stronger if such a family is ever added to the
  build matrix.
