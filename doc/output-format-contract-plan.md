# Plan: Output format contract (P1 reader-side) + event-stream conversion (P2-lite)

Status: plan — open items only, deferred (2026-10-18 decision: focus stays on
the functional test layer; revisit when functional work settles).

## Context

PeTar outputs are physically mixed: state-like bulk data (snapshots,
`data.status/lagr/core`, escapers, group files, `.single/.binary`, `object.*`)
is BINARY, while append-only event streams (`data.sse*`, `data.bse*`,
`data.interrupt`, `data.prof.*`) are ASCII, and `data.par*` is ASCII text.
The split is defensible (bulk → size/IO; streams → inspectability/append),
but the contract is undocumented and unenforced, which produced three
functional-layer incidents in 2026-10: the group-file glob written from the
help text instead of the real naming, the `object.1` → `object.1.origin`
rename going unnoticed, and the `data.interrupt` bogus pass (ASCII file read
with the binary reader under `strict_mismatch=False`). Design discussion
(2026-10-18) converged on: contract + reader auto-detection first; physical
format unification shelved; conversion support only where demand is real.

## Design decisions (agreed)

- **Contract over unification**: keep the purpose-based split; make it
  explicit, self-enforcing in readers, and asserted by tests.
- **P2 motivation re-scoped**: transfer coverage is NOT a bridge toward a
  format flip (P3 shelved); its only live value is on-demand
  compression/inspection of the large event streams.
- **P3 (single physical format) is out of scope** unless revisited; if ever
  revived, the full transfer bridge must exist first plus a dual-support
  release cycle.

## Open items

### P1 — reader-side contract (highest value, do first)

1. **Format map table**: one table in the SKILL (or its asset) listing every
   output file → format (ASCII/BINARY) → correct reader call. Acceptance: a
   user never guesses `fromfile` vs `loadtxt` again; functional README links it.
2. **Reader auto-detection** in `tools/analysis/base.py`: `fromfile()` sniffs
   the first bytes (printable text vs binary header) and dispatches to
   `loadtxt` automatically; explicit format still wins; mismatch warnings
   unchanged. Acceptance: reading any PeTar output with the "wrong" API still
   returns correct data; backward compatible for all existing callers.
3. **Functional format assertions**: extend the smoke readback to assert each
   artifact's expected format (text vs binary sniff) so renames/format flips
   fail loudly instead of surfacing as reader confusion.

### P2-lite — event-stream conversion (only if demand confirmed)

4. **Large event streams get a format choice**: `data.sse*`, `data.bse*`,
   `data.interrupt`, `data.prof.*` (Gyr-scale runs produce GB-size text
   logs). Either a writer switch (mirroring `--write-group-info 1|2`,
   snapshot `-i 0..3`) or transfer-tool support. Acceptance: round-trip
   (write A → convert B → parse equals direct parse) proven per stream; the
   Python readers handle both formats via the P1 auto-detect.

### Explicitly out of scope (revisit only on demand)

- Converters for `status/lagr/core/esc_single/esc_binary/group` (small files,
  conversion value marginal).
- Physical format unification (P3).

## Rough effort

P1: 1–2 sessions (base.py sniffing + tests; format map; assertions).
P2-lite: 1 session per stream family (writer switch + round-trip tests).
