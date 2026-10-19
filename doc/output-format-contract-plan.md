# Plan: Output format contract (P1 reader-side) + event-stream conversion (P2-lite)

Status: plan — open items only. P1 executed 2026-10-19 (see "Already in
place"); only P2-lite remains open, gated on confirmed demand.

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

## Already in place (do not re-implement)

Pre-existing at the 2026-10-18 re-baseline:

- **Format map**: SKILL.md "File-Type → Reader-Class Quick Reference" plus
  the format column in `assets/data-readback-patterns.md` (commit 63fac41).
- **Reader guards** (SDAR `tools/base.py`, commit eae5961): `fromfile()`
  raises on byte misalignment by default (`strict_mismatch`) with
  reader-kwargs hints; `loadtxt()` detects binary content and raises an
  advisory `ValueError`.
- **Functional readback**: the smoke calls the correct reader per artifact
  (`data.interrupt` via `loadtxt`, strict column mismatch fails the case —
  7dd238e) and verifies ASCII event-stream layouts in the restart flow
  (`ascii_rule`/`ascii_line_time`, a0115ef).

P1 executed 2026-10-19 (SDAR commit 9dac32b; PeTar side: the commit carrying
this plan update):

- **P1.1 `fromfile()` text sniff + auto-dispatch** (SDAR `tools/base.py`):
  `detectAsciiFile()` sniffs the first bytes BEFORE the alignment/
  `strict_mismatch` logic; ASCII content with binary-only kwargs (`offset`,
  `count`) raises an advisory pointing at `loadtxt(skiprows=...)`; ASCII with
  compatible kwargs dispatches to `loadtxt()` with a UserWarning (silent for
  explicit `format='ascii'`; `format='binary'` skips the sniff; `sep` text
  mode untouched). Delivered behavior: the dispatch is a pure
  `loadtxt(**kwargs)` passthrough — headered ASCII files still need
  caller-supplied `skiprows=1`, and a wrong column count then fails loudly
  inside loadtxt. Deployed via `make -C SDAR/tools` — the installed copy at
  `~/include/sdar` is what all consumers import (the source tree is not a
  package).
- **P1.2 Functional format assertions**: smoke `add_check(expect_format=...)`
  sniffs every non-empty readback artifact via `petar.detectAsciiFile` and
  fails the case on a format flip (empty files are exempt; sniffed format is
  recorded in check details); all ASCII event-stream call sites are marked
  `expect_format="ascii"`. A stale installed sdar package degrades to a
  recorded per-case failure instead of aborting the matrix run.
- **P1.3 README link**: `test/functional/README.md` item 6 documents the
  assertion and links the SKILL quick-reference table; the format-choice
  note in `assets/data-readback-patterns.md` describes the reverse guard and
  the `format=` override.

Validation: synthetic dispatch matrix (8 vectors) plus real-file checks — the
2026-10-18 `data.interrupt` vector (ASCII + `fromfile` +
`strict_mismatch=False`) now dispatches with records=1 instead of 2 garbage
records; full functional smoke 9/10 pass (`sevn` skipped: binary family
unavailable, pre-existing); simulated writer flip turns exactly the affected
check red. Independent review verdict: approve-with-nits, nits addressed.

## Open items

### P2-lite — event-stream conversion (only if demand confirmed)

1. **Large event streams get a format choice**: `data.sse*`, `data.bse*`,
   `data.interrupt`, `data.prof.*` (Gyr-scale runs produce GB-size text
   logs). Either a writer switch (mirroring `--write-group-info 1|2`,
   snapshot `-i 0..3`) or transfer-tool support. Constraint added
   2026-10-19: `petar.data.clear` now parses the ASCII event streams with
   per-stream time-column rules (96bb607) — any format switch must keep it
   working (binary variant relies on the P1.1 dispatch). Acceptance:
   round-trip (write A → convert B → parse equals direct parse) proven per
   stream, **and** `petar.data.clear` verified against both formats.

### Explicitly out of scope (revisit only on demand)

- Converters for `status/lagr/core/esc_single/esc_binary/group` (small files,
  conversion value marginal).
- Physical format unification (P3).

## Rough effort

P1: executed 2026-10-19 (~1 session, as estimated).
P2-lite: 1 session per stream family (writer switch + round-trip tests +
`petar.data.clear` compatibility), only after demand is confirmed.
