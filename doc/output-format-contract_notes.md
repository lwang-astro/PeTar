# Output format contract — design record

Status: reference note — documents the output-format contract, its
enforcement, and the shelved conversion/unification decisions. P1
(reader-side contract) executed 2026-10-19; P2-lite (event-stream
conversion) shelved 2026-10-26 as low-necessity; P3 (physical
unification) rejected. Supersedes the retired `output-format-contract-plan.md`.

## The contract

PeTar outputs are deliberately mixed: state-like bulk data (snapshots,
`data.status/lagr/core`, escapers, group files, `.single/.binary`,
`object.*`) is BINARY (size/IO); append-only event streams (`data.sse*`,
`data.bse*`, `data.interrupt`, `data.prof.*`) are ASCII
(inspectability/append); `data.par*` is ASCII text. The authoritative
file → format → reader mapping is the SKILL "File-Type → Reader-Class
Quick Reference" table plus `assets/data-readback-patterns.md`.

## Enforcement (executed 2026-10-19)

- Reader auto-detection: SDAR `tools/base.py` `fromfile()` sniffs the
  first bytes and dispatches ASCII tables to `loadtxt()` BEFORE the
  byte-alignment/`strict_mismatch` logic, so a wrong-reader call can no
  longer return garbage records (SDAR commit 9dac32b).
- Test-side assertion: the functional smoke `add_check(expect_format=...)`
  sniffs every non-empty readback artifact and fails the case on a writer
  format flip (PeTar commit 68deaa7).
- Known enforcement gap: only artifacts passing through the smoke readback
  are format-asserted; `data.par*`, `snap.lst`, `.bk` backups, and
  `.randseeds` are outside that path (par files are exercised only by the
  `petar.read.par` tool smoke, which asserts nothing about format), so a
  hypothetical par-writer flip would currently stay silent.

## Shelved decisions (recorded so they are not re-litigated blind)

- **P2-lite (event-stream conversion) — shelved 2026-10-26, demand judged
  low.** If ever revived: give the large streams (`data.sse*`, `data.bse*`,
  `data.interrupt`, `data.prof.*`) a format choice via a writer switch
  (mirroring `--write-group-info 1|2`, snapshot `-i 0..3`) or
  transfer-tool support. Acceptance per stream family: round-trip (write
  A → convert B → parse equals direct parse) AND `petar.data.clear`
  verified against both formats — its ASCII time-column rules (96bb607)
  are load-bearing, and a binary variant must ride the `fromfile()`
  auto-dispatch.
- **P3 (single physical format) — rejected unless demand changes.** Revive
  only after a full transfer bridge exists plus a dual-support release
  cycle; conversion coverage alone is not a bridge toward a format flip.
