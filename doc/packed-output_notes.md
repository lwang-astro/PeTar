# Packed snapshot container (`-w 4`) — implementation notes

Status: reference note — documents the packed snapshot output mechanism
(solver option `-w 4`, container/index layout, transactional commit
integration, restart semantics, and the `petar.pack` tool family) as
implemented in 1875e. Design decision: chosen over an HDF5 runtime container
(HDF5 has no journaling, so a monolithic file risks losing all prior
snapshots on an interrupted write; zero new dependencies; existing binary
record format and readers stay authoritative). HDF5 remains available as an
offline conversion (`petar.pack hdf5`).

## Motivation

The default output style (`-w 1`) writes one snapshot file per output time;
long runs and large run grids produce thousands-to-millions of files, making
directory access and copying slow and stressing filesystem inode quotas.
The packed style keeps the on-disk record format unchanged but concatenates
records into a few append-only volume files with a small line-per-snapshot
index.

## Container layout

- `[prefix].pack.<id>` (e.g. `data.pack.0000`): append-only concatenation of
  standalone BINARY snapshot records (one `FileHeader` + N particle records),
  byte-identical to what `-w 1` writes into `data.<N>`.
- `[prefix].pack.idx`: one line per packed record,
  `nfile time n_body volume nbytes offset` (`offset` is the record start
  within the volume file; headerless like the other output streams).
- Volume rolling: `--pack-volume-size <GiB>` starts a new volume when the
  next record would exceed the limit; `0` (default) = single unlimited
  volume.
- Per-snapshot sidecars are **not** packed and keep their `-w 1` names:
  `data.<N>.randseeds` (BSE/DSM; restart RNG continuity) and
  `data.<N>.galpy`.
- `data.meta` gains `snapshot_container = file|pack` so tooling can detect
  the container; all other schema keys are unchanged.
- `-w 4` requires a BINARY snapshot write format (`-i 0` or `-i 2`); other
  write formats abort at startup with an explicit message.
- No `data.snap.lst` is written at runtime in packed mode — the index
  replaces it (`petar.pack extract --all --snap-lst` regenerates the list
  from the index on demand).

## Transactional commit integration

- The snapshot **record bytes** are appended to the volume file at output
  time: the FDPS writer stages the complete record in
  `[prefix].pack.staging.<nfile>.tmp` (gathering particles to rank 0 under
  MPI), `OutputCommitManager::appendFile` moves it into the volume, and the
  staging file is removed — one staging file per record, keyed by `nfile`.
- The **index line is the commit record**. Index lines accumulate in
  `[prefix].pack.idx.tmp` (rank 0), which is registered with
  `OutputCommitManager` in `Append` mode from `registerTransactionalOutputs()`
  and committed with each output window alongside the status row. A packed
  snapshot becomes visible to readers exactly when its index line commits.
- Immediate append is load-bearing: the initial t=0 `output()` call runs
  *outside* a commit window, and `beginWindow()` discards uncommitted
  registrations — a batched staging file would be silently dropped, leaving
  every index row misaligned by one (found by the `std-packed` functional
  smoke, 2026-11-06; per-record staging alone did not fix it; final design
  = immediate append + index-line commit).
- Crash healing on startup: after recovering the rolling state from the
  last index line, the current volume is truncated to the committed byte
  count — a partial record append from an interrupted run is dropped before
  it can corrupt later offsets. Uncommitted staging/index tmp files are
  removed by the standard `--keep-tmp-on-startup` policy. A missing volume
  referenced by the index aborts loudly (recover via `petar.pack
  rebuild`/`extract` or remove the container).
- The index remains derivable from the volumes alone (fixed-size headers +
  known record size), so `petar.pack rebuild` regenerates it after idx-file
  loss.

## Restart semantics

- Restart input: `petar.pack extract <prefix> --last` materializes
  `data.<N>`, then the canonical restart form applies unchanged.
- Append restart (`-a 1`, default): the rolling state (volume id, committed
  bytes) is recovered from the last index line during `readParameters`
  (before the input snapshot is read); a missing index starts an empty
  container (snapshot numbering continues from the restart header). The
  stale-trailing-records warning fires later, in `initialParameters`, where
  the restart snapshot header (`file_header.nfile`) is valid.
- Restarting from an **older** snapshot leaves trailing packed records
  stale; the solver prints a warning suggesting
  `petar.pack trim <prefix> -t <time>` (the container analog of
  `petar.data.clear` for event streams). Trim truncates the last kept
  record's volume, deletes higher volumes, and rewrites the index.
- **Untrimmed stale-tail semantics**: restarting without trimming keeps the
  stale tail and re-appends from the restart nfile, producing *duplicate*
  nfile rows in the index (offsets stay self-consistent). Consumers then
  differ by design: `PackedSnapshots.read(nfile=…)` selects the FIRST
  (stale) match, `petar.pack extract --all` overwrites files so the LAST
  record wins, and `petar.pack rebuild` refuses the non-monotonic sequence
  with a trim hint — trim before restarting to keep the index unique.
- Append restart with a missing/empty index but **non-empty volume files**
  aborts loudly (no committed line references those bytes; they would
  silently misalign every appended record) — recover via
  `petar.pack rebuild` or remove the container.
- Fresh start (`-a 0`): all volumes and the index are removed.

## Readers and tools

- `petar.PackedSnapshots(prefix)` (`tools/analysis/data.py`): resolves the
  reader layout from `data.meta` (or explicit kwargs), then reads any
  snapshot by seek — `readSnapshot(nfile=…)`, `readSnapshot(time=…)`,
  `pack[-1]` — returning a `petar.Snapshot` object (`.particle` +
  `.header`, built via `Snapshot.fromRecord`) through the existing binary
  `fromfile(offset, count)` path; no unpacking needed.
- `petar.pack` (`tools/pack_manager.py`): `list`, `extract`
  (byte-identical `data.<N>`; `--snap-lst` rewrites the snapshot list from
  the full index, requires `--all`, idempotent), `trim`, `rebuild`,
  `hdf5` (offline single-file conversion, groups `snap_<nfile>`,
  optional compression; lazy `h5py` import).
- Because extracted records are byte-identical standalone snapshots, the
  whole readback contract (sidecar/auto-detect readers, `petar.movie`,
  `petar.data.process`) applies to them unchanged.

## Retained limitations

- BINARY snapshot write only (`-i 0`/`-i 2`).
- Volume rolling is decided per record from the staged size; a single
  record larger than the limit still occupies its own (oversized) volume.
- The `.randseeds` / `.galpy` per-snapshot sidecars still create one small
  file per output time in the corresponding solver families.
