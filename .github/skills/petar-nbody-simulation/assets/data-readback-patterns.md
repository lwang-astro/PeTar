# Python Data Readback Patterns

This file documents the verified patterns for reading PeTar output files from Python (`import petar`). These patterns have been iterated and validated through `test/functional/run_functional_smoke.py` and should be used when an agent needs to read or validate simulation outputs.

**Reference implementation**: `test/functional/run_functional_smoke.py` function `run_python_output_read_checks()`.

## Prerequisites

```python
import sys
sys.path.insert(0, "<repo_root>/tools")  # or ensure petar is in PYTHONPATH
import petar
import numpy as np
```

## Common Parameters

### Keyword arguments

Keyword arguments control column layout and must match the producing solver:

| Argument | Values | Purpose |
|----------|--------|---------|
| `interrupt_mode` | `none`, `merger`, `bse`, `bseEmp`, `mobse`, `dsm` | Must match solver `--with-interrupt` |
| `external_mode` | `none`, `galpy`, `agama` | Must match solver `--with-external`; adds `pot_ext` (and snapshot header c.m. offsets) |
| `spin_3d` | bool (`True`) | BSE-family `star.spin` width: 3 columns (current) vs 1 column. **Legacy outputs (before 2024-12): use `spin_3d=False`** |
| `use_mpfrc` | bool (`False`) | MPFR build: adds `pos_high` (3 columns) |
| `collect_sp_acc` | bool (`False`) | Solver built with superparticle acc collection: adds `acc_sp` (3 columns) |
| `N_particle` | int (`0`) | `Status` only: number of per-particle blocks appended after the global block |
| `N` | int | `GroupInfo` only: group member count encoded in the filename (`data.group.n<N>`) |

These apply to every particle-composing reader (`Particle`, `Status`, `SingleEscaper`/`BinaryEscaper`, `GroupInfo`, `Binary` members). Profile and Lagrangian readers take their own layout kwargs (`use_gpu`, `FDPS_version`, `old_version`, `group_count`; `mass_fraction`, `calc_energy`, `add_star_type`, `add_mass_range`) — see the class docstrings.

Pass these consistently to all readers. **Do not guess them by trial and error** — recover them with auto-detection (next section) and pass explicit kwargs only to override.

### Schema auto-detection (`petar.getDataKwargs` / `petar.Snapshot`)

The layout kwargs are recoverable from the data itself, in this precedence (highest wins): **explicit kwargs → `<prefix>.meta` sidecar → structure sniffing**.

**Concise construction (preferred)**: every public petar reader accepts a file path as its first constructor argument and auto-reads it with the auto-detected layout — `petar.Particle('data.5')`, `petar.Particle('data.5.single')`, `petar.Binary('data.5.binary', G=...)`, `petar.SingleEscaper('data.esc_single')`, `petar.BinaryEscaper('data.esc_binary')`, `petar.GroupInfo('data.group.n2', N=2)`, `petar.LagrangianMultiple('data.lagr')`, `petar.Core('data.core')`, `petar.BSEStatus('data.bse_status')`, and the SSE/BSE event readers (`petar.SSETypeChange('data.sse.type_change')`, `petar.BSEKick('data.bse.gw_kick')`, ...). Explicit kwargs passed alongside win over detection; `.npy` / ASCII / binary (with or without snapshot header) are dispatched automatically. No explicit kwargs are needed for default-build data. Exceptions kept explicit: `Status` (the optional per-particle blocks need `N_particle`) and `object.*` (needs the extra `time` member added first). `petar.Profile(path)` also auto-reads (see Pattern 8).

**`meta=` injection** (concise constructor and `petar.getDataKwargs`): `meta='<path>'` uses an explicit sidecar — the escape hatch for files copied away from their run directory; `meta={...}` passes ready layout kwargs as a group (no detection, no warnings); `meta=False` disables the sidecar (structure/explicit only). Default `'auto'` locates `<prefix>.meta` by the generic prefix walk (successively shorter stems: `data.5.single` → `data.5.meta`, `data.meta`) — no per-file-type name mapping; works for any run prefix and all stream files.

- **`petar.getDataKwargs(path)`** returns the kwargs dict for any run-family file (`data.5`, `data.5.single`, `data.esc_single`, ...). It reads the solver-written sidecar `<prefix>.meta` (petar-meta/1, written next to `data.par`; absent for older runs) and cross-validates it against the file structure — a disagreement raises `ValueError` (the sidecar was likely copied from another run). Without a sidecar it sniffs the structure: ASCII header item count or binary header size (24/72 B; 32/128 B for float128) pins `external_mode`; the column count (ASCII) or bytes-per-record with N from the header (binary) pins the rest.
- **`petar.Snapshot(path)`** returns an object with `.particle` and `.header` for raw snapshots in one call (the concise `petar.Particle(path)` returns particles only).
- **CLI tools auto-fill**: `petar.data.process`, `petar.get.object.snap`, `petar.format.transfer.post`, and `petar.movie` detect the schema from the data when `-i`/`-t` (and `--use-mpfrc`/`--spin-1d`) are omitted, print what was detected, and honor explicit flags over detection.

Detection uniqueness (verified over the kwargs space, `external_mode` pinned by the header):

| interrupt | no external | external |
|---|---|---|
| none | 20 | 21 |
| merger | 24 | 25 |
| dsm | 30 | 31 |
| bse-family, spin 1d (pre-2024-12) | 34 | 35 |
| bse-family, spin 3d | 36 | 37 |
| sevn | 42 | 43 |

Default builds (no MPFR, no `COLLECT_SP_ACC`) are always uniquely determined. Layouts that only differ in rare build flags (`use_mpfrc` vs `collect_sp_acc`, both +3 columns) can collide on column count: detection then resolves to the **lowest plausibility tier** (fewest rare flags) with a `UserWarning` listing the suppressed alternatives, and a same-tier collision raises `AmbiguousSchemaError` listing the candidates — never a silent guess. On that error, pass explicit kwargs or provide/regenerate `data.meta` (`petar.autoFillKwargs(path, kwargs)` is the programmatic/CLI helper). Sidecar cross-validation only checks structure-unique facts (header-pinned `external_mode`/`float_type`, and that the sidecar layout is among the raw candidates), so a correct rare-flag sidecar is never rejected by a tier-0-resolved sniff.

Caveat (binary, contrived N): when two header variants divide the file size exactly (N satisfying `H′−H = 8N·Δncol`, e.g. a float64/float128 pair), aliasing is possible and the tier-0 float64 candidate wins with a warning; the sidecar pins `float_type` authoritatively. ASCII is immune (no precision in text columns).

**When no candidate matches (no-match policy)**: `petar.getDataKwargs(path)` raises `AmbiguousSchemaError` listing the candidates (or reporting a schema gap). In the concise constructor the outcome depends on the file family: for **snapshot-family files** (snapshots, derived `.single/.binary/...`, run streams) a detection failure **raises** unless both `interrupt_mode` and `external_mode` are given explicitly — a default-kwargs fallback could misparse silently when the byte size happens to align (e.g. an isolated `.single` copied away from its run directory: 27360 B divides by both the 288 B bse and the 160 B default records). For **other file types** (e.g. hermite `.log` dumps, which carry no detectable schema) the constructor warns and proceeds; the read then surfaces the precise sdar error if the layout is wrong. Weak evidence gets strict treatment: the headerless divisibility scan requires a **unique** match (no plausibility-tier resolution) precisely because several layouts can divide one file size.

**Maintenance model (where the layout knowledge lives)**: no dedicated schema module exists - the knowledge is distributed to its owning classes. Reader classes declare their layout-affecting constructor arguments in `LAYOUT_PARAMS` (first value = default/most plausible), and the generic engine in `DictNpArrayMix` (sdar `tools/base.py`: `layoutCandidates` / `layoutTier` / `_layoutMatches` / `detectLayoutKwargs`, sidecar `findLayoutMeta`/`readLayoutMeta`/`writeLayoutMeta`, concise constructor `resolvePathInput`/`readAutoFile`/`enablePathInit`, `AmbiguousSchemaError`, `autoFillKwargs`) enumerates the declared combinations and resolves them against structural measurements. Petar adds `ParticleRecordLayoutMixin` (`tools/analysis/data.py`): snapshot structure detection (header knowledge in `PeTarDataHeader.asciiHeaderInfo`/`binaryHeaderVariants`), sibling-snapshot fallback, the strict no-silent-misparse policy, and sidecar-vs-structure cross-validation - mixed into the particle chain, `Binary`, escapers and `GroupInfo`, so every class detects with its OWN record layout. **Per-run (non-enumerable) layout parameters**: readers whose layout depends on data.process-time options declare `LAYOUT_META_PARAMS` (always list-valued: mass-fraction arrays, star-type / mass-range name lists) - carried by the sidecar with comma-encoded tokens, never enumerated as candidates. `petar.data.process` writes the producer sidecar `<prefix>.lagr.meta` next to the Lagrangian table (the prefix walk finds it before the solver's `data.meta`); `LagrangianMultiple('data.lagr')` then auto-reads with the exact mass fractions / energy columns. Old lagr files without the sidecar need explicit kwargs as before.
The mixin code itself contains NO parameter names; the remaining role knowledge is declared next to `LAYOUT_PARAMS` on the classes: `HEADER_PINNED_PARAMS` (header fact -> parameter, e.g. `'cm': 'external_mode'`), `LABEL_REFINE_PARAMS` (family representatives refined by run-dir probing, e.g. galpy/agama), `PINNING_PARAMS` (explicit values that exempt a detection failure from the strict raise). `petar.Snapshot(path)` bundles header + particles. **Adding or changing a layout parameter to a reader = update its `LAYOUT_PARAMS` (and the keys it builds) - nothing else.** The public per-class entry point `Reader.detectLayoutKwargs(ncols=..., record_size=...)` works for any declared class.

### Column mismatch: stop rules

When any reader warns or errors on column count/shape (`mismatches the number of columns`, `IndexError` on the column axis, dtype-alignment warnings):

1. Enumerate the kwargs table above against the producing solver, **including the version-dependent rows** (`spin_3d`, `old_version`-style flags). Better: call `petar.getDataKwargs(path)` first — if detection succeeds, its kwargs are the combination the file structure supports.
2. If no kwargs combination matches, **stop the analysis and report to the user** with the file pattern, reader class, kwargs tried, and error text. The file is likely from a PeTar version with a different output schema. Known legacy options: pre-2024-12 BSE-family outputs → `spin_3d=False`; pre-2020-09 BSE merger records (`BSEDynamicMerge`) → `less_output=True`; pre-2020-12 galpy snapshots → `petar.format.transfer -c`; pre-2020 group data → `petar.format.transfer -g`.
3. **Never write a custom parser or trim/slice columns as a workaround** for snapshot/status/escaper/group files — silent misalignment corrupts physics results. This restates `Snapshot Read-Mismatch Policy` in `SKILL.md`, which owns the full policy.

### Input/output format: choosing the right read method

`petar.data.process` supports three output formats, each requiring a different read method:

| `output_format` | File extension | Read method | Write method |
|-----------------|----------------|-------------|--------------|
| `binary` (default) | (none) | `.fromfile(path)` | `.tofile(f)` |
| `npy` | `.npy` | **`.load(path)`** | `.save(f)` |
| `ascii` | (none) | `.loadtxt(path)` | `.savetxt(f)` |

When reading post-processed files, **check the file extension** to determine the format. If the file ends with `.npy`, use `.load()` instead of `.fromfile()`. Otherwise, use the method matching the table above. Calling `.loadtxt()` on a BINARY file raises a clear `ValueError` (use `fromfile()`). The reverse is also guarded: `fromfile()` sniffs the first bytes and, on an ASCII text table, dispatches to `loadtxt()` with a `UserWarning` — even under `strict_mismatch=False` — so a wrong-reader call can no longer return garbage records. Passing binary-only kwargs (`offset`, `count`) with ASCII content raises an advisory instead. Use `format='binary'`/`'ascii'` in `fromfile()` to force a physical format and skip the sniffing (no warning); the shared sniff helper is `petar.detectAsciiFile(path)`.

### Unit conversion (critical for scientific correctness)

**Never use approximate hardcoded constants** (e.g., `semi * 206265`). These introduce systematic bias that can corrupt scientific results.

PeTar's physical constants (`src/astro_units.hpp`) are now **auto-generated from `astropy`** via `tools/generate_astro_units.py`. This means `petar` module constants and `astropy.units` derive from the same IAU/CODATA sources — there is no precision hierarchy between them.

**Preferred approach — `astropy.units`** (cleaner API, handles any conversion in one line):
```python
import astropy.units as u

# pc → AU
semi_au = (binary.semi * u.pc).to(u.AU).value

# pc/Myr → km/s
vel_kms = (binary.vel * u.pc / u.Myr).to(u.km / u.s).value

# pc → km, Rsun, ly, etc.
semi_km = (binary.semi * u.pc).to(u.km).value
semi_rsun = (binary.semi * u.pc).to(u.Rsun).value
```

**Alternative — `petar` module constants** (same values, use when cross-checking against C++ solver output):
```python
# Gravitational constant (Msun, pc, Myr) — from astropy via astro_units.hpp
petar.G_MSUN_PC_MYR          # = 0.004498502151469553

# pc/Myr ↔ km/s — from astropy via astro_units.hpp
petar.PCMYR_TO_KMS           # = 0.9777922216807892
petar.KMS_TO_PCMYR           # = 1.022712165045695
```

## Pattern 1: Lagrangian Radii (`data.lagr`)

```python
# current data.process writes the producer sidecar <prefix>.lagr.meta:
# the path constructor auto-recovers mass_fraction, calc_energy,
# add_star_type, add_mass_range, external_mode - no kwargs needed
lagr = petar.LagrangianMultiple("data.lagr")

# older runs (no sidecar): pass the per-run layout kwargs explicitly
lagr = petar.LagrangianMultiple(external_mode="<mode>",
                                mass_fraction=np.array([0.1, 0.3, 0.5, 0.7, 0.9]))
lagr.fromfile("data.lagr")
```

- `external_mode` gates `epot_ext` (external potential energy), which exists only when `calc_energy=True` — it does **not** add COM-offset columns to the Lagrangian table.
- Structure: `time` + three `Lagrangian` blocks (`single`, `binary`, `all`; binaries counted once at their c.m.), each holding `r`, `m`, `n`, `vel`, `sigma` (+ `epot`/`vr`/`epot_ext` with `calc_energy`); extra blocks per `add_star_type` / `add_mass_range` entry — star-type entries add a member named **exactly as given** (e.g. `BH__in__MS`), mass-range entries get the `mass_` prefix (`0.1_1` → `mass_0.1_1`; cross ranges → `mass_A__in__mass_B`).
- **Column semantics** (per `tools/analysis/lagrangian.py`): **6 columns** — mass fractions `[10%, 30%, 50%, 70%, 90%]` **plus the core radius as the last column** (Casertano & Hut 1985). The velocity keys (`abs`, `x`, `y`, `z`, `rad`, `tan`, `rot`) likewise carry 6 entries, last = core. With a custom `--mass-fraction` list the column count changes accordingly.
- **`all.m` is the *average stellar mass* per particle within each radius, not a cumulative mass** (cumulative mass is divided by the particle count). Multiply by `all.n` to recover the enclosed mass; under `mode="shell"` both refer to the shell.
- **Mislabel hazard without the sidecar**: a lagr table produced with non-default `--mass-fraction` of the *same length* reads successfully without the sidecar but silently labels the default fractions — keep `<prefix>.lagr.meta` next to the file or pass `mass_fraction=` explicitly (the functional `lagr.variants` step pins this contract).

## Pattern 2: Core Properties (`data.core`)

```python
core = petar.Core()
core.fromfile("data.core")
```

## Pattern 3: Run Status (`data.status`)

```python
status = petar.Status()
status.fromfile("data.status")
```

- Add `N_particle=<n>` when the run wrote per-particle blocks after the global block (e.g. single-particle orbit tests use `N_particle=1`); the blocks appear as `status.particles.p<i>`.
- Legacy BSE-family outputs (before 2024-12) additionally need `spin_3d=False` — see "Column mismatch: stop rules".

## Pattern 4: Escapers

```python
# Single escapers
esc_single = petar.SingleEscaper(interrupt_mode="<mode>", external_mode="<mode>")
esc_single.fromfile("data.esc_single")

# Binary escapers
esc_binary = petar.BinaryEscaper(interrupt_mode="<mode>", external_mode="<mode>")
esc_binary.fromfile("data.esc_binary")
```

## Pattern 5: Stellar Evolution Tables

```python
# SSE (single star evolution)
sse = petar.SSEType()
sse.loadtxt("data.sse")

sse_type_change = petar.SSETypeChange()
sse_type_change.loadtxt("data.sse.type_change")

sse_sn_kick = petar.SSESNKick()
sse_sn_kick.loadtxt("data.sse.sn_kick")

# BSE (binary evolution)
bse = petar.BSEType()
bse.loadtxt("data.bse")

bse_status = petar.BSEStatus()
bse_status.fromfile("data.bse_status")

bse_type_change = petar.BSETypeChange()
bse_type_change.loadtxt("data.bse.type_change")

bse_sn_kick = petar.BSEKick()
bse_sn_kick.loadtxt("data.bse.sn_kick")

bse_gw_kick = petar.BSEKick()
bse_gw_kick.loadtxt("data.bse.gw_kick")

bse_dynamic_merge = petar.BSEDynamicMerge()
bse_dynamic_merge.loadtxt("data.bse.dynamic_merge")

bse_binary_merge = petar.BSETypeChange()  # Note: uses BSETypeChange, not a dedicated class
bse_binary_merge.loadtxt("data.bse.binary_merge")
```

## Pattern 6: Raw Snapshots (`data.<N>`)

### Header offset selection

The snapshot header offset depends on the build features (external potential, MPFRC precision):

| Build | Offset |
|---|---|
| default | `petar.HEADER_OFFSET` |
| `--with-external=galpy\|agama` | `petar.HEADER_OFFSET_WITH_CM` |
| `--enable-mpfrc` | `petar.HEADER_OFFSET_F128` |
| `--enable-mpfrc` + external potential | `petar.HEADER_OFFSET_WITH_CM_F128` |

```python
offset = petar.HEADER_OFFSET_WITH_CM if external_mode != "none" else petar.HEADER_OFFSET
# MPFRC builds use the *_F128 variants instead.
```

The header itself (time, N, optional COM offset) is read with `petar.PeTarDataHeader`, which also
takes `external_mode=`. ASCII snapshots (`-u`-independent) have a one-line header — use
`skiprows=1` with `.loadtxt()`.

### Reading particles

**Default read (auto-detected layout, header included)** — the same call the
functional report shows as the primary method:

```python
snap = petar.Snapshot("data.1")          # layout from <prefix>.meta sidecar,
                                          # structure sniff fallback
particle, header = snap.particle, snap.header
```

**Legacy explicit read (manual offset, particles only)** — for column subsets
or pre-sidecar files:

```python
particle = petar.Particle(interrupt_mode="<mode>", external_mode="<mode>")
particle.fromfile("data.1", offset=offset)
```

### Reading single/binary component snapshots

After `petar.data.process`, each snapshot has split files. **Choose read method by file extension:**

**Binary format (default, no `.npy` suffix):**

```python
# Single particles
single = petar.Particle(interrupt_mode="<mode>", external_mode="<mode>")
single.fromfile("data.1.single")

# Binary particles (requires G for orbital element computation)
binary = petar.Binary(
    member_particle_type=petar.Particle,
    interrupt_mode="<mode>",
    external_mode="<mode>",
    G=petar.G_MSUN_PC_MYR,
)
binary.fromfile("data.1.binary")
```

**npy format (`.npy` suffix):** use `.load()` instead of `.fromfile()`:

```python
# Single particles
single = petar.Particle(interrupt_mode="<mode>", external_mode="<mode>")
single.load("data.1.single.npy")        # .load() not .fromfile()

# Binary particles
binary = petar.Binary(
    member_particle_type=petar.Particle,
    interrupt_mode="<mode>",
    external_mode="<mode>",
    G=petar.G_MSUN_PC_MYR,
)
binary.load("data.1.binary.npy")        # .load() not .fromfile()
```

## Pattern 6b: Packed Snapshot Container (`data.pack.<id>` + `data.pack.idx`)

Runs with `-w 4` store snapshots as records inside append-only volume files
(`data.pack.0000`, …) with a line-per-snapshot index (`data.pack.idx`).
A volume is NOT a standalone snapshot file — the structure sniffing cannot
identify records in it; always go through `petar.PackedSnapshots`, which
resolves the layout from the run sidecar (`<prefix>.meta`) and seeks records
by the index (the record byte range is exact; tail-alignment handling is
internal — do NOT add `offset=`/`strict_mismatch` yourself).

```python
import numpy as np
import petar

pack = petar.PackedSnapshots("data")     # prefix, not a file path
pack.index["time"]                       # all output times (discovery)

# selection forms (all return a petar.Snapshot: .particle + .header)
snap = pack.readSnapshot(time=3.0)   # FIRST snapshot at/after t=3.0
                                     # (falls back to the last row; NOT nearest)
snap = pack.readSnapshot(nfile=3)    # by snapshot id (== data.3)
snap = pack[-1]                      # last record; pack[i] by row

# range query: mask the index, then read by row
rows = np.where((pack.index["time"] >= 2.0) & (pack.index["time"] <= 5.0))[0]
snaps = [pack.readSnapshot(row=int(r)) for r in rows]

for i in range(len(pack)):          # or iterate all
    snap = pack[i]

pack.readHeader(row=0)               # header only (fast: 24/72 bytes)
```

Materialize standalone snapshots for tool chains (`petar.movie`,
`petar.data.process`) or external software — the extracted files are
byte-identical records and read with the standard `petar.Snapshot(path)`:

```bash
petar.pack extract data --last            # or --nfile <N> / --time <t>
petar.pack extract data --all --snap-lst  # all + data.snap.lst for tools
```

Restart from a packed run: extract first (see SKILL "Restart/resume" item 8);
`petar.pack trim <prefix> -t <time>` replaces `petar.data.clear` for the
container when restarting from an older snapshot.

## Pattern 7: Object Snapshots (`object.<N>`)

Object snapshots produced by `petar.get.object.snap` require pre-registering a `time` member before reading, because the object time-series format includes a time column not in the standard Particle layout:

```python
obj = petar.Particle(interrupt_mode="<mode>", external_mode="<mode>")
obj.addNewMember("time", np.array([], dtype=float))
obj.fromfile("object.1")
```

This pattern is required for all interrupt/external mode combinations.

## Pattern 8: Profile Tables (`data.prof.rank.*`)

Profiles are ASCII tables whose column count depends on the compile-time configuration (GPU, FDPS version, old_version flags). **Use the concise constructor**: it detects the layout from the numeric body's column count.

```python
prof = petar.Profile("data.prof.rank.0")
```

- Output streams are uniformly headerless (like status/lagr/core/esc). Files from older versions (2020-02 to 1868e) may carry a column-name first line — the generic path-constructor title skip handles it for any ASCII file, no prof-specific logic.
- The detected layout covers `use_gpu` / `FDPS_version` (<8.0 key set vs 8.0) / `old_version` / `group_count`. Rare-flag combinations that collide on column count resolve to the fewest-rare-flags candidate with a warning (same-tier collisions raise `AmbiguousSchemaError`); GPU builds additionally add the `gpu` member — verify with `prof.ncols` when in doubt.
- Manual matching (legacy pattern): construct `petar.Profile(**kwargs)` candidates and compare `prof.ncols` against the first data row's column count.

### Column-count mismatch fallback

**Scope: `data.profile` files only** — never generalize this column trim to snapshot/status/escaper/group files (see "Column mismatch: stop rules").

If a profile has additional diagnostic columns in later rows (after a consistent numeric prefix), parse the minimum common column count (uniformly headerless; legacy files with a title line need one extra skipped row):

```python
min_cols = 0
with path.open("r") as fh:
    for line in fh:
        tokens = line.strip().split()
        if not tokens:
            continue
        try:
            float(tokens[0])
        except ValueError:
            continue  # title line (legacy files)
        row_cols = len(tokens)
        min_cols = row_cols if min_cols == 0 else min(min_cols, row_cols)

if min_cols > 0:
    np.loadtxt(str(path), usecols=tuple(range(min_cols)))
```

## Pattern 9: GroupInfo (`data.group.n<N>`)

Groups files encode the member count `N` in the filename suffix (e.g., `data.group.n2` for 2-member groups). The `N` parameter is required:

```python
import re

match = re.search(r"\.n(\d+)$", group_path.name)
n_member = int(match.group(1))

group = petar.GroupInfo(
    N=n_member,
    interrupt_mode="<mode>",
    external_mode="<mode>",
)
group.fromfile(str(group_path))
```

## Pattern 10: DSM InterruptBinary (`data.interrupt`)

The DSM interrupt binary may include trailing padding bytes depending on compile-time layout. Validate by record count rather than strict dtype alignment:

```python
interrupt = petar.InterruptBinary(
    particle_type=petar.HardParticle,
    interrupt_mode="dsm",
)
interrupt.fromfile("data.interrupt", strict_mismatch=False)

# Check: records > 0 is a valid read; a trailing-byte warning is expected and non-blocking
```

## Warning Classification

When running read checks, warnings fall into two categories:

### Blocking (treat as read failure)

- `Binary file size is not aligned with dtype itemsize` — binary schema mismatch; raises `ValueError` by default with reader-kwargs hints
- `not aligned with dtype itemsize` — binary layout mismatch
- `File may be truncated` — incomplete snapshot
- `dtype mismatch` — column schema mismatch
- `column` mismatch (shape, column count, column dtype) — column layout mismatch

These indicate a mode-flag or format mismatch. Stop, correct flags, retry. If corrected flags still mismatch, follow **Column mismatch: stop rules** — report, do not work around.

### Non-blocking (record but do not fail)

- `ffmpeg` macro-block resize warnings from `petar.movie`
- Trailing-byte tolerance via `strict_mismatch=False` (reads complete records, warns once): DSM `data.interrupt`, and `petar.data` crash-recovery partial files

## Complete Workflow Pattern

Combining all patterns for a full functional-style readback:

```python
import petar, numpy as np, re, warnings
from pathlib import Path

def readback(case_dir, output_prefix, interrupt_mode, external_mode):
    case_dir = Path(case_dir)

    def check(label, path, loader_fn, command, strict=True):
        if not path.exists():
            return
        with warnings.catch_warnings(record=True) as rec:
            warnings.simplefilter("always")
            ok = True
            try:
                result = loader_fn(path)
            except Exception as e:
                ok = False
                result = str(e)
        # ... evaluate ok + warnings ...

    # Lagrangian
    check("lagr", case_dir / f"{output_prefix}.lagr",
          lambda p: petar.LagrangianMultiple(external_mode=external_mode).fromfile(str(p)),
          f"LagrangianMultiple(external_mode={external_mode})")

    # Core
    check("core", case_dir / f"{output_prefix}.core",
          lambda p: petar.Core().fromfile(str(p)),
          "Core()")

    # Status
    check("status", case_dir / f"{output_prefix}.status",
          lambda p: petar.Status().fromfile(str(p)),
          "Status()")

    # Escapers
    check("esc_single", case_dir / f"{output_prefix}.esc_single",
          lambda p: petar.SingleEscaper(interrupt_mode=interrupt_mode, external_mode=external_mode).fromfile(str(p)),
          f"SingleEscaper({interrupt_mode}, {external_mode})")

    check("esc_binary", case_dir / f"{output_prefix}.esc_binary",
          lambda p: petar.BinaryEscaper(interrupt_mode=interrupt_mode, external_mode=external_mode).fromfile(str(p)),
          f"BinaryEscaper({interrupt_mode}, {external_mode})")

    # SSE/BSE tables (if present)
    # Each entry is (suffix, class, method_name) — method_name must be 'loadtxt' or 'fromfile'
    sse_bse_entries = [
        ("sse", petar.SSEType, "loadtxt"),
        ("sse.type_change", petar.SSETypeChange, "loadtxt"),
        ("sse.sn_kick", petar.SSESNKick, "loadtxt"),
        ("bse", petar.BSEType, "loadtxt"),
        ("bse_status", petar.BSEStatus, "fromfile"),
        ("bse.type_change", petar.BSETypeChange, "loadtxt"),
        ("bse.sn_kick", petar.BSEKick, "loadtxt"),
        ("bse.gw_kick", petar.BSEKick, "loadtxt"),
        ("bse.dynamic_merge", petar.BSEDynamicMerge, "loadtxt"),
        ("bse.binary_merge", petar.BSETypeChange, "loadtxt"),
    ]
    for suffix, cls, method in sse_bse_entries:
        path = case_dir / f"{output_prefix}.{suffix}"
        if path.exists():
            loader = lambda p, c=cls, m=method: getattr(c(), m)(str(p))
            check(suffix, path, loader, f"{cls.__name__}().{method}(path)")

    # Snapshots from snap.lst
    snap_list = case_dir / f"{output_prefix}.snap.lst"
    if snap_list.exists():
        last_snap = [l.strip() for l in snap_list.read_text().splitlines() if l.strip()][-1]
        offset = petar.HEADER_OFFSET_WITH_CM if external_mode != "none" else petar.HEADER_OFFSET

        check(f"snapshot:{last_snap}", case_dir / last_snap,
              lambda p: petar.Particle(interrupt_mode=interrupt_mode, external_mode=external_mode).fromfile(str(p), offset=offset),
              f"Particle({interrupt_mode}, {external_mode})")

        check(f"snapshot.single:{last_snap}.single", case_dir / f"{last_snap}.single",
              lambda p: petar.Particle(interrupt_mode=interrupt_mode, external_mode=external_mode).fromfile(str(p)),
              f"Particle({interrupt_mode}, {external_mode})")

        check(f"snapshot.binary:{last_snap}.binary", case_dir / f"{last_snap}.binary",
              lambda p: petar.Binary(member_particle_type=petar.Particle, interrupt_mode=interrupt_mode, external_mode=external_mode, G=petar.G_MSUN_PC_MYR).fromfile(str(p)),
              f"Binary({interrupt_mode}, {external_mode})")

    # Object snapshots
    check("object.1", case_dir / "object.1",
          lambda p: (lambda obj: (obj.addNewMember("time", np.array([], dtype=float)), obj.fromfile(str(p)))[1])(petar.Particle(interrupt_mode=interrupt_mode, external_mode=external_mode)),
          f"Particle(addNewMember('time')).fromfile('object.1')")

    # GroupInfo files
    for gp in sorted(case_dir.glob(f"{output_prefix}.group.n*")):
        m = re.search(r"\.n(\d+)$", gp.name)
        if m:
            n = int(m.group(1))
            check(f"group:{gp.name}", gp,
                  lambda p, n=n: petar.GroupInfo(N=n, interrupt_mode=interrupt_mode, external_mode=external_mode).fromfile(str(p)),
                  f"GroupInfo(N={n}, {interrupt_mode}, {external_mode})")

    # Profiles
    for pp in sorted(case_dir.glob(f"{output_prefix}.prof.rank.*")):
        # ... use pattern 8 matching strategy ...
        pass

    # DSM interrupt (if mode=dsm)
    if interrupt_mode == "dsm":
        check(f"{output_prefix}.interrupt", case_dir / f"{output_prefix}.interrupt",
              lambda p: petar.InterruptBinary(particle_type=petar.HardParticle, interrupt_mode="dsm").fromfile(str(p)),
              "InterruptBinary(HardParticle, dsm)", strict=False)
```

## Pattern 11: hard.debug H4 diagnostic logs (`data.*_h4_<N>_<rank>.log`)

Output of `petar.hard.debug` (hard-integrator diagnostic log). Three rules
separate this from every other file type:

1. **Reader: `petar.HermiteData` — never `sdar.HermiteData`.** The PeTar build
   adds columns (stellar evolution / external-field modes) the SDAR reader does
   not know; using it fails with column-index errors
   (`IndexError: index N is out of bounds for axis 1`). The constructor must
   match the build that wrote the log:
   ```python
   d = petar.HermiteData(N_particle=<N>, N_sd=<n_groups>, time_measure=False,
                         interrupt_mode='bse', external_mode='galpy')
   d.loadtxt('data..._h4_59_0.log', skiprows=0)
   ```
2. **Column count changes mid-file when AR groups form/dissolve** (each group
   adds slowdown/orbit columns). Detect and split FIRST, then read each segment
   with the matching `N_sd`:
   ```bash
   awk 'NF>1{c[NF]++} END{for(k in c) print k, c[k]}' data..._h4_59_0.log
   # e.g. {2369: 303651, 3083: 1, 4154: 1} -> split by NF into separate files
   ```
   Reading mixed-column lines with one constructor fails or corrupts data.
3. **`time` resets (`diff(time) < 0`) mark hard_debug restart rounds.** Analyze
   per round (split at resets): energy references are re-initialized each
   round, so a global cumulative `de` shows fake spikes at round boundaries. A
   few restart rounds after a group event are normal hard_debug behavior.

**Threading caveat:** with the default OpenMP thread count the log contains
only thread 0's particle subset — a *partial* view that silently invalidates
per-particle analysis (particles appear frozen or missing). Rerun with
`OMP_NUM_THREADS=1` (also ~30x faster under gdb; many threads plus gdb makes
the debugging session pathologically slow).
