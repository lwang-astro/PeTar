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

Two keyword arguments control column layout and must match the producing solver:

| Argument | Values | Purpose |
|----------|--------|---------|
| `interrupt_mode` | `none`, `merger`, `base`, `bse`, `mobse`, `dsm` | Must match solver `--with-interrupt` |
| `external_mode` | `none`, `galpy`, `agama` | Must match solver `--with-external` |

Pass these consistently to all readers. See `Snapshot Read-Mismatch Policy` in `SKILL.md` for how to diagnose mismatches.

## Pattern 1: Lagrangian Radii (`data.lagr`)

```python
lagr = petar.LagrangianMultiple(external_mode="<mode>")
lagr.fromfile("data.lagr")
```

- `external_mode` controls whether COM-offset columns are present.

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

The snapshot header offset depends on whether an external potential was active:

```python
offset = petar.HEADER_OFFSET_WITH_CM if external_mode != "none" else petar.HEADER_OFFSET
```

### Reading particles

```python
particle = petar.Particle(interrupt_mode="<mode>", external_mode="<mode>")
particle.fromfile("data.1", offset=offset)
```

### Reading single/binary component snapshots

After `petar.data.process`, each snapshot has split files:

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

## Pattern 7: Object Snapshots (`object.<N>`)

Object snapshots produced by `petar.get.object.snap` require pre-registering a `time` member before reading, because the object time-series format includes a time column not in the standard Particle layout:

```python
obj = petar.Particle(interrupt_mode="<mode>", external_mode="<mode>")
obj.addNewMember("time", np.array([], dtype=float))
obj.fromfile("object.1")
```

This pattern is required for all interrupt/external mode combinations.

## Pattern 8: Profile Tables (`data.prof.rank.*`)

Profiles are ASCII tables whose column count depends on the compile-time configuration (GPU, FDPS version, old_version flags). The correct `Profile()` keyword arguments must be determined by matching against the actual column count in the file.

### Matching strategy

Try the following `Profile()` constructor keyword arguments in order until `ncols` matches the file's column count:

```python
candidates = [
    {},                                                    # default
    {"use_gpu": True},
    {"old_version": True},
    {"use_gpu": True, "old_version": True},
    {"FDPS_version": 7.0},
    {"use_gpu": True, "FDPS_version": 7.0},
    {"old_version": True, "FDPS_version": 7.0},
    {"use_gpu": True, "old_version": True, "FDPS_version": 7.0},
]
```

For each candidate, read the first data row to check `ncols`:

```python
first_row = np.loadtxt(path, skiprows=1, max_rows=1, ndmin=2)
ncols = first_row.shape[1]

for kwargs in candidates:
    prof = petar.Profile(**kwargs)
    if prof.ncols == ncols:
        prof.loadtxt(str(path), skiprows=1)
        break
```

### Column-count mismatch fallback

If a profile has additional diagnostic columns in later rows (after a consistent numeric prefix), parse the minimum common column count:

```python
min_cols = 0
with path.open("r") as fh:
    for i, line in enumerate(fh):
        if i == 0:
            continue
        tokens = line.strip().split()
        if not tokens:
            continue
        row_cols = len(tokens)
        min_cols = row_cols if min_cols == 0 else min(min_cols, row_cols)

if min_cols > 0:
    np.loadtxt(str(path), skiprows=1, usecols=tuple(range(min_cols)))
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
interrupt.fromfile("data.interrupt")

# Check: records > 0 is a valid read; dtype warnings are expected and non-blocking
```

## Warning Classification

When running read checks, warnings fall into two categories:

### Blocking (treat as read failure)

- `Binary file size is not aligned with dtype itemsize` — format/schema mismatch
- `not aligned with dtype itemsize` — binary layout mismatch
- `File may be truncated` — incomplete snapshot
- `dtype mismatch` — column schema mismatch

These indicate a mode-flag or format mismatch. Stop, correct flags, retry.

### Non-blocking (record but do not fail)

- `ffmpeg` macro-block resize warnings from `petar.movie`
- Trailing-byte warnings on `data.interrupt` in DSM mode (`strict_mismatch=False`)

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
    for suffix, cls, kwargs in [
        ("sse", petar.SSEType, {}),
        ("sse.type_change", petar.SSETypeChange, {}),
        ("sse.sn_kick", petar.SSESNKick, {}),
        ("bse", petar.BSEType, {}),
        ("bse_status", petar.BSEStatus, {}),
        ("bse.type_change", petar.BSETypeChange, {}),
        ("bse.sn_kick", petar.BSEKick, {}),
    ]:
        path = case_dir / f"{output_prefix}.{suffix}"
        if path.exists():
            loader = lambda p, c=cls, k=kwargs: c().loadtxt(str(p)) if hasattr(c(), 'loadtxt') else c().fromfile(str(p))
            check(suffix, path, loader, f"{cls.__name__}()")

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
