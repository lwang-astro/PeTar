# DSM Workflow Reference

Disk Star Merger (DSM) specific workflow. The hard rules are in `SKILL.md`; this file provides the detailed parameter reasoning and example walkthrough.

## Overview

DSM enables tidal disruption of stars by a central massive black hole or disk. The binary family must be compiled with `--with-interrupt=dsm`.

## IC Preparation

DSM requires an additional column in the raw IC: the particle type (star vs disk gas). The `petar.init` call must specify DSM mode:

```bash
petar.init -s dsm --type <type> --radius <radius> -v <vel_conversion> -f input <raw_ic>
```

| Flag | Meaning | Example |
|------|---------|---------|
| `-s dsm` | Enable DSM initialisation mode | Mandatory |
| `--type` | Particle type column format | `--type 2` (type in column 7) |
| `--radius` | Disk outer edge in simulation-length units (pc for `-u 1`) | `--radius 4.5092203040509496e-08` (0.2 au in pc) |

### Radius calculation

`--radius` sets the disk's outer edge in the simulation length unit. For `-u 1` (pc):
- 1 au = 4.84813681e-06 pc
- 0.2 au = 9.69627362e-07 pc
- 0.01 au = 4.84813681e-08 pc

Example: a disk with outer edge 0.2 au → `--radius 9.69627362e-07`.

## Runtime Parameters

DSM runs require two additional flags beyond standard settings:

```bash
petar -u 1 -t <t_end> -o <dt_out> -f <prefix> \
  --detect-interrupt 1 --dsm-new-star-mode 0 \
  -b <primordial_binary_count> input
```

| Flag | Purpose | Recommended |
|------|---------|-------------|
| `--detect-interrupt 1` | Enable DSM event detection | Required |
| `--dsm-new-star-mode 0` | DSM new-star creation mode; 0 = stable for smoke | Required for smoke |

## Post-Processing

DSM mode uses the same pipeline structure but with `-i dsm` flags:

```bash
petar.data.process --no-auto-resume -p <prefix> -i dsm -G 0.00449830997959438 data.snap.lst
petar.get.object.snap -i dsm -p object -f origin -m id 1 data.snap.lst
petar.format.transfer.post -d single -s binary -o npy -i dsm data.snap.last.lst
```

## DSM-Specific Outputs

### `data.interrupt`

This binary file records interrupted DSM events. Its layout depends on compile-time member ordering and may include trailing padding bytes. Validate by record count rather than strict dtype alignment:

```python
import petar
import numpy as np

interrupt = petar.InterruptBinary(
    particle_type=petar.HardParticle,
    interrupt_mode="dsm",
)
interrupt.fromfile("data.interrupt")
print("records:", int(interrupt.size), "ncols:", int(interrupt.ncols))
```

## Python Read Pattern

```python
import petar

# Particles (note: HEADER_OFFSET, not HEADER_OFFSET_WITH_CM)
particle = petar.Particle(interrupt_mode="dsm")
particle.fromfile("data.1", offset=petar.HEADER_OFFSET)

# Single/Binary particle snapshots
single = petar.Particle(interrupt_mode="dsm")
single.fromfile("data.1.single")
binary = petar.Binary(
    member_particle_type=petar.Particle,
    interrupt_mode="dsm",
    G=petar.G_MSUN_PC_MYR,
)
binary.fromfile("data.1.binary")
```

## Artifact Checklist

After a DSM run, expect:

| File | Contents |
|------|----------|
| `data.snap.lst` | Snapshot list |
| `data.N` | Raw snapshots (N = snapshot number) |
| `data.lagr` | Lagrangian radii |
| `data.core` | Core properties |
| `data.status` | Run status |
| `data.esc_single` | Escaped single particles |
| `data.esc_binary` | Escaped binaries |
| `data.interrupt` | DSM interrupted-binary records |
| `data.sse` / `data.bse` | Stellar evolution (if enabled) |
