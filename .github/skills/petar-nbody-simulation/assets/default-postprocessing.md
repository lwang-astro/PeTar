# PeTar Default Post-processing

Use these defaults when the user asks for a full workflow rather than only the launch command.

## Isolated cluster

```bash
petar.data.process -G 0.00449830997959438 data.snap.lst
```

## BSE / SSE

```bash
petar.data.process -i bse data.snap.lst
```

## Galpy

```bash
petar.data.process -t galpy --r-escape tidal data.snap.lst
```

## BSE + Galpy

```bash
petar.data.process -i bse -t galpy --r-escape tidal data.snap.lst
```

## Agama

```bash
petar.data.process -t agama --r-escape tidal -G 0.00449830997959438 data.snap.lst
```

## BSE + Agama

```bash
petar.data.process -i bse -t agama --r-escape tidal -G 0.00449830997959438 data.snap.lst
```

## Rule

If the user explicitly asks for run command only, omit post-processing.

## petar.movie Defaults by Scenario

When the user asks for a full workflow that includes visualization, use these `petar.movie` argument templates based on the scenario. Mode flags `-i`/`-t` must match the producing solver. Add `--n-cpu 1` for deterministic smoke/functional runs.

### Isolated cluster, Merger, DSM

```bash
petar.movie -i <interrupt_mode> -t <external_mode> -m x-y -R 10 \
  -b -L data.lagr --rlagr-min 0 --rlagr-max 5 data.snap.lst
```

### BSE

```bash
petar.movie -i bse -t <external_mode> -m x-y -R 10 \
  -c logtemp -b -L data.lagr --rlagr-min 0 --rlagr-max 5 data.snap.lst
```

### Galpy / Agama

```bash
petar.movie -i <interrupt_mode> -t galpy -m x-y,x-y -R 10,10000 \
  --cm-mode core,none --marker-scale 1,0.1 \
  -L data.lagr --rlagr-min 0 --rlagr-max 5 data.snap.lst
```

### BSE + Galpy / BSE + Agama

```bash
petar.movie -i bse -t galpy -m x-y,x-y -R 10,10000 \
  --cm-mode core,none --marker-scale 1,0.1 \
  -c logtemp,logtemp -b \
  -L data.lagr --rlagr-min 0 --rlagr-max 5 data.snap.lst
```

## petar.movie flag reference

`petar.movie` generates movies from `data.snap.lst` (auto-generated at runtime). Mode-consistency
flags must match the producing solver family, or snapshot read errors occur
(SKILL.md → "Snapshot Read-Mismatch Policy").

### Required mode flags

| Flag | Value | When |
|------|-------|------|
| `-i` | `none`, `merger`, `bse`, `mobse`, `dsm` | Match solver interrupt mode |
| `-t` | `none`, `galpy`, `agama` | Match solver external mode |
| `-s` | `ascii`, `binary`, `npy` | `npy` for `petar.data.process` output; `binary` for raw solver output |
| `--snapshot-type` | `origin`, `post` | `post` for processed snapshots; `origin` for raw solver output |
| `-G` | float | Gravitational constant; 0.00449830997959438 for Msun/pc/Myr; 1.0 for Henon |

### Quick examples

```bash
# Particle distribution
petar.movie -m x-y -R 10 data.snap.lst

# HR diagram
petar.movie -H data.snap.lst

# Combined panels for BSE simulation
petar.movie -m x-y -R 10 -H -b -i bse -s npy data.snap.lst
```

### Snapshot format matching

- Raw solver output (`data.*`): use `-s binary --snapshot-type origin`
- Post-processed output (`data.*.single`, `data.*.binary`): use `-s npy --snapshot-type post`

For full flag reference and comparison-mode usage, see `README.md`.
