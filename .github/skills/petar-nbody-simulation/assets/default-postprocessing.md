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
