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
