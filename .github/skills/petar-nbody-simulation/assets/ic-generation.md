# IC Generation Reference

Generator-based initial-condition workflow: availability, unit conversion, the complete star-cluster
parameter checklist, and `petar.init` invocation. The three always-applying counter-intuitive rules
(`mcluster -C 5`, `petar.init -f` argument order, unit confirmation) are owned by
`SKILL.md` → "IC and Unit Rules"; this file owns the workflow detail only.

## Generator availability

When generator-based star-cluster IC creation is required:

1. Prefer `mcluster_gpu` if available.
2. Else use `mcluster` if available.
3. If neither exists, stop and ask the user whether to proceed with installation guidance.

## mcluster to PeTar unit handling

- Output file: `-C 5` emits the headerless `test.dat.10` consumed below (full rule and
  N-verification: `SKILL.md` → "IC and Unit Rules").
- `mcluster` output unit depends on generator `-u`.
- Common case: `mcluster -u 1` gives mass `Msun`, position `pc`, velocity `km/s`.
- If target runtime is `petar -u 1`, convert velocity with:
  - `petar.init -v 1.022712165045695`
- Prefer normalizing units in `petar.init` unless the user explicitly requests a native-unit workflow.
- For any raw IC source, confirm or infer source mass/length/velocity units before `petar.init`.
- Recommended path: normalize the IC into the target PeTar unit system during `petar.init` so that,
  for `petar -u 1`, the final IC is in `Msun`, `pc`, `pc/Myr`.
- Advanced native-unit paths require a matching `-G` and consistent conversion of all unit-sensitive
  runtime and post-processing options; do not recommend this path unless the user explicitly asks
  to preserve native units.

## Star-Cluster Generation Inputs (Complete Checklist)

If IC must be generated for a star cluster (e.g., via mcluster), every parameter below must be
explicitly confirmed with the user. Do not assume defaults for any physics-defining parameter.
If the user does not volunteer a value, ask; do not proceed until all fields are resolved.

| # | Parameter | mcluster flag | Why it matters | Default if not asked |
|---|-----------|---------------|----------------|----------------------|
| 1 | Size scale | `-N` or `-M` | Total particle count or total mass | N/A — must be specified |
| 2 | Density profile model + parameters | `-P` | Plummer, King (W0), fractal, etc. | N/A — must be specified |
| 3 | Half-mass radius | `-R` | Controls density and dynamical timescale; directly affects evolution speed and escape rate | **Not safe to default** — ask user |
| 4 | Virial ratio (Q) | `-Q` | Q=0.5 = virial equilibrium; Q=0 = cold collapse; Q>0.5 = expanding. Drives initial dynamical phase | **Not safe to default** — ask user |
| 5 | IMF model + parameters | `-f` | Kroupa, Salpeter, top-heavy, etc. Determines stellar mass distribution | **Not safe to default** — ask user |
| 6 | Stellar mass range (min, max) | `-m` | Lower (typically 0.08 Msun) and upper mass limits | **Not safe to default** — confirm with user |
| 7 | Binary fraction + distribution model | `-b` | Fraction of binaries, period/semi-major axis distribution, mass-ratio distribution | Must be confirmed even for zero binaries |
| 8 | Mass segregation | `-S` | S=0 = no primordial segregation; S>0 = segregated | mcluster default (0) — confirm with user |
| 9 | Random seed | `-s` | Reproducibility. seed=0 = automatic (non-reproducible) | mcluster default (0) — confirm with user |
| 10 | External potential COM position | `-c` (via petar.init) | Cluster center-of-mass initial position (x,y,z) in simulation units | **Not safe to default** — ask user |
| 11 | External potential COM velocity | `-c` (via petar.init) | Cluster center-of-mass initial velocity (vx,vy,vz) in simulation units | **Not safe to default** — ask user |

**Enforcement rule**: this table is the single source of truth for star-cluster IC parameters —
`assets/minimal-question-sets.md` points here rather than restating it. Before composing any
mcluster command, iterate all applicable rows above; for each row without a user-provided value,
ask. Do not proceed to command generation until every applicable field is resolved.
