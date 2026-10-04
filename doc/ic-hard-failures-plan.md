# Reproduction Plan: Primordial-Binary t=0 Hard Failure & Neighbor-List Overflow

Two pre-existing experiment-branch runtime failures, discovered while validating
the SEVN integration (2026-10-04) but **reproducible with the bse binary alone**
— they are NOT SEVN defects. This file is the cross-session reproduction ticket:
each recipe below was verified end to end on 2026-10-04 with mcluster seed 42;
follow it verbatim to regenerate the failures and their dump artifacts.

## Environment snapshot (of the verified reproduction)

- PeTar `experiment` branch (commit 55d09bd or later), configured
  `--with-interrupt=bse --enable-mpi --enable-omp --enable-avx512`; any
  equivalent bse-family binary works (the failures are family-independent —
  sevn reproduces both identically with its own flags).
- `mcluster` with `-s` seed support; `petar.init` (bse template).
- Runtime: `OMP_NUM_THREADS=4 OMP_STACKSIZE=128M`, 1 MPI rank.

## Common IC generation (deterministic, seed 42)

```bash
# A: primordial-binary IC (Issue 1)
mcluster -N 1000 -b 0.95 -m 0.7 -m 80 -C 5 -u 1 -s 42 >mc.log 2>&1
petar.init -s bse -v kms2pcmyr -f input_bse test.dat.10

# B: no-binary IC (Issue 2) — regenerates test.dat.10, run AFTER input_bse exists
mcluster -N 1000 -m 0.7 -m 80 -C 5 -u 1 -s 42 >mc_nb.log 2>&1
petar.init -s bse -v kms2pcmyr -f input_nb test.dat.10
```

Note `-m 0.7 -m 80` (IMF bounds) only makes the massive-star population
prominent; the failures also occur with other IMFs that populate ~75 M☉ stars.

## Issue 1: primordial-binary IC aborts at t=0 (two failure stages)

**Stage 1 — default options**: SDAR list overflow assert.
```bash
OMP_NUM_THREADS=4 OMP_STACKSIZE=128M petar.mpi.omp.avx512.bse \
    -u 1 -b 500 --bse-metallicity 0.00142857 -t 0.2 -o 0.2 -s 0.5 input_bse
```
Landmarks: `Assertion! ../SDAR/src/Common/list.h:236: (num_<nmax_) fail!`,
dump `data.hard_dump_t0.000000_M0_O0_c0_s*`; before the abort the log prints
`Large cluster, n_ptcl=1000 n_group=336` — the whole system forms ONE hard
cluster (chain-linked via huge r_search, up to ~1300 pc for binary members).

**Stage 2 — with `--hermite-n-neighbor-max 2000`**: survives the neighbor
overflow, then hits the hard-step safeguard instead.
```bash
OMP_NUM_THREADS=4 OMP_STACKSIZE=128M petar.mpi.omp.avx512.bse \
    -u 1 -b 500 --bse-metallicity 0.00142857 --hermite-n-neighbor-max 2000 \
    -t 0.2 -o 0.2 -s 0.5 input_bse
```
Landmarks: `Large H4-AR step cluster found (dump): step: 1000286`, dump
`data.dump_large_step_h4_n1000_g357_t0.000000_M0_O0_c0_s*` (>10⁶ AR steps at
the very first hard step). Tree step `-s 0.1` does not help (same dump).

**Evidence collected**:
- Identical failure (both stages) with `petar.mpi.omp.avx512.sevn` (gold-table
  flags) on the same IC file → physics-neutral to the SEVN integration.
- The stuck system contains a ~2.1 M☉ binary member with |v| ~ 875 km/s (tight
  primordial binary) inside the whole-cluster hard group; r_search chain
  (`r_in~0.4, r_out~4` pc for that member) links all particles.
- Related family: `hard_large_energy` dumps also appear in these runs — see
  `doc/hard-large-energy-se-mass-loss.md` (same investigation thread).

**Replay**: the dumped system replays with the family-matched
`petar.mpi.omp.avx512.bse.hard.debug` (ASan build, `setarch $(uname -m) -R`
prefix, full `data.par*` set in a scratch dir).

**Open questions / next steps**:
1. Why does group analysis chain the entire N=1000 system into one hard cluster
   (n_group=336 → single 1000-member group)? Inspect r_search/r_out growth for
   massive stars and binary members at t=0.
2. Is the slowdown hierarchy of tight primordial binaries inside such a large
   group handled (AR slowdown tree), or does the 10⁶-step count indicate the
   slowdown is bypassed for the whole-cluster group?
3. Compare against a production setup that historically worked (e.g. Pal5 runs
   with `-b 0.95`-type ICs): which parameter (eps, r_out limits, tree dt) differs?

## Issue 2: neighbor-address-list overflow (default 300)

```bash
OMP_NUM_THREADS=4 OMP_STACKSIZE=128M petar.mpi.omp.avx512.bse \
    -u 1 --bse-metallicity 0.00142857 -t 0.2 -o 0.2 -s 0.5 input_nb
```
Landmarks: `Error: neighbor address list full. max size=300, current size=300`,
raised in `../SDAR/src/Hermite/neighbor.h` (`checkAndAddNeighborSingle`); the
trigger particle is the most massive star (~75 M☉, id 548, `r_in~1.03 pc,
r_out~10.31 pc, r_search~13.84 pc`) whose search radius covers >300 neighbours
in the cluster core. Identical with the sevn binary.

**Workaround (verified)**: `--hermite-n-neighbor-max 2000` — the run then
completes (2 Myr and 12 Myr verified). Option definition: `n_neighbor_max` in
`src/hard.hpp` (IOParams name `hermite-n-neighbor-max`, default 300).

**Open questions / next steps**:
1. Should the default auto-scale (e.g. with N, or with the maximum stellar
   r_out / cluster half-mass ratio) instead of aborting?
2. Or should r_search growth for very massive stars be capped (it also drives
   the Issue-1 whole-cluster chaining — the two issues share this root)?
3. Document the flag requirement for massive-star clusters until resolved
   (done for the SEVN sample; the bse sample may need it too).

## Cross-references

- SEVN-side tracker: `doc/sevn_integration_plan.md` (§3 dependency note and
  the sample's known-issues comment).
- Hard-dump debugging procedure (replay rules, filename taxonomy):
  `.github/skills/petar-nbody-simulation/assets/script-tools.md` →
  "Hard dump debugging".
