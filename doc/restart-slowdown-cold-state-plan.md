# Plan: restart leaves AR slowdown state cold — tight-binary grind after restart

Status: plan — open items only. Found 2026-10-19 by the functional `bse` case
(6.3 min vs seconds); investigation deferred.

## Symptom

Restarting from an intermediate snapshot makes the hard integration of
pre-existing tight binaries pathological, while the identical evolution
integrated continuously is cheap:

| | fresh run 0→10 | restart run (from data.3 → 10) |
|---|---|---|
| wall time | **1.4 s** | **5 m 38 s** (restart segment only) |
| `Warning: step count is signficiant large` (10⁶-step caps) | 0 | **20** |
| 121+124 Msun contact-binary merge trajectory | t≈3 → t≈6 | **identical** (t≈3 → t≈6) |
| final state / all functional assertions | pass | pass |

Physics, BSE evolution, and merge times are the same along both paths — only
the cost differs. The restarted process grinds the still-unmerged tight
binaries with repeated step-count-max caps, i.e. the AR slowdown machinery
that makes them cheap in continuous runs is not effective after restart.

## Ruled out (2026-10-19)

- Redundant restart CLI options (`-b`, `-o`, `-f`, `--write-group-info`
  re-passed): the minimal canonical `petar -p c.par -i 0 -t 10 c.3` is
  equally slow. (The functional runner now uses the minimal form regardless.)
- `petar.data.clear` interference: physics identical with/without it.
- BSE / stellar-evolution state reconstruction: merge trajectory identical.

## Working hypothesis

The AR slowdown / Hermite per-cluster state (slowdown timescale, slowdown
tree, per-member AR internals) is **not persisted in snapshots and not
re-derived at restart** — the continuous run keeps it warm, the restarted run
integrates pre-existing tight groups from a cold/naive state. Related area:
commit 5765150 ("reset Hermite state (acc/jerk/pot/time/dt) on cluster
copy-in") and the restart reconstruction path in `src/petar.hpp` (initial()
restart branch: "update group_data.cm to pcm data for search cluster after
restart").

## Reproduction (deterministic, ~6 min)

```bash
mkdir /tmp/bse_restart_grind && cd /tmp/bse_restart_grind
cp <repo>/test/out/functional_bse/ic.raw .
petar.select --require bse --optional mpi,omp,avx512,avx2
petar.init -s bse -v 1.022712165045695 -f input ic.raw

# A: fresh baseline — completes in ~1.4 s, zero step warnings
OMP_NUM_THREADS=1 OMP_STACKSIZE=128M petar -u 1 -t 10 -o 1 -f a -b 4 input &>a.log

# C: same IC, phase to 4, trim, restart — restart segment takes ~5.6 min, 20 warnings
OMP_NUM_THREADS=1 OMP_STACKSIZE=128M petar -u 1 -t 4 -o 1 -f c -b 4 input &>c1.log
petar.data.clear -t 3 c
OMP_NUM_THREADS=1 OMP_STACKSIZE=128M petar -p c.par -i 0 -t 10 c.3 &>c2.log
grep -c 'step count is signficiant' c2.log   # 20
```

Compare: `grep -c 'step count is signficiant' a.log` = 0; merge trajectory
identical (m>100 particle count per snapshot `a.N`/`c.N`).

## Open items (suggested order)

1. **Localize the grind**: diff the solver profile outputs (`a.prof.rank.0`
   vs `c.prof.rank.0` — AR step/AR_tsyn columns) per output interval to
   identify which clusters/windows hit the caps; correlate with the
   dump_large_step files the restart writes.
2. **Boundary probes** (cheap, isolates the trigger):
   - restart from `c.0` (t=0, binaries pristine): grind or not?
   - second restart from `c.6` (binaries merged): does speed recover?
   → grind tied specifically to *pre-existing tight groups crossing the
   restart* vs any restart.
3. **Inspect the restart reconstruction** in `src/petar.hpp` initial()
   restart branch + `src/hard.hpp`: what per-member AR/slowdown state is
   rebuilt from the snapshot `group_data`, what is left at defaults
   (slowdown timescale, slowdown tree, imf factors); compare with the
   in-run cluster-copy-in path that 5765150 fixed.
4. **Fix direction**: either persist the needed state in the snapshot format
   (bump/format-gated) or derive a warm start at restart from current binary
   properties (semi-major axis / period → slowdown timescale), avoiding the
   naive integration that walks into 10⁶-step caps.
5. **Verify**: restart-segment wall time returns to O(fresh-run); the
   functional `bse` case drops from ~6 min to seconds (it currently acts as
   the live detector for this issue; `data.dump_large_step_*` artifacts from
   the reproduction are regenerable via the recipe above).

## Cross-references

- Functional detector: `test/functional` `bse` case restart segment
  (universal restart flow, aggressive `clear_t=3`).
- Lessons-learned 2026-10-19 entry (same finding, prevention rules).
- Related src areas: `src/petar.hpp` restart branch, `src/hard.hpp` cluster
  copy-in / slowdown handling, commit 5765150.
