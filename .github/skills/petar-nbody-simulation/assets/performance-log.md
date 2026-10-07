# Performance Log (measured run configurations)

Single authoritative home for **measured** solver performance: per-step cost, per-Myr pace, and
launch-configuration comparisons from real runs. Consult this file **before** proposing a launch
configuration (threads/ranks) or a runtime estimate — do not guess from intuition or short-probe
wall clock. Append one row per production/validation run (or configuration probe) after it finishes.

## How to measure (rules)

- **Per-step cost**: solver's own profile output — `Wallclock time per step (local): Total` at the
  end of the run log, or `<prefix>.prof.rank.<n>` column `Total` (5th numeric column; header is
  wrapped over several lines). Divide by the tree step to get s/Myr: `Total/step / dt_soft`.
- **Pace check**: mtime deltas of committed output snapshots (`stat -c "%y %n" <prefix>.<t>*`)
  — independent of startup/IO and works post-hoc on any run.
- **Never** compare configurations by short-probe wall clock: startup (MPI init, IC read) and
  output IO dominate probes shorter than a few Myr and invert rankings (measured 2026-10-07:
  a 0.25 Myr probe wall was 4× its compute time).
- Record the environment: missing `UCX_VFS_ENABLE=n` inflated per-step wall ~8× on WSL/UCX
  (2026-10-07) by adding per-communication overhead — every launch must export it.

## Measured records

| Date | N | Binary family | Machine / env | Launch config | ms/step (prof) | s/Myr | Source | Notes |
|------|---|---------------|---------------|---------------|----------------|-------|--------|-------|
| 2026-09-19 | 1000 | plain (mpi.omp.avx2) | WSL DESKTOP-AMD | unknown (predates log) | ~7.4 (derived) | 1.94 | `data.*` mtime span 47→200 | original 200 Myr run, baseline pace |
| 2026-10-07 | 1000 | plain, restart @47 | WSL DESKTOP-AMD, `UCX_VFS_ENABLE=n` | direct, 1 OMP thread | 10.03 | 2.57 | `bench_direct1` prof | probe 47.00→47.25; recommended config |
| 2026-10-07 | 1000 | plain, restart @47 | same | `mpirun -n 1` / `-n 2` | 10.58 / 9.19 | 2.71 / 2.35 | bench prof | no meaningful gain at N=10³ |
| 2026-10-07 | 1000 | plain, restart @47 | same | `mpirun -n 4` / `-n 6` | 9.78 / 12.70 | 2.50 / 3.25 | bench prof | 6 ranks already worse; comm overhead grows |
| 2026-10-07 | 1000 | plain, restart @47 | same, **no** `UCX_VFS_ENABLE=n` | direct, 1 OMP thread | (wall ≈ 660) | ≈174 wall | probe wall | UCX per-step overhead ~8×; wall-based numbers unusable |
| 2026-10-07 | 1000 | plain, restart @47 | same | direct, OMP default (15 threads) | — | ≈436 wall | probe wall | 1500% CPU spin-wait; **slower than 1 thread** |
| 2026-10-07 | 1000 | plain, restart @47 | same | direct, 1 OMP thread (production) | — | 2.5 steady | `fix47.*` mtimes | full fix-verification 47→200 (restart at 167 past a transient NaN assert); |Σmv| conserved 47–66 vs original 55→238 |
| 2026-09-2x | 500 | plain | WSL DESKTOP-AMD | 1 / 4 / 8 OMP threads | 1.13 / 1.35 / 1.88 s per 2 Myr | — | prior measurement (moved from script-tools.md) | 4 threads +19%, 8 threads +66% vs 1 thread |
