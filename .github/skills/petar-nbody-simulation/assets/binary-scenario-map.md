# PeTar Binary Scenario Map

This note documents the default scenario inference used by the skill.

## Solver binaries

- `petar.mpi.omp.avx512`
  Scenario: isolated cluster
  Feature stack: base solver only

- `petar.mpi.omp.avx512.bse`
  Scenario: stellar evolution
  Feature stack: base + BSE/SSE

- `petar.mpi.omp.avx512.galpy`
  Scenario: external potential with Galpy
  Feature stack: base + Galpy

- `petar.mpi.omp.avx512.agama`
  Scenario: external potential with Agama
  Feature stack: base + Agama

- `petar.mpi.omp.avx512.bse.galpy`
  Scenario: stellar evolution with Galpy external potential
  Feature stack: base + BSE/SSE + Galpy

- `petar.mpi.omp.avx512.bse.agama`
  Scenario: stellar evolution with Agama external potential
  Feature stack: base + BSE/SSE + Agama

## Helper tools

- `*.hard.debug`
  Purpose: replay and analyse a **dump file produced by a run**. Not a debug build of the solver — it cannot start a simulation from a snapshot. See `script-tools.md` → "Hard dump debugging".

- `*.format.transfer`
  Purpose: solver-snapshot format conversion (including legacy `-c`/`-g` reads). Not a production simulation executable.

- `*.dump2test`
  Purpose: convert hard dump files into `petar.hard.test` input snapshots. Not a solver.

- `petar.hard.test`
  Purpose: standalone hard-integrator test driver on prepared input files. Not a solver.

## Standalone stellar-evolution tools

- `petar.bse` / `petar.mobse` / `petar.bseEmp`
  Valid executables only for the standalone-bse scenario (stellar evolution without N-body); not cluster solvers.

## Skill behavior

- If the binary suffix already determines the physics stack, the skill should not ask redundant enable/disable questions for the same feature.
- If the requested options conflict with the selected binary suffix, the skill should suggest the nearest matching solver binary.
