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
  Purpose: diagnostics and hard-integrator debugging
  Not a production simulation executable

- `*.format.transfer`
  Purpose: format conversion helper
  Not a production simulation executable

## Skill behavior

- If the binary suffix already determines the physics stack, the skill should not ask redundant enable/disable questions for the same feature.
- If the requested options conflict with the selected binary suffix, the skill should suggest the nearest matching solver binary.
