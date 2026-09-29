# Reconfigure and Rebuild Toolchain Notes

Read before reconfiguring or rebuilding the solver (Gate 4 version mismatch, feature-flag change,
or compiler/MPI stack switch). The scenario-flag derivation rules stay in `SKILL.md` → Gate 4;
this file owns the toolchain environment and verification procedure.

## Before running `./configure`

- Check `env | grep -E '^(CC|CXX|FC)='`. Compiler modules may export bare compiler paths
  (e.g., `CXX=/.../bin/g++`); these override autoconf detection and produce Makefiles whose `CXX`
  is a bare `g++` (no MPI wrapper) and whose `FCLIBS` may lose `-lgfortran`. If set, `unset CXX CC FC`
  in the shell (or pass `CXX=<mpicxx> FC=gfortran` explicitly to configure).
- When `./configure` runs on a cluster login node, auto-detected SIMD capabilities may differ from
  compute nodes; verify with explicit `--enable-avx2` / `--enable-avx512` flags.

## After configure, before building

- Verify the generated Makefile: `CXX=` must be an MPI wrapper (`mpicxx`/`mpic++`), and with stellar
  evolution enabled `FCLIBS` must contain `-lgfortran`. Symptom if skipped: undefined `MPI::` /
  `_gfortran_*` symbols at link.
- Re-running configure with different flags silently changes the generated `TARGET` names and
  physics flavor (e.g., dropping `--with-interrupt=bseEmp` downgrades to `bse`). Compare the new
  Makefile's `TARGET` line against the previously installed binary family before `make install`.

## Makefile hygiene

- All Makefiles are configure-generated (`AC_CONFIG_FILES`), **including the sub-Makefiles**
  (`bse-interface/`, `galpy-interface/`, `parallel-random/`). Never fix build problems by
  hand-editing them — the next configure run silently overwrites manual edits, and editor undo can
  regress them to stale content. Fix by re-running `./configure` with a corrected environment.

## Full clean rebuild after toolchain switch

After switching compiler/MPI module stacks: `make clean`, `make -C bse-interface clean`,
`make -C galpy-interface clean`, `make -C parallel-random clean`, and `rm -rf build`. Stale objects
from the old toolchain fail to link with misleading symbols (e.g., `MPI::` bindings that current
sources do not use).

## Toolchain-specific hazards

- Debug-family targets (`*.hard.debug`, `*.dump2test`) link ASan; the toolchain must provide a
  matching `libasan` at link and runtime (system gcc8 on this cluster lacks the `libasan` RPM; the
  gcc-14.2 module ships `libasan.so.8`).
- Linking site libraries built against a different MPI flavor (e.g., a system GSL needing
  `libmpi.so.40` while PeTar uses Intel MPI's `libmpi.so.12`) pulls two MPI runtimes into one
  process — treat the `ld` "may conflict" warning as a real hazard; prefer matching stacks or
  MPI-free library builds.
- Full incident detail and cluster-specific paths: `assets/lessons-learned.md`, "Build & Configure"
  entries of 2026-09-21/22.
