# TSV110 (HiSilicon Kunpeng-920) NEON benchmark and report

This directory contains the benchmark harness and the test report for the
128-bit NEON tree-force kernels of PeTar (`src/force_tsv110.hpp`), enabled by:

```shell
./configure --with-arch=tsv110
```

Contents:

- `REPORT.md` — full test report (in Chinese): hardware characterization,
  kernel/end-to-end benchmarks, figures and a point-by-point optimization list.
- `quad_newton/` — follow-up study (report + code + figures) on whether the
  quadrupole kernel should adopt the extra half Newton step used by the Fugaku
  kernel. Result: keep the cubic-only `rsqrt4` (no accuracy gain, ~20% slower).
- `figs/` — performance comparison figures used in the main report.
- `data/` — raw benchmark data (CSV/TXT) and the benchmark sources
  (micro-benchmarks, kernel benchmark driver, plotting script).

Post-report changes to `src/force_tsv110.hpp`:

- the EP-EP kernels now skip `mass <= 0` EPJ entries, matching the x86 SIMD and
  Fugaku kernels (measured cost ~1.000x, see `quad_newton/REPORT.md` §3.5);
- `NEON_QUAD_NEWTON` (default 0) selects the optional Fugaku-style extra half
  Newton step for the quadrupole kernels, with the cost/accuracy study in
  `quad_newton/`.

Key results on tsv110 (24 MPI ranks x 1 OpenMP thread, GCC 13.3, ~2.59 GHz):

| workload | baseline | NEON | speedup |
|---|---|---|---|
| kernel geomean (EP-EP / EP-SP quad / neighbor) | 1.00x | 3.97x / 4.65x / 2.35x | |
| demo 2000 stars, t=10 Myr | 74.4 s | 48.0 s | 1.55x |
| N=10000 single stars, t=2 Myr | 277.6 s | 123.5 s | 2.25x |

Numerical validation against the reference (`src/simd_test.cxx` style):
max relative force error ~2e-5 (tolerance 7e-3), neighbor counts identical.
