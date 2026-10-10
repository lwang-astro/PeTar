# Validation Report Design Rules

Rules for designing validation scenarios and their HTML reports. **T6
(`t6_pert_kernel_match_report.py`) is the reference implementation**; T1-T5 are
legacy style and will be migrated to these rules when next touched.

## 1. Report structure (paper style)

Numbered sections in this order:

- **§1 Background and purpose** — the design under test, why the question
  arises, and what the test establishes.
- **§2 Method**
  - *2.1 Quantities under comparison* — numbered equations (1), (2), ... for
    everything computed or compared.
  - *2.2 Metrics* — a table with columns *symbol | definition | interpretation*.
    Every physical quantity appearing anywhere in the report is defined here
    exactly once, with a symbol (χ, Σ, Φ, ε_E, ...), and later references use
    the symbol with its equation number.
  - *2.3 Runs* — a table *run | configuration | purpose*, plus the rationale
    for why structurally similar runs exist.
  - *2.4 Reproducibility* — software version (from the run log) and git commit;
    build commands; IC generation commands plus an IC parameter table (model,
    N, mass function, velocities, RNG and seed, units/constants, file format);
    the exact execution commands as recorded by the runner.
- **§3 Results** — one subsection per experiment; each table/figure is numbered
  and captioned; end with the automatic-checks table.

## 2. Tables and figures

- Captions must state **what is shown and how to read it**: which curve the
  check metric is taken from, where/why maxima arise, what bounded versus
  growing behaviour means.
- Figure encoding: **color distinguishes the model/run**, **line style the
  quantity type** (solid vs dashed); state the encoding in the caption. The
  legend sits in a dedicated strip *below* the plot whose height adapts to the
  entry count, so it can never overlap curves. Threshold lines are dotted gray,
  distinct from any curve dash pattern.
- Do not plot near-flat trends: a stationary quantity is reported as a number
  with a one-line statement, not a figure.
- Axes: state units, log/linear scale, and the meaning of the abscissa when it
  is not physical time (e.g. sample count).

## 3. Threshold semantics

- Distinguish **regression bands** (loose lower/upper bounds that detect
  structural breakage, e.g. a measurement going silent or bands decoupling)
  from **accuracy thresholds** (tight upper bounds on conservation metrics).
  State in the caption which kind each check is.
- Thresholds are properties of the IC: calibrate on the target IC and never
  copy numeric thresholds across ICs without re-measuring (the same mismatch
  metric measured 1.2e-3 on a production N=1000 IC and 0.24 on a compact
  engineered IC).

## 4. Test design principles

- **Guard that the tested path actually executes**: add checks that the
  configuration under test is present (e.g. tidal-tensor-active group =>
  artificial-particle count >= expected; single tree step => final time
  <= dt_soft). All-green on unexecuted paths is a silent failure.
- **Control chaotic noise by design, not statistics**: engineer the IC so the
  target configuration exists from the start (implanted triple, binary, ...)
  instead of relying on spontaneous rare events, which differ between
  realizations. When comparing kernels or code variants, run exactly **one
  tree step** so all runs replicate one trajectory; `-t` must equal dt_soft
  exactly (end times round UP to the step grid — `-t 0.004` runs two steps
  when dt_soft = 2^-8).
- **A/B pairs must be structurally identical** apart from the tested knob.
  Keep instrumented runs out of the A/B pair; prove measurement neutrality by
  the instrumented and clean runs ending bit-identically on the single step.
- **Parameter coupling**: `r_search_group <= r_in = r-ratio * r_out` is
  asserted; large group radii need a large `--r-ratio`. Prefer specifying
  `-r` only and letting dt_soft auto-derive; pin `-s` (dt_soft) only when the
  design demands it.
- **Portability**: build via `CXXFLAGS='<extra flags>' ./configure && make -j`
  (the Makefile appends `-D PROFILE` itself; bare `make` builds the configured
  family). The binary name is family-dependent — parse `TARGET` from the
  Makefile instead of hardcoding `build/petar.mpi.omp.avx2`.
- **Determinism**: fixed RNG seeds in the IC generator; plain family for
  plain-format ICs; 1 rank x 1 thread unless the test targets parallel paths.

## 5. Shared code

All generic helpers live in `report_common.py` (JSON loading, HTML escaping,
log/version parsing, command parsing, IC reading, metric cards, check rows,
and the plotting functions). Report scripts **must not redefine them**;
script-specific parsers stay in the script. Two plot styles:

- `svg_series_plot` — auto-palette per-series legend (legacy, T1-T5).
- `svg_curves_plot` — color = model, dash = quantity, bottom legend strip,
  optional translucent activity bands (current style; use for new reports).
