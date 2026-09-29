# Changeover Radius and Timestep Tuning Reference

Detailed physics-selection guidance for `-r` / `-s` / `--r-group` overrides. The hard constraints
(parameter definitions, `--r-ratio` ↔ `r_in` coupling, `r_group < r_in`, do-not-confuse rule, and the
auto-detection chain) live in `SKILL.md` → "Changeover Radius and Tree Time Step". Read this file
before overriding the defaults or when the system is not a virialized star cluster.

## Switching criteria (`-s` ↔ `-r` coupling)

Controlled by `--dt-soft-sigma-factor`:

- **Freefall-based** (default, `--dt-soft-sigma-factor=0`): `dt_soft = P(r_in) / nstep`, where `P(r_in)`
  is the binary period at semi-major axis `r_in`. `nstep` is set by `--dt-soft-kepler-nstep`
  (default 16, recommended 64 for high accuracy). Mass-dependent via per-particle mass-weighted `r_in`.
  See Wang et al. 2026, ApJ, 998, 233.
- **σ-based** (`--dt-soft-sigma-factor > 0`, recommended 0.2): `dt_soft = alpha * r_in / (sqrt(3) * sigma_3D)`.
  Classic criterion from Iwasawa et al. 2015. When `-r=0`, `dt_soft` is first auto-estimated as
  `2.6e-4 * G*M / sigma_3D^3`, then `r_out` and `r_in` are derived.

## Criterion selection (Wang et al. 2026)

- Low-σ / loose clusters, wide binaries, subvirial/fractal ICs → **freefall-based** (σ hard to
  measure, or σ-based gives too-small r_in).
- High-σ / dense clusters (N ≥ 10⁴) → **σ-based** (larger r_in for same dt_soft).
- Both similar at optimal dt_soft for virial equilibrium N≈10³. Transition: `dt_soft / t_ch ∝ 1/N`.

## Environment considerations

- **Star cluster** (default): σ well-defined; both criteria applicable; auto chain works.
- **Isolated binary / few-body scattering**: σ undefined → auto defaults fail. Set `-r`, `-s`
  explicitly. Use the freefall-based criterion directly: choose `-s` and let the code derive `r_out`,
  or choose `-r` and set `--dt-soft-kepler-nstep` for accuracy. SDAR handles close encounters by
  default (`--r-group` auto from `r_in`). Size `r_in` (and hence the auto `r_group = 0.8·r_in`) to
  exceed the binary apocenter `a(1+e)` with margin — an `r_group` smaller than the apocenter means
  the binary is never captured as an SDAR group. Gate 5 is complementary, not conflicting: run
  `petar.find.dt` with the explicit `-r`/`--r-ratio` passed inside `-a`, then halve the recommended
  step; if the σ-based auto base scan step is unreliable, set it explicitly with
  `petar.find.dt -s <base>`.
- **Stellar disk / DSM**: rotation-dominated, σ inapplicable. Use the freefall-based criterion.
  Set `-r` explicitly from the disk scale height or the encounter region of interest. `--r-group`
  should be smaller than the disk scale height to avoid false SDAR triggers from shear.
- **Unknown / mixed**: fall back to the **freefall-based criterion** (no σ dependency, works for any
  particle configuration). Set `-r` and `-s` explicitly when possible; avoid relying on auto
  defaults unless the system is clearly a virialized cluster.

## Related parameters

- **`--dt-soft-kepler-nstep`** (default 16.0): steps per orbit for the freefall criterion.
  ns=16 matches pure Hermite final error; ns=64 matches peak error (recommended for precision).
  Active only with the freefall criterion.
- **`--dt-soft-sigma-factor`** (default 0.0): 0 = freefall criterion; > 0 = σ-based criterion
  with this alpha value.
- **`--r-search-min`** (default 0.0, auto): hard-cluster neighbor search radius.
  Auto: `max(search_vel_factor * sigma_1D * dt_soft + r_out, 1.2 * r_out)`. Override when the
  auto-detected radius misses a specific binary population.
- **`--r-search-group`** (default -1.0, auto): group candidate search radius. Auto: `1.0 * r_in`.
  Set to 0 to disable SDAR. Controls which particles are considered for SDAR group membership.
- **`--r-group`** (default -1.0, auto): multi-body detection radius and tidal tensor box size.
  Auto: `0.8 * r_search_group`. Must satisfy `r_group < r_in` for SDAR to activate on binaries.
  - **Accuracy trade-off**: SDAR (LogH) is designed and most accurate for **2-body** (binary)
    systems. Multi-body SDAR groups (hierarchical triples, 3+ bodies) have significantly degraded
    accuracy — the LogH method loses its Kepler-solver property when a third body is present
    (Wang 2025, ApJ, 978, 65). The hybrid BlogH method improves accuracy for weakly perturbed
    triples but is not yet integrated into PeTar's runtime SDAR.
  - **r_group too large**: risks capturing multiple stars into one SDAR group → multi-body SDAR
    (lower accuracy, especially for secular evolution like Kozai–Lidov). In disk/DSM environments,
    excessive r_group may trigger false SDAR groups from shear — keep r_group smaller than the
    disk scale height.
  - **r_group too small**: misses real binaries → they remain in the P3T Hermite/leapfrog
    integrator (also less accurate and slower for tight binaries).
  - **Default auto value** (0.8 * r_search_group ≈ 0.8 * r_in) is a reasonable balance for typical
    star clusters. Override only when a specific binary population is being missed (increase
    r_group) or too many multi-body groups degrade accuracy (decrease r_group).

## Worked examples: fix r_in and vary changeover width

```bash
# r_in = 0.01 pc, r_ratio = 0.5 → r_out = 0.02 pc (narrow)
-r 0.02 --r-ratio 0.5

# r_in = 0.01 pc, r_ratio = 0.1 → r_out = 0.10 pc (wide)
-r 0.10 --r-ratio 0.1
```

Set `-s` and `-r` explicitly to take full manual control. Any parameter left at default continues
to auto-derive.
