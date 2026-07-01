# Lessons Learned — PeTar Agent Workflow

This file records mistakes, gotchas, and pitfalls discovered during PeTar agent sessions.
Each entry documents: what went wrong, why, and how to prevent recurrence.

**Lifecycle**: New entries are added by agents automatically after encountering issues.
Periodically reviewed → verified entries are elevated to `SKILL.md` as hard rules.

---

## Categories

- [Build & Configure](#build--configure)
- [Binary Selection](#binary-selection)
- [Simulation Execution](#simulation-execution)
- [Post-Processing](#post-processing)
- [Restart / Resume](#restart--resume)
- [Documentation](#documentation)

---

## Build & Configure

*(No entries yet)*

---

## Binary Selection

*(No entries yet)*

---

## Simulation Execution

### 2026-06-30: r-ratio 分析反了 r_in 的变化方向

**Mistake**: 用户请求 `--r-ratio 0.5` 与默认 0.1 对比，初步分析时直觉认为 r-ratio 增大 → r_in 减小（转换区变窄）→ 更少系统进入 SDAR。实际上 r_in = r_out × r_ratio，且 r_out 也会因 auto 调整而改变。最终 r_in 从 0.000838 pc 增大到 0.00532 pc（6.4×），**更多**系统进入 SDAR，含一个 4 体层次系统中的极紧双星导致 large step dump。

**Root cause**: 仅凭 r_ratio 数值直觉判断，没有先计算 r_in 绝对值再下结论。也没有注意到 r_out 会随 r_ratio 变化而 auto 调整。

**Prevention rule**: 分析 `--r-ratio` 或 `-r` 对 SDAR/Hermite 划分的影响时，必须先从 `data.par` 或运行输出中提取 r_out 和 r_ratio 的数值，**显式计算 r_in = r_out × r_ratio**，再对比 baseline 的 r_in 值。不要仅凭 r_ratio 比值做推断。

### 2026-06-30: Silent defaults for physics-defining IC parameters
**Mistake**: When user requested a Plummer N=500 cluster without specifying half-mass radius, virial ratio, IMF, mass range, or tidal filling, the agent silently used mcluster defaults (Rh=0.8 pc, Q=0.5, Kroupa IMF, 0.08–150 Msun) without asking — violating reproducibility and potentially producing unintended physics.
**Root cause**: The SKILL.md "Star-Cluster Generation Inputs" checklist was too brief (5 items) and did not enumerate all parameters that must be confirmed. The agent interpreted "do not proceed until the IC parameter set is complete" as satisfied by the 5 listed items, ignoring unlisted parameters.
**Prevention rule**: The definitive checklist in SKILL.md must enumerate every mcluster parameter that affects the initial physical state. Agents must iterate through all applicable rows and ask for any missing value. A "not safe to default" classification must be applied to Rh, Q, IMF, mass range, and COM phase-space.

---

### 2026-07-01: r-ratio 增大使总模拟时间增加而非减少

**Mistake**: 在 N=1000 无星族演化的 Plummer 星团中，增大 `--r-ratio` 从 0.1 到 0.5 使 `r_in` 扩大 2.5 倍（0.00335 → 0.00844）。直观上 dt_tree 可随之增大（s/r_in 从 0.29 到 0.46），每 Myr 步数减少到 1/4（1024 → 256）。但每步 wallclock 从 0.006 s 暴涨到 0.072 s（12 倍），导致总耗时反而增加约 3 倍（6.1 → 18.5 s/Myr）。

**Root cause**: `r-ratio` 增大使更多粒子落入 changeover 区域（r_in ~ r < r_out），这些粒子需要用 Hermite 直接积分和 group 搜索，每步计算量大幅增加。find.dt 的早期性能评估（data2 仅需 1.38 s/Myr）严重低估了实际开销，因为早期双星尚未形成、硬积分成本低。

**Prevention rule**: 评估 `--r-ratio` 对性能的影响时，不能只看步数减少。必须同时考虑：
1. `r_in = r_out × r_ratio` 的绝对值变化
2. 每步 wallclock 的变化趋势（从实际运行的 timing profile 提取）
3. find.dt 的早期评估可能严重低估长时间模拟的硬积分成本
4. 总耗时 ≈ (每 Myr 步数) × (每步 wallclock)，两者可能反向变化

---

## Post-Processing

*(No entries yet)*

---

## Restart / Resume

*(No entries yet)*

---

## Documentation

*(No entries yet)*
