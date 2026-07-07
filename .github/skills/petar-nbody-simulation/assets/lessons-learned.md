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

### 2026-07-07: SKILL.md Python 分析工具节缺少 Binary 类文档

**Mistake**: SKILL.md 的 "Python Data Analysis Tools" 节只记录了 `petar.Particle` 的用法，完全缺少 `petar.Binary` 类。当用户问如何读取双星数据（离心率、半长轴）时，agent 错误地用 `petar.Particle` + `offset=HEADER_OFFSET` 去读 `data.*.binary` 文件，并用错误的属性名（`.sma` 而非 `.semi`）。

**Root cause**: SKILL.md 只覆盖了原始快照读取这一种场景，没有涵盖 `petar.data.process` 后的 `.binary` 文件读取。agent 在没有查阅 `sample/data_analysis.ipynb` 或 `assets/data-readback-patterns.md` 的情况下凭记忆作答。

**Prevention rule**: SKILL.md 的 Python 分析工具节必须覆盖以下关键类：
1. `petar.Particle` — 原始快照（需指定 `interrupt_mode`/`external_mode` 和 `offset`）
2. `petar.Binary` — 后处理双星数据（`data.*.binary`，无 header，无需 offset）
   - 核心属性：`.semi`（半长轴）、`.ecc`（离心率）、`.p1`/`.p2`（分量星）
   - 单位转换：`-u 1` 下 `semi * 206265` → AU
3. agent 在回答 Python 读取相关问题时，必须先查阅 `assets/data-readback-patterns.md` 和 `sample/data_analysis.ipynb`，不要凭记忆推断 API。

### 2026-07-07: SKILL.md Python 节重构为映射表+硬规则

**Mistake**: SKILL.md 的 Python 分析工具节试图以内联教程的方式覆盖所有类的用法（只写了 Particle 和 Binary），但实际有十几种文件类型和对应的 reader 类（LagrangianMultiple, Core, Status, Escaper, SSEType/BSE*, GroupInfo, Profile 等），零散列举永远不可能完整。而且 agent 倾向于从 SKILL.md 的内联内容直接作答，不会主动查阅外部的 `data-readback-patterns.md`。

**Root cause**: 结构设计错误 — 教程式内容放在 SKILL.md 本身会产生"这里有完整信息"的错觉，抑制了 agent 查阅权威参考文件的意愿。

**Prevention rule**: SKILL.md 的 Python 节不应尝试做教程。正确结构是：
1. **硬性规则**（红色警语）：写任何分析代码前 MUST READ `assets/data-readback-patterns.md`
2. **快速映射表**：文件类型 → Reader 类 + 一行关键提示（如 offset、kwargs 需求）
3. **通用关键字参数表**：`interrupt_mode` / `external_mode` 等跨类共用的参数
4. 详细的构造函数签名、属性列表、偏移量选择 → 全部放在 `data-readback-patterns.md`，禁止在 SKILL.md 中重复
