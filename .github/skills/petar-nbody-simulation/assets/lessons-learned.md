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
- [Agent Workflow & Delegation](#agent-workflow--delegation)
- [Configuration Hygiene](#configuration-hygiene)

---

## Build & Configure

### 2026-09-13: 头文件依赖缺口使 make 静默不重编——"改了没编"浪费整轮 gdb 推理

**Mistake**: 修 SDAR `symplectic_integrator.h` 的 merger 能量记账后跑 PeTar `make install`，输出全部 "up to date"，直接 install 了 `build/` 里的**旧二进制**；随后 gdb 行断点打印的账本值与源码逻辑矛盾，推演多轮（宏未定义？时序错位？）后才发现二进制根本不含新代码。同轮 SDAR 侧 `make -C sample/AR` 对 `information.h` 同样静默跳过。另两次犯同一 gdb 错误：**行断点停在语句执行前**——在第 N 行断点打印的是 N-1 及之前语句的效果，用它验证"清零是否生效"必然读到旧值（需断在块之后的行）。

**Root cause**: PeTar Makefile 对 `../SDAR/src/*.h` 与 SDAR sample/AR 对自身 src 的头依赖规则都不完整（未列出的头文件不在目标的依赖链上），mtime 变化不触发重编；install 目标无条件拷贝 build/ 内容，掩盖了未编译。gdb 行断点语义（语句前停止）与"验证赋值效果"的直觉冲突。

**Prevention rule**:
1. 改 SDAR/PeTar 头文件后，**不要信任 make 的 up-to-date 判断**：`rm` 目标二进制（或 `make -B <target>`）强制重编，跑前 `md5sum` 确认二进制 mtime/md5 已变；
2. gdb 验证一段赋值代码是否生效，断点设在块**结束后**的行；打印值与源码预期矛盾时，第一反应先确认二进制含新代码（`md5sum`/反汇编一行），再推演时序；
3. 长期修复方向：两个 Makefile 补全头依赖（wildcard `$(SDAR_SRC)/*.h` 入依赖表），或改用 compile_commands/化构建。

### 2026-09-17: "最新版安装后仍复现"实为多目标 install 非原子——复现前先核对已安装二进制 mtime

**Mistake**: 用户 12:18 `make install` 后报 `petar.hard.debug` 读 dump 断言"最新版可复现"。排查发现 12:18 install 只更新了主 `petar` 等目标，`petar.hard.debug` 目标 12:28 才重编落地；12:18–12:24 两次"复现"用的都是含旧代码（`hard_debug.cxx` 漏加 `time_offset`，见 Restart/Resume 2026-09-17 条）的陈旧工具二进制。主会话最初按"最新代码仍有 bug"推演窗口边界越界，白耗数轮——直接重跑用户命令却一次通过，才暴露真凶。

**Root cause**: `make install` 一次安装多个目标（petar、hard.debug、dump2test…），编译耗时长的目标落地晚于 symlink 创建时刻；"安装时间戳"≠"全部目标都已更新"。用户与 agent 都以"刚跑过 make install"为准据断言二进制内容。

**Prevention rule**:
1. 调试"确定可复现"的问题前，`ls -la --time-style=full-iso $(which <bin>)` 与 `readlink -f`，确认**每个用到的**已安装二进制 mtime ≥ 相关源文件 mtime；不符则先重装再谈复现；
2. "直接重跑一次用户命令"应是与源码推演并列的第一步（本例中它一次就否定了"最新版可复现"前提）；
3. 源码已修但未 commit 时（git status 出现 M），先看 diff——工作树修复可能已存在，问题只剩"哪个二进制含它"。

### 2026-09-17: 功能冒烟套件 `--phase all` 会原地重 configure 源码树并劫持 petar 符号链接

**Mistake**: 为验证一处修改跑 `python3 test/functional/run_functional_smoke.py --phase all`。该套件每个 case 在**源码树原地** `./configure <case 专用>` + make install（std→merger→dsm→bse→galpy→bse-galpy→bse-agama 依次覆盖），最后一个 case 把 `petar`/`petar.hard.debug` 符号链接与 Makefile/config.status 全部留在 bse-agama 配置——生产 dsm.galpy.gasdrag 状态被静默顶掉。且各 case 的 petar 步骤因 IC ASCII 列数不匹配（20 vs 25）全部 fail，与本修改无关。

**Root cause**: 套件设计为独立全量构建验证，不区分"源码树当前配置"是否是用户生产状态；无恢复步骤。

**Prevention rule**:
1. 验证单点修改优先用目标二进制直接跑（如 petar.hard.debug 读 dump），不要默认拉全套件；
2. 若必须跑套件：先记录 `./configure --help` 口径下的生产 configure 参数（可从二进制名反推：`petar.mpi.omp.avx512.64b.dsm.galpy.gasdrag` → `--with-interrupt=dsm --with-external=galpy --with-external-hard=gasdrag --enable-64b`），跑完立即重新 configure + make install 恢复；
3. 套件 run-phase 失败先看失败机制（启动即读文件格式错 vs 运行中断言），IC 列数类失败与运行时修改无关，不要据此回滚代码。

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

### 2026-09-13: SDAR integrator lessons moved

*(The four 2026-09-13 SDAR entries — getH evaluator mismatch, cross-build-flag step
comparisons, sorted-cck assumptions, ds-floor landing kills — moved to
`SDAR/.github/skills/sdar-fewbody-integration/assets/lessons-learned.md`.)*

### 2026-09-14: petar.init `-f` 参数顺序写反，覆盖原初条件文件

**Mistake**: 运行 `petar.init -v kms2pcmyr -s merger --radius 1e-5 -f N100_b0.5.txt input` 时，本意是将 `N100_b0.5.txt`（mcluster 输出）转换为 PeTar 输入文件 `input`，但 `-f` 指定的是**输出文件**，位置参数 `input` 才是输入文件。结果把 PeTar 格式的数据写回了 `N100_b0.5.txt`，覆盖了原始初始条件；后续重新跑模拟时只能从被污染的文件读取，导致粒子数错误（N=1）和一系列格式断言失败。

**Root cause**: 误以为 `-f` 是输入文件参数，与常见工具（如 `mcluster -o`）的语义混淆，也没有在运行后核对 `petar.init` 的提示语 `Transfer "..." to PeTar input data file "..."`。

**Prevention rule**:
1. `petar.init -f <output> <input>`：`-f` 永远是**输出**（PeTar 输入快照），位置参数是**输入**（原始粒子表）。
2. 执行后必须核对提示语 `Transfer "<input>" to PeTar input data file "<output>"` 是否符合预期。
3. 对原始 IC 文件做只读保护：用 `chmod -w <raw_ic>` 或保留一个 `.orig` 备份，防止误覆盖。

### 2026-09-17: mcluster 默认输出 `test.txt` 带表头，被 `petar.init` 吞成幽灵粒子

**Mistake**: 用 mcluster 默认输出生成 N=500 IC，`petar.init` 报 `Skip rows: 0`，把 `#` 表头行解析成零质量、位于原点的粒子——快照变成 N=501，平均质量被稀释。

**Root cause**: `-C` 决定输出格式：`-C 5` 给无表头的 `test.dat.10`，默认给带表头的 `test.txt`。SKILL.md 只记了 `-u`，`-C` 语义仅存在于 sample 脚本注释。

**Prevention rule**: 星团 IC 一律 `-C 5`；转换后核对 `head -1 input` 的第二个字段等于请求粒子数（实测 `-N 10` 默认 → `0 11 0`，`-C 5` → `0 10 0`）。`Skip rows: 0` 两种情况都打印，不是判据。

### 2026-09-17: `petar.find.dt` 漏 `-a "-u 1"` → 以 G=1 解读物理单位 IC

**Mistake**: 对 `-u 1` 快照运行 `petar.find.dt -i 1 input` 而 `-a` 未带 `-u 1`，自引力被放大 ~222 倍、集群假坍缩、group 爆炸（N_all=1520），最终以 `neighbor address list full` 崩溃——错误信息完全不指向单位。

**Root cause**: `petar.find.dt` 会真实启动 solver，`-a` 内一切均按生产参数解释；但 SKILL.md 的示例只举 `--galpy-set`/`-b` 等场景选项，把 `-u` 衬成可省项。

**Prevention rule**: `-a` 须复述单位模式与全部单位定义选项（`petar.find.dt -a "-u 1 <scenario-opts>" -i 1 input`）。遇 `neighbor address list full` 先核对 `-u`/`-G` 与 IC 量级（virial 平衡下 K ≈ |U|/2），再怀疑树/精度参数。

### 2026-09-17: 小 N 用多线程反而更慢——线程数须随 N 缩放

**Mistake**: N=500、t=2 Myr 给 8 线程，比单线程慢 66%（1.88 s vs 1.13 s；4 线程已慢 19%）；`petar.find.dt -o 8` 同样被并行开销拖长。

**Root cause**: 域分解/树构建/barrier 的每步开销与 N 无关，力的计算量才随 N 增长；N~500 时前者压倒后者。SKILL.md 只问“用几线程”而无缩放指引，双方只能“选个好看的数字”。

**Prevention rule**: 按 N 缩放（≲10³ → 1 线程；~10⁴ → 数条；≳10⁵ → 全部核心），对 `petar.find.dt -o` 与生产 launch 一视同仁；用户坚持小系统多线程时先给实测代价再确认。细节：`assets/script-tools.md` → "Parallel Launch Heuristics"。

---

## Post-Processing

### 2026-09-14: 读 PeTar 二进制快照漏掉 `offset=HEADER_OFFSET` 且未指定 `interrupt_mode`

**Mistake**: 用 `petar.Particle().fromfile('data.0')` 读取标准孤立星团快照时，第一个粒子被读成 header（mass=0, id=101），其余 100 个粒子的质量、ID 全部错乱（负质量、巨大 ID），导致总质量为负、Lagrangian 半径计算完全错误。同时产生 `Binary file size ... not aligned with dtype itemsize` 警告，但 initially 被当作非致命警告忽略。

**Root cause**: `petar.Particle.fromfile` 默认从文件字节 0 开始按粒子 dtype 解析，而 PeTar 二进制快照前 `HEADER_OFFSET`（24 字节）是文件头（time, N, file_id）；漏掉偏移会把 header 解释成第一个粒子。此外，不同编译特性（interrupt/external 模式）会改变粒子 dtype 的列布局，必须用对应的 `interrupt_mode`/`external_mode` 构造 reader。

**Prevention rule**:
1. 读取 PeTar 原始二进制快照时，**必须**使用 `petar.Particle(interrupt_mode='...', external_mode='...').fromfile(fname, offset=petar.HEADER_OFFSET)`。
2. `interrupt_mode` 取值为 `none` / `merger` / `bse` / `mobse` / `bseEmp` / `dsm`，必须与编译和运行时一致；`external_mode` 取值为 `none` / `galpy` / `agama`。
3. 任何 `Binary file size not aligned` 警告都应视为潜在 dtype/偏移/模式不匹配的信号，先修正 reader 参数再忽略警告。
4. 分析前做快速 sanity check：总质量应为正、粒子数等于 header.n、ID 在合理正整数范围内。

### 2026-09-14: hard.debug 日志用错 reader、未分割混合列、多线程残缺视图

**Mistake**: 分析 `petar.hard.debug` 的 h4 日志（n59 弹弓问题排查）时连环三个错：(1) 先用 `sdar.HermiteData` 读取，列匹配直接报 IndexError——PeTar 构建额外输出恒星演化/外场列，SDAR reader 不认识；(2) 换对 reader 后又整文件直读——AR group 形成使中途列数变化（2369→3083→4154），单一构造参数读不了；(3) 默认多线程下日志只含 thread 0 的粒子子集，前几轮“粒子冻结/缺失”的结论全部基于残缺视图，险些误导根因判断（实际那些粒子在其他线程里正常积分）。

**Root cause**: h4 日志的列布局同时依赖构建特性（interrupt/external 模式）与运行时组结构（SD 列随组增减）；hard_debug 的 OpenMP 输出只写 thread 0；time 重置标志重启轮、每轮能量参考重置。

**Prevention rule**: 读 `data.*_h4_*.log` 前必读 `data-readback-patterns.md` Pattern 11：(1) 只用 `petar.HermiteData`，构造参数匹配构建；(2) 先 `awk 'NF>1{c[NF]++}'` 列数直方图，按列数分割后分段读取；(3) 按 time 重置分轮分析；(4) 需要完整粒子视图时 `OMP_NUM_THREADS=1` 重跑（gdb 调试时也必须单线程）。

### 2026-09-17: `data.lagr` 的 `m` 是平均质量而非累计质量

**Mistake**: 把 `lagr.all.m` 当各 Lagrangian 半径内的累计质量使用，得出的 enclosed-mass 剖面自相矛盾，一度怀疑 `petar.data.process` 的粒子计数。

**Root cause**: `data.lagr` 列语义无任何文档。每块 6 列 = 质量分数 `[10%,30%,50%,70%,90%]` + 末列核心半径（Casertano & Hut 1985）；`m` 是平均质量、`n` 是粒子数（`tools/analysis/lagrangian.py`：`mcum[rindex]` 后再除以 `nlagr`）。

**Prevention rule**: 按 `data-readback-patterns.md` Pattern 1 读写；enclosed mass = `m × n`（shell 模式下 `m`/`n` 为壳层值）。下结论前先做一次 `m × n` 与独立求和的自洽检查。

---

## Restart / Resume

### 2026-09-16: 重启改 `-s` 后 par 文件中的派生参数过期，触发 DSM 559 断言崩溃

**Mistake**: tsalm2p3（DSM 构建）重启时用 `-s 0.0625` 缩小树时间步，未同步处理 `data.par.hard` 中存储的 `hermite-dt-max`（= 旧 dt_soft/2 = 0.0625，首启自动推导后被写入 par 文件）。重启后 drift 半步 `time_interrupt_max` = 0.03125，而 H4 块步边界可达 0.0625，DSM 活跃星的 `time_interrupt` 越过 timax 使 `next_dt<0`，在 `disk_star_merger.hpp:559`（`calcMassChange`）断言崩溃。同理 `-r` 改变后 `r-search-group`/`r-group`/`r-search-min` 等派生半径也不会自动更新，需手动加 `--r-search-min 0 --r-group -1 --r-search-group -1`，既繁琐又无文档。

**Root cause**: PeTar 每次启动把**解析后的最终值**（含自动推导结果）覆盖写入 data.par/data.par.hard；重启读取后，仅当值等于哨兵值（0/-1）才走自动推导路径——哨兵已被派生值替换，推导永不触发。命令行显式给出的选项与 par 文件载入的值无法区分，"改 `-s` 时哪些依赖参数需要联动"没有机制保证。

**Prevention rule**:
1. 2026-09-16 起源码已实现重启自动重算（`petar.hpp` `initialParameters()`），且 `-s`/`-r` 互耦（2026-09-17，与首启一致）：只给 `-s` → `r_out` 连同 `hermite-dt-max`、`r-search-min`、`r-search-group`、`r-group`、`hermite-acc-offset-sq` 重算；只给 `-r` → `dt_soft` 连同上述半径参数重算；两者都给则互不重置；`-s 0`/`-r 0` 哨兵是推导请求不触发交叉重置；`--r-ratio`/`--r-search-group-safety` 给出 → 组半径重算。当前命令行显式给出的参数不被重置；关闭特性的 0 值（如 `--r-search-group 0`）永不重置。半径变化时 `update_changeover_flag` 置位，重启逐粒子更新 changeover。重启只需 `petar -p data.par -s <new> [snap]`。
2. 显式 `--hermite-dt-max` 不得超过一个 tree drift step（KDKDK4 = dt_soft/2），否则会错过时间中断（恒星演化/DSM 事件）引发断言；越界时代码会打印警告。
3. 组合重启命令时牢记：PeTar 每次启动都会用解析后的参数**覆盖写** data.par*——复现实验之间不要互相复制 par 文件；重启二进制快照默认需 `-i 0`（`i=2` 默认按 ASCII 读会报 "cannot read header" 中止）。

### 2026-09-17: `petar.hard.debug` 读 merger dump 触发 559 断言——`time_interrupt_max` 忘加 `time_offset`

**Mistake**: 用 `petar.hard.debug` 读 DSM merger 触发的 `data.dump_interrupt_t*` 报 `disk_star_merger.hpp:559`（`time_interrupt>=time_record`）断言。HEAD 版 `hard_debug.cxx` 把 `interaction.time_interrupt_max = hard_dump.time_end`——而 dump 头里 `time_end` 是**相对步长**（如 0.0625），粒子 `time_record/time_interrupt/last_mass_change_time` 却是**物理时间**（含 offset，如 6321.0）。`calcMassChange` 收到 `next_dt = 0.0625 - 6321.06 ≈ -6321` → `time_interrupt < time_record` → 断言。另确认 dump 文件格式：`time_offset(F64)@0, time_end(F64)@8, gcm..., n_ptcl(S32)@72`，每粒子记录 192 字节（`time_record@80, time_interrupt@88, star.last_mass_change_time@120`）。

**Root cause**: 与 `-s` 重启同根：`time_interrupt_max` 的不变量是"≥ modifyOneParticle 可能收到的任何 `_time_end`"，但各调用点各自拼装（主运行 `stat.time+dt_drift`、hard_debug `time_end`），只保证窗口通常够用，不保证必然。hard_debug 拼装时漏加 offset 是直接死因。

**Prevention rule**:
1. 2026-09-17 起 `hard_debug.cxx` 已修为 `time_offset + time_end`，且 `ar_interaction.hpp::modifyOneParticle` 在 DSM/BSE 两分支统一用 `max(time_interrupt_max, _time_end)` 作有效上限——**任何** `_time_end > timax` 的调用点失配（陈旧 par、H4 块步越界、未来新调用点）都不再触发断言，只会把下次检查推迟到窗口末端。
2. 处理 hard dump 时间字段时记住语义：`time_offset`=物理时刻，`time_end`=相对步宽；粒子时间字段恒为物理时间。凡是与粒子时间比较/相加的量必须先加 offset。

### 2026-09-17: restart 漏掉 `<prefix>.par.*` 伴随文件与 `-i` 读取格式，两次启动即中止

**Mistake**: 只复制 `data.par` 与 `data.60` 重启 → `Cannot open file data.par.hard`；补上后又 `cannot read header` + core dump。两次都误判为“restart 本来就坏”，实为两个独立手续缺失。

**Root cause**: 参数被拆成主文件 `data.par` + 按特性拆分的伴随文件（所有构建都有 `data.par.hard`；BSE 系有 `.par.bse`/`.mobse`/`.bseEmp`；DSM/Galpy/Agama 同）；而快照默认写出二进制，`-i` 默认 2 只读 ASCII。canonical 形式两者都没写，照抄必失败。

**Prevention rule**: 复制 `<prefix>.par` **及全部 `<prefix>.par.*`**；读二进制快照用 `-i 0`（或 `-i 3`）。归属：`assets/input-source-workflows.md` → "Restart / resume"。修正后以“重启区间能量与原始运行逐位一致（ΔE = 0.000e+00）”作为手续齐全的判据。

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

---

## Agent Workflow & Delegation

### 2026-09-12: 大 subagent 拆分用于“已推导内容”的落地，延迟高且增量低

**Mistake**: 方案 G 的设计与验证约束已在主会话完成推导后，仍按 conductor 默认模式把“写实施计划”整体委派给 Planner（Pro 档），随后又把含全文底稿的文档整合委派给 Documentation Maintainer。两次调用用户体感异常缓慢；事后评估 Planner 增量价值 ≈ 10%（几个代码锚点核实与 .sh 命令摘录），Doc Maintainer 实为机械编辑执行（4 处 replace），且主会话事后仍自行 grep 复核，产生双重验证成本。

**Root cause**: 委派判断基于任务“形式”（写计划→Planner、改文档→Doc Maintainer）而非“上下文状态”。subagent 无状态：主会话已读过的 ~7000 行代码/文档（symplectic_integrator.h 4420 行、ar.cxx 717 行、plan doc 731 行、NOTES 680 行）必须由 agent 从磁盘重读；Pro 档模型推理慢；agent 内部 read→grep→read→write→verify 多轮串行，每轮都是完整模型推理。已有的“well-specified mechanical task 不升档”规则未被应用——当 prompt 已包含成品内容的绝大部分时，任务已经是 well-specified mechanical。

**Prevention rule**: 委派门槛基于上下文状态而非任务形式：
1. 上下文已在主会话建立、只剩“表达出来”（写计划/文档、聚焦单文件编辑）→ conductor 直接做，不经 subagent；
2. 只需少量新事实（锚点/命令/行号）→ 主会话 grep/read 直接获取；
3. subagent 保留给：大量**新**上下文探索（内置 Explore）、长仿真与回归战役（Simulation Engineer，已吸收 Build & Test Maintainer——隔离长输出有真实价值）、根因未明的开放式推理（才升 Pro）；
4. 委派前自问：这个任务的增量是“获取新上下文”还是“表达已有上下文”？后者不委派。

落地（2026-09-12）：agent 套件 9 → 3 + 内置 Explore。Planner/Researcher/Implementer/Build & Test Maintainer/Validation Analyst/Documentation Maintainer 移除；职责分别并入 Developer（规划/聚焦编辑/文档/lessons）、Simulation Engineer（configure/Makefile/smoke/validation harness 接线）、Reviewer（验证层级选择/T1-T3/阈值判定）。

---

## Configuration Hygiene

### 2026-09-17: 规则重复与矛盾——同一事实最多被复述 6 次，且存在 4 处互斥表述

**Mistake**: PeTar + SDAR 配置文件系统性健康检查发现：`petar.select` 规则出现 6 处、确认门 4 处、版本一致性检查 2 处（含逐字相同的 bash 块）、"用工具而非空谈" 4 处；必读清单在 `SKILL.md` 与 `assets/minimal-question-sets.md` 之间环形互指；`petar-simulation-engineer.agent.md` 在同一文件内把 assist 输出字段列了两遍。另有 4 处矛盾：SDAR 的 `G` 常量（已记入 SDAR lessons 2026-09-17）、PeTar "不得假设任何默认值" vs "默认 `-u 1`/前缀 `data`"、SKILL "现代运行无需 gather" vs agent "MPI 输出需先 gather"、`petar.select` 规则会在 SDAR-only 任务中被误用。常驻的 `AGENTS.md` 还要求每个会话读 357 行的 `HANDOFF.md`，单条典型执行路径达 ~1690 行。此外 4 个 asset（`binary-scenario-map`、`default-postprocessing`、`input-source-workflows`、`prompt-starters`）从未被 `SKILL.md`/`AGENTS.md` 引用，属死重。

**Root cause**: 规则按"写作时最方便的位置"落笔——首次出现处 + 每个新读者视角各复述一次——缺少"每个事实只有一个权威位置"的约束。`SKILL_CONTENT_INDEX.md` 的覆盖表只记录搬运历史，不校验重复；矛盾片段各自孤立看都正确，从未被并列对照。HANDOFF 等临时性文件被误列入常驻必读。

**Prevention rule**:
1. **每个事实一个权威位置**：执行正确性规则 → `SKILL.md`；内部实现/构建细节 → 仓库 `AGENTS.md`；参考型细节（签名、表、模板）→ `assets/`。其余出现处一律改为带路径的指针，不复述。
2. **常驻文件只放指针**：`AGENTS.md` 不复制 agent 文件、SKILL 规则或 HANDOFF 内容；HANDOFF/`prompt-starters.md` 属 on-demand 维护材料，不得进入常驻必读路径。
3. **同步新增时对照同一常量的既有取值**：把"X 在此处说 A、在彼处说 B"当作缺陷修复，而不是各自保留。
4. 分层口径见 `SKILL_CONTENT_INDEX.md` 的 "Layering Contract"；改动 `SKILL.md` 结构时同步更新该索引。
5. 注意"看起来像常识"的规则里，**反直觉的领域规则必须保留**（`find.dt` 结果减半、`commands.log`、求解器前确认门、`petar.init` 参数顺序、`--r-ratio` 与 `r_in` 同向），而通用的良好行为（"善用工具"、"简洁作答"、"不存在则创建目录"）保留一处即可。

### 2026-09-17: 优化设计目标未落盘，跨 session 后部分丢失

**Mistake**: 首次系统性健康检查优化时，用户给出的设计目标（① 减轻规则长度、保证上下文高效；② 层级化——顶层只留必读规则，场景细节下放子集文件；③ 职责清晰——`SKILL.md` 覆盖 PeTar 模拟与数据分析且自足，`AGENTS.md`/subagent 只管代码开发且不复述 SKILL），以及两条简化判据（主流模型必然会做的规则应删；重复/矛盾规则应合并，PeTar+SDAR 共有的下放 SDAR），只存在于那次会话的对话里，**没有写进任何仓库文件**。后续 session 再改 SKILL 时，落点与体积控制全凭临场判断：新增内容把上一轮的压缩成果退回约一半，并制造出同一事实最多 6 处复述。

**Root cause**: "设计目标"与"如何更新这些文件"本身**不属于任何既有文件的职责范围**——`SKILL.md` 是执行时加载的（写进去等于违反目标 ②），agent 定义是调用时加载的，而本应承载它的 `.github/skills/README.md` 当时**未被任何文件引用**（同时是 Hygiene 条点名的死重）。于是"维护这些文件"这一活动没有自己的规范。

**Prevention rule**:
1. 定制文件（`AGENTS.md`、`.github/agents/*.md`、`.github/skills/**`，两仓）的维护规范归 `.github/skills/README.md`：Design Goals / Simplification Questions / Layering Rules / Update Rules / Health-Check Procedure。两个 `AGENTS.md` 各留一行指针；改动前先读。
2. 改动的验收条件包含**体积预算**：报告 `SKILL.md`（及 `AGENTS.md`）的行/字节增减；把参考型表格或实测数据加进 `SKILL.md` 视为缺陷——移入 `assets/`，只留一行加指针。
3. 新增规则后立即 grep 它，其余出现处一律改为 `<path> → "<section>"` 指针。
4. 设计目标类信息一旦确立即随手落盘到上述规范文件，不要依赖会话记忆。
