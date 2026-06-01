# PeTar Skill 开发交接说明（2026-03-27）

本文件用于在新电脑/新 VS Code 会话中快速恢复当前 skill 开发状态。

## 快速发布说明（2026-04-13）

- 已完成 DSM 一步到位同步：`SKILL.md`、`minimal-question-sets.md`、`prompt-starters.md` 的 interrupt 枚举均包含 `dsm`。
- 已修复 `assets/update_option_inventory.sh` 的 interrupt 映射，防止资产刷新时丢失 `dsm`。
- 已重新生成 `assets/option-matrix.md`，`--with-interrupt` 映射为 `base | bse | mobse | bseEmp | dsm`。
- 一致性回归已通过：`make skill-check` 全部 PASS。

## 1) 当前完成状态

已完成：
- 建立并持续增强 PeTar skill 主文件。
- 支持按二进制后缀自动推断场景（isolated、bse、galpy、agama、组合场景）。
- 支持按二进制 `-h` 动态校验参数。
- 明确区分 solver 与 helper tools（`*.hard.debug`、`*.format.transfer`）。
- 支持按输入来源自动切换流程：raw 输入、现有 snapshot、restart。
- 支持场景化默认后处理。
- 将 `Makefile.in` 中 `install_script_tool` 的脚本工具纳入 skill 能力面。

## 2) 核心文件

主规则文件：
- `.github/skills/petar-nbody-simulation/SKILL.md`

能力资产文件：
- `.github/skills/petar-nbody-simulation/assets/option-matrix.md`
- `.github/skills/petar-nbody-simulation/assets/script-tools.md`
- `.github/skills/petar-nbody-simulation/assets/binary-scenario-map.md`
- `.github/skills/petar-nbody-simulation/assets/minimal-question-sets.md`
- `.github/skills/petar-nbody-simulation/assets/input-source-workflows.md`
- `.github/skills/petar-nbody-simulation/assets/default-postprocessing.md`

维护脚本：
- `.github/skills/petar-nbody-simulation/assets/update_option_inventory.sh`
- `.github/skills/petar-nbody-simulation/assets/update_script_tool_help.sh`

辅助帮助归档目录：
- `.github/skills/petar-nbody-simulation/assets/tool-help/`

## 3) 新电脑恢复步骤

在仓库根目录执行：

```bash
# 1. 进入仓库
cd /path/to/PeTar

# 2. 刷新二进制能力矩阵（按本机安装实际状态生成）
.github/skills/petar-nbody-simulation/assets/update_option_inventory.sh

# 3. 刷新脚本工具 help 归档
.github/skills/petar-nbody-simulation/assets/update_script_tool_help.sh
```

建议再做一次快速检查：

```bash
# 检查当前主机可见的主求解器家族
ls -1 /home/lwang/bin/petar.mpi.omp.avx512* 2>/dev/null || true

# 检查两份关键资产已更新
head -n 40 .github/skills/petar-nbody-simulation/assets/option-matrix.md
head -n 40 .github/skills/petar-nbody-simulation/assets/script-tools.md
```

## 4) 在新会话里如何“接上上下文”

建议在 VS Code Chat 首条消息直接粘贴：

```text
请读取 .github/skills/petar-nbody-simulation/HANDOFF.md，
并基于 .github/skills/petar-nbody-simulation/SKILL.md 继续维护。
先执行 assets/update_option_inventory.sh 与 assets/update_script_tool_help.sh，
然后告诉我当前机器上 solver/helper 划分与下一步建议。
```

## 5) 当前约定（重要）

- `*.hard.debug` 与 `*.format.transfer` 视为配套工具，不作为主模拟执行程序。
- 参数支持以所选二进制 `-h` 为准，不凭空猜测。
- 用户给定二进制时，优先二进制驱动场景推断，并仅追问最小缺失输入。

## 6) 提交建议

通常建议提交：
- `SKILL.md`
- 规则类资产文档
- 两个更新脚本

按需提交（环境相关生成物）：
- `option-matrix.md`
- `assets/tool-help/*`

若团队希望“仓库不包含主机快照”，可只提交规则与脚本，不提交生成物。

## 7) 下一步可继续做的增强

- 增加“提问模板到工具链”的触发样例库。
- 增加常见场景命令片段模板（isolated、bse、bse+galpy、restart）。
- 增加自动化一致性检查脚本（检查文档与本机能力矩阵是否冲突）。

## 8) 本次更新摘要（2026-04-13）

已完成 SKILL 与资产文档的一致性同步，重点如下：

- `SKILL.md`：
	- Sources 补全样例脚本：`star_cluster_plummer_N1k_binaries.sh`、`star_cluster_plummer_N1k_GalpyMWPot.sh`、`star_cluster_plummer_N1k_AgamaMWPotHunter24.sh`。
	- 工具能力面补充：`petar.galpy.help`、`petar.external.galpy`、`petar.external.agama`。
	- 新增“外势势场图生成链路”规则：`petar.external.<galpy|agama> -p input.par -m pot_conf` -> `petar.external.pot.movie` / `petar.movie --ext-pot`。
	- Agama 与 BSE+Agama 默认后处理更新为与 sample 对齐：
		- `petar.data.process -t agama --r-escape tidal -G 0.00449830997959438 data.snap.lst`
		- `petar.data.process -i bse -t agama --r-escape tidal -G 0.00449830997959438 data.snap.lst`

- `assets/script-tools.md`：
	- 增补上述外势/帮助工具说明。
	- 增加 external potential map workflow 小节。

- `assets/default-postprocessing.md`：
	- Agama、BSE+Agama 默认流程升级到与 `SKILL.md` 一致。

- `assets/prompt-starters.md`：
	- 新增“外势势场图与可视化链路”模板，便于续接对话时快速触发完整命令链。

## 9) 一键自检清单（更新后快速验一致）

在仓库根目录执行以下命令，用于确认 `SKILL.md` 与资产文档没有漂移：

```bash
.github/skills/petar-nbody-simulation/assets/check_skill_consistency.sh
```

若需要查看每类检查对应的原始 `grep` 命令，再使用下面明细：

```bash
# 1) SKILL 中是否包含新增样例来源
grep -nE 'star_cluster_plummer_N1k_binaries\.sh|star_cluster_plummer_N1k_GalpyMWPot\.sh|star_cluster_plummer_N1k_AgamaMWPotHunter24\.sh' \
	.github/skills/petar-nbody-simulation/SKILL.md

# 2) SKILL 中是否包含外势工具与外势 map 规则
grep -nE 'petar\.galpy\.help|petar\.external\.galpy|petar\.external\.agama|petar\.external\.<galpy\|agama>' \
	.github/skills/petar-nbody-simulation/SKILL.md

# 3) SKILL 与默认后处理中 Agama 后处理是否一致
grep -nE 'petar\.data\.process -t agama --r-escape tidal -G 0\.00449830997959438|petar\.data\.process -i bse -t agama --r-escape tidal -G 0\.00449830997959438' \
	.github/skills/petar-nbody-simulation/SKILL.md \
	.github/skills/petar-nbody-simulation/assets/default-postprocessing.md

# 4) script-tools 是否包含 external map workflow
grep -nE 'petar\.external\.galpy|petar\.external\.agama|pot_conf|petar\.external\.pot\.movie|petar\.movie --ext-pot' \
	.github/skills/petar-nbody-simulation/assets/script-tools.md

# 5) prompt-starters 是否包含外势可视化模板
grep -nE '外势势场图与可视化链路|petar\.external\.agama|petar\.external\.pot\.movie|petar\.movie --ext-pot' \
	.github/skills/petar-nbody-simulation/assets/prompt-starters.md
```

若任意命令没有匹配结果，优先检查最近一次对 `SKILL.md` 的修改是否已同步到对应 `assets/*.md`。

## 10) 本次补充更新（2026-04-13，DSM 一步到位收尾）

已完成 DSM 相关规则层同步与自动资产生成器修复：

- `SKILL.md`
	- interrupt 枚举补全：`base | bse | mobse | bseEmp | dsm`。
	- 增加 DSM 场景最小提问节（含 `--dsm-*` 关键参数范式）。
	- Binary-to-option family mapping 增加 `.dsm` 条目。
	- Quick checklist 的 interrupt 选项补全 `dsm`。

- `assets/minimal-question-sets.md`
	- 新增 `## DSM` 最小提问清单。

- `assets/prompt-starters.md`
	- “最小提问清单（开场）”中的 interrupt 候选补全 `dsm`。

- `assets/update_option_inventory.sh`
	- 修复 `--with-interrupt` 相关映射，确保生成结果包含 `dsm`：
		- `solver_has_interrupt` 支持 `dsm`
		- `token_group_from_configure_feature` 支持 `dsm`
		- configure 映射文案更新为 `base | bse | mobse | bseEmp | dsm`

- `assets/option-matrix.md`
	- 已重新生成，`Requirement-Driven Solver Filters` 中 `--with-interrupt` 映射包含 `dsm`。

回归结果：

```bash
cd /home/lwang/code/PeTar
./.github/skills/petar-nbody-simulation/assets/update_option_inventory.sh
make skill-check
```

- `make skill-check` 全部通过。

## 11) 本次补充更新（2026-04-13，验证框架 T1-T4）

已在 `test/validation` 下建立可扩展验证框架，用于代码大改后的快速回归。

新增/更新文件：

- 入口与指标
	- `test/validation/run_validation.py`
	- `test/validation/metrics.py`
- 初始条件生成
	- `test/validation/make_ic.py`
		- `--case t1`: 高偏心 2 体
		- `--case t2`: 层级 3 体（切换测试）
		- `--case t3`: 层级 3 体（稳定/压力）
		- `--case t4`: 层级 3 体 + 背景粒子（tree+hard）
- 场景定义
	- `test/validation/scenarios/t1_high_ecc_changeover.json`
	- `test/validation/scenarios/t2_hermite_sdar_switch.json`
	- `test/validation/scenarios/t3_hierarchical_triple.json`
	- `test/validation/scenarios/t4_tree_hard_from_triple.json`
- 配置与文档
	- `test/validation/criteria.json`
	- `test/validation/README.md`

当前验证状态（本机已完成）：

- Dry-run（T1-T4）通过：
	- `python3 test/validation/run_validation.py --dry-run`
- 实跑通过：
	- `test/out/report.t1.json`
	- `test/out/report.t2.json`
	- `test/out/report.t3.json`
	- `test/out/report.t4.json`

关键实现说明：

- `run_validation.py` 支持场景批量执行、日志解析、JSON 报告输出。
- 指标提取：`max_abs_error_over_total`、`max_abs_error_pp`、`regex_counts`。
- 判定类型：`max_threshold`、`regex_count_max`、`convergence_ratio`。
- T1 收敛比判定为条件启用：
	- 当 `petar_bin_order2 == petar_bin_order4` 时自动 skip。
	- 当 fine error 低于数值地板时自动 skip。

跨机器恢复建议（验证框架）：

```bash
cd /path/to/PeTar

# 1) 脚本语法检查
python3 -m py_compile test/validation/run_validation.py test/validation/make_ic.py test/validation/metrics.py

# 2) 全量 dry-run
python3 test/validation/run_validation.py --dry-run

# 3) 轻量 smoke
python3 test/validation/run_validation.py --scenario test/validation/scenarios/t2_hermite_sdar_switch.json --report test/out/report.t2.json
python3 test/validation/run_validation.py --scenario test/validation/scenarios/t3_hierarchical_triple.json --report test/out/report.t3.json

# 4) 全量
python3 test/validation/run_validation.py
```

如需强制执行 2/4 阶对照（T1：当前主测 kdkdk4，对照 kdk）：

```bash
python3 test/validation/run_validation.py \
	--var petar_bin_order2=<kdk_binary> \
	--var petar_bin_order4=<kdkdk4_binary> \
	--var petar_bin_switch=<switch_test_binary>
```

下一步待开发：

- 增加 blogh A/B 占位场景（待 blogh 接入后直接启用）。
- 让 runner 直接读取 `criteria.json`，减少场景文件中的阈值重复。
- 增加多 seed 统计回归模式（均值/方差/KS）。

补充（2026-04-13，跨机续开发后的最新增量）：

- `test/validation/run_validation.py` 已支持 `--criteria`（默认 `test/validation/criteria.json`），并优先从 criteria 解析阈值：
	- `max_threshold`: 使用 `criteria[scenario][metric][run_id]`
	- `regex_count_max`: 使用 `criteria[scenario]["regex_count_max"][regex_key]`
	- `convergence_ratio`: 支持 `criteria_key` 从 `criteria[scenario]["convergence"]` 读取
- 新增 blogh A/B 占位场景：
	- `test/validation/scenarios/t3_blogh_ab_placeholder.json`
	- 默认变量：`petar_bin_blogh_a`、`petar_bin_blogh_b`（默认跟随 `petar_bin_switch`）
- `test/validation/README.md` 已修复 markdown 代码块并更新新场景与用法。

## 12) 最短恢复提示词（跨机器一条消息）

在新电脑的 VS Code Chat 里，直接发送下面这段：

```text
请读取 .github/skills/petar-nbody-simulation/HANDOFF.md 并按其中“验证框架 T1-T4”继续开发。
先执行：
1) python3 -m py_compile test/validation/run_validation.py test/validation/make_ic.py test/validation/metrics.py
2) python3 test/validation/run_validation.py --dry-run
3) python3 test/validation/run_validation.py --scenario test/validation/scenarios/t2_hermite_sdar_switch.json --report test/out/report.t2.json
然后汇报当前机器可用的 petar 二进制家族，并继续实现下一步：blogh A/B 占位场景。
```

如需立即做 2/4 阶强制对照（不是默认 skip），使用：

```text
请读取 .github/skills/petar-nbody-simulation/HANDOFF.md。
先检查 ./configure -h 中 --with-step-mode（确认 kdk/kdkdk/kdkdk4 可用与默认值），再执行 T1 对照并指定二进制：
--var petar_bin_order2=<kdk_binary>
--var petar_bin_order4=<kdkdk4_binary>
--var petar_bin_switch=<switch_test_binary>
然后给出 report.t1.json 的判定摘要。
```

指定 kdk/kdkdk4（二进制优先）做 T1 强制对照（kdkdk 仅临时兜底）建议模板：

```text
请读取 .github/skills/petar-nbody-simulation/HANDOFF.md，并执行 T1 强制 2/4 阶对照。
请先执行 ./configure -h 并确认 --with-step-mode（kdk/kdkdk/kdkdk4）;
优先确认这台机器上的 kdk 与 kdkdk4 二进制名（若缺失请明确报错并给替代建议），再执行：
python3 test/validation/run_validation.py \
	--scenario test/validation/scenarios/t1_high_ecc_changeover.json \
	--var petar_bin_order2=<kdk_binary> \
	--var petar_bin_order4=<kdkdk4_binary> \
	--var petar_bin_switch=<switch_test_binary> \
	--report test/out/report.t1.json
最后输出：
1) 实际使用的 order2/order4/switch 二进制名
2) report.t1.json 的 PASS/FAIL 摘要
3) 若 kdkdk4 缺失，可用 kdkdk 二进制临时对照，并说明结果可比性风险
```

英文版最短恢复提示词（便于对外协作）：

```text
Please read .github/skills/petar-nbody-simulation/HANDOFF.md and continue from the “Validation framework T1-T4” section.
Run the following first:
1) python3 -m py_compile test/validation/run_validation.py test/validation/make_ic.py test/validation/metrics.py
2) python3 test/validation/run_validation.py --dry-run
3) python3 test/validation/run_validation.py --scenario test/validation/scenarios/t2_hermite_sdar_switch.json --report test/out/report.t2.json
Then summarize the available petar binary families on this machine and continue with the next planned task.
```

英文版 T1 强制对照模板（primary: kdk vs kdkdk4; kdkdk as temporary fallback）：

```text
Please read .github/skills/petar-nbody-simulation/HANDOFF.md and run the forced T1 comparison.
First run ./configure -h and check --with-step-mode (kdk/kdkdk/kdkdk4), then verify whether kdk and kdkdk4 binaries are available on this machine.
After that, run:
python3 test/validation/run_validation.py \
	--scenario test/validation/scenarios/t1_high_ecc_changeover.json \
	--var petar_bin_order2=<kdk_binary> \
	--var petar_bin_order4=<kdkdk4_binary> \
	--var petar_bin_switch=<switch_test_binary> \
	--report test/out/report.t1.json
If kdkdk4 is unavailable, you may temporarily fall back to a kdkdk binary, but you must state the comparability risk explicitly.
Finally report:
1) the actual order2/order4/switch binaries used
2) the PASS/FAIL summary from report.t1.json
3) whether skip_if_vars_equal was triggered
```
