# Prompt Starters（新机器续接用）

以下提示词用于在新电脑的 VS Code Chat 中快速进入同一开发状态。

## 1) 完整恢复上下文

```text
请先阅读 .github/skills/petar-nbody-simulation/HANDOFF-2026-03-27.md。
然后执行：
1) .github/skills/petar-nbody-simulation/assets/update_option_inventory.sh
2) .github/skills/petar-nbody-simulation/assets/update_script_tool_help.sh
最后总结当前机器上的 solver/helper 划分和可用脚本工具。
```

## 2) 仅更新二进制能力

```text
请刷新并检查 .github/skills/petar-nbody-simulation/assets/option-matrix.md，
并把 SKILL.md 里的 solver 列表更新为与本机一致。
```

## 3) 仅更新脚本工具能力

```text
请依据 Makefile.in 里的 install_script_tool 和本机安装情况，
刷新 script-tools.md 与 tool-help 目录，并同步更新 SKILL.md 的工具章节。
```

## 4) 回归测试 skill 行为

```text
请对以下请求做“预期输出检查”，不需要实际跑模拟：
- 用本机可用的 bse+galpy solver（例如 petar.mpi.omp.avx2.bse.galpy）跑，给完整 workflow。
- 用已有 snapshot 做 restart，不要出现 petar.init。
- 用户误选 *.hard.debug 时，给出正确 solver 替代建议。
- 对 `petar.movie` / `petar.data.process` 人为给错 `-i/-t/-s/--snapshot-type` 时，应先识别读取错位警告，停止当前处理，再改参数重试。
- 对 `petar.get.object.snap` / `petar.galev.process` / `petar.format.transfer.post` 人为给错 `-s` 或快照类型相关参数时，也应按同样机制处理：先识别错位警告或 `utf-8` 错误，停止当前处理，再修正参数；若仍存在，则明确说明输出不可信。
- 对 `petar.data.gether`，默认不要自动加 `-g`；只有用户明确要求整合 group 文件时，才使用 `petar.data.gether -g <prefix>`，并提醒 group 合并文件可能较大。
```

## 6) 聊天模板回归

```text
请用三段式回复模板处理以下需求（先需求筛选，再性能推荐，最后特殊模块可选建议）：
- 需求: bse + galpy, 不需要 gasdrag, 不需要 pn
- 再给默认候选 solver 和命令
- 最后额外说明何时才建议 --enable-64b / --enable-mpfrc / --enable-gperf / --with-debug
```

## 7) 最小提问清单（开场）

```text
请先用“最小提问清单”收集缺失信息，并且只问缺项：
1) interrupt: off/base/bse/mobse/bseEmp
2) external(long-timescale): off/galpy/agama
3) external-hard(short-timescale): off/gasdrag
4) pn: off/pnhermite/pnsdar/pnall
5) 是否需要 GPU；若不需要则给最佳 CPU 候选
6) 核心运行参数: -u, -t, -o, 初始数据来源(raw/snapshot/restart)
7) 并行规模信息: 粒子数级别(10^3/10^4/更大) 与 初始双星比例(无/少/多)
8) 若 N~10^3 且双星无/少，默认 OMP_NUM_THREADS=1 且 mpiexec -n 1；若双星多可考虑多线程
9) 输出前缀(-f): 若用户未指定且非重启且当前目录无既有 data* 输出，则默认用 data；否则再改为 data.<tag>
若信息已足够，请直接给候选 binary 与可运行命令。
```

## 5) 准备提交清单

```text
请给出本次 skill 改动建议提交文件清单，并单独列出环境相关生成文件（可选提交）。
```

## 8) Python 分析代码直出（基于 data_analysis.ipynb）

```text
请参考 sample/data_analysis.ipynb，直接给我可运行 Python 代码（不是 shell 命令）来做以下分析：
- 读取 data.10 与 data.10.binary，并绘制 semi-ecc 分布
- 读取 data.lagr 并绘制 0.1/0.3/0.5/0.7/0.9 及 core 半径随时间演化
- 若是 BSE 场景，额外绘制 HR 图（温度-光度，按 stellar type 着色）
要求：
1) 明确写出 interrupt_mode/external_mode/G/snapshot_format 的匹配参数
2) 原始 snapshot 使用正确 header offset
3) 输出一段可直接复制运行的完整代码
```

## 9) Python 事件与对象时序分析

```text
请给我 Python 代码，读取并分析以下文件（基于 sample/data_analysis.ipynb 模式）：
- data.sse.type_change, data.sse.sn_kick
- data.bse.type_change, data.bse.sn_kick, data.bse.gw_kick, data.bse.dynamic_merge
- object.1 或 object.MS.MS.binary（来自 petar.get.object.snap）
任务：
1) 打印前 10 条关键事件表
2) 用 BSEMerge 汇总 merger candidate
3) 画 object.1 的 x-y 轨迹（颜色为 time）或 BMS 的 semi-time 散点图
并说明每一步读取时需要的类与关键参数。
```
