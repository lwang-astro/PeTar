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
- 用 petar.mpi.omp.avx512.bse.galpy 跑，给完整 workflow。
- 用已有 snapshot 做 restart，不要出现 petar.init。
- 用户误选 *.hard.debug 时，给出正确 solver 替代建议。
```

## 5) 准备提交清单

```text
请给出本次 skill 改动建议提交文件清单，并单独列出环境相关生成文件（可选提交）。
```
