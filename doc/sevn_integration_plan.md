# SEVN 集成:现状、缺陷与开发计划

状态:2026-10-04 首次真实编译打通(bse 家族 + sevn 双模式均通过,含全部附属工具);
**运行未完成**——恒星初始化阶段受 SEVN 库-表匹配问题阻塞(见 §3)。
本文档是 SEVN 路线的权威缺陷/待办清单,完成一项划掉一项并注明验证方式。

## 1. 已确认缺陷(需开发)

### D1 随机数种子未接入 PeTar 并行随机数体系 —— 高优先级
- **2026-10-04 处置**:已删除误导性死代码——`value3_`/`rand3_` 结构定义与
  `BSEManager::initial()` 中的 `value3_.idum` 写入(SEVN 库不读取、BSE 家族
  也不用,kick 经 Fortran `rand_f64()` 走 PeTar 并行随机);`-idum` CLI 选项
  整体删除(原 experiment 无此选项,系 SEVN PR 引入且无任何消费者);
  README 新增 "Random numbers in stellar evolution" 一节。
  **仍待开发**:下方 rseed 暴露与接线。
- 现象:`BSEManager::initial()` 的 SEVN 分支把种子写入本地 `struct {int idum;} value3_`
  (bse_interface.h SEVN 声明块),但 **SEVN 库并不读取该结构**(它用自己的
  `utilities` 随机数);同时 `construct_default_sevn_params()` 把 `rseed` 列入
  `unneeded_options` 直接删除,用户无法给 SEVN 传种子。
- 后果:① SEVN 模式下 SN 踢等随机过程**不可复现**,且每次运行种子不受控;
  ② `rand_manager`(data.par.randseeds、每线程种子)在 SEVN 模式照常输出,
  形成"看起来有种子管理、实际没人消费"的假象;③ MPI 多 rank 下 SEVN 内部
  随机序列可能各 rank 相同(未验证)。
- 方案:在 `default_params_sevn` 暴露 `rseed`,由 `BSEManager::initial()` 用
  `rand_manager` 的每 rank/线程种子填充。

### D2 GW 并合反冲(gw_kick)未接入 SEVN 路径 —— 高优先级
- 现象:实验分支的 GW 反冲链(`GWKick::calcKickVel/calcFinalMass`、
  `getCompactChiRandom`、`compactOspinToChi`,bse_interface.h L1855-1890 共用区)
  只在 `evolveBinary()` 的 **非 SEVN 分支**(L2284-2335 的 merger_event 块)被调用;
  SEVN 分支走 `evolve_bco`(Peters 1964)后由 `merge_SEVN` 内部处理并合,
  **不产生** PeTar 侧的 GW kick 速度、并合产物质量比、event_flag=6 事件。
- 后果:SEVN 模式下 BBH/BNS 并合失去实验分支的自洽反冲物理与
  `fout_bse_gw_kick` 事件流;后处理 `BSEKick`/GW 统计缺失。
- 方案:在 SEVN 分支的 `merge_SEVN` 调用点之后,复用共用区的 gw_kick 链
  (SEVN 的 RemnantType/自旋字段可映射);或在 SEVN 库侧确认其内建 kick 处理
  并对齐事件记录格式。

### D3 `spin_3d` 语义在 SEVN 路径下不一致 —— 中优先级
- C++:`StarParameter::ospin[3]`(实验 3D 自旋)在 SEVN 下只写 `ospin[0]`
  (evolve_sevn.cpp 全部标量用法,`ospin[1..2]` 恒 0)。
- Python:`SEVNStarParameter` 却默认 `spin_3d=True`(3 列),`spin_3d=False`
  才是 1 列——与 bse 家族方向相反(bse 数据 3D 是新格式,SEVN 数据 1D 是唯一格式)。
- 后果:快照列宽解释混乱风险;BH 的无量纲 chi(3D)在 SEVN 数据中根本不存在。
- 方案:文档+类默认值明确"SEVN 恒为标量自旋",或将 `spin_3d` 对 SEVN 忽略并断言。

### D4 事件文件流在 SEVN 下的产出未验证 —— 中优先级
- `fout_bse_gw_kick/binary_merge/hyperbolic_tde/binary_tide/gw_tide_merge` 对
  所有 BSE_BASE 模式(含 SEVN)统一打开;但写入点大多在非 SEVN 代码路径
  (D2 的 GW 块、Binary_merge/Tide 记录)。SEVN 运行下这些流预计**只开不写**
  (空文件在退出时被 tmp 机制清掉,行为无害但能力缺失)。
- `type_change/sn_kick`(.sevn/.sevnB)由 `log_sevn_event` 驱动,结构待运行验证。
- 方案:随 D2 一并梳理;对照 BSE 事件全集(type_change、sn_kick、gw_kick、
  dynamic_merge、binary_merge、hyperbolic_tde、binary_tde、tide、gw_tide_merge)
  逐项标注 SEVN 的对应行为(原生/缺失/等价物)。

### D5 TDE 判据在 SEVN 下的字段假设未验证 —— 中优先级
- `getTDERadius()`(bse_interface.h L2421,共用区)依赖 `kw`(SSE 类型)、
  `mt/r/mco`;SEVN 经 PhaseBSE 映射提供 kw,但 `mco` 的语义(SEVN 的
  MCO_SEVN vs BSE 的 mc)**未核对**;ar_interaction 的 TDE 中断调用链
  (L1590-1612)在 SEVN 构建下的实际触发未测。
- 方案:构造 CO+星双星近遇用例,验证 SEVN 模式 TDE 事件能否产生与记录。

### D6 其余 BSE→SEVN 能力对照(初步排查,待运行期证实)
| BSE 家族特性 | SEVN 现状 |
|---|---|
| SN kick(random walk kick 速度+方向) | SEVN 原生(库内),事件经 `.sevnB.sn_kick`;方向随机性受 D1 种子问题影响 |
| GW 并合反冲 + 产物质量 | **缺失**(D2) |
| 双曲/双星 TDE | 判据共用但未验证(D5) |
| mass transfer / CE / common envelope | SEVN 原生(evolv2_SEVN/merge_SEVN),语义不同于 BSE,事件经 BEvent 列 |
| 动态并合(hyperbolic merge)记录 | 待验证(依赖事件 flag 链在 SEVN 分支的走通) |
| 3D 自旋/chi | **缺失**(D3,仅标量) |
| metallicity 表格范围 | 受限:MIST 0.7–150 M☉、Z 目录离散(无 0.001 精确值);PARSEC 默认表 ≥2.2 M☉ |
| 随机数复现 | **缺失**(D1) |

## 2. 工具链残留事项(低优先级)
- `petar.select` 的 `--require sevn` 已支持,但 sevn 家族二进制的实际选择流未测
  (需 sevn 版安装后回归一次)。
- `petar.bse.get.init.binary` 的 sevn 版(`petar.sevn.get.init.binary`)内容未按
  SEVN 参数调整(仍偏 BSE 假设),使用前需审。

## 3. 未确定问题(阻塞运行验证)
1. **MZAMS 0.91 "out of range"**(MIST 表 0.7–150 明明覆盖;换了精确 Z、
   禁了 tabuse_* 仍报)——怀疑 SEVN 库版本 ↔ 表集版本不匹配或网格读取方式
   与预期不同。下一步:用 SEVN 仓库自带示例输入驱动 `sevn.x`(注意其 CLI
   参数需 `-参数 值` 成对,多余位置参数会触发 odd-number 报错);若官方程序
   也失败,升级/更换表集或联系 SEVN 开发者。
2. `sevn.x`/`sevnB.x` CLI 解析怪癖(odd number of option-value arguments),
   尚未找到其期望的调用形式,影响"对照排除"手段。
3. SEVN 库版本较新(cmake 4.4.3 构建),与 PeTar 上游 PR 作者开发时的版本差异
   未评估——`evolve_sevn.cpp` 对 SEVN API 的使用若有版本敏感点,需对照其
   changelog。

## 4. 建议实施顺序
1. §3-1 运行阻塞(否则一切运行期验证无从谈起);
2. D1 种子接入(可复现性是科学产出的前提);
3. D2 GW 反冲接入(核心物理缺口);
4. D4+D5 事件流与 TDE 验证(可与 3 并行);
5. D3 自旋语义、§2 工具残留。
