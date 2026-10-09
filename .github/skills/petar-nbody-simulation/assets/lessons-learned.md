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

### 2026-09-21: 本机缺 libasan RPM,asan 调试目标链接失败但 configure 通过

**Mistake**: `make`/`make install` 在链接 `petar.hard.debug` 与 `petar.dump2test` 时报 `cannot find /usr/lib64/libasan.so.5.0.0`,误以为 PeTar 代码问题排查多轮;实际是机器环境问题。另外先尝试 `LIBRARY_PATH=$HOME/.local/asan-shim make` 注入影子链接脚本,同样失败——Intel MPI 的 mpicxx 包装器不吃这个变量。

**Root cause**: 1) 本机 RHEL8 系统未安装 `libasan` RPM(`/usr/lib64/libasan*` 不存在),GCC 8 自带的 `libasan.so` 链接脚本写死 `INPUT(/usr/lib64/libasan.so.5.0.0)`;2) configure 的 `AX_CHECK_COMPILE_FLAG([-fsanitize=address])` 只测编译不测链接,asan 默认(auto 模式)启用,于是 configure 通过、make 阶段才炸;3) `TARGET` 含两个 asan 目标且 `install: $(TARGET)`,所以 `make install` 必挂;4) 共享盘 intelpython3 带 GCC 8 同版本(soname .so.5)的 libasan 可借用。

**Prevention rule**:
1. 该机 asan 目标链接失败(`cannot find /usr/lib64/libasan.so.5.0.0`)时,标准修复路径:确保 `~/.local/asan-shim/libasan.so` 为影子链接脚本 `INPUT ( /publicfs10/fs10-share/soft/share-soft/intel/2020/intelpython3/lib/libasan.so.5.0.0 )`,并在生成的(已 gitignore 的)Makefile `CXXFLAGS=` 行加 `-L$(HOME)/.local/asan-shim`;永久修复是管理员 `dnf install libasan`,装好后移除补丁;
2. 不要用 `LIBRARY_PATH` 注入——mpicxx(2021.5.0 包装器)会忽略;用命令行/Makefile 内 `-L` 才有效;
3. 运行 asan 调试二进制需 `LD_LIBRARY_PATH` 含借用 libasan 目录(intelpython3/lib),缺 soname .so.5 的库只能用它,不能拿 gcc-13/14 的 .so.8 混用;
4. configure 报告编译器支持某 flag ≠ 链接期可用;`AX_CHECK_COMPILE_FLAG` 的局限在解读 configure 结果时要记得。

---

### 2026-09-21: gcc/14.2.0 module 导出绝对路径 CC/CXX/FC,诱发四连坑;工具链切换必须全量清理

**Mistake**: 用户把登录环境从 oneAPI(Intel MPI + gcc8.5)切换到 `gcc/14.2.0 + openmpi/4.1.0-gcc14.2.0 + gsl/2.8` 后,`make` 连环失败,agent 最初按报错逐个修(先 asan,再 Fortran,再 MPI:: 符号),每修一个又冒出下一个;且第一次在长会话终端里测"module 是否接管 PATH"得到误判(残留 PATH 遮蔽,误以为 openmpi module 不导出工具链)。

**Root cause**: 多因叠加:
1. `gcc/14.2.0` module 会把 `CC/CXX/FC` 导出为**裸编译器绝对路径**(不带 MPI 包装),configure 直接继承 `CXX=<裸 g++14>` → 所有链接缺 `-lmpi/-lmpi_cxx`;
2. configure.ac 的 `AS_IF([test x"$FC" == xgfortran], [FCLIBS+=' -lgfortran'])` 只认**字面量** `gfortran`,module 导出的绝对路径不匹配 → `-lgfortran` 静默丢失;
3. `bse-interface/Makefile`(用户此前已提交)硬编码 `CXX=<裸 g++14>` 且 `FCLIBS` 为空,`petar.bse` 链接同时缺 `-lgfortran` 和 `-lmpi*`;
4. `parallel-random/rand.cxx` 使用 MPI-2 `MPI::` C++ 绑定,链接方必须带 `-lmpi_cxx -lmpi`(OpenMPI 4.1 仍提供 libmpi_cxx,已验证 63 个符号);
5. 全部旧 `.a`/`.o` 是旧工具链产物,make 头依赖不全静默跳过重编(同 2026-09-13 条目),陈旧 `randc.o` 报出源码里根本不存在的 `MPI::` 引用,误导排查方向。

**Prevention rule**:
1. 切换工具链/module 栈后,**必须 `make clean` + 子目录 clean + `rm -rf build` 全量重编**,不要信任增量;排查报错前先看 `.o`/`.a` 是否为旧产物(`strings`/`nm` 查它引用的符号源码里是否存在);
2. 测 module 行为要用 `env -i bash -c 'source ~/.bashrc; ...'` 模拟全新登录 shell,长会话终端的 PATH 残留会得出相反结论;
3. **根因修复(用户 2026-09-22 修正)**:gcc module 把 `CC/CXX/FC` 导出为裸编译器绝对路径,覆盖了 autoconf 对这三个 precious 变量的自动探测(这是 configure 生成 Makefile 的输入,才是问题实质)。正确做法是在 `.bashrc` 的 `module load` 之后 `unset CXX CC FC`,让 configure 自动探测恢复工作(生成 `CXX=mpic++`、`FCLIBS=-lgfortran`);configure 命令行显式传 `CXX=<mpicxx> FC=gfortran` 只是应急替代,不要作为常规规则;
4. 换 `--with-interrupt` 等风味参数重跑 configure 时,对比新旧 Makefile 的 MT_FLAGS/TARGET 名(如 bseEmp→bse 被静默降级),确认没丢原有选项;
5. gcc14 自带 `libasan.so.8`(位于其 install/lib64 且驱动会解析),gcc8 时代"系统缺 libasan RPM"的 shim(`~/.local/asan-shim` + Makefile `-L` 补丁)在新栈下不再需要;
6. configure.ac 的 `test x"$FC" == xgfortran` 字面量比较对绝对路径 FC 失效,上游可改成 basename 比较(待办,勿忘);
7. `bse-interface/Makefile` **同样由 configure 生成**(configure.ac `AC_CONFIG_FILES`,模板 `Makefile.in` 中 `FC=@FC@`/`CXX=@CXX@`/`FCLIBS=@FCLIBS@`),其中的裸编译器绝对路径与空 `FCLIBS` 就是污染环境下 configure 的输出,不是手写硬编码;干净环境重跑 configure 后即生成 `FC=gfortran`/`CXX=mpic++`/`FCLIBS=-lgfortran`。**推论:生成文件(bse-interface/Makefile、根 Makefile)永远靠重跑 configure 修复——手动编辑会被下次 configure 静默覆盖(实测发生过),编辑器 undo 也可能把生成文件退回到污染版本**;判断子目录 Makefile 性质先查 `AC_CONFIG_FILES` 和是否存在 `Makefile.in`,不要凭内容里的绝对路径臆断为手写文件。

### 2026-09-22: sanitizer 可用性必须用链接级测试判定——编译通过≠libasan 可链,且 auto 模式会先被重写成 asan

**Mistake**: configure 原用 `AX_CHECK_COMPILE_FLAG([-fsanitize=address])` 只测编译;oneAPI+gcc8.5(系统缺 libasan RPM)下编译通过、`make install` 在 `hard.debug`/`dump2test` 链接期才炸。首次修复实现还误杀了 auto 模式:auto 在 Linux 下先被重写为 `with_sanitize=asan`,显式失败分支随之触发。

**Root cause**: gcc8 的 `libasan.so` 链接脚本写死 `INPUT(/usr/lib64/libasan.so.5.0.0)`,缺 RPM 即链接失败;编译级测试覆盖不到链接期;`AS_CASE` 分支无法区分"用户显式要求"与"auto 解析结果"(原始值须用重写前记录的 `sanitize_mode_requested`)。

**Prevention rule**: 运行时库可用性检测一律用 `AC_LINK_IFELSE` 走 `$CXX`(同文件 `-no-pie` 检测写法);显式请求失败给可操作报错,auto 失败优雅降级(`auto-no-asan-runtime`,debug 工具不带 sanitizer 照编,`make install` 不再阻断)。三分支已实测:gcc14 链接通过→asan;oneAPI 链接失败→降级;显式 `--with-sanitize=asan`→报错。**加固(同日 review 后)**:链接测试改为手动编译并捕获驱动输出——Intel 经典编译器(icc/icpc)对未知旗标只发 `#10006: ignoring unknown option` 警告仍链接成功,会被链接级测试误判为 asan 可用(实际静默无 sanitizer 且编译警告刷屏);现 grep 该警告判 `no`,icc 栈如实降级 `auto-no-asan-runtime`(实测)。

### 2026-09-22: 站点 GSL 可能带异构 MPI 的 DT_NEEDED——ld "may conflict" 警告必须在 configure 期拦截并自动切静态库

**Mistake**: 本集群 `gsl/2.8` 的 `libgsl.so` 带 `NEEDED libmpi.so.40`(上游 GSL 无 MPI,系站点构建污染);oneAPI/Intel 主 MPI 是 `libmpi.so.12`,纯 oneAPI 运行环境直接 `library not found` 起不来,双 module 挂着则同进程双 MPI 运行时。

**Root cause**: 共享库的 DT_NEEDED 整体进入进程;站点 gsl 在 OpenMPI 环境下构建被污染,集群又无干净 gsl(系统 `/usr/lib64` 亦无)。

**Prevention rule**: configure 已加自动回退:`readelf` 提取 `libgsl.so` 的 `libmpi.so.N`,与 `$CXX` 试链程序的 libmpi soname 比对,失配且存在 `libgsl.a`/`libgslcblas.a` 时自动改用静态归档并打 NOTICE(静态归档无 DT_NEEDED;`nm` 已验证 gsl 无 MPI 符号引用)。oneAPI 2026 与 icc20 两次真实构建 `ldd` 干净。

### 2026-09-22: TARGET 名 token 变化先查 git log 再怀疑 configure 丢参——BTLogH 已无条件编译、目标名不再带 btlogh

**Mistake**: 换 MPI 栈重 configure 后新二进制名少了 `.btlogh`,误判为 configure 旗标丢失(`config.status --config | head -3` 截断多行输出加剧误判),几乎重查 configure 参数。

**Root cause**: commit `197c2b9`("build: compile BTLogH support unconditionally in PeTar")把 BTLogH 改为无条件编译、运行时 `--ar-g-func` 切换,目标名随之去掉 token——二进制名是构建产物命名约定,不是 configure 旗标清单。

**Prevention rule**: 重 configure 后对比新旧 `TARGET` 差异,先 `git log --oneline -- Makefile.in configure.ac` 查命名约定演进;读 `./config.status --config` 输出不要用 `head` 截断。

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

### 2026-09-22: 输入快照必须是命令行最后一个参数——PeTar 无条件取 argv[argc-1] 为文件名

**Mistake**: 把 `input.B` 放在选项中间（`petar -u 1 -i 0 input.B -t 0.125`），文件名被解析为 `0.125`，全体 rank 打不开输入文件 MPI_ABORT；用户 `run.sh` 参照命令同样写错。

**Root cause**: `petar.hpp` 的 `read()` 在有剩余位置参数时执行 `fname_inp.value = argv[argc-1]`——文件名永远取最后一个 token，所有选项必须在其前面。

**Prevention rule**: 组 `petar` 命令时输入快照永远放最后（`petar [options] <snapshot>`）；该 run.sh 已修正。（规则已同步入 SKILL.md。）

### 2026-09-22: BSCC-M9 跨节点 MPI 仅经典 Intel MPI 20.4.3 可用——新栈用户态 fabric 与节点驱动失配（rc_mlx5/UD 端点必崩）

**Mistake**: 按常规假设逐栈排查跨节点段错误，多轮作业后才定位到栈级问题：OpenMPI 4.1.0/4.1.8+UCX 1.18 在 `rc_mlx5` inter-node 端点建立即 SIGSEGV（UCX info 日志直接可见）；oneAPI 2022/2026（Intel MPI 2021.5/2021.18）崩溃或 `UD endpoint unhandled timeout` 挂起。`UCX_VFS_ENABLE=n` 全程生效、纯 MPI（1 线程）也崩、单节点全过——排除 VFS bug/OMP 混用/应用与构建问题。

**Root cause**: 集群 inter-node RDMA 用户态栈失配：较新 UCX/libfabric 无法在节点对间建立 RC/UD 连接，而经典 Intel MPI 2020 系自带 fabric 层正常（作业 5394248：跨节点 16×8 ×3 + 16×1 全过；5394302/5394303 完整 0.125 Myr + 快照读回验证通过）。

**Prevention rule**: 该集群跨节点一律 `mpi/intel/20.4.3`（`MPICXX=mpiicpc CC=icc ./configure`）+ `srun` + `I_MPI_PMI_LIBRARY=/opt/slurm/slurm/lib/libpmi2.so`；作业脚本自包含 module；模板 `~/data/N1k/sub.template.sh`，完整证据与工单 `~/data/N1k/ib-crossnode-mpi-failure-report.md`。OpenMPI 若必须用：显式 hostfile + `--oversubscribe`（plm slurm 下槽位推算偏少会拒启）+ `--mca plm slurm`；注意 `hostname` 排布检查不是 MPI 程序，不能证明 wire-up 正常。

### 2026-09-22: sbatch 默认继承提交 shell 的全部环境——同 soname 的 MPI 库会被残留 module 污染

**Mistake**: 会话 shell 残留 oneAPI 2026 的 `LD_LIBRARY_PATH`，`ldd` 与新提交作业把经典 Intel MPI 构建的二进制解析到 oneAPI 的 `libmpi.so.12`（两代 Intel MPI 共享 soname），差点得出错误测试结论。

**Root cause**: `sbatch` 默认导出提交环境；不同代 Intel MPI 的 `libmpi.so.12` 路径不同，`LD_LIBRARY_PATH` 顺序决定胜者。

**Prevention rule**: 提交前清理 shell 无关 module，或在作业脚本内 `module load` 后加 `ldd $(readlink -f $(command -v petar))` 路径断言（模板已内置 `case "$MPIPATH" in *intel/2020*)`）。

---

### 2026-09-23: BSCC-M9 上多个 PeTar 作业并发互相拖慢 10–40×——基准与生产须串行/错峰

**Mistake**: N=10⁶ 基准时 6 个配置作业同时提交（各自独占节点，互不同节点）。全部作业每步 2–7 s（同配置单独跑仅 0.15 s），两轮（none/cores 绑定）皆然，几乎废弃全部数据。各 rank Total 逐位一致（全体同步等待签名），计算相位（PP_single/Tree_Force）与单独跑完全一致，膨胀全在未计时残差（hard 部分同步等待），突发式波动。

**Root cause**: 共享资源争用（最可能 publicfs10 共享文件系统元数据：每步每作业 ~50 个文件追加写 × 6 作业；petar 自己的 Output 计时不覆盖该等待）。作业独占节点不能隔离共享 FS/基础设施；该集群按核调度、夜间批处理重载（421/522 节点占用）时风险更高。

**Prevention rule**: 该集群上基准作业一律串行执行（`sbatch --dependency=afterok` 链，模板 `~/data/N1m/submit_seq.sh`）；生产运行避免多作业同时写同一共享 FS；对比实验若出现"各 rank 同步变慢 + 计算相位不变 + 残差膨胀"签名，先怀疑并发/共享 FS 干扰，再用单作业复测。

### 2026-09-23: MPI+OMP 启动显式 `--cpu-bind=cores`——`none` 在多 NUMA AMD 节点慢 2–3×；OMP 每 rank 多线程已验证有效

**Mistake**: 沿用 SKILL 旧建议 `--bind-to none`（防单核限制），N=10⁶ 基准实测 32×8 配置下 none 比 cores 慢 2–3×（稳态 ~0.5 vs ~0.16 s/步）。历史上"设 OMP=8 每 MPI 最多 2 线程"的担忧在本栈（Intel MPI 20.4.3 + icc + libiomp）未复现：探针实测每 rank 9 线程（1主+8 OMP）、全节点 32×8=256 线程、线程分散不同 CPU；8×32（32 线程/rank）是 256 核最快配置（134 ms/步 vs 32×8 的 180 ms），线程失效时这不可能。

**Root cause**: `--cpu-bind=none` 下线程在 256 核多 NUMA 域间浮动，缓存/NUMA 亲和反复失效；显式 cores 绑定每 rank 固定 8/16/32 个互不重叠 CPU。线程创建与亲和限制本身无问题（旧故障应为当年 mpirun 默认绑定到单核所致，非 OMP 本身）。

**Prevention rule**: 该集群 MPI+OMP 一律 `srun --cpu-bind=cores`（配合 `-c <线程数>`）；`--cpu-bind=none` 仅作 cores 异常时的对照。验证线程是否真并行：作业内 `ps -eo pid,comm | awk '$2 ~ /^petar/'` 取 pid 后读 `/proc/<pid>/status` 的 Threads 与 Cpus_allowed_list + `ps -L` 逐线程 psr（勿用 /proc/<pid>/exe readlink，见同日条目）。N=10⁶ 规模每 rank 线程越多越快（8×32 ≥ 16×16 > 32×8）。

### 2026-09-23: BSCC-M9 跨节点"能跑"≠可用——N=10⁶ 同核数拆双节点即慢 6×（512 核比 64 核慢 5×）

**Mistake**: 基于 N1k 结论"经典 Intel MPI 20.4.3 跨节点可用"直接把 512 核基准排成 2 节点；实测 64×8 每步 1.9 s，比单节点 256 核慢 10×、比 64 核还慢 5×。对照实验（32×8 同 256 核拆双节点）确认：计算相位不变，全部膨胀在跨节点同步等待残差，每 rank 8 线程仅 ~200% CPU 占用（多数时间阻塞）。复测一致，非环境波动。

**Root cause**: N1k 的"可用"只验证了正确性（N=10³ 通信量极小）；该集群 inter-node RDMA 用户态栈故障（见 `~/data/N1k/ib-crossnode-mpi-failure-report.md`）在经典 Intel MPI 下表现为能通信但性能严重受损——N=10⁶ 每步大量跨节点集合通信将其放大为每步秒级等待。

**Prevention rule**: 该集群上 N=10⁶ 生产运行与基准以单节点（≤256 核）为上限；跨节点方案须等集群 IB 工单解决并重测（对照法：同核数 1 节点 vs 2 节点，比值即跨节点代价）。N1k 小规模跨节点验证结论不可外推到大 N 性能场景。确需跨节点时的应用层缓解：**低 rank 高线程布局**——2 节点 8×32（4 rank/节点×32 线程）实测 359–458 ms/步 vs 32×8 的 ~1000–1120 ms（快 ~2.3×；管理员扫描 2026-09-23，`~/data/N1m.opt.20260923/optimization.report.20260923.md`）；但计算相位不变、残差仍 ~0.32 s/步，比单节点同布局（134 ms）慢 ~3×——缓解非根治，rank 数决定跨节点通信量。

### 2026-09-23: `petar.find.dt` 外层勿用 srun 启动（会起 N 份）；profile `Total` 列=每步墙钟非累计

**Mistake**: 作业里写 `srun -n 32 petar.find.dt ...`——find.dt 是 bash 包装器，srun 起 32 份各自内部再 srun，互相争抢 allocation（发现于 finddt.log 出现 32 行 commander）。又因 N1k 目录的 profile 只有一行（单次输出），误判 `Total` 列为累计值，用相邻行差分得到负数/巨幅波动。

**Root cause**: 包装器脚本无 MPI 感知；profile `Total` 实为"自上次输出以来的每步墙钟"（输出块标题 "Wallclock time per step"），`N_steps` 列在此构建恒为常量 2 不可用；hard/Hermite 无独立计时列，其耗时 = Total − Σ(列出的相位) 残差。

**Prevention rule**: find.dt 在作业内直接 `petar.find.dt -r "srun -n N" ...`（-r 指定内部前缀，不用外层 srun）；解析 profile 主块行用 `Time ≈ k·DT 且严格递增` 过滤（尾部 FDPS 附加块的 Time 列会撞上小数值）；每步耗时直接取 Total（-o=DT 时每行恰一步），不要差分。

### 2026-09-23: 作业脚本内检测 rank 进程用 ps comm 匹配——`/proc/<pid>/exe` readlink 受 yama ptrace_scope 限制；sbatch --export 变量勿被位置参数覆盖

**Mistake**: 采样器用 `readlink /proc/<pid>/exe` 匹配 petar 进程，登录节点自测通过（假进程是 shell 子进程），作业内却 0 命中（srun 进程由 slurmd 派生，非采样 shell 后代）。另一次脚本里 `TAG=$1` 把 `--export=ALL,TAG=...` 传入的环境变量覆盖为空（sbatch 不传位置参数），探针静默失败一轮。

**Root cause**: yama ptrace_scope=1 下 `/proc/<pid>/exe` 的 readlink 需 PTRACE_MODE_READ，仅对调用者的后代进程开放；`ps -eo pid,comm` 读 /proc 公开字段无此限制（comm 截断为 "petar.mpi.omp."，用前缀匹配）。

**Prevention rule**: 作业内进程检测用 `ps -eo pid,comm | awk '$2 ~ /^petar/'`；`--export` 传参的脚本直接引用环境变量，禁止 `X=$1` 形式；登录节点自测进程检测逻辑时须用非后代进程（否则结果不可信）。

### 2026-09-29: `petar.hard.debug` 的 dump 文件名是位置参数；ASan+Intel MPI 必须 `setarch -R`

**Mistake**: 用 `--dump-filename <file>` 传 dump 文件名——该选项不在 getopt 表内被静默忽略，回退默认 `hard_dump` 后 "Error: filename hard_dump cannot be open!" 中止。首次直接运行又撞 "Shadow memory range interleaves with an existing memory mapping. ASan cannot proceed"。

**Root cause**: `IOParamsHardDebug::read` 的 getopt 串只有 `-m:n:p:h`，dump 文件名取 `argv[argc-1]`（位置参数，用法 `petar.hard.debug -p data.par <dumpfile>`）；ASan shadow 与 Intel MPI so 映射冲突是 ASLR 相关（ELF_ET_DYN_BASE）。

**Prevention rule**: 调用形式固定为 `setarch $(uname -m) -R petar.<family>.hard.debug -p <prefix> <dump文件>`；参数文件用 `<prefix>.par*`（需含 `hermite-acc-offset-sq`、`hermite-dt-max` 的持久化值）。重现结果与生产日志逐位可比（本次 N210k Pal5 dump `dE_SD/Etot_SD` 复现 -1.64123）。

### 2026-09-29: 多线程 stderr 的 message↔文件名配对会错位——按线程交错，须以文件系统为准

**Mistake**: 从 `output` 日志解析 1648 条 "Hard energy significant !"+"Dump file:" 相邻两行当作同一事件，据"最差 |ratio|=998"去 run/ 找对应 dump 文件——文件名不存在（n/t/O/c 均对不上）。

**Root cause**: OpenMP 多线程写同一 stderr，message 行与其 Dump 文件名行可被其他线程插入隔开；相邻行配对随机错位。

**Prevention rule**: 分析 dump 事件时以 run/ 实际文件清单为准（文件数=事件数可先核对），日志只用来取 ratio/dE_SD 分布；需要 message↔文件精确配对时按线程号聚合连续行并对照文件名校验。

### 2026-09-29: DATADUMP 的 backup 每 cluster-步只支撑一次写盘——先发事件会静默消耗掉后续 dump

**Mistake**: 给 `hard_large_energy` 加重复抑制后，冷坍缩测试首例事件 `allow=1` 却既无文件也无 "Dump file:" 打印，一度怀疑 map/关键键逻辑；更早还把 `dumpThread` 的 `dump_once_flag` 误当死参数（只读了函数前半）。

**Root cause**: `dumpThread` 尾部 `if (dump_once_flag) hard_dump[ith].backup_flag = false;`——backup 在 `driveForMultiCluster` 每 cluster 每 tree 步取一次，之后**第一个** DATADUMP 写盘并清零 backup_flag；同步内后续任何 DATADUMP 静默 no-op（无警告、无打印）。测试中 `dump_binary_merger` 在积分中途先触发，把步末能量检查的 dump 吃掉。

**Prevention rule**: 排查"DATADUMP 无产物"先确认本步该 cluster 是否已有更早 dump 事件（grep 同线程相邻 "Dump file:" 行）；读 `dumpThread` 必须读全——副作用在函数尾部。临时打印 key/allow 一轮即可定位此类问题。

### 2026-10-04: t=0 全簇硬系统连锁的真链路是 r_search 连通簇收集,不是 SDAR 组分析——先查 data.par 解析半径再归因

**Mistake**: `doc/ic-hard-failures-plan.md` 初版把 N=1000 单簇硬系统归因于 "SDAR group analysis 经巨大 r_search 链式合并"(n_group=336→单组),方向误导:真正把全簇连成一个连通硬系统的是上游 r_search 连通簇收集(`cluster_list.hpp` 的 DFS);SDAR 的 r-search-group(物理钳制 a_hs≈0.0068 pc)没有参与——无初始双星 IC 以 n_group=0 失败即为反证。

**Root cause**: `-s 0.5`(dt_soft=0.5 Myr)≫ 该簇自动值(~2.4e-4 Myr),自动反推 r_out_base=3.26 pc(`src/petar.hpp:3807` Kepler 约束,nstep=16,m_avg=2.38 M☉)> 簇维里半径(G·M/3σ²≈2.27 pc);每粒子 r_search≥r_search_min=5.13 pc,再叠加 (m/m_avg)^{1/3} 质量放大(73.6 M☉ ×3.16→r_out=10.24 pc)与双星成员 |v|·dt 项(783 pc),连通判据覆盖全簇→999/1000 粒子一簇。outA/outB 日志中 r_in/r_out/r_search 与公式逐位吻合。双星成员 783 pc 的成因是时序空窗:每树步 `calcRsearchAndGetMassBackupAndsetGroupDataToCM`(`src/hard.hpp:4152`)对已标记组员(bid≠0)本用**质心速度**算 r_search 并复制给成员,但首树步建组发生在刷新**之后**,原始双星成员当时 bid==0 全按单粒子自身速度刷新;建组后同树步只 floor r_search 到 r_out(`syncMemberChangeoverScale`,commit 2fa15b7)不覆盖为 CM 值,783 pc 判据随即进入 SDAR 扰动者判据(hermite_integrator.h:1579 取成员 max)。`-b` 仅影响 IC 载入时一次性的 σ-based r_search(petar.hpp:3937)、σ 统计的成员配对(3684)和 tree 容量(3245),不进入运行期半径公式。

**Prevention rule**: 归因 "整簇进硬积分" 类故障,第一步读 `data.par`/`data.par.hard` 的解析半径(r、r-ratio、r-search-min、r-search-group)与日志的 "Mean inner/outer changeover radius / Average mass / Velocity dispersion",用 r_virial≈G·M/(3σ²) 对照 r_out;检查 -s 是否被显式放大到自动值之上。HARD_DEBUG 版的 "Maximum rsearch particle" 打印直接给出肇事粒子及其半径。

### 2026-10-05: 破组 changeover 残留的验证必须用 data.process 成员判定——status 列与 Kepler 双星判据都会假阳性;修复量化必须有 A/B 基线

**Mistake**: 对症检查初版用快照 `status<0` 当组员判据,758 个"违规"全是活跃双星成员(输出时刻 status 不可靠);二版要求所有 data.process 双星带 pair 缩放,又把宽双星(SDAR 不建组、按设计带自身值)误判——两版都差点得出"修复无效"的错误结论。

**Root cause**: 成员身份的权威来源是 SDAR 分组;`status` 列在输出时刻被重置/不完整,`petar.data.process` 的双星是 Kepler 束缚判定,与 SDAR r_group 判组是两套标准。正确不变量:单星 r_in == 自身质量公式;双星成员接受"pair 值或自身值",匹配两者皆非才是陈旧。终版结果:基线(无修复)561 例陈旧、最高 +260%;修复版(sync+簇守卫+孤立单星守卫)0 例,唯一成员案例 0.12%(SE 漂移量级)。

**Prevention rule**: 写 changeover/成员不变量检查用 petar.data.process 的 single/binary 分类,不用原始快照 status 列;双星成员判据必须同时接受 pair 与自身两种合法值,容差 ≥2e-4;量化修复效果必须带无修复基线(stash→重编译→同 IC 重跑);修复版与基线中止于同子系统不同断言 = 既有缺陷而非回归。**第三类假阳性:层级组(n≥3)**——data.process 只识别双星,n3/n4 组成员按设计带组总质量 changeover,与 pair/自身期望都不匹配;严格判定需用 data.group.nN 成员表。MPI 验证(2 秩,bse,N=500 双星簇,t=2 Myr):干净完成、hard_connected=9+1(跨节点连通簇路径与 adr<0 远端分支实际执行,守卫经既有 send-back/correctForceForChangeOverUpdateOMP 机制传播),所有不变量标记均属三体假阳性(经 data.group.n3 确证,ids 47/48/371、449/450)。复现配方(mcluster -N 500 -b 0.95 -m 0.5 -m 50 -C 5 -u 1 -s 7 → petar.init bse 模板 → petar -u 1 -b 237 --bse-metallicity 0.02)记录于 2026-10-05 断言修复提交信息。

### 2026-10-05: 双星中断断言与 SEVN merge 无关——归因先 blame 失败行再 cpp -E 剥离对比;固定 --rand-seed 也不可复现

**Mistake**: BSE 双星两个断言(ar_interaction.hpp:1230 `dt>0` 与 :1342 零质量星)恰好在 SEVN 集成提交(55d09bd/1a12bab)合入次日复现,时间上高度可疑;且 ±hard.hpp 修复 A/B 表现为"换个位置崩",容易被读成"merge 引入新缺陷"。

**Root cause**: (1) 归因依据是时间邻近而非代码谱系:git blame 显示两条断言及其触发判据(isCallBSENeeded 的 Roche/trflow 分支)分别引入于 2021-01、2024-04、2024-07、2024-10,远早于 2026-10-04 的 merge;用 `g++ -E` 携带目标构建完整 -D 集(SEVN 未定义)对 bse_interface.h 做 merge 前后真实预处理对比,BSE 可见差异仅剩:readAscii 错误消息计数修正(fscanf 本来就按 12 列读)、#endif 后分号重排、static constexpr 重构、选项解析重构——无数值路径变化(手写 ifdef 剥离脚本连续两版都有 bug,不可靠)。(2) 运行对照方法错误:默认 `--rand-seed` 固定后仍不可复现——同一二进制、同 seed、同线程数(1 rank/2 threads)两次运行分别崩在 hard.hpp:3116(NaN vbk)与 ar_interaction.hpp:1342,BSE 中断处理的线程交错使每次运行对"崩在哪"重新采样,单次崩溃签名不能区分代码变体;而 pre-SEVN 二进制(worktree 1a89cd2)同输入 seed 1 即精确复现 Case-B 断言。

**Prevention rule**: ① 归因"最近合并"前先两步:blame 失败行与触发判据的引入日期;对涉事头文件用 cpp -E(带目标构建完整 -D、目标宏未定义一侧)做前后对比,不要肉眼扫 diff 或手写宏剥离脚本;② BSE 家族运行对照必须重复多次比较"崩溃集合"是否一致,`--rand-seed` 固定不构成可复现性保证(结论自包含:两条断言行均 2021/2024 引入,cpp -E 全 -D 对比 merge 前后 BSE 可见代码无数值路径变化,pre-SEVN worktree 二进制同输入复现 Case-B 断言);③ 对称性检查:一侧修了某判据缺陷(如 1a12bab 只给 SEVN 分支加了 dt>0 guard)后,必须检查其他 interrupt 分支(BSE/mobse/bseEmp)是否需要等价 guard。

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

### 2026-09-27: 列数不匹配穷举 kwargs 后必须停下上报——读取入口文档须自带终局规则，仅指针不够

**Mistake**: Pal5 会话读 2021 年 `data.status`（73 列 vs bse+galpy 期望 75 列）时，穷举 interrupt×external 组合仍不匹配后没有停下上报，而是两次自写 workaround（先 `np.loadtxt` 手工取列；被用户指出后改为裁剪到 global 块宽 + `readArray`）。SKILL.md「Snapshot Read-Mismatch Policy」明确禁止此做法，但未被遵守。

**Root cause**: 政策硬规则只存在于 SKILL.md 尾部；Python 读取的强制入口 `data-readback-patterns.md` 只有 kwargs 匹配说明加一行指向政策的指针（2026-06 已存在，事故仍发生），其 Warning Classification 只写 "Stop, correct flags, retry" 而无终局上报步骤，且 "Column-count mismatch fallback" 示范了裁列读取（Profile 专用但无作用域限定），直接诱导同类做法。

**Prevention rule**:
1. patterns 文档「Common Parameters」现自带「Column mismatch: stop rules」：穷举 kwargs（含版本相关行）→ 仍不匹配 → 停止并上报（文件模式、reader 类、kwargs、报错原文）；禁止自写 parser 或裁列；Profile fallback 已显式限定仅 `data.profile`。
2. 写读取代码前通读 patterns 文档含 stop rules，不能只看 kwargs 表。
3. SKILL.md 政策段已声明同样适用于旧版本输出（schema 代差按 reader/producing-solver 版本缺口上报，而非 mode flags 错误）。

### 2026-09-27: 旧格式列差先做 git 考古——pre-2024-12 输出的 +2 列是 star.spin 1D→3D，官方读法是 spin_3d=False

**Mistake**: 初步诊断把 2 列差怀疑为 `mass_bk`/`status`、`r_in`/`r_out`、`pot_*`（全部错误——这些字段 2021 已存在），未先做 git 历史 diff，险些立项开发"旧格式读取器"。

**Root cause**: 差异实际在嵌套块 `SSEStarParameter` 内部：`spin` 1D（SSE 块 10 列）→ 3D（12 列），C++ 侧 2024-12-23 引入（`5c96ac1`）；Python 侧 2026-02-03 已加 `spin_3d=False` 读取选项（`669b6cb`），但 `Status`/`Particle`/escaper/`GroupInfo` 的 docstring 与 patterns 文档 kwargs 表均未提及，按文档穷举 kwargs 永远试不到它。

**Prevention rule**:
1. 列数差鉴定：`git show <同期 commit>:tools/analysis/<file>.py` 对比 keys 列表；优先检查嵌套子块的版本开关（`SSEStarParameter`/`BSEBinaryEvent` 的 `spin_3d`、`BSEDynamicMerge` 的 `less_output`）。
2. 已知 legacy 读取项：pre-2024-12 BSE 系输出 → `spin_3d=False`；pre-2020-09 并合记录 → `less_output=True`；pre-2020-12 galpy 快照 → `petar.format.transfer -c`；pre-2020 group data → `-g`。
3. 新增影响列布局的 kwargs 时，同步补齐所有 kwargs 转发容器的 docstring 与 patterns 文档 kwargs 表，并检查 `petar.data` / `petar.get.object.snap` 等 CLI 是否透传（本次已补 `spin_3d`/`collect_sp_acc` 文档，CLI 增加 `--spin-1d`）。

**Update（同日）**: 崩溃点已在 SDAR 侧修复——`DictNpArrayMix.readArray` 列数不匹配由 warning 改为抛带类名、列数与 reader kwargs 的描述性 `ValueError`（含 legacy 提示），列过剩静默丢列一并消除；详见 SDAR lessons-learned 2026-09-27「Python Tools」条目。

### 2026-09-27: petar.data.process 三层误导错误的根因是错位读取宽容 + auto-resume 无可观测性

**Mistake**: 漏 `-i bse` 的 `petar.data.process` 崩在 `np.histogram` "bins must increase monotonically"，与根因（interrupt 模式不匹配）相距三层：错位读取仅警告 → 垃圾"双星" c.m. 除零 NaN → histogram 表面报错。崩溃残留的 realtime partial 又被后续重跑的 auto-resume 静默合并追加，多轮日志呈现"不同失败"，成为排查耗时主因。

**Root cause**: 二进制字节错位在 sdar `fromfile` 中只警告不终止（已于 SDAR 侧改为默认抛错）；恢复/去重/跳过全程不打印记录数与时间范围，跨 key（lagr/core/tidal/bse_status）时间集不一致也无任何提示。

**Prevention rule**:
1. 错位读取现第一步即抛带 kwargs 提示的 `ValueError`；`Lagrangian.calcOneSnapshot` 计算前校验质量 NaN/inf——错位但字节碰巧对齐时也能早失败。
2. `petar.data` 恢复时打印每 key 记录数与时间范围；输出 key 文件时间覆盖不一致时显式警告（多由不同 flags 的中断运行造成）。
3. 排查 data.process 崩溃先看第一条 warning 而非最后的 traceback；重跑前清理 `.parallel.*` partial 并核对恢复日志。

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

### 2026-09-17: 规则正文夹杂论证性叙述——agent 文件按调用加载同样吃上下文

**Mistake**: Model Allocation 的 Tier rules 以"规则+论证+示例"的散文体写入 agent 文件（15 行），每次 Developer 调用整体加载；其中论证性内容（钱包枚举理由、家族匹配充分性讨论、失效成本说理）在委派时刻不可执行，且"never escalate well-specified task"与同文件 Delegation Threshold 已有规则重复。

**Root cause**: 把"agent 定义只承载委派规则"理解为只约束**主题**（内容属于委派域即可），未约束**文体**（可执行规则 vs 论证叙述）；且 agent 文件不在"always-loaded"字面范围内，上下文效率要求未被显式适用于它。

**Prevention rule**: 规则正文 = 可执行指令 + 至多一句反直觉理由；论证、历史、示例归 lessons。写入后 grep 同文件相邻 section 去重。该文体约束已固化到 `.github/skills/README.md` Update Rules（agent 定义条 + 预算条）。
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

### 2026-09-17: frontmatter `model:` pin 与 doctrine 矛盾——模型分层只走委派参数，不走 pin

**Mistake**: 工作区（未提交）给 `petar-developer.agent.md` 和 `petar-reviewer.agent.md` 的 frontmatter 重新加了 `model: GLM-5.3 ...` pin，而同一文件的 "Model Allocation" 仍写着 "No agent pins a model in frontmatter"（2026-09-12 移除 pin 的决策记录）。pin 使用的 picker 显示名是否可解析未经验证——正是 2026-09-12 记录的静默回退失效模式，会产生"已强制强模型"的假信心。

**Root cause**: 模型分层的意图（conductor/reviewer 强、executor 高效）没有权威落点，实现时直接改了 frontmatter 而未对照同文件内已有的分配规则；frontmatter 编辑也不会触发阅读正文 doctrine。

**Prevention rule**: 模型分层的唯一权威位置是 `petar-developer.agent.md` 的 "Model Allocation"（委派时 `model` 参数：Reviewer 默认强档、Simulation Engineer 默认 Flash 档 + 根因未知时升档）。不改 frontmatter pin。若确要为"picker 直接调用"加 pin 强制，必须先验证 pin 真实生效，并在同一次改动中更新 Model Allocation 文本，消除矛盾。（同日三探测后修正：pin 机制验证可行，策略见"政策先于验证落盘"一条。）

---

### 2026-09-17: 委派 `model` 参数失效即报错并列出模型目录——Flash 名可探测解析，"never guess" 对委派参数不成立

**Mistake**: 在未经探测的情况下，把 frontmatter pin 的"静默回退"失效模式外推到委派时 `model` 参数上，写成 "never guess — an unresolvable name silently falls back"，并把 Flash 名解析设计为"问用户一次"。四探测验证（继承对照 / 假名 / 猜测名 / 目录精确名）证明：委派参数解析失败是**响亮报错并附完整可用模型目录**（工具层拒绝，零推理成本），不运行任何回退；精确到目录字符串的名字（`GLM-5.3-Flash (ZhiPu AI (Coding Plan)) (unify-chat-provider)`）成功让 subagent 跑在 Flash 上。

**Root cause**: 两种机制（frontmatter pin vs 委派参数）失效模式不同，但 doctrine 只基于 pin 的历史教训（2026-09-12）写作，未区分二者；且假设"委派时无模型目录"，实际目录可通过一次无效名探测免费获得。

**Prevention rule**: 委派参数允许试错——假名以零成本换回完整目录。但**选名不是纯技术决策**：目录混合 Coding Plan / Free / PayGo / copilot 路由 / BSCC 等互不相同的计费来源，同族匹配只是能力约束，不是计费约束；tier→名的映射必须由用户指定（或按其约束表选择），不得由 agent 擅自硬编码。frontmatter pin 的谨慎仍然成立（其静默回退未被推翻）。任何"X 会静默失败"的规则写入配置前，先用零成本探测验证 X 的真实失效模式。

---

### 2026-09-17: 政策先于验证落盘——"不走 pin"一日内两次反转

**Mistake**: 基于未验证的历史教训（2026-09-12 "pin 静默回退"）把 "No agent pins a model in frontmatter" 写成硬规则，并据此移除了用户手工添加的 pin。用户质疑后三探测（T1–T3）确立完整机制：委派参数 > frontmatter pin > 继承；pin 用目录精确字符串完全生效（委派路径亦然）；仅 pin 名失效时静默回退到继承——该降级对成本分层是安全的（回退=强档，成本回退而非正确性故障）。同一条规则一日内两次反转。

**Root cause**: 把历史教训当作当前环境的机制事实，未区分"当时的失效可能源于名字格式"与"机制本身不可用"；移除用户手工配置前未先用其精确值做零成本探测。

**Prevention rule**: 涉及机制行为的配置规则（pin/参数/优先级/失效模式）落盘前必须探测验证；移除用户手工添加的配置前，先验证该配置是否实际生效。验证后的最终策略：Flash 默认 = Simulation Engineer 的 frontmatter pin（用户指定，目录精确字符串，写入后探测一次确认生效）；升档 = 委派参数（覆盖 pin，失效响亮报错）；pin 静默回退可接受，供应商/订阅变更后重测 pin。

### 2026-09-17: AGENTS.md 里的条件指针未被遵循——执行者自己的 agent 定义需要一跳直达指针

**Mistake**: 用户在独立会话中让 Agent 修改 agent 规则文件，"改任何定制文件前先读 `.github/skills/README.md`" 的指针（写在两个 `AGENTS.md`）未被遵循，直接动手修改了文件。

**Root cause**: 指针距执行者两跳（agent 定义 → `AGENTS.md` → README）。`AGENTS.md` 虽常驻加载，但条件式指令（"改 X 前先读 Y"）未被主动对照时约束力是概率性的；而唯一有权修改定制文件的角色（PeTar Developer）自己的定义里没有该指针，Reviewer 的检查清单里也没有对应的验收项。

**Prevention rule**: 把指针直接钉进 `petar-developer.agent.md`（Conductor Rules 2 + Key Reference Files）与 `petar-reviewer.agent.md`（PeTar-Specific Checks 第 5 项），一跳可达。指针不是规则复述，不违反 "维护规则只住在 README" 的归属约束——Design Goal 4 允许 pointer 出现在一切需要的发现点。

### 2026-09-20: DSM 并合后零质量粒子与并合残骸精确重合——eps=0 下 0/0=NaN 毒化能量记账

**Mistake**: `calcMergerProperties` 把去活化的零质量粒子 `p0` 停放在 `-pm->pos`（关于 group c.m. 的镜像点）。对双成员 group（AR 在 group c.m. 系积分，质量加权和为零），`pm->pos` 可舍入到精确 0（等质量对必然），镜像点与残骸**逐位重合**；eps_sq=0 时 AR/Hermite 力与能量循环算出 `0/0 = NaN`，依次毒化 `epot_`/`etot_ref_`/`gt_kick_inv_`/`de_change_interrupt_`，group break 时经 `accumDESlowDownChangeBreakGroup` 进入 `energy_.de_cum`，最终触发 `hard.hpp` 的 `ASSERT(!std::isnan(energy.de))`（生产运行表现为并合事件后能量误差永久异常）。2026-09-16 的修复（case 2，残骸重锚定）只覆盖了另一条路径，未覆盖此重合。

**Root cause**: 815af98 用精确镜像替换旧的 `pm->pos*(1+1e-8)+1e-12` 停放时，未考虑双成员 group 中 `pm->pos == 0` 精确成立的情形；NaN 只在**精确**重合时出现（r≠0 时零质量对贡献恰为 0），因此只有个别并合触发——表现为随机、不可必然复现。Makefile 的 `HARD_SRC` 依赖表漏掉 `disk_star_merger.hpp`，头文件改动不会触发重编，也会让"修复未生效"假象。

**Prevention rule**: 任何"停放/重生"粒子的代码必须保证与所有粒子保持有限距离（对零质量粒子，任意 r>0 的相互作用严格为 0，只有 r==0 产生 0/0）；DSM 侧已改为重合时按并合前分离向量偏移 `1e-3*|dr0|`。调试此类"偶发 NaN"用条件 watchpoint（`condition N isnan(*(double*)addr)`，注意用裸地址而非符号名——跨帧求值会失败）。`disk_star_merger.hpp` 已补入 `HARD_SRC`（Makefile 与 Makefile.in 同步）；对仍可能未列入依赖表的头文件，改动后需 `touch` 其 includer 或先对照 `HARD_SRC` 检查。同函数内 `pm->dm` 累加顺序（先加后赋值 mass）已对齐其他质量变更点的模式（2026-09-20 修复，原顺序恒加 0）。

### 2026-09-28: Makefile DEBFLAGS 的 += 块在 `ifneq(debug_mode,no)` 内默认全灭——改编译宏必须 `make -pn` 验证生效值

**Mistake**: 欲在生产主程序开启 hard_assert 的 ASSERT（`-D HARD_DEBUG`），先把宏加进 DEBFLAGS 的 `+=` 块（Makefile 238-250 行），重编后二进制零断言字符串；`make -pn | grep ^DEBFLAGS` 才发现该块整体位于 `ifneq ($(debug_mode),no)` 内，默认 configure（debug_mode=no）下全部无效，实际生效的只有块外无条件的 `-D AR_DEBUG_DUMP -D HARD_DUMP`。此前"生产构建带 AR_DEBUG（kickVelIter 精确质量求和断言激活）"的推断同样基于这块死代码，均为错误。

**Root cause**: 条件赋值与无条件赋值混排的 Makefile 里，读片段推断生效 flags 不可靠；块内本就有 HARD_DEBUG（hard.debug 语义），易误以为生产目标继承。

**Prevention rule**: 增删编译宏后必须 `make -pn | grep ^DEBFLAGS` 验证生效值，并 `strings <binary> | grep '!ISNAN(dt)'` 确认断言表达式真正嵌入（`__FILE__`/表达式字符串是活断言的标志）。HARD_DEBUG 现已移入无条件区（`DEBFLAGS += -D HARD_DEBUG`，AR_DEBUG_DUMP/HARD_DUMP 旁）。A/B 实测（functional `std` 双并合用例，N=100, -u 1, -b 4）：断言版死于 `ASSERT(!ISNAN(integration_error_rel_abs))`，而无断言旧版同点位 `|Int_err/E|` 已是 nan——watchdog 类断言无误报，抓到的是真实损坏；同时实测确认负 dt 无条件缩减（2026-09-28 SDAR 提交）让新版越过旧版致死的 negative-dt streak abort。该用例当前仍红：interval 1 内出现真 NaN（旧版亦有，早于本周所有改动；疑与 2026-09-20 DSM 0/0 同族），待用 SDAR 恒等式残差 gdb 法定位。

### 2026-09-28: T3/T4 验证管线对 `-b 1` 场景三处断裂——petar.init 缺 `-s bse`、管线强制重选 plain 家族、Status reader 缺 interrupt_mode

**Mistake**: 为验证 SDAR ds 改动跑验证层矩阵：T2 直接全绿，但 T3/T4 管线全断——(1) 场景 setup 的 `petar.init -f input input.base` 未加 `-s bse`，而运行命令全部带 `-b 1`，输入列数在读取期即不匹配（"requiring 4, only obtain 2"）；(2) `run_validation.py` 的 `_select_or_build_petar` 每次调用 `petar.select --optional avx2,omp` 把用户事先选好的 bse 家族切回 plain（症状同上但更隐蔽——手动选好后管线内部又改掉）；(3) 分析端 `Status(N_particle=...)` 不带 `interrupt_mode`，bse 列使 dtype 错位。另：`-f` 前缀残留旧输出时 petar 自动重启，HARD_DEBUG 断言 `id>0` 拦截损坏 status（正确行为，但排障时需先清 `t3_*/t4_*` 残留）。

**Root cause**: T3/T4 场景定义于 bse 化运行命令，但管线从选二进制、生成输入到读输出全链路都默认 plain 特性集——三层没有任何一处校验特性一致性；套件长期只被 plain 家族跑过（T2），bse 路径从未端到端通过。

**Prevention rule**: (1) t3/t4 场景 setup 已补 `petar.init -s bse`；extract 配置新增 `"interrupt_mode": "bse"`，`extract_orbital_drift_from_status` 增加 interrupt_mode 参数透传给 Status reader；(2) 跑 T3/T4 需用场景模式显式指定二进制：`run_validation.py --scenario <json> --var petar_bin_switch=<bse二进制>`（管线包装器的自动选择与此类场景不兼容，待加 per-scenario 家族声明）；(3) A/B 基线构建时警惕 configure 状态漂移：验证管线自身的 `_select_or_build_petar` 会静默 `./configure` 重置根 Makefile 特性集（本次两次"基线"二进制因此作废）；(4) 排障 `-f` 前缀运行先清残留输出。

### 2026-09-29: 精确文本替换在 Markdown 行尾双空格（硬换行）上两次失败——先 `cat -A` 诊断再改编辑策略

**Mistake**: 批量下沉 SKILL.md 章节到资产文件时，两处大段 `replace_string_in_file` 连续失败（find.dt 工作流段、hard dump 段），第一次失败后仍按原思路整段重试，浪费两轮；根因是段落内某行以两个行尾空格（Markdown 硬换行 `for binary.␣␣$`）结尾，整段 oldString 无法逐字节匹配。

**Root cause**: Markdown 硬换行的行尾双空格在普通读取视图中不可见；整段替换对隐藏空白零容错。

**Prevention rule**: 整段精确替换失败一次后，立即 `sed -n 'A,Bp' <file> | cat -A` 检查目标区间的隐藏空白/制表符，再决定：把编辑拆成以问题行为边界的两段，或改用智能编辑工具；不要原样重试。另记录本次分层下沉预算（复查修正后终值，供下次健康检查对比）：SKILL.md 565 行/46.6 KB → 387 行/33.3 KB（-32%/-29%），新增资产 ic-generation / changeover-tuning / build-toolchain / controlled-experiments，DSM 硬规则并入 dsm-workflow，hard dump 调用规则并入 script-tools，场景问询表去重后单一家园为 minimal-question-sets.md。

### 2026-09-29: dsm-workflow `--radius` 示例值与换算节自相矛盾（差 20×）——单位示例必须与同文件换算表核对

**Mistake**: SKILL 分层下沉后的实测（5 个全新上下文子代理纸面演练，追踪读取路径）发现 `assets/dsm-workflow.md` 的 `--radius` 旗标表示例值 `4.5092203040509496e-08` 标注 "(0.2 au in pc)"，与同文件 "Radius calculation" 节的 0.2 au = 9.69627362e-07 pc 矛盾（示例值实为 ~0.0093 au）；照抄示例会把盘外缘缩小 20×。同批发现并修复：default-postprocessing 的 `petar.movie -i` 枚举漏 `dsm`；changeover-tuning 缺孤立双星 r_in 定尺寸规则与 Gate 5 优先级说明（已补）。

**Root cause**: 同一物理量在两节各自写数值、从未并列对照；分层重组时只搬运未做数值一致性核对。实测还暴露的已知残留（未修，待决）：`-G 0.00449830997959438`（default-postprocessing 模板）与 `petar.G_MSUN_PC_MYR = 0.004498502…`（data-readback）不一致；binary-scenario-map 无 DSM 条目；`petar.HardData`（主 debug.log）在 11 个 Pattern 中无覆盖。

**Prevention rule**: 搬运或新增任何带物理数值的示例时，grep 同文件/同主题文件中该量的其他出现并对照；数值不一致即缺陷，当次修掉。纸面演练（新上下文子代理 + 读取路径追踪 + "将执行的命令"对照规则原文）是验证 skill 路由与规则保留的低成本手段，重组/大改后应跑。

### 2026-09-29: OMP_NUM_THREADS 规则藏在资产层且为条件式——反直觉环境规则必须无条件进常载层

**Mistake**: 实际会话中小 N 模拟未设 `OMP_NUM_THREADS`，OpenMP 默认吃满全部核心，占用巨大计算资源且效率更低。"Parallel Launch Heuristics"（N≲10³→1 线程）与实测数据（N=500：4 线程 +19%、8 线程 +66%）早已在 `script-tools.md`，但 (1) 常载层 "Environment requirements" 只要求 `OMP_STACKSIZE`；(2) 资产层措辞是条件式——"set OMP_NUM_THREADS explicitly **when user asks for a concrete launch layout**"——用户不主动问就不设。

**Root cause**: 反直觉规则（OpenMP 默认全核对小 N 是负优化）停留在按需读取的资产层且带触发条件；主流模型默认不设置线程数，任何条件化都会让规则在"简单运行"场景静默失效。

**Prevention rule**: 执行安全类环境规则（线程数、栈大小）必须无条件写入 SKILL.md 常载层并进入 Gate 3 确认摘要的可视项；资产层只保留 N→线程数对照表与实测数据。已落实：Environment requirements 新增无条件硬规则、Gate 3 摘要第 2 项显式含 `OMP_NUM_THREADS`、script-tools 条件式措辞改为 "every launch command"。

### 2026-09-29: "not a production solver" 只说不是、没说是——hard.debug 被误用为主程序 debug 版

**Mistake**: 实际使用中 `petar.hard.debug` 被误当作"主程序的 debug 版"用于跑/调试模拟，而非其真实用途：重放运行产生的 dump 文件（输入是 dump 文件，不能从快照启动模拟）。binary-scenario-map 的 Helper 条目只写 "diagnostics and hard-integrator debugging; not a production simulation executable"——负面否定 + 模糊正面描述正是误用诱因；且 script-tools 工具清单缺少 `hard.debug`/`dump2test`/`hard.test`/`format.transfer`/`petar.bse` 系的逐条用途条目。

**Root cause**: 工具分类采用"不是什么"而非"是什么 + 输入形态"；清单以常用工具为主，helper 与 standalone 工具无正向用途描述，模型只能按名字猜测语义（"debug" 后缀天然诱导"调试版主程序"解读）。

**Prevention rule**: 每个已安装工具必须有正向用途 + 输入形态描述（输入是什么文件、能否从快照启动模拟），helper 条目同时保留"不可替代 solver"的否定句。已落实：Gate 4 显式声明 hard.debug 是 dump 重放工具、binary-scenario-map Helper/standalone 条目重写为正向用途、script-tools 新增 "Solver-adjacent binaries and standalone tools" 一节（用途经源码 usage 文本核实）。

### 2026-09-29: 新硬规则须与同主题既有规则并列对照——never-delete 规则资产层只留指针

**Mistake**: 独立复查（Reviewer 对照 skills README 契约，8 项必查）抓出两类问题：(1) OMP 硬规则补写时写 "N ≲ 10³ 用 1 线程" 未加 binaries-rich 限定，与同仓 "Parallel sizing quick rule"（双星多的小 N 建议 2–4 线程）及 sample 脚本相抵，且实测样本（N=500 含 1 个双星）不支持无条件推广；(2) ic-generation.md 下沉时把三条 never-delete 规则全文复制进资产并与 SKILL.md 双向声明所有权（"restated"/"owns the detail"），违反 "Never restate a fact in a second place"。另抓出 6 项建议（机制细节复述、quick rule 回声、controlled-experiments 复述 gates、lessons 预算数字未随微调回写、README 章节名指针错误）。全部已修。

**Root cause**: 增写规则时只对照触发事故本身，未 grep 主题关键词把同主题全部既有表述拉出来并列对照（矛盾/张力在孤立看各自都正确）；下沉时用"复制+互指"代替"指针+增量"，把单一家园契约软化成了双家园。

**Prevention rule**: (1) 新增或修改任何规则前，grep 主题关键词（如线程数、-C 5、hard.debug）把所有出现并列对照，矛盾或未限定的一般化即缺陷；(2) never-delete 规则的唯一全文在 SKILL.md 常载层，资产层只允许指针加不超过一行的增量说明；(3) 记录在 lessons 的预算/测量数字在后续微调后必须回写；(4) 结构性改动（分层、下沉、规则强化）完成后跑一次独立复查（Reviewer 按 skills README 健康检查程序）。

### 2026-10-01: 冒烟"回归"三重假象——VERSION 匹配跳过重建跑陈旧家族二进制；长寿组暴露成员 changeover 漂移；numpy2/ffmpeg 环境坑

**Mistake**: SDAR 统一判据落地后跑 `test/functional` 冒烟，bse-galpy 稳定 SIGABRT，一度判定为新判据回归并用 stash 双向重建做 A/B。结论部分正确但过程被三个无关因素污染：(1) `install_petar_major_versions.sh` 在 installed VERSION == source VERSION 时**跳过重建**——改源码不 bump VERSION 时，冒烟跑的是陈旧家族二进制（bse 案例用 `petar.mpi.omp.avx512.bse`，与手动 make install 的 avx512.bse.galpy 不同族），失败被错误归因；(2) 真实根因是**长寿组成员 changeover 漂移**：CM 的 r_in/r_out 每步按组总质量重算，成员只在成组时缩放一次，BSE 质量漂移累积突破 DEBUG 断言（1e-3）与 `r_search>=r_out` 不变量——旧判据频繁进出组掩盖，新判据（设计上）保持紧组更久从而暴露；(3) `np.in1d`（numpy2 移除）与 imageio 无 ffmpeg 后端让全部案例倒在工具步。

**Root cause**: 成员同步只存在于 `collectGroupMemberAdrAndSetMemberParametersIter`（成组时一次性 r_scale_next）；质量变化路径（`correctSoftPotMassChange` 只修能量簿记）无 changeover/r_search 刷新；诊断依赖行号匹配源码但二进制来自 VERSION 缓存。

**Prevention rule**: 源码级 A/B 或回归验证前必须核对运行二进制的 mtime/家族与源码一致（`ls -la` + 断言行号对照），VERSION 不 bump 时安装脚本会静默跳过重建；修复方式：`syncMemberChangeoverScale()`（每个组级 `pcm.changeover.setR` 后把 CM changeover 复制给成员并 floor r_search）+ 硬域入口 `updateWithRScale()` 后 floor `r_search>=r_out`（4 处组初始化位点 + 入口守卫）；成员与 CM 在同一边界一起跳变，与 CM 自身 setR 的即时性对称。已验证：双并合冒烟 IC 连续多轮通过、诊断断言零违例、Pal5 811 dump 回放与修复前逐位一致。

### 2026-10-02: "SE 质量损失"诊断被事件级回放证伪——large_energy 真凶是双曲近遇的组进出注入；修复=事件步重基线

**Mistake**: 对 Pal5-IMF 生产运行新出现的大量 `hard_large_energy` dump，交接计划文档诊断为"BSE 质量损失能量修正未接线"（dm 记账 active、消费者被注释），并给出改 dm 簿记的方向。事件级回放（gdb 断点 `evolveStar` + ADJUST_GROUP_DEBUG 组事件打印）三案例全部证伪：步内 12 次 evolveStar 全部 dm=0 或 ~1e-14，`.sse.*` 事件文件为空；每个案例的步内都发生组 form/break（40.5+40.5 M☉ BH 对在无束缚/临界束缚近遇的近心点瞬间 d<r_crit 成组、越界即释放，两例还有 form→break→re-form→break 抖动），每次状态改写注入 O(1e-3)×相遇动能（16.4 / 1.12 / 0.138）。`dE_mod`、`dE_change` 恒 ~0 因为过渡注入从未进簿记列——参考相对 dE 携带全部偏移并反复触发告警。

**Root cause**: 判据性诊断只看了"dE 在首子步即达终值且恒定"的形态与 dm 簿记的存在，未做事件级归因（组事件计数/断点验证）；过渡注入正是此前回文工作已测得的层（判据对称≠过渡对称），在稠密星团 BH 近遇上以 ~40 dump/小时的量级显形。

**Prevention rule**: large_energy 归因必须先做事件级验证：`grep "Find new group\|Break group"` + gdb 断点 `evolveStar` 看 dm，再谈修正方向；"dE 恒定于首子步"同时兼容参考偏移（组事件）与真实误差，不能单凭形态定罪。修复（`hard.hpp` integrateToTime 能量读出后）：事件步 `calcEnergySlowDown(true)` 重基线，注入计入 `dE_change` 簿记列，告警只看残余积分误差——三案例回放 dE 降至 ~1e-10、零触发，非事件步与事件后漂移的检测灵敏度保留。

### 2026-10-03: 重基线修复撤回——de_change_cum 语义是"物理变化专用"，算法注入混入即污染能量诊断

**Mistake**: 前条把组事件注入计入 `de_change_cum` 的"治标"修复，经用户审查撤回：`de_change_cum`/`de_sd_change_cum` 设计语义是**物理**能量变化（恒星演化质量损失、双星中断）；把算法误差折入等于让簿记列失去诊断意义、并掩盖真实误差信号。告警压力本身是"过渡层注入真实误差"的正确信号。

**Root cause**: 修复定位（事件步重基线）在诊断上有效（三案例 dE→1e-10）但语义错误；正确的根本方向是消除注入本身——SDAR 过渡层逐事件互逆 + 派生量纯态化（form/break 互为逆正则变换、辅助量由当前态经共享代码推导、切换对齐共同同步边界），可执行计划见 `SDAR/docs/transition_unification_plan.md`。

**Prevention rule**: 数值方案中的簿记列各有语义边界（物理变化 vs 算法残差 vs 参考重置），任何"把 X 挪到 Y 列"的修复必须先核对 Y 的设计语义；治标修复落地时必须在文档中标注其掩盖性质与撤回条件。

### 2026-10-03: SEVN 上游合并四连坑——非 SEVN 构建静默失效;依赖仓版本配对;自动合并丢 #ifdef 块;同名工具误导冒烟

**Mistake**: 审视 SEVN PR(#79) 合入 master 时,静态审阅只抓到 Makefile/tools 层问题;实际构建验证又连续暴露:(1) `bse_interface.h` 的 Fortran `merge_()` 调用丢分号、两处空 `#elif`、`../src/io.hpp` include 被挪进 `#ifdef SEVN`——非 SEVN 构建(bse/mobse/bseEmp)全部编译失败,而 PR 作者只测过 sevn 模式(preprocessor 直接剪掉 #else 分支,语法错误不可见);(2) `Makefile.in` 的 `se_modfe` 与 `bse-interface/Makefile.in` 的 `BSE_LIB_S` 两个变量名 typo 让 `--with-interrupt=bse` 静默退化成无 BSE 构建、测试工具丢失 `-lbse`——无任何报错;(3) 用共享 checkout 验证 PeTar master 时误配 SDAR experiment(`adr` 私有)与 FDPS master v8.0(profile API 改名),`adr`/`collect_sample_particle`/`IOParams<PS::S64>` 三类错误全是版本配对假象,PeTar master 实际配 FDPS ≤v7.1c + SDAR master;(4) 冒烟误用 `~/bin/petar.bse`——那是 bse-interface 的独立测试工具(bse_test.cxx),与主程序家族二进制同名。

**Root cause**: 上游 PR 只在新增特性模式下编译验证,`#ifdef` 新分支天然掩盖旧分支的语法破坏;Makefile 条件分支 typo 属"静默功能退化"类,configure/编译都不报错;多仓工作区(共享 ../SDAR、../FDPS checkout)里分支配对关系(哪边动 API 哪边适配)没有显式记录;`petar.bse` 一名二物从未被文档强调。

**Prevention rule**: (1) 审阅/合并 stellar-evolution 类 PR 后,至少在 bse 模式完整 configure+make+微运行(10 分钟内),不接受"只看 diff";(2) 验证某分支的 PeTar 前先确认依赖配对:PeTar master↔SDAR master+FDPS ≤8.0 前版本、PeTar experiment↔SDAR experiment+FDPS 8.0,用临时 worktree 钉住依赖版本(`git worktree add --detach ... <tag>`),不碰共享 checkout;(3) 大型合并解决冲突后,用脚本 diff 两边所有特性关键词(SEVN 等)出现的行集,审计自动合并是否吞掉 `#ifdef` 块,再数 `#if*/#endif` 配平;(4) 冒烟主程序用 `petar`(经 petar.select)或完整家族名(如 `petar.mpi.omp.avx512.bse.galpy`),`petar.bse` 永远指测试工具。

### 2026-10-04: SEVN 首次真编译挖出七处合并暗伤;运行瓶颈在 SEVN 库-表匹配,非集成层

**Mistake**: SEVN 模式长期只有静态审查+伪造安装 dry-run,从未真编译。装好真 SEVN 后首次编译连续暴露 7 处:① bse-interface 子目录用顶层相对路径找 SEVN 头文件(IO.h not found);② `IOParamsBSE::idum` 成员声明在合并中丢失但 SEVN 构造函数仍初始化它;③ SEVN 构造缺 `fname_par` 初始化(IOParams<string> 无默认构造直接编译错);④ 保存参数文件块缺 `#elif SEVN` 分支;⑤ `isCallBSENeeded` 的 SEVN 分支少一个 `}`(整个 BSEManager 类花括号失衡,导致后续文件所有报错全被误导成"constexpr 成员"之类);⑥ `petar.sevn` 链接缺 `-lrand`;⑦ initdata.sh sevn 模板 19 列应为 18 列。此前"宏配平(#if/#endif)检查"无法发现 ⑤——花括号平衡是另一维度。

**Root cause**: 上游 PR 只在自己环境测过;合并时我保留的"实验侧 if 块开括号 + 上游侧结尾"混搭在条件分支组合下漏了一个闭合;dry-run 与静态检查都不能替代真实编译。运行期:PeTar 侧输入读取(42 列/星 18 列)、表加载均通过,瓶颈在 SEVN 库内部——默认表是 PARSEC 大星表(2.2–600 M☉),低质量星必须换 MIST 表并禁用 `--tabuse_rhe/rco/envconv false`(MIST 表缺这仨文件),之后仍报 `MZAMS 0.91 out of range`(表网格 0.7–150 明明覆盖),SEVN 自带 sevn.x 也因 CLI 解析怪癖跑不起来,无法对照排除——库版本与表集的匹配问题,非 PeTar 集成层。

**Prevention rule**: ① SEVN 类条件合并后必须真编译(sevn + bse 两模式),检查包括:#if/#endif 配平 + 花括号净深度(`cpp -E -DSEVN ... | 数 {}`)两层;② 子目录 Makefile 引用仓库根相对路径时必须做 `$(patsubst ./bse-interface/%,./%,...))` 之类的上下文换算或用 @prefix@ 单独导出;③ SEVN 运行环境三要素先备齐再跑:MIST 表(`--tables`)+ `--tabuse_rhe/rco/envconv false` + 质量下限 0.7 M☉(MIST)/2.2 M☉(PARSEC 默认表);④ 首跑报"out of range"先查所用表集的质量/金属丰度覆盖,再怀疑集成代码。

### 2026-10-04: value3_/rand3_ 是上游 PR 携带的死种子路径——"看似在播种,实则无人消费";sevn.x 要求 argc 为奇数

**Mistake**: SEVN 分支的 `BSEManager::initial()` 一直在写 `value3_.idum`(含 MPI rank 偏移),代码形态与旧 BSE 的 COMMON /VALUE3/ 播种完全一致,审阅时自然认定"SEVN 种子已接入"。实际上:SEVN 库既不导出也不读取 `value3_`/`rand3_`(`nm libsevn_lib_static.a` 无此符号,头文件零引用),写入的是 PeTar 侧自建的孤儿结构;而 BSE 家族的 kick 早已改走 Fortran `rand_f64()`(PeTar 并行随机,`--rand-seed` 播种),`-idum` 选项同样无人消费。死代码让"随机数不可复现"这一真实缺陷被掩盖。另外 sevn.x 的 CLI 解析要求 argc(含 argv[0])为**奇数**——所有选项必须严格 `-名 值` 成对、不接受无值开关,否则直接报误导性的 "odd number of option-value arguments"(此时 argc 其实是偶数),无法用作对照排查。

**Root cause**: 上游 PR 按旧版 SEVN 的 BSE 式接口写种子代码,SEVN 库后来换成自带 `utilities` 随机数,接口蒸发但调用侧残留;静态审阅只看"写没写",没验证"谁在读"。CLI 坑在 SEVN 源码 `src/general/params.cpp` 的 `SEVNpar::load`(n%2==0 即报错)。

**Prevention rule**: ① 接线随机种子(或任何跨库状态)时必须验证消费端存在:查库符号导出(`nm`)与头文件引用,写而无读即删除,不留"安慰性"代码;② 判断某 CLI 选项是否仍有效,反向 grep 其 `.value` 的全部消费点,而不是看帮助文本;③ 排查 sevn.x 时所有参数写成 `-名 值` 成对形式(布尔也用 `-名 true/false`),argc 保持奇数。

### 2026-10-04: SEVN "MZAMS out of range" 是表覆盖稀疏不是库表错配——gold 表实测 0.7–80 M☉;运行验证须防已安装工具过期

**Mistake**: SEVN 运行长期被 `MZAMS 0.91 out of range` 阻塞,当时推断为"库版本↔表集版本不匹配"并写进了 SKILL/plan。实测(petar.sevn 逐质量探测)发现:非 gold MIST 表(AGBnotpedantic)在部分 Z 下 ~1 M☉ 以下覆盖稀疏(SEVN 的 `tables/tables_info.md` 明说未达 CO 燃烧的轨道被裁剪),报错信息引用的 0.7–150 是标称范围不反映实际网格;换 `SEVNtracks_MIST_AGBrobust_gold` 后 0.7 M☉ 起全通(0.65 失败,90+ 失败@Z=0.00142857)。另两个坑:① SEVN 安装器**不拷贝 tables/**,`--tables` 必须指向 SEVN 源码目录;② 端到端首测报"requiring 18, only obtain 17",根因是 `~/bin/petar.init` 还是修复前的 19 列旧版——仓库 tools/ 改过后必须重新 install,否则浪费排查时间在代码侧。

**Root cause**: 报错范围来自表元数据标称值,真实覆盖由轨道裁剪决定;错误推断("版本错配")未被实测证伪就写进了文档;已安装工具与仓库 HEAD 漂移没有检查习惯。

**Prevention rule**: ① SEVN 质量/Z 越界先做逐点探测(`petar.sevn --tables <集> -z <精确Z> 质量1 质量2 ...`),再读 `tables/tables_info.md` 的覆盖说明,gold 表优先;② 凡"工具行为怪异",先 `diff` 已安装副本与仓库版再查代码;③ 运行期新结论必须实测(本条:gold 表 + N1k 无双星全程跑通 + type_change 事件落盘)后才更新 SKILL/plan,推断性结论要标注"未验证"。

### 2026-10-04: SEVN dt=0 断言——unbound 对无条件调用违反主调 dt>0 不变量;修复验证必须核二进制 mtime>头文件 mtime

**Mistake**: N1k 无双星 SEVN 运行在 t=7.75 Myr(BH 形成后)触发 `ar_interaction.hpp:1230 (dt>0)` 断言。直觉归因方向全错:先怀疑 SEVN overshoot(dtmiss<0 使 time_record 超前)——独立工具逐星实测排除;真因是 dt **恰好为 0**(bit 级 `time_record == time_interrupt == time_now`,gdb 十六进制比对确认):被 SN 踢出的 BH 与 MS 星新形成**双曲对**(semi<0, ecc 1.88),两成员 SSE 状态已同步到 time_now,而 `isCallBSENeeded` SEVN 分支对 unbound 对(`_semi<=0`/`_ecc>=1`)无条件返回 true、不查 dt——主调方的 `ASSERT(dt>0)` 不变量被破坏。非 SEVN 分支无此规则,故 bse 免疫(对照实验存活只是巧合性证据)。修复:unbound 对仅在 `max(dt1,dt2)>0` 时调用,否则推迟到下一中断点。另一坑:首次修复后回放仍复现断言——worktree 头文件 mtime(17:33)晚于二进制(12:04,时钟跳变导致 make 未感知头已更新),编译的还是旧头;gdb 断点行号与源码不符是唯一线索。

**Root cause**: 上游 SEVN 分支新增"负 SMA/高心率即调用"判据时,没有与主调方 `dt = time_now − max(time_record) > 0` 的前置约定对齐;dt=0 调用本就无演化量可做。构建侧:文件 mtime 时钟跳变使 make 依赖失效。

**Prevention rule**: ① 时序断言类崩溃,先 gdb 十六进制比对涉事时间量(==0 与 <0 的区分直接决定归因方向),再用同 IC 的对照模式解释差异来源,不接受"对照组存活"作为充分证据;② 给 `isCallBSENeeded` 类判据函数加新调用条件时,必须核对主调方对该调用成立的全部隐含前提(此处是 dt>0);③ cp+make 之后、回放之前,核对 `binary_mtime > header_mtime`,不符则 touch 强制重建——gdb 断点行号与源码不符即为陈旧二进制的特征信号。

### 2026-10-04: SEVN 静态库默认带 OpenMP 符号——拉取其归档成员的工具只需链接期旗标;-Dopenmp=OFF 因 TLS 错配与线程安全被否决

**Mistake**: 主树 `make install` 在 `format.transfer` 链接失败(undefined `omp_get_thread_num`)。误判一:以为 worktree 曾成功构建该目标即"免疫"——实为其产物从未真正重链(文件缺失),旧规则强制链接同样失败;误判二:尝试以 `-Dopenmp=OFF` 重建 SEVN 求干净——SEVN 头文件在 `_OPENMP` 下把 `liststars` 等静态流声明为 `threadprivate`(IO.h),PeTar 的 OMP 构建编译 `evolve_sevn.cpp` 产生 TLS 引用,与无 OMP 库的非 TLS 定义错配无法链接;且 petar 主程序在 OMP 并行区(integrateGroupsOneStep→evolveStar)、`petar.sevn` 逐星并行均调用 SEVN,线程安全恰恰依赖其 openmp 构建的 threadprivate。终案:SEVN 保持默认 openmp=ON,拉取 SEVN 归档成员的目标(format.transfer,因显式链接 evolve_sevn.o)加**仅链接**旗标 `OMPLINKFLAGS`(不带 `-D PARTICLE_SIMULATOR_THREAD_PARALLEL`、不产生任何并行);不引用 SEVN 符号的目标(simd.test/tt.test)静态库成员不被拉入、无需任何旗标;链接集统一为 `SELIBS`。另:hard.debug 回放退出时 LeakSanitizer 报 ~6.5 kB 字符串持有(IOParamsContainer::readAscii/add_sevn_param)属静态生命周期噪声,非泄漏,`ASAN_OPTIONS=detect_leaks=0` 可静默。

**Root cause**: SEVN cmake `option(openmp ... ON)` 默认开启,日志用 `omp_get_thread_num` 随库分发;静态库成员拉取规则(无未解析符号不拉入)决定哪些目标受影响;SEVN 头文件的 `_OPENMP` 条件编译使库构建选项与调用方编译旗标强耦合。

**Prevention rule**: ① 判断某工具是否需要 OpenMP 链接旗标:看它是否显式链接 SEVN 对象(`SE_OBJS`)或引用 SEVN 符号——显式列出的静态库仅按需拉成员;② 改动外部库构建选项前,先 grep 其头文件的条件编译宏(`_OPENMP` 类)与 PeTar 侧编译旗标的耦合,并确认调用是否发生在 OpenMP 并行区;③ 链接期 `-fopenmp` ≠ 启用并行(零线程占用),与编译宏 `THREAD_PARALLEL` 必须区分使用。

### 2026-10-04: hard.debug 回放改 par 判据半径无效——判据状态内嵌在 dump 粒子里;判据参数效应归因先看事件卡在哪个分支

**Mistake**: 想用 `*.hard.debug` 回放做 r-group 参数扫描,直接改 `data.par.hard` 的 `r-group`/`r-search-group`(0.1×–10× 共 30× 范围),结果三档 dE 与有效 r_crit 逐位不变,一度误以为"生产事件与半径无关"。

**Root cause**: 有效 r_crit = 每粒子 changeover 半径(dump 内嵌)× r_group_over_in,回放工具从 dump 恢复粒子态时判据半径随之固化,par 文件的 r-group 行不重算它(×30 还会触发 r_search 界断言)。真正该看的:Pal5 三案形成点全部精确卡在半径边界(dr≈r_crit_eff)、κ_org 超门 7–13 个量级——事件是**半径钳制**而非 κ 门限,κ 分支(自适应放置器,样例实测校准在最优)被 r_group/r_in=0.00375 上限拦住。

**Prevention rule**: 回放只能复现不能重parametrize 判据;判据参数的生产效应要么改 dump 二进制里的粒子半径,要么跑真实段。归因任何判据行为先打印事件触发点的 dr/r_crit/κ_org,判断卡在半径分支还是 κ 分支,再谈参数。


### 2026-10-05: 防御性 clamp 照搬前先查危害是否已被根修——"单星路径有防护"不是加 clamp 的理由;isCallBSENeeded 的 SEVN 轮询是死代码

**Mistake**: 修复双星两断言后,我在 postProcess 又加了 `std::max(time_interrupt_max, time_now)` clamp,理由是"单星路径 modifyOneParticle 已有同款防护,注释列了两个危害"。用户追问后查证:两个危害都是已修复的 bug——hard dump 时间偏移(40b1ee4 已在 hard_debug 源头修复,clamp 只是残留 harden)、Hermite dt-max 超 drift 步宽(1bc65f9 已修)。Case-A dump 里 stale-max 也从未发生(触发是正常排程 ti=rec+dtstar 恰落在 time_now)。同轮还发现 isCallBSENeeded 的 SEVN 分支里 `getTimeStepBinary(...) < dt_sevn` 轮询自 SEVN 集成起就是死代码:束缚对走第一判据,else 分支要求束缚对与 !call_flag 矛盾;SEVN 下一事件时间本就由 postProcess 的 getTimeStepBinary 写进 time_interrupt,不需要检查时二次查询。均已移除。

**Root cause**: 把"别处有防护"当成"此处需要防护",没有验证危害在当前代码是否可达;把 2024 年注释列举的 hazard 当作活机制,而没查这些 hazard 对应的根修提交。SEVN 死代码则源于集成时照搬了不成立的假设(检查时轮询替代排程),从未做可达性推演。

**Prevention rule**: ① 加防御性 clamp/assert 豁免前,先回答:当前代码里这个危害可达吗?用 git log -S 找防护代码的引入提交,读它修的根因是否已在源头修复;若已根修,防护只是残留,照搬等于掩盖未来的真 bug。② 不变式类(time_interrupt_max >= time_now)设计上应严格成立,违反即 bug:暴露优于 clamp,clamp 只应存在于已知的工具性时间偏移场景(如 hard_debug 回放)且需注释指向根修提交。③ 评审继承的判据代码时做可达性推演:第一判据已覆盖的状态,后续分支是否可达;不可达即删。

### 2026-10-05: BSE 双星两断言修复——dt>0 不变式要求两成员都落后;kw15 是 BSE 合法返回通道必须在 assert 中豁免

**Mistake**: (1) Case A(dt>0 断言)最初归因为 Roche/trflow 非时间判据在 dt=0 时触发;实际 trflow.f 对已溢出恒星返回正值(一个轨道周期),dt=0 永远不会触发——真机制是新 group 形成时成员 time_record 混合(一个已在 time_now、另一个因 SSE 事件 dtmiss 落后),isCallBSENeeded 第一判据 `(_dt1>0 or _dt2>0)` 只要有一个落后就调用,catch-up 后二体区间恰为 0。(2) Case B(零质量星断言)假设是 merger 簿记泄漏;实际 hrdiag.f 三处(544/798/979)在无残骸超新星(PISN/碳点燃失败)时合法返回 kw=15/mt=0,evolv2 内部产生时只记 type-change 事件(event_flag≤2),断言在 postProcess 已有的零质量清理之前就 abort。两次都是 hard.debug dump 重放 + gdb watchpoint(bin_interrupt.status)定案,静态假设全部偏了。

**Root cause**: (1) 断言编码了错误的不变式("dt>0"实际要求两个成员都严格落后于 time_now),而混合 record 状态在新鲜 group/重组时常规出现;(2) 事件分类根因:evolv2.f 的 kw15 转变经 goto 135 本就写 marker 12(no remnant),真正的漏在 PeTar 侧 evolveBinary 的 kw15 早退分支——SSE(evolv1)产生的 kw15 无二体事件通道,该分支把状态记成 type 2 且仅伴星有事件才记录;修正为始终记 type 12,isNoRemnant→event_flag=5,断言门自然关闭,断言恢复严格(豁免移除)。消失处理最终为"每处自含":叶子 catch-up(catchup_modify==3→status=merge 且 modify_return=max(...,3))、叶子 postProcess mass-zero(既有)、树分支单成员注册点(branch==3→merge,带 destroy 保护);最初的集中提升块已删。中途教训:catch-up 注册曾写死 modify_return=2,把 3(消失信号)吞掉——父树节点永远看不到 3,hard.hpp 非 merge 路径的 ASSERT(!isUnused) 会炸;返回码语义必须透传(std::max),不得整体覆盖。诊断方法论:t=12.5 快照全净 + binary_merge 无新行即可排除 merger 泄漏,锁定窗口内 SSE/evolv2 直接产生。

**Prevention rule**: ① 排查 BSE 断言先取 dump 实际状态再验证触发过程:状态假设(成员 time_record/time_interrupt/kw/tphys)用 hard.debug -m 2 静态读回,触发序列用 -m 0 重放;不要从判据代码反推触发路径——反推两次都错(Roche 假设、merger 泄漏假设);② 对 BSE/SSE 返回值写 invariant assert 时必须核对 Fortran 侧全部产出通道(kw15 三处、mix.f 一处),"mt>0"类断言要豁免 kw15 或改为检查 kw 一致性;③ 中断状态机:凡设 status!=none 必须同时 setBinaryTreeAddress(SDAR checkParams 强制),编辑该区域后若新增状态赋值路径,grep 确认每个路径都带地址;④ 修复后必须重放原 abort dump 回归——本次重构(catch-up 前置)曾意外丢失入口 setBinaryTreeAddress,gdb watchpoint 立即暴露。

### 2026-10-05: 写自定义 dump/格式读取工具前先查 skill 工具清单——静态读回缺口当日由 hard.debug -m 2 关闭

**Mistake**: 排查 BSE 断言时,为回答"dump 备份态里各粒子的 time_record/time_interrupt/kw"这个静态问题,直接手写了 struct 布局式 dump reader(/tmp/read_hard_dump.cxx),而没先查 skill 的 script-tools.md——其中 "Hard dump debugging" 节已完整覆盖 hard.debug 的调用规则、setarch、par 拷贝、版本匹配与 dump 分类。用户质疑后改用 hard.debug,重放输出(group 形成序列、interrupt 上下文)才是定案证据;手写 reader 则因 PtclHard 布局与写方二进制不完全一致输出乱码。

**Root cause**: ① 拿到"读数据"类需求时惯性写一次性脚本,没走"先查 skill 工具清单"的流程;② 工具定位与字段覆盖分辨不足(实证 2026-10-05):hard.debug=全状态重放;dump2test=转换为 petar.hard.test 输入,仅 20 列动力学字段(mass/pos/vel/id/bin_stat/r_search/r_in/r_out),**丢弃 t_record、t_interrupt、dm/radius 与全部恒星参数(kw/tphys/mt...)**——即 BSE 中断调度诊断所需的簿记字段;"含簿记字段的 dump 静态读回"当时确无文档化工具(当日已由 hard.debug -m 2 关闭,勿再视为缺口),正确反馈路径是 skill/工具改进,而非布局脆弱的临时脚本。附:dump2test 亦为 ASan 构建,调用需 setarch $(uname -m) -R 前缀(当时未注明,已补入 script-tools.md)。

**Prevention rule**: ① 涉及 dump/快照/事件文件的读取分析,第一步 grep SKILL.md 与 script-tools.md 的工具清单;不覆盖再考虑自写,且自写前评估格式耦合度(内存布局式格式严禁临时脚本);② 临时脚本连续失败两次即停,换 skill 文档化的正式工具;③ 发现 skill 工具矩阵缺口(如 dump 静态读回)时,记入 lessons 并在维护窗口提议:给 hard.debug 加 --print 模式或在 script-tools.md 增加读回配方,由维护流程决定。(已实施 2026-10-05:hard.debug -m 2 静态读回,含簿记字段,列图见 script-tools.md;经 Case-A dump 验证与已知真值一致。)
### 2026-10-05: checkConsistence 质量备份断言在 STELLAR_EVOLUTION 下结构性无效——断言点位于 setMassBackup 刷新之前,固定容差挡不住合法质量损失

**Mistake**: pulsar 分支 debug 模式跑 binary_test 在 T=18.49 abort 于 `artificial_particles.hpp` checkConsistence 的 `abs(mass_cm_check-getMassBackup())<1e-3`,形态像质量簿记 bug。addr2line 定位真实触发点是 hard.hpp `driftClusterAndArtificialCMAndWriteBack` 中 AR 积分刚结束处:该 checkConsistence 执行在 BSE 步内质量演化之后、`setMassBackup()` 刷新(hard.hpp ~1914)之前,比较对象天然是"新成员质量和 vs 旧备份快照"。任何单步超 1e-3 的恒星演化质量损失(星风/SN/质量转移,pulsar 测试必含 SN)都必然触发——断言编码的不变式与代码自身数据流矛盾,不是 bug 信号。

**Root cause**: 备份同步是离散事件(每硬步末刷新;软侧 CM 质量由 calcRsearchAndGetMassBackupAndsetGroupDataToCM 在刷新后另行同步),BSE 恰在两个同步点之间改成员质量;1e-3 容差只是对单步质量损失的任意猜测,SN 超它数个量级。同构隐患:`correctOrbitalParticleForce` 内 `m_ob_tot` vs backup 的 1e-10 断言(orbit 采样粒子开启时同样会比较刷新前备份,当前默认关)。真正有效的簿记校验已存在(如 `ASSERT(getMassBackup()!=0.0)`、writeBack 的 dm 记账)。

**Prevention rule**: ① 写/保留备份类调试断言前先排数据流:断言点必须落在"生产者已写、消费者未改"窗口内,备份不变式只能在刷新点之后紧邻校验;② STELLAR_EVOLUTION 构建下任何固定容差的质量差断言都无效(单步质量损失无上界),应改用簿记类校验(非零、dm 记账);③ abort 栈先 addr2line 定位真实调用点再判断断言合理性,不要只看断言表达式。处置:删除该断言的 STELLAR_EVOLUTION 分支,保留非 SE 构建的 1e-10 版本(无质量演化时不变式成立)。

### 2026-10-06: 大质量天体"逃逸"合理性必须先做守恒律审计——快照总动量 P(t) 与双星硬化能量预算比 data.core 异常更早定罪

**Mistake**: run_N1k_200Myr(N=1000, W0=5, Rh=1 pc, 583 Msun, Kroupa IMF, 纯引力, 自动 changeover: r_out=0.079/r_in=0.008 pc)在 t≈47-48 出现"77.6 Msun 巨双星(id 184+377, 13% 簇质量)携 ~40 颗近邻星逃逸、data.core 随之漂移到 (-110,18,-31) pc"。若只看 data.core/data.esc 或逃逸速度量级(1.5-2.1 km/s,与该处 v_esc 同量级),会误判为物理弹射。

**Root cause**: 质量分层后巨双星沉入核内;changeover 半径按 (m/m_avg)^{1/3} 质量放大(77.6 Msun → ~5×)叠加 r_search 链式连通(见 2026-10-04 条目),把大块核心粘进**单个连通 Hermite 硬簇**(内嵌小 AR 组,非巨型 AR 组——`hard.debug` 重放 dump 证实:t≈11 为 111-134 粒子 Hermite 簇 + 1 个 2→3 成员 AR 组;日志 "N_members" 直方图统计的就是连通簇,来源 `petar.hpp:2399-2413` `clusterCount(connected_cluster_n_list)`);t=47-48 该系统硬积分失败:output 日志 E_total -303→-159(+144,占 |E0| 48%,注入走 de_change_cum/Modify 列 +143,Modify_single/Modify_group=0)、|L| 58→270、N_all 1000→1012(≈2-3 个 AR 组的人工粒子);相邻快照全粒子总动量 P: (9.2,-1.7,-2.3)→(-119,8.6,-54.4) Msun·pc/Myr(|ΔP|≈130 Msun·km/s,为事件前 46 Myr 累计漂移 9.6 的 14 倍),且簇剩余部分反冲方向与逃逸组同向(物理反冲必须反向);能量来源不足且动量无载体:全局注入 +144 中双星硬化仅贡献 ΔEb≈+31(Eb=G m1 m2/(2a):370.8→401.8;**勿用 Binary.mass(总质量)代 m1m2——Eb=G·m_tot/(2a) 无物理意义,首次分析即错成 21.4**);动量侧决定性:双星 CM 获得 |Δp|≈114 Msun·km/s,事件中的三体入侵者星 87 仅 0.242 M☉、事件后速度 0.9-2.8 km/s(动量 ~0.25-0.67,需 ~471 km/s 才能平衡),簇剩余部分反冲方向还与逃逸组同向——逃逸速度是在 AR 组 (87,184,377) form/break 状态回写瞬间凭空写入的(data.group.n2/n3 事件级记录:两次速度阶跃 0.68→1.40→1.77 km/s 精确落在 t=47.8203/48.0664 组解散时刻,期间双星内部 a/e 三位有效数字不变)。dump 重放(petar.mpi.omp.avx2.hard.debug -p data.par,版本比 run 新 9 天但 dE_SD=0.250278 与 run stderr 逐位一致)直接测得该硬系统单树步相对能量误差 ~1e-3(阈值 1e-4,PeTar 自身告警),步内发生 AR 组 form(2→3 成员),与 2026-10-02 组转换注入机制同源。

**Prevention rule**: 判定逃逸/喷出事件是否物理,按序审计:① 相邻快照算全粒子总动量 P(t)(孤立系统严格守恒,比能量灵敏;快照坐标为单一惯性系,初始 P=0 可作基准;output 日志的 "C.M." 行是势中心 _mode=3 不是动量 CM,勿混用);② output 日志事件步的 E_total/Modify/|L|err_cum 跳变;③ 能量预算:逃逸组动能 vs 全部双星硬化增量(事后无残留硬化双星=数值注入);④ 反冲方向核对。大 AR 组(>10 成员)是红旗:致密核+大质量 IMF 下自动 changeover 的质量放大项会让巨双星捕获整核,须显式收紧 -s/-r 并复跑。

### 2026-10-06: 巨双星逃逸 bug 已复现并定位到组 form/break 表示切换的动量失配——诊断用逐层动量探针,硬域内部守恒不代表全局守恒

**Mistake**: 首轮分析止于"守恒律崩坏+组事件时间吻合",未定位到代码;若只看硬积分器内部能量(dE_SD ~1e-3)会低估——实际全局注入达 O(100) Msun·km/s/事件。

**Root cause**: 从 data.47 重启复现(当前代码仍在:-261→+325 能量、|P| 55→306/3Myr;重启快照文件名沿用输入 nfile 编号,re47.48=t=47.25)。逐层动量探针结论:① SDAR Hermite 簇间跳变全部严格配对抵消(两簇经 soft 力交换,net=0);② 泄漏在 petar.hpp drift() 内部:soft 侧实粒子动量在组 form/break 事件单步跳 |ΔP| 50-156 Msun·pc/Myr,连环 form/break 不抵消(净 ~82-88/事件期),与快照全局跳变吻合。机理:组存在期间其动量由人工 CM 粒子(在 n_real 之外,快照/探针都不计)承载;break 时成员恢复质量回写,恢复的动量 ≠ (form 时隐藏的动量 + soft 系统对人工 CM 传递的冲量)。嫌疑位点:hard.hpp integrateToTime 末尾人工粒子更新用陈旧 cm_vel_org(pcm 重算前的漂移值)平移 bink.vel(~2088 行)、updateArtificialParticles 直接取 slowdown 表示下的树根速度、AR perturber 冲量与人工 CM soft kick 的双重/遗漏记账。

**Prevention rule**: ① 动量泄漏诊断必须分层:SDAR 簇内(calcEnergySlowDown 加 ptot_,已实现 getPtotTrue)→ soft 实粒子(petar.hpp drift() 入口,已实现 SoftP drift jump 打印)→ 快照全局;硬域守恒≠全局守恒,表示切换(质量置零成员+人工 CM)是盲区;② 重启复现时快照文件名编号连续自输入快照,审计时用日志时间对齐;③ 诊断探针改动保持未提交并在完成后恢复 configure(已恢复 bse.g);④ 修复方向:break/form 回写时强制 Σm·v 成员 == 人工 CM 动量(边界守恒校正),并修正 cm_vel_org/bink.vel 两处陈旧值;验证即用本探针重跑(净 |dP| 应归零)。

### 2026-10-06: 动量修复两次尝试的排除性结论——匹配路径本已严格守恒(dv~1e-16),退役侧校正反而有害;泄漏点仍在,下一步是步内相位二分

**Mistake**: 修复初版打在两个猜测点上:① retire 侧(组变化/解散时陈旧人工 CM 的动量移交成员)——实测把本来自洽的 form/break 对打错方向(全局 |P| 从 143 恶化到 216);② matched 侧(updateArtificialParticles 后锁定人工粒子组动量==成员 Σmv)——实测 dv 恒为 ~1e-16(机器精度),说明 bink 基础的更新本来就是成员 CM 的精确表示(SD=1 无 slowdown 时),锁定是无效但无害的防御性代码。

**Root cause**: 泄漏的 +85.9@t=47.623 单步跳变既不在 SDAR 簇内部(硬域跳变全部配对抵消)、也不在 matched 人工粒子更新(严格一致),且不是轻星表示翻转(0.24 M☉ 的星 87 不可能贡献 85.9)——只能来自 77.6 M☉ 双星自身的 sys 槽质量表示翻转(createGroup 置零 ↔ breakGroups 恢复)与陈旧人工粒子在移除前的共存窗口。另注意:Ptcl::DataCopy **不拷贝 group_data**(回写会剥离成员标记),任何在 system_soft 侧用 isMember()/getParticleCMAddress() 找成员的方案在 removeParticles 时刻都会空手而归——成员信息权威在 ptcl_hard_。

**Prevention rule**: ① 边界修复前先验证不变量是否已被现有代码满足(matched 锁定的 dv=1e-16 一测即知),别在未测量的猜点上动刀;② 下一步探针:在泄漏步(47.623)内做相位二分——kick 前后、createGroup 前后、drift 各子调用(kickCluster/kickCM/drift/writeBack)前后分别打印 P_soft_real,单步定位产生 +85.9 的确切调用;③ system_soft 侧任何成员识别必须走 ptcl_hard_,DataCopy 剥离 group_data 是结构性陷阱;④ 命令行 CXXFLAGS=... 会整体覆盖 Makefile 的 += 累积(SIMD 标志丢失导致虚假编译错误),加宏应直接追加进 Makefile。

### 2026-10-06: gdb 定位动量泄漏真因——不在 form/break,而是 AR 组隐藏期间每步连续注入 dP=m_组×Δv_cm;表示切换本身守恒;AR 扰动已差分化非双重计数

**Mistake**: 前两轮把修复打在 form/break 表示切换上,方向全错。gdb(-O1 -g,与 -O3 逐位一致)裸地址硬件 watchpoint + 逐步 dump 证实:①"SoftP drift jump"的 ±53/+85.9 大跳变是**表示切换假象**(成员质量 46.06→0、动量转入人工 CM;总质量 583.6→506 即证据;P_full=P_real+P_arti 在切换时刻守恒到 0.02);②真实泄漏是**连续的**:双星处于 AR 组表示期间(n_all=1012, t=47.623-47.7246)每个树步 P_full 增 1.3-8.9(逐步增长),累计 +85.5 与快照 |P| 55→143.7 精确吻合,组一解散泄漏立即停止;③每步注入 dP_full 与人工 CM 速度变化 Δv_cm **方向大小精确一致**(dP=m_组×Δv_cm,方向收敛 (-0.99,+0.14,-0.02))——隐藏组的 CM 被加速(隐含 8-58 pc/Myr²,核内近邻力量级)而**全系统无反冲**(实粒子每步 ΔP≈0.02);④单方受力候选逐一排除:AR perturber 双重计数不成立(ar_interaction.hpp:505-510 明确扣除 CM 净分量做差分化)、matched 人工粒子更新严格一致(dv=1e-16)、kickCM 与成员 kick 用同一 CM 加速度。

**Root cause**: 单向力在 Hermite 簇内组 CM 的受力路径上——组 CM 加速无反应。剩余嫌疑(下一步单一实验可钉死):`need_resolve_flag` 分支的 resolved-member 对力(calcAccJerkPairSingleGroupMember,singles 感受组的成员)与组 CM 侧对力(changeover 权重 k 在两种配对下不对称),或簇 CM 帧漂移 pcm.pos+=pcm.vel*_time_end 用重算前后的陈旧/新值不一致(hard.hpp:1984-2042)。逃逸速度来源定量闭合:∫m×a_注入dt = +85.5 M☉·km/s ≈ 77.6 M☉ × 1.1 km/s。

**Prevention rule**: ① 排查守恒破坏先分类"一次性跳变 vs 连续漂移"——表示切换类跳变(质量置零/恢复)必须用 P_full=real+arti 审计,只看 real 会满眼假泄漏;② gdb 硬件 watchpoint 必须用裸地址(`watch *(double*)addr`),带 this-> 的表达式离开作用域会被自动删除造成"零触发"假阴性;-O1 -g 与 -O3 逐位一致可放心降级调试;③ 逐步 dump(id/mass/vel 全粒子+人工)到文件离线分析是最高性价比方案:方向对齐检验(dP vs m·Δv_cm)一步区分"组被单方加速"与"表示错位";④ 下一步:在泄漏步内 break 于 Hermite 子步尾,比较组 CM 的 Δv 与其扰源邻居 Σm·Δv,或 A/B 补丁 calcAccJerkPairSingleGroupCM 强制对称验证泄漏归零。

### 2026-10-06: 联锁修复(失配停用人工粒子)已实施并验证非主通道——主泄漏窗口组员稳定且 matched;泄漏在 SDAR 积分内部,簇 CM 速度每步跳 ~0.045

**Mistake**: "TT 按设计在成员变化后停用"的门控实际只覆盖张量力输入,且原 group_id 符号翻转联锁已被注释废弃(ar_perturber.hpp:28-50);据此实施的失配联锁(ARTI_MISMATCH_INTERLOCK:失配人工组质量置零)在真实失配事件正确触发(2 次),但主泄漏窗口(47.625-47.72)组为稳定 n2、每步 matched 更新,联锁零触发,|P| 泄漏不变(55.2→143.7)。

**Root cause**(排除链收口):kick 体/createGroup/表示切换/张量 eval(T1 置零逐位不变)/失配陈旧人工粒子全部排除;泄漏发生在 drift 体内部(matted 更新把 AR 视图写入质量承载的人工 CM),CM 分解显示簇 CM 速度 ccm 每大步跳 0.045 pc/Myr(小步 ~0.01,随窗口渐增)——量级 v_orb×Δκ≈0.017 与 slowdown 因子更新同阶,方向恒定;窗口内组员无变化,故与 n2↔n3 churn 无关的每步机制只剩 slowdown 因子更新与 perturber 列表更新两类。--tt-switch 0 仍是最有效缓解(绕过整条人工粒子路径)。

**Prevention rule**: ① 表示级联锁修复方向正确但需先证主通道:实施前先用最小探针(如 CMVIEW 分解 bink/gcm/ccm 三分量)确认泄漏载体;② SDAR 组 cm.vel 的写入者(calcCenterOfMass 重算/AR 内部漂移/initial)可用裸地址 watchpoint 逐个抓——模板头文件行断点(hermite_integrator.h:2694)在 -O1 下可能不绑定,改用函数名 rbreak 或在调用方(非模板)行设断;③ 下一实验:watch SDAR 组 cm.vel 写入点,或对比 perturber 列表逐步差异(NB 数、成员)与泄漏步的相关性;④ ARTI_MISMATCH_INTERLOCK 宏保留(CXXFLAGS += -D ARTI_MISMATCH_INTERLOCK),是正确的防御加固,勿与主泄漏修复混淆。

### 2026-10-06: slowdown/perturber 检查收口——组转换重写动量精确中性(3174 事件 dP≤4e-14,CM 含入);泄漏在 SDorg 尖峰步的积分段内;联锁已按用户决定回退

**Mistake**: 联锁修复虽在真实失配事件正确触发,但主泄漏窗口组员稳定、matched 更新,无改善——按"无改善不改码"回退(质量置零亦有引入双表示/失重的风险)。另:geid 转换探针初版只测簇坐标系动量(恒≈0),不含簇 CM 项,是假排除;补上 `particles.cm.mass*cm.vel` 后才有效。

**Root cause**(最终状态):泄漏与 `kappa_org = kref·pert_in/pert_out` 的尖峰步逐一对应(SDorg 基线 0.0017 ↔ 尖峰 0.02→0.22 递增,阈值 groupedCriterion 1e-2 恰在中间;NB 12↔14 交替 = perturber 列表 churn);kappa(实际应用因子)恒 1。组 form/break 状态重写经 CM 含入审计确认动量中性;泄漏在**转换之间的 AR 积分段**内累积(geid P_pre 序列平滑漂移 +0.1-0.15/步而非跳变)。即:perturber 列表churn → pert_out 骤降 → kappa_org 尖峰(宽对成组翻转,判据 5a62e4e 2026-10-01,晚于原 run)→ 该步 AR 积分(含 slowdown 状态机更新)单向注入簇 CM 速度 +0.018/步。

**Prevention rule**: ① 探针测动量必须含 CM 项(簇/组坐标系动量恒≈0 是结构性盲区,两次踩坑);② 下一步唯一悬置实验:在每个 AR integrateToTime 段前后打印 CM 含入 P + 内层(非根)树节点的 κ 值——注入点应在 slowdown 因子更新路径(syncTreeSlowDownAndDs/calcBinaryTreeSlowDown 之后的状态调整);③ 判据阈值(kappa_org≥1e-2)与 perturber 列表更新均无滞回,churn 环境下必然翻转——修复可考虑双阈值滞回或列表变化平滑,但须先定位积分段内的注入函数;④ Makefile 追加行不触发重编(SDAR 头文件不在依赖),改 Makefile 后必须 make -B;configure 备份要在 configure 之前拷贝(时序错误曾把 plain 当 bse.g 恢复)。

### 2026-10-06: pert_out 振荡的具体原因已查明——星 282 在 r_neighbor_crit 质量放大边界(~0.41 pc)按刷新周期翻转;TT 因果链闭合

**Mistake**: 曾以"kappa 恒 1"与"转换重写中性"推断 TT 无关,忽略了 TT/人工粒子在**上游**驱动扰动列表churn——tt-switch 0 的显著改善必然意味着 TT 在链上(用户指出)。

**Root cause**(NBLIST 逐步实测):perturber 列表成员 12↔14 翻转的实体是**星 282(0.111 M☉,距组 0.3-0.4 pc)以精确的 tree-nstep-mklist=2 刷新周期进出列表**(374 偶发同步);**星 87(0.242 M☉)于 t=47.627 永久离开列表,恰为泄漏窗口起点**。组的 perturber 搜索半径按质量放大 r_out=0.0804×(77.65/0.584)^(1/3)≈0.41 pc,282 正悬边界——每次 createGroup/人工粒子刷新重估 r_neighbor_crit 时被纳入/踢出 → pert_out 骤降 → kappa_org=kref·pin/pout 尖峰(跨 1e-2 判据阈值)→ 尖峰步 AR 积分段注入簇 CM 速度 +0.018/步。TT-off(无人工粒子、成员保真实质量)时邻域判据结构不同,无翻转、无泄漏——因果链与全部 A/B 一致。

**Prevention rule**: ① "机制 X 无关"的推断必须检查其是否在**上游**驱动观测变量(tt-off 改善与"kappa 恒 1"两事实并存的原因:TT 驱动列表churn而非直接施力);② 质量放大 changeover(r∝m^(1/3))使大质量组的邻域边界(~0.4 pc)覆盖大量场星,边界翻转是结构性隐患——perturber 列表与 kappa 判据都无滞回;③ 修复方向排序:(a) 列表/判据滞回(防翻转,治标);(b) 定位尖峰步 AR 积分内的注入函数(治本,最后悬置实验:积分段前后 CM 含入 P + 内层 κ 打印);④ NBLIST 探针已删,复现方法:在 matched 更新处打印 neighbor_address 的 id 列表(NBAdr::Group 取 ->cm.id,Single 取 ->id)。

### 2026-10-06: pert_out 量级反解定位真摄动源;TT 关联定性——churn 双配置皆有,注入只在 TT-on 机器

**Mistake**: 两处归因先后被证伪:282 列表翻转(仅症状,0.11 M☉@0.35 pc 占 pert_out ~1% 不可能造成 12-130× 振荡);87 三体churn(重启实现泄漏窗口 n3 记录为零——那是原始 run 的实现)。

**Root cause**(量级反解):存储 Pout=2.7e6-4.7e7 反解出"0.24-7.8 M☉ 天体在 0.01-0.06 pc",而组自身 perturber 列表全在 1.5-1.9 pc(物理和仅几十)——pert_out 由**步内暂现组**(geid 每步 n_form:1 事件,以 group 型成员带合并质量进入列表)或等价近距离通道主导;公式 m·m/r³ 与 apo 基 pert_in 均无算术错误(pert_in=1454/0.0135³=5.9e8 与实测 6.15e8 吻合)。TT 关联:churn 在 TT on/off 都发生(tt-off 2988 次 n2 事件)但泄漏只在 TT-on——注入需要 TT-on 独有的耦合链(组型 perturber 列表项 / slowdown 状态与人工 CM 同步);TT-off 的 isolated 组无此耦合,churn 周期安全通过。slowdown 状态按内层轨道周期(~4 树步)更新,存储值在更新间陈旧——振荡的离散跳变源于此。

**Prevention rule**: ① 归因反转解:从存储量(Pout/kappa)反推所需的 m/r 组合,与实测列表逐项对照,不一致即存在未识别通道(本轮由此定位暂现组通道);② 探针打印列表须在事件发生时刻(步内)而非 drift 末端——暂现组在末端已解散,末端快照看不见;③ churn 事件计数要区分 run 实现(重启与原 run 的混沌实现不同,事件主体可能不同——87 三体是原 run 的,暂现对是重启的);④ 最后一环验证:在 checkAndAddNeighborGroup 处打印加入的 group 型成员的 (id, m, r),确认暂现组以近距离合并质量进入列表。

### 2026-10-06: 簇 CM 写入源已锁定——calcCenterOfMass@hard.hpp:2042 唯一写入;gdb VLA/-O0 工具链打通

**Mistake**: gdb 追 SDAR 内部状态三轮失败:(1) this-> 表达式 watchpoint 离开作用域被删;(2) 模板头文件行断点不绑定;(3) C99 VLA(hard_int_thread)在 -O1 无调试信息("no such vector element")。resolved 分支假设也被证伪:checkGroupResolve() 恒 resolved(CM 模式被注释,注释自述 "_kappa>3.0 seems cause oscillation")。

**Root cause**(本轮定位):积分期间簇 CM(h4_int.particles.cm.vel)的全部写入唯一来源 = `COMM::ParticleGroup::calcCenterOfMass` @ `driftClusterAndArtificialCMAndWriteBack` hard.hpp:2042(HARD_CHECK_ENERGY 门控的 drift 尾重算);重算前后 CM 速度差即代码内 `dcm_vel`,其动能记账进能量日志 Modify 列(事件时 +115 的指纹)——代码自知重算改变 CM 但只做能量记账。动力学消费者(成员回写/人工 CM)均用重算前的 cm_vel_org,形式自洽;泄漏如何经此进入动量承载表示待 -O0 轨迹复核。另:-O0 与 -O1/-O3 **不逐位一致**(FP 求值序),跨优化级对拍会错位泄漏步。

**Prevention rule**: ① gdb 追 VLA 内部对象必须 -O0 构建;裸地址 watchpoint + 非模板调用行(hard.hpp:3491 lambda 行)布防;COMM::List 计数字段是 num_ 非 size_;栈上 VLA 每 drift 复用同地址,watch 会串簇/串噪声——用 backtrace 的参数(_n_group 等)区分;② 跨优化级重跑会漂移泄漏步,先在同构建里重定位再布防;③ checkGroupResolve 的注释历史("kappa>3.0 oscillation")说明 resolved/CM 表示切换早有振荡前科,值得在修复设计时回看;④ 下一步闭环:SDCORR 加回 -O0 构建重定位泄漏步 → watch 2042 重算前后 dump 成员速度 → 验证 (重算值−cm_vel_org)×m 是否=每步泄漏。

### 2026-10-06: 【结案】真凶=潮汐张量 T3 偶阶项的成员平均非零——AR 每子步净推组 CM 无反作用;修复(扣除张量均值)已验证动量+能量双恢复

**Mistake**: 多轮将泄漏归因于表示切换/失配人工粒子/组转换重写等下游症状;用户坚持"TT perturbation 计算导致积分内动量不守恒"才引入决定性探针(Σm·acc_pert after removal+tensor):正常窗口 <1e-3,泄漏窗口 **29–39 pc/Myr²**,×77.65 M☉×dt=+4.7/树步=泄漏量,333,630 次 AR 子步持续作用。另:首轮 NETP-A 探针误置于成员循环内(测到中途部分和,~0.027 假信号),教训:聚合探针必须在完整循环后。

**Root cause**(数学+测量双确认):calcAccPert 对列表力做 CM 扣除(Σm·a 精确=0)之后加入潮汐张量;fit() 按构造 T1=0("assume input force already remove the c.m.")、T2 线性项 Σm·T2·x=0(CM 系精确抵消),但 **T3 二次项 Σm·T3·x²≠0(x² 为偶函数,±x 成员不抵消)**——张量给组成员留下非零质量加权平均力,无任何反作用通道。近距暂现体主导张量拟合时 T3∝r⁻⁴ 爆炸(0.01-0.06 pc 源 → 净力 11-30 pc/Myr²,实测 29-39 ✓)。kappa_org 尖峰/pert_out 振荡/NB churn 全是同一近距源的伴生症状。TT-off 无张量故无泄漏;drift 尾 calcCenterOfMass 重算只是测量结果。

**Fix**(ar_interaction.hpp,已验证):张量 eval 后同列表力一样扣除其质量加权平均(pot 同步修正)。验证:|P| 55→57.3(修复前→143.7,tt-off→56.7);能量 err_cum 29.9→3.2;TT-on 守恒性与 TT-off 等同。原始 200 Myr run 的逃逸事件(+144 能量/+130 动量、巨双星 1.5-2 km/s 假逃逸、data.core 崩坏)即此机制在无近距暂现体猝发时段的持续累积。

**Prevention rule**: ① 微分/差分化校验必须逐项做**偶宇称检验**:奇阶(线性)项在 CM 系自动抵消,偶阶(常数/二次)项的平均非零——任何"扣除均值"操作之后加入的场展开都要重新扣除均值;② 聚合探针(Σ/均值)必须放完整循环后,放循环内会测到部分和假信号;③ 守恒审计的最短路径是直接测"力通道的净输出"(Σm·a per call)而非层层追状态写入——写入者(calcCenterOfMass)往往只是测量者;④ 修复后验证必须同时查动量(|P| 序列)与能量(err_cum),二者同源时应同时恢复。

**修正语义对照论文复核(2026-10-06 续)**:PeTar_code.pdf §6.4.1 自身定义 A′=A−A(rcm)(测量点先扣 CM 值,零阶项由 CM 粒子单独交付)——但内部动力学需要的严格量是 A′−⟨A′⟩(⟨A′⟩=⟨T3x²⟩≠0 即泄漏),扣均值是把论文定义执行到底而非改变设计。**T3 差分贡献完整保留**:对成员减同一常数使相对力 (a_i−a_j) 严格不变;非等质量双星 x1≠x2(m1/m2 缩放),T3·x1²≠T3·x2² 的差分项(论文测得的轨道修正效应)原封不动;等质量时 T3 项纯均匀,扣除是恒等变换。净分量的物理内容是 O((a/L)²) 有限尺寸净力修正,可忽略且已由 CM 通道(A(rcm))覆盖。反作用(§6.4.2 sampling/伪粒子)与张量均值处理无关,修复恢复一致性。

**修复重构(用户方案,2026-10-06 再续)**:张量 eval 挪入循环 1(列表力之后、acc_pert_cm 累加之前),既有 CM 扣除自动覆盖张量均值,独立修正循环删除——与分离循环版逐位等价(|P| 序列完全一致)。**n_pert==0 分支勘误(用户纠正)**:该分支原本就有张量施加(else 内 SOFT_PERT 块,我此前误读为"完全不施加"),同样无均值扣除——即**两个分支都有同一泄漏**;修正:原块就地加均值扣除(非新增张量),我误加的重复块已删。教训:① python 字符串手术连吞三处结构性括号/宏(每次以不同编译错误显形),多 edits 必须从 git 原版一次性重放并逐锚点断言;② 读分支代码必须读到块尾(EXTERNAL_HARD 嵌套截断了我的视野);③ 未触发路径的错误重构不会被验证场景暴露(n_pert==0 在本测试中不执行)——重构后要静态核对所有分支的施加次数恰好一次。

### 2026-10-06: T5 验证场景落地——回归阈值必须用未修复二进制实测定标;3 体 hard 模式 TT 结构性不触发

**Mistake**: T5 场景初版三处想当然:(1) r3 命令 `--r-search-group 0.03` 触发断言 `r_search_group_over_in<=1.0` 直接 abort(0.03/0.024>1,r_in=r_ratio·r_out);(2) r1(hard 模式)期望 `tt_engaged>=1` 失败——3 体事件构型中 churner 轨道(0.01–0.04 pc)与并组半径必然重叠,系统恒为孤立 AR 组、TT 结构性关闭(r1≡r2 逐位一致本身是有效回归:证明无人工粒子路径不被修复扰动);(3) 初版动量阈值 0.02 用修复版结果(0.0014)直接外推,未实测未修复版——实测未修复仅 0.0102(<0.02,拦不住旧 bug)。

**Root cause**: ① r_search_group 上限是 r_in(非 r_out),树模式缩小 r_out 时必须同步缩 r_search;② hard 模式 TT 人工粒子需要"组外 perturber",3 体尺度重叠构型无此状态(T4 能触发因其尺度干净分离:双星 1.9e-3 ≪ perturber 0.15);③ 泄漏在该 IC 上早期一次近距通过即饱和(t=3 与 t=10 同为 ~0.0103),延时不能增大判别力,判别力只能来自阈值。

**Prevention rule**: ① 新增回归检查的阈值必须双向定标:同一命令分别跑修复/未修复二进制,取"修复版裕度≥3×、未修复版超限≥2×"的阈值(此处定 0.005:修复 0.0014 / 未修复 0.0103);② 定标未修复版的最短路径:`git show HEAD:file > file && make install`,测完从暂存副本恢复再重装(全程 ~2 分钟);③ 泄漏型 bug 的回归指标先测时间饱和性再决定跑多长;④ 场景命令先手动跑通再写入 JSON,断言失败会在 work 目录留下部分输出易误导诊断方向。

### 2026-10-06: 场景 setup 的输入格式与二进制家族强耦合——BSE 中断构建要求 petar.init -s bse(勘误:.g 是 --with-debug=g,非 galpy)

**Mistake**: T5 HTML 报告首跑失败:恢复用户 bse.g 选定后,`petar` 读 setup 生成的 input 直接 abort("Ptcl data reading r_search, id, and group_data fails! requiring 4, only obtain 2")。同一 setup + plain 家族全部通过,误以为 setup 与二进制无关;首版归因还把 .g 后缀误判为 galpy——实为 `--with-debug=g`(Makefile.in:280 `ifeq ($(debug_flag),g)` 添加 .g 后缀)。

**Root cause**: 输入粒子列数是**编译期特性**(interrupt/external 列),与运行时旗标(-b 0)无关:BSE 中断构建的粒子类额外读恒星演化列,petar.init 需 `-s bse`(运行也要 -b 1);galpy 外势构建需 `-t`(外势列+header offsets)。plain 格式 input 只能被 no-interrupt/no-external 家族读。

**Prevention rule**: ① 编写验证场景先确定目标二进制家族,setup 的 petar.init 旗标必须与之一致;纯引力测试选 plain 家族,报告脚本经 `petar.select --optional mpi,omp,avx2` 选择(不带 --require 即 plain 家族)并显式 `--var petar_bin_switch=<选中路径>`——不硬编码绝对路径、不依赖运行前谁被选中;注意 select 会**改变当前选定**(副作用同 T1-T3 管线);② 该 abort 的指纹是"requiring 4, only obtain 2"——列数不匹配即查家族特性列;③ 报告脚本(T2/T3/T5 范式)把 --petar 做成参数,空值走 petar.select,失败时报错而非静默回退(静默回退正是本轮 bse.g 失败被延迟发现的原因)。

### 2026-10-06: T4 m1 在 HEAD 上的既有 abort——r_search 断言的 1-ULP 重构舍入;乘法放缩改精确赋值修复

**Mistake**: T4 管线重启即挂(m1 t≈10.9 abort `pcm.r_search >= member.getRout()`),先入为主怀疑当日 TT 修复/重建引入——用 `git show HEAD:` 回退重建对照后同点 abort,证明是 HEAD 既有问题(Sep 28 旧报告全过 → 之后 f14606f/2fa15b7 的 r_search 不变量提交引入)。静态分析三轮推不出违例(逻辑上 r_out_max 扫描+放缩应恒覆盖),误判方向包括 NaN/能量爆炸(实际能量 1e-9 干净)。

**Root cause**(gdb -O0 帧内实测):BSE 质量流失使 mass_cm 微降 → member r_out 比 pcm 基准值大 1.3e-12 相对量 → r_ratio=1.0000000000013>1 → `updateWithRScale()` 用**乘法重构** `pcm_r_out·ratio`,舍入后比 member 值短 ~4e-15;平时 cushion 是 |v_cm|·dt·factor,但质心系下 v_cm≈6e-13 → cushion≈1e-16 失效 → 断言差 1 ULP 失败。**修法**:ratio>1 时改用两参 `setR(rin·ratio, r_out_max)` 把 r_out **精确赋值**为 r_out_max(无舍入路径),r_search=max(≥0 项+r_out,·) ≥ r_out_max 恒成立。验证:-O0 与 -O3 完整跑通 t=20(161 输出,无断言),T4 全 5 模式恢复。

**Prevention rule**: ① "重构应有值"类断言禁止经乘法/除法往返(`a·(b/a)≠b`,1-ULP 短缺);需要 ≥b 保证时直接赋值 b;② 断言的 cushion 项(速度·dt)可在质心系退化到 0,边界分析必须考虑 cushion=0 情形;③ 恢复旧验证先跑再改:同命令 HEAD 回退对照是排除"新改动引入"的最快手段(~2 分钟);④ abort 在 omp outlined 函数内时,断点帧 locals 要 `frame N`+`info locals` 取(N=bt 中 _omp_fn 帧),-batch 脚本注意 grep "Assertion" 无匹配返回码非零会误报失败。

### 2026-10-07: 模拟配置先查记录再推测——并行/时长/重启三连错被用户纠正;实测基准统一进 performance-log.md

**Mistake**: N1k 验证重启的启动配置连错三轮,均为"推测优先、查证滞后":① 提议 OMP 4 线程(SKILL 明文 N≲10³ 用 1 线程,script-tools 有 N=500 实测 4 线程 +19%,我引用"binaries-rich 2-4"例外却未按 quick rule 先 benchmark);② wall 实测显示"mpirun -n 4 快 8×"推荐 MPI——实为对照探针缺 `UCX_VFS_ENABLE=n` 的混淆(环境变量只设了 MPI 组),直跑+环境变量后差距消失;③ 用短探针 wall 外推时长(0.25 Myr 探针 wall 5 s 里计算仅 1.3 s,启动/IO 占 3/4),先后报出 2 h→52 min→10 min 三个错误估计,真值(prof: 2.6 s/Myr)自始至终可一步获得。另:重启命令最初漏 `-p data.par`(SKILL Restart/resume 节明文"Use -p before overrides"与"-i 0 for binary snapshot"),该节内容我未读即拼命令。

**Root cause**: 并行规模与时长属于"有实测数据可查"的问题,却被当作可直觉推断的问题;wall 计时在短探针上系统性失真(启动+IO+UCX 逐通信开销),且 A/B 对照未控制环境变量;SKILL 相关节(Restart/resume、Parallel Launch Heuristics、UCX 环境要求)都已覆盖,失败在"没有先读对应章节"。

**Prevention rule**: ① 提议任何多线程/多进程配置前,先查 `assets/performance-log.md` 实测记录,无记录则当场用 prof 输出(`Wallclock time per step: Total`)做基准,不得凭直觉;② 时长估算一律用 prof 每步 Total × 总步数(或输出文件 mtime 间隔),禁止短探针 wall 外推;③ A/B 对照必须逐变量控制环境(UCX/OMP/MPI 变量在两组间完全一致);④ 拼重启命令前重读 SKILL "Workflow Patterns → Restart/resume"(-p 前置、-i 匹配快照格式、.par.* 伴随文件);⑤ 每次正式/验证模拟后向 performance-log.md 追加一行实测(N、家族、启动配置、ms/step、s/Myr)——该文件是实测基准唯一权威家,script-tools/SKILL 只留指针。

### 2026-10-07: 全系统动量审计以快照 Σmv 为准(日志 C.M.: 行语义不符);fix47 长程验证通过 + t≈167 暂态 NaN 断言为未决鲁棒性问题

**Mistake**: 验证时先用日志 `C.M.:` 行算 |P|(全系统质心速度×质量),得 39–136 振荡、误判"修复可能失效";三轮对照(未修复 re47/中间版 re47k/新 fix47)该指标不分彼此,才起疑指标本身。

**Root cause**: status.hpp 的 `C.M.:` 行速度与快照 Σmv/M 不一致(fix47.50: 快照 59.5 vs 日志行 113.3),语义非全系统动量。快照直算 |Σmv| 后结论干净:原始 55→132→165→238(单调泄漏),修复版 47–66 有界振荡 153 Myr 无增长;巨双星(id184+377,77.65 M☉)原始 t=200 飞至 243 pc(假逃逸),修复版留在 22 pc 束缚;能量 t=50 处 29.9→5.3。另:fix47 在 t≈167 遇 `driveForOneClusterOMP` 的 `!isnan(pi.vel.x)` 断言 abort——从 fix47.167 快照重启即绕过(重启重建 AR 组/slowdown 内态,不逐位复现),未固定根因,属 hard 积分鲁棒性独立问题,待查。

**Prevention rule**: ① 动量守恒审计一律直接对快照算 |Σ m·v|(`petar.Particle.fromfile(f, offset=petar.HEADER_OFFSET)`,不含 BSE 列时默认参数即可),不信日志 C.M. 行;② 多版本对照若"不分彼此",先怀疑度量再怀疑物理;③ 快照重启不逐位复现崩溃(内态重建),NaN/断言类问题复现要靠 dump 回放(hard.debug)或同进程续跑;④ 长程验证中出现的一次性 NaN 断言要单独立案,不得因"重启绕过"而静默关闭。

### 2026-10-07: solver 跨运行位级不可复现是地址依赖(ASLR)所致——调试/验证 run 必须 setarch -R;NaN 为实现系综中的稀有事件

**Mistake**: fix47 在 t≈167 NaN abort 后,先后用"快照 166/167 重启"与"data.47 完整重跑"试图复现,全部干净通过,一度归因于"重启重建内态不可复现"并怀疑内存问题。完整重跑(同快照/同 par/同二进制/同 seed/单线程)在第一个输出 t=48 即与 fix47 字节级不同(C.M. 打印仍相同——差异在打印精度下),证明**同配置重跑本来就不同**。

**Root cause**: 判别实验:`setarch $(uname -m) -R`(关 ASLR)下连跑两次输出**逐位一致**;ASLR 开时每次进程堆布局不同→某处存在地址依赖行为(指针序/布局敏感路径,值键容器与值排序已排除)→FP 求和序等微差被混沌放大。因此每个进程是一次"实现采样",NaN 是该实现系综的稀有事件(fix47 命中、repfull 未命中);不可复现是常态而非异常。内存泄漏假设排除(泄漏不产生 FP 分歧);经典越界损坏未被证实(地址依赖在固定布局下是确定性的)。

**Prevention rule**: ① 需要复现/对拍/调试的 run 一律加 `setarch $(uname -m) -R` 前缀(与 ASan 构建并用时本就是仓库惯例);② 采样可复现实现系综:关 ASLR + 改变环境变量填充(如加无效 env pad)平移初始栈→不同布局且可复现;③ 判定"是否位级可复现"用输出快照 `cmp` 而非日志打印(打印精度会掩盖分歧);④ "重启不复现"不能作为排除 bug 的证据——完整重跑不可复现时先做 setarch -R 对照;⑤ 下次 NaN 出现时:先打印失败簇成员(id/质量/间距)再 abort 的一行诊断 + 上述可复现采样,即可闭环定位。

### 2026-10-07: NaN 已定位——致密核紧双星形成瞬间 Pert_In~1.8e9(摄动源近零距/重合类),form-break-重构震荡中 AR 发散;可复现采样协议生效

**Mistake**: NaN 一度按"稀有随机事件"处理;env-padding 采样显示 t=0 起跑 8 实现中 5 个死在同一 t=79(去重 5/6)——并非稀有,而是几乎每实现必遇的物理阶段(大质量星内落成紧双星)触发。

**Root cause**(diag1 t=57.7012,打印+dump 回放链):致密核形成紧双星(a=0.0055 pc,P=0.0043 Myr);其 Pert_In=**1.85e9**、Pert_Out=5.1e5——摄动强度超 AR 有效性 9 个量级,指向摄动源与组几乎重合(r→0,同 2026-09-20"停放粒子重合→0/0"家族);组 0.0007 Myr 内 form→break→3 体重构(semi<0,ecc=1.0047),窗口内 41 次 large-energy;最终 AR 积分 NaN(id1),NaN 经树力一步毒化全系统(检测点在下一漂移步:全表 vel NaN、binID=0、位置正常)。t=79/57/167 是不同实现的同一机制首次触发时刻。

**Prevention rule**: ① NaN/全局毒化的诊断顺序:漂移步检测点打印全表状态 → 定位首个 NaN 粒子 → 取最近 dump 回放(hard.debug)看 form/break 序列与 Pert_In/Out;② Pert_In>1e6 即"重合/近零距摄动"红旗,值得在 perturber 收集处加距离下限断言;③ 采样协议(setarch -R + env padding)把"不可复现 NaN"变成可复现枚举,是此类问题的标准工具;④ 下一步闭环:在组初始化处打印 perturber 到组 CM 距离分布,验证重合粒子来源(人工粒子/CM 表示/真实恒星),对症加最小距离保护或重合消解。

### 2026-10-08: NaN 根因闭环——轻星穿越巨紧双星内部(r~7.6e-4 pc < 双星间距),有限强时变摄动使 AR 时间变换积分内部发散;非重合粒子、非摄动力溢出

**Mistake**: 依 Pert_In=1.85e9 与 2026-09-20 先例推测"摄动源与组近零距重合"(r→0 除零),在 calcAccPert 加 r²<1e-12 诊断——6 个全 NaN 实现零命中,假设被否定;gdb 捕获路线因实现彩票(+gdb LINES/COLUMNS 环境泄漏、输出前缀长度参与布局抽样)效率过低放弃。

**Root cause**(双阈值诊断,r<0.01 距离层命中 + |acc_pert|>1e60 层零命中):0.132 M☉ 星(id 848)从 77.65 M☉ 紧双星(184+377,a~0.005 pc)**内部穿过**(距成员 7.6e-4–2.7e-3 pc,eps_sq=0);摄动力量级大但**始终有限**(changeover 加权在),NaN 产生于 **AR 时间变换辛积分器内部**——紧双星 slowdown 饱和 + 组成员间强时变摄动超出"弱摄动 Kepler"适用域,form-break-重构震荡即吸收失败的表现,最终一步发散。每实现必经的物理阶段(大质量星内落成紧双星后遭遇穿越)→ 触发率高(同窗 6/6)。

**Prevention rule**: ① 双阈值诊断(几何距离层 + 物理量溢出层)可一次实验同时证实/证伪两类假设——本例距离层命中+溢出层零命中直接把责任从"力计算"移到"积分器";② gdb 批处理调试需 `unset environment LINES/COLUMNS`(gdb 向 inferior 泄漏自身终端变量→改变实现),且输出前缀字符串长度也参与地址布局抽样;③ 修复方向候选(待设计):(a) 摄动者进入组半径内即强制吸收为组成员(n→3 精确积分)而非保持 perturber;(b) 对 r<r_in 的 perturber 力加强 changeover 截断;(c) 强摄动紧双星的 AR 步长/ slowdown 上限保护。

### 2026-10-08: 方案(a)检测层吸收未达(持久组旁路+ap_manager 记账),NaN 归责修正为簇 Hermite 单粒子;多轮"零命中"系陈旧二进制假象——诊断前必须 strings 验证

**Mistake**: 三重流程事故:① 强制吸收补丁多次构建因 Makefile 家族切换(恢复 bse.g 后 plain 目标无规则)与头文件不被 make 跟踪而**静默使用陈旧二进制**——fa/tt/herm/vel 系列"零命中"部分是假象(直到 strings|grep BLAME=0 才暴露);② 吸收真正生效时触发 ap_manager 成员计数一致性断言(吸收未向 artificial-particle 记账注册);③ NaN 机制叙述(AR/摄动/张量通道)全部基于陈旧二进制的"安静"证据,归责错误。

**Root cause**(strings 验证后的可靠证据):BLAME 诊断:首个 NaN 粒子 id=1 m=1.83 **type=single**(簇 Hermite 域);PERT-DIAG(848 穿越双星,r~7.6e-4,acc 有限)与 TT/GT/VEL/HERM 诊断均来自混合新旧构建,仅 -B+strings 验证过的可信:AR 子步无 dt 奇点、无近距对(r<1e-4)。责任修正:**簇 Hermite 对单粒子的 predictor-corrector 在极端有限力梯度下失稳**(动力学诱因仍可能是 848 穿越巨双星造成的环境)。

**Prevention rule**: ① 任何"诊断零命中"结论前必须 `strings <binary> | grep <DIAG标记>` 验证诊断确实在二进制里,再跑实验;② make 家族切换(configure/restore)后立即核对目标存在(`make -n`),头文件改动一律 -B;③ 强制吸收类改动必须同步 ap_manager 记账(createArtificialParticles 的成员注册路径),否则触发 hard.hpp 成员计数断言;④ 下一步:(i) 修近邻打印的 NaN 比较陷阱(dr2 非有限时按 |pos| 排序或跳过 NaN offender 自身),重跑 BLAME 拿单粒子的近邻归责;(ii) 对簇 Hermite 加 acc/jerk 量级诊断(|acc|>1e6 时打印对 id/距离)定位失稳对;(iii) 吸收方案重做时走 ap_manager 正规注册。

### 2026-10-08: 用户指路的 SDAR 侧判据改造(成员选择性保组)已实现但未命中 NaN;力求值通道全排除,嫌疑收窄到速度更新通道(jerk/块步长/changeover 修正)

**Mistake**: 方案(a)放 hard.hpp 检测层是错的方向(用户纠正):① 持久组旁路检测;② 吸收需 ap_manager 记账。SDAR 侧积分时调整机制里**吸收入口本就存在**(checkNewGroup 的 group case 用 perturber.r_min + groupedCriterion),真正的判据缺口在 **break 侧只评估根对**:848(双曲线)远离时根对判散→**整组解散**(含内层紧双星)→重组窗口+震荡。

**Root cause(本轮实验)**: 在 hermite_integrator.h 实现 anyInnerNodeGroupedIter(任一内层树节点过判据即保组)+ break 调用点接入——6 实现电池仍 6/6 NaN,form/break 震荡假设被否定(或非唯一通道)。经 strings 验证的力通道全面排查(单粒子 |acc|<1e6、无 r<1e-2 近对、AR 子步有限、TT/张量有界、κ 有钳制)全部干净,而 NaN 首发恒为**簇域单粒子**——生成点只能在速度更新通道:Hermite 校正器 jerk 项、Aarseth 块时间步、PeTar 侧 changeover 力修正(soft↔hard 混合,Force_corr 通道,本 session ULP bug 即在其邻域)、或 AR↔Hermite 写回。

**Prevention rule**: ① 下一步诊断目标:单粒子**每步 Δv** 打印(>1e2 即报,附 acc/jerk/dt)+ jerk 幅值——Δv 直接锁定爆掉的那一次更新及其通道;② changeover 力修正(forceCorr/Force_corr)是高嫌疑(历史上 ULP bug 同域),PeTar 侧 hard.hpp 力修正路径需同样排查;③ 判据侧改动(成员选择性保组)语义合理但未验证收益,去留待定(与诊断同树未提交);④ 复用本轮基建:strings 验证 + setarch -R + env padding 电池是标准 NaN 狩猎流水线。

### 2026-10-08: 【NaN 闭环】毒源=AR 积分内 TT 采样粒子位置 NaN(成员全程干净)→ 写回 soft 毒化树 → 全体单粒子 acc NaN → 全局 vel NaN;用户"入口遗留"检查是破局关键

**Mistake**: 多轮排查都盯着"力通道"(摄动力/张量值/Hermite acc/AR dt)——全部安静;直到按用户提示做**相位扫描**(硬积分入口/出口)发现两端口真实粒子都干净,才把检查移到 kickOne:soft acc 已非有限(1936 处扩散,首发即 BLAME 罪犯 id=1)→ 树被毒化 → 毒源必在喂树的**soft 系统人工粒子**:出口扫描人工粒子块,抓到 arti j=8/12(TT 采样,块布局 4TT+1orb+1CM)pos=(-nan,-nan,-nan) 而 vel 有限。

**Root cause(完整链)**:某组的**潮汐张量采样粒子在 AR 积分内位置变 NaN**(TT 采样点由内层双星 Kepler 解放置;848 俯冲驱动的强扰动使内层轨道根数退化(ecc→1/近心距→0)时 Kepler 放置产生 NaN)→ 写回 soft → 树矩 NaN → 所有单粒子 acc NaN → kick 全局毒化 vel → 漂移步断言 abort。这解释了此前一切"安静":真实成员/AR 成员检查干净(毒在人工粒子几何),力值检查干净(毒在位置不在力)。修复方向:TT 采样放置的退化保护(Kepler 解非有限时保持上一位置/钳制)+ 上游内层轨道退化处理(与 break-guard 方向同源)。

**Prevention rule**: ① 相位扫描法(每个相位边界全表查非有限)应从第一次就用,且**必须覆盖人工粒子**(soft 系统全量,不只 ptcl_hard_)——"谁喂树谁可能是毒源";② 力值诊断的阴性结果只排除力通道,几何(pos)通道需要独立检查;③ SDAR break-guard(成员选择性保组)实验未命中本 bug 已回退,但其方向(内层退化保护)与本根因同族,重做时可复用;④ 下一轮:instrument tidal_tensor.hpp 采样放置函数,捕获退化 Kepler 输入(打印根数),然后加保护并跑 6 实现电池+T4/T5 回归。

### 2026-10-08: 【NaN 结案】轨道采样粒子对双曲根数无解→NaN 毒树;最终修复=open 组在写回处 continue(人工粒子标记变更、真实重建)+updateArtificialParticles 入口 ASSERT

**Root cause(终版确认)**:NaN 生成点是 `orbit_manager.createSampleParticles`(由 `updateArtificialParticles` 调用):组根为双曲线(semi<0,ecc>1)时 `orbitToParticle/calcMeanAnomaly` 无椭圆采样解→采样粒子 pos NaN→写回 soft→树矩 NaN→全体单粒子 acc NaN→全局 vel NaN。入口确认:1901 孤立路径有 `semi>0` 守卫,**2089(driftClusterAndArtificialCMAndWriteBack)缺失**——故障全部经此进入。逐层证据:TT 生成输入正常(TTGEN 零)、人工粒子入口干净(ENTRY 零)、出口 NaN(EXIT 命中 m=25.93=77.65/3 轨道采样)、"after create" 命中(TTUPD,semi=-0.046)。

**Fix(用户设计,已验证)**:hard.hpp 写回循环中 `if (bink.semi <= 0.0) continue;`(声明前移修复编译序)——跳过刷新+成员质量备份,且 `group_arti_update_list` 保持 false → 计入 `sdar_n_groups_arti_change` → 人工粒子集**标记变更、下次调整真实重建**(非保留陈旧);artificial_particles.hpp `updateArtificialParticles` 内加 `ASSERT(_bin.semi>0)` 固化调用契约。验证:NaN 窗口电池 6/6(修复前 6/6 必死)、T5 回归 9/9、ASSERT 零触发、全程跑穿 t=10–15 暴力相(0.59 单事件 dE 与对照变体逐位同——固有物理,两方案一致;对照变体 6/6 CLEAN 200)。中间方案反证:仅跳过刷新但保留人工粒子为"当前"(guard-wrap)在 t≈15 因 n63 巨簇能量灾难死亡——open 组的陈旧人工表示不能复用,必须真实删除。

**Prevention rule**: ① open(双曲)系统的椭圆类采样/开普勒类计算一律先判据;人工粒子生命周期变更必须走正式删除/重建通道,不能"跳过更新但保留状态";② 调用契约用入口 ASSERT 固化,调用点守卫保持与既有范式一致(1901/2089 对称);③ 重负载多实现验证避免 9p 盘(localdata=E:\):dump 洪泛阶段 37 s/Myr(15× 慢,wchan=p9_client_rpc),改用 ext4 目录;④ 单事件 dE/Etot~0.6 的暴力相是修复后大质量核心固有活动,评估方案优劣要看事件特征跨方案一致性而非绝对值。

### 2026-10-09: 【百 dump 洪泛滥根因结案】形成瞬间成员 changeover 被组 CM 质量换算(×m^(1/3))改写→calcEnergy 对势权重全变→Epot 记账跳 9.34;纯监视器假阳性,动力学无恙

**Mistake**: 归因两轮错误:先"AR 积分发散"(dE 形成后冻结排除)、再"SD 替换注入"(gdb 对账真/SD 四能量逐位相等推翻)。

**Root cause**(n100 dump 回放 + gdb 逐层对账,铁证):近抛物线双曲对 {m4.6, m46} 形成时,`collectGroupMemberAdrAndSetMemberParametersIter` 把成员 changeover **替换为组 CM 的**(质量换算,m^(1/3):轻成员 r_in/r_out ×2.22 跃变,gdb 实测 0.0160/0.160→0.0356/0.356);`HermiteInteraction::calcEnergy` 的每对势能权重 = `calcPotWTwo(chi_i, chi_j, r)` — 物理位置不变、权重全变 → 簇 Epot 跳 +9.34、Ekin +0.31,监视器对冻结的旧参考报 dE=9.65 并**永久带偏**。量纲论证封死"真实积分误差"可能:窗口内位移 1.2e-4 pc,产生 9.34 Epot 需 46 M☉ 星在 0.001 pc 内(实际 0.026,差 4 个量级)。力侧有配套修正(correctForceChangeOverUpdate),**能量参考侧无重锚定** — 缺口即此。

**Prevention rule**: ① 任何"表示切换改 changeover/质量"的事件,力与能量两侧必须对称处理:力侧修正已有,**能量参考需同步重锚定或把权重差记账进 Modify 列**(与质量变化 de_change_modify 同范式);② 能量误差归因三问:形成/打断瞬间?(查表示切换记账)· 冻结不涨?(查参考偏置)· 量纲可能?(ΔE vs Δr·∇Φ);③ gdb 断点命令里调方法会返回垃圾值(getRin() 等),读私有字段(r_in_)可靠;VLA(ptmp)与 vec 索引同样不可靠,用 .x/.y/.z 或绕开;④ dump 洪泛≠积分病:先指纹去重统计(88% 不同指纹=成员 churn),再抽单 dump 回放对账,一步分辨真误差/假阳性。

### 2026-10-09: 【记账修复落地】changeover 切换差延迟入账(检查态精确求值,ID 寻址)——dump 洪泛 366→5,T5 9/9;两轮实现各踩一坑(位置漂移过记/陈旧索引 NaN)

**Mistake**: 实现两轮失败:① 即时记账(切换时刻求值)过记 1.16 — 双求值实验证明精确表示差 9.6526 需在**监视器同状态**(singles 预测位/成员写回位)下求值,切换时刻位置漂移经权重过渡带放大;② 索引版首跑 NaN,误诊为"adjustGroups 重排数组"改用 ID 寻址——用户质疑后实测(4737 条记账,索引移动 0、离开 0:SDAR 步中不重排 particles 数组),判别实验证明真因是**重合位置对 r=0 除零**,r>0 防护才是真修复,ID 寻址冗余已回退索引版。

**Fix(最终版)**:sync 只记(ID, 旧chi) 快照;三个 calcEnergySlowDown 评估点后 consumeBookedChangeoverDE():按 ID 重解索引(离簇跳过,r≤0 防护),同状态算新旧权重势差,经 SDAR 新公共接口 addDEChangeModifySingle 入 de_cum/de_modify_single(参考自动跟随);再调一次 calcEnergySlowDown 刷新参考。验证:n100 dump 回放 14→**0** 事件;6 实现电池 6/6 干净且 dump **366→5**(残 5 个为 n10-16 小系统真实 ~1e-2 误差,监视器恢复判别力);T5 回归 9/9。

**Prevention rule**: ① 表示切换类记账必须在**消费点的状态**下求值(与监视器同帧),切换时刻求值会引入过渡带位置漂移误差;② 除零防护(r>0)是重合位置对的硬性要求;SDAR 步中 adjust 不重排/不删除 particles 数组(实测),数组索引跨迭代安全——但下结论前先做双因子分离实验,双改动同时上会混淆归因(本轮 ID+r防护 即教训);③ 消费型账目配对设计:记录点最小快照+消费点完整求值,中间状态变化由 ID 解析吸收;④ 修复后残差事件要做性质抽检(确认是真实误差而非新假阳性)。

### 2026-10-09: 【设计回退结案】力路径审计证伪同步收益→记账修复整体废弃,回退为"成员终身保持自身质量 changeover"最简设计;6/6 电池 dump=0

**Mistake**: 记账修复(上条)方向正确但治标;用户追问"同步设计本身是否必要、pcm 微扰失配有无实证"后,逐力路径审计发现前提不成立。

**Root cause**: 全部力路径(AR calcAccPert、Hermite PairSingleSingle/PairSingleGroupMember、PairSingleGroupCM 逐成员、PairGroupCMSingle)用的都是**成员自身** changeover,CM 模式权重取 max(CM.rout, 成员.rout)=CM 的更大值——把成员 changeover 同步为 CM 值对力精度**零收益**,只制造形成/打断瞬间的权重跳变(监视器假阳性源头)与配套补丁链(打断重缩放×2、r_search 拷贝、记账)。此外 petar.hpp 的 CM 模式 r_search 设置处本就把 CM.rout 撑到成员最大值(上游固有),断言 pcm.r_search>=成员.r_out 恒成立——最简设计与之自洽。

**Final design**: ① 删除形成时赋值(collectGroupMemberAdrAndSetMemberParametersIter 不再改成员 changeover);② 删除两处 ex-member 重缩放(driveForOneClusterOMP 写回、findGroups 重置循环);③ syncMemberChangeoverScale 只拷 r_search 且加 max 下限 `std::max(pcm.r_search, 成员.r_out)`(c.m. 速度近零时 r_search 可小于成员 r_out);④ 写回路径 Ptcl::calcRSearch 天然含 +r_out 项、setRNeighbor/轨道粒子处各有下限,不变量全覆盖;⑤ SDAR 侧 addDEChangeModifySingle 公共接口一并移除。副作用:i_cluster_changeover_update_/changeover_update_flag 机构成为死代码(标志永不置真)——留作防御,因发送列表侧 MPI 轨道粒子修正仍活跃。验证:6 实现电池 6/6 干净 dump=**0**(优于记账版 5);T5 9/9;n100 dump 回放 dE=-1e-9 无事件;功能冒烟通过。

**Prevention rule**: ① 加补偿层前先审计补偿对象的**全部消费路径**——若力路径根本不用被同步的量,同步本身就是噪声源;② "配套修改链"(2fa15b7→f14606f)回退时必须整链回退,靠 git log 找齐依赖 commit;③ 最小设计优于记账修补:消除跳变源 vs 为跳变记账,前者使 dump 从 5→0 且删代码;④ 回退上游固有行为(形成时赋值是 master 既有代码)前,先 `git show <commit>^` 确认归属,避免误删上游机制。
