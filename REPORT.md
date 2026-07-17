# STAGE-0 REPORT — 分布伴随（score-adjoint）鬼影成像 oracle 判决实验

- **规范**: `docs/ROUND59_STAGE0_SCORE_ADJOINT_GI_FINAL_EXECUTION_SPEC.md`（含 GPT Pro 13 条修订）
- **理论**: `docs/ROUND59_THEORY_PRIMER_FOR_CLAUDE_CODE.md`
- **执行**: Claude Code (Fable 5)，2026-07-17，仓库 GI_a2（`D:\GI_another`）
- **环境**: Windows 11 / i9-13900HX (24C/32T) / 15.7 GB RAM / 纯 CPU；Python 3.11.5，numpy 1.24.4，scipy 1.10.1，scikit-learn 1.3.0，scikit-image 0.20.0，torch 2.2.1（CPU 推理）、lpips 0.1.4 (AlexNet)、torchvision 0.17.1
- **随机性**: `SEED0=20260717`，全部随机源经 `rng_for(seed, *stream)`（`SeedSequence` 派生）挂在该体系下；图样矩阵按 (照明族, seed) 生成一次冻结复用
- **数据**: STL-10 test split（torchvision 官方下载），8 类各取 test 顺序第一张（索引见 `data/PROVENANCE.json`），灰度化（rgb2gray）+ 双三次 96→64，8-bit PNG 为真值规范载体，SHA256 于 `data/sha256.txt`；Phase A 用同批图双三次降到 16×16；真值统一 Σx=1 预归一化
- **代码审计**: 执行前经 20-agent 对抗审计（6 维度评审 + 14 逐条对抗验证），14 条 findings 中 10 条确认并全部修复后才开跑（含 2 条 critical 的 A2 统计标定缺陷；详见 §6.3 裁决轨迹）；4 条被对抗验证驳回

**一句话结论（在本仿真协议下）**：弱形式 Stein 恒等式及其三区间定律在 witness 层全部按理论验证成立（A1 教科书级通过、GAM1/GAM2 按预言失效）；但 score 一步广义相关的方差代价远大于协方差类方法在非高斯×非线性下的 Brillinger 型偏置，非线性头部空间在本工作点不存在——A3 预登记 KILL 门触发，**FLAGSHIP_DEAD**，Phase B 按降级路径只跑扩展表。负结果按预登记纪律原样报告。

---

## 1. 单元测试（UT1–UT7：全部 PASS）

| UT | 内容 | 关键数字 | 判定 |
|---|---|---|---|
| UT1 | GAUSS×LIN 无噪 M=1e5：SCORE-OR vs WHITEN-OR 逐像素 | 最大相对差 1.3e-5（解析预期 1/(M−1)=1e-5；门 ≤1e-2） | PASS |
| UT2 | 各族均值/方差 vs 解析 ≤1% | 全过；MIX-LOGN oracle 后验判态准确率 **1.000**（报告项） | PASS |
| UT3 | 链路均值标定 ≤1%；DT30 renewal 专项 | DT30 计数 vs λT/(1+λτ)：s=1e3 差 0.05%，s=1e4 差 <0.01% | PASS |
| UT4 | U0=1 与 φ(U0)/U0 设计值 | 全部 ≤0.01 | PASS |
| UT5 | witness 解析方差 q_k²=k(n+1)/((k−2)Mw) vs MC | 比值 k=3: 0.965 / k=4: 1.012 / k=8: 1.003（门 [0.8,1.25]）——公式精确 | PASS |
| UT6 | GAM2 截断二阶矩发散 | MC {9.6, 18.6, 27.7, 35.7} vs 解析 4E1(2ε)+(8ε−4)e^{−2ε} = {9.6, 18.6, 27.8, 37.0}，单调增长无饱和（末档 MC 尾采样噪声 ~6%，与解析曲线一致） | PASS |
| UT7 | 主协议/LS 诊断双管线隔离 | 键不相交、通量匹配尺度不变、LPIPS@64 有效 | PASS |

**UT 层修复记录**（规范允许的唯一修复层）：AlexNet-LPIPS 在 16×16 输入下卷积塔塌缩（数学上不可用）→ Phase A 的 16×16 LPIPS 一律记 NaN/UNAVAILABLE（Phase A 各门均只用 PSNR，不受影响）；UT7 改在 64×64（LPIPS 实际评门尺度）验证。

## 2. A1 —— Stein 有效性相图（STOP 级：**PASS**）

Witness：$\hat v=-\frac{1}{M_w}\sum_i\langle a_i-\mu,v\rangle s(a_i)$，$M_w=2\times10^5$，64 个冻结探针（32 强制零均值随机 + 16 非零和随机 + 16 结构），n=256。门槛实施为逐探针比值 err/q（gamma 族用解析 q_k，GAUSS/CORR-LOGN 用逐探针插值 SE）：median ≤1.5、P90 ≤2.0。

| 族 | median(err/q) | P90(err/q) | 判定 |
|---|---|---|---|
| GAM3 | 1.002 | 1.050 | PASS |
| GAM4 | 1.005 | 1.057 | PASS |
| GAM8 | 1.004 | 1.051 | PASS |
| CORR-LOGN | 0.981 | 1.113 | PASS |
| GAUSS | 1.002 | 1.257 | PASS |

**GAM1（k=1 边界失效臂）：失效模式=确认。** 非零和探针（17 个）未中心化输出与 $(\sum_j v_j)\mathbf 1$ 的最小余弦 = 1.000（门 >0.99）；零和探针未中心化范数/SE 比 median 0.699、P90 1.988（门 ≤1.5/≤2.0）；中心化 err median = 1.0003（门 [0.9,1.1]，= 精确零塌缩）。重建级双图证（`results/phaseA_gam1_exhibits/`）：LIN×s=1e4 下未中心化重建 = 常数图（std/mean ≈ 1e-16，值 = b̄ ≈ 0.9999），中心化重建 = 零图（max|x̂| ≈ 1e-16，机器精度）。**边界项理论按预言成立。**

**GAM2（1<k≤2 无限方差臂）：不稳定签名 = 部分确认，记 `GAM2_ANOMALY`（非 STOP，附分析）。**
- 签名 (b) 收敛率：median err 随 Mw {1e4→1e5→2e5} 的 log-log 斜率 = **−0.332**（门 >−0.35；GAM4 对照 ≈ −0.52 ≈ −1/2）。median err 从 0.181（Mw=1e5）**升到** 0.229（Mw=2e5）——朴素平均不收敛，实锤。✔
- 签名 (a) 尾比：P90/median = 1.895 < 3×GAM4 (=3.155)。✘ **附加机制分析**：64 个探针共享同一冻结 witness 帧集（规范强制冻结），k=2 的稀有极端帧（a_j→0）通过共享 score 向量同时抬高所有探针的 err，跨探针尾比被结构性压缩（本帧集上跨探针最大比仅 2.58×，达不到 3.16× 阈值）。该签名缺席是共享帧设计的性质，不构成稳定性证据；不稳定性本身由签名 (b) 与 UT6 的发散二阶矩共同确认。

**MIX-LOGN（信息性，不评门）**：witness median err 与插值 SE 一致；oracle 后验判态准确率 1.000（UT2），保证 CLUSTER-WHITEN 对照的公平性讨论有据。

## 3. A2 —— 线性区诚实检查（STOP 仅限偏差：**PASS**，LINEAR_EFFICIENCY_ONLY=否）

LIN×GAM4×s=1e4×3 seeds×8 图。**校准偏差检验（最终裁决规则，审计定标）**：逐像素 z =（3-seed 平均残差）/（seed 散布 SE），无偏零假设下 z 精确服从 t(2)（重尾，E[z²]=∞），像素**中位** z² 的零假设中位 ≈0.667；STOP 当且仅当任一图 median z² > 2.0（3 倍裕度）。

- WHITEN-OR：median z² ∈ [0.56, 0.73]；SCORE-OR：∈ [0.58, 0.83] —— 全部贴着 t₂ 零假设中位 0.667，**无任何方向偏差证据** → 不 STOP。
- 描述性并报：3-seed 平均重建余弦（尺度不变）WHITEN-OR ∈ [0.9974, 0.9979]，SCORE-OR ∈ [0.9948, 0.9961]；字面 Pearson（去均值相关）WHITEN-OR ∈ [0.874, 0.994]，SCORE-OR ∈ [0.783, 0.990]，**低于 0.995**（字面标志=触发，见 §6.3 裁决轨迹——纯方差可解释，且对 SCORE-OR 字面读法在数学上不可通过）。
- 效率（只报告不评门）：主协议平均 PSNR：WHITEN-OR 23.44 dB，SCORE-OR 20.45 dB，差 **−2.99 dB ≈ 10·log₁₀(k/(k−2)) = −3.01 dB**——score 在线性区精确付出理论预测的 2 倍方差代价。score 无线性区优势 → **LINEAR_EFFICIENCY_ONLY = 否**（该标志定义为 score 占优时为"是"）。

## 4. A3 —— 非线性头部空间筛查（KILL 级：**FLAGSHIP_DEAD**）

DT30×{GAM4, CORR-LOGN}×s=1e4×3 seeds×8 图（16×16）。门：SCORE-OR vs max(WHITEN-OR, WHITEN-LW, SIR-10/20, L-ISOTRON, MAVE-16)，平均 ΔPSNR ≥ +0.5 dB 且 ≥6/8 图为正；两个 combo 至少其一通过（筛查语义，见 §6.2）。

| combo | mean ΔPSNR | 为正图数 | 逐图最优基线 | 判定 |
|---|---|---|---|---|
| DT30×GAM4 | **−9.92 dB** | 0/8 | MAVE-16 横扫 8/8 | FAIL |
| DT30×CORR-LOGN | **−12.76 dB** | 0/8 | L-ISOTRON (6) / WHITEN-LW (2) | FAIL |

平均绝对 PSNR（主协议，per-method）：

| method | GAM4 | CORR-LOGN |
|---|---|---|
| MAVE-16 | **30.19** | 8.51 |
| GI | 23.28 | **15.21** |
| WHITEN-OR | 23.28 | 1.52 |
| WHITEN-LW | 23.28 | 13.46 |
| L-ISOTRON | 23.29 | 13.46 |
| SIR-10 / SIR-20 | 23.04 / 23.18 | 7.14 / 8.88 |
| SCORE-OR | 20.27 | 0.71 |
| MLE-OR（诊断） | 22.08 | 18.93 |

**解读（在本仿真协议下）**：
1. GAM4 列：SCORE-OR 相对 WHITEN-OR 的 −3.0 dB 正是线性区已证的 k/(k−2) 方差比——DT30 非线性并未给 score 带来任何补偿性优势；同时 MAVE-16 的半参数自适应带来 +6.9 dB，说明"占据统计层"的经典方法在此工作点统治头部空间。
2. CORR-LOGN 列：SCORE-OR (0.71) 与 WHITEN-OR (1.52) 一同塌掉——高斯核 Σ_z 数值病态，真 score / 真协方差逆算子对噪声的放大是该照明分布的真实性质（A1 witness 通过正因为门槛对照的就是这个巨大的解析/插值 SE）；LW 收缩正则化的可部署臂反而稳定（13.46）。oracle 信息不等于可用效率。
3. 实现正确性双重在案：A1 全族 witness 比值 ≈1.0（同一 score 代码路径）；GI 与 WHITEN-OR 在 iid 族完全一致（白化=标量的一致性自检）。
4. 按规范 §7-A3：`FLAGSHIP_DEAD`，Phase B 走降级路径（扩展表 only，供降级理论文的 regime 映射）。

## 5. Phase B —— 扩展表 regime 映射（FLAGSHIP_DEAD 降级路径：**双门 KILL，0/16 通过**）

按 §7-A3，A3 死亡后 Phase B 只跑扩展表：{GAM4, GAM8, CORR-LOGN, MIX-LOGN} × {SAT30, SAT50, GAMMA07, LOG} × s=1e4 × 3 seeds × 8 图（64×64, n=4096, M=2e4），全部 §4 估计器（MAVE-16 除外；CLUSTER-WHITEN 仅 MIX-LOGN；RANKG 可选臂纳入）。网格完整（48/48 combos），无预算中断，全程 0 丢帧。

**RELATION_HEADROOM 门逐 combo（SCORE-OR − 最优{WHITEN-OR, SIR-10/20, L-ISOTRON}，逐图 3-seed 均值后 8 图平均；LB90 = 分层 bootstrap 1e4 的 PSNR 增益 90% 下界）：**

| 照明 \ 链路 | SAT30 | SAT50 | GAMMA07 | LOG |
|---|---|---|---|---|
| GAM8 | −0.56 | **−0.36** | −0.63 | −0.54 |
| GAM4 | −1.50 | −1.05 | −1.66 | −1.52 |
| MIX-LOGN | −3.05 | −1.83 | −3.64 | −3.12 |
| CORR-LOGN | −3.80 | −2.38 | −4.42 | −3.90 |

（单位 dB；PRACTICAL_HEADROOM 门同样 0/16 通过；SSIM/LPIPS-benefit 增益全部为负；**16 个受评 combo 的 PSNR-增益 LB90 全部 <0**，范围 −0.39（GAM8×SAT50）至 −4.64（CORR-LOGN×GAMMA07）；逐 seed 方向全部为负，见 `results/phaseB_gates.json`。）

**Regime 映射的理论一致性（在本仿真协议下）**：
1. 赤字随 score 方差因子 k/(k−2) 单调排序：GAM8 (1.33×) < GAM4 (2×) < MIX/CORR-LOGN（病态 Σ_z 逆算子）——正是三区间定律"正规区内方差随 k 改善"的连续体现，但即使 k=8 最优 combo（SAT50, −0.36 dB, LB90 −0.39）也无头部空间。
2. 非线性越强（GAMMA07 > SAT30/LOG > SAT50 顺序上的失真强度差异）并未产生任何 score 相对优势的翻转——恒等式对链路形状的普适性（A1 已证）没有转化为估计效率优势。
3. CLUSTER-WHITEN 的无监督聚类准确率 ≈ 0.50（三 seeds：0.506/0.502/0.500）：k-means 在 4096 维 log 图样上完全无法分辨相关长 2 px vs 3 px 的两态——而 oracle 后验判态准确率为 1.000（UT2）。这份对比是 score 方法"无标签状态识别"卖点在原理层的孤立正面证据，但它未能在任何重建指标门上兑现。
4. 诚实臂与信息上界（LIN/FGAIN/MLE-OR）按规范只报告不评门，数据齐备于 `results/phaseB_ext_metrics.csv`（本降级路径未跑核心表，`phaseB_core_metrics.csv` 不存在）。

图：`results/figures/B1_relation_headroom.png`、`B1_practical_headroom.png`（双门热图）、`A1_phase_diagram.png`、`GAM2_divergence.png`、`GAM2_convergence.png`、`DT30_calibration.png`。

## 6. 方法与裁决备注（协议登记）

### 6.1 实现等价性与模型定义（规范预授权范围内）
- **LOO 中心化**：$b_i-\bar b_{-i}=\frac{M}{M-1}(b_i-\bar b)$，按 LOO 公式实现；与普通中心化精确成比例、方向相同（规范头部备注，非"偏置修复"）。
- **核协方差 jitter**：GAUSS/CORR-LOGN 的高斯核矩阵数值秩亏，模型定义统一含对角 jitter 1e-8·σ²（采样、oracle score、WHITEN-OR 三处一致使用同一 Σ，自洽）。
- **DT30**：精确 non-paralyzable renewal（间隔 τ+Exp(λ) 累计计数，逐帧分块），非 Poisson(sat(u)) 冒充；UT3 双光子档标定 ≤0.05%。
- **LPIPS@16×16**：AlexNet 卷积塔在 16×16 塌缩，Phase A LPIPS 记 UNAVAILABLE（NaN）；Phase A 门均为 PSNR-only，不受影响。
- **FGAIN 增益**：g 按 (帧, 图) 在噪声流内独立抽取（与图样独立）。
- **CLUSTER-WHITEN**：无监督 k-means(2) 在 log 图样上（对数正态的自然特征化，仍无标签），簇内 LW 白化重建按帧数加权合并；聚类准确率随表报告。
- **L-Isotron 早停"似然"**：φ 未知噪声未知 → 用留出帧 MSE 作高斯似然代理（协议内文档化）。

### 6.2 门槛操作化登记（规范文字有歧义处的执行读法，均先于相应数据固化于代码）
- A1 通过门：逐探针比值 err/q 的 median/P90（gamma 用解析 q_k，其余用逐探针插值 SE）。
- A1 GAM1 三件套操作化：非零和未中心化 min 余弦 >0.99；零和 ‖v̂‖/SE median ≤1.5 且 P90 ≤2.0；中心化 err median ∈ [0.9,1.1]。
- A3 "至少一 combo 通过即 PASS"（筛查语义：检验头部空间是否存在于任一处）；两 combo 均报告。本轮两者皆大幅 FAIL，读法选择不影响结论。
- B1 bootstrap LB90>0 要求施加于 PSNR 增益（主 dB 门对应项）；SSIM/LPIPS 的 LB90 并报。
- B1 评门臂：非线性链路 = DT30 + 扩展数学臂；LIN/FGAIN 为诚实臂只报告；MLE-OR 只作信息上界。

### 6.3 A2 裁决轨迹（完整披露，供裁决/复核）
1. **跑前解析发现**：字面门"3-seed 平均重建 Pearson ≥0.995"对 SCORE-OR 在本工作点**纯方差即不可通过**（需图像中心化能量占比 ≥0.85，对 Σx=1 非负自然图像不可能；SCORE-OR 线性区方差 = k/(k−2)=2× WHITEN-OR）。0.995 的量级与尺度不变余弦读法吻合。
2. **第一版代码内登记**（跑数据前）：STOP 改按偏差证据 = min-图余弦 <0.995 或像素均值 χ²>2。
3. **对抗审计（跑 Phase A 前）确认两条 critical**：(i) χ² 均值门统计无效——3 seeds 下 z 精确 t(2)，E[z²]=∞，零假设下必触发（零假设模拟 P(触发)≈1.0）；(ii) 审计验证 agent 在冻结流上精确复现 A2：SCORE-OR 最差图（car）余弦 0.99479，恰为 2×方差预测值（无偏差信号）——min-余弦门同样被纯方差击穿。
4. **最终规则（Phase A 正式运行前固化）**：校准偏差检验 = 像素中位 z²（t₂ 零假设中位 0.667）> 2.0 任一图则 STOP。正式运行结果：全部 ≤0.834，无偏差。字面 Pearson / 余弦 / 原始均值 χ² 均按描述性并报于 `results/phaseA_gates.json`。
5. **诚实声明**：步骤 3 中审计复现已使 A2 数字先于"最终规则固化"而存在（执行者本人在固化规则时未读取除审计文本外的数字）。字面读法的不可通过性是解析事实（步骤 1），不依赖任何已见数据；若裁决方（用户）要求按字面 Pearson 执行，A2 应记 STOP——相应替代裁决块见 §8 附注。

### 6.4 其他
- **下溢帧**：全程 0 帧丢弃（各 (族,seed) 的 n_dropped 已持久化于 gates/meta JSON）。
- **隔离**：除只读拷贝两份 docs 外未读写 `E:\GAN_FCC_WORK`；未 import 旧项目代码。
- **gh CLI 本机不可用**：按规范未阻塞实验；仓库为本地 git（main 分支），请用户手动 `gh repo create GI_a2 --private --source D:\GI_another --push` 或添加 remote 后 push。
- **MLE-OR 优化中的 overflow RuntimeWarning**：L-BFGS-B 线搜索探测极端 x 时 λ 溢出，优化器自行回退；MLE-OR 仅诊断臂，结果有限且合理（表中 22.08/18.93）。

## 7. 运行时间（预算：UT+A ≤1.5 CPU h；B ≤6 CPU h）

| 阶段 | wall | 进程 CPU | 预算判定 |
|---|---|---|---|
| UT1–UT7（最终绿跑） | 371 s | ~0.2 h | ✓ |
| Phase A（A1+A2+A3） | 645 s | ~0.4 h（含 8 进程 MAVE 并行） | UT+A 合计 wall ≈ 0.28 h ✓ |
| Phase B 扩展表 + B 门 + 图 | 1462 s | 3.25 h（MKL 多线程，24 核） | wall 0.41 h / CPU 3.25 h，双读法均 ✓ |

预算行使用说明："CPU 小时" 按进程 CPU 时间与 wall 双报告；两种读法下均未超预算，未触发中断。全部本机 CPU 完成，未动用 Colab。

## 8. 裁决块

替代裁决附注（对应 §6.3 步骤 5）：若裁决方要求 A2 按字面 Pearson≥0.995 执行，则第 2 行改记
`STAGE0_SCOREGI_A2_LINEAR_HONESTY: STOP(literal-pearson; 方差一致、无偏差证据，min Pearson 0.783@dog)`，
后续 A3/B 行记 NOT_RUN，总裁决仍为 KILL 级终局（不改变科学结论的方向）。

```text
STAGE0_SCOREGI_UT: PASS
STAGE0_SCOREGI_A1_PHASE_DIAGRAM: PASS ; GAM1失效模式=确认, GAM2不稳定签名=异常(斜率签名-0.332确认不收敛; 尾比签名被共享帧设计结构性压缩, 附机制分析于phaseA_gates.json)
STAGE0_SCOREGI_A2_LINEAR_HONESTY: PASS(校准中位z²偏差检验, 全部≤0.834<2.0) ; LINEAR_EFFICIENCY_ONLY=否 ; 字面Pearson读法见§6.3裁决轨迹
STAGE0_SCOREGI_A3_SEPARATION: FLAGSHIP_DEAD (DT30×GAM4 −9.92dB 0/8; DT30×CORR-LOGN −12.76dB 0/8)
STAGE0_SCOREGI_B1_RELATION_HEADROOM: KILL (0/16 combos, 扩展表-only降级路径, 含DT30=否)
STAGE0_SCOREGI_B1_PRACTICAL_HEADROOM: KILL (0/16 combos, 含DT30=否)
STAGE0_SCOREGI_B2_REGIME: KILL (双门全灭, HONEST/THEORY_ONLY 判据无从满足)
STAGE0_SCOREGI_BOOTSTRAP: LB90>0 组合清单 = 空 (16/16 combo 的 PSNR增益 LB90 ∈ [−4.64, −0.39] dB, 全部<0)
STAGE0_SCOREGI_VERDICT: KILL
STAGE0_SCOREGI_RUNTIME: UT+PhaseA wall 0.28h; PhaseB wall 0.41h / process-CPU 3.25h (i9-13900HX 24C, 纯CPU)
```
