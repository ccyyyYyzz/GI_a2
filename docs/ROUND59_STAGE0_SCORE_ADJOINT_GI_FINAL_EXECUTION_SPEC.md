# ROUND59 最终执行规范：分布伴随鬼影成像 Stage-0（交 Claude Code 执行）

**日期：** 2026-07-17
**写入者：** Fable5（整合 GPT Pro `APPROVE_WITH_EDITS` 全部 13 条修订，见 `ROUND58_GPTPRO_SPEC_AUDIT_RESPONSE.md`）
**状态：** 本文件**取代** `ROUND58_STAGE0_SCORE_ADJOINT_GI_EXECUTION_SPEC.md`（作废，仅存档）
**执行者：** Claude Code

**旗舰主张（本轮要裁决的科学命题）：** 复杂散斑照明分布的 score 是非线性鬼影成像的正确统计伴随；传统 Gaussian GI 是其线性-score 特例；该伴随关系随散斑模式数 $k$ 呈三区间定律——

$$k=1:\ \text{边界项不消失，恒等式失效};\qquad 1<k\le2:\ \text{恒等式成立但朴素平均方差无限};\qquad k>2:\ \text{恒等式成立且 }M^{-1/2}\text{ 正规收敛}.$$

**两点技术备注（Fable5，执行时照做即可）：**
- 修订 10 的 leave-one-out 中心化满足精确恒等式 $b_i-\bar b_{-i}=\frac{M}{M-1}(b_i-\bar b)$，即 LOO 与普通中心化只差全局常数 $\frac{M}{M-1}$，方向完全相同。按 LOO 公式实现并在报告注明该等价性，禁止把它写成"修复了偏置"。
- 强形式（条件–边缘 score 差）与原 A4/C2 已整体移出本轮（GPT 修订 7 的打乱门致命性成立），另案处理，本规范不含任何强形式内容。

---

## 0. 硬性守则

1. **工作目录与仓库**：工作目录为 `D:\GI_another`，作为新 GitHub 仓库 **GI_a2** 的本地根。第一步：`git init` 并用 gh CLI 创建私有仓库（`gh repo create GI_a2 --private --source D:\GI_another --push`）；若 gh 不可用，`git init` 后提示用户手动添加 remote，不阻塞实验。提交纪律：docs 就位后、UT 通过后、每个 Phase 完成后各 commit+push 一次，提交信息带阶段标签（如 `stage0: UT pass`）；`.gitignore` 排除图样矩阵缓存与任何 >50MB 文件；CSV/JSON/figures/REPORT.md 必须入库。
2. **隔离**：除允许从 `E:\GAN_FCC_WORK\paper_work\REVIEW_2026-07-11\innovation_physics_reboot\FABLE5_MAILBOX\` **只读拷贝**本规范与理论手册（`ROUND59_THEORY_PRIMER_FOR_CLAUDE_CODE.md`）到 `D:\GI_another\docs\` 外，禁止读写 `E:\GAN_FCC_WORK` 下任何内容；不 import 旧项目代码。
3. **必读前置**：执行前先通读 `docs/ROUND59_THEORY_PRIMER_FOR_CLAUDE_CODE.md`（完整数学物理理论）；`GAM1/GAM2` 臂"表现差"是理论预测的失效模式，属于验收对象而非 bug，禁止"修复"。
4. **纯 CPU、无学习**：本轮禁止任何神经网络训练。Phase C 仅预告（§10），另行授权。
5. **预登记不可改**：门槛以本文件为准；结果无论好坏原样报告；禁止调参重跑后只报最好一组。
6. **失败即停**：STOP 级门失败即中止后续 Phase，写报告与裁决块，commit 后停止。
7. **预算**：UT+Phase A ≤ 1.5 CPU 小时；Phase B ≤ 6 CPU 小时。超预算中断并报告瓶颈。
8. **冻结随机性**：主种子 `SEED0=20260717`。同一 `(照明族, seed)` 的图样矩阵 $\{a_i\}$ 生成一次、冻结复用于所有链路与光子档（协方差与白化算子同样复用）；不得每链路重抽。

## 1. 环境与全局常量

- Python ≥3.10；`numpy scipy scikit-image scikit-learn`；LPIPS 用 `torch`(CPU)+`lpips`(AlexNet，记录版本)。lpips 不可用 ⇒ 照跑 PSNR/SSIM，LPIPS 标 `UNAVAILABLE`，相关门记"待补"，不得视为通过。
- 帧数 $M=20000$；witness 专用 $M_w=2\times10^5$（只抽图样、不算 bucket）。噪声种子：Phase A `{0,1,2}`；Phase B `{0,1,2,3,4}`。
- 图像域 `[0,1]`；照明均值定标使 $U_0=\mathbb E[\langle a,x\rangle]=1$（对每幅图用 $x$ 的和归一实现，UT4 复核）。
- **测试图像（8 张 64×64 灰度）**：优先 STL-10 test split，{cat, deer, dog, horse} 与 {airplane, car, ship, truck} 各取该类第一张，灰度化、双三次到 64×64，存 PNG+SHA256。下载不可行 ⇒ 退回 `skimage.data` 内置图并逐张声明替代。Phase A 用同批图下采样 16×16。真值统一预归一化 $\sum_j x_j=1$ 后再缩放回 `[0,1]` 展示（指标协议见 §5）。

## 2. 照明族（oracle score 全解析；采样细节固定）

| 代号 | 分布 | oracle score $s(a)$ | 角色 |
|---|---|---|---|
| `GAUSS` | $a=\mu+\Sigma^{1/2}z$，$\mu=\mathbf 1$，$\Sigma$=高斯核相关（相关长 2 px，σ=0.25） | $-\Sigma^{-1}(a-\mu)$ | **数学恒等式控制**（非物理臂：单帧含负像素概率 n=4096 时约 12%，报告中声明） |
| `GAM1` | iid $\mathrm{Gamma}(1,\text{rate}=1)$ | $-\mathbf 1$ | 边界失效展示臂（k=1） |
| `GAM2` | iid $\mathrm{Gamma}(2,\text{rate}=2)$ | $(k-1)/a_j-k$ | **预测不稳定臂**（恒等式成立、$\mathbb E[s^2]=\infty$） |
| `GAM3/GAM4/GAM8` | iid $\mathrm{Gamma}(k,\text{rate}=k)$，$k\in\{3,4,8\}$ | 同上，禁止钳位（改用 float64 与安全 log；若出现下溢帧则丢弃并计数报告，不得截断 $a$） | 有限方差通过臂（$\mathbb E[s_j^2]=k^2/(k-2)$，UT5 用） |
| `CORR-LOGN` | $z\sim\mathcal N(m,\Sigma_z)$，$a=e^z$；$\Sigma_z$=高斯核相关（相关长 2 px，$\sigma_{\ln}^2=0.25$），$m=-\sigma_{\ln}^2/2$（单位均值） | $-\operatorname{diag}(1/a)\big[\Sigma_z^{-1}(\log a-m)+\mathbf 1\big]$ | **空间相关×正值×非高斯主物理臂**（湍流型乘性涨落） |
| `MIX-LOGN` | 两状态各 0.5：CORR-LOGN 相关长 2 px vs 3 px（同边缘），逐帧抽状态、标签不给任何估计器 | log-domain 后验加权 $\sum_c w_c(a)s_c(a)$，$w_c\propto\pi_c p_c(a)$（Cholesky+三角解算） | 动态状态混合臂；配套经典公平控制 `CLUSTER-WHITEN`（§4.11），并报告无监督聚类准确率 |

## 3. bucket 链路与噪声

**噪声基线**：$b=\mathrm{Poisson}(s\cdot\phi(u))/s+\mathcal N(0,\sigma_r^2)$，$u=\langle a,x\rangle$，光子档 $s\in\{10^3,10^4\}$，$\sigma_r=2/s$。

| 代号 | 定义 | 角色 |
|---|---|---|
| `LIN` | $\phi(u)=u$ | 线性诚实臂 |
| `DT30` | **物理臂**：非延长死时间 renewal 计数精确模拟——到达率 $\lambda=s\,u$、窗口 $T=1$、死时间 $\tau=\dfrac{3}{7s}$（$\tau\lambda_{\rm mean}=3/7\Rightarrow$ 均值压缩 30%）。实现：接收间隔 iid $=\tau+\mathrm{Exp}(\lambda)$，累计求和数落在 $T$ 内的个数，$b=N/s+\mathcal N(0,\sigma_r^2)$。$\mathbb E[b\mid a]\approx u/(1+\tfrac37u)$ 单调单指标 ✓，噪声为 renewal 而非 Poisson ✓ | **旗舰物理非线性臂**（B2 的 HONEST 判据只认它） |
| `FGAIN` | $b=\mathrm{Poisson}(s\,g_i u)/s+\mathcal N(0,\sigma_r^2)$，$g_i\sim\mathrm{LogNormal}$，CV=0.2，与图样独立 | 动态标度诚实臂（预期 score 无偏置优势，如实报告） |
| `SAT30`/`SAT50` | $\phi(u)=u/(1+\alpha u)$，$\alpha$ 使 $\phi(U_0)/U_0=0.7/0.5$ | 数学应力臂（扩展表） |
| `GAMMA07` | $u^{0.7}$ | 数学应力臂（扩展表） |
| `LOG` | $\log(1+\alpha u)$，$\alpha$ 使 $\phi(U_0)/U_0=0.7$ | 数学应力臂（扩展表） |

## 4. 估计器（同数据、同帧、同光子、同后处理）

1. `GI`：$\frac1M\sum(b_i-\bar b)(a_i-\bar a)$
2. `DGI`：差分 GI 标准式（$b_i-\frac{\bar b}{\bar r}r_i$，$r_i=\sum_j a_{ij}$）
3. `CORR`：correspondence GI（bucket 上/下四分位图样均值差）
4. `SIR-10` / `SIR-20`：多片 sliced inverse regression——用 `WHITEN-LW` 的 $\hat\Sigma_{\rm LW}^{-1/2}$ 标准化图样，按 $b$ 分 10/20 片，片均值加权外积的主特征向量，反变换回原坐标
5. `WHITEN-LW`：$\hat\Sigma_{\rm LW}^{-1}\cdot\frac1M\sum(b_i-\bar b)(a_i-\bar a)$（sklearn LedoitWolf；每 `(照明,seed)` 一次复用）
6. `WHITEN-OR`：同上用真 $\Sigma$（各族解析）
7. `SCORE-OR`：$-\frac1M\sum_i(b_i-\bar b_{-i})\,s(a_i)$（LOO 中心化，见头部备注）；**GAM1 专属附加输出**未中心化版 $-\frac1M\sum b_i s(a_i)$（展示常数塌缩 vs 零塌缩双图证）
8. `L-ISOTRON`：Isotron/L-Isotron 冻结协议——迭代 $x^{t+1}=x^t+\eta\cdot\frac1M\sum_i\big(b_i-\hat\phi_t(\langle a_i,x^t\rangle)\big)a_i$，$\hat\phi_t$ 由 PAV isotonic 回归（Lipschitz 截断可选），$\eta=1/\mathrm{tr}(\hat\Sigma)$，最多 200 轮或验证集（10% 留出帧）似然容差 $10^{-6}$ 早停；3 个初始化（`WHITEN-LW`、`SIR-10`、随机）按留出似然选优，**不得用真值选择**
9. `MAVE-16`：MAVE（或 Poisson 局部似然 MADE）——仅 Phase A 16×16 运行（rMAVE 标准实现，OPG 初始化）
10. `MLE-OR`：已知真 $\phi$ 的 Poisson 似然 MLE（L-BFGS，非负投影，`WHITEN-LW` 初始化）——**信息上界诊断，不进任何头部空间门的"传统控制"集合**
11. `CLUSTER-WHITEN`（仅 MIX-LOGN）：k-means(2) 帧聚类 → 每簇 LW 白化重建 → 按簇合并；报告聚类准确率
12. `RANKG`（可选）：逐像素 rank-Gaussianize 后跑 `GI`（陈组 Gaussian-constraint 机制的内部对照）

## 5. 指标协议（消除真值泄漏）

**主协议（部署合法）**：$\hat x_+=\max(\hat x,0)$；通量匹配归一化 $\hat x_\star=\hat x_+\cdot\frac{\sum_jx_j}{\sum_j\hat x_{+,j}}$（真值已归一 $\sum x=1$；等价于双方单位和）。PSNR（data_range=真值最大值）、SSIM、LPIPS 均在 $\hat x_\star$ 上计算。
**附表诊断（`DIRECTION_ONLY_DIAGNOSTIC`）**：真值 LS 标量拟合版指标 + 尺度不变角误差 $\arccos\frac{\langle\hat x_+,x\rangle}{\|\hat x_+\|\|x\|}$ + Pearson 相关。两套并报，主门只看主协议。

## 6. 单元测试（先于一切 Phase，全过才继续）

- UT1：`GAUSS`×`LIN` 无噪、$M=10^5$：`SCORE-OR` 与 `WHITEN-OR` 逐像素相对差 ≤1e-2（高斯下恒等）。
- UT2：各族样本均值/方差 vs 解析值相对差 ≤1%；`MIX-LOGN` 两状态样本上 oracle 后验平均判别准确率报告（不设门，供公平性讨论）。
- UT3：Poisson 均值标定 $\mathbb E[b]$ vs $\phi(U_0)$ 相对差 ≤1%；**DT30 专项**：模拟均值计数 vs $\lambda T/(1+\lambda\tau)$ 相对差 ≤1%。
- UT4：$U_0=1$ 与各链路 $\phi(U_0)/U_0$ 标定复核（差 ≤0.01）。
- UT5：中心化 witness 的解析方差 $q_k^2=\frac{k(n+1)}{(k-2)M_w}$ 与 Monte Carlo 方差匹配（k∈{3,4,8}，比值 ∈[0.8,1.25]）。
- UT6：`GAM2` score 二阶矩随下截断阈值 $\epsilon\to0$ 发散（对 $\epsilon\in\{10^{-2},10^{-3},10^{-4},10^{-5}\}$ 展示单调增长无饱和）。
- UT7：主协议与 LS 诊断两套指标管线各自输出且互不污染（同一重建两套数字并存）。

## 7. Phase A（16×16，n=256；3 seeds）

**A1 —— Stein 有效性相图（STOP 级）**
探针族（64 个，冻结）：32 个强制零均值随机 $v$、16 个非零和随机 $v$、16 个结构向量（常数、横/竖条纹、棋盘、低频 Fourier 模）。witness 一律用**中心化探针** $\langle a_i-\mu,v\rangle$，$M_w=2\times10^5$：

$$\widehat v=-\frac1{M_w}\sum_i\langle a_i-\mu,v\rangle\,s(a_i),\qquad \mathrm{err}(v)=\|\widehat v-v\|_2/\|v\|_2.$$

门（分布特异理论/插值 SE）：对 gamma 族用 $q_k=\sqrt{\frac{k(n+1)}{(k-2)M_w}}$；对 `GAUSS`/`CORR-LOGN` 用插值 SE $\hat q=\sqrt{\widehat{\mathrm{Var}}[\langle a-\mu,v\rangle s(a)]/M_w}$（样本插值，逐 $v$）。判定：

- `GAM3/GAM4/GAM8/CORR-LOGN/GAUSS`：$\operatorname{median}(\mathrm{err})\le1.5q$ 且 $P_{90}(\mathrm{err})\le2q$ ⇒ PASS；任一族 FAIL ⇒ **STOP**（实现或理论有错）。
- `GAM1`：非零和探针的未中心化输出与 $(\sum_jv_j)\mathbf 1$ 余弦 >0.99；零和探针输出范数 ≤ 插值 SE 地板；中心化输出 err≈1。不满足此**失效模式** ⇒ **STOP**（边界条件理论有错，上报）。
- `GAM2`：展示不稳定签名——err 的 $P_{90}/\mathrm{median}$ 比值 ≥3× `GAM4` 同值，且 err 随 $M_w\in\{10^4,10^5,2\times10^5\}$ 不按 $M^{-1/2}$ 收敛（log-log 斜率 > −0.35）。签名缺席 ⇒ 记 `GAM2_ANOMALY` 并附分析（不 STOP）。

**A2 —— 线性区诚实检查（STOP 仅限偏差）**
`LIN`×`GAM4`×$s=10^4$：`WHITEN-OR` 与 `SCORE-OR` 各 3 seed 的**平均重建**与真值 Pearson 相关 ≥0.995（无偏数值确认）。两者 PSNR/方差/MSE 差**只报告不设胜负门**；若 score 在线性区占优，标 `LINEAR_EFFICIENCY_ONLY`，不计旗舰增益。仅当出现系统性方向偏差（平均重建相关 <0.995）⇒ **STOP** 排查。

**A3 —— 非线性头部空间筛查（KILL 级）**
`DT30`×{`GAM4`,`CORR-LOGN`}×$s=10^4$：

$$\mathrm{SCORE\text{-}OR}\ \ \text{vs}\ \ \max(\mathrm{WHITEN\text{-}OR},\ \mathrm{WHITEN\text{-}LW},\ \mathrm{SIR\text{-}10/20},\ \mathrm{L\text{-}ISOTRON},\ \mathrm{MAVE\text{-}16})$$

门：平均 PSNR ≥ +0.5 dB 且 ≥6/8 图为正（主协议指标）。不过 ⇒ `FLAGSHIP_DEAD`（仍跑完 Phase B 扩展表供降级论文，裁决块如实记）。

**（原 A4 已删除，无强形式内容。）**

## 8. Phase B（64×64，n=4096；5 seeds；A1/A2 过后执行）

**核心物理表（决策表）**：照明 {GAM4, GAM8, CORR-LOGN, MIX-LOGN} × 链路 {LIN, DT30, FGAIN} × 光子 {1e3, 1e4} × 5 seeds × 8 图 × 估计器 §4（`MAVE-16` 除外；`CLUSTER-WHITEN` 仅 MIX-LOGN）。
**扩展表（regime 映射）**：同照明 × {SAT30, SAT50, GAMMA07, LOG} × $s=10^4$ × 3 seeds。

**B1 —— 双头部空间门（KILL 级；仅非线性链路参与评门）**

- `RELATION_HEADROOM`：SCORE-OR vs $\max(\mathrm{WHITEN\text{-}OR},\ \mathrm{SIR\text{-}10/20},\ \mathrm{L\text{-}ISOTRON})$；
- `PRACTICAL_HEADROOM`：SCORE-OR vs $\max(\mathrm{WHITEN\text{-}LW},\ \mathrm{DGI},\ \mathrm{CORR},\ \mathrm{SIR\text{-}10/20},\ \mathrm{L\text{-}ISOTRON})$。

每门要求（主协议指标）：PSNR ≥+0.50 dB、SSIM ≥+0.010、LPIPS benefit ≥+0.010、三项同时改善 ≥6/8 图；在 ≥2 个非高斯照明 × ≥2 个非线性链路组合上通过，其中**必须含 DT30**。任一门全网格无通过组合 ⇒ `STAGE0_KILL`。`MLE-OR` 同表报告作信息上界，不入门。

**B2 —— regime 判据**：`HONEST_NONLINEARITY` 当且仅当 DT30（中等、物理化非线性）上双门通过；仅 SAT50/GAMMA07/LOG 通过 ⇒ `THEORY_ONLY`。

**统计协议**：先对每图取 5-seed 均值，再对 8 图做 paired hierarchical bootstrap（10^4 次重采样）；除均值门外，每个"通过组合"的增益 90% bootstrap lower bound 必须 >0；逐 seed 方向全部报告，禁止只报 pooled mean。预算允许时可扩 16 图（优先于加 seed）。

## 9. 输出契约

```
D:\GI_another\               # git 仓库 GI_a2 根
  docs/                      # 本规范 + ROUND59_THEORY_PRIMER_FOR_CLAUDE_CODE.md（只读拷贝）
  code/                      # 入口 run_all.py（UT → A1 → A2 → A3 → B）
  data/                      # 8 图 PNG + sha256.txt
  results/
    unit_tests.json
    phaseA_witness.csv       # (family, probe_type, v_idx, err, q_theory)
    phaseA_gates.json
    phaseA_gam1_exhibits/    # 未中心化常数塌缩图 + 中心化零塌缩图
    phaseB_core_metrics.csv  # (image, illum, link, photons, seed, method, PSNR, SSIM, LPIPS, angerr, pearson)
    phaseB_ext_metrics.csv
    phaseB_gates.json        # 含 bootstrap LB90
    figures/                 # A1 相图、B1 双门热图、GAM2 发散曲线、DT30 标定曲线
  REPORT.md
```

`REPORT.md` 结尾必须为：

```text
STAGE0_SCOREGI_UT: PASS | FAIL(...)
STAGE0_SCOREGI_A1_PHASE_DIAGRAM: PASS | STOP(...) ; GAM1失效模式=确认|异常, GAM2不稳定签名=确认|异常
STAGE0_SCOREGI_A2_LINEAR_HONESTY: PASS | STOP(...) ; LINEAR_EFFICIENCY_ONLY=是|否
STAGE0_SCOREGI_A3_SEPARATION: PASS | FLAGSHIP_DEAD
STAGE0_SCOREGI_B1_RELATION_HEADROOM: PASS(m/n combos, 含DT30=是|否) | KILL
STAGE0_SCOREGI_B1_PRACTICAL_HEADROOM: PASS(m/n combos, 含DT30=是|否) | KILL
STAGE0_SCOREGI_B2_REGIME: HONEST_NONLINEARITY | THEORY_ONLY
STAGE0_SCOREGI_BOOTSTRAP: LB90>0 组合清单
STAGE0_SCOREGI_VERDICT: PROCEED_TO_PHASE_C | THEORY_PAPER_ONLY | KILL
STAGE0_SCOREGI_RUNTIME: <CPU 时长>
```

## 10. Phase C 预告（本轮不执行；协议冻结，另行授权）

仅弱形式 learned marginal score：网络只见目标无关参考散斑（不见 bucket、不见物体、不见自然图像）；DSM σ 阶梯 {0.05,0.1,0.2} 并报告 σ-偏置；held-out Stein witness（A1 同门槛）；与经典 score 估计器（低维核估计）对照；错误散射状态 score 删除 ≥70% 增益；未见散射强度保留 ≥50% 增益；**最终必须在真实波光学散斑**（随机相位屏、角谱/Fresnel、有限 NA、多层动态散射——此处无解析 oracle，witness 为唯一合格证）上恢复 ≥80% 解析族头部空间。学习恢复 ≥80% oracle 增益、三指标 6/8 门沿用。

---

**给 Claude Code 的一句话总纲：** 顺序 UT → A1 → A2 → A3 → B，STOP 级失败立即停止；所有门槛与协议数值不可修改；`GAM1/GAM2` 是预测失效/不稳定臂——它们"表现差"是理论确认而非 bug，必须按 §7 的失效模式验收；结果无论好坏原样写入 `REPORT.md` 并以裁决块收尾。
