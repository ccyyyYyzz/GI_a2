# ROUND63 战役规范 D2.2（冻结候选版）：Dead-time-aware high-flux single-pixel imaging

**日期**: 2026-07-18（D2.2 第三次修订同日）
**状态**: D2.2 — 按 GPT 第六轮裁决（docs/ROUND63_GPT_ROUND6_RULING.md）：S1 试点证明旧协议被"正则化选择不可辨识 + 空间欠采样天花板"双重支配 → 主几何改 M/n=1、ν 网格九档、λ_TV 换 analytic_score_concentration 解析规则、删失分类 A–F 冻结；两步 S1（Pass A 校准 C₀ / Pass B 验证）后经一次不可变 commit 冻结为 F1。第四/五轮裁决的 AUDIT 拆分与描述性审计不变
**取代**: docs/ROUND63_HIGHFLUX_CAMPAIGN_SPEC_DRAFT.md（作废存档）
**产出**: Optics Express 投稿论文 *Dead-time-aware high-flux single-pixel imaging: operating beyond the conventional photon-counting regime with renewal-statistical reconstruction*

---

## 0. 预注册主假设与过门（逐字冻结）

> "We test whether renewal quasi-likelihood reconstruction permits operation
> in ρ̄ ∈ [0.3, 1] with a statistically significant reduction in total optical
> integration time relative to the conventional ρ̄ = 0.05 operating point."

- 安全参考点 **ρ₀ = 0.05**；**预注册主高通量点 ρ₁ = 0.6**；ρ = 1 为次级强负载点
- **两个工作点用同一 RQL 估计器**（隔离"换工作点收益"与"换模型收益"）
- 主质量目标（逐图）：Q_j\* = PSNR_rad(ρ=0.05, ν=2000, j) − 1 dB；固定 25/28/30 dB 为次级终点
- 主速度比：**S_j = T_opt,ρ=0.05(Q_j\*) / T_opt,ρ=0.6(Q_j\*)**，T_opt = M_physical·T；PSNR_rad–log T_opt 单调拟合插值；**删失按第六轮冻结分类 A–F**（T_min=M_physical·τ·ν_min, ν_min=5：B fast 触地 → S_gate=T_s/T_min 保守下界进全部门；C 双触地 → S_gate=1 进中位数但不计入 S>1 张数、标 unresolved below floor；D 反向 → S_gate=T_min/T_f<1 计非正；E fast 到 ν=2000 未达标 → 计失败不删；F 数据缺失 → ANALYSIS_FAILURE=0 并触发完整性审计；**bootstrap 每次抽样内重判删失类别**，不得预冻结逐图状态）
- **过门（三条同时）**：24 张确认自然图 **中位 S_j ≥ 3**；图像级分层 bootstrap **95% 下界 > 1**；**≥ 18/24 图 S_j > 1**
- 结果单位 = 图像（seed 图内配对平均；bootstrap 外层重采样图像、内层重采样 seed；不将 24×5 当 120 个独立样本）

## 1. 冻结边界（逐字条款）

> "S1 is exploratory and development-only. S1 images and seeds are disjoint
> from all S2 confirmatory images and seeds. S1 results shall not enter
> confirmatory confidence intervals or main quantitative tables."

> "Before the first S2 cell is executed, the following shall be frozen in one
> immutable commit: the primary hypothesis and pass rule; primary ρ₀, ρ₁;
> image and seed identifiers; all grids; reconstruction arms; TV-selection
> rule; quality targets; censoring rule; bootstrap procedure; figure subjects
> and crop coordinates; code, environment and input SHA256; all shard
> manifests; and the complete expected-cell table."

> "After S2 begins, changes are limited to demonstrable implementation
> defects established by outcome-blind unit tests. Every such change requires
> a version bump and rerunning all affected cells. No method, hyperparameter
> rule, grid point, image, seed, or failed/slow cell may be added, removed,
> or replaced based on observed reconstruction quality."

- **S1 白名单**：数值稳定；solver tolerance；memory/block size；manifest/shard 大小；无真值测试确定的路径覆盖；明确软件 bug
- **S1 黑名单**：主 claim；主高通量点；success gate；基线增删；最终图像/crop；按 PSNR 选 TV 规则；按初步结果删失配档
- **预算截断顺序（预冻结）**：先砍 jitter → 次级 PnP → 128² stress → **永不砍 S2-A、主失配锚点、S4**
- **SHA 硬门**：分片缺 SHA 即失败；代码 commit、环境 lockfile、图像、manifest、分析脚本全部进 provenance；合并前查 expected/重复/缺 cell
- 测试纪律：性能差异只作 smoke/report，禁作代码 PASS 门；真实 make_patterns→campaign→run_arm 集成测试为 S0 硬门

## 2. 物理模型（阻断项已修）

- 探测器：**τ = 50 ns 固定**（25/100 ns 只进"绝对量不缩放"敏感性）；non-paralyzable 主臂；帧起始 **active-start 主约定**：P(N≥m) = F_{Γ(m,λ)}(T − (m−1)τ)；delayed/continuous 为消融（continuous 修复后语义：carry = max(0, ready − T − g)，afterpulse 跨帧携带）
- 到达率 λ_i = Φ·⟨a_i, x⟩ + d（暗计数在死时间之前进入）；无量纲轴 ρ̄ = τ·E[λ]（**逐 pattern ρ 的 5/50/95 分位一并落盘**）、ν = T/τ
- **主观测 = 原始计数 N，σ_b = 0**；读出噪声只作次级消融（届时 exact 似然用离散卷积 p(b|λ)=Σ_m p_m·N(b;m,σ_b²)）
- **主估计器 RQL**（renewal quasi-likelihood）：min_{x≥0} (1/M)Σ[(T−N_iτ)λ_i − N_i log λ_i] + λ_TV·TV(x)——直接凸目标，无 IRLS；ceiling-count（N_iτ≥T）比例落盘为诊断；QMLE-IRLS 仅留 σ_b>0 扩展；QMLE-FULLGAUSS 为病灶消融臂
- exact 似然/Fisher：Poisson-logsf 表示 + logdiffexp + 解析 Ġ_m = e^{−z}z^m/(m−1)!（已实现 physics.exact_logpmf / exact_fisher_analytic，无有限差分、无截断）
- paralyzable 诚实定位：主文 = "paralyzable 生成器对 non-paralyzable RQL 的模型类别压力测试"；完整 2×2 需 matched-paralyzable RQL（若 S1 预算允许则补，否则如实声明）

## 3. 网格（GPT 授权最小充分集）

**S2-A 主战役（不可砍）**：64²，Bernoulli-50%，**M/n=1（M=4096，第六轮裁决：M/n≤0.5 的欠采样天花板使平滑自然图全光子轴动态范围仅 0.7 dB，不可辨识；M/n∈{0.25,0.5} 降入 S2-B 作为被报告的现象）**，ρ̄∈{0.05, 0.3, 0.6, 1, 2}，**ν∈{5, 10, 20, 50, 100, 200, 500, 1000, 2000}（九档；ν<20 论文称 short-window stress region，不得声称 RQL≈exact；固定 M 扫 dwell T=ντ）**，24 张确认自然图（STL-10 **test** split 索引 0..23，S1 不得接触），5 seeds。RQL 全网格；POISSON-LIN/SAT-POISSON/PRECORRECT 至少在 ρ̄∈{0.05, 0.6, 1}×全 ν；GI/DGI 选点+展示。M 扫描降级为次级"系统资源曲线"（附录）。
**S2-B 欠采样稳健**：M/n∈{0.25, 0.5, 1}×ρ̄∈{0.05, 0.6, 1}×ν∈{100, 500, 2000}，12 图、3 seeds。
**S2-C 图样/尺度锚点**：hadpair 与 gam4：ρ̄∈{0.05, 0.6, 1}×ν∈{100, 500, 2000}，M/n=0.5，8–12 图、3 seeds（互补对 M_physical=2M_signed 计费）；128²：Bernoulli，同锚点，仅 RQL/SAT/PRECORRECT/PnP-BM3D，**matrix-free 算子实现为前置硬门**；**PnP 冻结 = PnP-BM3D**。
**S3 失配（OAT + 三交互）**：固定 64²/bern50/M-n=0.5/ν=500；ρ̄∈{0.3, 0.6, 1}+参考 0.05；12 图、3 seeds。轴：τ 误差{−20,−10,0,+10,+20}%（+flat-field 标定部署版；联合 profile 只在 2 代表点）；暗{0,.05,.1,.25,.5}（已知/±10%/联合估计）；afterpulse{0,1,2,5,10}%；start modes{active,delayed,continuous}；guard{0,1,5}τ；jitter{0,5,10}% 附录。交互仅 ρ×τ-err、ρ×p_ap、continuous×afterpulse。**continuous 格子无独立 AUDIT（afterpulse 无界尾）：η\* 继承同 arm/image/seed/ρ/ν/A 的 active-start 格子，GOF_STATUS=GOF_NA_DEPENDENT，MODEL_FAIL_PREDICTIVE=NA（F1 冻结条款）。**
**S4 双拆**：(i) 标量 exact-Fisher map：ν=20..2000，ρ 自适应至 max(64, 2ρ\*(ν))（解析导数版）；(ii) 图像 exact-vs-RQL：8²/16²，ρ∈{.03,.1,.3,.6,1,2}×ν∈{20,100,500,2000}，3 seeds，**exact 与 RQL 同 TV 同选择器**。

## 4. 指标、图像与 λ_TV

- **主指标 = radiometric PSNR（PSNR_rad，不重标度）** + radiometric NRMSE + flux bias；flux-matched 改名 **shape-PSNR** 降为次级；SSIM/LPIPS 次级
- 图像：S1 开发集 = STL-10 **train** 0..15 + dev 结构靶（data/r63_images_dev/）；S2 确认集 = test 0..23 + 6 结构靶（结构靶只展示不进推断）；主图图名与 crop 坐标在冻结 commit 中登记
- **λ_TV 选择（F1 冻结 = 第六轮 `analytic_score_concentration`，取代第四轮 DEV 交叉拟合——后者被 S1 证明在桶式 SPI 多路复用劣势下结构性无功率、全网格坍缩到 η=1 完全平滑）**（实现 = code/round63/select_eta.py；裁决 docs/ROUND63_GPT_ROUND6_RULING.md）：**λ_TV,a = c_i·σ_{g,a}·√(2 ln n)**，σ_{g,a}=Φ√(κ_A·v_{s,a}/M)，v_{s,a}=各臂 score 在**精确 renewal 律**下的方差（physics.score_variance 枚举；ν=5,10 禁用 CLT 矩），κ_A=实测图案列能量；**c_i 两档**：Ĉ≤C₀→0.50，Ĉ>C₀→0.25；**Ĉ = clip[n(S_N²−V₀)₊/((Φμ′₀)²ω_A), 1, 64]**（DEV 原始计数 + 精确 PMF 矩，S_N²≤V₀→Ĉ=1）；**C₀ 为唯一开发集校准常数**（Pass A：endpoint-oracle regret，J=Q_{0.90}，tie 取大，冻结于 C0_FROZEN.json 后不可改）；全臂共用 c_i；逐 cell 落盘 C_hat/S_N2/V0/μ′₀/ω_A/c_used/σ_g/λ_TV；论文表述 "analytically noise-scaled TV with one development-calibrated concentration threshold"，禁称 training-free。每 cell 逻辑组（hadpair 对原子）冻结 hash 拆 **80% DEV / 20% AUDIT**（DEV 只算 Ĉ，AUDIT 只做描述性审计，均无超参搜索）；RNG 键全整数 (cell_key, seed, 63, 4, tag)
- **测量审计 = 纯描述性（第五轮终裁：无二元充分性门）**：结果盲探针证明 p≤0.025 门在该统计量上既无尺寸（正确模型 7/20 误报——refit 零假设条件于平滑 plugin 场景，对更粗糙真场景上尾反保守）也无功效（paralyzable/τ-err 0/10——探测器均值错配被场景尺度吸收，count-only 单工作点结构性不可辨识）。每 RQL cell 记录：AUDIT_STATUS∈{OK, UNDERPOWERED(<128 audit 组), NA_DEPENDENT(continuous)}、D_obs/D̄\*/sd/D_ratio、plugin_upper_rank q_D（**非标定 p 值**）、q_mean(+MEAN_RESIDUAL_WARN 描述性标记)、q_corr(+LOAD_CORR_WARN)、固定-λ̂ 下尾 q_low(+LEAKAGE_SUSPECT≤0.01)、residual-vs-load 图；**B_DIAG=39 全跑（禁早停——截断副本集的连续秩不可比）**；**任何审计统计量不得影响 η\*、重建、cell 取舍、战役门或确认推断**；探测器错配的后果由 S3 预注册 radiometric 指标直接量化。预注册措辞逐字条款见 ROUND63_GPT_ROUND5_RULING.md "Preregistered wording"。证据链：字面六步 10/10 误报 → 候选 A 8/10 → 第四轮门 7/20+0/10 → 删门（results/round63_gof_probe/）。truth-oracle λ 只作诊断，标 ORACLE—NOT DEPLOYABLE
- 公平性表述（论文逐字）：见 GPT 第三轮 REPLY §3.5 存档段落

## 5. 理论节（论文骨架）

- 主命题：**ρ\*(ν) ∼ 6^{1/3}ν^{1/3}**；formal second-order：**ρ\* = (6ν)^{1/3} − 2/3 + O(ν^{−1/3})**（"supported by exact numerical evaluation"）；峰值信息 J(ρ\*) = 1 − 0.8255·ν^{−1/3}
- 机制：missing-information identity **I_N = E[N] − ρ²·E[Var(R_ν|N)]**；边际平衡 1/ρ² = ρ/(6ν)；ridge 处 VarN 仍发散（禁用"方差 O(1)"表述）
- 分区：RQL deployment（逐 pattern ρ_95≤1）/ transition / exact-reference（"short-window & extreme-saturation reference mode"）/ information-decreasing（ρ≥ρ\*）；CLT 10% 信任边界 ρ_0.9 ≍ √(1.2ν)
- 表述纪律：**"We derive a finite-window count-information ridge…"**；无 "first"。**Grönberg 2018 逐式核对已完成（GPT 第四轮 Q2）：无 ridge/ν^{1/3}/identity——novelty "高置信可守"，但必须窄限定**（ideal nonparalyzable + active-start + scalar integrated count + exact finite-window FI + principal ridge + (ρ,ν) 渐近）；宽泛"最优通量"故事已拥挤（Wang 2011 / Bécares 2012 / Gupta / Rapp / Jorgensen 2026 / 2025 dToF——全部须防御性引用）；论文 novelty 段采用第四轮 digest 的安全版本文本；S5 残留硬项 = 对 2018 正文 PDF/KTH thesis 的最后人工关键词 sign-off
- Supplement 半页：*Why a full heteroscedastic Gaussian likelihood fails at high load*（ρ=1/2 符号翻转 + ceiling-row variance-collapse 非强制性 + radial objective 曲线 + 失败/修复图对）

## 6. 四联主图（修订版）

(a) exact count-FI 等高线 + principal ridge ρ\*(ν) + 10% discrepancy 边界 + 部署区着色 + CLT 虚线；(b) 预注册 1 自然图 + 1 文字/细线靶 × 5 列（安全全时/安全短时/高通量 naive/预校正/RQL），辐射显示范围统一，crop 冻结；(c) PSNR_rad–log T_opt 曲线族 + 目标线加速箭头 + bootstrap 区间；(d) empirical S_Q vs ρ̄ + exact-FI time-efficiency envelope + τ±10%/ap/背景带 + 失效区（CLT 不越 ridge 作上界）。

## 7. 分期与算力

S0'（本轮）：六硬门修复 + 规则三轮迭代（第 4/5/6 轮）→ **S1 两步**（Pass A = 唯一 C₀ 校准（c=.25/.50 端点 → regret → C0_FROZEN 单独 commit，此后端点/公式/目标不可改）；Pass B = 冻结 C₀ 下 M=4096 九档 ν 全臂验证：确认曲线不再平坦、存在可解析光子受限区、校准 Colab 分片时间）→ **F1 不可变冻结 commit**（含重生成 manifests/expected-cell 表/成本模型/SHA 账本/删失逻辑）→ S2/S3/S4（Colab pro2×3 + pro1×2 + 本地）→ S5 论文（OE opticajnl，Chen 体裁表述层，GPT 评审循环 ≥2 轮）。
预算：S1 ≤ 1 h 本地 wall；S2-A 为最大块（估 ~2–4 千 arm-fit，含 λ 选择 ~5× 开销——S1 校准后定分片）；超预算按 §1 截断顺序。
RNG：SEED0 体系，round63 流 (seed, 63, part, ...)；全输入输出 SHA256。
