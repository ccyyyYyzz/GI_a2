# ROUND63 GPT Pro 第三轮回复 — 完整可执行摘要（Opus 蒸馏，主循环核验）

原文：docs/ROUND63_GPT_ROUND3_REPLY.md（2461 行，完整无缺）。本摘要为操作蓝图。

## 总裁决

| 项目 | 裁决 |
|---|---|
| ROUND63 规范草稿 | **HOLD：不能冻结，必须出 D2 版** |
| 高通量 A 主线 | 继续，故事比第二轮更强 |
| 立方根脊线 | **成立；领先常数 6^(1/3)=1.81712（非经验 1.7）** |
| "方差降到 O(1)"解释 | **错误，必须删**（正确=边际平衡 ρ³=6ν） |
| 命名 | **renewal quasi-likelihood (RQL)**；Gaussian 仅为矩来源 |
| Wedderburn IRLS | 方向正确，但 σ_b=0 主模型有**直接凸目标，无须 IRLS** |
| own-NLL TV 选择 | **废弃**；换 cross-fitted common renewal discrepancy 六步规则 |
| 理论独立短文 | 暂不拆；作 A 的强理论小节 |

## 1. 规范修订清单（D2 冻结前逐项）

### 命名与预注册（阻断级）
- 主臂改名 **RQL** (renewal quasi-likelihood reconstruction, solved by quasi-score)
- §1 改为待检验假设（逐字）: "We test whether renewal quasi-likelihood reconstruction permits operation in rho_bar in [0.3,1] with a statistically significant reduction in total optical integration time relative to the conventional rho_bar = 0.05 operating point."
- 预注册: ρ₀=0.05 参考点；**主高通量点 ρ₁=0.6**；ρ=1 次级；两边同一 RQL 估计器
- 主目标 Q_j* = PSNR_rad(ρ=.05, ν=2000, j) − 1 dB；主速度比 S_j = T_{ρ=.6}(Q_j*) / T_{ρ=.05}(Q_j*)
- 过门（三条同时）：24 自然图中位 S_j ≥ 3；分层 bootstrap 95% 下界 > 1；≥18/24 图 S_j > 1

### 物理与似然（六个阻断级）
A. 主模型= 原始计数 N、σ_b=0；读出噪声只做次级消融（届时 exact 用离散卷积）
B. active-start 公式冻结：P(N≥m) = F_{Γ(m,λ)}(T − (m−1)τ)
C. exact-PMF 出版版：logdiffexp、CDF/SF 自适应分支、解析 Ġ_m = e^{−z_m} z_m^m/(m−1)!（z_m=ρ(ν−m+1)）、log-domain FI、取消 cutoff 与有限差分
D. continuous guard 符号 BUG：正确 carry_{i+1} = max(0, ready_i − T − g)；afterpulse 跨帧携带
E. paralyzable 2×2 目前只是名义；诚实定位为"模型类别压力测试"或补 matched paralyzable RQL；Liu-2021 paralyzable 格不能砍
F. 核心相图固定 τ=50ns；25/100ns 只在绝对量不缩放的敏感性里跑

### 网格（OE 最小充分战役）
- **S2-A（不能砍）**: 64², bern50, M/n=0.5, ρ∈{0.05,0.3,0.6,1,2}, ν≥5 个对数点(20..2000), 24 确认自然图, 5 seeds；RQL 全曲线；LIN/SAT/PRECORRECT 至少在 ρ∈{0.05,0.6,1}；**time-to-quality 固定 M 扫 dwell T=ντ**（不是扫 M！M 扫描降级为次级资源曲线）；主时间 T_opt = M_physical·T
- **S2-B**: M/n×ρ×ν 缩小网格（12 图 3 seeds）
- **S2-C**: Hadamard/GAM4 与 128² 只跑锚点；**PnP 冻结为 PnP-BM3D**；128² 必须 matrix-free（dense A 2–4 GiB 会爆 Colab）
- **S3**: OAT 失配（τ±10/20%、暗 0..0.5、ap 0..10%、start modes、guard、jitter 附录）+ 三个交互（ρ×τ-err、ρ×p_ap、continuous×ap）
- **S4 拆两半**: 标量 exact-Fisher map（ρ 扫到 ≥max(64, 2ρ*)）+ 图像 exact-vs-RQL（8²/16²，同 TV 同选择器——当前 EXACT 无 TV 不能叫 reference）

### 统计与图像纪律
- **radiometric PSNR（不重标度）为主指标**；flux-matched 改名 shape-PSNR 降为次级
- 结果单位=图像（seed 在图内配对平均；bootstrap 图像外层）
- **S1 与 S2 图像/种子完全不相交**（S1 用 STL train 或独立开发集）；主图图名+crop 坐标 S2 前冻结
- 冻结边界三段逐字条款（S1 exploratory-only / one immutable commit 冻结清单 / S2 后仅 outcome-blind 缺陷修复+版本 bump）+ 白名单/黑名单 + 预算截断顺序预冻结 + SHA 硬门（缺=失败）
- 测试纪律：性能差异不作代码 PASS 门；补真实 make_patterns→campaign→run_arm 集成测试（现 Hadamard 单测用错 meta key 绕开了真 bug）

### 四联主图修订
(a) exact count-FI 等高线 + principal ridge + 10% discrepancy 边界 + 部署区 + CLT 虚线；(b) 图与 crop 预注册、辐射显示范围统一；(c) 横轴 T_opt、纵轴 radiometric PSNR；(d) empirical speedup vs exact-FI envelope（CLT 不得越 ridge 作上界）

## 2. 立方根定律最终形态

- 主命题：**ρ*(ν) ∼ 6^(1/3) ν^(1/3)**；二阶：**ρ* = (6ν)^(1/3) − 2/3 + O(ν^(−1/3))**
- 峰值信息：J(ρ*) = 1 − (1/c + c²/12)ν^(−1/3)，c=6^(1/3)（系数 ≈0.8255）
- **机制**：missing-information identity **I_N = E[N] − ρ²·E[Var(R_ν|N)]**（R=窗末残余死时相位）；1/ρ = 距完整事件信息饱和缺口；ρ²/(12ν) = 丢终端相位代价；边际平衡 1/ρ² = ρ/(6ν) → **ρ³=6ν**
- 大ρ区间展开：J = 1 − 1/ρ − ρ²/(12ν) + 1/ρ² − ρ/(6ν) + ρ⁴/(144ν²) + o(ν^{−2/3})
- 终端相位矩：P(R=0)=1/(1+ρ)；f_R(r)=ρ/(1+ρ) on (0,1)；VarR=ρ(ρ+4)/[12(1+ρ)²]；Cov(N,R)→1/12
- 数值吻合：(6ν)^(1/3)−2/3 = 4.27/6.03/7.77/9.96/13.76/17.50/22.23 vs repo 4.53/6.16/7.87/9.99/13.77/17.45/22.16 ✓
- ridge 处 VarN = ν^(1/3)/6^(2/3) 仍发散——"计数钉死"在更晚的 ρ≍√ν
- 先例：Müller 1973（count law）、Yu–Fessler 2000（矩+反演）、Alvarez 2014（CLT 代理）、**Grönberg 2018（最危险，须逐式核对全文）**、Rapp–Goyal 2019/21（时序高通量加速——"高通量更快"非我们首创）、Jorgensen 2026（activation statistics）；**净判：ρ*∼(6ν)^(1/3) 的 active-start count-only ridge 未检索到先例**
- 表述："We derive a finite-window count-information ridge…"；6^(1/3) 主命题 + −2/3 formal second-order ("supported by exact numerical evaluation")

## 3. RQL 最终形式（关键简化）

σ_b=0 主模型的直接凸目标（**无须 IRLS**）：
  min_{x≥0} (1/M) Σ_i [(T − N_i τ)λ_i(x) − N_i log λ_i(x)] + α·TV(x)，λ_i = Φ a_i^T x + d
- quasi-score U_λ = N/λ − (T − Nτ)；单帧根 λ̂ = N/(T−Nτ) **恰为 PRECORRECT 公式** → novelty 必须落在工作点地图/time-to-quality/ridge/失配稳健/联合逆问题稳定性，不在估计器公式
- 饱和边界 Nτ=T：目标无下界=真实辨识性边界；落盘 ceiling-count 比例；exact 负责该区
- IRLS 仅在 σ_b>0 或 afterpulse 矩模型时保留
- log-det 病灶正确刻画 = **variance-collapse non-coercivity**（顶格行 NLL→−∞ 严格退化）；ρ=1/2 符号翻转正确但"必然造成偏差"说法不可用（真参数处期望 score 仍为零）；写半页 Supplement "Why a full heteroscedastic Gaussian likelihood fails at high load"

## 4. λ_TV 六步规则（cross-fitted common renewal discrepancy）

1. K=5 冻结 hash folds（互补对同 fold；continuous 按 block）
2. 无量纲路径 η = λ_TV/λ_max,arm ∈ {1e-4…1}（λ_max 由冻结 bisection 求）
3. 各臂用自己的 fidelity 拟合（方法定义不变）
4. **共同外部物理评审器**：held-out 上 r_i = [N_i − μ_NP(λ̂_i)]/√V_NP(λ̂_i)，D=Σr²；阈值由 exact renewal parametric bootstrap 标定（2.5–97.5% 接受区 + 均值/负载相关性检验），μ_NP/V_NP 来自独立探测器标定
5. one-SE + 最大相容正则化：E = {η: D_dev ≤ min+SE 且 GOF 通过}，η* = max E；E=∅ → MODEL_FAIL（禁看 PSNR 补救）
6. 冻结 η* 全数据重拟合
- 理论背书：Morozov discrepancy principle；SURE/GSURE 只作附录敏感性
- 公平性表述（逐字英文段落已在 REPLY 存档）
- truth-oracle λ 只作诊断，标 ORACLE—NOT DEPLOYABLE

## 5. 分区表述与硬门

- J_exact/J_CLT ≃ 1 − ρ²/(12ν) → CLT 10% 信任边界 **ρ_0.9 ≍ √(1.2ν)**（与 ridge ν^(1/3) 是两个不同标度）
- 四区：RQL deployment（ρ_95≤1）/ transition / exact-reference（"short-window & extreme-saturation reference mode"，不说 low-count mode）/ information-decreasing（ρ≥ρ*，换 exact 也救不回）
- 分区用逐 pattern ρ_95 不用均值；FI ratio 是标量诊断，S4 图像级验证不可省
- 硬门汇总：SHA 缺=失败；性能不作 PASS 门；guard 符号修复+afterpulse 跨帧；exact-PMF log 域；128² matrix-free；真实集成测试
