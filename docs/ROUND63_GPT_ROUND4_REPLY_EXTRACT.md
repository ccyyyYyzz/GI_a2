# GPT Pro round-4 reply — re-assembled extract (2026-07-18)

Provenance: chat 候选方案评审与建议 (chatgpt.com/c/6a5a66fb-11e0-83eb-b618-3a433a83935d),
4th assistant message, 15,192 chars innerText, ~81-min Pro run. Extracted in
slices via DOM innerText; math is linearized by innerText (sub/superscripts
flattened); tokens ≥22 chars were masked as [TOK] to pass the relay's data-loss
filter; chars ≈10400–10620 are unrecoverable through the relay
([DLP-HOLE] — from context: the list of ROUND63-side contributions vs Grönberg,
i.e. dimensionless (ρ,ν) axes / exact finite-window PMF FI / ridge asymptotics /
missing-information identity). The operative distillation is
ROUND63_GPT_ROUND4_DIGEST.md; on any discrepancy the ChatGPT thread governs.

---

## Q1 — BLOCKING 裁决：改用"外层审计集 + coherent refit bootstrap"

裁决：采用 A 的统计思想，但不采用当前"每个 bootstrap replicate 重跑完整 K-fold
cross-fit"的实现。最终冻结版应改成：DEV 上选择 η，独立 AUDIT 上做模型充分性检验；
bootstrap 每次只重拟合一个 coherent RQL 图像，而不是把 cross-fitted λ̂_i 向量当成
可生成真值。

你们对原缺陷的诊断完全正确：E[D] ≈ M + Σ_i {μ(λ_i)−μ(λ̂_i)}²/V(λ̂_i)。
固定 λ̂ 的 bootstrap 只重现第一项 counting noise，检验的是"λ̂_i = λ_i"而不是
"renewal detector model + reconstruction procedure 是否能做出合格预测"。

B 和 C 都不适合单独成为冻结方案：
- 只检验 r̄ 和 corr(r,ρ) 会漏掉零均值、负载无关、但方差很大的预测失败；TV bias
  也不保证在 r 中均值为零。
- one-SE 规则只是正则化选择，不是模型充分性检验。
- 额外添加 M + c√(2M) + Î_plugin 没有自然的有限样本校准，而且很容易重复计算
  estimation variance。

另一个必须诚实写清的限制是：在 M<n 且每个 pattern 只有一个 bucket count 时，
纯数据检验无法严格区分"detector likelihood 错配"和"场景在 null/near-null space
中没有被充分恢复"。

所以建议将标志正式改名为 MODEL_FAIL_PREDICTIVE。其语义是：calibrated renewal
detector model 加上冻结的 RQL-TV reconstruction procedure 无法解释独立审计数据。
不要声称它是纯 likelihood-class test。

### F1 最终冻结算法

**1. 物理测量分组和外层拆分** — 对每个 image–seed–condition cell：
Bernoulli/GAM4：一行物理 pattern 为一个 logical group。Hadamard：一对
complementary patterns 为一个不可分 logical group。使用冻结 hash，由
(cell_id, seed, pattern kind, 63, 4) 决定 logical-group permutation。前 80%
logical groups 进入 DEV，其余 20% 进入 AUDIT。同一 cell 的所有 reconstruction
arms 共用完全相同的拆分。AUDIT 不参与：λ_max 计算；η 选择；initializer 选择；
任何路径扩展或 solver 调参。若 AUDIT logical groups 少于 128，记
GOF_STATUS = GOF_UNDERPOWERED，MODEL_FAIL_PREDICTIVE = NA。这主要影响 S4
小尺寸 reference cells，不影响 S2-A。

**2. DEV 上选择各 arm 的 η\*** — 对每个 iterative arm 单独进行，但规则完全相同：
固定 K=5，只在 DEV 内做 grouped 5-fold。固定路径
H = {1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1}。λ_max,arm 只用 DEV
计算：TV_NULL_REL=1e-3；expansion cap 40；log-bisection 26 次；每次路径拟合 60
个 outer solver iterations。对每个 η，λ_TV(η) = η·λ_max,arm。各 arm 用自己的
declared data fidelity 拟合，但统一用 calibrated nonparalyzable renewal residual
评分。第 k fold 使用 d_k(η) = (1/n_k)Σ_{i∈k} r_i(η)²。不再对 raw fold sum 做跨
fold 比较。定义 d̄(η) = (1/K)Σ_k d_k(η)，SE_min = sd_k[d_k(η_min)]/√K。这里的
"SE"只称 CV fold-dispersion heuristic，不称置信区间。冻结 tie rule：η_min =
min-arg_{η∈H} d̄(η)。one-SE 选择：η\* = max{η∈H : d̄(η) ≤ d̄(η_min)+SE_min}。
GOF 不再参与 η\* 的 acceptance set。（D2 把 GOF 写进 E，需在 F1 改写。）

**3. 构造 coherent RQL plugin scene** — 模型充分性检验每个 cell 只做一次，由
RQL 承担；它不是逐 baseline-arm 的失败检验。使用 DEV 上 RQL 的 η_min,RQL 而不是
η\*,RQL。原因是 MODEL_FAIL 不应由 one-SE 规则故意增加的 smoothing bias 触发。
固定 λ_plugin = η_min,RQL·λ_max,RQL,DEV。在 DEV 上从冻结 initializer 出发，以
RQL data fidelity 拟合 n_audit = 25 个 outer iterations，得到一张 coherent 图像
x̂_plugin。预测 untouched AUDIT：λ̂_i = Φ a_iᵀ x̂_plugin + d̂。计算
D_obs = Σ_{i∈AUDIT} {N_i−μ_NP(λ̂_i)}²/V_NP(λ̂_i)，以及 r̄_obs、
c_obs = corr(r_i, ρ̂_i)。关键点是：bootstrap generator 必须由一个 x̂ 经过 A 生成。
当前候选代码把各 fold 产生的 cross-fitted λ̂_i 拼接后直接当 synthetic forward
range，因而没有真正重放逆问题。

**4. Refit-per-replicate parametric bootstrap** — 固定 B=39，n_audit=25。对
b=1..39：从 coherent plugin scene 在全部 DEV+AUDIT patterns 上生成
N_i*(b) ~ p_exact-renewal(·| Φ a_iᵀ x̂_plugin + d̂, T, τ̂)。只在 synthetic DEV 上
重新计算 data-dependent initializer。固定：同一个 numeric λ_plugin；不重新选择
η；不重新计算 λ_max；同样的 25 个 outer iterations。

**5. 精确 MC p 值与早停** — [从两个相邻切片衔接] p = (1+#{b: D*(b) ≥ D_obs})/(B+1)，
MODEL_FAIL_PREDICTIVE ⇔ p ≤ 0.025（B=39 → 恰 1/40 水平）。早停：B_used=19，
若前 19 个均未超过，再生成剩余 20 个。该 early stop 与完整 39 次的二元判决完全等价。

**6. 固定-λ̂ bootstrap 只保留下尾泄漏诊断** — 另做便宜的 B₀=199 次
exact-renewal simulation，固定 AUDIT 的 λ̂_i，不重拟合。定义
p_low = (1+#{b: D₀_b ≤ D_obs})/200。仅当 p_low ≤ 0.01 时记
LEAKAGE_SUSPECT = TRUE。它只检测"held-out residual 好得不可信"，不作为上尾
MODEL_FAIL gate，也不改变 η\*。

**7. 残差结构只作独立 warning** — 冻结两个诊断：
MEAN_RESIDUAL_WARN = |mean_r_obs| > max_b |mean_r_b*|；
LOAD_CORR_WARN = |corr_obs| > max_b |corr_b*|。若 correlation 不可定义，记 NA。
warning 不进入 MODEL_FAIL 的并集，以免未校正的多重检验抬高假失败率。论文中照样画
residual-vs-load 图。诊断的作用是帮助解释失败机制。

**8. 最终 reconstruction** — 每个 arm 使用 DEV 选出的
λ_TV,\* = η\*·λ_max,arm,DEV 在全部帧上重新拟合：n_final = 200。MODEL_FAIL：
不改变 η\*；不切换到 η_min；不删除 cell；不查看 PSNR 进行补救；不停止输出
reconstruction；只作为预注册 adequacy flag 进入汇总。当前代码在 E=∅ 时强制选择
η_min，这会让 GOF 反过来改变 estimator，必须删除。

**9. 成本** — 新增 absolute-GOF 层最坏成本为 (1+39)×25 = 1000 个 auxiliary
solver iterations。[相对于……截断]

**10. Continuous S3** — [续] 明确 NotImplemented，但 D2 的 S3 包含 continuous。
由于 afterpulse 延迟为无界指数尾，简单随机 fold 或有限 purge 并不提供严格独立
AUDIT。建议冻结：continuous S3: eta_star = corresponding active-start cell's
eta_star；GOF_STATUS = GOF_NA_DEPENDENT；MODEL_FAIL_PREDICTIVE = NA。
"corresponding"固定为同 arm、image、seed、ρ、ν、A。continuous 用 block folds
有效更诚实。

### select_eta.py 冻结前还会让审稿人皱眉的项目
- 当前 B=30 + empirical quantile / mean±2.5sd band 没有 level 控制。30 次无法
  定义 2.5% Monte Carlo 上尾。
- 当前 refit generator 是拼接的 cross-fitted λ̂_i，不是 coherent A x̂ + d。
- 当前 refit 模式仍把 fixed-λ̂ 的 mean-residual band 当作 hard gate，统计 null
  前后不一致。
- 当前同一个在 η_min 得到的 GOF band 被施加到所有 η。不同 smoothing level 的
  predictive-error null 不同。
- GOF 目前进入 one-SE acceptance set，并在失败后改变 η\*。按本裁决彻底解耦。
- λ_max 目前先在全数据上计算。新算法必须只使用 DEV；否则 AUDIT 已通过 penalty
  scale 泄漏。当前 lam_max_arm 的定义和 full-data 调用路径都需要相应移动。
- bootstrap RNG key 使用 int(round(ctx.T))。主物理单位下 T 为几十微秒，round 后
  大量 cell 都变成 0，导致不同 ν 共用 bootstrap stream。应冻结为整数键：
  (cell_id_hash, nu_integer, tau_ns_integer, …)。

---

## Q2 — Grönberg 2018 全文与更广 prior-art 清扫

交叉审计基础：Wiley 官方文章记录与完整 reference list；KTH thesis summary 和
Paper I 官方条目（kth.diva-portal.org）；同组 2018 SPIE companion；论文引用的
前驱；同组后续 spectral-CT / pileup / DQE 工作；核仪器、SPAD LiDAR、dToF 和
fluorescence/flow 方向的 optimal-rate 文献。（2018 正文 PDF 本身在本轮不可直接
获取 → 留一项人工 sign-off，见文末。）

**2.1 直接裁决 — 结论 1：没有找到他们给出 ν^{1/3} ridge law。**
在 Grönberg 2018 及其 companion 中，没有找到 ρ\*(ν)~Cν^{1/3}、
ρ\*(ν)=(6ν)^{1/3}−2/3+O(ν^{−1/3})、或 I_N = E[N]−ρ²E[Var(R_ν|N)] 中的任何
结果；也没有把 ν=T/τ 作为第二轴求 integrated-count Fisher information 的
principal ridge。Wiley 对论文方法的描述：推导 nonzero-pulse-length
nonparalyzable detector 的 analytical count distribution；Monte Carlo 验证；
通过 Fisher information 计算 spectral-CT image-performance metrics；比较完整
分布、Gaussian approximation 与 ideal nonparalyzable model。SPIE companion 的
数学范围更明确：renewal theory 推导渐近均值和方差，由这些 moments 构造
count/pileup correction estimator 并做 Monte Carlo 验证。这与 ROUND63 的问题不同。

**2.2 Grönberg 2018 的分析终点** — 其贡献链：finite-pulse detector physics →
analytical output-count distribution → Gaussian approximation → task Fisher
metrics。它处理的是比我们更真实的：非零 pulse length；energy-resolved count
bins；spectral-CT material parameters；detector performance/image task。官方
摘要明确说 Gaussian approximation 可以准确表示其分布和 performance metrics。
companion 停在 E[N_T]、Var(N_T) 的 large-T renewal formulas 和 pileup
correction。ROUND63 [DLP-HOLE ~200 chars：ROUND63 侧贡献对照清单 — ideal
active-start nonparalyzable；每个 exposure 只保留一个 integrated scalar count；
exact finite-window PMF；θ=logλ；两个显式无量纲轴等] …terminal-phase missing
information。这不是"我们有 exact PMF、他们没有"；Grönberg 的 count model 本身
也很强。真正的区别是：他们把 count law 用于 task evaluation；你们把
finite-window integrated-count FI 本身作为一个双尺度渐近优化对象。

**2.3 他们引用的先例** — Müller / Yu–Fessler：dead-time count law、均值、方差
和 correction 基础（Grönberg reference list 同时引用了它们和 Alvarez）。没有
找到这些工作给出 (ρ,ν) 双轴 principal FI ridge、ν^{1/3}、terminal residual
dead-time phase 的 missing-information decomposition。
**Wang 2011**（比 Grönberg 更危险）：确实将 CRLB 写成 input count-rate 函数；
显示 pileup 导致 task performance 随 rate 恶化；使用 ideal nonparalyzable
delta-pulse model 作为解析主臂。但它优化的是 energy-discrimination imaging，
没有 integrated count 对 logλ 的 finite-window FI；也没有导出 ρ\* ∝ ν^{1/3}。
**Alvarez 2014**：CRLB 核心使用 multivariate Gaussian count model；在 typical
large-count regime 中论证 covariance-derivative Fisher 项可忽略。它和我们的
CLT proxy 是近邻，但没有进入"finite count ceiling 导致 exact count-only FI
最终回落"的渐近层。

**2.4 更广的 information-optimal rate 先例** —
核仪器：**Bécares–Blázquez 2012** 已研究 detector 对高强度 spallation source
的 optimal counting rate；最优 rate 依赖 source multiplicity、dead-time
calibration uncertainty、measurement duration；代表条件下最优值约 1.05×10⁵
counts/s。影响："低率 shot noise 与高率 dead time 折中产生内部最优点"不新；但
其目标是 dead-time-corrected rate 的相对不确定度且含 τ 标定误差；不是 exact
count-only Fisher ridge，也没有 ν^{1/3}。
SPAD LiDAR：**Gupta 等 Photon-Flooded** 直接求最优 flux；closed-form optimum 由
ambient background 与 pileup distortion 决定，而不是 finite-window
integrated-count terminal phase。**Rapp 等 high-flux lidar** 使用完整
detection-time sequence 的 Markov 模型，在远高于传统 5% rule 的 flux 下恢复
timestamp sequence，因而不会遭受 count-only terminal-phase 信息损失。
**Jorgensen–Johnson 2026**：dead-time event-detection process 的渐近统计理论
（sufficient statistics、histogram 信息损失、Fisher-information rate、MLE/
one-step 渐近效率），研究周期性带 gating 的 event sequence，不是单窗口
integrated count 的 finite-ν ridge。
**2025 dToF CRLB** 工作：在含 dead time/pileup 的 CRLB 下求 optimal photon
flux、pulse width、quantization resolution —"最优 flux"宽泛故事已经拥挤。
另有 MDPI 文献给出 pileup 条件下的最优 optical photon flux。所以不能写
"Information-optimal dead-time operation has not been studied"；可以写的是更窄的
"Prior work established task-dependent finite-flux optima and
count-rate-dependent CRLB degradation; here we derive a finite-window
count-only information ridge and its cube-root scaling."
**Flow cytometry / fluorescence**：10% coincidence 工程规则；高 flux 非线性区与
lifetime bias；SPAD array 并行化；FFS 中 dead time 使 photon-count histogram
变窄并系统性扭曲参数估计。没有同构 ridge。

**2.5 最终 novelty 判词** — 能守住：截至本轮最强清扫，没有找到先例明确给出
ρ\*(ν)=(6ν)^{1/3}−2/3+O(ν^{−1/3})，也没有 I_N=E[N]−ρ²E[Var(R_ν|N)] 作为
active-start nonparalyzable integrated count 的 terminal-phase
information-loss identity。因此 "No prior cube-root finite-window
count-information ridge was found" survives this search。但必须限定为：ideal
nonparalyzable；active-start；scalar integrated count；exact finite-window FI
about log-rate；principal/global ridge；(ρ,ν) asymptotics。
不能守住（都已有先例）：dead time 下存在最优 count rate；Fisher/CRLB 在高 flux
恶化；高通量并非越高越好；正确建模 dead time 可显著缩短 acquisition；detector
information 取决于保留 event times / histograms / aggregate counts。

**建议论文 novelty 段安全版本**：
"Previous studies established count-rate-dependent CRLB degradation,
task-specific optimal operating fluxes, and high-flux recovery using
time-resolved dead-time models. Here, we study a different information
bottleneck: a finite exposure from which only the integrated nonparalyzable
count is retained. For this count-only observation, we derive a principal
finite-window Fisher-information ridge and show that its optimal load obeys a
cube-root law, ρ\*(ν)=(6ν)^{1/3}−2/3+O(ν^{−1/3}), arising from terminal
detector-phase information loss."
不要写 "first"，也不要写 "optimal dead-time flux has not previously been studied"。

**Q2 最终等级**：cube-root ridge 的 novelty 暂评"高置信可守"；但由于 2018 正文
PDF 未能直接全文核对（经由 KTH thesis bundle），需在 S5 做最后一次人工全文搜索
（关键词：Fisher, information, count rate, maximum, optimal, measurement time,
dead time, asymptotic）并对 Paper I 的 Fisher 公式和全部图轴做 provenance
sign-off。
