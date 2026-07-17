# ROUND63 战役规范（草稿，待 GPT Pro 第三轮审计后冻结）：死时间感知的高通量单像素成像

**日期**: 2026-07-17 晚
**状态**: DRAFT — 待 GPT Pro 规范审计（第三轮）后按修订冻结
**目标产出**: 一篇可直接投 **Optics Express** 的完整论文（仿真-only，陈组体裁 + 本程序的可复现性纪律）
**工作题目**: *Dead-time-aware high-flux single-pixel imaging: operating beyond the conventional photon-counting regime with renewal-statistical reconstruction*
**上游**: GPT Pro 两轮对抗评审（chat 6a5a66fb，2026-07-17）的全部修订已并入本规范；ROUND59/62 harness 沿用其可复用部分（图像协议、RNG 体系、审计纪律），物理参数化**重写**（见 §2）

**主张纪律（不可违反）**：
- 禁止 "first/首次"；禁止把方法叫 "exact renewal MLE"——正式名称 **Gaussian-renewal quasi-likelihood (QMLE)**，精确 renewal 似然只作小规模参考上界
- 不主打 "+X dB"，主打 **time-to-quality 加速**与**工作点相图**
- equal-time 与 equal-photon 双资源制度并报（后者为诚实负对照：固定入射光子数下死时间只丢事件，不得声称优于理想探测器）
- 复现 Liu-2021 式预校正作为基线（paralyzable 数据上用其原式；non-paralyzable 主臂用匹配的均值反演 λ̂=r/(1−rτ)）
- 逐条如实报告；速度加速若实测 < 理论投影（4–8× vs 5–10×），按实测写

---

## 1. 科学主张（预注册）

单光子单像素成像为回避死时间失真惯于在低负载区（ρ=λτ≲0.05）工作。我们证明（在本仿真协议下）：**同时建模死时间的非线性均值压缩与 renewal 亚泊松色散后，中等负载区 ρ∈[0.3, 1] 成为可用的高速工作区**——在光源功率可提升、采集时间受限的制度下，达到同等图像质量所需的光学积分时间比传统安全工作点缩短 S_Q 倍（理论上界 S_FI(ρ;ρ₀)=[ρ/(1+ρ)]/[ρ₀/(1+ρ₀)]；ρ₀=0.05→ρ=1 时 10.5×，实证预期 4–8×）。三项同时展示的能力：(i) 均值+色散联合重建；(ii) 采集时间加速相图；(iii) 标定误差/探测器类别/背景/afterpulsing 失配下的稳健工作区间。

**理论深挖交付物（与 GPT 第三轮确认）**：Gaussian-renewal 信息率 I_u/T = Φ/[u(1+ρ)] 随 Φ→∞ 饱和于 1/(u²τ)，但**精确离散计数信息在 Var N ≲ O(1)（N 被钉在 ν=T/τ）时必然衰减到 0**——CLT 近似在 λT/(1+ρ)³ < O(1) 处失效。论文相图叠加**精确 renewal Fisher 信息的最优脊线 ρ*(ν)**，把"高通量不是越高越好"的物理写成定量结论。这也自然划定 exact-likelihood 相对 QMLE 的价值区（低 ν / 高 ρ 角落）。

## 2. 物理模型与参数化（与 ROUND59 的关键差异）

**解耦参数化**（ROUND59 的 λ=su、τ=(3/7)/s 把 λτ 钉死，不可用于本战役）：
- 探测器：死时间 τ 固定为物理常数（主档 τ=50 ns；敏感档 25/100 ns，取真实 SPAD 规格范围）；non-paralyzable 为主臂，paralyzable 为消融臂
- 光源/物体：到达率 λ_i = Φ·⟨a_i, x⟩ + d（**暗计数/背景 d 在死时间之前进入**，与信号共同占用探测器）
- 帧窗 T 独立可变；无量纲主轴 **ρ̄ = τ·E[λ_i]**（同时报告逐 pattern ρ 的 5/50/95 分位）与 **ν = T/τ**
- bucket：b = N/s_ref + 读出噪声（读出档位按真实器件设 σ_r；s_ref 仅为数值缩放）
- 帧起始约定三种：active-start（首事件 Exp(λ)）/ delayed-start（ROUND59 约定，τ+Exp）/ continuous（探测器状态跨 pattern 延续，guard g/τ∈{0,1,5}）——主臂 active-start（标准物理），其余为消融
- afterpulsing：branching 模型（有效雪崩以 p_ap 产生指数延迟余脉冲，余脉冲亦触发死时间），p_ap∈{0,0.5,1,2,5,10}%
- 死时间抖动 CV(τ)∈{0,2,5,10}%（次级消融）

**似然（估计器侧）**：
- Gaussian-renewal QMLE：μ_N=λT/(1+λτ)，v_N=λT/(1+λτ)³（+读出方差；objective 按帧数归一）
- 精确 renewal（小规模参考）：P(N≥m)=Γ CDF(T−(m−1)τ−?；按 active-start 约定推导并与生成器逐约定对齐)，读出噪声用离散卷积 p(b|λ)=Σ_m P(N=m|λ)·N(b; m/s_ref, σ_r²)
- 全部 log-CDF-difference 数值稳定实现；优化器 L-BFGS-B（逐图独立、KKT 相对残差 ≤1e-6 出版口径、双求解器对照、多起点一致性）

## 3. 网格

**S2 核心相图（主图数据）**：ρ̄∈{0.01,0.03,0.05,0.1,0.2,0.3,3/7,0.6,1,2,3} × ν∈{20,50,100,500,2000}（主文固定 ν=500 切片，附录全相图）× 图样{Bernoulli-50%, Hadamard 互补对（按 2 次曝光计费）, GAM4 应力臂} × 分辨率{64²,128²} × M/n∈{0.25,0.5,1.0} × ≥5 seeds × 图像集（§4）
**S3 失配地图**：τ̂=(1+δ)τ, δ∈{−20,−10,−5,0,5,10,20}%（+flat-field 标定部署版 + τ 联合 profile 版）；paralyzable/non-paralyzable 2×2；d/E[信号]∈{0,0.01,0.05,0.1,0.25,0.5}（已知/±10% 误差/联合估计）；afterpulsing 档；帧起始三约定；CV(τ) 档
**S4 exact-vs-QMLE 验证套件**：8×8/16×16 × λτ∈{0.03,…,2} × T/τ∈{20,…,2000}，PSNR 差 / exact-NLL regret / 梯度夹角 / Fisher 主方向 / 多起点，预注册门：主区间（ρ∈[0.2,1], ν≥100）median[PSNR_exact−PSNR_quasi]≤0.25 dB 且 95% ≤0.5 dB，不过则 exact=低计数模式、QMLE=高计数快速模式如实分区

## 4. 图像与基线

**图像**：≥24 张自然图（STL-10 test 前 24 张按既有协议灰度化）+ 6 张结构靶（USAF 1951 数字合成、细线组、文字、稀疏点阵、低对比小结构、阶梯灰度）@64²/128²
**基线臂（同一 TV 先验、同一优化器、同一 pattern/到达流耦合随机数）**：GI/DGI；linear-Poisson MLE+TV；saturated-mean Poisson+TV（只修均值不用色散——**机制消融关键臂**）；Liu-2021 式预校正+TV；renewal-QMLE+TV（主臂）；exact-renewal MLE（小规模）；PnP-BM3D（或 DIP，二选一冻结）作为现代先验臂在各数据项上复用
**指标**：PSNR/SSIM/LPIPS + 无重标度辐射 NRMSE + 总通量偏差 + held-out predictive NLL + **主纵轴 S_Q = T_safe(Q)/T_highflux(Q)**（Q∈{25,28,30 dB}+图像自适应目标；PSNR–log T 单调拟合插值，达不到按 censored 处理）+ 墙钟敏感度 T_wall=M(T+t_switch)

## 5. 分期与算力

- **S0** 建库+UT（本地）：physics 核心（主循环亲写）、估计器套件（Opus）、pattern/图像/Colab 分片设施（Opus）、对抗审计 workflow 后冻结
- **S1** 试点（本地 i9）：64²、ρ̄∈{0.05, 3/7, 1}、ν=500、8 图、3 seeds——验证加速信号存在性与量级，校准全网格预算
- **S2/S3/S4** 全网格：**Colab pro2×3 会话 + pro1×2 会话 + 本地 = 6 路并行**，按 (阶段×ρ̄ 分片) 切分；每分片自包含（代码+冻结输入 SHA256 经 Google Drive），keep-alive+watchdog 按既有纪律；结果 CSV 回传合并，META-as-truth 断点
- **S5** 论文（OE opticajnl 模板，Chen 风格指南应用于表述层，主图四联）+ GPT 评审循环 ≥2 轮至可投

RNG：SEED0=20260717 体系沿用，stream 标签 63。全部冻结输入/输出 SHA256 落盘。

## 6. 主图（预设计）

四联横图：**(a)** 入射率-记录率曲线 + 安全区/高速区/失配敏感区着色 + 精确信息脊线 ρ*(ν)；**(b)** 同质量-更短时间图像阵列（1 自然图+1 文字靶 × 5 列：安全全时/安全短时/高通量 naive/高通量预校正/高通量 QMLE，细节放大框）；**(c)** PSNR–总积分时间曲线族（ρ̄ 四档）+30 dB 水平线上的加速箭头 + bootstrap 区间；**(d)** 速度相图：实证 S_Q vs ρ̄ + 理论 ρ/(1+ρ) 上界带 + τ±10%/afterpulse/背景带 + 失效区。
