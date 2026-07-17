# ROUND59 理论手册：分布伴随鬼影成像的完整数学物理背景（执行会话必读）

**日期：** 2026-07-17
**写入者：** Fable5
**读者：** 在 `D:\GI_another`（GitHub 仓库 GI_a2）执行 Stage-0 的 Claude Code 会话
**配套文件：** `ROUND59_STAGE0_SCORE_ADJOINT_GI_FINAL_EXECUTION_SPEC.md`（执行规范；先读本手册再读它）

---

## 0. 你在做什么

本项目为鬼影成像（ghost imaging, GI）寻找一个旗舰级理论候选。候选经过 30+ 轮敌意评审与两轮外部审计（Fable5 + GPT Pro）后存活，现在需要你执行预登记的 oracle Stage-0 判决实验。**候选的科学主张**：

> 复杂散斑照明分布的概率梯度场（score）是非线性鬼影成像的正确统计伴随参考；传统 Gaussian GI 是其线性-score 特例；该伴随关系随散斑模式数 $k$ 呈三区间定律（$k=1$ 边界失效 / $1<k\le2$ 方差无限 / $k>2$ 正规收敛）。

你的实验不是"让方法看起来好"，而是**裁决**：若关键门失败，如实报告 KILL——负结果同样是本项目的合格产出。三条纪律：(a) 门槛预登记、不可调参择优；(b) `GAM1/GAM2` 臂的失效是**理论预测**，是验收对象，禁止当 bug 修掉；(c) 失败即停、原样报告。

## 1. 物理设定与符号

**计算鬼影成像。** 一束结构化/随机光场照明未知物体，物体反射/透射的总光强由一个**无空间分辨的桶探测器**（bucket）收集。第 $i$ 帧：

- $a_i\in\mathbb R_{\ge0}^n$：物面照明强度图样（$n$ 像素，展平）；
- $x\in\mathbb R_{\ge0}^n$：物体反射率图（未知，待重建；本实验归一化 $\sum_jx_j=1$，使 $U_0=\mathbb E[\langle a,x\rangle]=1$）;
- $u_i=\langle a_i,x\rangle$：帧能量投影（单指标）；
- $b_i$：bucket 记录，$\mathbb E[b_i\mid a_i]=\phi(u_i)$，$\phi$ 为单调标量响应（可能未知、非线性），噪声为 Poisson/读出/renewal。

**为什么图样是随机散斑。** 陈文组（PolyU）场景：光经过动态复杂散射介质（漫射器、湍流、多层散射）后到达物面，图样是散斑。**参考臂**（分束器+相机）在不接触物体的情况下记录 $a_i$ 的样本——这就是"照明分布 $p_A$ 物理可采样"的仪器基础，也是学习 score 的合法数据来源（Phase C）。

**散斑统计的物理词典（本实验照明族的物理含义）：**

- 单模完全发展散斑：复场为圆高斯，强度逐点服从指数分布 = $\mathrm{Gamma}(k{=}1)$。这是转动漫射器 GI 的教科书照明。
- $k$ 模散斑：曝光时间 $T$ 内介质经历 $k\approx T/\tau_c$ 次去相关（或多偏振/多波长模式求和），强度趋于 $\mathrm{Gamma}(k)$，$k$ 即"散斑模式数"。**$k$ 是可从参考臂边缘矩直接估计的光学参数**。
- 弱起伏湍流：强度近似对数正态 → `CORR-LOGN` 族（空间相关、正值、非高斯）。
- 动态散射状态切换：照明 ensemble 是多状态混合 → `MIX-LOGN` 族。

## 2. 经典估计器谱系与两条精确定理

**相关 GI（1995 至今的标准做法）**：

$$\hat x_{\rm GI}\propto\frac1M\sum_i(b_i-\bar b)(a_i-\bar a_i)\ \xrightarrow{M\to\infty}\ \operatorname{Cov}(a,b).$$

**定理 1（白化的分布无关无偏性）。** 线性链路 $b=\langle a,x\rangle+\varepsilon$（$\mathbb E[\varepsilon\mid a]=0$）下，对**任意**有限二阶矩分布：

$$\operatorname{Cov}(a,b)=\mathbb E[(a-\mu)(a-\mu)^\top]x=\Sigma x\ \Longrightarrow\ \Sigma^{-1}\operatorname{Cov}(a,b)=x\ \text{精确成立}.$$

推论：**非高斯性本身不破坏线性 GI**（只影响方差），帧独立随机增益 $g_i$ 也只改全局尺度。这就是执行规范 A2"线性诚实检查"的理论依据：线性区 score 没有偏置优势可言，任何增益只能标 `LINEAR_EFFICIENCY_ONLY`。

**定理 2（Brillinger 1982，高斯特例）。** 若 $a$ 高斯，则对任意光滑 $\phi$：$\operatorname{Cov}(a,\phi(\langle a,x\rangle))=\mathbb E[\phi'(u)]\,\Sigma x$。即**经典 GI 能容忍未知非线性，靠的是高斯性**。当散斑非高斯且 $\phi$ 非线性时，$\operatorname{Cov}(a,b)$ 混入高阶累积量项，不再与 $x$ 成比例——这是候选要修复的精确失效点。

**其他传统臂的身份**：DGI = 用总照明能量 $r_i=\sum_ja_{ij}$ 做归一的差分相关；correspondence GI（`CORR`）= 按 bucket 高低选帧求图样均值差——**在统计学里这正是两片 sliced inverse regression**（详见 §7）。

## 3. 核心恒等式（弱形式）及其两个配套恒等式

定义照明分布的 score field：$s_A(a)=\nabla_a\log p_A(a)$。

**Stein/分部积分推导。** 设 $b=\phi(\langle a,x\rangle)+\varepsilon$，$\mathbb E[\varepsilon\mid a]=0$：

$$-\mathbb E[b\,s_A(a)]=-\int\phi(\langle a,x\rangle)\,\nabla p_A(a)\,da
=\underbrace{\int p_A(a)\,\nabla_a\phi(\langle a,x\rangle)\,da}_{=\ \mathbb E[\phi'(u)]\,x}\ -\ \text{边界项}.$$

**当且仅当边界项消失**（$p_A\phi$ 在支撑边界为零），得

$$\boxed{\,-\mathbb E[b\,s_A(a)]=\kappa_x\,x,\qquad\kappa_x=\mathbb E[\phi'(u)]>0\ (\phi\text{ 单调增}).}$$

即：不知道协方差、不知道 $\phi$、不知道完整分布，只要有正确的 score，一次广义相关就恢复图像方向；剩下的是无害全局尺度。**高斯特例**：$s_A(a)=-\Sigma^{-1}(a-\mu)$，恒等式退化为白化相关——经典 GI 是本候选的特例，这是论文级叙事的骨架。

**配套恒等式一（score 零均值）**：边界项消失时 $\mathbb E[s_A]=\int\nabla p_A=0$。这使 bucket 中心化合法（LOO 中心化 $b_i-\bar b_{-i}=\frac{M}{M-1}(b_i-\bar b)$，与普通中心化精确成比例，方向不变——按规范实现并注明等价性）。GAM1 的 $\mathbb E[s]=-\mathbf 1\ne0$ 正是边界违反的另一表现。

**配套恒等式二（witness，无物体证书）**：对任意冻结探针 $v$，

$$-\mathbb E[\langle a-\mu,v\rangle\,s_A(a)]=\int\nabla_a\langle a-\mu,v\rangle\,p_A\,da=v.$$

它不需要物体、不需要 bucket，直接检验"手里的 score 是不是这个分布的 score"。探针中心化（减 $\mu$）不改期望但大幅降 DC 方差。iid Gamma($k$, rate $k$) 下其蒙特卡洛误差有解析公式

$$\frac{\mathbb E\|\hat v-v\|_2^2}{\|v\|_2^2}=\frac{k(n+1)}{(k-2)M_w}\equiv q_k^2,$$

这是 A1 门槛 $\mathrm{median}\le1.5q_k$、$P_{90}\le2q_k$ 的来源（UT5 数值复核）。非 gamma 族用样本插值 SE。

## 4. 三区间定律（本论文的核心光学结论，推导必须吃透）

**区间一：$k=1$，边界项不消失，恒等式失效。** 指数密度 $p(a_j)=e^{-a_j}$ 在 $a_j=0$ 处取最大值 1。一维分部积分显式给出破缺项：

$$\int_0^\infty g\,p'\,da=[gp]_0^\infty-\int_0^\infty g'p\,da=-g(0)p(0)-\mathbb E[g'],$$

即 $\mathbb E[g\,s]=-\mathbb E[g']-g(0)p(0)$，比 Stein 恒等式多出 $-g(0)p(0)\neq0$。此时 $s_A(a)\equiv-\mathbf 1$（常向量，不含任何空间信息），于是：**未中心化**估计器 $-\frac1M\sum b_is(a_i)=\bar b\,\mathbf 1$ = 常数图；**中心化**估计器 $-\frac1M\sum(b_i-\bar b)s(a_i)\equiv 0$ = 精确零图（因 $\sum_i(b_i-\bar b)=0$）。两种塌缩都是理论预测（A1 的 GAM1 验收项），而同一设定下白化 GI 精确恢复 $x$——**score 不是普适更优，它有清晰的失效域**。

**区间二：$1<k\le2$，恒等式成立但估计不可用。** $k>1$ 时密度 $\propto a^{k-1}$ 在 0 处消失，边界项消失、期望恒等式恢复。但 $s_j=(k-1)/a_j-k$（rate $=k$）的二阶矩：用 $\mathbb E[a^{-1}]=\frac{k}{k-1}$、$\mathbb E[a^{-2}]=\frac{k^2}{(k-1)(k-2)}$，

$$\mathbb E[s_j^2]=(k-1)^2\mathbb E[a^{-2}]-2k(k-1)\mathbb E[a^{-1}]+k^2=\frac{k^2}{k-2}\ \xrightarrow{k\to2^+}\ \infty.$$

$k\le2$ 时朴素平均无 $M^{-1/2}$ 收敛（重尾、样本均值被 $a_j\to0$ 的帧主导）。**实现警告**：对 $a$ 钳位（如 `1e-6`）会把无限方差藏起来同时破坏 oracle score 的精确性——规范禁止钳位；GAM2 的验收是**展示**发散签名（UT6：二阶矩随截断阈值单调增长不饱和；A1：err 不按 $M^{-1/2}$ 收敛、$P_{90}/\mathrm{median}$ 比值异常大）。

**区间三：$k>2$，正规。** 有限方差、$M^{-1/2}$ 收敛、witness 误差满足 $q_k$ 公式。

**光学解读（写进论文的话）**：单帧单模散斑不可用；时间积分/多模散斑（$k>2$，即曝光超过两个去相关时间或有多个独立模式）进入正规区。有效性边界由一个**可测的光学参数**（散斑模式数）控制——这使定律成为光学结论而非纯统计注记。

## 5. 各照明族 oracle score（实现细节）

| 族 | score | 数值要点 |
|---|---|---|
| GAUSS | $-\Sigma^{-1}(a-\mu)$ | 预分解 $\Sigma$（Cholesky），解算不显式求逆 |
| GAM$k$ (iid, rate $k$) | $s_j=(k-1)/a_j-k$ | float64；**禁钳位**；下溢帧丢弃并计数报告 |
| CORR-LOGN | $-\operatorname{diag}(1/a)\big[\Sigma_z^{-1}(\log a-m)+\mathbf 1\big]$ | $z=\log a$；$\Sigma_z$ Cholesky 复用；$m=-\sigma_{\ln}^2/2$ 保单位均值 |
| MIX-LOGN | $\sum_cw_c(a)s_c(a)$，$w_c=\mathrm{softmax}_c[\log\pi_c+\log p_c(a)]$ | $\log p_c$ 经三角解算的二次型 + logdet；softmax 必须 log-domain 防下溢 |

混合 score 推导：$\nabla\log\sum_c\pi_cp_c=\frac{\sum_c\pi_cp_cs_c}{\sum_c\pi_cp_c}$。注意：混合权重是**后验软分配**——score 自动做无标签状态识别，这正是它对陈组动态散射的卖点；与之公平对照的经典臂是"无监督聚类+分簇白化"（`CLUSTER-WHITEN`）。

## 6. 链路物理：为什么 DT30 是旗舰非线性臂

**死时间物理。** 单光子探测器每记录一个光子后有死区 $\tau$（non-paralyzable：死区内到达的光子被丢弃且不延长死区）。记录事件间隔 iid $=\tau+\mathrm{Exp}(\lambda)$（死区 + 无记忆等待），故长窗均值计数率

$$R=\frac{\lambda}{1+\lambda\tau},$$

单调饱和于 $1/\tau$。取 $\tau=\frac{3}{7s}$ 使 $\tau\lambda_{\rm mean}=3/7$，均值响应压缩至 $1/(1+3/7)=0.7$（30%）。有效链路 $\mathbb E[b\mid a]\approx u/(1+\tfrac37u)$：单调单指标 ✓，$\mathbb E[\phi']>0$ ✓——弱形式恒等式只依赖条件均值结构，对 renewal 噪声照常成立（噪声只进方差）。

**为什么不能用 $\mathrm{Poisson}(\mathrm{sat}(u))$ 冒充**：死时间后的计数方差低于 Poisson（Fano<1），把饱和均值套 Poisson 噪声是物理错误的分布；审计（GPT 修订 12）要求旗舰非线性臂用**精确 renewal 模拟**。实现：逐帧生成间隔 $\tau+\mathrm{Exp}(\lambda_i)$ 的累计和，计数落入 $[0,T{=}1]$ 的个数；$\lambda=su$ 大时每帧约 $0.7s$ 个事件——**按帧分块**（如每块 1000 帧）控制内存（$s=10^4$ 时整批间隔矩阵约 1.2 GB，分块后 <100 MB）。UT3 验证模拟均值 vs $\lambda T/(1+\lambda\tau)$ 相对差 ≤1%。

数学应力臂（SAT30/SAT50/GAMMA07/LOG + Poisson）保留在扩展表：它们检验恒等式对链路形状的普适性，但不承担"物理诚实"的举证责任（B2 只认 DT30）。

## 7. 对照估计器的统计原理（实现时避免把预期现象当 bug）

- **SIR-10/20**（Li 1991）：把 $b$ 分片，片内标准化图样均值 $\mathbb E[\tilde a\mid b\in\text{slice}]$ 的加权外积取主方向。其无偏性依赖**线性条件**（设计分布椭圆，如高斯）；gamma 乘积与 lognormal **不是椭圆分布**，SIR 在这些臂上可以有系统偏差——这是预期现象、是 score 合法胜出的空间，报告即可，不要"修"。`CORR` 是它的两片粗化版。
- **L-Isotron**（Kakade et al. 2011）：交替 [PAV isotonic 拟合 $\hat\phi$] 与 [梯度更新 $x$]，是未知单调单指标的正统经典算法，也是"深度网络是否必要"的关键碰撞控制。协议已冻结（200 轮/早停/3 初始化/留出似然选优）——不得用真值挑初始化。
- **MAVE-16**（Xia et al. 2002）：局部线性平滑的最小平均方差估计，计算重，仅 16×16。它与 SIR、L-Isotron 一起构成"占据统计层"的强制对照——**候选的价值必须在超过它们的部分**。
- **MLE-OR**：已知真 $\phi$ 的似然极大化 = 信息上界诊断。它回答"score 方法离满信息还差多少"，**不进**任何头部空间门的对照集合（不公平地知道 $\phi$）。
- **WHITEN-LW vs WHITEN-OR**：$n=4096$、$M=2\times10^4$ 时样本协方差奇异，Ledoit–Wolf 收缩是可部署版本；WHITEN-OR（真 $\Sigma$）进 `RELATION_HEADROOM` 门以消除"oracle 信息不对称"（score 用了真分布，白化也必须给真 $\Sigma$）。

## 8. 指标协议为什么这样定

理论只恢复方向（正尺度规范自由）。旧版用"真值 LS 标量拟合"定尺度——审计判定这向部署输出注入真值信息、且掩盖幅度失真。主协议改为**通量匹配**：非负截断后把 $\hat x$ 的总和缩放到真值总和（真值已归一 $\sum x=1$）。物理上总通量可由 DC/均值照明标定近似测得，部署合法。真值 LS 版本降级为 `DIRECTION_ONLY_DIAGNOSTIC` 附表，另报尺度不变角误差与 Pearson 相关。**所有方法同一协议**——公平性高于一切。

## 9. 先例地图与措辞纪律（影响 REPORT 的每一句话）

统计层已被占据，禁止任何"首次"式表述：Härdle–Stoker 1989 平均导数估计（ADE）= "样本估计密度 score × 响应相关恢复未知链路方向"；Brillinger 1982 = 高斯特例；Yang et al. 2017 = 高维非高斯 Stein 单指标；Li 1991 SIR、Xia 2002 MAVE/OPG、Kakade 2011 Isotron = 直接算法先例。**可守的只有光学层**：(i) 参考臂散斑 ensemble 使 score 物理可采样可构造；(ii) Gaussian GI 与陈组 Gaussian-constraint 校正被统一为特例/投影；(iii) 三区间定律把有效边界交给可测光学参数（模式数 $k$）。REPORT 中的结论句必须写成"在本仿真协议下……"，禁止外推到未测场景。

## 10. 结果解读决策树（每个门的语义）

- **UT 失败** → 实现错误，修复后重跑（唯一允许"修"的层）。
- **A1 中 GAM3/4/8/CORR-LOGN/GAUSS 不过** → score 实现或恒等式推导有错 → STOP 上报（不是调门槛）。
- **A1 中 GAM1"意外通过"恒等式** → 边界条件理论有错 → STOP 上报（这是重大理论事件，不是好消息）。
- **A2 出现系统方向偏差** → 实现不对称 → STOP 排查；score 线性区效率优势 → 只标 `LINEAR_EFFICIENCY_ONLY`。
- **A3 不过** → 非线性头部空间不存在 → `FLAGSHIP_DEAD`，跑完扩展表供降级理论文，如实收尾。
- **B1 双门**：`RELATION_HEADROOM` 检验"伴随关系本身有信息优势"；`PRACTICAL_HEADROOM` 检验"可部署版本仍有优势"。任一全灭 → `STAGE0_KILL`。
- **B2**：只有 DT30（物理臂）通过才配 `HONEST_NONLINEARITY`；只在强数学臂通过 → `THEORY_ONLY`。
- **全过** → `PROCEED_TO_PHASE_C`：下一会话训练 learned score（只见参考散斑），最终在真实波光学散斑（相位屏+角谱传播+有限 NA+动态多层）上以 witness 为合格证复现头部空间。

## 11. 工程备忘

float64 全程；图样矩阵按 `(照明族, seed)` 生成一次冻结复用（含 $\hat\Sigma$、白化算子）；大数组分块（DT30 间隔矩阵、lognormal Cholesky 解算）；witness 用 $M_w=2\times10^5$ 且不算 bucket；LPIPS 缺失时如实标 `UNAVAILABLE` 并把相关门记"待补"；每 Phase 结束 commit+push；预算超限中断报告。随机性：`np.random.default_rng(SEED0+seed)`，任何新随机源必须挂在这个体系下。

## 12. 一句话总纲

你在检验一条恒等式的三区间定律和它的非线性头部空间：**让失效臂按预言失效，让通过臂在公平对照下自己说话**——两者同样是成果。
