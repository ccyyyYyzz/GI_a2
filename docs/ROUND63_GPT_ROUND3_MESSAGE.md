# 第三轮消息（发给 GPT Pro，chat 候选方案评审与建议）

第三轮。两件事：规范审计 + 一个理论深挖（有新数值结果）。仓库已更新，请先拉最新 main（github.com/ccyyyYyzz/GI_a2）。

【1. 规范审计】docs/ROUND63_HIGHFLUX_CAMPAIGN_SPEC_DRAFT.md 是按你第二轮全部修订写成的预注册战役规范草稿（解耦 (Phi,tau,T)、rho/nu 主轴、DMD 可实现图样含互补对双曝光计费、Liu-2021 预校正基线、六件套失配消融、equal-time/equal-photon 双制度、time-to-quality 主纵轴、四联主图）。请逐节审计：(a) 有无你第二轮要求被遗漏或写歪的；(b) 网格规模 vs "OE 可过审"最小充分集——哪些格子可以砍掉不伤论文，哪些绝不能砍；(c) 分期（S1 本地试点 -> S2/S3/S4 Colab 分片）有无预注册漏洞（试点后冻结协议再跑全网格，是否需要把"试点可调、全网格不可调"的边界写得更死？）。

【2. 理论深挖：最优负载脊线的立方根标度律】我们实现了精确 renewal 计数似然（code/round63/physics.py, exact_logpmf, active-start 约定 S_m=(m-1)tau+Gamma(m,lam)）并数值计算了精确计数 Fisher 信息（关于 log-lambda，每单位时间，tau 单位）：I_log(rho;nu)。结果（results/round63_theory/fisher_ridge.json）：
- CLT 代理 rho/(1+rho) 单调升；但精确信息在内点达峰后衰减（N 被钉在 nu=T/tau 时信息趋零）；
- 脊线数值：rho* = 4.53/6.16/7.87/9.99/13.77/17.45/22.16 @ nu=20/50/100/200/500/1000/2000；
- 拟合发现 rho* ≈ 1.7 * nu^(1/3)（log-log 斜率 ~0.345）。启发式推导：CLT 破裂发生在 Var N = rho*nu/(1+rho)^3 ~ O(1)，大 rho 下 (1+rho)^3/rho ≈ rho^2 => rho ~ (c*nu)^(1/3)... 但这给出的是"方差~O(1)边界"而非信息峰位置，两者常数不同。
请你：(a) 从精确 renewal 似然出发解析推导 rho*(nu) 的标度与常数（大 nu 渐近；信息峰的正确判据是什么——是 Var N ~ O(1) 还是别的量，比如 pmf 支撑宽度与均值灵敏度 dmu/dlogrho 的竞争？dmu/dlogrho = nu*rho/(1+rho)^2，pmf 宽度 ~ sqrt(rho*nu)/(1+rho)^{3/2}，信噪 = 灵敏度/宽度 => I_log ≈ nu*rho/(1+rho) ... 检验：这个比值平方 = nu*rho/(1+rho)（正是 CLT 公式），峰值来自离散性/边界效应的修正项——请把修正项算出来）；(b) 判断这个标度律 + 脊线图作为论文的一个理论小节（或独立短文？）的价值与最贴近先例（dead-time 计数统计的 Fisher 信息文献：Müller、Rapp&Goyal 等有没有已给出 count-only Fisher 信息的内点峰/标度律？若有请给出与我们的差异）；(c) QMLE 信任边界（exact/CLT 比值 0.9）数值在 rho=5(nu=20) 到 53(nu=2000)，全部高于部署区 rho in [0.3,1]——把"QMLE 有效区/exact 专属区"的分区写法给一个你认为审稿人最舒服的表述。

顺带：物理核心冒烟数据（均值/方差与 CLT 公式 0.2% 内、pmf TV 距离 0.003、paralyzable 均值符合 lam*T*exp(-lam*tau)、active/delayed 恰差 0.5 计数@rho=1）在 REPORT 之外的 commit 信息里，可信其已验证。
