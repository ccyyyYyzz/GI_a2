# GPT Pro 第三轮回复（原文存档，2026-07-17 深夜）

总裁决

我已按最新 main 核对了规范、physics.py、solvers.py、campaign/shard 设施和 Fisher 输出。结论如下：

项目 裁决
ROUND63 规范草稿 HOLD：当前不能冻结，需出 D2 版
高通量 A 主线 继续，故事比第二轮时更强
立方根脊线 成立；领先常数应为 6
1/3
=1.81712，不是经验的 1.7
“方差降到 O(1)”解释 错误，必须从规范删掉
当前 Gaussian-QMLE 命名 过时；改为 renewal quasi-likelihood / quasi-score
Wedderburn IRLS 修复 方向正确，但原始计数情形有更干净的直接凸目标，无须 IRLS
当前 TV 选择 必须替换；逐臂 own-NLL 选择会系统性污染公平比较
独立理论短文 暂不建议；先作为 A 的强理论小节。取得统一渐近证明后再考虑拆稿

最重要的两点是：

你们发现的内点信息峰和 ν
1/3
脊线是真结果，而且有清楚的统计机制。

当前规范和代码仍有若干会改变结果的实现缺口；这不是文字润色问题。

1. ROUND63 规范逐节审计
1.1 §0–§1：科学主张与命名
A. 当前名称已经写歪

规范仍把主臂称为 Gaussian-renewal QMLE，但现在实际主估计方程是 Wedderburn quasi-score：只使用均值和方差关系，不把数据假设成服从某个 Gaussian density。规范虽然禁止叫 exact MLE，却尚未反映这次方法变更。

建议正式名称改为：

renewal quasi-likelihood reconstruction，RQL

算法层面可写：

solved by quasi-score IRLS

不要再把主算法称为 Gaussian likelihood。Gaussian renewal approximation 可以作为推导矩函数的来源，但不是最终概率模型。

B. 预注册不能提前写“我们证明”

规范 §1 目前直接说“我们证明中等负载成为高速工作区”，并给出实证预期 4–8×。这在战役完成前应改成待检验假设，而不是结果陈述。

建议冻结成：

We test whether renewal quasi-likelihood reconstruction permits operation in
ρ
ˉ
​

∈[0.3,1] with a statistically significant reduction in total optical integration time relative to the conventional
ρ
ˉ
​

=0.05 operating point.

还缺一个明确的主假设、主高通量点和过门规则。否则跑完以后可以在 ρ=0.3,0.6,1 里选最好看的一个。

我建议：

安全参考点：ρ
0
​

=0.05；

预注册主高通量点：ρ
1
​

=0.6；

ρ=1 为次级强负载点；

两边使用同一个 RQL 估计器，避免把“换模型收益”和“换工作点收益”混在一起；

主目标

Q
j
⋆
​

=PSNR
rad
​

(ρ=.05,ν=2000,j)−1 dB;

主速度比

S
j
​

=
T
ρ=.6
​

(Q
j
⋆
​

)
T
ρ=.05
​

(Q
j
⋆
​

)
​

;

建议门：

24 张自然图的中位 S
j
​

≥3；

分层 bootstrap 的 95% 下界 >1；

至少 18/24 张图 S
j
​

>1。

固定 25/28/30 dB 可以继续作为次级终点。

C. 规范中的理论解释必须重写

规范目前写：

精确信息在 VarN≲O(1) 时衰减，因此得到立方根脊线。

这在代数上不成立。

因为大 ρ 时

VarN≃
(1+ρ)
3
νρ
​

≃
ρ
2
ν
​

.

令其为 O(1) 得到的是

ρ=O(ν
1/2
),

不是 ν
1/3
。

在真正的信息峰

ρ
⋆
≍ν
1/3

处，

VarN≍ν
1/3
,

仍然发散，计数分布还跨越许多个整数。立方根峰发生在“计数被钉死”之前。

正确机制是：

CLT 信息距离饱和值的剩余缺口为 1/ρ，而把完整检测时间序列压缩成一个总计数会丢失窗口末端 detector phase，其信息损失首项为 ρ
2
/(12ν)。二者边际收益平衡产生 ρ
3
=6ν。

后面给出完整推导。

1.2 §2：物理与似然规范

这一节有六个冻结前硬问题。

A. bucket 模型与实现不一致

规范写

b=N/s
ref
​

+readout noise,

并要求 exact likelihood 对读出噪声做离散卷积。

实际代码的主路径是原始数字计数 N，默认没有读出噪声；Gaussian 噪声只是可选 ablation。

对于 SPAD/PMT counter 论文，最干净的选择是：

主模型：直接观察原始计数 N，σ
b
​

=0；

additive readout noise 只做次级消融；

主 exact likelihood 不必背负卷积计算；

只有启用 σ
b
​

>0 的附录臂才用

p(b∣λ)=
m
∑
​

p
m
​

(λ)N(b;m,σ
b
2
​

).

否则规范承诺了一套当前代码完全没有实现的 likelihood。

B. active-start 公式应删除问号并冻结

正确公式是

P(N≥m)=F
Γ(m,λ)
​

(T−(m−1)τ),m≥1.

delayed-start 才是

F
Γ(m,λ)
​

(T−mτ).

当前 physics.py 的公式与 active-start 生成器是一致的，但规范仍保留了 −? 占位。

C. “稳定 log-CDF difference”尚未实现

当前代码实际是：

Python
Run
p = p_ge(N) - p_ge(N + 1)
log(max(p, 1e-300))

这是普通 CDF 相减再截断，不是稳定的 log-CDF difference。Fisher 又使用中心差分，并删除 p<10
−15
的项。

出版版本必须改成：

CDF/SF 分支自适应；

logdiffexp(logG_m, logG_m1)；

解析的 ∂
logρ
​

G
m
​

，不用有限差分；

FI 在 log domain 求和；

不用任意 10
−15
cutoff 决定高通量尾部。

当前 ridge 数字方向可信，但不能直接作为最终理论图。

D. continuous guard 有符号错误

当前实现为

carry=ready−T+g.

这意味着 guard 越长，下一帧反而越不 ready。正确应为

carry
i+1
​

=max(0,ready
i
​

−T−g).

此外，在 continuous 模式中，当前帧末之后发生的 afterpulse 被直接丢弃；它们应作为下一帧事件队列的一部分继续携带。

这两项会直接改变“跨帧 memory”消融的方向，必须在冻结前修。

E. paralyzable/non-paralyzable 2×2 目前只是名义存在

生成器可以生成 paralyzable 数据，PRECORRECT 也有 Lambert-W 反演；但 QMLE、SAT-POISSON 和 EXACT 仍一律使用 non-paralyzable 均值、方差和 count likelihood。

所以现在不能声称完成：

真 detector 假设 detector
nonpar nonpar
nonpar par
par nonpar
par par

你们只能声称：

paralyzable generator 对 non-paralyzable RQL 的模型类别压力测试。

要保留真正的 2×2，必须增加 matched paralyzable mean/variance/quasi-score，且 shared initializer 也不能继续无条件使用 non-paralyzable 预校正。

Liu-2021 的最近邻本来就是 paralyzable detector；忠实复现这一格不能砍。该文已在 DMD+PMT 单光子压缩成像中做了死时间预校正和实验验证。
PubMed

F. 改变绝对 τ 的主网格大多是重复计算

当 ρ=λτ、ν=T/τ 固定，并且没有绝对时间尺度的背景、afterpulse、switch overhead 或电子噪声时，把

τ=25,50,100 ns

分别跑一遍会产生相同的无量纲计数律。

因此：

核心相图固定 τ=50 ns 即可；

25/100 ns 只在保持暗计数 cps、afterpulse delay 秒数、DMD overhead 等绝对量不缩放的敏感性实验中有意义。

这是可以直接砍掉的大块重复网格。

1.3 §3：网格规模与最小充分集

当前 S2 核心网格是

11 ρ×5 ν×3 patterns×2 resolutions×3 M/n×5 seeds=4950

个 condition-seed 单元；再乘 30 张图就是 148,500 个 image-condition 单元，尚未乘方法数、TV 路径、fold、PnP 和 exact。

这不是 OE 的“充分”，而是会把战役变成算力驱动的全因子表。

更严重的是，128
2
、M/n=1 时：

Bernoulli dense float64 A：约 2 GiB；

Hadamard complementary raw A：约 4 GiB；

再加完整 Hadamard 矩阵、求解变量和 fold，标准 Colab 很容易爆内存。

当前 pattern 实现确实显式构造整个 dense matrix。

128² 必须采用：

matrix-free Hadamard；

on-the-fly Bernoulli block；

或分块 A@x / A.T@y operator。

我建议的 OE 最小充分战役
S2-A：主 time-to-quality 实验——绝不能砍
维度 冻结值
分辨率 64²
主图样 Bernoulli-50%
M/n 0.5
ρ {0.05,0.3,0.6,1,2}
ν 建议 {20,50,100,200,500,1000,2000}；最低不能少于 5 个对数尺度点
图像 24 张确认性自然图
seeds 5
全曲线主臂 RQL
全曲线比较臂 POISSON-LIN、SAT-POISSON、PRECORRECT 至少在 ρ={0.05,0.6,1}
GI/DGI 选定曲线点和展示图即可

主 time-to-quality 必须固定 M，改变 dwell T=ντ。

当前 S1 的 Part B 通过改变 M 得到时间曲线，这同时改变：

光学时间；

空间测量秩；

pattern 覆盖；

正则化难度。

因此它不能作为 headline acquisition-time speedup。

改变 M 的曲线可以保留为次级“总系统资源曲线”，但主速度结论应来自固定 pattern 数、只改变每幅 pattern 的曝光时间。

主时间成本必须定义为

T
opt
​

=M
physical
​

T,

而不是目前规范中的容易与单帧 T 混淆的 T
safe
​

/T
highflux
​

。

Hadamard pair 使用

M
physical
​

=2M
signed
​

.
S2-B：欠采样稳健性——缩小即可

M/n∈{0.25,0.5,1}；

ρ∈{0.05,0.6,1}；

ν∈{100,500,2000}；

12 张自然图；

3 seeds；

64² Bernoulli；

四个机制臂。

不需要对全部 ρ×ν 重跑。

S2-C：图样与尺度外推——选锚点

Hadamard/GAM4：

ρ∈{0.05,0.6,1}；

ν∈{100,500,2000}；

M/n=0.5；

8–12 张图；

3 seeds。

128²：

Bernoulli 为主；

ρ∈{0.05,0.6,1}；

ν∈{100,500,2000}；

M/n=0.5；

8–12 张图；

3 seeds；

只跑 RQL、SAT、PRECORRECT 和一个现代先验臂。

PnP 建议冻结为 PnP-BM3D，不要继续保留“BM3D 或 DIP 二选一”。“或”就是分析分叉。DIP 还会额外引入训练停止、初始化和网络结构选择。

S3：失配地图——必须 OAT，不要全笛卡尔积

固定：

64²；

Bernoulli；

M/n=0.5；

ν=500；

ρ∈{0.3,0.6,1}，另加 ρ=.05 参考；

12 图、3 seeds。

推荐：

失配 保留档位
τ 误差 −20,−10,0,10,20%；可砍 ±5%
暗计数比例 0,.05,.1,.25,.5
afterpulse 0,1,2,5,10%；可砍 0.5%
detector class 真正实现后保留 2×2
start mode active、delayed、continuous
guard continuous 下 0,1,5τ
dead-time jitter 0,5,10%；放附录

只额外跑三个交互：

ρ×τ-error；

ρ×p
ap
​

；

continuous × afterpulse。

“flat-field 标定版”要保留。“联合 profile τ”可以只在 1–2 个代表点做；不必全地图。

S4：必须拆成“标量理论”与“图像验证”

当前 S4 图像网格只到 ρ=2，但你们的 Fisher ridge 在 ρ=4.5–22；它无法验证理论脊线。

应拆成：

标量 exact-Fisher map

ν=20,…,2000；

ρ 自适应到至少

max(64,2ρ
⋆
(ν));

不跑图像重建。

图像 exact-vs-RQL

8²/16²；

ρ∈{.03,.1,.3,.6,1,2}；

ν∈{20,100,500,2000}；

3 seeds；

exact 与 RQL 使用相同 TV、相同无真值选择器。

当前 EXACT 是无 TV 的 L-BFGS-B，而 RQL 有 TV，不能被称作 reference upper bound。

1.4 哪些绝不能砍，哪些可以砍
绝不能砍

RQL 在安全点与高通量点的同估计器比较。

SAT-POISSON：它是“只修均值、不用 renewal 色散”的关键机制臂。

Liu-style PRECORRECT。

equal-time 与 equal-incident-photon 双制度。

radiometric 指标，不允许只用 flux-matched PSNR。

τ±10%、背景、afterpulse、detector class、continuous carry 的锚点消融。

exact-vs-RQL 小规模公平验证。

至少一组 128² 的可视展示。

exact Fisher principal ridge。

common truth-free TV 选择规则。

可以砍

64/128 × 三图样 × 三采样率的全因子积；

image grid 中的 ρ=.01,.03,3；

ρ=.03 只留在 exact/低计数验证；

每个失配轴的全组合；

τ=25/100 ns 的无量纲重复格；

±5% τ；

0.5% afterpulse；

jitter 全主图；

PnP 在每个网格点；

30 张图在每个 stress cell；

exact likelihood 在 64²/128²；

所有结构靶进入统计推断。结构靶适合作为展示，不应和自然图一起充当独立统计样本。

1.5 §4：图像、指标与基线
A. 试点图像与确认性图像不能重叠

当前 pilot8 使用 STL 前四张及文字、USAF、细线、低对比图，而最终集合又包含这些图。

这意味着 S1 看到 PSNR 和图像以后，可以调整：

TV 规则；

主工作点；

展示图；

crop；

甚至基线。

建议：

S1 使用 STL-10 train 或独立开发图集；

S2 确认性自然图完全不进入 S1；

结构靶可以另做一套开发版；

主图最终展示的图名和 crop 坐标在 S2 前冻结。

B. radiometric 指标应升为主指标

当前 campaign 的 main_metrics 是 flux-matched，之后才单独算 radiometric NRMSE。

你们刚发现的 full-Gaussian 病灶恰恰是亮度过估。若主 PSNR 在评价前重新缩放总通量，它会把最危险的错误擦掉。

建议：

主 time-to-quality：

PSNR
rad
​

——不对
x
^
做重标度；

同时报告 radiometric NRMSE 和 flux bias；

flux-matched PSNR 仅称 shape-PSNR，作为次级指标；

SSIM/LPIPS 仍可按归一图像计算，但不能替代 radiometry。

C. 结果单位应是图像，不是 image×seed

主统计建议：

图像为外层独立单位；

seed 在图像内配对平均；

bootstrap 先重采样图像，再在图像内部重采样 seed；

不把 24 图 × 5 seeds 当成 n=120 个独立样本。

1.6 §5：分期与预注册漏洞

答案是：必须把“试点可调、全网格冻结后不可调”的边界写得更死。

当前 S1 文件明确写：

verify the acceleration signal and calibrate its magnitude；protocol knobs may be tuned HERE。

这几乎允许在看到真值指标后改任何东西。

建议把以下文字直接写进 D2：

S1 is exploratory and development-only. S1 images and seeds are disjoint from all S2 confirmatory images and seeds. S1 results shall not enter confirmatory confidence intervals or main quantitative tables.

Before the first S2 cell is executed, the following shall be frozen in one immutable commit: the primary hypothesis and pass rule; primary ρ
0
​

,ρ
1
​

; image and seed identifiers; all grids; reconstruction arms; TV-selection rule; quality targets; censoring rule; bootstrap procedure; figure subjects and crop coordinates; code, environment and input SHA256; all shard manifests; and the complete expected-cell table.

After S2 begins, changes are limited to demonstrable implementation defects established by outcome-blind unit tests. Every such change requires a version bump and rerunning all affected cells. No method, hyperparameter rule, grid point, image, seed, or failed/slow cell may be added, removed, or replaced based on observed reconstruction quality.

再补四条：

允许 S1 调整的白名单

数值稳定；

solver tolerance；

memory/block size；

manifest/shard 大小；

经无真值测试确定的路径覆盖范围；

明确的软件 bug。

禁止 S1 调整的黑名单

主 claim；

主高通量点；

success gate；

基线增删；

最终图像/crop；

根据 PSNR 选择 TV 规则；

根据初步结果删掉“不好看”的 mismatch 档。

预算截断顺序提前冻结

先砍 jitter；

再砍次级 PnP；

再砍 128² stress；

永不砍主 S2-A、主失配锚点和 S4。

不允许跑到一半后按结果决定保留哪些格子。

SHA 必须真的是硬门

当前 shard runner 对缺 SHA 的输入只给 warning，并继续运行。
冻结阶段必须改成：

缺 SHA 即失败；

代码 commit、环境 lockfile、图像、manifest 和分析脚本都进入 provenance；

合并前检查 expected cell、重复 cell、缺 cell。

还有两个 S0 级测试问题

当前 unit test 把“QMLE 必须比 Poisson 高 1 dB”写成 PASS 条件。方法是否优于基线是科学结果，不是软件单元测试。

Hadamard 单测手工构造 meta["pairs"]，但真实 patterns.py 输出的是 meta["pair_indices"]；所以测试避开了实际 integration bug。

应改成：

软件 UT 只检验分布、梯度、objective、配对、accounting；

性能差异一律是 smoke/report，不作为代码 PASS 门；

增加真实 make_patterns → campaign → run_arm integration test。

1.7 §6：四联主图

四联结构可以保留，但需四处改写。

Panel a

不要只画 CLT 单调曲线。应画：

exact count-FI 等高线；

principal ridge ρ
⋆
(ν)；

10% information-discrepancy boundary；

deployment zone；

CLT 曲线用虚线标为 approximation。

有限 ν 下，主峰之后可能有离散计数产生的次级小波峰；应称 principal/global information ridge 和“decaying envelope”，不要声称峰后严格单调。

Panel b

图和 crop 必须冻结。使用：

一张预注册自然图；

一个文字/细线靶；

radiometric 显示范围统一；

不逐列独立自动拉伸。

Panel c

横轴是

T
opt
​

=M
physical
​

T,

不是 T 或 M 单独。

纵轴主用 radiometric PSNR。

Panel d

应画 empirical speedup 与 exact-FI time-efficiency envelope。CLT

ρ/(1+ρ)

只能在其信任区内作为理论近似；过了 ridge 不能继续作为“上界”。

2. 最优负载脊线的解析推导

下面把 τ 作为时间单位，令

τ=1,T=ν,ρ=λτ,θ=logρ.
2.1 exact count likelihood 与解析 Fisher

active-start 下，第 m 个记录时刻为

S
m
​

=(m−1)+Γ
m
​

,Γ
m
​

∼Gamma(m,rate=ρ).

令

z
m
​

=ρ(ν−m+1).

则

G
m
​

(ρ):=P(N≥m)=P(Γ
m
​

≤ν−m+1)=P{Pois(z
m
​

)≥m}.

计数 pmf 为

p
m
​

=G
m
​

−G
m+1
​

.

关于 θ=logρ 的导数可以完全解析：

G
˙
m
​

:=
∂θ
∂G
m
​

​

=ρ
∂ρ
∂G
m
​

​

=
(m−1)!
e
−z
m
​

z
m
m
​

​

,m≥1.

因此

p
˙
​

m
​

=
G
˙
m
​

−
G
˙
m+1
​

,

精确的 count-only Fisher information rate 是

J
ν
​

(ρ)=
ν
1
​

m
∑
​

G
m
​

−G
m+1
​

(
G
˙
m
​

−
G
˙
m+1
​

)
2
​

​

不需要对 ρ 做中心差分。

当前 fisher_exact 使用有限差分，fisher_ridge.py 只取 49 个对数网格点，再做三点二次插值，因此 JSON 数字只能算 preliminary。

2.2 关键机制：count-only 丢失终端 detector phase

定义窗口结束时剩余的 dead-time phase：

R
ν
​

={
S
N
​

+1−ν,
0,
​

窗口结束时 detector 仍 dead,
窗口结束时已 active.
​

窗口内总 active 时间为

L
ν
​

=ν−N+R
ν
​

.

如果保留完整记录时刻，关于 θ 的 score 是

U
full
​

=N−ρL
ν
​

.

完整记录的 Fisher 信息是

I
full
​

=E[N].

只观察总计数 N 时，Fisher identity 给出

U
N
​

=E[U
full
​

∣N].

利用 missing-information identity：

I
N
​

=E[N]−ρ
2
E[Var(R
ν
​

∣N)].
​

这是本问题最干净的物理解释：

CLT 把计数看成连续量，近似继承了完整事件记录的信息；精确的 count-only 观测却不知道窗口结束时 detector 还剩多少 dead time。这个终端相位不确定性被乘上 ρ
2
，最终压倒高通量的边际收益。

这也和最新的 dead-time event-detection 理论方向一致：Jorgensen–Johnson 强调，较丰富的 activation/event 统计包含传统 histogram 或压缩计数丢掉的信息，但他们研究的是周期门控事件序列，而不是这里的单窗口 integrated count ridge。
arXiv

2.3 终端相位矩与一阶修正

长窗口下，随机观察时刻的 R 有：

P(R=0)=
1+ρ
1
​

,

以及在 0<r<1 上的常密度

f
R
​

(r)=
1+ρ
ρ
​

.

因此

E[R]=
2(1+ρ)
ρ
​

,
VarR=
12(1+ρ)
2
ρ(ρ+4)
​

.

由

E[N]=ρE[L]=ρ(ν−E[N]+E[R])

得到

E[N]=
1+ρ
νρ
​

+
2(1+ρ)
2
ρ
2
​

+o(1).

如果暂时忽略 N 对终端相位的少量预测能力，即用

E[Var(R∣N)]≃VarR,

则

J
ν
​

(ρ)=
1+ρ
ρ
​

+
12ν(1+ρ)
2
ρ
2
(6−4ρ−ρ
2
)
​

+⋯.

在 1≪ρ≪
ν
​

下，

J
ν
​

(ρ)=1−
ρ
1
​

−
12ν
ρ
2
​

+O(
ρ
2
1
​

+
ν
ρ
​

).
​

这里：

1/ρ：离完整事件信息饱和值还剩的缺口；

ρ
2
/(12ν)：count-only 丢失终端 phase 的代价。

2.4 立方根标度与领先常数

令

ρ=cν
1/3
.

则

J
ν
​

(ρ)=1−(
c
1
​

+
12
c
2
​

)ν
−1/3
+o(ν
−1/3
).

最大化信息等价于最小化括号：

−
c
2
1
​

+
6
c
​

=0,

所以

c
3
=6.

因此

ρ
⋆
(ν)∼6
1/3
ν
1/3
=1.8171206ν
1/3
.
​

峰值信息为

J
ν
​

(ρ
⋆
)=1−
26
1/3
3
​

ν
−1/3
+o(ν
−1/3
).
​

所以你们经验拟合的 1.7 是有限 ν 下的有效系数，不是渐近常数。

2.5 二阶修正：为什么数值更接近 −2/3

终端 phase 并非完全不能由 N 预测。令

C
ρ
​

=
ν→∞
lim
​

Cov(N
ν
​

,R
ν
​

).

利用 score identity

∂
θ
​

E[R]=Cov(R,U
full
​

)=(1+ρ)C
ρ
​

−ρVarR,

可得

C
ρ
​

=
12(1+ρ)
3
ρ(ρ
2
+4ρ+6)
​

⟶
12
1
​

.

同时

VarN=
(1+ρ)
3
νρ
​

+O(1).

bivariate renewal local-CLT 的首个回归项给出

Var{E[R∣N]}≃
VarN
C
ρ
2
​

​

.

因此 information expansion 多出正修正

ν
ρ
2
​

VarN
C
ρ
2
​

​

.

在 ridge 标度下，形式二阶展开为

J
ν
​

(ρ)=1−
ρ
1
​

−
12ν
ρ
2
​

+
ρ
2
1
​

−
6ν
ρ
​

+
144ν
2
ρ
4
​

+o(ν
−2/3
).

令

ρ=aν
1/3
+b+O(ν
−1/3
),

逐阶解 ∂
ρ
​

J=0：

a
3
=6,b=−
3
2
​

.

所以更精细的预测是

ρ
⋆
(ν)=(6ν)
1/3
−
3
2
​

+O(ν
−1/3
).
​

这与当前数值非常吻合：

ν repo ρ
⋆
(6ν)
1/3
(6ν)
1/3
−2/3
20 4.53 4.93 4.27
50 6.16 6.69 6.03
100 7.87 8.43 7.77
200 9.99 10.63 9.96
500 13.77 14.42 13.76
1000 17.45 18.17 17.50
2000 22.16 22.89 22.23

repo 数值来自粗 ρ 网格与有限差分。

我的严格表述建议是：

6
1/3
ν
1/3
：可以作为主要渐近命题；

−2/3：作为 formal second-order proposition；

在论文定稿前补一个 uniform local-CLT/Edgeworth 证明，或诚实写成“second-order asymptotic supported by exact numerical evaluation”。

2.6 “方差边界”应如何改写

在 ridge：

VarN≃
(ρ
⋆
)
2
ν
​

=
6
2/3
ν
1/3
​

→∞.

所以信息峰不是 Gaussian pmf 宽度降到 1。

真正的判据是边际平衡：

继续增大通量的收益
dρ
d
​

(−
ρ
1
​

)
​

​

=
ρ
2
1
​

,
终端 phase 损失
dρ
d
​

(−
12ν
ρ
2
​

)
​

​

=−
6ν
ρ
​

.

令两者相等：

ρ
3
=6ν.

“计数被钉死”发生在更晚的

ρ≍
ν
​

区间。

2.7 理论先例与 novelty

最近邻按危险程度排序：

Müller 1973 已系统处理 non-extending dead time 的 renewal count law，包括有限观察窗与起点约定；因此 exact count pmf 本身不是新贡献。
ScienceDirect

Yu–Fessler 2000 已给出死时间计数的一、二阶矩以及强通量下的强度估计器；所以均值、方差与反演也不是新贡献。
deepblue.lib.umich.edu

Alvarez 2014 使用 renewal CLT 的均值、方差和 Gaussian/CRLB 分析 pileup 性能；这正是你们 CLT proxy 的理论近邻。
aapm.onlinelibrary.wiley.com

Grönberg–Danielsson–Sjölin 2018 是最危险的理论先例：他们推导了更一般的 nonzero pulse-length count distribution，并使用 Fisher information 评价成像性能。冻结 novelty 表述前必须逐式核对全文，尤其是其高通量 Fisher 图。当前可检索摘要没有给出这里的 active-start、count-only principal ridge 或 ν
1/3
标度，但不能仅凭摘要断言不存在。
aapm.onlinelibrary.wiley.com

Rapp–Ma–Dawson–Goyal 2019/2021 已明确证明：准确建模 dead-time detection-time 序列可在高通量下大幅减少 lidar acquisition time。他们保留时间分布和 Markov sequence，而不是每个 pattern 只保留一个 integrated count；因此不会得到同一个 count-only ridge，但“高通量可比传统 5% rule 更快”绝对不能由你们首创。
arXiv
+1

Jorgensen–Johnson 2026 已给出更一般 dead-time event-detection 的充分统计、Fisher rate 和渐近有效估计。他们还强调 activation statistics 中存在 histogram 丢掉的信息。你们的终端 phase identity 是这个思想在 constant-flux integrated-count 场景下的更具体结果。
arXiv

novelty 判断

目前我没有检索到文献明确给出：

ρ
⋆
(ν)∼(6ν)
1/3

这一 active-start、non-paralyzable、integrated-count-only Fisher ridge。

但论文应写：

We derive a finite-window count-information ridge…

不要写：

We are the first…

论文价值

这个理论结果非常值得进入 A：

它解释为什么“更高通量并非无限更好”；

它把 QMLE 的部署区与 exact 区分开；

它给四联主图的 Panel a 一个真正的理论骨架；

它使文章不再只是“换 likelihood 多几 dB”。

暂不建议独立短文。要拆稿，至少还需要：

对 ρ=O(ν
1/3
) 的统一渐近证明；

−2/3 二阶项的严格证明；

非整数 ν、stationary-start/delayed-start；

有限样本上下界；

最好再有 paralyzable 或随机 dead-time 扩展。

2.8 QMLE 信任边界的正确表述

当前 JSON 的 exact/CLT=0.9 边界为：

ν provisional ρ
0.9
​

20 4.97
50 8.59
100 12.37
200 17.83
500 25.69
1000 37.01
2000 53.32

大 ρ 下，

J
CLT
​

J
exact
​

​

≃1−
12ν
ρ
2
​

.

令损失为 10% 得

ρ
0.9
​

(ν)≍
1.2ν
​

.
​

因此有两个不同标度：

最优 count-information ridge：

ρ
⋆
≍ν
1/3
;

CLT 10% 信任边界：

ρ
0.9
​

≍ν
1/2
.

部署区 [0.3,1] 确实远低于所有当前边界。这是好消息。

但审稿表述不要说：

exact is the low-count mode.

正确是：

exact is the short-window / extreme-saturation reference mode.

建议划分：

区域 定义 表述
RQL deployment zone ρ
95
​

≤1, ν≥20 主生产算法；exact FI 与 quasi-score 信息局部一致
transition zone 1<ρ
95
​

<ρ
0.9
​

(ν) RQL 可用，但需 exact 诊断
exact-reference zone ρ
95
​

≥ρ
0.9
​

(ν) Gaussian/CLT moment approximation不再适合定量信息预测
information-decreasing zone ρ≥ρ
⋆
(ν) principal count-information envelope下降；换 exact likelihood 也不能恢复已被 count compression 丢失的信息

这里应使用逐 pattern 的 ρ
95
​

，不能只用均值
ρ
ˉ
​

。

另外，FI ratio 只是标量局部信息诊断，不能单独证明高维 RQL 重建与 exact 重建等价；S4 的图像级验证仍然必须保留。

3. 方法学裁决
3.1 Full-Gaussian log-det 病灶
(i) Wedderburn quasi-score 是否正确？

是，方向完全正确。

Wedderburn quasi-likelihood 只要求指定条件均值和方差，不要求存在相应的完整概率密度。
DOI

对于

μ(λ)=
1+λτ
λT
​

,V(λ)=
(1+λτ)
3
λT
​

,

quasi-score 是

U
λ
​

=
V(λ)
μ
′
(λ)
​

[N−μ(λ)].

因为

V
μ
′
​

=
λ
1+λτ
​

,

可化简为

U
λ
​

=
λ
N
​

−(T−Nτ).
​

其积分型负 quasi-likelihood 为

Q(λ;N)=(T−Nτ)λ−Nlogλ+const.
​

因此对于主模型 σ
b
​

=0，你们甚至不需要四轮 IRLS。可以直接解：

x≥0
min
​

M
1
​

i
∑
​

[(T−N
i
​

τ)λ
i
​

(x)−N
i
​

logλ
i
​

(x)]+αTV(x),
​

λ
i
​

(x)=Φa
i
⊤
​

x+d.

只要 N
i
​

τ<T，这是一个干净的凸 data fidelity。

这带来一个必须面对的 novelty 事实

单帧无正则根为

λ
^
i
​

=
T−N
i
​

τ
N
i
​

​

,

恰好就是 non-paralyzable pre-correction。

所以：

RQL 与 PRECORRECT 不是两个互不相干的思想；

PRECORRECT 是先逐帧反演，再做二次近似；

RQL 是把同一 quasi-score 直接嵌入联合图像逆问题，避免逐帧爆炸、clipping 和误差传播。

这进一步说明，A 的 novelty 不能落在“发明一个新 estimator 公式”，必须落在：

高通量 operating-point map；

time-to-quality；

count-Fisher ridge；

mismatch robustness；

联合逆问题比逐帧校正更稳定。

当前代码的 IRLS 仍有一个具体不一致：TV 选择时权重冻结在 x
0
​

，最终拟合却迭代更新权重。

直接使用上面的凸 quasi-likelihood 可以同时消掉：

selector/final objective 不一致；

四轮是否收敛的问题；

“Gaussian”命名；

frozen-weight 路径依赖。

只有 σ
b
​

>0 或更复杂 afterpulse moment model 时，IRLS 才作为通用实现保留。

饱和边界

若 N
i
​

τ=T，则线性项消失：

Q=−N
i
​

logλ
i
​

,

随 λ
i
​

→∞ 无下界。

这不是实现 bug，而是 moment quasi-likelihood 在计数顶格处失去辨识性的真实边界。应：

在主 deployment zone 避免此概率；

exact likelihood 负责该区；

或加物理强度上界；

把出现 ceiling count 的比例作为诊断落盘。

(ii) “log-det 在 ρ=1/2 翻转”是否值得写？

值得，但只能写成一个严谨的小节或 Supplement，不应成为主贡献。

你们算得没错：

dlogλ
dlogV
​

=
1+ρ
1−2ρ
​

.

所以在

ρ=
2
1
​

处符号翻转。当前代码注释也已记录这一点。

但不能简单写：

log-det 在 ρ>1/2 必然造成偏差。

在一个真正正确的 Gaussian 模型下，

E[(N−μ)
2
]=V,

残差项对 V
′
的期望会与 log-det 项抵消，真参数处期望 score 仍为零。

真正可证明的病灶是：

heteroscedastic Gaussian pseudo-likelihood 的 variance-collapse degeneracy。

例如顶格计数 N=ν、τ=1 时，ρ→∞：

μ=
1+ρ
νρ
​

=ν−
ρ
ν
​

+o(ρ
−1
),
V=
(1+ρ)
3
νρ
​

=
ρ
2
ν
​

+o(ρ
−2
).

所以

2V
(N−μ)
2
​

→
2
ν
​

,

而

2
1
​

logV=−logρ+O(1).

于是 full-Gaussian NLL

⟶−∞.

这是严格的非 coercive 退化。欠定逆问题里，一些信息贫方向可以把若干近饱和预测继续抬亮，同时通过缩小 V 获得伪奖励；你们观测到的 0.9 dB 和径向误差 1.5 与这个机制一致。

建议小节标题：

Why a full heteroscedastic Gaussian likelihood fails at high load

内容控制在半页：

ρ=1/2 符号翻转；

ceiling-row 非 coercivity 推导；

一张 radial objective 曲线；

Full-Gaussian 失败图；

RQL quasi-score 修复图。

3.2 TV 强度的无真值选择
先给根本裁决

你们观察到的现象不是 selector 写得不够聪明，而是欠定逆问题的基本限制：

数据预测风险和图像重建风险不是一回事。数据只约束 measurement-visible 子空间；强 TV 对 null/near-null 空间的改善可以提升 PSNR，却同时轻微恶化 held-out data fit。

因此没有任何纯数据规则能保证复现 truth-oracle PSNR 最优 λ
TV
​

。这和你们自己的 no-free-audit 刻画是一致的。

目标应改成：

选择与已标定测量统计相容的最强正则化，而不是宣称无真值地找到 PSNR-oracle 正则化。

为什么当前 own-held-out-NLL 必须废弃

当前选择器：

前 90% 拟合；

最后 10% 用该臂自己的 data term评分；

每臂又用自身 x
0
​

梯度范数缩放 λ 网格。

这造成两个偏差：

正确模型和错模型的 NLL 不在同一标尺上；

错模型的巨大梯度会自动获得更强 TV，正是你们看到的 POISSON-LIN 偶然占优。

“每臂都用同一程序”不等于公平，因为各自 objective 的单位、曲率和模型偏差不同。

SURE/GSURE 是否适合作为主规则？

不建议。

Eldar 的 generalized SURE 主要针对 exponential-family 模型；你们的 exact renewal pmf 是 Gamma-CDF difference，并不是可以直接套用的简单一参数指数族。要做 exact renewal Stein/SURE，几乎会再开一篇方法论文。
arXiv

projected GSURE 对 rectangular/ill-posed operator 估计的是可观测投影误差，不能评估真正不可见的 null-space 结构。
arXiv

对严重病态问题，SURE/GSURE 本身可能系统性偏向过小正则化；已有理论和数值研究明确报告这种现象。
arXiv

可以把 GSURE 作为附录 sensitivity，不应作为冻结主 selector。

3.3 推荐主规则：cross-fitted common renewal discrepancy

建议预注册一个统一规则：

选择通过独立 renewal 残差检验的最大 TV 强度。

Step 1：固定 folds

K=5；

用冻结 hash 分 fold，不用“最后 10%”；

Hadamard complementary pair 必须进入同一 fold；

continuous acquisition 必须按完整连续 block 分 fold，block 之间用预注册 reset/guard，不能逐帧随机打乱。

Step 2：使用无量纲正则化路径

不要共用数值 λ
TV
​

。

建议定义

η=
λ
max,arm
​

λ
TV
​

​

,

其中 λ
max,arm
​

是该臂达到 TV-null/近常数解所需的最小 penalty，通过冻结的 bisection 规则求得。

例如冻结：

η∈{10
−4
,3×10
−4
,10
−3
,…,1}.

更严格的替代是用相同的 achieved roughness：

R(x)=
∥x∥
1
​

TV(x)
​

.

数值 λ 可以不同，但先验强度路径必须可比。

Step 3：各臂保留自己的 reconstruction data fidelity

POISSON-LIN 用自己的 Poisson data term；

SAT 用自己的 mean-only data term；

PRECORRECT 用自己的 WLS；

RQL 用 renewal quasi-likelihood。

这保持方法定义不被改变。

Step 4：所有臂用同一个外部物理 evaluator

在 held-out fold 上统一计算：

r
i
​

(η)=
V
NP
​

(
λ
^
i
​

)
​

N
i
​

−μ
NP
​

(
λ
^
i
​

)
​

,
D(η)=
i∈V
∑
​

r
i
​

(η)
2
.

这里的 μ
NP
​

,V
NP
​

来自独立 detector calibration，不来自图像真值。

阈值不要机械使用 χ
m
2
​

，而应由冻结的 exact renewal parametric bootstrap 标定：

条件在 (ρ,ν) 和 held-out patterns；

生成 exact renewal counts；

得到 D 的 2.5%–97.5% 接受区；

同时标定 residual mean 和 residual–predicted-load correlation。

Step 5：one-SE + 最大相容正则化

先找 common held-out renewal deviance 最小的 η
min
​

。

定义候选集合：

E={η:D
dev
​

(η)≤D
dev
​

(η
min
​

)+SE
min
​

,D(η) 通过 GOF}.

最终选：

η
⋆
=maxE
​

即在数据无法拒绝的范围内选择最平滑解。

若 E=∅：

选择 common deviance 最小的值作为描述性结果；

标记 MODEL_FAIL；

不允许转而看 PSNR 调参数。

Step 6：用全部帧重拟合

折外选择完成后，用冻结的 η
⋆
在全数据上重拟合。

这本质上是带独立残差检验的 Morozov-style discrepancy principle。非线性逆问题中使用 discrepancy principle 选择正则化参数有成熟理论基础。
OSTI.gov

3.4 这套规则如何解决 POISSON-LIN 的“意外重 TV”问题

POISSON-LIN 可以继续因自身模型梯度得到某条重建路径，但它不能用自己的错 Poisson NLL宣布该路径“预测最好”。

所有臂都必须接受同一 exact-renewal calibrated residual test：

如果重 TV 恰好消除了某些结构，但仍与真实 detector statistics 相容，它可以被选择；

如果重 TV 只是在错模型下隐藏偏差，common residual 会显示：

残差均值偏移；

方差不匹配；

残差与预测 ρ
i
​

相关；

GOF 拒绝。

给错模型使用正确 detector calibration 做选择并不偏袒提出法，反而是在强化基线：所有方法都获得同一份可部署校准信息。

3.5 可直接放进论文的方法公平性表述

All reconstruction arms retained their declared data-fidelity models. Hyperparameters were selected without image ground truth using identical cross-validation folds, the same isotropic-TV functional, the same normalized regularization path, and a common detector-calibrated renewal discrepancy evaluator. Numerical TV coefficients were not forced to be equal because the data fidelities have different units and curvature. Among statistically admissible candidates, the strongest regularization was selected by a preregistered one-standard-error discrepancy rule.

还应附上：

truth-oracle λ 只作诊断，标记 ORACLE—NOT DEPLOYABLE；


[output truncated at 50000 of 51458 characters. Pass a larger max_chars (default 50000) to see more, or use read_page with a ref_id to focus on a smaller section.]