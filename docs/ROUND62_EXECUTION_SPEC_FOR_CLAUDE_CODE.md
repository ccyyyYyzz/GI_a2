# ROUND62 执行规范：MAVE 异常审计（T1-G0）× 头部空间-规模地图（T1-G1）× 块长交叉地图（P3）

**日期：** 2026-07-17
**写入者：** Fable5（依据 `ROUND62_TWO_TRACK_PLAN.md`）
**执行者：** Claude Code（GI_a2 会话，`D:\GI_another`，续用既有 harness）
**上游：** ROUND59 规范与理论手册仍然有效（照明族、链路、指标协议、种子体系全部沿用）；本文件只定义增量。

---

## 0. 硬性守则（增量）

1. **复用不重写**：沿用 `code/` 既有照明族、DT30 renewal、估计器、指标协议、`SEED0=20260717` 体系；新代码入 `code/round62/`，新结果入 `results/round62_g0/`、`round62_g1/`、`round62_p3/`。不改动、不覆盖 Phase A/B 任何已有结果。
2. **封闭盆地纪律**：任何 score/密度参考臂**禁止出现**（ROUND59 终局封闭）。
3. **三个 Part 相互独立**：G0 失败 ⇒ 跳过 G1，但 **P3 无论如何都跑**（Track 2 保底）。
4. **预登记不可改；结果原样报告；每 Part 完成即 commit。**
5. **预算**：Part 1 ≤ 3 CPU·h；Part 2 ≤ 10 CPU·h；Part 3 ≤ 4 CPU·h。超限中断报告。
6. **MAVE 全程禁止接触真值**；一切超参选择只允许用留出帧 MSE。所有 MAVE/OPG 超参（带宽、ridge、迭代数、子采样）逐项落盘。

## 1. Part 1 —— T1-G0：MAVE 异常审计（16×16，GAM4×DT30 为主，3 seeds）

审计对象：Phase A 中 MAVE-16 = 30.19 dB（超最强白化 +6.9 dB）。四项检验：

**(a) 置换零检验（硬门）。** 对每 seed 生成 3 个固定置换 $\pi$，在 $(a_i, b_{\pi(i)})$ 上重跑同一 MAVE 管线；另跑一个纯噪声臂（$b$ 换为方差匹配的独立噪声）。门：
- 未置换 − 置换 ≥ **15 dB**（8 图平均），且置换后 ≤ 12 dB；纯噪声臂同标准。
- 不过 ⇒ MAVE 数字含隐式平滑/正则化红利伪装的重建，`ARTIFACT`，T1 终止。

**(b) M-scaling 诊断（软门，只报告）。** M ∈ {5e3, 1e4, 2e4}：MAVE 相对 WHITEN-LW 的增益曲线。增益随 M **减小而增大** ⇒ 平滑先验签名，标记可疑并入报告。

**(c) MLE-OR 重做（硬门的对照组成部分）。** Phase B 的 MLE-OR 用了 Poisson 似然（对 DT30 的 renewal 噪声错配）且有溢出警告。重做：
- 主实现：高斯-renewal 近似似然，均值 $m(u)=\frac{\lambda T}{1+\lambda\tau}$、方差 $v(u)=\frac{\lambda T}{(1+\lambda\tau)^3}$，$\lambda=su$（Fano $=(1+\lambda\tau)^{-2}\approx0.49$ 于工作点）；
- 精确校验（子样 1000 帧）：$P(N\ge m)=\Gamma\text{CDF}(T-m\tau;\,\text{shape}=m,\,\text{rate}=\lambda)$（与 harness 的"全部间隔含 $\tau$"约定一致），$P(N=m)=P(N\ge m)-P(N\ge m+1)$，比对高斯近似的似然面一致性；
- 多起点（WHITEN-LW / L-ISOTRON / MAVE 输出）+ 梯度范数收敛诊断，报告 best-of。

**(d) 敏感性与增益来源分类（诊断）。** 带宽 ×{0.5, 1, 2}、ridge ×{0.1, 1, 10} 的 MAVE 稳定性；另跑 **GAM4×LIN**（16×16）：若 MAVE 在线性链路仍保有 ≥50% 增益 ⇒ `GENERAL_EFFICIENCY`（噪声加权/自适应故事）；若增益基本消失 ⇒ `NONLINEARITY_SPECIFIC`。两者都可发表，故事不同。

**G0 总门：** (a) 通过，且审计后（含重做 MLE、超参敏感性最差档）MAVE 相对 max(WHITEN-LW, L-ISOTRON, MLE-renewal) 仍 ≥ **+3.0 dB**（6/8 图）⇒ PASS，进 G1。否则 `T1_KILL`。

## 2. Part 2 —— T1-G1：头部空间-规模地图（G0 过后执行）

**半参数 oracle 类（预登记定义，防分叉花园）：** {OPG, rMAVE}，带宽 = rule-of-thumb ×{0.5,1,2}，ridge=1e-4，迭代封顶 20；每配置按**留出帧 MSE** 选最优成员与超参（禁真值选择）。$n=1024$ 及以上允许子采样 OPG（锚点 ≤ 5000，平均外积），近似方式落盘。

**网格：** $n\in\{256, 576, 1024\}$（4096 为可选扩展，预算内才跑，用子采样 OPG）× $M\in\{2\times10^4, 10^5\}$ × 照明 {GAM4, GAM8, CORR-LOGN} × 链路 {DT30, SAT30} × 3 seeds × 8 图。基线：WHITEN-LW、L-ISOTRON（协议与 Phase B 相同）、WHITEN-OR（解析族）。指标协议不变（通量匹配主协议）。

**G1 门：** $n=1024, M=10^5$ 下，{GAM4×DT30, GAM8×DT30} 至少其一：半参数类 − max(WHITEN-LW, L-ISOTRON) ≥ **+1.0 dB** 且 6/8 图。通过 ⇒ `T1_PROCEED_TO_G2`（G2 规模化构造器另立规范）。不过 ⇒ `T1_KILL`，输出头部空间随 $n$ 与 $M/n$ 的衰减曲线（併入 Track 2 论文消融节）。

**必交图：** headroom vs $n$（各照明×链路一条线）；headroom vs $M/n$ 合并散点；CORR-LOGN 行如实呈现（预期为负——regime 边界素材）。

## 3. Part 3 —— P3：块长交叉地图（Track 2，无条件执行；64×64）

**漂移模型（新增 nuisance）：** 帧间乘性增益 $g_t$，log-AR(1)：$\ln g_t=\rho\ln g_{t-1}+\sigma_\eta\varepsilon_t$，$\sigma_{\ln}=\sqrt{\ln(1+\mathrm{CV}^2)}$，$\sigma_\eta=\sigma_{\ln}\sqrt{1-\rho^2}$，均值校正使 $\mathbb E[g]=1$；$\rho=0.99$ 主档（相关时间 ~100 帧），$\rho=0.9$ 敏感档。bucket：$b_t=\mathrm{Poisson}(s\,g_t\,u_t)/s+\mathcal N(0,\sigma_r^2)$（**线性链路**，本实验只考漂移免疫）。

**估计器族（块中心化相关 = 块 multinomial 条件化的线性化）：**

$$\hat x_L=\frac1M\sum_{\text{blocks}}\sum_{i\in\text{blk}}(b_i-\bar b_{\text{blk}})(a_i-\bar a_{\text{blk}}),\qquad L\in\{2,4,8,16,32,64\},$$

外加两个端点臂：`RAW`（不中心化，失败展示）与 `GLOBAL`（全局中心化 = 标准 GI，$L=M$）。每臂各出两个读出：普通相关与 WHITEN-LW（$\hat\Sigma$ 用块中心化 $a$ 估计）。次级臂：比值归一 $b_i/\bar b_{\text{blk}}$（NDHSI 实践形态，报告但不评门）。秩代价 $M/L$ 自由度与漂移免疫的权衡即为地图。

**网格：** CV ∈ {0.02, 0.05, 0.1} × $\rho$ ∈ {0.99, 0.9} × 照明 {GAM4, MIX-LOGN} × $M=2\times10^4$ × $s=10^4$ × 3 seeds × 8 图。

**P3 门（ROUND61 预登记）：** 存在 $L\in[4,64]$ 在 ≥2 个 CV 档同时超过 `L=2` 与 `GLOBAL` 各 ≥ **+0.3 dB** 且 6/8 图 ⇒ `MIDDLE_L_DOMINATES`（交叉地图成为管线论文的设计贡献）。否则 `ENDPOINTS_ONLY`（管线第 1 级用 L=2，诚实引 NDHSI，照常推进）。**必交图：** PSNR–L 曲线族（按 CV×ρ 分面），标注交叉点。

## 4. 输出契约

```
results/round62_g0/  g0_audit.json, g0_metrics.csv, figures/(置换塌缩图, M-scaling, 带宽敏感性)
results/round62_g1/  g1_scale_map.csv, g1_gates.json, figures/(headroom-vs-n, headroom-vs-M/n)
results/round62_p3/  p3_crossover.csv, p3_gates.json, figures/(PSNR-L 曲线族)
REPORT_ROUND62.md    （结尾裁决块如下）
```

```text
ROUND62_T1_G0_AUDIT: PASS | ARTIFACT | KILL(...)
ROUND62_T1_G0_GAIN_SOURCE: NONLINEARITY_SPECIFIC | GENERAL_EFFICIENCY | ARTIFACT
ROUND62_T1_G1_SCALE_MAP: PASS(+X.XX dB @ n=1024) | KILL(衰减曲线要点) | NOT_RUN
ROUND62_P3_CROSSOVER: MIDDLE_L_DOMINATES(L=?, CV=?, ρ=?) | ENDPOINTS_ONLY
ROUND62_VERDICT: T1_PROCEED_TO_G2 | T1_DEAD_T2_PROCEED
ROUND62_RUNTIME: <各 Part CPU·h>
```

**给执行会话的一句话总纲：** 这一轮你在做两件事——用零检验和公平对照审判你自己上一轮跑出的 +6.9 dB 异常（它可能是发现也可能是伪影，两个答案都值钱），以及把 60 轮悬而未决的块长交叉问题一次跑清；门槛照旧不可改，坏消息照旧原样写。
