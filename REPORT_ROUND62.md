# ROUND62 REPORT — MAVE 异常审计（T1-G0）× 头部空间-规模地图（T1-G1）× 块长交叉地图（P3）

- **规范**: `docs/ROUND62_EXECUTION_SPEC_FOR_CLAUDE_CODE.md`（ROUND59 harness 全部沿用）
- **执行**: Claude Code (Fable 5)，2026-07-17，仓库 GI_a2；G1/图表脚本由 Opus 子代理编写并经主循环逐行核查
- **纪律**: 封闭盆地（全程无 score/密度参考臂）；MAVE 全程未接触真值；超参逐项落盘于 `results/round62_g0/g0_audit.json` 与 `results/round62_g1/`（未跑）；跑前 7-agent 对抗审计确认并修复 4 项（投影梯度诊断、P3/G0 预算中断工件标记、G0 预算守卫死代码）
- **算力**: 全部本机 CPU（i9-13900HX 24C）。Part 1 wall 417 s / process-CPU 141 s（预算 ≤3 CPU·h ✓ 双读法）；Part 3 见 §3

**一句话总结（在本仿真协议下）**：上一轮的 +6.9 dB MAVE 异常**不是伪影**（置换零检验分离 26.6 dB），但 (i) 它是**与非线性无关的一般效率优势**（LIN 增益 8.18 dB > DT30 增益 6.91 dB），(ii) 带**平滑先验签名**（增益随 M 增大缩水 8.23→6.91 dB），(iii) 有**超参脆弱轴**（ridge×10 摧毁至 12.94 dB），且 (iv) **被修正似然的参数 MLE 直接覆盖**（MLE-renewal 33.09 dB > MAVE 30.19 dB——Phase A 的 Poisson-MLE "上界" 22.08 dB 是噪声模型错配的假象）——按预登记总门判 **T1_KILL**（对读法稳健），G1 不跑；P3 交叉地图见 §3。

---

## 1. Part 1 —— T1-G0 MAVE 异常审计（16×16, GAM4, 3 seeds）

### (a) 置换零检验（硬门）：PASS（异常为真）

| 臂 | 8图×3seed 平均 PSNR | 门槛 |
|---|---|---|
| MAVE 未置换（= Phase A 流精确复现） | **30.19 dB** | — |
| MAVE 置换（3 个固定 π × 3 seeds 合并） | 3.56 dB（单臂最大 7.03） | ≤12 ✓ |
| MAVE 纯噪声（方差匹配独立噪声） | 3.82 dB | ≤12 ✓ |
| 分离 | 未置换−置换 = **26.63 dB**；未置换−噪声 = 26.37 dB | ≥15 ✓ |

破坏 (a_i, b_i) 配对后 MAVE 完全塌缩——30.19 dB 是真实的测量-图像信息重建，不含隐式平滑/正则化红利伪装。

### (b) M-scaling 诊断（软门）：平滑先验签名 = 触发

MAVE − WHITEN-LW 增益：M=5e3 → **8.23 dB**；1e4 → 7.61；2e4 → 6.91。增益随 M 单调缩水 ⇒ 有限样本效率优势成分显著（局部平滑在小样本更值钱），外推到大 M 的前景收窄——与 (d) 的结论互证。

### (c) MLE 重做（Gaussian-renewal 似然）：对照组就位，近似验证通过

- 主实现：均值 $u/(1+\tfrac37 u)$、方差 $u/[s(1+\tfrac37u)^3]+(2/s)^2$（Fano≈0.49 于工作点），L-BFGS-B、x≥0、三起点（WHITEN-LW / L-ISOTRON / MAVE）、按 NLL best-of。
- 精确校验（1000 帧）：$P(N\ge m)=\Gamma\mathrm{CDF}(T-m\tau;m,\lambda)$ 的精确 renewal 似然面与高斯近似沿 $\alpha x_{\rm true}$ 射线 **argmax 完全一致（均为 α=1.00）**，峰邻域偏差可忽略。
- 收敛诊断（投影梯度 inf-范数，审计修复后口径）：210 / 952 / 8749（三起点，seed 0）——未达严格一阶最优（maxiter 300 截断），如实报告。
- **审计关键发现：MLE-renewal = 33.09 dB（8图×3seed 均值），高出默认 MAVE (+2.90 dB)、高出白化 (+9.81 dB)。** Phase A/B 的 Poisson 版 MLE-OR（22.08 dB）确系噪声模型错配（renewal Fano≈0.49 vs Poisson Fano=1）所致的假上界；修正方差后，**正确设定的参数似然直接统治头部空间**——MAVE 的半参数优势真实但完全被正确似然覆盖，不存在残余的"MAVE 魔法"。

### (d) 敏感性与增益来源：GENERAL_EFFICIENCY，脆弱轴 = 过正则

带宽×{0.5,1,2} × ridge×{0.1,1,10}（seed 0 全网格，最差档再跑 seeds 1,2）：

| | r×0.1 | r×1 | r×10 |
|---|---|---|---|
| bw×0.5 | 29.95 | 30.18 | **12.94（最差）** |
| bw×1 | 30.17 | 30.16 | 21.10 |
| bw×2 | 30.25 | 30.19 | 28.87 |

7/9 配置稳定在 30±0.3 dB——带宽方向极稳；**唯一脆弱轴是 rMAVE H-矩阵的过正则（ridge×10）**。
增益来源：GAM4×LIN 下 MAVE−WHITEN-LW = **8.18 dB** ≥ 50%×DT30 增益（6.91 dB）⇒ **GENERAL_EFFICIENCY**——MAVE 的优势在线性链路上不缩反增，与非线性无关；这是"自适应局部平滑/噪声加权"的一般半参数效率故事，不支持任何非线性伴随叙事。

### G0 总门：KILL(margin<3dB)

预登记规则：(a) 通过 **且** 最差超参档 MAVE 相对 max(WHITEN-LW, L-ISOTRON, MLE-renewal) 仍 ≥+3.0 dB（≥6/8 图，3-seed 均值）。最差档 bw0.5_r10（12.94 dB）对最优基线逐图 Δ ∈ [−21.4, −17.3]，0/8 ⇒ **T1_KILL**。
**KILL 对门槛读法稳健**：即使按最宽的"默认配置"读法，MAVE (30.19) 对含 MLE-renewal (33.09) 的预登记对照集也是 **−2.90 dB**——远不及 +3.0 dB，0/8 图为正。异常真实、可复现、但 (i) 与非线性无关，(ii) 被正确似然覆盖，(iii) 有过正则脆弱轴：T1 按任何读法均不成立。

## 2. Part 2 —— T1-G1：NOT_RUN（G0 未过，规范 §0.3）

`code/round62/g1_scale_map.py` 已完成编写并通过审计（预登记 oracle 类 {OPG,rMAVE}×bw{0.5,1,2}、留出帧 MSE 选择、n∈{256,576,1024}×M∈{2e4,1e5} 网格、META 断点、预算守卫），未执行；如后续授权重跑可直接使用。

## 3. Part 3 —— P3 块长交叉地图（无条件执行；64×64, LIN + log-AR(1) 漂移）

网格完整：{GAM4, MIX-LOGN} × CV{0.02,0.05,0.1} × ρ{0.99,0.9} × 3 seeds × 8 图 × 臂{L2..L64, RAW, GLOBAL} × 读出{plain, whiten, ratio}；无预算中断。**裁决：`MIDDLE_L_DOMINATES`**（主档 ρ=0.99；预登记读法：逐读出族、逐照明族评门，任一成立即过——通过集合：GAM4-plain L∈{4,8,16}、GAM4-whiten L∈{4,8}、MIX-whiten L∈{4,8}；其中 **GAM4-plain L=4 在全部 3 个 CV 档通过**）。

**代表值（CV=0.05, ρ=0.99, 8图×3seed 均值 PSNR）：**

| 臂 | GAM4 plain / whiten | MIX-LOGN plain / whiten |
|---|---|---|
| GLOBAL | 4.86 / 4.86 | 17.07 / 3.67 |
| L2 | 8.00 / 8.62 | 16.98 / 5.84 |
| **L8** | **8.92 / 9.12** | **17.30** / 6.25 |
| RAW | 13.34†/ — | 13.35†/ 3.16 |

†RAW 的高 PSNR 是 **DC 基座假象**（未中心化的 b̄·ā 项给出近平坦图；通量匹配后 PSNR/SSIM 偏爱 DC 主导的自然图像）——RAW 仍是失败展示臂，不入任何比较。

**结构解读（在本仿真协议下）：**
1. **GAM4（iid 图样）= 教科书交叉**：漂移摧毁 GLOBAL（−3.1 dB vs L2），块中心化恢复，中段 L 支配（L=4~16 比 L2 高 +0.4~+1.5 dB、比 GLOBAL 高 +2.0~+4.1 dB）；最优 L* 随 CV 增大左移（CV=0.02 时 L*=16~32，CV=0.1 时 L*=4）——**秩代价 M/L vs 漂移免疫的权衡地图**正是预期形状，两个读出族一致。
2. **MIX-LOGN（相关图样）**：相关场自身的帧能量涨落远大于 CV≤0.1 的漂移 ⇒ GLOBAL 不塌（17.07）；块长的边际收益小（L8−L2 = +0.32 dB，对 GLOBAL +0.23 dB，仅 CV=0.1 档过门）。**whiten 读出对 MIX-LOGN 灾难性**（3.7–6.3 dB；LW 白化长相关对数正态场的病态放大，与 Stage-0 A3 的 CORR-LOGN 结论同源）——MIX-whiten 的门内通过是低绝对水平族内的比较，**部署建议以 plain/ratio 读出为准**。
3. **比值臂（NDHSI 实践形态，报告不评门）**：与 plain 数值几乎重合（GAM4 L8: 8.92 vs 8.92；MIX L8: 17.30 vs 17.30）——小漂移下除以块均值 ≈ 减去块均值；其价值在更大 CV/更强乘性扰动区间，本网格内无额外增益。
4. ρ=0.9 敏感档：同构但整体更低（去相关更快 → 有效漂移带宽升高），L* 进一步左移（多为 L*=4）；不改变裁决。

图：`results/round62_p3/figures/p3_psnr_vs_L.png`（分面曲线族+argmax 标注）；G0 图三张于 `results/round62_g0/figures/`。

**管线论文落点**：第 1 级采用块中心化 + **L∈[4,8]**（iid/弱相关照明），相关照明场景 plain/ratio 读出 + 块长不敏感；诚实引 NDHSI 作为比值归一先例。

## 4. 运行时间与算力

| Part | wall | process-CPU | 预算 |
|---|---|---|---|
| P1 (G0) | 417 s | 141 s 主进程（MAVE 在 8 个 joblib 子进程，合计约 0.5 CPU·h） | ≤3 CPU·h ✓ |
| P2 (G1) | 未跑（G0 KILL） | — | — |
| P3 | 1533 s | 12660 s = 3.52 CPU·h | ≤4 CPU·h ✓（wall 读法 0.43 h） |

全部本机 CPU（i9-13900HX），未动用 Colab/Drive（本机 24 核对此工作负载优于 Colab CPU/GPU 路径）。

## 5. 裁决块

```text
ROUND62_T1_G0_AUDIT: KILL(margin<3dB: 最差档 bw0.5_r10 0/8; 宽读法下默认档 −2.90 dB vs MLE-renewal; 置换零检验本身 PASS——异常为真)
ROUND62_T1_G0_GAIN_SOURCE: GENERAL_EFFICIENCY (LIN 增益 8.18 dB ≥ DT30 增益 6.91 dB; M-scaling 平滑先验签名触发)
ROUND62_T1_G1_SCALE_MAP: NOT_RUN (G0 未过, 规范 §0.3; 代码已就绪并通过审计)
ROUND62_P3_CROSSOVER: MIDDLE_L_DOMINATES(L=4–8; GAM4-plain L4 全 3 CV 档, ρ=0.99; GAM4 双读出 + MIX-whiten 同过; MIX 部署建议 plain/ratio 读出)
ROUND62_VERDICT: T1_DEAD_T2_PROCEED
ROUND62_RUNTIME: P1 wall 0.12 h (主进程 CPU 0.04 h + joblib 子进程 ~0.5 CPU·h); P3 wall 0.43 h / process-CPU 3.52 h; P2 未跑
```
