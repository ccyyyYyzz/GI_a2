# R25 ruling (GitHub issue #17, raw)

Title: R25 ruling: standardized maximin bank allocation and mixable-design requirement
Posted: 2026-07-22T01:48:49Z

---

# R25 — Final ruling: replace threshold routing by standardized maximin bank allocation

**Reference commit:** `d1b1c21`  
**Source question:** `docs/ROUND63_GPT_ROUND25_QUESTION.md`

**Scope:** pre-implementation amendment to the R24 DEV bridge only. T-A scenes remain frozen. T-C has not been built. No M1 result, paper, or campaign endpoint is changed.

> **R25 VERDICT: adopt the convex-allocation architecture, but reject the proposal in its literal form.**
>
> Retire the three runtime thresholds, hysteresis, and categorical `KNOB/LIBRARY/ABSTAIN` router. Replace them with one **standardized maximin estimator-risk allocation** over the library simplex, followed by a parameter-free tie-break toward the knob corner.
>
> Two amendments are non-negotiable:
>
> 1. optimize **worst relative regret**, not the raw worst utility; and
> 2. rebuild the library as **mixable approximate-design measures / guard-aware row-atom pools**. Arbitrarily taking `round(972 w_k)` rows from a jointly optimized 972-row bank is not valid and can destroy both information value and dose balance.
>
> The resulting method is a robust library **mixture**, not an abstaining router. The only hard fallback retained is operational: solver failure, exact-materialization failure, or a post-materialization guard violation returns the frozen knob bank and emits an explicit failure flag.

I recommend renaming the prototype from **Residual-Gated Library Illumination** to **Robust Library Mixture Illumination (RLMI)**. If the acronym `RGLI` must be retained for continuity, expand `G` as “guided,” not “gated.” The R23 rank-one residual remains a disclosure and diagnostic; it is no longer a runtime switch.

---

## Executive ruling

| claim | ruling |
|---|---|
| **P-A — convexity** | **Substantially true.** `log det` and `-tr(WX^{-1})` are concave in the simplex weights when `W\succeq0`, `X\succ0`, and every scenario matrix is affine in `w`. The pointwise minimum preserves concavity. “Milliseconds” is unproven: cost is dominated by `S` factorizations of the task matrix per iteration, not by `K=8`. |
| **P-B — continuity / emergent abstention** | **False as stated.** The optimal value is continuous, but the optimizer can be nonunique and discontinuous when the active worst scenario or support changes. Minimax mixing is diversification, not abstention, and it need not move toward `L0`. Use a lexicographic closest-to-`e0` tie-break; do not claim global continuity. |
| **P-C — guard inheritance** | **True only for continuous design measures, false for naive row subsets.** Linear guards and convex information guards inherit under a true convex mixture. A fraction of a whole exact bank is generally not that mixture. Exact constrained materialization and postcheck are mandatory. |
| **P-D — certificate** | **True.** The finite-scenario simplex problem has an exact KKT / least-favorable-scenario certificate. It certifies the declared local estimator-risk problem, not PSNR or the true scene. Recompute the realized regret after row materialization. |
| **P-E — deployment / bridge** | **Partly true after redesign.** Gates A, B, and D stay numerically unchanged. Gate C is replaced by a continuous allocation-informativeness gate below. The composite row sequence can be preloaded, but T-B must expose mixable row measures rather than monolithic banks. |

---

# Q1 — Verification and failure analysis

## 1. Convexity of the two proposed criteria

For posterior/bootstrap draw `s`, define

\[
X_s(w)=\mathcal H_{0,s}+\sum_{k=0}^{K-1}w_kF_{k,s},
\qquad w\in\Delta_K,
\]

where

- `F_{k,s}\succeq0` is fixed once the draw, bank rows, and source schedule are fixed;
- `\mathcal H_{0,s}\succ0` includes the pre-scan and frozen smoothed RQL-TV curvature;
- neither the regularization rule nor the source normalization is recomputed as a function of `w`.

Then

\[
D_s(w)=\log\det X_s(w)
\]

is concave. For a fixed task weight `W_s\succeq0`,

\[
R_s(w)=\operatorname{tr}\{W_sX_s(w)^{-1}\}
\]

is convex, because it is a matrix-fractional function; therefore

\[
U_s(w)=-R_s(w)
\]

is concave. These are standard facts from matrix convex analysis; see Boyd and Vandenberghe, *Convex Optimization*, Chapters 3 and 7: [official text](https://web.stanford.edu/~boyd/cvxbook/).

For concave functions `U_s`, the pointwise minimum is concave:

\[
\min_sU_s(\theta w+(1-\theta)v)
\ge
\theta\min_sU_s(w)+(1-\theta)\min_sU_s(v).
\]

Thus

\[
\max_{w\in\Delta_K}\min_sU_s(w)
\]

is a valid convex-optimization problem in the standard sense of maximizing a concave objective.

### Conditions that would break P-A

Concavity is lost or no longer established if any of the following is allowed inside the solve:

1. `W_s` depends on `w`;
2. `\lambda_{TV}` or the IRLS curvature is reselected after changing `w`;
3. bank power schedules are renormalized jointly after mixing;
4. `F_{k,s}` is recomputed from a reconstruction that itself depends on the proposed mixture;
5. a nonconcave image-domain surrogate is substituted for A-risk.

Freeze `H_0`, `W_s`, bank rows, and row-wise source settings before the simplex solve.

### Runtime caveat

`K=8` does not by itself imply millisecond runtime. A function/gradient evaluation requires up to `S=16` Cholesky factorizations or inverse solves on the declared task subspace. The R24 latency Gate D remains the deciding measurement; do not advertise subsecond runtime before benchmarking optimized batched linear algebra.

---

## 2. P-B: what dissolves and what does not

The threshold rule can be retired. The discontinuity claim cannot be promoted to a theorem.

By Berge’s maximum theorem, the optimal value is continuous under continuous scenario matrices and compact `\Delta_K`. The **argmax** is only upper hemicontinuous in general. It can jump when:

- two allocations tie;
- the least-favorable draw changes;
- an active bank enters or leaves the support;
- exact row rounding crosses an integer boundary.

Even if every scenario utility is strictly concave, the pointwise minimum need not be strictly concave. Integer materialization is necessarily piecewise constant.

### Parameter-free tie-break

Use a two-stage lexicographic definition:

1. solve the robust primary objective and obtain its exact optimal set `\mathcal W^*`;
2. choose
   \[
   \boxed{
   w^\dagger=\arg\min_{w\in\mathcal W^*}\|w-e_0\|_2^2.
   }
   \]

This secondary problem is strictly convex. It selects the knob corner when the robust objective is genuinely indifferent, without a `1.02` preference threshold or a temperature parameter.

It still does **not** guarantee global continuity. Describe the allocation as stable only where the robust optimum is unique and separated.

### Hedging is not abstention

A mixture commits measurement rows to several banks. It is robust diversification, not abstention. If all posterior draws share the same bias or every library member is poor, the minimax solution still returns the least-bad mixture. The only honest fallback is triggered by an implementation/physical-feasibility failure, not by a semantic claim that uncertainty has “emergently abstained.”

---

## 3. Failure modes of finite-draw minimax

Raw minimax introduces four important risks:

1. **one-draw domination:** a single rough, nearly singular, or bootstrap-pathological draw can determine all weights;
2. **scale domination:** an intrinsically difficult draw with large absolute A-risk can dominate even when its *relative* design sensitivity is ordinary;
3. **coverage failure:** draws centered on a biased pre-scan can all miss the true scene;
4. **model failure:** robustness over `x_s` does not cover wrong dead-time, background, motion, or reconstruction curvature.

`S=16` is an empirical scenario set, not a formal coverage guarantee. Classical scenario theory requires sample counts that depend on decision/support dimension and a declared violation probability; with roughly eight optimization variables, `N=16` gives essentially vacuous classical bounds at `10–20%` violation levels. Moreover, parametric-bootstrap draws conditional on one pre-scan are not automatically samples from the physical out-of-sample distribution. See Campi and Garatti’s exact scenario-feasibility result ([SIAM J. Optim. 2008](https://doi.org/10.1137/07069821X)) and the min-max scenario consistency analysis ([Ramponi 2017](https://doi.org/10.1137/16M109819X)).

Therefore the method may claim:

> robust over the declared finite posterior/bootstrap scenario set.

It may not claim distribution-free posterior coverage.

Mandatory no-gate disclosures are:

- active worst-scenario multipliers;
- number of active scenarios;
- leave-one-draw allocation spread
  \[
  \Delta_{\rm LOO}=\max_s\|w^\dagger-w^\dagger_{(-s)}\|_1;
  \]
- continuous versus materialized objective loss.

---

# Q2 — Correct robust objective

## 4. Reject raw maximin utility; adopt standardized maximin A-risk

The raw problem

\[
\min_w\max_sR_s(w)
\]

compares absolute risks on different scenario scales. The hardest draw can dominate simply because its baseline inverse problem is harder, not because the allocation is fragile.

Use **standardized maximin relative risk**. For each draw compute its bank-simplex oracle

\[
R_s^*=\min_{v\in\Delta_K}R_s(v).
\]

Then solve

\[
\boxed{
\begin{aligned}
\min_{w\in\Delta_K,t}\quad &t\\
\text{s.t.}\quad &R_s(w)\le tR_s^*,\qquad s=1,\ldots,S.
\end{aligned}}
\tag{R25.1}
\]

Here `t\ge1`; `t-1` is the worst relative A-risk regret over the sampled scenes. This is one parameter-free convex program.

This criterion is not a new generic design principle. It is the established standardized-maximin / maximin-efficiency construction; see Dette, “Designing Experiments with Respect to Standardized Optimality Criteria” ([JRSS B 1997](https://doi.org/10.1111/1467-9868.00056)) and Imhof and Wong, “A Graphical Method for Finding Maximin Efficiency Designs” ([Biometrics 2000](https://doi.org/10.1111/j.0006-341X.2000.00113.x)). The specialized content is the dead-time/RQL information matrix, nonnegative SPI bank, physical guards, and two-shot implementation.

Use `W_s=I` on the frozen reconstruction subspace for the first bridge, as in R24. `W_s` may vary by draw but must remain PSD and fixed during the solve.

### Log-det ablation

For D-optimality define

\[
D_s^*=\max_{v\in\Delta_K}\log\det X_s(v)
\]

and solve

\[
\min_{w\in\Delta_K}\max_s\{D_s^*-\log\det X_s(w)\}.
\]

The worst-case D-efficiency is `\exp(-t/r)`. This is an ablation; the estimator-aware A-risk criterion is the default.

## 5. Why not CVaR, expectation, or mean-variance?

### CVaR / soft-min

Lower-tail CVaR can preserve convexity, but it requires a tail level `\alpha`. At empirical level `\alpha=1/S`, CVaR is exactly the minimum and adds nothing. Any less extreme level is a new scientific robustness preference. Do not introduce it before the bridge.

### Posterior expectation

\[
\min_w\frac1S\sum_sR_s(w)
\]

is convex and is the natural Bayes action for the approximate posterior, but it is not conservative. Keep it as a no-gate comparison, not the default.

### Mean-variance

A mean-minus-variance utility generally destroys convexity and introduces a penalty coefficient. Reject it.

### Frozen default

The R25 default is therefore:

1. standardized maximin A-risk, Eq. (R25.1);
2. lexicographic closest-to-knob tie-break;
3. no runtime risk-level, hysteresis, or route-margin threshold.

---

## 6. Exact simplex certificate

Let

\[
r_s(w)=R_s(w)/R_s^*.
\]

The derivative is

\[
\frac{\partial r_s}{\partial w_k}
=-\frac{\operatorname{tr}\{W_sX_s(w)^{-1}F_{k,s}X_s(w)^{-1}\}}
{R_s^*}.
\]

Define the positive marginal value

\[
g_{s,k}(w)=
\frac{\operatorname{tr}\{W_sX_s(w)^{-1}F_{k,s}X_s(w)^{-1}\}}
{R_s^*}.
\]

At an optimum there exist scenario multipliers `\alpha_s\ge0`, simplex multipliers `\mu_k\ge0`, and `\lambda` such that

\[
\sum_s\alpha_s=1,
\qquad
\alpha_s\{r_s(w)-t\}=0,
\]

\[
\sum_s\alpha_sg_{s,k}(w)
\begin{cases}
=\lambda,&w_k>0,\\
\le\lambda,&w_k=0.
\end{cases}
\tag{R25.2}
\]

The `\alpha_s` are the least-favorable empirical scenario distribution. Equation (R25.2), primal feasibility, and complementary slackness give an exact certificate for the continuous scenario problem.

After row materialization report separately:

1. continuous KKT residual;
2. realized worst relative regret;
3. materialization regret `t_real-t_cont`;
4. all physical-guard residuals.

The certificate remains a local estimator-risk certificate, not an image-quality certificate.

---

# Q3 — Row mixing and the required T-B redesign

## 7. The strongest objection is valid

If `L_k` is merely a fixed, jointly optimized 972-row exact design, then

\[
F_k=\sum_{r=1}^{972}H_{k,r}
\]

does **not** imply that an arbitrary subset of `m` rows has information `(m/972)F_k`. Its whole-design value can come from a carefully balanced support and spectrum. A prefix or random fraction can be rank-deficient, dose-unbalanced, or spectrally poor.

Therefore the literal rule

```text
n_k = round(972*w_k)
take n_k rows from L_k
```

is rejected.

## 8. Required library representation

Each entry must be built and stored as a **mixable approximate design measure**

\[
\xi_k=\{(a_{k,r},p_{k,r})\}_r,
\qquad p_{k,r}\ge0,\quad\sum_rp_{k,r}=1,
\]

with per-row:

- physical binary/gray pattern;
- source setting/load;
- incident cost;
- pixel-dose vector;
- peak and detector-admission metadata.

For draw `s`, define the per-972-row information

\[
F_{k,s}=972\sum_rp_{k,r}H_s(a_{k,r}).
\]

The robust allocation induces one union measure

\[
\boxed{
\xi(w)=\sum_kw_k\xi_k.
}
\tag{R25.3}
\]

At this continuous-measure level, information, incident budget, and cumulative dose are affine in `w`. If every `\xi_k` satisfies the same linear guard, the mixture does too. Spectral LMIs inherit by convexity. A-risk caps inherit because `tr(WX^{-1})` is convex.

### Immediate T-B ruling

- `L0–L6` may continue only if their generators expose row atoms and normalized weights and support deterministic materialization at arbitrary requested counts.
- **Pause the current exact-bank freeze of `L7`.** Store the offline FW support and weights as `\xi_7`; do not treat one optimized 972-row realization as fractionally mixable.
- If a bank cannot be represented this way without losing its defining property, it is a nonmixable corner and cannot participate in the simplex. Keep it only as a standalone baseline.

## 9. Exact 972-row materialization

Merge identical physical rows across banks and form target expected counts

\[
p_a=972\,\xi(w)(a).
\]

Materialize one exact sequence by constrained rounding over the union support:

\[
\sum_an_a=972,\qquad n_a\in\mathbb Z_{\ge0},
\]

while enforcing:

- total incident budget;
- the per-pixel dose band;
- row-wise peak/admission constraints;
- any frozen family-allocation restriction.

Efficient rounding of approximate designs controls ordinary design-efficiency loss but does not by itself enforce all side constraints; see Pukelsheim and Rieder ([Biometrika 1992](https://doi.org/10.1093/biomet/79.4.763)). Exact/constrained design can be handled by mixed-integer conic methods ([Sagnol and Harman, Annals of Statistics 2015](https://www.zib.de/userpage/sagnol/publication/sh15_aos/)) or a balanced-sampling/dependent-rounding stage followed by deterministic exchange ([Deville and Tillé, Biometrika 2004](https://doi.org/10.1093/biomet/91.4.893)).

For this prototype, use the repository’s deterministic exchange machinery over the union atom pool; do **not** reintroduce online full-dictionary OED. The task is rounding a fixed `K=8` mixture, not optimizing a new scene-specific design.

Post-materialization:

1. recompute actual dose, budget, spectral and A-risk guards;
2. recompute all `X_s` from selected rows;
3. report continuous and realized relative regret;
4. if no compliant exact design exists, return `L0` and flag
   ```text
   MIXTURE_MATERIALIZATION_FAIL
   ```

This is an operational honesty fallback, not an uncertainty threshold.

## 10. Shared rows and ordering

Banks sharing a physical row do not create statistical correlation under the frozen conditionally independent active-start frame model. Merge duplicate atoms and add their target weights. Repeated exposures are legitimate, though possibly inefficient.

Randomize or balance the final block order with a frozen seed so source drift is not confounded with bank identity. If detector state carries across patterns, afterpulsing spans frames, or the scene changes between blocks, the additive-Gramian model is no longer exact; that is outside this first bridge and must be declared.

---

# Q4 — Gate C reinterpretation

## 11. Primary allocation diagnostic

For the **realized** 972-row design, let

\[
\hat w_{j,k}=n_{j,k}/972
\]

be the row fraction attributed to library component `k` for scene `j`. Define

\[
\boxed{
A_j=1-\hat w_{j,0}
=\tfrac12\|\hat w_j-e_0\|_1,
}
\tag{R25.4}
\]

which is the physical fraction of the acquisition moved away from the knob bank.

Use `A_j`, not entropy, for the gate. Entropy cannot distinguish useful departure from the knob from gratuitous spreading among redundant banks.

Report, but do not gate on,

\[
H_j=-\frac{\sum_k\hat w_{j,k}\log\hat w_{j,k}}{\log K},
\]

active-bank count, worst-scenario multipliers, and `\Delta_LOO`.

## 12. Frozen Gate C replacement

Keep the existing DEV-only ORACLE-LIB labels and the existing `1 dB`, `0.25 dB`, `80%`, and `75%` constants. Introduce no new threshold.

### C1 — rescue-needed scenes

Among scenes where ORACLE-LIB beats the knob by more than `1 dB`, require

\[
\boxed{
\operatorname{mean}_j A_j\ge0.80.
}
\]

### C2 — aligned controls

Among aligned controls where ORACLE-LIB advantage is below `0.25 dB`, require

\[
\boxed{
\operatorname{mean}_j\hat w_{j,0}\ge0.75.
}
\]

This is the exact continuous relaxation of the former categorical route fractions:

- binary bank selection gives `A_j\in\{0,1\}` and recovers the original `80%` rule;
- knob/abstain selection gives `w_{j,0}=1` and recovers the original `75%` rule.

Use materialized weights, not the continuous optimizer output.

If Gates A and B pass but this Gate C fails, the method may still be an effective **static composite design**, but it has not validated scene-adaptive allocation and must not lead M2 as an adaptive method.

Gate A, Gate B, and Gate D remain numerically unchanged.

---

# Q5 — Net architecture ruling

## 13. Adopt, with the amendments above

The final R25 Step 3 is:

1. produce `S=16` declared posterior/bootstrap scenario matrices;
2. compute each scenario’s bank-simplex A-risk oracle `R_s^*`;
3. solve the standardized maximin problem (R25.1);
4. among exact minimizers, choose the allocation closest to `e_0`;
5. materialize 972 rows from the union design measure with constrained rounding/exchange;
6. recompute guards, realized regret, and the simplex KKT certificate;
7. use `L0` only on solver/materialization/guard failure, with an explicit flag;
8. acquire once and reconstruct with RQL.

No runtime `0.01`, `1.02`, `2%`, hysteresis, or uncertainty-abstention threshold remains.

The R23 quantities `\varepsilon_1`, `m`, and `\varepsilon_1^2/(2m)` remain mandatory disclosures. They explain whether rank one was locally information-sufficient; they do not decide the allocation.

## 14. What is genuinely improved

The amended allocation is better than the R24 rulebook because it replaces:

- three unrelated route thresholds;
- categorical switching;
- ad hoc max-envelope comparison;
- nominal “abstention”;

with one established robust-design object and one exact KKT certificate.

It is also more honest. It does **not** promise a continuous optimizer, posterior coverage from 16 draws, or image-quality certification.

## 15. Brutal stop condition

The architecture is worth building **only if T-B is changed now**. A naive prefix/subset mixture of eight monolithic exact banks is mathematically inconsistent with the objective being optimized. If the library cannot expose mixable measures and a compliant exact 972-row materializer within Gate D latency, abandon the mixture architecture and return to the single-bank R24 router. Do not publish a convex continuous allocation whose physical row realization is a different design.

The DEV Fisher-to-image Gates A and B remain the ultimate kill test. Convex elegance does not rescue a surrogate that fails to improve images.