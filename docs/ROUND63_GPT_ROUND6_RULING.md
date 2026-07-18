# GPT round-6 frozen ruling (2026-07-18, 15m21s reply) — analytic λ rule + geometry/grid + censoring

Reconstructed from screenshots (the relay's DLP filter blocked several DOM
spans for tool extraction; formulas transcribed visually — on any discrepancy
the ChatGPT thread governs). THE operative spec for F1 alongside rounds 4–5.

## Summary table (verbatim)
| 项目 | 裁决 |
|---|---|
| R1 废弃 DEV discrepancy 选 λ_TV | 批准 |
| 全局 c vs 自适应 g | 极简自适应规则；不用原始 GI-TV，改用物理可解释的 bucket excess-variance 浓集度 |
| R2 主几何改为 M/n=1 | 批准 |
| R3 加入 ν=5,10 | 批准，但保留 ν=1000 |
| R4 原主门 | 保留，不降门槛 |
| floor censoring | 固定规则，24 图一张不删 |
| R5 重跑 S1 | 批准；分"校准 pass"和"冻结后验证 pass"两步 |

## 1. F1 frozen λ_TV rule — `analytic_score_concentration`

**1.1 Deleted from the production path entirely**: lam_max_arm; K=5 cross-fit
discrepancy path; one-SE η\*; per-image 9-η reconstructions. The DEV/AUDIT
80/20 split SURVIVES but: DEV is used only to compute the concentration
statistic; AUDIT only for the round-5 descriptive diagnostics; neither runs
any reconstruction hyperparameter search.

**1.2 Per-arm analytic noise scale.** Each data term f_a(x) =
(1/M)Σ q_a(N_i; λ_i); per-frame score s_a(N; λ) = ∂q_a/∂λ. E.g.
q_RQL = (T−Nτ)λ − N log λ, s_RQL = T − Nτ − N/λ. PRECORRECT's delta-method
weights w(N) keep their formula but the normalization constant becomes the
exact-PMF E[w_raw(N)] (not the sample mean). Score variance under the EXACT
active-start NP renewal law at λ̄ = Φ + d̂:
  v_{s,a} = Var_{N~P_NP(·|λ̄,T,τ̂)}[ s_a(N; λ̄) ]
(scalar enumeration over m, negligible cost; at ν=5,10 the CLT variance is
FORBIDDEN — exact PMF moments only). Pattern column-energy factor:
  κ_A = max_{1≤j≤n} (1/M) Σ_i a_ij²    (≈2 for bern50 {0,2})
Then (boxed, frozen):
  σ_{g,a} = Φ √(κ_A · v_{s,a} / M)
  λ_{TV,a} = c_i · σ_{g,a} · √(2 ln n)
For bern50 the RQL case reduces to σ_{g,RQL} ≈ (Φ/M)√(2M·V_N(λ̄)·(τ+1/λ̄)²)
(matches the session's ad-hoc scan formula). ALL arms share the same c_i —
no per-arm c tuning.

**1.3 Why NOT raw TV(x_GI)**: GI's TV mixes real object detail, (ρ,ν)
measurement noise, and pattern-energy fluctuation; at low photons noise makes
scenes look "rough" → picks SMALLER c exactly when more regularization is
needed. Freeze instead a statistic that directly measures the multiplex
disadvantage.

**1.4 Object concentration statistic (DEV raw counts only).**
  S_N² = (1/(M_dev−1)) Σ_{i∈DEV} (N_i − N̄)²
At λ̄ = Φ + d̂, from the exact active-start PMF: μ₀ = E[N], V₀ = Var N,
μ′₀ = ∂E[N]/∂λ. Pattern covariance mean scale:
  ω_A = (1/(M_dev n)) Σ_{i∈DEV} Σ_j (a_ij − ā_j)²
Then (boxed, frozen):
  Ĉ = clip[ n · (S_N² − V₀)₊ / ((Φ μ′₀)² ω_A), 1, 64 ]
Interpretation: for unit-mean iid patterns Var(aᵀx) ≃ Σx², so C = n Σx²;
uniform/diffuse ≈ 1; text/dot/concentrated structures larger. If S_N² ≤ V₀:
fix Ĉ = 1 (smooth-scene bin); NEGATIVE excess variance may not extrapolate.
Per cell log: C_hat, S_N2, V0_exact, muprime0_exact, omega_A, c_used,
sigma_grad_arm, lambda_tv_arm.

**1.5 Minimal two-bin c rule** (no continuous nets/regressions; boxed):
  c_i = 0.50 if Ĉ_i ≤ C₀ ;  0.25 if Ĉ_i > C₀
(preserves the two stable operating points seen in the session scans).

**1.6 C₀ calibration (Pass A, ONCE).** Candidates
C₀ ∈ {0, 1.5, 2, 3, 4, 6, 8, 12, 16, ∞}. Per DEV image j, endpoint-oracle
regret R_j(C₀) = mean_{ρ,ν,seed}[ max{P_{j,.25}, P_{j,.50}} − P_{j,c(Ĉ;C₀)} ]
with P = PSNR_rad. Objective J(C₀) = Q_{0.90}{R_j(C₀) : j}. Frozen choice
C₀* = argmin J; tie rule: within 0.02 dB take the LARGER C₀ (prefer stronger
regularization / the simpler global-0.5 limit). This calibration explicitly
uses DEV ground truth, exactly once, fully disjoint from the STL-10 test
confirmatory set. Paper wording: "analytically noise-scaled TV with one
development-calibrated concentration threshold" — never "training-free" /
"zero-calibration".

## 2–3. Geometry and grid
- Primary S2-A: **M/n = 1 (M = 4096 at 64²)**; M/n ∈ {0.25, 0.5} demoted to
  the S2-B robustness axis (undersampling ceiling = reported phenomenon).
- **ν grid: {5, 10, 20, 50, 100, 200, 500, 1000, 2000}** (9 tiers; 1000 kept).
- ν < 20 is called the **short-window stress region** in the paper; do NOT
  claim RQL ≈ exact likelihood there; RQL stays the current convex
  quasi-likelihood (no temporal-mixture estimator); its short-window
  approximation error is quantified by S4 (S4 exact-vs-RQL suite gains
  ν ∈ {5,10}).
- σ_g and Ĉ always use exact finite-window PMF moments.

## 4. Main gate (R4) + frozen floor-censoring taxonomy
Gate unchanged: 24-image median S_j ≥ 3; image-level stratified bootstrap 95%
lower bound > 1; ≥ 18/24 images with S_j > 1. Not lowered for anticipated
censoring. T_min = M_physical·τ·ν_min, ν_min = 5. Per image:
- **A** both crossings interior → normal S_j, all gates.
- **B** fast reaches floor, safe later: only S_j ≥ T_s/T_min is known →
  status=FAST_LEFT_CENSORED, S_gate = T_safe/T_min; main table shows
  S_j ≥ S_gate; the conservative bound enters median, 18/24 and bootstrap.
- **C** BOTH at floor (crossing outside observed range, ratio unidentifiable):
  status=BOTH_LEFT_CENSORED, S_gate = 1.0; enters the 24-image median; does
  NOT count toward the S_j>1 tally; not deleted/excluded; main table marks
  "unresolved below floor" — never claimed as a measured 1×.
- **D** safe at floor, fast interior (pathological direction):
  status=SAFE_LEFT_FAST_INTERIOR, S_gate = T_min/T_f < 1 (descriptive upper
  bound); counts as a NON-positive result.
- **E** fast never reaches Q* by ν=2000: counted as failure, NOT removed from
  the median.
- **F** safe fails / data missing (should be impossible by Q*_j definition):
  status=ANALYSIS_FAILURE, S_gate = 0.0; triggers a cell-completeness audit.
**Bootstrap**: every draw must (1) resample seeds within image, (2) re-form
both PSNR–time curves, (3) RE-JUDGE the censoring type, (4) regenerate
S_gate, (5) resample at the image level. Censoring statuses may NOT be frozen
per-image before bootstrapping.

## 5. R5 — final S1 rerun, two fixed passes
- **Pass A (the single C₀ calibration)**: run endpoints c = .25/.50 per §1.6,
  choose C₀*. Deliver c_calibration.csv, c_threshold_regret.json, C0_FROZEN,
  calibration image/seed list, code/input SHA256 → commit. Afterwards the c
  endpoints, C₀, formula, and calibration objective may not change.
- **Pass B (final protocol verification)**: frozen C₀*; M=4096; all 9 ν;
  ρ={0.05,0.6}; the 6 S1 dev images; seeds {0,1}; arms RQL/POISSON-LIN/
  SAT-POISSON/PRECORRECT/GI; 200 iterations; NO per-image search. Purpose
  ONLY: confirm curves are no longer TV-null flat; confirm an interpretable
  photon-limited region exists; calibrate Colab shard time; check numerical
  stability/completeness.

## F1 must update before freeze
select_rule = "analytic_score_concentration"; M_S2A = 4096;
NU_S2A = [5,10,20,50,100,200,500,1000,2000]; regenerate manifests,
expected-cell table, shard cost model, frozen SHA ledger, and the S1/S2
analysis censoring logic.

## Final framing (verbatim intent)
S1 proved the original protocol was dominated by regularization-selection
non-identifiability AND a spatial undersampling ceiling, hence could not test
the preregistered photon-time hypothesis; F1 moves the primary claim to the
identifiable M/n = 1 region and pins regularization strength via
detector-calibrated analytic score noise. After this change ROUND63 measures
"whether high flux saves exposure time" — not "whether cross-validation
likes a flat image".
