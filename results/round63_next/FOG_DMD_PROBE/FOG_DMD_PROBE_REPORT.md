# FOG-as-a-second-DMD — GPU feasibility probe (R38-aligned)

**Date:** 2026-07-23 · **Status:** DEV-GRADE exploration.
**Verdict arc:** (1) *pathwise blind* cold-start **F3 NOT MET** (§7); (2) two init-based fixes
**BOTH_FAIL** (§11) — the sharper diagnosis is that the drift is at the *medium re-estimation*
step (collusion with the free realizations `z_t`), not the initialization (§11c); (3) **RESURRECTION
(§12):** removing the realizations and estimating in the **covariance domain** (R38 eq 1.18) is
drift-free, injective, photon-robust, and **recovers the fixed null space** — verdict
**COVARIANCE_STICKS** — at a quantified temporal-sample price (crosses the 0.25 bar at T≈2048 vs the
pathwise T=128 that fails). The concept is **alive in statistical-identifiability form.**

## 0. Where it ran, and why (honest disclosure)

Ran **locally on an NVIDIA RTX 4060 Laptop GPU** (Anaconda `base`, Python 3.11.5,
torch 2.2.1+cu121, numpy 1.24.4), **float64** throughout. Colab was *not* used.
Rationale, per the standing "science over venue" fallback: the problem is small
(N=1024; even the E3 scale N=4096 fits this GPU), the decisive obstacle is the
**bilinear cold-start optimizer**, which needed many fast iterate–measure cycles to
develop — a loop that Colab's token-churn/rebind discipline would have throttled.
All heavy linear algebra is on-GPU. The metric is recomputed in float64 on CPU.

A key speed lever made the whole grid run in minutes: for uniform weights the scene
normal matrix factors as **AᵀA = (PᵀP) ⊙ (WᵀW)** (Hadamard of two N×N matrices),
turning the x-solve from O(T·M·N²) into O(T·N²). The Poisson-weighted x-step is
replaced by unweighted LS + Tikhonov (a dev-grade approximation; Poisson statistics
still enter through the noisy buckets and the Kalman R-diagonal in the medium E-step).

## 1. What this probe tested, and the R38 course-correction

Original idea ("fog as a second DMD"): a bucket detector, a **reused** known binary
pattern bank {P_i} (M rows), a static scene x (N pixels), and an **unknown**
spatially-structured medium field w_t fluctuating across epochs. Buckets
B_{i,t} = ⟨P_i ⊙ w_t, x⟩ (+ Poisson noise). Claim: medium fluctuation rotates the
fixed row space and makes the fixed-optics **null space** recoverable by **blind**
joint estimation of (x, {w_t}) — no reference arm.

Mid-probe, the GPT **R38 ruling** (docs/ROUND63_GPT_ROUND38_RULING_RAW.md) landed and
analytically reshaped the test. Adopted in full:

- **Oracle rank ≠ blind information.** The exact blind local phase boundary is
  `d<M` and `min(T,d+1)(M−d) ≥ N−g`. At **M=64, N=1024** the best capacity is
  `(31+1)(64−31)=1056`, only 3.1% above N (near-square, photon-hostile).
  ⇒ **M=64 runs relabelled as NEGATIVE CONTROLS** that must *reproduce* the theorem.
- **Primary cells: M=96**, d∈{32,48,64}, T=128 (χ_full ≈ 2.0–2.3), N=1024.
- **Normalization fix:** hold pixelwise field RMS **σ_f fixed across d**
  (per-coeff sd = σ_f·√(N/d)); positive **lognormal** medium w_t = exp(g_t − v/2),
  g_t = U z_t. Verified: DCT basis orthonormal to 1e-15, E[w]=1.0000, RMS(log w)=σ_f
  exactly (0.150 / 0.301) independent of d.
- **Staged, data-only cold start** (no truth access): Stage A nonneg row-space scene
  from temporal-mean buckets; Stage B centered spectral init
  B_c = B(I−11ᵀ/T) = P·diag(x)·U·Z_cᵀ → temporal subspace; Stage C tracking EM
  (OU Kalman smoother E-step + nonnegative LS+Tikhonov M-step), Tikhonov **homotopy**,
  multiple data-only random starts.
- **Photons:** primary headline ≤1e5; 1e6 as ceiling; 1e4/1e3 for the low-dose edge.
- **Metric: null-NMSE** = (‖null-component error‖ / ‖null_true‖)² after the scalar gauge;
  fixed-averaging baseline = 1.0. Report medians over ≥5 seeds, oracle-improvement
  fraction captured, and an OU-prior-weakened-4× sensitivity.

**Success signal (dev shadow of ruling bar F3):** blind cold-start null-NMSE **≤ 0.25
median** at any M=96 cell **AND ≥ 50%** of the oracle improvement captured.

## 2. Setup reproduction (sanity — matches the established toys exactly)

Using the original 32×32 scene, M=64, T=64, iid σ=0.30, signed-Gaussian P (the toy
setup), my common code reproduces the established numbers exactly:

| d_w | stacked rank | oracle | ALS warm (truth+15%) | ALS cold |
|----:|:---:|:---:|:---:|:---:|
| 16  | 1024/1024 | 0.033 | 0.558 | 1.000 |
| 64  | 1024/1024 | 0.000 | 0.133 | 0.991 |
| 256 | 1024/1024 | 0.000 | 0.100 | 1.002 |

(amplitude null-err; brief quoted oracle 0.033–0.000, warm 0.10–0.56, cold ~1.0). ✔

## 3. Why blind cold-start is hard — identifiability probe (M=64, T=256, iid)

Before the pivot I instrumented the failure directly (exp_identify.py):

- **ORACLE (true W):** null-err 0.161, data-residual 1.85e-5 → the data *do* determine
  the null when W is known.
- **SPECROT blind cold:** null-err 0.998 but **‖x̂‖ = 20.3** (vs oracle 8.6) — it fits
  the data with a **blown-up, wrong-direction** null solution.
- **ALS warm-started from the TRUE medium mapping:** *drifts* 0.161 → **0.549**, and its
  residual (1.13e-5) is **lower** than the oracle's — a worse-null solution fits *better*.
- **range-only x:** residual 1.41 (cannot fit).

Reading: under the plain least-squares blind objective the minimum is **not at truth**;
near-square geometry admits large-norm spurious minima. This is exactly the R38
"nonidentifiable-fiber / near-square Fisher edge" picture. The two fixes it implies —
**capacity margin (M=96)** and **physical regularization (nonnegativity + Tikhonov)** —
are what the primary probe applies.

## 4. Negative controls (M=64) — reproduce the theorem ✔

Lognormal σ_f=0.30, T=128, physical 0/1 patterns, 5-seed medians. null-NMSE:

| control | χ_full | oracle | blind cold (median) | warm-Q ceiling |
|:--|:--:|:--:|:--:|:--:|
| scalar gain (d=1) | 0.12 | **1.000** | 0.999 | — |
| d=16 | 0.80 | 0.000 | 0.698 (cap 30%) | 0.857 |
| d=31 (near-square) | 1.03 | 0.000 | 0.738 (cap 26%) | 0.996 |
| d=64 (=M) | 0.00 | 0.000 | 0.854 (cap 15%) | — |

All four behave as R38 predicts: scalar medium = zero-diversity degeneracy (oracle
itself stuck at 1.0); the oracle sees the null at d≥16 (full stacked rank) but **blind
fails everywhere at M=64**, worst at d=64 (no free-path profiled information) and at the
near-square χ≈1 edge the warm ceiling is ~1.0. **Blind M=64 dense-scene recovery is
dead — confirmed empirically, matching the analytic kill.**

## 5. Primary blind-capacity cells (M=96)

_[table filled from E1_results.json below]_

**Lognormal medium (ruling-mandated), physical 0/1 patterns, T=128, 5-seed medians.**

| cell | χ_full | photons | oracle | lin-approx oracle | **blind cold (median)** | captured | warm-Q ceiling |
|:--|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| M96 d32 σ_f0.15 | 2.06 | clean | 0.000 | 0.006 | **0.586** | 41% | 0.826 |
| M96 d32 σ_f0.15 | 2.06 | 1e+06 | 0.010 | 0.016 | **0.740** | 26% | — |
| M96 d32 σ_f0.15 | 2.06 | 1e+05 | 0.022 | 0.029 | **0.689** | 32% | — |
| M96 d32 σ_f0.15 | 2.06 | 1e+04 | 0.131 | 0.140 | **0.979** | 2% | — |
| M96 d32 σ_f0.3 | 2.06 | clean | 0.000 | 0.028 | **0.715** | 28% | 0.946 |
| M96 d32 σ_f0.3 | 2.06 | 1e+06 | 0.001 | 0.027 | **0.645** | 36% | — |
| M96 d32 σ_f0.3 | 2.06 | 1e+05 | 0.005 | 0.032 | **0.905** | 10% | — |
| M96 d32 σ_f0.3 | 2.06 | 1e+04 | 0.041 | 0.068 | **0.940** | 6% | — |
| M96 d48 σ_f0.15 | 2.30 | clean | 0.000 | 0.008 | **0.630** | 37% | 0.261 |
| M96 d48 σ_f0.15 | 2.30 | 1e+06 | 0.006 | 0.016 | **0.582** | 42% | — |
| M96 d48 σ_f0.15 | 2.30 | 1e+05 | 0.018 | 0.026 | **0.769** | 24% | — |
| M96 d48 σ_f0.15 | 2.30 | 1e+04 | 0.121 | 0.131 | **1.019** | -2% | — |
| M96 d48 σ_f0.3 | 2.30 | clean | 0.000 | 0.034 | **0.544** | 46% | 0.905 |
| M96 d48 σ_f0.3 | 2.30 | 1e+06 | 0.001 | 0.033 | **0.674** | 33% | — |
| M96 d48 σ_f0.3 | 2.30 | 1e+05 | 0.004 | 0.037 | **0.733** | 27% | — |
| M96 d48 σ_f0.3 | 2.30 | 1e+04 | 0.037 | 0.071 | **0.865** | 14% | — |
| M96 d64 σ_f0.15 | 2.03 | clean | 0.000 | 0.012 | **0.583** | 42% | 0.203 |
| M96 d64 σ_f0.15 | 2.03 | 1e+06 | 0.005 | 0.016 | **0.548** | 45% | — |
| M96 d64 σ_f0.15 | 2.03 | 1e+05 | 0.015 | 0.026 | **0.795** | 21% | — |
| M96 d64 σ_f0.15 | 2.03 | 1e+04 | 0.112 | 0.123 | **0.988** | 1% | — |
| M96 d64 σ_f0.3 | 2.03 | clean | 0.000 | 0.036 | **0.676** | 32% | 0.999 |
| M96 d64 σ_f0.3 | 2.03 | 1e+06 | 0.001 | 0.035 | **0.696** | 30% | — |
| M96 d64 σ_f0.3 | 2.03 | 1e+05 | 0.004 | 0.038 | **0.641** | 36% | — |
| M96 d64 σ_f0.3 | 2.03 | 1e+04 | 0.034 | 0.069 | **0.814** | 19% | — |

Best primary blind cold-start: **NMSE 0.544** at M96 d48 sf0.3 ph=clean (captures 46% of the oracle improvement).

### 5a. Model-matched linear-medium reference (isolates cold-start vs mismatch)

Same geometry, affine medium w=1+Uz that the estimator models exactly. The warm-Q column is now the true **profiling-tax ceiling** (no model mismatch).

| cell | χ_full | photons | oracle | lin-approx oracle | **blind cold (median)** | captured | warm-Q ceiling |
|:--|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| M96 d32 σ_f0.3 | 2.06 | clean | 0.000 | 0.000 | **0.782** | 22% | 0.178 |
| M96 d48 σ_f0.3 | 2.30 | clean | 0.000 | 0.000 | **0.525** | 48% | 0.084 |
| M96 d64 σ_f0.3 | 2.03 | clean | 0.000 | 0.000 | **0.590** | 41% | 0.100 |

Reading: with the model matched, the warm-start ceiling falls to 0.084–0.178 (vs 0.86–1.0 at M=64) — the capacity margin is real — yet blind **cold-start** stays 0.52–0.78. Cold-start is the wall.

### 5b. OU-prior sensitivity (medium OU τ=16; prior weakened 4×)

| cell | χ_full | photons | oracle | lin-approx oracle | **blind cold (median)** | captured | warm-Q ceiling |
|:--|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| M96 d48 σ_f0.3 | 2.30 | 1e+05 | 0.014 | 0.115 | **0.653** | 35% | — |
| M96 d48 σ_f0.3 | 2.30 | 1e+05 | 0.014 | 0.115 | **0.629** | 38% | — |

Weakening the OU precision 4× changes blind null-NMSE by -0.024 (0.653 → 0.629). The (small) recovery is **not** prior-imputed — it is likelihood-driven — but it remains far above the bar either way.


## 6. E2 — noise scaling curve (regularized)

_[filled from E2_noise_curve.json below; PNG: E2_noise_curve.png]_

Cell M=96, d=48, T=128, σ_f=0.3, lognormal; 5-seed medians. null-NMSE vs mean detected photons/bucket (baseline fixed-averaging = 1.0):

| photons/bucket | oracle NMSE | blind cold NMSE | baseline |
|:--:|:--:|:--:|:--:|
| 1e+03 | 0.245 | 1.061 | 1.000 |
| 1e+04 | 0.034 | 1.018 | 1.000 |
| 1e+05 | 0.007 | 0.859 | 1.000 |
| 1e+06 | 0.004 | 0.843 | 1.000 |
| 1e+07 | 0.003 | 0.795 | 1.000 |

The oracle degrades **gracefully** toward 1.0 as photons fall (regularized, no explosion): 0.003 at 1e7 → 0.245 at 1e3. The blind arm never approaches the F3 bar at any dose. See `E2_noise_curve.png`.


## 7. Verdicts

**E1 (blind cold-start recovery) — SUCCESS SIGNAL: NOT MET.**
Best primary blind cold-start null-NMSE = **0.544** (at M96 d48 sf0.3,
ph=clean), capturing **46%** of the oracle
improvement. The F3-shadow bar requires ≤0.25 **and** ≥50% captured; no cell meets it
(F3 pass = False). Across the whole M=96 grid, blind cold-start sits at NMSE ≈ 0.54–1.0.

- **Oracle viability holds** (null-NMSE ≈ 0 clean; ≤0.02 at ≤1e5 photons at d≤48): the data *do*
  contain the null information, and the M=96 capacity margin (χ≈2.0–2.3) is real — the
  *model-matched* warm-start profiling-tax ceiling (§5a) drops to **0.08–0.18** (vs 0.86–1.0 at
  M=64). (The lognormal-medium warm-Q in §5 is higher, 0.20–1.0, because the linear estimator
  mismatches the strongly nonlinear medium — a solver limitation, not a geometry limit.)
- **The wall is the bilinear cold-start optimizer.** Even where the warm ceiling reaches the bar
  (d=64, σ_f=0.15: warm-Q 0.203), the data-only staged cold start stalls at 0.58. This is R38's
  predicted kill-tree node: *F1/F2-type geometry OK, F3 cold-start fails → software-only blind
  deployment defeated.*

**E2 (noise scaling):** oracle degrades gracefully toward 1.0 under regularization; blind never
approaches the bar at any dose (1e3–1e7). Delivered as table above + `E2_noise_curve.png`.

**E3 (scale):** **not run** — gated on E1 showing life, which it does not.

**Bottom line for the direction:** consistent with the R38 ruling. The oracle "second-DMD" effect
is real and the M=96 geometry is healthy, but **reference-free blind cold-start recovery of the
fixed null space fails** at the dev-grade F3 bar. Under the ruling's kill tree this is a **kill of
the software-only blind deployment**; any surviving value is the oracle/known-medium-law theorem
(materials-bank) or a much lower-dimensional scene class — not the arbitrary-dense-scene flagship.


## 8. E3 — scale (64×64)

**Not run.** E3 was explicitly gated on "only if E1 shows life." E1 does not meet the
success signal, so escalating scene dimension would not change the verdict. The
infrastructure (fast Hadamard x-solve, GPU float64) is N=4096-ready if a future,
lower-dimensional-scene or higher-M regime revives the direction.

## 9. Failure modes observed

1. **Bilinear cold-start non-convexity is the dominant obstacle.** Naive cold ALS,
   cold joint-Adam, spectral-init ALS, null-kick symmetry breaking, and many random
   restarts all collapse to null-NMSE ≈ 1.0 at M=64 (exp_coldstart.py).
2. **Spurious large-norm minima** with wrong null direction (‖x̂‖≈20) — killed only by
   Tikhonov + nonnegativity, which is why the staged solver uses a ridge homotopy.
3. **Profiling tax is small at M=96 but cold-start still fails.** On the *model-matched* linear
   medium the warm-start reduced-ALS ceiling is **0.08–0.18** at M=96 (vs ~0.86–1.0 at M=64) —
   *below* the 0.25 bar — yet the data-only cold start stalls at 0.52–0.78. So the geometry and
   profiling tax are fine; the bilinear cold-start optimizer alone defeats the bar.
4. **Lognormal ⇄ linear model mismatch** inflates the linear estimator's medium fit at
   strong σ_f (coeff-sd up to 1.7); the *known-Z* linear-approx oracle stays small
   (≈0.03–0.05), so the mismatch is in blind *estimation*, not the model ceiling.
5. **Full-Z joint Adam is numerically unstable** (diverges even from truth); the
   subspace-reduced ALS (Z∈span of centered-SVD) is far better behaved.

## 11. E4 — cold-start fixes (operator-approved +2h extension)

Two attempts to break the pure-optimization wall on the best cell **M=96, d=48, σ_f=0.15,
T=128** (clean + 1e6 + 1e5 photons). Headline up front: **BOTH_FAIL**, and the extension
produced a *sharper* diagnosis of why — see §11c.

### 11a. E4a — classical rotation fix via multi-timescale medium + SOBI

Medium changed (declared, per R38 §1.1) so DCT components carry **scale-dependent** OU times
τ_k log-spaced 64→4 (coarse→fine; ρ 0.984→0.779), holding σ_f fixed. Then the lagged bucket
covariances R_b(ℓ)=H·diag(σ²ρ_kˡ)·Hᵀ are simultaneously diagonalizable *only* in the true
source coordinates, so a SOBI joint-diagonalization (whiten with R_b(0); Cardoso–Souloumiac
Jacobi joint-diag of lags 1,2,3,5,8,13,21,34) should fix the d×d rotation the SVD Stage-B
init cannot resolve. Recovered sources seed the reduced-ALS Q, then the usual likelihood
refinement. Single-τ (τ=16) is the control (all ρ_k equal ⇒ every R_b(ℓ)∝HHᵀ ⇒ joint-diag
degenerate ⇒ SOBI cannot help).

**Medium generated (ruling-declared), physical 0/1 patterns, 5-seed medians; null-NMSE. 'frozen-known' = solve x with the medium held at truth (no re-estimation).**

| medium | photons | oracle | frozen-known | true-Z-seed refine | baseline (random) | **SOBI-fix** | SOBI quality |
|:--|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| single | clean | 0.000 | 0.029 | 0.669 | 0.652 | **0.625** | 0.26 |
| single | 1e+06 | 0.034 | 0.062 | 0.782 | 0.679 | **0.678** | 0.24 |
| single | 1e+05 | 0.053 | 0.077 | 0.683 | 0.643 | **0.648** | 0.22 |
| multi | clean | 0.000 | 0.029 | 0.642 | 0.700 | **0.708** | 0.28 |
| multi | 1e+06 | 0.030 | 0.061 | 0.647 | 0.703 | **0.651** | 0.26 |
| multi | 1e+05 | 0.048 | 0.076 | 0.699 | 0.694 | **0.644** | 0.25 |

Best multi-τ SOBI-fix null-NMSE = **0.644** (ph=1e+05); the SOBI init never beats the random baseline by more than noise, and even the frozen-known medium (0.076) is unreachable once the medium is re-estimated (true-Z-seed 0.699).


Two independent reasons E4a fails:
1. **SOBI barely separates at d=48.** Validated on *clean* synthetic BSS: quality (mean
   per-source max|corr|) is 0.52 at d=16 but only **0.28–0.42 at d=48** even at T=2048 — the 48
   log-spaced timescales are too finely spaced (adjacent τ differ ~6%), so the joint
   diagonalization is near-degenerate. On the actual buckets SOBI quality is ~0.26–0.31 for
   *both* single- and multi-τ — the multi-timescale physics does not yield a usable rotation
   fix at this d.
2. **The init is not the bottleneck (§11c).** Even a *perfect* rotation (seeding the refinement
   with the true Z) drifts to ≈0.67 — indistinguishable from the random-seed baseline.

Multi-timescale media also *reduce* effective diversity (T_eff≈5.5 at T=128, like single-τ=16),
so the very structure that could enable SOBI trades away the diversity the mechanism needs.

**E4a verdict: ROTATION_FIX_FAILS.**

### 11b. E4b — DL amortized initializer (honest protocol)

A small MLP maps the **blind sufficient statistic** — the bucket covariance
`C = Bcᵀ Bc / T` (R38 §1.18/§2.5, the blind-SIM covariance mechanism) plus the Stage-A
row-space scene — to a scene estimate. Trained on **unlimited** draws from the *declared*
model (fixed P,U; lognormal iid σ_f=0.15 medium; mixed synthetic scene classes: blobs, bars,
smooth texture). Then DL init → **pure-likelihood** refinement (existing tracker, no learned
terms), evaluated on **fresh** draws (never trained on). Mandatory honesty checks:

| check | result | verdict |
|:--|:--:|:--|
| (iii) DL-init **alone**, in-distribution null-NMSE | **0.292** | misses the 0.25 bar even *with* the prior |
| (i) prior-ablation: after pure-likelihood refinement (seeded from the DL scene) | **0.440** | refinement drifts the DL init **worse** ⇒ the likelihood does **not** support the DL null |
| (ii) OOD (checkerboard / text — scene class held out of training) | **1.037** | catastrophic ⇒ the in-dist score was **scene-prior leakage**, not extraction |

All three checks fail. The net can only "recover" null content it has memorised the
distribution of; on unseen scene structure it is worse than the fixed-averaging baseline, and
the likelihood pulls *away* from its guess — the exact prior-imputation failure R38 §2.4 warns
about.

**E4b verdict: DL_INIT_FAILS.**

### 11c. The sharper diagnosis (why any init-based fix fails)

The decisive control, run across iid / single-τ / multi-τ media:

| refinement condition | null-NMSE |
|:--|:--:|
| oracle (medium fully known) | ~0.000 |
| **frozen known medium** (solve x at W=1+U·Z_true, no medium re-estimation) | **≈0.02–0.03** |
| refinement seeded from **TRUE Z**, medium re-estimated (linear or lognormal, reduced or full) | **≈0.67–0.71** |
| refinement from random / SOBI / DL init | ≈0.6–0.9 |

Freezing the true medium recovers the null (≈0.02); **re-estimating the medium from *any*
initialization — including the exact truth — drifts to a wrong-null solution that fits the
buckets as well or better.** This holds for the linear model, the model-matched lognormal
forward, and both the reduced (d×d) and full (T×d) parametrizations. It is the concrete,
empirical form of R38 §1.5 (local≠global; the blind likelihood admits secant/collision
solutions with spurious null content) and of the E1 "SPECROT found a small-residual, wrong-null
solution" result.

**Consequence:** cold-start recovery is *not* an initialization problem, so neither a classical
rotation-fixing init (E4a) nor a learned init (E4b) can rescue it — the blind likelihood
landscape itself does not place truth at (or below) its optimum when the medium is unknown.
The only thing that recovers the null is *knowing the medium* (a reference/calibration channel),
which is exactly the capability the flagship set out to avoid. This strengthens the E1 kill:
under R38's kill tree the direction is a **software-only-blind kill**, with residual value only
in the known-medium-law/oracle regime.

## 12. E5 — covariance-domain resurrection (operator-approved final round)

> **HEADLINE: the direction is RESURRECTED in statistical-identifiability form.** The E4 drift
> was collusion between the scene and the *free medium realizations* `z_t`. Remove `z_t` from the
> estimand entirely and work in the covariance domain (R38 eq 1.18) and the collusion vanishes:
> the estimator is injective, multi-start data-only runs **agree** (no drift), and it **recovers
> the fixed null space** — at a *quantified temporal-sample price*. Verdict: **COVARIANCE_STICKS**.

### 12a. E5a — covariance-domain primary estimator (decisive)

Never estimate `z_t`. With the declared law `K_w = σ_coeff²·U Uᵀ`, the model bucket covariance
factors as `G_P = P diag(x) K_w diag(x) Pᵀ = σ_coeff²·H(x)H(x)ᵀ`, `H(x)=P diag(x)U` — a
**rotation-invariant** object (the medium rotation that colluded in E4 is marginalized, not a free
parameter). The lag-ℓ centered-bucket covariances obey `S_ℓ ≈ r_ℓ G_P` for ℓ≥1 (`r_ℓ=ρ^ℓ`, OU
τ=16 declared) and are **shot-noise-free** (shot is white in time → only `S_0`). Estimator:
`min_{x≥0} ‖H(x)H(x)ᵀ − Ĝ/σ_coeff²‖²_F` from the Stage-A row-space init (Adam, multi-start).

**Three things all hold:**
1. **Injective / recovers the null.** Against the *exact* (noise-free) covariance target, cold data-only
   Adam recovers null-NMSE **0.054** — near-oracle, far below the 0.25 bar. The covariance *does*
   determine the fixed-null content.
2. **Sticks — no collusion.** Across 6 data-only random starts the null-NMSE agrees to std **0.017** (at T=128, clean) — every start lands in the same basin. There is no realization vector to trade against, so the E4 drift is structurally impossible here.
3. **Consistent — recovers the null given enough temporal samples.** The only obstacle at T=128 is
   finite-sample covariance-estimation noise (τ=16 ⇒ ~4 effective independent frames). Null-NMSE
   falls monotonically with T and crosses the bar:

| T (epochs) | 128 | 256 | 512 | 1024 | 2048 | 4096 |
|:--|:--:|:--:|:--:|:--:|:--:|:--:|
| blind null-NMSE | **0.721** | **0.630** | **0.436** | **0.340** | **0.240** | **0.204** |


The blind cell numbers (T=128, the spec) and photon-robustness (lags ℓ≥1 are shot-free, so the
covariance estimate is nearly photon-independent):

Blind cell (T=128, 5-seed medians; multi-start agreement std in parentheses):

| photons/bucket | blind null-NMSE | agreement std |
|:--:|:--:|:--:|
| clean | 0.797 | 0.017 |
| 1e+06 | 0.801 | 0.012 |
| 1e+05 | 0.767 | 0.013 |


Photon-robustness at T=1024 (lags ℓ≥1 are shot-free):

| photons/bucket | blind null-NMSE |
|:--:|:--:|
| clean | 0.340 |
| 1e+06 | 0.327 |
| 1e+05 | 0.343 |
| 1e+04 | 0.322 |


**Statistical-efficiency price.** The pathwise *oracle* (known medium) needs only T≈64 to reach
≈0; the covariance-*blind* estimator needs T≈2048 to cross 0.25 — roughly an order of
magnitude more frames. That is the honest cost of never measuring the medium *and* refusing to
estimate its realizations. But unlike every prior arm it is a **cost, not a wall**: the estimator
is consistent and collusion-free.

**E5a verdict: COVARIANCE_STICKS** — and recovers the fixed null space at a temporal-sample price,
photon-robustly. This is a live, statistical-identifiability resurrection of the "fog as a second
DMD" concept: *repeated bucket measurements with a fixed code bank, under a declared low-dimensional
medium law, recover fixed-null-space scene content with no reference arm and no realization
measurement — via second-order statistics rather than pathwise inversion.*

### 12b. E5b — sparsity-constrained pathwise

Sparse scene (few bright squares; T=128, single-τ OU). §11c-style table, null-NMSE:

| scene | oracle | frozen-known | true-Z-seed (no L1) | true-Z-seed (+L1) | baseline (+L1) |
|:--|:--:|:--:|:--:|:--:|:--:|
| dense (control) | 0.000 | 0.028 | 0.863 | 0.709 | 0.766 |
| **sparse** | 0.000 | 0.013 | 0.539 | **0.972** | 0.756 |

**SPARSITY_BREAKS_COLLUSION: NO** — with the L1/sparsity constraint the true-Z-seeded pathwise refinement goes 0.539 → 0.972 (no-L1 → +L1). The sparsity prior does not by itself restore truth as the pathwise minimum here.

### 12c. E5c — anchor frames

_E5c (anchor frames) not run — E5a already resurrects the direction; deprioritized per instruction._

### 12d. E5d — declared-law mismatch (R38 F5 kill test)

The sharpest remaining threat to E5a is F5: *if an exact simulator-matched declared law is
load-bearing, it is a simulation artifact and must be killed.* Data are generated with the TRUE
law; the estimator is fed a DECLARED (possibly wrong) law. Cell M=96, d=48, σ_f=0.15; 1e5 photons;
5-seed medians; degradation is relative to the matched baseline at the same T.

**2048 epochs** — matched baseline null-NMSE = **0.310**:

| mismatch axis | condition | null-NMSE | degradation |
|:--|:--|:--:|:--:|
| basis rotation (estimator U rotated; simulator U_true) | rot~5% (∠5%) | 0.316 | +2% |
|  | rot~10% (∠10%) | 0.312 | +1% |
|  | rot~20% (∠20%) | 0.318 | +3% |
| wrong τ (estimator assumes 16) | true_tau=8 | 0.285 | -8% |
|  | true_tau=32 | 0.377 | +22% |
| wrong σ_f (estimator assumes 0.15 → pure scale/gauge) | true_sf=0.1 | 0.313 | +1% |
|  | true_sf=0.2 | 0.331 | +7% |
| subspace dim (estimator d=48) | true_d=40 | 0.331 | +7% |
|  | true_d=56 | 0.320 | +3% |

**1024 epochs** — matched baseline null-NMSE = **0.376**:

| mismatch axis | condition | null-NMSE | degradation |
|:--|:--|:--:|:--:|
| basis rotation (estimator U rotated; simulator U_true) | rot~5% (∠5%) | 0.375 | -0% |
|  | rot~10% (∠10%) | 0.368 | -2% |
|  | rot~20% (∠20%) | 0.375 | -0% |
| wrong τ (estimator assumes 16) | true_tau=8 | 0.342 | -9% |
|  | true_tau=32 | 0.433 | +15% |

**E5d verdict: MISMATCH_ROBUST.** In the F5 key cells @T=2048 the degradation is rot-10% +1%, τ-mismatch (true 8 / 32) -8% / +22% (worst 22% < 25%). 
The estimator does **not** require the exact simulator law: τ and σ_f mismatches are absorbed by the scalar gauge (every lag stays ∝ G_P, so a wrong τ/σ_f only rescales the target), and moderate basis rotation degrades gracefully. **E5a graduates from signal to candidate.** The load-bearing assumption is the *subspace* (U's span and its dimension), not its exact realization — consistent with a declared-medium-class model, not an exact-simulator artifact.

Two honest caveats: (i) the τ=32 cell's +22% is **partly reduced medium diversity, not pure
declared-law mismatch** — a true τ=32 halves the effective independent-frame count vs τ=16
(T_eff ∝ (1−φ)/(1+φ)), so it degrades even a *matched* estimator; the estimator's wrong-τ
assumption *alone* is gauge-absorbed. (ii) The run was **externally terminated** (another
concurrent session cleared the GPU's python processes, exit 127) during the T=1024 block; the
**T=2048 block — which carries the verdict — completed in full**, and the missing T=1024 σ_f/dim
cells are redundant with their T=2048 counterparts. Numbers reconstructed from the run log into
`E5d_results.json`.

## 13. E7 — sample-complexity & optimal design (pre-R39 anchor)

The T≈2048 cost of the moment-matching estimator (E5a) is root-caused to the number of
**independent medium realizations** (T_eff ≈ T(1−φ)/(1+φ) ≈ 64 at the bar-crossing), which
design tricks cannot mint. E7 attacks the two levers that touch the *right* resource:
statistical efficiency of the estimator (E7a) and the exposure allocation (E7b).

### 13a. E7a — statistical efficiency: marginal-likelihood MLE

Because the centered model is linear-Gaussian state space (`z_t` OU, `b_t = H(x)z_t + noise`,
`H=P diag(x)U`), the Kalman recursion gives the **z-marginalized** likelihood exactly. We
maximize it in x by **EM** (E-step = RTS smoother returning E[z] *and* E[zzᵀ]; M-step =
closed-form nonneg x-solve). Crucially only the *posterior moments* enter — no `z` point
estimate — so, like E5a, there is nothing to collude with. Two diagnostics confirm the estimator
is sound (mini cell, T=512): from the **true** x, EM reaches null-NMSE **0.002** (the marginal
likelihood's global max is at truth, nll monotonically ascended); from Stage-A it sticks in a bad
local max (0.74). So the efficient recipe is **two-stage: moment-matching (E5a) for the global
basin → MLE/EM for efficiency** — at T=512 mini this sharpens 0.189 (moment) → **0.057** (MLE).

**Mini (16×16, M=48, d=24)** — null-NMSE, 1e5 photons (wall-clock/fit in parentheses):

| T | moment-matching (E5a) | **two-stage MLE** | t_moment | t_MLE |
|:--:|:--:|:--:|:--:|:--:|
| 128 | 0.348 | **0.098** | 11.8s | 3.4s |
| 256 | 0.252 | **0.062** | 11.8s | 6.8s |
| 512 | 0.174 | **0.035** | 11.8s | 13.3s |
| 1024 | 0.129 | **0.014** | 11.8s | 27.3s |

Crossing 0.25: moment at T≈512, **MLE at T≈128** (**4.0× reduction**).

**Full (32×32, M=96, d=48)** — null-NMSE, 1e5 photons (wall-clock/fit in parentheses):

| T | moment-matching (E5a) | **two-stage MLE** | t_moment | t_MLE |
|:--:|:--:|:--:|:--:|:--:|
| 512 | 0.467 | **0.149** | 12.3s | 16.0s |
| 1024 | 0.408 | **0.074** | 11.4s | 31.7s |

Crossing 0.25: moment at T≈None, **MLE at T≈512**.


### 13b. E7b — exposure economics

At **fixed total exposures** (m·T constant), does subsampling patterns per epoch — fewer patterns
m, proportionally more epochs T (more, cheaper epochs, hence more T_eff) — beat the full bank?
Constraint: blind capacity needs m > d, so m ≤ d is degenerate.

**Mini d=24**, fixed total exposures m·T = 196608:

| m (patterns/epoch) | T (epochs) | χ_full | MLE null-NMSE |
|:--:|:--:|:--:|:--:|
| 24 | 8192 | 0.00 | degenerate (m≤d) |
| 36 | 5461 | 1.17 | 0.018 |
| 48 | 4096 | 2.34 | 0.015 |

Exposure-optimal viable point: **m=48, T=4096** (null-NMSE 0.015). 
At fixed exposure, **m is roughly neutral** (36→48: 0.018→0.015). 
Below m=d the estimator is degenerate (no blind capacity), so subsampling cannot go arbitrarily cheap — the medium dimension d floors patterns-per-epoch.


**E7 verdict.** The marginal-likelihood MLE (moment-matching init → EM) is the efficiency win: on the mini cell it crosses 0.25 at T≈128 vs T≈512 for moment-matching (**4.0× sample-complexity reduction**); on the 32×32 cell the MLE crosses at T≈512. Exposure economics: patterns-per-epoch is floored by the medium dimension (m>d), so the exposure budget is best spent on epochs at m just above d, not on the full bank if d≪M. These numbers anchor the sample-complexity (E7a) and optimal-design (E7b) questions for R39.

## 14. E8 — beyond-modulator-band super-resolution (F4-immune)

The strongest form of the claim. The pattern bank is **band-limited to the DMD pixel limit**
(spatial frequencies fx,fy ≤ 3), but the medium's fine speckle **exceeds** that band
(1 ≤ i+j ≤ 6), and the scene carries deliberate content *outside* the pattern band. A
band-limited operator **cannot leave its band**, so the fresh-known-pattern comparator (the R38
F4 threat) is **physically walled at rel-err 1.000** on beyond-band content — no software trick
lets fresh band-limited patterns recover what they cannot address. The fixed-bank *fluctuation*
route can, because the unmeasured fine medium mixes beyond-band scene content down into the
measurable band. Cell: 16×16, band-limited patterns M=24, σ_f=0.30, τ=8, 1e5 photons; metric =
relative error on the 8 beyond-band scene frequencies (after scalar gauge); 3 seeds.

Beyond-band recovery — relative error on out-of-pattern-band scene content (1.000 = recovers nothing; lower is better):

| T (epochs) | fresh band-limited patterns (physics wall) | fixed-bank moment | **fixed-bank MLE** | oracle (medium known) |
|:--:|:--:|:--:|:--:|:--:|
| 512 | 1.000 | 0.650 | **0.505** | 0.389 |
| 1024 | 1.000 | 0.547 | **0.439** | 0.221 |
| 2048 | 1.000 | 0.480 | **0.481** | 0.091 |

**E8 verdict: super-resolution is REAL and F4-immune, but blind recovery is QUALITATIVE ONLY (best blind 0.44 > 0.3).**

- **Physics wall (F4-immune):** fresh band-limited patterns are pinned at **1.000** at every T — a band-limited operator cannot address beyond-band content, so no fresh-pattern or software route can compete. This is the comparator-proof core of the claim.
- **Physical headroom (oracle):** with the medium known, beyond-band error falls to **0.091** at T=2048 — the fine speckle genuinely mixes beyond-band scene content into the measurable band, and it is physically recoverable.
- **Blind gap:** the fixed-bank routes recover beyond-band content (moment 0.48, MLE best 0.44) but do **not** reach ≤0.3, and here the MLE does **not** beat moment (unlike E7). Mechanism: a pattern bank band-limited to fx,fy≤3 has **fluctuation rank (PB+1)²−1 = 15**, independent of how many patterns M are stored — while the medium band is **d=27 modes**. So the beyond-band geometry is inherently **d > M_eff**, R38's hardest profiling regime, where blind medium estimation is weakly determined and the MLE cannot add efficiency. Increasing M does not escape it (rank-capped by the band limit).

**For R39:** the *qualitative* super-resolution claim — fog reveals content the modulator physically cannot address — is airtight (1.000 physics wall + 0.09 oracle ceiling). The *quantitative* blind number is not yet ≤0.3 in this native d>M_eff geometry; closing the oracle→blind gap (0.09 → 0.44) is the open problem, not the physics.

## 10. Deliverables in this directory

- `FOG_DMD_PROBE_REPORT.md` (this file)
- `E1_results.json` (controls + primary + linear-reference + OU-sensitivity, 5-seed)
- `E2_noise_curve.json`, `E2_noise_curve.png`
- `E4a_results.json` (single/multi-τ × 3 photons, SOBI + attribution), `E4b_results.json`
- `E5_results.json` (covariance estimator: cell + T-scaling + photon-robustness + injectivity),
  `E5b_results.json` (sparsity-constrained pathwise)
- Solvers: `fog_common.py`, `fog_tracker.py` (incl. `cov_estimate`), `fog_sobi.py` (joint-diag)
- Runners: `run_grid.py` (E1), `run_e2.py` (E2), `run_e4a.py` (E4a), `run_e4b.py` (E4b),
  `run_e5.py` (E5a), `run_e5b.py` (E5b), `fill_report.py` / `fill_e4a.py` / `fill_e5.py`
- Diagnostics: `sanity_reproduce.py`, `exp_identify.py`, `exp_coldstart.py`, `exp_specrot.py`,
  `exp_diversity.py`, `exp_m96_diag.py`, `exp_ceiling.py`, `smoke96.py`
- No git commit (per instruction; coordinator reviews first).
