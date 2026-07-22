# PROBE A — RQL efficiency vs the jitter-capped information ceiling

**ROUND63-NEXT software-method hunt · Probe A (efficiency probe)**
Repo `D:\GI_another` · nu∈{2000,200}, c=0.05, loads ρ∈{1, 5.7, 22.25} · 2026-07-22 (Fable)

---

## VERDICT — THIS LANE IS DEAD (at the deployed operating point)

At the deployed terminal dwell **ν = 2000**, the frozen RQL reconstruction is
**information-efficient against the jitter-capped ceiling J(ρ;c=0.05)**:

- **Scalar (data term, per-frame CRB):** RQL is **98.9–99.5 % variance-efficient**
  and **≥ 95.9 % MSE-efficient** at every load. Certified data-term headroom on the
  range component is **< 0.05 dB (variance) / < 0.19 dB (MSE)**.
- **Image (range subspace):** in the data-constrained directions the deployed
  RQL+TV reconstruction sits **at the unbiased Fisher floor** (per-direction
  efficiency 0.96–1.38; max "actual-above-floor" headroom **≤ 0.20 dB**). In the
  ill-conditioned directions the TV prior already **beats** the unbiased CRB by
  100–400×. Null-space error is 1–10 % of total MSE and is *not* Fisher-bounded.

There is **no certified software headroom of a size worth a positive-result
campaign** on the range component at ν = 2000. A better data-fit estimator (or "DL
as tool" racing the certified ceiling) cannot recover more than ≈ 0.2 dB.

The **only** measurable inefficiency is at the shorter certification dwell
**ν = 200 with high load (ρ = 22.25): ≈ 0.62 dB (variance) / ≈ 1.7 dB (MSE)** — and
it is present **already at c = 0** (2.3 dB MSE), so it is the **finite-sample bias
of the closed-form dead-time inversion `λ̂ = N/(T−Nτ)`**, not a jitter/DL
phenomenon. It is recoverable by a one-line **bias-corrected / exact-MLE scalar
inversion**, not by a learned reconstruction.

---

## 1. Question

Chapter 1 proved a per-frame information ceiling J(ρ;c) (jitter-capped ridge law).
The frozen reconstruction in all campaigns is **RQL** (renewal quasi-likelihood +
isotropic TV, truth-free λ_TV). **How close does RQL get to the ceiling at the
deployed operating points?** If a large efficiency gap exists on the
*data-constrained (range) component*, a better estimator / DL has certified
headroom → a positive-result candidate. If RQL is near-efficient, the lane is dead.

## 2. Method

**Forward model** (`bridge_harness.simulate_counts_jitter`, identical draws to the
M1 campaign): active-start non-paralyzable renewal counting with i.i.d. lognormal
holds B ~ τ·LogNormal(mean τ, CV = c). The count law depends only on (ρ, ν, c);
τ cancels. Deployed jitter c = 0.05.

**RQL rate estimate** (scalar closed form; minimizer of the deployed per-frame RQL
objective Q(λ;N)=(T−Nτ)λ − N·log λ, identical to `physics.precorrect_rates`):
`λ̂ = N / max(T − Nτ, floor)`.

**Ceiling J(ρ;c)** = per-slot Fisher information about **log-rate**; per-frame
information = ν·J. Estimated with the FROZEN, validated score machinery
(`jitter_score_diag_colab_frozen.py`, imported verbatim): `fi_mi` (score-identity,
primary) and `fi_cs` (independent cross-check). The image J(ρ) lookup uses the
robust histogram-score `fi_hist` (positive & stable at high load, where the MI
identity's λ²·Var term needs millions of frames — this is the estimator whose
peak-extraction R16 adjudicated; here it is used only for the **value** of J at
fixed ρ, where it is well-behaved; fi_mi/fi_cs/fi_hist agree to ~1–2 % at N_MC=1e5).

**Scalar efficiency** (log-rate θ = log λ, the parameterization of J):
`eff = CRB_θ / Var(θ̂)`, `CRB_θ = 1/(ν·J)`, N_MC = 150 000 frames/point.

**Image Fisher floor (range component only).** Each row i measures the linear
functional ρ_i = a_i·x with per-frame info I_ρ(ρ_i) = ν·J(ρ_i;c)/ρ_i². The
information matrix is **I_x = Aᵀ diag(I_ρ) A** (design A = 52 pre-scan + 972 main
rows = the deployed `bridge_harness.reconstruct` path, scaled so mean main load =
ρ). Eigendecompose I_x = V diag(μ) Vᵀ; the **range** = top-r eigenvectors
(r = numerical rank of A), **null** = the rest. **floor = Σ_range 1/μ_k** is the
CRB trace bound on the range-space MSE of any *unbiased* estimator. Errors of the
deployed RQL+TV reconstruction (8 seeds) are decomposed in the **same** eigenbasis:
range_mse = Σ_range (v_k·e)², null_mse = ‖e‖² − range_mse. The **null (prior)
part is explicitly NOT bounded**. Efficiency_range = floor / range_mse_actual.

## 3. Pipeline validation (c = 0, against the exact analytic Fisher)

`fi_mi` reproduces `physics.exact_fisher_analytic` (exact deterministic-dead-time
per-frame log-rate information) to **< 0.02 %** at ν = 2000 and **< 0.12 %** at
ν = 200 — end-to-end confirmation of the ceiling estimator (consistent with the
R16 anchor ρ* = 22.2543 validated to < 0.1 %).

| ν | ρ | J_mi (MC) | J_exact | rel. diff |
|---|---|-----------|---------|-----------|
| 2000 | 1.00 | 0.500011 | 0.500010 | 0.000 % |
| 2000 | 5.70 | 0.849270 | 0.849264 | 0.001 % |
| 2000 | 22.25 | 0.935473 | 0.935449 | 0.003 % |
| 200 | 22.25 | 0.779882 | 0.779019 | 0.111 % |

## 4. Scalar efficiency (data-term, the certified core)

CRB and Var are per-frame, log-rate. `head` = 10·log₁₀(1/eff) dB (positive = RQL
above the bound = recoverable headroom).

| ν | ρ | c | J(ρ;c) | eff_var | eff_MSE | bias_θ | head_var | head_MSE |
|---|---|---|--------|---------|---------|--------|----------|----------|
| **2000** | 1.00 | 0.05 | 0.4988 | **0.9946** | 0.9945 | +0.0002 | +0.024 dB | +0.024 dB |
| **2000** | 5.70 | 0.05 | 0.7855 | **0.9889** | 0.9853 | +0.0015 | +0.049 dB | +0.064 dB |
| **2000** | 22.25 | 0.05 | 0.4228 | **0.9893** | 0.9593 | +0.0061 | +0.047 dB | +0.180 dB |
| 2000 | 1.00 | 0.00 | 0.5000 | 1.0013 | 1.0012 | +0.0003 | −0.006 dB | −0.005 dB |
| 2000 | 5.70 | 0.00 | 0.8493 | 1.0007 | 0.9961 | +0.0016 | −0.003 dB | +0.017 dB |
| 2000 | 22.25 | 0.00 | 0.9355 | 0.9897 | 0.9303 | +0.0059 | +0.045 dB | +0.314 dB |
| 200 | 1.00 | 0.05 | 0.4989 | 0.9883 | 0.9877 | +0.0026 | +0.051 dB | +0.054 dB |
| 200 | 5.70 | 0.05 | 0.7752 | 0.9733 | 0.9333 | +0.0169 | +0.118 dB | +0.300 dB |
| **200** | **22.25** | 0.05 | 0.3904 | **0.8665** | **0.6785** | +0.0640 | **+0.622 dB** | **+1.684 dB** |
| 200 | 22.25 | 0.00 | 0.7799 | 0.8922 | 0.5909 | +0.0605 | +0.495 dB | +2.285 dB |

**Reading.**
- The jitter **cap is real**: J at ρ = 22.25 falls from 0.935 (c = 0) to 0.423
  (c = 0.05) — jitter roughly halves the ceiling. RQL's efficiency is measured
  *against this lowered ceiling* and still holds at ~0.99.
- **At ν = 2000 the data term is efficient at every load** (eff_var 0.989–0.995).
  Headroom < 0.05 dB (variance), < 0.19 dB (MSE, worst at ρ = 22.25 where the
  closed form is mildly biased).
- The **ν = 200 / ρ = 22.25** corner is the only real gap (0.62 dB var / 1.68 dB
  MSE). It is **bias-dominated** (bias_θ = +0.064) and **exists at c = 0** (2.29 dB
  MSE), i.e. it is the low-count nonlinearity of `N/(T−Nτ)` (T−Nτ ≈ 8 of 200), not
  jitter. A bias correction / exact renewal MLE closes it — a trivial, non-DL fix.

## 5. Image range-space efficiency (ν = 2000, c = 0.05, 8 seeds)

Design rank = **1012** of 1024 (972 main rows are full-rank; +52 pre-scan → 12-dim
null space). Null-space error is small; **most error is in the range space**.

| scene | ρ | PSNR_rad | total MSE | range MSE | null MSE | null % | agg. floor/actual |
|-------|---|----------|-----------|-----------|----------|--------|-------------------|
| contour_1 | 1.00 | 19.79 | 3.52e-8 | 3.47e-8 | 4.5e-10 | 1.3 % | 206× |
| contour_1 | 5.70 | 20.68 | 2.87e-8 | 2.84e-8 | 3.2e-10 | 1.1 % | 217× |
| contour_1 | 22.25 | 19.01 | 4.21e-8 | 4.16e-8 | 5.2e-10 | 1.2 % | 190× |
| microtex_1 | 1.00 | 24.37 | 7.13e-8 | 6.91e-8 | 2.2e-9 | 3.1 % | 102× |
| microtex_1 | 5.70 | 24.91 | 6.28e-8 | 6.08e-8 | 2.0e-9 | 3.3 % | 99× |
| microtex_1 | 22.25 | 23.93 | 7.88e-8 | 7.67e-8 | 2.1e-9 | 2.6 % | 102× |
| control_2 | 1.00 | 23.41 | 7.98e-8 | 7.25e-8 | 7.3e-9 | 9.2 % | 99× |
| control_2 | 5.70 | 24.17 | 6.69e-8 | 6.00e-8 | 7.0e-9 | 10.4 % | 104× |
| control_2 | 22.25 | 21.99 | 1.11e-7 | 1.03e-7 | 7.4e-9 | 6.7 % | 83× |

The aggregate floor/actual is **83–217× (i.e. actual ≪ floor)**: the deployed
reconstruction is already **far inside** the unbiased CRB. This aggregate is
dominated by the ill-conditioned tail (the design's near-null directions), where
the unbiased CRB diverges — an artifact of `trace(I_x⁺)`, not a clean per-direction
statement. The **eigenband breakdown resolves it**: the 1012 range directions in
descending-μ quartiles, efficiency = floor/actual per band.

| scene | ρ | band 1 (best) | band 2 | band 3 | band 4 (worst) |
|-------|---|---------------|--------|--------|----------------|
| contour_1 | 5.70 | 1.19 | 1.31 | 2.06 | 417 |
| contour_1 | 22.25 | 1.19 | 1.38 | 2.33 | 377 |
| microtex_1 | 5.70 | **0.97** | 1.06 | 1.26 | 142 |
| microtex_1 | 22.25 | 1.02 | 1.15 | 1.52 | 168 |
| control_2 | 5.70 | 1.05 | **0.97** | 1.42 | 146 |
| control_2 | 22.25 | 1.00 | **0.96** | 1.21 | 125 |

**Reading.**
- **Well-measured directions (bands 1–2, the data-constrained half):** efficiency
  0.96–1.38. The deployed RQL+TV range error is **at the unbiased Fisher floor**.
  The worst "actual-above-floor" case is eff = 0.956 → **+0.20 dB** headroom;
  most bands are ≥ 1 (already below the unbiased bound). **No certified headroom.**
- **Ill-conditioned directions (band 4):** efficiency 125–417. The TV prior beats
  the unbiased CRB by two orders of magnitude; an unbiased estimator would be
  100–400× *worse* here. This is the null-adjacent subspace — **not** certified
  headroom for a better data-term estimator.

## 6. Where the gap lives

- **Data-term inefficiency: ≈ 0.** Scalar eff_var 0.989–0.995 at ν = 2000; image
  bands 1–2 at the floor. RQL's renewal quasi-score is the (asymptotically)
  efficient score for the count under the ceiling model, and the numbers confirm it.
- **TV-prior effect: helps, not headroom.** In the ill-conditioned directions the
  prior drives the deployed error 100–400× below the unbiased CRB. That residual
  (band 4 holds the largest *absolute* range error) could in principle be reduced
  further by a **different prior** (e.g. a learned/DL prior) — but that is the
  **null-space/prior lane, which Probe A explicitly does not count**, and the
  Fisher ceiling does **not** certify it (it says the opposite: unbiased error is
  larger there). Its benefit is scene-dependent and uncertified (TV is already
  near-optimal on the smooth `control` scene; a learned prior might help `microtex`,
  but speculatively).
- **Finite-sample bias:** the sole certified gap (ν = 200, high load) is the
  closed-form inversion's low-count bias — a cheap deterministic fix, not a lane.

## 7. Verdict

**Certified software headroom on the RANGE component ≈ 0.2 dB — the lane is DEAD**
at the deployed terminal operating point (ν = 2000, c = 0.05, ρ up to the ridge
22.25). Both independent measurements agree: the scalar per-frame data term is
98.9–99.5 % efficient (< 0.05 dB variance / < 0.19 dB MSE), and the image
range-space error is at the unbiased Fisher floor in the data-constrained
directions (≤ 0.20 dB). A better reconstruction *estimator* — including "DL as a
tool racing the certified ceiling" — has no certified room to run at the deployed
point.

Two honest footnotes, neither a positive-result campaign:
1. **ν = 200 / high load** carries ≈ 1.7 dB MSE of recoverable inefficiency, but it
   is the `N/(T−Nτ)` closed form's finite-sample bias (present at c = 0), closed by
   a bias-corrected / exact-renewal-MLE scalar inversion — deterministic, one-off.
2. **Ill-conditioned (band-4) residual** is the largest absolute range error and is
   a *prior*-improvement target (uncertified, scene-dependent), not a data-term /
   certified-ceiling target. If a future line wants DL, it would be a **learned
   prior for near-null directions**, not a better data-fit — and it should be
   chartered against a *design/prior* baseline, not the Fisher ceiling, which this
   probe shows RQL already saturates.

## 8. Provenance & honesty

- **Scripts** (this dir): `probe_a_scalar.py`, `probe_a_image.py`,
  `probe_a_image_bands.py`. Logs: `scalar_run.log`, `image_run.log`,
  `bands_run.log`. Data: `scalar_efficiency.json`, `image_efficiency.json`,
  `image_eigenbands.json`, `J_grid_c05_nu2000.npz`.
- **Read-only** on all existing results; ceiling estimators imported verbatim from
  the frozen `jitter_score_diag_colab_frozen.py`; reconstruction is the deployed
  `bridge_harness.reconstruct` (RQL + isotropic TV, truth-free λ_TV) via
  `solvers.run_arm("RQL", …)`. No RQL tuning; no load cherry-picking (all three
  loads and both ν reported, including the one corner where RQL underperforms).
- **MC budget:** scalar N_MC = 150 k/point; J-grid N_MC = 100 k/point (47 loads);
  image 8 seeds/cell (bands 6 seeds). Runtime ≈ 45 min CPU. The score-identity MI
  estimator is unstable at high ρ below ~10⁶ frames (its λ²·Var term); the image
  J-lookup therefore uses `fi_hist`, cross-checked against fi_mi/fi_cs (agree
  ~1–2 %) and against the exact analytic Fisher at c = 0.
- **Bound scope:** the floor bounds only the range (measured) subspace under the
  unbiased CRB; the 12-dim null space and the biased-estimator regime in
  ill-conditioned directions are stated as *not* Fisher-bounded.
