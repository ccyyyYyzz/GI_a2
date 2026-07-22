# R29 ruling (GitHub issue #21, raw)

Title: R29 ruling: software-method triage and drift-orthogonal scheduling campaign
Posted: 2026-07-22T09:50:27Z

---

# R29 — Independent software-method generation and campaign ruling

**Reference commit:** `61c7b98`  
**Source:** `docs/ROUND63_GPT_ROUND29_QUESTION.md`  
**Binding negative:** `results/round63_bridge/BRIDGE_VERDICT.md`

**Scope:** software-only methods. No detector modification, no new optical component, no hold-off intervention, no timing telemetry, and no hardware-side escape hatch. The frozen M1/bridge results remain unchanged.

> **R29 VERDICT.** The software space is not completely exhausted, but it is much narrower than it first appears. The bridge result rules out the most tempting lane—scene-adaptive Fisher/OED geometry driven by a noisy pre-scan—and the count-only identity rules out recovering physically discarded dead-time information in reconstruction.
>
> The only candidate presently strong enough to justify a confirmatory campaign is **drift-orthogonal pattern scheduling**: reorder a fixed, already-accepted pattern set so that its object score is orthogonal to the low-dimensional temporal gain nuisance. It is scene-independent, costs no photons, does not depend on a noisy router, and attacks the Peng–Chen failure mechanism directly.
>
> This independently generated candidate is the same scientific lane as the operator’s **S1**. That convergence increases confidence; it is not a reason to invent another layer.
>
> **There is no unconditional campaign GO today.** The running image-level S1 probe is the gate. If that probe fails the modest DEV bar specified below, I would rule that **no current software candidate clears the operator’s bar**. S2 and S3 are not substitute flagship methods. The only reserve worth a one-day micro-probe is the exact conditional-packet likelihood in §2.1.

---

# 1. What the bridge negative changes

The frozen bridge result is unusually informative:

- the finite-library image oracle gained only `+0.68 dB` median;
- every winning bank was the knob/scattered family `L0/L1`;
- structured geometry never won;
- the maximin pre-scan allocator lost `−8.3 dB` relative to the library oracle;
- even the true-scene, dose-relaxed Fisher/A-risk FW oracle was `−12.9 dB` below the simple library at image level.

The next method should therefore have all four properties:

1. **scene-independent or nearly so**—no fragile route from a low-SNR pre-scan;
2. **invariant/orthogonal rather than “more optimal” in a surrogate**;
3. **no extra physical resource unless explicitly charged**;
4. **the same reconstruction across compared schedules**, so an image gain cannot be attributed to a hidden estimator change.

This points toward nuisance cancellation and score orthogonalization, not another adaptive geometry search.

---

# 2. Independent candidate generation

The seven candidates below were generated from the frozen identities and failure mechanisms before merging them with S1–S3.

## 2.1 Conditional Packet Likelihood GI (`CPL-GI`)

### Mechanism

Partition the acquisition into short packets in which the dynamic-medium gain is approximately common. For packet `b`, patterns `m_{bj}` and photon counts satisfy

\[
Y_{bj}\mid a_b,x\sim\operatorname{Poisson}\!\left(a_b\phi\,m_{bj}^{\top}x\right),
\qquad j=1,\ldots,q.
\]

Conditioning on the packet total

\[
N_b=\sum_jY_{bj}
\]

eliminates the unknown multiplicative gain exactly:

\[
\boxed{
Y_b\mid N_b,x
\sim
\operatorname{Multinomial}\!\left(
N_b,
\frac{m_{b1}^{\top}x}{\sum_jm_{bj}^{\top}x},\ldots,
\frac{m_{bq}^{\top}x}{\sum_jm_{bj}^{\top}x}
\right).
}
\]

If the packet is constant-sum,

\[
\sum_jm_{bj}=c_b\mathbf 1,
\]

and image flux is normalized or separately anchored, the denominator is fixed. The conditional negative log-likelihood

\[
-\sum_{b,j}Y_{bj}\log(m_{bj}^{\top}x)
\]

is convex on `x>=0`, and a TV or other frozen prior can be added. Complementary pairs are the `q=2` special case.

### Why a positive is plausible

The method does not estimate the gain and then divide by a noisy estimate; it removes the block gain from the likelihood. It directly targets the hidden-gain covariance term in O1.

### Closest prior art

- Conditional Poisson/multinomial inference is classical; see [Palmgren, *Biometrika* 1981](https://doi.org/10.1093/biomet/68.2.563) and [Sei, *Information Geometry* 2024](https://doi.org/10.1007/s41884-022-00082-w).
- Adjacent complementary normalization for illumination fluctuations already exists in SPI: [“Illumination Temporal Fluctuation Suppression for Single-Pixel Imaging,” *Sensors* 2023](https://doi.org/10.3390/s23031478).
- Periodic assay-pattern correction is already established: [“Noise reduction in computational ghost imaging by interpolated monitoring,” *Optics Express* 2018](https://pubmed.ncbi.nlm.nih.gov/30118039/).
- The Chen line already contains temporal and learned corrections: [Xiao–Zhou–Chen 2022](https://doi.org/10.1364/OL.463897) and [Peng–Chen 2023](https://doi.org/10.1364/OL.499787).

The defensible residual is **not normalization**. It is an exact conditional likelihood for arbitrary constant-sum packets, with a coherence-matched packet construction and a convex image estimator.

### Kill risk

- gain varies materially inside the packet;
- dead-time/non-Poisson counting destroys the exact multinomial factorization;
- packet totals are too small, so conditioning throws away useful scale information;
- absolute radiometry remains unidentifiable without a stationarity/flux anchor;
- reviewers may regard the method as a statistically dressed version of differential Hadamard normalization.

**Ruling:** strong reserve candidate; not the first campaign because of the prior-art collision and block-gain assumption.

---

## 2.2 Nuisance-Orthogonal Score GI (`NOS-GI`)

### Mechanism

Let `eta` parameterize the gain path or a low-rank drift basis. From the joint model form the efficient object score

\[
\boxed{
S_x^{\perp}
=
S_x-I_{x\eta}I_{\eta\eta}^{-1}S_{\eta}.
}
\]

Estimate `eta` by a low-dimensional smoother, preferably cross-fitted, and reconstruct `x` from the penalized orthogonal estimating equation. The expected derivative of `S_x^perp` with respect to nuisance error vanishes at the model, so first-stage gain-estimation errors enter only at second order.

### Why a positive is plausible

A large fraction of the failure in dynamic scattering is leakage of gain drift into the object score. The method removes that leakage at estimation time rather than asking a neural network to learn it implicitly. When the schedule satisfies O4-A, `I_xeta=0` and the method automatically collapses to the ordinary reconstruction.

### Closest prior art

Efficient-score projection and Neyman orthogonality are classical semiparametric statistics; blind gain calibration is also mature. Relevant neighbors include [Ling–Strohmer, SIIMS 2018](https://doi.org/10.1137/16M1103634), [Cambareri–Jacques 2019](https://doi.org/10.1093/imaiai/iay004), and general higher-order orthogonalization such as [Bonhomme–Jochmans–Weidner 2024](https://arxiv.org/abs/2412.10304). Joint calibration and imaging is also standard in radio interferometry; e.g. [Roth et al., A&A 2023](https://doi.org/10.1051/0004-6361/202346851).

The residual novelty would be the exact dynamic-GI score, its O1 information ceiling, and the reduction to ordinary reconstruction under an O4-A schedule.

### Kill risk

- regularization bias dominates semiparametric efficiency;
- the gain model/basis is misspecified;
- `I_etaeta` is ill-conditioned;
- the efficient score amplifies noise even while removing bias;
- implementation becomes a joint nonconvex reconstruction under another name.

**Ruling:** mathematically strong, but too complex to be the first software-positive attempt.

---

## 2.3 Drift-Orthogonal Pattern Scheduling (`DOPS-GI`)

### Mechanism

Keep the pattern multiset, pattern count, dwell, dose, and reconstruction fixed. Change only the temporal permutation. For a declared gain basis/covariance `H`, choose a permutation `pi` minimizing the nuisance cross-information or the equivalent Schur-complement loss, for example

\[
\boxed{
\mathcal L(\pi)
=
\left\|
I_{xx}^{-1/2}
M_{\pi}^{\top}R^{-1}D_{\bar s}H
\left(H^{\top}D_{\bar s}R^{-1}D_{\bar s}H+\Lambda\right)^{-1/2}
\right\|_F^2.
}
\]

For balanced patterns, `D_bar s` can be frozen from the scene-independent carrier. Complementary patterns remain adjacent; pair blocks are then ordered to cancel low-frequency drift moments.

### Why a positive is plausible

It attacks the exact O4-A loss term. The repository already reports violation metrics of approximately `1770` for time order, `63.6` for random order, and `0` for the paired construction. Unlike RLMI, it has:

- no noisy pre-scan;
- no scene-specific routing;
- no Fisher-selected geometry;
- no additional photons;
- no change to the inverse operator as a set—only its relation to time.

### Closest prior art

The abstract idea is not new:

- trend-free run orders are classical: [Jacroux–Saha Ray, *Biometrika* 1990](https://doi.org/10.1093/biomet/77.1.187), [Atkinson–Donev, *Technometrics* 1996](https://doi.org/10.1080/00401706.1996.10484545), and [Tack–Vandebroek 2002](https://doi.org/10.1016/S0167-9473(02)00056-7);
- Hadamard ordering for SPI is established, though usually for early image quality rather than drift: [López-García et al., *Optics Express* 2022](https://doi.org/10.1364/OE.451656);
- adjacent complementary normalization already suppresses illumination fluctuations: [Sensors 2023](https://doi.org/10.3390/s23031478);
- Chen’s temporal correction addresses the same physical setting after acquisition: [Optics Letters 2022](https://doi.org/10.1364/OL.463897).

The narrow novelty must be:

> a multiplicative-gain, task-information schedule for programmable GI that minimizes the exact object–gain cross-information and is validated at image level without extra measurements.

### Kill risk

The bridge warning still applies: a huge Fisher-side improvement may yield little PSNR gain if the reconstruction is prior-limited or if pair differencing raises shot noise. The assumed drift basis/correlation time may also be wrong.

**Ruling:** best campaign candidate. This is the operator’s S1 lane.

---

## 2.4 Pilot-Interleaved State-Space Normalization (`PIS-GI`)

### Mechanism

Insert sparse all-ones or assay patterns, fit an OU/GP/Kalman gain path, interpolate it, normalize bucket measurements, and reconstruct with the existing method.

### Why a positive is plausible

It directly observes the nuisance carrier and usually works when drift is slow.

### Closest prior art

This is very close to interpolated monitoring CGI ([Optics Express 2018](https://pubmed.ncbi.nlm.nih.gov/30118039/)), ordinary pilot-symbol channel estimation, flat-field calibration, and Chen’s temporal-correction line.

### Kill risk

It spends measurements on calibration, interpolates badly when `t_c` is short, and is unlikely to survive novelty review.

**Ruling:** likely effective, insufficiently novel for the requested lead method.

---

## 2.5 Joint Gain-Marginal Reconstruction (`JGM-GI`)

### Mechanism

Model `log a_t` as OU/low-rank drift and jointly estimate/marginalize gain and image using Laplace–EM, variational smoothing, or alternating penalized likelihood. Use Asset-A identifiability conditions and the stationarity anchor to fix the gauge.

### Why a positive is plausible

It uses temporal correlation and all measurements, unlike framewise correction. It can be applied to Gaussian or Poisson readout.

### Closest prior art

Blind gain calibration and joint imaging/calibration are mature; see [Ling–Strohmer 2018](https://doi.org/10.1137/16M1103634), [Cambareri–Jacques 2019](https://doi.org/10.1093/imaiai/iay004), and joint Bayesian radio imaging/calibration [Roth et al. 2023](https://doi.org/10.1051/0004-6361/202346851). The Peng–Chen network methods already estimate/correct dynamic scaling factors empirically.

### Kill risk

Bilinear nonconvexity, gauge leakage, long runtime, prior sensitivity, and a high probability that the result looks like a standard latent-state model transplanted into GI.

**Ruling:** general, but not simple or aesthetically clean enough for the first campaign.

---

## 2.6 Gain-Covariance-Whitened GI (`GCW-GI`)

### Mechanism

Without estimating each gain value, use the gain covariance. Under a small-fluctuation model,

\[
\Sigma(x)
=
D_{Mx}R_aD_{Mx}+R_n.
\]

Alternate between image estimation and generalized least squares / likelihood reconstruction with `Sigma(x)^{-1/2}`-whitened patterns and buckets.

### Why a positive is plausible

The method uses the full correlation structure and may suppress the long-memory stripes that ordinary correlation or i.i.d. losses treat incorrectly.

### Closest prior art

This is generalized least squares with signal-dependent multiplicative covariance, adjacent to classical multiplicative-noise restoration and correlated-noise inverse problems. It is less crowded in GI than pilot correction, but the statistical idea is standard.

### Kill risk

The approximation is second-order; `Sigma` depends on the unknown image; strong gain fluctuations invalidate the Gaussianized model; iterative reweighting may lock onto a bad image.

**Ruling:** plausible backup, not as elegant or low-risk as scheduling.

---

## 2.7 Exact-Score One-Step RQL (`EOS-RQL`)

### Mechanism

Start from RQL and apply one exact-score / observed-information correction using the validated score machinery.

### Why a positive is plausible

If RQL is materially inefficient, a one-step estimator can close much of the gap at low cost.

### Closest prior art

Le Cam one-step efficient estimators for dead-time event detection are now explicit in [Jorgensen–Johnson 2026](https://arxiv.org/abs/2605.23210). Fisher scoring and one-step likelihood correction are classical.

### Kill risk

Probe A may show that RQL is already near the ceiling. Regularization, not likelihood efficiency, may dominate image error.

**Ruling:** merges with S2 and is unlikely to be a novel lead.

---

# 3. Ranking all independent candidates plus S1–S3

The numerical index below is a **subjective triage product**, not a calibrated probability:

\[
P(\text{positive})\times\text{novelty}\times\text{simplicity}\times\text{generality},
\]

with every factor scaled to `[0,1]`. It is used only to enforce the operator’s four-way criterion.

| rank | lane | P+ | novelty | simplicity | generality | product | ruling |
|---:|---|---:|---:|---:|---:|---:|---|
| **1** | **DOPS-GI / operator S1** | 0.82 | 0.62 | 0.95 | 0.85 | **0.410** | Only campaign-ready lane, conditional on Probe B. |
| **2** | **CPL-GI** | 0.78 | 0.55 | 0.88 | 0.68 | **0.257** | Best reserve; exact and elegant, but close to differential normalization. |
| **3** | **NOS-GI** | 0.68 | 0.78 | 0.55 | 0.85 | **0.248** | Strong theory, too complicated for first implementation. |
| **4** | **GCW-GI** | 0.70 | 0.45 | 0.68 | 0.82 | **0.176** | Reasonable backup, standard statistical core. |
| **5** | **PIS-GI** | 0.90 | 0.20 | 0.85 | 0.75 | **0.115** | Probably works; prior art makes it a poor flagship. |
| **6** | **JGM-GI** | 0.74 | 0.45 | 0.35 | 0.88 | **0.103** | General but nonconvex, slow, and crowded. |
| **7** | **operator S3: count-statistics self-calibration** | 0.65 | 0.25 | 0.82 | 0.72 | **0.096** | Useful infrastructure, not yet an image method; variance-to-mean dead-time estimation is old. |
| **8** | **EOS-RQL / operator S2** | 0.35 | 0.30 | 0.90 | 0.65 | **0.061** | Dies if RQL is already efficient; one-step estimation is established. |

## S1–S3 explicit rulings

### S1

**Promote to the top lane.** It is the only candidate that simultaneously has a strong mechanism signal, zero extra acquisition cost, no pre-scan fragility, and a defensible narrow novelty claim.

### S2

**Treat Probe A as a kill-only diagnostic.** If RQL is `>=90%` efficient under the declared image-relevant metric, close the lane. If a large gap exists, do not jump immediately to a learned model; first test the exact-score one-step update. A network can be a tool only after the analytic estimator fails.

### S3

**Do not make it the method lead.** Dead-time estimation from variance-to-mean/count statistics has clear prior art, including [Hashimoto–Ohya–Yamane 1996](https://doi.org/10.3327/jnst.33.863) and high-accuracy dead-time calibration work from [NIST](https://www.nist.gov/publications/calibrating-photon-counting-detectors-high-accuracy-background-and-deadtime-issues-0). S3 may become a useful software utility if it estimates both `tau` and hidden-hold variance robustly, but it must demonstrate an image-level benefit on deliberately miscalibrated systems before being called an imaging method.

---

# 4. Top pick: Drift-Orthogonal Pattern Scheduling

## 4.1 Frozen method statement

**Input:** a fixed accepted pattern multiset, a declared gain-drift basis or covariance, and the existing reconstruction.

**Output:** one scene-independent permutation of the exact same physical exposures.

The method computes

\[
\pi^*\in\arg\min_{\pi}\mathcal L(\pi),
\]

where `L` is the normalized object–gain cross-information loss in §2.3. For complementary illumination, the two members of each pair remain adjacent; the optimization acts on pair blocks. No scene pre-scan, learned router, or online optimization is used.

The method claim is:

> **pattern order is chosen to make the object score orthogonal to the declared temporal gain nuisance, reducing dynamic-medium information loss without changing patterns, photons, dwell, or reconstruction.**

Do not claim that randomization itself is new, that every drift is removed, or that the schedule is universal outside the declared gain class.

---

# 5. Preregistered confirmatory campaign for DOPS-GI

## 5.1 Pre-campaign DEV gate

The currently running Probe B is controlling. Launch the confirmatory campaign only if, on its frozen DEV cohort and hard drift cell, DOPS satisfies:

1. median radiometric-PSNR gain over the better of canonical and frozen-random order `>=1.0 dB`;
2. at least `5/6` DEV scenes positive;
3. no gain-free DEV scene loses more than `0.5 dB`.

If this gate fails, **do not launch a DOPS campaign and do not add an estimator to rescue it**.

## 5.2 Confirmatory cohort

Exactly **24 fresh 64×64 scenes**, six frozen strata × four instances:

1. glyph/text;
2. maze/barcode;
3. contour/thin outline;
4. low-amplitude microtexture;
5. smooth natural scenes;
6. textured natural scenes.

Five paired gain/noise seeds. No scene rejection, replacement, or post-generation filter.

## 5.3 Acquisition

- fixed complementary Hadamard pattern multiset;
- `2048` signed coefficients, realized as `4096` physical complementary exposures;
- identical dwell, photons, pattern content, and reconstruction for all schedule-only arms;
- every arm uses the same frozen regularization rule;
- the schedule is computed before any confirmatory scene is reconstructed.

## 5.4 Gain process

Primary hard condition:

\[
\log a_t=z_t-\tfrac12\operatorname{Var}(z_t),
\qquad
z_t=\varphi z_{t-1}+\epsilon_t,
\]

with stationary `CV(a)=0.20` and correlation length `t_c=64` physical frames.

No-gate secondary cells:

- `t_c in {16,256}`;
- `CV in {0.10,0.40}`;
- deterministic exponential drift matching the Peng–Chen model;
- gain-free control.

## 5.5 Arms

1. canonical Hadamard/sequency order;
2. one frozen random permutation;
3. adjacent complementary normalization baseline from the illumination-fluctuation literature;
4. published temporal-correction baseline, if the operator’s implementation is already validated;
5. **DOPS-GI**;
6. known-gain correction oracle, descriptive only.

All schedule-only arms use identical reconstruction. The learned Peng–Chen correction may be included as a no-gate context arm only if trained without confirmatory scenes.

## 5.6 Primary comparator

For image `j`, define the strong frozen baseline

\[
Q_{\mathrm{base},j}
=
\max\{Q_{\mathrm{random},j},Q_{\mathrm{pairnorm},j},Q_{\mathrm{tempcorr},j}\},
\]

using five-seed means. If temporal correction is not available before freeze, omit it everywhere rather than adding it after outcomes are seen.

Define

\[
\Delta Q_j
=Q_{\mathrm{DOPS},j}-Q_{\mathrm{base},j}.
\]

## 5.7 Primary image-level PASS

At the primary hard drift cell, all three must hold:

1. median `Delta Q >= 1.0 dB`;
2. family-stratified 10,000-replicate bootstrap 95% lower bound on the median `>0`;
3. at least `18/24` scenes have `Delta Q>0`.

## 5.8 Mandatory no-harm PASS

On the gain-free control:

1. median `Delta Q >= -0.10 dB`;
2. bootstrap lower bound `>-0.25 dB`;
3. no more than `2/24` scenes lose more than `0.50 dB`.

The method requires both the drift-rescue and no-harm conditions.

## 5.9 Secondary outputs

No gate:

- SSIM and CNR;
- the O4-A cross-information loss for every order;
- gain-model misspecification sweep;
- per-family effects;
- oracle-known-gain headroom;
- wall time to compute the offline permutation.

The mechanism metric does not rescue a failed image endpoint.

## 5.10 Decision tree

- DEV gate fails → DOPS campaign does not launch.
- Confirmatory primary or no-harm gate fails → software-method line fails; no estimator add-on or new cohort.
- Both pass → DOPS is the software method positive and may support the next manuscript.

---

# 6. Reserve lane if S1 fails

If Probe B fails, do **not** launch S2, S3, NOS-GI, or JGM-GI as multi-week campaigns.

The only justified follow-up is a one-day, six-scene **CPL-GI micro-probe**:

- Poisson photon-counting channel;
- complementary packets acquired within a declared `t_c`;
- compare exact conditional likelihood against the published adjacent-pair normalization and the operator’s temporal correction;
- kill unless median gain over the best comparator is `>=1 dB` and at least `5/6` scenes are positive.

If that micro-probe also fails, my ruling is:

> **No current software candidate clears the operator’s novelty × simplicity × generality × image-positive bar.**

At that point the correct move is not another reconstruction or router. The software layer for these frozen channels is saturated at the requested effect size, and the program should wait for a genuinely new observation model or application—not relabel an established correction.

---

# Final blunt assessment

The RLMI negative did not show that software is useless. It showed that **scene-adaptive surrogate optimization is the wrong kind of software here**.

The best remaining software method is almost embarrassingly simple:

> **use exactly the same patterns and photons, but order them so the object cannot masquerade as temporal gain drift.**

That simplicity is a strength. If it produces a preregistered image gain, it is a clean method paper. If it does not, the honest conclusion is that the remaining software ideas are either known corrections, indirect calibration utilities, or mathematically attractive estimators without enough image-level leverage.