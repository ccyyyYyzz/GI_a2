# Adversarial anomaly-fidelity probe — decisive feasibility test

**Date:** 2026-07-22 · **Mode:** strictly read-only on `E:\`; all writing under
`D:/GI_another/results/round63_next/ANOMALY_PROBE/`. No commits.
**Runtime:** ~66 s (three regimes, CPU). **Code:** `anomaly_probe.py` (self-contained),
`make_figures.py`. **Data:** `anomaly_probe_results.json`. **Figure:**
`anomaly_fidelity_3panel.png`.

---

## VERDICT (conclusion first)

**The improvement claim HAS LEGS for resolvable anomalies, through the certified
mechanism — but the effect is a preservation guarantee, not a background-quality
win, and the "prior invents plausible detail" concern is not the failure mode at
5% sampling; the failure mode is erasure.**

Per anomaly class, with numbers (clean regime = the operator's faithful
noiseless-design point, peak ~1e6 counts, measurement SNR ~87; real lane0 rate05
operator + real VQAE prior; 12 images × 2 anomalies × 3 seeds = 72 instances):

- **Resolvable anomalies (4–6 px blobs / 3×8 bars / 5 px rings) — FIDELITY ROUTE WORKS.**
  - The unconstrained DL prior (**A2**) **erases ~30% of the anomaly content the
    data actually measured**: recovery efficiency ρ/witnessed-fraction **0.72**
    (VQAE) / **0.67** (VQGAN) vs **1.00–1.05** for the information-guided arm (**A3**).
  - Anomaly **survival** (matched-filter detection > 5σ): **A3 83%** vs **A2 53%**
    (VQAE) / **28%** (VQGAN). The unconstrained prior removes resolvable anomalies
    in roughly **half to three-quarters** of cases; A3 in 17%.
  - Anomaly-region **bbox PSNR: A3 − A2 = +0.65 dB** (VQAE) / **+1.00 dB** (VQGAN),
    **A3 wins 83–89%** of paired cases. Modest because the *specified sharp 4–6 px
    shapes* are only ~46% witnessed by this operator; the gap **scales with the
    witnessed fraction** (corr 0.54–0.63; high-witnessed tertile **+0.80 / +1.45 dB**),
    so broader anomalies (the audit's σ3 bump reaches 96% witnessed) would show a
    larger fidelity advantage.
- **Sub-resolution anomalies (1–2 px points) — FIDELITY ROUTE DOES NOT APPLY;
  FLAGGING ROUTE WORKS PERFECTLY.** Witnessed fraction ~0.085 (null-dominated);
  **survival ≈ 0% for every arm** — no reconstruction can preserve them. But the
  range/null certificate flags them: the per-anomaly **witnessed fraction
  separates the two classes with AUC = 1.000** (min resolvable 0.286 > max
  sub-resolution 0.135, margin +0.15).
- **The certified mechanism is exact and legible.** Across all resolvable
  anomalies, **A3's recovered fraction ≈ its witnessed fraction** (linear fit
  slope **1.09 / 0.97**, intercept ≈ 0): the data-subspace lock recovers precisely
  the range component of the anomaly, no more and no less. A2 recovers *less than
  the witnessed fraction* — it overwrites even measured directions.
- **Background concession: none at this operating point.** A3 background PSNR ≥ A2
  (A3 − A2 = **+0.48 dB** VQAE / **+1.28 dB** VQGAN outside the anomaly boxes);
  A3 − A1 ≈ 0. So A3 **dominates** A2 here — better anomaly fidelity *and* equal/
  better background — because the operator's VQAE/VQGAN prior is itself lossy
  (~21 dB) and buys no background advantage from overwriting the data. This is
  *stronger* than the proposal's "preserve anomalies at a modest background cost"
  framing (caveat below).
- **Hallucination check on A2 — erasure YES, invention NO (at 5%).** A2 erases
  resolvable anomalies (above). It does **not** inject false anomaly-strength
  structure in anomaly-free regions: on 6 clean control images the blob detector
  fires **0 times for A2** vs **33 natural blob structures in the truth images**
  (the 5%-sampled prior is bandwidth-limited). The going-in "invents plausible
  detail" story is therefore only half right: at 5% sampling the DL prior's
  measured failure mode is **removing** atypical content, not fabricating it.

**Bottom line for the direction:** the improvement leg is real and certifiable for
resolvable anomalies via `ρ_A3 = witnessed-fraction`; the certificate leg
(witnessed fraction as an anomaly-flag) is perfect for the sub-resolution regime.
Both legs survive. The headline that needs rewording is "unconstrained DL *erases
or invents*" → at this sampling it is specifically **erasure**, and A3's advantage
is **anomaly preservation at zero background cost**, with the fidelity gap tunable
by the operator's spread-spectrum design (which sets the witnessed fraction).

---

## Design (frozen before running)

- **Operator.** Real lane0 **rate05** structured GI operator `A` (m=205, n=4096,
  sampling ratio **5.0%**), rebuilt verbatim in numpy (DC + low-freq DCT + low-
  sequency Hadamard + random ±1 rows) and **SHA-256 verified** against the frozen
  receipt: `62305c8e…e79d6` = match. Orthonormal row basis `Q` (numerical rank
  **200**) by SVD → range projector `P_R=QᵀQ`, null projector `P_N=I−P_R`.
  Verified `A·truth = y` to 1.7e-7 (the cache is **noiseless**, `noise_std=0`) and
  `A·(A⁺y)=y` to 3.9e-8 (min-norm range estimate is exactly data-consistent).
- **Prior.** Real lane0 **VQAE** (`priors/vqae.pt`, codebook 256 / z 64 / base 48),
  loaded with 0 missing / 0 unexpected keys; a vector-quantized autoencoder that
  projects onto the natural-image manifold (clean-image recon 20–29 dB — a genuine
  lossy manifold map). Robustness rerun with the more aggressive **VQGAN** prior
  (`priors/vqgan.pt`, same architecture, adversarially trained).
- **Images.** `truth` from the frozen `rate05_test_cache.pt` — STL-10 **test
  split** (`source_index ∈ [0,7999]`), disjoint from the priors' `train+unlabeled`
  training pool → never used for prior training. 20 deterministic indices
  (seed 20260722): **2 freeze** (`[915,1197]`, hyperparameters frozen here), **6
  anomaly-free controls** (`[1231,1649,2440,2759,3969,4200]`), **12 test**
  (`[4209,4225,4924,5030,5203,5243,5275,5764,6309,6391,6677,6700]`).
- **Anomaly planting (frozen, seed 90210).** Per image: **1 resolvable** +
  **1 sub-resolution**, non-overlapping.
  - Resolvable set (cycled): `blob5_bright` (5 px disk), `blob6_dark` (6 px dark
    disk), `bar_3x8_bright`, `ring_r5_bright`. Contrast rule: additive Δ=**0.5**
    (bright +, dark −), then clip to [0,1]; the realized post-clip change is the
    template `a`.
  - Sub-resolution set: `point_1px`, `point_2px`, additive Δ=**0.7**.
  - Planted in the **true image before measurement**, so the data genuinely
    contains the anomaly.
- **Measurement.** Honest **differential-bucket Poisson** counting channel: each
  zero-mean row is split `a=a⁺−a⁻`, buckets `B±~Poisson(N0·⟨a±,x⟩)`,
  `y=(B⁺−B⁻)/N0` (unbiased, `E[y]=Ax`). **Primary "clean" regime** N0=34757
  (peak ~1e6 counts, SNR ~87) — matches the operator's noiseless design and
  isolates the prior-vs-data split. **Secondary "noisy" stress** N0=348 (peak ~1e4,
  SNR ~8). 3 noise seeds {101,202,303}.
- **Arms (same `y`, same prior, same data estimate `x_range=A⁺y`; differ only in
  how range/null are filled).**
  - **A1 data-consistent baseline:** range = data, null = 0 (min-norm), box-fiber
    projected onto `{Ax=y}∩[0,1]ⁿ`. **No learned prior.**
  - **A2 unconstrained DL:** `x = VQAE(clip(x_range))` — the prior rewrites the
    **full image** including measured directions; **not** re-projected onto the data.
  - **A3 information-guided:** `x = box_fiber_project( x_range + P_N(prior−base) )`
    — **range locked to the data**, the same prior writes **only the null-space**
    component; exact `Ax=y` and box feasibility enforced. Reimplemented in-probe
    (their `apply_generator` geometry); their code untouched.
  A2 and A3 thus share the null component and **differ only in whether the measured
  (range) directions come from the prior (A2) or the data (A3)** — an exact
  isolation of the mechanism.
- **Frozen readout thresholds.** Survival = matched-filter response > **5σ**
  (σ = per-arm noise floor from control-image independent-seed reconstruction
  differences). Witnessed-fraction reported as a continuous certificate (class
  separation is threshold-free). Projection iters (50) and N0 frozen on the 2
  freeze images by background-quality/record-error only — **no tuning on anomaly
  outcomes.**

---

## Readouts (clean regime, VQAE prior unless noted; per anomaly class)

### 1 + 2. Anomaly-region fidelity & survival

| class | arm | bbox PSNR (dB) | recovered ρ | efficiency ρ/wf | survival |
|---|---|---|---|---|---|
| resolvable (wf 0.46) | A1 (data-only) | 16.26 | 0.465 | 1.01 | 83% |
| resolvable | **A2 (uncon. DL)** | 15.54 | **0.339** | **0.72** | **53%** |
| resolvable | **A3 (info-guided)** | 16.18 | **0.465** | **1.00** | **83%** |
| sub-res (wf 0.085) | A1 | 16.52 | 0.084 | 0.93 | 8% |
| sub-res | A2 | 16.05 | 0.027 | 0.19 | 0% |
| sub-res | A3 | 16.60 | 0.073 | 0.71 | 0% |

**Headline (paired, resolvable):** A3 − A2 bbox PSNR **+0.65 dB** (median +0.68),
**win 89%**, Δρ **+0.126**. With the **VQGAN** prior: **+1.00 dB** (median +1.20),
win 83%, Δρ +0.171, and A2 survival drops to **28%** (efficiency 0.67).

**Witnessed-fraction dependence (validates extrapolation):**

| wf tertile | mean wf | A3−A2 bbox (VQAE) | A3−A2 bbox (VQGAN) |
|---|---|---|---|
| low  | 0.36 | +0.32 dB | +0.18 dB |
| mid  | 0.48 | +0.82 dB | +1.38 dB |
| high | 0.55 | +0.80 dB | +1.45 dB |

Per shape (VQGAN), the erasure is sharpest for the least-natural anomalies:
`blob5_bright` A3 67% vs A2 0% survival; `ring_r5_bright` 67% vs 0%;
`bar_3x8_bright` 100% vs 33%; `blob6_dark` (highest wf) 100% vs 78%.

### 3. Background quality (PSNR outside anomaly boxes)

| arm | VQAE | VQGAN |
|---|---|---|
| A1 | 21.84 | 21.84 |
| A2 | 21.46 | 20.08 |
| A3 | 21.94 | 21.36 |

**A3 pays no background concession** (A3 − A2 = **+0.48 / +1.28 dB**; A3 − A1
≈ 0 / −0.48 dB). At this operating point A3 dominates A2.

### 4. Flagging correctness (certificate decision quality)

- `ρ_A3 = witnessed-fraction` linear fit: slope **1.09** (VQAE) / **0.97** (VQGAN),
  intercept ≈ 0 → A3 recovers exactly the witnessed (range) fraction.
- Witnessed fraction: resolvable mean 0.461 (0.286–0.578), sub-resolution mean
  0.085 (max 0.135). **Class AUC = 1.000**, margin **+0.150** — the certificate
  perfectly distinguishes "data can vouch" (resolvable) from "data cannot vouch"
  (sub-resolution). (A naive flag threshold of 0.5 over-flags 27/36 resolvable
  because these sharp shapes are only ~46% witnessed; any threshold in
  [0.135, 0.286] separates the classes exactly — the certificate is better read as
  the continuous witnessed fraction.)

### 5. Hallucination check on A2

- **Removal (erasure):** YES — resolvable survival A2 53% (VQAE) / 28% (VQGAN) vs
  A3 83%; A2 discards ~30% of witnessed anomaly content (efficiency 0.72 / 0.67).
- **Invention:** NO at 5% sampling — on anomaly-free controls the blob detector
  (threshold = 99.5th pct of natural truth-image blob responses) fires **0** for
  A2 vs **33** natural structures in the truth images and 1 for A1/A3. The prior is
  bandwidth-limited at this sampling; its measured failure mode is erasure, not
  fabrication. (In the noisy stress regime the *data-consistent* arms A1/A3 show
  ~230 noise-driven false blobs while A2 denoises them to 0 — the honest flip side:
  the prior's denoising is valuable under heavy noise, at the cost of anomaly
  erasure.)

---

## Honesty / limitations

- **Effect size is bounded by the specified anomaly set.** The prescribed sharp
  4–6 px shapes are only ~46% witnessed by this operator, so the *absolute* bbox
  PSNR gap is modest (+0.65–1.00 dB). The scale-free efficiency (0.72 vs 1.00) and
  survival (53% vs 83%) show the mechanism cleanly, and the gap grows with
  witnessed fraction; broader/smoother anomalies would show larger PSNR gaps. This
  is a genuine finding, not a strawman — the fidelity route's payoff is an
  operator-design (spread-spectrum) knob.
- **No background concession is partly an artifact of a lossy prior.** The
  operator's VQAE/VQGAN reconstruct clean images at only ~20–29 dB, so A2 (pure
  prior) has no background advantage to trade. A stronger generative prior
  (e.g. diffusion) could reintroduce a fidelity↔background tradeoff; **untested
  here** — the "modest background cost" framing may still hold against a
  higher-capacity prior.
- **Noisy regime confounds bbox PSNR.** At SNR ~8 the AC-coupled measurement is so
  noisy that A2's smoothing wins bbox PSNR by denoising (A3 − A2 = −5.3 dB) even
  though A2 recovers *less* anomaly (ρ 0.24 vs 0.43); detection survival collapses
  to 0% for all arms (below 5σ). The noise-robust recovered-fraction still shows
  A3 > A2. The operator's real regime is noiseless, so the clean regime is the
  faithful primary; the noisy case is a stress test, not the operating point.
- **Poisson is a differential-bucket surrogate**, chosen to keep the channel honest
  (unbiased, shot-noise-scaled); the probe's point is the prior-vs-data split, not
  deep noise physics.
- **Scope:** all-simulation; single operator/lane (lane0 rate05); STL-10 grayscale
  64×64; the prior is the operator's own VQAE/VQGAN. A2/A3 use a single decode-
  encode of the prior (not their full measurement-conditioned refiner pipeline);
  the isolation is deliberate (range-source is the only difference).

---

## Provenance / reproduction

- Run: `conda run -n pytorch python anomaly_probe.py` (writes
  `anomaly_probe_results.json`, `_recon_cache.npz`, `_plans.json`), then
  `python make_figures.py` (writes `anomaly_fidelity_3panel.png`).
- Read-only E: assets: operator config
  `E:/GAN_FCC_WORK/artifacts/gan_gi_journal_round59_recovery/lane0/content/gan_operator_assets/config_rate05.yaml`;
  priors `…/priors/vqae.pt`, `…/priors/vqgan.pt`; test cache
  `…/lane0/caches/rate05_test_cache.pt`. Operator SHA `62305c8e…e79d6` (match).
- Template for the operator + range/null split: the audit's `anomaly_split_test.py`
  (`…/GAN_LINE_AUDIT/`).
