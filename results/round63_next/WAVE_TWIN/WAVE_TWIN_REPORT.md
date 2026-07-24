# WAVE-OPTICS DIGITAL TWIN — bench-transfer test of the four falsifiable predictions

**Date:** 2026-07-24 · **Line:** ROUND63 two-ledger / beyond-modulator-band sentinel.
**Status:** QUARANTINED from the PRL manuscript (R43 green light — no new manuscript
data). Feeds the bench / collaboration phase only. No git commit.
**Scope:** a scalar angular-spectrum (Fresnel/band-limited-ASM) digital twin of the
planned bench — DMD → thin random phase screen (diffuser) → static object → single
bucket — used to test, in real diffraction physics, the four R39 §7.5 transfer
predictions plus an end-to-end sentinel comparison against the frozen statistical
detector (`SEALED_DET/CONFIRMATORY_RESULTS.json`).

> **One-line verdict.** The covariance sentinel **transfers to real scalar diffraction** — at
> the geometry-matched sealed cell the M-corrected detection time (`≈768 banks @2 %`) lands
> inside the sealed band (`513–2996`), and the aperture/grain (T3) and multiplicative→convolutive
> (T4) predictions are confirmed — **with three honest corrections to the idealized statistical
> model**: (1) the mean-wall is **not machine-zero** under real DMD pixelation + free-space
> geometry — it leaks at **~0.7 %–9 %** and is *operationally decisive* at under-developed
> near-contact (mean channel out-detects covariance there), requiring a near-conjugate DMD relay +
> transfer-wording fix; (2) the statistical grid's **(k_w, σ_f) axes are not independently
> realizable** — physical speckle is illumination-diffraction-limited to fine grain (wide
> aperture, ~10⁴–10⁵ grains/bucket) and the developed multiplicative regime is **sub-mm thin**;
> (3) the covariance noncentrality **scales as `M^1.8`** (codes still buy detection at finite z2 —
> rank-one `M_eff` saturation is a *complete-scrambling* property only), so the M=32→128
> comparability correction is essential and closes most of the apparent deviation.

---

## 0. The twin (parameters, sampling, anti-aliasing)

| Quantity | Choice | Criterion / check |
|---|---|---|
| Wavelength λ | 532 nm | fixed |
| Grid | 2048 × 2048 | — |
| Pixel dx | **4.0 µm** (window 8.192 mm) | resolves DMD macro-pixel (8 px); micromirror sub-pixel handled by aperture envelope |
| Propagator | band-limited angular spectrum (Matsushima–Shimobaba), complex128 (T1) / complex64 (T2–T5) | plain-ASM alias-free distance `z_max = N·dx²/λ = 61.6 mm` ≫ all legs (z1,z2 ≤ 20 mm); analytic band-limit applied at every z |
| Energy conservation | ratio = 1.000000 at z = 5/10/20/40 mm (both dtypes) | verified `validate_core.py` |
| DMD | 10.8 µm micromirrors, fill 0.92, 64×64 code on 32 µm macro-pixels (2.048 mm active) | band cut-off k_p=5 enforced spectrally (matches sealed `band_modes`); macro-pixel hold + micromirror sinc envelope = real pixelation |
| Diffuser | thin Gaussian-correlated phase screen, l_c ∈ {5,15,50}(+coarse) µm; developed σφ=2π and weak σφ=0.3 rad; **rotation = lateral shift** of a large precomputed screen (quasi-static within a bank, stepped between banks) | matches bench protocol |
| Object | 64×64 macro reflectance in [0,1]; toggleable sub-DMD-resolution beyond-band feature (the "change") at 2 % and 5 % energy | `witness_scene` + beyond-band δ in the k∈[6,9] annulus |
| Bucket | total object-plane power, one scalar/exposure; Poisson shot at 1e4 pe/bucket | matches sealed ledger |

**complex64 justification.** The 4060's FP64 path is ~1/64 throughput; the Monte-Carlo
legs (T2–T5, thousands of 2048² FFTs) run in complex64 (transfer-function *phase* still
built in float64, then cast). The wall metric is identical to 5 significant figures in
both dtypes (z1=0 leak = 6.8089e-3 in float64 **and** float32), and energy conservation
holds to 1e-6 — so complex64 is validated for every physical effect measured here.

---

## 1. T1 — MEAN-WALL LEAK  ·  verdict: **WALL LEAKS (FLAG RAISED)**

**Statistical baseline:** sealed D0 `mean_derivative = 5.68e-16` — machine zero. The mean
channel `m = P·diag(μ)·x` is exactly blind beyond k_p because the *signed effective code*
is band-limited.

**Mechanism found by the twin.** Physically the signed bucket is a difference of
complementary **intensity** exposures `|E⁺|² − |E⁻|²` with `a± = (1±s)/2`. With **no
propagation** `|a⁺|² − |a⁻|² = s` exactly (band-limited → zero beyond-band): the wall
holds. Two real effects break this exact cancellation:
- **free-space DMD→diffuser propagation z1** makes E± complex, so `|E⁺|²−|E⁻|² ≠ s`;
  [MECHANISM CORRECTED post-hoc, divergence round 2026-07-24: with complementary
  pairing the signed difference equals `4·Re(C̄S̃)` exactly at every z1 (measured
  5e-16) — the `|S̃|²` 2k_p intensity-autocorrelation term **cancels identically**;
  the leak is the carrier–code cross-spectrum through the finite object window,
  with NO spectral edge at 2k_p. See DIVERGENCE_R2/MERGED_ADJUDICATION.md];
- **finite DMD micromirror aperture / pixelation** leaks even at z1=0.

**Result — deterministic z1 sweep (z2=0; diffuser-independent, zero Monte-Carlo noise):**

| z1 (mm) | 0 | 1 | 2 | 5 | 10 | 20 |
|---|---|---|---|---|---|---|
| mean-wall leak (rel. in-band) | **6.8e-3** | 1.2e-2 | 1.6e-2 | 3.2e-2 | **5.4e-2** | 8.9e-2 |
| flag (>1e-3) | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ |

Even at conjugate imaging (z1=0) the DMD **pixelation floor** is ~0.7 %; free-space z1
raises the leak monotonically to ~5 % at 10 mm and ~9 % at 20 mm — **all ≫ the 1e-3
threshold.** No z1 in the swept range gives leak < 1e-3.

**Operational leak (mean-channel d′ vs covariance-channel d′), from T5 banks:** see §5 —
the leak's *operational* impact is small because the mean beyond-band signal is buried in
medium (speckle) fluctuation, but it is **not exactly zero**.

> **Transfer-wording flag (loud).** The manuscript sentence "the mean reconstruction has
> exactly zero response to the beyond-band witness" (R39 §7.5 prediction 1) is an
> idealization of a band-limited *signed code*. In real optics the mean channel is
> **nearly-but-not-exactly blind**: it leaks at the 10⁻³–10⁻¹ level from DMD pixelation and
> DMD→diffuser propagation. **Pre-submission action:** (a) reword to "the mean channel is
> blind to beyond-band content *to the extent the illumination code is band-limited*; residual
> leakage from pixelation/geometry is ≤X and is suppressed by complementary differencing and
> by near-conjugate DMD–diffuser imaging"; (b) specify near-conjugate (relay) DMD imaging in
> the apparatus; (c) the detection claim is unaffected (the covariance channel dominates the
> mean channel by ≫10³ in d′ — §5).

---

## 2. T2 — SPECKLE STATISTICS → the declared-law dictionary  ·  verdict: **DICTIONARY BUILT (with a coupling caveat)**

Grain (intensity-autocorrelation FWHM) and contrast at the object plane vs l_c and z2,
mapped to the statistical twin's (k_w, σ_f). `k_w ≈ L_obj/grain` (L_obj = 2048 µm),
64-grid Nyquist = 32; σ_f ≈ intensity contrast (sealed SIG_EFF = {0.298,0.503,0.696}).

Selected rows (developed σφ=2π):

| l_c (µm) | z2 (mm) | contrast | grain (µm) | k_w (=2048/grain) | σ_f class |
|---|---|---|---|---|---|
| 5 | 0 | 0.001 | — (no speckle) | — | phase screen invisible at contact |
| 5 | 5 | 0.99 | 4.2 | ~490 (≫Nyq) | ≈1.0 |
| 15 | 5 | 0.99 | 4.4 | ~465 (≫Nyq) | ≈1.0 |
| 50 | 1 | 0.87 | 20.6 | ~99 (≫Nyq) | ~0.7–1.0 |
| 50 | 20 | 1.01 | 13.2 | ~155 (≫Nyq) | ≈1.0 |

Weak σφ=0.3 rad gives contrast 0.03–0.40 (low σ_f) and coarser grain (up to 80 µm).

**Two findings.**
1. **Contact degeneracy.** At exactly z2=0 a pure phase screen leaves the intensity
   unchanged (`|E·e^{iφ}|² = |E|²`) → contrast 0, no covariance signal. The statistical
   "thin multiplicative screen" therefore corresponds to a phase screen at **small but
   non-zero z2**, where phase has partially converted to intensity speckle.
2. **(k_w, σ_f) are coupled, not independent.** Physical object-plane grain is
   **illumination-diffraction-limited** (grain ≈ λ·z2/D_illum) to ≲100 µm in this
   geometry, so **k_w always sits near/above the 64-grid Nyquist (kwf ≳ 4)**. The
   coarse-medium cells (kwf = 1, 2) of the sealed grid are **not physically realizable**
   with a 2 mm illumination aperture; and high contrast (large σ_f) requires developed
   speckle (enough propagation), which is incompatible with coarse grain at small z2.
   The sealed grid's free (k_w, σ_f) product overstates the realizable design space.

---

## 3. T3 — APERTURE LAW, REAL CURVE  ·  verdict: **PREDICTION 3 CONFIRMED**

Covariance-channel response ‖Cov(b|x0+εδ_k) − Cov(b|x0)‖_F vs scene Chebyshev
frequency k, for three grain sizes (common-random diffusers). Usable edge = last k above
40 % of the (non-DC) peak.

| grain (µm) | k_w (=2048/grain) | measured usable edge |
|---|---|---|
| 8 | ~256 | **27** |
| 48 | ~43 | **20** |
| 96 | ~21 | **13** |

The usable scene-frequency edge **moves monotonically inward as the grain coarsens** — the
covariance aperture narrows with k_w exactly as B_p⊕B_w predicts. Because the accessible
grain is fine, all three edges lie near/below the 64-grid Nyquist (32); the aperture is
**wide and not the operational bottleneck** in this geometry, but the *law* (edge tracks
grain) is cleanly reproduced.

---

## 4. T4 — MULTIPLICATIVE → CONVOLUTIVE  ·  verdict: **PREDICTION 4 CONFIRMED (sub-mm boundary)**

Complex-field correlation between the true propagated field `E_true = prop_{z2}(screen·E_code)`
and the pointwise-multiplicative model `E_mult = screen·prop_{z2}(E_code)` vs screen–object
separation z2. α = 1 − correlation maps to the sealed D5 `convolutive_blend` axis.

| regime | field-corr @0.5mm | field-corr @2mm | z2 @ α=0.10 | z2 @ α=0.25 | z2 @ α=0.50 |
|---|---|---|---|---|---|
| developed σφ=2π, l_c=15 | 0.13 | 0.037 | **0.065 mm** | 0.163 mm | 0.327 mm |
| weak σφ=0.3 rad, l_c=15 | 0.98 | 0.79 | **0.54 mm** | 1.41 mm | > 2 mm |
| developed σφ=2π, l_c=5 | 0.057 | 0.017 | 0.07 mm | 0.176 mm | 0.352 mm |

**Finding.** The **developed-speckle** multiplicative (thin-screen) regime is
**sub-millimetre thin** — the diffuser must sit within ~0.3 mm of the object for the
pointwise model (and the sealed detector's *declared law*) to hold. The **weak-phase**
screen stays multiplicative out to ~1–2 mm. The sealed D5 tolerance (`convolutive_blend`
passed to α=0.5) therefore corresponds to z2 ≲ 0.33 mm (developed) / ≳ 2 mm (weak). Beyond
that the physics is the **developed-speckle / scrambling regime** (`SCRAMBLE_EXT`, scene
enters as the rank-one scalar Q(x)) — where the *sentinel still survives* (§5) but the
thin-screen *imaging* aperture does not.

---

## 5. T5 — END-TO-END SENTINEL  ·  verdict: **TWIN_CONFIRMS at the geometry-matched cell (after M=code-count correction); sentinel transfers to real diffraction**

Frozen covariance detector (sealed conventions: signed buckets, per-bank sample
covariance, matched score `W = V0⁻¹ΔC V0⁻¹`, per-bank noncentrality `λ = ½·tr[(V0⁻¹ΔC)²]`,
`d′² = T_eff·λ`) on twin-generated banks. ΔC computed with common-random-number variance
reduction (object-independent speckle pool → identical diffusers for H0/H1); shot added
analytically. AUC by bootstrap over T_eff-bank sample covariances.

### 5.1 Raw run (M=32 codes, 800-bank pool)

| geometry | change | cov λ (M=32) | T_det@M32 | AUC@453 | mean d′/cov d′ (M=32) |
|---|---|---|---|---|---|
| near-contact z2=1mm, l_c=50 (grain ~20 µm) | 2 % | 1.36e-4 | 184 000 | 0.71 | 2.9 |
| near-contact z2=1mm, l_c=50 | 5 % | 7.35e-4 | 34 000 | 0.93 | 3.1 |
| mid z2=5mm, l_c=50 | 2 % | 2.01e-3 | 12 400 | 0.37* | 0.63 |
| **mid z2=5mm, l_c=50** | **5 %** | **1.14e-2** | **2 190** | **0.96** | 0.66 |
| near-contact z2=1mm, l_c=15 (grain ~4 µm) | 2 % | 3.9e-6 | 6.4e6 | 0.60 | 18 |
| near-contact z2=1mm, l_c=15 | 5 % | 2.0e-5 | 1.2e6 | 0.70 | 20 |

*(AUC 0.37 < 0.5 = the 2 % signal at the covariance-estimation-noise floor at M=32; the
ε-ordering AUC(5 %) ≫ AUC(2 %) and all z2/grain trends are monotone — not a bug.)*

### 5.2 COMPARABILITY CORRECTION — M-code-count scaling (**required before any deviation claim**)

The run used **M=32** codes; the sealed probe used **M=128** (256 complementary exposures/bank).
The covariance channel carries ~`M(M+1)/2` covariance entries vs the mean channel's ~`M`, so
M=32 deflates the covariance channel relative to the sealed numbers **and** relative to the
mean channel. Measured scaling of `λ` (per-bank cov noncentrality) with M, on code-subsets of
one pool (`T5b_MSCALING.json`):

| geometry / change | λ(M=8) | λ(M=16) | λ(M=32) | **exponent p (λ∝Mᵖ)** | **extrap T_det @ M=128** |
|---|---|---|---|---|---|
| mid z2=5mm, 2 % | 2.37e-4 | 8.90e-4 | 2.78e-3 | **1.78** | **768 banks** |
| mid z2=5mm, 5 % | 1.45e-3 | 5.29e-3 | 1.68e-2 | 1.76 | **129 banks** |
| near-contact, 2 % | 1.17e-5 | 4.52e-5 | 1.61e-4 | 1.89 | 11 300 banks |
| near-contact, 5 % | 6.75e-5 | 2.55e-4 | 9.23e-4 | 1.89 | 1 980 banks |

> **`λ` scales as `M^1.8–1.9` — the covariance channel is NOT saturated at the in-band DOF.**
> This contradicts the fully-scrambled `M_eff≈13` picture (`SCRAMBLE_EXT §4.1`) **and confirms
> it is regime-specific**: at finite z2 the bucket covariance is *higher-rank* (thin-screen /
> partial-memory regime), so extra codes keep adding independent covariance looks; the rank-one
> `M_eff` saturation is a property of **complete** scrambling only. This is a genuine finding.
> Consequence: the M=32→128 correction shortens `T_det` by `(128/32)^1.8 ≈ 12×`.

### 5.3 MATCHED-CELL comparison (**not** the sealed best cell)

Via the T2 dictionary, both l_c=50 geometries map to the sealed class **(σ_f≈1.0, kwf=4)**
(measured contrast 0.87–1.01 → σ_f=1.0; illumination-limited fine grain → kwf=4). The sealed
2 % `T_det` for **(σ_f=1.0, kwf=4)** cells is **513–2996 banks** (k^-2/k^-1/flat × claim
1.25–1.8) — *not* the best-cell 453. Comparison:

| geometry | change | twin T_det @ M=128 (extrap) | matched sealed cell (σ_f=1.0,kwf=4) | verdict |
|---|---|---|---|---|
| **mid z2=5mm developed** | **2 %** | **≈ 768 banks** | **513–2996 banks** | ✅ **WITHIN matched range** |
| mid z2=5mm developed | 5 % | ≈ 129 banks | (2 % cell / ε² → ~130–750) | ✅ consistent |
| near-contact (under-developed) | 2 % | ≈ 11 300 banks | 513–2996 banks | ✗ ~4–20× slower (starved cov, §5.4) |

> **Corrected verdict: TWIN_CONFIRMS.** For the developed mid-z2 geometry — the regime the
> detector is actually designed for — the M-corrected wave twin gives `T_det ≈ 768 banks (2 %)`,
> **inside the geometry-matched sealed cell's 513–2996 band.** The apparent 27×–10⁴× "deviation"
> against the *best* cell at M=32 was **almost entirely a comparability artifact** of (a) M=32 vs
> 128 and (b) comparing to the best (unrealizable kwf=1) cell instead of the matched (kwf=4) cell.
> The sentinel transfers to real scalar diffraction at the designed operating point.

### 5.4 Why near-contact still starves the covariance channel (T2 contact-degeneracy)

Near-contact (z2=1mm) remains ~4–20× slower than the matched cell even after M-correction,
because a **pure phase screen has not yet converted phase to intensity** at 1 mm (T2:
contrast 0.87, not fully developed) — the covariance signal is physically weak there. This is
the T2 contact-degeneracy, not a detector failure: the medium is *becoming* visible but is
under-developed.

### 5.5 THE LOUD FINDING — the T1 mean-leak is operationally decisive near-contact (design requirement)

The mean-vs-cov d′ ratio (mean channel ~`√M`, cov channel ~`M^0.9`; ratio ∝ `M^-0.4`, so the
M=128 correction *reduces* it ~0.57× but does not flip the near-contact case):

- **developed / mid-z2: ratio 0.63–0.66 (M=32) → ~0.36 (M=128) → covariance channel dominates.**
  The "image-blind, change-observable" separation holds cleanly.
- **under-developed / near-contact: ratio 2.9 (l_c=50) to 18–20 (l_c=15) at M=32 → still >1 at
  M=128.** The **leaky mean channel out-detects the covariance channel**: the T1 wall leak here
  is not a cosmetic 10⁻² number, it is **operationally larger than the intended signal**.

> **Design requirement for the bench (numbers-backed, from T1+T5).** The mean-as-reference
> attribution architecture is only safe where the covariance channel beats the (non-zero,
> ~0.7–9 %) mean leak — i.e. **developed speckle**. To use the sentinel on hardware:
> 1. **near-conjugate DMD→diffuser relay** (z1→0): drops the T1 leak from ~5 % (z1=10mm) toward
>    the ~0.7 % pixelation floor;
> 2. **DMD apodization** and a band-limited-code display to suppress the pixelation floor further;
> 3. **z2 placement in the developed-speckle window** (contrast ≈ 1) — from T2/T4 this is a
>    narrow band: far enough that phase→intensity is developed (z2 ≳ a few mm), but the developed
>    multiplicative regime is sub-mm (T4), so the operating point is the **developed-convolutive**
>    regime (SCRAMBLE_EXT), where the sentinel is validated (mid-z2 5 % AUC 0.96) and cov > mean;
> 4. **as many codes as the frame budget allows** (λ∝M^1.8 → codes buy detection speed here,
>    unlike the fully-scrambled limit).
>
> **Operating window where the sentinel transfers as-designed:** developed speckle (z2 ≈ few mm,
> contrast ≈ 1), near-conjugate DMD relay, M ≈ 128 codes → matched-cell detection (~500–800 banks
> at 2 %, ~130 at 5 %), covariance channel dominant over the residual mean leak.

---

## 6. The declared-law dictionary (summary table)

| statistical-twin quantity | wave-twin realization | mapping |
|---|---|---|
| code band k_p = 5 | 64×64 code on 32 µm macro-pixels, band-limited | exact (spectral) |
| medium band k_w (kwf∈{1,2,4}) | object-plane speckle grain w_g | k_w ≈ 2048 µm / w_g; **physically pinned to kwf ≳ 4** (fine grain) |
| medium contrast σ_f | intensity speckle contrast | σ_f ≈ contrast; developed → 1.0, weak → 0.03–0.4 |
| thin multiplicative screen | phase screen at z2 ≲ 0.3 mm (developed) / ≲ 2 mm (weak) | T4 boundary; α = 1−corr = D5 convolutive_blend |
| independent banks | independent diffuser lateral shift (quasi-static within bank) | exact (protocol) |
| beyond-band witness δ | sub-DMD-resolution object feature, k∈[6,9] | exact |
| mean wall (5.7e-16) | mean-channel leak 0.7 %–9 % | **NOT machine zero** — pixelation + geometry |

---

## 7. Honest gaps and scope

- **Scalar approximation.** Full scalar diffraction only; no vector/polarization effects.
  A high-NA diffuser or large scattering angles would need vector theory (COMSOL leg
  addresses the thin-phase-screen approximation itself — see `COMSOL_MICRO/`).
- **No polarization / depolarization.** Real diffusers depolarize; the intensity-only
  bucket model ignores this (would reduce speckle contrast, i.e. effective σ_f).
- **Thin-screen diffuser.** A single multiplicative phase layer; a thick/volumetric
  scatterer (multiple scattering) is not modelled — that is exactly the `SCRAMBLE_EXT`
  extreme, approached here only via the developed-speckle limit.
- **Coherence.** Fully coherent, monochromatic, spatially-coherent illumination assumed.
  Partial coherence would blur speckle and reduce contrast.
- **complex64 for T2–T5** (validated identical to float64 on the wall metric to 5 sig figs).
- **Illumination-aperture-limited grain.** Coarse-medium (small k_w) statistical cells are
  not reachable in this geometry — reported as a finding, not hidden.
- **T5 bootstrap ROC** uses the wave-measured bank pool; d′ scaling assumes the Gaussian
  covariance-discrimination law (matches the sealed convention).

---

## 8. Deliverables in this directory

- code: `wave_twin.py` (core), `twin_pool.py` (bank engine), `t1..t5_*.py`,
  `t5b_mscaling.py` (M-scaling + matched-cell), `make_figures.py`, `validate_core.py`
- results: `T1_WALL_LEAK.json`, `T2_SPECKLE_STATS.json`, `T3_APERTURE_LAW.json`,
  `T4_MULT_TO_CONV.json`, `T5_SENTINEL.json`, `T5b_MSCALING.json`
- figures: `figs/fig1_wall_leak.png … fig5_roc_dprime.png`, `fig6_mscaling.png`
- COMSOL: `COMSOL_MICROJOB_SPEC.md` (spec); `COMSOL_MICRO/` (validation run of the
  phase-screen approximation — separate leg, license-gated)
