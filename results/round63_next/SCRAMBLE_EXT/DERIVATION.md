# SCRAMBLE_EXT — the beyond-band sentinel through a FULLY SCRAMBLING dynamic medium

**Date:** 2026-07-24 · **Line:** ROUND63 two-ledger / beyond-modulator-band sentinel.
**Predecessor:** R39 thin multiplicative screen (`docs/ROUND63_GPT_ROUND39_RULING_RAW.md`
§1–2) and the validated detector `results/round63_next/FOG_DMD_PROBE64/`
(DETECTION_ROC.md, DETECTION_MAP_CORRECTED.md).
**Scope:** derive the mean-wall / covariance-aperture / sentinel structure when the thin screen
is replaced by a **strongly scattering dynamic medium** (multimode fibre, thick diffuser, fully
developed speckle). Companion sim: `scramble_toy.py`, results `SCRAMBLE_RESULTS.{md,json}`.
Write-scope: this directory only. No git commit.

---

## 0. Why this is the decisive extension

R39 proved a clean *made-possible-at-all* separation for a **thin multiplicative screen**
`w_t`: the band-limited modulator is provably blind in its **first moment** to frequencies
outside its own band `B_p`, while the unmeasured dynamic screen shifts beyond-band scene
frequencies into a **bucket-covariance** channel with aperture `A_cov = B_w ⊕ B_p`. The
`FOG_DMD_PROBE64` detector then showed this is not just estimable but *detectable*: a 2 %
beyond-band anomaly is caught in seconds on a covariance channel the mean route provably cannot
access, with a four-way `mean / beyond-band-cov / medium-amplitude` specificity split.

A thin screen is the *weakest* possible complex medium — a single multiplicative layer. The
open question that decides the tier of the whole program is:

> Does the wall / aperture / sentinel structure survive when the medium **fully scrambles** the
> field — the multimode-fibre / thick-scattering / developed-speckle regime that covers essentially
> **all** dynamic complex media?

This note answers it. The short version, stated up front and honestly:

- **The mean-wall gets STRONGER.** Under full scrambling the mean speckle envelope is *flat*
  (memory-effect–limited), so the first-moment channel collapses to the scene's **DC term
  alone** — band `B_mean = B_ME`, which → `{0}` in the fully mixing limit. The wall is no longer
  `B_p`; it is essentially total.
- **The covariance aperture becomes GIGANTIC in support but COLLAPSES in rank.** The
  speckle-grain correlation kernel has a band `2·B_grain` reaching the whole grid (wavelength-scale
  grains) — spectrally the wall is fully pierced. **But** the scene enters the fully-scrambled
  covariance only through a **single nonnegative scalar** `Q(x) = xᵀG x` (grain-weighted energy).
  The mode-coverage curse does not merely return — under *complete* scrambling it collapses the
  entire spatial aperture to **rank one**. Reconstruction of the beyond-band pattern does **not**
  survive; only its *energy footprint* does.
- **The sentinel SURVIVES — and in a strictly stronger form.** Detection rides on `Q`. Any
  beyond-band change `δ` lifts `Q` by `2 xᵀGδ + δᵀGδ`, and the quadratic term `δᵀGδ =
  Σ_k Ĝ(k)|δ̂(k)|² > 0` for **every** nonzero `δ`, including one entirely beyond band, DC-free, and
  `G`-orthogonal to `x`. There is **no blind direction for detection**, whereas the mean channel is
  exactly blind to all of it. The specificity split (mean = in-band only, covariance = beyond-band)
  ports intact, and the temporal-lag structure carries `Q` coherently across lags so the
  `FOG_DMD_PROBE64` LR-on-lag-covariances detector applies unchanged.

So the flagship claim splits cleanly on the honest boundary: **beyond-band *imaging* is a
thin-/thick-screen (finite-memory) capability; beyond-band *detection* (the sentinel) is
universal — it survives complete scrambling.**

---

## 1. Forward model — structured illumination through a dynamic scrambler onto a static scene

We commit to the cleanest physically defensible model of *structured illumination through a
dynamic strongly-scattering medium onto a static scene with bucket detection*.

**Geometry.** A DMD imprints a real, nonnegative, **band-limited** code `p_i` (Fourier support
`B_p`, cut-off `k_p` — the hardware pixel-pitch wall) on the illumination field at the input/pupil
plane, `A_i(a) = √I₀ · p_i(a)`. The field propagates through a dynamic medium whose state at epoch
`t` is a random complex transmission operator `T_t`. At the scene plane this produces a **speckle
field**
```
E_{i,t}(r) = Σ_a T_t(r,a) A_i(a).                                     (1.1)
```
A **static** scene `x(r) ≥ 0` (reflectance/transmittance) modulates the local intensity, and a
single-pixel **bucket** collects the total:
```
b_{i,t} = Σ_r S_{i,t}(r) x(r) + shot,     S_{i,t}(r) = |E_{i,t}(r)|².  (1.2)
```
This is the honest bucket model: the speckle **intensity** pattern `S_{i,t}` produced by code `i`
through medium state `t` multiplies the scene pointwise and is summed by the bucket. `S` is the
scrambling-medium analogue of the thin screen's `P diag(w_t)`.

**Medium statistics (fully developed speckle).** "Fully developed speckle / fully scrambling" is
exactly the statement that `T_t` has zero-mean **circular complex Gaussian** entries, so for each
`(i,t)` the field `E_{i,t}(·)` is a zero-mean circular complex Gaussian field over the medium
ensemble. Its law is then fixed by the mutual coherence
```
J_{ij}(r,r') := ⟨ E_{i,t}(r) E_{j,t}(r')* ⟩_t.                        (1.3)
```
The medium's known second-order statistics are
```
⟨ T_t(r,a) T_t(r',b)* ⟩ = C(r,r') · Γ(a,b),                          (1.4)
```
with two physically distinct kernels:

- **Output grain kernel** `C(r,r') = C(r−r')` (stationary): the field-coherence width is the
  **speckle grain** `w_g ~ λ/(2·NA)` — wavelength-scale, hence a **sharp** kernel with an
  **enormous** spectral band `B_grain` (cut-off `k_grain ≫ k_p`). This is the object that plays the
  role of the thin-screen medium spectrum `K_w`, but with a vastly larger band.
- **Input coupling** `Γ(a,b)`: how input channels mix. **Full scrambling** ⟺ `Γ(a,b) = τ(a) δ_{ab}`
  (independent input channels, no residual angular correlation). A **finite angular memory effect**
  ⟺ `Γ` with finite correlation width / finite rank `r_ME`.

Substituting (1.4) into (1.3):
```
J_{ij}(r,r') = C(r−r') · Σ_{a,b} Γ(a,b) A_i(a) A_j(b)*.              (1.5)
```

**The fully-scrambling factorization.** With `Γ(a,b)=τ(a)δ_{ab}`,
```
J_{ij}(r,r') = C(r−r') · O_{ij},   O_{ij} = I₀ Σ_a τ(a) p_i(a) p_j(a).  (1.6)
```
`O` is the (real, symmetric, PSD) **weighted code Gram matrix**; `O_{ii}` is the total power code
`i` delivers. The spatial part `C(r−r')` and the code part `O_{ij}` **separate**. This one line is
the source of every result below.

> **Generator identity (used by the sim).** The joint field `{E_{i,t}(r)}` is zero-mean complex
> Gaussian with the **Kronecker** covariance `C(r−r') ⊗ O_{ij}`. Hence it is generated exactly by
> `E_t = B·(Z_t·O^{1/2})`, with `Z_t` iid unit circular-Gaussian `(N×M)`, `O^{1/2}` the code-cov
> root, and `B` the spatial grain operator (`BB† = C`). This is the transmission-matrix model
> `E_{i,t}=B G_t A_i` with `G_t` iid, and it is what `scramble_toy.py` simulates — no approximation.

---

## 2. Derivation 1 — the mean channel and the (stronger) wall

The mean speckle intensity is the coherence at zero lag:
```
μ_i(r) := ⟨S_{i,t}(r)⟩ = J_{ii}(r,r) = C(0) · O_{ii}.                (2.1)
```
`C(0)` is the (constant) peak coherence and `O_{ii}` is a scalar per code. **`μ_i(r)` is
independent of `r`** — the mean speckle envelope is **flat**. The mean bucket is therefore
```
E[b_{i,t}] = Σ_r μ_i(r) x(r) = C(0) O_{ii} · Σ_r x(r) = C(0) O_{ii} · x̂(0).  (2.2)
```

### Theorem 1 (mean wall of the scrambling channel).
Under (1.1)–(1.6) with a fully scrambling medium, the first-moment channel measures **only the
scene DC term** `x̂(0)=Σ_r x(r)`, scaled by the known per-code power `C(0)O_{ii}`. Its spatial band
is `B_mean = {0}`. Consequently, for **every** perturbation `h` with `ĥ(0)=0` (any non-DC change,
in-band or beyond-band),
```
E[b_i(x+h)] = E[b_i(x)],   I_mean h = 0.                             (2.3)
```

**Proof.** (2.2) depends on `x` only through `x̂(0)`; differentiating, `∂E[b_i]/∂x = C(0)O_{ii}·𝟙`,
a rank-one operator onto the DC mode. Any `h ⊥ 𝟙` is annihilated. ∎

**Memory-effect refinement.** With finite memory, `C` and `Γ` give a smooth envelope
`μ_i(r) = (h_ME * |p_i|²)(r)` and the wall relaxes to `B_mean = B_ME ∩ 2B_p`, where `B_ME` is the
angular-memory band. In the strong-scattering regime `B_ME ≪ B_p`, so the wall is still far tighter
than the thin-screen wall `B_p`. **The mean-wall is strictly stronger under scrambling**, degrading
continuously from `B_p` (thin screen) to `{0}` (complete scramble). This is the exact,
mechanism-level reason the empirical `1.000` beyond-band baseline is even *more* robust here: the
mean route loses not just the beyond-band frequencies but essentially all spatial structure.

---

## 3. Derivation 2 — the covariance channel, the aperture law, and the rank collapse

For a circular complex Gaussian field the intensity covariance obeys the **Siegert relation**
```
Cov_t[S_{i,t}(r), S_{j,t}(r')] = |J_{ij}(r,r')|².                    (3.1)
```
Hence the same-epoch bucket covariance is
```
Cov_t[b_{i,t}, b_{j,t}] = Σ_{r,r'} x(r) x(r') |J_{ij}(r,r')|².        (3.2)
```
This is the covariance analogue of R39's `C(x)=Φ²P diag(x)K_w diag(x)Pᵀ`, but with the medium
kernel now the **speckle grain**. Using the factorization (1.6),
```
|J_{ij}(r,r')|² = |C(r−r')|² · O_{ij}²,
```
so the entire same-epoch bucket covariance collapses to
```
┌─────────────────────────────────────────────────────────────────┐
│  Cov_t[b_{i,t}, b_{j,t}] = O_{ij}² · Q(x),                        │
│  Q(x) = Σ_{r,r'} x(r) x(r') |C(r−r')|² = xᵀ G x,   G = |C|²-Toepl. │   (3.3)
└─────────────────────────────────────────────────────────────────┘
```

### 3.1 The aperture law (support) — gigantic
`G(r,r') = |C(r−r')|²` is a stationary PSD kernel; in Fourier
```
Q(x) = Σ_k Ĝ(k) |x̂(k)|²,   Ĝ(k) = FT{|C|²}(k) = (m̂ ⋆ m̂)(k) ≥ 0,   (3.4)
```
where `m̂` is the grain pupil (`Ĉ`). `Ĝ` is the **autocorrelation of the grain spectrum**, supported
on `|k| ≤ 2 k_grain`. Because the grain is wavelength-scale, `2 k_grain` reaches essentially the
**whole grid**:
```
Covariance support aperture:  A_scramble = 2·B_grain ⊇ B_p ∪ B_beyond.  (3.5)
```
Spectrally the wall is **fully pierced** — the grain kernel plays the role of `K_w` with a *huge*
band, exactly as anticipated. `Ĝ(k) > 0` out to `2 k_grain`, so **every** scene frequency, in-band
or beyond, contributes to `Q`.

### 3.2 The rank collapse (Fisher) — the mode-coverage curse returns in full
The scene enters (3.3) **only** through the scalar `Q(x)`. The full `M×M` same-epoch covariance is
the outer object `Q(x)·(O∘O)` with `O∘O` **known** (Hadamard square of the code Gram). Measuring all
`M(M+1)/2` entries yields `M(M+1)/2` noisy looks at the **same single number** `Q(x)`. Therefore the
covariance Fisher information for the scene is **rank one**:
```
I_xx^{cov} = κ · (∇Q)(∇Q)ᵀ = 4κ · (Gx)(Gx)ᵀ,   κ = precision of Q̂.   (3.6)
```
Only the one-dimensional scene direction `Gx` is identifiable from the covariance channel under
complete scrambling.

> **This is the honest heart of the extension.** In R39 the thin-screen covariance `C(x)` had a
> *rich* `M×M` structure whose entries carried **distinct** Fourier components of `x` across the
> aperture `B_w⊕B_p` (mode-coverage curse, but many recoverable modes). Complete scrambling
> **factorises** the code and space parts (`Γ` diagonal), collapsing all `M(M+1)/2` entries to a
> **single scalar**. The mode-coverage curse does not merely return — under full scrambling it is
> total: the recoverable scene aperture for *reconstruction* is **rank 1**.

### 3.3 The interpolation — how much imaging survives with a residual memory effect
For a memory effect of rank `r_ME` (i.e. `Γ = Σ_{m≤r_ME} λ_m u_m u_mᵀ`), (1.5) gives
```
J_{ij}(r,r') = Σ_m λ_m C(r−r') (u_mᵀA_i)(u_mᵀA_j)*,
|J_{ij}|² = Σ_{m,m'} λ_m λ_{m'} |C(r−r')|² (u_mᵀA_i)(u_mᵀA_j)*(u_{m'}ᵀA_i)*(u_{m'}ᵀA_j),
```
so the covariance becomes a rank-`O(r_ME²)` combination of quadratic scene functionals
`Q_{mm'}(x) = Σ_{r,r'} x(r)x(r')|C(r-r')|²(…)` — a **larger** recoverable scene subspace whose
dimension **grows with the number of memory-effect / open channels**. The recoverable aperture thus
**interpolates monotonically**:
```
   r_ME = 1  (complete scramble)  →  rank-1 scalar Q(x)          [imaging dead, detection alive]
   r_ME large (≈ thin phase screen) →  R39 rich aperture B_w⊕B_p [imaging alive]
```
**Beyond-band imaging survives exactly to the extent the medium retains memory-effect modes.** This
is the precise, testable boundary of the flagship claim, and it is physically correct: imaging
through scattering *is* the exploitation of the memory effect; when the memory effect is gone, so is
the image — but, as §4 shows, **not** the sentinel.

---

## 4. Derivation 3 — the sentinel (detection) survives, and is stronger

Static scene `x`; a **sub-resolution beyond-band change** `δ` with `‖δ‖=ε‖x‖`, spectral support
outside `B_mean` (take `δ̂(0)=0` and support beyond `B_p`). H0: scene `x`. H1: scene `x+δ`.

**Mean channel (blind).** By Theorem 1, `Δm_i = C(0)O_{ii}·δ̂(0) = 0`. The first-moment detector has
**zero** response to `δ` — the wall.

**Covariance channel (responds, unconditionally).** From (3.3) the observable scene functional
shifts by
```
ΔQ = Q(x+δ) − Q(x) = 2 xᵀG δ + δᵀG δ
   = 2 Σ_k Ĝ(k) x̂(k)* δ̂(k)  +  Σ_k Ĝ(k) |δ̂(k)|².                    (4.1)
```
Two terms, both live because `Ĝ(k)>0` out to `2k_grain ⊇` the beyond-band support of `δ`:

- **Quadratic term** `δᵀGδ = Σ_k Ĝ(k)|δ̂(k)|² > 0` for **any** nonzero `δ` — including one entirely
  beyond band, DC-free, and `G`-orthogonal to `x`. **This cannot be nulled.** It is a guaranteed,
  sign-definite, positive lift of `Q`. There is **no blind direction for detection.**
- **Coherent term** `2xᵀGδ` — generically nonzero, dominant (∝ ε) when `δ` is not `G`-orthogonal to
  `x`.

So the beyond-band change is detected **entirely** through the covariance channel while the mean
channel is exactly blind: the R39/`FOG_DMD_PROBE64` **specificity split survives intact**, and the
`δᵀGδ>0` floor makes it *stronger* than the thin-screen sentinel (which had estimable but not
sign-definite responses).

### 4.1 Detection Fisher and the d′ law
`Q` is estimated from the sample covariance across `T_eff` independent medium realisations. The
efficient (Fisher) estimator whitens the `M×M` covariance by `Σ₀⁻¹` before projecting onto the
signal direction `∂Σ/∂Q = O∘O` — the matched-score statistic `t = ⟨Σ₀⁻¹(O∘O)Σ₀⁻¹, Ĉov⟩`, identical
to the FOG `W=V₀⁻¹dV V₀⁻¹` detector. In the medium-dominated regime (shot below fluctuation,
`n≳10⁴`, cf. FOG §4.3):
```
d′  ≈  (ΔQ / Q) · √(M_eff · T_eff / 2) · [ r/(1+r) ],   r = n σ_S² · geom.   (4.2)
```
with `M_eff` the effective number of independent covariance looks and `r` the covariance visibility.

> **`M_eff` is capped by the in-band code DOF — a second face of the rank collapse.** Although
> `Σ = Q·(O∘O)+shot` is a *full-rank* `M×M` matrix (so whitening gives `M_eff` far above the naive
> projection's `M_eff≈1`), the informative code directions live in the band-limited illumination
> subspace of dimension `(k_p+1)²`. The numerics give `M_eff ≈ 13` for all duels, matching
> `(k_p+1)²−1 = 15`. Adding codes beyond the in-band DOF buys nothing: **independent medium banks,
> not codes or photons, are the currency.** This is the sample-complexity price of full scrambling.

The scene-signal ratio is
```
ΔQ/Q = (2xᵀGδ + δᵀGδ)/(xᵀGx).
 · δ ⟂_G x (worst case):  ΔQ/Q = δᵀGδ/xᵀGx ≈ ε²   →  d′ ~ ε² √(M_eff T_eff/2).
 · generic δ:             ΔQ/Q ≈ O(ε)             →  d′ ~ ε  √(M_eff T_eff/2).
```
The **guaranteed** worst-case (orthogonal, DC-free, beyond-band) detection scales as `ε²√T_eff` and
is recovered by aggregation and banks — never zero. This mirrors the FOG rate law: **independent
medium realisations, not photons, are the currency** (past `n≈10⁴`, buy banks not flux).

### 4.2 Temporal dimension — the lag ledger ports unchanged
Let the medium decorrelate with field-correlation `ρ_E(h)` over lag `h` (correlation time `τ_c`).
Within `τ_c` the speckle is frozen (one **bank** = all codes under one realisation); across banks,
independent. The lag-`h` bucket covariance is, by the temporal Siegert relation,
```
Cov_t[b_{i,t}, b_{j,t+h}] = O_{ij}² · Q(x) · ρ_S(h),   ρ_S(h) = |ρ_E(h)|².   (4.3)
```
The **scene functional `Q(x)` is identical across all lags**; only the profile `ρ_S(h)` — a pure
**medium-law fingerprint** — varies. Therefore:

- The `FOG_DMD_PROBE64` **LR-on-lag-covariances** detector applies verbatim: the scene signal `Q`
  modulates every lag coherently, so stacking lags `0…L` boosts `M_eff`.
- A medium **τ change** moves `ρ_S(h)` but leaves `Q` fixed → it is **separable** from a scene
  change on the lag-decay axis. A medium **amplitude (σ_S) change** rescales `Q`'s prefactor and is
  separable on the amplitude axis. This reproduces the validated **four-way discriminator**
  (`mean / beyond-band-cov / medium-amplitude / medium-τ`) — now with the mean axis even cleaner
  (pure DC) and the beyond-band axis carrying the sign-definite `δᵀGδ` floor.

The quasi-static stepped-diffuser acquisition of FOG (`T_eff` = independent banks, no wall-clock τ
penalty) carries over directly: detection time = `T_det · t_bank`.

---

## 5. What the toy sim must confirm (and where it can break the claim)

`scramble_toy.py` simulates (1.1)–(1.6) exactly via the Kronecker generator, on a 32×32 grid with
complex fields, a fully-scrambling medium, band-limited codes `B_p`, and a grain band `B_grain ≫ B_p`.
It is instructed to try to **falsify** the following, numerically:

1. **Mean-band collapse (Thm 1).** Empirical `μ_i(r)` flat (spatial contrast → 0 up to `1/√epochs`);
   `E[b_i]` responds only to `x̂(0)`; **null response** of the mean channel to a beyond-band, DC-free
   `δ` (mean-detector AUC ≈ 0.5).
2. **Covariance = `Q(x)` (Eq. 3.3).** Empirical `Var_t[b_i] = O_{ii}²Q(x)+shot` and
   `Cov_t[b_i,b_j] = O_{ij}²Q(x)` — i.e. the `M×M` covariance is rank-1 in matrix space up to shot;
   `ΔQ` predicted by (4.1) matches the measured shift.
3. **The sentinel duel.** H0 vs H1 (beyond-band sub-resolution `δ`), LR on the (lag-)covariance,
   `d′` and ROC at 2–3 bank budgets; confirm the covariance detector fires while the mean detector
   is at chance.

**Where it could break (honest failure modes to watch in the sim):**
- If the grain is *not* much finer than a scene pixel, `Ĝ(k)` may roll off before the beyond-band
  support of `δ`, weakening `δᵀGδ` — the sentinel would degrade for coarse speckle.
- The `x⊗x` coupling means two different `δ`'s with equal `‖δ‖_G` and equal low-`Q`-projection are
  **indistinguishable** in the aggregate — detection yes, *identification/reconstruction* no. This
  is the rank-1 curse of §3.2 and is expected, not a bug.
- With shot dominating (`n` small) the `δ²`-floor sinks below the shot variance and worst-case
  (orthogonal) detection stalls at finite `T_eff` — a finite-data CRB effect, not a physical floor
  (cf. R39 Thm 3). The sim reports the required `T_eff`.

The numbers, not this prose, decide. See `SCRAMBLE_RESULTS.md`.
