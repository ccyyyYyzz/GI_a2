# R44 — Wave-optics divergence: the finite-aperture wall, physical jet crossover, and the next paper

**Reference brief:** [`docs/ROUND63_GPT_ROUND44_QUESTION.md`](https://github.com/ccyyyYyzz/GI_a2/blob/60d5fc7/docs/ROUND63_GPT_ROUND44_QUESTION.md) at [`60d5fc7`](https://github.com/ccyyyYyzz/GI_a2/commit/60d5fc7d00e60be0c0e12598c30540fabad748ef).  
**Wave twin:** [`results/round63_next/WAVE_TWIN/`](https://github.com/ccyyyYyzz/GI_a2/tree/60d5fc7/results/round63_next/WAVE_TWIN).  
**Letter quarantine:** the twin contributes only the frozen transfer paragraph in §3. No T1–T5 result becomes a Letter endpoint, figure, or theorem test.

## Executive synthesis

The wave twin does not break Rank–Jet Separation. It reveals the exact distinction between an **algebraic wall** and a **finite-aperture physical wall**.

1. **Free-space propagation does not create new spatial support.** For a truly band-limited complementary field, the signed-intensity support is a difference set and remains compact at every propagation distance.
2. **The present DMD is not truly band-limited.** Its finite active aperture, sample-and-hold lattice, mirror aperture, finite object window, and macro-pixel readout produce sinc/Dirichlet tails and aliases. Those tails are generically nonzero on every open spectral annulus. Propagation changes their phase and destroys accidental cancellation; it does not create them.
3. Therefore **placing the witness beyond `2 k_p` removes the dominant ideal intensity core but does not, by itself, restore an exact physical wall** in the current bench. An exact support wall requires an explicit compact Fourier pupil plus alias/guard-band control, or a witness constructed in the measured physical mean operator's null space.
4. Near screen–object contact, phase-to-intensity fluctuation is linear in `z_2`; bucket covariance is `O(z_2^2)` and covariance Chernoff information is `O(z_2^4)`. A pre-existing mean leak is `O(z_2^0)`. This gives a fourth-root channel crossover and explains why the mean wins near contact while covariance wins after conversion develops.
5. The physically realized jet order is selected by the **lowest nonzero channel order**. A tiny first-order mean leak changes the strict local order of an ideal curvature-rescued direction from `m=2` to `m=1`; the `m=2` class is restored only after the physical leak tangent is removed by an optical wall or nuisance projection.
6. There is an exact Gaussian **mean–covariance Fisher decomposition**, but no conservation law of the form “mean information lost equals covariance information gained.” A weak-phase Fresnel quadrature identity exists before bucket compression; it does not survive as a constant Fisher sum.

The next paper should be the **memory-effect Rank–Jet phase diagram**: the end-to-end wave-optics family connecting the thin-screen and complete-scrambling endpoints, with an exact coherence-separation-rank bound and a measured transition from image recovery to change-only observability.

---

# 1. Exact leak-support derivation

## 1.1 Complementary fields: the intensity square cancels exactly

Write the physical complementary DMD amplitudes as

```text
a_+(r) = [c(r)+s(r)]/2,
a_-(r) = [c(r)-s(r)]/2,
```

where `c` is the common illumination/aperture field and `s` is the signed code. Let `L_z` be the linear optical map from the DMD to the plane at which the signed intensity is formed. Define

```text
U_z = L_z c,
S_z = L_z s.
```

Then

```text
E_+ = (U_z+S_z)/2,
E_- = (U_z-S_z)/2,
```

and the physical signed intensity is exactly

```text
D_z(r)
 = |E_+(r)|^2-|E_-(r)|^2
 = Re[U_z(r) S_z(r)^*].                                      (1.1)
```

The `|S_z|^2` intensity-autocorrelation term cancels. The leak is a **common-field/code-field heterodyne cross term**, not an uncancelled square of the signed code.

Taking a Fourier transform gives

```text
Dhat_z(q)
 = 1/2 ∫ [ Uhat_z(k+q) Shat_z(k)^*
          +Shat_z(k+q) Uhat_z(k)^* ] dk.                       (1.2)
```

Hence

```text
supp(Dhat_z)
 ⊆ [B_U - B_S] ∪ [B_S - B_U],                                 (1.3)
```

where `B_U=supp(Uhat_z)` and `B_S=supp(Shat_z)`.

## Theorem 1 — propagation-invariant difference-set wall

Let `L_z` be shift invariant with transfer function `H_z(k)`. If the input common and signed fields have compact Fourier supports `B_c` and `B_s`, then for every propagation distance

```text
supp(Dhat_z)
 ⊆ [(B_c∩B_H)-(B_s∩B_H)]
    ∪[(B_s∩B_H)-(B_c∩B_H)].                                   (1.4)
```

Free-space Fresnel or angular-spectrum propagation changes only the coefficients through the phase factor

```text
H_z(k+q) H_z(k)^*;
```

it cannot enlarge the support.

### Important special cases

1. **Infinite periodic plane-wave common mode.** `B_c={0}`, `B_s=B_p`. Then

   ```text
   supp(Dhat_z) ⊆ B_p                                         (1.5)
   ```

   for every `z`. Propagation alone cannot leak beyond the code band.

2. **Both complementary fields hard-pupil limited to a symmetric relay band `B_R`.** Then

   ```text
   supp(Dhat_z) ⊆ B_R-B_R.                                    (1.6)
   ```

   For a circular or Chebyshev cutoff `k_R`, this is the familiar `2 k_R` intensity band.

3. **Both fields merely described as “code-band-limited” but their common aperture is not.** The relevant bound is not `2 k_p`; it is the difference set of the *physical field supports*. The common-mode aperture is load bearing.

This corrects one phrase in the current T1 mechanism note: propagation does not manufacture a `2 k_p` tail. It changes the Fresnel phases of field components already supplied by the finite aperture and pixel lattice, thereby activating or amplifying their signed-intensity cross terms.

## 1.2 Why the present DMD has no exact open spectral wall

For a finite mirror array with mirror footprint `m(r)`, pitch `p`, and coefficients `s_n`, the signed field immediately after the DMD has the form

```text
s_D(r) = Σ_{n∈A} s_n m(r-np),
```

so

```text
shat_D(k) = mhat(k) Σ_{n∈A} s_n exp(-i k·np).                  (1.7)
```

The second factor is a finite trigonometric polynomial. The mirror transform `mhat` is a separable sinc. Neither is compactly supported. For an infinite periodic lattice, the spectrum is a comb of replicated code bands weighted by the mirror sinc; the finite active array broadens each order by a Dirichlet/sinc envelope. The common field has the same structure with `s_n=1`.

Therefore the current full-field supports satisfy, generically,

```text
B_U = B_S = R^2                                                (1.8)
```

up to isolated sinc zeros. Equation (1.3) then supplies no nontrivial exact annular exclusion.

There is a second alias path. Macro-cell integration followed by 64-grid sampling gives

```text
Jhat_macro(κ)
 = Σ_{g∈Λ*} hhat_cell(κ+g) Dhat_z(κ+g),                        (1.9)
```

so high diffraction orders can fold into the macro-grid baseband, weighted by the cell-box sinc. The “core plus pedestal” picture is therefore:

- a dominant zero-order core whose ideal support is at most `B_R-B_R` (often approximately `2 k_p`);
- a noncompact finite-aperture/mirror pedestal;
- lattice-difference orders and their macro-sampling aliases.

The T1 floor of `6.8×10^-3` at `z_1=0` is exactly the signature expected from this pedestal and finite-window sampling. The monotone increase to `8.9×10^-2` at `20 mm` is phase reweighting of existing tails, not support creation.

## Proposition 1 — finite-aperture no-open-wall principle

A nonzero field confined to a finite DMD aperture has an entire Fourier transform. By the Paley–Wiener principle, if that transform vanished on an open spectral set it would vanish identically. Thus a generic finite DMD field cannot possess an exact open Fourier stop band without a separate Fourier-plane filter.

The same warning applies to a spatially localized object change: a compactly supported defect is not strictly band limited. An exact “beyond-band witness” is an ideal periodic/apodized mode or a vector in the actual discrete mean-operator null space, not an arbitrary local feature described only by its nominal size.

## 1.3 Ruling on the `>2 k_p` witness

> **Moving the witness outside `2 k_p` is strongly recommended, but it does not restore an exact wall in the present unfiltered DMD geometry.**

It will remove the dominant ideal zero-order intensity core and should sharply lower the mean leak. The residual sinc/Dirichlet pedestal and aliases remain nonzero. Exact support exclusion is obtained only if all of the following are enforced:

1. a hard Fourier relay pupil with compact cutoff `B_R`;
2. witness support disjoint from `B_R-B_R`;
3. sampling above the intensity-band Nyquist, with no macro-order alias;
4. a guard/apodization region so the finite object window does not reintroduce low-frequency tails;
5. or, most robustly, construction of the test change in the measured physical mean Jacobian's null space.

For `k_p=5`, a practical twin/bench pair is:

- current witness `k=6…9`: tests the region where ideal intensity leakage is expected;
- guarded witness `k=12…16`: outside the `2 k_p=10` core, while remaining inside the measured covariance aperture for the 8–48 µm grain cases;
- a narrow `k=11…12` witness for the coarsest 96 µm case, whose measured covariance edge is about `13`.

A `4f` pupil with cutoff near `k_R≈5` and a witness beginning at `k≥12` gives the cleanest wall experiment.

## 1.4 Decisive ≤60-GPU-minute twin test

Run a deterministic shell-resolved wall scan with `M=8–16`, no diffuser (`z_2=0`), and `z_1∈{0,1,5,10,20} mm`. For each Chebyshev shell `k=1,…,31`, report

```text
L_k = ||J U_k||_F / ||J U_[0,k_p]||_F.                         (1.10)
```

Use seven nested optical variants:

1. periodic ideal macro-code, constant common field;
2. finite active aperture only;
3. macro sample-and-hold;
4. mirror sinc envelope;
5. full current twin;
6. full twin plus hard pupil `k_R∈{5,6,8}`;
7. hard pupil plus guarded/apodized object ROI.

### Frozen predictions

- Variant 1: machine-zero response for `k>k_p` at every `z_1`.
- Hard-pupil variants: machine-zero response for `k>2k_R`, subject to the ROI/sampling guard.
- Full current twin: a sharp drop beyond the core, but a nonzero tail on every shell; exact zeros only at unstable isolated comb/sinc zeros.
- Fresnel propagation changes amplitudes inside the allowed support but never moves the ideal support edge.

A full-twin value `L_k>10^-10` on any shell above `2k_p` falsifies the “frequency placement alone restores an exact wall” claim. A pupil-limited value below `10^-12` beyond `2k_R` validates the engineered exact wall.

---

# 2. Two-channel `(z_1,z_2)` phase diagram

## 2.1 Near-contact phase-to-intensity expansion

At the screen plane write one code field as

```text
ψ_t(r,0)=A(r) exp{i[θ(r)+φ_t(r)]}.
```

Under paraxial propagation the transport-of-intensity equation gives

```text
∂_z I_t|_0
 = -(1/k_0) ∇·[I_0 ∇(θ+φ_t)].                                  (2.1)
```

Separating deterministic and zero-mean random parts,

```text
I_t(r,z_2)
 = I_0(r)+z_2 L_0(r)+z_2 ξ_t(r)+O(z_2^2),
E[ξ_t]=0.                                                       (2.2)
```

For a scene `x`, the random bucket component is

```text
b_t-Eb_t = z_2 <x,ξ_t> + O(z_2^2),                              (2.3)
```

so

```text
Cov(b_i,b_j)
 = z_2^2 K_ij(x)+O(z_2^3).                                     (2.4)
```

This is the exact local reason for T2 contact degeneracy: a phase screen at `z_2=0` changes phase but not intensity.

For a scene path `x_ε=x+εδ`, the covariance change has the quadratic expansion

```text
ΔΣ
 = z_2^2 [ ε B_δ + ε^2 C_δ ]
   +O(z_2^3, ε^3 z_2^2),                                      (2.5)
```

where `B_δ=DK_x[δ]` and `C_δ=K(δ)`.

## 2.2 Mean and covariance Chernoff rates

Let the physical mean leak along `δ` be

```text
Δμ = ε ℓ(z_1,z_2;δ)+O(ε^2),                                   (2.6)
```

and let `R` be the shot/read covariance at contact. In the shot-dominated near-contact regime, the per-bank Bhattacharyya/Chernoff contributions are

```text
C_mean
 = (ε^2/8) ℓ^T R^-1 ℓ + o(ε^2),                               (2.7)

C_cov
 = (z_2^4/16)
   || R^-1/2 (ε B_δ+ε^2 C_δ) R^-1/2 ||_F^2
   +o(z_2^4).                                                   (2.8)
```

Thus:

- mean leakage is `O(z_2^0)` at contact;
- covariance distinguishability is `O(z_2^4)`;
- covariance `d'/sqrt(T)` is `O(z_2^2)`.

As covariance becomes medium dominated, the `z_2^4` growth saturates through the usual visibility factors `a/(1+a)`; it does not grow indefinitely.

## 2.3 Crossover law

Define

```text
L = ℓ^T R^-1 ℓ,
K_1 = ||R^-1/2 B_δ R^-1/2||_F^2,
K_2 = ||R^-1/2 C_δ R^-1/2||_F^2.
```

For a generic covariance direction (`B_δ≠0`), equal mean and covariance exponents occur at

```text
z_2,* = (2L/K_1)^(1/4).                                       (2.9)
```

For a first-order-orthogonal, curvature-rescued covariance direction (`B_δ=0`, `C_δ≠0`),

```text
z_2,*(ε)
 = [2L/(ε^2 K_2)]^(1/4)
 = ε^-1/2 (2L/K_2)^(1/4).                                     (2.10)
```

Consequences:

1. Larger `z_1` increases the leak coefficient `L` and moves the crossover outward.
2. Smaller changes move the curvature-channel crossover outward as `ε^-1/2`.
3. Near contact the mean must dominate whenever `L>0`; covariance can dominate only after phase-to-intensity conversion develops.
4. The T5 result brackets the present `M=32` crossover between `z_2=1` and `5 mm`. The absolute location remains provisional until the queued `M=128` scaling run lands, but the fourth-root structure does not depend on `M`.

## 2.4 Physical jet order

The statistical theorem classifies contact order in change amplitude `ε`. Real hardware adds a second contact variable, `z_2`. Near contact the channel-resolved bi-orders are

```text
mean leak, generic δ:           (m_ε,m_z)=(1,0)
covariance, generic δ:          (1,2)
covariance, curvature δ:        (2,2)                            (2.11)
```

because Chernoff distance scales as `ε^(2m_ε) z_2^(2m_z)`.

For a nominal curvature-rescued direction, the full physical law has

```text
C_phys(ε)
 = A(z_1,z_2) ε^2 + B(z_2) ε^4 + o(ε^4).                       (2.12)
```

If `A>0`, the strict local physical jet order is `m=1`: the smallest nonzero contact order wins. The finite-amplitude effective order is

```text
m_eff(ε)
 = [A+2B ε^2]/[A+B ε^2],                                      (2.13)
```

which crosses continuously from `1` to `2` near

```text
ε_cross = sqrt(A/B).                                           (2.14)
```

This is the physical version of the mixed-direction crossover already validated in `JET_TEST`.

> **A tiny mean leak does not kill curvature-rescued covariance detection, but it does change its asymptotic class unless the leak is removed from the quotient.**

The `m=2` class is physically restored by either:

- an engineered optical wall;
- a witness in the measured mean null space;
- or nuisance projection of the calibrated physical mean-leak tangent before the covariance score is formed.

This is a central theorem target for the next paper, not a write-back to the current Letter.

## 2.5 Phase regions

| Region | Physical condition | Dominant channel | Information statement |
|---|---|---|---|
| **Contact wedge** | `z_2→0` | mean leak | covariance absent as `z_2^4`; physical order selected by leak |
| **Conversion band** | `z_2≈z_2,*` | joint | mean and covariance scores both needed; best depth sensitivity |
| **Developed thin-screen band** | phase converted, memory retained | covariance | rich covariance aperture; multiplicative approximation only in the sub-mm/weak-mm range measured by T4 |
| **Scrambling band** | memory exhausted | covariance functional `Q` | reconstruction rank collapses toward one, but quotient jets remain finite after anchoring |

The present developed-screen thin multiplicative interval is only about `0.3 mm`; the weak-screen interval is roughly `1–2 mm`. The `z_2=5 mm` T5 point is already a scrambling/convolutive transfer point, not a validation of the thin multiplicative law.

## 2.6 Is there an exact sum rule?

### Exact statistical identity: yes, but it is additive—not conservative

For a Gaussian record `Y~N(μ_θ,Σ_θ)`, the mean score is linear in the centered record and the covariance score is centered quadratic. Their covariance is zero because Gaussian third central moments vanish. Therefore

```text
I_θθ
 = (∂θμ)^T Σ^-1 (∂θμ)
   +1/2 tr(Σ^-1 ∂θΣ Σ^-1 ∂θΣ).                                (2.15)
```

This is an exact **mean–covariance Fisher Pythagoras**. With nuisances, the same statement holds in the direct-sum score Hilbert space after one joint projection; it remains separately additive only when the nuisance tangent also respects the mean/covariance split.

### Conservation: no

There is no exact law

```text
I_mean(z_1,z_2)+I_cov(z_1,z_2)=constant.                        (2.16)
```

Unitary field propagation conserves optical energy, not object-weighted bucket Fisher information after finite aperture, shot noise, spatial integration, and nuisance profiling. Mean leakage gained is not covariance information lost.

### Weak-phase optical quadrature: a limited precursor

For a weak complex modulation `a+iφ`, coherent Fresnel propagation gives the contrast-transfer pair

```text
H_a(f,z)=2 cos χ,
H_φ(f,z)=2 sin χ,
χ=π λ z |f|^2,                                                  (2.17)
```

so

```text
|H_a|^2+|H_φ|^2=4.                                             (2.18)
```

This exact per-frequency quadrature identity explains phase-to-intensity conversion and the `z_2^2` fluctuation power near contact. Bucket integration, finite support, and noise prevent it from becoming an information conservation law.

---

# 3. Frozen Letter transfer paragraph

Use the following four sentences as the **only** manuscript write-back from `WAVE_TWIN`:

> **The exact mean-channel wall applies to a compactly band-limited signed illumination field; a physical complementary DMD pair has finite-aperture and pixel-lattice tails, so a scalar diffraction twin found residual mean leakage of about `0.7%` in a near-conjugate relay, increasing to about `9%` over a `20 mm` unrelayed DMD–screen gap. Near-conjugate imaging and apodization suppress this leakage, and placing the change beyond twice the relay passband removes the ideal intensity core, but an exact physical wall additionally requires a hard Fourier pupil and guarded sampling because the sinc/comb pedestal is noncompact. The declared multiplicative thin-screen model corresponds to screen–object spacings of order `≤0.3 mm` for developed speckle and roughly `1–2 mm` for a weak phase screen; beyond that the system crosses toward the scrambling regime, where Rank–Jet Separation survives although reconstruction rank and detection rate change. Diffraction also couples grain size and contrast, so neither the sealed grid's coarse-grain cells nor its absolute latency values are claimed to map one-to-one to a particular bench geometry.**

This paragraph is deliberately stronger than “the bench will work” and weaker than a new result section. It states the physical validity boundary without importing quarantined T1–T5 endpoints into the Letter's evidentiary chain.

---

# 4. Next-paper ruling

## 4.1 Verdict

> **The next paper should be the wave-optical memory-effect Rank–Jet phase diagram.**

Working title:

> **From Image to Sentinel across the Optical Memory Effect**

More theoretical alternative:

> **A Rank–Jet Phase Diagram for Scattered-Light Sensing**

This is the missing connecting family between the Letter's two exact endpoints. It can contain the T1 leak, T2 coupling, T3 aperture, T4 regime boundary, and T5 crossover without becoming an anthology because every measurement is plotted against one physical disorder coordinate.

## 4.2 Exact theorem the data should demand

Define the **coherence separation rank** `R` by a Schmidt decomposition of the code–space mutual-coherence tensor:

```text
J_ij(r,r')
 = Σ_{α=1}^R O_ij^(α) C_α(r,r').                               (4.1)
```

By the Siegert relation,

```text
Cov(b_i,b_j)
 = Σ_{α,β} O_ij^(α) O_ij^(β)* Q_αβ(x),                         (4.2)

Q_αβ(x)
 = ∬ x(r)x(r') C_α(r,r') C_β(r,r')* dr dr'.                    (4.3)
```

The marginal law therefore depends on the scene through at most `R^2` real quadratic functionals, giving the exact bound

```text
rank(I_x^cov) ≤ min(N,R^2),                                    (4.4)
```

with `R(R+1)/2` in the real-symmetric specialization.

- Complete scrambling is the separable endpoint `R=1`, hence rank one.
- A thin/memory-preserving screen has larger `R` and a richer scene Fisher space.
- If a positive-definite combination of the quadratic Hessians survives nuisance profiling, local quotient-jet order remains at most two even while reconstruction rank collapses.

This is the **Coherence-Rank / Jet-Order Separation Theorem**. The next paper should measure `R`, reconstruction rank, and jet order on the same physical `z_2` continuum.

## 4.3 Campaign sketch

### Physical sweeps

- `z_2`: dense logarithmic/linear grid from contact through `20 mm`, with extra points below `0.5 mm`;
- `z_1`: `{0,1,5,10} mm` or relay-equivalent pupil settings;
- screen RMS phase: `{0.3,1,2π}`;
- screen correlation length: `{5,15,50} µm`;
- illumination aperture: at least four diameters to expose the grain/contrast/étendue coupling;
- codes: `M={32,64,128}` with the queued subset scaling completed;
- photons: `10^4` primary, one low and one high audit;
- changes: generic, `G`-orthogonal, monotone-additive, and physical operator-null witnesses.

### Per-cell objects

1. physical mean-transfer support and pedestal;
2. mutual-coherence separation rank `R` by randomized SVD/sketching;
3. covariance-Jacobian eigenvalues and effective reconstruction rank;
4. generic and curvature quotient-jet orders;
5. mean and covariance Chernoff rates, crossover `z_2,*`, and CUSUM latency;
6. speckle grain, contrast, and effective grain count;
7. amplitude/timescale attribution leakage.

### One hero figure

A single horizontal disorder axis (`z_2` on the lower axis, measured coherence separation rank on the upper axis) with three aligned rows:

1. **support:** mean band contracts; covariance support expands;
2. **rank:** covariance scene rank rises/peaks/collapses toward one;
3. **jet:** generic/curvature orders and mean–covariance dominance remain finite, with the depth crossover marked.

The thin-screen and complete-scramble endpoints from the Letter appear as anchored limits, not separate acts.

### Compute budget

Use a staged hard cap of about `48 L4-equivalent GPU-hours`:

- `6 h` deterministic support/TIE/pupil sweep;
- `10 h` coherence-rank and Jacobian sketches;
- `24 h` selected-cell Monte Carlo and latency curves;
- `8 h` independent reruns and dtype/sampling audits.

A lower-resolution preflight must reproduce the rank bound and crossover before the 2048-grid Monte Carlo is released.

## 4.4 Ranking against the alternatives

1. **Memory-effect Rank–Jet family — first.** Highest theorem content, connects the Letter's endpoints, and turns all T1–T5 corrections into one physical phase diagram.
2. **Depth localization by channel crossover — second, as a corollary/application of the same paper.** Near contact `J_cov/J_mean∝z_2^4`; a calibrated ratio can range the screen/object separation or classify the scattering regime.
3. **Physical specificity recalibration — third.** Necessary for a bench and useful as a methods section, but not a new conceptual paper by itself.
4. **Leak-as-asset two-channel sentinel — fourth.** Potentially useful engineering, but weaker and more prior-art-adjacent than the rank–jet family. Fold it into the crossover/ranging section rather than lead with it.

---

# 5. Ranked innovations

Ranking is by **surprise × rigor × reach**, each scored from 1 to 5. All tests fit within `≤60 GPU-min` on the existing engines.

## 1. Finite-Aperture Wall Impossibility / Pupil-Hardened Wall

**Score:** `5 × 5 × 5 = 125`

**Mechanism.** A finite DMD aperture cannot have an exact open spectral stop band; the apparent wall is a compact-support idealization. A Fourier-plane pupil changes the mathematical class and restores a difference-set wall.

**Exact statement sketch.** Equations (1.3)–(1.8) plus Paley–Wiener: finite spatial support gives noncompact spectral support; hard pupil support `B_R` gives signed-intensity support `B_R-B_R` exactly.

**Decisive test.** The seven-variant shell sweep in §1.4. Require full-current tails above numerical zero and hard-pupil zero beyond `2k_R`.

**Nearest prior art.** Classical Paley–Wiener theory, coherent optical transfer functions, and DMD diffraction. The open square is using the theorem to distinguish an algebraic specificity wall from a physically realizable one in a complementary bucket detector.

**Potential named result:** **Finite-Aperture Wall Impossibility.**

## 2. Physical Jet Transmutation by Channel Leakage

**Score:** `5 × 5 × 5 = 125`

**Mechanism.** A weak first-order hardware leak adds `A ε^2` to an ideal curvature law `B ε^4`; the smallest nonzero order wins locally, while an efficient projection can restore the higher-order class.

**Exact statement sketch.** Equation (2.12) gives `m=1` for every `A>0`; Eq. (2.13) gives the finite-amplitude order and Eq. (2.14) the crossover. With `B∝z_2^4`, the crossover surface is a predictable function of `z_1,z_2,ε`.

**Decisive test.** On four `(z_1,z_2)` cells, sweep `ε` over two decades for a `G`-orthogonal direction. Fit the Chernoff slope before and after projecting out the measured mean-leak tangent. Prediction: unprojected slope tends to `2`; projected slope returns to `4`.

**Nearest prior art.** Weak-identification/singular-model contact orders and optical heterodyne leakage. No known wave-optics experiment treats hardware contamination as a change in quotient-jet universality class.

**Potential named result:** **Jet Transmutation.**

## 3. Coherence-Separation Rank Bound

**Score:** `5 × 5 × 5 = 125`

**Mechanism.** The number of separable code–space coherence modes, not nominal spectral support, bounds reconstruction dimension; intensity covariance squares that rank, while local curvature can remain complete.

**Exact statement sketch.** Equations (4.1)–(4.4): Schmidt rank `R` of the mutual-coherence tensor implies scene Fisher rank at most `R^2`; a positive quadratic Hessian combination guarantees quotient-jet order at most two.

**Decisive test.** On a reduced grid, compute randomized SVD of the mutual-coherence tensor at `z_2={0.1,0.5,1,2,5,10,20} mm`, compare measured covariance-Jacobian rank with `R^2`, and test generic/orthogonal jet slopes at three points.

**Nearest prior art.** Coherent-mode decompositions and the optical memory effect (Freund–Rosenbluh–Feng), but not an `R^2` reconstruction bound paired with a finite-jet theorem.

## 4. Étendue–Information Conductance Law

**Score:** `5 × 4 × 5 = 100`

**Mechanism.** The T5 slowdown is controlled by the number of independently integrated speckle grains. Optical étendue fixes that number and therefore the covariance visibility available per bank.

**Exact statement sketch.** For independent grain weights `w_g`,

```text
C_b^2 = Var(b)/E[b]^2
      = Σ_g w_g^2 /(Σ_g w_g)^2
      = 1/N_eff.                                                (5.1)
```

With `n` photoelectrons and a Gaussian covariance channel, the visibility is approximately

```text
r = n/N_eff,
I_logQ = 1/2 Σ_j [r_j/(1+r_j)]^2.                              (5.2)
```

Thus photons cease to help for `n≫N_eff`, while in the shot-limited regime information scales as `n^2/N_eff^2`. The physical resource is the pair `(étendue, independent banks)`.

**Decisive test.** Sweep illumination aperture and object area at fixed screen statistics; compute `N_eff` directly from the intensity weights and test whether all `T_det` curves collapse against `n/N_eff`.

**Nearest prior art.** Goodman's speckle-contrast/effective-speckle-number law and optical étendue conservation. The new object is its insertion into the exact covariance-information conductance and detection-latency law.

## 5. Mean–Covariance Channel-Ratio Ranging

**Score:** `4 × 4 × 5 = 80`

**Mechanism.** The mean leak is present at contact while covariance information grows as `z_2^4`; their ratio is a steep axial coordinate.

**Exact statement sketch.** In the near-contact shot-limited regime,

```text
J_cov/J_mean
 = [K_1/(2L)] z_2^4 + o(z_2^4).                               (5.3)
```

A calibrated ratio therefore estimates `z_2`, with a Fisher uncertainty obtained from the derivative `4/z_2` of its log.

**Decisive test.** Dense `z_2=0.05…5 mm` sweep for two `z_1` and two screen strengths; fit the quartic region and perform blind depth classification from held-out banks.

**Nearest prior art.** Transport-of-intensity and defocus diversity (Teague, DOI `10.1364/JOSA.73.001434`) plus speckle-contrast ranging. The distinction is a zero-dimensional two-moment range observable with no image.

## 6. Leak-Orthogonalized Efficient Sentinel

**Score:** `4 × 5 × 5 = 100`

**Mechanism.** Treat the measured physical mean leakage as a nuisance tangent, not as an unmodeled false alarm. Project the covariance score off that tangent in the joint mean–covariance score Hilbert space.

**Exact statement sketch.** If `s_μ` is the leak score and `s_Σ` the covariance-change score, the efficient signal is

```text
s_eff = s_Σ - Proj_span{s_μ,s_law} s_Σ,                        (5.4)
```

and retained information is `||s_eff||^2`. The wall need not be exact; specificity requires only that the covariance tangent not lie in the leak/law span.

**Decisive test.** Use the existing T5 bank pool and a calibrated physical mean Jacobian. Compare power and non-target false alarms before/after projection across `z_2=1,5 mm` and one mismatch cell.

**Nearest prior art.** Cox–Reid parameter orthogonality and semiparametric efficient scores. The open optics square is using the physical leak itself as the profiled tangent that restores the wall-based attribution claim.

## 7. Two-Kernel Jet Interferometry

**Score:** `5 × 4 × 4 = 80`

**Mechanism.** A second wavelength, pupil, polarization, or defocus supplies a noncollinear quadratic kernel. Directions first-order-orthogonal in one channel are generically first order in the other, converting an `ε^-4` class to `ε^-2` without reconstruction.

**Exact statement sketch.** For two anchored functionals

```text
Q_a=x^T G_a x,
Q_b=x^T G_b x,
```

joint first-order blindness requires

```text
x^T G_a δ = x^T G_b δ = 0.                                    (5.5)
```

For generic nonproportional `G_a,G_b`, this intersection has codimension two; a designed pair can eliminate it on a target subspace. The joint quotient jet is `m=1` unless both linear forms vanish.

**Decisive test.** Reuse the twin at two wavelengths or two illumination pupils. Construct directions orthogonal under one kernel, measure the joint Chernoff slope, and require conversion from `4` to `2` without an amplitude-gauge failure.

**Nearest prior art.** Multiwavelength speckle, phase diversity, and heterodyne interferometry. The new statement is deliberate engineering of the statistical contact order rather than image recovery.

## 8. Physical Mean-Null / Covariance-Max Witness

**Score:** `4 × 5 × 4 = 80`

**Mechanism.** Replace a nominal Fourier annulus by the exact finite-dimensional physical null of the mean operator while maximizing covariance detectability.

**Exact statement sketch.** With measured mean Jacobian `M_phys` and covariance Fisher `J_cov`, solve

```text
maximize_δ   δ^T J_cov δ
subject to   M_phys δ=0, ||δ||=1.                              (5.6)
```

The solution is the top eigenvector of the projected matrix `Π_ker(M) J_cov Π_ker(M)`. It is an exact discrete witness even when the optical spectrum has sinc tails.

**Decisive test.** Build `M_phys` from T1 and `J_cov` from a T5 geometry, solve (5.6), and compare the resulting mean/covariance `d'` with the current `k=6…9` and guarded `k>10` witnesses.

**Nearest prior art.** Generalized eigenmode discrimination and SPADE-like mode sorting. The novelty is a zero-dimensional physical-wall witness for a fluctuating-medium covariance channel.

---

# Final binding conclusions

1. **No exact `>2k_p` wall exists in the current unfiltered finite-DMD geometry.** There is an ideal compact-support theorem and a noncompact physical pedestal.
2. **Beyond `2k_p` remains the correct next witness**, because it removes the dominant core; add a hard pupil/guarded ROI or construct a physical mean-null witness for exactness.
3. **The channel crossover has a derived fourth-root law.** Near contact, mean leak is order zero in `z_2`; covariance Chernoff information is order four.
4. **Physical jet order is channel selected.** Any unprofiled first-order mean leak makes the strict local order one; projection or an optical wall restores curvature order two.
5. **No mean↔covariance information conservation law exists.** The exact object is an additive Gaussian Fisher decomposition plus a limited weak-phase quadrature identity.
6. **The Letter receives only the four frozen transfer sentences in §3.** All other R44 results remain quarantined.
7. **The next paper is the wave-optical memory-effect Rank–Jet phase diagram**, built around the coherence-separation-rank theorem and the measured support/rank/jet crossover.
