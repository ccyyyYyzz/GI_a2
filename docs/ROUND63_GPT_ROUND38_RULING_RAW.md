# R38 ruling (GitHub issue #29, raw)

Title: R38
Posted: 2026-07-23T19:26:58Z

---

# R38 — Fluctuating-medium measurement diversity: corrected phase boundary, Fisher prognosis, and one-shot ruling

**Reference request:** [`docs/ROUND63_GPT_ROUND38_QUESTION.md`](https://github.com/ccyyyYyzz/GI_a2/blob/9fc8208/docs/ROUND63_GPT_ROUND38_QUESTION.md) at commit [`9fc8208`](https://github.com/ccyyyYyzz/GI_a2/commit/9fc8208f1641f97d41d0d2bd63c5d9b9fc62fbd6).  
**Housekeeping:** R37 is superseded and receives no ruling. The DLGI v1 flagship kill remains binding.

## Executive verdict

> **CONDITIONAL GO to one deliberately hostile falsification probe; not a campaign or manuscript GO. The inversion is genuinely eye-lighting, but the present oxygen test is not a positive blind-recovery result.**

The key distinction is **oracle rank versus blind information**. If the medium realizations `w_t` are known, their spatial fluctuations can indeed rotate a repeated fixed pattern bank out of its original row space and expose fixed-system null modes. But when every `w_t` must be estimated from the same buckets, each epoch spends measurement dimensions on the medium. The exact local phase boundary is therefore much narrower than stacked-operator rank.

For the frozen toy `N=1024, M=64, T=64`:

- `d_w=16` has a full-rank **oracle** stack, but its blind profiled scene rank is at most `17(64-16)=816 < 1024`; full blind scene recovery is locally impossible.
- `d_w=64` and `256` have `d_w >= M`; a free per-epoch medium can generically absorb every bucket perturbation, so the likelihood supplies zero local scene information after profiling unless a temporal prior is doing the identification.
- only `d_w = 26,...,37` can even meet the generic blind local rank count at `M=64`.
- the best point is still almost square: `max_d (d+1)(64-d)=1056`, only `3.1%` above `N`. A random-matrix edge prognosis gives a worst-mode information factor of only about `2.3e-4` of the average. The observed `10^5–10^6` photons/bucket requirement is therefore not surprising; it is the expected near-phase-boundary penalty.

Thus:

> **If `M=64` and arbitrary 1024-pixel dense scenes are binding, kill the flagship now.** The geometry is too close to the blind phase transition for robust, photon-efficient all-mode recovery. The concept survives only by moving the decisive test into a genuinely overdetermined blind regime—roughly `M >= 80`, preferably `M≈96` for a twofold capacity margin—or by explicitly narrowing the claim to a lower-dimensional scene class.

The broad mechanism is also not new: blind structured illumination already owns “unknown random modulation reveals otherwise invisible spatial content.” The potentially new object is the **bucket-only configuration and its exact blind/Fisher boundary**.

---

# 1. Identifiability theorem

## 1.1 Model and load-bearing assumption

Let `P ∈ R^{M×N}` be the known fixed pattern matrix and

```text
w_t = 1 + U z_t,          U ∈ R^{N×d},
b_t = P diag(w_t) x,      t=1,...,T.                            (1.1)
```

The realizations `z_t` are unknown. The spatial subspace `U`, or at minimum its covariance/operator class, must be declared.

This distinction is non-negotiable:

- **unmeasured realizations** is compatible with the claim;
- **unknown arbitrary spatial subspace with no calibration or physical model** is not.

If `U` is arbitrary and unknown, the data contain products `x ⊙ u_k`; a smooth multiplicative envelope can migrate between scene and medium, and fixed-null content is not generally identifiable. “No reference arm” does not mean “no medium model.”

## 1.2 Gauge group

For the affine medium family

```text
W_U = 1 + range(U),
```

the exact multiplicative stabilizer is

```text
G_U = { h ∈ (R\{0})^N : diag(h)^(-1) W_U = W_U }.              (1.2)
```

Every `h ∈ G_U` generates the indistinguishable transformation

```text
x'   = h ⊙ x,
w'_t = w_t ⊘ h.                                                 (1.3)
```

The recoverable scene is the equivalence class `[x]` under this group, together with any source/PDE scale gauge.

Important cases:

1. **Known generic `U`, fixed mean field `1`, known global photon scale.** Generically `G_U={1}`; there is no automatic smooth-envelope gauge.
2. **Unknown global source/PDE scale.** One scalar scale gauge remains.
3. **Unknown static mean transmission `μ(r)`.** Only `μ ⊙ x` is identifiable without an additional normalization; this is the genuine smooth-envelope ambiguity.
4. **Special invariant bases.** Fourier/polynomial/smooth subspaces can have nontrivial or nearly nontrivial diagonal stabilizers. Near-stabilizers appear as tiny Fisher eigenvalues even when exact rank is full.

The probe must compute this stabilizer numerically for the realized `U`; merely setting `E[z_t]=0` in a simulator is not a gauge proof.

## 1.3 Exact local kernel theorem

Define

```text
H_x = P diag(x) U,                                               (1.4)
```

and let `r=rank(H_x)`. Let `Q_x` have rows spanning the left null space of `H_x`, so `Q_x H_x=0`. The Jacobian perturbation at epoch `t` is

```text
δb_t = P diag(w_t) δx + H_x δz_t.                               (1.5)
```

### Theorem 1 — profiled local identifiability

The joint Jacobian kernel consists exactly of perturbations satisfying

```text
G(x,Z) δx = 0,                                                  (1.6)

G(x,Z) = stack_t [ Q_x P diag(w_t) ],                           (1.7)
```

with

```text
δz_t = -H_x^† P diag(w_t)δx + ν_t,
ν_t ∈ ker(H_x).                                                  (1.8)
```

Therefore joint local identifiability up to the gauge tangent holds if and only if

```text
rank(H_x)=d
and
ker G = tangent(G_U · x).                                       (1.9)
```

### Proof

Equation (1.5) is solvable in `δz_t` exactly when its scene term lies in `range(H_x)`. Left multiplication by `Q_x` gives (1.6); the general solution is (1.8). No parameter count is being substituted for the kernel calculation. `square`

This is the production certificate. Stacked oracle rank is not.

## 1.4 Generic phase boundary

Let

```text
v_0 = 1,  v_k = U_:k,
C   = [1, Z] ∈ R^{T×(d+1)}.
```

Because

```text
Q_x P diag(w_t)
 = Q_x P diag(v_0) + sum_k z_tk Q_x P diag(v_k),                (1.10)
```

the profiled operator factors through only `rank(C)` temporal coefficients. Hence

```text
rank G
 <= min{N-g, rank(C)(M-d)}
 <= min{N-g, min(T,d+1)(M-d)},                                  (1.11)
```

where `g` is the gauge dimension.

Under generic `P,U,x,Z`, a trivial stabilizer, and nonvanishing relevant minors, equality holds. Thus the generic **local** phase boundary is

```text
d < M,
min(T,d+1)(M-d) >= N-g.                                         (1.12)
```

For `T >= d+1`, temporal diversity saturates algebraically and the condition becomes

```text
(d+1)(M-d) >= N-g.                                              (1.13)
```

The maximum over `d` is

```text
max_d (d+1)(M-d) = floor((M+1)^2/4).                            (1.14)
```

Therefore no medium dimension can make an arbitrary `N`-dimensional scene locally identifiable unless approximately

```text
M >= 2 sqrt(N-g) - 1.                                           (1.15)
```

For the fixed-system null component alone, conditional on a stable row-space estimate, the corresponding necessary generic capacity is

```text
min(T,d)(M-d) >= N-rank(P).                                     (1.16)
```

Equation (1.12), not `MT >= N+Td`, is the honest boundary. The latter ignores the fact that all epoch operators live in one `(d+1)`-dimensional affine family.

### Oxygen-test adjudication

At `N=1024,M=64,T=64`:

```text
d=16:  (d+1)(M-d)=816      -> blind full-scene impossible
d=31:  (d+1)(M-d)=1056     -> locally possible, extremely ill-conditioned
d=64:  M-d=0               -> zero likelihood scene information after profiling
d=256: rank(H_x)=M<d        -> medium and scene jointly nonidentifiable
```

The warm-start `d=64/256` reconstructions are solutions selected by proximity to truth or by implicit regularization on a nonidentifiable fiber. They are not evidence for blind recovery.

## 1.5 Local is not global

Full local rank does not prove global uniqueness. In lifted form, with

```text
c_t = [1;z_t],
X = x [c_1^T ... c_T^T],
```

the measurements are linear in a structured rank-one object. Global identifiability requires the measurement kernel to avoid the secant/difference variety of admissible lifted objects, modulo the gauge. The probe must therefore include an adversarial collision search and multi-start agreement; one converged warm solution is insufficient.

## 1.6 What OU correlation does—and does not do

A nondegenerate OU process has full-support innovations. Therefore:

> **An OU prior does not repair deterministic likelihood nonidentifiability.**

It can:

- provide a tracking initialization after the first epoch;
- improve finite-noise conditioning;
- shrink a Bayesian posterior;
- reduce the effective number of rapidly varying nuisance degrees of freedom when the innovations are genuinely small.

It cannot turn `d>=M` or a rank-deficient `G` into data identifiability. Any information appearing only after adding the OU precision is prior-supported and must be reported separately.

For conditioning and covariance estimation, use an effective independent-epoch count such as

```text
T_eff ≈ T (1-φ)/(1+φ),      φ=exp(-Δ/t_c),                      (1.17)
```

as a first-order diagnostic. Algebraic rank may remain full as `φ→1`, while useful diversity collapses numerically.

There is a second, distinct route: integrate out `z`. With known spatial covariance `K_w` and temporal correlation `R_ts`,

```text
Cov(b_t,b_s)
 = R_ts P diag(x) K_w diag(x) P^T + shot noise.                 (1.18)
```

This marginal covariance can identify `x` statistically even when pathwise joint recovery is not injective. It is also the natural spectral initializer. But this is a **known-medium-law** result, not recovery of arbitrary unknown realizations, and it is precisely the mechanism closest to blind-SIM theory.

---

# 2. Fisher and shot-noise accounting

## 2.1 Exact profiled Fisher matrix

Let

```text
Y_it ~ Poisson(μ_it),
μ_t = Φ_t P diag(w_t)x,
A_t = Φ_t P diag(w_t),
H_t = Φ_t P diag(x)U,
W_t = diag(1/μ_t).                                               (2.1)
```

The likelihood blocks are

```text
I_xx   = sum_t A_t^T W_t A_t,
I_xz_t = A_t^T W_t H_t,
I_z_tz_t = H_t^T W_t H_t.                                      (2.2)
```

Let `Q_z` be the OU prior precision, if used. The efficient scene information is

```text
J_x
 = I_xx - I_xz (I_zz + Q_z)^† I_zx.                            (2.3)
```

The likelihood-only result sets `Q_z=0`. The difference between the two is the prior contribution and must not be hidden.

For a basis `V_0` of `ker P`, the information that actually breaks the fixed wall is

```text
J_0 = V_0^T J_x V_0.                                            (2.4)
```

All campaign claims should be made from the eigenvalues of `J_0`, not from the rank of the oracle stacked matrix.

If `d>=M`, `H_t` is generically row-full. With free per-epoch `z_t` and `Q_z=0`, the Schur projection removes the entire measurement space at that epoch. This is the exact Fisher version of the phase-boundary failure above.

## 2.2 Small-fluctuation scaling

Since `P V_0=0`, at small spatial fluctuation

```text
A_t V_0
 ≈ Φ_t P diag(Uz_t)V_0.                                        (2.5)
```

Define the physically meaningful pixelwise medium RMS

```text
σ_f^2 = (1/N) tr(U Σ_z U^T).                                   (2.6)
```

For isotropic, well-normalized patterns, the **average oracle** null-mode information scales as

```text
λ_avg(J_0^oracle)
 ~ (PH/N) σ_f^2 × pattern/scene constant,                       (2.7)
```

where `PH` is the total detected photon budget across all buckets and epochs. The blind eigenvalues are (2.7) multiplied by profiling efficiencies in `[0,1]`.

Two consequences:

1. At fixed total `PH`, increasing `T` improves rank, isotropy and concentration only until the temporal medium span is filled. It does not create photons; after that point the information trace does not grow merely because the same dose was split into more epochs.
2. The useful signal is second order in fluctuation amplitude. Halving `σ_f` costs about four times the photons before profiling.

### Required normalization correction

The oxygen test states orthonormal `U` and per-coefficient `sd(z_k)=0.30`. Then

```text
σ_f = 0.30 sqrt(d/N),                                           (2.8)
```

so the actual pixelwise RMS values were approximately

```text
d=16:  0.0375
d=64:  0.075
d=256: 0.150.                                                    (2.9)
```

The apparent improvement with `d` partly increased the physical fluctuation energy. Every production sweep over `d` must hold `σ_f`, positivity and total transmission fixed. Otherwise geometry and signal amplitude are confounded.

Use a positive model—preferably a mean-normalized lognormal field—or a rigorously bounded linearized field. Gaussian `1+Uz_t` cannot silently generate negative transmission.

## 2.3 Recoverable-mode count

Let `λ_j` be the eigenvalues of the likelihood-only `J_0`, and let `a_j` be the scene coefficient in the corresponding mode. Define

```text
k_rec(γ)
 = # {j : a_j^2 λ_j >= γ^2}.                                   (2.10)
```

This is the number of null modes recoverable at Fisher SNR at least `γ`. A rank count corresponds only to `γ=0` and is not an imaging result.

Useful rank ceilings are

```text
rank(J_0^oracle)
 <= min{N-M, M min(d,T-1)},                                     (2.11)

rank(J_0^blind)
 <= min{N-M, min(T,d)(M-d)}                                    (2.12)
```

for generic full-row-rank `P`, after separating the fixed row-space component.

In a random-like regime let

```text
χ_0 = [min(T,d)(M-d)]/(N-M).                                   (2.13)
```

A Marchenko–Pastur edge prognosis gives approximately

```text
λ_min/λ_avg ≈ (1-χ_0^(-1/2))^2,      χ_0>1.                    (2.14)
```

This is a prognosis, not a theorem for the structured operator, but it explains the oxygen result. At `M=64,d=31,N=1024`,

```text
χ_full = 1056/1024 = 1.031,  edge ≈ 2.3e-4,
χ_0    = 1023/960  = 1.066,  edge ≈ 9.8e-4.                    (2.15)
```

The worst revealed modes need roughly `10^3–10^4` times the photons of an average mode even before model mismatch. A robust arbitrary-scene probe should not sit there.

For comparison, `M=96,d≈47` gives a full-scene capacity around `2352`, over twice `N`, and a much healthier random-edge prognosis. More generally:

```text
capacity ratio >=1.5 requires about M>=79 for N=1024,
capacity ratio >=2   requires about M>=91.                       (2.16)
```

That is why `M≈96` is the minimum honest primary setting for the decisive probe.

## 2.4 Profiling tax

Let `J_0^or` be the known-medium oracle information. Define the generalized eigenvalues

```text
J_0 v_j = η_j J_0^or v_j,       0<=η_j<=1.                     (2.17)
```

The `η_j` values are the exact multiplicative-information tax for estimating the medium. Report their full spectrum. Trace ratios alone can hide dead edge modes.

An OU prior may raise the apparent `η_j`; therefore every result must report:

- likelihood-only `η_j`;
- posterior/Bayesian `η_j`;
- the fraction of final precision attributable to the prior;
- sensitivity to `t_c`, covariance and basis misspecification.

If the claimed null recovery disappears when the temporal prior is weakened or misspecified, the method is not extracting a free physical operator; it is imputing one.

## 2.5 Crossover against averaging—and the missing comparator

Classic medium averaging estimates the mean operator and can improve the fixed row-space component. It can never recover a true `ker P` component when `E[w_t]=1`; its null information is exactly zero. Joint inference can therefore win by a made-possible-at-all effect in that subspace.

But it can lose elsewhere because photons and degrees of freedom are spent estimating `w_t`. The correct equal-budget endpoint is

```text
null-space MSE reduction
minus
row-space variance/profile penalty.                            (2.18)
```

There is also a stronger practical comparator:

> If a programmable DMD can display `MT` distinct known patterns at the same exposure and switching cost, fresh known patterns generate the same rank without blind estimation and should dominate.

The “second DMD for free” claim is defensible only in a declared regime where the unique code library or pattern-update rate is constrained, or where many medium realizations occur within one pattern-hold interval. If fresh known patterns are available at equal cost, this is not a universal advantage; it is a more difficult way to obtain diversity already available in software.

---

# 3. Optimal exploitation architecture

## 3.1 Recommended pipeline

A cold-start method should separate the easy fixed row-space problem from the new null-space problem.

### Stage A — mean/row-space estimate

Use the temporal mean and the known `P` to obtain a nonnegative row-space scene estimate. Do not ask a bilinear optimizer to discover this component and the fluctuation simultaneously.

### Stage B — centered spectral/covariance initialization

Center the epoch matrix:

```text
B_c = B (I-11^T/T)
     = P diag(x) U Z_c^T.                                      (3.1)
```

Use its singular subspace and the known/predeclared medium covariance to initialize `x` and the medium subspace coefficients. Equivalently, fit the lagged covariance (1.18), with the Poisson diagonal removed analytically.

This is the correct cold-start bridge and the direct point of contact with blind-SIM. It must be a named baseline, not quietly absorbed into the proposed method.

### Stage C — Poisson tracking EM / Gauss–Newton

With the spectral initializer:

- E step: OU Kalman/extended-Kalman/particle smoothing for `z_{1:T}` conditional on `x`;
- M step: nonnegative Poisson update for `x` conditional on the path;
- continuation: increase `d` or `σ_f` gradually;
- fixed trust region and multiple truth-independent starts;
- terminate by held-out likelihood and gauge-invariant residuals, not training loss alone.

A low-rank lifted or Burer–Monteiro formulation is a valuable initializer/diagnostic. A full convex lift is likely too large at `N=1024,T=64`; it should not be the production algorithm unless its memory and photon accounting are demonstrated.

## 3.2 Can pattern reordering convexify the cold start?

Not in the full Cartesian experiment. If every pattern `P_i` is measured under every epoch `t`, a permutation merely relabels rows and leaves the likelihood and bilinear geometry unchanged.

If only a subset of patterns is assigned to each medium state, a balanced randomized/Latin schedule can improve lifted incoherence and conditioning. Complement or anchor patterns can also stabilize scale. But:

- this does not change a failing dimension boundary;
- repeated anchors consume scene rows and are pilot-like resources;
- known scalar source modulation can be divided out and adds no spatial diversity;
- any schedule benefit must be compared with spending the same slots on fresh known patterns.

Call this **conditioning**, not convexification.

## 3.3 Cold start is a primary endpoint

Truth-plus-15% initialization is diagnostic only. No publication or campaign bar may read it. The primary estimator must start from data-only spectral/covariance initialization, with random starts and no access to true `x`, `z`, `w` or their realized subspace.

---

# 4. Prior-art fence

## 4.1 The mechanism-level killing paper

Mudry et al., **“Structured illumination microscopy using unknown speckle patterns,”** *Nature Photonics* 6, 312–315 (2012), DOI [10.1038/nphoton.2012.83](https://doi.org/10.1038/nphoton.2012.83), is the killing paper for the broad slogan.

It already shows that several uncontrolled, unknown random speckle illuminations can be jointly exploited to recover object frequencies invisible to the conventional passband, without illumination calibration. Therefore the following are forbidden:

- first use of unknown spatial modulation as imaging diversity;
- first recovery of otherwise invisible content from unknown speckles;
- first inversion of random illumination from nuisance into resource;
- first blind joint object/illumination reconstruction.

The theoretical continuation is also occupied:

- Labouesse et al., **“Joint Reconstruction Strategy for Structured Illumination Microscopy With Unknown Illuminations,”** DOI [10.1109/TIP.2017.2675200](https://doi.org/10.1109/TIP.2017.2675200).
- Idier et al., **“On the Super-Resolution Capacity of Imagers Using Unknown Speckle Illuminations,”** DOI [10.1109/TCI.2017.2771729](https://doi.org/10.1109/TCI.2017.2771729). Its covariance mechanism is particularly close to (1.18).

## 4.2 Bilinear/self-calibration ownership

The gauge, lifting, spectral initialization and sample-complexity machinery are classical:

- Ling & Strohmer, **“Self-calibration and biconvex compressive sensing,”** DOI [10.1088/0266-5611/31/11/115002](https://doi.org/10.1088/0266-5611/31/11/115002).
- Ling & Strohmer, **“Self-Calibration and Bilinear Inverse Problems via Linear Least Squares,”** DOI [10.1137/16M1103634](https://doi.org/10.1137/16M1103634).
- Ahmed, Recht & Romberg, **“Blind Deconvolution Using Convex Programming,”** DOI [10.1109/TIT.2013.2294644](https://doi.org/10.1109/TIT.2013.2294644).
- Kech & Krahmer, **“Optimal Injectivity Conditions for Bilinear Inverse Problems…,”** DOI [10.1137/16M1067469](https://doi.org/10.1137/16M1067469).

The new theorem cannot be advertised as generic bilinear identifiability. Its residual value is the structured bucket-GI specialization and the oracle-versus-profiled phase split.

## 4.3 Scattering as encoder and ghost-imaging neighbors

- Liutkus et al., **“Imaging With Nature: Compressive Imaging Using a Multiply Scattering Medium,”** DOI [10.1038/srep05552](https://doi.org/10.1038/srep05552), uses a scattering medium as a random encoder, but characterizes/calibrates a stable transmission matrix.
- Paniagua-Diaz et al., **“Blind Ghost Imaging,”** DOI [10.1364/OPTICA.6.000460](https://doi.org/10.1364/OPTICA.6.000460), removes direct knowledge of transmitted patterns by using correlated reflected patterns—there is still side/reference information.
- Meyers, Deacon & Shih, **“Turbulence-free ghost imaging,”** DOI [10.1063/1.3567931](https://doi.org/10.1063/1.3567931), concerns a reference-arm correlation architecture and turbulence immunity, not bucket-only blind null-space recovery. The rhetoric should not imply that all GI deliberately discards a universally available channel.
- Peng & Chen, **“Ghost Imaging Through Complex Scattering Media with Random Light Disturbance,”** DOI [10.1063/5.0252090](https://doi.org/10.1063/5.0252090), and Peng, Zhang & Chen, **“Super-Resolution Ghost Imaging Through Complex Scattering in Dynamic Media,”** DOI [10.1002/lpor.202502103](https://doi.org/10.1002/lpor.202502103), are dangerous keyword neighbors. Their reported architecture suppresses/corrects dynamic scattering using learned/deep-image priors; I found no explicit fixed-pattern null-space-diversity theorem or reference-free recovery of the spatial medium realization.

## 4.4 What appears open

I found no exact paper occupying the complete square:

> **a bucket-only detector with a reused known code bank, an uncontrolled and unmeasured dynamic spatial medium, no reference/wavefront/realization measurement, explicit recovery of directions in the fixed-code null space, and a blind-identifiability plus profiled-Fisher theorem.**

That is a defensible novelty candidate, not a certified novelty claim.

Freeze the claim wording as:

> **Under a declared low-dimensional spatial-statistical model of an unmeasured dynamic medium, repeated bucket measurements with a fixed programmed code bank can acquire scene directions lying in the fixed-bank null space. We characterize the blind identifiability and profiled-Fisher boundary and test reference-free cold-start recovery.**

Do not say “no calibration” if `U` or `K_w` was learned from separate data. Say **no realization measurement**. Do not say turbulence is always beneficial or free; it is beneficial only inside the Goldilocks dimension, fluctuation, cadence and photon regime.

---

# 5. Blunt judgment

## Eye-lighter score

The inversion itself is strong:

> a scalar medium is the exact zero-diversity degeneracy, while spatial fluctuation can turn a provable fixed null space into measured directions.

That is cleaner and more surprising than another reconstruction improvement. A fixed system moving from “mathematically impossible” to “identifiable” is the right kind of effect for the operator’s bar.

## Why it is not yet the sledgehammer

1. The broad mechanism is owned by blind-SIM.
2. The current oracle full-rank evidence is not blind evidence.
3. The no-reference claim still requires a known or accurately modeled medium subspace/covariance.
4. At `M=64,N=1024`, the only blind-identifiable dimensions lie in an almost-square, noise-hostile window.
5. The first toy already needs `10^5–10^6` photons/bucket before the blind profiling tax.
6. A DMD permitted to display fresh patterns is a simpler diversity generator and must be an equal-resource comparator.
7. Cold-start failure is currently complete.

### Binding judgment

> **The concept earns one life-or-death probe, not emotional investment. `M=64` dense-scene recovery is killed analytically. Continue only in a capacity-margin regime and kill the direction if the profiled Fisher spectrum or cold-start result misses the bars below.**

If the probe passes, this is a strong specialized paper candidate. It is not universal until it survives unknown-basis mismatch and beats the practical fresh-pattern comparator in a real fixed-code/fast-medium operating regime.

---

# 6. Frozen one-shot feasibility probe

## 6.1 Sealed banks

Create four disjoint banks before endpoint generation:

1. **Calibration/theory bank:** algorithm selection, threshold calibration, no final reporting.
2. **Confirmatory bank:** fresh scenes, pattern matrices, medium bases, OU paths and Poisson seeds; sealed until code and bars freeze.
3. **Mismatch bank:** basis rotation, covariance error, correlation-time error, lognormal/nonlinear medium and weak spatial-transfer residuals.
4. **Oracle bank:** the same confirmatory records with `w_t` revealed; diagnostic ceiling only.

No scene, pattern, medium path, seed or initialization may cross banks.

## 6.2 Primary geometry

Use `N=1024` and physical nonnegative or complementary patterns with full photon/time accounting.

### Negative controls

```text
M=64, d in {16,31,64}
scalar gain / d=1
```

These must reproduce the theorem: oracle diversity at `d=16`, blind rank failure at `16`, near-edge conditioning at `31`, and zero free-path profiled information at `64`.

### Primary blind-capacity cells

```text
M in {96,128}
d chosen for capacity ratios χ_full in approximately {1.5,2,3}
T in {64,128}, with T>=d+1 where applicable.                        (6.1)
```

Suggested `M=96` cells are `d∈{32,48,64}`. Do not select cells after seeing reconstruction results.

Scenes must include:

- nonnegative dense synthetic scenes with controlled fixed-null energy;
- fresh natural/bridge-style scenes;
- null-energy fractions spanning roughly 25%, 50% and 75%;
- no learned scene prior in the primary arm.

A low-dimensional scene-prior arm may be secondary, but cannot rescue the arbitrary-scene verdict.

## 6.3 Medium and photon grid

Hold the actual field RMS, mean transmission and positivity fixed across `d`:

```text
σ_f in {0.05,0.15,0.30};
T_eff/d in approximately {0.5,1,2,4};
mean detected photons/bucket in {10^3,10^4,10^5,10^6}.             (6.2)
```

Use a predeclared smooth basis family and a positive mean-normalized medium model. The primary headline budget is at most `10^5` detected photons/bucket. `10^6` is a stress/bench-ceiling cell, not a universal result.

## 6.4 Arms

1. **Fixed-bank medium averaging:** strongest classic reconstruction from the same repeated patterns and total resources.
2. **Fresh known-pattern comparator:** `MT` distinct patterns when the hardware/schedule permits it.
3. **Known-medium oracle:** exact `w_t`, same photons.
4. **Blind covariance/spectral initializer plus Poisson tracking EM:** proposed production arm.
5. **Blind-SIM-style covariance-only reconstruction:** mechanism-level prior-art baseline.
6. **Cold random ALS:** optimization baseline.
7. **Warm truth-neighborhood ALS:** diagnostic only; barred from every endpoint.

All arms receive identical detected photons, exposures, elapsed time and physical mask implementation. If the fresh-pattern arm is declared infeasible, the exact hardware/code-rate constraint must be frozen and demonstrated rather than asserted.

## 6.5 Binding bars

### Bar F0 — theorem and gauge

For every primary cell, the numerical Jacobian kernel must match Theorem 1, including the predicted negative controls. The affine stabilizer/gauge must be explicitly computed. The claimed production cells require

```text
χ_full >= 1.5
and
smallest non-gauge profiled singular value^2 / median >= 0.02.    (F0)
```

**Kill** on a rank mismatch, unremoved envelope gauge, or reliance on an unreported mean-field calibration.

### Bar F1 — oracle Fisher viability before blind reconstruction

At the primary `<=10^5` photons/bucket cell:

- at least 50% of the frozen scene’s null-space energy must lie in modes with oracle Fisher SNR `>=3`;
- the oracle null NMSE must be `<=0.15` on at least 90% of sealed replicates;
- the scalar-gain control must remain at null NMSE approximately one.

**Kill immediately** if the oracle fails. Do not run or tune a blind solver to rescue a photon-dead channel.

### Bar F2 — profiling tax and data dominance

On the modes claimed recoverable:

```text
median η_j >=0.25,
10th-percentile η_j >=0.10,                                     (F2a)
```

and at least half of the final precision must remain in the likelihood-only information when the OU prior is removed from the Fisher calculation.

Perturbing `t_c`, covariance amplitude and basis by the preregistered mismatch must change null NMSE by no more than 25% in the claimed domain.

**Kill** if the result is predominantly prior imputation or exact-basis calibration.

### Bar F3 — blind cold-start null recovery

Using no truth-dependent initialization:

- median null NMSE `<=0.25` (`>=6 dB` null recovery);
- at least 80% of sealed replicates have null NMSE `<=0.5`;
- the blind method captures at least 50% of the oracle improvement over fixed averaging;
- independent data-only starts agree modulo the declared gauge, and an adversarial collision search finds no equal-likelihood remote solution in the claimed cell.

**Kill** if only warm-start or selected-start runs succeed.

### Bar F4 — image-level and equal-resource effect

Against fixed-bank averaging at identical resources:

```text
median full-image gain >=3 dB,
90% scene-bootstrap lower bound >1 dB,
at least 80% of scenes improve.                                 (F4)
```

The gain must be concentrated neither in one artificial null witness nor in a learned scene prior.

Against the fresh-known-pattern arm, either:

- the proposed method is competitive under the same time/photon ledger; or
- the manuscript freezes a measured code-library/update constraint under which the fresh-pattern arm is physically unavailable and quantifies the resulting systems advantage.

If fresh patterns are feasible and dominate, **kill the universal/software-only flagship claim**; retain at most a fixed-code niche result.

### Bar F5 — physical model transfer

Repeat the primary cell with:

- physical nonnegative/complement masks;
- mean-normalized lognormal or otherwise positive transmission;
- `10–20%` spatial-basis rotation/mismatch;
- `2x` correlation-time mismatch;
- weak pattern-dependent residual transfer beyond pure multiplication.

Held-out bucket residuals must remain calibrated and performance degradation must stay below 25% inside the claimed domain.

**Kill** if signed-Gaussian patterns, negative transmission or exact simulator-matched `U` are load-bearing.

### Bar F6 — edge honesty and dose statement

Publish the phase map over `(M,d,T_eff,σ_f,PH)` including:

- scalar/no-diversity edge;
- `d>=M` nuisance-absorption edge;
- near-square Fisher edge;
- slow-medium `T_eff` edge;
- weak-fluctuation photon edge;
- fresh-pattern dominance region.

No `10^6`-photon result may be described as photon-efficient or universal.

## 6.6 Kill tree

1. **F0 fails:** theorem/model killed; archive.
2. **F0 passes, F1 fails:** geometric rank is noise-fiction; kill before blind inversion.
3. **F1 passes, F2 fails:** oracle-only/prior-driven phenomenon; materials-bank theorem only.
4. **F2 passes, F3 fails:** bilinear cold start defeats software-only deployment; kill method.
5. **F3 passes, F4 fails:** mathematically interesting but no image-level value; short theory note at most.
6. **F4 passes, F5 fails:** exact-model simulation artifact; no flagship.
7. **All pass but fresh patterns dominate whenever available:** narrow fixed-code application, not universal sledgehammer.
8. **All pass in a measured fixed-code/fast-medium regime:** paper GO and bench-transfer GO.

No repair round follows this probe. Thresholds, code and banks freeze before confirmatory generation.

---

# Final ruling

> **GO once, adversarially. The medium can be a second pattern generator in the oracle experiment, but the blind experiment has a Goldilocks dimension and a severe profiling/Fisher tax. The current `M=64` oxygen test sits outside or directly on that boundary and does not justify a campaign. Move to a capacity-margin setting, normalize physical fluctuation energy, require cold-start recovery at `<=10^5` photons/bucket, and include the fresh-pattern comparator. If any of those conditions fails, kill the direction immediately.**