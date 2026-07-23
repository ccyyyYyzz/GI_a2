# R39 ruling (GitHub issue #30, raw)

Title: R39
Posted: 2026-07-23T23:46:03Z

---

# R39 — Beyond-modulator-band covariance imaging: aperture law, information rate, and sealed-probe ruling

**Reference request:** [`docs/ROUND63_GPT_ROUND39_QUESTION.md`](https://github.com/ccyyyYyzz/GI_a2/blob/d6a2801/docs/ROUND63_GPT_ROUND39_QUESTION.md) at commit [`d6a2801`](https://github.com/ccyyyYyzz/GI_a2/commit/d6a280111481ab3ee7182d76fbee05d59863dc11).  
**Predecessor:** R38 / issue #29. The pathwise blind branch is dead; only the marginal covariance branch is in scope.

## Executive verdict

> **CONDITIONAL GO to exactly one publication-scale, sealed, simulation-only campaign. This is now the strongest surviving flagship candidate, but it is not promoted on the 16×16 result.**

The inversion is genuinely strong: a band-limited modulator is provably blind in its first-moment channel to frequencies outside its own band, while an unmeasured dynamic thin screen can shift those frequencies into a bucket-covariance channel whose support is the Minkowski sum of the pattern and medium bands. The mean-channel impossibility and covariance-channel possibility form a clean made-possible-at-all separation.

The central E8 question has a precise answer:

> **There is no universal nonzero blind-error floor such as 0.4 under an injective, correctly specified covariance law. The fundamental penalty is a slower information rate, not an asymptotic floor.**

If the profiled covariance Fisher matrix is positive definite on the beyond-band scene subspace, the blind marginal MLE has error proportional to `T_eff^(-1/2)` and is consistent. A true floor occurs only if a beyond-band scene derivative lies in the nuisance tangent of the medium law or in an exact gauge. The observed `0.44–0.48` stall is therefore one of two things:

1. a finite-`T` CRB caused by very small Fisher eigenvalues near the aperture edge; or
2. an estimator/initialization plateau.

It is not licensed as a physical floor until the exact profiled CRB is computed. The existing nonmonotone MLE curve, the continued improvement of moment matching, and the exact-covariance injectivity result make estimator immaturity a serious suspect. The sealed campaign must use a lifted covariance initializer followed by Fisher-scored marginal-likelihood refinement; pathwise ALS is retired.

The blind/oracle gap itself is real and fundamental. Once shot noise is below the fluctuation covariance, extra photons continue to improve a realization-aware oracle but no longer buy comparable blind covariance information. **Independent medium realizations, not photon flux, become the primary currency.** This explains why the oracle can reach `0.091` while the blind route remains much worse at the same `T`, without implying a nonzero blind floor.

**Immediate kill before launch:** if the 64×64 profiled-Fisher preflight predicts beyond-band NRMSE above the frozen limits in §6 at `T_eff ≤ 4096` and `10^4` detected photoelectrons per physical bucket, do not run the reconstruction campaign. Archive the aperture and wall theorems as a materials-bank note.

---

# 1. Model and exact covariance channel

Let `x ∈ R_+^N` be a static object on a finite periodic grid. A fixed bank of effective signed codes has rows `p_i^T`, collected in `P ∈ R^{M×N}`. Physical signed rows are implemented by complementary nonnegative exposures and charged twice in the ledger.

For medium epoch `t`,

```text
w_t = μ + δw_t,
E[δw_t] = 0,
Cov(δw_t,δw_s) = r_ts K_w,
```

and the bucket vector is

```text
Y_t = Φ P diag(w_t) x + ε_t.                                   (1.1)
```

In the declared Gaussian/linearized observation model, integrating out the unmeasured realization gives

```text
m(x)   = Φ P diag(μ) x,                                        (1.2)
C(x)   = Φ^2 P diag(x) K_w diag(x) P^T,                         (1.3)
Cov(Y_t,Y_s) = r_ts C(x) + 1{t=s} R_shot(x).                    (1.4)
```

Equation (1.3), not an estimated `w_t`, is the production object. It can also be written as a linear map of a lifted rank-one scene matrix:

```text
C(x) = L(xx^T),
L(X) = Φ^2 P [K_w ⊙ X] P^T.                                    (1.5)
```

This representation is the route around the pathwise nonidentifiability found in E4/E6.

---

# 2. Aperture-law theorem

Assume a stationary medium covariance on the discrete torus. Let `S_w(k)` be its spatial power spectrum, with support `B_w`, and let every code have Fourier support inside `B_p`.

For

```text
q_i(k;x) = Σ_ξ p̂_i(k-ξ) x̂(ξ),
```

Parseval gives

```text
C_ij(x)
 = Φ^2 Σ_{k∈B_w} S_w(k) q_i(k;x) overline{q_j(k;x)}.            (2.1)
```

Define the covariance aperture

```text
A_cov = B_w - B_p
      = { ξ : ξ = k-η for some k∈B_w, η∈B_p }.                 (2.2)
```

For symmetric bands this is the usual Minkowski sum `B_w ⊕ B_p`.

## Theorem 1 — exact covariance aperture

Under (1.1)–(1.4):

1. **Exact exclusion.** If a perturbation `h` has Fourier support disjoint from `A_cov`, then

   ```text
   C(x+h) = C(x)                                                 (2.3)
   ```

   for every scene `x`. Such a perturbation carries zero covariance Fisher information.

2. **Potential inclusion.** A frequency `ξ` can affect the covariance channel only if the coverage vector

   ```text
   r_ξ(i,k) = sqrt(S_w(k)) p̂_i(k-ξ),   k∈B_w                  (2.4)
   ```

   is nonzero. Thus bare spectral support is determined by (2.2).

3. **Local recoverability.** Let `D C_x` be the derivative of (1.3). On a declared scene model `X`, local recovery inside `A_cov` holds exactly when

   ```text
   ker(D C_x | T_x X) = T_x(G·x),                               (2.5)
   ```

   where `G` is the true gauge group. Bare support is not sufficient; the covariance Jacobian must have full rank after gauge projection.

4. **Global recoverability.** Global uniqueness up to gauge holds exactly when

   ```text
   ker L ∩ {xx^T - x'x'^T : x,x'∈X} = {0 modulo G}.             (2.6)
   ```

   This is a rank-one secant condition. It is the covariance analogue of phase-retrieval injectivity and must be checked numerically on the realized code bank.

### Proof of exclusion

If `ĥ` is supported outside `A_cov`, then for every `k∈B_w` and every code `i`,

```text
Σ_ξ p̂_i(k-ξ) ĥ(ξ) = 0.
```

Hence `q_i(k;x+h)=q_i(k;x)` on every medium frequency carrying nonzero power. Substitution into (2.1) proves (2.3). The local and global conditions are the standard tangent and secant conditions for the lifted map (1.5). `square`

## 2.1 Genericity and code requirements

For continuously randomized, incoherent Fourier coefficients inside `B_p`, Jacobian rank is algebraically generic once one full-rank witness exists. Necessary counting is

```text
M(M+1)/2 ≥ d_A - g,                                             (2.7)
```

where `d_A` is the number of real scene degrees of freedom claimed inside `A_cov` and `g` is the gauge dimension. This count is not sufficient for structured Fourier lattices; the committed collapse of lattice codes is therefore expected. The campaign must use band-limited incoherent random codes and certify the actual Jacobian/secant geometry.

## 2.2 Gauge

There is no automatic smooth-envelope gauge when `K_w`, the mean field `μ=1`, and the global photon scale are fixed. The genuine ambiguity appears when the medium has an unknown static envelope:

```text
w_t = g ⊙ (1+u_t),
K_w = diag(g) K_0 diag(g).                                      (2.8)
```

Then both mean and covariance depend on

```text
q = g ⊙ x,                                                       (2.9)
```

and only `q` is identifiable without an external normalization or a low-dimensional model for `g`. A covariance-only model also has the usual `x↔-x` ambiguity; nonnegativity and the in-band mean channel fix it. Unknown covariance amplitude creates a scalar `x`/variance gauge unless the mean channel or `E[w]=1` fixes scene scale.

The primary campaign fixes `E[w]=1` and profiles only a finite-dimensional law nuisance. Static-envelope uncertainty is a stress test, not silently ignored.

## 2.3 Fisher-weighted aperture, not support art

Define the overlap multiplicity

```text
m(ξ) = Σ_{k∈B_w} S_w(k) Σ_i |p̂_i(k-ξ)|^2.                    (2.10)
```

It is positive on the interior of `A_cov` and tends to zero near its edge for finite bands. The theorem therefore gives a **support aperture**; finite-data resolution is smaller.

Let `β` parameterize scene Fourier coefficients, `V` be the one-epoch bucket covariance including shot noise, and

```text
D_β = ∂ vech(V) / ∂β.
```

For `T_eff` independent-equivalent epochs, the covariance Fisher matrix is

```text
I_ββ = (T_eff/2) D_β^T
       (V^-1 ⊗_s V^-1) D_β.                                    (2.11)
```

After profiling all other scene coefficients and medium-law parameters, define the conditional information of frequency `ξ` by

```text
j_cond(ξ) = 1 / [(J_β^†)_ξξ].                                  (2.12)
```

The usable aperture at Fisher SNR `γ` is

```text
A_γ = { ξ : |x̂(ξ)|^2 j_cond(ξ) ≥ γ^2 }.                       (2.13)
```

The production claim is therefore:

> the support law reaches `k_p+k_w`, while the preregistered finite-data claim is a robust `1.8×` linear cutoff when `k_w=k_p`.

The outer `1.8–2.0×` rim is an edge audit, not part of the primary quality claim.

---

# 3. First-moment impossibility and information ordering

## Theorem 2 — band wall of the mean channel

The first-moment operator is

```text
m(x)=Φ P diag(μ)x.                                               (3.1)
```

If `μ=1` and every admissible code is supported in `B_p`, then for every perturbation `h` supported outside `B_p`,

```text
m(x+h)=m(x),
I_mean h = 0.                                                    (3.2)
```

This remains true for:

- an arbitrarily large number of fresh or reused codes;
- any temporal schedule or allocation;
- any estimator whose statistical model uses only the mean response.

Fresh code directions cannot leave a hardware band they do not possess. If the static mean field has band `B_μ`, replace the wall by `B_p-B_μ`.

This is the formal content of the observed `1.000` baseline. It is an information statement, not a failure of a particular reconstruction algorithm.

## 3.1 Which channel carries beyond-band information?

For the declared zero-mean Gaussian medium, the marginal bucket vector is Gaussian before the Poisson/readout approximation. The mean and covariance are then sufficient: the beyond-band information first appears in (1.3), and medium cumulants of order `k≥3` vanish exactly.

For weak non-Gaussian fluctuations, a `k`th cumulant can carry additional polyspectral combinations, but its signal scales as `σ_f^k`; a sample-cumulant information proxy scales no better than

```text
T_eff σ_f^(2k)                                                   (3.3)
```

before conditioning and combinatorial variance. Relative to second order, each extra order pays approximately `σ_f^2` in the weak-fluctuation regime. The committed third-order failure is therefore the expected boundary, not an implementation surprise.

The manuscript must not claim that covariance is universally sufficient for every medium law. It is sufficient under the declared Gaussian law and information-dominant under the tested weakly non-Gaussian laws.

---

# 4. Sample complexity and the blind/oracle gap

## 4.1 Exact profiled Fisher matrix

Let `η` collect all nuisances: in-band scene coefficients, covariance amplitude, temporal correlation time, spatial-band edge/slope, anisotropy, and any admitted envelope coefficients. For the full temporal covariance `Σ(x,η)`, the Gaussian Fisher blocks are

```text
I_uv = (1/2) tr[Σ^-1 Σ_u Σ^-1 Σ_v]
       + m_u^T Σ^-1 m_v.                                       (4.1)
```

For beyond-band coefficients `β`, the mean derivative is zero. The efficient information is

```text
J_B = I_ββ - I_βη I_ηη^† I_ηβ.                                (4.2)
```

Equivalently, after temporal whitening,

```text
J_B
 = (T_eff/2) D_B^T W^(1/2)
   [I - Π_{W^(1/2)D_η}]
   W^(1/2)D_B,                                                  (4.3)

W = V^-1 ⊗_s V^-1.                                              (4.4)
```

The projector makes the profiling tax explicit. If a scene derivative lies in the medium-law tangent, its eigenvalue is zero and that mode has a genuine floor/gauge. Otherwise all positive eigenvalues scale linearly with `T_eff`.

## Theorem 3 — no-floor criterion

Assume the marginal model is correctly specified, regular, globally identifiable up to the frozen gauge, and

```text
J_B/T_eff ≻ 0                                                   (4.5)
```

on the claimed beyond-band subspace. Then the marginal MLE is consistent and asymptotically normal,

```text
sqrt(T_eff)(β̂-β) -> N(0, J̄_B^-1),
J̄_B = J_B/T_eff.                                               (4.6)
```

Therefore the asymptotic beyond-band error floor is zero. A nonzero floor occurs only if (4.5) fails, the model is misspecified, or the optimizer does not reach the likelihood basin.

This resolves the logical status of E8: `0.44–0.48` is not a fundamental constant. The sealed Fisher calculation must determine whether it is the finite-`T=2048` CRB or an estimator gap.

## 4.2 Exact campaign sample-complexity objects

For a task weight `W_B`, the CRB-predicted relative MSE is

```text
NMSE_CRB(T_eff)
 = tr(W_B J̄_B^-1)
   / [T_eff ||W_B^(1/2)β||^2].                                 (4.7)
```

The required number of independent-equivalent medium states for target error `ε` is

```text
T_req(ε)
 = tr(W_B J̄_B^-1)
   / [ε^2 ||W_B^(1/2)β||^2].                                   (4.8)
```

If `λ_j` and `a_j` are the per-effective-epoch Fisher eigenvalues and scene coefficients, the fraction of modes recoverable at Fisher SNR `γ` is

```text
f_rec(γ,T_eff)
 = (1/d_B) # {j : T_eff a_j^2 λ_j ≥ γ^2}.                      (4.9)
```

Every publication-scale prognosis must report (4.7)–(4.9), the full eigenvalue spectrum, and the profiling efficiency

```text
η_prof,j = λ_j(J_B) / λ_j(I_ββ),                                (4.10)
```

rather than one reconstructed image.

## 4.3 Photon visibility law

Let `n` be mean detected photoelectrons per physical bucket and define the covariance visibility parameter

```text
r = n σ_f^2 × geometry factor.                                  (4.11)
```

In a normalized scalar approximation,

```text
I_blind / T_eff  ∝ [r/(1+r)]^2,
I_oracle / T     ∝ r.                                           (4.12)
```

Hence:

- **shot-dominated:** `r<<1`, blind covariance information scales as `n^2 σ_f^4`;
- **medium-dominated:** `r>>1`, blind information per independent state saturates, whereas a realization-aware oracle continues to gain with photon flux;
- blind/oracle relative efficiency is largest near `r≈1` and decreases again at very high flux.

This is the fundamental blind/oracle tax. It is a rate gap, not an error floor. It also explains the empirical observation that `10^4` photons is already essentially clean: beyond that point, buy more independent medium states, not more photons.

## 4.4 Temporal correlation

For a covariance statistic from a Gaussian OU/AR(1) sequence with `ρ(h)=φ^|h|`, a first-order effective sample count is

```text
T_eff,2
 ≈ T / [1+2Σ_{h≥1}ρ(h)^2]
 = T(1-φ^2)/(1+φ^2).                                           (4.13)
```

The production code must use the exact Toeplitz Fisher matrix, not substitute (4.13), but (4.13) is the physical planning quantity. In continuous time it approaches total duration divided by the medium correlation time.

## 4.5 Fixed versus fresh codes

For covariance-design atoms `a`, repeated allocation gives

```text
J_design = Σ_a n_a J_a.                                        (4.14)
```

Reusing a bank concentrates all independent medium realizations on the same set of covariance quadratic forms; fresh codes spread the same budget over many one-sample forms. This can favor repetition, but it is **not a universal theorem**. A non-iid exact likelihood can exploit fresh designs without first estimating every covariance entry.

Therefore the “concentration principle” remains a preregistered empirical/design claim. The fresh arm must receive best-effort GLS or the exact non-iid marginal likelihood with identical lifting and regularization. Only that comparison may support fixed-bank superiority.

## 4.6 E8 estimator ruling

The current MLE must be replaced by the following data-only pipeline:

1. **Lifted covariance GLS.** Restrict to the Fourier aperture and solve

   ```text
   X0 = argmin_{X>=0}
        || W_C^(1/2)[C_hat - L(X)] ||_F^2 + λ tr(X),             (4.15)
   ```

   with the in-band mean constraints, where `X=xx^T`. This is convex in the lifted variable; use a certified low-rank factorization if the full SDP is too large.

2. **Spectral extraction.** Take the leading rank-one factor, impose realness/nonnegativity, and anchor scale with the mean channel.

3. **Band continuation.** Expand the solved object band from `1.0 k_p` to `1.25, 1.5, 1.8, 2.0 k_p` rather than opening all fine modes at once.

4. **Fisher-scored marginal refinement.** Optimize the exact declared marginal likelihood with a trust-region or Riemannian Gauss–Newton method while jointly profiling the finite-dimensional medium-law nuisance.

5. **One-step efficiency audit.** Starting from truth plus a CRB-sized perturbation is allowed only as a diagnostic: one Fisher-scoring step must attain the predicted local covariance. It is not a reconstruction arm.

If lifted initialization plus exact marginal refinement still stalls well above the profiled CRB, the method—not the information channel—fails its flagship bar.

---

# 5. Thin-screen validity and mismatch boundary

The declared home is an intensity-multiplicative thin screen at or optically conjugate to the object plane. A screen separated far enough to create nonlocal Fresnel mixing does not obey (1.3); the correct covariance operator becomes convolutive or field-coherent.

There is no universal statement such as “20% convolution is safe.” Let `c(x)` be the modeled covariance vector and `Δc` the true-model discrepancy. The first-order beyond-band bias is

```text
δβ
 ≈ J_B^-1 D_B^T W [I-Π_η] Δc,                                  (5.1)
```

so the tolerated physical mismatch is set jointly by mismatch size and the smallest profiled Fisher eigenvalue. Near the aperture edge, even a small operator error can dominate.

Define the whitened mismatch

```text
ε_geom = ||W^(1/2)Δc|| / ||W^(1/2)c||.                          (5.2)
```

The operational validity boundary is reached when the bias predicted by (5.1) exceeds either half a standard error or the frozen 25% degradation bar. The sealed probe must measure this boundary rather than invent a universal screen distance.

Required mismatch axes:

- multiplicative/convolutive mixture `α ∈ {0.05,0.10,0.20,0.30}`;
- screen defocus/Fresnel propagation;
- medium cutoff and spectral slope errors `±20%`;
- anisotropy and band rotation;
- Gaussian versus bounded/lognormal fields with matched covariance;
- temporal correlation time `×0.5` and `×2`;
- `10–20%` static transmission envelope;
- detector gain/shot-noise mismatch and dead-time clipping.

Primary robustness is claimed only through `α=0.10`; `α=0.20` and `0.30` define the failure map. At every mismatch, held-out covariance residuals must flag model failure before image quality collapses.

---

# 6. Publication-scale sealed probe

## 6.1 Frozen geometry

```text
scene grid:             64×64, N=4096
fixed signed codes:     M=128
physical exposures:     complementary pair for each signed code
pattern Fourier band:   |kx|,|ky| <= 5  (121 real band degrees)
medium primary band:    same cutoff k_w=k_p
support aperture:       |kx|,|ky| <= 10  (2.0× cutoff)
primary usable claim:   frequencies <= 1.8 k_p
edge audit:             1.8–2.0 k_p
outside control:        >2.0 k_p
medium RMS:             σ_f in {0.15,0.30,0.50}; primary 0.30
T_eff grid:             {256,512,1024,2048,4096}
photons/physical bucket:{10^3,10^4,10^5}; primary ceiling 10^4
```

Use band-limited incoherent random codes. Fourier-lattice codes are a frozen negative control, not a production option.

The primary medium field is a bounded, mean-one, positive band-limited process. Rejection/truncation may enforce `w∈[0.2,1.8]`, but it may not alter the declared Fourier subspace. A lognormal field is a mismatch arm because exponentiation creates spectral harmonics.

At a 20 kHz DMD rate, one epoch of `2M=256` physical exposures lasts `12.8 ms`. `T_eff=4096` corresponds to approximately `52 s` and

```text
2 M T_eff n ≈ 1.05×10^10
```

detected photoelectrons at `n=10^4`, before correlation inflation. These are the physical units to report. A bench realization would step or rotate a fine diffuser once per completed code bank and keep it quasi-static within the bank.

## 6.2 Sealed banks

1. **Theory/unit bank:** Fourier impulses, exact same-mean pairs, and analytic covariance checks. No tuning.
2. **Design bank:** separate codes/scenes used to select `λ`, optimizer tolerances, and the fixed code bank.
3. **Fisher prognosis bank:** fresh scenes used only to compute the prelaunch profiled Fisher spectrum and `T_req`.
4. **Confirmatory scene bank:** at least 24 unseen natural/texture scenes with preregistered annular energy, plus 16 synthetic beyond-band witnesses.
5. **Mismatch bank:** separate scenes, codes, laws, and random seeds for §5.
6. **Outside-aperture bank:** frequencies above `2k_p`, used only to measure hallucination/leakage.

No object, code bank, medium seed, or optimization seed crosses banks. Commit all hashes and thresholds before opening the confirmatory bank. **No repair round.**

## 6.3 Arms

All deployable arms receive identical physical exposures, detected photons, duration, pattern-band limit, and object.

1. **FIXED + marginal MLE:** production lifted initializer plus exact marginal refinement.
2. **FIXED + moment GLS:** transparent covariance estimator and inverse.
3. **FRESH + mean:** strongest admissible first-moment estimator; exact wall arm.
4. **FRESH + covariance:** best-effort non-iid GLS/marginal implementation with the same lift, priors, and computation budget.
5. **ORACLE medium:** true `w_t`; nondeployable ceiling.
6. **CLASSIC averaging:** average the dynamic medium, then perform conventional band-limited reconstruction.
7. **LATTICE covariance:** structured-code negative control.

The fresh-covariance arm is not allowed to remain first-pass if that disadvantages it. Its weighting and law profiling must be frozen on the design bank.

## 6.4 Primary metrics

- beyond-band NRMSE on `1.0–1.8 k_p`;
- edge NRMSE on `1.8–2.0 k_p`;
- outside-aperture leakage above `2.0 k_p`;
- Fourier amplitude correlation and phase/sign error;
- full-image PSNR and a high-frequency task score;
- profiled CRB, `f_rec`, and empirical MSE/CRB ratio;
- multistart solution spread and lifted rank-one ratio;
- covariance prediction residual on held-out epochs;
- exact photon/time/pattern ledger.

## 6.5 Binding bars

### P0 — law and wall

- mean-channel Fisher is numerically zero above `k_p` (`<10^-10` relative);
- covariance Fisher is zero outside `2k_p`;
- at least 95% of frequencies inside `1.8k_p` have positive profiled information;
- empirical single-frequency recovery map agrees with the aperture prediction on at least 95% of probes;
- same-mean beyond-band scene pairs are indistinguishable to every mean arm.

**Kill** on any violation; it indicates a theorem/implementation defect.

### P1 — Fisher viability before reconstruction

At `N=4096`, `M=128`, `σ_f=0.30`, `n=10^4`, and `T_eff≤4096`:

- no exact profiled null on the primary annulus;
- 10th-percentile profiling efficiency `η_prof >= 0.10` and median `>=0.25`;
- CRB-predicted beyond-band NRMSE `<=0.25` on synthetic witnesses and median `<=0.35` on the natural-scene bank;
- at least 70% of primary beyond-band modes have Fisher SNR `>=3`.

**Kill before launch** if this bar fails. This is the publication-scale width-floor analogue.

### P2 — estimator efficiency and no plateau

- empirical beyond-band MSE is at most `2×` the profiled CRB median;
- from `T_eff=512` to `4096`, error follows the Fisher-predicted downward trend, with no two-step plateau exceeding 10%;
- lifted initialization and at least 10 perturbed cold starts converge to solutions with pairwise beyond-band correlation `>=0.98` and NRMSE standard deviation `<=0.03`;
- the lifted solution has leading-eigenvalue fraction `>=0.90`.

**Kill** if the production estimator remains at the E8 plateau while the CRB predicts substantially better performance.

### P3 — made-possible-at-all imaging effect

At the primary resource point:

- median beyond-band NRMSE `<=0.30` on synthetic witnesses;
- median `<=0.40` and 90th percentile `<=0.55` on natural scenes;
- mean/averaging arms remain `>=0.95` on beyond-band NRMSE;
- production gives at least a `2×` reduction in beyond-band error relative to the mean wall;
- on scenes with at least 20% energy in the primary annulus, full-image PSNR improves by at least `1.5 dB` median with in-band loss `<=0.2 dB`.

**Kill flagship** if only oracle recovery is strong or if blind recovery remains a weak `0.4–0.5` toy effect without full-image consequence.

### P4 — fixed-bank concentration claim

FIXED + marginal MLE must improve median beyond-band NRMSE by at least 15% over the fully optimized FRESH + covariance arm and may not lose on more than 20% of confirmatory scenes.

Failure kills the fixed-bank concentration claim. If the aperture method otherwise passes, the result is demoted to a general covariance-aperture method rather than promoted as the fixed-code flagship.

### P5 — mismatch and geometry

- degradation `<=25%` for `α<=0.10` partial convolution;
- degradation `<=25%` for each frozen law/band/τ mismatch in its primary range;
- held-out residual diagnostics detect `α>=0.20` or other out-of-domain laws with at least 90% sensitivity at 5% false-positive rate;
- outside-aperture leakage remains below 10% of target-annulus energy.

**Kill flagship** if a small undetected geometry error creates apparent superresolution.

### P6 — scale and equal-ledger audit

- all deployable arms are byte-audited for exposures, photons, duration, and admissible pattern band;
- 32×32 and 64×64 trends agree after dimensionless rescaling;
- a sealed 128×128 spot check preserves the aperture map and Fisher scaling on at least four scenes;
- no truth warm start, measured speckle realization, reference detector, transmission matrix, or hidden calibration enters a deployable arm.

**Kill** on any ledger or scaling failure.

## 6.6 Kill tree

1. **P0 fails:** theorem/channel bug — stop.
2. **P1 fails:** information too thin at publication scale — stop before reconstruction.
3. **P2 fails:** estimator cannot extract the available information — archive theorem, no flagship.
4. **P3 fails:** effect is not large enough — kill mainline.
5. **P4 fails:** remove fixed-bank concentration claim; flagship survives only if the operator accepts the more occupied general covariance method.
6. **P5 fails:** thin-screen domain too fragile — kill flagship.
7. **P6 fails:** simulation/accounting artifact — stop.

No post-unblinding threshold change, scene removal, estimator switch, or repair campaign is authorized.

---

# 7. Prior-art fence and sim-only credibility

## 7.1 The mechanistic collision must be ceded in the first page

The closest mechanistic twin is **Chaigne et al., Optica 3, 54–57 (2016)**, not the 2017 flow-fluctuation paper:

- Chaigne et al., “Super-resolution photoacoustic fluctuation imaging with multiple speckle illumination,” DOI [10.1364/OPTICA.3.000054](https://doi.org/10.1364/OPTICA.3.000054).

It used unknown dynamic speckle illumination, second-order photoacoustic fluctuations, and a spatial ultrasound array to exceed the acoustic passband. The 2017 Optica paper instead exploited flow-induced absorber fluctuations: DOI [10.1364/OPTICA.4.001397](https://doi.org/10.1364/OPTICA.4.001397).

The general covariance superresolution capacity is also established:

- Mudry et al., blind SIM, DOI [10.1038/nphoton.2012.83](https://doi.org/10.1038/nphoton.2012.83).
- Labouesse et al., joint blind-SIM reconstruction, DOI [10.1109/TIP.2017.2675200](https://doi.org/10.1109/TIP.2017.2675200).
- Idier et al., covariance capacity of unknown-speckle imagers, DOI [10.1109/TCI.2017.2771729](https://doi.org/10.1109/TCI.2017.2771729).
- Yeh, Tian & Waller, covariance-prior SIM, DOI [10.1364/BOE.8.000695](https://doi.org/10.1364/BOE.8.000695).
- Labouesse et al., variance-based twofold-resolution analysis, DOI [10.1103/PhysRevA.109.033525](https://doi.org/10.1103/PhysRevA.109.033525).
- Dertinger et al., SOFI, DOI [10.1073/pnas.0907866106](https://doi.org/10.1073/pnas.0907866106).

The paper must say plainly:

> Unknown-fluctuation superresolution and covariance bandwidth extension are established. Our candidate contribution is the zero-dimensional bucket configuration and its exact separation into an impossible first-moment channel and a possible covariance channel relative to a separate, band-limited programmable modulator.

## 7.2 Defensible claim wording

> **For a band-limited code bank observed only by a single bucket, an unmeasured dynamic thin screen with a declared covariance law creates a second-order synthetic aperture: all first-moment code designs are exactly blind beyond the modulator band, whereas the bucket covariance carries scene frequencies over the pattern–medium Minkowski aperture. We characterize the profiled information rate and recover those frequencies without measuring any screen realization.**

This is a novelty candidate, not a certified priority claim.

## 7.3 Forbidden claims

- first use of unknown speckles for superresolution;
- first fluctuation/covariance image beyond a passband;
- first scattering-assisted or random-illumination imaging;
- first single-pixel or computational ghost imaging;
- first synthetic aperture from a scattering medium;
- arbitrary-medium or arbitrary-distance validity;
- “turbulence gives free information” without integration-time and model costs;
- twofold usable resolution merely because support reaches `2k_p`;
- experiment, experimental feasibility, or biological validation.

## 7.4 Sim-only defense by architecture

The manuscript is credible without an experiment only if it is built as a theorem-and-falsification paper:

1. **Exact law and exact wall first.** The aperture theorem and first-moment impossibility precede all reconstructed images.
2. **Prelaunch Fisher kill.** The paper shows the full profiled spectrum and predicts the recoverable fraction before displaying outcomes.
3. **Sealed campaign.** All bars, banks, code hashes, budgets, and failures are committed before confirmation.
4. **Independent checks.** Analytic covariance, Monte Carlo covariance, marginal likelihood, and lifted GLS are separate implementations.
5. **No simulation laundering.** State “simulation” in the abstract and every apparatus schematic. Do not draw a completed bench.
6. **Full adverse record.** Include the dead pathwise branch, fresh-pattern in-band kill, third-order failure, and lattice-code failure in the supplement.
7. **Cross-scale and mismatch maps.** The 64×64 primary and 128×128 spot check must demonstrate that the effect is not a 16×16 discretization artifact.

## 7.5 Falsifiable experimental transfer, not dependency

The Discussion may specify one future test:

- a DMD projects complementary, low-pass random codes onto a static transmissive target;
- a fine rotating diffuser is optically conjugate to, or placed sufficiently near, the target plane;
- the diffuser is held fixed for one `2M`-exposure code bank and advanced between banks;
- one integrating photodiode records the buckets;
- an optional camera or wavefront sensor records the diffuser only as validation truth and is excluded from reconstruction;
- remove the diffuser, enlarge its grains beyond the code band, or move it out of the multiplicative regime as three falsification controls.

Predictions fixed from simulation:

1. the mean reconstruction has exactly zero response to the beyond-band witness;
2. covariance recovery follows `B_p⊕B_w` and fails outside it;
3. usable edge moves with measured grain size;
4. the method fails as the screen becomes nonlocal/convolutive.

The manuscript does not depend on this experiment being performed.

---

# 8. Blunt judgment and manuscript consequence

## 8.1 Operator-bar judgment

> **Eye-lighter if—and only if—the sealed blind result clears P1–P3.**

The phrase “fog as a second DMD” is not itself new; unknown-speckle superresolution is occupied. What can still be exquisite is the exact channel inversion:

```text
mean channel:       mathematically impossible beyond k_p
covariance channel: recoverable over k_p+k_w from the same 0-D record
```

That is one compact protagonist, not a method anthology. It creates information that the declared hardware mean channel cannot carry at all.

The Fisher accounting does **not** justify killing the direction today:

- the exact covariance map has already shown injectivity;
- shot noise is not the dominant price once `n σ_f^2 >> 1`;
- there is no asymptotic `0.4` floor if the profiled derivative is full rank;
- the current plateau has a credible lifted-initialization remedy.

But this is the last allowed optimism. If the 64×64 profiled CRB or blind sealed result fails the frozen bars, kill the flagship. A visually suggestive `0.45` reconstruction at 16×16 is not enough.

## 8.2 One-manuscript ruling

If P0–P6 pass, this becomes **the** flagship and replaces the R36 fail-branch manuscript. Do not split it into a letter plus a map paper.

Frozen `≤7 pp / ≤4 figure` seed:

1. **DEFINE — One bucket and two statistical channels.** Thin-screen model, 0-D record, declared law.
2. **LAW — The covariance aperture.** Exact Minkowski support plus Fisher-weighted usable aperture.
3. **WALL — What no admissible mean measurement can see.** First-moment impossibility and head-on Chaigne/Idier fence.
4. **METHOD — Statistics-only recovery beyond the modulator band.** Lifted marginal estimator, profiled Fisher, sealed 64×64 results, mismatch boundary.

Main figures:

1. physical concept and pattern/medium/aperture bands;
2. exact aperture map, mean wall, and Fisher-weighted edge;
3. equal-budget arm comparison and blind-versus-CRB/oracle scaling;
4. natural-scene outcomes, mismatch domain, and physical budget.

M1 appears only as a short design-window ceiling: detector bandwidth/dead-time limits the rate at which independent covariance samples can be acquired. All fixed-operator negatives, DLGI, pathwise fog failure, MOLT, and audit detail remain supplementary.

If the sealed probe fails, do not force a standalone simulation letter. Preserve the aperture and impossibility theorems in the materials bank and return to the already authorized law/boundary fallback.