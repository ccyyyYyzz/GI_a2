# R33 ruling (GitHub issue #25, raw)

Title: R33
Posted: 2026-07-22T13:53:05Z

---

# R33 — Dual-ledger imaging: exact Fisher partition, medium identifiability, and flagship ruling

**Reference request:** [`docs/ROUND63_GPT_ROUND33_QUESTION.md`](https://github.com/ccyyyYyzz/GI_a2/blob/833c60d/docs/ROUND63_GPT_ROUND33_QUESTION.md) at commit [`833c60d`](https://github.com/ccyyyYyzz/GI_a2/commit/833c60d7022703a2a7b94de0a442e8e99ae0fba7).  
**Predecessor evidence:** MOLT finishing verdict at [`cf3b3db`](https://github.com/ccyyyYyzz/GI_a2/commit/cf3b3db57d2bd7b05e96ea66521a57588d68c5b7).

## Executive verdict

> **CONDITIONAL GO to one decisive feasibility campaign; not yet a flagship GO.** The same unmodified bucket record can in principle support a scene estimate and a medium-dynamics estimate. However, the proposed scalar “two-ledger conservation law” is false, and two advertised interpretations must be reversed.

The theorem that survives is cleaner:

1. the joint observed Fisher matrix equals a block-diagonal complete-data ledger minus one positive-semidefinite posterior-score covariance;
2. scene and medium are coupled by a cross/confusion block;
3. after profiling either ledger, the **same canonical-confusion singular values** quantify its normalized information loss;
4. the O4-A moment condition annihilates that confusion and therefore protects **both** ledgers—it does not protect the scene by blinding the medium.

The scene missing-information term

```text
M^T E[Cov(a | Y)] M
```

is **not** the medium Fisher information and is not equal to the medium’s remaining uncertainty. It is a scene-space image of posterior path uncertainty. The medium ledger is formed from the posterior expectation of the medium-parameter score. They share the same posterior, but they are different statistical objects with different dimensions and units.

The present concept is therefore an elegant, potentially publishable dual-use measurement architecture, but not yet the operator’s universal sledgehammer. It earns that status only if the running probe shows, over a nontrivial `(T/t_c, Delta/t_c, CV, photon-SNR)` region, that pilot-free `t_c/CV` precision is within `1.5 x` of the pilot baseline while scene performance is noninferior under an exactly matched photon/time/pattern ledger.

---

# 1. Exact partition theorem

## 1.1 Separated hierarchical model

Let

```text
x      = scene parameter,
theta  = medium-process parameter, e.g. (t_c, CV),
a      = latent gain path,
Y      = observed bucket record.
```

Assume the separated hierarchy

```text
p(y,a | x,theta) = p(y | x,a) p_theta(a),                       (1.1)
```

so `theta` controls the latent medium law and `x` controls the optical observation conditional on the realized path. Assume the usual dominated-differentiability and finite-score conditions.

Define complete-data scores

```text
U_x     = grad_x log p(Y | x,a),
U_theta = grad_theta log p_theta(a).
```

Because `E[U_x | a]=0`, the complete-data cross information vanishes:

```text
E[U_x U_theta^T] = 0.                                           (1.2)
```

Write

```text
C_x     = E[U_x U_x^T],
C_theta = E[U_theta U_theta^T]
```

for the oracle/complete-data ledgers.

## 1.2 Theorem 1 — joint two-ledger missing-information identity

Let

```text
U = [U_x ; U_theta],
L = E_Y Cov(U | Y).
```

Then

```text
I_Y(x,theta)
  = diag(C_x,C_theta) - L                                      (1.3)
```

and, blockwise,

```text
I_xx         = C_x     - L_xx,
I_thetatheta = C_theta - L_thetatheta,
I_xtheta     =          - L_xtheta.                             (1.4)
```

Here

```text
L = [ L_xx       L_xtheta
      L_thetax   L_thetatheta ] >= 0.                           (1.5)
```

### Proof

Fisher’s identity gives

```text
score_x(Y)     = E[U_x | Y],
score_theta(Y) = E[U_theta | Y].                                (1.6)
```

Apply the law of total covariance to the stacked complete score:

```text
Cov(U)
  = Cov(E[U | Y]) + E Cov(U | Y).                               (1.7)
```

The first term on the right is the observed joint Fisher matrix. By (1.2), the left side is block diagonal with blocks `C_x` and `C_theta`. Rearrangement proves (1.3)–(1.5). `square`

This is Louis/Orchard–Woodbury missing-information theory specialized to the dual scene/medium hierarchy. The specialization, not the abstract identity, is the possible contribution.

## 1.3 Poisson GI specialization—and the exact correction to the teaser

For

```text
s = Mx,
Y_n | a,x ~ Poisson(a_n s_n),
```

the complete scene score is

```text
U_x = M^T D_s^{-1}Y - M^T a.
```

Therefore

```text
L_xx = M^T E_Y[Cov(a | Y,x,theta)] M,                            (1.8)
```

recovering the R26 HSGI identity.

But the medium score is

```text
U_theta = grad_theta log p_theta(a),
```

so

```text
I_thetatheta
  = C_theta
    - E_Y Cov(grad_theta log p_theta(a) | Y).                   (1.9)
```

Equations (1.8) and (1.9) are not equal. Even their dimensions generally differ. The strongest valid wording is:

> The scene-information loss and the medium-information deficit are two projections of the same posterior ambiguity over the latent path, coupled through the posterior cross-covariance of their complete scores.

Because `L >= 0`, the generalized Schur complement gives the useful inequality

```text
L_xtheta L_thetatheta^dagger L_thetax <= L_xx,                 (1.10)
```

with the usual range condition, and symmetrically in the other direction. Thus cross-confusion is bounded by the two diagonal missing-information ledgers; there is no equality in general.

---

# 2. The true complementarity law: canonical-confusion reciprocity

Let the observed joint Fisher matrix for one schedule `pi` be

```text
I(pi) = [ A  C
          C^T B ],                                              (2.1)
```

where

```text
A = I_xx       : scene information if theta were known,
B = I_thetatheta: medium information if x were known,
C = I_xtheta   : scene/medium score confusion.
```

Assume `A>0`, `B>0`; use pseudoinverses and identifiable subspaces otherwise.

The efficient information ledgers after treating the other parameter as nuisance are

```text
J_x     = A - C B^{-1} C^T,
J_theta = B - C^T A^{-1} C.                                    (2.2)
```

Define the whitened cross-information operator

```text
K = A^{-1/2} C B^{-1/2}.                                       (2.3)
```

## Theorem 2 — canonical-confusion ledger reciprocity

The efficient ledgers satisfy

```text
J_x
  = A^{1/2}(I - K K^T)A^{1/2},
J_theta
  = B^{1/2}(I - K^T K)B^{1/2}.                                 (2.4)
```

Consequently the nonzero normalized loss eigenvalues in the two ledgers are identical: they are the squared singular values

```text
kappa_1^2,...,kappa_r^2
```

of `K`. In particular,

```text
tr[A^{-1/2}(A-J_x)A^{-1/2}]
  = ||K||_F^2
  = tr[B^{-1/2}(B-J_theta)B^{-1/2}],                            (2.5)
```

and

```text
det(J_x)/det(A)
  = det(I-KK^T)
  = det(I-K^T K)
  = det(J_theta)/det(B)
  = product_i (1-kappa_i^2).                                   (2.6)
```

The determinant statement extends with pseudo-determinants on the common identifiable subspaces.

### Proof

Substitute (2.3) into the two Schur complements to obtain (2.4). The matrices `KK^T` and `K^T K` have the same nonzero eigenvalues, proving (2.5). Sylvester’s determinant identity gives

```text
det(I-KK^T)=det(I-K^TK),
```

proving (2.6). `square`

### Meaning

This is the honest “two-ledger” theorem:

> The same canonical confounding modes tax both efficient ledgers by the same normalized factors.

It is **reciprocity**, not additive conservation. There is no meaningful identity of the form

```text
scene information + medium information = constant,
```

because the ledgers have different coordinates, dimensions and units. Nor is `||K(pi)||` conserved along a schedule family. The schedule generally changes all three objects `(A(pi),B(pi),K(pi))`.

A particularly important consequence is the opposite of the proposed simple dial:

> If a schedule reduces `K` while preserving `A` and `B`, it improves both scene and medium inference simultaneously. There is no forced scene-versus-medium trade.

A genuine trade appears only when a schedule changes the oracle ledgers themselves—for example, increasing temporal excitation and hence `B` while sacrificing spatial diversity and hence `A`. The correct schedule object is therefore a Pareto surface such as

```text
maximize_pi  lambda logdet J_x(pi)
           + (1-lambda) logdet J_theta(pi)                      (2.7)
```

under the same time, photon and pattern-multiset constraints.

---

# 3. What O4-A really does

In the R26 local Gaussian path model

```text
y = M_pi x + D_s H beta + epsilon,
epsilon ~ N(0,R),
s = M_pi x,
```

the blocks are

```text
A = M_pi^T R^{-1} M_pi,
C = M_pi^T R^{-1} D_s H,
B = H^T D_s R^{-1} D_s H + Lambda_beta.                         (3.1)
```

Thus the O4-A moment condition

```text
M_pi^T R^{-1}D_sH = 0                                           (3.2)
```

is exactly `C=0`, hence `K=0`. Equations (2.4) then give

```text
J_x = A,
J_beta = B.                                                     (3.3)
```

So nuisance-orthogonal scheduling makes the path coefficients and the scene locally separable and lets **both** attain their oracle Fisher information. It does not blind the medium.

A paired/differenced analysis can nevertheless destroy medium information if it discards the common-mode/raw bucket component. That is a data-processing loss caused by the reduction, not a consequence of O4-A. The campaign must retain and analyze the full raw record.

Likewise, an “anti-paired” schedule is not automatically medium-optimal. It may increase the oracle temporal excitation `B`, but it can also increase cross-confusion `C`. It helps the medium only when the increase in `B` exceeds the Schur penalty `C^T A^{-1}C`.

Finally, (3.2) is a path-coefficient statement. For hyperparameters `theta=(t_c,CV)` after integrating out the path, it need not zero the full hyperparameter cross block. In the Gaussian marginal model

```text
Y | x,theta ~ N(s, Sigma),
Sigma = R + D_s C_theta D_s,
```

the exact Fisher formula is

```text
I_uv
 = (partial_u s)^T Sigma^{-1}(partial_v s)
   + 1/2 tr(Sigma^{-1} partial_u Sigma
            Sigma^{-1} partial_v Sigma).                       (3.4)
```

Since `Sigma` depends on `x` through `D_s`, the `x-theta` block can survive even when the local path-coefficient moment condition is met. The production schedule dial must therefore be evaluated with the full joint Fisher matrix for `(x,t_c,CV)`, not inferred from (3.2) alone.

---

# 4. Medium-channel identifiability

## 4.1 Exact local criterion

After fixing gauges, the medium parameter is locally identifiable in the presence of an unknown scene if and only if its efficient Fisher matrix is positive definite on the claimed parameter subspace:

```text
J_theta
  = I_thetatheta - I_thetax I_xx^dagger I_xtheta > 0.          (4.1)
```

This is the clean test for the feasibility probe. High gain-path correlation is not a substitute for (4.1).

Necessary physical/statistical conditions are:

1. **Scale gauge fixed.** The transformation `(x,a)->(c x,a/c)` leaves a multiplicative bucket model unchanged. Fix `E[a]=1`, mean log-gain zero, or an equivalent throughput normalization. Absolute mean transmission cannot otherwise be separated from scene scale.
2. **Scalar-gain model valid.** If the medium changes the spatial transfer operator rather than one common frame gain, `t_c` and `CV` belong to a higher-dimensional transfer process and the proposed two-parameter ledger is misspecified.
3. **Temporal span and cadence resolve the process.** The useful regime is roughly `Delta <= t_c << T=N Delta`.
4. **Nonzero fluctuation visibility.** At `CV=0`, `t_c` has no physical meaning; with finite bucket noise, both parameter informations vanish as the gain spectrum falls below the noise spectrum.
5. **Scene tangent cannot absorb the temporal law.** Pattern ordering and scene dimension must leave temporal residual contrasts after profiling `x`; this is exactly (4.1).

## 4.2 Oracle OU/AR(1) Fisher scaling

Let a normalized small log-gain/deviation process be sampled every `Delta` and obey a stationary Gaussian AR(1) law with

```text
Var(a_n-1) = CV^2,
phi = exp(-Delta/t_c),
eta = log CV.
```

For a directly observed high-SNR path and large `N`, the per-transition Fisher matrix for `(eta,phi)` is

```text
I/N = [ 2                         -2phi/(1-phi^2)
        -2phi/(1-phi^2)  (1+phi^2)/(1-phi^2)^2 ].               (4.2)
```

Profiling the other parameter yields

```text
Var(phi_hat) >= (1-phi^2)/N,                                   (4.3)

Var(log CV_hat)
  >= (1+phi^2)/[2N(1-phi^2)].                                  (4.4)
```

With `r=Delta/t_c`, transformation of (4.3) gives

```text
Var(t_c_hat)/t_c^2
  >= [exp(2r)-1]/(N r^2).                                      (4.5)
```

In the dense-sampling, long-span regime `Delta << t_c << T=NDelta`,

```text
sd(t_c_hat)/t_c  >= sqrt(2 t_c/T),                              (4.6)

sd(CV_hat)/CV    >= sqrt(t_c/(2T)).                             (4.7)
```

These are useful oracle floors. For example, even with the path directly observed, approximately 10% relative precision requires about

```text
T >= 200 t_c  for t_c,
T >=  50 t_c  for CV.                                          (4.8)
```

The limits are diagnostic:

- **fast drift, `t_c << Delta`:** `phi` is nearly zero, and (4.5) blows up exponentially; sampled buckets can retain variance/CV information but lose correlation-time information;
- **slow drift, `t_c >> T`:** one acquisition contains fewer than one effective correlation length, so stationary asymptotics fail and neither `t_c` nor ensemble CV is consistently estimable from that record;
- **small CV:** the path-oracle scale family in (4.2) hides the physical loss. With additive/shot noise, the gain spectrum disappears below the measurement spectrum and information tends to zero.

## 4.3 Bucket/noise penalty

After scene normalization, an approximate stationary residual has spectrum

```text
f_Y(omega)=f_a(omega;theta)+f_noise(omega).
```

The long-record Whittle information is

```text
[I_theta]_{rs}
  ~= N/(4pi) integral
      partial_r log f_Y(omega)
      partial_s log f_Y(omega) domega.                          (4.9)
```

Medium derivatives are attenuated by the spectral visibility

```text
W(omega)=f_a(omega)/[f_a(omega)+f_noise(omega)].                (4.10)
```

Thus each spectral contribution carries approximately `W(omega)^2`. In the weak-fluctuation regime `f_a proportional CV^2`, the information for correlation-shape and log-amplitude parameters can scale as `N CV^4` before scene profiling, so standard errors can grow as `1/(CV^2 sqrt(N))`. Unknown `x` then imposes the additional Schur penalty in (4.1).

This is why the committed gain-path correlation `0.998` is encouraging but not decisive. A smoother can correlate almost perfectly with a low-frequency path while producing biased or poorly covered estimates of `t_c` and `CV`. The probe must report parameter bias, standard error, interval coverage and efficient Fisher ratios—not path correlation alone.

---

# 5. Ruthless prior-art fence

## 5.1 Missing information and parameter orthogonality

The abstract decomposition and nuisance orthogonality are classical:

- Louis, “Finding the Observed Information Matrix When Using the EM Algorithm,” DOI [10.1111/j.2517-6161.1982.tb01203.x](https://doi.org/10.1111/j.2517-6161.1982.tb01203.x).
- Cox & Reid, “Parameter Orthogonality and Approximate Conditional Inference,” DOI [10.1111/j.2517-6161.1987.tb01422.x](https://doi.org/10.1111/j.2517-6161.1987.tb01422.x).

**Forbidden:** “first information partition between a parameter and nuisance” or “first nuisance-orthogonal design.”

**Possible residual:** the exact programmable-GI block specialization, the canonical-confusion reciprocity as the campaign’s organizing theorem, and its experimentally validated schedule consequences.

## 5.2 DWS, DCS and speckle correlometry — closest optical collision

Medium dynamics from temporal intensity correlations are foundational:

- Pine et al., “Diffusing-wave spectroscopy,” DOI [10.1103/PhysRevLett.60.1134](https://doi.org/10.1103/PhysRevLett.60.1134).
- Boas, Campbell & Yodh, “Scattering and Imaging with Diffusing Temporal Field Correlations,” DOI [10.1103/PhysRevLett.75.1855](https://doi.org/10.1103/PhysRevLett.75.1855).
- Durduran et al., “Diffuse optical measurement of blood flow, blood oxygenation, and metabolism in a human brain…,” DOI [10.1364/OL.29.001766](https://doi.org/10.1364/OL.29.001766).
- Bandyopadhyay et al., “Speckle-visibility spectroscopy,” DOI [10.1063/1.2037987](https://doi.org/10.1063/1.2037987).

DCS/DWS already use one or a few photodetection channels and intensity autocorrelation to estimate decorrelation/flow parameters. Therefore:

**Forbidden:** first use of a single-pixel intensity record to estimate medium dynamics; first optical measurement of `t_c`; first use of speckle fluctuations as a medium sensor.

**Surviving distinction:** a programmable GI sequence encodes an unknown scene while the same unmodified bucket stream is jointly decoded into a scene product and a medium-process product, with no pilot/reference channel and with a joint Fisher certificate quantifying confusion between the products.

## 5.3 Dynamic-medium GI correction and the pilot baseline

Relevant direct neighbors include:

- Yang et al., “Noise reduction in computational ghost imaging by interpolated monitoring,” DOI [10.1364/AO.57.006097](https://doi.org/10.1364/AO.57.006097).
- Xiao, Zhou & Chen, “High-resolution ghost imaging through complex scattering media via a temporal correction,” DOI [10.1364/OL.463897](https://doi.org/10.1364/OL.463897).
- Peng & Chen, “Learning-based correction with Gaussian constraints for ghost imaging through dynamic scattering media,” DOI [10.1364/OL.499787](https://doi.org/10.1364/OL.499787).
- Zhou, Xiao & Chen, “High-resolution self-corrected single-pixel imaging through dynamic and complex scattering media,” DOI [10.1364/OE.489808](https://doi.org/10.1364/OE.489808).

These establish pilot/interpolated monitoring, temporal correction, learning-based correction and dual-detector self-correction.

**Surviving distinction:** no inserted assay pattern, no second detector, and treating the recovered medium parameters as a certified science output rather than merely a correction variable. The pilot-interleaved method remains the mandatory strongest baseline.

## 5.4 Radio self-calibration solutions as science products

The broad “calibration nuisance becomes a medium measurement” idea is already mature in radio astronomy:

- Intema et al., SPAM ionospheric calibration, DOI [10.1051/0004-6361/200811094](https://doi.org/10.1051/0004-6361/200811094).
- Mevius et al., “Probing ionospheric structures using the LOFAR radio telescope,” DOI [10.1002/2016RS006028](https://doi.org/10.1002/2016RS006028).
- Loi et al., “Real-time imaging of density ducts between the plasmasphere and ionosphere,” DOI [10.1002/2015GL063699](https://doi.org/10.1002/2015GL063699).

**Forbidden:** first conversion of calibration/gain solutions into a physical-medium science product; first simultaneous instrument/propagation inference and imaging from astronomical data.

## 5.5 Adaptive-optics telemetry as turbulence instrumentation

AO systems already reuse wavefront-sensor/control telemetry to estimate atmospheric turbulence and wind:

- Laidlaw et al., “Optimizing the accuracy and efficiency of optical turbulence profiling using adaptive optics telemetry…,” DOI [10.1093/mnras/sty3285](https://doi.org/10.1093/mnras/sty3285).
- Laidlaw et al., “Automated wind velocity profiling from adaptive optics telemetry,” DOI [10.1093/mnras/stz3062](https://doi.org/10.1093/mnras/stz3062).

**Forbidden:** first reuse of imaging/control telemetry to identify the propagation medium with no dedicated monitor.

## 5.6 Blind calibration and joint signal/gain recovery

Simultaneous object and gain estimation is a crowded inverse-problem class:

- Ling & Strohmer, “Self-Calibration and Bilinear Inverse Problems via Linear Least Squares,” DOI [10.1137/16M1103634](https://doi.org/10.1137/16M1103634).
- Cambareri & Jacques, “Through the haze: a non-convex approach to blind gain calibration…,” DOI [10.1093/imaiai/iay004](https://doi.org/10.1093/imaiai/iay004).

**Forbidden:** first blind joint scene/gain recovery or first pilot-free bilinear calibration.

## 5.7 What is actually left

I found no exact prior instance of the full proposed square:

> **one ordinary programmable GI bucket record, with no pilot/reference measurement, simultaneously delivering an unknown-scene reconstruction and quantitative dynamic-medium hyperparameters, accompanied by an exact joint-information/confusion theorem and a schedule chosen to control the two efficient ledgers.**

That is a defensible **novelty candidate**, not a certified novelty claim. Its novelty is compositional and theorem-backed; none of its abstract ingredients is individually new.

---

# 6. Blunt flagship judgment

At present this is **competent-to-elegant integration, not yet the operator’s sledgehammer**.

Why it is not yet universal:

- the scalar per-frame gain approximation is a significant physical restriction;
- `t_c` is unidentifiable in both the undersampled-fast and one-realization-slow limits;
- medium information can vanish as `CV` falls below bucket noise;
- the core missing-information and orthogonality machinery is classical;
- DCS, radio selfcal science and AO telemetry already own the broad dual-use-sensor narrative.

Why it still has a real chance to become an eye-lighter:

- “one bucket stream, two certified products” is immediately understandable;
- the canonical-confusion reciprocity gives the story a compact mathematical spine;
- no extra photons/exposures/pilots is a genuinely strong systems fact **if** the schedule is only a permutation of the same pattern multiset and the full acquisition ledger is identical;
- near-pilot medium precision with no scene penalty would be a surprising empirical result in GI, not merely a reframing.

### Promotion bar

Promote to flagship only if the probe shows all of the following on a predeclared interior grid, not one favorable point:

1. identical photons, duration, exposures and pattern multiset for pilot-free schedule comparisons; only order and inference may change;
2. scalar-gain model and scale gauge pass held-out residual tests;
3. pilot-free RMSE/CRB for both `t_c` and `CV` is at most `1.5 x` the strongest pilot-interleaved baseline and reasonably close to the oracle monitor;
4. calibrated interval coverage is correct;
5. scene performance is noninferior—suggested bar: no more than `0.2 dB` PSNR loss and no more than 5% task/Fisher loss;
6. measured schedule behavior is predicted by `(A,B,K)`, including canonical-confusion singular values, rather than merely ordered post hoc;
7. edge regimes are shown as honest failure regions.

If those bars land, the sharpest honest headline is:

> **One bucket stream, two certified products: pilot-free ghost imaging reconstructs the scene and estimates the medium’s correlation time and fluctuation strength within 1.5× of pilot precision, with no extra measurements and no measurable scene penalty.**

A shorter campaign line is:

> **Ghost imaging becomes its own medium monitor.**

Do not use “information conservation,” “free medium information,” or “paired schedules trade scene information into medium information.”

---

# 7. Names

My ranking:

1. **MIRROR-GI — Medium Identification and Reconstruction from a Reused Observation Record**  
   Best flagship-facing name. It captures one record serving two scientific products.

2. **DLGI — Dual-Ledger Ghost Imaging**  
   Best literal technical name; safe and easy to search.

3. **DUET-GI — Dual-Use Estimation and Tomography in Ghost Imaging**  
   Best presentation name if “ledger” feels too accounting-heavy.

Recommended theorem name:

> **Canonical-Confusion Ledger Reciprocity Theorem**

The method name and theorem name should remain distinct. The first is the experimental architecture; the second is the exact information geometry.

---

# Final ruling

> **Continue the running feasibility probe, but freeze the corrected theorem now.** The bucket record may support two ledgers, yet there is no additive two-ledger conservation law. The exact invariant within each schedule is reciprocal normalized loss through the same canonical confounding modes. O4-A orthogonality decouples and protects both ledgers. The program becomes flagship-worthy only if that geometry yields broad, near-pilot medium precision at genuinely zero acquisition cost and with scene noninferiority. Otherwise preserve it as a strong theorem-backed dual-use GI method, not the universal mainline.