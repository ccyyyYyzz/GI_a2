# R42 innovation round (GitHub issue #33, raw)

Title: R42
Posted: 2026-07-24T02:50:10Z

---

# R42 — PRL innovation round: information jets, curvature-rescued detection, and the scrambling inversion

**Reference brief:** [`docs/ROUND63_GPT_ROUND42_QUESTION.md`](https://github.com/ccyyyYyzz/GI_a2/blob/ed7a1e0/docs/ROUND63_GPT_ROUND42_QUESTION.md) at commit [`ed7a1e0`](https://github.com/ccyyyYyzz/GI_a2/commit/ed7a1e0).  
**New endpoint:** `results/round63_next/SCRAMBLE_EXT/` at the same commit.  
**Round type:** pure co-theory and invention. No campaign adjudication is made here.

## Executive synthesis

The failed R41 arithmetic–harmonic “duality” points to the right methodological correction. A spectral inequality at fixed trace is not a physical principle when the trace and anisotropy co-vary across realizable media. The next object must be invariant under reparameterization and nuisance profiling, and it must classify the **statistical law itself**, not one eigenvalue summary evaluated after physics has already selected the spectrum.

The deepest surviving object is the **quotient information jet**: the order of first statistical distinguishability of a physical perturbation after minimizing over every admitted nuisance. It yields a PRL-shaped effect:

> **Curvature-rescued detection.** A perturbation can be invisible to the mean and to every first-order response, yet remain detectable through a positive second statistical variation. Generic changes then require observation time proportional to `epsilon^-2`; first-order-orthogonal changes are rescued at curvature order and require `epsilon^-4`; a true blind direction is one for which every quotient jet vanishes.

This is not a restatement of “detection is easier than imaging.” It is a classification of physical observability by **contact order**. It survives the physical covariation that killed R41 because it asks which derivative of the *profiled probability law* first leaves the nuisance orbit. Thin fog and complete scrambling are different ranks and supports, but the same quotient-jet language covers both.

There is also a necessary correction to the current fully-scrambling slogan. For `Q(x)=x^T Gx`,

```text
Delta Q = 2 x^T G delta + delta^T G delta.
```

`delta^T G delta>0` makes every nonzero direction locally visible at first or second order when `G>0` and the medium law is fixed. It does **not** imply global injectivity for arbitrary signed finite changes: the ellipsoid `Q(x+delta)=Q(x)` contains nontrivial points. Nor does zero-lag covariance distinguish `Q(x)` from a free medium-amplitude nuisance, because their covariance derivatives are collinear. The strong theorem is therefore a **local quotient theorem**, with a global sign-definite corollary on monotone additive-change cones or after adding noncollinear lag/wavelength anchors.

My PRL candidate sentence is:

> **Strong disorder can erase every image coordinate while leaving every local scene-change direction statistically observable; the price of observing a direction is set by the order of the first nonzero quotient information jet.**

---

# 1. Mandatory move 1 — invariant hunt, derivation first

## 1.1 Why the R41 rank-1 idea failed

For a positive Fisher spectrum `lambda_1,...,lambda_d`,

```text
[tr(J)/d] [tr(J^dagger)/d] >= 1
```

is exact, with equality at isotropy. But it supplies no physical ordering unless the trace is held fixed. In the actual 81-cell medium family, total supply and anisotropy co-vary; the measured `corr(P,T_det)=+0.505` is opposite the intended causal story, and lower-tail mass did not collapse the physical map (`R^2=0.085`).

The lesson is binding:

> **No algebraic comparison at an artificially fixed resource may be promoted unless that resource is either physically controlled or eliminated by an invariant quotient.**

I therefore reject an “imaging–detection uncertainty relation” as the central object. It is a useful conditional inequality, not the mechanism.

## 1.2 Common statistical model

Let one independent medium bank generate record `Y` with law

```text
P_{x,eta},
```

where `x` is the scene and `eta` contains the admitted medium law, throughput, temporal, and detector nuisances. Consider a physical scene path

```text
x_epsilon = x + epsilon delta.
```

The nuisance-free question is not whether a selected observable moves. It is whether the perturbed law leaves the nuisance orbit

```text
N_x = {P_{x,eta'} : eta' admissible}.
```

Define the **profiled Chernoff distance**

```text
C_*(epsilon;delta)
 = inf_{eta'} C(P_{x+epsilon delta,eta}, P_{x,eta'}),            (1.1)
```

and analogously the profiled Kullback–Leibler distance `D_*(epsilon;delta)`. The infimum makes the object invariant to nuisance coordinates, gain gauges, and changes of medium basis.

If

```text
C_*(epsilon;delta)
 = c_m(delta) epsilon^(2m) + o(epsilon^(2m)),
 c_m(delta)>0,                                                  (1.2)
```

call `m=m(delta)` the **quotient information-jet order**. Set `m=infinity` if the profiled distance vanishes identically near zero.

For `T` independent-equivalent banks, Chernoff theory gives

```text
-log P_e^*(T,epsilon) = T C_*(epsilon;delta) + o(T),             (1.3)
```

so fixed-power detection obeys

```text
T_req(epsilon,delta) proportional to epsilon^(-2m(delta)).      (1.4)
```

This order is the first candidate invariant. It is a property of the physical experiment after nuisance quotienting; it does not assume fixed trace, isotropy, or a chosen estimator.

## 1.3 Fluctuation-response invariant and equality case

For a differentiable law `P_epsilon` with efficient score `s_eff`, every centered statistic `A` satisfies

```text
[d/d epsilon E_epsilon A |_{0}]^2
 <= Var_0(A) I_eff,                                             (1.5)
```

by

```text
partial_epsilon E A = Cov(A,s_eff)
```

and Cauchy–Schwarz. Equality holds if and only if

```text
A - E A = c s_eff                                               (1.6)
```

almost surely. Thus the nuisance-profiled matched score is not merely convenient: it is the unique local statistic that converts fluctuation variance into maximal response. This is the statistical-optics transplant of the fluctuation-response inequality of Dechant and Sasa ([DOI 10.1073/pnas.1918386117](https://doi.org/10.1073/pnas.1918386117)).

The identity does **not** say that more physical fluctuation always helps. It says that, for a fixed law, the efficient score saturates the response–variance bound. The physical law still controls `I_eff`.

## 1.4 An exact information-conductance sum rule

Consider a zero-mean Gaussian covariance channel

```text
Sigma(q)=R+qH,
R>0, H>=0,
```

and parameter `theta=log q`. Let `a_j` be the eigenvalues of

```text
A=R^(-1/2) qH R^(-1/2),
```

and define visibility eigenchannels

```text
tau_j = a_j/(1+a_j) in [0,1).                                  (1.7)
```

The per-bank Fisher information is exactly

```text
I_theta
 = 1/2 tr[(Sigma^-1 qH)^2]
 = 1/2 sum_j tau_j^2.                                          (1.8)
```

### Proof

`Sigma=R^(1/2)(I+A)R^(1/2)`, so `Sigma^-1 qH` is similar to `(I+A)^-1 A`, whose eigenvalues are `tau_j`. The Gaussian covariance Fisher formula gives (1.8).

### Consequences

```text
0 <= I_theta < rank(H)/2.                                      (1.9)
```

Each covariance eigenchannel carries at most one-half nat of local Fisher information per independent bank for `log q`. The upper bound is approached only when every visible channel is medium dominated (`a_j -> infinity`). This is the exact saturation behind the empirical facts that photons and fluctuation contrast become inert while independent banks remain valuable.

Equation (1.8) is a Landauer-like **fluctuation information conductance**: the ordinary transmission probabilities are replaced by squared covariance visibilities. Its algebra is classical Gaussian Fisher theory; the new opportunity is its physical interpretation and verification across the thin-screen-to-scrambler family. The closest disordered-systems relatives are Landauer eigenchannels and universal conductance fluctuations ([DOI 10.1103/PhysRevLett.55.1622](https://doi.org/10.1103/PhysRevLett.55.1622)).

## 1.5 Chernoff exponent and statistical order

For the full-scramble scalar channel

```text
Q(x)=x^T Gx,
Sigma(x)=R+Q(x)H,
```

set

```text
c=x^T G delta,
d=delta^T G delta,
z(epsilon)=Q(x+epsilon delta)-Q(x)=2c epsilon+d epsilon^2.       (1.10)
```

Whiten `H` as

```text
K=Sigma_0^(-1/2) H Sigma_0^(-1/2).
```

The exact Gaussian KL divergence is

```text
D(P_epsilon || P_0)
 = 1/2 sum_j [z kappa_j - log(1+z kappa_j)]                     (1.11)
 = z^2 ||K||_F^2/4 + O(z^3).                                   (1.12)
```

The Bhattacharyya/Chernoff distance has leading term

```text
C(P_epsilon,P_0)
 = z^2 ||K||_F^2/16 + O(z^3).                                  (1.13)
```

Therefore:

```text
c != 0: C = [c^2 ||K||_F^2/4] epsilon^2 + o(epsilon^2), m=1;
c  = 0, d>0: C = [d^2 ||K||_F^2/16] epsilon^4 + o(epsilon^4), m=2.   (1.14)
```

The generic `epsilon^-2` and orthogonal `epsilon^-4` latency laws are not two empirical fits. They are two statistical observability classes.

## 1.6 Why there is no universal imaging–detection conservation law

The mean channel, covariance rank, `tr J`, lower spectral tail, and jet order are distinct invariants. The physical medium can change them in the same direction, opposite directions, or not at all. In particular, complete scrambling gives:

- mean support `{0}`;
- covariance support essentially full;
- covariance scene rank one;
- quotient jet order at most two for a fixed positive-definite `G`;
- a finite information conductance set by visible covariance channels.

No scalar sum of “image information plus detection information” is conserved. The right universality object is the tuple

```text
(support, score rank, quotient-jet order, conductance spectrum). (1.15)
```

I will call this the **support–rank–jet classification** below.

---

# 2. Mandatory move 2 — assumption-inversion sweep

| Unquestioned assumption | Inversion | What breaks | Can this program realize it? |
|---|---|---|---|
| Scene information enters primarily through a mean image | The mean can be exactly blind while covariance is informative | Mean-based resolution and sufficient-statistic intuition | **Yes:** exact band wall and full-scramble DC wall |
| Fluctuations are nuisance | Fluctuations are the transducer | Averaging them away destroys the signal-bearing layer | **Yes:** thin fog and full scramble |
| Temporal averaging always improves a measurement | Retaining bank-resolved fluctuations is necessary | Averaging first projects onto the blind mean channel | **Yes:** independent-bank covariance ledger |
| Resolution is an instrument-only property | Usable aperture is task and statistical-layer dependent | One hardware aperture no longer defines all capabilities | **Yes:** imaging dead, detection live on same spectrum |
| Stronger scattering always worsens sensing | It can destroy imaging while strengthening wall-based specificity and preserving detection | “More disorder is worse” monotonicity | **Yes, conditionally:** scrambling inversion; rates may worsen |
| Detection requires prior identifiability of the changed object | A functional can be detectable when the object is unreconstructable | Reconstruction-first pipelines and identifiability dogma | **Yes:** rank-one `Q(x)` sentinel |
| Zero linear response means impossibility | Positive curvature can rescue detectability | First-order Fisher tests and ordinary score diagnostics | **Yes:** `m=2`, `epsilon^-4` class |
| More photons always help | Beyond visibility saturation, only independent disorder realizations help | Photon-budget optimization | **Yes:** conductance sum rule (1.8) |
| More codes always help | Code information saturates at the illuminated subspace dimension | Pattern-count scaling | **Yes:** measured `M_eff` cap |
| Unknown realizations are unusable | A known law can be sufficient after marginalization | Pathwise calibration requirement | **Yes:** covariance MLE; pathwise route died |
| Formal spectral support implies recoverability | Support can be full while score rank is one | Aperture-as-resolution claims | **Yes:** complete scramble |
| A noise floor is an absolute floor | It may be a higher-order integration law | Fixed-duration declarations | **Yes:** CUSUM and `epsilon^-2/-4` laws |
| A positive quadratic term guarantees global detection | It guarantees local curvature, not global injectivity | Overstated “any finite change” language | **Only on monotone cones or with multiple anchors** |
| Nuisance profiling merely reduces SNR | A nuisance tangent can annihilate the entire scene channel | “Robustness tax” intuition | **Yes:** medium amplitude is exactly confounded with `Q` at zero lag |

The most valuable inversions are not slogans. They each change the order, rank, or quotient of the statistical experiment.

---

# 3. Mandatory move 3 — cross-field transplant sweep

## 3.1 Quantum metrology: Rayleigh’s curse and measurement-dependent visibility order

**The theorem to transplant.** Tsang, Nair, and Lu showed that direct imaging loses Fisher information for sub-Rayleigh separation while a mode-resolved measurement retains it ([DOI 10.1103/PhysRevX.6.031033](https://doi.org/10.1103/PhysRevX.6.031033)).

**Transplant.** Replace “which spatial optical measurement?” by “which statistical layer of the same scalar record?” The mean channel has `m=infinity` for beyond-band perturbations; the covariance channel has `m=1` generically and `m=2` on first-order-orthogonal directions. The analogue of SPADE is not another detector array but a covariance score that changes the statistical contact order.

**New statement sought.** A **statistical Rayleigh curse** occurs when a chosen data reduction pushes the quotient jet to higher order or infinity. A sufficient-statistic change can lower the order without changing hardware.

## 3.2 Stochastic thermodynamics: fluctuation-response inequalities

**The theorem to transplant.** Fluctuation-response and thermodynamic uncertainty relations bound response by fluctuation and information/dissipation; finite-time TURs follow from Cramér–Rao ([DOI 10.1073/pnas.1918386117](https://doi.org/10.1073/pnas.1918386117); [DOI 10.1103/PhysRevLett.125.140602](https://doi.org/10.1103/PhysRevLett.125.140602)).

**Transplant.** Equation (1.5) becomes an optical response bound. The efficient covariance score is the observable that saturates it. Equation (1.8) then says how much response capacity one independent disorder realization can supply.

**New statement sought.** The disorder bank is an information fuel, but its “cost” is not thermodynamic entropy production unless a physical actuation model is added. Do not claim an FDT in the equilibrium sense; claim an exact fluctuation–information response identity with equality conditions.

## 3.3 Disordered systems: Landauer channels and universal conductance

**The theorem to transplant.** Mesoscopic conductance decomposes into transmission eigenchannels, with disorder producing universal fluctuation laws ([DOI 10.1103/PhysRevLett.55.1622](https://doi.org/10.1103/PhysRevLett.55.1622)).

**Transplant.** The covariance visibilities `tau_j=a_j/(1+a_j)` are information-transmission eigenvalues and `I=(1/2)sum tau_j^2`. Complete scrambling can collapse the scene rank while leaving several record-space conductance channels for estimating the surviving scalar functional.

**New statement sought.** Disorder has two distinct ranks: **scene rank** and **information-conductance rank**. Their separation is the reason a rank-one scene functional can still be detected efficiently.

## 3.4 Random matrices: BBP transition for law-light sentinels

**The theorem to transplant.** A covariance spike separates from the sample bulk only above the Baik–Ben Arous–Péché threshold ([DOI 10.1214/009117905000000233](https://doi.org/10.1214/009117905000000233)).

**Transplant.** If the medium law is not known well enough for a matched score, a scene change is a covariance spike or scale perturbation in a high-dimensional bucket covariance. A law-blind largest-eigenvalue sentinel should exhibit a BBP threshold, whereas the matched score can retain power below that eigenvalue-separation threshold.

**New statement sought.** There are two detection phase transitions: physical quotient-jet visibility and finite-sample spectral emergence. Their separation quantifies the price of law ignorance.

## 3.5 Sequential analysis: Lorden/CUSUM optimality

**The theorem to transplant.** Lorden’s minimax change-detection asymptotics make worst-case delay inversely proportional to the per-sample KL information.

**Transplant.** Insert the quotient KL expansion. A first-order-visible change has delay `epsilon^-2`; a curvature-rescued change has delay `epsilon^-4`. The false-alarm threshold supplies only the logarithmic prefactor.

**New statement sought.** Sequential latency exponent equals twice the quotient-jet order. This is a falsifiable scaling law across optics, acoustics, radar, and any stochastic transducer.

## 3.6 Control/estimation duality: nonlinear observability rank conditions

**The theorem to transplant.** Hermann–Krener nonlinear observability uses successive Lie derivatives of outputs to reveal states invisible to the first differential ([DOI 10.1109/TAC.1977.1101601](https://doi.org/10.1109/TAC.1977.1101601)).

**Transplant.** Replace time derivatives of a deterministic output by derivatives of the profiled probability law. The quotient information jet is a stochastic observability jet. `m=2` means the state is invisible to the first statistical differential but visible to curvature.

**New statement sought.** A stochastic observability rank condition built from efficient score jets, with support–rank–jet as its optical specialization.

---

# 4. Mandatory move 4 — named-effect audit

## 4.1 Best candidate: Curvature-Rescued Detection

**Definition.**

> A physical change exhibits **curvature-rescued detection** when its profiled first statistical derivative vanishes but its profiled second derivative is nonzero. Its optimal local error exponent is quartic in change amplitude and its fixed-power observation time scales as the inverse fourth power.

This is general, exact, and not tied to bucket optics. The full-scramble orthogonal direction is the first clean instantiation.

## 4.2 Scrambling Inversion

**Provisional definition.**

> In a **scrambling inversion**, increasing disorder contracts reconstruction rank while preserving or strengthening a lower-dimensional change statistic and its invariance-based specificity.

This is visually powerful, but a true named effect requires the memory-effect interpolation to show a broad monotone family, not only two endpoints.

## 4.3 Blindness-Assisted Specificity

**Definition.**

> An exactly invariant observation channel acts as a physical negative control: activity in an orthogonal statistical layer can be attributed to a class of changes that the invariant channel cannot express.

This is likely Optica-shaped rather than PRL-shaped unless generalized into a theorem on orthogonal experiment quotients.

## 4.4 Fluctuation Information Conductance

**Definition.**

> The fluctuation information conductance of one disorder realization is `G_I=(1/2)sum_j tau_j^2`, where `tau_j` are covariance visibility eigenchannels; it saturates at one-half per visible channel.

This has exact equality conditions and broad reach. It is a strong second theorem, but by itself it is a reinterpretation of Gaussian Fisher algebra.

## 4.5 Disorder Clock

**Definition.**

> After photon and code visibility saturate, independent disorder realizations become the sole extensive resource and detection precision grows only as the square root of their number.

The phenomenon is real but too close to ordinary independent-sample scaling to carry PRL alone.

### Audit verdict

The named-effect candidate with the best combination of novelty, exact scaling, and cross-field reach is **Curvature-Rescued Detection**. “Scrambling inversion” is the best physical framing and may be the title-level effect if the memory-effect interpolation lands.

---

# 5. Mandatory move 5 — ranked innovation portfolio

Scores are `surprise × rigor × reach`, each from 1 to 5. The failed R41 spectral-duality line is not repeated.

## Rank 1 — Quotient Information Jets / Curvature-Rescued Detection

**Score:** `5 × 5 × 5 = 125`

**Mechanism.** Profile the physical law over nuisances and classify a change by the first nonzero derivative of its statistical distance from the nuisance orbit.

**Exact-statement sketch.** If `C_*(epsilon)=c epsilon^(2m)+o(epsilon^(2m))`, optimal fixed-sample and sequential resources scale as `epsilon^(-2m)`. In the full-scramble `Q` channel, `m=1` for `x^T G delta !=0`, `m=2` for `x^T G delta=0` and `delta^T G delta>0`, and `m=infinity` only after exact gauge/nuisance cancellation.

**Decisive test.** Construct generic and exactly `G`-orthogonal beyond-band changes over an amplitude decade; require exact KL slopes `2` and `4`, CUSUM delay slopes `-2` and `-4`, and the predicted crossover when the linear coefficient is small but nonzero. Repeat with a free medium-amplitude nuisance and verify complete annihilation unless a noncollinear lag/wavelength anchor is added.

**Nearest prior art.** Rayleigh-curse quantum metrology, nonregular statistical rates, and nonlinear observability. None appears to formulate nuisance-quotiented statistical jet order as a physical observability class for disordered sensing.

**PRL shape.** A general principle with two universal exponents and one striking optical realization.

## Rank 2 — Fluctuation Information Conductance

**Score:** `5 × 5 × 4 = 100`

**Mechanism.** Diagonalize the medium-to-shot covariance ratio into information eigenchannels; each independent bank carries `1/2 sum tau_j^2` Fisher information for log fluctuation scale.

**Exact-statement sketch.** Equation (1.8), with nuisance-projected generalization. Equality to `r/2` requires every one of `r` visible channels to be medium dominated.

**Decisive test.** Across thin fog, memory-effect interpolation, and full scrambling, compute the exact Fisher independently and require machine-precision agreement with the visibility sum rule; then collapse measured bank requirements against `G_I` rather than photons or nominal contrast.

**Nearest prior art.** Landauer conductance and classical Gaussian covariance Fisher information. The optical disorder-bank interpretation and scene-rank/conductance-rank split appear open.

**PRL shape.** A compact sum rule that explains why disorder states, not photons, are the asymptotic currency.

## Rank 3 — Support–Rank–Jet universality classes

**Score:** `5 × 4 × 5 = 100`

**Mechanism.** Classify a stochastic sensing medium by what frequencies can enter, how many scene directions survive, and at what statistical derivative order each remaining direction becomes observable.

**Exact-statement sketch.** Associate to a model the tuple `(A,r,m(delta))`: covariance support `A`, local scene-score rank `r`, and quotient-jet field `m(delta)`. Thin screens and complete scramblers occupy different endpoints of one memory-effect family.

**Decisive test.** Sweep a controlled transmission ensemble from diagonal/multiplicative to fully random using memory-effect rank `r_ME`; map support, covariance Jacobian rank, and jet-order distribution. Falsify the conjecture if rank does not scale as the predicted `O(r_ME^2)` or if a finite-memory model develops local blind directions inside the positive grain aperture.

**Nearest prior art.** Memory-effect imaging, observability theory, and universality classes in random media. The combined support–rank–jet classification is not standard.

**PRL shape.** A phase diagram of what disorder permits one to image, detect linearly, or detect only through curvature.

## Rank 4 — Minimal-anchor theorem for scene/medium attribution

**Score:** `4 × 5 × 5 = 100`

**Mechanism.** Detection survives nuisance profiling only when the scene-change covariance signature has a component outside the medium-law tangent.

**Exact-statement sketch.** For whitened derivative vectors `h_s` and nuisance matrix `H_eta`, efficient information is

```text
J_s = ||(I-Pi_eta)h_s||^2/2.
```

Zero-lag full-scramble scene energy and medium amplitude are exactly confounded. The minimum number and type of lag, wavelength, polarization, or mean anchors is the minimum augmentation making the projected signature nonzero.

**Decisive test.** Compute derivative signatures for scene `Q`, medium amplitude, and correlation time across lags. Verify zero efficient information with amplitude alone, then identify the smallest lag/wavelength set that gives a robust canonical angle and correct attribution.

**Nearest prior art.** Parameter orthogonality and system identifiability. The exact minimal-anchor result for a covariance-only disordered sentinel appears open.

**PRL shape.** Converts a hidden fatal gauge into a design theorem.

## Rank 5 — Multi-kernel covariance spectroscopy without imaging

**Score:** `5 × 4 × 4 = 80`

**Mechanism.** Probe the same scene with several grain kernels `G_l` generated by wavelength, numerical aperture, or polarization. The vector `Q_l=x^T G_l x` cannot reconstruct the scene in full scramble, but it can encode anomaly scale and spectral class.

**Exact-statement sketch.** A family of kernels separates two anomaly classes whenever their vectors of quadratic energies differ after nuisance quotienting. Local class identifiability is the rank of `[delta^T G_l delta]_l` and cross terms.

**Decisive test.** Use two or three controlled grain sizes and ask whether low-, mid-, and high-frequency anomalies become separable with no image reconstruction. Require class AUC above `0.95` under unknown amplitude and medium drift.

**Nearest prior art.** HBT spectroscopy, multiwavelength speckle methods, and SOFI. The “spectral diagnosis from several rank-one scrambled bucket functionals” conjunction appears open.

**PRL shape.** A capability leap: complete scrambling becomes a spectral analyzer rather than an imager.

## Rank 6 — BBP transition of a law-light covariance sentinel

**Score:** `4 × 5 × 4 = 80`

**Mechanism.** When the medium covariance law is unavailable, detect changes through sample-covariance spectral structure rather than a matched score.

**Exact-statement sketch.** In a high-dimensional limit `M/T -> c`, a scene-induced covariance spike becomes spectrally visible only above a BBP-type threshold; below it, eigenvalue-only detection loses power although a matched or full likelihood-ratio detector may remain informative.

**Decisive test.** Vary `M/T`, spike strength, and law mismatch; map the largest-eigenvalue transition and compare it with the matched-score Chernoff exponent. Kill the line if realistic bucket covariance is too structured for any universal transition.

**Nearest prior art.** BBP and high-dimensional covariance change detection. The disordered optical transducer supplies the physical spike but does not own the mathematics.

**PRL shape.** Gives a sharp law-knowledge phase transition and a reason the semiparametric detector may or may not work.

## Rank 7 — Fluctuation-response optimality of the matched sentinel

**Score:** `4 × 5 × 4 = 80`

**Mechanism.** Treat every implementable statistic as a response observable; the efficient score uniquely saturates the response-to-variance bound.

**Exact-statement sketch.** Equation (1.5), extended to vector alternatives and nuisance quotienting, yields an efficiency ratio equal to the squared correlation with the efficient score.

**Decisive test.** Compare raw covariance energy, eigenvalue, unwhitened projection, and efficient-score detectors. Their measured local efficiencies must equal their squared score correlations; the matched score must sit at one.

**Nearest prior art.** Score-test optimality and fluctuation-response inequalities. The value is compression and physical interpretation, not abstract novelty.

**PRL shape.** Strong supporting theorem, unlikely to carry the Letter alone.

## Rank 8 — Covariance heterodyne rank lift

**Score:** `5 × 3 × 5 = 75`

**Mechanism.** Add a deliberately weak, known unscattered or memory-preserving reference field. Cross-covariance terms between the reference and scrambled field become linear in scene modes, lifting the rank-one `Q(x)` collapse.

**Exact-statement sketch.** With field `E=E_s+alpha E_r`, fourth moments contain terms proportional to `alpha^2 J_s J_r^*`; if the reference coherence spans `r` modes, the covariance scene rank can grow from one to `O(r)` or `O(r^2)`.

**Decisive test.** Extend the exact Kronecker generator with a tunable coherent reference fraction. Measure the smallest reference amplitude at which beyond-band covariance Jacobian rank and Fisher mode coverage rise materially without destroying the wall-based sentinel.

**Nearest prior art.** Heterodyne holography and coherent reference imaging. The surprising target is not ordinary imaging but restoring rank through fourth-order cross terms while retaining a single bucket.

**PRL shape.** A shocking experiment if it works, but it adds hardware and needs a clean theorem.

## Rank 9 — Disorder-echo differential sentinel

**Score:** `5 × 3 × 4 = 60`

**Mechanism.** Use two statistically matched disorder ensembles or two time-reversed code sequences so medium-law drift cancels while scene-induced covariance changes add.

**Exact-statement sketch.** Construct a difference statistic whose nuisance derivative vanishes and whose scene derivative doubles, analogous to an echo. The gain is governed by the canonical angle between the scene and law tangents.

**Decisive test.** Simulate paired diffuser ensembles with common slow drift and independent fast disorder. Require at least a twofold latency improvement and robust cancellation of amplitude and correlation-time changes.

**Nearest prior art.** Coda-wave interferometry and spin-echo logic. The optical one-bucket realization is new-looking but the principle is familiar.

**PRL shape.** Excellent experiment; secondary theory unless the echo cancellation is exact.

---

# 6. Mandatory move 6 — moonshot derivation

## 6.1 The object: quotient Chernoff contact order

Let `P_{epsilon,eta}` be a dominated statistical model. `epsilon` is a scalar physical-change amplitude along scene direction `delta`; `eta` is nuisance. Define

```text
C_*(epsilon)
 = inf_{eta'} max_{0<=s<=1}
   [-log integral p_{epsilon,eta}^s p_{0,eta'}^(1-s)].          (6.1)
```

Assume the infimum is attained locally and

```text
C_*(epsilon)=c_m epsilon^(2m)+o(epsilon^(2m)), c_m>0.           (6.2)
```

### Theorem A — information-jet detection law

For `T` independent observations with equal priors,

```text
lim_{T->infinity} -T^-1 log P_e^*(T,epsilon)=C_*(epsilon).      (6.3)
```

Consequently, in the joint small-change/large-sample regime with target error exponent `h`,

```text
T_req(epsilon)
 = h/[c_m epsilon^(2m)] [1+o(1)].                               (6.4)
```

For quickest detection with false-alarm level `alpha`, if the profiled per-bank KL distance satisfies

```text
D_*(epsilon)=k_m epsilon^(2m)+o(epsilon^(2m)),                  (6.5)
```

then first-order Lorden/CUSUM asymptotics give

```text
E_delay
 = |log alpha|/[k_m epsilon^(2m)] [1+o(1)].                    (6.6)
```

### Why this is invariant

- Reparameterizing `eta` leaves the nuisance orbit unchanged.
- A gain gauge absorbed into `eta` is automatically quotiented.
- No statistic, estimator, or Fisher eigenbasis is selected.
- Physical covariation changes `c_m`; it does not change the definition of `m`.

## 6.2 Gaussian covariance specialization

Let

```text
P_epsilon=N(0,Sigma_epsilon),
Sigma_epsilon=Sigma_0+z(epsilon)H,
Sigma_0>0,                                                     (6.7)
```

and assume first that `H` is not a nuisance direction. Put

```text
A_epsilon=z(epsilon)K,
K=Sigma_0^-1/2 H Sigma_0^-1/2.                                 (6.8)
```

The KL divergence is

```text
D(P_epsilon||P_0)
 = 1/2[tr A_epsilon-log det(I+A_epsilon)]                      (6.9)
 = z^2 tr(K^2)/4-z^3 tr(K^3)/6+O(z^4).                        (6.10)
```

The Bhattacharyya distance is

```text
B(P_epsilon,P_0)
 = 1/2 log det(I+A/2)-1/4 log det(I+A)                         (6.11)
 = z^2 tr(K^2)/16+O(z^3).                                     (6.12)
```

For nearby regular distributions the optimizing Chernoff parameter tends to `1/2`, so (6.12) is the leading Chernoff coefficient.

## 6.3 Full-scramble scene functional

Take

```text
Q(x)=x^T Gx,
G>0,
Sigma(x)=R+Q(x)H.                                               (6.13)
```

Along `x_epsilon=x+epsilon delta`,

```text
z(epsilon)
 = 2c epsilon+d epsilon^2,
c=x^T Gdelta,
d=delta^T Gdelta>0 for delta!=0.                              (6.14)
```

### Case 1 — generic direction

If `c!=0`, then `z=2c epsilon+O(epsilon^2)`. Substituting into (6.10)–(6.12),

```text
D = c^2 ||K||_F^2 epsilon^2+o(epsilon^2),                      (6.15)
C = c^2 ||K||_F^2 epsilon^2/4+o(epsilon^2).                    (6.16)
```

Thus `m=1` and `T_req proportional to epsilon^-2`.

### Case 2 — first-order-orthogonal direction

If

```text
x^T Gdelta=0,
delta!=0,                                      (6.17)
```

then `z=d epsilon^2`. Hence

```text
D = d^2 ||K||_F^2 epsilon^4/4+o(epsilon^4),                    (6.18)
C = d^2 ||K||_F^2 epsilon^4/16+o(epsilon^4).                   (6.19)
```

Thus `m=2` and `T_req proportional to epsilon^-4`.

This is **curvature-rescued detection**: the ordinary score vanishes, but the Hessian of `Q` is `2G>0`, so the statistical law leaves the null at second order in every nonzero first-order-orthogonal direction.

### Equality and boundary conditions

1. **Fast-class equality (`m=1`).** `c!=0`.
2. **Curvature-class equality (`m=2`).** `c=0`, `d>0`, and `H` has nonzero efficient component.
3. **True local blindness.** `d=0` (possible if `G` has a nullspace at coarse grain) or the covariance signature is entirely nuisance.
4. **Global finite-change cancellation.** For signed changes, `Q(x+epsilon delta)=Q(x)` also at

   ```text
   epsilon_*=-2c/d,                                             (6.20)
   ```

   when this value is physically admissible. Therefore the theorem is local. On the monotone cone `x>=0`, `delta>=0`, and entrywise nonnegative `G`, `c,d>=0` and every `epsilon>0` gives `Delta Q>0`; only there is the sign-definite statement global.
5. **Unknown amplitude gauge.** If medium amplitude `a` produces

   ```text
   Sigma=R+a Q(x)H,
   ```

   then `partial_a Sigma=QH` and `partial_Q Sigma=aH` are collinear. Zero-lag covariance has zero efficient information for separating scene `Q` from `a`. A lag, wavelength, polarization, or mean anchor must create a noncollinear signature.

## 6.4 Linear nuisance projection

Let nuisance covariance derivatives form columns `D_eta` in the Gaussian covariance Hilbert space

```text
<A,B>_0=1/2 tr(Sigma_0^-1 A Sigma_0^-1 B).                      (6.21)
```

Let `H_perp=(I-Pi_eta)H`. Replace `K` by

```text
K_perp=Sigma_0^-1/2 H_perp Sigma_0^-1/2.                       (6.22)
```

Equations (6.15)–(6.19) remain valid with `||K_perp||_F`. If `H_perp=0`, every coefficient vanishes and the change is unidentifiable from that ledger. This is the exact profiling condition the current “universal sentinel” language must satisfy.

For nonlinear nuisance families the fully invariant definition remains (6.1); higher nuisance jets may cancel a second-order scene jet even when the linear tangent does not. The numerical experiment must profile the exact model, not only project first derivatives.

## 6.5 Falsifying numerical experiment

Use the exact `SCRAMBLE_EXT` Kronecker generator and no reconstruction.

### Bank A — fixed-law jet slopes

1. Generate positive-definite grain kernel `G` with full spectral support.
2. Construct:
   - `delta_g`: generic beyond-band/DC-free direction with `c!=0`;
   - `delta_o`: exact `G`-orthogonal projection satisfying `x^T Gdelta_o=0`;
   - `delta_n`: a direction concentrated near any small eigenvalues of `G`.
3. Sweep `epsilon` logarithmically over at least one decade and compute:
   - exact Gaussian KL and Chernoff distances;
   - Monte Carlo log-likelihood-ratio error;
   - CUSUM delay at fixed false alarm.
4. Frozen predictions:

```text
slope log KL vs log epsilon:       2 for delta_g, 4 for delta_o;
slope log delay vs log epsilon:   -2 for delta_g, -4 for delta_o;
coefficient ratio: equations (6.15) and (6.18) within 10%;
```

5. For a mixed direction with small `c`, verify crossover near

```text
epsilon_cross approximately 2|c|/d.                            (6.23)
```

### Bank B — nuisance quotient

1. Admit unknown medium amplitude. Prediction: zero-lag efficient information for `Q` collapses to numerical zero.
2. Add lags with known shape but unknown amplitude. If scene and amplitude remain proportional across every lag, the collapse must persist.
3. Add a second grain kernel/wavelength or an independently normalized amplitude anchor. Prediction: information returns exactly when the scene signature leaves the nuisance span.

### Bank C — global cancellation and cone

1. Choose signed `delta` with admissible `epsilon_*=-2c/d`; verify AUC returns to `0.5` at the iso-`Q` crossing despite local visibility.
2. Choose additive nonnegative defects; verify `Delta Q>0` for every tested amplitude.

### Kill conditions

Kill the PRL moonshot if any of the following occurs:

- exact profiled KL does not exhibit integer contact orders predicted by the quotient-jet calculation;
- the `epsilon^-4` class disappears after finite-sample correction rather than merely becoming expensive;
- nuisance profiling leaves unexplained information in the exactly collinear amplitude model;
- memory-effect interpolation produces a broad set of directions with all quotient jets zero inside the claimed grain aperture;
- the same phenomenon is found already named and proved in a directly equivalent optical/statistical experiment.

## 6.6 Why this is PRL-shaped

A general reader need not care about DMDs or ghost imaging to understand the principle:

> A physical state can be invisible to linear response yet observable through statistical curvature, and the order of that curvature fixes the universal resource exponent.

The optical system then supplies a severe realization: complete scrambling removes the image, the mean, and almost every scene coordinate, yet leaves a locally complete change detector whose hard directions obey a different asymptotic law. That is a physics result about observability in disorder, not only an imaging algorithm.

---

# 7. Suggested PRL seed if the moonshot survives

## Candidate titles

1. **Curvature-rescued detection through complete optical scrambling**
2. **Statistical observability beyond linear response in disordered light**
3. **Strong disorder erases images but not changes**
4. **Information jets of a scrambled optical measurement**

## One-sentence thesis

> **Complete optical scrambling collapses scene reconstruction to rank one, yet every local scene-change direction remains statistically observable at first or second order, producing universal `epsilon^-2` and `epsilon^-4` detection-time classes after nuisance quotienting.**

## The one figure a PRL referee should remember

A disorder axis from thin fog to full scrambling, with three aligned rows:

```text
mean support:        B_p ------------------------------> {0}
reconstruction rank: many modes -----------------------> 1
detection jet:       m=1/2 directions persist --------> m<=2 locally
```

Under it, one log–log panel shows the two latency slopes `-2` and `-4`, and one nuisance panel shows the amplitude-confounded ledger at zero followed by restored information after the minimal anchor. The visual inversion is immediate: **rank dies; jet visibility survives.**

---

# Final innovation verdict

The broadest correct object is not “detection beats estimation,” not a fixed-trace spectral inequality, and not an information-conservation slogan. It is:

> **Statistical observability has an order. Support says where a perturbation may enter; score rank says what can be reconstructed; the quotient information jet says how an otherwise invisible change first becomes distinguishable.**

The full-scramble endpoint is the decisive test bed because it separates these three objects maximally: total mean blindness, full covariance support, rank-one reconstruction, and first-/second-order local detectability. If the slope and nuisance tests above land, **Curvature-Rescued Detection** is the first candidate in this program that I would deliberately shape for PRL rather than merely hope PRL accepts.