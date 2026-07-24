# R45 — Deep mine first: symmetry walls, jet flow, and the certifiability hierarchy; then Paper 2

**Reference brief:** [`docs/ROUND63_GPT_ROUND45_QUESTION.md`](https://github.com/ccyyyYyzz/GI_a2/blob/8701970/docs/ROUND63_GPT_ROUND45_QUESTION.md) at [`8701970`](https://github.com/ccyyyYyzz/GI_a2/commit/87019708c2abe57ddb58d46129694e8f36fd4442).  
**Evidence index:** [`results/round63_next/R44_KILL_TESTS/INDEX.md`](https://github.com/ccyyyYyzz/GI_a2/blob/8701970/results/round63_next/R44_KILL_TESTS/INDEX.md).  
**Order lock:** Part B below is the derivation-first deep mine. The Paper-2 ruling appears only afterward in Part C.

---

# PART B — THE DEEP MINE

## B0. What the full battery has actually exposed

The 17-test record supports a deeper hierarchy than “mean channel versus covariance channel.” There are three logically distinct questions:

1. **Wall:** does a symmetry or conservation law make the full record law invariant to the perturbation? If yes, no amount of data helps.
2. **Jet:** after quotienting nuisances and broken symmetries, at what order does the law first move? This sets the physical sample-complexity exponent.
3. **Certificate:** even when the law moves, is the direction known, merely known to exist somewhere in a large tangent space, or required to be estimated? These cost parametrically different numbers of banks.

The central object produced by this dig is therefore the

> **Symmetry–Jet–Certifiability hierarchy:** symmetry fixes whether a direction exists in the statistical quotient; the first nonzero quotient jet fixes its amplitude exponent; the dimension of the unknown efficient tangent fixes the finite-sample certification penalty.

The strongest new scaling law is

```text
known direction        : T · lambda             ~ 1
blind existence test   : T · lambda / sqrt(p)   ~ 1
full direction estimate: T · lambda / p         ~ 1,
```

where `lambda = (1/2)||A_perp||_F^2` is the per-bank efficient covariance noncentrality after whitening and nuisance projection, and `p` is the dimension of the admitted efficient covariance tangent. This is derived in full in §B11.

It turns the “possible third ledger” into a precise statement:

> The estimator ledger is not another physical moment channel. It is the adaptation price paid when the informative tangent is unknown.

For KT6, `M=16` gives `p=M(M+1)/2=136`. Its signal was `0.0054` of the raw `p/T` Wishart floor. Even the optimal blind threshold is lower than the plug-in floor by only `sqrt(p)=11.66`, so the signal remains

```text
0.0054 sqrt(136) = 0.063
```

of the minimax blind-certification scale—about **16× too small**. The ranging kill is therefore fundamental at that ledger, not merely a defect of the plug-in estimator.

---

# B1. Ore vein 1 — a statistical Noether theorem for blind channels

## B1.1 Orbit–null theorem

Let `x` lie on a smooth parameter manifold `X`, and let the complete observed record have a dominated law `P_x` with density `p(y|x)`. Let a Lie group `G` act smoothly on `X`, with infinitesimal generator `X_xi(x)` for Lie-algebra element `xi`.

### Theorem 1 — Statistical Noether wall

If

```text
P_{g·x} = P_x    for every g in G,                              (1.1)
```

then every orbit tangent is an exact score-null direction:

```text
D_x log p(Y|x)[X_xi(x)] = 0          almost surely,             (1.2)
I_x X_xi(x) = 0.                                               (1.3)
```

Every `f`-divergence, Chernoff distance, and quotient information jet between two points on the same orbit is zero. No estimator, moment order, exposure count, or reconstruction prior can distinguish points on that orbit from this record.

### Proof

Differentiate `log p(y|exp(t xi)·x)=log p(y|x)` at `t=0` to obtain (1.2). Taking the score Gramian gives (1.3). Equality of the complete laws gives zero divergence at every finite displacement along the orbit. `□`

### Local converse

Suppose a smooth, constant-rank distribution `D_x` is integrable and

```text
D_x log p(y|x)[v] = 0
```

for every `v in D_x`, almost every `y`, and every `x` in a connected neighborhood. Then `p(y|x)` is constant along each integral leaf of `D`. The leaves are local blindness orbits. Thus, under ordinary smoothness and integrability, continuous local walls are exactly local gauge symmetries of the statistical experiment.

This is the honest Noether analogy. A symmetry of the **underlying wave equation** does not automatically create a wall. It creates a wall only when the **full measurement law** is invariant under the induced parameter action.

## B1.2 The two confirmed walls are two group actions

### Support wall — additive translation symmetry

For a band-limited mean operator `A`, every `h in ker A` generates

```text
x -> x + h,
P_{x+h}^{mean} = P_x^{mean}.                                   (1.4)
```

The group is the additive group of the unmeasured subspace. A hard field pupil restored this wall to machine precision (`7e-16`); the finite detector window then supplied a different breaking operator.

### Energy wall — unitary orbit symmetry

A lossless full-aperture bucket measures only the Hilbert-space norm:

```text
B(E) = ||E||_2^2.
```

For any unitary spatial transformation `U`,

```text
B(UE)=B(E).                                                     (1.5)
```

The phase-only wall confirmed by IT1 is the diagonal subgroup `U=diag(e^{i phi})`; the measured null `2.18e-16` is one orbit of the much larger unitary-mode wall. NA clipping replaces the identity collector by a projector `P_NA`, and

```text
||P_NA U E||^2
```

is no longer invariant. The wall then breaks at first order, exactly as IT1 measured.

## B1.3 Wall periodic table for an optical chain

| Wall class | Invariance / group | Exact null | Breaking operator | Restoration or readout hardware |
|---|---|---|---|---|
| **Support wall** | additive translations by `ker A` | scene modes outside the physical transfer range | larger pupil, aliasing, finite window, nonlinear moment | hard pupil + guarded sampling; measured-null codes |
| **Unitary-energy wall** | `E -> U E`, `U†U=I`, under total-energy collection | every lossless mode mixing, including phase-only changes | NA clipping, segmented or mode-selective collection | full aperture restores; mode sorter/analyzer intentionally breaks |
| **Global phase wall** | `E -> e^{i phi}E` under intensity-only detection | absolute optical phase | coherent reference | homodyne/heterodyne arm |
| **Parity/exchange wall** | `x(r)->x(-r)` when all codes, propagation, and collection are parity symmetric | mirror pair; at a parity-even baseline, odd tangent modes | asymmetric code, pupil, or off-axis reference | parity-diverse code pair |
| **Temporal-origin wall** | stationary law invariant under `t->t+t0` | absolute phase/origin of a stationary fluctuation process | synchronized clock or nonstationary drive | pilot timing/reference marker |
| **Polarization-unitary wall** | local/global `SU(2)` rotation under a polarization-insensitive total-intensity bucket | polarization rotation at fixed total intensity | analyzer, birefringent reference | Stokes-resolved bucket pair |
| **Fully scrambled orbit wall** | orthogonal/unitary transformations preserving the exposed scalar `Q(x)` | the full iso-`Q` orbit | second noncollinear kernel, wavelength, memory-effect mode | multi-kernel anchor |

Two genuinely new walls beyond the R44 pair follow immediately:

1. **Mode-unitary wall:** the full-aperture energy bucket is blind not merely to phase screens but to *any* unitary spatial-mode transformation. Phase-only blindness is one subgroup.
2. **Parity wall:** a reflection-symmetric chain cannot distinguish mirror-related scenes; at a symmetric baseline its first-order odd subspace is score-null. It is broken by one asymmetric code or pupil, not by more symmetric exposures.

The polarization row is the internal-state analogue and is structurally identical to a decoherence-free/noiseless subsystem: the detector observes the invariant norm, not the internal orientation.

## B1.4 Completeness and non-completeness

There is no finite universal list of walls independent of apparatus. The complete local list for a declared experiment is the integrable kernel distribution of its score map. The theorem supplies the classification rule:

```text
exact continuous wall = orbit/leaf on which the complete law is constant.
```

Discrete degeneracies—parity, conjugation, source–detector exchange—may produce globally indistinguishable pairs without a nonzero infinitesimal tangent at a generic point. They are walls, but not Lie-generated local nulls.

Closest conceptual relatives are gauge nonidentifiability, control-theoretic unobservable subspaces, and decoherence-free subspaces. In the latter, error generators annihilate a protected subspace; see Lidar, Chuang, and Whaley, [Phys. Rev. Lett. 81, 2594 (1998)](https://doi.org/10.1103/PhysRevLett.81.2594), and Zanardi and Rasetti, [Phys. Rev. Lett. 79, 3306 (1997)](https://doi.org/10.1103/PhysRevLett.79.3306). The optical novelty candidate is not the abstract kernel idea; it is **hardware-calibrated, task-selective blindness with a measured dimension–leak capacity curve**.

---

# B2. Ore vein 2 — the geometric jet flow

KT2 established the exact physical normal form

```text
C_phys(epsilon) = A epsilon^2 + B epsilon^4,                    (2.1)
```

for a curvature-rescued covariance direction contaminated by a first-order mean leak. IT1 established a second geometry coefficient that grows over six orders with propagation while retaining first-order amplitude slope.

The RG analogy is valid—but only for the asymptotic normal form, not as a universal flow of the entire wave system.

## B2.1 General two-monomial flow

Let the profiled Chernoff distance have two leading operators:

```text
C(epsilon,g)
 = A g^(2 a1) epsilon^(2 m1)
 + B g^(2 a2) epsilon^(2 m2),
0 < m1 < m2.                                                    (2.2)
```

Define

```text
y = (B/A) g^[2(a2-a1)] epsilon^[2(m2-m1)].                     (2.3)
```

The effective amplitude and geometry orders are

```text
m_eff = (1/2) d log C / d log epsilon
      = (m1 + m2 y)/(1+y),                                     (2.4)

a_eff = (1/2) d log C / d log g
      = (a1 + a2 y)/(1+y).                                     (2.5)
```

The flow equations are exact:

```text
d m_eff / d log epsilon
 = 2 (m_eff-m1)(m2-m_eff),                                     (2.6)

d m_eff / d log g
 = 2[(a2-a1)/(m2-m1)](m_eff-m1)(m2-m_eff).                     (2.7)
```

The fixed points are `m1` and `m2`. The lower-order wall-breaking perturbation is relevant in the small-amplitude limit; setting its coefficient to zero by an optical wall, certified code, or nuisance projection puts the experiment exactly on the higher-order fixed manifold.

## B2.2 R44 specialization

For a pre-existing mean leak and near-contact covariance conversion,

```text
m1=1, a1=0,
m2=2, a2=2,
```

because covariance Chernoff information is `O(z2^4 epsilon^4)`. Therefore

```text
dm/dlog epsilon = 2(m-1)(2-m),                                 (2.8)
dm/dlog z2      = 4(m-1)(2-m),                                 (2.9)
```

and the crossover invariant is

```text
y proportional to epsilon^2 z2^4,
epsilon z2^2 = constant.                                       (2.10)
```

Thus

```text
z2,* proportional to epsilon^(-1/2),                           (2.11)
```

which is the geometric version of KT2’s `epsilon_cross=sqrt(A/B)` law.

### Interpretation

- **Unprojected:** `A>0`; as `epsilon->0`, `m_eff->1`, matching the measured slopes `2.12–2.20` in the divergence.
- **Projected or pupil-hardened:** `A=0`; the system sits at `m=2`, matching slope `4.00` in all four cells.
- **Increasing conversion distance:** increases `y` and moves the finite-amplitude experiment from leak-dominated to curvature-dominated behavior.

This is a genuine logistic jet flow, with universal crossover exponents inside the two-term normal form. It is not a full renormalization group for diffraction: additional channels, saturation, finite apertures, and higher jets add operators and can change the beta function.

## B2.3 Named effect

> **Jet transmutation:** a lower-order symmetry-breaking channel changes the strict local observability class; calibrating and projecting that tangent restores the higher-order class without changing the underlying scene perturbation.

KT2 is an exact physical realization. The closest classical analogue is the change of observability rank under output injection in nonlinear control, following Hermann and Krener’s differential-observability framework ([IEEE TAC 22, 728 (1977)](https://doi.org/10.1109/TAC.1977.1101601)).

---

# B3. Ore vein 3 — the estimator-layer impossibility

The observed Wishart floor should not be promoted in its raw plug-in form. The plug-in `p/T` floor is the **full estimation** scale. An optimal direction-agnostic existence test can improve it to `sqrt(p)/T`, but not to the oracle `1/T` scale.

This gives three finite-sample tiers for a whitened covariance displacement `A_perp`:

```text
per-bank efficient signal lambda = (1/2)||A_perp||_F^2.
```

| Task | Direction knowledge | Required scale |
|---|---|---|
| matched detection | direction known | `T lambda ~ 1` |
| blind existence certification | direction unknown in `p`-dimensional tangent | `T lambda ~ sqrt(p)` |
| estimate/localize the direction | all `p` coordinates required | `T lambda ~ p` |

The second line is the actual minimax “third ledger.” It is derived, with lower and upper bounds, in §B11.

## B3.1 KT6 is below even the optimal blind tier

For `M=16`, `p=136`. KT6 measured

```text
lambda_signal / (p/T) = 0.0054.
```

Relative to the optimal blind threshold `sqrt(p)/T`, this is

```text
0.0054 sqrt(136) = 0.063.
```

The signal is therefore approximately `15.9×` too small even for the optimal unknown-direction quadratic test. No better plug-in covariance, shrinkage inverse, or matched-filter regularizer can turn the near-contact range observable into a certified blind bench result at the frozen `(M,T)`.

The correct conclusion is stronger than KT6’s original wording:

> Channel-ratio ranging is not merely plug-in-estimator-limited; at the tested bank count it lies below the minimax blind covariance-certification boundary.

## B3.2 What can and cannot be certified

Let `H` be the nuisance-profiled covariance tangent of dimension `p`.

- A specified scalar functional `⟨A,B⟩` is certifiable at the oracle rate if `B` is fixed independently of the data.
- Existence of an unspecified perturbation somewhere in `H` pays the `sqrt(p)` adaptation penalty.
- Recovering the perturbation direction or a full covariance difference pays `p`.
- If the physical direction lies in the nuisance tangent, `A_perp=0`, the problem returns to the wall ledger and no finite bank count helps.

This is the clean separation between **physical absence**, **channel weakness**, and **finite-sample adaptation**.

Closest statistical relatives are Gaussian covariance identity testing and high-dimensional mean testing. The nonrobust covariance-testing problem has sample complexity linear in vector dimension for fixed Frobenius separation; see the context summarized by Diakonikolas and Kane, [COLT 2021](https://proceedings.mlr.press/v134/diakonikolas21a.html). The new optical object is the coupling of that minimax threshold to a quotient jet and a calibrated physical wall.

---

# B4. Ore vein 4 — blindness certificates as a metrological primitive

## B4.1 Blindness-capacity theorem

Let a calibrated linear leak operator be

```text
L : C -> Y,
```

mapping a code vector to the forbidden mean-channel response. Let the code norm already be whitened by the desired in-band response metric. Write singular values in ascending order

```text
sigma_1 <= ... <= sigma_n.
```

Define the best worst-case leak achievable by a `d`-dimensional code subspace:

```text
beta_d(L)
 = min_{dim S=d} max_{c in S, ||c||=1} ||Lc||.                  (4.1)
```

### Theorem 2 — Blindness capacity

```text
beta_d(L) = sigma_d(L).                                        (4.2)
```

Equality is attained by the span of the bottom `d` right singular vectors. No other `d`-dimensional code family can guarantee a smaller worst-case leak.

For a generalized relative-leak certificate with desired-response metric `B†B`, replace singular values by the generalized eigenvalues of

```text
L†L v = lambda B†B v;                                          (4.3)
beta_d = sqrt(lambda_d).
```

This is precisely the mathematical object measured by IT4/IT4b. The single-plane conjugate result—80 dimensions all below `1e-5`—is a measured blindness-capacity point, not an arbitrary code-design success.

### Proof

Equation (4.2) is the min–max characterization of singular values. The bottom `d` right-singular subspace attains `sigma_d`; any `d`-dimensional subspace intersects the orthogonal complement of the bottom `d-1` singular vectors and therefore contains a unit vector with leak at least `sigma_d`. `□`

## B4.2 Robust calibration certificate

Suppose a calibration produces `Lhat` and an independently justified operator error bound

```text
||L-Lhat||_op <= delta.
```

If `S_d` is the bottom-`d` right-singular subspace of `Lhat`, then

```text
max_{c in S_d, ||c||=1} ||Lc||
 <= sigma_d(Lhat)+delta.                                       (4.4)
```

This converts one calibration run into a metrological **cannot-see guarantee**. Equation (4.4), with a fresh-data estimate of `delta`, is the missing ingredient that separates a numerical null space from a certified instrument property.

## B4.3 Noncomposability theorem

For several conditions or leak routes `L_1,...,L_k`, the simultaneous root-sum-square certificate is governed by the stacked operator

```text
L_stack = [w_1 L_1; ...; w_k L_k],
L_stack† L_stack = sum_i w_i^2 L_i†L_i.                        (4.5)
```

Therefore

```text
beta_d^2 = lambda_d^up(sum_i w_i^2 L_i†L_i).                   (4.6)
```

The certificate factors are **not multiplicative**. Multiplication occurs only for true serial contractions whose singular-vector chains align. Independent hardware modifications generally change the leak operator and its metric; they must be recalibrated jointly.

IT5 is the concrete counterexample:

```text
SVD-only single-z1 capacity at 1e-4 : 80 dimensions
DPSS + recalibrated SVD capacity    : 24 dimensions.
```

The DPSS route suppresses the targeted annular pedestal to `2.4e-4`, but it rotates/reshapes the operator and trades code-space capacity. The routes are alternatives, not factors.

## B4.4 Cross-field fence

- **Decoherence-free subspaces:** states annihilated by error generators; the closest structural cousin.
- **NMR/dynamical decoupling:** pulse averaging suppresses unwanted generators; Viola, Knill, and Lloyd, [PRL 82, 2417 (1999)](https://doi.org/10.1103/PhysRevLett.82.2417).
- **Null-space control:** controlled-invariant unobservable subspaces.
- **Differential privacy:** also certifies limited distinguishability, but probabilistically over neighboring datasets; the optical certificate is a deterministic hardware/operator guarantee and should not be called privacy.

The defensible new statement is:

> A physical sensor can publish a calibrated dimension–leak frontier describing the largest code space on which it is guaranteed not to respond to a declared nuisance.

---

# B5. Ore vein 5 — the unsaturation crossover and coherence-rank spectroscopy

KT4 showed that a single thin phase screen stays highly multimode (`R≈107–122`) and never approaches the `R=1` complete-scrambling endpoint. T5b showed an unsaturated code-count law close to `lambda proportional to M^1.8`.

The correct crossover variable is the number of measured code pairs versus the number of latent quadratic coherence modes.

Let

```text
p = M(M+1)/2
```

be the number of symmetric code-pair coordinates, and let

```text
q ≈ R^2
```

be the number of latent quadratic coherence products. For an isotropic random-feature sketch `X in R^{p×q}` with rows `N(0,I/q)`, the information Gramian is `S=X†X`. Exactly,

```text
E tr S   = p,
E tr S^2 = p + p(p+1)/q.                                      (5.1)
```

The participation-rank approximation therefore gives

```text
r_eff(M;R)
 ≈ [E tr S]^2 / E tr S^2
 = p q/(p+q+1).                                                 (5.2)
```

Consequences:

```text
p << q : r_eff ≈ p ≈ M^2/2      (unsaturated pair growth),
p >> q : r_eff ≈ q ≈ R^2        (coherence-limited saturation). (5.3)
```

The crossover occurs at

```text
p ≈ q  ->  M_* ≈ sqrt(2) R.                                    (5.4)
```

The local log-slope is approximately

```text
alpha(M)=d log r_eff/d log M ≈ 2q/(p+q).                       (5.5)
```

Thus `alpha≈1.8` means `p/q≈0.11`: the present system is still well before saturation, which is exactly why KT4’s `R^2` bound was nonbinding.

### What is theorem and what is model

- `rank <= min(p,q,N)` is algebraic and exact under the declared factorization.
- Equation (5.2) is a random-isotropic effective-rank law, not universal wave physics.
- The testable new proposal is **coherence-rank spectroscopy**: infer `q` or `R` from the code-count crossover without reaching `R=1` geometrically.

### Decisive test

Use shared wave pools and `M in {8,16,32,64,96,128,192}`. Fit the full curve (5.2), not one exponent. Independently estimate coherence participation rank from the field. Require:

```text
fitted sqrt(q) within 25% of the independent R proxy,
pre-crossover slope within 0.2 of 2,
post-crossover slope < 0.8 when p/q > 3.
```

Failure would kill the spectroscopy claim while leaving the algebraic rank bound intact.

---

# B6. Ore vein 6 — the task-dependent étendue sweet spot

The conductance identity

```text
I = (1/2) sum_j [r_j/(1+r_j)]^2                       (6.1)
```

alone does not determine an aperture optimum, because opening the aperture changes both the number of visible channels and the coupling of the target change to them.

Let `N=N_eff` be effective speckle/grain count, total photons per bank be `n`, the number of useful record channels scale as

```text
K(N) proportional to N^beta,
```

and the fractional task coupling scale as

```text
chi(N) proportional to N^(-gamma).
```

Under an equal-channel approximation `r=n/N`, the task information scales as

```text
I_task(N;n)
 = C N^(beta-2gamma) [n/(n+N)]^2.                              (6.2)
```

Write

```text
nu = beta-2gamma.
```

Then

```text
d log I_task/d log N
 = nu - 2N/(n+N).                                               (6.3)
```

An interior optimum exists only if `0<nu<2`, at

```text
N_* = [nu/(2-nu)] n.                                           (6.4)
```

### Important cases

1. **Global covariance-amplitude metrology:** `beta=1`, `gamma=0`, so `nu=1` and
   
   ```text
   N_*=n.
   ```
   Information is maximized at covariance/shot visibility matching `r=1`.
2. **Localized anomaly diluted by bucket averaging:** `gamma` can be large enough that `nu<=0`; then there is no interior sweet spot and the smallest étendue compatible with spatial coverage is best.
3. **Saturated high-photon regime:** if `N<<n`, `I_task proportional to N^nu`; KT5’s `T_det proportional to N_eff` implies a negative effective `nu` for its particular change coupling, not a universal amplitude-metrology law.

This resolves KT5’s partial result. The saturated-regime `25%` collapse is real; the global `64%` failure occurs because `beta` and `gamma` are task- and aperture-dependent. There is no universal scalar `n/N_eff` collapse without measuring the coupling exponent.

### Decisive test

For each aperture, separately estimate:

- `K(N)` from the record-space participation rank;
- `chi(N)` from the normalized derivative of the target covariance;
- then predict `N_*` from (6.4) before running the photon sweep.

The law survives only if the observed optimum/monotonic branch agrees without refitting `beta` and `gamma` to the latency curve itself.

---

# B7. Assumption-inversion sweep

| Standing assumption | Inversion | What breaks | Evidence/status |
|---|---|---|---|
| A wall is merely a limitation | A wall can be an engineered specificity resource | “more sensitivity is always better” | two exact walls; certified codes |
| Certification says what a sensor sees | Certification can guarantee what it cannot see | ordinary calibration vocabulary | IT4/IT4b capacity curve |
| Every useful code should illuminate the target channel strongly | Codes may be selected to lie in a nuisance-blind subspace | generic random-code doctrine | 80 conjugate-plane dimensions ≤`1e-5` |
| Wall-restoration mechanisms multiply | They generally rotate the operator and are noncomposable | independent-factor budgeting | DPSS + SVD: `80 -> 24` dimensions |
| Oracle physical information implies a bench observable | Unknown-direction adaptation may make it uncertifiable | oracle-wave-twin optimism | KT6 minimax deficit |
| More detector segments improve covariance SNR | Segmentation can destroy the very conservation wall supplying specificity | detector-array monotonicity | IT6 kill |
| A second arm repairs false alarms | The wall can starve that arm of target power | coincidence-veto intuition | IT7: 76% power loss |
| Banks alone are the currency | Tangent dimension is an additional currency | scalar sample-complexity laws | `1:sqrt(p):p` hierarchy |
| More codes always help | Pair growth saturates at latent coherence dimension | M-only scaling | `M^1.8` unsaturated; crossover open |
| More aperture/étendue always helps | It trades more modes against weaker per-mode contrast/coupling | throughput intuition | KT5 two-regime law |
| Exact statistical wall means exact bench wall | finite aperture and sampling change the function class | model-to-bench transfer | KT1 layered wall |
| Symmetry breaking is always harmful | controlled breaking selects the jet order and can create sensitivity | protection-only framing | IT1 clipping; KT2 projection |
| Higher-order observability is purely intrinsic to the scene change | lower-order hardware leakage is a relevant operator | static jet classification | jet transmutation |
| The bucket is just the detector element | collection aperture defines the conserved quantity and the wall | component-centric design | full vs clipped bucket |
| A nuisance projection only costs power | it can reduce tangent dimension and improve certifiability while retaining 86–99% target information | one-sided “profiling tax” intuition | KT3 |

---

# B8. Cross-field transplant sweep

## B8.1 Decoherence-free subspaces

**Transplant:** replace environmental error generators by calibrated optical leak operators. A code subspace is “blindness-free” with respect to the nuisance when those operators annihilate it.

**New optical statement:** the maximum protected dimension at leak tolerance `tau` is the count of generalized singular values below `tau`, with an out-of-sample bound (4.4). Unlike an abstract DFS, the generators are measured from the actual instrument and the protected subspace is allowed to retain a different sensing channel.

## B8.2 Noether/gauge theory

**Transplant:** a measurement symmetry, not merely a dynamical symmetry, produces a score-null orbit. The local converse identifies every integrable continuous score-null distribution with a local gauge leaf.

**New optical statement:** hard pupils and full-aperture collectors instantiate two inequivalent gauge groups—additive support translations and unitary-energy orbits—and their breaking operators carry different jet orders.

## B8.3 Renormalization-group language

**Transplant:** lower-order wall-breaking terms are relevant perturbations. The exact beta function (2.6) flows between integer contact-order fixed points.

**Boundary:** this is a normal-form RG, not a claim of universality for the full diffraction field. It survives only while two monomials dominate.

## B8.4 Le Cam minimax theory

**Transplant:** quotient the nuisance tangent, pass to the local Gaussian score experiment, and compare known-direction, unknown-direction, and full-estimation risks.

**New optical statement:** the adaptation dimension `p` multiplies the bank requirement by `sqrt(p)` or `p` without changing the physical jet exponent.

## B8.5 NMR and dynamical decoupling

**Transplant:** complementary exposures/PWM are toggling-frame controls; whole-period integration cancels the nuisance generator. The ordered-dither failure is the optical version of an incomplete decoupling cycle.

**New optical statement:** wall restoration is governed by the average leak operator, while noncommuting windows/codes require joint—not factorized—certification.

## B8.6 Control observability

**Transplant:** the measured mean leak operator gives a Kalman-like unobservable subspace; higher quotient jets play the role of higher Lie derivatives in nonlinear observability.

**New optical statement:** a direction may be unobservable to first order yet observable at second order after nuisance quotienting; the finite-sample adaptation tier is separate from observability rank.

## B8.7 Random matrices / BBP

**Transplant:** an unknown covariance direction becomes statistically visible only when its score norm emerges above a high-dimensional noise shell. The `sqrt(p)` blind threshold is the dense analogue of a spectral-emergence transition.

**New optical statement:** the KT6 oracle signal can be nonzero at every geometry yet remain below the blind-certification transition.

---

# B9. Named-effect audit

1. **Statistical Noether Wall**  
   A continuous symmetry of the complete record law generates an exact Fisher-null orbit.
2. **Blindness Capacity**  
   The `d`-dimensional cannot-see guarantee of a calibrated sensor equals the `d`th smallest generalized singular value of its leak operator.
3. **Jet Transmutation**  
   A lower-order symmetry-breaking channel changes the local observability class; projection or wall restoration returns the higher-order class.
4. **Certifiability Gap**  
   The same physical covariance displacement costs `1`, `sqrt(p)`, or `p` units of `T lambda` depending on whether its direction is known, only its existence is sought, or its direction is estimated.
5. **Noncomposable Blindness**  
   Independent suppression routes combine through the singular spectrum of a stacked operator, not by multiplying leak factors.
6. **Coherence-Rank Unsaturation**  
   Code-pair information grows approximately as `M^2` until the pair count reaches the latent quadratic coherence dimension.
7. **Étendue Matching**  
   For global covariance-amplitude metrology with equal modes, information is maximized near `N_eff=n`; other tasks shift or remove the optimum through their coupling exponent.

The best general-physics name is **Certifiability Gap**. The best optics-facing name is **Blindness Capacity**.

---

# B10. Ranked innovations

Scores are `surprise × rigor × reach`, each on a 1–5 scale.

## Rank 1 — Symmetry–Jet–Certifiability hierarchy

**Score:** `5 × 5 × 5 = 125`

**Mechanism:** quotient physical symmetries, identify the first nonzero jet, then pay the adaptation dimension of the remaining tangent.

**Exact statement:** §B11 gives the `1:sqrt(p):p` theorem and its jet/geometry specialization.

**Kill test:** the THREE-LEDGER SCALING GATE in Part C. Known-direction, blind-existence, and full-estimation contours must collapse on `T lambda`, `T lambda/sqrt(p)`, and `T lambda/p` across `M,T,p`.

**Nearest prior art:** Le Cam LAN, Gaussian mean/covariance testing, nonlinear observability. None alone couples the minimax hierarchy to a physical symmetry wall and quotient jet.

**Why broad physicists care:** a nonzero physical signal is not the same thing as a certifiable signal; the gap has a universal dimension law.

## Rank 2 — Blindness Capacity theorem

**Score:** `5 × 5 × 4 = 100`

**Mechanism:** calibrate the nuisance leak operator and use its bottom generalized singular subspace as the largest guaranteed blind code family.

**Exact statement:** `beta_d=sigma_d`, with calibration robustness `sigma_d(Lhat)+delta` and stacked-operator noncomposability.

**Kill test:** fresh calibration/test split with hardware-law perturbations. Require at least `95%` simultaneous coverage of the predicted leak envelope and no code exceeding `sigma_d(Lhat)+delta`.

**Nearest prior art:** decoherence-free subspaces, null-space control, robust SVD. The optical novelty is a publishable, measured cannot-see specification that preserves a second sensing channel.

## Rank 3 — Statistical Noether wall classification

**Score:** `5 × 4 × 5 = 100`

**Mechanism:** every continuous invariance of the full measurement law generates a score-null orbit; integrable score-null distributions are locally gauge leaves.

**Kill test:** implement parity and mode-unitary walls in the twin, break each with one controlled operator, and verify predicted contact order and restoration.

**Nearest prior art:** gauge nonidentifiability, Noether theory, control unobservability. The new value is the optical wall periodic table tied to measured breaking jets.

## Rank 4 — Geometric jet-flow normal form

**Score:** `5 × 5 × 4 = 100`

**Mechanism:** competing wall-breaking monomials generate an exact logistic beta function for effective contact order.

**Exact statement:** Eqs. (2.6)–(2.10).

**Kill test:** dense `(epsilon,z2)` grid after mean-leak projection and deliberate leak injection; all curves must collapse against `y proportional to epsilon^2 z2^4` with no refit of exponents.

**Nearest prior art:** crossover scaling and normal-form RG. The new physical object is the flow of statistical observability order.

## Rank 5 — Coherence-rank spectroscopy from code unsaturation

**Score:** `5 × 4 × 4 = 80`

**Mechanism:** code-pair features fill a latent `R^2` quadratic space, giving `r_eff≈pq/(p+q)` and a measurable slope crossover near `M≈sqrt(2)R`.

**Kill test:** shared-pool `M` sweep through and beyond the fitted crossover; independent coherence rank must agree within `25%`.

**Nearest prior art:** random-feature rank, Wishart participation ratios, mesoscopic channel counting. The novelty is measuring scattering coherence rank from detector-code scaling without transmission-matrix access.

## Rank 6 — Task-dependent étendue matching

**Score:** `4 × 5 × 4 = 80`

**Mechanism:** aperture changes channel count and task coupling; the optimum follows `N_*=[nu/(2-nu)]n`, not a universal `n/N` curve.

**Kill test:** estimate `beta` and `gamma` independently, predict the aperture optimum, then run a held-out photon/aperture sweep.

**Nearest prior art:** optical étendue, mode counting, Fisher-information matching. The task-coupling exponent is the missing ingredient.

## Rank 7 — Dual-wall self-calibration

**Score:** `4 × 4 × 4 = 64`

**Mechanism:** support-wall leakage responds to relay/pixel geometry; energy-wall leakage responds to collection clipping. Their two calibrated breaking coordinates identify which hardware symmetry failed.

**Exact sketch:** form two score coordinates from the support-wall and unitary-energy-wall breakers; their Jacobian with respect to `(z1,NA)` is locally invertible when the breaking vectors are noncollinear.

**Kill test:** joint `(z1,NA)` grid, estimate the 2×2 breaking Jacobian, and require condition number `<10` plus held-out geometry recovery within `10%`.

**Nearest prior art:** self-calibration and null metrology. The novelty is using two exact blind channels as orthogonal hardware diagnostics.

## Rank 8 — Certificate-transfer under drift

**Score:** `4 × 5 × 4 = 80`

**Mechanism:** Davis–Kahan/Weyl perturbation bounds convert calibration drift into a guaranteed loss of blind dimension.

**Exact sketch:** if `||L(z)-L(z0)||<=Delta(z)`, then a subspace certified at `sigma_d(z0)` remains blind to `sigma_d(z0)+Delta(z)`; principal-angle gaps control subspace rotation.

**Kill test:** calibrate at one plane and test held-out `z1`, fill, wavelength, and DMD-render perturbations. The predicted envelope must cover every observed leak.

**Nearest prior art:** robust subspace tracking. The optical consequence is a time-valid cannot-see certificate.

## Rank 9 — Wall-selective sequential sensing

**Score:** `4 × 4 × 4 = 64`

**Mechanism:** run separate CUSUMs on symmetry-breaking score coordinates; a wall-preserving event stays dark by theorem, while a wall-breaking event produces a class-specific delay exponent.

**Kill test:** simultaneous scene/hardware/medium events with controlled wall classes; require empirical delay ordering to match the quotient-jet orders and cross-trigger probability below `0.1` inside the certified domain.

**Nearest prior art:** multi-chart CUSUM and fault isolation. The new piece is theorem-defined event classes from physical invariances.

---

# B11. MOONSHOT DERIVATION — the Three-Ledger Certifiability Theorem

## B11.1 Model and nuisance quotient

Let one independent bank yield

```text
Y_t ~ N(0, Sigma_epsilon),   t=1,...,T.
```

Whiten the null covariance and project every admitted nuisance tangent in the Gaussian covariance metric. The remaining efficient covariance displacement is

```text
A_perp(epsilon,g)
 = Pi_perp[ Sigma_0^(-1/2)(Sigma_epsilon-Sigma_0)Sigma_0^(-1/2) ].   (11.1)
```

Let the efficient symmetric tangent space `H` have dimension `p`, with basis `B_1,...,B_p` normalized by

```text
(1/2) tr(B_j B_k)=delta_jk.                                    (11.2)
```

Write

```text
A_perp = sum_j theta_j B_j,
lambda = ||theta||^2 = (1/2)||A_perp||_F^2.                    (11.3)
```

The local score coordinate for one bank is

```text
s_j(Y) = (1/2)[Y^T B_j Y-tr(B_j)].                              (11.4)
```

At the null,

```text
E s=0,
Cov(s)=I_p.                                                     (11.5)
```

For small `A_perp`, `E_A s = theta+o(||theta||)`. Therefore the normalized aggregate score experiment is locally

```text
Z_T = T^(-1/2) sum_t s(Y_t)
 approximately N(sqrt(T) theta, I_p).                           (11.6)
```

This is the nuisance-profiled LAN reduction. It is the only reduction used below.

## B11.2 Ledger I — known-direction detection

If the unit direction `u=theta/||theta||` is known independently of the data, the matched score

```text
u = u^T Z_T                                                     (11.7)
```

obeys

```text
H0: nu ~ N(0,1),
H1: nu ~ N(sqrt(T lambda),1).                                   (11.8)
```

Hence fixed nontrivial power requires and is achieved by

```text
T lambda = O(1).                                                (11.9)
```

This is the oracle channel ledger. The matched score is the equality case of the fluctuation–response/Cauchy–Schwarz bound.

## B11.3 Ledger II — blind existence testing

Now suppose only that `theta` lies somewhere on the `p`-sphere `||theta||^2=lambda`; its direction is unknown. The invariant statistic is

```text
Q_T = ||Z_T||^2.                                                (11.10)
```

Under the null,

```text
Q_T ~ chi^2_p,
E_0 Q_T=p,
SD_0(Q_T)=sqrt(2p).                                             (11.11)
```

Under an alternative,

```text
Q_T ~ chi^2_p(T lambda),
E_1 Q_T=p+T lambda.                                             (11.12)
```

Therefore the norm test separates the hypotheses when

```text
T lambda >> sqrt(p).                                            (11.13)
```

### Minimax lower bound

In the Gaussian score experiment, mix the alternative uniformly over directions with norm `r=sqrt(T lambda)`. For two independent mixture directions `mu,mu'`, the chi-square divergence of the mixture from the null is

```text
chi^2+1 = E exp(mu^T mu').                                      (11.14)
```

For uniform sphere directions,

```text
E(mu^T mu')=0,
Var(mu^T mu')=r^4/p.                                            (11.15)
```

If

```text
r^4/p = T^2 lambda^2/p -> 0,                                   (11.16)
```

then the mixture chi-square divergence tends to zero, hence total variation tends to zero and every test has asymptotic sum of errors approaching one. Thus no direction-agnostic procedure succeeds below

```text
T lambda ~ sqrt(p).                                             (11.17)
```

The norm/score-energy test attains the same scale. This proves the blind-certification tier.

## B11.4 Ledger III — estimate or localize the direction

Equation (11.6) gives the efficient estimator

```text
theta_hat = Z_T/sqrt(T).                                        (11.18)
```

Its exact local risk is

```text
E||theta_hat-theta||^2 = p/T.                                  (11.19)
```

The Cramér–Rao bound gives the same lower bound for regular unbiased/local-minimax estimation. To obtain relative squared error below a constant requires

```text
T lambda ~ p.                                                   (11.20)
```

Thus

```text
known direction : blind existence : full direction
       1        :     sqrt(p)     :       p.                    (11.21)
```

## B11.5 Add the physical jet and geometry

Suppose the first nonzero efficient physical displacement has the quotient-jet form

```text
A_perp(epsilon,g)
 = kappa g^a epsilon^m B + o(g^a epsilon^m),
(1/2)||B||_F^2=1.                                               (11.22)
```

Then

```text
lambda(epsilon,g)
 = kappa^2 g^(2a) epsilon^(2m)+o(...).                          (11.23)
```

The three bank requirements become

```text
T_known  asymp kappa^-2 g^(-2a) epsilon^(-2m),
T_blind  asymp sqrt(p) T_known,
T_est    asymp p T_known.                                      (11.24)
```

The physical exponent `2m` is unchanged by adaptation. The estimator ledger multiplies the prefactor by `sqrt(p)` or `p`.

## B11.6 Equality and boundary cases

1. **`p=1`:** all three ledgers coincide. There is no adaptation penalty.
2. **Known tangent:** the matched score attains the oracle scale.
3. **Isotropic unknown direction:** the spherical mixture gives the lower bound and the score-norm test attains it; this is the minimax equality case.
4. **Full estimation:** the efficient score estimator attains `p/T` risk.
5. **Exact wall:** `B=0` at every order; `lambda=0` and no finite `T` succeeds.
6. **Nuisance collapse:** if the physical derivative lies in the nuisance tangent, projection makes `B=0`; a hardware or statistical anchor is required.
7. **Structured/sparse direction:** replace `p` by the metric entropy/effective dimension of the admitted alternative class; the dense `sqrt(p)` bound is not universal for sparse alternatives.
8. **Non-Gaussian banks:** the theorem survives under LAN with the efficient-score covariance replacing the Gaussian basis; constants and remainder conditions change.

## B11.7 Falsifying numerical experiment

Run the same nuisance-profiled covariance family at

```text
M in {8,16,32,64},
p_eff in {10,36,136,528} by controlled tangent restriction,
T on a logarithmic grid,
lambda on a logarithmic grid.
```

For each cell evaluate:

1. matched known-direction score;
2. direction-agnostic score norm / unbiased quadratic U-statistic;
3. full tangent estimator MSE.

### Frozen falsification bars

- 50%-power contours for the matched score collapse against `T lambda` with CV `<=20%`.
- blind contours collapse against `T lambda/sqrt(p)` with CV `<=25%` and fail to collapse against `T lambda` alone;
- relative estimation risk collapses against `T lambda/p` with CV `<=25%`;
- fitted `p` exponents are `0±0.15`, `1/2±0.15`, and `1±0.15` for the three ledgers;
- the KT6 operating point remains below blind power `0.2` at FPR `0.05` under the optimal blind statistic.

Any violation kills the optical Three-Ledger claim while leaving the individual wall and jet theorems intact.

## B11.8 Why this is the deepest candidate

It survives every failure mode that killed the earlier spectral-duality idea:

- it does not hold trace fixed;
- it allows physical coefficients to co-vary arbitrarily;
- it is invariant after nuisance quotienting;
- it distinguishes oracle, blind, and estimation tasks rather than comparing them through one scalar spectrum;
- it predicts the KT6 negative quantitatively;
- it composes directly with the validated quotient-jet exponent.

The candidate broad statement is:

> **A physical perturbation has three independent prices: symmetry decides whether it exists in the record, jet order decides how its signal scales, and tangent dimension decides whether that signal can be blindly certified.**

---

# PART C — PAPER-2 ARCHITECTURE, INFORMED BY THE DIG

## C1. Center-of-gravity ruling

> **Choose option (c), provisionally: Paper 2 should be a Certified Blindness / Three-Ledger paper, not a memory-effect phase diagram and not a catalog of wall-repair tricks.**

The one protagonist is:

> **Engineering and certifying what an optical sensor cannot see.**

The causal spine is

```text
SYMMETRY CREATES A WALL
    -> CALIBRATION MEASURES ITS BLINDNESS CAPACITY
    -> GEOMETRY/PROJECTION SELECTS THE JET ORDER
    -> FINITE TANGENT DIMENSION SETS WHAT CAN BE CERTIFIED.
```

This unifies all of the day’s positive pillars:

- hard-pupil field wall `7e-16`;
- lossless-energy wall `2.18e-16`;
- 80-dimensional conjugate-plane certified code family at `<=1e-5`;
- DPSS annular suppression `<=2.4e-4`, honestly noncomposable with the code route;
- leak projection restoring slope `4.00` while retaining `86–99%` target information;
- the KT6 estimator kill as the third-ledger boundary;
- the calibrated physical residual law as the bench-transfer object.

The certified-wall engineering option (a) becomes the experimental body of this paper. The Three-Ledger theorem supplies the conceptual spine and prevents the paper from reading like three unrelated repair methods.

## C1.1 Why not the memory-effect phase diagram now

KT4 killed the proposed `z2` disorder axis: a single thin screen remains highly multimode and never approaches `R=1`. The correct family requires screen stacking, thickness, or a volumetric/random-matrix propagation model. That is a real future paper, but it is not the next paper on the present evidence.

Its named gate is **COHERENCE-CROSSOVER**: the code-count sweep must exhibit the predicted unsaturation-to-saturation transition and agree with an independent coherence-rank estimate. Until then, the memory-effect phase diagram stays future work.

## C1.2 Binding gate before Paper 2 is frozen

The new center requires one compact decisive exam:

> **TLSG — THREE-LEDGER SCALING GATE.**

It is exactly the falsification experiment in §B11.7 plus a calibration-transfer arm for Eq. (4.4).

### TLSG pass bars

1. known/blind/estimation contours scale as `p^0`, `p^1/2`, and `p^1` within the frozen tolerances;
2. KT6 remains below optimal blind certifiability, confirming the negative for the right reason;
3. a fresh calibration/test split validates `sigma_d(Lhat)+delta` as a simultaneous leak envelope;
4. at the conjugate plane, at least `64` code dimensions remain certified below `1e-4` on fresh physical-law draws;
5. no repair round.

If TLSG passes, the Paper-2 architecture below is frozen. If it fails, fall back cleanly to option (a), a narrower Optica paper on two exact walls and three restoration routes; do not preserve the Three-Ledger headline by rhetoric.

---

# C2. Frozen theorem set after TLSG

## Main-text theorems

1. **Statistical Noether Wall Theorem** — measurement-law invariance generates exact score-null orbits; local integrable score-null distributions are gauge leaves.
2. **Blindness Capacity Theorem** — the dimension–leak frontier is the bottom generalized singular spectrum; calibration uncertainty gives a finite metrological envelope.
3. **Three-Ledger Certifiability Theorem** — `1:sqrt(p):p` known/blind/estimation scaling, composed with quotient-jet order.

## Main-text corollary

4. **Jet Transmutation** — a lower-order physical leak changes the strict class from `m=2` to `m=1`; an exact wall or efficient projection restores `m=2`.

## Supplement/future objects

- DPSS versus SVD noncomposability proof and full operator spectra;
- task-dependent étendue law;
- coherence-rank unsaturation model;
- IT7 coincidence-veto negative;
- segmented-bucket negative and PWM apparatus constraint.

Do not make the paper carry every killed shortcut.

---

# C3. Hero figure and section architecture

## Hero figure — “What the sensor cannot see”

One horizontal figure with three linked panels:

### Left — symmetry walls

Two exact examples:

```text
hard-pupil support orbit       -> 7e-16 field wall
full-aperture unitary orbit    -> 2.18e-16 energy wall
```

Each has one controlled breaker: finite window/relay leakage and NA clipping.

### Center — blindness capacity

Plot generalized singular value versus subspace dimension. Mark:

```text
80 dimensions <=1e-5 at fixed conjugate plane,
35 dimensions <=1e-4 jointly over four z1 planes,
DPSS route <=2.4e-4 but joint capacity 80->24.
```

The visual message is “blindness has capacity and is not automatically composable.”

### Right — certifiability ladder

Three parallel thresholds on one axis:

```text
known direction  T lambda ~ 1
blind existence T lambda ~ sqrt(p)
estimation      T lambda ~ p
```

Place KT6 below the blind line and the projected jet above its channel line. This panel is added only after TLSG passes.

## Printed sections

1. **Symmetry walls of an optical measurement** — theorem plus two exact walls.
2. **Measuring blindness capacity** — calibrated operator, dimension–leak frontier, robustness.
3. **Restoring the desired jet** — optical, code-space, and statistical routes; noncomposability.
4. **What finite data can certify** — Three-Ledger theorem and KT6 boundary.
5. **Wave-optics transfer and apparatus constraints** — calibrated leak law, DPSS/SVD choice, single bucket, PWM.

No chronology of 17 tests appears in the article.

---

# C4. Titles and venue

## Preferred title after TLSG

> **Engineering What an Optical Sensor Cannot See**

Alternatives:

1. **Certified Blindness in Computational Optics**
2. **Symmetry Walls and Certifiability in Optical Sensing**
3. **The Blindness Capacity of an Optical Measurement**
4. **From Exact Optical Walls to Certifiable Change Detection**

## Venue

> **Optica first; Physical Review Applied second. Do not lead with PRX on the present evidence.**

### Why Optica

- two exact optical walls of different physical type;
- a calibrated code-space capacity result;
- three physically distinct wall-restoration routes;
- a wave-optics/Maxwell transfer chain;
- a broad but clearly optical metrology thesis.

### Why not PRX yet

The abstract ingredients—gauge nulls, singular-value min–max, LAN/Le Cam scaling—are classical individually. The novelty is their synthesis into a physical optical metrology program. Without an experimental instrument or a second non-optical physical realization, PRX would be a reach and could obscure the strongest optics contribution.

PRX becomes arguable only if TLSG reveals a genuinely new finite-sample law beyond the standard dense Gaussian minimax hierarchy and the result is demonstrated in a second physical model. Do not promise that outcome.

### Physical Review Applied fallback

PRApplied is the right fallback if Optica judges the paper too control/statistics-heavy. It would emphasize the calibrated wall, hardware routes, and bench protocol while retaining the theorem set.

---

# C5. Campaign ruling beyond the completed battery

Authorize **only TLSG**, because it decides the center. Do not start:

- the screen-stacking/memory-effect campaign;
- a new D3/D5 hardware repair;
- channel-ratio ranging;
- segmented-bucket designs;
- new code optimization;
- higher cumulants.

If TLSG passes, writing starts. If it fails, the narrow wall-engineering paper starts. The architecture no longer needs another open-ended divergence round.

---

# C6. D3/D5 and IT7 disposition

The Paper-1 specificity dent remains exactly as disclosed. IT7 does not reopen it.

> The contact-arm veto is a valid medium-event discriminator but an invalid beyond-band repair: it reduced medium false alarms from `0.9965` to `0.056` while retaining only `24%` of scene power because the same wall that supplies specificity starves the contact arm of beyond-band signal.

Ruling:

- no strong-specificity wording returns;
- the semiparametric projection remains a quarantined future method until estimator-layer certification is solved;
- IT7 belongs in the Paper-2 supplement and bench manual as a theorem-consistent negative control;
- it does **not** merit a standalone methods note.

Its useful residual application—rejecting in-band/amplitude classes—is an apparatus diagnostic, not a paper headline.

---

# C7. Final Letter transfer-sentence ruling

Replace only the current generic leak sentence in End Matter C with the following calibrated wording:

> **In the unfiltered complementary-DMD twin, the algebraic mean wall became a reproducible calibrated residual \(L(z_1)=6.94\times10^{-3}+5.92\times10^{-3}(z_1/\mathrm{mm})^{0.90}\); the conjugate-plane floor is set by fill-factor and macro-pixel staircase effects, and the fitted floor, coefficient, and exponent were stable across code seeds and code counts with coefficients of variation \(2.2\%\), \(2.2\%\), and \(0.1\%\). Near-conjugate relay imaging suppresses the propagation term.**

Retain the already-written following statements that a hard Fourier pupil plus guarded sampling is required for an exact physical support wall and that the sealed latency/grid does not map one-to-one onto a particular bench.

## PWM line

> **Keep PWM temporal dithering out of the Letter transfer paragraph.**

It is renderer-specific apparatus discipline, not part of Rank–Jet Separation or the transfer validity boundary. Put it in the bench protocol and, at most, the future Paper-2 apparatus capsule:

```text
Binary gray rendering must use full-period PWM with whole-period bucket integration;
ordered spatial dithering is disallowed because it raises the conjugate-plane leak to 1.38e-2.
```

Paper 1 otherwise remains frozen.

---

# Final combined ruling

> **The deep object is a Symmetry–Jet–Certifiability hierarchy. Exact measurement symmetries create walls; the first nonzero quotient jet fixes the physical amplitude law; unknown tangent dimension imposes a separate `1:sqrt(p):p` certifiability ladder. Paper 2 should provisionally center on certified blindness—engineering and publishing what an optical sensor cannot see—subject to one decisive THREE-LEDGER SCALING GATE. The memory-effect phase diagram is deferred to a genuine thickness/stacking axis. Optica is the honest first venue; PR Applied is the fallback.**
