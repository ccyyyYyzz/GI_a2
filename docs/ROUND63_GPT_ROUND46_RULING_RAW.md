# R46 — The Curved Blindness Paradox: a camera that cannot image but notices every event

**Reference brief:** [`docs/ROUND63_GPT_ROUND46_QUESTION.md`](https://github.com/ccyyyYyzz/GI_a2/blob/352c300/docs/ROUND63_GPT_ROUND46_QUESTION.md) at [`352c300`](https://github.com/ccyyyYyzz/GI_a2/commit/352c300).  
**Round type:** unbound invention. This does not alter the frozen Letter, reopen a dead line, or constrain TLSG.

## Executive synthesis — the thing with 感觉

The strongest jump is not another sensor, estimator, or wall certificate. It is a dynamical paradox hiding inside the validated walls:

> **A one-number sensor can hide an entire `N−1`-dimensional world perfectly and still notice every straight-line event. To remain invisible while moving, the event must follow the curved blindness orbit and pay an exact centripetal-acceleration tax.**

This is the **Curved Blindness Paradox**. It turns the Noether-wall result from static metrology into a law of motion.

For the completely scrambled channel, the record depends on the scene only through

```text
Q(x)=x^T G x,    G>0.
```

Every ellipsoid `Q(x)=q` is an exact global blindness manifold: all of its points generate the same complete record law, even with infinite data. Yet every nonzero straight perturbation `x+epsilon delta` leaves that manifold at first or second order. Generic directions are detected with `T~epsilon^-2`; tangent directions are curvature-rescued with `T~epsilon^-4`. The only way to make a moving event exactly invisible is to steer it along the ellipsoid. The minimum required acceleration is not a metaphor; it is

```text
A_* = ||v||_G^2 / ||x||_G.
```

If the event cannot supply that normal acceleration, the unavoidable Chernoff leakage is exactly controlled by the acceleration deficit.

That gives a dinner-table sentence a physicist can retell:

> **Blind spots are curved roads, not directions. A state may be perfectly hidden, but motion stays hidden only if it pays the road’s curvature.**

It also yields two hard capabilities:

1. a **zero-knowledge camera** that can prove an event occurred while being mathematically incapable of forming the picture; and
2. a **Noether event cloak** with an exact speed/acceleration limit, rather than an informal “stealth subspace.”

---

# Top three — ranked by goosebump × rigor × sim reach

| Rank | Candidate | Goosebump | Rigor | Sim reach | Product |
|---|---|---:|---:|---:|---:|
| **1** | **Curved Blindness / Noether Event Cloaks** | 5 | 5 | 5 | **125** |
| **2** | **The Fog Quotient Channel: calibration-free communication through complete scrambling** | 5 | 4 | 5 | **100** |
| **3** | **The Law-Only Quadratic Computer** | 4 | 4 | 5 | **80** |

The first candidate is the only one I would currently call a new physical principle. The second is a powerful configuration-level result but sits inside noncoherent communication unless its image-erasing quotient property proves decisive. The third is exciting engineering but is dangerously close to random-feature and reservoir computing.

---

# 1. Rank 1 — Curved Blindness / Noether Event Cloaks

## 1.1 Dinner-table sentence

> **A camera can be mathematically unable to make a picture—even with infinite data—yet detect every local event; the only invisible motions are those that curve exactly along its Noether blindness orbit.**

## 1.2 What is actually new

The static ingredients have classical relatives:

- nonlinear observability and output-nulling trajectories;
- zero dynamics and stealth attacks in monitored systems;
- group-invariant statistics and maximal invariants;
- constrained motion on a curved manifold.

The candidate new theorem is their missing junction:

> **The second fundamental form of a statistical blindness orbit is the minimum control required to keep an event invisible, and any shortfall produces a universal fourth-order Chernoff leakage with a closed-form equality case.**

This is not Efron’s statistical curvature, which measures curvature of a family of probability laws and second-order estimator efficiency ([Efron 1975, DOI 10.1214/aos/1176343282](https://doi.org/10.1214/aos/1176343282)). It is curvature of the **state-space fiber on which the complete law is constant**, converted into a detection exponent. Nor is it merely the existence of output-nulling motions from geometric control ([Hermann–Krener 1977, DOI 10.1109/TAC.1977.1101601](https://doi.org/10.1109/TAC.1977.1101601)) or stealthy attacks in cyber-physical systems ([Pasqualetti–Dörfler–Bullo 2013, DOI 10.1109/TAC.2013.2266831](https://doi.org/10.1109/TAC.2013.2266831)). The new object is the **exact effort–leakage law** between a static null orbit and a moving physical event.

---

# 2. MOONSHOT DERIVATION — the Curvature Tax of Undetectability

## 2.1 General setting: a blind fiber of a complete statistical experiment

Let the physical state be `x in R^N`. The full record law depends on `x` through a smooth invariant

```text
h:R^N -> R^r,
P_x = P_{h(x)},                                                  (2.1)
```

after every admitted nuisance has already been quotient-projected. A level set

```text
M_c = {x : h(x)=c}                                               (2.2)
```

is an exact blindness manifold: all states on it generate the same complete law.

Let a physical event follow a twice differentiable trajectory

```text
x(t)=x+t v+(t^2/2)a+O(t^3),                                     (2.3)
```

with velocity `v` and acceleration `a`. Expanding the invariant,

```text
h(x(t))-h(x)
 = t Dh_x v
 + (t^2/2)[Dh_x a + D^2h_x[v,v]]
 + O(t^3).                                                       (2.4)
```

The event is first-order invisible exactly when

```text
Dh_x v=0,                                                       (2.5)
```

so `v` lies in the tangent space of the blind manifold. To stay invisible to second order it must additionally satisfy

```text
Dh_x a = -D^2h_x[v,v].                                         (2.6)
```

Equation (2.6) is the dynamical content missing from a static null-space certificate. A tangent direction is free; a tangent **motion** is not. The required normal acceleration is the second fundamental form of the blindness fiber.

## 2.2 Theorem 1 — minimum acceleration for second-order event cloaking

Let `W>0` define the physical control metric

```text
||a||_W^2 = a^T W a,
```

and assume `Dh_x` has full row rank. Set

```text
b = D^2h_x[v,v],
S = Dh_x W^-1 Dh_x^T.                                           (2.7)
```

Among all accelerations that maintain blindness through second order, the unique minimum-norm acceleration is

```text
a_* = -W^-1 Dh_x^T S^-1 b,                                     (2.8)
```

and the exact curvature tax is

```text
||a_*||_W^2 = b^T S^-1 b.                                      (2.9)
```

### Proof

Minimize `(1/2)a^TWa` subject to `Dh_x a=-b`. The Lagrangian

```text
L(a,lambda)=(1/2)a^TWa+lambda^T(Dh_x a+b)
```

gives `Wa+Dh_x^Tlambda=0`. Thus `a=-W^-1Dh_x^Tlambda`; imposing the constraint yields `S lambda=b`, which gives (2.8). Substitution gives (2.9). `square`

The equality condition is exact: every other second-order cloak has larger control norm by the `W`-orthogonal decomposition into the minimum normal acceleration plus an arbitrary tangent acceleration.

For a regular level set, a local curve satisfying `h(x(t))=h(x)` exists through every tangent direction. Equation (2.8) is the minimum initial acceleration of such a curve. Higher derivatives determine how the cloak is continued beyond second order.

## 2.3 Theorem 2 — residual curvature becomes fourth-order statistical leakage

Let `J_h` be the efficient Fisher matrix of the reduced law `P_h`. For a small parameter displacement `Delta h`, the profiled Chernoff distance obeys

```text
C(P_{h+Delta h},P_h)
 = (1/8) Delta h^T J_h Delta h + o(||Delta h||^2).               (2.10)
```

For a first-order-invisible velocity satisfying (2.5), define the residual normal acceleration

```text
r(a;v)=Dh_x a + D^2h_x[v,v].                                   (2.11)
```

Then

```text
C(P_{x(t)},P_x)
 = (t^4/32) r(a;v)^T J_h r(a;v) + o(t^4).                       (2.12)
```

Thus an imperfect event cloak does not merely “leak a little.” It enters a universal curvature class:

```text
T_req ~ t^-4.                                                   (2.13)
```

The coefficient is the Fisher-weighted squared residual of the second fundamental form after the available control has tried to cancel it.

### Proof

Under (2.5), equation (2.4) gives `Delta h=(t^2/2)r+O(t^3)`. Substitute into (2.10). `square`

This theorem is coordinate-free in substance: reparameterizing `h` changes both `r` and `J_h` covariantly, leaving the quadratic form invariant.

## 2.4 Exact specialization to the validated scrambled channel

Take

```text
h(x)=Q(x)=x^T Gx,   G>0,                                       (2.14)
```

with inner product and norm

```text
<u,v>_G=u^TGv,
||u||_G=sqrt(u^TGu).
```

Write

```text
r_0=||x||_G,
s=||v||_G.
```

The blind manifolds are ellipsoids `Q(x)=r_0^2`. A velocity is tangent exactly when

```text
<x,v>_G=0.                                                      (2.15)
```

Along (2.3),

```text
Q(x(t))-Q(x)
 = 2t<x,v>_G
 + t^2[s^2+<x,a>_G]
 + O(t^3).                                                       (2.16)
```

For a tangent event, exact second-order cloaking requires

```text
<x,a>_G=-s^2.                                                   (2.17)
```

By Cauchy–Schwarz,

```text
s^2 = |<x,a>_G| <= r_0 ||a||_G,
```

so every exact cloak must obey

```text
||a||_G >= A_* = s^2/r_0.                                      (2.18)
```

Equality holds only when

```text
a_* = -(s^2/r_0^2)x.                                           (2.19)
```

This is the exact **centripetal acceleration of blindness**.

## 2.5 Exact all-orders cloak and equality case

The lower bound is not merely local. Define

```text
omega=s/r_0,

gamma(t)
 = x cos(omega t) + (v/omega) sin(omega t).                     (2.20)
```

Because `<x,v>_G=0` and `||v/omega||_G=r_0`,

```text
Q(gamma(t))
 = r_0^2 cos^2(omega t)+r_0^2 sin^2(omega t)
 = Q(x)                                                         (2.21)
```

for every `t`. Moreover `gamma(0)=x`, `gamma'(0)=v`, and

```text
gamma''(0)=-(s^2/r_0^2)x,
||gamma''(0)||_G=s^2/r_0.                                      (2.22)
```

So the equality trajectory is a great circle of the `G`-ellipsoid. It is perfectly invisible to the complete `Q`-channel at all times, not merely to first or second order.

This is the event cloak.

## 2.6 Theorem 3 — acceleration-limited stealth speed and unavoidable leakage

Suppose the event can supply at most

```text
||a||_G <= A.                                                   (2.23)
```

The most negative possible radial acceleration is `<x,a>_G=-A r_0`. Therefore the minimum unavoidable second-order displacement of `Q` is

```text
kappa_min(A)
 = [s^2-A r_0]_+,                                               (2.24)
```

where `[z]_+=max(z,0)`. Equality is attained by choosing `a=-(A/r_0)x` below threshold, and by the great-circle acceleration once `A>=A_*`.

Hence the exact stealth speed limit is

```text
s <= sqrt(A r_0).                                               (2.25)
```

A faster tangent event cannot remain invisible with the available acceleration.

Let `I_Q` be the efficient Fisher information for the scalar `Q`. The best possible Chernoff leakage under the acceleration cap is

```text
C_min(t;v,A)
 = (I_Q/8)[s^2-A r_0]_+^2 t^4 + o(t^4).                        (2.26)
```

For the Gaussian covariance channel

```text
Sigma(Q)=R+QH,
K=Sigma^-1/2 H Sigma^-1/2,
I_Q=(1/2)||K||_F^2,                                             (2.27)
```

so

```text
C_min
 = (||K||_F^2/16)[s^2-A r_0]_+^2 t^4+o(t^4).                   (2.28)
```

This is the moonshot identity:

> **The fourth-order detection coefficient is the square of the unpaid curvature tax.**

It supplies a hard boundary between three regimes:

```text
radial velocity present             -> T~t^-2;
tangent but under-actuated motion   -> T~t^-4;
exact orbit-following cloak         -> C=0 for all T.            (2.29)
```

## 2.7 Why this resolves the apparent paradox

The validated Rank–Jet result said that every nonzero straight-line direction is locally visible at order one or two. The blindness-orbit theorem says that an `N−1`-dimensional set of exact invisible motions nevertheless exists.

There is no contradiction:

- a tangent **line** leaves a curved ellipsoid at second order;
- a tangent **trajectory** can remain on it by adding the exact normal acceleration.

The state is globally nonidentifiable, the tangent is first-order null, the straight event is curvature-visible, and the controlled orbit motion is exactly cloaked. These are four different statements that standard “observable/unobservable” language conflates.

## 2.8 Perfect privacy corollary

Let

```text
O_G={U:U^TGU=G}
```

be the `G`-orthogonal group. The complete record is invariant under `x->Ux` because `Q(Ux)=Q(x)`.

Write a random state as `(R,Omega)`, where `R=Q(X)` and `Omega` is its orbit coordinate. For any number of banks `T`,

```text
Omega -> R -> Y^T                                               (2.30)
```

is a Markov chain. Therefore

```text
I(Omega;Y^T | R)=0.                                             (2.31)
```

For `M` equiprobable scenes on the same `Q`-orbit, every classifier—neural, Bayesian, computationally unbounded—has success probability exactly `1/M`. With a Haar-uniform orbit prior, the posterior orientation remains Haar, and the normalized Bayes reconstruction error does not decrease with `T`.

Yet every nonzero straight local event is detected at finite jet order. Thus:

> **The sensor has perfect infinite-data image privacy and an empty straight-event blind set.**

This is the theorem behind the “camera that cannot see.”

## 2.9 Boundaries that must be stated

1. **Anchoring.** If medium covariance amplitude is a free nuisance collinear with `Q`, then `I_Q=0` and every coefficient above vanishes. A wavelength, lag-shape, polarization, or calibrated amplitude anchor remains mandatory.
2. **Positive definiteness.** If `G` has a kernel, directions in `ker G` are truly blind and require no cloak acceleration.
3. **Local versus global.** A signed straight path can cross the same ellipsoid again at the already validated iso-`Q` cancellation amplitude. The theorem classifies local motion and exact orbit paths, not arbitrary finite chords.
4. **Control model.** The acceleration tax is with respect to the chosen physical metric `G`; another actuator metric uses the general formula (2.8)–(2.9).
5. **Dual use.** The same theorem supports privacy-preserving sensing and adversarial evasion. Publication should lead with measurement geometry and privacy, not operational recipes for defeating a specific security system.

## 2.10 Decisive first simulation — `CURVED_BLINDNESS_TEST`

Use the existing `SCRAMBLE_EXT/JET_TEST` exact Gaussian `Q` engine. No new wave solver is required.

### Frozen arms

Choose one baseline `x`, one exactly `G`-tangent velocity `v`, and define `A_*=s^2/r_0`.

1. **STRAIGHT:** `gamma_0(t)=x+tv`.
2. **PARTIAL CLOAK:** initial accelerations `a_alpha=-alpha(s^2/r_0^2)x`, with `alpha in {0.25,0.5,0.75}`.
3. **EXACT CLOAK:** the great circle (2.20).
4. **RADIAL CONTROL:** a generic velocity with `<x,v>_G !=0`.
5. **WRONG-METRIC CONTROL:** acceleration optimized in Euclidean norm rather than the declared actuator metric.

### Frozen predictions

- STRAIGHT and every PARTIAL arm: Chernoff/KL slope `4.00±0.05`.
- Coefficient ratio relative to STRAIGHT:

  ```text
  C_alpha/C_0 = (1-alpha)^2                                   (2.32)
  ```

  within `5%`.
- EXACT CLOAK: `|Delta Q|/Q < 1e-12`, exact divergence below numerical floor, and Monte Carlo AUC in `[0.47,0.53]` for every tested `t` and bank count.
- RADIAL CONTROL: slope `2.00±0.05`.
- Acceleration-cap sweep: the fitted fourth-order coefficient follows

  ```text
  [1-A/A_*]_+^2                                                (2.33)
  ```

  with kink location within `5%` of `A_*`.
- Sequential delay under partial cloaking scales as `t^-4`; exact cloak never crosses except at the false-alarm rate.

### Kill conditions

Kill or sharply demote the moonshot if any occurs:

- the great-circle record is distinguishable despite fixed `Q` under the exact channel;
- the coefficient law misses `(1-alpha)^2` by more than `10%` after numerical convergence;
- the acceleration threshold misses `A_*` by more than `10%`;
- nuisance profiling introduces a lower-order term not removed by the declared anchor;
- the result cannot be reproduced in a second non-optical statistical model with a curved blind fiber.

Runtime should be well below one GPU-hour for exact divergences and under half a GPU-day including Monte Carlo and CUSUM.

---

# 3. Rank 2 — The Fog Quotient Channel

## 3.1 Dinner-table sentence

> **A completely scrambling medium can erase every spatial symbol and still carry a calibration-free message: encode in a ratio of quadratic energies, and the unknown fog cancels itself.**

## 3.2 Exact channel law

At the fully scrambled endpoint, whiten the informative covariance subspace into `k` effective real Gaussian looks per bank. Suppose a reference block and a message block share an unknown multiplicative medium scale `a`, while the sender controls the scene invariant `q=Q(x)`.

With `T_0,T_1` banks and degrees of freedom

```text
nu_0=kT_0,
nu_1=kT_1,
```

the sufficient energies obey

```text
S_0/(a q_0) ~ chi^2_{nu_0},
S_1/(a q_m) ~ chi^2_{nu_1}.                                    (3.1)
```

Therefore

```text
R=(S_1/nu_1)/(S_0/nu_0)
  = rho_m F_{nu_1,nu_0},
rho_m=q_m/q_0.                                                  (3.2)
```

The unknown medium amplitude `a` cancels exactly. In log coordinates,

```text
Z=log R = log rho_m + N,                                       (3.3)
```

where `N=log F`. For equal block sizes `nu_0=nu_1=nu`,

```text
E[N]=0,
Var(N)=2 psi_1(nu/2) ~ 4/nu.                                   (3.4)
```

Thus the exact Shannon problem is a one-dimensional additive-noise channel

```text
C_nu(L)=sup_{P_U, 0<=U<=L} I(U;U+N_nu),                         (3.5)
```

with `U=log rho`. It can be solved numerically without a transmission matrix, wavefront estimate, or medium-amplitude calibration.

The physically interesting version holds the mean bucket fixed while changing `Q`: symbols have equal DC/total flux but different quadratic texture energy. The mean channel then carries no symbol, while covariance does. This requires an isotropic or declared `G` law; arbitrary unknown `G` does not permit transmitter-side control of `Q` at fixed mean.

## 3.3 Decisive first simulation

Use `SCRAMBLE_EXT` with medium amplitude randomized over `40 dB` between codewords, but held common within each reference/message pair.

- alphabets `rho in {1,1.25,1.6,2.2,3.2}`;
- `T in {16,64,256,1024}`;
- absolute-energy decoder, pilot-estimated decoder, and exact ratio decoder;
- both global-scale symbols and constant-DC texture symbols;
- exact `F`-law transition matrix and Blahut–Arimoto capacity.

Pass only if:

1. empirical log-ratio noise matches the predicted `log F` law by a predeclared KS/CvM bar;
2. BER and mutual information are invariant to the unknown amplitude distribution;
3. constant-DC symbols remain invisible to the mean channel;
4. rates remain nontrivial under the measured `M_eff≈13` bank cost.

## 3.4 Prior art to beat

This is adjacent to noncoherent block-fading capacity, differential modulation, and energy-based massive-SIMO communication. Energy modulation without instantaneous channel state is established; for example, Manolakos, Chowdhury, and Goldsmith develop energy-based noncoherent modulation using channel statistics ([arXiv:1507.04978](https://arxiv.org/abs/1507.04978)). Communication through scattering media without a full transmission matrix is also established, including adaptive channel formation ([Hao, Martin-Rouault, and Cui 2014, DOI 10.1038/srep05874](https://doi.org/10.1038/srep05874)).

The surviving novelty is narrower:

> a **true image-erasing, `Q`-only covariance channel** through complete dynamic scrambling, with exact differential cancellation of the unknown medium scale and a mean-invariant texture alphabet.

Without the mean-invariant alphabet and the exact nuisance-cancelled capacity law, this is a renaming of noncoherent amplitude communication and should be killed as a flagship.

---

# 4. Rank 3 — The Law-Only Quadratic Computer

## 4.1 Dinner-table sentence

> **The fog can destroy the picture and still compute its quadratic kernel at light speed.**

## 4.2 Exact random-feature statement

Suppose `L` independently addressable medium laws—wavelengths, polarizations, lags, or reconfigurable scattering states—supply PSD kernels `G_l`, and the covariance readout returns

```text
q_l(x)=x^T G_l x.                                               (4.1)
```

For the canonical random rank-one case

```text
G_l=a_l a_l^T,
a_l~N(0,I),                                                     (4.2)
```

define the centered feature

```text
z_l(x)=(a_l^Tx)^2-||x||^2.                                     (4.3)
```

Isserlis’ theorem gives the exact kernel identity

```text
E[z_l(x)z_l(y)] = 2(x^Ty)^2.                                   (4.4)
```

Thus a linear electronic readout of law-only covariance scalars implements a degree-two polynomial kernel without measuring or reconstructing the input. With enough independent kernels, the empirical feature Gram matrix converges at the usual `L^-1/2` rate. More generally, any symmetric quadratic classifier

```text
x^TAx
```

can be written as the difference of two PSD forms using `A=A_+-A_-`; two appropriately engineered covariance laws compute it exactly.

The scientifically sharp claim would be:

> **A dynamic scattering medium can be used as a realization-free quadratic-feature machine: individual speckles are never calibrated, and only their stable law is read through a zero-dimensional covariance detector.**

## 4.3 Decisive first simulation

Extend the scrambling engine to `L in {1,4,16,64,256}` independent law kernels.

- Verify (4.4) on held-out vectors.
- Compare empirical kernel error with `L^-1/2`.
- Solve an XOR/parity or quadratic-boundary classification task where every linear bucket statistic is provably insufficient.
- Compare a camera-speckle random-feature arm, the `0-D` law-only covariance arm, and a digital random-feature baseline under matched photon and MAC budgets.
- Replace unconstrained random `G_l` by physically admissible Toeplitz grain kernels and measure the feature-rank collapse.

Kill the general claim if physically admissible kernels span only a very low-dimensional cone or if independent-bank cost overwhelms the optical compute advantage.

## 4.4 Prior art to beat

Random scattering and nonlinear disordered media are already major optical-computing platforms. Large-scale nonlinear photonic computing with disordered media has been demonstrated ([Nature Computational Science 2024](https://doi.org/10.1038/s43588-024-00644-1)), and nonlinear scattering tensors have been measured and proposed for optical information processing ([Moon et al., Nature Physics 2023, DOI 10.1038/s41567-023-02163-8](https://doi.org/10.1038/s41567-023-02163-8)). Optical random projections are also established.

The only defensible new square is:

```text
0-D bucket
× unmeasured realizations
× declared-law covariance readout
× exact quadratic-kernel theorem
× no input reconstruction.
```

If that square does not deliver a scaling or robustness advantage, the idea is reservoir/random-feature computing under a new noun and should be archived.

---

# 5. Physical zero-knowledge sensing — the strongest corollary of Rank 1

Although it is folded into Rank 1 rather than counted separately, this is probably the most powerful public framing.

## 5.1 Definition

A sensing transcript `Y` is **perfectly witness-private for property q** when

```text
P(Y|x)=P(Y|q(x))                                                (5.1)
```

for every admissible state. A simulator given only `q` can then generate an exactly distributed transcript. For any two witnesses `x,x'` with the same property,

```text
P(Y|x)=P(Y|x').                                                 (5.2)
```

This is an honest-verifier statistical zero-knowledge condition in the literal information-theoretic sense, although the system is not an interactive cryptographic proof.

In the scrambled channel, `q=Q(x)`. The transcript can certify a change in `Q`, while revealing zero information about the orbit coordinate. A second noncollinear kernel acts as an authorization key: it reduces the privacy-orbit dimension and reveals additional structure.

## 5.2 Dinner-table sentence

> **The fog can prove that something happened while revealing literally nothing about what the object looked like.**

## 5.3 Why privacy-optics precedent does not own it

Trainable optical privacy kernels and class-specific diffractive cameras already suppress or erase selected visual content; examples include Sepehri et al.’s trainable optical kernel ([arXiv:2106.14577](https://arxiv.org/abs/2106.14577)) and Bai et al.’s class-specific diffractive camera ([eLight 2022, DOI 10.1186/s43593-022-00021-3](https://doi.org/10.1186/s43593-022-00021-3)). Those are learned, distribution-dependent utility–privacy systems and are judged empirically against specified attacks.

The new theorem would be stronger and narrower:

- privacy holds against every algorithm and unlimited data;
- it is generated by an exact physical invariance, not training;
- the hidden dimension is known (`N−1` for one regular scalar invariant);
- local event completeness and its `epsilon^-2/epsilon^-4` rates are proved;
- the exact invisible trajectories and their minimum actuation cost are known.

---

# 6. Honest kill paragraph — what is secretly a known field

## 6.1 Event cloaking without the curvature theorem dies

“Design an undetectable event” by itself is output-nulling control, zero dynamics, or a stealthy attack. The control and cyber-physical-security literature already characterizes undetectable inputs and monitor limitations. The candidate survives only as the **Noether-orbit curvature theorem plus the exact acceleration-deficit Chernoff law**. Remove equations (2.8)–(2.28), and there is no new flagship.

## 6.2 The zero-knowledge camera without exact invariance dies

Privacy-preserving optical front ends, event-camera anonymization, and task-specific diffractive cameras are established. A learned system that merely makes reconstruction difficult is not new enough. The survivor is the exact Markov/simulator statement `I(Omega;Y|Q)=0` together with an empty straight-event blind set and a tunable symmetry-breaking anchor.

## 6.3 The fog communication idea is mostly noncoherent communication

If the sender merely scales brightness and the receiver energy-detects it, the result is classical noncoherent amplitude modulation. “No transmission matrix” is not enough. The new content must be the exact image-erasing quotient channel, constant-mean texture alphabet, and `F`-ratio nuisance cancellation. If those do not yield a measurable rate or privacy advantage, kill it.

## 6.4 The fog computer is mostly random features/reservoir computing

Scattering media as random linear maps, nonlinear speckle processors, reservoir computers, and optical kernels are occupied. The exact expectation (4.4) is classical Gaussian moment algebra. The candidate lives only if a zero-dimensional **law-only** covariance implementation reaches useful feature rank without realization calibration and beats matched digital/camera baselines.

## 6.5 A thermodynamic “first law of observability” dies for now

The validated record contains no conserved scalar that trades symmetry breaking, detection information, and certification cost. The ledgers are hierarchical, not additive: a wall can set information to zero; a jet sets an amplitude exponent; adaptation sets a dimension penalty. Calling that a first law would be metaphor, not physics. The Curvature Tax is the rigorous replacement: it is an effort–leakage law with equality conditions, not a pretend conservation law.

---

# 7. What to run first

> **Run only `CURVED_BLINDNESS_TEST` first.**

It is the smallest, most decisive, and most foundational experiment. It uses the exact engine that already validated Rank–Jet Separation, has parameter-free coefficient and threshold predictions, and can kill the entire frame in under a GPU-day.

A full pass would establish the following one-sentence result:

> **A symmetry can hide a state for free, but hiding motion requires paying the curvature of the blind orbit; the unpaid curvature is exactly the fourth-order information leaked to the detector.**

That is the first post-R42 idea in this program that I believe has genuine dinner-table retellability **and** a theorem sharp enough to survive the honesty machinery.

