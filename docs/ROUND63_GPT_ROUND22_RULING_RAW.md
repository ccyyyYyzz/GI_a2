# R22 ruling (GitHub issue #14, raw)

Title: R22 ruling: the controlled score-transfer object and limits of unification
Posted: 2026-07-19T21:58:12Z

---

# R22 — Compression ruling: the controlled score-transfer object and the limits of unification

**Audit target:** commit `4acfd34`, especially `docs/ROUND63_GPT_ROUND22_QUESTION.md`.

**Scope:** pure theory exploration only. Nothing here modifies either frozen manuscript, any campaign, or any frozen claim.

> **R22 VERDICT:** there is one honest object underneath the six results, but it is **not a new scalar law**. It is the **controlled statistical experiment**, represented locally by its **conditional-score Fisher Gramian**. The detector-side results are contractions of the complete-data score through a hidden-state observation channel; the design-side results are variational optimization of the resulting Gramian under a feasible control cone.
>
> This object compresses the system conceptually, but its abstract ingredients are classical: Fisher/Louis missing information, monotonicity under statistical transformations, Blackwell/Le Cam experiment comparison, and Kiefer–Wolfowitz optimal design. The publishable residual is therefore not “a new universal information geometry.” It would have to be a genuinely new theorem about **low-rank control of a dead-time score-transfer Gramian under physical constraints**.
>
> Of the three proposed candidates, **C2 is the best theorem candidate**, C3 is only a KKT/tangent-cone interpretation, and C1’s mechanism-by-mechanism product law is false in general. A fourth object, defined below, subsumes all three and is the only compression I recommend retaining.

---

## Executive ranking

Scores are deliberately harsh. “Compression” counts exact mathematical reductions, not loose narrative associations.

| candidate | theorem probability | compression over E1–E6 | novelty after prior art | R22 ruling |
|---|---:|---:|---:|---|
| **C1 — per-mechanism Friis/product law** | **1/5** | **2/6** | high *if true*, but the proposed form is false | Reject as a general scalar algebra. Replace by conditional-score channel composition; scalar products occur only in special invariant one-dimensional cases. |
| **C2 — control-rank/alignment hierarchy** | **4/5** | **4/6** | medium in the dead-time specialization; low at the abstract KKT level | Best of the submitted three. The correct object is a projected KKT residual, not `C_u` and not a scene-only Grassmann overlap. |
| **C3 — conjugate-plane obstruction** | **2/5** | **3/6** | low; mostly normal-cone geometry in different language | Demote to an interpretation. “Acts trivially on the binding plane” is sufficient in some cases, not necessary, and compensating controls immediately break it. |
| **C4 — controlled conditional-score Gramian** | **5/5** | **E1–E5 exactly at Fisher/local-design level; E6 mechanistically, not its dB values** | low as an abstract object; potentially medium for a new dead-time rank-control theorem | Adopt as the single essential object. It retires C1’s product algebra, C2’s proposed `A_k`, and C3’s plane metaphor. |

---

# Q1 — Candidate adjudication and prior-art boundary

## 1. C1: the proposed Friis algebra is not structural

### 1.1 What *is* structural

Let `X` be a complete latent record with score `S(X)`, and let a detector mechanism be a Markov channel `K` producing an observation `Y`. Define the conditional-score operator

\[
\mathsf C_K f(Y)=\mathbb E[f(X)\mid Y].
\]

For two cascaded channels

\[
X\longrightarrow Y\longrightarrow Z,
\]

the tower property gives

\[
\mathsf C_{X\to Z}=\mathsf C_{Y\to Z}\,\mathsf C_{X\to Y}.
\]

For one fixed scalar score direction `s`, define its Fisher contraction

\[
\eta_s(K)=\frac{\|\mathsf C_K s\|_2^2}{\|s\|_2^2}.
\]

Then the exact chain rule is

\[
\boxed{
\eta_s(K_2K_1)
=
\eta_s(K_1)\,
\eta_{\mathsf C_{K_1}s}(K_2).
}
\]

The second factor is evaluated on the **score already reshaped by the first channel**. It is not a mechanism constant. A mechanism-independent product occurs only when the relevant score lies in a common one-dimensional invariant/eigen-score subspace of every contraction. That is the exceptional case, not the general algebra.

This is the Fisher-local counterpart of information loss under garbling/statistical transformations, already situated inside the Blackwell/Le Cam and Chentsov traditions. Louis’s missing-information formula is the incomplete-data version ([Louis 1982](https://doi.org/10.1111/j.2517-6161.1982.tb01203.x)); Fisher monotonicity under Markov morphisms is part of information geometry ([Ay, Jost, Lê & Schwachhöfer 2015](https://doi.org/10.1007/s00440-014-0574-8)).

### 1.2 A third independent hold mechanism does **not** create a third factor

Take a renewal cycle

\[
C_\lambda=E_\lambda+\sum_{j=1}^{m} B_j,
\qquad E_\lambda\sim\mathrm{Exp}(\lambda),
\]

where the hidden delays `B_j` are independent of `lambda` and of one another. Let

\[
b=\sum_j\mathbb E B_j,
\qquad
v=\sum_j\operatorname{Var}B_j,
\qquad
\rho=\lambda b,
\qquad
c^2=v/b^2.
\]

The renewal local-normal calculation gives

\[
\mathbb E N_T\sim\frac{T}{b+1/\lambda},
\qquad
\operatorname{Var}N_T\sim
\frac{(v+1/\lambda^2)T}{(b+1/\lambda)^3},
\]

and hence, for `theta=log lambda`, normalized by `nu=T/b`,

\[
\boxed{
J_\infty(\rho,c)
=
\frac{\rho}{(1+\rho)(1+c^2\rho^2)}.
}
\]

Adding `B_2`, `B_3`, and so on changes only the **aggregate mean** `b` and **aggregate variance** `v`. It does not produce

\[
\prod_j(1+c_j^2\rho^2)^{-1}.
\]

The two denominator factors in E2 correspond to two cycle cumulants—mean cycle occupancy and cycle variance—not to one factor per physical mechanism. This is a direct refutation of the proposed Friis interpretation.

### 1.3 Dark counts also couple the factors

If signal arrivals of rate `lambda` and indistinguishable dark arrivals of rate `d` both trigger the same detector, the active-arrival rate is `r=lambda+d`. For information about `log lambda`, the long-window count information becomes

\[
\boxed{
J_{\log\lambda}
=
\left(\frac{\lambda}{\lambda+d}\right)^2
\frac{rb}{(1+rb)\{1+c^2(rb)^2\}}.
}
\]

There is a label-ambiguity multiplier, but the dead-time and jitter terms use the **total** load `rb`; changing `d` changes every factor. This is not a context-free dark-count noise figure.

### 1.4 Afterpulsing is worse for the product conjecture

Afterpulsing creates self-excitation and serial dependence rather than an iid additive cycle delay. The long-run Fisher rate depends on the conditional law of the memory state or, equivalently in suitable limits, a zero-frequency covariance/Poisson-equation object. Correlated stochastic processes have no generic Fisher subadditivity, superadditivity, or multiplicative decomposition; explicit counterexamples are known ([Radaelli et al. 2023](https://doi.org/10.1088/1367-2630/acd321)). Detector characterization studies also show higher-order afterpulse correlations rather than a single independent loss stage ([Wayne, Bienfang & Polyakov 2017](https://doi.org/10.1364/OE.25.020352)).

**C1 final:** retain only the operator chain rule. Retire the scalar “one factor per mechanism” program.

---

## 2. C2: a real theorem exists, but `C_u` is not the alignment functional

Let a finite-dimensional control vector be `z`, with concave local information objective

\[
f(z)=\Phi\{V(z)\},
\]

and convex feasible set `C`. A rank-`k` controller restricts `z` to an affine control subspace

\[
\mathcal M_k=z_0+\operatorname{range}(U_k).
\]

Let `z_k` maximize `f` over `C cap M_k`. The invariant diagnostic is not a scene-only scalar. It is the distance of the full gradient from the full feasible normal cone:

\[
\boxed{
\varepsilon_k
=
\operatorname{dist}_{*}
\bigl(\nabla f(z_k),N_C(z_k)\bigr).
}
\]

For differentiable concave `f`,

\[
\varepsilon_k=0
\quad\Longleftrightarrow\quad
z_k\text{ is also globally optimal over }C.
\]

If `f` is `m`-strongly concave on the relevant feasible segment, then

\[
\boxed{
0\le f(z^*)-f(z_k)
\le
\frac{\varepsilon_k^2}{2m}.
}
\]

Proof: choose the closest `n in N_C(z_k)`, write `r=grad f(z_k)-n`, apply strong concavity, use `n^T(z-z_k)<=0` for all feasible `z`, and maximize `r^T d-(m/2)||d||^2`.

This is the clean rank-control theorem. The global source multiplier is the rank-one manifold `rho_i=kappa u_i`; its exact scalar KKT condition

\[
\ell_i J'(\kappa u_i)=\beta
\]

is simply the condition `epsilon_1=0` for the load-allocation slice.

### Why `C_u` cannot be `A_1`

Bucket contrast contains only the dispersion of brightness `u_i`. It omits the information directions `q_i`, the current inverse information matrix `V^{-1}`, the kernel derivative `J'`, and active constraints. Two pattern families can have identical brightness vectors and identical `C_u`, while one has rank-one/collinear `q_i` and the other spans the task space orthogonally; their D-information and rank-one-control headroom are radically different.

Therefore:

- `C_u` is a useful first-order **photon-visibility proxy**;
- the correct alignment diagnostic is a projected KKT residual involving `ell_i=q_i^T V^{-1}q_i`;
- the contour failure is consistent with rank-one misalignment, but is not “exactly” an `A_1` theorem from the observed correlation `r=0.64`.

This territory is adjacent to classical general-equivalence/KKT optimal-design theory ([Kiefer 1974](https://doi.org/10.1214/aos/1176342810); [White 1973](https://doi.org/10.1093/biomet/60.2.345); [Molchanov & Zuyev 2004](https://doi.org/10.1051/ps:2003016)). The dead-time kernel and physical control manifold are the specialized content, not the existence of the residual condition.

---

## 3. C3: useful picture, false as an iff obstruction theorem

Write `D(z)` for cumulative object-plane dose and `V(z)` for the information operator. At a current design, a control direction `h` is first-order feasible and useful when

\[
D'(z)h\in T_{\mathcal D}(D(z)),
\qquad
\langle\nabla_V\Phi,V'(z)h\rangle>0,
\]

together with the other active constraints.

The proposed statement “a control survives iff it acts trivially where the binding constraint lives” is too strong.

A two-direction counterexample is immediate:

\[
D'h_1=v,
\qquad
D'h_2=-v,
\]

while both directions improve different information modes. Neither control is dose-trivial, but `h_1+h_2` is dose-neutral and beneficial. A direction can also reduce an active dose inequality (`D'h<0`) while improving information. Conversely, a dose-null direction can be useless if `V'h` lies in an already saturated information subspace.

The correct statement is only the normal-cone/KKT condition already contained in C2/C4. “Object plane versus measurement spectrum” remains a helpful explanatory drawing, not a separate theorem.

The repository’s own evidence makes this distinction unavoidable:

- the frozen one-dimensional mixture path collapsed;
- dose-only feasible directions retained headroom;
- full-stack feasible directions also retained headroom.

That is a statement about different **control tangent sets**, not a universal commuting-plane obstruction.

---

# Q2–Q3 — The one object that actually compresses the system

## 4. The controlled conditional-score Gramian

I introduce one fourth object because it strictly dominates C1–C3. The name is descriptive, not a priority claim.

Let

\[
\mathcal E_u=\{Q_{\theta,u}:\theta\in\Theta\},
\qquad
Q_{\theta,u}(dy)=\int K_u(dy\mid x)P_\theta(dx),
\]

be a **controlled statistical experiment**: `X~P_theta` is the complete latent optical/detector record, and `K_u` is the measurement/detector channel under control `u` (pattern, gain, detector mode, etc.).

Let the complete score operator be

\[
\mathsf S_\theta h
=
h^T\nabla_\theta\log p_\theta(X),
\qquad
h\in T_\theta\Theta,
\]

and let

\[
\mathsf C_u f(Y_u)=\mathbb E[f(X)\mid Y_u]
\]

be conditional expectation through the controlled observation channel.

Define the local information object

\[
\boxed{
\mathcal G_u
=
\mathsf S_\theta^*\mathsf C_u^*\mathsf C_u\mathsf S_\theta.
}
\]

Equivalently, for tangent directions `h,k`,

\[
\langle h,\mathcal G_u k\rangle
=
\mathbb E\!\left[
\mathbb E(h^TS\mid Y_u)
\mathbb E(k^TS\mid Y_u)
\right].
\]

This is simply the Fisher information matrix of the observed controlled experiment, written so that hidden-state loss and control geometry are visible in the same object.

## 5. Projection–design theorem

Under the usual dominated differentiability and integrability assumptions:

### (A) Missing information is Pythagoras in score space

\[
\boxed{
\mathcal I_{\rm complete}-\mathcal G_u
=
\mathbb E[\operatorname{Cov}(S\mid Y_u)]
\succeq0.
}
\]

For the dead-time log-rate score `S=N-lambda L`, this gives E1 immediately:

\[
I_N=\mathbb E N-\lambda^2\mathbb E\operatorname{Var}(L\mid N).
\]

### (B) Detector mechanisms compose as contractions, not independent scalar factors

For nested channels, the conditional-score operators compose. Fisher information is Loewner-monotone under further garbling. C1’s product law occurs only when the active score subspace is one-dimensional and invariant under each contraction.

### (C) Independent controlled exposures add Gramians

For conditionally independent exposures `u_1,...,u_M`,

\[
V=V_0+\sum_{i=1}^M\mathcal G_{u_i}.
\]

For an approximate design measure `xi`,

\[
V(\xi)=V_0+M\int\mathcal G_u\,d\xi(u).
\]

For local D-optimality,

\[
f(\xi)=\log\det V(\xi),
\]

and the directional sensitivity is

\[
\boxed{
d(u;\xi)=M\operatorname{tr}\{V(\xi)^{-1}\mathcal G_u\}.}
\]

Constraint prices give the Kiefer–Wolfowitz/KKT condition

\[
d(u;\xi)-\lambda^Tc(u)\le\eta,
\]

with equality on support atoms. E3 is the load-coordinate derivative of this one equation.

### (D) Control rank is a restriction of the same variational object

A rank-`k` controller restricts the admissible `u` or design coordinates to a `k`-dimensional manifold. Full sufficiency is exactly the vanishing of the full projected KKT residual; strong concavity yields the residual-squared gap bound above. E4 is a zero-feasible-ascent result on one restricted path. E5 is a constructive nonzero-ascent result in a larger tangent cone. E6 is the rank-one case.

## 6. Proof strategy

1. **Score projection:** Fisher’s identity gives the observed score as `E[S|Y]`; the law of total covariance gives part (A).
2. **Channel composition:** apply the tower property to nested observations; use contraction of conditional expectation in `L^2`.
3. **Design additivity:** conditional independence makes log likelihoods and Fisher Gramians additive.
4. **Variational derivative:** differentiate `Phi(V)` through the affine information map; for `log det`, `D Phi[V](H)=tr(V^{-1}H)`.
5. **Rank-gap bound:** combine the normal-cone optimality condition with strong concavity.

## 7. Corollary map to E1–E6

| frozen result | status under the object |
|---|---|
| **E1** | exact matrix projection identity; the live-time formula is its dead-time specialization |
| **E2** | scalar long-window evaluation of one score-transfer direction; the denominator product is a renewal mean/variance simplification, not a general cascade law |
| **E3** | KKT derivative of `Phi(V(xi))` along load coordinates |
| **E4** | restricted-path tangent-cone obstruction; not a collapse theorem for the full feasible class |
| **E5** | constructive evidence that the full projected KKT residual is nonzero for SCAT32 in those local cells |
| **E6** | rank-one control case; successes suggest small residual on many scenes, failures suggest misalignment. The measured PSNR/dB values remain empirical and are **not** Fisher-theorem corollaries |

## 8. First counterexample risks

This unification has hard boundaries.

1. **Cross-pattern memory.** If detector state carries across pattern changes, exposure information is not additive. The latent state must be included in one controlled Markov experiment; the design-measure concavity can disappear.
2. **Nonconcave load branch.** Beyond an information ridge, `J(rho)` is not globally concave/increasing. Water-filling statements require a declared branch or a global nonconvex analysis.
3. **Image-domain risk.** `log det` Fisher design is a local surrogate. It does not imply TV-regularized PSNR ordering. This prevents E6’s numerical image gains from becoming formal corollaries.
4. **Scene-only alignment.** No scalar depending only on `C_u` can certify rank-one sufficiency; the task subspace and current information matrix are indispensable.

---

# Q4 — Honest depth ruling

The six results are **already near their natural mathematical depth**. There is no credible new scalar principle from which all six numerical and image-domain statements fall out.

The controlled conditional-score Gramian provides a clean conceptual compression, but the abstraction itself is established mathematics:

- incomplete-data/observed information: [Louis 1982](https://doi.org/10.1111/j.2517-6161.1982.tb01203.x);
- Fisher geometry under sufficient statistics and Markov morphisms: [Ay et al. 2015](https://doi.org/10.1007/s00440-014-0574-8), building on Chentsov;
- statistical-experiment ordering: Blackwell and Le Cam;
- design measures and equivalence/KKT theory: [Kiefer 1974](https://doi.org/10.1214/aos/1176342810), [White 1973](https://doi.org/10.1093/biomet/60.2.345), [Molchanov & Zuyev 2004](https://doi.org/10.1051/ps:2003016);
- task-based imaging information operators and subspaces: Barrett & Myers, *Foundations of Image Science*;
- brightness-dependent multiplex advantage/disadvantage: [Oliver & Pike 1974](https://doi.org/10.1364/AO.13.000158) and, for modern single-pixel imaging under shot noise, [Scotté, Galland & Rigneault 2023](https://doi.org/10.1088/2515-7647/acc70b).

Therefore **do not write a paper whose novelty is the object itself**.

The one focused derivation attempt still worth making is:

> **A rank-constrained control theorem for hidden-state photon-counting experiments:** exact/full KKT sufficiency, computable projected-residual upper bounds, and lower bounds showing when rank-one global power control cannot approximate unrestricted load/geometry control under nonnegative illumination and physical dose/peak constraints.

That is C2, but formulated inside C4. It could genuinely bridge the two papers. C1 should not be pursued as a universal law; C3 should remain explanatory geometry.

### Retirement ledger

The compression requirement is satisfied as follows:

- the **conditional-score Gramian** retires the separate “missing-information principle” and the proposed scalar Friis algebra;
- its **channel-composition law** retires mechanism-by-mechanism factors and ad hoc detector-noise taxonomy;
- the **projected KKT residual**, derived from the same Gramian, retires the proposed `A_k` functional and the conjugate-plane obstruction metaphor;
- `C_u` returns to its proper role as a visibility proxy, not a universal alignment statistic.

---

# Q5 — Physics reach probes

## (a) Paralyzable and correlated-hold detectors

The controlled-score object survives; the simple live-time identity does not. In a paralyzable detector, arrivals during a blocked interval alter future state, so the complete score is no longer `N-lambda L` with a state path independent of hidden arrivals. One must take the full latent arrival/state trajectory as `X`; the observed information is still the conditional-score Gramian and still obeys the exact missing-information matrix identity. Long-window evaluation becomes a Markov-renewal or Markov-additive-functional problem, typically involving a Poisson equation or the zero-frequency covariance of the score. For finite-order correlated processes, asymptotic Fisher rates remain linear under suitable conditions, but correlations have no universal favorable or unfavorable algebra ([Radaelli et al. 2023](https://doi.org/10.1088/1367-2630/acd321)). Richer event records can preserve information discarded by integrated counts, as the dead-time event-detection theory of [Jorgensen & Johnson 2026](https://arxiv.org/abs/2605.23210) makes explicit. The exponent pair and two-factor renewal law must therefore be rederived, not transferred.

## (b) Sub-Poissonian or quantum illumination meeting dead time

For a fixed photodetection measurement, the same classical controlled experiment applies after the source state is mapped to a point process; the complete score simply comes from a non-Poisson latent arrival law. If the optical measurement itself is optimized, the upstream object is a quantum statistical experiment and the detector is a memoryful quantum instrument followed by a classical record; quantum Fisher information supplies an upper metric, while the observed click record is governed by the classical conditional-score Gramian. Dead time does not merely replace the Poisson variance by a Fano factor: it reshapes temporal correlations and can bias measured `g^(2)` and Mandel-Q even in two-detector configurations ([Ann et al. 2015](https://arxiv.org/abs/1503.06263)). Thus the operator unification survives, but E2’s scalar kernel and any product rule do not.

## (c) Detector arrays sharing one source

For `m` conditionally independent detector elements, the observed Gramian is a sum of block/element Gramians. Under ideal equal Poisson thinning, identical detectors, and no crosstalk, a total load `rho` split equally gives the scalar rate

\[
\boxed{J_m(\rho,c)=m\,J(\rho/m,c).}
\]

As `m to infinity`, `J_m to rho`, so spatial parallelism removes the single-element dead-time bottleneck and recovers the low-load Poisson information rate. A shared source constraint still couples the control allocation, and source fluctuations, optical crosstalk, shared quenching, or common electronics introduce off-diagonal score covariance. The global source multiplier is then a rank-one common-mode control, while per-element gains form higher-rank controls—exactly the rank hierarchy described above. Existing SPAD-array work already analyzes dead-time-dependent count distributions and information rates ([Sarbazi, Safari & Haas 2020](https://doi.org/10.1109/TCOMM.2020.2993374)); the residual opportunity is a task-based controlled-Gramian allocation theorem, not the fact that arrays parallelize dead time.

---

# Final blunt conclusion

There is a unification, but it is **architectural rather than miraculous**:

\[
\boxed{
\text{latent score}
\xrightarrow{\text{controlled observation channel}}
\text{conditional-score Gramian}
\xrightarrow{\text{variational design}}
\text{KKT/control-rank residual}.
}
\]

E1 and E2 evaluate the first arrow. E3–E5 optimize the second. E6 tests the rank-one restriction empirically.

That is the cleanest compression available. It is also mostly known mathematics. The next work is justified only if it proves a nontrivial **rank-control approximation or impossibility theorem** for physically constrained dead-time counting channels. Any broader “Friis law of detector mechanisms” or “conjugate-plane obstruction principle” would currently be numerology dressed as unification.