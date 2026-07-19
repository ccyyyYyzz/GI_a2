# R10 — Final ruling: dead-time-aware optimal experiment design

**Scope lock:** this ruling governs the follow-up method line only. It does not reopen or modify the R9 scope lock on the Optics Express manuscript or reinterpret either frozen Study-1/Study-2 verdict.

## Executive ruling

| Question | Frozen ruling |
|---|---|
| Q1 — equivalence theorem | The approximate-design program is valid after one essential correction: the information atom for image parameters is obtained by a chain rule from the exact Fisher information for **log rate**, not by using `g ∝ Φa` alone. A generalized Kiefer–Wolfowitz/KKT certificate applies if the coarse estimate is frozen, the atom set is compact, and the information is nonsingular on a declared estimable subspace. Under pointwise equal load, the renewal kernel factors out and does **not** change the optimal geometry. Genuine dead-time-aware design therefore requires controlled pattern-to-pattern load variation or a load-dependent feasible set. |
| Q2 — rounding | For an unregularized rank-`r` rank-one design, i.i.d. rounding has the exact falling-factorial determinant guarantee `((M)_r/M^r)^(1/r)` in D-efficiency; this is weak when `M≈r`. Use efficient weight rounding followed by exact exchange/reweighting, and report the achieved continuous-to-exact D-efficiency and spectral guard. Joint TV reconstruction does not invalidate information additivity under independent active-start exposures. |
| Q3 — inexact oracle | Guarantees attach only to an oracle that approximates the **Frank–Wolfe gap**, not the raw atom score. A multiplicative gap oracle with quality `δ` degrades the standard rate to `O(C_F/(δ²t))`. A learned proposer is computation-only unless its candidate is certified against an exact or upper-bounded dictionary oracle; the frozen confirmatory implementation must finish with an exact dictionary scan. |
| Q4 — objective and guards | Freeze local sequential **regularized D-optimality** as the primary design surrogate, with pre-scan information included as `V0`; do not claim it optimizes TV-PSNR. Use A-risk and relative spectral floors as mandatory guards. Exclude `p=2`; cap weight dispersion, peak irradiance, per-pixel dose imbalance, load tails, FW certificate slack, and exact-rounding loss. `Γ` and `C_u` remain descriptive only. |
| Q5 — adaptivity theorem | A clean `L`-factor advantage is provable only for a restricted union-of-scales model with minimum contrast and scale-separated dictionaries. The allocation lower bound and two-stage identification argument are standard; the new work is verifying those assumptions under nonnegative equal-resource patterns and the exact renewal channel. It is a follow-up theorem target, not a universal “adaptivity is necessary” claim. |
| Q6 — novelty | No closer prior was found that combines exact finite-window nonparalyzable renewal information, illumination OED, physical equal-resource constraints, and a generalized-equivalence certificate. The safe claim is this combination. Do not claim first OED for photon counting, first learned/adaptive SPI, first dead-time optimum, or a universal rescuability certificate. |
| Q7 — campaign | Use fresh confirmatory scenes, a common 52-pattern pre-scan charged to every arm, 972 designed/main patterns, a variable-load dead-time OED primary plus an equal-load OED ablation, fixed baselines, the existing Q90 speed gate, and a separate design-superiority secondary. The certificate is explicitly relative to a frozen structured dictionary. |

---

# Q1 — General equivalence theorem under the exact dead-time kernel

## 1. Correct information atom

The proposed formula needs a chain-rule correction. Let the coarse pre-scan reconstruction be frozen at `x̂`, and parameterize the locally estimable image perturbation as

\[
 x(\theta)=\hat x+B\theta,
 \qquad B\in\mathbb R^{n\times r},
\]

where the columns of `B` span the declared `r`-dimensional estimable/task subspace. For a physical nonnegative illumination atom `a`,

\[
 \lambda(a)=\Phi a^{\mathsf T}\hat x+d,
 \qquad \rho(a)=\tau\lambda(a).
\]

Let `J(ρ,ν)` be the exact **per-slot** Fisher information of the integrated renewal count for `log λ`, as derived in ROUND63. One exposure contains

\[
 I_{\log\lambda}(a)=\nu J(\rho(a),\nu).
\]

The score direction for `θ` is

\[
 q(a)=\frac{\partial\log\lambda(a)}{\partial\theta}
 =\frac{\Phi B^{\mathsf T}a}{\lambda(a)}.
\]

Therefore the correct rank-one local information atom is

\[
 \boxed{
 H(a;\hat x)=\nu J(\rho(a),\nu)\,q(a)q(a)^{\mathsf T}.
 }
\]

Equivalently,

\[
 H(a;\hat x)=
 \nu J(\rho(a),\nu)
 \frac{\Phi^2}{\lambda(a)^2}
 B^{\mathsf T}aa^{\mathsf T}B.
\]

Thus `g_i ∝ Φa_i` is incomplete unless the missing factor `sqrt(νJ)/λ` is absorbed into the definition of `g_i`.

If source flux, dark rate, or another calibration variable is uncertain, either augment `θ` with those nuisance parameters or form the effective image information by the Schur complement after assembling the joint information matrix. A D-optimal design that silently treats uncertain flux/dark as known can select patterns that are excellent only because they confound radiometry with structure.

## 2. Equal-load cancellation is a theorem, not a detail

If every feasible atom obeys the pointwise equal-load constraint

\[
 \lambda(a)=\bar\lambda
\]

at `x̂`, then `J(ρ(a),ν)`, `λ(a)^{-2}`, and `Φ` are constant over the atom set. Hence

\[
 H(a)=c(\bar\rho,\nu)B^{\mathsf T}aa^{\mathsf T}B,
\]

and the scalar `c` factors out of every homogeneous optimality criterion. The D-optimal geometry is then identical to ordinary local linear D-optimality. In that regime the method may be **dead-time evaluated**, but the exact renewal kernel does not alter the selected geometry.

Accordingly, freeze two distinct variants:

1. **OED-EQLOAD** — pointwise equal load; geometry-only ablation.
2. **OED-DT** — the gate-carrying method; pattern amplitude/load may vary within frozen bounds, while average load, total dose, peak irradiance, and per-pixel cumulative dose remain constrained. Here `J(ρ(a),ν)` changes the reduced sensitivity and genuinely affects the design.

Do not call an equal-load-only optimizer “dead-time-aware OED” without this caveat.

## 3. Design-measure formulation

Let `𝒜` be the compact feasible atom set after pointwise constraints, and let `ξ` be a probability measure over `𝒜`. Define

\[
 V(\xi)=V_0+\int_{\mathcal A}H(a)\,d\xi(a),
\]

where `V0 ≽ 0` is frozen pre-scan information plus a numerical ridge if needed. The local D-criterion is

\[
 \phi_D(\xi)=\log\det V(\xi).
\]

Allow linear design-measure constraints

\[
 \int e_j(a)d\xi(a)=b_j,
 \qquad
 \int c_k(a)d\xi(a)\le u_k,
\]

for average load, cumulative dose, family allocation, or other resources.

## 4. Generalized equivalence theorem

### Theorem R10.1 — constrained local renewal D-optimality

Assume:

1. `𝒜` is compact;
2. `H(a)`, `e_j(a)`, and `c_k(a)` are continuous;
3. the feasible measure set is nonempty and satisfies the relevant constraint qualification/Slater condition;
4. there exists a feasible `ξ` with `V(ξ) ≻ 0` on the declared `r`-dimensional subspace;
5. `x̂`, `B`, `V0`, detector parameters, and the atom dictionary are fixed during the optimization.

Then a feasible design `ξ*` maximizes `log det V(ξ)` if and only if there exist free multipliers `α_j` for the equality constraints and nonnegative multipliers `β_k` for the inequalities, satisfying complementary slackness, such that the reduced sensitivity

\[
 s_{\alpha,\beta}(a)
 =\operatorname{tr}\!\left[V(\xi^*)^{-1}H(a)\right]
 -\sum_j\alpha_j e_j(a)-\sum_k\beta_k c_k(a)
\]

obeys

\[
 \boxed{
 s_{\alpha,\beta}(a)\le\eta^*
 \quad\text{for all }a\in\mathcal A,
 }
\]

with equality for `ξ*`-almost every support atom. The constant is

\[
 \eta^*=
 \operatorname{tr}\!\left[V(\xi^*)^{-1}
 \{V(\xi^*)-V_0\}\right]
 -\sum_j\alpha_j b_j
 -\sum_k\beta_k\bar c_k^*,
\]

where `c̄_k* = ∫c_k dξ*`; active inequalities have `c̄_k*=u_k`.

In the unconstrained classical case with `V0=0`, this reduces to

\[
 \boxed{
 \max_{a\in\mathcal A}
 \operatorname{tr}[V(\xi^*)^{-1}H(a)]=r,
 }
\]

with equality on the support. If `H(a)=h(a)h(a)^T`, the sensitivity is simply

\[
 h(a)^T V(\xi^*)^{-1}h(a).
\]

With `V0>0`, the threshold is not `r`; it is

\[
 \operatorname{tr}[V^{-1}(V-V_0)]<r.
\]

### Proof sketch

The map `ξ ↦ V(ξ)` is affine and `log det` is concave on positive-definite matrices. Its directional derivative from `ξ` toward a feasible measure `ζ` is

\[
 D\phi_D(\xi;\zeta-\xi)
 =\int\operatorname{tr}[V(\xi)^{-1}H(a)]
 \,d(\zeta-\xi)(a).
\]

The measure-space KKT condition says this derivative, after subtracting the linear resource multipliers, is nonpositive in every feasible direction. Concentrating an infinitesimal mass at any atom gives the pointwise inequality; integrating the support equality against `ξ*` yields `η*`. Conversely, the pointwise inequality integrated against any feasible `ζ`, together with complementarity, gives a nonpositive directional derivative. Concavity then makes the condition sufficient for global optimality.

This is the Kiefer–Wolfowitz logic adapted to a nonlinear local model and linear resource restrictions. Classical sources are Kiefer’s optimum-design framework and the nonlinear extension of the general equivalence theorem; constrained designs require the KKT/reduced-sensitivity form rather than the bare `max sensitivity = r` statement. See [Kiefer (1959)](https://doi.org/10.1111/j.2517-6161.1959.tb00338.x) and [White (1973)](https://doi.org/10.1093/biomet/60.2.345).

## 5. Hidden failure modes

- **Noncompactness:** `a ≥ 0` alone is an unbounded cone. A finite peak-irradiance bound and dose/load constraints are mandatory for existence.
- **Rank deficiency:** if the dictionary does not span the chosen task subspace, unregularized `log det` is `−∞`. Restrict to the estimable subspace or include frozen `V0 ≻ 0`; do not hide a singular design with an arbitrary large ridge.
- **Changing `x̂` inside one convex solve:** updating `x̂` changes every atom `H(a)` and destroys the one-shot concavity statement. Adaptivity must be stagewise: freeze `x̂`, solve, acquire, then optionally relinearize in a new stage.
- **Nonnegativity:** the nonnegative atom constraint causes no theorem failure once the atom set is compact.
- **TV-MAP:** the theorem certifies the local information criterion, not unbiasedness or TV-reconstruction risk. The campaign must state D-optimality as a surrogate and validate the frozen image endpoint empirically.

---

# Q2 — Continuous-to-exact rounding

## 1. Standard determinant guarantee

For the unregularized rank-one case in an `r`-dimensional estimable subspace, write

\[
 M(\xi^*)=\mathbb E_{a\sim\xi^*}[h(a)h(a)^T].
\]

Draw `M` atoms independently from `ξ*` and form

\[
 \widehat M_M=\frac1M\sum_{m=1}^M h_mh_m^T.
\]

Cauchy–Binet gives the exact identity

\[
 \mathbb E\det(\widehat M_M)
 =\frac{(M)_r}{M^r}\det M(\xi^*),
\]

where `(M)_r=M(M−1)…(M−r+1)`. Therefore at least one exact multiset obeys

\[
 \boxed{
 \operatorname{Eff}_D(\widehat M_M;M^*)
 \ge
 \left(\frac{(M)_r}{M^r}\right)^{1/r}.
 }
\]

This is the clean universal multiplicative guarantee. It is also a warning: when `M=r=1024`, the bound is only about `0.369`. No theorem can turn a nearly square, unrestricted full-image design into a uniformly high-efficiency exact design without additional structure, prior information, or more measurements.

Volume-sampling and proportional-volume-sampling work provides related exact-design approximation guarantees; see [Nikolov, Singh & Tantipongpipat](https://doi.org/10.1287/moor.2021.1129). Classical approximate-to-exact apportionment is due to [Pukelsheim & Rieder](https://doi.org/10.1093/biomet/79.4.763).

## 2. High-probability spectral bound

If the leverage is bounded by

\[
 h(a)^TM(\xi^*)^{-1}h(a)\le L
\]

for every atom, matrix Chernoff yields, for `0<ε<1`,

\[
 P\{\widehat M_M\not\succeq(1-\epsilon)M(\xi^*)\}
 \le r\exp\!\left(-\frac{M\epsilon^2}{2L}\right).
\]

Thus

\[
 M\ge\frac{2L}{\epsilon^2}\log\frac r\delta
\]

suffices for the lower spectral approximation with probability at least `1−δ`; then

\[
 \det(\widehat M_M)^{1/r}
 \ge(1-\epsilon)\det(M^*)^{1/r}.
\]

At an unconstrained D-optimum `L=r`, so this bound is also weak when `M≈r`.

## 3. Frozen rounding implementation

Use the following exact-design pipeline:

1. fully corrective Frank–Wolfe/column generation to obtain support atoms and continuous weights;
2. Pukelsheim–Rieder efficient rounding of `M_main ξ*` to integer replicate counts;
3. exact reoptimization of the support weights subject to integer total count;
4. deterministic one-swap/one-replicate exchange until no dictionary atom improves the exact objective;
5. compute the final exact-design objective, spectral guard, A-risk guard, and continuous upper bound.

The confirmatory design is admissible only if

\[
 \boxed{\operatorname{Eff}_D^{\rm exact}\ge0.95}
\]

relative to the continuous optimum **and** all Q4 guards pass. If not, the method returns `ROUNDING_FAIL`; it is not repaired after seeing image endpoints.

## 4. Why joint reconstruction does not change the bound

The information matrices add because the active-start exposures are conditionally independent given the scene and calibration parameters. The fact that all rows enter one joint TV reconstruction is irrelevant to Fisher-information additivity. The standard exact-design logic therefore applies to repeated and nonrepeated rows alike.

It does **not** apply unchanged to continuous acquisition with detector-state carryover or unpurged afterpulse memory, because then the likelihood does not factor by pattern. Keep the OED campaign on independent active-start exposures.

---

# Q3 — Frank–Wolfe with an inexact inner atom oracle

## 1. Approximate the gap, not the raw score

For concave maximization of `F(ξ)` over a compact convex design set, define the exact Frank–Wolfe gap

\[
 g_t=\max_{s\in\mathcal D}
 \langle\nabla F(\xi_t),s-\xi_t\rangle.
\]

The oracle must return `s̃_t` satisfying

\[
 \boxed{
 \tilde g_t=
 \langle\nabla F(\xi_t),\tilde s_t-\xi_t\rangle
 \ge\delta g_t,
 \quad 0<\delta\le1.
 }
\]

A multiplicative approximation to `⟨∇F,s⟩` itself is not shift invariant and is not a valid substitute.

## 2. Convergence degradation

Let `C_F` be the Frank–Wolfe curvature constant for `−F`. With line search, the standard curvature recursion is

\[
 h_{t+1}
 \le(1-\delta\gamma_t)h_t+rac{C_F}{2}\gamma_t^2,
 \qquad h_t=F^*-F(\xi_t).
\]

Taking `γ_t=2/(δt+2)` gives, after the first step,

\[
 \boxed{
 F^*-F(\xi_t)
 \le
 \frac{2C_F}{\delta(\delta t+2)}
 =O\!\left(\frac{C_F}{\delta^2t}\right).
 }
\]

Thus a 90%-quality gap oracle (`δ=0.9`) costs at most the standard `1/δ²≈1.23` factor in the leading iteration bound; a weak learned oracle can be much more damaging.

For additive inexactness, Jaggi’s standard theorem is directly quotable: if the linear subproblem error at iteration `t` is at most a prescribed `O(C_F/(t+2))` term, the `O(1/t)` rate is retained with a multiplicative constant. See [Jaggi, ICML 2013](https://proceedings.mlr.press/v28/jaggi13.html).

## 3. Frozen learned-proposer rule

A learned proposer may:

- warm-start the atom search;
- nominate shape parameters;
- prioritize a subset of translates;
- expand the dictionary during **development**.

It does not inherit a guarantee merely because it performs well on images. For a theorem-backed run, one of the following must hold:

1. an exact structured-dictionary scan computes the true best reduced sensitivity and therefore certifies `δ`; or
2. a rigorous upper bound on the unknown best sensitivity is available, allowing a certified lower bound on `δ`.

Because the frozen dictionary consists of circulant translate families and all translations can be evaluated by FFT, the confirmatory solver must end every outer iteration with the exact dictionary oracle (`δ=1`). A learned proposer is therefore optional computation, not a scientific contribution and not part of the frozen guarantee.

If later work enlarges the atom set beyond exact enumeration, report the achieved certified `δ_t` sequence and use the minimum certified value in the convergence statement. No certificate against an unsearched continuous universe is allowed.

---

# Q4 — Endpoint-aligned objective and mandatory guards

## 1. Primary objective

Freeze **sequential pre-scan-regularized local D-optimality**:

\[
 \boxed{
 \max_{\xi}\;
 \log\det\left[V_0+\int H(a;\hat x)d\xi(a)\right].
 }
\]

Here `V0` is the exact local information contributed by the common pre-scan, plus only the minimal numerical ridge required for stable linear algebra. This choice is preferred because it is:

- concave in the design measure;
- certificateable by the generalized equivalence theorem;
- naturally sequential;
- sensitive to global conditioning, unlike `C_u` or `Γ`;
- free of a new learned image prior that could be tuned to the endpoint.

The paper must call it a **local information surrogate**. It is not mathematically equivalent to PSNR, Q90, or TV-MAP risk.

Do not freeze low-frequency weighted A-optimality: the failure classes include fine microtexture, so a low-frequency weight would encode the wrong task. Do not freeze a TV-Bayes-risk optimizer as the primary method: a nonsmooth, scene-dependent TV prior creates a circular design/reconstruction model and removes the clean certificate.

## 2. Endpoint-alignment audit

In addition to D-optimality, compute the unweighted local A-risk on the radiometric task subspace:

\[
 R_A(\xi)=\operatorname{tr}[V(\xi)^{-1}].
\]

It is a mandatory guard and reported secondary design metric, not the optimized objective. Final performance is still judged by the preregistered image endpoints.

## 3. Mandatory guards

The raw OED candidate is accepted only after all guards pass.

### G1 — exact physical information

Use the exact finite-window renewal `J(ρ,ν)` and the chain-rule atom `H(a)` from Q1. Include known nuisance parameters correctly; use a Schur complement if they are estimated.

### G2 — genuine dead-time dependence

The primary OED-DT design permits controlled variable load. Freeze:

- safe arm atom loads: `ρ(a) ∈ {0.025, 0.05, 0.10}`;
- high-flux arm atom loads: `ρ(a) ∈ {0.30, 0.60, 1.00}`;
- design-average load exactly `0.05` or `0.60`, respectively.

The equal-load OED-EQLOAD ablation fixes every atom at the mean load. It has no gate.

### G3 — peak irradiance

Normalize the physical atom so that

\[
 0\le a_j\le64,
\]

matching the Study-2 normalized `k=16` peak factor at `n=1024`. No optimized pattern may demand more peak irradiance than the established sparse benchmark.

### G4 — weight dispersion

For matched-intensity atoms, normalize positive weights to mean one on their support and enforce

\[
 \boxed{1/4\le w_j\le4.}
\]

The frozen dictionary contains powers `p∈{0,1}` only. `p=2` is excluded from confirmatory design because the dev evidence already demonstrated that increasing proxy contrast can worsen conditioning and endpoint performance.

### G5 — per-pixel cumulative dose

Across the full exact main design,

\[
 0.95\le
 \frac{M^{-1}\sum_i a_{ij}}
 {n^{-1}\sum_{\ell}M^{-1}\sum_i a_{i\ell}}
 \le1.05
 \quad\text{for every pixel }j.
\]

This prevents the optimizer from winning by repeatedly overexposing an easy region of `x̂`.

### G6 — relative spectral floor

Let `V_fix` be the information matrix of the globally frozen strongest fixed baseline under identical pre-scan and resource accounting. Require

\[
 \boxed{
 V_{\rm OED}\succeq0.5V_{\rm fix}
 }
\]

on the declared task subspace. If the raw optimum fails, form the largest convex mixture

\[
 \xi_{\alpha}=\alpha\xi_{\rm OED}+(1-\alpha)\xi_{\rm fix}
\]

that satisfies the floor. The mixing coefficient is determined from information matrices before any image endpoint is evaluated.

### G7 — A-risk floor

Require

\[
 \boxed{
 R_A(\xi_{\rm final})
 \le1.05R_A(\xi_{\rm fix}).
 }
\]

This is the direct guard against a high-determinant design with an unacceptable small-eigenvalue tail.

### G8 — continuous-design certificate

For the final continuous design, use the exact constrained reduced-sensitivity oracle and require

\[
 \boxed{
 g_{\rm FW}/r\le10^{-3}.
 }
\]

The certificate is explicitly relative to the frozen dictionary and resource constraints.

### G9 — exact-design loss

After rounding and exchange,

\[
 \boxed{
 \operatorname{Eff}_D^{\rm exact}\ge0.95.
 }
\]

### G10 — no proxy gates

`C_u`, `Γ`, pre-scan roughness, and predicted endpoint curves are descriptive. None may accept/reject an atom, a scene, or a confirmatory result except insofar as they enter the already frozen physical dictionary construction.

## 4. Computable “do not bother” verdict

Reject the proposed rule

> `max predicted Γ < 1 ⇒ certified do not bother`.

The `p=2` result has already shown that higher proxy contrast can coexist with worse endpoint performance. A certificate must use the optimized matrix objective and its dual upper bound.

Define the best possible local D-gain within the frozen dictionary using the FW upper bound:

\[
 U_D=\frac1r\left[
 \phi_D(\xi_t)+g_t-\phi_D(\xi_{\rm fix})
 \right].
\]

Also compute the best certified A-risk improvement available from the final candidate. The allowed verdict is:

```text
LOCAL_NO_MATERIAL_GAIN
```

only if

\[
 \exp(U_D)\le1.05
\]

and the certified A-risk improvement is below 5%. Its frozen wording is:

> Within the frozen atom dictionary, local linearization at the pre-scan estimate, and declared resource constraints, the continuous-design upper bound permits less than a 5% information-efficiency improvement over the fixed baseline.

It is not a universal impossibility result and must not be called “unrescuable scene.”

---

# Q5 — Adaptivity advantage in OED language

## 1. Rejected universal form

No theorem can say that adaptivity rescues every nonflat scene. Without a minimum contrast, `x=x0+εv` makes every finite-budget information measure vanish as `ε→0`. A full fixed raster dictionary also defeats any universal claim that every fixed ensemble has a nonflat blind direction.

The valid target is class-conditional.

## 2. Stylized theorem target

### Theorem target R10.2 — scale-adaptive allocation advantage

Let the scene class be a union of `L` scale classes

\[
 \mathcal X=\bigcup_{\ell=1}^L\mathcal X_\ell,
\]

with minimum structural amplitude `ε>0`. Let `𝒜_ℓ` be the scale-matched atom family. Suppose the exact-renewal local information obeys:

1. **matched information:** there is a design `ξ_ℓ*` supported on `𝒜_ℓ` whose per-measurement task information is at least `μ>0` for every scene in `𝒳_ℓ`;
2. **cross-scale leakage:** an atom assigned to another scale contributes at most `βμ` information to class `ℓ`, with `0≤β<1/L`;
3. **bounded atom information:** no atom contributes more than `κμ` to any class;
4. **coarse identifiability:** there is a frozen coarse action family for which the pairwise exact-renewal KL divergence between scale hypotheses is at least `d0>0` per coarse measurement.

Then:

### Fixed design lower bound

Any nonadaptive allocation of `M` measurements has some scale `ℓ` receiving at most `M/L` matched measurements. Its worst-scale information is bounded by

\[
 \boxed{
 I_{\rm fixed}^{\rm worst}
 \le\mu M\left(\frac{\kappa}{L}+\beta\right)
 }
\]

up to the precise task-normalization constants.

### Two-stage design

Use

\[
 m_0\ge\frac{\log L+\log(1/\delta)}{d_0}
\]

coarse measurements to select the scale, then allocate the remaining `M−m0` measurements according to `ξ_{\hat\ell}*`. Standard Chernoff/multi-hypothesis testing bounds give scale-identification error at most `δ`, and hence

\[
 \boxed{
 \mathbb E I_{\rm adaptive}
 \ge(1-\delta)\mu(M-m_0).
 }
\]

When `m0=o(M)` and `βL=o(1)`, the adaptive-to-fixed worst-case ratio is order `L`.

## 3. Proof sketch and provenance

- The fixed-design result is a pigeonhole/allocation argument plus the cross-scale leakage bound.
- The coarse-stage sample complexity is a standard multi-hypothesis Chernoff/KL argument; sequential experiment selection goes back to Chernoff’s work on sequential design.
- The adaptive sensing literature supplies both positive constructions and warnings that adaptivity is not universally superior; see [Arias-Castro, Candès & Davenport](https://arxiv.org/abs/1111.4646) and [Chernoff’s sequential design report](https://statistics.stanford.edu/technical-reports/sequential-design-experiments).
- The new scientific burden is to construct physically feasible nonnegative pattern families satisfying the matched/leakage assumptions under equal resource constraints, and to calculate `μ`, `β`, and `d0` from the exact renewal FI/KL rather than a Gaussian channel.

## 4. Dead-time-specific caveat

If every atom is pointwise equal-load, the scalar renewal kernel cancels from the geometry and the order-`L` advantage is a scale-matching result, not a dead-time-specific adaptivity theorem. To make dead time essential, at least one of the following must be shown:

- the optimal scale allocation changes with the load-dependent kernel;
- controlled variable-load actions improve the exact-renewal KL separation;
- the adaptive policy keeps informative patterns near the finite-window ridge while a fixed policy cannot under the same average resources.

This theorem is the anchor of the follow-up theory program, not a claim to place in the current OE paper.

---

# Q6 — Narrow prior-art claim

## 1. Crowded neighboring areas

The following are established and cannot support priority claims:

- classical and nonlinear optimal experimental design and equivalence certificates;
- Poisson/photon-counting Fisher-information design;
- Bayesian and goal-oriented OED in imaging inverse problems;
- learned sensing-matrix design for nonlinear inverse problems;
- adaptive or dictionary-trained illumination in SPI;
- optimal photon flux under SPAD pile-up/dead time;
- exact/asymptotic dead-time information analysis.

Representative neighboring sources include:

- [Feng et al., optimal illumination patterns for SPI](https://arxiv.org/abs/1806.01340);
- [Burger et al., sequentially optimized X-ray projections](https://doi.org/10.1088/1361-6420/ac01a4);
- [Ma, Xu & Maleki, sensing matrices for generalized nonlinear inverse problems](https://arxiv.org/abs/2111.03237);
- [Attia, Alexanderian & Saibaba, goal-oriented OED](https://doi.org/10.1088/1361-6420/aad210);
- [Jorgensen & Johnson, dead-time-constrained event information](https://arxiv.org/abs/2605.23210);
- recent Fisher-driven photon-efficient FLIM design under a Poisson model, [Sumaya-Martinez & Torres-Garcia](https://arxiv.org/abs/2601.00490).

## 2. Search verdict

No closer prior was found that combines all of the following:

1. illumination-pattern design for integrated-count single-pixel imaging;
2. atom information computed from the exact finite-window nonparalyzable renewal count law;
3. equal average detector load, equal dose/per-pixel dose, and peak-irradiance constraints;
4. generalized-equivalence/FW certification of the continuous design;
5. an exact-design rounding certificate;
6. a dictionary-relative upper-bound verdict on whether further local design gain is available.

This is a defensible narrow combination, not a priority proof.

## 3. Frozen novelty wording

> We formulate a local approximate-design framework for single-pixel illumination in which each pattern’s information is computed from the exact finite-window nonparalyzable renewal count law. The design is optimized under explicit detector-load, dose, per-pixel exposure, and peak-irradiance constraints, and is accompanied by generalized-equivalence and exact-rounding certificates. Existing photon-counting OED, adaptive SPI, and dead-time optimal-flux studies address neighboring components, but we found no prior framework combining these elements for integrated-count single-pixel illumination.

Do not use:

- “first OED for photon counting”;
- “first adaptive/learned SPI design”;
- “first information-optimal SPAD flux”;
- “first sparse illumination design”;
- “certified scene unrescuability”;
- “global optimum over all physically possible patterns.”

The certificate is always stated as relative to the frozen dictionary, local estimate, and resource model.

---

# Q7 — Method-campaign specification skeleton

## 1. Campaign identity and stop rule

This is a separate follow-up campaign. It must receive its own preregistration, seed namespace, manifests, expected-cell ledger, and immutable confirmatory freeze.

Development may tune implementation and select one global fixed comparator as described below. After confirmatory freeze, no scene, atom family, load grid, guard, endpoint, or gate may be changed. If OED fails, there is no third scene-class or endpoint redesign.

## 2. Geometry and detector constants

Freeze:

\[
 n=32^2=1024,
 \qquad M_{\rm total}=1024,
\]

\[
 M_{\rm pre}=52,
 \qquad M_{\rm main}=972.
\]

Thus the pre-scan is `52/1024 = 5.078125%` of physical patterns. Use:

- nonparalyzable active-start detector;
- `τ=50 ns`;
- `d=0` in the primary campaign;
- `σ_b=0`;
- `ρ̄∈{0.05,0.60}`;
- `ν∈{5,10,20,50,100,200,500,1000,2000}`;
- five measurement seeds per image/condition.

Total optical time is

\[
 T_{\rm opt}=M_{\rm total}\nu\tau.
\]

Design-computation wall time is reported separately and never hidden in optical time.

## 3. Fresh scenes

Use fresh confirmatory instances. Do not reuse any Study-2 confirmatory image as a decision unit.

Generate 24 images as six frozen families × four instances, retaining the same broad families so the new method is tested on the previously heterogeneous regime:

1. USAF/line-pair and chirp;
2. glyph/text;
3. spokes;
4. maze/barcode;
5. contours;
6. band-pass microtexture.

Freeze confirmatory seeds

\[
 s_{f,r}=633000+100f+r,
 \quad f=0,\ldots,5,
 \quad r=0,\ldots,3.
\]

Use six development-only images with seeds `632900+f`. Generator code, parameter manifests, images, and ROIs receive SHA256. No acceptance filter based on OED objective, `C_u`, `Γ`, PSNR, CNR, or visual quality is allowed.

The already frozen Study-2 cohort may be run as a no-gate external replication after the new verdict, but cannot enter confirmatory intervals.

## 4. Common pre-scan accounting

Every arm—OED and fixed—receives and pays for the **same** 52-pattern pre-scan at the same `(ρ̄,ν)` as its main acquisition.

- The pre-scan patterns are a frozen balanced multiscale set independent of the scene.
- All arms may include the pre-scan counts in their final reconstruction.
- Only adaptive/OED arms use the pre-scan to choose subsequent patterns.
- Fixed arms do not receive the 52 patterns back as extra main measurements.

This is the only fair accounting. Charging the pre-scan only to OED confounds method overhead with optical design; giving fixed methods 1024 main patterns while OED has 972 is permitted only because both still have exactly 1024 total physical exposures, including the common pre-scan.

## 5. Minimum frozen atom dictionary

The equivalence certificate is quoted only against this dictionary.

### Support families

Include all cyclic translates of:

1. scattered balanced supports, `k=16`;
2. scattered balanced supports, `k=32`;
3. solid `4×4` patch (`k=16`);
4. committed `Lblob6×6` compact 16-pixel support;
5. rectangles/bars `1×16`, `2×8`, `4×4`, `8×2`, `16×1`;
6. compact 32-pixel square/rectangle variants;
7. one annular/ring-like 16-pixel family for curved contours.

All translation families use the circulant/FFT evaluator.

### Intensity families

For each support, include:

- uniform weights `p=0`;
- matched-intensity weights `p=1` from the frozen pre-scan proxy, clipped and renormalized according to Q4-G4.

Exclude `p=2`.

### Load levels

Include the three Q4-G2 load levels for each atom shape. Average-load and per-pixel-dose constraints are imposed in the master design problem.

A learned proposer is not part of the minimum dictionary and carries no scientific claim.

## 6. Arms and baselines

### Gate-carrying arm

- **OED-DT:** variable-load, exact-renewal D-optimal design with all guards and exact dictionary certificate.

### Mandatory no-gate arms

1. **OED-EQLOAD:** same optimizer and dictionary, every atom fixed at `ρ̄`; isolates whether the renewal kernel changes design.
2. **SCAT16:** frozen scattered `k=16` benchmark.
3. **SCAT32:** best Study-2 ladder rung by fixed-dwell median.
4. **LBLOB16:** committed clustered compact-support rule.
5. **MATCH1:** frozen `p=1` matched-intensity heuristic on its declared support family.
6. **RASTER:** optional upper/control arm on the six mechanism representatives only.

Before confirmatory freeze, select **one global fixed comparator** `FIXED*` from `{SCAT32, LBLOB16, MATCH1}` using the six development images only. The selection score is the median fixed-dwell radiometric PSNR at `(ρ̄=0.60,ν=2000)`, averaged over development seeds. Tie within `0.05 dB` is resolved in order `SCAT32 → LBLOB16 → MATCH1`, favoring the least scene-dependent rule. `FIXED*` is then frozen for every confirmatory image; no per-scene oracle baseline is allowed.

All arms use the same RQL reconstruction, pooled initialization, analytic `λ_TV` rule with `c=0.50`, solver budget, and metrics.

## 7. Primary endpoint and gate

Use the existing Q90 time-to-quality machinery verbatim for OED-DT:

- seed-mean radiometric PSNR at each dwell;
- equal-weight PAVA isotonic curves;
- safe-range informativeness floor `0.50 dB`;
- Q90 target;
- frozen censoring taxonomy;
- 10,000 nested family-stratified bootstrap replicates.

Primary conditions remain:

\[
 \operatorname{median}S_{\rm gate}\ge3,
\]

\[
 \mathrm{LB}_{2.5\%}
 [\operatorname{median}S_{\rm gate}]>1,
\]

\[
 \#\{j:S_{{\rm gate},j}>1\}\ge18/24.
\]

Only OED-DT carries this gate. No baseline, OED-EQLOAD, or external cohort can rescue a failure.

## 8. Key secondaries

### Fixed-dwell safe-versus-fast gain

Keep the existing frozen endpoint at `ν=2000`:

\[
 \Delta Q_j^{\rm flux}
 =Q_{\rm OED-DT,0.60,j}-Q_{\rm OED-DT,0.05,j}.
\]

Use the existing practical bar:

- median `≥1.0 dB`;
- bootstrap lower bound `>0`;
- at least `18/24` positive.

It cannot rescue the primary.

### OED-versus-fixed design gain

At high flux and terminal dwell, define

\[
 \Delta Q_j^{\rm design}
 =Q_{\rm OED-DT,0.60,2000,j}
 -Q_{\rm FIXED^*,0.60,2000,j}.
\]

Preregister:

- median `≥0.50 dB`;
- family-stratified bootstrap lower bound `>0`;
- at least `18/24` positive images.

This is a key secondary and cannot rescue the Q90 primary. It is required because a method paper must show that OED improves over a strong globally frozen hand design, not merely that high flux helps under its own patterns.

## 9. Design diagnostics and fail states

Record for every designed cell:

- continuous `log det` objective;
- exact rounded objective and D-efficiency;
- A-risk;
- relative spectral floor;
- average and 5/50/95% pattern loads;
- peak irradiance;
- per-pixel cumulative dose spread;
- support size and replicate counts;
- exact constrained FW/KW gap;
- wall time;
- `C_u` and `Γ` as descriptive variables.

Frozen fail states are:

```text
DICTIONARY_RANK_FAIL
CONTINUOUS_CERTIFICATE_FAIL
SPECTRAL_GUARD_FAIL
A_RISK_GUARD_FAIL
DOSE_GUARD_FAIL
ROUNDING_FAIL
LOCAL_NO_MATERIAL_GAIN
```

A fail state is reported and the predeclared fallback is `FIXED*`; no atom family or threshold is changed after confirmatory freeze.

## 10. Development sequence before confirmatory freeze

Development may be used to:

1. verify the exact information atom against numerical Fisher derivatives;
2. validate the generalized-equivalence certificate on finite dictionaries;
3. verify FFT and brute-force atom scores agree;
4. validate continuous-to-exact rounding;
5. choose the single global `FIXED*` comparator by the frozen rule;
6. calibrate runtime and numerical tolerances;
7. test whether all guards are feasible.

Development may not change the scene families, endpoint form, primary thresholds, `p∈{0,1}`, peak cap, or basic load grid after the confirmatory freeze. Any pre-freeze amendment to an infeasible guard must be documented before any confirmatory seed is generated.

---

# Final method-line position

The unified scientific claim should be:

> Hand-designed sparse and matched patterns are atoms in a constrained experiment-design problem. Their usefulness is governed jointly by exact dead-time count information, spatial conditioning, and physical dose/peak constraints. A certified local OED can select and mix these atoms, while its generalized-equivalence gap quantifies the remaining design opportunity within the declared dictionary.

The method must **not** be sold as “maximize bucket contrast.” The dev results already falsify that proxy as a sufficient objective. The contribution is the replacement of one-dimensional contrast heuristics by a physically constrained information matrix with a certificate, plus a clean empirical test of whether that local information surrogate translates into the frozen image endpoint.

## Core references

- [Kiefer, Optimum Experimental Designs (1959)](https://doi.org/10.1111/j.2517-6161.1959.tb00338.x)
- [White, General Equivalence Theorem for nonlinear models (1973)](https://doi.org/10.1093/biomet/60.2.345)
- [Pukelsheim & Rieder, Efficient rounding of approximate designs (1992)](https://doi.org/10.1093/biomet/79.4.763)
- [Jaggi, Frank–Wolfe with approximate subproblems (2013)](https://proceedings.mlr.press/v28/jaggi13.html)
- [Sagnol & Harman, exact D-optimal design by mixed-integer conic optimization](https://arxiv.org/abs/1307.4953)
- [Nikolov, Singh & Tantipongpipat, volume-sampling design approximations](https://doi.org/10.1287/moor.2021.1129)
- [Arias-Castro, Candès & Davenport, limits of adaptive sensing](https://arxiv.org/abs/1111.4646)
- [Feng et al., adaptive illumination-pattern design in SPI](https://arxiv.org/abs/1806.01340)
- [Ma, Xu & Maleki, sensing-matrix design for generalized nonlinear inverse problems](https://arxiv.org/abs/2111.03237)
- [Jorgensen & Johnson, dead-time-constrained event information](https://arxiv.org/abs/2605.23210)
