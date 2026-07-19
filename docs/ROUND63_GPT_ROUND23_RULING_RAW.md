# R23 ruling (GitHub issue #15, raw)

Title: R23 ruling: rank-constrained control theorem for hidden-state photon counting
Posted: 2026-07-19T22:25:28Z

---

# R23 — Build the bridge: rank-constrained control for hidden-state photon-counting experiments

**Reference commit:** `f7723f6`  
**Source question:** `docs/ROUND63_GPT_ROUND23_QUESTION.md`

**Scope:** pure theory exploration for the next work. Both existing papers and all campaigns remain frozen. No endpoint, arm, or empirical claim is changed here.

> **R23 VERDICT.** A rigorous bridge exists, but it is narrower than the motivating slogan. On the monotone-concave branch of the dead-time kernel, rank-constrained illumination is a convex restricted-design problem. Its adequacy is governed by the **normal-cone marginal-information residual**, not by bucket contrast alone. This residual gives:
>
> 1. an exact KKT sufficiency theorem for every fixed control rank;
> 2. a computable strong-concavity upper bound on the unrestricted gap;
> 3. a finite-transfer lower bound and an explicit bright/uninformative versus dim/informative witness showing an `Omega(1)` rank-one gap;
> 4. an exact cluster-symmetry theorem and a local singular-value rank-decay theorem;
> 5. a first-order bridge in which bucket contrast appears only after leverage homogeneity is imposed.
>
> The genuinely new-looking specialization is the **dead-time elasticity alignment functional**
>
> \[
> \mathfrak M_1^2
> \simeq
> (\alpha-\alpha_J^*)^2 C_u^2+\sigma_\zeta^2,
> \qquad
> \alpha_J^*(\bar\rho)=-\frac{\bar\rho J''(\bar\rho)}{J'(\bar\rho)},
> \]
>
> not the generic KKT theorem. The latter is classical optimal-design/convex-analysis machinery specialized to the renewal-information kernel.

---

## Executive ruling

| deliverable | ruling |
|---|---|
| **T1 — KKT + upper gap** | **THEOREM.** Full proof below. The exact strong-concavity modulus is the minimum eigenvalue of the pulled-back log-det Hessian on the feasible difference space. A conservative closed-form lower bound is given. `1/lambda_max(V)^2` alone is not a modulus in gain/load coordinates unless the information-map derivative is injective. |
| **T2 — rank-one impossibility** | **THEOREM.** A finite feasible load-transfer with nonzero marginal-information contrast yields an explicit positive lower bound. A two-population witness gives an exact `Omega(1)` gap as brightness contrast grows. |
| **T3 — rank decay** | **THEOREM under exact cluster symmetry; THEOREM for the local quadratic model; CONJECTURE for global nonlinear continuation.** Current leverage values alone do not define exact clusters because leverage changes with the design. |
| **T4 — visibility proxy** | **ASYMPTOTIC PROPOSITION.** Under small bucket variation and a stable leverage law, the true rank-one residual contains `C_u` as a first-order factor. The corrected alignment object also contains the leverage–brightness slope and an orthogonal leverage residual. Identical `C_u` can still produce arbitrarily different control gaps. |
| **T5 — fences** | Beyond the monotone-concave branch, KKT loses sufficiency. With selectable geometry, the R18 `G_stack` quantity is a Frank–Wolfe/support-function gap, not literally `epsilon_infinity`, though both vanish at an optimum and are comparable under geometric radius bounds. Image-domain PSNR/TV claims are not corollaries. |

---

# 1. Natural coordinates and branch assumptions

The clean control coordinate is the **dimensionless load**, not the source gain itself. Define

\[
 b_i=\tau a_i^\top x>0,
 \qquad
 \rho_i=b_i\Phi_i,
 \qquad
 Q_i=q_iq_i^\top .
\]

A zero-brightness pattern has fixed `rho_i=0` and is removed from the controllable coordinate list. If the physical source uses a fixed base-gain profile `g`, a global multiplier gives

\[
 \Phi_i=\kappa g_i,
 \qquad
 \rho_i=\kappa s_i,
 \qquad
 s_i=b_i g_i.
\]

The local task information is

\[
 V(\rho)=V_0+\sum_{i=1}^{M}J(\rho_i)Q_i,
 \qquad
 f(\rho)=\log\det V(\rho),
\]

with `V_0 \succ 0`. For hidden iid hold-time coefficient of variation `c`, use

\[
 J(\rho)=\frac{\rho}{(1+\rho)(1+c^2\rho^2)}.
\]

Its first derivative is

\[
 J'(\rho)=
 \frac{1-c^2\rho^2(1+2\rho)}
 {(1+\rho)^2(1+c^2\rho^2)^2}.
\]

For `c>0`, let `rho_c` be the positive root of

\[
 c^2\rho_c^2(1+2\rho_c)=1.
\]

On the whole increasing branch `0 <= rho <= rho_c`, `J' >= 0` and `J''<0`. One direct check is to put `z=c^2rho^2`. The numerator of `J''/2` becomes

\[
 F(z)=z^2\left(3\rho+3+\rho^{-1}\right)
 -z\left(\rho+6+3\rho^{-1}\right)-1.
\]

The condition `J' >= 0` gives `z <= (1+2rho)^{-1}`, and throughout this interval

\[
 F'(z)
 \le
 -\frac{(\rho+1)(2\rho^2+5\rho+1)}{\rho(2\rho+1)}<0,
\]

so `F(z)<F(0)=-1`. For `c=0`, `J(rho)=rho/(1+rho)` is increasing and strictly concave on `rho>=0`.

The theorem below assumes a compact declared branch

\[
 0\le\rho_i\le\bar\rho_i<\rho_c
\]

(or `bar rho_i < infinity` when `c=0`) and therefore

\[
 \gamma_i:=\min_{0\le t\le\bar\rho_i}\{-J''(t)\}>0.
\]

The campaign constraints become a polyhedron in load coordinates. The incident budget is `1^T rho <= B`. Since

\[
 \Phi_i=\frac{\rho_i}{b_i},
 \qquad
 d_j=\sum_i\frac{a_{ij}}{b_i}\rho_i,
\]

the `+-delta` dose band is a finite set of homogeneous linear inequalities `U rho <= 0`; peak caps become coordinate upper bounds. Write the resulting compact convex set as `C`.

A rank-`k` controller is a fixed load-profile subspace

\[
 \mathcal S_k=\operatorname{range}(R_k)\subseteq\mathbb R^M,
\]

with all nonnegativity and resource restrictions still imposed through `C`. The global source multiplier is `S_1=span{s}`.

---

# 2. T1 — KKT sufficiency and the computable upper gap

## 2.1 Gradient and Hessian

Define the pattern leverage

\[
 \ell_i(\rho)=q_i^\top V(\rho)^{-1}q_i.
\]

Then

\[
 \boxed{\nabla_i f(\rho)=J'(\rho_i)\ell_i(\rho).}
\]

For any direction `h in R^M`, put

\[
 \dot V_\rho[h]=\sum_iJ'(\rho_i)h_iQ_i.
\]

Differentiating `log det` gives

\[
 \boxed{
 -D^2f(\rho)[h,h]
 =
 \left\|V^{-1/2}\dot V_\rho[h]V^{-1/2}\right\|_F^2
 +\sum_i[-J''(\rho_i)]\ell_i(\rho)h_i^2.
 }
 \tag{T1.1}
\]

Equivalently,

\[
 [\nabla^2 f]_{ij}
 =\delta_{ij}J''(\rho_i)\ell_i
 -J'(\rho_i)J'(\rho_j)
 (q_i^\top V^{-1}q_j)^2.
\]

Both terms in (T1.1) are nonnegative after sign reversal, so `f` is concave on the declared branch.

### Proof

Use

\[
 D\log\det V[A]=\operatorname{tr}(V^{-1}A),
\]

and

\[
 D^2\log\det V[A,A]
 =-\operatorname{tr}(V^{-1}AV^{-1}A).
\]

The second derivative of `V(rho+th)` is `sum_i J''(rho_i)h_i^2Q_i`. Substitution gives (T1.1).

---

## 2.2 Exact restricted KKT theorem

### Theorem T1-A — rank-`k` KKT sufficiency

Assume `ri(C) cap S_k` is nonempty. A point

\[
 \rho_k\in C\cap\mathcal S_k
\]

maximizes `f` over `C cap S_k` if and only if there exists

\[
 n_k\in N_C(\rho_k)
\]

such that

\[
 \boxed{
 R_k^\top\{\nabla f(\rho_k)-n_k\}=0.
 }
 \tag{T1.2}
\]

Equivalently,

\[
 \nabla f(\rho_k)\in N_C(\rho_k)+\mathcal S_k^\perp.
\]

This is a necessary **and sufficient** condition, not merely stationarity.

### Proof

Because `f` is concave and `C cap S_k` is convex, `rho_k` is optimal if and only if

\[
 \nabla f(\rho_k)^\top(\rho-\rho_k)\le0
 \quad \forall\rho\in C\cap\mathcal S_k,
\]

that is,

\[
 \nabla f(\rho_k)\in N_{C\cap\mathcal S_k}(\rho_k).
\]

The relative-interior qualification gives the normal-cone sum rule

\[
 N_{C\cap\mathcal S_k}(\rho_k)
 =N_C(\rho_k)+\mathcal S_k^\perp,
\]

which is exactly (T1.2).

---

## 2.3 Exact and conservative strong-concavity moduli

Let

\[
 \mathcal L=\operatorname{span}(C-C)
\]

be the feasible difference space. The exact checkable modulus is

\[
 \boxed{
 m_C=
 \inf_{\rho\in C}
 \inf_{\substack{h\in\mathcal L\\\|h\|_2=1}}
 \left[
 \left\|V^{-1/2}\dot V_\rho[h]V^{-1/2}\right\|_F^2
 +\sum_i[-J''(\rho_i)]\ell_i h_i^2
 \right].
 }
 \tag{T1.3}
\]

If `m_C>0`, `f` is `m_C`-strongly concave on `C` in the Euclidean load metric.

A simple conservative certificate is available. Let

\[
 \Lambda=\lambda_{\max}(V_0)
 +\sum_iJ(\bar\rho_i)\|q_i\|_2^2.
\]

Then `V(rho) <= Lambda I`, hence

\[
 \ell_i(\rho)\ge\frac{\|q_i\|_2^2}{\Lambda},
\]

and therefore

\[
 \boxed{
 m_C\ge
 m_{\rm diag}:=
 \min_i
 \frac{\gamma_i\|q_i\|_2^2}{\Lambda}.
 }
 \tag{T1.4}
\]

This bound is positive after zero-task directions are removed. If it is too conservative, (T1.3) is directly evaluable from the exact Hessian on the active feasible space.

**Important correction to a tempting shortcut.** `log det` is `1/Lambda^2` strongly concave as a function of the information matrix in Frobenius norm on `V <= Lambda I`. That number alone is **not** a Euclidean strong-concavity modulus in `rho`; the pullback `h -> dot V[h]` can have a null space. Formula (T1.3), or the kernel-curvature term in (T1.4), is required.

---

## 2.4 Full-optimality residual and upper gap

Define

\[
 \boxed{
 \varepsilon_k
 =\operatorname{dist}_2
 \left(\nabla f(\rho_k),N_C(\rho_k)\right).
 }
 \tag{T1.5}
\]

By Moreau decomposition,

\[
 \varepsilon_k
 =\left\|\Pi_{T_C(\rho_k)}\nabla f(\rho_k)\right\|_2.
\]

Thus `epsilon_k` is the magnitude of the feasible first-order ascent that remains after the rank-`k` problem has been solved.

### Theorem T1-B — residual-to-gap bound

If `f` is `m_C`-strongly concave on `C`, and `rho*` is the unrestricted maximizer over `C`, then

\[
 \boxed{
 0\le f(\rho^*)-f(\rho_k)
 \le\frac{\varepsilon_k^2}{2m_C}.
 }
 \tag{T1.6}
\]

### Proof

Let `n` be the Euclidean projection of `grad f(rho_k)` onto `N_C(rho_k)` and put `r=grad f-n`. Strong concavity gives, with `Delta=rho*-rho_k`,

\[
 f(\rho^*)-f(\rho_k)
 \le \nabla f(\rho_k)^\top\Delta
 -\frac{m_C}{2}\|\Delta\|_2^2.
\]

Because `n in N_C(rho_k)`, `n^T Delta <=0`. Hence

\[
 f(\rho^*)-f(\rho_k)
 \le \|r\|_2\|\Delta\|_2
 -\frac{m_C}{2}\|\Delta\|_2^2
 \le\frac{\|r\|_2^2}{2m_C}.
\]

Minimizing over `n` proves (T1.6).

---

## 2.5 Fully computable rank-one residual

Let the global-multiplier load profile be

\[
 \rho_1=\kappa s.
\]

Suppose the active face of the polyhedron is represented by

\[
 A_I\rho=b_I.
\]

Then

\[
 N_C(\rho_1)=\{A_I^\top\mu:\mu\ge0\},
\]

including active nonnegativity and peak rows. The exact residual is

\[
 \boxed{
 \varepsilon_1^2
 =\min_{\mu\ge0}
 \left\|h-A_I^\top\mu\right\|_2^2,
 \qquad
 h_i=J'(\rho_{1,i})\ell_i(V_1).
 }
 \tag{T1.7}
\]

This is a small nonnegative least-squares projection using only campaign-computable quantities. Once the active face is fixed, it is piecewise closed form. If `A_I` has full row rank and

\[
 \mu_0=(A_IA_I^\top)^{-1}A_Ih\ge0,
\]

then

\[
 \boxed{
 \varepsilon_1
 =\left\|
 \left[I-A_I^\top(A_IA_I^\top)^{-1}A_I\right]h
 \right\|_2.
 }
 \tag{T1.8}
\]

If some entries of `mu_0` are negative, the exact formula is obtained by the finite active-set NNLS reduction. There is no single scalar closed form valid across every possible combination of active dose and peak faces; the piecewise projection is the mathematically exact answer.

The rank-one restricted KKT condition itself is

\[
 s^\top(h-A_I^\top\widehat\mu)=0.
\]

Multiplying by `kappa` gives the requested shadow-price identity

\[
 \boxed{
 \sum_i\rho_{1,i}J'(\rho_{1,i})\ell_i
 =\widehat\mu^\top b_I.
 }
 \tag{T1.9}
\]

For the incident-budget row, the right side contributes `beta B`. Homogeneous dose-band rows contribute zero under radial scaling. Active peak caps contribute their cap times their shadow price.

### Budget-only interior corollary

If all coordinates are positive, no peak is active, and only `1^T rho=B` is active, then

\[
 N_C(\rho_1)=\{\beta\mathbf1:\beta\ge0\},
\]

and

\[
 \boxed{
 \varepsilon_1^2
 =\sum_i(h_i-\bar h)^2,
 \qquad
 \bar h=M^{-1}\sum_i h_i.
 }
 \tag{T1.10}
\]

Thus rank-one global power is unrestricted-load optimal if and only if the marginal information values

\[
 h_i=\ell_iJ'(\rho_i)
\]

are equal across all active patterns. This is the corrected alignment criterion. Bucket contrast `C_u` is not present in the exact condition.

The multiplier obtained from the one-dimensional rank-one KKT equation is

\[
 \widehat\beta_s
 =\frac{s^\top h}{s^\top\mathbf1},
\]

and the rank-orthogonal residual `h-beta_s 1` is an immediately computable upper bound on the exact projection residual in (T1.10).

---

# 3. T2 — the impossibility side

## 3.1 A finite-transfer lower bound

The upper bound in T1 is not enough to prove that extra control is necessary. The lower bound requires both a marginal-value mismatch and finite room to move inside the physical constraint set.

### Theorem T2-A — feasible-transfer lower bound

Let `rho_1 in C` be the rank-one optimum. Suppose a direction `d` and a number `a>0` satisfy

\[
 \rho_1+td\in C,
 \qquad 0\le t\le a.
\]

Put

\[
 \Delta_d=\nabla f(\rho_1)^\top d>0
\]

and assume along this segment

\[
 -d^\top\nabla^2 f(\rho_1+td)d\le L_d.
\]

Then

\[
 \boxed{
 f(\rho^*)-f(\rho_1)
 \ge
 \Psi(\Delta_d,L_d,a),
 }
\]

where

\[
 \boxed{
 \Psi(\Delta,L,a)=
 \begin{cases}
 \Delta^2/(2L),&a\ge\Delta/L,\\[3pt]
 a\Delta-(L/2)a^2,&a<\Delta/L.
 \end{cases}
 }
 \tag{T2.1}
\]

### Proof

Taylor's theorem with the curvature upper bound gives

\[
 f(\rho_1+td)
 \ge f(\rho_1)+t\Delta_d-\frac{L_d}{2}t^2.
\]

The unrestricted optimum is at least as good as every point on this feasible segment. Maximizing the quadratic over `0<=t<=a` yields (T2.1).

### Pairwise load-transfer form

For a budget-neutral transfer from pattern `i` to pattern `j`, let

\[
 d=e_j-e_i,
 \qquad
 \Delta_{ij}=h_j-h_i.
\]

Whenever `rho_1+t(e_j-e_i)` remains dose/peak feasible for `0<=t<=a_ij`, Theorem T2-A applies. The exact one-dimensional curvature is

\[
 \begin{aligned}
 L_{ij}
 =\sup_t\Bigl[&
 \|V^{-1/2}(J'_jQ_j-J'_iQ_i)V^{-1/2}\|_F^2\\
 &+[-J''_i]\ell_i+[-J''_j]\ell_j
 \Bigr].
 \end{aligned}
 \tag{T2.2}
\]

With `v_0=lambda_min(V_0)`, a coarse analytic bound is

\[
 L_{ij}\le
 \frac{(\bar J'_i\|q_i\|^2+\bar J'_j\|q_j\|^2)^2}{v_0^2}
 +\frac{\bar K_i\|q_i\|^2+\bar K_j\|q_j\|^2}{v_0},
\]

where `bar J'_i=max |J'|` and `bar K_i=max(-J'')` on the transfer interval.

Consequently, if a scene family admits constants

\[
 \Delta_{ij}\ge\Delta_0>0,
 \quad
 a_{ij}\ge a_0>0,
 \quad
 L_{ij}\le L_0<\infty,
\]

then rank-one control leaves a uniform `Omega(1)` log-det gap

\[
 f(\rho^*)-f(\rho_1)
 \ge\Psi(\Delta_0,L_0,a_0)>0.
\]

This is the precise impossibility statement: **misalignment alone is infinitesimal; misalignment plus feasible transfer reach produces a finite gap.**

---

## 3.2 Exact two-population witness

### Theorem T2-B — bright/uninformative versus dim/informative witness

Consider two patterns on the increasing branch, with budget

\[
 \rho_H+\rho_L\le B.
\]

Let their scene brightnesses obey `b_H/b_L=R`, so one global source gain enforces

\[
 \rho_H:\rho_L=R:1.
\]

Take a scalar task with

\[
 V_0=v>0,
 \qquad
 Q_H=0,
 \qquad
 Q_L=1.
\]

Then the unrestricted optimum allocates all load to the informative pattern,

\[
 \rho^*=(0,B),
\]

whereas the rank-one global multiplier gives

\[
 \rho_1=
 \left(\frac{RB}{R+1},\frac{B}{R+1}\right).
\]

The gap is exactly

\[
 \boxed{
 f(\rho^*)-f(\rho_1)
 =
 \log\frac{v+J(B)}{v+J(B/(R+1))}.
 }
 \tag{T2.3}
\]

Hence

\[
 \lim_{R\to\infty}
 [f(\rho^*)-f(\rho_1)]
 =\log\left(1+\frac{J(B)}{v}\right)>0.
\]

The rank-one gap is therefore `Omega(1)` even though the controller knows the scene brightness perfectly. The failure is that brightness and task value point in opposite directions.

### Nondegenerate version

Let

\[
 V_0=vI_2,
 \qquad
 Q_H=\eta e_1e_1^\top,
 \qquad
 Q_L=e_2e_2^\top,
\]

with `eta>0` small. The all-to-`L` feasible design gives determinant `v(v+J(B))`, while the global design has determinant at most

\[
 (v+\eta J(B))
 [v+J(B/(R+1))].
\]

Therefore

\[
 f(\rho^*)-f(\rho_1)
 \ge c_R-\log\left(1+\frac{\eta J(B)}{v}\right),
\]

where

\[
 c_R=\log\frac{v+J(B)}{v+J(B/(R+1))}.
\]

If

\[
 \eta\le\frac{v}{J(B)}(e^{c_R/2}-1),
\]

then

\[
 f(\rho^*)-f(\rho_1)\ge c_R/2.
\]

Thus the witness survives when the bright population is weakly informative rather than exactly null.

The contour family is an empirical analogue of this witness, not a proof that it obeys the two-population model.

---

## 3.3 The corrected alignment functional

For the full physical cone, the exact infinitesimal misalignment is

\[
 \boxed{
 \mathfrak M_1^{\rm exact}
 :=\varepsilon_1
 =\|\Pi_{T_C(\rho_1)}h\|_2,
 \qquad h_i=\ell_iJ'(\rho_i).
 }
 \tag{T2.4}
\]

A normalized score is

\[
 \operatorname{Align}_1
 =1-\frac{\varepsilon_1^2}{\|h\|_2^2},
\]

when `h != 0`. This depends on brightness through the achieved loads, on the task directions through `ell_i`, on the current information matrix, and on all active physical constraints. It replaces `C_u` as the rank-one adequacy diagnostic.

For finite impossibility, pair it with feasible reach through (T2.1). No scene-only scalar can supply that information.

---

# 4. T3 — rank decay

## 4.1 Exact cluster theorem

A theorem based only on `r` currently observed leverage values is false: two patterns can have equal leverage at one matrix `V` and diverge immediately after reweighting. Exact rank reduction requires symmetry of the whole design problem.

### Theorem T3-A — orbit/cluster sufficiency

Partition the patterns into `r` clusters `C_1,...,C_r`. Suppose the product permutation group that permutes indices within each cluster leaves invariant:

1. the feasible set `C`;
2. every information signature `Q_i` and load kernel within the cluster, so `f(P rho)=f(rho)` for each cluster permutation `P`.

Let

\[
 \mathcal S_r
 =\{\rho:\rho_i=\theta_s\text{ for }i\in C_s\}.
\]

Then there exists an unrestricted optimum in `C cap S_r`. If `f` is strictly concave on `C`, the unrestricted optimum is unique and belongs to `S_r`. Consequently, rank `r` cluster control is exact and

\[
 \varepsilon_r=0.
\]

### Proof

Let `rho*` be any optimum and average it over the finite permutation group:

\[
 \bar\rho=|G|^{-1}\sum_{P\in G}P\rho^*.
\]

Convexity and invariance of `C` give `bar rho in C cap S_r`. Concavity and invariance of `f` give

\[
 f(\bar\rho)
 \ge |G|^{-1}\sum_Pf(P\rho^*)
 =f(\rho^*).
\]

Thus `bar rho` is optimal. Strict concavity gives uniqueness.

This is a classical invariant-design argument. The dead-time specialization is the definition of the complete cluster signature: brightness/load law, task matrix, and all resource columns must share the symmetry.

---

## 4.2 Local numerical control rank

For a collection of scenes or operating conditions indexed by `s`, suppose the same active face and tangent space are valid in a neighborhood of a reference control. Let the local quadratic models be

\[
 \widehat f_s(\delta)
 =f_s(0)+r_s^\top\delta-rac12\delta^\top H\delta,
 \qquad H\succ0,
\]

with a common curvature `H`. Define the whitened Newton corrections

\[
 z_s=H^{-1/2}r_s
\]

and form

\[
 Z=[z_1,\ldots,z_S],
 \qquad
 \mathcal M=ZZ^\top
 =\sum_sz_sz_s^\top.
\]

### Theorem T3-B — local singular-value rank decay

Among all shared `k`-dimensional control subspaces, the minimum total quadratic optimality gap is

\[
 \boxed{
 \min_{\dim\mathcal U=k}
 \sum_s
 [\widehat f_s(\delta_s^*)-
 \max_{\delta\in\mathcal U}\widehat f_s(\delta)]
 =\frac12\sum_{j>k}\sigma_j(Z)^2.
 }
 \tag{T3.1}
\]

The optimal whitened control space is the span of the first `k` left singular vectors of `Z`. If

\[
 \operatorname{rank}(\mathcal M)=r,
\]

then rank `r` control is exact for every local quadratic model.

### Proof

In whitened coordinates `w=H^{1/2}delta`,

\[
 \widehat f_s
 =\text{const}
 -\frac12\|w-z_s\|_2^2
 +\frac12\|z_s\|_2^2.
\]

Restriction to a whitened subspace `W` replaces `z_s` by its orthogonal projection `P_Wz_s`, so the gap is `1/2||(I-P_W)z_s||^2`. Summing and applying Eckart–Young–Mirsky proves (T3.1).

### Cluster corollary

If the whitened correction vectors are combinations of `r` cluster-indicator profiles, then `rank(M)<=r`, so the local rank is at most the number of distinct control clusters.

### CONJECTURE T3-C — nonlinear continuation

If the active face is stable, the Hessian is uniformly Lipschitz, and both unrestricted and rank-`k` optima remain inside a fixed trust ball, then the nonlinear average gap equals the quadratic tail in (T3.1) up to a cubic Taylor remainder.

**Named obstacle:** the leverage `ell_i(V)`, normal cone, and active dose/peak faces all move with the control. Active-set bifurcations can change the relevant tangent space discontinuously, and scene-dependent Hessians need not share one whitening operator. Without those stability assumptions, “numerical rank of the current leverage profile” is not a global control-rank theorem.

---

# 5. T4 — how bucket contrast enters, and where it stops

Consider the budget-only interior case of (T1.10). Let the fixed-pattern bucket brightnesses be

\[
 u_i=a_i^\top x,
 \qquad
 u_i=\bar u(1+z_i),
 \qquad
 \frac1M\sum_i z_i=0,
 \qquad
 C_u^2=\frac1M\sum_i z_i^2.
\]

A global source multiplier produces

\[
 \rho_i=\bar\rho(1+z_i),
 \qquad
 \bar\rho=\kappa\tau\bar u.
\]

Write the leverage variation as

\[
 \ell_i=\bar\ell(1+\eta_i),
 \qquad
 \frac1M\sum_i\eta_i=0.
\]

Regress the leverage fluctuation on brightness:

\[
 \eta_i=\alpha z_i+\zeta_i,
 \qquad
 \sum_i z_i\zeta_i=0.
\]

Assume `max_i(|z_i|+|eta_i|)<=delta` and `J in C^3` on the relevant interval.

### Proposition T4-A — first-order visibility/alignment expansion

Let

\[
 \alpha_J^*(\bar\rho)
 =-\frac{\bar\rho J''(\bar\rho)}{J'(\bar\rho)}.
\]

Then

\[
 \boxed{
 \Pi_{\mathbf1^\perp}h
 =\bar\ell J'(\bar\rho)
 \left[(\alpha-\alpha_J^*)z+\zeta\right]
 +R,
 }
 \tag{T4.1}
\]

with

\[
 \|R\|_2\le K\sqrt M\,\delta^2,
\]

where `K` is computable from bounds on `bar ell`, `J'`, `J''`, and `J'''` over the declared branch. Hence

\[
 \boxed{
 \frac{\varepsilon_1^2}
 {M\bar\ell^2J'(\bar\rho)^2}
 =
 (\alpha-\alpha_J^*)^2C_u^2
 +\sigma_\zeta^2
 +O(\delta^3),
 }
 \tag{T4.2}
\]

where

\[
 \sigma_\zeta^2=M^{-1}\sum_i\zeta_i^2.
\]

### Proof

Taylor-expand

\[
 J'(\rho_i)
 =J'(\bar\rho)
 +\bar\rho J''(\bar\rho)z_i
 +O(z_i^2).
\]

Multiplying by `ell_i=bar ell(1+eta_i)`, centering, and using `eta=alpha z+zeta` gives (T4.1). Orthogonality of `z` and `zeta` gives (T4.2).

### Homogeneous-leverage corollary

If `eta_i=0`, then

\[
 \boxed{
 \varepsilon_1
 =\sqrt M\,\bar\ell
 |\bar\rho J''(\bar\rho)|C_u
 +O(\sqrt M\delta^2).
 }
 \tag{T4.3}
\]

Thus `C_u` is a monotone first-order proxy only when task leverage is homogeneous or changes according to a fixed law.

### Exact first-order alignment slope

The first-order residual vanishes when

\[
 \boxed{
 \alpha=\alpha_J^*(\bar\rho)
 =-\frac{\bar\rho J''(\bar\rho)}{J'(\bar\rho)}.
 }
 \tag{T4.4}
\]

For deterministic holds, `J=rho/(1+rho)`, so

\[
 \alpha_J^*(\rho)=\frac{2\rho}{1+\rho}.
\]

For random holds, `J'(rho)` tends to zero at the jitter-limited ridge while `J''` remains negative, so `alpha_J^*` diverges as the ridge is approached. This makes one-knob alignment increasingly fragile near the cap.

### Interpretation

The corrected rank-one alignment functional is therefore

\[
 \boxed{
 \mathfrak M_{1,\mathrm{local}}^2
 =(\alpha-\alpha_J^*)^2C_u^2+\sigma_\zeta^2.
 }
 \tag{T4.5}
\]

This is the mathematical bridge to the operating-map axis:

- `C_u` is the visibility amplitude;
- `alpha-alpha_J^*` measures whether leverage changes with brightness in the way the dead-time kernel requires;
- `sigma_zeta` measures task geometry invisible to brightness.

If `alpha` and `sigma_zeta` are stable across a restricted scene family, ordering by `C_u` is justified to first order. If they are not, identical `C_u` can correspond to zero or large rank-one gap. The R22 identical-contrast counterexample is exactly the `sigma_zeta` failure mode.

This proposition does **not** turn the paper-1 `Gamma=1` visibility proxy into an optimal-control boundary. It only explains why `C_u` can correlate with one-knob performance in a homogeneous task family and why that correlation breaks for contour-like geometry.

---

# 6. T5 — scope fences

## 6.1 Beyond the ridge/nonconcave branch

When the feasible loads cross a region where `J' < 0` or `J'' > 0`, `f` need not be concave in the control. The KKT condition remains necessary under standard constraint qualifications but is no longer sufficient; (T1.6) and the global cluster-symmetrization conclusion can fail. The local Hessian formulas remain valid, so a theorem can still be stated on a connected trust region where the pulled-back Hessian is uniformly negative definite. Outside such a certified basin, global branch selection is a genuinely nonconvex problem.

## 6.2 Geometry as part of the control

When patterns are selected from a dictionary, the natural variable is a design measure `xi` and

\[
 V(\xi)=V_0+M\int H(a,\rho)\,d\xi(a,\rho).
\]

The R18 quantity

\[
 G_{\rm stack}
 =\sup_{\xi\in C_{\rm stack}}
 \langle\nabla f(\xi_d),\xi-\xi_d\rangle
\]

is the constrained Frank–Wolfe/Kiefer–Wolfowitz support-function gap. It is **not literally**

\[
 \varepsilon_\infty
 =\operatorname{dist}(\nabla f,N_C).
\]

Both vanish exactly at a convex optimum. If the feasible-set diameter is `D`, then

\[
 G_{\rm stack}\le D\varepsilon_\infty.
\]

If the tangent set contains a ball of radius `r` in the projected-gradient direction, then

\[
 G_{\rm stack}\ge r\varepsilon_\infty.
\]

Without such geometric radii there is no equality. Thus `G_stack` is the measure-space analogue of the residual doctrine, but not the same numerical object.

## 6.3 Image-domain claims

Everything above concerns local Fisher information and `log det` on a declared task subspace. It does not imply ordering in PSNR, SSIM, TV reconstruction error, hallucination rate, or any nonlinear estimator risk. A Fisher-improving design can worsen a regularized reconstruction through bias, model mismatch, or task-subspace misspecification. Image-domain claims require separate analysis or experiment and are explicit non-corollaries.

---

# 7. Prior-art and novelty ruling

The components must be attributed aggressively.

## Classical or substantially known

1. **KKT/general equivalence.** Kiefer and Wolfowitz established the equivalence architecture for optimal designs; Kiefer generalized the theory, and White extended it to nonlinear models:
   - Kiefer, *Optimum Experimental Designs*, JRSS B 21 (1959), DOI: `10.1111/j.2517-6161.1959.tb00338.x`.
   - Kiefer & Wolfowitz, *The Equivalence of Two Extremum Problems* (1960), DOI: `10.4153/CJM-1960-030-4`.
   - Kiefer, *General Equivalence Theory for Optimum Designs* (1974), DOI: `10.1214/AOS/1176342810`.
   - White, *An extension of the General Equivalence Theorem to nonlinear models* (1973), DOI: `10.1093/biomet/60.2.345`.
   - Pukelsheim, *Optimal Design of Experiments*, SIAM, DOI: `10.1137/1.9780898719109`.

2. **Restricted/resource-constrained designs.** The idea that a design optimum and its equivalence theorem must be computed inside a restricted convex design class is established:
   - Cook & Thibodeau, *Marginally Restricted D-Optimal Designs* (1980), DOI: `10.1080/01621459.1980.10477478`.
   - Cook & Wong, *On the Equivalence of Constrained and Compound Optimal Designs* (1994), DOI: `10.1080/01621459.1994.10476794`.
   - Harman, Bachrata & Filova, *Construction of efficient experimental designs under multiple resource constraints* (2016), DOI: `10.1002/asmb.2117`.

3. **Normal-cone residual and strong-concavity bounds.** These are standard convex-analysis consequences; see Boyd & Vandenberghe, *Convex Optimization*, DOI: `10.1017/CBO9780511804441`, and Nesterov, *Introductory Lectures on Convex Optimization*, DOI: `10.1007/978-1-4419-8853-9`.

4. **Low-dimensional parameterizations and low-rank computation in OED.** Existing work uses reduced design parameterizations or low-rank forward structure computationally, but does not appear to state the dead-time rank-one adequacy/impossibility package above:
   - Ruthotto, Chung & Chung, *Optimal Experimental Design for Inverse Problems with State Constraints* (2018), DOI: `10.1137/17M1143733`.
   - Alexanderian & Saibaba, *Efficient D-Optimal Design of Experiments for Infinite-Dimensional Bayesian Linear Inverse Problems* (2018), DOI: `10.1137/17M115712X`.

5. **Cluster symmetry and singular-value rank decay.** Symmetrization belongs to invariant-design theory; the quadratic rank theorem is Eckart–Young–Mirsky after Hessian whitening. Neither abstract ingredient is new.

6. **Bright/dim multiplex disadvantage.** The physical principle that multiplexing can waste photon budget on bright or task-irrelevant components is established in optical multiplex analysis:
   - Oliver & Pike, *Multiplex advantage in the detection of optical images in the photon noise limit* (1974), DOI: `10.1364/AO.13.000158`.
   - Schechner, Nayar & Belhumeur, *Illumination Multiplexing within Fundamental Limits* (2007), DOI: `10.1109/CVPR.2007.383162`.
   - Scotte, Galland & Rigneault, *Photon-noise: is a single-pixel camera better than point scanning?* (2023), DOI: `10.1088/2515-7647/acc70b`.

7. **Task-space information.** The use of Fisher information and task subspaces for imaging assessment is established; see Barrett et al., *Objective assessment of image quality II* (1995), DOI: `10.1364/JOSAA.12.000834`.

## Narrow residual claim that may be worth a paper

Do **not** claim a new general rank-constrained OED theory. A defensible residual is:

> For a physically constrained, nonnegative photon-counting experiment with hidden dead-time state, derive and test a load-coordinate rank-control adequacy certificate, its finite-transfer impossibility bound, and the kernel-specific leverage–brightness elasticity `alpha_J^*=-rho J''/J'` that determines when one global power control approximates unrestricted pattern-wise gains.

The strongest original-looking object is (T4.5), together with the explicit two-sided rank-one bounds. T1 and T3 are the necessary classical scaffolding.

---

# 8. Final operative summary

1. **Use load coordinates.** They make the resource cone linear and expose the exact marginal value `h_i=ell_i J'(rho_i)`.
2. **Define adequacy by the normal-cone residual**, not by `C_u`:
   \[
   \varepsilon_k=\operatorname{dist}(\nabla f,N_C).
   \]
3. **Upper-bound the unrestricted loss** by `epsilon_k^2/(2m_C)` with the explicit pulled-back Hessian modulus (T1.3) or the conservative bound (T1.4).
4. **Prove necessity with finite transfer**, not an angle alone. The lower bound is (T2.1); the two-population witness gives an exact constant gap.
5. **Rank equals cluster count only under full problem symmetry.** Current leverage clustering alone is insufficient.
6. **The corrected bridge to bucket contrast is two-factor:**
   \[
   (\alpha-\alpha_J^*)C_u
   \quad\text{plus}\quad
   \sigma_\zeta.
   \]
7. **Do not identify `epsilon_infinity` with `G_stack`.** One is a normal-cone distance; the other is a support-function/Frank–Wolfe gap.
8. **Do not claim image-quality corollaries.** The theorem is local-information theory.

This is mathematically sufficient to begin a focused next paper. It also gives a clear falsifier: if the measured leverage–brightness slope and orthogonal residual do not predict the observed need for more than one gain profile, the proposed bridge fails rather than expanding into another layer of terminology.