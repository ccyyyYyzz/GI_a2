# R14 — Final ruling: unification architecture

**Scope:** this ruling is independent of the R13/M1 implementation hold. It governs the theory architecture of the follow-up method line only. It does not relax any R13 freeze blocker, reopen the OE manuscript, or alter the frozen Study-1/2 verdicts.

## Executive ruling

| Candidate | R14 ruling |
|---|---|
| **B — random-hold universality** | **Substantially confirmed, with a corrected law and corrected constant.** The exact missing-information identity extends to arbitrary hidden iid hold times. Random hold variance creates an extensive count-only information loss. The long-window count information is `J_inf(rho,c)=rho/[(1+rho)(1+c^2 rho^2)]`, so the unique jitter-limited optimum solves `c^2 rho^2(1+2rho)=1`. In the small-jitter limit, `rho*=2^(-1/3)c^(-2/3)-1/6+...`; the asymptotic coefficient is `2^(-1/3)=0.79370`, not `0.72`. Combining this with the deterministic finite-window boundary term gives the matched crossover `rho* ~ [2c^2+1/(6nu)]^(-1/3)`, not a literal `min`. The exponent pair `(1/3,-2/3)` is universal for a regular finite-variance iid hidden-hold class, not for every possible dead-time model. |
| **A — one master concave program** | **Approved as the framework theorem.** Every true interior optimum of the declared local D-optimal design problem satisfies one generalized KKT/reduced-sensitivity system. Load, support scale, and weight parameters are coordinate slices of that system. However, the empirically observed `k≈32` and `p≈1` endpoint optima are not automatically KKT theorems: they were discrete and judged by image endpoints, not necessarily by the master log-det objective. |
| **C — brightness/information alignment** | **Approved as a conditional proposition, not a universal approximation theorem.** Fixed global flux is exactly water-filling-optimal when `ell_i J'(kappa u_i)` is constant across active directions. A quantitative near-optimality bound follows from the projected KKT residual when `J` is strongly concave on the admitted load interval. Without this alignment, fixed flux can be arbitrarily bad. |
| **Paper-2 theory spine** | **B is the headline theorem; A is the unifying method framework; C is the explanatory proposition for RIDGE-FIXED/self-water-filling.** The deterministic detector is the zero-jitter critical manifold of B. A consumes the generalized kernel supplied by B, and C explains when the full variable-load optimization collapses approximately to one global source knob. |

---

# 1. Candidate B — random hold times and the information crossover

## 1.1 Model and observation

Let the detector be active at `t=0`. While active, the next registered event occurs after

\[
X_i\sim\operatorname{Exp}(\lambda).
\]

After each registered event the detector is blocked for an iid hidden hold time

\[
B_i\sim H,
\qquad
\mathbb E B_i=\tau,
\qquad
\operatorname{Var}B_i=\sigma_B^2=c^2\tau^2,
\]

independent of the Poisson arrivals and independent of `lambda`. The observed datum is only the integrated registered count

\[
N_T,
\qquad T=\nu\tau,
\qquad \rho=\lambda\tau.
\]

Let

\[
L_T=\int_0^T A_t\,dt
\]

be the total time for which the detector is active, where `A_t` is the active-state indicator. The realized `B_i`, the state path, and `L_T` are hidden from the count-only observer.

The mean-cycle statement in the R14 question is correct only asymptotically:

\[
\mathbb E N_T
=
\nu\frac{\rho}{1+\rho}+O(1).
\]

It is not an exact finite-window equality; the `O(1)` boundary term depends on the start convention and on `H`.

## 1.2 Exact missing-information identity

### Theorem R14.B1 — hidden-hold count-information identity

Assume the law of `B` does not depend on `lambda`. For `theta=log lambda`, the complete detector-path score is

\[
U_T
=
\frac{\partial}{\partial\theta}\log p_\theta(\text{path})
=
N_T-\lambda L_T.
\]

The full-path Fisher information is

\[
I_{\rm path}(\theta;T)
=
\mathbb E[U_T^2]
=
\mathbb E[N_T].
\]

The Fisher information retained by the scalar count is exactly

\[
\boxed{
I_N(\theta;T)
=
\mathbb E[N_T]
-
\lambda^2\,
\mathbb E\!\left[\operatorname{Var}(L_T\mid N_T)\right].
}
\]

### Proof

Conditional on the detector state path, registered events form a counting process with intensity `lambda A_t`; hence

\[
N_t-\lambda\int_0^tA_sds
\]

is a martingale whose predictable quadratic variation is `lambda L_T`. Therefore

\[
\mathbb E U_T^2=\lambda\mathbb E L_T=\mathbb E N_T.
\]

Fisher's identity gives the count score

\[
U_N=\mathbb E(U_T\mid N_T).
\]

The law of total variance then gives

\[
\operatorname{Var}(U_T)
=
\operatorname{Var}\{\mathbb E(U_T\mid N_T)\}
+
\mathbb E\{\operatorname{Var}(U_T\mid N_T)\}.
\]

Because `N_T` is fixed under the conditional variance,

\[
\operatorname{Var}(U_T\mid N_T)
=
\lambda^2\operatorname{Var}(L_T\mid N_T),
\]

which proves the identity.

This theorem is exact for every finite `T` and every iid hidden hold law `H`. The deterministic identity in the OE paper is the special case in which uncertainty in `L_T` is confined to the terminal detector phase.

### Scope failures

The identity requires amendment if:

- the hold law depends on `lambda` or the scene parameter;
- realized hold durations are themselves observed and used;
- timestamps/state transitions rather than one count are retained;
- holds are serially correlated;
- the counter is paralyzable and arrivals during a hold alter the state.

## 1.3 Long-window information theorem

Let one ordinary renewal cycle have duration

\[
C_i=B_i+X_i.
\]

Its mean and variance are

\[
\mu_C
=
\tau+\frac1\lambda
=
\tau\left(1+\frac1\rho\right),
\]

\[
v_C
=
\sigma_B^2+\frac1{\lambda^2}
=
\tau^2\left(c^2+\frac1{\rho^2}\right).
\]

Under the standard nonlattice renewal local-limit/LAN conditions—for example, `H` lognormal, gamma/exponential, or bounded with a finite `3+epsilon` moment—the renewal count obeys

\[
\mathbb E N_T
=
\frac{T}{\mu_C}+O(1)
=
\nu\frac{\rho}{1+\rho}+O(1),
\]

\[
\operatorname{Var}N_T
=
\frac{v_CT}{\mu_C^3}+O(1)
=
\nu\frac{\rho(1+c^2\rho^2)}{(1+\rho)^3}+O(1).
\]

Moreover,

\[
\frac{\partial}{\partial\theta}\mathbb E N_T
=
\nu\frac{\rho}{(1+\rho)^2}+O(1).
\]

The parameter derivative of the asymptotic variance contributes only `O(1)` Fisher information, whereas the mean shift contributes `O(nu)`. Consequently:

### Theorem R14.B2 — long-window count-information rate

\[
\boxed{
\frac{I_N(\log\lambda;T)}{\nu}
\longrightarrow
J_\infty(\rho,c)
=
\frac{\rho}{(1+\rho)(1+c^2\rho^2)}.
}
\]

The leading coefficient depends on `H` only through its mean and variance.

A compact proof uses the renewal local normal approximation:

\[
I_N
=
\frac{\{\partial_\theta\mathbb EN_T\}^2}
     {\operatorname{Var}N_T}
+O(1),
\]

which yields the displayed expression. A second proof combines the exact identity R14.B1 with the joint renewal CLT for `(N_T,L_T)`.

## 1.4 The extensive missing-information term

Subtracting `J_inf` from the complete-path information rate gives

\[
\boxed{
\frac{\lambda^2}{\nu}
\mathbb E\!\left[\operatorname{Var}(L_T\mid N_T)\right]
\longrightarrow
\frac{c^2\rho^3}
{(1+\rho)(1+c^2\rho^2)}.
}
\]

When `c^2 rho^2 << 1`,

\[
\frac{\lambda^2}{\nu}
\mathbb E\!\left[\operatorname{Var}(L_T\mid N_T)\right]
=
\frac{c^2\rho^3}{1+\rho}
+o(c^2\rho^2).
\]

The raw heuristic from the question,

\[
\lambda^2\sigma_B^2\mathbb E N_T,
\]

is therefore the correct first-order extensive term. It is not exact: conditioning on the total count partially reveals the accumulated hold-time fluctuation, producing the factor

\[
(1+c^2\rho^2)^{-1}.
\]

At large load and small jitter, the per-slot extensive loss is `c^2 rho^2` to leading order. Unlike the deterministic terminal-phase loss, it remains nonzero as `nu -> infinity`.

## 1.5 Exact long-window jitter optimum

For every fixed `c>0`, `J_inf(rho,c)` has a unique interior maximum. Differentiation gives

\[
\frac1\rho-
\frac1{1+\rho}-
\frac{2c^2\rho}{1+c^2\rho^2}=0,
\]

or equivalently

\[
\boxed{
c^2\rho_c^2(1+2\rho_c)=1.}
\]

This cubic, not a fitted empirical constant, is the long-window jitter-limited load law.

For small `c`,

\[
\boxed{
\rho_c
=
2^{-1/3}c^{-2/3}
-
\frac16
+
\frac{2^{1/3}}{36}c^{2/3}
+O(c^{4/3}).
}
\]

Thus the universal leading coefficient is

\[
\boxed{2^{-1/3}=0.793700526\ldots}
\]

and not `0.72`.

The current Monte-Carlo values are directionally consistent but are not yet a coefficient measurement:

| `c` | long-window cubic prediction | current MC peak |
|---:|---:|---:|
| 0.05 | 5.6860 | about 5.43 |
| 0.10 | 3.5247 | about 3.30 |
| 0.30 | 1.6191 | at or below 2, grid-edge limited |

The measured `0.72` is an effective finite-grid/finite-sample coefficient. The ratio test is less sensitive to common bias and therefore happened to agree more closely with `-2/3`.

### Lognormal holds

For

\[
B=\tau\exp\{\sigma Z-\sigma^2/2\},
\qquad
c^2=e^{\sigma^2}-1,
\]

Theorem R14.B2 and the same cubic apply. Lognormality changes only finite-window/higher-order corrections through skewness and higher cumulants.

### Exponential holds

An exponential hold has `c=1`, so

\[
J_\infty(\rho,1)
=
\frac{\rho}{(1+\rho)(1+\rho^2)}.
\]

The exact positive solution of

\[
\rho^2(1+2\rho)=1
\]

is

\[
\boxed{\rho_{\rm exp}=0.657298106\ldots .}
\]

Because `c=1` is not a small-jitter limit, it should not be approximated by `2^{-1/3}c^{-2/3}`.

Gamma, bounded two-point, and lognormal holds having the same small `c` share the leading coefficient; distribution shape enters at the next order.

## 1.6 Finite-window crossover: replace `min` by a scaling function

On the deterministic line, the OE derivation gives

\[
J_{\nu,0}(\rho)
=
1-rac1\rho-rac{\rho^2}{12\nu}
+\text{lower-order terms}.
\]

For small random hold variance, the long-window expansion gives

\[
J_{\infty,c}(\rho)
=
1-rac1\rho-c^2\rho^2
+\text{lower-order terms}.
\]

In the joint regime

\[
\nu\to\infty,
\qquad c\to0,
\qquad \rho\to\infty,
\]

with `c^2 nu=O(1)`, the matched leading expansion is

\[
\boxed{
J_{\nu,c}(\rho)
=
1-rac1\rho
-ho^2\left(c^2+\frac1{12\nu}\right)
+o\!\left(rac1\rho+c^2\rho^2+rac{\rho^2}{\nu}\right).
}
\]

Balancing marginal gain and loss gives

\[
\boxed{
\rho_*(\nu,c)
\sim
\left(2c^2+\frac1{6\nu}\right)^{-1/3}
=
(6\nu)^{1/3}
(1+12\nu c^2)^{-1/3}.
}
\]

This is the correct crossover architecture. The proposed

\[
\min\{(6\nu)^{1/3},\,c_0c^{-2/3}\}
\]

is only the two-limit envelope and has an artificial kink.

The scaling variable is

\[
\boxed{x=12\nu c^2.}
\]

The limits are:

- **deterministic/critical:** `x -> 0`,
  \[
  \rho_*\sim(6\nu)^{1/3};
  \]
- **jitter dominated near the critical manifold:** `x -> infinity` with `c -> 0`,
  \[
  \rho_*\sim2^{-1/3}c^{-2/3}.
  \]

The crossover occurs at

\[
c\asymp(12\nu)^{-1/2}.
\]

## 1.7 What is rigorous now, and what still needs a dedicated proof

### Provable immediately at manuscript-theorem rigor

1. exact identity R14.B1;
2. renewal mean/variance formulas;
3. long-window rate R14.B2 under a stated renewal local-limit/LAN assumption;
4. the exact cubic optimum;
5. the small-`c` expansion and constants;
6. the KKT results in Candidate A;
7. the alignment proposition in Candidate C.

### Dedicated proof required before calling the full crossover a theorem

The uniform two-parameter expansion

\[
J_{\nu,c}
=
1-\rho^{-1}-c^2\rho^2-\rho^2/(12\nu)+\cdots
\]

requires a triangular-array renewal local-limit/Edgeworth argument uniform in

\[
c\to0,
\quad \rho\asymp\nu^{1/3},
\quad c^2\nu=O(1).
\]

In particular, one must control:

- the start/terminal phase uniformly;
- third and fourth hold-time cumulants;
- differentiation of the count local limit with respect to `log lambda`;
- the remainder uniformly over the peak neighborhood.

Until that proof is complete, freeze the language as:

> The exact missing-information identity and long-window renewal limit imply the jitter exponent `-2/3`; matching this extensive loss to the deterministic terminal-phase loss predicts the crossover scaling function `(6nu)^(1/3)(1+12nu c^2)^(-1/3)`.

Do not yet write “we prove the full crossover law” solely from the current Monte Carlo.

## 1.8 Universality class and its boundary

The exponent pair `(1/3,-2/3)` is universal within the following class:

- iid hidden nonnegative holds;
- hold law independent of `lambda`;
- finite, nonzero mean and finite variance;
- a regular small-jitter family converging to deterministic hold in `L^{3+epsilon}`;
- count-only observation;
- active-start nonparalyzable renewal operation.

Within this class:

- `1/3` comes from balancing the `1/rho` saturation gap with the `rho^2/nu` terminal-boundary loss;
- `-2/3` comes from balancing the same `1/rho` gap with the extensive `c^2 rho^2` hidden-hold loss.

“Deterministic dead time is the critical line” may be used only with the more precise wording:

> zero hold-time variance is the critical manifold on which the optimal load continues to grow with window length; any fixed nonzero hidden hold variance produces a finite long-window optimum.

The universality statement does **not** cover:

- infinite-variance/heavy-tailed holds, where a stable-law exponent can replace `2`;
- correlated or drifting holds;
- holds observed by the electronics;
- paralyzable detectors;
- event-time observations;
- `lambda`-dependent recovery dynamics.

## 1.9 Prior-art verdict

A targeted search found substantial adjacent theory but no located statement of the above count-only Fisher crossover:

- renewal-CLT analyses of dead-time counters establish the mean/variance ingredients used here; Alvarez's pileup analysis is a close deterministic-count surrogate application: https://doi.org/10.1118/1.4898102
- Dubi and Atar prove diffusion-scale limits for counters with general random dead time, but for extendable/type-II counters and without this count-only Fisher optimum: https://doi.org/10.1016/j.net.2019.04.015
- He, Yang, Fang, and Widmann study Poisson-intensity estimation with recurrent dead-time gaps and residual lifetimes, but not the hidden-random-hold scalar-count ridge: https://www.nist.gov/publications/consistent-estimation-poisson-intensity-presence-dead-time
- Jorgensen and Johnson derive LAN and Fisher-information rates for richer periodic event-detection records with fixed dead time and gating, not the random-hold integrated-count crossover: https://arxiv.org/abs/2605.23210
- Müller, Yu–Fessler, Grönberg, and the nuclear-counting literature establish count laws, moments, and task-specific finite-rate optima, but no located source gives the exponent pair or the extensive/terminal crossover.

Safe novelty wording:

> For count-only nonparalyzable renewal detection with hidden random holds, we derive an exact missing-information identity and a long-window Fisher-information rate whose optimum is finite for every nonzero hold-time variance. Matching that extensive loss to the deterministic terminal-phase loss yields a two-parameter crossover between the `nu^(1/3)` and `cv^(-2/3)` regimes.

Do not claim that renewal CLTs, random dead-time limits, or finite-rate detector optima are new.

---

# 2. Candidate A — one constrained concave program

## 2.1 Master program

Let an atom be

\[
z=(\text{support family},\text{translation},p,\rho),
\]

with local chain-rule direction `q_z` and generalized hold-law kernel `J_H(rho,nu)`. Define

\[
H_z
=
\nu J_H(\rho_z,\nu)q_zq_z^T.
\]

For a design probability measure `xi`,

\[
V(\xi)
=
V_0+M\int H_zd\xi(z).
\]

The master problem is

\[
\boxed{
\max_{\xi}
\log\det V(\xi)
\quad\text{subject to}\quad
\int c_m(z)d\xi(z)\le C_m,
\quad \int d\xi=1.
}
\]

This remains concave in `xi` for deterministic or random holds because all detector physics is inside the positive-semidefinite atom `H_z`.

## 2.2 Reduced-sensitivity/KKT unification

At an optimum, with resource multipliers `beta_m`, define

\[
\psi(z)
=
M\operatorname{tr}[V^{-1}H_z]
-
\sum_m\beta_mc_m(z).
\]

The generalized equivalence condition is

\[
\psi(z)\le\eta
\]

for all feasible atoms, with equality on support atoms.

If a support atom is interior in a smooth scalar coordinate `y`, then

\[
\boxed{
\partial_y\psi(z)=0,
\qquad
\partial_y^2\psi(z)\le0.
}
\]

For `H=nu J q q^T`,

\[
\partial_y\operatorname{tr}(V^{-1}H)
=
\nu J'(\rho)\rho_y\,\ell
+
2\nu J(\rho)q_y^TV^{-1}q,
\]

where

\[
\ell=q^TV^{-1}q.
\]

### Load coordinate

If amplitude changes the load but leaves `q` invariant,

\[
\boxed{
M\nu J'(\rho)\ell
=
\text{marginal resource shadow price}.
}
\]

This is the information-water-filling equation.

### Support-scale coordinate

At fixed load,

\[
\boxed{
2M\nu J(\rho)q_s^TV^{-1}q
=
\text{marginal support/dose/peak price}.
}
\]

The apparent “conditioning penalty” is not a separate ad hoc term; it is the `V^{-1}` geometry inside the reduced sensitivity.

### Weight-power coordinate

Likewise,

\[
\boxed{
2M\nu J(\rho)q_p^TV^{-1}q
=
\text{marginal weight-dispersion/peak/dose price}.
}
\]

This explains why increasing proxy contrast from `p=1` to `p=2` can reduce the true objective: the change in direction can worsen global matrix conditioning even when `C_u` rises.

## 2.3 Theorem versus narrative

Candidate A is a theorem about the declared master problem, but not a theorem that every empirical optimum observed so far is one of its KKT points.

The observed `k≈32` and `p≈1` results were:

- selected from discrete grids;
- evaluated by PSNR/Q90 rather than necessarily by `log det`;
- produced under different fixed pattern families and constraints.

For discrete `k` or `p`, the correct KKT statement is a reduced-sensitivity inequality/comparison, not a derivative equal to zero.

Frozen wording:

> The master design problem places load, support scale, and intensity weighting in one constrained concave measure optimization. Its support atoms satisfy a common reduced-sensitivity condition; the empirical occupancy and weighting optima motivate these coordinates but are not themselves claimed as continuous KKT theorems.

Candidate A should be the framework section, not the headline theorem.

---

# 3. Candidate C — brightness/information alignment

## 3.1 Scalar allocation problem at a fixed design iterate

Fix candidate directions `q_i`, design frequencies `pi_i`, and the current matrix `V`. Let

\[
\ell_i=q_i^TV^{-1}q_i
\]

be the local directional value. Allocate scalar loads under

\[
\sum_i\pi_i\rho_i=B,
\qquad
\rho_i\in[\rho_{\min},\rho_{\max}],
\]

to maximize

\[
F(\rho)
=
\sum_i\pi_i\ell_iJ(\rho_i).
\]

On an interval where `J` is increasing and strictly concave, the water-filling optimum satisfies

\[
\boxed{
\ell_iJ'(\rho_i^*)=\beta
}
\]

for every interior coordinate, with the usual inequalities at the bounds.

A global fixed source gives

\[
\widetilde\rho_i=\kappa u_i,
\qquad
u_i=a_i^Tx
\]

(`u_i`, not `nu_i`, is the bucket brightness; notation corrected hereafter).

## 3.2 Exact brightness-alignment condition

The fixed-flux allocation is exactly optimal for this scalar subproblem if and only if there exists `beta` such that

\[
\boxed{
\ell_iJ'(\kappa u_i)=\beta
}
\]

for all interior active directions, with KKT-compatible inequalities for clipped directions.

This is the precise definition of brightness–information alignment.

It says neither “bright is always informative” nor “load should be uniform.” It says that the source-induced load ordering must match the inverse marginal-information demand of the directions.

## 3.3 Quantitative near-optimality bound

Assume

\[
-\ell_iJ''(\rho)\ge m>0
\]

throughout the admitted interval. Choose `beta` as the projection of the fixed-flux gradient onto the budget-normal direction and define residuals

\[
r_i
=
\ell_iJ'(\kappa u_i)-\beta.
\]

Then strong concavity gives

\[
\boxed{
0\le
F(\rho^*)-F(\widetilde\rho)
\le
\frac1{2m}
\sum_i\pi_i r_i^2.
}
\]

Thus a computable projected KKT residual—not `C_u` alone—certifies whether one global flux knob is close to the scalar water-filling optimum.

## 3.4 Power-law interpretation

If over the relevant interval

\[
J'(\rho)\asymp C\rho^{-s}
\]

and

\[
\ell_i\asymp u_i^\alpha,
\]

then fixed flux aligns when

\[
\boxed{\alpha\approx s.}
\]

For the deterministic long-window rising kernel

\[
J(\rho)=\frac\rho{1+\rho},
\qquad
J'(\rho)=\frac1{(1+\rho)^2},
\]

the local log-slope is

\[
s(\rho)
=-\frac{d\log J'}{d\log\rho}
=
\frac{2\rho}{1+\rho}.
\]

It ranges from `0` at very low load toward `2` at high pre-ridge load. Hence the often-invoked `ell proportional to u^2` alignment is only a high-load approximation, not a universal law.

With random holds, `J'` vanishes at the jitter cap and becomes negative beyond it; fixed flux can then overdrive the brightest directions, and explicit gain control becomes necessary.

## 3.5 Failure modes and impossibility of a universal constant factor

Without the residual/alignment condition there is no scene-independent constant-factor guarantee. Construct two directions with

\[
u_1/u_2\to\infty
\]

but

\[
\ell_1/\ell_2\to0.
\]

Fixed flux allocates most load to direction 1 while water filling allocates it to direction 2; the objective ratio can approach zero.

Other failures are:

- the pre-scan misses a low-amplitude but high-value direction;
- brightness is dominated by DC rather than the target task;
- a direction is bright but redundant with the current information matrix;
- load crosses the scalar ridge;
- peak/dose bounds clip the fixed-flux allocation;
- dark counts break simple scale invariance.

## 3.6 Interpretation of RIDGE-FIXED

RIDGE-FIXED is near-optimal only when its realized residuals

\[
\ell_iJ'(\kappa u_i)-\beta
\]

are small. The uniform-servo failure is then transparent: forcing equal `rho_i` discards a naturally favorable alignment between `u_i` and `ell_i` on sparse-support scenes.

Frozen wording:

> Fixed flux performs implicit information water filling when bucket brightness and local directional value satisfy the marginal-information alignment condition. Its adequacy is measured by a projected KKT residual; outside that condition, variable-load OED can improve without any universal bound on the fixed-flux approximation ratio.

Candidate C is a proposition/explanation, not the paper's primary theorem.

---

# 4. Frozen paper-2 theory architecture

## Theory spine

### Section I — random-hold count information

1. exact missing-information identity R14.B1;
2. long-window rate R14.B2;
3. exact cubic optimum;
4. deterministic and jitter exponents;
5. matched crossover prediction;
6. hardware warning: even small hidden hold variability can make deterministic-ridge operation severely information-decreasing.

### Section II — dead-time-aware OED master program

1. generalized kernel `J_H` enters the rank-one atoms;
2. the design measure remains a concave program;
3. generalized equivalence/KKT certificate;
4. load, scale, and weighting as coordinate stationarity conditions.

### Section III — self-water-filling proposition

1. exact alignment condition;
2. residual-based near-optimality bound;
3. uniform-servo refutation as a predicted failure mode;
4. when RIDGE-FIXED is sufficient and when OED-DT is needed.

This order produces one coherent chain:

\[
\text{detector recovery law}
\longrightarrow
J_H(\rho,\nu)
\longrightarrow
\text{constrained information allocation}
\longrightarrow
\text{fixed-flux special case}.
\]

## Claim discipline

Permitted after proving B1/B2:

> Hidden hold-time variability produces an extensive loss of count-only Fisher information and a finite long-window optimum.

> In the small-jitter regime the optimum scales as `cv^(-2/3)`.

Permitted only after the uniform two-parameter proof, or otherwise labeled a prediction:

> The full finite-window crossover follows `(6nu)^(1/3)(1+12nu cv^2)^(-1/3)`.

Not permitted:

- “all detector jitter obeys the same law”;
- “0.72 is the universal constant”;
- “every observed empirical optimum is a KKT point”;
- “fixed flux is universally within a constant of OED”;
- “dead-time randomness has not previously been studied.”

---

# 5. Refined numerical sweep for an exponent measurement

## 5.1 Primary estimator: use the exact missing-information score, not finite differences

For every simulated frame, record

\[
N_T,
\qquad L_T,
\qquad U_T=N_T-\lambda L_T.
\]

The count score is

\[
s(n)=\mathbb E(U_T\mid N_T=n),
\]

and

\[
I_N=\mathbb E\{s(N_T)^2\}.
\]

Estimate it by cross-fitted conditional-score regression over the integer count bins. Use independent simulation batches; for each held-out batch, estimate `s(n)` from the remaining batches. Pool only the two tails, deterministically from the outside inward, until every training bin has at least 500 frames. Subtract the estimated conditional-mean estimation variance from `s_hat(n)^2` before averaging.

Run an independent audit estimator from R14.B1:

\[
\widehat I_{\rm miss}
=
\overline N
-
\lambda^2
\sum_n\widehat p_n\widehat{\operatorname{Var}}(L_T\mid N_T=n).
\]

Require the two estimators to agree within

\[
\max\{1\%\text{ relative},\,2\text{ pooled standard errors}\}.
\]

If they do not, record `FI_ESTIMATOR_DISAGREEMENT`; do not choose the more favorable curve.

The deterministic `c=0` cells must additionally agree with the exact PMF table within 0.5% at every audited load.

## 5.2 Hold-time families

### Full lognormal sweep

Use

\[
c\in\{0,0.005,0.01,0.015,0.02,0.03,0.05,0.075,
0.10,0.15,0.20,0.30,0.50,1.00\}.
\]

### Matched-CV gamma sweep

Use gamma holds with shape `1/c^2` and scale `tau c^2` at

\[
c\in\{0.03,0.05,0.10,0.20,0.30,1.00\}.
\]

The `c=1` gamma cell is the exponential-hold test.

### Bounded-shape universality audit

Use the symmetric two-point law

\[
B=\tau(1\pm c)
\]

with equal probability at

\[
c\in\{0.05,0.10,0.30\}.
\]

This third family tests whether matched variance but different skewness/tails changes only higher-order terms.

## 5.3 Window lengths

Primary Monte Carlo:

\[
\boxed{\nu\in\{200,2000\}.}
\]

This gives a factor-10 window change, for which the deterministic ridge changes by `10^(1/3)=2.154...`, while jitter-dominated peaks should remain approximately fixed.

Add a limited long-window audit at

\[
\nu=20000,
\qquad
c\in\{0,0.01,0.02,0.05,0.10\}
\]

for lognormal holds only. This audit is not required for the first dev decision, but it cleanly separates finite-window bias from the long-window cubic.

The deterministic `c=0` exponent is measured from the existing exact-PMF ridge at

\[
\nu\in\{50,100,200,500,1000,2000,5000,10000\},
\]

not from Monte Carlo.

## 5.4 Load grid and peak extraction

Because the same exponential uniforms and hold draws can be reused across all loads of one `(H,c,nu,batch)` condition, use a dense common-random-number grid:

### Coarse grid

\[
\boxed{161\text{ log-spaced loads on }[0.15,48].}
\]

### Refinement grid

After the coarse maximum, evaluate

\[
41\text{ log-spaced loads on }
[0.65\widehat\rho_{\rm coarse},
 1.50\widehat\rho_{\rm coarse}].
\]

If the coarse maximum is at an endpoint, expand once to `[0.05,96]`; a second endpoint maximum is reported as `PEAK_NOT_BRACKETED`.

Fit a weighted cubic in `log rho` to the nine refinement points centered on the largest estimated `J`; select the interior local maximum. The fitted maximum must lie within one refinement-grid interval of the raw grid maximum. Otherwise use the raw maximum and mark `PEAK_FIT_UNSTABLE`.

Peak confidence intervals come from resampling whole independent simulation batches, rerunning the same extraction algorithm 2000 times.

## 5.5 Sample sizes and precision escalation

For every nonzero-CV condition:

1. **Stage 1:** 8 independent batches × 7,500 active-start frames = 60,000 frames.
2. **Stage 2:** for every condition entering an exponent or distribution-shape fit, add 24 batches × 7,500 frames, giving 240,000 total.
3. **Precision escalation:** if the 95% batch-bootstrap half-width of `log rho_peak` exceeds `0.025`, add exactly 16 further batches × 7,500 frames once, giving 360,000 total. No further discretionary sampling is allowed.

Use distinct RNG streams for hold times and exponential active waits, with common random numbers across `rho` only within the same batch. Report batch-level results so Monte-Carlo uncertainty is independently auditable.

## 5.6 Frozen exponent tests

### Test E1 — deterministic exponent

Regress

\[
\log\{\rho_*(\nu,0)+2/3\}
\]

on `log nu` using the exact-PMF ridge. Report the slope and its numerical-fit residual against `1/3`.

### Test E2 — jitter exponent

For each hold family, fit

\[
\log\rho_*=a-\gamma\log c
\]

using only cells satisfying

\[
12\nu c^2\ge10,
\qquad
c\le0.20,
\qquad
\rho_*>1.5,
\]

and having neither peak-bracketing nor estimator-disagreement flags. Use inverse-variance weights from the batch bootstrap. The target is

\[
\gamma=2/3.
\]

Do not estimate the exponent from the single ratio `c=0.05` versus `0.10` alone.

### Test E3 — long-window cubic

For every jitter-dominated cell report

\[
R_{\rm cubic}
=c^2\rho_*^2(1+2\rho_*)-1.
\]

The long-window theorem predicts `R_cubic -> 0` as `nu` grows.

### Test E4 — scaling collapse

Plot

\[
Y=\frac{\rho_*}{(6\nu)^{1/3}}
\]

against

\[
X=12\nu c^2.
\]

The matched prediction is

\[
Y=(1+X)^{-1/3}.
\]

Report the weighted RMS error in `log Y`, separately by hold family and by `nu`.

### Test E5 — distribution-shape universality

At matched `(c,nu)`, compare lognormal, gamma, and bounded holds. The leading finite-variance theorem predicts convergence to the same peak as `nu` increases and `c` decreases. Differences are reported as higher-order shape effects, not hidden by pooling families.

## 5.7 Hardware-relevant warning metric

For every `c`, report

\[
\mathcal L_{\rm detridge}(c,\nu)
=
1-
\frac{J(\rho_*(\nu,0),\nu,c)}
     {J(\rho_*(\nu,c),\nu,c)}.
\]

This directly measures the information lost by operating a jittered detector at the deterministic ridge. The current observation of roughly 60% loss at `c=0.05`, `nu=2000` is a strong dev signal, but its final uncertainty must come from the batch protocol above.

---

# 6. Final R14 decision

## Theory-spine decision

**GO with Candidate B as paper 2's theory spine**, using the following hierarchy:

1. **Theorem:** exact hidden-hold missing-information identity.
2. **Theorem:** long-window count-information rate and exact cubic optimum.
3. **Asymptotic theorem:** small-jitter `cv^(-2/3)` law with coefficient `2^(-1/3)`.
4. **Matched-asymptotic prediction until fully proved:** finite-window crossover scaling function.
5. **Framework theorem:** Candidate A generalized OED/KKT program.
6. **Conditional proposition:** Candidate C brightness-alignment residual bound.

The attractive but incorrect elements of the initial conjecture are frozen out:

- replace the literal `min` by the smooth crossover;
- replace fitted `c≈0.72` by the theoretical `2^(-1/3)` leading coefficient and the exact finite-CV cubic;
- do not call the exponent pair universal outside regular finite-variance hidden iid holds;
- do not claim every empirical ladder optimum is a master-program KKT point;
- do not give fixed flux a universal approximation ratio.

## One-sentence unification

> The recovery-time law determines a scalar count-information kernel; constrained OED allocates illumination directions and loads against that kernel; fixed global flux is the special case in which scene brightness already satisfies the kernel's water-filling KKT condition.

That is the simplest rigorous architecture connecting B, A, and C.