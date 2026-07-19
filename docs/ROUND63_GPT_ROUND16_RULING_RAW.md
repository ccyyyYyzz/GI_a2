# R16 — Final ruling: higher-order constant corrections

**Scope:** this ruling refines only the random-hold ridge theory from R14. It does not reopen the OE manuscript or modify the R13/R15 M1 freeze conditions.

**Evidence audited:** `docs/R14_PREREGISTERED_PREDICTIONS.md` at commit `bd731d5`, including the two independent 150k/400k Monte-Carlo peak estimates at `nu=2000`.

## Executive verdict

| Question | R16 ruling |
|---|---|
| Higher-order jitter-limited expansion | The exact long-window cubic gives further terms, but they are tiny. The next term after the previously quoted expansion is `-2^(2/3)c^(4/3)/324`; the `c^2` coefficient is zero. |
| Does lognormal skewness change the extensive coefficient? | **No.** The extensive `O(nu)` missing-information coefficient depends only on hold-time variance. For lognormal holds, the third central moment enters only the `O(1)` total-FI / `O(1/nu)` per-slot correction. Its induced peak shift is leftward but below `10^-4` at the four tested CVs. It cannot explain the replicated 5–8% displacement. |
| Composite crossover correction | The leading crossover law must receive a uniform constant correction. With `x=12 nu c^2`, the corrected matched expansion is `rho_* = [2c^2+1/(6nu)]^(-1/3) - (x+4)/[6(x+1)] + O(nu^(-1/3))`. This interpolates the deterministic `-2/3` and jitter-limited `-1/6` constants. |
| Corrected predictions at `nu=2000` | `rho_* ≈ 10.22, 5.65, 3.52, 2.17` for `c=0.02,0.05,0.10,0.20`. The correction absorbs part, but not all, of the replicated low bias. |
| Next suspect | The lognormal third moment is ruled out as the source. The **histogram/tail-floor FI estimator and peak interpolation** are the leading suspects. `dlog=0.01` should still be tested, but its Gaussian benchmark predicts only a 0.2–0.5% effect and the wrong sign for the observed left shift. |

---

# 1. Exact long-window expansion beyond R14

For fixed nonzero hold-time CV `c`, the R14 long-window count-information rate is

\[
J_\infty(\rho,c)=\frac{\rho}{(1+\rho)(1+c^2\rho^2)},
\]

and its unique maximizer solves

\[
c^2\rho_c^2(1+2\rho_c)=1.
\]

Let `epsilon=c^(2/3)`. Direct series solution of the cubic gives

\[
\boxed{
\rho_c=
2^{-1/3}c^{-2/3}
-\frac16
+\frac{2^{1/3}}{36}c^{2/3}
-\frac{2^{2/3}}{324}c^{4/3}
+0\,c^2
+\frac{2^{1/3}}{5832}c^{8/3}
+O(c^{10/3}).
}
\]

The newly exposed terms are far too small to create an 8% shift. For example, at `c=0.2` the `c^(4/3)` term is approximately `-5.7e-4` in `rho`; at `c=0.02` it is approximately `-2.7e-5`.

---

# 2. Where the lognormal third moment enters

## 2.1 Delayed-renewal notation

Set the mean hold time to one time unit. The first delay is

\[
D=X\sim\mathrm{Exp}(\rho),
\]

and ordinary inter-count cycles are

\[
C=B+X.
\]

Write

\[
\mathbb EB=1,\qquad \operatorname{Var}B=c^2,
\qquad h_3=\mathbb E[(B-1)^3].
\]

For a mean-one lognormal hold,

\[
\boxed{h_3=3c^4+c^6},
\]

while its standardized skewness is `(3+c^2)c`.

The leading renewal quantities are

\[
a(\rho)=\frac{\rho}{1+\rho},
\]

\[
b(\rho,c)=\frac{\rho(1+c^2\rho^2)}{(1+\rho)^3},
\]

where `E N_nu = nu a + O(1)` and `Var N_nu = nu b + O(1)`.

The asymptotic third-cumulant rate of the count is

\[
d=
\frac{3(c^2+\rho^{-2})^2}{(1+\rho^{-1})^5}
-
\frac{h_3+2\rho^{-3}}{(1+\rho^{-1})^4}.
\]

## 2.2 Finite-window mean and variance constants

Laplace expansion of the delayed-renewal first and second factorial moments gives

\[
\mathbb EN_\nu=\nu a+a_0+o(1),
\qquad
\operatorname{Var}N_\nu=\nu b+b_0+o(1),
\]

with

\[
\boxed{
a_0=\frac{(1+c^2)\rho^2}{2(1+\rho)^2},
}
\]

and

\[
\boxed{
b_0=
\frac{\rho^2}{12(1+\rho)^4}
\left[
15c^4\rho^2-12c^2\rho+18c^2
-8h_3\rho^2-8h_3\rho
+\rho^2+4\rho+18
\right].
}
\]

Thus the third moment first appears in a boundary/finite-window constant, not in the extensive variance rate `b`.

## 2.3 First finite-window Fisher correction

Let a dot denote differentiation with respect to `theta=log rho`. Under the standard local Edgeworth expansion for a nonlattice delayed renewal count,

\[
\frac{I_N(\theta;\nu)}{\nu}
=
\frac{\dot a^2}{b}
+\frac{C_1(\rho,c,h_3)}{\nu}
+o(\nu^{-1}),
\]

where

\[
\boxed{
C_1=
\frac{2\dot a\dot a_0}{b}
-
\frac{\dot a^2 b_0}{b^2}
+
\frac{(b\dot b-d\dot a)^2}{2b^4}.
}
\]

The final square is the non-Gaussian Edgeworth correction. It can be written compactly as

\[
\boxed{
\frac{(b\dot b-d\dot a)^2}{2b^4}
=
\frac{\rho^4(h_3\rho-2c^2)^2}
{2(1+c^2\rho^2)^4}.
}
\]

For deterministic holds (`c=h_3=0`), this formula reduces to

\[
C_1(\rho,0,0)
=-\frac{\rho^2(\rho^2+4\rho-6)}{12(1+\rho)^2}
=-\frac{\rho^2}{12}-\frac{\rho}{6}+O(1),
\]

recovering the finite-window boundary terms used in the deterministic ridge derivation.

## 2.4 Sign and magnitude of lognormal skewness

The extensive `O(nu)` information rate is exactly independent of `h_3`; only the mean and variance of `B` enter it. For a lognormal hold, the skewness contribution to `C_1` begins at

\[
\Delta C_1
=-2^{-2/3}c^{8/3}+O(c^{10/3})
\]

when evaluated on the jitter ridge. The corresponding peak displacement is

\[
\boxed{
\Delta\rho_{\rm skew}
=-\frac{11\,2^{2/3}}{12\nu}c^{4/3}
+O\!\left(\frac{c^2}{\nu}\right).
}
\]

It is leftward, but negligible. At `nu=2000`, the numerical shifts relative to a same-mean/same-variance symmetric hold are approximately

| `c` | `Delta rho_skew` |
|---:|---:|
| 0.02 | `-2.7e-6` |
| 0.05 | `-1.1e-5` |
| 0.10 | `-2.4e-5` |
| 0.20 | `-4.9e-5` |

Therefore the replicated end-point displacement is **not** a lognormal-third-moment correction to the extensive missing-information coefficient.

---

# 3. Uniform correction to the matched crossover

The leading R14 composite law was

\[
\rho_0(\nu,c)=
\left(2c^2+\frac{1}{6\nu}\right)^{-1/3}.
\]

To obtain the next uniform term, define

\[
t=\nu^{-1/3},
\qquad
x=12\nu c^2,
\qquad
\rho=t^{-1}y.
\]

Combining the exact long-window kernel, the delayed-renewal `C_1/nu` correction, and the known deterministic second boundary term gives the matched expansion

\[
J_{\nu,c}=1-tF_1(y,x)+t^2F_2(y,x)+O(t^3),
\]

with

\[
F_1=\frac1y+\frac{1+x}{12}y^2,
\]

\[
F_2=\frac1{y^2}+\frac{x-2}{12}y
+\frac{(1+x)^2}{144}y^4.
\]

The optimizer has

\[
y_0=\left(\frac6{1+x}\right)^{1/3},
\qquad
 y_1=-\frac{x+4}{6(x+1)}.
\]

Hence the corrected matched ridge is

\[
\boxed{
\rho_*(\nu,c)=
\left(2c^2+\frac1{6\nu}\right)^{-1/3}
-
\frac{12\nu c^2+4}{6(12\nu c^2+1)}
+O(\nu^{-1/3}).
}
\]

This constant correction has both required limits:

- `c=0`: `-(4/6)=-2/3`, the deterministic constant;
- `12 nu c^2 -> infinity`: `-1/6`, the jitter-limited cubic constant.

The lognormal third moment does not enter either `F_1` or `F_2`; it first appears in the next `O(t^3)` objective term.

## Error at `c=0.02`, `nu=2000`

Here

\[
x=12\nu c^2=9.6.
\]

The leading matched law gives

\[
\rho_0=10.4222.
\]

The uniform constant correction is

\[
-\frac{13.6}{63.6}=-0.21384,
\]

giving

\[
\rho_{\rm comp,2}=10.2083.
\]

The known next long-window terms add only about `+0.0026`. Direct numerical maximization of `J_inf+C_1/nu` plus the deterministic `rho^4/(144 nu^2)` term gives `10.216`.

Thus the leading matched formula had an internally predictable error of about **2%**, not 7.5%. After correction, the two observed peaks `9.645/9.610` remain approximately `5.6–5.9%` low.

The formal next remainder scale is `nu^(-1/3)=0.0794` in absolute load times a dimensionless coefficient. A coefficient of order one gives sub-percent relative error here. Explaining a further shift of about `-0.59` would require an unexpectedly large coefficient near `-7`, inconsistent with the exact long-window and deterministic limiting expansions.

---

# 4. Corrected numerical predictions at `nu=2000`

The R16 source-of-record prediction is the maximizer of

\[
J_{\rm R16}(\rho)=J_\infty(\rho,c)
+\frac{C_1(\rho,c,3c^4+c^6)}{\nu}
+\frac{\rho^4}{144\nu^2},
\]

with the final term retained because it contributes at the same order in the crossover scaling.

| `c` | R16 corrected `rho_*` | measured 150k | measured 400k | pooled residual vs prediction |
|---:|---:|---:|---:|---:|
| 0.02 | **10.22** | 9.645 | 9.610 | `-5.7%` |
| 0.05 | **5.65** | 5.154 | 5.700 | `-3.9%` |
| 0.10 | **3.52** | 3.752 | 3.399 | `+1.6%` |
| 0.20 | **2.17** | 2.106 | 2.051 | `-4.0%` |

The exponent result survives unchanged. The higher-order theory absorbs roughly two to four percentage points of the original end-point discrepancy, but it does not plausibly absorb all of it.

Frozen interpretation:

> The `-2/3` CV exponent is supported. The leading coefficient is consistent only to several percent under the current histogram-Fisher estimator; lognormal skewness does not account for the remaining displacement.

---

# 5. Estimator diagnosis

## 5.1 Which suspect comes first?

The primary suspect is the **histogram tail floor / sparse-tail derivative estimator**, followed by the local quadratic peak interpolation. The finite-difference step `dlog=0.01` is a secondary suspect.

At the corrected peaks, a `0.01` perturbation moves the asymptotic Gaussian mean by only

| `c` | standardized displacement `q` |
|---:|---:|
| 0.02 | 0.418 |
| 0.05 | 0.397 |
| 0.10 | 0.372 |
| 0.20 | 0.339 |

For a pure Gaussian location family, central differencing of the density has FI ratio

\[
\frac{I_h}{I}=\frac{\sinh(q^2)}{q^2},
\]

which gives only `+0.51%, +0.41%, +0.32%, +0.22%`, respectively. This is too small and tends in the wrong direction to explain a persistent leftward 5–6% peak shift.

By contrast, replacing rare-bin probabilities by `0.5/N` directly suppresses or distorts high-score tail contributions. Its effect is load dependent and can move a broad maximum. The replicated low end-point values are consistent with such a systematic estimator effect.

## 5.2 Decisive diagnostic: use the exact complete-score identity

The simulator already knows each frame's active time `L`. For `tau=1`, define the complete-path score

\[
U=N-\rho L.
\]

The count-only Fisher information is exactly

\[
I_N=\operatorname{Var}\{\mathbb E(U\mid N)\}
=\mathbb EN-\rho^2\mathbb E\{\operatorname{Var}(L\mid N)\}.
\]

This yields a Monte-Carlo estimator that uses **one load only**, with no finite difference and no PMF floor.

For bins with `k_n>=2`, use

\[
\widehat I_{\rm MI}
=\bar N-ho^2
\sum_n\frac{k_n}{N_{\rm MC}}s^2_{L\mid n},
\]

where `s^2` is the unbiased within-bin variance. As a cross-check, use the bias-corrected conditional-score form

\[
\widehat I_{\rm CS}
=\sum_n\frac{k_n}{N_{\rm MC}}
\left(\bar U_n^2-\frac{s^2_{U\mid n}}{k_n}\right).
\]

Pool only the extreme left and right tails whose individual bin counts are below 50, and report their total probability. The two estimators must agree within Monte-Carlo uncertainty.

## 5.3 Frozen diagnostic run

Run only the two replicated end points first:

- `nu=2000`;
- lognormal `c in {0.02,0.20}`;
- `2,000,000` frames per load;
- common load grids:
  - `c=0.02`: 31 points on `[8.5,11.5]`;
  - `c=0.20`: 29 points on `[1.75,2.45]`;
- record `(N,L)` for every frame;
- 200 block-bootstrap replicates for the peak location;
- peak extracted from a shape-preserving cubic interpolation and, separately, a three-point quadratic; both reported.

On the same simulated paths compute the existing histogram estimator for

\[
d\log\rho\in\{0.0025,0.005,0.01,0.02\}
\]

and tail pseudocount/floor levels

\[
\alpha/N,\qquad \alpha\in\{0,0.05,0.10,0.50,1.00\}.
\]

Also run the identical estimator on deterministic holds (`c=0`) around the exact `rho_*=22.2543` as an estimator calibration.

## 5.4 Decision rule

- If `I_MI/I_CS` peak near `10.22` and `2.17`, while the histogram peak moves with `alpha` or sample size, the residual is a tail-floor artifact.
- If the histogram peak moves quadratically with `dlog^2` but is stable to `alpha`, extrapolate to `dlog=0` and attribute it to finite differencing.
- If the score-identity estimator independently reproduces `9.62` and `2.08`, the discrepancy is physical/theoretical; then the next task is a genuinely uniform second-order renewal Edgeworth calculation, not another fitted constant.
- No new constant is to be fitted from the present eight peak values before this diagnostic.

---

# 6. Frozen claim discipline

Permitted after the current evidence:

> The jitter-limited load follows the predicted `c^(-2/3)` scaling; a uniform next-order correction interpolates the deterministic `-2/3` and random-hold `-1/6` constants.

> Hold-time skewness enters only the finite-window correction and is too small to explain the remaining Monte-Carlo peak displacement.

Not permitted yet:

- a new empirical coefficient replacing `2^(-1/3)`;
- a lognormal-specific `4/3` multiplier;
- claiming the residual proves a new universality subclass;
- treating histogram-Fisher peak locations as exact without the score-identity diagnostic.
