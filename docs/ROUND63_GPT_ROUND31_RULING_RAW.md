# R31 ruling (GitHub issue #23, raw)

PROVENANCE NOTE (2026-07-22): this issue body is the EXTRA-HIGH-tier
response (the first R31 submission). The PRO-tier response arrived as an
independent audit ADDENDUM (comment 0 on #23, archived as
ROUND63_GPT_ROUND31_PRO_APPENDIX_0.md) which confirms this ruling as
"substantively correct" and adds binding amendments: the correlated-
speckle determinant law, the mask/power-dependent C_eff definition, the
irreducible fixed-speckle floor, the n_eff^2 photon-dose law, fence
additions, and a narrowed GO. WHERE THEY DIFFER, THE PRO ADDENDUM
GOVERNS (operator rule).

Title: R31 theory ruling — conditional GO for microcell-occupancy Laplace tomography
Posted: 2026-07-22T12:25:33Z

---

# R31 theory ruling — conditional GO for microcell-occupancy Laplace tomography

**Reference request:** `docs/ROUND63_GPT_ROUND31_QUESTION.md` at commit [`b962ae903ae342fc57a90d94188221b912f8b5f0`](https://github.com/ccyyyYyzz/GI_a2/commit/b962ae903ae342fc57a90d94188221b912f8b5f0)

## Executive verdict

> **CONDITIONAL GO — the corrected theorem is real, but the proposed physical premise and novelty claim need two severe repairs before this can be called the flagship.**

The mathematical core survives in a narrower and cleaner form. A bank of binary microcells, operated with a controlled attenuation/power sweep, samples a **Laplace transform of the per-cell optical-load distribution**. If, and only if, the masked scene reaches each microcell as an **intensity-additive sum of independent positive modes**, the no-fire curve factorizes and its log-curvature exposes the masked power sums
\[
p_j(m\odot x)=\sum_p (m_p x_p)^j.
\]
For binary masks, \(p_1=Mx\) and \(p_2=M(x^{\circ 2})\). This second channel is generically outside the linear bucket span and separates almost every *fixed* linear-null pair.

Three corrections are non-negotiable:

1. **“Fully developed speckle” is not enough.** With ordinary coherent field addition,
   \[
   E_c=\sum_p h_{cp}\sqrt{v_p},\qquad h_{cp}\sim\mathcal{CN}(0,1),
   \]
   one has \(|E_c|^2\sim{\rm Exp}(\sum_p v_p)\), hence
   \[
   r(\rho)=P\{\text{no fire}\}=\frac{1}{1+\rho\sum_pv_p}.
   \]
   The response then depends only on the ordinary bucket and the entire \(x^2\) claim dies. The product formula requires incoherent/intensity-level addition, independent temporal/frequency/polarization modes, fluorescence-like emission, or another explicit mechanism that destroys the coherent cross terms **before** the microcell saturation.

2. **The stable observable is \(p_2\), not \(e_2\).** Although \(r^{-1}\) is an elementary-symmetric polynomial in the ideal exponential model, recovering
   \[
   p_2=p_1^2-2e_2
   \]
   is a cancellation problem when many pixels contribute. The production estimator must fit the curvature of \(\log r(\rho)\) directly.

3. **At 5–10% mask rate this is not exact dense-scene recovery.** For \(K\) binary masks in \(N\) pixels, the generic local rank is only
   \[
   \operatorname{rank}\!\begin{bmatrix}M\\M\operatorname{diag}(x)\end{bmatrix}
   =\min(N,2K).
   \]
   Thus the linear fiber dimension falls from \(N-K\) to \(N-2K\), but remains \(0.8N\)–\(0.9N\) at 10–5% sampling. Exact recovery at those rates requires a sparse or low-dimensional image prior; it is not supplied by saturation alone.

The prior-art verdict is equally narrow. Swept-efficiency on/off detection as a generating-function measurement, multiplexed click counting, occupancy saturation, and higher-order optical correlations are all established. The potentially publishable intersection is:

> **deliberately using an unaddressed SiPM's spatially heterogeneous microcell occupancy, under mask-indexed power sweeps, to obtain a second image equation \(M(x^{\circ2})\) and thereby reduce the ordinary single-pixel/ghost-imaging fiber.**

I found no exact prior instance of that complete spatial inverse construction. That is a **novelty candidate**, not yet a certified novelty claim.

My authorization is therefore:

- proceed with T1–T3 and one physically faithful forward-model probe;
- freeze \(p_1,p_2\) as the only primary information channels;
- demote \(p_3+\) to a stress-test/negative boundary;
- do not call this “zero-calibration”: it may avoid per-cell calibration, but still needs global PDE/power, dark-count, recovery, crosstalk, effective-mode-number, and \(C_{\rm eff}\) calibration;
- do not call it a universal null-space jailbreak unless the restricted-prior Jacobian is shown to close.

---

# 1. Corrected theorem seed

## 1.1 Detector and optical assumptions

Let \(v=m\odot x\in\mathbb R_+^N\). For microcell \(c=1,\ldots,C\), define the optical load

\[
T_c(v)=\sum_{p=1}^{N}v_pW_{cp}.
\]

For a gate at scalar power/efficiency \(\rho\), assume:

- **A1 — intensity addition:** the photon intensity at a cell is \(T_c(v)\), not the squared modulus of a coherent sum;
- **A2 — Poisson photo-arrivals:** conditional on \(T_c\), the number of primary avalanches is Poisson with mean \(\rho T_c\);
- **A3 — one binary opportunity per cell:** a cell records \(B_c=1\) iff at least one primary avalanche occurs during a gate; no refiring within the integration window;
- **A4 — no unmodelled coupling:** prompt/delayed crosstalk, afterpulsing, dark counts, recovery, PDE heterogeneity and readout saturation are absent or explicitly included as nuisance parameters;
- **A5 — exchangeable cells:** the row vectors \(W_c=(W_{c1},\ldots,W_{cN})\) are iid, or at least stationary/ergodic with a known effective sample size;
- **A6 — power scale:** \(\rho\) is known up to a separately identifiable global factor.

Write
\[
A_c(\rho)=1-B_c,\qquad
r(\rho;v)=E[A_c(\rho)]
\]
for the no-fire probability.

## 1.2 General Laplace theorem

### Theorem 1 — Masked microcell Laplace response

Under A1–A5,

\[
\boxed{
r(\rho;v)
=E_W\!\left[e^{-\rho W^\top v}\right]
=\mathcal L_W(\rho v),
}
\]
the multivariate Laplace transform of the positive random vector \(W\) evaluated along the ray \(\rho v\).

If the components \(W_p\) are independent, with marginal Laplace transforms
\[
L_p(s)=E[e^{-sW_p}],
\]
then
\[
\boxed{
r(\rho;v)=\prod_{p=1}^{N}L_p(\rho v_p).
}
\]

**Proof.** Conditional on \(W_c\), a cell does not fire iff its Poisson primary count is zero, so
\[
P(A_c=1\mid W_c)=e^{-\rho W_c^\top v}.
\]
Averaging gives the multivariate Laplace transform. Component independence factorizes the expectation. \(\square\)

### Corollary 1 — exponential one-mode law

If \(W_{cp}\stackrel{\rm iid}{\sim}{\rm Exp}(1)\), then
\[
\boxed{
r(\rho;v)=\prod_p(1+\rho v_p)^{-1},
\qquad
r(\rho;v)^{-1}=\prod_p(1+\rho v_p)
=\sum_{j=0}^N e_j(v)\rho^j.
}
\]

Thus the request's algebraic seed is correct under A1–A6.

### Corollary 2 — Gamma/multimode law

If
\[
W_{cp}\stackrel{\rm iid}{\sim}{\rm Gamma}(k,\text{scale}=1/k)
\]
with mean one and known \(k\), then
\[
\boxed{
r(\rho;v)=\prod_p\left(1+\frac{\rho v_p}{k}\right)^{-k}.
}
\]

The coefficients remain identifiable in principle, but the \(j\)th power-sum signal is suppressed by \(k^{-(j-1)}\). Unknown \(k\) is not benign: at second order it confounds \(p_2\) with \(k\).

## 1.3 The coherent-speckle counterexample

Suppose instead that the detector-cell field is a coherent sum
\[
E_c=\sum_p h_{cp}\sqrt{v_p},
\qquad h_{cp}\stackrel{\rm iid}{\sim}\mathcal{CN}(0,1).
\]
Then
\[
E_c\sim\mathcal{CN}\!\left(0,\sum_pv_p\right),
\]
so \(T_c=|E_c|^2\) is exponential with mean \(p_1=\sum_pv_p\). Therefore
\[
\boxed{
r(\rho;v)=\frac{1}{1+\rho p_1}.
}
\]

Every \(v\) with the same ordinary bucket has the same entire power-sweep curve. This is an exact kill, not a small correction.

The first experimental/simulation gate is therefore not “does the product fit better than Monte Carlo noise?” It is:

\[
H_{\rm product}:\quad
r(\rho;v)=\prod_p(1+\rho v_p)^{-1}
\]
versus
\[
H_{\rm coherent}:\quad
r(\rho;v)=(1+\rho p_1)^{-1},
\]
plus an intermediate correlated/multimode model. The comparison must use scenes with identical \(p_1\) but strongly different \(p_2\).

## 1.4 Power sums are the statistically natural coordinates

For the one-mode exponential law,

\[
\log r(\rho)
=-\sum_p\log(1+\rho v_p)
=\sum_{j\ge1}\frac{(-1)^j}{j}\rho^j p_j(v),
\]
where
\[
p_j(v)=\sum_pv_p^j.
\]

Hence
\[
\boxed{
p_j(v)=
\frac{(-1)^j}{(j-1)!}
\left.\frac{d^j}{d\rho^j}\log r(\rho)\right|_{\rho=0}.
}
\]

For Gamma shape \(k\),
\[
\log r(\rho)
=\sum_{j\ge1}\frac{(-1)^j}{j\,k^{j-1}}\rho^jp_j(v).
\]

The first identities are
\[
e_1=p_1,\qquad
e_2=\frac{p_1^2-p_2}{2},\qquad
e_3=\frac{p_1^3-3p_1p_2+2p_3}{6}.
\]

For a diffuse mask define
\[
n_{\rm eff}=\frac{p_1^2}{p_2}.
\]
Then \(p_2/p_1^2=1/n_{\rm eff}\), while \(e_2/p_1^2\approx1/2\). Extracting \(p_2\) by subtracting two order-one quantities is ill-conditioned by approximately \(n_{\rm eff}\). Therefore:

> **Fit \(\log r\) in power-sum coordinates. Do not estimate \(e_2\) and subtract.**

---

# 2. Quenched self-averaging and whether \(C=3600\) is enough

Let a fixed detector/speckle realization have

\[
r_C(\rho)=\frac1C\sum_{c=1}^C e^{-\rho T_c}.
\]

This is the curve approached after arbitrarily many repeated gates on the same fixed \(W\).

## 2.1 Exact finite-\(C\) covariance

For iid cells,

\[
E[r_C(\rho)]=r(\rho),
\]
and
\[
\boxed{
\operatorname{Cov}\!\left[r_C(\rho_i),r_C(\rho_j)\right]
=
\frac{r(\rho_i+\rho_j)-r(\rho_i)r(\rho_j)}{C}.
}
\]

In particular,
\[
\boxed{
\frac{\operatorname{sd}[r_C(\rho)]}{r(\rho)}
=
\frac1{\sqrt C}
\left[
\frac{r(2\rho)}{r(\rho)^2}-1
\right]^{1/2}.
}
\]

For bounded \(Y_c=e^{-\rho T_c}\in[0,1]\), Bernstein's inequality gives, with probability at least \(1-\delta\),

\[
|r_C-r|
\le
\sqrt{\frac{2\sigma_\rho^2\log(2/\delta)}{C}}
+
\frac{2\log(2/\delta)}{3C},
\]
where
\[
\sigma_\rho^2=r(2\rho)-r(\rho)^2.
\]
A union bound gives the same control over a finite sweep grid.

If cells are correlated, replace \(C\) by an empirically justified
\[
C_{\rm eff}
=
\frac{C}{1+\text{design-effect from optical/electrical correlations}}.
\]
The physical number \(C=3600\) is irrelevant if a speckle grain covers many microcells or crosstalk creates cell clusters.

## 2.2 Shot noise versus the quenched floor

Suppose level \(\rho_j\) is repeated for \(G_j\) gates on the same \(W\), and
\[
\widehat r_j=1-\frac{1}{CG_j}\sum_{g=1}^{G_j}S_{jg}.
\]
Then

\[
\boxed{
\operatorname{Var}(\widehat r_j)
=
\frac{r(2\rho_j)-r_j^2}{C}
+
\frac{r_j-r(2\rho_j)}{CG_j}.
}
\]

For different levels on the same fixed \(W\),
\[
\operatorname{Cov}(\widehat r_i,\widehat r_j)
=
\frac{r(\rho_i+\rho_j)-r_ir_j}{C}
\qquad(i\ne j),
\]
assuming independent photo-arrival noise across gates.

At \(G=1\), the two terms sum to the ordinary binomial variance \(r(1-r)/C\). As \(G\to\infty\), photon noise vanishes but the first, quenched term remains. Repeating pulses cannot calibrate away a finite random microcell ensemble.

## 2.3 Numerical scale at the useful operating point

Let \(t=\rho p_1\). In a diffuse scene with \(n_{\rm eff}\gg1\),

\[
\frac{r(2\rho)}{r(\rho)^2}
\approx
\exp\!\left(\frac{t^2}{n_{\rm eff}}\right).
\]

At \(C=3600\), \(t=3.5\), \(n_{\rm eff}=200\),

\[
\frac{\operatorname{sd}(r_C)}{r}\approx0.42\%.
\]

For a one-component scene at the same \(t\),

\[
r=(1+t)^{-1},\qquad
\frac{\operatorname{sd}(r_C)}{r}
=
\frac1{\sqrt C}
\sqrt{\frac{(1+t)^2}{1+2t}-1}
\approx2.1\%.
\]

A second, useful view is that \(p_2=\operatorname{Var}(T)\) in the exponential-component model. The asymptotic relative standard deviation of the sample variance across \(C\) cells is

\[
\boxed{
\frac{\operatorname{sd}(\widehat p_2)}{p_2}
\approx
\frac1{\sqrt C}
\sqrt{2+6\frac{p_4}{p_2^2}}.
}
\]

At \(C=3600\), this is about \(2.36\%\) for many equal contributors and \(4.71\%\) for one dominant contributor.

**Ruling on \(C=3600\):** it is adequate for a percent-level instrument curve under near-independent cells, but not a theorem that calibration is negligible. The campaign must estimate \(C_{\rm eff}\) and the effective mode number \(k\); otherwise the advertised self-averaging rate is unverified.

---

# 3. Exact and operational Fisher geometry

There are two different likelihoods. They must not be conflated.

## 3.1 Quenched, fixed-\(W\) Poisson-binomial likelihood

Let
\[
a_{cj}=e^{-\rho_jT_c},\qquad q_{cj}=1-a_{cj}.
\]
For one aggregate count \(S_j\),

\[
G_j(z;\theta)
=
\prod_{c=1}^C[a_{cj}+(1-a_{cj})z],
\qquad
P_j(s;\theta)=[z^s]G_j(z;\theta).
\]

The exact Fisher information is

\[
\boxed{
\mathcal I_{ab}^{\rm PB}
=
\sum_jG_j^{\rm rep}
\sum_{s=0}^C
\frac{\partial_aP_j(s)\,\partial_bP_j(s)}{P_j(s)},
}
\]
with
\[
\partial_aP_j(s)
=
[z^s]\,
G_j(z)
\sum_c
\frac{(1-z)\partial_a a_{cj}}
{a_{cj}+(1-a_{cj})z}.
\]

This is exact conditional on known \(T_c\). Aggregation loses information relative to labeled microcell outputs:
\[
\mathcal I(S)\preceq\sum_c\mathcal I(B_c),
\]
because the aggregate score is the conditional expectation of the labeled score given \(S\).

A crucial limitation follows: for a fixed unknown \(W\), a truncated vector \((e_1,e_2)\) does **not** determine the individual \(q_c\). Thus “the exact Poisson-binomial FIM of \(e_1,e_2\) with no cell calibration” is not a closed model. One must either:

- refresh \(W\) between gates and integrate it out;
- adopt the ensemble likelihood with the finite-\(C\) covariance above;
- or introduce a calibrated/random-effects model for \(W\).

## 3.2 Annealed/exchangeable binomial likelihood

If the cell weights are independently refreshed, or the ensemble approximation is accepted, then
\[
S_j\sim{\rm Binomial}(C,1-r_j).
\]
For any parameter vector \(\theta\),

\[
\boxed{
\mathcal I_{ab}
=
C\sum_jG_j^{\rm rep}
\frac{r_j}{1-r_j}
\,
\partial_a\log r_j\,
\partial_b\log r_j.
}
\]

In the one-mode power-sum expansion,
\[
\frac{\partial\log r}{\partial p_\ell}
=
\frac{(-1)^\ell\rho^\ell}{\ell},
\]
so the truncated information matrix is

\[
\boxed{
\mathcal I_{\ell m}
=
C\sum_jG_j^{\rm rep}
\frac{r_j}{1-r_j}
\frac{(-1)^{\ell+m}\rho_j^{\ell+m}}{\ell m}.
}
\]

In elementary-symmetric coordinates, with
\[
P(\rho)=\sum_ke_k\rho^k,\qquad r=P^{-1},
\]
\[
\frac{\partial\log r}{\partial e_k}=-\rho^kr,
\]
hence
\[
\mathcal I_{ab}^{(e)}
=
C\sum_jG_j^{\rm rep}
\frac{\rho_j^{a+b}r_j^3}{1-r_j}.
\]
This representation makes the high-order ill-conditioning visible.

## 3.3 Nuisance geometry

At minimum include:

- global power/PDE scale \(\alpha\), replacing \(\rho\to\alpha\rho\);
- dark no-fire factor, e.g. \(r_{\rm obs}=e^{-d}r\);
- Gamma/effective-mode parameter \(k\);
- crosstalk and afterpulse parameters;
- recovery/refiring parameter determined by pulse width relative to recovery time;
- an over-saturation/readout parameter if analog charge can exceed \(C\).

Without a low-power anchor or external power reference, \(\alpha\) and \(p_1\) are confounded. Without identifying \(k\), the second-order coefficient identifies \(p_2/k\), not \(p_2\).

The production Fisher matrix must be the Schur complement after nuisance elimination, not the optimistic \(p_2\)-only diagonal entry.

---

# 4. Optimal power sweep and the saturation ridge

Use dimensionless load
\[
t=\rho p_1,\qquad
\beta_2=\frac{p_2}{p_1^2}=\frac1{n_{\rm eff}}.
\]
For diffuse one-mode inputs,
\[
\log r=-t+\frac{\beta_2t^2}{2}+O(\beta_3t^3).
\]

## 4.1 If \(p_1\) is known

At \(\beta_2=0\), the ideal per-gate information for \(\beta_2\) is

\[
\boxed{
I_{\beta_2\beta_2}^{\rm gate}(t)
=
\frac{Ct^4}{4(e^t-1)}.
}
\]

Incident photons per gate are proportional to \(Ct\), so information per incident photon is proportional to
\[
\frac{t^3}{e^t-1}.
\]
The optimum solves
\[
\boxed{
t=3(1-e^{-t}),
\qquad
t_\star=2.821439\ldots
}
\]
with no-fire fraction
\[
r_\star\simeq e^{-t_\star}=0.0595.
\]

Thus the useful ridge is **moderate/deep occupancy with 6% cells still dark**, not the limit \(S/C\to1\).

If gate count rather than incident photons were the resource, the optimum solves
\[
t=4(1-e^{-t}),
\qquad
t=3.920690\ldots,
\]
but this is not the fair primary design because it spends more incident photons.

## 4.2 If \(p_1\) is nuisance: local two-point c-optimal design

Parameterize by \((\log p_1,\beta_2)\). Per incident photon, the local information contribution is proportional to

\[
M(t)
=
\frac1{e^t-1}
\begin{bmatrix}
t & -t^2/2\\
-t^2/2 & t^3/4
\end{bmatrix}.
\]

The numerical two-point c-optimal design for \(\beta_2\) is:

\[
\boxed{
t_L\to0,\qquad t_H=3.54403\ldots
}
\]
with incident-photon allocation approximately
\[
\boxed{
24.5\%\ \text{to the low-power anchor},\qquad
75.5\%\ \text{to the high-curvature level}.
}
\]

At the high point,
\[
r_H\approx e^{-3.544}=0.0289,
\qquad Cr_H\approx104
\quad(C=3600).
\]

The \(t_L\to0\) limit is an ideal local-design result, not a literal zero-light exposure. Use the smallest level that gives reliable power/dark calibration.

## 4.3 Recommended robust three-level sweep

The exact two-point design is fragile to \(k\), \(p_3\), dark counts and model error. Freeze a three-level discovery sweep:

1. **anchor:** \(t\in[0.05,0.2]\), for \(p_1\), scale and dark response;
2. **audit:** \(t\approx1\), to detect link-model failure;
3. **curvature ridge:** \(t\in[3.0,3.6]\), for \(p_2\).

Allocate most incident photons to the curvature level, but require at least 20–25% to the anchor until nuisance identifiability is demonstrated.

Hard guards:

- predicted and observed no-fire count \(Cr\ge50\), preferably \(\ge100\);
- no analog over-saturation/refiring evidence;
- fitted \(k\) and \(C_{\rm eff}\) stable across masks;
- leave-one-level-out prediction error below a predeclared tolerance;
- the product model must beat the coherent-collapse model on equal-\(p_1\), unequal-\(p_2\) controls.

## 4.4 Which orders are physically identifiable?

For normalized
\[
\beta_j=\frac{p_j}{p_1^j},
\]
and lower moments known, the ideal per-photon information for \(\beta_j\) is proportional to
\[
\frac{t^{2j-1}}{e^t-1}.
\]
Its optimum solves
\[
t=(2j-1)(1-e^{-t}).
\]

Therefore:

| order | per-photon optimum \(t\) | \(Cr\) at \(C=3600\) | ruling |
|---|---:|---:|---|
| \(p_2\) | 2.821 | 214 | primary, plausible |
| \(p_3\) | 4.965 | 25 | stress test only |
| \(p_4\) | 6.994 | 3.3 | practically nonregular |

For \(n_{\rm eff}=200\), \(\beta_2=0.005\). At \(t=3.54\),
\[
I_{\beta_2\beta_2}^{\rm gate}\approx4225,
\]
so an optimistic known-\(p_1\) lower bound for \(d'=3\) against zero curvature is about
\[
G\approx\frac{9}{I\beta_2^2}\approx85
\]
high-level gates per mask. At \(n_{\rm eff}=100\), the corresponding lower bound is about 21 gates.

For equal contributors, \(\beta_3=1/n_{\rm eff}^2\). At \(n_{\rm eff}=200\), even the ideal \(p_3\) ridge needs on the order of \(3\times10^5\) gates for \(d'=3\). Unknown nuisances make it worse.

For Gamma shape \(k\), the \(p_2\) effect is divided by \(k\), so gate demand scales approximately as \(k^2\); \(p_3\) demand scales as \(k^4\). A modest \(k=10\) turns the 85-gate ideal \(p_2\) example into roughly \(8.5\times10^3\) gates.

> **Freeze \(p_1,p_2\). Kill \(p_3+\) as a campaign claim unless the actual \(n_{\rm eff}\), \(k\), and measured Fisher information contradict this bound.**

---

# 5. Downstream identifiability: what the extra channel really buys

## 5.1 Binary-mask measurement map

For \(M\in\{0,1\}^{K\times N}\),

\[
y=Mx,\qquad z=M(x^{\circ2}).
\]

Define
\[
F_M(x)=
\begin{bmatrix}
Mx\\
M(x^{\circ2})
\end{bmatrix}.
\]

The Jacobian is

\[
\boxed{
J_F(x)=
\begin{bmatrix}
M\\
2M\operatorname{diag}(x)
\end{bmatrix}.
}
\]

### Theorem 2 — generic local rank

For continuous generic \(M\) and generic nonzero \(x\),

\[
\boxed{
\operatorname{rank}J_F(x)=\min(N,2K).
}
\]

**Proof witness.** Choose distinct positive \(x_j\), and let
\[
M_{ij}=x_j^{2(i-1)},\qquad i=1,\ldots,K.
\]
The top block has powers \(0,2,\ldots,2K-2\); the bottom block has powers \(1,3,\ldots,2K-1\). After row permutation, any \(2K\)-column minor is a Vandermonde matrix. Hence one rank-\(\min(N,2K)\) instance exists, and the corresponding determinant polynomial is nonzero; the deficient set has measure zero. \(\square\)

For a fixed binary \(M\), generic rank in \(x\) is governed by the union of two copies of the column matroid of \(M\). Random Bernoulli masks are expected to attain \(\min(N,2K)\) in the relevant regime, but the exact production statement should be a numerical rank/smallest-singular-value certificate for the realized mask bank.

### Corollary — local fiber dimension

At a regular interior point,
\[
\dim F_M^{-1}(F_M(x))
=
N-\operatorname{rank}J_F(x)
=
\max(0,N-2K).
\]

Thus:

| mask rate \(K/N\) | linear fiber \(N-K\) | joint local fiber \(N-2K\) |
|---:|---:|---:|
| 5% | 95% of \(N\) | 90% of \(N\) |
| 10% | 90% of \(N\) | 80% of \(N\) |

This is a real factor-of-two gain in local measurement rank, not exact dense-scene recovery.

## 5.2 Null-pair separation and its limit

Let \(x'=x+h\) with \(Mh=0\). The second channel changes by

\[
\boxed{
M[(x+h)^{\circ2}-x^{\circ2}]
=
M(2x\odot h+h^{\circ2}).
}
\]

For a **fixed nonzero** \(h\in\ker M\), this is generically nonzero as a function of \(x\). So almost every chosen linear-null pair is separated.

But when \(2K<N\), there remains a positive-dimensional family of \(h\)'s satisfying both equations. Therefore both statements are simultaneously true:

- a random prechosen null pair is usually separated;
- the full nonlinear map is still noninjective at 5–10% for arbitrary dense \(x\).

T3 must not stop after finding one separable null pair. It should estimate the smallest Fisher-weighted singular values of the restricted Jacobian and search adversarially for residual joint-fiber pairs.

## 5.3 Exact recovery under structure

Let the admissible images form a \(d\)-dimensional smooth manifold \(\mathcal M\), with tangent basis \(B_x\in\mathbb R^{N\times d}\). Local identifiability on \(\mathcal M\) holds iff

\[
\boxed{
\operatorname{rank}
\begin{bmatrix}
MB_x\\
2M\operatorname{diag}(x)B_x
\end{bmatrix}
=d.
}
\]

A generic dimension threshold is \(2K\ge d\). For a known support of size \(s\), replace \(d\) by \(s\); local recovery can therefore begin near \(2K\ge s\).

This is the mathematically honest 5–10% route: the extra channel can close a low-dimensional prior tangent space that the linear bucket alone leaves open. Global uniqueness still requires a separate secant/no-collision theorem or empirical adversarial certificate. Local rank alone is not a global recovery theorem.

## 5.4 Mask details that matter

- For gray nonnegative masks, the second channel is
  \[
  z=M^{\circ2}(x^{\circ2}),
  \]
  not \(M(x^{\circ2})\).
- A formal \(\pm1\) mask gives \(m_p^2=1\), so one nonlinear exposure has no spatially varying \(p_2\) weights. Physical complementary \(0/1\) exposures can recover a signed \(x^2\) equation by subtracting two separately estimated \(p_2\)'s, at double exposure cost.
- Saturated complement counts cannot be subtracted before inversion as though the detector were linear.
- Reconstruction should use the joint count likelihood. A two-stage estimate \((\widehat y,\widehat z)\) is useful diagnostically, but production should not square a noisy estimate of \(x\).

A suitable estimator is
\[
\widehat x
=
\arg\min_{x\ge0}
\left\{
-\log L(\{S_{kj}\}\mid x)
+\lambda R(x)
\right\},
\]
with the full sweep link and frozen nuisance model.

## 5.5 Fisher-weighted null-pair discriminability

For two scenes with the same linear buckets, under the annealed model,

\[
\boxed{
d'^2
\approx
\sum_{k,j}
CG_{kj}\frac{r_{kj}}{1-r_{kj}}
\left[
\log r_{kj}(x')-\log r_{kj}(x)
\right]^2.
}
\]

At leading second order,
\[
\Delta\log r_{kj}
\approx
\frac{\rho_j^2}{2k_{\rm eff}}
\Delta p_{2,k}.
\]

This is the correct analytic target for T3. The gate is not merely \(d'\ge3\) on one favorable pair; it must hold over a frozen bank stratified by \(n_{\rm eff}\), \(\Delta p_2\), mask brightness, and nuisance uncertainty.

---

# 6. Ruthless prior-art fence

## 6.1 Swept-efficiency on/off photon statistics — direct theorem collision

The use of no-click probabilities measured at multiple efficiencies/attenuations to infer photon-number statistics is established:

- Rossi, Olivares & Paris, “Photon statistics without counting photons,” *Phys. Rev. A* 70, 055801 (2004), [DOI 10.1103/PhysRevA.70.055801](https://doi.org/10.1103/PhysRevA.70.055801).
- Zambra et al., “Experimental reconstruction of photon statistics without photon counting,” *Phys. Rev. Lett.* 95, 063602 (2005), [DOI 10.1103/PhysRevLett.95.063602](https://doi.org/10.1103/PhysRevLett.95.063602).
- Sperling, Vogel & Agarwal, “True photocounting statistics of multiple on-off detectors,” *Phys. Rev. A* 85, 023820 (2012), [DOI 10.1103/PhysRevA.85.023820](https://doi.org/10.1103/PhysRevA.85.023820).
- Miatto, Safari & Boyd, “Explicit formulas for photon number discrimination with on/off detectors,” *Appl. Opt.* 57, 6750 (2018), [DOI 10.1364/AO.57.006750](https://doi.org/10.1364/AO.57.006750).
- Kovtoniuk, Bohmann & Semenov, “Nonclassical photocounting statistics with a single on-off detector,” *Phys. Rev. A* 113, 063726 (2026), [DOI 10.1103/pknl-24xd](https://doi.org/10.1103/pknl-24xd).

**Already known:** attenuation sweeps turn binary click/no-click data into a generating-function or photon-statistics measurement; multiple on/off elements yield click-count statistics.

**Surviving distinction:** these works infer the photon-number statistics of a field. The candidate here uses a known random spatial mixing model so that the field-statistics parameters become **mask-indexed scene power sums**, then solves a spatial inverse problem jointly with ordinary buckets.

**Forbidden claim:** “first use of a power/efficiency sweep to recover moments or a generating function.”

## 6.2 SiPM saturation and nonuniform illumination — direct hardware collision

The finite-cell occupancy law and its deviations are established:

- van Dam et al., “A Comprehensive Model of the Response of Silicon Photomultipliers,” *IEEE TNS* 57, 2254 (2010), [DOI 10.1109/TNS.2010.2053048](https://doi.org/10.1109/TNS.2010.2053048).
- Gruber et al., “Over saturation behavior of SiPMs at high photon exposure,” *NIM A* 737, 11 (2014), [DOI 10.1016/j.nima.2013.11.013](https://doi.org/10.1016/j.nima.2013.11.013).
- Weitzel et al., “Measurement of the response of Silicon Photomultipliers from single photon detection to saturation,” *NIM A* 936, 558 (2019), [DOI 10.1016/j.nima.2018.10.074](https://doi.org/10.1016/j.nima.2018.10.074).
- Kumar, Herzkamp & van Waasen, “Non-linearity Simulation of Digital SiPM Response for In-homogeneous Light,” *IEEE TNS* (2021), [DOI 10.1109/TNS.2021.3049675](https://doi.org/10.1109/TNS.2021.3049675).
- Moya-Zamanillo & Rosado, “Understanding the Nonlinear Response of SiPMs,” *Sensors* 24, 2648 (2024), [DOI 10.3390/s24082648](https://doi.org/10.3390/s24082648).

**Already known:** the textbook homogeneous response \(C(1-e^{-\mu/C})\), saturation correction, spatial nonuniformity as a systematic, recovery/refiring, crosstalk, afterpulsing, pulse-shape dependence and over-saturation.

**Surviving distinction:** deliberately choosing a structured, scene-dependent microcell-load distribution and treating its deviation from homogeneous occupancy as an image-bearing statistic.

**Forbidden claim:** “SiPM nonuniform saturation has not been modeled” or “no calibration is needed.” The honest phrase is **no per-cell addressing/calibration under an exchangeable random-mixing model**.

## 6.3 One-bit and quantized compressed sensing — adjacent inverse-problem family

- Boufounos & Baraniuk, “1-Bit Compressive Sensing,” CISS 2008, [DOI 10.1109/CISS.2008.4558487](https://doi.org/10.1109/CISS.2008.4558487).
- Plan & Vershynin, “One-Bit Compressed Sensing by Linear Programming,” *CPAM* 66, 1275 (2013), [DOI 10.1002/cpa.21442](https://doi.org/10.1002/cpa.21442).

**Difference:** one-bit CS quantizes each macroscopic linear projection \(a_k^\top x\). Here many latent on/off channels occur **inside one mask exposure**, and only their aggregate occupancy is read. The mean link can encode the distribution of microscopic loads, not merely a scalar quantization of \(m_k^\top x\).

**But:** this still belongs broadly to nonlinear/generalized compressed sensing. Do not claim that existing nonlinear sensing theory is irrelevant.

## 6.4 Occupancy sketches and balls-into-bins — mathematical cousin

- Flajolet & Martin, “Probabilistic counting algorithms for data base applications,” *JCSS* 31, 182 (1985), [DOI 10.1016/0022-0000(85)90041-8](https://doi.org/10.1016/0022-0000(85)90041-8).

Poissonized occupancy, empty-bin probabilities and generating functions are classical. A SiPM is physically an occupancy counter. The novelty cannot be the occupancy identity by itself.

**Surviving distinction:** the “item weights” are spatially modulated optical intensities tied to an unknown image, and repeated mask-indexed occupancy curves provide coupled inverse constraints.

## 6.5 HBT, thermal-light GI and higher-order GI — closest optical relative

Representative anchors:

- Gatti et al., “Ghost imaging with thermal light: comparing entanglement and classical correlation,” *PRL* 93, 093602 (2004), [DOI 10.1103/PhysRevLett.93.093602](https://doi.org/10.1103/PhysRevLett.93.093602).
- Shapiro, “Computational ghost imaging,” *PRA* 78, 061802 (2008), [DOI 10.1103/PhysRevA.78.061802](https://doi.org/10.1103/PhysRevA.78.061802).
- Chan, O'Sullivan & Boyd, “High-order thermal ghost imaging,” *Opt. Lett.* 34, 3343 (2009), [DOI 10.1364/OL.34.003343](https://doi.org/10.1364/OL.34.003343).
- Chen et al., “High-visibility, high-order lensless ghost imaging with thermal light,” *Opt. Lett.* 35, 1166 (2010), [DOI 10.1364/OL.35.001166](https://doi.org/10.1364/OL.35.001166).

**Classical \(g^{(2)}\)/higher-order GI:** estimates intensity correlation functions from products/coincidences across arms, pixels, or times. It can reveal second and higher optical moments and has long been used for ghost images.

**Candidate distinction here:** one unaddressed aggregate detector, no reference-arm coincidence record, and a controlled attenuation sweep estimates a Laplace/cumulant curve. The scene statistic is intended to be \(M(x^{\circ2})\), not the conventional image formed by correlating a bucket sequence with known reference patterns.

**Forbidden claims:** “first optical second-moment imaging,” “first higher-order ghost imaging,” or “first nonlinear ghost imaging.”

## 6.6 Nonlinear single-pixel/ghost imaging — name collision, different mechanism

Nonlinear GI already names schemes with nonlinear pattern generation, frequency conversion, or field-sensitive THz detection, e.g. Olivieri et al., “Time-Resolved Nonlinear Ghost Imaging,” *ACS Photonics* (2018), [DOI 10.1021/acsphotonics.8b00653](https://doi.org/10.1021/acsphotonics.8b00653).

Therefore **“nonlinear ghost imaging” is too broad and already occupied**. The manuscript must identify the specific microcell-occupancy/Laplace mechanism.

## 6.7 Goodman speckle boundary — mandatory optics citation

Goodman explicitly distinguishes incoherent intensity addition from coherent amplitude addition: J. W. Goodman, “Some fundamental properties of speckle,” *JOSA* 66, 1145 (1976), [DOI 10.1364/JOSA.66.001145](https://doi.org/10.1364/JOSA.66.001145).

This is not background decoration. It is the citation that determines whether the product theorem is physically applicable.

---

# 7. Names

My ranking:

1. **MOLT** — **Microcell-Occupancy Laplace Tomography**  
   Best technical name. It names the detector statistic and the mathematical object without claiming all nonlinear GI.

2. **SCOUT-GI** — **Saturation-Coded Occupancy Tomography for Ghost Imaging**  
   Best campaign-facing name, though “saturation-coded” should be defined as controlled occupancy, not arbitrary deep saturation.

3. **MOST-GI** — **Microcell-Occupancy Sweep Tomography for Ghost Imaging**  
   Less elegant mathematically, but immediately operational.

Recommended theorem name:

> **Masked Laplace–Power-Sum Theorem**

Do not name the central result the “elementary-symmetric-polynomial theorem”; the reciprocal-polynomial identity is exact, but power-sum/log-Laplace coordinates are what remain statistically and physically usable.

---

# 8. Frozen life-or-death gates for the running probe

## Gate P0 — optical mechanism

Use at least three equal-\(p_1\), separated-\(p_2\) scenes/masks. Compare:

1. independent intensity-additive product model;
2. coherent-amplitude collapse model;
3. Gamma/correlated intermediate model.

**Kill** if the product model cannot predict held-out power levels or if fitted curvature tracks only \(p_1\).

## Gate P1 — detector validity

At the proposed pulse width and integration window:

- no evidence of cell refiring/over-saturation;
- crosstalk/afterpulse correction frozen;
- \(C_{\rm eff}\) and \(k_{\rm eff}\) stable;
- no-fire fraction at the high level remains above the frozen guard;
- response reproducible under repeat sweeps and power order reversal.

**Kill** the “self-calibrating” narrative if nuisance drift creates curvature comparable to the scene-induced \(p_2\) effect.

## Gate P2 — Fisher effect size

For the frozen \(C=3600\) detector and equal incident-photon budget, estimate the nuisance-profiled Cramér–Rao bound for each mask.

Authorize \(p_2\) only if:

\[
\operatorname{median}\frac{|\Delta p_2|}
{\operatorname{sd}(\widehat p_2)}
\ge3
\]
over the frozen null-pair bank, with a lower-tail gate such as at least 80% of pairs above 2.

Do not use the optimistic known-\(p_1\), known-\(k\) FI as the final gate.

## Gate P3 — inverse-space teeth

Report:

- rank and singular spectrum of \(M\);
- rank and Fisher-weighted singular spectrum of
  \[
  \begin{bmatrix}
  M\\2M\operatorname{diag}(x)
  \end{bmatrix};
  \]
- the same restricted to the chosen prior tangent/active support;
- adversarial joint-fiber search, not only random null pairs;
- reconstruction against linear-GI under equal masks, equal incident photons, and equal total acquisition time.

A single favorable \(d'\ge3\) pair is a teaser, not a theorem.

## Gate P4 — order boundary

Fit \(p_1,p_2\) as primary. Add \(p_3\) only as a preregistered nested-model diagnostic.

**Kill \(p_3+\)** unless held-out likelihood and CRB show material, reproducible information after multiplicity and nuisance profiling.

---

# 9. Final strategic judgment

This is **not a renaming**, but neither is the current teaser yet the promised universal sledgehammer.

The strongest honest paper shape is:

> **A detector defect becomes a second sensing operator:** under verified intensity-additive speckle mixing, the SiPM occupancy sweep measures a masked Laplace transform whose log-curvature yields \(M(x^{\circ2})\). This adds one generic constraint per mask, doubles local measurement rank, and can close low-dimensional image-prior fibers without adding a detector array or reference arm.

The weakest—and currently forbidden—paper shape is:

> “Deep SiPM saturation automatically reveals all image moments, self-calibrates with 3600 cells, and exactly breaks the 5% GI null space.”

The mainline should continue only through the four gates above. If P0 survives and the nuisance-profiled T2/T3 result remains large, this is a genuinely unusual, theory-led imaging mechanism. If P0 collapses to coherent speckle, or if \(k_{\rm eff}\) makes the \(p_2\) curvature require thousands of gates per mask, kill it immediately and preserve the theorem as a conditional materials-bank result.