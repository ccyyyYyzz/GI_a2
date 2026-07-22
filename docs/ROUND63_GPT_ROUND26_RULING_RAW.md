# R26 ruling (GitHub issue #18, raw)

Title: R26 ruling: information limits and temporal design for dynamic-scattering GI
Posted: 2026-07-22T02:23:16Z

---

# R26 — New theory mainline: information limits and temporal design for dynamic-scattering ghost imaging

**Reference commit:** `ba0089c`  
**Source question:** `docs/ROUND63_GPT_ROUND26_QUESTION.md`

**Scope:** pure theory architecture. The two ROUND63 manuscripts and the RLMI bridge remain independent and untouched.

> **R26 VERDICT.** The dynamic-scattering program is worth pursuing, but its clean spine is narrower than the proposed five-object list.
>
> 1. **O1 is the first theorem.** The exact object is the conditional covariance of the complete object-score given the bucket record. In the Poisson multiplicative-gain model it simplifies to the posterior covariance of a hidden **gain-weighted illumination vector**, the precise analogue of hidden live time.
> 2. **O3 is not a separate foundational object.** It is O1 inserted into a van Trees/posterior Cramér–Rao inequality. This yields a rigorous estimator-independent Bayes-MSE ceiling, including neural estimators, but only under the declared stochastic-gain model and prior.
> 3. **O4 is the second theorem and the main novelty opportunity.** In a linearized correlated-gain model, the information available for the object is a Schur complement. Chopping, interleaving, and temporal randomization are optimal exactly when they make the object score orthogonal to the gain-drift score while preserving oracle object information.
> 4. **O2, as a universal two-axis “temporal ridge,” is not presently a theorem and is false without additional resource assumptions.** A decorrelation-time ratio alone does not create an interior optimum. Read noise, switching overhead, total duration, pattern rank, and the exact within-frame gain functional are load-bearing. O2 should be derived only after those are fixed.
> 5. **O5 composes at the conditional-score-operator level, not by multiplying scalar loss factors.** Gain drift and dead-time saturation interact; high-gain states are preferentially clipped, so the two hidden mechanisms do not separate.
>
> The recommended program is a **two-paper companion pair**: retain the existing identifiability manuscript as Part I; write one new Part II on information limits and temporal design. Do not merge them.

---

## Executive ranking

Scores are on a five-point scale and intentionally conservative.

| object | theorem probability | residual novelty | mathematical depth | Peng–Chen relevance | ruling |
|---|---:|---:|---:|---:|---|
| **O1 — scattering missing information** | **5** | **3** | **5** | **5** | **Prove first.** Generic identity is classical; the exact GI specialization and hidden-exposure corollaries are the real contribution. |
| **O4 — optimal temporal design** | **4** | **4** | **5** | **5** | **Second theorem.** Strongest new-design opportunity, but only under a stated gain basis/process and local model. |
| **O3 — correction ceiling** | **5** | **2** | **4** | **5** | Immediate O1 + van Trees corollary. Important yardstick, not an independent foundational invention. |
| **O2 — temporal ridge** | **2** | **2** | **4** | **5** | Demote to a model-dependent operating map. The universal two-dimensional conjecture is under-specified. |

Weighted order for the next work:

\[
\boxed{\mathrm{O1}\;>\;\mathrm{O4}\;>\;\mathrm{O3}\;>\;\mathrm{O2}.}
\]

O3 should be presented immediately after O1, despite ranking third as a distinct idea, because it costs almost no additional machinery and directly addresses supervised corrections.

---

# Q1 — Aggressive prior-art sweep

## 1. The Peng–Chen dynamic-scattering line already owns the correction problem

The immediate optical neighbors are not hypothetical. The group has developed:

- temporal correction for GI through complex scattering ([Xiao, Zhou & Chen, *Optics Letters* 2022](https://doi.org/10.1364/OL.463897));
- learning-based gain correction with Gaussian constraints ([Peng & Chen, *Optics Letters* 2023](https://doi.org/10.1364/OL.499787));
- dual-detector self-correction ([Zhou, Xiao & Chen, *Optics Express* 2023](https://doi.org/10.1364/OE.489808));
- the supervised dynamic-scaling-factor model that motivates this round ([Peng & Chen, *Applied Physics Letters* 2024](https://doi.org/10.1063/5.0213138));
- photon-limited correction with a single-photon detector ([Hao, Peng, Zhang & Chen, *APL Photonics* 2025](https://doi.org/10.1063/5.0293168)).

That literature establishes methods and experiments. It does **not** appear to establish an exact observed-versus-complete Fisher decomposition, an estimator-independent Bayes risk ceiling, or an information-optimal temporal schedule under correlated gain. Those are the defensible residual questions.

## 2. O1’s abstract identity is classical missing-data theory

The generic statement

\[
I_{\rm obs}=I_{\rm complete}-I_{\rm missing}
\]

is the Orchard–Woodbury/Louis missing-information principle, not a new information law. Louis derived the observed information from the complete-data score and its conditional covariance in latent-variable models ([Louis 1982](https://doi.org/10.1111/j.2517-6161.1982.tb01203.x)).

Therefore the novelty cannot be “we discover missing information for hidden gains.” It must be:

> the exact matrix specialization for programmable bucket imaging, its simplification to a hidden gain-weighted illumination vector in the Poisson channel, and the temporal-design consequences.

## 3. O3’s universal estimator bound is Bayesian estimation theory

Van Trees/posterior Cramér–Rao bounds already lower-bound MSE for arbitrary estimators under a prior. For time-correlated latent states, a canonical reference is Tichavský, Muravchik and Nehorai’s posterior CRB for nonlinear dynamical systems ([IEEE TSP 1998](https://doi.org/10.1109/78.668800)).

The new content is not the inequality. It is the information matrix inserted into it and the resulting dynamic-GI ceiling. A U-Net is covered only because it is one estimator using the same observation and prior—not because the inequality is network-specific.

## 4. Blind gain/object identifiability is crowded

The existing Asset-A manuscript must be positioned against:

- general bilinear identifiability and ambiguity groups ([Li, Lee & Bresler](https://arxiv.org/abs/1501.06120));
- blind self-calibration by linear least squares ([Ling & Strohmer, SIIMS 2018](https://doi.org/10.1137/16M1103634));
- nonconvex blind gain calibration with near-linear sample complexity ([Cambareri & Jacques](https://doi.org/10.1093/imaiai/iay004));
- optimal injectivity conditions for bilinear inverse problems ([Kech & Krahmer](https://doi.org/10.1137/16M1067469)).

Asset A may still contain a distinctive square-design obstruction, temporal-stationarity anchor, and schedule flip, but it should not advertise blind calibration itself as new.

## 5. Self-calibration and temporal coding already appear across imaging

- Astronomy jointly solves sky brightness and detector gains from dithered measurements; Fixsen, Moseley and Arendt explicitly encode calibration into the observing strategy ([Calibrating Array Detectors](https://arxiv.org/abs/astro-ph/0002260)), and Arendt, Fixsen and Moseley compare dithering strategies by a calibration figure of merit ([Dithering Strategies](https://arxiv.org/abs/astro-ph/0002258)).
- Parallel MRI jointly estimates image content and coil sensitivity maps by regularized nonlinear inversion ([Uecker et al. 2008](https://doi.org/10.1002/mrm.21691)).
- Microscopy flat-field methods such as [CIDRE](https://doi.org/10.1038/nmeth.3323) and [BaSiC](https://doi.org/10.1038/ncomms14836) estimate multiplicative illumination and temporal baseline components from image collections.

These precedents substantially narrow O4’s novelty. “Temporal randomization helps calibration” is known. The new claim must be the exact nuisance-information criterion and its specialization to programmable single-pixel measurements through a correlated scattering gain.

## 6. Speckle science already has decorrelation-time/exposure laws

Dynamic light scattering and diffusing-wave spectroscopy relate temporal field correlations to scatterer motion ([Pine et al. 1988](https://doi.org/10.1103/PhysRevLett.60.1134)). Speckle-visibility spectroscopy analytically connects finite exposure to field autocorrelation ([Bandyopadhyay et al. 2005](https://doi.org/10.1063/1.2037987)). Laser-speckle contrast imaging has explicit exposure-time optima: Yuan et al. found a task-specific optimum near 5 ms ([Applied Optics 2005](https://doi.org/10.1364/AO.44.001823)), while later analysis found best precision typically around \(T/\tau_c\in[2,10]\) under its own noise and speckle-statistics model ([Scientific Reports 2023](https://doi.org/10.1038/s41598-023-45303-z)).

Thus “an optimum near the decorrelation time” is not new. O2 must be specific to **object reconstruction from sequential programmable bucket measurements with a hidden scalar gain and a rank/identifiability cost**.

## 7. Fading-channel theory already owns coherence/training tradeoffs

Block-fading and noncoherent communications ask how much resource to spend learning a time-varying multiplicative channel. Hassibi and Hochwald show the tradeoff directly: too little training leaves the channel unknown; too much consumes the coherence block ([IEEE TIT 2003](https://authors.library.caltech.edu/records/y6rsb-e9851)). Noncoherent capacity theories likewise organize performance by coherence dimension and SNR.

This is a close conceptual prior for O2, not a fatal collision. Imaging estimates a fixed object through programmable linear functionals rather than communicates symbols, but the training/coherence architecture is already known.

## 8. O4 also sits inside classical optimal input design

Fisher-optimal input and schedule design for dynamic systems is classical; see Mehra’s optimal-input survey, Goodwin–Payne, and modern reviews such as Pronzato’s optimal experimental design/control framework. Dynamic-system design commonly uses Schur complements, \(D_s\)-optimality, and input spectra; for example, partial-parameter dynamic-system design and equivalence results are treated in [Qureshi & Ng](https://doi.org/10.1137/0320052).

The residual novelty is the **gain-drift nuisance geometry induced by programmable GI patterns**, not KKT certificates or Fisher-optimal scheduling in the abstract.

---

# Q2 — The first theorem: exact hidden-gain missing information

## 9. General theorem

Let \(X\) denote the object parameter, \(L\) the latent gain path, and \(Y\) the observed bucket record. Assume:

1. \(p(l)\) does not depend on \(x\);
2. \(p(y,l\mid x)=p(l)p(y\mid x,l)\) is dominated and differentiable in \(x\);
3. differentiation may pass under the integral;
4. the complete score is square-integrable.

Define the complete object score

\[
U_x(Y,L)=\nabla_x\log p(Y,L\mid x)
        =\nabla_x\log p(Y\mid x,L).
\]

### Theorem O1 — scattering missing-information identity

The observed score is

\[
\boxed{
\nabla_x\log p(Y\mid x)
=
\mathbb E\{U_x(Y,L)\mid Y,x\}.
}
\tag{O1.1}
\]

The complete and observed Fisher matrices satisfy

\[
\boxed{
I_Y(x)
=
I_{Y,L}(x)
-
\mathbb E_Y\!\left[
\operatorname{Cov}\{U_x(Y,L)\mid Y,x\}
\right].
}
\tag{O1.2}
\]

The missing term is positive semidefinite. It vanishes if and only if the complete object score is measurable from the observed record.

### Proof

Differentiate the marginal likelihood:

\[
\nabla_xp(y\mid x)
=
\int \nabla_xp(y,l\mid x)\,dl
=
\int p(y,l\mid x)U_x(y,l)\,dl.
\]

Division by \(p(y\mid x)\) gives (O1.1). Scores have mean zero, so the law of total covariance gives

\[
\operatorname{Cov}(U_x)
=
\operatorname{Cov}\{\mathbb E(U_x\mid Y)\}
+
\mathbb E\{\operatorname{Cov}(U_x\mid Y)\}.
\]

The first covariance is the observed Fisher matrix and the left side is the complete Fisher matrix, proving (O1.2).

This proof is classical Louis theory. The following programmable-GI specialization is the important part.

---

## 10. Exact Poisson bucket corollary

Let the rows of \(M\in\mathbb R_+^{N\times K}\) be \(m_n^\top\), let

\[
s=Mx,\qquad s_n>0,\qquad a_n=e^{l_n}>0,
\]

and conditionally on the entire correlated gain vector \(a\), let

\[
Y_n\mid a,x\sim\operatorname{Poisson}(a_ns_n)
\]

independently across frames. The prior on \(a\) may have arbitrary temporal correlation and need not be Gaussian.

The complete score is

\[
U_x
=
M^\top D_s^{-1}Y-M^\top a.
\tag{O1.3}
\]

Conditionally on \(a\),

\[
I_{Y\mid a}(x)
=
M^\top\operatorname{diag}\!\left(\frac{a_n}{s_n}\right)M.
\]

Therefore:

### Corollary O1-P — hidden exposure-vector identity

\[
\boxed{
I_Y(x)
=
M^\top\operatorname{diag}\!\left(\frac{\mathbb E a_n}{s_n}\right)M
-
M^\top
\mathbb E_Y\{\operatorname{Cov}(a\mid Y,x)\}
M.
}
\tag{O1.4}
\]

### Proof

Average the conditional Fisher matrix for the complete-information term. Given \(Y\), the first term of (O1.3) is fixed, hence

\[
\operatorname{Cov}(U_x\mid Y,x)
=
M^\top\operatorname{Cov}(a\mid Y,x)M.
\]

Insert this into (O1.2).

### What plays the role of live time?

For the dead-time channel, hidden active time is the latent exposure to the rate parameter. Here the corresponding object is

\[
\boxed{
G_M(a)=M^\top a=\sum_{n=1}^Na_nm_n,
}
\tag{O1.5}
\]

—the **gain-weighted illumination/exposure vector**. The information loss is exactly the posterior covariance of this hidden exposure vector:

\[
I_{\rm missing}
=
\mathbb E_Y\operatorname{Cov}\{G_M(a)\mid Y,x\}.
\]

This is the cleanest bridge from ROUND63’s live-time identity to dynamic-scattering GI.

Three consequences are immediate:

1. temporal correlation enters through the full posterior covariance of \(a\), not through one scalar variance;
2. pattern ordering matters because it rotates gain uncertainty through \(M^\top(\cdot)M\);
3. a second monitor detector or other gain side information can only reduce the missing term by conditioning more strongly.

---

## 11. Gaussian readout corollary

For analog bucket data

\[
Y\mid a,x\sim\mathcal N(D_aMx,R),
\qquad R\succ0,
\]

the complete score and information are

\[
U_x=M^\top D_aR^{-1}(Y-D_aMx),
\]

\[
I_{Y,a}(x)=
\mathbb E_a[M^\top D_aR^{-1}D_aM].
\]

Identity (O1.2) remains exact, but the missing term is

\[
\mathbb E_Y\operatorname{Cov}\!\left[
M^\top D_aR^{-1}(Y-D_aMx)\mid Y,x
\right],
\]

not merely \(M^\top\operatorname{Cov}(a\mid Y)M\). The unusually clean hidden-exposure formula (O1.4) is therefore a property of the conditionally Poisson multiplicative channel, not a universal gain law.

---

## 12. First counterexample risk

The first serious risk is not algebraic. It is physical:

> the medium may not produce one scalar gain per frame.

If the dynamic medium changes the spatial transfer operator, the true model is closer to

\[
y_n=\sum_j a_{nj}m_{nj}x_j+\epsilon_n
\quad\text{or}\quad
Y_n=\mathcal K_n(m_n\odot x)+\epsilon_n,
\]

and (O1.4) no longer applies. The general identity (O1.2) survives, but its attractive hidden-vector simplification and the low-dimensional gain prior disappear.

A second risk is gauge singularity: in a nonidentifiable square design, both complete and observed information can be singular. O1 remains true but cannot manufacture identifiability. Asset A’s thresholds must be assumed before interpreting inverse information or Bayes risk.

---

# O3 — The correction ceiling as an O1 corollary

## 13. Estimator-independent Bayesian bound

Let the object have differentiable prior \(\pi(x)\), with prior information

\[
I_\pi=
\mathbb E_\pi[
\nabla\log\pi(X)\nabla\log\pi(X)^\top].
\]

Under the usual van Trees boundary conditions, every estimator \(\widehat x(Y)\) satisfies

\[
\boxed{
\mathbb E[(\widehat x-X)(\widehat x-X)^\top]
\succeq
\left[I_\pi+\mathbb E_XI_Y(X)\right]^{-1}.
}
\tag{O3.1}
\]

Combining with O1 gives

\[
\boxed{
\mathbb E[(\widehat x-X)(\widehat x-X)^\top]
\succeq
\left[
I_\pi+
\mathbb E I_{Y,L}(X)
-
\mathbb E\operatorname{Cov}(U_X\mid Y,X)
\right]^{-1}.
}
\tag{O3.2}
\]

For a task matrix \(W\succeq0\),

\[
\mathbb E\|W^{1/2}(\widehat x-X)\|_2^2
\ge
\operatorname{tr}\!\left\{
W[I_\pi+\mathbb EI_Y(X)]^{-1}
\right\}.
\tag{O3.3}
\]

This covers a U-Net, an untrained network, a joint MAP estimator, or any other estimator that receives the same bucket record and side information.

### What may be claimed

> Under the declared gain prior and observation model, no estimator using the same measurements can achieve Bayes MSE below (O3.3).

### What may not be claimed

- The bound does not show that a particular network is close to optimal.
- It does not directly bound PSNR, SSIM, visibility, or CNR without an explicit transformation/loss argument.
- If a network receives additional calibration data, reference-detector data, or a mismatched training distribution, those must be included in the experiment and prior.
- A lower bound alone does not prove that “more training is pointless”; that conclusion requires the achieved risk to be close to the bound.

---

## 14. How Asset A’s static-correction ceiling can emerge

Use the small-log-gain linearization

\[
a\simeq\mathbf 1+\ell,
\qquad
Y=s+D_s\ell+\epsilon,
\qquad
s=Mx,
\]

with

\[
\ell\sim\mathcal N(0,C_\ell),
\qquad
\epsilon\sim\mathcal N(0,R).
\]

Suppose a static correction removes a declared drift subspace with projector \(P\). The residual gain is

\[
\ell_\perp=(I-P)\ell,
\qquad
C_\perp=(I-P)C_\ell(I-P)^\top.
\]

The marginalized model is

\[
Y\mid x\sim\mathcal N(s,\Sigma_x),
\qquad
\Sigma_x=R+D_sC_\perp D_s.
\tag{O3.4}
\]

Its exact Gaussian Fisher matrix is

\[
[I_Y]_{jk}
=(\partial_js)^\top\Sigma_x^{-1}(\partial_ks)
+rac12\operatorname{tr}\!\left(
\Sigma_x^{-1}\partial_j\Sigma_x
\Sigma_x^{-1}\partial_k\Sigma_x
\right).
\tag{O3.5}
\]

If \(C_\perp=v\bar C\) and \(v\to0\), the covariance-derivative term is \(O(v^2)\), while

\[
\boxed{
I_Y
=
M^\top R^{-1}M
-vM^\top R^{-1}D_s\bar C D_sR^{-1}M
+O(v^2).
}
\tag{O3.6}
\]

A scalarized risk or relative-MSE expansion of (O3.6) yields a term of the form \(vB_L\). Thus Asset A’s static-correction ceiling can be interpreted as the first-order scalarization of the exact missing-information/Bayes-risk architecture—but only under the small-gain Gaussian model, the declared correction projector, and the chosen task scalarization.

That is a theorem path, not an automatic identity between manuscripts.

---

# O2 — Temporal ridge: reject the universal form, retain a conditional program

## 15. Why \(T_f/t_c\) alone does not produce a ridge

Let a stationary continuous gain process have mean one and covariance \(C_a(u)\). If a frame records the frame-average gain

\[
\bar a_T=\frac1T\int_0^Ta(t)\,dt,
\]

then

\[
\boxed{
\operatorname{Var}(\bar a_T)
=
\frac{2}{T^2}\int_0^T(T-u)C_a(u)\,du.
}
\tag{O2.1}
\]

For an OU/exponential covariance

\[
C_a(u)=\sigma_a^2e^{-|u|/t_c},
\qquad r=T/t_c,
\]

\[
\boxed{
\operatorname{Var}(\bar a_T)
=
2\sigma_a^2\frac{r-1+e^{-r}}{r^2}.
}
\tag{O2.2}
\]

This decreases from \(\sigma_a^2\) at short exposure to approximately \(2\sigma_a^2/r\) at long exposure. Gain averaging alone therefore gives no interior optimum.

If instead the model takes one quasi-static gain sample per frame, changing \(T_f\) primarily changes the correlation between successive gains and the number of patterns collected—not the marginal gain variance. Again, no universal ridge follows.

An interior optimum requires a competing mechanism, such as:

- finite total acquisition time, so \(N\approx T_{\rm total}/(T_f+t_{\rm switch})\) decreases with frame duration;
- read noise or switching overhead at short frames;
- insufficient photons per frame;
- loss of object/gain identifiability when too few distinct patterns remain;
- within-frame model failure when the medium changes spatially rather than through one scalar average.

## 16. Correct dimensionless architecture

A minimally complete operating map generally needs at least

\[
r=\frac{T_f}{t_c},
\qquad
\eta=\lambda t_c
\quad\text{(photons per correlation time)},
\]

plus

\[
R=\frac{T_{\rm total}}{t_c},
\qquad
\omega=\frac{t_{\rm switch}}{t_c},
\]

and an object/design dimension such as \(N/K\). A read-noise-to-shot-noise ratio may be another independent parameter.

Therefore the requested two-axis map \((T_f/t_c,\text{photon budget})\) is not universal. The correct target is a **resource-conditioned temporal operating surface**. A ridge may emerge on a fixed slice of \((R,\omega,N/K,\text{noise model})\), but its constant and even its existence are model-dependent.

### Rigor label

- Exact: O1 and the frame-average covariance (O2.1).
- Long-window/model-specific: information rates after a precise gain process and resource constraint are fixed.
- **CONJECTURE:** a matched temporal ridge for a specified OU-gain, fixed-total-time, read-noise-plus-shot-noise GI model.
- Named obstacle: coupled growth/decline of pattern rank, gain posterior uncertainty, and per-frame photon information. There is no scalar renewal reduction analogous to dead time.

---

# O4 — Optimal temporal design as nuisance-score orthogonalization

## 17. Linearized drift model

Let the gain path be represented locally by

\[
\ell=H\beta+r,
\]

where the columns of \(H\in\mathbb R^{N\times p}\) are declared low-frequency or OU/KL drift modes, \(\beta\) has prior precision \(\Lambda_\beta\), and \(r\) is residual gain noise. At \(\beta=0\), use the local Gaussian model

\[
y=M_\pi x+D_{s_\pi}H\beta+\epsilon,
\qquad
s_\pi=M_\pi x,
\qquad
\epsilon\sim\mathcal N(0,R),
\tag{O4.1}
\]

where \(\pi\) denotes the temporal order/schedule of the programmable patterns.

The joint information matrix for \((x,\beta)\) is

\[
\mathcal I(\pi)=
\begin{bmatrix}
I_{xx} & I_{x\beta}\\
I_{\beta x} & I_{\beta\beta}
\end{bmatrix},
\]

with

\[
I_{xx}=M_\pi^\top R^{-1}M_\pi,
\]

\[
I_{x\beta}=M_\pi^\top R^{-1}D_{s_\pi}H,
\]

\[
I_{\beta\beta}=H^\top D_{s_\pi}R^{-1}D_{s_\pi}H+\Lambda_\beta.
\]

The efficient/Bayesian information for the object after treating drift as nuisance is the Schur complement

\[
\boxed{
I_{x\mid\mathrm{gain}}(\pi)
=
I_{xx}(
\pi)
-I_{x\beta}(\pi)I_{\beta\beta}(\pi)^{-1}I_{\beta x}(\pi).
}
\tag{O4.2}
\]

## 18. Theorem O4-A — nuisance-orthogonal schedule

Assume \(I_{\beta\beta}\succ0\). For schedules having the same oracle object information \(I_{xx}\),

\[
I_{x\mid\mathrm{gain}}\preceq I_{xx},
\]

with equality if and only if

\[
\boxed{
M_\pi^\top R^{-1}D_{M_\pi x}H=0.
}
\tag{O4.3}
\]

### Proof

The subtracted term in (O4.2) is of the form

\[
AA^\top,
\qquad
A=I_{x\beta}I_{\beta\beta}^{-1/2},
\]

and is therefore positive semidefinite. It vanishes exactly when \(I_{x\beta}=0\), proving the result.

### Meaning

The optimal temporal schedule does not merely make gains “stationary.” It makes the object score orthogonal, in the noise metric, to every admitted drift mode. Expanded by columns \(h_j\) of \(H\), the exact moment conditions are

\[
\boxed{
\sum_{n=1}^N h_j(t_n)\,
\frac{s_n}{\sigma_n^2}\,m_n=0,
\qquad j=1,\ldots,p,
}
\tag{O4.4}
\]

for diagonal \(R=\operatorname{diag}(\sigma_n^2)\).

This is the theorem behind chopping, interleaving, and randomization:

- **chopping/paired designs** are exact when their signed effective rows annihilate the low-order temporal moments in (O4.4);
- **interleaving** reduces the correlation between object-carrier directions and low-frequency gain modes;
- **random permutation** makes the cross block mean-zero under centered exchangeable rows and gives concentration around zero.

The first two are exact only when the pattern/carrier moment equations actually hold. “Alternating is always optimal” is false.

## 19. Randomization proposition and obstacle

Under bounded centered carrier vectors and a uniformly random permutation, standard matrix concentration can bound

\[
\|I_{x\beta}\|_{\rm op}
=
O_{\mathbb P}
\left(
B\sqrt{\frac{\log((K+p)/\delta)}{N}}
\right),
\]

which in turn bounds the Schur-complement loss by

\[
\|I_{xx}-I_{x\mid\mathrm{gain}}\|_{\rm op}
\le
\frac{\|I_{x\beta}\|_{\rm op}^2}
{\lambda_{\min}(I_{\beta\beta})}.
\]

**Status:** theorem under explicitly bounded, centered, permutation-sampling assumptions. It is not yet a theorem for arbitrary nonnegative GI patterns because \(s_n=m_n^\top x\) is scene-dependent and can destroy centering.

The named obstacle is the product \(s_nm_n\): temporal randomization randomizes a scene-dependent carrier, not a fixed zero-mean row sequence. Asset A’s stationary-carrier hypothesis is exactly the condition needed to close this proof.

## 20. Relation to the existing flip boundary

The paired schedule can spend measurements/references to reduce \(I_{x\beta}\); a randomized schedule preserves more distinct object rows but only suppresses the cross block statistically. The correct flip criterion is therefore equality of a task criterion applied to the two Schur complements,

\[
\Phi\{I_{x\mid\mathrm{gain}}^{\rm pair}\}
=
\Phi\{I_{x\mid\mathrm{gain}}^{\rm rand}\},
\]

under the same total-time and photon accounting.

Asset A’s existing \(\rho^*\) flip may be a closed-form instance of this equality. It should be rederived from (O4.2), not merely re-described as a heuristic schedule boundary.

O4’s measure-valued/KKT extension is then classical optimal design applied to atoms that carry both a spatial pattern and a temporal slot. The novel content is the correlated-gain information atom and the nuisance-orthogonality geometry.

---

# O5 — Composite gain drift plus dead-time counting

Let the complete latent record contain the gain path, photon arrivals, and detector state. Conditional-score operators compose through nested observation channels by the tower property, so the controlled-score Gramian architecture survives.

What does **not** survive is a scalar product law. If a frame’s incident rate is \(a_ns_n\) and the counter response is nonlinear, the score depends on

\[
\frac{\partial}{\partial x}\log p_{\rm count}(Y_n\mid a_ns_n),
\]

so the posterior uncertainty of \(a_n\) is weighted by the derivative and curvature of the dead-time count law. High-gain states are preferentially compressed or saturated. Consequently:

- gain uncertainty and dead-time loss are coupled;
- temporal carriers are clipped most strongly exactly when \(a_n\) is large;
- a correction trained for an analog linear bucket may become biased under photon-counting saturation;
- the best temporal schedule may shift because extreme-gain frames carry less marginal object information.

The exact O1 identity still holds if the complete score includes both hidden mechanisms. The first composite paper should derive one explicit Poisson-renewal specialization; it should not claim a mechanism-by-mechanism multiplicative efficiency law.

---

# Q3 — Honest scope and the physical validity boundary

## 21. Hypothesis MG: when one scalar gain per frame is physically valid

Write the dynamic medium’s intensity-transfer functional at time \(t\) as \(\mathcal K_t\). The scalar multiplicative model is a valid approximation on the experiment’s illuminated object subspace \(\mathcal S\) only if there exist a fixed functional \(\mathcal K_0\), a positive scalar \(a(t)\), and a small residual \(E_t\) such that

\[
\boxed{
\mathcal K_t|_{\mathcal S}
=
a(t)\mathcal K_0|_{\mathcal S}+E_t,
\qquad
\|E_tz\|
\le
\varepsilon a(t)\|\mathcal K_0z\|
\quad\forall z\in\mathcal S.
}
\tag{MG}
\]

After frame integration, \(a_n\) is either the quasi-static gain or the properly weighted frame average.

Plausible sufficient physical conditions are:

1. the bucket collects enough output speckles/modes that temporal medium variation mainly changes total throughput rather than the normalized spatial transfer;
2. the object and medium are static enough within a pattern exposure for a scalar frame average to be meaningful;
3. the pattern does not perturb the medium;
4. incoherent/intensity-linear propagation is an adequate description, or coherent cross terms average out at the bucket;
5. additive scattered background is separately modeled;
6. the same collection aperture, polarization, and wavelength distribution apply across patterns.

## 22. Where the model breaks

The scalar model fails when:

- the dynamic medium changes the spatial transmission matrix across object coordinates;
- only a few speckles are collected, producing pattern-dependent non-Gaussian fluctuations;
- coherent bucket integration retains interference cross terms that depend on the displayed pattern;
- the medium evolves substantially within a frame and cannot be represented by one scalar average;
- object motion and medium motion couple to the pattern order;
- polarization, wavelength, or path-length composition changes with time;
- multiple scattering creates an additive or convolutive component rather than common multiplicative throughput.

In those regimes the correct latent object is a vector/matrix-valued transfer process. O1 survives in its general form, but the low-dimensional correction theory and Poisson hidden-exposure corollary do not.

## 23. Claims this program must not make

- No universal image-quality bound from Fisher information alone.
- No claim that a neural network is close to the ceiling without measuring its gap to the bound.
- No universal two-dimensional temporal ridge.
- No claim that temporal randomization is always optimal.
- No claim that the scalar gain model describes arbitrary dynamic turbid media.
- No claim that identifiability thresholds imply statistical efficiency.
- No claim that an information-optimal schedule is automatically PSNR-optimal under nonlinear regularization.

The clean hierarchy is:

\[
\text{physical scalar-gain hypothesis}
\rightarrow
\text{identifiability}
\rightarrow
\text{observed information}
\rightarrow
\text{Bayes risk bound}
\rightarrow
\text{temporal design}.
\]

Each arrow requires assumptions; none may be skipped.

---

# Q4 — Program architecture

## 24. Do not merge the existing identifiability manuscript with the new theory

The existing manuscript asks:

> when can the object and gain class be separated at all?

The new work asks:

> once separation is possible, how much object information survives hidden gain, what is the estimator-independent risk floor, and how should patterns be ordered in time?

Those are companion questions, not sections of one oversized paper.

### Recommended two-paper dynamic-scattering program

## Part I — existing identifiability paper

Keep:

- square-design obstruction and gauge;
- local/uniform sample thresholds;
- stationarity anchor;
- temporal randomization for identifiability;
- finite-noise relMSE identity and static-correction term;
- paired-versus-randomized flip boundary.

Its conclusion should state that identifiability is necessary but does not determine efficiency.

## Part II — new information/design paper

Recommended spine:

1. **Physical hypothesis MG** and channel model.
2. **Exact gain-path missing-information theorem** O1.
3. **Poisson hidden-exposure corollary** and analog Gaussian counterpart.
4. **Universal correction ceiling** O3 via van Trees, with the small-gain derivation connecting to Asset A’s \(vB_L\) term.
5. **Nuisance-orthogonal temporal-design theorem** O4.
6. **Randomization/chopping corollaries** and rederivation of the flip boundary as a Schur-complement comparison.
7. **One model-specific temporal operating surface**, not a universal ridge claim.
8. **One reach section** on gain plus dead-time counting.

O2 should not appear in the title until a concrete resource-conditioned law has been derived. A safe working title is:

> **Information limits and temporal design for ghost imaging through dynamic scattering media**

## 25. Why this is one new paper, not two

O1 and O3 are one theorem package. O4 becomes meaningful only after O1 defines the information being optimized. Splitting “information ceiling” and “temporal design” would create two incomplete papers and recreate the fragmentation the operator is trying to eliminate.

The correct relationship is therefore:

\[
\boxed{
\text{existing identifiability manuscript}
\quad+\quad
\text{one information-and-design companion}.
}
\]

Do not merge. Do not create a separate “temporal ridge” paper at this stage.

---

# Final blunt assessment

The pivot is scientifically sound. Dynamic-scattering GI is closer to the target group’s experimental axis than another layer of dead-time OED, and the ROUND63 machinery transfers cleanly at the score level.

But the program will fail aesthetically and mathematically if it begins by announcing a new temporal ridge. The first theorem is simpler and stronger:

> **hidden dynamic gain removes exactly the conditional covariance of the complete object score; in the Poisson bucket channel, this is the posterior covariance of the gain-weighted illumination vector.**

That theorem immediately yields the correction ceiling and tells temporal design what it must suppress. The next real theorem is nuisance-score orthogonalization. Only after those are proved should an operating ridge be sought inside a fully specified time/photon/overhead model.

That order converts the five proposed objects into one chain rather than five partially connected claims.