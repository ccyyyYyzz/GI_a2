# R19 ruling (GitHub issue #11, raw)

Title: R19 ruling: post-unblinding analysis correction and proof completion
Posted: 2026-07-19T17:53:11Z

---

# R19 — Final ruling: post-unblinding analysis correction and proof completion

**Audit target:** commit `dbba1a5`, especially:

- `docs/ROUND63_GPT_ROUND19_QUESTION.md`;
- `docs2/01_FROM_AUDITOR_FINDINGS.md` and `docs2/03_FROM_AUDITOR_CLARIFICATION.md`;
- `results/round63_m1/AUDIT_2026-07-19/`;
- the immutable campaign tag `m1-freeze` at `6f00932`.

**Scope:** this ruling adjudicates analysis/specification discrepancies discovered after unblinding. It does not alter raw cells, the immutable tag, the frozen endpoint definitions, the scene cohort, any reconstruction, or any threshold. The correction must be additive and fully provenance-preserving.

> **R19 VERDICT: the elapsed optical-time endpoint is unambiguous and the frozen analyzer was nonconformant. Publish the spec-conformant corrected verdicts through a dated correction artifact. Preserve the original frozen outputs unchanged and show them only in the provenance/correction record.**
>
> The corrected scientific results are:
>
> - `RIDGE_OPERATING_PASS = TRUE`: median `1.86692 dB`, fixed-six-strata paired-seed bootstrap lower bound `1.41348975 dB`, `19/24` positive;
> - `RIDGE_SPEED_PASS = TRUE` on elapsed optical time `T_opt=M_total*nu*tau`: median `19.1270431`, lower bound `18.3284924`, `21/24` above one;
> - the corrected `nu*rho` analysis is retained only as a post-hoc incident-exposure sensitivity: median `0.2786469`, lower bound `0.2397469`, `1/24` above one, with no preregistered verdict.

## Executive ruling

| Question | R19 ruling |
|---|---|
| **Q1 — speed axis** | **Elapsed `T_opt` is the frozen axis.** The specification chain is explicit and not ambiguous. The `nu*rho` implementation is an analysis artifact. Correct the verdict additively; retain `nu*rho` only as a labeled post-hoc photon-exposure sensitivity. |
| **Q2 — bootstrap/PAVA/censoring** | **Use the full spec-conformant machinery.** Preserve all six strata, resample seeds pairwise within image, rebuild curves/target/censoring inside every replicate, use block-collapse equal-weight PAVA, and apply the complete frozen censoring taxonomy. Corrected numbers replace the nonconformant numbers in the abstract and main results; originals appear only in the dated correction/provenance table. |
| **Q3 — certificate semantics** | The campaign has **two empirical verdicts plus one descriptive structural analysis**. Remove every categorical `FULL_STACK_CERT_PASS`/gate/failure/confirmatory-certificate statement. Report `0 CERTIFIED / 299 COUNTEREXAMPLE / 181 NUMERICAL_UNRESOLVED` descriptively. |
| **Q4 — figures and paper-1 boundaries** | Comparator and wording corrections approved. Panel (e) must use `SCAT32-SAFE` for the Q90 comparison and elapsed `T_opt`; remove “preferred” and “+ certificate.” Limit the exact-likelihood uniqueness statement to the shown grid/TV evidence, delete the universal convex-regularization negative, and label `Gamma=1` only as a descriptive onset proxy. |
| **Q5 — certificate addendum** | Approved. Never synthesize unavailable fields. Mark them `N/A`, disclose the frozen `arisk_ratio=1.0` as an implementation constant rather than a measurement, report the true coefficient-range distribution, and regenerate only deterministic primal-screen fields with hashes. |
| **Q6 — proofs** | The LaTeX-ready proof blocks below are approved. The crossover proposition is rigorous **conditional on the stated uniform `C^2` triangular-array expansion**; it does not by itself remove the main text’s “matched-asymptotic prediction” label for the whole hold-law class. The alignment proposition requires a projected normal-cone residual at bounded coordinates; the interior formula in the current text is recovered as a special case. |

---

# 1. Q1 — Frozen speed axis and correction route

## 1.1 The endpoint is not ambiguous

The controlling documents agree:

1. `ROUND63_METHOD_SPEC_M1.md` defines
   \[
   T_{\rm opt}=M_{\rm total}\,\nu\,\tau.
   \]
2. The transferred R7 endpoint fits the isotonic curves and their crossings in `log T_opt`.
3. R17 and R18 explicitly retain the Q90 machinery unchanged.
4. The raw rows persist `optical_time_s=M*nu*tau`, verified to numerical precision.

By contrast, the analyzer used `log(nu*rho)`. Since `nu*rho=lambda T`, this is an incident-signal-exposure coordinate, not elapsed integration time. It changes the resource being measured and cannot supersede the explicit frozen definition merely because it was implemented in the tagged analyzer.

The correction is therefore an **implementation-to-specification correction**, not a post-hoc endpoint redesign.

## 1.2 Required additive correction artifact

Do not edit or replace:

```text
m1-freeze
results/round63_m1/M1_VERDICTS.json
results/round63_m1/M1_VERDICTS.md
any raw shard CSV
```

Create a dated additive bundle, for example:

```text
results/round63_m1/CORRECTION_2026-07-19/
    M1_VERDICTS_SPEC_CORRECTED_R19.json
    M1_VERDICTS_SPEC_CORRECTED_R19.md
    ANALYSIS_CORRECTION_DISCLOSURE.md
    scripts/
    outputs/
    SHA256SUMS
```

The bundle must identify:

- immutable input tag and raw-file hashes;
- the original analyzer outputs;
- D1, D2, and D3 separately;
- the independent auditor and date of discovery;
- the promoted independent implementation;
- the in-project verification result (`18/18` checked numbers, zero deep-diff leaves);
- software/environment and RNG tags;
- corrected source-of-record values.

Use an analysis-version field such as:

```text
analysis_version = R19_SPEC_CORRECTION_V1
```

## 1.3 Frozen corrected-result sentence

Use this sentence in the abstract/results, with normal journal rounding:

> **Using the preregistered elapsed optical-time coordinate, \(T_{\mathrm{opt}}=M_{\mathrm{total}}\nu\tau\), and the frozen Q90/PAVA/censoring procedure, the ridge time-to-quality endpoint passed: the median conservative speed ratio was \(19.13\times\), the family-stratified nested-bootstrap 2.5th-percentile bound was \(18.33\times\), and 21 of 24 scenes had \(S_{\mathrm{gate}}>1\).**

The exact source-of-record values are:

```text
median = 19.127043091646133
LB2.5  = 18.328492357080282
breadth = 21/24
```

## 1.4 Frozen correction-disclosure paragraph

> **Analysis correction.** After the campaign was unblinded, an independent peer audit found that the frozen scoring implementation did not reproduce three parts of the preregistered analysis. It placed the Q90 crossings on \(\log(\nu\rho)\), an incident-photon-exposure coordinate, rather than on the frozen elapsed optical time \(\log T_{\mathrm{opt}}\); it resampled family labels rather than retaining all six fixed strata with paired within-image seed resampling; and it used a nonstandard isotonic implementation while omitting the frozen censoring taxonomy. The immutable tag, raw cells, and original analyzer outputs were preserved unchanged. A dated correction artifact recomputed the endpoints from the same frozen raw cells using the preregistered axis and machinery; the promoted independent implementation and an in-project rerun agreed on all 18 audited numerical outputs with no deep differences. Scientific claims and tables use the spec-conformant corrected results, while the original outputs remain available in the provenance record.

Do not describe the correction as a new endpoint, alternative analysis choice, or threshold change.

## 1.5 Frozen `nu*rho` sensitivity wording

> As a post-hoc resource sensitivity, we repeated the corrected Q90 machinery on the incident-exposure coordinate proportional to \(
u\rho\). It yielded a median ratio of \(0.279\), a 2.5th-percentile bound of \(0.240\), and only 1 of 24 scenes above one. Thus ridge tracking greatly reduced elapsed dwell in this campaign but did not reduce incident exposure to Q90 under this coordinate. This sensitivity was not the preregistered speed endpoint and carries no PASS/FAIL verdict.

Do not call this result “optical time,” “the preregistered speed negative,” or a universal photon-efficiency theorem.

---

# 2. Q2 — Spec-conformant endpoint machinery

## 2.1 Equal-weight PAVA

For each seed-resampled image curve, run ordinary nondecreasing equal-weight PAVA by maintaining one record per current block:

```text
(block mean, block weight, block start/end)
```

When adjacent block means violate monotonicity, replace the two blocks by one block whose mean is their weight-averaged mean and whose weight is the sum. Expand the final block means back to the original dwell positions only after all merging is complete.

The test vector

```text
[3, 2, 1] -> [2, 2, 2]
```

is a required unit test. Duplicating pooled weights at each original position is not PAVA.

## 2.2 Per-image speed endpoint

For each image and seed draw:

1. Apply the **same five resampled seed indices** to `SCAT32-SAFE` and `RIDGE-SCAT32` at every dwell.
2. Average the resampled seed curves separately for safe and ridge.
3. Apply block-collapse PAVA to both curves.
4. Define the safe dynamic range from the fitted endpoints:
   \[
   R_j=\widetilde Q_{\rm safe}(T_{\max})-
       \widetilde Q_{\rm safe}(T_{\min}).
   \]
5. If `R_j < 0.50 dB`, set
   ```text
   SAFE_RANGE_UNINFORMATIVE
   S_gate = 1
   ```
6. Otherwise use
   \[
   Q_{90,j}=\widetilde Q_{\rm safe}(T_{\min})+0.9R_j
   \]
   and interpolate the first crossing linearly in `log T_opt`.
7. Apply the complete frozen taxonomy:
   - `NORMAL`;
   - `FAST_LEFT_CENSORED`;
   - `BOTH_LEFT_CENSORED`;
   - `SAFE_LEFT_FAST_INTERIOR`;
   - `FAST_RIGHT_CENSORED`;
   - `SAFE_RANGE_UNINFORMATIVE`;
   - `ANALYSIS_FAILURE`.

No image is deleted. `m1_microtexture_1` is `SAFE_RANGE_UNINFORMATIVE` because its fitted safe range is `0.27336 dB`; it is clamped to `S_gate=1`. This censoring step, not PAVA, changes breadth from 22 to 21.

## 2.3 Nested family-stratified bootstrap

Freeze `B=10,000` and the existing RNG tags:

```text
operating: tag 13
speed:     tag 14
```

Each bootstrap replicate must:

1. draw five seed indices with replacement **inside every image**;
2. use those paired indices for both compared arms;
3. rebuild the endpoint, including PAVA, target, crossing, and censoring;
4. within each of the six fixed families, sample four of that family’s four images with replacement;
5. concatenate the six four-image strata and take the 24-image median.

The family labels themselves are never resampled. Every replicate contains all six declared strata.

## 2.4 Corrected source-of-record numbers

### Ridge operating-point primary

```text
median deltaQ = 1.866920 dB
LB2.5         = 1.41348975 dB
positive      = 19/24
RIDGE_OPERATING_PASS = TRUE
```

### Ridge elapsed-time speed secondary

```text
median S_gate = 19.127043091646133
LB2.5         = 18.328492357080282
S_gate > 1    = 21/24
RIDGE_SPEED_PASS = TRUE
```

### Post-hoc `nu*rho` sensitivity

```text
median        = 0.27864687917406294
LB2.5         = 0.23974692258639926
above 1       = 1/24
no categorical verdict
```

The corrected numbers **replace** the nonconformant numbers in the abstract, main results, conclusions, figure labels, and cross-arm summary. The original numbers appear alongside the corrected values only in the dated correction/provenance table and remain in their immutable source files.

The “axis-only” lower bound `16.5778528`, which preserves other analyzer defects, is diagnostic only and is not a manuscript result.

---

# 3. Q3 — Certificate semantics after `FALLBACK_DESCRIPTIVE`

## 3.1 Frozen campaign structure

The campaign reports:

1. one corrected empirical operating-point verdict;
2. one corrected empirical elapsed-time speed verdict;
3. one **descriptive structural analysis** of the full-stack certificate machinery.

Delete every use of:

```text
FULL_STACK_CERT_PASS
certificate PASS
certificate FAIL
certificate gate
confirmatory certificate verdict
failed confirmatory certificate
```

The source-of-record status distribution remains:

```text
CERTIFIED             0
COUNTEREXAMPLE       299
NUMERICAL_UNRESOLVED 181
TOTAL                480
```

A `COUNTEREXAMPLE` means a frozen support-expanding primal screen found a full-stack-feasible design whose actual local log-determinant improvement exceeded the declared `0.01*r` threshold. A `NUMERICAL_UNRESOLVED` cell supplies neither certification nor refutation under the frozen numerical protocol. Neither status is a categorical campaign gate.

## 3.2 Frozen replacement wording

Replace “near-optimality is certified, not assumed” by:

> **The balanced design’s local near-optimality is evaluated rather than assumed.**

Immediately after the two frozen R18 mechanism paragraphs, add:

> **Before the immutable campaign tag, the required 24-cell development feasibility screen selected the preregistered `FALLBACK_DESCRIPTIVE` branch. Consequently, `FULL_STACK_CERT_PASS` was removed before freeze, and the post-freeze full-stack computation was specified only as a descriptive three-status analysis, not as a PASS/FAIL endpoint.**

Frozen results wording:

> Across the 480 frozen full-stack analysis cells, none was certified, 299 produced feasible counterexamples, and 181 were numerically unresolved. The counterexamples provide constructive evidence of local geometry headroom relative to deployed SCAT32 in those cells; the unresolved cells provide no finite upper bound. Because the pre-freeze descriptive branch had already been selected, this distribution is reported without a categorical certificate verdict.

Use “preregistered descriptive analysis” or “frozen descriptive analysis,” not “confirmatory certificate.”

---

# 4. Q4 — Figure and claim-boundary corrections

## 4.1 Act III panel (e)

For a Q90 time-to-quality panel:

- comparator: `SCAT32-SAFE`, not `SCAT32-060`;
- fast arm: `RIDGE-SCAT32`;
- horizontal resource: elapsed
  \[
  T_{\rm opt}=M_{\rm total}\nu\tau;
  \]
- caption: “elapsed optical integration time” or “elapsed dwell,” not “photon time”;
- remove “preferred,” “+ certificate,” and any certified-pattern recommendation.

A separate `nu*rho` panel is allowed only if explicitly labeled:

```text
post-hoc incident-exposure sensitivity
```

and visually separated from the preregistered elapsed-time endpoint.

## 4.2 Paper 1 — exact-likelihood scope

Use wording of this form:

> On the reported small exact-reference grid and under the frozen TV settings, the exact and approximate likelihood reconstructions agreed to the stated numerical tolerance, and the exact calculation added little runtime at those tested sizes. This empirical agreement does not establish uniqueness of the TV-regularized solution; convexity of the data term alone is insufficient when the design is underdetermined and the TV penalty is not strictly convex.

Delete any general statement that the exact likelihood “guarantees a unique optimum.”

## 4.3 Paper 1 — Gaussian pathology scope

Delete:

> no amount of convex-in-x regularization repairs the Gaussian pathology

Use:

> The exhibited heteroscedastic-Gaussian pseudo-likelihood is noncoercive along the displayed high-load direction, and the tested TV settings did not remove the resulting failure. Sufficiently coercive penalties can bound such a direction, so no universal claim is made over all convex regularizers.

## 4.4 Mechanism figure — `Gamma=1`

Freeze the label:

> **\(\Gamma=1\): descriptive photon-limited-onset proxy, marking bucket modulation comparable to count noise.**

The caption must add:

> It is not a predicted or validated boundary for high-flux image benefit.

Do not use “separates the region where high flux helps,” “critical benefit boundary,” or equivalent language.

---

# 5. Q5 — Dated certificate-disclosure addendum

## 5.1 Never synthesize missing diagnostics

The original frozen certificate rows remain untouched. The addendum must be status-aware.

### `COUNTEREXAMPLE`

The primal screen short-circuited the dual solve. Persist and report measured primal fields, including the feasible lower-bound gap and deterministically recoverable composition fields. Any dual, cut, complementarity, or coefficient field not computed before the short circuit is:

```text
N/A — dual not executed after primal counterexample
```

### `NUMERICAL_UNRESOLVED`

No successful certificate exists. Failed-attempt metadata and wall time are reportable, but unavailable bound, cut, PSD, complementarity, generalized-eigenvalue, and measured A-risk fields are:

```text
N/A — no numerically valid certificate
```

### `CERTIFIED`

There are no such M1 cells. Do not populate a synthetic row or template example in the campaign table.

## 5.2 Hard-coded `arisk_ratio=1.0`

The frozen value must be disclosed as:

```text
original_arisk_ratio = 1.0
original_arisk_ratio_source = IMPLEMENTATION_CONSTANT_NOT_MEASURED
corrected_arisk_ratio = N/A
```

Do not plot, average, or interpret the hard-coded value as an observed ratio. If a separate deployed-reference ratio equal to one by definition is desired, give it a distinct field name and state the definition explicitly.

## 5.3 Coefficient range

Replace the scalar shorthand “4000” by the actual frozen distribution:

```text
minimum = 2134.505
maximum = 5015.934
median  = approximately 4000
```

The addendum should calculate and archive the exact median and any quoted quantiles from the immutable rows rather than round all cells to one representative value.

## 5.4 Deterministic regeneration

Primal-screen composition diagnostics may be regenerated only when:

- all inputs are immutable frozen files;
- the exact frozen algorithm and ordering are used;
- no new solver option, random seed, threshold, or iteration count is introduced;
- the script, environment, input/output hashes, and deep-diff verification are archived.

Do not regenerate unavailable dual diagnostics by running a new or improved certificate solver. Such a run would be a new post-hoc analysis and must remain separately labeled.

---

# 6. Q6 — LaTeX-ready supplement proof blocks

## 6.0 Installation note

The current supplement loads `amsmath` and `amssymb` but not `amsthm`. Add:

```latex
\usepackage{amsthm}
```

or replace each `proof` environment below by an equivalent bold “Proof.” paragraph.

The crossover proposition below is an **optimizer theorem conditional on a uniform `C^2` renewal-information expansion**. It rigorously derives the crossover constant once that expansion holds. It does not, by itself, prove the triangular-array renewal expansion for every hidden-hold family; therefore the main text must retain the word “prediction” for the full universality-class crossover unless that separate uniform expansion is established.

## 6.1 Uniform crossover-constant proposition

```latex
\begin{proposition}[Uniform crossover constant from the triangular expansion]
\label{prop:s-uniform-crossover}
Let $t=\nu^{-1/3}$, let $x=12\nu c_v^2$, and write
$\rho=t^{-1}y$. Fix $X<\infty$ and a compact interval
$K\Subset(0,\infty)$ containing
\[
  y_0(x)=\left(\frac{6}{1+x}\right)^{1/3}
\]
for every $x\in[0,X]$. Assume that, uniformly for
$x\in[0,X]$ and $y\in K$, the finite-window count-information rate has the
$C^2$ expansion
\begin{equation}
  J_{\nu,c_v}(t^{-1}y)
  =1-tF_1(y,x)+t^2F_2(y,x)+R_t(y,x),
  \label{eq:s-uniform-expansion}
\end{equation}
where
\begin{align}
  F_1(y,x)&=\frac1y+\frac{1+x}{12}y^2,\\
  F_2(y,x)&=\frac1{y^2}+\frac{x-2}{12}y
             +\frac{(1+x)^2}{144}y^4,
\end{align}
and
\[
  \sup_{x\in[0,X],\,y\in K}
  \left(|R_t|+|\partial_yR_t|+|\partial_y^2R_t|\right)=O(t^3).
\]
Then, for all sufficiently large $\nu$, the principal maximizer in $K$ is
unique and satisfies, uniformly for $x\in[0,X]$,
\begin{equation}
  \rho_*(\nu,c_v)
  =\left(2c_v^2+\frac1{6\nu}\right)^{-1/3}
   -\frac{12\nu c_v^2+4}{6(12\nu c_v^2+1)}
   +O(\nu^{-1/3}).
  \label{eq:s-uniform-ridge}
\end{equation}
In particular, the constant correction tends to $-2/3$ on the deterministic
line and to $-1/6$ in the jitter-dominated matching limit.
\end{proposition}

\begin{proof}
Maximizing Eq.~\eqref{eq:s-uniform-expansion} is, to first order, equivalent
to minimizing $F_1$. Its derivatives are
\[
  \partial_yF_1=-\frac1{y^2}+\frac{1+x}{6}y,
  \qquad
  \partial_y^2F_1=\frac2{y^3}+\frac{1+x}{6}.
\]
Thus $F_1$ has the unique minimizer
$y_0=(6/(1+x))^{1/3}$ and
$\partial_y^2F_1(y_0,x)=(1+x)/2>0$. The uniform $C^2$ remainder and the
implicit-function theorem therefore give a unique nearby maximizer
$y_*(t,x)=y_0(x)+t y_1(x)+O(t^2)$, uniformly on $x\in[0,X]$.

Differentiate Eq.~\eqref{eq:s-uniform-expansion}, divide by $t$, and insert
this expansion:
\[
  0=-\partial_yF_1(y_*,x)+t\partial_yF_2(y_*,x)+O(t^2).
\]
Because $\partial_yF_1(y_0,x)=0$, the coefficient of $t$ yields
\[
  y_1(x)=
  \frac{\partial_yF_2(y_0,x)}{\partial_y^2F_1(y_0,x)}.
\]
Using $y_0^3=6/(1+x)$ gives
\[
  y_1(x)=-\frac{x+4}{6(x+1)}.
\]
Finally,
\[
  t^{-1}y_0
  =\nu^{1/3}\left(\frac6{1+12\nu c_v^2}\right)^{1/3}
  =\left(2c_v^2+\frac1{6\nu}\right)^{-1/3},
\]
while $x=12\nu c_v^2$. Since $\rho_*=t^{-1}y_*$, Eq.~\eqref{eq:s-uniform-ridge}
follows. Setting $x=0$ gives $y_1=-2/3$, and taking $x\to\infty$ in the
matching coefficient gives $y_1\to-1/6$.
\end{proof}
```

**Claim-discipline sentence to retain after this proposition:**

> Proposition S\(\cdot\) proves the optimizer expansion conditional on the stated uniform delayed-renewal information expansion. Establishing that expansion uniformly for a triangular family of hold laws is a separate renewal local-limit problem; absent that result, Eq. (crossover) remains a matched-asymptotic prediction rather than a universal finite-window theorem.

## 6.2 Proof of Theorem 1 — missing-information identity

```latex
\begin{proof}[Proof of Theorem~\ref{thm:identity}]
Let $A_t\in\{0,1\}$ be the predictable indicator that the detector is active,
and let
\[
  L_T=\int_0^T A_t\,dt
\]
be its total active time. Conditional on the detector-state path, registered
events have stochastic intensity $\lambda A_t$. Because the hold law is
independent of $\lambda$, all hold-law factors in the complete-path likelihood
are constant with respect to $\theta=\log\lambda$. Hence the complete-path
score is
\[
  U_T=\frac{\partial}{\partial\theta}\ell_T(\theta)
      =N_T-\lambda L_T.
\]
The process
\[
  M_t=N_t-\lambda\int_0^t A_s\,ds
\]
is a square-integrable martingale with predictable quadratic variation
$\langle M\rangle_T=\lambda L_T$. Therefore
\[
  \mathbb E(U_T)=0,
  \qquad
  I_{\rm path}(\theta;T)=\mathbb E(U_T^2)
  =\mathbb E\langle M\rangle_T
  =\lambda\mathbb E L_T
  =\mathbb E N_T.
\]
Let $\mathcal G_T=\sigma(N_T)$ be the sigma-field generated by the retained
scalar count. Fisher's identity gives the observed-count score
\[
  U_N=\mathbb E(U_T\mid\mathcal G_T).
\]
Applying the law of total variance to $U_T$ yields
\[
  \operatorname{Var}(U_T)
  =\operatorname{Var}\{\mathbb E(U_T\mid N_T)\}
   +\mathbb E\{\operatorname{Var}(U_T\mid N_T)\}.
\]
The first term is the Fisher information in $N_T$. Conditional on $N_T$, the
count itself is fixed, so
\[
  \operatorname{Var}(U_T\mid N_T)
  =\lambda^2\operatorname{Var}(L_T\mid N_T).
\]
Substitution of $\operatorname{Var}(U_T)=\mathbb E N_T$ proves
\[
  I_N(\theta;T)=\mathbb E N_T
  -\lambda^2\mathbb E\!\left[\operatorname{Var}(L_T\mid N_T)\right].
\]
\end{proof}
```

## 6.3 Proof of Theorem 2 — long-window information rate

```latex
\begin{proof}[Proof of Theorem~\ref{thm:rate}]
Rescale time by the mean hold time $\tau$, so that the observation window is
$\nu=T/\tau$, the arrival load is $\rho=\lambda\tau$, and
$\mathbb EB=1$, $\operatorname{Var}B=c_v^2$. After the first active waiting
time, the registered-count process is a delayed renewal process with ordinary
cycle
\[
  C=B+X,
  \qquad X\sim\operatorname{Exp}(\rho),
\]
where $B$ and $X$ are independent. Thus
\[
  \mu_C=\mathbb EC=1+\rho^{-1},
  \qquad
  \sigma_C^2=\operatorname{Var}C=c_v^2+\rho^{-2}.
\]
The first delay changes only $O(1)$ terms. Standard delayed-renewal mean and
variance asymptotics give, locally uniformly in $\rho$,
\begin{align}
  m_\nu:=\mathbb E N_\nu
    &=\frac{\nu}{\mu_C}+O(1)
      =\nu\frac{\rho}{1+\rho}+O(1),\\
  v_\nu:=\operatorname{Var}N_\nu
    &=\frac{\sigma_C^2\nu}{\mu_C^3}+O(1)
      =\nu\frac{\rho(1+c_v^2\rho^2)}{(1+\rho)^3}+O(1).
\end{align}
With $\theta=\log\rho=\log\lambda+\mathrm{const}$,
\[
  \dot m_\nu:=\partial_\theta m_\nu
  =\nu\frac{\rho}{(1+\rho)^2}+O(1).
\]

Under the stated nonlattice renewal local-limit/LAN assumptions, the scalar
renewal count has the same leading Fisher term as the normal experiment with
mean $m_\nu$ and variance $v_\nu$:
\[
  I_N(\theta;\nu)
  =\frac{\dot m_\nu^{\,2}}{v_\nu}+O(1).
\]
Equivalently, this follows from the local normal score expansion; the term
containing $\partial_\theta v_\nu$ and the first Edgeworth correction contribute
only $O(1)$ information, whereas the mean-shift term is $O(\nu)$. Therefore
\[
  \frac{I_N(\theta;\nu)}{\nu}
  \longrightarrow
  \frac{\rho}{(1+\rho)(1+c_v^2\rho^2)}
  =J_\infty(\rho,c_v).
\]
Finally, the complete-path information rate from
Theorem~\ref{thm:identity} is
$\lim_{\nu\to\infty}\mathbb EN_\nu/\nu=\rho/(1+\rho)$. Their difference is
\[
  \frac{\rho}{1+\rho}-J_\infty(\rho,c_v)
  =\frac{c_v^2\rho^3}{(1+\rho)(1+c_v^2\rho^2)},
\]
which is the claimed extensive missing-information rate.
\end{proof}
```

**Required assumption sentence immediately before or after this proof:**

> The local-limit step is invoked for a regular nonlattice delayed-renewal family with a finite \(3+\epsilon\) moment and parameter derivatives controlled on compact load intervals; lognormal, gamma/exponential, and bounded regular hold families satisfy the intended conditions. The theorem is not asserted for infinite-variance, correlated, observed, or load-dependent holds.

## 6.4 Proof of Theorem 3 — cubic optimum and small-jitter expansion

```latex
\begin{proof}[Proof of Theorem~\ref{thm:cubic}]
For fixed $c_v>0$,
\[
  J_\infty(\rho,c_v)
  =\frac{\rho}{(1+\rho)(1+c_v^2\rho^2)}.
\]
Its logarithmic derivative is
\[
  \frac{d}{d\rho}\log J_\infty
  =\frac1\rho-\frac1{1+\rho}
   -\frac{2c_v^2\rho}{1+c_v^2\rho^2}
  =\frac{1-c_v^2\rho^2(1+2\rho)}
  {\rho(1+\rho)(1+c_v^2\rho^2)}.
\]
The function $c_v^2\rho^2(1+2\rho)$ is strictly increasing from zero to
infinity, so there is exactly one stationary point, determined by
\[
  c_v^2\rho_c^2(1+2\rho_c)=1.
\]
Since $J_\infty\to0$ as $\rho\downarrow0$ and as $\rho\to\infty$, that
stationary point is the unique global maximum.

For the small-jitter expansion, put $\varepsilon=c_v^{2/3}$ and seek
\[
  \rho_c=a\varepsilon^{-1}+b+d\varepsilon+e\varepsilon^2+O(\varepsilon^3).
\]
Substitution into
$\varepsilon^3\rho_c^2(1+2\rho_c)=1$ and matching powers gives
\begin{align*}
  2a^3&=1,\\
  a^2(6b+1)&=0,\\
  6a^2d+6ab^2+2ab&=0,\\
  6a^2e+12abd+2ad+2b^3+b^2&=0.
\end{align*}
Hence
\[
  a=2^{-1/3},\qquad
  b=-\frac16,\qquad
  d=\frac{2^{1/3}}{36},\qquad
  e=-\frac{2^{2/3}}{324}.
\]
In particular,
\[
  \rho_c
  =2^{-1/3}c_v^{-2/3}-\frac16
   +\frac{2^{1/3}}{36}c_v^{2/3}
   +O(c_v^{4/3}),
\]
which proves the stated expansion.
\end{proof}
```

## 6.5 Alignment proposition — amended statement and proof

The current raw residual formula is correct when all fixed-flux loads are
interior. For manuscript rigor with box-constrained loads, replace it by the
projected normal-cone residual below. The existing interior formula follows by
setting the normal-cone term to zero.

```latex
\begin{proposition}[Alignment and projected-residual bound]
\label{prop:s-alignment}
Let $\pi_i>0$, $\sum_i\pi_i=1$, and consider
\[
  \max_{\rho}
  F(\rho)=\sum_i\pi_i\ell_iJ(\rho_i)
  \quad\text{subject to}\quad
  \sum_i\pi_i\rho_i=B,
  \qquad \rho_i\in[\underline\rho_i,\overline\rho_i].
\]
Assume $J\in C^2$ on the admitted interval and
$-\ell_iJ''(\rho)\ge m>0$ for every $i$ and every admitted $\rho$.
Let the feasible fixed-flux allocation be
$\widetilde\rho_i=\kappa u_i$.

The fixed-flux allocation is optimal if and only if there exist a scalar
$\beta$ and normal-cone elements
$s_i\in N_{[\underline\rho_i,\overline\rho_i]}(\widetilde\rho_i)$ such that
\[
  \ell_iJ'(\widetilde\rho_i)-\beta-s_i=0
\]
for every $i$. For any $\beta$, choose $s_i$ as a closest normal-cone element
and define the projected KKT residual
\[
  r_i=\ell_iJ'(\widetilde\rho_i)-\beta-s_i.
\]
Then
\begin{equation}
  0\le F(\rho^*)-F(\widetilde\rho)
  \le \frac1{2m}\sum_i\pi_i r_i^2,
  \label{eq:s-residual}
\end{equation}
where $\rho^*$ is the water-filling optimum. If every
$\widetilde\rho_i$ is interior, $s_i=0$ and
$r_i=\ell_iJ'(\kappa u_i)-\beta$, recovering the interior alignment condition.
\end{proposition}

\begin{proof}
The feasible set is convex and compact, and $F$ is strictly concave under the
curvature assumption, so the KKT conditions are necessary and sufficient and
the optimum is unique. The equality multiplier is $\beta$; the box constraints
contribute the normal-cone elements $s_i$. This proves the first statement.

Let $\Delta_i=\rho_i^*-\widetilde\rho_i$. Strong concavity of each
$f_i(\rho)=\ell_iJ(\rho)$ gives
\[
  f_i(\rho_i^*)-f_i(\widetilde\rho_i)
  \le f_i'(\widetilde\rho_i)\Delta_i-\frac m2\Delta_i^2.
\]
Write
$f_i'(\widetilde\rho_i)=\beta+s_i+r_i$ and sum with weights $\pi_i$.
The budget equality gives $\sum_i\pi_i\Delta_i=0$, while the normal-cone
property gives $s_i\Delta_i\le0$. Hence
\[
  F(\rho^*)-F(\widetilde\rho)
  \le\sum_i\pi_i\left(r_i\Delta_i-\frac m2\Delta_i^2\right).
\]
Completing the square term by term,
\[
  r_i\Delta_i-\frac m2\Delta_i^2
  \le\frac{r_i^2}{2m},
\]
which proves Eq.~\eqref{eq:s-residual}. The lower bound follows from optimality
of $\rho^*$.
\end{proof}
```

The main-text alignment sentence remains valid after replacing “residual” by
“projected KKT residual” whenever a realized load touches a frozen bound.

---

# 7. Submission gate after R19

Submission remains **HOLD** until all of the following are complete:

1. dated correction artifact written and hashed;
2. corrected analyzer unit tests, including PAVA `[3,2,1]`, all censoring states, and fixed-strata bootstrap;
3. corrected manuscript numbers and correction disclosure installed;
4. original frozen outputs preserved and linked in provenance;
5. certificate gate language removed and disclosure addendum installed;
6. panel (e), paper-1 scope sentences, and `Gamma=1` label corrected;
7. proof blocks installed and compiled without placeholders;
8. a second independent check verifies the corrected manuscript numbers against the correction JSON.

Once these items pass, the corrected elapsed-time speed verdict is the manuscript source of record; no further endpoint or threshold adjudication is required.