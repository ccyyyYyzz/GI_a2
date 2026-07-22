# Numerical verification of the Part II scattering theorems (O1-P and O4-A)

*Generated 2026-07-22 03:05:45Z | 3.11.5 | numpy 1.24.4 | runtime 64.2s | script `code/scatter/verify_o1_o4.py`*

Normative source: `paper3_scatter/main_scatter.tex` (Thm O1 / Cor O1-P, eqs. O1.2/O1.4; Thm O4-A, eqs. O4.1-O4.4; Sec. 6 randomization proposition). Plan sketch: `paper3_scatter/NOTES.md`.

## Verdicts

| Block | Claim verified | Verdict |
|---|---|---|
| 1 | O1-P hidden-exposure identity (O1.4) | **PASS** |
| 2 | O4-A Schur / orthogonality (O4.2-O4.4) | **PASS** |
| 3 | Randomization concentration + named obstacle (Sec. 6) | **PASS** |

## Block 1 - O1-P hidden-exposure identity

Exact instance: K=8 object, N=64 frames, nonnegative bounded patterns, positive object; conditionally-Poisson buckets with mean count 24.9. Gain a_n=exp(l_n), OU log-gain (sd=0.3, E[a]=1.0004) discretised onto a G=31-point Tauchen HMM grid so the posterior and all expectations are **exact** via forward-backward.

The exact posterior machinery is validated up front: the backward second-moment accumulator for `Cov(M^T a | Y)` matches a brute-force enumeration of all G**N state paths on a tiny HMM to relative error **1.0e-13**.

`I_Y(x)` computed two independent ways over a common batch of S=1500 bucket records (common random numbers):

* **(a)** the O1.4 identity `M^T diag(E a/s) M - M^T E_Y[Cov(a|Y)] M`, with the loss term `M^T E_Y[Cov(a|Y)] M = E_Y[Cov(G|Y)]` computed as the **exact per-record posterior covariance** of the gain-weighted exposure vector `G = M^T a` (eq. O1.5) via a backward second-moment accumulator on the HMM (no cancellation of large second moments);
* **(b)** brute-force Fisher `E_Y[ g g^T ]`, `g = grad_x log p(Y|x)` by central finite differences of the exact HMM marginal log-likelihood (step 0.0001).

The hidden gain removes a large information fraction here (missing trace / complete-Fisher trace = 0.86): a 30% correlated multiplicative gain is nearly degenerate with object scale, so this is a stringent regime for the identity.

| Quantity | Value | Target |
|---|---|---|
| Rel. Frobenius `||I_Y^a - I_Y^b|| / ||I_Y^a||` | 5.12e-02 | MC floor |
| Max per-entry rel. discrepancy | 1.22e-01 | MC floor |
| Median per-entry rel. discrepancy | 4.29e-02 | MC floor |
| Paired CRN max\|mean diff\| (per entry) | 9.802e-03 | within few x SE |
| **Paired max z-score** (\|mean\|/SE) | **2.19** | O(1) - **the key test** |
| FD observed score vs analytic O1.1 score (rel.) | 9.52e-10 | ~FD trunc. |

The decisive test is the paired common-random-number statistic: per record it forms `g g^T - (M^T diag(E a/s) M - Cov(G|Y))`, whose expectation is *exactly zero* under O1.4. Since both methods are unbiased for `I_Y`, the raw Frobenius gap is pure Monte-Carlo noise; the max z-score (largest per-entry `|mean|/SE` over the 64 matrix entries) tests the identity at the MC floor. A z-score of O(1) confirms O1.4; a systematic model error would show large z. The observed-score law O1.1 is separately confirmed to finite-difference precision.

**Loss term PSD and gauge checks:**

* every per-record `Cov(G|Y)` is PSD (min eigenvalue over sampled records = 2.81e-02), and the averaged missing term `M^T E_Y[Cov(a|Y)] M` has min eigenvalue 4.242e-02 - **PSD confirmed**.
* deterministic-gain limit (G=1 point-mass prior): missing term max-abs = 1.14e-12 (machine zero), and the FD Fisher matches the analytic complete Fisher `M^T diag(a0/s) M` to rel. 1.81e-02 - **the loss vanishes exactly when the gain is deterministic**.

**OU correlation-time sweep** (loss depends on the full posterior covariance, not a scalar variance - consequence (1) of Sec. 4):

| corr. time t_c (frames) | lag-1 corr | missing-term trace | missing-term op-norm |
|---|---|---|---|
| 0.5 | 0.135 | 6.0493 | 5.4248 |
| 2.0 | 0.607 | 7.0322 | 6.5286 |
| 8.0 | 0.882 | 7.2997 | 6.9994 |

**Pattern-ordering dependence** (consequence (2)): the loss rotates gain uncertainty through `M^T(.)M`, so it is order-dependent in general. For the i.i.d. patterns used here the rows are nearly exchangeable, so a random permutation barely moves the trace (identity 7.042 vs permuted 7.015) - the effect is real but small for unstructured patterns. The *strong* ordering effect is exhibited in Block 2(iii), where structured schedules on a fixed pattern set change the Schur loss by many orders of magnitude.

## Block 2 - O4-A Schur complement and nuisance-orthogonality

### (i) Schur formula O4.2 vs numerically inverted joint information

Across 24 random instances sweeping (N,K,p), the Schur complement `I_xx - I_xB I_BB^-1 I_Bx` equals the reinverted xx-block of `J^-1` to worst relative error **9.12e-16** (machine precision).

### (ii) Exact moment-satisfying chop and quadratic perturbation growth

A mirror-paired (time-symmetric) chopping schedule is constructed for p=2 **AC drift modes** H=[t, t^3] (odd in time). Mirror pairs annihilate every odd temporal moment exactly:

* cross block `|M^T R^-1 D_s H|_max = 0.00e+00` (machine zero),
* Schur-complement loss tr(I_xx - I_x|gain) = 0.00e+00 - **zero to machine precision**.

**Gauge caveat (honest, and a genuine finding).** The *constant* (DC) drift mode cannot be orthogonalised by any schedule: dotting its moment with x gives `sum_n s_n^2/sigma_n^2 = 2819.711 > 0` (here `I_xB[:,DC]·x = 2819.711`), which is strictly positive whenever the object is visible. A constant multiplicative gain is the scale/gauge direction, degenerate with object amplitude - exactly the gauge-singularity flagged in the manuscript's MG-scope remark and settled only by Part I identifiability. Hence O4.4 is satisfiable for centered/AC drift modes but never for the DC mode; the exact-zero demonstration uses the realisable (AC) modes, as the theory requires.

Breaking the mirror symmetry by shifting the +side times by delta grows the loss quadratically:

| delta | loss tr(I_xx - I_x\|gain) |
|---|---|
| 1e-01 | 2.241e+00 |
| 3e-02 | 2.312e-01 |
| 1e-02 | 2.614e-02 |
| 3e-03 | 2.361e-03 |
| 1e-03 | 2.625e-04 |
| 3e-04 | 2.363e-05 |
| 1e-04 | 2.626e-06 |

log-log slope = **1.984** (theory: 2 - loss is `I_xB I_BB^-1 I_Bx` with `I_xB = O(delta)`).

### (iii) Three schedules on the same pattern set

Same fixed multiset of matched-pair patterns (carrier direction drifts with index), only the time schedule changes; loss = tr(I_xx - I_x|gain), modes H=[t,t^3]:

| schedule | loss |
|---|---|
| naive ordered | 1.7699e+03 |
| random permutation (mean of 50) | 6.3593e+01 (min 2.65e+00, max 2.32e+02) |
| paired-chopped | 0.00e+00 |

Ordering: **ordered (1.77e+03) >> random (6.36e+01) > paired-chopped (0.00e+00)** - ratios ordered/random = 27.8, random/chopped = 6.4e+31. This is exactly the ruling's predicted ordering.

## Block 3 - Randomization concentration and the named obstacle

Averaged cross block `I_xB = (1/N) sum_n c_n h_{pi(n)}^T` (intensive normalisation, so the O(sqrt(log/N)) rate is a decay in N).

### Centered bounded carriers - concentration holds

| N | op-norm mean | op-norm q95 | predicted B·sqrt(log((K+p)/delta)/N) |
|---|---|---|---|
| 64 | 0.1196 | 0.1671 | 0.2710 |
| 256 | 0.0588 | 0.0855 | 0.1355 |
| 1024 | 0.0292 | 0.0421 | 0.0678 |

Fitted log-log slope of q95 vs N = **-0.497** (theory -0.5). The empirical op-norm tracks the predicted `sqrt(log/N)` curve: **concentration confirmed** for centered bounded carriers, as Proposition (Sec. 6) requires.

### Named obstacle - uncentered scene-dependent carrier (OPEN LEMMA)

With nonnegative patterns and a **bright structured scene**, the carrier `c_n = s_n m_n` is elementwise nonnegative, so its mean `c_bar` does not vanish. Random permutation cannot center the cross block: `E_pi[I_xB] = c_bar h_bar^T` (h_bar != 0 for the constant drift mode). The residual bias therefore does **not** decay with N:

| N | residual bias ||c_bar·h_bar^T|| | realised op-norm (mean) | ||c_bar|| |
|---|---|---|---|
| 64 | 30.9152 | 30.9286 | 29.2346 |
| 256 | 32.3555 | 32.3596 | 30.6710 |
| 1024 | 32.7614 | 32.7618 | 31.0741 |

**empirical face of the OPEN LEMMA (Remark rem:rand-obstacle): scene-dependent carrier s_n m_n destroys centering; closing the lemma requires the stationary-carrier hypothesis of Part I.** The residual bias is flat in N (it is set by the scene, not the sample size), which is the empirical signature that the randomization proposition is only a theorem under the bounded/centered hypothesis; for arbitrary nonnegative GI patterns the product `s_n m_n` breaks centering. This is carried in the manuscript as an explicit open-lemma Remark and must remain until the stationary-carrier hypothesis closes it.

## Honest caveats

* Block 1 (a)/(b) agreement is Monte-Carlo-limited (both estimators are unbiased for `I_Y`); the report quantifies the residual as a z-score against the MC standard error, not as a machine-precision identity. The exact per-record identity that *is* verified to FD precision is the observed-score law O1.1 (`grad log p(Y|x) = M^T(D_s^-1 Y - E[a|Y])`).
* The OU path is discretised onto a finite grid (Tauchen); every 'exact' statement is exact for the grid chain, which is the same generative model used to draw Y, so no model mismatch enters the comparison.
* Block 2's exact-zero chop uses AC (odd) drift modes; the DC/gauge mode is structurally non-orthogonalisable and is reported as such rather than forced to zero.
* Block 3 uses the intensive (1/N) normalisation of the information so the concentration rate is a decay; the absolute (summed) cross block grows like sqrt(N) as expected.

