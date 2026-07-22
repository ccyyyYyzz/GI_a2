# R31 Pro audit appendix (issue #23 comment 0, raw)

Posted: 2026-07-22T12:37:14Z

---

## Independent audit addendum — correlated-speckle determinant, fixed-pattern floor, and photon-dose law

I independently rederived the R31 mechanism after #23 was posted. The ruling in #23 is substantively correct. The amendments below sharpen the theorem and turn the Fisher gate into a genuine life-or-death test.

### A. The product law is one endpoint of a more general determinant law

Let `V = diag(v)`, and suppose one microcell receives `k` mutually incoherent complex-Gaussian modes:

```text
T_c = (1/k) Σ_{r=1}^k g_cr† V g_cr,
with g_cr iid ~ CN(0, Σ).
```

With conditional Poisson primary avalanches, the no-fire probability is the Gaussian quadratic-form Laplace transform

```text
q(ρ;v) = E[exp(-ρ T_c)] = det(I + (ρ/k) ΣV)^(-k).                 (A1)
```

The proof is the standard complex-Gaussian integral, applied independently to the `k` modes. For one mode, the reciprocal has the principal-minor expansion

```text
det(I + ρΣV)
  = Σ_{I subset [N]} ρ^|I| det(Σ_I) Π_{i in I} v_i.             (A2)
```

Therefore the plain elementary symmetric polynomials `e_j(v)` occur only when the scene-coordinate modes are independent and equally normalized, `Σ = I`. For diagonal known `Σ`, one obtains weighted `e_j(Σ_ii v_i)`. For non-diagonal `Σ`, the coefficients contain cross-coordinate principal minors and are not image `e_j` values.

The log curve identifies modal power sums:

```text
log q(ρ)
  = Σ_{j>=1} (-1)^j ρ^j tr[(ΣV)^j] / (j k^(j-1)).              (A3)
```

Thus the second-order observable is

```text
(1/k) tr[(ΣV)^2] = (1/k) Σ_{i,j} |Σ_ij|^2 v_i v_j,             (A4)
```

not `p_2(v)` unless the optical coherence matrix is diagonal in the pixel basis, or is known well enough to invert. An unknown coherence matrix is a scene-moment confounder.

The coherent-collapse counterexample in #23 is the rank-one quadratic-form endpoint. If a single coherent field is `E_c = a(v)† g_c`, the matrix determinant lemma gives

```text
q(ρ) = 1 / [1 + ρ a(v)† Σ a(v)].                              (A5)
```

For iid propagation coefficients and `a_i = sqrt(v_i)`, this is `q = 1/(1 + ρp_1)`: the complete sweep carries only the ordinary bucket.

**Binding wording:** “fully developed speckle” alone does not imply the product theorem. A single fully developed coherent mode gives the rank-one collapse. The image-power-sum theorem requires mutually incoherent or orthogonal scene-coordinate contributions before microcell saturation.

**P0 should fit three nested physical models**, not merely product versus coherent collapse:

1. rank-one coherent law;
2. diagonal/product law;
3. correlated determinant law with a low-rank or calibrated coherence matrix.

Use both (i) equal-`p_1`, unequal-`p_2` controls and (ii) equal-`p_1,p_2` spatial permutations that alter correlation-weighted cross terms. The latter separates a true pixelwise product law from a correlated determinant masquerading as curvature.

---

### B. `C_eff` is mask- and power-dependent, not one detector constant

For the quenched cell variables

```text
Y_c(ρ,v) = exp[-ρT_c(v)],
q_C      = C^(-1) Σ_c Y_c,
```

define exactly

```text
C_eff(ρ,v) = C^2 Var(Y_c) / Σ_{c,c'} Cov(Y_c,Y_c'),
Var(q_C)   = Var(Y_c) / C_eff.                                 (B1)
```

For iid cells, `C_eff = C`. Optical grain overlap, electrical crosstalk, and common gain drift reduce it. Because `Y_c` changes with both `ρ` and `v`, a scalar measured on one flat-field exposure is not generally transferable across the sweep or mask bank. P1 should report `C_eff(ρ,v)`, or a conservative lower envelope, on the actual discovery masks.

---

### C. There is an irreducible fixed-speckle detectability floor

Write

```text
t       = ρp_1,
β_2     = p_2/p_1^2 = 1/n_eff,
log q   = -t + (β_2 t^2)/2 + O(β_3 t^3).
```

For a fixed cell-mixing realization, the finite-`C` covariance in #23 gives, in the diffuse regime,

```text
Var(q_C)/q^2 ≈ [exp(β_2 t^2)-1]/C_eff ≈ β_2 t^2/C_eff.          (C1)
```

Even after infinitely many repeated gates, fitting the curvature therefore has the approximate quenched floor

```text
SE(β_2-hat)/β_2 ≈ (2/t) sqrt(n_eff/C_eff).                      (C2)
```

For two equal-`p_1` scenes whose `β_2` values differ by a fraction `f`, the best-case fixed-pattern discriminability is

```text
d'_quenched ≈ (f t/2) sqrt(C_eff/n_eff),
f required for d'=3 ≳ (6/t) sqrt(n_eff/C_eff).                 (C3)
```

At `C_eff = 3600`:

| `n_eff` | minimum fractional `β_2` difference for `d'=3`, `t=3` | same, `t=5` |
|---:|---:|---:|
| 200 | 47% | 28% |
| 1000 | 105% | 63% |
| 1800 | 141% | 85% |

The apparent benefit of `t=5` is not free: `q ≈ exp(-5) = 0.0067`, only about 24 un-fired cells at `C=3600`, before nuisance and model error. This is why the no-fire guard in #23 is essential.

Repeated gates beat shot noise but **do not beat this floor** when the cell mixing `W` is fixed. If the speckle or modal realization is independently refreshed between gates, the effective sample size can grow approximately as `G C_eff`; that is the cleanest rescue mechanism, but it must be physically demonstrated rather than silently assumed.

---

### D. The photon budget scales quadratically with scene diffuseness

Let

```text
B = C G t
```

be the mean number of primary-photoelectron opportunities delivered across all microcells and gates for one mask. Converting this to incident optical photons requires the corresponding PDE factor. In the optimistic known-`p_1`, one-mode, nuisance-free approximation,

```text
Var(β_2-hat)/β_2^2 ≈ 4 n_eff^2 [exp(t)-1] / (B t^3).            (D1)
```

For a fractional difference `f` and `d'=3`, minimizing over `t` gives `t = 2.821439...` and

```text
B_3σ, known p_1 ≳ 25.33 n_eff^2/f^2.                            (D2)
```

When `p_1` is a nuisance and the local two-point c-optimal sweep from #23 is used, the inverse-information constant is `5.2997/B`, yielding

```text
B_3σ, nuisance p_1 ≳ 47.70 n_eff^2/f^2.                         (D3)
```

For `n_eff=1000`, these bounds imply roughly `4.8e7` detected primary-photoelectron opportunities per mask for an order-one curvature difference, and `4.8e9` for a 10% difference, before `k>1`, dark counts, crosstalk, power uncertainty, or the quenched floor. The `n_eff^2` law is the central thin-effect warning for diffuse masks.

P2 should therefore report both:

1. the nuisance-profiled finite-dose CRB; and
2. the fixed-pattern floor above.

Passing only one is not sufficient.

---

### E. Prior-art fence additions

1. In the ideal intensity-additive exponential model,

   ```text
   E[T]=p_1,  Var(T)=p_2,  and speckle contrast squared = p_2/p_1^2.
   ```

   Thus `p_2` as a speckle-contrast statistic is classical. The novelty cannot be “speckle reveals a second moment.” The candidate novelty is converting an **unaddressed SiPM occupancy sweep** into mask-indexed Laplace/cumulant constraints and solving the joint spatial inverse problem.

2. The occupancy fence should include linear counting and frequency-moment sketches, not only Flajolet–Martin:
   - Whang, Vander-Zanden & Taylor, “A linear-time probabilistic counting algorithm for database applications,” DOI [10.1145/78922.78925](https://doi.org/10.1145/78922.78925).
   - Alon, Matias & Szegedy, “The space complexity of approximating the frequency moments,” DOI [10.1006/jcss.1997.1545](https://doi.org/10.1006/jcss.1997.1545).

3. The quantized-CS fence should include saturation-aware recovery:
   - Laska, Boufounos, Davenport & Baraniuk, “Democracy in action: Quantization, saturation, and compressive sensing,” DOI [10.1016/j.acha.2011.02.002](https://doi.org/10.1016/j.acha.2011.02.002).

These do not appear to contain the complete SiPM-mask-Laplace inverse construction, but they forbid claiming novelty for occupancy moments, saturation as information, or nonlinear compressed sensing in isolation.

---

### Binding amendment to the final verdict

The conditional GO survives, but I would freeze it more narrowly:

> **GO only for a verified incoherent/modal-mixing implementation in which the coherence operator is diagonal or calibrated, the realistic null-pair `Δβ_2` clears both the `n_eff^2` photon-dose law and the fixed-`C_eff` floor, and the joint Fisher-weighted prior tangent actually closes.**

Immediate kill conditions are:

- a rank-one coherent response at the microcells;
- an uncalibrated correlated determinant response whose curvature is not `M(x^2)`;
- realistic null-pair fractional differences below `(6/t)sqrt(n_eff/C_eff)` before the no-fire fraction enters the fragile `<1%` regime;
- or a nuisance-profiled dose requirement incompatible with equal-time/equal-photon GI.

This strengthens, rather than reverses, #23: the mechanism is not a renaming, but its physically useful domain is controlled by coherence rank, effective scene support, and whether the speckle ensemble is refreshed.