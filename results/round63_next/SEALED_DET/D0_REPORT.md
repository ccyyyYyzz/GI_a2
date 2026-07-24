# D0 — Mechanism and Engine Integrity (dry run) · **PASS**

**R41 §4.6 bar D0** · engine `run_d0.py` (sha `cba967444d177f51`) · data `D0_RESULTS.json` ·
RTX 4060, corrected-shot ledger · wall 35 s · **no confirmatory scene opened** (D0 is scene-independent
geometry + physics; run on the best/mid/floor reference cells + a spectrum/`k_w` sweep).

> **Verdict: D0_PASS.** All three sub-bars clear their frozen thresholds. Kill-tree node 1 is cleared;
> the probe may proceed to the sealed confirmatory run on the coordinator's order.

## D0.1 — Mean derivative to every beyond-band event (the exact wall)

`||P·U_β||_F / ||P||_F` for each claim shell:

| claim | 1.25·k_p | 1.5·k_p | 1.8·k_p | **max** | bar |
|---|---|---|---|---|---|
| rel. norm | 2.0e-16 | 3.5e-16 | 5.7e-16 | **5.68e-16** | ≤ 1e-10 |

The band-limited codes (`|k|≤k_p`) are Fourier-orthogonal to the beyond-band annulus (`|k|>k_p`), so a
beyond-band anomaly changes the **mean** bucket by exactly zero to machine precision — 14 orders of
magnitude inside the bar. This is the first-moment invariance wall that supplies the sentinel's
specificity (and makes FIXED-MEAN / FRESH-MEAN provably blind).

## D0.2 — Two independent Fisher implementations agree ≤ 10%

Engine A (einsum-trace + Schur profiling) vs Engine B (whitened vech-derivative + QR projection,
§3.1). Max relative difference over `{lam_mean, lam_max, tr, d'_energy-spread, d'_matched}`:

| cell (sf, spectrum, k_w/k_p, claim) | lam_mean | lam_max | tr | d'_ES | d'_MT |
|---|---|---|---|---|---|
| 1.0, k⁻², 1, 1.25 (best) | 1.5% | 1.6% | 1.5% | 0.7% | 0.8% |
| 0.3, flat, 1, 1.5 (mid)  | 1.1% | 1.6% | 1.1% | 0.6% | 0.8% |
| 0.3, flat, 1, 1.8 (floor)| 1.0% | 1.7% | 1.0% | 0.5% | 0.8% |
| 0.6, k⁻¹, 2, 1.5         | 0.7% | 0.6% | 0.7% | 0.4% | 0.3% |
| 1.0, flat, 4, 1.8        | 4.8% | 2.7% | 4.8% | 2.4% | 1.4% |

**Worst = 4.8%** (bar ≤ 10%). The two engines were written through fully distinct code paths
(Kronecker/vech whitening + QR vs. per-mode traces + Schur pinv), so the agreement is a genuine
cross-implementation check of `J_B`, `d'`, and `T_det`. The near-null `lam_min` at wide claims
(240 beyond-band dof, many modes at the numerical floor) is **excluded** from the gate: both engines
return ≈0 there and a relative comparison on a ~1e-8 eigenvalue is not physically meaningful (it is the
adversarial worst-mode direction, governed by `λ_min`, reported separately in D2/D5, not a Fisher-map
integrity failure).

## D0.3 — Monte-Carlo shot / covariance vs the `|P|x` ledger ≤ 10%

Three checks against the corrected ledger `V0 = C + R`, `R = diag(|P|x·mean(|P|x)/PHOT)`:

| cell | (a) shot var vs R | (b) C+R self-consistency | (c) physical end-to-end | medium/shot |
|---|---|---|---|---|
| 1.0, k⁻², 1, 1.25 | **1.15%** | **1.01%** | 12.5% † | 62.0 |
| 0.3, flat, 1, 1.5 | **1.15%** | 0.88% | 0.76% | 11.8 |
| 0.3, flat, 1, 1.8 | **1.15%** | 0.88% | 0.76% | 11.8 |

**Gate = max(a, b) = 1.15%** (bar ≤ 10%). PASS.

- **(a) Shot ledger** — frozen medium (`w=1`), full 256-complementary-exposure **physical Poisson**
  simulation: the empirical signed-bucket variance matches `R_i = |P|_i·x·mean(|P|x)/PHOT` to **1.15%**.
  This is the literal R41 D0.3 test and the direct confirmation of the shot correction (the `|P|x`
  throughput, not the signed mean `Px≈0`).
- **(b) C+R self-consistency** — the engine's own second-order generative law (`gen_records`) matches
  the analytic `V0` diagonal to **1.0%** (validates the covariance analytic to sampling noise, all
  contrasts).
- **(c) Physical end-to-end** † — Gaussian second-order medium + physical Poisson. At `sf∈{0.3,0.6}`
  the physical bucket covariance matches `V0` to **0.76%**. At `sf=1.0` (realized contrast 0.70) the
  residual is **12.5%**: this is the **linear-law / non-negativity boundary** — at that contrast the
  physical (non-negative) medium departs from the linearized second-order law by `exp(K_g)−1 ≈ 1.28·K_g`
  plus intensity clipping. This is **reported, not gated**, and is exactly why the lognormal / non-Gaussian
  medium is a **D5 mismatch axis** (Rank-8 non-Gaussian universality), not the primary Fisher-prognosis
  law (cf. `p1_fisher.py`: "lognormal harmonics are a mismatch arm, not the Fisher-prognosis law").

## Bottom line

The mechanism and both engines are sound: the mean-wall is exact (5.7e-16), the two Fisher engines
agree (4.8%), and the corrected `|P|x` shot ledger + the C+R covariance analytic are confirmed against
Monte Carlo (1.15% / 1.0%). The only >10% quantity is the sf=1.0 physical end-to-end covariance, a
documented high-contrast linear-law boundary that is itself a declared mismatch axis. **D0 = PASS.**
