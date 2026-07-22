# R34 ruling (GitHub issue #26, raw)

Title: R34
Posted: 2026-07-22T16:16:22Z

---

# R34 — Final campaign ruling: GO with preregistered simulation-calibrated inference

**Reference question:** [`docs/ROUND63_GPT_ROUND34_QUESTION.md`](https://github.com/ccyyyYyzz/GI_a2/blob/81c031b/docs/ROUND63_GPT_ROUND34_QUESTION.md) at commit [`81c031b`](https://github.com/ccyyyYyzz/GI_a2/commit/81c031b).  
**Frozen feasibility evidence:** [`results/round63_next/DUAL_LEDGER_PROBE/`](https://github.com/ccyyyYyzz/GI_a2/tree/95949da/results/round63_next/DUAL_LEDGER_PROBE) at [`95949da`](https://github.com/ccyyyYyzz/GI_a2/commit/95949da).

## One decision

> **Choose (a): CAMPAIGN GO with a preregistered simulation-calibrated interval construction. Bar 4 becomes a binding confirmatory endpoint. The feasibility probe remains scored 6/7; its failed intervals are not retroactively repaired.**

This is the right disposition because the failure is narrow and diagnosable. DLGI already passed the hard information and systems tests: identical acquisition ledgers; model/gauge residuals; pilot-free `t_c` and CV risk within `1.5 x` of the pilot arm and usually better; scene noninferiority; exact Canonical-Confusion Reciprocity; and honest edge mapping. The missing item is trustworthy finite-sample uncertainty for `t_c`. That is an inference-procedure deficiency, not evidence that the record lacks medium information.

Option (b) is too weak for the intended flagship claim. RMSE-only `t_c` estimates plus CV intervals would support a competent method note, but not “two certified products.” Option (c) is too harsh: archiving after six substantive passes because a skewed finite-sample statistic needs calibrated critical values would confuse a repairable statistics problem with a failed sensing mechanism.

The repair is allowed **once**, in a new confirmatory campaign, because it can be specified without seeing confirmatory outcomes. No third post-hoc interval repair is authorized.

---

# 1. The interval remedy that is authorized

A table of ad hoc quantiles indexed by the observed point estimate is **not** sufficient; selecting the calibration cell using the same noisy estimate can itself break coverage. Freeze a simulation-calibrated test inversion.

Parameterize

```text
eta = (log t_c, log CV).
```

For every candidate `eta_0` on the declared grid and photon-SNR stratum, compute the frozen full-pipeline profile statistic

```text
W(Y;eta_0)
 = 2 { ell_joint(eta_hat, x_hat) - ell_joint(eta_0, x_hat_eta0) }.
```

The exact likelihood approximation may remain the frozen log-Whittle/profile construction from `95949da`, but every calibration replicate must regenerate the gain path and bucket noise and rerun the **entire scene+medium pipeline**. Holding the reconstructed scene fixed is forbidden.

Before confirmatory data are generated:

1. simulate an independent calibration bank at every declared `(t_c,CV,SNR)` cell and frozen scene stratum;
2. estimate the 95% critical value `c_0.95(eta_0)` using a conservative Monte Carlo order statistic;
3. freeze a monotone upper-envelope interpolation in `(log t_c,log CV,log SNR)`;
4. commit the critical-value table, interpolation rule, code hash and random seeds;
5. form the confidence region by Neyman inversion:

```text
C_0.95(Y) = { eta_0 : W(Y;eta_0) <= c_0.95(eta_0) }.
```

Report the projections of this joint region as the `t_c` and CV intervals. This preserves the candidate-parameter indexing needed for model-based coverage; it does not choose a correction after seeing `eta_hat`.

These intervals must be called **simulation-calibrated, model-conditional intervals under the validated scalar-gain/OU model**. They are not distribution-free guarantees.

Coverage alone is not enough. A trivially huge interval would pass. The confirmatory campaign therefore freezes both a coverage bar and an informativeness bar.

---

# 2. Confirmatory campaign skeleton

## 2.1 Freeze and separation

Use three sealed objects:

- **calibration bank:** used only to construct `c_0.95(eta)`;
- **confirmatory bank:** fresh scenes and seeds used only for final endpoints;
- **edge bank:** used only to display failure boundaries.

No scene, gain path, Poisson seed or bootstrap seed may cross banks. Commit the calibration table before opening the confirmatory bank. No tuning, regridding, estimator switching or interval repair is permitted after that commit.

Use at least `5,000` full-pipeline calibration replicates per primary parameter/SNR cell and at least `1,000` independent confirmatory records per cell, balanced across a fresh cohort of at least `12` scenes not used in the DLGI feasibility probe. This gives useful precision for a 95% coverage endpoint and prevents another 204-record ambiguity.

## 2.2 Frozen primary grid

Retain the difficult cells; do not remove the locations that failed in the probe.

```text
t_c / Delta in {16, 32, 64}
CV          in {0.05, 0.15, 0.40}
photon scale in {0.5, 1, 2} x the frozen nominal 2200 counts/exposure
```

The primary claim therefore spans 27 cells. Keep `t_c/Delta = 2` as the fast-drift edge audit and add one slower-than-64 cell as the slow-drift audit. Edge cells cannot rescue or kill the interior result; they define the published failure map.

The scalar-per-frame gain model, `E[a]=1` gauge, exposure count, object dimension, pattern multiset and noise law are frozen exactly before generation.

## 2.3 Arms

1. **Pilot-free DLGI:** same ordinary bucket record supplies the scene and medium estimates.
2. **Pilot-interleaved baseline:** the frozen strongest 5% pilot design, with identical total exposures, photons and duration; pilot rows replace scene rows and their scene cost remains charged.
3. **Oracle monitor:** the gain path observed at every exposure; nondeployable information ceiling.
4. **Plain-linear camera:** reconstruction from the byte-identical pilot-free bucket record, with no medium correction, for the scene noninferiority endpoint.

For pilot-free schedule tests, use exactly

```text
{paired, random, ordered}
```

with the same pattern multiset, photons, exposure count and duration. Only temporal order may change. Preserve the raw record; no differencing step may discard the common-mode medium information.

## 2.4 Preregistered `(A,B,K)` predictions

Before generating confirmatory data, compute for every schedule and cell:

```text
A = scene oracle ledger,
B = medium oracle ledger,
K = A^(-1/2) C B^(-1/2),
J_x     = A^(1/2)(I-KK^T)A^(1/2),
J_theta = B^(1/2)(I-K^TK)B^(1/2).
```

Freeze the predicted schedule ordering for scene risk, medium risk and normalized determinant retention. The empirical ordering may not be selected after seeing results. The paired schedule is the primary acquisition because the feasibility probe found it best for both ledgers, exactly as the reciprocity theorem predicts.

---

# 3. Primary endpoints and binding kill bars

All bars below must pass. There is no averaging away a failed parameter cell.

## Bar C1 — calibrated `t_c` coverage

For every primary `(t_c,CV,SNR)` cell, the empirical coverage of the projected nominal 95% `t_c` interval must lie in

```text
[0.92, 0.98].
```

Report exact binomial intervals and multiplicity-adjusted diagnostics. The three previously failed cells—`(16,0.05)`, `(16,0.40)` and `(64,0.05)` at nominal SNR—are mandatory named stress cells.

**Kill:** any primary cell below `0.92`, or any unreported/nonfinite interval rate above 5%.

## Bar C2 — interval informativeness

At every primary cell:

- median log-scale interval width for `t_c` must be at most `1.5 x` that of the pilot baseline;
- median CV interval width must be at most `1.5 x` pilot;
- intervals must be bounded and connected in at least 95% of records.

Also report oracle-relative width. This prevents calibration from passing through vacuous intervals.

**Kill:** coverage passes only by producing intervals that fail these width/boundedness bars.

## Bar C3 — point precision

Pilot-free RMSE for both `t_c` and CV must remain no worse than `1.5 x` the strongest pilot arm at every primary cell. Report the oracle ratio and signed bias separately.

**Kill:** any cell exceeds `1.5 x`, or the previously observed advantage disappears because the calibration procedure changes the estimator rather than only its uncertainty statement.

## Bar C4 — scene noninferiority

Against the plain-linear reconstruction from the byte-identical record:

```text
PSNR loss <= 0.2 dB
and task/Fisher loss <= 5%
```

at every primary cell, with cohort-level confidence intervals. Pilot-free improvement is welcome but not required.

**Kill:** either noninferiority margin is crossed.

## Bar C5 — physical model and gauge

Repeat held-out innovation, residual-whiteness, scalar-gain adequacy and scale-gauge tests on the fresh scenes and all SNR strata.

**Kill:** systematic pattern-dependent residuals, gauge sensitivity, or a spatially varying medium operator large enough to invalidate the scalar ledger.

## Bar C6 — reciprocity and schedule prediction

Numerically verify

```text
det(J_x)/det(A) = det(J_theta)/det(B)
```

and the shared canonical-confusion spectrum to numerical tolerance. The preregistered `(A,B,K)` schedule ordering must predict the measured ordering for both ledgers, not merely correlate post hoc.

**Kill:** theorem implementation fails, or empirical schedule behavior contradicts the frozen prediction on more than a predeclared 10% of primary cells.

## Bar C7 — identical-ledger discipline and edge honesty

Byte-audit photons, exposures, duration and pattern multiset. Publish the fast- and slow-drift edge failures and the weak-CV/noise boundary.

**Kill:** hidden pilot/reference information, unequal resource accounting, or removal of a failed edge after unblinding.

Failure of any confirmatory bar ends the flagship claim. No additional repair round is authorized. The surviving output would be a shorter theorem-backed method paper with RMSE precision, CV intervals and an explicit `t_c` interval-calibration failure.

---

# 4. Claim ruling

Until the campaign passes, the authorized statement is:

> **DLGI extracts a competitive point estimate of medium correlation time and fluctuation strength from the same bucket record used for scene reconstruction, with no added acquisition and no detected scene penalty; model-calibrated `t_c` uncertainty remains under confirmation.**

Do not use “two certified products” yet.

If all confirmatory bars pass, restore and sharpen the R33 headline:

> **One bucket stream, two model-certified products: pilot-free ghost imaging reconstructs the scene and estimates the medium’s correlation time and fluctuation strength with calibrated intervals and no additional measurements, while matching or beating a pilot monitor under the same acquisition ledger.**

The words **model-certified** and the validated scalar-gain operating domain are mandatory. DLGI is not a universal replacement for DCS, a reference arm, or a monitor when the medium changes the spatial transfer operator.

---

# 5. Manuscript placement

> **Reserve DLGI as the lead of the METHOD act in the operator’s one-flagship manuscript; do not split it into a standalone short paper before confirmation.**

At 6/7 bars it is the strongest surviving method because it has a compact theorem, an immediately legible systems statement, zero acquisition overhead, near-pilot medium risk and scene noninferiority. It is also stronger as the culmination of the manuscript’s arc—first quantify hidden-channel information, then show the saturation boundary, finally reuse one record to deliver a second scientific product—than as an isolated “joint smoother” paper.

This lead position is conditional, not honorary. Passing the calibrated-interval campaign promotes DLGI to the METHOD act and licenses “two model-certified products.” Failing it demotes DLGI to a shorter theorem/method note or a substantial section in the materials bank; it must not lead the flagship on RMSE alone.

## Final ruling

> **GO once, under option (a). Freeze a simulation-calibrated Neyman construction before confirmation, require both coverage and useful width, rerun every R33 bar on fresh scenes and a photon-SNR grid, and permit no further repairs. DLGI conditionally leads the flagship METHOD act.**