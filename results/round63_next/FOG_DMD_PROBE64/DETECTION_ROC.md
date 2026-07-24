# DETECTION ROC — Monte-Carlo validation of the (corrected) analytic prognosis

**Date:** 2026-07-24 · **Engine:** `detection_roc.py` · **Data:** `DETECTION_ROC.json`,
`DETECTION_ROC_scatter.png`. Local RTX 4060 (~21 s). Fixed: `N=4096`, `M=128`, `n=1e4`.
**Shot term = physical (`|P|·x` throughput)** — the MC is compared against the **CORRECTED**
analytics (`DETECTION_MAP_CORRECTED`), the bug-free MC being the arbiter. Bank-acquisition:
`T_eff` = independent banks, **12.8 ms/bank**.

## Verdict

> **PROGNOSIS VALIDATED.** The empirical detector `d'` matches the corrected analytic `d'` to
> **0.98–1.04× across all 9 (cell × ε) points** (well inside the ±20 % bar), with AUC 0.93–1.00.
> The specificity theorem holds operationally: a **four-way discriminator** (mean, beyond-band-cov,
> amplitude-cov statistics) **cleanly separates H0 / scene-inband / scene-beyond / medium-change**,
> and the beyond-band detector **rejects both in-band scene changes and medium-law changes**
> (separation `d'≈0.1`). ROC is **invariant to a 10 % declared-basis rotation** (AUC 1.00→1.00,
> `d'` 5.30→5.31) — the trace-invariance prediction confirmed. Online cost: **~M²≈16 k MAC/bank**.

## 1. d' validation (empirical vs corrected analytic)

Detector = profiled matched-score statistic `t = ⟨Ŝ₀ − V₀, W⟩`, `W = V₀⁻¹ dV_a V₀⁻¹` (nuisance =
in-band scene + medium amplitude, profiled). `H1 = x + δ`, `δ` energy-spread beyond-band,
`‖δ‖=ε‖x‖`. 250 H0 + 250 H1 records/point; each record = `T_eff` independent banks. (`T_eff` capped
at 1400 for compute; the analytic `d'` is evaluated at the *tested* `T_eff`, so the comparison is
exact-at-that-T.)

| cell | ε | T_eff tested | d'_analytic | d'_empirical | ratio | AUC | T_det (banks) | T_det (s) |
|---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| **best** (k⁻²,k_w1,c1.25,σ0.3) | 1 % | 1400 | 3.46 | 3.61 | **1.04** | 0.998 | 2956 | 37.8 |
| | 2 % | 739 | 5.07 | 5.17 | 1.02 | 1.000 | 739 | 9.5 |
| | 5 % | 118 | 5.09 | 5.25 | 1.03 | 1.000 | 118 | 1.5 |
| **mid** (flat,k_w1,c1.5,σ0.3) | 1 % | 1400 | 2.47 | 2.44 | 0.99 | 0.960 | 5474 | 70.1 |
| | 2 % | 1369 | 5.04 | 4.94 | 0.98 | 1.000 | 1369 | 17.5 |
| | 5 % | 219 | 5.06 | 5.08 | 1.00 | 1.000 | 219 | 2.8 |
| **floor** (flat,k_w1,c1.8,σ0.3) | 1 % | 1400 | 2.22 | 2.23 | 1.00 | 0.934 | 7738 | 99.0 |
| | 2 % | 1400 | 4.32 | 4.27 | 0.99 | 0.999 | 1934 | 24.8 |
| | 5 % | 310 | 4.88 | 5.02 | 1.03 | 1.000 | 310 | 4.0 |

**Every ratio ∈ [0.98, 1.04].** The corrected Fisher engine is sound; the analytic `T_det` table is
trustworthy. (Had the *buggy* analytics been the target, the empirical `d'` would have sat ~1/1.38×
below them — the correction is exactly what closes the gap.)

## 2. Specificity — the four-way discriminator (best cell, ε=5 %)

Three statistics: `t_mean` (mean-channel χ²), `t_cov` (beyond-band matched, **profiled**), `t_amp`
(medium-amplitude matched). Class-separation `d'` between each population and H0:

| population ↓ / statistic → | t_mean | t_cov (beyond) | t_amp (medium) |
|---|:--:|:--:|:--:|
| scene **in-band** | **61.0** | 0.1 | — |
| scene **beyond-band** | −0.1 | **13.5** | 0.3 |
| **medium** (σ_f +20 %) | — | 0.1 | **60.8** |

Reading (see `DETECTION_ROC_scatter.png`): an **in-band** scene change lights only the **mean**
channel (`d'=61`) and leaves the beyond-band cov statistic dark (`0.1`); a **beyond-band** change
lights **only** the covariance statistic (`d'=13.5`) and is invisible to the mean (`−0.1`, the
specificity theorem in action); a **medium-law** change lights **only** the amplitude statistic
(`d'=61`) and is **rejected by the beyond-band detector** (`0.1`) because a σ_f shift lies in the
law-amplitude nuisance direction that the profiling removes. The three classes occupy three
orthogonal statistic axes — a clean **scene-inband / scene-beyond / medium** classifier seed.

## 3. Third class (medium-law change) — separability

`σ_f +20 %` (a covariance-scale change) is **profiled out** of the beyond-band detector
(`t_cov` separation 0.1 = H0-like) and separates strongly on the amplitude statistic
(`t_amp` separation **60.8**). So a medium drift will **not** false-alarm the beyond-band
anomaly detector, and is itself detectable on a dedicated amplitude channel. (A **τ** change would
separate via the *lag-decay* of `S_0..S_4`, a channel not exercised by the iid-bank MC — flagged
for the deployment detector.)

## 4. Robustness — ROC under 10 % declared-basis rotation

Best cell, ε=2 %, `T_eff=800`: detector filter `W` built from a **10 % basis-rotated** declared
medium law, applied to **true-law** data.

| detector filter | AUC | empirical d' |
|---|:--:|:--:|
| matched (true law) | 1.000 | 5.30 |
| **10 % rotated declared law** | 1.000 | 5.31 |

AUC and `d'` are unchanged — the empirical confirmation of the analytic **trace-invariance**: the
aggregate detection statistic does not care about a within-aperture rotation of the declared medium
subspace. Detection inherits (and exceeds) the E5d estimation robustness.

## 5. Online compute cost of the deployed detector

Running one-epoch (per-bank) covariance accumulation `Ŝ₀ ← Ŝ₀ + b_t b_tᵀ`: **M² = 16 384 MAC/bank**;
the final score `⟨Ŝ₀ − V₀, W⟩` is one more M² inner product. The M×M matched filter `W` is
**precomputed offline** from the declared law. At 12.8 ms/bank the detector is trivially real-time
(≈1.3 MFLOP/s of covariance arithmetic).

## Bottom line

The corrected analytic detection prognosis **survives Monte Carlo** (empirical `d'` within 4 %),
the mean-channel specificity is **operationally exact** (a beyond-band anomaly is invisible to the
mean and rejected-as-medium-safe), and the detector is **rotation-robust and cheap**. Combined with
`DETECTION_MAP_CORRECTED` (DETECTION_REGION_FOUND: 2 %-anomaly detection nearly universal, best
**≈6–10 s**; 1 % at the best cells; 0.5 % below floor), the demand-tier row is **validated and
honest** — ready for the R41 ruling and the parallel prior-art scout. No git commit.

## Files
`detection_roc.py`, `DETECTION_ROC.json`, `DETECTION_ROC_scatter.png`;
map: `detection_map.py`, `DETECTION_MAP_CORRECTED.json`, `DETECTION_MAP_CORRECTED.md`.
Superseded (buggy shot): `DETECTION_MAP.json`, `DETECTION_PROGNOSIS.md`.
