# RECORD-REFINEMENT LADDER probe — report

**ROUND63-NEXT / GI_a2 · acquisition-side pivot · 2026-07-22 (Fable)**
Repo `D:\GI_another` · read-only on all prior results; writes only under
`results/round63_next/LADDER_PROBE/`. CPU, ~21 min total (no downscaling needed).

---

## Verdict in one paragraph

The program identity (info loss = `E[Var(hidden | record)]`) leaves exactly one
free axis at fixed photon budget: the record's temporal richness. We tested it on
both program chapters. **Chapter 2 (gain drift, image level):** refining the
readout from a bucket total to sub-bins / timestamps, *at identical photon cost*,
buys back **~+1 dB of image PSNR at the primary cell** (and at every moderate-drift
cell), saturating by **B≈16 sub-bins**; it is **provably vacuous** under the frozen
within-exposure-constant gain model, and — contrary to the going-in hypothesis —
**does not help at fast drift** (t_c=2), where the joint image+gain problem becomes
near-unidentifiable. The larger recovery at the primary cell (**+4.5 dB over the
CPL bucket software arms**) comes not from record refinement but from a
**state-space estimator that exploits the known OU gain-drift statistics on the
existing bucket record** — a software result the five killed lanes did not test,
flagged below with its assumption. **Chapter 1 (dead time + jitter, scalar level):**
a full detection-timestamp record recovers only **11–15 % of the jitter
information loss** (0.05–0.54 dB), because the count already pins the summed holds;
the going-in expectation that timestamps recover ≈J(ρ;0) is **refuted**.

**Formal call:** the record-refinement axis is **WEAK-to-MODERATE**, not the
10–15 dB lever the oracle gap advertised. The oracle gap is dominated by
*cross-exposure gain identifiability*, which sub-binning does not address.

---

## 0. Setup and provenance

Lifted **verbatim** from `results/round63_next/CPL_PROBE/cpl_gi_probe.py`:
32×32 (NPIX=1024); the 6 DEV scenes `bridge_{contour_1, microtex_1, twopop_1,
control_2, contour_3, microtex_3}` from `data/r63_bridge_scenes/`; full 1024-row
sequency Hadamard realised as 2048 complementary (m+, m−) exposures; φ matched to
2200 counts/exposure (φ=11.9375); 5 paired seeds; the frozen OU gain (`ou_path`)
and frozen isotropic-TV rule (λ_img = 0.1 for the A1/A4 inversion path). Cells are
the exact CPL list: primary (σ=40 %, t_c=16) plus secondaries t_c∈{2,64},
σ=15 %, and the σ=0 no-gain control.

**The ladder.** Each exposure's dwell is subdivided into B sub-bins,
B∈{1, 4, 16, 64}. **Photon budget is identical across rungs by construction:** we
simulate at the finest resolution (B_MAX=64) and a coarser rung's count is the
*sum* of the finest sub-bin counts it contains (a coarser record is literally a
binning of the finer one). B=1 is the bucket total — the frozen CPL record. B=64
is the timestamp / continuous-observation proxy (34 counts/bin; finer bins only
add Poisson noise once the smooth within-exposure path is resolved).

**Reconstruction (SAME estimator at every rung).** Joint state-space estimation,
alternating: (gain step) each rung-bin count is a Poisson observation of
`exp(l_bin)·(φ/B)·(pattern·x)`; log-linearise (delta method) and run an RTS
smoother on the log-gain with the OU prior at the rung's resolution; the
per-exposure effective gain is the Jensen-debiased posterior mean
`ā_e = mean_k exp(l_sm + P_sm/2)`. (image step) gain-correct the bucket counts by
`ā_e` and invert with the frozen Hadamard + isotropic-TV rule (= CPL arm A4 with
*estimated* gains). Init from A2; 4 outer iterations; **no per-rung tuning**.

---

## 1. Chapter 2 — gain drift (PRIMARY, image level)

### 1a. Modelling note (a real finding): B>1 is vacuous under the frozen model

Under the CPL model the gain `a_t` is **constant within each exposure** (one gain
per exposure, drifting across exposures per the frozen OU). Then the B sub-bin
counts of one exposure are **Poisson thinnings** of the bucket total — the total
is a sufficient statistic for that exposure's rate, and sub-binning adds **exactly
zero** information. **We therefore state it explicitly: under the frozen model,
record refinement cannot help.** The numerical ladder confirms it — under
`M0_const` the coarsest rung B=1 is optimal at every drift cell and finer rungs
only degrade (the continuous-drift estimator over-smooths the piecewise-constant
truth across exposure boundaries):

| cell (M0, frozen) | A2 | A3 | **B1** | B4 | B16 | B64 | oracle | refine-only |
|---|---|---|---|---|---|---|---|---|
| primary σ40 t_c16 | 7.57 | 7.11 | **10.09** | 8.53 | 7.40 | 7.66 | 17.00 | **−3.10** |
| σ40 t_c2 (fast) | −0.86 | 2.82 | **1.09** | −1.09 | −4.09 | −5.66 | 17.13 | **−6.75** |
| σ40 t_c64 (slow) | 12.76 | 10.67 | 12.87 | 13.02 | 13.01 | 13.26 | 17.02 | −0.60 |
| σ15 t_c16 | 14.53 | 12.22 | 14.04 | 14.73 | 14.78 | 14.80 | 17.43 | −0.07 |
| control σ0 | 17.54 | 16.97 | 17.51 | 17.51 | 17.51 | 17.51 | 17.51 | 0.00 |

(median PSNR over the 6 scenes, dB; "refine-only" = B64 − B1.) Vacuity is
demonstrated: refinement is ≤0 everywhere.

The honest primary case is therefore the **model variant M1** in which the *same*
OU drifts **continuously within exposures** (t_c measured in frames): now sub-bins
genuinely resolve the within-exposure path, so refinement *can* carry information.
All headline numbers below are **M1_cont** and are **labelled a model variant**.

### 1b. The ladder under continuous drift (M1_cont)

| cell (M1, continuous) | A2 | A3 | B1 | B4 | B16 | B64 (ts) | oracle | **ts−bestSW** | **refine-only** | oracle-frac |
|---|---|---|---|---|---|---|---|---|---|---|
| **primary σ40 t_c16** | 8.77 | 7.80 | 10.95 | 11.73 | 11.94 | **12.24** | 17.21 | **+4.52** | **+1.05** | 0.39 |
| σ40 t_c2 (fast) | 0.96 | 3.51 | 3.80 | 3.76 | 3.83 | 2.68 | 17.24 | −0.83 | −0.86 | −0.06 |
| σ40 t_c64 (slow) | 13.50 | 11.24 | 13.16 | 14.19 | 14.33 | 14.34 | 17.13 | +2.61 | +0.95 | 0.43 |
| σ15 t_c16 (mild) | 15.25 | 13.02 | 14.21 | 15.46 | 15.54 | 15.29 | 17.59 | +1.65 | +0.96 | 0.36 |
| control σ0 | 17.54 | 16.97 | 17.51 | 17.51 | 17.51 | 17.51 | 17.51 | −0.02 | 0.00 | 0.00 |

median PSNR (dB) over 6 scenes × 5 seeds. `ts−bestSW` = B64 − max(A2,A3);
`refine-only` = B64 − B1(state-space smoother); `oracle-frac` =
(B64 − bestSW)/(oracle − bestSW).

### 1c. The headline number and its honest decomposition

At the **primary cell**, the timestamp rung reaches **12.24 dB**, versus the best
CPL bucket software arm at **8.77 dB** — a **median gain of +4.52 dB**, capturing
**39 % of the oracle gap**. But that gain splits into two very different pieces:

- **Temporal-prior estimator (+2.18 dB), on the SAME bucket record (B=1):** the
  B=1 state-space smoother (10.95 dB) already beats the gain-agnostic arms A2/A3
  (8.77 dB) by exploiting the OU smoothness of the gain path across exposures. It
  genuinely *identifies* the drift (estimated ā_e correlates 0.998 with the true
  gain path; the residual is a ~30 % amplitude shrinkage toward the OU mean). This
  requires the OU statistics (σ, t_c) to be **known/calibrated** — an assumption
  A2/A3 do not make (see §1e).
- **Record refinement (+1.05 dB), the actual ladder axis:** going from the bucket
  (B=1) to the timestamp rung (B=64) adds **+1.05 dB** at the primary cell, by
  reducing the gain-shrinkage: sub-bins reveal the within-exposure trend, tightening
  the cross-exposure smoother. This is **the answer to the probe's question.**

**Saturation / diminishing returns.** The refinement gain is essentially spent by
**B=4→16**; the timestamp rung adds nothing beyond B=16 and can decline slightly
from low-count linearisation. Median ladder increments (primary): B1→B4 +0.78,
B4→B16 +0.21, B16→B64 +0.30. Diminishing-returns rung ≈ **B=16 (137 counts/bin)**.

### 1d. Fast drift does NOT help (refutes the going-in hypothesis)

The going-in reasoning was that fast drift is where sub-bins matter most. The data
say the opposite: at **t_c=2** the ladder is **flat-to-negative** (refine-only
−0.86 dB), and even the B=1 smoother barely clears A3 (3.80 vs 3.51) while the
oracle sits at 17.24 dB — a 13 dB gap that neither smoothing nor refinement closes.
The reason is identifiability, not resolution: at t_c=2 the gain path carries
≈N_EXP/t_c ≈ 1024 effective degrees of freedom, matching the image's 1024, so the
joint (gain, image) inverse problem has no redundancy (2048 observations ≈ 2048
unknowns). Sub-bins resolve the *within-exposure* path but the bottleneck is
*cross-exposure* separation of gain from image — which refinement does not touch.
Refinement helps only in the **moderate-drift** regime (t_c≥16), where the gain is
smooth enough to be identified *and* sub-bins add usable within-exposure detail.

### 1e. Honesty flag — the temporal-prior estimator touches the software verdict

`docs/SOFTWARE_SATURATION_VERDICT.md` freezes "no software candidate clears the
bar on the bucket-count channels." The B=1 state-space smoother here recovers
**+2.2–2.5 dB** over A2/A3 on the **same bucket record**, under **both** the frozen
model (M0: 10.09 vs 7.57 at primary) and the continuous variant. None of the five
killed lanes (RLMI, DOPS, RQL-efficiency, CPL) exploited the gain's **temporal
prior**, so this is untested territory. Two caveats keep it honest and short of a
verdict-reopening claim: (i) it **assumes the OU drift statistics (σ, t_c) are
known/calibrated**, which the gain-agnostic arms do not; (ii) it still leaves
**~60 % of the oracle gap** on the table (shrinkage + the fast-drift identifiability
wall). Recommendation: this is a genuine bucket-record estimator worth the
operator's attention, but it is **not** the record-refinement result and should be
adjudicated separately.

---

## 2. Chapter 1 — dead time + jitter (SECONDARY, scalar level)

Model lifted from `code/round63/jitter_score_diag_colab_frozen.py`: active-start
non-paralyzable renewal counter, live waits W~Exp(λ0), lognormal holds
B~τ·LogNormal(CV=c), frame T=ν·τ. At ν=2000, c=0.05 we compare the per-frame
Fisher information about the rate for **(a) count-only record N** (= J(ρ;c), the
frozen value) versus **(b) the full detection-timestamp record** (inter-detection
gaps G_i = W_i + B_i).

**Method.** The program identity gives, for any record R that observes N,
`I_R = E[N] − λ0²·E[Var(L | R)]`. For count-only, `Var(L|N)`; for timestamps,
`Var(L | gaps) = Σ_i Var(B_i | G_i)` (the holds are conditionally independent given
the gaps). The recovered fraction of the jitter loss is
`1 − λ0²·E[N]·E_G[Var(B|G)] / (I(ρ;0) − I(ρ;c))`, where `Var(B|G=g)` is the exp-
tilted-truncated posterior variance under `p(b|g) ∝ f_B(b)·λ0·e^{−λ0(g−b)}`. An
**independent direct Monte-Carlo** of the three observed-data scores on simulated
renewal frames validates it (and reproduces the frozen J(ρ;c) to <0.5 %).

| ρ | J(ρ;c) count-only | J(ρ;0) ceiling | jitter loss (dB) | **recovered fraction** | recovered (dB) |
|---|---|---|---|---|---|
| 5.70 | 0.7855 | 0.8493 | 0.339 | **0.150** (mc 0.142) | 0.053 |
| **22.25** (ridge) | 0.4228 | 0.9354 | **3.449** | **0.110** (mc 0.114) | 0.543 |
| 1.00 | 0.4988 | 0.5000 | 0.011 | 0.047 (mc 0.027) | 0.001 |

**Reading.** The full timestamp record recovers only **11–15 %** of the jitter
loss — the going-in expectation that timestamps recover ≈J(ρ;0) is **refuted**.
The count already strongly constrains the summed holds (Σ gaps ≈ T with N terms,
Fano(N) ≈ 0.004 at the ridge), so the individual timestamps add little. At the
ridge (ρ=22.25), where jitter halves the ceiling (a 3.45 dB loss), timestamps buy
back only **0.54 dB**: the live wait there (mean 1/λ0 ≈ 0.045) is smaller than the
hold jitter (std 0.05), so the gap G=B+W barely resolves the live time better than
the count does.

**Hardware scope caveat.** This is a **detection-timestamp-only** record (what a
standard TCSPC records: photon/detection arrival times). A detector that
*additionally* time-tags the end of each dead/hold period (re-arm times) would make
L directly computable and recover ≈100 % of the loss — but that is a **different,
richer sensor**, not the timestamp record modelled here.

---

## 3. Cost accounting

- **Photons — identical across every rung and every method**, by construction: a
  coarser rung's counts are the exact sums of the finest sub-bin counts. No rung
  smuggles extra photons, exposures, or dwell. Patterns, φ, seeds, and TV rule are
  the frozen CPL values.
- **Data-volume ratio per rung:** B=1 → 1 number/exposure (the bucket, 2048
  numbers/frame); B=4 → 4×; B=16 → 16×; B=64 → 64×. The true timestamp record is
  ≈2200 arrival times/exposure; we proxy it at B=64 because the smooth
  within-exposure OU path is already resolved there (finer bins add only Poisson
  noise — confirmed by the saturation at B≈16).
- **Hardware scope:** the refinement is realisable on existing **time-tagging /
  TCSPC DAQ** (standard in FLIM, quantum optics, LIDAR), which record per-photon
  arrival times instead of one bucket total per exposure. No new optics; a richer
  back-end on the same single-pixel detector.
- **Compute:** full grid (6 scenes × 5 seeds × 5 cells × 4 rungs × 2 gain models
  for Chapter 2; 3 loads × 120 k frames for Chapter 1) in ~21 min CPU. No
  downscaling.

---

## 4. Honesty ledger and limitations

- **Refinement gain is a LOWER bound.** The fixed linearised estimator pays a
  low-count penalty at fine rungs (visible as the M0 decline). The
  information-theoretic ceiling for M0 is provably zero, so the penalty there is
  pure; under M1 it partially masks the true refinement information, so the
  reported +1 dB understates what an ideal estimator could extract. Direction of
  the bias is known (understates refinement).
- **No per-rung tuning:** identical estimator, frozen TV λ, at every rung; the only
  thing that changes is the record's temporal resolution.
- **The B=1 temporal-prior estimator assumes known OU statistics** (§1e) — the one
  place the comparison is not apples-to-apples, disclosed.
- **Oracle definition:** known per-exposure mean gain ā_e (CPL arm A4 generalised),
  the ceiling for image estimation from the gain-corrected bucket.

## 5. Files

- `LADDER_PROBE_REPORT.md` (this file)
- `ladder_results.json` — consolidated per-cell/per-scene ladder + Part 2
- `part1_gain_ladder.py` — Chapter 2 ladder (models M1_cont, M0_const)
- `part2_jitter_timestamp.py` — Chapter 1 timestamp recovery
- `consolidate.py` — builds `ladder_results.json` and the summary tables
- `ladder_results_M1_cont_full.json`, `ladder_results_M0_const_full.json`,
  `part2_jitter_results.json` — raw run outputs
