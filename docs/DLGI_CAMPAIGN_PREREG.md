# DLGI confirmatory campaign — frozen preregistration

**Program:** Dual-Ledger Ghost Imaging (DLGI) — one bucket record, a scene product
and a medium-dynamics product.
**Authority:** GPT R34 ruling (`docs/ROUND63_GPT_ROUND34_RULING_RAW.md`), building
on R33 (`docs/ROUND63_GPT_ROUND33_RULING_RAW.md`, the corrected reciprocity theorem).
**Feasibility evidence:** `results/round63_next/DUAL_LEDGER_PROBE/` (scored 6/7; the
one failed item was finite-sample `t_c` interval coverage, which this campaign
repairs *once* with a preregistered simulation-calibrated construction).
**Status:** protocol FROZEN. This document is committed before the confirmatory bank
is opened. No tuning, regridding, estimator switching, or interval repair is
permitted after the calibration table is committed.

This file is the normative protocol. The scaffold implementing it lives in
`code/round63_dlgi/`; sealed data and frozen artifacts in
`results/round63_dlgi_campaign/`.

---

## 1. Three sealed, disjoint scene banks (R34 §2.1)

Each bank is generated **blind** (no reconstruction, arm, or bucket metric touches
the images; only image-domain statistics). Every image is a pure deterministic
function of its integer seed (numpy `default_rng` only); regeneration is
byte-reproducible (`dlgi_scene_bank.py` self-checks determinism + all-distinct
across banks). **No scene, and no seed, crosses banks.** Seed bases are disjoint
from every prior repo cohort (630900 / 631000 / 631900 / 632000 / 632900 / 633000 /
650000, the last being the 6 DEV bridge scenes used by the feasibility probe).

| bank | role (and ONLY that role) | seed base | scenes |
|---|---|---|---|
| **calibration** | construct the critical-value table `c_0.95(eta)` | 710000 + 100·stratum + i | 12 (6 strata × 2) |
| **confirmatory** | final endpoints only (fresh scenes/seeds) | 700000 + 100·stratum + i | 12 (6 strata × 2) |
| **edge** | display the fast/slow-drift + weak-CV failure map | 720000 + 100·stratum + i | 6 (6 strata × 1) |

**Strata (6, canonical order):** `contour`, `texture`, `twopop` (two-population),
`smooth`, `binsparse` (binary-sparse), `natural`. `contour`/`texture`/`twopop`
reuse the frozen `detail24` and `gen_bridge_scenes` generators **verbatim** (import
only); `smooth` (low-pass Gaussian field), `binsparse` (sparse bright features), and
`natural` (1/f^β power-law field — a synthetic natural-image-statistics proxy,
**not** a photographic crop) are new generators archived in `dlgi_scene_bank.py`.
Cohort mean intensity is held in the bridge-scene band (~0.16–0.56) so per-scene
detected-count levels match the probe's regime under the frozen forward PHI.

Each bank writes `MANIFEST.json` with, per scene: `scene_id`, `stratum`, `seed`,
`gen_params`, `x_sha256` (content hash of `x.tobytes()`), `sha256_npz` (file bytes),
and image statistics. The manifests are the freeze fingerprint.

---

## 2. Frozen primary grid (R34 §2.2)

Δ = 1 exposure, so `t_c/Δ` equals `t_c` in exposure units (probe convention).

```
t_c/Δ       in {16, 32, 64}
CV          in {0.05, 0.15, 0.40}
photon scale in {0.5, 1, 2} × the frozen nominal 2200 counts/exposure
```

**Primary claim spans 27 cells.** Difficult cells are retained. The three cells that
FAILED intervals in the probe are **mandatory named stress cells** at nominal SNR:
`(16, 0.05)`, `(16, 0.40)`, `(64, 0.05)`.

**Edge cells (define the published failure map; cannot rescue or kill the interior):**
fast-drift `t_c/Δ = 2` and slow-drift `t_c/Δ = 128` (slower than 64), each across
the CV sweep at nominal SNR.

The scalar-per-frame gain model, `E[a]=1` gauge, exposure count (2048 complementary
exposures), object dimension (32×32, NPIX=1024), full 1024-row sequency-Hadamard
pattern multiset, and the Poisson noise law are frozen exactly before generation.
The point **estimator is bit-identical to the probe** — imported verbatim from
`results/round63_next/DUAL_LEDGER_PROBE/code/dl_common.py`; the campaign changes only
the uncertainty statement (§4). The photon-SNR axis multiplies the frozen PHI; at
photon scale 1 the forward is bit-identical to the probe (self-tested in
`dlgi_common.py`), and the `log(photon_scale)` offset is absorbed by the estimator's
gauge demean, leaving the medium estimate unchanged in output.

**Replicate counts:** ≥ 5000 full-pipeline calibration replicates per primary cell;
≥ 1000 independent confirmatory records per primary cell, balanced across the fresh
≥12-scene cohort.

---

## 3. Arms and schedules (R34 §2.3)

1. **Pilot-free DLGI** — the ordinary bucket record supplies both the scene and the
   medium estimates (`joint_dual_ledger`, frozen).
2. **Pilot-interleaved baseline** — the frozen strongest 5% pilot design; identical
   total exposures/photons/duration; pilot rows replace scene rows and their scene
   cost remains charged. Mandatory strongest baseline.
3. **Oracle monitor** — the gain path observed at every exposure; nondeployable
   information ceiling.
4. **Plain-linear camera** — reconstruction from the byte-identical pilot-free
   record with no medium correction (scene-noninferiority endpoint).

**Pilot-free schedules:** exactly `{paired, random, ordered}`, with the same pattern
multiset, photons, exposure count, and duration — only temporal order changes. The
raw record is preserved; no differencing step discards the common-mode medium
information. The **paired** schedule is the primary acquisition (best for both
ledgers in the probe, as the reciprocity theorem predicts).

*Naming reconciliation (flagged):* the probe implements schedules `{paired, split,
random}`; R34's `ordered` is the probe's deterministic non-paired `split` order
(complementary members maximally separated in time). The frozen `make_schedule` is
reused verbatim through this fixed rename.

---

## 4. Analysis plan — the authorized interval remedy (R34 §1)

Parameterize `eta = (log t_c, log CV)`. For every candidate `eta_0` on the frozen
grid and photon-SNR stratum, compute the frozen full-pipeline profile statistic

```
W(Y; eta_0) = 2 { ell_joint(eta_hat, x_hat) − ell_joint(eta_0, x_hat_eta0) }.
```

The exact-likelihood approximation is the **frozen log-Whittle/profile construction**
from the probe (`dl_bar4_final.whittle_nll_grid` + `_periodogram`, imported verbatim):
`ell_joint(eta, x_hat) = −whittle_nll(eta)` on the time-ordered log-gain residual
series `z(x_hat)` with the white floor `sv2` profiled; `x_hat_eta0 ≈ x_hat` (the
joint reconstruction), so `W` over the `(log t_c, log CV)` grid is the joint 2-DOF
profile-LR surface `W_grid = 2 (nll_grid − min nll_grid)`.

Every calibration replicate **regenerates the gain path and Poisson noise and reruns
the entire scene+medium pipeline** (`joint_dual_ledger`). Holding the reconstructed
scene fixed is forbidden.

Before confirmatory data are generated:

1. simulate the calibration bank at every declared `(t_c, CV, SNR)` cell (≥5000
   full-pipeline replicates);
2. `c_0.95(eta_0)` = **conservative Monte-Carlo order statistic** of
   `W(Y_rep; eta_0 = true)` — the `k`-th order statistic with
   `k = ⌈(n+1)·0.95⌉` (errs toward over-coverage);
3. freeze a **monotone upper-envelope interpolation** in `(log t_c, log CV, log SNR)`:
   the enclosing-box maximum of the knot critical values, with hold-flat
   extrapolation outside the calibration box (conservative between knots);
4. commit the critical-value table, interpolation rule, code hash, and seeds
   (`CALIBRATION_TABLE.json`) **before** opening the confirmatory bank;
5. form the confidence region by Neyman inversion
   `C_0.95(Y) = { eta_0 : W(Y; eta_0) ≤ c_0.95(eta_0) }`, and report the projections
   as the `t_c` and CV intervals.

These are **simulation-calibrated, model-conditional intervals under the validated
scalar-gain/OU model** — not distribution-free guarantees.

**Preregistered `(A,B,K)` predictions (R34 §2.4)** are frozen in
`predictions/ABK_PREDICTIONS.json` before any confirmatory data: `A = I_xx`,
`B = I_thetatheta`, `C = I_xtheta`, `K = A^{-1/2} C B^{-1/2}`,
`J_x = A^{1/2}(I−KK^T)A^{1/2}`, `J_theta = B^{1/2}(I−K^TK)B^{1/2}`, with the frozen
predicted ordering (full-record Fisher info schedule-invariant; paired best for both
ledgers at the estimator level; no forced scene-vs-medium trade). The empirical
ordering may not be selected after seeing results.

---

## 5. Binding kill bars (R34 §3, verbatim)

All bars must pass. There is no averaging away a failed parameter cell. Failure of
any confirmatory bar ends the flagship claim; no additional repair round is
authorized.

### Bar C1 — calibrated `t_c` coverage
For every primary `(t_c, CV, SNR)` cell, the empirical coverage of the projected
nominal 95% `t_c` interval must lie in **[0.92, 0.98]**. Report exact binomial
intervals and multiplicity-adjusted diagnostics. The three previously failed cells —
`(16, 0.05)`, `(16, 0.40)`, `(64, 0.05)` at nominal SNR — are mandatory named stress
cells.
**Kill:** any primary cell below 0.92, or any unreported/nonfinite interval rate
above 5%.

### Bar C2 — interval informativeness
At every primary cell: median log-scale `t_c` interval width ≤ 1.5× pilot baseline;
median CV interval width ≤ 1.5× pilot; intervals bounded and connected in ≥95% of
records. Also report oracle-relative width.
**Kill:** coverage passes only by producing intervals that fail these
width/boundedness bars.

### Bar C3 — point precision
Pilot-free RMSE for both `t_c` and CV ≤ 1.5× the strongest pilot arm at every primary
cell. Report the oracle ratio and signed bias separately.
**Kill:** any cell exceeds 1.5×, or the previously observed advantage disappears
because the calibration procedure changes the estimator rather than only its
uncertainty statement.

### Bar C4 — scene noninferiority
Against the plain-linear reconstruction from the byte-identical record:
`PSNR loss ≤ 0.2 dB` and `task/Fisher loss ≤ 5%` at every primary cell, with
cohort-level confidence intervals. Pilot-free improvement is welcome but not
required.
**Kill:** either noninferiority margin is crossed.

### Bar C5 — physical model and gauge
Repeat held-out innovation, residual-whiteness, scalar-gain adequacy, and scale-gauge
tests on the fresh scenes and all SNR strata.
**Kill:** systematic pattern-dependent residuals, gauge sensitivity, or a spatially
varying medium operator large enough to invalidate the scalar ledger.

### Bar C6 — reciprocity and schedule prediction
Numerically verify `det(J_x)/det(A) = det(J_theta)/det(B)` and the shared
canonical-confusion spectrum to numerical tolerance. The preregistered `(A,B,K)`
schedule ordering must predict the measured ordering for both ledgers, not merely
correlate post hoc.
**Kill:** theorem implementation fails, or empirical schedule behavior contradicts
the frozen prediction on more than a predeclared 10% of primary cells.

### Bar C7 — identical-ledger discipline and edge honesty
Byte-audit photons, exposures, duration, and pattern multiset. Publish the fast- and
slow-drift edge failures and the weak-CV/noise boundary.
**Kill:** hidden pilot/reference information, unequal resource accounting, or removal
of a failed edge after unblinding.

---

## 6. Claim ruling (R34 §4)

Until the campaign passes, the authorized statement is:

> DLGI extracts a competitive point estimate of medium correlation time and
> fluctuation strength from the same bucket record used for scene reconstruction,
> with no added acquisition and no detected scene penalty; model-calibrated `t_c`
> uncertainty remains under confirmation.

Do **not** use "two certified products" yet. If all confirmatory bars pass, the
sharpened R33 headline is restored:

> One bucket stream, two **model-certified** products: pilot-free ghost imaging
> reconstructs the scene and estimates the medium's correlation time and fluctuation
> strength with calibrated intervals and no additional measurements, while matching
> or beating a pilot monitor under the same acquisition ledger.

The words **model-certified** and the validated scalar-gain operating domain are
mandatory. DLGI is not a universal replacement for DCS, a reference arm, or a monitor
when the medium changes the spatial transfer operator. Manuscript placement: DLGI
conditionally leads the METHOD act of the one-flagship manuscript; passing this
campaign licenses that lead, failing it demotes DLGI to a shorter theorem/method note.

---

## 7. Commit-before-confirmatory discipline

1. Freeze + commit the three sealed bank manifests (§1) and this document.
2. Freeze + commit `predictions/ABK_PREDICTIONS.json` (§4).
3. Run the calibration bank → **commit `CALIBRATION_TABLE.json`** (the critical-value
   table, interpolation rule, code hash, seeds).
4. **Only then** generate the confirmatory bank records and evaluate bars C1–C7.
5. No tuning/regridding/estimator-switch/interval-repair after step 3. No third
   post-hoc interval repair is authorized.

The calibration launch (step 3) is a **separate GO decision**; tonight's scaffold
stops after the machinery + smoke.

---

## 8. Smoke verification (machinery unit-test)

`dlgi_neyman.py` smoke on cell `tc64_cv40` (nominal SNR; a cell where the probe's
profile likelihood already worked), held-out split (180 calibration replicates set
`c_0.95`, 180 disjoint confirmatory-bank replicates measure coverage):

- MC-calibrated critical value `c_0.95 = 39.9` vs asymptotic `χ²₂ = 5.99`
  (ratio 6.66× — the asymptotics badly under-shoot);
- **held-out pointwise joint coverage = 0.956** (target ~0.95, inside [0.92, 0.98]);
- asymptotic-`χ²₂` joint coverage = 0.628 (confirms the R34 diagnosis: the naive
  threshold under-covers; simulation calibration is necessary and sufficient here).

This verifies the Neyman machinery end-to-end. Full projected-interval coverage
(bar C1) requires the complete 27-knot calibration table (the separate GO).

---

## 9. File inventory

**Code** (`code/round63_dlgi/`): `dlgi_scene_bank.py` (banks + generators),
`dlgi_common.py` (config, seed hygiene, forward, frozen re-exports),
`dlgi_neyman.py` (W-statistic, `c_0.95`, `CritSurface`, inversion, smoke),
`dlgi_arms.py` (four arms, three schedules), `dlgi_abk.py` (frozen `(A,B,K)`
predictions), `dlgi_calibration_launcher.py` (cost/shard plan, cell/shard runners,
merge).

**Data/artifacts** (`results/round63_dlgi_campaign/`): `scene_banks/{calibration,
confirmatory,edge}/` (npz + png + `MANIFEST.json`), `predictions/ABK_PREDICTIONS.json`,
`smoke/smoke_tc64_cv40.json`, `calibration_plan/COST_AND_SHARD_PLAN.json`.
Produced by the separate GO: `CALIBRATION_TABLE.json`.

**Cost:** 27 × 5000 = 135,000 full-pipeline replicates at 1.26 s each = 47.2
core-hours → ~1.5 h on 32 cores (or ~8 h on the memory-safe worker count of this
16.9 GB box; each worker re-imports numpy/scipy at ~0.4 GB, so cap workers with
`--nproc`). Well under 12 h → **local multiprocess recommended**; a 5-route Colab
split (~5 h/route) is specified as an alternative.
