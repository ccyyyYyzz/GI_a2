# E6 -- Does measurement-design redundancy kill the bilinear collusion?

**Date:** 2026-07-23/24 · **Status:** DEV-GRADE, **STOPPED EARLY by operator** (compute
discipline — the two mini-scale predictions had already landed and E5 already resurrected the
direction via the covariance route). · **Ran on:** Google Colab **pro2 (Pro+), 2x NVIDIA L4**
(torch 2.11.0+cu128, numpy 2.0.2, float64), per COLAB_USAGE_GUIDE.md. VMs **released**
(assignments = NONE) and keep-alives killed at stop.

## Headline verdict

> **Of the arms that completed (A0 control, A2 chronology-randomization, A4 chronology+temporal-law),
> EVERY ONE DRIFTS on the true-seed test (null-NMSE 0.69-0.80 clean, ~0.79-0.90 @1e5). The
> collusion is NOT killed by chronology randomization (A2) or by temporal-law enforcement (A4)
> at d=48.** Blind cold-start fails on every arm (>=1.0, worse than fixed averaging); the
> covariance route fails in absolute terms at this T=128 cell (T_eff~4) but its RELATIVE ordering
> matches the coordinator's prediction. **The spatial route-redundancy arms A1 / A3 (Fourier
> overlap) were NOT RUN (stopped first) — the one design the pilot smoke hinted might resist the
> drift remains OPEN at full scale.**

The GI_a1 design principle (randomize chronology to break the gain-object mimicry) **does not
visibly generalise from d=1 to d=48 on the pathwise estimator** here: making the object's
temporal carrier white (certificate flatness 0.00 -> 0.56, verified) did not move the pathwise
medium-re-estimation minimum back onto truth.

## Cell (fixed across all arms)

N=1024 (32x32), M=96 patterns, d=48 DCT smooth medium basis, sigma_f=0.15 mean-normalized
**lognormal** medium (E[w]=1/pixel), OU tau=16, **T=128 medium blocks x S=96 slots/block =
12288 exposures fixed across arms**, clean + 1e5 photons/bucket, 5 seeds. Metric = null-NMSE
after the oracle scale gauge (fixed-averaging baseline = 1.0). Conventions copied from
fog_common.py / fog_tracker.py. Scene null fraction 0.68 (binary) / 0.72 (Fourier).

## Verdict table (arm x estimator x endpoint)

null-NMSE. TRUE-SEED = full-Z refinement seeded from the true medium path, medium re-estimated
(the collusion detector): STAY(<=0.10)=collusion killed, DRIFT(~0.7)=collusion present
(calibrated so A0 reproduces the E5 Sec.11c drift). BLIND = data-only cold start, 5-seed median;
agree = median within-seed std across data-only starts. COV = covariance-domain (E5a) blind.

| arm | photons | oracle | **PW true-seed** | verdict | PW blind (agree) | COV blind (agree) | cert SF | n_seed |
|---|---|---|---|---|---|---|---|---|
| A0 (ctrl, rand+cart) | clean | 0.001 | **0.73-0.75** | DRIFT | 0.805 | 0.833 | 0.000 | 1 + calib |
| A0 | 1e5 | 0.027 | **0.86-0.90** | DRIFT | -- | -- | 0.000 | calib |
| A2 (rand+**slot**) | clean | 0.001 | **0.694** | DRIFT | 1.073 (0.006) | 0.894 (0.018) | 0.560 | 5 |
| A2 | 1e5 | 0.032 | **0.821** | DRIFT | 1.074 | 0.885 | 0.553 | 1* |
| A4 (rand+slot+multiτ+law) | clean | 0.001 | **0.802** | DRIFT | 1.074 (0.006) | 0.883 (0.015) | 0.560 | 5 |
| A4 | 1e5 | 0.026 | **0.792** | DRIFT | 1.070 | 0.889 | 0.553 | 1* |
| **A1** (Fourier+cart) | -- | -- | **NOT_RUN_STOPPED** | -- | -- | -- | -- | 0 |
| **A3** (Fourier+slot) | -- | -- | **NOT_RUN_STOPPED** | -- | -- | -- | -- | 0 |

`*` 1e5 rows for A2/A4 are single-seed partials fetched from the live driver logs at the stop
(the 1e5 cells were mid-sweep). A0 rows come from the calibration run (calib_trueseed.py, same
ts config outer=30/zsteps=60) plus the first Colab launch's seed-0 (A0 is random-binary, so
unaffected by the Fourier-bank bug that forced the mid-run relaunch). A0 was re-queued last on
its session (priority A2>A3>A4>A0>A1) and the stop landed before its 5-seed grid re-ran.

### Per-arm one-liners
- **A0 (control):** DRIFTS to 0.73-0.90 — reproduces the E5 Sec.11c collusion. Baseline OK.
- **A2 (chronology randomization, the GI_a1 port):** DRIFTS to 0.69 clean / 0.82 @1e5. The
  carrier is fully whitened (SF 0.56 vs A0's 0.00) yet the pathwise medium re-estimation still
  collapses to a wrong-null solution. **Chronology randomization did not kill the d=48 collusion.**
  Blind is actually WORSE than A0 (1.07 vs 0.81): the slot schedule's per-block partial pattern
  coverage removes the full-rank per-block constraint A0's Cartesian schedule provides.
- **A4 (A2 + multi-timescale medium + temporal-law penalty):** DRIFTS to 0.80/0.79. Penalising
  OU-timescale mixing did not prevent the drift — the wrong-null solution respects the declared
  law well enough to survive the penalty.
- **A1, A3 (Fourier overlap):** NOT RUN. The pilot smoke (under-refined) hinted the Fourier arms
  resisted drift (TS ~0.06-0.09), but that was not confirmed at full refinement, and a confound
  (Fourier patterns may make the medium under-observable, so "stay" could be "medium can't move"
  rather than "collusion killed") was never resolved. **Open question.**

## Endpoint detail

**1. True-seed stability (collusion detector).** Every completed arm DRIFTS. Calibrated on A0:
at outer=30/zsteps=60/near-pure-likelihood the A0 refinement leaves truth for ~0.73 (clean)
climbing toward ~1.0 with more iterations — a genuine wrong-null minimum, matching Sec.11c.
A2 and A4 drift to the same regime at matched budget.

**2. Blind cold-start + multi-start agreement.** All arms fail (median 0.81-1.07). The
within-seed multi-start agreement is TINY (std 0.006 pathwise; 0.015-0.018 cov): the data-only
starts **agree on the same wrong solution**. This is the R38 Sec.1.5 signature — the blind
likelihood landscape places truth away from its optimum, so the failure is identifiability, not
optimizer luck. (Cov agreement is likewise tight -> the covariance estimator is injective/
collusion-free as in E5, it just needs far more temporal samples than T=128 supplies.)

**3. Whitening certificate vs outcome (slot arms).** The Fourier/Walsh spectral-flatness
certificate cleanly separates schedules (A0 Cartesian SF=1.9e-34 / periodic carrier;
A2/A4 slot SF~0.55 / white carrier). Across the 5 seeds the certificate **negatively correlates
with blind null-NMSE** (A2 corr -0.97, A4 corr -0.96): whiter realized carrier -> lower blind
NMSE. **Caveat:** the dynamic range is tiny (SF 0.553-0.569, blind 1.04-1.09), so the
correlation is a direction signal inside the failing regime, not evidence of recovery. Per-seed
data in `E6_whitening_certificate.json`.

## Prediction-vs-outcome scorecard (coordinator mini-scale, 16x16)

| prediction | full-scale (32x32) test | outcome |
|---|---|---|
| **P1**: Fourier-lattice patterns FAIL on the covariance estimator (~1.12) | A1/A3 x COV | **UNSCORED — A1/A3 NOT_RUN.** (Pilot smoke on the bugged rank-63 Fourier bank gave cov ~0.95, weakly consistent but not a clean test.) |
| **P2**: on COV, ordered/full-coverage BEATS randomized slot | A0 x COV (Cartesian=full pair-coverage=ordered comparator) vs A2 x COV (randomized slot=masked ~63%) | **DIRECTIONALLY CONFIRMED (weak):** A0 cov 0.833 < A2 cov 0.894 (both fail absolutely; A0 single-seed vs A2 5-seed). Ordered/full-coverage is better, as predicted. |
| Thesis: chronology randomization helps PATHWISE, neutral/harmful on COV | A2 pathwise vs A0 pathwise; A2 cov vs A0 cov | **PARTIAL:** on COV, A2 (randomized) is worse than A0 (ordered) — matches. On PATHWISE, A2 did NOT beat A0 (both drift ~0.7); at d=48 randomization did not deliver the d=1 GI_a1 benefit. |

## Compute spent vs saved

- **Spent (Colab):** 2x L4 on pro2, ~3 driver relaunches (one forced by the Fourier rank-63
  bug fix, one by the A4-rebase/priority-reorder). Completed on-VM: A2 {clean x5, 1e5 x1},
  A4 {clean x5, 1e5 x1}, plus first-launch A0 seed0 — roughly ~35 cell-solves x ~40 s ~= **~23
  GPU-minutes** across the two L4s. **Local (RTX 4060):** calibration + 2 smokes + 1 validation
  ~= **~15 min**.
- **Saved by the stop:** A0 {9 cells}, A1 {10}, A3 {10}, A2/A4 {8 remaining 1e5 seeds} ~= **~37
  cells (~25 GPU-min)** not run.

## Honest limitations / what stays open
- **A1/A3 (spatial route-redundancy) untested at full scale** — the single design most likely
  (per pilot) to behave differently. The confound (Fourier -> under-observable medium) was not
  resolved.
- **Block-constant medium** for slot arms: medium held constant within each 96-slot block, so
  A2/A4 keep A0's block-level temporal correlation (tau=16) and do NOT add independent medium
  realizations. This is the fair setting for the pathwise-collusion question but cannot test
  whether continuous slot-level drift would improve the covariance-route sample complexity.
- Reduced-Q ALS + unweighted-LS x-solve are the E4/E5 dev-grade approximations. T=128
  (T_eff~4) is far below the covariance route's ~10x sample requirement (E5), so absolute cov
  failure here is expected and not informative about the covariance direction (E5 already
  settled that positively).

## Deliverables in this directory
- `E6_REPORT.md` (this file), `E6_REPORT_TABLE.md` (auto verdict table).
- `E6_results.json` (merged A2+A4 5-seed cells), `E6_s1.json`/`E6_s2.json` (raw per-session),
  `E6_whitening_certificate.json` (per-seed certificate + outcome).
- `e6_s1.log` / `e6_s2.log` (full driver logs incl. the 1e5 partials).
- Code: `fog_e6.py` (self-contained model+estimators+certificate), `run_e6.py` (driver),
  `make_report.py`. Diagnostics: `scripts/design_diag.py` (bank rank / schedule redundancy),
  `scripts/calib_trueseed.py` (A0 drift calibration). VM ops: `vm/*.sh` + `vm/launch_s*.py`
  (create / deploy / launch / rebind / watchdog / release, all pro2 + --auth oauth2 + timeout).
- No git commit (per instruction; coordinator reviews).
