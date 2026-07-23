# DL-as-tool for DLGI: the analytic gauge fix, the learned temporal prior, and the certified ceiling

**EXPLORATORY / DEV-ONLY. Not preregistered. No commits. Local GPU only** (RTX 4060 Laptop, `pytorch` env py3.9 / torch2.1.0+cu121). All work on DEV assets: the 6 `data/r63_bridge_scenes/` bridge scenes + freshly-generated procedural training scenes with far-offset, disjoint seeds. The sealed confirmatory/edge banks under `results/round63_dlgi_campaign/confirmatory/` were never read or touched. Colab and the running fleet untouched. The frozen estimator `DUAL_LEDGER_PROBE/code/dl_common.py` was imported **read-only**; every variant here lives in an exploratory copy (`code/dl_tool_common.py`).

Program discipline honored throughout: DL is a **tool racing model-based / certified baselines**, never the headline; every learned component gets a held-out test and an explicit "where it loses."

---

## 0. Conclusions first

**E1 — the +3.5 dB corner deficit is a gauge bug, and physics (not learning) fixes most of it.**
Yesterday's fast-strong-drift corner (t_c=16, CV=0.40) deficit of DLGI vs in-family SCGI was diagnosed as a gain-*level* error. It is: at the corner the frozen DLGI under-estimates the gain level (`mean(â)/mean(a_true) = 0.757`), over-brightens the scene by 32%, while keeping gain-path correlation 0.996. **A one-line analytic fix — renormalize the estimated gain to the known physical gauge E[a]=1 — recovers +2.19 dB (D_geom +2.40 dB), closing 63% (69%) of the 3.48 dB deficit to SCGI with ZERO learning, and it does not touch the certified medium (t_c, CV) product.** The frozen estimator *already* carried the Jensen `+½V` term; it is negligible at 2200 counts (turning it off changes nothing), so the missing piece was purely the gauge. A tiny learned gain-path corrector buys at most **+0.96 dB more, only at the corner** (to 11.94, within 0.34 dB of SCGI), and **regresses the other three cells by 0.2–0.8 dB** — the learned layer does not robustly beat the free physics fix.

**E2 — a learned temporal prior is a genuine but targeted domain extension.**
The paper's validated domain is scalar OU. On three non-OU gain-path classes, the OU model's statistical adequacy provably breaks (regime-switching: innovation whiteness fails, 47% pass, p≈7e-3; heavy-tailed: whiteness holds but innovation excess-kurtosis +9.1; quasi-periodic: single-exponential t_c biased 37%). A small trained sequence model (TCN smoother + causal-predictor "learned transition density", each ≈105k params, both <1M) **earns +2.46 dB scene PSNR and restores innovation whiteness (44%→89%) on the regime-switching class — where the OU *autocovariance* is structurally wrong — while being no-harm on true OU (−0.07 dB).** It does **not** universally win: neutral on heavy-tail (−0.08 dB; the mis-specification is distributional, not autocovariance), and it **loses −0.92 dB on quasi-periodic**. The biased correlation-time readout for coherent drift is better fixed by a *learned autocovariance* (empirical, model-free: quasi-periodic tau rel-err 0.37→0.06) than by the trained smoother. Net: the domain extension matters most for the **measurement/certification** claim (you cannot certify OU on these media), and for imaging only where the autocovariance *structure* is wrong.

**E3 — the ceiling referees.** One figure (`figs/e3_certified_ceiling.png`) places every method on a performance-vs-certified-ceiling axis for the corner cell: floor → frozen DLGI → +physics fix → +learned tools all race the oracle known-gain ceiling, while the medium channel's Cramér–Rao floor certifies the measurement the learned correctors cannot.

**One-line takeaway:** model-based fixes first — a free gauge correction did 63% of the work the learned corrector was credited with; the learned layer earns real dB only where physics genuinely cannot write the model down (regime-switching drift), and even there it must be watched for the cells where it loses.

---

## 1. E1 — the Jensen/scale-bias fix, analytic first

### 1(a) Diagnosis: the deficit is a gain LEVEL error, and physics can name it
Decomposing the frozen DLGI reconstruction PSNR into a best-scalar-rescale (shape-limited) part and the residual level loss, over the 4 test cells (6 scenes × 3 seeds; `json/e1_diagnose.json`):

| cell | PSNR_D | +best-scale | level loss | gain ratio â/a | brightness x̂/x | gain corr |
|---|---|---|---|---|---|---|
| t_c=16, CV=0.15 | 15.40 | 15.56 | +0.16 | 0.963 | 1.040 | 0.988 |
| **t_c=16, CV=0.40** | **8.80** | **11.01** | **+2.22** | **0.757** | **1.323** | **0.996** |
| t_c=64, CV=0.15 | 17.23 | 17.41 | +0.18 | 1.001 | 0.999 | 0.994 |
| t_c=64, CV=0.40 | 14.02 | 14.33 | +0.31 | 0.947 | 1.060 | 0.998 |

At the fast-strong corner the gain is under-estimated by ~24% (scene 32% too bright) with correlation 0.996 — a pure multiplicative *level* error that grows sharply with CV, exactly the "Jensen/gauge term at high CV" signature. An oracle single-scalar rescale would recover +2.22 dB.

### 1(b) The analytic fix and what it isolates
In `code/dl_tool_common.py` I built a switchable copy of `medium_estimate`/`joint_dual_ledger` with three toggles: **gauge** (frozen counts-weighted vs uniform), **jensen** (the `+½V_sm` posterior-mean term), **renorm** (`none` / `arith` = impose E[a]=1 / `geom`). Head-to-head on the 4 cells (`json/e1_fix_headtohead.json`, `figs/e1_headtohead.png`):

| cell | floor | D0 frozen | **D_arith** | D_geom | D_noJ_arith | D_uni_arith | SCGI (L_mix) | oracle |
|---|---|---|---|---|---|---|---|---|
| t_c=16, CV=0.15 | 13.84 | 15.40 | 15.46 | 15.44 | 15.46 | 15.46 | 15.09 | 18.78 |
| **t_c=16, CV=0.40** | 6.16 | **8.80** | **10.98** | 11.20 | 10.98 | 10.98 | 12.28 | 18.39 |
| t_c=64, CV=0.15 | 16.46 | 17.23 | 17.19 | 17.13 | 17.19 | 17.19 | 15.73 | 18.66 |
| t_c=64, CV=0.40 | 11.30 | 14.02 | 14.23 | 14.09 | 14.23 | 14.23 | 12.67 | 18.48 |

**The key numbers.** Corner: D0 8.80 → **D_arith 10.98 (+2.19 dB), D_geom 11.20 (+2.40 dB)** — 63% / 69% of the 3.48 dB deficit to SCGI, essentially exhausting the +2.22 dB scale headroom, **without any learning**. Three diagnostics of *what* the fix is:
- `D_noJ_arith == D_arith`: the Jensen `+½V` term is **negligible** (obs noise ≈1/2200 → posterior variance tiny). The frozen code already had it; it was never the problem.
- `D_uni_arith == D_arith`: the counts-vs-uniform gauge choice is irrelevant once you renormalize. **The entire fix is the E[a]=1 gauge renormalization.**
- The medium (t_c, CV) rel-errors are **identical** for D0 and D_arith (renorm is a pure brightness gauge) — the certified science product is untouched. Off-corner cells move by ≤0.2 dB.

Adopt **D_arith** as the fix (arithmetic mean = 1 is the exact statement of the known E[a]=1 constraint, most robust across cells; D_geom is marginally better only at the extreme corner).

### 1(c) The learned residual corrector — marginal and non-robust
Because a material residual remains (1.3 dB to SCGI, 7.4 dB to oracle at the corner, all *shape*), I trained a tiny 1-D dilated-CNN (47.8k params) that refines the DLGI-smoothed **log-gain path** (the physics-adjacent quantity — it earns the local smoothing the fixed OU misses at fast-strong drift), on procedural scenes, certified on held-out bridge scenes (`json/e1c_corrector.json`, `checkpoints/e1c_corrector.pt`):

- In-distribution it cuts the gain-path log-MSE by **92%** (5.9e-2 → 4.7e-3 on the procedural val set).
- On the held-out **bridge** scenes it helps **only the fast-strong corner: +0.96 dB → 11.94** (within 0.34 dB of SCGI, and beating SCGI's need for a corruption-family assumption), but **regresses the other three cells by −0.37 / −0.79 / −0.20 dB** and increases gain-path RMSE in every cell.

**Honest verdict:** the physics fix earns the robust +2.2 dB; the learned layer adds marginal, corner-specific, non-robustly-transferring extra. The learned layer only earns the part physics cannot write down — and here that part is small and would need a CV-gated deployment to avoid the regressions.

---

## 2. E2 — learned temporal prior = domain extension

Domain: the paper's validated scalar-OU gain. I built three non-OU gain-path classes (`code/e2_common.py`), all mean(a)=1: **regime** (Markov switch between a fast/strong and slow/weak OU → 2-exponential autocovariance), **heavytail** (OU recursion with Student-t innovations → jump-like increments), **quasiper** (OU + a random-phase/period sinusoid → oscillating autocovariance), plus **ou** as control.

### 2(b) The frozen OU model misfits the non-OU classes
True-scene log-gain observation (isolates the temporal-gain model from reconstruction); 6 bridge scenes × 6 seeds (`json/e2_ou_misfit.json`):

| class | whiteness pass | Ljung-Box p | innov. kurtosis | CV rel-err | single-exp τ rel-err | learned-autocov τ rel-err |
|---|---|---|---|---|---|---|
| ou (control) | **100%** | 0.46 | −0.02 | 0.05 | 0.19 | 0.18 |
| regime | **47%** | 0.007 | +2.58 | 0.03 | 0.26 | 0.43 |
| heavytail | 97% | 0.44 | **+9.13** | 0.06 | 0.21 | 0.18 |
| quasiper | 75% | 0.34 | −0.06 | 0.03 | **0.37** | **0.06** |

The OU model's adequacy breaks in class-specific ways: **regime** fails innovation whiteness (the autocovariance is wrong); **heavytail** *passes* whiteness (Student-t innovations share the AR(1) autocovariance — correctly not flagged) but shows a strong distributional tell (excess kurtosis +9); **quasiper** biases the single-exponential correlation time by 37%. CV (fluctuation strength) is recovered robustly everywhere (prior-independent, noise-deconvolved). A **learned/empirical autocovariance** τ readout fixes the quasi-periodic bias (0.37→0.06) but is not better for regime/heavytail.

### 2(c–e) The learned temporal prior, honest recon loop
Two trained sequence models on the mixed non-OU family (procedural scenes, disjoint seeds): a **SmootherTCN** (acausal dilated TCN, 104.6k params → smoothed log-gain path) and a **PredictorTCN** (causal dilated TCN with Gaussian head, 104.6k params — a learned one-step transition density whose standardized residuals are the learned analog of the Kalman innovations). Checkpoints in `checkpoints/e2_*.pt`.

Scene/gain channel evaluated in the **honest full reconstruction loop** (A2-scene log-gain observation — no true-scene leak; the earlier true-scene isolation inflated PSNR *above* the oracle and is discarded). Model adequacy from the causal predictor. Per class, 6 bridge scenes × 6 seeds (`json/e2_recon_psnr.json`, `json/e2_learned_prior.json`, `figs/e2_learned_prior.png`):

| class | PSNR OU-model | PSNR learned | oracle | Δ(learned−OU) | gain-RMSE OU→learned | whiteness OU→learned |
|---|---|---|---|---|---|---|
| ou (control) | 12.80 | 12.73 | 18.65 | **−0.07** (no-harm) | 0.056→0.058 | 97%→89% |
| **regime** | 9.98 | **12.44** | 18.58 | **+2.46** | 0.037→0.035 | **44%→89%** |
| heavytail | 13.82 | 13.74 | 18.77 | −0.08 (neutral) | 0.027→0.031 | 100%→75% |
| quasiper | 15.48 | 14.55 | 18.77 | **−0.92** (loses) | 0.022→0.026 | 58%→42% |

**Readout.** The learned temporal prior is a **targeted** tool: it earns **+2.46 dB and restores innovation whiteness (44%→89%)** on the regime-switching class — the one class where the OU *autocovariance structure* is genuinely wrong — and is **no-harm on true OU (−0.07 dB)**. It is neutral on heavy-tail (whose mis-specification is distributional, not autocovariance, so a Gaussian smoother — learned or not — cannot exploit it) and it **loses on quasi-periodic (−0.92 dB, and whiteness 58%→42%)**: the trained model failed to capture the random-phase sinusoid better than the OU smoother's own tracking. The medium correlation-time bias of the coherent-drift class is instead fixed by the learned/empirical autocovariance (§2b, 0.37→0.06).

**Boundary (honest):** at the frozen 2200-count budget the log-gain observation is high-SNR (obs-noise std ≈0.02 vs gain std 0.1–0.4), so the temporal *prior* has little leverage on the scene channel except where the OU autocovariance is structurally wrong (regime). The domain extension's clearest value is the **measurement-adequacy** statement: you cannot certify a scalar-OU medium model on regime-switching / heavy-tailed / quasi-periodic drift, and a learned prior can restore adequacy for the autocovariance-structural case.

---

## 3. E3 — certified-ceiling demo

`figs/e3_certified_ceiling.png` (data `json/e3_ceiling.json`): the fast-strong corner cell, every method on a **scene-PSNR vs gain-path-fidelity** axis with the oracle known-gain arm as the ceiling (dashed) and the medium-channel Cramér–Rao sd floor (t_c 12.9%, CV 6.3%) as the annotated referee.

| method | scene PSNR | gain-path log-RMSE |
|---|---|---|
| no-corr floor | 6.16 | 0.396 |
| DLGI frozen | 8.80 | 0.281 |
| DLGI + E[a]=1 fix (physics) | 10.98 | 0.042 |
| DLGI fix + learned gain corrector | 11.94 | 0.053 |
| learned corrector (SCGI, in-family) | 12.28 | 0.095 |
| oracle known-gain (ceiling) | 18.39 | 0 |

The physics fix collapses the gain-path RMSE (0.281→0.042, the level error) and moves DLGI decisively toward the ceiling; the learned tools add the last increment of *shape*. The picture the program wants: **DL races the ceiling, the oracle/CRB referees.**

---

## 4. Honest caveats

- **Exploratory, not preregistered.** DEV scenes only, no blinding, no confirmatory split. Nothing here should enter the sealed campaign.
- **E1's fix is real and free**; the finding that most of yesterday's "SCGI beats DLGI at the corner" gap was a fixable gauge bug (not a learned-method advantage) is the load-bearing result. It leaves the medium-measurement thesis and the family-robustness argument untouched.
- **E1(c) and E2's learned components do not universally dominate** (E1c regresses 3/4 cells; E2's prior loses on quasi-periodic and its predictor regresses adequacy off-regime). Reported straight — the wins are class/corner-specific.
- **E2 uses two isolation protocols by design:** true-scene z for the clean model-adequacy statement (§2b), full recon-loop for the honest scene PSNR (§2c) after the true-scene version was found to leak (PSNR above oracle) — that discarded artifact is documented, not hidden.
- **Photon-budget dependence:** the temporal-prior's small scene-channel leverage is specific to the high-SNR 2200-count regime; at lower flux the OU misfit would cost more and the learned prior would likely earn more (untested here; flagged for a follow-up sweep).
- **Faithful-in-spirit** SCGI reference reused verbatim from yesterday's `DL_EXPLORATORY`; not code-identical to Peng & Chen.

---

## 5. Files (all under `results/round63_next/DL_TOOL_EXPLORATORY/`)

- `DL_TOOL_REPORT.md` — this report.
- `code/dl_tool_common.py` — configurable copy of the frozen estimator (gauge / Jensen / E[a]=1 renorm switches) + smoother-piece exposers. dl_common untouched.
- `code/e1_diagnose.py`, `code/e1_fix_headtohead.py`, `code/e1c_residual_corrector.py` — E1 (a)/(b)/(c).
- `code/e2_common.py` — non-OU generators + adequacy (Ljung-Box) + learned-autocovariance readout.
- `code/e2_ou_misfit.py` — E2(b). `code/e2_learned_prior.py` — E2(c–e) adequacy + true-scene (superseded PSNR). `code/e2_recon_psnr.py` — E2 honest recon-loop scene/gain channel.
- `code/e3_ceiling.py`, `code/make_figs.py` — E3 + E1/E2 summary figures.
- `figs/e1_headtohead.png`, `figs/e2_learned_prior.png`, `figs/e3_certified_ceiling.png`.
- `json/*.json` — all per-cell / per-class aggregates + per-record entries.
- `checkpoints/e1c_corrector.pt`, `e2_smoother.pt`, `e2_predictor.pt`, `e2_smoother_recon.pt` — trained models.

Reproduce (from `code/`, `pytorch` env): `python e1_diagnose.py` → `python e1_fix_headtohead.py` → `python e1c_residual_corrector.py` → `python e2_ou_misfit.py` → `python e2_learned_prior.py` → `python e2_recon_psnr.py` → `python e3_ceiling.py` → `python make_figs.py`.
