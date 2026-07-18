# ROUND63 S1 PILOT SUMMARY (round-7 Q90 endpoint, cohort=dev)

**S1 is exploratory and development-only** (spec D2 §1): dev images and seeds are disjoint from all S2 confirmatory images and seeds; these results do not enter confirmatory intervals or main tables. The analysis layer here is the round-7 FINAL endpoint reused verbatim for S2.

- rows: 1080   images: 6   seeds: [0, 1]   arms: RQL, POISSON-LIN, SAT-POISSON, PRECORRECT, GI   mode: pass_b   cohort: dev
- side 64, bern50, M 4096, nu [5, 10, 20, 50, 100, 200, 500, 1000, 2000], rho_bar [0.05, 0.6], tau 50 ns, sigma_b 0, active start
- select_iter 60, fista_iters 200, C0 inf, analytic score-concentration lam_TV rule + descriptive audit

## 1. Acquisition-speed endpoint S_gate — RQL (round-7 §2 Q90)

Per (image, rho): seed-mean PSNR_rad at each nu -> equal-weight PAVA isotonic fit vs log T_opt (raw curve kept in the JSON). Safe = rho=0.05; fast = rho=0.6. R_j = Qiso_safe(T_max) - Qiso_safe(T_min); if R_j<0.50 dB the safe range is uninformative (S_gate=1, not positive); else Q90_j = Qiso_safe(T_min) + 0.90*R_j and crossing times come by log-T interpolation on the isotonic curves. Censoring: NORMAL(=T_s/T_f)/FAST_LEFT_CENSORED(>=T_s/T_min)/BOTH_LEFT_CENSORED(=1, unresolved)/SAFE_LEFT_FAST_INTERIOR(<1)/FAST_RIGHT_CENSORED(fail=0)/SAFE_RANGE_UNINFORMATIVE(=1)/ANALYSIS_FAILURE(0). No image deleted.

| image | family | status | S_gate | T_safe | T_fast | R_j (dB) | Q90_j (dB) | DeltaQ_budget (dB) | seeds |
|---|---|---|---|---|---|---|---|---|---|
| dev_stl_00 | dev_stl_00 | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.37 | -- | 0.43 | 2 |
| dev_stl_01 | dev_stl_01 | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.26 | -- | 0.71 | 2 |
| dev_stl_02 | dev_stl_02 | NORMAL | 7.777 | 0.3101 | 0.03987 | 0.88 | 15.82 | 1.52 | 2 |
| dev_stl_03 | dev_stl_03 | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.33 | -- | 1.53 | 2 |
| dev_text | dev_text | NORMAL | 7.669 | 0.3413 | 0.04451 | 0.78 | 13.77 | 2.32 | 2 |
| dev_fine_lines | dev_fine_lines | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.11 | -- | 0.07 | 2 |
| **median S_gate** | | | 1.000 | | | | | | |

Status counts: NORMAL=2, SAFE_RANGE_UNINFORMATIVE=4

### Primary gate (§2.5) — DETAIL_SPEED_PASS = null (non-detail cohort: no positive gate, §4)

- median S_gate = 1.000  (>=3: False)
- bootstrap 2.5% LB of median S_gate = 1.000  (>1: False)
- images with S_gate>1 (excl. non-positive classes) = 2/6  (>=18: False)

### Fixed-budget quality gain (§3, key secondary)

DeltaQ_budget = RAW seed-mean PSNR_rad(fast) - (safe) at nu=2000 (NOT isotonic). Descriptive success: median>=1.0 dB AND LB2.5>0 AND >=18 images with DeltaQ>0.

- median DeltaQ = 1.12 dB  (>=1.0: True)
- bootstrap 2.5% LB of median DeltaQ = 0.09  (>0: True)
- images with DeltaQ>0 = 6/6  (>=18: False)   pass = null (non-detail cohort)

Bootstrap: B=500 draws, stream rng_for(0, 63, 6, 2), stratified by image (degenerate); finite draws S_gate=500 / DeltaQ=500.

## 2. Descriptive measurement audit by (rho, nu) — RQL cells

Continuous diagnostics only (round-5 ruling: no binary gate).

| rho | nu | n | D_ratio med | D_ratio max | leak |
|---|---|---|---|---|---|
| 0.05 | 5 | 12 | 1.004 | 1.091 | 0 |
| 0.05 | 10 | 12 | 1.033 | 1.111 | 1 |
| 0.05 | 20 | 12 | 0.984 | 1.008 | 0 |
| 0.05 | 50 | 12 | 0.989 | 1.112 | 0 |
| 0.05 | 100 | 12 | 0.989 | 1.101 | 0 |
| 0.05 | 200 | 12 | 1.046 | 1.115 | 0 |
| 0.05 | 500 | 12 | 1.021 | 1.185 | 0 |
| 0.05 | 1000 | 12 | 1.005 | 1.187 | 0 |
| 0.05 | 2000 | 12 | 0.994 | 1.412 | 0 |
| 0.6 | 5 | 12 | 0.986 | 1.045 | 0 |
| 0.6 | 10 | 12 | 0.994 | 1.100 | 0 |
| 0.6 | 20 | 12 | 0.989 | 1.059 | 0 |
| 0.6 | 50 | 12 | 1.019 | 1.135 | 0 |
| 0.6 | 100 | 12 | 1.022 | 1.162 | 0 |
| 0.6 | 200 | 12 | 1.020 | 1.248 | 0 |
| 0.6 | 500 | 12 | 1.025 | 1.478 | 0 |
| 0.6 | 1000 | 12 | 1.034 | 1.830 | 0 |
| 0.6 | 2000 | 12 | 1.046 | 2.721 | 0 |

## 3. Runtime calibration per (nu, arm)  [sizes Colab shards]

| nu | arm | n | wall mean (s) | wall max (s) | select mean (s) | select max (s) |
|---|---|---|---|---|---|---|
| 5 | GI | 24 | 0.12 | 0.17 | 0.00 | 0.00 |
| 5 | POISSON-LIN | 24 | 3.85 | 12.97 | 3.85 | 12.97 |
| 5 | PRECORRECT | 24 | 2.58 | 2.83 | 2.58 | 2.83 |
| 5 | RQL | 24 | 91.21 | 99.22 | 91.20 | 99.11 |
| 5 | SAT-POISSON | 24 | 6.43 | 15.91 | 6.43 | 15.91 |
| 10 | GI | 24 | 0.12 | 0.16 | 0.00 | 0.00 |
| 10 | POISSON-LIN | 24 | 8.92 | 16.10 | 8.92 | 16.10 |
| 10 | PRECORRECT | 24 | 2.66 | 2.98 | 2.66 | 2.98 |
| 10 | RQL | 24 | 97.91 | 110.13 | 97.89 | 110.13 |
| 10 | SAT-POISSON | 24 | 5.25 | 7.11 | 5.25 | 7.11 |
| 20 | GI | 24 | 0.12 | 0.16 | 0.00 | 0.00 |
| 20 | POISSON-LIN | 24 | 8.20 | 15.53 | 8.20 | 15.53 |
| 20 | PRECORRECT | 24 | 2.60 | 2.93 | 2.60 | 2.93 |
| 20 | RQL | 24 | 97.96 | 107.66 | 97.95 | 107.66 |
| 20 | SAT-POISSON | 24 | 5.69 | 7.81 | 5.68 | 7.80 |
| 50 | GI | 24 | 0.13 | 0.17 | 0.00 | 0.00 |
| 50 | POISSON-LIN | 24 | 10.25 | 17.68 | 10.25 | 17.68 |
| 50 | PRECORRECT | 24 | 2.73 | 3.49 | 2.72 | 3.48 |
| 50 | RQL | 24 | 103.58 | 115.75 | 103.57 | 115.75 |
| 50 | SAT-POISSON | 24 | 6.82 | 10.99 | 6.82 | 10.99 |
| 100 | GI | 24 | 0.13 | 0.15 | 0.00 | 0.00 |
| 100 | POISSON-LIN | 24 | 10.88 | 17.97 | 10.88 | 17.96 |
| 100 | PRECORRECT | 24 | 2.64 | 3.18 | 2.64 | 3.18 |
| 100 | RQL | 24 | 104.40 | 120.24 | 104.39 | 120.24 |
| 100 | SAT-POISSON | 24 | 7.23 | 10.96 | 7.23 | 10.96 |
| 200 | GI | 24 | 0.12 | 0.15 | 0.00 | 0.00 |
| 200 | POISSON-LIN | 24 | 11.70 | 20.23 | 11.69 | 20.23 |
| 200 | PRECORRECT | 24 | 2.54 | 2.82 | 2.54 | 2.82 |
| 200 | RQL | 24 | 106.07 | 125.81 | 106.06 | 125.63 |
| 200 | SAT-POISSON | 24 | 7.99 | 10.41 | 7.99 | 10.40 |
| 500 | GI | 24 | 0.12 | 0.15 | 0.00 | 0.00 |
| 500 | POISSON-LIN | 24 | 12.42 | 18.04 | 12.42 | 18.04 |
| 500 | PRECORRECT | 24 | 2.48 | 2.80 | 2.48 | 2.80 |
| 500 | RQL | 24 | 109.79 | 131.88 | 109.78 | 131.75 |
| 500 | SAT-POISSON | 24 | 8.53 | 14.89 | 8.53 | 14.89 |
| 1000 | GI | 24 | 0.12 | 0.15 | 0.00 | 0.00 |
| 1000 | POISSON-LIN | 24 | 12.82 | 18.11 | 12.82 | 18.11 |
| 1000 | PRECORRECT | 24 | 2.45 | 3.12 | 2.45 | 3.11 |
| 1000 | RQL | 24 | 112.33 | 136.39 | 112.32 | 136.23 |
| 1000 | SAT-POISSON | 24 | 10.81 | 13.54 | 10.81 | 13.54 |
| 2000 | GI | 24 | 0.12 | 0.16 | 0.00 | 0.00 |
| 2000 | POISSON-LIN | 24 | 12.82 | 17.67 | 12.82 | 17.67 |
| 2000 | PRECORRECT | 24 | 2.46 | 3.30 | 2.46 | 3.30 |
| 2000 | RQL | 24 | 115.43 | 138.59 | 115.42 | 138.59 |
| 2000 | SAT-POISSON | 24 | 10.28 | 14.28 | 10.27 | 14.28 |

