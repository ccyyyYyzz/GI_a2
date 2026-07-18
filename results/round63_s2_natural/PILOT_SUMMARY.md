# ROUND63 S1 PILOT SUMMARY (round-7 Q90 endpoint, cohort=natural)

**S1 is exploratory and development-only** (spec D2 §1): dev images and seeds are disjoint from all S2 confirmatory images and seeds; these results do not enter confirmatory intervals or main tables. The analysis layer here is the round-7 FINAL endpoint reused verbatim for S2.

- rows: 19800   images: 6   seeds: [0, 1]   arms: RQL, POISSON-LIN, SAT-POISSON, PRECORRECT, GI   mode: normal   cohort: natural
- side 64, bern50, M 4096, nu [5, 10, 20, 50, 100, 200, 500, 1000, 2000], rho_bar [0.05, 0.6], tau 50 ns, sigma_b 0, active start
- select_iter 60, fista_iters 200, C0 inf, analytic score-concentration lam_TV rule + descriptive audit

## 1. Acquisition-speed endpoint S_gate — RQL (round-7 §2 Q90)

Per (image, rho): seed-mean PSNR_rad at each nu -> equal-weight PAVA isotonic fit vs log T_opt (raw curve kept in the JSON). Safe = rho=0.05; fast = rho=0.6. R_j = Qiso_safe(T_max) - Qiso_safe(T_min); if R_j<0.50 dB the safe range is uninformative (S_gate=1, not positive); else Q90_j = Qiso_safe(T_min) + 0.90*R_j and crossing times come by log-T interpolation on the isotonic curves. Censoring: NORMAL(=T_s/T_f)/FAST_LEFT_CENSORED(>=T_s/T_min)/BOTH_LEFT_CENSORED(=1, unresolved)/SAFE_LEFT_FAST_INTERIOR(<1)/FAST_RIGHT_CENSORED(fail=0)/SAFE_RANGE_UNINFORMATIVE(=1)/ANALYSIS_FAILURE(0). No image deleted.

| image | family | status | S_gate | T_safe | T_fast | R_j (dB) | Q90_j (dB) | DeltaQ_budget (dB) | seeds |
|---|---|---|---|---|---|---|---|---|---|
| fine_lines | fine_lines | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.09 | -- | 0.07 | 5 |
| gray_staircase | gray_staircase | NORMAL | 4.939 | 0.313 | 0.06338 | 2.73 | 12.16 | 3.11 | 5 |
| low_contrast | low_contrast | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.00 | -- | -1.78 | 5 |
| sparse_dots | sparse_dots | NORMAL | 7.645 | 0.3155 | 0.04128 | 6.16 | 26.14 | 3.00 | 5 |
| stl_00 | stl_00 | NORMAL | 5.991 | 0.37 | 0.06177 | 0.77 | 13.22 | 1.25 | 5 |
| stl_01 | stl_01 | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.10 | -- | 0.29 | 5 |
| stl_02 | stl_02 | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.02 | -- | 0.86 | 5 |
| stl_03 | stl_03 | NORMAL | 5.711 | 0.3447 | 0.06036 | 1.33 | 16.23 | 1.25 | 5 |
| stl_04 | stl_04 | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.30 | -- | 1.42 | 5 |
| stl_05 | stl_05 | NORMAL | 7.209 | 0.3621 | 0.05023 | 0.57 | 14.77 | 1.25 | 5 |
| stl_06 | stl_06 | NORMAL | 6.987 | 0.3568 | 0.05107 | 0.93 | 13.23 | 1.51 | 5 |
| stl_07 | stl_07 | NORMAL | 6.761 | 0.344 | 0.05088 | 1.91 | 13.71 | 2.57 | 5 |
| stl_08 | stl_08 | NORMAL | 10.916 | 0.2978 | 0.02728 | 1.59 | 12.80 | 2.46 | 5 |
| stl_09 | stl_09 | NORMAL | 6.289 | 0.3677 | 0.05846 | 0.94 | 15.77 | 2.30 | 5 |
| stl_10 | stl_10 | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.00 | -- | 0.37 | 5 |
| stl_11 | stl_11 | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.41 | -- | 1.14 | 5 |
| stl_12 | stl_12 | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.43 | -- | 1.03 | 5 |
| stl_13 | stl_13 | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.28 | -- | 1.35 | 5 |
| stl_14 | stl_14 | NORMAL | 6.617 | 0.3131 | 0.04732 | 0.94 | 14.49 | 1.62 | 5 |
| stl_15 | stl_15 | NORMAL | 7.804 | 0.3301 | 0.0423 | 0.80 | 11.89 | 2.26 | 5 |
| stl_16 | stl_16 | NORMAL | 6.980 | 0.3547 | 0.05081 | 1.07 | 13.80 | 1.46 | 5 |
| stl_17 | stl_17 | NORMAL | 8.375 | 0.3617 | 0.04319 | 0.53 | 13.00 | 1.58 | 5 |
| stl_18 | stl_18 | NORMAL | 5.295 | 0.2284 | 0.04313 | 1.02 | 14.47 | 1.85 | 5 |
| stl_19 | stl_19 | NORMAL | 6.568 | 0.3557 | 0.05416 | 1.73 | 12.18 | 2.33 | 5 |
| stl_20 | stl_20 | NORMAL | 5.692 | 0.3342 | 0.05871 | 1.15 | 14.67 | 1.28 | 5 |
| stl_21 | stl_21 | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.34 | -- | 0.77 | 5 |
| stl_22 | stl_22 | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.36 | -- | 1.63 | 5 |
| stl_23 | stl_23 | NORMAL | 7.803 | 0.3186 | 0.04083 | 1.03 | 12.52 | 1.74 | 5 |
| text | text | NORMAL | 7.834 | 0.3393 | 0.04331 | 0.77 | 13.77 | 2.02 | 5 |
| usaf_bars | usaf_bars | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.38 | -- | 0.95 | 5 |
| **median S_gate** | | | 5.701 | | | | | | |

Status counts: NORMAL=18, SAFE_RANGE_UNINFORMATIVE=12

### Primary gate (§2.5) — DETAIL_SPEED_PASS = null (non-detail cohort: no positive gate, §4)

- median S_gate = 5.701  (>=3: True)
- bootstrap 2.5% LB of median S_gate = 1.000  (>1: False)
- images with S_gate>1 (excl. non-positive classes) = 18/30  (>=18: True)

### Fixed-budget quality gain (§3, key secondary)

DeltaQ_budget = RAW seed-mean PSNR_rad(fast) - (safe) at nu=2000 (NOT isotonic). Descriptive success: median>=1.0 dB AND LB2.5>0 AND >=18 images with DeltaQ>0.

- median DeltaQ = 1.44 dB  (>=1.0: True)
- bootstrap 2.5% LB of median DeltaQ = 1.18  (>0: True)
- images with DeltaQ>0 = 29/30  (>=18: True)   pass = null (non-detail cohort)

Bootstrap: B=10000 draws, stream rng_for(0, 63, 6, 2), stratified by image (degenerate); finite draws S_gate=10000 / DeltaQ=10000.

## 2. Descriptive measurement audit by (rho, nu) — RQL cells

Continuous diagnostics only (round-5 ruling: no binary gate).

| rho | nu | n | D_ratio med | D_ratio max | leak |
|---|---|---|---|---|---|
| 0.05 | 5 | 150 | 0.985 | 1.219 | 6 |
| 0.05 | 10 | 150 | 0.992 | 1.194 | 4 |
| 0.05 | 20 | 150 | 0.991 | 1.185 | 3 |
| 0.05 | 50 | 150 | 0.987 | 1.156 | 2 |
| 0.05 | 100 | 150 | 0.994 | 1.135 | 1 |
| 0.05 | 200 | 150 | 1.012 | 1.274 | 0 |
| 0.05 | 500 | 150 | 1.008 | 1.532 | 2 |
| 0.05 | 1000 | 150 | 0.996 | 2.186 | 1 |
| 0.05 | 2000 | 150 | 1.000 | 2.938 | 1 |
| 0.3 | 5 | 150 | 0.999 | 1.116 | 2 |
| 0.3 | 10 | 150 | 1.009 | 1.158 | 3 |
| 0.3 | 20 | 150 | 1.002 | 1.181 | 1 |
| 0.3 | 50 | 150 | 1.001 | 1.306 | 2 |
| 0.3 | 100 | 150 | 1.012 | 1.635 | 2 |
| 0.3 | 200 | 150 | 0.996 | 2.003 | 1 |
| 0.3 | 500 | 150 | 1.008 | 3.301 | 0 |
| 0.3 | 1000 | 150 | 1.012 | 5.455 | 1 |
| 0.3 | 2000 | 150 | 1.015 | 9.176 | 0 |
| 0.6 | 5 | 150 | 0.992 | 1.096 | 1 |
| 0.6 | 10 | 150 | 0.999 | 1.176 | 1 |
| 0.6 | 20 | 150 | 0.996 | 1.290 | 1 |
| 0.6 | 50 | 150 | 1.007 | 1.456 | 0 |
| 0.6 | 100 | 150 | 1.018 | 1.828 | 0 |
| 0.6 | 200 | 150 | 1.013 | 2.543 | 1 |
| 0.6 | 500 | 150 | 1.018 | 4.618 | 1 |
| 0.6 | 1000 | 150 | 1.003 | 8.265 | 0 |
| 0.6 | 2000 | 150 | 1.042 | 13.638 | 0 |
| 1 | 5 | 150 | 0.981 | 1.098 | 0 |
| 1 | 10 | 150 | 1.017 | 1.170 | 3 |
| 1 | 20 | 150 | 1.005 | 1.293 | 1 |
| 1 | 50 | 150 | 1.000 | 1.631 | 0 |
| 1 | 100 | 150 | 1.011 | 2.289 | 0 |
| 1 | 200 | 150 | 0.999 | 3.162 | 1 |
| 1 | 500 | 150 | 1.021 | 6.103 | 0 |
| 1 | 1000 | 150 | 1.013 | 10.038 | 0 |
| 1 | 2000 | 150 | 1.063 | 17.620 | 0 |
| 2 | 5 | 150 | 0.874 | 1.060 | 2 |
| 2 | 10 | 150 | 0.974 | 1.167 | 2 |
| 2 | 20 | 150 | 1.004 | 1.302 | 2 |
| 2 | 50 | 150 | 1.013 | 1.794 | 1 |
| 2 | 100 | 150 | 1.007 | 2.587 | 1 |
| 2 | 200 | 150 | 1.019 | 4.033 | 0 |
| 2 | 500 | 150 | 0.999 | 6.902 | 0 |
| 2 | 1000 | 150 | 1.038 | 7.359 | 0 |
| 2 | 2000 | 150 | 1.146 | 20.307 | 0 |

## 3. Runtime calibration per (nu, arm)  [sizes Colab shards]

| nu | arm | n | wall mean (s) | wall max (s) | select mean (s) | select max (s) |
|---|---|---|---|---|---|---|
| 5 | POISSON-LIN | 450 | 5.05 | 10.75 | 5.05 | 10.75 |
| 5 | PRECORRECT | 450 | 1.65 | 1.76 | 1.65 | 1.76 |
| 5 | RQL | 750 | 61.99 | 65.28 | 61.99 | 65.28 |
| 5 | SAT-POISSON | 450 | 3.16 | 5.94 | 3.16 | 5.94 |
| 10 | POISSON-LIN | 450 | 6.04 | 11.75 | 6.04 | 11.75 |
| 10 | PRECORRECT | 450 | 1.66 | 1.86 | 1.66 | 1.86 |
| 10 | RQL | 750 | 64.36 | 70.43 | 64.35 | 70.43 |
| 10 | SAT-POISSON | 450 | 3.68 | 6.41 | 3.68 | 6.41 |
| 20 | POISSON-LIN | 450 | 7.00 | 11.52 | 7.00 | 11.52 |
| 20 | PRECORRECT | 450 | 1.63 | 1.76 | 1.63 | 1.76 |
| 20 | RQL | 750 | 64.89 | 69.34 | 64.89 | 69.34 |
| 20 | SAT-POISSON | 450 | 4.15 | 6.99 | 4.15 | 6.99 |
| 50 | POISSON-LIN | 450 | 8.00 | 11.90 | 8.00 | 11.90 |
| 50 | PRECORRECT | 450 | 1.63 | 1.74 | 1.63 | 1.74 |
| 50 | RQL | 750 | 66.81 | 77.29 | 66.81 | 77.29 |
| 50 | SAT-POISSON | 450 | 4.87 | 11.50 | 4.86 | 11.50 |
| 100 | GI | 300 | 0.08 | 0.09 | 0.00 | 0.00 |
| 100 | POISSON-LIN | 450 | 8.43 | 11.74 | 8.43 | 11.74 |
| 100 | PRECORRECT | 450 | 1.62 | 1.77 | 1.62 | 1.77 |
| 100 | RQL | 750 | 68.49 | 78.59 | 68.49 | 78.59 |
| 100 | SAT-POISSON | 450 | 4.96 | 11.57 | 4.96 | 11.57 |
| 200 | POISSON-LIN | 450 | 8.83 | 11.48 | 8.83 | 11.48 |
| 200 | PRECORRECT | 450 | 1.61 | 1.74 | 1.61 | 1.74 |
| 200 | RQL | 750 | 69.63 | 78.64 | 69.62 | 78.64 |
| 200 | SAT-POISSON | 450 | 5.33 | 11.32 | 5.33 | 11.32 |
| 500 | GI | 300 | 0.08 | 0.09 | 0.00 | 0.00 |
| 500 | POISSON-LIN | 450 | 9.24 | 12.10 | 9.24 | 12.10 |
| 500 | PRECORRECT | 450 | 1.62 | 1.75 | 1.61 | 1.75 |
| 500 | RQL | 750 | 70.77 | 81.37 | 70.77 | 81.37 |
| 500 | SAT-POISSON | 450 | 5.12 | 11.58 | 5.12 | 11.58 |
| 1000 | POISSON-LIN | 450 | 9.35 | 11.91 | 9.35 | 11.91 |
| 1000 | PRECORRECT | 450 | 1.59 | 1.74 | 1.59 | 1.74 |
| 1000 | RQL | 750 | 75.30 | 84.57 | 75.30 | 84.57 |
| 1000 | SAT-POISSON | 450 | 6.03 | 11.52 | 6.03 | 11.52 |
| 2000 | GI | 300 | 0.08 | 0.08 | 0.00 | 0.00 |
| 2000 | POISSON-LIN | 450 | 6.94 | 11.62 | 6.94 | 11.62 |
| 2000 | PRECORRECT | 450 | 1.59 | 1.71 | 1.59 | 1.71 |
| 2000 | RQL | 750 | 80.13 | 91.86 | 80.13 | 91.86 |
| 2000 | SAT-POISSON | 450 | 7.16 | 11.48 | 7.16 | 11.48 |

