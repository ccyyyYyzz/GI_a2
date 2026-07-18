# ROUND63 S1 PILOT SUMMARY (round-7 Q90 endpoint, cohort=detail)

**S1 is exploratory and development-only** (spec D2 §1): dev images and seeds are disjoint from all S2 confirmatory images and seeds; these results do not enter confirmatory intervals or main tables. The analysis layer here is the round-7 FINAL endpoint reused verbatim for S2.

- rows: 15840   images: 6   seeds: [0, 1]   arms: RQL, POISSON-LIN, SAT-POISSON, PRECORRECT, GI   mode: normal   cohort: detail
- side 64, bern50, M 4096, nu [5, 10, 20, 50, 100, 200, 500, 1000, 2000], rho_bar [0.05, 0.6], tau 50 ns, sigma_b 0, active start
- select_iter 60, fista_iters 200, C0 inf, analytic score-concentration lam_TV rule + descriptive audit

## 1. Acquisition-speed endpoint S_gate — RQL (round-7 §2 Q90)

Per (image, rho): seed-mean PSNR_rad at each nu -> equal-weight PAVA isotonic fit vs log T_opt (raw curve kept in the JSON). Safe = rho=0.05; fast = rho=0.6. R_j = Qiso_safe(T_max) - Qiso_safe(T_min); if R_j<0.50 dB the safe range is uninformative (S_gate=1, not positive); else Q90_j = Qiso_safe(T_min) + 0.90*R_j and crossing times come by log-T interpolation on the isotonic curves. Censoring: NORMAL(=T_s/T_f)/FAST_LEFT_CENSORED(>=T_s/T_min)/BOTH_LEFT_CENSORED(=1, unresolved)/SAFE_LEFT_FAST_INTERIOR(<1)/FAST_RIGHT_CENSORED(fail=0)/SAFE_RANGE_UNINFORMATIVE(=1)/ANALYSIS_FAILURE(0). No image deleted.

| image | family | status | S_gate | T_safe | T_fast | R_j (dB) | Q90_j (dB) | DeltaQ_budget (dB) | seeds |
|---|---|---|---|---|---|---|---|---|---|
| detail_chirp_0 | chirp | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.12 | -- | 0.12 | 5 |
| detail_chirp_1 | chirp | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.16 | -- | 0.30 | 5 |
| detail_chirp_2 | chirp | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.04 | -- | 0.27 | 5 |
| detail_chirp_3 | chirp | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.09 | -- | 0.12 | 5 |
| detail_contour_0 | contour | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.45 | -- | 0.83 | 5 |
| detail_contour_1 | contour | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.37 | -- | 0.91 | 5 |
| detail_contour_2 | contour | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.13 | -- | 0.53 | 5 |
| detail_contour_3 | contour | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.25 | -- | 0.63 | 5 |
| detail_glyph_0 | glyph | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.11 | -- | 0.25 | 5 |
| detail_glyph_1 | glyph | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.35 | -- | 0.34 | 5 |
| detail_glyph_2 | glyph | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.14 | -- | 0.35 | 5 |
| detail_glyph_3 | glyph | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.06 | -- | 0.34 | 5 |
| detail_maze_0 | maze | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.01 | -- | 0.17 | 5 |
| detail_maze_1 | maze | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.01 | -- | 0.18 | 5 |
| detail_maze_2 | maze | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.02 | -- | 0.12 | 5 |
| detail_maze_3 | maze | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.01 | -- | 0.14 | 5 |
| detail_microtexture_0 | microtexture | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.02 | -- | 0.05 | 5 |
| detail_microtexture_1 | microtexture | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.00 | -- | 0.04 | 5 |
| detail_microtexture_2 | microtexture | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.01 | -- | 0.21 | 5 |
| detail_microtexture_3 | microtexture | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.02 | -- | 0.16 | 5 |
| detail_spokes_0 | spokes | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.27 | -- | 0.64 | 5 |
| detail_spokes_1 | spokes | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.24 | -- | 0.62 | 5 |
| detail_spokes_2 | spokes | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.29 | -- | 1.91 | 5 |
| detail_spokes_3 | spokes | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.30 | -- | 0.45 | 5 |
| **median S_gate** | | | 1.000 | | | | | | |

Status counts: SAFE_RANGE_UNINFORMATIVE=24

### Primary gate (§2.5) — DETAIL_SPEED_PASS = False

- median S_gate = 1.000  (>=3: False)
- bootstrap 2.5% LB of median S_gate = 1.000  (>1: False)
- images with S_gate>1 (excl. non-positive classes) = 0/24  (>=18: False)

### Fixed-budget quality gain (§3, key secondary)

DeltaQ_budget = RAW seed-mean PSNR_rad(fast) - (safe) at nu=2000 (NOT isotonic). Descriptive success: median>=1.0 dB AND LB2.5>0 AND >=18 images with DeltaQ>0.

- median DeltaQ = 0.28 dB  (>=1.0: False)
- bootstrap 2.5% LB of median DeltaQ = 0.23  (>0: True)
- images with DeltaQ>0 = 24/24  (>=18: True)   pass = False

Bootstrap: B=10000 draws, stream rng_for(0, 63, 6, 2), stratified by family; finite draws S_gate=10000 / DeltaQ=10000.

## 2. Descriptive measurement audit by (rho, nu) — RQL cells

Continuous diagnostics only (round-5 ruling: no binary gate).

| rho | nu | n | D_ratio med | D_ratio max | leak |
|---|---|---|---|---|---|
| 0.05 | 5 | 120 | 0.988 | 1.195 | 4 |
| 0.05 | 10 | 120 | 0.991 | 1.205 | 2 |
| 0.05 | 20 | 120 | 0.993 | 1.185 | 2 |
| 0.05 | 50 | 120 | 0.994 | 1.148 | 2 |
| 0.05 | 100 | 120 | 0.998 | 1.119 | 1 |
| 0.05 | 200 | 120 | 1.007 | 1.109 | 0 |
| 0.05 | 500 | 120 | 1.011 | 1.127 | 3 |
| 0.05 | 1000 | 120 | 1.008 | 1.172 | 1 |
| 0.05 | 2000 | 120 | 1.026 | 1.261 | 0 |
| 0.3 | 5 | 120 | 0.994 | 1.111 | 0 |
| 0.3 | 10 | 120 | 1.004 | 1.134 | 3 |
| 0.3 | 20 | 120 | 1.000 | 1.119 | 0 |
| 0.3 | 50 | 120 | 1.007 | 1.161 | 1 |
| 0.3 | 100 | 120 | 1.021 | 1.113 | 0 |
| 0.3 | 200 | 120 | 1.023 | 1.177 | 0 |
| 0.3 | 500 | 120 | 1.046 | 1.253 | 0 |
| 0.3 | 1000 | 120 | 1.074 | 1.337 | 0 |
| 0.3 | 2000 | 120 | 1.133 | 1.699 | 0 |
| 0.6 | 5 | 120 | 0.997 | 1.090 | 0 |
| 0.6 | 10 | 120 | 0.996 | 1.131 | 1 |
| 0.6 | 20 | 120 | 0.999 | 1.121 | 1 |
| 0.6 | 50 | 120 | 1.014 | 1.137 | 0 |
| 0.6 | 100 | 120 | 1.020 | 1.142 | 0 |
| 0.6 | 200 | 120 | 1.035 | 1.192 | 0 |
| 0.6 | 500 | 120 | 1.077 | 1.335 | 0 |
| 0.6 | 1000 | 120 | 1.078 | 1.511 | 0 |
| 0.6 | 2000 | 120 | 1.135 | 1.666 | 0 |
| 1 | 5 | 120 | 0.975 | 1.088 | 0 |
| 1 | 10 | 120 | 1.020 | 1.160 | 3 |
| 1 | 20 | 120 | 1.011 | 1.146 | 1 |
| 1 | 50 | 120 | 1.009 | 1.145 | 0 |
| 1 | 100 | 120 | 1.029 | 1.176 | 0 |
| 1 | 200 | 120 | 1.042 | 1.243 | 1 |
| 1 | 500 | 120 | 1.078 | 1.367 | 0 |
| 1 | 1000 | 120 | 1.101 | 1.475 | 0 |
| 1 | 2000 | 120 | 1.200 | 1.861 | 0 |
| 2 | 5 | 120 | 0.867 | 1.024 | 1 |
| 2 | 10 | 120 | 0.973 | 1.130 | 0 |
| 2 | 20 | 120 | 1.013 | 1.161 | 2 |
| 2 | 50 | 120 | 1.018 | 1.169 | 1 |
| 2 | 100 | 120 | 1.032 | 1.208 | 0 |
| 2 | 200 | 120 | 1.067 | 1.260 | 0 |
| 2 | 500 | 120 | 1.065 | 1.394 | 0 |
| 2 | 1000 | 120 | 1.136 | 1.579 | 0 |
| 2 | 2000 | 120 | 1.298 | 2.115 | 0 |

## 3. Runtime calibration per (nu, arm)  [sizes Colab shards]

| nu | arm | n | wall mean (s) | wall max (s) | select mean (s) | select max (s) |
|---|---|---|---|---|---|---|
| 5 | POISSON-LIN | 360 | 5.08 | 10.90 | 5.08 | 10.90 |
| 5 | PRECORRECT | 360 | 1.65 | 1.75 | 1.65 | 1.75 |
| 5 | RQL | 600 | 62.41 | 66.12 | 62.41 | 66.12 |
| 5 | SAT-POISSON | 360 | 3.13 | 6.23 | 3.13 | 6.23 |
| 10 | POISSON-LIN | 360 | 5.96 | 11.32 | 5.96 | 11.32 |
| 10 | PRECORRECT | 360 | 1.65 | 1.78 | 1.65 | 1.78 |
| 10 | RQL | 600 | 64.51 | 69.07 | 64.50 | 69.07 |
| 10 | SAT-POISSON | 360 | 3.63 | 6.32 | 3.63 | 6.32 |
| 20 | POISSON-LIN | 360 | 7.16 | 11.51 | 7.16 | 11.51 |
| 20 | PRECORRECT | 360 | 1.66 | 1.78 | 1.66 | 1.78 |
| 20 | RQL | 600 | 66.08 | 71.87 | 66.08 | 71.87 |
| 20 | SAT-POISSON | 360 | 4.24 | 6.44 | 4.24 | 6.44 |
| 50 | POISSON-LIN | 360 | 8.03 | 11.48 | 8.03 | 11.48 |
| 50 | PRECORRECT | 360 | 1.64 | 1.80 | 1.64 | 1.80 |
| 50 | RQL | 600 | 67.65 | 73.54 | 67.65 | 73.54 |
| 50 | SAT-POISSON | 360 | 4.73 | 6.74 | 4.73 | 6.74 |
| 100 | GI | 240 | 0.08 | 0.09 | 0.00 | 0.00 |
| 100 | POISSON-LIN | 360 | 8.46 | 11.55 | 8.46 | 11.55 |
| 100 | PRECORRECT | 360 | 1.61 | 1.74 | 1.61 | 1.74 |
| 100 | RQL | 600 | 68.80 | 73.40 | 68.80 | 73.40 |
| 100 | SAT-POISSON | 360 | 4.88 | 6.30 | 4.88 | 6.30 |
| 200 | POISSON-LIN | 360 | 8.90 | 11.96 | 8.90 | 11.96 |
| 200 | PRECORRECT | 360 | 1.62 | 1.73 | 1.62 | 1.73 |
| 200 | RQL | 600 | 70.07 | 75.30 | 70.07 | 75.30 |
| 200 | SAT-POISSON | 360 | 5.35 | 6.78 | 5.35 | 6.77 |
| 500 | GI | 240 | 0.08 | 0.09 | 0.00 | 0.00 |
| 500 | POISSON-LIN | 360 | 9.30 | 12.80 | 9.30 | 12.80 |
| 500 | PRECORRECT | 360 | 1.63 | 1.77 | 1.63 | 1.77 |
| 500 | RQL | 600 | 71.45 | 80.91 | 71.45 | 80.91 |
| 500 | SAT-POISSON | 360 | 5.96 | 8.80 | 5.96 | 8.80 |
| 1000 | POISSON-LIN | 360 | 9.32 | 11.49 | 9.32 | 11.48 |
| 1000 | PRECORRECT | 360 | 1.59 | 1.72 | 1.59 | 1.72 |
| 1000 | RQL | 600 | 75.94 | 85.59 | 75.94 | 85.59 |
| 1000 | SAT-POISSON | 360 | 5.55 | 8.95 | 5.55 | 8.95 |
| 2000 | GI | 240 | 0.08 | 0.09 | 0.00 | 0.00 |
| 2000 | POISSON-LIN | 360 | 8.84 | 11.72 | 8.84 | 11.72 |
| 2000 | PRECORRECT | 360 | 1.60 | 2.43 | 1.60 | 2.43 |
| 2000 | RQL | 600 | 80.12 | 91.11 | 80.12 | 91.11 |
| 2000 | SAT-POISSON | 360 | 7.08 | 8.71 | 7.08 | 8.71 |

