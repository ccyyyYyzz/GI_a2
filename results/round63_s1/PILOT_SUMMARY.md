# ROUND63 S1 DEV PILOT SUMMARY

**S1 is exploratory and development-only** (spec D2 §1): dev images and seeds are disjoint from all S2 confirmatory images and seeds; these results do not enter confirmatory intervals or main tables.

- rows: 840   images: 6   seeds: [0, 1]   arms: RQL, POISSON-LIN, SAT-POISSON, PRECORRECT, GI
- side 64, bern50, M 2048, nu [20, 50, 100, 200, 500, 1000, 2000], rho_bar [0.05, 0.6], tau 50 ns, sigma_b 0, active start
- select_iter 60, fista_iters 200, F1 discrepancy lam_TV rule + descriptive audit (spec D2 §4)

## 1. Time-to-quality speedup S_j (dev)

S_j = T_opt(rho=0.05, Q*_j) / T_opt(rho=0.6, Q*_j), T_opt = M_physical*T. Q*_j = PSNR_rad[RQL, rho=0.05, nu=2000] - 1 dB (same target for every arm). 'cens' = target never reached; '--' = missing.

| image | Q*_j (dB) | RQL S_j | POISSON-LIN S_j | SAT-POISSON S_j | PRECORRECT S_j |
|---|---|---|---|---|---|
| dev_stl_00 | 17.54 | 1.000 | cens | 1.000 | cens |
| dev_stl_01 | 14.91 | 1.000 | cens | 1.000 | cens |
| dev_stl_02 | 13.99 | 1.000 | cens | 1.000 | cens |
| dev_stl_03 | 13.87 | 1.000 | 1.000 | 1.000 | cens |
| dev_text | 13.40 | 8.689 | cens | 8.878 | 8.014 |
| dev_fine_lines | 5.11 | 1.000 | 1.000 | 1.000 | cens |
| **median** | | 1.000 | 1.000 | 1.000 | 8.014 |

Censored (image count per arm): RQL=0, POISSON-LIN=4, SAT-POISSON=0, PRECORRECT=5

## 2. Descriptive measurement audit by (rho, nu) — RQL cells

Continuous diagnostics only (round-5 ruling: no binary gate).

| rho | nu | n | D_ratio med | D_ratio max | leak |
|---|---|---|---|---|---|
| 0.05 | 20 | 12 | 0.966 | 1.093 | 0 |
| 0.05 | 50 | 12 | 1.004 | 1.086 | 0 |
| 0.05 | 100 | 12 | 1.013 | 1.109 | 0 |
| 0.05 | 200 | 12 | 0.999 | 1.115 | 0 |
| 0.05 | 500 | 12 | 0.967 | 1.122 | 0 |
| 0.05 | 1000 | 12 | 1.072 | 1.190 | 0 |
| 0.05 | 2000 | 12 | 0.993 | 1.339 | 0 |
| 0.6 | 20 | 12 | 0.974 | 1.128 | 0 |
| 0.6 | 50 | 12 | 1.068 | 1.194 | 0 |
| 0.6 | 100 | 12 | 1.025 | 1.143 | 0 |
| 0.6 | 200 | 12 | 1.052 | 1.238 | 0 |
| 0.6 | 500 | 12 | 1.037 | 1.837 | 0 |
| 0.6 | 1000 | 12 | 1.014 | 2.264 | 0 |
| 0.6 | 2000 | 12 | 1.082 | 2.923 | 0 |

## 3. Runtime calibration per (nu, arm)  [sizes Colab shards]

| nu | arm | n | wall mean (s) | wall max (s) | select mean (s) | select max (s) |
|---|---|---|---|---|---|---|
| 20 | GI | 24 | 0.05 | 0.05 | 0.00 | 0.00 |
| 20 | POISSON-LIN | 24 | 105.74 | 133.93 | 105.74 | 133.93 |
| 20 | PRECORRECT | 24 | 58.77 | 63.91 | 58.77 | 63.91 |
| 20 | RQL | 24 | 156.23 | 173.06 | 156.14 | 171.89 |
| 20 | SAT-POISSON | 24 | 121.86 | 132.22 | 121.86 | 132.22 |
| 50 | GI | 24 | 0.05 | 0.06 | 0.00 | 0.00 |
| 50 | POISSON-LIN | 24 | 112.38 | 131.30 | 112.37 | 131.30 |
| 50 | PRECORRECT | 24 | 57.80 | 63.36 | 57.80 | 63.36 |
| 50 | RQL | 24 | 160.81 | 173.98 | 160.72 | 172.85 |
| 50 | SAT-POISSON | 24 | 124.62 | 132.10 | 124.62 | 132.10 |
| 100 | GI | 24 | 0.05 | 0.06 | 0.00 | 0.00 |
| 100 | POISSON-LIN | 24 | 114.07 | 137.86 | 114.07 | 137.86 |
| 100 | PRECORRECT | 24 | 58.11 | 65.87 | 58.10 | 65.86 |
| 100 | RQL | 24 | 162.21 | 174.28 | 162.12 | 174.28 |
| 100 | SAT-POISSON | 24 | 125.06 | 133.47 | 125.06 | 133.47 |
| 200 | GI | 24 | 0.05 | 0.06 | 0.00 | 0.00 |
| 200 | POISSON-LIN | 24 | 116.50 | 134.55 | 116.50 | 134.55 |
| 200 | PRECORRECT | 24 | 53.88 | 62.76 | 53.88 | 62.76 |
| 200 | RQL | 24 | 150.48 | 174.07 | 150.39 | 174.07 |
| 200 | SAT-POISSON | 24 | 128.15 | 136.13 | 128.15 | 136.13 |
| 500 | GI | 24 | 0.05 | 0.06 | 0.00 | 0.00 |
| 500 | POISSON-LIN | 24 | 132.94 | 144.08 | 132.94 | 144.08 |
| 500 | PRECORRECT | 24 | 52.29 | 56.83 | 52.29 | 56.83 |
| 500 | RQL | 24 | 148.29 | 182.03 | 148.20 | 181.02 |
| 500 | SAT-POISSON | 24 | 113.20 | 133.63 | 113.20 | 133.63 |
| 1000 | GI | 24 | 0.05 | 0.06 | 0.00 | 0.00 |
| 1000 | POISSON-LIN | 24 | 109.28 | 142.93 | 109.28 | 142.93 |
| 1000 | PRECORRECT | 24 | 43.28 | 55.24 | 43.28 | 55.24 |
| 1000 | RQL | 24 | 130.56 | 183.38 | 130.47 | 182.30 |
| 1000 | SAT-POISSON | 24 | 94.93 | 136.55 | 94.93 | 136.56 |
| 2000 | GI | 24 | 0.05 | 0.05 | 0.00 | 0.00 |
| 2000 | POISSON-LIN | 24 | 112.61 | 140.55 | 112.61 | 140.55 |
| 2000 | PRECORRECT | 24 | 42.88 | 52.18 | 42.88 | 52.18 |
| 2000 | RQL | 24 | 125.01 | 176.54 | 124.91 | 176.54 |
| 2000 | SAT-POISSON | 24 | 88.48 | 133.92 | 88.47 | 133.91 |

