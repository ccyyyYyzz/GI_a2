# ROUND63 S1 DEV PILOT SUMMARY (round-6 protocol)

**S1 is exploratory and development-only** (spec D2 §1): dev images and seeds are disjoint from all S2 confirmatory images and seeds; these results do not enter confirmatory intervals or main tables.

- rows: 1080   images: 6   seeds: [0, 1]   arms: RQL, POISSON-LIN, SAT-POISSON, PRECORRECT, GI   mode: pass_b
- side 64, bern50, M 4096, nu [5, 10, 20, 50, 100, 200, 500, 1000, 2000], rho_bar [0.05, 0.6], tau 50 ns, sigma_b 0, active start
- select_iter 60, fista_iters 200, C0 inf, analytic score-concentration lam_TV rule + descriptive audit

## 1. Time-to-quality speedup S_gate — RQL (round-6 §4 censoring)

T_opt = M_physical*nu*tau; T_min at the shortest nu tier. Q*_j = PSNR_rad[RQL, rho=0.05, nu=2000] - 1 dB. S_gate per the frozen taxonomy A(normal)/B(fast-left-censored, lower bound)/C(both-left-censored, =1, unresolved below floor)/D(safe-left/fast-interior, <1)/E(fast-right-censored, fail=0)/F(analysis failure=0).

| image | status | S_gate | T_safe | T_fast | Q*_j (dB) | seeds |
|---|---|---|---|---|---|---|
| dev_stl_00 | BOTH_LEFT_CENSORED | 1.000 | 0.001024 | 0.001024 | 17.34 | 2 |
| dev_stl_01 | SAFE_LEFT_FAST_INTERIOR | 0.756 | 0.001024 | 0.001355 | 15.23 | 2 |
| dev_stl_02 | BOTH_LEFT_CENSORED | 1.000 | 0.001024 | 0.001024 | 14.90 | 2 |
| dev_stl_03 | BOTH_LEFT_CENSORED | 1.000 | 0.001024 | 0.001024 | 14.12 | 2 |
| dev_text | BOTH_LEFT_CENSORED | 1.000 | 0.001024 | 0.001024 | 12.85 | 2 |
| dev_fine_lines | BOTH_LEFT_CENSORED | 1.000 | 0.001024 | 0.001024 | 5.17 | 2 |
| **median S_gate** | | 1.000 | | | | |

Status counts: BOTH_LEFT_CENSORED=5, SAFE_LEFT_FAST_INTERIOR=1

Gate stats (S1 dev, descriptive — NOT a confirmatory gate): median S_gate=1.000; images with S_gate>1 (excl. C/E/F)=0/6; image-level stratified bootstrap (B=2000) 2.5% lower bound of the median = 0.878.

Reference gate check (median>=3: False; bootstrap_lower>1: False) — reported for orientation only; S1 never enters the confirmatory decision.

Secondary arms (point-estimate median S_gate): POISSON-LIN=0.245 (n>1=0), SAT-POISSON=1.000 (n>1=0), PRECORRECT=0.000 (n>1=1)

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

