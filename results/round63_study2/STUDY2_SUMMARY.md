# ROUND63 STUDY-2 SUMMARY — contrast-dead-time phase diagram

Primary confirmatory arm k=16 on the DETAIL-32 cohort (six families x four instances = 24), reconstructed with RQL through the frozen production path; round-7 Q90 acquisition-speed endpoint, family-stratified nested bootstrap (ruling §5/§6).

- rows: 10800   images: 24   rho_safe=0.05   rho_fast=0.6   nu_budget=2000.0

## Primary gate (ruling §6) — DETAIL_SPEED_PASS = False

- median S_gate = 6.796  (>=3: True)
- bootstrap 2.5% LB of median S_gate = 2.835  (>1: True)
- images with S_gate>1 = 16/24  (>=18: False)

## Fixed-budget quality gain (ruling §7, secondary)

- median DeltaQ = 1.16 dB  (>=1.0: True)
- bootstrap 2.5% LB of median DeltaQ = 0.72  (>0: True)
- images with DeltaQ>0 = 19/24  (>=18: True)   pass = True

Censoring status counts: FAST_RIGHT_CENSORED=2, NORMAL=17, SAFE_RANGE_UNINFORMATIVE=5

## Per-image endpoints (RQL, Q90)

| image | family | status | S_gate | T_safe | T_fast | R_j (dB) | Q90_j (dB) | DeltaQ_budget (dB) | seeds |
|---|---|---|---|---|---|---|---|---|---|
| detail32_glyph_0 | glyph | NORMAL | 32.516 | 0.07768 | 0.002389 | 2.98 | 10.88 | 1.11 | 5 |
| detail32_glyph_1 | glyph | NORMAL | 9.278 | 0.02088 | 0.002251 | 4.23 | 14.34 | 0.80 | 5 |
| detail32_glyph_2 | glyph | NORMAL | 21.339 | 0.07423 | 0.003478 | 3.54 | 11.75 | 0.90 | 5 |
| detail32_glyph_3 | glyph | NORMAL | 8.653 | 0.02169 | 0.002506 | 4.29 | 14.41 | 0.70 | 5 |
| detail32_chirp_0 | chirp | NORMAL | 5.013 | 0.09432 | 0.01882 | 1.23 | 7.17 | 2.30 | 5 |
| detail32_chirp_1 | chirp | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.00 | -- | 0.35 | 5 |
| detail32_chirp_2 | chirp | NORMAL | 1.713 | 0.08805 | 0.05141 | 2.42 | 8.14 | 2.88 | 5 |
| detail32_chirp_3 | chirp | NORMAL | 10.775 | 0.08003 | 0.007427 | 3.14 | 9.71 | 2.14 | 5 |
| detail32_spokes_0 | spokes | NORMAL | 13.232 | 0.07777 | 0.005877 | 2.64 | 9.27 | 1.21 | 5 |
| detail32_spokes_1 | spokes | NORMAL | 11.094 | 0.06092 | 0.005491 | 2.92 | 9.90 | 1.46 | 5 |
| detail32_spokes_2 | spokes | NORMAL | 10.737 | 0.0812 | 0.007562 | 3.93 | 10.86 | 4.38 | 5 |
| detail32_spokes_3 | spokes | NORMAL | 13.695 | 0.07038 | 0.005139 | 2.79 | 9.29 | 1.40 | 5 |
| detail32_maze_0 | maze | NORMAL | 17.712 | 0.07522 | 0.004247 | 7.93 | 14.53 | 7.62 | 5 |
| detail32_maze_1 | maze | NORMAL | 4.691 | 0.06769 | 0.01443 | 6.86 | 13.07 | 3.14 | 5 |
| detail32_maze_2 | maze | NORMAL | 9.227 | 0.08182 | 0.008868 | 13.06 | 20.89 | 9.70 | 5 |
| detail32_maze_3 | maze | NORMAL | 8.580 | 0.07372 | 0.008592 | 16.70 | 24.84 | 8.24 | 5 |
| detail32_contour_0 | contour | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.04 | -- | 1.01 | 5 |
| detail32_contour_1 | contour | NORMAL | 3.022 | 0.0796 | 0.02634 | 2.10 | 10.61 | 0.02 | 5 |
| detail32_contour_2 | contour | FAST_RIGHT_CENSORED | 0.000 | 0.04777 | -- | 1.43 | 8.48 | -0.32 | 5 |
| detail32_contour_3 | contour | NORMAL | 0.723 | 0.04604 | 0.06368 | 2.76 | 9.77 | 1.35 | 5 |
| detail32_microtexture_0 | microtexture | FAST_RIGHT_CENSORED | 0.000 | 0.09554 | -- | 1.38 | 14.23 | -0.98 | 5 |
| detail32_microtexture_1 | microtexture | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.00 | -- | -1.64 | 5 |
| detail32_microtexture_2 | microtexture | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.00 | -- | -0.11 | 5 |
| detail32_microtexture_3 | microtexture | SAFE_RANGE_UNINFORMATIVE | 1.000 | -- | -- | 0.32 | -- | -1.49 | 5 |

