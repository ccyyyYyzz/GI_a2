# Compute-Cost Forecast — 27-cell OA plan vs the 6 GPU-hour ceiling (R41 §4.3, §4.8) · **WITHIN**

**Engine** `mc_plan.py` (sha `2cb9396374c28cc8`) · data `COST_FORECAST.json` · measured on RTX 4060.

> **Total forecast: 1.21 GPU-hours** (ceiling **6.0**) · storage **0.001 GB** (ceiling 64 GB) ·
> **126 980 records**. Comfortably within both ceilings, with >4× headroom. Exceeding the ceiling at
> run time requires **stopping and reporting**, never silently trimming cells or records.

## Measured throughput (this box)

| generator | measured | planning value (conservative) |
|---|---|---|
| FIXED-COV Gaussian bank generator | ~0.78–2.0e6 record-banks/s (T-dependent; drops at T=4096) | **8.0e5 record-banks/s** |
| FRESH-COV-OPT per-bank (batched GPU, chunk 48) | **4.5 ms/bank** (measured) | 9.5 ms/bank (2× margin) |
| `setup_cell` (per-cell filter build) | 0.27 s | 0.27 s |

## The 27-cell orthogonal-array MC subset (R41 §4.3)

The exhaustive analytic map is the full **81 cells** (`3 σ_f × 3 spectrum × 3 k_w/k_p × 3 claim`). Monte
Carlo confirmation uses a precommitted **27-cell `3^(4-1)` orthogonal array**: `claim_level =
(σ_f_level + spectrum_level + k_w_level) mod 3`. **Every level of every factor appears 9 times and every
factor pair is balanced** (verified in `sealed_common.oa_grid_27`). Each primary point: **≥500 H0 + 500
H1** records; specificity and mismatch: **≥250 records/class**.

Primary-grid `T_eff` (per (cell, ε), `= min(T_det, 4096)`): min 84, median 2596, max 4096; 48 of 108
points hit the `T_eff=4096` cap (mostly ε=0.5% and 1% floor cells — the honestly-expensive small anomalies).

## Cost breakdown

| block | scope | record-banks | forecast (s) |
|---|---|---|---|
| **D0** mechanism | two-engine Fisher + physical Poisson shot/cov (few cells) | — | 300 |
| **D1/D2** primary | 27 OA cells × 4 ε × (500+500) records | 2.616e8 | **327** |
| **D3** specificity | 3 base cells × 6 classes × 250 @ T=2048 | 9.22e6 | 12 |
| **D5** mismatch | 2 base cells × 13 axis-levels × (250+250) @ T=800 | 1.12e7 | 14 |
| **D4** fresh-cov | analytic over 27 cells + MC validation at 6 points (40 rec) | 3.50e5 streams | 58 + **3324** |
| setup + overhead | per-cell filter builds + misc | — | 329 |
| **TOTAL** | | | **4364 s = 1.21 h** |

The dominant cost is the **D4 FRESH-COV-OPT MC validation** (3324 s ≈ 55 min at the conservative 9.5
ms/bank; ~26 min at the measured 4.5 ms/bank). The D4 **primary** comparison is analytic
(`fresh_cov_fisher`, `lam_eff = 0.1168 ≈ fixed 0.1117`, latency ratio 0.956 → retain-concentration
branch) and is cheap (58 s); the MC only validates the fresh arm at best/mid/floor × {2%,5%}.

## Storage

Records are reduced on the fly to **scalar score statistics** (per-record `t_beyond`, `t_amp`, `t_lag`,
`t_mean`) — never per-bank frames or M×M covariances. Total scalar storage **0.001 GB**. (If the full
M×M sample covariances were retained for audit it would be 16.6 GB — still under the 64 GB ceiling — but
they are **not** retained; only summaries + this forecast JSON.)

## Ceiling discipline

At 1.21 forecast GPU-hours the plan sits at ~20% of the 6-hour ceiling; even a 2× throughput shortfall
(≈2.4 h) or doubling the FRESH-COV MC scope stays within budget. If the live run projects over 6.0
GPU-hours or 64 GB, the run **halts and reports** — the OA subset and per-point record counts are
frozen and are not to be reduced to fit.
