# ROUND63 shard manifest summary (F1 freeze artifact, D2.3)

Generator: `code/round63/make_manifests.py` (deterministic; no wall-clock, no RNG). Target shard size: 6.00 h.

Cost model: S1 runtime calibration `results/round63_s1_passB/runtime_calibration.json` (smoke=False, side=64, scale x1, 45 (nu,arm) keys) — REAL 64^2/M=4096 walls INCLUDING the analytic lam_TV selection + descriptive audit; M-scaled by (M/4096); analytic fallback (RQL 300s / iter 150s / linear 5s @ nu=500, x(nu/500)^0.15), hadpair x2. 128^2 est is a nominal (M-scaled, blocked) figure — the matrix-free 128^2 operator is a prerequisite gate and its real per-fit cost is substantially higher than this lower bound.

Cohorts (GPT round-7 ruling §1): **S2A_DETAIL** = DETAIL-24 primary confirmatory positive-speed cohort (imageset detail24, 24 targets); **S2A_NATURAL** = STL-test NATURAL-24 secondary regime-boundary cohort (imageset conf, 30 names, NO positive gate). Both run the identical 64^2/bern50/M=4096 physics grid (rho in {0.05,0.3,0.6,1,2} x nu in {5,10,20,50,100,200,500,1000,2000} x 5 seeds).


## Default (runnable) shards

| shard | stage | pri | cells | est h | account | blocked |
|---|---|---|---|---|---|---|
| S2A_DETAIL_00 | S2A_DETAIL | 1 | 8 | 5.76 | pro1 |  |
| S2A_DETAIL_01 | S2A_DETAIL | 1 | 7 | 5.34 | pro2 |  |
| S2A_DETAIL_02 | S2A_DETAIL | 1 | 7 | 5.78 | pro1 |  |
| S2A_DETAIL_03 | S2A_DETAIL | 1 | 7 | 5.93 | pro2 |  |
| S2A_DETAIL_04 | S2A_DETAIL | 1 | 6 | 5.30 | pro1 |  |
| S2A_DETAIL_05 | S2A_DETAIL | 1 | 6 | 5.55 | pro2 |  |
| S2A_DETAIL_06 | S2A_DETAIL | 1 | 7 | 5.59 | pro1 |  |
| S2A_DETAIL_07 | S2A_DETAIL | 1 | 9 | 5.79 | pro2 |  |
| S2A_DETAIL_08 | S2A_DETAIL | 1 | 8 | 5.41 | pro1 |  |
| S2A_DETAIL_09 | S2A_DETAIL | 1 | 8 | 5.60 | pro2 |  |
| S2A_DETAIL_10 | S2A_DETAIL | 1 | 8 | 5.82 | pro1 |  |
| S2A_DETAIL_11 | S2A_DETAIL | 1 | 7 | 5.30 | pro2 |  |
| S2A_DETAIL_12 | S2A_DETAIL | 1 | 8 | 5.77 | pro1 |  |
| S2A_DETAIL_13 | S2A_DETAIL | 1 | 7 | 5.35 | pro2 |  |
| S2A_DETAIL_14 | S2A_DETAIL | 1 | 7 | 5.64 | pro1 |  |
| S2A_DETAIL_15 | S2A_DETAIL | 1 | 7 | 5.89 | pro2 |  |
| S2A_DETAIL_16 | S2A_DETAIL | 1 | 6 | 5.23 | pro1 |  |
| S2A_DETAIL_17 | S2A_DETAIL | 1 | 6 | 5.47 | pro2 |  |
| S2A_DETAIL_18 | S2A_DETAIL | 1 | 6 | 5.63 | pro1 |  |
| S2A_DETAIL_19 | S2A_DETAIL | 1 | 8 | 5.76 | pro2 |  |
| S2A_DETAIL_20 | S2A_DETAIL | 1 | 7 | 5.34 | pro1 |  |
| S2A_DETAIL_21 | S2A_DETAIL | 1 | 7 | 5.78 | pro2 |  |
| S2A_DETAIL_22 | S2A_DETAIL | 1 | 7 | 5.92 | pro1 |  |
| S2A_DETAIL_23 | S2A_DETAIL | 1 | 6 | 5.30 | pro2 |  |
| S2A_DETAIL_24 | S2A_DETAIL | 1 | 6 | 5.55 | pro1 |  |
| S2A_DETAIL_25 | S2A_DETAIL | 1 | 7 | 5.58 | pro2 |  |
| S2A_DETAIL_26 | S2A_DETAIL | 1 | 9 | 5.79 | pro1 |  |
| S2A_DETAIL_27 | S2A_DETAIL | 1 | 8 | 5.41 | pro2 |  |
| S2A_DETAIL_28 | S2A_DETAIL | 1 | 8 | 5.60 | pro1 |  |
| S2A_DETAIL_29 | S2A_DETAIL | 1 | 8 | 5.82 | pro2 |  |
| S2A_DETAIL_30 | S2A_DETAIL | 1 | 7 | 5.30 | pro1 |  |
| S2A_DETAIL_31 | S2A_DETAIL | 1 | 2 | 1.54 | pro2 |  |
| S2A_NATURAL_00 | S2A_NATURAL | 2 | 6 | 5.29 | pro1 |  |
| S2A_NATURAL_01 | S2A_NATURAL | 2 | 6 | 5.73 | pro2 |  |
| S2A_NATURAL_02 | S2A_NATURAL | 2 | 6 | 5.95 | pro1 |  |
| S2A_NATURAL_03 | S2A_NATURAL | 2 | 5 | 5.19 | pro2 |  |
| S2A_NATURAL_04 | S2A_NATURAL | 2 | 5 | 5.29 | pro1 |  |
| S2A_NATURAL_05 | S2A_NATURAL | 2 | 5 | 5.47 | pro2 |  |
| S2A_NATURAL_06 | S2A_NATURAL | 2 | 5 | 5.68 | pro1 |  |
| S2A_NATURAL_07 | S2A_NATURAL | 2 | 5 | 5.83 | pro2 |  |
| S2A_NATURAL_08 | S2A_NATURAL | 2 | 6 | 5.39 | pro1 |  |
| S2A_NATURAL_09 | S2A_NATURAL | 2 | 7 | 5.66 | pro2 |  |
| S2A_NATURAL_10 | S2A_NATURAL | 2 | 7 | 5.86 | pro1 |  |
| S2A_NATURAL_11 | S2A_NATURAL | 2 | 6 | 5.21 | pro2 |  |
| S2A_NATURAL_12 | S2A_NATURAL | 2 | 6 | 5.29 | pro1 |  |
| S2A_NATURAL_13 | S2A_NATURAL | 2 | 6 | 5.51 | pro2 |  |
| S2A_NATURAL_14 | S2A_NATURAL | 2 | 6 | 5.67 | pro1 |  |
| S2A_NATURAL_15 | S2A_NATURAL | 2 | 6 | 5.49 | pro2 |  |
| S2A_NATURAL_16 | S2A_NATURAL | 2 | 6 | 5.56 | pro1 |  |
| S2A_NATURAL_17 | S2A_NATURAL | 2 | 6 | 5.72 | pro2 |  |
| S2A_NATURAL_18 | S2A_NATURAL | 2 | 5 | 5.14 | pro1 |  |
| S2A_NATURAL_19 | S2A_NATURAL | 2 | 5 | 5.22 | pro2 |  |
| S2A_NATURAL_20 | S2A_NATURAL | 2 | 5 | 5.35 | pro1 |  |
| S2A_NATURAL_21 | S2A_NATURAL | 2 | 5 | 5.56 | pro2 |  |
| S2A_NATURAL_22 | S2A_NATURAL | 2 | 5 | 5.77 | pro1 |  |
| S2A_NATURAL_23 | S2A_NATURAL | 2 | 5 | 5.88 | pro2 |  |
| S2A_NATURAL_24 | S2A_NATURAL | 2 | 6 | 5.29 | pro1 |  |
| S2A_NATURAL_25 | S2A_NATURAL | 2 | 6 | 5.73 | pro2 |  |
| S2A_NATURAL_26 | S2A_NATURAL | 2 | 6 | 5.95 | pro1 |  |
| S2A_NATURAL_27 | S2A_NATURAL | 2 | 5 | 5.18 | pro2 |  |
| S2A_NATURAL_28 | S2A_NATURAL | 2 | 5 | 5.29 | pro1 |  |
| S2A_NATURAL_29 | S2A_NATURAL | 2 | 5 | 5.47 | pro2 |  |
| S2A_NATURAL_30 | S2A_NATURAL | 2 | 5 | 5.68 | pro1 |  |
| S2A_NATURAL_31 | S2A_NATURAL | 2 | 5 | 5.83 | pro2 |  |
| S2A_NATURAL_32 | S2A_NATURAL | 2 | 6 | 5.39 | pro1 |  |
| S2A_NATURAL_33 | S2A_NATURAL | 2 | 7 | 5.66 | pro2 |  |
| S2A_NATURAL_34 | S2A_NATURAL | 2 | 7 | 5.86 | pro1 |  |
| S2A_NATURAL_35 | S2A_NATURAL | 2 | 6 | 5.21 | pro2 |  |
| S2A_NATURAL_36 | S2A_NATURAL | 2 | 6 | 5.29 | pro1 |  |
| S2A_NATURAL_37 | S2A_NATURAL | 2 | 6 | 5.51 | pro2 |  |
| S2A_NATURAL_38 | S2A_NATURAL | 2 | 6 | 5.67 | pro1 |  |
| S2A_NATURAL_39 | S2A_NATURAL | 2 | 3 | 2.89 | pro2 |  |
| S2B_00 | S2B | 3 | 40 | 5.84 | pro1 |  |
| S2B_01 | S2B | 3 | 14 | 3.15 | pro2 |  |
| S2C_00 | S2C | 3 | 14 | 5.75 | pro1 |  |
| S2C_01 | S2C | 3 | 15 | 5.86 | pro2 |  |
| S2C_02 | S2C | 3 | 25 | 5.21 | pro1 |  |
| S3_00 | S3 | 2 | 32 | 5.99 | pro1 |  |
| S3_01 | S3 | 2 | 32 | 5.99 | pro2 |  |
| S3_02 | S3 | 2 | 32 | 5.99 | pro1 |  |
| S3_03 | S3 | 2 | 24 | 4.49 | pro2 |  |
| S3JIT_00 | S3JIT | 4 | 12 | 2.25 | pro1 |  |

## Blocked shards (EXCLUDED from the default run set)

| shard | stage | pri | cells | est h | account | blocked_by |
|---|---|---|---|---|---|---|
| S2C128_00 | S2C128 | 3 | 7 | 5.56 | pro1 | matrixfree_128 |
| S2C128_01 | S2C128 | 3 | 7 | 5.61 | pro2 | matrixfree_128 |
| S2C128_02 | S2C128 | 3 | 7 | 5.65 | pro1 | matrixfree_128 |
| S2C128_03 | S2C128 | 3 | 6 | 4.98 | pro2 | matrixfree_128 |
| S3CONT_00 | S3CONT | 2 | 18 | 3.37 | pro1 | continuous_eta_inheritance |

## Per-stage totals

| stage | pri | shards | cells | rows | est h | blocked |
|---|---|---|---|---|---|---|
| S2A_DETAIL | 1 | 32 | 225 | 15840 | 174.88 |  |
| S2A_NATURAL | 2 | 40 | 225 | 19800 | 218.59 |  |
| S3 | 2 | 4 | 120 | 2880 | 22.45 |  |
| S3CONT | 2 | 1 | 18 | 432 | 3.37 | continuous_eta_inheritance |
| S2B | 3 | 2 | 54 | 2592 | 8.99 |  |
| S2C | 3 | 3 | 54 | 2592 | 16.82 |  |
| S2C128 | 3 | 4 | 27 | 972 | 21.79 | matrixfree_128 |
| S3JIT | 4 | 1 | 12 | 288 | 2.25 |  |

## Account parity (frozen: even shard index -> pro1, odd -> pro2; default shards)

| account | shards | est h |
|---|---|---|
| pro1 | 42 | 231.12 |
| pro2 | 40 | 212.85 |

## Notes

- S3 (+S3JIT/S3CONT) stays at the historic mismatch anchor **64^2 / bern50 / M=2048 / nu=500 (M/n=0.5)** per spec §3 — it is NOT lifted to the S2-A M=4096 grid.
- S2B drops its M=4096 rows: they are a strict subset of the S2A_NATURAL anchor (same physics, subset of images/arms/seeds); S2B therefore runs only M in {1024, 2048}.
- No cell carries C0/c_force; the frozen `C0_FROZEN.json` (C0*=inf, Pass A) is authoritative at run time.

## S4 (EXCLUDED from manifests)

- S4(i) scalar exact-Fisher map (nu=20..2000, rho adaptive to max(64, 2 rho*(nu))): standalone analytic-derivative script, not a cell grid.
- S4(ii) 8^2/16^2 exact-vs-RQL image split (rho in {.03,.1,.3,.6,1,2} x nu in {20,100,500,2000}, 3 seeds): BLOCKED on the EXACT+TV arm (same TV + selector as RQL); emit once that arm exists.
- S4/Supplement optional: nu in {1,2} exact Bernoulli/trinomial vs RQL + ceiling-count fraction (outside the deployment zone; no nu<5 result enters either confirmatory endpoint).
