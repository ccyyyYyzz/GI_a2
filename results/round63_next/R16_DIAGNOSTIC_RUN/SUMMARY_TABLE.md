# R16 frozen diagnostic - summary table (2026-07-22)

N_MC=2,000,000 frames/load; nu=2000; 200-replicate frame-bootstrap peak CIs;
estimators/peaks byte-identical to the audited prep driver. Raw (N,L) per frame
in raw_nl_cv*.npz. Point curves + per-load J values in r16_part_{A,B}.json.
NO VERDICT DECLARED HERE - 5.4 adjudication belongs to the coordinator.

## Peak locations (point | bootstrap median [95% CI])

| channel | estimator | point peak (quad/cubic) | boot median [95% CI] | R16 pred | measured (hist, 150k/400k) |
|---|---|---|---|---|---|
| cv=0.02 | J_MI (score) | 10.1098 / 10.0998 | 10.0998 [10.0998, 10.3002] (sd 0.0956) | 10.22 | 9.62 |
| cv=0.02 | J_CS (score x-check) | 10.2155 / - | 10.2046 [9.4214, 10.7802] (sd 0.3431) | 10.22 | 9.62 |
| cv=0.02 | J_HIST (audited, dlog=0.01, a=0.5) | 9.2876 / 9.2998 | 9.2998 [9.2998, 10.9075] (sd 0.6747) | 10.22 | 9.62 |
| cv=0.20 | J_MI (score) | 2.1520 / 2.1500 | 2.1500 [2.1251, 2.1999] (sd 0.0155) | 2.17 | 2.08 |
| cv=0.20 | J_CS (score x-check) | 2.0981 / - | 2.0997 [2.0227, 2.1952] (sd 0.0435) | 2.17 | 2.08 |
| cv=0.20 | J_HIST (audited, dlog=0.01, a=0.5) | 2.0556 / 2.0500 | 2.2001 [2.0500, 2.2513] (sd 0.0876) | 2.17 | 2.08 |
| c=0 cal (exact 22.2543) | J_MI (score) | 22.2372 / 22.2497 | 22.2497 [21.9994, 22.7496] (sd 0.2681) | 22.2543 (exact) | - |
| c=0 cal (exact 22.2543) | J_CS (score x-check) | 22.5181 / - | 22.5216 [21.0603, 24.5000] (sd 0.8266) | 22.2543 (exact) | - |
| c=0 cal (exact 22.2543) | J_HIST (audited, dlog=0.01, a=0.5) | 24.2093 / 24.2502 | 22.2497 [20.0000, 24.2502] (sd 1.5454) | 22.2543 (exact) | - |

## Histogram-peak sensitivity sweeps (three-point quadratic peak)

- cv=0.02: alpha sweep (dlog=0.01): a=0.0:9.2880 a=0.05:9.2872 a=0.1:9.2873 a=0.5:9.2876 a=1.0:9.2877
- cv=0.02: dlog sweep (alpha=0.5): dlog=0.0025:10.1036 dlog=0.005:9.2902 dlog=0.01:9.2876 dlog=0.02:9.2890
- cv=0.20: alpha sweep (dlog=0.01): a=0.0:2.0553 a=0.05:2.2307 a=0.1:2.2349 a=0.5:2.0556 a=1.0:2.0555
- cv=0.20: dlog sweep (alpha=0.5): dlog=0.0025:2.2404 dlog=0.005:2.2333 dlog=0.01:2.0556 dlog=0.02:2.0538
- c=0 cal: alpha sweep (dlog=0.01): a=0.0:20.0000 a=0.05:24.2181 a=0.1:24.2152 a=0.5:24.2093 a=1.0:24.2071
- c=0 cal: dlog sweep (alpha=0.5): dlog=0.0025:20.5705 dlog=0.005:22.9792 dlog=0.01:24.2093 dlog=0.02:24.2252

## R16 5.2 agreement condition + tail handling

- cv=0.02: max|J_MI-J_CS|/J = 0.294% ; max tail_mass(bins<50) = 0.000069
- cv=0.20: max|J_MI-J_CS|/J = 0.226% ; max tail_mass(bins<50) = 0.000211
- c=0 cal: max|J_MI-J_CS|/J = 0.216% ; max tail_mass(bins<50) = 0.000040

## 5.4 pattern evidence (facts only; no verdict)

1. Score estimators vs the R16 predictions: cv=0.02 MI cubic boot CI [10.0998, 10.3002]
   contains 10.22 and excludes 9.62; cv=0.20 MI cubic boot CI [2.1251, 2.1999]
   contains 2.17 and excludes 2.08. CS agrees (wider CIs, same containment).
2. Score estimator vs measured values: it does NOT reproduce 9.62 / 2.08 at either end point.
3. c=0 calibration: J_MI peak 22.2372/22.2497 vs exact 22.2543 (-0.08%/-0.02%),
   boot CI [21.9994, 22.7496] contains the exact value. Same-draws J_HIST point peak
   24.2093/24.2502 (+8.8%) with boot CI [20.0000, 24.2502] spanning nearly the grid.
4. Histogram-peak stability: boot sd 0.675 (cv=0.02), 0.088 (cv=0.20), 1.545 (c=0)
   vs MI sd 0.096 / 0.016 / 0.268 - the histogram peak location is 5-10x noisier on
   these flat curves. Across sample size the audited-estimator peak moved:
   9.645 (150k) -> 9.610 (400k) -> 9.288 point / CI [9.27,10.93] (2M) at cv=0.02.
5. alpha (tail-floor) lever at 2M frames: inert at cv=0.02 and c=0 (peak flat to 4th
   decimal for a in {0,...,1}); non-monotone +-0.18 jitter at cv=0.20. Tail mass
   (bins<50) is <= 2.1e-4 everywhere, so the floor rarely binds at this scale.
6. dlog sweep is not a smooth dlog^2 drift; the peak jumps between locations
   (cv=0.02: 10.10 at dlog=0.0025 vs 9.29 at >=0.005; c=0: 20.57/22.98/24.21/24.23).

Caveat: MI point peaks sit ~1% below the R16 predictions (10.10 vs 10.22; 2.15 vs 2.17)
with CIs containing the predictions; the c=0 anchor bounds MI's own bias at <0.1%.
