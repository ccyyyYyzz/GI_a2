# QUANTILE-CALIBRATION RESCUE PROBE

**Status: POST-HOC DESCRIPTIVE EXPLORATION.** Not a preregistered endpoint, carries no PASS/FAIL verdict, and makes no paper edit. It reruns the frozen forward machinery (tag `m1-freeze` @ 6f00932) at `nu=2000` on the 5 frozen seeds and reports display/analysis only, on the post-unblinding-display precedent. n is small (5 negative scenes + 2 positive controls); the single-knob variants change only one scalar.

## 0. Frozen protocol (fixed before any outcome was read)

1. **Scenes.** The five scenes with negative frozen dQ from the R19 corrected `per_image_delta` (`M1_VERDICTS_SPEC_CORRECTED_R19.json`): `contour_0` (-3.612), `contour_1` (-5.794), `contour_2` (-4.738), `spokes_0` (-0.765), `spokes_1` (-5.085); plus two positive controls `maze_1` (+6.436, best case) and `chirp_3` (+1.539, median). *Deviation from the brief:* the brief's parenthetical said "microtexture negative", but all four microtexture instances are POSITIVE in the frozen table; the exact IDs were pulled from the per-image table as instructed, so the second spokes negative (`spokes_1`) is the fifth scene.
2. **Variants** (same single global multiplier, different calibration target). The three deployed SCAT32 modes share ONE 972-row multiset and differ only in the scalar load `rho_bar`; the finite-window ridge argmax at `nu=2000` is `rho_ridge = 22.2544` (= `rho_R` = 22.2545, clip does not bind):
   - **(a) FROZEN** — predicted MEAN load = `rho_R` (deployed policy; validates the harness).
   - **(b) Q90** — predicted q90 load = `rho_R`  (`rho_bar = rho_R * mean(pred)/q90(pred)`).
   - **(c) Q75** — predicted q75 load = `rho_R`.
   - **(d) MEAN-CLIPPED** — *NOT single-knob; upper reference only.* baseline multiplier + per-pattern amplitude servo clipping each predicted load at `rho_R` (the per-pattern servo R17 deliberately removed).
3. **Per scene x variant.** `nu=2000`, 5 frozen seeds, seed-mean `PSNR_rad` minus the frozen `SCAT32-060` terminal (the frozen dQ definition); plus achieved load quantiles and the fraction of the 972 patterns past the ridge argmax.
4. **Analysis.** Does (b)/(c) turn any negative dQ positive, and by how much? Do the controls survive? Is the improvement monotone in the past-ridge-fraction reduction (delta-dQ vs reduction correlation)?
5. **Stop rule.** If the harness cannot reproduce a frozen baseline dQ within 0.05 dB, STOP and report the discrepancy.

## 1. Harness validation (stop-rule gate)

Baseline (FROZEN) dQ reproduced through the literal `campaign.run_cell` vs the R19 corrected `per_image_delta`:

| scene | reproduced dQ | frozen dQ | \|delta\| |
|---|---|---|---|
| `m1_contour_0` | -3.61218 | -3.61218 | 0.00000 |
| `m1_contour_1` | -5.79446 | -5.79446 | 0.00000 |
| `m1_contour_2` | -4.73830 | -4.73830 | 0.00000 |
| `m1_spokes_0` | -0.76524 | -0.76524 | 0.00000 |
| `m1_spokes_1` | -5.08472 | -5.08472 | 0.00000 |
| `m1_maze_1` | +6.43642 | +6.43642 | 0.00000 |
| `m1_chirp_3` | +1.53884 | +1.53884 | 0.00000 |

**Max |delta| = 0.00000 dB < 0.05 dB -> PASS.** The MEAN-CLIPPED replica reproduces the unclipped `run_cell` baseline to < 5e-5 dB on every scene (seed 0), so the clipped path is trustworthy.

## 2. Scene x variant (nu=2000, seed-mean)

Achieved load quantiles are on the actual per-pattern dead-time load `load_i = rho_bar * (rows @ x_true)_i` (mean = `rho_R` = 22.25 for FROZEN, matching the frozen `rho_bar`); `past%` = fraction of the 972 patterns above the ridge argmax 22.25.

| scene | variant | PSNR_rad | dQ | delta-dQ | load q50 / q90 / q95 / max | past-ridge % |
|---|---|---|---|---|---|---|
| **m1_contour_0** (neg) | FROZEN | 8.857 | -3.612 | -- | 21.9 / 28.2 / 30.1 / 42.5 | 47.2% |
|  | Q90 | 10.118 | -2.352 | +1.261 | 16.6 / 21.4 / 22.9 / 32.2 | 6.7% |
|  | Q75 | 8.683 | -3.786 | -0.174 | 19.0 / 24.5 / 26.1 / 36.8 | 22.1% |
|  | MEAN_CLIPPED | 7.247 | -5.222 | -1.610 | 20.2 / 25.0 / 27.1 / 33.3 | 26.4% |
| **m1_contour_1** (neg) | FROZEN | 7.765 | -5.794 | -- | 21.9 / 28.2 / 30.2 / 39.7 | 46.3% |
|  | Q90 | 7.718 | -5.842 | -0.047 | 16.4 / 21.2 / 22.7 / 29.8 | 6.6% |
|  | Q75 | 7.714 | -5.845 | -0.051 | 19.0 / 24.6 / 26.3 / 34.6 | 23.2% |
|  | MEAN_CLIPPED | 7.115 | -6.445 | -0.650 | 20.1 / 24.8 / 26.2 / 32.7 | 26.2% |
| **m1_contour_2** (neg) | FROZEN | 10.008 | -4.738 | -- | 21.9 / 29.0 / 31.1 / 43.6 | 47.4% |
|  | Q90 | 9.957 | -4.788 | -0.050 | 16.2 / 21.5 / 23.0 / 32.2 | 7.1% |
|  | Q75 | 10.025 | -4.721 | +0.017 | 18.7 / 24.8 / 26.5 / 37.2 | 22.0% |
|  | MEAN_CLIPPED | 10.416 | -4.330 | +0.408 | 19.9 / 25.3 / 26.9 / 31.4 | 28.3% |
| **m1_spokes_0** (neg) | FROZEN | 9.769 | -0.765 | -- | 21.5 / 30.6 / 33.1 / 52.1 | 46.4% |
|  | Q90 | 9.777 | -0.757 | +0.008 | 16.4 / 23.4 / 25.3 / 39.9 | 14.3% |
|  | Q75 | 9.773 | -0.761 | +0.004 | 18.9 / 26.9 / 29.1 / 45.8 | 31.2% |
|  | MEAN_CLIPPED | 9.906 | -0.628 | +0.138 | 20.0 / 27.3 / 29.7 / 43.2 | 33.7% |
| **m1_spokes_1** (neg) | FROZEN | 9.752 | -5.085 | -- | 21.5 / 31.0 / 33.6 / 44.7 | 46.1% |
|  | Q90 | 9.769 | -5.068 | +0.017 | 16.7 / 24.1 / 26.1 / 34.7 | 15.7% |
|  | Q75 | 9.761 | -5.075 | +0.009 | 18.7 / 26.9 / 29.1 / 38.8 | 27.8% |
|  | MEAN_CLIPPED | 10.031 | -4.805 | +0.279 | 20.0 / 27.0 / 28.8 / 38.2 | 34.3% |
| **m1_maze_1** (ctrl) | FROZEN | 35.249 | +6.436 | -- | 20.7 / 37.8 / 43.0 / 77.2 | 43.3% |
|  | Q90 | 35.272 | +6.459 | +0.023 | 15.8 / 28.9 / 32.9 / 59.0 | 23.6% |
|  | Q75 | 35.060 | +6.247 | -0.189 | 18.1 / 33.1 / 37.7 / 67.5 | 33.3% |
|  | MEAN_CLIPPED | 35.748 | +6.935 | +0.499 | 19.2 / 33.1 / 37.3 / 54.6 | 37.0% |
| **m1_chirp_3** (ctrl) | FROZEN | 10.916 | +1.539 | -- | 22.0 / 28.8 / 30.4 / 41.2 | 47.9% |
|  | Q90 | 10.857 | +1.479 | -0.060 | 16.9 / 22.2 / 23.4 / 31.8 | 9.8% |
|  | Q75 | 10.923 | +1.546 | +0.007 | 19.3 / 25.3 / 26.7 / 36.2 | 25.1% |
|  | MEAN_CLIPPED | 11.079 | +1.701 | +0.162 | 20.4 / 25.2 / 27.0 / 31.9 | 30.1% |

## 3. Rescue outcome

- **No negative dQ is turned positive by any variant.** The best single-knob move is Q90 on `contour_0` (-3.612 -> -2.352, +1.26 dB) but it stays negative; the other four negatives move by < 0.06 dB under Q90/Q75. The genie upper-reference MEAN-CLIPPED also fails to rescue any negative (it *worsens* `contour_0`/`contour_1` by 0.65-1.61 dB and helps the rest by at most +0.41 dB).
- **Both positive controls survive** every variant: `maze_1` +6.44 -> +6.25..+6.94, `chirp_3` +1.54 -> +1.48..+1.70. The variants do not destroy successes.
- The Q90/Q75 knobs *do* pull the load distribution back under the ridge as designed (Q90 cuts the past-ridge fraction from ~46-48% to ~7-16%), but quality barely responds.

## 4. Mechanism correlation

delta-dQ (single-knob variant minus baseline) vs the reduction in past-ridge fraction, over the 14 (scene x {Q90,Q75}) points:

- **Pearson r = 0.405**, **Spearman rho = 0.169** — weak, non-monotone. `contour_1` and `contour_2` receive nearly the same Q90 past-ridge reduction as `contour_0` (0.40 vs 0.41) yet gain ~0 dB while `contour_0` gains +1.26 dB.
- Pulling patterns back under the finite-window ridge is therefore **not** the controlling variable for the contour/spokes deficit. This matches the independent `CONTOUR_DIAGNOSTIC.md` finding that the loss is a reconstruction/alignment effect separated by bucket contrast `C_u`, not a dead-time-saturation effect (its per-pattern load test was unpopulated in the frozen shards; this probe supplies that test and it comes back negative).

See `QUANTILE_RESCUE_PROBE.png` (per-scene dQ bars + mechanism scatter) and `QUANTILE_RESCUE_PROBE.json` (per-seed values, load quantiles).

## 5. Honest summary (exploratory, n small, single-knob)

Calibrating the same single global RIDGE multiplier to a load quantile (q75 or q90) instead of the mean successfully pulls the per-pattern load distribution back under the finite-window ridge, but it does **not** rescue any of the five failing scenes: the largest single-knob gain is +1.26 dB (`contour_0`, still net-negative) and four of five negatives move by under 0.06 dB, while both positive controls are preserved. The rescue is weak and non-monotone in the past-ridge-fraction reduction (Pearson 0.41 / Spearman 0.17), and even the non-deployable per-pattern-servo upper reference fails to flip the sign, so the contour/spokes deficit is not primarily a past-ridge dead-time effect and the cheapest single-knob quantile recalibration does not fix it. This is an exploratory post-hoc probe (n = 5 negatives + 2 controls, single-knob variants, nu=2000 only) and changes no frozen result or verdict.
