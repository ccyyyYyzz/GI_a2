# Speed-ratio decomposition of the corrected elapsed-time speedup

**POST-HOC DESCRIPTIVE.** Not a preregistered endpoint, carries no PASS/FAIL verdict.
Reuses only the frozen `m1-freeze @ 6f00932` campaign CSVs and the frozen kernel
table; **zero new simulation**. All figures are reproducible from
`results/round63_next/SPEED_DECOMPOSITION.json`.

## Inputs (all frozen)

| role | path |
|---|---|
| speed verdict of record | `results/round63_m1/CORRECTION_2026-07-19/M1_VERDICTS_SPEC_CORRECTED_R19.json` |
| imaging shards | `results/round63_m1/shards/M1_{SCAT32-SAFE,RIDGE-SCAT32,SCAT32-060}_00.csv` |
| finite-nu kernel | `results/round63_theory/fisher_ridge.csv` (`I_exact_log`=J_ν(ρ), `I_clt_log`=f(ρ)=ρ/(1+ρ), `ratio`=J/f) |
| ridge requested load ρ_R(ν) | ridge design npz / `ridge_targets.json` (global per dwell) |

Definitions: `detected_rate = mean_counts / optical_time_s` (per row, campaign
normalization); `T_opt = M·ν·τ` (M=1024, τ=50 ns); quality = `PSNR_rad`; Q90 target
= `safe_pava[0] + 0.9·(safe_pava[-1] − safe_pava[0])`; `S_gate = T_safe_cross /
T_ridge_cross` on the frozen elapsed-optical-time axis. Arm is read from the
`shard_id` prefix, never the CSV `arm` column (which is `RQL` for every row).

## The identity

For every NORMAL scene, with detected photons to the crossing `N = rate · T`:

```
S_gate  =  T_safe_cross / T_ridge_cross
        =  ( rate_ridge_cross / rate_safe_cross )  ×  ( N_safe_cross / N_ridge_cross )
        =  detected-rate ratio (at the crossing)   ×  residual factor
```

The residual factor `N_safe_cross / N_ridge_cross` is the ratio of detected photons
each arm needs to reach the **same** Q90 target — it absorbs per-detected-photon
informativeness and estimator/censoring efficiency. Crossings and all
interpolations are linear in `log(T_opt)` with the identical bracket and fraction
the frozen quality crossing used, so the decomposition is exact per scene
(`rate_ratio × residual − S_gate = 0` to machine precision).

**Validation.** The recomputed `S_gate` reproduces the frozen `per_image_s_gate`
for all 21 NORMAL scenes (max abs diff 4.6e-7); the recomputed all-24 median is
19.127043, matching the frozen `median_S_gate` bit-for-bit.

## Step 1 — Theoretical anchor vs empirical detected-rate ratio

Detected fraction f(ρ)=ρ/(1+ρ). Operating loads: SAFE ρ=0.05 (fixed);
RIDGE terminal requested ρ_R(2000)=22.2545 (achieved ≈22.2545 at ν=2000, no clip).

| quantity | value |
|---|---|
| f(0.05) | 0.047619 |
| f(22.2545) | 0.956989 |
| **theoretical terminal detected-rate ratio** f(22.25)/f(0.05) | **20.097** |
| **empirical terminal ratio** median rate_R/rate_S at ν=2000 | **20.085** (Q1–Q3 20.03–20.13) |

The empirical terminal detected-rate ratio matches the theoretical f-ratio to
0.06 %. This is the "detected-photon-rate ratio between the two operating points"
that the median S_gate (19.13) sits near. Across dwells, `mean_counts` tracks
ν·f(ρ_bar) to within 6 % everywhere (worst case the maze scenes, slightly
sub-Poisson), so `detected_rate` faithfully equals f(ρ)/(M·τ) — no `mean_counts`
anomaly survives a systematic |Δ|/pred>0.25 screen.

Per-dwell empirical vs theoretical detected-rate ratio (median across 24 scenes):

| ν | empirical rate_R/rate_S | theory f(ρ_ach)/f(0.05) |
|---|---|---|
| 5 | 4.18 | 4.16 |
| 10 | 6.86 | 6.90 |
| 20 | 10.29 | 10.35 |
| 50 | 14.67 | 14.85 |
| 100 | 17.30 | 17.38 |
| 200 | 18.95 | 19.01 |
| 500 | 19.54 | 19.58 |
| 1000 | 19.87 | 19.87 |
| 2000 | 20.09 | 20.10 |

(exact per-dwell/per-scene numbers in the JSON). The empirical rate ratio ramps
with ν because the ridge schedule ramps ρ_ach from ≈0.24 (ν=5) to 22.25 (ν=2000);
it saturates at the terminal 20.1 only at the top of the ladder.

## Step 2 — Per-scene decomposition

| statistic | rate-ratio at crossing | residual factor |
|---|---|---|
| median | 16.271 | 1.253 |
| Q1 / Q3 | 15.960 / 16.673 | 1.119 / 1.323 |
| min / max | 12.852 / 19.650 | 0.150 / 3.591 |

Residual factor median by family:

| family | median residual | n |
|---|---|---|
| chirp | 1.158 | 4 |
| contour | 0.679 | 2 (contour_0, contour_3) |
| glyph | 1.197 | 4 |
| maze | 1.761 | 4 |
| microtexture | 1.046 | 3 |
| spokes | 1.340 | 4 |

The maze family carries the largest residual (~1.76): its ridge reconstructions
reach Q90 with roughly half the detected photons SAFE needs — geometric/edge
structure is where the high-load ridge policy is most photon-efficient. The two
NORMAL contour scenes are the opposite: `contour_0` has residual 0.15 (ridge needs
~7× more photons, crosses late at ν≈543, unclipped), the low-contrast regime where
the ridge advantage nearly vanishes (S=2.94). `spokes_0` is the high tail
(residual 3.59, S=46), crossing earliest at bracket (20,50).

## Step 3 — Clip accounting and attribution of 20.1 → 19.13

The ridge trust-clip is active at ν∈{5,10,20,50,100} (achieved==requested exactly
for ν≥200). **20 of 21** NORMAL ridge Q90 crossings land at a clipped dwell — 19 in
bracket (50,100), one in (20,50); only `contour_0` crosses unclipped at (500,1000).
At the crossing the ridge runs at ρ_ach median ≈3.8 (clip ratio ρ_ach/ρ_req median
≈0.55), so its detected-rate ratio **there** is only 16.27, i.e. 0.81× the terminal
20.1.

Finite-ν Fisher-information cost of the clip, J(ρ_ach)/J(ρ_req), at the clipped
dwells that have kernel rows (ν=5,10 are absent from `fisher_ridge.csv`):

| ν | clip ratio (median) | J(ρ_ach)/J(ρ_req) median | min |
|---|---|---|---|
| 20 | 0.490 | 0.754 | 0.658 |
| 50 | 0.491 | 0.884 | 0.819 |
| 100 | 0.610 | 0.969 | 0.931 |

Because J(ρ) saturates near the ridge, halving the load at ν=50–100 costs only
~3–12 % of the finite-ν Fisher information even though it halves the incident load.

**Attribution.** The clip does *not* produce a net speed shortfall on NORMAL scenes:
it drops the detected-rate ratio to 16.27 (×0.81 of terminal), but the residual
factor (×1.253, ridge reaches Q90 with ~25 % fewer detected photons) very nearly
cancels it, leaving the NORMAL-scene median S = **20.384** — slightly *above* the
terminal 20.1. The move to the reported all-24 median **19.127** is entirely the
censoring step: the three hard scenes `m1_contour_1`→0, `m1_contour_2`→0
(FAST_RIGHT_CENSORED) and `m1_microtexture_1`→1 (SAFE_RANGE_UNINFORMATIVE) pull the
median down from 20.38 (21 NORMAL) to 19.13 (24). So of the 20.1→19.13 gap, the
clip explains ≈0 (offset by photon efficiency) and the censoring taxonomy explains
essentially all of it.

## Step 4 — The c=0 signature (information per detected photon)

At ν=2000, J_ν(ρ)/f(ρ) — finite-window Fisher information normalized by detected
fraction — is ≈1 at low load and falls as counting starts to outrun learning:

| ρ | J_ν(ρ) | f(ρ) | J/f |
|---|---|---|---|
| 0.05 | 0.04772 | 0.047619 | 1.00001 |
| 0.60 | 0.37525 | 0.375000 | 1.00005 |
| 5.00 | 0.83215 | 0.833333 | 0.99865 |
| 11.00 | 0.91082 | 0.916667 | 0.99389 |
| 22.25 | 0.93512 | 0.956989 | 0.97726 |
| 40.00 | 0.90981 | 0.975610 | 0.93267 |

J/f first drops below 0.99 at ρ≈14.4, below 0.95 at ρ≈34.0, below 0.90 at ρ≈50.0
(ν=2000). Operating points: SAFE (ρ=0.05) and SCAT32-060 (ρ=0.60) sit at J/f≈1.000
— every detected photon carries essentially its full CLT information; the RIDGE
terminal (ρ=22.25) already sits at J/f=0.977, i.e. ~2.3 % of each detected photon's
count no longer converts to learning. This is the finite-window ridge: "where
counting stops equaling learning." It also explains the sign of the residual/rate
trade-off — the ridge crossing lives at ρ_ach≈3.8 (J/f≈0.999), *below* the terminal
ridge, so at the crossing the ridge is still in the near-lossless per-photon regime,
which is consistent with the residual factor exceeding 1.

## Step 5 — Summary (for adaptation into paper-2 discussion; the paper is not edited here)

1. The corrected elapsed-time speedup factorizes exactly, per scene, as
   `S_gate = (detected-rate ratio at the Q90 crossing) × (detected photons to reach
   the same target, safe/ridge)`, validated against the frozen `per_image_s_gate` to
   4.6e-7.
2. The empirical detected-photon-rate ratio between the two terminal operating
   points is 20.09, matching the closed-form f(22.25)/f(0.05)=20.10 to 0.06 %, and
   this rate ratio is the quantity the median speedup sits near.
3. However S_gate is not that terminal ratio mechanically: nearly every ridge
   crossing lands at a clipped mid-ladder dwell where the detected-rate ratio is only
   16.3, and a residual factor of 1.25 (the ridge reaches Q90 with ~25 % fewer
   detected photons) restores it, giving a NORMAL-scene median of 20.4.
4. The trust clip, whose finite-ν Fisher-information cost is only ~3–12 % at the
   relevant dwells because J(ρ) saturates near the ridge, therefore does not create a
   net speed penalty — it is offset by per-photon efficiency — and the reported
   median (19.13) differs from the NORMAL median (20.38) purely through censoring of
   three low-contrast/degenerate scenes.
5. The information-per-detected-photon curve J_ν/f stays at ≈1 up to ρ≈14 and only
   then declines (0.977 at the ρ=22.25 terminal ridge), so the speedup is bought in
   the near-lossless counting regime, with the finite-window loss beginning exactly
   where the ridge peaks.

## Caveats

- POST-HOC DESCRIPTIVE; n=24 scenes (21 NORMAL after the frozen censoring taxonomy);
  no verdict, no threshold, no new data.
- Crossing and rate/photon interpolation are linear in `log(T_opt)`; a different
  interpolation of the rate curve would perturb the rate-ratio/residual split
  (their product S_gate is interpolation-invariant and matches the frozen value).
- `fisher_ridge.csv` lacks rows for ν=5,10, so the clip's Fisher-info cost is
  quantified only at ν∈{20,50,100}; the two lowest clipped dwells are reported by
  clip ratio only.
- Requested ρ_R(ν) is the global ridge schedule; achieved ρ is the per-scene,
  per-seed `rho_bar` averaged over the 5 seeds, consistent with how the frozen
  curves are built.
- The near-equality 16.3 × 1.25 ≈ 20.1 (crossing rate ratio × residual ≈ terminal
  rate ratio) is an emergent cancellation of two independent mechanisms (clip/schedule
  vs photon efficiency), not an identity; it should be described as such.
