# Contour-family diagnostic — why the primary ridge gain goes negative on `contour`

**Status: post-unblinding descriptive diagnostic (M1, tag `m1-freeze` @ 6f00932).**
Source of record: `results/round63_m1/shards/M1_{RIDGE-SCAT32,SCAT32-060}_00.csv`
(imaging arms; arm identity = `shard_id` prefix, estimator = `RQL`).
Aggregation: seed-mean radiometric PSNR (`PSNR_rad`) per image at `nu=2000`,
then per-family; `dQ = PSNR_rad(RIDGE-SCAT32) - PSNR_rad(SCAT32-060)`.
Primary condition only (`nu=2000`), matching the frozen primary endpoint.

## 1. The fact to explain

Per-family median `dQ` at `nu=2000` (matches `M1_VERDICTS.json`):

| family | median dQ (dB) | sign |
|---|---|---|
| maze | +6.369 | + |
| glyph | +3.582 | + |
| microtexture | +1.487 | + |
| chirp | +1.405 | + |
| spokes | +0.715 | + |
| **contour** | **−4.175** | **−** |

`contour` is the only family with a negative median. It is also the **only**
family whose mean `SCAT32-060` radiometric PSNR **exceeds** the ridge arm's:

| family | RIDGE mean PSNR_rad | 060 mean PSNR_rad | RIDGE − 060 |
|---|---|---|---|
| chirp | 10.97 | 9.42 | +1.54 |
| contour | **13.47** | **15.97** | **−2.50** |
| glyph | 24.31 | 20.46 | +3.85 |
| maze | 32.16 | 25.97 | +6.19 |
| microtexture | 18.58 | 17.14 | +1.44 |
| spokes | 13.99 | 13.77 | +0.22 |

## 2. Per-instance contour detail (`nu=2000`, seed-mean)

| instance | C_u | RIDGE PSNR_rad | 060 PSNR_rad | dQ (dB) | RIDGE rad_nrmse | 060 rad_nrmse | Γ(RIDGE) |
|---|---|---|---|---|---|---|---|
| contour_0 | 0.213 | 8.857 | 12.469 | −3.61 | 0.517 | 0.337 | 9.31 |
| contour_1 | 0.217 | 7.765 | 13.559 | −5.79 | 0.562 | 0.288 | 9.51 |
| contour_2 | 0.245 | 10.008 | 14.746 | −4.74 | 0.509 | 0.295 | 10.71 |
| contour_3 | 0.331 | 27.259 | 23.107 | +4.15 | 0.099 | 0.160 | 14.50 |

Three of four contour instances lose under the ridge; the fourth (contour_3)
gains as strongly as the other families. The ridge policy hits the same
achieved mean load on every scene (`rho_bar = 22.2545 = rho_R(2000)`; the
global safety clip does not bind), so the family split is **not** a load or
clip effect.

## 3. What is NOT the cause (no instrument/guard artifact)

For every contour cell in both arms, at `nu=2000`:

- `dark_frac = 0`, `tau_err = 0` (no dark counts, no timing-estimator error);
- predicted count-ceiling fraction `p_ceil ≈ 0` at all operating loads
  (frozen kernel `oed_design_v3.kernel_eval`; the ridge load 22.25 < cap 24,
  and `P(N=nu)` underflows to 0);
- no `RIDGE_GUARD_CLIPPED` event (achieved mean load equals the ridge target,
  so the global clip did not bind);
- in the confirmatory certificate shards `MU_CAP_ACTIVE = False` and
  `dep_feasible = True` for the contour cells;
- the trust/audit columns `audit_status`, `leak_suspect`, `cnr`, `d_ratio`,
  `q_d`, `q_mean` are **unpopulated** in the frozen shards (the once-per-cell
  descriptive audit wrote no value), so no trust or leakage flag is recorded
  as firing on these cells.

The contour deficit is therefore a genuine reconstruction-quality effect, not
a saturation, clipping, dark-count, or trust-flag artifact.

## 4. The covariate that separates win from loss: bucket contrast C_u

Within contour, `dQ` is **monotone** in the bucket-brightness contrast `C_u`
(a per-image quantity, identical across arms): the three low-contrast
instances (C_u 0.213–0.245) all lose; the one high-contrast instance
(C_u 0.331) wins. Across all 24 confirmatory scenes, `dQ` correlates with
`C_u` at Pearson **+0.639**, and the losing scenes are concentrated at low
C_u (the two other negative-dQ scenes outside contour are spokes_1 C_u 0.297
dQ −5.08 and spokes_0 C_u 0.296 dQ −0.77).

This is consistent with the manuscript's alignment condition (main §3.6 /
supplement S2): a single global flux multiplier performs implicit information
water-filling only when bucket brightness `u_i` is aligned with directional
value `ell_i`. Low-contrast scenes are exactly where that alignment is
weakest, so the ridge over-drives the scene without radiometric benefit — the
ridge-arm reconstruction carries a larger radiometric NRMSE (0.51–0.56) than
the balanced 0.60 design (0.29–0.34) on the three losing instances, and the
ordering reverses on the high-contrast instance (0.099 vs 0.160).

## 5. What C_u does NOT fully explain (hypothesis, flagged)

Contrast alone is not the whole story: at essentially matched low contrast,
chirp_0 (C_u 0.215) **gains** +2.36 dB while contour_0 (C_u 0.213) **loses**
−3.61 dB. A family-geometry factor therefore compounds the low-contrast
misalignment. The plausible mechanism — **hypothesis, not established by
these data** — is that thin, sparse contour supports concentrate the
informative signal on a few edge pixels whose radiometric scale the high-flux
ridge reconstruction degrades, whereas broadband families (chirp) tolerate
the extra flux. The per-pattern load-quantile columns that would test the
per-support over-drive directly (`q_d`, `q_mean`) are unpopulated in the
frozen shards, so this remains a hypothesis.

## 6. Most defensible one-sentence characterization (for the discussion)

The single negative-median family is `contour` (median dQ −4.18 dB), where the
balanced 0.60 design outperforms the ridge on three of four instances; the
separating covariate is bucket-brightness contrast (losing instances C_u
0.21–0.25 vs the winning instance's 0.33; all-scene dQ–C_u correlation +0.64),
consistent with the alignment condition under which a single global flux
multiplier is least effective when bucket brightness is weakly aligned with
directional value, with no guard, trust, dark-count, or ceiling flag firing on
these cells.
