# JET_TEST — decisive falsification of the R42 moonshot

**Quotient Information Jets / Curvature-Rescued Detection** — R42 PRL candidate.
Executes `docs/ROUND63_GPT_ROUND42_RULING_RAW.md` §6.5 verbatim (Banks A/B/C + kill conditions), on the exact `SCRAMBLE_EXT` Kronecker channel `Σ(x)=R+Q(x)H`, `Q(x)=xᵀGx`.

## VERDICT: **MOONSHOT_SURVIVES**

Runtime 59 s. All scored checks:

| check | pass |
|---|---|
| kl_order_generic_is_2 | ✅ |
| kl_order_ortho_is_4 | ✅ |
| orders_integer | ✅ |
| coef_generic_within_10pct | ✅ |
| coef_ortho_within_10pct | ✅ |
| crossover_within_25pct | ✅ |
| delay_slope_generic_near_-2 | ✅ |
| delay_slope_ortho_near_-4 | ✅ |
| mc_generic_m_near_1 | ✅ |
| mc_ortho_m_near_2 | ✅ |
| amplitude_collapse_to_zero | ✅ |
| collapse_persists_lags | ✅ |
| info_returns_with_anchor | ✅ |
| auc_returns_to_half_at_star | ✅ |
| local_visible_either_side | ✅ |
| cone_dQ_all_positive | ✅ |
| no_broad_blind_set | ✅ |

Kill conditions (any true ⇒ killed):

| kill condition | triggered |
|---|---|
| non_integer_contact_order | no |
| eps4_class_vanished_not_expensive | no |
| info_survives_collinear_amplitude | no |
| broad_in_aperture_blind_set | no |

## Bank A — fixed-law jet slopes

`||K||_F² = 0.000346`.  Directions (norm = ‖x‖):

| dir | c = xᵀGδ | d = δᵀGδ | ε_cross=2|c|/d |
|---|---|---|---|
| delta_g | +29.78 | 164.3 | 0.3626 |
| delta_o | +1.998e-15 | 137.4 | 2.908e-17 |
| delta_n | -3.066e-06 | 44.19 | 1.387e-07 |
| delta_mix | +0.6872 | 137.5 | 0.009995 |

| dir | KL slope (pred %s) | Bhatt slope | Chernoff* slope | coef ratio (meas/pred) |
|---|---|---|---|---|
| delta_g | **2.038** (2) | 2.037 | 2.037 | 1.0107 |
| delta_o | **4.000** (4) | 4.000 | 4.000 | 1.0000 |
| delta_n | **4.000** (4) | 4.000 | 4.000 | 0.9999 |

- **delta_g** (generic, c≠0): exact-KL contact order 2, coefficient `c²‖K‖²` matched to 1.07%.
- **delta_o** (exactly G-orthogonal, c=0): exact-KL contact order 4, coefficient `d²‖K‖²/4` matched to 0.00%. This is curvature-rescued detection: the score vanishes, the Hessian 2G>0 lifts the law at 2nd order.
- **delta_n** (concentrated on G's smallest eigenvalues): d>0 so jet order stays ≤2 (slope 4.00); no local blindness inside the aperture.

**Mixed-direction crossover (eq 6.23):** predicted ε_cross = 0.009995, empirical (local KL slope→3) = 0.009997, ratio 1.000.

## Bank A — Monte-Carlo detector (real Kronecker generator, finite-sample)

Paired matched-score duel; d′ ∝ εᵐ√T at fixed T. Slope of log d′ vs log ε gives m:

| dir | m (=slope) | pred | rows (ε: d′@T) |
|---|---|---|---|
| delta_g | **0.95** | 1 | 0.050:3.35, 0.035:1.76, 0.025:1.73, 0.018:1.14 |
| delta_o | **2.05** | 2 | 0.160:3.78, 0.120:2.65, 0.090:1.41, 0.060:0.53 |

The ε⁻⁴ class is **expensive, not vanishing**: delta_o d′ still grows as ε²√T on real generator banks (finite-sample CRB effect, not a floor).

## Bank A — CUSUM detection delay (sequential)

Fixed false-alarm (ARL₀≈500). Delay ∝ 1/D_*(ε) ⇒ slope −2m.

| dir | delay slope | pred | (ε: delay) |
|---|---|---|---|
| delta_g | **-2.16** | -2 | 0.060:1435, 0.030:6407 |
| delta_o | **-3.92** | -4 | 0.160:906, 0.090:8623 |

## Bank B — nuisance quotient (unknown medium amplitude)

Efficient info J_θ (nuisance-profiled) as a fraction of the naive I_θ:

| configuration | J_θ / I_θ (naive) | verdict |
|---|---|---|
| single zero-lag, unknown a | 0.000e+00 | **collapse → 0** |
| + proportional lags (known shape) | 0.000e+00 | **collapse persists** |
| + 2nd wavelength (kernel disk_r4) | 0.018 | **info returns** |
| + amplitude anchor (finite prior) | 0.200 | **info returns** |

Zero-lag scene/amplitude scores are exactly collinear (both ∝ H), so the profiled information annihilates the whole scene channel (J_θ/I_θ = 0.0e+00 ≈ numerical zero). The fractional covariance response c/Q differs between wavelengths (0.1589 vs 0.1114, gap 0.04751), which is exactly what restores identifiability.

## Bank C — global cancellation and monotone cone

**Signed cancellation:** ε_* = −2c/d = 0.3626. At ε_* the covariance is identical (z=0.00e+00), so AUC = **0.481** (chance) — despite clear local visibility either side (AUC 0.026 at 0.25ε_*, 1.000 at 1.8ε_*; |AUC−0.5| = 0.474 / 0.500). Detection is a *local* quotient statement; iso-Q surfaces are globally invisible.

**Monotone cone:** grain kernel G entrywise min = 2.890e-09 (nonneg=True); with x,δ≥0, ΔQ>0 at every tested amplitude (True). The sign-definite corollary holds globally on the cone.

## Kill-condition scan

- Contact orders integer (2, 4): PASS.
- No broad in-aperture blind set: scanned 400 beyond-band directions, min d = 2.544e-01 (median 2.648e-01), fraction blind = 0.0e+00.
- Collinear-amplitude residual info: 0.0e+00 (≈0).
- ε⁻⁴ class expensive-not-vanishing: MC delta_o m=2.05 (>0).

Figures: `jet_slopes.{png,pdf}` (THE PRL log-log figure), `bank_b_nuisance.{png,pdf}`.
