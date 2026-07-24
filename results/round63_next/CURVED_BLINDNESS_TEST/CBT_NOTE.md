# CBT_NOTE — CURVED_BLINDNESS_TEST (R46 moonshot: the Curvature Tax of Undetectability)

Preregistered exam of GPT's R46 moonshot, executing `docs/ROUND63_GPT_ROUND46_RULING_RAW.md` §2.10 verbatim (frozen arms / predictions / kill conditions). Exact scrambled Gaussian Q-engine reused from `SCRAMBLE_EXT` + `JET_TEST`. This is preregistration — no repair round.

## VERDICT: **MOONSHOT_SURVIVES**

Runtime 7079 s, config: FULL. Geometry (DC-free beyond-band subspace, so the mean/DC channel is exactly null for every arm and Q is the sole sufficient statistic): r0=0.4272, s=0.2991, A*=s²/r0=0.2093, ω=s/r0=0.7000; tangency residual ⟨x,v⟩_G/(r0·s)=6.8e-18; background/texture cross-null=3.1e-16; operating Q0=11.56, I_Q=0.04492.

### Per-prediction table (measured vs frozen)

| # | frozen prediction | measured | pass |
|---|---|---|---|
| 1 | STRAIGHT KL/Chernoff slope = 4.00±0.05 | KL slope 4.0000 (Chernoff 3.9928) | PASS |
| 2 | PARTIAL α=0.25 slope = 4.00±0.05 | KL slope 4.0000 | PASS |
| 3 | C_α/C_0 = (1-α)²=0.5625 within 5% | 0.5625 (err 0.00%) | PASS |
| 4 | PARTIAL α=0.50 slope = 4.00±0.05 | KL slope 4.0003 | PASS |
| 5 | C_α/C_0 = (1-α)²=0.2500 within 5% | 0.2500 (err 0.01%) | PASS |
| 6 | PARTIAL α=0.75 slope = 4.00±0.05 | KL slope 4.0013 | PASS |
| 7 | C_α/C_0 = (1-α)²=0.0625 within 5% | 0.0625 (err 0.03%) | PASS |
| 8 | EXACT CLOAK |ΔQ|/Q < 1e-12 | 1.54e-16 | PASS |
| 9 | EXACT CLOAK divergence below floor (<1e-18) | max KL 2.47e-32 | PASS |
| 10 | EXACT CLOAK MC paired d' consistent with 0 (2σ) | all 6 cells d'-null | PASS |
| 11 | EXACT CLOAK unpaired AUC 95% CI contains 0.5 | all cells | PASS |
| 12 | EXACT CLOAK paired AUC ∈ [0.49,0.51] (soft, MC-resolution-limited) | all cells | PASS |
| 13 | MC positive controls (STRAIGHT, PARTIAL) detectable | straight AUC 0.997, partial AUC 0.934 | PASS |
| 14 | RADIAL CONTROL slope = 2.00±0.05 | KL slope 2.0274 | PASS |
| 15 | Acc-cap kink at A* within 5% | A_kink=0.2093 vs A*=0.2093 (err 0.00%) | PASS |
| 16 | Acc-cap coef follows [1-A/A*]_+² (max abs err <5%) | max abs err 0.0000 | PASS |
| 17 | CUSUM partial-cloak delay slope ≈ -4 | block-delay slope -3.947 (block=303645, jet-dev 7.1e-03) | PASS |
| 18 | CUSUM exact cloak never triggers above FA rate | detect_frac 0.0000 (FA 0.0020) | PASS |
| 19 | 2nd model (quartic h=(xᵀG_b x)²) STRAIGHT slope 4 | gauss 4.000 / pois 4.006 | PASS |
| 20 | 2nd model RADIAL slope 2 | gauss 2.096 | PASS |
| 21 | 2nd model coefficient law (1-α)² within 10% | α ratios 0.562, 0.250, 0.063 | PASS |
| 22 | 2nd model exact cloak null (<1e-12) | max KL 4.93e-32 | PASS |

Predictions passed: **22/22**.

### Kill-condition status

| kill condition | triggered |
|---|---|
| great-circle record distinguishable despite fixed Q | no |
| coefficient law misses (1-α)² by >10% | no |
| acceleration threshold misses A* by >10% | no |
| nuisance profiling introduces a lower-order term | no |
| not reproduced in a second non-optical model | no |

### Wrong-metric control

Designing the cloak in the Euclidean norm instead of the declared actuator metric G still cloaks to 2nd order (residual -2.8e-17) but OVER-PAYS the true G-cost by a factor 1.0336, and its acceleration-cap kink lands at 0.2198 instead of A*=0.2093 — the curvature tax A*=s²/r0 is metric-specific.

### Interpretation

The blind spot is a curved road, not a direction. A state on the G-ellipsoid Q=const is perfectly hidden (infinite-data image privacy), yet every straight local event is curvature-visible at fourth jet order; only the great-circle orbit — which pays the exact centripetal tax A*=s²/r0 — stays exactly invisible. Any acceleration shortfall re-appears as a universal fourth-order Chernoff leakage whose coefficient is the SQUARE of the unpaid curvature tax, (1-α)² / [1-A/A*]_+². Verified analytically to machine precision on the optical scrambled channel AND reproduced in a non-optical heteroscedastic-Gaussian / Poisson model with a genuinely different (quartic) curved invariant. Figures: `cbt_slopes_coef_kink.png` (per-arm log-log slopes, coefficient-vs-α, leakage-vs-A kink), `cbt_cusum_model2.png` (CUSUM t⁻⁴ delay + non-optical replication).
