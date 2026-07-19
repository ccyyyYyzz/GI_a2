# Dose-only headroom — DEV constructive evidence (R18 §1.2/§6, frozen set)

**Development-only constructive analysis; not a confirmatory endpoint or
population estimate.**

Frozen replication set (declared before reading aggregates, R18 §6): all six
DEV family images (632900+ instances), seed 0, the two anchor cells
(nu=2000, b=0.60) and (nu=200, b=0.05). Every result reported without
selection. Solver: dose-only primal probe over D_cert within C_dose
(incident budget + ±5% per-pixel dose), multiplicative ascent + projection.

## Source-of-record fields (per cell)

DOSE_ONLY_SOLVER_STATUS = FEASIBLE_CONSTRUCTION (all 12 cells);
DOSE_ONLY_DUAL_UPPER_IF_AVAILABLE: glyph fast 26.7 (loose), glyph safe 79.1
(loose) / 106.2 (tightened, CG residual 3.5e4 unconverged); other cells: not
computed.

| scene | nu | b | SCAT32 logdet | probe logdet | DOSE_ONLY_PRIMAL_GAP_PER_R | DOSE_ONLY_D_EFFICIENCY_LOWER | D_gain mass | load q50 | budget used | dose dev |
|---|---|---|---|---|---|---|---|---|---|---|
| m1_dev_glyph | 2000 | 0.60 | 3646.425 | 3801.577 | 0.7758 | 2.17x | 0.000 | 0.300 | 0.3000 | 0.038 |
| m1_dev_glyph | 200 | 0.05 | 2788.834 | 2969.346 | 0.9026 | 2.47x | 0.000 | 0.025 | 0.0250 | 0.038 |
| m1_dev_chirp | 2000 | 0.60 | 3646.994 | 3819.047 | 0.8603 | 2.36x | 0.000 | 0.300 | 0.3000 | 0.038 |
| m1_dev_chirp | 200 | 0.05 | 2788.149 | 2887.359 | 0.4960 | 1.64x | 0.000 | 0.025 | 0.0250 | 0.038 |
| m1_dev_maze | 2000 | 0.60 | 3633.386 | 3750.598 | 0.5861 | 1.80x | 0.000 | 0.300 | 0.3000 | 0.038 |
| m1_dev_maze | 200 | 0.05 | 2790.259 | 2892.347 | 0.5104 | 1.67x | 0.000 | 0.025 | 0.0250 | 0.038 |
| m1_dev_spokes | 2000 | 0.60 | 3640.443 | 3820.941 | 0.9025 | 2.47x | 0.000 | 0.300 | 0.3000 | 0.038 |
| m1_dev_spokes | 200 | 0.05 | 2800.765 | 2958.083 | 0.7866 | 2.20x | 0.000 | 0.025 | 0.0250 | 0.038 |
| m1_dev_contour | 2000 | 0.60 | 3648.639 | 3829.357 | 0.9036 | 2.47x | 0.000 | 0.300 | 0.3000 | 0.038 |
| m1_dev_contour | 200 | 0.05 | 2772.045 | 2891.034 | 0.5949 | 1.81x | 0.000 | 0.025 | 0.0250 | 0.038 |
| m1_dev_microtexture | 2000 | 0.60 | 3667.825 | 3851.978 | 0.9208 | 2.51x | 0.000 | 0.300 | 0.3000 | 0.038 |
| m1_dev_microtexture | 200 | 0.05 | 2810.246 | 2953.895 | 0.7182 | 2.05x | 0.000 | 0.025 | 0.0250 | 0.038 |

Gap/r range: 0.496 .. 0.921 across all six families and both anchors. All
twelve constructions are budget- and dose-feasible (verified at every
recorded point).

## Mandatory caveat (R18 §1.2, recorded verbatim)

The probes used for this table initialize positive mass ONLY on the lowest
D_load column; the multiplicative update `xi = xi * sqrt(score/mean_score)`
is support-preserving, so a zero column can never become positive. Therefore
the following are valid: a geometry-only, lowest-load feasible construction
exists; that construction already gives the reported lower-bound headroom;
gain coupling is not necessary to refute the dose-only certificate. The
following are NOT established: that the dose-only optimum assigns zero mass
to D_gain; that the optimum prefers the lowest palette load; that the
optimum intentionally spends only half the photon budget. These numbers are
lower bounds on available headroom (constructive existence), not estimates
of the dose-only optimum, and no DOSE_ONLY_*_PASS verdict exists.

Approved future-work sentence (R18 §6): "The dose-only construction shows
that uniform exposure still permits scene-adapted geometry selection.
Optimizing and validating that geometry under relaxed or task-specific
conditioning constraints is left to a separate follow-up; no such arm or
endpoint is added to the present campaign."
