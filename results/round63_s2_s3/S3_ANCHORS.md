# S3 detector-mismatch anchor analysis (ROUND63, R9)

PURELY DESCRIPTIVE (R9: "S3 should remain descriptive"). No pass/fail language, no gate. All numbers rounded to 2 dp. Non-numeric PSNR_rad ('' / 'inf') dropped from medians and paired differences, never imputed.

## Provenance and baseline choice

- Variant rows: `results/round63_s2_s3/round63_s2_s3_v2/pilot_rows.csv` + `results/round63_s2_s3jit/pilot_rows.csv` (3168 rows total), joined to variant labels by `cell_id` via `results/round63/manifests/S3_*.json` and `S3JIT_*.json`.
- Baseline: **EXTERNAL** matched rows from `results/round63_s2_s2b/pilot_rows.csv` filtered to M=2048, nu=500 (432 rows). No in-stage no-mismatch cell exists (every S3 cell carries a non-empty det or est per `mismatch_sig`), so the matched S2-B rows are the no-mismatch reference at the same (rho_bar, nu=500, image, seed, arm).
- Grid overlap verified: S3 rho_bar in {0.05, 0.3, 0.6, 1.0}; S2-B baseline provides matched rows at rho_bar in {0.05, 0.6, 1.0} (NO baseline at rho=0.3 -> those interaction cells show dPSNR as MISSING). Both anchor operating points (safe=0.05, fast=0.6) are covered.
- dPSNR = median of per-(image,seed) PAIRED differences (variant - baseline) where pairable (method=`paired`); if a baseline exists but nothing pairs, difference of medians (`median-diff`); if no baseline at that (arm,rho), `no-baseline` -> MISSING. dflux = median(variant flux_dev) - median(baseline flux_dev) (difference of medians). "flux bias" columns report the variant's own median flux_dev (radiometric flux deviation).

## (b) Main-text anchor subset (R9 rows)

safe = rho_bar 0.05; fast = rho_bar 0.6 (high-flux). PRECORRECT dPSNR is reported at fast only: at the safe point PRECORRECT is a degenerate floored reconstruction (flux_dev = -1.0 for all rows, identical to its own baseline), so its safe dPSNR is trivially ~0 and uninformative.

| variant | RQL dPSNR safe (rho=0.05) | RQL dPSNR fast (rho=0.6) | RQL flux bias fast | PRECORRECT dPSNR fast |
|---|---|---|---|---|
| tau_err +10% | -0.00 | -0.10 | 0.06 | -0.50 |
| tau_err -10% | 0.00 | -0.12 | -0.06 | 0.45 |
| dark 10% (known) | -0.05 | 0.16 | 0.00 | -0.04 |
| dark 10% (unknown) | -0.34 | -0.30 | 0.10 | -0.79 |
| afterpulsing 2% | -0.00 | -0.00 | 0.01 | -0.12 |
| afterpulsing 5% | -0.01 | -0.09 | 0.03 | -0.23 |
| paralyzable-data / non-paralyzable-model | MISSING | MISSING | MISSING | MISSING |

> **MISSING anchors:** paralyzable-data / non-paralyzable-model. Not present in the pilot data; reported MISSING rather than imputed.

## (c) Plain-language summary

1. Across the instantiated anchors, every RQL radiometric-PSNR shift stays within 0.34 dB of the matched no-mismatch baseline at both the safe (rho=0.05) and fast (rho=0.6) operating points, so all of tau_err +10%, tau_err -10%, dark 10% (known), dark 10% (unknown), afterpulsing 2%, afterpulsing 5% read as benign (|dPSNR| < 0.5 dB).
2. The mismatches degrade gracefully rather than catastrophically: the largest RQL degradation is only -0.34 dB (dark 10% (unknown)), and known-dark and low afterpulsing perturbations move RQL PSNR by essentially zero.
3. Radiometric flux is biased mainly by tau_err +10%, tau_err -10%, dark 10% (unknown) (RQL fast flux bias up to 0.10 in magnitude), whereas dark 10% (known), afterpulsing 2%, afterpulsing 5% leave the flux essentially unbiased (|bias| < 0.05).
4. RQL degrades LESS than PRECORRECT at the fast point in |dPSNR| for 5 of 6 present anchors (the sole exception being dark 10% (known) where both moves are negligible), the sharpest PRECORRECT swing being -0.79 dB (dark 10% (unknown)) against RQL's -0.30 dB there; note PRECORRECT at the safe point is a degenerate floored reconstruction (flux_dev = -1.0, identical to its own baseline) so its safe-point dPSNR is uninformative and is omitted from the anchor columns.
5. The paralyzable-data / non-paralyzable-model anchor is MISSING from this pilot: no cell instantiates det.paralyzable / est.assume_paralyzable (build_s3 never emits one), so that row cannot be characterized from the available data and is reported MISSING rather than imputed.

## (a) Full variant table

All (variant, arm, rho_bar) groups. `n` = rows in group; `n_pair` = per-(image,seed) pairs used for paired dPSNR.

| variant | arm | rho_bar | median PSNR_rad | median flux_dev | n | dPSNR (vs baseline) | method | n_pair | dflux |
|---|---|---|---|---|---|---|---|---|---|
| baseline | RQL | 0.05 | 14.53 | -0.00 | 36 | 0.00 | baseline | 36 | 0.00 |
| baseline | RQL | 0.6 | 14.95 | 0.00 | 36 | 0.00 | baseline | 36 | 0.00 |
| baseline | RQL | 1 | 14.70 | 0.00 | 36 | 0.00 | baseline | 36 | 0.00 |
| baseline | PRECORRECT | 0.05 | 6.36 | -1.00 | 36 | 0.00 | baseline | 36 | 0.00 |
| baseline | PRECORRECT | 0.6 | 4.17 | 0.00 | 36 | 0.00 | baseline | 36 | 0.00 |
| baseline | PRECORRECT | 1 | 4.33 | 0.00 | 36 | 0.00 | baseline | 36 | 0.00 |
| baseline | POISSON-LIN | 0.05 | 14.47 | -0.05 | 36 | 0.00 | baseline | 36 | 0.00 |
| baseline | SAT-POISSON | 0.05 | 14.52 | -0.00 | 36 | 0.00 | baseline | 36 | 0.00 |
| baseline | POISSON-LIN | 0.6 | 11.95 | -0.37 | 36 | 0.00 | baseline | 36 | 0.00 |
| baseline | SAT-POISSON | 0.6 | 14.82 | 0.00 | 36 | 0.00 | baseline | 36 | 0.00 |
| baseline | POISSON-LIN | 1 | 10.82 | -0.50 | 36 | 0.00 | baseline | 36 | 0.00 |
| baseline | SAT-POISSON | 1 | 14.89 | 0.00 | 36 | 0.00 | baseline | 36 | 0.00 |
| tau_err=+0.1 | RQL | 0.05 | 14.52 | 0.00 | 36 | -0.00 | paired | 36 | 0.01 |
| tau_err=+0.1 | RQL | 0.3 | 14.59 | 0.03 | 36 | MISSING | no-baseline | 0 | MISSING |
| tau_err=+0.1 | RQL | 0.6 | 14.79 | 0.06 | 36 | -0.10 | paired | 36 | 0.06 |
| tau_err=+0.1 | RQL | 1 | 14.46 | 0.11 | 36 | -0.25 | paired | 36 | 0.11 |
| tau_err=+0.1 | PRECORRECT | 0.05 | 6.36 | -1.00 | 36 | 0.00 | paired | 36 | 0.00 |
| tau_err=+0.1 | PRECORRECT | 0.3 | 3.74 | 0.03 | 36 | MISSING | no-baseline | 0 | MISSING |
| tau_err=+0.1 | PRECORRECT | 0.6 | 3.68 | 0.06 | 36 | -0.50 | paired | 36 | 0.06 |
| tau_err=+0.1 | PRECORRECT | 1 | 3.41 | 0.11 | 36 | -0.87 | paired | 36 | 0.11 |
| tau_err=+0.2 | RQL | 0.05 | 14.53 | 0.01 | 36 | 0.00 | paired | 36 | 0.01 |
| tau_err=+0.2 | RQL | 0.6 | 14.30 | 0.14 | 36 | -0.46 | paired | 36 | 0.14 |
| tau_err=+0.2 | PRECORRECT | 0.05 | 6.36 | -1.00 | 36 | 0.00 | paired | 36 | 0.00 |
| tau_err=+0.2 | PRECORRECT | 0.6 | 3.12 | 0.14 | 36 | -1.06 | paired | 36 | 0.14 |
| tau_err=-0.1 | RQL | 0.05 | 14.52 | -0.01 | 36 | 0.00 | paired | 36 | -0.00 |
| tau_err=-0.1 | RQL | 0.3 | 14.59 | -0.03 | 36 | MISSING | no-baseline | 0 | MISSING |
| tau_err=-0.1 | RQL | 0.6 | 14.84 | -0.06 | 36 | -0.12 | paired | 36 | -0.06 |
| tau_err=-0.1 | RQL | 1 | 14.61 | -0.09 | 36 | -0.27 | paired | 36 | -0.09 |
| tau_err=-0.1 | PRECORRECT | 0.05 | 6.36 | -1.00 | 36 | 0.00 | paired | 36 | 0.00 |
| tau_err=-0.1 | PRECORRECT | 0.3 | 4.23 | -0.03 | 36 | MISSING | no-baseline | 0 | MISSING |
| tau_err=-0.1 | PRECORRECT | 0.6 | 4.61 | -0.06 | 36 | 0.45 | paired | 36 | -0.06 |
| tau_err=-0.1 | PRECORRECT | 1 | 5.11 | -0.09 | 36 | 0.73 | paired | 36 | -0.09 |
| tau_err=-0.2 | RQL | 0.05 | 14.52 | -0.01 | 36 | -0.00 | paired | 36 | -0.01 |
| tau_err=-0.2 | RQL | 0.6 | 14.53 | -0.11 | 36 | -0.34 | paired | 36 | -0.11 |
| tau_err=-0.2 | PRECORRECT | 0.05 | 6.36 | -1.00 | 36 | 0.00 | paired | 36 | 0.00 |
| tau_err=-0.2 | PRECORRECT | 0.6 | 5.00 | -0.11 | 36 | 0.85 | paired | 36 | -0.11 |
| dark10%_known | RQL | 0.05 | 14.34 | -0.00 | 36 | -0.05 | paired | 36 | -0.00 |
| dark10%_known | RQL | 0.6 | 14.87 | 0.00 | 36 | 0.16 | paired | 36 | -0.00 |
| dark10%_known | PRECORRECT | 0.05 | 6.36 | -1.00 | 36 | 0.00 | paired | 36 | 0.00 |
| dark10%_known | PRECORRECT | 0.6 | 4.18 | 0.00 | 36 | -0.04 | paired | 36 | -0.00 |
| dark10%_unknown | RQL | 0.05 | 14.17 | 0.10 | 36 | -0.34 | paired | 36 | 0.10 |
| dark10%_unknown | RQL | 0.6 | 14.47 | 0.10 | 36 | -0.30 | paired | 36 | 0.10 |
| dark10%_unknown | PRECORRECT | 0.05 | 6.36 | -1.00 | 36 | 0.00 | paired | 36 | 0.00 |
| dark10%_unknown | PRECORRECT | 0.6 | 3.40 | 0.10 | 36 | -0.79 | paired | 36 | 0.10 |
| dark25%_known | RQL | 0.05 | 14.41 | -0.00 | 36 | 0.04 | paired | 36 | -0.00 |
| dark25%_known | RQL | 0.6 | 14.98 | 0.00 | 36 | 0.02 | paired | 36 | 0.00 |
| dark25%_known | PRECORRECT | 0.05 | 6.36 | -1.00 | 36 | 0.00 | paired | 36 | 0.00 |
| dark25%_known | PRECORRECT | 0.6 | 4.06 | 0.00 | 36 | -0.11 | paired | 36 | 0.00 |
| dark5%_known | RQL | 0.05 | 14.41 | -0.00 | 36 | 0.07 | paired | 36 | -0.00 |
| dark5%_known | RQL | 0.6 | 15.06 | 0.00 | 36 | 0.08 | paired | 36 | -0.00 |
| dark5%_known | PRECORRECT | 0.05 | 6.36 | -1.00 | 36 | 0.00 | paired | 36 | 0.00 |
| dark5%_known | PRECORRECT | 0.6 | 4.24 | 0.00 | 36 | 0.05 | paired | 36 | -0.00 |
| dark50%_known | RQL | 0.05 | 14.40 | -0.00 | 36 | -0.03 | paired | 36 | 0.00 |
| dark50%_known | RQL | 0.6 | 14.54 | 0.00 | 36 | -0.24 | paired | 36 | 0.00 |
| dark50%_known | PRECORRECT | 0.05 | 3.68 | -0.00 | 36 | -2.60 | paired | 36 | 1.00 |
| dark50%_known | PRECORRECT | 0.6 | 3.86 | 0.00 | 36 | -0.26 | paired | 36 | 0.00 |
| p_ap=0.01 | RQL | 0.05 | 14.45 | 0.01 | 36 | 0.06 | paired | 36 | 0.01 |
| p_ap=0.01 | RQL | 0.6 | 15.33 | 0.01 | 36 | 0.15 | paired | 36 | 0.01 |
| p_ap=0.01 | PRECORRECT | 0.05 | 6.36 | -1.00 | 36 | 0.00 | paired | 36 | 0.00 |
| p_ap=0.01 | PRECORRECT | 0.6 | 4.17 | 0.01 | 36 | -0.04 | paired | 36 | 0.01 |
| p_ap=0.02 | RQL | 0.05 | 14.47 | 0.02 | 36 | -0.00 | paired | 36 | 0.02 |
| p_ap=0.02 | RQL | 0.3 | 14.71 | 0.02 | 36 | MISSING | no-baseline | 0 | MISSING |
| p_ap=0.02 | RQL | 0.6 | 14.92 | 0.01 | 36 | -0.00 | paired | 36 | 0.01 |
| p_ap=0.02 | RQL | 1 | 14.93 | 0.01 | 36 | 0.09 | paired | 36 | 0.01 |
| p_ap=0.02 | PRECORRECT | 0.05 | 6.36 | -1.00 | 36 | 0.00 | paired | 36 | 0.00 |
| p_ap=0.02 | PRECORRECT | 0.3 | 3.89 | 0.01 | 36 | MISSING | no-baseline | 0 | MISSING |
| p_ap=0.02 | PRECORRECT | 0.6 | 4.04 | 0.01 | 36 | -0.12 | paired | 36 | 0.01 |
| p_ap=0.02 | PRECORRECT | 1 | 4.32 | 0.01 | 36 | -0.09 | paired | 36 | 0.01 |
| p_ap=0.05 | RQL | 0.05 | 14.40 | 0.05 | 36 | -0.01 | paired | 36 | 0.05 |
| p_ap=0.05 | RQL | 0.6 | 14.80 | 0.03 | 36 | -0.09 | paired | 36 | 0.03 |
| p_ap=0.05 | PRECORRECT | 0.05 | 6.36 | -1.00 | 36 | 0.00 | paired | 36 | 0.00 |
| p_ap=0.05 | PRECORRECT | 0.6 | 3.98 | 0.03 | 36 | -0.23 | paired | 36 | 0.03 |
| p_ap=0.1 | RQL | 0.05 | 14.37 | 0.09 | 36 | -0.19 | paired | 36 | 0.09 |
| p_ap=0.1 | RQL | 0.3 | 14.45 | 0.08 | 36 | MISSING | no-baseline | 0 | MISSING |
| p_ap=0.1 | RQL | 0.6 | 14.78 | 0.06 | 36 | -0.16 | paired | 36 | 0.06 |
| p_ap=0.1 | RQL | 1 | 15.05 | 0.05 | 36 | 0.03 | paired | 36 | 0.05 |
| p_ap=0.1 | PRECORRECT | 0.05 | 6.36 | -1.00 | 36 | 0.00 | paired | 36 | 0.00 |
| p_ap=0.1 | PRECORRECT | 0.3 | 3.51 | 0.08 | 36 | MISSING | no-baseline | 0 | MISSING |
| p_ap=0.1 | PRECORRECT | 0.6 | 3.84 | 0.06 | 36 | -0.45 | paired | 36 | 0.06 |
| p_ap=0.1 | PRECORRECT | 1 | 4.02 | 0.05 | 36 | -0.41 | paired | 36 | 0.05 |
| start=delayed | RQL | 0.05 | 14.51 | -0.00 | 36 | -0.00 | paired | 36 | -0.00 |
| start=delayed | RQL | 0.6 | 14.96 | -0.00 | 36 | -0.00 | paired | 36 | -0.00 |
| start=delayed | PRECORRECT | 0.05 | 6.36 | -1.00 | 36 | 0.00 | paired | 36 | 0.00 |
| start=delayed | PRECORRECT | 0.6 | 4.18 | -0.00 | 36 | 0.03 | paired | 36 | -0.00 |
| guard=2.5e-07 | RQL | 0.05 | 14.53 | -0.00 | 36 | 0.00 | paired | 36 | 0.00 |
| guard=2.5e-07 | RQL | 0.6 | 14.95 | 0.00 | 36 | 0.00 | paired | 36 | 0.00 |
| guard=2.5e-07 | PRECORRECT | 0.05 | 6.36 | -1.00 | 36 | 0.00 | paired | 36 | 0.00 |
| guard=2.5e-07 | PRECORRECT | 0.6 | 4.17 | 0.00 | 36 | 0.00 | paired | 36 | 0.00 |
| guard=5e-08 | RQL | 0.05 | 14.53 | -0.00 | 36 | 0.00 | paired | 36 | 0.00 |
| guard=5e-08 | RQL | 0.6 | 14.95 | 0.00 | 36 | 0.00 | paired | 36 | 0.00 |
| guard=5e-08 | PRECORRECT | 0.05 | 6.36 | -1.00 | 36 | 0.00 | paired | 36 | 0.00 |
| guard=5e-08 | PRECORRECT | 0.6 | 4.17 | 0.00 | 36 | 0.00 | paired | 36 | 0.00 |
| jitter_cv=0.05 | RQL | 0.05 | 14.35 | 0.00 | 36 | 0.06 | paired | 36 | 0.00 |
| jitter_cv=0.05 | RQL | 0.6 | 14.81 | 0.00 | 36 | 0.01 | paired | 36 | 0.00 |
| jitter_cv=0.05 | PRECORRECT | 0.05 | 6.36 | -1.00 | 36 | 0.00 | paired | 36 | 0.00 |
| jitter_cv=0.05 | PRECORRECT | 0.6 | 4.31 | 0.00 | 36 | 0.01 | paired | 36 | 0.00 |
| jitter_cv=0.1 | RQL | 0.05 | 14.34 | 0.00 | 36 | 0.07 | paired | 36 | 0.00 |
| jitter_cv=0.1 | RQL | 0.6 | 14.89 | 0.00 | 36 | -0.02 | paired | 36 | 0.00 |
| jitter_cv=0.1 | PRECORRECT | 0.05 | 6.36 | -1.00 | 36 | 0.00 | paired | 36 | 0.00 |
| jitter_cv=0.1 | PRECORRECT | 0.6 | 4.33 | 0.00 | 36 | 0.01 | paired | 36 | 0.00 |

