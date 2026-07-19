# M1 verdicts — spec-conformant correction (R19)

Analysis version: `R19_SPEC_CORRECTION_V1`  
Date: 2026-07-19  
Immutable input: tag `m1-freeze` @ `6f00932` (`6f00932710dff5ea3afefc0faee54f231e9b915a`)  
Ruling: `docs/ROUND63_GPT_ROUND19_RULING_RAW.md` (R19, GitHub issue #11)

This is the human-readable companion to the corrected source of record `M1_VERDICTS_SPEC_CORRECTED_R19.json`. It is an additive, provenance-preserving implementation-to-specification correction. The immutable tag, raw cells, frozen endpoint definitions, scene cohort, reconstructions, and thresholds are unchanged. The endpoints below were recomputed from the same frozen raw cells using the preregistered elapsed optical-time axis `T_opt = M_total * nu * tau` and the full spec-conformant Q90 / equal-weight-PAVA / censoring-taxonomy / nested family-stratified bootstrap machinery.

## Corrected verdicts (source of record)

**RIDGE_OPERATING_PASS = True**

- median dQ = 1.866920 dB (bar >= 1.0)
- LB2.5 = 1.41348975 (bar > 0)
- n_pos = 19/24 (bar >= 18)

**RIDGE_SPEED_PASS = True** (elapsed optical time `T_opt = M_total * nu * tau`)

- median S_gate = 19.127043091646133 (bar >= 3)
- S LB2.5 = 18.328492357080282 (bar > 1)
- n_S>1 = 21/24 (bar >= 18)

**Post-hoc nu*rho incident-exposure sensitivity** (no PASS/FAIL verdict)

- median = 0.27864687917406294
- LB2.5 = 0.23974692258639926
- above 1 = 1/24

Bootstrap: B=10000, RNG operating tag 13, speed tag 14, `numpy.random.default_rng([20260717, 63, tag])`.

## Original vs corrected (provenance)

The **Original** column reproduces the frozen nonconformant-implementation outputs from `results/round63_m1/M1_VERDICTS.json` (preserved unchanged). These original values — in particular the speed `median 0.2759082546618772 / LB2.5 0.17225913049832775 / 1-of-24 FAIL` on the `log(nu*rho)` axis and the operating `LB2.5 0.11989000000000072` from the non-frozen bootstrap — are shown here **only** for provenance and must not be used as scientific results.

| Metric | Original (frozen, nonconformant implementation) | Corrected (R19 spec-conformant) |
|---|---|---|
| Operating median dQ (dB) | 1.8669199999999986 | 1.86692 |
| Operating LB2.5 | 0.11989000000000072 (non-frozen bootstrap, D2) | 1.41348975 |
| Operating n_pos | 19/24 | 19/24 |
| Operating verdict | PASS | PASS |
| Speed axis | `log(nu*rho)` incident exposure (D1) | elapsed `log T_opt = log(M_total*nu*tau)` |
| Speed median S | 0.2759082546618772 | 19.127043091646133 |
| Speed LB2.5 | 0.17225913049832775 | 18.328492357080282 |
| Speed n>1 | 1/24 | 21/24 |
| Speed verdict | FAIL | PASS |
| nu*rho role | reported as the speed endpoint | post-hoc incident-exposure sensitivity, no verdict (median 0.27864687917406294, LB2.5 0.23974692258639926, 1/24) |

Discrepancies D1 (axis), D2 (bootstrap), D3 (PAVA/target/censoring) are detailed in `ANALYSIS_CORRECTION_DISCLOSURE.md`. A diagnostic-only 'axis-only' speed LB `16.577852783743403` (breadth 22/24) replaces only the axis but preserves D2/D3; it is not a manuscript result.

## Per-image RIDGE_OPERATING terminal quality gain

dQ = mean over 5 frozen seeds of (RIDGE-SCAT32 − SCAT32-060) PSNR_rad at nu=2000.

| Image | dQ (dB) | dQ > 0 |
|---|---:|:--:|
| m1_chirp_0 | 2.361180 | yes |
| m1_chirp_1 | 1.270420 | yes |
| m1_chirp_2 | 1.005020 | yes |
| m1_chirp_3 | 1.538840 | yes |
| m1_contour_0 | -3.612180 | no |
| m1_contour_1 | -5.794460 | no |
| m1_contour_2 | -4.738300 | no |
| m1_contour_3 | 4.151300 | yes |
| m1_glyph_0 | 4.734300 | yes |
| m1_glyph_1 | 6.005160 | yes |
| m1_glyph_2 | 2.430660 | yes |
| m1_glyph_3 | 2.212580 | yes |
| m1_maze_0 | 6.329740 | yes |
| m1_maze_1 | 6.436420 | yes |
| m1_maze_2 | 5.583480 | yes |
| m1_maze_3 | 6.408240 | yes |
| m1_microtexture_0 | 1.529520 | yes |
| m1_microtexture_1 | 1.252680 | yes |
| m1_microtexture_2 | 1.525840 | yes |
| m1_microtexture_3 | 1.448760 | yes |
| m1_spokes_0 | -0.765240 | no |
| m1_spokes_1 | -5.084720 | no |
| m1_spokes_2 | 4.530260 | yes |
| m1_spokes_3 | 2.195000 | yes |

Median dQ = 1.866920 dB; positive in 19/24 images.

Per-family median dQ: chirp 1.4046, contour -4.1752, glyph 3.5825, maze 6.3690, microtexture 1.4873, spokes 0.7149.

## Per-image RIDGE_SPEED (elapsed optical time)

S_gate = Q90 time-to-quality ratio (SCAT32-SAFE crossing / RIDGE-SCAT32 crossing) on `T_opt = M_total * nu * tau`, with the frozen censoring taxonomy applied.

| Image | S_gate | status | Q90 target (dB) |
|---|---:|---|---:|
| m1_chirp_0 | 17.800203 | NORMAL | 6.914794 |
| m1_chirp_1 | 20.459721 | NORMAL | 7.389890 |
| m1_chirp_2 | 18.431468 | NORMAL | 7.261996 |
| m1_chirp_3 | 19.261646 | NORMAL | 7.404890 |
| m1_contour_0 | 2.940635 | NORMAL | 9.099688 |
| m1_contour_1 | 0.000000 | FAST_RIGHT_CENSORED | 9.584150 |
| m1_contour_2 | 0.000000 | FAST_RIGHT_CENSORED | 10.442378 |
| m1_contour_3 | 18.784826 | NORMAL | 15.186056 |
| m1_glyph_0 | 18.992440 | NORMAL | 17.970960 |
| m1_glyph_1 | 21.500155 | NORMAL | 18.454512 |
| m1_glyph_2 | 18.312756 | NORMAL | 9.343498 |
| m1_glyph_3 | 20.593897 | NORMAL | 9.230494 |
| m1_maze_0 | 21.440543 | NORMAL | 13.789082 |
| m1_maze_1 | 26.814231 | NORMAL | 18.569958 |
| m1_maze_2 | 26.844183 | NORMAL | 15.161556 |
| m1_maze_3 | 29.012398 | NORMAL | 15.054804 |
| m1_microtexture_0 | 18.142114 | NORMAL | 14.435742 |
| m1_microtexture_1 | 1.000000 | SAFE_RANGE_UNINFORMATIVE | N/A |
| m1_microtexture_2 | 18.014590 | NORMAL | 14.874426 |
| m1_microtexture_3 | 15.355503 | NORMAL | 15.629622 |
| m1_spokes_0 | 46.150824 | NORMAL | 7.989012 |
| m1_spokes_1 | 20.384410 | NORMAL | 9.191172 |
| m1_spokes_2 | 21.208912 | NORMAL | 11.207580 |
| m1_spokes_3 | 22.034111 | NORMAL | 8.462386 |

Median S_gate = 19.127043091646133; S_gate > 1 in 21/24 images. `m1_microtexture_1` is SAFE_RANGE_UNINFORMATIVE (fitted safe range 0.27336 dB < 0.50 dB) and is clamped to S_gate=1 (target undefined); this censoring step sets breadth to 21/24. `m1_contour_1` and `m1_contour_2` are FAST_RIGHT_CENSORED (S_gate=0). All 24 images retained.

## Per-image post-hoc nu*rho incident-exposure sensitivity

Not the preregistered endpoint; no PASS/FAIL verdict. Reported for transparency only.

| Image | S_gate (nu*rho) | status |
|---|---:|---|
| m1_chirp_0 | 0.209262 | NORMAL |
| m1_chirp_1 | 0.310040 | NORMAL |
| m1_chirp_2 | 0.278868 | NORMAL |
| m1_chirp_3 | 0.248876 | NORMAL |
| m1_contour_0 | 0.010351 | NORMAL |
| m1_contour_1 | 0.000000 | FAST_RIGHT_CENSORED |
| m1_contour_2 | 0.000000 | FAST_RIGHT_CENSORED |
| m1_contour_3 | 0.305285 | NORMAL |
| m1_glyph_0 | 0.229447 | NORMAL |
| m1_glyph_1 | 0.280244 | NORMAL |
| m1_glyph_2 | 0.247760 | NORMAL |
| m1_glyph_3 | 0.273390 | NORMAL |
| m1_maze_0 | 0.293458 | NORMAL |
| m1_maze_1 | 0.414036 | NORMAL |
| m1_maze_2 | 0.425202 | NORMAL |
| m1_maze_3 | 0.474813 | NORMAL |
| m1_microtexture_0 | 0.189512 | NORMAL |
| m1_microtexture_1 | 1.000000 | SAFE_RANGE_UNINFORMATIVE |
| m1_microtexture_2 | 0.204535 | NORMAL |
| m1_microtexture_3 | 0.152605 | NORMAL |
| m1_spokes_0 | 1.401020 | NORMAL |
| m1_spokes_1 | 0.278426 | NORMAL |
| m1_spokes_2 | 0.282317 | NORMAL |
| m1_spokes_3 | 0.311989 | NORMAL |

Median = 0.27864687917406294; above 1 in 1/24 images. Ridge tracking greatly reduced elapsed dwell but did not reduce incident exposure to Q90 under this coordinate.

## Full-stack certificate — descriptive (not a gate)

Branch = FALLBACK_DESCRIPTIVE (selected before the immutable tag). Reported descriptively; no categorical certificate verdict.

- 0 CERTIFIED / 299 COUNTEREXAMPLE / 181 NUMERICAL_UNRESOLVED (480 cells)
- descriptive certified fraction = 0.0
- by anchor: {"nu200_b0.05": {"NUMERICAL_UNRESOLVED": 43, "COUNTEREXAMPLE": 77}, "nu200_b0.6": {"COUNTEREXAMPLE": 67, "NUMERICAL_UNRESOLVED": 53}, "nu2000_b0.05": {"COUNTEREXAMPLE": 84, "NUMERICAL_UNRESOLVED": 36}, "nu2000_b0.6": {"COUNTEREXAMPLE": 71, "NUMERICAL_UNRESOLVED": 49}}

Per-status field availability, the hard-coded `arisk_ratio=1.0` disclosure, and the true coefficient-range distribution are in `results/round63_m1/CERT_DISCLOSURE_ADDENDUM_2026-07-19.md`.

## Provenance and reproduction

- Immutable input: tag `m1-freeze` @ `6f00932`; raw shard sha256 list in `M1_VERDICTS_SPEC_CORRECTED_R19.json` (`provenance.raw_shard_sha256`).
- Promoted independent implementation: `scripts/m1_independent_recompute.py`, `scripts/m1_frozen_and_axis_only.py`; outputs in `outputs/`.
- Verification: 18/18 audited numbers matched, zero deep-diff leaves (`results/round63_m1/AUDIT_2026-07-19/RERUN_REPORT.md`).
- Auditor: independent peer auditor, docs2/01, 2026-07-19 (implementation written before reading `code/round63/m1_score.py`).
- Environment: Windows 11, Python 3.11.5, numpy 1.24.4, pandas 1.5.3.
