# M1 FREEZE CHECKLIST LEDGER (R13 Sec 10 as amended by R15)

Outcome-blind selftest, 2026-07-19 10:50. 13/14 hard items PASS.

FIXED* selection rule (frozen, R10): Before confirmatory freeze, select one global fixed comparator FIXED* from {SCAT32, LBLOB16, MATCH1} using the six development images only. The selection score is the median fixed-dwell radiometric PSNR at (rho=0.60, nu=2000), averaged over development seeds. Tie within 0.05 dB is resolved in order SCAT32 -> LBLOB16 -> MATCH1, favoring the least scene-dependent rule. [ROUND63_GPT_ROUND10_RULING_RAW.md Sec 7 / line 801]

## Amended Box 11 wording (R15 Sec 5, verbatim)

```
[ ] relaxed_and_full_constrained_certificates

PASS iff, for BOTH the safe and fast OED policies:

(a) the dose-relaxed simplex+incident-budget reference problem, after the
    frozen support-restricted exact reweighting, satisfies
    RELAXED_REFERENCE_KW_UPPER_BOUND / r <= 1e-3 under a full-dictionary scan;

(b) the physically dose-feasible continuous design is primal-feasible for
    the incident budget and every +/-5% per-pixel dose inequality and satisfies
    FULL_CONSTRAINED_KW_UPPER_BOUND / r <= 1e-2, where the certificate includes
    one budget multiplier and all upper/lower per-pixel dose multipliers and
    the final maximum is scanned over the entire frozen dictionary;

(c) the full constrained certificate passes the independent active-dose toy
    verification and records all primal, dual, complementarity, and dictionary-
    scan residuals.
```

## Amended Box 13 wording (R15 Sec 7, verbatim)

```
[ ] final_exact_rows_full_guard_suite

PASS iff, for BOTH safe and fast policies, the final returned main-pattern
matrix contains exactly 972 physical rows and is either the direct rounded OED
design or the materialized COMPLIANT_VIA_MIXTURE design selected by the frozen
m=972,...,0 search. On that exact matrix, all of the following must pass:

- incident_sum_rho <= 972*budget_mean + 1e-9;
- predicted and realized detected-budget fields are emitted;
- per-pixel cumulative-dose max relative deviation <= 0.05;
- physical peak <= 1536;
- design-weighted ceiling probability <= 0.01;
- exact D-efficiency >= 0.95 relative to the certified full dose-constrained
  continuous design;
- A-risk <= 1.05 times the load-matched FIXED* reference;
- minimum generalized eigenvalue versus the load-matched FIXED* reference
  >= 0.5;
- exactly 52 common pre-scan rows plus 972 returned main rows are accounted.

If no integer mixture candidate, including m=0, passes all guards, Box 13 FAILS.
```

## Results

| # | item | verdict | detail |
|---|---|---|---|
| 1 | prescan_row_average_all_ones | PASS | max|mean-1|=0.00e+00 |
| 2 | prescan_pixel_dose_exactly_equal | PASS | dose=52 per pixel |
| 3 | V0_actual_row_loads_no_oracle | PASS | row-load sd=9.226e-02 |
| 4 | dictionary_manifest_complete | PASS | 13 families, sha 710aa49f43d0.. |
| 5 | fixed_star_r10_psnr_selection_committed | PASS | FIXED*=SCAT32 scores(dB)={'LBLOB16': 11.908, 'MATCH1': 9.181, 'SCAT32': 17.498} (tie window 0.05, order SCAT32->LBLOB16->MATCH1) |
| 6 | all_arms_52_plus_972 | PASS | SCAT16:972, SCAT32:972, LBLOB16:972, MATCH1:972, RIDGE-FIXED:972 (+52 pre-scan) |
| 7 | nine_dwell_kernel_trust_bias_cert | PASS | 18 palettes certified |
| 8 | toy_certificates_verified | PASS | unconstrained exact; active-dose agreement 1.8e-09 <= 1e-8 |
| 9 | fixed_dose_972_exact_rows | PASS | SCAT32-based (information-matched to FIXED*), dose dev=0.0400 |
| 10 | relaxed_and_full_constrained_certificates | FAIL | fast: relref=2.01e+00 full=inf; safe: relref=3.51e+00 full=inf |
| 11 | final_exact_rows_full_guard_suite | PASS | fast: PASS m=0; safe: PASS m=0 (order sha e040663b1446..) |
| 12 | analyzer_emits_three_verdicts_on_dev_mock | PASS | positive all PASS, negative all FAIL, all strings on both |
| 13 | no_hard_check_uses_confirmatory_scenes_or_endpoints | PASS | hard checks read: synthetic scene, m1_dev scenes (FIXED* rule per frozen R10 DEV-PSNR selection - sanctioned), kernel tables, design matrices; zero confirmatory loads |
| 14 | manifests_and_expected_cells_regenerated | PASS | 15120 expected cells, 40 shards, FIXED*=SCAT32 propagated |

HOLD REMAINS (failures above; trajectories in the selftest log).

wall time: 1644.2 s
