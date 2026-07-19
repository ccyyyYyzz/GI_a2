# M1 FREEZE CHECKLIST LEDGER (R13 Section 10)

Outcome-blind selftest, 2026-07-19 05:48. 14/16 hard items PASS.

| # | item | verdict | detail |
|---|---|---|---|
| 1 | prescan_row_average_all_ones | PASS | max|mean-1|=0.00e+00 |
| 2 | prescan_pixel_dose_exactly_equal | PASS | dose=52 per pixel, spread=0.00e+00 |
| 3 | V0_actual_row_loads_no_oracle | PASS | row loads spread sd=9.226e-02; V0 differs from oracle version |
| 4 | dictionary_manifest_complete | PASS | 13 families, global sha 710aa49f43d0.. |
| 5 | fixed_star_scores_committed | PASS | FIXED*=MATCH1  scores={'SCAT32': 3591.6, 'LBLOB16': 3736.3, 'MATCH1': 3736.8} |
| 6 | all_arms_52_plus_972 | PASS | SCAT16:972, SCAT32:972, LBLOB16:972, MATCH1:972, RIDGE-FIXED:972 (+52 pre-scan) |
| 7 | nine_dwell_kernel_trust_bias_cert | PASS | 18 (nu, policy) palettes certified |
| 8 | cert_convention_toy_verified | PASS | G_rel=0.0e+00 at the hand-solved optimum |
| 9 | fixed_dose_972_exact_rows | PASS | dose dev=0.0376 (<=0.05), 972 exact rows |
| 10 | safe_and_fast_oed_through_final_solver | PASS | both policies solved |
| 11 | relaxed_kw_upper_bound_1e-3 | FAIL | dose-feasible: fast=2.73e+00 safe=2.86e+00 (/r); relaxed-solver: fast=2.94e-02 safe=2.46e-02 -- if the dose-feasible values plateau above 1e-3 while the relaxed solver converges, the binding per-pixel dose band is the cause (R13 Sec 6.1 tolerance amendment evidence) |
| 12 | mixture_materialized_exact_rows | PASS | machinery demonstrated: exact 972-row OED+FIXED_DOSE stack; in-pipeline trigger not needed (dose passed directly) |
| 13 | final_exact_rows_full_guard_suite | FAIL | fast: budget=True dose=False peak=True ceil=True effD=True arisk=False spectral=False; A5 bound=7.00 |
| 14 | analyzer_emits_three_verdicts_on_dev_mock | PASS | all three verdict strings emitted on BOTH cohorts; positive all PASS, negative all FAIL |
| 15 | no_hard_check_uses_endpoint_psnr_or_conf_scenes | PASS | hard checks touch: synthetic scene, m1_dev mocks, kernel tables, design matrices; zero PSNR reads, zero m1_ confirmatory loads |
| 16 | manifests_and_expected_cells_regenerated | PASS | 15120 expected cells, 40 shards |

R13 verdict rule: GO only when every box passes on an immutable commit. Current run: HOLD REMAINS (failures above).

wall time: 720.6 s
