# M1 REFREEZE CHECKLIST LEDGER (R17 §7 / amendment §E)

Outcome-blind selftest, 2026-07-19 13:05. 9/9 items PASS.

Architecture: R17 (issue #9). Deployed design = balanced exact 972-row SCAT32 multiset (sha ff133328f00b..). Modes: SCAT32-SAFE 0.05 / SCAT32-060 0.60 / RIDGE-SCAT32 rho_R(nu) policy. Context: SCAT16, LBLOB16 (descriptive). D_cert sha 908cfccbd249..

Interpretations I1-I5 logged in m1_runner.py docstring.

| # | item (R17 Sec 7) | verdict | detail |
|---|---|---|---|
| 1 | oed_dt_and_m0_pass_semantics_removed | PASS | prod tokens absent; ARMS=['SCAT32-SAFE', 'SCAT32-060', 'RIDGE-SCAT32', 'SCAT16', 'LBLOB16']; retired-arm loader raises: arm 'OED-DT' RETIRED_BY_R17 (issue #9); ; m=0 -> ADAPTIVE_COLLAPSE_UNDER_GUARDS (is_pass=False) |
| 2 | deployed_scat32_and_modes_frozen_hashed | PASS | 972x1024 sha=ff133328f00b.. dose_dev=0.0400; modes SAFE=0.05 060=0.60 RIDGE=rho_R(nu) policy; prescan sha=0b93e4cf9915.. |
| 3 | three_verdicts_emitted_pos_and_neg | PASS | RIDGE_OPERATING/RIDGE_SPEED/DOSE_SAFE_CERT on both cohorts; gain toy 1.8e-09; negative toy G=2000000001.0 |
| 4 | d_cert_complete_and_sha_frozen | PASS | L=11 (load 6 + gain 5); admissible: 151344 load + 133120 gain atoms; sha=908cfccbd249de22.. (23s) |
| 5 | expanded_class_cert_toys | PASS | base toy 1.8e-09, gain-atom toy 1.8e-09, dose band ACTIVE in both |
| 6 | dev_expanded_class_feasibility_runtime | PASS | certified finite bounds (walls: [30.0, 24.3] s); gate bar 1e-2 UNCHANGED; G_full/r: [26.6761383383072, 79.10859500908491] |
| 7 | accounting_52_plus_972_and_ridge_disclosures | PASS | 5 arms x (52+972); 9-dwell calibration + full disclosure block (achieved=22.255 at nu=2000, ratio=37.09x) |
| 8 | manifests_and_expected_cells_r17 | PASS | 5880 cells (480 cert of 480 required), 43 shards; stale design frozen_inputs: none |
| 9 | no_confirmatory_data_opened | PASS | confirmatory-imageset guard raises pre-freeze; selftest touched only the synthetic scene and m1_dev (632900+) instances |

**REFREEZE_READY** -- every Sec-7 item passes; one immutable revised m1-freeze tag may be created and all arms launched together (R17 Sec 5).

wall time: 83.8 s
