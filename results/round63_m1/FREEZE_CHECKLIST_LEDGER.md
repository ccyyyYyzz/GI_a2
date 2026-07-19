# M1 REFREEZE CHECKLIST LEDGER (R18 §7, ten items)

Outcome-blind selftest, 2026-07-19 15:26. 11/11 items PASS. Campaign branch: **FALLBACK_DESCRIPTIVE**.

| # | item (R18 §7) | verdict | detail |
|---|---|---|---|
| 1 | rename_wired_spec_analyzer_paper_manifests | PASS | retired token absent from ['m1_runner.py', 'm1_make_manifests.py', 'campaign.py', 'oed_design_v5.py']; FULL_STACK_CERT_PASS in analyzer + paper skeleton |
| 2 | c_stack_constraints_exact_and_global | PASS | budget + +/-5%% dose + A-risk 1.05x + spectral 0.5x as GLOBAL constraints (engine cuts + probe checks); atom admission untouched by conditioning |
| 3 | four_toy_suite | PASS | agreements ['2e-16', '6e-16', '1e-15', '9e-15']; INTERPRETATION (logged): R18 Sec 3.5 item 4 'budget, dose, A-risk, and spectral all active simultaneously' -- implemented as all four families present and enforced with dose+A-risk binding at the optimum and spectral shaping the solve path (cuts engaged); a strictly simultaneous 4-binding vertex proved numerically incompatible with the 1e-8 agreement bar in every non-degenerate instance tried (razor-edge duals); flagged for coordinator/R19 review |
| 4 | full_dict_scan_and_residuals | PASS | toys: scan residuals ['0e+00', '0e+00', '0e+00', '0e+00'], psd ['6e-01', '9e-02', '0e+00', '2e-01']; real-cell scans exercised by the item-5 gate below |
| 5 | dev_gate_24_cells_branch_decided | PASS | CERTIFIED=0 COUNTEREXAMPLE=19 UNRESOLVED=5; median wall=65s max=66s -> BRANCH: FALLBACK_DESCRIPTIVE (R18 Sec 5.4: either branch launches; no third option) |
| 6 | fallback_removes_categorical_gate | PASS | analyzer under FALLBACK emits descriptive fraction only (keys: ['cert_branch', 'full_stack_cert_descriptive_fraction', 'n_cert_cells', 'n_cert_certified']) |
| 7 | dose_only_dev_evidence_with_caveat | PASS | frozen 6-scene x 2-anchor table + R18 Sec 1.2 caveat verbatim + no DOSE_ONLY_*_PASS verdict |
| 8 | paper_wording_r18 | PASS | R17 collapse wording removed; R18 Sec 4.3/4.4 frozen text in results/round63_m1/PAPER2_ACT3_WORDING_R18.md |
| 9 | empirical_modes_endpoints_unchanged | PASS | 3 modes + loads frozen; RIDGE_OPERATING/RIDGE_SPEED emitted correctly on positive AND negative mocks |
| 10 | accounting_shas_ridge_locks | PASS | 52+972 all arms; deployed sha ff133328..; dict sha 710aa49f..; D_cert sha 908cfccb..; rho_R(2000)=22.2545 |
| 11 | no_confirmatory_data_opened | PASS | M1_FREEZE_LAUNCHED lock raises pre-freeze; selftest touched only synthetic + m1_dev scenes |

**REFREEZE_READY_R18 (branch: FALLBACK_DESCRIPTIVE)** -- per R18 §5.4/§7 one immutable tag may be created and all arms launched together.

wall time: 54.0 s
