# CLAIM–SOURCE MATRIX — PRL Letter (ROUND63 / R43 FINAL LOCK)

**Purpose (R43 §7.3 step 1):** every sentence-level numerical claim that will appear in the
Letter — drawn from R43's frozen wording (§4.1–4.6) and the three figure specs (§6.2–6.4) —
mapped to exactly one committed frozen artifact: value → file → field/table → commit hash.

**Binding spec:** `docs/ROUND63_GPT_ROUND43_RULING_RAW.md`.
**Green-light protocol (R43 §7):** figures and claims are built ONLY from committed frozen
artifacts. No data generation. Any claim without a committed source is flagged as a GAP below
and is NOT filled by regeneration.

**Repo HEAD at build time:** `18891c65` (branch `main`).

## Frozen source commits

| Artifact | Frozen commit | Working-tree state | Extracted copy used by figure scripts |
|---|---|---|---|
| `results/round63_next/JET_TEST/JET_TEST.json` | **`1bf29f1`** | clean (restored; GAP-1 RESOLVED) | `paper_prl/figures/_frozen_sources/JET_TEST.committed.json` |
| `results/round63_next/FOG_DMD_PROBE64/CORRECTION_NOTE.md` + `P1_results_CORRECTED.json` | `099ce7f` | clean (== committed) | (End Matter C disclosure 2 numbers) |
| `results/round63_next/SCRAMBLE_EXT/SCRAMBLE_RESULTS.json` | `ed7a1e0` | clean (== committed) | `paper_prl/figures/_frozen_sources/SCRAMBLE_RESULTS.committed.json` |
| `results/round63_next/SEALED_DET/CONFIRMATORY_RESULTS.json` | `b37c841` | clean (== committed) | `paper_prl/figures/_frozen_sources/CONFIRMATORY_RESULTS.committed.json` |
| `results/round63_next/SCRAMBLE_EXT/DERIVATION.md` | `ed7a1e0` | clean | (prose/equations) |
| `docs/ROUND63_GPT_ROUND42_RULING_RAW.md` | HEAD (`18891c65`) | tracked | (theorem equations) |

The `CONFIRMATORY_RESULTS.json` internal **analysis-freeze** commit is `5910277` (field
`freeze_commit`); the artifact itself was committed at `b37c841`. Both are recorded.

---

## A. R43 §4.1 — Detection capability (thin-screen sealed test)

| # | Claim (frozen value) | File | Field / table | Commit |
|---|---|---|---|---|
| A1 | `2%` beyond-band detectable in **77/81** cells | CONFIRMATORY | `D2.eps2_cells` = 77 (`D2.eps2_pass_cells`=77); total `D2.eps5_cells` = 81 | b37c841 |
| A2 | best cell required **453** banks | CONFIRMATORY | `D2.best_Tdet_2pct` = 453.0 (= `D2.eval.measured.best_cell_Tdet`); map cell {sf1.0, k^-2, kwf1, claim1.25} `T_det`=453 | b37c841 |
| A3 | Monte Carlo power lower bound **0.990** | CONFIRMATORY | `D2.mc_power_lcb` = 0.9899 (= `D2.eval.measured.eps2_mc_lb`); pooled n `D2.mc_pooled_n`=13500 | b37c841 |
| A4 | **459-versus-1013**-bank latency advantage (fixed vs fresh) | CONFIRMATORY | `D4.rows[0]` (cell "best"): `fixed_latency`=459.0, `fresh_latency`=1013.0 | b37c841 |
| A5 | online scoring consumed **0.16%** of each bank interval | CONFIRMATORY | `D6.eval.measured.online_frac_of_bank` = 0.0016; `online_ms`=0.0203 | b37c841 |
| A6 | four-class balanced accuracy **0.916** at frozen 1% threshold | CONFIRMATORY | `D3.balanced_accuracy` = 0.916; threshold `D3.threshold_1pct_FA`=0.0510 | b37c841 |
| A7 | in-band false alarms **0.020** | CONFIRMATORY | `D3.fa.inband` = 0.02 | b37c841 |
| A8 | medium-amplitude / medium-timescale false alarms **0.096** and **0.084** | CONFIRMATORY | `D3.fa.amplitude` = 0.096, `D3.fa.timescale` = 0.084 | b37c841 |

## B. R43 §4.2 — Robustness domain

| # | Claim (frozen value) | File | Field / table | Commit |
|---|---|---|---|---|
| B1 | primary domain = aperture-preserving basis rotations up to **10%** (**FA=0.032**, latency inflation **19%**) | CONFIRMATORY | `D5.rows` axis `basis_rotation` level 0.1: `nontarget_fa`=0.032, `Tdet_inflation`=0.192 | b37c841 |
| B2 | at **20%** rotation and spectral-slope **-1**: non-target FA rose to **0.052** and **0.076** | CONFIRMATORY | `D5.rows` basis_rotation 0.2 `nontarget_fa`=0.052; spectral_slope -1 `nontarget_fa`=0.076 | b37c841 |
| B3 | robustness frozen to ≤10% rotation regime | CONFIRMATORY | `D5.eval.first_primary_fail`="basis_rotation=0.2"; `kill_tree.verdict`="narrow the physical domain or kill robustness wording" | b37c841 |

## C. R43 §4.3 — D3 disclosure (specificity killed, power intact)

| # | Claim (frozen value) | File | Field / table | Commit |
|---|---|---|---|---|
| C1 | target TPR **0.988** | CONFIRMATORY | `D3.tpr` = 0.988 (`per_class_recall.beyond`=0.98) | b37c841 |
| C2 | beyond score centered on non-targets, **\|d'\|≤0.02** | CONFIRMATORY | `D3.beyond_nontarget_dp` = [-0.0192, -0.0170, -0.0045]; max\|·\|=0.0192 | b37c841 |
| C3 | balanced accuracy **0.916** | CONFIRMATORY | `D3.balanced_accuracy` = 0.916 (dup A6) | b37c841 |
| C4 | medium tails exceeded frozen **0.05** bar (**0.096** and **0.084**) | CONFIRMATORY | `D3.fa.amplitude`/`timescale`; bar in checks `fa_amplitude_le_5`=false, `fa_timescale_le_5`=false | b37c841 |
| C5 | weakest class-specific separation **1.79** (not 5) | CONFIRMATORY | `D3.intended_dp` = [1.7928, 43.06, 4.78]; min=1.79; check `intended_scores_dp_ge_5`=false | b37c841 |

## D. R43 §4.4 — D5 disclosure

| # | Claim (frozen value) | File | Field / table | Commit |
|---|---|---|---|---|
| D1 | D5 preserved AUC **1.000** on every tested axis | CONFIRMATORY | `D5.rows[*].auc` all = 1.0; `D5.eval.base_auc`=1.0, `base_dprime`=6.91 | b37c841 |
| D2 | bounded latency inflation on every axis | CONFIRMATORY | `D5.rows[*].Tdet_inflation` (range −0.035…+0.249) | b37c841 |
| D3 | **0.052** at 20% rotation, **0.076** at spectral slope −1 | CONFIRMATORY | (dup B2) | b37c841 |

## E. R43 §4.5 — Generator-chunking integrity (qualitative + provenance)

| # | Claim | File | Field / table | Commit |
|---|---|---|---|---|
| E1 | post-freeze refactor streamed MC records in GPU-memory chunks; no change to hypotheses, thresholds, bank counts, scenes, estimators, endpoints; altered RNG → different still-blinded realization set | CONFIRMATORY | `D7.attestation`; `freeze_commit`="5910277"; `banks.confirmatory`=12 | b37c841 |
| E2 | shot-model integrity disclosure (End Matter C disclosure 2): signed-code shot error inflated covariance Fisher "≈1.9×"; witness `f_rec` "0.050 → 0.033" | FOG_DMD_PROBE64 | `CORRECTION_NOTE.md` ("Fisher over-stated ≈ 1.9×"; table row "witness f_rec 0.050 → **0.033**"); `P1_results_CORRECTED.json:f_rec_snr3`=0.03333 | 099ce7f |
| E3 | bench-transfer validity boundary (End Matter C, R44-frozen wording): mean leak "≈0.7%" (near-conjugate) → "≈9%" (20 mm gap); thin-screen spacing "≤0.3 mm" (developed) / "1–2 mm" (weak) | WAVE_TWIN | `T1_WALL_LEAK.json` (z1 sweep 6.8e-3 → 8.9e-2); `T4_MULT_TO_CONV.json` (α=0.5 at z2=0.327 mm developed; weak holds to ~1.4–2 mm); frozen wording = R44 §3 (issue #35) | 60d5fc7 |

## F. R43 §4.6 / §6.3 — Quotient-jet theorem statement + measured validation

**Theorem equations (analytical — no data):**

| # | Statement | File | Location | Commit |
|---|---|---|---|---|
| F0a | `C_*(ε;δ)=c_m(δ) ε^{2m}+o(ε^{2m})`, `c_m>0` ⇒ `T_req ≍ ε^{-2m}` | ROUND42 ruling | §6.1 eqs (6.2)(6.4); §1.1 eqs (1.2)–(1.4) | HEAD 18891c65 |
| F0b | `Σ(x)=R+[xᵀGx]H`, `G>0`; generic dir `m=1`, first-order-orthogonal `m=2`; no nonzero local blind dir though scene Fisher rank = 1 | ROUND42 §6.2 (6.13)(6.21–6.23); DERIVATION §3–4 (3.3)(3.6)(4.1); R43 §2 (2.4) | HEAD / ed7a1e0 |

**Measured validation (from FROZEN JET_TEST commit `1bf29f1` — see GAP-1):**

| # | Claim (frozen value) | Field in `JET_TEST.committed.json` | Commit |
|---|---|---|---|
| F1 | Exact KL orders **2.038** and **4.000** | `bank_A.per_dir.delta_g.kl_slope`=2.0381; `delta_o.kl_slope`=3.9999 (delta_n 4.0000) | 1bf29f1 |
| F2 | coefficient agreement (generic 1.07%, ortho 0.00%) | `per_dir.delta_g.coef_ratio`=1.0107; `delta_o.coef_ratio`=1.0000 | 1bf29f1 |
| F3 | Monte Carlo response exponents **0.95** and **2.05** | `bank_A_montecarlo.delta_g.dprime_vs_eps_slope`=0.9547; `delta_o`=2.0455 | 1bf29f1 |
| F4 | CUSUM delay slopes **-2.16** and **-3.92** | `cusum.delta_g.delay_slope`=-2.1583; `delta_o.delay_slope`=-3.9162 | 1bf29f1 |
| F5 | predicted crossover (ratio 1.000) at `ε_cross` | `bank_A.crossover.eps_cross_pred`=0.009995, `eps_cross_emp`=0.009997, `ratio`=1.0003 | 1bf29f1 |
| F6 | nuisance collapse `J_θ=0` (amplitude gauge) | `bank_B.single_zero_lag.I_theta_eff`=0.0 (`proportional_lags.I_theta_eff`=0.0) | 1bf29f1 |
| F7 | anchor restoration `J_θ>0` (finite prior recovers 20%) | `bank_B.amplitude_anchor.frac_recovered`=0.20 (2nd-wavelength 0.018) | 1bf29f1 |
| F8 | iso-`Q` cancellation → AUC ≈ chance at `ε_*` | `bank_C.cancellation.eps_star`=0.3626, `auc_at_star`=0.4807; local visible either side (auc 0.026 / 1.0; JSON `auc_left`=0.02638) | 1bf29f1 |
| F9 | monotone-cone corollary (ΔQ>0 all tested) | `bank_C.cone.all_positive`=true; `G_entrywise_nonneg`=true | 1bf29f1 |
| F10 | no broad in-aperture blind set | `kill_scan.frac_blind`=0.0; `d_min`=0.2544 over 400 dirs | 1bf29f1 |
| F11 | verdict MOONSHOT_SURVIVES, 17/17 checks, 0/4 kills | `verdict`, `checks` (all true), `kills` (all false) | 1bf29f1 |

## G. R43 §6.2 — Figure 1 (hero) content

| # | Element | File | Field / location | Commit |
|---|---|---|---|---|
| G1 | mean spatial support `B_p → {0}` | DERIVATION Thm 1 (`B_mean={0}`); SCRAMBLE `valA_mean_null.delta_mean_rel`=6.560e-19 | ed7a1e0 |
| G2 | reconstruction score rank many → **1** (99.99% energy in one direction) | DERIVATION §3.2 (3.6); SCRAMBLE `valB_rank1.frac_cov_energy_in_OhO_direction`=0.99992 | ed7a1e0 |
| G3 | quotient jet order finite `m=1/2` → finite `m=1/2` (survives) | JET_TEST `kl_slope` 2.038/4.000; `kill_scan.frac_blind`=0 | 1bf29f1 |
| G4 | `Cov(b)=Q(x)(O∘O)+R`, `Q=xᵀGx` | DERIVATION eq (3.3); SCRAMBLE `valB_cov_is_Q` (`Q_x_theory`=187.4, `var_over_O2Q_ratio_median`=1.019) | ed7a1e0 |
| G5 | two directional classes `u=2xᵀGδ`, `v=2δᵀGδ` | R43 §2.4; SCRAMBLE `valB_dQ.terms` {cross, quad}; per_delta ortho/coherent | ed7a1e0 |
| G6 | amplitude-anchor condition (icon) | R43 §2 boundary #2; JET_TEST `bank_B` (gauge collapse `J=0` → anchor `J>0`) | 1bf29f1 |

## H. R43 §6.4 — Figure 3 (two optical realizations + honest boundary)

| # | Element | File | Field / location | Commit |
|---|---|---|---|---|
| H1 | thin-screen 2% power: **77/81**, best **453** banks, MC LCB **0.990** | CONFIRMATORY | `D2` (dup A1–A3) | b37c841 |
| H2 | fixed-vs-fresh latency **459 vs 1013** | CONFIRMATORY | `D4.rows[0]` (dup A4) | b37c841 |
| H3 | calibration boundary: bal_acc **0.916**, medium FA **0.096/0.084**, ≤**10%** rotation | CONFIRMATORY | `D3`/`D5` (dup A6, A8, B1) | b37c841 |
| H4 | complete-scrambling mean AUC at chance | SCRAMBLE | `valC_mean_detector.coherent.auc`=0.5016, `ortho.auc`=0.5045 | ed7a1e0 |
| H5 | covariance rank-one functional (**99.99%**) | SCRAMBLE | `valB_rank1.frac_cov_energy_in_OhO_direction`=0.99992 | ed7a1e0 |
| H6 | generic **5%** duel: **d'=4.1**, AUC **0.997** | SCRAMBLE | `detection.coherent.eps_0.05."8192"`: `d_prime`=4.098, `auc`=0.997 | ed7a1e0 |
| H7 | mean-null magnitude **6.6e-19** | SCRAMBLE | `valA_mean_null.delta_mean_rel`=6.560e-19 | ed7a1e0 |

## I. Shared model configuration (Methods / captions)

| # | Value | File | Field | Commit |
|---|---|---|---|---|
| I1 | grid 32×32 (`N`=1024), code band `PB`=3, grain band `K_GRAIN`=13, `M`=36 codes, `PHOT`=1e4 | JET_TEST / SCRAMBLE `config` (identical) | 1bf29f1 / ed7a1e0 |
| I2 | scrambling aperture: `Ghat_nonzero_bins`=1024/1024 (full grid), `C0`=0.5166, grain_band ≥ pattern_band | SCRAMBLE `aperture` | ed7a1e0 |
| I3 | sealed run = 12 confirmatory banks, 2.76/6.0 GPU-h, within ceiling | CONFIRMATORY `banks.confirmatory`=12, `total_gpu_hours`=2.764, `within_ceiling`=true | b37c841 |

---

# GAPS AND DISCREPANCIES (flagged, NOT regenerated)

## GAP-1 — RESOLVED (was CRITICAL): JET_TEST working tree diverged from the frozen commit

**RESOLUTION (coordinator, option (a) executed):** `git checkout` restored
`results/round63_next/JET_TEST/` to the frozen `1bf29f1` state (working tree verified clean);
the post-freeze re-run was quarantined to `results/round63_next/JET_TEST_POSTFREEZE_RERUN/`
(@ `4e235a5`, provenance-only, used nowhere in the Letter). Original narrative retained below
for the audit trail.

`results/round63_next/JET_TEST/` is **modified in the working tree** (git status `M`) relative
to the frozen commit **`1bf29f1`** referenced by R43 §Reference and §Decisive-exams. Modified
files: `JET_TEST.json`, `JET_TEST_REPORT.md`, `jet_slopes.{png,pdf}`, `bank_b_nuisance.{png,pdf}`,
`run_full.log`. The working-tree copy is an **uncommitted re-run** of `jet_test.py`.

The **committed frozen** values (`1bf29f1`) match R43 §4.6 exactly; the uncommitted re-run does not:

| Quantity | Frozen `1bf29f1` (== R43 §4.6, USED) | Uncommitted working-tree re-run (NOT used) |
|---|---|---|
| KL order, generic (`delta_g`) | **2.038** | 2.015 |
| KL order, orthogonal (`delta_o`) | 4.000 | 4.000 |
| MC exponent, generic | **0.95** | 1.17 |
| MC exponent, orthogonal | **2.05** | 1.99 |
| CUSUM slope, generic | **-2.16** | -2.19 |
| CUSUM slope, orthogonal | **-3.92** | -3.98 |
| generic coef ratio | 1.011 (1.07%) | 1.003 (0.27%) |
| iso-Q cancellation AUC at ε* | 0.481 | 0.516 |

**Resolution applied (per R43 §7 / §7.2 "disclose, do not regenerate"):** Figure 2 and all §4.6
claims are built from the committed frozen JSON at `1bf29f1`, extracted verbatim to
`paper_prl/figures/_frozen_sources/JET_TEST.committed.json`. The uncommitted re-run is **not
used** anywhere in the Letter.

**Both value sets support the identical qualitative conclusion** — integer contact orders 2 and 4,
MC exponents ≈1 and ≈2, CUSUM slopes ≈−2 and ≈−4, 17/17 checks pass, empty blind set. Only
second-decimal values differ; no new science is present in either.

**Action requested of coordinator (do NOT let me regenerate):** choose one —
(a) **RECOMMENDED** — `git checkout -- results/round63_next/JET_TEST/` to restore the frozen
artifact and remove the stray re-run (a re-run of `jet_test.py` after freeze is a §7.2-forbidden
regeneration); or
(b) formally re-freeze the re-run as a new commit and update R43 §4.6 decimals to 2.015 / 1.17 /
1.99 / -2.19 / -3.98. Figure 2 would then be rebuilt from the new commit.

Until the coordinator rules, the Letter cites the R43-consistent frozen values.

## GAP-2 — Non-blocking: no committed source computes a single scalar "T_req crossover" figure

R43 §6.3 asks for "a compact crossover marker at the predicted `ε_cross`". The committed JET_TEST
provides `bank_A.crossover.eps_cross_pred`=0.009995 and the empirical local-slope crossover
`eps_cross_emp`=0.009997 (ratio 1.000). This is sufficient for the marker; recorded here only to
note the marker is the **mixed direction** (`delta_mix`, `c`=0.687, `d`=137.5), not `delta_g`.

## GAP-3 — Non-blocking: End-Matter theorem proofs are Phase 2

The full Rank–Jet proof and quotient-jet derivation (R43 End Matter A/B, Supplement S1) live in
`docs/ROUND63_GPT_ROUND42_RULING_RAW.md` §6 and `SCRAMBLE_EXT/DERIVATION.md` §2–4. These are
prose/equation sources, not numerical claims, and are deferred to Phase 2 per the task scope. No
gap in evidence; flagged only so the coordinator knows the proof text is not yet transcribed into
the manuscript.

**No numerical claim in the Letter (R43 §4.1–4.6, figure specs, or End Matter C disclosures)
lacks a committed source.** Coverage extended to the shot-model disclosure numbers via row E2
(hostile-read finding M4). GAP-1 is resolved (working tree restored to the frozen commit;
re-run quarantined).
