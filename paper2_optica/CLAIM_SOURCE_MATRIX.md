# CLAIM–SOURCE MATRIX — Paper 2 (Optica) / ROUND63 R45 §C FROZEN

**Paper:** *Engineering What an Optical Sensor Cannot See* (venue **Optica**; PR Applied fallback).
**Binding spec:** `docs/ROUND63_GPT_ROUND45_RULING_RAW.md` **Part C** (§C1–C7), architecture
**FROZEN** by the **TLSG_PASS** gate (`results/round63_next/TLSG/`, all 7 bars).

**Purpose (paper-1 §7.3 step-1 discipline):** every sentence-level numerical claim that will
appear in Paper 2 — drawn from the R45 §C2 theorem set and the §C3 section/figure architecture —
mapped to exactly one **committed frozen artifact**: value → file → field → commit hash. Figures
and claims are built ONLY from committed frozen artifacts extracted verbatim to
`paper2_optica/figures/_frozen_sources/`. **No data generation.** Any claim without a committed
source is flagged as a **GAP-n** below and is NOT filled by regeneration.

**Repo HEAD at build time:** `35d858e` (branch `main`; ≥ the required `35d858e` = TLSG_PASS commit).
**Nothing here touches `paper_prl/` (paper 1 frozen) or any sealed artifact.**

---

## Frozen source commits

| Artifact | Frozen commit | Working-tree state | Extracted copy used by figure scripts |
|---|---|---|---|
| `.../TLSG/TLSG_RESULTS.json` | **`35d858e`** (HEAD) | clean (== committed) | `_frozen_sources/TLSG_RESULTS.committed.json` |
| `.../R44_KILL_TESTS/KT1_shell_sweep.json` | `19338a4` | clean | `_frozen_sources/KT1_shell_sweep.committed.json` |
| `.../R44_KILL_TESTS/KT2_jet_transmutation.json` | `05f1029` | clean | `_frozen_sources/KT2_jet_transmutation.committed.json` |
| `.../R44_KILL_TESTS/KT3_leak_orthogonalized_sentinel.json` | `2430166` | clean | `_frozen_sources/KT3_leak_orthogonalized_sentinel.committed.json` |
| `.../R44_KILL_TESTS/KT6_channel_ratio_ranging.json` | `635ff05` | clean | `_frozen_sources/KT6_channel_ratio_ranging.committed.json` |
| `.../R44_KILL_TESTS/IT1_two_wall_jet_invariance.json` | `f13a323` | clean | `_frozen_sources/IT1_two_wall_jet_invariance.committed.json` |
| `.../R44_KILL_TESTS/IT2_floor_decomposition.json` | `acfcda0` | clean | (prose; leak-floor decomposition) |
| `.../R44_KILL_TESTS/IT3_leak_law_stability.json` | `6fb5149` | clean | `_frozen_sources/IT3_leak_law_stability.committed.json` |
| `.../R44_KILL_TESTS/IT4_svd_certified_blind_codes.json` | `33b8ef2` | clean | (§3 route-2 numbers) |
| `.../R44_KILL_TESTS/IT4b_spectrum_trade.json` | `5666d1e` | clean | `_frozen_sources/IT4b_spectrum_trade.committed.json` |
| `.../R44_KILL_TESTS/IT5_dpss_concentration.json` | `c8a0487` | clean | `_frozen_sources/IT5_dpss_concentration.committed.json` |
| `.../R44_KILL_TESTS/IT6_segmented_bucket_p2.json` | `dade877` | clean | (supplement negative) |
| `.../R44_KILL_TESTS/IT7_contact_coincidence_veto.json` | `e1546ac` | clean | (supplement negative) |
| `.../R44_KILL_TESTS/IT8_dither_render_leak.json` | `6fb5149` | clean | (apparatus PWM constraint) |
| `.../WAVE_TWIN/T1_WALL_LEAK.json` | `0c3b0d9` | clean | (§5 mean-leak floor) |
| `.../WAVE_TWIN/T5b_MSCALING.json` | `0c3b0d9` | clean | (§5 M^1.8 + matched-cell) |
| `.../WAVE_TWIN/WAVE_TWIN_REPORT.md` | `0c3b0d9` | clean | (§5 COMSOL 15.5%) |
| `.../JET_TEST/JET_TEST.json` | **`1bf29f1`** | clean (paper-1 frozen) | `_frozen_sources/JET_TEST.committed.json` |
| `.../SCRAMBLE_EXT/SCRAMBLE_RESULTS.json` | **`ed7a1e0`** | clean (paper-1 frozen) | `_frozen_sources/SCRAMBLE_RESULTS.committed.json` |
| `docs/ROUND63_GPT_ROUND45_RULING_RAW.md` | HEAD (`35d858e`) | tracked | (theorem statements + jet-flow eqs §B2/§B11) |

All cited artifact directories are **clean in the working tree** (`git status --porcelain` empty for
`R44_KILL_TESTS/`, `TLSG/`, `WAVE_TWIN/`, `JET_TEST/`, `SCRAMBLE_EXT/`).

---

## Theorem set (analytical — no data). R45 §C2 / §B

| # | Statement | Source | Location |
|---|---|---|---|
| T0-1 | **Statistical Noether Wall Theorem** — `P_{g·x}=P_x ∀g∈G ⇒ D_x log p[X_ξ]=0 a.s., I_x X_ξ=0`; integrable score-null distributions are local gauge leaves | R45 ruling | §B1.1 Thm 1, eqs (1.1)–(1.3); local converse §B1.1 |
| T0-2 | **Blindness Capacity Theorem** — `β_d(L)=σ_d(L)`; generalized `L†L v=λ B†B v`, `β_d=√λ_d`; robust `≤σ_d(L̂)+δ`; stacked-operator noncomposability `β_d²=λ_d^↑(Σ w_i² L_i†L_i)` | R45 ruling | §B4.1 Thm 2 eqs (4.1)–(4.4); §B4.3 (4.5)–(4.6) |
| T0-3 | **Three-Ledger Certifiability Theorem** — `known : blind : full = 1 : √p : p` in `T·λ`, `λ=½‖A_⊥‖_F²`; composed with quotient jet `T ≍ κ^-2 g^-2a ε^-2m` | R45 ruling | §B11 (LAN reduction 11.1–11.6; ledgers 11.2–11.5) |
| T0-4 | **Jet Transmutation** (corollary) — two-monomial flow `dm_eff/dlogε=2(m_eff−m1)(m2−m_eff)`; lower-order leak → `m=1`, projection/wall restores `m=2`; `z2,*∝ε^{-1/2}` | R45 ruling | §B2.1 eqs (2.4)–(2.11); §B2.3 named effect |

---

## §1 — Symmetry walls of an optical measurement

### A. Two exact walls (measured score-null orbits)

| # | Claim (frozen value) | File | Field | Commit |
|---|---|---|---|---|
| A1 | **Support (field) wall**: hard Fourier pupil makes the signed-intensity field exact to `2kR`; max field leak beyond `2kR` = **6.6×10⁻¹⁶** (machine zero), all pupil variants, ∀z1 | KT1 | `verdicts.P2a_field_wall.per_variant.V6_kR5.max_field_leak_beyond_2kR` = 6.591e-16 (V6_kR6 6.44e-16, V6_kR8 6.37e-16) | 19338a4 |
| A2 | ideal periodic-code detector wall machine-zero for `k>kp`, ∀z1 | KT1 | `verdicts.P1.value` = 2.070e-15 (thr 1e-12) | 19338a4 |
| A3 | **Energy (unitary) wall**: phase-only object change invisible to a lossless full-aperture bucket at **every order**; null = **2.18×10⁻¹⁶** ∀(ε,z2) | IT1 | `verdict.max_full_aperture_null` = 2.179772703870535e-16 | f13a323 |

### B. The two controlled breakers

| # | Claim (frozen value) | File | Field | Commit |
|---|---|---|---|---|
| B1 | **NA clipping breaks the energy wall at first order** (ε¹): clip-bucket slopes **0.97–1.03**, monotone in z2 | IT1 | `verdict.clip_slopes` = [1.033, 0.982, 0.970, 1.000]; `P_z2_trend`="MONOTONE" | f13a323 |
| B2 | **Finite object window breaks the support wall**: full-twin noncompact pedestal, detector leak beyond `2kp` = **5.3×10⁻²** (field 1.6×10⁻¹); frequency placement alone gives no wall | KT1 | `verdicts.P3.value` = 0.053154; `field_leak_beyond_2kp` = 0.16356 | 19338a4 |
| B3 | naive separable Tukey guard reduces the detector tail only ~8–40×, never to `1e-12` → an exact detector wall needs a concentration-grade (DPSS) guard | KT1 | `KT1_NOTE.md` (V7 guard 8–40×; 1.5e-3 kR5 / 1.6e-4 kR8) | 19338a4 |

---

## §2 — Measuring blindness capacity (Blindness Capacity Theorem, measured)

| # | Claim (frozen value) | File | Field | Commit |
|---|---|---|---|---|
| C1 | **Single-plane capacity**: at the fixed relay-conjugate plane (z1=0), **80 code dimensions ≤ 1×10⁻⁵** (and ≤1×10⁻⁴) — a measured blindness-capacity point | IT4b | `results.per_z1_vs_joint.z1_0_only_d_at_1e5` = 80; `z1_0_only_d_at_1e4` = 80 | 5666d1e |
| C2 | **Joint (four-plane) capacity**: **35 dimensions ≤ 1×10⁻⁴** jointly over z1∈{0,2,5,10} mm; d@1e-5 = 8; single-plane gain **+45** | IT4b | `results.joint_spectrum.{d_at_1e4=35, d_at_1e5=8, d_at_1e3=80}`; `per_z1_vs_joint.single_z1_gain`=45 | 5666d1e |
| C3 | full ascending generalized-leak spectrum `σ_d` vs `d` (d=1…120): `σ_1`=1.43e-6 → `σ_{64}`=4.99e-4 → `σ_{120}`=1.13e-1 (the plotted capacity curve) | IT4b | `results.joint_spectrum.rel_leak_vs_d` (30-point map) | 5666d1e |
| C4 | d=64 subspace: max leak **4.99×10⁻⁴**, operational RMS random-code leak **1.68×10⁻⁴** (3× below the worst-direction max); ridge-stable (d@1e-4=35 across ridge 0/1e-8/1e-6) | IT4b | `verdict.{joint_d64_max_leak=4.99e-4, joint_d64_rms_random_combo=1.68e-4}`; `results.ridge_sanity` | 5666d1e |
| C5 | certified-blind codes retain a second sensing channel: `λ_random/λ_blind` = **1.20** (≪3×), and 219× lower operational mean leak (2.58e-4 vs 5.64e-2) | IT4 | `measured.{lam_ratio_rand_over_blind=1.202, meanleak_blind_codes=2.58e-4, meanleak_random_codes=5.64e-2}` | 33b8ef2 |
| C6 | **fresh-draw robustness (TLSG bar 7)**: ≥64 conjugate-plane dims certified <1e-4 on fresh physical-law draws; min = **80** (calib d@1e-4 = d@1e-5 = 80) | TLSG | `bars[6].measured` ("fresh d@1e-4 min = 80"); arm B | 35d858e |
| C7 | **calibration robustness (TLSG bar 6, Eq. 4.4)**: fresh split, `σ_{64}(L̂)+δ` envelope, coverage = **1.000**; σ₆₄=3.4e-8, δ=1.07e-2, worst fresh S₆₄ leak 1.99e-6 (six orders below envelope) | TLSG | `bars[5].measured`; arm B | 35d858e |

---

## §3 — Restoring the desired jet (three routes + noncomposability)

| # | Claim (frozen value) | File | Field | Commit |
|---|---|---|---|---|
| D1 | **Route 1 — optical projection (jet transmutation)**: projecting out the measured mean-leak tangent restores the covariance-only slope to **4.00** in all four cells; unprojected total leak-dominated slope **2.12–2.20** (m=1); orthogonality residual **~1×10⁻¹⁸** | KT2 | `verdict.{projected_cov_slopes=[4,4,4,4], unprojected_lowEps_slopes=[2.136,2.2,2.117,2.176], orthogonality_residuals}` (max 5.09e-18, min 1.45e-18) | 05f1029 |
| D2 | generic-direction covariance control stays sub-quartic (slopes 1.58–2.28) | KT2 | `verdict.generic_cov_slopes` = [1.576, 2.284, 1.576, 2.285] | 05f1029 |
| D3 | **Route 2 — code-space (statistical projection)**: projecting the covariance score off the leak/law span drives nuisance overlaps to **0** (amplitude 0.36→0, medium 0.20→0) while retaining **86–99%** of target information | KT3 | `verdict.{clean_nuisance_overlap_amp_after=[0,0,0,0], clean_nuisance_overlap_medium_after=[0,0,0,0], retained_info_frac=[0.857,0.956,0.919,0.989]}` | 2430166 |
| D4 | **Route 3 — hardware DPSS concentration**: DPSS render+object window holds the beyond-2kp witness leak ≤ **2.4×10⁻⁴** uniformly in z1 (max 2.38e-4), 349× below no-window (8.3e-2); Tukey ~4e-2 (rim-kill) | IT5 | `verdict.{dpss_max_witness_leak=2.378e-4, none_max_witness_leak=8.30e-2, tukey_max_witness_leak=3.99e-2, dpss_improvement_over_none=349.29}` | c8a0487 |
| D5 | **Noncomposability**: DPSS does NOT stack with the SVD code route — single-z1 certified dimension drops **80→24** at 1e-4 (routes are alternatives, not factors) | IT5 | `combined_single_z1_SVD.{no_window_d_at_1e4=80, dpss_d_at_1e4=24}` | c8a0487 |
| D6 | Tukey-window noncomposability (2nd instance): apodized joint d@1e-4 = **24** (no-apod 35); Tukey **93× worse** at d=64 (rim effect) — apod rel-leak 4.66e-2 vs no-apod 4.99e-4 | IT4b | `results.window_apodized.{joint_apod_d_at_1e4=24, apod_rel_leak_at_d64=0.04655, noapod_rel_leak_at_d64=4.99e-4}` (ratio 93.3×) | 5666d1e |

---

## §4 — What finite data can certify (Three-Ledger + KT6 boundary)

**All from the frozen TLSG gate (`35d858e`) unless noted.**

| # | Claim (frozen value) | File | Field | Commit |
|---|---|---|---|---|
| E1 | **Ledger separation `1 : √p : p`**: fitted p-exponents **−0.004 / 0.444 / 1.010** (known / blind / estimation; frozen bars 0±0.15, ½±0.15, 1±0.15) | TLSG | `bars[3].measured` = "−0.004 / 0.444 / 1.010" | 35d858e |
| E2 | matched-score 50%-power contours collapse vs `T·λ`: CV = **0.125** (≤0.20) | TLSG | `bars[0].measured` = "CV = 0.1250" | 35d858e |
| E3 | blind contours collapse vs `T·λ/√p` (CV = **0.124**) AND fail vs `T·λ` (CV = **0.628**) | TLSG | `bars[1].measured` = "CV_scaled=0.1242; CV_raw=0.6283" | 35d858e |
| E4 | estimation risk collapses vs `T·λ/p`: CV = **0.035**, mean risk·`T·λ/p` = **0.985** (≈ CRB `p/T`) | TLSG | `bars[2].measured` = "CV=0.0348 (mean risk·T·λ/p=0.985)" | 35d858e |
| E5 | **KT6 boundary**: the near-contact operating point sits **below** optimal blind certifiability — blind power = **0.051** at FPR 0.05 (10×: 0.056; 100×: 0.134); p=136, `T·λ`=0.0866 | TLSG | `bars[4].measured` | 35d858e |
| E6 | KT6 oracle law confirmed: near-contact `J_cov/J_mean ∝ z2⁴`, exponents **3.99–4.08**; but blind signal/floor ≈ **5.4×10⁻³ ≪ 1** (below the Wishart `M²/T` floor) → ESTIMATOR_KILLED | KT6 | `verdict.{oracle_near_contact_exponents=[4.081,4.054,4.022,3.989], blind_median_signal_to_floor_near≈0.0054, worst=5.25e-3, P_blind_estimable="ESTIMATOR_KILLED"}` | 635ff05 |
| E7 | **projected jet above its channel line** (hero right panel): the projected m=2 covariance channel (D1, slope 4.00) is a known-direction (Ledger I) observable — contrasts with KT6 below the blind line | KT2 + TLSG | KT2 `projected_cov_slopes`=[4,4,4,4] (05f1029); ledger geometry TLSG (35d858e) | 05f1029 / 35d858e |

---

## §5 — Wave-optics transfer and apparatus constraints

| # | Claim (frozen value) | File | Field | Commit |
|---|---|---|---|---|
| F1 | **Calibrated leak law** (End-Matter-C wording, R45 §C7): `L(z1)=6.94×10⁻³ + 5.92×10⁻³·(z1/mm)^{0.90}`; floor, coefficient, exponent stable across code seeds & counts with CVs **2.2% / 2.2% / 0.1%** | IT3 | `verdict.{L_pix_mean=6.94e-3, a_mean=5.92e-3, p_mean=0.9036, L_pix_CV=0.0219, a_CV=0.022, p_CV=0.001}` | 6fb5149 |
| F2 | conjugate-plane floor = fill-factor + macro-pixel staircase (instrument constant `L_pix`=7.15e-3): floor moves 6% (fill), 50% (macro-staircase); pupil does not reduce it (kR=8 → 1.29e-1) | IT2 | `verdict.{floor_move_fill=0.056, floor_move_macro=0.501, L_pix=7.15e-3, pupil_reduces_floor=false}` | acfcda0 |
| F3 | **COMSOL full-Maxwell validation**: thin-phase-screen speckle contrast **0.467** vs FEM (COMSOL 6.3 Wave Optics, 2D TE) **0.395** → **15.5%** agreement; VALIDATED (developed fine-grain cell) | WAVE_TWIN | `WAVE_TWIN_REPORT.md` §9 table (0.467 / 0.395 / 15.5%); `COMSOL_MICRO/` executed | 0c3b0d9 |
| F4 | **Code-count law**: covariance `λ ∝ M^{1.8–1.9}` (unsaturated at finite z2); fitted exponent **1.78** (mid z2=5mm, 2%) / **1.89** (near-contact) | WAVE_TWIN | `T5b_MSCALING.json` `rows[*].scaling_exponent_p` = 1.78 / 1.76 / 1.89 | 0c3b0d9 |
| F5 | **Matched-cell transfer**: M=128 extrapolated `T_det` ≈ **768 banks @2%** for the geometry-matched developed cell — inside the sealed band **[513, 2996] banks** | WAVE_TWIN | `T5b_MSCALING.json` `rows[0].extrap_Tdet_M128`=767.7; `matched_sealed_2pct_range`=[513,2996] | 0c3b0d9 |
| F6 | **Mean-wall leak is nearly-but-not-exactly blind on a real DMD**: deterministic beyond-band leak grows 6.8e-3 (z1=0 pixelation floor) → 5.4e-2 (z1=10mm) → 8.9e-2 (z1=20mm) | WAVE_TWIN | `T1_WALL_LEAK.json` `{conjugate_z1_0_pixelation_floor=6.81e-3, nominal_z1_10mm_leak=5.43e-2, ...z1=20mm 8.92e-2}` | 0c3b0d9 |
| F7 | **PWM apparatus constraint** (R45 §C7): binary rendering must use full-period PWM + whole-period bucket integration; Bayer **ordered spatial dither doubles** the conjugate-plane leak (7.0e-3 → **1.38×10⁻²** ≥ 1e-2) | IT8 | `verdict.{dither_leak_z1_0=1.384e-2, analog_leak_z1_0=7.01e-3}` | 6fb5149 |

---

## Supplement / future objects (R45 §C2, §C6). Theorem-consistent negatives.

| # | Claim (frozen value) | File | Field | Commit |
|---|---|---|---|---|
| G1 | **Segmented-bucket negative** (IT6): the `P²` law is real (exp 1.68–2.02, P=4 buys 12–17×) but segmenting **breaks the mean wall** (per-segment d′ 0.8→1.8 at z2=5) — single-bucket stands | IT6 | `IT6_segmented_bucket_p2.json` verdict (KILL) | dade877 |
| G2 | **Coincidence-veto negative** (IT7, §C6): contact arm exactly medium-blind (AUC **0.5009**), medium FA repaired 17.8× (0.9965→0.056), but only **24%** scene power retained — invalid beyond-band repair | IT7 | `measured.{contact_AUC_H0_vs_medium=0.5009, single_arm...FA_medium=0.9965}` | e1546ac |
| G3 | task-dependent étendue law (future): `N_* = [ν/(2−ν)] n`; KT5 N_eff grows 46→811, T_det∝N_eff (saturated CV 0.25), universal collapse partial (normCV 0.64) | KT5 | `KT5_etendue_neff_collapse.json` verdict (PARTIAL) | 752114a |
| G4 | coherence-rank unsaturation (future, COHERENCE-CROSSOVER gate): rank ≤ min(N,R²) holds ∀z2; thin screen stays R≈110 (no R→1 collapse) → memory-effect phase diagram deferred | KT4 | `KT4_coherence_rank_sketch.json` verdict (PARTIAL) | 494b46b |

---

## Cross-paper anchors (JET_TEST `1bf29f1` / SCRAMBLE `ed7a1e0`) — cited where the jet/scrambling law is invoked

| # | Claim (frozen value) | File | Field | Commit |
|---|---|---|---|---|
| H1 | integer contact orders exact: KL slopes **2.038 / 4.000** (generic / first-order-orthogonal); coef ratio 1.07% / 0.00% | JET_TEST | `bank_A.per_dir.{delta_g.kl_slope=2.0381, delta_o.kl_slope=3.9999}` | 1bf29f1 |
| H2 | gauge collapse `J_θ=0` → finite-prior anchor recovers 20% (amplitude-anchor condition) | JET_TEST | `bank_B.{single_zero_lag.I_theta_eff=0.0, amplitude_anchor.frac_recovered=0.20}` | 1bf29f1 |
| H3 | iso-Q cancellation → AUC ≈ chance at ε*; no broad in-aperture blind set (frac_blind=0) | JET_TEST | `bank_C.cancellation.{eps_star=0.3626, auc_at_star=0.4807}`; `kill_scan.frac_blind=0.0` | 1bf29f1 |
| H4 | complete-scrambling law `Cov(b)=Q(x)(O∘O)+R`, `Q=xᵀGx`; 99.99% covariance energy in one direction (rank-one) | SCRAMBLE | `valB_rank1.frac_cov_energy_in_OhO_direction`=0.99992 | ed7a1e0 |

---

## Shared model configuration (Methods / captions)

| # | Value | File | Field | Commit |
|---|---|---|---|---|
| I1 | KT-battery grid 2048 (IT1/KT1), complex128 / float64; M∈{8,12,16,32}, n_bank 200–700 per test | KT1/KT2/IT1 | `params` blocks | (per test) |
| I2 | TLSG arm A: M∈{8,16,32,64} × p_eff∈{10,36,136,528}, λ=0.15, N_MC=1500, FPR=0.05, seed 20260724 | TLSG | `partA_params` | 35d858e |
| I3 | TLSG arm B: d_cert=64, n_cal=n_test=24, conjugate z1=0, perturb fill±5%/pitch±3%/relay defocus 0–0.5mm | TLSG | `partB_params` | 35d858e |

---

# GAPS AND DISCREPANCIES (flagged, NOT regenerated)

## GAP-1 — Non-blocking: R45 §C3 "7e-16" field wall vs committed 6.6e-16
The R45 §C1/§C3 hero text names the hard-pupil support wall as `7e-16`. The **committed** KT1
summary statistic (`verdicts.P2a_field_wall`) is the max field leak beyond `2kR` over the three
hard-pupil variants = **6.591×10⁻¹⁶** (largest single per-cell field leak across all variants/z1 is
7.70×10⁻¹⁶). Both are machine-zero. **Resolution:** figures and text cite the committed **6.6×10⁻¹⁶**
(A1); R45's `7e-16` is a loose round of the same machine-zero wall. No regeneration; no science change.

## GAP-2 — Non-blocking: hero "35@1e-4 joint" reads a d@1e-4 spectrum point, not a fresh-draw min
R45 §C3 center panel marks "35 dimensions ≤1e-4 jointly over four z1 planes." The committed source
is the **single deterministic joint spectrum** (`IT4b joint_spectrum.d_at_1e4=35`, C2), not a
fresh-draw ensemble. The **fresh-draw** certification (TLSG bar 7, C6) is single-plane (min 80). Both
are cited with their correct provenance; the hero labels the joint-35 point as a spectrum value and
the 80 point as the fresh-draw / single-plane capacity. No conflict, flagged for caption precision.

## GAP-3 — Non-blocking: "80→24" noncomposability has two committed instances
The 80→24 single-plane capacity collapse is committed in **two** places: IT5
`combined_single_z1_SVD` (SVD 80 → DPSS 24, D5) and IT4b `window_apodized` (joint 35 → Tukey-apod 24,
D6). The hero/§3 cite the IT5 single-z1 SVD→DPSS instance (matches R45 §C3 "DPSS route ≤2.4e-4 but
joint capacity 80→24"); the IT4b Tukey instance (93× rim-kill, D6) is the second, independent
noncomposability witness. Recorded so the two "→24" numbers are not conflated.

## GAP-4 — Non-blocking: Optica length target vs guideline
R45 sets no page budget; the task targets ~10–12 pp main. Optica **Research Article** guideline is
**6–8 pages** (Letters 4 pp; mini-reviews 8–12 pp), abstract **~100 words** (see BUILD_STATUS.md).
Flagged as a length-budget decision for the writing phase, not an evidence gap.

**No numerical claim in the Paper-2 architecture (R45 §C2 theorem set, §C3 sections/hero, §C7
transfer sentences) lacks a committed source.** All four GAPs are non-blocking (rounding / caption
precision / length budget); none requires data generation.
