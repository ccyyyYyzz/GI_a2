# R23 bridge-theorem falsifier regression — POST-HOC DESCRIPTIVE

**Label: POST-HOC DESCRIPTIVE. Zero new simulation. Not preregistered, not a
verdict input.** Computed 2026-07-20 on the frozen M1 campaign data
(m1-freeze @ 6f00932, R19-corrected outcomes). Authority for the functionals:
`docs/ROUND63_GPT_ROUND23_RULING_RAW.md` §2.5 (T1.7–T1.10), §3.3 (T2.4),
§5 (T4.1–T4.5). Raw numbers: `FALSIFIER_REGRESSION.json` (same directory).

**R24 context note (2026-07-22).** After this analysis was specified, the R24
ruling (`docs/ROUND63_GPT_ROUND24_RULING_RAW.md`, issue #16) redefined the
deployment gate as a **robust, estimator-aware finite-library router with
abstention (RGLI)** — not the raw R23 residual. This regression therefore
stands as the M1-data validation of the *raw* alignment-diagnostic concept;
the deployed gate will use the richer R24 robust form. The degeneracies found
below (branch-edge collapse of the raw functionals) are directly relevant
input to that re-architecture.

---

## 1. Verdict summary

| prediction | statement | verdict |
|---|---|---|
| **P1** | \|corr(dQ, Align1_exact)\| > \|corr(dQ, C_u)\| | **FAIL** — \|r\| = 0.16 (seed 0) / 0.20 (5-seed mean) vs 0.64 for C_u |
| **P2** | per-scene cert gap correlates positively with misalignment (eps1 or M1sq_local) | **NOMINAL PASS, CONFOUNDED** — Pearson +0.50/+0.51, Spearman +0.54/+0.55; but partial correlation given the xhat-predicted contrast Cu_t4 collapses to +0.08 / −0.16 |
| **P3** | contour scenes occupy the top misalignment ranks; load dispersion adds explanatory power for contour | **FAIL** — contour ranks 10–16/24 by eps1 (5–17 by M1sq), mid-pack; q95/q50 adds ΔR² ≈ +0.06 overall and does not separate contour_3 (+4.15 dB) from the negative contour trio |

**Bottom line (by the ruling's own falsifier criterion, §8):** at the M1 ridge
operating point, the measured leverage–brightness slope and orthogonal
residual do **not** predict the observed per-scene need for more than one gain
profile. The raw R23 bridge functionals fail on this data — but with a
structural reason (§5 below) that the ruling itself anticipated: the
production ridge sits at the finite-nu kernel argmax, **outside the declared
monotone-concave branch**, where alpha* diverges and the alignment residual
degenerates. The negative result is a scope finding, not merely noise.

---

## 2. Methods — exact code paths

All quantities were rebuilt with the FROZEN campaign machinery, run from repo
root with `code/` and `code/round63/` on `sys.path`, `M1_FREEZE_LAUNCHED=1`
(post-unblinding descriptive use of the confirmatory pre-scan path — noted per
the R17 §5 rule; the campaign is unblinded and concluded).

1. **Deployed design**: `m1_runner.deployed_scat32()` — the frozen balanced
   972-row SCAT32 multiset, sha `ff133328f00b0867` (from
   `results/round63_m1/designs/scat32_deployed.npz`).
2. **Pre-scan estimate xhat**: `m1_runner.prescan_estimate(x_true, image,
   seed, 0.60, 2000)` — the I2 common per-(image, seed) realization at the
   (0.60, nu=2000) anchor, stream `rng_for(seed, 63, 9, img)`; truth of record
   via `campaign._images(32, "all", imageset="m1")` (`data/r63_images_m1/`).
   Primary = seed 0; all diagnostics also averaged over the 5 frozen seeds.
3. **Ridge loads rho_i**: inline replica of
   `m1_runner.ridge_scat32_calibration` at nu=2000 (identical arithmetic,
   **no cache write** into the frozen `designs/` dir): global multiplier
   kappa = rho_R / mean(rows @ xhat) with rho_R(2000) = 22.2545 from
   `oed_design_v4.ridge_target4(2000)` (clip reason NONE), then the frozen
   `oed_design_v5._kernel_grids` guard bisection (never triggered: 0/24
   scenes clipped; achieved mean load = 22.2545 in all scenes).
   Verified to reproduce the cached
   `results/round63_m1/designs/ridge_scat32_m1_*_0.npz` calibrations.
4. **Kernel J_nu(rho)**: `oed_design_v3.kernel_eval(rho, 2000)["J_exact"]` —
   the finite-nu exact active-start non-paralyzable renewal kernel. This is
   the object the certificate uses (`oed_design_v4._J` is the same quantity
   through the committed table `results/round63_theory/fig_a_sweep.npz`;
   agreement checked < 1e-5 relative at rho ∈ {0.5, 5, 22}).
   Derivatives: 500-node log-grid on rho ∈ [1e-3, 90] (dt = 0.0229 in
   ln rho), J′ = np.gradient in t = ln rho divided by rho; **alpha\*** uses a
   direct 5-point central stencil on the exact kernel with h = 2e-3 in
   ln rho (documented; the grid gradient is too coarse at the ridge where
   J′ → 0). Finite-nu argmax: rho ≈ 22.32 (grid) / 22.21 (dense scan) — the
   production ridge target 22.2545 sits essentially **at** the argmax.
5. **Leverages ell_i**: exactly the `oed_design_v5.cert_deployed_rows`
   assembly (V0 + per-row rank-one updates): task subspace B (r=200) and
   eps0 from `oed_design_v4.subspace_from_fixedstar(info_matrix_full(rows,
   xhat, 2000, rho_bar, P=balanced_prescan_52()))` (the I4 frozen rule);
   V0 = `oed_design_v4.V0_prescan(P52, xhat, 2000, rho_prescan, B)`;
   V = V0 + Σ_i 2000·J(rho_i)·q_i q_iᵀ + eps0·I with
   q_i = Bᵀ a_i/(a_i·xhat); ell_i = q_iᵀ V⁻¹ q_i. The certificate's
   sensitivities are d_i ∝ nu·J(rho_i)·ell_i; the T1 marginal value used here
   is h_i = ell_i·J′(rho_i).
6. **Which-V0 ambiguity — resolved by sensitivity, not by guessing.** The
   certificate freezes the subspace at the cell's own budget corner (I4), but
   the ridge deployment has no single frozen corner. Three variants were run:
   (a) **primary** `subspace+V0 at rho_bar=0.60` (the SCAT32-060 cert corner,
   the arm dQ is measured against), (b) `0.05` corner, (c) subspace at the
   ridge multiplier itself with V0 at 0.60. **All correlations agree to
   < 0.05** (see JSON `variants`), so no STOP was needed — the choice is
   immaterial at this operating point.
7. **Outcomes (frozen, not recomputed)**: dQ_j from
   `results/round63_m1/CORRECTION_2026-07-19/M1_VERDICTS_SPEC_CORRECTED_R19.json`
   → `RIDGE_OPERATING.per_image_delta` (24 scenes); per-scene mean
   certificate counterexample gap = mean `primal_gap_lower_per_r` over the
   299 COUNTEREXAMPLE rows of `results/round63_m1/shards/M1_CERT_*.csv`,
   parsed with the `m1_score.py` overflow-merge reader (reader logic inlined
   — importing `m1_score` executes its unblinding script); campaign C_u from
   the frozen `M1_SCAT32-060_00.csv` shard (scene-constant column).
   Baseline validation: corr(dQ, C_u) = **+0.6392** Pearson — reproduces the
   +0.64 reference exactly.

### Diagnostics per scene (ruling references)

- **eps1 / Align1** (T1.7–T1.10, T2.4): h_i = ell_i·J′(rho_i);
  s = rho/kappa (the deployed load profile = rows @ xhat);
  eps1 = ‖h − (hᵀs/sᵀs)s‖₂; Align1 = 1 − eps1²/‖h‖².
  **Budget-interior form: projection onto span{s}^⊥ only — the dose-band /
  peak cone (the full NNLS in T1.7) is ignored in this first pass, as
  specified.** The T1.10 constant-vector variant (span{1}) is also reported
  (`eps1_const`, `Align1_const`).
- **T4 local functional** (T4.5): b_i = a_i·xhat, z_i = b_i/mean(b) − 1,
  C_u^{T4} = std(b)/mean(b) (note: this is the *xhat-predicted* contrast,
  distinct from the campaign C_u which uses the true image);
  eta_i = ell_i/mean(ell) − 1; alpha = OLS slope of eta on z;
  sigma_zeta = std of residuals; alpha\* = −rhobar·J″(rhobar)/J′(rhobar) at
  the achieved mean load (unclipped, as instructed);
  M1sq_local = (alpha − alpha\*)²·C_u² + sigma_zeta².
- **Load dispersion**: q95/q50 of rho_i; fraction of rho_i beyond the
  finite-nu argmax (22.32).

Runtime: 315 s total (kernel grid 21 s, 120 pre-scans 1 s, 3 variants × 24
scenes × leverage assembly). Script: session scratchpad `falsifier.py`
(analysis-only; writes nothing outside `results/round63_next/`).

---

## 3. Per-scene table (primary variant, seed 0; dQ/C_u/gap are frozen outcomes)

| image | family | dQ (dB) | C_u | cert gap | eps1 ×1e5 | Align1 | Cu_t4 | alpha | alpha\* | sigma_zeta | M1sq ×1e-10 | q95/q50 | frac>argmax |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| m1_chirp_0 | chirp | +2.36 | 0.215 | 1.583 | 2.15 | 0.018 | 0.289 | −2.01 | −7.10e5 | 0.57 | 4.22 | 1.51 | 0.47 |
| m1_chirp_1 | chirp | +1.27 | 0.267 | 1.606 | 1.69 | 0.020 | 0.270 | −1.91 | −7.10e5 | 0.53 | 3.68 | 1.50 | 0.47 |
| m1_chirp_2 | chirp | +1.01 | 0.253 | 1.384 | 1.99 | 0.023 | 0.302 | −1.94 | −7.10e5 | 0.56 | 4.60 | 1.57 | 0.46 |
| m1_chirp_3 | chirp | +1.54 | 0.234 | 1.310 | 0.86 | 0.023 | 0.224 | −1.54 | −7.10e5 | 0.46 | 2.53 | 1.40 | 0.48 |
| m1_contour_0 | contour | **−3.61** | 0.213 | 1.239 | 1.08 | 0.021 | 0.239 | −1.59 | −7.11e5 | 0.46 | 2.89 | 1.41 | 0.48 |
| m1_contour_1 | contour | **−5.79** | 0.217 | 1.347 | 1.18 | 0.026 | 0.254 | −1.72 | −7.10e5 | 0.50 | 3.25 | 1.49 | 0.47 |
| m1_contour_2 | contour | **−4.74** | 0.245 | 1.674 | 1.24 | 0.028 | 0.268 | −1.72 | −7.10e5 | 0.48 | 3.63 | 1.51 | 0.46 |
| m1_contour_3 | contour | +4.15 | 0.331 | 1.626 | 1.49 | 0.024 | 0.276 | −1.78 | −7.10e5 | 0.52 | 3.85 | 1.50 | 0.46 |
| m1_glyph_0 | glyph | +4.73 | 0.371 | 1.390 | 1.63 | 0.017 | 0.260 | −1.82 | −7.10e5 | 0.51 | 3.42 | 1.49 | 0.48 |
| m1_glyph_1 | glyph | +6.01 | 0.363 | 1.416 | 1.50 | 0.019 | 0.252 | −1.86 | −7.10e5 | 0.49 | 3.20 | 1.43 | 0.49 |
| m1_glyph_2 | glyph | +2.43 | 0.292 | 1.382 | 1.17 | 0.027 | 0.261 | −1.70 | −7.10e5 | 0.49 | 3.43 | 1.44 | 0.47 |
| m1_glyph_3 | glyph | +2.21 | 0.291 | 1.418 | 1.62 | 0.018 | 0.267 | −1.80 | −7.10e5 | 0.51 | 3.58 | 1.49 | 0.46 |
| m1_maze_0 | maze | +6.33 | 0.339 | 1.262 | 1.10 | 0.030 | 0.259 | −1.68 | −7.10e5 | 0.48 | 3.40 | 1.50 | 0.45 |
| m1_maze_1 | maze | +6.44 | 0.528 | 1.306 | 0.77 | 0.028 | 0.229 | −1.44 | −7.10e5 | 0.42 | 2.64 | 1.44 | 0.47 |
| m1_maze_2 | maze | +5.58 | 0.461 | 1.186 | 0.82 | 0.019 | 0.219 | −1.35 | −7.10e5 | 0.41 | 2.42 | 1.39 | 0.48 |
| m1_maze_3 | maze | +6.41 | 0.458 | 1.241 | 1.05 | 0.026 | 0.246 | −1.63 | −7.10e5 | 0.45 | 3.06 | 1.47 | 0.47 |
| m1_microtexture_0 | microtexture | +1.53 | 0.197 | 1.374 | 1.47 | 0.025 | 0.276 | −1.81 | −7.10e5 | 0.50 | 3.83 | 1.51 | 0.46 |
| m1_microtexture_1 | microtexture | +1.25 | 0.177 | 1.366 | 1.84 | 0.019 | 0.279 | −1.91 | −7.10e5 | 0.52 | 3.92 | 1.54 | 0.46 |
| m1_microtexture_2 | microtexture | +1.53 | 0.200 | 1.406 | 1.58 | 0.019 | 0.271 | −1.74 | −7.10e5 | 0.50 | 3.71 | 1.49 | 0.47 |
| m1_microtexture_3 | microtexture | +1.45 | 0.197 | 1.462 | 1.61 | 0.027 | 0.287 | −1.93 | −7.10e5 | 0.57 | 4.17 | 1.57 | 0.46 |
| m1_spokes_0 | spokes | −0.77 | 0.296 | 1.381 | 0.79 | 0.024 | 0.226 | −1.43 | −7.10e5 | 0.43 | 2.57 | 1.40 | 0.47 |
| m1_spokes_1 | spokes | **−5.08** | 0.297 | 1.333 | 0.88 | 0.025 | 0.234 | −1.50 | −7.10e5 | 0.43 | 2.76 | 1.42 | 0.47 |
| m1_spokes_2 | spokes | +4.53 | 0.350 | 1.443 | 1.02 | 0.017 | 0.229 | −1.51 | −7.10e5 | 0.44 | 2.65 | 1.41 | 0.47 |
| m1_spokes_3 | spokes | +2.20 | 0.309 | 1.448 | 0.97 | 0.025 | 0.239 | −1.61 | −7.10e5 | 0.46 | 2.87 | 1.42 | 0.48 |

(All 24 scenes: achieved mean load 22.2545, guard clip never active. cert gap
= mean over that scene's COUNTEREXAMPLE cells, n = 8–19 of 20 per scene, 299
total. Full precision + 5-seed means + eps1_const columns in the JSON.)

---

## 4. Correlations and prediction verdicts (n = 24 — CIs are wide; a Pearson r
of ±0.40 has a 95% CI spanning roughly ±0.38 around it; these are
descriptive, not inferential)

### P1 — Align1 vs C_u as a dQ predictor: **FAIL**

| | Pearson | Spearman |
|---|---:|---:|
| corr(dQ, Align1_exact), seed 0 | −0.156 | −0.111 |
| corr(dQ, Align1_exact), 5-seed mean | −0.198 | −0.123 |
| corr(dQ, C_u) (baseline) | **+0.639** | **+0.701** |

Align1 is uniformly tiny (0.017–0.030 across all 24 scenes — h is ~98%
orthogonal to the rank-one direction everywhere) and carries essentially no
between-scene signal. The span{1} variant (T1.10 exact corollary form) is
equally uninformative. C_u remains by far the stronger predictor.

### P2 — cert gap vs misalignment: **NOMINAL PASS, CONFOUNDED**

| | Pearson | Spearman |
|---|---:|---:|
| corr(cert gap, eps1), seed 0 | +0.501 | +0.542 |
| corr(cert gap, M1sq_local), seed 0 | +0.514 | +0.545 |
| (5-seed means: +0.527 / +0.539 Pearson) | | |

The signs are as predicted. But: corr(M1sq_local, Cu_t4) = **+0.999** and
corr(eps1, Cu_t4) = +0.907 — at this operating point both misalignment
functionals are near-deterministic transforms of the xhat-predicted bucket
contrast (because alpha\* is scene-constant to 0.04% and dominates alpha, and
h's magnitude profile is contrast-driven). Partial correlations given Cu_t4:
eps1 → **+0.08**, M1sq → **−0.16**. The honest reading: the cert gap
correlates with predicted bucket contrast (+0.52); the "misalignment" framing
adds nothing beyond that on this data. (Curiosity flag: Cu_t4 is
*anti*-correlated with the true-image campaign C_u, r = −0.54 — the 4×4-block
pre-scan smoothing inverts the contrast ordering; and corr(cert gap,
campaign C_u) = −0.34.)

### P3 — contour at top misalignment ranks + load-dispersion power: **FAIL**

- Contour misalignment ranks (1 = most misaligned of 24): eps1 → 16, 13, 12,
  10; M1sq → 17, 14, 9, 5. Mid-pack, not top. The most misaligned scenes are
  chirp_0/2 and microtexture_1 — all solidly dQ-positive.
- Within contour, misalignment points the wrong way: contour_3 (the +4.15 dB
  outlier) has *higher* eps1 and M1sq than the three negative contour scenes.
- Load-dispersion covariates: R²(dQ ~ C_u) = 0.409; adding q95/q50 →
  0.468 (ΔR² = +0.059); adding frac>argmax → 0.419 (ΔR² = +0.010). Modest at
  best, and the load-dispersion columns are themselves nearly scene-constant
  (q95/q50 ∈ [1.39, 1.57]; frac>argmax ∈ [0.45, 0.49]) because every scene is
  calibrated to the same mean load through near-identical smoothed xhat load
  profiles.

---

## 5. Why the raw functionals degenerate here (structural finding)

The production ridge target rho_R(2000) = 22.2545 sits **at the finite-nu
kernel argmax** (dense-scan argmax ≈ 22.21–22.32, clip reason NONE — the
guards never bound). Consequences, all anticipated by the ruling's own
fences:

1. **J′(rhobar) ≈ 0 and changes sign across the pattern population**: every
   scene deploys ~46–49% of its 972 loads beyond the argmax. h = ell·J′ is a
   mixed-sign, near-zero-sum vector, hence almost orthogonal to the positive
   load profile s — Align1 ≈ 0.02 for every scene, with no discriminative
   variance. The R23 branch assumption (declared compact branch with J′ ≥ 0,
   J″ < 0) is **violated at the deployed operating point**; T1 concavity and
   the T1.6 gap bound do not formally apply there (ruling §6.1).
2. **alpha\* diverges**: alpha\* = −rhobar·J″/J′ = −7.10e5 (sign and
   magnitude are numerically fragile since J′ ≈ 0 — reported unclipped as
   instructed). Ruling §5: "alpha_J\* diverges as the ridge is approached.
   This makes one-knob alignment increasingly fragile near the cap." At the
   cap itself, (alpha − alpha\*)² is a scene-constant ~5e11 and M1sq_local
   collapses to alpha\*²·C_u², i.e. a monotone transform of predicted
   contrast — the two-factor structure (slope mismatch + orthogonal residual)
   is unobservable here.
3. Therefore this dataset **cannot cleanly test the bridge in its intended
   regime** (interior of the monotone-concave branch). What it does show: at
   the operating point the campaign actually used, the raw diagnostics
   reduce to contrast proxies and add no predictive power over C_u — which
   is itself the strongest argument for the R24 re-architecture (robust
   gate, abstention, estimator-aware routing) over raw-residual gating.

---

## 6. Honest caveats

- **n = 24 scenes**; all correlations are descriptive with wide CIs. No
  multiplicity control; three predictions × several functionals were
  examined. Nothing here is a preregistered endpoint.
- **POST-HOC**: computed after unblinding, on outcomes already known
  (dQ = R19-corrected, cert gaps descriptive under FALLBACK_DESCRIPTIVE).
  The confirmatory pre-scan path was re-invoked descriptively under
  `M1_FREEZE_LAUNCHED=1` (deterministic streams — identical realizations to
  the campaign's, no new randomness).
- **Machinery approximations, as specified**: (i) eps1 uses the
  budget-interior span{s} projection, ignoring the dose/peak normal cone
  (full T1.7 NNLS not attempted in this pass); (ii) J′/J″ by numerical
  differentiation (500-node log-grid; 5-point stencil h = 2e-3 for alpha\*)
  of the exact kernel, not closed form; (iii) leverage/subspace at the frozen
  0.60 cert corner (primary) — but variant sensitivity shows < 0.05 movement
  in every correlation, so the which-V0 question did not require a STOP;
  (iv) per-scene cert gap averages over that scene's own COUNTEREXAMPLE cell
  mix (8–19 of 20 cells; the NUMERICAL_UNRESOLVED remainder carries no gap
  value) — mixing nu ∈ {200, 2000} and b ∈ {0.05, 0.60} anchors, whereas the
  alignment diagnostics are at the ridge/nu=2000 point only.
- **The T4 C_u is the xhat-predicted contrast**, not the campaign's
  true-image C_u; they are anti-correlated here (r = −0.54) due to pre-scan
  smoothing. Statements about "contrast" must say which one.
- The five-seed averaging changes nothing (pre-scan estimates are
  low-variance after 4×4 smoothing; all seed-0 vs seed-mean correlations
  agree to < 0.05).
