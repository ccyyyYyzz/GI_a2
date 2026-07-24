# R44 Kill-Test Battery — index & verdicts

Bench-phase / next-paper kill tests for GPT-Pro R44 (`docs/ROUND63_GPT_ROUND44_RULING_RAW.md`)
and the internal divergence R2 (`results/round63_next/DIVERGENCE_R2/`). Run on Colab L4
(pro2), float64 for machine-zero walls / float64 covariance algebra. **Quarantined**:
nothing here touches `paper_prl/` or the sealed artifacts. Each test = standalone
`.py` + `<name>.json` (parameters, frozen predictions, measured values, verdict) +
`<name>_NOTE.md`.

## Wave-1 (R44 KT1–KT6) — COMPLETE

| Test | R44 | Score | Verdict |
|---|---|---|---|
| **KT1** shell sweep | §1.4 / #1 | 125 | **PASS** — ideal wall exact (2e-15); hard pupil = exact FIELD wall to 2kR (7e-16 ∀z1); finite object window reintroduces a z1-independent detector tail (naive Tukey guard 8–40× only; DPSS needed = IT5); full twin noncompact pedestal → frequency placement alone gives no wall |
| **KT2** jet transmutation | §5.2 / #2 | 125 | **PASS** — first-order-null direction exists (residual 1e-18); projected cov slope **4.00**, unprojected 2.12–2.20 (leak-dominated m=1), generic sub-quartic; mean-leak projection restores m=2 |
| **KT3** leak-orthogonalized sentinel | §5.6 / #6 | 100 | **PASS** — projection drives nuisance covariance overlaps to **0** (amp 0.36→0, medium 0.20→0), target info retained 86–99%; operational AUC honestly estimator-limited (oracle-W anti-orders at small λ) |
| **KT4** coherence-rank sketch | §5.3 / #3 | 125 | **PARTIAL** — rank ≤ min(N,R²) holds ∀z2 + orthogonal jet slope stays 4.0 (Jet-Order Separation); BUT bound non-binding (rank capped by code count M(M+1)/2, not R²≈12000) and R stays ~110 (no R→1 collapse on a thin screen) |
| **KT5** étendue/N_eff | §5.4 / #4 | 100 | **PARTIAL** — N_eff grows with aperture 46→811 (PASS); T_det ∝ N_eff (saturated CV 0.25); universal collapse vs n/N_eff partial (normCV 0.64) |
| **KT6** channel-ratio ranging | §5.5 / #5 | 80 | **ORACLE-PASS / BLIND-ESTIMATOR_KILLED** — oracle z2⁴ law confirmed (exponents 3.99–4.08); blind λ_cov ~185× under the Wishart M²/T floor (signal/floor 5.4e-3 ≪ 1) — binding caveat honored |

## Wave-2 (internal divergence queue) — COMPLETE (IT1–IT9)

| Test | Source | Verdict |
|---|---|---|
| **IT1** two-wall / jet invariance | higher #1 | **PASS — SECOND EXACT WALL** (energy-conservation type): phase-only object change invisible to a lossless full-aperture bucket at every order (null 2.18e-16 ∀z2); NA-clipped bucket responds at ε¹ (slopes 0.97–1.03), monotone in z2 |
| **IT2** floor decomposition | nextpaper #9 | **DECOMPOSES** — z1=0 floor moves 50% under macro-staircase, 6% under fill → L_pix=7.2e-3 real instrument constant; leak law L_pix + 6.1e-3·z1^0.91 (near-linear) |
| **IT3** leak-law stability | stranger #8 | **REPRODUCIBLE** — (L_pix,a,p) stable to CV 2.2/2.2/0.1% across 3 seeds × M{8,32}; leak law 6.94e-3 + 5.92e-3·z1^0.90 quotable |
| **IT4** SVD certified-blind codes | finer #3 (sledgehammer gate) | **PARTIAL/near-PASS** — mean-blind subspace exists, covariance retained (λ ratio 1.20 ≪ 3×), 219× lower operational mean leak; misses strict 64-dim ≤1e-4 gate by ~5× (a ~48-dim subspace clears it) |
| **IT4b** dimension/leak trade | coordinator follow-up | joint d@1e-4=35 (fundamental); d=64 RMS random-combo 1.68e-4 (3× below gated max); **single-z1 relay-conjugate d@1e-4=80** (+45); Tukey window is a KILL (rim effect); ridge-stable |
| **IT5** DPSS concentration window | finer #4 | **PASS** — DPSS render+object window holds beyond-2kp witness leak ≤2.4e-4 uniformly in z1 (349× below no-window; Tukey ~4e-2 rim-kill). The surviving window route. Does NOT stack with SVD (d@1e-4 80→24) |
| **IT6** segmented-bucket P² | stranger #7 | **KILL** (mandatory audit) — P² law real (exp 1.68–2.02, P=4 buys 12–17×) BUT segmenting breaks the mean wall: per-segment mean/cov d′ 0.8→1.8 (z2=5) / 11.6 (z2=1). Single-bucket stands |
| **IT7** contact coincidence veto | stranger #6 (D3/D5 FA repair) | **KILL** by power cost — contact arm exactly medium-blind (AUC 0.5009), medium FA repaired 17.8× (0.9965→0.056), but only 24% scene power retained (beyond-band change too weak in the mean channel). Fall back to semiparametric repair |
| **IT8** dither-render leak | finer #5 | **PWM CONSTRAINT NEEDED** — Bayer ordered-dither doubles the z1=0 leak (7e-3→1.38e-2 ≥1e-2); PWM temporal dithering + whole-period integration is a hard apparatus-spec line |
| **IT9** carrier-cross grid-independence | finer #2 residual | **CONFIRMED** — D=Re[C·conj(S)] exact to 5e-16 at N=2048 AND N=4096 (\|S\|² cancels, no 2kp edge); beyond-2kp leak 9.58e-3 identical (0.0% diff) → report §1 correction freezable |

## Bearing on the three 125-score R44 theorems
- **#1 Finite-Aperture Wall Impossibility / Pupil-Hardened Wall** — CONFIRMED (KT1):
  no exact open wall for a finite DMD; a hard pupil restores an exact field wall;
  exact detector wall needs a concentration-grade guard.
- **#2 Jet Transmutation** — CONFIRMED (KT2): mean-leak projection restores the m=2 class.
- **#3 Coherence-Separation Rank Bound** — bound + jet-order separation hold (KT4), but
  the twin geometry cannot reach R→1 / make R² binding → the memory-effect rank line
  still needs corrected bookkeeping + a genuine scrambling endpoint before it can anchor
  the next paper.
- New: a **SECOND exact wall** of energy-conservation type (IT1) and a viable
  **certified-blind code** route (IT4/IT4b), strongest at a fixed relay-conjugate plane.
