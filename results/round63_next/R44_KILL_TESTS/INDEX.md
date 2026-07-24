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

## Wave-2 (internal divergence queue) — priority subset done

| Test | Source | Verdict |
|---|---|---|
| **IT1** two-wall / jet invariance | higher #1 | **PASS — SECOND EXACT WALL** (energy-conservation type): phase-only object change invisible to a lossless full-aperture bucket at every order (null 2.18e-16 ∀z2); NA-clipped bucket responds at ε¹ (slopes 0.97–1.03), magnitude monotone in z2 |
| **IT4** SVD certified-blind codes | finer #3 (sledgehammer gate) | **PARTIAL/near-PASS** — mean-blind subspace exists, covariance retained (λ ratio 1.20 ≪ 3×), 219× lower operational mean leak; misses strict 64-dim ≤1e-4 gate by ~5× (a ~48-dim subspace clears it) |
| **IT4b** dimension/leak trade | coordinator follow-up | joint d@1e-4=35 (fundamental); d=64 RMS random-combo 1.68e-4 (3× below gated max, conservative); **single-z1 relay-conjugate d@1e-4=80** (+45); Tukey window is a KILL (rim effect); ridge-stable |

## Not run this session (deferred, WAVE-2 secondary/confirmatory)
IT7 (contact-arm coincidence veto — D3/D5 FA repair candidate; sketch: contact arm at
z2≈0 is intensity-blind to medium events by the T2 contact identity but sees scene
changes, so a 2D coincidence with the z2=5mm covariance arm vetoes medium false alarms —
worth building next), IT2 (floor decomposition), IT3 (leak-law seed/M stability), IT5
(DPSS concentration — the live window route after IT4b killed naive Tukey), IT6
(segmented-bucket P²), IT8 (dither-render leak), IT9 (carrier-cross grid-independence).

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
