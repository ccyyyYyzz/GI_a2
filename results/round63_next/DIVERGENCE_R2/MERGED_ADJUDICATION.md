# Merged adjudication — internal 5-lens divergence (wf_d590903e) × GPT R44 (issue #35)

Date 2026-07-24. Coordinator ruling merging the two independent divergence engines
run over the wave-twin T1–T5 findings. Internal: 30 ideas → 9 survivors after
adversarial verify (66 agents; full record `synthesis_full.json`). R44: 8 ranked
innovations with derivations. Both quarantined from the PRL Letter (only already-
executed R44 §3 transfer paragraph touched the manuscript).

## Cross-validation between the two engines (independent, convergent)

- **Same mechanism, found twice**: internal measured the exact complementary
  cancellation `D = 4·Re(C̄S̃)` (5e-16) = R44's derived `D_z = Re[U_z S_z*]`
  (eq 1.1). The leak is a carrier–code heterodyne cross term; **NO spectral edge
  at 2k_p**. The twin report §1's "2k_p intensity-autocorrelation tail" wording
  is corrected by both engines.
- **Divergence on witness relocation**: R44 recommends >2k_p placement (removes
  the *ideal* core). Internal **measured** it: buys ≤3.5× at z1=0, ~1× at
  z1≥5 mm — operationally near-useless because the carrier pedestal dominates.
  RULING: the Letter's frozen sentence stands (it claims only ideal-core
  removal + states the pupil requirement); bench protocol does NOT rely on
  witness relocation.
- **Pupils**: hard FIELD pupil (relay) = GOOD, measured leak z1-independent to
  3 digits (internal) and required for exactness (R44 Thm 1 case 2).
  Code-band apodization at the pupil = MEASURED KILL (rim effect, 4–20× worse).
  Naive separable windows (rect/Tukey/Gauss) = measured kills. Only
  DPSS/Slepian-grade concentration remains live (internal test 5).

## Three independent wall-restoration routes (the emerging next-paper spine)

1. **Optical**: hard Fourier relay pupil → exact difference-set wall
   (R44 Thm 1; Colab KT1 variant 6/7 decides).
2. **Code-space**: SVD of the measured linear leak operator L(z1) → certified
   mean-blind code subspace, bench-measurable as a calibration step
   (internal #1 sledgehammer; internal test 4 gates it).
3. **Statistical**: leak tangent as profiled nuisance; efficient score
   projection (R44 #6; Colab KT3 decides). Jet-transmutation restoration
   (R44 #2, KT2) is the theory statement of this route.

Plus a possible **second exact wall** of a different conservation type:
phase-only object change invisible to a lossless full-aperture bucket at every
order (internal higher #1; internal test 1 arm ii).

## Estimator-layer caveat (the internal round's missing-lens finding — binding)

Most internal kills were **inferential, not physical**: blind λ_cov estimation
sits 600–9000× under the Wishart M²/T noise floor (killed range-gating and
channel-switching as bench observables); frozen plug-in W anti-orders the
covariance channel at small λ; CRN/oracle ΔC quantities are twin-only, not
bench-observable. CONSEQUENCE for the Colab battery: **KT6 (channel-ratio
ranging) must report the blind-estimable version alongside the oracle version;
an oracle-only pass is NOT a pass.** KT4/KT5 unaffected (deterministic /
oracle-legitimate quantities).

## Master kill-test queue (Colab battery, wave 1 + wave 2)

Wave 1 (already dispatched): R44 KT1–KT6 with the KT6 caveat above.
Wave 2 (internal queue, ~200 GPU-min total, decision trees in
`synthesis_full.json → immediate_tests`):
IT1 two-wall + jet-order invariance (35m) · IT2 floor decomposition/leak-law
fit (10m) · IT3 leak-law seed/M stability (8m) · IT4 **SVD certified-blind
codes — the sledgehammer gate** (50m) · IT5 DPSS concentration (20m) ·
IT6 segmented-bucket P² with per-segment leak audit (15m) · IT7 contact-arm
coincidence veto — the D3/D5 FA repair candidate (40m) · IT8 dither-render
leak (12m) · IT9 carrier-cross grid-independence confirm (10m).

## Next-paper center of gravity — DEFERRED to test results

R44 ranks the memory-effect Rank–Jet phase diagram first (coherence-rank
theorem rank ≤ R²; first evidence = T5b's unsaturated M^1.8 scaling).
Internal ranks null-space wall engineering first and DEMOTES memory-effect
(registration inconsistency: r_ME vs r_ME² exponent clash; rank capped at 24
by the scene basis — reopen only with corrected bookkeeping). RULING: run the
tests; the center is chosen in an R45 round with KT/IT results on the table.
One-good-paper discipline applies: ONE center, the other becomes a section.

## Corrections applied to in-repo documents

- WAVE_TWIN_REPORT.md §1 mechanism sentence: 2k_p-tail attribution → carrier-
  tail (this commit). Letter untouched.
