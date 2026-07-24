# R45 — Deep mine first, then rule the architecture: the full kill-test record is in

The R44 battery + the internal divergence queue have run (Colab, float64, all
committed under `results/round63_next/R44_KILL_TESTS/` with per-test JSON + notes
+ `INDEX.md`). Paper 1 (the Letter) is untouched and submit-ready. This is a
COMBINED round by operator instruction: an R42-class maximal-deliberation mining
pass over the new evidence FIRST, and only then the paper-2 architecture ruling —
so the dig can reshape the frame instead of being fenced by it. Do not rush the
frame; where the dig opens something that needs a decisive test before the
architecture can be fixed, SAY SO and name the test — provisional rulings with a
named gate are allowed.

## Part A — the frozen verdict table (quote-grade)

| Test | Object | Verdict | Key numbers |
|---|---|---|---|
| KT1 | Finite-aperture wall impossibility / pupil-hardened wall (#1, 125) | **PASS (layered)** | ideal-code wall 2.1e-15; hard-pupil FIELD wall exact 7e-16 ∀z1,kR; full-twin pedestal 5.3e-2; finite object window reintroduces z1-independent detector tail |
| KT2 | Jet transmutation (#2, 125) | **PASS** | projected slope 4.00 ×4 cells; unprojected 2.12–2.20; first-order-null residual 1e-18 |
| KT3 | Leak-orthogonalized sentinel (#6, 100) | **PASS** | nuisance overlaps →0 exactly ×4; target info retained 86–99%; oracle-W small-λ anti-ordering flagged |
| KT4 | Coherence-rank bound (#3, 125) | **PARTIAL** | bound never violated ∀z2; orthogonal jet slope 4.0 ∀z2; BUT rank capped by M(M+1)/2=136 not R²≈12000 (non-binding), and R stays ~107–122 — a single thin screen NEVER collapses toward R=1 on the z2 axis |
| KT5 | Étendue/N_eff law (#4, 100) | **PARTIAL** | N_eff 46→811 with aperture; saturated-regime collapse CV 0.25; universal single-variable collapse fails (normCV 0.64) — two-regime law stands |
| KT6 | Channel-ratio ranging (#5, 80) | **ORACLE-PASS / BLIND-ESTIMATOR_KILLED** | z2⁴ exponents 3.93–4.08 ×4; blind λ_cov ~185–220× under the Wishart M²/T floor |
| IT1 | Second exact wall + jet invariance | **PASS** | energy-conservation wall: phase-only change invisible to lossless full-aperture bucket, null 2.18e-16 ∀z2,ε; NA-clipped slope 0.97–1.03 (ε¹); clip response climbs 6 orders over z2 |
| IT4 | SVD certified-blind codes (sledgehammer gate) | **PARTIAL** (original strict gate) | bottom-64 joint max 4.99e-4; λ ratio 1.20 (83% power kept); 219× operational leak cut |
| IT4b | Dimension/leak trade re-scoring | **PASS (operational gate)** | joint d@1e-4=35; d=64 RMS-of-drawn-codes 1.68e-4; **single-z1 relay-conjugate: 80 dims ALL ≤1e-5**; Tukey window KILL (93× worse, rim); ridge-stable |
| IT7 | Contact-arm coincidence veto (D3/D5 repair) | **KILL (power)** | hardware identity holds (contact medium-AUC 0.5009); medium FA repaired 0.9965→0.056 (17.8×); but coincidence keeps only 24% scene power — the beyond-band change barely leaks into the contact arm BY the wall itself; fallback = quarantined semiparametric repair; veto remains valid for in-band/amplitude classes (bench manual only) |
| IT5 | DPSS concentration window | **PASS, non-composable** | witness leak ≤2.38e-4 uniform in z1 (349× over bare, 168× over killed Tukey); BUT DPSS reshapes the leak operator — SVD-certified subspace shrinks 80→24 dims under it: window and code routes are ALTERNATIVES, not multipliers (partial answer to the certificate-composability vein) |
| IT6 | Segmented-bucket P² law | **KILL (mandatory audit)** | P² gain real (exponents 1.68–2.02; P=4 buys 12–17×) but segmentation breaks whole-bucket complementary cancellation — per-segment mean/cov d′ 1.4–11.9 floods attribution; consistent with IT1 (segmenting = per-segment NA clipping → ε¹ wall breaking); bench stays single-bucket |
| IT2/IT3 | Leak-law metrology | **QUOTABLE** | floor decomposes (fill 6% + macro-staircase 50%): L_pix=7.15e-3 is a real DMD instrument constant; law leak(z1)=6.94e-3+5.92e-3·z1^0.90, stable across seeds×M (CV ≤2.2%) — the CALIBRATED-LAW transfer sentence is unlocked |
| IT8 | Dither-render leak | **PWM REQUIRED** | ordered-dither render leaks 1.38e-2 beyond-band at z1=0 (≥1e-2 bar): PWM temporal dithering + whole-period bucket integration = hard apparatus-spec line |
| IT9 | Carrier-cross grid independence | **CONFIRMED** | D=Re[C·S̄] machine-zero on 2048² AND 4096²; pedestal leak grid-difference 0.0% — the twin-report §1 mechanism correction is frozen with confidence |

Plus: wave twin complete (T1–T5b; TWIN_CONFIRMS at matched cell 768∈[513,2996];
λ∝M^1.8; COMSOL full-Maxwell 15.5% agreement) and the internal divergence
synthesis (30→9 survivors; `results/round63_next/DIVERGENCE_R2/`).

## Part B — THE DIG (do this FIRST; R42 protocol; maximal deliberation)

Mandatory ore veins — examine each, derive before proposing, then go beyond:

1. **A Noether theorem for blind channels?** Two exact walls of different
   conservation types now exist: support-type at the modulator (band-limited code
   algebra; physically restored by a hard pupil, field-exact 7e-16) and
   energy-type at the collector (phase-only change invisible to a lossless
   full-aperture bucket, 2.18e-16; broken at ε¹ by NA clipping). Each has its own
   geometric softening knob (defocus z1 vs collection NA). Is EVERY
   symmetry/conservation of an optical measurement chain the generator of an
   exact statistical blindness? Enumerate the walls of a scalar-diffraction chain
   (energy, support, parity, reciprocity, time-reversal, polarization…), each
   with its exact null, breaking operator, and restoration hardware. Derive at
   least two NEW walls or prove the list complete.
2. **The geometric jet flow.** IT1's clipped response climbs SIX ORDERS
   monotonically over z2 at fixed slope 1; KT2's m_eff(ε) crosses 1→2 at
   ε_cross; your R44 §2.4 gave bi-orders (m_ε, m_z). Is there a genuine
   two-variable jet FLOW m_eff(ε, geometry) with universal crossover exponents —
   RG-like? Derive the flow equation or kill the analogy.
3. **The estimator-layer impossibility (a possible THIRD ledger).** Blind λ_cov
   sits 185–9000× under the Wishart M²/T floor; plug-in matched filters
   anti-order at small λ; oracle quantities are twin-only. Derive the exact
   second-kind impossibility: WHICH functionals of the record a bench can never
   blindly certify at finite (M,T) — a minimax statement whose equality cases
   define the certifiable set. Physical wall / channel demand / estimator floor
   would then be three distinct ledgers.
4. **Blindness certificates as a metrological primitive.** The leak operator is
   exactly linear in the code; one calibration run measures it; its null space
   gives 80-dim ≤1e-5 certified codes. What is the general theory of measurable
   CANNOT-see guarantees? Compare honestly: decoherence-free subspaces, NMR/
   dynamical decoupling, null-space control, differential privacy. What exact new
   statement does optics get? Is the certificate COMPOSABLE (pupil × code ×
   projection = product of leak factors — derive or refute multiplicativity)?
5. **The unsaturation crossover.** λ∝M^1.8 at finite z2 vs exact M_eff≈R
   saturation at complete scrambling. Derive M_eff(M;R) interpolating the two
   regimes (plausibly a crossover at M*~R) — this would make the coherence-rank
   object measurable as a scaling exponent WITHOUT the R→1 endpoint geometry
   that KT4 showed a single thin screen cannot reach.
6. **The étendue sweet spot.** From I=½Σ[r/(1+r)]², r=n/N_eff: is there a
   maximum-information étendue at fixed photon budget (an optimal aperture the
   bench should sit at)? Derive the two-variable scaling function; KT5 measured
   the saturated-regime collapse (CV 0.25) and the global failure (CV 0.64).

Then the R42 moves: assumption-inversion sweep (≥8 — e.g. "walls protect the
scene from the medium", "certification is per-code", "the bucket is the
detector", "banks are the currency", "blindness is a defense rather than a
product"); cross-field transplants (decoherence-free subspaces, Noether/gauge,
RG, Le Cam minimax, NMR decoupling, control observability); named-effect audit;
**≥7 ranked ideas** (surprise × rigor × reach, each with mechanism, exact-statement
sketch, decisive ≤60-GPU-min kill test, nearest prior art); and THE MOONSHOT:
derive your deepest candidate in full — assumptions, identity/inequality,
equality cases, falsification numbers. R42's moonshot survived 17/17; that is
the bar.

## Part C — architecture ruling (AFTER the dig, informed by it)

1. **Center of gravity for paper 2** — now a THREE-way: (a) certified-wall
   engineering (two exact walls, three restoration routes with certified
   numbers, jet transmutation as the theory); (b) the memory-effect Rank–Jet
   phase diagram REBUILT on a corrected disorder axis (screen stacking /
   thickness — KT4 measured that single-thin-screen z2 does NOT scramble); or
   (c) something the dig just surfaced (e.g. the wall-classification theorem or
   the third-ledger impossibility as the new protagonist). One-good-paper
   discipline binds: ONE center; everything else becomes sections. Provisional
   ruling with a named decisive gate is allowed if the dig demands a test first.
2. **The theorem set, hero figure, title candidates, venue** (Optica vs PR
   Applied vs PRX — argue against yourself honestly), and any campaign beyond
   the battery.
3. **D3/D5 disposition** given IT7's KILL (hardware repair dead for the
   beyond-band class by the wall's own physics; 17.8× medium-FA repair real but
   76% power cost): rule whether the specificity dent stays as-disclosed with the
   semiparametric repair quarantined-future, and whether the IT7 diagnostics
   merit a methods note.
4. **Transfer-paragraph conditionals — now unlocked.** IT2/IT3 passed
   (instrument-constant decomposition + CV ≤2.2% stability), so the internal
   synthesis's CALIBRATED-LAW sentence is eligible: "the mean channel is blind up
   to a calibrated residual obeying leak(z1)=L_pix+a·z1^p with L_pix=6.9e-3 set
   by DMD pixelation (fill + macro-staircase), p≈0.90, suppressed by
   near-conjugate relay imaging" — rule the final wording (this touches ONLY the
   already-written End Matter C transfer paragraph's leak sentence; everything
   else in the Letter stays frozen). Also rule whether IT8's PWM line belongs in
   the Letter's transfer paragraph or the bench protocol only.

Constraints unchanged: paper 1 frozen; sim-only until publication; dead lines
dead; honesty machinery symmetric (your ranging application died at its own
floor calculation this week — the system works). Deliver as ONE GitHub issue
titled R45. Take maximal deliberation on Part B before touching Part C.
