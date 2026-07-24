# READ2 — HOSTILE REFEREE REPORT

**Manuscript:** *Engineering What an Optical Sensor Cannot See* (Optica Research Article; PR Applied fallback)
**Repo/commit:** `github.com/ccyyyYyzz/GI_a2`, branch `main`, HEAD `0f9f48d`
**Reviewer stance:** maximally hostile; goal is to reject. Worked independently of any READ1 report.
**Materials read:** `paper2_optica/main.tex` + `main.pdf`, `supplement.tex` + `supplement.pdf`, `CLAIM_SOURCE_MATRIX.md`; cross-examined against `paper_prl/main.tex`/`supplement` (frozen companion Letter) and the committed JSON under `results/round63_next/`.

**Provenance spot-check (traced to the committed artifact, not just the matrix row):** KT1 `6.591e-16` (field wall), KT1 `P1=2.070e-15`, KT1 `P3=0.053154`/`0.16356`; IT1 `2.1798e-16`, IT1 clip slopes `[1.033,0.982,0.970,1.000]`; KT2 projected `[4,4,4,4]` / unproj `[2.136,2.20,2.117,2.176]` / generic `[1.576,2.284,1.576,2.285]`; IT4b `d@1e-5=80`, `joint d@1e-4=35`, `σ1=1.425e-6`/`σ64=4.99e-4`/`σ120=0.1125`, apod `d@1e-4=24`; IT4 `lam_ratio=1.202`, `219×`; IT5 `dpss=2.378e-4`/`none=8.30e-2`/`349.29`, `80→24`; KT6 `[4.081,4.054,4.022,3.989]`, `s/floor≈0.0054`, `ESTIMATOR_KILLED`; IT3 `6.94e-3/5.92e-3/0.904`, CVs `2.2/2.2/0.1%`; IT8 `1.384e-2/7.01e-3`; TLSG bars 1–7; SCRAMBLE `0.99992`; JET_TEST `2.038/4.000`, anchor `0.20`. **Every number reproduced exactly.** The arithmetic and provenance discipline are genuinely clean; the manuscript's weaknesses are therefore *not* fabrication — they are novelty, theorem scope, overclaim, selective reporting, and prior-art collapse. That is where I attack.

**Severity tally:** 2 REJECT-LEVEL, 7 MAJOR, 6 MINOR.

---

## REJECT-LEVEL FINDINGS

### R1 — The three "theorems" are classical results renamed; the paper's own proofs admit it
**Location:** Thm 1 (Statistical Noether wall, main §1 / S1); Thm 2 (Blindness capacity, §2 / S2); Thm 3 (Three-ledger, §4 / S4); jet corollary (§3 / S3).
**Attack:** The load-bearing contribution is "a hierarchy of three theorems." Under scrutiny each is a textbook result in optics vocabulary, and the manuscript itself supplies the citation that kills it.

- **Thm 2 `β_d(L)=σ_d(L)` is Courant–Fischer.** The proof sketch (§2) states verbatim that Eq. (2) "is the Courant–Fischer min–max characterization of σ_d," the robust envelope Eq. (3) "is Weyl's inequality," and the noncomposability Eq. (5) is the block-Gram identity `L_stack†L_stack = Σ w_i² L_i†L_i`. A "Blindness Capacity Theorem" that is Courant–Fischer + Weyl + a block matrix multiply is not a theorem; it is nomenclature. The physical content — "the least-observable code subspace is the bottom singular subspace of the leak operator" — is exactly **truncated-SVD / Backus–Gilbert null-space resolution analysis** (also the noise-subspace result underlying MUSIC), a half-century old in linear inverse theory. No delta over Backus–Gilbert is stated.
- **Thm 1 (Noether wall)** is the classical fact *group invariance of a dominated likelihood ⇒ score vanishes along orbit tangents ⇒ Fisher information singular there* (van der Vaart, invariance chapter — which the paper cites as [Vaart]). Worse, the theorem's headline clause — "no estimator, moment order, exposure count, or reconstruction prior can distinguish them" — is a **tautology**: it follows from `P_{g·x}=P_x` (identical laws are identical). The only non-trivial content is the infinitesimal score identity, which is the cited classical result.
- **Thm 3 (`1:√p:p`)** is the standard triad: matched filter needs SNR `O(1)`; the omnibus/χ² detection boundary against a p-dimensional alternative scales as `√p` (Ingster–Suslina high-dimensional detection boundary; the supplement's own "uniform spherical mixture / χ²-divergence `T²λ²/p`" argument is Ingster's second-moment method); estimation risk `p/T` is Cramér–Rao. None of Ingster, Donoho–Jin, or Barrett–Myers is cited.

**Evidence:** main.tex L251–257 (Courant–Fischer/Weyl admission), L166–178 & S1 L93–102 (invariance⇒score-null), L394–403 & S4 L331–349 (Ingster/CRB argument).
**Fix that would neutralize:** Demote all three to "Propositions (classical; restated for the optical setting)"; add explicit deltas over Backus–Gilbert/TSVD, van der Vaart invariance, and Ingster/Barrett–Myers; and relocate the paper's genuine novelty (a *measured hardware* dimension–leak curve with a calibration-transfer envelope) to the front as the sole contribution. As written, the novelty claim ("a hierarchy of three theorems… turns invisibility into a metrological primitive," L107–108) is unsupportable and, for a venue that gates on novelty, is a standalone rejection basis.

### R2 — "Full-Maxwell validation" cherry-picks the one favorable statistic; the same committed artifact records a 65% disagreement
**Location:** main §5 L465–466 ("A full-Maxwell check anchors the model: the thin-phase-screen speckle contrast 0.467 agrees with a … finite-element … contrast 0.395 to 15.5%"); supplement S6 L452–454; matrix F3 ("VALIDATED").
**Attack:** The manuscript reports **only** the ensemble-mean contrast agreement (15.5%) and calls the model "VALIDATED / anchored." The committed source `WAVE_TWIN/COMSOL_MICRO/COMSOL_developed_results.json` `summary` block contains, side by side with `contrast_rel_diff=0.155`, the value **`grain_rel_diff=0.6523`** — the speckle grain size (arguably the more diagnostic speckle statistic) disagrees by 65% (FEM mean 1.347 µm vs thin-screen 3.875 µm). The per-realization contrast disagreements are 11.6%–**60.9%** (seed 5: FEM 0.782 vs thin-screen 1.258), over **n = 8** realizations. A referee reads this as: the reduced model agrees with full Maxwell on one ensemble-averaged scalar to 15.5% while disagreeing on grain size by 65%, and the manuscript quotes the favorable number and suppresses the unfavorable one that sits in the same JSON dictionary.
**Evidence:** `COMSOL_developed_results.json` → `summary.contrast_rel_diff=0.155`, `summary.grain_rel_diff=0.6523`, `rows[0..7]`; recomputed per-seed |Δ|/FEM = [35.9, 16.5, 13.1, 22.0, 11.6, 60.9, 26.9, 15.9]%.
**Fix that would neutralize:** Report both statistics ("contrast agrees to 15.5%; grain size to only 65%, n=8"), downgrade "validation/anchors" to "order-of-magnitude consistency of ensemble contrast," and drop "VALIDATED." Until then this is selective reporting supporting a validation claim — an integrity flag, not a rounding quibble.

---

## MAJOR FINDINGS

### M1 — Theorem 3 is a Gaussian-covariance-model result, but the STATEMENT does not say so
**Location:** Thm 3 main L375–392; the Gaussian assumption appears only in supplement S4.1 L304 ("One independent bank yields `Y_t∼N(0,Σ_ε)`") and in the score `s_j=½[YᵀB_jY−tr B_j]`.
**Attack:** The main-text statement says "efficient covariance displacement," "score experiment," "locally `N(√Tθ,I_p)`" — never "Gaussian." The entire LAN reduction, the `½‖A_⊥‖_F²` signal, and the χ²/norm-test machinery are specific to a zero-mean Gaussian covariance family. A reader takes the `1:√p:p` ladder as model-free; it is not. Scope buried in methods is exactly what a hostile referee penalizes.
**Fix:** Add "for a Gaussian covariance model (LAN reduction)" to the theorem statement and to the abstract's ledger sentence.

### M2 — The calibrated leak law is published verbatim in the companion Letter; §5 reuses it without attribution
**Location:** paper2 §5 Eq. (7), main.tex L453; identical to paper_prl End Matter C, main.tex L423.
**Attack:** `L(z_1)=6.94×10⁻³+5.92×10⁻³(z_1/mm)^{0.90}`, the same three CVs (2.2/2.2/0.1%), and the same fill-factor/macro-pixel interpretation appear in **both** manuscripts. The introduction concedes the Letter "owns … the algebraic mean wall" (L100–102), yet §5 re-presents that wall's calibrated realization as this paper's own wave-optics transfer result, and **no `[COMPANION]` citation appears on Eq. (7) or anywhere in §5** (companion cites are only at L86, L98, L300). For Optica this is a redundant-publication / salami flag: the two papers share the same bench artifact (IT3/WAVE_TWIN) with overlapping text.
**Fix:** Cite `[COMPANION]` at Eq. (7), state explicitly which paper owns the leak law, and reduce §5 to the parts genuinely new to Paper 2 (COMSOL, `M^{1.8}` scaling, PWM/dither).

### M3 — TLSG exponents `−0.004/0.444/1.010` carry no CI, no SE, rest on 4 distinct p-values and a single seed
**Location:** main §4 L408–409, abstract L82; supplement Table S2 bar 4; TLSG bar4.
**Attack:** The ladder's headline numerical confirmation is three exponents quoted to 3 decimals with **zero uncertainty**. Fit is `log(contour) vs log(p)` over only **four distinct p** (10, 36, 136, 528; 8 cells). The blind exponent is **0.444 — 11% below the theoretical 0.5** — and is presented as confirmation of ½ with no error bar; only the frozen ±0.15 gate makes it "pass." Everything is a **single seed** (20260724), one sealed run, no bootstrap: the reported CVs (0.125, 0.124, 0.035) are cross-cell dispersions, not sampling uncertainty on the exponents. Reproduced here: blind exponent 0.444 (R²=0.983 on 8 pts), matched −0.004 (R²=0.003), estimation 1.010 (R²=0.9995).
**Fix:** Report bootstrap CIs on all three exponents over ≥2 seeds; state n and the p-grid in the caption; acknowledge the blind estimate sits low.

### M4 — The "measured projected slope 4.00" is analytically forced, and Fig 2 is a reconstruction from 2 fitted coefficients per cell
**Location:** §3 Route 1 L330–334, Fig 2 caption L291–300 ("reconstructed from the committed fitted coefficients"; "the four measured projected slopes are 4.00").
**Attack:** `KT2 cov_only_slope_fit = 4.0` **exactly** in all four cells, while the generic-direction fits in the same file are non-integer (1.576, 2.284). Exactly 4.000 for a "measured" log-log slope is the signature of a value that is *imposed*, not measured: once the ε² mean-leak coefficient A is projected out, the covariance channel is a pure `B ε⁴` monomial whose log-log slope is 4 by construction. Fig 2 (left) is explicitly "reconstructed from the committed fitted coefficients" — i.e., a two-parameter (A,B) analytic curve drawn with the visual authority of data. The claim "restores the covariance-only slope to 4.00" is therefore near-circular.
**Fix:** Plot the raw C(ε) sweep points under the reconstructed curve; report the slope as the *analytic consequence* of zeroing A, not as an independent measurement; give the fit residual.

### M5 — Missing prior-art pillars: physical-layer security, compressed-sensing certifiability, high-dimensional detection, ideal-observer detectability
**Location:** Discussion L490–494 cites only decoherence-free subspaces, dynamical decoupling, nonlinear-control unobservability.
**Attack:** "Engineering a subspace a sensor cannot see" with "a declared leak" and "the largest guaranteed-invisible code space" is the **wiretap/physical-layer-security** framing (β_d is a leakage/secrecy bound); none cited. "Certifiability ladder" invites **compressed-sensing certificates** (RIP, null-space property); none cited. The `√p` blind boundary is **Ingster**'s; the decision-theoretic "cannot-see" is **Barrett–Myers** ideal-observer detectability; neither cited. A knowledgeable referee will say each pillar is an existing framework and demand the explicit delta.
**Fix:** Add these citations and a one-sentence delta each; the honest delta ("hardware-calibrated, measured capacity curve") is defensible but must be stated against the right baselines.

### M6 — "Exactly invisible" / "exact wall" in abstract vs a 6.94×10⁻³ bench residual (13 orders larger)
**Location:** abstract L75–84 ("makes a scene subspace exactly invisible"; "an exact wall … restores"); Fig 1 L144 ("both machine zero"); vs §5 L452–453.
**Attack:** The abstract and intro lead with exactness and machine-zero walls (6.6×10⁻¹⁶, 2.18×10⁻¹⁶); the paper's own §5 shows that on the unfiltered bench the wall is `6.94×10⁻³` — thirteen orders of magnitude larger. The analytic theorem is exact; the *physical realization the paper offers* is not, and the reader meets the caveat only in §5. Promoting machine-precision numerics to "exact" is precisely the overclaim vector flagged for this review.
**Fix:** Qualify the abstract: "exact under the idealized band-limited model; on a bench the wall becomes a calibrated ~10⁻³ residual (§5)."

### M7 — Mode-unitary/parity/polarization walls are asserted as established but never measured
**Location:** §1 L205–207 ("blind not only to phase screens but to any unitary spatial-mode transformation"); Table S1 (six-row "wall periodic table," four rows unmeasured, labeled "immediate consequences").
**Attack:** Only the phase-only (diagonal) subgroup is measured (null 2.18×10⁻¹⁶). The full mode-unitary, parity, global-phase, and polarization-`SU(2)` rows are theoretical extrapolations presented in a "periodic table" that reads as a catalog of demonstrated walls. "Blind to any unitary spatial-mode transformation" is an unmeasured universal.
**Fix:** Mark the four unmeasured rows as conjectural corollaries, not confirmed walls; soften the universal.

---

## MINOR FINDINGS

### m1 — "Maximum field leak … 6.6×10⁻¹⁶ across all … z1" is contradicted by the matrix's own 7.70×10⁻¹⁶
**Location:** §1 L197–199; CLAIM_SOURCE_MATRIX GAP-1 ("largest single per-cell field leak … is 7.70×10⁻¹⁶"). The quoted 6.591e-16 is the max over per-variant summary maxima, not the true maximum. Both machine zero, but the word "maximum" is factually wrong by the authors' own admission. **Fix:** say "per-variant maximum" or quote 7.7×10⁻¹⁶.

### m2 — "Coverage 1.000" stated in §2 without its n=24 denominator
**Location:** §2 L278–282 ("coverage 1.000"); denominator (24 fresh draws) appears only in Fig 3 caption L442. Coverage 1.000 on 24 draws has a Wilson 95% lower bound ≈ 0.86; presenting "1.000" as the certifiable property overstates it. **Fix:** append "(24/24 fresh draws)" at first use.

### m3 — Fig 1 overlays single-plane and joint-four-plane capacity points on one "capacity curve"
**Location:** Fig 1 caption L146–151; GAP-2/GAP-3. The plotted σ_d curve is the joint spectrum, yet the "80 dims ≤10⁻⁵" (single-plane) and "35 dims ≤10⁻⁴" (joint) points are from different operators/conditions marked on the same axis. Condition conflation the authors themselves flag. **Fix:** separate axes or clearly label the operator for each point.

### m4 — Not submission-ready; reproducibility gated by a placeholder URL and commercial software
**Location:** L507–512, S6 L492–493. Authors, affiliation, funding, and Data-availability `[URL]` are all placeholders. Engine scripts *are* committed in-repo (good: `tlsg_A/B_*.py`, `IT*/KT*.py`, `run_comsol_micro.py`), but external reproduction depends on an unspecified URL to what appears to be a private repo, and the "full-Maxwell" check requires a COMSOL 6.3 license. **Fix:** public archival DOI (Zenodo) with seeds/hashes; note the COMSOL license dependency.

### m5 — Length/abstract exceed Optica Research Article format
**Location:** GAP-4; abstract ~150 words / target 10–12 pp vs Optica guideline 6–8 pp, ~100-word abstract. Desk-reject risk on format before science is read. **Fix:** cut to guideline or redirect to PR Applied (already the stated fallback).

### m6 — Theorem 1's "no estimator can distinguish" is a tautology, not a result
**Location:** Thm 1 L166–178. Restating that identical laws are indistinguishable adds no content beyond the score identity. **Fix:** state only the score-null/Fisher-singularity consequence as the mathematical content.

---

## VERDICT

**Would I reject at Optica? YES — reject (or transfer to PR Applied) in the present form.** The numbers are honest and fully traceable, and the measured hardware dimension–leak curve with a calibration-transfer envelope is a real, publishable object. But the paper is *sold* as "a hierarchy of three theorems," and on inspection all three are classical results renamed in optics vocabulary — the manuscript's own proof sketches identify them as Courant–Fischer, Weyl, invariance⇒singular-Fisher, and the Ingster/Cramér–Rao boundary. For a novelty-gated flagship that is fatal as framed. Compounding it: a "full-Maxwell validation" that quotes 15.5% and hides a 65% grain-size disagreement from the same file, a leak law shared verbatim with the companion Letter and reused uncited, and headline exponents with no error bars from four p-points and one seed.

**Single most dangerous finding: R1** — the three central "theorems" collapse to textbook linear algebra and classical detection theory, admitted in the paper's own proofs. It removes the stated basis for contribution and cannot be patched by a wording change; it forces a reframing of the entire paper around its one genuine novelty (the measured, calibration-transferred capacity curve), at which point the venue question (Optica flagship vs PR Applied) reopens. **R2** (the COMSOL selective-reporting) is the most dangerous *integrity* finding and is the one an editor is most likely to act on first.
