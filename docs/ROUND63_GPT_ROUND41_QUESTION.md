# R41 — Admitting the demand tier: beyond-band anomaly detection into the aperture atlas

Housekeeping since R40. Three developments, all committed or committing:

1. **Shot-model bug found and fixed (multi-engine discipline paid off).** The
   living_region_map/p1_fisher engine used R_shot ∝ clamp(Px,1e-12) — near
   zero for zero-mean signed codes — inflating Fisher ~1.9x. Physical shot =
   |P|·x throughput (medium/shot ≈ 12–15). Direction of correction: ALL kill
   verdicts HARDEN (f_rec drops), all positive headlines shrink. P1-FAIL and
   REGION_EMPTY stand, corrected artifacts + CORRECTION_NOTE archived. The
   manuscript quotes corrected numbers only.
2. **Pocket reconciled and demoted further.** The 8x map-vs-demo CRB gap
   decomposes exactly (T_eff 4096 vs 200: 4.53x; shot bug 1.88x; scene set
   1.16x). Honest pocket headline: blind natural-median 0.30 needs ~2028
   banks ≈ 26 s (witness: 112 banks ≈ 1.4 s); blind ≈ CRB (0.98–0.99x — the
   lifted convex estimator was built, validated convex, and did NOT beat the
   already-CRB-efficient moment pipeline; the E8 plateau was the finite-T
   CRB, not an estimator gap; Theorem 3 confirmed AND practically attained).
   Bank-acquisition model clarified: quasi-static diffuser per 2M-exposure
   bank, T_eff = #banks, 12.8 ms/bank at 20 kHz, no tau penalty.
3. **The demand axis was attacked (the operator's divergence mandate) and it
   is GREEN.** Detection needs only aggregate Fisher; estimation died on
   mode coverage. Corrected 81-cell detection map + Monte-Carlo ROC:

   - DETECTION_REGION_FOUND (corrected): eps=2% beyond-band anomaly
     detectable (d'>=5) in 95% of cells, best 467 banks = 6.0 s; eps=5%:
     100% of cells (worst-mode 91%), 1.9 s; eps=1%: best cells only,
     ~1869 banks = 23.9 s; eps=0.5%: below floor everywhere. Best cell
     moved to sigma_f=1.0 — with physical shot, detection is
     throughput-limited so DEEPER fluctuation HELPS (opposite to
     estimation's sigma_f-inertness): the worse the medium, the better
     the detector.
   - ROC validation: real profiled matched-score covariance-LR detector,
     250+250 records/point, empirical d'/analytic = 0.98–1.04 across all
     9 (cell x eps) points; AUC 0.93–1.00. Prognosis engine
     Monte-Carlo-sound.
   - Specificity operational (four-way): scene-inband events → t_mean
     d'=61 (beyond-detector −0.1); scene-beyond → t_cov 13.5 (others
     0.1); medium-law change (sigma_f+20%) → t_amp 61 (beyond 0.3). The
     beyond-band detector fires ONLY for beyond-band scene change and
     rejects both in-band scene changes and medium changes — the exact
     mean-wall theorem operating as a zero-false-alarm discriminator.
     Three orthogonal statistic axes (t_mean, t_cov, t_amp) = a
     scene-inband / scene-beyond / medium classifier seed; the DLGI
     reciprocity machinery is the natural attribution theory. A tau-type
     medium change (lag-decay signature) was NOT exercised by the
     iid-bank MC — flagged for the sealed probe.
   - Robustness: 10% declared-basis rotation → AUC 1.000→1.000, d'
     5.30→5.31. Mechanism: the aggregate statistic reads tr(J_B), which
     is invariant under within-aperture orthogonal reshuffles of the
     medium basis — a free lemma ("trace-invariance robustness") we ask
     you to formalize.
   - Deployment: running MxM covariance + one quadratic form, ~16k
     MAC/bank, real-time trivial; no reconstruction anywhere.

4. **Hostile prior-art scout (detection square): NOVEL_SQUARE_OPEN,
   crowded neighborhood.** Zero optics hits on the two pillars
   (mean-invariance wall as specificity discriminator; unmeasured
   fluctuating medium as transducer). Must buy off head-on: Grace–Guha–
   Dutton arXiv:2308.07262 (SPADE sub-diffraction change detection —
   different channel entirely), SAR Coherent Change Detection (covariance
   catches sub-resolution change — the classical math family), coda-wave
   interferometry (cross-domain spirit twin), second-order speckle SR +
   blind-SIM (the covariance-band family, already fenced in R39/R40).
   Scout's frozen claim sentence is on file.

Also: R40 figure layer is built (Fig 1 master atlas passed the coordinator
gate, three-pass self-critique recorded; 10 frozen source tables, RULE ZERO;
pocket/branch variants prepared; committed @ 063ef05).

## R41 asks

1. **Admission ruling.** Does the detection capability enter the flagship
   as (a) the demand-tier row of Fig 1 + a constructive Section (the atlas
   predicts the one affordable demand and the paper delivers it), (b) a
   compressed paragraph + supplement, or (c) a separate later paper
   (violating the one-paper rule — argue only if you believe it)? If (a):
   re-cut the R40 page map (what shrinks to make ~1 page of room; we
   suspect Section 3 or 5) and the Fig 1 third-row design (demand tiers:
   estimation ✕ / detection ✓, with the same support/supply geometry).
2. **The frozen claim wording** for the detection capability, integrating
   the scout's sentence with the corrected numbers (>=2% energy, 6 s,
   95% of physical space; 1% best-cells 24 s; 0.5% out of reach; the
   sigma_f=1.0 deep-fluctuation advantage) and the four-way specificity.
   Also the forbidden-claims list for this capability.
3. **The trace-invariance lemma.** Formal statement + proof sketch of the
   aggregate-detection robustness (tr(J_B) invariance under within-aperture
   basis rotations of the declared law; when does it fail — non-orthogonal
   reshuffles, aperture-boundary changes?).
4. **Sealed detection probe spec** (the paper's constructive demo must be
   preregistered-grade): banks, arms (incl. fresh-pattern mean arm as the
   zero-information control, best-effort fresh-covariance arm per your P4
   logic — does F4-immunity hold for detection? the wall says the MEAN arm
   is blind; the fresh-COVARIANCE arm is not obviously blind for
   detection — rule the honest comparator), bars (d' vs analytic within
   20%; specificity confusion matrix; tau-change class; mismatch axes),
   kill tree, compute budget (local GPU scale).
5. **The pocket sixth sentence, final.** With the corrected numbers
   (blind ≈ CRB efficiency achieved; natural-median 0.30 at ~26 s;
   witness 1.4 s), which R40 branch wording applies — (A) achieved-as-
   calibration-point with corrected honest magnitudes, or a hybrid
   wording? Freeze the sentence.
6. **Disclosure wording** for the shot-model correction (found by the
   mandated multi-engine cross-check; hardens all kills; corrected
   artifacts only are quoted) — one referee-proof sentence.
7. **Title/thesis impact.** Does the green demand tier change the frozen
   thesis/title (the demand-relative extension: "usability is
   demand-relative — the same aperture that cannot support imaging
   supports detection")? Rule whether "The information aperture of a
   bucket record" stands or a variant is stronger, given the paper now
   contains a certified-positive mean layer, a dead imaging claim, AND a
   validated detection capability.
8. **Venue impact.** With a working, validated, novel-square detection
   capability inside: does the OE branch reopen / Optica become arguable,
   or does JOSA A/TCI remain the honest first target?

Deliver as a GitHub issue titled R41. All artifacts committed at 063ef05
(+ corrected detection/P1/map artifacts committing next).
