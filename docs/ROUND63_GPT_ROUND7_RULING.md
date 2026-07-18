# GPT round-7 FINAL F1 ruling (2026-07-18) — DETAIL-24 primary cohort + Q90 target

Reconstructed from screenshots (DLP-blocked tool extraction); on discrepancy
the ChatGPT thread governs. Operative for F1 together with rounds 4–6.

## Decision (verbatim opening)
Adopt P1, replace the defective Q_safe,max − 1 dB target with a normalized
safe-curve target, and adopt P2 only as a preregistered KEY SECONDARY
endpoint — not a co-primary veto. Two explicitly separate scene regimes:
- **DETAIL-24**: procedurally generated diagnostic targets; the confirmatory
  cohort for the positive acquisition-speed claim.
- **NATURAL-24**: the frozen STL-test naturals; a secondary regime-boundary
  cohort with NO positive pass gate.
"High-flux reconstruction is beneficial when the inverse problem is
photon-limited" — the scientifically correct narrowing.

## 1. DETAIL-24
**1.1 Composition** — exactly 24 images = six fixed families × four
independent instances:
1. Bitmap glyph/stroke targets (random strings from a repository-embedded
   5×7 alphanumeric bitmap alphabet; no system fonts)
2. Line-pair and spatial-frequency chirps (h/v/diagonal line pairs,
   widths/gaps 1–4 px)
3. Radial-spoke targets (Siemens-star-like, frozen spoke counts + center
   offsets)
4. Orthogonal maze/barcode targets (1–2 px connected stroke networks,
   frozen fill-factor ranges)
5. Thin-contour composites (ellipses/polygons/curved+straight outlines on a
   NONZERO background)
6. Band-pass microtexture targets (seeded Rademacher fields through fixed
   band-pass kernels, incl. moderate band)
**1.2** Implementation testing may use exactly SIX development instances
(one per family), seeds 630900..630905 — only to verify deterministic
generation, numerical stability, endpoint implementation, runtime/manifest
accounting. They must NOT change: c=0.50, generator parameter ranges,
endpoint constants, pass thresholds, family composition. The 24 confirmatory
instances must not be reconstructed before F1.
**1.3** Inference is conditional on the declared procedural target
distribution; the paper must NOT present DETAIL-24 as representative of
arbitrary natural imagery. Bootstrap stratified by family (each replicate
holds 4 sampled images from each of the six families).

c = 0.50 remains PERMANENTLY frozen (Pass A chose C0*=inf under the frozen
tie rule); it must NOT be recalibrated on the new cohort (no Pass-A rerun).

## 2. Endpoint
**2.1 Retire** Q*_j = Q_safe(ν_max) − 1 dB (falls below the observable knee
whenever the safe curve has < 1 dB dynamic range — the Pass-B defect).
**2.2 Monotone curves**: average PSNR_rad over the FIVE measurement seeds at
each ν; fit nondecreasing Q̃_{ρ,j}(log T_opt) by ordinary equal-weight PAVA
isotonic regression; preserve/report the raw seed-mean curve separately
(isotonic used only for endpoint extraction). Grid stays
ν ∈ {5,10,20,50,100,200,500,1000,2000}, M = n = 4096.
**2.3 New per-image target**: safe-side observable dynamic range
R_j = Q̃_{0.05,j}(T_max) − Q̃_{0.05,j}(T_min).
If R_j < 0.50 dB: status = SAFE_RANGE_UNINFORMATIVE, S_gate = 1.0, not
counted positive. Otherwise the target is (boxed)
  **Q90,j = Q̃_{0.05,j}(T_min) + 0.90 · R_j**
("time to recover 90% of the quality improvement the safe operating point
achieves over the tested budget"; the 10% margin avoids anchoring at the
noisiest terminal sample). Crossing times by linear interpolation in
log T_opt.
**2.4 Conservative censoring** — round-6 taxonomy augmented by
SAFE_RANGE_UNINFORMATIVE:
- both interior: S_j = T_safe,j / T_fast,j
- fast at floor, safe later: S_gate,j = T_safe,j / T_min (report S_j ≥ S_gate)
- both at floor: S_gate = 1, unresolved and not positive
- safe floor / fast interior: S_gate = T_min / T_fast,j < 1
- fast never reaches: S_gate = 0
- missing/incomplete safe curve: S_gate = 0 + analysis-completeness failure
No image is deleted because of censoring.
**2.5 Primary gate — constants UNCHANGED**: median_{j=1..24} S_gate,j ≥ 3
AND LB_{2.5%}[median S_gate] > 1 AND #{j: S_gate,j > 1} ≥ 18/24 ⇒
DETAIL_SPEED_PASS. Bootstrap = **10,000 nested family-stratified
replicates**: (1) within each image resample the five seeds with
replacement and REBUILD both curves, the target and the censoring; (2)
within each family resample its 4 images with replacement; (3) compute the
24-image median. Rationale kept verbatim: ideal information-rate ratio
between ρ=0.6 and 0.05 ≈ 7.9×, so a median ≥3× stays a meaningful,
conservative engineering claim.

## 3. Fixed-budget quality gain — key secondary, not co-primary
ΔQ^budget_j = Q_fast,j(ν=2000) − Q_safe,j(ν=2000) using the RAW five-seed
mean PSNR_rad at the terminal grid point (not isotonic). Call it
**fixed-budget quality gain**; do NOT call it asymptotic "ceiling lift" or
claim the safe arm can never reach that quality (with longer dwell it may
catch up; the claim is specifically about the frozen acquisition-time
budget). Preregistered descriptive success rule: median ΔQ ≥ 1.0 dB AND
LB_{2.5%}[median ΔQ] > 0 AND #{j: ΔQ_j > 0} ≥ 18/24, same 10k bootstrap.
Key confirmatory secondary; does not veto or rescue DETAIL_SPEED_PASS. If
speed fails and quality-gain passes: claim fixed-budget quality improvement
only, not speed.

## 4. NATURAL-24 (STL test) — secondary regime-boundary cohort
No positive gate. Report its curves, censoring classes, representative
reconstructions. Manuscript interpretation (verbatim): "For smooth or
prior-limited natural scenes, the safe operating point may reach its
recoverable spatial-information ceiling with very few detected photons,
leaving little photon-rate acceleration to measure. High flux is
advantageous primarily in the detail-resolving, photon-limited regime."
The natural cohort must be shown prominently enough that DETAIL-24 cannot
be mistaken for a universal-image claim.

## 5. Exact-likelihood patch at ν ≤ 2 — DO NOT add to production
Confirmatory grid begins at ν=5; switching estimators at short windows
would complicate the operating-point comparison; ν ≤ 2 is an identified RQL
boundary, not part of the deployment claim; at ν=1 the ceiling-count
fraction is a physical identifiability warning, not merely an optimizer
defect (the RQL objective loses its linear identifying term when Nτ ≥ T).
Include ν = 1, 2 only as an optional S4/Supplement reference: exact
Bernoulli/trinomial likelihood, RQL comparison, ceiling-count fraction, and
the statement that this is outside the production deployment zone. No ν < 5
result enters either confirmatory endpoint.

## 6. Pooled-flux initializer — CONFIRMED whitelist-legal
(Changes x0 for all arms identically.) Wording: do NOT describe it as an
exact finite-window estimator of the active-start mean; describe it as a
stable pooled moment initialization.

## 7. F1 lock conditions (before F1)
1. Freeze DETAIL-24 generator code, parameter manifest, seeds and image
   SHA256. 2. Freeze c=0.50; do not rerun Pass A. 3. Replace the old −1 dB
endpoint with Q90. 4. Implement PAVA isotonic fitting and the R_j < 0.50 dB
uninformative rule. 5. Implement the exact censoring taxonomy and the
10,000-draw nested family-stratified bootstrap. 6. Add the fixed-budget
secondary endpoint. 7. Regenerate manifests and expected-cell tables with
M = 4096, ν = {5,10,20,50,100,200,500,1000,2000}. 8. Run only an
implementation smoke on the six nonconfirmatory generator seeds — no
threshold or generator modification may follow from its PSNR results.
9. Freeze the commit and launch S2.
**After these changes, no further endpoint, image-class, regularization,
geometry or grid revision is permitted. If DETAIL-24 fails, that is the
campaign's confirmatory result.**
