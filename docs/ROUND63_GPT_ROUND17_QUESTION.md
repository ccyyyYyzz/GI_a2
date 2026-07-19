# ROUND63 — GPT Round 17 consultation (pre-freeze STOP: OED-DT gate arm is vacuous under the frozen guard stack)

Status: **m1-freeze is HELD.** The freeze-9 selftest reached 13/14 green, but the one
"green" that matters most (box 11/13, final exact-row full-guard suite) passed in the
degenerate sense: the accepted mixture has **m = 0 OED rows** for BOTH policies —
the deployed OED-DT arm would be 100% its own FIXED_DOSE fallback and the primary
gate would be empty. Per the standing stop rule, **no confirmatory scene has been
opened**; everything below is DEV/selftest evidence. Full durable evidence file:
`results/round63_m1/R17_CONSULTATION_EVIDENCE.md` (same commit as this question).

We ask R17 to adjudicate the architecture BEFORE launch. This is the last blocker;
all other boxes (1–9, 12–14) pass on the corrected pipeline, including
FIXED\* = SCAT32 by the frozen R10 DEV-PSNR rule (17.498 dB, no tie).

---

## 1. The frozen guard stack (R10 G1–G6 + R11 + R13/R15 amendments, as implemented)

The OED-DT arm accepts a 972-row design only if ALL of:

- **budget**: incident-energy sum ≤ dwell-matched cap (equality allowed, tol 1e-9);
- **dose (G5)**: per-pixel incident dose within ±5% of uniform (band ±0.05 relative);
- **ceiling**: per-row load ≤ ρ_R (clip NONE);
- **effD**: effective-dimension floor vs the balanced reference;
- **A-risk**: worst-case per-pixel amplification ≤ 1.05× the FIXED\* (SCAT32) reference;
- **spectral**: minimum trust spectrum ≥ 0.5× the FIXED\* reference.

Acceptance rule (`oed_design_v5.mixture_search`): descending scan over m
(OED rows retained, fallback rows filling the rest), accept the FIRST m where the
full suite passes; `alpha = m/972`.

## 2. The finding (freeze-9, ledger 2026-07-19 10:50)

- Accepted point, safe policy: `mixture=PASS m=0 alpha=0.0`, all guards True,
  `incident=48.60 effD=2.2498 dose=0.0400 arisk=1.046x mu_min=0.651`.
- Pure-OED endpoint, fast policy (m=972, rounded+repaired):
  `checks={budget:T, dose:F, ceiling:T, effD:T, arisk:F, spectral:F}`,
  `inc=583.2000 cap=583.2000 diff=-1.14e-13, effD=1.195, dose=0.1891,
  arisk=1.95x, mu=0.126`.
- The scan walked from m=972 to m=0: **at every m ∈ [1,972] at least one of
  {dose ≤ 0.05, A-risk ≤ 1.05×, spectral ≥ 0.5×} failed.**

Three convergent measurements of one structural fact:

1. RELAXED_REFERENCE gap at the dose-feasible continuous design: relref = 2.01/r
   (fast), 3.51/r (safe) — the dose-feasible OED sits far below the dose-relaxed
   optimum.
2. Pure SCAT32 FIXED_DOSE **beats** the dose-feasible continuous OED reference
   outright (safe: eff_D 2.25 vs 1.195-class rounded designs).
3. The rounded OED design fails dose/A-risk/spectral against the load-matched
   SCAT32 reference at every OED fraction above zero.

Physics reading: the R11 variable-load OED earns its Fisher advantage exactly by
concentrating incident load (in ρ and per-pixel dose). The guard stack — a ±5%
per-pixel dose band plus A-risk/spectral anchored to the strongest *dose-balanced*
baseline — forbids precisely that concentration. The dose-safe feasible set
collapses onto balanced-like designs, inside which SCAT32 is already (near-)optimal.
The vacuous arm is not a solver bug; the pure-OED endpoint's dose = 0.189 (3.8× the
band) and A-risk = 1.95× (vs 1.05× cap) show the frozen safety region and the OED
advantage region are **disjoint** at these operating points.

## 3. Box-10 interaction (full constrained certificate)

The remaining red box shares the same root. Diagnosis chain is closed:
freeze-8 `inf` = HiGHS time-limit on the dense degenerate 2×1024 μ block →
pixel-subset fix → freeze-9 `inf` = genuine LP unboundedness (μ-pixel with no
covering atom in the candidate set ⇒ t → −∞). Fix implemented (per-pixel coverage
seeding + finite μ cap; restricting the dual cone only loosens the bound, so
validity is preserved); verified on the R15 active-dose toy at 1.8e-9 ≤ 1e-8.
Real-scale re-validation was deliberately halted by the STOP. **The certificate's
target depends on R17's ruling** (certify against the dose-relaxed optimum vs the
dose-feasible optimum vs the balanced class).

## 4. Options we ask R17 to rank (and freeze wording for the winner)

- **(a) Widen the dose band for k-sparse atoms.** The A5 row-granularity analysis
  showed the rounding granule is comparable to the ±5% band at these loads, so part
  of the band violation is quantization, not physics. But the pure-OED endpoint is
  at dose 0.189 — no defensible widening (e.g. granule-scaled band) reaches it, and
  A-risk 1.95× stands regardless. Widening enough to admit OED (~2×) would gut the
  safety narrative.
- **(b) Re-anchor A-risk/spectral** to a load-matched but non-dose-balanced
  reference, or gate them on the continuous (pre-rounding) design. Same objection:
  anchors were chosen for the safety story; moving them post-hoc to rescue the
  method inverts the burden of proof.
- **(c) Accept the structural negative and INVERT the arm's claim (our
  recommendation to consider seriously).** The finding "under honest dose-safety
  guards, variable-load OED collapses onto the balanced fallback" is itself a
  clean, publishable result — and the box-10 machinery then carries the headline
  *as a certificate of near-optimality of the simple design*: SCAT32 at
  ridge-zone global load is within a certified gap of the best dose-safe design.
  The M1 arm keeps its 972-row budget but its preregistered claim becomes
  (i) primary: RIDGE-FIXED / SCAT32-at-ridge-load vs FIXED\* at matched dwell
  (the operating-map confirmation), and (ii) the certificate statement
  G_full ≤ ε for the deployed balanced design. No safety anchor moves; the
  "certified design" act of paper 2 becomes *certified near-optimality* instead
  of *certified adaptive win*. This also unifies with the R11 water-filling
  refutation (uniform servo −7 dB): pattern-level adaptivity was already dead;
  this closes the design-level version under safety, with a certificate.
- **(d) Split launch.** Freeze and launch the five non-gated arms now
  (OED-EQLOAD, SCAT16, SCAT32, LBLOB16, RIDGE-FIXED — none touches the guard
  stack), and re-derive the gated arm after R17. Gets confirmatory dwell curves
  ~1 day earlier; cost: two freeze events and a spec amendment note. Note
  OED-EQLOAD (un-gated OED) already provides the "what OED does without safety
  guards" contrast inside the campaign.

## 5. Questions for R17 (operative ruling requested, R13/R15 style)

- **Q1.** Is the m=0 outcome to be treated as (i) a spec defect to amend pre-launch
  (options a/b), or (ii) a genuine structural negative to be reported and
  re-architected around (option c)? We believe the three convergent measurements
  + the 1.95× A-risk endpoint support (ii), but we ask for the external ruling
  before touching a frozen guard.
- **Q2.** If (c): freeze the exact primary-endpoint wording for the re-architected
  arm — what is gated, against what comparator, with what pass bar — and the
  certificate statement's normative form (target class, ε tolerance, and whether
  the certificate is confirmatory or descriptive).
- **Q3.** Acceptance-rule audit: is descending-scan/first-pass the right mixture
  rule at all, or should the re-architected arm drop the mixture entirely
  (fallback-only + certificate)? Any residual role for the `alpha` disclosure?
- **Q4.** Box-10 target under the chosen architecture: certify the deployed design
  against the dose-FEASIBLE constrained optimum (our proposal: yes — that is the
  scientifically meaningful gap) with the coverage-seeded dual construction; state
  the acceptance tolerance (current toy verification 1.8e-9; proposed real-scale
  bar ≤ 1e-2 per R15).
- **Q5.** Launch policy: given (c) or (d), may the five non-gated arms launch under
  the current 13/14 ledger while the gated arm's re-derivation completes, or must
  the campaign wait for one coherent freeze? Preregistration integrity is the
  binding constraint; wall-clock is secondary but the operator has asked for
  results as soon as discipline allows.
- **Q6.** Paper-2 impact check: under (c), Act III's headline becomes "certified
  near-optimality of a simple dose-safe design + operating map," with the
  adaptive-OED negative reported alongside the servo refutation as a unified
  "adaptivity does not survive safety constraints" section. Confirm this framing
  is scientifically sound and OE-appropriate, or propose the correct framing.

Ruling delivery: as before, please post the R17 ruling as a GitHub issue on
`ccyyyYyzz/GI_a2` (DLP blocks chat extraction). Reference this document and
`results/round63_m1/R17_CONSULTATION_EVIDENCE.md` at this commit.
