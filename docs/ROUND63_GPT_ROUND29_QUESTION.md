# R29 — Independent software-method generation (think first, then judge ours)

Mandate change from the operator (verbatim constraints):
1. **Software only.** No hardware innovation of any kind. Shelved: hold-off
   clamping, timing telemetry, source self-gating, comparator circuits.
   "Software" = anything implementable as code or settings on the existing
   simulated optical system: illumination pattern CONTENT and ORDER (SLM is
   programmable), reconstruction, calibration, estimation, scheduling.
   No new physical components, no detector modifications.
2. **No manuscript is written until a software METHOD positive exists.**
   The M1 ridge-tracking positives (+1.87 dB / 19.13×, frozen) do not meet
   the operator's bar for the method act — the demand is a genuinely
   effective, deployable software method, theory-guided, validated at
   IMAGE level in a preregistered campaign with kill bars.
3. Judging criteria: 创新 (novelty), 简洁 (simplicity), 普适 (generality),
   美观 (elegance). Deep learning admissible as TOOL only (racing a
   certified ceiling), never as the claimed innovation. Target problem
   space: ghost imaging broadly, incl. dynamic/scattering media and
   high-quality reconstruction.

## The frozen walls (do not propose things these kill)

- **Count-only impossibility (proven)**: no post-processing of counts
  alone recovers the dead-time-jitter information loss (E1 + CRB
  argument). Reconstruction cannot ADD information in ch.1.
- **Bridge negative (2026-07-22, results/round63_bridge/BRIDGE_VERDICT.md
  @ 026c239)**: within the ±5% dose-uniform class at the hard cell, the
  8-bank finite-library oracle clears the knob baseline by only
  +0.68 dB median (bar was 1.0); k* was ALWAYS L0/L1 — structured
  geometry never beat scattered; a robust maximin allocator on noisy
  pre-scans HARMED (−8.3 dB median allocation loss); the dose-relaxed
  Fisher/A-risk FW oracle sits −12.9 dB BELOW the simple library at image
  level (surrogate-to-image gap, measured). Lesson binding on every new
  candidate: Fisher-side superiority must be demonstrated at image level.
- Prior negatives from this program (do not re-propose): CSL, CSNT
  deployment, CFPC, CBFSA, CR-MFD structural probes; runtime-threshold
  triage mechanisms (operator finds them unnatural/rigid); third papers.

## The assets a method may lean on

- E1 identity (info loss = E[Var(hidden | record)]) and its ch.2 twin O1
  (I_Y(x) = MᵀD M − MᵀE[Cov(a|Y)]M, hidden gain-weighted illumination
  vector G = Mᵀa) — both proven, numerically verified.
- Jitter-capped ridge law + operating maps (blind-verified; R16 closure).
- O4-A optimal-scheduling moment conditions (MᵀR⁻¹D_sH = 0) — proven;
  numerical: time-ordered schedule violation metric 1770, random 63.6,
  paired 0.0 (results/scatter_verify/).
- Score-identity estimator machinery validated to <0.1% at an exact
  anchor (results/round63_next/R16_DIAGNOSTIC_RUN/).
- Certified-ceiling discipline: any learned/heuristic component races a
  computable information bound, so claims stay honest.
- Range-null accountability machinery (two-ledger audit) from the
  operator's companion program: null-space content of reconstructions is
  attributable/auditable — a hallucination audit for DL-as-tool.

## What we already generated (OUR candidates — judge them AFTER
## generating your own; do not anchor on them)

- **S1 Pattern scheduling under gain drift/dynamic media**: pattern ORDER
  is free software; O4-A constructs drift-immune schedules; Probe B
  (image-level ordered/random/paired, 5 paired seeds, frozen drift model)
  is running now.
- **S2 Estimator-efficiency gap**: measure RQL reconstruction efficiency
  against the certified ceiling at the deployed operating points; if a
  large gap exists, close it (possibly DL-as-tool with the ceiling as
  yardstick + null-space audit). Probe A running; if ≥90% efficient the
  lane dies.
- **S3 TDC-free self-calibration**: identify (τ, c) purely from counting
  statistics (mean + Fano vs drive), no timing hardware; enables
  operating maps on uncharacterized systems.

## The R29 ask

1. **Generate independently 5–8 software-method candidates** meeting the
   mandate (software-only, image-level positive plausible, theory-guided,
   novel vs prior art). Think from the identities/assets, from the
   bridge's mechanistic lessons (e.g. WHY did the allocator harm? the
   answer may itself define a method), and from the Chen-group problem
   space (dynamic scattering media, high-quality reconstruction). For
   each: mechanism, why theory predicts a positive, image-level test
   design sketch, the closest prior art you can name, and the kill risk.
2. **Rank all candidates (yours + S1–S3)** by
   P(preregistered positive) × novelty × simplicity × generality.
   Be ruthless about prior-art renamings — this program has killed 15+
   candidates for that before.
3. For your **top pick**: a preregistered campaign sketch (endpoints,
   kill bars, frozen comparators, scene/drift cohort, seeds) in the style
   of the M1/bridge campaigns — small enough to run in days on
   CPU/5 Colab routes.
4. Flag explicitly if you believe NO candidate clears the bar — a
   "the software space at this layer is saturated" verdict is admissible
   and would redirect the program (that itself is valuable).

Deliver as a GitHub issue on ccyyyYyzz/GI_a2 referencing this document.
Reference commit will be stated in the chat message. Depth matters more
than speed this round: this decision sets the next multi-week campaign.
