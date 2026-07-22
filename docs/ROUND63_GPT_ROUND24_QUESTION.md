# ROUND63 — GPT Round 24 (method round: from theory to a genuinely effective, deployable scheme — critique ours, propose better)

Operator directive (verbatim in spirit): "我要一个真正有效的方案。现在的理论
很深刻,但理论没有指引,别人看都不看" — deep theory without an effective
method will not even be read. The program stays at TWO papers; if a real
method exists it becomes the lead of paper 2 (the global power knob is its
rank-one special case), validated by its own preregistered campaign (M2).
This round happens BEFORE any implementation: critique our candidate hard,
and try to give a BETTER one. Nothing is frozen yet.

## 1. What "effective" means here (the acceptance bar)

A scheme another lab can adopt next month: inputs = detector parameters
(τ, jitter c), a cheap pre-scan, a pattern budget; outputs = operating
setting + illumination design + a per-scene guarantee; no per-lab
retraining, no fragile tuning; degrades gracefully to the one-knob policy
when adaptation is unnecessary. Plugin-like. Every theoretical component we
have is admissible material; none is sacred.

## 2. Our candidate: RGAI (Residual-Gated Adaptive Illumination)

Loop: (0) optional hardware step — clamp the hold-off register so c→0 by
construction (the anti-jitter menu's cheapest item); (1) 52-pattern
balanced pre-scan → x̂; (2) closed-form diagnostics from R23-T1: leverage
ℓ_i, marginal value h_i = ℓ_i J'(ρ_i), rank-one residual ε₁ (and the
alignment decomposition (α−α*_J)²C_u² + σ_ζ²); (3) ROUTE by failure mode:
ε₁ ≈ 0 → knob-only to the (clamped/jitter-aware) ridge (zero added cost;
the 19–21/24 scenes); σ_ζ/misalignment large → a few Frank–Wolfe geometry
steps over the admissible dictionary inside the full safety cone (engine =
our R18 support-expanding primal, which constructively found 3.7–4.8×
feasible D-eff headroom on every DEV family); C_u tiny → report the photon
tax C_u^{-2}, optionally lengthen dwell; (4) acquire + RQL. Self-certifying:
T1's bound f* − f_k ≤ ε_k²/(2m) ships WITH each output ("the method knows
when it isn't needed").

## 3. Known risks we see (be harsher if you can)

- **R1 (make-or-break): the Fisher→image gap.** R23 T5.3 is explicit:
  D-information gains are not PSNR corollaries. Our rescue probe already
  showed load-only control cannot fix the misaligned scenes (even a
  per-pattern genie failed); geometry steps are the untested limb. If
  FW-designed geometry does not lift image quality on misaligned scenes,
  RGAI collapses to the knob we already have.
- **R2: acquisition-time compute.** FW ascent took 60–240 s in our probes;
  a deployable loop wants seconds. Options: fewer steps, warm starts,
  precomputation (see Q2b).
- **R3: pre-scan estimate quality.** ℓ_i and the routing depend on x̂ from
  52 patterns; diagnostic noise → misrouting. Needs a stability story.
- **R4: estimator interaction.** λ_TV selection under nonstandard designs;
  the reconstruction layer may eat the design gains (our contour evidence
  hints the estimator is implicated).
- **R5: governance.** M1 verdicts are frozen and stand as reported; a new
  method requires its own preregistered campaign (M2) with fresh scenes.

## 4. Questions for R24

- **Q1 — kill or crown RGAI.** Attack the loop: failure modes we have not
  listed, prior art that already does this (adaptive/sequential OED for
  imaging, plug-and-play adaptive sensing, any published residual-gated
  design), and whether the self-certifying angle is genuinely novel or
  known in sequential-design literature.
- **Q2 — propose BETTER, if better exists.** Directions we explicitly
  invite you to weigh against RGAI (pick, combine, or replace):
  (a) estimator co-design instead of (or with) illumination design — put
  the adaptation in the reconstruction, keep illumination fixed;
  (b) **library routing instead of online optimization**: precompute a
  small library of k designs (from T3 cluster structure or scene-class
  priors), route by the diagnostic — one matrix multiply at runtime, no
  online solve; (c) amortization: train a tiny network to map pre-scan →
  design ONCE (offline, on synthetic scenes), keep the certificate as a
  post-hoc CHECK rather than the optimizer (deep learning as plumbing,
  certificates as guarantees — matches our earlier doctrine); (d) two-shot
  batch adaptation (pre-scan → one design → done) vs sequential; (e) any
  scheme of your own design. Criteria: universality across our three
  failure modes, deployability (R2/R3), guarantee strength, and survival
  probability against R1.
- **Q3 — the go/no-go experiment.** Specify the DEV image-level probe that
  decides R1 BEFORE any campaign: scenes (misaligned DEV instances),
  arms, dwell, seeds, the exact kill criterion (e.g. "if the adapted
  design's median PSNR gain over the better of SCAT32-060/RIDGE-SCAT32 on
  misaligned DEV scenes is < X dB, stop"), and what X should be.
- **Q4 — M2 sketch + paper-2 integration.** If the probe passes: the
  minimal preregistered M2 (scenes, endpoints, bars) that would let the
  method lead paper 2 without touching M1's frozen verdicts; and the
  paper-2 architecture with the method as the spine (law → knob as rank-one
  case → RGAI/your-better-scheme → M2 verdicts), within OE scope.
- **Q5 — deployment reality check.** DMD/SLM pattern-set switching costs,
  hold-off clamp availability on commercial SPAD modules, pre-scan photon
  budget accounting — anything that breaks the "plugin" claim in a real
  lab.

House rules: brutal candor; prior-art citations aggressively (if this
method exists under another name, we want the citation, not a rediscovery);
if RGAI (or your improvement) is NOT worth building, say so and why.
Ruling as a GitHub issue on ccyyyYyzz/GI_a2 referencing this document.
