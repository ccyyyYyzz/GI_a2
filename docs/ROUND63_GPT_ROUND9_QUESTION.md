# ROUND63 — GPT round-9 consultation request (post Study-2 verdict)

Repo: github.com/ccyyyYyzz/GI_a2. State at time of writing: tag
`study2-scored` (commit 9073425), paper filled through the verdict
(ab95c93). This document is the complete briefing; the questions are at the
end. Please answer as a formal ruling (R9) the way rounds 4–8 were ruled.

## 1. Frozen Study-2 outcome (no redesign follows — this is context, not a
question)

k=16 primary, DETAIL-32 cohort (24 scenes), frozen analyzer, 10k nested
family-stratified bootstrap:

- **DETAIL_SPEED_PASS = FALSE**, failed on breadth ONLY: 16/24 scenes with
  S_gate>1 vs preregistered 18/24. Effect-size criteria PASSED: median
  S_gate = 6.796 (bar ≥3), bootstrap LB2.5 = 2.835 (bar >1).
- Fixed-budget secondary PASSED all three criteria: median +1.16 dB
  (bar ≥1), 19/24 positive (bar ≥18), LB2.5 = +0.72 (bar >0). Per the
  frozen hierarchy it cannot rescue the primary.
- Censoring: NORMAL=17, SAFE_RANGE_UNINFORMATIVE=5, FAST_RIGHT_CENSORED=2.
- Family texture: glyph 4/4, maze 4/4 (S up to 33×, ΔQ up to +9.7 dB),
  spokes 4/4, chirp 3/4, contour 1/4, microtexture 0/4.
- Descriptive mechanism check (computed after unblinding, no decision
  role): per-family gate success is EXACTLY rank-ordered by measured C_u:
  maze .53 (4/4) > glyph .41 (4/4) > spokes .38 (4/4) > chirp .28 (3/4) >
  contour .27 (1/4) > microtexture .17 (0/4).
- Ladder fixed-budget medians (ν=2000): k=512 +0.65 (24/24), k=32 +4.17
  (24/24), k=16 +1.16 (19/24), k=1 +3.64 (21/24) — non-monotone; contrast
  vs conditioning (κ_A = n/k) compete; k=32 is the sweet spot.

Current paper wording (main_oe.tex, committed): abstract reports "formally
negative on scene breadth (16/24 vs 18/24) while the acceleration, where
measurable, was large (median 6.8×, LB 2.8×) and the fixed-budget
secondary passed"; results §6.3 reports all three criteria, the family
composition, and the C_u rank-order as an explicitly-labeled descriptive
observation; §6.4 uses "every rung positive, dense rung smallest,
non-monotone across sparse rungs" wording.

## 2. Questions

### Q1 — Verdict presentation calibration (paper-blocking)
Given the outcome above: (a) does the current title ("When high flux helps
single-pixel imaging: a contrast–dead-time phase diagram") remain the
right claim level, or does 16/24 force a further hedge? (b) Is the
"phase boundary passes through the confirmatory cohort, gate outcome
monotone in measured C_u" framing review-safe as a descriptive claim, or
does it need a preregistration-status caveat beyond the one sentence we
have? (c) Any wording traps you see in abstract/§6.3 as summarized above?
(d) Should the k=32 sweet-spot observation be promoted into the abstract
(currently only ladder medians are there)?

### Q2 — Scene-side contrast: scale-matched patch-sparse patterns
The 8 failing scenes are exactly the low-C_u ones. Mechanism: with k
SCATTERED on-pixels, bucket variance of fine-scale texture is √k-averaged
away (C_u ≈ σ_pix/(√k·mean)). Proposal to evaluate (idea-space, likely
future work): keep k fixed but make the on-support a CONTIGUOUS patch
whose size matches the scene correlation length; generalize to a frozen
multi-scale family (wavelet-like ladder of patch scales). Questions:
(a) is the √k-averaging analysis correct as stated? (b) prior art check —
does patch/block/superpixel sparse illumination for CONTRAST (not for
resolution/compressive reasons) exist in SPI literature? (c) does this
belong in THIS paper as a discussion paragraph only, or is any inclusion
beyond one sentence a scope error?

### Q3 — Is adaptivity provably necessary?
Conjecture: for every FIXED pattern ensemble satisfying the equal-load
constraints, there exists a structured (non-flat) scene with Γ < 1,
whereas an oracle-adaptive designer (allowed a coarse pre-scan) achieves
Γ > 1 on every structured scene. (a) Is this provable / already known
(adversarial scene vs fixed dictionary)? (b) Our earlier linear-Gaussian
no-adaptation lemma says adaptivity gains nothing asymptotically in that
layer — the dead-time channel is nonlinear, so the lemma's premise fails;
is "adaptivity strictly helps under dead time" a defensible theorem
target, and of what difficulty? (c) If provable, is it a THIS-paper
theorem or the anchor of the follow-up?

### Q4 — Optical common-mode cancellation (the Study-3 / collaboration
pitch)
Nonnegative intensity projection forces DC ≫ modulation on low-contrast
scenes; dead time taxes total counts, so the DC eats the budget.
Interferometric proposal: interfere the bucket field with a reference
shaped to the mean scene; balanced two-port output — "difference" port
(near-zero rate, all structure) and "sum" port (DC monitor). Dead-time
budget at the difference port is then spent almost entirely on structure;
effective C_u → large for ANY structured scene. (a) Physical soundness
check — what breaks first (coherence through scattering media, phase
stability, shot noise of the cancellation itself, reference-shaping
error)? (b) Photon statistics at the difference port: still Poisson to
adequate approximation for our renewal analysis? (c) Prior-art check:
does "dark-field"/null-interferometric single-pixel imaging exist? How
does this relate to complex-field SPI (Hao & Chen, Opt. Lett. 2025)?
(d) Is this the right anchor for the PolyU collaboration pitch, and does
it merit a forward-looking paragraph in this paper's discussion?

### Q5 — Remaining scope triage for the OE manuscript
Pending: NATURAL-24 cohort analysis (Colab tail finishing), S2B/S2C
ensembles+scale, S3 mismatch map, S4 exact-likelihood small-scale split.
For an OE-length two-study paper, which of these are load-bearing
(reviewer-expected) vs supplement vs droppable? Specifically: is S4
(EXACT+TV arm at 8²/16²) still worth its cost now that the two-study
structure carries the argument?

Please rule R9 on Q1–Q5. Q1 and Q5 block the manuscript; Q2–Q4 are
idea-space (no execution is planned before your ruling and user sign-off).
