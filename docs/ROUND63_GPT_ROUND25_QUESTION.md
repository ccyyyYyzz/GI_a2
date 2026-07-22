# ROUND63 — GPT Round 25 (router redesign: dissolve the threshold gate into a convex bank allocation)

Context: T-A (16 bridge scenes) is frozen; T-B (the 8-bank library) is in
construction; T-C (the router) is NOT yet built — this is the moment to fix
the architecture's least principled part. The operator's critique of the
R24 Step-3 gate, verbatim in spirit: "鲁棒分诊很不自然,很刻板" — the
three hand-set thresholds (U1/r ≤ 0.01, 1.02 knob preference, 2% margin),
the S=16 max-envelope, the hysteresis, and the categorical
KNOB/LIBRARY/ABSTAIN switch feel like a rulebook approximating something
that should be a single principled object. We agree and propose the
replacement below. Attack it.

## 1. Proposal: minimax-convex allocation over the bank simplex

Replace R24 §2.3 Step 3 entirely. With preloaded banks L_0..L_{K-1}
(K=8), per-bank predicted Gramians F_k(x) at the candidate operating
points, and S posterior/bootstrap draws {x_s} around the pre-scan
estimate, solve

    w* = argmax_{w in Delta_K}  min_{s=1..S}  Phi( H_0(x_s) + sum_k w_k F_k(x_s) )

with Phi = logdet (or -tr[W(.)^{-1}] to keep R24's estimator-aware risk),
then realize the design by drawing round(972·w_k) rows from bank k
(largest-remainder rounding).

Claimed properties (verify or refute each):
- P-A: the inner objective is concave in w for Phi = logdet (linear map of
  w into PSD cone; min over s preserves concavity), so this is a small
  convex program solvable in milliseconds at K=8.
- P-B: thresholds/hysteresis/abstain DISSOLVE: the solution is a
  continuous function of the data (fixes R24-F2 discontinuity by
  construction); under draw disagreement the minimax solution hedges
  toward banks robust across draws — the "abstain to robust default"
  behavior emerges as a property of the solution rather than a rule;
  knob-only emerges as the corner w = e_0 when evidence warrants.
- P-C: guard inheritance by convexity: each bank satisfies the ±5% dose
  band and budget; proportional row-mixing gives dose = convex combination
  → in-band; budget linear; peak per-row → mixture-safe up to rounding
  perturbation (to be re-verified post-materialization per F5).
- P-D: a sharper certificate: exact KKT verification over the
  K-simplex (per-draw active sets), replacing the fuzzy "route margin".
- P-E: deployment unchanged: banks preloaded; a mixture is an interleaved
  subsequence of preloaded patterns; Gate A–D of the bridge unchanged
  except the RGLI arm's internals.

## 2. Questions

- **Q1.** Verify/refute P-A..P-E. In particular: (i) does -tr[W(.)^{-1}]
  (the R24 estimator-aware risk) preserve concavity in w? (it is convex in
  the information matrix argument being inverted — check the composition
  direction carefully and give the correct statement); (ii) does the
  minimax-over-draws formulation have failure modes the max-envelope gate
  does not (e.g. one adversarial draw dominating; S small; draws not
  covering the true scene)?
- **Q2.** The correct robust objective: minimax over S draws vs
  CVaR/soft-min vs posterior expectation with a variance penalty — which
  is right HERE, given S=16 parametric-bootstrap draws of a 52-pattern
  pre-scan, and given that we want continuity AND conservatism without
  new tuning knobs? If a temperature/level parameter is unavoidable, is
  there a principled default (e.g. minimax exactly = level 1/S)?
- **Q3.** Row-mixing subtleties: rounding 972·w_k across 8 banks;
  correlation structure when banks share pattern families; whether mixing
  degrades any single bank's design coherence (e.g. L7's FW bank was
  optimized as a WHOLE — does taking a fraction of its rows preserve its
  value, or must banks be constructed as mixable atoms from the start?).
  This is the strongest objection we see — rule on it and, if needed,
  prescribe the bank-construction change (e.g. build banks as exchangeable
  row-families rather than jointly-optimized sets).
- **Q4.** Disclosure and bridge integration: what replaces the
  KNOB/LIBRARY/ABSTAIN route fractions in Gate C (routing informativeness)
  — e.g. weight entropy, distance of w* from e_0, per-draw allocation
  spread? Freeze the Gate C reinterpretation so the bridge remains
  decision-valid without new tuning.
- **Q5.** Net ruling: adopt the convex allocation as RGLI Step 3, keep the
  threshold gate, or a hybrid (allocation + one honesty flag when the
  minimax value is nearly w-flat)? One recommendation, with the same
  brutal-candor standard as R24.

Ruling as a GitHub issue on ccyyyYyzz/GI_a2 referencing this document.
