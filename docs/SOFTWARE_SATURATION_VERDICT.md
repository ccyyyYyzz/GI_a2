# Software-layer saturation verdict (2026-07-22) — program-level, frozen

Per R29 §6 (issue #21), the pre-written conclusion is now in force:

> **No current software candidate clears the operator's
> novelty × simplicity × generality × image-positive bar.**
> The software layer for these frozen channels is saturated at the
> requested effect size (≥1 dB deployable image gain).

## The five preregistered certificates behind it (all frozen, one day)

| # | lane | mechanism tested | verdict | governing number |
|---|---|---|---|---|
| 1 | RLMI bridge (Gate A) | scene-adaptive geometry from noisy pre-scan, dose-uniform class | FAIL | +0.68 dB median < 1.0 (@ 026c239) |
| 2 | RLMI allocator (Gate B/C) | robust maximin routing | HARM | −8.3 dB allocation loss |
| 3 | DOPS scheduling | O4-A moment-orthogonal pattern order | FAIL | +0.039 dB median < 1.0 (@ c592a4e) |
| 4 | Estimator efficiency (Probe A) | RQL vs certified Fisher ceiling | SATURATED | ~99% efficient; ≤0.2 dB headroom (@ 69207b3) |
| 5 | CPL-GI (Probe C / reserve) | exact conditional-likelihood gain elimination | KILL | +0.412 dB < 1.0 (@ 73b33ab) |

Supporting descriptive facts: the dose-relaxed FW oracle sits −12.9 dB
below simple designs at image level (surrogate-image gap); the known-gain
oracle sits 10–15 dB above ALL software arms (the headroom exists but is
reachable only by KNOWING the hidden state, i.e., by a richer record).

## Why (one sentence, the identity's own explanation)

Information lost = E[Var(hidden state | record)]; the bucket-count record
is too thin to identify the hidden state, so every software post-hoc
trick is bounded by that conditional variance — and five independent
mechanisms measured the bound as binding.

## What this does NOT say

- It does not kill the science already frozen (M1 +1.87 dB / 19.13×,
  ridge law, identities, ceilings — all stand).
- It does not kill software methods under a RICHER observation model
  (e.g., time-tagged records that existing photon-counting hardware
  already produces) or in a DIFFERENT application/channel.
- It does not kill prior-side work: Probe A explicitly locates the
  remaining image error in near-null directions — prior territory,
  uncertified by Fisher, where accountability (range–null audit) rather
  than "recovery" is the honest deliverable.

## Frozen consequences

1. No further method campaigns on the bucket-count channels at this
   effect size; no relabeled corrections.
2. Manuscript status: per the operator's standing order (no writing
   without a method positive), the flagship consolidation remains ON HOLD
   pending the operator's strategic decision on program direction.
3. The next direction decision is the OPERATOR's (scope fork): richer
   record (software on time-tagged data), application pivot, the
   accountability/audit line, or accepting the bracketing manuscript.
