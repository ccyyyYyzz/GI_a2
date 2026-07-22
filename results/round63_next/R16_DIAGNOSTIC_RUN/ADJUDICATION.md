# R16 §5.4 adjudication (Fable, 2026-07-22) — FROZEN

**Verdict: ESTIMATOR_ARTIFACT.** The residual −8% deviation between the
measured ridge locations (9.62 / 2.08, audited histogram estimator at
150k–400k frames) and the R16 predictions (10.22 / 2.17) is attributable
to peak extraction in the audited histogram estimator on ultra-flat
objective curves, not to physics missing from the jitter-capped ridge law.
The R16 predictions are CONFIRMED within Monte-Carlo uncertainty at the
frozen scale (N_MC = 2,000,000 frames/load, ν = 2000).

## Evidence chain (all from the same common-random-number draws)

1. **Exact anchor.** At c = 0 the ridge location is known exactly
   (ρ* = 22.2543). The score-identity estimator J_MI lands on it to
   <0.1% (22.2497, CI [21.999, 22.750]). The audited histogram estimator
   at the SAME draws misses by **+8.8%** (24.21 at dlog = 0.01) with a CI
   spanning nearly the whole grid — direct demonstration, against known
   truth, that the audited peak extraction carries an O(8%) location
   artifact of exactly the size under adjudication.
2. **CI separation.** The score CIs contain the predictions and exclude
   the measured values at both channels: cv=0.02 J_MI
   [10.100, 10.300] ∋ 10.22, ∌ 9.62; cv=0.20 J_MI [2.125, 2.200] ∋ 2.17,
   ∌ 2.08. The independent cross-check J_CS agrees (point 10.216 ≈ 10.22;
   MI–CS agreement ≤0.29% of J everywhere, satisfying §5.2).
3. **Bin-width smoking gun.** The histogram peak at cv=0.02 flips from
   10.10 (dlog = 0.0025 — agreeing with the score estimators) to 9.29
   (dlog ≥ 0.005): the artifact spans precisely the disputed deviation
   and vanishes as bin width → 0.
4. **Instability signature.** Histogram peak bootstrap sd is 5–10× the
   score estimator's (0.675 vs 0.096 at cv=0.02), and the audited peak
   drifted with sample size (9.645 @150k → 9.610 @400k → 9.288 @2M) —
   behavior of interpolation noise on a flat curve, not of a physical
   ridge displacement.
5. **Floor lever inert.** The α tail-floor changes nothing at this scale
   (tail mass ≤ 2.1×10⁻⁴; peak flat to the 4th decimal in α), eliminating
   the other hypothesized artifact channel.

## Honest bounds

The score point peaks sit ~1% below the predictions (10.10 vs 10.22;
2.15 vs 2.17) with CIs containing them; the c=0 anchor bounds the score
estimator's own bias at <0.1%. The confirmation is therefore at the
CI level with ≲1% residual consistent with MC noise; no claim of
agreement finer than that is made.

## Paper consequence

The paper-2 flag "the ~−8% residual displacement of the measured ridge
remains under diagnosis" is RESOLVED: the displacement was an artifact of
the audited histogram estimator's peak extraction; the jitter-capped
ridge law's predicted optima stand as stated. Wording installed in
paper2 by the coordinator (frozen-wording custody).

## Provenance

Runner `code/round63/jitter_score_diag_colab_frozen.py` (estimator
mathematics byte-identical to the audited prep driver; diffs = frozen
header, §5.3 bootstrap loop, heartbeat). Executed on two pro2 L4
sessions, PART_DONE both, VMs released. Full per-frame (N, L) records
are LOCAL-ONLY (603 MB):

```
raw_nl_cv0.0.npz   0f97eed0f0bcced3b11b1447442e917528436382cd8a5ec7f4cee3c8ea0f9d6d
raw_nl_cv0.02.npz  6516887cf6532dc9c4b4ba6121fb5e876496bbb5e6f30a3b4d4f70c2a8630d7c
raw_nl_cv0.2.npz   88b82ae21a16516a0f61aa5b75e1e7a9fdb559fd3efc8a1ec52626ba3e00c26a
```

Committed: SUMMARY_TABLE.md, r16_frozen_summary_combined.json,
r16_part_A.json, r16_part_B.json, run_part{A,B}.log, this adjudication.
