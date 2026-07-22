# R28 — Gate A oracle-arm definition under a discovered dose-band/FW incompatibility (fast ruling requested)

Context: all R27 launch conditions were executed — k_eff≥32 library rebuilt,
8/8 banks pass exact-972 admission with witnesses (independently verified,
manifest @ 2c11130), allocator caching verified output-equivalent, smoke-3
clean for ORACLE-LIB and RLMI (all genuine materialization). Smoke-4 on the
online FW arms then surfaced a finding that makes frozen Gate A unevaluable
as written. Grid remains HOLD. This round needs exactly one decision.

## The finding

Per the R27-§1.1-consistent implementation ruling (disclosed below), the
online primal_probe dictionary for TRUE-X-FW/XHAT-FW was restricted to the
same k_eff≥32 super-atom dictionary as L7, and each per-cell FW design is
materialized by the proven T-B witness construction (coverage-aware
selection + deterministic Sinkhorn nominal amplitudes, all frozen guards
enforced post-hoc). The construction WORKS — but the result:

| scene | arm | PSNR_rad | mat | dose_dev |
|---|---|---|---|---|
| twopop_0 | SCAT32-060 baseline | 22.78 | — | — |
| twopop_0 | TRUE-X-FW (materialized) | **9.66** | OK | 0.021 |
| twopop_0 | XHAT-FW (materialized) | **9.52** | OK | 0.047 |
| control_0 | ORACLE-LIB (best = L7) | **34.60** (c=0.05) | OK, all 8 banks genuine | 0.021–0.040 |
| control_0 | SCAT32-060 | 30.25 | — | — |

Sanity checks isolating the cause:
- Witness-materializing a KNOWN-GOOD SCAT32 support → 20.8 dB (−2 dB from
  tuned SCAT). The construction is sound.
- Re-biasing row powers toward the FW weights recovers 27.8 dB but breaks
  the dose band (0.075 > 0.05). The FW design's value LIVES IN nonuniform
  per-pixel dose — inherently incompatible with the ±5% uniformity guard.
- Fragility: the same FW design materializes at support 7278 (fw_iters=10)
  but falls back at support 11000 (fw_iters=40) — sensitivity to support
  size, second symptom of the same incompatibility.

So under the frozen deployment class (±5% dose band = the safety/resource
constraint every arm honors), an unstructured per-pixel FW oracle is not a
usable reachability witness: dose-uniformizing its support yields badly
conditioned designs (9.5 dB ≪ 22.8 baseline), while the information
weighting that made it an oracle is exactly what the class forbids. Gate A
as frozen (R24 §3.4: TRUE-X-FW ≥1 dB median gain over the fixed composite
baseline, 9/12, LB>0.5; fail kills the method line and M2) would fail at
−13 dB for reasons that have nothing to do with the presence or absence of
image-level adaptivity headroom — which ORACLE-LIB demonstrates at +4.3 dB
under EXACT materialization in the same class.

## Reframing observation (for your consideration, not a claim)

The ±5% band quotients out "concentrate dose on informative pixels" —
trivial, often unsafe adaptivity. What remains admissible is GEOMETRY
adaptivity: which structured banks to mix. ORACLE-LIB is then precisely
the reachability ceiling of the admissible class, and RLMI is its
deployable estimator; the FW arms' collapse is itself evidence that
per-pixel weighting does not survive the class. Under this reading the
method story sharpens rather than weakens.

## Options to rule on (we recommend (c))

(a) Re-anchor Gate A on ORACLE-LIB (same thresholds: ≥1 dB median over the
    composite baseline, 9/12, LB>0.5), demote both FW arms to descriptive
    gap diagnostics under the band (their honest materialized numbers
    reported; no gate reads them).
(b) Exempt the nondeployable TRUE-X-FW oracle from the dose band
    (disclosure-only), preserving Gate A's original Fisher-headroom
    meaning — touches the frozen guard profile of a gating arm.
(c) BOTH ROLES SPLIT: Gate A re-anchored on ORACLE-LIB as the
    deployable-class reachability gate (thresholds unchanged); the FW arms
    run UNCONSTRAINED by the band as descriptive Fisher-ceiling arms
    (clearly labeled nondeployable; band replaced by disclosure of their
    realized dose profile; no gate reads them). Gate B unchanged (RLMI
    ≥60% of ORACLE-LIB gain + no-harm) — note A and B then form a coherent
    pair: A = is the class ceiling ≥1 dB real; B = does RLMI capture ≥60%
    of it. Gate C, PLUGIN_LATENCY_PASS, four-gap decomposition unchanged.
    The gap-2 (Fisher→image) diagnostic becomes FW-unconstrained minus
    ORACLE-LIB.

If (c): confirm the exact frozen wording for Gate A's new definition and
whether the twopop_0 exposure during smoke-4 (materialization validation
only; no gate statistics computed; 14 stress scenes still blind) requires
any scene replacement or is acceptable as-disclosed.

## Disclosures for ratification (made under R27 execution, documented in
docs/ROUND63_BRIDGE_BUILD_PLAN.md)

1. Bank admission witnesses use each bank's FROZEN per-atom nominal power
   schedule (Sinkhorn amplitudes within w∈[1/4,4], peak≤1536) — read as
   compliant with your §2 text "under the bank's frozen nominal powers"
   (L0 itself is grandfathered with a ridge schedule). Strict equal-weight
   would instead exclude L1/L5 (they stall at 0.0535).
2. Online FW dictionary restricted to k_eff≥32 (conservative direction:
   weakens the oracle only). Made Gate A evaluable in principle; the
   incompatibility above is what remained.

Grid is otherwise fully green: 8/8 banks admitted, equivalence verified,
unit suites 13/13 + 8/8, cost ~10 CPU-hr across 5 routes. Speed matters —
one decision unblocks launch.

Deliver as a GitHub issue on ccyyyYyzz/GI_a2 referencing this document.
