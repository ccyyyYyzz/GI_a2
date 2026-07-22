# DEV bridge verdict (2026-07-22) — FROZEN

**LIBRARY_REACHABILITY_FAIL → RLMI/M2 STOP (R28 §1; R24 hard stop — no
rescue redesigns, no re-runs, no scene swaps).** The method line ends here
as a preregistered negative. The flagship manuscript ships with the
law/knob arc as primary and this bridge as its honest closing act.

Grid: 320/320 cell-groups, 1920 arm records, 5 local shards, 117 min,
zero skips/errors. Evaluation: frozen `bridge_gates.evaluate_all` @
merged main (commit 384315a) via `evaluate_gates_run.py` (schema-flatten
adapter only — field promotion, no value changes). Full output:
`GATE_VERDICT.json`.

## Gate A — LIBRARY_REACHABILITY_PASS = FALSE (the deciding gate)

Hard cell (c=0.05, ν=2000), 12 stress scenes, 5 paired seeds:

- median ΔQ^A = **+0.680 dB** — REQUIRED ≥ 1.0 → **FAIL**
- positive scenes 11/12 — required ≥9 → pass
- 90% bootstrap LB **+0.543 dB** — required >0.50 → pass

Two of three conditions passed; the median missed the frozen bar. The
declared 8-bank library helps almost everywhere but not enough: the
finite-library image oracle's advantage over the composite knob baseline
is real, scene-dependent, and small.

Per-scene ΔQ^A (five-seed means; k* = winning bank):

| scene | ΔQ^A (dB) | k* |
|---|---|---|
| contour_3 | +2.051 | L1 |
| contour_1 | +1.256 | L1 |
| twopop_3 | +1.052 | L0 |
| microtex_3 | +0.998 | L1 |
| contour_0 | +0.995 | L0 |
| microtex_1 | +0.697 | L1 |
| contour_2 | +0.663 | L0 |
| microtex_0 | +0.583 | L1 |
| twopop_1 | +0.502 | L1 |
| microtex_2 | +0.498 | L1 |
| twopop_2 | +0.002 | L0 |
| twopop_0 | −0.085 | L1 |

Notable: k* is ALWAYS L0 (the deployed SCAT32 knob design) or L1
(scattered-k32) — none of the structured banks (rings, bars, multiscale,
compact pairs, FW-geometry) ever wins a scene. Within the dose-uniform
class, geometry diversity beyond scattered supports bought nothing.

Leave-smoke-exposed-out sensitivity (disclosure-only): median +0.697,
same FAIL — the verdict is robust to the smoke exposure.

## Gate B / Gate C (moot given A, recorded for completeness)

- Gate B FAIL: RLMI median −7.45 dB vs base; 3 controls below −1 dB
  (control_1 −19.8, control_3 −26.1). The deployable allocator not only
  failed to capture the (sub-bar) oracle gain — it actively harmed on
  5/12 stress scenes and 3/4 controls.
- Gate C FAIL: mean A_j = 0.53 on rescue scenes (bar 0.80), mean ŵ₀ = 0.50
  on aligned controls (bar 0.75) — the allocator deviated from the knob
  about half the time under real pre-scan noise, and those deviations are
  where the harm lives. (Contrast: in noise-free smoke it returned w=e₀
  on controls. The S=16 bootstrap scenario machinery did not protect
  against pre-scan-noise-driven misallocation at ν=2000.)
- PLUGIN_LATENCY_PASS = FALSE (median 13.5 s, p95 30.5 s) — moot; claim
  label was already batch/two-shot per R27/R28.

## Four-quantity decomposition (R28 §4)

1. in-class finite-library reachability: **+0.680 dB** (Gate A, above)
2. deployable allocation loss: **−8.33 dB** (RLMI − ORACLE-LIB)
3. plug-in FW loss (descriptive, dose-relaxed): −1.73 dB (XHAT − TRUE-X)
4. dose-relaxation/library contrast (descriptive): **−12.91 dB**
   (TRUE-X-FW − ORACLE-LIB) — even with the dose band removed, the local
   Fisher/A-risk FW surrogate materializes to images far below the simple
   library designs; with TRUE-X − base = −11.31 dB (nondeployable label).

`no_gate_reads_fw = True` verified on the final evaluation.

## Scientific reading (for the manuscript's closing act)

Within the safe deployment class (±5% per-pixel dose uniformity), at the
hard operating cell, the preregistered conclusion is:

1. **The global ridge knob captures nearly everything this class offers.**
   The best of 8 exactly-materialized geometry banks clears the knob
   composite by a median of only 0.68 dB (11/12 positive, LB 0.54) —
   below the 1 dB deployment bar we preregistered as "worth adapting for".
2. **Structured geometry never beat scattered.** k* ∈ {L0, L1} on all 12
   scenes.
3. **Per-pixel information weighting fails at image level in this regime
   even when the dose constraint is removed** (q4 = −12.9 dB): the
   Fisher-style surrogate asks for designs whose image-domain
   reconstructions are poor — the surrogate-to-image gap, measured.
4. **A robust maximin allocator on noisy pre-scans caused harm** where the
   oracle margin was thin (allocation loss −8.3 dB median) — adaptivity
   at ν=2000 pre-scan SNR is not merely unhelpful, it is risky.

Together with the proven software-side ceiling (count-only correction
impossibility), this completes the story: at high flux under dead time +
jitter, neither software post-correction nor safe-class illumination
adaptation clears the deployment bar — **the certified global operating
law (ridge tracking, +1.87 dB / 19.13×, M1) is the practical optimum of
this regime**, now bracketed from both sides by preregistered evidence.

## Decision consequences (frozen tree)

- M2: DEAD (never preregistered, never run).
- RLMI: reported as the preregistered negative + descriptive diagnostics;
  does not lead any manuscript section as a deployable method.
- Flagship one-paper consolidation proceeds: Acts I–III (identity → law →
  campaign) + Act IV = this bracketing negative (its strength is the
  completed impossibility-both-sides argument).
