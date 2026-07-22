# DEV bridge build plan (Fable-owned; R24 issue #16 is the normative spec)

**R27 amendment layer (issue #19, docs/ROUND63_GPT_ROUND27_RULING_RAW.md —
overrides anything below that conflicts):**

1. **k_eff >= 32 for every physical main-acquisition atom.** All 16-support
   mixture entries rebuilt as support-32 counterparts; multiscale minimum
   physical support 32; L7 re-solved on the restricted dictionary; L0
   unchanged. 16-pixel motifs only inside manifest-frozen disjoint-pair
   super-atoms (one simultaneous 32-pixel exposure, scored with its actual
   load/dose/peak/Gramian). Post-allocation pairing forbidden.
2. **Bank admission before allocation (R25 §9 frozen replacement, R27 §2):**
   every bank must demonstrate, BEFORE the grid, at least one exact 972-row
   realization satisfying budget + ±5% dose band + peak/admission — witness
   stored in the manifest. Failing banks are excluded (NONMIXABLE), not
   rescued. The γ∈[0.94,1.06] trim is polish, not a feasibility certificate.
   On materializer failure at runtime: return L0 + MIXTURE_MATERIALIZATION_FAIL;
   never loosen band/row-count/peak.
3. **Gate D renamed PLUGIN_LATENCY_PASS — claim classification only.**
   Gates A–C = science go/no-go. D-FAIL does not block M2: the method is
   then described as batch/two-shot adaptive illumination (never
   real-time/subsecond/plugin) with median/p95/max/CPU/threading/cache-state
   disclosed. No substitute threshold.
4. **Output-equivalent caching authorized pre-grid** (H0_s/W_s/F_{k,s}, R*_s
   oracles, cross-arm reuse on identical pre-scan state, batched Cholesky,
   warm starts, exact Woodbury): must match reference on all phase-1.5 cells
   within |t_fast−t_ref|≤1e-8, ||w_fast−w_ref||₁≤1e-6, KKT ≤ frozen tol,
   identical 972-row realization after tie-break. r=200 and S=16 must NOT
   be reduced.
5. **RIDGE-SCAT32 degradation at c=0.05 is physics, not a bug** (R27 §4) —
   preserved; composite baseline absorbs it per R24 §3.3.
6. **Launch HOLD until**: rebuilt library + manifests regenerated; every bank
   passed exact-972 admission or excluded; optimized allocator passed
   equivalence; phase-1.5 smoke reruns clean.

**R28 amendment layer (issue #20, docs/ROUND63_GPT_ROUND28_RULING_RAW.md —
supersedes Gate A and the FW-arm guard profile below):**

1. **Gate A → LIBRARY_REACHABILITY_PASS** (finite-library image-reachability
   gate): ΔQ^A_j = ORACLE-LIB − max(SCAT32-060, RIDGE-SCAT32), per-scene
   oracle bank = argmax over 5-seed means (min-index tie-break, one bank
   per scene); PASS iff median ≥ 1.0 dB, 9/12 positive, 90% scene-bootstrap
   LB > 0.50 dB. Fail ⇒ RLMI/M2 stop (statement scoped to the DECLARED
   library, not the physical class). Terminology: "finite-library image
   oracle" — never "class ceiling".
2. **Gate B** capture denominator = ORACLE-LIB (median ΔQ_RLMI ≥ 0.60 ×
   median ΔQ_ORACLE-LIB); all else unchanged.
3. **FW arms = NONDEPLOYABLE_DOSE_RELAXED_FW_DIAGNOSTIC**: ±5% dose band
   removed (ONLY that); all other guards + k≥32 dictionary retained; full
   dose-profile disclosure; no gate reads them.
4. **Four quantities** (§4): in-class finite-library reachability;
   deployable allocation loss; plug-in FW loss; dose-relaxation/library
   contrast (never labeled a pure gap).
5. **SMOKE_EXPOSED_PREGRID = {twopop_0, control_0}** — retained, disclosed;
   post-grid leave-out sensitivity (no-gate).
6. Both implementation disclosures **ratified** (§6): pattern-and-power
   banks (nominal schedules are bank identity); k≥32 online dictionary.
7. Launch authorized once §7 checks pass; bank manifest NOT regenerated.

**Fable implementation ruling 2026-07-22 (post-smoke-3, within R27 — no
guard touched; §6-ratified by R28):** the ONLINE FW designs (TRUE-X-FW / XHAT-FW oracle arms)
are subject to the same arithmetic obstruction R27 §1.1 identified for
banks; without a fix they fall back to L0 and Gate A becomes vacuous.
Remedy = R27's own: (1) online primal_probe dictionary restricted to the
k_eff>=32 physical super-atom dictionary (identical to L7's) — conservative,
can only weaken the oracle, cannot inflate Gate A; (2) per-cell exact-972
realization via T-B's witness construction (coverage-aware selection +
Sinkhorn nominal-amplitude balancing within the peak guard), deterministic
given (measure, dictionary), all guards verified post-hoc per the R25 §9
replacement text; (3) L0 + MIXTURE_MATERIALIZATION_FAIL remains the honest
fallback. Option "exempt oracle from the dose guard" REJECTED (frozen-guard
change). Nominal-power-schedule witnesses ruled COMPLIANT for banks (R27 §2
verbatim presupposes "the bank's frozen nominal powers"; L0 grandfathered
with a ridge schedule); strict equal-weight NOT required. Disclosed here
for the post-grid GPT report.

Command change 2026-07-22: the build is owned directly by the main agent
and executed by decomposed Opus subagents. The docs2/ letters 13–14 are
superseded as assignments (their technical content, including the
on-branch-evaluation requirements from the falsifier addendum, is folded
here). Output root: `results/round63_bridge/`. Scenes:
`data/r63_bridge_scenes/`.

## Interface contract (fixed now so tasks parallelize)

- **Scene** files: `data/r63_bridge_scenes/{scene_id}.npz` with field `x`
  (32x32 float, [0,1]) + `{scene_id}.png` preview; manifest
  `data/r63_bridge_scenes/BRIDGE_SCENES.json`: list of
  {scene_id, group (contour|twopop|microtex|control), gen_params, sha256};
  seed base 650000+ (disjoint from 632900 DEV and 633000 confirmatory).
  FROZEN (committed) before any arm evaluation.
- **Bank** files: `results/round63_bridge/library/L{k}.npz` with fields
  `rows` (972 x 1024 nonneg), `meta` (json string: family, params,
  guard-check results, sha256 of rows); index
  `library/LIBRARY_MANIFEST.json`. Guards re-verified AFTER
  materialization (R24 F5): budget, peak (w in [1/4,4] & max<=1536),
  +-5% dose band, admission, exact row count.
- **Arm outputs**: `results/round63_bridge/cells/{scene}_{arm}_{nu}_{c}_{seed}.json`
  (PSNR_rad, PSNR, achieved loads q5/50/95/max, incident, detected,
  route metadata for RGLI, wall times). Kill-gate evaluation reads only
  these.

## Task decomposition

- **T-A (scenes)**: 16 scenes per R24 §3.1. Two-population witnesses per
  R23 T2-B recipe (bright smooth low-leverage region + dim task-carrying
  fine structure; parameters documented). Blind generation: NO
  reconstruction may be run during generation; freeze+commit before arms.
- **T-B (library)**: L0 = frozen deployed SCAT32 + ridge power schedule;
  L1 scattered; L2 clustered; L3 compact k=32; L4 bar mixture;
  L5 radial/ring/contour; L6 multiscale blocks — ALL rebuilt at
  k_eff>=32 per R27 amendment 1; L7 = offline FW bank from
  the R18 support-expanding solver run on PROTOTYPE scenes = the six
  m1_dev (632900) instances ONLY (never the 16 bridge scenes — no
  leakage).
- **T-C (allocator; SUPERSEDED BY R25 — issue #17, docs/
  ROUND63_GPT_ROUND25_RULING_RAW.md is normative)**: the method is now
  **RLMI (Robust Library Mixture Illumination)**. Step 3 = the R25 §13
  eight-step pipeline: S=16 declared bootstrap scenario matrices →
  per-scenario bank-simplex A-risk oracles R*_s → standardized maximin
  relative-regret program (R25.1, parameter-free convex) → lexicographic
  closest-to-knob tie-break → exact 972-row materialization from the
  union design measure by constrained rounding/deterministic exchange
  (§9; budget/dose/peak/family constraints enforced; NO online
  full-dictionary OED) → post-materialization guard + realized-regret +
  KKT recomputation (R25.2 certificate, least-favorable scenario
  multipliers) → L0 fallback ONLY on solver/materialization/guard
  failure with explicit flag (MIXTURE_MATERIALIZATION_FAIL etc.) →
  acquire once + RQL. NO runtime thresholds (0.01/1.02/2%/hysteresis all
  retired). ε1, m, ε1²/2m remain mandatory DISCLOSURES (on-branch
  evaluation + BRANCH_VIOLATION flag per the falsifier addendum), never
  allocation decisions. Banks consumed as mixable measures ξ_k
  (T-B corrected schema: atom pool + normalized weights + per-row dose/
  cost/peak metadata); duplicate atoms merged with weight addition;
  final block order randomized with a frozen seed.
- **T-D (arms + harness)**: six arms (SCAT32-060, RIDGE-SCAT32, TRUE-X FW
  oracle, XHAT FW, RGLI, ORACLE-LIB), resources per R24 §3.2 (52+972,
  tau=50ns, nu in {200,2000}, c in {0,0.05}, 5 paired seeds; adaptive
  incident budgets capped by RIDGE's; pre-scan photons charged). Verify
  the c=0.05 jittered-hold forward path end-to-end on one smoke cell
  BEFORE the grid. Gates A, B numerically unchanged per R24 §3.4; Gate D
  renamed PLUGIN_LATENCY_PASS per R27 amendment 3 (claim label only);
  Gate C REPLACED by R25 §12: allocation informativeness via realized
  A_j = 1 − ŵ_j0 (C1: mean A_j ≥ 0.80 on rescue-needed scenes; C2: mean
  ŵ_j0 ≥ 0.75 on aligned controls; same 1 dB / 0.25 dB labels, no new
  thresholds; entropy/active-bank count reported, not gated). Four-gap
  decomposition per §3.5; gates are decision rules — no tuning, no
  re-runs, no scene swaps. If A+B pass but C fails: the method may be a
  static composite design but must NOT lead M2 as adaptive.
- **Compute**: local smoke; full grid sharded to Colab (pro2 x3 sessions
  first, pro1 added for speed) after the cost projection.

## Decision tree (frozen)

Gate A fail → method line + M2 dead; paper 2 ships as the jitter/knob
paper. Gates A–C pass → M2 preregistration (R24 §4) → RLMI leads paper 2;
PLUGIN_LATENCY_PASS then only classifies the claim (PASS: subsecond
plugin wording allowed; FAIL: batch/two-shot adaptive wording + full
latency distribution disclosed). No third option; no rescue redesigns
(R24 hard stop).
