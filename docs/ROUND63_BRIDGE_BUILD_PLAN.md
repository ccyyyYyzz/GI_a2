# DEV bridge build plan (Fable-owned; R24 issue #16 is the normative spec)

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
  L1 scattered k=16; L2 clustered k=16; L3 compact k=32; L4 bar mixture;
  L5 radial/ring/contour; L6 multiscale blocks; L7 = offline FW bank from
  the R18 support-expanding solver run on PROTOTYPE scenes = the six
  m1_dev (632900) instances ONLY (never the 16 bridge scenes — no
  leakage).
- **T-C (router)**: RGLI per R24 §2.3 steps 0–5 with prototype constants
  AS GIVEN (U1/r<=0.01; 1.02 knob preference; 2% margin; S=16 frozen
  bootstrap draws; ABSTAIN + OUT_OF_LIBRARY_OR_UNCERTAIN). Plus the
  falsifier-derived requirements: all diagnostics evaluated ON-BRANCH at
  each candidate bank's own loads; BRANCH_VIOLATION disclosure when the
  R23 bound hypothesis fails at a candidate's operating point (bound
  reported as UNAVAILABLE, never vacuous).
- **T-D (arms + harness)**: six arms (SCAT32-060, RIDGE-SCAT32, TRUE-X FW
  oracle, XHAT FW, RGLI, ORACLE-LIB), resources per R24 §3.2 (52+972,
  tau=50ns, nu in {200,2000}, c in {0,0.05}, 5 paired seeds; adaptive
  incident budgets capped by RIDGE's; pre-scan photons charged). Verify
  the c=0.05 jittered-hold forward path end-to-end on one smoke cell
  BEFORE the grid. Gates A–D + four-gap decomposition per §3.4–3.5;
  gates are decision rules — no tuning, no re-runs, no scene swaps.
- **Compute**: local smoke; full grid sharded to Colab (pro2 x3 sessions
  first, pro1 added for speed) after the cost projection.

## Decision tree (frozen)

Gate A fail → method line + M2 dead; paper 2 ships as the jitter/knob
paper. Gates A–D pass → M2 preregistration (R24 §4) → RGLI leads paper 2.
No third option; no rescue redesigns (R24 hard stop).
