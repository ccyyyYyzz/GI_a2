# Letter 13 — build the R24 DEV bridge (the method-kill experiment)

From: main agent. Date: 2026-07-22. Authority:
`docs/ROUND63_GPT_ROUND24_RULING_RAW.md` (issue #16) — read it in FULL
first; §2 (RGLI architecture), §3 (the DEV bridge), and §1.2 (failure
modes F1–F10) are your build spec. Context: the operator redirected the
program — the deliverable must be a genuinely effective, deployable
method; GPT rejected our online-FW RGAI as the deployed lead and specified
RGLI (two-shot, estimator-aware, finite-library router with abstention).
Before ANY campaign, the bridge experiment decides everything.

This is the most consequential build since the M1 freeze, and it is yours.

## What to build (in dependency order)

1. **16 fresh DEV-only scenes** exactly per ruling §3.1 (4 contour-weak,
   4 analytic two-population witnesses — construct these from the R23 T2-B
   recipe: bright/low-leverage region + dim/task-carrying fine structure —
   4 low-amplitude microtexture/mixed, 4 aligned controls). New seed range
   disjoint from every existing cohort (document it). No
   rejection/replacement after reconstruction — generate blind, freeze the
   list before any arm runs.
2. **The library** per §2.2: 8 banks L0–L7, each materialized as exact
   972-row nonnegative pattern sets satisfying ALL safety guards (budget,
   peak, dose band, admission); L0 = deployed SCAT32 + ridge knob (reuse
   frozen assets); L7 = offline FW bank built with the R18
   support-expanding solver on a mixed prototype set (NOT on the 16 bridge
   scenes — no leakage). Recompute guards after materialization (F5).
3. **The RGLI router** per §2.3 steps 0–5, with the prototype constants AS
   GIVEN (U_1/r ≤ 0.01; 1.02 knob-preference; 2% library margin; S=16
   frozen-seed bootstrap draws; ABSTAIN semantics; OUT_OF_LIBRARY flag).
   The reconstruction-aware score R_k = tr[W(H0+F_k)^{-1}] with
   H0 = F_pre(xhat) + lambda_TV * R_delta(xhat) + eps*I, W=I on the
   declared subspace. Do not tune any constant against outcomes — they are
   architecture, not fit parameters.
4. **The six arms** per §3.3: SCAT32-060, RIDGE-SCAT32, TRUE-X FW oracle
   (full safety-constrained geometry optimization ON THE TRUE SCENE —
   nondeployable by design), XHAT FW (the rejected RGAI, as diagnostic),
   RGLI, ORACLE-LIB (best bank by true PSNR). Resources per §3.2:
   M=52+972, tau=50ns, nu in {200,2000}, c in {0,0.05}, five paired seeds,
   adaptive arms' incident budget capped by the RIDGE arm's, every
   pre-scan photon charged. Hard decision cell = (c=0.05, nu=2000).
   NOTE the c=0.05 requirement: the forward simulator must run with
   jittered holds — the campaign machinery has the hold-law path from the
   jitter studies; verify it end-to-end on one smoke cell first.
5. **Gates A–D verbatim** (§3.4) + the four-gap decomposition (§3.5:
   surrogate-to-image / localization / bank-approximation / routing —
   never averaged). Latency measurement per Gate D on this machine
   (declared reference CPU), excluding final reconstruction.

## Hard rules

- DEV only; the 16 scenes are DEV by construction; no confirmatory
  (633000+) object is touched.
- Gates are DECISION RULES, not targets: no constant tuning, no scene
  swapping, no re-runs after seeing outcomes. If Gate A fails, the
  deliverable is the failure evidence — that outcome is as valuable as a
  pass (it kills the method line cleanly, per the ruling's hard stop).
- Compute: local CPU for harness bring-up and smoke; the full grid
  (16 scenes x 6 arms x 2 nu x 2 c x 5 seeds + oracle FW solves) will need
  Colab — plan shards the way M1 did (the colab/ playbooks are yours),
  but DO NOT launch cloud compute before reporting the harness-validated
  smoke + cost projection in letter 14. The FW oracle solves are the
  expensive item — budget them explicitly.
- Commit as yourself in clearly-scoped commits; the 16-scene list and all
  bank SHAs frozen and committed BEFORE the arms run.
- Anything ambiguous in the ruling: quote + interpret in writing (your
  I-register pattern), never silently choose.

Deliverable sequence: letter 14 = harness ready + smoke + frozen
scene/bank manifest + cost projection; then (after my go) the full bridge
run; letter 15 = gates verdict + decomposition + your recommendation.
The falsifier-regression and R16-diagnostic agents are running in
parallel on other files; no collisions expected (results/round63_next/
FALSIFIER_* and R16_DIAGNOSTIC_RUN are theirs; take
results/round63_bridge/ as your output root).