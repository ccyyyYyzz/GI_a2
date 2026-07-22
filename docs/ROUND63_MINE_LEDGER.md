# ROUND63 mine ledger — anomaly excavation tracker (2026-07-20)

Status of every item from the full-journey reflection (出乎意料/奇怪/忽视/避开),
post-R23. Papers remain two; digs are post-hoc/descriptive unless promoted.

## ACTIVE DIGS

1. **Falsifier regression** — COMPLETE, verdict NEGATIVE for the raw
   bridge (2026-07-20/22): P1 FAIL (Align1 ~0.02 uniform, no between-scene
   signal), P2 nominally passes but is contrast-in-disguise
   (corr(M1sq, Cu)=0.999; partials collapse), P3 FAIL (contour mid-pack).
   ROOT CAUSE IS SCOPE, not concept: the production ridge sits AT the
   finite-nu kernel argmax, ~47% of loads beyond it, J' sign-splits, and
   alpha* is scene-constant (-7.1e5) — the R23 declared monotone-concave
   branch is violated at the deployed operating point, exactly the
   fragility R23 T4 predicted. CONSEQUENCE (pre-agreed gate): the R23
   alignment page does NOT enter paper 2; the material stays next-line.
   The failure directly supports R24's re-architecture (estimator-aware
   robust router + abstention instead of raw-residual gating), and yields
   a build requirement: diagnostics must be evaluated on-branch (at
   candidate-bank operating points), never at the deployed ridge alone.
   → results/round63_next/FALSIFIER_REGRESSION.md.
2. **R16 2M-frame score-identity diagnostic** (the −8% endpoint bias;
   now load-bearing for the elasticity map too, since alpha*_J = −ρJ''/J').
   → Colab pro2 r63diag_a (Opus, running); runbook frozen protocol.
3. **Speed-ratio decomposition** (the c=0 signature): S_gate vs
   detected-rate ratio (~20.1); per-scene deviations; clip-cost accounting
   J(11)/J(22)≈0.96 at low-nu crossings; the derived statement that at c=0
   long-window, information per detected photon is load-invariant (=1) and
   the FINITE-window ridge is exactly where counting stops equaling
   learning. → results/round63_next/SPEED_DECOMPOSITION.md (Opus, running).
4. **Quantile-calibration rescue probe** — COMPLETE, verdict NEGATIVE
   (2026-07-20): harness reproduced frozen dQ to 0.00000 dB; Q90/Q75
   calibration pulls the past-ridge load fraction from 46–48% to 7–16%
   with essentially no quality response; no failing scene rescued; even
   the per-pattern MEAN-CLIPPED genie fails (worsens contour_0/1 by
   0.65–1.61 dB). The dispersion-past-ridge hypothesis is REFUTED for the
   confirmatory negatives (which are contour_0/1/2 + spokes_0/1;
   microtexture is positive — earlier brief premise corrected). Deficit
   lives in geometry/task alignment, NOT load allocation: no reweighting
   of light (rank-one OR per-pattern) helps these scenes; only different
   PATTERNS could. Sharpens the falsifier's role: if the alignment
   functional predicts these scenes, the closed story is "load control
   cannot help misaligned scenes (probe) + geometry headroom exists
   (certificate) + the functional says which is which (falsifier)".
   → results/round63_next/QUANTILE_RESCUE_PROBE.md.

## QUEUED (batch into the next GPT round)

5. **Source-RIN into the law**: does source intensity noise enter as a
   c²-like cycle-variance term or differently? (R22's C1 refutation warns
   against casual factor-adding — needs a real derivation.) Shortest
   bridge from the law to real hardware.
6. **Photon-information exchange rate box**: formal one-paragraph statement
   of dig 3's identity for paper-2 discussion (if R24 approves wording).
7. **Arrays preemption paragraph**: J_m = m·J(ρ/m,c) + why single-detector
   operation still matters (cost/crosstalk/shared-source; the rank
   hierarchy reappears as common-mode vs per-element gains).
8. **Telemetry removal of the jitter cap** (2026-07-20, from the
   core-limitation analysis): E1 localizes the loss as E[Var(L|N)];
   observing the hidden state (TDC event timestamps + re-arm/quench
   telemetry) drives the conditional variance to ~0, restoring monotone
   J → ρ/(1+ρ) and erasing the jitter cap up to TDC precision. Prior-art
   fence to check: nuclear live-time clocks (Gedcke–Hale) measure L for
   rate CORRECTION; the Fisher-restoration/optimal-load claim appears
   unwritten. Natural bench proposal (SPAD + TDC + quench monitor:
   "the ridge flattened by telemetry"). Needs: prior-art sweep + one
   short GPT derivation round before any claim.

## PARKED (mine, next work)

8. **181 NUMERICAL_UNRESOLVED recondition**: R23 predicts load-coordinate
   LP fixes the uniform ~65s failures; would sharpen 0/299/181 to a
   post-hoc sensitivity. Next-paper material; frozen distribution stands.
9. **Estimator-layer science** (lambda_TV × load interaction): image-domain
   non-corollary fence (R23 T5.3); needs its own study.
10. **Sub-Poissonian illumination × dead time** (R22 reach b): quantum
    interface; archived.
11. **T3-C nonlinear cluster continuation**: conjecture with named obstacle.

## CLOSED BY UNDERSTANDING (no further dig)

- 19.13≈20.1 near-tautology → reframed as the c=0 channel signature
  (dig 3 formalizes).
- Contour mechanism duality (dispersion-past-ridge vs misalignment) →
  SETTLED EMPIRICALLY by the rescue probe (2026-07-20): the dispersion
  mechanism is refuted; misalignment carries the weight. The earlier
  "ridge-top flatness unifies both" synthesis was half wrong — the
  flatness does make clipping cheap (speed decomposition confirmed
  3–12% clip cost), but the contour deficit is not a load effect at all.
- 两边对比度 → the three failure modes span the alignment functional:
  amplitude (C_u→0, microtexture), elasticity mismatch ((α−α*)², contour),
  brightness-invisible geometry (σ_ζ, R22 counterexample).
- Absolute-quality honesty, microtexture relative-success framing,
  GPT-connector phantom commits, PSNR_rad vs PSNR: final-pass wording
  notes, not science digs.
