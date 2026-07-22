# R34 — DLGI campaign adjudication at 6/7 promotion bars (fast ruling)

Reference: results/round63_next/DUAL_LEDGER_PROBE/ @ 95949da
(DUAL_LEDGER_PROBE_REPORT.md incl. the frozen bar-4 final section,
bar4_coverage.json, t3_reciprocity.json).

## The frozen probe outcome against YOUR R33 promotion bars

- Bar 1 identical-ledger comparisons: PASS (byte-identical records).
- Bar 2 model/gauge held-out residuals: PASS (0.89-0.97; gauge 3e-15).
- Bar 3 pilot-free precision <= 1.5x pilot: PASS everywhere; usually
  BETTER than pilot (t_c RMSE ratio <= 0.92, CV <= 1.15); at fast drift
  0.01x (pilots alias t_c <= 40).
- Bar 4 calibrated interval coverage: **FAIL after two honest repairs**
  (bias-corrected bootstrap then profile likelihood on logs; CV covered
  0.985-1.0 everywhere; t_c covered on 3/6 interior cells, the misses
  0.848/0.868/0.892 vs the 0.92 floor; mechanism = right-skew + bias at
  the slow-drift CRB floor and the high-CV cell).
- Bar 5 scene noninferiority: PASS (-0.04 to +1.76 dB; superior 7/9).
- Bar 6 (A,B,K)-predicted schedule behavior: PASS (reciprocity verified
  1e-15; paired best for BOTH ledgers as your theorem predicts).
- Bar 7 edge honesty: PASS (both t_c edges mapped as failure regions).

## The one decision

Per your all-bars rule DLGI is not auto-promoted. Rule on disposition:

(a) **Campaign GO with a frozen interval remedy**: preregister a
    SIMULATION-CALIBRATED interval procedure (per-cell coverage-corrected
    quantiles computed on the declared grid BEFORE the campaign, applied
    mechanically during it), with coverage re-verified as a campaign
    endpoint. Bar 4 becomes a campaign gate instead of a pre-gate.
(b) **Campaign GO with claim narrowing**: drop CI certification for t_c;
    claim RMSE-based precision (the bar-3 result) + calibrated CIs for
    CV only; intervals for t_c reported as descriptive.
(c) **No campaign**: 6/7 insufficient; DLGI archives as a strong
    theorem-backed method note.

If (a) or (b): freeze the confirmatory campaign skeleton in your R33 §6
style — fresh scene cohort, the (t_c, CV, photon-SNR) interior grid,
arms (pilot-free / pilot / oracle / plain-linear), primary endpoints and
kill bars, schedule set {paired, random, ordered} with (A,B,K)
predictions preregistered, and the equal-ledger discipline. Also state:
does DLGI at this evidence level lead the METHOD act of the operator's
one-flagship-manuscript (with the M1 ridge positives and the saturation
certificates as earlier acts), or is it a standalone shorter paper?
Blunt, fast — this is the last ruling of a long day.

Deliver as a GitHub issue on ccyyyYyzz/GI_a2 titled R34.
