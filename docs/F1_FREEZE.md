# F1 IMMUTABLE FREEZE — ROUND63 dead-time-aware high-flux SPI campaign

**Date**: 2026-07-18. **This commit is the F1 freeze** (spec D2.3 §1; GPT
round-7 §7 lock conditions 1–9 all satisfied). Ledger:
results/round63/F1_LEDGER.json (payload self-hash inside; every frozen input
surface enumerated with sha256).

Frozen by this commit (verbatim clause, D2.3 §1): the primary hypothesis and
pass rule; primary ρ₀=0.05, ρ₁=0.6; image and seed identifiers (DETAIL-24
generator code + params manifest + seeds 631000..631023 + image SHA256;
NATURAL-24 = STL test 0..23; 5 confirmatory seeds 0..4); all grids
(M=n=4096, ν∈{5,10,20,50,100,200,500,1000,2000}, ρ̄∈{0.05,0.3,0.6,1,2} +
S2B/S2C/S3 per D2.3 §3); reconstruction arms; the analytic_score_concentration
λ_TV rule with C₀=∞ ⇒ global c=0.50 (C0_FROZEN.json, Pass-A commit 9f1740f);
the Q90 quality target + SAFE_RANGE_UNINFORMATIVE rule + censoring taxonomy +
10,000-draw nested family-stratified bootstrap; the fixed-budget quality-gain
secondary endpoint; figure subjects and crop coordinates
(docs/ROUND63_FIGURE_PREREG.md); code + environment (env_lock.txt) + input
SHA256; all shard manifests (82 default + 5 blocked); and the complete
expected-cell table (45,396 rows).

Lock trail: rule iterations rounds 3→7 with outcome-blind probes at every
step (docs/ROUND63_GPT_ROUND{4,5,6,7}_*.md; results/round63_gof_probe/;
results/round63_s1_passA/, _passB/); Pass-A C₀ freeze 9f1740f; lock-8
implementation smoke on the six dev generator seeds (structural PASS,
results/round63_detail_devsmoke/).

**After this commit** (round-7 §7, verbatim): no further endpoint,
image-class, regularization, geometry or grid revision is permitted. Changes
are limited to demonstrable implementation defects established by
outcome-blind unit tests, each requiring a version bump and rerunning all
affected cells. **If DETAIL-24 fails, that is the campaign's confirmatory
result.**
