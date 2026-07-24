# QUARANTINE: post-freeze JET_TEST re-run (NOT used anywhere)

The JET_TEST agent's buffered process finalized a second run AFTER the
coordinator's freeze commit 1bf29f1. Per the R43 green-light protocol
(§7.2: no post-freeze regeneration), the FROZEN artifacts at 1bf29f1 are
authoritative and are the sole source for the manuscript (see
paper_prl/CLAIM_SOURCE_MATRIX.md). This directory preserves the re-run
for provenance only. Both value sets support identical qualitative
conclusions (e.g. KL orders 2.015 vs frozen 2.038 / 4.000 both -> order 2/4;
MC 1.17/1.99 vs 0.95/2.05 both -> m=1/2; CUSUM -2.19/-3.98 vs -2.16/-3.92).
