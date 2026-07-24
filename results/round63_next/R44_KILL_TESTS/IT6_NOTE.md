# IT6 — Segmented-bucket P² law (internal divergence stranger #7)

**Verdict: KILL by the mandatory leak audit — the P² covariance gain is real and large,
but segmenting the bucket breaks the mean-channel wall (per-segment mean leak floods the
witness). Single-bucket bank counts stand; the bench stays simple.**
`IT6_segmented_bucket_p2.py` + `.json`. M=16, 700 banks (P·M=256 ≪ 700, well sampled),
~7 min on one L4.

## Method
Re-reduce the SAME twin speckle kernels (no re-propagation) into P∈{1,4,16} spatial
segments; covariance detection noncentrality λ over the P·M buckets; scaling exponent of
λ vs P; **mandatory** per-segment beyond-band mean-leak audit (mean d′ / cov d′).

## Measured
| cell | change | exponent (λ∝Pᵖ) | P=4 buys | mean/cov d′ @ P=1 → P=4 → P=16 |
|---|---|---|---|---|
| z2=5mm | 2% | 1.68 | 11.8× | 0.80 → 1.81 → 1.43 |
| z2=5mm | 5% | 1.74 | 13.7× | 0.88 → 1.84 → 1.44 |
| z2=1mm | 2% | 1.96 | 15.3× | 3.25 → 11.6 → 7.53 |
| z2=1mm | 5% | 2.02 | 17.3× | 3.56 → 11.9 → 7.61 |

- **The P² law is real** — exponents 1.68–2.02 (≈2, super-linear), P=4 buys **12–17×**
  covariance noncentrality: a quadrant/multi-cell diode would sharply raise per-bank
  covariance information. Both the exponent (≥1.5) and P=4-gain (≥6×) criteria PASS.
- **But the mandatory mean-leak audit KILLS it** — at P=1 (single bucket) the mean leak
  is below the covariance signal (ratio 0.8 at z2=5mm: the mean-channel wall holds), but
  **segmenting to P≥4 breaks it**: the per-segment mean d′ exceeds the covariance d′
  (ratio 1.8 at z2=5mm, up to 11.9 near contact). The whole-bucket complementary
  cancellation that suppresses the mean leak is incomplete per-segment, so a segmented
  detector reintroduces exactly the mean leakage the two-ledger design cancels.

## Reading
The segmented-bucket P² gain is genuine but **cannot be used as-is**: it comes packaged
with a mean-leak flood (per-segment d′ 1.4–11.9× the covariance signal) that would
contaminate attribution. This is the coordinator's explicit mandatory-audit KILL branch
("segment windowing floods the witness annulus"). A segmented detector would need to be
paired with the per-segment mean-blind code subspace (IT4) or a concentration guard (IT5)
to be viable — a next-paper engineering direction, not a drop-in bench upgrade.
Single-bucket bank counts stand. Nothing here touches the Letter or sealed artifacts.
