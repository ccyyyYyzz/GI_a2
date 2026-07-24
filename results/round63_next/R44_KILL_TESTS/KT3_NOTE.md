# KT3 — Leak-Orthogonalized Efficient Sentinel (R44 §5.6, innovation #6, score 100)

**Verdict: PASS on the guaranteed R44 #6 statement — projecting the covariance score
off the leak/law span provably removes the nuisance covariance overlap (→0) while
retaining the bulk of target information; specificity is restored without an exact
wall.** `KT3_leak_orthogonalized_sentinel.py` + `.json`. M=32, 700 banks, T5
convention (z1=10mm, developed speckle l_c=50 σ=2π), z2∈{1,5}mm, ~23 min on one L4.

## Construction
Score objects live in covariance space with the Gaussian metric `<A,B>_V=tr(V⁻¹AV⁻¹B)`.
Because the signed flux is exactly linear in the scene, a global amplitude law is
**exactly** `dC_amp=((1+a)²−1)C0`; the nuisance span is `{C0(amp law), diag(shot),
dC_medium(l_c 50→45, σ 2π→1.8π)}`. The efficient signal is
`s_eff = dC_target − Proj_span{...} dC_target`.

## Measured (2 geometries × 2 changes)
| Quantity | z2=1(2%) | z2=1(5%) | z2=5(2%) | z2=5(5%) |
|---|---|---|---|---|
| retained info frac | 0.857 | 0.956 | 0.918 | 0.989 |
| clean amp overlap before→after | 0.36→**0** | 0.18→**0** | 0.25→**0** | 0.07→**0** |
| clean medium overlap before→after | 0.20→**0** | 0.11→**0** | 0.11→**0** | 0.03→**0** |
| operational AUC target change | −0.011 | −0.036 | −0.052 | −0.004 |

- **P_clean_nuisance_overlap_zeroed: PASS** — the projection drives the amplitude and
  medium covariance overlaps to exactly 0 while target info is retained 86–99%.
- **P_target_info_retained: PASS**.

## Estimator-layer honesty (per coordinator)
The **operational** AUC (oracle-W applied to shot-noisy blindly-estimated sample
covariances) is reported for reference but is **estimator-limited**: near contact
(z2=1mm, small λ) the oracle matched filter **anti-orders** — target AUC 0.24/0.17
(<0.5) — exactly the coordinator's "frozen plug-in W anti-orders the covariance
channel at small λ". Target-AUC change under projection is small (−0.01 to −0.05), so
projection does not itself cost operational power; the small-λ anti-ordering is the
independent estimator effect. The covariance-space specificity statement (overlaps→0)
is the robust, bench-agnostic result; oracle dC/W are twin-only. Nothing here touches
the Letter or sealed artifacts.
