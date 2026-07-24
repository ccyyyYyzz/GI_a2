# KT4 — Coherence-separation rank sketch (R44 §5.3/§4.2, innovation #3, score 125)

**Verdict: PARTIAL — the rank ≤ min(N,R²) bound HOLDS at every z2 and the
Jet-Order Separation is clean (orthogonal quotient-jet slope stays 4.0), but the bound
is non-binding in this regime and R does NOT collapse toward 1.** M=16, 250 banks,
z1=10mm, developed speckle, ~8 min on one L4. `KT4_coherence_rank_sketch.py` + `.json`.

## Measured (z2 ∈ {0.1,0.5,1,2,5,10,20} mm)
| z2 (mm) | rank(I_cov) (PR) | count | R (coherence PR) | R² | min(N,R²) | bound |
|---|---|---|---|---|---|---|
| 0.1 | 99.8 | 136 | 107.4 | 11528 | 361 | ✓ |
| 1.0 | 101.5 | 136 | 108.2 | 11704 | 361 | ✓ |
| 5.0 | 110.2 | 136 | 114.4 | 13081 | 361 | ✓ |
| 10.0 | 111.5 | 136 | 120.3 | 14481 | 361 | ✓ |
| 20.0 | 102.0 | 136 | 122.0 | 14877 | 361 | ✓ |

- **Rank bound: PASS** — `rank(I_cov) ≤ min(N_scene, R²)` holds at every z2 (rank ~100–112
  ≪ min(361, R²≈12000)). Never violated.
- **Jet-Order Separation: PASS** — the first-order-orthogonal (curvature) quotient-jet
  slope is **4.0 at every tested z2** (generic 1.6–2.3): the m=2 curvature class
  survives regardless of z2/rank, exactly R44's "quotient jets remain finite even while
  reconstruction rank collapses."

## Honest caveats (the coordinator's contested-bookkeeping flag, confirmed)
1. **The R² bound is non-binding here.** The covariance-Fisher rank is capped by the
   probe's code-pair count `M(M+1)/2 = 136` (count = 136 at every z2), **not** by R²
   (≈12000) nor N_scene (361). So the `≤R²` inequality holds only vacuously in this
   regime — the code count is the real limit.
2. **R does not collapse toward 1.** The coherence rank R stays ~107–122 (a
   thin developed phase screen remains highly multimode); it even rises slightly with
   z2. The R44 "complete-scrambling R=1" endpoint is **not** reached on this
   single-thin-screen z2 continuum — it would require memory exhaustion (thicker/stacked
   scattering). R here is a field-coherence participation-ratio proxy and is
   n_bank-sensitive (M=8/120-bank smoke gave R~73–79), so absolute R values are
   estimator-dependent; the qualitative "R large, non-collapsing" conclusion is robust.

## Reading
The two structural claims of the Coherence-Rank / Jet-Order Separation theorem that
*are* testable on this geometry hold cleanly (bound never violated; curvature jet-order
survives). But this twin geometry cannot exhibit the R→1 collapse or make R² the
binding constraint, so it does not by itself certify the memory-effect rank narrative —
consistent with the internal divergence's call to reopen that line only with corrected
bookkeeping and a genuine scrambling endpoint. Nothing here touches the Letter or
sealed artifacts.
