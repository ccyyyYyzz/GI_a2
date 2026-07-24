# COHERENCE_CROSSOVER — R45 named gate for the memory-effect line (R45 §B5)

**This is the R45-NAMED GATE for the memory-effect line (the gate, not the campaign).**

**VERDICT: PASS (primary bar).** In the decisive reduced-R geometry the code-count crossover
independently recovers the field coherence rank: crossover-fitted `sqrt(q_fit) = 17.99` vs the
independent KT4 coherence rank `R = 18.9` — **4.8% disagreement**, well inside the frozen 25%
bar. The memory-effect line's named gate opens.

Engine: WAVE_TWIN `wave_twin.py` + `twin_pool.py`, shared bank pool. r_eff(M) = participation
ratio of the whitened covariance-Fisher Jacobian (KT4 method); independent R = participation
ratio of the object-plane field mutual-coherence tensor (KT4 method). Run: `coherence_crossover.py`
on **Colab L4** (complex64, WAVE_DTYPE=32), 861 s.

## Frozen bar (declared before running)

Primary: `|sqrt(q_fit) − R_KT4| / R_KT4 < 0.25`. Corroborating (R45 §B5): pre-crossover slope
within 0.2 of 2; post-crossover slope < 0.8 where p/q > 3.

## Two configurations

**(A) Named-gate reference — z2=5 mm, l_c=50 µm, developed (σ_φ=2π).** KT4 `R = 123` (matches
KT4's ~110–122). Since `R² ≈ 15000` exceeds the resolvable rank (bank/scene budget), the sweep
saturates at the structural cap (`r_eff → ~1189`, not R²) and the crossover is **not reachable** —
exactly KT4's finding. Reported as the reference; not the decisive test.

**(B) Reduced-R decisive geometry — z2=2 mm, l_c=250 µm, σ_φ=0.6, full illumination**,
auto-calibrated to `R_KT4 = 18.9` (target 20–40), n_bank=800 (> R²≈357 so the R² saturation is
resolvable), scene basis ns=1681.

| M | p=M(M+1)/2 | r_eff (measured) | regime |
|---|---|---|---|
| 4 | 10 | 9.3 | pre-crossover |
| 8 | 36 | 30.7 | pre-crossover |
| 12 | 78 | 62.6 | rising (~p) |
| 16 | 136 | 100.1 | near M*≈25 |
| 24 | 300 | 172.0 | crossover |
| 32 | 528 | 219.5 | saturating |
| 48 | 1176 | 257.2 | saturating |
| 64 | 2080 | 277.3 | saturated → q |
| 96 | 4656 | 289.0 | saturated → q |

Fit `r_eff = pq/(p+q+1)` → `q_fit = 323.5`, `sqrt(q_fit) = 17.99`, `M* = √2·R = 25.4`.

## Verdict detail

- **Primary bar PASS:** `sqrt(q_fit)=17.99` vs `R_KT4=18.9` → rel err **0.048 < 0.25**. Two
  independent estimators of the coherence rank — the code-count crossover and the field
  mutual-coherence participation ratio — agree within 5%, with no transmission-matrix access.
- **Corroborating post-slope PASS:** post-crossover slope 0.182 < 0.8 (clean R² saturation).
- **Corroborating pre-slope soft (reported honestly):** pre-crossover slope **1.61**, below the
  [1.8, 2.2] band. Cause: even the smallest M in the sweep already has `p` not `≪ q` (at M=8,
  p/q≈0.11), so the finite-q term pulls the local slope below the ideal 2 before a clean
  unsaturated decade can develop. A cleaner slope-2 leg would need smaller M or larger R; it does
  not affect the primary agreement result, which is decisive.

## Consequence

The named COHERENCE-CROSSOVER gate **opens**: coherence-rank spectroscopy recovers the scattering
coherence rank from detector-code scaling. Per R45, this reopens the memory-effect line as a
genuine future thickness/stacking campaign (the algebraic `rank ≤ min(N,R²)` bound was never in
question; this validates the *measurable* spectroscopy claim). Bench-phase only; nothing here
touches the Letter or sealed artifacts.
