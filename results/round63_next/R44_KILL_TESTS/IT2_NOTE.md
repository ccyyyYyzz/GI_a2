# IT2 — Floor decomposition + leak-law fit (internal divergence nextpaper #9)

**Verdict: DECOMPOSES — the z1=0 pixelation floor is a real DMD instrument constant
(macro-pixel staircase dominated), quotable with a leak law, not just a bracket.**
`IT2_floor_decomposition.py` + `.json`. float64, M=12, ~6 s on one L4.

## Floor decomposition (z1=0)
| variant | leak | Δ vs baseline |
|---|---|---|
| baseline (fill=0.92, mp=8) | **7.15e-3** | — |
| fill=1.0 (no mirror gap) | 7.55e-3 | +6% |
| mp=16 (2× macro staircase) | 3.56e-3 | **−50%** |

The floor **moves 50% when the macro-pixel staircase is coarsened** (and 6% with
fill-factor) → it **decomposes** into physical DMD terms. `L_pix = 7.15e-3` is a real
instrument constant dominated by the macro-pixel sample-and-hold staircase, with a
minor fill-factor contribution — **quotable**, not a twin grid artifact.

## Leak law (extended z1 grid)
leak(z1) = 7.15e-3 → 9.38e-2 over z1 = 0 → 20 mm, fit
**leak(z1) = L_pix + 6.10e-3 · z1^0.91** (nearly linear in z1, p≈0.91), consistent with
the sealed T1 bracket [6.8e-3, 8.9e-2]. The near-linear (not quadratic) growth confirms
the R44 caveat that no linear-vs-quadratic small-z1 claim should be made beyond this fit.

## Caveat
The relay-pupil sub-variant in the script probes shells 6–9, which lie *within* the
2k_R passband for kR∈{5,6,8}, so it is not the right probe for the pupil benefit (KT1
already established the pupil walls energy *beyond* 2k_R, not within) — those pupil
numbers are not load-bearing here. The decomposition result (macro-staircase 50%, fill
6%) is the finding.

## Reading
The transfer paragraph can quote `L_pix ≈ 7.2e-3` as a DMD instrument constant plus the
near-linear leak law leak(z1) = L_pix + 6.1e-3·z1^0.91, rather than only the bracket.
Nothing here touches the Letter or sealed artifacts.
