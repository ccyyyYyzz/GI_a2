# IT4b — certified-blind-code dimension/leak trade (coordinator follow-up to IT4)

**Verdict: the IT4 "64-dim at 5e-4" reading is BOTH geometry-fundamental AND
gate-statistic-conservative — and the operationally-relevant single-z1 case is far
more forgiving.** `IT4b_spectrum_trade.py` + `.json`. Deterministic float64,
ridge-free, 16 s on one L4.

## 1. Full joint spectrum (z1∈{0,2,5,10}mm stacked)
Relative-leak generalized eigenvalue crossings: **d@1e-5 = 8, d@1e-4 = 35,
d@1e-3 = 80**. So the strict joint-blind subspace is ~35 dims (not 64) — this part
is geometry-fundamental (the joint-over-4-z1 constraint is what shrinks it).

## 2. Gate-statistic variants (the max-statistic is conservative)
| d | max leak | median | RMS of random code in subspace |
|---|---|---|---|
| 32 | **8.6e-5** | 1.8e-5 | 4.4e-5 |
| 48 | 1.43e-4 | 5.5e-5 | 7.2e-5 |
| 64 | 4.99e-4 | 8.6e-5 | **1.68e-4** |
| 96 | 3.5e-3 | 2.0e-4 | 1.3e-3 |

A **32-dim** subspace is cleanly blind on the worst-case max (8.6e-5 ≤ 1e-4). At
d=64 the operationally-relevant RMS leak of random codes drawn from the subspace is
**1.68e-4 — 3× below the 4.99e-4 worst-direction max** that IT4 gated on. So codes
actually drawn from the d=64 subspace leak ~1.7e-4, not 5e-4.

## 3. Per-z1 vs joint — the decisive operational point
**Single-z1 (z1=0, relay-conjugate bench): d@1e-4 = 80 (and d@1e-5 = 80)** vs
**joint d@1e-4 = 35 (+45 gain)**. A bench operating at a fixed conjugate plane has a
**2.3× larger** certified-blind subspace, essentially perfectly blind (≤1e-5) in 80
dims. The joint-over-all-z1 requirement is what costs the dimension; a real bench does
not need it.

## 4. Window-apodized route (naive window is a KILL — confirms coordinator)
V7-style Tukey-apodized ROI makes the spectrum **WORSE**: joint d@1e-4 drops 35→24,
d@1e-5 drops 8→3, and the d=64 leak jumps 4.99e-4 → **4.66e-2 (93× worse)**. The
naive separable window's rim effect *hurts* — routes do **not** multiply for naive
windows (DPSS/Slepian remains the only live window route, = IT5).

## 5. Sanity
Operator build is deterministic float64; bottom-5 eigenvalues (1.43e-6, 1.43e-6,
1.91e-6, 4.23e-6, 7.48e-6) and d@1e-4=35 are **identical** for B-ridge 0, 1e-8, 1e-6.

## Reading
The certified-blind-code sledgehammer (IT4) is **more viable than the d=64/5e-4
reading suggested**: (i) a 32-dim subspace is cleanly ≤1e-4; (ii) operational RMS leak
is 3× below the gated max; (iii) the operationally-relevant single-z1 conjugate bench
gives 80 blind dims. The only failed route is naive window apodization (rim-effect
kill). Nothing here touches the Letter or sealed artifacts.
