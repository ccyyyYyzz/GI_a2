# IT3 — Leak-law stability (internal divergence stranger #8)

**Verdict: REPRODUCIBLE instrument constant — (L_pix, a, p) stable to ≤2.2% across 3
code seeds × M∈{8,32}. The leak law is quotable.** `IT3_leak_law_stability.py` + `.json`.
float64, ~34 s on one L4.

## Fits leak(z1) = L_pix + a·z1^p
| M | seed | L_pix | a | p |
|---|---|---|---|---|
| 8 | 10 | 6.99e-3 | 5.97e-3 | 0.904 |
| 8 | 21 | 6.71e-3 | 5.72e-3 | 0.904 |
| 8 | 33 | 6.97e-3 | 5.93e-3 | 0.902 |
| 32 | 10 | 7.15e-3 | 6.10e-3 | 0.904 |
| 32 | 21 | 6.78e-3 | 5.78e-3 | 0.904 |
| 32 | 33 | 7.05e-3 | 6.00e-3 | 0.903 |

- **CV(L_pix) = 2.2%, CV(a) = 2.2%, CV(p) = 0.1%**, all points within 2.2% — well inside
  the 15%/10% stability gate.

## Reading
The pixelation-floor + leak law is a **reproducible instrument constant**:
`leak(z1) = 6.94e-3 + 5.92e-3·z1^0.90` (mean over seeds/M). The exponent p ≈ 0.90 is
extremely stable (0.1% CV), confirming a near-linear (not quadratic) z1 growth. The law
is the carrier-cross `4·Re(C̄S̃)` factorization (IT9), **not** a 2k_p edge. The transfer
paragraph may use "reproducible instrument constant" wording rather than a bare bracket.
Nothing here touches the Letter or sealed artifacts.
