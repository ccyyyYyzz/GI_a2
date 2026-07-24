# IT4 — SVD certified-blind codes (internal divergence finer #3, THE SLEDGEHAMMER GATE)

**Verdict: PARTIAL (near-PASS) — a certified mean-blind code subspace exists with
near-perfect covariance retention; it misses the strict 64-dim ≤1e-4 gate by ~5×
but delivers a 219× operational mean-leak reduction.** `IT4_svd_certified_blind_codes.py`
+ `IT4_svd_certified_blind_codes.json`. float64, ~4 min on one L4.

## Method
The mean-channel signed intensity is **exactly linear** in the code
(`D=Re[U_c·conj(U_s)]`), so the beyond-band mean leak is a linear operator L_be(z1).
Built L_be(z1) and the in-band response L_in(z1) for z1∈{0,2,5,10} mm (n_in=120
in-band DOF, n_be=240 beyond-band annulus), stacked over z1, and solved the
generalized eigenproblem `min v^T(L_beᵀL_be)v / v^T(L_inᵀL_in)v` — the bottom
eigenvectors are code directions jointly mean-blind at every z1.

## Gate (coordinator) vs measured
| Quantity | Gate | Measured | |
|---|---|---|---|
| bottom-64 joint leak (max) | ≤1e-4 PASS / ≥1e-3 KILL | **4.99e-4** (median **8.59e-5**) | between gates → PARTIAL |
| covariance λ (blind vs random) | within 3× | **λ_blind 2.64e-3, λ_rand 3.17e-3, ratio 1.20** | **retained** |
| full-spectrum rel leak | — | min **1.4e-6**, max 0.113 | huge blind subspace exists |
| operational mean leak / code | — | **blind 2.6e-4 vs random 5.6e-2 (219× lower)** | strong positive |

## Reading
The sledgehammer is **real**: a code subspace exists that is jointly mean-blind at
all four propagation distances and retains ~83% of the covariance detection power
(λ ratio 1.20 ≪ 3×), cutting the per-code mean leak **219×** (5.6e-2 → 2.6e-4). The
bottom of the spectrum reaches 1.4e-6 relative leak. The only miss is the strict
64-dimensional ≤1e-4 threshold: the 64th direction sits at 5e-4 (median 8.6e-5), so
a ~48-dimensional blind subspace would clear 1e-4 cleanly. This is **not a KILL**
(far from the 1e-3 code-independence threshold) — certified-blind codes are a viable
next-paper core and a bench calibration step (route 2, code-space wall restoration),
pending a subspace-dimension/leak trade sweep. Nothing here touches the Letter or
sealed artifacts.
