# KT2 — Physical Jet Transmutation by Channel Leakage (R44 §5.2, innovation #2, score 125)

**Verdict: R44 innovation #2 (Jet Transmutation) CONFIRMED across all four (z1,z2)
cells.** `KT2_jet_transmutation.py` + `KT2_jet_transmutation.json`. M=16, 350 banks,
developed speckle (l_c=50, σ=2π), ~6 min on one L4.

## Exact structure exploited
The noiseless signed flux `f=(Gp−Gm)·x` is **linear** in the scene, so the bucket
covariance is **exactly quadratic** in x and `dC(ε)=ε·B_δ+ε²·C_δ` with **no higher
terms**. The first-order-orthogonal (curvature-rescued) direction is the minimum
right-singular vector of the whitened first-order response operator over the
beyond-band annulus; its orthogonality residual s_min/s_max ≈ **1–5×10⁻¹⁸** (a
genuine, exact first-order-null direction exists — as the coherence-rank bound
predicts). The mean-leak coefficient is `A=ℓ₁ᵀV⁻¹ℓ₁` with `ℓ₁=J_mean·δ` (the T1 wall
leak); the curvature coefficient is `B=½tr[(V⁻¹C_δ)²]`.

## Frozen prediction vs measured (4 cells: (10,1),(10,5),(5,1),(5,5) mm)

| Quantity | Predicted | Measured (4 cells) | Verdict |
|---|---|---|---|
| projected (covariance-only) slope | → 4 (m=2 restored) | **4.00, 4.00, 4.00, 4.00** | **PASS** |
| unprojected total slope, low-ε | → 2 (m=1, leak-dominated) | **2.14, 2.20, 2.12, 2.18** | **PASS** |
| unprojected total slope, high-ε (>ε_cross) | → 4 (crossover) | 3.06, 3.69, 3.29, 3.65 | crossover captured |
| generic direction cov slope (control) | ~2, sub-quartic (no rescue) | **1.58, 2.28, 1.58, 2.29** | **PASS** |

The full physical law `C_phys(ε)=A ε²+B ε⁴` (R44 eq 2.12) and the effective order
`m_eff(ε)=(A+2Bε²)/(A+Bε²)` crossing 1→2 near `ε_cross=√(A/B)` are reproduced
exactly: unprojected slope rises from 2 (below ε_cross) toward 4 (above); projecting
out the measured mean-leak tangent removes the `Aε²` term and restores the pure `ε⁴`
(slope 4) curvature class at every ε.

## Bearing on the 125-score theorem
A tiny first-order hardware mean leak (present here as the beyond-band T1 leak
`J_mean·δ≠0`) changes the strict local contact order of a curvature-rescued covariance
direction from m=2 to m=1 — and an efficient projection of the calibrated mean-leak
tangent restores the m=2 class. This is the theory statement of the statistical
wall-restoration route (its operational realization is KT3). Oracle contact-order
quantities are the legitimate object of a jet-order theorem (unlike KT6 detection,
which carries the blind-estimator caveat). Nothing here touches the Letter or sealed
artifacts.
