# KT5 — Étendue / N_eff collapse (R44 §5.4, innovation #4, score 100)

**Verdict: N_eff grows with illumination étendue (PASS) and controls detection latency
(T_det ∝ N_eff, ~25% scatter); the precise universal collapse vs n/N_eff is PARTIAL.**
`KT5_etendue_neff_collapse.py` + `.json`. M=16, 300 banks, z1=10mm, z2=5mm, ~6 min/L4.

## Method
Illumination aperture diameter swept over 5 diameters (0.40–1.00 of the DMD active
region, a circular stop on the illumination). N_eff read directly from the intensity
bucket contrast `N_eff = E[flux]²/Var(flux)` (= 1/C_b², R44 eq 5.1). Photon count n
swept {1e3…1e5}; per-bank covariance noncentrality λ(n) with shot ∝1/n folded into V0;
T_det(n) = 25/λ(n).

## Measured
| aperture frac | N_eff | T_det (n=1e4) | T_det/N_eff |
|---|---|---|---|
| 0.40 | 46.4 | 316 | 6.82 |
| 0.55 | 90.3 | 560 | 6.20 |
| 0.70 | 214.2 | 1308 | 6.11 |
| 0.85 | 369.8 | 2299 | 6.22 |
| 1.00 | 810.6 | 4180 | 5.16 |

- **P_Neff_grows: PASS** — N_eff rises 46 → 811 (~17×) as the aperture opens 2.5×
  (super-area growth: opening the aperture both enlarges the illuminated object area and
  shrinks the speckle grain, both raising the grain count). Étendue sets N_eff.
- **Étendue-conductance law holds (T_det ∝ N_eff)** — T_det/N_eff = 5.2–6.8
  (saturated-regime CV **0.25**): at fixed photon budget, detection latency scales with
  the effective grain count, i.e. the resource is (étendue, independent banks) and
  photons stop helping once n ≫ N_eff.
- **P_Tdet_over_Neff_collapses: PARTIAL** — the full universal collapse of T_det/N_eff
  onto a single function of r=n/N_eff has normalized CV **0.64** (vs raw CV 0.81 at
  fixed n): the rescaling helps but is not a clean single-curve collapse, because the
  covariance-channel latency depends on the per-aperture speckle structure and the
  change coupling, not on r=n/N_eff alone.

## Reading
The core R44 §5.4 claim — étendue fixes the number of independently integrated speckle
grains N_eff, which governs covariance detection latency — is confirmed: N_eff ∝ étendue
and T_det ∝ N_eff (25% scatter). The stronger claim of an exact universal `T_det/N_eff =
f(n/N_eff)` collapse is only partially supported on this twin (the covariance channel
carries structure beyond the scalar r). Honest PARTIAL. Nothing here touches the Letter
or sealed artifacts.
