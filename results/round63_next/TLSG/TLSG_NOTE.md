# TLSG — Three-Ledger Scaling Gate

**VERDICT: TLSG_PASS — all 7 frozen bars pass. The Paper-2 Certified-Blindness /
Three-Ledger architecture FREEZES.**

Preregistered exam for R45 (issue #36; ruling §B11.7 + §C1.2). Thresholds frozen by the
committed R45 ruling; runner built, run ONCE, verdict — no repair round. Seed 20260724.
`tlsg_A_three_ledger.py` (arm A) + `tlsg_B_calibration_transfer.py` (arm B) →
`TLSG_partA.json`, `TLSG_partB.json`, combined `TLSG_RESULTS.json`. Colab pro2, 2× L4.

## Per-bar table (frozen vs measured)

| # | Arm | Bar | Frozen | Measured | Verdict |
|---|---|---|---|---|---|
| 1 | A §B11.7 | matched 50%-power contours collapse vs T·λ | CV ≤ 0.20 | **CV = 0.125** | **PASS** |
| 2 | A §B11.7 | blind contours collapse vs T·λ/√p AND fail vs T·λ | CV_scaled ≤ 0.25 AND CV_raw > 0.25 | **CV_scaled = 0.124**, CV_raw = 0.628 | **PASS** |
| 3 | A §B11.7 | relative estimation risk collapses vs T·λ/p | CV ≤ 0.25 | **CV = 0.035** (mean risk·T·λ/p = 0.985) | **PASS** |
| 4 | A §B11.7 | fitted p-exponents (matched/blind/estimation) | 0±0.15, ½±0.15, 1±0.15 | **−0.004 / 0.444 / 1.010** | **PASS** |
| 5 | A §B11.7 | KT6 below blind power 0.2 at FPR 0.05 (optimal blind stat) | power < 0.2 | **0.051** (10×: 0.056, 100×: 0.134) | **PASS** |
| 6 | B §C1.2 | ≥95% coverage of σ_d(L̂)+δ envelope, no code exceeding | coverage ≥ 0.95 | **coverage = 1.000** (σ₆₄=3.4e-8, δ=1.07e-2, worst fresh S₆₄ leak 2.0e-6) | **PASS** |
| 7 | B §C1.2 | ≥64 conjugate-plane dims certified < 1e-4 on fresh draws | ≥ 64 dims | **min = 80** (calib d@1e-4 = d@1e-5 = 80) | **PASS** |

## Arm A — three-ledger scaling (§B11.7)
Nuisance-profiled Gaussian covariance family Y~N(0,Σ_ε), Σ_ε = I + Σθ_jB_j, efficient
symmetric tangent basis (½)tr(B_jB_k)=δ_jk, λ=‖θ‖². Grid M∈{8,16,32,64} × p_eff∈{10,36,
136,528} (controlled tangent restriction, p_eff ≤ M(M+1)/2), N_MC=1500, λ=0.15, full
quadratic-form score simulation (not the LAN reduction assumed a priori). The three
ledgers separate exactly as `1 : √p : p`:
- **matched** (Ledger I) known-direction score → 50%-power at T·λ ≈ 2.6, **p⁰** (exponent −0.004);
- **blind** (Ledger II) score-norm ‖Z_T‖² → 50%-power at T·λ/√p ≈ 2.85, **p^½** (exponent 0.444; slightly below ½ by the expected finite-p χ² offset, inside tolerance) — and it demonstrably **fails** to collapse vs T·λ (CV 0.628);
- **estimation** (Ledger III) θ̂=Z_T/√T → risk·T·λ = p exactly (CRB), **p¹** (exponent 1.010).
The KT6 near-contact operating point (signal/floor 5.25e-3, p=136) maps to the optimal
blind statistic's detectability T·λ/√(2p) → blind power 0.051 ≈ FPR ≪ 0.2 (robust to 10×/100×),
confirming the KT6 negative for the right reason (below optimal blind certifiability).

## Arm B — calibration transfer (§C1.2, Eq. 4.4)
Whitened leak operator L_w = L_be B^{-1/2} at the conjugate plane (z1=0); calibration L̂_w
(nominal hardware) → bottom-64 subspace S₆₄, σ₆₄ = 3.4e-8. Fresh calibration/test split
with hardware-law perturbations (fill ±5%, mirror pitch ±3%, relay defocus 0–0.5 mm);
δ = 1.07e-2 estimated on the 24-draw calibration ensemble. On 24 disjoint FRESH draws:
**100% simultaneous coverage** of the σ₆₄+δ envelope (worst S₆₄ leak 2.0e-6, six orders
below the envelope), and **all draws keep 80 code dimensions below 1e-4** (≥ the required
64) at the conjugate plane — the single-plane blindness-capacity point survives fresh
physical-law draws.

## Consequence
Every frozen bar passes; per the binding R45 verdict logic the **Paper-2
Certified-Blindness / Three-Ledger architecture freezes** (Statistical Noether Wall +
Blindness Capacity + Three-Ledger Certifiability, with Jet Transmutation as corollary).
Bench-phase only; nothing here touches the Letter or sealed artifacts.
