# Study-2 no-gate modes — descriptive summary (pre-primary)

Generated 2026-07-18 (Berlin), after completion of the five non-confirmatory
modes (robustness / controls / fluxmap / separation / scattering). The k=16
primary was still running when this file was written; nothing here feeds any
gate decision (all five modes are preregistered NO-GATE per round-8 ruling).
All numbers: RQL arm, PSNR_rad, means over completed rows.

## 1. Fluxmap (1 rep/family × 3 seeds; ν=2000 fixed budget)

| k | C_u | ρ=0.05 | 0.1 | 0.3 | 0.6 | 1.0 | 2.0 | Γ(0.6, 2000) |
|---|-----|--------|-----|-----|-----|-----|-----|---------------|
| 512 | 0.040 | 8.67 | 8.81 | 9.01 | 9.32 | 9.51 | 9.74 | 1.09 |
| 32  | 0.222 | 11.46 | 12.57 | 14.58 | 15.94 | 16.63 | 17.36 | 6.07 |
| 16  | 0.316 | 11.85 | 12.46 | 13.05 | **13.41** | 13.32 | 13.18 | 8.66 |
| 1   | 1.270 | 17.12 | 18.77 | 19.64 | **19.82** | 19.61 | 19.36 | 34.79 |

- C_u ladder behaves as designed (0.04 → 0.22 → 0.32 → 1.27); Γ crosses ≫1
  for all sparse rungs ⇒ photon-limited structural regime reached.
- Fixed-budget high-flux gain grows with sparsity ladder as predicted for
  k=512→32; k=16 and k=1 show an **interior ρ optimum near 0.6–1.0**
  (quality declines again at ρ=2), i.e. the phase-diagram boundary is visible
  descriptively.
- Note (honest): k=32 exceeds k=16 in absolute PSNR at long dwell —
  conditioning (κ≈n/k) opposes contrast; the ladder is non-monotone in
  absolute quality. Gate remains k=16 per freeze; this is reported as a
  descriptive observation only.

## 2. Mechanism separation (k=16, 6 reps × 3 seeds) — the closest-prior killer

| flux mode | ν | safe (ρ=0.05) | fast (ρ=0.6) | Δ | S_det safe/fast |
|---|---|---|---|---|---|
| normalized (Study-2 mechanism) | 20 | 8.57 | 8.84 | +0.27 | 0.0009 / 0.0073 |
| normalized | 200 | 8.49 | 10.30 | **+1.81** | 0.0093 / 0.0723 |
| normalized | 2000 | 10.26 | 11.36 | **+1.10** | 0.0927 / 0.7227 |
| unnormalized (prior-art mitigation) | 20 | 2.74 | 4.47 | +1.73 | 0.0002 / 0.0024 |
| unnormalized | 200 | 4.62 | 7.48 | +2.86 | 0.0025 / 0.0235 |
| unnormalized | 2000 | 7.45 | 10.57 | +3.12 | 0.0246 / 0.2343 |

- The Study-2 claim survives the separation: **at equal mean detector load and
  equal dose** (normalized), high flux still gains +1.1…+1.8 dB at ν≥200.
- The unnormalized route (what Zhang/Liu/Han-style sparsity does) has larger
  *relative* deltas but strictly worse *absolute* quality at every ν
  (e.g. 7.45 vs 10.26 at safe, ν=2000): it wins by starving the detector.
  This is exactly the defensible-difference table of the round-8 ruling.

## 3. Fixed-budget Δ(ρ0.6−ρ0.05) at ν=2000, per-image (24 confirmatory images)

| arm | k | mean Δ | median Δ | positive |
|---|---|---|---|---|
| robustness (no gate) | 32 | **+5.10 dB** | **+4.17 dB** | 24/24 |
| control | 512 | +0.79 | +0.65 | 24/24 |
| control | 1 | +2.59 | +3.64 | 21/24 |

(The k=16 primary fixed-budget secondary is computed only by the frozen
analyzer after the primary grid completes.)

## 4. Dynamic scattering (k=16, CV(α)=0.20, 6 reps × 3 seeds)

| ν | safe | fast | Δ |
|---|---|---|---|
| 20 | 8.71 | 9.37 | +0.66 |
| 200 | 8.71 | 10.55 | +1.84 |
| 2000 | 9.96 | 12.28 | **+2.32** |

High-flux advantage survives (indeed grows under) dynamic scattering —
directly relevant to the Hao/Chen turbid-water setting; reported in S_det
units in the paper panel (f).

## Paper-panel mapping

- Panel (c): §1 C_u/Γ ladder — data complete.
- Panel (e): §3 fixed-budget by occupancy — awaits k=16 primary for its bar.
- Panel (f): §4 scattering — data complete.
- Mechanism-separation paragraph: §2 — data complete.
