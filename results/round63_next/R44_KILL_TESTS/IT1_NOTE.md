# IT1 — Two-wall + jet-order invariance (internal divergence higher #1, WAVE-2 first)

**Verdict: PASS — a SECOND EXACT WALL (energy-conservation type) is confirmed.**
`IT1_two_wall_jet_invariance.py` + `.json`. float64, 2 s on one L4.

## Arm (ii) — second exact wall
A phase-only object change `t=exp(iεφ)` is lossless (|t|=1). A **lossless full-aperture
bucket** integrates all the light, so by Parseval + unitary propagation
`B_full = ∫_all |prop(E·t)|² = ∫|E|²` — **independent of ε at every order**.

| z2 (mm) | full-aperture null (max rel) | clip slope | clip response @2% |
|---|---|---|---|
| 0.5 | **2.18e-16** | 1.033 | 3.7e-15 |
| 1.0 | **2.18e-16** | 0.982 | 5.2e-15 |
| 5.0 | **2.18e-16** | 0.970 | 3.1e-13 |
| 20.0 | **2.18e-16** | 1.000 | 1.3e-9 |

- **P_second_exact_wall: PASS** — full-aperture null is machine-zero (2.18e-16) at
  every z2 and every ε. This is a genuine second exact wall of a *different
  conservation type* than the pupil/code walls: energy conservation makes a
  phase-only object modulation invisible to a lossless bucket at every order.
- **P_clip_slopes_physical: PASS** — the NA-clipped bucket responds at O(ε¹)
  (slopes 0.97–1.03, all in [0.9,1.15]): clipping breaks the wall at first order,
  the "eps^1 restoration under clipping" prediction.
- **z2 trend MONOTONE** — the clip response grows ~6 orders of magnitude from z2=0.5mm
  (3.7e-15, essentially null near contact) to z2=20mm (1.3e-9): phase-to-intensity
  conversion develops with propagation, so the clip bucket "sees" the phase change
  only after enough lateral energy redistribution. This is the geometric jet-order
  interpolation seed (a flagship follow-up candidate).

## Bearing
Confirms the merged adjudication's conjectured "possible second exact wall of a
different conservation type: phase-only object change invisible to a lossless
full-aperture bucket at every order." Two-wall symmetry (the code/pupil support wall
of KT1 + this energy-conservation wall) is real and enters the next-paper frame.
The near-contact clip-null (≤1e-14 at z2≤1mm) also means an NA-limited bench near
conjugate is nearly phase-blind. Nothing here touches the Letter or sealed artifacts.
