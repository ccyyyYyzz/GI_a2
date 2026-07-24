# KT1 — Seven-variant shell sweep (R44 §1.4, innovation #1, score 125)

**Verdict: R44 innovation #1 (Finite-Aperture Wall Impossibility / Pupil-Hardened
Wall) CONFIRMED, with the exact field-vs-detector wall distinction resolved.**

Deterministic, no diffuser (z2=0), M=12 codes, z1∈{0,1,5,10,20} mm, float64
(complex128), 36 s on one L4. `KT1_shell_sweep.py` + `KT1_shell_sweep.json`.

Shell-resolved mean-channel leak `L_k = ||J U_k||_F / ||J U_[0,kp]||_F` over
Chebyshev shells k=1..31, seven nested optical variants. We report BOTH a
detector-level `L_k` (block-summed Jacobian over the finite 512 object window) and
a FIELD-level spectral leak (|FFT(D)| energy beyond a shell, no object window), which
is what cleanly separates the two walls.

## Frozen predictions vs measured

| Pred | Statement | Measured | Verdict |
|---|---|---|---|
| P1 | variant 1 (ideal periodic code) machine-zero for k>kp, every z1 | max L_k = **2.1e-15** | **PASS** |
| P2a | hard Fourier pupil → FIELD-level D energy beyond 2kR < 1e-12 | max = **6.6e-16** (all kR∈{5,6,8}, all z1) | **PASS** |
| P2b | detector-level wall beyond 2kR | **finite-object-window-limited** (below) | see note |
| P3 | full twin nonzero tail (>1e-10) above 2kp | det **5.3e-2**, field **1.6e-1** | **PASS** |
| P4 | propagation never moves the ideal edge; field pupil → z1-independent leak | V1 edge=5 ∀z1; V6 leak rel-spread 1.3e-3 (kR5), 5.0e-3 (kR6), 2.0e-2 (kR8) | **PASS (kR5,6) / PARTIAL (kR8)** |

## The physics the ladder exposes

- **V1 ideal / V2 aperture / V3 hold**: machine-zero (~4e-16) beyond kp at z1=0 —
  at contact `D = c·s = s` is exactly band-limited regardless of aperture/hold.
  V2/V3 leak (~2e-2) appears only at z1>0 (propagation makes fields complex and
  breaks the complementary cancellation).
- **V4 mirror sinc / V5 full twin**: break the wall even at z1=0 (floor ~3.8e-3),
  = the DMD pixelation floor (cf. sealed T1 6.8e-3). V5≈V4: the mirror-sinc model is
  complete; the "full twin" adds no ingredient beyond it.
- **V6 hard pupil**: the signed-intensity FIELD is band-limited to 2kR **exactly**
  (6.6e-16, machine-zero) at every z1 and kR — the pupil changes the mathematical
  class (Paley–Wiener → compact difference-set support, R44 Thm 1 case 2). The
  detector-level leak beyond 2kR is **z1-independent** (1.19e-2 kR5 → 6.2e-3 kR8,
  invariant across z1 to 3 digits for kR5/6) — reproducing the coordinator's
  measured "hard field pupil → z1-independent leak" finding.
- **The finite OBJECT window is the residual wall-breaker**: block-summing D over
  the finite 512 active region (= the physical bucket integrating a finite object)
  reconvolves the pupil-limited spectrum → a detector tail beyond 2kR that the
  exact field pupil cannot remove. A naive separable **Tukey guard (V7) reduces it
  ~8–40×** (kR5 1.5e-3, kR8 1.6e-4) **but never to 1e-12** → an exact DETECTOR wall
  needs a concentration-grade (DPSS/Slepian) guard. This is exactly R44 binding
  conclusions #1/#2 and the coordinator's "naive separable windows are kills"; the
  DPSS route is queued as IT5.

## Bearing on the 125-score theorem

Confirms all three legs: (i) a finite DMD field has **no exact open spectral wall**
(P3 noncompact pedestal); (ii) a **hard Fourier pupil restores an exact difference-set
field wall** to 2kR (P2a machine-zero ∀z1); (iii) **exactness at the detector
additionally requires guarded sampling** (P2b + V7), not achievable with a naive
window. Frequency placement alone does **not** restore a wall (P3). Nothing here
touches the Letter or sealed artifacts.
