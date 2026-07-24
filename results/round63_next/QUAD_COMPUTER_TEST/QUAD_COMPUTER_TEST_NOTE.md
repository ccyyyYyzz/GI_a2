# QUAD_COMPUTER_TEST — the law-only quadratic computer (R46 §4.2–4.3)

**VERDICT: KILL_GENERAL_CLAIM (preregistered negative).** The idealized random-law quadratic
kernel is exact and beats every linear statistic, but *physically admissible* circulant/Toeplitz
grain kernels collapse the quadratic feature cone from `d(d+1)/2 = 300` to `~13` (power-spectrum
only) and cannot compute a generic quadratic boundary. Per R46 §6.4 the general "law-only
quadratic computer" flagship is **killed**; it survives only as a power-spectrum feature machine.

Engine: SCRAMBLE_EXT law-only covariance convention. Run: `quad_computer_test.py`, full scale
`d=24`, numpy/CPU, 1.4 s (Colab GPU offers nothing for this job, so it was run locally at full
scale — deterministic and reproducible).

## Frozen bars (declared before running) — measured

| Bar | Frozen | Measured | Verdict |
|---|---|---|---|
| 1 kernel identity `E[z z]=2(x·y)²` (Isserlis) | median diagonal rel err < 0.05 at L=1e5 | **0.0085** (Frobenius 0.022) | **PASS** |
| 2 feature-Gram error ~ `L^{-1/2}` | log-log slope ∈ [−0.6,−0.4] | **−0.477** | **PASS** |
| 3 quadratic boundary a linear stat can't solve | law-only & digital > 0.85; best linear < 0.62; \|law−dig\|<0.05 | **law 0.903, digital 0.903, camera-speckle 0.903, linear 0.495** | **PASS** |
| 4 admissible-cone KILL test | general claim KILLED if collapse-ratio < 0.25 AND admissible fails generic quadratic | **cone rank: random-rank1 300, circulant 13 (ratio 0.043); admissible generic-quadratic acc 0.549, power-spectrum acc 0.955** | **KILL_GENERAL** |

## Reading

- **Idealized version is exact and non-trivial.** With unconstrained random rank-one laws
  `G_l = a_l a_lᵀ`, the centered feature `z_l(x)=(a_lᵀx)²−‖x‖²` satisfies the exact Isserlis
  identity `E[z_l(x)z_l(y)] = 2(xᵀy)²` (diagonal rel err 0.85%), the empirical feature Gram
  converges as `L^{-1/2}` (slope −0.477), and the law-only covariance features solve a generic
  quadratic boundary (90.3%) that the best linear bucket statistic provably cannot (49.5% ≈
  chance) — matching a digital random-feature and a camera-speckle random-feature arm exactly.
  This is classical Gaussian moment algebra (R46 §6.4), so it was expected to pass.

- **Physics kills the general claim.** The decisive test replaces unconstrained `G_l` with
  physically admissible circulant/Toeplitz grain kernels (the real scattering medium). For a
  circulant `G_l`, `x^T G_l x = Σ_k g_l(k)|x̂(k)|²` is a *linear readout of the power spectrum*.
  The quadratic-form cone therefore collapses from the full `d(d+1)/2 = 300` dimensions to
  `~d/2 = 13` (Hermitian-symmetric power-spectrum DoF) — a collapse ratio of **0.043**. The
  admissible arm accordingly **fails** a generic quadratic boundary (54.9% ≈ chance) while
  perfectly solving a power-spectrum-expressible one (95.5%).

- **Bank cost also overwhelms.** Spanning the full random cone needs `~d(d+1)/2` independent law
  kernels = independent medium banks; at that point the "0-D law-only" advantage over a matched
  digital baseline is gone (R46 §6.4 second kill condition).

**Consequence.** Per the R46 §4.3 / §6.4 kill conditions the general law-only quadratic-computer
flagship is retired: physically admissible kernels span only a very low-dimensional (power-
spectrum) cone. The idealized exact-kernel result stands as classical random-feature/reservoir
computing under a new noun. Bench-phase only; nothing here touches the Letter or sealed artifacts.
