# SCRAMBLE_EXT — numerical ledger (fully scrambling dynamic medium)

**Date:** 2026-07-24 · **Engine:** `scramble_toy.py` (torch-CUDA, RTX 4060, 598 s) ·
**Data:** `SCRAMBLE_RESULTS.json`, figure `SCRAMBLE_DUEL.png`. **Derivation:** `DERIVATION.md`.
Grid `32×32` (`N=1024`), pattern band `k_p=3`, grain band `k_grain=13` (`≫ k_p`), `M=36`
band-limited codes, `n=10⁴` photoelectrons/bucket. Exact Kronecker speckle generator
`E=B(Z O^{1/2})` (fully-developed speckle, iid circular-Gaussian transmission, grain low-pass).
Numbers only. No git commit.

---

## Verdict

> **The beyond-band SENTINEL survives complete scrambling; beyond-band IMAGING does not.**
> Under a fully scrambling medium the **mean-wall gets stronger** (the speckle envelope is flat —
> the first-moment channel sees only the scene DC term, mean-null verified to `6.6e-19`), the
> **covariance channel responds to any beyond-band change** through a single scalar `Q(x)=xᵀGx`
> whose spectral aperture spans the **whole grid** (`Ĝ` nonzero on `1024/1024` bins), and the
> `δᵀGδ>0` floor guarantees a positive response with **no blind direction**. The four-way
> specificity split ports intact (mean detector at chance, AUC `0.50`). **But** the same-epoch
> covariance collapses to **rank one** (`99.99 %` of its Frobenius energy lies in the single
> `Q·(O∘O)` direction), so *reconstruction* of the beyond-band pattern dies — full scrambling is a
> **detect-but-cannot-image** regime. Detection is real but bank-limited: a **physically-generic
> (coherent) 5 % beyond-band change is caught at `d′=4.1`, AUC `0.997` in `8192` banks (≈1.7 min)**;
> the spectrally-orthogonal worst case rides the `ε²` floor and needs `~10⁶` banks (a rate gap, not
> a floor).

---

## 1. Aperture — support gigantic, rank one

| quantity | value | meaning |
|---|---|---|
| pattern band `k_p` | 3 | DMD/modulator wall `B_p` |
| grain band `k_grain` | 13 | speckle-grain spectrum `B_grain ≫ B_p` |
| `Ĝ(k)` nonzero bins | **1024 / 1024** | covariance support aperture `= 2·B_grain` reaches the **whole grid** |
| `C(0)` | 0.517 | peak grain coherence (passband fraction) |

The grain kernel `Ĝ=FT{|C|²}` is nonzero on **every** frequency — the covariance channel's *support*
pierces the wall completely (Eq. 3.5). What collapses is not the support but the **rank** (§3).

## 2. Validation A — the mean-wall is stronger (Theorem 1)

| test | result | expectation |
|---|---|---|
| mean speckle envelope `μ_i(r)` spatial CV | **0.0156** (4 codes) | = MC floor `1/√n_env = 0.0158` → **flat** |
| envelope non-DC power fraction | **2.4e-4** | pure Monte-Carlo residual (→0) |
| mean-null `‖Δm‖/‖m‖` for DC-free beyond-band `δ` | **6.6e-19** | machine zero → mean channel **exactly blind** |

The empirical mean envelope is flat to the Monte-Carlo noise floor: `E[b_i]` depends on the scene
**only** through `x̂(0)`, so every non-DC change — in-band or beyond — is invisible to the first
moment. The wall is no longer `B_p`; it is `{0}`. (The observed `1.000` beyond-band baseline of R39
is here promoted to a total mean blindness.)

## 3. Validation B — the covariance is the single scalar `Q(x)` (Eq. 3.3, rank collapse)

| test | result | expectation |
|---|---|---|
| `Var[b_i]/(O_ii² Q(x))` median over codes | **1.019** | = 1 (Siegert + factorisation) |
| fraction of `Cov` Frobenius energy in the `O∘O` direction | **0.99992** | = 1 → covariance is `Q·(O∘O)`, **rank one in matrix space** |
| top eigenvalue fraction of sample `Cov` | 0.826 | leading `O∘O` mode dominates |
| `ΔQ` predicted vs measured (Eq. 4.1) | **rel err 2.4e-15** | exact |

The entire `M×M` same-epoch bucket covariance is one scalar `Q(x)=xᵀGx` times the known Hadamard
square `O∘O` of the code Gram — confirmed to `10⁻⁴`. The scene's reconstruction Fisher from this
channel is therefore **rank one** (Eq. 3.6): only the direction `Gx` is identifiable. **The
mode-coverage curse returns in its most extreme form — beyond-band imaging does not survive full
scrambling.**

### The sentinel floor (`δᵀGδ > 0`, no blind direction)

| `δ` type | cross `2xᵀGδ` | quad `δᵀGδ` | `ΔQ/Q` (ε=0.05) | regime |
|---|---|---|---|---|
| **ortho** (spectrally disjoint from `x`) | 0.0031 | **0.345** | 1.86e-3 `≈ ε²` | worst case: only the guaranteed floor responds |
| **coherent** (`x`'s own beyond-band texture) | **2.978** | 0.411 | 1.81e-2 `≈ 0.9 ε` | generic physical change: linear-in-ε response |

For the orthogonal `δ` the coherent term nearly vanishes (the two are DCT-disjoint) and detection
rides **only** `δᵀGδ = 0.345 > 0` — the sign-definite floor that cannot be nulled. For a generic
(coherent) change the `2xᵀGδ ∝ ε` term dominates and is `~10×` larger.

## 4. Validation C — the sentinel duel (H0 vs H1, efficient matched-score LR)

Detector: efficient matched-score `t = ⟨Σ₀⁻¹(O∘O)Σ₀⁻¹, Ĉov⟩` (the FOG `W=V₀⁻¹dV V₀⁻¹` statistic;
whitening by `Σ₀⁻¹` de-correlates the code looks). Mean-channel detector at chance for **both** `δ`
types (AUC **0.505 / 0.502**) — the wall holds operationally.

| `δ` | ε | `T_eff`=512 | 2048 | 8192 | 32768 | 65536 | `T_req(d′=5)` |
|---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| **coherent** | 0.05 | d′0.99 A.762 | d′2.22 A.943 | **d′4.10 A.997** | — | — | **1.18e4 banks ≈ 152 s** |
| **coherent** | 0.02 | d′0.38 A.605 | d′0.85 A.724 | d′1.60 A.887 | — | — | 7.82e4 banks ≈ 1000 s |
| **ortho** (worst) | 0.05 | — | — | d′0.36 A.601 | d′0.99 A.763 | d′1.23 A.809 | 1.02e6 banks ≈ 3.6 h |

(A = AUC; TPR@FPR5 %: coherent ε=0.05 reaches **0.99** at 8192 banks; ortho reaches 0.42 at 65536.)

**`d′ ∝ √T_eff` in every case** (clean fits, see figure), so the `0.44–0.48`-style stalls of the
thin-screen line are confirmed here as **finite-T rate effects, not floors** (R39 Thm 3): the fitted
slope `a` extrapolates a finite `T_req` for any target `d′`.

### Effective looks `M_eff` — capped by the in-band code DOF

From `d′ = (ΔQ/Q)·√(M_eff·T_eff/2)` and the fitted slopes:

| duel | `M_eff` (from slope) |
|---|:--:|
| coherent ε=0.05 | 12.9 |
| coherent ε=0.02 | 12.2 |
| ortho ε=0.05 | 14.3 |

All land at `M_eff ≈ 13`, matching the **in-band code degrees of freedom** `(k_p+1)²−1 = 15`. This is
the precise cost of the rank collapse: the `M=36` codes give only `~13` independent covariance looks
because the informative code directions live in the band-limited illumination subspace. Adding codes
beyond the in-band DOF buys nothing — **independent medium banks, not codes or photons, are the
currency** (the FOG rate law, reconfirmed under full scrambling). The naive `⟨Ĉov,O∘O⟩` projection
gets only `M_eff≈1`; the efficient whitened detector recovers the full `~13`.

---

## 5. Honest ledger — what survives, what breaks

**Survives (upgrades to universal — thin screen → complete scramble):**
- **Mean-wall.** Strictly stronger: band `{0}` vs `B_p`. Mean-null `6.6e-19`. (§2)
- **Covariance response to any beyond-band change.** `Ĝ` spans the whole grid; `δᵀGδ>0` guarantees a
  sign-definite positive response for **every** nonzero `δ` — including one beyond band, DC-free, and
  `G`-orthogonal to `x`. **No blind direction for detection.** (§3)
- **Specificity split.** Mean detector at chance (AUC 0.50); the beyond-band change lives entirely on
  the covariance channel. Temporal-lag structure carries `Q` coherently across lags (analytic §4.2),
  so the FOG four-way `mean / beyond-cov / medium-amp / medium-τ` discriminator ports unchanged.
- **Detectability of a physically-generic change.** Coherent 5 % beyond-band anomaly: `d′=4.1`,
  AUC 0.997, TPR@FPR5 %=0.99 in **8192 banks (≈1.7 min)**; `T_req(d′=5)≈1.2e4` banks.

**Breaks (does NOT upgrade — remains a finite-memory capability):**
- **Beyond-band imaging / reconstruction.** Full scrambling collapses the covariance to rank one
  (`0.99992` energy in one direction); the recoverable scene aperture is a **single scalar** `Q(x)`.
  You learn that the beyond-band **energy footprint** changed, not **what** it changed to. Two `δ`'s
  with equal `‖δ‖_G` and equal `Q`-projection are indistinguishable. Imaging survives only with a
  residual memory effect (the analytic rank-`O(r_ME²)` interpolation, §3.3 — untested here).

**Sample complexity (the cost, stated plainly):**
- Detection is `√T`-limited with `M_eff ≈` in-band code DOF (`~13`), prefactor `ΔQ/Q`.
- Coherent (linear-in-ε): `T_req(d′=5) ≈ 1.2e4 (ε=5%) / 7.8e4 (ε=2%)` banks — minutes to ~15 min.
- Orthogonal worst case (`ε²`): `T_req ≈ 1.0e6` banks — hours. Real, but a **rate gap, not a floor**.

**Where it could still break (flagged, untested):**
- Coarse speckle (`k_grain ≲ k_p`): `Ĝ` would roll off before `δ`'s support, weakening `δᵀGδ`.
  Here `k_grain=13 ≫ k_p=3`, the strong-scattering regime; the claim is scoped to fine speckle.
- The memory-effect interpolation (§3.3) is analytic only. A shift-memory-effect transmission model
  is the clean next probe to map exactly how much imaging returns as `r_ME` grows.

---

## 6. Bottom line for the program

The flagship claim splits on an honest boundary the numbers draw sharply:

- **Beyond-band DETECTION (the sentinel) is universal.** It survives complete scrambling — the
  strongest-possible dynamic complex medium — with the wall *stronger* and the response *sign-
  definite*. This is the Nature-Photonics/PRL-tier claim, now demonstrated end-to-end (mean-wall +
  covariance `Q` + duel + ROC) through full speckle, not just a thin screen.
- **Beyond-band IMAGING is not universal.** It is a memory-effect (finite-scattering) capability;
  complete scrambling reduces it to a rank-one energy read-out. The manuscript must state this — it
  is the correct and defensible scope, and it is exactly what R39's aperture law predicts in the
  `Γ→diagonal` limit.

Files: `DERIVATION.md`, `scramble_toy.py`, `SCRAMBLE_RESULTS.json`, `SCRAMBLE_DUEL.png`.
