> ⚠️ **SUPERSEDED (buggy shot term).** This map used `R_shot ∝ signed Px` (near-zero → Fisher
> ~1.9× inflated, d' ~1.38× / T_det ~1.9× optimistic). Use **`DETECTION_MAP_CORRECTED.md`** (physical
> `|P|·x` throughput shot) and **`DETECTION_ROC.md`** (MC-validated) instead. Kept for the record only.

# DETECTION LIVING-REGION MAP — the demand-axis divergence

**Date:** 2026-07-24 · **Engine:** `detection_map.py` + `detection_robust.py` (reuse the validated
P1 / living-region profiled-Fisher analytics; **NO reconstruction**) · **Data:** `DETECTION_MAP.json`,
`DETECTION_ROBUST.json`. Ran locally on the RTX 4060 (~1 min). Fixed: `N=4096`, `M=128`, `k_p=5`,
`n=1e4`, same 81-cell grid (σ_f × spectrum-shape × k_w/k_p × claim-shell).

## Verdict

> **DETECTION_REGION_FOUND.** The estimation map is `REGION_EMPTY` — it fails **only** the
> mode-coverage bar (`f_rec ≥ 0.70`), while *aggregate* information passes everywhere. **Detection
> needs only the aggregate**, and the map is green: a beyond-band anomaly of **energy ≥ 2 % is
> detectable at `d'≥5` within 4096 epochs in 100 % of the 81 cells** (even worst-case/adversarial
> direction at ≥5 %); a **1 % anomaly is detectable (energy-spread) in 81 % of cells** (100 % at
> claim ≤ 1.5 k_p). Best cell: **≈1335 epochs** for a 1 % anomaly (`d'≈8.8`). The **mean/row channel
> is exactly blind** to every such anomaly (specificity theorem, `‖P U_β‖ ~ 1e-13`), so this is a
> genuine **covariance-only, F4-immune demand-tier capability** — the fresh-pattern comparator
> cannot detect it. Detection **inherits and exceeds** the estimation robustness: a 10 % medium-law
> basis rotation changes `d'` by **+0 %** (aggregate Fisher trace is invariant to 0.04 %).

## 1. Detection model and statistic

`H0`: declared scene `x`. `H1`: `x + δ`, with `δ` supported **strictly beyond the modulator band
and inside the covariance aperture** (`δ ∈ U_β`, `k_p<|k|≤k_claim≤2k_p`), energy `ε = ‖δ‖/‖x‖`.
The lag-covariance likelihood ratio (Gaussian) has deflection
`d'(T_eff)² = T_eff · δᵀ J̄_B δ`, where `J̄_B` is the **per-epoch profiled beyond-band Fisher**
(`J_B = T_eff·J̄_B` scales exactly linearly, so `J̄_B` is `T_eff`-independent — a fixed
per-cell geometry object). Two directions:
- **least-favorable (worst-mode):** `δ` aligned with the min-eigenvector → `d'² = T_eff·λ_min·ε²‖x‖²`;
- **energy-spread (average):** `δ` isotropic → `d'² = T_eff·λ_mean·ε²‖x‖²`, `λ_mean = tr(J̄_B)/d_β`.

`T_det(ε) = 25 / (λ·ε²‖x‖²)` = epochs for **strong detection `d'≥5`**.

## 2. The specificity theorem (mean-channel null) — numerically exact

A beyond-band `δ` changes the **mean** bucket by `dm = P·δ = P·U_β·c`. Because the codes are
band-limited (`|k|≤k_p`) and `U_β` is beyond-band (`|k|>k_p`), Fourier orthogonality gives
`P·U_β = 0` **exactly**:

| claim | ‖P U_β‖_F | max\|·\| | rel. to ‖P‖ |
|---|---|---|---|
| 1.25 k_p | 4.4e-14 | 3.6e-15 | 2.0e-16 |
| 1.5 k_p | 7.1e-14 | 6.7e-15 | 3.2e-16 |
| 1.8 k_p | 1.2e-13 | 8.9e-15 | 5.7e-16 |

The mean/row channel — and therefore any fresh-known-pattern or classical-averaging route — carries
**zero** detection information about a beyond-band anomaly, and generates **zero** in-band false
alarm from it. All detectability lives in the **covariance** channel, and it is **specific**. This
is the demand-tier analogue of the E8 super-resolution physics wall: comparator-proof (F4-immune).

## 3. The detection region (T_det for d'≥5)

Fraction of the 81 cells with `T_det ≤ 4096`:

| anomaly energy ε | energy-spread (average) | worst-mode (adversarial) |
|:--:|:--:|:--:|
| 0.5 % | 0 % | 0 % |
| **1 %** | **81 %** | 15 % |
| **2 %** | **100 %** | 67 % |
| 5 % | 100 % | **100 %** |

By claim shell (energy-spread, 1 %):

| claim | avg-direction pass (1 %) | worst-mode pass (1 %) | T_det_avg(1 %) range | d_β |
|---|:--:|:--:|:--:|:--:|
| 1.25 k_p | 100 % | 44 % | 1335 – 3214 | 48 |
| 1.5 k_p | 100 % | 0 % | 1810 – 3693 | 104 |
| 1.8 k_p | 44 % | 0 % | 2950 – 5136 | 240 |

Reading: the **energy-spread** anomaly (the physically relevant "does something beyond-band exist")
is detectable at 1 % across **all narrow/mid claims** and everywhere at ε≥2 %. The **worst-mode**
(an adversary hiding energy in the near-null modes) is only reachable at the narrow claim (1.25 k_p)
for a 1 % anomaly, and everywhere by ε=5 % — the same "some modes are near-null" structure that
sinks *estimation* mode-coverage, softened for *detection* because detection sums over modes. The
best cell (**k⁻² spectrum, k_w=1·k_p, claim 1.25, σ_f=0.3**) detects a 1 % anomaly in **≈1335
epochs** at `d'≈8.8`; the same cell was the estimation map's best residue (T_req(0.30)≈197).

## 4. Detection vs estimation — the divergence

The estimation map is `REGION_EMPTY` on the **mode-coverage** flagship bar (witness `f_rec` never
reaches 0.70) yet passes the **aggregate** bars everywhere (`NRMSE_CRB` and `T_req(0.30)≤4096` in
all 81 cells, `T_req(0.30)≈200–1000`). **Detection consumes only the aggregate**, so it lives where
imaging's flagship fails:

- A **2–5 % beyond-band anomaly is universally detectable** (all 81 cells) at cost comparable to or
  below imaging: `T_det_avg(2%) ≈ 330–1280`, `T_det_avg(5%) ≈ 55–210` epochs — vs `T_req(0.30)≈200–1000`.
- A **1 % anomaly** costs `T_det_avg ≈ 1300–5100` (3–13× imaging) because 1 % ≪ the 30 % imaging
  tolerance — a smaller signal, honestly more expensive, but still within budget in 81 % of cells.
- A **0.5 % anomaly** is below the floor at this `(n=1e4, T_eff≤4096)` budget everywhere.

So the honest claim is **"can't image it, but can detect it"**: the modulator-band wall blocks
imaging beyond `~1.1–1.25×`, but the covariance channel *detects* a `1.25–1.8×` beyond-band anomaly
of ≥1–2 % energy, with the mean channel provably blind.

## 5. Robustness — detection inherits (and exceeds) estimation robustness

Recomputing the profiled Fisher (hence `d'`) at 3 representative cells under E5d-style declared-law
mismatch:

| cell | d'_avg(1 %) | +10 % basis rotation | τ × 2 (true = 2× assumed) |
|---|:--:|:--:|:--:|
| k⁻², k_w=1, claim 1.25 | 8.76 | 8.76 (**+0 %**) | ×0.707 (**−29 %**) |
| flat, k_w=1, claim 1.5 | 5.72 | 5.72 (**+0 %**) | ×0.707 (−29 %) |
| flat, k_w=1, claim 1.8 | 4.82 | 4.82 (**+0 %**) | ×0.707 (−29 %) |

- **Basis rotation (10 %):** the medium covariance `K_w` changes by 14 % (Frobenius) and the Fisher
  eigenvalues *redistribute*, yet the **aggregate detection information `λ_mean = tr(J̄_B)/d_β`
  changes by only 0.04 %**. Detection is far more robust than estimation (E5d showed +1–3 % for the
  estimation ceiling) precisely because it reads the **trace**, not individual modes — the rotation
  is a within-aperture reshuffle that conserves aggregate detectability.
- **τ × 2:** a pure `T_eff` rescaling (the engine works in independent-epoch count; a 2× slower
  medium halves `T_eff`, so `d' ∝ √T_eff` drops 29 %). This is a *diversity* cost at fixed
  wall-clock, identical in kind to E5d's τ finding — not a law-specification fragility.

## 6. Decision

**DETECTION_REGION_FOUND** — the demand-tier row is green and comparator-proof. Where the *supply*
(imaging) map is a tombstone on mode coverage, the *demand* (detection) map lives: a beyond-band
anomaly of ≥1–2 % energy is detectable at `d'≥5` within a 4096-epoch budget across the grid, using a
covariance channel that the mean/fresh-pattern route provably cannot access (specificity theorem),
and with aggregate detectability invariant to declared-law basis rotation. This is the natural
constructive demonstration for the flagship's Fig 1 demand tier — **subject to the R41 ruling and
the parallel hostile prior-art scout** (covariance-based change/anomaly detection and blind-SIM
statistics are the adjacent occupied territory to clear).

### Honest limits (stated plainly)
- The **0.5 % anomaly is undetectable** within 4096 epochs at `n=1e4` everywhere — there is a
  detection floor, not an unlimited demand-tier win.
- The **worst-mode (adversarial) 1 % detection** only lives at the narrow 1.25 k_p claim; the wide
  1.8 k_p claim (240 beyond-band dof) has near-null worst modes that no aggregate softening rescues
  at 1 %.
- These are **CRB/deflection prognoses** (analytic, Gaussian LR), consistent with the MC-validated
  engine (living-region §5); they are a go/no-go, not a built detector.

## Files
`DETECTION_MAP.json` (81-cell grid: `d'` and `T_det` for both directions × 4 energies, specificity),
`DETECTION_ROBUST.json` (3-cell rotation/τ column), `detection_map.py`, `detection_robust.py`.
No git commit.
