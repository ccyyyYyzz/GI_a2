> 🛑 **SUPERSEDED — do not quote these numbers.** The CRB/f_rec here use a buggy shot term
> (`clamp(Px)` on zero-mean signed codes ≈ 0 → Fisher inflated ~1.9×; see `CORRECTION_NOTE.md`).
> Use **`LIVING_REGION_MAP_CORRECTED.md`**. Verdict (REGION_EMPTY) is unchanged and hardened.
> Kept in place for provenance only.

# LIVING-REGION MAP — where (if anywhere) is the beyond-band method P1-viable?

**Date:** 2026-07-24 · **Engine:** `living_region_map.py` (reuses the P1 exact profiled-Fisher
analytics; NO reconstruction) · **Data:** `LIVING_REGION_MAP.json`. Ran locally on the RTX 4060
(~2 min). Fixed: `N=4096`, `M=128`, `k_p=5`, `n=1e4`, `T_eff≤4096`.

## Verdict

> **REGION_EMPTY — 0 / 81 cells reach P1-grade viability** under the ruling's strict (synthetic-
> witness) `f_rec ≥ 0.70` bar. The tombstone map. **σ_f is NOT the lever** the naive `σ_f⁴` scaling
> suggested: at `n=1e4` the channel is medium-dominated (eq 4.12 saturation), so f_rec is nearly
> flat across `σ_f`, and the bounded mean-one field caps `σ_f=1.0` at an effective contrast of only
> **0.70** anyway. The binding wall is **mode coverage (f_rec)**, not aggregate error, not the
> profiling tax, and not shot noise.
>
> **One honest residue** (not a flagship): a **≈1.25× extension on natural scenes** with a `k⁻²`
> (realistic decaying) medium spectrum is information-feasible — `f_rec(natural)` reaches **0.79**,
> `NRMSE_CRB ≈ 0.10`, `T_req(0.30) ≈ 200` epochs. It fails the worst-case witness bar, is a modest
> 1.25× (not the 1.8× flagship), and sits in occupied prior art. **Decision: law/boundary fallback,
> not R40.**

---

## 1. The bounded-field truncation deflates the σ_f lever (MC-measured)

The bounded mean-one positive field (`w = clip(1+g)`) cannot realize high nominal contrast — the
clip caps it:

| σ_f nominal | bounds | **σ_f effective** | clip fraction |
|---|---|---|---|
| 0.3 | [0.2, 1.8] | 0.298 | 0.8% |
| 0.6 | [0.2, 1.8] | 0.503 | 18.2% |
| 1.0 | [0.05, 1.95] | **0.696** | 34.2% |

So the "`σ_f=1.0` fully-developed-speckle big lever" realizes only ≈0.70 effective contrast even
after loosening the bounds to `[0.05, 1.95]`. Combined with medium-dominated saturation (§3), the
lever is inert here.

## 2. The binding wall is f_rec, not error — ceiling per claim

`f_rec` = fraction of the claimed beyond-band annulus at Fisher SNR ≥ 3 (the mode-coverage bar).
Best cell in the whole grid for each claim (across all σ_f / shape / k_w):

| claim | max f_rec (witness) | max f_rec (natural) | NRMSE_CRB witness | T_req(0.30) |
|---|---|---|---|---|
| 1.25 k_p | **0.50** | **0.79** | 0.07–0.10 | 197–482 |
| 1.5 k_p | 0.34 | 0.71 | 0.08–0.11 | 309–576 |
| 1.8 k_p | 0.22 | 0.44 | 0.13–0.15 | 744–995 |

Reading: **the NRMSE_CRB bars pass in 78/81 cells** (witness ≤ 0.25 AND natural ≤ 0.35), and
**T_req(0.30) ≤ 4096 in all 81 cells** — the *aggregate* beyond-band error target is reachable.
But **witness f_rec never reaches 0.70** (max 0.50). The support aperture is fine; a *minority of
strong near-k_p modes* carries the aggregate error while the rest of even the 1.25× annulus is
sub-SNR-3. This is the exact "support ≠ Fisher-weighted aperture" gap (§2.3), now mapped across
physical space and found everywhere.

## 3. Why the axes don't rescue it

- **σ_f (the "big lever"):** f_rec at (flat, k_w=1, claim 1.25) is 0.44 / 0.40 / 0.35 for
  σ_f = 0.3 / 0.6 / 1.0 — **flat-to-slightly-decreasing**. At `n=1e4` we are deep in the medium-
  dominated regime (`nσ_f² ≫ 1`): extra contrast does not buy blind covariance information per
  independent state (ruling eq 4.12). The lever is inert; the currency is independent medium states
  (T_eff), not contrast or photons.
- **Spectrum shape:** `k⁻²` (realistic decay) is the **best** for aggregate error — it concentrates
  medium power in low modes near the aperture, giving NRMSE 0.07 and the lowest T_req (~200) — yet
  witness f_rec stays 0.44–0.48. (Note the reversal vs G3: at the *aggressive* geometry decaying
  killed recovery; at the *frozen* geometry with a *narrow* claim, decaying helps.)
- **Medium band k_w:** widening (1→4×k_p) nudges witness f_rec up modestly (0.44→0.50 at claim
  1.25) but dilutes per-mode power (total variance fixed over more modes), leaving NRMSE ~flat. An
  honest wash — never crosses 0.70.
- **Profiling tax:** `η_prof` p10 ≈ 0.97–0.99 everywhere — the in-band mean channel pins the
  nuisance so profiling costs almost nothing. Profiling is *not* the problem.

## 4. The natural-scene modest residue (15 cells, all sub-P1)

Relaxing f_rec to natural scenes (energy concentrated in the recoverable near-modes) and requiring
the NRMSE bars: **15 cells** clear `f_rec(natural) ≥ 0.70`, all at **claim ≤ 1.5 k_p** (mostly
1.25) and **all with `k⁻²` spectrum**. Most defensible:

> **σ_f=0.3, k⁻² spectrum, k_w=1·k_p, claim 1.25·k_p:** f_rec(nat) 0.79, f_rec(witness) 0.46,
> NRMSE_CRB 0.07 (witness) / 0.10 (natural), **T_req(0.30) ≈ 197 epochs** (≈ 2.5 s at 20 kHz DMD).

This is a genuine, small, honest effect — a ~1.25× covariance-aperture extension on natural imagery
at very modest integration. It is **not** a flagship: it fails the worst-case witness f_rec bar, is
1.25× (not the ruling's 1.8× claim), and covariance superresolution is occupied prior art
(Chaigne/Idier/Mudry). Per the ruling kill tree it does not revive the mainline.

## 5. MC cross-check of the Fisher engine (ruling §7.4 independence)

| cell | analytic vs MC (Gaussian) | analytic vs MC (bounded) |
|---|---|---|
| σ_f=0.3, flat, k_w=1 | 0.122 | 0.122 |
| σ_f=1.0, flat, k_w=1 | 0.122 | 0.312 |
| σ_f=0.6, k⁻², k_w=2 | 0.087 | 0.222 |

The **Gaussian** rel-error is a constant ≈0.09–0.12 **independent of σ_f** — it is the MC-estimation
floor at `T_mc=4000` (a 128×128 covariance from 4000 epochs), not an engine error; a wrong analytic
`C(x)` would drift with σ_f. **Engine validated.** The **bounded** rel-error grows 0.12 → 0.31 as
σ_f 0.3 → 1.0, independently quantifying the truncation deviation of the realized speckle covariance
from the nominal band-limited law at high contrast (consistent with the clip diagnostic in §1).

## 6. Decision

**REGION_EMPTY → law/boundary fallback, not R40.** No point in the swept physical space clears the
ruling's P1 flagship bar; the beyond-band *method* is Fisher-rate-limited to a thin ≈1.1–1.25× rim
regardless of σ_f, spectrum, band width, or claim shell. The **aperture-law + first-moment-
impossibility theorems** (proven, and P0-clean: no null, wall exact) are the durable contribution →
archive as the materials-bank / law-boundary note. If the operator wants a *small* positive result,
the §4 natural-scene 1.25× `k⁻²` pocket is the only honest candidate — a modest covariance-aperture
demonstration, explicitly sub-P1 and prior-art-adjacent, not the flagship.

---

## Appendix — full 81-cell grid

f_rec at Fisher SNR≥3, T_eff=4096, n=1e4. `f_rec wit` = synthetic witness (the P1 bar);
`f_rec nat` = natural-scene median. P1 requires witness f_rec ≥ 0.70 (met by **none**).

| sigma_f (eff) | shape | k_w/k_p | claim | f_rec wit | f_rec nat | NRMSE_w | NRMSE_nat | eta_p10 | T_req(0.30) |
|---|---|---|---|---|---|---|---|---|---|
| 0.3 (0.298) | flat | 1 | 1.25 | 0.44 | 0.67 | 0.10 | 0.13 | 0.98 | 409 |
| 0.3 (0.298) | flat | 1 | 1.5 | 0.31 | 0.53 | 0.10 | 0.18 | 0.98 | 494 |
| 0.3 (0.298) | flat | 1 | 1.8 | 0.05 | 0.41 | 0.14 | 0.31 | 0.99 | 891 |
| 0.3 (0.298) | flat | 2 | 1.25 | 0.40 | 0.56 | 0.10 | 0.16 | 0.99 | 457 |
| 0.3 (0.298) | flat | 2 | 1.5 | 0.34 | 0.56 | 0.11 | 0.21 | 0.99 | 542 |
| 0.3 (0.298) | flat | 2 | 1.8 | 0.09 | 0.42 | 0.14 | 0.34 | 0.99 | 926 |
| 0.3 (0.298) | flat | 4 | 1.25 | 0.50 | 0.58 | 0.10 | 0.16 | 0.99 | 482 |
| 0.3 (0.298) | flat | 4 | 1.5 | 0.34 | 0.51 | 0.11 | 0.21 | 0.99 | 576 |
| 0.3 (0.298) | flat | 4 | 1.8 | 0.13 | 0.34 | 0.15 | 0.35 | 0.99 | 993 |
| 0.3 (0.298) | k^-1 | 1 | 1.25 | 0.44 | 0.62 | 0.09 | 0.12 | 0.98 | 347 |
| 0.3 (0.298) | k^-1 | 1 | 1.5 | 0.28 | 0.58 | 0.10 | 0.17 | 0.98 | 442 |
| 0.3 (0.298) | k^-1 | 1 | 1.8 | 0.08 | 0.42 | 0.14 | 0.30 | 0.99 | 852 |
| 0.3 (0.298) | k^-1 | 2 | 1.25 | 0.31 | 0.73 | 0.09 | 0.14 | 0.98 | 366 |
| 0.3 (0.298) | k^-1 | 2 | 1.5 | 0.19 | 0.55 | 0.10 | 0.19 | 0.99 | 464 |
| 0.3 (0.298) | k^-1 | 2 | 1.8 | 0.07 | 0.40 | 0.14 | 0.33 | 0.99 | 884 |
| 0.3 (0.298) | k^-1 | 4 | 1.25 | 0.38 | 0.73 | 0.09 | 0.14 | 0.98 | 374 |
| 0.3 (0.298) | k^-1 | 4 | 1.5 | 0.19 | 0.57 | 0.10 | 0.19 | 0.99 | 476 |
| 0.3 (0.298) | k^-1 | 4 | 1.8 | 0.07 | 0.39 | 0.14 | 0.34 | 0.99 | 918 |
| 0.3 (0.298) | k^-2 | 1 | 1.25 | 0.46 | 0.79 | 0.07 | 0.10 | 0.99 | 197 |
| 0.3 (0.298) | k^-2 | 1 | 1.5 | 0.33 | 0.70 | 0.08 | 0.14 | 0.98 | 309 |
| 0.3 (0.298) | k^-2 | 1 | 1.8 | 0.22 | 0.42 | 0.13 | 0.28 | 0.99 | 744 |
| 0.3 (0.298) | k^-2 | 2 | 1.25 | 0.38 | 0.71 | 0.07 | 0.10 | 0.99 | 202 |
| 0.3 (0.298) | k^-2 | 2 | 1.5 | 0.32 | 0.58 | 0.08 | 0.15 | 0.99 | 316 |
| 0.3 (0.298) | k^-2 | 2 | 1.8 | 0.15 | 0.42 | 0.13 | 0.30 | 0.99 | 764 |
| 0.3 (0.298) | k^-2 | 4 | 1.25 | 0.38 | 0.67 | 0.07 | 0.11 | 0.99 | 203 |
| 0.3 (0.298) | k^-2 | 4 | 1.5 | 0.28 | 0.59 | 0.08 | 0.15 | 0.99 | 320 |
| 0.3 (0.298) | k^-2 | 4 | 1.8 | 0.15 | 0.44 | 0.13 | 0.31 | 0.99 | 777 |
| 0.6 (0.503) | flat | 1 | 1.25 | 0.40 | 0.65 | 0.10 | 0.14 | 0.98 | 410 |
| 0.6 (0.503) | flat | 1 | 1.5 | 0.28 | 0.53 | 0.10 | 0.19 | 0.98 | 495 |
| 0.6 (0.503) | flat | 1 | 1.8 | 0.04 | 0.39 | 0.14 | 0.31 | 0.98 | 894 |
| 0.6 (0.503) | flat | 2 | 1.25 | 0.40 | 0.58 | 0.10 | 0.16 | 0.99 | 458 |
| 0.6 (0.503) | flat | 2 | 1.5 | 0.32 | 0.57 | 0.11 | 0.21 | 0.99 | 543 |
| 0.6 (0.503) | flat | 2 | 1.8 | 0.08 | 0.40 | 0.14 | 0.34 | 0.99 | 927 |
| 0.6 (0.503) | flat | 4 | 1.25 | 0.50 | 0.58 | 0.10 | 0.16 | 0.99 | 482 |
| 0.6 (0.503) | flat | 4 | 1.5 | 0.33 | 0.53 | 0.11 | 0.21 | 0.99 | 576 |
| 0.6 (0.503) | flat | 4 | 1.8 | 0.12 | 0.35 | 0.15 | 0.35 | 0.99 | 994 |
| 0.6 (0.503) | k^-1 | 1 | 1.25 | 0.44 | 0.65 | 0.09 | 0.13 | 0.98 | 349 |
| 0.6 (0.503) | k^-1 | 1 | 1.5 | 0.26 | 0.56 | 0.10 | 0.17 | 0.98 | 443 |
| 0.6 (0.503) | k^-1 | 1 | 1.8 | 0.05 | 0.44 | 0.14 | 0.30 | 0.98 | 854 |
| 0.6 (0.503) | k^-1 | 2 | 1.25 | 0.31 | 0.77 | 0.09 | 0.14 | 0.98 | 366 |
| 0.6 (0.503) | k^-1 | 2 | 1.5 | 0.17 | 0.56 | 0.10 | 0.19 | 0.99 | 464 |
| 0.6 (0.503) | k^-1 | 2 | 1.8 | 0.07 | 0.41 | 0.14 | 0.33 | 0.99 | 885 |
| 0.6 (0.503) | k^-1 | 4 | 1.25 | 0.35 | 0.75 | 0.09 | 0.14 | 0.98 | 374 |
| 0.6 (0.503) | k^-1 | 4 | 1.5 | 0.18 | 0.57 | 0.10 | 0.19 | 0.99 | 476 |
| 0.6 (0.503) | k^-1 | 4 | 1.8 | 0.07 | 0.40 | 0.14 | 0.34 | 0.99 | 919 |
| 0.6 (0.503) | k^-2 | 1 | 1.25 | 0.48 | 0.77 | 0.07 | 0.10 | 0.98 | 198 |
| 0.6 (0.503) | k^-2 | 1 | 1.5 | 0.34 | 0.67 | 0.08 | 0.14 | 0.98 | 310 |
| 0.6 (0.503) | k^-2 | 1 | 1.8 | 0.15 | 0.43 | 0.13 | 0.28 | 0.98 | 746 |
| 0.6 (0.503) | k^-2 | 2 | 1.25 | 0.42 | 0.71 | 0.07 | 0.11 | 0.99 | 202 |
| 0.6 (0.503) | k^-2 | 2 | 1.5 | 0.28 | 0.58 | 0.08 | 0.15 | 0.99 | 317 |
| 0.6 (0.503) | k^-2 | 2 | 1.8 | 0.13 | 0.42 | 0.13 | 0.30 | 0.98 | 765 |
| 0.6 (0.503) | k^-2 | 4 | 1.25 | 0.44 | 0.71 | 0.07 | 0.11 | 0.99 | 204 |
| 0.6 (0.503) | k^-2 | 4 | 1.5 | 0.27 | 0.55 | 0.08 | 0.16 | 0.99 | 320 |
| 0.6 (0.503) | k^-2 | 4 | 1.8 | 0.14 | 0.44 | 0.13 | 0.31 | 0.99 | 778 |
| 1.0 (0.696) | flat | 1 | 1.25 | 0.35 | 0.67 | 0.10 | 0.14 | 0.98 | 411 |
| 1.0 (0.696) | flat | 1 | 1.5 | 0.26 | 0.51 | 0.10 | 0.19 | 0.98 | 497 |
| 1.0 (0.696) | flat | 1 | 1.8 | 0.03 | 0.40 | 0.14 | 0.31 | 0.98 | 897 |
| 1.0 (0.696) | flat | 2 | 1.25 | 0.40 | 0.58 | 0.10 | 0.16 | 0.99 | 458 |
| 1.0 (0.696) | flat | 2 | 1.5 | 0.31 | 0.60 | 0.11 | 0.21 | 0.99 | 543 |
| 1.0 (0.696) | flat | 2 | 1.8 | 0.07 | 0.37 | 0.14 | 0.34 | 0.99 | 928 |
| 1.0 (0.696) | flat | 4 | 1.25 | 0.50 | 0.58 | 0.10 | 0.16 | 0.99 | 482 |
| 1.0 (0.696) | flat | 4 | 1.5 | 0.33 | 0.52 | 0.11 | 0.21 | 0.99 | 576 |
| 1.0 (0.696) | flat | 4 | 1.8 | 0.11 | 0.33 | 0.15 | 0.35 | 0.99 | 995 |
| 1.0 (0.696) | k^-1 | 1 | 1.25 | 0.44 | 0.67 | 0.09 | 0.13 | 0.97 | 350 |
| 1.0 (0.696) | k^-1 | 1 | 1.5 | 0.25 | 0.57 | 0.10 | 0.17 | 0.98 | 445 |
| 1.0 (0.696) | k^-1 | 1 | 1.8 | 0.04 | 0.42 | 0.14 | 0.30 | 0.98 | 857 |
| 1.0 (0.696) | k^-1 | 2 | 1.25 | 0.31 | 0.75 | 0.09 | 0.14 | 0.98 | 366 |
| 1.0 (0.696) | k^-1 | 2 | 1.5 | 0.16 | 0.57 | 0.10 | 0.19 | 0.99 | 465 |
| 1.0 (0.696) | k^-1 | 2 | 1.8 | 0.06 | 0.40 | 0.14 | 0.33 | 0.99 | 887 |
| 1.0 (0.696) | k^-1 | 4 | 1.25 | 0.33 | 0.75 | 0.09 | 0.14 | 0.98 | 374 |
| 1.0 (0.696) | k^-1 | 4 | 1.5 | 0.18 | 0.57 | 0.10 | 0.19 | 0.99 | 476 |
| 1.0 (0.696) | k^-1 | 4 | 1.8 | 0.07 | 0.40 | 0.14 | 0.34 | 0.99 | 919 |
| 1.0 (0.696) | k^-2 | 1 | 1.25 | 0.44 | 0.79 | 0.07 | 0.10 | 0.98 | 199 |
| 1.0 (0.696) | k^-2 | 1 | 1.5 | 0.33 | 0.71 | 0.08 | 0.14 | 0.98 | 310 |
| 1.0 (0.696) | k^-2 | 1 | 1.8 | 0.12 | 0.43 | 0.13 | 0.28 | 0.98 | 747 |
| 1.0 (0.696) | k^-2 | 2 | 1.25 | 0.40 | 0.69 | 0.07 | 0.11 | 0.98 | 203 |
| 1.0 (0.696) | k^-2 | 2 | 1.5 | 0.26 | 0.59 | 0.08 | 0.16 | 0.98 | 318 |
| 1.0 (0.696) | k^-2 | 2 | 1.8 | 0.12 | 0.42 | 0.13 | 0.30 | 0.98 | 766 |
| 1.0 (0.696) | k^-2 | 4 | 1.25 | 0.38 | 0.71 | 0.07 | 0.11 | 0.99 | 204 |
| 1.0 (0.696) | k^-2 | 4 | 1.5 | 0.24 | 0.58 | 0.08 | 0.16 | 0.99 | 321 |
| 1.0 (0.696) | k^-2 | 4 | 1.8 | 0.13 | 0.44 | 0.13 | 0.31 | 0.99 | 779 |
