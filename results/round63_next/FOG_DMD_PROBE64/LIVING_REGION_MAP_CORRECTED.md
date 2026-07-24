# LIVING-REGION MAP — CORRECTED (physical shot model)

> **This supersedes `LIVING_REGION_MAP.md`.** Regenerated with the physically-correct shot term
> (nonneg photon throughput `|P|·x`) after the shot-model bug named in `CORRECTION_NOTE.md`. Only the
> shot term changed; every other choice is identical. **The manuscript quotes ONLY these corrected
> numbers.** Engine: `living_region_map.py` (corrected), data: `LIVING_REGION_MAP_CORRECTED.json`.

## Verdict (unchanged, hardened)

> **REGION_EMPTY — 0/81 cells reach P1-grade viability**, now with **max witness f_rec = 0.542** and
> **max natural f_rec = 0.667** (both < 0.70). The correction **eliminates even the modest natural-
> scene pocket** the buggy map reported (15 cells → **0**): with the correct shot, no cell clears the
> mode-coverage bar under the strict *or* the relaxed (natural-scene) criterion. σ_f remains inert
> (medium-dominated + truncation-capped). The tombstone is now complete.

## What changed (bug direction confirmed)

The buggy shot (`clamp(Px,1e-12)·mean(Px)/n ≈ 0` for zero-mean signed codes) let `V=C+R` go rank-
deficient (M=128 > db=120) and `V⁻¹` explode in the 8 null directions → **Fisher inflated ~1.9×**.
Correcting it:

| quantity (over 81 cells) | buggy | corrected | direction |
|---|---|---|---|
| max witness f_rec | 0.50 | 0.542 | ~flat, still ≪ 0.70 |
| **max natural f_rec** | **0.79** | **0.667** | **drops below 0.70 → natural pocket gone** |
| mean witness f_rec | 0.258 | 0.215 | drops (kills harden) |
| median NRMSE_nat | 0.210 | 0.392 | **×1.63 (rises)** |
| cells NRMSE_nat ≤ 0.35 | 78 | 46 | fewer pass |
| natural-pocket cells (f_rec_nat ≥ 0.70) | 15 | **0** | eliminated |
| cells P1-viable | 0 | 0 | REGION_EMPTY holds |

**Pocket cell** (σ_f=0.3, k⁻², k_w=1, claim 1.25) — the manuscript's cited cell:

| metric | buggy | corrected |
|---|---|---|
| f_rec (witness) | 0.458 | **0.25** |
| f_rec (natural median) | 0.792 | **0.521** |
| NRMSE_CRB (witness) | 0.066 | **0.099** |
| NRMSE_CRB (natural median) | 0.097 | **0.183** |
| T_req(0.30) (witness) | 197 | **444** |

The natural-median blind acquisition headline (physical shot, blind ≈ CRB) is **T_eff ≈ 2028 banks
≈ 26 s** to reach 0.30 NRMSE (see `CRB_RECONCILIATION.json` / `POCKET_DEMO.md` §6); the witness-based
`T_req` (444) is lower only because the witness carries stronger, fuller annular content.

## Decision (unchanged)

**REGION_EMPTY → law/boundary fallback.** The correction strengthens, not weakens, the conclusion:
no living region exists at flagship grade, and the modest natural-scene residue the buggy map hinted
at does not survive the correct shot. The durable contribution remains the aperture-law + first-
moment-impossibility theorems plus the honest (corrected) information-rate characterization.

---

## Appendix — full 81-cell grid (CORRECTED, physical shot)

f_rec at Fisher SNR≥3, T_eff=4096, n=1e4. P1 requires witness f_rec ≥ 0.70 (met by **none**;
max 0.542).

| sigma_f (eff) | shape | k_w/k_p | claim | f_rec wit | f_rec nat | NRMSE_w | NRMSE_nat | eta_p10 | T_req(0.30) |
|---|---|---|---|---|---|---|---|---|---|
| 0.3 (0.298) | flat | 1 | 1.25 | 0.25 | 0.38 | 0.12 | 0.23 | 0.98 | 712 |
| 0.3 (0.298) | flat | 1 | 1.5 | 0.19 | 0.36 | 0.14 | 0.28 | 0.99 | 860 |
| 0.3 (0.298) | flat | 1 | 1.8 | 0.03 | 0.21 | 0.19 | 0.44 | 0.99 | 1569 |
| 0.3 (0.298) | flat | 2 | 1.25 | 0.48 | 0.29 | 0.17 | 0.35 | 0.99 | 1299 |
| 0.3 (0.298) | flat | 2 | 1.5 | 0.23 | 0.20 | 0.18 | 0.46 | 0.99 | 1560 |
| 0.3 (0.298) | flat | 2 | 1.8 | 0.06 | 0.14 | 0.25 | 0.82 | 0.99 | 2778 |
| 0.3 (0.298) | flat | 4 | 1.25 | 0.29 | 0.10 | 0.28 | 0.56 | 0.98 | 3468 |
| 0.3 (0.298) | flat | 4 | 1.5 | 0.11 | 0.09 | 0.30 | 0.80 | 0.98 | 4209 |
| 0.3 (0.298) | flat | 4 | 1.8 | 0.03 | 0.03 | 0.41 | 1.61 | 0.98 | 7658 |
| 0.3 (0.298) | k^-1 | 1 | 1.25 | 0.27 | 0.52 | 0.12 | 0.22 | 0.98 | 636 |
| 0.3 (0.298) | k^-1 | 1 | 1.5 | 0.19 | 0.36 | 0.14 | 0.28 | 0.99 | 841 |
| 0.3 (0.298) | k^-1 | 1 | 1.8 | 0.11 | 0.22 | 0.20 | 0.49 | 0.99 | 1824 |
| 0.3 (0.298) | k^-1 | 2 | 1.25 | 0.27 | 0.52 | 0.14 | 0.28 | 0.99 | 856 |
| 0.3 (0.298) | k^-1 | 2 | 1.5 | 0.18 | 0.30 | 0.16 | 0.38 | 0.99 | 1152 |
| 0.3 (0.298) | k^-1 | 2 | 1.8 | 0.13 | 0.17 | 0.24 | 0.75 | 0.99 | 2599 |
| 0.3 (0.298) | k^-1 | 4 | 1.25 | 0.31 | 0.38 | 0.16 | 0.34 | 0.98 | 1231 |
| 0.3 (0.298) | k^-1 | 4 | 1.5 | 0.18 | 0.23 | 0.19 | 0.47 | 0.99 | 1699 |
| 0.3 (0.298) | k^-1 | 4 | 1.8 | 0.14 | 0.11 | 0.30 | 1.02 | 0.99 | 4062 |
| 0.3 (0.298) | k^-2 | 1 | 1.25 | 0.25 | 0.52 | 0.10 | 0.18 | 0.98 | 444 |
| 0.3 (0.298) | k^-2 | 1 | 1.5 | 0.26 | 0.41 | 0.13 | 0.28 | 0.98 | 768 |
| 0.3 (0.298) | k^-2 | 1 | 1.8 | 0.34 | 0.22 | 0.23 | 0.61 | 0.99 | 2520 |
| 0.3 (0.298) | k^-2 | 2 | 1.25 | 0.27 | 0.54 | 0.10 | 0.20 | 0.98 | 501 |
| 0.3 (0.298) | k^-2 | 2 | 1.5 | 0.25 | 0.38 | 0.14 | 0.31 | 0.99 | 882 |
| 0.3 (0.298) | k^-2 | 2 | 1.8 | 0.31 | 0.19 | 0.26 | 0.80 | 0.99 | 3009 |
| 0.3 (0.298) | k^-2 | 4 | 1.25 | 0.25 | 0.52 | 0.11 | 0.21 | 0.98 | 554 |
| 0.3 (0.298) | k^-2 | 4 | 1.5 | 0.24 | 0.35 | 0.15 | 0.34 | 0.99 | 989 |
| 0.3 (0.298) | k^-2 | 4 | 1.8 | 0.32 | 0.17 | 0.28 | 0.88 | 0.99 | 3498 |
| 0.6 (0.503) | flat | 1 | 1.25 | 0.27 | 0.58 | 0.11 | 0.19 | 0.98 | 557 |
| 0.6 (0.503) | flat | 1 | 1.5 | 0.20 | 0.42 | 0.12 | 0.24 | 0.99 | 674 |
| 0.6 (0.503) | flat | 1 | 1.8 | 0.03 | 0.32 | 0.16 | 0.38 | 0.99 | 1224 |
| 0.6 (0.503) | flat | 2 | 1.25 | 0.50 | 0.48 | 0.14 | 0.27 | 0.99 | 837 |
| 0.6 (0.503) | flat | 2 | 1.5 | 0.23 | 0.35 | 0.15 | 0.34 | 0.99 | 1000 |
| 0.6 (0.503) | flat | 2 | 1.8 | 0.07 | 0.17 | 0.20 | 0.58 | 0.99 | 1755 |
| 0.6 (0.503) | flat | 4 | 1.25 | 0.46 | 0.23 | 0.19 | 0.38 | 0.99 | 1578 |
| 0.6 (0.503) | flat | 4 | 1.5 | 0.17 | 0.16 | 0.20 | 0.52 | 0.99 | 1904 |
| 0.6 (0.503) | flat | 4 | 1.8 | 0.04 | 0.11 | 0.27 | 0.95 | 0.98 | 3406 |
| 0.6 (0.503) | k^-1 | 1 | 1.25 | 0.31 | 0.48 | 0.10 | 0.18 | 0.98 | 489 |
| 0.6 (0.503) | k^-1 | 1 | 1.5 | 0.18 | 0.38 | 0.12 | 0.23 | 0.98 | 636 |
| 0.6 (0.503) | k^-1 | 1 | 1.8 | 0.06 | 0.30 | 0.17 | 0.39 | 0.99 | 1314 |
| 0.6 (0.503) | k^-1 | 2 | 1.25 | 0.33 | 0.52 | 0.12 | 0.22 | 0.98 | 597 |
| 0.6 (0.503) | k^-1 | 2 | 1.5 | 0.19 | 0.41 | 0.13 | 0.29 | 0.99 | 785 |
| 0.6 (0.503) | k^-1 | 2 | 1.8 | 0.10 | 0.20 | 0.19 | 0.54 | 0.99 | 1658 |
| 0.6 (0.503) | k^-1 | 4 | 1.25 | 0.25 | 0.42 | 0.13 | 0.26 | 0.99 | 753 |
| 0.6 (0.503) | k^-1 | 4 | 1.5 | 0.22 | 0.33 | 0.15 | 0.34 | 0.99 | 1007 |
| 0.6 (0.503) | k^-1 | 4 | 1.8 | 0.15 | 0.23 | 0.22 | 0.67 | 0.99 | 2231 |
| 0.6 (0.503) | k^-2 | 1 | 1.25 | 0.27 | 0.67 | 0.08 | 0.15 | 0.98 | 318 |
| 0.6 (0.503) | k^-2 | 1 | 1.5 | 0.21 | 0.50 | 0.11 | 0.22 | 0.99 | 526 |
| 0.6 (0.503) | k^-2 | 1 | 1.8 | 0.24 | 0.30 | 0.18 | 0.44 | 0.99 | 1524 |
| 0.6 (0.503) | k^-2 | 2 | 1.25 | 0.19 | 0.54 | 0.09 | 0.16 | 0.98 | 346 |
| 0.6 (0.503) | k^-2 | 2 | 1.5 | 0.24 | 0.48 | 0.11 | 0.24 | 0.99 | 579 |
| 0.6 (0.503) | k^-2 | 2 | 1.8 | 0.23 | 0.28 | 0.20 | 0.55 | 0.99 | 1724 |
| 0.6 (0.503) | k^-2 | 4 | 1.25 | 0.21 | 0.54 | 0.09 | 0.17 | 0.98 | 370 |
| 0.6 (0.503) | k^-2 | 4 | 1.5 | 0.21 | 0.43 | 0.12 | 0.25 | 0.99 | 625 |
| 0.6 (0.503) | k^-2 | 4 | 1.8 | 0.27 | 0.25 | 0.20 | 0.58 | 0.99 | 1912 |
| 1.0 (0.696) | flat | 1 | 1.25 | 0.29 | 0.60 | 0.10 | 0.17 | 0.98 | 506 |
| 1.0 (0.696) | flat | 1 | 1.5 | 0.20 | 0.52 | 0.12 | 0.22 | 0.99 | 612 |
| 1.0 (0.696) | flat | 1 | 1.8 | 0.04 | 0.36 | 0.16 | 0.36 | 0.99 | 1110 |
| 1.0 (0.696) | flat | 2 | 1.25 | 0.54 | 0.48 | 0.12 | 0.24 | 0.99 | 700 |
| 1.0 (0.696) | flat | 2 | 1.5 | 0.21 | 0.38 | 0.14 | 0.30 | 0.99 | 835 |
| 1.0 (0.696) | flat | 2 | 1.8 | 0.06 | 0.23 | 0.18 | 0.49 | 0.99 | 1456 |
| 1.0 (0.696) | flat | 4 | 1.25 | 0.48 | 0.31 | 0.16 | 0.32 | 0.99 | 1131 |
| 1.0 (0.696) | flat | 4 | 1.5 | 0.21 | 0.24 | 0.17 | 0.42 | 0.99 | 1362 |
| 1.0 (0.696) | flat | 4 | 1.8 | 0.06 | 0.15 | 0.23 | 0.74 | 0.98 | 2416 |
| 1.0 (0.696) | k^-1 | 1 | 1.25 | 0.33 | 0.62 | 0.10 | 0.17 | 0.98 | 441 |
| 1.0 (0.696) | k^-1 | 1 | 1.5 | 0.18 | 0.51 | 0.11 | 0.21 | 0.98 | 569 |
| 1.0 (0.696) | k^-1 | 1 | 1.8 | 0.04 | 0.29 | 0.16 | 0.36 | 0.99 | 1152 |
| 1.0 (0.696) | k^-1 | 2 | 1.25 | 0.27 | 0.54 | 0.11 | 0.20 | 0.98 | 516 |
| 1.0 (0.696) | k^-1 | 2 | 1.5 | 0.18 | 0.36 | 0.12 | 0.26 | 0.99 | 671 |
| 1.0 (0.696) | k^-1 | 2 | 1.8 | 0.08 | 0.27 | 0.17 | 0.46 | 0.99 | 1377 |
| 1.0 (0.696) | k^-1 | 4 | 1.25 | 0.29 | 0.52 | 0.12 | 0.22 | 0.99 | 615 |
| 1.0 (0.696) | k^-1 | 4 | 1.5 | 0.17 | 0.44 | 0.13 | 0.29 | 0.99 | 811 |
| 1.0 (0.696) | k^-1 | 4 | 1.8 | 0.11 | 0.23 | 0.20 | 0.55 | 0.99 | 1732 |
| 1.0 (0.696) | k^-2 | 1 | 1.25 | 0.23 | 0.65 | 0.08 | 0.13 | 0.98 | 277 |
| 1.0 (0.696) | k^-2 | 1 | 1.5 | 0.21 | 0.49 | 0.10 | 0.19 | 0.98 | 451 |
| 1.0 (0.696) | k^-2 | 1 | 1.8 | 0.16 | 0.32 | 0.17 | 0.38 | 0.98 | 1235 |
| 1.0 (0.696) | k^-2 | 2 | 1.25 | 0.23 | 0.58 | 0.08 | 0.15 | 0.98 | 297 |
| 1.0 (0.696) | k^-2 | 2 | 1.5 | 0.19 | 0.51 | 0.10 | 0.22 | 0.98 | 486 |
| 1.0 (0.696) | k^-2 | 2 | 1.8 | 0.20 | 0.29 | 0.17 | 0.46 | 0.99 | 1362 |
| 1.0 (0.696) | k^-2 | 4 | 1.25 | 0.21 | 0.56 | 0.08 | 0.15 | 0.98 | 312 |
| 1.0 (0.696) | k^-2 | 4 | 1.5 | 0.19 | 0.50 | 0.11 | 0.22 | 0.99 | 516 |
| 1.0 (0.696) | k^-2 | 4 | 1.8 | 0.20 | 0.29 | 0.18 | 0.49 | 0.99 | 1478 |
