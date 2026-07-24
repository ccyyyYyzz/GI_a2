# ⚠️ CORRECTION NOTE — shot-model bug in the Fisher/CRB engines (read first)

**Date:** 2026-07-24 · **Scope:** the profiled-Fisher CRB numbers in `p1_fisher.py` and
`living_region_map.py` (and any artifact quoting them). **No verdict flips. All kill verdicts
HARDEN; all positive headlines SHRINK. The manuscript quotes ONLY the `_CORRECTED` artifacts.**

## The bug

Both Fisher engines set the per-bucket shot term as

```python
R_shot = clamp(P @ x, 1e-12) * mean(P @ x) / n        # WRONG for signed codes
```

For the **zero-mean signed codes** used here (`P` band-limited, DC excluded), the signed mean
bucket `P @ x ≈ 0`, so `R_shot ≈ 1e-3` — **~26,000× too small** (measured: 1.09e-3 vs the correct
28.7). The one-epoch covariance `C(x) = P·diag(x)·K_w·diag(x)·Pᵀ` is **rank-deficient (rank 120 of
M=128, because M > db=120)**; with `R_shot ≈ 0` the 8 null eigenvalues of `V = C + R` collapse to
the `1e-9` numerical floor, so `V⁻¹` **explodes in those 8 directions** and inflates the covariance
Fisher `½·tr(V⁻¹ V_u V⁻¹ V_v)`. Net effect: **Fisher over-stated ≈ 1.9×** (natural median), i.e. the
CRB was **under-stated** and f_rec **over-stated**.

The physically correct shot is the **nonneg photon throughput** (each signed code is a complementary
pair of nonneg exposures; the detected photons scale with `|P|·x`, not the signed difference):

```python
flux = abs(P) @ x
R_shot = flux * mean(flux) / n                        # CORRECT (medium/shot = 14.7, "1e4 clean")
```

This is the shot model `pocket_demo.py` and the gates already used, so **only `p1_fisher.py` and
`living_region_map.py` were affected.** The fix is a one-line change in each (applied in place;
corrected artifacts written to `_CORRECTED` filenames; originals preserved for provenance).

## Direction of the correction (kills harden, positives shrink)

| quantity | direction under correction |
|---|---|
| Fisher information | ↓ ~1.9× (was inflated) |
| f_rec (mode coverage) | **↓** — every P1/REGION kill **hardens** |
| NRMSE_CRB, T_req | **↑ ~1.6–1.9×** — every positive headline **shrinks** |

## Affected committed artifacts and the corrected replacements

| artifact | status | corrected replacement |
|---|---|---|
| `P1_results.json` | SUPERSEDED | `P1_results_CORRECTED.json` |
| `LIVING_REGION_MAP.json` | SUPERSEDED | `LIVING_REGION_MAP_CORRECTED.json` |
| `LIVING_REGION_MAP.md` | SUPERSEDED (banner added) | `LIVING_REGION_MAP_CORRECTED.md` |
| `RISK_GATES_REPORT.md` (quotes P1 CRB/f_rec) | banner added | corrected numbers below |
| `DETECTION_MAP.*` (detection agent) | affected — **that agent regenerates its own** | — |
| `POCKET_DEMO.md/.json`, `G1/G2/G3_results.json`, gates | **unaffected** (already physical/empirical shot) | — |

## No verdict flips — corrected numbers

**P1 (prelaunch Fisher kill):**

| | buggy | corrected |
|---|---|---|
| witness f_rec (bar ≥0.70) | 0.050 | **0.033** |
| witness NRMSE_CRB | 0.140 | 0.185 |
| natural-median NRMSE_CRB (bar ≤0.35) | 0.312 (pass) | **0.443 (FAIL)** |
| failing bars | {f_rec} | **{f_rec, nrmse_nat}** |
| **verdict** | P1_FAIL | **P1_FAIL (hardened)** |

**Living-region map (81 cells):**

| | buggy | corrected |
|---|---|---|
| P1-viable cells | 0/81 | **0/81** |
| max witness f_rec | 0.50 | 0.542 (still ≪ 0.70) |
| max natural f_rec | 0.79 | **0.667 (< 0.70)** |
| natural-scene pocket (f_rec_nat ≥ 0.70) | 15 cells | **0 cells (eliminated)** |
| median NRMSE_nat | 0.210 | 0.392 (×1.63) |
| pocket cell T_req(0.30), witness | 197 | **444** |
| **verdict** | REGION_EMPTY | **REGION_EMPTY (hardened; natural residue gone)** |

**Pocket demo** (`pocket_demo.py`, already correct shot): **POCKET_ACHIEVED** on the efficiency
criterion — the moment estimator is CRB-efficient at all T; the E8 plateau was the finite-T CRB, not
an estimator gap. **Unchanged by this correction.**

## Corrected manuscript headline

Blind natural-median beyond-band NRMSE reaches **0.30 at T_eff ≈ 2028 banks ≈ 26 s** (20 kHz,
M=128 complementary; one medium realization per bank ⇒ T_eff = banks, no τ penalty); synthetic
witness ≈ **112 banks ≈ 1.4 s**. (`CRB_RECONCILIATION.json`, `POCKET_DEMO.md` §5–6.) The buggy
`T_req ≈ 197` was the *witness* value *with the shot bug* — doubly optimistic.

**Bottom line:** the correction only strengthens the program's conclusion — REGION_EMPTY and P1_FAIL
harden, the modest natural pocket disappears, and the honest pocket headline lengthens to ~26 s.
Nothing flips. Decision stands: **law/boundary fallback.**
