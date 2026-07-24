# DETECTION LIVING-REGION MAP — CORRECTED (physical shot term)

**Date:** 2026-07-24 · **Supersedes** `DETECTION_PROGNOSIS.md` / `DETECTION_MAP.json` (buggy shot).
**Engine:** `detection_map.py` (shot fixed) · **Data:** `DETECTION_MAP_CORRECTED.json`. Same 81-cell
grid, N=4096, M=128, n=1e4.

## The correction

The shared Fisher engine used `R_shot = clamp(Px,1e-12)·mean(Px)/n` with the **signed** bucket mean
`Px` — near-zero for zero-mean signed codes, which nulled the shot term and inflated the Fisher
~1.9×. The **physically correct** shot uses the nonneg **photon throughput** `|P|·x` (signed codes
are realized as complementary nonneg pairs, so the light on the detector is `|P|·(w⊙x)`):
`R_shot = diag(|P|·x · mean(|P|·x)/n)`. At the pocket cell the realized **medium/shot ratio is
≈11.8** (shot is a real ~8% of the diagonal, not negligible). Net effect: analytic `d'` ≈ **1.38×
optimistic**, `T_det` ≈ **1.9× optimistic** in the old map. **All numbers below are corrected.**

## Acquisition model (reconciliation)

Bank acquisition: a quasi-static diffuser is held for one code bank (256 exposures @ 20 kHz =
**12.8 ms/bank**), then **stepped** to an independent realization for the next bank. So
`T_eff = number of independent banks` (**no τ / OU wall-clock penalty**), and detection time is
`T_det · 12.8 ms`. All `T_det` below are in **banks**; seconds shown alongside.

## Verdict

> **DETECTION_REGION_FOUND — survives the correction, region shrinks.** A **2 % beyond-band anomaly
> is detectable (energy-spread, `d'≥5`) within 4096 banks in 95 % of the 81 cells** (best **≈467
> banks = 6.0 s**); a **1 % anomaly only at the best cells (30 %)** (best **1869 banks = 23.9 s**);
> a **5 % anomaly universally, even worst-case/adversarial in 91 % of cells** (best worst-mode **146
> banks = 1.9 s**). `0.5 %` is below the detection floor everywhere. The mean/row channel is exactly
> blind to every such anomaly (`‖P U_β‖~1e-13`), so this remains a **covariance-only, F4-immune**
> demand-tier capability.

## Region fractions (T_det ≤ 4096 banks)

| anomaly ε | energy-spread (avg) | worst-mode | (buggy avg, for comparison) |
|:--:|:--:|:--:|:--:|
| 0.5 % | 0 % | 0 % | 0 % |
| 1 % | **30 %** | 2 % | (81 %) |
| **2 %** | **95 %** | 52 % | (100 %) |
| 5 % | **100 %** | 91 % | (100 %) |

By claim shell (energy-spread):

| claim | avg pass (1 %) | avg pass (2 %) | T_det_avg(2 %) range (banks) |
|---|:--:|:--:|:--:|
| 1.25 k_p | 56 % | 96 % | 467 – 5715 |
| 1.5 k_p | 33 % | 96 % | 644 – 6558 |
| 1.8 k_p | 0 % | 93 % | 1083 – 9125 |

Best cells (corrected, banks / seconds):
- **2 % anomaly:** k⁻², k_w=1, claim 1.25, σ_f=1.0 → **467 banks = 6.0 s** (d'≈? at 4096; T_det is the operative number).
- **1 % anomaly:** same family, σ_f=1.0 → **1869 banks = 23.9 s**.
- **5 % worst-mode:** **146 banks = 1.9 s**.

Note the best cell now sits at **σ_f=1.0** (not 0.3): once shot is physical, the higher realized
field contrast (σ_f,eff=0.70) helps the covariance SNR — the opposite of the estimation map's
σ_f-inert finding, because detection is throughput/shot-limited here, not mode-coverage-limited.

## What survives, plainly

- **The qualitative demand-tier claim survives:** you cannot image beyond `~1.1–1.25×`, but you can
  **detect** a `1.25–1.8×` beyond-band anomaly of **≥2 % energy in seconds** (6–120 s across the
  green cells), on a covariance channel the fresh-pattern/mean route provably cannot access.
- **The aggressive 1 % claim does not survive as universal** — it lives only at the best (k⁻²,
  narrow-claim, high-σ_f) cells. The honest headline is **"2 %-anomaly detection, nearly universal,
  seconds"**, with 1 % as a best-cell result and 0.5 % out of reach.
- Robustness (from `DETECTION_ROBUST.json`, recomputed under the corrected shot in the ROC run):
  aggregate detectability is invariant to 10 % medium-law basis rotation; τ×2 is a bank-count
  (diversity) cost only.

Full 81-cell table with both directions × 4 energies (banks + seconds): `DETECTION_MAP_CORRECTED.json`.
The Monte-Carlo ROC validation of these corrected `d'` values is in `DETECTION_ROC.md`. No git commit.
