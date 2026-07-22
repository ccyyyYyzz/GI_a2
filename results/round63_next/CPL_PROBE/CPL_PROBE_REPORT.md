# CPL-GI micro-probe — report

**Lane:** R29 §6 reserve (the LAST software-method candidate), activated by the recorded
`DEV_GATE_FAIL` of DOPS-GI (`results/round63_next/SCHEDULE_PROBE/DEV_GATE_VERDICT.json`,
Probe B). Governing spec: `docs/ROUND63_GPT_ROUND29_RULING_RAW.md` §2.1 (method) + §6 (kill bar).

**Verdict: `CPL_KILL`.** At the frozen primary declared-coherence cell, CPL-GI's median
PSNR gain over the best of the naive/ratio baselines is **+0.412 dB** (5/6 scenes positive) —
below the **≥ 1.0 dB** kill bar. Per R29 §6, this being the last reserve candidate, the frozen
program conclusion is: **no current software candidate clears the operator's
novelty × simplicity × generality × image-positive bar.**

- Run: 2026-07-22 10:18:43Z · Python 3.11.5 · NumPy 1.24.4 · CPU · **110.3 s** wall (target ≤ 1–2 h).
- Code: `results/round63_next/CPL_PROBE/cpl_gi_probe.py` · results: `cpl_probe_results.json` · gate: `CPL_GATE_VERDICT.json`.
- Read-only on all inputs; all writes confined to `results/round63_next/CPL_PROBE/`. Not committed.

---

## 1. Frozen design (declared before the drift cells were run)

**Method (CPL-GI, R29 §2.1).** Complementary single-pixel acquisition. Hadamard rows
`h ∈ {−1,+1}^1024` (full 1024-row, **sequency order**, row 0 = DC) realized as physical pattern
pairs `m+ = (1+h)/2`, `m− = (1−h)/2` (both in {0,1}; constant-sum `m+ + m− = 1`). Each Hadamard
coefficient is one **packet of q = 2 adjacent complementary exposures** projected back-to-back →
**2048 physical exposures**. Bucket counts under multiplicative gain drift:
`Y+_b ~ Poisson(a_{2b}·φ·m+_bᵀx)`, `Y−_b ~ Poisson(a_{2b+1}·φ·m−_bᵀx)`. Because `m+ + m− = 1`,
conditioning on the packet total `N_b = Y+_b + Y−_b` removes the packet gain **exactly when the
gain is constant within the packet**: `Y+_b | N_b ~ Binomial(N_b, m+_bᵀx / 1ᵀx)`. CPL reconstructs
by maximizing the exact conditional (binomial) log-likelihood over `x ≥ 0`, plus a TV penalty.

The two packet frames are simulated as **adjacent samples of the OU path** (`a_{2b}`, `a_{2b+1}`),
i.e. the within-packet gain is *not* forced equal — the method's approximation error is included
honestly and is precisely what the primary cell stresses.

- **Scenes:** the 6 DEV scenes `bridge_contour_1, bridge_microtex_1, bridge_twopop_1, bridge_control_2, bridge_contour_3, bridge_microtex_3` (32×32, values in [0,1]), from `data/r63_bridge_scenes/`. SHA-256 (first 16 hex) recorded in `cpl_probe_results.json`.
- **Patterns:** full 1024-row sequency Hadamard, 2048 complementary exposures. **Identical exposures, photons, and pattern content for all arms; only the processing differs** (data simulated once per (scene, seed, cell) and fed to every arm).
- **Gain model:** the same frozen OU log-gain as Probe B — exact AR(1), `μ = −σ²/2` (so `E[a]=1`), lifted verbatim from `results/round63_next/SCHEDULE_PROBE/probe_b.py` (`ou_path`).
- **Photon scale:** `φ = 11.9375`, set so the mean detected counts per exposure ≈ the Probe B bucket scale (~2200 counts/exposure; mean `m±ᵀx = S/2 = 184.29` over the 6 scenes).
- **Seeds:** 5 paired gain/noise seeds (`rng(1000+seed)` gain, `rng(5000+seed)` Poisson — identical offsets to Probe B), shared across arms.

**Cells (frozen now).** PRIMARY (gated) `σ = 40%, t_c = 16 frames` — the declared-coherence cell:
within-packet correlation `exp(−1/16) ≈ 0.94`, so the packet assumption holds by design.
Secondary descriptive cells (no gate): `t_c = 2` (assumption-stressed), `t_c = 64` (easy),
`σ = 15% at t_c = 16` (milder), and `σ = 0` (gain-free no-harm control).

**Arms** (identical data per seed; processing-only differences):
- **A1 naive differential** — `c_b = Y+_b − Y−_b`; standard Hadamard inversion (`x = Hᵀc/(1024·φ)`), then TV denoise. Drift-vulnerable baseline.
- **A2 adjacent-pair ratio normalization** (Sensors-2023-style) — `r_b = (Y+_b − Y−_b)/(Y+_b + Y−_b)`, scaled by the global mean packet total `N̄/φ`, then the same inversion + TV denoise as A1.
- **A3 CPL-GI** (the method) — exact conditional binomial MLE + TV. Flux pinned by a firm quadratic anchor to `Ŝ = N̄/φ` (the stationarity flux anchor, `E[a]=1`); **the same `Ŝ` is A2's global scale**, so any flux bias affects A2 and A3 identically.
- **A4 oracle** (descriptive only) — divide each exposure by the TRUE `a_t`, then A1. Headroom reference; never a comparator.

**Frozen TV / λ rule.** Smoothed isotropic TV, `ε = 1e-3`, identical form for every arm.
The λ **rule** is identical for all arms: *per-arm λ chosen to maximize mean control-cell (σ=0)
PSNR over the shared dimensionless grid `{0, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1}`, frozen before any
drift cell is evaluated.* Because the binomial NLL's data-term curvature differs from the
unit-curvature denoise term, A3's grid is mapped through a single analytic curvature scale
`H_SCALE = 75.07` (mean NLL Hessian diagonal at the uniform anchored image, averaged over the
control scenes/seeds) so that the same grid means the same smoothness/data trade-off for every arm.
Selected: **A1 = A2 = A4 → λ_img = 0.1; A3 → λ_img = 0.03.** (Control-cell PSNR-vs-λ curves are in §5.)

**Kill bar (R29 §6, frozen).** At the PRIMARY cell: median PSNR gain of CPL over the **best of
{A1, A2}** per scene (five-seed means) **≥ 1.0 dB**, AND **≥ 5/6** scenes positive. No-harm at
σ=0: CPL vs best comparator median **≥ −0.25 dB**. Anything less = KILL.

---

## 2. Numerical sanity (pre-grid)

1. **TV gradient** vs finite differences: max abs error `2.4e-8` (analytic gradient correct).
2. **No-drift, no-noise, λ = 0** (exact expected counts, `a ≡ 1`): every arm recovers every scene at **283–306 dB**, and A3's achieved flux equals `Ŝ` to `≤ 2.2e-16` relative error — the exact conditional MLE + flux anchor is correctly implemented and identifiable. (A3 flux-anchor residual stays `≤ 1.3e-3` across all noisy drift runs, i.e. the scale is effectively hard-pinned.)

---

## 3. PRIMARY cell (σ=40%, t_c=16) — the gate

Five-seed mean PSNR (dB); `dQ = CPL − best{A1,A2}`. At this cell `best{A1,A2} = A2` for every scene.

| scene | A1 naive | A2 ratio | **A3 CPL** | A4 oracle | dQ (CPL−best) |
|---|---:|---:|---:|---:|---:|
| bridge_contour_1  | 1.59 | 2.18 | **3.25** | 16.29 | **+1.070** |
| bridge_microtex_1 | 10.79 | 11.34 | **10.94** | 19.81 | **−0.405** |
| bridge_twopop_1   | 7.46 | 8.26 | **8.37** | 23.08 | **+0.101** |
| bridge_control_2  | 9.24 | 9.95 | **10.30** | 17.47 | **+0.350** |
| bridge_contour_3  | 1.52 | 2.10 | **3.18** | 16.82 | **+1.081** |
| bridge_microtex_3 | 6.35 | 6.84 | **7.32** | 16.55 | **+0.474** |

- **median dQ = +0.412 dB** — condition (a) `≥ 1.0 dB` **FAILS**.
- **n_positive = 5/6** — condition (b) `≥ 5/6` passes.
- no-harm (σ=0) median **+0.213 dB** — condition (c) `≥ −0.25 dB` passes.

**→ VERDICT `CPL_KILL`** (fails the primary margin condition).

Two facts frame the result. (i) CPL's only clear wins are the low-texture **contour** scenes
(+1.07, +1.08 dB); it *loses* on `microtex_1` (−0.41). (ii) The known-gain **oracle sits 10–15 dB
above every software arm** at this cell — the drift residual CPL cannot remove dominates the error
budget, and the exact likelihood captures only a sliver of the recoverable headroom.

---

## 4. Secondary cells (descriptive, no gate)

Median `dQ = CPL − best{A1,A2}` (five-seed means), frozen protocol vs the λ = 0 (no-TV) view:

| cell | within-pkt corr | median dQ (frozen λ) | n_pos | median dQ (λ=0) | n_pos |
|---|---|---:|---:|---:|---:|
| **primary** σ40, t_c16 | 0.94 | **+0.412** | 5/6 | +1.585 | 6/6 |
| sec σ40, t_c2 (stressed) | 0.61 | +3.880 | 6/6 | +4.559 | 6/6 |
| sec σ40, t_c64 (easy) | 0.98 | **−0.224** | 3/6 | +0.444 | 6/6 |
| sec σ15, t_c16 (mild) | 0.94 | **−0.296** | 3/6 | +0.411 | 6/6 |
| control σ0 (no drift) | — | +0.213 | 5/6 | +0.117 | 6/6 |

Under the frozen protocol CPL's advantage is **not robust across regimes**:

- **Large only where the comparator has collapsed** — at `t_c = 2` (severe high-frequency drift) CPL leads by +3.9 dB, but there the *absolute* PSNRs are ~ −6 to +6 dB (A2 negative on the contour scenes): unusable images either way.
- **Marginal at the primary** declared-coherence cell (+0.41 dB).
- **Negative where imaging is actually feasible** — at the easy (`t_c = 64`, +14–16 dB regime) and mild (`σ = 15%`) cells CPL is **−0.22 / −0.30 dB** vs the ratio baseline, losing on the higher-SNR texture scenes.

A deployable drift correction should help where imaging works, not only where it is already broken.

---

## 5. Why the gain evaporates under regularization (the R29 §2.1 kill risk, realized)

The **raw** (λ = 0) exact-MLE advantage over ratio-normalization is a robust **+1.585 dB** at the
primary cell (all six scenes +1.49 to +1.78 dB). Control-calibrated TV regularization — which every
deployable pipeline uses, and which the comparators use too — **shrinks this to +0.412 dB**, because
TV independently cleans the drift/shot-noise artifacts that CPL's likelihood was correcting. The
comparators gain *more* from TV than CPL does:

- `bridge_contour_1`: A2 improves 1.03 → 2.18 dB (+1.15) under TV while A3 improves 2.52 → 3.25 (+0.73); gap 1.49 → 1.07.
- `bridge_microtex_1`: A2 improves 8.56 → 11.34 (+2.78) and **overtakes** A3 (10.21 → 10.94); gap +1.65 → −0.41.

Control-cell PSNR-vs-λ curves (mean over 6 scenes × 5 seeds), used to freeze the per-arm λ:

| λ_img | 0 | 3e-3 | 1e-2 | 3e-2 | 1e-1 | 3e-1 |
|---|---:|---:|---:|---:|---:|---:|
| A1 | 16.69 | 16.86 | 17.24 | 18.16 | **18.86** | 15.17 |
| A2 | 16.63 | 16.80 | 17.18 | 18.11 | **18.87** | 15.18 |
| A3 | 16.92 | 17.25 | 17.94 | **19.13** | 16.53 | 12.83 |

This is exactly the kill risk flagged in R29 §2.1 ("reviewers may regard the method as a
statistically dressed version of differential Hadamard normalization") materializing at image
level: once a standard image prior is applied on equal terms, the exact conditional likelihood does
not add a deployable margin over the published ratio normalization at the declared-coherence cell.

The λ = 0 view is reported for full transparency: **the verdict is sensitive to regularization** —
without any prior CPL would clear the 1.0 dB bar (+1.585). But the kill bar is defined for the
regularized reconstruction with a λ frozen on the control cell before any drift outcome is seen
(R29 §6 / task spec), and under that pre-registered protocol CPL fails. No softening is applied in
either direction: the governing number is +0.412 dB → KILL.

---

## 6. No-harm (σ = 0 control) and honesty notes

- No-harm median **+0.213 dB** ≥ −0.25 (passes), but per-scene it is heterogeneous: CPL gains on `bridge_control_2` (+2.44) and the microtexture scenes (+0.42, +0.70) yet **loses −2.06 dB on `bridge_twopop_1`**. The loss is an estimator/regularization interaction (A3's lighter λ_img=0.03 + flux-anchored joint reconstruction underperforms strongly-smoothed Hadamard inversion on that high-contrast two-population scene), not a drift effect.
- "Expected-under-null" outcome, recorded plainly (per task): at the declared-coherence cell A2 (the published normalization) already captures essentially all of the deployable drift robustness; the exact likelihood adds **< 0.5 dB** median over it. The interesting question — whether exact likelihood beats ratio normalization by a *deployable* margin — is answered **no** under a fair, pre-registered, regularized protocol.

---

## 7. Conclusion

CPL-GI **fails** its frozen R29 §6 kill bar at the primary declared-coherence cell
(median +0.412 dB < 1.0 dB; the raw un-regularized margin of +1.6 dB does not survive an
equal-terms image prior; the advantage reverses on the higher-SNR feasible-imaging cells; and the
known-gain oracle shows 10–15 dB of headroom that CPL does not capture). Verdict: **`CPL_KILL`**.

As the last software reserve candidate, this closes the R29 triage. Following R29 §6 verbatim:

> **No current software candidate clears the operator's novelty × simplicity × generality ×
> image-positive bar.** The software layer for these frozen channels is saturated at the requested
> effect size; the correct next move is a genuinely new observation model or application, not a
> relabeled established correction.
