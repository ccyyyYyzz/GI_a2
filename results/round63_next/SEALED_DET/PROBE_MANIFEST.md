# SEALED DETECTION PROBE — PROBE MANIFEST (R41 §4, frozen)

**Campaign:** R41 sealed detection confirmation — "a bucket detector can detect what it cannot image."
**Location:** `results/round63_next/SEALED_DET/` · **Date built:** 2026-07-24 · **No git commit.**
**Status:** machinery built + calibrated + frozen; **D0 + calibration-side D1 + simplex dry run PASS**.
The confirmatory (sealed) run awaits the coordinator's explicit freeze order. **No confirmatory scene,
code, medium, anomaly, or seed has been opened.**

---

## 1. Frozen physical ledger (R41 §4.2) — every choice load-bearing

| quantity | frozen value | source |
|---|---|---|
| scene grid | 64×64, `N=4096` | §4.2 |
| modulator band | `k_p=5` (Chebyshev `max(\|kx\|,\|ky\|)≤5`) | §4.2 |
| signed effective codes | `M=128` | §4.2 |
| physical realization | **256 complementary nonnegative exposures/bank** | §4.2 |
| pattern rate | 20 kHz | §4.2 |
| bank duration | **12.8 ms/bank** (256/20 kHz) | §4.2 |
| bank cap | `T_eff=4096` = **52.43 s** | §4.2 |
| photon level | **1e4 detected photoelectrons / physical bucket** | §4.2 |
| shot term | `R = diag(\|P\|x · mean(\|P\|x)/PHOT)` — **the corrected `\|P\|x` throughput** | CORRECTION_NOTE |
| medium | quasi-static within a bank; **independent between banks** (`T_eff`=indep. banks) | §4.2 |
| medium contrasts | nominal `{0.3,0.6,1.0}`, realized `sig_eff={0.298,0.503,0.696}` | §4.2 |
| spectra | `{flat, k^-1, k^-2}` | §4.2 |
| medium bands | `k_w/k_p ∈ {1,2,4}` | §4.2 |
| claim shells | `{1.25,1.5,1.8}·k_p` → beyond-band annulus upper freq `{6,7,9}` | §4.2 |
| beyond-band dof | `{1.25→48, 1.5→104, 1.8→240}` | derived |
| in-band nuisance dof | 121 | derived |
| anomaly magnitudes | `ε ∈ {0.005,0.01,0.02,0.05}` | §4.2 |
| medium correlation | `τ=8` banks, `ρ=exp(-1/τ)=0.8825` (lag/τ attribution only) | §4.2 |

**Mandatory everywhere:** the physical complementary-exposure shot term `|P|x`. The signed-difference
mean `Px (≈0)` is **never** used as photon throughput (this was the corrected bug; see CORRECTION_NOTE
and §6 disclosure). Shot-ledger validated to **1.15%** in D0.3.

## 2. Machinery (frozen scripts + sha256-16)

| file | role | sha256(16) |
|---|---|---|
| `sealed_common.py` | frozen ledger, corrected-shot Fisher (Engine A), Engine B, 4 tangents, Gaussian + physical-Poisson generators, 81/27-cell grids, frozen thresholds | `15abc6535409133a` |
| `arms.py` | the 7 arms (FIXED-COV, FRESH-COV-OPT full strength, FIXED/FRESH-MEAN, TRUE-LAW ORACLE, CROSSFIT-LAW, AMPLITUDE, LAG) | `0d64172ecf955c3e` |
| `simplex.py` | attribution simplex — 4 tangent spaces, canonical correlations, 5-class LDA | `8eedde19bc313a0a` |
| `sealed_banks.py` | five sealed banks + committed hashes + MANIFESTs | `08caf86a245823d7` |
| `bars.py` | D0–D7 evaluator + kill tree + fresh-vs-fixed 1.20× branch | `52afa6e5bc7e9efc` |
| `mc_plan.py` | 27-cell OA plan + compute-cost forecast | `2cb9396374c28cc8` |
| `run_d0.py` | D0 dry run (mechanism/engine integrity) | `cba967444d177f51` |
| `run_d1_cal.py` | calibration-side D1 + simplex Gram | `b1a93d24e0089b42` |

## 3. Two independent Fisher engines (D0.2)

- **Engine A** (`setup_cell`): einsum traces `I[a,b]=½ tr(V⁻¹ dV_a V⁻¹ dV_b)` + Schur-complement
  nuisance profiling (`J_B = I_bb − I_be I_ee⁺ I_eb`).
- **Engine B** (`fisher_engine_B`): the whitened vech-derivative construction of §3.1,
  `D̃_B = (I−Π_η) G D_B` with `GᵀG = V⁻¹`, `J_B = D̃_Bᵀ D̃_B` (QR projection out of the nuisance span).
- Agreement (D0.2): **max 4.8%** over `{lam_mean, lam_max, tr, d'_energy-spread, d'_matched}` across 5
  cells (bar ≤ 10%). Near-null `lam_min` at wide claims is excluded (both ≈0; not a physical mode).

## 4. Seven arms (R41 §4.5)

1. **FIXED-COV** — repeated sealed code bank `P_FIXED` (seed 10); running M×M covariance + profiled
   efficient beyond-band score `make_W`. Production detector.
2. **FRESH-COV-OPT** — fresh band-limited code bank each bank; **exact code-conditioned covariance
   score** `s_t=⟨r_tr_tᵀ−V_t, W_t⟩` per bank (NOT a coordinate-changing sample covariance). Full
   strength, GPU-batched (~4.5 ms/bank). Analytic `lam_eff=0.1168 ≈ fixed 0.1117` → **latency ratio
   0.956** (retain-concentration branch).
3. **FIXED-MEAN** — mean Mahalanobis on the repeated bank. Blind beyond-band (`d'≈0`, the wall).
4. **FRESH-MEAN** — strongest fresh-pattern mean detector; mean deflection **2.6e-17** (exact wall).
5. **TRUE-LAW ORACLE** — covariance score with the generating law (nondeployable ceiling).
6. **CROSSFIT-LAW** — declared law estimated from held-out baseline banks (deployable).
7. **AMPLITUDE / LAG** — lag-0 amplitude filter + lag-1 cross-bank filter (medium diagnostics).

## 5. Attribution simplex (R41 §9 Rank 2)

Four tangent spaces in the whitened joint observation space [mean ⊕ lag-0 cov ⊕ lag-1 cov],
joint dim 16 640: `mean/in-band`, `cov/beyond-band`, `cov-amplitude`, `cov-lag`.
- **Key specificity orthogonality (calibration):** beyond↔amplitude cc = **0.039**, beyond↔lag cc =
  **0.026**, both `< 0.10`. The efficient scores are mutually orthogonal by construction (max
  off-diagonal efficient cc = 0 → asymptotically independent Gaussian coordinates).
- **Five-class LDA balanced accuracy = 0.925 at T=2048 banks** (≥ 0.90, within the 4096 cap).

## 6. Five sealed banks — committed hashes (`banks/BANKS_INDEX.json`)

| bank | scenes | role | scene-set sha256(12) |
|---|---|---|---|
| calibration | 6 | filters/thresholds/cross-fit; **only bank opened** | `1deab0dbf1c1` |
| confirmatory | 12 (6 natural + 6 synthetic) | all primary endpoints — **SEALED** | `87da1b866ca3` |
| specificity | 3 + 6 event classes | H0/in-band/beyond/amplitude/τ×{0.5,2}/mixed — **SEALED** | `23105ff55c8a` |
| mismatch | 3 + 7 D5 axes | rotation/slope/band-edge/τ/shot/envelope/convolution — **SEALED** | `737e4ef22726` |
| oracle | 2 | exact mean-wall / true-law ceiling — **SEALED** | `4c7a25547ea3` |

**26 of 26 scenes disjoint** across banks (deterministic; no scene/seed crosses partitions). Salt bases
`{840000,841000,842000,843000,844000}`, disjoint from FOG_DMD_PROBE64 (`800000/810000/820000/830000`)
and DLGI (`700/710/720`). **Thresholds frozen before generation** (`sealed_common.BARS`).

## 7. Frozen bars (R41 §4.6) — `sealed_common.BARS`

| bar | frozen thresholds |
|---|---|
| D0 | mean-deriv ≤ 1e-10; two-engine ≤ 10%; MC shot/ledger ≤ 10% |
| D1 | `d'_emp/d'_an ∈ [0.80,1.20]`; median\|ARE\| ≤ 10%; none outside ±30% |
| D2 | ε2% ≥ 77/81 cells + MC LB > 0.90; best-cell T_det ≤ 600; ε5% all 81 + ≥90% worst-mode; ε1% best ≤ 2048; ε0.5% edge audit |
| D3 | target TPR ≥ 90% @ 1% FA; off-target FA ≤ 5% each; bal-acc ≥ 90%; beyond non-target \|d'\| ≤ 0.5; intended scores d' ≥ 5; simplex off-diag cc < 0.10 |
| D4 | FIXED retained if latency ≤ 1.20× FRESH; else switch to FRESH production |
| D5 | AUC loss ≤ 0.05; T_det inflation ≤ 25%; non-target FA ≤ 5% |
| D6 | online ≤ 0.1× bank (1.28 ms); memory ≤ 10 MB |
| D7 | one sealed run, no repair after unblinding |

## 8. Compute ceiling (R41 §4.8)

Forecast **1.21 GPU-hours** (ceiling 6.0) · storage **0.001 GB scalars** (ceiling 64 GB) ·
126 980 records. See `COST_FORECAST.md`. Exceeding the ceiling requires stopping and reporting.

## 9. Dry-run verdicts (calibration only; confirmatory sealed)

- **D0 = PASS** (`D0_REPORT.md`): wall 5.7e-16, two-engine 4.8%, shot+self-consistency 1.15%.
- **D1 (calibration) = PASS** (`SIMPLEX_CALIBRATION.md`): median ratio 1.022, median\|ARE\| 4.0%.
- **Simplex = PASS**: beyond ⟂ amplitude/lag (0.039/0.026 < 0.10); 5-class 0.925 @ T=2048.

**Kill-tree node 1 (D0) cleared.** The probe is ready for the sealed confirmatory run on order.
