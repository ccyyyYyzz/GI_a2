# MOLT probe — finishing run report

**MOLT = Microcell-Occupancy Laplace Tomography** (working name *saturation
jailbreak*). A SiPM bucket detector is C≈3600 binary microcells whose per-cell
Poisson (≥1 photon) saturation acts *before* the spatial sum, so a swept-power
"no-fire" curve encodes not just the ordinary linear bucket `p1 = Mx` but the
masked **quadratic** functional `p2 = M(x∘x)` — a second sensing operator from a
detector defect.

Rulings governing this run: `docs/ROUND63_GPT_ROUND31_RULING_RAW.md` (R31),
its Pro addendum `..._PRO_APPENDIX_0.md` (determinant law, `C_eff`, quenched
floor, `n_eff²` dose law), and — landed during this run and **binding** —
`docs/ROUND63_GPT_ROUND32_RULING_RAW.md` (R32), which adopts **equal wall-clock
time** as the primary frame, **retires dense-mask MOLT permanently**, and freezes
six go/no-go bars. This report is evaluated against those six bars.

---

## Conclusion (first)

> **MOLT does not die. Dense-mask MOLT is dead.** As a *sparse-support,
> dose-insensitive, switching-limited, equal-time* method the mathematics is
> fully verified and four of six binding bars pass outright; the two economic
> bars that touch hardware (wall-clock overhead, equal-time noisy reconstruction)
> are **conditional** on operating in the MHz-gated / slow-DMD regime and remain
> honest photon-for-rank costs. The matched-photon economics fail by orders of
> magnitude and are reported as a failure, exactly as R32 requires.

### Frozen thesis (R32, verbatim)

> Sparse MOLT trades photon dose for sensing rank: a microcell-occupancy sweep on
> support-targeted sparse masks supplies a mask-indexed quadratic operator in
> addition to ordinary linear buckets. At an expected dose floor of order `30×`
> the original linear campaign for 10%-precision `p2` when `n_eff ≤ 32`, the
> repeated gates can fit within one mask dwell in a MHz-gated, switching-limited
> system, allowing known or certified low-dimensional support fibers to close when
> linear sensing alone remains nonidentifiable.

---

## R32 six-bar verdict table

| # | Bar (R32 §4) | Result | Verdict | One-line evidence |
|---|---|---|---|---|
| 1 | `n_eff ≤ 32` on real sparse masks × scenes; report `C_eff`, mode number `k` | per-mask `n_eff ∈ [6.0, 19.0]` (all ≤32); `C_eff(iid)=3750≈C=3600`; `k=1` (fully-developed speckle) | **PASS** | `sparse_arm.json`, `quenched_floor.json` |
| 2 | 10% `p2` map at **≤ ×40** linear detected-photon budget/mask (expect ~×30) | production 3-level CRB **max ×11.3** over masks; dose-law boundary `n_eff=32 → ×28.8` (known-`p1`) confirms the ~×30 landing | **PASS** | `bars_dose.json` |
| 3 | favorable order-one `Δp2/p2` null pair `d'≥3` at **≤ ×4** dose (expect ~2.6) | dose-law `f=1` nuisance-`p1` **×2.64** (matches ruling ×2.6); aligned measured pair ×4.52 | **PASS** | `bars_dose.json` |
| 4 | ≤10% **end-to-end** wall-clock overhead incl. mask/settle/power/readout | ≤10% in **2/12** regimes: only switching-limited (1 MHz gates + 1 kHz DMD → **4.0%**). 1 MHz+10 kHz DMD → 29.5%; 100 kHz gates → 22–85% | **CONDITIONAL** | `timing_endtoend.json`, `fig_timing_overhead.png` |
| 5 | P0 mechanism identifiable from data at realistic noise | 3×3 model-recovery **99.7%** (coherent/product/correlated each ≥99%); equal-`p1`/equal-`p2` permutation moves `tr[(ΣV)²]` 198%, correlated `d'=3` at ×2.2, product exactly blind | **PASS** | `p0_modelfit.json`, `fig_p0_mechanism.png` |
| 6 | R32.1 Jacobian full-rank+conditioned on support **and** noisy joint recon materially beats strongest equal-time linear comparator | closure window **exact** (`rank[M_D|_S; M_S|_S; 2M_S|_S diag x_S]=s` for `40<s≤100`, linear stuck at 40); noiseless recon ~297 dB; **but** noisy equal-time recon only ties linear (+0.8 dB) | **CONDITIONAL** | `certificate_r321.json` |

Frame note (R32): the primary comparison is **equal wall-clock time**; the
**matched-photon** result below is a *mandatory adverse audit reported as a
failure*, not the headline.

---

## Honest bottom line (one paragraph)

The MOLT mathematics is real and now fully verified end to end: the masked
Laplace/power-sum identity recovers `p2` at ~5×10⁻⁷ relative error noiselessly
(committed T1), the finite-`C` covariance and saturation ridge `t*=2.8214` are
exact (committed T0b), the joint Jacobian doubles local rank `51→102` (committed
P3), and the `n_eff²` dose law is confirmed to a constant — the dense-vs-sparse
`B_sat` ratio equals `(n_eff ratio)²` to within 6–18% across six scenes. On
**matched photons** the method fails by orders of magnitude and we say so: dense
50% masks need `×1.6×10⁵ – ×1.1×10⁶` for a 10% `p2` map and `×2.84×10⁵` for a
null-pair `d'=3`; moving to `k=32`-support sparse masks (`n_eff` 6–19) cuts this
~100× via the `n_eff²` law, landing the decisive null-pair economics at **×126**
matched-photon — still a failure, but the production ridge design reaches a 10%
`p2` map at only **×11** and a favorable order-one null pair at **×2.6**, both
inside the R32 dose bars. On **matched time** the same photons are delivered in
`G=B/(C·t*)` gates at 0.1–1 MHz (`C·t*=10 157 ≈ B_LIN`, so one linear-campaign
photon budget ≈ one gate); a ×30 photon multiple is ~30 gates ≈ **30 µs/mask at
1 MHz**, but the *end-to-end* schedule including DMD switching, settling, power
changes and readout meets a ≤10% overhead budget **only** in the MHz-gated /
slow-DMD corner (4% at 1 MHz+1 kHz, but 30–85% once gates slow to 100 kHz or the
DMD reaches 22 kHz) — so "negligible time" is false in general and true only in
the declared regime. The mechanism is identifiable from data (99.7% three-model
recovery; an equal-`p1`/equal-`p2` spatial permutation exposes `tr[(ΣV)²]`, not
`p2`, as the correlated-speckle observable, which product/coherent are blind to).
The known-support jailbreak is exact **noiselessly** (committed T4: joint 296.7 dB
vs linear 27.8 dB at `s=80`, **+268.9 dB**) and the R32.1 hybrid certificate
closes the fiber geometrically for `K_D<s≤K_D+2K_S` — but the fiber-resolving
quadratic rows are Fisher-weak, so at the *strict equal-time* 10% `p2` dose the
noisy joint reconstruction only **ties** the strongest (averaging) linear
comparator: the win is provable RANK/identifiability, not an equal-time SNR win.
Net: MOLT is a photon-for-rank trade that is publishable only as the narrow
sparse/support-certified/equal-time method R32 authorizes.

---

## Claim fence (R32, verbatim)

**Authorized** — MOLT creates a second sensing operator by deliberately using
microcell occupancy saturation; sparse support-targeted masks improve the
economics through the `n_eff²` dose law; the joint operator can double local rank
for a common mask bank and can close known or certified low-dimensional prior
fibers; in switching-limited, dose-insensitive systems the photon-expensive sweep
may add little wall-clock time; the method is a **photon-for-rank trade**.

**Forbidden** — universal or dense-scene null-space jailbreak; exact recovery at
5% sampling without a support/prior certificate; photon efficiency, dose
neutrality, zero-cost sensing or zero calibration; "negligible time" based only on
MHz gate arithmetic without end-to-end timing; hiding the matched-photon failure;
extrapolating the `s=80` noiseless witness to unknown-support noisy recovery.

---

## Evidence by bar (finishing-run detail)

### Bar 1 — sparse `n_eff`, `C_eff`, mode number `k`
Six 32-support masks × scenes: `n_eff` = 15.0, 17.8, 17.9, 6.0, 19.0, 13.4 (all
≤32; ≤32 by construction since ≤32 active pixels). Measured effective sample
size `C_eff = Var(Y_c)/Var(q_C)`: **iid speckle → 3750 ≈ C=3600** (ratio 1.04);
block-correlated grains → `C/4`, `C/16` recovered (852, 220) — the honest
degradation channel. Mode number `k=1` (fully-developed speckle = one mode); a
bench `k>1` divides the `p2` signal by `k` (gate demand ~`k²`; committed P0:
`k3/k1=9.6`, `k10/k1=108.6`).

### Bar 2 — 10% `p2` map economics
Two independent estimates agree. Closed-form dose law (Pro D1) `B₁₀ = 281.3·n_eff²`
(known-`p1`, `t*=2.8214`), and the committed annealed-FIM three-level design give
per-mask multiples **×0.8–11.3** (max 11.3), well inside ×40. The ruling's ~×30
landing is the `n_eff=32` boundary (dose law → ×28.8 known-`p1`). *(Reported
adverse: the conservative flat 12-level discovery sweep in committed T2 costs
×330–3720/mask — a wasteful design that spends photons on information-poor anchor
levels; the production ridge design above is the correct economics.)*

### Bar 3 — favorable null-pair discrimination
Dose law (Pro D2) `B_{3σ}=25.33·n_eff²/f²`. For a **favorable order-one**
`f=|Δp2|/p2=1` pair at mean `n_eff=23.5`: known-`p1` **×1.40**, nuisance-`p1`
**×2.64** (matches ruling ~×2.6), both ≤×4. A box-constrained aligned null
perturbation only reached `f_max=0.25` (measured ×4.52), consistent with the
finding that interior-clamped null directions cap achievable `Δp2` — the dose law
at the bar's stipulated order-one `f` is the correct evaluation.

### Bar 4 — end-to-end wall-clock overhead
Model (µs): DMD switch `1e6/R_dmd`, settle 10, power switch 10 (2 per 3-level
sweep), readout 5/level; gate `1e6/R_gate`. Equal DMD-slot budget `K_D+K_S` so
switching is common-mode; MOLT adds the ridge-heavy sweep on the `K_S` sparse
slots.

| gates \ DMD (K_S=51) | 1 kHz | 10 kHz | 22 kHz |
|---|---|---|---|
| **100 kHz** | 22.1% | 72.8% | 84.6% |
| **1 MHz** | **4.0%** | 29.5% | 47.7% |

Halving the sparse bank (`K_S=25`) roughly halves overhead (2.7% at 1 MHz+1 kHz).
Meets ≤10% only in the switching-limited corner; **"negligible time" is false**
once gates slow to 100 kHz or the DMD reaches ≥10 kHz. Per-regime reporting is
mandatory. (`fig_timing_overhead.png`.)

### Bar 5 — P0 mechanism identifiability
Three nested physical models (all `q = det(I+aA)⁻¹`): rank-one **coherent**
(`A=√v √vᵀ`, `q=1/(1+a p1)`), diagonal **product** (`A=diag v`), low-rank
**correlated determinant** (block-coherent `A`, `q=∏_b 1/(1+a s_b)`). Fitting each
(one global-scale nuisance) to noisy sweeps from each generative model gives a
near-diagonal 3×3 confusion (300/300, 300/300, 297/300 → **99.7%**). Permutation
control: an equal-`p1`/equal-`p2` clustered-vs-spread permutation leaves
product & coherent curves **exactly invariant** (diff ~1e-16) but moves
`tr[(ΣV)²]=Σ_b s_b²` by 198%, detectable under the correlated model (`d'=2.0` at
`B=1e4`, `d'=3` at ×2.2) — so under correlated speckle the second-order observable
is `tr[(ΣV)²]`, **not** `p2`, unless `Σ` is diagonal/calibrated. (`fig_p0_mechanism.png`.)

### Bar 6 — R32.1 hybrid support certificate + equal-time reconstruction
Hybrid dense linear bank `M_D` (`K_D=40`) + support-targeted sparse swept bank
`M_S` (`K_S=30`), sparse `p1` equations **retained**. Closure window for
`rank[M_D|_S ; M_S|_S ; 2 M_S|_S diag x_S]`:

| `s` | rank(M_D\|_S) | rank(joint) | σ_min (unwt) | σ_min (Fisher-wt) | closes |
|---|---|---|---|---|---|
| 30 | 30 | 30 | 6.5e-1 | 6.1e-1 | both |
| 40 | 40 | 40 | 1.2e0 | 1.0e0 | both |
| 55 | 40 | 55 | 7.8e-1 | 5.5e-1 | **joint only** |
| 70 | 40 | 70 | 4.4e-1 | 9.9e-2 | **joint only** |
| 85 | 40 | 85 | 2.3e-1 | 2.5e-2 | **joint only** |
| 100 | 40 | 100 | 9.5e-3 | 9.0e-4 | **joint only** |
| 115 | 40 | **100** | 0 | 0 | neither |

The joint operator closes exactly `K_D < s ≤ K_D+2K_S = 100`; the Fisher-weighted
`σ_min` collapses toward the row-budget edge — the quadratic rows are
geometrically full-rank but carry `~1/n_eff` the information/photon (the
photon-for-rank cost, made visible). **Noiseless** joint reconstruction on the
support = 296–306 dB where the linear arm is underdetermined, confirming the fiber
closes. **Noisy equal-time**: the strongest linear comparator gets the *same DMD
slots re-spent as averaging* (not photon-starved); at the strict 10% `p2` dose the
joint recon only **ties** it (+0.8 dB at `s=80`) because box-constrained null
perturbations give sub-order-one `Δp2` → O(1) resolution SNR. Verdict:
**geometry+identifiability PASS, equal-time noisy material win NOT established** →
CONDITIONAL.

---

## Adverse audit panels (matched-photon; reported as failures per R32)

| Panel | Number | Source |
|---|---|---|
| dense 50% `p2` map (10%) | ×1.6×10⁵ – ×1.1×10⁶ / mask | committed T2 |
| dense null-pair `d'=3` | median `d'=0.0056`, ×2.84×10⁵ | committed P2 |
| dense null-space break (T3) | flat ×27 436 / optimal ×9 311 of the linear campaign | committed T3 |
| **sparse** null-pair `d'=3` | median `d'=0.267`, ×126 (×2254 cheaper than dense) | `sparse_arm.json` |
| known-support witness `s=80` | joint 296.7 dB vs linear 27.8 dB = **+268.9 dB** *(noiseless)* | committed T4 |

---

## Files (all under `results/round63_next/SATURATION_JAILBREAK_PROBE/`)

Finishing-run code (imports/extends the committed framework, does not rewrite it):
`code/sparse_arm.py` (bars 1 econ + dual-resource), `code/bars_dose.py` (bars 2–3),
`code/p0_modelfit.py` (bar 5), `code/quenched_floor.py` (quenched floor + refresh +
`C_eff`), `code/timing_endtoend.py` (bar 4), `code/certificate_r321.py` (bar 6),
`code/update_results.py`.
Results: `sparse_arm.json`, `bars_dose.json`, `p0_modelfit.json`,
`quenched_floor.json`, `timing_endtoend.json`, `certificate_r321.json`, and the
merged `jailbreak_results.json` (`FINISHING_RUN` + `R32_SIX_BAR_VERDICT`).
New figures: `figs/fig_quenched_floor.png`, `figs/fig_timing_overhead.png`.
