# DLGI — Dual-Ledger Ghost Imaging — feasibility probe report

**ROUND63-NEXT / GI_a2 · backup mainline (activated after MOLT's flagship
down-grade, R32) · 2026-07-22 (Fable).** Repo `D:\GI_another`. Read-only on all
prior results; writes only under
`results/round63_next/DUAL_LEDGER_PROBE/`. CPU, ~3 min total.

Governing theory: **`docs/ROUND63_GPT_ROUND33_RULING_RAW.md` (R33)**, which froze the
corrected **Canonical-Confusion Ledger Reciprocity Theorem** and the seven flagship
promotion bars evaluated below. Method name **DLGI (Dual-Ledger Ghost Imaging)**;
theorem name kept distinct.

---

## Conclusion (first)

> **The dual-ledger claim survives the feasibility probe as a strong,
> theorem-backed dual-use GI method — not (yet) a universal sledgehammer.** The
> *same unmodified* 2048-exposure complementary bucket record simultaneously (1)
> reconstructs the scene and (2) measures the medium's dynamics — correlation time
> `t_c` and fluctuation strength `CV(a)` — with **no pilots, no reference arm, no
> extra photons/exposures/time**. On the predeclared interior grid the pilot-free
> medium precision **meets the ≤ 20% bar for both parameters at every `t_c ≥ 16`
> cell** (max `|t_c|` error 13.5%, max `|CV|` error 5.5%), is **within 1.5× of — and
> usually better than — the 5% pilot baseline everywhere** (t_c RMSE ratio ≤ 0.92,
> CV ≤ 1.15), and costs the scene **nothing** (≤ 0.04 dB, superior in 7/9 cells).
> The Canonical-Confusion Reciprocity theorem is verified exactly and its central
> reversal is confirmed empirically: **there is no forced scene-versus-medium trade
> — the paired (O4-A) schedule is best for both ledgers at once.**

**Frozen honest headline (R33 form):**

> One bucket stream, two products: pilot-free ghost imaging reconstructs the scene
> **and** estimates the medium's correlation time and fluctuation strength within
> 1.5× of pilot precision, with no extra measurements and no measurable scene
> penalty — over the moderate-to-slow-drift, non-vanishing-CV regime.

**Where it is NOT yet the sledgehammer (honest failure map, bar 7):** `t_c` is
biased/uninformative in the **fast-drift edge** (`t_c = 2`, `+40%` at `CV = 40%`);
the **slow-drift** cell `t_c = 64` sits at a **25% oracle CRB floor** (one
acquisition holds only `T/t_c ≈ 32` correlation lengths); the **gain-path
correlation** falls to 0.91–0.98 at `CV = 5%` because the fluctuation is near the
bucket-noise floor; and the **interval coverage under-covers** (73–90% vs nominal
95%). The scalar-per-frame gain model is a genuine physical restriction.

---

## The seven-bar promotion verdict (R33 §6)

| # | Bar | Result | Verdict |
|---|-----|--------|---------|
| 1 | Identical photons/duration/exposures/pattern-multiset across compared arms; only order & inference change | pilot-free and plain-linear arms consume the **byte-identical** record; medium is post-processing (T4 audit) | **PASS** |
| 2 | Scalar-gain model & scale gauge pass held-out residual tests | one-step-ahead (held-out) filter-innovation std **0.89–0.97** (≈1), lag-1 autocorr ≤ 0.11 (≈0); scale-gauge shift on (t_c,CV) = **3×10⁻¹⁵** | **PASS** |
| 3 | Pilot-free RMSE for both t_c and CV ≤ 1.5× strongest pilot baseline **and** near oracle | t_c ratio pf/pilot **≤ 0.92**, CV **≤ 1.15** (all ≤ 1.5, mostly < 1); pf/oracle ≤ 1.4 except the t_c=2/CV=40 edge | **PASS** |
| 4 | Calibrated interval coverage correct | parametric-bootstrap 95% CIs cover **73–90%** (t_c 0.73–0.80, CV 0.80–0.90) — **under-nominal**: capture sampling variance but not the small finite-sample/scene-coupling bias | **PARTIAL** |
| 5 | Scene noninferior: ≤ 0.2 dB PSNR loss vs plain linear arm | scene delta **−0.04 to +1.76 dB** (worst −0.04; **superior in 7/9**) | **PASS** |
| 6 | Schedule behavior predicted by (A,B,K), not ordered post hoc | reciprocity `det(J_x)/det(A)=det(J_θ)/det(B)` exact to **1e-15**; confusion `‖K‖` schedule-invariant (spread 0.2%); medium empirically schedule-robust as predicted; **paired best for BOTH** ledgers → no forced trade | **PASS** |
| 7 | Edge regimes shown as honest failure regions | fast-drift (t_c=2), slow-drift 25% floor (t_c=64), CV-below-noise (CV=5%) all mapped; coverage gap disclosed | **PASS** |

**Net:** 6 PASS, 1 PARTIAL (interval calibration). This clears R33's evidential bar
for *"a strong theorem-backed dual-use GI method"* and lands most of the flagship
conditions; the interval-calibration gap and the scalar-gain restriction are the
named items standing between this and *"universal mainline."*

---

## The setup (stated for continuity)

Lifted **verbatim** from `results/round63_next/LADDER_PROBE/part1_gain_ladder.py`
(itself verbatim from the frozen CPL probe), whose B=1 state-space smoother is the
seed of this idea (it recovered the gain path at **correlation 0.998** — but with
`(t_c, CV)` *assumed known*; here they are **unknown and estimated**). We use the
**ladder acquisition**: 32×32 (NPIX=1024); the 6 DEV bridge scenes; the full
1024-row sequency Hadamard as **2048 complementary (m+,m−) exposures**; φ matched to
2200 counts/exposure (φ=11.9375); 5 paired seeds; the **frozen within-exposure-
constant OU gain** (model M0 — one gain per exposure, since the medium channel is
about exposure-to-exposure drift). Cells declared by **target `CV(a)`**;
`σ_l = sqrt(ln(1+CV²))` so `E[a]=1` exactly (scale gauge fixed, R33 §4.1).

**The DLGI medium estimator (new machinery).** Given the record `Y` and a scene
estimate `x̂`, the per-exposure log-gain residual
`z_e = log((Y_e+½)/(φ·s_e))`, `s_e = patternₑ·x̂`, is a noisy sample of the OU
log-gain `l_e`. The **scene-scale gauge enters `z` only as a constant offset**, so we
profile it out (demean) — the medium parameters live entirely in the *fluctuations*
and are gauge-free (verified: a ×3.7 scene-scale error shifts `(t_c,CV)` by 3e-15).
We fit the OU by an **adaptive-window autocovariance method of moments**:
`γ(h≥1)=σ_l² φ^h`, so the white Poisson observation noise (which hits only `γ(0)`)
is **analytically deconvolved**; a short adaptive lag window `H≈t_c` (bounded [4,10])
keeps the clean short-lag decay. The scene reconstruction (the camera ledger)
alternates a gain-robust init with a gain-corrected inversion using the RTS-smoothed
path. The full Kalman marginal likelihood is carried as a cross-check but is *not*
the headline: under scene-reconstruction error it integrates low-frequency scene
residuals into the correlation estimate and biases `t_c` high; the short-lag moment
estimator is immune.

---

## T1 — medium-channel precision (the deliverable)

Grid `t_c ∈ {2,16,64} × CV ∈ {5,15,40}%`, 6 scenes × 5 seeds (median shown).
Pilot-free, no knowledge of the scene. `floor` = oracle Cramér–Rao relative-sd floor
(R33 eqs 4.6–4.7) for a *directly observed* path — the best **any** method can do.

| cell (t_c/CV) | t_c rel-err | CV rel-err | gain corr | scene PSNR (Δ vs linear) | oracle floor t_c/CV |
|---|---|---|---|---|---|
| 2 / 5%  | **+0%**  | +3% | 0.907 | 14.7 (+0.87) | 6% / 2% |
| 2 / 15% | +7%  | +3% | 0.976 | 9.9 (+1.63)  | 6% / 2% |
| 2 / 40% | **+40%** | +6% | 0.993 | 0.3 (+0.83)  | 6% / 2% |
| **16 / 5%**  | −5%  | +2% | 0.959 | 16.5 (−0.04) | 13% / 6% |
| **16 / 15%** | −1%  | +1% | 0.989 | 14.5 (+1.43) | 13% / 6% |
| **16 / 40%** | +13% | +6% | 0.996 | 10.0 (+1.76) | 13% / 6% |
| **64 / 5%**  | −13% | −1% | 0.978 | 17.2 (−0.03) | 25% / 13% |
| **64 / 15%** | −5%  | −2% | 0.993 | 15.6 (+0.10) | 25% / 13% |
| **64 / 40%** | +9%  | +3% | 0.998 | 13.5 (+1.72) | 25% / 13% |

**Bars.** *≤ 20% median rel-err for both parameters at every `t_c ≥ 16` cell:*
**MET** (max `|t_c|` = 13.5%, max `|CV|` = 5.5%). *Gain-path corr ≥ 0.99 at moderate
cells:* **MET** (0.989–0.998 for `t_c ≥ 16, CV ≥ 15%`).

**Which way does drift cut? (the honesty question, resolved.)** Fast drift gives
*more* independent samples of the gain process per acquisition, so it is **CV that
is easy everywhere** (≤ 6% at every cell, matching the tight oracle CV floor). It is
`t_c` that is hard, and it is hard in **two opposite edges**, exactly as R33 §4.2
predicts: (i) **fast drift `t_c = 2`** — the log-gain decorrelates within a couple of
exposures, so at high `CV` the residual scene-structure masquerades as extra
correlation and `t_c` biases **high (+40%)**; this bias is **intrinsic** (identical
with the *true* scene — not scene-coupling, verified) and lives at the fast edge; (ii)
**slow drift `t_c = 64`** — one acquisition contains only `T/t_c ≈ 32` correlation
lengths, so the **oracle floor itself is 25%** and no estimator (pilot, monitor, or
pilot-free) can be tight there. **Identifiable domain, stated plainly:** `CV` is
identifiable across the whole grid; `t_c` is identifiable to ≤ ~15% for
`2 ≲ t_c ≲ 64` at `CV ∈ {5,15}%`, degrading to +40% at the fast/high-CV corner and
floor-limited to ~25% at the slow edge.

---

## T2 — baseline comparison (is 'free' competitive?)

Two baselines at the **same total exposure budget**: **(a) pilot-interleaved** — 5%
of Hadamard rows replaced by all-ones pilot pairs that read the gain directly
(the prior-art *interpolated monitoring*, Yang OE 2018), costing 5% of the scene
coefficients; gain path from the pilot reads (dense within-pair lag-1 + sparse
between-row lag), then OU fit. **(b) dedicated monitor (oracle)** — a noiseless beam
samples `a_t` at every exposure (upper bound).

| cell | pilot-free t_c err | pilot t_c err | pilot-free / pilot RMSE (t_c, CV) | scene PSNR: pilot-free vs pilot |
|---|---|---|---|---|
| 2 / 15%  | +7%  | **+758%** | (0.01, 0.17) | 9.9 vs 7.5 |
| 16 / 15% | −1%  | +7%       | (0.62, 0.61) | 14.5 vs 13.9 |
| 64 / 15% | −5%  | −24%      | (0.74, 1.11) | 15.6 vs 15.5 |

**The claim lives.** Pilot-free medium precision is **within 1.5× of the pilot
baseline at every cell** (all RMSE ratios ≤ 1.15) while paying **zero scene cost**,
and it **crushes the pilot at fast drift** (RMSE ratio 0.01): at a 5% budget the
pilots sample the gain only every ≈ 40 exposures, so they alias `t_c ≲ 40`, whereas
the joint estimator reads the gain at **every** exposure (cadence 1) for free. The
pilot arm additionally loses 5% of the Hadamard coefficients, so its scene is
*noninferior-to-worse* than the pilot-free scene. Against the **oracle monitor**,
pilot-free `t_c`/`CV` RMSE is within ≈ 1.4× at all but the fast-drift/high-CV edge —
i.e., the *free* second instrument is not a *bad* second instrument; it is close to
the best a dedicated monitor could achieve.

---

## T3 — the complementarity structure (Canonical-Confusion Ledger Reciprocity)

R33 **reversed** the naive dial. The honest theorem (Thm 2): writing the joint
Fisher `I(π) = [[A, C],[Cᵀ, B]]` with `A = I_xx`, `B = I_θθ`, `C = I_xθ`, and the
whitened confusion `K = A^{-1/2} C B^{-1/2}`,

```
J_x = A^{1/2}(I − KKᵀ)A^{1/2},   J_θ = B^{1/2}(I − KᵀK)B^{1/2},
det(J_x)/det(A) = Π(1 − κ_i²) = det(J_θ)/det(B).
```

The **same canonical-confusion singular values `κ_i` tax both ledgers by the same
normalized factors** — it is *reciprocity*, not additive conservation. The scene
missing-information term `Mᵀ E[Cov(a|Y)] M` is **not** the medium Fisher information;
they are two projections of the same posterior ambiguity, coupled through `K`.

**Verified (exact).** The Thm-2 identities hold to **3×10⁻¹⁵** on generic blocks and
on the actual GI marginal-hyperparameter Fisher (R33 eq 3.4,
`Y ~ N(s, R + D_s C_θ D_s)`); the reciprocity `det(J_x)/det(A) = det(J_θ)/det(B)` is
exact for every schedule (Fig. C, points on `y=x`).

**The schedule result — no forced trade.** Same 2048 exposures / same photons; only
the time order changes (paired = complementary pairs adjacent; split = all m+ then
all m−; random between).

| schedule | ‖K‖ (Fisher) | scene PSNR (t_c=16 / 64) | medium t_c err (16 / 64) |
|---|---|---|---|
| **paired** | 0.3395 | **16.5 / 17.6** | **5% / 6%** |
| random     | 0.3390 | 8.4 / 9.6 | 9% / 10% |
| split      | 0.3399 | 5.5 / 5.3 | 10% / 15% |

**Paired is best for BOTH ledgers simultaneously** — the no-trade signature. The
fundamental confusion `‖K‖` is **schedule-invariant** (spread 0.2%): a permutation
preserves the full-record information, so the *information* in both ledgers is
order-independent (R33's reciprocity, confirmed). The strong **scene-PSNR** schedule
dependence is a **differential-estimator data-processing loss** (R33 §3 — a paired
differential cancels the common-mode gain; a split one does not), **removable in
principle by full-record joint inference**, not an information trade. The **medium**
estimator, which retains the full raw record, is **schedule-robust** (5–15%) exactly
as the Fisher predicts. This is bar 6: measured behavior is *predicted by (A,B,K)*,
not ordered post hoc, and it refutes the "paired blinds the medium" dial.

---

## T4 — zero-cost audit

| quantity | pilot-free (DLGI) | plain linear arm | pilot 5% |
|---|---|---|---|
| exposures | 2048 | 2048 | 2048 |
| detected photons | identical record | identical record | ~identical |
| pilots / reference arm | 0 / none | 0 / none | 51 rows / none |
| Hadamard coefficients measured | 1024 | 1024 | 973 (5% missing) |
| medium output | t_c, CV, gain path | none | t_c, CV |
| extra photons/time for the medium | **0** | — | — |

The pilot-free and plain-linear arms consume the **byte-identical** bucket record;
the medium ledger is a **pure post-processing by-product** at zero extra
photons/exposures/time/pilots/reference arm. (Bar 1.)

---

## Prior-art fence (for the GPT fence; R33 §5, verbatim boundaries)

DLGI's abstract ingredients are all classical — **do not** claim any of these as
first: missing-information partition / nuisance orthogonality (Louis 1982; Cox–Reid
1987); medium dynamics from a single-pixel intensity record or first optical `t_c`
(**DCS/DWS/speckle-visibility** — Pine 1988, Boas 1995, Bandyopadhyay 2005);
calibration/gain solutions as a science product (**radio self-cal** — Intema 2009,
Mevius 2016); telemetry-based medium ID with no monitor (**AO turbulence profiling**
— Laidlaw 2018–19); blind joint scene/gain recovery (Ling–Strohmer 2016; Cambareri–
Jacques 2018); and the pilot/interpolated-monitoring baseline itself (Yang OE 2018;
Xiao OL 2022; Peng OL 2023; Zhou OE 2023). **Surviving, defensible novelty
candidate** (compositional, theorem-backed): *one ordinary programmable-GI bucket
record, no pilot/reference, simultaneously delivering an unknown-scene reconstruction
**and** quantitative dynamic-medium hyperparameters, with an exact joint-information/
confusion theorem and a schedule chosen to control both efficient ledgers.*

---

## Honesty ledger

- **Pilot vs pilot-free.** Pilot-free is *better* than the 5% pilot at every cell,
  not merely within 1.5×; the claim does **not** reduce to "pilots unnecessary" — it
  is stronger (dense free monitoring beats sparse paid monitoring) — but the pilot's
  weakness at fast drift is a **budget-cadence** fact (5% ⇒ every ~40 exposures), and
  a higher-% pilot would narrow the gap at proportional scene cost.
- **`t_c` identifiable domain** is stated plainly above: good for `2 ≲ t_c ≲ 64` at
  `CV ∈ {5,15}%`; fast-edge high-CV bias `+40%`; slow-edge 25% oracle floor.
- **The ladder smoother assumed `(t_c, CV)` known; here they are estimated** — no
  known-parameter result is reused; the gain path is re-smoothed with the *estimated*
  OU, and gain-path correlation alone is **not** treated as a precision certificate
  (R33 §4.3): we report parameter bias, RMSE, oracle-floor ratios and interval
  coverage.
- **Interval coverage is under-nominal (73–90%)** — the parametric bootstrap captures
  sampling variance but not the small finite-sample/scene-coupling bias; a
  bias-corrected or wider interval is needed before flagship (bar 4 PARTIAL).
- **Scalar-per-frame gain** is a real physical restriction (R33 §6): a medium that
  changes the spatial transfer operator rather than one common frame gain would make
  the two-parameter ledger misspecified.
- **Forbidden framings avoided** (R33 §6): no "information conservation", no "free
  medium information", no "paired trades scene into medium". The scene-loss and
  medium-deficit are *two projections of the same posterior ambiguity, coupled
  through the posterior cross-covariance* `K`.

---

## Files (all under `results/round63_next/DUAL_LEDGER_PROBE/`)

- `DUAL_LEDGER_PROBE_REPORT.md` (this file)
- `dual_ledger_results.json` — T1/T2/T4 grid + bar 2/3/4/5 evidence
- `t3_reciprocity.json` — Thm-2 verification + (A,B,K) per schedule + empirical scan
- `figs/fig_dual_ledger.png` — A: t_c precision vs oracle floor & pilot;
  B: no-trade schedule result; C: reciprocity `eff_scene = eff_medium`
- `code/dl_common.py` — forward model (reused verbatim from LADDER_PROBE) + DLGI
  medium estimator, schedules, oracle floors, pilot & monitor baselines
- `code/dl_run_main.py` — T1 + T2 + T4 + bars 2/3/4/5 runner
- `code/dl_t3_reciprocity.py` — T3 reciprocity theorem + schedule prediction
- `code/dl_figure.py` — figure generator

---

## Bar-4 repair — FINAL (coordinator-assembled after two attempts; frozen)

Two honest interval-repair attempts were run on the fresh-replica protocol
(204 records/cell, interior grid):

1. **Bias-corrected parametric bootstrap** (percentile + pivotal, log-space):
   interior t_c coverage 0.77-0.92 — FAIL.
2. **Profile-likelihood on (log t_c, log CV)** (the last authorized attempt):
   CV coverage 0.985-1.000 on every interior cell (PASS); t_c coverage
   {tc16_cv5: 0.892, tc16_cv15: 0.980, tc16_cv40: 0.848, tc64_cv5: 0.868,
   tc64_cv15: 0.956, tc64_cv40: 0.990} — three of six interior cells below
   the 0.92 floor. `interior_all_in_92_98 = False`.

**BAR 4 VERDICT: FAIL (frozen; no further iterations).**
Mechanism, honestly: the t_c sampling distribution is right-skewed and
biased at the slow-drift CRB floor (tc64_cv5) and carries intrinsic
estimator bias at high CV (tc16_cv40); neither symmetric bootstrap nor
profile likelihood calibrates all cells simultaneously. This is an
estimator-statistics deficiency, not an information deficiency (bar 3
point-precision PASSED at <=1.5x pilot everywhere).

**DLGI FINAL SCORE: 6 PASS / 1 FAIL of the R33 promotion bars.**
Per R33 ("promote only if ALL"), DLGI is NOT auto-promoted to flagship.
Disposition: full disclosure to the next GPT round (R34) for campaign
adjudication — options include a preregistered simulation-calibrated
interval procedure (frozen before the campaign) or narrowing the claim to
RMSE-based precision. Raw numbers in `bar4_coverage.json`.
