# FLAGSHIP MATERIALS MAP — one-paper consolidation (GI_a2 / ROUND63)

**Purpose.** Single exhaustive catalog of every source the writing phase will
consume. INVENTORY ONLY — no manuscript prose here. Consumed by the operator +
the R35 architecture ruling to produce the frozen outline.

**Repo:** `D:\GI_another` (cwd = repo root). **Compiled:** 2026-07-22.
**No commits made by this map.**

**Ruled arc (R34 §5, the frozen story spine):**
- **Act I** — hidden-state identities + jitter-capped ridge law (quantify hidden-channel information)
- **Act II** — M1 certified campaign (global power-control policy, +1.87 dB / 19.13×)
- **Act III** — saturation boundary: five preregistered negatives, double-sided bracketing (the impossibility-both-sides argument)
- **Act IV** — DLGI (CONDITIONAL: leads METHOD act only if the R34 calibrated-interval campaign passes; else demoted to note/materials bank)
- **Venue:** Optica-first — lean main text + deep supplement.

**Provenance anchors referenced throughout:**
- Paper-1 Study-1 numbers = frozen verdict `results/round63_s2_detail/pilot_summary.json` @ commit `07ce17d`; Study-2 no-gate numbers = mode CSVs @ `b0a70ad`.
- Paper-2 M1 numbers = `results/round63_m1/CORRECTION_2026-07-19/M1_VERDICTS_SPEC_CORRECTED_R19.json` (R19-corrected).
- Act-III certificates: bridge @ `026c239`/`384315a`; DOPS @ `c592a4e`; Probe A @ `69207b3`; CPL @ `73b33ab`.
- GitHub rulings archive: `github.com/ccyyyYyzz/GI_a2`, issues #1–#26 = R-rounds; raw copies in `docs/ROUND63_GPT_ROUND*_RULING_RAW.md`.

---

## 1. PAPER-1 — `paper/main_oe.tex` (750 ln) + `paper/supplement.tex` (270 ln)

**Line:** "When high flux helps single-pixel imaging: a contrast–dead-time
operating map." Two-study operating-map paper. **Frozen title (R9):** "operating
map", NOT "phase diagram". Contains BOTH the Act-I ridge theory AND a
self-contained two-study campaign that is largely OE-derivative-specific.

### 1a. main_oe.tex section-by-section

| Sec (lines) | What it is | Key frozen numbers / wording | Destination |
|---|---|---|---|
| Abstract (37–70) | Two-study operating-map summary | Study-1 fixed-dwell **+0.28 dB**; Study-2 primary formally negative (16/24), fixed-dwell **+1.16 dB**; ridge $\rho^*(\nu)=(6\nu)^{1/3}-2/3$ | rewrite for flagship |
| Intro (72–141) | Field-first Chen-style; low-flux convention; ridge reframe; sparse prior-art; contributions | Para 3 + Para 4 "safe novelty sentence" are FROZEN verbatim (round-4/5/8) | Act I intro seed |
| Fig 1 mechanism (143–160) | SPI chain + dead-time + operating map schematic | $\rho^*\approx7.8$ ($\nu$=100), $22.2$ ($\nu$=2000); $\Gamma=1$ = descriptive proxy only | Act I hero candidate |
| §2 Physical model (162–191) | Nonparalyzable renewal counting, exact PMF Eq (1), $\rho_i=\lambda_i\tau$ | Exact PMF `P(N≥m)=P(Pois(z_m)≥m)`; deployment uses per-pattern $\rho_{95}$ not mean | Act I methods/supp |
| §3 RQL reconstruction (193–239) | Convex renewal quasi-likelihood Eq (2); truth-free reg. rule; descriptive audit | $\hat\lambda=N/(T-N\tau)$ = classical dead-time correction; $\lambda_{TV}=c\,\sigma_g\sqrt{2\ln n}$, **c=0.50 frozen**; ceiling $N\tau\ge T$ | Act I/II methods (RQL is shared estimator) |
| §4 Ridge theory (241–378) | **Act I CORE.** Missing-info identity Eq (3); cube-root law Eqs (4–6); deployment zoning; relation-to-prior-work | $I_N=E[N]-\rho^2 E[Var(R_\nu|N)]$; $\rho^*=(6\nu)^{1/3}-2/3+O(\nu^{-1/3})$; $J(\rho^*)=1-0.8255\nu^{-1/3}$; exact $\rho^*$ table {4.53,6.16,7.87,9.99,13.77,17.45,22.16} at $\nu$={20,50,100,200,500,1000,2000}; CLT trust $\rho_{0.9}\asymp\sqrt{1.2\nu}$ | **Act I main** (identity + law) |
| §4.4 Relation to prior work (327–378) | R12 frozen wording, verbatim w/ bibkeys | Müller/Yu-Fessler/Wang/Alvarez/Bécares/Grönberg fences; **Grönberg full text paywalled** (S5 hard sign-off item) | Act I supp / related work |
| §5 Study 1 dense (380–469) | Preregistered dense-Bernoulli study; formally negative primary + small positive secondary | 24/24 `safe_range_uninformative` (S_gate=1); fixed-dwell **+0.28 dB** median, 24/24 positive, LB **+0.23 dB**; natural-30 cohort median **+1.44 dB**, 29/30, LB **+1.18 dB**; classical GI 4.8–6.3 dB; $C_u\sim1/\sqrt n\approx1.6\%$ at $64^2$ | **mostly CUT** (OE two-study framing); mechanism ($C_u$) → Act I/III bridge |
| §6 Study 2 ladder (471–619) | Equal-load equal-dose occupancy ladder; contrast-conditioning trade-off | Primary formally negative (16/24 <18); median $S_{gate}$=**6.796**, LB **2.835**; fixed-dwell **+1.16 dB**, 19/24, LB +0.72; family contrast order maze .53→microtex .17; $C_u$={.040,.222,.316,1.270}, $\Gamma$={1.09,6.07,8.66,34.79}; $k$=32 median **+4.13 dB** | **CUT / materials bank** (superseded by M1 SCAT32); contrast mechanism reusable |
| §6.5 Mechanism separation (605–619) | Normalized vs unnormalized flux — not rate starvation | normalized +1.81/+1.10 dB ($\nu$=200/2000); unnormalized 7.45 vs 10.26 dB (worse absolute) | supp (mechanism defense) |
| §6.6 Scattering (621–629) | Dynamic scattering CV(α)=0.20 | +0.66/+1.84/+2.32 dB at $\nu$=20/200/2000; $S_{det}$ axis (Hao 0.096–0.44) | supp candidate |
| §6.7 Detector mismatch (631–648) | Descriptive robustness at anchor $M/n$=0.5,$\nu$=500 | all perturbations ≤0.34 dB; dead-time miscal ≤0.12 dB; dark 0.16/0.34 dB; afterpulse ≤0.09 dB | supp |
| §6.8 Exact-vs-RQL (650–669) | $8^2$/$16^2$ exact likelihood verification, 864 cells | median gap ≤0.01 dB in deployment zone, −0.05 dB extreme corner; individual cells ±9 dB | Act I/II supp (surrogate justification) |
| §7 Discussion+Limitations (671–737) | Mixed scientific position; 3 scope boundaries; 3 extensions; coherent dark-port extension | peak irradiance rises as $n/k$ (64× at gate occupancy); $C_u^2=\frac{n-k}{k(n-1)}(n\|x\|_2^2-1)$ | Act III discussion seed; cite companion |
| Availability/Ack (739–744) | \SPH author/repo/funding | **USER DECISION placeholders** | flagship front matter |

### 1b. supplement.tex sections (all tables machine-generated by `code/round63/figs/supp_s*.py`)

| Sec | Content | Key frozen number | Destination |
|---|---|---|---|
| S1 (47–117) | Gaussian-likelihood non-coercivity derivation | log-det sign reverses at $\rho=1/2$; $\mu=\nu\rho/(1+\rho)$, $\sigma^2=\nu\rho/(1+\rho)^3$; unbounded-below at ceiling | Act I supp (RQL justification) — KEEP |
| S2 (120–142) | Full NATURAL-30 cohort table | low_contrast $\Delta Q$=−1.78 dB; fine_lines +0.07 dB | cut/materials |
| S3 (145–164) | Detector-mismatch full tables | dark50%-known PRECORRECT un-flooring −2.60 dB anomaly | supp (thin) |
| S4 (167–180) | Exact-vs-RQL 12×4 grid, 864 pairs | up to 9.23 dB per-cell both directions | Act I supp |
| S5 (183–203) | No-gate occupancy ladder | medians +0.64/+4.13/+1.16/+3.56 dB (k=512/32/16/1) | cut/materials |
| S6 (206–227) | Dynamic-scattering detail | +0.66/+1.84/+2.32 dB | supp candidate |
| S7 (230–248) | Ensemble/undersampling S2B/S2C | separately-provenanced regenerated artifact | cut |
| S8 (251–268) | Descriptive measurement-audit distributions | $d_{ratio}$ quantiles; not calibrated p-values | cut |

**Verdict on paper-1:** §4 (ridge theory) + S1 = **Act I backbone, reuse
heavily**. The two-study campaign (§5, §6) is **OE-derivative-specific** — the M1
campaign (paper-2) is the stronger, certified campaign, so Study-1/Study-2 are
CUT to the materials bank except: the missing-info identity, cube-root law,
$C_u$ contrast mechanism, and the exact-vs-RQL surrogate justification.

---

## 2. PAPER-2 — `paper2/main_m1.tex` (450 ln) + `paper2/supplement_m1.tex` (770 ln)

**Line:** "Jitter-capped high-flux single-pixel imaging with global power
control." This is the **Act I law extension (jitter cap) + Act II campaign**.
Restructured R20/R21 to one storyline; R19-corrected numbers.

### 2a. main_m1.tex section-by-section

| Sec (lines) | What it is | Key frozen numbers / wording | Destination |
|---|---|---|---|
| Abstract (29–59) | Jitter cap + campaign summary | $c_v^2\rho^2(1+2\rho)=1$; cap scales $c_v^{-2/3}$; **+1.87 dB** 19/24; **37.1×** dose, **2.6×** counts; speedup **19.13×** LB **18.33×** 21/24 | Act I+II abstract seed |
| §1 Intro (61–97) | Deterministic→random recovery; $c_v$=0.05 moves optimum 22.3→5.7, forfeits ~55% info | power-for-time framing; count not timestamps; simulation-only | Act I/II intro |
| Fig 1 hero (99–116) | Jitter-cap mechanism; det vs jitter $J(\bar\rho)$ at $\nu$=2000 | $\bar\rho^*\approx22$ (det) → $\approx5.7$ ($c_v$=0.05); ~55% forfeit | **Act I hero candidate** |
| §2 Random-hold info (118–167) | **Thm 1** missing-info identity; **Thm 2** long-window rate | `I_N=E[N_T]−λ²E[Var(L_T|N_T)]`; `J_∞=ρ/[(1+ρ)(1+c_v²ρ²)]` | **Act I main** |
| §3 Jitter-capped load (169–230) | **Thm 3** exact optimum; small-jitter expansion; finite-window crossover; MC verification | `c_v²ρ_c²(1+2ρ_c)=1`; `ρ_c=2^{-1/3}c_v^{-2/3}−1/6+...`; crossover `ρ_*∼(6ν)^{1/3}(1+12νc_v²)^{-1/3}`; pooled slope **−0.658** vs −2/3 (1.3%) | **Act I main** |
| §4 Power-control policy (232–262) | Global multiplier policy; $\rho_R(2000)$=22.25; resource interpretation | one global source multiplier, no per-pattern servo; power-for-time explicit | **Act II main** |
| §5 Campaign (264–311) | 24 scenes × 6 families × 5 seeds × 9 dwells; SCAT32; arms; endpoints; **R19 correction subsec** | primary `RIDGE_OPERATING_PASS`; ≥1 dB / LB>0 / ≥18-24; elapsed $S_{gate}$≥3/LB>1/≥18; geometry 0/299/181 no PASS/FAIL | **Act II main** |
| §5 Analysis correction (294–311) | R19 peer-audit correction record | 3 fixes: axis $\log(\nu\rho)$→$\log T_{opt}$; family-strata bootstrap; PAVA+censoring; 18/18 agreed; bundle `CORRECTION_2026-07-19/` | **Act II provenance (FROZEN wording)** |
| §6 Imaging results (313–381) | Fixed-dwell + elapsed-time pass; resource table | **+1.87 dB** LB **1.41** 19/24; **37.1×** dose **2.6×** counts; **19.13×** LB **18.33** 21/24; family +6.37 (maze) to −4.18 dB (contour); contrast corr +0.64; post-hoc incident ratio 0.279 no verdict | **Act II main** (Table 1 = operating-results) |
| §7 Geometry headroom (383–399) | Descriptive scene-adapted geometry | dose-uniform 12/12 D-eff LB 1.64–2.51×; full-cond 0/480 certified, 299 counterexample, 181 unresolved | Act II supp / next-paper seed |
| §8 Discussion+Limitations (401–438) | Detector-specific rule; open geometry axis; 3 limitations; immediate bench test | 37.1× peak-power limit; iid-hold assumption; sim-only; bench: $c_v$=0.05→peak $\rho$=5.7 | Act II discussion |

### 2b. supplement_m1.tex sections

| Sec (lines) | Content | Frozen number/wording | Destination |
|---|---|---|---|
| S1 (47–75) | Framework theorem (R14 Cand A generalized-KKT); discrete-optima caveat verbatim | information-water-filling; $k\approx32$/$p\approx1$ motivate NOT KKT theorems | Act II supp |
| S2 (78–159) | Alignment proposition (R14 Cand C); projected-KKT residual bound Eq (S-residual) | $0\le F(\rho^*)-F(\tilde\rho)\le\frac1{2m}\sum\pi_i r_i^2$; uniform servo loses ~7 dB sparse | Act II supp |
| S3 (162–425) | Full design machinery: r=200 subspace, full-stack certificate $\mathcal C_{stack}$, guards-as-constraints, frozen dictionary, budget accounting, dose-uniform construction | $\epsilon_0=10^{-9}$; Eff$_D\ge e^{-0.01}$=0.99005; $D_{cert}$ SHA `908cfccbd249de22`; dev-replication table (12 cells, D-eff 1.64–2.51×) | Act II deep supp |
| S4 (428–476) | Monte-Carlo verification detail | `jitter_sfi_v2` 100k frames; peaks {22.297,10.739,5.173,3.862,2.153,1.607}; pooled slope **−0.658**; score-identity resolves +8.8% histogram bias; k=4/3 refuted | Act I supp (verification) |
| S5 (479–540) | Campaign accounting + R19 correction record table | 1024 exposures=52 pre-scan+972 main; Table S-correction (3 discrepancies) | Act II supp (provenance) |
| S6 (543–665) | Secondary-arm + certificate full results; per-dwell table; context arms | per-dwell $\Delta Q$ table (Table S-perdwell); SCAT32-060 17.12 > SCAT16 14.18 > LBLOB16 9.63 dB; 0/299/181 dist. | Act II supp |
| S7 Proofs (667–765) | Thm 1–3 proofs; uniform crossover prop; alignment proof | martingale score + total variance; delayed-renewal LAN | Act I supp (proofs) |

**Verdict on paper-2:** This is the **strongest existing draft** — its Act-I
theorems (jitter cap) and Act-II campaign are the flagship's Acts I and II
almost as-is. The R19 correction subsection wording is FROZEN (see Register).

---

## 3. ACT III — SATURATION BOUNDARY (five certificates → one boundary table)

Master synthesis: `docs/SOFTWARE_SATURATION_VERDICT.md` (82 ln, R29 §6 frozen).
**Program-level frozen conclusion (verbatim):** *"No current software candidate
clears the operator's novelty × simplicity × generality × image-positive bar.
The software layer for these frozen channels is saturated at the requested effect
size (≥1 dB deployable image gain)."*
**Why (identity's own explanation, verbatim):** *"Information lost = E[Var(hidden
state | record)]; the bucket-count record is too thin to identify the hidden
state, so every software post-hoc trick is bounded by that conditional variance —
and five independent mechanisms measured the bound as binding."*

### 3a. THE SINGLE ACT-III BOUNDARY TABLE (one line per certificate)

| # | Lane | Mechanism tested | Verdict | Governing number | Bar | Source files |
|---|---|---|---|---|---|---|
| 1 | RLMI bridge (Gate A) | scene-adaptive geometry from noisy pre-scan, dose-uniform library | **FAIL** | median ΔQ^A **+0.680 dB** < 1.0 (11/12 pos, boot LB +0.543) | ≥1 dB | `results/round63_bridge/BRIDGE_VERDICT.md` + `GATE_VERDICT.json` @026c239/384315a |
| 2 | RLMI allocator (Gate B/C) | robust maximin routing on noisy pre-scans | **HARM** | RLMI median **−7.45 dB** vs base; allocation loss **−8.33 dB**; controls to −26.1 dB | no-harm | same (Gate B/C) |
| 3 | DOPS scheduling | O4-A moment-orthogonal pattern order (drift) | **FAIL** | median dQ **+0.039 dB** < 1.0 (5/6 pos) @c592a4e | ≥1 dB median | `results/round63_next/SCHEDULE_PROBE/DEV_GATE_VERDICT.json` + `PROBE_B_REPORT.md` |
| 4 | Estimator efficiency (Probe A) | RQL vs certified jitter-capped Fisher ceiling | **SATURATED** | ~98.9–99.5% var-efficient; ≤0.2 dB range headroom @69207b3 | headroom | `results/round63_next/EFFICIENCY_PROBE/PROBE_A_REPORT.md` |
| 5 | CPL-GI (reserve) | exact conditional-likelihood gain elimination | **KILL** | median dQ **+0.412 dB** < 1.0 (5/6 pos) @73b33ab | ≥1 dB | `results/round63_next/CPL_PROBE/CPL_GATE_VERDICT.json` + `CPL_PROBE_REPORT.md` |

**Supporting descriptive facts (from master doc + bridge q4):** dose-relaxed FW
oracle sits **−12.91 dB** below simple designs at image level (surrogate-image
gap); known-gain oracle sits **10–15 dB above ALL software arms** (headroom
exists but reachable only by KNOWING hidden state = richer record).

### 3b. Per-certificate detail (for the supplement / methods appendix)

- **Bridge (cert 1/2):** grid 320/320 cell-groups, 1920 arm records, 117 min, zero skips. Decision cell (c=0.05,ν=2000), 12 stress scenes, 5 seeds. k* always ∈ {L0 (SCAT32 knob), L1 (scattered-k32)} — structured banks never win. Four-quantity decomposition: (1) +0.680, (2) −8.33, (3) −1.73 (plug-in FW), (4) −12.91 dB. Verdict = LIBRARY_REACHABILITY_FAIL → RLMI/M2 STOP. **M2: DEAD (never preregistered/run).**
- **DOPS (cert 3):** OU tc=2, 4 scenes × 3 drifts × 5 seeds. Fisher-side paired-over-random 63.6→0.0 does NOT survive to image (+0.10 dB). Survivor by-product: **anti-naive-ordering +1.34 dB pooled / +1.96 dB hard cell** (= "interleave, randomly is enough" — prior-art-adjacent). Dominant lever: correlation time vs acquisition length (tc=2 « N=1944).
- **Probe A (cert 4):** ν=2000 data term efficient at every load; only real gap = ν=200/ρ=22.25 (0.62 dB var / 1.7 dB MSE) but present at c=0 → finite-sample bias of `N/(T−Nτ)`, fixed by bias-corrected/exact-MLE, not DL. Band-4 (ill-conditioned) residual = prior-improvement target, uncertified, NOT Fisher-ceiling target. Pipeline validated to <0.02% vs exact analytic Fisher.
- **CPL (cert 5):** raw λ=0 margin +1.585 dB evaporates to +0.412 under equal-terms TV prior; reverses on feasible-imaging cells (tc=64: −0.224; σ15: −0.296); oracle 10–15 dB above. Kill-risk (R29 §2.1 "statistically dressed differential Hadamard normalization") realized.

### 3c. Ladder-probe scope amendment (adjudicated 2026-07-22, in master doc §Adjudicated amendment)

- Record refinement per se: WEAK (~1 dB, saturates B≈16); VACUOUS under within-exposure-constant gain; HARMFUL at fast drift (tc=2). Ch.1 detection-timestamps recover only 11–15% of jitter loss (count pins summed holds; Fano≈0.004 at ridge) — observe-route REFUTED for detection-only tagging.
- FLAG upheld: temporal-prior B=1 BUCKET estimator recovers **+2.2–2.5 dB** over killed A2/A3 arms; ADMISSIBLE only as ch.2 DATA TERM of certified-imaging direction, NOT a standalone novel method.

**Act-III scientific reading (bridge doc, for closing act):** *"neither software
post-correction nor safe-class illumination adaptation clears the deployment bar
— the certified global operating law (ridge tracking, +1.87 dB / 19.13×, M1) is
the practical optimum of this regime, now bracketed from both sides by
preregistered evidence."*

---

## 4. ACT IV — DLGI (Dual-Ledger Ghost Imaging) — CONDITIONAL

Method name **DLGI**; theorem name **Canonical-Confusion Ledger Reciprocity
Theorem** (must stay distinct — R33 §7). Backup mainline activated after MOLT's
R32 downgrade.

### 4a. Feasibility evidence

`results/round63_next/DUAL_LEDGER_PROBE/DUAL_LEDGER_PROBE_REPORT.md` +
`dual_ledger_results.json` + `t3_reciprocity.json` + `bar4_coverage.json` +
`figs/fig_dual_ledger.png`. Governing theory = R33; campaign ruling = R34.

**Seven-bar verdict (R33 §6): 6 PASS / 1 FAIL** (final, frozen — bar 4 intervals FAIL after two repair attempts).

| # | Bar | Result | Verdict |
|---|---|---|---|
| 1 | Identical photons/exposures/multiset; only order+inference change | byte-identical record; medium = post-processing | PASS |
| 2 | Scalar-gain model + gauge pass held-out residual tests | innovation std 0.89–0.97; lag-1 ≤0.11; gauge shift 3×10⁻¹⁵ | PASS |
| 3 | Pilot-free RMSE ≤1.5× pilot both params, near oracle | tc ratio ≤0.92, CV ≤1.15 | PASS |
| 4 | Calibrated interval coverage | tc 0.73–0.90 (profile-LR: 3/6 interior <0.92) | **FAIL** |
| 5 | Scene noninferior ≤0.2 dB vs plain linear | Δ −0.04 to +1.76 dB; superior 7/9 | PASS |
| 6 | Schedule behavior predicted by (A,B,K) | reciprocity det exact to 1e-15; ‖K‖ schedule-invariant 0.2%; paired best for BOTH | PASS |
| 7 | Edge regimes honest failure map | fast tc=2 (+40%), slow tc=64 (25% floor), CV-below-noise mapped | PASS |

**T1 grid (pilot-free medium precision, tc×CV):** at every tc≥16 cell, max |tc| err **13.5%**, max |CV| err **5.5%** (bar ≤20% MET). Fast edge tc=2/CV40 biases tc **+40%**; slow tc=64 at **25%** oracle floor.
**T3 no-forced-trade:** paired schedule best for BOTH ledgers (paired ‖K‖=0.3395, scene PSNR 16.5/17.6, tc err 5%/6%) > random > split. Reciprocity `det(J_x)/det(A)=det(J_θ)/det(B)` exact 3×10⁻¹⁵.
**Bar-4 mechanism (frozen):** tc sampling distribution right-skewed/biased at slow-drift CRB floor + intrinsic bias at high CV; estimator-statistics deficiency NOT information deficiency.

### 4b. R33 theorem statements (FROZEN — mark verbatim; see Register)

- **Theorem 1** (joint two-ledger missing-info identity): `I_Y(x,θ)=diag(C_x,C_θ)−L`, blockwise `I_xx=C_x−L_xx`, `I_θθ=C_θ−L_θθ`, `I_xθ=−L_xθ`, `L≥0`.
- **Theorem 2** (canonical-confusion ledger reciprocity): with `K=A^{-1/2}CB^{-1/2}`, `J_x=A^{1/2}(I−KK^T)A^{1/2}`, `J_θ=B^{1/2}(I−K^TK)B^{1/2}`, and `det(J_x)/det(A)=∏(1−κ_i²)=det(J_θ)/det(B)`.
- **What O4-A does (§3):** moment condition `M_π^T R^{-1}D_sH=0` ⟺ `C=0` ⟺ `K=0` → both ledgers attain oracle Fisher. It does NOT blind the medium.
- **Medium identifiability (§4.1):** `J_θ=I_θθ−I_θx I_xx^† I_xθ > 0`. Oracle floors: `sd(t_c)/t_c ≥ √(2t_c/T)`, `sd(CV)/CV ≥ √(t_c/2T)`.

### 4c. R33 FORBIDDEN claims (must never appear)

- No "information conservation"; no "free medium information"; no "paired schedules trade scene information into medium".
- No "scene information + medium information = constant" (ledgers have different coordinates/dimensions/units).
- FORBIDDEN novelty firsts (prior-art fence): first info partition/nuisance-orthogonal design; first single-pixel/optical medium-dynamics or t_c measurement (DCS/DWS own it); first calibration-as-science-product (radio selfcal); first telemetry medium ID (AO); first blind joint scene/gain recovery.

### 4d. R34 authorized claim wording (until campaign passes — FROZEN, see Register)

> *"DLGI extracts a competitive point estimate of medium correlation time and fluctuation strength from the same bucket record used for scene reconstruction, with no added acquisition and no detected scene penalty; model-calibrated t_c uncertainty remains under confirmation."*
> Do NOT use "two certified products" yet.

**R34 campaign bars (C1–C7, binding kill bars if the confirmatory campaign runs):**
- Freeze simulation-calibrated Neyman construction on η=(log t_c, log CV) BEFORE confirmation; repair allowed ONCE.
- Grid: t_c∈{16,32,64} × CV∈{0.05,0.15,0.40} × photon∈{0.5,1,2}× = 27 cells; ≥5000 calibration + ≥1000 confirmatory records/cell; ≥12 fresh scenes.
- C1 coverage ∈[0.92,0.98] every cell (named stress: (16,.05),(16,.40),(64,.05)); C2 width ≤1.5× pilot; C3 RMSE ≤1.5× pilot; C4 scene noninferiority ≤0.2 dB & ≤5% Fisher; C5 model/gauge; C6 reciprocity+schedule; C7 identical-ledger + edge honesty.
- If all pass, restore R33 headline ("two model-certified products"); "model-certified" + validated scalar-gain domain MANDATORY.

**R34 §5 placement (FROZEN):** DLGI **conditionally leads the flagship METHOD
act**; do NOT split into standalone short paper before confirmation. Passing →
METHOD act + "two model-certified products"; failing → shorter theorem/method
note or materials-bank section.

**GAP:** the R34 confirmatory campaign does NOT yet exist — Act IV is contingent.

---

## 5. MATERIALS BANK / SUPPLEMENT CANDIDATES

### 5a. MOLT (Microcell-Occupancy Laplace Tomography) — RETIRED, materials-bank only

- **Theory (R31 `docs/ROUND63_GPT_ROUND31_RULING_RAW.md`, 982 ln + PRO appendix `..._PRO_APPENDIX_0.md`, 203 ln):** Masked Laplace–Power-Sum Theorem. `r(ρ;v)=∏(1+ρv_p)^{-1}`, `log r=Σ(-1)^j ρ^j p_j/j`; second channel `p_2=M(x∘x)`; saturation ridge `t*=2.821439`; Jacobian rank `min(N,2K)` doubles local rank 51→102. PRO adds: correlated-speckle determinant law `q=det(I+(ρ/k)ΣV)^{-k}`, `n_eff²` photon-dose law `B_{3σ}≳25.33 n_eff²/f²` (known-p1) / `47.70 n_eff²/f²` (nuisance), quenched floor `SE(β_2)/β_2≈(2/t)√(n_eff/C_eff)`.
- **Probe (`results/round63_next/SATURATION_JAILBREAK_PROBE/JAILBREAK_PROBE_REPORT.md`):** R32 six-bar = 4 PASS / 2 CONDITIONAL. Frozen thesis: "sparse MOLT trades photon dose for sensing rank." Noiseless witness s=80: joint 296.7 dB vs linear 27.8 dB (+268.9 dB); matched-photon fails by orders of magnitude (dense ×1.6e5–1.1e6; sparse ×126).
- **Revival calc (`results/round63_next/MOLT_REVIVAL_CALC/REVIVAL_CALC.md`):** **STAY DEAD.** Best honest-dose margin 4.13× (need 10×); binary trick `x²=x` kills classes a/b; DL-tangent route dies (linear closes d≤102). Decision: do NOT open GPT round; MOLT retired.
- **Destination:** ONE-PARAGRAPH materials-bank mention at most (the R32-authorized narrow claim), OR cut. Its DOIs feed the defensive bibliography (§8).

### 5b. Anomaly probe (`results/round63_next/ANOMALY_PROBE/ANOMALY_PROBE_REPORT.md`)

- Range–null accountability / anomaly-fidelity: A3 (info-guided, range-locked) recovers exactly witnessed fraction (slope 1.09/0.97); A2 (unconstrained DL) **erases ~30%** of measured anomaly content (survival A3 83% vs A2 53% VQAE / 28% VQGAN). Sub-resolution: witnessed fraction flags with **AUC=1.000**. Hallucination check: at 5% sampling failure mode is **erasure NOT invention** (0 false blobs A2). Governed by R30 (certified information-guided imaging fence).
- **Destination:** materials bank / next-paper (accountability line, R30). Real operator/prior (lane0 rate05, VQAE/VQGAN). Not part of the frozen 4-act arc unless the operator forks to the audit line.

### 5c. Ladder probe (`results/round63_next/LADDER_PROBE/LADDER_PROBE_REPORT.md`)

- **Ch.2 temporal-prior result:** B=1 state-space smoother recovers **+2.2–2.5 dB** over A2/A3 on same bucket record (assumes known OU stats). This is the estimator that becomes DLGI's ch.2 DATA-TERM lineage (feeds Act IV estimator provenance).
- **Ch.1 timestamp refutation:** detection-timestamps recover only 11–15% of jitter loss (Act-I aside — "the count already pins the summed holds"). Nice Act-I supplement footnote.
- **Destination:** Act IV estimator lineage (ch.2) + Act I aside (ch.1 timestamp refutation, supp).

---

## 6. FIGURES INVENTORY

### 6a. paper/figs (paper-1; owned generators `code/round63/figs/*.py`; README `FIG_MECHANISM_P1_README.md`)

| File | Thumbnail description | Reuse candidate |
|---|---|---|
| `fig_mechanism_p1.pdf` | Fig 1: SPI chain + dead-time/ridge + Γ=1 schematic (neutral) | **Act I hero** (mechanism) |
| `fig_a_ridge_map.pdf` (in `results/round63_theory/`) | exact count-Fisher info per unit time vs load, peak on ρ*(ν) | **Act I law figure** |
| `fig_b_masks.pdf` | representative masks k∈{512,32,16,1} | Act III / materials |
| `fig_c_ladder.pdf` | measured $C_u$ + Γ vs occupancy | Act III mechanism |
| `fig_d_dwell.pdf` | radiometric PSNR + CNR vs dwell | cut/materials |
| `fig_e_gate.pdf` | per-image fixed-dwell gain + per-scene speed | cut/materials |
| `fig_f_scattering.pdf` | dynamic scattering vs $S_{det}$ | supp candidate |
| `fig_p1_families.pdf` | six-family empirical acceleration | Act III / materials |
| `fig_s1b/s1c/s1d`, `fig_s2_gallery`, `fig_nat_pair(_16)` | Study-1 images/curves/gallery, natural pair | cut/supp |

Style (M10/SF7/SF8): sans-serif, neutral gray=safe, blue=ridge/high-flux, vermillion=dead-time/negative. Sizes approximate; Optica resize pending (M9).

### 6b. paper2/figs (paper-2; generator `code/round63/figs/fig_paper2_r20.py`; README `ACTIII_README.md`)

| File | Thumbnail description | Reuse candidate |
|---|---|---|
| `fig_mechanism.pdf` | Fig 1 hero: chain + det-vs-jitter + jitter cap | **Act I hero (jitter)** |
| `fig_jitter_validation.pdf` | Fig 2: peak-vs-$c_v$ (slope −0.658) + retention warning | **Act I verification** |
| `fig_speed_results.pdf` | recon + fixed-dwell ΔQ + elapsed speed (2-resource) | **Act II results** |
| `fig_audit_supp.pdf` | SUPP: guard path + certificate status (0/299/181) | Act II supp |
| `actiii_b.pdf` | SUPP: dose-only headroom (D-eff 1.6–2.5×) | Act II supp |
| `actiii_a/c/d/e.pdf` | RETIRED/CUT (merged/folded); kept for provenance | — |

### 6c. Probe figures (results/round63_next/*/figs and roots)

| File | Description | Reuse |
|---|---|---|
| `DUAL_LEDGER_PROBE/figs/fig_dual_ledger.png` | A: tc precision vs oracle/pilot; B: no-trade schedule; C: reciprocity eff_scene=eff_medium | **Act IV hero candidate** |
| `SATURATION_JAILBREAK_PROBE/figs/fig_dprime_budget.png`, `fig_p0_mechanism.png`, `fig_p3_rank_doubling.png`, `fig_quenched_floor.png`, `fig_timing_overhead.png`, `fig_t3_scenes/x/xprime/delta.png` | MOLT mechanism/dose/rank/timing panels | materials bank |
| `ANOMALY_PROBE/anomaly_fidelity_3panel.png` | 3-panel anomaly fidelity | next-paper (audit line) |
| `SCHEDULE_PROBE`, `CPL_PROBE` | (no figs dir — JSON/report only) | Act III table only |

**GAP:** No flagship hero figure is designed. Existing heroes are per-paper
(paper-1 mechanism, paper-2 jitter, DLGI 3-panel). A unified 4-act arc figure
does not exist.

---

## 7. GPT ROUND LEDGER (audit-trail citation — "all rulings archived")

All at `github.com/ccyyyYyzz/GI_a2` issues #1–#26; raw in `docs/ROUND63_GPT_ROUND*`.

| Round | Subject (one line) | File |
|---|---|---|
| R3 | 总裁决 — outer AUDIT split + coherent refit bootstrap (digest) | `..._ROUND3_DIGEST.md` |
| R4 | Q1 outer AUDIT split + coherent refit bootstrap (digest) | `..._ROUND4_DIGEST.md` |
| R5 | no binary adequacy gate | `..._ROUND5_RULING.md` |
| R6 | summary-table ruling | `..._ROUND6_RULING.md` |
| R7 | decision (verbatim opening) | `..._ROUND7_RULING.md` |
| R8 | study2-freeze framing/intent | `..._ROUND8_RULING.md` |
| R9 | Final ruling (operating-map title, supplement structure S1–S8) | `..._ROUND9_RULING_RAW.md` |
| R10 | dead-time-aware optimal experiment design | `..._ROUND10_RULING_RAW.md` |
| R11 | ridge-zone load architecture | `..._ROUND11_RULING_RAW.md` |
| R12 | Grönberg check + X.5 frozen wording | `..._ROUND12_RULING_RAW.md` |
| R13 | Final M1 freeze audit (r=200, guards) | `..._ROUND13_RULING_RAW.md` |
| R14 | unification architecture (Cand A framework thm, Cand C alignment) | `..._ROUND14_RULING_RAW.md` |
| R15 | certificate amendment + freeze authorization | `..._ROUND15_RULING_RAW.md` |
| R16 | higher-order constant corrections (2M-frame diagnostic) | `..._ROUND16_RULING_RAW.md` |
| R17 | OED-DT retirement / SCAT32 amendment | `..._ROUND17_RULING_RAW.md` |
| R18 | full-stack near-optimality certificate | `..._ROUND18_RULING_RAW.md` |
| R19 | **peer-audit correction adjudication** (both verdicts PASS; +1.87 dB / 19.13×) | `..._ROUND19_RULING_RAW.md` |
| R20 | reader-side presentation audit + revision plan (figs restructure) | `..._ROUND20_REVIEW_RAW.md` |
| R21 | abstract relocation + correction-disclosure placement | `..._ROUND21_RULING_RAW.md` |
| R22 | controlled score-transfer object; limits of unification | `..._ROUND22_RULING_RAW.md` |
| R23 | rank-constrained control for hidden-state photon counting (bridge) | `..._ROUND23_RULING_RAW.md` |
| R24 | deployable residual-gated illumination router | `..._ROUND24_RULING_RAW.md` |
| R25 | standardized maximin bank allocation | `..._ROUND25_RULING_RAW.md` |
| R26 | information limits + temporal design for dynamic-scattering GI | `..._ROUND26_RULING_RAW.md` |
| R27 | dose-feasible bank rebuild + latency claim split (PLUGIN_LATENCY) | `..._ROUND27_RULING_RAW.md` |
| R28 | finite-library Gate A + dose-relaxed FW diagnostics | `..._ROUND28_RULING_RAW.md` |
| R29 | **software-method triage** (DOPS campaign; saturation-verdict basis) | `..._ROUND29_RULING_RAW.md` |
| R30 | certified information-ledger imaging fence + DEV gates (anomaly probe) | `..._ROUND30_RULING_RAW.md` |
| R31 | **MOLT** conditional GO (+PRO appendix binding amendments) | `..._ROUND31_RULING_RAW.md` |
| R32 | MOLT → time-limited dose-expensive sparse-support method | `..._ROUND32_RULING_RAW.md` |
| R33 | **DLGI** exact Fisher partition + reciprocity theorem + flagship ruling | `..._ROUND33_RULING_RAW.md` |
| R34 | **DLGI final campaign GO** (calibrated intervals; METHOD-act lead) | `..._ROUND34_RULING_RAW.md` |
| R35 | (question posed — architecture ruling PENDING) | `..._ROUND35_QUESTION.md` |

---

## 8. DEFENSIVE CITATION SEED (deduplicated; tag = claim it fences)

### 8a. Already in the manuscripts (paper/refs.bib ∪ paper2/refs.bib)

GI/SPI core: `duarte2008singlepixel`, `shapiro2008computational`, `edgar2019principles`, `harwit1979hadamard`, `wuttig2005multiplex`.
Photon-limited/scattering (Hao): `hao2025photonlimited`, `hao2025randommedia`, `hao2025complexfield`.
Dead-time/pileup/lidar: `rapp2019deadtime`, `rapp2021highflux`, `gupta2019photonflooded`, `gupta2019asynchronous`, `liu2021photoncounting`, `jorgensen2026deadtime`, `wang2011pileup`, `becares2012deadtime`, `gronberg2018deadtime`, `muller1973deadtime`, `yu2000mean`, `alvarez2014pileup`.
Sparse-SPI pileup: `liu2019pileup`, `han2025sparsehadamard`.
Optimization/OED: `beck2009fast`, `rudin1992nonlinear`, `kiefer1959optimum`, `white1973equivalence`, `jaggi2013frankwolfe`, `pukelsheim1992rounding`, `nikolov2022volume`, `sagnol2015exact`.
Adaptive/sensing (paper2 only): `ariascastro2013adaptive`, `attia2018goaloriented`, `burger2021sequential`, `feng2018adaptive`, `maxu2021sensing`.
Cross-ref: `companion` (each paper cites the other; MERGE means this key is removed on consolidation).

### 8b. NEW DOIs from R31/R31-PRO fences (MOLT — needed only if MOLT keeps its paragraph)

| DOI | Work | Fences |
|---|---|---|
| 10.1103/PhysRevA.70.055801 | Rossi–Olivares–Paris, photon stats without counting | "first efficiency-sweep generating-function" |
| 10.1103/PhysRevLett.95.063602 | Zambra et al. reconstruction w/o counting | same |
| 10.1103/PhysRevA.85.023820 | Sperling–Vogel–Agarwal, multi on-off click stats | on/off click statistics |
| 10.1364/AO.57.006750 | Miatto–Safari–Boyd, on/off number discrimination | on/off discrimination |
| 10.1103/pknl-24xd | Kovtoniuk et al., nonclassical single on-off | on/off |
| 10.1109/TNS.2010.2053048 | van Dam et al., SiPM response model | SiPM saturation law |
| 10.1016/j.nima.2013.11.013 | Gruber et al., SiPM over-saturation | over-saturation |
| 10.1016/j.nima.2018.10.074 | Weitzel et al., SiPM single-photon→saturation | saturation |
| 10.1109/TNS.2021.3049675 | Kumar et al., digital SiPM inhomogeneous light | nonuniform illumination |
| 10.3390/s24082648 | Moya-Zamanillo–Rosado, SiPM nonlinearity | nonlinearity |
| 10.1109/CISS.2008.4558487 | Boufounos–Baraniuk, 1-bit CS | one-bit/quantized CS |
| 10.1002/cpa.21442 | Plan–Vershynin, 1-bit CS by LP | one-bit CS |
| 10.1016/0022-0000(85)90041-8 | Flajolet–Martin, probabilistic counting | occupancy sketch |
| 10.1145/78922.78925 | Whang et al., linear-time counting | frequency-moment sketch |
| 10.1006/jcss.1997.1545 | Alon–Matias–Szegedy, frequency moments | frequency moments |
| 10.1016/j.acha.2011.02.002 | Laska et al., quantization/saturation/CS | saturation-aware recovery |
| 10.1103/PhysRevLett.93.093602 | Gatti et al., thermal-light GI | higher-order GI |
| 10.1103/PhysRevA.78.061802 | Shapiro, computational GI | GI (already ≈shapiro2008) |
| 10.1364/OL.34.003343 | Chan–O'Sullivan–Boyd, high-order thermal GI | higher-order GI |
| 10.1364/OL.35.001166 | Chen et al., high-order lensless GI | higher-order GI |
| 10.1021/acsphotonics.8b00653 | Olivieri et al., time-resolved nonlinear GI | "nonlinear GI" name collision |
| 10.1364/JOSA.66.001145 | Goodman, speckle fundamentals | coherent-vs-incoherent addition (MANDATORY optics cite) |

### 8c. NEW DOIs from R33 fence (DLGI — needed for Act IV)

| DOI | Work | Fences |
|---|---|---|
| 10.1111/j.2517-6161.1982.tb01203.x | Louis, observed info in EM | missing-info partition (no "first") |
| 10.1111/j.2517-6161.1987.tb01422.x | Cox–Reid, parameter orthogonality | nuisance orthogonality (no "first") |
| 10.1103/PhysRevLett.60.1134 | Pine et al., diffusing-wave spectroscopy | first optical t_c / medium dynamics |
| 10.1103/PhysRevLett.75.1855 | Boas–Campbell–Yodh, diffusing temporal field correlations | DCS |
| 10.1364/OL.29.001766 | Durduran et al., diffuse optical blood flow | DCS |
| 10.1063/1.2037987 | Bandyopadhyay et al., speckle-visibility spectroscopy | speckle as medium sensor |
| 10.1364/AO.57.006097 | Yang et al., interpolated monitoring GI | pilot baseline (mandatory strongest baseline) |
| 10.1364/OL.463897 | Xiao–Zhou–Chen, temporal correction through scattering | dynamic-medium GI correction |
| 10.1364/OL.499787 | Peng–Chen, learning-based Gaussian-constraint correction | learning-based correction |
| 10.1364/OE.489808 | Zhou–Xiao–Chen, self-corrected SPI through scattering | dual-detector self-correction |
| 10.1051/0004-6361/200811094 | Intema et al., SPAM ionospheric calibration | calibration-as-science (no "first") |
| 10.1002/2016RS006028 | Mevius et al., LOFAR ionosphere | radio selfcal science |
| 10.1002/2015GL063699 | Loi et al., density ducts imaging | radio selfcal |
| 10.1093/mnras/sty3285 | Laidlaw et al., AO turbulence profiling | AO telemetry medium ID |
| 10.1093/mnras/stz3062 | Laidlaw et al., AO wind profiling | AO telemetry |
| 10.1137/16M1103634 | Ling–Strohmer, self-calibration/bilinear | blind joint scene/gain (no "first") |
| 10.1093/imaiai/iay004 | Cambareri–Jacques, blind gain calibration | pilot-free bilinear calibration |

---

## 9. FROZEN-WORDING REGISTER (may NOT be paraphrased — reuse verbatim)

### 9a. R19 verdict sentences (paper-2 §5–6, §S5)
- Fixed-dwell: "ridge control improved median fixed-dwell radiometric quality over the matched $0.60$ comparator by **+1.87 dB**, with a family-stratified lower bound of **1.41 dB** and **19 of 24** scenes positive."
- Resource: "The same operating point used **37.1×** the incident dose and **2.6×** the detected counts."
- Speed: "the median conservative speed ratio was **19.13×**, the family-stratified nested-bootstrap 2.5th-percentile bound was **18.33×**, and **21 of 24** scenes had $S_{gate}>1$."
- Exact frozen values (§S6): $S_{gate}$ median = **19.127043091646133**, LB = **18.328492357080282**.
- Frame: "These are fixed-dwell and elapsed-time benefits bought with power, not photon-efficiency gains."
- R19 correction record (main §5 subsec + Table S-correction): the 3 discrepancies (axis $\log(\nu\rho)$→$\log T_{opt}$; family-label→six-fixed-strata; nonstandard isotonic→PAVA+censoring) and "agreed on all 18 audited numerical outputs." Source-of-record `M1_VERDICTS_SPEC_CORRECTED_R19.json` + `ANALYSIS_CORRECTION_DISCLOSURE.md`.

### 9b. Paper-1 frozen novelty wording (round-4/5/8, PLACEHOLDER_LEDGER)
- NO "first" anywhere. Ridge claim scoped to: ideal nonparalyzable + active-start + scalar integrated count + exact finite-window FI + principal ridge + (ρ,ν) asymptotics.
- Intro Para 4 "safe novelty sentence" (verbatim, already inlined): "…spatial occupancy is varied under equal mean detector load and equal incident dose, revealing a contrast-controlled boundary between multiplex-limited and photon-limited high-flux single-pixel imaging."
- Panel (d) may NOT claim coincidence with ρ*; discussion MUST carry the peak-irradiance + bench-validation limitation sentence.
- R12 relation-to-prior-work paragraph (§4.4) is frozen verbatim; Grönberg full-text sign-off outstanding (S5 hard item).

### 9c. R33 theorem names & forbidden claims
- Theorem name: **Canonical-Confusion Ledger Reciprocity Theorem** (distinct from method name DLGI).
- MOLT theorem name: **Masked Laplace–Power-Sum Theorem** (NOT "elementary-symmetric-polynomial theorem").
- Forbidden (DLGI): no "information conservation", no "free medium information", no "paired trades scene into medium", no additive conservation law. Correct wording: "two projections of the same posterior ambiguity, coupled through the posterior cross-covariance K."
- Forbidden novelty "firsts" list (R33 §5 + R31 §6) — see §4c and §5a.

### 9d. R34 authorized DLGI claim (until campaign passes — verbatim)
> "DLGI extracts a competitive point estimate of medium correlation time and fluctuation strength from the same bucket record used for scene reconstruction, with no added acquisition and no detected scene penalty; model-calibrated t_c uncertainty remains under confirmation."
- "Do not use 'two certified products' yet." On pass, restore: "One bucket stream, two model-certified products…" ("model-certified" + validated scalar-gain domain MANDATORY).

### 9e. Act-III master verdict (SOFTWARE_SATURATION_VERDICT.md — verbatim)
> "No current software candidate clears the operator's novelty × simplicity × generality × image-positive bar. The software layer for these frozen channels is saturated at the requested effect size (≥1 dB deployable image gain)."
> Identity explanation: "Information lost = E[Var(hidden state | record)]; the bucket-count record is too thin to identify the hidden state…and five independent mechanisms measured the bound as binding."
- Bridge closing: "…the certified global operating law (ridge tracking, +1.87 dB / 19.13×, M1) is the practical optimum of this regime, now bracketed from both sides by preregistered evidence."

### 9f. MOLT authorized/forbidden fence (R32, verbatim — if MOLT kept)
- Authorized: "photon-for-rank trade"; sparse support-targeted masks; second sensing operator via microcell occupancy; may add little wall-clock in switching-limited systems.
- Forbidden: universal/dense null-space jailbreak; exact 5% recovery without support/prior certificate; photon efficiency / dose neutrality / zero-cost / zero-calibration; "negligible time" on MHz arithmetic alone; hiding matched-photon failure.

---

## 10. GAP LIST (needed by the flagship, does not yet exist)

| # | Gap | Blocking | Notes |
|---|---|---|---|
| G1 | **Unified notation table** | writing phase | Paper-1 uses ρ,ν,$C_u$,Γ,J; Paper-2 adds $c_v$,$\rho_R$,$\rho_c$,K,A,B; DLGI adds t_c,CV,κ. No single symbol table reconciles the three acts. |
| G2 | **DLGI R34 confirmatory campaign** | R34 GO issued, not run | 27-cell grid, ≥5000 calib + ≥1000 confirm records/cell, ≥12 fresh scenes, frozen Neyman intervals. Act IV lead is CONDITIONAL on this passing. Until then, only the "under confirmation" claim is licensed. |
| G3 | **Flagship hero figure not designed** | architecture ruling (R35) | Existing heroes are per-paper (paper1 mechanism, paper2 jitter cap, DLGI 3-panel). A single figure carrying the 4-act arc (identity → law → saturation boundary → dual-ledger) does not exist. |
| G4 | **Act-III single boundary figure** | writing phase | The five-certificate table (§3a) has no companion figure; a double-sided-bracket schematic (software ceiling vs illumination ceiling) is implied by the arc but unbuilt. |
| G5 | **Merged bibliography** | writing phase | Two refs.bib + ~39 new DOIs (§8b/c) to dedupe; `companion` cross-ref key must be removed on merge; Grönberg full-text sign-off still open. |
| G6 | **Optica template port (M9)** | pre-submission | Both papers written in `article` class; figures at approximate 8.6/17.8 cm; opticajnl.cls port + font/resize pass pending for all figures. |
| G7 | **Author block / repo URL / funding** | USER DECISION | \SPH placeholders in all four .tex files (paper1 main, paper2 main + supp). |
| G8 | **Act-I/II unification prose** | R35 outline | Paper-1 ridge law (deterministic) and Paper-2 jitter cap (random-hold) are currently two separate theorem sets; the flagship needs them as ONE Act-I identity→law progression with the crossover Eq. as the bridge. |
| G9 | **Decision on cut scope** | R35 | Paper-1 Study-1/Study-2 campaign is OE-derivative; must be adjudicated CUT vs materials-bank vs retained-as-contrast-mechanism. §5 anomaly/ladder/MOLT placements also need the operator's fork decision (arc-only vs audit-line). |
| G10 | **R16 2M-frame diagnostic** | open, flagged in paper-2 | Score-identity diagnostic resolved the histogram bias but the full 2M-frame confirmatory run is still noted as open in paper-2 §S4; carry-forward flag. |

---

*End of map. Source files are absolute-pathed under `D:\GI_another\`. This
document is inventory only; it prescribes no prose and makes no architecture
decision — those belong to the operator + R35.*
