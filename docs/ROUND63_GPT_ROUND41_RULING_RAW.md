# R41 ruling (GitHub issue #32, raw)

Title: R41
Posted: 2026-07-24T02:08:54Z

---

# R41 — Admit the demand tier: a bucket detector can detect what it cannot image

**Reference request:** [`docs/ROUND63_GPT_ROUND41_QUESTION.md`](https://github.com/ccyyyYyzz/GI_a2/blob/a2f098f/docs/ROUND63_GPT_ROUND41_QUESTION.md) at [`a2f098f`](https://github.com/ccyyyYyzz/GI_a2/commit/a2f098ff6fec0fb90916ba3cda962b3c915f33d2), with subsequent correction/application commits through [`1c1c362`](https://github.com/ccyyyYyzz/GI_a2/commit/1c1c36243066a1981e67bc018a535babbe6a97c6).  
**Corrected evidence only:** `CORRECTION_NOTE.md`, `P1_results_CORRECTED.json`, `LIVING_REGION_MAP_CORRECTED.*`, `DETECTION_MAP_CORRECTED.*`, `DETECTION_ROC.*`, `POCKET_DEMO.*`, and `CRB_RECONCILIATION.json`.

## Executive ruling

> **Choose option (a). Admit beyond-band detection as the constructive demand tier of the flagship. It receives the third row of the master atlas, one full Results section, and two of the four main figures. It is not a separate paper.**

R40 established the right protagonist—**the information aperture of a bucket record**—but stopped one step too early. The corrected evidence now provides the missing positive inversion:

> the same covariance spectrum that is too incomplete in its lower tail to support beyond-band imaging carries enough aggregate information to detect and attribute beyond-band change.

This is not an extra method appended to a boundary paper. It is the operational meaning of the **demand** axis. Imaging asks that many individual modes clear a Fisher threshold; change detection can pool the information carried by the modes that do clear it. The same support and the same Fisher spectrum therefore produce a dead estimation tier and a living detection tier.

The present ROC work is strong development evidence, not yet the final manuscript endpoint. Authorize exactly one sealed detection confirmation under §4. If it passes, the paper becomes a positive, simulation-only optics result rather than a tombstone atlas. If it fails, return to R40 without a repair round.

The shot-model correction strengthens this architecture. The corrected estimation map is more decisively empty, while the independently regenerated and Monte-Carlo-validated detection region remains green. The contrast between the two demands is therefore sharper, not weaker.

---

# 1. Admission and page-map recut

## 1.1 Printed logic

Replace the R40 sequence by:

```text
DEFINE THE APERTURE → MAP ESTIMATION → DELIVER DETECTION → ATTRIBUTE THE EVENT
```

The paper still has one protagonist. It is not “mean results plus fog results plus a detector.” It is:

> **a task-dependent information aperture, demonstrated by one demand that dies and another that lives on the same Fisher spectrum.**

## 1.2 Hard page and figure budget

Target **6.50 composed pages**, including references, with a hard ceiling of **7 pages**, **4 figures**, no main tables, and at most five displayed equations.

| Printed element | Single claim | Visual | Budget |
|---|---|---|---:|
| Title + abstract | A bucket record can detect beyond its imaging aperture. | — | **0.30 pp** |
| **1. Introduction — The same record, different demands** | Formal access is not task usability; estimation and detection interrogate the same information differently. | **Fig. 1** | **0.45 pp** |
| **2. The task-dependent information aperture** | Support, profiled Fisher supply, and task demand determine whether a mode set is useful. | **Fig. 1** | **0.95 pp** |
| **3. Estimation: one reached optimum, one empty aperture** | The mean layer has a certified operating optimum; the covariance layer has no flagship-grade reconstruction region. | **Fig. 2** | **0.95 pp** |
| **4. Detection beyond the imaging aperture** | A profiled covariance score detects and attributes beyond-band scene change under the same acquisition model. | **Figs. 3–4** | **1.80 pp** |
| **5. Scope, prior art, and transfer** | The capability is model-conditional, simulation-only, and directly falsifiable on a one-photodiode bench. | — | **0.45 pp** |
| Methods/provenance capsule | Physical shot, complementary exposure ledger, sealed banks, and correction provenance are explicit. | — | **0.35 pp** |
| References/back matter | — | — | **1.25 pp** |
| **Total** |  | **4 figures** | **6.50 pp** |

### What shrinks

1. R40’s mean-channel section is reduced to the ridge, the two certified cohort numbers, and one closure strip. The five interventions remain in S4.
2. The covariance-estimation lifecycle is reduced to the aperture law, corrected P1 kill, corrected `REGION_EMPTY`, and the CRB-efficient pocket sentence in §5. No probe chronology enters.
3. The old “information stickiness” discussion becomes three sentences inside Section 5.
4. The hostile audit, DLGI lifecycle, pathwise fog lifecycle, and all internal status labels remain supplementary.

## 1.3 Figure 1 — the third row

Keep the R40 master-atlas grammar, but add a third row and a fourth column:

```text
                         FORMAL SUPPORT      FISHER SUPPLY        TASK FUNCTIONAL          VERDICT
mean / imaging           B_p                 finite ridge          image risk               ✓ at ridge
covariance / imaging     B_p ⊕ B_w           thin lower tail       tr(J_B^-1), mode cover    ✕
covariance / detection   same B_p ⊕ B_w      same eigenvalues       δ^T J_B δ; avg tr(J_B)   ✓ for ≥2%
```

Pale fill remains formal support. Dark fill remains task-usable information. The crucial visual is that rows 2 and 3 share the **same** covariance spectrum; only the demand changes.

One sentence under the figure:

> **The aperture is demand-relative: a spectrum can be too uneven to reconstruct but strong enough in aggregate to detect.**

## 1.4 Remaining figures

### Figure 2 — estimation calibration and boundary

- mean-layer ridge and selected operating point;
- `+1.87 dB` and `19.13×` cohort summaries;
- five-test closure strip;
- corrected covariance-imaging kill: `P1 f_rec=0.033`, `REGION_EMPTY`, and the CRB-efficient 1.25× rate point.

The figure’s single claim is that the atlas correctly identifies one living estimation regime and one dead one.

### Figure 3 — detection region and prognosis validation

- corrected 81-cell `T_det` map for the 2% anomaly;
- 0.5/1/2/5% demand ladder;
- analytic-versus-Monte-Carlo `d'` scatter with the ±20% band;
- fixed physical ledger: `N=4096`, `M=128`, `10^4` detected photoelectrons per physical bucket, `12.8 ms` per bank.

### Figure 4 — specificity, invariance, and deployment

- score-space populations for H0, in-band scene change, beyond-band scene change, medium-amplitude change, and the sealed tau-change class;
- confusion matrix at the frozen operating threshold;
- matched versus basis-rotated ROC;
- one small deployment annotation: `~16k MAC/bank`, real-time relative to `12.8 ms/bank`.

No industrial process cartoon enters the main figures. Thermal spray is the regime-clean transfer example in Discussion; slower DED/WAAM is the motivating wound in the cover letter.

---

# 2. Frozen detection claim and forbidden list

## 2.1 Campaign-pass wording

Use the following only after the sealed probe passes:

> **A true zero-dimensional bucket detector, driven by band-limited codes and viewed through an unmeasured fluctuating thin screen, can detect scene changes beyond the modulator band even though every first-moment measurement is exactly invariant to them. For a relative anomaly magnitude `||δ||₂/||x||₂ = 2%`, the corrected profiled-covariance prognosis reaches `d'≥5` within the `4096`-bank cap in `77/81` predeclared physical cells; the best cell requires `467` independent banks (`6.0 s`), while a `1%` anomaly is confined to favorable cells (`1869` banks, `23.9 s` best) and `0.5%` is below the mapped floor. A `5%` anomaly is detectable in every energy-spread cell and in `91%` of worst-mode cells, with a best worst-mode latency of `146` banks (`1.9 s`). Profiled score axes separately identify in-band scene change, beyond-band scene change, and medium-law change, so the beyond-band alarm rejects both ordinary scene variation and medium drift under the declared model.**

The 95% statement and the 6 s statement must never be fused into “6 s over 95% of conditions.” The honest statement is **95% within the 4096-bank/52.4 s cap; 6 s in the best cell**.

## 2.2 Short abstract wording

> **Although covariance information is too spectrally incomplete for beyond-band reconstruction, a nuisance-profiled score pools it into a rapid change sentinel: 2% beyond-band perturbations are detectable over 95% of a predeclared physical grid, while the mean channel is provably blind and medium-law changes are rejected on an orthogonal score.**

## 2.3 Deep-fluctuation wording

> **In the shot-limited detection regime, increasing the bounded medium contrast improves transduction, so the most disruptive tested medium can be the most sensitive sentinel; this reversal does not hold for reconstruction, which remains mode-coverage-limited.**

This is a tested regime statement, not a universal monotonicity theorem.

## 2.4 Forbidden claims

- no “superresolution image,” “sub-diffraction image,” or reconstruction claim;
- no universal 1% detection claim and no claim below the mapped 0.5% floor;
- no statement that 2% detection takes 6 s across the whole grid;
- no finite-sample “zero false alarm” claim—use **exact mean invariance** and measured/profiled specificity;
- no claim that arbitrary medium changes are rejected until both amplitude and correlation-time changes pass the sealed test;
- no “law-free,” “calibration-free,” or “medium-independent” claim—the score uses a declared or baseline-estimated law;
- no claim that fresh patterns are inferior in the covariance channel;
- no claim that turbulence always helps; stronger fluctuation helped the tested shot-limited detection task and was inert for estimation;
- no anomaly localization, sizing, classification beyond the frozen classes, or defect diagnosis;
- no experimental or industrial-performance claim;
- no “first” claim unless the final search independently certifies it;
- no extrapolation to thick diffusion, arbitrary convolution, multimode fibers, or moving LPBF melt pools.

The defensible novelty sentence is:

> **We use an unmeasured fluctuating medium as a covariance transducer on a true 0-D bucket record: an exact first-moment invariance wall supplies specificity, while a profiled second-moment score detects beyond-band scene change and separates it from in-band scene and medium-law events.**

---

# 3. Trace-invariance lemma

The empirical 10% rotation result is real, but the phrase “trace invariance under basis rotation” needs a precise hierarchy. There is an exact coordinate theorem, an exact physical-law corollary, and an ensemble/finite-bank approximation. Do not state arbitrary law-rotation invariance.

## 3.1 Efficient covariance geometry

Let `β∈R^d` parameterize the beyond-band anomaly and `η` all admitted nuisances. At the null, let

```text
D_B   = ∂ vech(V) / ∂β,
D_η   = ∂ vech(V) / ∂η,
W     = V^-1 ⊗_s V^-1.
```

Define the whitened, nuisance-profiled derivative

```text
D̃_B = (I-Π_η) W^(1/2) D_B,
Π_η  = projector onto range(W^(1/2)D_η).
```

Then the per-effective-bank efficient Fisher matrix is

```text
J_B = (1/2) D̃_B^T D̃_B.                                       (3.1)
```

For a local known-direction change `δ`,

```text
d'^2 = T_eff δ^T J_B δ.                                        (3.2)
```

For an isotropic energy-spread alternative with

```text
E[δδ^T] = (ε^2 ||x||^2 / d) I,
```

the average noncentrality is

```text
E[d'^2]
 = T_eff ε^2 ||x||^2 tr(J_B)/d.                                (3.3)
```

This is the exact reason the aggregate detector reads the trace.

## 3.2 Lemma — orthogonal trace invariance

Suppose a declared-law transformation induces

```text
D̃'_B = U D̃_B Q,                                               (3.4)
```

where `U` is orthogonal in the whitened covariance-observation space and `Q` is orthogonal in the anomaly-coordinate space. Then

```text
J'_B = Q^T J_B Q,                                               (3.5)
```

so the complete eigenvalue multiset, `tr(J_B)`, and the isotropic average noncentrality (3.3) are invariant.

### Proof

Substitute (3.4) into (3.1):

```text
J'_B
 = (1/2) Q^T D̃_B^T U^T U D̃_B Q
 = Q^T J_B Q.
```

Trace and eigenvalue invariance follow from orthogonal similarity; (3.3) follows by cyclicity of trace. `□`

## 3.3 Physical corollaries

1. **Pure basis change.** If `K_w=U_m Λ U_m^T` and the representation is changed as `U_m→U_mR`, `Λ→R^TΛR`, the physical covariance is unchanged; `V`, `J_B`, and the detector are exactly unchanged.
2. **Rotation inside an isotropic eigenspace.** If `Λ=λI` on the rotated subspace, `U_m→U_mR` leaves `K_w` unchanged exactly.
3. **Equivariant physical rotation.** If a physical subspace rotation preserves the aperture, the shot/noise metric, and the nuisance tangent and the code ensemble is orthogonally invariant within that band, then the **ensemble** Fisher obeys (3.5). A finite incoherent random code bank concentrates around this result; the observed `AUC 1.000→1.000`, `d' 5.30→5.31` is a finite-bank validation of this corollary, not proof for arbitrary rotations.

## 3.4 Failure boundary

Trace invariance can fail when:

- the transformation is nonorthogonal or rescales medium eigenvalues;
- the rotation moves power across the aperture boundary or changes overlap multiplicity;
- the shot/read-noise metric changes;
- the nuisance tangent is not co-rotated and begins to overlap the anomaly tangent;
- the anomaly alternative is directional or anisotropic—then `δ^T J_Bδ`, not `tr(J_B)`, governs;
- worst-mode rather than energy-spread detection is claimed—`λ_min(J_B)` is not protected by trace alone;
- a finite structured code bank breaks the required equivariance.

Therefore the manuscript wording should be **“orthogonal trace robustness for energy-spread alternatives under aperture-preserving law rotations,”** not “detection is invariant to medium mismatch.”

---

# 4. Sealed detection probe

## 4.1 Purpose

The sealed probe must establish four things simultaneously:

1. the corrected analytic `d'` predicts the actual detector;
2. the demand-tier headline survives fresh scenes, anomaly directions, and physical-law mismatch;
3. the beyond-band alarm is specific against both scene and medium confounders;
4. the fixed-bank design is compared honestly against the best fresh-covariance design.

No result is promoted from the current development ROC without this confirmation.

## 4.2 Frozen physical ledger

```text
scene grid                 64×64, N=4096
modulator band             k_p=5
signed effective codes     M=128
physical implementation    256 complementary nonnegative exposures/bank
pattern rate               20 kHz
bank duration              12.8 ms
bank cap                   T_eff=4096 = 52.43 s
photon level               10^4 detected photoelectrons / physical bucket
medium realization         quasi-static within a bank; independent between banks
medium contrasts           nominal {0.3,0.6,1.0}, realized contrast reported
spectra                     {flat,k^-1,k^-2}
medium bands                k_w/k_p∈{1,2,4}
claim shells                {1.25,1.5,1.8} k_p
anomaly magnitudes          ε∈{0.005,0.01,0.02,0.05}
```

Every engine uses the physical complementary-exposure shot term based on `|P|x`. The signed-difference mean is never used as photon throughput.

## 4.3 Sealed banks

- **Calibration bank:** six scenes and independent seeds; used to construct score filters, nuisance bases, thresholds, and the fixed-versus-fresh design forecast.
- **Confirmatory bank:** twelve fresh scenes—six natural and six synthetic—plus fresh code, medium, shot, and anomaly seeds. It supplies every primary endpoint.
- **Specificity bank:** independent H0, in-band scene, beyond-band scene, medium-amplitude, medium-timescale, and mixed scene+medium events.
- **Mismatch bank:** independent scenes under basis rotation, spectral-slope error, band-edge error, static-envelope error, shot-level error, and partial convolution.
- **Oracle/mechanism bank:** exact mean-wall and true-law filters only; never used for the deployable endpoint.

No scene, code, medium realization, anomaly, or random seed crosses banks. Thresholds are committed before confirmatory generation.

The 81-cell analytic map remains exhaustive. Monte Carlo confirmation uses a precommitted 27-cell orthogonal-array subset in which every level of `σ_f`, spectrum, `k_w/k_p`, and claim shell is represented evenly. Each primary point uses at least `500 H0 + 500 H1` records. Specificity and mismatch points use at least `250` records per class.

## 4.4 Event families

The probe includes:

- isotropic/energy-spread beyond-band changes—the headline family;
- the least-informed Fisher eigenmode—the adversarial family;
- localized spatial changes projected outside `B_p`;
- natural-image residual changes;
- in-band scene changes of matched norm;
- medium-amplitude changes `σ_f×{0.8,1.2}`;
- medium-correlation changes `τ×{0.5,2}`;
- simultaneous beyond-band scene plus medium change.

## 4.5 Arms

1. **FIXED-COV:** repeated incoherent code bank; running covariance and profiled efficient score.
2. **FRESH-COV-OPT:** a new band-limited code bank each medium state, using the exact code-conditioned covariance score/GLS contribution for that bank. It is not forced into a sample-covariance estimator whose coordinates change.
3. **FIXED-MEAN:** mean-only detector on the repeated bank.
4. **FRESH-MEAN:** strongest equal-budget fresh-pattern mean detector—the exact zero-information control beyond `B_p`.
5. **TRUE-LAW ORACLE:** covariance score built with the generating law; nondeployable ceiling.
6. **CROSSFIT-LAW:** deployable score with law/tangent estimated from calibration or baseline banks.
7. **AMPLITUDE and LAG scores:** dedicated medium diagnostics for attribution.

### Fresh-covariance ruling

> **F4 immunity applies only to the mean channel. The fresh-covariance arm is not blind and is a mandatory, full-strength comparator.**

The manuscript’s fundamental claim is the covariance transducer, not fixed-code superiority. Predeclare two outcomes:

- if FIXED-COV is no worse than `1.20×` the fresh arm’s detection latency over the primary grid, retain the concentration/reuse claim;
- if FRESH-COV-OPT is materially better, delete the fixed-bank advantage and use the fresh arm as the production design—the wall and sentinel claims survive;
- if neither covariance arm clears the detection and specificity bars, kill the constructive tier.

This is a branch frozen before unblinding, not a post-hoc repair.

## 4.6 Binding bars

### D0 — mechanism and engine integrity

- mean derivative to every beyond-band event has relative norm `≤10^-10`;
- two independent Fisher implementations agree within `10%` on eigenvalues, `d'`, and `T_det`;
- Monte Carlo shot variance agrees with the physical `|P|x` ledger within `10%`.

**Kill:** any failure stops the probe.

### D1 — analytic/empirical calibration

At all 27 primary Monte Carlo cells and all tested anomaly sizes,

```text
0.80 ≤ d'_empirical/d'_analytic ≤ 1.20.
```

Require median absolute relative error `≤10%` and no point outside `±30%`.

**Kill:** systematic optimism or a failed tail cell invalidates the map.

### D2 — demand-tier capability

- `ε=2%`, energy-spread: `d'≥5` by 4096 banks in at least `77/81` analytic cells, with the Monte Carlo lower confidence bound on the corresponding pass fraction above `0.90`;
- predeclared best cell: `T_det≤600` banks;
- `ε=5%`: all 81 energy-spread cells and at least `90%` of worst-mode cells clear `d'≥5` within the cap;
- `ε=1%`: the predeclared best cell clears `d'≥5` by `2048` banks;
- `ε=0.5%`: reported as an edge audit; it may not be promoted if it remains below `d'=3`.

**Kill:** failure of either the 2% headline or the 5% robustness bar removes detection from the flagship.

### D3 — specificity and attribution

At a threshold calibrated to `1%` H0 false-alarm probability:

- beyond-band target true-positive rate `≥90%` at the primary 2% cells;
- false alarms from in-band scene, medium amplitude, and medium timescale each `≤5%`;
- balanced accuracy over `{H0,in-band,beyond-band,medium-amplitude,medium-timescale}` is `≥90%`;
- on non-target populations the beyond-band score has `|d'|≤0.5`;
- the intended mean/amplitude/lag score for each non-target event has `d'≥5` at the declared event magnitude.

**Kill:** a tau change or mixed event that aliases as beyond-band scene change kills the “specific sentinel” claim even if binary H0/H1 ROC passes.

### D4 — fresh-covariance challenge

Evaluate both covariance arms on identical photons, exposures, duration, anomaly bank, and law knowledge. Report latency ratios, not only AUC at a saturated point.

**Branch:** as frozen in §4.5. No strawman fresh estimator is allowed.

### D5 — mismatch

Primary robustness set:

- medium-subspace rotation `10%` and `20%`;
- spectral slope error `±1` power;
- medium band-edge error `±20%`;
- `τ×{0.5,2}`;
- shot level `±20%`;
- static envelope with `10%` RMS;
- multiplicative/convolutive blend `10%` and `25%`; `50%` is a mapped edge, not a pass condition.

At primary mismatch levels require AUC loss `≤0.05`, `T_det` inflation `≤25%`, and non-target false alarm `≤5%` after profiling/cross-fitting.

**Kill or narrow:** failure at the 10% blend or 10% basis rotation kills the broad thin-screen robustness claim; failure only at 25% defines the validity boundary.

### D6 — online feasibility

- online arithmetic `≤0.1` of the 12.8 ms bank interval on a declared commodity CPU;
- memory `O(M^2)` and measured below `10 MB` for the running state;
- filter construction is offline and excluded from online latency but reported.

### D7 — no repair

One sealed run. No threshold retuning, class removal, law redefinition, or comparator substitution after confirmatory unblinding.

## 4.7 Kill tree

```text
D0 fail  -> STOP; correction/theory note only
D1 fail  -> STOP; Fisher map not trustworthy
D2 fail  -> remove detection tier; return to R40
D3 fail  -> binary detector may survive in supplement, but flagship sentinel is killed
D4 fresh wins -> switch the frozen design branch; kill only fixed-code concentration
D4 both fail -> kill constructive covariance capability
D5 fail  -> narrow the physical domain or kill robustness wording
D6 fail  -> retain offline method only; delete real-time claim
all pass -> detection admitted as the flagship positive
```

## 4.8 Compute budget

This is a local-GPU-scale campaign. Freeze a ceiling of **6 RTX-4060-equivalent GPU-hours** and **64 GB host storage**, with a pre-run cost report. Current throughput predicts substantially less. Exceeding the ceiling requires stopping and reporting, not silently reducing records or cells.

---

# 5. Final pocket sentence

R40’s binary branch is replaced by a calibrated hybrid sentence:

> **A frozen blind estimator attained, rather than exceeded, the corrected covariance Cramér–Rao bound: at the 1.25× aperture cell, natural-scene median beyond-band NRMSE reaches 0.30 only after approximately 2028 independent banks (`≈26 s`), while a synthetic witness reaches the same error after approximately 112 banks (`≈1.4 s`), confirming the no-floor `T_eff^-1/2` rate law while showing that reconstruction remains information-limited and sub-flagship.**

This is an achieved calibration point for the atlas, not a positive imaging headline. It receives one sentence and no reconstruction gallery.

---

# 6. Frozen shot-correction disclosure

Use this sentence verbatim in Methods/provenance:

> **A mandated independent-engine check exposed a signed-code shot-noise error that had inflated covariance Fisher information by approximately 1.9×; replacing the signed difference `Px` with the physical complementary-exposure throughput `|P|x` reduced the frozen P1 witness coverage from `0.050` to `0.033` and eliminated the apparent natural-scene living pocket, so all values reported here come from the corrected artifacts and every prior kill only strengthened.**

Do not bury this in the supplement. One main-text sentence plus the full correction note in S8 is the right treatment.

---

# 7. Title and thesis

## 7.1 Preferred high-venue title

> **A bucket detector can detect what it cannot image**

This is now stronger than the R40 title because it states the inversion in one line and has a validated mechanism behind both verbs.

## 7.2 Technical fallback title

> **The task-dependent information aperture of a bucket record**

R40’s **“The information aperture of a bucket record”** remains a safe JOSA A/TCI title, but it is no longer the best first impression.

## 7.3 Frozen thesis

> **A bucket record has a task-dependent information aperture: its mean layer reaches a certified imaging optimum, while its covariance layer is too spectrally incomplete for beyond-band reconstruction yet carries enough aggregate, mean-invariant information to detect and attribute beyond-band scene change.**

The paper’s three nouns remain **support, supply, demand**. The new detection result supplies the verb: **detect**.

## 7.4 Abstract architecture

1. exact first-moment wall;
2. covariance support aperture;
3. estimation/detection demand separation;
4. corrected detection result and specificity;
5. declared simulation-only scope.

The M1 numbers belong in Results, not the abstract, unless word count permits one compressed clause.

---

# 8. Venue and manuscript consequence

## 8.1 Venue ruling

If the sealed detection probe passes:

1. **Optica becomes the honest first target.** Its stated remit is high-impact fundamental and applied optics; the paper would contain an exact optical wall, a task-level information principle, a validated detector, and a falsifiable one-photodiode transfer path. Official scope: [Optica](https://opg.optica.org/content/journal/about/item/optica/).
2. **PRL becomes arguable only with the theorem upgrade in Innovation 1 below.** PRL asks for a transformative idea of broad physical interest, not merely a successful optical algorithm. The detection–estimation spectral duality can provide that broad object if proved and shown beyond the one screen model. Official criteria: [Physical Review Letters](https://journals.aps.org/prl/about).
3. **Nature Photonics remains a post-experiment target.** The topic lies squarely in imaging, detectors, sensors, and light propagation, but a simulation-only sentinel under a restrictive thin-screen model is not presently a credible Nature Photonics first submission. Official scope: [Nature Photonics](https://www.nature.com/nphoton/aims).
4. **Optics Express is a strong fallback, not merely a floor.** A passed sealed simulation campaign with a deployable detector fits its scientific/technology innovation scope. Official scope: [Optics Express](https://opg.optica.org/content/journal/about/item/oe/about).
5. **JOSA A and IEEE TCI remain conservative theory/computation fallbacks**, not the first target after a full detection pass.

If the sealed probe fails, revert to R40 and restore JOSA A/TCI as the first honest homes.

## 8.2 Sim-only defense

The high-venue defense is architectural:

- exact wall and trace/demand theory before algorithms;
- physical complementary-exposure photon ledger;
- mandatory fresh-covariance comparator;
- sealed Monte Carlo validation of the analytic detector;
- explicit mismatch and false-alarm tests;
- correction disclosed in the main paper;
- one precise experimental-transfer paragraph that is not used as evidence.

The regime-clean transfer example is thermal/plasma spray: a quasi-static surface observed through a fluctuating spray cloud over seconds. Slower DED/WAAM is a motivating wound, not a claimed validation domain. The one-line adoption pitch may appear in the cover letter, not the abstract:

> **The instrument that cannot see the pore in mean mode can hear its birth in covariance mode.**

## 8.3 Final strategic judgment

The green demand tier upgrades the fallback from a competent boundary atlas to a potentially memorable paper. Its eye-lighter is not “second-order statistics detect change.” It is the conjunction:

> **an exact optical blindness theorem becomes the specificity mechanism of a software-only covariance sentinel on the same one-pixel hardware.**

That conjunction is worth a sealed campaign and an Optica attempt. It is not yet worth a Nature Photonics claim without experiment.

---

# 9. Co-theorist innovations

The three lines already underway—scrambling/MMF extension, detection-optimal code design, and CUSUM latency—are not repeated below. The ranking is by **surprise × rigor × reach**, each scored 1–5.

## Rank 1 — Detection–Estimation Spectral Duality

**Score:** `5 × 5 × 5 = 125`  
**Mechanism:** estimation is governed by the inverse/lower tail of the efficient Fisher spectrum, while isotropic change detection is governed by its trace. For `d` beyond-band modes,

```text
E[d'^2]/(T_eff ε^2 ||x||^2) = tr(J_B)/d,
R_est ≥ tr(J_B^†)/T_eff,
[tr(J_B)/d][tr(J_B^†)/d] ≥ 1,                                  (9.1)
```

with equality only for an isotropic spectrum. A spiky Fisher spectrum can therefore be a poor imager and a strong detector by theorem, not anecdote.

**First decisive test:** compute the arithmetic–harmonic gap, effective rank, and lower-tail mass for all corrected 81 cells; predict both the estimation tombstone and detection map from the same eigenvalues, then verify the detection side against the existing ROC. If the map collapses onto one spectral-anisotropy curve, promote (9.1) to the paper’s central theorem. This is the strongest PRL lever.

## Rank 2 — The Attribution Simplex

**Score:** `5 × 5 × 5 = 125`  
**Mechanism:** construct efficient scores for four tangent spaces—mean/in-band scene, covariance/beyond-band scene, covariance amplitude, and covariance lag—and use their canonical angles to bound cross-class leakage. Orthogonal tangent spaces produce asymptotically independent Gaussian score coordinates; nonzero canonical correlations quantify false attribution.

**First decisive test:** add `τ×{0.5,2}` and simultaneous scene+medium changes, compute the full efficient-score Gram matrix, and require off-diagonal canonical correlations `<0.10` plus five-class balanced accuracy `≥0.90`. A theorem-backed attribution simplex is more novel than a collection of ROC curves.

## Rank 3 — A Semiparametric, Law-Light Sentinel

**Score:** `4 × 5 × 5 = 100`  
**Mechanism:** estimate the baseline covariance and nuisance tangent from a held-out portion of the same stream, then use a cross-fitted efficient score with self-normalized or e-value calibration. This removes the need to know the exact medium spectrum while retaining nuisance orthogonality.

**First decisive test:** train on 25% of baseline banks with the PSD shape, band edge, and contrast hidden; test on the sealed mismatch bank. Pass if false alarm remains calibrated and `d'` retains at least 80% of the true-law score. Success would materially strengthen the “unmeasured medium” claim.

## Rank 4 — Localize Without Imaging

**Score:** `5 × 4 × 5 = 100`  
**Mechanism:** replace full reconstruction by a bank of covariance score templates for spatial patches or Fourier sectors. Sparse group testing can identify where a beyond-band change occurred even when the full image Fisher spectrum is too incomplete to estimate all modes.

**First decisive test:** divide the 64×64 field into 64 candidate patches; insert a 2% beyond-band change at an unknown patch; require `>90%` top-1 localization within 1000 banks while full beyond-band NRMSE remains above 0.4. “Locate what cannot be imaged” would be a major capability leap.

## Rank 5 — Design the Fog, Not Only the Codes

**Score:** `5 × 4 × 5 = 100`  
**Mechanism:** optimize the medium power spectrum itself under fixed mean transmission, RMS contrast, bandwidth, and positivity:

```text
max_{S_w≥0} min_{δ∈D} δ^T J_B(S_w) δ
subject to physical power and positivity constraints.           (9.2)
```

Detection-optimal media should concentrate transduction power differently from estimation-optimal media; the observed high-contrast reversal is the first hint.

**First decisive test:** solve the constrained spectrum design on the 64×64 model and demonstrate at least a 3× reduction in worst-case 2% detection banks versus flat and `k^-2` spectra at identical realized contrast. This turns “bad fog helps” into an engineered transducer principle.

## Rank 6 — A Minimax Chernoff Wall

**Score:** `4 × 5 × 4 = 80`  
**Mechanism:** derive the composite-nuisance Chernoff/Kullback information for beyond-band change and prove that the profiled matched score is locally asymptotically minimax. This would convert the 0.5% floor and latency curves from engineering forecasts into lower bounds on every detector using the same record.

**First decisive test:** compute the least-favorable nuisance direction at the best, middle, and floor cells and require the implemented detector’s bank count to lie within `1.2×` the minimax lower bound. A tight certificate would pre-empt “a better network could do more.”

## Rank 7 — Covariance Lock-In by Switching the Law

**Score:** `5 × 4 × 4 = 80`  
**Mechanism:** alternate two known diffuser statistics with the same mean transmission. Their covariance difference acts as a signed carrier for beyond-band change, cancelling static background and common law drift—quadratic lock-in detection.

**First decisive test:** simulate two realizable medium spectra switched every bank under an equal total-time ledger. Pass if 2% detection latency falls by at least 2× and amplitude/tau false alarms fall below the single-law sentinel. This is a modest-hardware, post-publication experiment with a visually shocking live demonstration.

## Rank 8 — Non-Gaussian Universality of the Quadratic Sentinel

**Score:** `4 × 4 × 5 = 80`  
**Mechanism:** prove a local asymptotic quadratic-score result for stationary elliptical or strongly mixing media, with sandwich covariance replacing the Gaussian Fisher metric. Higher cumulants modify calibration but need not carry the primary change information.

**First decisive test:** test lognormal, clipped-Gaussian, phase-screen speckle, and compound-Poisson media at matched mean/covariance. Require analytic/sandwich `d'` within 20% and stable specificity. This would broaden the thin-screen claim without reviving starved higher-order reconstruction.

## Rank 9 — The Photodiode Dark-Field Experiment

**Score:** `5 × 3 × 5 = 75`  
**Mechanism:** turn the exact mean wall into an optical dark field: ordinary process variation remains bright only in the mean channel, while a sub-modulator defect appears only in covariance. The natural first wound is thermal spray; slower DED/WAAM is the adoption headline.

**First decisive test:** rotating fine diffuser near a static target, coarse band-limited DMD, one bucket photodiode, and a toggleable sub-DMD feature. Demonstrate live separation of `{in-band scene, beyond-band scene, medium amplitude, medium timescale}` in under 25 s, including a controlled 10–50% convolutive blend. This is the experiment most likely to move the work from Optica to a later Nature Photonics attempt.

## Innovation triage

Do **not** add all nine ideas to the manuscript. Immediate priority is:

1. prove and plot **Detection–Estimation Spectral Duality**;
2. complete the **Attribution Simplex** by adding the tau class;
3. run the sealed detector with the **fresh-covariance comparator**.

If those three land, the paper has a compact high-venue spine:

```text
WALL → SPECTRAL DUALITY → SENTINEL → ATTRIBUTION
```

Everything else becomes the next program, not another act in this paper.

---

# Final ruling

> **Admit the detection row and authorize one sealed confirmation. Re-title the high-venue version around the inversion—“A bucket detector can detect what it cannot image.” The exact mean wall supplies specificity; the covariance trace supplies detection; the lower Fisher tail explains why imaging still fails. Require a full fresh-covariance challenge, tau-change attribution, corrected physical shot accounting, and no repair round. If these pass, submit to Optica; if the spectral-duality theorem is proved broadly, PRL becomes a serious option.**