# R40 ruling (GitHub issue #31, raw)

Title: R40
Posted: 2026-07-24T00:47:00Z

---

# R40 — Final fallback architecture: the information aperture of a bucket record

**Reference request:** [`docs/ROUND63_GPT_ROUND40_QUESTION.md`](https://github.com/ccyyyYyzz/GI_a2/blob/af3ae1f/docs/ROUND63_GPT_ROUND40_QUESTION.md) at commit [`af3ae1f`](https://github.com/ccyyyYyzz/GI_a2/commit/af3ae1f2b8655c18708ca7cb95c5a6570d2ad7b1).  
**Frozen evidence:** [`1f1e4fe`](https://github.com/ccyyyYyzz/GI_a2/commit/1f1e4fe60c49597ec33629fc4521e14cffdfda71), especially `P1_results.json` and `LIVING_REGION_MAP.md`.  
**Pocket status at ruling time:** neither `POCKET_ACHIEVED` nor `ESTIMATOR_GAP_PERSISTS` has landed on `main` or issue #30. Both editorial branches below are therefore binding.

## Executive ruling

> **Choose option (a), but collapse it into one object: the Fisher-weighted information aperture of a bucket record. The mean and covariance channels are two rows of one atlas, not two acts and not a chronology of campaigns.**

The R36 fallback is now too narrow because it omits the durable result of the fog-DMD arc: the exact distinction between **formal support** and **usable information**. A full two-layer program summary would repeat the R35 bloat failure. The correct paper is the middle path:

> one master information criterion, evaluated on the first two moment layers of the same bucket record.

For moment layer `r`, scene mode `j`, and profiled Fisher eigenvalue `lambda_j^(r)`, define the task supply

```text
S_j^(r) = |a_j|^2 lambda_j^(r).
```

A mode is usable at task threshold `gamma` only when

```text
S_j^(r) >= gamma^2.                                             (R40.1)
```

Everything in the main article has one role relative to (R40.1):

- the **support law** says where `lambda_j` may be nonzero;
- the **profiled Fisher atlas** says how large it actually is;
- the **task demand** says whether the mode matters;
- the certified campaigns test whether the predicted living and dead regions are real.

This is the anti-anthology device. The paper is not “our program found many things.” It is:

> **formal access to a mode is not the same as usable information about it.**

The mean layer supplies the positive calibration point: a finite operating ridge that was reached experimentally in simulation and survived hostile closure. The covariance layer supplies the exact counterpoint: a larger Minkowski support aperture whose flagship-grade usable region is empty because most modes never clear the Fisher demand.

---

# 1. Frozen thesis and titles

## 1.1 One-sentence thesis

> **A bucket record has a Fisher-weighted information aperture: its mean layer can be driven to a certified optimum but cannot escape the modulator band, while its covariance layer enlarges formal support only where mode-by-mode information supply exceeds the task demand.**

This is the manuscript contract. Every paragraph must advance one of three nouns: **support, supply, demand**.

Do not claim a “complete information atlas” without qualification. The paper treats the first two moment layers under declared detector and medium models; it does not prove completeness over every possible statistic or optical measurement.

## 1.2 Titles — ranked

1. **The information aperture of a bucket record**  
   Preferred. One protagonist, concept-first, and broad enough to contain both moment layers without sounding like a review.

2. **From formal support to usable information in bucket imaging**  
   Most explicit statement of the discovery.

3. **What a bucket record can and cannot reveal**  
   Strongest editorial title; use only if the abstract immediately defines the optical setting.

4. **A two-layer information atlas of bucket imaging**  
   Safest descriptive title, but slightly more programmatic.

5. **Moment-channel limits of bucket imaging**  
   Most conservative technical fallback.

Do not use “fog as a second DMD,” “superresolution,” “complete supply-demand ledger,” R-round labels, or method acronyms in the title.

---

# 2. Printed architecture: 6.5-page target, three figures

Hard constraints:

- **7 composed pages maximum**, including references and back matter;
- **3 main figures**, not 4 unless typesetting proves one indispensable;
- **no main tables**;
- **5 displayed equations maximum**;
- **80–105 abstract words**;
- no campaign chronology, internal tree, status-code vocabulary, or “Act I/II” language.

## 2.1 Page map

| Printed element | Single claim | Main visual | Hard budget |
|---|---|---|---:|
| Title + abstract | Support, supply, and demand define what a bucket record can reveal. | — | **0.30 pp** |
| **1. Introduction — One record, two information layers** | The same bucket stream has a mean channel and a covariance channel with different information apertures. | **Fig. 1** | **0.55 pp** |
| **2. The Fisher-weighted information aperture** | Equation (R40.1) unifies hidden-state loss, band support, profiling, and task recoverability. | **Fig. 1** | **0.95 pp** |
| **3. Mean channel: a living optimum with closed headroom** | The mean layer has a finite, reached operating optimum and no material residual count-only scene headroom in the declared class. | **Fig. 2** | **1.05 pp** |
| **4. Covariance channel: broad support, thin usable aperture** | Medium covariance expands formal support, but the publication-scale Fisher map leaves no flagship-grade living region and only one modest natural-scene pocket. | **Fig. 3** | **1.45 pp** |
| **5. Information stickiness, scope, and transfer** | Extra moment order or formal aperture cannot defeat a mode-wise supply deficit; the result is a falsifiable simulation-only boundary map. | — | **0.55 pp** |
| Methods and provenance capsule | Models, ledgers, preregistration, simulation-only status, and reproducibility are explicit. | — | **0.35 pp** |
| References/back matter | — | — | **1.25 pp** |
| **Total** |  | **3 figures** | **6.45 pp** |

The remaining 0.55 page is contingency, not an invitation to add another result.

## 2.2 Section execution

### Section 1 — Introduction

Exactly three paragraphs:

1. A bucket detector compresses each patterned exposure to one scalar, but the scalar has more than one statistical layer.
2. The mean layer can be optimized but has exact geometric walls; the covariance layer can cross a support wall, yet formal support may still fail the task mode by mode.
3. State the thesis, the certified positive, the prelaunch Fisher kill, and the atlas contribution. End by saying that the paper maps both the living and dead regions before any bench claim is made.

Do not open with project history, ghost-imaging applications, or a contribution list.

### Section 2 — The Fisher-weighted information aperture

Keep only the minimum common mathematics:

1. the observed-versus-complete information identity in one line;
2. the profiled Fisher object;
3. the mean support wall and covariance support aperture in one paired statement;
4. the mode-level criterion (R40.1).

The hidden-state identity is the bridge between layers, not a separate paper inside the paper. Louis/Orchard–Woodbury lineage is credited in the first theory paragraph. Full proofs and detector-specific specializations move to S1–S2.

### Section 3 — Mean channel: a living optimum with closed headroom

One paragraph establishes the finite ridge and global control. One paragraph reports the certified positive:

- `+1.87 dB`, lower bound `1.41 dB`, `19/24` positive;
- `19.13x` elapsed-time result, lower bound `18.33x`, `21/24`.

One final paragraph states that five preregistered interventions and a hostile resurrection audit found no material count-only headroom. The five interventions are not named in five subsections. DLGI appears only in the frozen disclosure sentence in §6 below.

### Section 4 — Covariance channel: broad support, thin usable aperture

This section contains exactly the six-sentence fog-DMD summary frozen in §5.2 and no lifecycle chronology. Its logic is:

```text
support law -> exact mean wall -> P1 preflight kill -> 81-cell map -> pocket branch.
```

The pathwise blind route, eight probes, six design gates, fresh-pattern in-band duel, cumulant tests, and every reconstruction attempt are supplement-only evidence for this one boundary conclusion.

### Section 5 — Information stickiness, scope, and transfer

Three short paragraphs:

1. **Information stickiness:** new formal support, more photons, more contrast, and higher cumulants do not help when the mode-wise profiled Fisher spectrum stays below demand.
2. **Scope:** mean-layer results apply to the declared integrated-count detector class; covariance-layer results apply to the stationary thin-screen, multiplicative, law-known setting. Neither is universal over all optical measurements.
3. **Transfer:** specify a rotating fine diffuser near the object, a band-limited DMD, and a bucket detector as a falsifiable future bench test. The article does not depend on that experiment.

---

# 3. Figure architecture

## Figure 1 — Hero: one bucket record, two information apertures

A single master atlas, not a pipeline cartoon.

Use a two-row layout:

```text
                         FORMAL SUPPORT          FISHER-USABLE REGION
mean channel             modulator band          finite ridge; closed headroom
covariance channel       pattern ⊕ medium band   thin rim; REGION_EMPTY at P1
```

Horizontal position is spatial-frequency/task mode. Pale fill denotes formal support; dark fill denotes modes satisfying (R40.1). Place:

- one green marker for the certified mean-channel operating point;
- one red cross for the killed `1.8x` covariance claim;
- one amber marker for the `1.25x` natural-scene pocket;
- one common legend: **support is not supply**.

A referee should understand in ten seconds that the paper is a map, not a method anthology.

## Figure 2 — Mean-layer certificate

One scientific sentence:

> **The mean layer has a reachable optimum and no material residual headroom in the declared count-only class.**

Panels:

1. hidden-state/jitter information ridge with the selected operating point;
2. the two confirmed cohort statistics (`+1.87 dB`, `19.13x`);
3. a compact closure strip showing all five preregistered interventions below the materiality threshold, with the harmful arm visually capped.

No method names, solver labels, audit codes, or gallery mosaics.

## Figure 3 — Covariance-layer law and boundary map

This is the **only** main figure for the entire fog-DMD lifecycle.

Panels:

1. spectral law: `B_p` versus `B_p ⊕ B_w`, next to the exact first-moment wall;
2. empirical aperture validation: `18.76x` in/out separation;
3. P1 and the 81-cell Fisher atlas, showing `REGION_EMPTY`, the mode-coverage bottleneck, and the natural-scene pocket;
4. branch inset:
   - `POCKET_ACHIEVED`: one blind reconstruction inset;
   - `ESTIMATOR_GAP_PERSISTS`: no reconstruction gallery; use a hollow amber marker labelled **Fisher-feasible, algorithmically unrealized**.

Do not give the pocket a standalone figure, title phrase, abstract sentence, or conclusion headline.

---

# 4. Pocket ruling

## 4.1 If `POCKET_ACHIEVED`

The blind 1.25x demonstration **strengthens the paper**, but only as a calibration point for the atlas:

> the atlas predicts one small living room and the estimator walks into it.

Include one inset in Fig. 3 and one sentence in Section 4. The claim remains “a modest natural-scene pocket,” not superresolution leadership. Do not compare the 1.25x number rhetorically with SIM’s approximately 2x regime; the main comparison is against the paper’s own exact mean-channel wall and Fisher prediction.

The pocket may be called achieved only under the already frozen pocket verdict. No new post-hoc visual-quality threshold is invented for the manuscript.

## 4.2 If `ESTIMATOR_GAP_PERSISTS`

Do **not** omit the result and do not turn it into a method claim. Mark the pocket with a hollow point in Fig. 3 and state:

> The atlas identifies a Fisher-feasible 1.25x natural-scene pocket, but the frozen blind lifted-GLS reconstruction did not realize it; this is an estimator-reach gap, not a contradiction of the no-floor theorem.

Full reconstruction diagnostics go to S9. This branch is still publishable because the protagonist is the information atlas. Hiding the failed pocket attempt would weaken trust more than the one-sentence disclosure.

---

# 5. Fog-DMD main-text boundary

## 5.1 Single-figure rule

Only Fig. 3 enters the main text. The full lifecycle receives **six sentences total** outside its caption.

## 5.2 Frozen six-sentence summary

Use the first five sentences verbatim; choose the sixth by the pocket verdict.

> The second moment of the same bucket record obeys a distinct support law: a medium with spatial band `B_w` expands covariance sensitivity from the pattern band `B_p` to `B_p ⊕ B_w`, whereas the first moment remains exactly blind outside `B_p`. At `64×64`, the support theorem and an `18.76×` in/out separation confirmed that this formal covariance aperture exists. The preregistered publication-scale campaign was then stopped before reconstruction because its profiled-Fisher gate recovered only `5%` of the worst-case witness modes at SNR at least three, below the frozen `70%` requirement, even though both aggregate-NRMSE forecasts passed. An exhaustive `81`-cell physical sweep found no flagship-grade region: profiling cost was negligible, and the universal bottleneck was mode coverage rather than photons, fluctuation contrast, or nominal medium bandwidth. The sole residue was a `1.25×` natural-scene pocket under a `k^-2` medium spectrum, with predicted NRMSE `0.10` and about `197` independent epochs. **[Choose one:]** (A) A frozen blind lifted-GLS run realized this predeclared pocket, providing a modest positive check of the atlas without changing the headline claim. (B) A frozen blind lifted-GLS run did not realize this Fisher-feasible pocket, leaving an explicit estimator-reach gap rather than a superresolution result.

No additional E-number, F-number, arm name, or development chronology appears in the main article.

## 5.3 Supplement map

- **S1:** joint missing-information identity, regularity conditions, and hidden-state specializations;
- **S2:** dead-time ridge derivation, scaling laws, R16 artifact closure;
- **S3:** M1 campaign, operating maps, R19 correction and audit bundle;
- **S4:** five mean-channel closure certificates and the hostile resurrection audit;
- **S5:** DLGI reciprocity theorem, confirmatory lifecycle, C2 width kill and v2 prognosis;
- **S6:** covariance model, aperture law, first-moment wall, no-floor theorem and cumulant ordering;
- **S7:** pathwise E1–E8 lifecycle, global collusion and marginalization mandate;
- **S8:** R39 risk gates, exact profiled-Fisher engine and P1 prelaunch kill;
- **S9:** 81-cell living-region map, Monte Carlo validation and pocket demonstration;
- **S10:** prior-art fence, reproducibility manifests, physical-unit ledger and experimental-transfer protocol.

Demotion is not deletion. The supplement is the full scientific record; the main text is the argument.

---

# 6. Frozen disclosure sentences

These sentences are main-text-safe and should be used without euphemism.

## 6.1 R39 P1 prelaunch kill

> **The beyond-band reconstruction campaign was never launched: its preregistered profiled-Fisher preflight recovered only 5% of the witness modes at SNR at least three against a 70% requirement, so the frozen kill tree stopped the study before confirmatory reconstruction despite passing both aggregate-NRMSE forecasts.**

## 6.2 Fresh-pattern in-band dominance

> **An equal-budget fresh-pattern control dominated all nine in-band cells—about eightfold lower scene error while also estimating the medium law—so we withdraw any claim that uncontrolled medium fluctuations are advantageous when new modulator directions are physically available.**

## 6.3 DLGI C2 kill

> **DLGI retained its reciprocity theorem and competitive point estimates, but its certified second-product claim was killed at the interval-width endpoint; the v2 preflight then found 0 of 27 viable cells in the calibration-hull domain, so no repair campaign was launched.**

Each appears once. Full context and numbers live in the supplement.

---

# 7. Anti-anthology rule and venue risk

## 7.1 Binding structural rule

The article may never enumerate the program inventory. Every result must be attached to one cell of Fig. 1 and one of the three questions:

1. **Where can information exist?**
2. **How much usable information is supplied?**
3. **Which task demand does that supply meet?**

The paper does not narrate discovery order. It does not say “we next tried.” It does not use campaign names as section headings. The master atlas is shown before any individual result, so every later number is evidence for a predeclared coordinate rather than another story.

The top structural risk is no longer “two papers stitched together.” It is:

> **“This is an internal postmortem recast as a journal article.”**

The answer is architectural:

- one theorem object, (R40.1);
- one master figure;
- one living mean-layer region;
- one dead flagship covariance region and one honest pocket;
- all chronology and forensic detail outside the main text.

## 7.2 Sim-only desk-reject risk by venue

The following are editorial judgments, not acceptance guarantees, based on the journals’ stated scopes.

| Venue | Desk-risk for this fallback | Ruling |
|---|---|---|
| **JOSA A** | **Low–moderate.** Its scope explicitly includes coherence/statistical optics, computational sensing and imaging, image science, and numerical methods. A theorem-plus-boundary atlas is native to the journal. | **Recommended first choice if the pocket estimator fails.** |
| **IEEE TCI** | **Low–moderate.** Its scope covers fundamental computational-imaging theory, model-based inversion, nontraditional sensing, dynamic information acquisition, and reproducibility. | Strong second choice; make the mode-wise Fisher atlas general enough to outlive the GI example. |
| **Optics Express** | **Moderate–high.** Its scope is broad but emphasizes scientific and technology innovation. A simulation-only boundary map with no realized covariance positive can be read as insufficiently concrete. | Try first only if `POCKET_ACHIEVED` and the optical law/wall figure is exceptionally clean; otherwise do not lead here. |
| **Physical Review Applied** | **High.** The journal seeks significant applications-based physical insight. `REGION_EMPTY`, no experiment, and a modest pocket make the applied advance vulnerable at editorial screening. | Do not use as the first fallback venue. |

Official scope references: [JOSA A](https://opg.optica.org/content/journal/about/item/josaa/), [IEEE TCI](https://signalprocessingsociety.org/publications-resources/ieee-transactions-computational-imaging/about-tci), [Optics Express](https://opg.optica.org/content/journal/about/item/oe/about), and [Physical Review Applied](https://journals.aps.org/prapplied/about).

### Submission order

- If `POCKET_ACHIEVED`: **Optics Express → JOSA A → IEEE TCI**. Skip PR Applied unless an editor signals interest.
- If `ESTIMATOR_GAP_PERSISTS`: **JOSA A → IEEE TCI**. Do not spend the first submission on OE or PR Applied.

The one-paper rule remains binding; these are transfer destinations, not split manuscripts.

---

# 8. Simulation-only credibility architecture

The paper is credible without a bench only if the architecture makes the nature of each claim unmistakable:

1. theorem statements are exact under declared models;
2. the P1 gate is presented before any reconstruction result;
3. the full physical sweep is exhaustive over the frozen axes rather than a selected gallery;
4. the Fisher engine has an independent Monte Carlo validation;
5. all positive and negative comparators are equal-ledger;
6. source code, manifests, scene banks and figure tables are public at submission;
7. the experimental-transfer paragraph is a falsifiable protocol, not evidence.

State “simulation-only” in the abstract or final Methods sentence. Do not use an optical-bench rendering that could be mistaken for completed hardware.

---

# 9. Writing order

Do not draft the Introduction first.

1. **Lock the pocket branch.** Insert either the achieved inset or the hollow estimator-gap marker; no other architecture changes.
2. **Build Fig. 1 first.** The manuscript is not allowed to proceed until the master atlas is intelligible without a caption.
3. **Build Figs. 2–3 from frozen source tables.** No prose values are typed by hand; every plotted number has a machine-readable source.
4. **Write captions and one-sentence claims.** If a panel needs more than one claim, split or delete it.
5. **Write Sections 2–4 around the figures.** Theory first, then the living mean layer, then the covariance boundary.
6. **Write the Discussion and Introduction last.** They may summarize only what the figures already establish.
7. **Assemble S1–S10 at full depth.** Preserve all negative results, correction history, preregistration gates, and code provenance there.
8. **Run the disclosure audit.** Check that P1, fresh-pattern dominance, and DLGI C2 appear exactly once in main text and are not softened elsewhere.
9. **Port to the target journal template and enforce 7 pages at 100% print scale.** Any overflow is solved by cutting, not by shrinking figures or captions.
10. **Freeze a submission-ready repository release.** Experimental contact begins only after the manuscript, supplement, code, data manifests, and transfer protocol are complete.

---

# Final ruling

> **The fallback flagship is a three-figure information-atlas paper, not the R36 mean-only paper and not a two-layer program review. Its protagonist is the gap between formal support and Fisher-usable information. M1 is the living calibration point; the covariance aperture, P1 prelaunch kill, REGION_EMPTY map, and modest pocket define the second layer. If the pocket is achieved, it is one inset; if it is not, the estimator gap is disclosed rather than hidden. Submit to JOSA A or TCI unless the pocket gives Optics Express a genuine positive optical anchor.**