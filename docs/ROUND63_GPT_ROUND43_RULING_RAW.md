# R43 ruling (GitHub issue #34, raw)

Title: R43
Posted: 2026-07-24T06:26:34Z

---

# R43 — Final architecture: rank collapse without observability collapse

**Reference request:** [`docs/ROUND63_GPT_ROUND43_QUESTION.md`](https://github.com/ccyyyYyzz/GI_a2/blob/198632e/docs/ROUND63_GPT_ROUND43_QUESTION.md) at [`198632e`](https://github.com/ccyyyYyzz/GI_a2/commit/198632e31e67e7123006565fbd25e7c6bf0dfbe7).  
**Decisive exams:** `JET_TEST/` at [`1bf29f1`](https://github.com/ccyyyYyzz/GI_a2/commit/1bf29f137e3f05f0175d2d5e87455cd26f1b6d3b) and `SEALED_DET/` at [`b37c841`](https://github.com/ccyyyYyzz/GI_a2/commit/b37c841adc4910edb480f8af6d873f37ac49a872).

## Executive ruling

> **Choose a third cut: a JETS-centered PRL Letter sharpened into a Rank–Jet Separation paper.**
>
> **Strong optical disorder can collapse reconstruction to rank one without creating a local blind direction: after nuisance quotienting, every nonzero scene change first leaves the statistical law at order one or two, which fixes `epsilon^-2` or `epsilon^-4` detection-time scaling.**

This is stronger and cleaner than either proposed center:

- It is broader than a sentinel paper because it identifies a new invariant of a statistical experiment: the **order of local observability** after nuisance quotienting.
- It is more physical than an abstract jets paper because complete optical scrambling supplies an extreme realization: the mean sees only DC, covariance imaging collapses to rank one, yet local change observability survives.
- It is not weakened by the D3/D5 specificity dents. Those dents concern finite-threshold attribution calibration in one thin-screen implementation. They do not touch the quotient-jet theorem, its exact contact orders, the scrambling inversion, or the certified detection-power result.

The paper's protagonist is therefore not the old information atlas, M1, DLGI, or even the detector itself. It is:

> **reconstruction rank and change observability are different invariants.**

The sealed detector is the severe optical instantiation and engineering reality check. Its power endpoints enter; its false-alarm leakage enters with equal prominence; it does not carry the title or the theorem.

---

# 1. The center of gravity

## 1.1 Preferred title

> **Strong Optical Disorder Erases Images but Preserves Local Change Observability**

This is the best PRL title because it states the physical paradox without introducing project acronyms or a new technical noun in the title.

Acceptable alternatives:

1. **Rank Collapse without Local Blindness in Scrambled Light**
2. **Curvature-Rescued Detection through Complete Optical Scrambling**
3. **Statistical Observability beyond the Imaging Limit**

Do not use `DLGI`, `fog-DMD`, `information atlas`, `sentinel`, `quotient jet`, or any R-number in the title. `Curvature-rescued detection` may be named inside the article after the effect is defined.

## 1.2 Frozen thesis

> **After profiling all admitted nuisances, the first nonzero statistical jet of a physical perturbation fixes its observability class: in an amplitude-anchored scrambled covariance channel, reconstruction Fisher rank can collapse to one while generic changes remain first-order observable and first-order-orthogonal changes are rescued by curvature, producing `T_req ~ epsilon^-2` and `T_req ~ epsilon^-4` laws.**

## 1.3 Abstract seed

> Strong disorder is expected to erase spatial information. We show that reconstruction rank and local change observability are distinct invariants. After nuisance profiling, the first nonzero Chernoff jet of a perturbation fixes its detection-time exponent. In an amplitude-anchored fully scrambled optical channel, the scene Fisher matrix collapses to rank one, yet every nonzero local change remains visible: generic directions require `T~epsilon^-2`, while first-order-orthogonal directions are curvature-rescued with `T~epsilon^-4`. Exact divergences, Monte Carlo tests, sequential detection, and a sealed bucket-optics experiment validate the result.

## 1.4 Why this is PRL-shaped

The broad statement is not “we built a better optical alarm.” It is:

> **identifiability dimension can vanish while local hypothesis distinguishability remains complete at a higher contact order.**

That is relevant to disordered sensing, singular statistical models, sequential change detection, inverse problems, and quantum-inspired superresolution—not only to ghost imaging. The optical system gives an unusually clean physical realization and an exact falsification chain.

---

# 2. The sharper theorem: Rank–Jet Separation

This theorem should be named in the Letter. It is the conceptual bridge between the R42 jets and the scrambling three-row inversion.

Let the marginal record law depend on the scene through `r` smooth functionals

```text
q(x) = (q_1(x),...,q_r(x)),
P_x = P_{q(x),eta},
```

with every admitted nuisance `eta` profiled out. Let `J_q` be the efficient Fisher matrix in functional space. Then

```text
I_x = Dq_x^T J_q Dq_x,
rank(I_x) <= r.                                                  (2.1)
```

Thus `r=1` permits at most rank-one reconstruction information, regardless of how many detector samples are acquired.

For a physical direction `delta`, write the quotient expansion

```text
q(x+epsilon delta)-q(x)
 = epsilon u(delta) + (epsilon^2/2) v(delta) + o(epsilon^2),     (2.2)
```

where `u` and `v` are projected away from the nuisance tangent. Under regularity of the reduced law:

```text
u != 0             -> C_*(epsilon;delta) = c_1 epsilon^2 + o(epsilon^2),
                       T_req ~ epsilon^-2;

u  = 0, v != 0     -> C_*(epsilon;delta) = c_2 epsilon^4 + o(epsilon^4),
                       T_req ~ epsilon^-4.                       (2.3)
```

For the fully scrambled covariance channel

```text
Sigma(x)=R+Q(x)H,
Q(x)=x^T Gx,
```

with medium amplitude anchored and `G>0`,

```text
u(delta)=2x^T G delta,
v(delta)=2delta^T G delta.                                      (2.4)
```

Hence every `delta != 0` is in exactly one of two local classes:

- `x^T G delta != 0`: first-order, `m=1`;
- `x^T G delta = 0`: curvature-rescued, `m=2`, because `delta^T G delta>0`.

Therefore:

> **The reconstruction Fisher matrix is rank one, yet the local blind set is empty.**

This is the **Rank–Jet Separation Theorem**. It is the exact structural result the two completed exams jointly reveal.

### Necessary boundaries

1. The statement is **local**. Signed finite changes can return to the same `Q` ellipsoid at `epsilon_*=-2x^TGdelta/(delta^TGdelta)`; JET_TEST correctly recovered AUC `≈0.5` there.
2. With an unanchored covariance-amplitude nuisance, the scene and amplitude derivatives are collinear and the quotient information collapses to zero. A wavelength, lag-shape, polarization, or finite calibration anchor is required.
3. Global sign-definiteness holds on the monotone cone when `x,delta>=0` and the kernel is entrywise nonnegative; it is not claimed for arbitrary signed changes.

The main article should state all three in compressed form. They make the theorem stronger, not weaker, because its exact domain is known and experimentally falsified at its boundaries.

---

# 3. One-paper discipline

## 3.1 Formal legitimacy versus strategic ruling

A simultaneous PRL Letter plus a longer Physical Review companion is formally legitimate. APS explicitly permits joint submissions when the longer article contains considerably more or different information and substantially improves understanding ([PRL author guidance](https://journals.aps.org/prl/authors); [editorial on joint submissions](https://journals.aps.org/prl/edannounce/10.1103/PhysRevLett.118.240001)).

**Do not use that option here.**

It would violate the operator's one-paper discipline in substance even if it complied with APS policy:

- the putative companion would mostly unpack the same atlas and sentinel lifecycle rather than introduce a second independent physical result;
- it would revive the anthology structure R35 correctly killed;
- it would invite referees to treat the Letter as an advertisement for a methods paper;
- it would split the clean `rank collapse / jet survival` protagonist from its evidence.

## 3.2 Binding publication unit

> **Submit one self-contained PRL Letter, with concise End Matter and a focused Supplemental Material deposit. No simultaneous companion.**

PRL currently permits a `3750`-word core—about four journal pages—and up to two pages of End Matter for specialist material. The Letter must remain convincing without the supplement. The formal supplement should contain only material necessary to audit this Letter:

- complete quotient-jet and Rank–Jet proofs;
- JET_TEST methods and all 17 checks;
- SCRAMBLE_EXT derivation and numerical generator;
- SEALED_DET protocol, D0–D6 tables, mismatch boundaries, and chunking audit;
- the minimum thin-screen derivation needed to connect the two disorder endpoints.

The full M1, DLGI, MOLT, five-certificate, living-region, and project-history archives remain public in the repository and are cited through the data/provenance statement. They do **not** become a hundred-page PRL supplement.

A later experimental paper is legitimate because it would contain new physical evidence. A later memory-effect phase-diagram paper is legitimate only if it establishes new physics between the two endpoints. Neither is a companion to this submission.

---

# 4. Final frozen wording under D3/D5

## 4.1 Detection capability

> **In a preregistered thin-screen test, the nuisance-profiled covariance alarm met every detection-power and deployment endpoint: `2%` beyond-band changes were predicted detectable in `77/81` cells, the best cell required `453` banks and achieved a Monte Carlo power lower bound of `0.990`, repeated banks retained a `459`-versus-`1013`-bank latency advantage over fresh covariance codes, and online scoring consumed `0.16%` of each bank interval. At the frozen `1%` alarm threshold, four-class balanced accuracy was `0.916` and in-band false alarms were `0.020`, while medium-amplitude and medium-timescale events produced false-alarm rates of `0.096` and `0.084`; the result is therefore a power-certified beyond-band alarm with measured attribution leakage, not a fully specific sentinel.**

The main Letter may shorten this to one sentence plus a figure caption, but it may not replace the measured false-alarm numbers with “specific,” “zero-false-alarm,” or “medium-immune.”

## 4.2 Robustness domain

> **All frozen performance and calibration bars passed in the primary mismatch domain of aperture-preserving medium-basis rotations up to `10%` (`FA=0.032`, latency inflation `19%`). At `20%` rotation and a spectral-slope mismatch of `-1`, non-target false alarms rose to `0.052` and `0.076`; these points define the measured robustness boundary, and no broader law-mismatch claim is made.**

## 4.3 D3 disclosure

> **D3 did not weaken detection power, but it killed the strong specificity claim: target TPR was `0.988`, the beyond score remained centered on non-targets (`|d'|<=0.02`), and balanced accuracy was `0.916`, yet medium-amplitude and timescale tails exceeded the frozen `0.05` false-alarm bar (`0.096` and `0.084`) and the weakest class-specific separation was `1.79` rather than `5`.**

## 4.4 D5 disclosure

> **D5 preserved AUC `1.000` and bounded latency inflation on every tested mismatch axis, but false-alarm calibration failed outside the primary domain—`0.052` at `20%` rotation and `0.076` at spectral slope `-1`—so robustness is frozen to the `<=10%` rotation regime and the broader wording is withdrawn.**

## 4.5 Generator-chunking integrity disclosure

> **After protocol freeze and before unblinding, the record generator was refactored solely to stream Monte Carlo records in GPU-memory chunks; an audited diff found no changes to hypotheses, thresholds, bank counts, scenes, estimators, or endpoints, although altered random-number consumption produced a different, still-blinded realization set. All reported decisions apply the frozen analysis unchanged to that disclosed chunked stream.**

This belongs in Methods/End Matter and the repository provenance note, not in the abstract.

## 4.6 Quotient-jet theorem statement for the main text

> **Let `C_*(epsilon;delta)=inf_eta' C(P_{x+epsilon delta,eta},P_{x,eta'})` be the Chernoff distance after nuisance profiling. If `C_*=c_m(delta)epsilon^(2m)+o(epsilon^(2m))` with `c_m>0`, then fixed-error sample complexity obeys `T_req asymp epsilon^(-2m)`. In the amplitude-anchored scrambled covariance law `Sigma(x)=R+[x^TGx]H`, `G>0`, generic directions have `m=1`, first-order-orthogonal directions have `m=2`, and no nonzero local direction is blind, even though the scene Fisher matrix has rank one.**

Follow immediately with the measured validation:

> **Exact KL orders were `2.038` and `4.000`, Monte Carlo response exponents were `0.95` and `2.05`, and CUSUM delay slopes were `-2.16` and `-3.92`; the predicted crossover, nuisance collapse, anchor restoration, iso-`Q` cancellation, and monotone-cone corollary all passed their frozen tests.**

---

# 5. Final venue ruling and sequencing

## 5.1 First submission: Physical Review Letters

> **Submit to PRL first, with the Rank–Jet Separation cut.**

The case is now materially stronger than in R41. PRL asks for a concise, complete result combining impact, innovation, and broad interest; its current Letter limit is `3750` words, approximately four core pages, with up to two pages of End Matter ([official author guidance](https://journals.aps.org/prl/authors)). The manuscript now has:

- a reparameterization- and nuisance-invariant theoretical object;
- an exact theorem separating reconstruction dimension from change observability;
- integer contact orders and universal latency exponents;
- complete-scrambling and thin-screen optical endpoints;
- exact-KL, Monte Carlo, CUSUM, gauge, cancellation, and cone tests;
- a sealed practical instantiation with its failed calibration bars disclosed.

Simulation-only is not a conceptual defect for this cut because the paper is a theorem plus controlled numerical falsification, not an unverified experimental-performance claim.

## 5.2 Second submission if PRL declines

If PRL declines at editorial screening for breadth or because the result is judged too optical, recut the **same single paper** for **Optica**:

- lead with the beyond-band bucket alarm;
- retain Rank–Jet Separation as the theoretical spine;
- expand the thin-screen power map, fixed-versus-fresh comparison, and measured specificity boundary;
- use the R41 `<=7 pp / 4 figures` architecture.

Optica explicitly targets high-impact fundamental and applied optics and photonics ([journal scope](https://opg.optica.org/content/journal/about/item/optica/)).

If PRL sends the paper for review and the reports affirm the physics but judge it too specialized, an APS transfer to Physical Review Applied is acceptable if offered; otherwise proceed to Optica. **Do not maintain simultaneous active submissions and do not create a companion before the first decision.**

Nature Photonics is not the first submission for a simulation-only result. It becomes credible only after a physical demonstration or a substantially broader scrambling-channel experiment.

## 5.3 PRL justification paragraph

Use a version of the following in the required impact statement:

> Strong disorder is expected to destroy spatial sensing. We prove and numerically validate a different outcome: in a scrambled covariance experiment, reconstruction Fisher rank collapses to one while every nonzero local scene change remains observable at a finite nuisance-profiled contact order. The first nonzero Chernoff jet fixes universal `epsilon^-2` or `epsilon^-4` detection-time laws, verified by exact divergences, Monte Carlo statistics, and sequential detection. A sealed bucket-optics experiment then demonstrates the effect under severe band limitation. The result separates identifiability dimension from observability order and applies broadly to disordered sensing, change detection, and nonregular statistical experiments.

---

# 6. Final PRL page map and figures

Target **3.55–3.75 core journal pages**, below the official approximately four-page/`3750`-word ceiling, plus at most two pages of End Matter. Use **three main figures**.

## 6.1 Core page map

| Core element | Single purpose | Target |
|---|---|---:|
| Title + abstract | State rank collapse without local observability collapse. | `0.20 pp` |
| **Opening: the disorder paradox** | Strong scattering erases images; ask whether it must erase change observability. | `0.40 pp` |
| **Quotient information jets** | Define profiled contact order and derive `T_req~epsilon^-2m`. | `0.80 pp` |
| **Rank–Jet Separation** | Prove rank `<=r`; specialize to rank-one `Q=x^TGx`, `m=1/2`. | `0.55 pp` |
| **Scrambling inversion** | Mean support `->{0}`, covariance scene rank `->1`, local jets survive. | `0.55 pp` |
| **Decisive tests** | Exact/MC/CUSUM slopes, crossover, nuisance collapse/anchor, cancellation/cone. | `0.60 pp` |
| **Sealed optical instantiation** | Power-certified detector plus measured D3/D5 boundaries. | `0.35 pp` |
| Conclusion | General principle and exact scope. | `0.20 pp` |
| **Total** |  | **`3.65 pp`** |

No M1 result, five-certificate strip, DLGI lifecycle, 81-cell estimation tombstone, reconstruction gallery, or full prior-art atlas enters the PRL core.

## 6.2 Figure 1 — Hero: rank collapses; jets survive

A horizontal disorder axis from thin multiplicative screen to complete scrambling, with three rows:

```text
mean spatial support        B_p  -------------------------->  {0}
reconstruction score rank   many -------------------------->   1
quotient jet order          finite m=1/2 ------------------> finite m=1/2
```

At the scrambling endpoint, show the factorization

```text
Cov(b)=Q(x)(O o O)+R,
Q=x^TGx,
```

and the two directional classes. Include the amplitude-anchor condition as a visible icon/label, not a caption caveat.

**Ten-second conclusion:** disorder can erase coordinates without erasing local change observability.

## 6.3 Figure 2 — The PRL evidence figure: `jet_slopes`

Use the existing log-log asset as the center and add only the minimum insets:

1. exact KL/Chernoff contact orders `2` and `4` with coefficient agreement;
2. Monte Carlo `d'` slopes `1` and `2`;
3. CUSUM delay slopes `-2` and `-4`;
4. a compact crossover marker at the predicted `epsilon_cross`;
5. a two-symbol nuisance inset: amplitude gauge `J=0`, anchor restores `J>0`.

Do not turn this into a six-panel dashboard. The log-log exponent hierarchy is the figure.

## 6.4 Figure 3 — Two optical realizations and the honest boundary

One full-width figure with three compact parts:

- **Thin screen:** sealed `2%` power result (`77/81`, best `453` banks, MC LCB `0.990`) and fixed-versus-fresh latency (`459` versus `1013`).
- **Calibration boundary:** balanced accuracy `0.916`, medium-event false alarms `0.096/0.084`, and the `<=10%` rotation domain. Show this as an honest threshold graphic, not fine print.
- **Complete scrambling:** mean AUC at chance, covariance response/rank-one functional, and the generic `5%` duel (`d'=4.1`, AUC `0.997`) as the opposite physical endpoint.

The figure's single sentence is:

> **The order law survives both weak and complete scrambling; practical power survives, while attribution calibration has a measured domain.**

## 6.5 End Matter and Supplemental Material

**End Matter A:** proof of Rank–Jet Separation and the quadratic specialization.  
**End Matter B:** nuisance quotient, amplitude-gauge collapse, and minimal anchor.  
**End Matter C:** model/protocol capsule and integrity disclosures.

Focused Supplemental Material:

- **S1:** general Chernoff/KL jet derivation and equality cases;
- **S2:** JET_TEST design, all `17/17` checks, raw exponent tables;
- **S3:** full-scramble Gaussian derivation, generator, rank-one result, and duel;
- **S4:** sealed thin-screen protocol, D0–D6, fixed/fresh comparator, D3/D5 failures, chunking audit;
- **S5:** thin-screen covariance aperture and the minimum prior-art/provenance needed to connect the endpoints.

All older program material remains in the cited repository rather than being copied into the PRL supplement.

---

# 7. Writing green-light protocol

> **GREEN LIGHT. No further campaign, probe, parameter sweep, estimator repair, threshold calibration, mismatch test, or physics extension may enter this manuscript.**

The memory-effect interpolation and any other in-flight line are quarantined as future work. They do not delay writing and cannot alter the R43 title, theorem, figures, or claims.

## 7.1 Allowed work

- typeset and simplify existing derivations;
- render figures only from committed frozen tables and images;
- perform deterministic consistency checks of equations, labels, and source-table joins;
- verify references and priority wording;
- run exact reproduction checks only when they are expected to reproduce a frozen artifact and cannot create a new endpoint;
- audit code/provenance and write the data-availability statement.

## 7.2 Forbidden work

- new random seeds or scene banks;
- a D3/D5 threshold repair;
- additional specificity classes;
- a memory-effect campaign;
- an experimental emulation presented as evidence;
- a new title-driving result discovered after this freeze;
- replacing a failed endpoint with a newly invented endpoint.

If a factual or code error is found, stop and disclose/correct it. Do not compensate by launching a replacement campaign.

## 7.3 Build order

1. **Freeze the claim-source matrix:** every sentence-level numerical claim mapped to one committed artifact.
2. **Build Figure 1 first:** if rank/jet separation is not legible in ten seconds, revise the figure, not the story.
3. **Lock Figures 2–3** from frozen sources and print-test at final PRL width.
4. **Write the theorem core** and End Matter before Introduction or Discussion.
5. **Write the 3.65-page core** to the PRL word budget.
6. **Build S1–S5** as an audit trail, not an alternate manuscript.
7. **Run two final reads:** one general-physics read for PRL interest and one hostile statistical-optics read for scope/priority.
8. **Submit.** No “one more test” after the hostile read.

---

# 8. Innovation unlocked by the two exams

## 8.1 Rank–Jet Separation is the sharper cut

The two candidate manuscripts in the question still treat theory and detector as separate pillars. The exams reveal one deeper object:

> **the rank of parameter reconstruction and the order of directional observability can separate arbitrarily.**

In the fully scrambled endpoint, the separation is maximal:

```text
reconstruction rank = 1,
local observable directions = all delta != 0,
jet order = 1 generically, 2 on the orthogonal class.           (8.1)
```

This is more surprising than “detect what cannot be imaged.” It says why that sentence is possible and provides a classification that transfers beyond optics.

## 8.2 The Disorder Observability Triangle

A statistical sensing experiment should be classified by three independent objects:

```text
support: where the law can respond;
rank:    how many parameter coordinates can be reconstructed;
jet:     at what contact order a chosen physical direction becomes detectable. (8.2)
```

Strong disorder can increase support, decrease rank, and leave jet order finite. No single scalar information budget captures this. Figure 1 is the first physical realization of this triangle.

This is an outlook object, not another main theorem. It should appear in the final paragraph as the generalization invited by the Letter.

## 8.3 Cone-complete observability

The iso-`Q` cancellation test prevents a global signed-change claim. The monotone-cone result supplies a clean positive corollary:

> if `x>=0`, `delta>=0`, and `G` is entrywise nonnegative with `G>0`, then `Q(x+delta)-Q(x)>0` for every nonzero `delta`.

Thus local completeness becomes global on physically monotone addition processes. This is a useful bridge to deposition, coating growth, contamination, and damage accumulation without overstating arbitrary signed-change detection.

These three innovations are already supported by the frozen derivations and exams. They require no new campaign.

---

# Final lock

> **One PRL Letter. One protagonist: Rank–Jet Separation. Three figures. No companion. No further data generation.**
>
> The core result is that strong disorder can erase images without erasing local change observability; quotient information jets determine the price. The sealed detector proves that the effect is not algebra alone, and its D3/D5 failures define—not undermine—the honest engineering boundary.