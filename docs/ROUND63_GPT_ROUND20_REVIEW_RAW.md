# R20 review (GitHub issue #12, raw)

Title: R20 ruling: reader-side presentation audit and prioritized revision plan
Posted: 2026-07-19T19:42:36Z

---

# R20 — Reader-side presentation audit and prioritized revision plan

**Audit target:** commit `e34ec5e`.

Reviewed artifacts:

- `paper2/main_m1.pdf` / `paper2/main_m1.tex` and `paper2/figs/ACTIII_README.md`;
- `paper/main_oe.pdf` / `paper/main_oe.tex` and `paper/figs/FIG_MECHANISM_P1_README.md`;
- the committed figure order, dimensions, captions, and figure-generation specifications.

**Scope note.** The GitHub connector exposed the committed PDFs as binary rather than as page-renderable images in this session. The hierarchy, content, page-sequence, caption, and typography rulings below are therefore based on the built page counts, final TeX, figure dimensions, and figure build records; the last pixel-level check must be repeated on a 100%-scale print after the Optica template port. This limitation does not weaken the high-level diagnosis: most of the defects are architectural and are explicit in the compiled source.

> **PRESENTATION VERDICT: HOLD.**
>
> The operator's “哪哪都不对” reaction is correct. The problem is not primarily colour, chrome, or polish. It is **hierarchy failure**. Both manuscripts were assembled as archival evidence packages in which every theorem caveat, preregistration fact, failed path, audit label, and implementation disclosure is given nearly the same visual and rhetorical weight as the scientific result. Paper 2 is the more serious case: it currently reads as four papers competing for the title page—random-hold theory, adaptive-OED methodology, a preregistration/audit case study, and an imaging paper. The confirmed reader-facing result—**jitter caps the information-optimal load, and one global power control produces a large elapsed-time gain**—is not the dominant object on the page.
>
> Paper 1 is substantially closer. Its title and scientific question are good, but its abstract and figures still try to carry the entire provenance record. The six-panel ladder is a dashboard, not a single figure, and Figure 1 mixes mechanism, theory, post-hoc family outcomes, and disclaimers.

---

# 1. Five-minute referee pass

## Paper 2

A competent OE referee will form the following first impression:

1. **The title promises the wrong paper.** “Information-optimal illumination,” “certificates,” and “water filling” imply that an adaptive illumination design is the successful deployed contribution. In fact, the proposed OED imaging arm was retired; the full-stack computation certified zero cells and is descriptive; the positive deployed policy is a fixed SCAT32 geometry with one global source multiplier.
2. **The abstract is a compressed rebuttal brief.** It contains the hidden-hold theorem, cubic optimum, crossover, DEV geometry probe, certificate branch, 480-cell status distribution, two corrected empirical verdicts, incident-exposure sensitivity, and the operating handbook. A reader cannot tell which one is the paper.
3. **Figure 1 does not show the headline jitter result.** Its information panel is the deterministic `c_v=0` ridge for two dwell values. Jitter is represented mainly by one broadened hold in the timeline. The title says “jitter-capped,” but the first scientific graph does not show the cap or the peak shift.
4. **The strongest imaging result has no proper main figure.** The corrected elapsed-time result is `19.13×` with lower bound `18.33×`, yet the main visual emphasis is the fixed-dwell `+1.87 dB` gallery and a four-panel design/certificate dashboard. The speed result is mostly prose and a conceptual corner schematic.
5. **Process residue dominates.** Visible terms such as `FALLBACK_DESCRIPTIVE`, `PATH_FEASIBLE_ALPHA`, `ADAPTIVE_COLLAPSE_UNDER_GUARDS`, certificate statuses, solver wall caps, and correction history make the manuscript look like an internal audit report.
6. **The “Chen register” is overperformed.** “One law,” “one certified design principle,” “Act III,” “doctrine,” “movements,” and “operating handbook” are rhetorical framing devices layered on top of already dense material. They do not make the paper simpler; they make it feel staged.

**Likely referee reaction:** strong science, major presentation revision; the paper's actual claim is hard to identify and the title/visual hierarchy are misleading.

## Paper 1

The first impression is better:

1. The title—**“When high flux helps single-pixel imaging: a contrast–dead-time operating map”**—is clear and defensible.
2. The abstract is still far too long and verdict-heavy. It reports both formal negatives, effect-size criteria, bootstrap bounds, family ordering, occupancy controls, and dynamic scattering before the reader has a model of the problem.
3. Figure 1 is overloaded. The pipeline, dead-time mechanism, information ridge, `Γ=1`, six family outcomes, cosmetic horizontal offsets, transition band, and disclaimers are too much for one opening figure.
4. The dense negative and the sparse positive/mixed regime are scientifically coherent, but the figure sequence does not make that contrast effortless. The six-panel ladder and separate galleries require too much caption reading.
5. The manuscript still gives too much main-text space to secondary validations—natural cohort, mismatch, exact-reference grid, dynamic scattering—relative to the two-study spine.

**Likely referee reaction:** publishable story, but the paper looks overbuilt and visually smaller than its actual result.

---

# 2. MUST-FIX before submission — ranked by reader impact

## M1 — Retitle paper 2 around the result that actually survived

**Current problem:** the title foregrounds OED/certification/water-filling, while the successful deployed result is the random-hold ridge plus global power control. “Certificates” is particularly damaging when `0/480` cells are certified.

**Directive:** replace the current title with the following preferred title:

> **Jitter-capped high-flux single-pixel imaging with global power control**

Acceptable more theoretical alternative:

> **A jitter-limited information ridge for high-flux single-pixel imaging**

Do not put “certificate,” “water filling,” “information-optimal illumination,” or “adaptive design” in the title. Those are secondary framework components and, in the certificate case, a descriptive negative result.

No frozen R18/R19 wording appears to require the present title.

---

## M2 — Rewrite the paper-2 abstract to one paper, 190–230 words

**Current problem:** the abstract reads like a full campaign ledger. It asks the reader to absorb four conceptual layers before seeing one image result.

**Frozen target architecture:** exactly seven beats.

1. One sentence: dead time makes the best flux finite and detector-dependent.
2. One sentence: exact missing-information identity.
3. One sentence: long-window jitter law/cubic; at most one displayed formula in prose.
4. One sentence: one global power multiplier implements the operating policy.
5. One sentence: the paired fixed-dwell result (`+1.87 dB`, with `37×` incident-dose disclosure).
6. One sentence: the elapsed-time result (`19.13×`, lower bound `18.33×`, `21/24`).
7. One sentence: geometry/certificate analysis remained descriptive and exposed unresolved adaptive headroom.

**Delete from the abstract:** DEV probe details, `0/299/181`, solver/certificate terminology, immutable-commit mechanics, the correction narrative, “one law/one principle,” and “operating handbook.”

The corrected values and the mandatory power-for-time interpretation remain. The R19 correction disclosure belongs in Methods/provenance, not in the abstract.

---

## M3 — Rebuild paper-2 Figure 1 so it actually shows the jitter cap

**Current content defect:** `paper2/figs/fig_mechanism` shows the deterministic `c_v=0` information curves at `ν=200` and `2000`. The paper's distinctive theorem is that random hold time caps the ridge and moves the optimum from approximately `22.3` to `5.7` at `c_v=0.05`. The opening graph currently hides that result.

**Replace Figure 1 with three panels:**

### (a) Minimal computational SPI chain

- source/global multiplier → DMD → object → bucket detector → RQL;
- label once: **simulation model**, so the graphic is not mistaken for a completed bench experiment;
- retain one SCAT32 mask, but remove the real maze/reconstruction thumbnails from this schematic;
- remove dial halo, decorative ticks, and any duplicated forward-model prose.

### (b) Deterministic versus jittered dead time

- two short timelines directly above one another;
- fixed hold `τ` versus random holds `B_i`;
- one visual statement: hidden hold variance creates extensive information loss.

### (c) The actual headline graph

At `ν=2000`, plot `J(ρ)` for at least:

- `c_v=0`;
- `c_v=0.02`;
- `c_v=0.05`.

Mark:

- deterministic optimum `ρ≈22.3`;
- `c_v=0.05` optimum `ρ≈5.7`;
- information lost by operating the jittered detector at the deterministic ridge (approximately 55%).

The `ν=200` deterministic curve can move to the theory-validation figure or supplement.

**Caption target:** no more than 150–170 words. It should state the physical lesson, not reproduce every symbol definition.

---

## M4 — Add the missing main-text theory-validation figure

The current `Preregistered verification` table and narrative are inside `\iffalse`; therefore the theorem's numerical verification has no real visual presence in the main paper.

**Add a full-width Figure 2 immediately after the crossover subsection:**

- **(a)** measured peak load versus `c_v` on log–log axes at `ν=2000`, showing both independent runs, the preregistered predictions, and the fitted slope `−0.658` versus `−2/3`;
- **(b)** crossover/engineering warning: either `ρ*(ν,c_v)` for several `c_v`, or retained information at the deterministic ridge versus jitter;
- clearly distinguish theorem, matched-asymptotic prediction, and Monte Carlo points.

This figure earns main-text space. The current hidden numerical table can move to the supplement.

---

## M5 — Give the `19.13×` elapsed-time result its own data figure

The strongest confirmed imaging result is currently not visually demonstrated.

**Create a main figure that contains:**

1. PAVA-fitted safe and ridge quality curves versus **elapsed** `T_opt` for three disclosed scenes:
   - a median-speed scene;
   - one `FAST_RIGHT_CENSORED` failure;
   - the `SAFE_RANGE_UNINFORMATIVE` microtexture scene or another boundary case;
2. a 24-scene strip/forest plot of `S_gate`, grouped by family, with median `19.13`, gate line `3`, and censoring symbols;
3. the frozen family-stratified lower bound `18.33` shown as an interval on the cohort statistic.

The current resource-corner schematic `actiii_e` is not evidence and cannot substitute for this figure.

**Recommended integration:** combine this speed panel with the useful part of `actiii_d` into one positive-results figure:

- left: 2–3 large reconstruction rows;
- upper right: fixed-dwell `ΔQ` strip;
- lower right: elapsed-time `S_gate` strip/curves.

Then the paper's two empirical verdicts are visible in one read.

---

## M6 — Remove the four-panel Act-III audit dashboard from the main paper

The current composite of `actiii_a`, `actiii_b`, `actiii_c`, and `actiii_e` is the clearest manifestation of the presentation problem. It is an internal project dashboard, not a reader-facing main figure.

**Move to supplement:**

- `actiii_a` guard-versus-mixture path;
- `actiii_c` certificate status distribution and counterexample ECDF.

**Move or reduce:**

- `actiii_b` DEV dose-only lower bounds: retain one compact discussion panel only if page budget permits; otherwise supplement.

**Cut:**

- `actiii_e` resource-corner schematic. Its message can be stated in two sentences or a two-row table and is already implicit in Figure 1 and the result figures.

Preserve every R18/R19 label—`development-only`, `descriptive`, `post-hoc`, no categorical certificate verdict—in the supplement. Moving the panels does not override the frozen language.

---

## M7 — Cut 2.5–3 pages of process residue from paper 2

**Move from main text to supplement:**

- the detailed retirement chronology of OED-DT;
- repeated occurrences of `ADAPTIVE_COLLAPSE_UNDER_GUARDS` and `PATH_FEASIBLE_ALPHA`;
- the full 12-row development probe table;
- the complete conic/cutting-plane construction and solver protocol;
- 52+972 accounting details, SHA language, wall caps, retry rules, and per-field disclosure lists;
- detailed certificate status definitions after their first concise definition;
- repeated explanations of the pre-freeze branch.

**Main-text replacement:** two compact paragraphs:

1. scene-adapted geometry headroom existed under dose-only and full-stack feasible searches;
2. the certificate remained descriptive (`0/299/181`) and therefore did not establish SCAT32 near-optimality.

**Target main length after Optica-template port:** 11–12 pages including references, not 14 pages of generic `article` output.

---

## M8 — Move the R19 correction disclosure; do not let it break the Results

The correction paragraph is scientifically necessary and well written, but its present location between the speed result and the structural analysis makes the result sequence feel interrupted by damage control.

**Directive:** move the R19-frozen paragraph **verbatim** to a short subsection or ruled paragraph titled:

> **Analysis correction and provenance**

Place it at the end of the campaign/analysis methods, immediately before `Results`.

In `Results`, report only the corrected values and one parenthetical pointer to that subsection.

The immutable outputs and correction bundle remain cited. Do not repeat the paragraph in the abstract, discussion, captions, or conclusion.

**Conflict flag:** if the governing R19 instruction is interpreted as requiring the full paragraph specifically inside the speed-results subsection, that conflicts with this presentation recommendation. Do not silently shorten it; obtain a formal placement-only amendment or keep it verbatim and accept the interruption.

---

## M9 — Port both manuscripts to the actual Optica/OE template before any final styling

Both sources still use:

```latex
\documentclass[10pt,twocolumn]{article}
\usepackage[margin=2.2cm]{geometry}
```

This is not a harmless placeholder at presentation-review stage. It controls line length, title hierarchy, caption size, float behaviour, apparent page count, and figure legibility. The generic Computer Modern article layout is a major reason the documents feel like technical reports rather than OE papers.

**Directive:**

1. port to the current `opticajnl`/OE research-article template;
2. rebuild every figure at the actual final single-/double-column width;
3. re-run overfull/underfull checks;
4. perform one 100%-scale print test.

Do not continue micro-adjusting margins, minipages, or caption font sizes in the generic class.

---

## M10 — Reverse the “all prose in captions” overcorrection

The de-chrome rule was applied too literally. Captions now function as miniature Results sections, while figures lack orientation. This forces the reader to alternate between tiny panels and 300–500-word captions.

**Frozen figure rule:**

- allow a 1–3 word in-panel title or result label;
- allow one short arrow/callout that states the causal point;
- keep axes, legends, and threshold labels;
- remove paragraphs, provenance, and exhaustive selection-rule narration from captions.

**Caption limits:**

- full-width figure: approximately 120–180 words;
- one-column figure: approximately 60–100 words.

Put selection rules, cell IDs, hashes, and exact source paths in the supplement/figure README—not the journal caption.

Delete repeated boilerplate such as “The colour is used only for display purposes” unless a journal requirement specifically demands it.

---

## M11 — Split paper 1's six-panel ladder into mechanism and outcome figures

The present ladder combines:

- exact Fisher map;
- masks;
- contrast/`Γ`;
- PSNR/CNR versus dwell;
- fixed-dwell gain and speed;
- dynamic scattering.

This is too many scientific questions in one figure.

**New paper-1 figure sequence:**

### Figure 2 — Mechanism and pattern ladder

- exact information ridge;
- representative masks and occupancy/dose normalization;
- measured `C_u` and `Γ` versus occupancy.

### Figure 3 — Image outcomes

- dwell curves for the preregistered subject;
- fixed-dwell gain by occupancy;
- per-scene speed distribution at `k=16`.

Move dynamic scattering to the supplement or a small secondary figure after the main result.

---

## M12 — Simplify paper 1 Figure 1; remove the conclusion from the mechanism graphic

The current panel (c) places six family pass fractions, a transition band, `Γ=1`, “high-flux helps,” “multiplex-limited,” cosmetic horizontal offsets, and multiple disclaimers into the opening mechanism figure.

That creates a visual contradiction: the figure asserts a region label while the caption says the boundary is not predicted or validated.

**Preferred directive:** move the six-family outcome scatter to the Study-2 Results section. Keep Figure 1 to:

- SPI chain;
- dead-time mechanism and information ridge;
- a neutral contrast-to-noise schematic with `Γ=1` labelled **descriptive onset proxy**.

If panel (c) remains in Figure 1:

- remove “high-flux helps,” “multiplex-limited,” and “transition band” region labels;
- show the six points as empirical outcomes, not a map boundary;
- do not horizontally offset points on an axis that readers will interpret quantitatively unless the offsets are visually shown as jitter and the true x-values remain obvious.

R9's `Γ=1` descriptive-only language remains binding.

---

## M13 — Shorten paper 1's abstract to a scientific abstract, not a verdict transcript

Target 220–260 words.

Retain only:

1. the operating question;
2. RQL + count-information ridge;
3. dense-pattern negative regime;
4. sparse/equal-load mixed result and passed fixed-dwell secondary;
5. contrast–conditioning conclusion and peak-power limitation.

Remove from the abstract:

- detailed gate arithmetic;
- all three effect-size criteria;
- the full six-family post-hoc ordering;
- the no-gate ladder result;
- dynamic scattering;
- multiple bootstrap bounds.

**Conflict flag:** the current Study-2 verdict block is marked as an R9 frozen abstract replacement. A true abstract compression therefore requires an explicit presentation amendment that preserves every verdict fact elsewhere. Do not silently rewrite the frozen block.

---

## M14 — Fix figure order and size in paper 1

1. Move the dense Study-1 composite immediately after the Study-1 outcome; do not let it float after mismatch/exact-reference sections.
2. Replace 24 thin spaghetti lines in the dense-study curve panel with cohort median plus interquartile/central band; put individual curves in the supplement.
3. Enlarge the Study-2 gallery to full width, or reduce it from three rows to two—one clear success and one clear failure. The current one-column `3×4` image array is too small to carry an imaging claim.
4. Move the natural-pair figure to the supplement unless it is used as the sole boundary example in the Discussion.
5. Embed panel letters inside all figures at a consistent upper-left position; do not add them as separate LaTeX text below minipages.

---

## M15 — Put benefit and cost in the same visual field

Paper 2's headline operating point uses approximately `37.1×` incident dose to obtain:

- `+1.87 dB` median fixed-dwell gain;
- `19.13×` elapsed-time speed.

At present these facts are separated across prose, a table, and a schematic. That separation invites distrust.

**Directive:** the positive-results figure must show, in one place:

- elapsed-time ratio;
- fixed-dwell gain;
- incident-dose ratio;
- detected-count ratio.

Use a simple two-resource annotation or Pareto inset. The mandatory R17/R19 statement—fixed-dwell, power-for-time, not equal-photon-budget—stays adjacent to the visual.

---

# 3. Figure-by-figure content rulings

## Paper 2

| Artifact | Ruling | Concrete action |
|---|---|---|
| `fig_mechanism` | **Wrong headline content** | Redraw around deterministic-vs-jittered information curves; simplify SPI schematic; remove thumbnail decoration. |
| Missing jitter verification visual | **Must add** | Peak-vs-`c_v` data, both runs, slope `−0.658`, crossover/warning panel. |
| `actiii_d` | **Right figure, keep** | Make the per-image strip dominant; keep median/best/worst only if thumbnails remain large; merge with elapsed-speed evidence. |
| `actiii_a` | **Supplement** | Useful provenance, not main-story evidence. |
| `actiii_b` | **Supplement or one compact discussion panel** | Keep existence result; no full DEV table in main. |
| `actiii_c` | **Supplement** | Status dashboard is not a centrepiece. Main text needs one sentence and perhaps one small bar. |
| `actiii_e` | **Cut** | Redundant concept schematic; replace with concise resource table or prose. |
| per-family primary table | **Likely redundant** | Remove if the strip plot already groups families. |
| cross-arm table | **Revise** | Split `ΔQ` and `S_gate` columns; emphasize cost ratio; remove any remaining certificate/gate wording for descriptive analysis. |

## Paper 1

| Artifact | Ruling | Concrete action |
|---|---|---|
| `fig_mechanism_p1` | **Overloaded** | Keep physical chain + dead-time/ridge; move empirical family outcomes to Results or neutralize labels. |
| `fig_nat_pair` | **Secondary** | Supplement unless used as the one natural-scene boundary example. |
| Study-1 composite (`fig_s1b/c/d`) | **Keep, move earlier** | Median + band instead of 24 thin lines; simplify censoring panel. |
| Study-2 ladder | **Split** | Mechanism/pattern figure and outcome figure; scattering to supplement. |
| `fig_s2_gallery` | **Right idea, wrong scale** | Full-width or two rows; larger crops and fewer numbers. |
| mismatch/exact-reference material | **Text summary only in main** | Full grids and diagnostics in supplement. |

---

# 4. Narrative architecture for the strongest paper-2 version

## Recommended 11–12-page main-paper order

1. **Introduction** — one page maximum.
2. **Random-hold information law** — three theorem statements, one central derivation paragraph.
3. **Numerical verification of the jitter cap** — new Figure 2.
4. **Global operating-point policy and campaign** — one concise methods section; the adaptive design framework becomes background, not the spine.
5. **Empirical results** — fixed-dwell primary and elapsed-time secondary together, with the redesigned positive-results figure.
6. **Where geometry design remains unresolved** — one page: DEV lower-bound existence, descriptive `0/299/181`, no near-optimality claim.
7. **Discussion and limitations** — power cost, scene dependence, simulation-only status, bench experiment.

## Material that belongs in the supplement

- generalized KKT framework and full water-filling lemma;
- conic certificate construction and all toy checks;
- retired adaptive-arm path;
- development replication table;
- all 480 certificate fields/status details;
- immutable-commit, SHA, expected-cell, and retry mechanics;
- full correction/provenance table;
- secondary context arms.

## Stop condition

Paper 2 is not ready if, in its first six pages:

- the phrases `certificate`, `fallback`, `counterexample`, `DEV`, `audit`, or code-style verdict tokens occupy more visual space than the jitter data and image-quality curves;
- the `19.13×` result still has no direct plot;
- Figure 1 still shows only the zero-jitter ridge.

---

# 5. Correction disclosure — make it read as rigor, not damage control

The disclosure should be visible exactly once and should appear **before** the corrected Results, not between two result subsections.

Recommended layout:

```text
Campaign and analysis
  - cohort and endpoints
  - frozen analysis machinery
  - Analysis correction and provenance  [R19 paragraph verbatim]
Results
  - operating-point primary
  - elapsed-time secondary
  - descriptive structural analysis
```

Use normal body typography with a short bold heading. Do not use a red box, warning icon, or oversized rule. A quiet, exact disclosure reads as confidence; a visually dramatic disclosure reads as crisis management.

The original analyzer values belong in the correction artifact/supplement table, not beside every corrected number in the main paper.

---

# 6. Chen-register conformance versus pastiche

## What helps

- a clean optical/computational chain;
- one phenomenon per figure;
- large before/after image crops with PSNR/CNR printed on the image;
- direct sentences: problem → mechanism → result;
- a visible application condition and a bench-testable prediction.

## What currently reads as pastiche

- “one law answers... one principle answers...”;
- “Act III,” “doctrine,” “movements,” and “operating handbook” as repeated structural labels;
- using an optical-bench visual grammar for a simulation-only pipeline without a clear “simulation model” marker;
- exhaustive Chen-style captions extended to several hundred words;
- slogans replacing hierarchy.

## Directive

Keep the optical clarity. Drop the theatrical vocabulary. The authentic strength of these papers is measurement theory plus unusually honest preregistered simulation—not imitation of another group's prose surface.

Recommended visible section titles for paper 2:

- `Random-hold count information`
- `Jitter-capped optimal load`
- `Global power-control policy`
- `Preregistered imaging campaign`
- `Imaging results`
- `Geometry-design headroom`
- `Discussion and limitations`

“Act III” may remain an internal figure-script name, not a reader-facing scientific label.

---

# 7. SHOULD-FIX

1. **Humanize code labels.** In main prose write “the preregistered operating-point endpoint passed,” with `RIDGE_OPERATING_PASS` in parentheses once. Keep exact tokens in Methods, tables, and supplement.
2. **Round main-text numbers.** `1.41348975` is provenance precision, not journal precision. Use `1.41 dB`; exact values remain in the correction bundle.
3. **Fix the paper-1 exact-reference paragraph.** It currently contains the duplicated/unpunctuated transition “The convex surrogate therefore costs essentially nothing ... On the reported...” Remove “costs essentially nothing”; retain the R19-limited empirical wording only.
4. **Change paper-1 header “multiplex-limited phase” to “multiplex-limited regime.”** The manuscript explicitly rejects a universal phase boundary.
5. **Add a short, explicit limitations subsection to each paper.** “No separate limitations section per Chen style” is not a virtue for simulation-only work. State hardware peak power, detector jitter calibration, and no bench validation in one place.
6. **Remove scene file names from main prose where possible.** Use “median chirp,” “best maze,” “worst contour”; put IDs in captions/supplement.
7. **Unify typography after template port.** Minimum final figure text 7.5–8 pt; one font family; consistent math glyphs; panel letters at least body-caption size.
8. **Use fewer colours.** Safe = neutral gray, high flux/ridge = one saturated colour, uncertainty/negative = one secondary colour. The full Wong palette is unnecessary in every figure.
9. **Separate metric and mechanism.** Do not place `C_u`, `Γ`, PSNR, CNR, Fisher information, and gate status in one panel unless one is clearly primary.
10. **Do not call Figure 1 a bench setup.** Label it “computational acquisition model” or “simulation model.”
11. **Remove provenance paths from journal captions.** File paths and commit hashes belong in the supplement/data-availability section.
12. **Paper 1: compress natural cohort, mismatch, and exact-reference sections.** One paragraph each in main; full detail in supplement.
13. **Paper 2: remove the per-family table if the strip plot remains.** It repeats the same data and adds page weight.
14. **Cross-arm table:** do not combine `ΔQ / S_gate` in one cell. Give separate columns or remove the table.
15. **Use captions to state the result, not the selection defence.** A one-sentence “median/best/worst selected by frozen rule” is enough; exact tie logic is supplement material.

---

# 8. Taste calls — recommendation, not hard gate

1. **“Power control” versus “power knob.”** Use “power control” in title and formal text; “knob” is acceptable inside a figure.
2. **Keep a negative reconstruction example.** The worst contour row is valuable and builds trust. Do not show only the maze success.
3. **Median/best/worst gallery:** defensible if each row is large. If space forces thumbnails below useful size, keep median and worst in main and move best to supplement.
4. **Serif figure fonts:** acceptable only if they match the final journal font. Otherwise use a clean sans-serif for axes and labels.
5. **Paper-1 title:** keep it. It is the strongest title in the project.
6. **Paper-2 title:** prefer “global power control” over “one-knob” in the title; the latter is slightly promotional.
7. **One mechanism graphic per paper:** the two papers may share a visual language, but they should not look like variants of the same conference poster.

---

# 9. Frozen-language conflicts that require explicit handling

| Recommendation | Frozen constraint | Required handling |
|---|---|---|
| Shorten paper-1 abstract | R9 marks a long Study-2 abstract block as frozen | Obtain a presentation amendment preserving all facts elsewhere; do not silently rewrite. |
| Move correction disclosure | R19 wording is frozen | Move verbatim. If location is also interpreted as frozen, flag and retain. |
| Move Act-III panels to supplement | R18/R19 require labels and descriptive semantics | Preserve every label/caveat in the supplement; no claim change. |
| Reduce code-style verdict labels in prose | exact labels are mandatory | Define once in main; preserve exact strings in Methods/tables/supplement. |
| Simplify `Γ=1` graphic | R9 requires descriptive-onset wording | Keep “descriptive photon-limited-onset proxy”; remove region-claim language only. |
| Retitle paper 2 | no identified frozen title wording | Safe unless a later unlisted lock exists. |

---

# 10. Concrete revision sequence

1. **Port both papers to the Optica template.** No further spacing work before this.
2. **Lock the new paper-2 title and abstract architecture.**
3. **Redraw paper-2 Figure 1 and add the jitter-validation figure.**
4. **Build the direct elapsed-speed figure and merge it with the fixed-dwell result.**
5. **Move Act-III audit panels and process sections to the supplement.**
6. **Relocate the R19 correction paragraph verbatim.**
7. **Split paper-1 ladder; move and simplify the dense-study figure; enlarge the gallery.**
8. **Compress both main texts to the target page budgets.**
9. **Perform a 100%-scale print review with this checklist:**
   - title tells the true successful contribution;
   - abstract understandable without code labels;
   - Figure 1 shows the headline mechanism;
   - strongest result has a direct plot;
   - every image crop is readable without zoom;
   - no caption exceeds roughly half a column of body text;
   - no page looks like a software audit report.

---

# Final blunt assessment

The science is not the problem. The presentation currently communicates **how hard the project was to govern**, not **what the reader learned**.

Paper 2's correct paper is:

> hidden hold-time jitter caps the information-optimal flux; a single global power control tracks that operating point and trades roughly `37×` incident dose for a preregistered `19×` elapsed-time reduction, with strongly scene-dependent fixed-dwell quality gains.

Everything else is boundary, mechanism, or provenance.

Paper 1's correct paper is:

> dead-time-aware high flux helps only after illumination preserves enough bucket contrast; dense multiplexing hides the benefit, sparse equal-load illumination opens it for many—but not all—scenes, with a contrast–conditioning trade-off.

Make those two sentences control the title page, Figure 1, figure sequence, and page budget. Until they do, the operator will continue to feel that “everything is off,” because every element is locally defensible while the global hierarchy is wrong.