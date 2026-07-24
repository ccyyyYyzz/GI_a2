# BUILD STATUS — PRL Letter, Phases 1–2 (R43 FINAL LOCK)

**Paper:** *Strong Optical Disorder Erases Images but Preserves Local Change Observability*
(Rank–Jet Separation). Binding spec: `docs/ROUND63_GPT_ROUND43_RULING_RAW.md`.
**Scope of this build:** R43 §7.3 steps 1–3 — claim-source matrix + the three main figures,
built ONLY from committed frozen artifacts. No data generation (§7 green light honored).
**Repo HEAD:** `18891c65` (branch `main`).

---

## Done (Phase 1 deliverables)

| # | Deliverable | File(s) | Status |
|---|---|---|---|
| 1 | Claim–source matrix | `paper_prl/CLAIM_SOURCE_MATRIX.md` | ✅ complete; every R43 §4.1–4.6 + figure-spec number mapped to a committed artifact (value→file→field→commit) |
| 2 | Figure 1 (hero / the gate) | `figures/fig1_hero.{py,pdf,png}` | ✅ 3 design passes + self-critique (log below) |
| 3 | Figure 2 (jet_slopes exponent hierarchy) | `figures/fig2_jet_slopes.{py,pdf,png}` | ✅ rebuilt at PRL width from committed `JET_TEST` (1bf29f1) |
| 4 | Figure 3 (two optical realizations + honest boundary) | `figures/fig3_optical_realizations.{py,pdf,png}` | ✅ from committed `CONFIRMATORY` (b37c841) + `SCRAMBLE` (ed7a1e0) |
| 5 | This status file | `paper_prl/BUILD_STATUS.md` | ✅ |

**Provenance discipline.** All figure scripts read ONLY from
`figures/_frozen_sources/*.committed.json`, extracted verbatim from the frozen commits
(`JET_TEST`@1bf29f1, `SCRAMBLE`@ed7a1e0, `CONFIRMATORY`@b37c841). Each script header records
its source path and commit. A deterministic consistency check confirmed all **22 hard-coded
label literals** match the committed JSON exactly (KL 2.04/4.00, coef 1.07%/0.00%, MC 0.95/2.05,
CUSUM −2.16/−3.92, 77/81, 453 banks, LCB 0.990, 459/1013, bal-acc 0.916, TPR 0.988, FA
0.096/0.084, duel d′ 4.1 / AUC 0.997, rank-one 99.99%, mean-null 6.6e-19).

**Publication-grade compliance (R43 §6 + task):**
- vector PDF + 600-dpi PNG for all three figures ✅
- `pdf.fonttype=42` verified (each PDF embeds TrueType `FontFile2` subsets) ✅
- colorblind-safe Okabe-Ito palette (`figures/prl_style.py`) ✅
- widths: Fig 1 & 3 double-column (180 mm), Fig 2 double-column; all axis labels/ticks/titles ≥ 7 pt ✅
- secondary micro-annotations (schematic m=1/2 tags, sub-captions) sit at 6.5–6.6 pt — standard for
  double-column sub-labels; can be enlarged or moved to the caption if PRL production enforces a hard 7-pt floor.

---

## Figure 1 self-critique log (R43 §7.3 step 2: "if not legible in ten seconds, revise the figure")

**Pass 1** — three rows (mean support disk→X'd {0}; rank many-axis ellipse→rank-1 line; jet
slope-ladder→boxed "SURVIVES") + disorder arrow + scrambling-law callout + drawn anchor icon.
*Critique:* (i) right-state labels "{0} (DC only)" and "rank 1 / 99.99%" bled into the callout box,
overlapping its equations; (ii) two grey brace connectors cut through text — clutter; (iii) the
collapsed rank-1 line was vermillion = same color as the surviving jet m=2 class — a semantic clash.

**Pass 2** — moved the right-state column left (clean gutter before the box), relocated right
labels to sit under/over their glyphs, dropped the braces, recolored both collapsed glyphs
neutral-dark so the *dying* channels read as "lost" and only the *surviving* jet row keeps vivid
blue+vermillion. *Critique:* collisions gone and the dead-vs-alive contrast now works, but the two
disorder regimes were not visually grouped and the callout did not read as a governing law.

**Pass 3 (final)** — added two regime background bands (thin-screen grey | complete-scrambling
light-blue) to structure the axis, gave the callout a blue accent header + rule and a green
"blind set empty" pill, and balanced the callout's vertical rhythm. Bumped the smallest secondary
labels toward the 7-pt floor. *Result:* the 10-second read is immediate — support→{0} and rank→1
collapse (dark), while the quotient jet order (m = 1 or 2) survives (green), governed by
`Cov(b)=Q(x)(O∘O)+R`, `Q=xᵀGx`, with the amplitude-anchor condition shown as an icon.

All R43 §6.2 required elements are present: horizontal disorder axis; the three rows with the exact
transitions `B_p→{0}`, `many→1`, `finite m=1/2 → finite m=1/2`; the endpoint factorization; the two
directional classes; and the amplitude-anchor icon (not a caption caveat).

---

## Gaps flagged (NOT regenerated — R43 §7 absolute)

**GAP-1 (CRITICAL, needs a coordinator ruling).** `results/round63_next/JET_TEST/` is **modified in
the working tree** relative to the frozen commit `1bf29f1` that R43 references. The uncommitted
re-run's numbers differ from R43 §4.6 (KL 2.015 vs frozen 2.038; MC 1.17/1.99 vs 0.95/2.05; CUSUM
−2.19/−3.98 vs −2.16/−3.92; iso-Q AUC 0.516 vs 0.481). **Resolution applied:** Figure 2 and all
§4.6 claims use the **committed frozen** JSON (`1bf29f1`), extracted to
`figures/_frozen_sources/JET_TEST.committed.json`; the working-tree re-run is not used anywhere.
Both value sets support the identical qualitative conclusion (integer orders 2/4, MC ≈1/2, CUSUM
≈−2/−4, empty blind set). **Coordinator action:** either `git checkout -- results/round63_next/JET_TEST/`
to restore the frozen artifact (recommended — a post-freeze `jet_test.py` re-run is a §7.2-forbidden
regeneration), or formally re-freeze the re-run and update R43 §4.6 decimals. Full table in
`CLAIM_SOURCE_MATRIX.md` GAP-1.

**GAP-2 (non-blocking).** Fig 2's crossover marker uses the committed `bank_A.crossover` block
(mixed direction, pred `ε_cross`=0.009995 vs empirical 0.009997, ratio 1.000). Recorded only to note
the marker is the *mixed* direction, not `delta_g`.

**GAP-3 (non-blocking).** End-Matter theorem proofs (Rank–Jet, quotient-jet) are prose/equation
sources in `docs/ROUND63_GPT_ROUND42_RULING_RAW.md` §6 and `SCRAMBLE_EXT/DERIVATION.md` §2–4; not yet
transcribed into the manuscript. No missing evidence — this is Phase-2 scope.

**No numerical claim in R43 §4.1–4.6 or the figure specs lacks a committed source.**

---

## Phase 2 — LaTeX build (R43 §7.3 steps 4–6) — COMPLETE

Built in R43 §7.3 order: theorem core + End Matter **before** Introduction, then the 3.65-page
core, then S1–S5 as an audit trail. TeXLive 2024 (`pdflatex`), REVTeX 4.2 (`revtex4-2.cls`).

| # | Deliverable | File | Status |
|---|---|---|---|
| 6 | Main Letter (core + End Matter A/B/C) | `paper_prl/main.tex` → `main.pdf` | ✅ compiles clean |
| 7 | Supplemental Material S1–S5 | `paper_prl/supplement.tex` → `supplement.pdf` | ✅ compiles clean |

### Compiled page counts (100% scale, PRL two-column)

- **`main.pdf` = 6 pages.** Core = **pages 1–4** (page 4 is Figs 2–3, top-floated); End Matter A/B/C
  = **page 5 + 3 lines of page 6**; references = **rest of page 6**. → **core ≤ 4 pp, End Matter ≈ 1.1 pp**,
  both inside the R43 budget (4-page core + 2-page End Matter).
- **`supplement.pdf` = 5 pages** (S1–S5, Tables S1–S5).

### Word count (texcount)

- **Core** (abstract + body + 3 figure captions, excl. End Matter & refs): **1379 text + 269 caption
  = 1648 words**. PRL-formula estimate incl. 3 full-width figures + 5 display equations
  ≈ **2.7–2.9k words**, well under the **3750** ceiling.
- End Matter text: **604 words**. Supplement: **1474 text + 238 caption + 5 tables**.
- **Abstract: ≈ 588 characters** (adapted from R43 §1.3 seed; ≤ 600 char PRL limit).

### Compile health

- **No overfull hboxes, no undefined references/citations** in either document (both fixed:
  eq. `jet`/`classes` split to fit column; one supplement inline limit promoted to a display).
- **One residual underfull hbox** in `main.tex` (badness **1122**, sealed-instantiation paragraph) —
  a single mildly-loose line, far below LaTeX's 10000 problem threshold and routine in REVTeX; not a
  layout defect.

### Deviations from the R43 §6.1 page map

- **Core composes to ~3.0 pp of text + 1 pp of figures (Figs 2–3 float to p.4).** The per-part
  proportions (opening 0.40 / jets 0.80 / Rank–Jet 0.55 / inversion 0.55 / tests 0.60 / sealed 0.35 /
  conclusion 0.20) are followed in **content and ordering**; the composed text is **lighter** than
  3.65 pp because the three double-column figures consume ~1 full page. This is *under* budget, not
  over — acceptable, and leaves headroom for the user's author/affiliation/acknowledgment blocks.
- **Exactly 5 displayed equations in the core** (`jet`, `rank`, `Q`+ΔQ, `classes`, `mean`) — at the
  R43 cap.
- **Frozen wording used verbatim** where it fits: §4.6 validation sentence (KL 2.038/4.000, MC
  0.95/2.05, CUSUM −2.16/−3.92), §4.2 robustness-domain sentence, §4.1 capability compressed with the
  D3/D5 numbers preserved in body + Fig. 3 caption (no "specific"/"zero-false-alarm"/"medium-immune").
  End Matter C carries both integrity disclosures **verbatim** (R43 §4.5 chunking; R41 §6 shot
  correction).
- **End Matter placement:** after the core + acknowledgments, before the shared bibliography (so both
  core and End Matter resolve the same numbered references). A `\clearpage` starts End Matter on a
  fresh page per convention.

### Every number traces to CLAIM_SOURCE_MATRIX.md

All sentence-level values in both documents were re-checked field-by-field against the committed
frozen JSON (`JET_TEST`@1bf29f1, `SCRAMBLE_EXT`@ed7a1e0, `CONFIRMATORY`@b37c841). No new analysis; no
softening or strengthening of frozen wording. GAP-1 resolution unchanged: the Letter uses the
**committed frozen** JET_TEST values; the uncommitted working-tree re-run is not used anywhere.

---

## TODO — user-level items (clearly-marked placeholders in the source)

These are the only non-automatable blanks. Each is a literal placeholder token in `main.tex` /
`supplement.tex`:

1. **`\author{[AUTHORS]}`** — author block (both files).
2. **`\affiliation{[AFFILIATION]}`** — institutional affiliation(s) (both files).
3. **`[FUNDING/ACKNOWLEDGMENTS]`** — the `acknowledgments` environment in `main.tex` (funding,
   grant numbers, thanks).
4. **Data-availability statement + URL** — add a short statement pointing to the GI_a2 repository
   (`github.com/ccyyyYyzz/GI_a2`) with the frozen commit hashes; recommended placement is the end of
   End Matter C or a dedicated back-matter note. *(Not yet inserted — needs the user's preferred
   public URL/DOI wording.)*
5. **Reference verification** — the 11 main + 5 supplement `\bibitem` entries are best-effort
   bibliographic reconstructions of the R41/R43 head-on prior art (Feng/Lee memory effect, Lee–Stone
   UCF, Tsang PRX, Mudry blind-SIM, Idier, Chaigne, Snieder coda-wave, Dechant, Lorden, Hermann–Krener,
   BBP). Volume/page/year should be verified against the user's Zotero before submission.

## Two-reads checklist for Phase 3 (R43 §7.3 step 7)

1. **General-physics interest read** — does the rank-vs-jet inversion land in ten seconds from
   title + abstract + Fig. 1 alone, without ghost-imaging background?
2. **Hostile statistical-optics read** — scope/priority: are the three theorem boundaries (locality,
   amplitude anchor, monotone cone) stated where a referee will look, and is the D3/D5 boundary
   impossible to misread as a specificity claim? Confirm the prior-art fence (S5) pre-empts
   memory-effect / SPADE / coda-wave equivalence objections.

Phase 2 timebox (~4 h) met; both documents compile from source with a single `pdflatex ×2` pass each.

---

## Phase 3 — pre-submission fix-pass (R43 §7.3 step 7) — COMPLETE

Both pre-submission reads (`REVIEW_READ1_GENERAL.md` F1–F10; `REVIEW_READ2_HOSTILE.md`
B1, M1–M5, m1–m7) were adjudicated by the coordinator and every accepted repair applied.
No new data was generated: figures were re-rendered from
`figures/_frozen_sources/*.committed.json` via the existing `fig*.py` scripts with
**label/layout/caption changes only** — plotted data untouched. `CLAIM_SOURCE_MATRIX.md`
was not modified (coordinator-owned; m1/M4 handled there).

### Edits applied (finding ID → file:line, post-edit)

| ID | Fix | Location |
|---|---|---|
| B1 / A1 | Abstract: "sealed bucket-optics test" → "sealed **simulated** single-pixel optical test"; also F2 gloss pass (see below) | `main.tex:31–40` |
| B1 | "sealed single-pixel optical **experiment**" → "sealed, **fully simulated** single-pixel optical **model**" | `main.tex:89` |
| B1 / item 3 | "thin-screen **experiment on a true** single-pixel bucket record" → "**simulated** thin-screen **study of the** single-pixel bucket record" | `main.tex:289` |
| B1 / item 4 | Fig. 3 caption: inserted "**simulated**" before sealed thin-screen detector | `main.tex:272–286` |
| M1 / item 5 | empty-blind-set headline qualified: "…empty **on the aperture where $G\succ0$ (grain band covering the scene band)**" | `main.tex:171` |
| M1 / item 6 | "no blind direction for detection **within the informative aperture**" | `main.tex:221–222` |
| M1 / item 7 | boundaries paragraph: added premise sentence ($G\succ0$ needs grain modulus-squared band to cover scene band; change beyond $2k_{\rm grain}$ ⇒ $u=v=0$, genuinely blind, Suppl. case iii) | `main.tex:184–189` |
| M2 / item 8 | "mean wall supplies **specificity**" → "mean wall **rejects DC and flux changes**"; noun "specificity" now absent from theory section (only D3 "…not a fully **specific** sentinel" remains, permitted) | `main.tex:223` |
| M5 / item 12 | Fig. 2 caption: CUSUM slopes flagged as two-point confirmations of the densely sampled exact-KL exponents, not independent fits | `main.tex` (Fig. 2 caption) |
| m4 / item 13 | End Matter A: added "(so $d=v/2$ …)" reconciling $v=2\delta^\top G\delta$ with $d=\delta^\top G\delta$ | `main.tex:366` |
| m6 / item 15 (+F10) | "**unique** statistic that **saturates**" → "**attains**"; sentence split; triple name-drop lightened | `main.tex:119–121` |
| m2 / item 17 | Outlook: Lorden cited for "sequential change detection"; Snieder moved to "disordered and medium-change sensing" | `main.tex` (Outlook) |
| m3 / item 18 | "measured covariance energy … single direction $Gx$" → "**scene-dependent** covariance energy … **rank-one structure $O\!\circ\!O$**" (body + Fig. 1 caption aligned) | `main.tex:216`, `:55` (Fig.1 cap) |
| F1 / item 20 | Fig. 1 `figure*` moved to the front matter (after `\maketitle`); now typesets at **top of page 2**, co-visible register with abstract (was p.3) | `main.tex` (front matter) |
| F2 / item 21 | Abstract gloss pass (see deviations) | `main.tex:31–40` |
| F3 / item 22 | sealed section: "cells (**medium-parameter grid points**)" gloss added; "banks" already glossed at `:73` ("one independent medium realization"); rot-20%/slope−1 FA + latency-19% **kept** (part of frozen §4.2 robustness sentence) | `main.tex:292` |
| F4 / item 23 | Fig. 3 caption rewritten to **lead with the physics panel** (complete scramble: image dies, detection lives); `fig3` panel order swapped (physics‖thin-screen‖calibration), data untouched | `main.tex:272–286`, `fig3_optical_realizations.py` |
| F5 / item 24 | `fig1_hero.py`: row label + banner "quotient jet order" → "**local change order $m$**"; added "$O\!\circ\!O$: code overlaps; $R$: noise floor" gloss in law box; regenerated | `fig1_hero.py` |
| F6 | body: "exact Gaussian divergences" → "exact Gaussian **statistical (KL and Chernoff)** divergences" (abstract had no char budget for this gloss) | `main.tex:252` |
| F7 / item 19 | `fig2_jet_slopes.py`: main-panel slope label "2.04" → "**2.038**"; regenerated | `fig2_jet_slopes.py` |
| F8 / item 25 | Eq. (5) glossed: "$O_{ii}$ the code weight and $\hat{x}(0)$ the scene DC component" | `main.tex:206` |
| F9 / item 26 | added pointer "—the reduced law is analytic in $Q$ (End Matter A)—" | `main.tex:158–159` |
| M3 / items 9–10 | Supplement S5 prior-art fence rebuilt: Mudry/Idier re-attributed to **blind-SIM/speckle-SR** (not memory-effect); memory-effect imaging (Bertolotti/Katz/Freund), Grace–Guha quantum change detection, SAR-CCD (Preiss&Stacy), and Lee–Stone/UCF fenced; 6 new bibitems | `supplement.tex:417–433`, `:455–461` |
| m5 / item 14 | Supplement S2: Siegmund (Springer 1985) cited for ARL₀→threshold; bibitem added | `supplement.tex:143`, `:461` |
| M5 / item 12 | Supplement S2: CUSUM two-point-confirmation sentence added | `supplement.tex:143–146` |
| m7 / item 16 | Supplement S1: exponent $m$ robust to unknown nuisance; constant is a composite-hypothesis/least-favorable quantity for which profiled Chernoff is a bound | `supplement.tex:84–89` |
| item 11 | prior-art bibitems added to **supplement only** (main refs kept lean; no main-text sentence required them) | — |
| m1, M4 | **NOT applied here** — `CLAIM_SOURCE_MATRIX.md` edits, coordinator-owned (per task) | — |

### Frozen-wording deviations (adjudicated)

R43 §4.1–4.6 frozen sentences were preserved verbatim **except** the following surgical
edits, all explicitly adjudicated in the fix list. No frozen numerical value or claim was
altered; the §4.6 validation sentence (KL 2.038/4.000, MC 0.95/2.05, CUSUM −2.16/−3.92) and
the §4.2 robustness sentence (FA 0.032, latency 19%, 0.052, 0.076) are untouched in the body.

1. **Abstract** (deviates from the R43 §1.3 *seed*, not a §4 frozen sentence): reworded to
   (a) label the sealed test "**simulated**" [A1]; (b) gloss the jet concept —
   "first nonzero Chernoff jet…detection-time exponent" → "the lowest order at which a change
   perturbs the statistics (its jet order) fixes the detection time" [F2]; (c) "nuisance
   profiling" → "projecting out nuisances" [F2]; (d) **dropped** "amplitude-anchored" and
   "first-order-" (before "orthogonal") for the 600-char budget [F2, permitted]; (e) "Exact
   divergences" gloss and explicit T/ε gloss **omitted** — no char budget (priority: sim-label
   > jet gloss > rest). Final rendered length **598 chars** (≤ 600). "**All results are
   numerical.**" was **not** added — no char headroom (task: "if chars allow").
2. **§4.1 lead-in framing** (not the frozen numbers sentence): "experiment on a true
   single-pixel bucket record" → "simulated … study of the single-pixel bucket record" [B1/3];
   added "cells (medium-parameter grid points)" gloss [F3].
3. **Fig. 3 caption** rewritten to lead with physics + "simulated" inserted + D2 capability
   numbers (77/81, 453, 0.990, 459/1013) **de-duplicated** out of the caption (they remain in
   the body's frozen §4.1 sentence and on-figure) [F4, item 22].
4. **Theory sentence** `:223`: "supplies specificity" → "rejects DC and flux changes" [M2].
5. **§67 / §271 capability framing**: "optical experiment" → "fully simulated optical model" [B1].

### New compile stats (post-fix, `pdflatex ×2` each)

| Check (task item 27) | Result |
|---|---|
| Overfull hboxes (both docs) | **0 / 0** — PASS |
| Underfull hboxes (main) | **0** — PASS |
| Undefined refs / citations (both docs) | **none** — PASS (7 new supplement bibitems all resolve) |
| Core ≤ 4 pp | **PASS** — core = pages 1–4 (text p.1–3; Figs 2–3 float to p.4); End Matter starts p.5 |
| Main total ≤ 6–7 pp | **PASS** — `main.pdf` = **6 pages** |
| Abstract ≤ 600 chars | **PASS** — 598 (source-rendered) |
| texcount estimate ≤ 3750 | **PASS** — core ≈ 2.7 k PRL-formula words (core text 1367 + captions 286 + 5 eqns + 3 full-width figs); full-doc texcount Sum = 2639 |
| Exactly 5 displayed core equations | **PASS** — eq:jet(1), eq:rank(2), eq:Q(3), eq:classes(4), eq:mean(5); End Matter uses inline math only |
| Fig. 1 on page 1–2 | **PASS** — top of **page 2** |
| Supplement | 6 pages (was 5; +1 from the extended S5 fence + 6 bibitems), 0 overfull, no undefined |

**Not committed** (per task). All figure PDFs/PNGs regenerated (fig1, fig2, fig3). One
finding intentionally not fully realized: the abstract's "All results are numerical." clause
and the divergence/T-ε glosses were dropped for the 600-char limit (priority-ordered, as
instructed) — the sim-label and jet gloss, the higher-priority items, are in.
