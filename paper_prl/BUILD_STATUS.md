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
