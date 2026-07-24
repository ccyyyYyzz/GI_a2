# BUILD STATUS — Paper 2 (Optica): *Engineering What an Optical Sensor Cannot See*

**Program:** ROUND63 · Certified-Blindness / Three-Ledger line.
**Binding spec:** `docs/ROUND63_GPT_ROUND45_RULING_RAW.md` **Part C** (architecture FROZEN).
**Gate:** `results/round63_next/TLSG/` → **TLSG_PASS** (all 7 bars) → §C architecture frozen.
**Repo HEAD at build:** `35d858e` (branch `main`).
**Working dir:** `paper2_optica/` — nothing here touches `paper_prl/` (paper 1 frozen) or sealed artifacts.

---

## Phase log

### Phase 1 — claim matrix + hero figure (THIS PHASE) — COMPLETE

| Deliverable | Status | File |
|---|---|---|
| Claim–source matrix (every quoted number → committed artifact → hash) | **DONE** | `CLAIM_SOURCE_MATRIX.md` |
| Frozen-source extraction (verbatim committed blobs) | **DONE** | `figures/_frozen_sources/*.committed.json` (11 files) |
| Figure style module (Optica; reuse of paper-1 `prl_style.py`) | **DONE** | `figures/o2_style.py` |
| Fig 1 hero "What the sensor cannot see" (3 linked panels) | **DONE** | `figures/fig1_hero.{py,pdf,png}` |
| Fig 2 jet transmutation (spec + built) | **DONE** | `figures/fig2_jet_transmutation.{py,pdf,png}` |
| Fig 3 three routes + calibration transfer (spec + built) | **DONE** | `figures/fig3_routes_transfer.{py,pdf,png}` |
| Build status + Optica limits | **DONE** | this file |

**Discipline held:** all figures read ONLY committed frozen JSONs extracted verbatim to
`figures/_frozen_sources/` (same protocol as paper 1). **No data generation.** All cited
artifact directories verified clean in the working tree.

### Phase 2 — prose (COMPLETE)

| Deliverable | Status | File |
|---|---|---|
| Main article (5 sections, 3 theorems + jet-transmutation corollary, proof sketches) | **DONE** | `main.tex` → `main.pdf` (**6 pp**) |
| Supplement (full proofs S1–S4, negatives S5, methods+provenance S6, TLSG gate table) | **DONE** | `supplement.tex` → `supplement.pdf` (**6 pp**) |

**Compile (pdflatex, TeX Live 2024):** main and supplement both compile clean —
**0 undefined references, 0 missing citations, 0 overfull hbox** (worst residual box
eliminated; only a few benign underfull vboxes remain from full-width float pages).
All three committed figures embed (`fig1_hero`, `fig2_jet_transmutation`,
`fig3_routes_transfer`). All 11 main-text references are cited (no dangling entries).

**Section map (R45 §C3, organized by theorem — no 17-test chronology):**
1. Symmetry walls — Statistical Noether Wall Thm + support wall 6.6e-16 (A1) / energy wall 2.18e-16 (A3) + two breakers (B1 NA-clip ε¹ 0.97–1.03; B2 window 5.3e-2) + wall periodic table (supp S1).
2. Blindness capacity — Thm β_d=σ_d (4.2/4.3), robust envelope (4.4); 80@1e-5 single-plane (C1), 35@1e-4 joint (C2), σ_d curve (C3), C4/C5, fresh-draw min 80 + coverage 1.000 (C6/C7).
3. Restoring the jet — Jet Transmutation corollary + flow eq; Route 1 slope 4.00 (D1), Route 2 86–99% (D3), Route 3 DPSS ≤2.4e-4 (D4); noncomposability (4.6) with **both** →24 instances kept distinct (D5 IT5 SVD→DPSS; D6 IT4b Tukey 93×).
4. What finite data can certify — Three-Ledger Thm 1:√p:p; exponents −0.004/0.444/1.010 (E1), CVs (E2–E4), KT6 boundary blind power 0.051 + z2⁴ (E5/E6), projected jet above channel line (E7).
5. Wave-optics transfer — leak law (F1), floor decomposition (F2), COMSOL 15.5% (F3), M^1.8 + matched-cell 768 (F4/F5), real-DMD leak (F6), PWM constraint (F7).

**GAP flags — all honored in prose/captions (see CLAIM_SOURCE_MATRIX.md):**
- **GAP-1**: field wall cited as committed **6.6×10⁻¹⁶** everywhere (not R45's loose 7e-16). Hero + §1 + supp table all read 6.6e-16.
- **GAP-2**: joint **35@1e-4** (IT4b spectrum) and fresh-draw single-plane **80** (TLSG bar 7) are labelled with distinct provenance in the hero caption, §2, and Fig-3 caption.
- **GAP-3**: the two "→24" instances are separated — §3/supp S2 name IT5 SVD→DPSS (D5) and IT4b Tukey rim-kill 93× (D6) as different operators/splits.
- **GAP-4** (length): Optica Research Article guideline is 6–8 pp. The scaffold typesets at **6 pp** in a dense 10 pt / 0.8 in two-column layout; under `opticajnl` the same body expands toward the upper end of 6–8 pp. No padding was added — the anti-hype PRL house voice is kept. Abstract is ~100 words (Optica target), not the mistaken ≤35.

**Honesty machinery (symmetric):** IT6 segmentation-flood and IT7 coincidence-veto
(24% scene power) appear in supplement **S5** with the same KILL-verdict prominence as
the positives, and S5 is referenced from the main-text Discussion.

**Class-file status:** `opticajnl.cls` is **NOT installed** in this TeX Live 2024
toolchain (`kpsewhich opticajnl.cls` → empty; `opticameeting` also absent). Per the
task fallback, both files are **clean two-column (main) / one-column (supplement)
`article`-class scaffolds** with Optica-matching typography (Times, full-width
title/abstract, numbered sections, spanning `figure*`, ~100-word abstract). A TODO
header in each `.tex` records that the **final submission is a mechanical class swap**
to `\documentclass[10pt,twocolumn]{opticajnl}` + `\journal{opticajnl}` with **no body
edits required**. Same pdflatex toolchain that built `paper_prl`.

**Discipline held:** zero new numbers; every numerical claim traces to a
CLAIM_SOURCE_MATRIX.md row. `paper_prl/`, `results/`, `docs/`, the figures, and the
matrix were **not modified**. Companion Letter cited as `[COMPANION]`; user-reserved
placeholders (`[AUTHORS]`, `[AFFILIATION]`, `[FUNDING/ACKNOWLEDGMENTS]`, `[URL]`) left
intact exactly as in paper 1.

---

## Frozen source provenance (see CLAIM_SOURCE_MATRIX.md for the full map)

| Source | Commit | Role |
|---|---|---|
| `TLSG/TLSG_RESULTS.json` | `35d858e` | three-ledger scaling + calibration transfer (gate) |
| `R44_KILL_TESTS/KT1` | `19338a4` | support/field wall 6.6e-16 + breaker |
| `R44_KILL_TESTS/IT1` | `f13a323` | energy/unitary wall 2.18e-16 + NA-clip ε¹ |
| `R44_KILL_TESTS/KT2` | `05f1029` | jet transmutation (slopes 2→4) |
| `R44_KILL_TESTS/KT3` | `2430166` | statistical projection (overlaps→0, 86–99%) |
| `R44_KILL_TESTS/IT4,IT4b` | `33b8ef2`,`5666d1e` | blindness-capacity spectrum (80/35/24 dims) |
| `R44_KILL_TESTS/IT5` | `c8a0487` | DPSS route 2.4e-4, noncomposability 80→24 |
| `R44_KILL_TESTS/KT6` | `635ff05` | third-ledger boundary (oracle z⁴, blind killed) |
| `R44_KILL_TESTS/IT2,IT3,IT8` | `acfcda0`,`6fb5149`,`6fb5149` | leak-law + PWM apparatus constraint |
| `R44_KILL_TESTS/IT6,IT7` | `dade877`,`e1546ac` | supplement negatives |
| `WAVE_TWIN/*` | `0c3b0d9` | COMSOL 15.5%, M^1.8, matched-cell 768∈[513,2996] |
| `JET_TEST/`, `SCRAMBLE_EXT/` | `1bf29f1`,`ed7a1e0` | cross-paper jet/scrambling anchors (paper-1 frozen) |

---

## Section architecture (R45 §C3) — frozen

1. **Symmetry walls of an optical measurement** — Statistical Noether Wall theorem + two exact walls (§1 / matrix A,B).
2. **Measuring blindness capacity** — Blindness Capacity theorem, dimension–leak frontier, robustness (§2 / matrix C).
3. **Restoring the desired jet** — optical / code-space / statistical routes + noncomposability (§3 / matrix D).
4. **What finite data can certify** — Three-Ledger theorem + KT6 boundary (§4 / matrix E).
5. **Wave-optics transfer and apparatus constraints** — leak law, DPSS/SVD choice, single bucket, PWM (§5 / matrix F).

**No chronology of the 17 tests appears in the article** (R45 §C3). Supplement carries IT6/IT7
negatives, DPSS-vs-SVD noncomposability spectra, étendue and coherence-rank future objects (§C2).

### Frozen theorem set (R45 §C2)
1. Statistical Noether Wall Theorem · 2. Blindness Capacity Theorem (`β_d=σ_d`; `+δ` robustness) ·
3. Three-Ledger Certifiability Theorem (`1:√p:p`) · corollary 4. Jet Transmutation.

### Hero figure self-critique (R45 §C3)
The hero delivers the story in ~10 s to a non-specialist optics reader: a left-to-right causal
spine (**symmetry → capacity → certifiability**) over three quantitative panels — two machine-zero
walls with their breakers; a measured dimension–leak capacity curve (80 dims ≤1e-5, does not
compose); and the `1:√p:p` ladder with KT6 below the blind line and the projected jet above its
channel line. One-line takeaways sit under the two data panels. Verdict: **PASS** (see final report).

---

## Budgets & FORMAT LIMITS (Optica — verified against author resources, 2026-07)

**Venue:** Optica (flagship). Fallback: Physical Review Applied (R45 §C4).

| Item | Optica limit | Source |
|---|---|---|
| **Abstract** | **~100 words** (explicit summary: problem, methods, major results/conclusions) | Optica Publishing Group Journal Style Guide |
| **Research Article length** | **6–8 published pages** | Optica submission guidelines |
| Letter (Optica Express-style short form) | 4 pages | Optica submission guidelines |
| Mini-review | 8–12 pages | Optica submission guidelines |
| Length metric | published pages (figures/tables count toward the page total); LaTeX via `opticameeting`/`opticajnl` class on Overleaf | Optica author portal |
| Figures | vector PDF/EPS or ≥300–600 dpi raster; colorblind-safe encouraged | Optica style guide |

**Corrected assumption:** the task asked whether Optica's abstract limit is ≤35 words — it is **not**;
the Optica abstract target is **~100 words**. The ~10–12 pp main-text target in the task **exceeds** the
Optica **Research Article** guideline of **6–8 pp** (10–12 pp is mini-review length). This is a
**Research Article** (original results), so the recommendation is to target **~8 pp main** with the
supplement carrying the negatives and extra spectra; a 10–12 pp submission is possible but invites an
over-length flag. Logged as **GAP-4** in the matrix (length-budget decision for Phase 2, not an
evidence gap).

Sources:
- [Optica Publishing Group Journal Style Guide](https://opg.optica.org/content/author/portal/item/style-optica-styleguide/)
- [Author Resources: Journal Style Guide](https://opg.optica.org/submit/style/style_traditional_journals.cfm)
- [Author Resources (Optica review)](https://opg.optica.org/submit/review/optica.cfm)

---

## Title candidates (R45 §C4)
Preferred: **Engineering What an Optical Sensor Cannot See.** Alternatives: *Certified Blindness in
Computational Optics*; *Symmetry Walls and Certifiability in Optical Sensing*; *The Blindness Capacity
of an Optical Measurement*; *From Exact Optical Walls to Certifiable Change Detection*.

---

## Open items (for Phase 2 / user)
- **[GAP-4]** Length-budget decision: Optica Research Article 6–8 pp vs task's 10–12 pp target (recommend ~8 pp + supplement).
- **[GAP-1]** Hero cites the committed KT1 field wall **6.6e-16** (R45 §C3 rounds to 7e-16) — confirm caption wording.
- **[GAP-2]** Caption precision: joint-35 (spectrum) vs fresh-draw-80 (single plane) provenance kept distinct.
- **[GAP-3]** Two committed "→24" noncomposability instances (IT5 SVD→DPSS; IT4b Tukey) — do not conflate.
- **USER slots** (deferred, per program convention): author block, affiliation, repo URL wording, funding/URIS acknowledgment, ORCID.
- Theorem-proof transcription from R45 §B1/§B4/§B11 into main text + supplement (analytical; no data).
- Abstract draft (≤~100 words) + cover letter (Phase 2).
