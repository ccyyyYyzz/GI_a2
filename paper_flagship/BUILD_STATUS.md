# BUILD STATUS — flagship manuscript skeleton

**Repo:** `D:\GI_another` · **Dir:** `paper_flagship/` · **Compiled:** 2026-07-22
**Architecture:** R35 frozen (LIMIT → REACH → CLOSE → REUSE). No commits.

Fill-state legend: **skeleton** = structure + TODO pointers only · **drafted** =
prose installed, authorized by R35 · **frozen** = verbatim custody-controlled
string installed (see `FROZEN_REGISTER.md`).

---

## 1. Compile status

| Artifact | pdflatex | Pages | Errors | Notes |
|---|---|---:|---:|---|
| `main.tex` | ✔ 2 passes | **5** | **0** | 5 overfull hboxes (placeholder boxes + italic transition quotes; cosmetic, expected in skeleton) |
| `NOTATION.tex` | ✔ 2 passes | 3 | 0 | standalone; longtable notation reconciliation |
| `supplement/supplement.tex` | ✔ 2 passes | 1 | 0 | headings only; content is TODO |

- Toolchain: TeX Live 2024 `pdflatex`. **No `opticajnl.cls` installed / in-repo**
  → skeleton uses the same `article` two-column preamble as `paper/` + `paper2/`.
  **Optica port is gap G6** (PORT-TODO in `main.tex` header).
- Skeleton is **bibtex-free**: citations carried as `\CITE{}` markers, not
  `\cite`, until the merged bibliography (gap G5) is built. `\bibliography` is
  commented out.
- 5 composed skeleton pages ≪ the 9.0-page hard cap; page budget below tracks the
  target once prose lands.

---

## 2. Per-section fill state (main.tex)

| § | Printed section | R35 page budget | Fill state | What is installed now |
|---|---|---:|---|---|
| Title | — | — | **frozen** (pass) + comments | Pass title #1 active; alternates #2–#6 + fail #7/#8 + forbidden-word list as comments |
| Abstract | — | ~100 words | **skeleton** | `\SPH` body wired to the 6-beat sheet; §1.1 thesis (pass) + §1.2 (fail) in comment blocks |
| 1 | One hidden channel, two consequences | 1.25 | **drafted** (partial) | Frozen opening + 2nd sentence installed; frozen Intro→§2 transition; paras 1–3 finish = TODO |
| 2 | Hidden-state information sets a finite ridge | 1.45 | **skeleton + frozen transitions** | 4 equation-import TODOs (paper2 Thm 1–3 + paper-1 ridge + crossover); 2 frozen transitions installed |
| 3 | One global control reaches the ridge | 1.45 | **frozen results installed** | R19 fixed-dwell + resource + speed + frame sentences **verbatim**; beats 1/4/5 = TODO; frozen end-of-M1 + §3→§4 transitions |
| 4 | Closing the integrated-count scene channel | 0.95 | **drafted** (3 paras) | All three closure paragraphs drafted from frozen certs (+0.68/−8.3/+0.04/≤0.2/+0.41); frozen scope + pivot verbatim |
| 5 | One record, two scientific outputs | 1.75 | **skeleton (PASS-BRANCH)** | Block-Fisher Eq. installed; DCS 3-sentence priority **drafted**; R34 claim + 7 confirmatory beats = `\CAMPAIGN{}` TODO |
| 6 | Boundary, scope, experimental transfer | 0.60 | **skeleton** | MOLT one-sentence slot + bench-falsification slot + 3-sentence conclusion = TODO; §5→Disc transition in comment |
| M | Methods and provenance capsule | 0.40 | **skeleton** | 6-point stub with materials-map pointers |
| — | References / back matter | 1.15 | **blocked (G5)** | merged bib not built |
| **Σ** | | **9.00** | | skeleton renders 5 pp |

Branch discipline: every PASS-BRANCH element (title, abstract, §5, Fig. 1 Stage 3,
Fig. 5) tagged `%% PASS-BRANCH` with an adjacent `%%>> FAIL-BRANCH` block →
degradation is a mechanical edit (R35 §8).

---

## 3. Figures

| Fig | Section | State | Source now | Blocking |
|---|---|---|---|---|
| 1 hero | §1 | placeholder box | 3 per-paper heroes exist (paper1 mech, paper2 jitter, DLGI 3-panel) | **G3** unified 3-stage arc unbuilt; Stage 3 gated on **G2** |
| 2 ridge | §2 | placeholder box | `paper2/figs/fig_mechanism.pdf`, `fig_jitter_validation.pdf`, `fig_a_ridge_map.pdf` | panel (a) schematic new |
| 3 imaging | §3 | placeholder box | `paper2/figs/fig_speed_results.pdf` (near-direct reuse) | relabel internal tokens |
| 4 closure | §4 | placeholder box | 5 frozen cert numbers; no figure exists | **G4** double-sided schematic unbuilt |
| 5 dual (PASS) | §5 | placeholder box | `DUAL_LEDGER_PROBE/figs/fig_dual_ledger.png` (feasibility only) | **G2** confirmatory grid; panels b/c/d `\CAMPAIGN{}`-gated |

Full specs: `figures/FIGURE_SPECS.md`. **Five figures, no sixth** (hard rule).

---

## 4. Supplement (S1–S10, R35 §9)

All ten stubs present with content manifests + materials-map pointers; master
`supplement/supplement.tex` compiles.

| Note | Title | State | Port source |
|---|---|---|---|
| S1 | Unified notation + hidden-state experiment | stub | `NOTATION.tex` (G1) + paper2 S7 |
| S2 | Detector hidden-state specialization | stub | paper2 Thm1–2/S7 + paper-1 Eq.(missing) + paper-1 S1 |
| S3 | Ridge asymptotics + verification | stub | paper2 S4 + main_m1 §3 + paper-1 §4.2 (+G10 flag) |
| S4 | M1 protocol + provenance | stub | paper2 §5 + S5 + `CORRECTION_2026-07-19/` (frozen R19) |
| S5 | Complete M1 results + operating maps | stub | paper2 S6 + main_m1 §6 + paper-1 supp_tables |
| S6 | Count-only closure campaigns | stub | 5 cert reports/JSON (acronyms live HERE) |
| S7 | Reciprocity + DLGI estimator | stub | R33 Thm1–2 + DUAL_LEDGER_PROBE + LADDER_PROBE |
| S8 | Interval construction + C1–C7 campaign | stub | **BLOCKED G2** — R34 spec only |
| S9 | Detector-side escape (MOLT) | stub | R31/R31-PRO + JAILBREAK_PROBE + REVIVAL_CALC |
| S10 | Reproducibility | stub | provenance anchors + round ledger + env |

Note: R35 §9 enumerates S1–**S10** (the task brief said "S1–S9"; S10
Reproducibility is included to match the authoritative R35 list).

---

## 5. Dependency ledger — what waits on C1–C7 (gap G2)

**The R34 confirmatory campaign does not yet exist.** Until it runs and all seven
bars clear, the entire REUSE movement is contingent (R35 §5 opening).

| Depends on C1–C7 | Manifestation in skeleton | If PASS | If FAIL (R35 §8) |
|---|---|---|---|
| §5 entire section | `%% PASS-BRANCH` + all numbers `\CAMPAIGN{}` | keep §5; restore "two model-certified products" | delete §5 → 0.45-pp Discussion subsection "What lies beyond the count-only boundary" |
| Abstract beats 5–6 | `\SPH` body notes | keep pass abstract | switch to fail thesis (A2) |
| Title | pass #1 active | keep | swap to #7/#8 |
| Fig. 1 Stage 3 | `%% PASS-BRANCH` in hero box | keep reuse stage | remove Stage 3 (hero ends at boundary) |
| Fig. 5 | placeholder (PASS) | build from campaign grid | drop figure → 4 figures total |
| Page budget | 9.0 pp / 5 figs | as-is | 7.5–8.0 pp / 4 figs |

C1–C7 bars (kill bars, materials map §4d): C1 coverage ∈[0.92,0.98] every cell ·
C2 width ≤1.5× pilot · C3 RMSE ≤1.5× pilot · C4 scene noninferiority ≤0.2 dB &
≤5% Fisher · C5 model/gauge · C6 reciprocity+schedule · C7 identical-ledger +
edge honesty. Grid: t_c{16,32,64}×CV{0.05,0.15,0.40}×photon{0.5,1,2}× = 27 cells;
≥5000 calib + ≥1000 confirm/cell; ≥12 fresh scenes; frozen Neyman; repair once.

**Probe caveat:** the DLGI *feasibility* probe already FAILED bar 4 (interval
coverage) — 6 PASS / 1 FAIL (materials map §4a). The frozen-Neyman campaign is
the confirmatory retest; a second failure sends the manuscript to the fail branch.

---

## 6. Other open gaps (from materials map §10) blocking the writing phase

| Gap | Blocks | State |
|---|---|---|
| G1 unified notation | writing | **addressed** — `NOTATION.tex` built (collisions flagged; freeze pending) |
| G2 DLGI C1–C7 campaign | §5 / Fig. 5 / abstract | **open** — not run |
| G3 flagship hero figure | Fig. 1 | open — placeholder only |
| G4 Act-III closure figure | Fig. 4 | open — placeholder only |
| G5 merged bibliography | references / all `\CITE` | open — 2 bibs + ~39 DOIs to dedupe; `companion` key to drop; Grönberg sign-off |
| G6 Optica template port | pre-submission | open — `article` class placeholder + PORT-TODO |
| G7 author/repo/funding | front matter | open — `\SPH` sites |
| G8 Act-I/II unification prose | §2 | open — crossover Eq. is the bridge; TODO installed |
| G9 cut-scope adjudication | §4/materials bank | open — paper-1 Study-1/2 CUT decision |
| G10 R16 2M-frame diagnostic | S3 flag | open — carry-forward |

---

## 7. Next actions (writing phase, not this skeleton)

1. Run the C1–C7 campaign (G2) → unblocks §5, Fig. 5, abstract, branch decision.
2. Draft §1 paras 1–3, §2 equations (import + unify, G8), §3 beats 1/4/5, §5 beats.
3. Build the five figures (G3/G4 new; Figs 2/3/5 adapt existing generators).
4. Merge + dedupe bibliography (G5); clear Grönberg sign-off.
5. Fill `\SPH` front-matter (G7); Optica port (G6) at submission prep.
