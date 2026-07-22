# BUILD STATUS — flagship manuscript skeleton (R36)

**Repo:** `D:\GI_another` · **Dir:** `paper_flagship/` · **Compiled:** 2026-07-22
**Architecture:** **R36 frozen** (DEFINE → DECODE → CERTIFY, one DLGI
protagonist). R36 (issue #28, `79a4508`) **overrides R35** (issue #27). No commits.

Fill-state legend: **skeleton** = structure + TODO pointers only · **drafted** =
prose installed, authorized by the ruling · **frozen** = verbatim custody-
controlled string installed (see `FROZEN_REGISTER.md`).

---

## 0. What changed from the R35 build (restructure log)

The operator rejected the R35 four-movement version as bloated. R36 rebuilt it:
- **Spine:** LIMIT→REACH→CLOSE→REUSE (5 sections + Methods, 9.0 pp, 5 figs) →
  **DEFINE→DECODE→CERTIFY (5 sections + Methods, 6.80 pp, 4 figs)**.
- **Protagonist:** the whole program → **DLGI + reciprocity theorem only**.
- **M1/ridge (+1.87 dB / 19.13×):** main-text Results → **stripped from the body
  (binding cut)**; frozen sentences **moved verbatim to Supplement S4**; one
  restoration line left in `main.tex` §2 tagged `%% EDITOR-REQUEST-ONLY`.
- **Five closure campaigns:** a full section + figure → **one frozen four-sentence
  boundary paragraph in the Introduction** (the only main-text trace) + Supp S6.
- **Fail branch:** inline `%% PASS/FAIL` tags → **separate file
  `main_failbranch_notes.tex`** (protagonist switch makes the branches
  structurally different; inline tagging would defeat the slimming).
- **Supplement:** renumbered R35 S1–S10 → **R36 S1–S10** (crosswalk §4 below).
- **Figures:** R35 five-figure specs → **R36 four-figure architecture** (old specs
  archived at the bottom of `figures/FIGURE_SPECS.md`).

---

## 1. Compile status

| Artifact | pdflatex | Pages | Errors | Notes |
|---|---|---:|---:|---|
| `main.tex` | ✔ 2 passes | **4** | **0** | slimmer than the R35 skeleton (was 5); overfull boxes cosmetic |
| `NOTATION.tex` | ✔ 2 passes | 3 | 0 | unchanged; now front matter of the supplement (R36 §6 has no notation note) |
| `supplement/supplement.tex` | ✔ 2 passes | 1 | 0 | headings only; content TODO. (One transient error — a literal `\input` inside prose — fixed.) |

- Toolchain TeX Live 2024 `pdflatex`; **no `opticajnl.cls`** → `article` two-column
  preamble + Optica PORT-TODO (**gap G6**).
- Bibtex-free skeleton: `\CITE{}` markers, `\bibliography` commented (**gap G5**).
- 4 skeleton pages ≪ the 7.0-page hard cap; the 6.80-page target is tracked below.

---

## 2. Per-section fill state (main.tex) — R36 §1 page-and-claim map

| § | Printed section (R36) | Budget | Fill state | Installed now |
|---|---|---:|---|---|
| Abstract | — | 80–100 w | **skeleton** | `\SPH` body wired to the DEFINE/DECODE/CERTIFY beats; binding-cut + pre-verdict wording noted |
| 1 | One record, two products | 0.80 | **drafted + frozen** | R36.2 opening installed; **R36.3 frozen boundary paragraph verbatim**; **R36.1 manuscript-contract sentence verbatim**; frozen Intro→theory transition; para 1/3 finish = TODO |
| 2 | Two ledgers of one hidden-state experiment (DEFINE) | 1.05 | **equations installed + frozen transition** | Eqs (1)–(4) installed (hierarchy, partition, Schur, K+reciprocity); Louis credit + O4-A = TODO; dead-time sibling cut + `%% EDITOR-REQUEST-ONLY` line; frozen theory→DLGI transition |
| 3 | Pilot-free scene-and-medium inference (DECODE) | 0.75 | **equation installed + frozen transition** | Eq (5) obs+gauge installed; residual/inference/arms = TODO; word-cap routing noted; frozen architecture→evidence transition |
| 4 | Two products under one acquisition ledger (CERTIFY) | 2.10 | **skeleton** | R34 authorized claim (frozen) referenced; 3-claim structure; all numbers `\CAMPAIGN{}`; frozen evidence→scope transition |
| 5 | Scope and experimental transfer | 0.55 | **skeleton** | 3 one-paragraph TODOs (physical scope / prior-art boundary / bench transfer) |
| M | Methods and provenance capsule | 0.40 | **skeleton** | 6-point stub; word-cap guard installed |
| — | References / back matter | 1.15 | **blocked (G5)** | merged bib not built |
| **Σ** | | **6.80** | | skeleton renders 4 pp |

Equation budget: **5 displayed equations used / 5 max** (Eqs 1–4 in §2, Eq 5 in
§3). Word caps: Kalman/RTS/Whittle/bootstrap each appear **exactly once on a
printed line** — all inside the single §3 word-cap TODO (0 in final manuscript
prose, which is still to be written); Methods is guarded against repeating them.
Verified `grep -vE '^\s*%'` → 1 each. Satisfies R36 §3.

---

## 3. Figures (R36 four-figure architecture)

| Fig | Section | State | Source now | Blocking |
|---|---|---|---|---|
| 1 hero | §1 | placeholder box | `paper/figs/fig_mechanism_p1.pdf` (left) + `DUAL_LEDGER_PROBE` (right) | **G3**; "model-certified" wording gated on **G2**; **binding cut: no M1 numbers** |
| 2 reciprocity | §2 | placeholder box | R33 Thm 2 + `DUAL_LEDGER_PROBE/t3_reciprocity.json` | schematic unbuilt; measured schedule ordering **G2**-gated |
| 3 precision+calibration | §4 | placeholder box | `DUAL_LEDGER_PROBE/figs/fig_dual_ledger.png` panel A (feasibility) | **G2** confirmatory grid; all `\CAMPAIGN{}` |
| 4 scene+domain | §4 | placeholder box | `DUAL_LEDGER_PROBE/` noninferiority + edge (feasibility) | **G2** grid + edge bank; all `\CAMPAIGN{}` |

Full specs + R35→R36 mapping + archive: `figures/FIGURE_SPECS.md`. **Four figures,
no fifth.** §3 and §5 carry no new visual.

---

## 4. Supplement crosswalk (R35 S1–S10 → R36 S1–S10)

R36 §6 mapping. Files renamed; master `supplement/supplement.tex` updated; all
content manifests rehomed. Notation table (gap G1) is now supplement front matter
(R36 §6 defines no notation note).

| R36 note | File | Holds | From R35 |
|---|---|---|---|
| S1 | `S1_regularity_partition.tex` | regularity + joint missing-info proof | old S1(notation)+S2(detector) proof half |
| S2 | `S2_reciprocity_proof.tex` | reciprocity proof, pseudoinverse, determinants | old S7 (proof half) |
| S3 | `S3_dlgi_model_estimator.tex` | DLGI model, gauge, identifiability, estimator, calibration | old S7 (estimator half), expanded |
| S4 | `S4_deadtime_sibling.tex` | **dead-time sibling: ridge law + M1 +1.87 dB/19.13× + maps** (binding-cut destination; **frozen R19 sentences verbatim**) | old S2+S3+S4+S5 merged |
| S5 | `S5_auxiliary_derivations.tex` | auxiliary hidden-state + schedule derivations (**+ folded MOLT**) | old S9 (MOLT) folded in |
| S6 | `S6_closure_certificates.tex` | five closure certificates at full depth | old S6 (role unchanged) |
| S7 | `S7_priorart_fence.tex` | prior-art fence (DCS/DWS, pilot, self-cal, AO/radio) | NEW dedicated note |
| S8 | `S8_confirmatory_protocol.tex` | R34 C1–C7, full grid, calibration + multiplicity | old S8 (role unchanged) |
| S9 | `S9_edge_sensitivity.tex` | edge maps, residual/model checks, sensitivity | NEW (old S9 slot was MOLT → moved to S5) |
| S10 | `S10_provenance.tex` | provenance, correction record, code, seeds, repro | old S10 (role unchanged) |

Note on MOLT: R36 §6 gives MOLT no dedicated note but the executive keeps it "in
S1–S10", so it is folded into **S5** as a labeled auxiliary subsection.
Acronyms/audit labels appear only in **S6** and **S8** (R36 §2, §6).

---

## 5. Fail branch (R36 §7) — separate variant, not inline tags

`main_failbranch_notes.tex` (not compiled, not `\input` by `main.tex`) holds the
protagonist-switch spec: fail thesis (R36.6), title #7/#8, **≤6.2 pp / ≤3
figures**, section spine restoring ridge + M1 to the main text, the three
**RETIRED** R35 transitions reinstated (legitimate here), and DLGI demoted to the
supplement with its failed endpoint disclosed. Rationale for the file-not-tags
mechanism: under R36 the two branches differ *structurally* (different thesis,
order, figures, cap), so inline tagging would fold two papers into one file and
defeat the slimming. Promote to `main_law.tex` on any failed C1–C7 bar.

---

## 6. Dependency ledger — what waits on C1–C7 (gap G2)

**The R34 confirmatory campaign does not yet exist.** The entire CERTIFY movement
(§4) is contingent (R36 §4; pass section appears in full only after all bars pass).

| Depends on C1–C7 | Manifestation | If PASS | If FAIL |
|---|---|---|---|
| §4 entire section | 3-claim skeleton, all numbers `\CAMPAIGN{}` | open with one-line headline; restore "model-certified" | switch to `main_failbranch_notes.tex` (protagonist switch) |
| Abstract | `\SPH` body, "two inferred products" | keep; upgrade wording | fail thesis (R36.6) |
| Fig. 1 wording | "two inferred products" | "two model-certified products" | recast to law/boundary hero |
| Figs. 3 & 4 | placeholders, `\CAMPAIGN{}` | build from S8 grid | replaced by fail-branch fig set |
| Title | "One bucket record measures both scene and medium" | keep | swap to #7/#8 |

C1–C7 bars (materials map §4d): C1 coverage ∈[0.92,0.98]/cell · C2 width ≤1.5×
pilot · C3 RMSE ≤1.5× pilot · C4 scene noninferiority ≤0.2 dB & ≤5% Fisher · C5
model/gauge · C6 reciprocity+schedule · C7 identical-ledger + edge honesty. Grid
27 cells; ≥5000 calib + ≥1000 confirm/cell; ≥12 fresh scenes; frozen Neyman;
repair once. **Probe caveat:** the feasibility probe already FAILED bar 4
(coverage) — 6 PASS / 1 FAIL (materials map §4a).

---

## 7. Gap ledger (materials map §10 + R36 deltas)

| Gap | Blocks | State under R36 |
|---|---|---|
| G1 unified notation | writing | **addressed** — `NOTATION.tex` built; now supplement front matter; collisions still to freeze |
| G2 DLGI C1–C7 campaign | §4 / Figs 3–4 / abstract / branch | **open** — not run |
| G3 flagship hero figure | Fig. 1 | open — placeholder (recut to hidden-channel-as-instrument, no M1 numbers) |
| G4 closure figure | — | **downgraded** — closure is now text (R36.3) + Supp S6; a closure figure is only needed in the fail branch |
| G5 merged bibliography | references / `\CITE` | open — 2 bibs + ~39 DOIs; drop `companion`; Grönberg sign-off |
| G6 Optica template port | pre-submission | open — `article` + PORT-TODO |
| G7 author/repo/funding | front matter | open — `\SPH` sites |
| **G8 Act-I/II unification prose** | — | **RESOLVED / MOOT** — see §8 |
| G9 cut-scope adjudication | supplement | **largely resolved** — R36 demotes M1/ridge/closure to supplement wholesale; the cut is now architectural, not open |
| G10 R16 2M-frame diagnostic | S4 flag | open — carry-forward into S4 (dead-time sibling) |

---

## 8. Risk list (updated for R36)

**G8 (Act-I/II seam) is GONE — confirmed.** R36 removes the ridge law and jitter
cap from the main text entirely (binding cut → Supplement S4). There is no longer
an Act-I/II progression in the body to unify: §2 is a single hidden-state
partition + reciprocity theorem, with the dead-time material demoted wholesale.
The two-theorem-set seam that G8 flagged cannot appear in a paper that carries
neither theorem set in its body.

Top risks now:

1. **G2 — CERTIFY rests on an unrun campaign with a failing precedent.** Under
   R36 this is *more* concentrated than before: §4 is the paper's payload and the
   feasibility probe already failed the coverage bar (bar 4). A second failure
   triggers the full protagonist switch to the fail branch.
2. **The "just a Kalman estimator on its own OU simulation" attack (R36 §5).** R36
   names this the top structural risk after slimming. The architecture answers it
   by ordering (principle → theorem → demoted estimator → external baselines →
   falsifiability), so the *writing* must hold that order exactly; any leak of
   filter detail into §3, or of a self-referential result into Figs 3–4, reopens it.
3. **G1 — frozen-glyph collisions.** Still stands and is now sharper: `θ`
   (log-rate vs medium hyperparameters), `c_v` vs `CV`, `J` (scalar) vs `J_x,J_θ`
   (ledgers), `K` vs `k`, `A,B` blocks. R19/R33/R34 strings are custody-frozen
   using these glyphs; the dead-time `c_v/ρ` glyphs now live only in S4, which
   slightly *reduces* main-text collision surface but the S4↔body consistency
   must still be policed.
4. **G5 — bibliography integrity.** Unchanged: merge two bibs + ~39 DOIs, drop the
   `companion` self-citation key, clear the Grönberg full-text sign-off.
5. **6.80-page cap with 0.2 pp slack (R36 §1).** §4 alone is 2.10 pp carrying
   three claims + two figures; the boundary paragraph, reciprocity, and DCS
   sentence all compete for the 0.80 pp Introduction. Very little room before a
   section overruns its budget.

---

## 9. Next actions (writing phase)

1. Run the C1–C7 campaign (G2) → unblocks §4, Figs 3–4, abstract, branch decision.
2. Draft §1 paras 1/3, §2 equation readings + Louis credit + O4-A, §3 architecture
   prose (spend the four word-caps once), §4 three claims, §5 three paragraphs.
3. Build four figures (G3 hero new; Fig 2 schematic new; Figs 3–4 from S8 grid).
4. Merge + dedupe bibliography (G5); clear Grönberg sign-off.
5. Fill `\SPH` front matter (G7); Optica port (G6) at submission prep.
6. Freeze the notation collisions (G1) before first full draft.
