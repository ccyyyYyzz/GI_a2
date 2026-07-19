# Optics Express submission checklist — both manuscripts

Covers the two R20/R21-overhauled manuscripts:

- **Paper 1** — `paper/main_oe.tex` — *When high flux helps single-pixel imaging:
  a contrast–dead-time operating map* (+ `paper/supplement.tex`).
- **Paper 2** — `paper2/main_m1.tex` — *Jitter-capped high-flux single-pixel
  imaging with global power control* (+ `paper2/supplement_m1.tex`).

Cover letters: `paper/COVER_LETTER_OE.md`, `paper2/COVER_LETTER_OE.md`.

Status legend: `[ ]` not done · `[~]` partially done, note in text · `[USER]`
requires a user decision.

---

## 1. opticajnl (Optica) template port

Both mains and both supplements are currently written against the standard
`article` class (`\documentclass[10pt,twocolumn]{article}`) and both mains
carry an in-file porting note stating the Optica port is a deferred task
(`paper/main_oe.tex:2`, `paper2/main_m1.tex:2–3`). The section structure of both
mains already follows the OE research-article shape, so the port is mechanical.

- [USER] **Download the opticajnl template.** Deferred pending user approval —
  do not fetch the template class/style files until the user authorizes the
  download. Target: current `opticajnl.cls` + `opticajnl.bst` from the Optica
  Publishing Group author-tools / Overleaf template.
- [ ] **Class swap (×4 files).** Replace
  `\documentclass[10pt,twocolumn]{article}` (paper1 main, paper2 main and their
  supplements) with the Optica journal class invocation (`\documentclass[...]{opticajnl}`
  with the OE options). Remove the manual `geometry` margin package and the
  manual two-column float tuning in `paper2/main_m1.tex:11–14`
  (`\dbltopfraction` etc.); the class governs these.
- [ ] **Bibliography style.** Both mains currently use `\bibliographystyle{unsrt}`
  (`paper/main_oe.tex:745`, `paper2/main_m1.tex:447`). Swap to the Optica style
  (`opticajnl`/`osajnl`) at port.
- [ ] **Backmatter blocks — reorder into Optica order.** Optica requires, in
  this order at the end of the main text:
  **Funding → Acknowledgments → Disclosures → Data availability → Supplemental
  document.** Current state:

  | Block | Paper 1 (`paper/main_oe.tex`) | Paper 2 (`paper2/main_m1.tex`) | Remaining |
  |---|---|---|---|
  | Funding | folded into Acknowledgments SPH (line 742) | folded into Acknowledgments SPH (line 445) | split into a dedicated `\section*{Funding}`; fill from [USER] |
  | Acknowledgments | `\section*{Acknowledgments}` line 741 (SPH) | `\section*{Acknowledgments}` line 444 (SPH) | fill from [USER]; move after Funding |
  | Disclosures | **missing** | **missing** | add `\section*{Disclosures}` (e.g. "The authors declare no conflicts of interest.") — [USER] confirm |
  | Data availability | `\section*{Code and data availability}` line 737 | `\section*{Code and data availability}` line 440 | rename to Optica "Data availability" (or keep code+data combined per OE policy); fill repo URL [USER] |
  | Supplemental document | referenced in-text as "supplement"/"Supplement" throughout | referenced in-text as "supplement" throughout | add explicit `\section*{Supplemental document}` pointing to "Supplement 1" (see §3) |
- [ ] **Optica-mandated ordering nuance:** Disclosures and Data availability must
  appear **before** the references; Supplemental document line goes with the
  backmatter. Verify against the current OE author guide at port time.

---

## 2. Figure requirements

All figures already exist as vector PDF plus PNG previews
(`paper/figs/*.pdf`, `paper2/figs/*.pdf`). Per-figure build records:
`paper/figs/FIG_MECHANISM_P1_README.md`, `paper2/figs/ACTIII_README.md`.

- [ ] **Vector PDF only** for the submitted build (drop the PNG twins from the
  LaTeX include path; keep them only as previews).
- [ ] **Fonts fully embedded** in every PDF — verify with
  `pdffonts figs/*.pdf` (every font must show `emb yes`).
- [ ] **Minimum text size 7.5 pt at final printed width.** Optica single-column
  ≈ 3.25 in (≈ 8.2 cm), double-column ≈ 6.75 in (≈ 17 cm). Several current
  figures were sized against the `article` two-column geometry; after the class
  swap changes the column width, **re-run a resize/regeneration pass** and
  re-check the smallest tick/label against 7.5 pt at the new final width.
- [ ] **Figure inventory** (confirm each survives the R20 rebuild and the resize
  pass):
  - Paper 1: `fig_mechanism_p1`, `fig_b_masks`, `fig_c_ladder`, `fig_d_dwell`,
    `fig_e_gate`, `fig_f_scattering`, `fig_p1_families`, `fig_nat_pair`,
    `fig_nat_pair_16`, and supplement figs `fig_s1b_images`, `fig_s1c_curves`,
    `fig_s1d_svsrho`, `fig_s2_gallery`.
  - Paper 2: `fig_mechanism`, `fig_speed_results`, `fig_jitter_validation`,
    `actiii_a…e`, `fig_audit_supp`.
- [ ] **R20/R21 pixel-level check** (from `docs/ROUND63_GPT_ROUND20_REVIEW_RAW.md`
  §Scope note): the last 100%-scale print inspection of hierarchy, captions, and
  typography must be repeated **after** the template port, on the ported build.

---

## 3. Supplemental document

- [ ] Compile each supplement standalone as the Optica **"Supplement 1"** PDF
  (`paper/supplement.tex` → Supplement 1 for paper 1;
  `paper2/supplement_m1.tex` → Supplement 1 for paper 2). Each paper has its own
  "Supplement 1".
- [ ] Add the Optica "Supplemental document" backmatter line in each main text
  referencing "Supplement 1" (see §1).
- [ ] Optica requires supplementary material to be **cited in the main text**;
  both mains already reference the supplement extensively, so ensure at least
  one explicit "Supplement 1" citation survives the wording port.
- [~] Supplement author blocks must mirror the main author block once filled —
  `paper2/supplement_m1.tex:35` is an SPH mirror; `paper/supplement.tex:31` is
  an **empty** `\author{}` that must be synced to the paper-1 author block.

---

## 4. arXiv package

- [ ] Assemble a self-contained source tarball per paper: `main_*.tex`,
  `supplement*.tex`, `refs.bib`, the compiled `.bbl`, and the `figs/*.pdf`
  actually included (exclude PNG previews and build junk — see §7).
- [ ] Verify the tarball compiles from clean (`latexmk -pdf`) with no missing
  assets and **zero surviving `\SPH{}` placeholders** (campaign rule: zero red
  placeholders at submission).
- [USER] **arXiv license choice** — select the arXiv distribution license
  (e.g. CC BY 4.0 vs arXiv non-exclusive). No license statement currently exists
  in either source; decide before upload and keep it consistent with the Optica
  copyright/open-access choice.
- [ ] Post the arXiv preprint and update the companion `@misc{companion}` entry
  (`paper2/refs.bib:388`) with the paper-1 arXiv ID; add a matching companion
  entry in `paper/refs.bib` pointing at paper 2 (see reported inconsistency).

---

## 5. USER slots to fill (file:line)

The seven author/URL/funding decision points, each a red `\SPH{}`:

1. `paper/main_oe.tex:31` — Paper 1 author block (affiliation + collaborators).
2. `paper/main_oe.tex:739` — Paper 1 repo URL wording (Data availability).
3. `paper/main_oe.tex:742` — Paper 1 Funding / Acknowledgments text.
4. `paper2/main_m1.tex:23` — Paper 2 author block.
5. `paper2/main_m1.tex:442` — Paper 2 repo URL wording (Data availability).
6. `paper2/main_m1.tex:445` — Paper 2 Funding / Acknowledgments text.
7. `paper2/supplement_m1.tex:35` — Paper 2 supplement author block (mirror of #4).

Related sync/verification items (not among the seven but must not ship unfilled):

- `paper/supplement.tex:31` — empty `\author{}`, sync to #1.
- `paper/main_oe.tex:245` — §X.5 relation-to-prior-work SPH pending the Grönberg
  full-text sign-off (S5 hard reference item); resolve before submission.
- Suggested-reviewer lists in both cover letters (`[USER: suggested reviewers]`).

---

## 6. ORCID and funding-ID reminders

- [USER] Register/confirm an **ORCID iD for every author** and enter them in the
  Optica submission system (Prism) — ORCID is collected at submission, not in the
  `.tex`.
- [USER] Supply the **funder name(s) and grant/award numbers** for the Funding
  block (§1). Optica cross-checks these against the Crossref Open Funder Registry
  and the Funding statement must match the metadata entered in Prism exactly.
- [USER] Confirm the **corresponding author** and their institutional email.

---

## 7. Repository hygiene before public release

The repo currently carries LaTeX build junk alongside source. Sweep before the
repository is made public / the arXiv tarball is cut:

- [ ] Remove or `.gitignore` build artifacts under `paper/` and `paper2/`:
  `*.aux`, `*.bbl`, `*.blg`, `*.log`, `*.out`, `*.fls`, `*.fdb_latexmk`,
  `main_*.txt`, `supplement*.txt`, and `tmp_compile_editorial_*.log`
  (e.g. `paper/tmp_compile_editorial_main_oe.log`).
- [ ] Confirm the two stray untracked paths in `git status`
  (`E：/` directory and `water_lens_repro/`) are intended before any public
  push; exclude if not part of the release.
- [ ] Verify the immutable freeze tags referenced by the manuscripts still resolve
  (`f1-freeze`, `study2-freeze`, `study2-scored`, `m1-freeze`,
  `paper1-freeze-v1`) and that the Data availability URL points at the intended
  public mirror.
- [ ] Final zero-placeholder scan across both mains and supplements:
  `grep -rn "\\SPH" paper/ paper2/` must return only the `\newcommand` definition
  lines.

---

## Cross-paper consistency note (found during drafting)

- The companion cross-reference is currently **one-directional**: paper 2 cites
  paper 1 via `\cite{companion}` with a `@misc{companion}` entry
  (`paper2/refs.bib:388`), but **paper 1 has no reciprocal bib entry or citation
  to paper 2** — it only alludes to jitter qualitatively in its Discussion
  (`paper/main_oe.tex:695`). Add a companion entry + citation in paper 1 (or
  decide deliberately to keep it one-way) before submission so the two papers
  cross-reference symmetrically.
- The deterministic optimum at ν = 2000 is written as ≈22 in one place and
  ≈22.3 in another within paper 2 (`main_m1.tex:111` vs `:78`); both round the
  same value (6·2000)^(1/3) − 2/3 = 22.23 and agree with paper 1's ridge law, so
  this is cosmetic rounding, not a claim conflict — optionally unify the wording.
