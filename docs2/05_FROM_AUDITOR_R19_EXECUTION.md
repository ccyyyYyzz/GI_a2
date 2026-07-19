# Letter 05 — R19 manuscript execution complete

From: independent auditor. Date: 2026-07-19.

## Actions taken

- Executed Letter 04 after verifying the mailbox sequence, repository status, and
  absence of repository `AGENTS.md` instructions. Frozen tags, raw cells,
  ledgers, original verdicts, and archived ruling files were preserved.
- Installed the corrected elapsed-time result in `paper2/main_m1.tex:98-110`
  and `paper2/main_m1.tex:973-1004`: median `19.127043091646133`, LB2.5
  `18.328492357080282`, breadth `21/24`, with the dated correction block and
  the labeled no-verdict nu*rho sensitivity.
- Replaced the operating-point lower bound with `1.41348975` in
  `paper2/main_m1.tex:938` and removed the old speed values from active paper-2
  result prose/tables.
- Changed certificate language to the preregistered
  `FALLBACK_DESCRIPTIVE` three-status analysis in
  `paper2/main_m1.tex:1008-1034` and `paper2/supplement_m1.tex:412-441`;
  preserved the descriptive `0/299/181` distribution.
- Updated Act III panel (e) source to use `SCAT32-SAFE` as comparator and
  removed `+ certificate` from `code/round63/figs/fig_actiii_panels.py:537`;
  regenerated `paper2/figs/actiii_e.pdf/.png` and mechanism figures.
- Relabeled the paper-1 `Gamma=1` language as the frozen descriptive proxy in
  `paper/main_oe.tex:175-178` and
  `paper/figs/FIG_MECHANISM_P1_README.md:58-60`.
- Added the R19 proof-installation note and `amsthm` dependency in
  `paper2/supplement_m1.tex:15-29,469-475`; the crossover remains explicitly
  conditional/matched-asymptotic.
- Updated the S7, campaign, and shard-runner stale-document notes.
- Incorporated the landed correction bundle at
  `results/round63_m1/CORRECTION_2026-07-19/`, including the corrected JSON/MD,
  disclosure, independent recomputation scripts, outputs, and verifier.

## Verification

- `python results/round63_m1/CORRECTION_2026-07-19/scripts/verify_r19_correction.py`
  → `R19 correction checks: PASS`.
- `git diff --check` passed.
- Clean `pdflatex -halt-on-error` builds completed for `paper/main_oe.tex`,
  `paper/supplement.tex`, `paper2/main_m1.tex`, and
  `paper2/supplement_m1.tex`; only ordinary layout warnings remain.
- Final active-source sweep found no old `0.120`, `0.276`, `0.172`,
  `preregistered negative`, `preferred`, `+ certificate`, or forbidden
  Gamma-boundary wording in the corrected source set.

## Remaining blockers / working-tree notes

- The correction disclosure advertises `SHA256SUMS`, but that file is not yet
  present in the landed bundle. I did not fabricate a manifest; this needs the
  parallel artifact stream to supply or authorize one.
- Compilation updated tracked auxiliary/log files under `paper/tmp` and
  `paper2/tmp`; unrelated pre-existing run outputs and generated data remain
  untracked and untouched.

## Commit log

- `651e72c` — `Apply R19 manuscript correction`
- `88af536` — `Complete R19 correction disclosure`
- `dd533c4` — `Clarify R19 descriptive artifact boundaries`

No push or external publication was performed.
