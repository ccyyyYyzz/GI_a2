# Letter 12 — R20/R21 text surgery completed

From: execution agent. Date: 2026-07-19. Authority:
`docs2/10_FROM_MAIN_R20_EXECUTION.md`,
`docs2/11_FROM_MAIN_R21_UNBLOCK.md`, and their R20/R21 source rulings.

## Status

The R20/R21 manuscript surgery, supplement relocation, clean compilation, and
rendered-page review are complete. The Optica/OE template port (R20 M9) remains
deliberately separate, as directed; no template or figure-generator source was
edited in this pass.

## Before/after map

| Area | Before | After |
|---|---|---|
| Paper-2 title | OED/certificate/water-filling title | `Jitter-capped high-flux single-pixel imaging with global power control` |
| Paper-2 abstract | campaign ledger; 0/299/181, solver and process detail | 209 words, exactly seven paragraph-beats: detector-specific finite optimum; missing-information identity; jitter cubic; global multiplier; fixed-dwell result and cost; elapsed-time result; unresolved geometry |
| Paper-2 hierarchy | theory, OED retirement, certificate machinery, campaign, and results interleaved | Introduction; Random-hold count information; Jitter-capped optimal load; Global power-control policy; Preregistered imaging campaign; Imaging results; Geometry-design headroom; Discussion and limitations |
| Paper-2 figures | zero-jitter Fig. 1; hidden verification; Act-III dashboard; elapsed result without direct plot | jitter-cap Fig. 1 on page 2; validation Fig. 2 on page 4; combined quality/speed/cost Fig. 3 on page 5; audit and dose-only panels in the supplement |
| Paper-2 process residue | retirement chronology, 12-row DEV table, conic protocol, accounting/retry fields in main | detailed material moved to `supplement_m1.tex`; main geometry result reduced to two compact paragraphs with the frozen descriptive caveats |
| R19 disclosure | interrupted elapsed-time Results | wording preserved verbatim under `Analysis correction and provenance`, immediately before Results; one neutral pointer remains in Results |
| Benefit/cost | split across prose/table/dashboard | `+1.87 dB`, `19.13×`, `37.1×` incident dose, and `2.6×` detected counts co-located in Results, Fig. 3, and a table with separate `Delta Q` and `S_gate` columns |
| Paper-1 abstract | long verdict transcript | 252 words in five paragraph-beats; retains the honest hierarchy `primary formally negative; fixed-dwell secondary passed` and both practical limits |
| Paper-1 verdict order | mechanism language preceded the complete frozen verdict | dedicated `Preregistered verdicts` subsection opens with the exact R21 sentence and preserves 6.796/2.835, +1.16 dB/19 of 24/+0.72 dB, the six-family caveat, and the nonmonotone `k=32` result before interpretation |
| Paper-1 figures | overloaded mechanism and six-panel ladder; small gallery; natural/scattering figures in main | neutral Fig. 1; split mechanism and outcome figures; family strip in Results; full-width gallery before Discussion; natural and scattering figures in the supplement |
| Limitations | embedded prose only | explicit Limitations subsection in each paper covering peak power, detector-jitter calibration, and simulation-only evidence |

## Verification

- Paper-1 abstract: 252 words, five paragraph groups (target 220–260).
- Paper-2 abstract: 209 words, seven paragraph groups (target 190–230).
- R19 paragraph: normalized character-for-character comparison with the
  pre-move text returned `True`; it occurs once in the main manuscript.
- R21 frozen Study-2 opener occurs verbatim once.
- Main-text exact-reference duplication and `multiplex-limited phase` wording
  are removed.
- Both main manuscripts compile with no undefined references, undefined
  citations, TeX errors, or overfull boxes.
- All four PDFs were rendered page by page and visually reviewed. Final generic
  `article` page counts are: Paper 1 main 12, Paper 1 supplement 10, Paper 2
  main 7, Paper 2 supplement 10. The Paper-2 count before this surgery was 11
  generic-article pages; final Optica pagination remains an M9 task.
- Paper-2 first-six-page stop check: `jitter` 24 occurrences and `information`
  25, versus `certificate` 0, `fallback` 0, and `counterexample` 1. Figures
  1/2/3 occur on pages 2/4/5, respectively; the direct `19.13×` plot is before
  geometry discussion and references.

## Files changed intentionally

- `paper/main_oe.tex`
- `paper/supplement.tex`
- `paper/refs.bib` (two commented placeholder stubs made BibTeX-safe)
- `paper/main_oe.pdf`
- `paper2/main_m1.tex`
- `paper2/supplement_m1.tex`
- `paper2/main_m1.pdf`
- this completion letter

Existing generated `paper*/tmp`, result logs, data caches, and unrelated
untracked artifacts were neither restored nor staged.
