# Letter 10 — R20 presentation overhaul: manuscript surgery is yours

From: main agent. Date: 2026-07-19. Authority:
`docs/ROUND63_GPT_ROUND20_REVIEW_RAW.md` (issue #12, commit 90f3b67) — read it
in full first; it is unusually concrete. The operator's "哪哪都不对" was
adjudicated as **hierarchy failure**, and I concur with essentially every
item. Execution split: you do the TEXT surgery (this letter); a figure agent
is concurrently doing the FIGURE rebuilds (M3/M4/M5/M6-recomposition/M10
balance/M11/M12 — its files: paper*/figs/*, code/round63/figs/*; do not edit
those); the Optica template port (M9) is a separate step blocked on a user
download approval and will follow your surgery.

Your work order (R20 item numbers):

1. **M1 — retitle paper 2**: "Jitter-capped high-flux single-pixel imaging
   with global power control". Sweep title-dependent references.
2. **M2 — paper-2 abstract rebuild**: exactly the seven frozen beats,
   190–230 words; the deletion list is explicit. Corrected values + the
   power-for-time sentence stay; correction narrative OUT of the abstract.
3. **M6/M7 — process-residue relocation**: move the listed blocks to the
   supplement (retirement chronology, PATH_FEASIBLE_ALPHA /
   ADAPTIVE_COLLAPSE mentions beyond first definition, 12-row DEV table,
   conic construction + solver protocol, accounting/SHA/wall-cap/retry
   details, repeated branch explanations); replace with the two compact
   main-text paragraphs specified in M7. Preserve every R18/R19 frozen
   label/caveat in the supplement (conflict table §9 rules this explicitly).
4. **M8 — move the R19 correction disclosure VERBATIM** to a "Analysis
   correction and provenance" subsection at the end of campaign/analysis
   methods, immediately before Results; Results keep corrected values + one
   parenthetical pointer. Quiet typography (§5): no box, no drama.
5. **§6 section titles**: adopt the seven recommended reader-facing titles;
   "Act III"/doctrine/handbook/movement vocabulary out of reader-facing
   text (internal script names may stay in code/READMEs).
6. **M15 — benefit+cost co-location**: ensure the results prose/table puts
   19.13×, +1.87 dB, incident-dose 37.1×, detected-count ratio in one
   place (the figure agent handles the visual counterpart; your table per
   SHOULD-FIX 14: split ΔQ and S_gate columns, emphasize cost ratio).
7. **SHOULD-FIX text items**: 1 (humanize code labels — token once in
   parentheses, exact strings live in Methods/tables/supplement), 2 (round
   main-text numbers: 1.41 dB, 19.13×, 18.33×), 5 (add a short explicit
   Limitations subsection to EACH paper: peak power, jitter calibration,
   simulation-only), 6 (scene IDs out of main prose), 11 (provenance paths
   out of journal captions), 12 (paper-1: compress natural-cohort /
   mismatch / exact-reference to one paragraph each, detail to supplement),
   13 (drop the per-family table if the strip remains), 15 (captions state
   results; selection defence to supplement).
8. **Paper-1 items that need no amendment**: M14.1 (Study-1 composite
   placement), SF3 (delete "costs essentially nothing" duplication), SF4
   ("multiplex-limited phase" → "regime").
9. **HOLD paper-1 abstract (M13)** until the R21 presentation amendment
   lands (I am requesting it from GPT now — R9 froze that block; §9 row 1).
   Also HOLD any interpretation question on disclosure placement — R21
   covers it.
10. **Coordination**: a submission-prep agent may have left uncommitted
    back-matter edits (Optica-ordered Funding/Disclosures/Data-availability
    blocks) and cover-letter drafts; absorb compatible pieces, supersede
    conflicts, note what you did.
11. **Stop condition** (R20 §4): in paper 2's first six pages,
    certificate/fallback/counterexample/DEV/audit tokens must not outweigh
    jitter data + image-quality curves; verify before reporting.

Compile clean after each major step; commit as yourself; letter 11 with the
before/after map and page counts. The figure agent's new PDFs will land
under the existing filenames plus new ones (fig_jitter_validation,
fig_speed_results) — reference them where R20 places them; if a figure is
not ready at your compile time, use a draft box placeholder and note it.