# R21 ruling (GitHub issue #13, raw)

Title: R21 ruling: abstract relocation and correction-disclosure placement
Posted: 2026-07-19T19:51:46Z

---

# R21 — Presentation amendments: abstract relocation and correction-disclosure placement

**Reference commit:** `aac9d53`.

**Scope:** presentation location and hierarchy only. No result, endpoint, number, label, or scientific interpretation is changed.

> **R21 VERDICT: both amendments are approved.** R9 froze the Study-2 verdict hierarchy and factual content, not its permanent residence in the abstract. R19 froze the correction paragraph's wording and disclosure function, not its location inside the Results section.

## A1 — Paper-1 abstract relocation

### Ruling

**Approved.** The R9 Study-2 verdict block may be removed from the abstract and relocated to a dedicated main-text Results/verdict subsection. This satisfies R9's intent provided that the relocated text preserves the complete hierarchy and facts verbatim or semantically exactly, and appears before any mechanism interpretation that might soften the negative primary.

The dedicated Results subsection must begin with the R9-required sentence:

> **The preregistered Study-2 primary endpoint was formally negative because its breadth criterion failed: 16 of 24 scenes accelerated, below the required 18 of 24.**

It must then preserve all of the following:

1. median `S_gate = 6.796` and family-stratified bootstrap lower bound `2.835`;
2. fixed-dwell secondary median `+1.16 dB`, `19/24` positive, lower bound `+0.72 dB`;
3. the post-unblinding six-family ordering by mean measured bucket contrast, with the frozen caveat that it was not preregistered, entered no decision rule, and does not estimate a universal critical contrast;
4. the no-gate occupancy result: fixed-dwell gain was nonmonotone and largest **among the tested rungs** at `k=32` (`3.1%` occupancy), revealing a contrast–conditioning trade-off—not an “optimal sparsity” claim;
5. the Study-1 formal negative/uninformative result and its small positive fixed-dwell secondary remain in the Study-1 Results subsection, so both preregistered negative outcomes remain prominent in the main paper.

The relocation is not permission to paraphrase away “formally negative,” combine the primary and secondary into a rescued positive verdict, or move the complete verdict record to the supplement.

### Frozen five-beat abstract skeleton

Target: **220–260 words**. Use five short paragraphs or five clearly separated sentence groups:

1. **Operating question.** State the low-flux convention and ask when deliberately entering the nonlinear dead-time regime improves single-pixel imaging.
2. **Method and law.** Introduce the convex renewal quasi-likelihood reconstruction and the finite-window count-information ridge; give only the cube-root law and the role of bucket contrast.
3. **Dense-pattern regime.** State that the dense-pattern time-to-quality endpoint was uninformative/formally negative, while a small positive fixed-dwell gain remained. One headline number is sufficient; detailed gate arithmetic stays in Results.
4. **Sparse equal-load regime.** State that the Study-2 acquisition-speed primary remained formally negative on breadth, whereas the fixed-dwell secondary passed and the response exposed a contrast–conditioning trade-off. Do not reproduce the full R9 verdict block here.
5. **Conclusion and boundary.** State the operating-map conclusion and the two practical limits: sparse equal-dose operation requires increased peak irradiance, and the evidence is simulation-only pending bench validation.

Remove from the abstract:

- detailed gate arithmetic and all effect-size criteria;
- multiple bootstrap bounds;
- the complete six-family ordering;
- the no-gate rung-by-rung result;
- dynamic-scattering details;
- provenance or freeze-process language.

The abstract must still use the honest hierarchy “primary formally negative; fixed-dwell secondary passed.”

## A2 — R19 correction-disclosure placement

### Ruling

**Approved.** R19 §1.4 binds the correction paragraph's **wording**, not its exact location. Move that paragraph **verbatim** to a quiet subsection titled:

> **Analysis correction and provenance**

Place it at the end of the campaign/analysis methods, immediately before `Results`. Keep normal body typography—no warning box, icon, oversized rule, or repeated disclosure.

In the elapsed-time Results subsection, retain the corrected values and add exactly one neutral pointer, for example:

> *(The post-unblinding implementation-to-specification correction and its provenance are reported in the preceding “Analysis correction and provenance” subsection.)*

The full paragraph must not be repeated in the abstract, Results, Discussion, captions, or Conclusion. The immutable original outputs and dated correction bundle remain cited in the relocated subsection.

This placement preserves R19's purpose: readers encounter the complete disclosure before the corrected results, while the Results sequence is no longer interrupted by process narration.

## Authorization

A1 and A2 are now presentation-authorized. No further ruling is needed to execute these two relocations, provided the mandatory R9 facts and the R19 paragraph remain exactly within the boundaries above.