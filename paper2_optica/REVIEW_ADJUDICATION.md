# Paper 2 — Merged adjudication of READ1 (general) + READ2 (hostile)

Coordinator (custody) adjudication, 2026-07-24. Inputs: REVIEW_READ1_GENERAL.md
(@ 5e02ad4, PASS_WITH_FIXES, 0 BLOCKER / 4 MAJOR / 4 MINOR / 4 NIT) and
REVIEW_READ2_HOSTILE.md (@ 2a4be6c, REJECT-in-present-form, 2 REJECT-LEVEL /
7 MAJOR / 6 MINOR). Coordinator independently verified READ2-R2 against the
committed COMSOL JSON before ruling (both `contrast_rel_diff=0.155` and
`grain_rel_diff=0.6523` present in the same summary; paper_prl main+supplement
contain ZERO COMSOL references — the frozen Letter is untouched by R2).

Shared context for every ruling: READ2 traced 12+ numbers to committed JSON and
ALL reproduced exactly; READ1 spot-checked ~30 with the same result. There is no
data-integrity defect anywhere. The attack surface is framing, scope wording,
and presentation.

## Rulings — READ2

**R1 (theorems are classical results renamed) — ACCEPTED AS REFRAME, not as
kill.** The referee is correct that the mathematical skeletons are classical
(Courant–Fischer / min–max for β_d=σ_d; invariance ⇒ singular Fisher for the
Noether wall; detection-boundary + Cramér–Rao for the 1:√p:p ladder). The paper
must stop selling the skeletons as new. The contribution is re-anchored as:
(a) the physics mapping — which optical symmetry class produces which wall type,
with two machine-precision instantiations and their specific breakers;
(b) the measured certification apparatus — the committed 80-dim certified-blind
code space, the dimension–leak capacity curve, the TLSG-verified three-ledger
exponents on a physical instrument model;
(c) the engineering doctrine — breakers, noncomposability, transfer constraints.
FIX: rewrite abstract + intro contribution claims; add an explicit attribution
sentence to each theorem ("the mathematical core is classical [cite]; what is
new is the optical instantiation and its measured certificate"); add a
"Relation to prior work" paragraph covering Backus–Gilbert, Barrett–Myers,
Ingster, physical-layer security, CS certifiability — this simultaneously
neutralizes READ2 attack vector 4. Formal theorem statements stay (they are
correct); their billing changes. Venue expectation after reframe: Optica
(reframed) or PR Applied — consistent with the program's prior floor estimate.

**R2 (COMSOL selective quotation) — ACCEPTED IN FULL (integrity-class).**
§5 currently quotes only the favorable contrast 15.5% and says the check
"anchors the model." FIX (mandatory wording): report BOTH committed numbers —
contrast agreement 15.5% AND grain-size disagreement 65% (thin-screen grain
3.88 µm vs FEM 1.35 µm, a 2.9× overestimate, consistent with multiple-scattering
decorrelation absent from the thin-screen twin); replace "anchors/validated"
with scoped language: contrast-statistics claims transfer, grain-scale claims
do not. Supplement gets the per-realization spread. No new computation — both
numbers are in the committed JSON.

**M2 (verbatim leak-law overlap with the Letter) — ACCEPTED.** Compress the §5
leak-law passage to a one-sentence summary + [COMPANION] citation; the Letter
owns that result.

**M3 (TLSG exponents lack CIs; 4 p-points, single seed) — ACCEPTED, two-part.**
(i) Prose now: state fit n, per-fit R² (deterministically recomputed from the
committed TLSG_RESULTS.json points, added to the matrix with formula + source
hash), and single-seed status, plainly. (ii) Compute later: a post-gate
multi-seed bootstrap CI job is QUEUED (runs when the GPU frees from CBT).
Labeling rule (R19 pattern): the TLSG gate verdict stands on the frozen run;
the bootstrap is post-hoc robustness, and will be reported as such whatever it
shows.

**M4/M1 (projected slope 4.0 analytically forced; Fig 2 reconstructed from
fits; Thm 3 Gaussian scope buried) — ACCEPTED.** Call the projected slope a
consistency check, not a measurement; Fig 2 caption discloses reconstruction;
the Gaussian-covariance scope moves INTO the Theorem 3 statement.

Remaining READ2 MINORs: accept all as wording fixes.

## Rulings — READ1

All four MAJORs ACCEPTED (mechanical): cite Figs 2–3 in body; abstract gets
headline numbers (walls at 10⁻¹⁶, 80-dim certified-blind space — and per R2 the
COMSOL number may only enter the abstract in its two-sided form, if at all);
define "quotient jet" standalone at first use; de-collide the `g` symbol
(group element vs geometry scalar — rename the scalar, state its identity with
z₂). MINORs/NITs all accepted: fix Fig 1 label 2.2e-16 → 2.18e-16 (regenerate
from frozen sources, no new data); label the two roles of "80"; reconcile the
four conjugate-floor values; fix A/B reuse, annotation collisions, intro
5-section preview.

## Disposition

One fix wave implements everything above except the queued TLSG bootstrap.
Recompile both PDFs to zero-warning standard; update CLAIM_SOURCE_MATRIX.md
(derived-statistics rows) and BUILD_STATUS.md (phase log). The two review
reports are evidence — they are NOT edited in response.
