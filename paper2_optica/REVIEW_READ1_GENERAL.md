# REVIEW_READ1_GENERAL — Optica Research Article, constructive general referee (read 1 of 2)

**Persona:** an Optica referee with an optics / computational-imaging background but **no**
prior exposure to this program's 17-test history and **without the companion Letter in hand**,
reading the manuscript cold for coherence, self-containment, and internal consistency.
**Scope:** `paper2_optica/main.tex` + `main.pdf` (6 pp) + `supplement.tex` + `supplement.pdf`
(6 pp) + `figures/fig1_hero`, `fig2_jet_transmutation`, `fig3_routes_transfer`.
Cross-checked against `CLAIM_SOURCE_MATRIX.md`, `BUILD_STATUS.md`, and (claim-boundary only)
`paper_prl/main.pdf`. **Report-only; no manuscript file edited.** READ2 handles the hostile pass.

---

## (a) VERDICT: PASS_WITH_FIXES (accept-trajectory)

This reads as a genuine Optica Research Article, not a bag of stapled results. The five-section
arc the venue wants — **symmetry walls → blindness capacity → restoring the jet → what finite
data can certify → wave-optics transfer** — is real, and it is carried by *prose*, not just by
section titles: every section opens with a one-clause recap that positions it against the
previous one ("An exact wall is a limiting case"; "The exact walls of Section 2 are the
endpoints of a continuous mechanism"; "A wall says a change is invisible; blindness capacity
says how large a subspace is guaranteed invisible; the jet order says how a visible change
scales. A fourth, independent price…"; "The theorems above are algebraic; a bench is not").
The hero figure lands at the **top of page 2**, co-visible with the abstract and intro (this
is the exact layout defect the companion Letter's READ1 flagged, and it is fixed here). I
spot-checked ~30 numerical claims across main + supplement against the claim–source matrix and
**every one traces correctly**; the three known GAP items (6.6e-16 wall, 35@1e-4-vs-80, the two
"→24" instances) are all handled with correct provenance. There is **no BLOCKER**: nothing is
factually wrong, the LaTeX compiles clean, the honest negatives (S5) carry KILL verdicts with
the same detail as the positives and are referenced from the Discussion.

What holds it back from a clean PASS are four MAJOR issues, none touching the science: (1) **Figs
2 and 3 are never cited by number anywhere in the body** — only Fig 1 is called out; (2) the
**abstract states no quantitative result at all** — it is entirely qualitative where an Optica
abstract should carry a headline number; (3) the load-bearing coined term **"quotient
(information) jet" is used in Theorem 1 and Theorem 3 without a self-contained definition**,
silently outsourcing it to the companion Letter a standalone reader may not have; and (4) the
symbol **`g` is overloaded** across two theorems (group element `g∈G` in Thm 1 vs a gain/geometry
scalar `g` in §4/Thm 3, never physically identified). All four are local fixes.

---

## (b) FINDINGS

### Finding 1 — Figures 2 and 3 are never referenced by number in the body (MAJOR)
**Q3/Q7.** `grep` of `main.tex` for figure cross-references returns exactly two hits, both
`\ref{fig:hero}` (L268, L426). **`fig:jet` (Fig. 2) and `fig:transfer` (Fig. 3) are never
`\ref`-ed in the running text.** Fig. 2 is the jet-transmutation figure for §4, whose Routes 1–3
prose (L330–361) never says "Fig. 2"; Fig. 3 is the routes/calibration-transfer figure for §6,
whose prose (L449–477) never says "Fig. 3." In the typeset PDF this is visible as two full-width
floats (p. 4 top, p. 6 top) that a reader arrives at with no textual anchor telling them when to
look or what they add. Optica production requires every display item to be called out in text in
order; a copyeditor will bounce this, and a reader loses the figure↔claim binding.
*Repair (do not apply):* add an in-text callout at the natural place — e.g. at the head of the
Routes paragraph "(Fig.~\ref{fig:jet})", and in §6 at the leak-ladder / calibration-transfer
sentences "(Fig.~\ref{fig:transfer})". No new content, one `\ref` each.

### Finding 2 — Abstract carries no quantitative result (MAJOR)
**Q2.** The abstract (L75–84, ~98 words, correctly at the Optica ~100-word target) states each
theorem's *form* well — score-null orbit, `β_d=σ_d`, the `1:√p:p` ladder, jet transmutation — but
contains **zero numbers**. A reader deciding whether to read on gets no magnitude: not the
machine-zero walls (`~10⁻¹⁶`), not the measured capacity (80 code dimensions certified invisible),
not the wave-optics validation (15.5%). The companion Letter's abstract, by contrast, gives the
reader `T~ε⁻²`/`ε⁻⁴`. For an Optica Research Article the abstract should name at least one headline
result. The closing sentence "This opens calibrated, publishable optical blindness" is also the
vaguest place to spend the last 6 words.
*Repair:* fold one or two concrete numbers into the existing sentences, e.g. "…a measured singular
spectrum (80 code dimensions certified invisible below 10⁻⁵ at a fixed plane)…" and "A scaling
gate and a full-Maxwell wave-optics twin (15.5% contrast agreement) validate the program."

### Finding 3 — "quotient (information) jet" is undefined; self-containment gap vs the companion (MAJOR)
**Q4.** The noun **"jet"** is load-bearing and never glossed in this article. It appears as
"nuisance-profiled statistical jet" (intro, L101), "quotient information jet" inside
**Theorem 1** (L173), the corollary title "Jet transmutation" (L317), and "Composed with a
quotient jet" inside **Theorem 3** (L387). The main text only ever defines the *effective
amplitude order* `m_eff` and the leading orders `m₁,m₂` (§4) — it never tells a cold reader what
a "jet" *is*. Verified against `paper_prl/main.pdf`: the companion Letter is where "the quotient
information-jet order `m=m(δ)`" is actually defined. An Optica article must stand on its own; a
reader without the Letter meets an undefined coined term in the paper's very first theorem
statement.
*Repair:* one clause at first use, e.g. "…the first nonzero *jet* — the lowest order at which a
scene change perturbs the nuisance-profiled measurement statistics — fixes the price…", and a
parenthetical in Theorem 1 tying "quotient information jet" to that definition.

### Finding 4 — Symbol `g` is overloaded across two theorems (MAJOR)
**Q4.** `g` is used with two incompatible meanings in adjacent formal statements:
- **Theorem 1 / §2:** `g∈G`, a **group element** (`P_{g·x}=P_x for every g∈G`, L167–170).
- **§4 / Theorem 3:** `g` is a continuous **gain/geometry scalar** (`C(ε,g)=A g^{2a₁}ε^{2m₁}+…`,
  L306; `A_⊥=κ g^{a} ε^{m} B`, L388), and its physical identity is never stated in the main text —
  only Supplement S3 reveals that in the near-contact case it is effectively `z₂`.
A referee reading the theorems in sequence will trip on the collision, and "what is `g`?" is a
genuine unanswered question in §4/§5.
*Repair:* rename the gain/geometry variable (e.g. `γ` or reuse `z₂` explicitly where the physics
is near-contact), and add one clause at L306 identifying it ("a geometric gain, here the
contact distance `z₂`"). Keep `g` for the group element only.

### Finding 5 — Fig. 1 graphic contradicts its own caption on the energy-wall value (MINOR)
**Q3/Q5.** The Fig. 1 graphic labels the energy (unitary) wall **"null 2.2×10⁻¹⁶"** (left panel),
while the caption (L144) and body (L205) both say **"2.18×10⁻¹⁶"** (matrix A3: 2.1798e-16). This
is an internal inconsistency *inside the same figure* — the one figure-vs-text drift I found in
the whole manuscript. (The support wall "6.6×10⁻¹⁶" is consistent everywhere, GAP-1 correctly
honored.)
*Repair:* regenerate `fig1_hero` with the label reading `2.18×10⁻¹⁶` (or make the caption say
`≈2.2×10⁻¹⁶`) so figure and caption agree to the same precision.

### Finding 6 — The value "80" plays two distinct roles a reader can conflate (MINOR)
**Q5 (GAP-2).** `80` appears both as the **single-plane deterministic capacity** (80 dims certified
`≤10⁻⁵`, C1/IT4b, §3 L261) and as the **fresh-draw minimum** (fresh `d@10⁻⁴` min = 80, C6/TLSG bar
7, §3 L282 and Fig. 3 right). The provenance is technically correct and the hero caption does
label the joint-35 and single-plane-80 points distinctly — but the same integer standing for two
different certificates (deterministic ≤1e-5 vs fresh-draw ≤1e-4) is a trap for a careful reader,
especially since Fig. 3's caption says "at least 80 code dimensions remain certified below 10⁻⁴"
(the fresh-draw min) while the fig-3 bar's "80 dims" is the *deterministic ≤10⁻⁵* capacity.
*Repair:* one half-sentence disambiguation on first collision, e.g. "(the fresh-draw minimum of
80 coincides numerically with, but is a stronger statement than, the single-plane deterministic
capacity)"; or annotate the two 80s with their thresholds wherever both appear.

### Finding 7 — Conjugate-plane floor appears as four slightly different numbers in §6 (MINOR)
**Q5.** Within §6 the conjugate-plane leak floor is quoted as **6.94×10⁻³** (leak-law intercept,
Eq. 7, L453), **7.15×10⁻³** (instrument constant `L_pix`, L459), **6.8×10⁻³** (real-DMD pixelation
floor, L461), and **7.0×10⁻³** (analog PWM baseline, L477). Each is correctly sourced (IT3, IT2,
WAVE_TWIN, IT8 respectively), but presented back-to-back with no unifying remark they read as
measurement drift on "the same" quantity.
*Repair:* add a clause noting these are consistent independent estimates of one ~`7×10⁻³`
instrument floor from different tests/generators, so the reader sees agreement rather than scatter.

### Finding 8 — `A` and `B` are reused across sections (MINOR)
**Q4.** Beyond `g` (Finding 4): `A` = bandlimited mean operator (§2, L194), the leading Chernoff
coefficient (§4, L306), and (as `A_⊥`) the efficient displacement (§5). `B` = the energy-bucket
functional `B(E)=‖E‖₂²` (§2, L203), the second Chernoff coefficient (§4, L306), the certificate
metric `B†B` (Thm 2), and the normalized tangent operator in `A_⊥=κg^aε^mB` (§5, L388). Contexts
differ enough to disambiguate on close reading, but the reuse is dense.
*Repair:* at minimum a symbol table (Supplement) or one-word local reminders on reuse; ideally
distinct glyphs for the Chernoff coefficients (e.g. `c₁,c₂`).

### Finding 9 — Fig. 2 right panel: annotation text collides with the fixed-point label (MINOR)
**Q3.** In the typeset Fig. 2 (p. 4), the orange annotations "projected cov slope = 4.00 (all 4
cells)" / "unprojected low-ε slope = 2.12–2.20" overprint the gray dotted-line label "m = 1
leak-dominated fixed pt" at the bottom of the axes; the text is partly illegible.
*Repair:* reposition the annotations (or the gray label) in `fig2_jet_transmutation.py` so they
do not overlap; move the measured-slope note into clear space or into the caption.

### Finding 10 — Fig. 3 middle: rotated axis label overlaps the right-panel tick labels (NIT)
**Q3.** In Fig. 3 the shared-axis seam places the rotated y-label "bottom-64 subspace leak" over
the right panel's `10⁻¹…10⁻⁵` tick labels (visible p. 6). Cosmetic but it muddies an otherwise
clean two-panel figure.
*Repair:* add horizontal padding between panels or move the label in `fig3_routes_transfer.py`.

### Finding 11 — "seventeen-test verification battery" is asserted but never enumerated (NIT)
**Q7.** The intro (L130) promises "the seventeen-test verification battery… in the Supplement,"
but the Supplement is organized as S1–S6 by theorem and never labels or counts a set of 17; a
reader who tries to reconcile the number cannot. Minor, because it is a forward pointer, not a
claim the argument leans on.
*Repair:* either drop the specific count ("the full verification battery") or add a one-line
roster in S6 mapping the 17 test IDs.

### Finding 12 — Intro preview lists four items for five content sections (NIT)
**Q1.** The enumerated preview (L110–123) has four items mapping to §2–§5; §6 (wave-optics
transfer) is introduced only in the following sentence ("a wave-optics transfer chain (Section
6)"). Coherent, but a reader using the list as a table of contents will not find the fifth
section in it.
*Repair:* optional — either add a fifth bullet for the transfer/validation section or a half-clause
signaling that §6 is the physical-realization capstone.

---

## Question-by-question summary (persona answers)

1. **Coherent Optica arc, or stapled results?** Coherent. The five-section spine is real and is
   carried by genuine inter-section transitions, not just headings (Finding 12 is the only seam).
2. **Abstract (~100 words, theorems + headline numbers)?** ~98 words, states all three theorems +
   corollary content clearly, but **no numbers at all** (Finding 2), and the closing clause is
   vague.
3. **Figures land, captions self-contained, no caption↔body contradiction?** Fig. 1 lands
   perfectly (p. 2, co-visible with abstract) and is well cross-referenced; captions are
   self-contained. But **Figs 2 & 3 are never cited in text** (Finding 1), Fig. 1's graphic
   contradicts its own caption on 2.18 vs 2.2 (Finding 5), and two figures have annotation
   collisions (Findings 9, 10).
4. **Every symbol defined before use; consistent with the companion?** Mostly, but "quotient
   jet" is undefined (Finding 3), `g` is overloaded across theorems (Finding 4), and `A`/`B` are
   reused (Finding 8). Shared macros (`β_d`, `σ_d`, `ε`) are consistent with the Letter; the
   claim boundary ("that Letter owns the Rank–Jet Separation theorem and the algebraic mean
   wall") is stated correctly and does not contradict `paper_prl`.
5. **Internal consistency of ~20+ numerical claims vs the matrix?** Spot-checked ~30 (walls,
   capacity spectrum, routes, ledger exponents/CVs, KT6, leak law, COMSOL, M-scaling, negatives)
   — **all trace correctly.** The three GAP items are handled with correct provenance. Only
   residual issues are presentational: the 2.18/2.2 figure drift (Finding 5), the double-role of
   "80" (Finding 6), and the four floor values in §6 (Finding 7).
6. **Honest negatives (S5) symmetric and referenced from main text?** Yes. IT6 and IT7 each carry
   an explicit "Verdict: KILL" with the same detail as the positives, plus a deferred-future-work
   paragraph; the Discussion (L500–503) references S5 by name. Genuinely symmetric.
7. **Readability defects / assumed history?** The main ones are the two orphaned figures (Finding
   1) and the undefined "jet" (Finding 3). The "seventeen-test" pointer (Finding 11) is the only
   place that gestures at unshared history, and it is only a pointer.

---

### Severity tally
- BLOCKER: 0
- MAJOR: 4 (Findings 1, 2, 3, 4)
- MINOR: 4 (Findings 5, 6, 7, 8)
- NIT: 4 (Findings 9, 10, 11, 12)

### Summary verdict (as an Optica referee would frame it)
This is an accept-trajectory manuscript. The scientific content is sound, every number is
traceable to a committed artifact, the honest-negatives machinery is symmetric, and — unlike many
theory-into-optics syntheses — the five-part story is genuinely narrated rather than assembled.
I would return it for a **minor-to-moderate revision**. The four MAJOR items are all
presentational/self-containment fixes an author can complete in an afternoon: cite Figures 2 and
3 in the text, put one or two headline numbers in the abstract, define "jet" once, and
de-collide the `g` symbol. Clearing those, plus the four MINOR consistency/notation polish items,
would leave a clean, self-contained Optica Research Article. Nothing here requires new
computation or threatens the claims.
