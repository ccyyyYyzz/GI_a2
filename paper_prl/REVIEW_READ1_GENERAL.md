# REVIEW_READ1_GENERAL — PRL general-physics interest read (R43 §7.3 step 7, read 1 of 2)

**Persona:** busy general physicist (condensed-matter / stat-mech), no ghost-imaging or
computational-imaging background, flipping through a PRL issue.
**Scope:** `paper_prl/main.tex` + `main.pdf` (6 pp) + `figures/fig1_hero.png`,
`fig2_jet_slopes.png`, `fig3_optical_realizations.png`. Report-only; no repo file edited.

---

## (a) VERDICT: PASS_WITH_FIXES

The central inversion lands from the title and Fig. 1, the opening hook is genuinely
physics-first, and the five displayed equations are the right five. Nothing here is a
BLOCKER: no factual error, and the "disorder erases the image but not local change" thesis
is legible. But four MAJOR general-interest issues keep it from a clean PASS — an abstract
that clots with coined jargon exactly where the ten-second reader is deciding whether to
care; a hero figure that floats to page 3 and is never co-visible with the abstract; a
momentum valley on page 2 that ends in a detector-calibration number-wall; and a back half
whose register drifts from "broad physical principle" toward "single-pixel imaging-methods
paper." All four are fixable without touching the science.

---

## (c) TEN-SECOND TAKEAWAY I ACTUALLY GOT (verbatim, before reading the body)

> "In a strongly scattering medium you can't form an image any more — the reconstruction
> collapses to basically one degree of freedom (rank 1) — but you can still *detect* that
> the scene changed, and there's a clean two-tier scaling law for how long that takes:
> T ~ ε⁻² for ordinary changes, ε⁻⁴ for special 'orthogonal' ones. So the *rank* of what you
> can reconstruct and the *observability* of a local change are two different invariants, and
> disorder kills the first without killing the second."

That is the intended inversion, and it landed. What I did **not** get in ten seconds, and had
to go into the body for: what a "jet" is, what "nuisance profiling" buys me, what
"amplitude-anchored" means, and what a "bucket-optics test" is. Those four terms sit in the
abstract and gate the read (see Finding 2).

---

## (b) FINDINGS

### Finding 1 — Hero figure floats to page 3; never co-visible with the abstract (MAJOR)
**Q1/Q3.** The ten-second test is *title + abstract + Fig. 1*. In the actual typeset PDF the
abstract is on page 1 (bottom of the abstract block, `main.tex` L31–40) but Fig. 1
(`\begin{figure*}[t]`, L165) lands at the **top of page 3**. Page 2 carries no float at all.
A physicist flipping through sees the abstract and must turn two pages to reach the one figure
that is engineered to deliver the payoff in ten seconds. R43 §7.3 step 2 makes Fig. 1's
ten-second legibility the load-bearing test; the float placement silently defeats it for a
browse reader.
*Repair (do not apply):* force Fig. 1 to the top of page 1 or 2 (e.g. move the `figure*`
block ahead of "The disorder paradox", or `[t!]`/manual placement), so it is co-visible with
the abstract. Figs. 2–3 clustering on page 4 is fine; Fig. 1 is the one that must migrate
forward.

### Finding 2 — Abstract clots with coined/undefined jargon at the decision point (MAJOR)
**Q2.** The abstract's middle sentences (L34–39) are where a non-imaging physicist decides to
read on, and they are the densest with undefined terms. Terms used before they are defined:
- **"nuisance profiling"** / **"nuisance-profiled"** (L34) — statisticians know "nuisance
  parameters," but a condensed-matter reader will not parse "profiling" as a projection.
- **"first nonzero Chernoff jet"** (L34–35) — *"jet"* is a coined contact-order noun, never
  glossed; "Chernoff" is unexplained. This is the single most opaque phrase in the abstract.
- **"amplitude-anchored"** (L36) — this is boundary 2 of the theorem; in the abstract it is a
  bare modifier with no referent.
- **"first-order-orthogonal"** (L37) — guessable (orthogonal to the gradient) but undefined.
- **"Exact divergences"** (L38) — a stat-mech reader reads "divergences" as singularities, not
  as KL/Chernoff statistical divergences. Genuinely ambiguous (see Finding 6).
- **"sealed bucket-optics test"** (L39) — *"bucket"* is a single-pixel-detector term; a general
  physicist will not know it.
- Symbols **T** and **ε** appear (L37–38) with no in-abstract gloss (detection time /
  perturbation amplitude).

Exact quoted sentence: *"After nuisance profiling, the first nonzero Chernoff jet of a
perturbation fixes its detection-time exponent. In an amplitude-anchored scrambled channel
the scene Fisher matrix collapses to rank one..."*
*Repair:* keep sentences 1–2 (they are excellent and jargon-free), then gloss the terms in
plain physics, e.g. "…the lowest order at which a scene change perturbs the measurement
statistics (its *jet order*) fixes the detection time"; replace "bucket-optics" with
"single-pixel optical" and drop or gloss "amplitude-anchored." Note the abstract is otherwise
disciplined — "banks," "beyond-band," and "DMD" do **not** appear in it (they do appear in the
body/figures); the fix is small and local.

### Finding 3 — Momentum valley: page 2 is all-text and dies in the sealed-detector number-wall (MAJOR)
**Q5.** Attention first sags at the cross-reference catalog (L96–99, Finding 10) but survives.
It **dies** in the "Sealed optical instantiation" paragraph (L271–287), and the layout makes
it worse: page 2 carries *no figure* (Fig. 1 is on p3, Figs. 2–3 on p4), so the reader crosses
"Decisive tests" (L232–248) and "Sealed optical instantiation" back-to-back as an unbroken
column of numbers. The sealed paragraph reads, in one breath: *"77/81 cells; the best cell
required 453 banks and achieved a Monte Carlo power lower bound of 0.990; a repeated code bank
retained a latency advantage of 459 versus 1013 banks … four-class balanced accuracy was 0.916
… false-alarm rates of 0.096 and 0.084 … aperture-preserving medium-basis rotations up to 10%
(FA=0.032, latency inflation 19%). At 20% rotation and a spectral-slope mismatch of −1,
non-target false alarms rose to 0.052 and 0.076…"* For a reader with no imaging background,
"banks," "cells," "beyond-band," "aperture-rotation," and "FA" are opaque, and the density of
calibration figures reads as a detector-engineering report with no physical narrative. This is
the paragraph where a general reader stops.
*Repair:* cut the sealed paragraph to 2–3 sentences of physics (power-certified beyond-band
alarm; measured attribution leakage as the honest boundary) and push the calibration numbers
into End Matter C / Fig. 3 caption. Consider letting a figure float onto page 2 to break the
column.

### Finding 4 — Register drift: the back half reads as an imaging-methods paper (MAJOR)
**Q7, brutal answer.** The first ~two pages are unambiguously "a broad physical principle
instantiated in optics": the title, the abstract's "distinct invariants," the opening
(disorder / speckle / memory effect), the reparameterization-invariant Rank–Jet theorem, and
the Outlook's "Disorder Observability Triangle" transferring to disordered sensing, change
detection, singular models, and quantum-inspired superresolution. That is PRL register.
But from "Scrambling inversion" onward the paper spends its most quantitative real estate on
single-pixel-detector methodology — "single-pixel bucket," "banks," "cells," "code band,"
"grain band," "photons per bucket," "four-class balanced accuracy," "false-alarm rates,"
"fixed-versus-fresh latency," "aperture-rotation domain" — and **Fig. 3 is a
detector-calibration dashboard** (a banks heatmap, a false-alarm bar chart against a "frozen
5% bar"). Combined with End Matter C's two integrity disclosures (a shot-noise error that
"inflated covariance Fisher information by approximately 1.9×"; a generator refactor), the back
half reads like the validation section of a computational-imaging methods paper bolted onto a
statistics theorem. A hostile general physicist would say the theorem could be stated with no
optics at all, and the optics that *is* here is mostly detector engineering. It does not fail —
the principle is still the star — but the balance leans too far into single-pixel-detector
calibration.
*Repair:* keep the optical *endpoint* prominent (it is the right PRL move) but demote the
calibration/leakage numbers to End Matter/SM, and reframe Fig. 3 so its physics panel
(complete-scramble: image dies, detection lives) leads and the calibration panel supports,
rather than the reverse weighting a methods reader would infer.

### Finding 5 — "jet" and "quotient" are undefined inside Fig. 1 (self-containment gap) (MINOR)
**Q3.** Fig. 1 is otherwise strongly self-contained: the three rows (mean support Bp→{0} "DC
only"; reconstruction Fisher rank "many"→"1 (Gx), 99.99% energy in one mode"; local change jet
"finite"→"finite — SURVIVES"), the "complete-scrambling law" box, the amplitude-anchor icon,
and the bottom banner ("Strong disorder erases coordinates … but not local change
observability: the quotient jet order stays finite") together tell the law without the body.
The gap: the row label **"local change quotient jet order"** and the banner's **"quotient jet
order"** use two coined nouns ("quotient," "jet") that a reader meeting the figure first cannot
decode. `O∘O` and `R` in the law box are also unglossed (minor).
*Repair:* in the figure, replace or gloss "quotient jet order" with something like "local
change order (m)" and let the m=1/m=2 lines carry it; the concept is already visually shown by
the two lines, so the coined word is not doing needed work in the graphic.

### Finding 6 — "Exact divergences" is ambiguous to a condensed-matter reader (MINOR)
**Q2.** Abstract L38 and body ("forming exact Gaussian divergences," L236–237) use "divergences"
in the information-theoretic sense (KL/Chernoff). A stat-mech reader's default reading of
"exact divergences" is *singularities*, which momentarily derails the sentence.
*Repair:* write "exact statistical divergences (KL and Chernoff)" on first use, or "exact
Gaussian KL/Chernoff distances."

### Finding 7 — Fig. 2 plotted slope label disagrees with the text exponent (MINOR)
**Q6/consistency.** Fig. 2 plots "slope 2.04" and its subpanel title reads
"single-look d′~ε^m (m = 0.95, 2.05)," while the body (L238) and Fig. 2 caption (L221) state
"exact contact orders 2.038 and 4.000." 2.04 vs 2.038 is only rounding, but a careful reader
cross-checking figure against text notices the mismatch.
*Repair:* make the plotted label read "slope 2.038" (or state "≈2.04") so figure and text
agree to the same precision.

### Finding 8 — Eq. (5) is the most optics-internal of the five displays (MINOR)
**Q6.** The five numbered equations each earn display: Eq. (1) T_req≍ε⁻²ᵐ (central law),
Eq. (2) the rank bound, Eq. (3) Σ=R+Q(x)H with ΔQ, Eq. (4) the m=1/m=2 dichotomy — these four
are followable by a general reader. Eq. (5), 𝔼[b_i]=C(0)O_ii x̂(0) (L191), introduces
imaging-specific notation (O_ii, x̂(0), C(0)) and is the one display a non-imaging reader
cannot parse from the main text. It earns its place for the "mean sees only DC" point but is
the weakest of the five for this audience.
*Repair:* optional — keep it, but state in one clause that O_ii is the code weight and x̂(0)
the scene DC component, so the equation is readable in place.

### Finding 9 — The Q → m → T chain has one unproven hand-wave in the main text (MINOR)
**Q6.** The chain Q=xᵀGx → ΔQ=2xᵀGδ+δᵀGδ (Eq. 3) → m=1/2 (Eq. 4) → T∝ε⁻²ᵐ (Eq. 1) *is*
followable from the main text without the supplement, with one exception: L136,
*"Because a regular reduced law inherits the contact order of Q, every direction falls into
exactly one of two local classes,"* is the single step a general reader must take on faith
(the justification lives in End Matter A). It is stated, so this is minor.
*Repair:* add a half-sentence pointer ("— the reduced law is analytic in Q, End Matter A —")
so the reader knows the gap is closed, not skipped.

### Finding 10 — Cross-reference catalog is the first (survivable) momentum sag (MINOR)
**Q5, secondary.** L96–99: *"The quotient jet is a stochastic observability order, the
statistical counterpart of successive-derivative observability rank conditions [10], and the
efficient profiled score is the unique statistic that saturates the fluctuation–response bound
[8]."* Three cross-field name-drops (Hermann–Krener observability, Dechant fluctuation–response)
in one sentence is where attention first dips on page 1. It recovers, so minor, but it is a
candidate to lighten.
*Repair:* split into two sentences or move the fluctuation–response tie to a footnote/EM.

---

## Question-by-question summary (persona answers)

1. **Ten-second test:** Yes, the inversion lands — from the title and Fig. 1's banner. It
   stumbles on (i) the abstract's coined jargon (Finding 2) and (ii) the layout: Fig. 1 is on
   page 3, so the browse reader never sees title+abstract+Fig.1 together (Finding 1).
2. **Abstract jargon:** terms used before defined — *nuisance profiling, (Chernoff) jet,
   amplitude-anchored, first-order-orthogonal, exact divergences (ambiguous), bucket-optics*,
   plus undefined symbols T, ε. "banks/beyond-band/DMD" are correctly kept out of the abstract.
3. **Fig. 1 self-contained:** Yes, almost fully — three rows + law box + banner tell the law.
   Only gap: the coined nouns "quotient/jet" are unglossed in the graphic (Finding 5).
4. **Opening hook:** PASS. The opening is physics-first (scattering, speckle, memory effect →
   "received wisdom … We show that this conflates two logically independent quantities"), frames
   observability as an invariant with an order, and the Feng / Lee–Stone references land in the
   target reader's own (mesoscopic) backyard. Good.
5. **Momentum:** dies in "Sealed optical instantiation" (L271–287), amplified by page 2 being
   all-text (Finding 3); first sag at the L96–99 catalog (Finding 10).
6. **Five equations:** the right five; the Q→m→T chain is followable without the supplement
   (one stated hand-wave, Finding 9); Eq. (5) is the most optics-internal (Finding 8).
7. **PRL register:** broad principle for ~2 pages, then drifts toward a single-pixel
   imaging-methods paper in the sealed section and Fig. 3 (Finding 4). Passes, but the balance
   over-invests in detector calibration.

---

### Severity tally
- BLOCKER: 0
- MAJOR: 4 (Findings 1, 2, 3, 4)
- MINOR: 6 (Findings 5, 6, 7, 8, 9, 10)
