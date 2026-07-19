# ROUND63 — GPT Round 22 consultation (pure theory exploration: is there ONE essential object under this system?)

Status framing: the two manuscripts are presentation-complete and FROZEN for
this round — nothing in this consultation modifies them. This is exploratory
theory for the NEXT work. The operator's charge, verbatim in spirit: the
system now feels "many but incomplete, connected but not unified — which is
ugly and cumbersome" (多而不完整,相连而不统一). The request is COMPRESSION:
find the single object from which the existing results fall out as
corollaries — or rule honestly that no such object exists at this depth.
Do not propose additional components; every new concept must retire at
least two existing ones.

## 1. The six results that want to be corollaries

(All frozen and evidenced in this repo; refs: R14/R16 theory chain, R11,
R17–R19, `results/round63_m1/`.)

- **E1 (exact identity).** I_N = E[N] − ρ²·E[Var(L|N)]: information lost to
  a dead-time detector equals the conditional variance of its hidden live
  time. Mechanism-level, exact.
- **E2 (factorized long-window law).** J∞(ρ,c) = ρ/[(1+ρ)(1+c²ρ²)]: the
  loss factorizes as (saturation) × (hidden-hold noise); optimum
  c²ρ²(1+2ρ)=1; universal exponents (1/3, −2/3) for regular finite-variance
  iid holds; verified out-of-sample (slope −0.658, constants 5–10%).
- **E3 (water-filling).** The variable-load relaxed OED optimum follows an
  information water-filling; naive per-pattern load equalization destroys
  brightness–information alignment (−7 dB, refuted preregistered).
- **E4 (collapse localization).** Under the deployed guard stack the
  adaptive design space collapses onto the balanced family — and the audit
  chain localized the operative restriction to the CONDITIONING anchors
  (A-risk/spectral), not dose uniformity (R18 primal constructions:
  dose-only class holds 2.2–2.5× D-eff headroom).
- **E5 (headroom existence).** Even the FULL safety stack retains 3.7–4.8×
  feasible D-efficiency headroom via scene-adapted geometry selection
  (support-expanding primal, every DEV family; 0/480 cells certified
  near-optimal on the confirmatory cohort).
- **E6 (empirical map + rank-one success).** One global power multiplier
  (the rank-one control) delivers +1.87 dB fixed-dwell and 19.13× elapsed
  time on 19–21/24 scenes; its failures (contour) correlate with bucket
  contrast / alignment (r = +0.64); the operating map's axes are (C_u, ρ).

## 2. Three candidate unifying structures (all currently UNREFINED — refute freely)

- **C1 — a product law for cascaded hidden-state losses ("Friis algebra").**
  Is the factorization in E2 structural rather than accidental? Conjecture:
  for a counting channel with several independent hidden-state
  nonidealities (dead time, hold jitter, afterpulsing, dark counts…), the
  long-window information efficiency factorizes (or composes via a simple
  algebra) with one factor per mechanism, each a conditional-variance term
  of E1 type. If true: detector design becomes noise-figure algebra, E1 is
  the generator, E2 its first two terms, and the universality question
  becomes a classification of factor exponents (a table of detector
  universality classes). Please derive the composition rule for at least
  one added mechanism (e.g. afterpulsing as correlated holds, or dark
  counts) — or exhibit a counterexample showing the factorization is a
  two-mechanism accident.
- **C2 — a control-rank hierarchy with an alignment criterion.** The global
  knob is the rank-one element of the control space (one scalar scaling
  all pattern energies). Conjecture: there is an alignment functional
  A_k(scene, patterns) such that rank-k control captures the achievable
  gain iff the information spectrum is A_k-aligned with brightness; E6's
  contrast axis is the first-order coefficient (C_u ≈ A_1), the contour
  failure is exactly A_1-misalignment, and E5's headroom is the rank-∞
  residual. If formalizable, paper 1's empirical map becomes the
  first-order term of paper 2's theory — the missing bridge. Please
  propose the correct functional (candidates: the overlap between the
  brightness vector and the leading task-subspace eigenvectors of the
  per-pattern information; something Grassmannian) and the statement of
  the rank-k sufficiency theorem, or show why no clean criterion exists.
- **C3 — a conjugate-plane obstruction principle.** Dose uniformity
  constrains the OBJECT plane (an L¹-type budget on cumulative exposure);
  conditioning anchors constrain the MEASUREMENT spectrum (spectral/trace
  bounds on the information operator); geometry selection acts nontrivially
  on the measurement spectrum while acting trivially on object-plane
  exposure — hence it survives dose but is killed by conditioning (E4),
  while per-pattern load control acts on BOTH and dies twice (E3+E4).
  Conjecture: a control helps iff it acts nontrivially on the plane where
  the information functional lives and trivially on the plane where the
  binding constraint lives — a commutation/obstruction statement. Please
  either give this a precise operator-theoretic form (which pairs of
  constraint classes and control groups commute?) or demote it to a
  heuristic with a counterexample.

## 3. The adjudication we ask of R22

- **Q1.** Rank C1/C2/C3 by (i) probability of being a real theorem,
  (ii) compression power over E1–E6 (how many existing results become
  corollaries), (iii) novelty against prior art (renewal/point-process
  information theory, Barrett–Myers task-based assessment, OED classics,
  Fellgett/multiplex analyses — cite what already exists; we would rather
  learn this is known than rediscover it).
- **Q2.** For the top candidate: sketch the actual mathematics — the
  object, the theorem statement, the proof strategy, the first
  counterexample risk. Depth over breadth; one page of real derivation
  beats ten of taxonomy.
- **Q3.** Is there a FOURTH object we have not seen that compresses more?
  (E.g. an information-geometric formulation; a single variational
  principle whose Euler–Lagrange conditions produce E2, E3, and the ridge
  simultaneously; a semigroup/channel-composition view.) One candidate
  maximum, and only if it beats C1–C3 on compression.
- **Q4.** The honest option: if E1–E6 are already at their natural depth
  and further "unification" would be numerology, say so bluntly and
  identify which ONE of the three candidates is still worth a focused
  derivation attempt as the next paper's core.
- **Q5.** Physics reach check, one paragraph each: does the chosen
  structure say anything nontrivial for (a) paralyzable/correlated-hold
  detectors, (b) sub-Poissonian (quantum) illumination meeting dead time,
  (c) detector arrays (many buckets sharing a source)? These are reach
  probes, not scope commitments.

Constraint reminder: the deliverable is UNDERSTANDING, not a campaign. No
experiment design, no endpoints, no preregistration machinery in this
round — mathematics and literature only. Ruling as a GitHub issue on
ccyyyYyzz/GI_a2 referencing this document.
