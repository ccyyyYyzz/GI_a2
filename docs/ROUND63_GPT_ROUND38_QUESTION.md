# R38 — Fluctuating-medium measurement diversity ("fog as a second DMD"): identifiability theorem + Fisher accounting + ruling

Housekeeping: R37 (v2 instrument-fix prereg, @ 4ee707e) was committed but
never sent; it is SUPERSEDED — the width prognosis it requested was run
locally first and returned 0/27 (calibration-hull domain) / 2/27 (honest
wide domain), with narrow intervals achievable only via assumption
clamping (bounded fraction 0.00 everywhere) and the gauge fix WORSENING
the C4 task metric in all 3 failed cells. DLGI v2 is dead at prognosis;
no campaign was launched; the v1 FLAGSHIP_KILL stands. Do not answer R37.

A hostile whole-archive re-audit (38 findings → 7 candidates → 0
survivors under adversarial verification) then confirmed: no honest,
novel, universal, ≥1 dB software gain exists on the FIXED-operator
bucket record. The operator rejected the audit/accountability line as a
deliverable and demanded divergent rethinking. R38 is that rethink. It
attacks the premise every prior wall shares: A IS FIXED.

## The proposal

Computational GI/SPI: known pattern set {P_i}, i=1..M (DMD, reused every
epoch), STATIC scene x (N pixels), and an UNKNOWN spatially-structured,
temporally-fluctuating medium between source and bucket. Per coherence
epoch t = 1..T the effective measurement row is P_i ⊙ w_t (elementwise),
with w_t the medium's spatial transmission field:

    B_{i,t} = <P_i ⊙ w_t, x> (+ shot noise).

Claim under test: the medium fluctuation is a FREE measurement-diversity
generator. The stacked operator {P_i ⊙ w_t} escapes the fixed row space;
scene content in the fixed system's null space (A P_0 = 0 — provably
invisible forever at fixed optics) becomes recoverable by JOINT
scene+medium estimation, with NO reference arm, NO wavefront sensor, NO
hardware change. The scalar-gain model of the entire prior program
(a_t · P_i — row direction unchanged) is precisely the unique degenerate
fluctuation that generates zero diversity: our 27-cell world was the
worst case, and the walls we proved there do not bind here because the
RECORD itself is richer.

Inversion pitch (the operator's aesthetic bar): the field's celebrated
"turbulence-free GI" proves correlation GI is IMMUNE to turbulence —
i.e. it discards exactly this information by design. We claim the
turbulence is the resource. Biological echo: fixational eye movements —
evolution jitters the operator on purpose.

## Oxygen test (local, clean linear, dev-grade, done tonight)

N=1024 (32×32), M=64 known Gaussian patterns, T=64 epochs, w_t = 1 + U
z_t with U an orthonormal smooth spatial basis (d_w columns), z_t ~ N(0,
σ_w²), σ_w = 0.30, scene with strong fixed-A null content.

| case | stacked rank | oracle (w known) null-err | ALS warm (truth+15%) | ALS cold |
|---|---|---|---|---|
| scalar gain (d_w=1) | 64 (row space unchanged) | 1.000 (nothing) | — | — |
| spatial d_w=16 | **1024/1024 full** | **0.033** | 0.558 | 1.000 (fail) |
| spatial d_w=64 | full | **0.000** | 0.133 | 0.991 (fail) |
| spatial d_w=256 | full | **0.000** | 0.100 | 1.002 (fail) |

Reading: (i) the record CONTAINS the full scene including 100% of the
fixed-system null space — geometry confirmed; (ii) blind extraction is
the obstacle: cold-start alternating LS fails (classic bilinear
hardness), warm-start largely succeeds; (iii) the naive counting
"M > d_w" phase transition is WRONG in detail — d_w=256 makes the blind
problem globally underdetermined (N + T·d_w > M·T) yet warm-init sticks
near truth; the correct identifiability statement is what we ask you
for. Note the medium is temporally correlated (OU-like) in reality, so
only epoch 1 is truly cold: tracking (the DLGI medium-estimation
machinery, whose point estimates passed 26/27 cells at 1.0–1.5× oracle)
is the natural warm-start ladder.

## R38 asks

1. **Identifiability theorem.** The correct generic-identifiability
   condition for joint (x, {z_t}) recovery from {B_{i,t}}: counting
   (M·T ≥ N + T·d_w − gauges?), the gauge group (global scale; static
   smooth envelope x·h vs w_t/h confounding — state exactly what is
   recoverable: x up to a d_w-smooth multiplicative envelope?),
   genericity conditions on {P_i} and U, and the role of temporal
   correlation of z_t (does an OU prior on z restore identifiability
   where per-epoch counting fails?). State the honest phase boundary in
   (M, d_w, T, N).
2. **Fisher/noise accounting.** Under Poisson shot noise with total
   photon budget PH split over T epochs: the Fisher information of the
   stacked fluctuating record for the null-component of x. Does the
   extra rank survive realistic budgets, or is the smallest singular
   value of the effective operator so small that the "revealed" modes
   carry sub-noise information? Give the scaling of recoverable-mode
   count vs PH, σ_w, d_w, T — the analogue of the width-floor analysis
   that killed DLGI certification, run BEFORE we fall in love. Include
   the multiplicative-noise penalty: w_t must itself be estimated, so
   the Schur/profiling tax applies; when does joint estimation beat
   simply averaging the medium out (classic GI)?
3. **Optimal exploitation.** Given the physics (OU-correlated z_t,
   known correlation time), the right estimation architecture: tracking
   EM / spectral init / convex lifting (Ahmed–Recht–Romberg style) —
   and whether a DELIBERATE weak modulation of the source (still
   software: pattern reordering across epochs, no new hardware) can
   convexify the cold start.
4. **Prior-art fence.** Closest families we already see: imaging with
   nature (Liutkus 2014, static medium + calibration), blind
   ptychography, self-calibration/bilinear CS (Gribonval; Ling &
   Strohmer), speckle-correlation imaging (Bertolotti 2012, Katz 2014),
   turbulence-free GI (Meyers), pseudothermal GI itself (medium IS the
   pattern source but measured via reference). A hostile scout is
   running in parallel; from your knowledge: name the killing paper if
   one exists, and if not, state the defensible novelty square exactly
   ("no reference arm + known codes + unknown dynamic medium + joint
   estimation + identifiability/Fisher theorem + turbulence-as-resource
   framing"?).
5. **Blunt judgment against the operator's bar** (one surprising
   inversion, dB-scale or made-possible-at-all effects, robust,
   universal, software-only, honest): eye-lighter or boring? If the
   Fisher accounting in (2) predicts the revealed modes are
   noise-drowned at realistic budgets, say KILL now and save the
   program the campaign.
6. **If GO: the preregistered probe spec.** Frozen bars for a decisive
   feasibility probe (sim first, bench-transferable): what must be
   demonstrated (null-mode recovery at what SNR margin vs the
   fixed-system impossibility baseline; blind, not oracle; cold-start
   included), sealed-bank structure, and the kill tree — in the R34
   style that has served this program well.

Deliver as a GitHub issue titled R38. The machinery is warm; the
operator is waiting.
