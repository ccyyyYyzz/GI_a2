# ROUND63 — GPT round-14: the unification round (theory architecture)

Repo: github.com/ccyyyYyzz/GI_a2. Purpose: the user has asked for a
SIMPLER, MORE BEAUTIFUL, MORE UNIVERSAL theoretical architecture for the
method line. Three candidate unifications, in decreasing order of
ambition. Please derive/refute at R11/R12 rigor, then rule which becomes
paper-2's theory spine. Fire after R13; file as issue "R14 ruling:
unification architecture".

## Candidate B — Universality class of the cube-root ridge (most exciting)

Generalize the detector: after each registered event the detector blocks
for a RANDOM hold time B ~ H with mean τ_B and cv_B = σ_B/τ_B
(nonparalyzable renewal; deterministic B = our current model, cv_B = 0).
Sketch of what we believe happens (please verify or correct):

1. Mean cycle is τ_B + 1/λ regardless of H, so the occupancy term
   E[N]/ν = ρ/(1+ρ) and the complete-data information E[N] (one unit per
   registered event from the exponential wait) are H-INDEPENDENT.
2. The missing information changes structurally: with live time
   L = T − Σ_{i≤N} B_i − (terminal), the conditional variance
   Var(λL | N) acquires an EXTENSIVE term λ²·N·σ_B² (each hidden hold
   time contributes), not just the boundary term ρ²·Var(R_ν)/ν of the
   deterministic case.
3. Consequence (conjecture): per-slot info
   J ≈ ρ/(1+ρ) − a·ρ²·cv_B²·(occupancy) − b·ρ²/(12ν), so the marginal
   cost picks up a ν-INDEPENDENT piece ∝ cv_B². Balancing against the
   saturation gain 1/ρ² gives a crossover:
   **ρ*(ν, cv_B) ≈ min{ (6ν)^{1/3}, c·cv_B^{-2/3} }** (up to constants):
   the deterministic dead time is the CRITICAL LINE (cv → 0) on which
   the ridge grows unboundedly as ν^{1/3}; ANY hold-time randomness caps
   the ridge at a jitter-limited ceiling.
4. If true, this is a universality-class statement: exponent 1/3 in ν on
   the critical line, exponent −2/3 in cv off it, constants
   model-dependent. It also makes an immediately TESTABLE prediction in
   our own simulator: the S3JIT machinery already implements
   tau_jitter_cv ∈ {0.05, 0.10} — predicted ceilings c·cv^{-2/3} ≈ small
   integers times 10, i.e. within our reachable grid. We can run a
   dev-legal jitter-ridge sweep within hours of your derivation.

**NUMERICAL CONFIRMATION (added before dispatch; scratchpad
jitter_scalar_fi.py, results/round63_study2/jitter_scalar_fi.log):**
Monte-Carlo scalar J(ρ, ν=2000, cv) for lognormal hold times, common
random numbers, 60k frames, validated against the exact table at cv=0
(J(2)=0.676 vs exact 0.667+corr):

| cv | peak location | J at peak |
|---|---|---|
| 0 | ≈22–24 (exact ridge 22.25) | 0.948 |
| 0.05 | ≈5.43 | 0.796 |
| 0.1 | ≈3.30 | 0.701 |
| 0.3 | ≤2 (grid edge) | 0.499 |

Exponent check: cap(0.05)/cap(0.1) = 1.65 vs 2^{2/3} = 1.587 (−2/3 law
within MC noise); constant fit c ≈ 0.72 from both points. Additional
finding: beyond the cap J FALLS steeply — at ρ = 24 with cv = 0.05, J =
0.38 vs deterministic 0.95, i.e. ~60% information loss from 5% jitter at
the deterministic ridge. Hardware-relevant headline.

Questions: (a) verify/correct the extensive-term derivation and the
crossover law, giving the constants for exponential and lognormal H —
your derived c should match the measured 0.72 for lognormal;
(b) is the ν^{1/3} critical exponent universal across ALL hold-time
distributions with cv = 0 (i.e., deterministic is the only cv=0 case —
so state it as: the exponent pair (1/3, −2/3) characterizes the family);
(c) does any literature (renewal information theory, queueing estimation)
anticipate this crossover?

## Candidate A — One concave program, all optima as shadow prices

Claim: every interior optimum we observed — flux ridge ρ*, occupancy
sweet spot k≈32, weight power p≈1, gain-range guard — is a KKT
stationarity condition of the single master program
max_ξ log det[V0 + ∫ J(ρ(a),ν)·g(a)g(a)ᵀ dξ(a)] over the constrained
design measure, sliced along different coordinates of the atom
parameterization. Please formalize: state the master program once, then
derive each observed optimum as its directional stationarity along
(load), (support scale), (weight power) respectively, with the
conditioning penalty appearing through the log det curvature. Rule
whether this is a THEOREM (each slice's optimum = restricted KKT) or a
narrative (in which case say what breaks).

## Candidate C — Brightness-information alignment (self-water-filling)

Empirical: under fixed flux, per-pattern load ρ_i ∝ u_i auto-allocates
load toward bright patterns, and bright patterns carry the structural
information on nonneg scenes; the uniform servo REFUTATION showed this
allocation is hard to beat. Conjecture: for nonneg scenes and nonneg
k-sparse patterns, the fixed-flux allocation achieves within a constant
factor of the water-filling optimum whenever directional value scales as
u^α with α ∈ [some range]; outside it, explicit gain control wins by
more than any constant. Please formalize the alignment condition and its
failure modes (this would explain WHY RIDGE-FIXED with one global knob
is near-optimal and delimit when OED-DT's variable load must win).

## Ruling requested

R14: (1) verdicts + derivations on B/A/C; (2) which becomes the theory
spine of paper 2 (our guess: B as the headline theorem, A as the
framework section, C as a proposition explaining RIDGE-FIXED); (3) what
of B is provable TONIGHT at the rigor of the ridge derivation vs what
needs a dedicated effort; (4) the dev-legal jitter-ridge sweep design
that would test B's crossover with the least compute.
