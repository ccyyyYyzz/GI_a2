# ROUND63 — GPT Round 23 (build the bridge: the rank-constrained control theorem)

Charge: R22 (issue #14) identified exactly one derivation still worth making:

> A rank-constrained control theorem for hidden-state photon-counting
> experiments: exact/full KKT sufficiency, computable projected-residual
> upper bounds, and lower bounds showing when rank-one global power control
> cannot approximate unrestricted load/geometry control under nonnegative
> illumination and physical dose/peak constraints.

The operator's instruction, verbatim in spirit: "build the bridge." This
round asks for the MATHEMATICS — theorem statements and proofs at
manuscript rigor where achievable, clearly-labeled conjectures where not
(R11/R22 house style). No campaign design, no endpoints. Papers remain
frozen; this is the core of the next work.

## 1. The concrete setting (fix notation once, use ours)

- Scene x ∈ R^n_{≥0}; patterns a_i ∈ R^n_{≥0} (nonnegative illumination);
  per-pattern bucket rate λ_i = Φ_i a_i^T x with per-pattern source gains
  Φ_i ≥ 0 (the control); load ρ_i = τλ_i.
- Dead-time kernel: per-pattern count information J(ρ_i) with the frozen
  long-window form J(ρ) = ρ/[(1+ρ)(1+c²ρ²)] (declared monotone-concave
  branch ρ ≤ ρ̄_br as needed — state the branch hypothesis explicitly,
  per R22 counterexample-risk 2).
- Information matrix V(Φ) = V_0 + Σ_i J(ρ_i(Φ_i)) q_i q_i^T on the task
  subspace (q_i the information directions; ℓ_i(V) = q_i^T V^{-1} q_i the
  leverages); objective f(Φ) = log det V(Φ).
- Constraint cone C from the campaign: Φ ≥ 0; incident budget
  Σ_i ρ_i ≤ B; per-pixel dose band |d_j/d̄ − 1| ≤ δ with
  d_j(Φ) = Σ_i Φ_i a_ij; peak caps. (A-risk/spectral MAY be added as a
  variant; start with budget+dose+peak.)
- Rank-k control manifolds: M_1 = {Φ_i = κ u_i : κ ≥ 0} (the global
  multiplier over a fixed base u); M_k = span of k fixed gain profiles;
  and the "geometry" variant where the pattern SET is also selectable from
  a dictionary (connect to the R18 C_stack machinery if useful).

## 2. Deliverables (in priority order)

- **T1 (KKT sufficiency + gap upper bound).** State and prove: (i) the
  first-order condition characterizing the rank-k optimum on C ∩ M_k;
  (ii) the projected-residual ε_k = dist_*(∇f(Φ_k), N_C(Φ_k)) and the gap
  bound f(Φ*) − f(Φ_k) ≤ ε_k²/(2m) with an EXPLICIT, checkable strong-
  concavity modulus m for log det on our feasible set (log det is not
  globally strongly concave — derive m from spectrum bounds
  m ≥ 1/λ_max(V)² on the segment, or the correct local statement; handle
  the J-branch hypothesis). Make ε_1 fully computable for M_1: the scalar
  condition Σ_i ℓ_i J'(ρ_i) ρ_i = (prices), and the residual vector in
  closed form so it can be evaluated per scene from campaign quantities.
- **T2 (the impossibility side — the actual bridge).** Prove a lower
  bound: construct (or characterize) scene/task families where rank-one
  control leaves an Ω(1) gap. Target form: if brightness u and the
  leverage-weighted information profile are misaligned in a precise sense
  (e.g. the vectors {ℓ_i J'(ρ_i)ρ_i} and the achievable ascent directions
  under M_1 span an angle bounded away from the full cone's), then
  f(Φ*) − f(Φ_1) ≥ c(misalignment) > 0, with c explicit. A two-population
  scene (bright/uninformative vs dim/informative patterns) is the natural
  witness — our contour family is the empirical instance. State exactly
  which functional of (u, {q_i}, V_0) controls the gap: this is the
  corrected "alignment functional" replacing C_u (R22 §2).
- **T3 (rank decay).** How does the optimal gap decay in k? Conjecture or
  theorem: for scenes whose leverage profile has r distinct "clusters,"
  rank-r control achieves ε_r = 0 (or ≤ tolerance); relate r to a
  numerical rank of an explicit matrix (e.g. the moment matrix of
  {(u_i, ℓ_i J'(ρ_i))}). Even a clean two-cluster theorem is valuable.
- **T4 (the visibility-proxy corollary — paper-1 bridge).** Derive the
  first-order expansion in which C_u appears as a proxy of the true
  residual: under stated homogeneity assumptions (e.g. q_i alignment
  homogeneous across patterns), show ε_1 reduces to a monotone function of
  bucket contrast — recovering the operating map's contrast axis as the
  aligned special case, and quantifying when the proxy breaks (the R22
  identical-C_u counterexample as the boundary).
- **T5 (scope fences).** One paragraph each: what survives on the
  nonconcave J branch (beyond the ridge); what changes when the pattern
  dictionary (geometry) is part of the control (connect to the R18
  support-function certificate: is ε_∞ exactly the G_stack object?); and
  the explicit statement that image-domain (PSNR/TV) claims are NOT
  corollaries (R22 risk 3) so the theorem's empirical companion stays
  honestly scoped.

## 3. House rules

Full proofs where achievable; label anything less as CONJECTURE with the
obstacle named. Cite prior art aggressively (Kiefer–Wolfowitz equivalence,
restricted-design/marginal designs, Pukelsheim, low-rank control in OED if
it exists — if T2 is already in the literature we want the citation, not a
rediscovery). LaTeX-ready statements welcome. Ruling as a GitHub issue on
ccyyyYyzz/GI_a2 referencing this document.
