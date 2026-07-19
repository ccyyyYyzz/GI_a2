# ROUND63 — GPT round-10 consultation (unified method: dead-time-aware OED)

Repo: github.com/ccyyyYyzz/GI_a2. Follows R9 (issue #1). Scope: the
FOLLOW-UP method line only (paper 2 / PolyU collaboration package); nothing
here touches the OE manuscript, whose R9 scope lock stands.

## 0. New evidence since R9 (all dev-legal, frozen production path, committed)

Dev pilots on the failing families (results/round63_study2/dev_pilot_*):

- Clustered compact support "Lblob6x6" (near-full-rank translate family,
  16 px): dev microtexture S_gate 0→33.2, dev contour 0→26.3, dev chirp
  4.1→37.5; all fixed-dwell deltas positive. solid4x4 KILLS microtexture
  (safe-range-uninformative) — low-pass transfer, as spectral overlap
  predicts.
- Matched-intensity weighting (A = (n/k)·B·diag(w), w=(x̂/mean)^p from a
  4×4-blur pre-scan proxy, load-renormalized): on bandpass microtexture
  USELESS (pre-scan cannot resolve the structure; w≈const — corrected
  theory), but on its proper class (synthetic low-amplitude coarse blobs,
  pixel contrast 0.13) p=1 flips SAFE_RANGE_UNINFORMATIVE → S_gate=19.8.
  p=2 fails AGAIN despite higher design C_u (0.312 vs 0.209): C_u proxy is
  insufficient; over-weighting pays a conditioning/mismatch cost.

## 1. The unified formulation we intend to build on

Replace the hand-designed "ladder" (support shapes / weightings) by
optimal experiment design over the illumination:

  A* = argmax_{A ∈ C} Ψ( I(A; x̂) ),
  I(A; x̂) = Σ_i J(ρ_i(a_i), ν) · g_i g_iᵀ  (up to parameterization),

where J(ρ,ν) is OUR exact per-slot renewal count information (the ridge
kernel from the OE paper), g_i ∝ ∂λ_i/∂x = Φ a_i, x̂ is a coarse pre-scan
estimate (5% dose), and C = {a ≥ 0; equal mean detector load a·x̂ = const;
per-pixel cumulative dose bounds; peak irradiance bound}. Over design
MEASURES ξ on patterns, I(ξ) = ∫ J(ρ(a)) g(a)g(a)ᵀ dξ(a) is linear in ξ,
so Ψ=log det gives a concave problem → Kiefer–Wolfowitz style equivalence
certificate; Frank–Wolfe outer loop; inner atom search over a structured
dictionary (multi-scale translate families × weight powers), with the
circulant/FFT trick evaluating all translates of a shape in O(n log n).
A learned proposer may serve as an INEXACT inner oracle (computation only,
not a claimed contribution).

## 2. Questions

### Q1 — Equivalence theorem under the dead-time kernel
State precisely the KW/D-optimal equivalence condition for our setting.
Note under EQUAL-LOAD constraint ρ(a)=ρ̄ is constant on C, so J factors
out; if per-pattern load is allowed to vary, J(ρ(a)) varies but linearity
in ξ persists. (a) Any hidden failure of the classical theorem here
(noncompact atom set, rank deficiency, the nonneg cone)? (b) Correct
handling of the reconstruction being penalized-MAP (TV) rather than an
unbiased estimator — is a D-optimality criterion still the right design
objective, or should we design against a task-weighted/Bayes criterion
(see Q4)? (c) Please draft the theorem statement + proof sketch at the
rigor level of our R11/R12 collaborations.

### Q2 — Rounding bound
Continuous design ξ* → M=1024 discrete patterns: state the sharpest
standard bound (multiplicative info loss) applicable, with citation, and
any modification needed because our patterns enter row-wise into one
joint reconstruction rather than as independent replicates.

### Q3 — Inexact inner oracle
FW with multiplicative-error LMO: state the convergence degradation
theorem we can quote so that a learned proposer slots in without touching
the guarantees. Citation preferred.

### Q4 — Design objective aligned with the preregistered endpoint
Our confirmatory endpoint is Q90 time-to-quality on PSNR_rad under a TV
reconstruction, not CRLB of an unbiased estimator. Options: (i) plain
D-opt as honest surrogate + empirical validation; (ii) weighted A-opt
targeting low frequencies; (iii) Bayes-risk design with a TV-plausible
prior. Rule on which to freeze for the method campaign, and on the
mandatory GUARD you'd impose given the p=2 lesson (weight-dispersion cap?
condition-number cap? certificate-slack floor?).

### Q5 — Adaptivity-necessity theorem, OED language
Recast the R9-Q3 program: for the union-of-scales class with minimum
contrast ε at unknown scale, (a) fixed designs suffer worst-case info
M/L; (b) two-stage (O(log L) coarse + matched) attains constant fraction
of oracle; (c) ratio ~L. Please draft the cleanest statement + sketch,
flagging which steps are standard adaptive-sensing imports vs new (the
dead-time kernel weighting).

### Q6 — Prior-art narrow claim
OED for Poisson/count channels, learned sensing matrices, adaptive SPI
all exist. Proposed narrow claim: "illumination experiment design whose
information kernel is the exact dead-time renewal information J(ρ,ν),
under equal-load/equal-dose constraints, with a KW certificate and a
computable rescuability verdict (max predicted Γ < 1 ⇒ certified
don't-bother)." Verify no closer prior; adjust wording.

### Q7 — Method-campaign spec skeleton (we run experiments after your
ruling; dev-level first, confirmatory freeze later)
Please rule on: (a) fresh confirmatory scene instances (new seeds
633000+) vs reuse; (b) baselines to preregister (scattered k=16; best
fixed ladder arm; OED design) and which carries a gate; (c) endpoints
(same Q90 + fixed-dwell machinery verbatim?); (d) the pre-scan budget
accounting rule (5% off the top of BOTH arms or only the OED arm);
(e) required guards frozen into the design rule (from Q4); (f) minimum
dictionary the certificate must be quoted against.

Please rule R10 on Q1–Q7 in the operative style, with frozen wording
where needed, and file the complete ruling directly as a NEW GitHub issue
on ccyyyYyzz/GI_a2 titled "R10 ruling: dead-time OED unified method"
(the DLP blocked in-chat extraction last round; the issue IS the delivery
channel).
