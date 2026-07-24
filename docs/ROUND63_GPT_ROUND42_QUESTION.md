# R42 — Pure innovation round at PRL grade: find the right deep object

This round has NO adjudication asks. It is a maximal-depth creation round.
The operator's instruction: aim at PRL deliberately — the demand itself
raises the ceiling. Take your maximal deliberation. Derive before you
propose. Depth over coverage.

## A. Where the program stands (honest, compressed)

PROVEN (all corrected-shot, MC-validated where empirical): the exact
first-moment wall (any band-limited code set, any schedule: zero response
beyond B_p; operationally verified to 1e-16 and used as a specificity
discriminator); the covariance aperture law (support = B_p ⊕ B_w, 18.76×
in/out at 64×64); the demand-threshold dissociation (the SAME efficient
Fisher spectrum yields dead imaging — mode-coverage f_rec 0.033 vs 0.70 —
and live detection — 2% anomaly at d'≥5 in 95% of an 81-cell physical
grid, 467 banks = 6.0 s best; MC ROC 0.98–1.04 of analytic); the
streaming CUSUM sentinel (latency 1.1/5.4/17.5/54 s for 5/2/1/0.5% at
1 FA/hr; Siegmund-validated; ~1/ε² law); the 0.5% "floor" broken (it was
integration time: ε_min 0.37–0.41% at 16384 banks); four-way attribution
(in-band scene / beyond-band scene / medium-amplitude separate on
orthogonal scores, d' 61/13.5/61; tau class pending in the sealed
probe); trace robustness (10% law rotation → AUC 1.000→1.000; exact
orthogonal-similarity lemma with stated failure boundary); blind
CRB-efficiency (the moment pipeline attains the corrected CRB 0.98–0.99×;
the convex lift was built, validated convex, and did not beat it — the
E8 plateau WAS the finite-T CRB); M1 (certified +1.87 dB / 19.13×
mean-channel operating optimum).

DEAD (each at a frozen bar or decisive test): beyond-band IMAGING at
flagship grade (REGION_EMPTY across the physical space; natural pocket
eliminated at f_rec level after the shot correction; the CRB-efficient
1.25× pocket needs ~26 s for NRMSE 0.30 — calibration point only);
third-order cumulants (starved, σ_f^2k scaling); detection-optimal code
design (inert, ~3% — 128 band-limited codes already span the ~121-dim
band); pathwise blind estimation (global collusion, marginalization
mandated); AND — post-mortem you must absorb — YOUR R41 Rank-1
spectral-duality story: P=[tr(J)/d][tr(J⁺)/d] ≥ 1 holds exactly 81/81
(range 1.05–2.64), but the qualitative claim "spiky spectrum = strong
detector" is WRONG-SIGNED on the physical grid (corr(P, T_det) = +0.505):
anisotropy CO-VARIES with total supply across physically-reachable
media, so spiky cells are worse at both tasks; "spiky helps detection"
requires fixed trace, which physics never grants. The genuinely
non-circular collapse (f_rec vs lower-tail mass) failed (R²=0.085). The
duality is demoted to a lemma. The REAL mechanism of the dissociation is
demand thresholding: imaging needs each mode above SNR γ; detection
pools the trace. Learn from this: the right deep object must survive
the covariation structure of PHYSICAL parameter space, not just algebra.

JUST LANDED (results/round63_next/SCRAMBLE_EXT/, derivation + torch sim):
the FULLY-SCRAMBLING extension. Under developed speckle
(zero-mean circular-Gaussian field, ⟨TT*⟩ = C(r−r')·δ), three exact
results, all numerically confirmed: (i) the mean wall becomes TOTAL —
E[b] sees only the scene DC (band {0}; null verified 6.6e-19);
(ii) covariance support pierces the entire spectrum (1024/1024 bins) but
COLLAPSES TO RANK ONE — the scene enters only through Q = xᵀGx (99.99%
of covariance energy in one direction) → beyond-band imaging does NOT
survive full scrambling; (iii) the SENTINEL survives and upgrades:
ΔQ = 2xᵀGδ + δᵀGδ with δᵀGδ > 0 for ANY nonzero δ — a sign-definite
response with NO blind direction, while the mean stays exactly blind
(5% coherent change: d'=4.1, AUC 0.997 in 8192 banks ≈ 1.7 min; d'∝√T
throughout). Currency: independent banks (M_eff ≈ 13 ≈ in-band code
DOF). VERDICT: detection is UNIVERSAL across the medium family — from
thin fog (quantitative 81-cell maps) to total scrambling (exact
structural results), the wall only strengthens and the sentinel never
goes blind; imaging degrades from thin-rim to rank-one. The
memory-effect interpolation (rank vs r_ME) is the flagged connecting
family between the two anchored endpoints.

IN FLIGHT: the sealed detection probe (D0–D7, fresh-covariance
full-strength comparator, tau attribution); the memory-effect
interpolation probe.

## B. What PRL-shaped means here (the target definition)

One broad physical principle; exact or near-exact; surprising to a
general physics reader (not only to imaging specialists); instantiated
cleanly in our bucket/coded-illumination system but NOT confined to it;
carrying one decisive falsifiable prediction that our engines can test
this week. Named-effect potential is a plus (PRL rewards effects with
names). The current best candidate sentence — "an exact optical
blindness theorem becomes the specificity mechanism of a software-only
covariance sentinel" — is Optica-shaped. Your job: find what is
PRL-shaped above or beneath it.

## C. Mandatory thinking protocol (do these moves, in order, showing work)

1. **The invariant hunt.** Before proposing anything: derive. Search the
   record's structure for the deepest unnamed invariant / conservation /
   exact tradeoff. Candidate veins you must examine and then go beyond:
   (a) a fluctuation–dissipation-type relation — our transducer IS a
   fluctuation; is there an exact FDT-like identity linking the medium's
   fluctuation spectrum to the information current through the bucket
   (information gained per unit fluctuation power per unit time), with
   the mean channel as the "dissipationless" (information-free) part?
   (b) a Stein/Chernoff-exponent formulation — detection latency laws as
   large-deviation exponents of the record; is there an exact exponent
   duality between the two moment layers? (c) an uncertainty-type
   inequality between what a record can localize (image) and what it can
   notice (detect) — with our data as its first instantiation; (d) a
   universality statement — which of our laws survive ANY stationary
   medium law and ANY code bank (the aperture ⊕ law looks
   universality-class-shaped)?
2. **The assumption inversion sweep.** List ≥8 deep unquestioned
   assumptions of statistical optics / computational imaging (e.g.,
   "information about a scene enters through the mean field", "temporal
   averaging improves measurement", "resolution is a property of the
   instrument", "noise floors are floors"). For each: state what exactly
   breaks if inverted, and whether our machinery can realize the
   inversion.
3. **The cross-field transplant sweep.** For each of: quantum metrology
   (SPADE/Rayleigh-curse is our sibling — what is OUR curse-and-cure
   statement?), stochastic thermodynamics (Landauer/Sagawa-Ueda
   information-heat bounds), disordered-systems statistical mechanics
   (self-averaging, universality), random-matrix theory (our aperture
   edge), sequential analysis (SPRT optimality), control/estimation
   duality — name the ONE theorem there whose transplant into the bucket
   record would yield an exact new statement, and sketch the transplant.
4. **The named-effect audit.** Among our verified phenomena — the
   dissociation (detect what cannot be imaged), the deep-fog advantage
   (contrast helps detection, inert for estimation), the mean-wall
   specificity, the 1/ε² latency law, anisotropy-blindness of detection
   — which can be elevated to a general named effect with a universal
   scaling law and a one-line definition? Draft the definition.
5. **Synthesis.** ≥7 ranked ideas (surprise × rigor × reach, each
   scored), each with: the one-line mechanism; the exact-statement
   sketch (assumptions → claim); the decisive kill test runnable on our
   engines this week; your best guess of the nearest prior art and why
   it does not own the statement; and the PRL-shape verdict (what a
   general physicist finds surprising).
6. **The moonshot derivation.** Pick your single deepest candidate and
   DERIVE it in the issue — actual mathematics, not a proposal: state
   assumptions, derive the identity/inequality, exhibit the equality
   conditions, and give the numerical experiment that would falsify it.
   If the derivation fails midway, show where and say what is missing.

## D. Constraints

Honesty machinery applies to you symmetrically (your Rank-1 died at its
own decisive test — that is the system working). Ideas may break our
current frame; the aperture atlas is not sacred. Sim-only until the
manuscript is done. The bucket/coded-illumination system is the home
instantiation; generality is the prize. Do not repeat the nine R41
ideas except where a new derivation materially upgrades one.

Deliver as a GitHub issue titled R42. Take maximal thinking time.
