# Paper-2 Act III wording — R18-frozen (ruling §4; binding)

## Removed (R17 sentences that do not survive)

- R17 sentence 2 ("Under the declared dose-safety constraints, per-pattern
  load adaptation provided no certified material D-information advantage
  over the balanced design...") — DOES NOT SURVIVE AS WRITTEN; dose-only
  geometry adaptation has a constructive advantage (R18 §1.1).
- R17 final headline ("Strict uniform-dose safety collapses the tested
  adaptive design space...") — REPLACED (below).
- Every figure-caption and discussion claim that uniform dose itself
  collapses adaptive concentration — DELETED (R18 §4.5).

## Frozen four-beat Act III structure (R18 §4.2)

1. Information rewards scene-adapted geometry (dose-relaxed OED + the
   dose-only constructive probe expose substantial local D-information
   headroom).
2. Uniform dose is not enough to eliminate that headroom (a balanced
   cumulative exposure can coexist with scene-adapted measurement
   directions).
3. Conditioning safeguards define the deployable class (the A-risk and
   spectral anchors test whether the geometry headroom survives
   reconstruction-trust requirements).
4. Global operating-point control is orthogonal and survives (ridge
   tracking changes one global source multiplier; evaluated by the
   empirical power-for-time endpoints).

## Frozen paragraphs (R18 §4.3)

> A development-only constructive probe showed that the ±5% per-pixel dose
> constraint alone did not eliminate scene-adapted geometry selection: a
> feasible low-load geometry design achieved a substantially larger local
> determinant than SCAT32 in the audited cells. This was a lower-bound
> construction, not an estimate of the dose-only optimum.

> The confirmatory certificate therefore targeted the full deployed stack:
> the frozen dictionary and local pre-scan linearization, incident budget,
> ±5% per-pixel dose band, A-risk cap, and spectral floor. Global ridge
> tracking was evaluated separately as the time-limited, power-available
> control.

Conditional paragraph — permitted ONLY if FULL_STACK_CERT_PASS passes on
480/480 confirmatory cells (forbidden under the fallback branch):

> Within the frozen full-stack class at the preregistered anchor cells,
> balanced SCAT32 was certified within 1% local D-efficiency of the best
> admissible design. The useful deployed positive control was the global
> detector operating point; geometry optimization under relaxed
> conditioning remains a separate opportunity.

If the certificate fails or is removed pre-freeze (R18 §5.4 fallback):
report the distribution of CERTIFIED / COUNTEREXAMPLE /
NUMERICAL_UNRESOLVED cells without any categorical collapse claim.

## Frozen final headline (R18 §4.4)

> Uniform dose controls cumulative exposure but does not eliminate
> scene-adapted geometry selection. The operative restriction is the full
> deployed conditioning stack: balanced SCAT32 is tested for certified
> near-optimality within that frozen class, while time-limited image gains
> are tested separately by global ridge tracking.

After a successful certificate only, "is tested for" may be replaced by
"is certified within 1% local D-efficiency." No stronger substitution is
permitted.

## Figure/handbook edits (R18 §4.5)

- Act III (a): retain the path-specific guard-collapse DEV plot.
- Add/replace one DEV panel with the constructive dose-only primal lower
  bound (results/round63_m1/R18_GAP_PROBE_REPLICATION.md), explicitly
  marked development-only.
- The confirmatory certificate panel is labeled FULL_STACK_CERT_PASS.
- The handbook says "full deployed conditioning class," not "expanded
  dose-safe class."
