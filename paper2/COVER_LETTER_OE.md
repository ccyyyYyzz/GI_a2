# Cover letter — Optics Express (Paper 2)

**Manuscript:** *Jitter-capped high-flux single-pixel imaging with global power control*

**Authors:** [USER: author list]
**Affiliation(s):** [USER: affiliation block]
**Corresponding author:** [USER: name, email]

---

Dear Editor,

We submit for your consideration the manuscript *Jitter-capped high-flux
single-pixel imaging with global power control* for publication in Optics
Express as a research article.

Dead time makes the flux that maximizes an integrated photon count's information
finite; this manuscript shows that hidden recovery-time jitter makes that
optimum depend on the detector itself, and turns that fact into a simple,
deployable operating policy. Starting from an exact missing-information identity
for the count-only nonparalyzable channel — uncertainty in the unobserved active
time subtracts information from the complete detector path — we derive a
long-window law in which the retained information depends on the hold-time
variance and its unique optimum obeys c_v²ρ²(1+2ρ) = 1, so the information-optimal
load is capped and falls as c_v^(−2/3). The consequence is practical: at ν = 2000
a coefficient of variation c_v = 0.05 moves the optimum from ≈22.3 to ≈5.7, and
driving the jittered detector at the deterministic ridge forfeits about half of
the available count information. We implement the resulting detector-specific
operating point with a single global source multiplier calibrated from a common
pre-scan, so the illumination geometry never changes from pattern to pattern.
Both preregistered empirical verdicts pass. In a paired simulation campaign the
policy improves fixed-dwell radiometric quality by a median +1.87 dB with 19 of
24 scenes positive; we state plainly that this power-for-time comparison spends
37.1× the incident dose and 2.6× the detected counts of its matched comparator
and is therefore not an equal-photon-budget result. Against the safe-flux
reference the median elapsed-time speedup is 19.13×, with a family-stratified
lower bound of 18.33× and 21 of 24 scenes accelerated. A separate scene-adapted
geometry search is reported as descriptive only: it locates feasible headroom
but does not resolve the scope of adaptive geometry, and we do not claim it as a
positive result.

The work is well matched to the scope of Optics Express. It concerns
computational and single-pixel imaging with photon-counting detectors, the
dead-time and recovery-time physics of real single-photon devices, and a
source-power control policy that a bench operator can apply directly through one
global multiplier. The central deliverable — an analytically derived,
out-of-sample-verified jitter cap tied to a one-knob global power policy and a
paired image-level test — targets practitioners who must decide where to run a
jittering photon-counting single-pixel system.

We have built the study to be checkable rather than merely asserted. All
simulation code, preregistration documents, per-cell manifests, correction
records, and evidence bundles are openly released under version control with
immutable freeze tags, so that each reported number traces to a fixed artifact.
The manuscript was hardened through more than a dozen rounds of external
adversarial review (R4–R21). After the campaign was unblinded, an independent
peer audit found that the frozen scoring implementation departed from the
preregistered analysis in three respects; we recomputed the affected endpoints
from the same frozen raw cells using the preregistered axis and machinery, an
independent reimplementation and an in-project rerun agreed on all 18 audited
outputs, and the corrected, spec-conformant values are the ones reported while
the original outputs remain in the provenance record. We regard this transparent
analysis-correction trail, together with the preregistration and the public
repository, as a strength of the submission, and we note that all image-level
evidence is simulation-based and is offered as motivation for, not a substitute
for, optical-bench validation.

A companion manuscript, *When high flux helps single-pixel imaging: a
contrast–dead-time operating map*, establishes the deterministic count-information
ridge that this paper extends to the jittered detector; the two are
self-contained but complementary, and we are happy for them to be considered
together or by overlapping referees.

We confirm that this manuscript is original, is not under consideration
elsewhere, and that all authors approve the submission. We suggest the following
referees as qualified and without conflict of interest: [USER: suggested
reviewers]. We thank you for your consideration.

Sincerely,

[USER: corresponding author name, on behalf of all authors]
