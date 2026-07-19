# Cover letter — Optics Express (Paper 1)

**Manuscript:** *When high flux helps single-pixel imaging: a contrast–dead-time operating map*

**Authors:** [USER: author list]
**Affiliation(s):** [USER: affiliation block]
**Corresponding author:** [USER: name, email]

---

Dear Editor,

We submit for your consideration the manuscript *When high flux helps
single-pixel imaging: a contrast–dead-time operating map* for publication in
Optics Express as a research article.

Photon-counting single-pixel imaging is conventionally run at low flux so that
detector dead time removes few arrivals, at a direct cost in acquisition time.
This manuscript asks the operating question that convention leaves unanswered:
when does deliberately driving a count-only bucket detector into the nonlinear
dead-time regime actually improve the image, and when does it merely add
saturated photons? We answer it with an exact finite-window count-information
analysis. A convex renewal quasi-likelihood reconstruction models the integrated
nonparalyzable count without per-scene tuning, and its Fisher information has a
principal ridge whose optimal load follows a cube-root law,
ρ\*(ν) = (6ν)^(1/3) − 2/3 + O(ν^(−1/3)); the bucket-energy contrast of the
illumination patterns then determines how much of that scalar count information
reaches scene structure. We test this in two preregistered simulation studies
whose endpoints, gates, scene lists, and display subjects were frozen in
immutable commits before any confirmatory reconstruction existed, and we report
the outcomes exactly as they fell — including the formally negative primaries.
In the dense-pattern study the time-to-quality primary was uninformative and
therefore formally negative, leaving only a small, uniformly positive
fixed-dwell benefit (median +0.28 dB); in the equal-load occupancy study the
acquisition-speed primary again failed its breadth criterion while the
fixed-dwell secondary passed (median +1.16 dB), with a nonmonotone response
across occupancies. Together these results define a contrast–dead-time operating
map: high flux helps only where the illumination preserves enough bucket
modulation for the scene and the operator stays sufficiently conditioned. We
believe this honest primary-negative, secondary-positive structure is exactly
what a deployable operating map should look like, and that its value lies in
telling practitioners where the high-flux regime is worth entering rather than
in overclaiming a universal gain.

The work sits squarely within the scope of Optics Express. It concerns
computational and single-pixel imaging, photon-counting detection, and the
dead-time physics of single-photon avalanche diodes and photomultipliers, and it
connects directly to the journal's active literature on high-flux operation,
dead-time compensation, and count-rate correction in single-photon lidar and
photon-counting imaging. The deliverable is an operating map with a computable
engagement criterion, aimed at practitioners running dwell-time-limited and
photon-budget-limited single-pixel systems deep in the dead-time regime.

We have built this study to be checkable rather than merely stated. All
simulation code, the frozen preregistration documents, per-cell manifests, and
evidence bundles are openly released under version control with immutable
freeze tags, so that every reported number traces to a fixed artifact. The
manuscript was also hardened through more than a dozen rounds of external
adversarial review (R4–R21) and an independent audit of the analysis pipeline;
where that scrutiny surfaced a discrepancy we corrected it against the frozen
raw data and preserved the original outputs in the provenance record. We regard
this transparency — preregistered negatives reported as such, an auditable
analysis trail, and a public repository — as a strength of the submission, and
we note explicitly that all evidence here is simulation-only and is presented as
a quantitative map that motivates, rather than substitutes for, optical-bench
validation.

A companion manuscript, *Jitter-capped high-flux single-pixel imaging with
global power control*, extends the deterministic ridge analyzed here to the case
of hidden recovery-time jitter and its immediately deployable global power
policy; the two papers are self-contained but complementary, and we are happy
for them to be considered together or by overlapping referees.

We confirm that this manuscript is original, is not under consideration
elsewhere, and that all authors approve the submission. We suggest the following
referees as qualified and without conflict of interest: [USER: suggested
reviewers]. We thank you for your consideration.

Sincerely,

[USER: corresponding author name, on behalf of all authors]
