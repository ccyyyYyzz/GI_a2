# R30 — Prior-art fence + campaign design for CERTIFIED INFORMATION-GUIDED IMAGING

The direction assembled after R29's saturation conclusion. All evidence
frozen in-repo. This round: ruthless fence, headline-claim sharpening,
DEV-gate design. The operator's bar is unchanged
(创新 × 简洁 × 普适 × 美观 + a genuine visible improvement, software only).

## The assembled package (each part with measured evidence)

**Thesis:** in compressed-sampling photon-counting imaging, ship every
reconstruction with an exact, physics-grounded, per-direction split of
"what the data determined" vs "what the prior invented" — and use that
split (a) to make the reconstruction ITSELF better where it matters
(anomaly/atypical-content fidelity: the prior may speak only where data
is silent), and (b) as a GT-free third-party auditable certificate.

Components:

1. **Chassis (exists, operator's GAN_FCC line — audited today,
   results/round63_next/GAN_LINE_AUDIT/)**: exact range/null gauge
   geometry; generator structurally confined to the null space; exact
   dual projector (~1e-8 in 0.6 s/img); truth-blind attainable-record
   closure. Known state: hard NOISELESS Euclidean consistency (needs the
   Fisher-weighted upgrade); FOHI orthogonalization is a measured no-op
   (discarded); measurement certificate in the frozen eval used a
   stalling Dykstra (the working exact-dual solver simply wasn't wired).
2. **Physics ledger (this program, proven)**: E1/O1 missing-information
   identities for real counting channels (dead time, jitter, gain drift)
   → the consistency weights and the certificate's actuarial table are
   EXACT, not model-trusting heuristics.
3. **Ledger-completeness certificates (this program, measured)**: five
   saturation certificates (docs/SOFTWARE_SATURATION_VERDICT.md) — at
   the deployed points the data ledger is provably fully spent (RQL
   ~99% efficient), so "beyond the ledger = prior" is a measured fact.
4. **Data term with real dB (measured today, LADDER_PROBE)**: for drift
   channels, the gain-marginalized temporal-prior likelihood recovers
   **+2.2–2.5 dB over published-style corrections on the same bucket
   record**. Adjudicated: admissible as the data term (it is your R29
   JGM lane — prior-art-heavy alone), novelty carried by the
   integration, not the smoother.
5. **Anomaly-fidelity feasibility (measured today, GAN_LINE_AUDIT Q4)**:
   at 5%/10% sampling, localized anomalies land 0.43/0.60 of their
   energy in the measured subspace (≥3px: 30–96%; sub-3px: 65–87% null).
   So: fidelity route viable for resolvable anomalies; attribution/
   flagging route for sub-resolution ones; the split is an
   operator-design knob we control.
6. **Audit protocol (operator's prior theory line)**: no-free-audit
   theorem, impossibility quadrilateral, GT-free post-commitment
   challenge (catches the barrier witness 100% at k≥4). Internal
   adjudication already on record: the structural inequalities are
   prior art (Backus–Gilbert lineage); the novel slice is the
   OPERATIONALIZATION (GT-free challenge + witness + attribution).

**The improvement claim to be tested (DEV gate):** unconstrained DL
priors improve average PSNR by betting typicality, and thereby ERASE
atypical content (documented fastMRI-style failure); information-guided
reconstruction (data-subspace locked with Fisher weights, prior confined
to null space, attribution shipped) wins anomaly-region fidelity
decisively while conceding little background quality — plus flags
sub-resolution anomalies as unwitnessed instead of silently inventing
or deleting them.

## R30 asks

1. **Fence, component by component, against the sharpest prior art you
   can name** (with DOIs): (a) information/noise-WEIGHTED data
   consistency vs hard data-consistency layers in DL-MRI/unrolled nets;
   (b) per-direction data-vs-prior attribution certificates vs conformal
   prediction for imaging, Bayesian deep-imaging UQ, posterior-sampling
   uncertainty maps; (c) adversarial anomaly-fidelity evaluation vs
   fastMRI+ pathology metrics and hallucination-measurement literature;
   (d) exact counting-channel ledger (dead time/jitter/drift Fisher
   accounting) inside a reconstruction+certificate system — we believe
   THIS is the unclaimed square; check it; (e) GT-free post-commitment
   audit challenges in imaging. For each: what survives as narrow
   novelty, what must be cited defensively, what must NOT be claimed.
2. **The headline**: state the sharpest honest one-sentence claim of the
   package (the "eye-lighting" formulation), and the claim most likely
   to get the paper killed if we overreach.
3. **DEV gate design** for the anomaly-fidelity probe now running
   (compressed GI on the operator's real lane0 rate05 operator, arms:
   data-consistent baseline / unconstrained DL prior / information-
   guided+certificate): frozen bars for "improvement is real"
   (anomaly-region fidelity margin, background-quality concession limit,
   flagging correctness for sub-resolution anomalies), in the style of
   the M1/bridge gates.
4. **Scope ruling**: one flagship system paper (operator's standing
   one-good-paper order) — which pieces of the ROUND63 dead-time line
   (identity, ridge law, M1 verdicts, saturation certificates) belong
   INSIDE this paper as load-bearing sections vs get cited as companion
   evidence? The operator's two prior manuscript lines (GI_a1
   identifiability; FOHI draft as quarry) exist; how do they relate
   without violating the one-paper order?
5. If you believe the package still fails the operator's bar, say so
   plainly and say why.

Deliver as a GitHub issue on ccyyyYyzz/GI_a2. Depth over speed.
