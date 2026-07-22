# Theory mainline: Hidden-State Ghost Imaging (program charter, 2026-07-22)

Operator mandate: direction = ghost imaging broadly (Chen-group problem
space, not limited to dynamic media); deep learning admissible as a TOOL,
never the claimed innovation; a lean toward high-quality reconstruction is
welcome; judging criteria = innovation, simplicity, universality, beauty.
Full autonomy under the standing operating rules.

## The program in one sentence

Every real-world GI degradation — detector dead time, hold-time jitter,
dynamic-media gain drift, speckle decorrelation, source RIN, detector
nonlinearity (SiPM cell exhaustion) — is a HIDDEN-STATE CHANNEL, and one
exact identity governs them all: information lost = the conditional
(co)variance of the hidden state given the record (the E1/Louis object;
R22's controlled conditional-score Gramian is the architecture).

## The self-similar pipeline, per channel

1. **Identity specialization** — what plays the role of "live time".
2. **Operating map / ridge** — the channel's optimal operating point and
   its law (dead-time instance: rho* ~ (2c^2+1/6nu)^(-1/3), DONE).
3. **Correction-ceiling theorem** — a Fisher/van-Trees bound on ANY
   estimator including any network: the DL yardstick ("no depth crosses
   this line; below it, we show how close one can get"). This is where
   deep learning enters as a tool racing toward a certified ceiling.
4. **Optimal design** — measure-valued design with KKT certificates
   (R25 machinery); temporal scheduling for drifting channels.

High-quality-reconstruction axis: the ceiling theorems DEFINE achievable
quality; estimators (classical RQL-TV and learned) approach it; the
range-null audit (the operator's GAN_FCC line) discloses which part of an
output is data-determined vs prior-bet — quality WITH accountability.

## Instantiation status

| channel | identity | map/ridge | ceiling | design | status |
|---|---|---|---|---|---|
| dead time + jitter | DONE (papers 1–2) | DONE | partial (E1 cap) | R23-R25 | papers in submission prep |
| gain drift (dynamic media) | R26 O1 (pending) | R26 O2 | R26 O3 | R26 O4 + GI_a1 SRHT | architecture round out |
| speckle decorrelation | open | open | open | open | after R26 |
| SiPM cell exhaustion | mine item 10 | open | open | open | parked |
| source RIN | mine item 5 | — | — | — | queued |

Chen-portfolio mapping: their DL-GI papers are empirical instances of
"network approaching a channel ceiling"; their SCGI line is the gain-drift
channel; the program supplies the theory spine their empirical line lacks.

Governance: same as ROUND63 — GPT Pro architecture/proof rounds via GitHub
issues, brutal-candor house rules, conjectures labeled, prior art
aggressive. The two ROUND63 papers and the RLMI bridge run to completion
independently and become the dead-time chapter's evidence base.
