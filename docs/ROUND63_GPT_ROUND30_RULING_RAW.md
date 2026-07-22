# R30 ruling (GitHub issue #22, raw)

Title: R30 ruling: certified information-ledger imaging fence and DEV gates
Posted: 2026-07-22T12:00:53Z

---

# R30 — Prior-art fence and DEV ruling for certified information-guided imaging

**Reference commit:** `d5170d5`  
**Source:** `docs/ROUND63_GPT_ROUND30_QUESTION.md`  
**Binding program record:** `docs/SOFTWARE_SATURATION_VERDICT.md`

**Scope:** software-only reconstruction and audit. The frozen ROUND63 campaigns and their verdicts remain unchanged. No manuscript is authorized by this ruling; the running anomaly-fidelity DEV probe must first pass the gates below.

> **R30 VERDICT: CONDITIONAL GO TO THE DEV DECISION, NOT YET A PAPER GO.**
>
> The package contains a potentially publishable method, but only in a narrow form. None of its pieces is individually new: weighted data consistency, null-space networks, likelihood-informed subspaces, Fisher-based prior-dominance maps, Bayesian/conformal uncertainty, pathology-aware reconstruction tests, and held-out measurement checks all have strong prior art. The residual contribution that may survive is the **system-level coupling** of:
>
> 1. a detector-channel observed-information ledger that includes dead time, hidden hold jitter, and marginalized gain drift;
> 2. a reconstruction architecture that locks the likelihood-supported component and confines the learned prior to the declared information-null complement;
> 3. an explicit per-direction data-versus-prior provenance output; and
> 4. a post-commitment, fresh-measurement challenge that audits witnessable claims without deployment ground truth.
>
> Even that combination is not enough by itself. The running DEV test must show a visible anomaly-fidelity gain over **both** an unconstrained generative prior and the strongest ordinary Euclidean range/null method. If it only produces a nicer heat map, or if its gain comes entirely from the already-known gain-marginalized likelihood, the method does not clear the operator’s novelty × simplicity × generality × image-positive bar.

---

# Executive ruling

| item | ruling |
|---|---|
| **Information-weighted consistency** | Established in broad form. The narrow residual is the channel-specific observed-information metric and its use in a structurally separated data/prior reconstruction. |
| **Range/null confinement** | Established. Deep null-space learning and GAN null-space constructions are direct precedents. It cannot carry the novelty claim. |
| **Data-versus-prior attribution** | Very crowded. Bhadra et al. already define generalized measurement/null hallucination maps; Sampson–Melchior already compare likelihood and prior Fisher information to flag prior-dominated regions. The surviving distinction must be exact architectural provenance plus a real counting-channel ledger, not a generic uncertainty map. |
| **Anomaly-fidelity evaluation** | Established as a need and as an evaluation paradigm. The residual is the operator-resolved/null-energy stratification and the paired anomaly-response endpoint. |
| **Counting-channel ledger inside the system** | The strongest narrow novelty candidate. I found no prior system combining dead-time/jitter/gain-drift observed information with null-confined generation, per-direction provenance, and a post-commitment audit. Phrase as “we formulate and evaluate,” never as an unsupported priority claim. |
| **GT-free post-commitment challenge** | Potentially distinctive operationalization. Held-out measurements, self-supervised measurement splitting, posterior checks, and data consistency are old. The novelty is commit-first, fresh randomized challenge acquisition for a declared reconstruction claim. It verifies witnessability/consistency, not truth. |
| **Current status** | **Not yet above the bar.** Await the frozen DEV gates. Failure against the Euclidean null-space comparator is a hard no-go for the flagship method. |

---

# 1. Component-by-component prior-art fence

## 1.1 Information/noise-weighted data consistency

### Closest prior art

Physics-guided and unrolled networks have enforced hard or soft measurement consistency for years:

- Schlemper et al., deep cascaded MRI reconstruction, DOI [10.1109/TMI.2017.2760978](https://doi.org/10.1109/TMI.2017.2760978).
- Hammernik et al., learned variational network with learned data-term weights, DOI [10.1002/mrm.26977](https://doi.org/10.1002/mrm.26977).
- Aggarwal, Mani and Jacob, MoDL with numerical data-consistency blocks, DOI [10.1109/TMI.2018.2865356](https://doi.org/10.1109/TMI.2018.2865356).
- Hammernik et al., systematic study of data-consistency layers, DOI [10.1002/mrm.28827](https://doi.org/10.1002/mrm.28827).

Structural null-space restriction is also established:

- Schwab, Antholzer and Haltmeier, deep null-space learning, DOI [10.1088/1361-6420/aaf14a](https://doi.org/10.1088/1361-6420/aaf14a).
- Wang et al., GAN-prior null-space learning for consistent super-resolution, DOI [10.1609/aaai.v37i3.25372](https://doi.org/10.1609/aaai.v37i3.25372).
- Quan et al., cooperative range/null learning from incomplete measurements, DOI [10.1109/TPAMI.2024.3359087](https://doi.org/10.1109/TPAMI.2024.3359087).

Weighted likelihoods, generalized least squares, Poisson likelihoods, and Fisher scoring are classical. The gain-marginalized state-space likelihood is effective, but by the frozen adjudication it is a data term, not the novelty.

### What survives

The defensible residual is:

> a **channel-derived observed-information lock** for photon-counting bucket data, where the metric is computed from the declared dead-time/jitter/gain-drift likelihood and the learned prior is structurally prevented from altering the declared data-supported component.

This is narrower than “physics-guided reconstruction” and narrower than “Fisher-weighted data consistency.”

### Critical mathematical qualification

For a linear forward operator `A` and strictly positive row weights `W`,

\[
\ker(W^{1/2}A)=\ker(A),
\qquad
\operatorname{range}(A^TW^{1/2})=\operatorname{range}(A^T).
\]

Thus Fisher weighting does **not** create a new algebraic null space. It changes the metric, the stable right inverse, and the relative strength of measured directions. This is precisely why the Euclidean range/null chassis is a mandatory comparator: if the proposed method does not beat it on anomaly fidelity, the counting physics improved bookkeeping but did not create an effective reconstruction method.

### Must not be claimed

Do not claim:

- first data-consistent learned reconstruction;
- first null-space generative prior;
- first information-weighted inverse solver;
- that Fisher weighting turns weakly observed directions into known truth;
- that the projector is globally exact for a nonlinear, plug-in gain-marginalized model.

Use the term:

> **channel-exact under the declared counting model, and local/plug-in at the reconstructed scene.**

If the gain posterior covariance is numerically approximated, state that approximation separately; do not label the entire deployed ledger exact.

---

## 1.2 Per-direction data-versus-prior attribution

### Closest prior art

This is the most crowded component.

- Cui et al. identify likelihood-informed and prior-dominated subspaces in Bayesian inverse problems, DOI [10.1088/0266-5611/30/11/114015](https://doi.org/10.1088/0266-5611/30/11/114015).
- Spantini et al. derive optimal low-rank likelihood-informed approximations, DOI [10.1137/140977308](https://doi.org/10.1137/140977308), and goal-oriented variants, DOI [10.1137/16M1082123](https://doi.org/10.1137/16M1082123).
- Bhadra et al. decompose estimates into generalized measurement and null components and introduce hallucination maps; their measurement-space map is ground-truth-free, while null-space hallucination assessment requires truth, DOI [10.1109/TMI.2021.3077857](https://doi.org/10.1109/TMI.2021.3077857).
- Repetti et al. test whether particular structures are supported under a Bayesian imaging model (BUQO), DOI [10.1137/18M1173629](https://doi.org/10.1137/18M1173629).
- Narnhofer et al. combine posterior variance with conformal prediction to obtain pixelwise error bounds, DOI [10.1137/23M1546129](https://doi.org/10.1137/23M1546129).
- Luo et al. perform posterior MRI reconstruction with uncertainty maps, DOI [10.1002/mrm.29624](https://doi.org/10.1002/mrm.29624).
- Angermann, Göppel and Haltmeier combine data-consistent null-space networks with learned uncertainty maps, arXiv [2304.06955](https://arxiv.org/abs/2304.06955).
- Most dangerously, Sampson and Melchior explicitly compare Fisher information from the likelihood and the learned prior to flag prior-dominated image regions, arXiv [2306.13272](https://arxiv.org/abs/2306.13272).

The last paper is a direct collision with any broad “Fisher data-versus-prior attribution” claim.

### What survives

The package must distinguish **provenance** from **uncertainty**:

1. the output is architecturally decomposed,
   \[
   \hat x=\hat x_{\rm data}+\hat x_{\rm prior};
   \]
2. the prior component satisfies a numerically certified channel-null constraint;
3. every direction is accompanied by the declared detector information and by the amount contributed by each component;
4. the detector ledger includes dead time, hidden hold variance, and/or marginalized temporal gain rather than a Gaussian surrogate alone.

That is an attribution ledger. It is not a calibrated posterior probability and not a conformal interval.

### Must not be claimed

Do not call the map:

- a truth map;
- a pixelwise probability of hallucination;
- a calibrated uncertainty map unless coverage is separately proved;
- proof that a data-supported feature is correct;
- proof that a prior-supported feature is false.

A data-supported direction may still be noisy or model-mismatched. A null-space feature may be correct by luck or prior knowledge. The certificate states **source and witnessability**, not truth.

---

## 1.3 Adversarial anomaly-fidelity evaluation

### Closest prior art

- fastMRI+ provides extensive expert pathology annotations for reconstruction research, DOI [10.1038/s41597-022-01255-z](https://doi.org/10.1038/s41597-022-01255-z).
- Cheng et al. and Calivá et al. construct adversarial small abnormalities and robustly train MRI reconstruction models to reduce false negatives; journal DOI [10.59275/j.melba.2021-df47](https://doi.org/10.59275/j.melba.2021-df47).
- Bhadra et al. directly study hallucinations and prior-induced structures, DOI [10.1109/TMI.2021.3077857](https://doi.org/10.1109/TMI.2021.3077857).
- The Hallucination Index is an explicit metric for generative reconstruction errors, DOI [10.1007/978-3-031-72117-5_42](https://doi.org/10.1007/978-3-031-72117-5_42).
- The TrueCT challenge evaluates reconstruction with disease-specific detectability and lesion-volume metrics, DOI [10.1002/mp.17619](https://doi.org/10.1002/mp.17619).
- Antun et al. establish severe instabilities of learned inverse maps, DOI [10.1073/pnas.1907377117](https://doi.org/10.1073/pnas.1907377117).
- A clinical review explicitly warns that deep MRI reconstruction may obscure small lesions, DOI [10.1148/rg.220133](https://doi.org/10.1148/rg.220133).

### What survives

The distinctive evaluation is not “we test anomalies.” It is:

> anomaly insertion/deletion stratified by the anomaly’s **operator-measured versus information-null energy**, with a paired anomaly-response metric and a reconstruction-produced attribution/challenge record.

This directly tests the method’s declared mechanism: resolvable anomaly energy should be preserved by the data lock; sub-resolution content cannot be recovered honestly and must be reported as unwitnessed rather than silently erased or invented.

### Must not be claimed

Do not claim clinical lesion performance from synthetic GI anomalies. Use “localized atypical structures” or “anomaly phantoms,” not “diagnostic lesions,” unless a clinically validated cohort and reader study are added later.

---

## 1.4 Exact counting-channel ledger inside reconstruction and certification

### Closest prior art

Exact or accurate count likelihoods, Poisson inverse reconstruction, and Fisher analysis are mature. Examples include:

- Bardsley and Goldes, Poisson MAP reconstruction, DOI [10.1137/080726884](https://doi.org/10.1137/080726884).
- Grönberg, Danielsson and Sjölin, exact finite-frame count distribution and Fisher performance for a nonzero-pulse-length detector, DOI [10.1002/mp.13063](https://doi.org/10.1002/mp.13063).
- The data/prior Fisher comparison of Sampson–Melchior, arXiv [2306.13272](https://arxiv.org/abs/2306.13272).

### Prior-art verdict

I did **not** locate a published system that combines all of the following:

1. a dead-time/jitter or hidden-gain observed-information ledger for integrated photon counts;
2. a generative reconstruction structurally confined to the channel-information complement;
3. an explicit per-direction data/prior provenance output;
4. anomaly-fidelity evaluation stratified by measured/null energy; and
5. a fresh post-commitment measurement challenge.

That is the narrow unclaimed square most likely to survive. Phrase it as a specific formulation and system, not as “the first certified imaging framework.”

### Exactness fence

Use three separate labels:

- **model-exact channel identity:** the E1/O1 missing-information relation;
- **channel-exact local ledger:** Fisher/observed information at a declared scene and calibrated detector model;
- **numerically exact architectural split:** projector and prior-leak residual at machine tolerance.

Do not collapse these into “exact truth attribution.”

---

## 1.5 GT-free post-commitment audit challenge

### Closest prior art

Held-out measurement use is established:

- SSDU partitions acquired measurements into consistency and loss sets, DOI [10.1002/mrm.28378](https://doi.org/10.1002/mrm.28378).
- Multi-mask SSDU repeats this splitting, DOI [10.1002/nbm.4798](https://doi.org/10.1002/nbm.4798).
- BUQO tests support for claimed structures, DOI [10.1137/18M1173629](https://doi.org/10.1137/18M1173629).
- Posterior and conformal imaging checks are represented by DOI [10.1137/23M1546129](https://doi.org/10.1137/23M1546129).
- Measurement self-consistency has a long history; for example SPIRiT, DOI [10.1002/mrm.22428](https://doi.org/10.1002/mrm.22428).

### What survives

The operational slice is:

1. commit the reconstruction, attribution ledger, and any claimed region before the challenge;
2. generate fresh challenge patterns from a frozen hash/seed rule after commitment;
3. acquire and charge those measurements;
4. test the committed image against the fresh record;
5. return a witness/consistency result, not a truth verdict.

I did not find this exact commit-first challenge protocol in reconstruction literature. It is closer to active validation/property testing than to ordinary data consistency.

### Must not be claimed

Do not say that the challenge verifies null-space content or proves an image correct. Without fresh measurements that intersect a disputed direction, the no-free-audit obstruction remains. The challenge only tests what the added measurements make witnessable.

---

# 2. Frozen headline discipline

## 2.1 Sharpest honest one-sentence claim

Before the DEV result, use the method claim:

> **We reconstruct photon-counting images with a channel-derived information ledger that locks measurement-supported directions to the counting likelihood, confines the generative prior to the information-null complement, and reports—then optionally challenges—which structures were witnessed by data and which were supplied by the prior.**

After, and only after, a positive anomaly gate, the eye-lighting result sentence may be:

> **Every generated detail comes with a receipt: the method preserves measurement-resolvable anomalies that an unconstrained generative prior erases, while marking sub-resolution content as unwitnessed instead of silently inventing or deleting it.**

Use “channel-derived” in the title/abstract unless every deployed Fisher quantity is genuinely exact rather than plug-in or numerically approximated.

## 2.2 The overclaim most likely to kill the paper

The fatal sentence would be any variant of:

> “Our certificate determines whether every reconstructed structure is real and eliminates hallucinations without ground truth.”

That contradicts the null-space impossibility result and collides directly with existing hallucination/UQ literature. The strongest permitted semantics are:

- source attribution;
- measurement support under the declared channel;
- consistency with fresh challenge measurements;
- calibrated uncertainty only if separately established.

---

# 3. Frozen DEV gate for the running anomaly probe

## 3.1 Decision status

The current probe is DEV-only and may decide whether a confirmatory campaign is worth freezing. It cannot be used to tune the same thresholds after outcomes are inspected.

No scene, anomaly, prior checkpoint, operator, sampling mask, or seed may be replaced after this ruling. Any already exposed smoke cases remain in the analysis and receive a leave-exposed-out sensitivity only.

## 3.2 Minimum cohort and inference unit

A decisive run requires at least:

- **12 resolvable anomaly cases**, defined before reconstruction by support diameter `>=3 px` and declared measured-subspace energy fraction `eta_data >= 0.40`;
- **12 sub-resolution/unwitnessed cases**, defined before reconstruction by `eta_data <= 0.25` or null-energy fraction `>=0.65`;
- **12 anomaly-free controls**.

If the running probe has fewer independent base scenes, it is a smoke/proof-of-concept and cannot authorize a flagship campaign.

Use the base scene as the outer inference unit. Multiple anomaly shapes or contrasts on the same base scene are nested, not independent. Use five paired pattern/noise seeds per case. The primary operating point is the frozen lane-0 `rate05`, 5% operator; 10% sampling is a secondary replication and cannot rescue a 5% failure.

## 3.3 Mandatory arms

All learned arms must share the identical gain-marginalized counting data term and identical acquired record.

1. `PHYSICS_ONLY`: gain-marginalized likelihood plus the frozen nonlearned prior/regularizer.
2. `DL_UNCONSTRAINED`: the VQGAN/generative prior without structural range/null confinement.
3. `NULL_EUCLIDEAN`: the strongest working GAN_FCC range/null chassis with ordinary Euclidean/noiseless consistency and the exact dual projector.
4. `LEDGER_LOCK`: Fisher/observed-information-weighted data component, prior structurally confined to the declared information-null complement, attribution ledger, and exact dual solve.

The proposed method must beat `NULL_EUCLIDEAN`; otherwise the new counting ledger is an audit layer attached to known null-space learning, not an effective reconstruction method.

## 3.4 Paired anomaly-response metric

For each base image `x_j^0`, form an anomaly-paired image

\[
x_j^1=x_j^0+a_j,
\]

using identical measurement and noise seeds. For method `m`, define the reconstructed anomaly response

\[
\hat a_{m,j}=\hat x_{m,j}^1-\hat x_{m,j}^0.
\]

On the frozen anomaly ROI `Omega_j` plus its one-pixel dilation, define

\[
AF_{m,j}
=20\log_{10}
\frac{\|a_j\|_{2,\Omega_j}}
{\|\hat a_{m,j}-a_j\|_{2,\Omega_j}+10^{-12}}.
\]

This is the primary anomaly-fidelity score. Also report signed amplitude recovery

\[
r_{m,j}=\frac{\langle\hat a_{m,j},a_j\rangle_{\Omega_j}}
{\|a_j\|_{2,\Omega_j}^2}
\]

and anomaly localization overlap, but neither substitutes for `AF`.

For each scene define the strong comparator

\[
AF_{\mathrm{comp},j}
=
\max\{AF_{\mathrm{PHYSICS}},AF_{\mathrm{DL}},AF_{\mathrm{NULL\text{-}EUC}}\},
\]

where the maximum is taken after five-seed averaging, never per seed.

## 3.5 Gate 0 — exact implementation integrity

All must pass before image gates are read:

- relative dual/projector residual `<=1e-7`;
- prior leakage into the declared data component
  \[
  \|P_D\hat x_{\rm prior}\|_2/(\|\hat x_{\rm prior}\|_2+10^{-12})\le10^{-6};
  \]
- decomposition closure
  \[
  \|\hat x-(\hat x_{\rm data}+\hat x_{\rm prior})\|_2/\|\hat x\|_2\le10^{-7};
  \]
- identical counting likelihood, calibration, and regularization rule across arms;
- no ground truth, anomaly mask, or anomaly class enters reconstruction or attribution generation.

Failure is an implementation failure, not a negative scientific result.

## 3.6 Gate A — resolvable-anomaly improvement

On the resolvable cohort, define

\[
\Delta AF_j=AF_{\mathrm{LEDGER},j}-AF_{\mathrm{comp},j}.
\]

PASS iff all hold:

1. median `Delta AF >= 2.0 dB`;
2. at least `75%` of independent base scenes have `Delta AF > 0`;
3. the 90% scene-bootstrap lower bound on the median exceeds `+0.50 dB`.

### Required novelty subgate

Against `NULL_EUCLIDEAN` alone, require:

- median anomaly-fidelity gain `>=1.0 dB`;
- 90% bootstrap lower bound `>0`.

If the main gate passes but this subgate fails, the positive belongs to ordinary null-space confinement, not to the counting-information ledger. Do not launch the claimed flagship method in that case.

## 3.7 Gate B — background no-harm

Let `Omega_j^c` exclude the anomaly ROI and its two-pixel dilation. Against the stronger of `DL_UNCONSTRAINED` and `NULL_EUCLIDEAN`, define background radiometric-PSNR difference `Delta B_j`.

PASS iff all hold:

1. median `Delta B >= -0.50 dB`;
2. 90% scene-bootstrap lower bound `>-1.0 dB`;
3. no more than `20%` of independent scenes lose more than `1.5 dB`.

The method cannot win by preserving an anomaly while visibly degrading the rest of the image.

## 3.8 Gate C — attribution and sub-resolution flagging

The certificate must emit a frozen scalar witnessability score for every challenged direction/ROI; the score definition and threshold must be committed before outcomes.

Use true operator energy only to label DEV cases:

- `WITNESSED`: `eta_data >=0.50`;
- `UNWITNESSED`: `eta_data <=0.25`;
- the interval `(0.25,0.50)` is reported but excluded from binary classification.

PASS iff all hold:

1. AUROC for `WITNESSED` versus `UNWITNESSED` is `>=0.90`;
2. sensitivity for flagging `UNWITNESSED` is `>=90%` at the frozen threshold;
3. specificity on `WITNESSED` cases is `>=80%`;
4. the score is computed without ground truth or the injected anomaly mask.

If the implementation can calculate the score only after being handed the true anomaly direction, it is an offline analysis—not a deployable certificate—and Gate C fails.

## 3.9 Gate D — post-commitment challenge

For every audited image:

1. cryptographically commit the reconstruction, ledger, and automatically proposed challenge region(s);
2. derive the challenge seed from the commit hash;
3. acquire exactly `k=4` fresh patterns under the same detector model;
4. charge their dwell and photons separately;
5. return a witness/consistency decision without updating the displayed reconstruction.

On the sub-resolution and deliberately prior-only challenge cohort, PASS iff:

- challenge sensitivity `>=95%`;
- false-alarm rate on clean/consistent controls `<=5%`;
- no challenge is selected, repeated, or discarded after observing its measurements.

This gate validates the operational audit. It does not validate all null-space content.

## 3.10 Overall DEV verdict

```text
CERTIFIED_LEDGER_DEV_PASS
    = Gate0 AND GateA AND novelty_subgate AND GateB AND GateC AND GateD
```

All gates are non-rescuing. A 10% sampling success cannot rescue a 5% primary failure. Global PSNR, SSIM, LPIPS, challenge latency, and the previously measured `+2.2–2.5 dB` data-term gain are mandatory disclosures but not substitutes for the anomaly gate.

### Kill consequences

- **Gate A fail:** no effective reconstruction method; stop the flagship.
- **Novelty subgate fail:** known null-space method is doing the work; no counting-ledger method claim.
- **Gate B fail:** anomaly preservation is bought by unacceptable image degradation; stop.
- **Gate C fail:** attribution is not operationally meaningful; drop the certificate thesis.
- **Gate D fail:** post-commitment auditing does not work robustly; it cannot appear as a system claim.

---

# 4. One-flagship-paper scope ruling

## 4.1 Paper authorization

No manuscript work begins until `CERTIFIED_LEDGER_DEV_PASS = TRUE`. A partial pass is not enough for the proposed one-system paper.

If it passes, write **one flagship system paper**, not a merger of every ROUND63 result.

## 4.2 Load-bearing main-paper content

1. **Problem:** learned priors can erase or invent atypical structure under compressed photon-counting measurements.
2. **One theorem/object:** the channel-derived observed-information ledger and the exact architectural data/prior decomposition.
3. **Method:** gain-marginalized counting likelihood; information-weighted data lock; null-confined generative completion; attribution output.
4. **Audit:** post-commitment four-pattern challenge with its no-free-audit scope boundary.
5. **Primary evidence:** paired anomaly insertion/deletion, resolvable versus unwitnessed strata, background no-harm, and visible galleries.
6. **Limitations:** model calibration, local Fisher geometry, synthetic anomaly scope, challenge photon cost, and impossibility of truth verification in the untouched null space.

## 4.3 ROUND63 content to include only compactly

The following are load-bearing but should occupy one compact section/figure, not become a second paper inside the paper:

- E1/O1 as the mathematical source of the counting ledger;
- the measured RQL `~99%` efficiency and software-saturation table, establishing that the data-side likelihood has already spent the available record;
- the gain-marginalized temporal-prior likelihood as the data term, with its `+2.2–2.5 dB` effectiveness but no standalone novelty claim.

## 4.4 Content to cite as companion evidence

Do not reproduce full derivations of:

- the cube-root ridge and jitter-capped operating law;
- M1’s `+1.87 dB / 19.13x` verdicts;
- the RLMI, DOPS, CPL, and other negative campaigns.

Cite them as companion evidence that motivates the “ledger fully spent” premise. A one-paragraph saturation summary is enough.

## 4.5 Relationship to the operator’s other manuscripts

- **GI_a1 identifiability/audit line:** companion theory. Import only the exact no-free-audit statement and the post-commitment challenge theorem needed to make the system self-contained. Do not merge the full identifiability manuscript.
- **FOHI draft:** quarry only. The measured no-op orthogonalization must not appear as a contribution or named method.
- **GAN_FCC:** chassis and strongest baseline. Cite its prior lineage honestly; the flagship novelty must be the counting-ledger integration and the anomaly/audit result, not the existence of a null-space GAN.

---

# 5. Final no-go test

The package does **not yet** clear the operator’s bar. It has a credible architecture, a strong data term, exact numerical geometry, and a plausible narrow novelty square—but the central image-level result is still being measured.

My final decision rule is blunt:

> **Build the paper only if the Fisher/counting ledger adds at least 1 dB anomaly fidelity beyond the Euclidean null-space chassis, the complete method gains at least 2 dB over the strongest comparator on resolvable anomalies, background loss stays within the frozen no-harm band, and the GT-free flag/challenge gates pass.**

If the proposed method merely matches the Euclidean null-space network, or only supplies an attribution heat map, the package is an elegant integration of known ideas rather than the genuinely effective software method the operator requested. In that event, stop rather than rename it “certified.”