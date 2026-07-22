# R24 ruling (GitHub issue #16, raw)

Title: R24 ruling: deployable residual-gated illumination method
Posted: 2026-07-22T00:47:41Z

---

# R24 — Method ruling: from residual-gated OED to a deployable illumination router

**Reference commit:** `1db245b`  
**Source question:** `docs/ROUND63_GPT_ROUND24_QUESTION.md`

**Scope:** pre-implementation method architecture for the next paper-2 line. M1 and both existing papers remain untouched. No campaign is authorized by this ruling; the DEV bridge in §4 must pass first.

> **R24 VERDICT: do not build RGAI as written as the lead method.**
>
> The residual gate is valuable, but the proposed online Frank–Wolfe branch is the wrong deployment architecture. It is slow, it certifies a local Fisher surrogate rather than image quality, and it amplifies pre-scan/model error exactly on the scenes for which it is invoked. The R18 support-expanding primal is strong evidence that geometry headroom exists; it is not evidence that online D-optimal geometry will improve RQL images.
>
> **Conditional GO:** re-architect the method as a **two-shot, estimator-aware, finite-library router with abstention**, called below **Residual-Gated Library Illumination (RGLI)**. Use online FW only as an offline library generator and DEV oracle. The global power knob is library element 0 and is selected whenever a robust rank-one sufficiency test says adaptation is unnecessary. Geometry adaptation is selected from a small preloaded safety-compliant bank using a reconstruction-aware local risk, not D-optimality alone.
>
> **Hard stop:** if the true-scene OED oracle fails the image-level bridge test in §4, do not launch M2 and do not lead paper 2 with an adaptive method. That failure would mean the Fisher-to-image gap is structural in the present architecture.

## Executive ruling

| question | ruling |
|---|---|
| **Q1 — RGAI** | **Reject as the deployed lead; retain as diagnostic/oracle.** Adaptive OED, uncertainty-driven acquisition, amortized policies, and FW/equivalence certificates are established components. The narrow potentially new element is using a dead-time-aware rank-one residual as a route/abstention gate under physical illumination constraints. “Self-certifying imaging” is not a defensible phrase: the bound certifies the local information surrogate, not PSNR or RQL risk. |
| **Q2 — better method** | **Adopt RGLI:** 52-pattern pre-scan → local RQL-TV covariance → robust residual gate → choose one preloaded pattern bank and one power schedule → acquire once → RQL. No online full-dictionary solve. Online FW is an oracle and library-construction tool only. |
| **Q3 — DEV go/no-go** | Run the two-stage bridge in §4. Require a true-scene OED image gain first, then a deployable RGLI gain and no-harm control. Any hard gate failure kills M2. |
| **Q4 — M2** | Only after the bridge passes. Use a fresh 24-scene stress/general cohort, a strong composite fixed baseline, one rescue primary and one no-harm requirement. M1 remains the frozen rank-one special-case evidence. |
| **Q5 — hardware** | The plugin claim is conditional. DMD banks must be preloaded; switch/blanking time must enter wall-clock accounting; programmable hold-off is not universal and does not imply `c=0`; every pre-scan and discriminator photon is charged. |

---

# 1. Q1 — Attack on RGAI

## 1.1 The make-or-break defect is worse than R1 states

R23 gives

\[
 f(\Phi^*)-f(\Phi_k)\le \frac{\varepsilon_k^2}{2m}
\]

for a **local information objective** under a declared model, scene linearization, control set, and monotone-concave kernel branch. It does not bound

\[
\operatorname{PSNR}(\hat x),\quad
\operatorname{SSIM}(\hat x),\quad
\text{TV bias},\quad
\text{or image-domain risk}.
\]

Therefore the proposed phrase “the method knows when it is not needed” is currently too strong. What it knows is narrower:

> under the plug-in local information model, the rank-one controller is within a computable information gap of the full declared control class.

That can be a routing feature. It is not an image certificate.

The conservative modulus `m` is also likely to be smallest exactly on weakly conditioned scenes, so the bound may become largest or vacuous where the gate is supposed to help. Report both `epsilon_k` and the actual bound; never hide a tiny `m` behind a normalized residual.

## 1.2 Unlisted failure modes

### F1 — plug-in certification is not robust certification

`epsilon_1`, the leverage profile, and the atom ordering are evaluated at `xhat`. A 52-pattern pre-scan can alter the active face, change the normal-cone projection, and reorder the best atoms. A certificate at `xhat` need not cover the true scene. The route needs an uncertainty envelope or an abstention state.

### F2 — route discontinuity

A hard threshold on `epsilon_1` can flip the entire pattern bank under an arbitrarily small pre-scan perturbation. This creates unstable decisions even if the final reconstruction is stable. Use a margin/hysteresis rule and abstain to the robust fixed bank when the decision is not separated.

### F3 — local OED can overfocus on the wrong inverse-problem geometry

D-optimality rewards volume in the declared task subspace. RQL-TV image quality may instead be dominated by a few edge modes, radiometric bias, or a regularization-null direction. The contour failure is consistent with this problem. The design objective must see the reconstruction curvature, not only the detector Fisher matrix.

### F4 — online FW is least reliable when invoked

The large-residual scenes are the scenes on which the local linearization is most suspect. Running more FW steps does not repair model error; it can confidently optimize the wrong local objective.

### F5 — finite-row realization and calibration

A continuous design or a sparse FW measure must be rounded into physical rows. Rounding changes dose balance, peak power, spectral conditioning, and the residual. The deployed certificate must be recomputed after materialization and after measured pattern-intensity calibration.

### F6 — policy staleness

For dynamic scattering or scene motion, the pre-scan and main acquisition do not describe the same object. An adaptive pattern bank can be worse than a scene-independent bank because it is stale. The method needs a change detector or an explicit static-scene assumption.

### F7 — detector-model drift

`tau`, effective hold jitter, PDE, background, and afterpulsing depend on temperature and operating point. A single scalar `c` calibrated once is not enough if the detector is driven far from the calibration regime.

### F8 — selection-policy likelihood

A deterministic two-stage policy is statistically legitimate when reconstruction conditions on the acquired patterns and pre-scan. A sequential policy that repeatedly adapts to noisy counts is more delicate operationally and analytically; the pattern-selection history must be retained as part of the data record. This is another reason to prefer one post-scan batch over per-pattern adaptation.

### F9 — source and modulator constraints

Arbitrary per-pattern gains require a fast, calibrated AOM/EOM or equivalent source control. Many SPI benches can change binary patterns rapidly but cannot settle an independently chosen optical power on every row with the same speed and accuracy.

### F10 — out-of-dictionary scenes

A residual can indicate that rank one is inadequate without identifying a useful geometry inside the available bank. A method needs an `OUT_OF_LIBRARY` state rather than always returning an adaptive pattern set.

## 1.3 Prior-art ruling

The broad loop “pre-scan, estimate uncertainty/information, choose the next measurements, reconstruct” is crowded:

- Info-Greedy Sensing chooses sequential compressed measurements by conditional information gain ([Braun, Pokutta & Xie](https://arxiv.org/abs/1407.0731)).
- Active MRI selects measurements during inference to reduce reconstruction uncertainty ([Zhang et al., CVPR 2019](https://openaccess.thecvf.com/content_CVPR_2019/html/Zhang_Reducing_Uncertainty_in_Undersampled_MRI_Reconstruction_With_Active_Acquisition_CVPR_2019_paper.html)).
- Edge-promoting adaptive Bayesian design chooses new tomographic views from a TV-type posterior approximation ([Helin, Hyvönen & Puska](https://epubs.siam.org/doi/10.1137/21M1409330)).
- AdaSense uses posterior samples from a diffusion model to choose future measurements in CS, MRI, and CT ([Elata, Michaeli & Elad, ECCV 2024](https://www.ecva.net/papers/eccv_2024/papers_ECCV/html/10059_ECCV_2024_paper.php)).
- DAD amortizes sequential Bayesian experimental design into a one-pass policy ([Foster et al., ICML 2021](https://proceedings.mlr.press/v139/foster21a.html)); Step-DAD adds test-time refinement ([Hedman et al., ICML 2025](https://proceedings.mlr.press/v267/hedman25a.html)); PASOA jointly updates the posterior and the sequential design ([Iollo et al., ICML 2024](https://proceedings.mlr.press/v235/iollo24a.html)).
- SPI-specific precedents include adaptive illumination from image dictionaries ([Feng et al.](https://arxiv.org/abs/1806.01340)), adaptive pattern ordering with an image dictionary, learning-based adaptive Fourier under-sampling ([Huang et al., Optics Letters 2023](https://pubmed.ncbi.nlm.nih.gov/37262260/)), adaptive regional SPI ([Jiang et al., Optics Express 2017](https://pubmed.ncbi.nlm.nih.gov/28788943/)), and adaptive real-time SPI ([Zhu et al., Optics Letters 2024](https://doi.org/10.1364/OL.514934)).

The certificate component is also classical at the abstract level. The general equivalence theorem gives optimality certificates for design measures ([Pukelsheim, Ch. 7](https://epubs.siam.org/doi/10.1137/1.9780898719109.ch7)); cutting-plane equivalence methods are longstanding ([Gribik & Kortanek](https://epubs.siam.org/doi/10.1137/0132021)); Frank–Wolfe gaps are standard primal-dual stopping certificates, including stochastic estimators of the gap ([Négiar et al., ICML 2020](https://research.google/pubs/stochastic-frank-wolfe-for-constrained-finite-sum-minimization/)).

### Narrow novelty that may survive

I did not find an existing method that combines all of the following:

1. an exact/renewal dead-time information kernel;
2. a rank-one global-flux policy as the default controller;
3. a per-scene normal-cone residual used to **route or abstain** between rank-one and higher-rank illumination;
4. a safety-compliant nonnegative SPI pattern bank;
5. a post-route information-gap disclosure.

That is a plausible narrow method claim. Do not call the generic FW/equivalence bound novel, and do not call the output “self-certifying image quality.” Use:

> **residual-gated, dead-time-aware illumination with a local information-sufficiency certificate.**

---

# 2. Q2 — Better architecture: Residual-Gated Library Illumination (RGLI)

## 2.1 Why library routing beats online FW here

Online FW is attractive mathematically and poor operationally:

- current latency is 60–240 s;
- it depends most strongly on the noisy local scene estimate;
- it optimizes the surrogate most aggressively;
- every result then needs rounding and guard repair.

A finite bank can be built once using the full solver, rounded once, calibrated once, and preloaded on the controller. Runtime becomes one reconstruction plus a small number of matrix evaluations. This is much more likely to survive R2, R3, and R5.

Library routing is not novel by itself; the novelty must be the **dead-time residual gate plus robust estimator-aware routing and abstention**.

## 2.2 Offline bank

Construct a public detector-normalized library of `K=8` exact pattern banks for the first prototype:

1. `L0`: balanced SCAT32 with one global ridge-aware power multiplier;
2. `L1`: scattered `k=16`;
3. `L2`: compact/clustered `k=16`;
4. `L3`: compact `k=32`;
5. `L4`: horizontal/vertical bar mixture;
6. `L5`: radial/ring/contour bank;
7. `L6`: multiscale block bank;
8. `L7`: a geometry bank produced offline by the full safety-constrained FW solver on a mixed prototype set.

The exact contents may change during DEV, but the deployment principle should not: a small, finite, physically materialized bank, not a live full-dictionary solve.

Every bank must already satisfy:

- nonnegative/binary implementability;
- exact physical row count;
- incident budget and peak constraints;
- per-pixel dose guard;
- declared spectral/A-risk floors where used;
- detector-load admission and ceiling guards.

The source scaling is recomputed from calibrated `(tau,c,nu)`; the spatial rows do not require per-lab retraining.

## 2.3 Runtime algorithm

### Step 0 — detector characterization, not a mandatory clamp

Input calibrated `tau`, effective hold-law `c`, PDE/background, afterpulse setting, and source/modulator limits. If the module exposes a programmable hold-off, it may be used as a detector-conditioning option, but the effective hold law must be remeasured afterward. Do **not** assume `c -> 0`.

### Step 1 — charged balanced pre-scan

Acquire the 52 balanced rows. They count toward optical time, incident photons, detected counts, and reconstruction. Produce `xhat` with the same RQL pipeline used at the end.

### Step 2 — local reconstruction-aware state

Form the smoothed local RQL-TV curvature

\[
\mathcal H_0
=F_{\rm pre}(\hat x)+\lambda_{\rm TV}R_\delta(\hat x)+\epsilon I,
\]

where `F_pre` uses the dead-time kernel and `R_delta` is the frozen differentiable TV/IRLS curvature used only for local design scoring.

For each candidate bank `L_k`, compute its predicted future information `F_k` and the local task-weighted A-risk

\[
\mathcal R_k(\hat x)
=\operatorname{tr}\!\left[W_{\rm task}
(\mathcal H_0+F_k)^{-1}\right].
\]

For the first implementation use `W_task=I` on the declared reconstruction subspace. This is still a local surrogate, but it is closer to image MSE than log-det volume and directly incorporates the reconstruction curvature. It is the same estimator-aware idea used in adaptive Bayesian inverse-problem design; the dead-time kernel and route architecture are the specialized contribution.

### Step 3 — robust gate

Generate a small deterministic uncertainty set around the pre-scan estimate—for example, `S=16` Laplace/parametric-bootstrap draws using frozen seeds. For every draw compute:

- the rank-one residual and its R23 information-gap bound;
- `R_k` for all library banks.

Define

\[
U_1=\max_s\frac{\varepsilon_1(x^{(s)})^2}{2m(x^{(s)})},
\qquad
\bar R_k=\max_s\mathcal R_k(x^{(s)}).
\]

Route as follows:

1. **KNOB:** choose `L0` if `U_1/r <= 0.01` and
   \[
   \bar R_0\le1.02\min_k\bar R_k.
   \]
2. **LIBRARY:** otherwise choose the bank with smallest `bar R_k` only if its robust advantage over the second-best bank is at least `2%`.
3. **ABSTAIN:** if the uncertainty set changes the winner, the margin is below `2%`, or no bank satisfies all guards, use the robust fixed `L0` bank and flag `OUT_OF_LIBRARY_OR_UNCERTAIN`.

These constants are prototype values, not campaign thresholds. Their role is to make the architecture unambiguous before implementation.

### Step 4 — one main batch

Project the selected preloaded bank once. No per-pattern sequential adaptation. Reconstruct using all pre-scan and main measurements with RQL.

### Step 5 — output disclosures

Ship:

- selected route and bank;
- calibrated `(tau,c,nu)` and achieved loads;
- robust library-risk margin;
- `epsilon_1`, `m`, and the R23 information-gap bound;
- guard status;
- `OUT_OF_LIBRARY_OR_UNCERTAIN` if applicable;
- optical time, wall time, incident photons, and detected counts.

Call the bound a **local information-sufficiency certificate**. It is not a PSNR certificate.

## 2.4 Role of online FW and learned amortization

- **Online FW:** DEV oracle and offline bank generator only. It measures available surrogate headroom and library coverage.
- **Tiny network:** optional runtime ablation after the deterministic router works. It may predict the bank index, but the analytic scorer must verify or override it. DAD and Step-DAD already establish the amortization idea; a network is plumbing, not the method claim.
- **Estimator co-design:** included through the local RQL-TV A-risk. Do not jointly retrain a reconstruction network; that would destroy the no-retraining/plugin claim and invite domain-shift failure.
- **Sequential acquisition:** reject as the core. Two-shot acquisition—balanced pre-scan then one selected bank—is the deployable architecture.

## 2.5 Comparative ruling

| criterion | online-FW RGAI | RGLI |
|---|---|---|
| three failure modes | theoretically broad, practically fragile | bank-limited but explicit abstention |
| R1 survival | low–medium; pure D-opt | medium–high; estimator-aware scoring + direct DEV test |
| runtime | unacceptable without amortization | seconds or subsecond with preloaded banks/low-rank updates |
| pre-scan robustness | hard route on one plug-in scene | robust envelope + abstention |
| guarantee | local D-gap only | local risk ranking within library + D-gap disclosure |
| hardware | online pattern generation/loading | one preloaded bank switch |
| novelty | components largely known | narrow dead-time residual-gated routing claim is plausible |

**Recommendation:** build RGLI, not production RGAI. Keep RGAI as the oracle that can kill or validate the library concept.

---

# 3. Q3 — Decisive DEV image-level bridge

This is a method-kill experiment, not a campaign.

## 3.1 Cohort

Generate **16 fresh DEV-only `32 x 32` scenes**, disjoint from every confirmatory set:

- 4 contour/outline scenes with weak bucket contrast;
- 4 analytic two-population witnesses: a bright broad region that contributes little task leverage plus a dim fine structure carrying the task;
- 4 low-amplitude microtexture/mixed-scale scenes;
- 4 aligned controls from maze/glyph families.

No rejection or replacement after reconstruction.

## 3.2 Physics and resources

- `M_total=1024=52+972`;
- `tau=50 ns`;
- `nu in {200,2000}`;
- `c in {0,0.05}`;
- five paired pattern/noise seeds;
- identical elapsed dwell and pattern count;
- each adaptive arm has incident budget no larger than the corresponding RIDGE-SCAT32 arm and obeys the same peak/dose guards;
- every pre-scan photon is charged.

The hard bridge decision is made at `(c=0.05, nu=2000)`. Other cells diagnose generality.

## 3.3 Arms

1. `SCAT32-060`;
2. `RIDGE-SCAT32`;
3. **TRUE-X FW oracle:** full safety-constrained geometry optimization using the true DEV scene; this is nondeployable and tests Fisher-to-image validity;
4. **XHAT FW:** original RGAI geometry using the 52-pattern estimate;
5. **RGLI:** the proposed library router;
6. **ORACLE-LIB:** best library bank selected by true radiometric PSNR; this tests bank coverage separately from routing.

For each scene define the strong fixed comparator

\[
Q_{\rm base,j}
=\max\{Q_{\rm SCAT32-060,j},Q_{\rm RIDGE-SCAT32,j}\}.
\]

All gains below are paired radiometric-PSNR gains over this comparator.

## 3.4 Hard gates

### Gate A — does Fisher-guided geometry improve images at all?

On the 12 misalignment-stress scenes, TRUE-X FW must satisfy all three:

1. median gain `>= 1.0 dB`;
2. at least `9/12` scenes positive;
3. 90% scene-bootstrap lower bound on the median `> 0.50 dB`.

**If Gate A fails, kill RGAI, RGLI, and M2.** The local Fisher geometry is not an effective image method in this setting.

### Gate B — can the deployable router retain the gain?

RGLI must satisfy all five:

1. median gain on the 12 stress scenes `>= 1.0 dB`;
2. at least `9/12` positive;
3. 90% scene-bootstrap lower bound `> 0`;
4. median gain at least `60%` of the TRUE-X FW median gain;
5. among the four aligned controls, median loss `>= -0.25 dB` and no more than one scene below `-1.0 dB`.

**If Gate B fails, do not launch M2.** A library or router redesign would be a new method line, not a minor implementation fix.

### Gate C — is routing informative?

Using the ORACLE-LIB image gain as a DEV-only label:

- among scenes where ORACLE-LIB beats the knob by `>1 dB`, RGLI must choose a non-knob bank in at least `80%`;
- among aligned controls where ORACLE-LIB advantage is `<0.25 dB`, RGLI must choose KNOB or ABSTAIN in at least `75%`.

This gate tests the route, not just the bank.

### Gate D — latency

On the declared reference CPU, excluding final reconstruction but including pre-scan feature extraction and all bank scoring:

- median route latency `<1.0 s`;
- 95th percentile `<2.0 s`.

If this fails, the plugin claim fails even if PSNR improves.

## 3.5 Required decomposition

Report the four gaps separately:

1. TRUE-X FW vs fixed baseline — surrogate-to-image bridge;
2. XHAT FW vs TRUE-X FW — pre-scan/localization loss;
3. ORACLE-LIB vs TRUE-X FW — finite-bank approximation loss;
4. RGLI vs ORACLE-LIB — routing loss.

This decomposition identifies the failure layer. Do not average them into one “adaptive gain.”

---

# 4. Q4 — Minimal M2 if and only if DEV passes

## 4.1 Cohort

Exactly 24 fresh `32 x 32` scenes:

- **12 stress:** four contour, four two-population, four low-amplitude microtexture/mixed-scale;
- **12 general:** two fresh instances from each of the six established structural families.

Five paired seeds. No acceptance filter.

## 4.2 Conditions

Primary detector condition:

- calibrated `c=0.05`;
- `nu=2000`;
- `M_total=1024`, including pre-scan;
- identical elapsed dwell;
- RGLI incident budget capped by the RIDGE-SCAT32 budget;
- actual incident/detected photons and peak dose disclosed.

Secondary condition: `nu=200`. A `c=0` arm is a no-gate detector ablation.

## 4.3 Confirmatory arms

1. RGLI;
2. RIDGE-SCAT32;
3. SCAT32-060;
4. one preregistered robust fixed-library bank.

Online FW and true-scene oracles are DEV/supplement diagnostics only and never confirmatory arms.

The primary comparator is the preregistered paired composite

\[
Q_{\rm fixed,j}=\max\{Q_{\rm RIDGE,j},Q_{0.60,j}\}.
\]

This is deliberately hard: the method must beat the better fixed operating policy on each scene.

## 4.4 Primary rescue endpoint

On the 12 stress scenes at `nu=2000`, define

\[
\Delta Q_j=Q_{\rm RGLI,j}-Q_{\rm fixed,j}.
\]

PASS iff all hold:

- median `Delta Q >= 1.0 dB`;
- family/stratum-preserving bootstrap lower bound `>0`;
- at least `9/12` scenes positive.

## 4.5 Mandatory no-harm endpoint

On the 12 general scenes:

- median `Delta Q >= -0.25 dB`;
- bootstrap lower bound `>-0.50 dB`;
- no more than `2/12` scenes lose more than `1 dB`.

The method claim requires both rescue and no-harm. The no-harm endpoint cannot be replaced by an average over all 24 scenes.

## 4.6 Secondary disclosures

- route fractions: KNOB / LIBRARY / ABSTAIN;
- route latency and pattern-bank switching latency;
- robust risk margin and R23 information-gap bound;
- out-of-library rate;
- incident/detected photon ratios;
- `nu=200` replication;
- post-hoc correlation between predicted local risk gain and actual PSNR gain, explicitly not a gate.

## 4.7 Paper-2 architecture if M2 passes

1. **Problem:** the global ridge knob is powerful but scene-dependent.
2. **Law:** jitter-aware optimal load and the rank-one residual theorem.
3. **Method:** RGLI; the knob is its rank-one/default bank.
4. **DEV bridge:** true-scene oracle proves the Fisher-to-image bridge; decomposition shows where deployability is lost or retained.
5. **M2:** rescue plus no-harm verdicts.
6. **M1:** compact frozen evidence for the rank-one special case, not the new paper’s lead campaign.
7. **Deployment:** power, switching, calibration, abstention, and limitations.

If M2 fails, paper 2 must remain the jitter/global-control paper. Do not rescue the method with a new cohort, learned router, or altered image endpoint.

---

# 5. Q5 — Deployment reality check

## 5.1 DMD switching is not free

TI specifies maximum 1-bit rates of `32,552 Hz` for DLP7000 and `23,148 Hz` for DLP9500 ([DLP7000](https://www.ti.com/product/DLP7000), [DLP9500](https://www.ti.com/product/DLP9500)). Ideal pattern times are therefore about `30.7 us` and `43.2 us`.

Consequences:

- 52 pre-scan patterns require at least about `1.6–2.2 ms` of DMD display time;
- 1024 patterns require at least about `31–44 ms`, before detector dwell, blanking, source settling, USB/PCIe transfer, or reconstruction;
- at `tau=50 ns, nu=200`, optical dwell is only `10 us/pattern`, so DMD switching dominates wall time;
- at `nu=2000`, dwell is `100 us/pattern`, so switching remains a substantial overhead.

The method must report

\[
T_{\rm wall}=\sum_i(T_i+t_{\rm switch,i}+t_{\rm blank,i}+t_{\rm source,i})
+t_{\rm route},
\]

not optical dwell alone. Library banks should be preloaded. Loading newly generated FW patterns from the host during an acquisition can erase the claimed speedup.

TI support also documents pattern-transition transients in some DLP7000 configurations; optical blanking/guard intervals must be characterized rather than assumed negligible.

## 5.2 SLM class matters

The proposed scheme is plausible for a high-speed binary DMD. It is much less attractive for liquid-crystal SLMs with frame rates orders of magnitude lower. A two-shot bank switch survives; per-pattern sequential adaptation does not.

## 5.3 Hold-off control is not universal

Commercial detector behavior varies:

- ID Quantique’s ID230 exposes adjustable dead time from `2 us` to `100 us` ([official product page](https://www.idquantique.com/quantum-detection-systems/products/id230/)).
- Excelitas SPCM-AQRH modules are sold with model-dependent fixed dead times around `24–42 ns`, rather than a universal runtime hold-off register ([official product page](https://prod.excelitas.com/product/spcm-aqrh)).

A programmable dead time is therefore not a plugin assumption. Moreover, imposing a deterministic minimum hold does not prove the effective hold distribution has `c=0`; avalanche recovery, afterpulsing, electronics, and gating can leave state dependence. Longer hold-off also changes mean dead time, maximum count rate, and afterpulse probability. Recalibrate the effective `H`, then recompute the ridge.

Remove “clamp so `c -> 0` by construction” from the core method. Replace it with:

> if the detector exposes programmable hold-off, choose and characterize a detector setting jointly with the illumination policy.

## 5.4 Pre-scan accounting

The 52-pattern pre-scan is about `5.1%` of a 1024-pattern acquisition, but its photon fraction is not necessarily 5.1% if its load differs from the main bank. Charge:

- every pre-scan exposure;
- every tie-break/discriminator exposure if later added;
- all source power used during blanking or settling where physically relevant;
- all pre-scan counts in the final reconstruction.

Do not quote a method speedup that excludes pre-scan optical time or decision latency.

## 5.5 Additional bench requirements

- calibrate actual row-wise optical power, not just binary masks;
- retain the selected pattern history and measured source setting in the data file;
- verify DMD switching transients and blanking;
- enforce peak irradiance and sample photodamage limits;
- monitor scene change between pre-scan and main batch;
- expose an ABSTAIN/fixed-bank fallback in the user interface;
- distinguish detector timing jitter from hold-time variability—these are not the same `c`.

---

# 6. Final blunt assessment

**RGAI is worth keeping as a scientific oracle, not yet worth deploying.** Its attractive sentence—“compute the residual, then run a few FW steps if needed”—hides the central unresolved issue: the residual and the FW objective live in Fisher space, while the claimed product lives in reconstructed-image space.

**RGLI is worth building conditionally** because it makes three disciplined concessions:

1. it routes with an estimator-aware local risk;
2. it moves expensive optimization offline into a preloaded finite bank;
3. it abstains when the pre-scan cannot support a stable decision.

The self-certifying angle is only defensible in this narrower form:

> the method reports when the rank-one policy is locally information-sufficient, when a banked geometry is predicted to reduce local reconstruction risk, and when neither conclusion is stable enough to trust.

The DEV oracle gate is non-negotiable. If true-scene OED cannot buy at least a robust `1 dB` image gain on the misaligned scenes, stop. No amount of routing, amortization, or neural plumbing can rescue an objective that does not cross into image quality.