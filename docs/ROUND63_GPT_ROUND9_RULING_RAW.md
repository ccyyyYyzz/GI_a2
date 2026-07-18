# R9 — Final ruling

| Question | Frozen ruling |
|---|---|
| Q1 verdict presentation | Keep “When high flux helps…,” but replace **phase diagram** with **operating map**. Primary remains formally negative. Contrast ordering is post hoc descriptive evidence, not a validated boundary. |
| Q2 patch-sparse patterns | The \(k^{-1/2}\) argument is correct only for uniformly scattered supports. Patch/multiscale patterns are a strong follow-up, but discussion-only here. |
| Q3 adaptivity | The stated necessity conjecture is false or vacuous without a restricted scene class. A class-conditional adaptive advantage is a follow-up theorem, not part of this paper. |
| Q4 interferometric cancellation | Physically sound only as **pre-detection optical displacement/nulling**, not electronic balanced subtraction. Strong PolyU collaboration anchor; one forward-looking paragraph here. |
| Q5 scope | S3 core and S4 exact-vs-RQL are load-bearing. NATURAL-24 is secondary external validity. S2B/S2C are supplement-if-complete, otherwise droppable. No Study 3 experiments in this manuscript. |

---

# Q1 — Verdict presentation calibration

## 1. Frozen title

Use:

> **When high flux helps single-pixel imaging: a contrast–dead-time operating map**

Do not retain **phase diagram**. The paper has compelling regime structure, but it did not preregister or estimate a universal critical \(C_u\), and the Study-2 primary endpoint is formally negative because the breadth criterion failed. “Operating map” preserves the strong conditional claim without implying a thermodynamic-style phase boundary established by a confirmatory threshold.

The opening “When high flux helps…” remains justified: Study 1 gives a dense-pattern null regime, while Study 2 gives large acquisition-speed effects in 16/24 scenes and a fixed-dwell secondary that passes all three criteria.

## 2. Frozen hierarchy

The first sentence of the Study-2 verdict must be:

> **The preregistered Study-2 primary endpoint was formally negative because its breadth criterion failed: 16 of 24 scenes accelerated, below the required 18 of 24.**

Only after that sentence may the paper report:

- median \(S_{\rm gate}=6.796\);
- bootstrap lower bound \(2.835\);
- fixed-dwell gain \(+1.16\) dB;
- 19/24 fixed-dwell-positive scenes;
- fixed-dwell bootstrap lower bound \(+0.72\) dB.

Do not use:

- “near miss”;
- “essentially passed”;
- “positive except for breadth”;
- “where measurable” without defining the denominator;
- any subgroup gate based on \(C_u\);
- any logical combination in which the secondary rescues the primary.

## 3. Contrast-ordering claim

The family ordering is scientifically valuable but was computed after unblinding and contains only six family-level points. It therefore supports a mechanism interpretation, not a confirmatory boundary estimate.

Replace “the phase boundary passes through the confirmatory cohort” with the following frozen wording:

> **After unblinding, a descriptive family-level analysis showed that the six family success fractions were monotonically ordered by their corresponding mean measured bucket contrasts. This six-point ordering was not preregistered, did not enter any decision rule, and does not estimate a universal critical contrast; it is consistent with the confirmatory cohort spanning a contrast-controlled transition.**

The phrase “predicted boundary” should not be used. The theory predicted a count-information ridge and a qualitative contrast limitation, not a numerical scene-side \(C_u\) cutoff.

## 4. Frozen abstract replacement

Replace the current Study-2 result passage with:

> In Study 2, the preregistered acquisition-speed endpoint was formally negative because its breadth criterion failed (16/24 scenes with \(S_{\mathrm{gate}}>1\), versus the required 18/24), although both effect-size criteria passed (median \(6.80\times\); family-stratified bootstrap lower bound \(2.84\times\)). The preregistered fixed-dwell secondary passed (median \(+1.16\) dB at the terminal dwell, 19/24 scenes positive, lower bound \(+0.72\) dB). A post hoc descriptive family-level analysis found that the six success fractions were monotonically ordered by measured bucket contrast, consistent with—but not establishing—a contrast-controlled transition. Across preregistered no-gate occupancy controls, the fixed-dwell gain was nonmonotone and largest among the tested rungs at \(3.1\%\) occupancy, exposing a contrast–conditioning trade-off.

This reports the result forcefully without converting a failed primary into a positive claim.

## 5. Promote \(k=32\) to the abstract?

**Yes, as a no-gate descriptive result.**

Use:

> **Among the tested no-gate occupancy controls, the largest median fixed-dwell gain occurred at \(k=32\) (\(3.1\%\) occupancy), rather than at the sparsest pattern, revealing competing benefits of bucket contrast and sensing-operator conditioning.**

Do not call \(k=32\):

- “the optimum”;
- “the optimal sparsity”;
- “the sweet spot” without “among the tested rungs”;
- a confirmatory result.

The result is useful precisely because it prevents an oversimplified “more sparse is always better” story.

## 6. Mandatory manuscript corrections

There are four hard textual errors in the current draft:

1. **Delete “grew monotonically along the contrast ladder.”**  
   The ladder is explicitly nonmonotone: \(k=32\) exceeds \(k=16\) and raster under the frozen summary convention. The current abstract contradicts the results section.

2. Replace every occurrence of:
   - “fixed-photon-budget gain”;
   - “fixed incident photon budget” when comparing \(\rho=0.6\) with \(\rho=0.05\);

   with:

   > **fixed-dwell gain**  
   > or  
   > **fixed optical-integration-time gain**

   At equal \(\nu\), the high-\(\rho\) arm uses greater incident flux. Equal incident dose was enforced **across occupancy \(k\) at a fixed \((\rho,\nu)\)**, not between safe and high-flux operating points.

3. The current results text says the failures include “half of the contour family.” The frozen outcome is contour \(1/4\) successful, hence **three of four contour instances failed**.

4. Resolve the numerical-source mismatch before submission. The R9 briefing lists ladder medians \(+4.17\) and \(+3.64\) dB, while the cab9197 manuscript and figure convention use \(+4.13\) and \(+3.56\) dB. Use the frozen analyzer’s source-of-record convention everywhere, including abstract, captions and supplement.

---

# Q2 — Scale-matched patch-sparse illumination

## 1. Is the \(k^{-1/2}\) analysis correct?

Yes, under the specific scattered-support model.

Let \(x\ge0\), \(\sum_jx_j=1\), and let a pattern select exactly \(k\) pixels uniformly without replacement, with the Study-2 equal-load scaling

\[
a_j=\frac nk\,\mathbf 1_{\{j\in S\}}.
\]

Then \(E[u]=1\), where \(u=a^\top x\), and

\[
\boxed{
C_u^2
=
\operatorname{Var}(u)
=
\frac{n-k}{k(n-1)}
\left(n\|x\|_2^2-1\right).
}
\]

Defining the pixel-level coefficient of variation

\[
C_{\rm pix}^2=n\|x\|_2^2-1,
\]

for \(k\ll n\),

\[
\boxed{
C_u\approx \frac{C_{\rm pix}}{\sqrt{k}}.
}
\]

Thus the scattered support averages approximately \(k\) independent spatial samples.

This does **not** hold universally for a contiguous patch. For a stationary scene with covariance \(R_x(\Delta)\), a patch average has variance proportional to

\[
\frac{1}{k^2}
\sum_{p,q\in {\cal P}}R_x(p-q),
\]

or spectrally,

\[
\int S_x(\omega)|H_{\cal P}(\omega)|^2\,d\omega.
\]

Positive within-patch correlation reduces the effective number of independent samples, increasing bucket contrast relative to scattered support. But a patch larger than the relevant feature or correlation scale can still average away fine texture. The correct proposal is therefore a **multiscale support bank**, not one supposedly optimal patch size.

## 2. Prior-art risk

Local, foveated, wavelet-tree, spatially variant-resolution and multiscale Hadamard illumination already exist in SPI/ghost imaging. Their stated purposes include adaptive resolution, ROI allocation, reduced sampling and improved reconstruction efficiency.

The defensible future-work distinction is narrower:

> local or multiscale supports optimized for **bucket-contrast preservation under equal detector load, equal dose and dead-time constraints**.

Do not claim invention of patch, block, wavelet or multiscale single-pixel patterns.

## 3. Current-paper scope

Discussion paragraph only. No new simulation, figure or theorem should enter this manuscript.

Frozen wording:

> For uniformly scattered fixed-\(k\) supports, \(C_u^2=\frac{n-k}{k(n-1)}(n\|x\|_2^2-1)\), giving the observed \(k^{-1/2}\) spatial-averaging law. A local support instead averages only the effective number of independent correlation cells within the patch and may therefore preserve modulation from scale-matched structure. Local, foveated and multiscale single-pixel patterns have previously been developed for resolution and measurement efficiency; optimizing their spatial scale for equal-load dead-time contrast is a distinct future design problem.

---

# Q3 — Is adaptivity provably necessary?

## 1. The conjecture as written is rejected

The statement

> every fixed ensemble has a structured nonflat scene with \(\Gamma<1\), while an adaptive oracle achieves \(\Gamma>1\) for every structured scene

is either false or vacuous.

First, without a lower bound on structural amplitude, take

\[
x_\epsilon=x_{\rm flat}+\epsilon v.
\]

For every fixed or adaptive finite-dose strategy,

\[
C_u,\Gamma\rightarrow0
\quad\text{as}\quad
\epsilon\rightarrow0.
\]

No adaptive method can uniformly guarantee \(\Gamma>1\) for every merely nonflat scene.

Second, a sufficiently complete fixed ensemble—raster being the simplest example—need not possess a nonflat spatial null direction. Whether it has poor statistical power depends on dose, minimum contrast and loss, not fixedness alone.

Third, the failure of a linear-Gaussian no-adaptation lemma’s assumptions does not imply that adaptivity strictly helps. Adaptive sensing theory contains both strong-gain examples for structured sparse classes and lower bounds showing little or no improvement in other regimes.

## 2. Defensible theorem target

A useful follow-up must freeze:

- a restricted scene class;
- a minimum structural contrast;
- nonnegative and equal-dose pattern constraints;
- a peak-irradiance constraint;
- an exact renewal observation channel;
- a precise loss or testing problem.

The most promising class is a union of unknown spatial scales:

\[
{\cal X}=\bigcup_{\ell=1}^{L}{\cal X}_{\ell},
\]

where every \(x\in{\cal X}_{\ell}\) has a specified minimum contrast at correlation scale \(\ell\).

A plausible theorem program is:

1. any fixed design distributing measurements over all \(L\) scales has worst-case structural Chernoff/Fisher information of order \(M/L\);
2. a two-stage design uses \(m_0=O(\log L)\) coarse measurements to identify the active scale;
3. the remaining \(M-m_0\) measurements use matched local patterns and attain a constant fraction of the scale-oracle information;
4. under appropriate separation conditions, the adaptive-to-nonadaptive ratio scales as \(L\).

That would establish a class-conditional adaptive advantage. It would not establish universal necessity.

## 3. Paper placement

This is the anchor of a follow-up, not a theorem for the present paper.

Frozen discussion sentence:

> The present results do not establish that adaptive illumination is necessary. Although dead-time nonlinearity makes the statistically efficient pattern scale scene-dependent, a strict adaptive advantage requires a specified scene class, resource constraint and loss. A natural follow-up is to seek a minimax or Chernoff-information separation between fixed multiscale illumination and a two-stage scale-estimation-plus-matched-pattern policy.

---

# Q4 — Interferometric common-mode cancellation

## 1. Physical ruling

The central idea is sound only if cancellation occurs **optically before photon counting**.

Let

\[
E_s=E_0+\delta E
\]

be the bucket field and \(E_r\) a shaped reference. A beamsplitter produces output rates proportional to

\[
\lambda_\pm\propto
\frac12|E_s\pm E_r|^2.
\]

With \(E_r=E_0\), the dark port has

\[
\lambda_-\propto\frac12|\delta E|^2,
\]

so common-mode photons are diverted before reaching the SPAD and do not consume its dead-time budget.

However, an exact null has a major drawback: sensitivity to \(\delta E\) is quadratic and loses its sign/quadrature. A deliberately biased displacement

\[
E_r=E_0-b
\]

gives

\[
\lambda_-\propto
|b+\delta E|^2
=
|b|^2+2\operatorname{Re}(b^*\delta E)+|\delta E|^2.
\]

The bias restores a linear structural term but consumes part of the dead-time budget. The key design variable is therefore **optimal displacement bias**, not maximum cancellation.

Electronic subtraction of two already-recorded bright SPAD outputs does not solve the problem: each detector has already incurred dead-time loss before subtraction.

## 2. What breaks first?

In order:

1. **Complex-field matching.**  
   The reference must match the common-mode field’s amplitude, phase, polarization and spatial mode, not merely its mean intensity.

2. **Scattering-mode mismatch.**  
   Through a diffuse or multimode scattering medium, a single reference mode cannot generally null the collected field. Wavefront shaping or single-mode collection is required.

3. **Phase drift.**  
   Dynamic scattering changes the complex common mode. The reference must be updated faster than its coherence time.

4. **Reference error.**  
   A small phase or amplitude error leaves residual DC that can again dominate the dark port.

5. **Bright-port saturation.**  
   The sum port should probably use an analog or high-dynamic-range monitor rather than another SPAD.

6. **Vanishing absolute information.**  
   \(C_u\) can become arbitrarily large as the dark-port mean approaches zero while the absolute structural photon rate and Fisher information approach zero. Therefore “\(C_u\to\infty\) for every structured scene” is not a performance guarantee.

## 3. Photon statistics

For a coherent input field and deterministic linear optics, each output port is conditionally a coherent field, so its ideal arrival process remains Poisson and a dead-time renewal model is appropriate.

Under time-varying speckle or reference drift within an exposure, the rate itself is random. The marginal process is then a Cox or mixed-Poisson process rather than a homogeneous Poisson process; the present exact renewal likelihood applies only conditionally or when exposures are short compared with the field’s coherence time.

## 4. Prior art

Interferometric, phase-shifting, homodyne and complex-field single-pixel imaging are established. Prior systems include compressive holography with a Mach–Zehnder interferometer, full-wavefront SPI by homodyne detection, spatial-coherence interferometric SPI, multicore-fiber interferometric SPI, and recent common-path interferometric single-pixel photothermal imaging.

Hao and Chen’s 2025 work retrieves complex fields through scattering media with single-pixel measurements, which provides a particularly relevant route to estimating the reference field, but it is not a dead-time-optimized optical nulling receiver.

The defensible Study-3 contribution is:

> **dead-time-aware optical-displacement single-pixel photon counting, with null depth optimized for structural Fisher information under SPAD saturation.**

Do not pitch it as the first interferometric SPI system or the first dark-port single-pixel detector.

## 5. Collaboration ruling

**Approved as the strongest PolyU collaboration pitch.** It directly joins:

- Wen Chen’s complex-field/scattering expertise;
- your dead-time operating-point theory;
- a hardware question that cannot be settled by simulation alone.

Recommended bench sequence:

1. static coherent amplitude/phase target without scattering;
2. static diffuser with a shaped reference;
3. optimize displacement bias versus dead time;
4. dynamic scattering with reference tracking.

Current-paper inclusion: one forward-looking discussion paragraph only.

Frozen wording:

> A coherent extension would use a coarse complex-field estimate to synthesize an optical displacement that suppresses the predicted common mode before photon counting. Unlike electronic balanced subtraction, a dark-port receiver removes photons before dead-time loss. Exact nulling, however, gives quadratic residual sensitivity and demands spatial-mode matching and phase stability; the relevant operating variable is therefore an optimized displacement bias rather than maximal null depth. Testing this receiver first under static and then dynamic scattering is a natural optical-bench extension.

---

# Q5 — Remaining scope triage

## Frozen priority order

\[
\boxed{
\text{S4 exact reference}
>
\text{S3 core mismatch}
>
\text{NATURAL-24 summary}
>
\text{stop}.
}
\]

## 1. S4 exact-vs-RQL: load-bearing

**Run it. Do not drop it.**

The exact Fisher ridge validates a scalar information calculation; it does not validate that the RQL image reconstruction tracks exact-likelihood reconstruction under TV regularization. Because the paper’s central practical method is RQL, the \(8^2/16^2\) EXACT+TV comparison closes an obvious reviewer question.

Placement:

- main text: one compact table or one sentence reporting the exact–RQL gap in the deployment zone and its growth in the short-window/extreme-load corner;
- supplement: full grid, optimization diagnostics and multistart checks.

No new gate is needed.

## 2. S3 mismatch: load-bearing core

A simulation-only detector paper requires a credible mismatch section.

Main-text anchors should include:

- \(\hat\tau\) error \(\pm10\%\);
- calibrated background/dark rate at the 10% level;
- afterpulsing at 2–5%;
- paralyzable data reconstructed with the nonparalyzable model.

Report radiometric degradation and flux bias. Put the complete OAT/interactions/start-mode/jitter grids in the supplement.

S3 should remain descriptive. Do not introduce a late robustness pass gate.

## 3. NATURAL-24: finish, but secondary

Finish the analysis because it is already near completion and provides external-validity honesty.

Main text needs only:

- the proportion of natural scenes that are safe-range uninformative;
- median fixed-dwell gain;
- one representative image pair;
- a sentence stating that smooth/prior-limited scenes occupy the low-benefit regime.

Put all curves and per-image tables in the supplement.

## 4. S2B and S2C: not load-bearing

- If already substantially completed, place them in the supplement as geometry/ensemble sensitivity.
- If not completed, they may be dropped without weakening the central argument.
- Do not delay submission for 128², broad pattern-family sweeps or additional sampling-rate grids.
- Do not promote an S2B/S2C result into the abstract.

The occupancy ladder already demonstrates the key contrast–conditioning trade-off.

## 5. Manuscript scope lock

### Main paper

1. finite-window exact information ridge;
2. RQL reconstruction;
3. Study 1 dense-pattern negative;
4. Study 2 formally negative primary plus passed fixed-dwell secondary and mixed scene response;
5. occupancy ladder and contrast–conditioning interpretation;
6. compact S3 mismatch summary;
7. compact exact-vs-RQL validation;
8. limitations and collaboration outlook.

### Supplement

- full NATURAL-24;
- complete S3;
- complete S4;
- no-gate ladder tables;
- dynamic-scattering detail;
- S2B/S2C only if already available;
- descriptive residual audit.

### Explicitly excluded

- patch-sparse execution;
- adaptive illumination theorem;
- interferometric Study 3;
- new image classes;
- new endpoints or gates;
- another attempt to obtain a formally positive primary result.

The paper’s final scientific position is not “Study 2 passed.” It is:

> **Dense illumination is multiplex-limited; sparse equal-load illumination produces large and statistically supported benefits for many, but not sufficiently many, preregistered scenes to pass the breadth gate. The scene response is nevertheless structured by bucket contrast, while the tested occupancy ladder exposes an independent conditioning penalty.**

That mixed verdict, combined with the exact ridge and fixed-dwell secondary, is a defensible and more interesting operating-map paper than a post hoc all-positive story.