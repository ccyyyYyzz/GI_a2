# GAN_FCC completion-line audit vs. the information-guided certified-reconstruction proposal

**Date:** 2026-07-22 · **Mode:** strictly read-only on `E:\`; all writing under this directory.
**Scope:** connect SIDE 1 (operator's GAN_FCC completion line, worktree
`E:/GAN_FCC_WORK/active_code/completion_gan_round18`, branch
`codex/gan-gi-journal-poc-20260718`, commit 2697699) to SIDE 2 (HSGI
information-guided certified reconstruction, `D:/GI_another`).

---

## VERDICT (conclusion first)

**CHASSIS on the plumbing, COUSIN on the science.** The GAN line is the natural
implementation base for the *infrastructure* the certified/information-guided
direction needs — it already ships, tested and hash-audited, three of the
proposal's four load-bearing components:

1. **Hard data-consistency = subspace lock.** Exact projection onto the
   measurement-attainable set (`src/gauge_geometry.py`), not a soft penalty.
2. **Prior containment to data-silent directions.** *Every* learned correction is
   null-space-projected before use (`apply_generator`, `null_project_flat`), so
   the prior can only write where the data is silent — this is exactly "let an
   aggressive prior speak only where data is silent."
3. **Record realizability / certified provenance.** Round-61's truth-blind
   nearest box-attainable-record closure + SHA-chained receipts.
4. **A range/null "two-ledger" decomposition primitive** used throughout.

But the proposal's *distinctive scientific content is entirely absent*, and the
GAN line's geometry actively fights it:

- **(a) information-WEIGHTED consistency (the Fisher ledger).** The whole GAN
  stack is **exact, noiseless, Euclidean**. There is no noise model, no Fisher
  weighting, no covariance metric. The range/null split is Euclidean-orthogonal;
  a Fisher metric makes range and null **oblique**, so this is a genuine reframe
  of the geometry, not a parameter. It was *contemplated* (the Round-61 note's
  "tolerance tube / calibrated covariance norm") but deliberately **not built**.
- **(b) counting-channel physics.** The forward model is a deterministic
  `y = A x` with structured DCT+Hadamard+random rows and `noise_std = 0`. There
  is no Poisson/renewal/dead-time channel — the entire SIDE-2 counting-channel
  Fisher ledger has **no counterpart** here.
- **(c) attribution certificate.** The *primitives* exist (orthonormal row basis
  `Q`, singular values, an LMMSE **per-pixel posterior-std map**), but nothing is
  *shipped* as a per-direction data-vs-prior certificate.
- **(d) anomaly-fidelity evaluation.** None. The eval is mean PSNR/SSIM/LPIPS
  over 6740 natural images. It is in fact the **opposite** of anomaly fidelity: it
  rewards the LPIPS gain that a texture prior buys — the very "prior paints
  plausible detail" mode the proposal wants to catch.

So the operator's **"接不太上" intuition is right on the science** and **too
pessimistic on the plumbing**. The certified direction should *inherit the gauge
geometry + projection + provenance chassis* and *discard the FOHI/GAN
perceptual refiner* (which is the erasing prior, and whose headline mechanism is
in any case a measured no-op — see Q1/Q2). The information-weighting, the
counting physics, the attribution certificate, and the anomaly-fidelity
benchmark are all new work: the chassis saves months of already-audited
geometry/provenance code, but the *paper* is essentially new.

**Q4 makes the anomaly-fidelity claim conditionally feasible:** at 5–10% sampling
with *this* operator, a localized anomaly's *measured (range) energy fraction*
ranges from ~13% (2px point / thin stroke) to ~96% (≥3px broad bump), mean
**0.43 @ 5% / 0.60 @ 10%**. Not the pessimistic "~5–10%" — for features ≥3px the
data sees a substantial fraction, so a data-consistent reconstruction *can*
preserve them (fidelity route has legs); for sub-3px features it is
null-dominated (attribution/flagging route required). The split is a knob of the
operator, which the GAN line's spread-spectrum DCT/Hadamard rows set favorably.

---

## Q1 — MECHANISM INVENTORY (what the pipeline enforces, and where)

### Forward model / geometry
- Operator built by `gan_high_quality_gi.build_structured_operator_rows`
  (`src/dc_balanced.py`): row 0 = DC (constant), then low-frequency **DCT** rows,
  low-sequency **Hadamard** rows, then **random ±1 zero-mean** rows; each non-DC
  row zero-mean, unit-norm. Images are **64×64 grayscale STL-10** (n = 4096).
  rate05 → m = 205 (1+128+56+20, seed 772101); rate10 → m = 410
  (1+257+112+40, seed 772011). I **reproduced both operators in pure numpy and
  the rate05 SHA-256 matches the frozen receipt exactly**
  (`62305c8e…e79d6` → True), so the geometry below is the real one.
- `GaugeGeometry` (`src/gauge_geometry.py`) stores the thin factorization
  `A = U_r · diag(s) · Q` with **`Q` orthonormal row basis**. It provides the
  two-ledger primitives: `row_project_flat` (range = `Qᵀ Q`),
  `null_project_flat` (I − range), `intrinsic_record(y) = (y U_r)/s` (whitened
  record coords), `affine_project_flat` (put a vector on the fiber `Qx=z`).

### Is consistency EXACT or noise-weighted? — **EXACT / hard, Euclidean.**
The forward cache is **noiseless** (`GhostMeasurementOperator … noise_std=0.0`);
there is **no noise-weighting anywhere**. The terminal step is a hard Euclidean
projection onto the measurement-and-box-feasible set
`F_y = {x∈[0,1]ⁿ : Ax = y}`, in three interchangeable exact solvers:
- `project_box_fiber_q` — Dykstra alternating projection (used by the eval
  pipeline via `project_predictions`);
- `project_box_fiber_exact_dual` — semismooth-Newton on the low-dim dual
  (fast, ~1e-8 record error; used by the dual_v2 benchmark);
- `project_box_fiber_feasible_ray` — closed-form maximal feasible retraction.
The metric is the plain Euclidean `‖x−p‖₂` in image space; the only "whitening"
(division by `s` in `intrinsic_record`) is a coordinate convenience, **not** an
information weighting of the consistency constraint.

### Where the generator acts — **null space only.**
`FiberResidualPhaseGenerator` (`src/fiber_residual_phase_gan.py`) is an
identity-initialized local rotation of a supplied residual; it does **not** touch
the measured component. `apply_generator` calls
`geometry.null_project_flat(raw_correction)` — every correction is forced into
`null(A)` before it is added to the base and box-clipped. So the learned prior is
**structurally contained to the data-silent subspace**; the range component is
frozen by the record.

### What the discriminator sees.
`ConditionalHighPassDiscriminator` sees `concat(base, high_pass(image))` — a
patch discriminator on the **high-frequency texture** of the candidate,
conditioned on the stable GI estimate. It is a perceptual-texture critic in the
null space, not a data-consistency term. (In the frozen selected method the
adversarial weight is 0 — the discriminator was found unnecessary.)

### What FOHI added, and HOW it degenerated.
Claimed mechanism (`paper/fohi/sections/02_principle.tex`): with `c_S∈null(A)` the
VQAE "structural" coordinate and `c_G∈null(A)` the VQGAN coordinate,
`u = P_N H_{0.12,0.03} P_N (c_G − c_S)` (high-passed null-space difference), then
the **rank-one "fiber-orthogonalization"** `v = u − (⟨u,c_S⟩/‖c_S‖²) c_S` so that
`v ⟂ c_S`, and `x̂ = Proj_{F_y}(x_B + c_S + 0.5 v)`. The orthogonalization is the
paper's headline geometric novelty.

**It is a measured no-op.** Across *all* six frozen held-out cells (3 lanes × 2
rates) the diagnostics record `beta ≡ 0`, `parallel_energy_fraction ≡ 0`, and
`orthogonal_innovation_norm == innovation_norm`, and FOHI's output is
**byte-identical to the plain "fixed" high-pass arm** (`fohi_vs_fixed` = 0 in
PSNR/SSIM/LPIPS; the pairwise record certificate `fohi_minus_fixed` = 0.0,
passed). So `v = u` exactly — the rank-one removal subtracts nothing.

**Root cause (traced through the code):** the "structural" coordinate `c_S`
collapses to ≈0. In `prepare_split` the `vqae_control` arm's source direction is
built from `x_A` — the *same* reconstruction that builds `base` — so its
`direction` input is ~0, and the control adapter (which never left its
identity init: `weight≈0.10`, `rotation_rms≈8.5e-8`) emits a correction of
RMS ≈ 8.5e-8. `fiber_orthogonal_innovation` (`src/fiber_orthogonal_innovation.py`)
then hits its float32 degeneracy guard: `‖c_S‖² ≈ 3e-11` is far below
`eps₃₂·n ≈ 5e-4`, so `beta` is forced to 0. The named novelty therefore never
engages; "FOHI" reduces to a plain high-pass null-space innovation. (The claim
audit `paper/fohi/notes/final_claim_audit.md` P2-1 flags the near-zero threshold;
the measured all-zero result confirms it is not a corner case but the operating
regime.) A separate, more serious claim defect is P0-1: the frozen projection
targeted a **clipped-anchor record**, not raw `y` — which is what Round-61's
canonical closure was built to fix.

### What Round-61 "canonical attainable target" does — **confirmed: it builds
exactly-feasible measurement records.** `project_record_to_box_attainable` solves
`min_{x∈[0,1]ⁿ} ½‖Ax − y‖²` by monotone restarted FISTA (float64), returning the
box-feasible witness `x_c`, its record `y* = A x_c` (the unique Euclidean
projection of `y` onto the attainable set `C = A·[0,1]ⁿ`), and the discrepancy
`ρ = ‖y − y*‖`. This is truth-blind (uses only `A`, `y`, precision budget), and
its purpose is to remove the float32 "empty-fiber" incompatibility so the terminal
projection has a nonempty target (`theory/ROUND61_MEASUREMENT_COMPATIBILITY_NOTE.md`).
FOHI Proposition 3 (any `h` with `Ah=0` leaves `y*` unchanged) is the exact
statement that null-space innovations add no measurement information — i.e. the
formal version of "the prior only writes where the data is silent."

---

## Q2 — EMPIRICAL STATE

### Headline numbers (Round 59 held-out, 6740 STL-10 test images)
Status: **`QUALITY_EFFECT_POSITIVE_BUT_MEASUREMENT_CERTIFICATE_FAILED`, final claim
NO** (`aggregate_diagnostic/round59_diagnostic.json`). The measured effect of the
selected method vs. the matched **structural** baseline:

| rate | PSNR Δ (dB) | SSIM Δ | LPIPS Δ | crossed 95% CI (PSNR) |
|---|---|---|---|---|
| 5%  | +0.078…+0.082 (mean **+0.0805**) | +0.0045 | −0.0096 | [0.0790, 0.0820] |
| 10% | +0.206…+0.238 (mean **+0.222**)  | +0.0108 | −0.0159 | [0.219, 0.224] |

Absolute quality (Round-61, lane0): rate05 PSNR **23.30 dB** / SSIM 0.674 / LPIPS
0.285; rate10 PSNR **25.53 dB**. All 18 lane×rate×metric means favorable; six
crossed-CI and six Bonferroni directional bounds favorable.

### Solid vs. known-broken
- **Solid:** the *sign and size* of the effect (a small but consistent null-space
  refinement worth ~0.08 dB @5% / ~0.22 dB @10%); the operator/lane provenance;
  the crossed-image bootstrap (Round-59 fixed the nested→crossed bug flagged in
  the claim audit P1-4).
- **Known-broken / not-claimable:**
  - **FOHI ≡ fixed** — the orthogonalization novelty is a no-op (Q1). The entire
    "fiber-orthogonal" contribution is vacuous on the data; what remains is a
    plain high-pass null-space texture add worth ~0.08–0.22 dB.
  - **Measurement certificate fails.** The eval targets raw `y` with the Dykstra
    projector, which **never converges**: `max_relative_record_error ≈ 3.0–4.5e-4`
    at 4096 iterations, `all_converged=false`, for all arms — because the cached
    float32 `y` is not exactly attainable (empty fiber) *and* Dykstra has a slow
    feasibility tail. Positive quality but no final claim.
  - The claim-audit P0-1 (clipped-anchor vs raw-`y`), P1-1/2 (adapter coords not
    raw VQAE/VQGAN), P1-8 ("sampling ratio" is a signed-row budget, not a photon
    budget), and the all-simulation caveat remain.

### What the projection benchmark (dual_v2) established
`projection_benchmark/lane0_rate05_limit64_dual_v2.json` (64 images): the
**exact-dual** semismooth-Newton projector converges to
`max_relative_record_error ≈ 9.8e-9`, box-violation 0, complementarity ~5.5e-17
in **64 iterations, ~0.6 s/image**; the closed-form **ray** retraction is
machine-exact-feasible (~1.9e-15) in ~0.03 s but suboptimal (squared distance
72.6 vs the dual's 7.9; mean retention 0.64). **Reconciliation:** the feasible
box-fiber projection *is solved* to ~1e-8 by the exact-dual solver — the Round-59
certificate failure is a *solver-wiring* issue (the eval used Dykstra, not the
dual) *plus* the raw-`y` incompatibility that Round-61's `y*` closure removes. So
"exact hard consistency" is demonstrably achievable at ~0.6 s/image; the frozen
eval just didn't wire it into the terminal projection. (The Round-61 eval still
shows the Dykstra tail at 3.77e-4 against `z*`, i.e. the fix is designed but not
yet plumbed into the metric path.)

---

## Q3 — CONNECTION MAP (component by component)

### Already exists in the GAN line
- **Hard data consistency = subspace lock → YES.** Exact box-fiber projection
  (`project_box_fiber_exact_dual` / `_q` / `_feasible_ray`) onto `{Ax=y}∩[0,1]ⁿ`.
- **Null-space-only generation = prior containment → YES.** Every learned
  correction is `null_project`ed (`apply_generator`). The prior physically cannot
  alter the measured component.
- **Canonical targets = record realizability → YES.** Round-61
  `project_record_to_box_attainable` builds the truth-blind nearest attainable
  record `y*`, with SHA-chained receipts and a float32 precision gate.
- **Range/null two-ledger primitive → YES.** `row_project_flat` /
  `null_project_flat` and the SVD factorization are the audit substrate.

### Missing (the proposal's actual novelty)
**(a) noise/information-WEIGHTED consistency (Fisher ledger vs binary
projection).** Absent. Concretely, how hard to extend `gauge_geometry`'s geometry
to a weighted metric:
- *Easy part:* the **soft** closure already minimizes `½‖Ax−y‖²` by FISTA
  (`project_record_to_box_attainable`) — swapping to a diagonal-`W` weighted
  residual `½‖Ax−y‖²_W` (with `W` = per-bucket counting-Fisher information) is a
  few-line change to the residual/gradient. The Round-61 note itself sketches this
  as the "tolerance-tube / calibrated covariance-norm" variant.
- *Hard part:* the **exact-equality attribution machinery** is built entirely on
  an **orthonormal** `Q` and Euclidean orthogonality. Under a Fisher metric `M`,
  range and null become **`M`-orthogonal (oblique)**; `intrinsic_record`,
  `affine_project_flat`, and the exact-dual box-fiber projector all assume
  `QQᵀ=I`. Producing a *weighted range/null certificate* requires the generalized
  SVD/oblique projectors `A(AᵀWA)⁺AᵀW` etc. — a real re-derivation, moderate
  effort but not a knob. Net: **the soft-consistency change is small; the
  weighted two-ledger *audit* is a rebuild.**

**(b) counting-channel physics.** The GAN line assumes a **deterministic,
noiseless** `y = A x` (structured DCT/Hadamard/random rows). No Poisson/renewal,
no dead-time, no jitter — none of SIDE-2's channel model. To carry the Fisher
ledger you must *first* give the GAN line a stochastic counting observation model
(today it has none).

**(c) attribution / certificate outputs.** No per-direction data-vs-prior content
is shipped. The nearest existing objects: `GaugeEmpiricalAnchor` computes an LMMSE
**per-pixel posterior-std map** (`normalized_posterior_map`) — a data-uncertainty
proxy — and Probe A (SIDE 2) already does the **eigenband range/null MSE
decomposition** (`PROBE_A_REPORT.md`). The certificate is one wrapper away from
these primitives, but does not exist as an output.

**(d) anomaly-fidelity evaluation.** None, and the current objective is
adversarial to it: the method is *selected* on a joint PSNR↑/SSIM↑/**LPIPS↓**
gate, i.e. it rewards perceptually-plausible null-space texture — the exact
"prior invents plausible detail" behaviour the proposal targets as a failure. No
planted-anomaly set, no worst-case/region metric exists.

---

## Q4 — ANOMALY ENERGY SPLIT (numerical test)

Test: `anomaly_split_test.py` (+ `anomaly_split.json`). Reproduced the **real**
lane0 operators (rate05 SHA verified `62305c8e…`=match; rate10 rebuilt from
`config_rate10.yaml`), computed an orthonormal row basis `Q` by SVD
(rank 200 @5%, 393 @10%), planted 10 localized anomalies on the native 64×64
grid (Gaussian bumps σ=1–3, 2–5px squares, 5px strokes) at random positions, and
measured `range_fraction = ‖Q a‖²/‖a‖²`.

| anomaly (≈scale) | range @5% | null @5% | range @10% |
|---|---|---|---|
| point 2px / thin stroke 5px | **0.13–0.19** | 0.81–0.87 | 0.22–0.33 |
| square 3px | 0.27 | 0.73 | 0.43 |
| gauss σ1.0–1.2 | 0.32–0.42 | 0.58–0.68 | 0.55–0.68 |
| square 5px | 0.52 | 0.48 | 0.75 |
| gauss σ1.5 | 0.56 | 0.44 | 0.82 |
| gauss σ2.0 | 0.76 | 0.24 | 0.95 |
| gauss σ3.0 (≈6px FWHM) | **0.96** | 0.04 | 0.9999 |
| **mean (10)** | **0.426** | 0.574 | **0.598** |

DC contributes negligibly (0.001–0.028), so this is genuine spread-spectrum
coverage from the low-frequency DCT/Hadamard rows, **not** a mean artifact.

**Reading.** The pessimistic "~5–10% range" scenario is **false for this
operator.** Because a spatially-localized bump has broad spectral support, the
structured low-frequency rows capture 30–96% of a ≥3px anomaly at 5% sampling
(50–99% at 10%). Consequences for the proposal:
- For **resolvable anomalies (≳3px)** the data *sees* most of the anomaly →
  a data-consistent reconstruction **can** preserve it → the **fidelity route has
  legs** (the DL prior's failure to preserve it is then a certifiable error, not
  an unmeasurable one).
- For **sub-resolution anomalies (≤2px / thin strokes)** the anomaly is
  null-dominated (65–87% null) → range-locking alone cannot preserve it → these
  require the **attribution/flagging route** (declare "prior-determined, not
  data-supported here"), not a fidelity guarantee.
- The split is an **operator-design knob** (spread-spectrum ↑ range fraction).
  The GAN line's DCT/Hadamard-heavy design is favorable; a random-speckle GI
  operator would push more anomaly energy into null — a lever the proposal should
  own explicitly.

---

## Q5 — VERDICT (three-way, with the minimal delta)

**(i) CHASSIS — YES for the infrastructure layer.** Reuse, largely as-is:
`src/gauge_geometry.py` (range/null projectors, exact box-fiber projectors,
LMMSE anchor + posterior-std), the Round-61 canonical-record closure, and the
SHA-receipt/truth-blind provenance discipline. These already implement
subspace-lock + null-space prior-containment + certified realizable record —
three of the proposal's four pillars — and are audited to ~1e-8 record error at
~0.6 s/image. This saves the certified direction its entire geometry+provenance
build.

**Minimal delta to reach the information-guided certified method:**
1. **Add a counting observation model** to the GAN forward path (Poisson/renewal
   with dead-time/jitter), porting SIDE-2's channel; today `noise_std=0`. *[new;
   moderate]*
2. **Weight the consistency** — change the FISTA closure residual to `‖Ax−y‖²_W`
   with `W` = per-bucket counting-Fisher weights (the Round-61 "tolerance-tube").
   *[small for the solver]*
3. **Weighted two-ledger audit** — generalize `intrinsic_record` /
   `affine_project_flat` / exact-dual projector from orthonormal-Euclidean to the
   `M = AᵀWA` oblique metric (GSVD / oblique projectors). *[real re-derivation;
   moderate]*
4. **Ship a per-direction attribution certificate** — wrap the (weighted) SVD
   directions + LMMSE posterior-std into a data-vs-prior content ledger per
   output. *[small–moderate; primitives exist, incl. Probe-A eigenbands]*
5. **Anomaly-fidelity benchmark** — planted-anomaly test set + worst-case/region
   fidelity metric; `anomaly_split_test.py` is a seed. *[new; small]*
6. **Drop** the FOHI orthogonalization (degenerate no-op) and **repurpose, not
   reuse,** the perceptual GAN refiner: as the "aggressive prior where data is
   silent" it must be gated by the attribution certificate, **not** by an LPIPS
   texture gate (the current selection objective is the failure mode the proposal
   indicts).

**(ii) COUSIN — YES on the headline science.** The operator's "接不太上" is correct
here: noiseless-exact-Euclidean consistency and Fisher-weighted consistency are
different animals; the counting-channel physics lives only on SIDE 2 (the GAN
line has no forward-noise model at all); and the GAN refiner's whole purpose —
perceptually pleasing null-space texture selected on LPIPS — is the *opposite* of
anomaly preservation. Merging the *method claims* (not the plumbing) would be
forced.

**(iii) INDEPENDENT — the FOHI/GAN *result* is a dead end for this direction, but
two things are worth importing either way:** (a) the **exact box-fiber projection
+ certified-record provenance** machinery (chassis, above); (b) the empirical
**range/null decomposition tooling** and the Q4 finding that operator design sets
the anomaly range-fraction — which tells the proposal *which* anomalies are a
fidelity target vs a flagging target, and gives it a design lever.

**Bottom line:** inherit the chassis, discard the refiner and its degenerate
novelty, and build the four missing scientific pieces (counting model,
Fisher-weighted consistency, attribution certificate, anomaly-fidelity eval).
The infrastructure reuse is substantial and real; the *paper* is new.

---

### Evidence index (all read-only)
- Geometry / projectors: `…/completion_gan_round18/src/gauge_geometry.py`
- Generator / discriminator: `…/src/fiber_residual_phase_gan.py`,
  `train_fiber_residual_phase_gan.py`, `train_vqae_centered_residual_adapter.py`
- FOHI orthogonalization + degeneracy guard: `…/src/fiber_orthogonal_innovation.py`;
  eval `…/diagnose_fiber_orthogonal_highpass_innovation.py`
- Operator builder: `…/gan_high_quality_gi.py` (`build_structured_operator_rows`),
  `…/src/dc_balanced.py`; configs `…/artifacts/gan_gi_journal_round59_recovery/lane0/…/config_rate05.yaml`,`…/config_rate10.yaml`
- Round-61 theory: `…/theory/ROUND61_MEASUREMENT_COMPATIBILITY_NOTE.md`
- Empirical: `…/experiments/gan_gi_journal_round59/aggregate_diagnostic/ROUND59_DIAGNOSTIC.md` + `round59_diagnostic.json`;
  `…/round61/lanes/lane0/rate05/fohi/summary.json`;
  `…/round61/canonical_targets/lane0_rate05/CANONICAL_ATTAINABLE_TARGET_RECEIPT.json`;
  `…/round61/projection_benchmark/lane0_rate05_limit64_dual_v2.json`
- Claim audit trail: `…/paper/fohi/notes/final_claim_audit.md` (P0-1, P2-1),
  `project_truth.md`, `decision_log.md`
- SIDE 2: `D:/GI_another/docs/SOFTWARE_SATURATION_VERDICT.md`,
  `docs/THEORY_MAINLINE.md`, `results/round63_next/EFFICIENCY_PROBE/PROBE_A_REPORT.md`
- This audit's test: `anomaly_split_test.py`, `anomaly_split.json` (this dir)
