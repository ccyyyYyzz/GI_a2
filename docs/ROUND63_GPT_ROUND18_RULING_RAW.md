# R18 ruling (GitHub issue #10, raw)

Title: R18 ruling: full-stack certificate localization and final refreeze
Posted: 2026-07-19T12:00:02Z

---

# R18 — Final ruling: full-stack certificate localization and final refreeze

**Audit target:** commit `c67ba07`, especially:

- `docs/ROUND63_GPT_ROUND18_QUESTION.md`;
- `results/round63_m1/R18_GAP_PROBE_LOG.txt`;
- `code/round63/dev_gap_probe.py`;
- the R17 architecture frozen in issue #9 and `docs/ROUND63_METHOD_SPEC_M1_R17_AMENDMENT.md`.

**Scope:** this ruling amends only the still-unfrozen M1 follow-up campaign. No confirmatory scene, pre-scan count, or endpoint has been opened. It does not reopen paper 1 or alter the R14–R16 detector theory.

> **CURRENT FREEZE VERDICT: HOLD, with one final conditional refreeze path.**
>
> Adopt option **(iii)**, but rename and sharpen it. The dose-only certificate is retired as a categorical gate because a concrete feasible design already refutes near-optimality of SCAT32 in that class. The confirmatory certificate is redefined over the **full deployed conditioning stack**—incident budget, ±5% dose, A-risk, and spectral floor—while the dose-only class is retained as a no-gate development analysis that localizes where the design headroom lives.
>
> The empirical ridge primary and Q90 speed secondary remain unchanged. One coherent refreeze is still mandatory.

## Executive ruling

| Question | R18 ruling |
|---|---|
| **Q1 — certified class** | **Adopt option (iii), amended.** Replace `DOSE_SAFE_CERT_PASS` with `FULL_STACK_CERT_PASS`, targeting the complete deployed stack. Keep a no-gate `DOSE_ONLY_HEADROOM_DEV` analysis. The current constructive design proves that dose uniformity alone does not collapse geometry selection. |
| **Q2 — dual construction** | Use a **global conic/multiplier certificate**, not atom-level admission. A-risk is represented by a Schur-complement LMI and the spectral guard by an LMI. A certified cutting-plane outer approximation is permitted. The hard bar remains `G_stack/r <= 1e-2`. |
| **Q3 — paper wording** | The R17 uniform-dose-collapse sentences do **not** survive verbatim. The corrected architecture is: dose uniformity controls exposure but leaves geometry headroom; the conditioning anchors are the operative restriction; global ridge tracking remains a separate power-for-time control. |
| **Q4 — compute/failures** | Every certificate cell is classified as `CERTIFIED`, `COUNTEREXAMPLE`, or `NUMERICAL_UNRESOLVED`. All 480 must be `CERTIFIED` for `FULL_STACK_CERT_PASS`. Freeze a 420 s per-cell maximum with one deterministic rescaled retry. No imputation or cell deletion. |
| **Q5 — DEV finding** | **Allowed in the paper as explicitly DEV-only descriptive evidence and future-work motivation.** It adds no arm or endpoint. However, the current probe does not establish that the optimum uses zero gain atoms or intentionally leaves half the budget unused; its multiplicative update cannot activate initially zero columns. |

---

# 1. Preliminary audit of the new evidence

## 1.1 The near-optimality refutation is decisive

For the two audited DEV cells, the probe produced concrete probability designs over the frozen `D_cert` satisfying the incident budget and ±5% dose band, with log-determinant improvements over deployed SCAT32 of

\[
\frac{\Delta\log\det V}{r}=0.7758
\quad\text{and}\quad
0.9026.
\]

Therefore the corresponding D-efficiency ratios are at least

\[
\exp(0.7758)\approx2.17,
\qquad
\exp(0.9026)\approx2.47.
\]

A single feasible construction is sufficient to refute the proposed `G_full/r <= 10^{-2}` near-optimality claim over `C_dose`. The loose or failed dual calculations are immaterial to that logical conclusion.

Accordingly, `DOSE_SAFE_CERT_PASS` is structurally incompatible with the audited dose-only target class and must not be frozen as an expected-failure gate under its old name.

## 1.2 Mandatory caveat on the reported composition

The current `primal_probe` initializes positive mass **only** on the lowest `D_load` column and sets every other load/gain column exactly to zero. Its multiplicative update

```python
xi = xi * sqrt(score / mean_score)
```

is support-preserving: a zero column can never become positive. Thus the following are valid:

- a geometry-only, lowest-load feasible construction exists;
- that construction already gives the reported lower-bound headroom;
- gain coupling is not necessary to refute the dose-only certificate.

The following are **not** established by the current probe:

- that the dose-only optimum assigns zero mass to `D_gain`;
- that the optimum prefers the lowest palette load;
- that the optimum intentionally spends only half the photon budget.

Freeze the paper wording as **constructive existence**, not optimal-composition inference. Any later composition claim requires a support-expanding primal method—positive initialization on all admitted columns or Frank–Wolfe atom injection—with the algorithm frozen before inspecting its result.

---

# 2. Q1 — Certified-class ruling

## 2.1 Retire the misleading verdict name

Delete:

```text
DOSE_SAFE_CERT_PASS
```

Replace it with:

```text
FULL_STACK_CERT_PASS
```

The new name is mandatory because A-risk and spectral conditioning—not dose alone—are part of the certified target class.

## 2.2 Full deployed-stack class

For each local pre-scan estimate `xhat`, dwell `nu`, and photon-budget corner `b`, retain the frozen expanded dictionary

\[
\mathcal D_{\rm cert}=\mathcal D_{\rm load}\cup\mathcal D_{\rm gain}.
\]

Let

\[
V(\xi)=V_{\rm pre}+M_{\rm main}\sum_{a\in\mathcal D_{\rm cert}}
\xi_a H_a+\epsilon_0 I
\]

on the already frozen `r=200` task subspace. Let `V_fix` be the load-matched deployed SCAT32 information matrix for the same pre-scan, dwell, budget, subspace, and numerical ridge.

Freeze

\[
\boxed{
\mathcal C_{\rm stack}=\left\{\xi:\begin{array}{l}
\xi\ge0,\quad\mathbf1^T\xi=1,\\
c^T\xi\le b,\\
U^+\xi\le0,\quad U^-\xi\le0,\\
\operatorname{tr}\{V(\xi)^{-1}\}
\le1.05\operatorname{tr}(V_{\rm fix}^{-1}),\\
V(\xi)\succeq0.5V_{\rm fix}
\end{array}\right\}.
}
\]

The pointwise peak, ceiling, support, weight, and RQL-trust restrictions remain atom-admission guards. The A-risk and spectral conditions are **global design constraints** and must not be approximated by dropping individual atoms.

This class is comparator- and subspace-relative by construction. The manuscript must state that the certificate preserves the predeclared conditioning of SCAT32 on the directions SCAT32 resolves; it is not a global optimum over all images, subspaces, risk tradeoffs, or physical patterns.

## 2.3 Confirmatory verdict

Run the certificate at the same 480 frozen cells:

- 24 confirmatory scenes;
- five pre-scan/noise seeds;
- `nu in {200,2000}`;
- `b in {0.05,0.60}`.

Freeze:

```text
FULL_STACK_CERT_PASS
```

PASS iff every cell is:

1. deployed-design feasible for budget, dose, A-risk, and spectral constraints;
2. numerically resolved as `CERTIFIED`;
3. certified at
   \[
   G_{\rm stack}/r\le10^{-2}.
   \]

This remains a confirmatory structural secondary and cannot rescue either empirical ridge verdict.

## 2.4 Dose-only companion

Retain the dose-only class only as a no-gate development analysis with source-of-record fields:

```text
DOSE_ONLY_PRIMAL_GAP_PER_R
DOSE_ONLY_D_EFFICIENCY_LOWER
DOSE_ONLY_DUAL_UPPER_IF_AVAILABLE
DOSE_ONLY_SOLVER_STATUS
```

The current values are lower bounds on available headroom, not estimates of the global dose-only optimum.

No confirmatory `DOSE_ONLY_*_PASS` verdict is permitted.

---

# 3. Q2 — Full-stack certificate construction

## 3.1 Normative support-function certificate

Let `V_d` denote deployed SCAT32 information in one cell and define the local atom sensitivity

\[
d_a=M_{\rm main}\operatorname{tr}(V_d^{-1}H_a).
\]

Let

\[
\bar d_d=\operatorname{tr}\left[
V_d^{-1}(V_d-V_{\rm pre}-\epsilon_0 I)
\right].
\]

By concavity of `log det`, for every feasible design `xi`,

\[
\log\det V(\xi)-\log\det V_d
\le d^T\xi-\bar d_d.
\]

Define

\[
\boxed{
G_{\rm stack}
=\sup_{\xi\in\mathcal C_{\rm stack}}
\{d^T\xi-\bar d_d\}.
}
\]

Then

\[
\sup_{\xi\in\mathcal C_{\rm stack}}
\left[\log\det V(\xi)-\log\det V_d\right]
\le G_{\rm stack},
\]

and

\[
\operatorname{Eff}_D(V_d;\mathcal C_{\rm stack})
\ge\exp(-G_{\rm stack}/r).
\]

The existing `10^{-2}` bar therefore retains its exact interpretation:

\[
G_{\rm stack}/r\le0.01
\quad\Longrightarrow\quad
\operatorname{Eff}_D\ge e^{-0.01}=0.99005.
\]

## 3.2 Exact conic representation

The spectral constraint is the LMI

\[
V(\xi)-0.5V_{\rm fix}\succeq0.
\]

The A-risk constraint is represented exactly by introducing a symmetric matrix `W`:

\[
\begin{bmatrix}
V(\xi)&I\\
I&W
\end{bmatrix}\succeq0,
\qquad
\operatorname{tr}W
\le1.05\operatorname{tr}(V_{\rm fix}^{-1}).
\]

The Schur complement makes this equivalent to the declared A-risk cap.

The resulting support-function problem is a conic linear program. Its dual contains:

- one nonnegative budget multiplier;
- upper/lower pixel-dose multipliers;
- a PSD matrix price for the spectral LMI;
- a PSD block-matrix price and trace multiplier for the A-risk LMI;
- a full-dictionary maximum reduced sensitivity.

This multiplier/conic formulation is approved. **Restricted atom admission is rejected** for A-risk or spectral enforcement because neither condition is an atomwise property.

## 3.3 Approved scalable implementation

A direct SDP is the mathematical reference, but a certified cutting-plane outer approximation is permitted:

1. **Spectral cuts.** For every discovered violating generalized eigenvector `z`, add
   \[
   z^T\{V(\xi)-0.5V_{\rm fix}\}z\ge0.
   \]
2. **A-risk cuts.** Since
   \[
   h(V)=\operatorname{tr}(V^{-1})
   \]
   is convex with gradient `-V^{-2}`, at a trial matrix `V_k` add the necessary tangent inequality
   \[
   h(V_k)-\operatorname{tr}\{V_k^{-2}(V-V_k)\}
   \le1.05h(V_{\rm fix}).
   \]
3. Solve the resulting linear master with budget/dose prices and column generation.
4. Scan the entire frozen dictionary for reduced-sensitivity violations.
5. Add the most violated atom/cut and repeat.

Finite cut sets define an **outer relaxation** of `C_stack`; therefore a small master upper bound remains conservative. Do not call a large bound a counterexample unless a feasible primal design is constructed.

## 3.4 Mandatory primal discriminator

Before spending the full dual budget, run a frozen support-expanding primal search over `C_stack`.

- If it finds a feasible design with actual
  \[
  [\log\det V(\xi)-\log\det V_d]/r>10^{-2},
  \]
  classify the cell as `COUNTEREXAMPLE` and the certificate gate is refuted for that cell.
- A primal lower bound never certifies near-optimality; it only refutes it.

## 3.5 Toy suite

Before refreeze, require four independent finite-dictionary toys:

1. active budget + active dose, with A-risk/spectral inactive;
2. active A-risk only;
3. active spectral floor only;
4. budget, dose, A-risk, and spectral all active simultaneously.

For each toy:

- solve the primal independently by exhaustive gridding or an independent conic solver;
- require primal–dual objective agreement `<=1e-8`;
- require minimum PSD residual `>=-1e-9`;
- require normalized complementarity residual `<=1e-8`;
- require the full finite-dictionary reduced-score scan residual `<=1e-8`.

The all-active toy must contain at least one design that is dose-feasible but rejected by A-risk or spectral conditioning, matching the real architectural distinction.

---

# 4. Q3 — Frozen manuscript correction

## 4.1 Status of the R17 permitted sentences

The three R17 §6.2 sentences do not all survive verbatim.

### Sentence 1

Old:

> Within the frozen nonnegative illumination dictionary, local pre-scan linearization, incident budget, and ±5% per-pixel dose constraint, the balanced SCAT32 design was tested against a full constrained dual upper bound.

**Amend.** Replace “incident budget and ±5% dose constraint” with the complete conditioning stack.

### Sentence 2

Old:

> Under the declared dose-safety constraints, per-pattern load adaptation provided no certified material D-information advantage over the balanced design...

**Does not survive as written.** Dose-only geometry adaptation has a constructive advantage. The sentence may be retained only in the narrower, contingent full-stack form below.

### Sentence 3

Old:

> The useful surviving control was the global operating point rather than spatial redistribution of incident dose.

**Amend.** Global operating-point control remains the deployed positive arm, but geometry selection is now known to retain dose-feasible headroom outside the conditioning stack. Do not present the two as an exhaustive either/or.

## 4.2 Corrected Act III beats

Freeze the new four-beat structure:

1. **Information rewards scene-adapted geometry.** Dose-relaxed OED and the dose-only constructive probe expose substantial local D-information headroom.
2. **Uniform dose is not enough to eliminate that headroom.** A balanced cumulative exposure can coexist with scene-adapted measurement directions.
3. **Conditioning safeguards define the deployable class.** The A-risk and spectral anchors test whether the geometry headroom survives reconstruction-trust requirements.
4. **Global operating-point control is orthogonal and survives.** Ridge tracking changes one global source multiplier and is evaluated by the empirical power-for-time endpoints.

## 4.3 Frozen manuscript wording

Use the following, contingent on the corresponding results:

> A development-only constructive probe showed that the ±5% per-pixel dose constraint alone did not eliminate scene-adapted geometry selection: a feasible low-load geometry design achieved a substantially larger local determinant than SCAT32 in the audited cells. This was a lower-bound construction, not an estimate of the dose-only optimum.

> The confirmatory certificate therefore targeted the full deployed stack: the frozen dictionary and local pre-scan linearization, incident budget, ±5% per-pixel dose band, A-risk cap, and spectral floor. Global ridge tracking was evaluated separately as the time-limited, power-available control.

If `FULL_STACK_CERT_PASS` passes:

> Within the frozen full-stack class at the preregistered anchor cells, balanced SCAT32 was certified within 1% local D-efficiency of the best admissible design. The useful deployed positive control was the global detector operating point; geometry optimization under relaxed conditioning remains a separate opportunity.

If it fails, the last paragraph is forbidden. Report the distribution of certified, counterexample, and unresolved cells without a categorical collapse claim.

## 4.4 Corrected final headline

Replace the R17 final sentence

> Strict uniform-dose safety collapses the tested adaptive design space...

with:

> **Uniform dose controls cumulative exposure but does not eliminate scene-adapted geometry selection. The operative restriction is the full deployed conditioning stack: balanced SCAT32 is tested for certified near-optimality within that frozen class, while time-limited image gains are tested separately by global ridge tracking.**

After a successful certificate, “is tested for” may be replaced by “is certified within 1% local D-efficiency.” No stronger substitution is permitted.

## 4.5 Figure and handbook edits

- Act III panel (a): retain the path-specific guard-collapse DEV plot.
- Add or replace one DEV panel with the **constructive dose-only primal lower bound**, explicitly marked development-only.
- The confirmatory certificate panel must be labeled `FULL_STACK_CERT_PASS`, not `DOSE_SAFE_CERT_PASS`.
- The handbook must say “full deployed conditioning class,” not “expanded dose-safe class.”
- Delete the figure-caption and discussion claims that uniform dose itself collapses adaptive concentration.

---

# 5. Q4 — Numerical and compute rules

## 5.1 Per-cell status taxonomy

Every one of the 480 cells receives exactly one terminal status:

```text
CERTIFIED
COUNTEREXAMPLE
NUMERICAL_UNRESOLVED
```

Definitions:

- `CERTIFIED`: deployed design and certificate class are primal-feasible; all residual checks pass; `G_stack/r <= 1e-2`.
- `COUNTEREXAMPLE`: a concrete full-stack-feasible design has actual logdet improvement per dimension `>1e-2`.
- `NUMERICAL_UNRESOLVED`: neither conclusion was obtained within the frozen numerical protocol.

`FULL_STACK_CERT_PASS` requires 480/480 `CERTIFIED`. Either of the other statuses makes the categorical gate fail. `NUMERICAL_UNRESOLVED` means “not certified,” not “scientifically refuted.”

All cells remain in the source-of-record table; no retry based on image quality and no imputation are allowed.

## 5.2 Frozen compute budget

Per certificate cell, after the pre-scan estimate and fixed context are built:

1. support-expanding primal counterexample screen: **60 s**;
2. primary full-stack dual/cutting-plane solve: **180 s**;
3. one deterministic rescaled retry, only if unresolved: **180 s**.

Freeze

```text
CERT_CELL_WALL_CAP = 420 seconds
```

The retry must use the same class and threshold, with only frozen numerical equilibration:

- Cholesky whitening/scaling by the load-matched SCAT32 matrix;
- objective scores scaled by `r`;
- budget coefficients scaled by `b`;
- dose coefficients scaled by deployed mean dose.

No third solver attempt and no threshold change are permitted.

## 5.3 Solver-failure disclosure

Every cell records:

```text
solver_status_primary
solver_status_retry
wall_seconds
coefficient_dynamic_range
n_dictionary_scans
n_arisk_cuts
n_spectral_cuts
full_dictionary_scan_residual
primal_gap_lower_per_r
dual_gap_upper_per_r
min_generalized_eigenvalue
arisk_ratio
MU_CAP_ACTIVE
```

A HiGHS/SDP unknown status, time limit, nonfinite objective, or coefficient range above `1e10` after frozen scaling leads to the deterministic retry. If the retry fails, the cell is `NUMERICAL_UNRESOLVED`.

If a finite dual-multiplier cap is active, certification is accepted only after a cap ×100 stability rerun changes the bound by at most `1e-4*r` and all residuals still pass; otherwise the cell is unresolved.

## 5.4 Prefreeze feasibility gate

The current `9/9 REFREEZE_READY` ledger is superseded by R18.

Before tagging, run the full-stack protocol on the complete six-image DEV family set, seed `0`, at all four anchor combinations

\[
\nu\in\{200,2000\},
\qquad b\in\{0.05,0.60\},
\]

for **24 DEV certificate cells**.

The categorical certificate remains in the campaign only if:

- 24/24 are `CERTIFIED`;
- no full-stack primal counterexample exceeds `1e-2`;
- median certificate wall time is `<=120 s`;
- maximum certificate wall time is `<=420 s`.

If this DEV gate fails, remove `FULL_STACK_CERT_PASS` before the immutable tag and launch the empirical ridge campaign with the full-stack and dose-only analyses descriptive only. Do not create another certificate class, move another anchor, or alter the empirical endpoints. This is the final predetermined stop rule.

---

# 6. Q5 — Use of the constructive DEV result

The finding may appear in paper 2 because:

- it was generated before confirmatory freeze;
- it changes no confirmatory scene, arm, endpoint, or gate;
- it explains why the R17 dose-only class was rejected;
- it defines a clear future problem rather than rescuing a failed outcome.

Mandatory labeling:

> **Development-only constructive analysis; not a confirmatory endpoint or population estimate.**

Report:

- exact DEV image and seed identifiers;
- the two `(nu,b)` cells;
- budget and dose feasibility;
- logdet lower-bound gap and D-efficiency lower bound;
- solver initialization and support restriction;
- the fact that the construction proves existence only.

The current single-scene numbers may be shown as such. If the in-flight replication is included, freeze its set before reading the aggregate: all six DEV images, the already declared seeds, and the same two anchor cells; report every result without selection.

Approved future-work sentence:

> The dose-only construction shows that uniform exposure still permits scene-adapted geometry selection. Optimizing and validating that geometry under relaxed or task-specific conditioning constraints is left to a separate follow-up; no such arm or endpoint is added to the present campaign.

This does not violate R17’s no-redesign rule because it creates no new confirmatory method.

---

# 7. Final refreeze authorization

Commit `c67ba07` is **not** authorized for `m1-freeze`.

A new immutable candidate receives GO only after one outcome-blind ledger shows:

1. `DOSE_SAFE_CERT_PASS` is removed and `FULL_STACK_CERT_PASS` is wired through spec, analyzer, paper skeleton, manifests, and expected-cell ledger;
2. `C_stack` includes budget, ±5% dose, A-risk, and spectral constraints exactly as frozen above;
3. the conic/cutting-plane certificate passes the four-toy suite;
4. the full-dictionary scan and all residual checks are implemented;
5. the 24-cell DEV full-stack feasibility/runtime gate passes, or the categorical certificate is removed under the predetermined fallback branch;
6. the dose-only finding is stored only as no-gate DEV evidence with the support-preserving initialization caveat;
7. the paper’s uniform-dose-collapse wording is removed and the R18 wording inserted;
8. the three empirical operating modes, endpoints, bootstrap rules, and power-for-time disclosures remain unchanged;
9. all 52+972 accounting, dictionaries, SHAs, ridge guards, and confirmatory-access locks still pass;
10. no confirmatory image, pre-scan count, certificate, or endpoint has been opened.

Only then: create one coherent immutable tag and launch all arms together.

## Final scientific position

The R18 finding improves the architecture rather than weakening it:

> **Dose uniformity controls photon placement but not measurement geometry. Scene-adapted geometry retains substantial development-stage information headroom under dose alone; the confirmatory structural question is whether the predeclared conditioning safeguards remove that headroom. Independently, global ridge tracking tests the time-limited, power-available operating-point gain.**

This is narrower, better localized, and more defensible than attributing the collapse to uniform dose itself.