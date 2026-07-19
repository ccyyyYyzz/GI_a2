# R11 — Final ruling: ridge-zone load architecture

**Scope:** this issue amends only the load architecture in R10 issue #2. Every other R10 ruling, including its implementation errata, objective, dictionary, pre-scan accounting, rounding, oracle, guards, fresh-scene rule, endpoints, and manuscript separation, remains in force.

The post-R10 dev result is scientifically material: the original R10 high-flux palette and exact-average-0.60 anchor explored only the foothills of the exact count-information map. The fixed-pattern gains of approximately +6.7 dB at `nu=200` and +12.2 dB at `nu=2000`, with peak locations near the principal ridge, justify adding a ridge-tracking operating-point demonstration and allowing ridge atoms inside OED-DT. The uniform per-pattern servo is rejected; its loss on sparse-support scenes is the expected consequence of erasing the coupling between directional value and pattern brightness.

## Executive ruling

| Amendment item | Frozen ruling |
|---|---|
| (a) RIDGE-FIXED | **Approved as a no-gate arm**, with one naming/guard correction: it is a fixed-support, globally scaled **ridge-tracking policy** whose target changes with `nu`; it is not a per-pattern load servo. |
| (b) OED-DT loads and budget | **Extend the load palette into the ridge zone.** Replace the exact average-load equality by a hard **total incident-signal budget inequality** equivalent to the original mean-0.60 resource. Do **not** use total detected counts as the sole fairness budget; saturation makes that cost degenerate at high load. Report incident and detected budgets jointly. |
| (c) endpoints | The frozen Q90 endpoint, PAVA fitting, censoring taxonomy, and nested bootstrap transfer unchanged. A fast curve whose load tracks `rho_R(nu)` is a policy-level time-to-quality curve. Fixed-dwell results require explicit incident-dose, detected-count, load-tail, and ceiling-count disclosure. |
| (d) water filling | **Confirmed.** Variable-load OED-DT already implements information water filling through its reduced-sensitivity/KKT conditions. No separate uniform-servo arm is permitted in the confirmatory campaign. |
| (e) new guards | Add an exact-ridge target definition, RQL-trust guard, exact ceiling-probability guard, post-rounding cell-level ceiling guard, kernel numerical-certification gate, and mandatory realized ceiling/load disclosure. |

---

# 1. Exact ridge target used by every new arm

Let `J_exact(rho,nu)` be the exact per-slot Fisher information of the integrated nonparalyzable count for `log lambda`, and let `J_q(rho,nu)` be the corresponding quasi-score/Godambe information used to diagnose RQL trust.

For an integer dwell `nu`, define the active-start ceiling probability

\[
 p_{\rm ceil}(\rho,\nu)
 =P_{\rho,\nu}\{N\tau\ge T\}
 =P_{\rho,\nu}\{N\ge\nu\}.
\]

Freeze the scalar admissible set

\[
 \mathcal R_{\rm adm}(\nu)=
 \left\{\rho:\ 0<\rho\le24,\quad
 p_{\rm ceil}(\rho,\nu)\le0.01,\quad
 \frac{J_{\rm exact}(\rho,\nu)}{J_q(\rho,\nu)}\ge0.90
 \right\}.
\]

The production ridge target is

\[
 \boxed{
 \rho_R(\nu)=
 \min\arg\max_{\rho\in\mathcal R_{\rm adm}(\nu)}
 J_{\rm exact}(\rho,\nu).
 }
\]

The outer `min` is the frozen tie rule: if multiple loads are numerically tied, choose the smaller load. The cap `rho<=24` is absolute and is not tuned after dev results. At long dwell it leaves the observed ridge near `rho≈22`; at very short dwell the ceiling/RQL guards may clip the unconstrained information ridge. Every cell records:

```text
rho_ridge_exact_unconstrained
rho_R_production
ridge_clip_reason in {NONE, RHO_CAP, CEILING, RQL_TRUST}
J_exact_at_target
J_q_at_target
p_ceiling_scalar
```

Use the **exact numerical principal ridge**, not the asymptotic expression, to set `rho_R`. The cube-root formula remains the analytic interpretation and initialization.

If `R_adm(nu)` is empty, record `RIDGE_TARGET_FAIL`; no analyst-selected replacement load is allowed.

---

# 2. RIDGE-FIXED arm

## 2.1 Definition

Add the arm named `RIDGE-FIXED`, defined as:

- fixed support family: the frozen strongest fixed family `LBLOB16`;
- the same common 52-pattern pre-scan and 972 main exposures as R10;
- the same reconstruction, analytic TV rule, pooled initializer, solver budget, and scene set as the fixed baselines;
- one **global** source multiplier per dwell;
- no per-pattern normalization or servo.

For each dwell, choose the global multiplier so that the pre-scan-predicted mean main-pattern signal load equals `rho_R(nu)`:

\[
 \frac1{M_{\rm main}}
 \sum_{i=1}^{M_{\rm main}}
 \hat\rho_i
 =\rho_R(\nu),
 \qquad
 \hat\rho_i=\tau\Phi_R a_i^T\hat x.
\]

When the fixed translate family has exact unit mean under the frozen normalization, this calibration is scene-independent; otherwise it uses only the common pre-scan estimate. It never uses truth.

The load **changes with dwell**, so the precise interpretation is a fixed-support, ridge-tracking operating policy. Retain the requested label `RIDGE-FIXED` for code/table continuity, but define it this way in the paper.

## 2.2 Cell-level guard without destroying self-water-filling

After the single global scaling, evaluate every predicted main-pattern load. Do not equalize them. If necessary, reduce the **single global multiplier** by deterministic bisection until

\[
 \frac1{M_{\rm main}}\sum_i p_{\rm ceil}(\hat\rho_i,\nu)\le0.01,
 \qquad
 \max_i p_{\rm ceil}(\hat\rho_i,\nu)\le0.05,
\]

and

\[
 \min_i\frac{J_{\rm exact}(\hat\rho_i,\nu)}{J_q(\hat\rho_i,\nu)}\ge0.90.
\]

This is a global safety clip, not a servo. Record the originally requested target and the achieved mean load. If the achieved mean is below `0.90*rho_R`, mark `RIDGE_GUARD_CLIPPED`; the result remains descriptive.

## 2.3 Status and comparisons

`RIDGE-FIXED` carries **no confirmatory gate** and cannot rescue any R10 primary endpoint. Its mandatory comparisons are:

1. `LBLOB16 @ rho=0.60` — isolates the global operating-point change;
2. `SCAT16 @ rho=0.60` — connects to the original fixed sparse baseline;
3. safe `rho=0.05` — permits a descriptive Q90 speed curve.

At every dwell report:

- PSNR_rad and the frozen secondary metrics;
- predicted and realized incident signal counts;
- predicted and realized detected counts;
- dose ratio relative to `rho=0.60` and `rho=0.05`;
- source/peak irradiance ratio;
- `rho_5`, `rho_50`, `rho_95`, and `rho_max`;
- predicted and realized ceiling-count fractions;
- mean exact information `mean_i J_exact(rho_i,nu)`.

The equal-dwell gain is an **operating-power-for-time** result, not an equal-photon-efficiency result.

---

# 3. OED-DT variable-load palette

For each structured base direction and dwell, form candidate atoms at the frozen load palette

\[
 \boxed{
 \mathcal L_\nu=
 \operatorname{unique}\left\{
 0.30,\ 0.60,\ 1.00,\
 \min(3.00,\rho_R),\
 \tfrac12\rho_R,\
 \rho_R
 \right\}.
 }
\]

Values are sorted and deduplicated to relative tolerance `1e-10`. A candidate is removed if it violates:

- the R10 finite gain range;
- the R10 zero-support rule;
- the R10 peak-irradiance or weight-dispersion guards;
- `rho<=24`;
- `p_ceiling(rho,nu)<=0.05`;
- `J_exact/J_q>=0.90`.

The half-ridge level is mandatory. The dev sweep at `nu=200` peaked below the scalar ridge, and a palette jumping directly from `3` to approximately `10` or `22` is unnecessarily coarse. Zero design mass already represents starving/omitting an atom; no explicit zero-load exposure is required.

For each base direction `b` and target load `r`, its amplitude is chosen from the pre-scan estimate so that

\[
 \tau\lambda(b,r;\hat x)=r,
\]

subject to the frozen gain and zero-support guards. The exact information atom remains the R10 chain-rule atom

\[
 H(b,r;\hat x)=
 \nu J_{\rm exact}(r,\nu)
 \frac{\Phi^2}{\lambda(b,r)^2}
 B^Taa^TB.
\]

---

# 4. Resource architecture: incident budget is primary, detected counts are dual-reported

## 4.1 Reject a detected-count-only budget

Do **not** replace the R10 average-load constraint by total detected counts alone. For a nonparalyzable detector,

\[
 E[N\mid\rho]\simeq\frac{\nu\rho}{1+\rho},
\]

whose marginal cost tends to zero at high load. A detected-count-only constraint can therefore permit arbitrarily large incident power while appearing nearly cost-free. That is not a defensible design-comparison resource.

## 4.2 Gate-carrying OED-DT budget

Replace the exact equality `mean rho = 0.60` by the hard **incident-signal budget inequality**

\[
 \boxed{
 \sum_{m=1}^{972}\nu\rho_{s,m}
 \le972\,\nu\,(0.60),
 }
\]

where `rho_s` excludes calibrated dark/background load. In normalized design-measure notation this is

\[
 \int\rho_s(a)\,d\xi(a)\le0.60.
\]

This preserves the original R10 total source/dose resource while allowing a subset of high-value atoms to enter the ridge zone and forcing compensation by low-load or omitted atoms. It is the correct budget for testing information water filling. The optimizer need not spend the entire budget if all remaining feasible photons would fall beyond the useful information region or violate a guard.

The R10 per-pixel cumulative-dose bound, total-dose accounting, peak irradiance bound, and 52+972 exposure accounting remain unchanged.

The statement that an average-0.60 **budget** forbids ridge atoms is incorrect; only a pointwise/uniform load constraint forbids them. The variable-load measure can allocate positive mass at `rho_R` while respecting the same total incident resource.

## 4.3 Mandatory dual reporting

For every arm report both

\[
 B_{\rm inc}=\sum_i\lambda_{s,i}T
\]

and

\[
 B_{\rm det}^{\rm pred}=\sum_iE[N_i],
 \qquad
 B_{\rm det}^{\rm obs}=\sum_iN_i.
\]

Also report total and maximum per-pixel dose. `B_det` may be used as an electronics-throughput guard if a hardware limit is frozen before confirmation, but it is never the sole fairness constraint.

`RIDGE-FIXED` deliberately uses a larger incident budget than the mean-0.60 arms and is therefore no-gate. OED-DT retains the original resource budget and remains the gate-carrying method.

No additional ridge-budget OED gate is added in R11; that would mix operating-point gain with design-allocation gain.

---

# 5. Q90 endpoint and censoring at ridge loads

The existing endpoint transfers unchanged:

- same safe curve at `rho=0.05`;
- same Q90 target derived from the frozen safe curve;
- same equal-weight PAVA fit in log optical time;
- same crossing interpolation;
- same censoring taxonomy;
- same nested scene/seed bootstrap and frozen gate constants for gate-carrying OED-DT.

For OED-DT, the fast arm is a fixed resource-constrained design policy whose per-pattern loads vary. For `RIDGE-FIXED`, the fast arm is a dwell-dependent policy `rho_R(nu)`. A policy-level time-to-quality curve is valid, but the manuscript must state that the source operating point is retuned as dwell changes.

The terminal comparison remains named

> **fixed-dwell gain**

or

> **fixed optical-integration-time gain**.

It must never be called equal-photon-budget gain when comparing `rho=0.05`, `0.60`, or ridge operation.

## Intention-to-treat handling

No cell is deleted after observing ceiling counts. If a gate-carrying RQL reconstruction fails numerically or becomes non-identifiable because of realized ceiling events, use the frozen failure/censoring path (`S_gate=0` where applicable) and retain the cell. `RIDGE-FIXED` is no-gate but must still report the failure.

---

# 6. Information water filling is already inside OED-DT

No separate servo arm is added.

For a base direction `b`, let

\[
 \ell_V(b)=q(b)^TV^{-1}q(b)
\]

be its current local directional value. In the zero-dark idealization, gain changes the scalar load but not `q(b)`, so under an incident-budget multiplier `beta` the reduced utility of assigning load `rho` is

\[
 U_b(\rho)=
 \nu J_{\rm exact}(\rho,\nu)\ell_V(b)
 -\beta\nu\rho
 -\text{other active resource costs}.
\]

An interior load satisfies

\[
 \boxed{
 \ell_V(b)\,J_{\rm exact}'(\rho,\nu)
 =\beta+\text{marginal guard costs}.
 }
\]

Thus high-value directions receive larger loads, generally up to but not beyond the useful ridge region; low-value or blank directions receive low load or zero design mass. This is information water filling. Uniform per-pattern load is optimal only in the special case that directional values and marginal resource costs are equal.

With dark counts, gain caps, per-pixel dose multipliers, and nonzero pre-scan uncertainty, the scalar closed form changes, but the R10 constrained LMO/KKT system performs the same allocation numerically.

The refuted uniform servo may be shown as dev-level motivation in a methods supplement, but it is not a confirmatory baseline and cannot be used to define the final method.

---

# 7. New ridge-zone guards

All R10 guards remain. Add the following hard pre-confirmation guards.

## G7 — exact kernel numerical certification

For every frozen `(nu,rho)` node used by `RIDGE-FIXED` or the OED palette:

1. compute the exact PMF in log space through the ceiling count;
2. require absolute PMF mass error `<=1e-10`;
3. require absolute expected-score error `<=1e-9`;
4. require relative agreement between analytic score-sum Fisher information and an independent derivative check `<=1e-6`;
5. require finite, nonnegative `J_exact` and `J_q`;
6. verify the principal maximum by a dense log-load bracket followed by deterministic one-dimensional refinement;
7. if two maxima differ by at most `1e-10` in `J`, choose the lower load.

Failure gives `KERNEL_CERT_FAIL`; the affected load cannot be silently replaced.

## G8 — predicted ceiling/RQL admissibility

For an OED atom:

\[
 p_{\rm ceil}(\rho,\nu)\le0.05,
 \qquad
 J_{\rm exact}/J_q\ge0.90.
\]

After exact rounding, require the design-weighted ceiling probability

\[
 \sum_jw_jp_{\rm ceil}(\rho_j,\nu)\le0.01.
\]

Failure gives `RQL_LOAD_GUARD_FAIL` before confirmatory data are generated.

For `RIDGE-FIXED`, use the global-reduction procedure in Section 2 rather than per-pattern equalization.

## G9 — realized stress disclosure

Every cell records

```text
ceiling_count_fraction_observed
ceiling_count_fraction_predicted
rho_5, rho_50, rho_95, rho_max
incident_budget
predicted_detected_budget
observed_detected_budget
source_gain_relative_to_rho_0p6
source_gain_relative_to_safe
```

If the realized ceiling fraction exceeds `0.05`, mark `RQL_CEILING_STRESS`. The cell remains in all frozen analyses; the label is disclosure, not a post hoc exclusion rule.

## G10 — budget and gain integrity

After continuous optimization and after exact rounding separately, verify:

- incident budget `<=972*nu*0.60` for gate-carrying OED-DT;
- all R10 per-pixel dose and peak constraints;
- all candidate gains inside the frozen finite range;
- no nonzero load assigned to a zero-support candidate;
- exact rounded D-efficiency and certificate conditions from R10.

Any violation is a design failure, not an invitation to retune the load palette.

---

# 8. Frozen campaign amendment

The R10 campaign is amended only as follows:

1. Add no-gate `RIDGE-FIXED` at the exact admissible ridge target for every dwell.
2. Extend OED-DT atoms using the frozen palette in Section 3.
3. Change OED-DT's exact mean-0.60 equality to the incident-budget inequality in Section 4.
4. Retain OED-EQLOAD at pointwise `rho=0.60` as the geometry-only ablation.
5. Keep OED-DT as the only gate-carrying designed arm.
6. Keep the common 52-pattern pre-scan and 972 main-exposure accounting exactly as corrected in R10 erratum E2.
7. Keep all R10 endpoints and gates unchanged.
8. Add G7-G10 and the mandatory dual resource disclosures.
9. Do not add a uniform-servo confirmatory arm.
10. Do not reopen the OE manuscript or the R9 scope lock.

## Frozen interpretation

The follow-up method line now separates two questions cleanly:

- **Operating-point question:** how much is gained by globally moving a fixed sensing geometry from `rho=0.60` toward the exact finite-window ridge? Answered descriptively by `RIDGE-FIXED`.
- **Design-allocation question:** under the original total incident resource, can a dead-time-aware optimizer spend photons selectively on high-value directions near the ridge and starve low-value directions? Answered confirmatorily by budgeted variable-load OED-DT.

This separation is mandatory. The +6 to +12 dB dev result is strong evidence for the operating-point arm, but it must not be presented as equal-dose proof of the OED allocation method.