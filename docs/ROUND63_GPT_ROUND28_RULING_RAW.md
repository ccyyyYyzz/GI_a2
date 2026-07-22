# R28 ruling (GitHub issue #20, raw)

Title: R28 ruling: finite-library Gate A and dose-relaxed FW diagnostics
Posted: 2026-07-22T06:44:55Z

---

# R28 — Fast ruling: finite-library reachability gate and dose-relaxed FW diagnostics

**Reference commit:** `fa568ca`  
**Source:** `docs/ROUND63_GPT_ROUND28_QUESTION.md`

**Scope:** pre-grid amendment to the DEV bridge only. No M1 result, paper, or confirmatory campaign is changed.

> **R28 VERDICT: adopt option (c), with two precision amendments.**
>
> 1. Gate A becomes a **finite-library image-reachability gate**, not a claim about the entire dose-feasible design class.
> 2. The FW oracle arms relax **only the ±5% per-pixel dose band**. They retain the same row count, incident-budget cap, peak/load/admission constraints, `k_eff >= 32` physical dictionary, dwell, and reconstruction. They are descriptive, nondeployable surrogate diagnostics and no gate reads them.
>
> The amendment is justified because the original TRUE-X-FW gate no longer tests its intended proposition inside the frozen deployment class: the FW objective’s useful direction is removed by dose uniformization, while an exactly materialized in-class library design already demonstrates image-level headroom. Allowing the broken oracle to kill the method line would be a category error.

---

## 1. Frozen replacement for Gate A

At the hard bridge cell

```text
c = 0.05
nu = 2000
12 frozen misalignment-stress scenes
5 paired seeds
```

define the fixed composite comparator

\[
Q_{\mathrm{base},j}
=\max\{Q_{\mathrm{SCAT32-060},j},Q_{\mathrm{RIDGE-SCAT32},j}\},
\]

where each `Q` is the five-seed mean radiometric PSNR.

For each scene, define one finite-library oracle bank by

\[
k_j^*=\min\arg\max_{k\in\{0,\ldots,7\}}
\overline{Q}_{j,k},
\]

with the bar denoting the mean over the five frozen paired seeds and the smallest bank index breaking an exact tie. The oracle value is

\[
Q_{\mathrm{ORACLE-LIB},j}=\overline{Q}_{j,k_j^*},
\qquad
\Delta Q^{A}_j=Q_{\mathrm{ORACLE-LIB},j}-Q_{\mathrm{base},j}.
\]

Every candidate bank must be its genuine, exact-972, postchecked materialization under the common deployment guards. No bank may be selected separately for each noise seed.

### `LIBRARY_REACHABILITY_PASS`

PASS iff all three original numerical conditions hold:

1. \(\operatorname{median}_{j=1}^{12}\Delta Q^A_j\ge1.0\) dB;
2. at least `9/12` stress scenes have \(\Delta Q^A_j>0\);
3. the frozen 90% scene-bootstrap lower bound on the median exceeds `0.50 dB`.

**Decision consequence:** if this gate fails, the declared finite bank does not contain a robust image-rescuing design and RLMI/M2 stop. It does **not** claim that no design exists elsewhere in the full physical class.

Do not call ORACLE-LIB “the class ceiling.” Use:

> **finite-library image oracle**  
> or  
> **declared-library reachability reference**.

---

## 2. Gate B becomes the coherent capture gate

All existing RLMI image-gain and no-harm conditions remain unchanged. The capture condition is frozen relative to the new Gate-A reference:

\[
\operatorname{median}\Delta Q_{\mathrm{RLMI}}
\ge0.60\,
\operatorname{median}\Delta Q_{\mathrm{ORACLE-LIB}}.
\]

If any implementation still uses TRUE-X-FW in this denominator, replace it before the grid. Gate A asks whether the declared library contains useful geometry; Gate B asks whether the deployable allocator captures enough of that achievable library gain.

Gate C, `PLUGIN_LATENCY_PASS`, the baseline definition, scene cohort, thresholds, and all other bridge constants remain frozen.

---

## 3. FW arms: descriptive dose-relaxed diagnostics only

Run TRUE-X-FW and XHAT-FW with:

- the same `k_eff >= 32` super-atom dictionary;
- exactly 972 main rows;
- the same total incident-budget cap;
- the same row-wise peak, detector-load, ceiling, nonnegativity, and admission guards;
- the same dwell and RQL reconstruction;
- **no ±5% cumulative per-pixel dose-band constraint**.

Record the full realized per-pixel dose profile, maximum deviation, source-power distribution, budget, peak, achieved loads, and materialization status. Label both arms:

```text
NONDEPLOYABLE_DOSE_RELAXED_FW_DIAGNOSTIC
```

No Gate A–C statistic may read these arms. A dose-relaxed FW failure or success cannot kill or rescue RLMI.

Do not call these arms a physical-class ceiling. They are a diagnostic of what the local Fisher/A-risk objective asks for when the dose-uniformity restriction is removed.

---

## 4. Revised decomposition

Retire the old interpretation in which every difference was called a Fisher-to-image or finite-bank gap. Report four quantities separately:

1. **in-class finite-library reachability**  
   `ORACLE-LIB − fixed composite baseline` — Gate A;
2. **deployable allocation loss**  
   `RLMI − ORACLE-LIB` — Gate B/capture analysis;
3. **plug-in FW loss**  
   `XHAT-FW − TRUE-X-FW`, both dose-relaxed and descriptive;
4. **dose-relaxation/library contrast**  
   `TRUE-X-FW − ORACLE-LIB`, descriptive only.

Quantity 4 conflates the removal of the dose band, the FW surrogate, and restriction to the finite bank. It must not be labeled a pure “finite-bank approximation gap” or a pure “Fisher-to-image gap.”

The direct descriptive surrogate-to-image quantity `TRUE-X-FW − fixed baseline` may also be shown, with the same nondeployable label.

---

## 5. Smoke-exposed scenes

**Do not replace `twopop_0`, `control_0`, or any other frozen scene.** Replacement after seeing smoke values would be less defensible than retaining them.

Freeze the following provenance:

```text
SMOKE_EXPOSED_PREGRID = {twopop_0, control_0}
```

- retain them in every originally assigned DEV calculation;
- do not change any threshold, family count, or bootstrap rule;
- report a no-gate leave-smoke-exposed-out sensitivity after the full grid;
- the official gates continue to use the complete frozen cohort.

This is a DEV bridge, not M2 confirmatory data. The exposure is acceptable because it was a materialization/architecture diagnostic, no cohort member was replaced, and the amendment is fully recorded before the grid.

---

## 6. Ratification of the two implementation disclosures

### 6.1 Frozen nominal-power witnesses

**Ratified.** A bank may use a manifest-frozen per-atom nominal power schedule, including the committed Sinkhorn amplitudes in `[1/4,4]` under the peak cap, provided that the schedule is part of the bank identity rather than a post-hoc feasibility trim.

The identical physical powers must be used in:

- exact-972 admission;
- dose, budget, and peak accounting;
- every scenario Gramian;
- materialization;
- acquisition simulation;
- final reconstruction metadata.

Strict equal row power is not required. The manuscript/method description must say **pattern-and-power bank**, not imply that only binary geometry varies. L0’s frozen ridge schedule remains an explicit special member of the same convention.

### 6.2 `k_eff >= 32` online dictionary

**Ratified.** The restriction is physically conservative and remains in force for both the deployed library and the dose-relaxed FW diagnostics. It does not by itself make the FW arm deployable; that status is determined by the dose-band relaxation above.

---

## 7. Launch authorization

The T-B bank manifest does **not** need regeneration: the eight bank identities and exact-admission witnesses are unchanged. Regenerate only the bridge specification/analyzer provenance needed to encode this R28 amendment.

Before launch:

1. replace Gate A in code and tests by `LIBRARY_REACHABILITY_PASS` exactly as above;
2. verify Gate B’s 60% denominator is ORACLE-LIB;
3. remove the dose band only from the two FW diagnostic arms and add the mandatory label/disclosures;
4. update the four-gap names;
5. add the smoke-exposure provenance flag and leave-out sensitivity output;
6. rerun the unit suite and one outcome-blind smoke confirming genuine ORACLE-LIB/RLMI materialization and that no gate reads an FW arm.

After those checks pass, the held bridge grid is authorized to launch. No further scene, bank, threshold, or endpoint redesign is authorized by this ruling.