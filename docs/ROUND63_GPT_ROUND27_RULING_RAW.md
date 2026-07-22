# R27 ruling (GitHub issue #19, raw)

Title: R27 ruling: dose-feasible bank rebuild and latency claim split
Posted: 2026-07-22T04:59:03Z

---

# R27 — Fast ruling: exact-dose bank rebuild and latency-claim split

**Reference commit:** `57f0abf`  
**Source:** `docs/ROUND63_GPT_ROUND27_QUESTION.md`

**Scope:** pre-grid amendment to the DEV bridge only. No confirmatory cell has run. T-A scenes and Gates A–C remain frozen except for the bank identities forced by the physical-feasibility correction below.

> **R27 VERDICT.**
>
> 1. **Collision 1 is a real specification defect, not a rounding-algorithm failure.** Keep the common ±5% dose band and the exact 972-row budget. Rebuild the RLMI dictionary so every *physical* atom has effective support at least 32. A paired pair of disjoint 16-pixel supports is permitted only when it is defined and optimized from the outset as one simultaneous 32-pixel super-atom. Do not relax the guard and do not keep trying cube/dependent rounding against an arithmetic obstruction.
> 2. **Collision 2 does not kill the image-science bridge or M2.** It kills the subsecond “plugin/real-time” claim. Gates A–C remain the scientific go/no-go conditions; the original Gate D becomes a claim-classification gate named `PLUGIN_LATENCY_PASS`. A method that passes A–C but fails D may proceed only as a **two-shot batch adaptive illumination method**, with the full latency distribution disclosed.
> 3. Output-equivalent caching/vectorization is authorized before the grid. Reducing `r=200` or `S=16` is not authorized for the production bridge merely to hit the latency target.

---

## 1. Collision 1 — ruling on the dose band

### 1.1 Adopt option (b), with option (d) allowed only as a bank-construction primitive

The observed failure is combinatorial. At `M_main=972`, a uniform `k=16` design has mean pixel coverage

\[
972\,16/1024=15.1875.
\]

Even the best equal-weight integer coverage must use 15 and 16 hits; the 16-hit level is about `+5.35%` above the mean. Thus the ±5% band can be infeasible before any optimizer or rounding heuristic is considered. Per-row trims cannot generally repair this because each trim moves all illuminated pixels together.

The amendment is therefore:

\[
\boxed{k_{\rm eff}\ge 32\quad\text{for every physical main-acquisition atom admitted to RLMI.}}
\]

- Replace scattered-16, compact-16, ring-16, and other 16-support mixture entries by support-32 counterparts.
- Rebuild multiscale entries with minimum physical support 32.
- Re-solve `L7` over a dictionary already restricted to physically admissible support-32-or-larger atoms.
- `L0` remains the existing exact SCAT32 design.

A **paired-super-atom** is permitted to preserve some 16-scale spatial structure:

\[
a^{\rm pair}=a^{(1)}+a^{(2)},
\qquad
\operatorname{supp}(a^{(1)})\cap\operatorname{supp}(a^{(2)})=\varnothing,
\qquad
|\operatorname{supp}(a^{\rm pair})|=32.
\]

It is one simultaneous physical exposure, one dwell, and one information atom—not two temporally paired rows. The pair must be fixed jointly in the bank manifest and scored with its actual 32-pixel load, dose, peak power, and Gramian. Post-allocation pairing of arbitrary 16-pixel rows is forbidden. Weighted/FW atoms are not paired after the fact; `L7` must be reoptimized on the physical super-atom dictionary.

This is option **(b)** plus a disciplined version of **(d)**. It is not a waiver for k=16 rows.

### 1.2 Reject options (a) and (c)

- **Cube/dependent rounding:** may improve a feasible integer problem but cannot defeat the one-hit granularity floor. It remains an implementation option after the bank has passed the exact-feasibility admission test; it is not the remedy.
- **Per-bank dose-band relaxation:** rejected. A bank-specific ±7% rule would weaken the common safety/resource class to admit a desired geometry, would not rescue banks whose structured floors are much larger than 7%, and would make the simplex guards noncommon.

The exact main row count remains `972`; do not redesign the 52/972 accounting in this amendment.

---

## 2. Exact replacement for R25 §9

Replace the opening of the materializer specification by the following frozen text.

> **Bank admission before allocation.** Continuous dose feasibility is necessary but not sufficient. Every library measure must expose physical atoms with effective support `k_eff >= 32` and must demonstrate, before the bridge grid, at least one exact 972-row realization satisfying the common incident budget, ±5% per-pixel dose band, row-wise peak/admission constraints, and any frozen family restrictions. A 16-pixel motif may enter only inside a manifest-frozen, simultaneously projected, disjoint-pair super-atom of effective support 32. A bank that fails this exact-count admission test is excluded from the simplex rather than rescued by a bank-specific guard.
>
> **Exact mixture materialization.** Merge identical physical atoms and form the union target measure. Solve for nonnegative integer multiplicities `n_a` with `sum_a n_a=972` while enforcing the common guards. Dependent rounding, exchange, or mixed-integer repair may be used only inside this feasible support class. The row-power trim `gamma in [0.94,1.06]` is a calibration/polish variable, not a certificate of feasibility: the integer support/multiplicity design must first satisfy the dose band under the bank’s frozen nominal powers, after which any trim is applied and every guard and scenario matrix is recomputed.
>
> **Failure.** If no compliant exact design is obtained, return `L0` and emit `MIXTURE_MATERIALIZATION_FAIL`; do not loosen the dose band, row count, or peak guard.

The post-materialization requirements from R25 remain unchanged: recompute actual dose/budget/peak guards, every scenario matrix, realized standardized-maximin regret, and materialization loss.

### 2.1 Manifest consequence

**Yes: regenerate the complete T-B manifest before the grid.** The bank identities, physical rows, source normalization, expected Gramians, exact-feasibility witnesses, hashes, and ORACLE-LIB labels all change. Required pre-grid sequence:

1. rebuild `L1–L7` under `k_eff>=32`;
2. run exact-972 bank admission for every entry;
3. regenerate bank and union-atom hashes;
4. rebuild all cached `F_{k,s}` machinery;
5. rerun the outcome-blind T-B/T-C unit suite and phase-1.5 smoke;
6. only then launch the frozen bridge grid.

T-A scenes and image-level gate constants do not change.

---

## 3. Collision 2 — latency decision tree

### 3.1 Clarification of R24 Gate D

The original sentence—“if this fails, the plugin claim fails even if PSNR improves”—is controlling. It did **not** logically imply that a scientifically effective two-shot method must be killed.

Freeze the split as follows:

- **Gates A–C:** method-science go/no-go. If any fails, no M2.
- **Gate D, renamed `PLUGIN_LATENCY_PASS`:** claim-classification gate:
  - PASS: median route latency `<1.0 s` and 95th percentile `<2.0 s`; the method may be described as subsecond/plugin-like.
  - FAIL: M2 may still launch if A–C pass, but the method must be described as **batch/two-shot adaptive illumination**, not real-time, immediate, or subsecond plugin operation. Report median, 95th percentile, maximum, CPU, threading, and cache state.

There is no substitute latency threshold. The measured `~9.6 s` is neither hidden nor converted into a new pass bar.

### 3.2 Pre-grid engineering that is authorized

The following are output-equivalent implementation changes and are authorized before the grid:

- cache `H_{0,s}`, `W_s`, and every per-bank scenario contribution `F_{k,s}` for the identical scene/seed/physics condition;
- cache the per-scenario simplex-oracle values `R_s^*`;
- reuse those caches across allocator calls and bridge arms that share the identical pre-scan state;
- batched/parallel Cholesky or triangular solves;
- deterministic warm starts;
- exact Woodbury/low-rank updates where they reproduce the same matrix objective;
- compiled/vectorized kernels.

Before use, the optimized path must agree with the reference implementation on all phase-1.5 cells within:

```text
|t_fast - t_ref| <= 1e-8
||w_fast - w_ref||_1 <= 1e-6
KKT residual <= the existing frozen tolerance
identical exact 972-row realization after the lexicographic tie-break
```

If the optimum is numerically nonunique, the frozen closest-to-`e0` tie-break is applied before comparing weights.

### 3.3 Changes that are not authorized as engineering

Do not reduce production

```text
r = 200
S = 16
```

to chase Gate D. Both alter the declared estimator-risk problem; `S=16` already provides only finite-scenario robustness, and a smaller task subspace can change bank ordering. `S=8`, reduced-`r`, mean-risk, or other fast variants may be reported later as no-gate latency/accuracy ablations, but they cannot replace the bridge method without a new pre-grid method amendment and rerun of all relevant DEV gates.

---

## 4. Jitter-cap smoke result

No corrective action is authorized for the observed `RIDGE-SCAT32` degradation at `c=0.05`. The `c=0` ridge load is intentionally wrong for a jittered detector; the composite fixed comparator protects the bridge decision exactly as specified. Preserve this result as physics, not as a calibration bug.

---

## Launch authorization

The grid remains **HOLD** until:

1. the support-32-or-larger T-B library and manifests are regenerated;
2. every bank passes exact-972 admission or is excluded before allocation;
3. the optimized/cached allocator passes the equivalence tests above;
4. the phase-1.5 smoke reruns cleanly.

After those conditions, the bridge grid may launch. A remaining `PLUGIN_LATENCY_PASS=FALSE` does not block the launch; it irrevocably removes the subsecond/plugin claim unless a later output-equivalent implementation actually passes the original timing bar.