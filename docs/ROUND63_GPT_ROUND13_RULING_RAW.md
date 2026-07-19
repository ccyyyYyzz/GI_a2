# R13 — Final M1 freeze audit

**Audit target:** `docs/ROUND63_METHOD_SPEC_M1.md` and `code/round63/oed_design_v3.py` at commit `6bdfc3c`.

**Scope:** this ruling governs only the follow-up M1 method campaign. It does not reopen the frozen OE manuscript, Study-1/2 outcomes, or the R10–R12 claim discipline.

> **FREEZE VERDICT: HOLD / NO-GO.**
>
> The scientific architecture is coherent, and the photon-budgeted-versus-time-budgeted framing is approved. However, the current M1 draft and v3 implementation are not yet freezeable. Several deviations are not cosmetic: the pre-scan called “balanced” is strongly dose-unbalanced and its information calculation uses an unimplementable row-by-row oracle normalization; the implemented atom dictionary is smaller than the R10 minimum dictionary; `RIDGE-FIXED` still emits 1024 main rows rather than 972; the current G5 mixture is not materialized as an exact design; and the advertised KW/rounding/A-risk guards are not yet computed for the final exact design.
>
> An immutable `m1-freeze` commit is authorized only after every hard item in Section 10 below passes outcome-blind self-tests.

## Executive ruling

| Item | R13 ruling |
|---|---|
| A1 — trust-ratio direction | **Direction approved, implementation amended.** The efficiency direction must be approximate-information / exact-Fisher, not the vacuous reciprocal. Replace the present `J_q` terminology and add a finite-window RQL bias guard; details below. |
| A2 — peak guard | **Shape-only as the sole guard is rejected.** R11 necessarily supersedes the old `a_j<=64` load-scaled cap, but M1 must retain both a shape-concentration guard and an explicit finite physical-peak guard. |
| A3 — `M_main=972` | **Approved.** This is the R10 erratum and must hold for every arm. The present `ridge_fixed_design` violates it and must be fixed. |
| A4 — ceiling probability | **Approved with scope.** For integer `nu` and the active-start model, `p_ceil=P(N>=nu)=P(Pois(rho)>=nu)` is exact. Do not reuse this identity for noninteger windows or other start conventions. |
| A5 — dose-cap exclusion | **Mechanism approved, “theorem at 9.6” rejected.** There is an exact row-specific necessary bound; `9.6` is the p=0, k=16, mean-dose≈0.60 specialization, not a universal dictionary theorem. |
| G5 repair / mixture | **Not freezeable as implemented.** Best-so-far repair plus strict-descent polish is an acceptable deterministic heuristic, but `COMPLIANT-VIA-MIXTURE` is valid only after an exact 972-row mixed design is actually constructed and all guards are recomputed. A continuous dose-vector mixture is not an executable design. |
| Certificate | **Current convention not approved as a KW certificate.** It can be converted into a rigorous conservative relaxed upper bound by using the additive dual-gap formula with budget slack. The current ratio gap and `tol=0.1` do not satisfy the frozen R10 certificate guard. |
| §12(a) family mapping | **Approved with naming amendment.** Reuse the six generators; call the first family `chirp/line-pair` only to the extent implemented, and do not label it USAF if no USAF target is generated. |
| §12(b) per-arm safe reference | **Approved after exact policy definition.** Fixed arms reuse the same main pattern multiset at safe and fast loads; OED-DT has a separately frozen safe-load OED solve and fast-load OED solve. No arm borrows another arm’s safe curve. |
| §12(c) 52-row pre-scan | **Rejected in current form.** Replace it with the balanced, scene-independent construction in Section 2. |
| §12(d) analyzer naming | **Approved after restoring the omitted design-gain secondary and exact source-of-record names.** |
| §12(e) ridge values | **Approved for `nu=200,2000` as reported.** They do not certify the remaining seven dwell values. Full-grid G7/trust/ceiling certification is mandatory before freeze. |
| §12(f) photon economics | **Approved with a nonuniversality caveat.** OED-DT and RIDGE-FIXED are conjugate photon-budgeted and time-budgeted corners. Absence of ridge atoms in one self-test is a valid KKT outcome, not a universal theorem that photon-budgeted OED never selects them. |

---

# 1. Ambiguity ledger A1–A5

## A1 — RQL trust: correct direction, corrected quantity

The literal R11 ratio `J_exact/J_q>=0.90` is vacuous if `J_q` is a valid Godambe information, because Godambe information cannot exceed Fisher information. The intended direction is therefore

\[
\eta_{\rm RQL}=\frac{J_{\rm surrogate}}{J_{\rm exact}}\ge0.90.
\]

However, the current `J_q` calculation is not a fully rigorous finite-window Godambe information of the production RQL score: the active-start RQL score is only asymptotically centered. Do not retain the statement “Godambe information of any estimating equation is below Fisher” while using an uncentered finite-window score.

Freeze the following two diagnostics instead.

### A1.1 Exact mean-information efficiency

Using the exact PMF, compute

\[
J_{\rm mean}(\rho,\nu)
=\frac{\{\partial_{\log\rho}\,\mathbb E_\rho N\}^{2}}
       {\nu\,\operatorname{Var}_\rho N}.
\]

By the score/Cauchy–Schwarz information inequality,

\[
0\le \eta_{\rm mean}
=\frac{J_{\rm mean}}{J_{\rm exact}}\le1.
\]

Use

\[
\boxed{\eta_{\rm mean}\ge0.90}
\]

as the information-efficiency part of the trust guard.

### A1.2 Finite-window RQL location bias

For repeated scalar counts from true load `rho`, the population root of the production RQL moment equation is

\[
\rho^\dagger(\rho,\nu)
=\frac{\mathbb E_\rho N}{\nu-\mathbb E_\rho N}.
\]

Record

\[
b_{\log\rho}=\left|\log\frac{\rho^\dagger}{\rho}\right|
\]

and require

\[
\boxed{b_{\log\rho}\le0.05.}
\]

Rename the fields to `J_mean`, `mean_info_efficiency`, and `rql_logload_bias`; do not call the current quantity exact Godambe information. Recompute every ridge/palette trust decision after this amendment. The reported long-window ridge locations should remain unchanged if the claimed ≥98.5% efficiency is genuine.

## A2 — split the peak guard into shape and physical components

R10’s original load-scaled cap `max_j a_j<=64` is incompatible with the later R11 ridge architecture. That conflict must be recorded as an explicit R11 supersession, not hidden as an “interpretation.”

Freeze two guards:

### G3-SHAPE

For a support of size `k`, positive weights have support mean one and satisfy

\[
\frac14\le w_j\le4,
\qquad k\ge16.
\]

The unit-load shape peak is therefore bounded by

\[
\max_j a^{\rm unit}_j\le\frac nk\,4.
\]

The implementation must use the actual support size `k`; the current hard-coded `16` is not valid for the required `k=32` dictionary families.

### G3-PHYSICAL

In the Study-2 relative-irradiance convention, freeze

\[
\boxed{\max_{i,j} a^{\rm physical}_{ij}\le1536.}
\]

This is the peak of a uniform `k=16` atom at the approved absolute load cap `rho=24` (`64×24`). It permits the `p=0` ridge arm but prevents matched-intensity ridge atoms from silently demanding four times more peak power. Every cell must additionally report its realized peak ratio.

If a different hardware peak cap is desired, it must be set now from a declared source specification; it cannot be removed entirely. A finite rho cap plus a shape cap establishes compactness, but it does not by itself preserve the physical-resource claim made in R10.

## A3 — exposure accounting

Approved:

\[
M_{\rm total}=1024=52+972.
\]

This must be asserted at the returned-pattern-matrix level for every arm. The current `ridge_fixed_design` constructs all `n=1024` translations as main rows and therefore violates A3. It must return exactly 972 main rows, with a frozen deterministic row multiset and the same final guard checks as other fixed arms.

## A4 — ceiling probability

For integer `nu=T/tau` under the active-start convention,

\[
P(N\ge\nu)
=P\{\Gamma(\nu,\rho)\le1\}
=P\{\operatorname{Pois}(\rho)\ge\nu\}.
\]

The implementation is correct. Record that its scope is integer `nu`, active start, ideal nonparalyzable counting. The tiny probabilities at `nu=200,2000` do not imply tiny probabilities at `nu=5,10`; all dwell values require independent certification.

## A5 — exact dose-cap mechanism, amended wording

For a final exact row `r` with load `rho_r` and unit-load direction `g_r`, define

\[
d_j=\frac1M\sum_s\rho_s g_{sj},
\qquad
\bar d=\frac1n\sum_j d_j.
\]

G5 implies `d_j<=1.05 dbar`. Because every contribution is nonnegative, a necessary condition for row `r` is

\[
\boxed{
\rho_r\le
\frac{1.05\,M\,\bar d}{\max_j g_{rj}}.
}
\]

For a uniform `k=16` row with `max g=64`, `M=972`, and `dbar≈0.60`, this gives approximately `9.57`. For matched weights, bright supports, other `k`, or a different realized mean dose, the bound changes. Freeze the label

```text
A5_ROW_FEASIBLE_LOAD_BOUND
```

and compute it from each realized cell. Do not write that G5 universally forbids every row above 9.6.

---

# 2. Mandatory replacement of the 52-row pre-scan

The current pre-scan is not balanced. Under a uniform scene and support-size normalization:

- pixels in the selected 4×4 checkerboard blocks receive `64+16+4=84` relative dose units;
- pixels in the omitted checkerboard blocks receive `16+4=20`;
- the ratio is `4.2:1`.

More seriously, `build_V0` rescales each pre-scan row by `1/(a^T xhat)` so that every row has exactly the arm load. The pre-scan cannot use `xhat` to set its own row intensity because `xhat` is produced by that pre-scan. This is a circular/oracle physical model.

Freeze this scene-independent 52-row replacement:

1. Partition the 8×8 grid of 4×4 blocks into 32 deterministic pairs. Use row-major block index `q=0,...,31` paired with `q+32`. Each fine row illuminates both blocks, support size 32, at amplitude `n/32=32`.
2. Use all sixteen 8×8 block rows, amplitude `n/64=16`.
3. Use the four 16×16 quadrant rows, amplitude `n/256=4`.

Then every pixel receives cumulative amplitude

\[
32+16+4=52,
\]

and

\[
\boxed{\frac1{52}\sum_{m=1}^{52}a_m=\mathbf1.}
\]

Therefore the mean pre-scan load is exactly the declared global load for every sum-normalized scene, without an oracle or row servo.

For local information after the pre-scan, use the actual frozen physical rows:

\[
\hat\rho_m=\bar\rho\,a_m^T\hat x,
\qquad
q_m=\frac{B^Ta_m}{a_m^T\hat x},
\]

and

\[
V_0=\sum_m\nu J_{\rm exact}(\hat\rho_m,\nu)q_mq_m^T.
\]

Do not replace all row kernels by `J(bar rho,nu)`.

Accounting convention:

- safe pre-scans use global mean load 0.05;
- budgeted fast pre-scans use global mean load 0.60;
- `RIDGE-FIXED` uses the common 0.60 fast pre-scan and retunes only its 972 main rows toward the ridge;
- all 52 counts enter final reconstruction for every arm.

The replacement patterns, exact balance identities, rank, and SHA256 must be frozen before any confirmatory scene is opened.

---

# 3. Dictionary and arm compliance

The v3 solver currently imports `v1._default_shapes`, which contains six 16-pixel families only. It does not implement the minimum R10 dictionary against which the certificate was promised.

Before freeze, the exact searchable dictionary must include all cyclic translates of:

1. scattered `k=16`;
2. scattered `k=32`;
3. solid 4×4;
4. `Lblob6×6` 16-pixel support;
5. bars/rectangles 1×16, 2×8, 4×4, 8×2, 16×1;
6. compact 32-pixel square/rectangle variants;
7. the frozen 16-pixel ring/annulus family;
8. `p=0` and guarded `p=1` weights for every applicable support.

The implementation must support variable `k`; concatenating every geometry into one fixed-width `(G,16)` array is not compliant.

Restore the R10 arm bookkeeping:

- `MATCH1` remains a mandatory no-gate arm;
- `FIXED*` must be selected once from `{SCAT32, LBLOB16, MATCH1}` using the frozen DEV rule and committed score table;
- if that table selected `LBLOB16`, state `FIXED*=LBLOB16` explicitly rather than merely calling it “strongest”;
- optional raster may remain omitted.

The spectral and A-risk guards use `FIXED*`, not an undocumented substitute. A separate dose-balancing fallback may use another predeclared exact design, but it must not be mislabeled as the G6 comparator.

---

# 4. Safe and fast policy definitions

Section 6’s “each arm’s own safe reference” is approved only with these exact conventions.

## Fixed arms

`SCAT16`, `SCAT32`, `LBLOB16`, `MATCH1`, and the selected `FIXED*` use the same frozen 972-row main multiset for their safe and fast curves. Only the global source/load changes.

## OED-DT

The endpoint evaluates the full adaptive policy:

- safe OED-DT: load palette `{0.025,0.05,0.10}` and hard mean incident-load budget `<=0.05`;
- fast OED-DT: the R11 palette and hard mean incident-load budget `<=0.60`.

Both use the same dictionary, pre-scan accounting, guards, and algorithm. Their optimized pattern multisets may differ. The paper must say that the Q90 result compares two resource-constrained OED policies, not a flux-only perturbation of one fixed matrix.

## OED-EQLOAD

Use one pointwise-equal-load design policy at 0.05 and one at 0.60; no gate. Because the kernel factors from geometry under exact equal load, any geometry differences caused solely by numerical noise must be reported.

## RIDGE-FIXED

Safe curve: the frozen 972-row LBLOB16 multiset at load 0.05. Fast curve: the same support multiset under the dwell-dependent global ridge policy. It remains no-gate.

No arm may use another arm’s safe curve, and no cross-arm Q90 ratio may enter `METHOD_SPEED_PASS`.

---

# 5. G5 rounding and exact mixture

The following parts are approved:

1. largest-remainder/apportionment as an initial integer design;
2. deterministic budget demotion;
3. keeping the best state seen along a nonmonotone prescribed-repair trajectory;
4. a strictly descending dose-penalty polish.

Do not call the original exchange trajectory “provably divergent.” The dev trajectory was nonmonotone and failed to converge; that is sufficient motivation for the polish.

The current `COMPLIANT-VIA-MIXTURE` status is not valid because only the continuous dose vectors are mixed. `Aout` and the returned atom counts remain the noncompliant design.

Freeze this rule:

1. A fallback design `FIXED_DOSE` must be precomputed as an **exact 972-row multiset** that itself passes incident budget, G5, physical peak, and row-count checks.
2. If polish fails, construct an exact integer mixture of the rounded OED multiset and `FIXED_DOSE`; the mixture coefficient is quantized through actual row counts, not applied only to information/dose vectors.
3. Recompute from the materialized rows: incident and detected budgets, per-pixel dose, `V_exact`, D-efficiency, spectral floor, A-risk, ceiling guards, and peak irradiance.
4. `COMPLIANT-VIA-MIXTURE` is permitted only if the returned matrix itself passes every guard.
5. Otherwise return `DOSE_GUARD_FAIL` or `ROUNDING_FAIL`; do not run confirmatory reconstruction with the noncompliant matrix.

There is a real lattice warning: with 972 uniform `k=16` rows, the mean support coverage is `972×16/1024=15.1875`; coverage 16 is `5.35%` above the mean. Thus the continuous uniform scatter16 measure currently used as fallback cannot automatically be realized within a ±5% exact band. `FIXED_DOSE` must solve the exact feasibility problem, for example by frozen variable loads/weights or another predeclared balanced row family. Do not relax the 5% band after seeing endpoints.

The reported `alpha<=0.2548` is therefore only a continuous feasibility diagnostic until the exact mixture exists.

---

# 6. Task subspace, objective ridge, and certificates

M1 currently has no declared estimable subspace, although R10’s theorem and guards require one. The full 1024-dimensional objective is stabilized with a design-dependent numerical ridge, and repeated exact rows can leave the unregularized matrix rank deficient. A certificate cannot be quoted against an undeclared, ridge-dependent target.

Freeze the least-disruptive convention already anticipated by the v2 G6 code:

\[
\boxed{r=200.}
\]

For each local cell, let `B` be the top 200 eigenvectors of the `FIXED*` information matrix under the same common pre-scan, dwell, and resource convention. Freeze `B` before the OED solve. Project `V0`, every atom, `FIXED*`, A-risk, and spectral comparisons into this same subspace. Reconstruction remains full 1024 pixels.

Use a fixed numerical ridge

\[
\epsilon_0
=10^{-9}\frac{\operatorname{tr}(B^TV_{\rm FIXED*}B)}{r}
\]

chosen before optimization and held constant. Do not recompute `epsilon` from the current design trace; a design-dependent ridge changes the objective derivative and invalidates the present sensitivity calculation.

## 6.1 Conservative relaxed certificate

Omitting per-pixel-dose multipliers can still yield a rigorous conservative certificate because removing those constraints enlarges the feasible set. Replace the current ratio gap by the additive budget-dual bound.

Let

\[
d_a=\operatorname{tr}(V^{-1}M_{\rm main}H_a),
\qquad c_a=\rho_a,
\]

let `theta>=0`, and let the budget be `c^T xi<=b`. Define

\[
\boxed{
G_{\rm rel}
=\max_a(d_a-\theta c_a)
-\sum_a\xi_a(d_a-\theta c_a)
+\theta\{b-c^T\xi\}.
}
\]

For a dose-feasible current design this upper-bounds the objective gap to the optimum of the larger simplex+budget problem and therefore also upper-bounds the gap to the full dose-constrained optimum.

Freeze:

\[
\boxed{G_{\rm rel}/r\le10^{-3}.}
\]

Call this field `RELAXED_KW_UPPER_BOUND`, not a full per-pixel KKT/equivalence certificate. The exact frozen dictionary must be scanned when computing the maximum. The present `max(s)/mean(s)-1`, without the slack term and with `tol=0.1`, is not accepted.

If per-pixel dual multipliers are later implemented, a full constrained equivalence certificate may be reported in addition.

## 6.2 Exact-design guards

After rounding, polish, and any materialized mixture, require all of the original R10 guards on the final returned rows:

\[
\operatorname{Eff}_D^{\rm exact}
=\exp\left\{\frac{\log\det V_{\rm exact}-\log\det V_{\rm cont}}{r}\right\}
\ge0.95,
\]

\[
\operatorname{tr}(V_{\rm exact}^{-1})
\le1.05\operatorname{tr}(V_{\rm FIXED*}^{-1}),
\]

and

\[
V_{\rm exact}\succeq0.5V_{\rm FIXED*}
\]

on `B`.

The current v3 output does not compute these three final exact-design gates. Add explicit fail states:

```text
DICTIONARY_RANK_FAIL
CONTINUOUS_CERTIFICATE_FAIL
SPECTRAL_GUARD_FAIL
A_RISK_GUARD_FAIL
DOSE_GUARD_FAIL
ROUNDING_FAIL
```

No `ALL HARD CHECKS PASS` message is valid until these are evaluated on the returned exact matrix.

---

# 7. Full-grid ridge and kernel certification

The values

```text
rho_R(200)  = 10.0369
rho_R(2000) = 22.2543
```

with no clip are approved as v3 development facts. They do not certify the production grid.

Before freeze, produce a source-of-record table for every

\[
\nu\in\{5,10,20,50,100,200,500,1000,2000\}
\]

containing:

- unconstrained exact ridge;
- production ridge;
- clip reason;
- every safe and fast palette level;
- exact PMF mass/score/derivative checks;
- `J_exact`, `J_mean`, efficiency, and RQL log-load bias;
- `p_ceil`;
- admission verdict.

Every node used by any arm must be G7-certified. The fact that `p_ceil` is astronomically small at `nu=200,2000` cannot be generalized to `nu=5,10`.

---

# 8. Items §12(a)–(f)

## (a) Family mapping

Approved. Reuse the committed glyph, chirp, spokes, maze, contour, and microtexture generators with fresh `633000+` seeds. In manifests and the paper use the actual generator names. “USAF” is not permitted unless a USAF object is actually generated; `spokes` may be described as a radial resolution target.

## (b) Each arm’s safe reference

Approved exactly as Section 4 of this ruling. The primary is a policy-level safe-versus-fast comparison within OED-DT. Add a sentence that design adaptation and flux allocation are both components of the OED policy; the separate high-flux OED-versus-`FIXED*` secondary isolates design value.

## (c) 52-row implementation

Not approved. Use Section 2’s balanced, scene-independent physical pre-scan and recompute all V0 terms.

## (d) Analyzer naming and cohort

Freeze the fresh 24 scenes as the only confirmatory decision units and use these source-of-record names:

```text
METHOD_SPEED_PASS          # primary Q90 gate, OED-DT only
METHOD_FIXED_DWELL_PASS    # safe-vs-fast terminal secondary
METHOD_DESIGN_PASS         # OED-DT vs FIXED* at fast terminal dwell
```

The latter restores the R10 key secondary: median design gain ≥0.50 dB, family-stratified lower bound >0, and at least 18/24 positive. Neither secondary rescues the primary. `RIDGE-FIXED` and other arms remain descriptive.

## (e) Ridge numbers and certificate convention

Ridge numbers approved for the two tested dwells. The current certificate convention is replaced by Section 6’s additive relaxed upper bound. Quote its exact formula and the fact that it is relative to the frozen dictionary and local pre-scan estimate.

## (f) Photon-economics framing

Approved frozen wording:

> Under a linear incident-photon budget, the relevant scalar efficiency is `J_exact(rho,nu)/rho`. In the v3 development cell this efficiency decreased across the admitted palette, and the active budget together with the per-pixel dose band led the optimizer to use only loads 0.3 and 0.6. This is a cell-specific constrained-design outcome, not a universal prohibition on ridge-zone allocation: sufficiently valuable directions may still receive larger loads in other cells. OED-DT tests the photon-budgeted allocation corner, whereas RIDGE-FIXED tests the time-limited, power-available ridge-tracking corner. The two arms expose conjugate resource regimes of the same dead-time information map.

The self-test requirement “ridge atoms must be selected” is correctly demoted to a load-histogram diagnostic. Conversely, absence of ridge atoms must not be advertised as a theorem solely from one self-test.

---

# 9. M1 specification edits required before freeze

The next draft must explicitly contain:

1. the balanced physical pre-scan and nonoracle V0 formula;
2. the complete frozen dictionary;
3. the DEV selection record for `FIXED*` and the mandatory `MATCH1` arm;
4. safe and fast OED load palettes/budgets;
5. the two-part peak guard;
6. the exact trust-efficiency and bias guards;
7. the declared `r=200` task subspace and fixed ridge;
8. the additive relaxed certificate and `10^-3` threshold;
9. exact materialization of any mixture;
10. D-efficiency, A-risk, spectral, dose, budget, ceiling, and physical-peak gates on the final returned rows;
11. all three analyzer verdict names;
12. the full nine-dwell kernel/ridge certification table.

---

# 10. Freeze authorization checklist

R13 changes from **HOLD** to **GO** only when one outcome-blind self-test/ledger demonstrates all of the following:

- [ ] pre-scan row average equals the all-ones vector to machine precision;
- [ ] pre-scan per-pixel cumulative dose is exactly equal;
- [ ] V0 is computed from actual row loads, with no `xhat`-dependent acquisition scaling;
- [ ] dictionary manifest contains every R10 family, every translation, both powers, support sizes, and SHA256;
- [ ] `FIXED*` selection scores and tie rule are committed;
- [ ] every arm returns exactly 52+972 physical rows;
- [ ] safe and fast OED policies both run through the final solver;
- [ ] all nine dwell values and all admitted loads pass kernel/trust/bias certification;
- [ ] the final continuous design passes `RELAXED_KW_UPPER_BOUND/r<=10^-3`;
- [ ] any fallback mixture is materialized as an exact 972-row matrix;
- [ ] final exact rows pass incident budget, detected-budget reporting, G5, physical peak, ceiling, D-efficiency, A-risk, and spectral guards;
- [ ] analyzer self-test emits `METHOD_SPEED_PASS`, `METHOD_FIXED_DWELL_PASS`, and `METHOD_DESIGN_PASS` on a DEV-only mock cohort;
- [ ] no hard check depends on endpoint PSNR or confirmatory scenes;
- [ ] expected-cell ledger and manifests are regenerated from the corrected spec.

Once those boxes pass and the resulting commit is immutable, **GO: tag `m1-freeze` and launch.** Until then, the current `ALL HARD CHECKS PASS` label is superseded by this audit and must not be used as freeze authorization.