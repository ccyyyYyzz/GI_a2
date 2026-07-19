# R15 — Final ruling: certificate amendment and freeze authorization

**Audit target:** commit `55448e3`, especially `results/round63_m1/FREEZE_CHECKLIST_LEDGER.md` and `code/round63/oed_design_v4.py`.

**Scope:** this ruling amends only the R13 certificate/final-exact-design boxes and resolves the three implementation questions raised in R15. All other R10–R14 rulings remain in force.

> **CURRENT FREEZE VERDICT: HOLD, with a precise conditional GO path.**
>
> Adopt remedy **(a)**: a full constrained dual/KKT certificate that includes the per-pixel dose-band multipliers. Retain a separately solved relaxed problem as a solver-reference check. Reject the proposed transfer-bound remedy (b) as the freeze certificate. The next outcome-blind selftest is decisive if it implements the exact Box 11 and Box 13 wording below.
>
> One additional blocker was exposed by this audit: the committed `FIXED*` table uses an information score and a lexicographic tie rule, whereas R10 froze a DEV radiometric-PSNR selection rule with a different tie rule. Neither `MATCH1` nor `LBLOB16` may be pinned from the current table. The original R10 selection must be run, then every `FIXED*`-dependent object and guard must be regenerated.

## Executive ruling

| Item | R15 ruling |
|---|---|
| Remedy (a): full dose-constrained KKT certificate | **ADOPTED.** Include the budget multiplier and all upper/lower per-pixel dose multipliers. Hard gate: `FULL_CONSTRAINED_KW_UPPER_BOUND / r <= 1e-2`. A value `<=1e-3` is reported as a tight certificate but is not a second moving gate. |
| Remedy (b): relaxed certificate + transfer bound `T<=3.5` | **REJECTED as the freeze certificate.** The transfer factor is not a standard, objective-invariant optimality certificate for the full constrained problem and is too loose to replace the KKT bound. |
| Remedy (c): tolerance-only relaxation | **REJECTED.** The dose-feasible value `G_rel/r≈2.8` is a distance to the wrong feasible set, not a meaningful stationarity tolerance. |
| Support-restricted exact reweighting | **APPROVED and REQUIRED as a reference-solver upgrade.** It strengthens the relaxed lower bound, but a full-dictionary scan is still mandatory and it cannot replace the full constrained certificate. |
| `FIXED*` | **Current MATCH1 selection is not valid under R10.** Re-run the frozen PSNR rule. Pin whichever arm wins that rule; do not default to `MATCH1` or `LBLOB16` from the information-score table. |
| `m=0` FIXED_DOSE tolerance | **APPROVED.** Include the exact `m=0` endpoint and use `incident_sum <= cap + 1e-9`. |
| Box 13 mixture fallback | **APPROVED only as an exact, materialized, all-guards search.** Search the OED-row count down to zero and accept the largest OED fraction for which every final guard passes. |

---

# 1. Why the R13 relaxed certificate must be amended

The R13 bound

\[
G_{\rm rel}=\min_{\theta\ge0}
\left[\max_a(d_a-\theta c_a)-d^T\xi+\theta b\right]
\]

is a valid upper bound for the larger simplex-plus-incident-budget problem. It is therefore conservative for a dose-feasible design, but it cannot generally tend to zero when the omitted per-pixel dose constraints are binding.

The committed trajectories establish exactly this geometry:

- the same pairwise-FW machinery reaches `G_rel/r≈2.9e-2` and `2.5e-2` on the relaxed fast/safe problems and then becomes bit-stationary;
- dose projection monotonically lowers the physical dose violation while increasing the relaxed bound;
- the dose-feasible values plateau near `2.73` and `2.86`.

That is evidence that the relaxed optimum lies outside the dose polytope. It does **not** justify weakening `G_rel` until the dose-feasible design passes. It requires certifying the actual constrained problem.

---

# 2. Full constrained equivalence/dual certificate

## 2.1 Linear form of the dose band

For dictionary atom `a`, let

\[
z_a\in\mathbb R^n
\]

be its per-exposure physical dose vector in the design-measure convention, and let

\[
\bar z_a=\frac1n\mathbf1^Tz_a.
\]

With frozen dose half-width

\[
\delta=0.05,
\]

define the upper and lower dose coefficients

\[
u^+_{ja}=z_{aj}-(1+\delta)\bar z_a,
\qquad
u^-_{ja}=(1-\delta)\bar z_a-z_{aj}.
\]

The physical dose band is exactly the set of linear inequalities

\[
\sum_a\xi_a u^+_{ja}\le0,
\qquad
\sum_a\xi_a u^-_{ja}\le0,
\qquad j=1,\ldots,n.
\]

The incident budget remains

\[
\sum_a c_a\xi_a\le b,
\qquad c_a=\rho_a.
\]

## 2.2 Certificate formula

At a feasible continuous design `xi`, let

\[
d_a=\operatorname{tr}\!\left[V(\xi)^{-1}M_{\rm main}H_a\right].
\]

For multipliers

\[
\theta\ge0,
\qquad
\mu_j^+\ge0,
\qquad
\mu_j^-\ge0,
\]

define the reduced sensitivity

\[
s_a(\theta,\mu^+,\mu^-)
=d_a-\theta c_a
-\sum_j\mu_j^+u^+_{ja}
-\sum_j\mu_j^-u^-_{ja}.
\]

The full constrained upper bound is

\[
\boxed{
G_{\rm full}(\xi)
=
\min_{\theta,\mu^+,\mu^-\ge0}
\left\{
\max_{a\in\mathcal D}s_a
-d^T\xi+\theta b
\right\}.
}
\]

The dose constraints have zero right-hand side, so they do not add a separate constant term. Equivalently, this can be written as the maximum reduced-sensitivity excess plus all multiplier-weighted primal slacks.

For every feasible design `zeta`, concavity gives

\[
\log\det V(\zeta)-\log\det V(\xi)
\le G_{\rm full}(\xi).
\]

The multiplier minimization is a convex piecewise-linear problem. It may be solved as an epigraph LP or by certified column generation, but the final maximum must be evaluated over the **entire frozen dictionary**, not only the active support.

## 2.3 Frozen constants and interpretation

Freeze the hard gate

\[
\boxed{
G_{\rm full}/r\le10^{-2}.
}
\]

This is not a weak tolerance: it certifies D-efficiency relative to the full dose-constrained continuous optimum of at least

\[
\exp(-0.01)=0.99005.
\]

Also report

```text
FULL_CONSTRAINED_KW_TIGHT = (G_full/r <= 1e-3)
```

but do not alter the hard gate after the freeze according to whether the tighter target is reached.

Mandatory numerical checks:

- budget primal violation `<=1e-10` in design-measure units;
- dose-band primal excess `<=1e-8`;
- all dual multipliers nonnegative to `1e-10`;
- a full-dictionary reduced-sensitivity scan after the multiplier solve;
- any numerical negative bound may be set to zero only if it is no smaller than `-1e-10*r`; otherwise certificate failure;
- report budget and dose complementarity residuals;
- add one finite toy problem with an **active dose constraint**, solve it by brute force or an independent convex solver, and require objective/dual-bound agreement `<=1e-8`.

Field name:

```text
FULL_CONSTRAINED_KW_UPPER_BOUND
```

Do not call the budget-only `G_rel` a certificate for the dose-feasible optimum.

---

# 3. Relaxed reference solve and support-restricted exact reweighting

The support-restricted reweighting upgrade is approved and becomes a mandatory solver-reference check.

Freeze the support construction as:

1. all atoms with final relaxed weight `xi_a > 1e-12`;
2. union with the top **300** atoms from the final full-dictionary budget-adjusted sensitivity scan;
3. deterministic atom-index tie breaking.

On that fixed support, solve the simplex-plus-incident-budget weight problem to high precision. Require:

- support KKT residual divided by `r` `<=1e-8`;
- relative objective change over the last accepted update `<=1e-12`;
- then recompute the additive relaxed dual bound with a **full-dictionary** scan.

Freeze the relaxed-reference machinery check:

\[
\boxed{
\texttt{RELAXED_REFERENCE_KW_UPPER_BOUND}/r\le10^{-3}.
}
\]

This check certifies that the design/certificate implementation can solve the easier relaxed problem. It is not the certificate of the physical dose-constrained method.

The same support-reweighting machinery may be used as a warm start for the full constrained solve, but the full constrained column-generation loop must add outside-support atoms whenever the full reduced-sensitivity scan demands them.

---

# 4. Rejection of the proposed transfer-bound option

Do not freeze remedy (b) as stated.

A scalar condition such as

```text
T <= 3.5
```

between a relaxed solution and a projected or dose-feasible solution is not, without a new theorem, a dictionary-wide upper bound on the constrained log-determinant optimum. Its numerical meaning depends on the chosen objective normalization and projection path, and a factor of 3.5 would permit a materially large uncertainty in design quality.

It may be reported descriptively after the full constrained certificate exists, but it cannot authorize the campaign.

---

# 5. Exact amended Box 11 wording

Replace R13 checklist Box 11 with:

```text
[ ] relaxed_and_full_constrained_certificates

PASS iff, for BOTH the safe and fast OED policies:

(a) the dose-relaxed simplex+incident-budget reference problem, after the
    frozen support-restricted exact reweighting, satisfies
    RELAXED_REFERENCE_KW_UPPER_BOUND / r <= 1e-3 under a full-dictionary scan;

(b) the physically dose-feasible continuous design is primal-feasible for
    the incident budget and every +/-5% per-pixel dose inequality and satisfies
    FULL_CONSTRAINED_KW_UPPER_BOUND / r <= 1e-2, where the certificate includes
    one budget multiplier and all upper/lower per-pixel dose multipliers and
    the final maximum is scanned over the entire frozen dictionary;

(c) the full constrained certificate passes the independent active-dose toy
    verification and records all primal, dual, complementarity, and dictionary-
    scan residuals.
```

Fail state remains:

```text
CONTINUOUS_CERTIFICATE_FAIL
```

No transfer factor and no tolerance-only exception is permitted.

---

# 6. Final exact design and the materialized mixture fallback

## 6.1 Approved rounding path

The following sequence is approved:

1. deterministic apportionment/largest remainder;
2. deterministic incident-budget demotion;
3. best state retained from the prescribed dose-repair trajectory;
4. strictly descending dose-penalty polish;
5. if any final guard still fails, an exact materialized mixture search with the frozen `FIXED_DOSE` design.

The prescribed exchange trajectory is not required to be monotone. Keeping its best visited state, followed by a strict-descent polish, is acceptable because every subsequent certificate and guard is recomputed on the actual returned rows.

## 6.2 Exact mixture rule

Search integer OED row count

\[
m=972,971,\ldots,0
\]

in descending order. The `m=0` endpoint must be included explicitly.

For each `m`:

- construct an actual 972-row matrix containing exactly `m` frozen-order OED rows and `972-m` frozen-order `FIXED_DOSE` rows;
- recompute the information matrix and **every** final guard from those rows;
- accept the first, hence largest-OED-fraction, candidate satisfying all guards.

The frozen order/selection rule for both row sources must be committed and hashed. A continuous information-matrix or dose-vector mixture is not sufficient.

Use the approved floating-point comparison

```python
incident_sum <= cap + 1e-9
```

including at `m=0`.

The status

```text
COMPLIANT_VIA_MIXTURE
```

is valid only when the materialized exact matrix passes the complete guard suite. Report

```text
mixture_oed_rows = m
mixture_alpha = m / 972
```

A mixture may be heavily dominated by `FIXED_DOSE`; that is an honest result. No minimum alpha is added after seeing endpoints.

## 6.3 D-efficiency reference

For the exact-row gate, use the corresponding certified **full dose-constrained continuous design** as the `logdet V_cont` reference, not the dose-relaxed optimum.

---

# 7. Exact amended Box 13 wording

Replace R13 checklist Box 13 with:

```text
[ ] final_exact_rows_full_guard_suite

PASS iff, for BOTH safe and fast policies, the final returned main-pattern
matrix contains exactly 972 physical rows and is either the direct rounded OED
design or the materialized COMPLIANT_VIA_MIXTURE design selected by the frozen
m=972,...,0 search. On that exact matrix, all of the following must pass:

- incident_sum_rho <= 972*budget_mean + 1e-9;
- predicted and realized detected-budget fields are emitted;
- per-pixel cumulative-dose max relative deviation <= 0.05;
- physical peak <= 1536;
- design-weighted ceiling probability <= 0.01;
- exact D-efficiency >= 0.95 relative to the certified full dose-constrained
  continuous design;
- A-risk <= 1.05 times the load-matched FIXED* reference;
- minimum generalized eigenvalue versus the load-matched FIXED* reference
  >= 0.5;
- exactly 52 common pre-scan rows plus 972 returned main rows are accounted.

If no integer mixture candidate, including m=0, passes all guards, Box 13 FAILS.
```

The corresponding fail state is whichever first failed guard applies; do not label a design `OK` merely because dose, budget, and peak pass.

---

# 8. `FIXED*` selection ruling

The current committed selection is not compliant with R10.

R10 froze:

1. six development images only;
2. each candidate reconstructed through the frozen production path at `rho=0.60`, `nu=2000`;
3. average radiometric PSNR over the frozen development measurement seeds for each image;
4. score = median of those six image means;
5. ties within `0.05 dB` resolved in the order

```text
SCAT32 -> LBLOB16 -> MATCH1
```

The committed table instead uses a local information/log-eigenvalue score and a lexicographic tie rule. Those values may remain as a design diagnostic, but they cannot select the confirmatory comparator.

Therefore:

- do **not** pin `MATCH1` from `3736.8 > 3736.3`;
- do **not** default to `LBLOB16` from earlier wording;
- run the original R10 PSNR rule and pin its winner;
- if the winner is `MATCH1`, then `MATCH1` is accepted without further objection;
- after the winner is fixed, regenerate `B`, `eps0`, `V_FIX`, all A-risk/spectral references, Boxes 11 and 13, manifests, and the ledger.

This means the present ledger's Box 5 must be reopened. The current `14/16` count is not the final count until this correction is made.

---

# 9. Freeze authorization

The next immutable candidate commit receives **GO** only if one outcome-blind ledger shows:

1. the corrected R10 PSNR-based `FIXED*` selection and all dependent matrices regenerated;
2. amended Box 11 PASS for safe and fast;
3. amended Box 13 PASS for safe and fast;
4. every other R13 box still PASS after regeneration;
5. no confirmatory scene or endpoint is read by any hard check.

If those conditions hold:

> **GO: create the immutable `m1-freeze` commit/tag and launch.**

Until that rerun, commit `55448e3` remains **HOLD**.