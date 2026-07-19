# M1 method spec — R17 amendment (OPERATIVE; supersedes conflicting parts of ROUND63_METHOD_SPEC_M1.md)

Authority: GPT R17 ruling, GitHub issue #9, archived verbatim at
`docs/ROUND63_GPT_ROUND17_RULING_RAW.md` (posted 2026-07-19T10:36:40Z, audit
target commit 7ba4183). Where this amendment is silent, the R10–R16 spec
stands. Where they conflict, this amendment wins. No confirmatory scene,
pre-scan count, or endpoint has been opened as of this writing.

## A. Retired semantics (R17 §1, §3)

1. OED-DT is **retired** as a deployed/gated imaging arm. Its
   `METHOD_SPEED_PASS` gate is withdrawn *pre-freeze as infeasible* — it is
   never to be reported as a confirmatory failure.
2. Frozen rule: `m = 0 ⇒ ADAPTIVE_COLLAPSE_UNDER_GUARDS`; `m = 0 ≠ PASS`.
   The production semantics `COMPLIANT_VIA_MIXTURE`, `mixture PASS`, and
   `OED-DT with alpha=0` are deleted from all production code paths.
3. The descending exact scan survives only as the descriptive DEV diagnostic
   `PATH_FEASIBLE_ALPHA = max{m/972 : frozen ordered mixture at m passes}`,
   explicitly path-specific (frozen row ordering, not a subset max), reported
   in the methods-development supplement only. It cannot alter patterns,
   endpoints, certificates, or inclusion.
4. The dose-relaxed and dose-constrained OED solutions are retained as
   **no-endpoint design diagnostics** (they feed the Act III figure), never
   labeled the proposed imaging method.

## B. Arm roster (revised architecture)

Deployed spatial design: the balanced exact 972-row **SCAT32** multiset,
identical in all operating-point comparisons; the common 52-row balanced
pre-scan at global mean load 0.60 remains charged to every arm.

Endpoint-carrying operating modes (R17 §2.1; name `RIDGE-FIXED` retired):

- **SCAT32-SAFE** — global mean load 0.05.
- **SCAT32-060** — global mean load 0.60.
- **RIDGE-SCAT32** — one global source multiplier per dwell, calibrated from
  the common pre-scan estimate so the predicted mean main-pattern load hits
  the exact production ridge ρ_R(ν), followed only by the frozen global
  safety clip. No per-pattern servo.

Context arms (no gate, descriptive only, disclosed as such in the spec and
ledger): **SCAT16** and **LBLOB16** at load 0.60, retained so the cross-arm
occupancy-ladder tables in paper 2 keep their comparison columns. This
retention is a pre-freeze bookkeeping choice made before any confirmatory
observation; it adds no endpoint and will be listed explicitly in the
refreeze audit for external confirmation. The former OED-EQLOAD imaging arm
is dropped (its role is covered by the §A.4 design diagnostics).

## C. Endpoints (all frozen verbatim from R17 §2.2–2.3, §4.5)

### C.1 Primary — `RIDGE_OPERATING_PASS`

At terminal dwell ν=2000, per confirmatory image j:
ΔQ_ridge_j = mean-over-5-frozen-seeds PSNR_rad(RIDGE-SCAT32, ν=2000, j)
− same(SCAT32-060, ν=2000, j), with identical 52-row pre-scan at load 0.60,
identical 972 SCAT32 multiset, paired pattern/noise seeds, same
reconstruction path, same optical dwell; only the single global
main-acquisition source multiplier differs.

PASS iff ALL of: median_j ΔQ_ridge_j ≥ 1.0 dB; LB_2.5%[median ΔQ_ridge] > 0
(existing 10,000-replicate nested family-stratified bootstrap);
#{j: ΔQ_ridge_j > 0} ≥ 18/24. Intention-to-treat: any numerical/guard
failure counts as nonpositive.

Mandatory wording: “This is a fixed-dwell, power-for-time operating-point
test. It is not an equal-photon-budget or photon-efficiency comparison.”
Mandatory disclosures per cell: incident-dose ratio, detected-count ratio,
achieved mean load, load quantiles, physical peak, ceiling fraction,
`RIDGE_GUARD_CLIPPED` flag.

### C.2 Secondary — `RIDGE_SPEED_PASS`

Nine-dwell Q90 time-to-quality machinery (existing PAVA/censoring), comparing
SCAT32-SAFE (load 0.05) vs the dwell-dependent RIDGE-SCAT32 policy
(ρ_R(ν) per dwell). Bars: median S_gate ≥ 3; bootstrap lower bound > 1;
≥ 18/24 images with S_gate > 1. Cannot rescue the primary; primary cannot
rescue it.

### C.3 Confirmatory structural secondary — `DOSE_SAFE_CERT_PASS`

Certificate of near-optimality of the deployed SCAT32 design over the
expanded dose-safe class (R17 §4):

- **Dictionary**: D_cert = D_load ∪ D_gain, SHA256-frozen BEFORE any
  confirmatory scene is opened. D_load = existing frozen 13-family,
  all-translation, p∈{0,1} target-load-normalized dictionary. D_gain =
  certificate-only family: for every frozen base geometry and weight
  pattern, physical rows a_{g,γ} = γ·a_g^base, γ ∈ {0.2, 0.5, 1, 2, 5},
  base normalized by the scene-independent physical convention of the fixed
  patterns; predicted load EMERGES as ρ_{g,γ} = ρ̄·γ·(a_g^base)ᵀx̂ — never
  rescaled to a target load. Existing support/weight/physical-peak/ceiling/
  RQL-trust admission guards apply to both parts.
- **Target class** C_dose(x̂, ν, b): probability design measures over D_cert
  with ∫ρ_s dξ ≤ b and per-pixel dose within ±5% of mean. A-risk and
  spectral are NOT in the dual target class (larger class ⇒ stronger
  certificate); the deployed SCAT32 rows still pass all physical and
  conditioning disclosures separately. Same balanced pre-scan V0, declared
  r=200 task subspace, fixed numerical ridge as the campaign.
- **Bound**: R15 full dose-constrained KKT dual G_full (one budget
  multiplier, all ± pixel-dose multipliers, coverage-seeded column
  generation, final scan over ALL of D_cert). Retain: active-dose toy
  agreement ≤ 1e-8; primal feasibility; full-dictionary scan residual;
  dual/complementarity reporting; finite μ cap disclosed as
  `MU_CAP_ACTIVE`.
- **Scope**: every actual confirmatory pre-scan realization — 24 scenes ×
  5 seeds × ν∈{200, 2000} × b∈{0.05, 0.60} = **480 cells**.
- **Gate**: PASS iff EVERY cell is primal-feasible, numerically certified,
  and G_full/r ≤ 1e-2 (⇒ certified local D-efficiency ≥ exp(−0.01) =
  0.99005). If any cell fails: report the full distribution and fraction;
  the categorical near-optimality headline is forbidden.

The two resource corners are conjugate and mutually non-rescuing (R17 §4.1):
photon-budgeted corner → certificate; time-limited/power-available corner →
`RIDGE_OPERATING_PASS`. A ridge-budget certificate may appear later as
no-gate analysis only.

## D. Launch policy (R17 §5)

Split launch REJECTED. One revised spec, one expected-cell ledger, one
immutable tag, then one launch of all arms together.

Allowed pre-freeze: exact kernels and ridge tables independent of
confirmatory scenes; dictionary construction + SHA checks; DEV-only
certificate feasibility and runtime tests; scene-independent fixed pattern
matrices; mock analyzer tests. Forbidden: loading any 633000+ image into a
reconstruction or certificate; generating confirmatory pre-scan counts;
running any imaging arm on the confirmatory cohort.

## E. Refreeze authorization checklist (R17 §7 — selftest must map 1:1)

GO only after one outcome-blind ledger shows:
1. original OED-DT primary and `m=0 PASS` semantics removed;
2. exact SCAT32 deployed matrix and three operating modes frozen and hashed;
3. `RIDGE_OPERATING_PASS`, `RIDGE_SPEED_PASS`, `DOSE_SAFE_CERT_PASS`
   emitted correctly on positive AND negative mock cohorts/certificate toys;
4. D_cert = D_load ∪ D_gain complete and SHA-frozen;
5. expanded-class certificate passes the active-dose toy and
   full-dictionary scan tests;
6. DEV-only expanded-class feasibility confirms the certificate pipeline
   terminates without changing the 1e-2 bar;
7. common 52+972 accounting and all ridge disclosures pass;
8. all manifests and expected-cell ledgers regenerated from the revised
   architecture;
9. no confirmatory image, pre-scan count, or endpoint opened.

## F. Paper-2 Act III (R17 §6 — wording constraints binding)

Four-beat structure: (1) unconstrained information wants concentration
(dose-relaxed OED illustration); (2) uniform-dose safety removes that degree
of freedom (collapse within the frozen class); (3) simplicity is certified,
not assumed (expanded-class dual certificate); (4) global operating-point
control survives (ridge tracking, empirical primary). R11 + R17 unified as:
naive per-pattern equalization destroyed brightness–information alignment;
adaptive concentration conflicts with dose/conditioning safety; the balanced
global-flux policy is the robust middle.

Permitted and forbidden sentences: use the verbatim lists in
`ROUND63_GPT_ROUND17_RULING_RAW.md` §6.2. Final claim (frozen): “Strict
uniform-dose safety collapses the tested adaptive design space toward a
balanced pattern family, for which near-optimality can be certified;
high-flux benefit is then recovered by tuning the global detector operating
point.” Act III figure: (i) guard-vs-alpha DEV panel (labeled development
evidence); (ii) confirmatory G_full/r distributions; (iii) ridge-vs-0.60
fixed-dwell image + quality gains; (iv) photon-budgeted vs time-limited
resource-corner schematic.
