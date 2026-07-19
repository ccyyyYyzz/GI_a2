# R17 ruling (GitHub issue #9, raw)

Title: R17 ruling: safety-constrained OED collapse and coherent refreeze
Posted: 2026-07-19T10:36:40Z

---

# R17 — Final ruling: safety-constrained OED collapse and coherent refreeze

**Audit target:** commit `7ba4183`, especially `docs/ROUND63_GPT_ROUND17_QUESTION.md` and `results/round63_m1/R17_CONSULTATION_EVIDENCE.md`.

**Scope:** this ruling governs the still-unfrozen M1 follow-up campaign only. No confirmatory scene has been opened. It does not reopen the paper-1/OE scope lock or alter R14–R16 theory results.

> **FREEZE VERDICT: HOLD. Do not launch any confirmatory arm yet.**
>
> The `m=0` outcome has two distinct meanings. Allowing `m=0` to count as a passing OED-DT arm is a specification defect because the deployed arm contains no OED rows. The fact that every nonzero mixture fails is, however, a genuine structural negative for the **current target-load-normalized OED architecture under the declared guard stack**. The guards must not be widened or re-anchored post hoc. The campaign may be re-architected around a balanced fixed design and a near-optimality certificate, but the certificate target class must first be expanded to include bounded **gain-coupled/fixed-flux atoms**; otherwise the claimed collapse would be an artifact of the anti-brightness load-normalization parameterization.

## Executive ruling

| Question | R17 ruling |
|---|---|
| **Q1 — defect or structural negative?** | **Both, at different levels.** `m=0 => OED-DT PASS` is a procedural defect and is retired. The no-nonzero-mixture result is accepted as a structural negative for the frozen load-targeted atom class. It is not yet a universal negative for dose-safe adaptivity; a broader gain-coupled certificate class is mandatory before that wording is used. Options (a) and (b) are rejected. |
| **Q2 — inverted claim and endpoints** | Adopt option (c), amended. The new empirical primary is a paired **SCAT32 ridge-operating-point** test, not an adaptive-OED gate. A separate confirmatory certificate tests whether balanced SCAT32 is within 1% D-efficiency of the best dose-safe design in the expanded frozen class at the original photon-budgeted corners. The certificate cannot rescue the empirical primary, and the empirical primary cannot rescue the certificate. |
| **Q3 — mixture rule** | Remove the mixture from every deployed/gated arm. `m=0` is `ADAPTIVE_COLLAPSE_UNDER_GUARDS`, never PASS. The descending scan may remain only as a frozen, path-specific descriptive diagnostic `PATH_FEASIBLE_ALPHA`; it is not an estimator, fallback, or certificate. |
| **Q4 — certificate target/tolerance** | Certify SCAT32 against the continuous dose-safe optimum over the union of the existing target-load atoms and a new bounded gain-coupled atom family. Use the full dose KKT certificate, full-dictionary scan, and `G_full/r <= 1e-2`; this implies at least `exp(-0.01)=0.99005` local D-efficiency. The certificate is a confirmatory structural secondary. |
| **Q5 — launch policy** | **One coherent freeze is mandatory.** None of the five non-gated arms may touch the confirmatory scenes before the revised spec, endpoints, target-class SHA, manifests, and decisive ledger are frozen together. Scene-independent kernels and DEV-only computation may continue. |
| **Q6 — paper framing** | Scientifically sound and OE-appropriate after narrowing: under the declared uniform-dose and detector constraints, spatial/per-pattern load adaptation collapses within the frozen class, while global ridge tracking remains useful. Do not write that adaptivity is universally impossible or that all safety constraints destroy all adaptive sensing. |

---

# 1. Q1 — What the `m=0` result establishes

## 1.1 The arm semantics are defective

The current accepted matrix at `m=0` consists of 972 fallback rows and zero OED rows. It therefore cannot be named, gated, or interpreted as OED-DT. Freeze the rule:

```text
m = 0  => ADAPTIVE_COLLAPSE_UNDER_GUARDS
m = 0  != PASS
```

The original `METHOD_SPEED_PASS` gate for OED-DT is withdrawn **before confirmatory freeze as infeasible**, not reported later as a confirmatory failure.

## 1.2 The observed incompatibility is real within the implemented architecture

The evidence is stronger than a rounding accident:

- the pure OED exact design has dose deviation `0.1891`, A-risk `1.95x`, and generalized spectral floor `0.126`;
- every exact mixture with `1 <= m <= 972` violates at least one of dose, A-risk, or spectral safety;
- the dose-feasible design is far from the dose-relaxed optimum;
- balanced SCAT32/FIXED_DOSE outperforms the current dose-feasible continuous reference.

Do not widen the ±5% dose band enough to admit a 19% deviation, and do not move the A-risk/spectral anchors after seeing this result. Options (a) and (b) would erase the preregistered safety question.

## 1.3 The negative is not yet universal

The present atom construction chooses a per-pattern amplitude to hit a target predicted load. On dim supports this introduces large physical amplitude, so the dose vector is anti-aligned with brightness. R11 and R14-C already showed that a fixed global source naturally preserves brightness–information alignment.

Therefore the current evidence supports only:

> Within the target-load-normalized variable-load atomization and frozen guard stack, no nonzero exact adaptive mixture is admissible in the audited development cell.

Before the manuscript says that dose-safe adaptivity collapses, the certificate target class must include gain-coupled/fixed-flux atoms as specified below. This is the final escape-hatch test; no further scene or endpoint redesign follows it.

---

# 2. Q2 — Re-architected empirical primary

## 2.1 Deployed design

Freeze the balanced exact 972-row SCAT32/FIXED_DOSE multiset as the deployed spatial design. Use the same pattern multiset in all operating-point comparisons. The common 52-row balanced pre-scan remains charged to every arm.

The previous LBLOB-based name `RIDGE-FIXED` is retired for the primary. Use:

```text
SCAT32-SAFE
SCAT32-060
RIDGE-SCAT32
```

`RIDGE-SCAT32` uses one global source multiplier per dwell, calibrated from the common pre-scan estimate so that the predicted mean main-pattern load targets the exact production ridge `rho_R(nu)`, followed only by the frozen global safety clip. There is no per-pattern servo.

## 2.2 Primary endpoint

At terminal dwell `nu=2000`, define for each confirmatory image `j`

\[
\Delta Q^{\rm ridge}_j =
\overline{\operatorname{PSNR}_{\rm rad}}
(\mathrm{RIDGE\!\! -\!SCAT32},\nu=2000,j)
-
\overline{\operatorname{PSNR}_{\rm rad}}
(\mathrm{SCAT32\! -\!060},\nu=2000,j),
\]

where the bar is the mean over the five frozen measurement seeds.

Both sides use:

- the identical 52-row pre-scan at global mean load `0.60`;
- the identical 972 SCAT32 main-pattern multiset;
- paired pattern/noise seeds and the same reconstruction path;
- the same optical dwell.

Only the single global main-acquisition source multiplier differs.

Freeze the primary verdict name and gate:

```text
RIDGE_OPERATING_PASS
```

PASS iff all three hold:

\[
\operatorname{median}_{j=1}^{24}\Delta Q^{\rm ridge}_j\ge1.0\ {\rm dB},
\]

\[
\mathrm{LB}_{2.5\%}
[\operatorname{median}\Delta Q^{\rm ridge}]>0,
\]

and

\[
\#\{j:\Delta Q^{\rm ridge}_j>0\}\ge18/24.
\]

Use the existing 10,000-replicate nested family-stratified bootstrap. Any numerical/guard failure remains in the intention-to-treat cohort and counts as nonpositive.

Mandatory wording:

> This is a fixed-dwell, power-for-time operating-point test. It is not an equal-photon-budget or photon-efficiency comparison.

Report the incident-dose ratio, detected-count ratio, achieved mean load, load quantiles, physical peak, ceiling fraction, and any `RIDGE_GUARD_CLIPPED` flag.

## 2.3 Q90 time-to-quality secondary

Retain the nine-dwell Q90 machinery as a key secondary:

```text
RIDGE_SPEED_PASS
```

It compares `SCAT32-SAFE` at global load `0.05` with the dwell-dependent `RIDGE-SCAT32` policy using the existing PAVA/censoring machinery and the existing bars:

- median `S_gate >= 3`;
- bootstrap lower bound `>1`;
- at least `18/24` images with `S_gate>1`.

It cannot rescue `RIDGE_OPERATING_PASS` and vice versa.

## 2.4 Original adaptive arm

OED-DT is no longer a deployed or gated reconstruction arm. A dose-relaxed OED solution and a dose-constrained OED solution may be retained as no-endpoint design diagnostics, but neither may be labeled the proposed imaging method.

---

# 3. Q3 — Mixture audit

## 3.1 Remove the fallback from deployment

The production design is exactly SCAT32. Do not generate it by running an OED arm and falling back. `FIXED_DOSE` is a baseline/design in its own right, not an OED repair.

Delete these production semantics:

```text
COMPLIANT_VIA_MIXTURE
mixture PASS
OED-DT with alpha=0
```

## 3.2 Residual diagnostic

The descending exact scan may remain as:

```text
PATH_FEASIBLE_ALPHA = max{m/972 : the frozen ordered mixture at m passes}
```

with two restrictions:

1. it is explicitly relative to the frozen row ordering and is not a global maximum over all subsets;
2. it is descriptive only and cannot alter patterns, endpoints, certificates, or inclusion.

The existing DEV value may be reported in a methods-development supplement. It is not itself confirmatory evidence. No confirmatory mixture scan is required once the stronger near-optimality certificate below is run.

---

# 4. Q4 — Confirmatory near-optimality certificate

## 4.1 Why the certificate and ridge arm use different resource corners

Do not claim that the current certificate proves ridge-load SCAT32 optimality. The structural collapse was observed at the original photon-budgeted policies `b=0.05` and `b=0.60`. The architecture contains two conjugate results:

1. **photon-budgeted corner:** balanced SCAT32 is tested for certified near-optimality among dose-safe designs;
2. **time-limited/power-available corner:** global ridge tracking is tested empirically by `RIDGE_OPERATING_PASS`.

A ridge-budget certificate may be reported later as no-gate analysis, but it is not required and cannot replace the two frozen claims.

## 4.2 Expanded certificate dictionary

Let the existing frozen 13-family, all-translation, `p in {0,1}` dictionary be `D_load`. It contains the current target-load-normalized safe/fast atoms.

Add a certificate-only family `D_gain`. For every frozen base geometry and weight pattern, form physical rows

\[
a_{g,\gamma}=\gamma\,a_g^{\rm base},
\qquad
\gamma\in\{0.2,0.5,1,2,5\},
\]

where `a_g^base` is normalized by the scene-independent physical convention used by the fixed patterns; its predicted load is allowed to emerge as

\[
\rho_{g,\gamma}=\bar\rho\,\gamma\,a_g^{\rm base\,T}\hat x.
\]

It is **not** rescaled to hit a target load. This finite family includes fixed-flux self-water-filling and bounded per-pattern gain control. The final target dictionary is

\[
\boxed{D_{\rm cert}=D_{\rm load}\cup D_{\rm gain}.}
\]

Apply the existing support, weight, physical-peak, ceiling, and RQL-trust admission guards to both parts. Freeze and record a SHA256 for `D_cert` before any confirmatory scene is opened.

If SCAT32 is not near-optimal against this expanded class, the manuscript may not claim that adaptivity collapses under safety constraints.

## 4.3 Target feasible class

For local pre-scan estimate `xhat`, dwell `nu`, and incident budget `b`, define `C_dose(xhat,nu,b)` as all probability design measures over `D_cert` satisfying:

\[
\int \rho_s(a)\,d\xi(a)\le b,
\]

and every per-pixel dose inequality

\[
-0.05\le d_j(\xi)/\bar d(\xi)-1\le0.05.
\]

Pointwise peak, ceiling, and trust restrictions are enforced through atom admission. Use the same balanced pre-scan information `V0`, declared `r=200` task subspace, and fixed numerical ridge as the campaign.

A-risk and spectral constraints are **not** included in the dual target class. This makes the certified class larger; a near-optimality result over this larger budget+dose class is stronger. The deployed SCAT32 rows must still pass all physical and conditioning disclosures separately.

## 4.4 Full constrained certificate

At the exact SCAT32 design measure, compute the R15 full dose-constrained bound

\[
G_{\rm full}=\min_{\theta,\mu^+,\mu^-\ge0}
\left[
\max_{a\in D_{\rm cert}}s_a(\theta,\mu^+,\mu^-)
-d^T\xi+\theta b
\right],
\]

with one incident-budget multiplier, every upper/lower pixel-dose multiplier, coverage-seeded column generation, and a final scan over **all of `D_cert`**.

Freeze:

\[
\boxed{G_{\rm full}/r\le10^{-2}.}
\]

This certifies local D-efficiency of at least

\[
\exp(-0.01)=0.99005
\]

relative to the continuous optimum of the expanded dose-safe class.

Retain:

- active-dose toy agreement `<=1e-8`;
- primal dose/budget feasibility;
- full-dictionary scan residual;
- dual/complementarity reporting;
- the finite `mu` cap as a conservative restriction, with `MU_CAP_ACTIVE` disclosed.

## 4.5 Confirmatory scope and verdict

Run the certificate on every actual pre-scan realization for:

- all 24 confirmatory scenes;
- all five measurement seeds;
- `nu in {200,2000}`;
- `b in {0.05,0.60}`.

This is 480 local certificate cells.

Freeze:

```text
DOSE_SAFE_CERT_PASS
```

PASS iff **every** one of the 480 cells is primal-feasible, numerically certified, and satisfies `G_full/r <= 1e-2`.

This is a confirmatory structural secondary. It authorizes the sentence “SCAT32 is certified within 1% D-efficiency of the best design in the frozen dose-safe class at the two preregistered dwell anchors.” It cannot rescue a failed `RIDGE_OPERATING_PASS`.

If some cells fail, report the complete distribution and fraction, but do not use the categorical near-optimality headline.

---

# 5. Q5 — Launch policy

Reject split launch.

The primary arm, comparator semantics, certificate target class, analyzer names, expected cells, and paper claim have changed. Running any confirmatory fixed arm now would expose outcomes before the final architecture is frozen and would create pressure to retain or drop arms based on observed images.

Allowed before the coherent freeze:

- exact kernels and ridge tables independent of confirmatory scenes;
- dictionary construction and SHA checks;
- DEV-only certificate feasibility and runtime tests;
- scene-independent fixed pattern matrices;
- mock analyzer tests.

Not allowed:

- loading any `633000+` confirmatory image into a reconstruction or certificate;
- generating confirmatory pre-scan counts;
- running the five non-gated imaging arms on the confirmatory cohort.

One revised spec, one expected-cell ledger, one immutable tag, then one launch.

---

# 6. Q6 — Paper-2 Act III framing

## 6.1 Scientifically sound framing

The proposed inversion is sound after the certificate-class amendment. The strongest architecture is:

1. **Unconstrained information wants concentration.** The dose-relaxed OED illustrates the information incentive to concentrate spatial and pattern load.
2. **Uniform-dose safety removes that degree of freedom.** Under the declared budget and ±5% dose constraints, the current adaptive architecture collapses to the balanced design.
3. **Simplicity is certified, not assumed.** The expanded-class dual certificate quantifies how little local D-information remains available beyond SCAT32.
4. **Global operating-point control survives.** Ridge tracking changes source power globally without spatially redistributing dose and is tested by the empirical primary.

R11 and R17 should not be summarized as “all adaptivity is bad.” R11 refuted naive per-pattern load equalization because it destroyed brightness–information alignment; R17 shows that the opposite extreme—adaptive concentration—conflicts with strict dose/conditioning safety. The balanced global-flux policy occupies the robust middle.

## 6.2 Frozen manuscript wording

Permitted:

> Within the frozen nonnegative illumination dictionary, local pre-scan linearization, incident budget, and ±5% per-pixel dose constraint, the balanced SCAT32 design was tested against a full constrained dual upper bound. Global ridge tracking was evaluated separately as the time-limited, power-available control.

> Under the declared dose-safety constraints, per-pattern load adaptation provided no certified material D-information advantage over the balanced design at the preregistered anchor dwells.

> The useful surviving control was the global operating point rather than spatial redistribution of incident dose.

Not permitted:

- “adaptivity is impossible under dead time”;
- “safety constraints always eliminate adaptive sensing”;
- “SCAT32 is globally optimal over all physical patterns”;
- “the m=0 DEV result is a confirmatory theorem”;
- “the ridge arm is photon efficient.”

## 6.3 OE suitability

This is OE-appropriate if both pieces are shown in the main text:

- an image-level ridge operating-point result with honest power/dose disclosure;
- a certificate/constraint figure showing the dose-relaxed optimum, the safety frontier, and the SCAT32 upper-bound gap.

A clean Act III figure should include:

1. dose deviation, A-risk, and spectral guard versus path-specific `alpha` on a DEV representative, explicitly labeled development evidence;
2. confirmatory `G_full/r` distributions for the expanded class;
3. ridge-versus-0.60 fixed-dwell image and quality gains;
4. a schematic distinguishing photon-budgeted and time-limited resource corners.

The final claim is not “adaptive OED won.” It is:

> **Strict uniform-dose safety collapses the tested adaptive design space toward a balanced pattern family, for which near-optimality can be certified; high-flux benefit is then recovered by tuning the global detector operating point.**

---

# 7. Revised freeze authorization

The next candidate commit receives GO only after one outcome-blind ledger shows:

1. the original OED-DT primary and `m=0 PASS` semantics are removed;
2. the exact SCAT32 deployed matrix and three operating modes are frozen and hashed;
3. `RIDGE_OPERATING_PASS`, `RIDGE_SPEED_PASS`, and `DOSE_SAFE_CERT_PASS` are emitted correctly on positive and negative mock cohorts/certificate toys;
4. `D_cert = D_load union D_gain` is complete and SHA-frozen;
5. the expanded-class certificate passes the active-dose toy and full-dictionary scan tests;
6. DEV-only expanded-class feasibility confirms the certificate pipeline can terminate without changing `1e-2`;
7. the common 52+972 accounting and all ridge disclosures pass;
8. all manifests and expected-cell ledgers are regenerated from the revised architecture;
9. no confirmatory image, pre-scan count, or endpoint has been opened.

Only then: create one immutable revised M1 freeze tag and launch all arms together.