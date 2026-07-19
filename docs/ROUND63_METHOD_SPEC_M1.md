# M1 — Method-campaign specification (assembles R10 + R11; R13 audit target)

Status: DRAFT for GPT round-13 final audit. Freezes upon the R13 ruling +
one immutable commit (tag m1-freeze). After that commit: no scene, atom
family, load grid, guard, endpoint, or gate may change; if OED-DT fails,
there is no third redesign (R10 Q7 §1).

## 1. Identity

Separate follow-up campaign for the method line (paper 2). Own seed
namespace, manifests, expected-cell ledger. Nothing here reopens paper-1
(frozen at tag paper1-freeze-v1) or the Study-1/2 verdicts.

## 2. Geometry and constants (R10 Q7 §2, verbatim)

n = 32² = 1024; M_total = 1024 = M_pre (52) + M_main (972); nonparalyzable
active-start, τ = 50 ns, d = 0, σ_b = 0; ρ̄ ∈ {0.05, 0.60} for the
budget-anchored arms; ν ∈ {5,10,20,50,100,200,500,1000,2000}; 5
measurement seeds per (image, condition). T_opt = M_total·ν·τ.
Design-computation wall time reported separately, never inside T_opt.

## 3. Scenes (R10 Q7 §3)

24 fresh confirmatory instances: six families × four instances, seeds
s_{f,r} = 633000 + 100f + r; six dev instances 632900+f. **Family mapping
decision for R13 blessing:** the six committed detail32 generators
(glyph, chirp, spokes, maze, contour, microtexture) are reused verbatim
as the "same broad families"; R10's family list names "USAF/line-pair and
chirp" as family 1 — we map it to the chirp generator and note that
spokes serves as the resolution-target analog (per the frozen ladder
addendum precedent). No new generator code; no acceptance filter of any
kind. Generator params + images + ROIs receive SHA256 manifests.

## 4. Pre-scan (R10 Q7 §4, verbatim)

Every arm receives and pays for the same 52-pattern frozen balanced
multiscale pre-scan at its own (ρ̄, ν): the deterministic 52-row set
implemented in oed_design_v2/v3 (32 even-parity 4×4 blocks + 16 8×8
blocks + 4 quadrants). All arms may fold pre-scan counts into their final
reconstruction; only adaptive arms use it to choose patterns. Fixed arms
run 972 main patterns; every arm totals exactly 1024 physical exposures.

## 5. Load architecture (R11, verbatim)

- Exact ridge target ρ_R(ν): min-argmax of J_exact over
  R_adm(ν) = {ρ ≤ 24, p_ceil(ρ,ν) ≤ 0.01, J_exact/J_q ≥ 0.90}; recorded
  fields per cell: rho_ridge_exact_unconstrained, rho_R_production,
  ridge_clip_reason, J_exact_at_target, J_q_at_target, p_ceiling_scalar.
- OED-DT palette L_ν = unique{0.30, 0.60, 1.00, min(3, ρ_R), ρ_R/2, ρ_R};
  atom admission adds p_ceil ≤ 0.05 and J_exact/J_q ≥ 0.90 to G3/G4.
- Budget: hard total incident-signal inequality Σ ρ_i ≤ M_main·0.60
  (the original mean-0.60 resource); incident AND detected budgets
  reported jointly.
- RIDGE-FIXED: LBLOB16 supports, ONE global multiplier per dwell
  calibrated on the pre-scan estimate to mean load ρ_R(ν), then the R11
  bisection safety clip; requested vs achieved load + RIDGE_GUARD_CLIPPED
  recorded. No per-pattern servo anywhere (R11-d prohibition).

## 6. Arms

| arm | role | gate |
|---|---|---|
| OED-DT | variable-load exact-kernel D-optimal (v3 solver, KW certificate vs the frozen dictionary) | **primary gate** |
| OED-EQLOAD | same optimizer, all atoms at ρ̄ | no gate (kernel ablation) |
| SCAT16 | frozen scattered k=16 | no gate (Study-2 bridge) |
| SCAT32 | best ladder rung by fixed-dwell median | no gate |
| LBLOB16 | committed clustered support | no gate (strongest fixed baseline; G6 reference) |
| RIDGE-FIXED | ridge-tracking policy | no gate (demonstration) |

Reconstruction: RQL + frozen analytic λ_TV rule + pooled init, identical
across arms (Study-2 machinery verbatim). Primary comparison: each arm's
own fast curve vs the common safe (ρ̄=0.05) reference of its own pattern
set, per the transferred Q90 endpoint.

## 7. Endpoints and gate (R10 Q7 + R11-c, transferred verbatim)

Q90 time-to-quality with PAVA, censoring taxonomy, 10,000-replicate
nested family-stratified bootstrap; gate on OED-DT ONLY: median S ≥ 3,
LB2.5 > 1, ≥ 18/24. Fixed-dwell secondary (median ≥ 1 dB at terminal
dwell, LB > 0, ≥ 18/24) cannot rescue. Mandatory per-cell disclosure:
incident/detected photons (predicted+realized), dose ratios, source/peak
irradiance ratio, ρ_5/ρ_50/ρ_95/ρ_max, ceiling fractions
(predicted+realized), mean J_exact. Fixed-dwell results additionally
carry the R11-c incident-dose/ceiling disclosure block.

## 8. Guards (R10 Q4 G1–G6 + R11 additions, all verbatim)

G1 chain-rule atoms on exact J table; G2 superseded by R11 load palette;
G3 peak ≤ 64; G4 weights ∈ [1/4,4], p ∈ {0,1}; G5 per-pixel dose ±5%
enforced in-optimizer + verified post-rounding; G6 spectral floor
V_OED ⪰ 0.5·V_LBLOB16 with α-mixture fallback; R11: ρ ≤ 24 cap,
p_ceil caps, RQL-trust ratio, kernel numerical certification, cell-level
post-rounding ceiling guard.

## 9. Analysis plan

Frozen analyzer = pilot_s1 endpoint machinery (cohort = the 24 fresh
scenes), executed once per arm after its grid completes; OED-DT verdict
named METHOD_SPEED_PASS. Descriptive cross-arm table (no new gates).

## 10. Compute plan

Local 2–3 CPU processes (no GPU anywhere; pure numpy). Colab pro1
available now, pro2 after ~2 h (user grant); shard by (arm × scene
block) with the proven v22 CLI infra if wall-clock requires. Watchdogs +
30-min digests as standard.

## 11. Learned proposer (no-claim ablation; R10 Q7 §5 compliant)

A small CPU-trained CNN proposer for the inner FW oracle, trained by
imitation of the v3 solver's own exhaustive atom sweeps on DEV scenes
only. Reported solely as a certificate-gap ablation (dictionary-only vs
+proposer); carries no gate, no scientific claim; dropped without
prejudice if it underperforms. Runs only during campaign idle time.

## 12. Items for R13 to bless or amend

(a) §3 family mapping (chirp covers "USAF/line-pair and chirp");
(b) §6 "each arm's own safe reference" reading of the endpoint transfer;
(c) the 52-row pre-scan implementation as committed in v2/v3 (even-parity
choice documented there);
(d) analyzer naming/cohort conventions;
(e) v3 selftest results: ρ_R(200) = 10.0369 (clip NONE), ρ_R(2000) =
22.2543 (clip NONE) — both G7-certified; RQL-trust ratio holds even at
the ν=2000 ridge. Constrained-gap certificate convention as implemented
in v3 (to be quoted verbatim in the freeze).
(f) **Photon-economics discovery (v3 selftest, needs R13 blessing of the
framing):** under the R11 total-incident-budget inequality the per-photon
information efficiency J(ρ)/ρ falls monotonically (0.55 at ρ=0.6 vs
0.036 at ρ=22), so the OED-DT optimizer RATIONALLY avoids ridge-zone
atoms — ridge operation is a TIME-constrained regime, not a
photon-constrained one. Proposed framing: OED-DT (photon-budgeted) and
RIDGE-FIXED (time-budgeted) are the two conjugate corners of the same
resource trade; the campaign measures both. The "ridge atoms selected"
selftest check is demoted to a diagnostic accordingly.
