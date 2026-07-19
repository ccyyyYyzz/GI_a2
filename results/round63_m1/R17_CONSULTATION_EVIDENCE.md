# R17 consultation evidence — OED-DT mixture degenerates to alpha = 0

Status: freeze HELD by the vacuous-arm stop rule (coordinator directive,
2026-07-19). Assembled from the freeze-9 selftest (ledger 2026-07-19 10:50,
13/14 PASS) before scratchpad log rotation; quoted lines were read from
scratchpad/freeze9.log in-session.

## 1. Box-11 "PASS m=0" semantics (the gating question)

In oed_design_v5.mixture_search, `m` counts OED rows retained and the
accepted matrix is `vstack(oed_rows[:m], fixed_dose_rows[...])` with
`mixture_alpha = m/972`. Therefore:

> **PASS m=0, alpha=0.0 = the accepted 972-row design contains ZERO OED
> rows and 972 FIXED_DOSE rows (100% fallback) — for BOTH policies.**

Quoted freeze-9 evidence:

```
[fast] ... mixture=PASS m=0 alpha=0.0
[safe] ... mixture=PASS m=0 alpha=0.0  (935s)
    guards: {'budget': True, 'dose': True, 'ceiling': True, 'effD': True,
             'arisk': True, 'spectral': True}
    incident=48.60 effD=2.2498 dose=0.0400 arisk=1.046x mu_min=0.651
[PASS] final_exact_rows_full_guard_suite -- fast: PASS m=0; safe: PASS m=0
```

At the pure-OED endpoint (fast policy, m=972, rounded + repaired design):

```
[mix] m=972 checks={'budget': True, 'dose': False, 'ceiling': True,
                    'effD': True, 'arisk': False, 'spectral': False}
      inc=583.2000 cap=583.2000 diff=-1.14e-13
      effD=1.195 dose=0.1891 arisk=1.95x mu=0.126
```

The descending m-scan accepted the FIRST m at which every guard passed and
walked all the way to m=0: for every m in [1, 972] at least one of
{dose <= 0.05, A-risk <= 1.05x FIXED*, spectral >= 0.5 FIXED*} failed.
The OED-DT arm as deployed would be vacuous (identical in kind to the
balanced fallback), which would empty the primary gate.

## 2. Three convergent measurements of the same structural fact

1. RELAXED_REFERENCE bounds at the dose-feasible designs: 2.01/r (fast),
   3.51/r (safe) — the dose-feasible continuous OED sits far below the
   dose-relaxed optimum.
2. Pure SCAT32-based FIXED_DOSE beats the dose-feasible continuous
   reference outright (safe: eff_D = 2.25 relative to it).
3. The rounded OED design fails dose/A-risk/spectral against the
   LOAD-MATCHED SCAT32 reference at every OED fraction above zero.

Interpretation for R17: under the frozen guard stack (G5 +/-5 percent
per-pixel dose band + A-risk/spectral anchored to the strongest balanced
baseline), the feasible set is so close to dose-uniform, baseline-like
designs that the R11/R13 variable-load OED has no room to differ from the
baseline it is gated against. Options R17 should adjudicate BEFORE any
confirmatory scene is opened:

  (a) widen the dose band for k-sparse atoms (the A5 row bound showed the
      rounding granule is comparable to the band at these loads);
  (b) re-anchor A-risk/spectral to a load-matched but not dose-balanced
      reference, or gate them on the continuous (pre-rounding) design;
  (c) accept the negative finding: dose-safe, guard-compliant OED-DT
      cannot beat balanced scattering under this physics -- and redesign
      or drop the gate arm rather than launch it vacuously.

## 3. Box-10 (full constrained certificate) status

Root-cause chain, fully diagnosed:
  freeze-8: G_full = inf via LP_FAIL "time limit" (HiGHS ground on the
  dense degenerate 2x1024 mu block);
  pixel-subset LP fix -> freeze-9: LP_FAIL "unbounded" (HiGHS Status 10):
  a mu-pixel with no covering atom in the LP candidate set makes
  u-_ja = (1-delta)*zbar_a > 0 for every candidate atom, so raising that
  mu- lowers every reduced sensitivity simultaneously and t -> -inf.
Fix implemented (oed_design_v5.full_constrained_cert): per-pixel coverage
seeding of the candidate set + a finite mu cap (restricting the dual cone
only loosens the bound, so validity is preserved). The R15 active-dose toy
verifies the formula at the constrained optimum to 1.8e-09 (<= 1e-8).
Real-scale re-validation was still running when the vacuous-arm stop rule
fired; it shares its root with Section 2 (the dose-feasible design's
distance from every optimum) and should be finished under whatever
architecture R17 selects.

## 4. What is NOT in question

Boxes 1-9, 12-14 all PASS on the corrected pipeline, including: balanced
52-row pre-scan (exact identities), non-oracle V0, 13-family variable-k
dictionary manifest, FIXED* = SCAT32 by the frozen R10 DEV-PSNR rule
(17.498 dB vs LBLOB16 11.908 vs MATCH1 9.181 -- no tie), 52+972 accounting
for every arm, nine-dwell kernel/trust/bias certification, analyzer
verdict emission on positive AND negative DEV mocks, outcome-blindness,
manifests + 15120-cell expected-cell ledger with FIXED* propagated.
