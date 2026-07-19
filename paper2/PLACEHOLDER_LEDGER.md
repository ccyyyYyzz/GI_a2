# Placeholder ledger — paper2 (M1 method line: main_m1.tex + supplement_m1.tex)

Campaign rule (inherited from the 2026-07-13 audit discipline): every `\SPH{}`
in the manuscript pair is indexed here with its fill condition; the pair is
submittable only when this ledger is EMPTY (zero `\SPH` in both tex files).

## Structure note (editorial restructure, 2026-07-19)

The manuscript was restructured to ONE storyline (user decision): Act 1 = the
law (main §2.1–2.2), Act 2 = preregistered verification (main §2.3), Act 3 =
the certified design principle + campaign + operating handbook (main §3–6).
`supplement_m1.tex` was created; the following moved there VERBATIM (moved,
not reworded):

- S1 — framework theorem (R14 Candidate A full generalized-KKT text, incl.
  the frozen discrete-optima caveat; main keeps a one-paragraph version with
  the caveat, §3.4);
- S2 — alignment proposition (R14 Candidate C full statement, residual
  bound, power-law reading, failure modes; main keeps one paragraph built on
  the R14 frozen wording, §3.5);
- S3 — design machinery in full: r=200 subspace + fixed ridge, the full
  relaxed-KW certificate equation, the complete guards paragraph
  (peak/weight/dose/spectral/kernel-trust/A-risk + fail states), the frozen
  dictionary listing, and the R11 §4.3 budget dual-reporting (main keeps one
  paragraph + the machinery-at-a-glance table, §3.2);
- S4 — Monte-Carlo verification detail (grid specs, both-run peak lists,
  ν=200 consistency, k=4/3 full texture, E1–E5 extension note; main keeps
  the pooled result + out-of-sample story + ONE k=4/3 sentence, §2.3);
- S5 — secondary-arm full results (SPH placeholders; main keeps one combined
  SPH subsection §5.4);
- S6 — proofs (SPH; main §2 points to it without its own SPH).

New in main: the 55%-trap opening hook (§1 ¶1), the operating-handbook table
(Table 4, Discussion), and the bench-testable SPCM prediction as the final
Discussion paragraph. All frozen wordings verified present: intro novelty ¶
(R10 Q6), R14 claim-discipline sentences (crossover prediction label,
universality scope, critical-manifold statement), R13 §8(f)
photon-economics wording, A's discrete-optima caveat, C's frozen wording,
Theorems 1–3, Table 1 + provenance.

## R17 re-architecture (2026-07-19, ruling issue #9)

The M1 campaign was re-architected per `docs/ROUND63_GPT_ROUND17_RULING_RAW.md`
and `docs/ROUND63_METHOD_SPEC_M1_R17_AMENDMENT.md` (§F). Structural changes to
the manuscript pair (claims + machinery only — NO campaign outcome filled,
every outcome slot stays `\SPH`):

- **Arm roster.** OED-DT retired as the proposed/gated method (its
  `METHOD_SPEED_PASS` gate withdrawn *pre-freeze as infeasible*, never a
  confirmatory failure); OED-EQLOAD arm dropped; `RIDGE-FIXED` renamed
  `RIDGE-SCAT32`; MATCH1 dropped as an arm (kept only in the deployed-design
  occupancy-rung selection record). New 5-arm roster over the common balanced
  972-row SCAT32 design: **SCAT32-SAFE / SCAT32-060 / RIDGE-SCAT32** +
  **SCAT16 / LBLOB16** descriptive context columns. OED survives only as
  dose-relaxed / dose-constrained *design diagnostics* feeding the Act III
  figure.
- **Endpoints.** Old `METHOD_SPEED_PASS` / `METHOD_DESIGN_PASS` /
  `METHOD_FIXED_DWELL_PASS` removed. New three mutually non-rescuing verdicts:
  primary `RIDGE_OPERATING_PASS` (paired fixed-dwell ν=2000, RIDGE-SCAT32 vs
  SCAT32-060; median ΔQ≥1.0 dB / LB2.5>0 / ≥18/24; mandatory power-for-time
  disclosure), secondary `RIDGE_SPEED_PASS` (nine-dwell Q90, SCAT32-SAFE @0.05
  vs dwell-dependent RIDGE-SCAT32; median S_gate≥3 / LB>1 / 18/24),
  confirmatory structural secondary `DOSE_SAFE_CERT_PASS` (480-cell full
  dose-constrained dual G_full/r≤1e-2 over expanded D_cert=D_load∪D_gain,
  D-eff≥0.99005).
- **Act III (§3–6) four-beat** (R17 §6.1): (1) information wants concentration
  (§3.1); (2) uniform-dose safety collapses it — DEV `ADAPTIVE_COLLAPSE_UNDER_
  GUARDS`, labeled development evidence, not a theorem (§3.2, new); (3)
  simplicity certified by the expanded-class dual (§3.3); (4) global operating
  point survives (§3.4 conjugacy + §3.6 alignment). R11+R17 unified: naive
  equalization killed alignment; adaptive concentration conflicts with dose
  safety; balanced global-flux is the robust middle (§3.6). Four-panel Act III
  figure planned (Fig. `fig:actiii`, §5) with `\SPH` data panels.
- **Certificate machinery** (§3.3 main, S3 supp): relaxed-KW G_rel/r≤1e-3
  replaced by the full dose-constrained dual G_full/r≤1e-2 over the expanded
  D_cert; D_gain gain-coupled family (γ∈{0.2,0.5,1,2,5}, emergent load) added
  to the dictionary listing; A-risk/spectral disclosed separately, not in the
  dual target class.

**Frozen R17 wordings verified present verbatim (PDF+tex checked):** final
claim (ruling §6.3), power-for-time disclosure (§2.2, appears in §4 endpoints
and §5.1 results), and the three §6.2 permitted wordings. **Five NOT-PERMITTED
sentences verified absent.**

## Previously FILLED (unchanged by the restructure)

- R14 theory content (all seven original theory placeholders) — Section 2 is
  R14-RULED; matched crossover stays labeled a PREDICTION per claim
  discipline.
- Sweep verdict, three stages: grid out-of-sample hits (3%/6%; commits
  8a627bb pre-data, 1d8d7aa verdict), pooled zoom+Colab exponent −0.658 vs
  −2/3 (1.3%), constant within 5–10% of 2^(−1/3) (~5% below, noise-edge),
  ≈−8% end-point low bias = open higher-order note; k=4/3 preregistered
  5e51088 / refuted 1a32758. Refined ~55% loss at the deterministic ridge.

## Open placeholders

Classes: **M1** (campaign outputs after tag `m1-freeze`), **R14-SUPP**
(supplement proofs / crossover proof), **USER** (author, repo URL, funding).

### main_m1.tex (R17 architecture)

| # | line | Location | Placeholder | Class |
|---|------|----------|-------------|-------|
| 1 | ~56  | title page | author block | USER |
| 2 | ~100 | abstract | ridge operating-point + speed + certificate verdicts + cross-arm numbers | M1 |
| 3 | ~185 | intro (contrib iii) | campaign outcomes (incl. any negative primary) | M1 |
| 4 | ~296 | §2.2 crossover ¶ | uniform two-parameter crossover proof, or prediction label retained | R14-SUPP |
| 5 | ~766 | §5.1 | `RIDGE_OPERATING_PASS`: median ΔQ_ridge (dB), LB2.5, count>0, PASS/FAIL + mandatory disclosures | M1 |
| 6 | ~775 | §5.2 | `RIDGE_SPEED_PASS`: median S_gate, LB, count>1, PASS/FAIL | M1 |
| 7 | ~780 | §5.3 | `DOSE_SAFE_CERT_PASS`: G_full/r distribution over 480 cells, fraction ≤1e-2, PASS/FAIL | M1 |
| 8 | ~785 | §5.4 | context arms (SCAT16/LBLOB16) + dose-relaxed/-constrained OED design diagnostics + DEV collapse summary | M1 |
| 9 | ~793 | §5 Fig `fig:actiii` | Act III four data panels (a guard-vs-α DEV / b G_full/r dist / c ridge-vs-0.60 images+gains / d resource-corner schematic) | M1 |
| 10 | ~821–825 | §5 cross-arm table | 5 arms × descriptive+paired cols = 21 cells | M1 |
| 11 | ~845 | §6 discussion | frozen-conclusion contingency verdicts (RIDGE_OPERATING/SPEED/DOSE_SAFE_CERT calls) | M1 |
| 12 | ~958 | §7 | repo URL wording | USER |
| 13 | ~959 | §7 | funding / acknowledgments | USER |

Note (site formerly #5, now retired as a placeholder): the deployed-design
occupancy-rung selection (SCAT32 over LBLOB16 and a matched baseline; DEV
radiometric-PSNR scores 17.498 / 11.908 / 9.181 dB; R15-remedied) is **FILLED
inline** in §4 as prose, not an `\SPH` — it justifies why the deployed balanced
design is the k=32 rung. Source of record: `results/round63_m1/FREEZE_CHECKLIST_LEDGER.md`.

### supplement_m1.tex (R17 architecture)

| # | line | Location | Placeholder | Class |
|---|------|----------|-------------|-------|
| S-1 | ~35  | title page | author block (mirrors main; fills with same user decision) | USER |
| S-2 | ~304 | S5.1 | ridge operating-point + nine-dwell speed full results (per-dwell gains, loads, clip flags, ceilings, S_gate curves) | M1 |
| S-3 | ~311 | S5.2 | near-optimality certificate distribution (per-cell G_full/r, feasibility, D-eff, dual/complementarity, MU_CAP_ACTIVE) | M1 |
| S-4 | ~317 | S5.3 | context arms + dose-relaxed/-constrained OED design diagnostics + DEV collapse (PATH_FEASIBLE_ALPHA) | M1 |
| S-5 | ~325 | S6 | proofs of Theorems 1–3 + Proposition S1; uniform crossover proof or prediction label | R14-SUPP |

**Counts.** 18 open sites (main 13 + supplement 5). Rendered `\SPH` commands:
main 33 (12 single + 21 table cells), supplement 5, total 38. By class:
M1 ×32 (main 29 + supp 3), R14-SUPP ×2 (main 1 + supp 1), USER ×4 (main 3 +
supp 1). The former FIXED*/occupancy-rung site is FILLED inline (no `\SPH`).

## Fill order (suggested)

1. **Supplement proofs** → S-5 and main #4 (R14 §1.7 "provable immediately"
   list; crossover keeps its prediction label unless the dedicated uniform
   proof lands).
2. **m1-freeze → run → analyze** (revised R17 architecture) → main #2, #3,
   #5–#11; supplement S-2–S-4. Verdicts fill only after the coherent R17
   refreeze tag and outcome-blind ledger (ruling §7 / amendment §E); the three
   verdicts (`RIDGE_OPERATING_PASS`, `RIDGE_SPEED_PASS`, `DOSE_SAFE_CERT_PASS`)
   are mutually non-rescuing.
3. **USER** → main #1, #12, #13; supplement S-1.

## Frozen-wording constraints (carried; verified surviving the restructure)

- No "first"; the novelty is the *combination* (R10 Q6 frozen wording,
  verbatim in intro ¶4). The certificate is always stated relative to the
  frozen dictionary, local pre-scan estimate, and declared resource model.
- Photon-economics framing is non-universal (R13 §8(f) wording verbatim in
  §3.3): OED-DT and RIDGE-FIXED are conjugate corners; absence of ridge
  atoms in a self-test is a KKT outcome, not a theorem.
- **R14 claim discipline (frozen):** permitted — extensive loss + finite
  long-window optimum; small-jitter `c_v^{-2/3}` scaling; "exponent
  confirmed at 1.3%, constant consistent within 5–10%". The finite-window
  crossover `(6ν)^{1/3}(1+12ν c_v²)^{-1/3}` stays labeled a PREDICTION until
  the uniform two-parameter proof. NOT permitted: "all detector jitter obeys
  the same law"; "0.72 is the universal constant"; "every observed empirical
  optimum is a KKT point"; "fixed flux is universally within a constant of
  OED"; "dead-time randomness has not previously been studied".
- Universality claims scoped to the regular finite-variance iid hidden-hold
  class only (exclusions listed verbatim in §2.2).
- Handbook Γ row is descriptive (companion's computable engagement
  criterion), never a confirmatory gate — G10 discipline.
- **R17 claim discipline (frozen):** the DEV `m=0` collapse is *development
  evidence* (`ADAPTIVE_COLLAPSE_UNDER_GUARDS`), explicitly NOT a confirmatory
  theorem; OED-DT never presented as a confirmatory failure (its gate was
  withdrawn pre-freeze as infeasible). Three §6.2 permitted wordings + the
  §6.3 final claim + the §2.2 power-for-time disclosure appear VERBATIM; the
  five NOT-PERMITTED sentences ("adaptivity is impossible under dead time";
  "safety constraints always eliminate adaptive sensing"; "SCAT32 is globally
  optimal over all physical patterns"; "the m=0 DEV result is a confirmatory
  theorem"; "the ridge arm is photon efficient") must never appear. Outcome
  claims (permitted §6.2 #2/#3, final claim) are stated only as
  campaign-licensed conclusions **contingent on** the `\SPH` verdicts, never
  as asserted results.
- Bibliography: `paper2/refs.bib` = copy of `paper/refs.bib` + flagged
  METHOD-LINE group; entries marked "VERIFY at submission" need a live
  Crossref pass. R14 §1.9 names Dubi & Atar and He et al. (NIST) as
  adjacent prior art — add + verify if the prior-art sentence is expanded.
