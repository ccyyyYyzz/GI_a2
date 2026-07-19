# Placeholder ledger — paper2/main_m1.tex (M1 method-line skeleton)

Campaign rule (inherited from the 2026-07-13 audit discipline): every `\SPH{}`
in `main_m1.tex` is indexed here with its fill condition; the manuscript is
submittable only when this ledger is EMPTY (zero `\SPH` in the tex).

## FILLED since the previous ledger revision (R14 landed 2026-07-19)

The R14 unification ruling (`docs/ROUND63_GPT_ROUND14_RULING_RAW.md`, commit
8a627bb) filled all seven theory-content placeholders. Section 3 is now
R14-RULED, not conjectural:

- abstract theory sentence — exact identity + `J_inf(ρ,c_v)` + cubic optimum
  + 2^{-1/3}c_v^{-2/3} law + smooth crossover (min{} conjecture SUPERSEDED);
- intro contribution (iii) — the ruled theory summary + spine assignment;
- §3 intro — spine: B headline theorems / A framework / C conditional
  proposition; deterministic detector = zero-jitter critical manifold;
- §3.1 — Theorems 1–3 stated (identity, long-window rate, exact cubic with
  small-jitter expansion; coefficient 2^{-1/3} = 0.794, NOT the fitted 0.72);
  matched crossover `(6ν)^{1/3}(1+12ν c_v²)^{-1/3}` stated as a
  matched-asymptotic PREDICTION per the frozen claim discipline;
  universality class scoped to regular finite-variance iid hidden holds;
- Table 1 — exact-law prediction column added (5.69 / 3.50 / 1.63 vs measured
  5.43 / 3.30 / ≤2; 5–6%, grid-limited) + registered out-of-sample rows
  (c_v=0.02 → ≈10.4, c_v=0.2 → ≈2.30; `docs/R14_PREREGISTERED_PREDICTIONS.md`,
  commit 8a627bb);
- §3.2 — Candidate A as framework theorem, with the frozen caveat that the
  empirical k≈32 / p≈1 optima are NOT claimed as continuous KKT theorems;
- §3.3 — Proposition 1 (exact alignment condition `ℓ_i J'(κu_i)=β` +
  projected-KKT-residual bound); "without alignment fixed flux can be
  arbitrarily bad" (no universal constant factor);
- discussion hardware warning — exact cubic + coefficient now cited, SPH
  removed.

The former §3.1 site "M1: confirmatory jitter-ridge crossover" was subsumed
into the §5.7 refined-sweep verdict site (#13 below).

## Open placeholders

Four fill classes remain:

- **R14-SUPP** — the supplement carrying proofs at manuscript rigor (R14 §1.7
  lists what is provable immediately; the uniform two-parameter crossover
  proof needs a dedicated triangular-array argument or the prediction label
  stays).
- **M1** — campaign outputs after tag `m1-freeze` (spec
  `docs/ROUND63_METHOD_SPEC_M1.md`; freeze audit R13 §10 checklist).
- **SWEEP** — the refined jitter sweep (jitter_sfi_v2; estimator, hold
  families, exponent tests E1–E5 frozen in R14 §5; predictions registered at
  commit 8a627bb BEFORE completion).
- **USER** — author/affiliation, repository URL wording, funding text.

| # | tex line | Location | Placeholder | Fills from | Class |
|---|----------|----------|-------------|-----------|-------|
| 1 | ~48  | title page | author block | user decision | USER |
| 2 | ~90  | abstract | primary/secondary verdicts + cross-arm numbers (frozen) | M1 campaign | M1 |
| 3 | ~167 | intro (contrib iv) | campaign outcomes (incl. any negative primary) | M1 campaign | M1 |
| 4 | ~336 | §3 intro | supplement: proofs of Theorems 1–3 + Proposition 1 | supplement writing | R14-SUPP |
| 5 | ~427 | §3.1 crossover ¶ | uniform two-parameter crossover proof, or prediction label retained | dedicated proof effort | R14-SUPP |
| 6 | ~623 | §Campaign (arms) | committed FIXED* selection record {SCAT32,LBLOB16,MATCH1} | M1 (DEV selection, R13 §3) | M1 |
| 7 | ~684 | §5.1 Results | METHOD_SPEED_PASS: median S, LB, count>1, PASS/FAIL | M1 campaign | M1 |
| 8 | ~688 | §5.2 Results | METHOD_DESIGN_PASS: median design gain (dB), LB, count | M1 campaign | M1 |
| 9 | ~692 | §5.3 Results | METHOD_FIXED_DWELL_PASS: median terminal gain (dB), count, LB | M1 campaign | M1 |
| 10 | ~696 | §5.4 Results | RIDGE-FIXED per-dwell gain, realized loads, ceiling fractions | M1 campaign | M1 |
| 11 | ~700 | §5.5 Results | OED-EQLOAD kernel ablation (geometry diff + noise disclosure) | M1 campaign | M1 |
| 12 | ~705 | §5.6 Results | relaxed-KW gap, exact D-efficiency, A-risk, spectral margin | M1 campaign | M1 |
| 13 | ~709 | §5.7 Results | refined-sweep verdict: 1/3 and −2/3 exponent tests, cubic residual, scaling collapse vs Eq.(crossover), shape universality, registered out-of-sample rows (0.02→10.4, 0.2→2.30) | refined sweep (E1–E5) | SWEEP |
| 14 | ~730–735 | §5 cross-arm table | 6 arms × 7 columns = 42 cells | M1 campaign | M1 |
| 15 | ~792 | §7 Reproducibility | repo URL wording | user decision | USER |
| 16 | ~793 | §7 Reproducibility | funding / acknowledgments | user decision | USER |

**Counts.** 16 placeholder sites (site 14 holds 42 `\SPH{M1}` cells, so the
raw `\SPH` command count in the tex is 57). By class: R14-SUPP ×2, M1 ×10,
SWEEP ×1, USER ×3.

## Fill order (suggested)

1. **Supplement** → sites 4, 5 (write the supplement with the R14 §1.7
   "provable immediately" list; the crossover either gets the dedicated
   uniform proof or keeps its prediction label — R14 claim discipline).
2. **Refined sweep completes** → site 13 (E1–E5 verdicts; out-of-sample rows
   scored against the registered predictions).
3. **m1-freeze → run → analyze** → sites 2, 3, 6, 7–12, 14. Site 6 (FIXED*)
   fills at freeze from the DEV selection record, before confirmatory scenes
   open.
4. **USER** → sites 1, 15, 16 (author block, repo URL, funding).

## Frozen-wording constraints (carried from paper-1 discipline + R14)

- No "first"; the novelty is the *combination* (R10 Q6 frozen wording, inlined
  in intro ¶4). The certificate is always stated relative to the frozen
  dictionary, local pre-scan estimate, and declared resource model.
- Photon-economics framing is non-universal: OED-DT (photon-budgeted) and
  RIDGE-FIXED (time-budgeted) are conjugate corners; absence of ridge atoms in
  a self-test is a KKT outcome, not a theorem (R13 §8(f)).
- **R14 claim discipline (frozen):** permitted now — extensive loss + finite
  long-window optimum; small-jitter `c_v^{-2/3}` scaling. The full
  finite-window crossover `(6ν)^{1/3}(1+12ν c_v²)^{-1/3}` stays labeled a
  PREDICTION until the uniform two-parameter proof. NOT permitted: "all
  detector jitter obeys the same law"; "0.72 is the universal constant";
  "every observed empirical optimum is a KKT point"; "fixed flux is
  universally within a constant of OED"; "dead-time randomness has not
  previously been studied".
- Universality claims scoped to the regular finite-variance iid hidden-hold
  class only (not heavy-tailed, correlated, observed, paralyzable, event-time,
  or λ-dependent-recovery cases).
- Bibliography: `paper2/refs.bib` is a COPY of `paper/refs.bib` (GROUPs 1–7
  unchanged, already Crossref-verified) plus a flagged METHOD-LINE group.
  Entries marked "VERIFY at submission" (kiefer1959optimum pages,
  nikolov2022volume vol/pages, feng2018adaptive & maxu2021sensing & burger2021
  authors) need the same live Crossref pass before submission. R14 §1.9 also
  names Dubi & Atar (10.1016/j.net.2019.04.015) and He et al. (NIST) as
  adjacent prior art — add + verify if the prior-art paragraph is expanded.
