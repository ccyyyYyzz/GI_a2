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

### main_m1.tex

| # | line | Location | Placeholder | Class |
|---|------|----------|-------------|-------|
| 1 | ~56  | title page | author block | USER |
| 2 | ~94  | abstract | primary/secondary verdicts + cross-arm numbers | M1 |
| 3 | ~173 | intro (contrib iii) | campaign outcomes (incl. any negative primary) | M1 |
| 4 | ~284 | §2.2 crossover ¶ | uniform two-parameter crossover proof, or prediction label retained | R14-SUPP |
| 5 | ~591 | §4 arms | committed FIXED* selection record {SCAT32,LBLOB16,MATCH1} | M1 |
| 6 | ~653 | §5.1 | METHOD_SPEED_PASS: median S, LB, count>1, PASS/FAIL | M1 |
| 7 | ~657 | §5.2 | METHOD_DESIGN_PASS: median design gain (dB), LB, count | M1 |
| 8 | ~661 | §5.3 | METHOD_FIXED_DWELL_PASS: median terminal gain (dB), count, LB | M1 |
| 9 | ~665 | §5.4 | secondary arms + diagnostics summary (full tables → supp S5) | M1 |
| 10 | ~682–687 | §5 cross-arm table | 6 arms × 7 columns = 42 cells | M1 |
| 11 | ~797 | §7 | repo URL wording | USER |
| 12 | ~798 | §7 | funding / acknowledgments | USER |

### supplement_m1.tex

| # | line | Location | Placeholder | Class |
|---|------|----------|-------------|-------|
| S-1 | ~35  | title page | author block (mirrors main; fills with same user decision) | USER |
| S-2 | ~277 | S5.1 | RIDGE-FIXED full results (per-dwell gains, loads, clip flags, ceilings) | M1 |
| S-3 | ~282 | S5.2 | OED-EQLOAD kernel ablation full results | M1 |
| S-4 | ~287 | S5.3 | certificate/rounding diagnostics per designed cell | M1 |
| S-5 | ~293 | S6 | proofs of Theorems 1–3 + Proposition S1; uniform crossover proof or prediction label | R14-SUPP |

**Counts.** 17 sites (main 12 + supplement 5); raw `\SPH` commands: main 53
(11 single + 42 table cells), supplement 5, total 58. By class: M1 ×11,
R14-SUPP ×2, USER ×4.

## Fill order (suggested)

1. **Supplement proofs** → S-5 and main #4 (R14 §1.7 "provable immediately"
   list; crossover keeps its prediction label unless the dedicated uniform
   proof lands).
2. **m1-freeze → run → analyze** → main #2, #3, #5–#10; supplement S-2–S-4.
   Main #5 (FIXED*) fills at freeze from the DEV selection record —
   NOTE: R15 (commit af6d0b6) invalidated the previous FIXED* selection;
   the record must be the R15-remedied one.
3. **USER** → main #1, #11, #12; supplement S-1.

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
- Bibliography: `paper2/refs.bib` = copy of `paper/refs.bib` + flagged
  METHOD-LINE group; entries marked "VERIFY at submission" need a live
  Crossref pass. R14 §1.9 names Dubi & Atar and He et al. (NIST) as
  adjacent prior art — add + verify if the prior-art sentence is expanded.
