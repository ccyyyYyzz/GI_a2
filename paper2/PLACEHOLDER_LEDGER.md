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
sentences verified absent.** *(Superseded in part by R18 below: the R17 final
claim and §6.2 sentences 1–3 were later removed/amended per ruling #10; the
power-for-time disclosure and the five NOT-PERMITTED sentences carry forward
unchanged.)*

## R18 amendment (2026-07-19, ruling issue #10)

Per `docs/ROUND63_GPT_ROUND18_RULING_RAW.md` (§4 Q3 manuscript correction,
§6 Q5 DEV finding binding) and `docs/ROUND63_METHOD_SPEC_M1_R18_AMENDMENT.md`
§6–§7. No campaign outcome filled; every outcome slot stays `\SPH`.

- **Verdict rename.** `DOSE_SAFE_CERT_PASS` → `FULL_STACK_CERT_PASS`
  everywhere (text, arms-table caption, §4 endpoints, §5.3, Act III figure,
  disclosure ledger, supplement S3/S5). Certified class = the FULL deployed
  conditioning stack C_stack: budget, ±5% dose, A-risk cap 1.05×, spectral
  floor 0.5×V_fix — A-risk/spectral are now GLOBAL design constraints of the
  class, LMI-represented and PSD-priced in the conic dual via a certified
  cutting-plane outer approximation; the R17-era "excluded from the dual
  target class" phrasing is deleted. Comparator-/subspace-relativity sentence
  (ruling §2.2) added in main §3.4 and supplement S3.2.
- **Uniform-dose-collapse claims removed** (ruling §4.1): R17 §6.2 sentence 2
  deleted; sentences 1 and 3 replaced by the amended forms; R17 final claim
  replaced by the §4.4 headline; figure-caption and discussion claims that
  uniform dose itself collapses adaptive concentration deleted. The DEV m=0
  collapse is now attributed to the JOINT guard stack (dose + A-risk +
  spectral), never to dose alone (§3.2 reworked).
- **New Act III four beats** (ruling §4.2): (1) information rewards
  scene-adapted geometry (§3.1 + dose-relaxed OED + probe); (2) uniform dose
  is not enough to eliminate that headroom (§3.3 devprobe); (3) conditioning
  safeguards define the deployable class (§3.4 full-stack certificate);
  (4) global operating-point control is orthogonal and survives (§3.5/§3.7).
- **New DEV-evidence subsection** (main §3.3 `sec:devprobe`, ruling §6):
  mandatory label verbatim; DEV IDs `m1_dev_glyph` seed 0; cells
  (ν,b)=(2000,0.60) and (200,0.05); gap/r 0.776 / 0.903, D-eff ≥ 2.17× /
  2.47× (DEV numbers allowed in prose; source of record
  `results/round63_m1/R18_GAP_PROBE_LOG.txt`); support-preserving-
  initialization caveat (existence only — composition NOT established);
  approved future-work sentence verbatim; `\SPH` slot reserved for the
  frozen six-image replication table (in flight, class R18-DEV).
- **Certificate cell taxonomy** (ruling §5.1): per-cell terminal status
  CERTIFIED / COUNTEREXAMPLE / NUMERICAL_UNRESOLVED; 480/480 CERTIFIED for
  the gate; 420 s per-cell wall cap (60 s primal screen + 180 s solve + one
  deterministic rescaled retry; no third attempt / threshold change /
  imputation / deletion) — methods sentence in §4 endpoints + disclosure-
  ledger rows (solver statuses, wall seconds, dynamic range, scan/cut
  counts, residuals, gaps per r, min generalized eigenvalue, A-risk ratio,
  MU_CAP_ACTIVE). Dose-only no-gate DEV fields listed in supplement S5.3.
- **Handbook** step 4: "full deployed conditioning class" (was "expanded
  dose-safe class"). Act III figure now five panels (added the constructive
  dose-only primal-lower-bound DEV panel; certificate panel labeled
  FULL_STACK_CERT_PASS).

**R18 wording verification (tex-normalized + compiled):** present verbatim —
§4.3 P1 (dose-only probe), §4.3 P2 (full-stack targeting), §4.3 conditional
pass-paragraph (inside the `\SPH` guard with the failure branch noted), §4.4
headline ("is tested for" form — swap to "is certified within 1% local
D-efficiency" ONLY after a passed certificate), §6 mandatory DEV label, §6
future-work sentence, and the R17 power-for-time disclosure (×2). Absent —
R17 §6.2 sentences 1/2/3 as written, R17 final claim, token
DOSE_SAFE_CERT_PASS, "excluded from the dual target class",
uniform-dose-collapse phrasings, and the five R17 NOT-PERMITTED sentences.

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
(supplement proofs / crossover proof), **R18-DEV** (pre-freeze dose-only
replication, frozen set, in flight), **USER** (author, repo URL, funding).

### main_m1.tex (R18 architecture)

| # | line | Location | Placeholder | Class |
|---|------|----------|-------------|-------|
| 1 | ~56  | title page | author block | USER |
| 2 | ~103 | abstract | ridge operating-point + speed + full-stack certificate verdicts + cross-arm numbers | M1 |
| 3 | ~190 | intro (contrib iii) | campaign outcomes (incl. any negative primary) | M1 |
| 4 | ~301 | §2.2 crossover ¶ | uniform two-parameter crossover proof, or prediction label retained | R14-SUPP |
| 5 | ~510 | §3.3 devprobe | frozen six-image dose-only replication table (in flight; set frozen before any aggregate is read) | R18-DEV |
| 6 | ~849 | §5.1 | `RIDGE_OPERATING_PASS`: median ΔQ_ridge (dB), LB2.5, count>0, PASS/FAIL + mandatory disclosures | M1 |
| 7 | ~858 | §5.2 | `RIDGE_SPEED_PASS`: median S_gate, LB, count>1, PASS/FAIL | M1 |
| 8 | ~864 | §5.3 | `FULL_STACK_CERT_PASS`: G_stack/r + terminal-status distribution over 480 cells, fraction ≤1e-2, wall-time vs 420 s cap, PASS/FAIL | M1 |
| 9 | ~871 | §5.4 | context arms (SCAT16/LBLOB16) + dose-relaxed/-constrained OED design diagnostics + DEV collapse summary | M1 |
| 10 | ~879 | §5 Fig `fig:actiii` | Act III five data panels (a guard-vs-α DEV / b dose-only primal lower bound DEV / c G_stack/r + status dist, FULL_STACK_CERT_PASS / d ridge-vs-0.60 images+gains / e resource-corner schematic) | M1 |
| 11 | ~913–917 | §5 cross-arm table | 5 arms × descriptive+paired cols = 21 cells | M1 |
| 12 | ~937 | §6 discussion | frozen-conclusion contingency verdicts (RIDGE_OPERATING/RIDGE_SPEED/FULL_STACK_CERT calls) | M1 |
| 13 | ~944 | §6 discussion | conditional R18 §4.3 pass-paragraph (insert verbatim ONLY if FULL_STACK_CERT_PASS passes; failure branch: report status distribution, no categorical claim; headline "is tested for" → "is certified within 1% local D-efficiency" only on PASS) | M1 |
| 14 | ~1063 | §7 | repo URL wording | USER |
| 15 | ~1064 | §7 | funding / acknowledgments | USER |

Note (retired placeholder, unchanged): the deployed-design occupancy-rung
selection (SCAT32 over LBLOB16 and a matched baseline; DEV radiometric-PSNR
scores 17.498 / 11.908 / 9.181 dB; R15-remedied) is **FILLED inline** in §4 as
prose, not an `\SPH`. Source of record:
`results/round63_m1/FREEZE_CHECKLIST_LEDGER.md`. The R18 dose-only probe
numbers (gap/r 0.776 / 0.903, D-eff ≥ 2.17× / 2.47×, `m1_dev_glyph` seed 0)
are likewise FILLED inline in §3.3 as development evidence (ruling §6 allows
them in prose); source of record `results/round63_m1/R18_GAP_PROBE_LOG.txt`.

### supplement_m1.tex (R18 architecture)

| # | line | Location | Placeholder | Class |
|---|------|----------|-------------|-------|
| S-1 | ~35  | title page | author block (mirrors main; fills with same user decision) | USER |
| S-2 | ~372 | S5.1 | ridge operating-point + nine-dwell speed full results (per-dwell gains, loads, clip flags, ceilings, S_gate curves) | M1 |
| S-3 | ~379 | S5.2 | full-stack certificate distribution (per-cell G_stack/r, terminal statuses, feasibility, D-eff, dual/complementarity, cut/scan counts, wall-time vs 420 s cap, MU_CAP_ACTIVE) | M1 |
| S-4 | ~388 | S5.3 | context arms + dose-relaxed/-constrained OED design diagnostics + DEV collapse (PATH_FEASIBLE_ALPHA); dose-only DEV source-of-record fields listed in prose | M1 |
| S-5 | ~404 | S6 | proofs of Theorems 1–3 + Proposition S1; uniform crossover proof or prediction label | R14-SUPP |

**Counts (pre-unblinding).** 20 open sites (main 15 + supplement 5). Rendered
`\SPH` commands: main 35 (14 single + 21 table cells), supplement 5, total 40.
By class: M1 ×33 (main 30 = 9 single + 21 cells; supp 3), R14-SUPP ×2 (main 1 +
supp 1), R18-DEV ×1 (main), USER ×4 (main 3 + supp 1).

## FILLED at unblinding (2026-07-19, tag `m1-freeze` @ 6f00932)

Verdicts unblinded; all M1 outcome sites and the R18-DEV replication table
filled from the frozen analyzer output and stated CSV aggregations. Compile
verified: `pdflatex`+`bibtex` ×2 on both docs exit 0, 0 undefined, 0 overfull
hbox/vbox. Forbidden certificate-pass paragraph absent; R18 frozen sentences
(§4.3 P1/P2, §4.4 "is tested for" headline, DEV label, future-work, ×2
power-for-time) verified present.

**Sources of record:** `M1_VERDICTS.json`/`.md`, `CERT_BRANCH.json`, the 5
imaging shards `shards/M1_<ARM>_00.csv` (5,400 rows; arm = `shard_id` prefix,
`arm` col = estimator `RQL`), the 480 `shards/M1_CERT_*.csv` cert rows, and
`R18_GAP_PROBE_REPLICATION.md`. Contour diagnostic:
`results/round63_m1/CONTOUR_DIAGNOSTIC.md`.

| # | site | status | fill + provenance |
|---|------|--------|-------------------|
| 2 | abstract | **FILLED (M1)** | primary 1.87 dB / 19-24 / LB 0.12; speed 0.28 neg; cert 0/299/181; 37× dose. `M1_VERDICTS.json` |
| 3 | intro (iii) | **FILLED (M1)** | primary PASS / speed preregistered-negative / cert returns counterexamples |
| 5 | §3.3 devprobe | **FILLED (R18-DEV)** | Table `tab:devreplication`, all 12 rows verbatim from `R18_GAP_PROBE_REPLICATION.md` (gap/r 0.496–0.921, D-eff 1.64×–2.51×); no selection |
| 6 | §5.1 primary | **FILLED (M1, R19 corrected)** | `RIDGE_OPERATING_PASS`=PASS; median ΔQ 1.866920, LB2.5 1.41348975, 19/24 (`CORRECTION_2026-07-19/M1_VERDICTS_SPEC_CORRECTED_R19.json`). |
| 7 | §5.2 speed | **FILLED (M1, R19 corrected)** | `RIDGE_SPEED_PASS`=PASS on elapsed `T_opt`; median S 19.127043091646133, LB 18.328492357080282, 21/24; `nu*rho` remains a post-hoc sensitivity with no verdict. |
| 8 | §5.3 cert | **FILLED (M1, descriptive)** | No categorical certificate endpoint. 0/299/181 over 480 cells + by-anchor; status distribution is reported descriptively with unavailable fields kept N/A. |
| 9 | §5.4 context | **FILLED (M1)** | context ladder mean PSNR_rad @(0.60,ν2000): SCAT32-060 17.12 > SCAT16 14.18 > LBLOB16 9.63 (arm-mean over 24×5); OED diagnostics + DEV collapse referenced to Fig/supp |
| 11 | cross-arm table | **FILLED (M1, R19 corrected)** | @ν2000: corrected ΔQ/S_gate for RIDGE (+1.87/19.127); `nu*rho` sensitivity is disclosed separately. |
| 12 | §6 verdicts | **FILLED (M1)** | three verdict calls: primary PASS, speed FAIL, cert FAIL |
| 13 | §6 conditional | **FILLED (M1, failure branch)** | cert FAILED → pass-paragraph NOT inserted; report 0/299/181 distribution, no categorical/collapse claim; headline kept "is tested for" (NOT swapped) |
| S-2 | S5.1 | **FILLED (M1, R19 corrected)** | per-dwell median ΔQ table plus corrected elapsed-time endpoint; S_gate 19.127043091646133/18.328492357080282/21-24 |
| S-3 | S5.2 | **FILLED (M1)** | full-stack cert distribution: statuses, anchors, wall, dep_feasible, MU_CAP, CE primal_gap 1.19–1.87, solver statuses (299 SKIPPED_COUNTEREXAMPLE, 181 LP_FAIL), d_cert_sha 908cfccbd249de22; descriptive branch |
| S-4 | S5.3 | **FILLED (M1)** | context ladder + OED diagnostics + ADAPTIVE_COLLAPSE_UNDER_GUARDS / PATH_FEASIBLE_ALPHA |

**Ambiguous / unpopulated source columns (disclosed, not guessed):** in the
frozen imaging shards the columns `q_d`, `q_mean` (per-pattern dose/load
quantiles), `audit_status`, `leak_suspect`, `cnr`, `d_ratio`, `k_occupancy`,
`LPIPS` are **empty** for all arms; a physical-peak column is absent. These
were reported as unpopulated/absent rather than fabricated (ceiling frac,
mean J_ex, and load-quantile-adjacent disclosures were instead derived from
the frozen renewal kernel `code/round63/oed_design_v3.kernel_eval` at the
achieved mean load, documented in-caption). `V0_exact` is the RQL
reconstruction variance (identical across arms at matched load), **not** the
paper's `J_ex`; mean `J_ex` was therefore computed from the frozen kernel, not
read from `V0_exact`. The incident-dose ratio uses `S_inc` normalized to the
`SCAT32-060` comparator (documented in `tab:crossarm` caption).

## Remaining open after unblinding pass (7 rendered `\SPH`)

| class | sites |
|-------|-------|
| USER ×4 | main #1 author, #14 repo URL, #15 funding; supp S-1 author |
| R14-SUPP ×2 | main #4 uniform crossover proof (prediction label retained); supp S-5 proofs of Thm 1–3 + Prop S1 |
| pending-figures ×1 | main #10 Act III five-panel figure `fig:actiii` (figure agent owns `paper2/figs/`; currently a boxed `\SPH` placeholder, **not open science** — all underlying panel numbers are filled in text/tables) |

No M1 or R18-DEV science placeholder remains open.

## Fill order (suggested)

1. **Supplement proofs** → S-5 and main #4 (R14 §1.7 "provable immediately"
   list; crossover keeps its prediction label unless the dedicated uniform
   proof lands).
2. **R18-DEV replication** → main #5, as soon as the frozen six-image
   dose-only replication completes (set frozen before aggregates are read;
   every result reported without selection).
3. **m1-freeze → run → analyze** (R18 architecture) → main #2, #3, #6–#13;
   supplement S-2–S-4. Verdicts fill only after the coherent R18 refreeze tag
   and outcome-blind ledger (ruling #10 §7 / amendment §8); the three verdicts
   (`RIDGE_OPERATING_PASS`, `RIDGE_SPEED_PASS`, `FULL_STACK_CERT_PASS`) are
   mutually non-rescuing. Main #13 fills per its conditional rule; if the
   24-cell DEV feasibility gate removes the categorical certificate before the
   tag (ruling §5.4 fallback), #8/#13 fill with the descriptive-only branch.
4. **USER** → main #1, #14, #15; supplement S-1.

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
- **R17→R18 claim discipline (frozen; R18 ruling #10 supersedes the R17
  wording set):** the DEV `m=0` collapse is *development evidence*
  (`ADAPTIVE_COLLAPSE_UNDER_GUARDS`) under the JOINT guard stack, explicitly
  NOT a confirmatory theorem and never attributed to uniform dose alone;
  OED-DT never presented as a confirmatory failure (gate withdrawn pre-freeze
  as infeasible). VERBATIM-present set (R18): §4.3 P1 (dose-only probe
  existence), §4.3 P2 (full-stack targeting), §4.3 conditional pass-paragraph
  (SPH-guarded, failure branch noted), §4.4 headline ("is tested for" form;
  swap to "is certified within 1% local D-efficiency" ONLY after a passed
  certificate — no stronger substitution), §6 mandatory DEV label, §6
  future-work sentence, R17 §2.2 power-for-time disclosure. MUST-NOT-appear
  set: the R17 §6.2 sentences as written, the R17 final claim, token
  DOSE_SAFE_CERT_PASS, "excluded from the dual target class", any claim that
  uniform dose itself collapses adaptive concentration, plus the five R17
  NOT-PERMITTED sentences ("adaptivity is impossible under dead time";
  "safety constraints always eliminate adaptive sensing"; "SCAT32 is globally
  optimal over all physical patterns"; "the m=0 DEV result is a confirmatory
  theorem"; "the ridge arm is photon efficient"). The certificate is always
  stated as comparator- and subspace-relative (preserves SCAT32's predeclared
  conditioning on the directions SCAT32 resolves; not a global optimum over
  images, subspaces, risk tradeoffs, or physical patterns). DEV probe numbers
  are existence-only lower bounds (support-preserving initialization caveat);
  no dose-only verdict may ever be created.
- Bibliography: `paper2/refs.bib` = copy of `paper/refs.bib` + flagged
  METHOD-LINE group; entries marked "VERIFY at submission" need a live
  Crossref pass. R14 §1.9 names Dubi & Atar and He et al. (NIST) as
  adjacent prior art — add + verify if the prior-art sentence is expanded.
