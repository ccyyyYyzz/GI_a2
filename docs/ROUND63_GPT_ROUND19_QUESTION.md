# ROUND63 — GPT Round 19 consultation (post-unblinding audit: implementation/spec discrepancies requiring adjudication before any manuscript claim changes)

Status: the M1 campaign ran to completion at tag `m1-freeze` (6f00932): 43
shards / 5,880 cells, sha-reconciled, provenance verified. The frozen
analyzer emitted: `RIDGE_OPERATING_PASS = TRUE` (median 1.867 dB, LB2.5
0.120, 19/24) and `RIDGE_SPEED_PASS = FALSE` (median S 0.276), certificate
descriptive 0/299/181 (`results/round63_m1/M1_VERDICTS.json|.md`).

An independent peer audit (full letter: `docs2/01_FROM_AUDITOR_FINDINGS.md`,
clarification `docs2/03_FROM_AUDITOR_CLARIFICATION.md`; reproduction
artifacts + in-project rerun verification under
`results/round63_m1/AUDIT_2026-07-19/`) then found that the frozen
analyzer deviates from the frozen specification in three ways. Per the
correction discipline: **no frozen artifact was modified, no manuscript
claim has been changed, and we bring the discrepancies here before any
interpretation flips.** The audit also independently confirmed: raw data
completeness/uniqueness, tag integrity, no confirmatory leakage pre-tag,
and all point estimates.

## 1. The discrepancies (evidence chain)

**D1 — Speed-endpoint axis (CRITICAL).** The controlling spec defines
elapsed optical time `T_opt = M_total·ν·τ` (`docs/ROUND63_METHOD_SPEC_M1.md`
§ endpoint; R7 fits Q90 crossings in log T_opt,
`docs/ROUND63_GPT_ROUND7_RULING.md`), and R17 §2.3 / R18 both freeze the
Q90 machinery "unchanged." The implemented analyzer instead fit crossings
in `log(ν·ρ)` (`code/round63/m1_runner.py`, `_q90_time`) — an
incident-photon-exposure coordinate (ν·ρ = λT), not elapsed time. Every raw
`optical_time_s` in the shards equals `M·ν·τ` exactly, so the raw data
carries the spec axis.

Recomputation table (auditor's independent implementation, written before
reading our scorer; verified by in-project rerun):

| analysis | median | LB2.5 | breadth | vs frozen bars (3 / >1 / 18-24) |
|---|---:|---:|---:|---|
| frozen output on ν·ρ (published) | 0.2759 | 0.1723 | 1/24 | FAIL |
| axis-only correction (elapsed ∝ ν) | 19.1270 | 16.5779 | 22/24 | PASS |
| full spec-conformant (axis + bootstrap + PAVA + censoring) | 19.1270 | 18.3285 | 21/24 | PASS |
| corrected machinery kept on ν·ρ | 0.2786 | 0.2397 | 1/24 | FAIL |

The 22→21 breadth difference is entirely the censoring taxonomy:
`m1_microtexture_1` has SAFE dynamic range 0.273 dB < 0.50 dB ⇒
`SAFE_RANGE_UNINFORMATIVE`, S_gate clamped to 1 (corrected PAVA contributes
zero images; median unaffected).

**D2 — Bootstrap conformance (MAJOR).** R7/R8 freeze a nested bootstrap
with paired within-image seed resampling and image resampling **within the
six fixed families**. The implemented `_nested_boot_lb` averages seeds
first and resamples the family labels themselves, so whole strata can
vanish (245/10,000 replicate medians ≤ 0 exactly reproduces the published
operating LB 0.1199). The spec-conformant bootstrap gives operating LB
**1.4135** (PASS unchanged) and elapsed-speed LB **18.3285**.

**D3 — PAVA + censoring conformance (MAJOR).** The reimplemented `_pava`
duplicates pooled weights rather than collapsing blocks (wrong fit for
[3,2,1]); the analyzer omits `SAFE_RANGE_UNINFORMATIVE` and the full frozen
censoring taxonomy of the source endpoint (`pilot_s1.py`). Corrected
machinery is what produces the table's "full spec-conformant" row.

## 2. Questions for R19 (operative ruling requested)

- **Q1 (the adjudication).** Rule on D1: is elapsed `T_opt` the frozen
  speed-endpoint axis (our and the auditor's reading of the spec chain),
  making the published FAIL an implementation artifact to be corrected via
  a dated correction artifact — with the ν·ρ result retained as an
  explicitly post-hoc, descriptive incident-photon-exposure sensitivity
  (an honest photon-cost negative)? Or does the discrepancy render the
  speed endpoint ambiguous, in which case BOTH results are reported
  descriptively and no preregistered speed verdict is claimed? If the
  former: freeze the manuscript wording for (i) the corrected verdict
  sentence, (ii) the mandatory correction-disclosure paragraph (what was
  implemented, what the spec froze, when it was caught, by whom), and
  (iii) the demoted ν·ρ sentence.
- **Q2.** Rule on D2/D3: confirm the spec-conformant machinery (fixed six
  strata, paired seed resampling, block-collapse PAVA, full censoring
  taxonomy) is the frozen analysis, and authorize publishing the corrected
  operating LB 1.4135 (median/breadth unchanged) through the same
  correction-artifact route. State whether the corrected numbers replace
  the originals in the abstract/results with the correction note, or
  appear alongside them.
- **Q3.** Confirm the C2 manuscript semantics: the campaign reports **two
  empirical verdicts + one descriptive structural analysis**; no
  PASS/FAIL/gate/confirmatory language for the certificate anywhere; the
  three-status 0/299/181 distribution stands; freeze the replacement
  sentence for "near-optimality is certified" → "is evaluated rather than
  assumed," and the added sentence after the R18 §4.3 frozen paragraphs
  noting the DEV gate selected the descriptive branch before freeze.
- **Q4.** Confirm the smaller wording boundaries: Act III panel (e)
  comparator must be SCAT32-SAFE (the actual Q90 comparator), no
  "preferred"/"+ certificate" labels; paper 1's uniqueness claim limited to
  the TV/grid evidence shown; "no amount of convex-in-x regularization
  repairs the Gaussian pathology" deleted (false for coercive penalties);
  the mechanism-figure Γ=1 curve labeled a descriptive photon-limited-onset
  proxy per R9 (no predicted benefit boundary).
- **Q5.** Certificate disclosure addendum semantics (auditor M3): the 299
  counterexample rows short-circuit before some §5.3 fields persist; the
  181 unresolved rows carry a hard-coded `arisk_ratio=1.0`. Approve a dated
  addendum marking unavailable diagnostics as N/A (never synthetic),
  reporting the true coefficient-range distribution (2134.5–5015.9,
  median ≈ 4000), and regenerating primal-screen composition diagnostics
  deterministically where possible.
- **Q6 (proofs, R11 pattern).** Please supply, at manuscript rigor, the two
  outstanding supplement proof blocks: (i) the uniform crossover-constant
  proposition (main text site, R14/R16 chain), and (ii) Theorems 1–3 +
  alignment lemma of the paper-2 supplement. We will verify and install
  verbatim with attribution per the established workflow.

Ruling delivery: as a GitHub issue on `ccyyyYyzz/GI_a2`, referencing this
document and the audit directory at this commit.
