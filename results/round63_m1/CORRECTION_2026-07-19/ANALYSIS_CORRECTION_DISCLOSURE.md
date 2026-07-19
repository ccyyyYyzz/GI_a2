# Analysis correction disclosure (R19)

Analysis version: `R19_SPEC_CORRECTION_V1`
Date: 2026-07-19
Immutable input: tag `m1-freeze` @ `6f00932`
Ruling: `docs/ROUND63_GPT_ROUND19_RULING_RAW.md` (R19, GitHub issue #11), section 1.4

## Frozen correction-disclosure paragraph (verbatim, ruling section 1.4)

The following paragraph is reproduced verbatim from the R19 ruling and is the frozen text
to be installed in the manuscript. It must not be paraphrased, and the correction must not be
described as a new endpoint, an alternative analysis choice, or a threshold change.

> **Analysis correction.** After the campaign was unblinded, an independent peer audit found that the frozen scoring implementation did not reproduce three parts of the preregistered analysis. It placed the Q90 crossings on $\log(\nu\rho)$, an incident-photon-exposure coordinate, rather than on the frozen elapsed optical time $\log T_{\mathrm{opt}}$; it resampled family labels rather than retaining all six fixed strata with paired within-image seed resampling; and it used a nonstandard isotonic implementation while omitting the frozen censoring taxonomy. The immutable tag, raw cells, and original analyzer outputs were preserved unchanged. A dated correction artifact recomputed the endpoints from the same frozen raw cells using the preregistered axis and machinery; the promoted independent implementation and an in-project rerun agreed on all 18 audited numerical outputs with no deep differences. Scientific claims and tables use the spec-conformant corrected results, while the original outputs remain available in the provenance record.

## Technical detail of the three discrepancies

The three discrepancies below correspond one-to-one, in order, to the three clauses of the
verbatim paragraph. Line references are to the analyzer/spec files as recorded by the
independent auditor (`results/round63_m1/AUDIT_2026-07-19/2026-07-19_M1_INDEPENDENT_AUDIT.md`).

### D1 - wrong speed horizontal axis (critical)

- **Frozen behavior.** The frozen analyzer `_q90_time` fit the Q90 crossings on `log(nu*rho)`
  (`m1_runner.py:457-463`). Because `nu*rho = lambda*T`, this is an incident-photon-exposure
  (photon-budget) coordinate, not elapsed integration time.
- **Specification.** The frozen spec defines elapsed optical time
  `T_opt = M_total * nu * tau` (`ROUND63_METHOD_SPEC_M1.md:16-20`); R17 transfers the nine-dwell
  Q90 machinery unchanged (`ROUND63_METHOD_SPEC_M1_R17_AMENDMENT.md:71-77`); R18 again leaves the
  empirical endpoints/bootstrap rules unchanged (`ROUND63_METHOD_SPEC_M1_R18_AMENDMENT.md:7-9`);
  the transferred R7 endpoint fits PAVA in `log T_opt` (`ROUND63_GPT_ROUND7_RULING.md:49-62`);
  every raw row persists `optical_time_s = M*nu*tau` to absolute error below `1e-16`.
- **Effect.** The wrong axis flips RIDGE_SPEED from PASS on the frozen elapsed-time endpoint to
  FAIL on a different resource endpoint.
- **Resolution.** Restore elapsed `T_opt` as the confirmatory speed endpoint; retain `nu*rho`
  only as a clearly labeled post-hoc incident-exposure-to-quality sensitivity with no verdict.

### D2 - non-frozen bootstrap (high)

- **Frozen behavior.** The frozen scoring first averaged seeds (`m1_score.py:89-105`), then
  `_nested_boot_lb` resampled the six family labels themselves followed by images
  (`m1_runner.py:422-439`); it never resampled seeds nor rebuilt an endpoint. This reproduces the
  published operating LB `0.11989`, and 245 of its 10,000 medians are `<= 0` because whole declared
  strata can disappear.
- **Specification.** Retain all six fixed strata; resample the five seeds pairwise within each
  image; rebuild curves, target, crossing, and censoring inside every replicate; then sample four
  of the four images with replacement within each of the six fixed families; concatenate the six
  four-image strata and take the 24-image median. Family labels are never resampled
  (`ROUND63_GPT_ROUND7_RULING.md:38-41,74-78`; `ROUND63_GPT_ROUND8_RULING.md:50-53`).
- **Effect.** The operating verdict stays PASS, but the reported LB is not the preregistered
  estimator. The frozen fixed-strata nested protocol gives operating LB `1.41348975`.
- **Resolution.** Use the frozen nested family-stratified bootstrap, retaining RNG tags 13
  (operating) / 14 (speed) and `B = 10000`.

### D3 - incomplete PAVA / Q90 target / censoring (high)

- **Frozen behavior.** The frozen `_pava` duplicated pooled weights rather than collapsing blocks
  (`m1_runner.py:442-454`); for `[3, 2, 1]` it returns approximately `[2.145898, 2.145898, 2.145898]`
  instead of `[2, 2, 2]`. The safe range was derived from raw endpoints and failed crossings were
  reduced to zero (`m1_runner.py:497-504`) without the frozen statuses.
- **Specification.** Ordinary equal-weight block-collapse PAVA (required unit test
  `[3,2,1] -> [2,2,2]`); PAVA-derived safe range and Q90 target; the `0.50 dB`
  `SAFE_RANGE_UNINFORMATIVE` rule; and the complete censoring taxonomy (`NORMAL`,
  `FAST_LEFT_CENSORED`, `BOTH_LEFT_CENSORED`, `SAFE_LEFT_FAST_INTERIOR`, `FAST_RIGHT_CENSORED`,
  `SAFE_RANGE_UNINFORMATIVE`, `ANALYSIS_FAILURE`) (`ROUND63_GPT_ROUND7_RULING.md:49-71`).
- **Effect.** `m1_microtexture_1` is correctly uninformative (fitted safe range `0.27336 dB`) and
  must enter with `S_gate = 1`. This censoring step (not PAVA) changes breadth from 22 (axis-only)
  to 21 under the full taxonomy.
- **Resolution.** Use a tested block-collapse PAVA, derive Q90 from its safe endpoints, and
  implement every censoring status explicitly on every bootstrap replicate.

## Corrected source-of-record values

| Endpoint | Median | LB2.5 | Breadth | Verdict |
|---|---:|---:|---|---|
| RIDGE_OPERATING (terminal dQ) | 1.866920 dB | 1.41348975 | 19/24 > 0 | PASS |
| RIDGE_SPEED (elapsed `T_opt`) | 19.127043091646133 | 18.328492357080282 | 21/24 > 1 | PASS |
| Post-hoc `nu*rho` sensitivity | 0.27864687917406294 | 0.23974692258639926 | 1/24 > 1 | no verdict |

The original nonconformant frozen outputs (operating LB `0.11989000000000072`; speed
`median 0.2759082546618772 / LB2.5 0.17225913049832775 / 1-of-24 FAIL` on `log(nu*rho)`) are
preserved unchanged in `results/round63_m1/M1_VERDICTS.json` and are quoted for provenance in
`M1_VERDICTS_SPEC_CORRECTED_R19.json` (`provenance.original_analyzer_outputs_frozen`). The
"axis-only" speed LB `16.577852783743403` (breadth 22/24) preserves D2/D3 and is diagnostic only,
not a manuscript result.

## File pointers

Corrected artifacts (this bundle, `results/round63_m1/CORRECTION_2026-07-19/`):

- `M1_VERDICTS_SPEC_CORRECTED_R19.json` - corrected source of record (numbers, per-image tables,
  full provenance block, raw shard sha256 list).
- `M1_VERDICTS_SPEC_CORRECTED_R19.md` - human-readable companion with the original-vs-corrected
  provenance table.
- `scripts/m1_independent_recompute.py`, `scripts/m1_frozen_and_axis_only.py` - promoted
  independent implementation (copied byte-verbatim from `results/round63_m1/AUDIT_2026-07-19/`,
  originally from `E:\GAN_FCC_WORK\.codex_tmp\`).
- `outputs/m1_independent_output.json`, `outputs/m1_frozen_and_axis_only_output.json` - promoted
  auditor outputs; `outputs/*_RERUN.json` - in-project reruns (all copied byte-verbatim from the
  audit directory).
- `SHA256SUMS` - hashes over every file in this bundle.

Immutable / preserved (unchanged):

- `results/round63_m1/M1_VERDICTS.json`, `results/round63_m1/M1_VERDICTS.md` - original frozen
  analyzer outputs (tag `m1-freeze` @ `6f00932`).
- `results/round63_m1/shards/*.csv` - raw shard cells (read-only; sha256 listed in the corrected JSON).

Audit and verification record (referenced, not copied):

- `results/round63_m1/AUDIT_2026-07-19/2026-07-19_M1_INDEPENDENT_AUDIT.md` - auditor findings.
- `results/round63_m1/AUDIT_2026-07-19/RERUN_REPORT.md` - promotion + in-repo rerun verification
  (18/18 audited numbers matched; zero deep-diff leaves).
- `results/round63_m1/AUDIT_2026-07-19/SHA256_MANIFEST.txt` - byte-for-byte promotion verification.

Certificate semantics:

- `results/round63_m1/CERT_DISCLOSURE_ADDENDUM_2026-07-19.md` - per-status field availability, the
  hard-coded `arisk_ratio=1.0` disclosure, and the true coefficient-range distribution.
