# 2026-07-19 M1 independent, read-only verdict recomputation

Scope: raw M1 imaging shard CSVs under `D:\GI_another\results\round63_m1\shards`.
The repository and frozen results were not modified. The independent implementation
was written before reading `code/round63/m1_score.py`.

## Reproduction artifacts

- Full spec-correct implementation: `E:\GAN_FCC_WORK\.codex_tmp\m1_independent_recompute.py`
- Full output: `E:\GAN_FCC_WORK\.codex_tmp\m1_independent_output.json`
- Independent frozen-semantics reproduction plus axis-only sensitivity:
  `E:\GAN_FCC_WORK\.codex_tmp\m1_frozen_and_axis_only.py`
- Its output: `E:\GAN_FCC_WORK\.codex_tmp\m1_frozen_and_axis_only_output.json`

Commands:

```powershell
python E:\GAN_FCC_WORK\.codex_tmp\m1_independent_recompute.py
python E:\GAN_FCC_WORK\.codex_tmp\m1_frozen_and_axis_only.py
```

## Input integrity

- 5,400 imaging rows = 5 arms x 24 images x 9 dwells x 5 seeds.
- Every `(derived_arm,image,nu)` group has exactly five rows and seed set
  `{0,1,2,3,4}`; zero duplicate `cell_id`; zero nonfinite `PSNR_rad`.
- The CSV `arm` column is always `RQL`; arm identity must be derived from the
  `shard_id` prefix. The independent parser does so at
  `m1_independent_recompute.py:19-24`.
- Every raw `optical_time_s` equals `M*nu*50e-9` to absolute error below `1e-16`.
  Validation output is at `m1_independent_output.json:3-44`.

## Exact numerical results

| Analysis | Median | LB2.5 | Breadth | Verdict |
|---|---:|---:|---:|---|
| Frozen output, RIDGE_OPERATING | 1.8669199999999986 dB | 0.11989000000000072 dB | 19/24 > 0 | PASS |
| Frozen output, RIDGE_SPEED (`nu*rho`) | 0.2759082546618772 | 0.17225913049832775 | 1/24 > 1 | FAIL |
| Axis-only sensitivity: replace `nu*rho` by elapsed `nu`, keep all other frozen-code behavior | 19.127043091645646 | 16.577852783743403 | 22/24 > 1 | PASS |
| Full spec-correct RIDGE_OPERATING | 1.86692 dB | 1.41348975 dB | 19/24 > 0 | PASS |
| Full spec-correct RIDGE_SPEED, elapsed `Topt` proportional to `nu` | 19.127043091646133 | 18.328492357080282 | 21/24 > 1 | PASS |
| Full corrected machinery, descriptive `nu*rho` sensitivity | 0.27864687917406294 | 0.23974692258639926 | 1/24 > 1 | FAIL |

Thresholds: operating median >= 1 dB, LB > 0, breadth >= 18/24; speed median
>= 3, LB > 1, breadth >= 18/24.

## Findings

### Critical: confirmatory speed endpoint uses the wrong horizontal axis

The frozen M1 spec defines `Topt=M_total*nu*tau`
(`ROUND63_METHOD_SPEC_M1.md:16-20`), and R17 transfers the existing nine-dwell
Q90 machinery unchanged (`ROUND63_METHOD_SPEC_M1_R17_AMENDMENT.md:71-77`). R18
again says the empirical endpoints and bootstrap rules are unchanged
(`ROUND63_METHOD_SPEC_M1_R18_AMENDMENT.md:7-9`). The transferred R7 endpoint fits
PAVA in log `Topt` (`ROUND63_GPT_ROUND7_RULING.md:49-62`).

However, `_q90_time` fits on `log(nu*rho)` (`m1_runner.py:457-463`). This is an
incident-exposure/photon-budget coordinate, not elapsed optical time. It flips
RIDGE_SPEED from PASS on the frozen elapsed-time endpoint to FAIL on a different
resource endpoint. The axis-only sensitivity already flips all three bars to PASS.

The manuscript is internally contradictory: it defines optical time as
`M_total*nu*tau` (`paper2/main_m1.tex:750-758`) but later calls `nu*rho`
"optical (photon) time" and labels its FAIL preregistered
(`paper2/main_m1.tex:941-958`; supplement lines 375-382). Recommended repair:
restore elapsed-time RIDGE_SPEED as the confirmatory endpoint, retain `nu*rho`
only as clearly post-hoc/descriptive photon-exposure-to-quality sensitivity, and
remove "preregistered negative" from the latter.

### High: bootstrap is not the frozen nested family-stratified bootstrap

The frozen protocol keeps all six strata and samples four images within each
family after paired within-image resampling of the five seeds, rebuilding curves,
target, and censoring (`ROUND63_GPT_ROUND7_RULING.md:38-41,74-78`;
`ROUND63_GPT_ROUND8_RULING.md:50-53`). Current scoring first averages seeds
(`m1_score.py:89-105`) and `_nested_boot_lb` then resamples the six family labels
themselves, followed by images (`m1_runner.py:422-439`). It never resamples seeds
or rebuilds an endpoint.

The current protocol exactly reproduces the published operating LB 0.11989; 245
of its 10,000 medians are <= 0 because whole declared strata can disappear. The
frozen fixed-strata nested protocol gives LB 1.41348975. Operating remains PASS,
but the reported LB is not the preregistered estimator. Recommended repair: take
raw paired seed arrays, resample five indices per image, rebuild the statistic,
then sample four images within every one of the six fixed families; retain tags
13/14 and 10,000 replicates.

### High: PAVA, Q90 target, and censoring transfer are incomplete/incorrect

The frozen rule requires ordinary equal-weight PAVA, PAVA-derived safe range and
Q90 target, the 0.50 dB `SAFE_RANGE_UNINFORMATIVE` rule, and the full censoring
taxonomy (`ROUND63_GPT_ROUND7_RULING.md:49-71`). Current `_pava` duplicates pooled
weights rather than collapsing blocks (`m1_runner.py:442-454`); for `[3,2,1]` it
returns approximately `[2.145898,2.145898,2.145898]` instead of `[2,2,2]`.
The analyzer derives the safe range from raw endpoints and reduces all failed
crossings to zero (`m1_runner.py:497-504`) without the frozen statuses.

In these data `m1_microtexture_1` is correctly uninformative and must enter with
`S_gate=1`; the frozen analyzer assigns a non-taxonomic value. This explains the
axis-only breadth difference (22 under current behavior versus 21 under the full
frozen taxonomy). Recommended repair: use a tested block-collapse PAVA, derive
Q90 from its safe endpoints, and implement every censoring status explicitly on
every bootstrap replicate.

### Medium: ITT behavior is not robustly enforceable

The current raw imaging data have no missing/nonfinite quality cells, so ITT does
not change this recomputation. Operating maps a missing terminal value to zero
(`m1_runner.py:487-490`). Speed can propagate NaN through `_q90_time`; `NaN` is
truthy in the `Ts and Tf` expression (`m1_runner.py:500-504`), so it is not
explicitly clamped to nonpositive. Imaging shard CSVs also contain no explicit
guard-failure status for the scorer to enforce. Recommended repair: emit a frozen
per-cell terminal status/guard flag, validate unique seeds, and explicitly map any
failed image endpoint to `S_gate=0` before medians/bootstrap.

### Informational: arm mapping and completeness are correct for the present data

The raw `arm` column denotes estimator (`RQL`), not M1 operating arm. Current
assembly correctly derives M1 arm from `shard_id` (`m1_score.py:74-80`), and the
present raw files satisfy the 5-seed design. Add assertions for the exact unique
seed set, not only list length, to prevent duplicate-seed substitution.

## Severity conclusion

The published RIDGE_OPERATING categorical PASS is stable, but its LB is computed
with a non-frozen bootstrap. RIDGE_SPEED is not stable: the frozen elapsed-time
endpoint is a strong PASS, whereas the published FAIL is a `nu*rho` photon-budget
sensitivity labeled as preregistered optical time. This requires a versioned
analysis correction and manuscript correction; frozen raw cells need not be
altered.
