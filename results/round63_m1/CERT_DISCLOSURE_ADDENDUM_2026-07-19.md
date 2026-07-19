# Certificate disclosure addendum (R19, section 5)

Date: 2026-07-19  
Immutable input: tag `m1-freeze` @ `6f00932`  
Ruling: `docs/ROUND63_GPT_ROUND19_RULING_RAW.md` (R19, GitHub issue #11), section 5  
Companion to: `results/round63_m1/CORRECTION_2026-07-19/` (dated correction bundle)

The original frozen certificate rows in `results/round63_m1/shards/M1_CERT_*.csv` remain untouched. This addendum is status-aware: it never synthesizes an unavailable diagnostic. It records which fields are genuinely available per terminal status, marks the rest `N/A`, discloses the hard-coded `arisk_ratio=1.0` as an implementation constant, and reports the true coefficient-range distribution computed from the immutable rows.

Campaign structure: the pre-freeze 24-cell development feasibility screen selected the preregistered `FALLBACK_DESCRIPTIVE` branch, so `FULL_STACK_CERT_PASS` was removed before freeze and the post-freeze full-stack computation is a descriptive three-status analysis, not a PASS/FAIL endpoint.

Frozen status distribution over the 480 full-stack analysis cells (24 confirmatory scenes x 5 seeds x 2 nu x 2 b anchors; source of record):

- CERTIFIED: 0
- COUNTEREXAMPLE: 299
- NUMERICAL_UNRESOLVED: 181
- TOTAL: 480

## Per-status field availability

Field names are the frozen `M1_CERT_*.csv` columns. `avail` = measured and reportable; `N/A` = never computed under the frozen protocol for that status (must not be synthesized). The single distinct raw value is shown in parentheses where the column is constant across the status group.

| Field | CERTIFIED (0 rows) | COUNTEREXAMPLE (299 rows) | NUMERICAL_UNRESOLVED (181 rows) |
|---|---|---|---|
| `status` | none exist | avail (`COUNTEREXAMPLE`) | avail (`NUMERICAL_UNRESOLVED`) |
| `dual_gap_upper_per_r` | -- | **N/A** - dual not executed after primal counterexample (`nan`) | **N/A** - no finite upper bound (`inf`) |
| `primal_gap_lower_per_r` | -- | avail - feasible lower-bound gap (range [1.184952, 1.867441]) | **N/A** - no feasible counterexample (`-inf`) |
| `solver_status_primary` | -- | avail (`SKIPPED_COUNTEREXAMPLE`) | avail (`LP_FAIL: infeasible (HiGHS Status 8)`) |
| `solver_status_retry` | -- | **N/A** - retry not executed (blank) | avail (`LP_FAIL: infeasible (HiGHS Status 8)`) |
| `wall_seconds` | -- | avail (range ~62.8-67.0 s) | avail (range ~63.2-67.3 s) |
| `coefficient_dynamic_range` | -- | **N/A** - dual not executed after primal counterexample (blank) | avail - see distribution below |
| `n_dictionary_scans` | -- | **N/A** (blank) | avail (`0`) - no successful solve |
| `n_arisk_cuts` | -- | **N/A** (blank) | avail (`0`) - no successful solve |
| `n_spectral_cuts` | -- | **N/A** (blank) | avail (`0`) - no successful solve |
| `full_dictionary_scan_residual` | -- | **N/A** (blank) | **N/A** - no numerically valid certificate (`inf`) |
| `min_generalized_eigenvalue` (PSD residual) | -- | **N/A** (blank) | **N/A** - PSD residual `nan` |
| `arisk_ratio` | -- | **N/A** (blank) | implementation constant `1.0` - NOT a measurement (see below) |
| `MU_CAP_ACTIVE` | -- | avail (`False`) | avail (`False`) |
| `dep_feasible` | -- | avail (`True`) | avail (`True`) |
| `dep_dose_dev` | -- | avail (~0.039991) | avail (~0.039991) |
| `d_cert_sha` | -- | avail (`908cfccbd249de22`) | avail (`908cfccbd249de22`) |

Notes.

- **CERTIFIED**: there are no such M1 cells. No synthetic row or template example is populated in the campaign table.
- **COUNTEREXAMPLE** (299 rows): the primal support-expanding screen short-circuited the dual solve. The measured primal fields (`primal_gap_lower_per_r`, `wall_seconds`) and the deterministic feasibility/composition fields (`dep_feasible`, `dep_dose_dev`, `d_cert_sha`, `MU_CAP_ACTIVE`) are reported. Every dual/cut/complementarity/coefficient field not computed before the short circuit is `N/A - dual not executed after primal counterexample`.
- **NUMERICAL_UNRESOLVED** (181 rows): no successful certificate exists (primal LP infeasible, no feasible counterexample and no finite dual upper bound). Failed-attempt metadata (`solver_status_primary/retry`) and `wall_seconds` are reportable; the scan/cut counts are `0`; and the unavailable bound, cut, PSD, complementarity, generalized-eigenvalue, and measured A-risk fields are `N/A - no numerically valid certificate`.

## Hard-coded `arisk_ratio = 1.0`

The `arisk_ratio` column on the 181 NUMERICAL_UNRESOLVED rows is the hard-coded constant `1.0` written in the dual-side output branch of `full_stack_cert_cell`, not a measured A-risk ratio. In `code/round63/oed_design_v5.py`, the CERTIFIED/NUMERICAL_UNRESOLVED output dict ends with:

```python
# code/round63/oed_design_v5.py, full_stack_cert_cell(), lines ~1783-1796
    out.update({
        "status": status,
        "dual_gap_upper_per_r": Gpr,
        ...
        "min_generalized_eigenvalue": eng.get("psd_residual", float("nan")),
        "arisk_ratio": 1.0})        # <-- hard-coded implementation constant
```

A genuinely measured A-risk ratio does exist, but as a distinct quantity: inside the primal screen composition dict (`stack_primal_probe` / `full_stack_cert_cell`, line ~1700) the code computes `"arisk_ratio": arb / h_fix` from the recovered design (`arb = trace(inv(Vb))`). That primal-composition value is a different field and is not the CSV `arisk_ratio` column. The CSV column value of `1.0` must therefore not be plotted, averaged, or interpreted as an observed ratio. Disclosure fields:

```text
original_arisk_ratio        = 1.0
original_arisk_ratio_source = IMPLEMENTATION_CONSTANT_NOT_MEASURED
corrected_arisk_ratio       = N/A
```

## True coefficient-range distribution

The scalar shorthand "4000" is replaced by the actual frozen distribution. The `coefficient_dynamic_range` field is populated on the 181 NUMERICAL_UNRESOLVED rows (it is a dual-engine quantity, `eng["dyn_range"]`); it is `N/A` (blank) on the 299 COUNTEREXAMPLE rows, which short-circuit before the dual engine runs. The distribution below is therefore computed over the 181 rows where the field is present, parsed directly from the immutable `M1_CERT_*.csv` shards with the overflow-merge reader from `code/round63/m1_score.py`:

- N (populated / total cells): 181 / 480
- minimum = 2134.5052083333335
- 25th percentile = 2145.998697916667
- median  = 3999.999999999999  (approximately 4000)
- 75th percentile = 4000.0
- maximum = 5015.934472647808
- mean    = 3120.8673634371034

Rounded for text: minimum 2134.505, median approximately 4000, maximum 5015.934. A large fraction of rows sit exactly at the 4000 cap (median and Q3 both at 4000), which is why the single-value "4000" shorthand was misleading and is replaced by the range and quantiles above. The exact per-row values remain in the immutable shards.

## Deterministic regeneration - possible and deferred

Primal-screen composition diagnostics (for example the `stack_primal_probe` composition fields `budget_used`, `dose_dev`, `gain_mass`, load quantiles, the primal `arb/h_fix` A-risk ratio, and `mu_min`, computed at `oed_design_v5.py` lines ~1684-1701) are deterministically recoverable from the immutable frozen inputs using the exact frozen algorithm and ordering, with no new solver option, random seed, threshold, or iteration count. Such a regeneration is **possible** and is **deliberately deferred**; it is **not run** in this addendum. If it is ever produced, the script, environment, input/output hashes, and deep-diff verification must be archived alongside it.

Unavailable **dual** diagnostics (dual gap, cuts, PSD/complementarity residuals, generalized eigenvalue, measured A-risk) must not be regenerated by running a new or improved certificate solver; that would be a new post-hoc analysis and would have to be separately labeled, not folded into this disclosure.
