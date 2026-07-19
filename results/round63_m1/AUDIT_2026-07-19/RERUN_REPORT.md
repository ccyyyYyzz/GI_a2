# M1 independent audit — promotion + in-repo rerun report

Date: 2026-07-19
Prepared by: peer-audit promotion step (read-only outside `results/round63_m1/AUDIT_2026-07-19/`).

## Purpose

Promote the peer auditor's independent M1 recomputation (scripts + outputs +
narrative) from its out-of-repo scratch location into a dated, hashed repo
artifact, and verify that both scripts reproduce the auditor's reported numbers
when rerun inside the project environment against the repo's frozen raw shards.

## Provenance of promoted files

Copied byte-for-byte (sha256-verified, see `SHA256_MANIFEST.txt`, all `OK`) from
`E:\GAN_FCC_WORK\.codex_tmp\` into this directory:

| Promoted file | sha256 (source == copy) |
|---|---|
| `m1_independent_recompute.py` | `0901ddaa092cf037b42394488bdc1693cee51b4573e36fdc2fc5abe70f7ea3b2` |
| `m1_independent_output.json` | `57bb3ca8cdb0487b9950dbd1df483dac9d1ef21119d8d039f1b72dceb10cbbcf` |
| `m1_frozen_and_axis_only.py` | `e71637efef519b65ee8084675897a67ffe3a2f6c99a2001b7f69e8c199b783ad` |
| `m1_frozen_and_axis_only_output.json` | `f9fd195f6a24ccd0f93b09125ac3efddef46d5fc17249d0c0797ad8a3a6b1cb0` |
| `2026-07-19_M1_INDEPENDENT_AUDIT.md` | `c8d5c863a3858c990944a4f0d2ca1c443d0b74808c19ea64fbc0bd703129d648` |

## Environment

- OS: Windows 11 (win32)
- Python: 3.11.5 (Anaconda, MSC v.1916 64-bit) — `python` on PATH (`D:\Anacondar\anaconda3\python`)
- numpy: 1.24.4
- pandas: 1.5.3
- Working directory for the rerun: `D:\GI_another`
- Input data (read-only): `D:\GI_another\results\round63_m1\shards\*.csv`
  (loader excludes `M1_CERT_*`; 5,400 imaging rows = 5 arms x 24 images x 9 dwells x 5 seeds)

Determinism note: both scripts draw all randomness from
`numpy.random.default_rng([20260717, 63, seed_tag])` (PCG64 via SeedSequence,
version-stable), so bootstrap results are reproducible bit-for-bit.

## Path handling (why the verbatim copies run unmodified)

No `_rerun` edit was required; the verbatim copies were run as-is:

- `m1_independent_recompute.py` hardcodes `ROOT = Path(r"D:\GI_another")` (absolute)
  and reads `ROOT/results/round63_m1/shards`, so its input is independent of cwd
  and of the script's own location.
- `m1_frozen_and_axis_only.py` imports its sibling via
  `Path(__file__).resolve().parent / "m1_independent_recompute.py"`, so from the
  promoted location it loads the promoted copy (which in turn reads the same
  frozen shards).

## Commands

```powershell
# cwd = D:\GI_another, python 3.11.5 on PATH
python results\round63_m1\AUDIT_2026-07-19\m1_independent_recompute.py   > results\round63_m1\AUDIT_2026-07-19\m1_independent_output_RERUN.json
python results\round63_m1\AUDIT_2026-07-19\m1_frozen_and_axis_only.py    > results\round63_m1\AUDIT_2026-07-19\m1_frozen_and_axis_only_output_RERUN.json
```

## Wall times (exit code 0, empty stderr for both)

| Script | Wall time | Return code |
|---|---:|---:|
| `m1_independent_recompute.py` (3 x 10,000-replicate bootstraps) | 199.32 s | 0 |
| `m1_frozen_and_axis_only.py` (2 x 10,000-replicate scalar bootstraps) | 10.07 s | 0 |

## Diff table — auditor findings-table numbers vs in-repo rerun

Format: `median / LB2.5 / breadth`. All values reproduce exactly (bit-for-bit).

| # | Analysis | Source | Auditor (promoted output) | Rerun (in-repo) | Match |
|---|---|---|---|---|:--:|
| 1 | Frozen output, RIDGE_OPERATING | `m1_frozen_and_axis_only_output.json` | 1.8669199999999986 / 0.11989000000000072 / 19 | 1.8669199999999986 / 0.11989000000000072 / 19 | OK |
| 2 | Frozen output, RIDGE_SPEED (`nu*rho`) | `m1_frozen_and_axis_only_output.json` | 0.2759082546618772 / 0.17225913049832775 / 1 | 0.2759082546618772 / 0.17225913049832775 / 1 | OK |
| 3 | Axis-only sensitivity (elapsed `nu`) | `m1_frozen_and_axis_only_output.json` | 19.127043091645646 / 16.577852783743403 / 22 | 19.127043091645646 / 16.577852783743403 / 22 | OK |
| 4 | Full spec-conformant RIDGE_OPERATING | `m1_independent_output.json` | 1.86692 / 1.41348975 / 19 | 1.86692 / 1.41348975 / 19 | OK |
| 5 | Full spec-conformant RIDGE_SPEED (elapsed `Topt`) | `m1_independent_output.json` | 19.127043091646133 / 18.328492357080282 / 21 | 19.127043091646133 / 18.328492357080282 / 21 | OK |
| 6 | Full corrected machinery, `nu*rho` descriptive | `m1_independent_output.json` | 0.27864687917406294 / 0.23974692258639926 / 1 | 0.27864687917406294 / 0.23974692258639926 / 1 | OK |

Cross-check against the task's expected 6-dp targets (all pass):
1.866920/0.119890/19, 0.275908/0.172259/1, 19.127043/16.577853/22,
1.866920/1.413490/19, 19.127043/18.328492/21, 0.278647/0.239747/1.

### Mismatches

None. Every reported number in the auditor's findings table reproduces to full
printed precision. A recursive deep-diff over the entire JSON structure of both
outputs (not just the six headline rows — every per-image delta, S_gate, status,
target, and every bootstrap summary field) found zero differing leaves between
the auditor's promoted outputs and the reruns.

### Note on saved-file byte encoding (not a numeric difference)

The rerun JSON files are byte-identical in content to the promoted auditor JSON
files after line-ending normalization, but their raw sha256 differ because the
auditor's saved `*_output.json` use CRLF line endings while the redirect-captured
`*_RERUN.json` use LF:

| File pair | raw sha256 differ? | LF-normalized content sha256 |
|---|:--:|---|
| `m1_independent_output.json` vs `..._RERUN.json` | yes (CRLF vs LF) | identical: `db7c1f7de526d884…` |
| `m1_frozen_and_axis_only_output.json` vs `..._RERUN.json` | yes (CRLF vs LF) | identical: `9eabdcc31942265a…` |

The size deltas (279 bytes over 279 lines; 21 bytes over 21 lines) are exactly the
per-line CR bytes. No numeric field differs.

## Integrity statement

No frozen artifact was modified. All writes in this step were confined to the new
directory `results/round63_m1/AUDIT_2026-07-19/`. The frozen raw shard CSVs under
`results/round63_m1/shards/` were read only. `M1_VERDICTS.json`, `M1_VERDICTS.md`,
`CERT_BRANCH.json`, any ledger, any `.tex`, and all code under `code/round63/`
were not touched. The five promoted source files were copied verbatim
(sha256-verified) and were not edited. No git commit was made.

## Files created in this directory

- `m1_independent_recompute.py` (verbatim copy)
- `m1_independent_output.json` (verbatim copy)
- `m1_frozen_and_axis_only.py` (verbatim copy)
- `m1_frozen_and_axis_only_output.json` (verbatim copy)
- `2026-07-19_M1_INDEPENDENT_AUDIT.md` (verbatim copy)
- `SHA256_MANIFEST.txt` (source vs copy hash verification)
- `m1_independent_output_RERUN.json` (in-repo rerun output)
- `m1_frozen_and_axis_only_output_RERUN.json` (in-repo rerun output)
- `RERUN_REPORT.md` (this file)
