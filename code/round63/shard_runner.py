"""ROUND63 shard runner — spec §5 (S2/S3/S4 全网格分片执行器).

One shard = one self-contained unit of the campaign grid: a manifest listing a
list of grid cells (spec §3), the frozen inputs those cells consume (verified by
sha256 before anything runs), and an output CSV.  Shards are farmed across the
6 parallel routes (Colab pro2 x3 + pro1 x2 + local) and their CSVs merged.

Manifest JSON schema
--------------------
{
  "shard_id":   "S2_rho0.05_a",          # required, string
  "stage":      "S2",                     # S2 | S3 | S4 (spec §3)
  "frozen_inputs": [                      # sha256-verified before any cell
      {"path": "data/r63_images/64/stl_00.png", "sha256": "<hex>"},
      ...                                 # path may be absolute or repo-relative
  ],
  "cells": [                              # each cell = dict of grid coordinates
      {"cell_id": "c0001", "rho_bar": 0.05, "nu": 500, "pattern": "bern50",
       "side": 64, "mn": 0.5, "seed": 0, ...},
      ...
  ],
  "output_csv": "results/round63/S2_shard_a.csv",   # required
  "meta_json":  "results/round63/S2_shard_a_meta.json"  # optional (derived)
}

Execution contract
------------------
For each not-yet-done cell the runner calls run_cell(cell) -> list[dict], each
dict one CSV row of metrics (spec §4).  run_cell is imported from
code/round63/campaign.py, which does NOT exist yet; the import lives inside
main() and fails cleanly (exit 3) so this skeleton is runnable/testable now.

META-as-truth + sanitize (copied from code/phase_b.py): the meta JSON records
which cells are fully appended; it is written only AFTER a cell's rows are
flushed+fsync'd.  On restart, CSV rows whose cell_id is not marked complete
(crash mid-append) are dropped, so pivots are never biased by partial rows and
the shard resumes cleanly.  Budget guard via --wall-budget-s.
"""
import argparse
import csv
import hashlib
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)                       # sibling campaign.py
sys.path.insert(0, os.path.dirname(HERE))      # gi_core (for campaign's use)
ROOT = os.path.dirname(os.path.dirname(HERE))

CELL_FIELD = "cell_id"                          # completion key (sanitize pivot)
DEFAULT_WALL_BUDGET_S = 6 * 3600


def _abspath(p):
    return p if os.path.isabs(p) else os.path.normpath(os.path.join(ROOT, p))


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_frozen_inputs(frozen):
    """Hard gate: every declared frozen input must exist and match its sha256.
    Returns (ok, report_lines).  Entries without a sha256 are warned, not failed
    (nothing to verify against)."""
    ok = True
    report = []
    for ent in frozen:
        if isinstance(ent, str):
            ent = {"path": ent}
        path = _abspath(ent["path"])
        want = (ent.get("sha256") or "").lower().strip()
        if not os.path.exists(path):
            ok = False
            report.append("MISSING frozen input: %s" % path)
            continue
        if not want:
            report.append("WARN no sha256 declared for %s (not verified)" % path)
            continue
        got = _sha256(path)
        if got != want:
            ok = False
            report.append("SHA MISMATCH %s\n    want %s\n    got  %s"
                          % (path, want, got))
        else:
            report.append("ok %s" % path)
    if not frozen:
        report.append("no frozen inputs declared")
    return ok, report


def _fieldnames_from(rows):
    """CSV header: cell_id, shard_id first, then remaining keys in first-seen order."""
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    lead = [k for k in (CELL_FIELD, "shard_id") if k in keys]
    return lead + [k for k in keys if k not in lead]


def append_rows(csv_path, fieldnames, rows):
    """Append dict rows; write header if new; flush + fsync (durable checkpoint)."""
    if not rows:
        return
    new = not os.path.exists(csv_path)
    with open(csv_path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, restval="",
                           extrasaction="raise")
        if new:
            w.writeheader()
        w.writerows(rows)
        f.flush()
        os.fsync(f.fileno())


def sanitize_csv(csv_path, done_cells):
    """Drop rows whose cell_id is not marked complete in the meta (crash mid-append).
    Mirrors code/phase_b.py:sanitize_csv, keyed on the cell_id column."""
    if not os.path.exists(csv_path):
        return
    with open(csv_path, newline="") as f:
        rd = csv.DictReader(f)
        header = rd.fieldnames
        if not header or CELL_FIELD not in header:
            return
        kept, dropped = [], 0
        for r in rd:
            if r.get(CELL_FIELD) in done_cells:
                kept.append(r)
            else:
                dropped += 1
    if dropped:
        print("[sanitize] %s: dropped %d orphan rows" % (csv_path, dropped),
              flush=True)
        with open(csv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=header, restval="",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(kept)


def _load_meta(meta_json, shard_id, stage):
    meta = {"shard_id": shard_id, "stage": stage, "cells_done": {},
            "fieldnames": None, "wall_accum_s": 0.0,
            "aborted_over_budget": False}
    if os.path.exists(meta_json):
        with open(meta_json) as f:
            meta = json.load(f)
        meta.setdefault("cells_done", {})
        meta.setdefault("fieldnames", None)
        meta.setdefault("wall_accum_s", 0.0)
        meta["aborted_over_budget"] = False
    return meta


def _save_meta(meta_json, meta):
    tmp = meta_json + ".tmp"
    with open(tmp, "w") as f:
        json.dump(meta, f, indent=2)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, meta_json)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="ROUND63 shard runner (spec §5): sha-verify frozen inputs, "
                    "run manifest cells via campaign.run_cell, checkpoint CSV+meta.")
    ap.add_argument("--manifest", required=True, help="path to shard manifest JSON")
    ap.add_argument("--wall-budget-s", type=float, default=DEFAULT_WALL_BUDGET_S,
                    help="wall-clock budget (accumulated across resumes); abort "
                         "cleanly when exceeded (default %d)" % DEFAULT_WALL_BUDGET_S)
    args = ap.parse_args(argv)

    with open(_abspath(args.manifest)) as f:
        man = json.load(f)
    shard_id = man["shard_id"]
    stage = man.get("stage", "?")
    cells = man["cells"]
    frozen = man.get("frozen_inputs", [])
    out_csv = _abspath(man["output_csv"])
    meta_json = _abspath(man.get("meta_json",
                                 os.path.splitext(out_csv)[0] + "_meta.json"))
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)

    t0 = time.time()
    cpu0 = time.process_time()
    print("[shard %s] stage=%s cells=%d out=%s"
          % (shard_id, stage, len(cells), out_csv), flush=True)

    # 1) frozen-input sha256 gate (before any compute)
    ok, report = verify_frozen_inputs(frozen)
    for line in report:
        print("[shard %s] %s" % (shard_id, line), flush=True)
    if not ok:
        print("[shard %s] FROZEN-INPUT SHA GATE FAILED - aborting before any cell"
              % shard_id, flush=True)
        return 4

    # 2) META-as-truth: load, sanitize orphan CSV rows, resume
    meta = _load_meta(meta_json, shard_id, stage)
    done = set(meta["cells_done"])
    sanitize_csv(out_csv, done)
    wall_prev = float(meta["wall_accum_s"])
    fieldnames = meta.get("fieldnames")

    # 3) campaign is intentionally absent for now -> clean, testable error path
    try:
        from campaign import run_cell
    except Exception as e:  # ImportError today; AttributeError if symbol missing
        print("[shard %s] campaign.run_cell unavailable: %r" % (shard_id, e),
              flush=True)
        print("[shard %s] shard_runner is a runnable SKELETON. Implement "
              "code/round63/campaign.py exposing run_cell(cell)->list[dict rows]. "
              "Frozen-input sha256 gate PASSED and %d/%d cells were already "
              "complete; nothing was computed this invocation."
              % (shard_id, len(done), len(cells)), flush=True)
        return 3

    # 4) cell loop with per-cell durable checkpoint + budget guard
    aborted = False
    for cell in cells:
        cid = cell[CELL_FIELD]
        if cid in done:
            print("[shard %s] cell %s already done - skip" % (shard_id, cid),
                  flush=True)
            continue
        if wall_prev + (time.time() - t0) > args.wall_budget_s:
            aborted = True
            meta["aborted_over_budget"] = True
            print("[shard %s] WALL BUDGET %.0fs EXCEEDED - aborting (resumable)"
                  % (shard_id, args.wall_budget_s), flush=True)
            break
        tc = time.time()
        rows = list(run_cell(cell))
        for r in rows:                              # inject provenance columns
            r.setdefault("shard_id", shard_id)
            r[CELL_FIELD] = cid
        if fieldnames is None and rows:
            fieldnames = _fieldnames_from(rows)
            meta["fieldnames"] = fieldnames
        append_rows(out_csv, fieldnames or [], rows)
        meta["cells_done"][cid] = {"wall_s": round(time.time() - tc, 2),
                                   "n_rows": len(rows)}
        meta["wall_accum_s"] = round(wall_prev + time.time() - t0, 2)
        _save_meta(meta_json, meta)
        print("[shard %s] cell %s done (%d rows, %.1fs)"
              % (shard_id, cid, len(rows), time.time() - tc), flush=True)

    meta["wall_s"] = round(time.time() - t0, 2)
    meta["wall_accum_s"] = round(wall_prev + time.time() - t0, 2)
    meta["process_cpu_s"] = round(time.process_time() - cpu0, 2)
    _save_meta(meta_json, meta)
    print("[shard %s] finished: %d/%d cells done, wall %.1fs, abort=%s"
          % (shard_id, len(meta["cells_done"]), len(cells), meta["wall_s"],
             aborted), flush=True)
    return 2 if aborted else 0


if __name__ == "__main__":
    sys.exit(main())
