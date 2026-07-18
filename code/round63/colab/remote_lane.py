#!/usr/bin/env python3
"""ROUND63 remote lane wrapper -- runs ONE campaign shard on the Colab VM.

This is the per-lane worker launched (six-up) by ``session_driver.sh`` inside a
Colab session. It:

  1. pins the four BLAS thread-count env vars to 1 **before** any numpy import
     (they are set in this process's environment and inherited by the
     shard_runner subprocess, which imports numpy fresh);
  2. runs ``code/round63/shard_runner.py`` for the given manifest as a
     subprocess, resolving paths against ``--bundle-root`` (== repo root layout
     inside the bundle);
  3. writes a heartbeat JSON every ``--heartbeat-interval`` seconds carrying the
     wall epoch, shard id, and progress (cells done + rows), obtained by tailing
     the shard's output CSV (the ``cell_id`` column injected by shard_runner);
  4. exits with the shard_runner return code, so a non-zero exit means the shard
     did not fully complete.

Resumability is inherited from shard_runner (META-as-truth): if the shard's
``*_meta.json`` and CSV are restored into the bundle before this runs (e.g. after
a VM recycle), shard_runner skips completed cells and continues. This wrapper
adds no state of its own beyond the heartbeat file, so re-invoking it on the same
shard simply resumes.

shard_runner exit codes (propagated): 0 all cells done; 2 wall-budget abort
(resumable); 3 campaign.run_cell unavailable; 4 frozen-input sha256 mismatch.
"""
import argparse
import csv
import json
import os
import subprocess
import sys
import time

# --- (1) BLAS pin: MUST precede any numpy import in this process or a child ---
_BLAS_VARS = ("OMP_NUM_THREADS", "MKL_NUM_THREADS",
              "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS")
for _v in _BLAS_VARS:
    os.environ[_v] = "1"


def _abspath(bundle_root, p):
    return p if os.path.isabs(p) else os.path.normpath(os.path.join(bundle_root, p))


def _resolve_manifest(bundle_root, manifest):
    """Accept a full path, a repo-relative path, a bare filename, or a shard id
    (with or without .json) under results/round63/manifests/."""
    cands = []
    if os.path.isabs(manifest):
        cands.append(manifest)
    else:
        cands.append(os.path.join(bundle_root, manifest))
        base = manifest if manifest.endswith(".json") else manifest + ".json"
        cands.append(os.path.join(bundle_root, "results", "round63",
                                  "manifests", base))
    for c in cands:
        if os.path.exists(c):
            return os.path.normpath(c)
    raise SystemExit("manifest not found; tried:\n  " + "\n  ".join(cands))


def _read_manifest(path):
    with open(path) as f:
        return json.load(f)


def _count_progress(csv_path):
    """Tail the shard CSV -> (distinct cells done, data rows). Robust to a CSV
    that does not exist yet or is mid-write."""
    if not csv_path or not os.path.exists(csv_path):
        return 0, 0
    cells = set()
    rows = 0
    try:
        with open(csv_path, newline="") as f:
            rd = csv.DictReader(f)
            if not rd.fieldnames or "cell_id" not in rd.fieldnames:
                # header not flushed yet, or unexpected shape
                return 0, 0
            for r in rd:
                rows += 1
                cid = r.get("cell_id")
                if cid:
                    cells.add(cid)
    except (OSError, csv.Error):
        return len(cells), rows
    return len(cells), rows


def _write_heartbeat(path, shard, csv_path, returncode, pid, n_cells_total):
    cells_done, rows = _count_progress(csv_path)
    hb = {
        "epoch": int(time.time()),
        "shard": shard,
        "cells_done": cells_done,
        "cells_total": n_cells_total,
        "rows": rows,
        "pid": pid,
        "returncode": returncode,
    }
    tmp = path + ".tmp"
    with open(tmp, "w", newline="\n") as f:
        json.dump(hb, f)          # single-line JSON (bash-greppable)
    os.replace(tmp, path)
    return hb


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Run one ROUND63 shard on the Colab VM with BLAS pinned and "
                    "a heartbeat (per-lane worker for session_driver.sh).")
    ap.add_argument("--manifest", required=True,
                    help="shard manifest: full path, repo-relative path, bare "
                         "filename, or shard id (resolved under the bundle's "
                         "results/round63/manifests/).")
    ap.add_argument("--bundle-root", required=True,
                    help="path to the unpacked bundle root (== repo root layout).")
    ap.add_argument("--wall-budget-s", type=float, default=21600.0,
                    help="accumulated wall-clock budget passed to shard_runner "
                         "(default %(default)s = 6 h).")
    ap.add_argument("--heartbeat-dir", default=None,
                    help="dir for <shard>.hb.json (default: "
                         "<bundle-root>/results/round63/heartbeats).")
    ap.add_argument("--heartbeat-interval", type=float, default=60.0,
                    help="seconds between heartbeats (default %(default)s).")
    ap.add_argument("--log-dir", default=None,
                    help="dir for <shard>.log (default: "
                         "<bundle-root>/results/round63/logs).")
    ap.add_argument("--python", default=sys.executable,
                    help="python interpreter for shard_runner (default: this one).")
    args = ap.parse_args(argv)

    bundle_root = os.path.abspath(args.bundle_root)
    manifest_path = _resolve_manifest(bundle_root, args.manifest)
    man = _read_manifest(manifest_path)
    shard = man.get("shard_id", os.path.splitext(os.path.basename(manifest_path))[0])
    n_cells_total = len(man.get("cells", []))
    out_csv = _abspath(bundle_root, man.get("output_csv",
                       "results/round63/shards/%s.csv" % shard))

    hb_dir = args.heartbeat_dir or os.path.join(
        bundle_root, "results", "round63", "heartbeats")
    log_dir = args.log_dir or os.path.join(
        bundle_root, "results", "round63", "logs")
    os.makedirs(hb_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    hb_path = os.path.join(hb_dir, "%s.hb.json" % shard)
    log_path = os.path.join(log_dir, "%s.log" % shard)

    runner = os.path.join(bundle_root, "code", "round63", "shard_runner.py")
    if not os.path.exists(runner):
        raise SystemExit("shard_runner.py not found at %s" % runner)

    cmd = [args.python, runner, "--manifest", manifest_path,
           "--wall-budget-s", "%r" % float(args.wall_budget_s)]
    env = dict(os.environ)  # already carries the pinned BLAS vars

    print("[lane %s] launching: %s" % (shard, " ".join(cmd)), flush=True)
    print("[lane %s] cwd-independent (shard_runner resolves via __file__); "
          "csv=%s" % (shard, out_csv), flush=True)

    _write_heartbeat(hb_path, shard, out_csv, None, None, n_cells_total)
    interval = max(1.0, float(args.heartbeat_interval))
    with open(log_path, "a", buffering=1) as logf:
        logf.write("\n===== lane %s start %s =====\n"
                   % (shard, time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())))
        proc = subprocess.Popen(cmd, env=env, stdout=logf, stderr=subprocess.STDOUT)
        _write_heartbeat(hb_path, shard, out_csv, None, proc.pid, n_cells_total)
        while True:
            try:
                rc = proc.wait(timeout=interval)
                break
            except subprocess.TimeoutExpired:
                hb = _write_heartbeat(hb_path, shard, out_csv, None, proc.pid,
                                      n_cells_total)
                print("[lane %s] heartbeat cells=%d/%d rows=%d"
                      % (shard, hb["cells_done"], n_cells_total, hb["rows"]),
                      flush=True)
        _write_heartbeat(hb_path, shard, out_csv, rc, proc.pid, n_cells_total)
        logf.write("===== lane %s exit rc=%d =====\n" % (shard, rc))

    done, rows = _count_progress(out_csv)
    print("[lane %s] shard_runner rc=%d (%d/%d cells, %d rows) -> %s"
          % (shard, rc, done, n_cells_total, rows,
             {0: "COMPLETE", 2: "BUDGET-ABORT/resumable",
              3: "campaign-missing", 4: "sha-mismatch"}.get(rc, "FAIL")),
          flush=True)
    return rc


if __name__ == "__main__":
    sys.exit(main())
