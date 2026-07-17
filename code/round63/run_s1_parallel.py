"""ROUND63 S1 pilot — local parallel launcher.

Runs the full S1 grid as a pool of single-cell pilot_s1.py subprocesses
(BLAS pinned to 1 thread each, so N workers scale ~linearly on the i9), then
merges the per-worker CSVs into results/round63_s1/pilot_rows.csv and runs the
standard --analyze-only pass there. Resume-safe at job granularity: a job whose
worker CSV already contains all its expected rows is skipped by pilot_s1's own
per-row resume logic (we simply re-invoke; it exits quickly).

Usage:
  python run_s1_parallel.py                 # full D2 grid, 6 workers
  python run_s1_parallel.py --workers 8
  python run_s1_parallel.py --dry-run       # print the job table only
"""
import argparse
import csv
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
PY = sys.executable
OUT_FINAL = os.path.join(REPO, "results", "round63_s1")

RHO_LIST = [0.05, 0.6]
NU_LIST = [20, 50, 100, 200, 500, 1000, 2000]


def build_jobs():
    """One job per (rho, nu) cell pair — 14 jobs, each its own out dir."""
    jobs = []
    for rho in RHO_LIST:
        for nu in NU_LIST:
            tag = "r%s_n%d" % (str(rho).replace(".", "p"), nu)
            jobs.append({
                "tag": tag, "rho": rho, "nu": nu,
                "out": os.path.join(REPO, "results", "round63_s1_workers", tag),
            })
    return jobs


def run_job(job, extra_args):
    env = dict(os.environ)
    for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
              "NUMEXPR_NUM_THREADS"):
        env[v] = "1"
    cmd = [PY, os.path.join(HERE, "pilot_s1.py"),
           "--rho-list", str(job["rho"]), "--nu-list", str(job["nu"]),
           "--out", job["out"], "--no-analysis"] + extra_args
    t0 = time.time()
    p = subprocess.run(cmd, cwd=HERE, env=env, capture_output=True, text=True)
    dt = time.time() - t0
    tail = "\n".join((p.stdout or "").strip().splitlines()[-3:])
    return job["tag"], p.returncode, dt, tail, (p.stderr or "")[-2000:]


def merge_and_analyze(jobs, extra_args):
    os.makedirs(OUT_FINAL, exist_ok=True)
    merged = os.path.join(OUT_FINAL, "pilot_rows.csv")
    rows, header = [], None
    for job in jobs:
        src = os.path.join(job["out"], "pilot_rows.csv")
        if not os.path.exists(src):
            print("[merge] MISSING %s" % src, flush=True)
            continue
        with open(src, newline="") as f:
            r = csv.reader(f)
            h = next(r)
            if header is None:
                header = h
            elif h != header:
                raise SystemExit("[merge] header mismatch in %s" % src)
            rows.extend(list(r))
    if header is None:
        raise SystemExit("[merge] nothing to merge")
    seen = set()
    dedup = []
    for row in rows:
        key = tuple(row[:8])
        if key in seen:
            continue
        seen.add(key)
        dedup.append(row)
    with open(merged, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(dedup)
    print("[merge] %d rows (%d after dedup) -> %s"
          % (len(rows), len(dedup), merged), flush=True)
    cmd = [PY, os.path.join(HERE, "pilot_s1.py"), "--analyze-only",
           "--out", OUT_FINAL] + extra_args
    p = subprocess.run(cmd, cwd=HERE, text=True)
    return p.returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--merge-only", action="store_true")
    args, extra = ap.parse_known_args()

    jobs = build_jobs()
    if args.dry_run:
        for j in jobs:
            print(j["tag"], "->", j["out"])
        return 0
    if args.merge_only:
        return merge_and_analyze(jobs, extra)

    print("[s1-par] %d jobs, %d workers, extra=%s"
          % (len(jobs), args.workers, extra), flush=True)
    t0 = time.time()
    failed = []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(run_job, j, extra): j for j in jobs}
        done = 0
        for fut in as_completed(futs):
            tag, rc, dt, tail, err = fut.result()
            done += 1
            print("[s1-par] %2d/%d %-10s rc=%d %5.0fs | %s"
                  % (done, len(jobs), tag, rc, dt,
                     tail.replace("\n", " | ")), flush=True)
            if rc != 0:
                failed.append(tag)
                print("[s1-par] STDERR(%s): %s" % (tag, err), flush=True)
    print("[s1-par] sweep wall=%.0fs failed=%s"
          % (time.time() - t0, failed or "none"), flush=True)
    if failed:
        return 1
    return merge_and_analyze(jobs, extra)


if __name__ == "__main__":
    sys.exit(main())
