"""Stage-0 entry point — spec §9: UT -> A1 -> A2 -> A3 -> B (+ gates + figures).

Usage: python code/run_all.py [--phase ut|a|b|gates|figures|all]
STOP-level failures abort subsequent phases (spec §0.6); exit codes:
0 ok, 1 UT fail, 2 Phase A STOP, 3 budget abort.
"""
import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable


def run(script, tag):
    t0 = time.time()
    print("\n===== %s =====" % tag, flush=True)
    r = subprocess.call([PY, os.path.join(HERE, script)])
    print("===== %s done rc=%d (%.1fs) =====" % (tag, r, time.time() - t0), flush=True)
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", default="all",
                    choices=["ut", "a", "b", "gates", "figures", "all"])
    args = ap.parse_args()
    t0 = time.time()
    if args.phase in ("ut", "all"):
        rc = run("ut.py", "UNIT TESTS")
        if rc != 0:
            print("UT FAILED - stopping (fix layer).", flush=True)
            return 1
    if args.phase in ("a", "all"):
        rc = run("phase_a.py", "PHASE A")
        if rc != 0:
            print("PHASE A STOP - aborting Phase B per spec.", flush=True)
            run("figures.py", "FIGURES (partial)")
            return 2
    if args.phase in ("b", "all"):
        rc = run("phase_b.py", "PHASE B")
        if rc == 2:
            print("PHASE B budget abort.", flush=True)
    if args.phase in ("gates", "b", "all"):
        run("gates_b.py", "B GATES")
    if args.phase in ("figures", "b", "all", "gates"):
        run("figures.py", "FIGURES")
    print("TOTAL wall %.1fs" % (time.time() - t0), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
