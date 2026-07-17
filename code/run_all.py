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


def run(script, tag, extra=()):
    t0 = time.time()
    print("\n===== %s =====" % tag, flush=True)
    r = subprocess.call([PY, os.path.join(HERE, script)] + list(extra))
    print("===== %s done rc=%d (%.1fs) =====" % (tag, r, time.time() - t0), flush=True)
    return r


def a3_flagship_dead():
    import json

    p = os.path.join(os.path.dirname(HERE), "results", "phaseA_gates.json")
    if not os.path.exists(p):
        return False
    with open(p) as f:
        g = json.load(f)
    a3 = g.get("A3")
    return bool(a3) and not a3.get("A3_PASS", True)


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
    budget_abort = False
    if args.phase in ("b", "all"):
        extra = []
        if a3_flagship_dead():
            extra = ["--ext-only"]
            print("A3 = FLAGSHIP_DEAD -> Phase B runs the EXTENDED TABLE ONLY "
                  "(spec §7-A3 downgraded-paper path).", flush=True)
        rc = run("phase_b.py", "PHASE B", extra)
        if rc == 2:
            budget_abort = True
            print("PHASE B budget abort — gates will be evaluated on the "
                  "PARTIAL grid and flagged grid_complete=false.", flush=True)
    if args.phase in ("gates", "b", "all"):
        run("gates_b.py", "B GATES")
    if args.phase in ("figures", "b", "all", "gates"):
        run("figures.py", "FIGURES")
    print("TOTAL wall %.1fs" % (time.time() - t0), flush=True)
    return 3 if budget_abort else 0


if __name__ == "__main__":
    sys.exit(main())
