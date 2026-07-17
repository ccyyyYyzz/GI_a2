"""ROUND63 S1 pilot (local) — spec §5. Purpose: (i) verify the acceleration
signal exists and calibrate its magnitude before freezing the full-grid
protocol; (ii) calibrate per-cell runtime for the Colab shard budget.

Part A  equal-time rho scan: side=64, bern50, nu=500, M=2048,
        rho_bar in {0.05, 0.2, 3/7, 1, 2}, seeds {0,1,2}, pilot-8 images,
        arms {GI, POISSON-LIN, SAT-POISSON, PRECORRECT, QMLE}.
Part B  M scan for time-to-quality: rho in {0.05 (safe), 1 (fast)} x
        M in {512, 1024, 2048, 4096, 8192}; safe arm = POISSON-LIN@0.05
        (its home regime) + QMLE@0.05; fast arm = QMLE@1; seeds {0,1,2}.
        S_Q extracted in analysis by PSNR-vs-logM monotone fit.

Pilot results are S1-exploratory by preregistration: protocol knobs may be
tuned HERE; the S2+ full grid is frozen afterwards and never retuned.
"""
import csv
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from campaign import run_cell

ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(ROOT, "results", "round63_s1")

ARMS_A = ["GI", "POISSON-LIN", "SAT-POISSON", "PRECORRECT", "QMLE"]


def cells():
    cs = []
    for rho in [0.05, 0.2, 3.0 / 7.0, 1.0, 2.0]:
        for seed in (0, 1, 2):
            cs.append(dict(part="A", side=64, pattern="bern50", rho_bar=rho,
                           nu=500, M=2048, seed=seed, arms=ARMS_A,
                           images="pilot8"))
    for rho, arms in ((0.05, ["POISSON-LIN", "QMLE"]), (1.0, ["QMLE"])):
        for M in (512, 1024, 2048, 4096, 8192):
            if rho == 0.05 and M == 2048:
                pass  # also covered in Part A but arms differ; keep for pairing
            for seed in (0, 1, 2):
                cs.append(dict(part="B", side=64, pattern="bern50", rho_bar=rho,
                               nu=500, M=M, seed=seed, arms=arms,
                               images="pilot8"))
    return cs


def main():
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, "s1_pilot.csv")
    done = set()
    if os.path.exists(path):
        with open(path) as f:
            for r in csv.DictReader(f):
                done.add((r["part"], float(r["rho_bar"]), int(r["M"]),
                          int(r["seed"])))
    t0 = time.time()
    writer = None
    f = open(path, "a", newline="")
    try:
        for cell in cells():
            key = (cell["part"], float(cell["rho_bar"]), cell["M"], cell["seed"])
            if key in done:
                continue
            t_c = time.time()
            rows = run_cell(cell)
            for r in rows:
                r["part"] = cell["part"]
            if writer is None:
                new = os.path.getsize(path) == 0
                writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
                if new:
                    writer.writeheader()
            writer.writerows(rows)
            f.flush()
            os.fsync(f.fileno())
            print("[s1] %s rho=%.3g M=%d seed=%d done (%.1fs, total %.0fs)"
                  % (cell["part"], cell["rho_bar"], cell["M"], cell["seed"],
                     time.time() - t_c, time.time() - t0), flush=True)
    finally:
        f.close()
    print("S1 PILOT DONE wall=%.0fs" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
