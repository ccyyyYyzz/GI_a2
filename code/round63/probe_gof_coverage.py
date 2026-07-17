"""ROUND63 outcome-blind GOF-coverage probe (round-4 adjudication evidence).

Question: does the digest-§4 MODEL_FAIL gate have honest coverage? Two gates
are compared on data whose generator IS the assumed model (non-paralyzable
renewal, exact detector params — so MODEL_FAIL firing is by definition a
false alarm), plus a misspecified-arm contrast (POISSON-LIN, dead time
ignored at rho=0.6 — MODEL_FAIL SHOULD fire; the repaired gate must not have
bought coverage by giving up power):

  fixed — literal digest §4: parametric bootstrap band at the FIXED
          cross-fitted lam_hat (counting noise only);
  refit — candidate A: refit-per-replicate bootstrap at eta_min (null includes
          estimation variance + smoothing bias).

Outcome-blind: no PSNR, no truth-based quantity enters any verdict; x_true is
used only to SIMULATE the data. One select_eta_and_fit(gof_mode="refit") call
per (dataset, arm) — its info carries both bands, so the fixed-gate verdict is
re-derived post hoc from the same D-curve (identical fits, no double compute).

Writes results/round63_gof_probe/probe_results.json + PROBE_SUMMARY.md.
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

from gi_core.utils import rng_for
from physics import Detector, simulate_counts
from solvers import ArmContext
from select_eta import select_eta_and_fit

OUT = os.path.join(os.path.dirname(HERE), "..", "results", "round63_gof_probe")
OUT = os.path.abspath(OUT)

SIDE = 32
N = SIDE * SIDE
TAU, NU, RHO = 1.0, 500.0, 0.6
T = NU * TAU
R_CORRECT = 5          # correct-model datasets per M regime
R_CONTRAST = 2         # misspecified-arm (POISSON-LIN) datasets per M regime
M_LIST = (1500, 512)   # over- and under-determined (M/n = 1.46 / 0.5)
B_REFIT = 25
SELECT_ITER = 40


def phantom():
    img = np.full((SIDE, SIDE), 0.03)
    yy, xx = np.mgrid[0:SIDE, 0:SIDE]
    img[(yy - 9) ** 2 + (xx - 9) ** 2 <= 25] = 1.0
    img[20:28, 4:12] = 0.55
    img[6:9, 16:30] = 1.0
    return img.ravel() / img.sum()


def fixed_gate_verdict(info):
    """Re-derive the literal-§4 (fixed-band) MODEL_FAIL from a refit-mode info:
    same D-curve/one-SE threshold, gate = the recorded fixed-lam_hat band."""
    Db, mrb = info["gof"]["D_band"], info["gof"]["mean_r_band"]
    D = info["D_curve"]
    mr = info["mean_r_curve"]
    thr = info["one_se_threshold"]
    E = [i for i in range(len(D))
         if D[i] <= thr and Db[0] <= D[i] <= Db[1] and mrb[0] <= mr[i] <= mrb[1]]
    return (not E), E


def run_one(arm, x_true, M, ds):
    rng_p = rng_for(0, 63, 1, 1, M, ds)             # per-dataset pattern stream
    A = 2.0 * (rng_p.random((M, N)) < 0.5)
    u = A @ x_true
    det = Detector(tau=TAU, dark=0.0, start_mode="active")
    Phi = RHO / (TAU * float(u.mean()))
    b, _ = simulate_counts(u, Phi, T, det, rng_for(0, 63, 3, 9, M, ds))
    ctx = ArmContext(Phi=Phi, det=det, T=T, side=SIDE, n_iter=150,
                     select_iter=SELECT_ITER, pattern_kind="bern50",
                     meta={"kind": "bern50"})
    t0 = time.time()
    _, info = select_eta_and_fit(arm, A, b, ctx, K=5, seed=ds,
                                 gof_mode="refit", B_refit=B_REFIT)
    dt = time.time() - t0
    fail_fixed, E_fixed = fixed_gate_verdict(info)
    rb = info["gof"]["refit"]
    row = {
        "arm": arm, "M": M, "dataset": ds, "runtime_s": round(dt, 1),
        "D_at_eta_min": info["gof"]["D_at_eta_min"],
        "fixed_band": [round(v, 1) for v in info["gof"]["D_band"]],
        "refit_gate": [round(v, 1) for v in rb[0]["gate"]],
        "refit_D_mean": round(rb[0]["mean"], 1),
        "refit_D_sd": round(rb[0]["sd"], 1),
        "MODEL_FAIL_fixed": bool(fail_fixed),
        "MODEL_FAIL_refit": bool(info["MODEL_FAIL"]),
        "eta_star_refit": info["eta_star"],
        "eta_star_fixed": (info["eta_grid"][max(E_fixed)] if E_fixed
                           else info["eta_min"]),
        "overfit_flag": info["overfit_flag"],
        # candidate-B evidence: is the MEAN residual honest under the null and
        # displaced under misspecification (the D band's power without its
        # estimation-variance fragility)?
        "mean_r_at_eta_min": round(info["gof"]["mean_r_at_eta_min"], 4),
        "mean_r_fixed_band": [round(v, 4) for v in info["gof"]["mean_r_band"]],
        "mean_r_refit_gate": [round(v, 4) for v in rb[1]["gate"]],
        "corr_r_rho_at_eta_min": round(
            info["corr_r_rho_curve"][info["i_min"]], 4),
        # candidate-D evidence: discrepancy RATIO D_obs / refit-null mean
        "D_ratio_refit": round(
            info["gof"]["D_at_eta_min"] / rb[0]["mean"], 3),
    }
    print("  %-11s M=%-5d ds=%d D=%.0f fixed=[%.0f,%.0f]->%s "
          "refit_gate=[%.0f,%.0f]->%s eta*=%.3g (%.0fs)"
          % (arm, M, ds, row["D_at_eta_min"], row["fixed_band"][0],
             row["fixed_band"][1],
             "FAIL" if fail_fixed else "pass",
             row["refit_gate"][0], row["refit_gate"][1],
             "FAIL" if row["MODEL_FAIL_refit"] else "pass",
             row["eta_star_refit"], dt), flush=True)
    return row


def main():
    os.makedirs(OUT, exist_ok=True)
    x_true = phantom()
    rows = []
    for M in M_LIST:
        print("[probe] correct-model RQL, M=%d (M/n=%.2f)" % (M, M / N),
              flush=True)
        for ds in range(R_CORRECT):
            rows.append(run_one("RQL", x_true, M, ds))
        print("[probe] misspecified contrast POISSON-LIN, M=%d" % M, flush=True)
        for ds in range(R_CONTRAST):
            rows.append(run_one("POISSON-LIN", x_true, M, ds))

    def rate(arm, key, M=None):
        sel = [r for r in rows if r["arm"] == arm and (M is None or r["M"] == M)]
        return (sum(r[key] for r in sel), len(sel))

    summary = {}
    for M in M_LIST:
        summary["RQL_M%d" % M] = {
            "false_alarm_fixed": "%d/%d" % rate("RQL", "MODEL_FAIL_fixed", M),
            "false_alarm_refit": "%d/%d" % rate("RQL", "MODEL_FAIL_refit", M),
        }
        summary["PLIN_M%d" % M] = {
            "detect_fixed": "%d/%d" % rate("POISSON-LIN", "MODEL_FAIL_fixed", M),
            "detect_refit": "%d/%d" % rate("POISSON-LIN", "MODEL_FAIL_refit", M),
        }
    out = {"config": {"side": SIDE, "tau": TAU, "nu": NU, "rho": RHO,
                      "M_list": list(M_LIST), "B_refit": B_REFIT,
                      "select_iter": SELECT_ITER, "R_correct": R_CORRECT,
                      "R_contrast": R_CONTRAST},
           "rows": rows, "summary": summary}
    with open(os.path.join(OUT, "probe_results.json"), "w") as f:
        json.dump(out, f, indent=1)

    lines = ["# GOF coverage probe (outcome-blind, round-4 adjudication)", "",
             "Correct-model false-alarm rate (MODEL_FAIL on data whose",
             "generator IS the assumed renewal model) and misspecified-arm",
             "detection rate (POISSON-LIN at rho=0.6, dead time ignored):", ""]
    for k, v in summary.items():
        lines.append("- **%s**: %s" % (k, json.dumps(v)))
    lines += ["", "Full rows in probe_results.json. Verdict logic: the honest",
              "gate must have LOW false alarm on the first block and HIGH",
              "detection on the second, in BOTH M regimes."]
    with open(os.path.join(OUT, "PROBE_SUMMARY.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("[probe] summary:", json.dumps(summary), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
