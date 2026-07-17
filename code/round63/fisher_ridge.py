"""ROUND63 theory deliverable: the optimal-load ridge rho*(nu).

Computes exact-count Fisher information about log-lambda per unit time,
I_log(rho; nu) = lambda^2 * I_exact(lambda; T=nu*tau, tau) * tau / T,
versus the CLT (Gaussian-renewal) proxy rho/(1+rho) (per-tau, same units).
Outputs the ridge rho*(nu), the QMLE trust boundary (exact/CLT ratio 0.9),
and the full map: results/round63_theory/fisher_ridge.csv + .json.

Scale-free choice: information about log-lambda makes curves comparable
across flux; tau is the only time unit.
"""
import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from physics import fisher_exact, fisher_gauss_renewal

ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(ROOT, "results", "round63_theory")

TAU = 1.0  # scale-free: everything in units of tau
NUS = [20, 50, 100, 200, 500, 1000, 2000]
RHOS = np.geomspace(0.01, 64.0, 49)


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = []
    ridge = {}
    for nu in NUS:
        T = nu * TAU
        ie_list, ig_list = [], []
        for rho in RHOS:
            lam = rho / TAU
            ie = fisher_exact(lam, T, TAU) * lam ** 2 * TAU / T
            ig = fisher_gauss_renewal(lam, T, TAU) * lam ** 2 * TAU / T
            ie_list.append(ie)
            ig_list.append(ig)
            rows.append([nu, float(rho), ie, ig, ie / max(ig, 1e-300)])
        ie_arr = np.array(ie_list)
        i_pk = int(np.argmax(ie_arr))
        # quadratic refine in log-rho
        if 0 < i_pk < len(RHOS) - 1:
            lr = np.log(RHOS[i_pk - 1:i_pk + 2])
            y = ie_arr[i_pk - 1:i_pk + 2]
            c = np.polyfit(lr, y, 2)
            rho_star = float(np.exp(-c[1] / (2 * c[0])))
        else:
            rho_star = float(RHOS[i_pk])
        ratio = ie_arr / np.maximum(np.array(ig_list), 1e-300)
        below = np.nonzero(ratio < 0.9)[0]
        trust = float(RHOS[below[0]]) if below.size else float("inf")
        ridge[str(nu)] = {"rho_star": rho_star,
                          "I_at_peak": float(ie_arr[i_pk]),
                          "I_asymptote_CLT": 1.0,
                          "qmle_trust_rho_ratio0.9": trust}
        print("nu=%5d  rho*=%7.2f  I_peak=%.4f  trust(0.9)=%s"
              % (nu, rho_star, ie_arr[i_pk],
                 ("%.2f" % trust) if np.isfinite(trust) else "none<=64"),
              flush=True)
    with open(os.path.join(OUT, "fisher_ridge.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["nu", "rho", "I_exact_log", "I_clt_log", "ratio"])
        w.writerows(rows)
    with open(os.path.join(OUT, "fisher_ridge.json"), "w") as f:
        json.dump({"ridge": ridge,
                   "definition": "I about log-lambda per unit time, tau units; "
                                 "CLT proxy = rho/(1+rho)"}, f, indent=2)
    print("DONE fisher_ridge", flush=True)


if __name__ == "__main__":
    main()
