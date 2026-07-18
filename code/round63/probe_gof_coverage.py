"""ROUND63 outcome-blind coverage probe for the F1-FROZEN adequacy audit
(select_eta.py, GPT round-4 rule: AUDIT split + coherent refit bootstrap).

Two questions, both outcome-blind (x_true only ever SIMULATES data; no PSNR or
truth-based quantity enters any verdict):

  NULL      generator == the calibrated model (non-paralyzable renewal, exact
            same detector params). MODEL_FAIL_PREDICTIVE is an exact MC test at
            level 1/40 = 0.025, so over R null runs we expect ~R/40 false
            alarms (0-1 at R=10 per regime).
  DETECT    generator violates the calibrated model in the two S3-relevant
            ways: (a) paralyzable true detector audited as non-paralyzable
            (~12% per-slot mean shift at rho=0.6); (b) tau_err: true tau 20%
            larger than calibrated. The flag SHOULD fire.

History: the round-3 literal rule failed this probe 10/10 (fixed band), our
candidate-A refit failed 8/10 (results/round63_gof_probe/probe_results.json @
commit 109b718); the F1 rule is the round-4 repair. M regimes chosen so the
AUDIT split stays powered (M=1500 -> 300 audit groups; M=768 -> 154 >= 128).

ROLE AFTER ROUND 5: the binary gate was DELETED based on this probe's first
run (7/20 null false alarms, 0/10 detection — see PROBE_F1_SUMMARY.md and
docs/ROUND63_GPT_ROUND5_RULING.md). The probe now records the DESCRIPTIVE
distributions of the audit statistics under null/misspecified generators —
supplement evidence for why the audit is reported without a gate.

Writes results/round63_gof_probe/probe_f1_results.json + PROBE_F1_SUMMARY.md.
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

OUT = os.path.abspath(os.path.join(os.path.dirname(HERE), "..",
                                   "results", "round63_gof_probe"))

SIDE = 32
N = SIDE * SIDE
TAU, NU, RHO = 1.0, 500.0, 0.6
T = NU * TAU
M_LIST = (1500, 768)
R_NULL = 10
R_DETECT = 5


def phantom():
    img = np.full((SIDE, SIDE), 0.03)
    yy, xx = np.mgrid[0:SIDE, 0:SIDE]
    img[(yy - 9) ** 2 + (xx - 9) ** 2 <= 25] = 1.0
    img[20:28, 4:12] = 0.55
    img[6:9, 16:30] = 1.0
    return img.ravel() / img.sum()


def run_one(scenario, x_true, M, ds):
    """scenario: 'null' | 'paralyzable' | 'tau_err'. Estimator/audit always
    calibrated NP with tau=TAU; the TRUE generator varies per scenario."""
    rng_p = rng_for(0, 63, 1, 1, M, ds)
    A = 2.0 * (rng_p.random((M, N)) < 0.5)
    u = A @ x_true
    if scenario == "paralyzable":
        det_true = Detector(tau=TAU, dark=0.0, start_mode="active",
                            paralyzable=True)
    elif scenario == "tau_err":
        det_true = Detector(tau=1.2 * TAU, dark=0.0, start_mode="active")
    else:
        det_true = Detector(tau=TAU, dark=0.0, start_mode="active")
    det_est = Detector(tau=TAU, dark=0.0, start_mode="active")
    Phi = RHO / (TAU * float(u.mean()))
    scen_id = {"null": 0, "paralyzable": 1, "tau_err": 2}[scenario]
    b, _ = simulate_counts(u, Phi, T, det_true,
                           rng_for(0, 63, 3, 9, scen_id, M, ds))
    ctx = ArmContext(Phi=Phi, det=det_est, T=T, side=SIDE, n_iter=150,
                     pattern_kind="bern50", meta={"kind": "bern50"})
    cell_key = (1, int(RHO * 1000), int(NU), M, SIDE)
    t0 = time.time()
    _, info = select_eta_and_fit("RQL", A, b, ctx, cell_key=cell_key, seed=ds)
    dt = time.time() - t0
    au = info["audit"]
    row = {"scenario": scenario, "M": M, "dataset": ds,
           "runtime_s": round(dt, 1), "eta_star": info["eta_star"],
           "AUDIT_STATUS": au["AUDIT_STATUS"],
           "D_obs": round(au.get("D_obs", float("nan")), 1),
           "D_star_mean": round(au.get("D_star_mean", float("nan")), 1),
           "D_star_sd": round(au.get("D_star_sd", float("nan")), 2),
           "D_ratio": round(au.get("D_ratio", float("nan")), 4),
           "q_d": au.get("plugin_upper_rank"),
           "q_mean": au.get("q_mean"),
           "q_corr": au.get("q_corr"),
           "LEAKAGE_SUSPECT": au.get("LEAKAGE_SUSPECT"),
           "MEAN_RESIDUAL_WARN": au.get("MEAN_RESIDUAL_WARN"),
           "LOAD_CORR_WARN": au.get("LOAD_CORR_WARN"),
           "mean_r_obs": round(au.get("mean_r_obs", float("nan")), 4)}
    print("  %-12s M=%-5d ds=%d D_ratio=%-7s q_d=%-7s q_mean=%-7s leak=%s "
          "mwarn=%s (%.0fs)"
          % (scenario, M, ds, row["D_ratio"], row["q_d"], row["q_mean"],
             row["LEAKAGE_SUSPECT"], row["MEAN_RESIDUAL_WARN"], dt),
          flush=True)
    return row


def main():
    os.makedirs(OUT, exist_ok=True)
    x_true = phantom()
    rows = []
    for M in M_LIST:
        print("[probe-F1] NULL (correct model), M=%d" % M, flush=True)
        for ds in range(R_NULL):
            rows.append(run_one("null", x_true, M, ds))
    for scen in ("paralyzable", "tau_err"):
        print("[probe-F1] DETECT %s, M=%d" % (scen, M_LIST[0]), flush=True)
        for ds in range(R_DETECT):
            rows.append(run_one(scen, x_true, M_LIST[0], ds))

    def stats(scen, M=None):
        sel = [r for r in rows if r["scenario"] == scen
               and (M is None or r["M"] == M)]
        drs = sorted(r["D_ratio"] for r in sel)
        return {"n": len(sel),
                "D_ratio_med": drs[len(drs) // 2] if drs else None,
                "extreme_rank_qd": "%d/%d" % (
                    sum(1 for r in sel if r["q_d"] is not None
                        and r["q_d"] <= 0.025 + 1e-12), len(sel)),
                "mean_warn": "%d/%d" % (
                    sum(bool(r["MEAN_RESIDUAL_WARN"]) for r in sel), len(sel))}

    summary = {"null_M1500": stats("null", 1500),
               "null_M768": stats("null", 768),
               "paralyzable": stats("paralyzable"),
               "tau_err_+20%": stats("tau_err"),
               "note": "descriptive distributions only (round-5: no gate); "
                       "extreme_rank_qd is the incidence the DELETED gate "
                       "would have flagged"}
    out = {"config": {"side": SIDE, "tau": TAU, "nu": NU, "rho": RHO,
                      "M_list": list(M_LIST), "R_null": R_NULL,
                      "R_detect": R_DETECT, "rule":
                      "f1_audit_split_coherent_bootstrap"},
           "rows": rows, "summary": summary}
    with open(os.path.join(OUT, "probe_f1_results.json"), "w") as f:
        json.dump(out, f, indent=1)
    lines = ["# F1 audit-statistic distribution probe (outcome-blind)", "",
             "Round-5: NO binary gate; these are descriptive distributions.",
             ""]
    lines += ["- **%s**: %s" % (k, json.dumps(v)) for k, v in summary.items()]
    lines += ["", "History: literal round-3 rule 10/10 false alarms; "
              "candidate-A refit 8/10 (probe_results.json); round-4 gate "
              "7/20 null + 0/10 detection -> gate deleted (round-5 ruling)."]
    with open(os.path.join(OUT, "PROBE_F1_SUMMARY.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("[probe-F1] summary:", json.dumps(summary), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
