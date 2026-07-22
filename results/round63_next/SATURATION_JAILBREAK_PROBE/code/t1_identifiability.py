"""
T1 -- IDENTIFIABILITY (exact / noiseless).

Verify the product formula and recover Sigma v_p^2 (= p_2) from the swept-power
survival curve on 6 masked bridge scenes.  BAR: relative error < 1% noiseless.

Two estimators reported:
  * log/power-sum linear fit  (L=-logQ = a p1 - a^2 p2/2 + ...)
  * elementary-symmetric poly fit (1/Q = 1 + a e1 + a^2 e2 + ...; p2=e1^2-2e2)
Plus a QUENCHED-vs-ANNEALED check at C=3600 (the no-calibration claim).
"""
import json
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import saturation_core as sc
import gi_operator as op

OUT = os.path.dirname(HERE)
DATA = "data/r63_bridge_scenes"
rng = np.random.default_rng(651_101)

SCENES = ["bridge_contour_0", "bridge_twopop_1", "bridge_microtex_2",
          "bridge_control_0", "bridge_contour_2", "bridge_twopop_3"]

M = op.build_operator()                    # 51 x 1024
L = 10                                      # sweep levels
ORDER = 4                                   # truncation order for the fit

report = {"scenes": {}, "settings": dict(L=L, order=ORDER, C_quenched=3600,
          occ_lo=0.05, occ_hi=3.0, mask_density=0.5, K=op.K)}
rel_log, rel_es, rel_quench = [], [], []

print("scene                 p2_true      p2_log     rel_log    p2_ES      rel_ES"
      "    rel_quench(C=3600)")
for name in SCENES:
    x = op.load_scene(os.path.join(DATA, name + ".npz"))
    m = M[0]                                # one representative 50%-density mask
    v = m * x                               # masked scene
    p_true = sc.power_sums(v, 4)
    p1t, p2t = p_true[1], p_true[2]

    a_lv = sc.default_a_levels(v, L=L, occ_lo=0.05, occ_hi=3.0)
    Q = sc.Q_exp(a_lv, v)                   # NOISELESS ensemble survival

    flog = sc.fit_powersums(a_lv, Q, order=ORDER)
    fes = sc.fit_es_poly(a_lv, Q, order=ORDER)
    e_log = abs(flog["p2"] - p2t) / p2t
    e_es = abs(fes["p2"] - p2t) / p2t

    # quenched: one fixed realised W at C=3600, annealed-model fit -> bias
    C = 3600
    W = sc.draw_W(C, len(v), rng)
    Qq = sc.empirical_survival(a_lv, v, W)
    fq = sc.fit_powersums(a_lv, Qq, order=ORDER)
    e_q = abs(fq["p2"] - p2t) / p2t

    rel_log.append(e_log); rel_es.append(e_es); rel_quench.append(e_q)
    report["scenes"][name] = dict(
        n_active=int((v > 0).sum()), p1_true=float(p1t), p2_true=float(p2t),
        a_lo=float(a_lv[0]), a_hi=float(a_lv[-1]), a_max_times_vmax=float(a_lv[-1] * v.max()),
        p2_log=float(flog["p2"]), relerr_log=float(e_log),
        p2_ES=float(fes["p2"]), relerr_ES=float(e_es),
        p2_quenched=float(fq["p2"]), relerr_quenched_C3600=float(e_q))
    print("%-20s %10.3f %10.3f %9.2e %10.3f %9.2e   %9.2e"
          % (name, p2t, flog["p2"], e_log, fes["p2"], e_es, e_q))

report["summary"] = dict(
    primary_estimator="log/power-sum (correct for real >order-pixel scenes)",
    max_relerr_log=float(np.max(rel_log)),
    max_relerr_ES=float(np.max(rel_es)),
    ES_note="ES-poly valid only for P<=order (unit-test regime); diverges for ~500-px scenes",
    median_relerr_log=float(np.median(rel_log)),
    max_relerr_quenched_C3600=float(np.max(rel_quench)),
    median_relerr_quenched_C3600=float(np.median(rel_quench)),
    BAR="relerr < 1% noiseless",
    PASS_noiseless=bool(np.max(rel_log) < 0.01))

print("\nNOISELESS (annealed) max rel err:  log=%.2e  ES=%.2e   BAR=1e-2  -> %s"
      % (report["summary"]["max_relerr_log"], report["summary"]["max_relerr_ES"],
         "PASS" if report["summary"]["PASS_noiseless"] else "FAIL"))
print("QUENCHED C=3600 (fixed W, annealed-model fit) max rel err = %.2e (median %.2e)"
      % (report["summary"]["max_relerr_quenched_C3600"],
         report["summary"]["median_relerr_quenched_C3600"]))

with open(os.path.join(OUT, "t1_identifiability.json"), "w") as fh:
    json.dump(report, fh, indent=2, default=float)
print("wrote t1_identifiability.json")
