#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
EXPLORATORY / DEV-ONLY.  NOT preregistered.  No commits.  CPU.

E2 step (a)+(b): build the non-OU gain-path classes and show the FROZEN OU-model
DLGI misfits them -- one-step-ahead innovation whiteness FAILS (Ljung-Box p small,
innovation variance off 1) and the single-exponential (t_c, CV) summary is biased,
where on true OU it passes.  This is the held-out model-adequacy failure the learned
temporal prior must then repair.

Evaluated on the 6 held-out BRIDGE scenes x several gain seeds per class, using the
TRUE scene for z (isolates the temporal-gain model adequacy from recon error).
"""
import os, sys, json, time
from datetime import datetime, timezone
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import e2_common as E
import dl_tool_common as T
import dl_common as DL

OUT = "D:/GI_another/results/round63_next/DL_TOOL_EXPLORATORY"
T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

CLASSES = ["ou", "regime", "heavytail", "quasiper"]
SEEDS = list(range(6))
SCENES = DL.SCENES


def relerr(h, t): return float(abs(h - t) / max(abs(t), 1e-9))


def main():
    rng = np.random.default_rng(7_000_000)
    rows = []
    for cls in CLASSES:
        recs = []
        for sc in SCENES:
            x = DL.scene_x[sc]
            for seed in SEEDS:
                r2 = np.random.default_rng(hash((cls, sc, seed)) % (2**32))
                a_time, meta = E.sample_class(r2, cls)
                Y = E.simulate(x, a_time, r2)
                zc, R, valid, mu_hat = E.z_from_true_scene(Y, x)
                tc_hat, sig_hat = T.mom_autocov(zc, valid)
                if not np.isfinite(tc_hat):
                    tc_hat, sig_hat = 1e3, 1e-3
                innov = E.ou_innovations(zc, R, valid, tc_hat, sig_hat)
                istd = float(innov.std())
                ikurt = float(((innov - innov.mean()) ** 4).mean() / max(innov.var() ** 2, 1e-12) - 3.0)
                Q, p, K = E.ljung_box(innov, K=20)
                tru = E.true_summary(a_time)
                cv_hat = float(T.cv_of_sigma_l(sig_hat))
                # learned-autocovariance medium readout (integrated tau, model-free)
                tau_emp, cv_emp = E.emp_autocov_tau(zc, valid)
                recs.append(dict(
                    p=p, istd=istd, ikurt=ikurt, passed=E.adequacy_pass(p, istd),
                    tc_hat=float(tc_hat), cv_hat=cv_hat,
                    tc_true=tru["tc"], tau_true=tru["tau"], cv_true=tru["cv"],
                    tau_emp=tau_emp,
                    tc_relerr=relerr(tc_hat, tru["tc"]), cv_relerr=relerr(cv_hat, tru["cv"]),
                    # single-exp OU t_c -> its implied integrated tau, vs empirical tau, vs true tau
                    tauOU_relerr=relerr((1 + np.exp(-1 / max(tc_hat, 1e-3))) / (1 - np.exp(-1 / max(tc_hat, 1e-3))), tru["tau"]),
                    tauEMP_relerr=relerr(tau_emp, tru["tau"]) if np.isfinite(tau_emp) else float("nan")))
        def med(k): return float(np.nanmedian([r[k] for r in recs]))
        frac_pass = float(np.mean([r["passed"] for r in recs]))
        row = dict(cls=cls, n=len(recs), frac_pass=frac_pass,
                   med_p=med("p"), med_istd=med("istd"), med_ikurt=med("ikurt"),
                   med_tc_hat=med("tc_hat"), med_tc_true=med("tc_true"),
                   med_cv_hat=med("cv_hat"), med_cv_true=med("cv_true"),
                   med_tc_relerr=med("tc_relerr"), med_cv_relerr=med("cv_relerr"),
                   med_tauOU_relerr=med("tauOU_relerr"), med_tauEMP_relerr=med("tauEMP_relerr"))
        rows.append(row)
        log(f"{cls:10s} pass={frac_pass*100:5.1f}%  p={row['med_p']:.4f}  istd={row['med_istd']:.3f}  ikurt={row['med_ikurt']:+.2f}  "
            f"tc_hat={row['med_tc_hat']:5.1f}/tc_true={row['med_tc_true']:5.1f}(rel{row['med_tc_relerr']:.2f})  "
            f"cv rel{row['med_cv_relerr']:.2f}  | tau_relerr OU={row['med_tauOU_relerr']:.2f} EMP={row['med_tauEMP_relerr']:.2f}")

    meta = dict(exploratory=True, preregistered=False,
                utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
                note="E2(b) OU-model adequacy failure on non-OU classes; true-scene z; Ljung-Box K=20",
                classes=CLASSES, seeds=SEEDS, scenes=SCENES, runtime_s=time.time() - T0)
    json.dump(dict(meta=meta, rows=rows), open(os.path.join(OUT, "json", "e2_ou_misfit.json"), "w"), indent=2)
    log(f"saved e2_ou_misfit.json ({time.time()-T0:.1f}s)")


if __name__ == "__main__":
    main()
