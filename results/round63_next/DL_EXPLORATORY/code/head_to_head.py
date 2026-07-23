#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
EXPLORATORY / DEV-ONLY.  Head-to-head: SCGI-style learned correction (Arm L) vs our
model-based DLGI (Arm D) vs no-correction floor (Arm 0), on the 6 DEV bridge scenes
x OU drift cells, all arms sharing the IDENTICAL records.  No commits.

Arms (all reconstructions use the SAME arm_A4_gaincorr + frozen TV inversion; only the
per-exposure gain estimate differs):
  Arm 0        : gain = 1            (no correction floor)
  Arm L (mix)  : learned correction, trained on exp+OU  (in-family; SCGI's best case)
  Arm L (exp)  : learned correction, trained on exp ONLY (OUT-OF-FAMILY on OU -> item 4)
  Arm D        : joint_dual_ledger   (pilot-free OU/Kalman joint estimator, frozen)
  oracle       : true per-exposure gain (ceiling, reference only)

Medium products (t_c, CV):
  Arm D            -> its native medium estimate (the certified 'product')
  Arm L corrections-> (t_c,CV) fit to SCGI's implicit correction factors ('their trash')
  oracle-monitor   -> OU fit to the TRUE gain path (best-possible reference)
"""
import os, sys, json, time, platform
from datetime import datetime, timezone
import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import scgi_common as S
import dl_common as DL

OUT = "D:/GI_another/results/round63_next/DL_EXPLORATORY"
CKPT = os.path.join(OUT, "checkpoints")
T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
SCENES = DL.SCENES                              # the 6 DEV bridge scenes
TC_TEST = [16.0, 64.0]
CV_TEST = [0.15, 0.40]
SEEDS = [0, 1, 2]

def load_model(name):
    ck = torch.load(os.path.join(CKPT, f"{name}.pt"), map_location=DEVICE)
    net = S.UNet(c=ck["channels"]).to(DEVICE)
    net.load_state_dict(ck["state_dict"]); net.eval()
    return net, np.asarray(ck["baseline"], dtype=np.float64)

def relerr(hat, tru):
    return float(abs(hat - tru) / max(abs(tru), 1e-9))

def main():
    net_mix, base_mix = load_model("scgi_mix")
    net_exp, base_exp = load_model("scgi_exp")
    log(f"models loaded  device={DEVICE}")

    records = []
    for tc in TC_TEST:
        for cv in CV_TEST:
            sig = S.sigma_l_of_cv(cv)
            for sc in SCENES:
                x = DL.scene_x[sc]
                for seed in SEEDS:
                    # ONE record, paired schedule, shared by ALL arms
                    Y, a_time, slot = DL.simulate_record(sc, seed, tc, sig, schedule="paired")
                    ap_t, am_t = a_time[0::2], a_time[1::2]

                    # --- scene PSNR per arm ---
                    x0 = S.reconstruct_from_gain(Y, np.ones(S.N_EXP))          # Arm 0
                    g_mix = S.apply_correction(net_mix, base_mix, Y, DEVICE)   # Arm L mix
                    xL_mix = S.reconstruct_from_gain(Y, g_mix)
                    g_exp = S.apply_correction(net_exp, base_exp, Y, DEVICE)   # Arm L exp (OOD)
                    xL_exp = S.reconstruct_from_gain(Y, g_exp)
                    jd = DL.joint_dual_ledger(Y, slot=slot)                    # Arm D
                    xD = jd["x_hat"]
                    x_or = np.clip(DL.arm_A4_gaincorr(Y[0::2], Y[1::2], ap_t, am_t, DL.LAM_TV), 0.0, None)

                    # --- medium products ---
                    medD = jd["med"]                                          # DLGI product
                    medL_mix = S.medium_from_gain(g_mix)                      # SCGI trash (mix)
                    medL_exp = S.medium_from_gain(g_exp)                      # SCGI trash (exp/OOD)
                    med_or = DL.oracle_monitor_estimate(a_time)              # best-possible ref

                    rec = dict(
                        tc=tc, cv=cv, scene=sc, seed=seed,
                        psnr=dict(
                            arm0=DL.psnr(x0, x),
                            L_mix=DL.psnr(xL_mix, x),
                            L_exp=DL.psnr(xL_exp, x),
                            D=DL.psnr(xD, x),
                            oracle=DL.psnr(x_or, x)),
                        gain_corr=dict(                                       # recovered-gain fidelity
                            L_mix=DL.gain_path_corr(g_mix, a_time),
                            L_exp=DL.gain_path_corr(g_exp, a_time),
                            D=DL.gain_path_corr(jd["med"]["a_hat"], a_time)),
                        med=dict(
                            D=dict(tc=medD["tc_hat"], cv=medD["cv_hat"]),
                            L_mix=dict(tc=medL_mix["tc_hat"], cv=medL_mix["cv_hat"]),
                            L_exp=dict(tc=medL_exp["tc_hat"], cv=medL_exp["cv_hat"]),
                            oracle=dict(tc=med_or["tc_hat"], cv=med_or["cv_hat"])),
                    )
                    records.append(rec)
            log(f"  cell tc={tc:.0f} cv={cv:.2f} done ({len([r for r in records if r['tc']==tc and r['cv']==cv])} recs)")

    # -------- aggregate per cell --------
    def cell_recs(tc, cv): return [r for r in records if r["tc"] == tc and r["cv"] == cv]
    cells = []
    for tc in TC_TEST:
        for cv in CV_TEST:
            rs = cell_recs(tc, cv)
            def mean_psnr(a): return float(np.mean([r["psnr"][a] for r in rs]))
            def med_relerr(arm, key, tru):
                return float(np.median([relerr(r["med"][arm][key], tru) for r in rs]))
            cells.append(dict(
                tc=tc, cv=cv, n=len(rs),
                psnr={a: mean_psnr(a) for a in ["arm0", "L_mix", "L_exp", "D", "oracle"]},
                dpsnr_L_minus_D=float(np.mean([r["psnr"]["L_mix"] - r["psnr"]["D"] for r in rs])),
                dpsnr_D_minus_arm0=float(np.mean([r["psnr"]["D"] - r["psnr"]["arm0"] for r in rs])),
                gain_corr={a: float(np.mean([r["gain_corr"][a] for r in rs]))
                           for a in ["L_mix", "L_exp", "D"]},
                tc_relerr={a: med_relerr(a, "tc", tc) for a in ["D", "L_mix", "L_exp", "oracle"]},
                cv_relerr={a: med_relerr(a, "cv", cv) for a in ["D", "L_mix", "L_exp", "oracle"]},
            ))

    meta = dict(
        exploratory=True, preregistered=False, note="DEV-ONLY exploratory; no commits; "
        "faithful-in-spirit SCGI reimplementation (not code-identical).",
        utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
        python=platform.python_version(), torch=torch.__version__, device=DEVICE,
        PHI=DL.PHI, N_EXP=DL.N_EXP, schedule="paired",
        scenes=SCENES, tc_test=TC_TEST, cv_test=CV_TEST, seeds=SEEDS,
        n_records=len(records), runtime_s=time.time() - T0)
    out = dict(meta=meta, cells=cells, records=records)
    fn = os.path.join(OUT, "head_to_head_results.json")
    json.dump(out, open(fn, "w"), indent=2)
    log(f"saved {fn}  ({len(records)} records, {time.time()-T0:.1f}s)")

    # -------- console summary --------
    log("")
    log("=== SCENE PSNR (mean over 6 scenes x 3 seeds), dB ===")
    log(f"{'cell':16s} {'arm0':>7s} {'L_exp':>7s} {'L_mix':>7s} {'D':>7s} {'oracle':>7s}  {'L_mix-D':>8s}")
    for c in cells:
        p = c["psnr"]
        log(f"tc={c['tc']:>2.0f} cv={c['cv']:.2f}    {p['arm0']:7.2f} {p['L_exp']:7.2f} "
            f"{p['L_mix']:7.2f} {p['D']:7.2f} {p['oracle']:7.2f}  {c['dpsnr_L_minus_D']:+8.2f}")
    log("")
    log("=== MEDIUM t_c median rel-err  (D=product, L=trash) ===")
    for c in cells:
        t = c["tc_relerr"]
        log(f"tc={c['tc']:>2.0f} cv={c['cv']:.2f}   D={t['D']:.3f}  L_mix={t['L_mix']:.3f}  "
            f"L_exp={t['L_exp']:.3f}  oracle={t['oracle']:.3f}")
    log("=== MEDIUM CV median rel-err ===")
    for c in cells:
        t = c["cv_relerr"]
        log(f"tc={c['tc']:>2.0f} cv={c['cv']:.2f}   D={t['D']:.3f}  L_mix={t['L_mix']:.3f}  "
            f"L_exp={t['L_exp']:.3f}  oracle={t['oracle']:.3f}")

if __name__ == "__main__":
    main()
