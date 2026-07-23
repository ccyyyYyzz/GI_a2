#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
EXPLORATORY / DEV-ONLY.  Off-grid sensitivity probe.  The head-to-head test cells
(t_c in {16,64}, CV in {0.15,0.40}) sit EXACTLY on the SCGI-MIX OU training grid
(t_c in {8,16,32,64}, CV in {5,15,40}%), giving the learned corrector a home-field
advantage for the medium POINT estimate.  Here we test cells OFF the training grid
to separate genuine within-OU generalisation from grid-node recall.  DLGI (Arm D)
is model-based / training-free, so it should be unaffected by grid position.

  interior     : (t_c=24, CV=0.25)  -- between grid nodes (interpolation)
  extrap_tc    : (t_c=128, CV=0.15) -- t_c beyond training max (64)
  extrap_cv    : (t_c=16,  CV=0.55) -- CV beyond training max (0.40)
No commits.
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

CELLS = [dict(name="interior", tc=24.0, cv=0.25),
         dict(name="extrap_tc", tc=128.0, cv=0.15),
         dict(name="extrap_cv", tc=16.0, cv=0.55)]
SEEDS = [0, 1, 2]

def load_model(name):
    ck = torch.load(os.path.join(CKPT, f"{name}.pt"), map_location=DEVICE)
    net = S.UNet(c=ck["channels"]).to(DEVICE)
    net.load_state_dict(ck["state_dict"]); net.eval()
    return net, np.asarray(ck["baseline"], dtype=np.float64)

def relerr(hat, tru): return float(abs(hat - tru) / max(abs(tru), 1e-9))

def main():
    net_mix, base_mix = load_model("scgi_mix")
    net_exp, base_exp = load_model("scgi_exp")
    log(f"models loaded  device={DEVICE}")
    out_cells = []
    for cell in CELLS:
        tc, cv = cell["tc"], cell["cv"]; sig = S.sigma_l_of_cv(cv)
        recs = []
        for sc in DL.SCENES:
            x = DL.scene_x[sc]
            for seed in SEEDS:
                Y, a_time, slot = DL.simulate_record(sc, seed, tc, sig, schedule="paired")
                x0 = S.reconstruct_from_gain(Y, np.ones(S.N_EXP))
                g_mix = S.apply_correction(net_mix, base_mix, Y, DEVICE)
                g_exp = S.apply_correction(net_exp, base_exp, Y, DEVICE)
                xL_mix = S.reconstruct_from_gain(Y, g_mix)
                jd = DL.joint_dual_ledger(Y, slot=slot)
                x_or = np.clip(DL.arm_A4_gaincorr(Y[0::2], Y[1::2], a_time[0::2], a_time[1::2], DL.LAM_TV), 0.0, None)
                medL_mix = S.medium_from_gain(g_mix)
                medL_exp = S.medium_from_gain(g_exp)
                med_or = DL.oracle_monitor_estimate(a_time)
                recs.append(dict(scene=sc, seed=seed,
                    psnr=dict(arm0=DL.psnr(x0, x), L_mix=DL.psnr(xL_mix, x),
                              D=DL.psnr(jd["x_hat"], x), oracle=DL.psnr(x_or, x)),
                    gain_corr=dict(L_mix=DL.gain_path_corr(g_mix, a_time),
                                   D=DL.gain_path_corr(jd["med"]["a_hat"], a_time)),
                    med=dict(D=dict(tc=jd["med"]["tc_hat"], cv=jd["med"]["cv_hat"]),
                             L_mix=dict(tc=medL_mix["tc_hat"], cv=medL_mix["cv_hat"]),
                             L_exp=dict(tc=medL_exp["tc_hat"], cv=medL_exp["cv_hat"]),
                             oracle=dict(tc=med_or["tc_hat"], cv=med_or["cv_hat"]))))
        def mp(a): return float(np.mean([r["psnr"][a] for r in recs]))
        def mre(arm, key, tru): return float(np.median([relerr(r["med"][arm][key], tru) for r in recs]))
        c = dict(name=cell["name"], tc=tc, cv=cv, n=len(recs),
                 psnr={a: mp(a) for a in ["arm0", "L_mix", "D", "oracle"]},
                 gain_corr={a: float(np.mean([r["gain_corr"][a] for r in recs])) for a in ["L_mix", "D"]},
                 tc_relerr={a: mre(a, "tc", tc) for a in ["D", "L_mix", "L_exp", "oracle"]},
                 cv_relerr={a: mre(a, "cv", cv) for a in ["D", "L_mix", "L_exp", "oracle"]},
                 records=recs)
        out_cells.append(c)
        log(f"  {cell['name']:10s} tc={tc:.0f} cv={cv:.2f}: PSNR L_mix={c['psnr']['L_mix']:.2f} D={c['psnr']['D']:.2f} "
            f"| tc_relerr D={c['tc_relerr']['D']:.3f} L_mix={c['tc_relerr']['L_mix']:.3f} "
            f"| cv_relerr D={c['cv_relerr']['D']:.3f} L_mix={c['cv_relerr']['L_mix']:.3f}")
    meta = dict(exploratory=True, preregistered=False, utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
                python=platform.python_version(), torch=torch.__version__, device=DEVICE,
                note="off-grid sensitivity; SCGI-MIX OU training grid tc in {8,16,32,64}, CV in {5,15,40}%",
                seeds=SEEDS, runtime_s=time.time() - T0)
    json.dump(dict(meta=meta, cells=out_cells), open(os.path.join(OUT, "off_grid_results.json"), "w"), indent=2)
    log(f"saved off_grid_results.json ({time.time()-T0:.1f}s)")

if __name__ == "__main__":
    main()
