#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
EXPLORATORY / DEV-ONLY.  NOT preregistered.  No commits.

E3: certified-ceiling demo for ONE cell (the fast-strong corner t_c=16, CV=0.40 --
the corner where the model-based DLGI was weakest and the learned tools help most).
Every method is placed on a (performance vs certified-ceiling) axis:
  * y = scene PSNR (performance)
  * x = gain-path fidelity (log-gain RMSE; ->0 is the ceiling)
  * ceiling = the ORACLE known-gain arm (star); the CRB sd floor on the MEDIUM
    parameters (t_c, CV) is the referee that certifies which methods actually
    MEASURE the medium vs merely correct it.
The picture: DL and model-based tools RACE the oracle ceiling; the oracle/CRB REFEREE.

All methods use the SAME full-recon protocol (E1's) on the 6 bridge scenes x 3 seeds.
"""
import os, sys, json, time
from datetime import datetime, timezone
import numpy as np
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "D:/GI_another/results/round63_next/DL_EXPLORATORY/code")
import dl_tool_common as T
import dl_common as DL
import scgi_common as S
from e1c_residual_corrector import Corrector1D, N_EXP

OUT = "D:/GI_another/results/round63_next/DL_TOOL_EXPLORATORY"
CKPT = os.path.join(OUT, "checkpoints")
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

TC, CV = 16.0, 0.40
SEEDS = [0, 1, 2]


def load_scgi():
    ck = torch.load(os.path.join(CKPT.replace("DL_TOOL_EXPLORATORY", "DL_EXPLORATORY"), "scgi_mix.pt"), map_location=DEVICE)
    net = S.UNet(c=ck["channels"]).to(DEVICE); net.load_state_dict(ck["state_dict"]); net.eval()
    return net, np.asarray(ck["baseline"], dtype=np.float64)


def load_corr():
    net = Corrector1D().to(DEVICE)
    net.load_state_dict(torch.load(os.path.join(CKPT, "e1c_corrector.pt"), map_location=DEVICE)["state_dict"]); net.eval()
    return net


def corr_deploy(net, Y):
    Yp = Y[0::2]; Ym = Y[1::2]
    x_a2 = np.clip(DL.arm_A2(Yp, Ym, DL.LAM_A2), 1e-9, None)
    p = T.gain_pieces(Y, x_a2)
    prec = p["valid"] * np.clip(1.0 / (p["R"] + 1e-9), 0, None); prec = prec / (prec[p["valid"] > 0].mean() + 1e-9)
    feat = np.stack([p["ms"], p["zc"] * p["valid"], p["valid"], np.log1p(prec)], axis=0).astype(np.float32)
    with torch.no_grad():
        dl = net(torch.from_numpy(feat[None]).to(DEVICE)).cpu().numpy().reshape(N_EXP)
    a_hat = T.a_hat_from_l(p["ms"] + dl, p["Vs"], p["mu_hat"], p["order"], jensen=True, renorm="arith")
    return np.clip(DL.arm_A4_gaincorr(Yp, Ym, a_hat[0::2], a_hat[1::2], DL.LAM_TV), 0.0, None), a_hat


def logrmse(a_hat, a_time):
    return float(np.sqrt(np.mean((np.log(np.maximum(a_hat, 1e-9)) - np.log(np.maximum(a_time, 1e-9))) ** 2)))


def main():
    net_scgi, base_scgi = load_scgi()
    net_corr = load_corr()
    sig = T.sigma_l_of_cv(CV)
    methods = ["floor", "D0", "D_arith", "D_arith+corr", "SCGI", "oracle"]
    acc = {m: dict(psnr=[], rmse=[], corr=[], tc=[], cv=[]) for m in methods}
    for sc in DL.SCENES:
        x = DL.scene_x[sc]
        for seed in SEEDS:
            Y, a_time, slot = DL.simulate_record(sc, seed, TC, sig, schedule="paired")
            Yp, Ym = Y[0::2], Y[1::2]
            ap_t, am_t = a_time[0::2], a_time[1::2]
            # floor
            xf = S.reconstruct_from_gain(Y, np.ones(N_EXP)); acc["floor"]["psnr"].append(DL.psnr(xf, x))
            acc["floor"]["rmse"].append(logrmse(np.ones(N_EXP), a_time)); acc["floor"]["corr"].append(0.0)
            acc["floor"]["tc"].append(np.nan); acc["floor"]["cv"].append(np.nan)
            # D0 frozen
            jd0 = T.joint_dual_ledger_cfg(Y, slot=slot, renorm="none")
            acc["D0"]["psnr"].append(DL.psnr(jd0["x_hat"], x)); acc["D0"]["rmse"].append(logrmse(jd0["med"]["a_hat"], a_time))
            acc["D0"]["corr"].append(DL.gain_path_corr(jd0["med"]["a_hat"], a_time))
            acc["D0"]["tc"].append(jd0["med"]["tc_hat"]); acc["D0"]["cv"].append(jd0["med"]["cv_hat"])
            # D_arith
            jda = T.joint_dual_ledger_cfg(Y, slot=slot, renorm="arith")
            acc["D_arith"]["psnr"].append(DL.psnr(jda["x_hat"], x)); acc["D_arith"]["rmse"].append(logrmse(jda["med"]["a_hat"], a_time))
            acc["D_arith"]["corr"].append(DL.gain_path_corr(jda["med"]["a_hat"], a_time))
            acc["D_arith"]["tc"].append(jda["med"]["tc_hat"]); acc["D_arith"]["cv"].append(jda["med"]["cv_hat"])
            # D_arith + learned corrector
            xC, aC = corr_deploy(net_corr, Y)
            acc["D_arith+corr"]["psnr"].append(DL.psnr(xC, x)); acc["D_arith+corr"]["rmse"].append(logrmse(aC, a_time))
            acc["D_arith+corr"]["corr"].append(DL.gain_path_corr(aC, a_time))
            mc = S.medium_from_gain(aC); acc["D_arith+corr"]["tc"].append(mc["tc_hat"]); acc["D_arith+corr"]["cv"].append(mc["cv_hat"])
            # SCGI learned corrector
            g = S.apply_correction(net_scgi, base_scgi, Y, DEVICE)
            xS = S.reconstruct_from_gain(Y, g); acc["SCGI"]["psnr"].append(DL.psnr(xS, x)); acc["SCGI"]["rmse"].append(logrmse(g, a_time))
            acc["SCGI"]["corr"].append(DL.gain_path_corr(g, a_time))
            ms = S.medium_from_gain(g); acc["SCGI"]["tc"].append(ms["tc_hat"]); acc["SCGI"]["cv"].append(ms["cv_hat"])
            # oracle
            xo = np.clip(DL.arm_A4_gaincorr(Yp, Ym, ap_t, am_t, DL.LAM_TV), 0.0, None)
            acc["oracle"]["psnr"].append(DL.psnr(xo, x)); acc["oracle"]["rmse"].append(0.0); acc["oracle"]["corr"].append(1.0)
            acc["oracle"]["tc"].append(np.nan); acc["oracle"]["cv"].append(np.nan)

    crb = DL.oracle_floor(TC)
    summ = {m: dict(psnr=float(np.mean(v["psnr"])), rmse=float(np.mean(v["rmse"])),
                    corr=float(np.mean(v["corr"])),
                    tc=float(np.nanmean(v["tc"])) if np.any(np.isfinite(v["tc"])) else None,
                    cv=float(np.nanmean(v["cv"])) if np.any(np.isfinite(v["cv"])) else None)
            for m, v in acc.items()}
    for m in methods:
        s = summ[m]; log(f"{m:14s} PSNR={s['psnr']:6.2f}  gain-logRMSE={s['rmse']:.4f}  corr={s['corr']:.4f}  "
                        f"tc={s['tc'] if s['tc'] else float('nan'):.2f}  cv={s['cv'] if s['cv'] else float('nan'):.3f}")
    log(f"CRB sd floor (t_c={TC:.0f}): tc_rel={crb['sd_tc_rel']:.3f}  cv_rel={crb['sd_cv_rel']:.3f}")
    json.dump(dict(meta=dict(exploratory=True, tc=TC, cv=CV, seeds=SEEDS, crb=crb,
                   utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ")), summary=summ),
              open(os.path.join(OUT, "json", "e3_ceiling.json"), "w"), indent=2)

    # ---------------- figure (single panel: performance vs certified ceiling) ------
    fig, axA = plt.subplots(1, 1, figsize=(8.6, 6.0))
    order = ["floor", "D0", "D_arith", "D_arith+corr", "SCGI", "oracle"]
    colors = {"floor": "#8c8c8c", "D0": "#d1495b", "D_arith": "#edae49",
              "D_arith+corr": "#66a182", "SCGI": "#2e86ab", "oracle": "#111111"}
    labels = {"floor": "no-corr floor", "D0": "DLGI frozen", "D_arith": "DLGI + E[a]=1 fix (physics)",
              "D_arith+corr": "DLGI fix + learned gain corrector", "SCGI": "learned corrector (SCGI, in-family)",
              "oracle": "oracle known-gain (ceiling)"}
    # Panel A: performance vs gain-fidelity ceiling
    orc = summ["oracle"]["psnr"]
    axA.axhline(orc, ls="--", c="k", lw=1.2, alpha=0.7)
    # label offsets (x-mult, y-add, ha) to avoid overlap
    off = {"floor": (1.0, 0.55, "left"), "D0": (1.0, 0.55, "center"),
           "D_arith": (1.0, -0.85, "center"), "D_arith+corr": (0.72, 0.55, "center"),
           "SCGI": (1.35, -0.9, "center"), "oracle": (1.9, 0.0, "left")}
    for m in order:
        s = summ[m]; xx = max(s["rmse"], 1.5e-3)
        axA.scatter(xx, s["psnr"], s=160, c=colors[m], edgecolor="k", zorder=3, marker=("*" if m == "oracle" else "o"))
        fx, fy, ha = off[m]
        axA.annotate(labels[m], (xx, s["psnr"]), (xx * fx, s["psnr"] + fy), fontsize=8.4, ha=ha)
    axA.set_xscale("log"); axA.invert_xaxis()
    axA.set_ylim(5.2, 19.8)
    axA.set_xlabel("gain-path log-RMSE   (lower = closer to certified gain ceiling  →)")
    axA.set_ylabel("scene PSNR  (dB)")
    axA.set_title(f"E3  DL races the ceiling, the oracle/CRB referees   (t_c={TC:.0f}, CV={CV:.2f})\nEXPLORATORY, DEV-only", fontsize=10.5)
    axA.grid(alpha=0.25)
    # the referee: annotate the medium-channel CRB certificate in prose on the ceiling
    axA.annotate("", xy=(summ["D_arith"]["rmse"], summ["D_arith"]["psnr"] + 0.2),
                 xytext=(summ["D0"]["rmse"], summ["D0"]["psnr"] + 0.2),
                 arrowprops=dict(arrowstyle="->", color="#555", lw=1.4))
    axA.text(0.13, 9.9, "physics fix\n(+2.2 dB, no learning)", fontsize=8, color="#555", ha="center", style="italic")
    axA.text(0.986, 0.02,
             "referee: the medium (t_c, CV) is read from the SAME record with a\n"
             f"Cramér–Rao certificate (sd floor t_c {crb['sd_tc_rel']*100:.0f}%, CV {crb['sd_cv_rel']*100:.0f}%);\n"
             "the learned correctors' medium point carries no such certificate.",
             transform=axA.transAxes, fontsize=7.6, ha="right", va="bottom",
             bbox=dict(boxstyle="round", fc="#f4f4f4", ec="#bbb", alpha=0.9))
    fig.tight_layout()
    fp = os.path.join(OUT, "figs", "e3_certified_ceiling.png")
    fig.savefig(fp, dpi=140); log(f"saved {fp}")


if __name__ == "__main__":
    main()
