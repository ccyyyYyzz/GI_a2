#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
EXPLORATORY / DEV-ONLY.  NOT preregistered.  No commits.  Local GPU.

E2 scene/gain channel, done HONESTLY in the full reconstruction loop (A2-scene log-
gain observation, NO true-scene leak -- the earlier true-scene isolation inflated
PSNR above the oracle because dividing Y by a gain formed from s_true re-injects the
scene).  Here the learned smoother is trained AND deployed on realistic A2-scene
features, matching how DLGI actually operates.

Arms per class {ou, regime, heavytail, quasiper}, 6 bridge scenes x seeds:
  OU-model DLGI  : D_arith  (single-exponential mom_autocov + Kalman RTS + E[a]=1)
  learned-prior  : same pipeline, Kalman smoother REPLACED by a trained TCN smoother
  oracle         : true per-exposure gain (ceiling)
Readouts: scene PSNR, gain-path log-RMSE, gain-path corr.  No-harm required on ou.
"""
import os, sys, json, time, platform
from datetime import datetime, timezone
import numpy as np
import torch
import torch.nn as nn

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "D:/GI_another/results/round63_next/DL_EXPLORATORY/code")
import e2_common as E
import dl_tool_common as T
import dl_common as DL
import scgi_common as S
from e2_learned_prior import SmootherTCN   # reuse arch

OUT = "D:/GI_another/results/round63_next/DL_TOOL_EXPLORATORY"
CKPT = os.path.join(OUT, "checkpoints")
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)
N_EXP = DL.N_EXP
CLASSES = ["ou", "regime", "heavytail", "quasiper"]


def a2_features(Y):
    """A2-scene (recon-loop) features: [zc, valid, prec] + gauge/order/Vs from the
    gain-robust init x_A2 (exactly what D_arith's smoother sees)."""
    Yp, Ym = Y[0::2], Y[1::2]
    x_a2 = np.clip(DL.arm_A2(Yp, Ym, DL.LAM_A2), 1e-9, None)
    p = T.gain_pieces(Y, x_a2)
    prec = p["valid"] * np.clip(1.0 / (p["R"] + 1e-9), 0, None)
    prec = prec / (prec[p["valid"] > 0].mean() + 1e-9)
    feat = np.stack([p["zc"] * p["valid"], p["valid"], np.log1p(prec)], axis=0).astype(np.float32)
    return feat, p, x_a2


def gen_data(n, seed0):
    rng = np.random.default_rng(5_500_000 + seed0)
    F = np.empty((n, 3, N_EXP), np.float32); L = np.empty((n, N_EXP), np.float32); V = np.empty((n, N_EXP), np.float32)
    for i in range(n):
        x = S.make_scene(rng); cls = CLASSES[i % len(CLASSES)]
        a_time, _ = E.sample_class(rng, cls)
        Y = E.simulate(x, a_time, rng)
        feat, p, _ = a2_features(Y)
        l = np.log(np.maximum(a_time, 1e-12)); lm = (l * p["valid"]).sum() / max(p["valid"].sum(), 1.0)
        F[i] = feat; L[i] = (l - lm).astype(np.float32); V[i] = p["valid"]
    return F, L, V


def train_smoother():
    cache = os.path.join(OUT, "json", "e2_recon_data.npz")
    if os.path.exists(cache):
        d = np.load(cache); F, L, V, Fv, Lv, Vv = d["F"], d["L"], d["V"], d["Fv"], d["Lv"], d["Vv"]
        log(f"loaded cached recon data {F.shape}")
    else:
        log("generating recon-loop training data (A2-scene features; arm_A2 per record)...")
        F, L, V = gen_data(3600, 0); log(f"train {F.shape} done")
        Fv, Lv, Vv = gen_data(900, 700000); log("val done")
        np.savez(cache, F=F, L=L, V=V, Fv=Fv, Lv=Lv, Vv=Vv)
    dev = DEVICE
    def t(a): return torch.from_numpy(a).to(dev)
    F, L, V, Fv, Lv, Vv = t(F), t(L), t(V), t(Fv), t(Lv), t(Vv)
    sm = SmootherTCN().to(dev); npar = sum(p.numel() for p in sm.parameters()); log(f"SmootherTCN params={npar}")
    opt = torch.optim.Adam(sm.parameters(), lr=2e-3)
    sched = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, factor=0.5, patience=5)
    def wmse(pr, tg, v): return ((pr - tg) ** 2 * v).sum() / v.sum().clamp_min(1.0)
    n = F.shape[0]; bs = 64; best = 1e9; best_state = None; bad = 0
    for ep in range(80):
        sm.train(); perm = torch.randperm(n, device=dev)
        for i in range(0, n, bs):
            idx = perm[i:i + bs]; pred = sm(F[idx])
            loss = wmse(pred, L[idx], V[idx]) + 0.01 * ((pred[:, 1:] - pred[:, :-1]) ** 2).mean()
            opt.zero_grad(); loss.backward(); opt.step()
        sm.eval()
        with torch.no_grad():
            vl = wmse(sm(Fv), Lv, Vv).item()
        sched.step(vl)
        if vl < best - 1e-7:
            best = vl; best_state = {k: v.detach().cpu().clone() for k, v in sm.state_dict().items()}; bad = 0
        else:
            bad += 1
        if ep % 10 == 0 or bad == 0:
            log(f"  [sm-recon] ep{ep:3d} val_wMSE={vl:.5e} (best {best:.5e})")
        if bad >= 12:
            log(f"  early stop ep{ep}"); break
    sm.load_state_dict(best_state); sm.eval()
    torch.save(dict(state_dict=sm.state_dict(), params=npar, best_val=best), os.path.join(CKPT, "e2_smoother_recon.pt"))
    return sm


def recon_learned(sm, Y):
    feat, p, _ = a2_features(Y)
    with torch.no_grad():
        l_hat = sm(torch.from_numpy(feat[None]).to(DEVICE)).cpu().numpy().reshape(N_EXP)
    a_hat = T.a_hat_from_l(l_hat, p["Vs"], p["mu_hat"], p["order"], jensen=False, renorm="arith")
    x = np.clip(DL.arm_A4_gaincorr(Y[0::2], Y[1::2], a_hat[0::2], a_hat[1::2], DL.LAM_TV), 0.0, None)
    return x, a_hat


def logrmse(a, at): return float(np.sqrt(np.mean((np.log(np.maximum(a, 1e-9)) - np.log(np.maximum(at, 1e-9))) ** 2)))


def evaluate(sm):
    SEEDS = list(range(6)); rows = []
    for cls in CLASSES:
        recs = []
        for sc in DL.SCENES:
            x = DL.scene_x[sc]
            for seed in SEEDS:
                r2 = np.random.default_rng(hash((cls, sc, seed, "e2recon")) % (2**32))
                a_time, _ = E.sample_class(r2, cls); Y = E.simulate(x, a_time, r2)
                jd = T.joint_dual_ledger_cfg(Y, slot=None, renorm="arith")          # OU-model D_arith
                x_ou = jd["x_hat"]; a_ou = jd["med"]["a_hat"]
                x_lp, a_lp = recon_learned(sm, Y)                                    # learned prior
                x_or = np.clip(DL.arm_A4_gaincorr(Y[0::2], Y[1::2], a_time[0::2], a_time[1::2], DL.LAM_TV), 0.0, None)
                recs.append(dict(psnr_ou=DL.psnr(x_ou, x), psnr_lp=DL.psnr(x_lp, x), psnr_or=DL.psnr(x_or, x),
                                 rmse_ou=logrmse(a_ou, a_time), rmse_lp=logrmse(a_lp, a_time),
                                 corr_ou=DL.gain_path_corr(a_ou, a_time), corr_lp=DL.gain_path_corr(a_lp, a_time)))
        def m(k): return float(np.mean([r[k] for r in recs]))
        row = dict(cls=cls, n=len(recs), psnr_ou=m("psnr_ou"), psnr_lp=m("psnr_lp"), psnr_or=m("psnr_or"),
                   dpsnr=m("psnr_lp") - m("psnr_ou"), rmse_ou=m("rmse_ou"), rmse_lp=m("rmse_lp"),
                   corr_ou=m("corr_ou"), corr_lp=m("corr_lp"))
        rows.append(row)
        log(f"{cls:10s} PSNR ou={row['psnr_ou']:6.2f} lp={row['psnr_lp']:6.2f} or={row['psnr_or']:6.2f} "
            f"(lp-ou {row['dpsnr']:+.2f})  RMSE {row['rmse_ou']:.3f}->{row['rmse_lp']:.3f}  corr {row['corr_ou']:.4f}->{row['corr_lp']:.4f}")
    return rows


def main():
    ck = os.path.join(CKPT, "e2_smoother_recon.pt")
    if os.path.exists(ck):
        sm = SmootherTCN().to(DEVICE); sm.load_state_dict(torch.load(ck, map_location=DEVICE)["state_dict"]); sm.eval(); log("loaded recon smoother")
    else:
        sm = train_smoother()
    rows = evaluate(sm)
    meta = dict(exploratory=True, preregistered=False, utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
                python=platform.python_version(), torch=torch.__version__, device=DEVICE,
                note="E2 scene/gain channel in the honest recon loop (A2-scene z, no leak)")
    json.dump(dict(meta=meta, rows=rows), open(os.path.join(OUT, "json", "e2_recon_psnr.json"), "w"), indent=2)
    log("saved e2_recon_psnr.json")


if __name__ == "__main__":
    main()
