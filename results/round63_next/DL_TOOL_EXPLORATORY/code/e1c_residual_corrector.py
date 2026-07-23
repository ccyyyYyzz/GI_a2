#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
EXPLORATORY / DEV-ONLY.  NOT preregistered.  No commits.  Local GPU.

E1 step (c): the analytic E[a]=1 fix (D_arith) closes ~2.2 dB of the fast-strong
corner deficit but leaves a residual (~1.3 dB to SCGI, ~7.4 dB to the oracle) that
is SHAPE (path-tracking) not scale.  Question: how much MORE does a tiny learned
layer buy OVER the analytic fix, with the oracle known-gain arm as the ceiling?

Design (honest "DL races the certified gain ceiling"):
  * The learned layer is a small 1-D temporal CNN that REFINES the DLGI-smoothed
    log-gain path (the physics-adjacent quantity), not the scene -- it earns exactly
    the part the fixed scalar-OU smoother cannot write down (optimal local smoothing
    at fast-strong drift, where the OU prior over-smooths sharp gain excursions).
  * Input (time order, 2048): [ms (OU-smoothed demeaned log-gain), zc (raw demeaned
    log-gain observation), validity, counts-precision].  Output: a correction dl to
    ms.  Corrected gain rebuilt with the SAME Jensen + E[a]=1 renorm -> scene recon.
  * Trained on PROCEDURAL scenes + OU gain, disjoint seeds; CERTIFIED on held-out
    procedural records (gain-path RMSE / corr) AND evaluated on the 6 held-out BRIDGE
    scenes at the 4 test cells (the honest out-of-distribution scene test).
  * Ceiling = oracle (true gain).  No-harm check: must not degrade the low-CV cells.
"""
import os, sys, json, time, platform
from datetime import datetime, timezone
import numpy as np
import torch
import torch.nn as nn

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "D:/GI_another/results/round63_next/DL_EXPLORATORY/code")
import dl_tool_common as T
import dl_common as DL
import scgi_common as S

OUT = "D:/GI_another/results/round63_next/DL_TOOL_EXPLORATORY"
CKPT = os.path.join(OUT, "checkpoints")
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

N_EXP = DL.N_EXP
IMG_H, IMG_W = 32, 64   # 2048 = 32x64, row-major = time order (paired schedule)

# ---------------- procedural record generation (disjoint seeds) ----------------
def make_features(Y, a_time):
    """One record -> (feat [4,2048], target [2048], valid [2048], mu_hat, order, Vs)."""
    Yp = Y[0::2]; Ym = Y[1::2]
    x_a2 = np.clip(DL.arm_A2(Yp, Ym, DL.LAM_A2), 1e-9, None)      # realistic init
    p = T.gain_pieces(Y, x_a2)                                    # slot=paired -> order=identity
    ms, zc, R, valid, mu_hat = p["ms"], p["zc"], p["R"], p["valid"], p["mu_hat"]
    prec = valid * np.clip(1.0 / (R + 1e-9), 0, None)
    prec = prec / (prec[valid > 0].mean() + 1e-9)                 # ~1 on valid
    feat = np.stack([ms, zc * valid, valid, np.log1p(prec)], axis=0).astype(np.float32)
    l_true = np.log(np.maximum(a_time, 1e-12))
    target = (l_true - mu_hat).astype(np.float32)                 # corrected ms should hit this
    return feat, target, valid.astype(np.float32), mu_hat, p["order"], p["Vs"]


def gen_dataset(n, seed0):
    rng = np.random.default_rng(3_000_000 + seed0)
    F = np.empty((n, 4, N_EXP), np.float32)
    Tg = np.empty((n, N_EXP), np.float32)
    V = np.empty((n, N_EXP), np.float32)
    for i in range(n):
        x = S.make_scene(rng)
        tc = float(np.exp(rng.uniform(np.log(6.0), np.log(80.0))))   # broad OU, fast-strong emphasis
        cv = float(rng.uniform(0.05, 0.55))
        sig = T.sigma_l_of_cv(cv)
        a_time = T.ou_path(rng, sig, tc, N_EXP)
        s = T.s_of_x(x)
        Y = rng.poisson(np.maximum(a_time * T.PHI * s, 0.0)).astype(np.float64)
        f, tg, v, *_ = make_features(Y, a_time)
        F[i] = f; Tg[i] = tg; V[i] = v
    return F, Tg, V


# ---------------- tiny 1-D temporal CNN corrector ----------------
class Corrector1D(nn.Module):
    """Dilated 1-D CNN on the length-2048 gain path.  <~120k params.  Predicts a
    correction dl added to ms; output re-gauged downstream by the E[a]=1 renorm."""
    def __init__(self, cin=4, c=48):
        super().__init__()
        def blk(ci, co, d):
            return nn.Sequential(nn.Conv1d(ci, co, 5, padding=2 * d, dilation=d),
                                 nn.BatchNorm1d(co), nn.ReLU(inplace=True))
        self.net = nn.Sequential(blk(cin, c, 1), blk(c, c, 2), blk(c, c, 4),
                                 blk(c, c, 8), blk(c, c, 1))
        self.out = nn.Conv1d(c, 1, 1)
    def forward(self, x):           # x: (B,4,L)
        return self.out(self.net(x)).squeeze(1)   # (B,L) correction dl


def train():
    log("generating training data (procedural scenes, disjoint seeds)...")
    ntr, nva = 3200, 800
    cache = os.path.join(OUT, "json", "e1c_data.npz")
    if os.path.exists(cache):
        d = np.load(cache)
        Ftr, Ttr, Vtr, Fva, Tva, Vva = d["Ftr"], d["Ttr"], d["Vtr"], d["Fva"], d["Tva"], d["Vva"]
        log(f"loaded cached data {Ftr.shape}")
    else:
        Ftr, Ttr, Vtr = gen_dataset(ntr, 0)
        log(f"train set {Ftr.shape} done")
        Fva, Tva, Vva = gen_dataset(nva, 900000)
        log(f"val set {Fva.shape} done")
        np.savez(cache, Ftr=Ftr, Ttr=Ttr, Vtr=Vtr, Fva=Fva, Tva=Tva, Vva=Vva)
    # tensors
    def t(a): return torch.from_numpy(a)
    Ftr, Ttr, Vtr = t(Ftr).to(DEVICE), t(Ttr).to(DEVICE), t(Vtr).to(DEVICE)
    Fva, Tva, Vva = t(Fva).to(DEVICE), t(Tva).to(DEVICE), t(Vva).to(DEVICE)
    net = Corrector1D().to(DEVICE)
    npar = sum(p.numel() for p in net.parameters())
    log(f"Corrector1D params={npar}")
    opt = torch.optim.Adam(net.parameters(), lr=2e-3)
    sched = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, factor=0.5, patience=6)
    ms_tr = Ftr[:, 0, :]; ms_va = Fva[:, 0, :]   # current estimate = channel 0
    n = Ftr.shape[0]; bs = 64; best = 1e9; best_state = None; patience = 16; bad = 0

    def wmse(pred_l, tgt, v):
        e = (pred_l - tgt) ** 2 * v
        return e.sum() / v.sum().clamp_min(1.0)

    # baseline (uncorrected) val error, for reference
    with torch.no_grad():
        base_val = wmse(ms_va, Tva, Vva).item()
    log(f"val MSE uncorrected (D_arith ms vs truth): {base_val:.5e}")
    for ep in range(120):
        net.train(); perm = torch.randperm(n, device=DEVICE)
        for i in range(0, n, bs):
            idx = perm[i:i + bs]
            dl = net(Ftr[idx])
            pred = ms_tr[idx] + dl
            # smoothness prior on the correction (drift is smooth)
            sm = ((dl[:, 1:] - dl[:, :-1]) ** 2).mean()
            loss = wmse(pred, Ttr[idx], Vtr[idx]) + 0.02 * sm
            opt.zero_grad(); loss.backward(); opt.step()
        net.eval()
        with torch.no_grad():
            vpred = ms_va + net(Fva)
            vloss = wmse(vpred, Tva, Vva).item()
        sched.step(vloss)
        if vloss < best - 1e-7:
            best = vloss; best_state = {k: v.detach().cpu().clone() for k, v in net.state_dict().items()}; bad = 0
        else:
            bad += 1
        if ep % 10 == 0 or bad == 0:
            log(f"  ep{ep:3d} val_MSE={vloss:.5e}  (best {best:.5e})")
        if bad >= patience:
            log(f"early stop ep{ep}"); break
    net.load_state_dict(best_state)
    torch.save(dict(state_dict=net.state_dict(), params=npar, base_val=base_val, best_val=best,
                    note="EXPLORATORY gain-path residual corrector; D_arith refinement"),
               os.path.join(CKPT, "e1c_corrector.pt"))
    log(f"saved corrector.  val MSE {base_val:.4e} (uncorrected) -> {best:.4e} (corrected), "
        f"variance-reduction {100*(1-best/base_val):.1f}%")
    return net


# ---------------- deploy: DLGI D_arith + learned gain refinement ----------------
def deploy_recon(net, Y, a_time):
    Yp = Y[0::2]; Ym = Y[1::2]
    x_a2 = np.clip(DL.arm_A2(Yp, Ym, DL.LAM_A2), 1e-9, None)
    p = T.gain_pieces(Y, x_a2)
    ms, zc, R, valid, mu_hat, Vs, order = (p["ms"], p["zc"], p["R"], p["valid"],
                                           p["mu_hat"], p["Vs"], p["order"])
    prec = valid * np.clip(1.0 / (R + 1e-9), 0, None); prec = prec / (prec[valid > 0].mean() + 1e-9)
    feat = np.stack([ms, zc * valid, valid, np.log1p(prec)], axis=0).astype(np.float32)
    with torch.no_grad():
        dl = net(torch.from_numpy(feat[None]).to(DEVICE)).cpu().numpy().reshape(N_EXP)
    ms_corr = ms + dl
    a_hat = T.a_hat_from_l(ms_corr, Vs, mu_hat, order, jensen=True, renorm="arith")
    x = np.clip(DL.arm_A4_gaincorr(Yp, Ym, a_hat[0::2], a_hat[1::2], DL.LAM_TV), 0.0, None)
    return x, a_hat


def evaluate(net):
    TC_TEST = [16.0, 64.0]; CV_TEST = [0.15, 0.40]; SEEDS = [0, 1, 2]
    cells = []
    for tc in TC_TEST:
        for cv in CV_TEST:
            sig = T.sigma_l_of_cv(cv); recs = []
            for sc in DL.SCENES:
                x = DL.scene_x[sc]
                for seed in SEEDS:
                    Y, a_time, slot = DL.simulate_record(sc, seed, tc, sig, schedule="paired")
                    jd = T.joint_dual_ledger_cfg(Y, slot=slot, renorm="arith")     # D_arith
                    xDa = jd["x_hat"]; a_Da = jd["med"]["a_hat"]
                    xC, a_C = deploy_recon(net, Y, a_time)                          # + corrector
                    x_or = np.clip(DL.arm_A4_gaincorr(Y[0::2], Y[1::2], a_time[0::2], a_time[1::2], DL.LAM_TV), 0.0, None)
                    recs.append(dict(
                        psnr_Darith=DL.psnr(xDa, x), psnr_corr=DL.psnr(xC, x), psnr_oracle=DL.psnr(x_or, x),
                        corr_Darith=DL.gain_path_corr(a_Da, a_time), corr_C=DL.gain_path_corr(a_C, a_time),
                        rmse_Darith=float(np.sqrt(np.mean((np.log(np.maximum(a_Da,1e-9))-np.log(np.maximum(a_time,1e-9)))**2))),
                        rmse_C=float(np.sqrt(np.mean((np.log(np.maximum(a_C,1e-9))-np.log(np.maximum(a_time,1e-9)))**2)))))
            def m(k): return float(np.mean([r[k] for r in recs]))
            c = dict(tc=tc, cv=cv, psnr_Darith=m("psnr_Darith"), psnr_corr=m("psnr_corr"),
                     psnr_oracle=m("psnr_oracle"), gain_dB=m("psnr_corr")-m("psnr_Darith"),
                     corr_Darith=m("corr_Darith"), corr_C=m("corr_C"),
                     rmse_Darith=m("rmse_Darith"), rmse_C=m("rmse_C"))
            cells.append(c)
            log(f"tc={tc:>2.0f} cv={cv:.2f}  D_arith={c['psnr_Darith']:6.2f} +corr={c['psnr_corr']:6.2f} "
                f"(+{c['gain_dB']:+.2f}dB)  oracle={c['psnr_oracle']:6.2f}  | gain-RMSE {c['rmse_Darith']:.3f}->{c['rmse_C']:.3f}")
    return cells


def main():
    ck = os.path.join(CKPT, "e1c_corrector.pt")
    net = Corrector1D().to(DEVICE)
    if os.path.exists(ck):
        net.load_state_dict(torch.load(ck, map_location=DEVICE)["state_dict"]); log("loaded corrector ckpt")
    else:
        net = train()
    cells = evaluate(net)
    meta = dict(exploratory=True, preregistered=False,
                utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
                python=platform.python_version(), torch=torch.__version__, device=DEVICE,
                note="E1(c) learned gain-path residual corrector over the analytic D_arith fix")
    json.dump(dict(meta=meta, cells=cells), open(os.path.join(OUT, "json", "e1c_corrector.json"), "w"), indent=2)
    log("saved e1c_corrector.json")


if __name__ == "__main__":
    main()
