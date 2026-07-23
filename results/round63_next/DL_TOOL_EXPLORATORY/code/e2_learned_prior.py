#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
EXPLORATORY / DEV-ONLY.  NOT preregistered.  No commits.  Local GPU.

E2 steps (c)+(d)+(e): a small trained sequence model as the LEARNED TEMPORAL PRIOR,
replacing the frozen scalar-OU assumption, tested as a domain extension beyond the
paper's validated OU domain.

Two <1M-param models, trained on the MIXED non-OU family (ou+regime+heavytail+
quasiper), procedural scenes, disjoint seeds:
  * SMOOTHER (acausal dilated TCN): [zc, valid, prec] -> smoothed demeaned log-gain
    path l_hat.  Plugged into an estimator copy for the scene/path channel.
  * PREDICTOR (causal dilated TCN, Gaussian head): the LEARNED TRANSITION DENSITY --
    one-step-ahead (mu_t, var_t) for z_t given the strict past.  Its standardized
    residuals are the learned analog of the OU-Kalman innovations, so the SAME
    whiteness test certifies model adequacy.

Isolation: E2 forms the log-gain observation z from the TRUE scene (identically for
every arm), so the comparison is purely about the TEMPORAL GAIN MODEL, not scene
reconstruction transfer (E1 already stress-tested the full recon loop).

Readouts per class {ou, regime, heavytail, quasiper} x arm {OU-model DLGI,
learned-prior DLGI, oracle}: scene PSNR, gain-path corr, log-gain RMSE, medium
(t_c / CV / tau) rel-err, innovation-whiteness pass rate.  No-harm required on OU.
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

OUT = "D:/GI_another/results/round63_next/DL_TOOL_EXPLORATORY"
CKPT = os.path.join(OUT, "checkpoints")
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)
N_EXP = DL.N_EXP
CLASSES = ["ou", "regime", "heavytail", "quasiper"]


# ------------------------------- feature extraction -------------------------------
def features_true_scene(Y, x_true):
    zc, R, valid, mu_hat = E.z_from_true_scene(Y, x_true)
    prec = valid * np.clip(1.0 / (R + 1e-9), 0, None)
    prec = prec / (prec[valid > 0].mean() + 1e-9)
    feat = np.stack([zc * valid, valid, np.log1p(prec)], axis=0).astype(np.float32)
    return feat, zc, R, valid, mu_hat


# ------------------------------- models -------------------------------
class DilBlock(nn.Module):
    def __init__(self, ci, co, d, causal=False):
        super().__init__()
        self.causal = causal; self.d = d; self.k = 5
        pad = (self.k - 1) * d if causal else (self.k - 1) * d // 2
        self.pad = pad
        self.conv = nn.Conv1d(ci, co, self.k, dilation=d, padding=(0 if causal else pad))
        self.bn = nn.BatchNorm1d(co); self.act = nn.ReLU(inplace=True)
    def forward(self, x):
        if self.causal:
            x = nn.functional.pad(x, (self.pad, 0))     # left-only pad
        return self.act(self.bn(self.conv(x)))


class SmootherTCN(nn.Module):
    """Acausal dilated TCN: 3ch -> smoothed demeaned log-gain path (1ch)."""
    def __init__(self, cin=3, c=64):
        super().__init__()
        ds = [1, 2, 4, 8, 16, 1]
        chs = [cin] + [c] * len(ds)
        self.blocks = nn.ModuleList([DilBlock(chs[i], chs[i + 1], ds[i], causal=False) for i in range(len(ds))])
        self.out = nn.Conv1d(c, 1, 1)
    def forward(self, x):
        for b in self.blocks:
            x = b(x)
        return self.out(x).squeeze(1)


class PredictorTCN(nn.Module):
    """Causal dilated TCN: strict-past features -> (mu, logvar) for z_t (learned
    transition/emission density)."""
    def __init__(self, cin=3, c=64):
        super().__init__()
        ds = [1, 2, 4, 8, 16, 32]
        chs = [cin] + [c] * len(ds)
        self.blocks = nn.ModuleList([DilBlock(chs[i], chs[i + 1], ds[i], causal=True) for i in range(len(ds))])
        self.mu = nn.Conv1d(c, 1, 1); self.lv = nn.Conv1d(c, 1, 1)
    def forward(self, x):
        for b in self.blocks:
            x = b(x)
        return self.mu(x).squeeze(1), self.lv(x).squeeze(1).clamp(-12, 6)


# ------------------------------- data generation -------------------------------
def gen_data(n, seed0):
    rng = np.random.default_rng(4_000_000 + seed0)
    F = np.empty((n, 3, N_EXP), np.float32)     # smoother/predictor input features
    Ltrue = np.empty((n, N_EXP), np.float32)    # demeaned true log-gain (smoother target)
    Z = np.empty((n, N_EXP), np.float32)        # zc (predictor target)
    V = np.empty((n, N_EXP), np.float32)
    for i in range(n):
        x = S.make_scene(rng)
        cls = CLASSES[i % len(CLASSES)]
        a_time, _ = E.sample_class(rng, cls)
        Y = E.simulate(x, a_time, rng)
        feat, zc, R, valid, mu_hat = features_true_scene(Y, x)
        l = np.log(np.maximum(a_time, 1e-12))
        lm = (l * valid).sum() / max(valid.sum(), 1.0)
        F[i] = feat; Ltrue[i] = (l - lm).astype(np.float32); Z[i] = zc.astype(np.float32); V[i] = valid
    return F, Ltrue, Z, V


# ------------------------------- training -------------------------------
def train_models():
    cache = os.path.join(OUT, "json", "e2_data.npz")
    if os.path.exists(cache):
        d = np.load(cache); F, L, Z, V, Fv, Lv, Zv, Vv = (d["F"], d["L"], d["Z"], d["V"],
                                                          d["Fv"], d["Lv"], d["Zv"], d["Vv"])
        log(f"loaded cached data {F.shape}")
    else:
        log("generating E2 training data (procedural scenes x 4 classes)...")
        F, L, Z, V = gen_data(4000, 0)
        Fv, Lv, Zv, Vv = gen_data(1000, 800000)
        np.savez(cache, F=F, L=L, Z=Z, V=V, Fv=Fv, Lv=Lv, Zv=Zv, Vv=Vv)
        log(f"data {F.shape} generated")
    dev = DEVICE
    def t(a): return torch.from_numpy(a).to(dev)
    F, L, Z, V = t(F), t(L), t(Z), t(V)
    Fv, Lv, Zv, Vv = t(Fv), t(Lv), t(Zv), t(Vv)
    n = F.shape[0]; bs = 64

    # ---- smoother ----
    sm = SmootherTCN().to(dev); npar = sum(p.numel() for p in sm.parameters())
    log(f"SmootherTCN params={npar}")
    opt = torch.optim.Adam(sm.parameters(), lr=2e-3)
    sched = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, factor=0.5, patience=5)
    def wmse(pred, tgt, v): return ((pred - tgt) ** 2 * v).sum() / v.sum().clamp_min(1.0)
    best = 1e9; best_state = None; bad = 0
    for ep in range(80):
        sm.train(); perm = torch.randperm(n, device=dev)
        for i in range(0, n, bs):
            idx = perm[i:i + bs]
            pred = sm(F[idx])
            sm_pen = ((pred[:, 1:] - pred[:, :-1]) ** 2).mean()
            loss = wmse(pred, L[idx], V[idx]) + 0.01 * sm_pen
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
            log(f"  [sm] ep{ep:3d} val_wMSE={vl:.5e} (best {best:.5e})")
        if bad >= 12:
            log(f"  [sm] early stop ep{ep}"); break
    sm.load_state_dict(best_state)
    torch.save(dict(state_dict=sm.state_dict(), params=npar, best_val=best), os.path.join(CKPT, "e2_smoother.pt"))

    # ---- predictor (learned transition density) ----
    pr = PredictorTCN().to(dev); npar2 = sum(p.numel() for p in pr.parameters())
    log(f"PredictorTCN params={npar2}")
    opt = torch.optim.Adam(pr.parameters(), lr=2e-3)
    sched = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, factor=0.5, patience=5)
    def shift_right(a):    # inp[t] = a[t-1]; predict target[t] from strict past
        z = torch.zeros_like(a); z[:, :, 1:] = a[:, :, :-1]; return z
    def nll(mu, lv, tgt, v):
        e = (tgt - mu) ** 2 * torch.exp(-lv) + lv
        return (0.5 * e * v).sum() / v.sum().clamp_min(1.0)
    Fp = shift_right(F); Fvp = shift_right(Fv)
    best = 1e9; best_state = None; bad = 0
    for ep in range(80):
        pr.train(); perm = torch.randperm(n, device=dev)
        for i in range(0, n, bs):
            idx = perm[i:i + bs]
            mu, lv = pr(Fp[idx]); loss = nll(mu, lv, Z[idx], V[idx])
            opt.zero_grad(); loss.backward(); opt.step()
        pr.eval()
        with torch.no_grad():
            mu, lv = pr(Fvp); vl = nll(mu, lv, Zv, Vv).item()
        sched.step(vl)
        if vl < best - 1e-6:
            best = vl; best_state = {k: v.detach().cpu().clone() for k, v in pr.state_dict().items()}; bad = 0
        else:
            bad += 1
        if ep % 10 == 0 or bad == 0:
            log(f"  [pr] ep{ep:3d} val_NLL={vl:.5f} (best {best:.5f})")
        if bad >= 12:
            log(f"  [pr] early stop ep{ep}"); break
    pr.load_state_dict(best_state)
    torch.save(dict(state_dict=pr.state_dict(), params=npar2, best_val=best), os.path.join(CKPT, "e2_predictor.pt"))
    return sm, pr


# ------------------------------- deploy / evaluate -------------------------------
def smoother_gain(sm, feat, mu_hat):
    with torch.no_grad():
        l_hat = sm(torch.from_numpy(feat[None]).to(DEVICE)).cpu().numpy().reshape(N_EXP)
    a = np.exp(mu_hat + l_hat); a = a / max(a.mean(), 1e-30)   # E[a]=1 renorm
    return a, l_hat


def predictor_innov(pr, feat, zc, valid):
    inp = np.zeros_like(feat); inp[:, 1:] = feat[:, :-1]        # strict-past shift
    with torch.no_grad():
        mu, lv = pr(torch.from_numpy(inp[None]).to(DEVICE))
    mu = mu.cpu().numpy().reshape(N_EXP); s = np.exp(0.5 * lv.cpu().numpy().reshape(N_EXP))
    r = (zc - mu) / np.maximum(s, 1e-6)
    return r[valid.astype(bool)]


def relerr(h, t): return float(abs(h - t) / max(abs(t), 1e-9))


def evaluate(sm, pr):
    SEEDS = list(range(6)); rows = []
    for cls in CLASSES:
        recs = []
        for sc in DL.SCENES:
            x = DL.scene_x[sc]
            for seed in SEEDS:
                r2 = np.random.default_rng((hash((cls, sc, seed, "e2eval")) % (2**32)))
                a_time, _ = E.sample_class(r2, cls)
                Y = E.simulate(x, a_time, r2)
                feat, zc, R, valid, mu_hat = features_true_scene(Y, x)
                tru = E.true_summary(a_time)

                # OU-model DLGI arm (single-exponential + Kalman smoother + E[a]=1)
                tc_hat, sig_hat = T.mom_autocov(zc, valid)
                if not np.isfinite(tc_hat): tc_hat, sig_hat = 1e3, 1e-3
                phi = np.exp(-1.0 / max(tc_hat, 1e-3))
                ms, Vs = T.kalman_rts(zc, R, valid.astype(bool), phi, max(sig_hat, 1e-6))
                a_ou = np.exp(mu_hat + ms + 0.5 * Vs); a_ou = a_ou / max(a_ou.mean(), 1e-30)
                x_ou = np.clip(DL.arm_A4_gaincorr(Y[0::2], Y[1::2], a_ou[0::2], a_ou[1::2], DL.LAM_TV), 0.0, None)
                innov_ou = E.ou_innovations(zc, R, valid, tc_hat, sig_hat)
                Q, p_ou, _ = E.ljung_box(innov_ou, 20)
                tau_emp, cv_emp = E.emp_autocov_tau(zc, valid)

                # learned-prior DLGI arm
                a_lp, l_hat = smoother_gain(sm, feat, mu_hat)
                x_lp = np.clip(DL.arm_A4_gaincorr(Y[0::2], Y[1::2], a_lp[0::2], a_lp[1::2], DL.LAM_TV), 0.0, None)
                r_lp = predictor_innov(pr, feat, zc, valid)
                Q2, p_lp, _ = E.ljung_box(r_lp, 20)

                # oracle
                x_or = np.clip(DL.arm_A4_gaincorr(Y[0::2], Y[1::2], a_time[0::2], a_time[1::2], DL.LAM_TV), 0.0, None)

                def rmse(a): return float(np.sqrt(np.mean((np.log(np.maximum(a, 1e-9)) - np.log(np.maximum(a_time, 1e-9))) ** 2)))
                recs.append(dict(
                    psnr_ou=DL.psnr(x_ou, x), psnr_lp=DL.psnr(x_lp, x), psnr_or=DL.psnr(x_or, x),
                    corr_ou=DL.gain_path_corr(a_ou, a_time), corr_lp=DL.gain_path_corr(a_lp, a_time),
                    rmse_ou=rmse(a_ou), rmse_lp=rmse(a_lp),
                    p_ou=p_ou, istd_ou=float(innov_ou.std()), pass_ou=E.adequacy_pass(p_ou, float(innov_ou.std())),
                    p_lp=p_lp, istd_lp=float(r_lp.std()), pass_lp=E.adequacy_pass(p_lp, float(r_lp.std())),
                    # medium: OU single-exp tau vs empirical-autocov tau vs true tau
                    tauOU_relerr=relerr((1 + phi) / (1 - phi), tru["tau"]),
                    tauEMP_relerr=relerr(tau_emp, tru["tau"]) if np.isfinite(tau_emp) else np.nan,
                    cv_ou_relerr=relerr(T.cv_of_sigma_l(sig_hat), tru["cv"]),
                    cv_emp_relerr=relerr(cv_emp, tru["cv"]) if np.isfinite(cv_emp) else np.nan))
        def m(k): return float(np.nanmean([r[k] for r in recs]))
        def md(k): return float(np.nanmedian([r[k] for r in recs]))
        row = dict(cls=cls, n=len(recs),
                   psnr_ou=m("psnr_ou"), psnr_lp=m("psnr_lp"), psnr_or=m("psnr_or"),
                   dpsnr_lp_ou=m("psnr_lp") - m("psnr_ou"),
                   corr_ou=m("corr_ou"), corr_lp=m("corr_lp"),
                   rmse_ou=m("rmse_ou"), rmse_lp=m("rmse_lp"),
                   pass_ou=float(np.mean([r["pass_ou"] for r in recs])),
                   pass_lp=float(np.mean([r["pass_lp"] for r in recs])),
                   med_p_ou=md("p_ou"), med_p_lp=md("p_lp"),
                   tauOU_relerr=md("tauOU_relerr"), tauEMP_relerr=md("tauEMP_relerr"),
                   cv_ou_relerr=md("cv_ou_relerr"), cv_emp_relerr=md("cv_emp_relerr"))
        rows.append(row)
        log(f"{cls:10s} PSNR ou={row['psnr_ou']:6.2f} lp={row['psnr_lp']:6.2f} or={row['psnr_or']:6.2f} "
            f"(lp-ou {row['dpsnr_lp_ou']:+.2f})  RMSE {row['rmse_ou']:.3f}->{row['rmse_lp']:.3f}  "
            f"white ou={row['pass_ou']*100:3.0f}% lp={row['pass_lp']*100:3.0f}%  "
            f"tau_relerr OU={row['tauOU_relerr']:.2f} EMP={row['tauEMP_relerr']:.2f}")
    return rows


def main():
    smck = os.path.join(CKPT, "e2_smoother.pt"); prck = os.path.join(CKPT, "e2_predictor.pt")
    if os.path.exists(smck) and os.path.exists(prck):
        sm = SmootherTCN().to(DEVICE); sm.load_state_dict(torch.load(smck, map_location=DEVICE)["state_dict"]); sm.eval()
        pr = PredictorTCN().to(DEVICE); pr.load_state_dict(torch.load(prck, map_location=DEVICE)["state_dict"]); pr.eval()
        log("loaded E2 model ckpts")
    else:
        sm, pr = train_models()
        sm.eval(); pr.eval()
    rows = evaluate(sm, pr)
    meta = dict(exploratory=True, preregistered=False,
                utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
                python=platform.python_version(), torch=torch.__version__, device=DEVICE,
                note="E2(c-e) learned temporal prior (smoother TCN + causal predictor TCN); true-scene z isolation",
                classes=CLASSES)
    json.dump(dict(meta=meta, rows=rows), open(os.path.join(OUT, "json", "e2_learned_prior.json"), "w"), indent=2)
    log("saved e2_learned_prior.json")


if __name__ == "__main__":
    main()
