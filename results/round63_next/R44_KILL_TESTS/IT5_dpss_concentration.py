#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""IT5 -- DPSS/Slepian concentration window (internal divergence finer #4).

After IT4b showed the naive Tukey window is a KILL (rim effect, 93x worse), the ONLY
live window route is a concentration-optimal (DPSS/Slepian) render kernel + object window.
Test: mean-channel leak into a beyond-2kp witness annulus k in [12,14] with DPSS render +
DPSS object window, deterministic z1 sweep {0,1,5,10,20} mm.  Compare none / Tukey / DPSS.
PASS = DPSS leak <=1e-3 uniform in z1; KILL = the window-alias floor is not
concentration-limited (DPSS no better than Tukey).

Also report the COMBINED route: rebuild the single-z1 (z1=0) SVD leak operator (IT4b) with
the DPSS object window and report d@1e-4 vs the no-window 80 -- does concentration grow the
certified-blind subspace?
"""
import os
os.environ.pop("WAVE_DTYPE", None)
import time, json
import numpy as np
import torch
from scipy.signal.windows import dpss
import wave_twin as wt

t0 = time.time()
DEV = wt.DEV; N = wt.N; NM = wt.NMACRO; MP = wt.MACRO_PX; KP = wt.KP
ACT = wt.DMD_PIX; c = wt._cx; a0 = c - ACT // 2
M = int(os.environ.get("M", "12"))
ZLIST = [0.0, 1.0, 5.0, 10.0, 20.0]
NW = float(os.environ.get("NW", "4"))

U_in = torch.tensor(wt.band_modes(1, KP), device=DEV, dtype=torch.float64)
U_wit = torch.tensor(wt.band_modes(12, 14), device=DEV, dtype=torch.float64)
U_be = torch.tensor(wt.band_modes(KP + 1, 9), device=DEV, dtype=torch.float64)   # for SVD operator

# ---- windows over the 512 object region, embedded in the 2048 grid ----
def embed(W512):
    full = np.zeros((N, N)); full[a0:a0 + ACT, a0:a0 + ACT] = W512
    return torch.tensor(full, device=DEV, dtype=torch.float64)

w_dpss1 = dpss(ACT, NW, Kmax=1)[0]; w_dpss1 = w_dpss1 / w_dpss1.max()
DPSS2 = embed(np.outer(w_dpss1, w_dpss1))
def _tukey(alpha=0.25):
    w = np.ones(ACT); e = int(alpha * (ACT - 1) / 2)
    for i in range(e + 1):
        w[i] = 0.5 * (1 + np.cos(np.pi * (2.0 * i / (alpha * (ACT - 1)) - 1))); w[-1 - i] = w[i]
    return np.outer(w, w)
TUK2 = embed(_tukey())
NONE2 = embed(np.ones((ACT, ACT)))

codes = wt.signed_codes(M=M, seed=10)

def render(code, win_render):
    fld = wt.dmd_field(code).to(torch.complex128)
    if win_render is not None:
        fld = fld * win_render
    return fld

def block_sum_64(D, win_obj):
    if win_obj is not None:
        D = D * win_obj
    reg = D[a0:a0 + ACT, a0:a0 + ACT].reshape(NM, MP, NM, MP)
    return reg.sum(dim=(1, 3))

def jac(win_render, win_obj, z1):
    rows = []
    for i in range(M):
        ap = render((1.0 + codes[i]) * 0.5, win_render)
        am = render((1.0 - codes[i]) * 0.5, win_render)
        if z1 > 0:
            ap = wt.propagate(ap, z1); am = wt.propagate(am, z1)
        D = (ap.real ** 2 + ap.imag ** 2) - (am.real ** 2 + am.imag ** 2)
        rows.append(block_sum_64(D, win_obj).reshape(-1))
    return torch.stack(rows)

def wit_leak(J):
    return float(torch.linalg.norm(J @ U_wit) / (torch.linalg.norm(J @ U_in) + 1e-30))

results = {"test": "IT5_dpss_concentration", "ref": "internal divergence finer #4",
           "params": dict(M=M, z1_mm=ZLIST, witness_shells="12-14", NW=NW),
           "decision": "DPSS witness leak <=1e-3 uniform => concentration route lives; "
                       "KILL if DPSS no better than Tukey (alias floor not concentration-limited)",
           "leak_vs_z1": {}}

configs = [("none", None, None), ("tukey_render+obj", TUK2, TUK2), ("dpss_render+obj", DPSS2, DPSS2)]
for tag, wr, wo in configs:
    row = {}
    for z1mm in ZLIST:
        J = jac(wr, wo, z1mm * 1e-3)
        row["z1_%g" % z1mm] = wit_leak(J)
    results["leak_vs_z1"][tag] = row
    print("[%s] leak(k12-14) " % tag + " ".join("z%g:%.2e" % (z, row["z1_%g" % z]) for z in ZLIST) + " [%.0fs]" % (time.time() - t0))

# ---- COMBINED: single-z1 (z1=0) SVD leak operator with DPSS object window; d@1e-4 ----
def leak_operator(win_obj):
    Uc = wt.dmd_field(np.ones((NM, NM))).to(torch.complex128)   # z1=0
    Uin_np = np.asarray(wt.band_modes(1, KP)); nin = Uin_np.shape[1]
    Ube_np = np.asarray(wt.band_modes(KP + 1, 9)); nbe = Ube_np.shape[1]
    Uin_t = torch.tensor(Uin_np, device=DEV, dtype=torch.float64)
    Ube_t = torch.tensor(Ube_np, device=DEV, dtype=torch.float64)
    Lbe = torch.zeros(nbe, nin, device=DEV, dtype=torch.float64)
    Lin = torch.zeros(nin, nin, device=DEV, dtype=torch.float64)
    for s in range(0, nin, 16):
        cols = Uin_np[:, s:s + 16].T.reshape(-1, NM, NM)
        Us = torch.stack([wt.dmd_field(cols[k]).to(torch.complex128) for k in range(cols.shape[0])])
        D = (Uc[None] * Us.conj()).real
        if win_obj is not None:
            D = D * win_obj
        reg = D[..., a0:a0 + ACT, a0:a0 + ACT].reshape(cols.shape[0], NM, MP, NM, MP)
        mu = reg.sum(dim=(2, 4)).reshape(cols.shape[0], -1)
        Lbe[:, s:s + 16] = (mu @ Ube_t).T; Lin[:, s:s + 16] = (mu @ Uin_t).T
    A = Lbe.T @ Lbe; B = Lin.T @ Lin
    wB, UB = torch.linalg.eigh(0.5 * (B + B.T)); wB = torch.clamp(wB, min=wB.max() * 1e-14)
    Bisq = UB @ torch.diag(wB.rsqrt()) @ UB.T
    ev = torch.linalg.eigvalsh(0.5 * (Bisq @ A @ Bisq + (Bisq @ A @ Bisq).T))
    rel = torch.sqrt(torch.clamp(ev, min=0.0)).cpu().numpy()
    return int((rel <= 1e-4).sum()), int((rel <= 1e-5).sum())

d_none, d5_none = leak_operator(None)
d_dpss, d5_dpss = leak_operator(DPSS2)
results["combined_single_z1_SVD"] = dict(no_window_d_at_1e4=d_none, no_window_d_at_1e5=d5_none,
                                         dpss_d_at_1e4=d_dpss, dpss_d_at_1e5=d5_dpss)
print("[combined z1=0 SVD] d@1e-4 none=%d dpss=%d | d@1e-5 none=%d dpss=%d" % (d_none, d_dpss, d5_none, d5_dpss))

# ---- verdict ----
dpss_leaks = list(results["leak_vs_z1"]["dpss_render+obj"].values())
tuk_leaks = list(results["leak_vs_z1"]["tukey_render+obj"].values())
none_leaks = list(results["leak_vs_z1"]["none"].values())
dpss_ok = max(dpss_leaks) <= 1e-3
dpss_beats_tukey = max(dpss_leaks) < max(tuk_leaks)
results["verdict"] = dict(
    dpss_max_witness_leak=max(dpss_leaks), tukey_max_witness_leak=max(tuk_leaks), none_max_witness_leak=max(none_leaks),
    dpss_improvement_over_none=round(max(none_leaks) / (max(dpss_leaks) + 1e-30), 2),
    combined_dpss_subspace_gain=d_dpss - d_none,
    P_dpss_concentration_lives="PASS" if dpss_ok else ("PARTIAL" if dpss_beats_tukey else "KILL"),
    overall=("PASS -- DPSS render+object window holds the beyond-2kp witness leak <=1e-3 uniformly "
             "in z1 (max %.2e); the hardware concentration route lives and %s the certified-blind "
             "subspace (single-z1 d@1e-4 %d->%d)." % (max(dpss_leaks),
             "grows" if d_dpss > d_none else "does not grow", d_none, d_dpss)
             if dpss_ok else
             "DPSS does not reach <=1e-3 (max %.2e); %s -- report honestly" %
             (max(dpss_leaks), "still beats Tukey" if dpss_beats_tukey else "no better than Tukey (alias floor not concentration-limited = KILL)")))
results["elapsed_s"] = round(time.time() - t0, 1)
with open(os.environ.get("OUT", "IT5_dpss_concentration.json"), "w") as f:
    json.dump(results, f, indent=2)
print("VERDICT:", results["verdict"]["overall"])
print("saved  elapsed %.1fs" % (time.time() - t0))
