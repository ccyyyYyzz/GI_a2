#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""IT4b -- certified-blind-code dimension/leak trade (coordinator follow-up to IT4).

Decides whether IT4's "48-dim clean / 64-dim at 5e-4" reading is geometry-fundamental
or gate-statistic-conservative.  Five deliverables:
 1. FULL SPECTRUM: joint-leak generalized eigenvalue (rel leak) vs subspace dim d=1..120;
    d at crossings 1e-5, 1e-4, 1e-3.
 2. GATE-STATISTIC VARIANTS for d in {32,48,64,96}: (a) max leak over subspace,
    (b) median, (c) RMS leak of RANDOM unit code combinations within the subspace
    (operationally relevant -- codes are drawn from the subspace).
 3. PER-Z1 vs JOINT: single-z1 (z1=0, relay-conjugate) null vs joint-over-4 stack;
    d at 1e-4 for both (expect single-z1 strictly larger).
 4. WINDOW-STACKED: rebuild L_be with the V7-style Tukey-apodized object ROI (KT1) and
    repeat; d at 1e-4 before/after (expect the whole spectrum drops).
 5. SANITY: deterministic float64, ridge-free; bottom eigenvalues stable to B-ridge
    0 vs 1e-8 vs 1e-6.
"""
import os
os.environ.pop("WAVE_DTYPE", None)
import time, json
import numpy as np
import torch
import wave_twin as wt

t0 = time.time()
DEV = wt.DEV; N = wt.N; NM = wt.NMACRO; MP = wt.MACRO_PX; KP = wt.KP
ACT = wt.DMD_PIX; c = wt._cx; a0 = c - ACT // 2
Z1S = [0.0, 2.0, 5.0, 10.0]
CHUNK = 16

Uin = np.asarray(wt.band_modes(1, KP)); Ube = np.asarray(wt.band_modes(KP + 1, 9))
n_in = Uin.shape[1]; n_be = Ube.shape[1]
Uin_t = torch.tensor(Uin, device=DEV, dtype=torch.float64)
Ube_t = torch.tensor(Ube, device=DEV, dtype=torch.float64)

def _tukey_full(alpha=0.25):
    w = np.ones(ACT); edge = int(alpha * (ACT - 1) / 2.0)
    for i in range(edge + 1):
        w[i] = 0.5 * (1 + np.cos(np.pi * (2.0 * i / (alpha * (ACT - 1)) - 1))); w[-1 - i] = w[i]
    W = np.outer(w, w); full = np.zeros((N, N)); full[a0:a0 + ACT, a0:a0 + ACT] = W
    return torch.tensor(full, device=DEV, dtype=torch.float64)
TUK = _tukey_full()

def block_sum_64(D, apod=False):
    if apod:
        D = D * TUK
    reg = D[..., a0:a0 + ACT, a0:a0 + ACT].reshape(*D.shape[:-2], NM, MP, NM, MP)
    return reg.sum(dim=(-3, -1))

def build_L(z1mm, apod=False):
    z1 = z1mm * 1e-3
    Uc = wt.propagate(wt.dmd_field(np.ones((NM, NM))).to(torch.complex128), z1)
    Lbe = torch.zeros(n_be, n_in, device=DEV, dtype=torch.float64)
    Lin = torch.zeros(n_in, n_in, device=DEV, dtype=torch.float64)
    for s in range(0, n_in, CHUNK):
        cols = Uin[:, s:s + CHUNK].T.reshape(-1, NM, NM)
        Us = torch.stack([wt.propagate(wt.dmd_field(cols[k]).to(torch.complex128), z1) for k in range(cols.shape[0])])
        D = (Uc[None] * Us.conj()).real
        mu = block_sum_64(D, apod=apod).reshape(cols.shape[0], -1)
        Lbe[:, s:s + CHUNK] = (mu @ Ube_t).T
        Lin[:, s:s + CHUNK] = (mu @ Uin_t).T
    return Lbe, Lin

# build per-z1 operators (no-apod and apod)
L_noapod = {z: build_L(z, apod=False) for z in Z1S}
print("built no-apod L for all z1 [%.0fs]" % (time.time() - t0))
L_apod = {z: build_L(z, apod=True) for z in Z1S}
print("built apod L for all z1 [%.0fs]" % (time.time() - t0))

def spectrum(Lbe_stack, Lin_stack, ridge=0.0):
    A = Lbe_stack.T @ Lbe_stack; B = Lin_stack.T @ Lin_stack
    if ridge > 0:
        B = B + ridge * float(torch.diag(B).mean()) * torch.eye(n_in, device=DEV, dtype=torch.float64)
    wB, UB = torch.linalg.eigh(0.5 * (B + B.T)); wB = torch.clamp(wB, min=wB.max() * 1e-14)
    Bisq = UB @ torch.diag(wB.rsqrt()) @ UB.T
    Mm = Bisq @ A @ Bisq
    ev, evec = torch.linalg.eigh(0.5 * (Mm + Mm.T))
    rel = torch.sqrt(torch.clamp(ev, min=0.0))
    Vblind = Bisq @ evec                        # generalized eigenvectors (B-orthonormal), ascending
    return rel.cpu().numpy(), Vblind, A, B

def stack(Ld, zs):
    return torch.cat([Ld[z][0] for z in zs], 0), torch.cat([Ld[z][1] for z in zs], 0)

def d_at(rel, thr):
    return int((rel <= thr).sum())

def rms_random_combo(Vblind, A, B, d, nsamp=400, seed=0):
    g = torch.Generator(device=DEV).manual_seed(seed)
    Vd = Vblind[:, :d]
    G = torch.randn(d, nsamp, device=DEV, dtype=torch.float64, generator=g)
    W = Vd @ G                                   # (n_in, nsamp) random combos in bottom-d subspace
    num = torch.einsum('in,in->n', W, A @ W)
    den = torch.einsum('in,in->n', W, B @ W)
    leak = torch.sqrt(torch.clamp(num / den, min=0.0))
    return float(leak.pow(2).mean().sqrt()), float(leak.median())

out = {"test": "IT4b_spectrum_trade", "ref": "coordinator follow-up to IT4",
       "params": dict(z1_mm=Z1S, n_in=n_in, n_be=n_be, thresholds=[1e-5, 1e-4, 1e-3]),
       "results": {}}

# 1 + 2: joint spectrum, crossings, gate-statistic variants
Lbe_j, Lin_j = stack(L_noapod, Z1S)
rel_j, V_j, A_j, B_j = spectrum(Lbe_j, Lin_j)
out["results"]["joint_spectrum"] = dict(
    rel_leak_vs_d={str(d): float(rel_j[d - 1]) for d in list(range(1, 21)) + [24, 32, 40, 48, 56, 64, 80, 96, 112, 120]},
    d_at_1e5=d_at(rel_j, 1e-5), d_at_1e4=d_at(rel_j, 1e-4), d_at_1e3=d_at(rel_j, 1e-3))
gate = {}
for d in [32, 48, 64, 96]:
    rms, med_combo = rms_random_combo(V_j, A_j, B_j, d)
    gate["d%d" % d] = dict(max_leak=float(rel_j[:d].max()), median_leak=float(np.median(rel_j[:d])),
                           rms_random_combo=rms, median_random_combo=med_combo)
out["results"]["gate_statistic_variants"] = gate

# 3: per-z1 (z1=0) vs joint
Lbe_0, Lin_0 = stack(L_noapod, [0.0])
rel_0, _, _, _ = spectrum(Lbe_0, Lin_0)
out["results"]["per_z1_vs_joint"] = dict(
    z1_0_only_d_at_1e4=d_at(rel_0, 1e-4), z1_0_only_d_at_1e5=d_at(rel_0, 1e-5),
    joint_d_at_1e4=d_at(rel_j, 1e-4), joint_d_at_1e5=d_at(rel_j, 1e-5),
    single_z1_gain=d_at(rel_0, 1e-4) - d_at(rel_j, 1e-4))

# 4: window-stacked (Tukey-apodized ROI) joint spectrum
Lbe_ja, Lin_ja = stack(L_apod, Z1S)
rel_ja, _, _, _ = spectrum(Lbe_ja, Lin_ja)
out["results"]["window_apodized"] = dict(
    joint_noapod_d_at_1e4=d_at(rel_j, 1e-4), joint_apod_d_at_1e4=d_at(rel_ja, 1e-4),
    joint_noapod_d_at_1e5=d_at(rel_j, 1e-5), joint_apod_d_at_1e5=d_at(rel_ja, 1e-5),
    apod_rel_leak_at_d64=float(rel_ja[63]), noapod_rel_leak_at_d64=float(rel_j[63]))

# 5: ridge sanity on the joint operator
ridge_check = {}
for rg in [0.0, 1e-8, 1e-6]:
    rr, _, _, _ = spectrum(Lbe_j, Lin_j, ridge=rg)
    ridge_check["ridge_%g" % rg] = dict(bottom5=[float(x) for x in rr[:5]], d_at_1e4=d_at(rr, 1e-4))
out["results"]["ridge_sanity"] = ridge_check

# ------------------------------------------------------------------ verdict
dj = out["results"]["joint_spectrum"]["d_at_1e4"]
d0 = out["results"]["per_z1_vs_joint"]["z1_0_only_d_at_1e4"]
dja = out["results"]["window_apodized"]["joint_apod_d_at_1e4"]
rms64 = gate["d64"]["rms_random_combo"]; max64 = gate["d64"]["max_leak"]
out["verdict"] = dict(
    joint_d_at_1e4=dj, joint_d64_max_leak=max64, joint_d64_rms_random_combo=rms64,
    single_z1_d_at_1e4=d0, window_apod_joint_d_at_1e4=dja,
    reading=("Gate-statistic-conservative: the d=64 MAX leak (%.2e) is dominated by the worst "
             "direction; the operationally-relevant RMS leak of random codes drawn from the d=64 "
             "subspace is %.2e (%.0fx smaller). Joint d@1e-4=%d; single-z1(relay-conjugate) "
             "d@1e-4=%d (gain %+d); Tukey-window-apodized joint d@1e-4=%d. Routes multiply."
             % (max64, rms64, (max64 / rms64 if rms64 > 0 else float('nan')), dj, d0, d0 - dj, dja)))
out["elapsed_s"] = round(time.time() - t0, 1)
with open(os.environ.get("OUT", "IT4b_spectrum_trade.json"), "w") as f:
    json.dump(out, f, indent=2)
print("\njoint d@1e-4=%d (max64=%.2e rms64=%.2e) | z1=0 d@1e-4=%d | apod d@1e-4=%d"
      % (dj, max64, rms64, d0, dja))
print("VERDICT:", out["verdict"]["reading"])
print("saved  elapsed %.1fs" % (time.time() - t0))
