#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""TLSG PART A -- Three-Ledger Scaling experiment (R45 sec B11.7).

Nuisance-profiled Gaussian covariance family: one bank Y ~ N(0, Sigma_eps), Sigma_0=I
(WLOG after whitening).  Efficient symmetric tangent basis B_1..B_p with (1/2)tr(B_jB_k)
=delta_jk; Sigma_eps = I + sum_j theta_j B_j, theta=sqrt(lambda) u, lambda=||theta||^2
=(1/2)||A_perp||_F^2.  Score s_j(Y)=(1/2)[Y^T B_j Y - tr B_j]; E_0 s=0, Cov_0(s)=I_p;
E_A s=theta.  Z_T = T^{-1/2} sum_t s(Y_t) ~ N(sqrt(T) theta, I_p).

Three statistics per cell (M, p_eff):
  1. matched known-direction score nu = u^T Z_T  (Ledger I; power at T lambda ~ 1, p^0)
  2. direction-agnostic score norm Q = ||Z_T||^2  (Ledger II; power at T lambda ~ sqrt(p))
  3. full-tangent estimator theta_hat = Z_T/sqrt(T), risk E||theta_hat-theta||^2 = p/T (Ledger III; T lambda ~ p)

FROZEN BARS (R45 B11.7, verbatim -- recorded per bar):
  - matched 50%-power contours collapse vs T*lambda, CV<=20%;
  - blind contours collapse vs T*lambda/sqrt(p), CV<=25% AND fail to collapse vs T*lambda alone;
  - relative estimation risk collapses vs T*lambda/p, CV<=25%;
  - fitted p-exponents = 0+-0.15, 1/2+-0.15, 1+-0.15;
  - the KT6 operating point stays below blind power 0.2 at FPR 0.05 under the optimal blind statistic.
Full simulation of the actual quadratic-form scores (NOT the LAN reduction assumed a priori).
"""
import os, time, json
import numpy as np
import torch
from scipy.stats import ncx2, chi2

t0 = time.time()
DEV = "cuda" if torch.cuda.is_available() else "cpu"
RDT = torch.float64
LAM = float(os.environ.get("LAM", "0.15"))
N_MC = int(os.environ.get("N_MC", "1500"))
CHUNK = int(os.environ.get("CHUNK", "32"))
SEED = int(os.environ.get("SEED", "20260724"))
FPR = 0.05
# cells: each p in {10,36,136,528} at >=2 M; each M in {8,16,32,64} at >=2 p; p_eff<=M(M+1)/2
CELLS = [(8, 10), (16, 10), (8, 36), (32, 36), (16, 136), (64, 136), (32, 528), (64, 528)]
rng = np.random.default_rng(SEED)

def sym_basis(M):
    """standard orthonormal symmetric basis under (1/2)tr(AB): sqrt2*E_ii and (E_ij+E_ji)."""
    G = []
    for i in range(M):
        A = np.zeros((M, M)); A[i, i] = np.sqrt(2.0); G.append(A)
    for i in range(M):
        for j in range(i + 1, M):
            A = np.zeros((M, M)); A[i, j] = 1.0; A[j, i] = 1.0; G.append(A)
    return np.stack(G)                       # (M(M+1)/2, M, M), (1/2)tr(G_a G_b)=delta_ab

def make_cell(M, p):
    G = sym_basis(M); dim = G.shape[0]
    Q, _ = np.linalg.qr(rng.standard_normal((dim, p)))    # random p-dim subspace (orthonormal cols)
    B = np.einsum('ap,amn->pmn', Q[:, :p], G)             # (p,M,M) orthonormal tangent basis
    u = rng.standard_normal(p); u /= np.linalg.norm(u)
    theta = np.sqrt(LAM) * u
    A_perp = np.einsum('p,pmn->mn', theta, B)
    Sig = np.eye(M) + A_perp
    assert np.linalg.eigvalsh(Sig).min() > 1e-6, "Sigma_eps not PSD"
    trB = np.einsum('pmm->p', B)
    return (torch.tensor(B, device=DEV, dtype=RDT), torch.tensor(trB, device=DEV, dtype=RDT),
            torch.tensor(u, device=DEV, dtype=RDT), torch.tensor(theta, device=DEV, dtype=RDT),
            torch.tensor(np.linalg.cholesky(Sig), device=DEV, dtype=RDT), torch.eye(M, device=DEV, dtype=RDT))

def scores_ZT(chol, B, trB, Tgrid, gen):
    """Z_T at each T in Tgrid, for n_mc replicates. chol: (M,M) cholesky of Sigma."""
    M = B.shape[1]; p = B.shape[0]; Tmax = int(Tgrid[-1])
    ZT = torch.zeros(N_MC, len(Tgrid), p, device=DEV, dtype=RDT)
    for s in range(0, N_MC, CHUNK):
        c = min(CHUNK, N_MC - s)
        Zr = torch.randn(c, Tmax, M, device=DEV, dtype=RDT, generator=gen)
        Y = Zr @ chol.T                                   # (c,Tmax,M) ~ N(0,Sigma)
        # scores s_j = 1/2 (Y^T B_j Y - trB_j): einsum in two steps
        BY = torch.einsum('pmk,ctk->ctpm', B, Y)          # (c,Tmax,p,M)
        quad = torch.einsum('ctm,ctpm->ctp', Y, BY)       # (c,Tmax,p)
        sc = 0.5 * (quad - trB[None, None])               # (c,Tmax,p)
        cs = torch.cumsum(sc, dim=1)
        for it, T in enumerate(Tgrid):
            ZT[s:s + c, it] = cs[:, T - 1] / np.sqrt(T)
        del Zr, Y, BY, quad, sc, cs
    return ZT

def power_curve(stat_null, stat_alt):
    thr = torch.quantile(stat_null, 1 - FPR)
    return float((stat_alt > thr).float().mean()), float(thr)

def find_T50(Tgrid, powers):
    """interpolate T where power crosses 0.5 (log-linear in T)."""
    p = np.array(powers); Tg = np.array(Tgrid, float)
    idx = np.where(p >= 0.5)[0]
    if len(idx) == 0:
        return None
    i = idx[0]
    if i == 0:
        return float(Tg[0])
    x0, x1 = np.log(Tg[i - 1]), np.log(Tg[i]); y0, y1 = p[i - 1], p[i]
    lt = x0 + (0.5 - y0) * (x1 - x0) / (y1 - y0 + 1e-12)
    return float(np.exp(lt))

results = {"test": "TLSG_partA_three_ledger", "R45_ref": "sec B11.7",
           "params": dict(lambda_=LAM, n_mc=N_MC, cells_M_p=CELLS, FPR=FPR, seed=SEED),
           "cells": []}

for (M, p) in CELLS:
    B, trB, u, theta, chol, I = make_cell(M, p)
    Tmax = int(min(1600, max(60, 8 * 2.33 * np.sqrt(p) / LAM)))
    Tgrid = sorted(set(int(round(x)) for x in np.unique(np.geomspace(4, Tmax, 16).astype(int)) if x >= 2))
    gen = torch.Generator(device=DEV).manual_seed(SEED + M * 1000 + p)
    ZT_null = scores_ZT(chol * 0 + I, B, trB, Tgrid, gen)   # null: Sigma=I -> chol=I
    ZT_alt = scores_ZT(chol, B, trB, Tgrid, gen)
    matched_pw = []; blind_pw = []; rel_risk = []
    for it, T in enumerate(Tgrid):
        zn = ZT_null[:, it]; za = ZT_alt[:, it]
        mp, _ = power_curve(zn @ u, za @ u)                # matched
        bp, _ = power_curve((zn ** 2).sum(1), (za ** 2).sum(1))  # blind score-norm
        risk = float(((za / np.sqrt(T) - theta[None]) ** 2).sum(1).mean())
        matched_pw.append(mp); blind_pw.append(bp); rel_risk.append(risk / LAM)
    T50m = find_T50(Tgrid, matched_pw); T50b = find_T50(Tgrid, blind_pw)
    cell = dict(M=M, p_eff=p, Tgrid=Tgrid, matched_power=matched_pw, blind_power=blind_pw,
                rel_risk=rel_risk, Tlam_50_matched=(T50m * LAM if T50m else None),
                Tlam_50_blind=(T50b * LAM if T50b else None),
                risk_times_T=[rel_risk[it] * LAM * Tgrid[it] for it in range(len(Tgrid))])  # = p if risk=p/T
    results["cells"].append(cell)
    print("[M=%d p=%d] Tlam50 matched=%s blind=%s | blind/sqrt(p)=%s | mean(risk*T*lam)=%.1f(p=%d) [%.0fs]"
          % (M, p, ("%.2f" % (T50m * LAM)) if T50m else "NA", ("%.2f" % (T50b * LAM)) if T50b else "NA",
             ("%.2f" % (T50b * LAM / np.sqrt(p))) if T50b else "NA",
             float(np.mean(cell["risk_times_T"])), p, time.time() - t0))

# ------------------------------------------------------------------ bar evaluation
def cv(vals):
    v = np.array([x for x in vals if x is not None and np.isfinite(x)])
    return float(np.std(v) / (np.mean(v) + 1e-30)) if len(v) else None

matched_c = [c["Tlam_50_matched"] for c in results["cells"]]
blind_c = [c["Tlam_50_blind"] for c in results["cells"]]
ps = [c["p_eff"] for c in results["cells"]]
blind_over_sqrtp = [b / np.sqrt(pp) if b else None for b, pp in zip(blind_c, ps)]

# Bar 1: matched collapse vs T*lambda
bar1_cv = cv(matched_c)
# Bar 2: blind collapse vs T*lambda/sqrt(p) AND fails vs T*lambda
bar2_cv_scaled = cv(blind_over_sqrtp)
bar2_cv_raw = cv(blind_c)      # should be HIGH (fails to collapse vs T*lambda)
# Bar 3: relative risk collapses vs T*lambda/p -> relative_risk*(T*lambda/p)=1; use risk*T*lam/p ~ 1
risk_ratio = []
for c in results["cells"]:
    for it in range(len(c["Tgrid"])):
        if c["Tgrid"][it] * LAM > 1:      # avoid tiny-T transient
            risk_ratio.append(c["rel_risk"][it] * (c["Tgrid"][it] * LAM / c["p_eff"]))
bar3_cv = cv(risk_ratio)
# Bar 4: p-exponents. matched: log(Tlam50m) vs log p -> 0; blind: -> 0.5; estimation: log(mean risk*T) vs log p -> 1
def slope(xp, yv):
    x = np.log(np.array(xp)); y = np.log(np.clip(np.array(yv), 1e-9, None))
    return float(np.polyfit(x, y, 1)[0])
exp_matched = slope(ps, [c if c else np.nan for c in matched_c])
exp_blind = slope(ps, [c if c else np.nan for c in blind_c])
est_pTerm = [float(np.mean(c["risk_times_T"])) for c in results["cells"]]   # = p if CRB
exp_est = slope(ps, est_pTerm)

# Bar 5: KT6 operating point blind power < 0.2 at FPR 0.05 under optimal blind statistic
KT6_SIGNAL_TO_FLOOR = 5.25e-3     # worst near-contact from committed KT6 M=16 run
KT6_P = 136                       # M=16 covariance-score tangent dim = M(M+1)/2
def blind_power(Tlam, p):
    thr = chi2.ppf(1 - FPR, df=p)
    return float(1 - ncx2.cdf(thr, df=p, nc=Tlam))
# detectability index of the optimal blind (score-norm) test = T*lambda/sqrt(2p); KT6 signal/floor maps to this
kt6_Tlam = KT6_SIGNAL_TO_FLOOR * np.sqrt(2 * KT6_P)
kt6_blind_power = blind_power(kt6_Tlam, KT6_P)
kt6_blind_power_10x = blind_power(10 * kt6_Tlam, KT6_P)
kt6_blind_power_100x = blind_power(100 * kt6_Tlam, KT6_P)

bars = {
    "bar1_matched_collapse_vs_Tlambda": dict(CV=bar1_cv, threshold=0.20, contours=matched_c,
                                             PASS=bool(bar1_cv is not None and bar1_cv <= 0.20)),
    "bar2_blind_collapse_vs_Tlambda_over_sqrtp": dict(CV_scaled=bar2_cv_scaled, threshold=0.25,
        CV_raw_vs_Tlambda=bar2_cv_raw, contours_over_sqrtp=blind_over_sqrtp,
        PASS=bool(bar2_cv_scaled is not None and bar2_cv_scaled <= 0.25 and bar2_cv_raw > 0.25)),
    "bar3_risk_collapse_vs_Tlambda_over_p": dict(CV=bar3_cv, threshold=0.25, mean_risk_ratio=float(np.mean(risk_ratio)),
                                                PASS=bool(bar3_cv is not None and bar3_cv <= 0.25)),
    "bar4_p_exponents": dict(matched=round(exp_matched, 3), blind=round(exp_blind, 3), estimation=round(exp_est, 3),
        targets=[0.0, 0.5, 1.0], tol=0.15,
        PASS=bool(abs(exp_matched) <= 0.15 and abs(exp_blind - 0.5) <= 0.15 and abs(exp_est - 1.0) <= 0.15)),
    "bar5_KT6_below_blind_0p2": dict(kt6_signal_to_floor=KT6_SIGNAL_TO_FLOOR, kt6_p=KT6_P, kt6_Tlambda=kt6_Tlam,
        blind_power=round(kt6_blind_power, 4), blind_power_10x=round(kt6_blind_power_10x, 4),
        blind_power_100x=round(kt6_blind_power_100x, 4), threshold=0.2,
        PASS=bool(kt6_blind_power < 0.2 and kt6_blind_power_10x < 0.2)),
}
results["bars"] = bars
results["partA_all_pass"] = bool(all(b["PASS"] for b in bars.values()))
results["elapsed_s"] = round(time.time() - t0, 1)
print("\n=== PART A BARS ===")
for k, b in bars.items():
    print(" ", k, "PASS" if b["PASS"] else "FAIL")
print("exponents matched/blind/est = %.3f / %.3f / %.3f (targets 0/0.5/1)" % (exp_matched, exp_blind, exp_est))
print("KT6 blind power = %.4f (10x %.4f, 100x %.4f)" % (kt6_blind_power, kt6_blind_power_10x, kt6_blind_power_100x))
with open(os.environ.get("OUT", "TLSG_partA.json"), "w") as f:
    json.dump(results, f, indent=2)
print("partA_all_pass:", results["partA_all_pass"], " elapsed %.1fs" % (time.time() - t0))
