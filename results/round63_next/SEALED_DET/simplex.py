#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
THE ATTRIBUTION SIMPLEX (R41 sec 9, Rank 2) -- corrected-shot engine.

Constructs efficient scores for FOUR tangent spaces and uses their canonical angles to bound
cross-class leakage, then builds the five-class classifier used by bar D3.

  tangent 1  mean / in-band scene       : lights the MEAN channel (dm = P U_in) + a lag-0 cov signature
  tangent 2  covariance / beyond-band    : the headline detector direction (lag-0 cov; mean = 0, wall)
  tangent 3  covariance-amplitude        : medium sigma_f direction (lag-0 cov ~ C)
  tangent 4  covariance-lag              : medium correlation-time (tau) direction (LAG-1 cross-bank cov)

Geometry: a whitened JOINT observation space with three blocks -- mean (metric V^-1), lag-0 covariance
(Gaussian metric 1/2 V^-1 (x)_s V^-1), lag-1 cross-covariance (same whitening, weak-correlation
approx). Each tangent = span of its whitened derivative columns across the three blocks.

Reported diagnostics (R41 Rank-2 targets):
  * RAW canonical correlations between the four tangent subspaces (intrinsic class separability);
  * the efficient-score construction (each tangent projected off the union of the other three) --
    by construction the efficient scores are mutually orthogonal (canonical corr = 0), the formal
    "asymptotically independent Gaussian score coordinates" statement;
  * off-diagonal canonical correlations must be < 0.10 (the intrinsic-separability target);
  * five-class balanced accuracy >= 0.90 (operational, on the specificity populations).
"""
import numpy as np
import torch

import sealed_common as sc

DEV = sc.DEV
DT = sc.DT
M = sc.M


def _whiten_cov(G, dV, idx, wt):
    """Whitened half-vech of symmetric dV: <w(A),w(B)> = 1/2 tr(Vinv A Vinv B) with G^T G = Vinv."""
    WdW = G @ dV.cpu().numpy() @ G.T
    return (WdW[idx[0], idx[1]] * wt) / np.sqrt(2.0)


def tangent_blocks(cell):
    """Return whitened derivative-column matrices for the four tangents in the joint observation
    space [mean (M) ; lag-0 cov (L) ; lag-1 cov (L)]. Each is (M+2L, d_family)."""
    Vinv_np = cell["Vinvd"].cpu().numpy()
    G = np.linalg.cholesky(Vinv_np + 1e-15 * np.eye(M)).T          # G^T G = Vinv
    idx = np.tril_indices(M)
    wt = np.where(idx[0] == idx[1], 1.0, np.sqrt(2.0))
    L = idx[0].shape[0]
    rho = cell["rho"]
    drho = cell["drho_dlogtau"]

    def emb(mean_block, dV0, dV1):
        m = np.zeros(M) if mean_block is None else mean_block
        c0 = np.zeros(L) if dV0 is None else _whiten_cov(G, dV0, idx, wt)
        c1 = np.zeros(L) if dV1 is None else _whiten_cov(G, dV1, idx, wt)
        return np.concatenate([G @ m, c0, c1])

    # tangent 1: in-band scene -> mean (P U_in) + lag-0 cov + lag-1 cov (rho*)
    Min = cell["Min"].cpu().numpy()                                # (M, d_eta)
    T1 = []
    for k in range(cell["Vk_in"].shape[0]):
        T1.append(emb(Min[:, k], cell["Vk_in"][k], rho * cell["Vk_in"][k]))
    T1 = np.array(T1).T

    # tangent 2: beyond-band scene -> mean = 0 (wall) + lag-0 cov + lag-1 cov (rho*)
    T2 = []
    for k in range(cell["Vk_b"].shape[0]):
        T2.append(emb(None, cell["Vk_b"][k], rho * cell["Vk_b"][k]))
    T2 = np.array(T2).T

    # tangent 3: amplitude -> lag-0 cov (C) + lag-1 cov (rho*C)
    T3 = emb(None, cell["Vk_amp"][0], rho * cell["Vk_amp"][0])[:, None]

    # tangent 4: lag / tau -> lag-0 = 0 (amplitude fixed) + lag-1 cov (drho * C)
    T4 = emb(None, None, cell["Vk_lag1"][0])[:, None]

    return dict(inband=T1, beyond=T2, amplitude=T3, lag=T4, dim=T1.shape[0])


def canonical_correlations(Ua, Ub):
    """Canonical correlations between two subspaces spanned by columns of Ua, Ub (singular values
    of Qa^T Qb, Q = orthonormal basis). Returns sorted-descending array."""
    Qa, _ = np.linalg.qr(Ua)
    Qb, _ = np.linalg.qr(Ub)
    s = np.linalg.svd(Qa.T @ Qb, compute_uv=False)
    return np.clip(s, 0, 1)


def simplex_gram(cell):
    """Full 4x4 canonical-correlation structure (max off-diagonal canonical correlation per pair) +
    the efficient-score orthogonality check. Returns a dict for SIMPLEX_CALIBRATION.md."""
    T = tangent_blocks(cell)
    fams = ["inband", "beyond", "amplitude", "lag"]
    U = {f: T[f] for f in fams}
    # RAW pairwise canonical correlations
    raw = {}
    max_off = 0.0
    for i in range(4):
        for j in range(i + 1, 4):
            cc = canonical_correlations(U[fams[i]], U[fams[j]])
            raw[f"{fams[i]}__{fams[j]}"] = dict(max_cc=float(cc.max()), mean_cc=float(cc.mean()),
                                                n_modes=int(len(cc)))
            max_off = max(max_off, float(cc.max()))
    # EFFICIENT scores: project each tangent off the union of the other three; report residual
    # canonical correlation (should collapse toward 0 -> asymptotically independent coordinates)
    eff = {}
    max_off_eff = 0.0
    for i in range(4):
        others = np.concatenate([U[fams[j]] for j in range(4) if j != i], axis=1)
        Qo, _ = np.linalg.qr(others)
        Ri = U[fams[i]] - Qo @ (Qo.T @ U[fams[i]])          # efficient (nuisance-profiled) tangent
        # residual energy fraction retained after profiling (efficiency of the class score)
        retained = float(np.linalg.norm(Ri) / (np.linalg.norm(U[fams[i]]) + 1e-30))
        # canonical corr of the efficient score with each OTHER raw tangent (leakage after profiling)
        leak = 0.0
        for j in range(4):
            if j == i:
                continue
            if np.linalg.norm(Ri) < 1e-12:
                continue
            cc = canonical_correlations(Ri, U[fams[j]])
            leak = max(leak, float(cc.max()))
        eff[fams[i]] = dict(efficiency_retained=round(retained, 4), max_leak_cc=round(leak, 5))
        max_off_eff = max(max_off_eff, leak)
    return dict(raw_canonical_correlations=raw, max_offdiag_raw_cc=round(max_off, 5),
                efficient_scores=eff, max_offdiag_efficient_cc=round(max_off_eff, 6),
                joint_dim=int(T["dim"]),
                target_max_offdiag=sc.BARS["D3"]["simplex_canon_corr_max"])


# ------------------------------------------------------------------ five-class classifier (bar D3)
def five_class_populations(cell, T_eff, n_per, eps_beyond, rng0):
    """Generate the five specificity populations and score each with the FOUR statistics
    (t_mean, t_beyond, t_amp, t_lag). Returns a dict class -> (n_per, 4) feature array."""
    import arms
    xn = np.sqrt(cell["xnorm2"])
    # beyond-band target direction
    c_b = sc.energy_spread_delta(cell["db"], eps_beyond, xn, seed=int(rng0) + 7)
    delta_b = cell["Phi_b"].cpu().numpy() @ c_b
    # in-band scene direction, matched norm
    Uin1 = sc.band_modes(1, sc.KP)
    rin = np.random.default_rng(int(rng0) + 8)
    c_in = rin.standard_normal(Uin1.shape[1]); c_in /= np.linalg.norm(c_in)
    delta_in = Uin1 @ (c_in * eps_beyond * xn)
    m0 = (cell["Pt"] @ cell["x"]).cpu().numpy()
    Vinv = cell["Vinvd"].cpu().numpy()

    def feats(Sc, Mb, Lag1):
        t_beyond = arms.fixed_cov_score(cell, Sc, c_b)
        t_amp = arms.amplitude_score(cell, Sc)
        t_lag = arms.lag_score(cell, Lag1)
        dmv = Mb.cpu().numpy() - m0[None, :]
        t_mean = T_eff * np.einsum("ai,ij,aj->a", dmv, Vinv, dmv)
        return np.stack([t_mean, t_beyond, t_amp, t_lag], axis=1)

    pops = {}
    specs = {
        "H0":       dict(kind="H0", delta=None, sf=1.0, tau=1.0),
        "inband":   dict(kind="inband", delta=delta_in, sf=1.0, tau=1.0),
        "beyond":   dict(kind="beyond", delta=delta_b, sf=1.0, tau=1.0),
        "amplitude": dict(kind="H0", delta=None, sf=1.2, tau=1.0),
        "lag":      dict(kind="H0", delta=None, sf=1.0, tau=2.0),
    }
    for i, (cls, sp) in enumerate(specs.items()):
        Sc, Mb, Lag1 = sc.gen_records(cell, n_per, T_eff, sp["kind"], sp["delta"],
                                      sf_scale=sp["sf"], tau_scale=sp["tau"],
                                      rng=int(rng0) + 100 + 10 * i, want_lag=True)
        pops[cls] = feats(Sc, Mb, Lag1)
    return pops


def balanced_accuracy(pops):
    """LDA (pooled-within-class-covariance Mahalanobis nearest-centroid) in the 4-D score space ->
    per-class recall + balanced accuracy over the five classes. LDA is the natural attribution
    classifier: it accounts for the different scales/correlations of the four score statistics."""
    classes = list(pops.keys())
    K = len(classes)
    X = np.concatenate([pops[c] for c in classes], axis=0)
    y = np.concatenate([[i] * len(pops[c]) for i, c in enumerate(classes)])
    mu = X.mean(0); sd = X.std(0) + 1e-12
    Xs = (X - mu) / sd
    cents = np.array([Xs[y == i].mean(0) for i in range(K)])
    # pooled within-class covariance
    Sw = np.zeros((Xs.shape[1], Xs.shape[1]))
    for i in range(K):
        Xi = Xs[y == i] - cents[i]
        Sw += Xi.T @ Xi
    Sw = Sw / len(y) + 1e-6 * np.eye(Xs.shape[1])
    Swi = np.linalg.inv(Sw)
    diff = Xs[:, None, :] - cents[None, :, :]                       # (n, K, d)
    d = np.einsum("nkd,de,nke->nk", diff, Swi, diff)               # Mahalanobis
    pred = d.argmin(1)
    recalls = {}
    for i, c in enumerate(classes):
        m = y == i
        recalls[c] = float((pred[m] == i).mean())
    bal = float(np.mean(list(recalls.values())))
    conf = np.zeros((K, K), int)
    for t, p in zip(y, pred):
        conf[t, p] += 1
    return dict(balanced_accuracy=round(bal, 4), per_class_recall={k: round(v, 3) for k, v in recalls.items()},
                classes=classes, confusion=conf.tolist())


if __name__ == "__main__":
    import time
    t0 = time.time()
    cell = sc.setup_cell(**sc.BEST_CELL)
    g = simplex_gram(cell)
    print("=== SIMPLEX GRAM (best cell) ===", flush=True)
    print(f"joint dim = {g['joint_dim']}", flush=True)
    for pair, v in g["raw_canonical_correlations"].items():
        print(f"  raw cc  {pair:28s} max={v['max_cc']:.4f} mean={v['mean_cc']:.4f} ({v['n_modes']} modes)", flush=True)
    print(f"  max off-diagonal RAW cc       = {g['max_offdiag_raw_cc']:.4f}", flush=True)
    print(f"  max off-diagonal EFFICIENT cc = {g['max_offdiag_efficient_cc']:.6f} (target < {g['target_max_offdiag']})", flush=True)
    for f, v in g["efficient_scores"].items():
        print(f"    eff[{f:10s}] retained={v['efficiency_retained']:.3f} leak_cc={v['max_leak_cc']:.5f}", flush=True)
    print("=== FIVE-CLASS CLASSIFIER (best cell, T=1200, eps_beyond=2%) ===", flush=True)
    pops = five_class_populations(cell, 1200, 200, 0.02, rng0=1000)
    ba = balanced_accuracy(pops)
    print(f"  balanced accuracy = {ba['balanced_accuracy']:.3f}", flush=True)
    print(f"  per-class recall  = {ba['per_class_recall']}", flush=True)
    print(f"[simplex selftest done {time.time()-t0:.0f}s]", flush=True)
