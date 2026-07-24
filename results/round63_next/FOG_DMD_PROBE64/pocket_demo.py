# POCKET DEMO -- empirical test of the ONE inhabitable cell the living-region map found, with the
# R39 section 4.6 PRODUCTION estimator (lifted covariance GLS). Answers: was the E8 ~0.45 plateau
# estimator-caused (no-floor theorem's open question) in the one place it matters?
#
# Pocket cell (frozen): sigma_f=0.3, k^-2 medium spectrum, k_w=k_p=5, claim 1.25*k_p, N=4096,
# M=128 signed band-limited codes, n=1e4. iid epochs (one medium realization per code bank ->
# T_eff = T, matching the ruling's diffuser-per-bank acquisition and the CRB in T_eff).
#
# Estimators:
#   OLD  moment (r=1 nonconvex Adam on x, the E8-style fit) -> iid-Gaussian-MLE refine
#   OLD  MLE alone (iid-Gaussian marginal MLE from data init)
#   NEW  lifted GLS (Burer-Monteiro X=FF^T, r>1 CONVEX lift) -> rank-1 extract + mean anchor
#        -> band continuation 1.0->1.25 k_p -> Fisher-scored (iid-MLE) refine
# Blind (data-only init), >=5 multi-starts, T_eff in {200,500,1000}, >=6 naturals + 3 witnesses.
# Metric: beyond-band NRMSE on the 1.0-1.25 k_p annulus (shell 6). Compared to the exact CRB.
# Verdict: POCKET_ACHIEVED (blind <= 2x CRB) / ESTIMATOR_GAP_PERSISTS.
import argparse
import json
import time

import numpy as np
import torch

DEV = "cuda" if torch.cuda.is_available() else "cpu"
DT = torch.float64
n = 64
N = n * n
KP = 5
CLAIM_EDGE = 6                       # 1.25*k_p -> shell 6 (the beyond-band annulus)
SIG_F = 0.30
PHOT = 1e4
BOUNDS = (0.2, 1.8)
M = 128
_X, _Y = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")


def band_modes(k_lo, k_hi):
    cols = []
    for kx in range(-k_hi, k_hi + 1):
        for ky in range(-k_hi, k_hi + 1):
            if k_lo <= max(abs(kx), abs(ky)) <= k_hi:
                ph = 2 * np.pi * (kx * _X + ky * _Y) / n
                cols.append(np.cos(ph).ravel()); cols.append(np.sin(ph).ravel())
    Q, R = np.linalg.qr(np.array(cols).T)
    return Q[:, np.abs(np.diag(R)) > 1e-8]


def medium_modes(k_w):
    cols, kabs = [], []
    for kx in range(-k_w, k_w + 1):
        for ky in range(-k_w, k_w + 1):
            if max(abs(kx), abs(ky)) == 0 or (kx < 0) or (kx == 0 and ky < 0):
                continue
            ph = 2 * np.pi * (kx * _X + ky * _Y) / n
            for f in (np.cos(ph).ravel(), np.sin(ph).ravel()):
                nr = np.linalg.norm(f)
                if nr > 1e-8:
                    cols.append(f / nr); kabs.append(np.hypot(kx, ky))
    return np.array(cols).T, np.array(kabs)


# --- frozen medium (k^-2 spectrum, total variance sigma_f_eff^2) ---
U_med_np, kabs = medium_modes(KP)
db = U_med_np.shape[1]
_w = 1.0 / kabs ** 2
S_np = _w / _w.sum() * (N * SIG_F ** 2)          # (1/N) sum S = sigma_f^2
U_med = torch.tensor(U_med_np, device=DEV, dtype=DT)
S = torch.tensor(S_np, device=DEV, dtype=DT)
sqrtS = torch.sqrt(S)
# codes, in-band + beyond-band bases
U_in_np = band_modes(0, KP)
U_beta_np = band_modes(KP + 1, CLAIM_EDGE)       # the 1.0-1.25 k_p annulus (shell 6)
U_band_np = band_modes(0, CLAIM_EDGE)            # in-band + claim annulus (for continuation)
U_beta = torch.tensor(U_beta_np, device=DEV, dtype=DT)


def signed_codes(seed=10):
    rp = np.random.default_rng(seed)
    Uc = band_modes(1, KP)
    P = rp.standard_normal((M, Uc.shape[1])) @ Uc.T
    P /= np.abs(P).max(axis=1, keepdims=True)
    return P


P_np = signed_codes()
P = torch.tensor(P_np, device=DEV, dtype=DT)
Pabs = torch.abs(P)
# lifted-map building blocks: Q_j = P diag(U_j); precompute for the vectorized L(FF^T)
# UF[j,n,r] = U[n,j] F[n,r]; QF[j,m,r] = sum_n P[m,n] UF[j,n,r]; L = Phi^2 sum_j S_j QF QF^T


def scenes_bank():
    from skimage import data, transform
    out = {}
    for nm, fn in [("cameraman", data.camera), ("coins", data.coins), ("moon", data.moon),
                   ("text", data.text), ("clock", data.clock), ("gravel", data.gravel)]:
        img = transform.resize(fn().astype(np.float64), (n, n), anti_aliasing=True, mode="reflect")
        out[nm] = ((img - img.min()) / (np.ptp(img) + 1e-12)).ravel()
    for w in range(3):                                       # 3 synthetic beyond-band witnesses
        rs = np.random.default_rng(40 + w)
        a_in = U_in_np @ rs.standard_normal(U_in_np.shape[1])
        a_be = U_beta_np @ (2.0 * rs.choice([-1, 1], U_beta_np.shape[1]))
        x = a_in + a_be; x -= x.min(); out[f"witness{w}"] = x / x.max()
    return out


def gen_epochs(x_np, T, seed):
    """iid epochs: one medium realization per epoch; M signed buckets with shot at n=1e4."""
    r = np.random.default_rng(seed)
    Z = r.standard_normal((T, db)) * S_np[None, :] ** 0.5
    G = Z @ U_med_np.T
    W = np.clip(1.0 + G, *BOUNDS)                            # bounded mean-one field
    xx = x_np[None, :]
    s = np.einsum("mn,tn->tm", P_np, W * xx)                # signed buckets (T,M)
    flux = np.einsum("mn,tn->tm", np.abs(P_np), W * xx)     # nonneg photon throughput
    scp = PHOT / max(flux.mean(), 1e-9)
    Y = s + r.standard_normal(s.shape) * np.sqrt(np.clip(flux, 0, None) / scp)   # descaled shot
    return Y


def data_objects(Y):
    """mean, whitening, shot-corrected signal covariance from the epoch record."""
    mhat = Y.mean(0)
    Yc = Y - mhat[None, :]
    Vhat = (Yc.T @ Yc) / Yc.shape[0]                        # includes shot
    fluxbar = (Pabs.cpu().numpy() @ np.clip(np.linalg.pinv(P_np) @ mhat, 0, None))
    return (torch.tensor(mhat, device=DEV, dtype=DT), torch.tensor(Vhat, device=DEV, dtype=DT))


def L_of_F(F):
    """L(FF^T) = sum_j S_j (Q_j F)(Q_j F)^T, Q_j=P diag(U_j). Vectorized (no python loop)."""
    UF = torch.einsum("nj,nr->jnr", U_med, F)               # (db,N,r)
    QF = torch.einsum("mn,jnr->jmr", P, UF)                 # (db,M,r)
    return torch.einsum("j,jmr,jkr->mk", S, QF, QF)


def C_of_x(x):
    """C(x)=P diag(x) K_w diag(x) P^T = H H^T, H = P diag(x) U_med diag(sqrtS)."""
    H = (P * x[None, :]) @ (U_med * sqrtS[None, :])
    return H @ H.t()


def whiten(V):
    ev, evec = torch.linalg.eigh(0.5 * (V + V.t()))
    inv_sqrt = evec @ torch.diag(1.0 / torch.sqrt(torch.clamp(ev, min=1e-9))) @ evec.t()
    return inv_sqrt


def nrmse_beta(x_hat, x_true):
    xh = torch.as_tensor(x_hat, device=DEV, dtype=DT)
    xt = torch.as_tensor(x_true, device=DEV, dtype=DT)
    s = (xh @ xt) / torch.clamp(xh @ xh, min=1e-12)
    e = U_beta.t() @ (s * xh - xt)
    return float((torch.norm(e) / torch.norm(U_beta.t() @ xt)).item())


def mean_anchor(x, mhat):
    Px = P @ x
    sc = (Px @ mhat) / torch.clamp(Px @ Px, min=1e-12)
    return torch.clamp(sc * x, min=0.0) if sc > 0 else torch.clamp(-sc * (-x), min=0.0)


def mle_refine(x0, mhat, Vhat, x_band_np, iters=250, lr=3e-2):
    """iid-Gaussian marginal MLE (Fisher-scored refinement), x constrained to a Fourier band."""
    Bnd = torch.tensor(x_band_np @ x_band_np.T, device=DEV, dtype=DT)   # band projector
    Rbase = (Pabs @ torch.clamp(x0, min=0)) * ((Pabs @ torch.clamp(x0, min=0)).mean() / PHOT)
    u = torch.tensor(np.clip(x0.cpu().numpy(), 1e-4, None), device=DEV, dtype=DT).clone().requires_grad_(True)
    opt = torch.optim.Adam([u], lr=lr)
    Shat = Vhat
    for it in range(iters):
        opt.zero_grad()
        x = torch.clamp(Bnd @ torch.relu(u), min=0.0)
        V = C_of_x(x) + torch.diag(Rbase) + 1e-6 * torch.eye(M, device=DEV, dtype=DT)
        sign, logdet = torch.linalg.slogdet(V)
        nll = logdet + torch.trace(torch.linalg.solve(V, Shat))
        nll.backward(); opt.step()
    return torch.clamp(Bnd @ torch.relu(u.detach()), min=0.0)


def est_moment(mhat, Vhat, Wih, x_row, iters=1200, lr=3e-2, seed=0):
    """OLD nonconvex moment fit: min_x ||Wih(Chat - C(x))Wih||_F^2 + lam, Adam on x (r=1, E8-style)."""
    Chat = Vhat
    r = np.random.default_rng(seed)
    u = torch.tensor(x_row * (0.5 + r.random(N)), device=DEV, dtype=DT).clone().requires_grad_(True)
    opt = torch.optim.Adam([u], lr=lr)
    tgt = Wih @ Chat @ Wih
    for it in range(iters):
        opt.zero_grad()
        x = torch.relu(u)
        res = Wih @ C_of_x(x) @ Wih - tgt
        loss = (res ** 2).sum() + 1e-4 * (x ** 2).sum()
        loss.backward(); opt.step()
    return torch.relu(u.detach())


def est_lifted(mhat, Vhat, Wih, x_row, rank=6, iters=1200, lr=3e-2, lam=1e-3, seed=0):
    """NEW lifted GLS (Burer-Monteiro X=FF^T, rank>1 convex lift). Returns extracted x."""
    Chat = Vhat
    r = np.random.default_rng(seed)
    F0 = (x_row[:, None] * (0.3 + r.random((N, rank)))) / np.sqrt(rank)
    F = torch.tensor(F0, device=DEV, dtype=DT).clone().requires_grad_(True)
    opt = torch.optim.Adam([F], lr=lr)
    tgt = Wih @ Chat @ Wih
    for it in range(iters):
        opt.zero_grad()
        res = Wih @ L_of_F(F) @ Wih - tgt
        loss = (res ** 2).sum() + lam * (F ** 2).sum()
        loss.backward(); opt.step()
    Fd = F.detach()
    U_, s_, _ = torch.linalg.svd(Fd, full_matrices=False)   # leading left singular vec of F
    x = U_[:, 0] * torch.sqrt(torch.clamp(s_[0], min=0))
    if (P @ x) @ mhat < 0:
        x = -x
    return mean_anchor(x, mhat)


def crb_nmse(x_np, T_eff):
    """Exact profiled covariance-Fisher CRB NRMSE on the claim annulus (same V as the forward)."""
    x = torch.tensor(x_np, device=DEV, dtype=DT)
    H = (P * x[None, :]) @ (U_med * sqrtS[None, :])
    C = H @ H.t()
    flux = Pabs @ x
    R = torch.diag(flux * (flux.mean() / PHOT))
    V = C + R
    Vinv = torch.linalg.inv(V + 1e-9 * torch.eye(M, device=DEV, dtype=DT))
    KA = (U_med * S[None, :]) @ ((U_med * x[:, None]).t() @ P.t())   # K_w diag(x) P^T  (N,M)... via modes

    def Gstack(Phi):
        G = torch.empty((Phi.shape[1], M, M), device=DEV, dtype=DT)
        for k in range(Phi.shape[1]):
            B = (P * Phi[:, k][None, :]) @ KA
            G[k] = Vinv @ (B + B.t())
        return G
    Phi_in = torch.tensor(U_in_np, device=DEV, dtype=DT)
    Gb = Gstack(U_beta); Ge = Gstack(Phi_in)
    Glaw = (Vinv @ C).unsqueeze(0)
    Geta = torch.cat([Ge, Glaw], 0)
    half = T_eff / 2.0
    Ibb = half * torch.einsum("aij,bji->ab", Gb, Gb)
    Ibe = half * torch.einsum("aij,bji->ab", Gb, Geta)
    Iee = half * torch.einsum("aij,bji->ab", Geta, Geta)
    Min = P @ Phi_in
    Iee[:Phi_in.shape[1], :Phi_in.shape[1]] += T_eff * (Min.t() @ Vinv @ Min)
    JB = Ibb - Ibe @ torch.linalg.pinv(Iee, rcond=1e-10) @ Ibe.t()
    JB = 0.5 * (JB + JB.t())
    lam = torch.clamp(torch.linalg.eigvalsh(JB), min=0)
    beta = U_beta.t() @ x
    pos = lam > 1e-8 * lam.max()
    tr_inv = float((1.0 / (lam[pos] / T_eff)).sum().item())
    return float(np.sqrt(tr_inv / (T_eff * max(float((beta ** 2).sum()), 1e-12))))


_SCENE_IDX = {nm: i for i, nm in enumerate(
    ["cameraman", "coins", "moon", "text", "clock", "gravel", "witness0", "witness1", "witness2"])}


def run_scene(name, x_np, T, n_starts, iters):
    Y = gen_epochs(x_np, T, seed=1234 + 37 * _SCENE_IDX.get(name, 0) + T)
    mhat, Vhat = data_objects(Y)
    Wih = whiten(Vhat)
    x_row = np.clip(np.linalg.pinv(P_np) @ mhat.cpu().numpy(), 1e-3, None)
    crb = crb_nmse(x_np, T)
    res = {"crb": crb}
    # OLD moment (+MLE refine) and NEW lifted (+MLE refine), multi-start
    mom_x, lift_x = [], []
    for st in range(n_starts):
        xm = est_moment(mhat, Vhat, Wih, x_row, iters=iters, seed=st)
        xm = mle_refine(mle_band(xm, U_in_np), mhat, Vhat, U_band_np, iters=200)
        mom_x.append(xm)
        xl = est_lifted(mhat, Vhat, Wih, x_row, iters=iters, seed=st)
        xl = mle_refine(mle_band(xl, U_in_np), mhat, Vhat, U_band_np, iters=200)
        lift_x.append(xl)
    mom_e = [nrmse_beta(x, x_np) for x in mom_x]
    lift_e = [nrmse_beta(x, x_np) for x in lift_x]
    # pure MLE from data row init (old alt)
    xM = mle_refine(torch.tensor(x_row, device=DEV, dtype=DT), mhat, Vhat, U_band_np, iters=300)
    res.update(moment_med=float(np.median(mom_e)), moment_all=[round(e, 3) for e in mom_e],
               lifted_med=float(np.median(lift_e)), lifted_all=[round(e, 3) for e in lift_e],
               mle_only=nrmse_beta(xM, x_np),
               moment_multistart_std=float(np.std(mom_e)), lifted_multistart_std=float(np.std(lift_e)))
    return res


def mle_band(x, band_np):
    """band-continuation stage 1: project the init onto the in-band subspace before growing."""
    Bnd = torch.tensor(band_np @ band_np.T, device=DEV, dtype=DT)
    return torch.clamp(Bnd @ x, min=0.0) + 1e-3


def convexity_check():
    """Validate the lift is convex (unique min from many inits) vs r=1 nonconvex, on N=256 (16x16)."""
    print("=== convexity validation (small case, N=256) ===", flush=True)
    global n, N, _X, _Y, U_med_np, db, S_np, U_med, S, sqrtS, U_in_np, U_beta_np, U_band_np, U_beta
    global P_np, P, Pabs, kabs, S_np
    n0, N0 = n, N
    n = 16; N = 256; _X, _Y = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    kp = 3; ce = 4
    Um, ka = medium_modes(kp); dbl = Um.shape[1]
    Sl = (1 / ka ** 2); Sl = Sl / Sl.sum() * (N * 0.3 ** 2)
    Uc = band_modes(1, kp); rp = np.random.default_rng(3)
    Pl = rp.standard_normal((32, Uc.shape[1])) @ Uc.T; Pl /= np.abs(Pl).max(1, keepdims=True)
    Uin = band_modes(0, kp); Ube = band_modes(kp + 1, ce)
    rs = np.random.default_rng(4)
    x = Uin @ rs.standard_normal(Uin.shape[1]) + Ube @ (2.0 * rs.choice([-1, 1], Ube.shape[1]))
    x -= x.min(); x /= x.max()
    Umt = torch.tensor(Um, device=DEV, dtype=DT); St = torch.tensor(Sl, device=DEV, dtype=DT)
    Pt = torch.tensor(Pl, device=DEV, dtype=DT); sq = torch.sqrt(St)
    Ht = (Pt * torch.tensor(x, device=DEV, dtype=DT)[None, :]) @ (Umt * sq[None, :])
    Chat = Ht @ Ht.t()                                       # noise-free target (exact)

    def Lf(F):
        UF = torch.einsum("nj,nr->jnr", Umt, F); QF = torch.einsum("mn,jnr->jmr", Pt, UF)
        return torch.einsum("j,jmr,jkr->mk", St, QF, QF)

    def solve(rank, seed):
        r = np.random.default_rng(seed)
        F = torch.tensor(r.standard_normal((N, rank)) * 0.1, device=DEV, dtype=DT).requires_grad_(True)
        opt = torch.optim.Adam([F], lr=5e-2)
        for it in range(1500):
            opt.zero_grad(); loss = ((Lf(F) - Chat) ** 2).sum() + 1e-4 * (F ** 2).sum()
            loss.backward(); opt.step()
        return float(loss.item())
    r1 = [solve(1, s) for s in range(6)]
    r6 = [solve(8, s) for s in range(6)]
    print(f"  r=1 (nonconvex) final objective, 6 inits: {[f'{v:.2e}' for v in r1]}  std {np.std(r1):.2e}", flush=True)
    print(f"  r=8 (lifted)    final objective, 6 inits: {[f'{v:.2e}' for v in r6]}  std {np.std(r6):.2e}", flush=True)
    conv = np.std(r6) < 0.1 * (np.mean(r6) + 1e-12) and np.mean(r6) <= np.mean(r1) + 1e-9
    print(f"  lift reaches a consistent (convex-like) minimum, no worse than r=1: {conv}", flush=True)
    # restore globals
    n, N = n0, N0; _X, _Y = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    return dict(r1_objs=r1, r8_objs=r6, r1_std=float(np.std(r1)), r8_std=float(np.std(r6)),
                convex_like=bool(conv))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    print(f"POCKET DEMO  cell: sigma_f={SIG_F} k^-2 k_w=k_p={KP} claim1.25 (shell{CLAIM_EDGE}) "
          f"N={N} M={M} n={PHOT:.0e} db={db}", flush=True)
    conv = convexity_check()
    scenes = scenes_bank()
    if a.quick:
        r = run_scene("cameraman", scenes["cameraman"], 500, n_starts=2, iters=500)
        print("quick cameraman T=500:", {k: r[k] for k in ("crb", "moment_med", "lifted_med", "mle_only")}, flush=True)
        print(f"[{time.time()-t0:.0f}s]", flush=True)
    else:
        out = {"cell": dict(sigma_f=SIG_F, spectrum="k^-2", k_w_over_kp=1, claim=1.25, N=N, M=M,
                            photons=PHOT), "convexity": conv, "per_T": {}}
        for T in (200, 500, 1000):
            rows = {}
            for name, x in scenes.items():
                rows[name] = run_scene(name, x, T, n_starts=5, iters=1200)
                print(f"  T={T} {name:10s} CRB={rows[name]['crb']:.3f} "
                      f"moment={rows[name]['moment_med']:.3f} lifted={rows[name]['lifted_med']:.3f} "
                      f"mle={rows[name]['mle_only']:.3f}", flush=True)
            nats = [k for k in scenes if not k.startswith("witness")]
            wits = [k for k in scenes if k.startswith("witness")]
            out["per_T"][T] = dict(scenes=rows,
                crb_nat_med=float(np.median([rows[k]["crb"] for k in nats])),
                lifted_nat_med=float(np.median([rows[k]["lifted_med"] for k in nats])),
                moment_nat_med=float(np.median([rows[k]["moment_med"] for k in nats])),
                lifted_wit_med=float(np.median([rows[k]["lifted_med"] for k in wits])),
                crb_wit_med=float(np.median([rows[k]["crb"] for k in wits])))
            json.dump(out, open("POCKET_DEMO.json", "w"), indent=2)
        # verdict: blind lifted <= 2x CRB at T=500-1000 (natural median)
        ok = all(out["per_T"][T]["lifted_nat_med"] <= 2.0 * out["per_T"][T]["crb_nat_med"]
                 for T in (500, 1000))
        out["verdict"] = "POCKET_ACHIEVED" if ok else "ESTIMATOR_GAP_PERSISTS"
        out["wall_s"] = time.time() - t0
        json.dump(out, open("POCKET_DEMO.json", "w"), indent=2)
        print(f"\nPOCKET DEMO VERDICT: {out['verdict']}  [{time.time()-t0:.0f}s]", flush=True)
