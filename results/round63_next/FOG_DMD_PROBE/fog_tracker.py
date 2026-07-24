# fog_tracker.py -- GPU blind joint (x, {z_t}) estimators for the fog-as-2nd-DMD probe.
#
# Model:  B_{t} = P (w_t (.) x) + shot noise,   w_t = 1 + U z_t,   z_t OU-correlated.
# Given x, the per-epoch obs is linear-Gaussian in z_t:
#     r_t = B_t - P x = G z_t + v_t,   G = P diag(x) U  (M x d_w),  v_t ~ N(0, R_t)
# with OU dynamics z_t = rho z_{t-1} + innov,  innov ~ N(0, (1-rho^2) sig_w^2 I).
# => the E-step over z is EXACTLY an RTS Kalman smoother (time-invariant G, per-epoch R_t).
#     The M-step over x is a weighted ridge given {w_t}.  This is tracking MAP-EM.
# Engine B is a joint-Adam optimizer over (x, Z) for comparison.
import numpy as np
import torch

# float64: this null-recovery metric measures small null-space components and is precision
# sensitive (float32 degraded even the oracle). Vectorized build_H + chunked einsum keep the
# fp64 cost acceptable on the RTX 4060 for N<=4096.
DTYPE = torch.float64


def _to(t, dev, dtype=DTYPE):
    return torch.as_tensor(np.asarray(t), device=dev, dtype=dtype)


# ---------------- x M-step (clean/unweighted): Hadamard-factored, O(T N^2) ----------------
# For uniform weights,  A^T A = sum_t (P (.) W_t)^T (P (.) W_t) = (P^T P) (.) (W^T W).
# PtP is precomputed once; PtBt = P^T B^T is constant across ALS iters.
def solve_x_clean(PtP, PtBt, W, lam_x):
    N = PtP.shape[0]
    WtW = W.t() @ W                          # (N,N)
    AtA = PtP * WtW
    Atb = (PtBt * W.t()).sum(dim=1)          # (N,)
    ridge = lam_x * torch.trace(AtA) / N
    AtA = AtA + ridge * torch.eye(N, device=PtP.device, dtype=PtP.dtype)
    return torch.linalg.solve(AtA, Atb)


# ---------------- x M-step: weighted ridge given W (chunked, memory-safe) ----------------
def solve_x_weighted(P, W, B, wts, lam_x, dev, chunk=64):
    """P:(M,N) W:(T,N) B:(T,M) wts:(T,M) -> x:(N,).  Solves (A^T Wgt A + lam I)x = A^T Wgt b."""
    M, N = P.shape
    T = W.shape[0]
    AtA = torch.zeros((N, N), device=dev, dtype=P.dtype)
    Atb = torch.zeros((N,), device=dev, dtype=P.dtype)
    for s in range(0, T, chunk):
        e = min(s + chunk, T)
        A = P[None, :, :] * W[s:e, None, :]            # (c,M,N)
        w = wts[s:e]                                    # (c,M)
        Aw = A * w[:, :, None]                          # (c,M,N)
        AtA += torch.einsum('cmn,cmk->nk', Aw, A)
        Atb += torch.einsum('cmn,cm->n', Aw, B[s:e])
    ridge = lam_x * torch.trace(AtA) / N
    AtA += ridge * torch.eye(N, device=dev, dtype=P.dtype)
    x = torch.linalg.solve(AtA, Atb)
    return x


# ---------------- z E-step: RTS Kalman smoother (time-invariant G) ----------------
def rts_smoother(G, r, R_diag, rho, sig_w, dev, prec_scale=1.0):
    """G:(M,d) r:(T,M) obs residuals  R_diag:(T,M) obs noise var  -> Zsmooth:(T,d).
    OU: z_t = rho z_{t-1} + innov, innov var q=(1-rho^2)sig_w^2; z_0 ~ N(0, sig_w^2 I)."""
    M, d = G.shape
    T = r.shape[0]
    var_mult = 1.0 / max(prec_scale, 1e-6)                 # weaken prior -> larger variances
    # rho may be scalar or per-component (d,): diagonal OU transition A = diag(rho)
    if torch.is_tensor(rho):
        rho_v = rho.to(device=dev, dtype=G.dtype)
    else:
        rho_v = torch.as_tensor(np.broadcast_to(np.asarray(rho, float), (d,)).copy(),
                                device=dev, dtype=G.dtype)
    q = (1.0 - rho_v ** 2) * sig_w ** 2 * var_mult        # (d,)
    Qmat = torch.diag(q)
    rr = rho_v[:, None] * rho_v[None, :]                  # (d,d) for cov propagation
    I = torch.eye(d, device=dev, dtype=G.dtype)
    Gt = G.t()
    zf = torch.zeros((T, d), device=dev, dtype=G.dtype)
    Pf = torch.zeros((T, d, d), device=dev, dtype=G.dtype)
    zp = torch.zeros((T, d), device=dev, dtype=G.dtype)
    Pp = torch.zeros((T, d, d), device=dev, dtype=G.dtype)
    z_pred = torch.zeros(d, device=dev, dtype=G.dtype)
    P_pred = (sig_w ** 2 * var_mult) * I
    for t in range(T):
        zp[t] = z_pred; Pp[t] = P_pred
        Rt = torch.diag(R_diag[t])
        S = G @ P_pred @ Gt + Rt
        K = torch.linalg.solve(S, G @ P_pred).t()
        innov = r[t] - G @ z_pred
        zf[t] = z_pred + K @ innov
        Pf[t] = (I - K @ G) @ P_pred
        z_pred = rho_v * zf[t]
        P_pred = rr * Pf[t] + Qmat
    zs = torch.zeros((T, d), device=dev, dtype=G.dtype)
    zs[-1] = zf[-1]
    for t in range(T - 2, -1, -1):
        C = (Pf[t] * rho_v[None, :]) @ torch.linalg.inv(Pp[t + 1])   # Pf[t] A^T Pp^{-1}
        zs[t] = zf[t] + C @ (zs[t + 1] - zp[t + 1])
    return zs


def spectral_init_Z(P_np, U_np, B_np, sig_w, dev='cuda'):
    """Data-driven Z init from SVD of mean-removed bucket fluctuations.
    Model: b_{i,t} = B_{i,t}-mean_t = (U^T diag(P_i) x) . z_t  => b = H(x) Z^T (rank d_w).
    SVD gives Z up to a rotation; rescaled to std ~ sig_w. Breaks the trivial z=0 saddle."""
    dev = torch.device(dev if torch.cuda.is_available() else 'cpu')
    B = _to(B_np, dev); d = U_np.shape[1]
    b = B - B.mean(dim=0, keepdim=True)          # (T,M) mean over epochs removed
    b = b.t()                                    # (M,T)
    Uu, S, Vh = torch.linalg.svd(b, full_matrices=False)
    k = min(d, S.shape[0])
    Zc = (Vh[:k].t() * S[:k])                    # (T,k) temporal factor up to rotation/scale
    Zc = Zc / (Zc.std() + 1e-12) * sig_w         # rescale to plausible medium magnitude
    Z = torch.zeros((B.shape[0], d), device=dev, dtype=B.dtype)
    Z[:, :k] = Zc
    return Z.detach().cpu().numpy()


def solve_spectral_rot(P_np, U_np, B_np, wts_np, sig_w, lam_x, d_keep=None,
                       n_als=60, n_restart=10, dev='cuda', seed0=0, refine_em=0,
                       rho=0.0, R_diag_np=None):
    """SPECTRAL-ROTATION blind solver.  The mean-removed fluctuation matrix b=B-mean_t
    is rank d_w with b = H(x) Z^T, so the TRUE temporal coeffs Z live in the known subspace
    spanned by the top right singular vectors Zhat (T x d).  Recovery reduces to a d_w x d_w
    mixing Q (Z = Zhat Q) plus x -- a tiny well-conditioned bilinear:
        x-step (given Q): weighted ridge (linear).
        Q-step (given x): CLOSED FORM  (Zhat^T Zhat) Q (H^T H) = Zhat^T (B-Px) H,  H=[P(U[:,k](.)x)]_k.
    Random Q restarts pick the lowest data residual.  Optional OU-Kalman EM refinement after."""
    dev = torch.device(dev if torch.cuda.is_available() else 'cpu')
    P = _to(P_np, dev); U = _to(U_np, dev); B = _to(B_np, dev); wts = _to(wts_np, dev)
    M, N = P.shape; d = U.shape[1]; T = B.shape[0]
    dk = d if d_keep is None else min(d_keep, d)
    clean = bool(np.allclose(wts_np, wts_np.flat[0]))   # uniform weights -> fast path
    PtP = P.t() @ P
    PtBt = P.t() @ B.t()                            # (N,T)
    # temporal subspace from SVD of fluctuations
    b = (B - B.mean(dim=0, keepdim=True))          # (T,M)
    Uu, S, Vh = torch.linalg.svd(b.t(), full_matrices=False)  # b^T = (M,T)
    Zhat = Vh[:dk].t()                              # (T,dk) orthonormal temporal basis
    ZtZ = Zhat.t() @ Zhat                           # (dk,dk) ~ I
    ZtZ_inv = torch.linalg.inv(ZtZ + 1e-9 * torch.eye(dk, device=dev, dtype=P.dtype))
    def build_H(x):
        # H[:,k] = P (U[:,k] (.) x)   -> (M,d)   (single matmul, vectorized over k)
        return P @ (U * x[:, None])

    def x_from_Q(Q):
        Z = Zhat @ Q                                # (T,d)
        W = torch.clamp(1.0 + Z @ U.t(), min=0.05)
        if clean:
            x = solve_x_clean(PtP, PtBt, W, lam_x)
        else:
            x = solve_x_weighted(P, W, B, wts, lam_x, dev)
        return x, W

    best = None
    for r in range(n_restart):
        torch.manual_seed(seed0 + r)
        Q = 0.1 * torch.randn(dk, d, device=dev, dtype=P.dtype)
        for it in range(n_als):
            x, W = x_from_Q(Q)
            H = build_H(x)                          # (M,d)
            HtH = H.t() @ H + 1e-8 * torch.eye(d, device=dev, dtype=P.dtype)
            HtH_inv = torch.linalg.inv(HtH)
            y = B - (P @ x)[None, :]                # (T,M)
            RHS = Zhat.t() @ y @ H                  # (dk,d)
            Q = ZtZ_inv @ RHS @ HtH_inv             # closed-form Q-step
        x, W = x_from_Q(Q)
        pred = torch.einsum('mn,tn->tm', P, W * x[None, :])
        resid = (wts * (pred - B) ** 2).sum().item()
        if best is None or resid < best[0]:
            best = (resid, x.detach().clone(), Q.detach().clone())
    resid, x, Q = best
    # optional OU-Kalman EM refinement seeded from the spectral solution
    if refine_em > 0 and R_diag_np is not None:
        Zc = (Zhat @ Q).detach().cpu().numpy()
        x = _to(track_em(P_np, U_np, B_np, wts_np, rho, sig_w, R_diag_np, lam_x,
                         n_em=refine_em, anneal=False, Z_init=Zc, dev=str(dev)), dev)
    return x.detach().cpu().numpy(), float(resid)


def kalman_em(P_np, U_np, B_np, R_diag_np, rho, sig_w, lam_x, n_em=30, dev='cuda',
              x_init=None, return_nll=False):
    """E7a -- MARGINAL-LIKELIHOOD MLE of x via EM. The centered model is linear-Gaussian
    state space: z_t OU, b_t = H(x) z_t + v_t, H = P diag(x) U. The Kalman recursion gives the
    z-marginalized likelihood exactly. EM: E-step = RTS smoother returning E[z_t] AND E[z_t z_t^T]
    (the posterior SECOND moment -- NOT a point estimate -- is what marginalizes z and blocks the
    pathwise collusion); M-step = closed-form nonneg x-solve using those moments. Uses ALL temporal
    structure optimally (vs E5a's few lag moments) -> statistically efficient. Returns x (numpy)."""
    dev = torch.device(dev if torch.cuda.is_available() else 'cpu')
    P = _to(P_np, dev); U = _to(U_np, dev); B = _to(B_np, dev); R_diag = _to(R_diag_np, dev)
    M, N = P.shape; d = U.shape[1]; T = B.shape[0]
    b = B - B.mean(dim=0, keepdim=True)                     # (T,M) centered
    PtP = P.t() @ P; Ptb = P.t() @ b.t()                    # (N,N), (N,T)
    q = (1.0 - rho ** 2) * sig_w ** 2
    Id = torch.eye(d, device=dev, dtype=P.dtype)
    # init x
    if x_init is None:
        x = torch.clamp(torch.linalg.pinv(P) @ B.mean(0), min=1e-3)
    else:
        x = _to(x_init, dev)
    for it in range(n_em):
        H = P @ (U * x[:, None])                            # (M,d)
        Ht = H.t()
        zf = torch.zeros((T, d), device=dev, dtype=P.dtype)
        Pf = torch.zeros((T, d, d), device=dev, dtype=P.dtype)
        zp = torch.zeros((T, d), device=dev, dtype=P.dtype)
        Pp = torch.zeros((T, d, d), device=dev, dtype=P.dtype)
        z_pred = torch.zeros(d, device=dev, dtype=P.dtype); P_pred = (sig_w ** 2) * Id
        nll = 0.0
        for t in range(T):
            zp[t] = z_pred; Pp[t] = P_pred
            S = H @ P_pred @ Ht + torch.diag(R_diag[t])     # (M,M)
            e = b[t] - H @ z_pred
            K = torch.linalg.solve(S, H @ P_pred).t()       # (d,M)
            zf[t] = z_pred + K @ e
            Pf[t] = P_pred - K @ (H @ P_pred)
            if return_nll:
                sign, logdet = torch.linalg.slogdet(S)
                nll = nll + logdet + e @ torch.linalg.solve(S, e)
            z_pred = rho * zf[t]; P_pred = rho ** 2 * Pf[t] + q * Id
        # backward RTS (means + covariances)
        zs = zf.clone(); Ps = Pf.clone()
        for t in range(T - 2, -1, -1):
            C = Pf[t] * rho @ torch.linalg.inv(Pp[t + 1])   # (d,d)
            zs[t] = zf[t] + C @ (zs[t + 1] - zp[t + 1])
            Ps[t] = Pf[t] + C @ (Ps[t + 1] - Pp[t + 1]) @ C.t()
        Mbar = (Ps + torch.einsum('ti,tj->tij', zs, zs)).sum(dim=0)   # sum_t E[z z^T] (d,d)
        # M-step (fast, clean-style; obs-noise handled in the E-step): x-solve
        LHS = PtP * (U @ Mbar @ U.t())                      # (N,N) = (P^T P) (.) (U Mbar U^T)
        Uzs = U @ zs.t()                                    # (N,T)
        RHS = (Uzs * Ptb).sum(dim=1)                        # (N,)
        ridge = lam_x * torch.trace(LHS) / N
        x = torch.clamp(torch.linalg.solve(LHS + ridge * torch.eye(N, device=dev, dtype=P.dtype), RHS), min=0.0)
    if return_nll:
        return x.detach().cpu().numpy(), float(nll.item())
    return x.detach().cpu().numpy()


def cov_estimate(P_np, U_np, B_np, coeff_sd, rho, lags=(1, 2, 4, 8, 16, 32), n_starts=6,
                 steps=5000, lr=1.5e-2, dev='cuda', x_row=None, seed0=0, rankproj=False,
                 That_override=None):
    """E5a COVARIANCE-DOMAIN estimator (R38 eq 1.18). Never estimates z_t -> nothing to collude.
    With the declared law K_w = coeff_sd^2 U U^T, the model covariance factors as
        G_P = P diag(x) K_w diag(x) P^T = coeff_sd^2 H(x) H(x)^T,   H(x)=P diag(x) U,
    which is ROTATION-INVARIANT (the medium rotation is marginalized, not a free parameter).
    Lag-l bucket covariances S_l = (1/T) sum_t b_t b_{t+l}^T (b_t centered) obey S_l ~ r_l G_P for
    l>=1 (r_l=rho^l), and are SHOT-NOISE-FREE (shot is white in time -> only S_0). Target Ghat is
    the SNR-weighted combine of S_l/r_l over l>=1. Estimator: min_{x>=0} ||H(x)H(x)^T - That||_F^2,
    Adam multi-start (data-only). Returns per-start x + objectives for the drift/agreement test."""
    dev = torch.device(dev if torch.cuda.is_available() else 'cpu')
    P = _to(P_np, dev); U = _to(U_np, dev); B = _to(B_np, dev)
    M, N = P.shape; d = U.shape[1]; T = B.shape[0]
    b = B - B.mean(dim=0, keepdim=True)                      # (T,M)
    num = torch.zeros((M, M), device=dev, dtype=P.dtype); wsum = 0.0
    for l in lags:
        if l >= T - 4:
            continue
        r_l = rho ** l
        if r_l < 1e-3:
            continue
        Sl = (b[:T - l].t() @ b[l:]) / (T - l)
        Sl = 0.5 * (Sl + Sl.t())
        w_l = r_l ** 2                                       # SNR weight (downweight noisy far lags)
        num = num + w_l * (Sl / r_l); wsum += w_l
    Ghat = num / max(wsum, 1e-12)                            # ~ coeff_sd^2 H H^T
    That = Ghat / (coeff_sd ** 2)                            # ~ H_true H_true^T
    if That_override is not None:
        That = _to(That_override, dev)
    if rankproj:                                             # denoise: G_P is rank <= d
        w, V = torch.linalg.eigh(That); w = torch.clamp(w, min=0); w[:-d] = 0
        That = (V * w) @ V.t()
    Tnorm = (That ** 2).sum() + 1e-12
    xr0 = _to(x_row, dev) if x_row is not None else torch.clamp(torch.linalg.pinv(P) @ B.mean(0), min=0)
    starts, objs = [], []
    for s in range(n_starts):
        torch.manual_seed(seed0 + s)
        if s == 0:
            u = torch.log(torch.clamp(xr0, min=1e-3)).clone()   # Stage-A init
        else:
            u = torch.log(torch.clamp(xr0 * torch.rand(N, device=dev, dtype=P.dtype) * 2, min=1e-3))
        u.requires_grad_(True)
        opt = torch.optim.Adam([u], lr=lr)
        for it in range(steps):
            opt.zero_grad()
            x = torch.nn.functional.softplus(u)
            H = P @ (U * x[:, None])                         # (M,d)
            resid = H @ H.t() - That
            loss = (resid ** 2).sum() / Tnorm
            loss.backward(); opt.step()
        x = torch.nn.functional.softplus(u).detach()
        starts.append(x.cpu().numpy()); objs.append(float(loss.item()))
    return dict(per_start_x=starts, objs=objs)


def refine_lognormal_reduced(P_np, U_np, B_np, z_init, coeff_sd, lam_x=1e-4,
                             outer=40, qsteps=40, qlr=5e-3, dev='cuda'):
    """Model-MATCHED reduced refinement for lognormal media. Medium coeffs constrained to the
    centered-SVD temporal subspace via a d x d mixing Q (z = Zh Q; ~d^2 params, stable, unlike
    full-z Adam). Forward model w=exp(Uz - v/2) matches the generator, so the true Q is a genuine
    residual minimum (no linear-mismatch drift). Alternates closed-form nonneg x-solve with a few
    Adam steps on Q (x detached -> no backprop through the solve, ~100x faster). Returns x (numpy)."""
    dev = torch.device(dev if torch.cuda.is_available() else 'cpu')
    P = _to(P_np, dev); U = _to(U_np, dev); B = _to(B_np, dev)
    M, N = P.shape; d = U.shape[1]; T = B.shape[0]
    PtP = P.t() @ P; PtBt = P.t() @ B.t(); Ut = U.t()
    var_n = (coeff_sd ** 2) * (U ** 2).sum(dim=1)
    Bc = B - B.mean(dim=0, keepdim=True)
    _, _, Vh = torch.linalg.svd(Bc.t(), full_matrices=False); Zh = Vh[:d].t()
    zt = _to(z_init, dev); zt = (zt - zt.mean(0, keepdim=True)) / (zt.std(0, keepdim=True) + 1e-12) * coeff_sd
    Q = torch.linalg.lstsq(Zh, zt).solution.clone().requires_grad_(True)
    opt = torch.optim.Adam([Q], lr=qlr)
    x = None
    for it in range(outer):
        with torch.no_grad():
            W = torch.exp((Zh @ Q) @ Ut - 0.5 * var_n[None, :])
            x = torch.clamp(solve_x_clean(PtP, PtBt, W, lam_x), min=0.0)
        for _ in range(qsteps):                       # Q gradient steps, x fixed (detached)
            opt.zero_grad()
            W = torch.exp((Zh @ Q) @ Ut - 0.5 * var_n[None, :])
            pred = torch.einsum('mn,tn->tm', P, W * x[None, :])
            loss = ((pred - B) ** 2).sum()
            loss.backward(); opt.step()
    with torch.no_grad():
        W = torch.exp((Zh @ Q) @ Ut - 0.5 * var_n[None, :])
        x = torch.clamp(solve_x_clean(PtP, PtBt, W, lam_x), min=0.0)
    return x.detach().cpu().numpy()


def refine_from_z(P_np, U_np, B_np, wts_np, z_init, lam_x=1e-4, iters=150,
                  nonneg=True, dev='cuda'):
    """PURE-LIKELIHOOD tight refinement of (x, Q) from a given medium-coefficient init z_init
    (T x d, any rotation/scale).  Reduced ALS in the centered-SVD temporal subspace: x-step =
    nonneg LS + small Tikhonov, Q-step = closed form. No homotopy, no learned terms. Used to
    carry a good SOBI/DL initialization to the likelihood optimum. Returns x (numpy)."""
    dev = torch.device(dev if torch.cuda.is_available() else 'cpu')
    P = _to(P_np, dev); U = _to(U_np, dev); B = _to(B_np, dev)
    M, N = P.shape; d = U.shape[1]; T = B.shape[0]
    PtP = P.t() @ P; PtBt = P.t() @ B.t(); Ut = U.t()
    Bc = B - B.mean(dim=0, keepdim=True)
    _, _, Vh = torch.linalg.svd(Bc.t(), full_matrices=False)
    Zh = Vh[:d].t()
    ZtZ_inv = torch.linalg.inv(Zh.t() @ Zh + 1e-9 * torch.eye(d, device=dev, dtype=P.dtype))
    zt = _to(z_init, dev)
    zt = zt / (zt.std(dim=0, keepdim=True) + 1e-12) * (B.std().item() * 0 + 1.0)  # scale-free
    Q = ZtZ_inv @ (Zh.t() @ zt)
    x = None
    for it in range(iters):
        W = torch.clamp(1.0 + (Zh @ Q) @ Ut, min=0.05)
        x = solve_x_clean(PtP, PtBt, W, lam_x)
        if nonneg:
            x = torch.clamp(x, min=0.0)
        H = P @ (U * x[:, None]); HtH = H.t() @ H + 1e-8 * torch.eye(d, device=dev, dtype=P.dtype)
        y = B - (P @ x)[None, :]
        Q = ZtZ_inv @ (Zh.t() @ y @ H) @ torch.linalg.inv(HtH)
    return x.detach().cpu().numpy()


def staged_blind(P_np, U_np, B_np, wts_np, rho, sig_w, R_diag_np, lam_x,
                 dev='cuda', seeds=(0, 1, 2, 3, 4), n_als=120,
                 lam_sched=(3e-1, 3e-2, 3e-3), nonneg=True, refine_em=25,
                 ou_prior_scale=1.0, d_keep=None, sobi_z=None):
    """Staged data-only cold-start blind solver (ruling 3.1). ZERO truth access.
      Stage A: nonnegative row-space scene from temporal-mean buckets (known P).
      Stage B: centered spectral init  B_c = B - mean_t = P diag(x) U Z_c^T  -> temporal
               subspace Zhat (SVD); reduce to a d x d mixing Q (Z = Zhat Q).
      Stage C: continuation ALS with nonnegativity + Tikhonov homotopy (high->low ridge,
               which forbids the blown-up ||x|| spurious minima), then OU-Kalman EM polish.
    Multi-start over `seeds` (data-only random Q); pick lowest data residual. ou_prior_scale
    multiplies the OU precision (set <1 to weaken the temporal prior for sensitivity tests).
    Returns dict(x=best, resid, per_seed_x=list, Zhat)."""
    dev = torch.device(dev if torch.cuda.is_available() else 'cpu')
    P = _to(P_np, dev); U = _to(U_np, dev); B = _to(B_np, dev); wts = _to(wts_np, dev)
    M, N = P.shape; d = U.shape[1]; T = B.shape[0]
    dk = d if d_keep is None else min(d_keep, d)
    # x M-step uses the fast Hadamard-factored unweighted-LS + Tikhonov solve (O(T N^2)).
    # Poisson statistics enter via the noisy buckets B and the Kalman R_diag (E-step); using
    # unweighted LS here is a dev-grade approximation that keeps the whole grid ~minutes.
    PtP = P.t() @ P; PtBt = P.t() @ B.t()
    Ut = U.t()
    def xsolve(W, lam):
        x = solve_x_clean(PtP, PtBt, W, lam)
        return torch.clamp(x, min=0.0) if nonneg else x
    # Stage A: nonneg row-space estimate from mean buckets
    Bbar = B.mean(dim=0)                                   # (M,)
    Pinv = torch.linalg.pinv(P)
    xA = torch.clamp(Pinv @ Bbar, min=0.0)
    # Stage B: centered spectral temporal subspace
    Bc = B - B.mean(dim=0, keepdim=True)                   # (T,M)
    _, _, Vh = torch.linalg.svd(Bc.t(), full_matrices=False)
    Zhat = Vh[:dk].t()                                     # (T,dk)
    ZtZ = Zhat.t() @ Zhat
    ZtZ_inv = torch.linalg.inv(ZtZ + 1e-9 * torch.eye(dk, device=dev, dtype=P.dtype))
    # E4a: SOBI rotation-fix seed. Q_sobi maps the SVD temporal basis onto the recovered
    # independent-source coordinates -> a near-truth Q init (up to perm/sign/scale) instead
    # of a random one.  scale each source to the declared per-coeff sd.
    Q_sobi = None
    if sobi_z is not None:
        Zs = _to(sobi_z, dev)
        Zs = Zs / (Zs.std(dim=0, keepdim=True) + 1e-12) * sig_w
        Q_sobi = ZtZ_inv @ (Zhat.t() @ Zs)                # (dk,d) s.t. Zhat Q_sobi ~ Zs

    def als_from(Q):
        x = xA.clone()
        for lam in lam_sched:
            for it in range(n_als // len(lam_sched)):
                Z = Zhat @ Q
                W = torch.clamp(1.0 + Z @ Ut, min=0.05)
                x = xsolve(W, lam)
                H = P @ (U * x[:, None])                   # (M,d)
                HtH = H.t() @ H + 1e-8 * torch.eye(d, device=dev, dtype=P.dtype)
                y = B - (P @ x)[None, :]
                Q = ZtZ_inv @ (Zhat.t() @ y @ H) @ torch.linalg.inv(HtH)
        return x, Q

    def data_resid(x, Q):
        W = torch.clamp(1.0 + (Zhat @ Q) @ Ut, min=0.05)
        pred = torch.einsum('mn,tn->tm', P, W * x[None, :])
        return float((wts * (pred - B) ** 2).sum().item()), W

    per_seed = []
    best = None
    for sd in seeds:
        torch.manual_seed(1000 + sd)
        if Q_sobi is not None:
            Q0 = Q_sobi + (0.0 if sd == seeds[0] else 0.03) * torch.randn(dk, d, device=dev, dtype=P.dtype)
        else:
            Q0 = 0.1 * torch.randn(dk, d, device=dev, dtype=P.dtype)
        x, Q = als_from(Q0)
        # Stage C: OU-Kalman EM polish seeded from spectral solution
        if refine_em > 0 and R_diag_np is not None:
            Zc = (Zhat @ Q).detach().cpu().numpy()
            xp = track_em(P_np, U_np, B_np, wts_np, rho, sig_w, R_diag_np, lam_sched[-1],
                          n_em=refine_em, anneal=False, Z_init=Zc, dev=str(dev),
                          nonneg=nonneg, ou_prior_scale=ou_prior_scale)
            x = _to(xp, dev)
            # recover Q consistent with polished x for residual accounting
            H = P @ (U * x[:, None]); HtH = H.t() @ H + 1e-8 * torch.eye(d, device=dev, dtype=P.dtype)
            y = B - (P @ x)[None, :]
            Q = ZtZ_inv @ (Zhat.t() @ y @ H) @ torch.linalg.inv(HtH)
        r, _ = data_resid(x, Q)
        xn = x.detach().cpu().numpy()
        per_seed.append(xn)
        if best is None or r < best[0]:
            best = (r, xn)
    return dict(x=best[1], resid=best[0], per_seed_x=per_seed,
                Zhat=Zhat.detach().cpu().numpy())


def track_em(P_np, U_np, B_np, wts_np, rho, sig_w, R_diag_np, lam_x,
             n_em=15, anneal=True, dev='cuda', x_init=None, Z_init=None,
             null_kick=0.0, nonneg=False, ou_prior_scale=1.0, verbose=False):
    """Blind tracking MAP-EM. If Z_init given, start with M-step (x given that medium).
    nonneg clamps the scene estimate >=0 each M-step; ou_prior_scale multiplies the OU
    precision (set <1 to weaken the temporal prior for the sensitivity test)."""
    dev = torch.device(dev if torch.cuda.is_available() else 'cpu')
    P = _to(P_np, dev); U = _to(U_np, dev); B = _to(B_np, dev); wts = _to(wts_np, dev)
    R_diag = _to(R_diag_np, dev)
    M, N = P.shape; d = U.shape[1]; T = B.shape[0]
    PtP = P.t() @ P; PtBt = P.t() @ B.t()
    def xsolve(W):
        x = solve_x_clean(PtP, PtBt, W, lam_x)   # fast unweighted LS + Tikhonov (see staged_blind)
        return torch.clamp(x, min=0.0) if nonneg else x
    if Z_init is not None:
        W = torch.clamp(1.0 + _to(Z_init, dev) @ U.t(), min=0.05)
        x = xsolve(W)
    elif x_init is None:
        # cold init: x = fixed-A weighted solution (W = 1 everywhere)
        x = xsolve(torch.ones((T, N), device=dev, dtype=P.dtype))
    else:
        x = _to(x_init, dev)
    if null_kick > 0:
        g = torch.randn(N, device=dev, dtype=P.dtype)
        Pinv_range = torch.linalg.pinv(P) @ (P @ g)      # range part of g
        gnull = g - Pinv_range
        x = x + null_kick * torch.linalg.norm(x) / torch.linalg.norm(gnull) * gnull
    for it in range(n_em):
        # optional sig_w annealing: start smaller (trust prior/continuity), relax to true
        if anneal:
            frac = min(1.0, (it + 1) / max(1, n_em // 2))
            sig_it = sig_w * (0.4 + 0.6 * frac)
        else:
            sig_it = sig_w
        # E-step: smoother for z given current x
        G = P @ (U * x[:, None])                          # (M,d)
        r = B - (P @ x)[None, :]                          # (T,M)
        Zs = rts_smoother(G, r, R_diag, rho, sig_it, dev, prec_scale=ou_prior_scale)
        W = 1.0 + Zs @ U.t()
        W = torch.clamp(W, min=0.05)
        # M-step: x given W
        x = xsolve(W)
    return x.detach().cpu().numpy()


# ---------------- Engine B: joint Adam ----------------
def joint_adam(P_np, U_np, B_np, wts_np, rho, sig_w, lam_x,
               steps=4000, lr=2e-2, dev='cuda', z_scale_init=0.0, seed=0, Z_init=None):
    dev = torch.device(dev if torch.cuda.is_available() else 'cpu')
    torch.manual_seed(seed)
    P = _to(P_np, dev); U = _to(U_np, dev); B = _to(B_np, dev); wts = _to(wts_np, dev)
    M, N = P.shape; d = U.shape[1]; T = B.shape[0]
    x = torch.zeros(N, device=dev, dtype=P.dtype, requires_grad=True)
    if Z_init is not None:
        Z = _to(Z_init, dev).clone().requires_grad_(True)
    else:
        Z = (z_scale_init * torch.randn(T, d, device=dev, dtype=P.dtype)).requires_grad_(True)
    opt = torch.optim.Adam([x, Z], lr=lr)
    q = (1.0 - rho ** 2) * sig_w ** 2 if rho > 0 else sig_w ** 2
    for step in range(steps):
        opt.zero_grad()
        W = 1.0 + Z @ U.t()                               # (T,N)
        pred = torch.einsum('mn,tn->tm', P, W * x[None, :])  # (T,M)
        misfit = (wts * (pred - B) ** 2).sum()
        # OU prior on Z
        if rho > 0:
            innov = Z[1:] - rho * Z[:-1]
            ou = (innov ** 2).sum() / q + (Z[0] ** 2).sum() / (sig_w ** 2)
        else:
            ou = (Z ** 2).sum() / (sig_w ** 2)
        loss = 0.5 * misfit + 0.5 * lam_x * (x ** 2).sum() + 0.5 * ou
        loss.backward()
        opt.step()
    return x.detach().cpu().numpy()


def joint_adam_reg(P_np, U_np, B_np, wts_np, rho, sig_w, lam_x,
                   steps=6000, lr=2e-2, dev='cuda', seed=0, Z_init=None, z_scale_init=0.0,
                   nonneg_x=True, nonneg_w=True, lognormal=False):
    """Physically-regularized joint solver: nonnegativity on the scene x and medium w
    (both are intensities), Tikhonov lam_x on x, OU prior on z.  These priors penalize the
    blown-up spurious minima (||x|| explodes) that the plain-LS blind objective admits."""
    dev = torch.device(dev if torch.cuda.is_available() else 'cpu')
    torch.manual_seed(seed)
    P = _to(P_np, dev); U = _to(U_np, dev); B = _to(B_np, dev); wts = _to(wts_np, dev)
    M, N = P.shape; d = U.shape[1]; T = B.shape[0]
    # x parametrized via softplus for nonnegativity (else raw)
    xr = torch.zeros(N, device=dev, dtype=P.dtype, requires_grad=True)
    if Z_init is not None:
        Z = _to(Z_init, dev).clone().requires_grad_(True)
    else:
        Z = (z_scale_init * torch.randn(T, d, device=dev, dtype=P.dtype)).requires_grad_(True)
    opt = torch.optim.Adam([xr, Z], lr=lr)
    q = (1.0 - rho ** 2) * sig_w ** 2 if rho > 0 else sig_w ** 2
    # lognormal-matched medium model: w_t = exp(U z_t - v/2), v = per-pixel var of Uz.
    if lognormal:
        var_n = (sig_w ** 2) * (U ** 2).sum(dim=1)          # per-pixel var of g=Uz (sig_w=coeff sd)
    for step in range(steps):
        opt.zero_grad()
        x = torch.nn.functional.softplus(xr) if nonneg_x else xr
        if lognormal:
            W = torch.exp(Z @ U.t() - 0.5 * var_n[None, :])
        else:
            W = 1.0 + Z @ U.t()
            if nonneg_w:
                W = torch.clamp(W, min=0.05)
        pred = torch.einsum('mn,tn->tm', P, W * x[None, :])
        misfit = (wts * (pred - B) ** 2).sum()
        if rho > 0:
            innov = Z[1:] - rho * Z[:-1]
            ou = (innov ** 2).sum() / q + (Z[0] ** 2).sum() / (sig_w ** 2)
        else:
            ou = (Z ** 2).sum() / (sig_w ** 2)
        loss = 0.5 * misfit + 0.5 * lam_x * (x ** 2).sum() + 0.5 * ou
        loss.backward(); opt.step()
    x = (torch.nn.functional.softplus(xr) if nonneg_x else xr).detach().cpu().numpy()
    return x


def solve_oracle(P_np, W_np, B_np, wts_np, lam_x, dev='cuda'):
    dev = torch.device(dev if torch.cuda.is_available() else 'cpu')
    P = _to(P_np, dev); W = _to(W_np, dev); B = _to(B_np, dev); wts = _to(wts_np, dev)
    if bool(np.allclose(wts_np, wts_np.flat[0])):
        x = solve_x_clean(P.t() @ P, P.t() @ B.t(), W, lam_x)
    else:
        x = solve_x_weighted(P, W, B, wts, lam_x, dev)
    return x.detach().cpu().numpy()
