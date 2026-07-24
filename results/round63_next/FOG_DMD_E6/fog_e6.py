# fog_e6.py -- E6: "does measurement-design redundancy kill the bilinear collusion?"
# Self-contained (no import from FOG_DMD_PROBE, which another agent may be writing).
# Conventions COPIED to match fog_common.py / fog_tracker.py exactly so numbers compare.
#
# Model (per medium block tau, per slot j):
#   b[tau,j] = < P[idx[tau,j]] (.) w_tau , x >  + Poisson shot,   w_tau = lognormal, E[w]=1
# Unified block representation for ALL arms:
#   P   : (M,N) pattern bank
#   idx : (T,S) long -- which pattern fires in block tau, slot j  (S = slots per block = 96)
#   W   : (T,N) block-constant medium field   (dev-grade: coherence ~ one block)
#   b   : (T,S) buckets
# Cartesian arms (A0/A1): idx[tau] = arange(M) for every tau  (all patterns every block).
# Slot arms (A2/A3): idx[tau] = a balanced random draw (each pattern equally often overall).
import numpy as np
import torch

DTYPE = torch.float64


def _dev(dev):
    return torch.device(dev if (dev != 'cuda' or torch.cuda.is_available()) else 'cpu')


def _to(t, dev, dtype=DTYPE):
    return torch.as_tensor(np.asarray(t), device=dev, dtype=dtype)


# ============================ scene / basis / medium ============================
def make_scene(n_side=32):
    yy, xx = np.mgrid[0:n_side, 0:n_side] / n_side
    x = (np.exp(-((xx - .3) ** 2 + (yy - .4) ** 2) / .02)
         + .8 * np.exp(-((xx - .7) ** 2 + (yy - .6) ** 2) / .01))
    s = n_side // 32
    x[10 * s:22 * s, 15 * s:18 * s] += 1.0
    x[5 * s:8 * s, 5 * s:26 * s] += .7
    x = (x / x.max()).ravel()
    return x


def dct_basis(d, n_side, exclude_dc=True):
    """Orthonormal 2D DCT-II low-frequency smooth basis, d columns (N x d), DC excluded.
    Columns ordered by radial frequency (matches fog_common.dct_basis)."""
    n = n_side
    k = np.arange(n)[:, None]; i = np.arange(n)[None, :]
    D = np.cos(np.pi * (2 * i + 1) * k / (2 * n))
    D *= np.sqrt(2.0 / n); D[0] *= 1 / np.sqrt(2)
    modes = []
    for p in range(n):
        for q in range(n):
            if exclude_dc and p == 0 and q == 0:
                continue
            modes.append((p * p + q * q, p, q))
    modes.sort()
    cols = []
    for _, p, q in modes[:d]:
        phi = np.outer(D[p], D[q])
        cols.append(phi.ravel())
    U = np.array(cols).T
    Q, _ = np.linalg.qr(U)
    return Q


def lognormal_medium(U, T, sigma_f, tau, rng):
    """Positive mean-normalized lognormal medium, fixed pixelwise RMS sigma_f (indep of d).
    Returns W:(T,N), Z:(T,d), rho, coeff_sd.  (matches fog_common.lognormal_medium)"""
    N, d = U.shape
    coeff_sd = sigma_f * np.sqrt(N / d)
    if tau is None:
        Z = coeff_sd * rng.standard_normal((T, d)); rho = 0.0
    else:
        rho = np.exp(-1.0 / tau)
        Z = np.zeros((T, d)); Z[0] = coeff_sd * rng.standard_normal(d)
        a = np.sqrt(1.0 - rho ** 2) * coeff_sd
        for t in range(1, T):
            Z[t] = rho * Z[t - 1] + a * rng.standard_normal(d)
    G = Z @ U.T
    var_n = (coeff_sd ** 2) * (U ** 2).sum(axis=1)
    W = np.exp(G - 0.5 * var_n[None, :])
    return W, Z, rho, coeff_sd


def tau_schedule(d, tau_coarse=64.0, tau_fine=4.0):
    return tau_coarse * (tau_fine / tau_coarse) ** (np.arange(d) / max(1, d - 1))


def lognormal_medium_mt(U, T, sigma_f, tau_arr, rng):
    """Multi-timescale lognormal medium (A4): each DCT component has its own OU tau_k."""
    N, d = U.shape
    coeff_sd = sigma_f * np.sqrt(N / d)
    rho = np.exp(-1.0 / np.asarray(tau_arr, dtype=float))
    a = np.sqrt(1.0 - rho ** 2) * coeff_sd
    Z = np.zeros((T, d)); Z[0] = coeff_sd * rng.standard_normal(d)
    for t in range(1, T):
        Z[t] = rho * Z[t - 1] + a * rng.standard_normal(d)
    G = Z @ U.T
    var_n = (coeff_sd ** 2) * (U ** 2).sum(axis=1)
    W = np.exp(G - 0.5 * var_n[None, :])
    return W, Z, rho, coeff_sd


# ============================ pattern banks ============================
def make_binary_bank(M, N, rng):
    """A0/A2 random 0/1 physical patterns (matches fog_common 'binary')."""
    return (rng.random((M, N)) < 0.5).astype(float)


def make_fourier_overlap_bank(M, n_side, spacing=3, band=4.0):
    """A1/A3 Fourier-lattice bank with deliberately OVERLAPPING sidebands.
    Physical [0,1] cosine/sine patterns P = 0.5*(1 + cos(2pi(kx*X+ky*Y)/n + phi)).
    Frequencies on a lattice with `spacing`, kept strictly BELOW Nyquist (|kx|,|ky| < n/2)
    to avoid aliasing that would collapse the bank rank; each freq appears as cos AND sin
    (phi in {0, pi/2}) -> M patterns. `band` is the medium spatial half-bandwidth
    (radial DCT-index ~ sqrt(d/(pi/4))/2 ~= 4 cycles for d=48). Overlap of adjacent
    passbands = max(0, 2*band - spacing)/(2*band).  spacing<band -> >50% overlap (>=2 routes
    per scene frequency, the ptychography overlap principle). Returns (P (M,N), meta)."""
    n = n_side
    yy, xx = np.mgrid[0:n, 0:n]
    lim = (n // 2) - 1                                       # strictly below Nyquist (<=15 for n=32)
    freqs = []
    for kx in range(0, lim + 1, spacing):
        for ky in range(-lim, lim + 1, spacing):
            if kx == 0 and ky <= 0:                          # half-plane; exclude DC + axis dups
                continue
            freqs.append((kx * kx + ky * ky, kx, ky))
    freqs.sort()
    nfreq = M // 2
    assert len(freqs) >= nfreq, f"need {nfreq} freqs, only {len(freqs)} below Nyquist at spacing={spacing}"
    chosen = freqs[:nfreq]
    cols, used = [], []
    for _, kx, ky in chosen:
        ph = 2 * np.pi * (kx * xx + ky * yy) / n
        cols.append((0.5 * (1 + np.cos(ph))).ravel())          # cos phase
        cols.append((0.5 * (1 + np.cos(ph - np.pi / 2))).ravel())  # sin phase
        used.append((kx, ky))
    P = np.array(cols[:M])
    overlap = max(0.0, 2 * band - spacing) / (2 * band)
    meta = dict(spacing=spacing, band=band, overlap_frac=float(overlap),
                n_freq=len(used), kmax=float(np.sqrt(chosen[-1][0])),
                rank=int((np.linalg.svd(P, compute_uv=False) > 1e-9).sum()))
    return P, meta


# ============================ metric ============================
def projectors(P):
    Pi = np.linalg.pinv(P)
    rangeP = lambda v: Pi @ (P @ v)
    return Pi, rangeP


def null_metric(x_hat, x_true, rangeP, null_true, nrm0):
    """EXACT metric from the toys: oracle scale gauge, then null-component error.
    null-NMSE = ne**2."""
    x_hat = np.asarray(x_hat, dtype=np.float64)
    s = np.dot(x_hat, x_true) / np.dot(x_hat, x_hat)
    x_hat = s * x_hat
    ne = np.linalg.norm((x_hat - rangeP(x_hat)) - null_true) / nrm0
    te = np.linalg.norm(x_hat - x_true) / np.linalg.norm(x_true)
    return float(ne ** 2), float(te)


def null_fraction(x_true, rangeP):
    nt = x_true - rangeP(x_true)
    return float(np.linalg.norm(nt) ** 2 / np.linalg.norm(x_true) ** 2)


# ============================ schedule ============================
def cartesian_idx(T, M):
    """A0/A1 schedule: every block shows all M patterns in fixed order 0..M-1."""
    return np.broadcast_to(np.arange(M)[None, :], (T, M)).copy()


def slot_idx(T, S, M, rng):
    """A2/A3 schedule: balanced random assignment. Each pattern appears exactly T*S/M
    times overall, order fully randomized, then reshaped to (T,S) blocks.
    Returns idx (T,S) and the flat acquisition-order pattern sequence (T*S,)."""
    reps = (T * S) // M
    assert reps * M == T * S, "T*S must be a multiple of M for a balanced assignment"
    pool = np.repeat(np.arange(M), reps)
    rng.shuffle(pool)
    return pool.reshape(T, S), pool.copy()


# ============================ forward + noise ============================
def forward_clean(P, idx, W, x):
    """pred[tau,j] = < P[idx[tau,j]] (.) W[tau] , x >.  numpy, (T,S)."""
    T, S = idx.shape
    Wx = W * x[None, :]                       # (T,N)
    Pg = P[idx]                               # (T,S,N)
    return np.einsum('tsn,tn->ts', Pg, Wx)


def add_photons(b_clean, photons, rng):
    """Poisson shot at `photons` mean detected counts/bucket. Returns (b_noisy, R_diag).
    photons=None -> clean (tiny R for the Kalman)."""
    if photons is None:
        return b_clean.copy(), np.full_like(b_clean, 1e-8 * (b_clean.mean() ** 2 + 1e-12))
    mb = b_clean.mean()
    scale = photons / max(mb, 1e-12)
    lam = np.clip(b_clean * scale, 0, None)
    counts = rng.poisson(lam).astype(float)
    b_noisy = counts / scale
    R = np.clip(counts, 1.0, None) / (scale ** 2)   # Var(counts)=lam ~ counts, back-scaled
    return b_noisy, R


# ============================ certificate ============================
def spectral_flatness(seq):
    """Wiener-entropy spectral flatness of a real sequence (mean removed).
    exp(mean(log P)) / mean(P) in (0,1]; white->1, tonal/periodic->0."""
    s = np.asarray(seq, float)
    s = s - s.mean()
    F = np.abs(np.fft.rfft(s)) ** 2
    F = F[1:]                                   # drop DC bin
    F = np.clip(F, 1e-30, None)
    return float(np.exp(np.mean(np.log(F))) / np.mean(F))


def walsh_flatness(seq):
    """Walsh(Hadamard)-domain flatness of the carrier (GI_a1 uses the Walsh spectrum).
    Pads to next power of two; same Wiener-entropy statistic on |Walsh coeff|^2."""
    s = np.asarray(seq, float); s = s - s.mean()
    n = 1 << int(np.ceil(np.log2(len(s))))
    v = np.zeros(n); v[:len(s)] = s
    # fast Walsh-Hadamard transform (natural order)
    h = 1
    while h < n:
        for i in range(0, n, h * 2):
            a = v[i:i + h].copy(); b = v[i + h:i + 2 * h].copy()
            v[i:i + h] = a + b; v[i + h:i + 2 * h] = a - b
        h *= 2
    F = v ** 2; F = F[1:]; F = np.clip(F, 1e-30, None)
    return float(np.exp(np.mean(np.log(F))) / np.mean(F))


def carrier_certificate(P, idx_flat, x):
    """Temporal carrier c_n = <P[idx_flat_n], x> laid out in acquisition order; return
    both Fourier and Walsh spectral flatness. High flatness = object 'stopped pretending
    to be time' (GI_a1). Uses the TRUE scene (a schedule property, not an estimate)."""
    c = P[idx_flat] @ x
    return dict(sf_fourier=spectral_flatness(c), sf_walsh=walsh_flatness(c),
                carrier_std=float(c.std()), carrier_mean=float(c.mean()))


# ============================ block-bucket matrix (for cov + spectral init) ============================
def block_bucket_matrix(idx, b, M):
    """Aggregate slot buckets into a (T,M) block x pattern matrix + observation mask.
    Bmat[tau,i] = mean of b over slots in block tau firing pattern i (0 if none);
    mask[tau,i] = 1 if pattern i observed in block tau. Cartesian -> mask all 1."""
    T, S = idx.shape
    Bmat = np.zeros((T, M)); cnt = np.zeros((T, M))
    for tau in range(T):
        np.add.at(Bmat[tau], idx[tau], b[tau])
        np.add.at(cnt[tau], idx[tau], 1.0)
    mask = (cnt > 0).astype(float)
    Bmat = np.where(cnt > 0, Bmat / np.clip(cnt, 1, None), 0.0)
    return Bmat, mask


# ============================ ESTIMATOR (i): staged pathwise ============================
def _solve_x(P, idx, W, b, lam_x):
    """x M-step given medium W: LS over all (T*S) slot equations + Tikhonov, nonneg clamp.
    All torch. P:(M,N) idx:(T,S) long W:(T,N) b:(T,S)."""
    T, S = idx.shape; N = P.shape[1]
    Pg = P[idx]                                  # (T,S,N)
    A = Pg * W[:, None, :]                        # (T,S,N)
    A = A.reshape(T * S, N)
    AtA = A.t() @ A
    Atb = A.t() @ b.reshape(-1)
    ridge = lam_x * torch.trace(AtA) / N
    AtA = AtA + ridge * torch.eye(N, device=P.device, dtype=P.dtype)
    x = torch.linalg.solve(AtA, Atb)
    return torch.clamp(x, min=0.0)


def pathwise_solve(P_np, U_np, idx_np, b_np, coeff_sd, rho, lam_x=1e-4,
                   outer=45, qsteps=25, qlr=6e-3, dev='cuda',
                   z_seed=None, n_starts=5, seed0=0, ou_lambda=0.0, tau_arr=None):
    """Reduced-Q lognormal pathwise EM (matches refine_lognormal_reduced structure,
    generalized to per-slot patterns). Medium z_tau = Zhat[tau] @ Q lives in the centered
    block-bucket SVD temporal subspace (d x d params: stable, unlike full-Z Adam).
    Forward w=exp(Uz - v/2) matches the generator.  Alternates closed-form nonneg x-solve
    with Adam Q-steps (x detached).
      z_seed given -> single warm start from that medium path (TRUE-SEED stability test).
      z_seed None  -> n_starts random data-only starts (BLIND cold-start), pick lowest resid.
    ou_lambda>0 adds the A4 temporal-law penalty (per-component AR(1) mixing rejection).
    Returns dict(per_start_x, objs, x_best)."""
    dev = _dev(dev)
    P = _to(P_np, dev); U = _to(U_np, dev); idx = torch.as_tensor(idx_np, device=dev, dtype=torch.long)
    b = _to(b_np, dev); Ut = U.t()
    M, N = P.shape; d = U.shape[1]; T, S = idx.shape
    var_n = (coeff_sd ** 2) * (U ** 2).sum(dim=1)
    # centered block-bucket SVD -> temporal subspace Zhat (T,d)
    Bmat_np, mask_np = block_bucket_matrix(idx_np, b_np, M)
    Bmat = _to(Bmat_np, dev); mask = _to(mask_np, dev)
    col_mean = (Bmat * mask).sum(0) / torch.clamp(mask.sum(0), min=1)
    Bfill = torch.where(mask > 0, Bmat, col_mean[None, :])
    Bc = Bfill - Bfill.mean(0, keepdim=True)
    _, _, Vh = torch.linalg.svd(Bc.t(), full_matrices=False)
    dk = min(d, Vh.shape[0])                       # temporal-subspace rank (dk<=d if T or M small)
    Zhat = Vh[:dk].t()                             # (T,dk)
    # Stage-A row-space init
    Pinv = torch.linalg.pinv(P)
    xA = torch.clamp(Pinv @ (Bfill.mean(0)), min=0.0)

    rho_v = _to(np.broadcast_to(np.asarray(rho, float), (d,)).copy(), dev) if not np.isscalar(rho) \
        else torch.full((d,), float(rho), device=dev, dtype=P.dtype)

    def run(Q0):
        Q = Q0.clone().detach().requires_grad_(True)
        opt = torch.optim.Adam([Q], lr=qlr)
        x = xA.clone()
        for it in range(outer):
            with torch.no_grad():
                W = torch.exp((Zhat @ Q) @ Ut - 0.5 * var_n[None, :])
                x = _solve_x(P, idx, W, b, lam_x)
            for _ in range(qsteps):
                opt.zero_grad()
                Z = Zhat @ Q
                W = torch.exp(Z @ Ut - 0.5 * var_n[None, :])
                Pg = P[idx]
                pred = torch.einsum('tsn,tn->ts', Pg, W * x[None, :])
                loss = ((pred - b) ** 2).sum()
                if ou_lambda > 0:                 # A4: reject component-timescale mixing
                    innov = Z[1:] - rho_v[None, :] * Z[:-1]
                    q = torch.clamp((1 - rho_v ** 2), min=1e-4) * (coeff_sd ** 2)
                    loss = loss + ou_lambda * (innov ** 2 / q[None, :]).sum()
                loss.backward(); opt.step()
        with torch.no_grad():
            W = torch.exp((Zhat @ Q) @ Ut - 0.5 * var_n[None, :])
            x = _solve_x(P, idx, W, b, lam_x)
            Pg = P[idx]
            pred = torch.einsum('tsn,tn->ts', Pg, W * x[None, :])
            resid = float(((pred - b) ** 2).sum().item())
        return x.detach().cpu().numpy(), resid

    per_start, objs = [], []
    if z_seed is not None:
        # TRUE-SEED: project the true medium path onto Zhat -> Q0
        Zs = _to(z_seed, dev)
        Zs = (Zs - Zs.mean(0, keepdim=True))
        Q0 = torch.linalg.lstsq(Zhat, Zs).solution
        x, r = run(Q0); per_start.append(x); objs.append(r)
    else:
        for s in range(n_starts):
            torch.manual_seed(1000 + seed0 + s)
            Q0 = 0.1 * torch.randn(dk, d, device=dev, dtype=P.dtype)
            x, r = run(Q0); per_start.append(x); objs.append(r)
    ibest = int(np.argmin(objs))
    return dict(per_start_x=per_start, objs=objs, x_best=per_start[ibest])


# ============ ESTIMATOR (i) variant: FULL-Z refinement seeded at truth (true-seed test) ============
def pathwise_fullz_refine(P_np, U_np, idx_np, b_np, coeff_sd, rho, z_init,
                          lam_x=1e-4, outer=18, zsteps=40, zlr=4e-3, dev='cuda',
                          ou_lambda=1.0, tau_arr=None):
    """FULL per-block medium z (T x d), seeded at z_init (=truth for the collusion test),
    then re-estimated. Alternates closed-form nonneg x-solve with Adam z-steps (x detached),
    lognormal forward w=exp(Uz - v/2) + OU prior. This is the faithful 'refinement seeded
    from the true medium, medium re-estimated' probe of the E5 report Sec.11c: STAY(<=0.1)
    means the design killed the collusion; DRIFT(~0.7) means it did not. Returns x (numpy)."""
    dev = _dev(dev)
    P = _to(P_np, dev); U = _to(U_np, dev)
    idx = torch.as_tensor(idx_np, device=dev, dtype=torch.long)
    b = _to(b_np, dev); Ut = U.t()
    M, N = P.shape; d = U.shape[1]; T, S = idx.shape
    var_n = (coeff_sd ** 2) * (U ** 2).sum(dim=1)
    rho_v = _to(np.broadcast_to(np.asarray(rho, float), (d,)).copy(), dev) if not np.isscalar(rho) \
        else torch.full((d,), float(rho), device=dev, dtype=P.dtype)
    q = torch.clamp((1 - rho_v ** 2), min=1e-4) * (coeff_sd ** 2)
    z = _to(z_init, dev).clone()
    Pg = P[idx]
    for it in range(outer):
        with torch.no_grad():
            W = torch.exp(z @ Ut - 0.5 * var_n[None, :])
            x = _solve_x(P, idx, W, b, lam_x)
        z = z.detach().requires_grad_(True)
        opt = torch.optim.Adam([z], lr=zlr)
        for _ in range(zsteps):
            opt.zero_grad()
            W = torch.exp(z @ Ut - 0.5 * var_n[None, :])
            pred = torch.einsum('tsn,tn->ts', Pg, W * x[None, :])
            loss = ((pred - b) ** 2).sum()
            innov = z[1:] - rho_v[None, :] * z[:-1]
            loss = loss + ou_lambda * ((innov ** 2 / q[None, :]).sum() + (z[0] ** 2 / (coeff_sd ** 2)).sum())
            loss.backward(); opt.step()
    with torch.no_grad():
        W = torch.exp(z @ Ut - 0.5 * var_n[None, :])
        x = _solve_x(P, idx, W, b, lam_x)
    return x.detach().cpu().numpy()


# ============================ ESTIMATOR (ii): covariance-domain ============================
def _masked_lagcov(Bmat, mask, l):
    """Masked lag-l covariance across blocks: S[i,j]=<Bc_i(tau) Bc_j(tau+l)> over blocks where
    both observed. Bmat:(T,M) mask:(T,M) torch. Returns (M,M)."""
    T, M = Bmat.shape
    # center each column over observed blocks
    csum = (Bmat * mask).sum(0); cnt = torch.clamp(mask.sum(0), min=1)
    cmean = csum / cnt
    Bc = (Bmat - cmean[None, :]) * mask
    A = Bc[:T - l]; Bb = Bc[l:]; mA = mask[:T - l]; mB = mask[l:]
    num = A.t() @ Bb
    den = torch.clamp(mA.t() @ mB, min=1.0)
    return num / den


def cov_solve(P_np, U_np, idx_np, b_np, coeff_sd, rho, lags=(1, 2, 4, 8, 16, 32),
              n_starts=5, steps=4000, lr=1.5e-2, dev='cuda', seed0=0, rankproj=True):
    """E5a covariance-domain estimator, generalized to the block schedule (masked lag-cov).
    Never estimates z_tau. Target That ~ H_true H_true^T with H(x)=P diag(x) U.
    min_{x>=0} ||H(x)H(x)^T - That||_F^2, Adam multi-start (data-only). Collusion-free by
    construction (medium marginalized). Returns dict(per_start_x, objs, x_best)."""
    dev = _dev(dev)
    P = _to(P_np, dev); U = _to(U_np, dev)
    M, N = P.shape; d = U.shape[1]; T, S = idx_np.shape
    Bmat_np, mask_np = block_bucket_matrix(idx_np, b_np, M)
    Bmat = _to(Bmat_np, dev); mask = _to(mask_np, dev)
    num = torch.zeros((M, M), device=dev, dtype=P.dtype); wsum = 0.0
    for l in lags:
        if l >= T - 4:
            continue
        r_l = float(np.mean(rho) ** l) if not np.isscalar(rho) else float(rho) ** l
        if r_l < 1e-3:
            continue
        Sl = _masked_lagcov(Bmat, mask, l)
        Sl = 0.5 * (Sl + Sl.t())
        w_l = r_l ** 2
        num = num + w_l * (Sl / r_l); wsum += w_l
    Ghat = num / max(wsum, 1e-12)
    That = Ghat / (coeff_sd ** 2)
    if rankproj:
        w, V = torch.linalg.eigh(That); w = torch.clamp(w, min=0); w[:-d] = 0
        That = (V * w) @ V.t()
    Tnorm = (That ** 2).sum() + 1e-12
    Pinv = torch.linalg.pinv(P)
    xr0 = torch.clamp(Pinv @ Bmat.mean(0), min=1e-3)
    per_start, objs = [], []
    for s in range(n_starts):
        torch.manual_seed(seed0 + s)
        if s == 0:
            u = torch.log(torch.clamp(xr0, min=1e-3)).clone()
        else:
            u = torch.log(torch.clamp(xr0 * torch.rand(N, device=dev, dtype=P.dtype) * 2, min=1e-3))
        u.requires_grad_(True)
        opt = torch.optim.Adam([u], lr=lr)
        for it in range(steps):
            opt.zero_grad()
            x = torch.nn.functional.softplus(u)
            H = P @ (U * x[:, None])
            resid = H @ H.t() - That
            loss = (resid ** 2).sum() / Tnorm
            loss.backward(); opt.step()
        x = torch.nn.functional.softplus(u).detach()
        per_start.append(x.cpu().numpy()); objs.append(float(loss.item()))
    ibest = int(np.argmin(objs))
    return dict(per_start_x=per_start, objs=objs, x_best=per_start[ibest])


# ============================ oracle (known medium) ============================
def oracle_solve(P_np, idx_np, W_np, b_np, lam_x=1e-4, dev='cuda'):
    dev = _dev(dev)
    P = _to(P_np, dev); idx = torch.as_tensor(idx_np, device=dev, dtype=torch.long)
    W = _to(W_np, dev); b = _to(b_np, dev)
    x = _solve_x(P, idx, W, b, lam_x)
    return x.detach().cpu().numpy()
