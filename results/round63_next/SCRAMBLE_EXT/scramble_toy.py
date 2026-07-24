#!/usr/bin/env python
# scramble_toy.py -- FULLY SCRAMBLING dynamic-medium extension of the beyond-band sentinel.
#
# Simulates the coherent scrambling channel of SCRAMBLE_EXT/DERIVATION.md exactly:
#   E_{i,t} = B (G_t A_i)   (fully-developed speckle, iid circular-Gaussian T, grain kernel B)
#   S_{i,t}(r) = |E_{i,t}(r)|^2,   b_{i,t} = sum_r S_{i,t}(r) x(r) + Poisson shot.
# Generated via the exact Kronecker identity  Cov[E_i(r),E_j(r')] = C(r-r') O_ij :
#   E = B ( Z O^{1/2} ),  Z ~ CN(0,1)^{N x M},  O = A^T A (code Gram),  B = FFT lowpass (grain).
#
# Validates numerically, trying to FALSIFY the derivation:
#   A. Mean-band collapse (Thm 1): mu_i(r) flat; E[b_i] ~ x_hat(0); mean-null to beyond-band delta.
#   B. Covariance = Q(x) (Eq 3.3): Var[b_i]=O_ii^2 Q; Cov[b_i,b_j]=O_ij^2 Q (rank-1); dQ matches.
#   C. The sentinel duel: H0 vs H1 (beyond-band sub-resolution delta), matched LR on the
#      bank-covariance, d' + ROC at 3 budgets; mean detector at chance.
#
# WRITE-SCOPE: results/round63_next/SCRAMBLE_EXT/ only.  numpy/scipy, CPU.  No git commit.
import json
import os
import time

import numpy as np
import torch
from scipy.fftpack import dct, idct

HERE = os.path.dirname(os.path.abspath(__file__))
RNG = np.random.default_rng(20260724)
DEV = "cuda" if torch.cuda.is_available() else "cpu"

# ----------------------------------------------------------------- config
N_SIDE = 32
N = N_SIDE * N_SIDE
PB = 3                    # pattern/modulator band (DCT fx,fy <= PB): the hardware wall B_p
K_GRAIN = 13              # grain pupil radius (FFT) -> B_grain >> B_p (fine speckle)
M = 36                   # number of band-limited codes
PHOT = 1e4               # mean detected photoelectrons per bucket (FOG-matched)
CHUNK = 128              # banks per vectorised batch


# ----------------------------------------------------------------- DCT helpers (scene/codes)
def _D2(a):
    return dct(dct(a.reshape(N_SIDE, N_SIDE), axis=0, norm="ortho"), axis=1, norm="ortho")


def _I2(A):
    return idct(idct(A, axis=0, norm="ortho"), axis=1, norm="ortho").ravel()


def dct_img(coeffs):
    C = np.zeros((N_SIDE, N_SIDE))
    for (i, j, v) in coeffs:
        C[i, j] = v
    return _I2(C)


# ----------------------------------------------------------------- codes (band-limited, nonneg)
def make_codes(M_, pb=PB, seed=1):
    rp = np.random.default_rng(seed)
    A = np.zeros((N, M_))
    for m in range(M_):
        Cp = np.zeros((N_SIDE, N_SIDE))
        Cp[:pb + 1, :pb + 1] = rp.standard_normal((pb + 1, pb + 1))
        Cp[0, 0] = 0.0
        f = _I2(Cp)
        f /= np.abs(f).max()
        A[:, m] = 0.5 + 0.45 * f          # nonneg DMD code, illumination field amplitude
    return A


# ----------------------------------------------------------------- grain operator B and kernel
def grain_pupil(k_grain=K_GRAIN):
    """FFT low-pass disk of radius k_grain -> B = ifft2(mask*fft2(.)); C = ifft2(mask) (PSF)."""
    ky = np.fft.fftfreq(N_SIDE) * N_SIDE
    KX, KY = np.meshgrid(ky, ky, indexing="ij")
    mask = (np.hypot(KX, KY) <= k_grain).astype(float)
    return mask


def grain_kernel_and_Q(mask):
    """Return g(r)=|C(r)|^2 image (C=PSF=ifft2 mask) and a closure Q(x)=x^T G x, G_{rr'}=g(r-r')."""
    C_psf = np.fft.ifft2(mask)                     # (N_SIDE,N_SIDE) complex PSF, origin at [0,0]
    g = np.abs(C_psf) ** 2                          # |C|^2 kernel (real, nonneg)
    Ghat = np.fft.fft2(g)                           # circulant symbol

    def Q(x_vec):
        X = x_vec.reshape(N_SIDE, N_SIDE)
        gx = np.fft.ifft2(Ghat * np.fft.fft2(X)).real     # circular conv (g (*) x)
        return float(np.sum(X * gx))

    return g, Ghat, C_psf, Q


# ----------------------------------------------------------------- scene + beyond-band deltas
def _beyond_band_project(v_img, pb=PB):
    """Project an N-vector onto DCT frequencies strictly beyond the pattern band (i+j>pb),
    excluding DC. Returns an N-vector (beyond-band, DC-free)."""
    Cc = _D2(v_img)
    mask = np.zeros_like(Cc)
    for i in range(N_SIDE):
        for j in range(N_SIDE):
            if (i + j) > pb:
                mask[i, j] = 1.0
    mask[0, 0] = 0.0
    return _I2(Cc * mask)


def make_scene_and_delta(eps, pb=PB, seed=7):
    """Static scene x plus two sub-resolution, beyond-band, DC-free perturbations:
       - 'ortho'    : high-DCT spikes spectrally DISJOINT from x  (worst case; only the eps^2
                       floor delta^T G delta responds -> the guaranteed-but-weak sentinel).
       - 'coherent' : amplified beyond-band texture OF x itself   (generic physical change;
                       the coherent term 2 x^T G delta ~ eps responds -> the strong sentinel)."""
    rs = np.random.default_rng(seed)
    inb = [(i, j, rs.standard_normal()) for i in range(pb + 1) for j in range(pb + 1)
           if (i + j) > 0]
    base_hi = [(5, 2, 1.4), (2, 6, -1.2)]            # beyond-band scene texture (in BOTH H0,H1)
    x = dct_img(inb + base_hi)
    x = x - x.min(); x = x / x.max()
    x = 0.12 + 0.76 * x                              # headroom in [0.12, 0.88]
    nx = np.linalg.norm(x)
    # ortho: high freqs NOT present in x
    hi_freqs = [(8, 3), (3, 9), (7, 7), (11, 1), (1, 10), (6, 8)]
    d_o = dct_img([(i, j, rs.choice([-1.0, 1.0])) for (i, j) in hi_freqs])
    d_o = d_o / np.linalg.norm(d_o) * (eps * nx)
    # coherent: the scene's own beyond-band content, amplified (correlates with x)
    d_c = _beyond_band_project(x)
    d_c = d_c / np.linalg.norm(d_c) * (eps * nx)
    return x, {"ortho": d_o, "coherent": d_c}, hi_freqs


# ----------------------------------------------------------------- bank generator (Kronecker, GPU)
def gen_banks_bstream(nbanks, A, O_half, mask, x0, x1, phot=PHOT, chunk=512, seed=0,
                      rng=None):
    """Generate nbanks fully-scrambled speckle realisations via the exact Kronecker identity
    E = B (Z O^{1/2}); return paired bucket streams b0 (scene x0), b1 (scene x1), shape
    (nbanks, M) -- SAME speckle for both (perfectly paired H0/H1). torch-CUDA, Poisson shot
    with mean `phot`."""
    if rng is not None and seed == 0:
        seed = int(rng.integers(1 << 31))
    M_ = A.shape[1]
    dev = DEV
    g = torch.Generator(device=dev).manual_seed(int(seed))
    Oh = torch.tensor(O_half, dtype=torch.complex64, device=dev)
    mk = torch.tensor(mask, dtype=torch.complex64, device=dev)
    tx0 = torch.tensor(x0, dtype=torch.float32, device=dev)
    tx1 = torch.tensor(x1, dtype=torch.float32, device=dev)
    b0 = torch.empty((nbanks, M_), dtype=torch.float32, device=dev)
    b1 = torch.empty((nbanks, M_), dtype=torch.float32, device=dev)
    done = 0
    while done < nbanks:
        nb = min(chunk, nbanks - done)
        Z = torch.randn(nb, N, M_, dtype=torch.complex64, device=dev, generator=g)
        W = Z @ Oh                                     # (nb,N,M) code-correlated white field
        W2 = W.reshape(nb, N_SIDE, N_SIDE, M_).permute(0, 3, 1, 2)   # (nb,M,ns,ns)
        E = torch.fft.ifft2(torch.fft.fft2(W2, dim=(-2, -1)) * mk, dim=(-2, -1))  # grain
        S = (E.real ** 2 + E.imag ** 2).reshape(nb, M_, N)           # speckle intensity
        b0[done:done + nb] = torch.einsum("bmn,n->bm", S, tx0)
        b1[done:done + nb] = torch.einsum("bmn,n->bm", S, tx1)
        done += nb
    # Poisson shot: scale so mean bucket = phot photoelectrons; realise per stream
    if phot < 1e17:
        sc0 = phot / torch.clamp(b0.mean(), min=1e-12)
        sc1 = phot / torch.clamp(b1.mean(), min=1e-12)
        b0 = torch.poisson(torch.clamp(b0 * sc0, min=0), generator=g) / sc0
        b1 = torch.poisson(torch.clamp(b1 * sc1, min=0), generator=g) / sc1
    return b0.cpu().numpy().astype(np.float64), b1.cpu().numpy().astype(np.float64)


# ----------------------------------------------------------------- covariance detector
def sample_cov(bstream):
    T = bstream.shape[0]
    D = bstream - bstream.mean(0)
    return (D.T @ D) / (T - 1)


def efficient_weight(OhO, Q0, shot_diag):
    """Fisher-efficient matched-score weight  W = Sigma0^{-1} (dSigma/dQ) Sigma0^{-1},
    dSigma/dQ = OhO, Sigma0 = Q0*OhO + diag(shot). This is the FOG 'W = V0^{-1} dV V0^{-1}'
    detector -- whitening by Sigma0^{-1} de-correlates the M code looks and recovers M_eff ~ M,
    where the naive projection <Cov,OhO> gets only M_eff ~ 1."""
    Sig0 = Q0 * OhO + np.diag(shot_diag)
    Sig0inv = np.linalg.inv(Sig0)
    return Sig0inv @ OhO @ Sig0inv


def score_stat(bstream, W):
    """Efficient matched-score statistic t = <W, Cov(bstream)> (monotone in Q_hat)."""
    return float(np.sum(W * sample_cov(bstream)))


def auc_dprime(s0, s1):
    s0 = np.asarray(s0); s1 = np.asarray(s1)
    # AUC via Mann-Whitney
    allv = np.concatenate([s0, s1])
    order = allv.argsort()
    ranks = np.empty_like(order, float)
    ranks[order] = np.arange(1, len(allv) + 1)
    r1 = ranks[len(s0):].sum()
    auc = (r1 - len(s1) * (len(s1) + 1) / 2) / (len(s0) * len(s1))
    m0, m1 = s0.mean(), s1.mean()
    v = 0.5 * (s0.var(ddof=1) + s1.var(ddof=1))
    dp = (m1 - m0) / np.sqrt(v + 1e-30)
    return float(auc), float(dp)


def roc_points(s0, s1, fprs=(0.01, 0.05, 0.10)):
    s0 = np.sort(np.asarray(s0)); s1 = np.asarray(s1)
    out = {}
    for f in fprs:
        thr = np.quantile(s0, 1 - f)
        out[f] = float((s1 > thr).mean())
    return out


# =================================================================== main
def main():
    t0 = time.time()
    _sm = os.environ.get("SMOKE", "0") == "1"
    rep = {"config": dict(N_SIDE=N_SIDE, N=N, PB=PB, K_GRAIN=K_GRAIN, M=M, PHOT=PHOT, smoke=_sm)}

    A = make_codes(M)
    O = A.T @ A
    w, V = np.linalg.eigh(O)
    O_half = (V * np.sqrt(np.clip(w, 0, None))) @ V.T
    OhO = O * O                                       # Hadamard square (known)
    mask = grain_pupil()
    g, Ghat, C_psf, Q = grain_kernel_and_Q(mask)
    C0 = C_psf[0, 0].real                             # C(0) = mask.sum()/N

    # ---- aperture facts
    kg_support = int((np.abs(Ghat) > 1e-9 * np.abs(Ghat).max()).sum())
    rep["aperture"] = dict(
        k_pattern=PB, k_grain=K_GRAIN,
        Ghat_nonzero_bins=kg_support, Ghat_total_bins=N,
        grain_band_ge_pattern=bool(K_GRAIN > PB),
        C0=float(C0),
    )

    # ---- VALIDATION A: mean flatness + mean-null (needs empirical mean speckle envelope)
    eps_main = 0.05
    x, deltas, hi_freqs = make_scene_and_delta(eps_main)
    d = deltas["ortho"]                                # use the strict beyond-band worst case here
    xd = x + d
    rep["scene"] = dict(eps=eps_main, min_x=float(x.min()), min_xd=float(xd.min()),
                        dc_delta_ortho=float(deltas["ortho"].sum()),
                        dc_delta_coherent=float(deltas["coherent"].sum()),
                        norm_ratio=float(np.linalg.norm(d) / np.linalg.norm(x)),
                        hi_freqs=hi_freqs)

    # empirical mean speckle intensity envelope for a few codes (many epochs)
    n_env = 512 if _sm else 4000
    codes_probe = list(range(min(4, M)))
    g_env = torch.Generator(device=DEV).manual_seed(101)
    Oh = torch.tensor(O_half, dtype=torch.complex64, device=DEV)
    mk = torch.tensor(mask, dtype=torch.complex64, device=DEV)
    acc = torch.zeros((len(codes_probe), N), dtype=torch.float32, device=DEV)
    done = 0
    while done < n_env:
        nb = min(512, n_env - done)
        Z = torch.randn(nb, N, M, dtype=torch.complex64, device=DEV, generator=g_env)
        W = Z @ Oh
        W2 = W.reshape(nb, N_SIDE, N_SIDE, M).permute(0, 3, 1, 2)
        E = torch.fft.ifft2(torch.fft.fft2(W2, dim=(-2, -1)) * mk, dim=(-2, -1))
        S = (E.real ** 2 + E.imag ** 2).reshape(nb, M, N)
        acc += S[:, codes_probe, :].sum(0)
        done += nb
    env = (acc / n_env).cpu().numpy().astype(np.float64)   # (4, N) empirical mu_i(r)
    # spatial contrast (CV) of the mean envelope, and its non-DC spectral fraction
    cv = (env.std(1) / env.mean(1))
    env_spec_nonDC = []
    for k in range(env.shape[0]):
        E2 = _D2(env[k])
        tot = (E2 ** 2).sum()
        nondc = tot - E2[0, 0] ** 2
        env_spec_nonDC.append(float(nondc / tot))
    # theory: CV -> sqrt(Var of estimator)/mean ~ contrast/sqrt(n_env); flat means CV small
    rep["valA_mean_flatness"] = dict(
        codes=codes_probe,
        envelope_CV=[float(c) for c in cv],
        envelope_nonDC_power_frac=env_spec_nonDC,
        note="CV ~ 1/sqrt(n_env) if truly flat; nonDC frac ~ residual MC noise",
        est_MC_CV=float(1.0 / np.sqrt(n_env)),
    )
    # analytic mean-null: Delta m_i = C0 * O_ii * sum(delta) = 0 for DC-free delta
    Oii = np.diag(O)
    dm = C0 * Oii * d.sum()
    m_x = C0 * Oii * x.sum()
    rep["valA_mean_null"] = dict(
        delta_mean_rel=float(np.linalg.norm(dm) / np.linalg.norm(m_x)),
        dc_of_delta=float(d.sum()),
        note="mean channel sees only x_hat(0); DC-free beyond-band delta -> exact null",
    )

    # ---- VALIDATION B: Var[b_i] = O_ii^2 Q(x); rank-1 covariance; dQ matches
    Qx = Q(x)
    dQ_both = {}
    for name, dd in deltas.items():
        Gd = _apply_G(dd, Ghat)
        cross = float(2 * (x @ Gd)); quad = float(dd @ Gd)
        dQ_both[name] = dict(cross=cross, quad=quad, dQ_pred=cross + quad,
                             dQ_meas=float(Q(x + dd) - Qx),
                             dQ_over_Q=float((cross + quad) / Qx))
    # keep the ortho one for the headline printout / back-compat fields
    Qxd = Q(xd)
    dQ_pred = dQ_both["ortho"]["dQ_pred"]
    dQ_terms = dict(cross=dQ_both["ortho"]["cross"], quad=dQ_both["ortho"]["quad"])
    # empirical covariance over many banks (no shot for the clean structural check)
    n_val = 800 if _sm else 6000
    b0v, _ = gen_banks_bstream(n_val, A, O_half, mask, x, x, phot=1e18,
                               rng=np.random.default_rng(202))  # huge phot ~ shot-free
    mu = b0v.mean(0); Dv = b0v - mu
    Covv = (Dv.T @ Dv) / (n_val - 1)
    var_emp = np.diag(Covv)
    Qhat_per_code = var_emp / (Oii ** 2)
    rep["valB_cov_is_Q"] = dict(
        Q_x_theory=float(Qx),
        Qhat_from_var_median=float(np.median(Qhat_per_code)),
        Qhat_from_var_iqr=[float(np.percentile(Qhat_per_code, 25)), float(np.percentile(Qhat_per_code, 75))],
        var_over_O2Q_ratio_median=float(np.median(var_emp / (Oii ** 2 * Qx))),
        note="empirical Var[b_i]/O_ii^2 should equal Q(x) across all codes",
    )
    # rank-1 structure: fraction of covariance frobenius energy along the OhO direction
    proj = np.sum(Covv * OhO) / np.sum(OhO * OhO)
    resid = Covv - proj * OhO
    rank1_frac = 1.0 - (np.linalg.norm(resid) ** 2) / (np.linalg.norm(Covv) ** 2)
    # eigenvalue spectrum of the sample covariance (leading eigenvalue dominance)
    evals = np.linalg.eigvalsh(Covv)[::-1]
    rep["valB_rank1"] = dict(
        frac_cov_energy_in_OhO_direction=float(rank1_frac),
        top_eig_frac=float(evals[0] / evals.sum()),
        top3_eig_frac=float(evals[:3].sum() / evals.sum()),
        note="fully-scrambled cov collapses to Q*(O.O): rank-1 in matrix space",
    )
    rep["valB_dQ"] = dict(
        dQ_predicted=float(dQ_pred), dQ_measured=float(Qxd - Qx),
        rel_err=float(abs(dQ_pred - (Qxd - Qx)) / abs(Qxd - Qx + 1e-30)),
        terms=dQ_terms, per_delta=dQ_both,
        quad_term_positive=bool(dQ_terms["quad"] > 0),
        note="dQ = 2 xGd + dGd; quad dGd>0 for ANY delta (no blind direction). "
             "ortho: coherent term ~0 -> dQ/Q~eps^2 (worst case). "
             "coherent: 2xGd~eps -> dQ/Q~eps (strong sentinel).",
    )

    # ---- VALIDATION C: the sentinel duel (H0 vs H1), efficient matched-score LR on bank covariance
    # Two delta types, budgets matched to each regime; d'(T)=a*sqrt(T) fit -> T_req(d'=5).
    if _sm:
        duel_cfg = {"coherent": dict(budgets=[64, 256, 1024], eps=[0.05]),
                    "ortho":    dict(budgets=[256, 1024], eps=[0.05])}
        R = 20
    else:
        duel_cfg = {"coherent": dict(budgets=[512, 2048, 8192], eps=[0.05, 0.02]),
                    "ortho":    dict(budgets=[8192, 32768, 65536], eps=[0.05])}
        R = 100                                        # records per hypothesis
    detection = {}
    raw_scores = {}                                    # (dtype,eps,B) -> (Q0 list, Q1 list) for fig
    for dtype, cfg in duel_cfg.items():
        budgets = cfg["budgets"]
        detection[dtype] = {}
        for eps in cfg["eps"]:
            xe, dd, _ = make_scene_and_delta(eps)
            de = dd[dtype]
            xde = xe + de
            Eb = C0 * Oii * xe.sum()
            scale = PHOT / max(Eb.mean(), 1e-12)
            shot_diag = Eb / scale
            W = efficient_weight(OhO, Q(xe), shot_diag)   # Fisher-efficient matched-score weight
            Q0 = {b: [] for b in budgets}
            Q1 = {b: [] for b in budgets}
            for rec in range(R):
                _doff = {"coherent": 0, "ortho": 100000}[dtype]
                rng_r = np.random.default_rng(9000 + 31 * rec + _doff + int(1000 * eps))
                b0, b1 = gen_banks_bstream(max(budgets), A, O_half, mask, xe, xde,
                                           phot=PHOT, rng=rng_r)
                for B_ in budgets:
                    Q0[B_].append(score_stat(b0[:B_], W))
                    Q1[B_].append(score_stat(b1[:B_], W))
            res = {}
            dps = []
            for B_ in budgets:
                auc, dp = auc_dprime(Q0[B_], Q1[B_])
                roc = roc_points(Q0[B_], Q1[B_])
                res[str(B_)] = dict(T_eff=B_, d_prime=dp, auc=auc, tpr_at_fpr=roc,
                                    Q0_mean=float(np.mean(Q0[B_])), Q1_mean=float(np.mean(Q1[B_])))
                dps.append((B_, dp))
                raw_scores[(dtype, eps, B_)] = (list(Q0[B_]), list(Q1[B_]))
            # fit d' = a*sqrt(T) (through the origin) -> a, and extrapolate T_req(d'=5)
            Ts = np.array([t for t, _ in dps], float)
            Dp = np.array([max(d_, 1e-6) for _, d_ in dps])
            a = float(np.sum(np.sqrt(Ts) * Dp) / np.sum(Ts))       # LS slope of d' vs sqrt(T)
            res["fit"] = dict(a_dprime_per_sqrtT=a,
                              T_req_dprime5=float((5.0 / max(a, 1e-9)) ** 2),
                              dQ_over_Q=float(dQ_both[dtype]["dQ_over_Q"]
                                              if abs(eps - eps_main) < 1e-9 else np.nan))
            detection[dtype][f"eps_{eps}"] = res

    # ---- mean-channel detector (should be at chance for DC-free beyond-band delta, BOTH types)
    xe, dd05, _ = make_scene_and_delta(0.05)
    Bmean = 512
    tmpl = C0 * Oii * xe.sum()                          # H0 mean template
    rep["valC_mean_detector"] = {}
    for dtype in ["ortho", "coherent"]:
        xde = xe + dd05[dtype]
        ms0, ms1 = [], []
        for rec in range(80):
            rng_r = np.random.default_rng(4400 + rec + (0 if dtype == "ortho" else 500))
            b0, b1 = gen_banks_bstream(Bmean, A, O_half, mask, xe, xde, phot=PHOT, rng=rng_r)
            ms0.append(float(np.sum((b0.mean(0) - tmpl) ** 2)))
            ms1.append(float(np.sum((b1.mean(0) - tmpl) ** 2)))
        auc_m, dp_m = auc_dprime(ms0, ms1)
        rep["valC_mean_detector"][dtype] = dict(budget=Bmean, auc=auc_m, d_prime=dp_m)
    rep["valC_mean_detector"]["note"] = ("beyond-band DC-free delta invisible to mean channel "
                                         "-> AUC~0.5 for BOTH delta types (the wall)")
    auc_m = rep["valC_mean_detector"]["ortho"]["auc"]
    dp_m = rep["valC_mean_detector"]["ortho"]["d_prime"]

    rep["detection"] = detection
    rep["runtime_sec"] = round(time.time() - t0, 1)

    with open(os.path.join(HERE, "SCRAMBLE_RESULTS.json"), "w") as f:
        json.dump(rep, f, indent=2)

    # ---- summary figure: d' vs sqrt(T) with fits + ROC curves (coherent) + envelope demo
    try:
        _make_figure(rep, raw_scores, env, C_psf, Ghat)
    except Exception as e:  # figure is a nicety; never fail the run on it
        print("  (figure skipped: %s)" % e)

    # ---- console summary
    print(f"[done in {rep['runtime_sec']}s]")
    print("APERTURE: k_pattern=%d  k_grain=%d  Ghat nonzero bins=%d/%d (grain band >> pattern band)"
          % (PB, K_GRAIN, kg_support, N))
    print("VAL A mean flatness: envelope CV=%s (MC floor %.4f); nonDC power frac=%s"
          % (["%.4f" % c for c in cv], 1 / np.sqrt(n_env), ["%.2e" % e for e in env_spec_nonDC]))
    print("VAL A mean-null: |dm|/|m| = %.2e (DC of delta = %.2e)"
          % (rep["valA_mean_null"]["delta_mean_rel"], d.sum()))
    print("VAL B cov=Q: median Var/(O^2 Q) = %.3f ; rank-1 energy frac = %.4f ; top-eig frac = %.3f"
          % (rep["valB_cov_is_Q"]["var_over_O2Q_ratio_median"],
             rank1_frac, evals[0] / evals.sum()))
    print("VAL B dQ: predicted %.4e  measured %.4e  (quad term dGd=%.3e > 0)"
          % (dQ_pred, Qxd - Qx, dQ_terms["quad"]))
    print("VAL C mean detector AUC = %.3f  d'=%.2f  (expect ~0.5 / ~0, both delta types)" % (auc_m, dp_m))
    for dtype, cfg in duel_cfg.items():
        for eps in cfg["eps"]:
            block = detection[dtype][f"eps_{eps}"]
            print("  [%s] dQ/Q(eps=.05)=%.2e  T_req(d'=5)=%.3g banks (=%.1f s @12.8ms)" %
                  (dtype, dQ_both[dtype]["dQ_over_Q"], block["fit"]["T_req_dprime5"],
                   block["fit"]["T_req_dprime5"] * 0.0128))
            for B_ in cfg["budgets"]:
                r = block[str(B_)]
                print("     DUEL eps=%.2f T_eff=%6d : d'=%.2f  AUC=%.3f  TPR@FPR5%%=%.2f"
                      % (eps, B_, r["d_prime"], r["auc"], r["tpr_at_fpr"][0.05]))


def _apply_G(v, Ghat):
    X = v.reshape(N_SIDE, N_SIDE)
    return np.fft.ifft2(Ghat * np.fft.fft2(X)).real.ravel()


def _roc_curve(s0, s1, npts=200):
    s0 = np.asarray(s0); s1 = np.asarray(s1)
    thr = np.linspace(min(s0.min(), s1.min()), max(s0.max(), s1.max()), npts)
    fpr = np.array([(s0 > t).mean() for t in thr])
    tpr = np.array([(s1 > t).mean() for t in thr])
    return fpr, tpr


def _make_figure(rep, raw_scores, env, C_psf, Ghat):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.1))

    # panel 0: mean speckle envelope (flat) + grain PSF inset info
    im = ax[0].imshow(env[0].reshape(N_SIDE, N_SIDE), cmap="magma")
    ax[0].set_title("mean speckle envelope $\\mu_i(r)$ (flat)\nCV=%.3f (MC floor %.3f)"
                    % (env[0].std() / env[0].mean(), 1 / np.sqrt(4000)))
    ax[0].set_xticks([]); ax[0].set_yticks([])
    plt.colorbar(im, ax=ax[0], fraction=0.046)

    # panel 1: d' vs sqrt(T_eff) with sqrt-law fits
    colors = {"coherent": "#1b7837", "ortho": "#762a83"}
    for dtype in rep["detection"]:
        for ek, block in rep["detection"][dtype].items():
            eps = float(ek.split("_")[1])
            Ts = [block[k]["T_eff"] for k in block if k != "fit"]
            Dp = [block[k]["d_prime"] for k in block if k != "fit"]
            Ts, Dp = np.array(sorted(Ts)), np.array([d for _, d in sorted(zip(Ts, Dp))])
            a = block["fit"]["a_dprime_per_sqrtT"]
            lbl = "%s $\\epsilon$=%.2f (T$_{req}$=%.1e)" % (dtype, eps, block["fit"]["T_req_dprime5"])
            ax[1].plot(np.sqrt(Ts), Dp, "o", color=colors[dtype],
                       alpha=0.5 if eps < 0.05 else 1.0)
            tt = np.linspace(0, np.sqrt(max(Ts)), 50)
            ax[1].plot(tt, a * tt, "-", color=colors[dtype],
                       alpha=0.5 if eps < 0.05 else 1.0, label=lbl)
    ax[1].axhline(5, ls=":", color="k", lw=0.8)
    ax[1].set_xlabel("$\\sqrt{T_{eff}}$  (independent banks)"); ax[1].set_ylabel("detector d'")
    ax[1].set_title("sentinel duel: d' $\\propto \\sqrt{T_{eff}}$\n(covariance channel; mean channel d'$\\approx$0)")
    ax[1].legend(fontsize=6.5, loc="upper left")

    # panel 2: ROC curves, coherent eps=0.05, all budgets
    for (dtype, eps, B_), (q0, q1) in raw_scores.items():
        if dtype != "coherent" or abs(eps - 0.05) > 1e-9:
            continue
        fpr, tpr = _roc_curve(q0, q1)
        auc = rep["detection"]["coherent"]["eps_0.05"][str(B_)]["auc"]
        ax[2].plot(fpr, tpr, label="T$_{eff}$=%d (AUC %.2f)" % (B_, auc))
    ax[2].plot([0, 1], [0, 1], "k:", lw=0.8)
    ax[2].set_xlabel("false-positive rate"); ax[2].set_ylabel("true-positive rate")
    ax[2].set_title("ROC — coherent beyond-band $\\delta$, $\\epsilon$=5%")
    ax[2].legend(fontsize=7, loc="lower right")

    fig.suptitle("Beyond-band sentinel through a FULLY SCRAMBLING medium "
                 "(mean-wall + covariance $Q(x)$ response)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(HERE, "SCRAMBLE_DUEL.png"), dpi=115)
    plt.close(fig)


if __name__ == "__main__":
    main()
