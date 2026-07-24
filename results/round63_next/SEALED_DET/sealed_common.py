#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
SEALED DETECTION PROBE -- shared frozen physical ledger + Fisher/score machinery (R41 sec 4.2-4.5).

ONE module holding every load-bearing physics choice so no gate script can drift:
  * frozen physical ledger (R41 sec 4.2): 64x64, k_p=5, M=128 signed effective codes realized as
    256 complementary nonneg exposures/bank, 20 kHz -> 12.8 ms/bank, T_eff cap 4096 = 52.43 s,
    1e4 detected photoelectrons per PHYSICAL bucket via the |P|x throughput shot term (CORRECTION_NOTE);
  * corrected shot everywhere: R_shot = diag(|P|x * mean(|P|x)/PHOT)  (NEVER the signed mean Px);
  * corrected-shot profiled beyond-band covariance Fisher J_B (two independent engines for D0);
  * lag-augmented covariance model (lag-0 V0=C+R, lag-1 V1=rho*C) so the medium-correlation (tau/lag)
    tangent is a genuine 4th direction distinct from the amplitude direction;
  * the four attribution tangent spaces (mean/in-band, cov/beyond-band, cov-amplitude, cov-lag)
    with efficient (nuisance-profiled) scores and the whitened Gram geometry (Rank-2 simplex);
  * profiled matched filters (make_W beyond-band, W_amp amplitude, W_lag lag) for the arms;
  * the 81-cell analytic grid (3 sigma_f x 3 spectrum x 3 k_w/k_p x 3 claim) and the precommitted
    27-cell 3^(4-1) orthogonal-array MC subset (claim = (sf+spec+kw) mod 3);
  * physical complementary-exposure Poisson MC engine (D0 shot-vs-ledger and end-to-end validation).

Reuses the VALIDATED dev engine (FOG_DMD_PROBE64/detection_roc.py, bold_push.py): same band_modes,
medium_modes, power_spectrum, signed_codes, corrected-shot setup, and Gaussian bank generator, so
every number here is directly comparable to the corrected DETECTION_MAP/ROC evidence chain.

WRITE-SCOPE: read-only config/engine; gate scripts write only under SEALED_DET/.
"""
import hashlib
import numpy as np
import torch

DEV = "cuda" if torch.cuda.is_available() else "cpu"
DT = torch.float64

# ------------------------------------------------------------------ frozen ledger (R41 sec 4.2)
n = 64
N = n * n                                  # 4096
KP = 5                                     # modulator (pattern) band cutoff, Chebyshev |k|<=k_p
M = 128                                    # signed effective codes
N_COMPLEMENTARY = 256                      # 256 complementary nonneg exposures / bank (2 per signed code)
PATTERN_RATE_HZ = 20_000.0                 # 20 kHz
MS_PER_BANK = N_COMPLEMENTARY / PATTERN_RATE_HZ * 1000.0   # 12.8 ms/bank
T_CAP = 4096                               # bank cap = 52.43 s
T_CAP_SEC = T_CAP * MS_PER_BANK / 1000.0   # 52.43 s
PHOT = 1e4                                 # detected photoelectrons / PHYSICAL bucket
TAU = 8.0                                  # nominal medium correlation time (banks); rho = exp(-1/tau)
RHO = float(np.exp(-1.0 / TAU))
BANKS_PER_HOUR = 3600.0 / (MS_PER_BANK / 1000.0)
BANKS_PER_DAY = 24.0 * BANKS_PER_HOUR

# axis levels (R41 sec 4.2)
SIGMA_F_LEVELS = [0.3, 0.6, 1.0]           # nominal medium contrasts (realized reported via sig_eff)
SPECTRA = ["flat", "k^-1", "k^-2"]
KW_OVER_KP = [1, 2, 4]
CLAIMS = [1.25, 1.5, 1.8]                  # claim shells (x k_p)
CLAIM_KHI = {1.25: 6, 1.5: 7, 1.8: 9}      # beyond-band annulus upper Chebyshev freq per claim
EPS_LEVELS = [0.005, 0.01, 0.02, 0.05]     # anomaly magnitudes ||delta||/||x||

# realized (per-pixel RMS) lognormal-OU field contrast for each nominal sigma_f (matches bold_push)
SIG_EFF = {0.3: 0.298, 0.6: 0.503, 1.0: 0.696}

_X, _Y = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")


# ------------------------------------------------------------------ Fourier subspaces
def band_modes(k_lo, k_hi):
    """Orthonormal REAL Fourier modes with Chebyshev freq max(|kx|,|ky|) in [k_lo,k_hi]."""
    cols = []
    for kx in range(-k_hi, k_hi + 1):
        for ky in range(-k_hi, k_hi + 1):
            if k_lo <= max(abs(kx), abs(ky)) <= k_hi:
                ph = 2 * np.pi * (kx * _X + ky * _Y) / n
                cols.append(np.cos(ph).ravel())
                cols.append(np.sin(ph).ravel())
    Q, R = np.linalg.qr(np.array(cols).T)
    return Q[:, np.abs(np.diag(R)) > 1e-8]


def medium_modes(k_w):
    """Zero-mean medium Fourier modes up to Chebyshev freq k_w (DC excluded), with radial |k|."""
    cols, kabs = [], []
    for kx in range(-k_w, k_w + 1):
        for ky in range(-k_w, k_w + 1):
            if max(abs(kx), abs(ky)) == 0:
                continue
            if (kx < 0) or (kx == 0 and ky < 0):
                continue
            ph = 2 * np.pi * (kx * _X + ky * _Y) / n
            for f in (np.cos(ph).ravel(), np.sin(ph).ravel()):
                nr = np.linalg.norm(f)
                if nr > 1e-8:
                    cols.append(f / nr)
                    kabs.append(np.hypot(kx, ky))
    return np.array(cols).T, np.array(kabs)


def power_spectrum(kabs, shape, se):
    """Medium power spectrum with total field variance N*se^2 (se = realized per-pixel RMS)."""
    w = {"flat": np.ones_like(kabs), "k^-1": 1.0 / kabs, "k^-2": 1.0 / kabs ** 2}[shape]
    return w / w.sum() * (N * se ** 2)


# in-band scene nuisance modes (includes DC) and per-claim beyond-band annulus
U_in_np = band_modes(0, KP)
d_eta = U_in_np.shape[1]
BETA_np = {c: band_modes(KP + 1, CLAIM_KHI[c]) for c in CLAIMS}


def signed_codes(seed=10):
    """M signed band-limited random codes (physical complementary nonneg pairs; the signed effective
    row enters the Fisher). Band-limited to 1<=|k|<=k_p (DC excluded -> zero-mean signed)."""
    rp = np.random.default_rng(seed)
    Uc = band_modes(1, KP)
    P = rp.standard_normal((M, Uc.shape[1])) @ Uc.T
    P /= np.abs(P).max(axis=1, keepdims=True)
    return P


# the frozen fixed sealed code bank (seed 10) -- FIXED-COV / FIXED-MEAN arms reuse this
P_FIXED = signed_codes(10)


# ------------------------------------------------------------------ scenes (k_p=5 convention)
def witness_scene(seed=4):
    """In-band random content + deliberate beyond-band structure -> nonneg image in [0,1]."""
    rs = np.random.default_rng(seed)
    a_in = U_in_np @ rs.standard_normal(U_in_np.shape[1])
    Ube = band_modes(KP + 1, 9)
    a_be = Ube @ (2.0 * rs.choice([-1, 1], Ube.shape[1]))
    x = a_in + a_be
    x -= x.min()
    return x / x.max()


def natural_scene(name):
    """Deterministic skimage natural image, downsampled to 64x64, normalized to [0,1]."""
    from skimage import data, transform
    src = getattr(data, name)().astype(np.float64)
    if src.ndim == 3:
        src = src.mean(axis=2)
    img = transform.resize(src, (n, n), anti_aliasing=True, mode="reflect")
    return ((img - img.min()) / (np.ptp(img) + 1e-12)).ravel()


# canonical witness used to BUILD per-cell filters/Fisher (scene-independent geometry object);
# endpoints are re-evaluated on the sealed-bank scenes.
XW_np = witness_scene(4)


# ------------------------------------------------------------------ corrected-shot cell setup
def setup_cell(sf, shape, kwf, claim, x_np=None, code_P=None, want_mc=True, declared=None):
    """Corrected-shot lag-augmented covariance geometry for ONE physical cell.

    Returns V0=C+R (lag-0), V1=rho*C (lag-1 medium), Vinv, per-mode dV stacks, the profiled
    beyond-band Fisher J_B (per effective bank), amplitude/lag score filters, and the four
    attribution tangent derivative blocks (mean/in-band, cov/beyond-band, cov-amplitude, cov-lag).

    `declared` optionally overrides the medium (Um, S, rho) used to BUILD the detector filters
    while the DATA is generated from (sf,shape,kwf); this is the law-mismatch / cross-fit path.
    """
    x = XW_np if x_np is None else np.asarray(x_np, float)
    xt = torch.tensor(x, device=DEV, dtype=DT)
    P_np = P_FIXED if code_P is None else code_P
    Pt = torch.tensor(P_np, device=DEV, dtype=DT)
    Pabs = torch.abs(Pt)

    # --- true generating medium (for DATA) ---
    Um_np, kabs = medium_modes(kwf * KP)
    Um = torch.tensor(Um_np, device=DEV, dtype=DT)
    S = torch.tensor(power_spectrum(kabs, shape, SIG_EFF[sf]), device=DEV, dtype=DT)
    rho = RHO

    # --- declared medium (for FILTERS); defaults to the true medium ---
    if declared is None:
        Um_d, S_d, rho_d = Um, S, rho
    else:
        Umd_np, kabs_d = declared.get("Um", (Um_np, kabs))
        Um_d = torch.tensor(Umd_np, device=DEV, dtype=DT)
        S_d = torch.tensor(power_spectrum(kabs_d, declared.get("shape", shape),
                                          SIG_EFF[declared.get("sf", sf)]), device=DEV, dtype=DT)
        rho_d = float(declared.get("rho", rho))

    def build(Um_, S_, rho_):
        A = Pt * xt[None, :]                              # P diag(x)   (M,N)
        KA = Um_ @ (S_[:, None] * (Um_.t() @ A.t()))      # Kw A^T      (N,M)
        C = A @ KA                                        # medium bucket cov (M,M)
        flux = Pabs @ xt                                  # |P| x  physical throughput
        R = flux * (flux.mean() / PHOT)                   # CORRECTED shot diagonal
        V0 = C + torch.diag(R)
        V1 = rho_ * C                                     # lag-1 (medium correlates; shot white)
        Vinv = torch.linalg.inv(V0 + 1e-9 * torch.eye(M, device=DEV, dtype=DT))
        return A, KA, C, R, V0, V1, Vinv

    A, KA, C, R, V0, V1, Vinv = build(Um, S, rho)                 # DATA-side (true)
    Ad, KAd, Cd, Rd, V0d, V1d, Vinvd = build(Um_d, S_d, rho_d)   # FILTER-side (declared)

    Phi_b = BETA_np[claim]
    Phi_b_t = torch.tensor(Phi_b, device=DEV, dtype=DT)
    db = Phi_b.shape[1]
    Phi_eta = torch.tensor(U_in_np, device=DEV, dtype=DT)

    def dV0_stack(Phi, KA_, C_ref):
        """dV0/dcoeff for each column of Phi: d/dphi of C = P diag(phi) Kw diag(x) P^T + transpose."""
        d = Phi.shape[1]
        out = torch.empty((d, M, M), device=DEV, dtype=DT)
        for k in range(d):
            Bk = (Pt * Phi[:, k][None, :]) @ KA_
            out[k] = Bk + Bk.t()
        return out

    # derivative stacks on the DECLARED medium (efficient score / Fisher use declared law)
    Vk_b = dV0_stack(Phi_b_t, KAd, Cd)          # beyond-band scene (lag-0)
    Vk_in = dV0_stack(Phi_eta, KAd, Cd)         # in-band scene (lag-0)
    Vk_amp = Cd.unsqueeze(0)                     # amplitude: dV0/dlogamp = C
    # lag tangent: d V1 / d log tau = (d rho / d log tau) C = rho*(1/tau) C ... acts on the LAG-1 block
    drho_dlogtau = rho_d * (1.0 / TAU)          # rho=exp(-1/tau); d rho/d log tau = rho/tau
    Vk_lag0 = torch.zeros(1, M, M, device=DEV, dtype=DT)   # lag-0 insensitive to tau (amplitude fixed)
    Vk_lag1 = (drho_dlogtau * Cd).unsqueeze(0)             # lag-1 amplitude changes with tau

    # ---- profiled beyond-band Fisher (per effective bank), ENGINE A: einsum-trace + Schur ----
    Vk_eta = torch.cat([Vk_in, Vk_amp], 0)      # nuisance = in-band scene + medium amplitude
    Gb = torch.einsum("ij,ajk->aik", Vinvd, Vk_b)
    Ge = torch.einsum("ij,ajk->aik", Vinvd, Vk_eta)
    Ibb = 0.5 * torch.einsum("aij,bji->ab", Gb, Gb)
    Ibe = 0.5 * torch.einsum("aij,bji->ab", Gb, Ge)
    Iee = 0.5 * torch.einsum("aij,bji->ab", Ge, Ge)
    Min = Pt @ Phi_eta                          # dm/d(in-band) mean channel (per bank)
    Iee[:d_eta, :d_eta] += (Min.t() @ Vinvd @ Min)
    JB = Ibb - Ibe @ torch.linalg.pinv(Iee, rcond=1e-10) @ Ibe.t()
    JB = 0.5 * (JB + JB.t())
    lamJ = torch.clamp(torch.linalg.eigvalsh(JB), min=0)

    def make_W(c_np):
        """Profiled matched filter for a beyond-band delta with coeffs c (nuisance profiled out)."""
        c = torch.tensor(np.asarray(c_np, float), device=DEV, dtype=DT)
        a_eta = -torch.linalg.solve(Iee + 1e-12 * torch.eye(Iee.shape[0], device=DEV, dtype=DT),
                                    Ibe.t() @ c)
        dV = torch.einsum("a,aij->ij", c, Vk_b) + torch.einsum("a,aij->ij", a_eta, Vk_eta)
        return Vinvd @ dV @ Vinvd

    W_amp = Vinvd @ Cd @ Vinvd                   # medium-amplitude matched filter (lag-0)
    # lag score: applied to the lag-1 sample cross-covariance; matched to Cd direction
    W_lag = Vinvd @ Cd @ Vinvd                   # same spatial shape, but read on the LAG-1 statistic

    out = dict(
        sf=sf, shape=shape, kwf=kwf, claim=claim, db=int(db), rho=rho, drho_dlogtau=drho_dlogtau,
        xnorm2=float((xt * xt).sum().item()),
        lam_mean=float(lamJ.mean().item()), lam_min=float(lamJ.min().item()),
        lam_max=float(lamJ.max().item()), tr=float(lamJ.sum().item()),
        JB=JB, lamJ=lamJ, Phi_b=Phi_b_t, make_W=make_W, W_amp=W_amp, W_lag=W_lag,
        # tangent-derivative blocks for the attribution simplex (efficient-score Gram)
        Vk_b=Vk_b, Vk_in=Vk_in, Vk_amp=Vk_amp, Vk_lag0=Vk_lag0, Vk_lag1=Vk_lag1,
        Min=Min, Iee=Iee, Ibe=Ibe, Ibb=Ibb)
    if want_mc:
        out.update(dict(Pt=Pt, Pabs=Pabs, Um=Um, S=S, C=C, R=R, V0=V0, V1=V1, Vinv=Vinv,
                        V0d=V0d, Vinvd=Vinvd, Cd=Cd, x=xt))
    return out


# ------------------------------------------------------------------ ENGINE B (independent Fisher)
def fisher_engine_B(cell):
    """Independent profiled beyond-band Fisher via the whitened vech-derivative construction
    (R41 sec 3.1): D~_B = (I - Pi_eta) W^{1/2} D_B with W = V^{-1} (x)_s V^{-1}; J_B = 1/2 D~_B^T D~_B.

    Computed through a DIFFERENT code path (explicit Frobenius-whitening of each dV followed by a QR
    projection out of the nuisance span) than ENGINE A (einsum traces + Schur complement). Agreement
    within 10% on eigenvalues/d'/T_det is the D0 cross-check.
    """
    Vinv_np = cell["Vinvd"].cpu().numpy()
    idx_i, idx_j = np.tril_indices(M)
    wt = np.where(idx_i == idx_j, 1.0, np.sqrt(2.0))   # sqrt(2) off-diagonal so vech preserves <.,.>_F
    # G with G^T G = Vinv (upper Cholesky); then <vech(G A G^T), vech(G B G^T)>_wt = tr(Vinv A Vinv B),
    # the Gaussian covariance Fisher metric (matched to ENGINE A's 1/2 tr(Vinv dV_a Vinv dV_b)).
    G = np.linalg.cholesky(Vinv_np + 1e-15 * np.eye(M)).T

    def whiten(dV):
        WdW = G @ dV.cpu().numpy() @ G.T
        return (WdW[idx_i, idx_j] * wt) / np.sqrt(2.0)

    Db = np.stack([whiten(cell["Vk_b"][k]) for k in range(cell["Vk_b"].shape[0])], axis=1)   # (L, db)
    Din = np.stack([whiten(cell["Vk_in"][k]) for k in range(cell["Vk_in"].shape[0])], axis=1)
    Damp = whiten(cell["Vk_amp"][0])[:, None]
    L = Db.shape[0]
    db = Db.shape[1]

    # mean-channel nuisance derivative (in-band scene), whitened by the same G (G^T G = Vinv);
    # beyond-band has an EXACT zero mean derivative (the wall), so its mean block is zero.
    Min = cell["Min"].cpu().numpy()                    # (M, d_eta)
    Dmean_in = G @ Min                                 # (M, d_eta)
    dme = Dmean_in.shape[1]

    # augment cov (length L) and mean (length M) channels into one whitened observation space
    Db_full = np.concatenate([Db, np.zeros((M, db))], axis=0)                          # (L+M, db)
    Deta_cov = np.concatenate([Din, Damp], axis=1)                                     # (L, d_eta+1)
    Deta_full = np.concatenate([
        np.concatenate([Deta_cov, np.zeros((M, Deta_cov.shape[1]))], axis=0),          # cov nuisance
        np.concatenate([np.zeros((L, dme)), Dmean_in], axis=0)                         # mean nuisance
    ], axis=1)                                                                         # (L+M, d_eta+1+d_eta)

    Q, _ = np.linalg.qr(Deta_full)
    R_perp = Db_full - Q @ (Q.T @ Db_full)             # (I - Pi_eta) D_B
    # whiten() already embeds the Gaussian 1/2 (via /sqrt(2) per vector), so <Db_a,Db_b> = Ibb[a,b];
    # the profiled Gram is J_B directly (no extra 1/2).
    JB_B = R_perp.T @ R_perp
    JB_B = 0.5 * (JB_B + JB_B.T)
    lam = np.clip(np.linalg.eigvalsh(JB_B), 0, None)[::-1]
    return dict(JB=JB_B, lam=lam, lam_mean=float(lam.mean()), lam_min=float(lam.min()),
                lam_max=float(lam.max()), tr=float(lam.sum()))


# ------------------------------------------------------------------ Gaussian bank generator
def gen_records(cell, n_rec, T_eff, kind="H0", delta_np=None, sf_scale=1.0, tau_scale=1.0,
                rng=0, want_lag=False, bank_rho=None):
    """n_rec records, each T_eff banks. Returns lag-0 sample covariance (and mean, optional lag-1).

    kind: H0 / beyond / inband; amplitude via sf_scale; medium-timescale via tau_scale.
    bank_rho: cross-bank medium correlation. Default None -> PRIMARY frozen ledger = INDEPENDENT banks
    (bank_rho=0, fully vectorized: matches "quasi-static within a bank; independent between banks").
    For the LAG / tau-attribution path pass bank_rho = exp(-1/(TAU*tau_scale)) (or want_lag=True which
    sets it) so the lag-1 cross-bank covariance E[Lag1] = bank_rho*C carries the tau signature.
    Records are processed in PARALLEL (batched over n_rec); only the OU path loops over T_eff banks.
    """
    Um, x = cell["Um"], cell["x"]
    S = cell["S"] * (sf_scale ** 2)               # amplitude change scales the power spectrum
    R = cell["R"]
    if bank_rho is None:
        bank_rho = float(np.exp(-1.0 / (TAU * tau_scale))) if (want_lag or tau_scale != 1.0) else 0.0
    xt = x if delta_np is None else x + torch.tensor(delta_np, device=DEV, dtype=DT)
    A2 = cell["Pt"] * xt[None, :]
    H = (Um.t() @ A2.t()).t()                     # (M, dm): P diag(x+delta) U_med  (medium loading)
    sdS = torch.sqrt(S)
    mfull = cell["Pt"] @ xt
    rr = torch.sqrt(R)
    dm = Um.shape[1]
    g = torch.Generator(device=DEV)
    g.manual_seed(int(rng))

    # chunk over RECORDS so the peak (chunk, T_eff, max(dm,M)) intermediates fit the 8 GB card
    # (shared with other jobs). Budget ~0.6 GB of working set per chunk.
    per_rec_bytes = T_eff * max(dm, M) * 8 * 7
    rchunk = max(1, min(n_rec, int(0.6e9 / max(per_rec_bytes, 1))))
    Scov = torch.empty(n_rec, M, M, device=DEV, dtype=DT)
    Mbar = torch.empty(n_rec, M, device=DEV, dtype=DT)
    Lag1 = torch.empty(n_rec, M, M, device=DEV, dtype=DT) if want_lag else None
    cf = np.sqrt(1 - bank_rho ** 2) if bank_rho > 1e-9 else None
    done = 0
    while done < n_rec:
        nb = min(rchunk, n_rec - done)
        if bank_rho <= 1e-9:
            Z = torch.randn(nb, T_eff, dm, device=DEV, dtype=DT, generator=g)
        else:
            Z = torch.empty(nb, T_eff, dm, device=DEV, dtype=DT)
            Z[:, 0] = torch.randn(nb, dm, device=DEV, dtype=DT, generator=g)
            for t in range(1, T_eff):
                Z[:, t] = bank_rho * Z[:, t - 1] + cf * torch.randn(nb, dm, device=DEV, dtype=DT, generator=g)
        med = (Z * sdS[None, None, :]) @ H.t()        # (nb, T_eff, M)
        del Z
        shot = torch.randn(nb, T_eff, M, device=DEV, dtype=DT, generator=g) * rr[None, None, :]
        B = mfull[None, None, :] + med + shot
        del med, shot
        Bc = B - B.mean(1, keepdim=True)
        Mbar[done:done + nb] = B.mean(1)
        del B
        Scov[done:done + nb] = torch.einsum("rti,rtj->rij", Bc, Bc) / T_eff
        if want_lag:
            Lag1[done:done + nb] = torch.einsum("rti,rtj->rij", Bc[:, :-1], Bc[:, 1:]) / (T_eff - 1)
        del Bc
        done += nb
    if want_lag:
        return Scov, Mbar, Lag1
    return Scov, Mbar


# ------------------------------------------------------------------ PHYSICAL Poisson MC (D0)
def physical_banks(cell, T_eff, delta_np=None, sf_scale=1.0, seed=0, freeze_medium=False,
                   medium_model="gaussian"):
    """Generate T_eff PHYSICAL banks via 256 complementary nonneg Poisson exposures each (independent
    banks, the frozen ledger). Each signed code P_i = P_pos - P_neg (nonneg); the detector sees photon
    counts through P_pos and P_neg separately; the signed effective bucket b_i = (n_pos - n_neg)/scp.
    Ground-truth physical ledger; used in D0 to check the |P|x shot term and the C + R covariance.

    medium_model: "gaussian" = the engine's linearized 2nd-order law w = 1 + U_med(sqrt(S) z), so the
    empirical bucket covariance matches V0 = C + R exactly (the D0.3 ledger law; p1_fisher.py:
    "lognormal harmonics are a mismatch arm, not the Fisher-prognosis law"). "lognormal" = the physical
    nonneg realization w = exp(g - var/2) (a D5 / non-Gaussian mismatch axis, NOT the primary law).
    Returns the signed bucket matrix B (T_eff, M) in physical (photon-normalized) units.
    """
    rng = np.random.default_rng(seed)
    P = P_FIXED if "Pt" not in cell else cell["Pt"].cpu().numpy()
    Ppos = np.clip(P, 0, None)
    Pneg = np.clip(-P, 0, None)
    x = cell["x"].cpu().numpy() if delta_np is None else \
        cell["x"].cpu().numpy() + np.asarray(delta_np, float)
    Um_np = cell["Um"].cpu().numpy()
    S_np = (cell["S"].cpu().numpy()) * (sf_scale ** 2)
    dm = Um_np.shape[1]
    sdS = np.sqrt(S_np)
    var_n = (S_np[None, :] * (Um_np ** 2)).sum(axis=1)
    flux_mean = (np.abs(P) @ x).mean()
    scp = PHOT / flux_mean            # PHOT detected photoelectrons/bucket at the mean throughput
    B = np.empty((T_eff, M))
    for t in range(T_eff):
        z = sdS * rng.standard_normal(dm)                       # independent bank
        if freeze_medium:
            w = np.ones(N)
        elif medium_model == "gaussian":
            w = 1.0 + Um_np @ z                                 # engine 2nd-order law: Cov(w)=U diag(S) U^T
        else:
            w = np.exp(Um_np @ z - 0.5 * var_n)                # lognormal (nonneg mismatch axis)
        wx = np.clip(w, 0, None) * x                           # physical intensity is nonneg
        n_pos = rng.poisson(np.clip(Ppos @ wx, 0, None) * scp)
        n_neg = rng.poisson(np.clip(Pneg @ wx, 0, None) * scp)
        B[t] = (n_pos - n_neg) / scp
    return B


# ------------------------------------------------------------------ cell grids
def full_grid_81():
    """The exhaustive 81-cell analytic grid (3 sigma_f x 3 spectrum x 3 k_w/k_p x 3 claim)."""
    cells = []
    for sf in SIGMA_F_LEVELS:
        for shape in SPECTRA:
            for kwf in KW_OVER_KP:
                for claim in CLAIMS:
                    cells.append(dict(sf=sf, shape=shape, kwf=kwf, claim=claim))
    return cells


def oa_grid_27():
    """Precommitted 27-cell 3^(4-1) orthogonal-array MC subset: claim level = (i+j+k) mod 3.
    Every level of sigma_f, spectrum, k_w/k_p, and claim appears 9 times; every factor PAIR balanced.
    """
    cells = []
    for i, sf in enumerate(SIGMA_F_LEVELS):
        for j, shape in enumerate(SPECTRA):
            for k, kwf in enumerate(KW_OVER_KP):
                cl = CLAIMS[(i + j + k) % 3]
                cells.append(dict(sf=sf, shape=shape, kwf=kwf, claim=cl))
    return cells


# best / mid / floor reference cells (from the corrected DETECTION_MAP)
BEST_CELL = dict(sf=1.0, shape="k^-2", kwf=1, claim=1.25)
MID_CELL = dict(sf=0.3, shape="flat", kwf=1, claim=1.5)
FLOOR_CELL = dict(sf=0.3, shape="flat", kwf=1, claim=1.8)


# ------------------------------------------------------------------ frozen thresholds (R41 sec 4.6)
BARS = dict(
    D0=dict(mean_deriv_rel_max=1e-10, two_engine_rel_max=0.10, mc_shot_rel_max=0.10),
    D1=dict(ratio_lo=0.80, ratio_hi=1.20, median_are_max=0.10, hard_lo=0.70, hard_hi=1.30),
    D2=dict(dp_strong=5.0, eps2_min_cells=77, eps2_total_cells=81, eps2_mc_lb_min=0.90,
            best_cell_Tdet_max=600, eps5_worstmode_frac=0.90, eps1_best_Tdet_max=2048,
            eps05_dp_floor=3.0),
    D3=dict(fa_prob=0.01, target_tpr_min=0.90, offtarget_fa_max=0.05, bal_acc_min=0.90,
            nontarget_dp_max=0.5, intended_dp_min=5.0, simplex_canon_corr_max=0.10),
    D4=dict(fixed_retain_ratio=1.20),
    D5=dict(auc_loss_max=0.05, Tdet_inflation_max=0.25, nontarget_fa_max=0.05),
    D6=dict(online_frac_of_bank=0.10, mem_mb_max=10.0),
)

# anomaly noncentrality convention: d'^2 = T_eff * (c^T JB c); strong detection d'>=5 -> 25.
DP_STRONG2 = 25.0


def energy_spread_delta(cell_geo_db, eps, xnorm, seed):
    """Isotropic beyond-band delta with ||delta|| = eps*||x||; returns coeffs c (in U_beta)."""
    rng = np.random.default_rng(seed)
    c = rng.standard_normal(cell_geo_db)
    c = c / np.linalg.norm(c) * eps * xnorm
    return c


def worstmode_delta(cell, eps, xnorm):
    """Least-informed (min-eigenvector) beyond-band delta -- the adversarial family."""
    evals, evecs = torch.linalg.eigh(cell["JB"])
    c = evecs[:, 0].cpu().numpy()               # min-eigenvalue direction
    return c / np.linalg.norm(c) * eps * xnorm


def sha256_arr(a):
    return hashlib.sha256(np.ascontiguousarray(a, np.float64).tobytes()).hexdigest()


if __name__ == "__main__":
    print(f"SEALED_DET frozen ledger: n={n} N={N} k_p={KP} M={M} exposures/bank={N_COMPLEMENTARY}")
    print(f"  {MS_PER_BANK:.1f} ms/bank, T_cap={T_CAP} banks = {T_CAP_SEC:.2f} s, "
          f"PHOT={PHOT:.0e}/bucket, tau={TAU} rho={RHO:.4f}")
    print(f"  in-band eta dof={d_eta}; beyond-band dof per claim: "
          f"{ {c: BETA_np[c].shape[1] for c in CLAIMS} }")
    print(f"  81-cell grid size={len(full_grid_81())}, 27-cell OA size={len(oa_grid_27())}")
    # OA balance self-check
    oa = oa_grid_27()
    for key, levels in [("sf", SIGMA_F_LEVELS), ("shape", SPECTRA),
                        ("kwf", KW_OVER_KP), ("claim", CLAIMS)]:
        counts = {lv: sum(1 for c in oa if c[key] == lv) for lv in levels}
        print(f"    OA balance {key}: {counts}")
    print(f"  DEV={DEV}")
