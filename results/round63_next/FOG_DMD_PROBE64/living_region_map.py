# LIVING-REGION MAP -- where (if anywhere) in physical parameter space does the beyond-band method
# reach P1-grade Fisher viability? (post-P1-kill divergence, R39 ruling section 4 analytics; NO
# reconstruction.) Reuses the p1_fisher.py exact profiled-Fisher engine, swept over:
#   sigma_f in {0.3, 0.6, 1.0}     (1.0 = fully developed speckle; bounded field, clip documented)
#   spectrum shape in {flat, k^-1, k^-2}   (per-mode power S_w(|k|) within the medium band)
#   medium band k_w in {1,2,4} x k_p       (total variance fixed -> wider band dilutes per-mode power)
#   claim shell in {1.25, 1.5, 1.8} x k_p
# Fixed: N=4096, M=128, k_p=5, n=1e4, T_eff<=4096.
# Per cell: f_rec(gamma=3,T_eff=4096), NRMSE_CRB(witness + natural median), eta_prof p10/median,
#           T_req(0.30) (the price tag, finite even where f_rec fails).
# Verdict: REGION_FOUND (passing cells + most defensible) or REGION_EMPTY (tombstone).
# Plus a Monte-Carlo covariance CROSS-CHECK of the analytic Fisher engine on representative cells.
import json
import time

import numpy as np
import torch

DEV = "cuda" if torch.cuda.is_available() else "cpu"
DT = torch.float64
n = 64
N = n * n
KP = 5
M = 128
PHOT = 1e4
t0 = time.time()
_X, _Y = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")


def band_modes(k_lo, k_hi):
    """Orthonormal REAL Fourier modes, Chebyshev freq max(|kx|,|ky|) in [k_lo,k_hi] (QR)."""
    cols = []
    for kx in range(-k_hi, k_hi + 1):
        for ky in range(-k_hi, k_hi + 1):
            m = max(abs(kx), abs(ky))
            if k_lo <= m <= k_hi:
                ph = 2 * np.pi * (kx * _X + ky * _Y) / n
                cols.append(np.cos(ph).ravel()); cols.append(np.sin(ph).ravel())
    Q, R = np.linalg.qr(np.array(cols).T)
    return Q[:, np.abs(np.diag(R)) > 1e-8]


def medium_modes(k_w):
    """Per-frequency orthonormal real Fourier modes for the medium band (Chebyshev |k|<=k_w),
    with the radial |k| of each mode (for the spectrum-shape power weighting). Half-plane
    enumeration avoids cos/sin duplicates; columns normalized to unit norm."""
    cols, kabs = [], []
    for kx in range(-k_w, k_w + 1):
        for ky in range(-k_w, k_w + 1):
            if max(abs(kx), abs(ky)) == 0:
                continue
            if (kx < 0) or (kx == 0 and ky < 0):        # half-plane (drop conjugate duplicate)
                continue
            ph = 2 * np.pi * (kx * _X + ky * _Y) / n
            for f in (np.cos(ph).ravel(), np.sin(ph).ravel()):
                nrm = np.linalg.norm(f)
                if nrm > 1e-8:
                    cols.append(f / nrm); kabs.append(np.hypot(kx, ky))
    U = np.array(cols).T
    return U, np.array(kabs)


def power_spectrum(kabs, shape, sig_eff):
    w = {"flat": np.ones_like(kabs), "k^-1": 1.0 / kabs, "k^-2": 1.0 / kabs ** 2}[shape]
    S = w / w.sum() * (N * sig_eff ** 2)                # (1/N) sum S = sig_eff^2 (pixel variance)
    return S


def effective_contrast(sig_f, lo, hi, ns=2_000_000, seed=0):
    """Realized contrast + clip fraction of the bounded mean-one field w=clip(1+g), g~N(0,sig_f^2).
    Pointwise (band-limited g is ~Gaussian per pixel by CLT), so depends only on (sig_f,bounds)."""
    r = np.random.default_rng(seed)
    g = sig_f * r.standard_normal(ns)
    w = np.clip(1 + g, lo, hi)
    return float(np.std(w)), float(np.mean((1 + g < lo) | (1 + g > hi)))


def signed_codes(seed=10):
    rp = np.random.default_rng(seed)
    Uc = band_modes(1, KP)
    P = (rp.standard_normal((M, Uc.shape[1])) @ Uc.T)
    P /= np.abs(P).max(axis=1, keepdims=True)
    return P


P_np = signed_codes()
Pt = torch.tensor(P_np, device=DEV, dtype=DT)
U_in_np = band_modes(0, KP)                             # in-band scene nuisance
Phi_eta = torch.tensor(U_in_np, device=DEV, dtype=DT)
d_eta = U_in_np.shape[1]
# claim-shell beyond-band bases (fixed by k_p)
CLAIMS = {1.25: 6, 1.5: 7, 1.8: 9}
BETA = {c: torch.tensor(band_modes(KP + 1, hi), device=DEV, dtype=DT) for c, hi in CLAIMS.items()}


def witness_scene(seed=4):
    rs = np.random.default_rng(seed)
    a_in = U_in_np @ rs.standard_normal(U_in_np.shape[1])
    Ube = band_modes(KP + 1, 9)
    a_be = Ube @ (2.0 * rs.choice([-1, 1], Ube.shape[1]))
    x = a_in + a_be; x -= x.min(); return x / x.max()


def natural_scene(name):
    from skimage import data, transform
    src = {"cameraman": data.camera, "coins": data.coins, "moon": data.moon}[name]()
    img = transform.resize(src.astype(np.float64), (n, n), anti_aliasing=True, mode="reflect")
    return ((img - img.min()) / (np.ptp(img) + 1e-12)).ravel()


SCENES = {"witness": witness_scene(), "cameraman": natural_scene("cameraman"),
          "coins": natural_scene("coins"), "moon": natural_scene("moon")}
SCENES_T = {k: torch.tensor(v, device=DEV, dtype=DT) for k, v in SCENES.items()}


def fisher_metrics(x, Um, S, Phi_beta, T_eff=4096, n_ph=PHOT):
    """Exact profiled covariance Fisher on the beyond-band subspace for one (scene, medium) cell."""
    A = Pt * x[None, :]
    UmA = Um.t() @ A.t()                                # (db,M)
    KA = Um @ (S[:, None] * UmA)                        # K_w A^T  (N,M)
    C = A @ KA
    m = Pt @ x                                            # signed mean bucket (~0); kept for the eta mean term
    flux = torch.abs(Pt) @ x                              # CORRECTED shot: physical nonneg photon throughput |P|x
    R = flux * (flux.mean() / n_ph)                       # var(b_i)=flux_i*mean(flux)/n (signed complementary-pair shot)
    V = C + torch.diag(R)
    Vinv = torch.linalg.inv(V + 1e-9 * torch.eye(M, device=DEV, dtype=DT))

    def Gstack(Phi):
        d = Phi.shape[1]
        G = torch.empty((d, M, M), device=DEV, dtype=DT)
        for k in range(d):
            B = (Pt * Phi[:, k][None, :]) @ KA
            G[k] = Vinv @ (B + B.t())
        return G

    Gb = Gstack(Phi_beta)
    Ge = Gstack(Phi_eta)
    Glaw = (Vinv @ C).unsqueeze(0)
    Geta = torch.cat([Ge, Glaw], 0)
    half = T_eff / 2.0
    Ibb = half * torch.einsum("aij,bji->ab", Gb, Gb)
    Ibe = half * torch.einsum("aij,bji->ab", Gb, Geta)
    Iee = half * torch.einsum("aij,bji->ab", Geta, Geta)
    Min = Pt @ Phi_eta
    Iee[:d_eta, :d_eta] += T_eff * (Min.t() @ Vinv @ Min)
    JB = Ibb - Ibe @ torch.linalg.pinv(Iee, rcond=1e-10) @ Ibe.t()
    JB = 0.5 * (JB + JB.t())
    evals, evecs = torch.linalg.eigh(JB)
    lam = torch.clamp(evals, min=0)
    lam_Ibb = torch.clamp(torch.linalg.eigvalsh(0.5 * (Ibb + Ibb.t())), min=0)
    lam_s = torch.sort(lam, descending=True).values
    lIbb_s = torch.sort(lam_Ibb, descending=True).values
    eta_prof = (lam_s / torch.clamp(lIbb_s, min=1e-30)).cpu().numpy()
    beta_coeff = Phi_beta.t() @ x
    a_eig = evecs.t() @ beta_coeff
    snr2 = (a_eig ** 2) * lam                           # = T_eff a^2 lam_bar (lam already ~ T_eff)
    f_rec = float((snr2 >= 9.0).float().mean().item())
    beta_norm2 = float((beta_coeff ** 2).sum().item())
    pos = lam > 1e-8 * lam.max()
    lam_bar = lam / T_eff
    tr_inv_bar = float((1.0 / lam_bar[pos]).sum().item())
    nmse_crb = float(np.sqrt(tr_inv_bar / (T_eff * max(beta_norm2, 1e-12))))
    t_req = tr_inv_bar / (0.30 ** 2 * max(beta_norm2, 1e-12))
    return dict(f_rec=f_rec, nmse_crb=nmse_crb, t_req_0p30=t_req,
                eta_p10=float(np.percentile(eta_prof, 10)), eta_med=float(np.median(eta_prof)),
                min_eig_bar=float((lam.min() / T_eff).item()),
                exact_null=bool((lam.min() / T_eff).item() <= 1e-12))


# ---- effective contrast per sigma_f (bounded-field clip diagnostic) ----
BOUNDS = {0.3: (0.2, 1.8), 0.6: (0.2, 1.8), 1.0: (0.05, 1.95)}
SIGS = [0.3, 0.6, 1.0]
SHAPES = ["flat", "k^-1", "k^-2"]
KWS = {1: KP, 2: 2 * KP, 4: 4 * KP}
sig_eff = {}
print("=== bounded-field truncation (clip) diagnostic ===", flush=True)
for sf in SIGS:
    lo, hi = BOUNDS[sf]
    se, cf = effective_contrast(sf, lo, hi)
    sig_eff[sf] = se
    print(f"  sigma_f_nominal={sf:.1f} bounds[{lo},{hi}] -> sigma_f_eff={se:.3f}  clip-frac={cf:.3f}", flush=True)

# precompute medium bases + spectra
MEDU = {kwf: medium_modes(kw) for kwf, kw in KWS.items()}
for kwf, (U, ka) in MEDU.items():
    MEDU[kwf] = (torch.tensor(U, device=DEV, dtype=DT), ka, U.shape[1])
    print(f"  medium band k_w={kwf}x k_p: db={MEDU[kwf][2]} modes", flush=True)

# ============================ the grid ============================
print("\n=== LIVING-REGION grid (81 cells): P1 viability @ T_eff=4096, n=1e4 ===", flush=True)
grid = []
for sf in SIGS:
    for shape in SHAPES:
        for kwf in (1, 2, 4):
            Um, kabs, db = MEDU[kwf]
            S = torch.tensor(power_spectrum(kabs, shape, sig_eff[sf]), device=DEV, dtype=DT)
            for claim in (1.25, 1.5, 1.8):
                Phi_beta = BETA[claim]
                wit = fisher_metrics(SCENES_T["witness"], Um, S, Phi_beta)
                nats = [fisher_metrics(SCENES_T[k], Um, S, Phi_beta)
                        for k in ("cameraman", "coins", "moon")]
                nat_nmse = float(np.median([r["nmse_crb"] for r in nats]))
                nat_frec = float(np.median([r["f_rec"] for r in nats]))
                p1 = (not wit["exact_null"] and wit["eta_p10"] >= 0.10 and wit["eta_med"] >= 0.25
                      and wit["nmse_crb"] <= 0.25 and nat_nmse <= 0.35 and wit["f_rec"] >= 0.70)
                cell = dict(sigma_f=sf, sigma_f_eff=round(sig_eff[sf], 3), shape=shape,
                            k_w_over_kp=kwf, claim=claim, db=db,
                            f_rec_witness=round(wit["f_rec"], 3), f_rec_nat_med=round(nat_frec, 3),
                            nmse_crb_witness=round(wit["nmse_crb"], 3), nmse_crb_nat_med=round(nat_nmse, 3),
                            eta_p10=round(wit["eta_p10"], 3), eta_med=round(wit["eta_med"], 3),
                            t_req_0p30=round(wit["t_req_0p30"], 1), exact_null=wit["exact_null"],
                            P1_pass=bool(p1))
                grid.append(cell)
                flag = "  <-- P1 PASS" if p1 else ""
                print(f"  sf{sf}({sig_eff[sf]:.2f}) {shape:5s} kw{kwf} claim{claim}: "
                      f"f_rec={wit['f_rec']:.2f} NRMSE_w={wit['nmse_crb']:.2f} NRMSE_nat={nat_nmse:.2f} "
                      f"eta10={wit['eta_p10']:.2f} Treq={wit['t_req_0p30']:.0f}{flag}", flush=True)

passing = [c for c in grid if c["P1_pass"]]
print(f"\n=== {len(passing)}/81 cells reach P1-grade viability ===", flush=True)

# ============================ MC cross-check (section 7.4) ============================
print("\n=== MC covariance cross-check of the analytic Fisher engine ===", flush=True)


def mc_covariance(x_np, Um_np, S_np, sf, lo, hi, T_mc=4000, bounded=True, seed=1):
    """Empirical one-epoch bucket covariance from simulated medium epochs (Gaussian or bounded)."""
    r = np.random.default_rng(seed)
    db = Um_np.shape[1]
    sd = np.sqrt(S_np)
    Y = np.empty((T_mc, M))
    for t in range(T_mc):
        g = Um_np @ (sd * r.standard_normal(db))
        w = np.clip(1 + g, lo, hi) if bounded else (1 + g)
        Y[t] = P_np @ (w * x_np)
    return np.cov(Y.T)


xc = SCENES["witness"]
check_cells = [(0.3, "flat", 1), (1.0, "flat", 1), (0.6, "k^-2", 2)]
xcheck = []
for (sf, shape, kwf) in check_cells:
    Um_np, kabs = medium_modes(KWS[kwf])
    S_np = power_spectrum(kabs, shape, sig_eff[sf])
    lo, hi = BOUNDS[sf]
    # analytic one-epoch covariance signal C(x)
    A = P_np * xc[None, :]
    C_an = A @ (Um_np @ (S_np[:, None] * (Um_np.T @ A.T)))
    V_g = mc_covariance(xc, Um_np, S_np, sf, lo, hi, bounded=False)   # engine check (Gaussian)
    V_b = mc_covariance(xc, Um_np, S_np, sf, lo, hi, bounded=True)    # truncation check (bounded)
    err_g = float(np.linalg.norm(V_g - C_an) / np.linalg.norm(C_an))
    err_b = float(np.linalg.norm(V_b - C_an) / np.linalg.norm(C_an))
    xcheck.append(dict(cell=f"sf{sf}_{shape}_kw{kwf}", rel_err_gaussian_engine=round(err_g, 3),
                       rel_err_bounded_truncation=round(err_b, 3)))
    print(f"  {sf} {shape} kw{kwf}: analytic-vs-MC(Gaussian) rel-err {err_g:.3f} "
          f"(engine ok if small); analytic-vs-MC(bounded) rel-err {err_b:.3f} (truncation effect)",
          flush=True)

# ============================ verdict ============================
if passing:
    best = min(passing, key=lambda c: (c["k_w_over_kp"], c["claim"], -c["sigma_f"]))  # most defensible
    verdict = "REGION_FOUND"
else:
    best = None
    verdict = "REGION_EMPTY"
out = dict(map="living_region", fixed=dict(N=N, M=M, k_p=KP, photons=PHOT, T_eff=4096),
           sigma_f_eff=sig_eff, n_cells=len(grid), n_pass=len(passing),
           grid=grid, passing_cells=passing, most_defensible=best,
           mc_cross_check=xcheck, verdict=verdict, wall_s=time.time() - t0)
json.dump(out, open("LIVING_REGION_MAP_CORRECTED.json", "w"), indent=2)  # physical-shot corrected artifact
print(f"\nLIVING-REGION VERDICT: {verdict}  ({len(passing)}/81 P1-viable)  [{time.time()-t0:.0f}s]", flush=True)
if best:
    print(f"  most defensible: {best}", flush=True)
