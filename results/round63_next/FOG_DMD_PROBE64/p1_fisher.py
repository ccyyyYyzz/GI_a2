# P1 -- PRELAUNCH PROFILED-FISHER PROGNOSIS at the R39 FROZEN geometry (ruling section 4 + 6.5).
#
# THE go/no-go BEFORE any reconstruction campaign. Computes the exact profiled covariance Fisher
# on the beyond-band scene subspace (eq 4.1-4.4) and evaluates bar P1 (section 6.5):
#   - no exact profiled null on the primary annulus (1.0-1.8 k_p);
#   - profiling efficiency eta_prof: p10 >= 0.10, median >= 0.25 (eq 4.10);
#   - CRB-predicted beyond-band NRMSE <= 0.25 synthetic / <= 0.35 natural median (eq 4.7);
#   - >= 70% of primary beyond-band modes at Fisher SNR >= 3 (eq 4.9).
# IF P1 FAILS -> STOP, archive the theorem, no reconstruction (kill-tree node 2).
#
# FROZEN geometry (ruling 6.1): N=4096 (64x64), M=128 signed band-limited codes, pattern band
# |kx|,|ky|<=k_p=5 (real Fourier box, ~121 dof), medium primary band k_w=k_p=5 (SAME cutoff),
# support aperture <=10 (=2 k_p), primary claim <= 1.8 k_p, sigma_f=0.30, n=1e4 photons/bucket,
# T_eff grid {256..4096}. Medium = flat band-limited stationary covariance (bounded field's
# 2nd-order law; lognormal harmonics are a mismatch arm, not the Fisher-prognosis law).
#
# Basis: REAL Fourier lattice modes on the 64x64 torus (matches the ruling's k-space geometry and
# the 121-dof in-band count), not DCT. Codes/medium/scene all live in these bands.
import json
import time

import numpy as np
import torch

DEV = "cuda" if torch.cuda.is_available() else "cpu"
DT = torch.float64
n = 64
N = n * n
KP = 5                       # pattern & medium cutoff (k_w = k_p)
K_CLAIM = 9                  # 1.8 * k_p (primary usable claim edge)
K_SUPPORT = 10               # 2.0 * k_p (support aperture)
SIG_F = 0.30
PHOT = 1e4                   # primary ceiling (ruling: 1e4 already essentially clean)
M = 128
TAU = 8.0
RHO = float(np.exp(-1.0 / TAU))
t0 = time.time()

_X, _Y = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")


def band_modes(k_lo, k_hi):
    """Orthonormal REAL Fourier modes with Chebyshev freq max(|kx|,|ky|) in [k_lo,k_hi]."""
    cols = []
    for kx in range(-k_hi, k_hi + 1):
        for ky in range(-k_hi, k_hi + 1):
            m = max(abs(kx), abs(ky))
            if k_lo <= m <= k_hi:
                ph = 2 * np.pi * (kx * _X + ky * _Y) / n
                cols.append(np.cos(ph).ravel())
                cols.append(np.sin(ph).ravel())
    A = np.array(cols).T
    Q, R = np.linalg.qr(A)
    keep = np.abs(np.diag(R)) > 1e-8
    return Q[:, keep]


# subspaces
U_med = band_modes(1, KP)                 # medium band (DC excluded; zero-mean fluctuation)
db = U_med.shape[1]
U_in = band_modes(0, KP)                  # in-band scene nuisance (includes DC)
U_beta = band_modes(KP + 1, K_CLAIM)      # primary beyond-band annulus (1.0-1.8 k_p)
d_beta = U_beta.shape[1]
coeff_sd2 = SIG_F ** 2 * N / db
Kw = coeff_sd2 * (U_med @ U_med.T)        # flat-spectrum band-limited medium covariance (N x N)
print(f"P1 frozen geometry: N={N} k_p={KP} claim<= {K_CLAIM} support<= {K_SUPPORT} "
      f"M={M} sigma_f={SIG_F} n={PHOT:.0e}", flush=True)
print(f"  dof: medium db={db}  in-band eta={U_in.shape[1]}  beyond-band beta(1.0-1.8kp)={d_beta}", flush=True)


def signed_codes(seed=10):
    """M signed band-limited random codes (physical complementary nonneg pairs; the signed
    effective row is what enters the Fisher). Band-limited to |k|<=k_p."""
    rp = np.random.default_rng(seed)
    Uc = band_modes(1, KP)                # exclude DC so codes are zero-mean (signed)
    C = rp.standard_normal((M, Uc.shape[1]))
    P = C @ Uc.T                          # (M, N)
    P /= np.abs(P).max(axis=1, keepdims=True)
    return P


P = signed_codes()


def witness_scene(seed=4):
    rs = np.random.default_rng(seed)
    a_in = U_in @ rs.standard_normal(U_in.shape[1])
    a_be = U_beta @ (2.0 * rs.choice([-1, 1], d_beta))
    x = a_in + a_be
    x = x - x.min(); return x / x.max()


def natural_scene(name):
    from skimage import data, transform
    src = {"cameraman": data.camera, "coins": data.coins, "moon": data.moon}[name]()
    img = transform.resize(src.astype(np.float64), (n, n), anti_aliasing=True, mode="reflect")
    return ((img - img.min()) / (np.ptp(img) + 1e-12)).ravel()


def fisher_prognosis(x, T_eff, tag):
    """Exact profiled covariance Fisher on the beyond-band subspace at scene x."""
    xt = torch.tensor(x, device=DEV, dtype=DT)
    Pt = torch.tensor(P, device=DEV, dtype=DT)
    Kwt = torch.tensor(Kw, device=DEV, dtype=DT)
    A = Pt * xt[None, :]                                  # P diag(x)  (M,N)
    C = A @ Kwt @ A.t()                                   # covariance signal (M,M)
    m = Pt @ xt                                           # signed mean bucket (~0); kept for the eta mean term
    flux = torch.abs(Pt) @ xt                             # CORRECTED shot: physical nonneg photon throughput |P|x
    R_shot = torch.diag(flux * (flux.mean() / PHOT))      # var(b_i)=flux_i*mean(flux)/n (signed complementary-pair shot)
    V = C + R_shot
    Vinv = torch.linalg.inv(V + 1e-9 * torch.eye(M, device=DEV, dtype=DT))
    KA = Kwt @ A.t()                                      # (N,M)

    def Dmat(U):                                          # stack dV/dcoeff for each column of U
        d = U.shape[1]
        Gs = torch.empty((d, M, M), device=DEV, dtype=DT)
        Ut = torch.tensor(U, device=DEV, dtype=DT)
        for k in range(d):
            phi = Ut[:, k]
            B = (Pt * phi[None, :]) @ KA                  # (M,M)
            Vk = B + B.t()
            Gs[k] = Vinv @ Vk
        return Gs, Ut

    Gb, _ = Dmat(U_beta)                                  # beyond-band derivatives
    # nuisance eta = in-band scene modes + law amplitude (dV/dlogamp = C)
    Gin, Uint = Dmat(U_in)
    Glaw = (Vinv @ C).unsqueeze(0)
    Geta = torch.cat([Gin, Glaw], 0)
    # covariance Fisher blocks: I[a,b] = (T_eff/2) tr(G_a G_b)
    half = T_eff / 2.0
    Ibb = half * torch.einsum("aij,bji->ab", Gb, Gb)
    Ibe = half * torch.einsum("aij,bji->ab", Gb, Geta)
    Iee = half * torch.einsum("aij,bji->ab", Geta, Geta)
    # mean-channel Fisher for eta (beyond-band mean deriv is exactly 0 -> Theorem 2)
    Min = Pt @ Uint                                       # (M, n_in)  = dm/deta
    Iee[: Uint.shape[1], : Uint.shape[1]] += T_eff * (Min.t() @ Vinv @ Min)
    # profiled (efficient) beyond-band Fisher  J_B = Ibb - Ibe Iee^+ Ieb
    Iee_pinv = torch.linalg.pinv(Iee, rcond=1e-10)
    JB = Ibb - Ibe @ Iee_pinv @ Ibe.t()
    JB = 0.5 * (JB + JB.t())

    lam_JB = torch.linalg.eigvalsh(JB).clamp(min=0).cpu().numpy()[::-1]
    lam_Ibb = torch.linalg.eigvalsh(0.5 * (Ibb + Ibb.t())).clamp(min=0).cpu().numpy()[::-1]
    eta_prof = lam_JB / np.maximum(lam_Ibb, 1e-30)       # eq 4.10 (sorted eigenvalue ratios)
    Jbar = JB / T_eff                                    # per-effective-epoch
    lam_bar = lam_JB / T_eff
    # scene beyond-band coefficients a_j in the JB eigenbasis
    _, Vec = torch.linalg.eigh(JB)
    beta_coeff = (torch.tensor(U_beta.T @ x, device=DEV, dtype=DT))   # scene coeffs in U_beta
    a_in_eig = (Vec.t() @ beta_coeff).cpu().numpy()[::-1]
    fisher_snr2 = T_eff * (a_in_eig ** 2) * lam_bar      # eq 4.9 argument
    f_rec = float(np.mean(fisher_snr2 >= 9.0))           # gamma=3
    # CRB-predicted beyond-band NRMSE (eq 4.7), W_B = I, on the non-null subspace
    pos = lam_JB > 1e-8 * lam_JB.max()
    tr_inv = float(np.sum(1.0 / lam_JB[pos]))            # tr(JB^-1) on positive subspace
    beta_norm2 = float(np.sum(beta_coeff.cpu().numpy() ** 2))
    nmse_crb = float(np.sqrt(tr_inv / max(beta_norm2, 1e-12)))
    min_eig_bar = float(lam_bar.min())
    return dict(tag=tag, T_eff=T_eff, d_beta=int(d_beta),
                min_eig_Jbar=min_eig_bar, n_positive=int(pos.sum()),
                eta_prof_p10=float(np.percentile(eta_prof, 10)),
                eta_prof_median=float(np.median(eta_prof)),
                f_rec_snr3=f_rec, nmse_crb=nmse_crb,
                exact_null=bool(min_eig_bar <= 1e-12))


# evaluate at the primary resource point (T_eff=4096, n=1e4) on witness + naturals
print("\n[P1] profiled-Fisher prognosis at T_eff=4096, n=1e4, sigma_f=0.30", flush=True)
rows = []
xw = witness_scene()
rows.append(fisher_prognosis(xw, 4096, "witness"))
for nm in ["cameraman", "coins", "moon"]:
    try:
        rows.append(fisher_prognosis(natural_scene(nm), 4096, f"nat_{nm}"))
    except Exception as e:
        print(f"   (skip {nm}: {e})", flush=True)
for r in rows:
    print(f"  {r['tag']:12s} null={r['exact_null']} min-eig(Jbar)={r['min_eig_Jbar']:.2e}  "
          f"eta_prof p10={r['eta_prof_p10']:.3f} med={r['eta_prof_median']:.3f}  "
          f"f_rec(SNR>=3)={r['f_rec_snr3']:.2f}  NRMSE_CRB={r['nmse_crb']:.3f}", flush=True)

# T_eff sweep on the witness (information rate)
print("\n[P1] T_eff sweep (witness) -- NRMSE_CRB and f_rec vs independent medium states", flush=True)
sweep = []
for Te in [256, 512, 1024, 2048, 4096]:
    r = fisher_prognosis(xw, Te, f"witness_T{Te}")
    sweep.append(r)
    print(f"  T_eff={Te:5d}  NRMSE_CRB={r['nmse_crb']:.3f}  f_rec(SNR>=3)={r['f_rec_snr3']:.2f}  "
          f"eta_prof med={r['eta_prof_median']:.3f}", flush=True)

# ============================ P1 bar (ruling 6.5) ============================
wit = [r for r in rows if r["tag"] == "witness"][0]
nat = [r for r in rows if r["tag"].startswith("nat_")]
nat_nmse_med = float(np.median([r["nmse_crb"] for r in nat])) if nat else float("nan")
checks = dict(
    no_exact_null=(not wit["exact_null"]),
    eta_p10_ge_0p10=(wit["eta_prof_p10"] >= 0.10),
    eta_median_ge_0p25=(wit["eta_prof_median"] >= 0.25),
    nrmse_synth_le_0p25=(wit["nmse_crb"] <= 0.25),
    nrmse_nat_median_le_0p35=(np.isfinite(nat_nmse_med) and nat_nmse_med <= 0.35),
    frec_ge_0p70=(wit["f_rec_snr3"] >= 0.70),
)
P1_pass = all(checks.values())
verdict = "P1_PASS -> proceed to reconstruction" if P1_pass else "P1_FAIL -> STOP (no reconstruction; archive theorem)"
out = dict(bar="P1_fisher_prognosis", frozen_geometry=dict(
    N=N, k_p=KP, k_claim=K_CLAIM, k_support=K_SUPPORT, M=M, sigma_f=SIG_F, photons=PHOT,
    medium_db=int(db), beyond_band_dof=int(d_beta)),
    primary_point=rows, T_eff_sweep=sweep, natural_nrmse_crb_median=nat_nmse_med,
    checks=checks, P1_pass=bool(P1_pass), verdict=verdict, wall_s=time.time() - t0)
json.dump(out, open("P1_results_CORRECTED.json", "w"), indent=2)  # physical-shot corrected artifact
print("\n[P1] bar checks:", flush=True)
for k, v in checks.items():
    print(f"  {'PASS' if v else 'FAIL'}  {k}", flush=True)
print(f"\nP1 VERDICT: {verdict}  [{time.time()-t0:.0f}s]", flush=True)
