# DETECTION robustness column (E5d analog): does detectability inherit the estimation robustness
# to declared-law mismatch? At 3 representative cells recompute the per-epoch profiled Fisher (hence
# d') with (a) a 10% medium-subspace basis rotation, (b) a tau x2 mismatch (pure T_eff rescaling).
# Self-contained (does not import the grid module, so the 81-cell grid is not re-run).
import json, numpy as np, torch
DEV = "cuda" if torch.cuda.is_available() else "cpu"; DT = torch.float64
n = 64; N = n * n; KP = 5; M = 128; PHOT = 1e4
_X, _Y = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")

def band_modes(k_lo, k_hi):
    cols = []
    for kx in range(-k_hi, k_hi + 1):
        for ky in range(-k_hi, k_hi + 1):
            if k_lo <= max(abs(kx), abs(ky)) <= k_hi:
                ph = 2 * np.pi * (kx * _X + ky * _Y) / n
                cols.append(np.cos(ph).ravel()); cols.append(np.sin(ph).ravel())
    Q, R = np.linalg.qr(np.array(cols).T); return Q[:, np.abs(np.diag(R)) > 1e-8]

def medium_modes(k_w):
    cols, kabs = [], []
    for kx in range(-k_w, k_w + 1):
        for ky in range(-k_w, k_w + 1):
            if max(abs(kx), abs(ky)) == 0: continue
            if (kx < 0) or (kx == 0 and ky < 0): continue
            ph = 2 * np.pi * (kx * _X + ky * _Y) / n
            for f in (np.cos(ph).ravel(), np.sin(ph).ravel()):
                nrm = np.linalg.norm(f)
                if nrm > 1e-8: cols.append(f / nrm); kabs.append(np.hypot(kx, ky))
    return np.array(cols).T, np.array(kabs)

def power_spectrum(kabs, shape, sig_eff):
    w = {"flat": np.ones_like(kabs), "k^-1": 1.0 / kabs, "k^-2": 1.0 / kabs ** 2}[shape]
    return w / w.sum() * (N * sig_eff ** 2)

def signed_codes(seed=10):
    rp = np.random.default_rng(seed); Uc = band_modes(1, KP)
    P = (rp.standard_normal((M, Uc.shape[1])) @ Uc.T); P /= np.abs(P).max(axis=1, keepdims=True); return P

P_np = signed_codes(); Pt = torch.tensor(P_np, device=DEV, dtype=DT)
U_in_np = band_modes(0, KP); Phi_eta = torch.tensor(U_in_np, device=DEV, dtype=DT); d_eta = U_in_np.shape[1]
CLAIMS = {1.25: 6, 1.5: 7, 1.8: 9}
BETA = {c: torch.tensor(band_modes(KP + 1, hi), device=DEV, dtype=DT) for c, hi in CLAIMS.items()}
def witness_scene(seed=4):
    rs = np.random.default_rng(seed); a_in = U_in_np @ rs.standard_normal(U_in_np.shape[1])
    Ube = band_modes(KP + 1, 9); a_be = Ube @ (2.0 * rs.choice([-1, 1], Ube.shape[1]))
    x = a_in + a_be; x -= x.min(); return x / x.max()
XW = torch.tensor(witness_scene(), device=DEV, dtype=DT)

def jbar_lam(x, Um, S, Phi_beta, T_ref=4096.0, n_ph=PHOT):
    A = Pt * x[None, :]; KA = Um @ (S[:, None] * (Um.t() @ A.t())); C = A @ KA
    m = Pt @ x; R = torch.clamp(m, min=1e-12) * (m.mean() / n_ph); V = C + torch.diag(R)
    Vinv = torch.linalg.inv(V + 1e-9 * torch.eye(M, device=DEV, dtype=DT))
    def Gstack(Phi):
        G = torch.empty((Phi.shape[1], M, M), device=DEV, dtype=DT)
        for k in range(Phi.shape[1]):
            B = (Pt * Phi[:, k][None, :]) @ KA; G[k] = Vinv @ (B + B.t())
        return G
    Gb = Gstack(Phi_beta); Ge = Gstack(Phi_eta); Glaw = (Vinv @ C).unsqueeze(0)
    Geta = torch.cat([Ge, Glaw], 0); half = T_ref / 2.0
    Ibb = half * torch.einsum("aij,bji->ab", Gb, Gb); Ibe = half * torch.einsum("aij,bji->ab", Gb, Geta)
    Iee = half * torch.einsum("aij,bji->ab", Geta, Geta)
    Min = Pt @ Phi_eta; Iee[:d_eta, :d_eta] += T_ref * (Min.t() @ Vinv @ Min)
    JB = Ibb - Ibe @ torch.linalg.pinv(Iee, rcond=1e-10) @ Ibe.t(); JB = 0.5 * (JB + JB.t())
    return torch.clamp(torch.linalg.eigvalsh(JB), min=0).cpu().numpy() / T_ref, float((x * x).sum().item())

def dprime_avg(lam, xnorm2, eps=0.01, T=4096):
    return float(np.sqrt(max(T * float(lam.mean()) * eps * eps * xnorm2, 0)))

def rotate_cols(U, eps, seed):
    rng = np.random.default_rng(seed); R = rng.standard_normal(U.shape)
    R = R - U @ (U.T @ R); R = R / (np.linalg.norm(R, axis=0, keepdims=True) + 1e-12)
    Ue, _ = np.linalg.qr(U * np.sqrt(1 - eps ** 2) + R * eps)
    ang = np.arccos(np.clip(np.linalg.svd(U.T @ Ue, compute_uv=False), -1, 1))
    return Ue, float(np.sin(ang).mean())

def T_eff_ratio(tt, ta):
    f = lambda t: (1 - np.exp(-1.0 / t)) / (1 + np.exp(-1.0 / t)); return f(tt) / f(ta)

sig_eff = {0.3: 0.298}   # from bounded-field diagnostic (matches living_region_map)
KWS = {1: KP, 2: 2 * KP, 4: 4 * KP}
CELLS = [(0.3, "k^-2", 1, 1.25), (0.3, "flat", 1, 1.5), (0.3, "flat", 1, 1.8)]
rows = []
for (sf, shape, kwf, claim) in CELLS:
    Um_np, kabs = medium_modes(KWS[kwf]); S = torch.tensor(power_spectrum(kabs, shape, sig_eff[sf]), device=DEV, dtype=DT)
    Um = torch.tensor(Um_np, device=DEV, dtype=DT)
    lam0, xn2 = jbar_lam(XW, Um, S, BETA[claim]); d0 = dprime_avg(lam0, xn2)
    Ur_np, ang = rotate_cols(Um_np, 0.10, seed=7); Ur = torch.tensor(Ur_np, device=DEV, dtype=DT)
    lamr, _ = jbar_lam(XW, Ur, S, BETA[claim]); dr = dprime_avg(lamr, xn2)
    r_tau = T_eff_ratio(16.0, 8.0); fac = float(np.sqrt(r_tau))
    rows.append(dict(cell=f"sf{sf}_{shape}_kw{kwf}_c{claim}", dprime_avg_1pct=round(d0, 2),
                     rot_angle=round(ang, 3), dprime_avg_1pct_rot=round(dr, 2),
                     rot_degradation_pct=round((dr / d0 - 1) * 100, 1),
                     tau2x_Teff_ratio=round(r_tau, 3), tau2x_dprime_factor=round(fac, 3),
                     tau2x_degradation_pct=round((fac - 1) * 100, 1)))
    print(f"{rows[-1]['cell']}: d'avg(1%) {d0:.2f} | +10%rot(ang {ang*100:.0f}%) {dr:.2f} "
          f"({(dr/d0-1)*100:+.0f}%) | tau x2 d'-factor {fac:.3f} ({(fac-1)*100:+.0f}%)", flush=True)
json.dump(dict(robustness=rows), open("DETECTION_ROBUST.json", "w"), indent=2)
print("robustness -> DETECTION_ROBUST.json", flush=True)
