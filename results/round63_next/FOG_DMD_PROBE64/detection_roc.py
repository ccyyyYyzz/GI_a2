# DETECTION ROC -- Monte-Carlo validation of the (CORRECTED-shot) analytic d' prognosis.
# Builds the actual lag-covariance LR / profiled matched-score detector and measures empirical
# ROC/AUC + empirical d', compared to the CORRECTED analytic (bug-free MC; corrected target).
# Bank-acquisition model: T_eff = number of independent banks (iid medium states), 12.8 ms/bank.
import json, time, numpy as np, torch
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
DEV = "cuda" if torch.cuda.is_available() else "cpu"; DT = torch.float64
n = 64; N = n * n; KP = 5; M = 128; PHOT = 1e4; MS_PER_BANK = 12.8
_X, _Y = np.meshgrid(np.arange(n), np.arange(n), indexing="ij"); t0 = time.time()

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
                nr = np.linalg.norm(f)
                if nr > 1e-8: cols.append(f / nr); kabs.append(np.hypot(kx, ky))
    return np.array(cols).T, np.array(kabs)
def power_spectrum(kabs, shape, se):
    w = {"flat": np.ones_like(kabs), "k^-1": 1.0 / kabs, "k^-2": 1.0 / kabs ** 2}[shape]
    return w / w.sum() * (N * se ** 2)
def signed_codes(seed=10):
    rp = np.random.default_rng(seed); Uc = band_modes(1, KP)
    P = (rp.standard_normal((M, Uc.shape[1])) @ Uc.T); P /= np.abs(P).max(axis=1, keepdims=True); return P
P_np = signed_codes(); Pt = torch.tensor(P_np, device=DEV, dtype=DT); Pabs = torch.abs(Pt)
U_in_np = band_modes(0, KP); Phi_eta = torch.tensor(U_in_np, device=DEV, dtype=DT); d_eta = U_in_np.shape[1]
CLAIMS = {1.25: 6, 1.5: 7, 1.8: 9}
BETA_np = {c: band_modes(KP + 1, hi) for c, hi in CLAIMS.items()}
BETA = {c: torch.tensor(v, device=DEV, dtype=DT) for c, v in BETA_np.items()}
def witness_scene(seed=4):
    rs = np.random.default_rng(seed); a_in = U_in_np @ rs.standard_normal(U_in_np.shape[1])
    Ube = band_modes(KP + 1, 9); a_be = Ube @ (2.0 * rs.choice([-1, 1], Ube.shape[1]))
    x = a_in + a_be; x -= x.min(); return x / x.max()
XW_np = witness_scene(); XW = torch.tensor(XW_np, device=DEV, dtype=DT)
sig_eff = {0.3: 0.298, 1.0: 0.696}

def setup_cell(sf, shape, kwf, claim):
    """Corrected-shot V0, Vinv, per-mode Vk stack, Fisher blocks; returns builder for the profiled
    matched filter W(c) and the analytic per-epoch profiled Fisher J̄_B."""
    Um_np, kabs = medium_modes({1: KP, 2: 2 * KP, 4: 4 * KP}[kwf])
    S_np = power_spectrum(kabs, shape, sig_eff[sf]); Um = torch.tensor(Um_np, device=DEV, dtype=DT)
    S = torch.tensor(S_np, device=DEV, dtype=DT)
    x = XW; A = Pt * x[None, :]; KA = Um @ (S[:, None] * (Um.t() @ A.t())); C = A @ KA
    flux = Pabs @ x; R = flux * (flux.mean() / PHOT); V0 = C + torch.diag(R)
    Vinv = torch.linalg.inv(V0 + 1e-9 * torch.eye(M, device=DEV, dtype=DT))
    Phi_b = BETA[claim]; db = Phi_b.shape[1]
    def Vks(Phi):
        d = Phi.shape[1]; out = torch.empty((d, M, M), device=DEV, dtype=DT)
        for k in range(d):
            Bk = (Pt * Phi[:, k][None, :]) @ KA; out[k] = Bk + Bk.t()
        return out
    Vk_b = Vks(Phi_b); Vk_e = Vks(Phi_eta); Vk_law = C.unsqueeze(0)   # dV/dlogamp = C
    Vk_eta = torch.cat([Vk_e, Vk_law], 0)
    Gb = torch.einsum("ij,ajk->aik", Vinv, Vk_b); Ge = torch.einsum("ij,ajk->aik", Vinv, Vk_eta)
    Ibb = 0.5 * torch.einsum("aij,bji->ab", Gb, Gb); Ibe = 0.5 * torch.einsum("aij,bji->ab", Gb, Ge)
    Iee = 0.5 * torch.einsum("aij,bji->ab", Ge, Ge)
    Min = Pt @ Phi_eta; Iee[:d_eta, :d_eta] += (Min.t() @ Vinv @ Min)   # per-epoch mean-channel
    JB = Ibb - Ibe @ torch.linalg.pinv(Iee, rcond=1e-10) @ Ibe.t(); JB = 0.5 * (JB + JB.t())
    lamJ = torch.clamp(torch.linalg.eigvalsh(JB), min=0)
    def make_W(c_np):                     # profiled matched filter for beyond-band delta coeffs c
        c = torch.tensor(c_np, device=DEV, dtype=DT)
        a_eta = -torch.linalg.solve(Iee + 1e-12 * torch.eye(Iee.shape[0], device=DEV, dtype=DT), Ibe.t() @ c)
        dV = torch.einsum("a,aij->ij", c, Vk_b) + torch.einsum("a,aij->ij", a_eta, Vk_eta)
        return Vinv @ dV @ Vinv
    Wamp = Vinv @ C @ Vinv                # medium-amplitude matched filter (law direction)
    return dict(Um=Um, S=S, V0=V0, Vinv=Vinv, C=C, R=R, JB=JB, lamJ=lamJ, db=db,
                Phi_b=Phi_b, make_W=make_W, Wamp=Wamp, x=x)

def gen_records(cell, n_rec, T_eff, kind="H0", delta_np=None, S_scale=1.0, rng=None):
    """n_rec records, each T_eff iid banks. kind: H0 / beyond / inband. Returns sample cov + mean per rec."""
    Um, S, V0, C, R, x = cell["Um"], cell["S"] * S_scale, cell["V0"], cell["C"], cell["R"], cell["x"]
    m0 = Pt @ x
    xt = x if delta_np is None else x + torch.tensor(delta_np, device=DEV, dtype=DT)
    A2 = Pt * xt[None, :]; H = (Um.t() @ A2.t()).t()          # (M, db): P diag(x+delta) U_med
    sdS = torch.sqrt(S); mfull = Pt @ xt                      # mean under this scene (beyond: =m0)
    rr = torch.sqrt(R)
    Scov = torch.empty((n_rec, M, M), device=DEV, dtype=DT); Mbar = torch.empty((n_rec, M), device=DEV, dtype=DT)
    g = torch.Generator(device=DEV); g.manual_seed(int(rng))
    for i in range(n_rec):
        Z = torch.randn(T_eff, Um.shape[1], device=DEV, dtype=DT, generator=g)
        B = mfull[None, :] + (Z * sdS[None, :]) @ H.t() + torch.randn(T_eff, M, device=DEV, dtype=DT, generator=g) * rr[None, :]
        Bc = B - B.mean(0, keepdim=True); Scov[i] = (Bc.t() @ Bc) / T_eff; Mbar[i] = B.mean(0)
    return Scov, Mbar, m0

def auc(a, b):   # P(b>a)
    a = np.asarray(a); b = np.asarray(b); return float(np.mean(b[:, None] > a[None, :]))

def dprime_emp(t0v, t1v):
    return float((np.mean(t1v) - np.mean(t0v)) / (np.std(t0v) + 1e-12))

CELLS = [("best", 0.3, "k^-2", 1, 1.25), ("mid", 0.3, "flat", 1, 1.5), ("floor", 0.3, "flat", 1, 1.8)]
EPS = [0.01, 0.02, 0.05]; NREC = 250; T_CAP = 1400
results = {"validation": [], "specificity": {}, "third_class": {}, "robustness": {}}

def energy_spread_delta(cell, eps, seed):
    """delta in the beyond-band subspace, isotropic (energy-spread), ||delta|| = eps*||x||."""
    rng = np.random.default_rng(seed); c = rng.standard_normal(cell["db"]); c = c / np.linalg.norm(c)
    c = c * eps * float(torch.linalg.norm(cell["x"]).item())
    return c, (cell["Phi_b"].cpu().numpy() @ c)   # coeffs, pixel-space delta

for (tag, sf, shape, kwf, claim) in CELLS:
    cell = setup_cell(sf, shape, kwf, claim)
    lam = cell["lamJ"].cpu().numpy(); lam_mean = lam.mean()
    for eps in EPS:
        # analytic (corrected) T_det + choose test T_eff
        xn2 = float((cell["x"] ** 2).sum().item()); e2x = eps * eps * xn2
        Tdet = 25.0 / (lam_mean * e2x); T_eff = int(min(round(Tdet), T_CAP))
        c_np, delta_px = energy_spread_delta(cell, eps, seed=100 + int(eps * 1000))
        dprime_an = float(np.sqrt(T_eff * float(c_np @ (cell["JB"].cpu().numpy() @ c_np)) / 1.0))  # =sqrt(Teff*c^T Jbar c)
        # Note JB stored is per-epoch (built with half=1/2, mean term per-epoch) -> multiply by T_eff
        W = cell["make_W"](c_np)
        V0W = float((cell["V0"] * W).sum().item())
        Sc0, Mb0, m0 = gen_records(cell, NREC, T_eff, "H0", None, 1.0, rng=11 + int(eps * 1e4))
        Sc1, Mb1, _ = gen_records(cell, NREC, T_eff, "beyond", delta_px, 1.0, rng=22 + int(eps * 1e4))
        t0v = (torch.einsum("aij,ij->a", Sc0, W)).cpu().numpy() - V0W
        t1v = (torch.einsum("aij,ij->a", Sc1, W)).cpu().numpy() - V0W
        d_emp = dprime_emp(t0v, t1v); A = auc(t0v, t1v)
        results["validation"].append(dict(cell=tag, eps=eps, T_eff=T_eff, T_det_analytic=round(Tdet, 0),
            dprime_analytic=round(dprime_an, 2), dprime_empirical=round(d_emp, 2),
            ratio=round(d_emp / max(dprime_an, 1e-9), 3), auc=round(A, 4),
            T_det_sec=round(Tdet * MS_PER_BANK / 1000, 1)))
        print(f"[{tag} eps{eps*100:.0f}%] T_eff={T_eff} d'_an={dprime_an:.2f} d'_emp={d_emp:.2f} "
              f"ratio={d_emp/max(dprime_an,1e-9):.2f} AUC={A:.3f}", flush=True)

print(f"[validation block done {time.time()-t0:.0f}s]", flush=True)

# ==================== specificity: two/three-channel discriminator ====================
# Best cell, eps=5%, four populations: H0 / scene-inband / scene-beyond / medium(sigma_f+20%).
cell = setup_cell(0.3, "k^-2", 1, 1.25)
Tsp = 800; eps_sp = 0.05
c_b, delta_b = energy_spread_delta(cell, eps_sp, seed=555)
W_b = cell["make_W"](c_b); V0Wb = float((cell["V0"] * W_b).sum().item())
W_amp = cell["Wamp"]; V0Wa = float((cell["V0"] * W_amp).sum().item())
Vinv = cell["Vinv"]; m0 = Pt @ cell["x"]
# in-band delta (lights the MEAN channel), same energy
rng = np.random.default_rng(556); Uin1 = band_modes(1, KP)
c_in = rng.standard_normal(Uin1.shape[1]); c_in /= np.linalg.norm(c_in)
delta_in = Uin1 @ (c_in * eps_sp * float(torch.linalg.norm(cell["x"]).item()))
def three_stats(Sc, Mb):
    tcov = (torch.einsum("aij,ij->a", Sc, W_b)).cpu().numpy() - V0Wb
    tamp = (torch.einsum("aij,ij->a", Sc, W_amp)).cpu().numpy() - V0Wa
    dm = Mb - m0[None, :]; tmean = (Tsp * torch.einsum("ai,ij,aj->a", dm, Vinv, dm)).cpu().numpy()
    return tmean, tcov, tamp
pops = {}
for kd, (kind, dpx, ssc, sd) in {"H0": ("H0", None, 1.0, 1), "inband": ("inband", delta_in, 1.0, 2),
                                  "beyond": ("beyond", delta_b, 1.0, 3), "medium": ("H0", None, 1.44, 4)}.items():
    Sc, Mb, _ = gen_records(cell, 200, Tsp, kind, dpx, ssc, rng=700 + sd)
    pops[kd] = three_stats(Sc, Mb)
def sep(a, b):  # class-separation d' on a stat: (mean_b - mean_a)/pooled_std
    a = np.asarray(a); b = np.asarray(b); return float((b.mean() - a.mean()) / (0.5 * (a.std() + b.std()) + 1e-12))
results["specificity"] = dict(
    tmean_inband_vs_H0=round(sep(pops["H0"][0], pops["inband"][0]), 1),
    tmean_beyond_vs_H0=round(sep(pops["H0"][0], pops["beyond"][0]), 1),
    tcov_beyond_vs_H0=round(sep(pops["H0"][1], pops["beyond"][1]), 1),
    tcov_inband_vs_H0=round(sep(pops["H0"][1], pops["inband"][1]), 1),
    tcov_medium_vs_H0=round(sep(pops["H0"][1], pops["medium"][1]), 1),
    tamp_medium_vs_H0=round(sep(pops["H0"][2], pops["medium"][2]), 1),
    tamp_beyond_vs_H0=round(sep(pops["H0"][2], pops["beyond"][2]), 1))
print("specificity (class separation d'):", results["specificity"], flush=True)
# scatter plot: t_mean vs t_cov (3 scene populations) + t_amp for medium
fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
col = {"H0": "0.5", "inband": "#1f77b4", "beyond": "#d62728", "medium": "#2ca02c"}
for k in ["H0", "inband", "beyond", "medium"]:
    ax[0].scatter(pops[k][0], pops[k][1], s=10, alpha=0.6, c=col[k], label=k)
ax[0].set_xlabel("mean-channel statistic  t_mean"); ax[0].set_ylabel("beyond-band cov statistic  t_cov")
ax[0].set_title("scene-inband ↔ mean;  scene-beyond ↔ cov"); ax[0].legend(fontsize=8); ax[0].set_xscale("symlog")
for k in ["H0", "beyond", "medium"]:
    ax[1].scatter(pops[k][2], pops[k][1], s=10, alpha=0.6, c=col[k], label=k)
ax[1].set_xlabel("medium-amplitude statistic  t_amp"); ax[1].set_ylabel("beyond-band cov statistic  t_cov")
ax[1].set_title("medium ↔ amplitude (profiled out of t_cov)"); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig("DETECTION_ROC_scatter.png", dpi=130)
print("saved DETECTION_ROC_scatter.png", flush=True)

# ==================== robustness: best-cell ROC under 10% declared-basis rotation ====================
def rotate_cols(U, eps, seed):
    r = np.random.default_rng(seed); R = r.standard_normal(U.shape); R = R - U @ (U.T @ R)
    R = R / (np.linalg.norm(R, axis=0, keepdims=True) + 1e-12)
    Ue, _ = np.linalg.qr(U * np.sqrt(1 - eps ** 2) + R * eps)
    ang = np.arccos(np.clip(np.linalg.svd(U.T @ Ue, compute_uv=False), -1, 1)); return Ue, float(np.sin(ang).mean())
cellb = setup_cell(0.3, "k^-2", 1, 1.25)
Um_np, kabs = medium_modes(KP); Ur_np, ang = rotate_cols(Um_np, 0.10, 7)
# build a DECLARED (rotated-law) cell for the detector filter; DATA generated with the TRUE law
Sdec = power_spectrum(kabs, "k^-2", sig_eff[0.3])
cell_dec = setup_cell(0.3, "k^-2", 1, 1.25)  # true law object (for make_W we override Um below)
# recompute make_W with the rotated declared medium
Ur = torch.tensor(Ur_np, device=DEV, dtype=DT); Sr = torch.tensor(Sdec, device=DEV, dtype=DT)
xb = cellb["x"]; Ab = Pt * xb[None, :]; KAr = Ur @ (Sr[:, None] * (Ur.t() @ Ab.t())); Cr = Ab @ KAr
fluxb = Pabs @ xb; Rr = fluxb * (fluxb.mean() / PHOT); V0r = Cr + torch.diag(Rr)
Vinvr = torch.linalg.inv(V0r + 1e-9 * torch.eye(M, device=DEV, dtype=DT))
def Vks_r(Phi):
    out = torch.empty((Phi.shape[1], M, M), device=DEV, dtype=DT)
    for k in range(Phi.shape[1]):
        Bk = (Pt * Phi[:, k][None, :]) @ KAr; out[k] = Bk + Bk.t()
    return out
Vk_b_r = Vks_r(cellb["Phi_b"]); Vk_e_r = Vks_r(Phi_eta); Vk_eta_r = torch.cat([Vk_e_r, Cr.unsqueeze(0)], 0)
Gb_r = torch.einsum("ij,ajk->aik", Vinvr, Vk_b_r); Ge_r = torch.einsum("ij,ajk->aik", Vinvr, Vk_eta_r)
Ibb_r = 0.5 * torch.einsum("aij,bji->ab", Gb_r, Gb_r); Ibe_r = 0.5 * torch.einsum("aij,bji->ab", Gb_r, Ge_r)
Iee_r = 0.5 * torch.einsum("aij,bji->ab", Ge_r, Ge_r); Minr = Pt @ Phi_eta
Iee_r[:d_eta, :d_eta] += (Minr.t() @ Vinvr @ Minr)
c_r, delta_r = energy_spread_delta(cellb, 0.02, seed=999)
a_eta_r = -torch.linalg.solve(Iee_r + 1e-12 * torch.eye(Iee_r.shape[0], device=DEV, dtype=DT),
                              Ibe_r.t() @ torch.tensor(c_r, device=DEV, dtype=DT))
dV_r = torch.einsum("a,aij->ij", torch.tensor(c_r, device=DEV, dtype=DT), Vk_b_r) + torch.einsum("a,aij->ij", a_eta_r, Vk_eta_r)
W_rot = Vinvr @ dV_r @ Vinvr; V0Wr = float((cellb["V0"] * W_rot).sum().item())   # applied to TRUE-law data cov
W_match = cellb["make_W"](c_r); V0Wm = float((cellb["V0"] * W_match).sum().item())
Tr = 800
Sc0, Mb0, _ = gen_records(cellb, 250, Tr, "H0", None, 1.0, rng=31)
Sc1, Mb1, _ = gen_records(cellb, 250, Tr, "beyond", delta_r, 1.0, rng=32)
for lab, W, V0W in [("matched", W_match, V0Wm), ("rot10", W_rot, V0Wr)]:
    t0v = (torch.einsum("aij,ij->a", Sc0, W)).cpu().numpy() - V0W
    t1v = (torch.einsum("aij,ij->a", Sc1, W)).cpu().numpy() - V0W
    results["robustness"][lab] = dict(auc=round(auc(t0v, t1v), 4), dprime=round(dprime_emp(t0v, t1v), 2))
results["robustness"]["rot_angle"] = round(ang, 3)
print("robustness (best cell, eps2%):", results["robustness"], flush=True)

# ==================== online compute cost ====================
results["online_cost"] = dict(
    covariance_update_flops_per_bank=int(M * M),            # b b^T outer product
    final_score_flops=int(M * M),                           # <Sc, W>
    note="running M x M covariance accumulation: ~M^2=16384 MAC/bank + one M^2 inner product; "
         "the M x M matched filter W is precomputed offline from the declared law.",
    ms_per_bank=MS_PER_BANK)
results["third_class"] = dict(
    note="medium-law change (sigma_f +20%) sits in the law-amplitude nuisance direction, which is "
         "PROFILED OUT of the beyond-band detector -> t_cov(medium)~H0; it separates on the amplitude "
         "statistic t_amp. tau change would separate via lag-decay (not exercised in the iid-bank MC).",
    tcov_medium_vs_H0=results["specificity"]["tcov_medium_vs_H0"],
    tamp_medium_vs_H0=results["specificity"]["tamp_medium_vs_H0"])
json.dump(results, open("DETECTION_ROC.json", "w"), indent=2)
print(f"[all done {time.time()-t0:.0f}s]", flush=True)
