# DETECTION LIVING-REGION MAP -- the demand-axis divergence (pure Fisher arithmetic, NO
# reconstruction). Reuses the validated profiled-Fisher engine (p1_fisher / living_region_map).
#
# Task model: H0 = declared scene x; H1 = x + delta, delta supported STRICTLY beyond the modulator
# band and INSIDE the covariance aperture (delta in U_beta, k_p<|k|<=k_claim<=2k_p). Energy
# eps = ||delta||/||x||. LR test on the lag covariances (Gaussian approx). Detection deflection is
#   d'(T_eff)^2 = T_eff * delta^T Jbar_B delta            (Jbar_B = per-epoch profiled beyond-band
# Fisher; J_B = T_eff*Jbar_B scales exactly linearly, so Jbar_B is T_eff-independent).
# Worst-mode (least-favorable) delta -> lam_min(Jbar_B); energy-spread delta -> lam_mean=tr/ d_beta.
#   T_det(eps) = 25 / (lam * eps^2 * ||x||^2)   (epochs for d'>=5, strong detection).
# Specificity theorem (mean-channel null): dm/d(beta) = P U_beta = 0 exactly (band orthogonality)
#   -> the mean/row channel carries ZERO detection info; detection is covariance-only and specific.
import json, time, numpy as np, torch
DEV = "cuda" if torch.cuda.is_available() else "cpu"; DT = torch.float64
n = 64; N = n * n; KP = 5; M = 128; PHOT = 1e4; TAU = 8.0
t0 = time.time(); _X, _Y = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")

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
BETA_np = {c: band_modes(KP + 1, hi) for c, hi in CLAIMS.items()}
BETA = {c: torch.tensor(v, device=DEV, dtype=DT) for c, v in BETA_np.items()}

def witness_scene(seed=4):
    rs = np.random.default_rng(seed)
    a_in = U_in_np @ rs.standard_normal(U_in_np.shape[1])
    Ube = band_modes(KP + 1, 9); a_be = Ube @ (2.0 * rs.choice([-1, 1], Ube.shape[1]))
    x = a_in + a_be; x -= x.min(); return x / x.max()

def natural_scene(name):
    from skimage import data, transform
    src = {"cameraman": data.camera, "coins": data.coins, "moon": data.moon}[name]()
    img = transform.resize(src.astype(np.float64), (n, n), anti_aliasing=True, mode="reflect")
    return ((img - img.min()) / (np.ptp(img) + 1e-12)).ravel()

SCENES = {"witness": witness_scene(), "cameraman": natural_scene("cameraman"),
          "coins": natural_scene("coins"), "moon": natural_scene("moon")}
SCENES_T = {k: torch.tensor(v, device=DEV, dtype=DT) for k, v in SCENES.items()}

Pabs = torch.abs(Pt)   # nonneg photon-throughput matrix for the physically-correct shot term

def jbar_B(x, Um, S, Phi_beta, T_ref=4096.0, n_ph=PHOT):
    """Per-epoch profiled beyond-band Fisher Jbar_B (=J_B/T_eff), its eigenvalues, and ||x||^2.
    CORRECTED shot: variance = throughput |P|*x scaled to n photons (NOT the near-zero signed Px)."""
    A = Pt * x[None, :]; UmA = Um.t() @ A.t(); KA = Um @ (S[:, None] * UmA); C = A @ KA
    flux = Pabs @ x; R = flux * (flux.mean() / n_ph); V = C + torch.diag(R)   # physical shot
    Vinv = torch.linalg.inv(V + 1e-9 * torch.eye(M, device=DEV, dtype=DT))
    def Gstack(Phi):
        d = Phi.shape[1]; G = torch.empty((d, M, M), device=DEV, dtype=DT)
        for k in range(d):
            B = (Pt * Phi[:, k][None, :]) @ KA; G[k] = Vinv @ (B + B.t())
        return G
    Gb = Gstack(Phi_beta); Ge = Gstack(Phi_eta); Glaw = (Vinv @ C).unsqueeze(0)
    Geta = torch.cat([Ge, Glaw], 0); half = T_ref / 2.0
    Ibb = half * torch.einsum("aij,bji->ab", Gb, Gb)
    Ibe = half * torch.einsum("aij,bji->ab", Gb, Geta)
    Iee = half * torch.einsum("aij,bji->ab", Geta, Geta)
    Min = Pt @ Phi_eta; Iee[:d_eta, :d_eta] += T_ref * (Min.t() @ Vinv @ Min)
    JB = Ibb - Ibe @ torch.linalg.pinv(Iee, rcond=1e-10) @ Ibe.t(); JB = 0.5 * (JB + JB.t())
    lam = torch.clamp(torch.linalg.eigvalsh(JB), min=0).cpu().numpy() / T_ref   # per-epoch
    beta_coeff = (Phi_beta.t() @ x).cpu().numpy()
    return lam, float(np.sum(beta_coeff ** 2)), float((x * x).sum().item())

EPS = [0.005, 0.01, 0.02, 0.05]; T_CAP = 4096

def det_row(lam_bar, xnorm2, T_ref=4096):
    lam_min = float(lam_bar.min()); lam_mean = float(lam_bar.mean())
    d_beta = len(lam_bar)
    out = dict(d_beta=d_beta, lam_min_bar=lam_min, lam_mean_bar=lam_mean,
               dprime_worst={}, dprime_avg={}, Tdet_worst={}, Tdet_avg={})
    for e in EPS:
        e2x = e * e * xnorm2
        out["dprime_worst"][str(e)] = float(np.sqrt(max(T_ref * lam_min * e2x, 0)))
        out["dprime_avg"][str(e)] = float(np.sqrt(max(T_ref * lam_mean * e2x, 0)))
        out["Tdet_worst"][str(e)] = float(25.0 / (lam_min * e2x)) if lam_min * e2x > 0 else float("inf")
        out["Tdet_avg"][str(e)] = float(25.0 / (lam_mean * e2x)) if lam_mean * e2x > 0 else float("inf")
    return out

# ---- specificity theorem: mean-channel null P U_beta = 0 (band orthogonality) ----
spec = {}
for c, Ub in BETA_np.items():
    PU = P_np @ Ub
    spec[str(c)] = dict(frob_norm=float(np.linalg.norm(PU)), max_abs=float(np.abs(PU).max()),
                        rel_to_P=float(np.linalg.norm(PU) / np.linalg.norm(P_np)))
print("=== specificity theorem: mean-channel response to beyond-band delta ===", flush=True)
for c, s in spec.items():
    print(f"  claim {c}k_p: ||P U_beta||_F = {s['frob_norm']:.2e}  max|.| = {s['max_abs']:.2e} "
          f"(rel to ||P|| = {s['rel_to_P']:.2e})  -> mean channel is EXACTLY blind", flush=True)

# ---- effective contrast per sigma_f (reuse bounded-field diagnostic) ----
def eff_contrast(sf, lo, hi, ns=2_000_000, seed=0):
    g = sf * np.random.default_rng(seed).standard_normal(ns); return float(np.std(np.clip(1 + g, lo, hi)))
BOUNDS = {0.3: (0.2, 1.8), 0.6: (0.2, 1.8), 1.0: (0.05, 1.95)}
SIGS = [0.3, 0.6, 1.0]; SHAPES = ["flat", "k^-1", "k^-2"]; KWS = {1: KP, 2: 2 * KP, 4: 4 * KP}
sig_eff = {sf: eff_contrast(sf, *BOUNDS[sf]) for sf in SIGS}
MEDU = {kwf: medium_modes(kw) for kwf, kw in KWS.items()}
MEDU = {kwf: (torch.tensor(U, device=DEV, dtype=DT), ka) for kwf, (U, ka) in MEDU.items()}

# ============================ the 81-cell detection grid ============================
print("\n=== DETECTION grid (81 cells): T_det for d'>=5, witness H0 ===", flush=True)
grid = []
for sf in SIGS:
    for shape in SHAPES:
        for kwf in (1, 2, 4):
            Um, kabs = MEDU[kwf]; S = torch.tensor(power_spectrum(kabs, shape, sig_eff[sf]), device=DEV, dtype=DT)
            for claim in (1.25, 1.5, 1.8):
                lam, bnorm2, xnorm2 = jbar_B(SCENES_T["witness"], Um, S, BETA[claim])
                d = det_row(lam, xnorm2)
                cell = dict(sigma_f=sf, sigma_f_eff=round(sig_eff[sf], 3), shape=shape,
                            k_w_over_kp=kwf, claim=claim, **d)
                grid.append(cell)
                print(f"  sf{sf} {shape:5s} kw{kwf} c{claim}: d'avg(1%)={d['dprime_avg']['0.01']:5.1f} "
                      f"Tdet_avg(1%)={d['Tdet_avg']['0.01']:8.0f} Tdet_worst(1%)={d['Tdet_worst']['0.01']:10.0f}",
                      flush=True)

# ============================ verdict ============================
def frac_le(rows, key, eps="0.01", cap=T_CAP):
    return float(np.mean([r[key][eps] <= cap for r in rows]))
verdict = dict(
    frac_avg_1pct_Tdet_le_4096=frac_le(grid, "Tdet_avg"),
    frac_worst_1pct_Tdet_le_4096=frac_le(grid, "Tdet_worst"),
    frac_avg_2pct_Tdet_le_4096=frac_le(grid, "Tdet_avg", "0.02"),
    frac_avg_0p5pct_Tdet_le_4096=frac_le(grid, "Tdet_avg", "0.005"))
best = sorted(grid, key=lambda r: r["Tdet_avg"]["0.01"])[:5]
region = "DETECTION_REGION_FOUND" if verdict["frac_avg_1pct_Tdet_le_4096"] > 0 else "DETECTION_REGION_EMPTY"
MS_PER_BANK = 12.8   # bank-acquisition: T_eff = # independent banks, each 256 exposures @20kHz
def secs(banks): return round(banks * MS_PER_BANK / 1000.0, 1)
# medium/shot ratio at the best cell (diagnostic of the corrected shot)
_Um, _ka = MEDU[1]; _S = torch.tensor(power_spectrum(_ka, "k^-2", sig_eff[0.3]), device=DEV, dtype=DT)
_A = Pt * SCENES_T["witness"][None, :]; _C = _A @ (_Um @ (_S[:, None] * (_Um.t() @ _A.t())))
_flux = Pabs @ SCENES_T["witness"]; _R = _flux * (_flux.mean() / PHOT)
med_shot_ratio = float((_C.diagonal().mean() / _R.mean()).item())
print(f"\ncorrected medium/shot ratio at pocket cell (k^-2,kw1,c1.25,sf0.3): {med_shot_ratio:.1f}", flush=True)

json.dump(dict(map="detection_living_region_CORRECTED", shot_model="physical: R=|P|x*mean(|P|x)/n",
               medium_shot_ratio_pocket=med_shot_ratio, ms_per_bank=MS_PER_BANK,
               fixed=dict(N=N, M=M, k_p=KP, photons=PHOT, tau=TAU, T_cap=T_CAP, eps=EPS),
               specificity_meanchannel_null=spec, sigma_f_eff=sig_eff,
               grid=grid, verdict_fractions=verdict,
               best_cells=[{k: c[k] for k in ("sigma_f", "shape", "k_w_over_kp", "claim")} |
                           {"Tdet_avg_1pct": round(c["Tdet_avg"]["0.01"], 0),
                            "Tdet_avg_1pct_sec": secs(c["Tdet_avg"]["0.01"]),
                            "dprime_avg_1pct_at4096": round(c["dprime_avg"]["0.01"], 1)} for c in best],
               region_verdict=region, wall_s=time.time() - t0),
          open("DETECTION_MAP_CORRECTED.json", "w"), indent=2)
print(f"\n=== {region} (CORRECTED shot) ===", flush=True)
for e in ["0.005", "0.01", "0.02", "0.05"]:
    fa = frac_le(grid, "Tdet_avg", e); fw = frac_le(grid, "Tdet_worst", e)
    print(f"  eps={float(e)*100:.1f}%: avg {fa*100:.0f}%  worst {fw*100:.0f}%", flush=True)
print(f"  best cell: Tdet_avg(1%)={best[0]['Tdet_avg']['0.01']:.0f} banks = {secs(best[0]['Tdet_avg']['0.01'])} s "
      f"at {best[0]['shape']} kw{best[0]['k_w_over_kp']} c{best[0]['claim']} sf{best[0]['sigma_f']}", flush=True)
print(f"[detection grid CORRECTED done {time.time()-t0:.0f}s]", flush=True)
