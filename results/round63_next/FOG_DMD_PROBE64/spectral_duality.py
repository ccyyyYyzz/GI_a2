# SPECTRAL DUALITY -- R41 Rank-1 innovation (the PRL lever). Detection-Estimation Spectral Duality.
#
# CLAIM (ruling section 9, Rank 1): estimation is governed by the inverse/lower tail of the efficient
# Fisher spectrum J_B; isotropic (energy-spread) detection by its trace. For d beyond-band modes,
#   E[d'^2]/(T_eff eps^2 ||x||^2) = tr(J_B)/d                                 (detection reads the trace)
#   R_est >= tr(J_B^+)/T_eff                                                  (estimation reads the tail)
#   P := [tr(J_B)/d][tr(J_B^+)/d] >= 1, equality iff isotropic spectrum       (arithmetic-harmonic gap)
# A spiky spectrum is a poor imager and a strong detector BY THEOREM (AM-HM gap), not anecdote.
#
# THIS SCRIPT reuses the CORRECTED physical-|P|x-shot profiled-Fisher engine (identical to
# detection_map.py::jbar_B / living_region_map.py::fisher_metrics) and, for all 81 corrected cells:
#   1. full per-bank efficient Fisher spectrum lam_bar of J_B on the claim annulus; arithmetic mean
#      tr(J)/d, harmonic mean d/tr(J^+), duality product P, effective rank (trJ)^2/tr(J^2),
#      lower-tail mass (scene-weighted [map f_rec convention] AND unweighted) at T_eff=4096.
#   2. estimation verdict (f_rec, NRMSE_CRB) from LIVING_REGION_MAP_CORRECTED.json and detection
#      verdict (T_det(2%) energy-spread) from DETECTION_MAP_CORRECTED.json.
#   3. COLLAPSE TEST: does f_rec collapse onto a function of lower-tail mass, and T_det(2%) onto a
#      function of tr(J)/d, across all 81 cells? (R^2 + Spearman, publication plots.)
#   4. verify P>=1 in every cell; most isotropic / spikiest cells; qualitative high-P prediction.
#   5. honest ledger: COLLAPSE_TIGHT (R^2>=0.7 both) or COLLAPSE_LOOSE.
#
# WRITE-SCOPE: outputs only under FOG_DMD_PROBE64/. No git. NO reconstruction.
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
T_REF = 4096
GAMMA2 = 9.0                     # Fisher SNR >= 3  ->  SNR^2 >= 9 (the map's f_rec convention)
t0 = time.time()
_X, _Y = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")


# ============================ corrected Fisher engine (verbatim geometry) ============================
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
            if max(abs(kx), abs(ky)) == 0:
                continue
            if (kx < 0) or (kx == 0 and ky < 0):
                continue
            ph = 2 * np.pi * (kx * _X + ky * _Y) / n
            for f in (np.cos(ph).ravel(), np.sin(ph).ravel()):
                nrm = np.linalg.norm(f)
                if nrm > 1e-8:
                    cols.append(f / nrm); kabs.append(np.hypot(kx, ky))
    return np.array(cols).T, np.array(kabs)


def power_spectrum(kabs, shape, sig_eff):
    w = {"flat": np.ones_like(kabs), "k^-1": 1.0 / kabs, "k^-2": 1.0 / kabs ** 2}[shape]
    return w / w.sum() * (N * sig_eff ** 2)


def signed_codes(seed=10):
    rp = np.random.default_rng(seed); Uc = band_modes(1, KP)
    P = (rp.standard_normal((M, Uc.shape[1])) @ Uc.T)
    P /= np.abs(P).max(axis=1, keepdims=True)
    return P


P_np = signed_codes()
Pt = torch.tensor(P_np, device=DEV, dtype=DT)
Pabs = torch.abs(Pt)
U_in_np = band_modes(0, KP)
Phi_eta = torch.tensor(U_in_np, device=DEV, dtype=DT)
d_eta = U_in_np.shape[1]
CLAIMS = {1.25: 6, 1.5: 7, 1.8: 9}
BETA = {c: torch.tensor(band_modes(KP + 1, hi), device=DEV, dtype=DT) for c, hi in CLAIMS.items()}


def witness_scene(seed=4):
    rs = np.random.default_rng(seed)
    a_in = U_in_np @ rs.standard_normal(U_in_np.shape[1])
    Ube = band_modes(KP + 1, 9)
    a_be = Ube @ (2.0 * rs.choice([-1, 1], Ube.shape[1]))
    x = a_in + a_be; x -= x.min(); return x / x.max()


X_WIT = torch.tensor(witness_scene(), device=DEV, dtype=DT)


def spectrum(x, Um, S, Phi_beta, T_eff=T_REF, n_ph=PHOT):
    """Full efficient Fisher spectrum of J_B (CORRECTED |P|x shot). Returns per-bank eigenvalues
    lam_bar (=J_B eigenvalues / T_eff, sorted ascending, clamped >=0), the scene's beyond-band
    projection energy per J_B eigenmode a_eig^2 (map f_rec convention), ||beta||^2, ||x||^2."""
    A = Pt * x[None, :]
    UmA = Um.t() @ A.t()
    KA = Um @ (S[:, None] * UmA)
    C = A @ KA
    flux = Pabs @ x
    Rd = flux * (flux.mean() / n_ph)                    # physical complementary-exposure shot
    V = C + torch.diag(Rd)
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
    lam = torch.clamp(evals, min=0)                     # T_eff-scaled eigenvalues of J_B
    beta_coeff = Phi_beta.t() @ x
    a_eig = evecs.t() @ beta_coeff
    lam_bar = (lam / T_eff).cpu().numpy()               # per-bank
    a_eig2 = (a_eig ** 2).cpu().numpy()
    beta_norm2 = float((beta_coeff ** 2).sum().item())
    xnorm2 = float((x * x).sum().item())
    return lam_bar, a_eig2, beta_norm2, xnorm2


# ============================ spectral scalars per cell ============================
def spectral_scalars(lam_bar, a_eig2, beta_norm2, xnorm2, d_beta):
    """AM, HM, duality product P, effective rank, lower-tail mass (scene-weighted + unweighted),
    reconstructed f_rec/NRMSE_CRB (engine self-check)."""
    lam = np.asarray(lam_bar, float)                    # per-bank
    d = int(d_beta)
    pos = lam > 1e-12 * lam.max()
    rank = int(pos.sum())
    lam_p = lam[pos]
    AM = float(lam.sum() / d)                           # tr(J)/d  (all d modes; nulls contribute 0)
    tr_Jplus = float((1.0 / lam_p).sum())               # tr(J^+) on the positive subspace
    tr_Jplus_over_d = tr_Jplus / d
    HM = 1.0 / tr_Jplus_over_d if tr_Jplus_over_d > 0 else 0.0
    P = AM * tr_Jplus_over_d                            # duality product  >= (rank/d)^2 ; = AM/HM
    eff_rank = float(lam.sum() ** 2 / max((lam ** 2).sum(), 1e-300))
    cond = float(lam.max() / lam_p.min()) if rank > 0 else float("inf")
    # lower-tail mass at T_eff=4096, Fisher SNR>=3 threshold (SNR^2>=9)
    lam_JB = lam * T_REF                                # back to J_B scale for the SNR test
    snr2_scene = a_eig2 * lam_JB                        # map's exact f_rec convention
    frac_below_scene = float(np.mean(snr2_scene < GAMMA2))
    f_rec_recon = float(np.mean(snr2_scene >= GAMMA2))
    a_flat2 = beta_norm2 / d                            # unweighted: flat scene energy per mode
    snr2_unw = a_flat2 * lam_JB
    frac_below_unw = float(np.mean(snr2_unw < GAMMA2))
    nmse_crb_recon = float(np.sqrt(tr_Jplus / (T_REF * max(beta_norm2, 1e-12))))
    return dict(d_beta=d, rank=rank, AM_trJ_over_d=AM, HM=HM, trJplus_over_d=tr_Jplus_over_d,
                duality_P=P, eff_rank=eff_rank, eff_rank_norm=eff_rank / d, cond=cond,
                lower_tail_scene=frac_below_scene, lower_tail_unw=frac_below_unw,
                f_rec_recon=f_rec_recon, nmse_crb_recon=nmse_crb_recon)


# ============================ effective contrast (bounded-field) ============================
def eff_contrast(sf, lo, hi, ns=2_000_000, seed=0):
    g = sf * np.random.default_rng(seed).standard_normal(ns)
    return float(np.std(np.clip(1 + g, lo, hi)))


BOUNDS = {0.3: (0.2, 1.8), 0.6: (0.2, 1.8), 1.0: (0.05, 1.95)}
SIGS = [0.3, 0.6, 1.0]
SHAPES = ["flat", "k^-1", "k^-2"]
KWS = {1: KP, 2: 2 * KP, 4: 4 * KP}
sig_eff = {sf: eff_contrast(sf, *BOUNDS[sf]) for sf in SIGS}
MEDU = {kwf: medium_modes(kw) for kwf, kw in KWS.items()}
MEDU = {kwf: (torch.tensor(U, device=DEV, dtype=DT), ka) for kwf, (U, ka) in MEDU.items()}


# ============================ load the corrected verdict maps ============================
HERE = "."
LRM = json.load(open("LIVING_REGION_MAP_CORRECTED.json"))
DET = json.load(open("DETECTION_MAP_CORRECTED.json"))


def key(c):
    return (round(float(c["sigma_f"]), 3), c["shape"], int(c["k_w_over_kp"]), float(c["claim"]))


LRM_BY = {key(c): c for c in LRM["grid"]}
DET_BY = {key(c): c for c in DET["grid"]}


# ============================ the 81-cell spectral sweep ============================
print("=== SPECTRAL DUALITY: full J_B spectrum on all 81 corrected cells ===", flush=True)
rows = []
for sf in SIGS:
    for shape in SHAPES:
        for kwf in (1, 2, 4):
            Um, kabs = MEDU[kwf]
            S = torch.tensor(power_spectrum(kabs, shape, sig_eff[sf]), device=DEV, dtype=DT)
            for claim in (1.25, 1.5, 1.8):
                lam_bar, a_eig2, bnorm2, xnorm2 = spectrum(X_WIT, Um, S, BETA[claim])
                sc = spectral_scalars(lam_bar, a_eig2, bnorm2, xnorm2, len(lam_bar))
                k = (round(sf, 3), shape, kwf, claim)
                le = LRM_BY[k]
                de = DET_BY[k]
                row = dict(sigma_f=sf, sigma_f_eff=round(sig_eff[sf], 3), shape=shape,
                           k_w_over_kp=kwf, claim=claim, **sc,
                           # estimation verdict (LIVING_REGION_MAP_CORRECTED)
                           f_rec=le["f_rec_witness"], f_rec_nat_med=le["f_rec_nat_med"],
                           nmse_crb=le["nmse_crb_witness"], nmse_crb_nat_med=le["nmse_crb_nat_med"],
                           t_req_0p30=le["t_req_0p30"],
                           # detection verdict (DETECTION_MAP_CORRECTED)
                           lam_mean_bar_det=de["lam_mean_bar"], lam_min_bar_det=de["lam_min_bar"],
                           Tdet_avg_2pct=de["Tdet_avg"]["0.02"], Tdet_worst_2pct=de["Tdet_worst"]["0.02"],
                           dprime_avg_2pct=de["dprime_avg"]["0.02"], dprime_worst_2pct=de["dprime_worst"]["0.02"],
                           beta_norm2=bnorm2, xnorm2=xnorm2)
                rows.append(row)
    print(f"  ... sigma_f={sf} done  [{time.time()-t0:.0f}s]", flush=True)

# ---- engine self-check: recomputed f_rec / AM must match the frozen maps ----
frec_err = np.array([abs(r["f_rec_recon"] - r["f_rec"]) for r in rows])
am_relerr = np.array([abs(r["AM_trJ_over_d"] - r["lam_mean_bar_det"]) / max(r["lam_mean_bar_det"], 1e-30)
                      for r in rows])
print(f"\n[self-check] |f_rec_recon - f_rec_map|: max={frec_err.max():.4f} median={np.median(frec_err):.4f}",
      flush=True)
print(f"[self-check] rel |AM - lam_mean_bar(det map)|: max={am_relerr.max():.2e} "
      f"median={np.median(am_relerr):.2e}", flush=True)


# ============================ statistics ============================
def r2_linear(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    b, a = np.polyfit(x, y, 1)
    yh = b * x + a
    ss_res = np.sum((y - yh) ** 2); ss_tot = np.sum((y - y.mean()) ** 2)
    return float(1 - ss_res / ss_tot), float(b), float(a)


def spearman(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    rx -= rx.mean(); ry -= ry.mean()
    return float((rx @ ry) / np.sqrt((rx @ rx) * (ry @ ry)))


f_rec = np.array([r["f_rec"] for r in rows])
nmse = np.array([r["nmse_crb"] for r in rows])
lt_scene = np.array([r["lower_tail_scene"] for r in rows])
lt_unw = np.array([r["lower_tail_unw"] for r in rows])
AM = np.array([r["AM_trJ_over_d"] for r in rows])
HMv = np.array([r["HM"] for r in rows])
trJp = np.array([r["trJplus_over_d"] for r in rows])
lam_min = np.array([r["lam_min_bar_det"] for r in rows])
Pdual = np.array([r["duality_P"] for r in rows])
effr = np.array([r["eff_rank_norm"] for r in rows])
Tavg = np.array([r["Tdet_avg_2pct"] for r in rows])
Tworst = np.array([r["Tdet_worst_2pct"] for r in rows])

# ---- COLLAPSE A (estimation): f_rec vs lower-tail mass ----
r2_frec_unw, *_ = r2_linear(lt_unw, f_rec)
r2_frec_scene, *_ = r2_linear(lt_scene, f_rec)
sp_frec_unw = spearman(lt_unw, f_rec)
sp_frec_scene = spearman(lt_scene, f_rec)
# NRMSE_CRB vs tr(J^+)/d  (CRB analog, log-log)
r2_nmse_trJp, *_ = r2_linear(np.log10(trJp), np.log10(nmse))
sp_nmse_trJp = spearman(trJp, nmse)

# ---- COLLAPSE B (detection): T_det(2%) energy-spread vs tr(J)/d ----
r2_tdet_AM, slope_tdet, icpt_tdet = r2_linear(np.log10(AM), np.log10(Tavg))
sp_tdet_AM = spearman(AM, Tavg)
# discriminant controls: energy-spread detection should NOT read the tail; worst-mode SHOULD
r2_tdet_lammin, *_ = r2_linear(np.log10(lam_min), np.log10(Tavg))
r2_tworst_lammin, *_ = r2_linear(np.log10(lam_min), np.log10(Tworst))
r2_tworst_AM, *_ = r2_linear(np.log10(AM), np.log10(Tworst))
# estimation should NOT read the trace
r2_frec_AM, *_ = r2_linear(AM, f_rec)

# ---- duality inequality ----
P_min_all = float(Pdual.min()); P_max_all = float(Pdual.max())
n_ge1 = int(np.sum(Pdual >= 1.0 - 1e-9))
i_iso = int(np.argmin(Pdual)); i_spiky = int(np.argmax(Pdual))

# ---- qualitative high-P prediction ----
sp_P_frec = spearman(Pdual, f_rec)        # expect negative (spiky -> bad imaging)
sp_P_nmse = spearman(Pdual, nmse)         # expect positive (spiky -> worse NRMSE)
sp_P_Tavg = spearman(Pdual, Tavg)         # ruling: high-P -> good detection -> lower T_det (negative)
sp_effr_frec = spearman(effr, f_rec)      # isotropy (high eff-rank) -> better imaging (positive)

print("\n=== COLLAPSE A (estimation: f_rec vs lower-tail mass) ===", flush=True)
print(f"  f_rec vs lower_tail_UNWEIGHTED : R^2={r2_frec_unw:.3f}  Spearman={sp_frec_unw:+.3f}", flush=True)
print(f"  f_rec vs lower_tail_SCENE(=1-f_rec identity): R^2={r2_frec_scene:.3f}  "
      f"Spearman={sp_frec_scene:+.3f}", flush=True)
print(f"  NRMSE_CRB vs tr(J^+)/d (log-log): R^2={r2_nmse_trJp:.3f}  Spearman={sp_nmse_trJp:+.3f}", flush=True)
print(f"  [control] f_rec vs tr(J)/d (should be POOR): R^2={r2_frec_AM:.3f}", flush=True)
print("\n=== COLLAPSE B (detection: T_det(2%) energy-spread vs tr(J)/d) ===", flush=True)
print(f"  log T_det_avg vs log tr(J)/d : R^2={r2_tdet_AM:.3f}  slope={slope_tdet:+.3f}  "
      f"Spearman={sp_tdet_AM:+.3f}", flush=True)
print(f"  [control] T_det_avg vs lam_min (should be POOR): R^2={r2_tdet_lammin:.3f}", flush=True)
print(f"  [dual]    T_det_WORST vs lam_min (worst-mode reads tail): R^2={r2_tworst_lammin:.3f}", flush=True)
print(f"  [control] T_det_WORST vs tr(J)/d (should be weaker): R^2={r2_tworst_AM:.3f}", flush=True)
print("\n=== DUALITY INEQUALITY  P=[tr(J)/d][tr(J^+)/d] >= 1 ===", flush=True)
print(f"  cells with P>=1: {n_ge1}/81   P range [{P_min_all:.2f}, {P_max_all:.1f}]", flush=True)
print(f"  most isotropic: {rows[i_iso]['shape']} kw{rows[i_iso]['k_w_over_kp']} "
      f"c{rows[i_iso]['claim']} sf{rows[i_iso]['sigma_f']}  P={Pdual[i_iso]:.2f} "
      f"eff_rank/d={effr[i_iso]:.2f} f_rec={f_rec[i_iso]:.2f} Tdet2%={Tavg[i_iso]:.0f}", flush=True)
print(f"  spikiest:       {rows[i_spiky]['shape']} kw{rows[i_spiky]['k_w_over_kp']} "
      f"c{rows[i_spiky]['claim']} sf{rows[i_spiky]['sigma_f']}  P={Pdual[i_spiky]:.1f} "
      f"eff_rank/d={effr[i_spiky]:.2f} f_rec={f_rec[i_spiky]:.2f} Tdet2%={Tavg[i_spiky]:.0f}", flush=True)
print("\n=== qualitative high-P prediction (Spearman) ===", flush=True)
print(f"  P vs f_rec       = {sp_P_frec:+.3f}  (expect <0: spiky -> bad imaging)", flush=True)
print(f"  P vs NRMSE_CRB   = {sp_P_nmse:+.3f}  (expect >0)", flush=True)
print(f"  P vs T_det(2%)   = {sp_P_Tavg:+.3f}  (ruling expect <0: spiky -> strong detector)", flush=True)
print(f"  eff_rank/d vs f_rec = {sp_effr_frec:+.3f}  (expect >0: isotropy -> better imaging)", flush=True)

# ============================ verdict ============================
def key_str(r):
    return f"{r['shape']} kw{r['k_w_over_kp']}x c{r['claim']} sf{r['sigma_f']}"


TIGHT = (r2_frec_unw >= 0.7) and (r2_tdet_AM >= 0.7)
verdict = "COLLAPSE_TIGHT" if TIGHT else "COLLAPSE_LOOSE"
print(f"\n=== SPECTRAL-DUALITY VERDICT: {verdict} "
      f"(estimation R^2={r2_frec_unw:.2f}, detection R^2={r2_tdet_AM:.2f}) ===", flush=True)

stats = dict(
    collapse_estimation=dict(
        frec_vs_lowertail_unweighted=dict(R2=r2_frec_unw, spearman=sp_frec_unw),
        frec_vs_lowertail_scene_identity=dict(R2=r2_frec_scene, spearman=sp_frec_scene,
            note="scene-weighted lower-tail mass == 1 - f_rec by construction (map convention); "
                 "reported as an identity/consistency check, not an independent collapse"),
        nmse_crb_vs_trJplus_loglog=dict(R2=r2_nmse_trJp, spearman=sp_nmse_trJp),
        control_frec_vs_trace=dict(R2=r2_frec_AM,
            note="estimation should NOT collapse on the trace")),
    collapse_detection=dict(
        Tdet_avg_2pct_vs_trace_loglog=dict(R2=r2_tdet_AM, slope=slope_tdet, spearman=sp_tdet_AM),
        control_Tdet_avg_vs_lammin=dict(R2=r2_tdet_lammin,
            note="energy-spread detection should NOT read the worst mode"),
        dual_Tdet_worst_vs_lammin=dict(R2=r2_tworst_lammin,
            note="worst-mode detection DOES read the lower tail"),
        control_Tdet_worst_vs_trace=dict(R2=r2_tworst_AM)),
    duality_inequality=dict(n_cells_P_ge_1=n_ge1, P_min=P_min_all, P_max=P_max_all,
        most_isotropic=key_str(rows[i_iso]), spikiest=key_str(rows[i_spiky]),
        P_most_isotropic=float(Pdual[i_iso]), P_spikiest=float(Pdual[i_spiky])),
    qualitative_highP=dict(spearman_P_frec=sp_P_frec, spearman_P_nmse=sp_P_nmse,
        spearman_P_Tdet_avg_2pct=sp_P_Tavg, spearman_effrank_frec=sp_effr_frec),
    engine_selfcheck=dict(frec_abs_err_max=float(frec_err.max()),
        frec_abs_err_median=float(np.median(frec_err)),
        AM_vs_detmap_relerr_max=float(am_relerr.max())))

out = dict(analysis="detection_estimation_spectral_duality", ruling="R41 Rank 1",
           shot_model="physical |P|x complementary-exposure (CORRECTED engine)",
           fixed=dict(N=N, M=M, k_p=KP, photons=PHOT, T_eff=T_REF, gamma=3.0),
           lower_tail_definition=dict(
               scene_weighted="frac_i( a_i^2 * lam_JB,i < 9 ), a_i^2=scene proj energy on J_B "
                              "eigenmode i (== 1 - f_rec, map convention)",
               unweighted="frac_i( (||beta||^2/d) * lam_JB,i < 9 ), flat energy-spread per mode"),
           verdict=verdict, statistics=stats, n_cells=len(rows), grid=rows,
           wall_s=time.time() - t0)
json.dump(out, open("SPECTRAL_DUALITY.json", "w"), indent=2)
print(f"[wrote SPECTRAL_DUALITY.json]  [{time.time()-t0:.0f}s]", flush=True)


# ============================ publication plots ============================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans", "axes.linewidth": 0.8,
                     "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
                     "figure.dpi": 150})
# Okabe-Ito colorblind-safe palette, one colour per claim shell
CLAIM_COL = {1.25: "#0072B2", 1.5: "#E69F00", 1.8: "#009E73"}
KW_MARK = {1: "o", 2: "s", 4: "^"}
claims = np.array([r["claim"] for r in rows])
kws = np.array([r["k_w_over_kp"] for r in rows])


def scatter_by_group(ax, x, y):
    for c in (1.25, 1.5, 1.8):
        for kw in (1, 2, 4):
            m = (claims == c) & (kws == kw)
            ax.scatter(x[m], y[m], c=CLAIM_COL[c], marker=KW_MARK[kw], s=34,
                       edgecolors="black", linewidths=0.4, alpha=0.9, zorder=3)


def legend_handles():
    h = [Line2D([0], [0], marker="o", color="w", markerfacecolor=CLAIM_COL[c],
                markeredgecolor="k", markersize=7, label=f"claim {c}$\\,k_p$") for c in (1.25, 1.5, 1.8)]
    h += [Line2D([0], [0], marker=KW_MARK[kw], color="w", markerfacecolor="0.6",
                 markeredgecolor="k", markersize=7, label=f"$k_w/k_p={kw}$") for kw in (1, 2, 4)]
    return h


# ---- Figure: estimation collapse  f_rec vs lower-tail mass ----
fig, ax = plt.subplots(figsize=(4.6, 4.0))
scatter_by_group(ax, lt_unw, f_rec)
xs = np.linspace(lt_unw.min(), lt_unw.max(), 100)
b, a = np.polyfit(lt_unw, f_rec, 1)
ax.plot(xs, b * xs + a, "-", color="0.25", lw=1.4, zorder=2,
        label=f"linear fit  $R^2$={r2_frec_unw:.2f}")
ax.set_xlabel("lower-tail mass  (frac. modes below Fisher SNR 3, energy-spread)")
ax.set_ylabel("$f_{\\mathrm{rec}}$  (frac. beyond-band modes recoverable)")
ax.set_title("Estimation does NOT collapse onto the spectral lower tail", fontsize=9.5)
ax.text(0.03, 0.06, f"$R^2$={r2_frec_unw:.2f}\nSpearman={sp_frec_unw:+.2f}",
        transform=ax.transAxes, va="bottom", ha="left", fontsize=8.5,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", lw=0.6))
ax.legend(handles=legend_handles() + [Line2D([0], [0], color="0.25", lw=1.4, label=f"fit $R^2$={r2_frec_unw:.2f}")],
          fontsize=6.6, loc="upper right", framealpha=0.9, ncol=2)
ax.grid(True, ls=":", lw=0.5, alpha=0.6)
fig.tight_layout()
fig.savefig("duality_collapse_frec.pdf"); fig.savefig("duality_collapse_frec.png", dpi=300)
plt.close(fig)

# ---- Figure: detection collapse  T_det(2%) vs tr(J)/d ----
fig, ax = plt.subplots(figsize=(4.6, 4.0))
scatter_by_group(ax, AM, Tavg)
ax.set_xscale("log"); ax.set_yscale("log")
xs = np.logspace(np.log10(AM.min()), np.log10(AM.max()), 100)
ax.plot(xs, 10 ** icpt_tdet * xs ** slope_tdet, "-", color="0.25", lw=1.4, zorder=2)
ax.axhline(T_REF, color="#D55E00", ls="--", lw=1.0, zorder=1)
ax.text(AM.min(), T_REF * 1.15, "4096-bank cap", color="#D55E00", fontsize=7.5, va="bottom")
ax.set_xlabel("tr$(J_B)/d$  (per-bank, energy-spread Fisher supply)")
ax.set_ylabel("$T_{\\mathrm{det}}$ (banks) for $d'\\!\\geq\\!5$,  $\\varepsilon=2\\%$")
ax.set_title("Isotropic detection collapses onto the trace", fontsize=9.5)
ax.text(0.03, 0.06, f"$R^2$={r2_tdet_AM:.3f}\nslope={slope_tdet:+.2f}\nSpearman={sp_tdet_AM:+.2f}",
        transform=ax.transAxes, va="bottom", ha="left", fontsize=8.5,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", lw=0.6))
ax.legend(handles=legend_handles(), fontsize=6.6, loc="upper right", framealpha=0.9, ncol=2)
ax.grid(True, which="both", ls=":", lw=0.5, alpha=0.6)
fig.tight_layout()
fig.savefig("duality_collapse_tdet.pdf"); fig.savefig("duality_collapse_tdet.png", dpi=300)
plt.close(fig)
print(f"[wrote duality_collapse_frec/tdet .pdf/.png]  [{time.time()-t0:.0f}s]", flush=True)


# ============================ markdown report ============================
def fnum(v, p=3):
    return f"{v:.{p}f}"


order = np.argsort(Pdual)
lines = []
A = lines.append
A("# Detection-Estimation Spectral Duality (R41 Rank-1, the PRL lever)\n")
A(f"**Verdict: {verdict}** — estimation collapse $R^2$={fnum(r2_frec_unw,2)}, "
  f"detection collapse $R^2$={fnum(r2_tdet_AM,3)}.\n")
A("Corrected physical-|P|x-shot profiled Fisher engine, all 81 cells, witness scene, "
  "$T_{\\mathrm{eff}}$=4096, $n$=10$^4$ e$^-$/bucket. No reconstruction.\n")
A("## The claim\n")
A("Estimation is governed by the **inverse/lower tail** of the efficient Fisher spectrum $J_B$; "
  "isotropic (energy-spread) detection by its **trace**:\n")
A("```\nE[d'^2]/(T_eff eps^2 ||x||^2) = tr(J_B)/d        (detection reads the trace)\n"
  "R_est >= tr(J_B^+)/T_eff                          (estimation reads the tail)\n"
  "P := [tr(J_B)/d][tr(J_B^+)/d] >= 1,  = iff isotropic\n```\n")
A("## 1. Engine self-check\n")
A(f"- recomputed $f_{{rec}}$ vs frozen map: max abs err {fnum(frec_err.max(),4)}, "
  f"median {fnum(np.median(frec_err),4)}.\n")
A(f"- recomputed tr$(J)/d$ vs detection-map `lam_mean_bar`: max rel err {am_relerr.max():.1e}. "
  "The spectrum is the same object that produced both frozen maps.\n")
A("\n## 2. The two collapses\n")
A("| collapse | predictor | $R^2$ | Spearman |\n|---|---|---:|---:|\n")
A(f"| **Estimation** $f_{{rec}}$ | lower-tail mass (unweighted/energy-spread) | "
  f"{fnum(r2_frec_unw,3)} | {sp_frec_unw:+.3f} |\n")
A(f"| Estimation NRMSE$_{{CRB}}$ | tr$(J^+)/d$ (log-log) | {fnum(r2_nmse_trJp,3)} | {sp_nmse_trJp:+.3f} |\n")
A(f"| **Detection** $T_{{det}}(2\\%)$ energy-spread | tr$(J)/d$ (log-log) | "
  f"{fnum(r2_tdet_AM,3)} | {sp_tdet_AM:+.3f} |\n")
A(f"| Detection $T_{{det}}(2\\%)$ worst-mode | $\\lambda_{{min}}$ (log-log, dual) | "
  f"{fnum(r2_tworst_lammin,3)} | — |\n")
A("\n**Discriminant controls** (the axes do not cross-talk):\n\n")
A("| test | $R^2$ | reads |\n|---|---:|---|\n")
A(f"| $f_{{rec}}$ vs tr$(J)/d$ | {fnum(r2_frec_AM,3)} | estimation does NOT read the trace |\n")
A(f"| $T_{{det}}$avg vs $\\lambda_{{min}}$ | {fnum(r2_tdet_lammin,3)} | isotropic detection does NOT read the tail |\n")
A(f"| $T_{{det}}$worst vs tr$(J)/d$ | {fnum(r2_tworst_AM,3)} | worst-mode detection does NOT read the trace |\n")
A(f"\nScene-weighted lower-tail mass is $\\equiv 1-f_{{rec}}$ by construction "
  f"($R^2$={fnum(r2_frec_scene,3)}); reported as an identity check, not an independent collapse.\n")
A("\n## 3. Duality inequality  P = [tr(J)/d][tr(J^+)/d] >= 1\n")
A(f"- Holds in **{n_ge1}/81** cells. Range $P\\in[{fnum(P_min_all,2)},\\,{fnum(P_max_all,1)}]$; "
  "all spectra full-rank (no exact null), so the AM-HM gap is exact.\n")
A(f"- **Most isotropic** ($P$={fnum(Pdual[i_iso],2)}): {key_str(rows[i_iso])} — "
  f"eff-rank/d={fnum(effr[i_iso],2)}, $f_{{rec}}$={fnum(f_rec[i_iso],2)}, "
  f"$T_{{det}}(2\\%)$={Tavg[i_iso]:.0f} banks.\n")
A(f"- **Spikiest** ($P$={fnum(Pdual[i_spiky],1)}): {key_str(rows[i_spiky])} — "
  f"eff-rank/d={fnum(effr[i_spiky],2)}, $f_{{rec}}$={fnum(f_rec[i_spiky],2)}, "
  f"$T_{{det}}(2\\%)$={Tavg[i_spiky]:.0f} banks.\n")
A("\n## 4. Qualitative high-P prediction (Spearman across 81 cells)\n")
A("| relation | Spearman | prediction |\n|---|---:|---|\n")
A(f"| $P$ vs $f_{{rec}}$ | {sp_P_frec:+.3f} | spiky -> worse imaging (expect <0) |\n")
A(f"| $P$ vs NRMSE$_{{CRB}}$ | {sp_P_nmse:+.3f} | spiky -> worse NRMSE (expect >0) |\n")
A(f"| $P$ vs $T_{{det}}(2\\%)$ | {sp_P_Tavg:+.3f} | ruling expected <0 (spiky=strong detector) — "
  "**CONTRADICTED**, sign is positive |\n")
A(f"| eff-rank/d vs $f_{{rec}}$ | {sp_effr_frec:+.3f} | isotropy -> better imaging (expect >0) |\n")
A("\n## Honest ledger — does (9.1) earn the central-theorem slot?\n")
A("**No. It stays a lemma.** The ruling's promotion test is whether the maps *collapse onto one "
  "spectral-anisotropy curve*. They do not, for three concrete reasons:\n\n")
A(f"1. **The two exact collapses are definitional, not predictive.** $T_{{det}}(2\\%)$ vs tr$(J)/d$ "
  f"($R^2$={fnum(r2_tdet_AM,3)}) and NRMSE$_{{CRB}}$ vs tr$(J^+)/d$ ($R^2$={fnum(r2_nmse_trJp,3)}) are "
  "both perfect because the detection deflection and the CRB are *built from these very spectral "
  "functionals* (eq 3.3 and the profiled CRB), with $\\|x\\|^2$ fixed across cells. They confirm the "
  "algebra of the duality but are not independent evidence. The independent detection check is the "
  "sealed Monte-Carlo probe (R41 section 4), not this analytic map.\n")
A(f"2. **The one genuinely non-circular estimation collapse is loose.** $f_{{rec}}$ — a "
  "scene-weighted *count* of illuminated modes — does not collapse onto any scene-independent "
  f"lower-tail mass ($R^2$={fnum(r2_frec_unw,3)}, Spearman={sp_frec_unw:+.2f}). $f_{{rec}}$ is not a "
  "spectral invariant; it depends on which modes the witness happens to occupy. The spectral "
  "estimation quantity that *is* invariant is tr$(J^+)$, and that collapse is the definitional one "
  "in point 1.\n")
A(f"3. **The anisotropy product $P$ does not deliver the promised task trade-off.** $P$ correctly "
  f"anti-predicts imaging ($P$ vs $f_{{rec}}$ Spearman {sp_P_frec:+.2f}; $P$ vs NRMSE {sp_P_nmse:+.2f}), "
  f"but has the **wrong sign for detection** ($P$ vs $T_{{det}}(2\\%)$ = {sp_P_Tavg:+.2f}; the ruling "
  "predicted spiky $\\to$ *stronger* detector, i.e. negative). Concretely the most isotropic cell "
  f"detects in {Tavg[i_iso]:.0f} banks while the spikiest needs {Tavg[i_spiky]:.0f} — spiky cells are "
  "worse at **both** tasks. The reason: across this physical grid, anisotropy and total Fisher "
  "supply co-vary — spiky cells also have a smaller trace — so 'spiky = strong detector' (true only "
  "at *fixed* trace) never materialises.\n\n")
A("**Consequence for the flagship.** The dead-estimation / live-detection dissociation is real, but "
  "it is driven by the **demand threshold** (imaging requires every mode above SNR 3; detection "
  "pools the trace), not by spectral anisotropy $P$. Eq (9.1) is an exact per-cell AM-HM identity "
  "worth one lemma/remark; it is not the paper's organising principle and is not a PRL-grade central "
  "object on this evidence. Promote the demand-threshold mechanism, keep (9.1) as support.\n")
A("\n## 5. Per-cell table (sorted by duality product P, spiky -> isotropic)\n")
A("| cell | $d$ | tr$(J)/d$ | $\\lambda_{min}$ | $P$ | eff-rank/d | low-tail | "
  "$f_{rec}$ | NRMSE | $T_{det}2\\%$avg |\n|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|\n")
for i in order[::-1]:
    r = rows[i]
    A(f"| {key_str(r)} | {r['d_beta']} | {r['AM_trJ_over_d']:.4f} | {r['lam_min_bar_det']:.4f} | "
      f"{r['duality_P']:.1f} | {r['eff_rank_norm']:.2f} | {r['lower_tail_unw']:.2f} | "
      f"{r['f_rec']:.2f} | {r['nmse_crb']:.2f} | {r['Tdet_avg_2pct']:.0f} |\n")
A("\n## Figures\n")
A("- `duality_collapse_frec.pdf/png` — estimation $f_{rec}$ vs lower-tail mass.\n")
A("- `duality_collapse_tdet.pdf/png` — detection $T_{det}(2\\%)$ vs tr$(J)/d$ (log-log).\n")
open("SPECTRAL_DUALITY.md", "w", encoding="utf-8").writelines(lines)
print(f"[wrote SPECTRAL_DUALITY.md]  total wall {time.time()-t0:.0f}s", flush=True)

