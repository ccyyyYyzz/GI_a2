# BOLD PUSH -- three attacks on the detection sentinel's limits (corrected-shot engine, MC-validated).
#  A1: detection-optimal code design (edge-concentrated / greedy spectral profile vs flat).
#  A2: 0.5% floor assault (optimal codes + longer T <=16384 + matched-target priors) -> eps_min(T).
#  A3: CUSUM streaming sentinel -> detection latency vs eps at ARL0 = 1/hour, 1/day.
import json, time, numpy as np, torch
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
DEV = "cuda" if torch.cuda.is_available() else "cpu"; DT = torch.float64
n = 64; N = n * n; KP = 5; M = 128; PHOT = 1e4; MS_PER_BANK = 12.8
_X, _Y = np.meshgrid(np.arange(n), np.arange(n), indexing="ij"); t0 = time.time()
BANKS_PER_HOUR = 3600.0 / (MS_PER_BANK / 1000.0)      # 281250
BANKS_PER_DAY = 24 * BANKS_PER_HOUR

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

Uc_np, kc = medium_modes(KP)                            # code Fourier modes (|k|<=k_p) + radial freq
def codes_with_profile(w_of_k, seed=10):
    """M signed band-limited codes with spectral weight profile w(|k|) over the in-band modes.
    Physical nonneg realization is 0.5+0.45*P/max|P|; the signed P (zero-mean) enters the Fisher."""
    rp = np.random.default_rng(seed); wsc = np.sqrt(np.asarray(w_of_k))
    C = rp.standard_normal((M, Uc_np.shape[1])) * wsc[None, :]
    P = C @ Uc_np.T; P /= np.abs(P).max(axis=1, keepdims=True); return P

U_in_np = band_modes(0, KP); Phi_eta = torch.tensor(U_in_np, device=DEV, dtype=DT); d_eta = U_in_np.shape[1]
CLAIMS = {1.25: 6, 1.5: 7, 1.8: 9}
BETA = {c: torch.tensor(band_modes(KP + 1, hi), device=DEV, dtype=DT) for c, hi in CLAIMS.items()}
def witness_scene(seed=4):
    rs = np.random.default_rng(seed); a_in = U_in_np @ rs.standard_normal(U_in_np.shape[1])
    Ube = band_modes(KP + 1, 9); a_be = Ube @ (2.0 * rs.choice([-1, 1], Ube.shape[1]))
    x = a_in + a_be; x -= x.min(); return x / x.max()
XW = torch.tensor(witness_scene(), device=DEV, dtype=DT)
sig_eff = {0.3: 0.298, 0.6: 0.503, 1.0: 0.696}

def setup_cell(P_np, sf, shape, kwf, claim, want_filter=False):
    """Corrected-shot profiled beyond-band Fisher for a given code bank P. Returns eigen-summary
    (lam_mean/min, tr), ||x||^2, and (optionally) V0/Vinv/Vk for the MC detector."""
    Pt = torch.tensor(P_np, device=DEV, dtype=DT); Pabs = torch.abs(Pt)
    Um_np, kabs = medium_modes({1: KP, 2: 2 * KP, 4: 4 * KP}[kwf])
    Um = torch.tensor(Um_np, device=DEV, dtype=DT)
    S = torch.tensor(power_spectrum(kabs, shape, sig_eff[sf]), device=DEV, dtype=DT)
    x = XW; A = Pt * x[None, :]; KA = Um @ (S[:, None] * (Um.t() @ A.t())); C = A @ KA
    flux = Pabs @ x; R = flux * (flux.mean() / PHOT); V0 = C + torch.diag(R)
    Vinv = torch.linalg.inv(V0 + 1e-9 * torch.eye(M, device=DEV, dtype=DT))
    Phi_b = BETA[claim]; db = Phi_b.shape[1]
    def Vks(Phi):
        out = torch.empty((Phi.shape[1], M, M), device=DEV, dtype=DT)
        for k in range(Phi.shape[1]):
            Bk = (Pt * Phi[:, k][None, :]) @ KA; out[k] = Bk + Bk.t()
        return out
    Vk_b = Vks(Phi_b); Vk_e = Vks(Phi_eta); Vk_eta = torch.cat([Vk_e, C.unsqueeze(0)], 0)
    Gb = torch.einsum("ij,ajk->aik", Vinv, Vk_b); Ge = torch.einsum("ij,ajk->aik", Vinv, Vk_eta)
    Ibb = 0.5 * torch.einsum("aij,bji->ab", Gb, Gb); Ibe = 0.5 * torch.einsum("aij,bji->ab", Gb, Ge)
    Iee = 0.5 * torch.einsum("aij,bji->ab", Ge, Ge); Min = Pt @ Phi_eta
    Iee[:d_eta, :d_eta] += (Min.t() @ Vinv @ Min)
    JB = Ibb - Ibe @ torch.linalg.pinv(Iee, rcond=1e-10) @ Ibe.t(); JB = 0.5 * (JB + JB.t())
    lam = torch.clamp(torch.linalg.eigvalsh(JB), min=0)
    out = dict(lam_mean=float(lam.mean().item()), lam_min=float(lam.min().item()),
               lam_max=float(lam.max().item()), tr=float(lam.sum().item()), db=db,
               xnorm2=float((x * x).sum().item()))
    if want_filter:
        def make_W(c_np):
            c = torch.tensor(c_np, device=DEV, dtype=DT)
            a_eta = -torch.linalg.solve(Iee + 1e-12 * torch.eye(Iee.shape[0], device=DEV, dtype=DT), Ibe.t() @ c)
            dV = torch.einsum("a,aij->ij", c, Vk_b) + torch.einsum("a,aij->ij", a_eta, Vk_eta)
            return Vinv @ dV @ Vinv
        out.update(dict(Pt=Pt, Pabs=Pabs, Um=Um, S=S, V0=V0, Vinv=Vinv, JB=JB, Phi_b=Phi_b,
                        make_W=make_W, x=x, R=R))
    return out

BEST = (0.3, "k^-2", 1, 1.25); MED = (0.3, "flat", 1, 1.5)
RES = {}
def save(): json.dump(RES, open("BOLD_PUSH.json", "w"), indent=2)

# ============================ ATTACK 1: detection-optimal code design ============================
print("=== ATTACK 1: detection-optimal code design ===", flush=True)
def tr_for_profile(w_of_k, cell_spec, nseed=3):
    trs, lmins, lmeans = [], [], []
    for s in range(nseed):
        P = codes_with_profile(w_of_k(kc), seed=10 + s)
        r = setup_cell(P, *cell_spec)
        trs.append(r["tr"]); lmins.append(r["lam_min"]); lmeans.append(r["lam_mean"])
    return float(np.mean(trs)), float(np.mean(lmins)), float(np.mean(lmeans))

PROFILES = {"flat": lambda k: np.ones_like(k),
            "k^1": lambda k: k, "k^2": lambda k: k ** 2, "k^3": lambda k: k ** 3,
            "edge_exp1": lambda k: np.exp(1.0 * k), "edge_exp2": lambda k: np.exp(2.0 * k),
            "low_k^-2": lambda k: 1.0 / k ** 2}
RES["attack1"] = {"best_cell": BEST, "median_cell": MED, "profiles": {}}
for cellname, spec in [("best", BEST), ("median", MED)]:
    RES["attack1"]["profiles"][cellname] = {}
    base_tr, base_lmin, base_lmean = tr_for_profile(PROFILES["flat"], spec)
    for pname, pf in PROFILES.items():
        tr, lmin, lmean = tr_for_profile(pf, spec)
        RES["attack1"]["profiles"][cellname][pname] = dict(
            tr=round(tr, 6), lam_min=round(lmin, 8), lam_mean=round(lmean, 8),
            dprime_gain_agg=round(np.sqrt(tr / base_tr), 3), dprime_gain_worst=round(np.sqrt(max(lmin, 0) / max(base_lmin, 1e-30)), 3))
    save()
    best_agg = max(RES["attack1"]["profiles"][cellname], key=lambda p: RES["attack1"]["profiles"][cellname][p]["tr"])
    print(f"  [{cellname}] best-aggregate profile: {best_agg} "
          f"gain d'x{RES['attack1']['profiles'][cellname][best_agg]['dprime_gain_agg']} "
          f"(worst-mode x{RES['attack1']['profiles'][cellname][best_agg]['dprime_gain_worst']})", flush=True)

# greedy per-shell refinement on the best cell (5 Chebyshev shells)
shell = np.clip(np.round(kc).astype(int), 1, KP)
def shellw_profile(wsh): return lambda k: wsh[np.clip(np.round(k).astype(int), 1, KP) - 1]
w = np.ones(KP)
cur_tr, _, _ = tr_for_profile(shellw_profile(w), BEST)
for it in range(12):
    improved = False
    for s in range(KP):
        for fac in (1.5, 1 / 1.5):
            wt = w.copy(); wt[s] *= fac; wt = wt / wt.sum() * KP
            tr, _, _ = tr_for_profile(shellw_profile(wt), BEST, nseed=2)
            if tr > cur_tr * 1.001:
                w, cur_tr, improved = wt, tr, True
    if not improved: break
base_tr_best = RES["attack1"]["profiles"]["best"]["flat"]["tr"]
RES["attack1"]["greedy_best"] = dict(shell_weights=[round(float(x), 3) for x in (w / w.sum() * KP)],
                                     tr=round(cur_tr, 6), dprime_gain_agg=round(np.sqrt(cur_tr / base_tr_best), 3))
save()
print(f"  greedy per-shell (best cell): weights={RES['attack1']['greedy_best']['shell_weights']} "
      f"gain d'x{RES['attack1']['greedy_best']['dprime_gain_agg']}", flush=True)
# pick the OPTIMAL profile (best of parametric + greedy) for downstream attacks
opt_name = max(PROFILES, key=lambda p: RES["attack1"]["profiles"]["best"][p]["tr"])
OPT_PROFILE = shellw_profile(w) if RES["attack1"]["greedy_best"]["tr"] > RES["attack1"]["profiles"]["best"][opt_name]["tr"] else PROFILES[opt_name]
P_OPT = codes_with_profile(OPT_PROFILE(kc), seed=10)
P_FLAT = codes_with_profile(PROFILES["flat"](kc), seed=10)
print(f"[attack1 done {time.time()-t0:.0f}s]", flush=True)

# ============================ ATTACK 2: the 0.5% floor assault ============================
print("\n=== ATTACK 2: 0.5% floor assault -- frontier eps_min(T) ===", flush=True)
cellsA2 = {"best": setup_cell(P_OPT, *BEST, want_filter=True), "median": setup_cell(P_FLAT, *MED, want_filter=True)}
Tgrid = [256, 512, 1024, 2048, 4096, 8192, 16384]
RES["attack2"] = {"note": "eps_min(T)=5/(||x|| sqrt(T*lam)); banks independent so multi-lag GLS is moot",
                  "T_grid": Tgrid, "frontier": {}}
for nm, cell in cellsA2.items():
    xn = np.sqrt(cell["xnorm2"]); fr = {"energy_spread": [], "matched_target": []}
    for T in Tgrid:
        fr["energy_spread"].append(round(5.0 / (xn * np.sqrt(T * cell["lam_mean"])), 5))
        fr["matched_target"].append(round(5.0 / (xn * np.sqrt(T * cell["lam_max"])), 5))
    RES["attack2"]["frontier"][nm] = fr
    def T_for_eps(eps, lam): return 25.0 / (eps * eps * cell["xnorm2"] * lam)
    RES["attack2"]["frontier"][nm + "_T_for_0.5pct"] = dict(
        energy_spread=round(T_for_eps(0.005, cell["lam_mean"]), 0),
        matched_target=round(T_for_eps(0.005, cell["lam_max"]), 0),
        energy_spread_sec=round(T_for_eps(0.005, cell["lam_mean"]) * MS_PER_BANK / 1000, 1),
        matched_target_sec=round(T_for_eps(0.005, cell["lam_max"]) * MS_PER_BANK / 1000, 1))
    print(f"  [{nm}] eps_min(16384): spread {fr['energy_spread'][-1]*100:.2f}% matched {fr['matched_target'][-1]*100:.2f}% ; "
          f"0.5% at T_spread={RES['attack2']['frontier'][nm+'_T_for_0.5pct']['energy_spread']:.0f} banks "
          f"({RES['attack2']['frontier'][nm+'_T_for_0.5pct']['energy_spread_sec']:.0f}s) "
          f"T_matched={RES['attack2']['frontier'][nm+'_T_for_0.5pct']['matched_target']:.0f} banks", flush=True)
save()

def gen_S(cell, n_rec, T_eff, delta_np, seed):
    Um, S, x, R = cell["Um"], cell["S"], cell["x"], cell["R"]
    xt = x if delta_np is None else x + torch.tensor(delta_np, device=DEV, dtype=DT)
    A2 = cell["Pt"] * xt[None, :]; H = (Um.t() @ A2.t()).t(); sdS = torch.sqrt(S); mfull = cell["Pt"] @ xt
    rr = torch.sqrt(R); g = torch.Generator(device=DEV); g.manual_seed(int(seed))
    Sc = torch.empty((n_rec, M, M), device=DEV, dtype=DT)
    for i in range(n_rec):
        Z = torch.randn(T_eff, Um.shape[1], device=DEV, dtype=DT, generator=g)
        B = mfull[None, :] + (Z * sdS[None, :]) @ H.t() + torch.randn(T_eff, M, device=DEV, dtype=DT, generator=g) * rr[None, :]
        Bc = B - B.mean(0, keepdim=True); Sc[i] = (Bc.t() @ Bc) / T_eff
    return Sc

cell = cellsA2["best"]; lam_use = cell["lam_max"]
T_val = int(min(round(25.0 / (0.005 ** 2 * cell["xnorm2"] * lam_use)), 4000))
eps_val = 5.0 / (np.sqrt(cell["xnorm2"]) * np.sqrt(T_val * lam_use))
evals, evecs = torch.linalg.eigh(cell["JB"]); c_top = evecs[:, -1].cpu().numpy()
c_np = c_top / np.linalg.norm(c_top) * eps_val * np.sqrt(cell["xnorm2"])
delta_px = cell["Phi_b"].cpu().numpy() @ c_np; W = cell["make_W"](c_np); V0W = float((cell["V0"] * W).sum().item())
Sc0 = gen_S(cell, 150, T_val, None, 41); Sc1 = gen_S(cell, 150, T_val, delta_px, 42)
t0v = torch.einsum("aij,ij->a", Sc0, W).cpu().numpy() - V0W; t1v = torch.einsum("aij,ij->a", Sc1, W).cpu().numpy() - V0W
d_emp = float((t1v.mean() - t0v.mean()) / (t0v.std() + 1e-12)); A = float(np.mean(t1v[:, None] > t0v[None, :]))
RES["attack2"]["mc_validation"] = dict(cell="best_matched_target", eps=round(eps_val, 5), T_eff=T_val,
    dprime_analytic=round(np.sqrt(T_val * float(c_np @ (cell["JB"].cpu().numpy() @ c_np))), 2),
    dprime_empirical=round(d_emp, 2), auc=round(A, 4))
print(f"  MC-validate (matched target, eps={eps_val*100:.2f}%, T={T_val}): d'_emp={d_emp:.2f} AUC={A:.3f}", flush=True)
save()
plt.figure(figsize=(6.5, 4.5))
for nm, ls in [("best", "-"), ("median", "--")]:
    fr = RES["attack2"]["frontier"][nm]
    plt.plot(Tgrid, [e * 100 for e in fr["energy_spread"]], ls, color="#1f77b4", label=f"{nm} energy-spread")
    plt.plot(Tgrid, [e * 100 for e in fr["matched_target"]], ls, color="#d62728", label=f"{nm} matched-target")
plt.axhline(0.5, color="green", ls=":", label="0.5% floor"); plt.axhline(1.0, color="0.6", ls=":")
plt.xscale("log", base=2); plt.yscale("log"); plt.xlabel("T (banks)"); plt.ylabel("eps_min(T) (%, d'=5)")
plt.title("Detection frontier (corrected shot, optimal codes)"); plt.legend(fontsize=7); plt.grid(alpha=0.3)
plt.tight_layout(); plt.savefig("BOLD_PUSH_frontier.png", dpi=130); print("saved BOLD_PUSH_frontier.png", flush=True)
print(f"[attack2 done {time.time()-t0:.0f}s]", flush=True)

# ============================ ATTACK 3: CUSUM streaming sentinel ============================
print("\n=== ATTACK 3: CUSUM streaming sentinel -- latency vs eps ===", flush=True)
def siegmund_solve_h(delta1, ARL0):
    def ARL0_of_b(b): return (np.exp(delta1 * b) - delta1 * b - 1.0) / (delta1 ** 2 / 2.0)
    lo, hi = 0.0, 5.0 + 3.0 * np.log(ARL0) / max(delta1, 1e-6)
    for _ in range(90):
        mid = 0.5 * (lo + hi)
        if ARL0_of_b(mid) < ARL0: lo = mid
        else: hi = mid
    return 0.5 * (lo + hi)
def siegmund_ARL1(delta1, b): return (np.exp(-delta1 * b) + delta1 * b - 1.0) / (delta1 ** 2 / 2.0)

cell = cellsA2["best"]; xn2 = cell["xnorm2"]
RES["attack3"] = {"cell": "best", "banks_per_hour": BANKS_PER_HOUR, "ms_per_bank": MS_PER_BANK, "latency": {}}
for eps in [0.005, 0.01, 0.02, 0.05]:
    delta1 = eps * np.sqrt(xn2 * cell["lam_mean"]); row = {}
    for arl_name, ARL0 in [("1/hour", BANKS_PER_HOUR), ("1/day", BANKS_PER_DAY)]:
        b = siegmund_solve_h(delta1, ARL0); L = siegmund_ARL1(delta1, b)
        row[arl_name] = dict(latency_banks=round(float(L), 0), latency_sec=round(float(L) * MS_PER_BANK / 1000, 1),
                             threshold_b=round(float(b), 2), delta1_per_bank=round(float(delta1), 4))
    RES["attack3"]["latency"][f"{eps*100:.1f}%"] = row
    print(f"  eps={eps*100:.1f}%: @1/hour={row['1/hour']['latency_banks']:.0f} banks ({row['1/hour']['latency_sec']:.0f}s)"
          f" | @1/day={row['1/day']['latency_banks']:.0f} banks ({row['1/day']['latency_sec']:.0f}s)", flush=True)
save()

def per_bank_scores(cell, n_banks, delta_np, W, V0W, seed):
    Um, S, x, R = cell["Um"], cell["S"], cell["x"], cell["R"]
    xt = x if delta_np is None else x + torch.tensor(delta_np, device=DEV, dtype=DT)
    A2 = cell["Pt"] * xt[None, :]; H = (Um.t() @ A2.t()).t(); sdS = torch.sqrt(S); mfull = cell["Pt"] @ xt
    rr = torch.sqrt(R); g = torch.Generator(device=DEV); g.manual_seed(int(seed))
    Z = torch.randn(n_banks, Um.shape[1], device=DEV, dtype=DT, generator=g)
    B = mfull[None, :] + (Z * sdS[None, :]) @ H.t() + torch.randn(n_banks, M, device=DEV, dtype=DT, generator=g) * rr[None, :]
    return (torch.einsum("ti,tj,ij->t", B, B, W) - V0W).cpu().numpy()
eps_mc = 0.02
c_es = np.random.default_rng(7).standard_normal(cell["db"]); c_es = c_es / np.linalg.norm(c_es) * eps_mc * np.sqrt(xn2)
delta_px = cell["Phi_b"].cpu().numpy() @ c_es; Wc = cell["make_W"](c_es); V0Wc = float((cell["V0"] * Wc).sum().item())
s0 = per_bank_scores(cell, 40000, None, Wc, V0Wc, 101); mu0, sd0 = s0.mean(), s0.std()
d1_emp = (per_bank_scores(cell, 40000, delta_px, Wc, V0Wc, 102).mean() - mu0) / sd0
k_ref = d1_emp / 2.0
def cusum_alarms(z, k, h):
    """all run-lengths in stream z (CUSUM with reset at alarm)."""
    rls = []; g = 0.0; last = 0
    for i in range(len(z)):
        g = g + z[i] - k
        if g < 0: g = 0.0
        if g > h: rls.append(i - last + 1); last = i + 1; g = 0.0
    return rls
def cusum_first(z, k, h):
    g = 0.0
    for i in range(len(z)):
        g = max(0.0, g + z[i] - k)
        if g > h: return i + 1
    return len(z)
ARL0_mc = 500.0                         # tractable ARL0 to validate the Siegmund formula
h_use = siegmund_solve_h(d1_emp, ARL0_mc)
z0 = (per_bank_scores(cell, 40000, None, Wc, V0Wc, 200) - mu0) / sd0   # long H0 stream
arl0_emp = float(np.mean(cusum_alarms(z0, k_ref, h_use)))
lat = [cusum_first((per_bank_scores(cell, 3000, delta_px, Wc, V0Wc, 300 + tr) - mu0) / sd0, k_ref, h_use) for tr in range(150)]
L_mc = float(np.median(lat)); L_pred = siegmund_ARL1(d1_emp, h_use)
RES["attack3"]["mc_validation"] = dict(eps=eps_mc, ARL0_target=ARL0_mc, ARL0_empirical=round(arl0_emp, 0),
    delta1_empirical=round(float(d1_emp), 4), threshold_h=round(float(h_use), 2),
    latency_mc_median_banks=round(L_mc, 0), latency_siegmund_banks=round(float(L_pred), 0),
    ratio=round(L_mc / max(L_pred, 1e-9), 2))
print(f"  MC-validate CUSUM (eps=2%): ARL0 target {ARL0_mc:.0f} vs empirical {arl0_emp:.0f}; "
      f"latency MC={L_mc:.0f} vs Siegmund={L_pred:.0f} (ratio {L_mc/max(L_pred,1e-9):.2f})", flush=True)
save()
plt.figure(figsize=(6.5, 4.5)); epss = [0.005, 0.01, 0.02, 0.05]
for arl_name, col in [("1/hour", "#1f77b4"), ("1/day", "#d62728")]:
    lat_b = [RES["attack3"]["latency"][f"{e*100:.1f}%"][arl_name]["latency_banks"] for e in epss]
    plt.plot([e * 100 for e in epss], lat_b, "o-", color=col, label=f"ARL0 = {arl_name}")
plt.xscale("log"); plt.yscale("log"); plt.xlabel("anomaly energy eps (%)"); plt.ylabel("detection latency (banks)")
plt.gca().secondary_yaxis("right", functions=(lambda b: b * MS_PER_BANK / 1000, lambda s: s * 1000 / MS_PER_BANK)).set_ylabel("latency (s)")
plt.title("CUSUM streaming sentinel latency (best cell)"); plt.legend(); plt.grid(alpha=0.3)
plt.tight_layout(); plt.savefig("BOLD_PUSH_latency.png", dpi=130); print("saved BOLD_PUSH_latency.png", flush=True)
print(f"[all attacks done {time.time()-t0:.0f}s]", flush=True)
