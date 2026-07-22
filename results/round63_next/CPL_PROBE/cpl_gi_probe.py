#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
CPL-GI micro-probe  --  ROUND63 R29 section 6 reserve lane.

Activated by the recorded DEV_GATE_FAIL of DOPS-GI (Probe B). This is the LAST
software-method candidate. Kill bar (R29 section 6, frozen BEFORE any result):
at the PRIMARY cell (sigma=40%, t_c=16 frames, the declared-coherence cell),
CPL-GI (A3) must achieve

    (a) median PSNR gain >= 1.0 dB over the BEST of {A1 naive, A2 ratio-norm}
        per scene (five-seed means), AND
    (b) >= 5/6 DEV scenes positive,
and no-harm at sigma=0 (CPL vs best comparator median >= -0.25 dB).
Anything less = CPL_KILL.

THE METHOD (CPL-GI). Complementary single-pixel acquisition. Hadamard rows
h in {-1,+1}^1024 (full 1024-row, sequency order) realized as physical pattern
pairs m+ = (1+h)/2, m- = (1-h)/2 (constant-sum: m+ + m- = 1). Each Hadamard
coefficient = one packet of q=2 adjacent complementary exposures projected
back-to-back. Bucket counts under multiplicative gain drift:
    Y+_b ~ Poisson(a_{2b}  * phi * m+_b^T x)
    Y-_b ~ Poisson(a_{2b+1}* phi * m-_b^T x)
Because m+ + m- = 1, conditioning on the packet total N_b = Y+_b + Y-_b
removes the packet gain EXACTLY when gain is constant within the packet:
    Y+_b | N_b ~ Binomial(N_b, m+_b^T x / 1^T x).
CPL reconstruction = exact conditional (binomial) MLE over x>=0 + TV penalty.
With the image flux separately anchored (1^T x = Shat), the conditional NLL
    -sum_b [ Y+_b log(m+_b^T x) + Y-_b log(m-_b^T x) ]      (DC packet inert)
is convex on x>=0 (- log of a positive linear form). Flux is pinned by a firm
quadratic anchor to Shat = mean_b(N_b)/phi (the stationarity flux anchor,
E[a]=1); the same Shat is the global scale used by the ratio arm A2, so any
flux bias affects A2 and A3 identically.

ARMS (identical simulated data per (scene,seed,cell); processing-only diffs):
  A1 naive differential  : c_b = Y+_b - Y-_b ; standard Hadamard inversion,
                           then TV denoise (drift-vulnerable baseline).
  A2 ratio normalization : r_b = (Y+_b - Y-_b)/(Y+_b + Y-_b), scaled by global
                           mean packet total (Sensors-2023-style); same inversion
                           + TV denoise as A1.
  A3 CPL-GI (the method) : exact conditional binomial MLE + TV.
  A4 oracle (descriptive): divide each exposure by TRUE a_t, then A1. Headroom.

Read-only on all inputs; writes only results/round63_next/CPL_PROBE/. CPU.
"""
import os, json, time, platform, hashlib
from datetime import datetime, timezone
import numpy as np
from scipy.optimize import minimize

T0 = time.time()
REPO = "D:/GI_another"
OUT = os.path.join(REPO, "results", "round63_next", "CPL_PROBE")
os.makedirs(OUT, exist_ok=True)
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

# ------------------------------------------------------------------ frozen scenes
SCENES = ["bridge_contour_1", "bridge_microtex_1", "bridge_twopop_1",
          "bridge_control_2", "bridge_contour_3", "bridge_microtex_3"]
scene_x = {}; scene_sha = {}
for sc in SCENES:
    p = os.path.join(REPO, "data", "r63_bridge_scenes", sc + ".npz")
    scene_x[sc] = np.load(p)["x"].astype(np.float64).ravel()   # (1024,), in [0,1]
    scene_sha[sc] = sha(p)
NPIX = 1024
SIDE = 32

# ------------------------------------------------------------------ sequency Hadamard
def sylvester(n):
    H = np.array([[1.0]])
    while H.shape[0] < n:
        H = np.block([[H, H], [H, -H]])
    return H
def sequency_order(H):
    # number of sign changes per row -> ascending = Walsh/sequency order
    changes = (H[:, 1:] != H[:, :-1]).sum(axis=1)
    return np.argsort(changes, kind="stable")
Hsyl = sylvester(NPIX)
order = sequency_order(Hsyl)
H = Hsyl[order].copy()                 # (1024,1024) +/-1, sequency order, row0 = DC (all +1)
assert np.allclose(H @ H.T, NPIX * np.eye(NPIX)), "Hadamard rows not orthogonal"
assert np.all(H[0] == 1.0), "row0 not DC"
Mp = (1.0 + H) / 2.0                   # m+ patterns as rows (1024,1024), in {0,1}
Mm = (1.0 - H) / 2.0                   # m- patterns
MpT = Mp.T.copy(); MmT = Mm.T.copy()
HT = H.T.copy()
# sign-change count of DC row = 0; all non-DC Walsh rows are balanced (512/512)
row_ones = Mp.sum(axis=1)
assert row_ones[0] == NPIX and np.all(row_ones[1:] == NPIX // 2), "unexpected pattern balance"

# ------------------------------------------------------------------ photon scale phi
# Target: mean detected counts per exposure ~= Probe B bucket scale (~2200/frame,
# DOSE=4.0 in probe_b.py). For complementary patterns mean(m+/-^T x)=S/2 per scene.
TARGET_COUNTS = 2200.0
mean_member = np.mean([scene_x[sc].sum() / 2.0 for sc in SCENES])   # mean over scenes of S/2
PHI = TARGET_COUNTS / mean_member
log(f"phi (photon scale) = {PHI:.4f}  (target {TARGET_COUNTS:.0f} counts/exposure; mean S/2={mean_member:.2f})")

# ------------------------------------------------------------------ frozen OU gain (lifted from probe_b.py ou_path)
def ou_path(rng, sigma_l, tc, N):
    phi = np.exp(-1.0 / tc); mu = -0.5 * sigma_l ** 2
    eps = rng.standard_normal(N); l = np.empty(N)
    l[0] = mu + sigma_l * eps[0]; sd = sigma_l * np.sqrt(max(1 - phi ** 2, 0.0))
    for n in range(1, N):
        l[n] = mu + phi * (l[n - 1] - mu) + sd * eps[n]
    return np.exp(l)

N_EXP = 2 * NPIX          # 2048 complementary exposures
SEEDS = list(range(5))    # 5 paired gain/noise seeds shared across arms

# frozen cells: primary (gated) + secondary (descriptive, no gate)
PRIMARY = dict(name="primary_sig40_tc16", sigma=0.40, tc=16.0, gate=True)
CELLS = [
    PRIMARY,
    dict(name="sec_sig40_tc2",  sigma=0.40, tc=2.0,   gate=False),   # assumption-stressed
    dict(name="sec_sig40_tc64", sigma=0.40, tc=64.0,  gate=False),   # easy
    dict(name="sec_sig15_tc16", sigma=0.15, tc=16.0,  gate=False),   # milder drift
    dict(name="control_sig0",   sigma=0.0,  tc=16.0,  gate=False),   # gain-free no-harm
]
CONTROL = "control_sig0"

# ------------------------------------------------------------------ smoothed isotropic TV (value + gradient)
TV_EPS = 1e-3
def tv_val_grad(x):
    X = x.reshape(SIDE, SIDE)
    dh = np.zeros_like(X); dv = np.zeros_like(X)
    dh[:, :-1] = X[:, 1:] - X[:, :-1]
    dv[:-1, :] = X[1:, :] - X[:-1, :]
    mag = np.sqrt(dh * dh + dv * dv + TV_EPS * TV_EPS)
    val = float(mag.sum())
    gh = dh / mag; gv = dv / mag
    G = np.zeros_like(X)
    G[:, :-1] += -gh[:, :-1]; G[:, 1:] += gh[:, :-1]
    G[:-1, :] += -gv[:-1, :]; G[1:, :] += gv[:-1, :]
    return val, G.ravel()

# ------------------------------------------------------------------ arms
def hadamard_inverse(coeff):
    # coeff_b ~ (H x)_b ; H H^T = 1024 I -> x = H^T coeff / 1024
    return HT @ coeff / NPIX

def rof_denoise(x0, lam):
    # min_{x>=0} 0.5||x-x0||^2 + lam*TV(x)
    x0c = np.clip(x0, 0.0, None)
    if lam <= 0.0:
        return x0c
    def f(x):
        tv, tvg = tv_val_grad(x)
        r = x - x0
        return 0.5 * float(r @ r) + lam * tv, r + lam * tvg
    res = minimize(f, x0c, jac=True, method="L-BFGS-B",
                   bounds=[(0.0, None)] * NPIX,
                   options=dict(maxiter=500, ftol=1e-12, gtol=1e-8))
    return np.clip(res.x, 0.0, None)

def arm_A1(Yp, Ym, lam):
    c = Yp - Ym
    x0 = hadamard_inverse(c) / PHI            # radiometric (assumes E[a]=1)
    return rof_denoise(x0, lam)

def arm_A2(Yp, Ym, lam):
    N = Yp + Ym
    r = np.where(N > 0, (Yp - Ym) / np.maximum(N, 1e-12), 0.0)
    scale = N.mean() / PHI                    # global mean packet total (radiometric) = Shat
    x0 = hadamard_inverse(r * scale)
    return rof_denoise(x0, lam)

def arm_A4(Yp, Ym, a, lam):
    # oracle known-gain correction: divide each exposure by TRUE gain, then A1
    ap = a[0::2]; am = a[1::2]
    c = Yp / ap - Ym / am
    x0 = hadamard_inverse(c) / PHI
    return rof_denoise(x0, lam)

# A3 flux-anchor strength (numerical; frozen after sanity verification of flux pinning)
RHO_MULT = 200.0
EPS_LOG = 1e-6
def arm_A3(Yp, Ym, lam3, rho=None):
    # Scale-invariant conditional NLL: -sum[Y+ log(m+^T x) + Y- log(m-^T x)] + Ntot*log(1^T x).
    # This is homogeneous degree 0 (invariant to x -> c x), so its radial gradient is 0 and
    # the flux anchor alone fixes the scale to Shat. The DC packet (m-_0=0, Y-_0=0 always) is
    # automatically inert. On the flux slice {1^T x = Shat} the NLL is convex (- log of a
    # positive linear form); the firm anchor confines the solution to that slice.
    N = Yp + Ym; Ntot = float(N.sum())
    Shat = N.mean() / PHI
    if rho is None:
        xu = (Shat / NPIX) * np.ones(NPIX)
        up = np.maximum(Mp @ xu, EPS_LOG); um = np.maximum(Mm @ xu, EPS_LOG)
        Hdiag = MpT @ (Yp / (up * up)) + MmT @ (Ym / (um * um))
        rho = RHO_MULT * float(Hdiag.mean()) / NPIX
    ones = np.ones(NPIX)
    def f(x):
        S = float(x.sum())
        up = np.maximum(Mp @ x, EPS_LOG); um = np.maximum(Mm @ x, EPS_LOG)
        nll = -float(Yp @ np.log(up)) - float(Ym @ np.log(um)) + Ntot * np.log(max(S, EPS_LOG))
        gnll = -(MpT @ (Yp / up)) - (MmT @ (Ym / um)) + (Ntot / max(S, EPS_LOG)) * ones
        d = S - Shat
        anch = 0.5 * rho * d * d
        ganch = rho * d * ones
        tv, tvg = tv_val_grad(x)
        return nll + anch + lam3 * tv, gnll + ganch + lam3 * tvg
    x_init = np.clip(arm_A2(Yp, Ym, 0.0), EPS_LOG, None)
    s0 = x_init.sum()
    if s0 > 0:
        x_init = x_init * (Shat / s0)
    res = minimize(f, x_init, jac=True, method="L-BFGS-B",
                   bounds=[(0.0, None)] * NPIX,
                   options=dict(maxiter=1000, ftol=1e-12, gtol=1e-8))
    x = np.clip(res.x, 0.0, None)
    return x, Shat, abs(x.sum() / Shat - 1.0)

def psnr(xh, x):
    mse = float(np.mean((xh - x) ** 2))
    return 10.0 * np.log10(1.0 / mse) if mse > 0 else 99.0

# ------------------------------------------------------------------ data simulation (shared across arms)
def simulate(sc, seed, sigma, tc):
    x = scene_x[sc]
    gain_rng = np.random.default_rng(1000 + seed)
    if sigma <= 0.0:
        a = np.ones(N_EXP)
    else:
        a = ou_path(gain_rng, sigma, tc, N_EXP)
    mp = Mp @ x; mm = Mm @ x
    lam_p = a[0::2] * PHI * mp
    lam_m = a[1::2] * PHI * mm
    noise_rng = np.random.default_rng(5000 + seed)
    Yp = noise_rng.poisson(lam_p).astype(np.float64)
    Ym = noise_rng.poisson(lam_m).astype(np.float64)
    return Yp, Ym, a

# ================================================================== SANITY
def sanity_checks():
    log("SANITY 1: TV gradient finite-difference check")
    rng = np.random.default_rng(0)
    xt = np.clip(rng.standard_normal(NPIX) * 0.3 + 0.5, 0, None)
    v0, g0 = tv_val_grad(xt)
    h = 1e-6; err = 0.0
    for k in rng.choice(NPIX, 20, replace=False):
        xp = xt.copy(); xp[k] += h; xm = xt.copy(); xm[k] -= h
        num = (tv_val_grad(xp)[0] - tv_val_grad(xm)[0]) / (2 * h)
        err = max(err, abs(num - g0[k]))
    log(f"  max |num-analytic| TV grad = {err:.2e}")
    assert err < 1e-4, "TV gradient mismatch"

    log("SANITY 2: no-drift no-noise, lam=0 -> all arms recover the scene")
    facts = {}
    for sc in SCENES:
        x = scene_x[sc]
        a = np.ones(N_EXP)
        mp = Mp @ x; mm = Mm @ x
        Yp = a[0::2] * PHI * mp     # exact expected counts (noiseless)
        Ym = a[1::2] * PHI * mm
        p1 = psnr(arm_A1(Yp, Ym, 0.0), x)
        p2 = psnr(arm_A2(Yp, Ym, 0.0), x)
        x3, Sh, fe = arm_A3(Yp, Ym, 0.0)
        p3 = psnr(x3, x)
        p4 = psnr(arm_A4(Yp, Ym, a, 0.0), x)
        facts[sc] = dict(A1=p1, A2=p2, A3=p3, A4=p4, A3_flux_err=fe)
        log(f"  {sc:18s} A1={p1:6.1f} A2={p2:6.1f} A3={p3:6.1f} A4={p4:6.1f} dB  (A3 flux err={fe:.2e})")
        assert min(p1, p2, p3, p4) > 55.0, f"noiseless recovery failed for {sc}"
    return facts

# ================================================================== LAMBDA SELECTION (control cell, sigma=0)
GRID = [0.0, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1]   # dimensionless TV weight (image domain)
# For A3, NLL data-term curvature differs from the unit-curvature denoise term; map
# the same image-domain grid through the analytic curvature scale H_SCALE (one frozen
# scalar) so the SAME grid means the same smoothness/data tradeoff for every arm.
def a3_curvature_scale():
    # mean diagonal NLL Hessian at the uniform anchored image, control-cell counts,
    # averaged over the 6 control scenes x 5 seeds -> one frozen scalar.
    vals = []
    for sc in SCENES:
        for seed in SEEDS:
            Yp, Ym, a = simulate(sc, seed, 0.0, 16.0)
            N = Yp + Ym; Shat = N.mean() / PHI
            xu = (Shat / NPIX) * np.ones(NPIX)
            up = np.maximum(Mp @ xu, EPS_LOG); um = np.maximum(Mm @ xu, EPS_LOG)
            Hdiag = MpT @ (Yp / (up * up)) + MmT @ (Ym / (um * um))
            vals.append(float(Hdiag.mean()))
    return float(np.mean(vals))

def select_lambdas():
    H_SCALE = a3_curvature_scale()
    log(f"A3 curvature scale H_SCALE = {H_SCALE:.4f} (maps image-domain grid -> NLL TV weight)")
    # evaluate each arm at each grid lambda on the control cell (sigma=0)
    perf = {a: {g: [] for g in GRID} for a in ("A1", "A2", "A3")}
    for sc in SCENES:
        x = scene_x[sc]
        for seed in SEEDS:
            Yp, Ym, a = simulate(sc, seed, 0.0, 16.0)
            for g in GRID:
                perf["A1"][g].append(psnr(arm_A1(Yp, Ym, g), x))
                perf["A2"][g].append(psnr(arm_A2(Yp, Ym, g), x))
                x3, _, _ = arm_A3(Yp, Ym, g * H_SCALE)
                perf["A3"][g].append(psnr(x3, x))
    sel = {}; curves = {}
    for a in ("A1", "A2", "A3"):
        means = {g: float(np.mean(perf[a][g])) for g in GRID}
        best_g = max(means, key=means.get)
        sel[a] = best_g
        curves[a] = means
        log(f"  {a} control-PSNR by lambda: " + "  ".join(f"{g:.0e}:{means[g]:.2f}" for g in GRID)
            + f"   -> selected lambda_img={best_g:.0e}")
    # A4 uses A1's selected lambda (same denoise convention)
    sel["A4"] = sel["A1"]
    return sel, curves, H_SCALE

# ================================================================== GRID RUN
def run_cell(cell, lam_sel, H_SCALE, lam0=False):
    sigma = cell["sigma"]; tc = cell["tc"]
    lamA1 = 0.0 if lam0 else lam_sel["A1"]
    lamA2 = 0.0 if lam0 else lam_sel["A2"]
    lamA4 = 0.0 if lam0 else lam_sel["A4"]
    lam3 = 0.0 if lam0 else lam_sel["A3"] * H_SCALE
    per_scene = {}
    flux_errs = []
    for sc in SCENES:
        x = scene_x[sc]
        pv = {a: [] for a in ("A1", "A2", "A3", "A4")}
        for seed in SEEDS:
            Yp, Ym, a = simulate(sc, seed, sigma, tc)
            pv["A1"].append(psnr(arm_A1(Yp, Ym, lamA1), x))
            pv["A2"].append(psnr(arm_A2(Yp, Ym, lamA2), x))
            x3, Sh, fe = arm_A3(Yp, Ym, lam3); flux_errs.append(fe)
            pv["A3"].append(psnr(x3, x))
            pv["A4"].append(psnr(arm_A4(Yp, Ym, a, lamA4), x))
        means = {a: float(np.mean(pv[a])) for a in pv}   # five-seed means
        base = max(means["A1"], means["A2"])             # best of {A1,A2}
        per_scene[sc] = dict(
            A1=means["A1"], A2=means["A2"], A3=means["A3"], A4=means["A4"],
            best_comparator=base,
            dQ_cpl_vs_bestcomp=means["A3"] - base,
            dQ_cpl_vs_A2=means["A3"] - means["A2"],
            seed_vals={a: pv[a] for a in pv},
        )
    dQ = [per_scene[sc]["dQ_cpl_vs_bestcomp"] for sc in SCENES]
    summary = dict(
        sigma=sigma, tc=tc,
        per_scene_dQ={sc: per_scene[sc]["dQ_cpl_vs_bestcomp"] for sc in SCENES},
        median_dQ=float(np.median(dQ)),
        n_positive=int(sum(1 for d in dQ if d > 0)),
        max_A3_flux_err=float(np.max(flux_errs)),
    )
    return per_scene, summary

# ================================================================== MAIN
sanity = sanity_checks()
lam_sel, lam_curves, H_SCALE = select_lambdas()
log(f"selected lambdas (image domain): {lam_sel}")

results = {}
log("=== primary + secondary cells (frozen selected lambda) ===")
for cell in CELLS:
    ps, summ = run_cell(cell, lam_sel, H_SCALE, lam0=False)
    results[cell["name"]] = dict(cell=cell, per_scene=ps, summary=summ)
    log(f"  {cell['name']:18s} median dQ(CPL vs best[A1,A2])={summ['median_dQ']:+.3f} dB  "
        f"n_pos={summ['n_positive']}/6  (A3 flux err<={summ['max_A3_flux_err']:.1e})")

# robustness: lambda=0 (no TV) descriptive view for every cell
log("=== robustness: lambda=0 (no TV) ===")
results_lam0 = {}
for cell in CELLS:
    ps, summ = run_cell(cell, lam_sel, H_SCALE, lam0=True)
    results_lam0[cell["name"]] = dict(cell=cell, per_scene=ps, summary=summ)
    log(f"  {cell['name']:18s} median dQ={summ['median_dQ']:+.3f} dB  n_pos={summ['n_positive']}/6")

# ------------------------------------------------------------------ VERDICT (primary cell)
prim = results[PRIMARY["name"]]["summary"]
ctrl = results[CONTROL]["summary"]
# no-harm: CPL vs best comparator on gain-free control (median dQ >= -0.25 dB)
ctrl_ps = results[CONTROL]["per_scene"]
noharm_dQ = [ctrl_ps[sc]["dQ_cpl_vs_bestcomp"] for sc in SCENES]
noharm_median = float(np.median(noharm_dQ))

cond_a = prim["median_dQ"] >= 1.0
cond_b = prim["n_positive"] >= 5
cond_c = noharm_median >= -0.25
verdict = "CPL_PASS" if (cond_a and cond_b and cond_c) else "CPL_KILL"

gate = dict(
    generated_utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
    method="CPL-GI exact conditional binomial MLE + TV",
    primary_cell=PRIMARY["name"],
    kill_bar="median dQ(CPL vs best[A1,A2]) >= 1.0 dB AND >=5/6 scenes positive; no-harm sigma=0 median >= -0.25 dB",
    per_scene_dQ=prim["per_scene_dQ"],
    median=prim["median_dQ"],
    n_positive=prim["n_positive"],
    no_harm=dict(control_median_dQ=noharm_median,
                 per_scene={sc: ctrl_ps[sc]["dQ_cpl_vs_bestcomp"] for sc in SCENES}),
    conditions=dict(
        a=dict(name="median dQ >= 1.0 dB", value=prim["median_dQ"], threshold=1.0, passed=bool(cond_a)),
        b=dict(name="#scenes dQ>0 >= 5 of 6", value=prim["n_positive"], threshold=5, passed=bool(cond_b)),
        c=dict(name="no-harm sigma=0 median dQ >= -0.25 dB", value=noharm_median, threshold=-0.25, passed=bool(cond_c)),
    ),
    verdict=verdict,
)
with open(os.path.join(OUT, "CPL_GATE_VERDICT.json"), "w", encoding="utf-8") as f:
    json.dump(gate, f, indent=2)

# ------------------------------------------------------------------ full results json
meta = dict(
    utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
    python=platform.python_version(), numpy=np.__version__, platform=platform.platform(),
    runtime_s=time.time() - T0,
    scenes={sc: scene_sha[sc] for sc in SCENES},
    config=dict(
        Npix=NPIX, side=SIDE, N_exposures=N_EXP, patterns="full 1024-row sequency Hadamard, complementary pairs",
        phi=PHI, target_counts_per_exposure=TARGET_COUNTS, seeds=SEEDS,
        cells=[c["name"] for c in CELLS], primary_cell=PRIMARY["name"], control_cell=CONTROL,
        ou_mu_rule="-0.5*sigma^2 (E[a]=1)", ou_source="lifted from probe_b.py ou_path (exact AR(1))",
        tv="smoothed isotropic, eps=%.0e" % TV_EPS, tv_eps=TV_EPS,
        lambda_grid=GRID, lambda_selected_image_domain=lam_sel, lambda_control_curves=lam_curves,
        a3_curvature_scale=H_SCALE, a3_flux_anchor_rho_mult=RHO_MULT,
        lambda_rule="per-arm lambda maximizing mean control-cell (sigma=0) PSNR over the shared grid; "
                    "frozen before any drift cell; A3 grid mapped through analytic curvature scale H_SCALE.",
    ),
    sanity_noiseless_lam0=sanity,
)
RES = dict(meta=meta, results=results, results_lambda0=results_lam0, gate=gate)
with open(os.path.join(OUT, "cpl_probe_results.json"), "w", encoding="utf-8") as f:
    json.dump(RES, f, indent=2)
log(f"VERDICT = {verdict}  (median dQ={prim['median_dQ']:+.3f} dB, n_pos={prim['n_positive']}/6, "
    f"no-harm median={noharm_median:+.3f} dB)")
log(f"TOTAL runtime {time.time()-T0:.1f}s")
