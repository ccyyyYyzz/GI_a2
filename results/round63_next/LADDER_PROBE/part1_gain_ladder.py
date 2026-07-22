#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
PART 1 of the RECORD-REFINEMENT LADDER probe -- CHAPTER 2 (gain drift), PRIMARY.
ROUND63-NEXT / GI_a2.  Image-level.

QUESTION.  At the CPL primary cell the known-gain oracle sits 10-15 dB above ALL
bucket-count software arms (results/round63_next/CPL_PROBE/).  The program
identity leaves ONE free axis at fixed photon budget: the RECORD's temporal
richness.  Same photons, same patterns, same dwell -- but the readout keeps finer
temporal structure (sub-bins of each exposure, down to the timestamp limit)
instead of collapsing to one bucket total per exposure.  How much of the oracle
gap comes back as the record refines, AT IMAGE LEVEL, at identical photon cost?

SETUP (lifted VERBATIM from results/round63_next/CPL_PROBE/cpl_gi_probe.py):
  32x32, NPIX=1024; 6 DEV scenes; full 1024-row sequency Hadamard realised as
  2048 complementary exposures (m+, m- pairs); PHI matched to 2200 counts/exposure;
  5 paired seeds; frozen OU gain (sigma_l, t_c) and frozen isotropic-TV rule.
  Cells: primary (sigma=40%, t_c=16) + secondaries (t_c=2, t_c=64, sigma=15%,
  sigma=0 control) -- the exact CPL CELLS list.

THE LADDER.  Per exposure, subdivide the dwell into B sub-bins, B in {1,4,16,64}.
Counts at coarser rungs are literally SUMS of the finest (B_MAX=64) sub-bin
counts -> IDENTICAL photon budget across the ladder (a coarser record is a
binning of the finer one).  B=1 = the bucket total (the frozen CPL record).
B=64 = timestamp/continuous-observation proxy.

MODELLING NOTE (stated as a finding).  Two gain models:
  (M0) frozen WITHIN-EXPOSURE-CONSTANT gain (the CPL model): a is constant across
       an exposure's sub-bins.  Then the sub-bin counts are Poisson THINNINGS of
       the bucket total; the total is a sufficient statistic and B>1 provably adds
       ZERO information vs B=1.  We verify this numerically (the ladder is flat).
  (M1) CONTINUOUS-DRIFT gain (the honest case at fast drift): the SAME OU drifts
       continuously in time, t_c measured in frames; sub-bins genuinely resolve
       the within-exposure path.  This is the honest PRIMARY variant and is
       labelled as a model variant throughout.

RECONSTRUCTION (SAME estimator at every rung, only the record differs).  Joint
state-space estimation, alternating:
  (gain step)  given x, each rung-bin count is a Poisson observation of
      exp(l_bin) * (PHI/B) * (pattern^T x); log-linearise (delta method) and run
      an RTS smoother on the log-gain with the OU prior at the rung's resolution.
  (image step) given the smoothed per-exposure mean gain a_bar_e, gain-correct the
      bucket counts and invert with the frozen Hadamard + isotropic-TV rule
      (= CPL arm A4 with ESTIMATED gains).
Baselines per cell: (i) B=1 bucket software arms A2 (ratio-norm, rerun in-frame)
and A3 (CPL, rerun in-frame) + the B=1 state-space smoother; (ii) known-gain
oracle A4 (ceiling).  Frozen TV lambda; NO per-rung tuning.

Read-only on inputs; writes only this dir.  CPU.
"""
import os, json, time, platform, hashlib, argparse
from datetime import datetime, timezone
import numpy as np
from scipy.optimize import minimize

T0 = time.time()
REPO = "D:/GI_another"
OUT = os.path.join(REPO, "results", "round63_next", "LADDER_PROBE")
os.makedirs(OUT, exist_ok=True)
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

# ------------------------------------------------------------------ frozen scenes (verbatim)
SCENES = ["bridge_contour_1", "bridge_microtex_1", "bridge_twopop_1",
          "bridge_control_2", "bridge_contour_3", "bridge_microtex_3"]
scene_x = {}; scene_sha = {}
for sc in SCENES:
    p = os.path.join(REPO, "data", "r63_bridge_scenes", sc + ".npz")
    scene_x[sc] = np.load(p)["x"].astype(np.float64).ravel()
    scene_sha[sc] = sha(p)
NPIX = 1024; SIDE = 32

# ------------------------------------------------------------------ sequency Hadamard (verbatim)
def sylvester(n):
    H = np.array([[1.0]])
    while H.shape[0] < n:
        H = np.block([[H, H], [H, -H]])
    return H
def sequency_order(H):
    changes = (H[:, 1:] != H[:, :-1]).sum(axis=1)
    return np.argsort(changes, kind="stable")
Hsyl = sylvester(NPIX); order = sequency_order(Hsyl)
H = Hsyl[order].copy()
assert np.allclose(H @ H.T, NPIX * np.eye(NPIX))
assert np.all(H[0] == 1.0)
Mp = (1.0 + H) / 2.0; Mm = (1.0 - H) / 2.0
MpT = Mp.T.copy(); MmT = Mm.T.copy(); HT = H.T.copy()

# ------------------------------------------------------------------ photon scale phi (verbatim)
TARGET_COUNTS = 2200.0
mean_member = np.mean([scene_x[sc].sum() / 2.0 for sc in SCENES])
PHI = TARGET_COUNTS / mean_member

# ------------------------------------------------------------------ frozen OU (verbatim ou_path) + continuous variant
def ou_path(rng, sigma_l, tc, N):
    """Exposure-resolution OU (the frozen within-constant path)."""
    phi = np.exp(-1.0 / tc); mu = -0.5 * sigma_l ** 2
    eps = rng.standard_normal(N); l = np.empty(N)
    l[0] = mu + sigma_l * eps[0]; sd = sigma_l * np.sqrt(max(1 - phi ** 2, 0.0))
    for n in range(1, N):
        l[n] = mu + phi * (l[n - 1] - mu) + sd * eps[n]
    return np.exp(l)

def ou_path_fine(rng, sigma_l, tc, N_exp, B):
    """Continuous-drift OU discretised at B sub-bins per exposure.  Same stationary
    std sigma_l and correlation time tc (in EXPOSURES); per sub-bin step = 1/B
    exposure.  Returns log-gain reshaped (N_exp, B)."""
    M = N_exp * B
    phi = np.exp(-1.0 / (tc * B)); mu = -0.5 * sigma_l ** 2
    sd = sigma_l * np.sqrt(max(1 - phi ** 2, 0.0))
    eps = rng.standard_normal(M); l = np.empty(M)
    l[0] = mu + sigma_l * eps[0]
    for n in range(1, M):
        l[n] = mu + phi * (l[n - 1] - mu) + sd * eps[n]
    return l.reshape(N_exp, B)

N_EXP = 2 * NPIX; B_MAX = 64
RUNGS = [1, 4, 16, 64]                 # 64 = timestamp/continuous proxy
SEEDS = list(range(5))

PRIMARY = dict(name="primary_sig40_tc16", sigma=0.40, tc=16.0)
CELLS = [
    PRIMARY,
    dict(name="sec_sig40_tc2",  sigma=0.40, tc=2.0),
    dict(name="sec_sig40_tc64", sigma=0.40, tc=64.0),
    dict(name="sec_sig15_tc16", sigma=0.15, tc=16.0),
    dict(name="control_sig0",   sigma=0.0,  tc=16.0),
]

# ------------------------------------------------------------------ frozen TV (verbatim) + selected lambda
TV_EPS = 1e-3
def tv_val_grad(x):
    X = x.reshape(SIDE, SIDE)
    dh = np.zeros_like(X); dv = np.zeros_like(X)
    dh[:, :-1] = X[:, 1:] - X[:, :-1]; dv[:-1, :] = X[1:, :] - X[:-1, :]
    mag = np.sqrt(dh * dh + dv * dv + TV_EPS * TV_EPS)
    val = float(mag.sum()); gh = dh / mag; gv = dv / mag
    G = np.zeros_like(X)
    G[:, :-1] += -gh[:, :-1]; G[:, 1:] += gh[:, :-1]
    G[:-1, :] += -gv[:-1, :]; G[1:, :] += gv[:-1, :]
    return val, G.ravel()
def hadamard_inverse(coeff):
    return HT @ coeff / NPIX
def rof_denoise(x0, lam):
    x0c = np.clip(x0, 0.0, None)
    if lam <= 0.0:
        return x0c
    def f(x):
        tv, tvg = tv_val_grad(x); r = x - x0
        return 0.5 * float(r @ r) + lam * tv, r + lam * tvg
    res = minimize(f, x0c, jac=True, method="L-BFGS-B", bounds=[(0.0, None)] * NPIX,
                   options=dict(maxiter=500, ftol=1e-12, gtol=1e-8))
    return np.clip(res.x, 0.0, None)
LAM_TV = 0.1                       # frozen A1/A4 selected lambda (image domain)
LAM_A2 = 0.1; LAM_A3 = 0.03       # frozen selected lambdas

# ------------------------------------------------------------------ CPL software arms (verbatim, bucket record)
def arm_A2(Yp, Ym, lam):
    N = Yp + Ym
    r = np.where(N > 0, (Yp - Ym) / np.maximum(N, 1e-12), 0.0)
    scale = N.mean() / PHI
    return rof_denoise(hadamard_inverse(r * scale), lam)
EPS_LOG = 1e-6; RHO_MULT = 200.0
def arm_A3(Yp, Ym, lam3, rho=None):
    N = Yp + Ym; Ntot = float(N.sum()); Shat = N.mean() / PHI
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
        d = S - Shat; anch = 0.5 * rho * d * d; ganch = rho * d * ones
        tv, tvg = tv_val_grad(x)
        return nll + anch + lam3 * tv, gnll + ganch + lam3 * tvg
    x_init = np.clip(arm_A2(Yp, Ym, 0.0), EPS_LOG, None); s0 = x_init.sum()
    if s0 > 0:
        x_init = x_init * (Shat / s0)
    res = minimize(f, x_init, jac=True, method="L-BFGS-B", bounds=[(0.0, None)] * NPIX,
                   options=dict(maxiter=1000, ftol=1e-12, gtol=1e-8))
    return np.clip(res.x, 0.0, None)
def arm_A4_gaincorr(Yp, Ym, ap, am, lam):
    """Oracle / gain-corrected inversion: divide each exposure by its (mean) gain."""
    c = Yp / np.maximum(ap, 1e-9) - Ym / np.maximum(am, 1e-9)
    return rof_denoise(hadamard_inverse(c) / PHI, lam)

def psnr(xh, x):
    mse = float(np.mean((xh - x) ** 2))
    return 10.0 * np.log10(1.0 / mse) if mse > 0 else 99.0

# ------------------------------------------------------------------ RTS smoother on log-gain (scalar OU)
def rts_smoother(z, R_obs, valid, mu, phi, Q, sigma_l):
    """Kalman forward + RTS backward for scalar state l_t: l_t = mu + phi(l_{t-1}-mu)+w,
    w~N(0,Q); obs z_t = l_t + v, v~N(0,R_t) where valid_t marks a usable observation.
    Stationary prior N(mu, sigma_l^2).  Vectorised loops in numpy (O(M))."""
    M = z.shape[0]
    a_pred = np.empty(M); P_pred = np.empty(M)
    a_filt = np.empty(M); P_filt = np.empty(M)
    a_pred[0] = mu; P_pred[0] = sigma_l ** 2
    for t in range(M):
        if t > 0:
            a_pred[t] = mu + phi * (a_filt[t - 1] - mu)
            P_pred[t] = phi * phi * P_filt[t - 1] + Q
        if valid[t]:
            S = P_pred[t] + R_obs[t]; K = P_pred[t] / S
            a_filt[t] = a_pred[t] + K * (z[t] - a_pred[t])
            P_filt[t] = (1.0 - K) * P_pred[t]
        else:
            a_filt[t] = a_pred[t]; P_filt[t] = P_pred[t]
    a_sm = np.empty(M); P_sm = np.empty(M)
    a_sm[-1] = a_filt[-1]; P_sm[-1] = P_filt[-1]
    for t in range(M - 2, -1, -1):
        C = phi * P_filt[t] / max(P_pred[t + 1], 1e-30)
        a_sm[t] = a_filt[t] + C * (a_sm[t + 1] - mu - phi * (a_filt[t] - mu))
        P_sm[t] = P_filt[t] + C * C * (P_sm[t + 1] - P_pred[t + 1])
    return a_sm, np.maximum(P_sm, 0.0)

# ------------------------------------------------------------------ joint state-space ladder estimator
def ladder_estimate(Y_rung, s_hint_scenexnpix, B, sigma_l, tc, n_outer=4):
    """Y_rung: (N_EXP, B) rung-bin counts.  Returns reconstructed image x.
    SAME structure at every B; only the record resolution changes."""
    mu = -0.5 * sigma_l ** 2
    phi = np.exp(-1.0 / (tc * B)) if sigma_l > 0 else 0.0
    Q = (sigma_l ** 2) * (1 - phi ** 2) if sigma_l > 0 else 1e-12
    Ye_tot = Y_rung.sum(axis=1)                          # bucket totals (identical photons)
    Yp = Ye_tot[0::2]; Ym = Ye_tot[1::2]
    if sigma_l <= 0:                                      # control: no gain, TV-denoise
        return np.clip(arm_A4_gaincorr(Yp, Ym, np.ones(NPIX), np.ones(NPIX), LAM_TV), 0.0, None)
    x = np.clip(arm_A2(Yp, Ym, 0.0), 1e-9, None)         # gain-robust init
    # exposure->pattern map: exposure 2b uses Mp[b], 2b+1 uses Mm[b]
    for it in range(n_outer):
        s_e = np.empty(N_EXP)
        s_e[0::2] = Mp @ x; s_e[1::2] = Mm @ x
        # gain step: per rung-bin Poisson obs of exp(l)*mu_bin
        mu_bin = (PHI / B) * s_e[:, None]                # (N_EXP, 1) -> broadcasts to (N_EXP,B)
        with np.errstate(divide="ignore", invalid="ignore"):
            z = np.log((Y_rung + 0.5) / np.maximum(mu_bin, 1e-12))
        R = 1.0 / (Y_rung + 0.5)
        # DC m- exposures (s_e=0) carry no gain info; broadcast per-exposure validity
        valid = np.broadcast_to(mu_bin > 1e-9, z.shape)
        zf = z.reshape(-1); Rf = R.reshape(-1); vf = np.ascontiguousarray(valid).reshape(-1)
        l_sm, P_sm = rts_smoother(zf, Rf, vf, mu, phi, Q, sigma_l)
        l_sm = l_sm.reshape(N_EXP, B); P_sm = P_sm.reshape(N_EXP, B)
        # MMSE estimate of the per-exposure mean gain a_bar_e = mean_k E[exp(l_k)|data]
        # = mean_k exp(l_sm + P_sm/2)  (uniform Jensen debias, same rule at every rung)
        abar = np.mean(np.exp(l_sm + 0.5 * P_sm), axis=1)
        # image step: gain-correct + frozen Hadamard/TV inversion
        x = np.clip(arm_A4_gaincorr(Yp, Ym, abar[0::2], abar[1::2], LAM_TV), 0.0, None)
    return x

# ------------------------------------------------------------------ simulate finest record, coarsen to rungs
def simulate_finest(sc, seed, sigma, tc, model):
    """model in {'M0_const','M1_cont'}.  Returns (Y_fine (N_EXP,B_MAX), abar_true
    (N_EXP,)) at IDENTICAL expected photon budget."""
    x = scene_x[sc]
    s_e = np.empty(N_EXP); s_e[0::2] = Mp @ x; s_e[1::2] = Mm @ x
    grng = np.random.default_rng(1000 + seed)
    if sigma <= 0:
        a_fine = np.ones((N_EXP, B_MAX))
    elif model == "M0_const":
        a_e = ou_path(grng, sigma, tc, N_EXP)            # frozen exposure-res path
        a_fine = np.repeat(a_e[:, None], B_MAX, axis=1)  # constant within exposure
    else:                                                # M1_cont
        a_fine = np.exp(ou_path_fine(grng, sigma, tc, N_EXP, B_MAX))
    lam_fine = a_fine * (PHI / B_MAX) * s_e[:, None]
    nrng = np.random.default_rng(5000 + seed)
    Y_fine = nrng.poisson(lam_fine).astype(np.float64)
    abar_true = a_fine.mean(axis=1)
    return Y_fine, abar_true

def coarsen(Y_fine, B):
    g = B_MAX // B
    return Y_fine.reshape(N_EXP, B, g).sum(axis=2)

# ------------------------------------------------------------------ one (scene,seed,cell) across the ladder
def run_cell_scene_seed(sc, seed, cell, model):
    sigma = cell["sigma"]; tc = cell["tc"]; x = scene_x[sc]
    Y_fine, abar_true = simulate_finest(sc, seed, sigma, tc, model)
    Yb = Y_fine.sum(axis=1); Yp = Yb[0::2]; Ym = Yb[1::2]     # bucket totals
    out = {}
    # B=1 software baselines (bucket record)
    out["A2"] = psnr(arm_A2(Yp, Ym, LAM_A2), x)
    out["A3"] = psnr(arm_A3(Yp, Ym, LAM_A3), x)
    # known-gain oracle (ceiling): true per-exposure mean gain
    out["oracle"] = psnr(arm_A4_gaincorr(Yp, Ym, abar_true[0::2], abar_true[1::2], LAM_TV), x)
    # ladder rungs (same estimator, coarsened records)
    for B in RUNGS:
        Yr = coarsen(Y_fine, B)
        xr = ladder_estimate(Yr, x, B, sigma, tc)
        out[f"B{B}"] = psnr(xr, x)
    return out

# ================================================================== MAIN
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenes", type=int, default=len(SCENES))
    ap.add_argument("--seeds", type=int, default=len(SEEDS))
    ap.add_argument("--cells", default="all")     # 'all' or comma names
    ap.add_argument("--model", default="M1_cont") # M1_cont primary; M0_const for vacuity check
    ap.add_argument("--tag", default="full")
    args = ap.parse_args()
    scenes = SCENES[:args.scenes]; seeds = SEEDS[:args.seeds]
    cells = CELLS if args.cells == "all" else [c for c in CELLS if c["name"] in args.cells.split(",")]
    log(f"PART 1 ladder  model={args.model}  scenes={len(scenes)} seeds={len(seeds)} "
        f"cells={[c['name'] for c in cells]}  rungs={RUNGS} B_MAX={B_MAX}")
    log(f"PHI={PHI:.4f}  N_EXP={N_EXP}")

    results = {}
    for cell in cells:
        cres = {}
        for sc in scenes:
            per_seed = {k: [] for k in ["A2", "A3", "oracle"] + [f"B{B}" for B in RUNGS]}
            for seed in seeds:
                o = run_cell_scene_seed(sc, seed, cell, args.model)
                for k in per_seed:
                    per_seed[k].append(o[k])
            cres[sc] = {k: float(np.mean(v)) for k, v in per_seed.items()}
            cres[sc]["_seed_vals"] = {k: per_seed[k] for k in per_seed}
        results[cell["name"]] = cres
        # cell summary
        def med(k): return float(np.median([cres[sc][k] for sc in scenes]))
        best_b1_sw = [max(cres[sc]["A2"], cres[sc]["A3"], cres[sc]["B1"]) for sc in scenes]
        d_ts = [cres[sc]["B64"] - max(cres[sc]["A2"], cres[sc]["A3"]) for sc in scenes]
        d_ts_vs_b1sm = [cres[sc]["B64"] - cres[sc]["B1"] for sc in scenes]
        frac_oracle = []
        for sc in scenes:
            base = max(cres[sc]["A2"], cres[sc]["A3"])
            gap = cres[sc]["oracle"] - base
            frac_oracle.append((cres[sc]["B64"] - base) / gap if gap > 1e-6 else 0.0)
        log(f"[{cell['name']:18s} {args.model}] median: A2={med('A2'):.2f} A3={med('A3'):.2f} "
            f"B1={med('B1'):.2f} B4={med('B4'):.2f} B16={med('B16'):.2f} B64={med('B64'):.2f} "
            f"oracle={med('oracle'):.2f}")
        log(f"    median dB(B64 - best_B1_software)={float(np.median(d_ts)):+.3f}  "
            f"refine-only(B64 - B1_smoother)={float(np.median(d_ts_vs_b1sm)):+.3f}  "
            f"oracle-frac={float(np.median(frac_oracle)):.3f}")
    meta = dict(utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
                python=platform.python_version(), numpy=np.__version__,
                model=args.model, PHI=PHI, N_EXP=N_EXP, B_MAX=B_MAX, rungs=RUNGS,
                scenes=scenes, seeds=seeds, cells=[c["name"] for c in cells],
                scene_sha={sc: scene_sha[sc] for sc in scenes},
                lam_tv=LAM_TV, target_counts=TARGET_COUNTS, runtime_s=time.time() - T0)
    fn = os.path.join(OUT, f"ladder_results_{args.model}_{args.tag}.json")
    json.dump(dict(meta=meta, results=results), open(fn, "w"), indent=2)
    log(f"saved {fn}  ({time.time()-T0:.1f}s)")

if __name__ == "__main__":
    main()
