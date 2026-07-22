"""PROBE A -- IMAGE-level range-space efficiency of the deployed RQL+TV
reconstruction vs the Fisher information floor, under the jittered-hold model.

Design: L0 bank 972 main rows + 52 balanced pre-scan rows (the deployed
bridge_harness path), scaled so mean MAIN operating load == target rho.
Forward model: bridge_harness.simulate_counts_jitter (c=0.05).
Reconstruction: RQL + isotropic TV (deployed, truth-free lam_tv) AND RQL with
no TV (lam_tv=0, pure data term) -- both via solvers.run_arm("RQL", ...).

Fisher floor (range component only; null space NOT bounded):
   Per-measurement info about the load rho_i = a_i.x is I_rho(rho_i) =
   nu*J(rho_i;c)/rho_i^2  (J = per-slot log-rate info from fi_mi).
   Information matrix  I_x = A^T diag(I_rho) A   (n_pix x n_pix).
   Eigendecompose I_x = V diag(mu) V^T.  Range = top-r eigenvectors
   (r = numerical rank of A_full); null = remaining.
   floor = sum_{range} 1/mu_k   (CRB trace bound on the range-space MSE of any
   UNBIASED estimator).  Errors are decomposed in the SAME eigenbasis:
   range_mse = sum_{range} (v_k.e)^2 ,  null_mse = ||e||^2 - range_mse.
Efficiency_range = floor / range_mse_actual.
"""
import json, os, sys, time
import numpy as np

ROOT = r"D:\GI_another"
sys.path.insert(0, os.path.join(ROOT, "code", "round63"))
import bridge_harness as bh
import jitter_score_diag_colab_frozen as jd
import physics, solvers

TAU = bh.TAU
SIDE = bh.SIDE
NU = 2000.0
C = 0.05
SCENES = ["bridge_contour_1", "bridge_microtex_1", "bridge_control_2"]
TARGET_LOADS = [1.0, 5.7, 22.25]
N_SEED = 8
OUT = os.path.join(ROOT, "results", "round63_next", "EFFICIENCY_PROBE")

# ---- J(rho;c) lookup (per-slot log-rate info) via the robust histogram-score
# estimator fi_hist (positive & stable at high load, where fi_mi's lam^2*Var term
# needs millions of frames; fi_mi/fi_cs/fi_hist agree to ~1-2% at N_MC=1e5). -----
DLOG = 0.01
ALPHA = 0.5
def build_J_grid(nu, c, grid_nmc=100_000):
    rhos = np.unique(np.concatenate([
        np.geomspace(0.05, 80.0, 45), np.array(TARGET_LOADS)]))
    mults = [np.exp(-DLOG), 1.0, np.exp(DLOG)]
    J = np.empty(rhos.shape[0])
    for i, rho in enumerate(rhos):
        _, _, cnt = jd.simulate_multi(float(rho), c, int(nu), grid_nmc,
                                      seed=20260722 + i, mults=mults)
        J[i] = jd.fi_hist(cnt[mults[0]], cnt[1.0], cnt[mults[2]],
                          int(nu), DLOG, ALPHA)
    assert np.all(J > 0), "non-positive J on grid: %s" % rhos[J <= 0]
    return rhos, J

def I_rho_of(loads, rhos, J, nu):
    """Per-frame Fisher info about the load rho at each operating load."""
    lr = np.log(np.clip(loads, rhos[0], rhos[-1]))
    Jv = np.interp(lr, np.log(rhos), J)
    rho_c = np.clip(loads, rhos[0], rhos[-1])
    return nu * Jv / (rho_c ** 2)          # I_rho = nu*J / rho^2

# ---- reconstruction (deployed RQL+TV, and RQL no-TV) -------------------------
def recon(A_full, N_full, nu, lam_tv, n_iter=200):
    det_est = physics.Detector(tau=TAU, dark=0.0)
    ctx = solvers.ArmContext(Phi=1.0 / TAU, det=det_est, T=nu * TAU, side=SIDE,
                             sigma_b=0.0, n_iter=n_iter, select_iter=60,
                             pattern_kind="bridge_main",
                             meta={"n_physical_rows": A_full.shape[0]},
                             lam_tv=lam_tv)
    xh, info = solvers.run_arm("RQL", A_full, N_full.astype(np.float64), ctx)
    return np.asarray(xh, dtype=np.float64).ravel(), float(info.get("lam_tv", np.nan))

def eig_decomp(A_full, I_rho):
    Ix = A_full.T @ (I_rho[:, None] * A_full)      # n x n  Fisher matrix
    mu, V = np.linalg.eigh(Ix)                      # ascending
    # numerical rank from A_full itself (design geometry)
    sv = np.linalg.svd(A_full, compute_uv=False)
    tol = sv.max() * max(A_full.shape) * np.finfo(float).eps
    rank = int((sv > tol).sum())
    order = np.argsort(mu)[::-1]                    # descending
    mu = mu[order]; V = V[:, order]
    rng_idx = np.arange(rank)                       # top-`rank` = range
    floor = float(np.sum(1.0 / mu[rng_idx]))
    return mu, V, rank, floor

def decomp_err(e, V, rank):
    coeff = V.T @ e
    range_sse = float(np.sum(coeff[:rank] ** 2))
    total_sse = float(e @ e)
    null_sse = total_sse - range_sse
    return total_sse, range_sse, null_sse

def run_cell(scene, target, rhos, J, nu, c):
    x = bh.load_scene(scene)
    rows = bh._l0_rows()
    mult = target / float((rows @ x).mean())
    main = rows * mult
    P = bh.prescan_matrix_load()
    A_full = np.vstack([P, main])
    loads_all = A_full @ x
    I_rho = I_rho_of(loads_all, rhos, J, nu)
    mu, V, rank, floor = eig_decomp(A_full, I_rho)
    n_pix = x.shape[0]
    out = {"scene": scene, "target_load": target, "nu": nu, "c": c,
           "n_rows": int(A_full.shape[0]), "rank": rank,
           "n_null": int(n_pix - rank),
           "main_load_q": [float(q) for q in np.percentile(main @ x, [5, 50, 95])],
           "floor_range_sse": floor,
           "floor_range_mse": floor / n_pix,
           "mu_min_range": float(mu[rank - 1]), "mu_max": float(mu[0])}
    for tag, lam_tv, n_iter in [("TV", None, 200), ("noTV", 0.0, 600)]:
        tot = np.zeros(N_SEED); rng_ = np.zeros(N_SEED); nul = np.zeros(N_SEED)
        lam_used = []
        for s in range(N_SEED):
            rp = bh._rng(s, 1, scene, int(nu), int(round(c * 1000)))
            rm = bh._rng(s, 2, scene, tag, int(nu), int(round(c * 1000)))
            N_pre = bh.simulate_counts_jitter(P @ x, nu, c, rp)
            N_main = bh.simulate_counts_jitter(main @ x, nu, c, rm)
            N_full = np.concatenate([N_pre, N_main])
            xh, lam = recon(A_full, N_full, nu, lam_tv, n_iter=n_iter)
            e = xh - x
            t, r, n = decomp_err(e, V, rank)
            tot[s], rng_[s], nul[s] = t, r, n
            lam_used.append(lam)
        out[tag] = {
            "lam_tv_mean": float(np.mean(lam_used)),
            "total_mse": float(tot.mean() / n_pix),
            "range_mse": float(rng_.mean() / n_pix),
            "null_mse": float(nul.mean() / n_pix),
            "range_sse": float(rng_.mean()),
            "null_sse": float(nul.mean()),
            "eff_range": floor / float(rng_.mean()),
            "range_headroom_dB": 10.0 * np.log10(float(rng_.mean()) / floor),
            "psnr_rad": 10.0 * np.log10(float(x.max()) ** 2 / float(tot.mean() / n_pix)),
        }
    return out

def main():
    t0 = time.time()
    print(f"[img] building J grid nu={NU} c={C} ...", flush=True)
    rhos, J = build_J_grid(NU, C)
    print(f"[img] J grid done ({time.time()-t0:.0f}s); "
          f"J range [{J.min():.4f},{J.max():.4f}] over rho[{rhos[0]:.3f},{rhos[-1]:.1f}]",
          flush=True)
    np.savez(os.path.join(OUT, "J_grid_c05_nu%d.npz" % int(NU)), rhos=rhos, J=J)
    cells = []
    for scene in SCENES:
        for target in TARGET_LOADS:
            r = run_cell(scene, target, rhos, J, NU, C)
            cells.append(r)
            tv = r["TV"]; nt = r["noTV"]
            print(f"[img] {scene:20s} rho={target:6.2f} rank={r['rank']} "
                  f"floor_mse={r['floor_range_mse']:.3e} | "
                  f"TV: rng={tv['range_mse']:.3e} null={tv['null_mse']:.3e} "
                  f"eff={tv['eff_range']:.3f} | noTV: rng={nt['range_mse']:.3e} "
                  f"eff={nt['eff_range']:.3f}  ({time.time()-t0:.0f}s)", flush=True)
    with open(os.path.join(OUT, "image_efficiency.json"), "w") as f:
        json.dump({"nu": NU, "c": C, "n_seed": N_SEED, "cells": cells}, f, indent=2)
    print("[img] DONE", flush=True)

if __name__ == "__main__":
    main()
