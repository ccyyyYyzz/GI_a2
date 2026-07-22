"""PROBE A image supplement: per-eigenband decomposition of the Fisher floor and
the deployed RQL+TV range error, to localize where the range-space MSE and its
CRB floor live (well-measured bulk vs ill-conditioned tail).  Reuses the saved
J grid.  Runs one representative operating point per scene (rho=22.25, the ridge)
plus rho=5.7, averaging the error over a few seeds.
"""
import json, os, sys, time
import numpy as np

ROOT = r"D:\GI_another"
sys.path.insert(0, os.path.join(ROOT, "code", "round63"))
import bridge_harness as bh
import physics, solvers

TAU = bh.TAU; SIDE = bh.SIDE
NU = 2000.0; C = 0.05
OUT = os.path.join(ROOT, "results", "round63_next", "EFFICIENCY_PROBE")
SCENES = ["bridge_contour_1", "bridge_microtex_1", "bridge_control_2"]
LOADS = [5.7, 22.25]
N_SEED = 6
NBAND = 4                      # quartile bands over the range eigen-spectrum

g = np.load(os.path.join(OUT, "J_grid_c05_nu%d.npz" % int(NU)))
RHOS, JG = g["rhos"], g["J"]

def I_rho_of(loads):
    lr = np.log(np.clip(loads, RHOS[0], RHOS[-1]))
    Jv = np.interp(lr, np.log(RHOS), JG)
    rc = np.clip(loads, RHOS[0], RHOS[-1])
    return NU * Jv / rc ** 2

def recon(A, N, lam_tv, n_iter=200):
    det = physics.Detector(tau=TAU, dark=0.0)
    ctx = solvers.ArmContext(Phi=1.0 / TAU, det=det, T=NU * TAU, side=SIDE,
                             sigma_b=0.0, n_iter=n_iter, select_iter=60,
                             pattern_kind="bridge_main",
                             meta={"n_physical_rows": A.shape[0]}, lam_tv=lam_tv)
    xh, info = solvers.run_arm("RQL", A, N.astype(np.float64), ctx)
    return np.asarray(xh, float).ravel(), float(info.get("lam_tv", np.nan))

def main():
    rows = bh._l0_rows(); P = bh.prescan_matrix_load()
    cells = []
    for scene in SCENES:
        x = bh.load_scene(scene)
        for target in LOADS:
            main = rows * (target / float((rows @ x).mean()))
            A = np.vstack([P, main])
            Irho = I_rho_of(A @ x)
            Ix = A.T @ (Irho[:, None] * A)
            mu, V = np.linalg.eigh(Ix)
            sv = np.linalg.svd(A, compute_uv=False)
            rank = int((sv > sv.max() * max(A.shape) * np.finfo(float).eps).sum())
            o = np.argsort(mu)[::-1]; mu = mu[o]; V = V[:, o]
            rng_mu = mu[:rank]
            # average deployed-TV error coeffs^2 over seeds
            c2 = np.zeros(A.shape[1]); lam_used = []
            for s in range(N_SEED):
                Np = bh.simulate_counts_jitter(P @ x, NU, C, bh._rng(s, 1, scene, int(NU), int(round(C*1000))))
                Nm = bh.simulate_counts_jitter(main @ x, NU, C, bh._rng(s, 2, scene, "TV", int(NU), int(round(C*1000))))
                xh, lam = recon(A, np.concatenate([Np, Nm]), None)
                lam_used.append(lam)
                coeff = V.T @ (xh - x)
                c2 += coeff ** 2
            c2 /= N_SEED
            # bands over range directions (index 0..rank-1, descending eigenvalue)
            edges = np.linspace(0, rank, NBAND + 1).astype(int)
            bands = []
            for bi in range(NBAND):
                sl = slice(edges[bi], edges[bi + 1])
                floor_b = float(np.sum(1.0 / rng_mu[sl]))
                err_b = float(np.sum(c2[sl]))
                bands.append({"idx": [int(edges[bi]), int(edges[bi+1])],
                              "mu_lo": float(rng_mu[sl][-1]), "mu_hi": float(rng_mu[sl][0]),
                              "floor_sse": floor_b, "err_sse": err_b,
                              "eff": floor_b / err_b if err_b > 0 else float('inf')})
            cell = {"scene": scene, "target_load": target, "rank": rank,
                    "lam_tv_mean": float(np.mean(lam_used)),
                    "floor_total": float(np.sum(1.0 / rng_mu)),
                    "range_err_total": float(np.sum(c2[:rank])),
                    "null_err_total": float(np.sum(c2[rank:])),
                    "bands": bands}
            cells.append(cell)
            print(f"{scene:20s} rho={target:5.2f} rank={rank} floor={cell['floor_total']:.3e} "
                  f"rngerr={cell['range_err_total']:.3e} nullerr={cell['null_err_total']:.3e}", flush=True)
            for b in bands:
                print(f"   band {b['idx']} mu[{b['mu_lo']:.2e},{b['mu_hi']:.2e}] "
                      f"floor={b['floor_sse']:.3e} err={b['err_sse']:.3e} eff={b['eff']:.3f}", flush=True)
    json.dump({"nu": NU, "c": C, "n_seed": N_SEED, "cells": cells},
              open(os.path.join(OUT, "image_eigenbands.json"), "w"), indent=2)
    print("BANDS DONE", flush=True)

if __name__ == "__main__":
    main()
