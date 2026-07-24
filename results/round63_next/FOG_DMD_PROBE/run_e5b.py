# E5b -- SPARSITY-CONSTRAINED PATHWISE. Does a structural sparsity constraint break the fiber
# collusion that made the E4 pathwise refinement drift away from truth? Sparse scene (few bright
# squares -> sparse in pixels & gradient); reduced-Q ALS pathwise refinement with L1 soft-threshold
# + continuation. Key test: does true-Z-seeded refinement now STAY at truth? Same table as sec11c.
import time, json, numpy as np, torch
from fog_common import *
import fog_tracker as ft
def log(*a): print(*a, flush=True)
dev = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
n = 32; N = n * n; M = 96; d = 48; sf = 0.15; TAU = 16.0

def sparse_scene(g):
    X = np.zeros((n, n))
    for _ in range(g.integers(3, 6)):
        r0 = g.integers(0, n - 5); c0 = g.integers(0, n - 5)
        X[r0:r0 + g.integers(2, 5), c0:c0 + g.integers(2, 5)] = g.uniform(.5, 1)
    return (X / X.max()).ravel()

def refine_sparse(P, U, B, z_init, l1_sched=(3e-2, 1e-2, 3e-3, 1e-3), iters=240, dev='cuda'):
    dev = torch.device(dev)
    P = torch.tensor(P, device=dev, dtype=torch.float64); U = torch.tensor(U, device=dev, dtype=torch.float64)
    B = torch.tensor(B, device=dev, dtype=torch.float64)
    Md, Nn = P.shape; dd = U.shape[1]; T = B.shape[0]
    PtP = P.t() @ P; PtBt = P.t() @ B.t(); Ut = U.t()
    Bc = B - B.mean(0, keepdim=True); _, _, Vh = torch.linalg.svd(Bc.t(), full_matrices=False); Zh = Vh[:dd].t()
    ZtZ_inv = torch.linalg.inv(Zh.t() @ Zh + 1e-9 * torch.eye(dd, device=dev, dtype=torch.float64))
    zt = torch.tensor(z_init, device=dev, dtype=torch.float64)
    Q = ZtZ_inv @ (Zh.t() @ (zt / (zt.std(0, keepdim=True) + 1e-12)))
    x = None
    for lam1 in l1_sched:
        for it in range(iters // len(l1_sched)):
            W = torch.clamp(1.0 + (Zh @ Q) @ Ut, min=0.05)
            x = ft.solve_x_clean(PtP, PtBt, W, 1e-5)
            x = torch.clamp(x - lam1, min=0.0)                # nonneg soft-threshold (L1 prox)
            H = P @ (U * x[:, None]); HtH = H.t() @ H + 1e-8 * torch.eye(dd, device=dev, dtype=torch.float64)
            y = B - (P @ x)[None, :]
            Q = ZtZ_inv @ (Zh.t() @ y @ H) @ torch.linalg.inv(HtH)
    return x.detach().cpu().numpy()

def build(seed):
    g = np.random.default_rng(6060 + seed)
    x_true = sparse_scene(g)
    P = make_patterns(M, N, "binary", g); Pi, rangeP = projectors(P)
    nt = x_true - rangeP(x_true); nrm0 = np.linalg.norm(nt)
    U = dct_basis(d, n); W, Z, rho, csd = lognormal_medium(U, T, sf, TAU, g)
    B = np.einsum('mn,tn->tm', P, W * x_true)
    return dict(x_true=x_true, P=P, U=U, W=W, Z=Z, B=B, csd=csd, rho=rho, rangeP=rangeP, nt=nt, nrm0=nrm0,
                nullfrac=np.linalg.norm(nt) ** 2 / np.linalg.norm(x_true) ** 2)

T = 128
def nm(x, c): return null_metric(x, c['x_true'], c['rangeP'], c['nt'], c['nrm0'])[0] ** 2
rows = []
log("E5b: sparse scene, T=128, single-tau OU. Does sparsity keep true-Z-seed at truth?")
for kind in ['dense_ref', 'sparse']:
    orac, froz, tz_plain, tz_sparse, base_sparse = [], [], [], [], []
    for sd in range(5):
        if kind == 'dense_ref':
            # dense scene control (reuses E4-style scene) to contrast
            g = np.random.default_rng(6060 + sd); xt = make_scene(n)
            P = make_patterns(M, N, "binary", g); Pi, rangeP = projectors(P)
            nt = xt - rangeP(xt); nrm0 = np.linalg.norm(nt); U = dct_basis(d, n)
            W, Z, rho, csd = lognormal_medium(U, T, sf, TAU, g); B = np.einsum('mn,tn->tm', P, W * xt)
            c = dict(x_true=xt, P=P, U=U, W=W, Z=Z, B=B, csd=csd, rho=rho, rangeP=rangeP, nt=nt, nrm0=nrm0)
        else:
            c = build(sd)
        orac.append(nm(ft.solve_oracle(c['P'], c['W'], c['B'], np.ones_like(c['B']), 1e-6, str(dev)), c))
        Wl = np.clip(1.0 + c['Z'] @ c['U'].T, 0.05, None)
        froz.append(nm(ft.solve_oracle(c['P'], Wl, c['B'], np.ones_like(c['B']), 1e-6, str(dev)), c))
        # true-Z-seed WITHOUT sparsity (plain reduced ALS) vs WITH sparsity
        tz_plain.append(nm(refine_sparse(c['P'], c['U'], c['B'], c['Z'], l1_sched=(0,), iters=200, dev=str(dev)), c))
        tz_sparse.append(nm(refine_sparse(c['P'], c['U'], c['B'], c['Z'], dev=str(dev)), c))
        rnd = np.random.default_rng(sd).standard_normal((T, d))
        base_sparse.append(nm(refine_sparse(c['P'], c['U'], c['B'], rnd, dev=str(dev)), c))
    row = dict(kind=kind, oracle=float(np.median(orac)), frozen_known=float(np.median(froz)),
               trueZ_noSparse=float(np.median(tz_plain)), trueZ_sparse=float(np.median(tz_sparse)),
               baseline_sparse=float(np.median(base_sparse)))
    rows.append(row)
    log(f"[{kind:9s}] oracle {row['oracle']:.3f} | frozen {row['frozen_known']:.3f} | "
        f"true-Z(noL1) {row['trueZ_noSparse']:.3f} | true-Z(+L1) {row['trueZ_sparse']:.3f} | "
        f"baseline(+L1) {row['baseline_sparse']:.3f}")

json.dump(dict(rows=rows), open('E5b_results.json', 'w'), indent=2)
sp = [r for r in rows if r['kind'] == 'sparse'][0]
broke = sp['trueZ_sparse'] < 0.5 * sp['trueZ_noSparse']
log(f"\nSPARSITY_BREAKS_COLLUSION: {'YES' if broke else 'NO'} "
    f"(true-Z: {sp['trueZ_noSparse']:.3f} noL1 -> {sp['trueZ_sparse']:.3f} +L1)")
