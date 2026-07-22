"""
T4 -- RECONSTRUCTION TEASER (descriptive; NO bar).

NOTE ON GATING: T4 is specified "only if T1-T3 pass".  T2 and T3 do NOT pass the
photon-cost bars, so this is run only as a NOISELESS UPSIDE CEILING: assuming the
51 quadratic constraints q_i = sum_p m_ip x_p^2 were measured EXACTLY (ignoring
the T2/T3 photon cost), how much would they improve reconstruction over linear-
only TV at 5% sampling?  This answers "is the quadratic side-information worth
anything at all, best case?" -- it does NOT claim it is measurable cheaply.

Solver: projected gradient on  TV_eps(x) + lam_lin||Mx-b||^2 + lam_q||Q(x)-q||^2,
x clamped to [0,1].  Linear-only sets lam_q=0.
"""
import json
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import saturation_core as sc  # noqa
import gi_operator as op

OUT = os.path.dirname(HERE)
DATA = "data/r63_bridge_scenes"
side = op.SIDE

M = op.build_operator()
SCENES = ["bridge_twopop_1", "bridge_contour_0", "bridge_microtex_2"]


def tv_and_grad(x2d, eps=1e-3):
    dx = np.zeros_like(x2d); dy = np.zeros_like(x2d)
    dx[:, :-1] = x2d[:, 1:] - x2d[:, :-1]
    dy[:-1, :] = x2d[1:, :] - x2d[:-1, :]
    mag = np.sqrt(dx ** 2 + dy ** 2 + eps ** 2)
    tv = mag.sum()
    gx = dx / mag; gy = dy / mag
    g = np.zeros_like(x2d)
    g[:, :-1] -= gx[:, :-1]; g[:, 1:] += gx[:, :-1]
    g[:-1, :] -= gy[:-1, :]; g[1:, :] += gy[:-1, :]
    return tv, g


def reconstruct(b, q, use_quad, iters=4000, lr=0.05, lam_lin=50.0, lam_q=50.0):
    x = np.full(side * side, 0.5)
    for it in range(iters):
        r_lin = M @ x - b
        g = 2 * lam_lin * (M.T @ r_lin)
        if use_quad:
            Qx = M @ (x ** 2)
            r_q = Qx - q
            g += 2 * lam_q * (2 * x * (M.T @ r_q))
        tv, gtv = tv_and_grad(x.reshape(side, side))
        g += 0.02 * gtv.ravel()
        x = np.clip(x - lr / (1 + 0.001 * it) * g, 0.0, 1.0)
    return x


def psnr(a, b):
    mse = np.mean((a - b) ** 2)
    return 99.0 if mse == 0 else 10 * np.log10(1.0 / mse)


report = {"dense": {}, "note": "NOISELESS upside ceiling only; not a photon-cost claim"}

# (a) DENSE scenes: linear-only TV vs linear+quad TV -> expect ~0 (fiber still open)
d_psnr = []
print("(a) DENSE scenes (5% sampling, TV):")
print("scene                PSNR_linear  PSNR_lin+quad   deltaPSNR")
for name in SCENES:
    x = np.clip(op.load_scene(os.path.join(DATA, name + ".npz")), 0, 1)
    b = M @ x; q = M @ (x ** 2)
    xl = reconstruct(b, q, use_quad=False)
    xq = reconstruct(b, q, use_quad=True)
    pl, pq = psnr(x, xl), psnr(x, xq)
    d_psnr.append(pq - pl)
    report["dense"][name] = dict(psnr_linear=float(pl), psnr_lin_quad=float(pq),
                                 delta_psnr=float(pq - pl))
    print("%-20s %10.3f  %12.3f  %10.3f" % (name, pl, pq, pq - pl))
report["dense_mean_delta_psnr"] = float(np.mean(d_psnr))
print("  dense mean delta PSNR = %.3f dB (fiber still ~90%% of N -> no dense gain)"
      % report["dense_mean_delta_psnr"])

# (b) LOW-DIM PRIOR: KNOWN support of size s with K < s <= 2K.  Linear-only is
# underdetermined (s>K); linear+quad is over-determined (2K>=s) -> recovers.
print("\n(b) KNOWN-SUPPORT sparse recovery (the honest positive teaser, noiseless):")
rng = np.random.default_rng(4242)
rec = {}
for s in (40, 80, 100, 120):
    xt = np.zeros(op.P)
    supp = rng.choice(op.P, size=s, replace=False)
    xt[supp] = rng.uniform(0.2, 1.0, size=s)
    b = M @ xt; q = M @ (xt ** 2)
    Ms = M[:, supp]
    # linear-only least-norm on support
    xl_s = np.linalg.lstsq(Ms, b, rcond=None)[0]
    # joint Gauss-Newton on support: solve [Ms; Ms diag] -> use LS iterations
    xs = np.full(s, 0.5)
    for _ in range(200):
        r1 = Ms @ xs - b
        r2 = Ms @ (xs ** 2) - q
        Jg = np.vstack([Ms, 2 * Ms * xs[None, :]])
        res = np.concatenate([r1, r2])
        dxs = np.linalg.lstsq(Jg, res, rcond=None)[0]
        xs = np.clip(xs - dxs, 0, 1.5)
    xl_full = np.zeros(op.P); xl_full[supp] = xl_s
    xq_full = np.zeros(op.P); xq_full[supp] = xs
    pl, pq = psnr(xt, xl_full), psnr(xt, xq_full)
    rec[f"s{s}"] = dict(support=s, psnr_linear_only=float(pl),
        psnr_joint=float(pq), linear_determined=bool(s <= op.K),
        joint_determined=bool(s <= 2 * op.K))
    print("  s=%3d (K=%d,2K=%d): PSNR linear-only=%6.2f  joint(quad)=%6.2f  %s"
          % (s, op.K, 2 * op.K, pl, pq,
             "quad RECOVERS" if pq > pl + 10 else "both fail" if pq < 15 else ""))
report["known_support"] = rec

with open(os.path.join(OUT, "t4_reconstruction.json"), "w") as fh:
    json.dump(report, fh, indent=2, default=float)
print("wrote t4_reconstruction.json")
