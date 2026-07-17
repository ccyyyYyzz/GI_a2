"""ROUND63 solver-suite unit tests (fast; whole file < 120 s).

Operating point (spec §2/§3 pilot cell): 32x32 phantom, Bernoulli-50% patterns,
M=3000 frames, non-paralyzable detector, mean load rho = tau*E[lambda] = 3/7,
nu = T/tau = 500 (=> ~150 recorded counts/frame, clearly sub-Poisson). Phi is
pinned so the EMPIRICAL mean load equals 3/7 exactly.

Checks (PASS/FAIL printed per check):
  (0) grad/div adjoint + tv_prox sanity (solver building blocks)
  (a) tv_fista converges: monotone objective, small gradient-mapping residual
  (b) QMLE beats POISSON-LIN by > 1 dB flux-matched PSNR at this operating point
  (c) SAT-POISSON sits between POISSON-LIN and QMLE
  (d) hadamard_pair_combine round-trips exactly on a hand-built 8x8 Hadamard case
  (e) precorrect_rates + variance-weighted-LS arm runs and is finite

Run:  D:/Anacondar/anaconda3/python.exe code/round63/test_solvers.py
(the script's own dir is auto on sys.path, so `import physics/solvers` resolve;
gi_core is not needed here.)
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import physics
import solvers as S

SEED0 = 20260717     # frozen RNG base (mirrors gi_core.config; utils.rng_for)
SIDE = 32
N = SIDE * SIDE
M = 3000
TAU = 1.0            # dimensionless dead time (only tau*lambda and T/tau matter)
NU = 500
T = NU * TAU
RHO = 3.0 / 7.0


def rng_for(seed, *stream):
    """Local mirror of gi_core.utils.rng_for (avoids adding code/ to path)."""
    from numpy.random import default_rng
    return default_rng([SEED0 + int(seed)] + [int(t) for t in stream])


def make_phantom(side):
    """High-contrast sum-normalized phantom (disc + square + bright bar + lines).
    High contrast -> real per-frame pattern-energy spread -> the dead-time mean
    compression differs frame-to-frame, which is what QMLE models and POISSON-LIN
    cannot (a pure global scale would be erased by flux matching)."""
    img = np.full((side, side), 0.03, dtype=np.float64)
    yy, xx = np.mgrid[0:side, 0:side]
    disc = (yy - 9) ** 2 + (xx - 9) ** 2 <= 25
    img[disc] = 1.0
    img[20:28, 4:12] = 0.55
    img[6:9, 16:30] = 1.0            # bright bar
    img[16:30, 22] = 0.85            # thin line
    img[24, 15:28] = 0.7
    x = img.ravel()
    return x / x.sum()


def flux_psnr(xhat, x_true, side):
    """MAIN-protocol flux-matched PSNR (mirrors gi_core.metrics)."""
    xp = np.maximum(xhat, 0.0)
    s = xp.sum()
    xs = xp * (x_true.sum() / s) if s > 0 else np.zeros_like(xp)
    dr = float(x_true.max())
    mse = float(np.mean((xs - x_true) ** 2))
    return np.inf if mse <= 0 else 10.0 * np.log10(dr ** 2 / mse)


def bern_patterns(M, n, rng):
    return (rng.random((M, n)) < 0.5).astype(np.float64)


PASS = []


def report(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    PASS.append(bool(ok))
    print("[%s] %s %s" % (tag, name, detail), flush=True)


# ----------------------------------------------------------------------
def check_building_blocks():
    rng = rng_for(0, 63, 1)
    m, n = 9, 7
    x = rng.standard_normal((m, n))
    p = rng.standard_normal((m, n))
    q = rng.standard_normal((m, n))
    gx, gy = S._grad(x)
    lhs = np.sum(gx * p) + np.sum(gy * q)
    rhs = -np.sum(x * S._div(p, q))
    adj_err = abs(lhs - rhs) / max(abs(lhs), 1e-12)

    # tv_prox: w=0 is nonneg projection; large w -> near-flat, lower TV, feasible
    b = rng.standard_normal((16, 16))
    p0 = S.tv_prox(b, 0.0)
    zero_ok = np.allclose(p0, np.maximum(b, 0.0))
    pw = S.tv_prox(b, 5.0)
    tv_drop = S.tv_value(pw, 16) < S.tv_value(np.maximum(b, 0.0), 16)
    feasible = pw.min() >= -1e-9
    ok = adj_err < 1e-10 and zero_ok and tv_drop and feasible
    report("(0) grad/div adjoint + tv_prox", ok,
           "adj_err=%.2e tv_drop=%s feasible=%s" % (adj_err, tv_drop, feasible))


def simulate_cell():
    """Build the shared (A, b, x_true, Phi, det) pilot cell."""
    x_true = make_phantom(SIDE)
    A = bern_patterns(M, N, rng_for(0, 63, 10))
    u = A @ x_true
    ubar = float(u.mean())
    det = physics.Detector(tau=TAU, dark=0.0, start_mode="active")
    Phi = RHO / (TAU * ubar)                       # pins empirical mean rho = 3/7
    b, Ncount = physics.simulate_counts(u, Phi, T, det, rng_for(0, 63, 20))
    return dict(x_true=x_true, A=A, b=b, N=Ncount, Phi=Phi, det=det,
                lam=Phi * u, rho_emp=float(np.mean(Phi * u) * TAU),
                mean_counts=float(Ncount.mean()))


def check_fista_converges(cell):
    ctx = S.ArmContext(Phi=cell["Phi"], det=cell["det"], T=T, side=SIDE)
    x0 = S.init_gi_flux(cell["A"], cell["b"], cell["Phi"], 0.0, T, TAU)
    fg = S._phys_factory("RQL", ctx)(cell["A"], cell["b"])
    # meaningful regularization (grid-max of the auto-scaled selection grid)
    _, g0 = fg(x0)
    lam = 1e-2 * float(np.linalg.norm(g0))
    n_max = 400
    x, info = S.tv_fista(fg, x0, lam, n_iter=n_max, side=SIDE)
    hist = np.asarray(info["obj_hist"])
    # MFISTA is monotone up to prox round-off
    max_rise = float(np.max(np.diff(hist))) if hist.size > 1 else 0.0
    monotone = max_rise <= 1e-8 * max(1.0, abs(hist[0]))
    total_drop = hist[0] - hist[-1]
    # objective converged: the rel-change < 1e-7 stop fired before the iter cap
    stopped = info["n_iter"] < n_max
    # stationarity residual (gradient mapping) collapsed >=10x vs the initializer
    gmap_ratio = info["grad_map_norm"] / max(info["grad_map_norm_x0"], 1e-30)
    ok = monotone and total_drop > 0 and stopped and gmap_ratio < 0.1
    report("(a) tv_fista converges", ok,
           "iters=%d(<%d=%s) obj %.4f->%.4f max_rise=%.1e gmap %.2e->%.2e (x%.3f)"
           % (info["n_iter"], n_max, stopped, hist[0], hist[-1], max_rise,
              info["grad_map_norm_x0"], info["grad_map_norm"], gmap_ratio))


def check_mechanism_ablation(cell):
    """(b) QMLE > POISSON-LIN by >1 dB and (c) SAT-POISSON in between. All arms:
    same solver, same x0, same budget, per-arm truth-free lam_tv selection."""
    A, b, x_true = cell["A"], cell["b"], cell["x_true"]
    # lam_tv chosen per-arm by the truth-free held-out-NLL rule; the grid
    # auto-calibrates to each arm's ||grad f_data(x0)||_2 (lam_grid_scale=None).
    x0 = S.init_gi_flux(A, b, cell["Phi"], 0.0, T, TAU)
    ctx = S.ArmContext(Phi=cell["Phi"], det=cell["det"], T=T, side=SIDE,
                       lam_tv=None, n_iter=180, x0=x0)
    psnr, lam = {}, {}
    for arm in ("POISSON-LIN", "SAT-POISSON", "QMLE"):
        x, info = S.run_arm(arm, A, b, ctx)
        psnr[arm] = flux_psnr(x, x_true, SIDE)
        lam[arm] = info["lam_tv"]
    gap = psnr["QMLE"] - psnr["POISSON-LIN"]
    ok_b = gap > 1.0
    report("(b) QMLE beats POISSON-LIN by >1 dB", ok_b,
           "QMLE=%.2f LIN=%.2f gap=%.2f dB (lam_tv=%.2g/%.2g)"
           % (psnr["QMLE"], psnr["POISSON-LIN"], gap, lam["QMLE"], lam["POISSON-LIN"]))
    # SAT sits between (with a small tolerance for solver/regularizer noise)
    lo, hi = psnr["POISSON-LIN"], psnr["QMLE"]
    tol = 0.35
    ok_c = (lo - tol) <= psnr["SAT-POISSON"] <= (hi + tol)
    report("(c) SAT-POISSON between LIN and QMLE", ok_c,
           "LIN=%.2f SAT=%.2f QMLE=%.2f" % (lo, psnr["SAT-POISSON"], hi))
    return psnr


def check_hadamard_roundtrip():
    """(d) hand-built 8x8 Hadamard complementary-pair round-trip."""
    H = np.array([[1, 1, 1, 1, 1, 1, 1, 1],
                  [1, -1, 1, -1, 1, -1, 1, -1],
                  [1, 1, -1, -1, 1, 1, -1, -1],
                  [1, -1, -1, 1, 1, -1, -1, 1],
                  [1, 1, 1, 1, -1, -1, -1, -1],
                  [1, -1, 1, -1, -1, 1, -1, 1],
                  [1, 1, -1, -1, -1, -1, 1, 1],
                  [1, -1, -1, 1, -1, 1, 1, -1]], dtype=np.float64)
    assert np.allclose(H @ H.T, 8 * np.eye(8))
    # realize each +/-1 row as a 0/1 pattern and its complement
    A = np.zeros((16, 8), dtype=np.float64)
    for k in range(8):
        a_plus = (H[k] + 1.0) / 2.0
        A[2 * k] = a_plus
        A[2 * k + 1] = 1.0 - a_plus
    x_true = rng_for(0, 63, 30).random(8) + 0.1
    b = A @ x_true                              # noiseless linear buckets
    meta = {"pairs": np.array([[2 * k, 2 * k + 1] for k in range(8)])}
    A_signed, b_signed = S.hadamard_pair_combine(A, b, meta)
    signed_ok = np.allclose(A_signed, H) and np.allclose(b_signed, H @ x_true)
    x_rec = (H.T @ b_signed) / 8.0              # H^T H = 8 I -> exact inverse
    recon_ok = np.allclose(x_rec, x_true, atol=1e-10)
    # the dispatcher routes hadpair -> the differenced signed system; the
    # (centered) GI arm on it is a biased estimator, so require a strong
    # correlation, not exact recovery (the exact round-trip is checked above).
    ctx = S.ArmContext(Phi=1.0, det=physics.Detector(tau=TAU), T=T, side=None,
                       pattern_kind="hadpair", meta=meta)
    x_gi, ginfo = S.run_arm("GI", A, b, ctx)
    corr = np.corrcoef(x_gi, x_true)[0, 1]
    ok = signed_ok and recon_ok and ginfo["hadpair"] and corr > 0.8
    report("(d) hadamard pair combine round-trips", ok,
           "signed=%s recon_max_err=%.2e dispatch_gi_corr=%.4f"
           % (signed_ok, float(np.max(np.abs(x_rec - x_true))), corr))


def check_precorrect(cell):
    """(e) precorrect_rates + WLS arm runs and is finite (via the dispatcher,
    truth-free lam_tv selection); also exercise the paralyzable Lambert-W
    inversion path."""
    A, b, x_true = cell["A"], cell["b"], cell["x_true"]
    x0 = S.init_gi_flux(A, b, cell["Phi"], 0.0, T, TAU)
    ctx = S.ArmContext(Phi=cell["Phi"], det=cell["det"], T=T, side=SIDE,
                       lam_tv=None, n_iter=180, x0=x0)
    x, info = S.run_arm("PRECORRECT", A, b, ctx)
    finite = np.all(np.isfinite(x)) and np.isfinite(info["final_obj"])
    ps = flux_psnr(x, x_true, SIDE)
    # paralyzable inversion path: recorded rate r=lam*exp(-lam*tau) -> lam via W0
    lam_true = 0.4
    r = lam_true * np.exp(-lam_true * TAU)
    lam_par = S._precorrect_rates_paralyzable(np.array([r * T]), T, TAU)[0]
    par_ok = np.isfinite(lam_par) and abs(lam_par - lam_true) < 1e-6
    ok = finite and ps > 0 and par_ok
    report("(e) precorrect_rates + WLS arm", ok,
           "PSNR=%.2f lam_tv=%.3g finite=%s lambertW_inv_err=%.1e"
           % (ps, info["lam_tv"], finite, abs(lam_par - lam_true)))


def check_exact_reference():
    """(f) small-scale EXACT reference arm runs (smooth exact renewal NLL, nonneg
    L-BFGS-B, no TV) + the dispatcher rejects unknown arms. 8x8 / 300 frames so
    the per-frame Gamma-CDF likelihood stays fast."""
    side, n, m, T8 = 8, 64, 300, 200.0
    img = np.zeros((side, side))
    img[2:6, 2:6] = 1.0
    x_true = img.ravel()
    x_true = x_true / x_true.sum()
    A = bern_patterns(m, n, rng_for(0, 63, 77))
    u = A @ x_true
    det = physics.Detector(tau=TAU, dark=0.0, start_mode="active")
    Phi = RHO / (TAU * float(u.mean()))
    b, _ = physics.simulate_counts(u, Phi, T8, det, rng_for(1, 63, 77))
    ctx = S.ArmContext(Phi=Phi, det=det, T=T8, side=side, exact_maxiter=40)
    x, info = S.run_arm("EXACT", A, b, ctx)
    ps = flux_psnr(x, x_true, side)
    exact_ok = (np.all(np.isfinite(x)) and info["reference_only"]
                and np.isfinite(info["fun"]) and ps > 0)
    try:
        S.run_arm("BOGUS", A, b, ctx)
        guard_ok = False
    except ValueError:
        guard_ok = True
    report("(f) EXACT reference arm + dispatch guard", exact_ok and guard_ok,
           "PSNR=%.2f nll=%.4f finite=%s unknown_arm_raises=%s"
           % (ps, info["fun"], np.all(np.isfinite(x)), guard_ok))


def main():
    t0 = time.time()
    np.seterr(over="ignore", invalid="ignore")
    check_building_blocks()
    cell = simulate_cell()
    print("cell: rho_emp=%.4f mean_counts=%.1f Phi=%.4g"
          % (cell["rho_emp"], cell["mean_counts"], cell["Phi"]), flush=True)
    check_fista_converges(cell)
    check_mechanism_ablation(cell)
    check_hadamard_roundtrip()
    check_precorrect(cell)
    check_exact_reference()
    dt = time.time() - t0
    n_pass = sum(PASS)
    print("\n==== %d/%d checks PASS  (%.1fs) ===="
          % (n_pass, len(PASS), dt), flush=True)
    return 0 if n_pass == len(PASS) else 1


if __name__ == "__main__":
    sys.exit(main())
