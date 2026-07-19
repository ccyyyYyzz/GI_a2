"""DEV gap discrimination probe (coordinator directive, R18 prep).

For the two DEV certificate cells: (1) tightened dual bound (full mu-pixel
set, MU_CAP 1e9, 25 CG rounds, 240s LP budget); (2) PRIMAL probe: best
feasible design in C_dose over D_cert via multiplicative ascent + dose/
budget projection, tracking the best FEASIBLE log-det vs deployed SCAT32;
(3) budget-binding checks for the gamma atoms. DEV scenes only.
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, r"D:/GI_another/code/round63")
sys.path.insert(0, r"D:/GI_another/code")
os.chdir(r"D:/GI_another")

import campaign                              # noqa: E402
import m1_runner as m1                       # noqa: E402
import oed_design_v4 as v4                   # noqa: E402
import oed_design_v5 as v5                   # noqa: E402

DEV_IMG = "m1_dev_glyph"
SEED = 0


def deployed_logdet(ctx, rows_rel, b):
    xhat = ctx["xhat"]
    loads_rel = rows_rel @ xhat
    rows_abs = rows_rel * (b / float(loads_rel.mean())) * (1 - 1e-12)
    loads = rows_abs @ xhat
    V = ctx["V0"].copy()
    for s in range(rows_abs.shape[0]):
        idx = np.nonzero(rows_abs[s])[0]
        q = ctx["B"][idx].T @ (rows_abs[s, idx] / loads[s])
        V += ctx["nu"] * v4._J(float(loads[s]), ctx["nu"]) * np.outer(q, q)
    V[np.diag_indices(ctx["r"])] += ctx["eps0"]
    return float(np.linalg.slogdet(V)[1]), float(loads.mean())


def primal_probe(ctx, b, iters=80):
    """Multiplicative ascent + dose/budget projection over the FULL D_cert;
    returns the best FEASIBLE (budget + dose band) log-det + composition."""
    xi = np.where(ctx["ALLOW"], 1.0, 0.0)
    lo = int(np.argmax(ctx["ALLOW"].any(axis=0)))
    xi = np.zeros_like(xi)
    xi[ctx["ALLOW"][:, lo], lo] = 1.0
    xi /= xi.sum()
    xi = v5._dose_project(ctx, xi, b, v5.DELTA)
    best = (-np.inf, None, -1)
    for it in range(iters):
        V = v5.assemble(ctx, xi)
        Vinv = np.linalg.inv(V)
        D = v5.dvals(ctx, Vinv)
        pos = np.where(ctx["ALLOW"], np.maximum(D, 0.0), 0.0)
        zb = float((xi * pos).sum())
        if zb <= 0:
            break
        xi = xi * np.sqrt(pos / zb)
        xi[~ctx["ALLOW"]] = 0.0
        xi /= xi.sum()
        xi = v5._dose_project(ctx, xi, b, v5.DELTA)
        load = float((xi * ctx["C"]).sum())
        dv = v5.dose_of(ctx, xi)
        dev = float(np.abs(dv / dv.mean() - 1.0).max())
        if load <= b + 1e-9 and dev <= v5.DELTA:
            ld = float(np.linalg.slogdet(v5.assemble(ctx, xi))[1])
            if ld > best[0]:
                best = (ld, xi.copy(), it)
        if it % 20 == 0:
            print("      [probe] it=%d load=%.4f dev=%.4f best_ld=%s"
                  % (it, load, dev,
                     ("%.3f" % best[0]) if np.isfinite(best[0]) else "-"),
                  flush=True)
    return best


def main():
    dev_imgs = campaign._images(32, "all", imageset="m1_dev")
    x_true = dev_imgs[DEV_IMG]
    rows, sha = m1.deployed_scat32()
    P = m1.prescan_matrix()
    for (nu, b) in ((2000.0, 0.60), (200.0, 0.05)):
        print("\n===== cell (%s, s%d, nu=%g, b=%g) =====" % (
            DEV_IMG, SEED, nu, b), flush=True)
        t0 = time.time()
        xhat = m1.prescan_estimate(x_true, DEV_IMG, SEED, b, nu,
                                   per_cell=True)
        V_full = v4.info_matrix_full(rows, xhat, int(nu), b, P=P)
        B, eps0, _tr = v4.subspace_from_fixedstar(V_full)
        ctx = v5.setup_ctx_cert(xhat, nu, b, B, eps0, 32)
        ld_dep, mload = deployed_logdet(ctx, rows, b)
        print("  setup %.0fs; deployed SCAT32: logdet=%.4f mean_load=%.6f"
              % (time.time() - t0, ld_dep, mload), flush=True)
        # ---- gamma-atom budget audit ---------------------------------- #
        Ll = ctx["L_load"]
        Cg = ctx["C"][:, Ll:]
        Ag = ctx["ALLOW"][:, Ll:]
        print("  gamma-atom loads (admissible): min=%.4f med=%.4f max=%.4f"
              " (budget b=%g); gamma=5 admissible: %d atoms, load range"
              " [%.3f, %.3f]"
              % (Cg[Ag].min(), np.median(Cg[Ag]), Cg[Ag].max(), b,
                 int(Ag[:, -1].sum()),
                 Cg[Ag[:, -1], -1].min() if Ag[:, -1].any() else -1,
                 Cg[Ag[:, -1], -1].max() if Ag[:, -1].any() else -1),
              flush=True)
        # ---- 1) tightened dual ---------------------------------------- #
        t1 = time.time()
        out = v5.cert_deployed_rows(ctx, rows, b, verbose=True,
                                    kw_rounds=25, kw_npix=1024,
                                    kw_mucap=1e9, kw_lptime=240)
        print("  TIGHT dual: status=%s G_full/r=%s theta=%.4g "
              "MU_CAP_ACTIVE=%s n_active_mu=%s wall=%.0fs"
              % (out["status"],
                 ("%.4f" % (out["G_full"] / ctx["r"]))
                 if np.isfinite(out.get("G_full", np.inf)) else "inf",
                 out.get("theta", -1), out.get("MU_CAP_ACTIVE"),
                 out.get("n_active_mu"), time.time() - t1), flush=True)
        # ---- 2) primal probe ------------------------------------------ #
        t2 = time.time()
        ld_best, xi_best, it_best = primal_probe(ctx, b)
        gap_primal = (ld_best - ld_dep) / ctx["r"] if np.isfinite(ld_best) \
            else float("nan")
        print("  PRIMAL probe: best feasible logdet=%.4f (it=%d, %.0fs); "
              "SCAT32=%.4f; primal gap/r = %.4f"
              % (ld_best, it_best, time.time() - t2, ld_dep, gap_primal),
              flush=True)
        if xi_best is not None:
            mass_load = float(xi_best[:, :Ll].sum())
            mass_gain = float(xi_best[:, Ll:].sum())
            loads_used = ctx["C"][xi_best > 1e-6]
            print("  probe composition: D_load mass=%.3f D_gain mass=%.3f; "
                  "support loads q5/q50/q95 = %.3f/%.3f/%.3f; "
                  "budget used=%.4f (cap %g)"
                  % (mass_load, mass_gain,
                     np.percentile(loads_used, 5),
                     np.percentile(loads_used, 50),
                     np.percentile(loads_used, 95),
                     float((xi_best * ctx["C"]).sum()), b), flush=True)
    print("\n[probe] DONE", flush=True)


if __name__ == "__main__":
    main()
